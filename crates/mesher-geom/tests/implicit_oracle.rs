//! The implicit predicates against an exact rational reference: each point's
//! coordinates divided out in `BigRational` (and checked to lie where its
//! definition says), then the textbook formulas. Random, degenerate and nearly
//! degenerate cases; 100,000 of each in release builds (2,000 in debug).
//! Also the filter rate: the f64 filter must settle at least 95% of random
//! queries.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::needless_range_loop
)]

use mesher_geom::implicit::{
    self, Kind, Point, Stage, incircle_coplanar, incircle_projected_staged, insphere_staged,
    orient2d_projected_staged, orient3d_staged,
};
use mesher_geom::{Point3, Sign};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, ToPrimitive, Zero};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

type Q = BigRational;
type V = [Q; 3];

fn cases() -> usize {
    if let Some(n) = std::env::var("IMPLICIT_CASES")
        .ok()
        .and_then(|s| s.parse().ok())
    {
        return n;
    }
    if cfg!(debug_assertions) {
        2_000
    } else {
        100_000
    }
}

/// Deterministic xorshift.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn point(&mut self, centre: Point3, size: f64) -> Point3 {
        [0, 1, 2].map(|k| centre[k] + self.range(-size, size))
    }
}

impl Rng {
    /// A nonzero state for chunk `chunk` of a test seeded with `seed`
    /// (SplitMix64).
    fn seeded(seed: u64, chunk: usize) -> Rng {
        let mut z = seed.wrapping_add((chunk as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15));
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        Rng((z ^ (z >> 31)) | 1)
    }
}

/// Run `cases()` cases of `f` on all cores in 64 fixed chunks, each with its
/// own seeded generator, so the cases do not depend on the core count. The
/// results come back in case order.
fn run<T: Send>(seed: u64, f: impl Fn(&mut Rng, usize) -> T + Sync) -> Vec<T> {
    const CHUNKS: usize = 64;
    let n = cases();
    let threads = std::thread::available_parallelism()
        .map_or(4, std::num::NonZero::get)
        .min(CHUNKS);
    let next = AtomicUsize::new(0);
    let done = Mutex::new(Vec::new());
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| {
                loop {
                    let chunk = next.fetch_add(1, Ordering::Relaxed);
                    if chunk >= CHUNKS {
                        break;
                    }
                    let mut rng = Rng::seeded(seed, chunk);
                    let out: Vec<T> = (n * chunk / CHUNKS..n * (chunk + 1) / CHUNKS)
                        .map(|i| f(&mut rng, i))
                        .collect();
                    done.lock().unwrap().push((chunk, out));
                }
            });
        }
    });
    let mut done = done.into_inner().unwrap();
    done.sort_by_key(|(c, _)| *c);
    done.into_iter().flat_map(|(_, o)| o).collect()
}

fn q(x: f64) -> Q {
    Q::from_float(x).unwrap()
}

fn v(p: Point3) -> V {
    p.map(q)
}

fn sub(a: &V, b: &V) -> V {
    [0, 1, 2].map(|k| &a[k] - &b[k])
}

fn add(a: &V, b: &V) -> V {
    [0, 1, 2].map(|k| &a[k] + &b[k])
}

fn scale(a: &V, s: &Q) -> V {
    [0, 1, 2].map(|k| &a[k] * s)
}

fn dot(a: &V, b: &V) -> Q {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

fn cross(a: &V, b: &V) -> V {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

fn sign(x: &Q) -> Sign {
    if x.is_positive() {
        Sign::Positive
    } else if x.is_negative() {
        Sign::Negative
    } else {
        Sign::Zero
    }
}

fn normal(p: &[Point3; 3]) -> (V, V) {
    let [a, b, c] = p.map(v);
    (cross(&sub(&b, &a), &sub(&c, &a)), a)
}

/// Solve a 3×3 system by Gaussian elimination with exact pivoting; `None`
/// when singular.
fn solve(mut m: [[Q; 3]; 3], mut r: [Q; 3]) -> Option<V> {
    for col in 0..3 {
        let pivot = (col..3).find(|&i| !m[i][col].is_zero())?;
        m.swap(col, pivot);
        r.swap(col, pivot);
        for i in 0..3 {
            if i != col && !m[i][col].is_zero() {
                let f = &m[i][col] / &m[col][col];
                for k in 0..3 {
                    let t = &f * &m[col][k];
                    m[i][k] -= t;
                }
                let t = &f * &r[col];
                r[i] -= t;
            }
        }
    }
    Some([0, 1, 2].map(|k| &r[k] / &m[k][k]))
}

/// The exact coordinates of a point, checked against its definition.
fn exact(p: &Point) -> V {
    match *p.kind() {
        Kind::Explicit(a) => v(a),
        Kind::Lnc { a, b, t } => {
            let t = q(t);
            let one_minus = Q::from_integer(BigInt::from(1)) - &t;
            add(&scale(&v(a), &t), &scale(&v(b), &one_minus))
        }
        Kind::Bpt { tri, weights } => {
            let w = weights.map(q);
            let total = &w[0] + &w[1] + &w[2];
            let s = add(
                &add(&scale(&v(tri[0]), &w[0]), &scale(&v(tri[1]), &w[1])),
                &scale(&v(tri[2]), &w[2]),
            );
            scale(&s, &(Q::from_integer(BigInt::from(1)) / total))
        }
        Kind::Lpi { line, plane } => {
            let (n, r) = normal(&plane);
            let (p0, p1) = (v(line[0]), v(line[1]));
            let d = sub(&p1, &p0);
            let lambda = dot(&n, &sub(&r, &p0)) / dot(&n, &d);
            let x = add(&p0, &scale(&d, &lambda));
            assert!(dot(&n, &sub(&x, &r)).is_zero(), "LPI off its plane");
            assert!(
                cross(&sub(&x, &p0), &d).iter().all(Zero::is_zero),
                "LPI off its line"
            );
            x
        }
        Kind::Tpi { planes } => {
            let rows = planes.map(|p| normal(&p));
            let m = [0, 1, 2].map(|i| rows[i].0.clone());
            let r = [0, 1, 2].map(|i| dot(&rows[i].0, &rows[i].1));
            let x = solve(m, r).expect("TPI planes are independent");
            for (n, a) in &rows {
                assert!(dot(n, &sub(&x, a)).is_zero(), "TPI off a plane");
            }
            x
        }
    }
}

fn ref_orient3d(p: [&Point; 4]) -> Sign {
    let [a, b, c, d] = p.map(exact);
    sign(&dot(&cross(&sub(&b, &a), &sub(&c, &a)), &sub(&d, &a)))
}

fn det(m: &[Vec<Q>]) -> Q {
    if m.len() == 1 {
        return m[0][0].clone();
    }
    let mut total = Q::zero();
    for j in 0..m.len() {
        let minor: Vec<Vec<Q>> = m[1..]
            .iter()
            .map(|r| {
                r.iter()
                    .enumerate()
                    .filter(|&(k, _)| k != j)
                    .map(|(_, x)| x.clone())
                    .collect()
            })
            .collect();
        let term = &m[0][j] * det(&minor);
        if j % 2 == 0 {
            total += term;
        } else {
            total -= term;
        }
    }
    total
}

/// The classic in-sphere determinant on the exact coordinates, with the
/// convention flipped as `predicates::insphere` does; with `check`, and a
/// non-flat `abcd`, also checked against the side of the exact circumcentre.
fn ref_insphere(p: [&Point; 5], check: bool) -> Sign {
    let x = p.map(exact);
    let e = &x[4];
    let rows: Vec<Vec<Q>> = x[..4]
        .iter()
        .map(|a| {
            let d = sub(a, e);
            vec![d[0].clone(), d[1].clone(), d[2].clone(), dot(&d, &d)]
        })
        .collect();
    let s = match sign(&det(&rows)) {
        Sign::Positive => Sign::Negative,
        Sign::Negative => Sign::Positive,
        Sign::Zero => Sign::Zero,
    };
    let orientation = if check {
        ref_orient3d([p[0], p[1], p[2], p[3]])
    } else {
        Sign::Zero
    };
    if check && orientation != Sign::Zero {
        // |c - a|² = |c - q|² for q in b, c, d: 2(q - a)·c = |q|² - |a|².
        let a = &x[0];
        let two = Q::from_integer(BigInt::from(2));
        let m = [1, 2, 3].map(|i| sub(&x[i], a).map(|c| &c * &two));
        let r = [1, 2, 3].map(|i| dot(&x[i], &x[i]) - dot(a, a));
        let centre = solve(m, r).unwrap();
        let radius = dot(&sub(a, &centre), &sub(a, &centre));
        let dist = dot(&sub(e, &centre), &sub(e, &centre));
        let inside = sign(&(radius - dist));
        let expected = match (inside, orientation) {
            (Sign::Zero, _) => Sign::Zero,
            (i, o) if i == o => Sign::Positive,
            _ => Sign::Negative,
        };
        assert_eq!(
            s, expected,
            "the reference determinant disagrees with the circumcentre"
        );
    }
    s
}

fn kept(axis: usize) -> [usize; 2] {
    [(axis + 1) % 3, (axis + 2) % 3]
}

fn ref_orient2d(p: [&Point; 3], axis: usize) -> Sign {
    let [a, b, c] = p.map(exact);
    sign(&cross(&sub(&b, &a), &sub(&c, &a))[axis])
}

fn ref_incircle_projected(p: [&Point; 4], axis: usize) -> Sign {
    let [i, j] = kept(axis);
    let x = p.map(exact);
    let d = &x[3];
    let rows: Vec<Vec<Q>> = x[..3]
        .iter()
        .map(|a| {
            let (u, w) = (&a[i] - &d[i], &a[j] - &d[j]);
            let lift = &u * &u + &w * &w;
            vec![u, w, lift]
        })
        .collect();
    sign(&det(&rows))
}

/// In the circumcircle of `abc`, for `d` in their plane: by the exact
/// circumcentre in that plane.
fn ref_incircle_coplanar(p: [&Point; 4]) -> Sign {
    let [a, b, c, d] = p.map(exact);
    let n = cross(&sub(&b, &a), &sub(&c, &a));
    assert!(dot(&n, &sub(&d, &a)).is_zero(), "d is not coplanar");
    if n.iter().all(Zero::is_zero) {
        return Sign::Zero;
    }
    let two = Q::from_integer(BigInt::from(2));
    let m = [
        sub(&b, &a).map(|x| &x * &two),
        sub(&c, &a).map(|x| &x * &two),
        n.clone(),
    ];
    let r = [
        dot(&b, &b) - dot(&a, &a),
        dot(&c, &c) - dot(&a, &a),
        dot(&n, &a),
    ];
    let centre = solve(m, r).unwrap();
    let radius = dot(&sub(&a, &centre), &sub(&a, &centre));
    let dist = dot(&sub(&d, &centre), &sub(&d, &centre));
    sign(&(radius - dist))
}

/// A random point of a random kind near `centre`.
fn random_point(rng: &mut Rng, centre: Point3, size: f64) -> Point {
    loop {
        let made = match rng.below(5) {
            0 => Point::explicit(rng.point(centre, size)),
            1 => Point::lnc(rng.point(centre, size), rng.point(centre, size), rng.unit()),
            2 => Point::bpt(
                [0; 3].map(|_| rng.point(centre, size)),
                [0; 3].map(|_| rng.range(0.05, 1.0)),
            ),
            3 => Point::lpi(
                [0; 2].map(|_| rng.point(centre, size)),
                [0; 3].map(|_| rng.point(centre, size)),
            ),
            _ => Point::tpi([0; 3].map(|_| [0; 3].map(|_| rng.point(centre, size)))),
        };
        if let Ok(p) = made {
            // Keep implicit points near the configuration: a line nearly
            // parallel to a plane meets it far away.
            let a = p.approx();
            if (0..3).all(|k| (a[k] - centre[k]).abs() <= 4.0 * size) {
                return p;
            }
        }
    }
}

/// A random centre: near the origin, or far from it.
fn centre(rng: &mut Rng) -> (Point3, f64) {
    let far = rng.below(3) == 0;
    let c = if far {
        rng.point([0.0; 3], 1e4)
    } else {
        [0.0; 3]
    };
    (c, if rng.below(2) == 0 { 1.0 } else { 1e-3 })
}

#[derive(Default)]
struct Tally {
    filter: usize,
    total: usize,
}

impl Tally {
    fn count(&mut self, stage: Stage) {
        if stage != Stage::Explicit {
            self.total += 1;
            if stage == Stage::Filter {
                self.filter += 1;
            }
        }
    }
    fn rate(&self) -> f64 {
        self.filter as f64 / self.total.max(1) as f64
    }
}

#[test]
fn random_queries_agree_and_the_filter_settles_them() {
    let stages = run(0x9e37_79b9_7f4a_7c15, |rng, i| {
        let (c, size) = centre(rng);
        let p: Vec<Point> = (0..5).map(|_| random_point(rng, c, size)).collect();
        let axis = rng.below(3);
        let (s, st0) = orient3d_staged(&p[0], &p[1], &p[2], &p[3]);
        assert_eq!(
            s,
            ref_orient3d([&p[0], &p[1], &p[2], &p[3]]),
            "orient3d {p:?}"
        );
        let (s, st1) = insphere_staged(&p[0], &p[1], &p[2], &p[3], &p[4]);
        assert_eq!(
            s,
            ref_insphere([&p[0], &p[1], &p[2], &p[3], &p[4]], i % 8 == 0),
            "insphere {p:?}"
        );
        let (s, st2) = orient2d_projected_staged(&p[0], &p[1], &p[2], axis);
        assert_eq!(
            s,
            ref_orient2d([&p[0], &p[1], &p[2]], axis),
            "orient2d {p:?}"
        );
        let (s, st3) = incircle_projected_staged(&p[0], &p[1], &p[2], &p[3], axis);
        assert_eq!(
            s,
            ref_incircle_projected([&p[0], &p[1], &p[2], &p[3]], axis),
            "incircle {p:?}"
        );
        [st0, st1, st2, st3]
    });
    let [mut o3, mut is, mut o2, mut ic] = [0; 4].map(|_| Tally::default());
    for [a, b, c, d] in stages {
        o3.count(a);
        is.count(b);
        o2.count(c);
        ic.count(d);
    }
    eprintln!(
        "filter settled: orient3d {:.4}, insphere {:.4}, orient2d {:.4}, incircle {:.4}",
        o3.rate(),
        is.rate(),
        o2.rate(),
        ic.rate()
    );
    for (name, t) in [
        ("orient3d", &o3),
        ("insphere", &is),
        ("orient2d", &o2),
        ("incircle", &ic),
    ] {
        assert!(
            t.rate() >= 0.95,
            "{name}: the filter settled only {:.4}",
            t.rate()
        );
    }
}

/// The same point in several implicit forms, all exactly equal to an
/// explicit `x` on the triangle `abc` (or anywhere, for LNC and TPI).
fn forms_of(rng: &mut Rng, x: Point3, size: f64) -> Vec<Point> {
    let mut out = vec![Point::explicit(x).unwrap()];
    let other = rng.point(x, size);
    // t = 1 gives a exactly; t = 0 gives b.
    out.push(Point::lnc(x, other, 1.0).unwrap());
    out.push(Point::lnc(other, x, 0.0).unwrap());
    out.push(Point::bpt([x, other, rng.point(x, size)], [2.5, 0.0, 0.0]).unwrap());
    // A line through x meeting a plane through x.
    if let Ok(p) = Point::lpi(
        [x, rng.point(x, size)],
        [x, rng.point(x, size), rng.point(x, size)],
    ) {
        out.push(p);
    }
    // Three planes through x.
    if let Ok(p) = Point::tpi([0; 3].map(|_| [x, rng.point(x, size), rng.point(x, size)])) {
        out.push(p);
    }
    out
}

#[test]
fn degenerate_queries_agree() {
    let kinds = run(0x2545_f491_4f6c_dd1d, |rng, i| {
        let (c, size) = centre(rng);
        match i % 4 {
            // Coplanar: the fourth point lies exactly in plane abc.
            0 => {
                let tri = [0; 3].map(|_| rng.point(c, size));
                let d = match rng.below(4) {
                    0 => Point::bpt(tri, [0; 3].map(|_| rng.range(-0.5, 1.5))),
                    1 => Point::lnc(tri[0], tri[1], rng.range(-1.0, 2.0)),
                    2 => Point::lpi([rng.point(c, size), rng.point(c, size)], tri),
                    _ => Point::tpi([
                        tri,
                        [0; 3].map(|_| rng.point(c, size)),
                        [0; 3].map(|_| rng.point(c, size)),
                    ]),
                };
                let Ok(d) = d else { return None };
                // Use implicit forms of the triangle's corners too.
                let a = forms_of(rng, tri[0], size);
                let a = &a[rng.below(a.len())];
                let (b, cc) = (Point::from(tri[1]), Point::from(tri[2]));
                let s = implicit::orient3d(a, &b, &cc, &d);
                assert_eq!(s, ref_orient3d([a, &b, &cc, &d]));
                assert_eq!(s, Sign::Zero, "coplanar {a:?} {d:?}");
                let axis = rng.below(3);
                let s = implicit::orient2d_projected(a, &b, &cc, axis);
                assert_eq!(s, ref_orient2d([a, &b, &cc], axis));
                // In-circle in the plane, with d coplanar.
                let s = incircle_coplanar(a, &b, &cc, &d);
                assert_eq!(
                    s,
                    ref_incircle_coplanar([a, &b, &cc, &d]),
                    "incircle_coplanar"
                );
                let s = implicit::incircle_projected(a, &b, &cc, &d, axis);
                assert_eq!(s, ref_incircle_projected([a, &b, &cc, &d], axis));
                Some(0)
            }
            // Cospherical: points of an integer sphere, scaled by a power of
            // two and moved by an exact offset, in implicit forms.
            1 => {
                const SPHERE: [[f64; 3]; 10] = [
                    [3.0, 4.0, 0.0],
                    [0.0, 5.0, 0.0],
                    [5.0, 0.0, 0.0],
                    [0.0, 0.0, 5.0],
                    [-3.0, -4.0, 0.0],
                    [0.0, 3.0, -4.0],
                    [4.0, 0.0, 3.0],
                    [-5.0, 0.0, 0.0],
                    [0.0, -4.0, 3.0],
                    [-4.0, 0.0, -3.0],
                ];
                let s = [0.25, 1.0, 64.0][rng.below(3)];
                let off = [0, 1, 2].map(|k| (c[k] * 8.0).round() / 8.0);
                let pts: Vec<Point> = (0..5)
                    .map(|_| {
                        let x = SPHERE[rng.below(SPHERE.len())].map(|t| t * s);
                        let x = [0, 1, 2].map(|k| x[k] + off[k]);
                        let forms = forms_of(rng, x, s);
                        forms[rng.below(forms.len())]
                    })
                    .collect();
                let got = implicit::insphere(&pts[0], &pts[1], &pts[2], &pts[3], &pts[4]);
                assert_eq!(
                    got,
                    ref_insphere([&pts[0], &pts[1], &pts[2], &pts[3], &pts[4]], true)
                );
                assert_eq!(got, Sign::Zero, "cospherical");
                Some(1)
            }
            // Coincident: an implicit point equal to another point of the query.
            2 => {
                let x = rng.point(c, size);
                let forms = forms_of(rng, x, size);
                let a = &forms[rng.below(forms.len())];
                let b = &forms[rng.below(forms.len())];
                let p = [random_point(rng, c, size), random_point(rng, c, size)];
                let s = implicit::orient3d(a, &p[0], b, &p[1]);
                assert_eq!(s, ref_orient3d([a, &p[0], b, &p[1]]));
                assert_eq!(s, Sign::Zero, "coincident");
                let e = random_point(rng, c, size);
                let s = implicit::insphere(a, &p[0], &p[1], &e, b);
                assert_eq!(s, ref_insphere([a, &p[0], &p[1], &e, b], i % 8 == 2));
                Some(2)
            }
            // On a segment: an LNC on ab with a, b and another point.
            _ => {
                let a = rng.point(c, size);
                let b = rng.point(c, size);
                let m = Point::lnc(a, b, rng.unit()).unwrap();
                let axis = rng.below(3);
                let pa = Point::from(a);
                let pb = Point::from(b);
                let s = implicit::orient2d_projected(&pa, &m, &pb, axis);
                assert_eq!(s, ref_orient2d([&pa, &m, &pb], axis));
                assert_eq!(s, Sign::Zero, "collinear");
                let d = random_point(rng, c, size);
                let s = implicit::orient3d(&pa, &m, &pb, &d);
                assert_eq!(s, Sign::Zero);
                Some(3)
            }
        }
    });
    let mut zeros = [0usize; 4];
    for k in kinds.into_iter().flatten() {
        zeros[k] += 1;
    }
    eprintln!("degenerate cases: {zeros:?}");
}

#[test]
fn coincidence_and_collinearity_agree() {
    let kinds = run(0x510e_527f_ade6_82d1, |rng, i| {
        let (c, size) = centre(rng);
        let x = rng.point(c, size);
        let forms = forms_of(rng, x, size);
        let a = forms[rng.below(forms.len())];
        let b = if i % 2 == 0 {
            forms[rng.below(forms.len())]
        } else {
            random_point(rng, c, size)
        };
        let same = exact(&a) == exact(&b);
        assert_eq!(implicit::coincident(&a, &b), same, "{a:?} {b:?}");
        let m = random_point(rng, c, size);
        let on = Point::lnc(x, m.approx(), rng.unit()).unwrap();
        let line = [a, Point::from(m.approx()), on];
        let (pa, pb, pc) = (exact(&line[0]), exact(&line[1]), exact(&line[2]));
        let straight = cross(&sub(&pb, &pa), &sub(&pc, &pa))
            .iter()
            .all(Zero::is_zero);
        assert_eq!(implicit::collinear(&line[0], &line[1], &line[2]), straight);
        same
    });
    eprintln!("coincident pairs: {}", kinds.iter().filter(|&&s| s).count());
}

#[test]
fn nearly_degenerate_queries_agree() {
    let runs = run(0x6a09_e667_f3bc_c909, |rng, i| {
        let mut seen: Vec<(&'static str, Stage)> = Vec::new();
        let (c, size) = centre(rng);
        let tilt = 10f64.powi(-(rng.below(14) as i32) - 2);
        // Three nearly parallel planes: a base plane and two tilted copies.
        let base = [0; 3].map(|_| rng.point(c, size));
        let tilted = |rng: &mut Rng| -> [Point3; 3] {
            let mut p = base;
            for q in &mut p {
                for k in 0..3 {
                    q[k] += tilt * size * rng.range(-1.0, 1.0);
                }
            }
            p
        };
        let planes = [base, tilted(rng), tilted(rng)];
        let Ok(t) = Point::tpi(planes) else {
            return seen;
        };
        let others: Vec<Point> = (0..4).map(|_| random_point(rng, c, size)).collect();
        let (s, st) = orient3d_staged(&t, &others[0], &others[1], &others[2]);
        assert_eq!(s, ref_orient3d([&t, &others[0], &others[1], &others[2]]));
        seen.push(("orient3d", st));
        // A point nearly in the plane of three others: a line nearly
        // parallel to plane abc, meeting the nearly parallel plane.
        let tri = [0; 3].map(|_| rng.point(c, size));
        let near = Point::lpi(
            [rng.point(c, size), rng.point(c, size)],
            [
                tri[0],
                tri[1],
                [0, 1, 2].map(|k| tri[2][k] + tilt * size * rng.range(-1.0, 1.0)),
            ],
        );
        if let Ok(near) = near {
            let (a, b, cc) = (
                Point::from(tri[0]),
                Point::from(tri[1]),
                Point::from(tri[2]),
            );
            let (s, st) = orient3d_staged(&a, &b, &cc, &near);
            assert_eq!(s, ref_orient3d([&a, &b, &cc, &near]));
            seen.push(("near plane", st));
            if i % 4 == 0 {
                let (s, _) = insphere_staged(&t, &a, &b, &cc, &near);
                assert_eq!(s, ref_insphere([&t, &a, &b, &cc, &near], i % 32 == 0));
            }
        }
        seen
    });
    let mut stages = std::collections::BTreeMap::new();
    for key in runs.into_iter().flatten() {
        *stages.entry(key).or_insert(0usize) += 1;
    }
    eprintln!("nearly degenerate, settled by stage: {stages:?}");
}

#[test]
fn constructors_reject_what_does_not_exist() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    assert!(Point::explicit([f64::NAN, 0.0, 0.0]).is_err());
    assert!(Point::lnc(a, b, f64::INFINITY).is_err());
    assert_eq!(
        Point::bpt([a, b, c], [1.0, -0.5, -0.5]),
        Err(implicit::ImplicitError::ZeroWeight)
    );
    // A line parallel to the plane z = 0, and one lying in it.
    assert_eq!(
        Point::lpi([[0.0, 0.0, 1.0], [1.0, 1.0, 1.0]], [a, b, c]),
        Err(implicit::ImplicitError::NoLinePlaneIntersection)
    );
    assert!(Point::lpi([a, b], [a, b, c]).is_err());
    // A degenerate plane.
    assert!(Point::lpi([[0.0, 0.0, 1.0], a], [a, b, [2.0, 0.0, 0.0]]).is_err());
    // Two parallel planes, and three planes through one line.
    let z1 = [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]];
    let x0 = [a, [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(
        Point::tpi([[a, b, c], z1, x0]),
        Err(implicit::ImplicitError::NoThreePlaneIntersection)
    );
    let y0 = [a, b, [0.0, 0.0, 1.0]];
    let diagonal = [a, b, [0.0, 1.0, 1.0]];
    assert!(Point::tpi([[a, b, c], y0, diagonal]).is_err());
    // The corner of the unit cube's three coordinate planes.
    let corner = Point::tpi([[a, b, c], y0, x0]).unwrap();
    assert_eq!(corner.approx(), [0.0, 0.0, 0.0]);
    // A negative W is fine: the predicates account for its sign.
    let p = Point::bpt([a, b, c], [-1.0, -1.0, -1.0]).unwrap();
    let up = Point::from([0.0, 0.0, 1.0]);
    assert_eq!(
        implicit::orient3d(&Point::from(a), &Point::from(b), &Point::from(c), &up),
        Sign::Positive
    );
    assert_eq!(
        implicit::orient3d(&p, &Point::from(b), &Point::from(c), &up),
        Sign::Positive
    );
}

#[test]
fn approximate_positions_are_close() {
    let mut rng = Rng(0x3c6e_f372_fe94_f82b);
    for _ in 0..1000 {
        let (c, size) = centre(&mut rng);
        let p = random_point(&mut rng, c, size);
        let x = exact(&p);
        let a = p.approx();
        for k in 0..3 {
            let exact_k = x[k].to_f64().unwrap();
            assert!(
                (a[k] - exact_k).abs() <= 1e-9 * (size + exact_k.abs()),
                "{p:?}: {a:?}"
            );
        }
    }
}

/// Milestone 14's angle predicates against the rational reference: random
/// queries (with the filter rate), and exactly right angles and right
/// dihedrals built on an integer lattice, with the points in implicit forms.
#[test]
fn angle_queries_agree() {
    let ref_dot = |p: &Point, q: &Point, r: &Point| {
        let (p, q, r) = (exact(p), exact(q), exact(r));
        sign(&dot(&sub(&p, &q), &sub(&r, &q)))
    };
    let ref_dihedral = |a: &Point, b: &Point, c: &Point, d: &Point| {
        let (a, b, c, d) = (exact(a), exact(b), exact(c), exact(d));
        let e = sub(&c, &d);
        sign(&dot(&cross(&sub(&b, &d), &e), &cross(&sub(&a, &d), &e)))
    };
    let stages = run(0x3c6e_f372_fe94_f82b, |rng, _| {
        let (c, size) = centre(rng);
        let p: Vec<Point> = (0..4).map(|_| random_point(rng, c, size)).collect();
        let (s, st0) = implicit::dot_sign_staged(&p[0], &p[1], &p[2]);
        assert_eq!(s, ref_dot(&p[0], &p[1], &p[2]), "dot {p:?}");
        let (s, st1) = implicit::dihedral_sign_staged(&p[0], &p[1], &p[2], &p[3]);
        assert_eq!(
            s,
            ref_dihedral(&p[0], &p[1], &p[2], &p[3]),
            "dihedral {p:?}"
        );
        [st0, st1]
    });
    let [mut dt, mut dh] = [0; 2].map(|_| Tally::default());
    for [a, b] in stages {
        dt.count(a);
        dh.count(b);
    }
    eprintln!(
        "filter settled: dot {:.4}, dihedral {:.4}",
        dt.rate(),
        dh.rate()
    );
    for (name, t) in [("dot", &dt), ("dihedral", &dh)] {
        assert!(
            t.rate() >= 0.95,
            "{name}: the filter settled only {:.4}",
            t.rate()
        );
    }

    // Exactly right: integer vectors, scaled by a power of two and moved by
    // an offset on an eighth grid, so every coordinate is exact in f64.
    let int3 = |rng: &mut Rng| [0; 3].map(|_| (rng.below(9) as f64) - 4.0);
    let icross = |a: Point3, b: Point3| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let nonzero = |a: Point3| a.iter().any(|&x| x != 0.0);
    let zeros = run(0x7137_4491_23ef_65cd, |rng, i| {
        let (c, _) = centre(rng);
        let s = [0.25, 1.0, 64.0][rng.below(3)];
        let o = [0, 1, 2].map(|k| (c[k] * 8.0).round() / 8.0);
        let at = |rng: &mut Rng, u: Point3| {
            let x = [0, 1, 2].map(|k| o[k] + s * u[k]);
            let forms = forms_of(rng, x, s);
            forms[rng.below(forms.len())]
        };
        if i % 2 == 0 {
            // (p − q) ⊥ (r − q): r − q = (p − q) × w.
            let u = int3(rng);
            let v = icross(u, int3(rng));
            if !nonzero(u) || !nonzero(v) {
                return None;
            }
            let (p, q, r) = (at(rng, u), at(rng, [0.0; 3]), at(rng, v));
            let got = implicit::dot_sign(&p, &q, &r);
            assert_eq!(got, ref_dot(&p, &q, &r));
            assert_eq!(got, Sign::Zero, "right angle");
            Some(0)
        } else {
            // Faces a c d and b c d at a right dihedral: a − d and b − d
            // perpendicular to each other and to the edge c − d.
            let e = int3(rng);
            let u = icross(e, int3(rng));
            let v = icross(e, u);
            if !nonzero(e) || !nonzero(u) || !nonzero(v) {
                return None;
            }
            let (a, b, cc, d) = (at(rng, u), at(rng, v), at(rng, e), at(rng, [0.0; 3]));
            let got = implicit::dihedral_sign(&a, &b, &cc, &d);
            assert_eq!(got, ref_dihedral(&a, &b, &cc, &d));
            assert_eq!(got, Sign::Zero, "right dihedral");
            Some(1)
        }
    });
    let mut count = [0usize; 2];
    for k in zeros.into_iter().flatten() {
        count[k] += 1;
    }
    eprintln!("right angles {}, right dihedrals {}", count[0], count[1]);
    assert!(count[0] > 0 && count[1] > 0);
}

/// The projection test on a flat triangle: an edge leaning over the
/// triangle projects into it; one leaning away, or over a side only, does
/// not.
#[test]
fn edge_projection_into_a_triangle() {
    let p = |x: [f64; 3]| Point::explicit(x).unwrap();
    let (a, b, c) = (p([0.0; 3]), p([4.0, 0.0, 0.0]), p([0.0, 4.0, 0.0]));
    assert!(implicit::edge_projects_into_triangle(
        &a,
        &p([1.0, 1.0, 0.5]),
        &b,
        &c
    ));
    assert!(!implicit::edge_projects_into_triangle(
        &a,
        &p([-1.0, -1.0, 0.5]),
        &b,
        &c
    ));
    assert!(!implicit::edge_projects_into_triangle(
        &a,
        &p([1.0, -1.0, 0.5]),
        &b,
        &c
    ));
}
