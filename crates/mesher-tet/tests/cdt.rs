//! The exact Delaunay core and segment recovery of the exact constrained
//! kernel (`mesher_tet::cdt`, Milestone 12).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_geom::Point3;
use mesher_geom::implicit::Point;
use mesher_tet::cdt::{ExactDelaunay, Inserted, recover_segments};
use mesher_tet::delaunay;

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
    fn point(&mut self) -> Point3 {
        [self.unit(), self.unit(), self.unit()]
    }
}

fn explicit(p: &[Point3]) -> Vec<Point> {
    p.iter().map(|&q| Point::explicit(q).unwrap()).collect()
}

#[test]
fn random_points_match_the_float_kernel() {
    for seed in 1..=5u64 {
        let mut rng = Rng(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15));
        let p: Vec<Point3> = (0..300).map(|_| rng.point()).collect();
        let exact = ExactDelaunay::new(explicit(&p), 7).unwrap();
        exact.check().unwrap();
        exact.check_delaunay().unwrap();
        // In general position the Delaunay tetrahedralization is unique.
        let float = delaunay(&p).unwrap();
        let key = |mut t: [u32; 4]| {
            t.sort_unstable();
            t
        };
        let mut a: Vec<[u32; 4]> = exact.tets().into_iter().map(key).collect();
        let mut b: Vec<[u32; 4]> = float.tets.into_iter().map(key).collect();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b, "seed {seed}");
    }
}

#[test]
fn a_cospherical_grid_is_valid_and_delaunay() {
    let mut p = Vec::new();
    for x in 0..5 {
        for y in 0..5 {
            for z in 0..5 {
                p.push([f64::from(x), f64::from(y), f64::from(z)]);
            }
        }
    }
    let d = ExactDelaunay::new(explicit(&p), 3).unwrap();
    d.check().unwrap();
    d.check_delaunay().unwrap();
    // The tets fill the cube exactly.
    let volume: f64 = d
        .tets()
        .iter()
        .map(|t| mesher_tet::quality::float_volume(t.map(|v| p[v as usize])))
        .sum();
    assert!((volume - 64.0).abs() < 1e-9, "{volume}");
    assert!(d.duplicates().is_empty());
}

#[test]
fn implicit_points_are_inserted_exactly() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let corners: Vec<Point3> = (0..8)
        .map(|i| {
            [
                f64::from(i & 1),
                f64::from((i >> 1) & 1),
                f64::from((i >> 2) & 1),
            ]
        })
        .collect();
    let mut pts = explicit(&corners);
    for k in 0..200 {
        let p = match k % 4 {
            // On an edge of the cube, exactly.
            0 => Point::lnc(corners[0], corners[1 << (k % 3)], rng.unit()).unwrap(),
            // On a face, exactly.
            1 => Point::bpt(
                [corners[0], corners[1], corners[3]],
                [rng.unit() + 0.01, rng.unit() + 0.01, rng.unit() + 0.01],
            )
            .unwrap(),
            // Three planes meeting inside.
            2 => {
                let q = rng.point();
                // Three planes through q, each spanned by two axis offsets.
                let mut plane = |k: usize| {
                    let (mut b, mut c) = (q, q);
                    b[(k + 1) % 3] += 0.3 + rng.unit();
                    c[(k + 2) % 3] += 0.3 + rng.unit();
                    [q, b, c]
                };
                Point::tpi([plane(0), plane(1), plane(2)]).unwrap()
            }
            _ => Point::explicit(rng.point()).unwrap(),
        };
        pts.push(p);
    }
    let mut d = ExactDelaunay::new(pts, 11).unwrap();
    d.check().unwrap();
    d.check_delaunay().unwrap();
    // A point equal to a corner, built implicitly, is a duplicate.
    let again = Point::lnc(corners[7], corners[0], 1.0).unwrap();
    assert!(matches!(d.insert(again).unwrap(), Inserted::Duplicate(_)));
    let new = d.insert(Point::explicit([0.5, 0.5, 2.0]).unwrap()).unwrap();
    assert!(matches!(new, Inserted::New(_)));
    d.check().unwrap();
    d.check_delaunay().unwrap();
}

type Surface = (Vec<Point3>, Vec<[u32; 3]>);

fn cube() -> Surface {
    let p: Vec<Point3> = (0..8)
        .map(|i| {
            [
                f64::from(i & 1),
                f64::from((i >> 1) & 1),
                f64::from((i >> 2) & 1),
            ]
        })
        .collect();
    let quads = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let f = quads
        .iter()
        .flat_map(|&[a, b, c, d]| [[a, b, c], [a, c, d]])
        .collect();
    (p, f)
}

fn uv_sphere(segments: u32, rings: u32, centre: Point3, round: bool) -> Surface {
    let r = 0.8;
    let snap = |x: f64| if round { (x * 1e6).round() / 1e6 } else { x };
    let at = |q: Point3| [0, 1, 2].map(|k| snap(q[k] + centre[k]));
    let mut p = vec![at([0.0, 0.0, r])];
    for i in 1..rings {
        let t = f64::from(i) / f64::from(rings) * std::f64::consts::PI;
        for j in 0..segments {
            let a = f64::from(j) / f64::from(segments) * std::f64::consts::TAU;
            p.push(at([
                r * libm::sin(t) * libm::cos(a),
                r * libm::sin(t) * libm::sin(a),
                r * libm::cos(t),
            ]));
        }
    }
    let south = p.len() as u32;
    p.push(at([0.0, 0.0, -r]));
    let ring = |i: u32, j: u32| 1 + i * segments + j % segments;
    let mut f = Vec::new();
    for j in 0..segments {
        f.push([0, ring(0, j), ring(0, j + 1)]);
    }
    for i in 0..rings - 2 {
        for j in 0..segments {
            let (a, b) = (ring(i, j), ring(i, j + 1));
            let (c, d) = (ring(i + 1, j), ring(i + 1, j + 1));
            f.push([a, c, d]);
            f.push([a, d, b]);
        }
    }
    for j in 0..segments {
        f.push([south, ring(rings - 2, j + 1), ring(rings - 2, j)]);
    }
    (p, f)
}

/// Two boxes sharing an edge's line at a small angle: segments meeting at a
/// vertex at about 3°, the case concentric shells exist for.
fn sharp_wedge() -> Surface {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.05, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 0.05, 1.0],
    ];
    let f = vec![
        [0, 2, 1],
        [3, 4, 5],
        [0, 1, 4],
        [0, 4, 3],
        [1, 2, 5],
        [1, 5, 4],
        [2, 0, 3],
        [2, 3, 5],
    ];
    (p, f)
}

/// Every subsegment is an edge, and each input segment's subsegments chain
/// from one end to the other.
fn check_recovered(name: &str, (p, f): &Surface) -> usize {
    let r = recover_segments(p, f, 0, 5).unwrap_or_else(|e| panic!("{name}: {e}"));
    r.dt.check().unwrap_or_else(|e| panic!("{name}: {e}"));
    for s in &r.subsegments {
        assert!(r.dt.has_edge(s.ends[0], s.ends[1]), "{name}: {s:?} missing");
    }
    for seg in mesher_tet::cdt::segments::segments_of(f) {
        let chain: Vec<_> = r.subsegments.iter().filter(|s| s.segment == seg).collect();
        assert_eq!(chain[0].ends[0], seg.0, "{name}");
        assert_eq!(chain[chain.len() - 1].ends[1], seg.1, "{name}");
        for w in chain.windows(2) {
            assert_eq!(w[0].ends[1], w[1].ends[0], "{name}");
            assert_eq!(w[0].t[1], w[1].t[0], "{name}");
        }
    }
    eprintln!(
        "{name:<20} {:>5} triangles, {:>5} subsegments, {:>5} Steiner points",
        f.len(),
        r.subsegments.len(),
        r.steiner
    );
    r.steiner
}

#[test]
fn segments_are_recovered_exactly() {
    // The cube's corners are cospherical: the Delaunay tetrahedralization may
    // choose the other diagonal of a face, which is then split (a flip would
    // recover it without a point; not done yet).
    assert!(check_recovered("cube", &cube()) <= 6);
    check_recovered("uv sphere", &uv_sphere(16, 10, [0.0; 3], false));
    // The rounded, off-centre sphere the current kernel cannot recover.
    check_recovered(
        "rounded uv sphere",
        &uv_sphere(16, 10, [1.0, 1.0, 1.5], true),
    );
    check_recovered("sharp wedge", &sharp_wedge());
    check_recovered("fine uv sphere", &uv_sphere(48, 24, [0.3, -0.2, 0.1], true));
}
