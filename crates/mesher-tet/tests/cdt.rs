//! The exact Delaunay core and segment recovery of the exact constrained
//! kernel (`mesher_tet::cdt`, Milestone 12).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_geom::Point3;
use mesher_geom::implicit::Point;
use mesher_tet::cdt::{
    ExactDelaunay, Inserted, mesh_region, recover_faces, recover_segments, round_region,
};
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
    // choose the other diagonal of a face. Hull edge removal recovers them,
    // some only after removing an edge around them first.
    assert_eq!(check_recovered("cube", &cube()), 0);
    // Flips recover every sphere segment without a point.
    assert_eq!(
        check_recovered("uv sphere", &uv_sphere(16, 10, [0.0; 3], false)),
        0
    );
    // The rounded, off-centre sphere the current kernel cannot recover.
    assert_eq!(
        check_recovered(
            "rounded uv sphere",
            &uv_sphere(16, 10, [1.0, 1.0, 1.5], true),
        ),
        0
    );
    assert!(check_recovered("sharp wedge", &sharp_wedge()) <= 1);
    assert_eq!(
        check_recovered("fine uv sphere", &uv_sphere(48, 24, [0.3, -0.2, 0.1], true)),
        0
    );
}

#[test]
fn protected_insertion_keeps_protected_edges() {
    let mut rng = Rng(0x5851_f42d_4c95_7f2d);
    let p: Vec<Point3> = (0..200).map(|_| rng.point()).collect();
    let mut d = ExactDelaunay::new(explicit(&p), 9).unwrap();
    // Protect every edge at the vertex nearest the centre, then insert points
    // all around it, inside the hull.
    let hub = (0..p.len())
        .min_by(|&i, &j| {
            let r = |q: Point3| (0..3).map(|k| (q[k] - 0.5).powi(2)).sum::<f64>();
            r(p[i]).total_cmp(&r(p[j]))
        })
        .unwrap() as u32;
    let protected: Vec<u32> = d.neighbours(hub);
    assert!(!protected.is_empty());
    let is = |x: u32, y: u32| {
        (x == hub && protected.contains(&y)) || (y == hub && protected.contains(&x))
    };
    for _ in 0..100 {
        let q = [0, 1, 2].map(|k| p[hub as usize][k] + 0.2 * (rng.unit() - 0.5));
        d.insert_protected(Point::explicit(q).unwrap(), &is)
            .unwrap();
        d.check().unwrap();
    }
    for &v in &protected {
        assert!(d.has_edge(hub, v), "edge {hub}-{v} lost");
    }
    // A point outside the hull is refused.
    let outside = d.insert_protected(Point::explicit([2.0, 2.0, 2.0]).unwrap(), &is);
    assert!(outside.is_err());
}

#[test]
fn edge_removal_keeps_the_tetrahedralization_valid() {
    let never = |_: u32, _: u32| false;
    // Interior edges of a random tetrahedralization.
    let mut rng = Rng(0x1405_7b7e_f767_814f);
    let p: Vec<Point3> = (0..200).map(|_| rng.point()).collect();
    let mut d = ExactDelaunay::new(explicit(&p), 13).unwrap();
    let mut removed = 0;
    for _ in 0..300 {
        let x = (rng.next() % 200) as u32;
        let n = d.neighbours(x);
        let y = n[(rng.next() % n.len() as u64) as usize];
        if d.remove_edge(x, y, &never) {
            assert!(!d.has_edge(x, y));
            removed += 1;
        }
        d.check().unwrap();
    }
    assert!(removed > 50, "{removed} of 300 removed");
    // A grid's hull: every face diagonal lies in a hull plane, and removing
    // it flips the square's other diagonal in. Many points around it are
    // coplanar, so some need longer flip sequences than removal tries.
    let mut g = Vec::new();
    for x in 0..3 {
        for y in 0..3 {
            for z in 0..3 {
                g.push([f64::from(x), f64::from(y), f64::from(z)]);
            }
        }
    }
    let mut d = ExactDelaunay::new(explicit(&g), 3).unwrap();
    let at = |x: u32, y: u32, z: u32| 9 * x + 3 * y + z;
    let mut flipped = 0;
    for a in 0..2 {
        for b in 0..2 {
            for (u, v) in [
                (at(0, a, b), at(0, a + 1, b + 1)),
                (at(0, a + 1, b), at(0, a, b + 1)),
            ] {
                if d.has_edge(u, v) {
                    let removed = d.remove_edge(u, v, &never);
                    assert_eq!(removed, !d.has_edge(u, v));
                    flipped += usize::from(removed);
                    d.check().unwrap();
                }
            }
        }
    }
    assert!(flipped >= 1, "{flipped} of 4 removed");
}

/// Every subface is a face, and each triangle's subfaces cover it: their
/// areas sum to its area.
fn check_faces(name: &str, (p, f): &Surface) -> usize {
    let r = recover_faces(p, f, 0, 5).unwrap_or_else(|e| panic!("{name}: {e}"));
    r.dt.check().unwrap_or_else(|e| panic!("{name}: {e}"));
    let area = |t: [Point3; 3]| {
        let u = [0, 1, 2].map(|k| t[1][k] - t[0][k]);
        let w = [0, 1, 2].map(|k| t[2][k] - t[0][k]);
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
    };
    for (ti, t) in f.iter().enumerate() {
        let whole = area(t.map(|v| p[v as usize]));
        let parts: f64 = r
            .subfaces
            .iter()
            .filter(|s| s.triangle == ti as u32)
            .map(|s| area(s.v.map(|v| r.dt.points()[v as usize].approx())))
            .sum();
        assert!(
            (parts - whole).abs() <= 1e-9 * whole.max(1.0),
            "{name}: triangle {ti}"
        );
    }
    eprintln!(
        "{name:<20} {:>5} subfaces, {:>4} missing after segments, {:>4} triangles crossed, {} + {} Steiner points",
        r.subfaces.len(),
        r.initially_missing,
        r.missing.len(),
        r.steiner,
        r.steiner_on_faces
    );
    r.missing.len()
}

/// Schönhardt's twisted prism: a triangle over a triangle turned by 30°,
/// each side quad split along its reflex diagonal. Its inside cannot be
/// tetrahedralized without Steiner points.
fn schonhardt() -> Surface {
    let mut p = Vec::new();
    for (z, turn) in [(0.0, 0.0), (1.0, std::f64::consts::FRAC_PI_6)] {
        for k in 0..3 {
            let a = turn + f64::from(k) * std::f64::consts::TAU / 3.0;
            p.push([libm::cos(a), libm::sin(a), z]);
        }
    }
    let mut f = vec![[0, 2, 1], [3, 4, 5]];
    for k in 0..3u32 {
        let (a, b) = (k, (k + 1) % 3);
        let (c, d) = (a + 3, b + 3);
        // The quad a b d c, split along a-d.
        f.push([a, b, d]);
        f.push([a, d, c]);
    }
    (p, f)
}

/// A large tetrahedron with a small one just above the middle of its base
/// (a void inside it, faces inward) and one just below (outside), a tenth
/// of their size from it: the Delaunay edges between the small ones cross
/// the base.
fn sandwich() -> Surface {
    let tet = |c: Point3, r: f64, flip: bool| -> Vec<Point3> {
        let z = if flip { -r } else { r };
        vec![
            [c[0] - r, c[1] - r, c[2]],
            [c[0] + r, c[1] - r, c[2]],
            [c[0], c[1] + r, c[2]],
            [c[0], c[1], c[2] + z],
        ]
    };
    let mut p = tet([0.0, 0.0, 0.0], 4.0, false);
    p.extend(tet([0.1, 0.05, 0.02], 0.2, false));
    p.extend(tet([-0.05, 0.1, -0.02], 0.2, true));
    let mut f = Vec::new();
    // Reversed: the void, and the mirrored one below to face outward.
    for (base, flip) in [(0u32, false), (4, true), (8, true)] {
        let [a, b, c, d] = [base, base + 1, base + 2, base + 3];
        let faces = [[a, c, b], [a, b, d], [b, c, d], [c, a, d]];
        for t in faces {
            f.push(if flip { [t[0], t[2], t[1]] } else { t });
        }
    }
    (p, f)
}

#[test]
fn faces_are_recovered_by_flips() {
    assert_eq!(check_faces("cube", &cube()), 0);
    assert_eq!(
        check_faces("uv sphere", &uv_sphere(16, 10, [0.0; 3], false)),
        0
    );
    assert_eq!(
        check_faces(
            "rounded uv sphere",
            &uv_sphere(16, 10, [1.0, 1.0, 1.5], true)
        ),
        0
    );
    assert_eq!(check_faces("sharp wedge", &sharp_wedge()), 0);
    assert_eq!(
        check_faces("fine uv sphere", &uv_sphere(48, 24, [0.3, -0.2, 0.1], true)),
        0
    );
}

#[test]
fn faces_crossed_by_edges_are_recovered() {
    let (p, f) = sandwich();
    let r = recover_faces(&p, &f, 0, 5).unwrap();
    assert!(r.initially_missing > 0, "the base should start crossed");
    assert_eq!(check_faces("sandwich", &(p, f)), 0);
}

#[test]
fn an_untetrahedralizable_prism_is_recovered_with_steiner_points() {
    let (p, f) = schonhardt();
    assert_eq!(check_faces("schonhardt", &(p.clone(), f.clone())), 0);
    let r = recover_faces(&p, &f, 0, 5).unwrap();
    assert!(r.steiner + r.steiner_on_faces > 0);
}

/// A large tetrahedron with a twisted triangular prism hovering just above
/// the middle of its base (a void) and one just below (mirrored, outside),
/// placed by `rng`: the edges between them cross the base where flips
/// cannot clear them.
fn twisted_sandwich(rng: &mut Rng) -> Surface {
    let mut p: Vec<Point3> = vec![
        [-4.0, -4.0, 0.0],
        [4.0, -4.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 4.0],
    ];
    let mut f: Vec<[u32; 3]> = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
    for below in [false, true] {
        let gap = 0.005 + 0.05 * rng.unit();
        let h = 0.1 + 0.4 * rng.unit();
        let turn = 0.3 + 0.4 * rng.unit();
        let r = 0.2 + 0.2 * rng.unit();
        let a0 = rng.unit() * std::f64::consts::TAU;
        let (cx, cy) = (0.2 * rng.unit() - 0.1, 0.2 * rng.unit() - 0.1);
        let z = |d: f64| if below { -(gap + d) } else { gap + d };
        let base = p.len() as u32;
        for (dz, tw) in [(0.0, 0.0), (h, turn)] {
            for k in 0..3 {
                let a = a0 + tw + f64::from(k) * std::f64::consts::TAU / 3.0;
                p.push([cx + r * libm::cos(a), cy + r * libm::sin(a), z(dz)]);
            }
        }
        let mut t = vec![[0u32, 2, 1], [3, 4, 5]];
        for k in 0..3u32 {
            let (a, b) = (k, (k + 1) % 3);
            t.push([a, b, b + 3]);
            t.push([a, b + 3, a + 3]);
        }
        // Reversed: above, a void; below, mirrored, so outward again.
        f.extend(t.iter().map(|q| [q[0] + base, q[2] + base, q[1] + base]));
    }
    (p, f)
}

/// A star-shaped sphere: a UV sphere whose vertices sit at random radii
/// between `lo` and 1.
fn star_sphere(segments: u32, rings: u32, lo: f64, rng: &mut Rng) -> Surface {
    let mut radius = || lo + (1.0 - lo) * rng.unit();
    let mut p = vec![[0.0, 0.0, radius()]];
    for i in 1..rings {
        let t = f64::from(i) / f64::from(rings) * std::f64::consts::PI;
        for j in 0..segments {
            let a = f64::from(j) / f64::from(segments) * std::f64::consts::TAU;
            let r = radius();
            p.push([
                r * libm::sin(t) * libm::cos(a),
                r * libm::sin(t) * libm::sin(a),
                r * libm::cos(t),
            ]);
        }
    }
    let south = p.len() as u32;
    p.push([0.0, 0.0, -radius()]);
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

/// A large tetrahedron with spiky star-shaped voids hovering just above its
/// base (a void) and just below it (outside).
fn spiky_sandwich(rng: &mut Rng) -> Surface {
    let mut p: Vec<Point3> = vec![
        [-4.0, -4.0, 0.0],
        [4.0, -4.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 4.0],
    ];
    let mut f: Vec<[u32; 3]> = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
    for below in [false, true] {
        let (q, t) = star_sphere(10, 7, 0.05, rng);
        let gap = 0.002 + 0.02 * rng.unit();
        let s = 0.3;
        let (cx, cy) = (0.2 * rng.unit() - 0.1, 0.2 * rng.unit() - 0.1);
        let base = p.len() as u32;
        let zmin = q.iter().map(|v| v[2]).fold(f64::INFINITY, f64::min);
        let zmax = q.iter().map(|v| v[2]).fold(f64::NEG_INFINITY, f64::max);
        for v in &q {
            let z = if below {
                -gap - s * (zmax - v[2])
            } else {
                gap + s * (v[2] - zmin)
            };
            p.push([cx + s * v[0], cy + s * v[1], z]);
        }
        // Above, reversed: a void.
        f.extend(t.iter().map(|q| {
            let q = q.map(|v| v + base);
            if below { q } else { [q[0], q[2], q[1]] }
        }));
    }
    (p, f)
}

fn seeded(seed: u64) -> Rng {
    Rng(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1)
}

/// Facets flips cannot recover go to the cavity stage: refilled from
/// Delaunay tetrahedralizations of its halves, or, where expanding the
/// cavity fails, by gift-wrapping (with symbolic perturbation).
#[test]
fn facet_cavities_are_refilled() {
    for seed in 1..=5 {
        let s = twisted_sandwich(&mut seeded(seed));
        let r = recover_faces(&s.0, &s.1, 0, 5).unwrap();
        assert!(r.crossed[1] > 0, "seed {seed}: flips should leave the base");
        assert_eq!(r.missing, Vec::<u32>::new(), "seed {seed}");
        assert_eq!(r.steiner_on_faces, 0, "seed {seed}");
        check_region(&format!("twisted sandwich {seed}"), &s);
    }
    let s = spiky_sandwich(&mut seeded(127));
    let r = recover_faces(&s.0, &s.1, 0, 5).unwrap();
    assert!(r.wrapped > 0, "a cavity should need gift-wrapping");
    assert_eq!(r.missing, Vec::<u32>::new());
    assert_eq!(r.steiner_on_faces, 0);
    check_region("spiky sandwich", &s);
}

/// The disjoint-cavity case (Diazzi et al., Section 3.5): growing a facet
/// cavity reaches a tet across the facet's plane, so the retetrahedralized
/// parts above and below it would overlap, as TetGen's face recovery
/// assumes they cannot. It is refused, and the facet is recovered another
/// way, with the region still exactly the inside.
#[test]
fn a_cavity_reaching_across_its_plane_is_refused() {
    let s = star_sphere(24, 14, 0.01, &mut seeded(161));
    let r = recover_faces(&s.0, &s.1, 0, 5).unwrap();
    assert!(r.across_plane > 0, "an expansion should reach across");
    assert_eq!(r.missing, Vec::<u32>::new());
    check_region("star sphere", &s);
}

/// A unit cube with a cube of side 0.4 removed from its middle (the inner
/// shell faces inward).
fn hollow_cube() -> Surface {
    let (mut p, mut f) = cube();
    let (inner, faces) = cube();
    let n = p.len() as u32;
    p.extend(inner.iter().map(|q| q.map(|x| 0.3 + 0.4 * x)));
    f.extend(faces.iter().map(|t| [t[0] + n, t[2] + n, t[1] + n]));
    (p, f)
}

/// The inside tets fill exactly the volume the surface encloses, and the
/// boundary is one outward face per subface, covering the surface.
fn check_region(name: &str, (p, f): &Surface) {
    let r = mesh_region(p, f, 0, 5).unwrap_or_else(|e| panic!("{name}: {e}"));
    let at = |v: u32| r.dt().points()[v as usize].approx();
    let det = |a: Point3, b: Point3, c: Point3| {
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
    };
    let enclosed: f64 = f
        .iter()
        .map(|t| det(p[t[0] as usize], p[t[1] as usize], p[t[2] as usize]) / 6.0)
        .sum::<f64>()
        .abs();
    let volume: f64 = r
        .tets
        .iter()
        .map(|t| mesher_tet::quality::float_volume(t.map(at)))
        .sum();
    // The divergence theorem over the boundary, outward, gives it again.
    let through_boundary: f64 = r
        .boundary
        .iter()
        .map(|(t, _)| det(at(t[0]), at(t[1]), at(t[2])) / 6.0)
        .sum();
    let tol = 1e-9 * enclosed.max(1.0);
    assert!(
        (volume - enclosed).abs() <= tol,
        "{name}: {volume} against {enclosed}"
    );
    assert!(
        (through_boundary - enclosed).abs() <= tol,
        "{name}: boundary {through_boundary}"
    );
    assert_eq!(r.boundary.len(), r.recovered.subfaces.len(), "{name}");
    for &(face, triangle) in &r.boundary {
        let s = r.recovered.subfaces.iter().find(|s| {
            let (mut a, mut b) = (s.v, face);
            a.sort_unstable();
            b.sort_unstable();
            a == b
        });
        assert_eq!(s.map(|s| s.triangle), Some(triangle), "{name}");
    }
    eprintln!(
        "{name:<20} {:>6} tets inside, {:>5} boundary faces",
        r.tets.len(),
        r.boundary.len()
    );
}

#[test]
fn regions_are_the_inside_of_the_surface() {
    check_region("cube", &cube());
    check_region("hollow cube", &hollow_cube());
    check_region("uv sphere", &uv_sphere(16, 10, [0.0; 3], false));
    check_region(
        "rounded uv sphere",
        &uv_sphere(16, 10, [1.0, 1.0, 1.5], true),
    );
    check_region("sharp wedge", &sharp_wedge());
    check_region("schonhardt", &schonhardt());
    check_region("sandwich", &sandwich());
    check_region("fine uv sphere", &uv_sphere(48, 24, [0.3, -0.2, 0.1], true));
}

/// Rounded to f64, every tet is strictly positive (exactly, on the rounded
/// coordinates), the volume is the enclosed volume to rounding, and the
/// boundary is one outward face per subface.
fn check_rounded(name: &str, (p, f): &Surface) {
    let region = mesh_region(p, f, 0, 5).unwrap_or_else(|e| panic!("{name}: {e}"));
    let subfaces = region.recovered.subfaces.len();
    let m = round_region(region, p, f).unwrap_or_else(|e| panic!("{name}: {e}"));
    let explicit: Vec<Point> = m
        .points
        .iter()
        .map(|&q| Point::explicit(q).unwrap())
        .collect();
    for t in &m.tets {
        let [a, b, c, d] = t.map(|v| &explicit[v as usize]);
        assert_eq!(
            mesher_geom::implicit::orient3d(a, b, c, d),
            mesher_geom::Sign::Positive,
            "{name}: {t:?}"
        );
    }
    let det = |a: Point3, b: Point3, c: Point3| {
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
    };
    let enclosed: f64 = f
        .iter()
        .map(|t| det(p[t[0] as usize], p[t[1] as usize], p[t[2] as usize]) / 6.0)
        .sum::<f64>()
        .abs();
    let volume: f64 = m
        .tets
        .iter()
        .map(|t| mesher_tet::quality::float_volume(t.map(|v| m.points[v as usize])))
        .sum();
    assert!(
        (volume - enclosed).abs() <= 1e-9 * enclosed.max(1.0),
        "{name}: {volume} against {enclosed}"
    );
    assert_eq!(m.boundary.len(), subfaces, "{name}");
    eprintln!(
        "{name:<20} {:>6} tets, {:>4} points rounded, {} turned over, {} fixes",
        m.tets.len(),
        m.rounded,
        m.inverted_before,
        m.untangles.len()
    );
}

#[test]
fn rounded_meshes_are_valid() {
    check_rounded("cube", &cube());
    check_rounded("hollow cube", &hollow_cube());
    check_rounded(
        "rounded uv sphere",
        &uv_sphere(16, 10, [1.0, 1.0, 1.5], true),
    );
    check_rounded("sharp wedge", &sharp_wedge());
    check_rounded("schonhardt", &schonhardt());
    check_rounded("sandwich", &sandwich());
    check_rounded("fine uv sphere", &uv_sphere(48, 24, [0.3, -0.2, 0.1], true));
    check_rounded("spiky sandwich", &spiky_sandwich(&mut seeded(127)));
    check_rounded("star sphere", &star_sphere(24, 14, 0.01, &mut seeded(161)));
    for seed in 1..=5 {
        check_rounded("twisted sandwich", &twisted_sandwich(&mut seeded(seed)));
    }
}
