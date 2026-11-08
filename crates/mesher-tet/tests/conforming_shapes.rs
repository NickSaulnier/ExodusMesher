//! The conforming mesher on shapes that exercise curvature, genus, concavity,
//! cavities and coplanar degeneracy: it must succeed, keep every tet positive,
//! reproduce the enclosed volume to rounding, and conform to the surface.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeMap;
use std::time::Instant;

use mesher_geom::predicates::orient3d;
use mesher_geom::vec3::{cross, dot, sub};
use mesher_geom::{Point3, Sign};
use mesher_tet::{ConformingParams, conforming_mesh, dihedral_range};

type Surface = (Vec<Point3>, Vec<[u32; 3]>);

fn icosphere(n: u32) -> Surface {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut p: Vec<Point3> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut f: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..n {
        let mut mid: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        let mut next = Vec::new();
        for [a, b, c] in f {
            let mut m = |i: u32, j: u32| {
                *mid.entry((i.min(j), i.max(j))).or_insert_with(|| {
                    let (x, y) = (p[i as usize], p[j as usize]);
                    p.push([0, 1, 2].map(|k| (x[k] + y[k]) / 2.0));
                    (p.len() - 1) as u32
                })
            };
            let (ab, bc, ca) = (m(a, b), m(b, c), m(c, a));
            next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        f = next;
    }
    for q in &mut p {
        let r = dot(*q, *q).sqrt();
        *q = q.map(|x| x / r);
    }
    (p, f)
}

fn torus(nu: u32, nv: u32) -> Surface {
    let mut p = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let u = f64::from(i) / f64::from(nu) * std::f64::consts::TAU;
            let v = f64::from(j) / f64::from(nv) * std::f64::consts::TAU;
            let a = 2.0 + 0.75 * v.cos();
            p.push([a * u.cos(), a * u.sin(), 0.75 * v.sin()]);
        }
    }
    let id = |i: u32, j: u32| (i % nu) * nv + (j % nv);
    let mut f = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            f.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
            f.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
        }
    }
    (p, f)
}

/// Axis-aligned box as 12 triangles; `flip` reverses them (for a cavity).
fn boxed(lo: Point3, hi: Point3, offset: u32, flip: bool) -> Surface {
    let p: Vec<Point3> = (0..8)
        .map(|i| {
            [
                if i & 1 == 0 { lo[0] } else { hi[0] },
                if i & 2 == 0 { lo[1] } else { hi[1] },
                if i & 4 == 0 { lo[2] } else { hi[2] },
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
        .map(|t: [u32; 3]| {
            let t = t.map(|v| v + offset);
            if flip { [t[0], t[2], t[1]] } else { t }
        })
        .collect();
    (p, f)
}

fn l_prism() -> Surface {
    let foot = [
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ];
    let mut p: Vec<Point3> = foot.iter().map(|q| [q[0], q[1], 0.0]).collect();
    p.extend(foot.iter().map(|q| [q[0], q[1], 1.0]));
    // Bottom (outward -z) and top (outward +z) as fans that respect the notch.
    let bottom = [[0, 3, 1], [1, 3, 2], [0, 5, 4], [0, 4, 3]];
    let mut f: Vec<[u32; 3]> = bottom.to_vec();
    f.extend(bottom.iter().map(|t| [t[0] + 6, t[2] + 6, t[1] + 6]));
    for s in 0..6u32 {
        let (a, b) = (s, (s + 1) % 6);
        f.push([a, b, b + 6]);
        f.push([a, b + 6, a + 6]);
    }
    (p, f)
}

fn grid_cube(n: i64) -> Surface {
    let mut ids: BTreeMap<[i64; 3], u32> = BTreeMap::new();
    let (mut p, mut f) = (Vec::new(), Vec::new());
    let (x, y, z) = ([1, 0, 0], [0, 1, 0], [0, 0, 1]);
    for (o, u, v) in [
        ([0, 0, 0], y, x),
        ([0, 0, n], x, y),
        ([0, 0, 0], x, z),
        ([0, n, 0], z, x),
        ([0, 0, 0], z, y),
        ([n, 0, 0], y, z),
    ] {
        let mut id = |i: i64, j: i64| {
            let q = [0, 1, 2].map(|k| o[k] + i * u[k] + j * v[k]);
            *ids.entry(q).or_insert_with(|| {
                p.push(q.map(|c| c as f64));
                (p.len() - 1) as u32
            })
        };
        for i in 0..n {
            for j in 0..n {
                let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
                f.push([a, b, c]);
                f.push([a, c, d]);
            }
        }
    }
    (p, f)
}

fn enclosed(p: &[Point3], f: &[[u32; 3]]) -> f64 {
    f.iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| p[v as usize]);
            dot(a, cross(b, c)) / 6.0
        })
        .sum()
}

fn check(name: &str, (p, f): &Surface, size: f64) {
    let start = Instant::now();
    let m = conforming_mesh(p, f, &ConformingParams::with_size(size))
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    let elapsed = start.elapsed();
    let volume: f64 = m
        .tets
        .iter()
        .map(|t| {
            let [a, b, c, d] = t.map(|v| m.points[v as usize]);
            assert_eq!(
                orient3d(a, b, c, d),
                Sign::Positive,
                "{name}: non-positive tet"
            );
            dot(sub(b, a), cross(sub(c, a), sub(d, a))) / 6.0
        })
        .sum();
    let expected = enclosed(p, f);
    assert!(
        (volume - expected).abs() <= 1e-9 * expected.abs(),
        "{name}: volume {volume} vs {expected}"
    );
    let (lo, hi) = dihedral_range(&m.points, &m.tets).unwrap();
    eprintln!(
        "{name:<12} {:>7} tets {:>6} pts (+{} boundary, +{} interior Steiner) dihedral {lo:5.1}..{hi:5.1} slivers<10°: {} [{elapsed:?}]",
        m.tets.len(),
        m.points.len(),
        m.boundary_steiner,
        m.interior_steiner,
        m.slivers_left
    );
    // The boundary covers the surface: the same total area, and every
    // boundary triangle lies in the plane of the input triangle it maps to.
    let area = |a: Point3, b: Point3, c: Point3| {
        dot(cross(sub(b, a), sub(c, a)), cross(sub(b, a), sub(c, a))).sqrt() / 2.0
    };
    let total: f64 = f
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| p[v as usize]);
            area(a, b, c)
        })
        .sum();
    let covered: f64 = m
        .boundary
        .iter()
        .map(|&(t, _)| {
            let [a, b, c] = t.map(|v| m.points[v as usize]);
            area(a, b, c)
        })
        .sum();
    assert!(
        (covered - total).abs() <= 1e-9 * total,
        "{name}: boundary area {covered} vs {total}"
    );
    for &(t, o) in &m.boundary {
        let [a, b, c] = f[o as usize].map(|v| p[v as usize]);
        let n = cross(sub(b, a), sub(c, a));
        let scale = dot(n, n).sqrt();
        for q in t.map(|v| m.points[v as usize]) {
            let d = dot(sub(q, a), n) / scale;
            assert!(
                d.abs() <= 1e-9,
                "{name}: boundary point {d} off its input plane"
            );
        }
    }
    if m.flattened > 0 {
        eprintln!("{name:<12} flattened {} rounding-flat slivers", m.flattened);
    }
    {
        assert_eq!(m.slivers_left, 0, "{name}: slivers left");
    }
}

#[test]
fn shapes_mesh_conformally() {
    check("cube", &boxed([0.0; 3], [1.0; 3], 0, false), 0.2);
    // Planar faces off the axes: Steiner points on them are off their planes
    // by rounding.
    let (p, f) = grid_cube(3);
    let (s1, c1, s2, c2) = (0.3f64.sin(), 0.3f64.cos(), 0.7f64.sin(), 0.7f64.cos());
    let p = p
        .into_iter()
        .map(|[x, y, z]| {
            let (x, y) = (c1 * x - s1 * y, s1 * x + c1 * y);
            let (y, z) = (c2 * y - s2 * z, s2 * y + c2 * z);
            [x, y, z]
        })
        .collect();
    check("tilted cube", &(p, f), 0.6);
    check("grid cube", &grid_cube(3), 0.6);
    check("icosphere", &icosphere(2), 0.25);
    check("torus", &torus(24, 12), 0.3);
    check("L-prism", &l_prism(), 0.25);
    let (mut p, mut f) = boxed([0.0; 3], [3.0; 3], 0, false);
    let (q, g) = boxed([1.0; 3], [2.0; 3], 8, true);
    p.extend(q);
    f.extend(g);
    check("hollow cube", &(p, f), 0.5);
}

/// A cube whose first triangle is split around a point 0.2% off one edge: the
/// tiny input angles must not make refinement run away, at any placement.
#[test]
fn sliver_input_terminates() {
    for (scale, offset) in [(0.001, 1000.0), (1.0, 0.0), (10000.0, -50000.0)] {
        let (mut p, mut f) = boxed([0.0; 3], [1.0; 3], 0, false);
        let [a, b, c] = f[0];
        let (pa, pb, pc) = (p[a as usize], p[b as usize], p[c as usize]);
        let mid = [0, 1, 2].map(|k| (pa[k] + pb[k]) / 2.0);
        p.push([0, 1, 2].map(|k| mid[k] + 0.002 * (pc[k] - mid[k])));
        let m = p.len() as u32 - 1;
        f[0] = [a, b, m];
        f.push([b, c, m]);
        f.push([c, a, m]);
        for q in &mut p {
            *q = q.map(|x| x * scale + offset);
        }
        let start = Instant::now();
        let r = conforming_mesh(&p, &f, &ConformingParams::with_size(0.104 * scale));
        let elapsed = start.elapsed();
        match &r {
            Ok(m) => eprintln!(
                "sliver x{scale} +{offset}: {} tets, {} pts [{elapsed:?}]",
                m.tets.len(),
                m.points.len()
            ),
            Err(e) => eprintln!("sliver x{scale} +{offset}: {e} [{elapsed:?}]"),
        }
        assert!(elapsed.as_secs() < 30, "took {elapsed:?}");
    }
}

/// The same cube at several scales: power-of-two scales mesh identically,
/// and rounding at the others leaves no slivers either.
#[test]
fn cube_quality_does_not_depend_on_scale() {
    let mut reference = None;
    for scale in [1.0, 1e-6, 1e7, 1.0 / 1_048_576.0, 1e-3, 1e-5] {
        let (mut p, f) = boxed([0.0; 3], [1.0; 3], 0, false);
        for q in &mut p {
            *q = q.map(|x| x * scale);
        }
        let h = 0.06 * 3f64.sqrt() * scale;
        let m = conforming_mesh(&p, &f, &ConformingParams::with_size(h)).unwrap();
        let (lo, hi) = dihedral_range(&m.points, &m.tets).unwrap();
        eprintln!(
            "scale {scale:e}: {} tets, dihedral {lo:.2}..{hi:.2}",
            m.tets.len()
        );
        assert_eq!(m.slivers_left, 0, "scale {scale:e}: {lo:.2}");
        if scale == 1.0 {
            reference = Some(m.tets.len());
        }
        if scale == 1.0 / 1_048_576.0 {
            assert_eq!(Some(m.tets.len()), reference);
        }
    }
}
