//! Isosurface stuffing on analytic shapes at random positions, sizes and cell
//! sizes: the dihedral-angle bounds hold, no tet is degenerate, and the mesh
//! conforms (every face is shared by two tets or lies on a closed, manifold
//! boundary).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeMap;

use mesher_geom::predicates::orient3d;
use mesher_geom::{Aabb, Point3, Sign};
use mesher_tet::stuffing::{DIHEDRAL_MAX_DEGREES, DIHEDRAL_MIN_DEGREES};
use mesher_tet::{Stuffed, StuffingParams, dihedral_range, stuff, stuff_graded};
use proptest::prelude::*;

#[derive(Debug, Clone)]
enum Shape {
    Sphere { c: Point3, r: f64 },
    Ellipsoid { c: Point3, axes: Point3 },
    Box { c: Point3, half: Point3 },
    Torus { c: Point3, big: f64, small: f64 },
}

impl Shape {
    fn level(&self, p: Point3) -> f64 {
        match *self {
            Shape::Sphere { c, r } => {
                ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt() - r
            }
            Shape::Ellipsoid { c, axes } => {
                (0..3)
                    .map(|k| ((p[k] - c[k]) / axes[k]).powi(2))
                    .sum::<f64>()
                    - 1.0
            }
            Shape::Box { c, half } => (0..3)
                .map(|k| (p[k] - c[k]).abs() - half[k])
                .fold(f64::NEG_INFINITY, f64::max),
            Shape::Torus { c, big, small } => {
                let q = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
                let ring = (q[0] * q[0] + q[1] * q[1]).sqrt() - big;
                (ring * ring + q[2] * q[2]).sqrt() - small
            }
        }
    }

    fn bounds(&self) -> Aabb {
        let (c, e) = match *self {
            Shape::Sphere { c, r } => (c, [r; 3]),
            Shape::Ellipsoid { c, axes } => (c, axes),
            Shape::Box { c, half } => (c, half),
            Shape::Torus { c, big, small } => (c, [big + small, big + small, small]),
        };
        Aabb::of_points(&[
            [c[0] - e[0], c[1] - e[1], c[2] - e[2]],
            [c[0] + e[0], c[1] + e[1], c[2] + e[2]],
        ])
    }
}

fn shape() -> impl Strategy<Value = Shape> {
    let c = prop::array::uniform3(-3.0f64..3.0);
    prop_oneof![
        (c.clone(), 0.5f64..2.0).prop_map(|(c, r)| Shape::Sphere { c, r }),
        (c.clone(), prop::array::uniform3(0.3f64..2.0))
            .prop_map(|(c, axes)| Shape::Ellipsoid { c, axes }),
        (c.clone(), prop::array::uniform3(0.3f64..1.5))
            .prop_map(|(c, half)| Shape::Box { c, half }),
        (c, 1.0f64..2.0, 0.3f64..0.8).prop_map(|(c, big, small)| Shape::Torus { c, big, small }),
    ]
}

/// Points on the shape's surface: sign changes along grid lines through its
/// bounds, bisected.
fn surface_seeds(shape: &Shape) -> Vec<Point3> {
    let b = shape.bounds();
    let n = 24;
    let mut seeds = Vec::new();
    for a in 0..3 {
        let (u, v) = ((a + 1) % 3, (a + 2) % 3);
        for i in 0..=n {
            for j in 0..=n {
                let mut p = [0.0; 3];
                p[u] = b.min[u] + (b.max[u] - b.min[u]) * f64::from(i) / f64::from(n);
                p[v] = b.min[v] + (b.max[v] - b.min[v]) * f64::from(j) / f64::from(n);
                let at = |t: f64| {
                    let mut q = p;
                    q[a] = b.min[a] - 0.1 + (b.max[a] - b.min[a] + 0.2) * t;
                    q
                };
                let steps = 4 * n;
                for k in 0..steps {
                    let (t0, t1) = (
                        f64::from(k) / f64::from(steps),
                        f64::from(k + 1) / f64::from(steps),
                    );
                    let (l0, l1) = (shape.level(at(t0)), shape.level(at(t1)));
                    if (l0 < 0.0) != (l1 < 0.0) {
                        let (mut lo, mut hi) = (t0, t1);
                        for _ in 0..40 {
                            let mid = 0.5 * (lo + hi);
                            if (shape.level(at(mid)) < 0.0) == (l0 < 0.0) {
                                lo = mid;
                            } else {
                                hi = mid;
                            }
                        }
                        seeds.push(at(lo));
                    }
                }
            }
        }
    }
    seeds
}

fn check(shape: &Shape, h: f64) -> Result<(), TestCaseError> {
    let out = stuff(
        |p| shape.level(p),
        &shape.bounds(),
        &StuffingParams::with_cell_size(h),
    )
    .unwrap();
    check_mesh(shape, h, &out)
}

fn check_mesh(shape: &Shape, h: f64, out: &Stuffed) -> Result<(), TestCaseError> {
    prop_assert_eq!(out.dropped, 0);
    prop_assert!(!out.tets.is_empty());
    for t in &out.tets {
        let [a, b, c, d] = t.map(|v| out.points[v as usize]);
        prop_assert_eq!(orient3d(a, b, c, d), Sign::Positive);
    }
    let (lo, hi) = dihedral_range(&out.points, &out.tets).unwrap();
    prop_assert!(
        lo >= DIHEDRAL_MIN_DEGREES,
        "min dihedral {lo} for {shape:?} at h = {h}"
    );
    prop_assert!(
        hi <= DIHEDRAL_MAX_DEGREES,
        "max dihedral {hi} for {shape:?} at h = {h}"
    );

    // Conformity: faces shared by at most two tets, with opposite orientation.
    let mut faces: BTreeMap<[u32; 3], Vec<[u32; 3]>> = BTreeMap::new();
    for t in &out.tets {
        for f in [
            [t[0], t[1], t[3]],
            [t[1], t[2], t[3]],
            [t[0], t[3], t[2]],
            [t[0], t[2], t[1]],
        ] {
            let mut k = f;
            k.sort_unstable();
            faces.entry(k).or_default().push(f);
        }
    }
    let mut boundary_edges: BTreeMap<(u32, u32), i32> = BTreeMap::new();
    for (key, uses) in &faces {
        prop_assert!(
            uses.len() <= 2,
            "face {key:?} shared by {} tets",
            uses.len()
        );
        if uses.len() == 1 {
            let f = uses[0];
            for (u, v) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
                *boundary_edges.entry((u.min(v), u.max(v))).or_default() +=
                    if u < v { 1 } else { -1 };
            }
        }
    }
    // A closed, consistently oriented boundary: every boundary edge is
    // traversed once in each direction.
    prop_assert!(
        boundary_edges.values().all(|&n| n == 0),
        "boundary is not closed"
    );
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn bounds_hold_and_meshes_conform(s in shape(), cells in 6.0f64..18.0) {
        let b = s.bounds();
        let extent = (0..3).map(|k| b.max[k] - b.min[k]).fold(0.0, f64::max);
        check(&s, extent / cells)?;
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    /// The graded background (fine at the surface, up to two levels coarser
    /// inside) keeps the same bounds and conformity.
    #[test]
    fn graded_bounds_hold_and_meshes_conform(
        s in shape(),
        cells in 8.0f64..20.0,
        levels in 1u32..=2,
    ) {
        let b = s.bounds();
        let extent = (0..3).map(|k| b.max[k] - b.min[k]).fold(0.0, f64::max);
        let h = extent / cells;
        let out = stuff_graded(
            |p| s.level(p),
            &b,
            &surface_seeds(&s),
            &StuffingParams::with_cell_size(h),
            levels,
        )
        .unwrap();
        check_mesh(&s, h, &out)?;
    }
}

/// Grading keeps the surface as fine as the uniform lattice of the same
/// leaf size (the same volume, to the lattice's accuracy) with fewer tets.
#[test]
fn graded_sphere_is_coarser_inside() {
    let s = Shape::Sphere {
        c: [0.1, -0.2, 0.3],
        r: 1.0,
    };
    let h = 0.05;
    let uniform = stuff(
        |p| s.level(p),
        &s.bounds(),
        &StuffingParams::with_cell_size(h),
    )
    .unwrap();
    let graded = stuff_graded(
        |p| s.level(p),
        &s.bounds(),
        &surface_seeds(&s),
        &StuffingParams::with_cell_size(h),
        2,
    )
    .unwrap();
    check_mesh(&s, h, &graded).unwrap();
    let volume = |m: &Stuffed| -> f64 {
        m.tets
            .iter()
            .map(|t| {
                let [a, b, c, d] = t.map(|v| m.points[v as usize]);
                let (u, v, w) = (
                    [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                    [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
                    [d[0] - a[0], d[1] - a[1], d[2] - a[2]],
                );
                (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                    + u[2] * (v[0] * w[1] - v[1] * w[0]))
                    / 6.0
            })
            .sum()
    };
    let (vu, vg) = (volume(&uniform), volume(&graded));
    assert!((vu - vg).abs() / vu < 1e-9, "{vu} vs {vg}");
    assert!(
        graded.tets.len() * 2 < uniform.tets.len(),
        "graded {} vs uniform {}",
        graded.tets.len(),
        uniform.tets.len()
    );
}
