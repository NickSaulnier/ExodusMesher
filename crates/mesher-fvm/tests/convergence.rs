//! The downstream check: a manufactured Poisson solution on meshes of the unit
//! cube from both kernels; the error must fall as the mesh is refined.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::f64::consts::PI;

use mesher_fvm::{Scheme, solve_manufactured, solve_with};
use mesher_geom::{Aabb, Point3};
use mesher_tet::{ConformingParams, StuffingParams, conforming_mesh, fv_quality, stuff};

fn exact(p: Point3) -> f64 {
    (PI * p[0]).sin() * (PI * p[1]).sin() * (PI * p[2]).sin()
}

fn source(p: Point3) -> f64 {
    3.0 * PI * PI * exact(p)
}

fn cube() -> (Vec<Point3>, Vec<[u32; 3]>) {
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
    let t = quads
        .iter()
        .flat_map(|&[a, b, c, d]| [[a, b, c], [a, c, d]])
        .collect();
    (p, t)
}

#[test]
fn conforming_meshes_converge() {
    let (p, t) = cube();
    let mut errors = Vec::new();
    for size in [0.25, 0.125, 0.0625] {
        let m = conforming_mesh(&p, &t, &ConformingParams::with_size(size)).unwrap();
        let r = solve_manufactured(&m.points, &m.tets, &exact, &source).unwrap();
        let two = solve_with(&m.points, &m.tets, &exact, &source, Scheme::TwoPoint).unwrap();
        let q = fv_quality(&m.points, &m.tets);
        eprintln!(
            "conforming h={size}: {} cells, L2 {:.3e} corrected ({} passes), {:.3e} two-point; non-orthogonality max {:.1}° mean {:.1}°",
            r.cells,
            r.l2_error,
            r.correction_passes,
            two.l2_error,
            q.max_non_orthogonality_degrees,
            q.mean_non_orthogonality_degrees
        );
        errors.push(r.l2_error);
    }
    // With the non-orthogonal correction the error falls by a clear factor
    // per halving of the element size.
    assert!(
        errors.windows(2).all(|w| w[1] < 0.6 * w[0]),
        "error did not fall: {errors:?}"
    );
}

#[test]
fn stuffed_meshes_converge() {
    let level = |q: Point3| {
        (0..3)
            .map(|k| (q[k] - 0.5).abs() - 0.5)
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let bounds = Aabb::of_points(&[[0.0; 3], [1.0; 3]]);
    let mut errors = Vec::new();
    for h in [0.1, 0.05] {
        let m = stuff(level, &bounds, &StuffingParams::with_cell_size(h)).unwrap();
        let r = solve_manufactured(&m.points, &m.tets, &exact, &source).unwrap();
        eprintln!("stuffing h={h}: {} cells, L2 {:.3e}", r.cells, r.l2_error);
        errors.push(r.l2_error);
    }
    assert!(
        errors[1] < 0.7 * errors[0],
        "error did not fall: {errors:?}"
    );
}
