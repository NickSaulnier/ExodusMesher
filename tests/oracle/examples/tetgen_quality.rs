//! TetGen's finite-volume quality on triangle OBJ files, for comparing with
//! `mesher-tet`'s `--example improve` (Milestone 15):
//!
//! ```text
//! cargo run --release -p mesher-oracle --example tetgen_quality -- a.obj b.obj ...
//! ```
//!
//! TetGen meshes each file with `pzAfq` and `-a` the regular-tet volume of
//! the improve example's element size (0.03 of the bounding-box diagonal).
//! Prints a line per file in the improve example's format (the same
//! quality twice, as before and after), so its summaries read both: tets,
//! the smallest dihedral angle, tets under 1° and under 10°, tets flat in
//! floating point, the largest non-orthogonality and skewness, and both
//! over wall-adjacent cells; or the error. Test-only: TetGen is AGPL and is
//! never shipped.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::time::Instant;

use mesher_geom::predicates::orient3d;
use mesher_geom::{Point3, Sign};
use mesher_tet::{dihedral_angles, fv_quality};
use tritet::Tetgen;

fn describe(points: &[Point3], tets: &[[u32; 4]]) -> String {
    let q = fv_quality(points, tets);
    let (mut under_one, mut under_ten) = (0, 0);
    for t in tets {
        let min = dihedral_angles(t.map(|v| points[v as usize]))
            .into_iter()
            .fold(180.0, f64::min);
        under_one += usize::from(min < 1.0);
        under_ten += usize::from(min < 10.0);
    }
    format!(
        "{} tets, min {:.3}, {under_one} <1, {under_ten} <10, {} flat, nonorth {:.2}, skew {:.3}, wall nonorth {:.2}, wall skew {:.3}",
        tets.len(),
        q.min_dihedral_degrees,
        q.flat_in_float,
        q.max_non_orthogonality_degrees,
        q.max_skewness,
        q.wall_max_non_orthogonality_degrees,
        q.wall_max_skewness
    )
}

fn tetgen(
    points: &[Point3],
    triangles: &[[u32; 3]],
    max_volume: f64,
) -> Result<(Vec<Point3>, Vec<[u32; 4]>), String> {
    let mut tg = Tetgen::new(points.len(), Some(vec![3; triangles.len()]), None, None)
        .map_err(str::to_owned)?;
    for (i, q) in points.iter().enumerate() {
        tg.set_point(i, 0, q[0], q[1], q[2])
            .map_err(str::to_owned)?;
    }
    for (i, t) in triangles.iter().enumerate() {
        for (m, &v) in t.iter().enumerate() {
            tg.set_facet_point(i, m, v as usize)
                .map_err(str::to_owned)?;
        }
    }
    tg.generate_mesh(false, false, Some(max_volume), None)
        .map_err(str::to_owned)?;
    let out: Vec<Point3> = (0..tg.out_npoint())
        .map(|i| [tg.out_point(i, 0), tg.out_point(i, 1), tg.out_point(i, 2)])
        .collect();
    let tets = (0..tg.out_ncell())
        .map(|c| {
            let mut t = [0, 1, 2, 3].map(|m| tg.out_cell_point(c, m) as u32);
            let q = t.map(|v| out[v as usize]);
            if orient3d(q[0], q[1], q[2], q[3]) == Sign::Negative {
                t.swap(0, 1);
            }
            t
        })
        .collect();
    Ok((out, tets))
}

fn main() {
    for path in std::env::args().skip(1) {
        let text = std::fs::read_to_string(&path).unwrap();
        let mut points = Vec::new();
        let mut triangles = Vec::new();
        for line in text.lines() {
            let mut w = line.split_whitespace();
            match w.next() {
                Some("v") => {
                    let c: Vec<f64> = w.take(3).map(|x| x.parse().unwrap()).collect();
                    points.push([c[0], c[1], c[2]]);
                }
                Some("f") => {
                    let v: Vec<u32> = w
                        .map(|x| x.split('/').next().unwrap().parse::<u32>().unwrap() - 1)
                        .collect();
                    triangles.push([v[0], v[1], v[2]]);
                }
                _ => {}
            }
        }
        let (lo, hi) = points.iter().fold(
            ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
            |(lo, hi), p: &[f64; 3]| {
                (
                    [0, 1, 2].map(|k| lo[k].min(p[k])),
                    [0, 1, 2].map(|k| hi[k].max(p[k])),
                )
            },
        );
        let diagonal = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
        let h = 0.03 * diagonal;
        let max_volume = 8.0 / (9.0 * libm::sqrt(3.0)) * h * h * h;
        let started = Instant::now();
        let result = tetgen(&points, &triangles, max_volume);
        let seconds = started.elapsed().as_secs_f64();
        match result {
            Ok((p, t)) if !t.is_empty() => {
                let q = describe(&p, &t);
                println!("{path}\tok\t{q}\t{q}\tTetGen\t{seconds:.2}\t0");
            }
            Ok(_) => println!("{path}\tfailed\tno tetrahedra\t{seconds:.2}"),
            Err(e) => println!("{path}\tfailed\t{e}\t{seconds:.2}"),
        }
    }
}
