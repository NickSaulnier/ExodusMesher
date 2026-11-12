//! Worst-cell improvement on triangle OBJ files, for measuring (Milestone 15):
//!
//! ```text
//! cargo run --release -p mesher-tet --example improve -- a.obj b.obj ...
//! ```
//!
//! Meshes each file as the pipeline does (the conforming kernel at 0.03 of
//! the diagonal, graded by the size field; `IMPROVE_KERNEL=exact` for the
//! exact kernel rounded to f64), then improves it (`IMPROVE_TARGET` degrees,
//! 15 by default). Prints, per file, before and after: tets, the smallest
//! dihedral angle, tets under 1° and under 10°, tets flat in floating point,
//! the largest non-orthogonality and skewness; then the operations and the
//! seconds for meshing and improvement; or the error. The OBJ is read
//! naively (`v` and triangular `f` lines); use the welded OBJ files the
//! Thingi10k harness writes.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::sync::Arc;
use std::time::Instant;

use mesher_geom::Point3;
use mesher_tet::cdt::{mesh_region, round_region};
use mesher_tet::{
    ConformingParams, ImproveParams, SizeField, SizingParams, conforming_mesh, dihedral_angles,
    fv_quality, improve,
};

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
        "{} tets, min {:.3}, {under_one} <1, {under_ten} <10, {} flat, nonorth {:.2}, skew {:.3}",
        tets.len(),
        q.min_dihedral_degrees,
        q.flat_in_float,
        q.max_non_orthogonality_degrees,
        q.max_skewness
    )
}

fn main() {
    let exact = std::env::var("IMPROVE_KERNEL").is_ok_and(|s| s == "exact");
    let target = std::env::var("IMPROVE_TARGET")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(15.0);
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
        let started = Instant::now();
        let mesh = if exact {
            mesh_region(&points, &triangles, 0, 5)
                .map_err(|e| e.to_string())
                .and_then(|r| round_region(r, &points, &triangles).map_err(|e| e.to_string()))
                .map(|m| (m.points, m.tets))
        } else {
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
            let field = SizeField::from_surface(&points, &triangles, &SizingParams::for_size(h));
            let params = ConformingParams {
                size: (!field.is_uniform()).then(|| Arc::new(field)),
                ..ConformingParams::with_size(h)
            };
            conforming_mesh(&points, &triangles, &params)
                .map_err(|e| e.to_string())
                .map(|m| (m.points, m.tets))
        };
        let t_mesh = started.elapsed().as_secs_f64();
        let (p, t) = match mesh {
            Ok(m) => m,
            Err(e) => {
                println!("{path}\tfailed\t{e}\t{t_mesh:.2}");
                continue;
            }
        };
        let before = describe(&p, &t);
        let started = Instant::now();
        let r = improve(
            &p,
            &t,
            &ImproveParams {
                target_degrees: target,
                collapses: std::env::var("IMPROVE_COLLAPSES").map_or(true, |s| s != "0"),
                splits: std::env::var("IMPROVE_SPLITS").map_or(true, |s| s != "0"),
                ..ImproveParams::default()
            },
        );
        let t_improve = started.elapsed().as_secs_f64();
        let after = describe(&r.points, &r.tets);
        let s = r.stats;
        println!(
            "{path}\tok\t{before}\t{after}\t{} flips, {} edge removals, {} moves, {} collapses, {} splits, {} passes, poor {} -> {}, min {:.4} -> {:.4}\t{t_mesh:.2}\t{t_improve:.2}",
            s.flips,
            s.edge_removals,
            s.moves,
            s.collapses,
            s.splits,
            s.passes,
            s.poor_before,
            s.poor_after,
            s.min_before,
            s.min_after
        );
    }
}
