//! Conforming Delaunay refinement of the chamfered surface, for measuring
//! (Milestone 14):
//!
//! ```text
//! cargo run --release -p mesher-tet --example refine -- a.obj b.obj ...
//! ```
//!
//! Filters (by `CHAMFER_TAU` times the diagonal, 1e-4 by default) and
//! chamfers the input, then recovers every segment and planar face (and the
//! box edges and sides; `REFINE_STAGE=segments` stops after the segments) of
//! the chamfered surface
//! and the outer box's edges in a Delaunay tetrahedralization. Prints, per
//! file: triangles, chamfered vertices in use, subsegments, Steiner points
//! on segments, seconds for chamfering and for recovery, the shortest
//! subsegment over the diagonal, subfaces, Steiner points on faces and
//! subfaces set aside; then tets, Steiner points in tets and tets set aside
//! (`REFINE_STAGE=faces` stops after the faces); or the error.
//! `REFINE_MAX_POINTS` caps the vertex count (default 2,000,000).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::time::Instant;

use mesher_tet::chamfer::chamfer_filtered;
use mesher_tet::refine::{
    recover_faces_conforming, recover_segments_conforming, refine_conforming,
};

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
        let tau = std::env::var("CHAMFER_TAU")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(1e-4)
            * diagonal;
        let max_points = std::env::var("REFINE_MAX_POINTS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(2_000_000);
        let started = Instant::now();
        let surface = match chamfer_filtered(&points, &triangles, tau) {
            Ok(s) => s,
            Err(e) => {
                println!("{path}\tfailed\tchamfer: {e}");
                continue;
            }
        };
        let chamfer_seconds = started.elapsed().as_secs_f64();
        let used = surface
            .chamfered
            .faces
            .iter()
            .flat_map(|f| f.polygon.iter().copied())
            .collect::<std::collections::BTreeSet<u32>>()
            .len();
        let started = Instant::now();
        // REFINE_STAGE: segments, faces or tets (the default).
        let stage = std::env::var("REFINE_STAGE").unwrap_or_else(|_| "tets".to_owned());
        let result = match stage.as_str() {
            "segments" => {
                recover_segments_conforming(&surface, &points, &triangles, 0.2, max_points)
            }
            "faces" => recover_faces_conforming(&surface, &points, &triangles, 0.2, max_points),
            _ => refine_conforming(&surface, &points, &triangles, 0.2, max_points, 2.0),
        };
        match result {
            Ok(r) => {
                let at = |v: u32| r.dt.points()[v as usize].approx();
                let length = |s: &[u32; 2]| {
                    let (p, q) = (at(s[0]), at(s[1]));
                    (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>().sqrt()
                };
                let shortest = r
                    .subsegments
                    .iter()
                    .map(length)
                    .fold(f64::INFINITY, f64::min);
                println!(
                    "{path}\tok\t{}\t{used}\t{}\t{}\t{chamfer_seconds:.2}\t{:.2}\tshortest {:.1e}\t{} subfaces, {} Steiner points on faces, {} set aside	{} tets, {} Steiner points in tets, {} set aside",
                    triangles.len(),
                    r.subsegments.len(),
                    r.steiner_segments,
                    started.elapsed().as_secs_f64(),
                    shortest / diagonal,
                    r.subfaces.len(),
                    r.steiner_faces,
                    r.stuck_subfaces,
                    r.dt.tets().len(),
                    r.steiner_tets,
                    r.stuck_tets
                );
            }
            Err(e) => println!(
                "{path}\tfailed\t{}\trefinement: {e}\t{:.2}",
                triangles.len(),
                started.elapsed().as_secs_f64()
            ),
        }
        let report = mesher_tet::prof::report();
        if !report.is_empty() {
            eprint!("{report}");
        }
    }
}
