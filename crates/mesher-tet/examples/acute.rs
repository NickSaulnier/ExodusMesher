//! Acute vertices and edges of triangle OBJ files, for measuring
//! (Milestone 14):
//!
//! ```text
//! cargo run --release -p mesher-tet --example acute -- a.obj b.obj ...
//! ```
//!
//! Prints, per file: vertices, acute vertices, edges, acute edges, flat
//! edges, the smallest local feature size relative to the bounding-box
//! diagonal, and seconds; or the error. The OBJ is read naively (`v` and
//! triangular `f` lines); use the welded OBJ files the Thingi10k harness
//! writes.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::time::Instant;

use mesher_tet::chamfer::classify;

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
        let started = Instant::now();
        match classify(&points, &triangles) {
            Ok(a) => {
                let seconds = started.elapsed().as_secs_f64();
                let lfs = a
                    .vertex_lfs
                    .iter()
                    .copied()
                    .chain(a.edges.values().map(|e| e.lfs))
                    .fold(f64::INFINITY, f64::min);
                println!(
                    "{path}\t{}\t{}\t{}\t{}\t{}\t{:.1e}\t{seconds:.2}",
                    points.len(),
                    a.acute_vertices().count(),
                    a.edges.len(),
                    a.acute_edges().count(),
                    a.edges.values().filter(|e| e.flat).count(),
                    lfs / diagonal,
                );
            }
            Err(e) => println!("{path}\tfailed\t{e}"),
        }
    }
}
