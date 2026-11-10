//! Rounding the exact kernel's mesh to f64 on triangle OBJ files, for
//! measuring (Milestone 13):
//!
//! ```text
//! cargo run --release -p mesher-tet --example cdt_round -- a.obj b.obj ...
//! ```
//!
//! Prints, per file: triangles, tets, constructed points rounded, inside
//! tets rounding turned over and those flat in floating point, fixes (vertex moves, 2-3 flips, edge
//! removals), seconds for the exact mesh and for rounding, and the volume
//! of the rounded mesh against the volume the surface encloses; or the
//! error. The OBJ is read naively (`v` and triangular `f` lines); use the
//! welded OBJ files the Thingi10k harness writes.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::time::Instant;

use mesher_tet::cdt::{Untangle, mesh_region, round_region};

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
        let started = Instant::now();
        let region = match mesh_region(&points, &triangles, 0, 5) {
            Ok(r) => r,
            Err(e) => {
                println!("{path}\tfailed\t{}\texact: {e}", triangles.len());
                continue;
            }
        };
        let exact_seconds = started.elapsed().as_secs_f64();
        let started = Instant::now();
        let result = round_region(region, &points, &triangles);
        let round_seconds = started.elapsed().as_secs_f64();
        let det = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
            a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
        };
        match result {
            Ok(m) => {
                let enclosed: f64 = triangles
                    .iter()
                    .map(|t| {
                        det(
                            points[t[0] as usize],
                            points[t[1] as usize],
                            points[t[2] as usize],
                        ) / 6.0
                    })
                    .sum::<f64>()
                    .abs();
                let volume: f64 = m
                    .tets
                    .iter()
                    .map(|t| mesher_tet::quality::float_volume(t.map(|v| m.points[v as usize])))
                    .sum();
                let count =
                    |f: &dyn Fn(&Untangle) -> bool| m.untangles.iter().filter(|u| f(u)).count();
                println!(
                    "{path}\tok\t{}\t{}\t{}\t{}/{}\t{} moves, {} 2-3 flips, {} edge removals\t{exact_seconds:.2}\t{round_seconds:.2}\tvolume off by {:.1e}",
                    triangles.len(),
                    m.tets.len(),
                    m.rounded,
                    m.inverted_before,
                    m.flat_before,
                    count(&|u| matches!(u, Untangle::Move { .. })),
                    count(&|u| matches!(u, Untangle::Flip23 { .. })),
                    count(&|u| matches!(u, Untangle::RemoveEdge { .. })),
                    (volume - enclosed).abs() / enclosed,
                );
            }
            Err(e) => println!(
                "{path}\tfailed\t{}\trounding: {e}\t{exact_seconds:.2}\t{round_seconds:.2}",
                triangles.len()
            ),
        }
    }
}
