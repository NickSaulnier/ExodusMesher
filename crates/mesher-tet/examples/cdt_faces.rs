//! Segment and face recovery of the exact kernel on triangle OBJ files, for
//! measuring (Milestone 12):
//!
//! ```text
//! cargo run --release -p mesher-tet --example cdt_faces -- a.obj b.obj ...
//! ```
//!
//! Prints, per file: triangles, subfaces, Steiner points on segments,
//! triangles still crossed, seconds, the subfaces missing after segment
//! recovery, the triangles crossed after segment recovery, flips and the
//! first cavities, and the Steiner points on faces; or the error. The OBJ is read naively (`v` and triangular `f` lines); use the
//! welded OBJ files the Thingi10k harness writes.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::time::Instant;

use mesher_tet::cdt::recover_faces;

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
        // CDT_MAX_POINTS caps the vertex count (0: the default budget).
        let max_points = std::env::var("CDT_MAX_POINTS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let result = recover_faces(&points, &triangles, max_points, 5);
        let seconds = started.elapsed().as_secs_f64();
        match result {
            Ok(r) => println!(
                "{path}	ok	{}	{}	{}	{}	{seconds:.2}	{} subfaces missing after segments; triangles crossed {:?}; {} Steiner points on faces",
                triangles.len(),
                r.subfaces.len(),
                r.steiner,
                r.missing.len(),
                r.initially_missing,
                r.crossed,
                r.steiner_on_faces,
            ),
            Err(e) => println!("{path}	failed	{}	{e}	{seconds:.2}", triangles.len()),
        }
        let report = mesher_tet::prof::report();
        if !report.is_empty() {
            eprint!("{report}");
        }
    }
}
