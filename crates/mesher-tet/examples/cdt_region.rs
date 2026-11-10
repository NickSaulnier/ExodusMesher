//! Region marking by the exact kernel on triangle OBJ files, for measuring
//! (Milestone 12):
//!
//! ```text
//! cargo run --release -p mesher-tet --example cdt_region -- a.obj b.obj ...
//! ```
//!
//! Prints, per file: triangles, tets inside, boundary faces, Steiner points
//! on segments and faces, seconds, and the inside tets' volume against the
//! volume the surface encloses (relative difference); or the error. The OBJ
//! is read naively (`v` and triangular `f` lines); use the welded OBJ files
//! the Thingi10k harness writes.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::time::Instant;

use mesher_tet::cdt::mesh_region;

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
        let result = mesh_region(&points, &triangles, max_points, 5);
        let seconds = started.elapsed().as_secs_f64();
        let det = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
            a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
        };
        match result {
            Ok(r) => {
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
                let at = |v: u32| r.dt().points()[v as usize].approx();
                let volume: f64 = r
                    .tets
                    .iter()
                    .map(|t| mesher_tet::quality::float_volume(t.map(at)))
                    .sum();
                println!(
                    "{path}	ok	{}	{}	{}	{}	{}	{seconds:.2}	volume off by {:.1e}",
                    triangles.len(),
                    r.tets.len(),
                    r.boundary.len(),
                    r.recovered.steiner,
                    r.recovered.steiner_on_faces,
                    (volume - enclosed).abs() / enclosed,
                );
            }
            Err(e) => println!("{path}	failed	{}	{e}	{seconds:.2}", triangles.len()),
        }
        let report = mesher_tet::prof::report();
        if !report.is_empty() {
            eprint!("{report}");
        }
    }
}
