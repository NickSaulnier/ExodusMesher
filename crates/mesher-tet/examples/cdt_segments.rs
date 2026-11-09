//! Segment recovery of the exact kernel on triangle OBJ files, for measuring
//! (Milestone 12):
//!
//! ```text
//! cargo run --release -p mesher-tet --example cdt_segments -- a.obj b.obj ...
//! ```
//!
//! Prints, per file: triangles, subsegments, Steiner points, seconds and how
//! many segments the first tetrahedralization lacked, or the
//! error. The OBJ is read naively (`v` and triangular `f` lines); use the
//! welded OBJ files the Thingi10k harness writes.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::time::Instant;

use mesher_tet::cdt::recover_segments;

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
        let result = recover_segments(&points, &triangles, 0, 5);
        let seconds = started.elapsed().as_secs_f64();
        match result {
            Ok(r) => println!(
                "{path}\tok\t{}\t{}\t{}\t{seconds:.2}\tinitially missing {} of {}, {} rounds",
                triangles.len(),
                r.subsegments.len(),
                r.steiner,
                r.initially_missing,
                mesher_tet::cdt::segments::segments_of(&triangles).len(),
                r.rounds
            ),
            Err(e) => println!("{path}\tfailed\t{}\t{e}\t{seconds:.2}", triangles.len()),
        }
    }
}
