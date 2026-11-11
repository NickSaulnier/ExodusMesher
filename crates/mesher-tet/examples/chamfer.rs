//! Chamfering triangle OBJ files, for measuring (Milestone 14):
//!
//! ```text
//! cargo run --release -p mesher-tet --example chamfer -- a.obj b.obj ...
//! ```
//!
//! Filters the input by `CHAMFER_TAU` times the bounding-box diagonal
//! (1e-4 by default), chamfers it and simplifies its bridges, then prints,
//! per file: triangles, triangles filtered out, rounds of filtering and
//! chamfering, acute vertices and edges of the filtered surface, vertices of
//! the chamfered surface, the smallest angle between two of its segments
//! meeting at a vertex (degrees, from rounded positions; flat input edges
//! are not segments), and seconds; or the error. The OBJ is read
//! naively (`v` and triangular `f` lines).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use mesher_tet::chamfer::chamfer_filtered;

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
        // CHAMFER_TAU: the filtering tolerance relative to the diagonal.
        let tau = std::env::var("CHAMFER_TAU")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(1e-4)
            * diagonal;
        let started = Instant::now();
        let filtered = match chamfer_filtered(&points, &triangles, tau) {
            Ok(f) => f,
            Err(e) => {
                println!("{path}\tfailed\t{e}");
                continue;
            }
        };
        let seconds = started.elapsed().as_secs_f64();
        let (acute, c) = (&filtered.acute, &filtered.chamfered);
        let used = c
            .faces
            .iter()
            .flat_map(|f| f.polygon.iter().copied())
            .collect::<BTreeSet<u32>>()
            .len();
        // Segments: face edges, except those lying along a flat input edge.
        let mut segments: BTreeSet<(u32, u32)> = BTreeSet::new();
        for f in &c.faces {
            let k = f.polygon.len();
            for i in 0..k {
                let (u, v) = (f.polygon[i], f.polygon[(i + 1) % k]);
                if c.is_segment(u, v, acute) {
                    segments.insert((u.min(v), u.max(v)));
                }
            }
        }
        let mut around: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
        for &(u, v) in &segments {
            around.entry(u).or_default().push(v);
            around.entry(v).or_default().push(u);
        }
        let at = |v: u32| c.points[v as usize].approx();
        let mut smallest = 180.0f64;
        for (&v, near) in &around {
            if near.len() > 64 {
                continue;
            }
            let here = at(v);
            for (i, &a) in near.iter().enumerate() {
                for &b in &near[i + 1..] {
                    let (p, q) = (at(a), at(b));
                    let u = [0, 1, 2].map(|k| p[k] - here[k]);
                    let w = [0, 1, 2].map(|k| q[k] - here[k]);
                    let d: f64 = (0..3).map(|k| u[k] * w[k]).sum();
                    let nu: f64 = (0..3).map(|k| u[k] * u[k]).sum::<f64>().sqrt();
                    let nw: f64 = (0..3).map(|k| w[k] * w[k]).sum::<f64>().sqrt();
                    let angle = libm::acos((d / (nu * nw)).clamp(-1.0, 1.0)).to_degrees();
                    if angle.is_finite() {
                        smallest = smallest.min(angle);
                    }
                }
            }
        }
        println!(
            "{path}\t{}\t{}\t{}\t{}\t{}\t{used}\t{smallest:.3}\t{seconds:.2}",
            triangles.len(),
            filtered.kept.iter().filter(|k| !**k).count(),
            filtered.rounds,
            acute.acute_vertices().count(),
            acute.acute_edges().count(),
        );
    }
}
