//! Chamfering triangle OBJ files, for measuring (Milestone 14):
//!
//! ```text
//! cargo run --release -p mesher-tet --example chamfer -- a.obj b.obj ...
//! ```
//!
//! Prints, per file: triangles, acute vertices and edges, vertices of the
//! chamfered surface, its faces, the smallest angle between two of its
//! segments meeting at a vertex (degrees, from rounded positions; flat input
//! edges are not segments), and seconds; or the error. The OBJ is read
//! naively (`v` and triangular `f` lines).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use mesher_tet::chamfer::{ChamferVertex, chamfer, classify};

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
        let acute = match classify(&points, &triangles) {
            Ok(a) => a,
            Err(e) => {
                println!("{path}\tfailed\t{e}");
                continue;
            }
        };
        let c = match chamfer(&points, &triangles, &acute) {
            Ok(c) => c,
            Err(e) => {
                println!("{path}\tfailed\t{e}");
                continue;
            }
        };
        let seconds = started.elapsed().as_secs_f64();
        // Segments: face edges, except those lying along a flat input edge.
        let on_edge = |v: u32| -> Vec<(u32, u32)> {
            match c.vertices[v as usize] {
                ChamferVertex::OnEdge { edge, .. } => vec![edge],
                ChamferVertex::Input(x) => acute
                    .edges
                    .keys()
                    .filter(|e| e.0 == x || e.1 == x)
                    .copied()
                    .collect(),
                ChamferVertex::InFace { .. } => Vec::new(),
            }
        };
        let mut segments: BTreeSet<(u32, u32)> = BTreeSet::new();
        for f in &c.faces {
            let k = f.polygon.len();
            for i in 0..k {
                let (u, v) = (f.polygon[i], f.polygon[(i + 1) % k]);
                let shared: Vec<(u32, u32)> = on_edge(u)
                    .into_iter()
                    .filter(|e| on_edge(v).contains(e))
                    .collect();
                if shared
                    .iter()
                    .any(|e| acute.edges.get(e).is_some_and(|i| i.flat))
                {
                    continue;
                }
                segments.insert((u.min(v), u.max(v)));
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
            "{path}\t{}\t{}\t{}\t{}\t{}\t{smallest:.3}\t{seconds:.2}",
            triangles.len(),
            acute.acute_vertices().count(),
            acute.acute_edges().count(),
            c.vertices.len(),
            c.faces.len(),
        );
    }
}
