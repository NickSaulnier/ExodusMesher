//! Wraps of triangle OBJ files, for measuring (Milestone 17):
//!
//! ```text
//! cargo run --release -p mesher-tet --example wrap -- a.obj b.obj ...
//! ```
//!
//! α is `WRAP_ALPHA` times the bounding-box diagonal (default 0.03, the
//! default element size), δ is α over `WRAP_OFFSET` (default 30). Prints, per
//! file: input triangles, wrap triangles and points (on Voronoi edges,
//! projected), stuck gates, non-manifold edges and vertices found and cells
//! relabelled to remove them, the wrap's volume over the input's signed
//! volume, the largest distance from an input vertex to the wrap over the
//! diagonal, seconds, and `valid` when every edge has two triangles, no two
//! cross and every input vertex is inside (else the counts); then whether
//! the input is open (an edge not used as often in each direction), the
//! winding guard's cap points and blocked cells (`WRAP_GUARD=1` turns it
//! on), and, on a 32³ grid over the input's box, the solid the input's
//! winding number encloses (magnitude at least 1/2) against the wrap's
//! inside: their intersection over union, the fraction of the solid kept,
//! and the wrap's extra over the solid; or the error. `WRAP_OUT=<dir>`
//! writes each wrap as `<dir>/<name>.wrap.obj`. The OBJ is read naively (`v` and polygon `f`
//! lines, fanned).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::io::Write;
use std::time::Instant;

use mesher_geom::vec3::{cross, dot};
use mesher_geom::{Aabb, Bvh, FastWindingNumber, TriangleSoup, faces_intersect};
use mesher_tet::{WrapParams, wrap};

fn env(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn main() {
    let relative_alpha = env("WRAP_ALPHA", 0.03);
    let ratio = env("WRAP_OFFSET", 30.0);
    let out_dir = std::env::var("WRAP_OUT").ok();
    println!(
        "file\ttriangles\twrap\tpoints\tedge\tprojected\tstuck\tnonmanifold\trelabelled\tvolume_ratio\tdistance\tseconds\tcheck\topen\tcap\tblocked\tiou\tkept\textra"
    );
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
                    for k in 1..v.len().saturating_sub(1) {
                        triangles.push([v[0], v[k], v[k + 1]]);
                    }
                }
                _ => {}
            }
        }
        let diagonal = Aabb::of_points(&points).diagonal();
        let alpha = relative_alpha * diagonal;
        let params = WrapParams {
            alpha,
            offset: alpha / ratio,
            max_points: 5_000_000,
            winding_guard: std::env::var("WRAP_GUARD").is_ok_and(|v| v == "1"),
        };
        let started = Instant::now();
        match wrap(&points, &triangles, params) {
            Ok(w) => {
                let seconds = started.elapsed().as_secs_f64();
                let signed = |p: &[[f64; 3]], t: &[[u32; 3]]| -> f64 {
                    t.iter()
                        .map(|t| {
                            let [a, b, c] = t.map(|v| p[v as usize]);
                            dot(a, cross(b, c)) / 6.0
                        })
                        .sum()
                };
                let ratio_v = signed(&w.points, &w.triangles) / signed(&points, &triangles).abs();
                let bvh = Bvh::new(TriangleSoup::new(&w.points, &w.triangles).unwrap());
                let distance = triangles
                    .iter()
                    .flatten()
                    .map(|&v| bvh.nearest(points[v as usize]).unwrap().distance_squared)
                    .fold(0.0f64, f64::max)
                    .sqrt();
                // The guarantees: two triangles per edge, none crossing, the
                // input's vertices inside.
                let mut uses = std::collections::BTreeMap::new();
                for t in &w.triangles {
                    for k in 0..3 {
                        let (a, b) = (t[k], t[(k + 1) % 3]);
                        *uses.entry((a.min(b), a.max(b))).or_insert(0u32) += 1;
                    }
                }
                let bad_edges = uses.values().filter(|&&n| n != 2).count();
                let corners = |t: [u32; 3]| t.map(|v| w.points[v as usize]);
                let crossing = bvh
                    .self_overlaps()
                    .into_iter()
                    .filter(|&(i, j)| {
                        let (a, b) = (w.triangles[i as usize], w.triangles[j as usize]);
                        faces_intersect(a, corners(a), b, corners(b))
                    })
                    .count();
                let fwn =
                    FastWindingNumber::new(TriangleSoup::new(&w.points, &w.triangles).unwrap());
                let outside = triangles
                    .iter()
                    .flatten()
                    .filter(|&&v| (fwn.winding_number(points[v as usize]) - 1.0).abs() > 1e-3)
                    .count();
                let check = if bad_edges + crossing + outside == 0 {
                    "valid".to_string()
                } else {
                    format!("edges {bad_edges} crossing {crossing} outside {outside}")
                };
                // Open: some edge not used as often in each direction.
                let mut balance = std::collections::BTreeMap::new();
                for t in &triangles {
                    for k in 0..3 {
                        let (a, b) = (t[k], t[(k + 1) % 3]);
                        *balance.entry((a.min(b), a.max(b))).or_insert(0i64) +=
                            if a < b { 1 } else { -1 };
                    }
                }
                let open = balance.values().any(|&n| n != 0);
                // The winding-number solid against the wrap's inside.
                let input_fwn =
                    FastWindingNumber::new(TriangleSoup::new(&points, &triangles).unwrap());
                let b = Aabb::of_points(&points);
                let n = 32;
                let (mut both, mut solid, mut wrapped) = (0usize, 0usize, 0usize);
                for i in 0..n {
                    for j in 0..n {
                        for k in 0..n {
                            let q = [i, j, k].map(|x| (x as f64 + 0.5) / n as f64);
                            let p = [0, 1, 2].map(|a| b.min[a] + q[a] * (b.max[a] - b.min[a]));
                            let s = input_fwn.winding_number(p).abs() >= 0.5;
                            let w = fwn.winding_number(p) >= 0.5;
                            solid += usize::from(s);
                            wrapped += usize::from(w);
                            both += usize::from(s && w);
                        }
                    }
                }
                let union = solid + wrapped - both;
                let iou = both as f64 / union.max(1) as f64;
                let kept = both as f64 / solid.max(1) as f64;
                let extra = (wrapped - both) as f64 / solid.max(1) as f64;
                println!(
                    "{path}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{ratio_v:.4}\t{:.4e}\t{seconds:.3}\t{check}\t{open}\t{}\t{}\t{iou:.3}\t{kept:.3}\t{extra:.3}",
                    triangles.len(),
                    w.triangles.len(),
                    w.points.len(),
                    w.edge_points,
                    w.projected_points,
                    w.stuck,
                    w.nonmanifold_edges + w.nonmanifold_vertices,
                    w.relabelled,
                    distance / diagonal,
                    w.cap_points,
                    w.blocked
                );
                if let Some(dir) = &out_dir {
                    let name = std::path::Path::new(&path)
                        .file_stem()
                        .unwrap()
                        .to_string_lossy()
                        .into_owned();
                    let mut f = std::fs::File::create(
                        std::path::Path::new(dir).join(format!("{name}.wrap.obj")),
                    )
                    .unwrap();
                    for p in &w.points {
                        writeln!(f, "v {} {} {}", p[0], p[1], p[2]).unwrap();
                    }
                    for t in &w.triangles {
                        writeln!(f, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1).unwrap();
                    }
                }
            }
            Err(e) => println!("{path}\terror\t{e}"),
        }
    }
}
