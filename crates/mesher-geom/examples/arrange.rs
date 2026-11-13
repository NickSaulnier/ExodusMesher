//! The intersections of OBJ files' triangle soups (Milestone 16, stage 1):
//!
//! ```text
//! cargo run --release -p mesher-geom --example arrange -- a.obj b.obj ...
//! ```
//!
//! Prints a line per file: triangles (and zero-area ones left out),
//! intersecting pairs (coplanar ones), new points, segments, curve ends
//! (points where an odd number of segments meet: none when every
//! intersection curve closes), and seconds; or the error. `ARRANGE_DUMP=1`
//! lists each curve end with the triangles holding it, on stderr.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeMap;
use std::time::Instant;

use mesher_geom::Point3;
use mesher_geom::arrangement::intersections;
use mesher_geom::predicates::plane_axis;

fn main() {
    println!(
        "file\ttriangles\tdegenerate\tpairs\tcoplanar\tnew_points\tsegments\topen_ends\tseconds"
    );
    for path in std::env::args().skip(1) {
        let text = std::fs::read_to_string(&path).unwrap();
        let mut points: Vec<Point3> = Vec::new();
        let mut triangles: Vec<[u32; 3]> = Vec::new();
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
        let total = triangles.len();
        triangles.retain(|t| {
            let [a, b, c] = t.map(|v| points[v as usize]);
            plane_axis(a, b, c).is_some()
        });
        let degenerate = total - triangles.len();
        let started = Instant::now();
        match intersections(&points, &triangles) {
            Ok(x) => {
                let seconds = started.elapsed().as_secs_f64();
                let mut segments: Vec<[u32; 2]> = x
                    .constraints
                    .iter()
                    .flat_map(|c| c.segments.iter().copied())
                    .collect();
                segments.sort_unstable();
                segments.dedup();
                let mut degree: BTreeMap<u32, usize> = BTreeMap::new();
                for s in &segments {
                    for &v in s {
                        *degree.entry(v).or_insert(0) += 1;
                    }
                }
                let open = degree.values().filter(|&&d| d % 2 == 1).count();
                if std::env::var("ARRANGE_DUMP").is_ok() {
                    for (&v, &d) in &degree {
                        if d % 2 == 0 {
                            continue;
                        }
                        let p = &x.points[v as usize];
                        eprintln!("end {v} degree {d} at {:?} {:?}", p.approx(), p.kind());
                        for (t, c) in x.constraints.iter().enumerate() {
                            let segs: Vec<_> =
                                c.segments.iter().filter(|s| s.contains(&v)).collect();
                            if !segs.is_empty() || c.points.contains(&v) {
                                eprintln!("   triangle {t} {:?}: segments {segs:?}", triangles[t]);
                            }
                        }
                    }
                }
                println!(
                    "{path}\t{}\t{degenerate}\t{}\t{}\t{}\t{}\t{open}\t{seconds:.3}",
                    triangles.len(),
                    x.pairs,
                    x.coplanar_pairs,
                    x.points.len() - points.len(),
                    segments.len(),
                );
            }
            Err(e) => println!("{path}\tfailed\t{e}"),
        }
    }
}
