//! The arrangement of OBJ files' triangle soups (Milestone 16, stages 1 and
//! 2):
//!
//! ```text
//! cargo run --release -p mesher-geom --example arrange -- a.obj b.obj ...
//! ```
//!
//! Prints a line per file: triangles (and zero-area ones left out),
//! intersecting pairs (coplanar ones), pieces, unbalanced edges of the
//! input and of the arrangement (an edge is balanced when it is used as
//! often in each direction: every edge of a closed oriented surface is, and
//! the arrangement of one must keep that), the largest error in a
//! triangle's area covered by its pieces (over the squared diagonal), and
//! seconds; or the error. `ARRANGE_DUMP=1` lists, on stderr, the triangles
//! whose pieces miss more than 1e-12 of that. `ARRANGE_RESOLVE=1` adds a
//! line for stage 3 (closed inputs): faces, patches, faces dropped and
//! kept, unbalanced edges, the volume against the input's, and, rounded to
//! f64, the triangles left flat and the pairs that intersect again;
//! `ARRANGE_PROPAGATE=1` also finds the windings from one cast per piece
//! and says on stderr where that differs.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeMap;
use std::time::Instant;

use mesher_geom::Point3;
use mesher_geom::arrangement::{
    ArrangementError, arrangement, intersections, resolve, resolve_with,
};
use mesher_geom::predicates::plane_axis;

fn area(p: Point3, q: Point3, r: Point3) -> f64 {
    let u = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
    let v = [r[0] - p[0], r[1] - p[1], r[2] - p[2]];
    let c = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
}

fn unbalanced(triangles: &[[u32; 3]]) -> usize {
    let mut net: BTreeMap<(u32, u32), i64> = BTreeMap::new();
    for t in triangles {
        for k in 0..3 {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            *net.entry((u.min(v), u.max(v))).or_insert(0) += if u < v { 1 } else { -1 };
        }
    }
    net.values().filter(|&&n| n != 0).count()
}

fn main() {
    println!(
        "file\ttriangles\tdegenerate\tpairs\tcoplanar\tpieces\tunbalanced_in\tunbalanced_out\tarea_error\tseconds"
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
        match arrangement(&points, &triangles) {
            Ok(a) => {
                let seconds = started.elapsed().as_secs_f64();
                let mut covered = vec![0.0; triangles.len()];
                for (t, &o) in a.triangles.iter().zip(&a.origin) {
                    let [x, y, z] = t.map(|v| a.points[v as usize].approx());
                    covered[o as usize] += area(x, y, z);
                }
                // Against the square of the diagonal: a sliver's pieces,
                // measured at rounded positions, lose its tiny area to
                // cancellation.
                let (lo, hi) = points.iter().fold(
                    ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
                    |(lo, hi), p| {
                        (
                            [0, 1, 2].map(|k| lo[k].min(p[k])),
                            [0, 1, 2].map(|k| hi[k].max(p[k])),
                        )
                    },
                );
                let d2: f64 = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum();
                let error = triangles
                    .iter()
                    .zip(&covered)
                    .map(|(t, c)| {
                        let [x, y, z] = t.map(|v| points[v as usize]);
                        (c - area(x, y, z)).abs() / d2
                    })
                    .fold(0.0, f64::max);
                if std::env::var("ARRANGE_DUMP").is_ok() {
                    for (o, (t, c)) in triangles.iter().zip(&covered).enumerate() {
                        let [x, y, z] = t.map(|v| points[v as usize]);
                        let want = area(x, y, z);
                        if (c - want).abs() / d2 > 1e-12 {
                            let pieces: Vec<[u32; 3]> = a
                                .triangles
                                .iter()
                                .zip(&a.origin)
                                .filter(|(_, g)| **g as usize == o)
                                .map(|(t, _)| *t)
                                .collect();
                            eprintln!(
                                "triangle {o} {t:?} at {x:?} {y:?} {z:?}: area {want:.3e}, covered {c:.3e}, pieces {pieces:?}"
                            );
                            let mut new: Vec<u32> = pieces
                                .iter()
                                .flatten()
                                .copied()
                                .filter(|&v| v as usize >= points.len())
                                .collect();
                            new.sort_unstable();
                            new.dedup();
                            for v in new {
                                let q = &a.points[v as usize];
                                let b = q.bounds();
                                eprintln!(
                                    "   point {v} approx {:?} bounds width {:?} {:?}",
                                    q.approx(),
                                    b.map(|r| r[1] - r[0]),
                                    q.kind()
                                );
                            }
                        }
                    }
                }
                println!(
                    "{path}\t{}\t{degenerate}\t{}\t{}\t{}\t{}\t{}\t{error:.1e}\t{seconds:.3}",
                    triangles.len(),
                    a.pairs,
                    a.coplanar_pairs,
                    a.triangles.len(),
                    unbalanced(&triangles),
                    unbalanced(&a.triangles),
                );
            }
            Err(e) => println!("{path}\tfailed\t{e}"),
        }
        if std::env::var("ARRANGE_RESOLVE").is_err() {
            continue;
        }
        // Stage 3, and a preview of stage 4: the boundary rounded to f64.
        let started = Instant::now();
        match resolve(&points, &triangles) {
            Ok(r) => {
                let seconds = started.elapsed().as_secs_f64();
                let rp: Vec<Point3> = r.points.iter().map(|q| q.approx()).collect();
                let proper: Vec<[u32; 3]> = r
                    .triangles
                    .iter()
                    .copied()
                    .filter(|t| {
                        let [a, b, c] = t.map(|v| rp[v as usize]);
                        plane_axis(a, b, c).is_some()
                    })
                    .collect();
                let flat = r.triangles.len() - proper.len();
                let crossing = intersections(&rp, &proper).map_or(usize::MAX, |x| x.pairs);
                if std::env::var("ARRANGE_PROPAGATE").is_ok() {
                    match resolve_with(&points, &triangles, false) {
                        Ok(by_edges) if by_edges.triangles == r.triangles => {}
                        Ok(_) => eprintln!("{path}: around edges gives another surface"),
                        Err(e) => eprintln!("{path}: around edges failed: {e}"),
                    }
                }
                let vol = |pts: &[Point3], tris: &[[u32; 3]]| -> f64 {
                    tris.iter()
                        .map(|t| {
                            let [a, b, c] = t.map(|v| pts[v as usize]);
                            (a[0] * (b[1] * c[2] - b[2] * c[1])
                                - a[1] * (b[0] * c[2] - b[2] * c[0])
                                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                                / 6.0
                        })
                        .sum()
                };
                println!(
                    "{path}\tresolved\tfaces {}\tpatches {}\tdropped {}\tkept {}\tunbalanced {}\tvolume {:.6e} from {:.6e}\trounded: flat {flat}, crossing pairs {crossing}\t{seconds:.3}",
                    r.faces,
                    r.patches,
                    r.dropped,
                    r.triangles.len(),
                    unbalanced(&r.triangles),
                    vol(&rp, &r.triangles),
                    vol(&points, &triangles),
                );
            }
            Err(ArrangementError::Open { edges }) => {
                println!("{path}\topen\t{edges} unbalanced edges");
            }
            Err(e) => println!("{path}\tresolve failed\t{e}"),
        }
    }
}
