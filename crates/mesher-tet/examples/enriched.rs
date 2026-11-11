//! Chamfer-and-refine end to end, for measuring (Milestone 14):
//!
//! ```text
//! cargo run --release -p mesher-tet --example enriched -- a.obj b.obj ...
//! ```
//!
//! Filters (by `CHAMFER_TAU` times the diagonal, 1e-4 by default) and
//! chamfers the input, refines it (radius-edge bound 2), makes the enriched
//! CDT with the input and rounds it to f64. Prints, per file: triangles;
//! refined points, those kept and those left out (balls, strips, near
//! filtered triangles, split points on bridges); the kernel's Steiner points on segments and faces;
//! tets inside; tets rounding turned over or left flat; the smallest face
//! angle and dihedral angle in degrees (and how many tets have a face angle
//! under 1°); seconds for chamfering, refinement, the CDT and rounding; or
//! the error. `REFINE_MAX_POINTS` caps refinement (default 2,000,000).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::time::Instant;

use mesher_tet::cdt::round_region;
use mesher_tet::chamfer::chamfer_filtered;
use mesher_tet::refine::{enriched_cdt, refine_conforming};

fn angle(p: [f64; 3], q: [f64; 3], r: [f64; 3]) -> f64 {
    let u = [0, 1, 2].map(|k| q[k] - p[k]);
    let w = [0, 1, 2].map(|k| r[k] - p[k]);
    let d = (0..3).map(|k| u[k] * w[k]).sum::<f64>();
    let n = |x: [f64; 3]| (0..3).map(|k| x[k] * x[k]).sum::<f64>().sqrt();
    libm::acos((d / (n(u) * n(w))).clamp(-1.0, 1.0)).to_degrees()
}

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
        let tau = std::env::var("CHAMFER_TAU")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(1e-4)
            * diagonal;
        let max_points = std::env::var("REFINE_MAX_POINTS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(2_000_000);
        let n = triangles.len();
        let started = Instant::now();
        let surface = match chamfer_filtered(&points, &triangles, tau) {
            Ok(s) => s,
            Err(e) => {
                println!("{path}\tfailed\t{n}\tchamfer: {e}");
                continue;
            }
        };
        let t_chamfer = started.elapsed().as_secs_f64();
        let started = Instant::now();
        let refined = match refine_conforming(&surface, &points, &triangles, 0.2, max_points, 2.0) {
            Ok(r) => r,
            Err(e) => {
                println!("{path}\tfailed\t{n}\trefinement: {e}");
                continue;
            }
        };
        let t_refine = started.elapsed().as_secs_f64();
        let started = Instant::now();
        let enriched = match enriched_cdt(&refined, &surface, &points, &triangles, tau, 0, 7) {
            Ok(e) => e,
            Err(e) => {
                println!(
                    "{path}\tfailed\t{n}\tenriched CDT: {e}\t{:.2}",
                    started.elapsed().as_secs_f64()
                );
                continue;
            }
        };
        let t_cdt = started.elapsed().as_secs_f64();
        let recovered = &enriched.region.recovered;
        let (steiner, steiner_faces) = (recovered.steiner, recovered.steiner_on_faces);
        let inside = enriched.region.tets.len();
        let started = Instant::now();
        let rounded = match round_region(enriched.region, &points, &triangles) {
            Ok(r) => r,
            Err(e) => {
                println!("{path}\tfailed\t{n}\trounding: {e}");
                continue;
            }
        };
        let t_round = started.elapsed().as_secs_f64();
        let at = |v: u32| rounded.points[v as usize];
        let mut min_face = f64::INFINITY;
        let mut under_one = 0usize;
        for t in &rounded.tets {
            let p = t.map(at);
            let mut worst = f64::INFINITY;
            for f in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
                let [a, b, c] = f.map(|k| p[k]);
                worst = worst
                    .min(angle(a, b, c))
                    .min(angle(b, c, a))
                    .min(angle(c, a, b));
            }
            if worst < 1.0 {
                under_one += 1;
            }
            min_face = min_face.min(worst);
        }
        let (min_dihedral, _) = mesher_tet::quality::dihedral_range(&rounded.points, &rounded.tets)
            .unwrap_or((0.0, 0.0));
        println!(
            "{path}\tok\t{n}\t{} refined, {} kept, {:?} left out\t{steiner}+{steiner_faces} Steiner\t{inside} tets\t{}+{} turned or flat\tface {min_face:.3}\tdihedral {min_dihedral:.3}\t{under_one} under 1°\t{t_chamfer:.2}\t{t_refine:.2}\t{t_cdt:.2}\t{t_round:.2}",
            refined.dt.points().len(),
            enriched.kept,
            enriched.cleared,
            rounded.inverted_before,
            rounded.flat_before,
        );
        let report = mesher_tet::prof::report();
        if !report.is_empty() {
            eprint!("{report}");
        }
    }
}
