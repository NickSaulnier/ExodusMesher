//! Profile the constrained kernel on one triangle OBJ, as the pipeline runs
//! it (element size 0.03 of the diagonal, graded by the size field):
//!
//! ```text
//! cargo run --release -p mesher-tet --features profile --example profile -- model.obj [--uniform] [--remesh] [--size 0.03] [--max-points N]
//! ```
//!
//! Prints the outcome and the time spent in each timed function (inclusive:
//! nested scopes count in their callers too). The OBJ is read naively (`v`
//! and triangular `f` lines only); use the welded OBJ the Thingi10k harness
//! writes.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::sync::Arc;
use std::time::Instant;

use mesher_tet::{
    ConformingParams, RemeshParams, SizeField, SizingParams, conforming_mesh, prof, remesh,
};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let path = args
        .first()
        .expect("usage: profile model.obj [--uniform] [--remesh] [--size F]");
    let flag = |f: &str| args.iter().any(|a| a == f);
    let size: f64 = args
        .iter()
        .position(|a| a == "--size")
        .and_then(|i| args.get(i + 1))
        .map_or(0.03, |s| s.parse().unwrap());
    // Stop early (as a budget failure) to profile models that run long.
    let max_points: usize = args
        .iter()
        .position(|a| a == "--max-points")
        .and_then(|i| args.get(i + 1))
        .map_or(0, |s| s.parse().unwrap());

    let text = std::fs::read_to_string(path).unwrap();
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
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let diag = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
    let h = size * diag;
    eprintln!(
        "{path}: {} points, {} triangles, h = {h:.4e}",
        points.len(),
        triangles.len()
    );

    let started = Instant::now();
    let field = if flag("--uniform") {
        SizeField::uniform(h)
    } else {
        SizeField::from_surface(&points, &triangles, &SizingParams::for_size(h))
    };
    eprintln!("size field: {:?} in {:?}", field, started.elapsed());
    let mut groups = Vec::new();
    if flag("--remesh") {
        let t = Instant::now();
        let r = remesh(
            &points,
            &triangles,
            &[],
            &|p| field.at(p),
            &RemeshParams::default(),
        )
        .unwrap();
        eprintln!(
            "remeshed to {} triangles in {:?}",
            r.triangles.len(),
            t.elapsed()
        );
        groups = (0..r.triangles.len() as u32).collect();
        points = r.points;
        triangles = r.triangles;
    }
    let params = ConformingParams {
        size: (!field.is_uniform()).then(|| Arc::new(field)),
        groups,
        max_points,
        ..ConformingParams::with_size(h)
    };
    let t = Instant::now();
    let result = conforming_mesh(&points, &triangles, &params);
    let elapsed = t.elapsed();
    match &result {
        Ok(m) => eprintln!(
            "meshed: {} tets, {} points (+{} boundary, +{} interior Steiner) in {elapsed:?}",
            m.tets.len(),
            m.points.len(),
            m.boundary_steiner,
            m.interior_steiner
        ),
        Err(e) => eprintln!("failed after {elapsed:?}: {e}"),
    }
    eprint!("{}", prof::report());
}
