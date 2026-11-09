//! Profile the constrained kernel on one triangle OBJ, as the pipeline runs
//! it (element size 0.03 of the diagonal, graded by the size field):
//!
//! ```text
//! cargo run --release -p mesher-tet --features profile --example profile -- model.obj [--uniform] [--remesh] [--size 0.03] [--max-points N]
//! ```
//!
//! Prints the outcome, where the slivers are, every tet flat in floating
//! point (with its vertices when `FLAT_DETAIL` is set), and the time spent in
//! each timed function (inclusive: nested scopes count in their callers too). The OBJ is read naively (`v`
//! and triangular `f` lines only); use the welded OBJ the Thingi10k harness
//! writes.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::sync::Arc;
use std::time::Instant;

use mesher_tet::{
    ConformingMesh, ConformingParams, RemeshParams, SizeField, SizingParams, conforming_mesh,
    dihedral_angles, prof, remesh,
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
    // As the pipeline does: a remeshed surface keeps its side sets.
    let mut side_sets = Vec::new();
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
        side_sets = vec![0; r.triangles.len()];
        points = r.points;
        triangles = r.triangles;
    }
    let params = ConformingParams {
        size: (!field.is_uniform()).then(|| Arc::new(field)),
        groups,
        side_sets,
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
    if let Ok(m) = &result {
        slivers(m, points.len());
    }
    eprint!("{}", prof::report());
}

/// Where the slivers are: tets below 10° (and 1°) by how many of their faces
/// are on the boundary, how many vertices are, and how many are input
/// vertices (indices below `inputs`).
fn slivers(m: &ConformingMesh, inputs: usize) {
    let key = |mut f: [u32; 3]| {
        f.sort_unstable();
        f
    };
    let faces: std::collections::HashSet<[u32; 3]> =
        m.boundary.iter().map(|&(f, _)| key(f)).collect();
    let on: std::collections::HashSet<u32> = m.boundary.iter().flat_map(|&(f, _)| f).collect();
    let mut table: std::collections::BTreeMap<(usize, usize, usize), (usize, usize)> =
        std::collections::BTreeMap::new();
    for t in &m.tets {
        let min = dihedral_angles(t.map(|v| m.points[v as usize]))
            .into_iter()
            .fold(180.0f64, f64::min);
        if min >= 10.0 {
            continue;
        }
        let bf = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
            .iter()
            .filter(|f| faces.contains(&key(f.map(|k| t[k]))))
            .count();
        let bv = t.iter().filter(|v| on.contains(v)).count();
        let iv = t.iter().filter(|&&v| (v as usize) < inputs).count();
        let e = table.entry((bf, bv, iv)).or_default();
        e.0 += 1;
        if min < 1.0 {
            e.1 += 1;
        }
    }
    eprintln!(
        "slivers by (boundary faces, boundary vertices, input vertices): count below 10°, below 1°"
    );
    for ((bf, bv, iv), (n, n1)) in table {
        eprintln!("  faces {bf} vertices {bv} input {iv}: {n:>6} {n1:>6}");
    }
    // Tets positive exactly but flat in floating point, which the volume
    // gate rejects.
    for t in &m.tets {
        let p = t.map(|v| m.points[v as usize]);
        let v = mesher_tet::quality::float_volume(p);
        if v > 0.0 {
            continue;
        }
        let bf = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
            .iter()
            .filter(|f| faces.contains(&key(f.map(|k| t[k]))))
            .count();
        let bv = t.iter().filter(|v| on.contains(v)).count();
        let iv = t.iter().filter(|&&v| (v as usize) < inputs).count();
        let longest = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            .iter()
            .map(|&(i, j)| {
                (0..3)
                    .map(|k| (p[i][k] - p[j][k]).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .fold(0.0f64, f64::max);
        let shortest = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            .iter()
            .map(|&(i, j)| {
                (0..3)
                    .map(|k| (p[i][k] - p[j][k]).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .fold(f64::INFINITY, f64::min);
        eprintln!(
            "  flat in float: tet {t:?} volume {v:e}, boundary faces {bf}, boundary vertices {bv}, input vertices {iv}, edges {shortest:.3e}..{longest:.3e}"
        );
        if std::env::var_os("FLAT_DETAIL").is_some() {
            for &w in t {
                let around = m.tets.iter().filter(|u| u.contains(&w)).count();
                eprintln!(
                    "    vertex {w}: {:?} on boundary {} input {} in {around} tets",
                    m.points[w as usize],
                    on.contains(&w),
                    (w as usize) < inputs
                );
            }
        }
    }
}
