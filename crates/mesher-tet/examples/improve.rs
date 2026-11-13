//! Worst-cell improvement on triangle OBJ files, for measuring (Milestone 15):
//!
//! ```text
//! cargo run --release -p mesher-tet --example improve -- a.obj b.obj ...
//! ```
//!
//! Meshes each file as the pipeline does (the conforming kernel at 0.03 of
//! the diagonal, graded by the size field; `IMPROVE_KERNEL=exact` for the
//! exact kernel rounded to f64), then improves it (`IMPROVE_TARGET` degrees,
//! 15 by default; `IMPROVE_NONORTH` and `IMPROVE_SKEW` for the face
//! targets, `IMPROVE_FV=0` to turn them off; `IMPROVE_SPLITS=0`,
//! `IMPROVE_COLLAPSES=0` and `IMPROVE_INSERT=0` to leave those operations
//! out). `IMPROVE_REMESH=1` remeshes the surface to the size field first,
//! as the pipeline's fallback does. Prints, per file, before and
//! after: tets, the smallest dihedral angle, tets under 1° and under 10°,
//! tets flat in floating point, the largest non-orthogonality and skewness,
//! and both over wall-adjacent cells; then the operations and the
//! seconds for meshing and improvement; or the error. `IMPROVE_WORST=1`
//! describes the worst non-orthogonal face after improvement on stderr
//! (its cells' boundary vertices, smallest angle, the free vertex's height
//! and the domain's thickness there). With `--features profile` the timing
//! scopes go to stderr. The OBJ is read
//! naively (`v` and triangular `f` lines); use the welded OBJ files the
//! Thingi10k harness writes.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::sync::Arc;
use std::time::Instant;

use mesher_geom::Point3;
use mesher_tet::cdt::{mesh_region, round_region};
use mesher_tet::{
    ConformingParams, ImproveParams, RemeshParams, SizeField, SizingParams, conforming_mesh,
    dihedral_angles, fv_quality, improve, remesh,
};

fn describe(points: &[Point3], tets: &[[u32; 4]]) -> String {
    let q = fv_quality(points, tets);
    let (mut under_one, mut under_ten) = (0, 0);
    for t in tets {
        let min = dihedral_angles(t.map(|v| points[v as usize]))
            .into_iter()
            .fold(180.0, f64::min);
        under_one += usize::from(min < 1.0);
        under_ten += usize::from(min < 10.0);
    }
    format!(
        "{} tets, min {:.3}, {under_one} <1, {under_ten} <10, {} flat, nonorth {:.2}, skew {:.3}, wall nonorth {:.2}, wall skew {:.3}",
        tets.len(),
        q.min_dihedral_degrees,
        q.flat_in_float,
        q.max_non_orthogonality_degrees,
        q.max_skewness,
        q.wall_max_non_orthogonality_degrees,
        q.wall_max_skewness
    )
}

/// The worst non-orthogonal face: its angle, whether it is on the
/// boundary, and per cell the vertices on the boundary and the smallest
/// dihedral angle.
fn worst_face(points: &[Point3], tets: &[[u32; 4]]) -> String {
    use std::collections::BTreeMap;
    const OUT: [[usize; 3]; 4] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];
    let mut faces: BTreeMap<[u32; 3], Vec<(usize, usize)>> = BTreeMap::new();
    for (i, t) in tets.iter().enumerate() {
        for (k, f) in OUT.iter().enumerate() {
            let mut key = f.map(|j| t[j]);
            key.sort_unstable();
            faces.entry(key).or_default().push((i, k));
        }
    }
    let mut on_wall = vec![false; points.len()];
    for (key, cells) in &faces {
        if cells.len() == 1 {
            for &v in key {
                on_wall[v as usize] = true;
            }
        }
    }
    let p = |t: usize| tets[t].map(|v| points[v as usize]);
    let c = |t: usize| {
        let q = p(t);
        [0, 1, 2].map(|k| (q[0][k] + q[1][k] + q[2][k] + q[3][k]) / 4.0)
    };
    let sub = |a: Point3, b: Point3| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let mut worst = (0.0, String::new());
    for (key, cells) in &faces {
        let (i, k) = cells[0];
        let f = OUT[k].map(|m| p(i)[m]);
        let fc = [0, 1, 2].map(|m| (f[0][m] + f[1][m] + f[2][m]) / 3.0);
        let (u, w) = (sub(f[1], f[0]), sub(f[2], f[0]));
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let d = if cells.len() == 2 {
            sub(c(cells[1].0), c(i))
        } else {
            sub(fc, c(i))
        };
        let dot = n[0] * d[0] + n[1] * d[1] + n[2] * d[2];
        let l =
            (n.iter().map(|x| x * x).sum::<f64>() * d.iter().map(|x| x * x).sum::<f64>()).sqrt();
        let angle = libm::acos((dot / l).clamp(-1.0, 1.0)).to_degrees();
        if angle > worst.0 {
            let describe = |t: usize| {
                let fixed = tets[t].iter().filter(|&&v| on_wall[v as usize]).count();
                let min = dihedral_angles(p(t)).into_iter().fold(180.0, f64::min);
                // The free vertex's height over the plane of the others, and
                // the domain's thickness there: along the normal of that
                // plane, the nearest boundary face on each side.
                let q = p(t);
                let Some(free) = (0..4).find(|&m| !on_wall[tets[t][m] as usize]) else {
                    // Every vertex on the wall: its boundary faces (by the
                    // opposite vertex), the cell's dihedral angle between
                    // each two of them (the domain's angle at that edge),
                    // and, across each other face, the opposite vertex's
                    // height over it, over its longest edge.
                    let opposite: Vec<usize> = (0..4)
                        .filter(|&m| {
                            let mut key = OUT[m].map(|j| tets[t][j]);
                            key.sort_unstable();
                            faces.get(&key).is_some_and(|c| c.len() == 1)
                        })
                        .collect();
                    let angles = dihedral_angles(p(t));
                    // dihedral_angles' order: edges 01, 02, 03, 12, 13, 23;
                    // the faces opposite a and b meet at the edge of the
                    // other two vertices.
                    let edge_index = |a: usize, b: usize| {
                        let mut e: Vec<usize> = (0..4).filter(|&m| m != a && m != b).collect();
                        e.sort_unstable();
                        match (e[0], e[1]) {
                            (0, 1) => 0,
                            (0, 2) => 1,
                            (0, 3) => 2,
                            (1, 2) => 3,
                            (1, 3) => 4,
                            _ => 5,
                        }
                    };
                    let mut between = Vec::new();
                    for (x, &a) in opposite.iter().enumerate() {
                        for &b in &opposite[x + 1..] {
                            between.push(format!("{:.1}°", angles[edge_index(a, b)]));
                        }
                    }
                    let q = p(t);
                    let heights: Vec<String> = (0..4)
                        .filter(|m| !opposite.contains(m))
                        .map(|m| {
                            let o: Vec<Point3> = OUT[m].iter().map(|&j| q[j]).collect();
                            let (u, w) = (sub(o[1], o[0]), sub(o[2], o[0]));
                            let n = [
                                u[1] * w[2] - u[2] * w[1],
                                u[2] * w[0] - u[0] * w[2],
                                u[0] * w[1] - u[1] * w[0],
                            ];
                            let nl = n.iter().map(|x| x * x).sum::<f64>().sqrt();
                            let d = sub(q[m], o[0]);
                            let h = (n[0] * d[0] + n[1] * d[1] + n[2] * d[2]).abs() / nl;
                            let e = (0..3)
                                .map(|k| {
                                    let x = sub(o[(k + 1) % 3], o[k]);
                                    (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt()
                                })
                                .fold(0.0, f64::max);
                            format!("{:.3}", h / e)
                        })
                        .collect();
                    // Each boundary face's smallest angle, and the cell's
                    // vertices.
                    let face_angles: Vec<String> = opposite
                        .iter()
                        .map(|&m| {
                            let o: Vec<Point3> = OUT[m].iter().map(|&j| q[j]).collect();
                            let angle = |a: Point3, b: Point3, c: Point3| {
                                let (u, w) = (sub(b, a), sub(c, a));
                                let d = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
                                let l = (u.iter().map(|x| x * x).sum::<f64>()
                                    * w.iter().map(|x| x * x).sum::<f64>())
                                .sqrt();
                                libm::acos((d / l).clamp(-1.0, 1.0)).to_degrees()
                            };
                            let a = angle(o[0], o[1], o[2])
                                .min(angle(o[1], o[2], o[0]))
                                .min(angle(o[2], o[0], o[1]));
                            format!("{a:.1}°")
                        })
                        .collect();
                    let corners: Vec<String> = q
                        .iter()
                        .map(|c| format!("({:.4}, {:.4}, {:.4})", c[0], c[1], c[2]))
                        .collect();
                    return format!(
                        "{fixed} wall vertices, min dihedral {min:.2}, {} boundary faces at {} (smallest angles {}), interior faces' h/e {}, at {}",
                        opposite.len(),
                        if between.is_empty() {
                            "-".to_owned()
                        } else {
                            between.join(" ")
                        },
                        face_angles.join(" "),
                        heights.join(" "),
                        corners.join(" ")
                    );
                };
                let o: Vec<Point3> = (0..4).filter(|&m| m != free).map(|m| q[m]).collect();
                let (u, w) = (sub(o[1], o[0]), sub(o[2], o[0]));
                let mut n = [
                    u[1] * w[2] - u[2] * w[1],
                    u[2] * w[0] - u[0] * w[2],
                    u[0] * w[1] - u[1] * w[0],
                ];
                let nl = n.iter().map(|x| x * x).sum::<f64>().sqrt();
                n = n.map(|x| x / nl);
                let dotp = |a: Point3, b: Point3| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                let height = dotp(sub(q[free], o[0]), n);
                let edge = (0..3)
                    .map(|m| dotp(sub(o[(m + 1) % 3], o[m]), sub(o[(m + 1) % 3], o[m])).sqrt())
                    .fold(0.0, f64::max);
                let oc = [0, 1, 2].map(|k| (o[0][k] + o[1][k] + o[2][k]) / 3.0);
                let mut thick = f64::INFINITY;
                for (key, cells) in &faces {
                    if cells.len() != 1 || key.iter().all(|v| tets[t].contains(v)) {
                        continue;
                    }
                    let f = key.map(|v| points[v as usize]);
                    // Ray oc + s n against triangle f (Moller-Trumbore).
                    for dir in [n, n.map(|x| -x)] {
                        let (e1, e2) = (sub(f[1], f[0]), sub(f[2], f[0]));
                        let pv = [
                            dir[1] * e2[2] - dir[2] * e2[1],
                            dir[2] * e2[0] - dir[0] * e2[2],
                            dir[0] * e2[1] - dir[1] * e2[0],
                        ];
                        let det = dotp(e1, pv);
                        if det.abs() < 1e-300 {
                            continue;
                        }
                        let tv = sub(oc, f[0]);
                        let a = dotp(tv, pv) / det;
                        let qv = [
                            tv[1] * e1[2] - tv[2] * e1[1],
                            tv[2] * e1[0] - tv[0] * e1[2],
                            tv[0] * e1[1] - tv[1] * e1[0],
                        ];
                        let b = dotp(dir, qv) / det;
                        let s = dotp(e2, qv) / det;
                        if a >= 0.0 && b >= 0.0 && a + b <= 1.0 && s > 0.0 {
                            thick = thick.min(s);
                        }
                    }
                }
                format!(
                    "{fixed} wall vertices, min dihedral {min:.2}, height/edge {:.4}, thickness/edge {:.4}",
                    height.abs() / edge,
                    thick / edge
                )
            };
            let mut cells_text: Vec<String> = cells.iter().map(|&(t, _)| describe(t)).collect();
            // The free vertex of the face, and its star: per tet, the
            // boundary vertices, the smallest angle and the free vertex's
            // height over the opposite face, over that face's longest edge.
            if let Some(&v) = key.iter().find(|&&v| !on_wall[v as usize]) {
                let mut star = String::new();
                for (t, tet) in tets.iter().enumerate() {
                    let Some(m) = tet.iter().position(|&w| w == v) else {
                        continue;
                    };
                    let q = p(t);
                    let o: Vec<Point3> = (0..4).filter(|&k| k != m).map(|k| q[k]).collect();
                    let (u, w) = (sub(o[1], o[0]), sub(o[2], o[0]));
                    let n = [
                        u[1] * w[2] - u[2] * w[1],
                        u[2] * w[0] - u[0] * w[2],
                        u[0] * w[1] - u[1] * w[0],
                    ];
                    let nl = n.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let d = sub(q[m], o[0]);
                    let h = (n[0] * d[0] + n[1] * d[1] + n[2] * d[2]).abs() / nl;
                    let e = (0..3)
                        .map(|k| {
                            let x = sub(o[(k + 1) % 3], o[k]);
                            (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt()
                        })
                        .fold(0.0, f64::max);
                    let fixed = tet.iter().filter(|&&w| on_wall[w as usize]).count();
                    let min = dihedral_angles(q).into_iter().fold(180.0, f64::min);
                    star += &format!(" ({fixed}w {min:.1}° h/e {:.3})", h / e);
                }
                cells_text.push(format!("star of {v}:{star}"));
            }
            worst = (
                angle,
                format!(
                    "{angle:.2} {} face with {} wall vertices [{}]",
                    if cells.len() == 1 {
                        "boundary"
                    } else {
                        "interior"
                    },
                    key.iter().filter(|&&v| on_wall[v as usize]).count(),
                    cells_text.join("; ")
                ),
            );
        }
    }
    worst.1
}

fn main() {
    let exact = std::env::var("IMPROVE_KERNEL").is_ok_and(|s| s == "exact");
    let var = |name: &str, default: f64| {
        std::env::var(name)
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(default)
    };
    let defaults = ImproveParams::default();
    let target = var("IMPROVE_TARGET", defaults.target_degrees);
    let fv = std::env::var("IMPROVE_FV").map_or(true, |s| s != "0");
    let non_orthogonality = if fv {
        var("IMPROVE_NONORTH", defaults.max_non_orthogonality_degrees)
    } else {
        90.0
    };
    let skewness = if fv {
        var("IMPROVE_SKEW", defaults.max_skewness)
    } else {
        f64::INFINITY
    };
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
        let mesh = if exact {
            mesh_region(&points, &triangles, 0, 5)
                .map_err(|e| e.to_string())
                .and_then(|r| round_region(r, &points, &triangles).map_err(|e| e.to_string()))
                .map(|m| (m.points, m.tets))
        } else {
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
            let h = 0.03 * diagonal;
            // `IMPROVE_SIZE_FLOOR=k`: the size field's floor h/k (default 4);
            // `IMPROVE_BOUNDARY_FLOOR=k`: the kernel's h/k (default 16);
            // `IMPROVE_SMALL_FEATURE=f`: remeshing's small-feature merge
            // threshold, a fraction of the target size (default 0.1).
            let floor = |name: &str, default: f64| {
                std::env::var(name)
                    .ok()
                    .and_then(|s| s.parse::<f64>().ok())
                    .unwrap_or(default)
            };
            let sizing = SizingParams {
                min_size: h / floor("IMPROVE_SIZE_FLOOR", 4.0),
                ..SizingParams::for_size(h)
            };
            let field = SizeField::from_surface(&points, &triangles, &sizing);
            // `IMPROVE_REMESH=1`: the surface remeshed to the size field
            // first, each triangle its own facet, as the pipeline's
            // fallback does.
            let remeshed = std::env::var("IMPROVE_REMESH")
                .is_ok_and(|s| s == "1")
                .then(|| {
                    remesh(
                        &points,
                        &triangles,
                        &vec![0; triangles.len()],
                        &|q| field.at(q),
                        &RemeshParams {
                            small_feature: floor("IMPROVE_SMALL_FEATURE", 0.1),
                            ..RemeshParams::default()
                        },
                    )
                });
            let (surface_points, surface, groups) = match remeshed {
                Some(Ok(r)) => {
                    let n = r.triangles.len() as u32;
                    (r.points, r.triangles, (0..n).collect())
                }
                Some(Err(e)) => {
                    println!("{path}	failed	remesh: {e}	0");
                    continue;
                }
                None => (points.clone(), triangles.clone(), Vec::new()),
            };
            let params = ConformingParams {
                size: (!field.is_uniform()).then(|| Arc::new(field)),
                groups,
                min_boundary_size: h / floor("IMPROVE_BOUNDARY_FLOOR", 16.0),
                ..ConformingParams::with_size(h)
            };
            conforming_mesh(&surface_points, &surface, &params)
                .map_err(|e| e.to_string())
                .map(|m| (m.points, m.tets))
        };
        let t_mesh = started.elapsed().as_secs_f64();
        let (p, t) = match mesh {
            Ok(m) => m,
            Err(e) => {
                println!("{path}\tfailed\t{e}\t{t_mesh:.2}");
                continue;
            }
        };
        let before = describe(&p, &t);
        let started = Instant::now();
        let r = improve(
            &p,
            &t,
            &ImproveParams {
                target_degrees: target,
                max_non_orthogonality_degrees: non_orthogonality,
                max_skewness: skewness,
                collapses: std::env::var("IMPROVE_COLLAPSES").map_or(true, |s| s != "0"),
                splits: std::env::var("IMPROVE_SPLITS").map_or(true, |s| s != "0"),
                insertions: std::env::var("IMPROVE_INSERT").map_or(true, |s| s != "0"),
                ..ImproveParams::default()
            },
        );
        let t_improve = started.elapsed().as_secs_f64();
        eprint!("{}", mesher_tet::prof::report());
        if std::env::var("IMPROVE_WORST").is_ok_and(|s| s == "1") {
            eprintln!("{path}: worst face {}", worst_face(&r.points, &r.tets));
        }
        let after = describe(&r.points, &r.tets);
        let s = r.stats;
        println!(
            "{path}\tok\t{before}\t{after}\t{} flips, {} edge removals, {} moves, {} collapses, {} splits, {} insertions, {} passes, poor {} -> {}, min {:.4} -> {:.4}\t{t_mesh:.2}\t{t_improve:.2}",
            s.flips,
            s.edge_removals,
            s.moves,
            s.collapses,
            s.splits,
            s.insertions,
            s.passes,
            s.poor_before,
            s.poor_after,
            s.worst_before,
            s.worst_after
        );
    }
}
