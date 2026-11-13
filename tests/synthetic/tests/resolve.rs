//! Milestone 16: repair rung 2. Parts whose shells overlap or that
//! intersect themselves are resolved to the boundary of the solid they
//! enclose; its components and Euler characteristic, as the rung's record
//! gives them, match a brute-force count on voxels classified by the
//! input's winding number (the surface's Euler characteristic is twice the
//! solid's).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeSet;

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::repair::{RepairedPart, repair};
use mesher_core::{Action, DefectCode, JobConfig, Outcome, Record, RecordLog};
use mesher_geom::TriangleSoup;
use mesher_geom::winding::winding_number_brute;
use mesher_io::parse_obj;
use mesher_synthetic::{Scene, Solid};

fn bar(lo: [f64; 3], hi: [f64; 3]) -> Solid {
    Solid::cuboid(lo, hi)
}

/// The scenes, each one part, and whether rung 2 must apply.
fn scenes() -> Vec<(&'static str, Scene)> {
    let mut pushed = Solid::icosphere(1);
    pushed.points[0] = pushed.points[0].map(|c| -1.5 * c);
    let cube = Solid::unit_cube;
    vec![
        (
            "shell overlap",
            Scene::default().with(
                Some("p"),
                cube().union(&bar([0.5, 0.25, 0.25], [1.5, 1.25, 1.25])),
            ),
        ),
        (
            "bar through a cube",
            Scene::default().with(
                Some("p"),
                cube().union(&bar([-0.5, 0.3, 0.35], [1.5, 0.7, 0.65])),
            ),
        ),
        (
            "common floor",
            Scene::default().with(
                Some("p"),
                cube().union(&bar([0.4, 0.3, 0.0], [1.6, 0.8, 0.7])),
            ),
        ),
        (
            "touching face to face",
            Scene::default().with(
                Some("p"),
                cube().union(&bar([1.0, 0.25, 0.2], [1.8, 0.75, 0.9])),
            ),
        ),
        (
            // Four bars overlapping at the corners: a solid with a hole.
            "square frame",
            Scene::default().with(
                Some("p"),
                bar([0.0, 0.0, 0.0], [4.0, 1.0, 1.0])
                    .union(&bar([0.0, 3.0, 0.0], [4.0, 4.0, 1.0]))
                    .union(&bar([0.0, 0.0, 0.0], [1.0, 4.0, 1.0]))
                    .union(&bar([3.0, 0.0, 0.0], [4.0, 4.0, 1.0])),
            ),
        ),
        ("pushed-through vertex", Scene::single(pushed)),
    ]
}

/// Repair the scene: its records and its one repaired part.
fn repaired(scene: &Scene) -> (Vec<Record>, RepairedPart) {
    let obj = scene.to_obj();
    let mut log = RecordLog::new();
    let model = parse_obj(obj.as_bytes(), &mut log, &mut |_| {});
    let diagnosis = diagnose(&model, &Tolerances::default(), &mut log).unwrap();
    let config = JobConfig::default();
    let selected = select_parts(&diagnosis, &config, &mut log).unwrap();
    let mut overlay = diagnosis.overlay.clone();
    check_part_pairs(&model, &diagnosis, &selected, &mut log, &mut overlay).unwrap();
    let mut parts = repair(
        &model,
        &diagnosis,
        &selected,
        &Tolerances::default(),
        &mut log,
    )
    .unwrap();
    assert_eq!(parts.len(), 1);
    (log.records().to_vec(), parts.remove(0))
}

/// Components (26-connected) and Euler characteristic of the union of the
/// closed voxels whose centres the scene wraps (winding number at least a
/// half), on an `n`³ grid over its box.
fn voxels(scene: &Scene, n: usize) -> (usize, i64) {
    let mut points = Vec::new();
    let mut tris = Vec::new();
    for s in &scene.solids {
        let o = points.len() as u32;
        points.extend(s.1.points.iter().copied());
        tris.extend(s.1.tris.iter().map(|t| t.map(|v| v + o)));
    }
    let soup = TriangleSoup::new(&points, &tris).unwrap();
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    // A margin off the faces' coordinates, so no centre lies on a face.
    let margin = 0.0173;
    let lo = lo.map(|x| x - margin);
    let hi = hi.map(|x| x + margin * 1.31);
    let cell = [0, 1, 2].map(|k| (hi[k] - lo[k]) / n as f64);
    let id = |i: usize, j: usize, k: usize| (i * n + j) * n + k;
    let mut solid = vec![false; n * n * n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let c = [
                    lo[0] + (i as f64 + 0.5) * cell[0],
                    lo[1] + (j as f64 + 0.5) * cell[1],
                    lo[2] + (k as f64 + 0.5) * cell[2],
                ];
                solid[id(i, j, k)] = winding_number_brute(&soup, c) >= 0.5;
            }
        }
    }
    // Euler characteristic of the cubical complex: lattice vertices,
    // edges and faces touching a solid voxel, minus edges, minus voxels.
    let occupied = |i: i64, j: i64, k: i64| -> bool {
        let n = n as i64;
        (0..n).contains(&i)
            && (0..n).contains(&j)
            && (0..n).contains(&k)
            && solid[id(i as usize, j as usize, k as usize)]
    };
    let (mut v, mut e, mut f, mut c) = (0i64, 0i64, 0i64, 0i64);
    let m = n as i64;
    for i in 0..=m {
        for j in 0..=m {
            for k in 0..=m {
                // The cells around lattice point (i, j, k) by offsets.
                let any = |di: &[i64], dj: &[i64], dk: &[i64]| {
                    di.iter().any(|&a| {
                        dj.iter()
                            .any(|&b| dk.iter().any(|&d| occupied(i + a, j + b, k + d)))
                    })
                };
                if any(&[-1, 0], &[-1, 0], &[-1, 0]) {
                    v += 1;
                }
                // Edges from (i, j, k) along each axis.
                if i < m && any(&[0], &[-1, 0], &[-1, 0]) {
                    e += 1;
                }
                if j < m && any(&[-1, 0], &[0], &[-1, 0]) {
                    e += 1;
                }
                if k < m && any(&[-1, 0], &[-1, 0], &[0]) {
                    e += 1;
                }
                // Faces at (i, j, k) normal to each axis.
                if j < m && k < m && any(&[-1, 0], &[0], &[0]) {
                    f += 1;
                }
                if i < m && k < m && any(&[0], &[-1, 0], &[0]) {
                    f += 1;
                }
                if i < m && j < m && any(&[0], &[0], &[-1, 0]) {
                    f += 1;
                }
                if occupied(i, j, k) {
                    c += 1;
                }
            }
        }
    }
    // Components, 26-connected.
    let mut seen = vec![false; n * n * n];
    let mut components = 0;
    for start in 0..n * n * n {
        if !solid[start] || seen[start] {
            continue;
        }
        components += 1;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(x) = stack.pop() {
            let (i, j, k) = ((x / (n * n)) as i64, ((x / n) % n) as i64, (x % n) as i64);
            for a in -1..=1 {
                for b in -1..=1 {
                    for d in -1..=1 {
                        if occupied(i + a, j + b, k + d) {
                            let y = id((i + a) as usize, (j + b) as usize, (k + d) as usize);
                            if !seen[y] {
                                seen[y] = true;
                                stack.push(y);
                            }
                        }
                    }
                }
            }
        }
    }
    (components, v - e + f - c)
}

#[test]
fn resolved_surfaces_match_a_brute_force_count() {
    for (name, scene) in scenes() {
        let (records, part) = repaired(&scene);
        let summary = part.summary;
        assert_eq!(summary.rung, 2, "{name}");
        assert!(summary.gate_passed, "{name}");
        let (components, euler) = voxels(&scene, 48);
        assert_eq!(summary.components_after, components as u64, "{name}");
        assert_eq!(summary.euler_after, 2 * euler, "{name}");
        // Every finding of the intersection codes is dealt with.
        let codes = [DefectCode::IntSelfIntersection, DefectCode::IntShellOverlap];
        let found: BTreeSet<DefectCode> = records
            .iter()
            .map(|r| r.code)
            .filter(|c| codes.contains(c))
            .collect();
        assert!(!found.is_empty(), "{name}");
        for rec in records.iter().filter(|r| codes.contains(&r.code)) {
            assert_eq!(rec.action, Some(Action::Resolve), "{name}: {rec:?}");
            assert_eq!(rec.outcome, Some(Outcome::Applied), "{name}: {rec:?}");
        }
    }
}
