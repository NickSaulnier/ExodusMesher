//! Exit criteria for Milestones 3 and 4, at several scales and offsets:
//! - diagnosis: every injected defect is detected with its exact location, and
//!   clean shapes produce no findings (prints precision and recall per code);
//! - repair: every repairable finding has a matching repair record, the
//!   repaired parts pass the surface gate (or fail it exactly where expected),
//!   no guard trips, and the face and vertex counts before and after repair are
//!   fully explained by the repair records;
//! - stuffing (Milestone 6): every repaired part of a closed or nearly closed
//!   case is meshed with positive tets inside the proven dihedral bounds, and
//!   accepted meshes are valid and match the enclosed volume;
//! - the exact constrained kernel (Milestone 12): every repaired part of a
//!   closed case that passes the surface gate has every segment and triangle
//!   recovered in exact form, and its inside region has the enclosed volume;
//!   rounded to f64 (Milestone 13), every tet is positive, exactly and in
//!   floating point, and the volume still the enclosed one.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::vec_init_then_push
)]

use std::collections::{BTreeMap, BTreeSet};

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::repair::{RepairedPart, repair};
use mesher_core::report::MAX_LOCATION_ITEMS;
use mesher_core::volume::{MeshOutcome, mesh_parts};
use mesher_core::{
    Action, DefectCode, JobConfig, Kernel, Outcome, Record, RecordLog, Severity, Stage,
};
use mesher_io::parse_obj;
use mesher_synthetic::{Scene, Solid};
use mesher_tet::stuffing::{DIHEDRAL_MAX_DEGREES, DIHEDRAL_MIN_DEGREES};

#[derive(Debug, Clone, PartialEq)]
enum Key {
    Faces(Vec<u32>),
    Vertices(Vec<u32>),
    /// At least one record of the code; locations checked elsewhere or not fixed.
    Any,
}

#[derive(Debug, Clone)]
struct Case {
    name: &'static str,
    scene: Scene,
    parts: Vec<&'static str>,
    expect: Vec<(DefectCode, Key)>,
    parts_count: Option<usize>,
    /// Also run scaled and translated copies.
    variants: bool,
    /// Codes the surface gate should fail on after repair (empty: it passes).
    gate_fails: Vec<DefectCode>,
    /// Codes whose repair should be recorded as failed.
    repair_fails: Vec<DefectCode>,
    /// The surface is open on purpose (stuffing is not expected to apply).
    open: bool,
    /// The conformal kernel should refuse the repaired surface (shells that
    /// touch at duplicated vertices), so stuffing must take over.
    conformal_fails: bool,
}

fn case(name: &'static str, scene: Scene, expect: Vec<(DefectCode, Key)>) -> Case {
    Case {
        name,
        scene,
        parts: Vec::new(),
        expect,
        parts_count: None,
        variants: true,
        gate_fails: Vec::new(),
        repair_fails: Vec::new(),
        open: false,
        conformal_fails: false,
    }
}

fn faces(f: impl IntoIterator<Item = u32>) -> Key {
    let mut v: Vec<u32> = f.into_iter().collect();
    v.sort_unstable();
    Key::Faces(v)
}

fn verts(f: impl IntoIterator<Item = u32>) -> Key {
    let mut v: Vec<u32> = f.into_iter().collect();
    v.sort_unstable();
    Key::Vertices(v)
}

/// Records list at most `MAX_LOCATION_ITEMS` entries (the smallest ones).
fn capped(key: &Key) -> Key {
    let cut = |v: &Vec<u32>| v.iter().copied().take(MAX_LOCATION_ITEMS).collect();
    match key {
        Key::Faces(v) => Key::Faces(cut(v)),
        Key::Vertices(v) => Key::Vertices(cut(v)),
        Key::Any => Key::Any,
    }
}

fn quad() -> Solid {
    Solid {
        points: vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        tris: vec![[0, 1, 2], [0, 2, 3]],
    }
}

fn cases() -> Vec<Case> {
    use DefectCode::*;
    let cube = Solid::unit_cube;
    let mut out = Vec::new();

    // Clean shapes.
    out.push(Case {
        parts_count: Some(1),
        ..case("clean cube", Scene::single(cube()), vec![])
    });
    out.push(case(
        "clean grid cube",
        Scene::single(Solid::grid_cube(3)),
        vec![],
    ));
    out.push(case(
        "clean icosphere",
        Scene::single(Solid::icosphere(2)),
        vec![],
    ));
    out.push(case(
        "clean torus",
        Scene::single(Solid::torus(2.0, 0.75, 24, 12)),
        vec![],
    ));
    out.push(Case {
        parts_count: Some(1),
        ..case(
            "clean hollow cube (cavity joins its component)",
            Scene::single(Solid::grid_cube(3).union(&Solid::cuboid([1.0; 3], [2.0; 3]).flipped())),
            vec![],
        )
    });
    out.push(Case {
        parts_count: Some(2),
        ..case(
            "clean separate cubes",
            Scene::single(cube().union(&Solid::cuboid([3.0; 3], [4.0; 3]))),
            vec![],
        )
    });
    out.push(Case {
        parts_count: Some(2),
        ..case(
            "touching boxes as two parts",
            Scene::default().with(Some("left"), cube()).with(
                Some("right"),
                Solid::cuboid([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]),
            ),
            vec![(PartSharedInterface, Key::Any)],
        )
    });

    // GEO_
    for (name, offset) in [
        ("duplicate vertex, exact", 0.0),
        ("duplicate vertex, within tolerance", 1e-9),
    ] {
        let mut s = cube();
        let d = s.diagonal();
        let p = s.points[0];
        s.points.push([p[0] + offset * d, p[1], p[2]]);
        let t = s.tris.iter().position(|t| t.contains(&0)).unwrap();
        for v in &mut s.tris[t] {
            if *v == 0 {
                *v = 8;
            }
        }
        out.push(case(
            name,
            Scene::single(s),
            vec![(GeoDuplicateVertex, verts([0, 8]))],
        ));
    }
    {
        let mut s = cube();
        s.tris.push([0, 0, 1]);
        out.push(case(
            "degenerate face, repeated vertex",
            Scene::single(s),
            vec![(GeoDegenerateFace, faces([12]))],
        ));
    }
    {
        let mut s = Solid::grid_cube(2);
        let line = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]].map(|p| s.find(p).unwrap());
        s.tris.push(line);
        let n = s.tris.len() as u32 - 1;
        out.push(case(
            "degenerate face, collinear",
            Scene::single(s),
            vec![(GeoDegenerateFace, faces([n]))],
        ));
    }
    {
        // A collapsed sliver at a T-junction: the front face uses the long edge
        // p-q while the bottom face uses p-m and m-q.
        let mut s = Solid::grid_cube(2);
        let pts = s.points.clone();
        let id = |p: [f64; 3]| s.find(p).unwrap();
        let (p, m, q) = (
            id([0.0, 0.0, 0.0]),
            id([1.0, 0.0, 0.0]),
            id([2.0, 0.0, 0.0]),
        );
        let (r, a, b) = (
            id([1.0, 0.0, 1.0]),
            id([0.0, 0.0, 1.0]),
            id([2.0, 0.0, 1.0]),
        );
        s.tris.retain(|t| {
            !t.iter()
                .all(|&v| pts[v as usize][1] == 0.0 && pts[v as usize][2] <= 1.0)
        });
        s.tris.extend([[p, q, r], [p, r, a], [q, b, r], [p, m, q]]);
        let n = s.tris.len() as u32 - 1;
        out.push(case(
            "degenerate face at a T-junction",
            Scene::single(s),
            vec![(GeoDegenerateFace, faces([n]))],
        ));
    }
    for (name, reversed) in [
        ("duplicate face", false),
        ("duplicate face, reversed", true),
    ] {
        let mut s = cube();
        let [a, b, c] = s.tris[0];
        s.tris.push(if reversed { [a, c, b] } else { [a, b, c] });
        out.push(case(
            name,
            Scene::single(s),
            vec![(GeoDuplicateFace, faces([0, 12]))],
        ));
    }

    // TOP_
    {
        let mut s = cube();
        let removed = s.tris.remove(0);
        out.push(case(
            "hole, one triangle",
            Scene::single(s),
            vec![(TopBoundaryEdge, verts(removed))],
        ));
    }
    {
        let mut s = Solid::grid_cube(2);
        let b = s.tris.remove(1);
        let a = s.tris.remove(0);
        let rim: BTreeSet<u32> = a.iter().chain(&b).copied().collect();
        out.push(case(
            "hole, one quad",
            Scene::single(s),
            vec![(TopBoundaryEdge, verts(rim))],
        ));
    }
    {
        let ico = Solid::icosphere(2);
        let bowl = Solid {
            tris: ico
                .tris
                .iter()
                .copied()
                .filter(|t| t.iter().all(|&v| ico.points[v as usize][2] >= 0.0))
                .collect(),
            points: ico.points.clone(),
        };
        // Too large to fill, and stuffing's boundary is too far from it, so
        // repair rung 3, the last resort, wraps it (Milestone 17), with the
        // winding guard closing it across the rim.
        out.push(Case {
            gate_fails: vec![TopBoundaryEdge],
            repair_fails: vec![TopBoundaryEdge],
            open: true,
            ..case(
                "hole, open bowl (too large to fill)",
                Scene::single(bowl.clone()),
                vec![(TopBoundaryEdge, Key::Any)],
            )
        });
        // Open and intersecting: rung 2 refuses an open part; stuffing or,
        // failing that, rung 3 meshes it.
        out.push(Case {
            gate_fails: vec![TopBoundaryEdge, IntShellOverlap],
            repair_fails: vec![TopBoundaryEdge, IntShellOverlap],
            open: true,
            ..case(
                "open bowl with a box through its wall",
                Scene::default().with(
                    Some("bowl"),
                    bowl.union(&Solid::cuboid([0.7, -0.2, 0.1], [1.3, 0.2, 0.5])),
                ),
                vec![(TopBoundaryEdge, Key::Any), (IntShellOverlap, Key::Any)],
            )
        });
    }
    {
        let s = cube().union(&Solid::cuboid([1.0, 1.0, 0.0], [2.0, 2.0, 1.0]));
        let e = [[1.0, 1.0, 0.0], [1.0, 1.0, 1.0]].map(|p| s.find(p).unwrap());
        out.push(Case {
            conformal_fails: true,
            ..case(
                "non-manifold edge (cubes sharing an edge)",
                Scene::single(s),
                vec![(TopNonmanifoldEdge, verts(e))],
            )
        });
    }
    {
        let s = cube().union(&Solid::cuboid([1.0; 3], [2.0; 3]));
        let v = s.find([1.0; 3]).unwrap();
        out.push(Case {
            parts_count: Some(1),
            conformal_fails: true,
            ..case(
                "non-manifold vertex (cubes sharing a corner)",
                Scene::single(s),
                vec![(TopNonmanifoldVertex, verts([v]))],
            )
        });
    }
    {
        let mut s = cube();
        s.tris[0].swap(1, 2);
        out.push(case(
            "inconsistent orientation, one face",
            Scene::single(s),
            vec![(TopInconsistentOrientation, faces([0]))],
        ));
    }
    {
        let mut s = Solid::icosphere(1);
        for t in &mut s.tris[..5] {
            t.swap(1, 2);
        }
        out.push(case(
            "inconsistent orientation, five faces",
            Scene::single(s),
            vec![(TopInconsistentOrientation, faces(0..5))],
        ));
    }
    out.push(case(
        "inverted cube",
        Scene::single(cube().flipped()),
        vec![(TopInvertedShell, faces(0..12))],
    ));
    {
        let torus = Solid::torus(2.0, 0.75, 24, 12).flipped();
        let n = torus.tris.len() as u32;
        out.push(case(
            "inverted torus",
            Scene::single(torus),
            vec![(TopInvertedShell, faces(0..n))],
        ));
    }
    {
        let outer = Solid::grid_cube(3);
        let n = outer.tris.len() as u32;
        let s = outer.union(&Solid::cuboid([1.0; 3], [2.0; 3]));
        out.push(case(
            "inverted cavity (inner shell facing outward)",
            Scene::default().with(Some("hollow"), s),
            vec![(TopInvertedShell, faces(n..n + 12))],
        ));
    }
    out.push(case(
        "open sheet",
        Scene::single(quad()),
        vec![(TopOpenSheet, faces([0, 1]))],
    ));
    {
        let mut sheet = quad();
        for p in &mut sheet.points {
            p[2] = 5.0;
        }
        out.push(case(
            "open sheet beside a solid",
            Scene::default().with(Some("mixed"), cube().union(&sheet)),
            vec![(TopOpenSheet, faces([12, 13]))],
        ));
    }

    // INT_
    {
        let mut s = Solid::icosphere(1);
        s.points[0] = s.points[0].map(|c| -1.5 * c);
        // Rung 2 resolves it: the gate passes on the resolved surface.
        out.push(case(
            "self-intersection (vertex pushed through)",
            Scene::single(s),
            vec![(IntSelfIntersection, Key::Any)],
        ));
    }
    let shifted = || Solid::cuboid([0.5, 0.25, 0.25], [1.5, 1.25, 1.25]);
    out.push(case(
        "shell overlap within one part",
        Scene::default().with(Some("pair"), cube().union(&shifted())),
        vec![(IntShellOverlap, Key::Any)],
    ));
    {
        // A cube and a hooked prism sharing an edge, the hook reaching back
        // into the cube: one shell that intersects itself, until rung 1
        // splits the edge and the intersections are between two shells.
        // Rung 2 resolves them all (Thingi10k 44104 kept its diagnosis
        // errors). The union still has the edge, split again, so stuffing
        // meshes it.
        let hook = Solid::prism(
            &[
                [1.5, 1.0],
                [0.5, 0.25],
                [2.0, 0.25],
                [2.0, 1.5],
                [1.0, 1.5],
                [1.0, 1.0],
            ],
            0.0,
            1.0,
            0,
        );
        let s = cube().union(&hook);
        let e = [[1.0, 1.0, 0.0], [1.0, 1.0, 1.0]].map(|p| s.find(p).unwrap());
        out.push(Case {
            conformal_fails: true,
            ..case(
                "self-intersection across a non-manifold edge",
                Scene::single(s),
                vec![
                    (TopNonmanifoldEdge, verts(e)),
                    (IntSelfIntersection, Key::Any),
                ],
            )
        });
    }
    // Rung 2 on more ways shells meet: a bar through a cube, boxes
    // overlapping on a common floor (coplanar faces), and boxes touching
    // face to face (their interface goes).
    out.push(case(
        "bar through a cube",
        Scene::default().with(
            Some("cross"),
            cube().union(&Solid::cuboid([-0.5, 0.3, 0.35], [1.5, 0.7, 0.65])),
        ),
        vec![(IntShellOverlap, Key::Any)],
    ));
    out.push(case(
        "boxes overlapping on a common floor",
        Scene::default().with(
            Some("floor"),
            cube().union(&Solid::cuboid([0.4, 0.3, 0.0], [1.6, 0.8, 0.7])),
        ),
        vec![(IntShellOverlap, Key::Any)],
    ));
    out.push(case(
        "boxes touching face to face",
        Scene::default().with(
            Some("touch"),
            cube().union(&Solid::cuboid([1.0, 0.25, 0.2], [1.8, 0.75, 0.9])),
        ),
        vec![(IntShellOverlap, Key::Any)],
    ));

    // PART_
    out.push(case(
        "overlapping parts",
        Scene::default()
            .with(Some("a"), cube())
            .with(Some("b"), shifted()),
        vec![(PartOverlap, Key::Any)],
    ));
    out.push(Case {
        parts_count: Some(2),
        ..case(
            "part inside another part",
            Scene::single(
                Solid::cuboid([0.0; 3], [3.0; 3]).union(&Solid::cuboid([1.0; 3], [2.0; 3])),
            ),
            vec![(PartOverlap, Key::Any)],
        )
    });
    out.push(Case {
        parts: vec!["nope"],
        ..case(
            "selecting a missing part",
            Scene::single(cube()),
            vec![(PartEmptySelection, Key::Any)],
        )
    });

    // QUAL_
    {
        let mut s = cube();
        let [a, b, c] = s.tris[0];
        let (pa, pb, pc) = (
            s.points[a as usize],
            s.points[b as usize],
            s.points[c as usize],
        );
        let mid = [0, 1, 2].map(|k| (pa[k] + pb[k]) / 2.0);
        s.points
            .push([0, 1, 2].map(|k| mid[k] + 0.002 * (pc[k] - mid[k])));
        let m = s.points.len() as u32 - 1;
        s.tris[0] = [a, b, m];
        s.tris.push([b, c, m]);
        s.tris.push([c, a, m]);
        out.push(case(
            "sliver face",
            Scene::single(s),
            vec![(QualSliverFace, faces([0]))],
        ));
    }
    for (name, scale) in [("tiny model", 1e-6), ("huge model", 1e7)] {
        out.push(Case {
            variants: false,
            ..case(
                name,
                Scene::single(cube()).transformed(scale, 0.0),
                vec![(QualSuspectScale, Key::Any)],
            )
        });
    }
    out
}

struct Run {
    records: Vec<Record>,
    parts: usize,
    repaired: Vec<RepairedPart>,
    meshes: MeshOutcome,
}

/// Coarser than the default, to keep the suite fast.
const SUITE_ELEMENT_SIZE: f64 = 0.06;

fn run(scene: &Scene, parts: &[&str]) -> Run {
    let obj = scene.to_obj();
    let mut log = RecordLog::new();
    let model = parse_obj(obj.as_bytes(), &mut log, &mut |_| {});
    assert!(
        log.records().is_empty(),
        "parse records: {:?}",
        log.records()
    );
    let diagnosis = diagnose(&model, &Tolerances::default(), &mut log).unwrap();
    let config = JobConfig {
        parts: parts.iter().map(|s| (*s).to_owned()).collect(),
        ..JobConfig::default()
    };
    let mut repaired = Vec::new();
    let mut meshes = MeshOutcome::default();
    if let Some(selected) = select_parts(&diagnosis, &config, &mut log) {
        let mut overlay = diagnosis.overlay.clone();
        check_part_pairs(&model, &diagnosis, &selected, &mut log, &mut overlay).unwrap();
        repaired = repair(
            &model,
            &diagnosis,
            &selected,
            &Tolerances::default(),
            &mut log,
        )
        .unwrap();
        let mesh_config = JobConfig {
            element_size: SUITE_ELEMENT_SIZE,
            ..config.clone()
        };
        meshes = mesh_parts(&model, &repaired, &mesh_config, &mut log).unwrap();
    }
    Run {
        records: log.records().to_vec(),
        parts: diagnosis.parts.len(),
        repaired,
        meshes,
    }
}

/// Milestone 6: stuffing on every repaired part of a closed or nearly closed case.
/// Parts meshed, accepted, and meshed by the conformal Delaunay kernel.
#[derive(Default)]
struct MeshTally {
    meshed: u32,
    accepted: u32,
    delaunay: u32,
}

/// Check every part's volume mesh: valid, inside the kernel's quality bounds
/// (proven for stuffing; 10 degrees for Delaunay unless the input itself has
/// a sliver face) and of the right volume (exact for a conformal mesh).
fn check_meshes(c: &Case, run: &Run, label: &str, failures: &mut Vec<String>) -> MeshTally {
    let mut tally = MeshTally::default();
    if c.open {
        return tally;
    }
    let sliver_input = c
        .expect
        .iter()
        .any(|(code, _)| *code == DefectCode::QualSliverFace);
    for r in &run.repaired {
        let Some((_, s)) = run.meshes.summaries.iter().find(|(i, _)| *i == r.part) else {
            failures.push(format!("{label}: part {} was not meshed", r.name));
            continue;
        };
        if s.tets == 0 {
            failures.push(format!("{label}: part {} produced no tets", r.name));
        }
        let (lo, hi) = match s.kernel {
            Kernel::Stuffing => (DIHEDRAL_MIN_DEGREES, DIHEDRAL_MAX_DEGREES),
            Kernel::Delaunay | Kernel::Exact if sliver_input => (0.0, 180.0),
            Kernel::Delaunay | Kernel::Exact => (10.0, 170.0),
        };
        if s.kernel == Kernel::Delaunay {
            tally.delaunay += 1;
        }
        if c.conformal_fails && s.kernel != Kernel::Stuffing {
            failures.push(format!("{label}: part {} was meshed conformally", r.name));
        }
        if s.min_dihedral_degrees < lo || s.max_dihedral_degrees > hi {
            failures.push(format!(
                "{label}: part {} dihedral range {:.2}..{:.2}",
                r.name, s.min_dihedral_degrees, s.max_dihedral_degrees
            ));
        }
        if !s.accepted {
            eprintln!(
                "  not accepted: {label} part {} (hausdorff {:.3} of the diagonal)",
                r.name, s.hausdorff_relative
            );
        }
        if s.accepted {
            tally.accepted += 1;
        }
        // The enclosed volume is the region's volume only for a surface that
        // passed the gate (overlapping shells count their overlap twice).
        if s.accepted && c.gate_fails.is_empty() {
            let error =
                (s.volume - s.surface_volume).abs() / s.surface_volume.abs().max(f64::MIN_POSITIVE);
            let allowed = match s.kernel {
                Kernel::Stuffing => 0.15,
                Kernel::Delaunay | Kernel::Exact => 1e-9,
            };
            if error > allowed {
                failures.push(format!(
                    "{label}: part {} ({:?}) volume off by {:.3e}",
                    r.name, s.kernel, error
                ));
            }
        }
    }
    for m in &run.meshes.meshed {
        let mesh = mesher_core::volume::assemble(std::slice::from_ref(m));
        let problems = mesh.validate(5);
        if !problems.is_empty() {
            failures.push(format!(
                "{label}: part {} mesh invalid: {problems:?}",
                m.name
            ));
        }
    }
    for code in [
        DefectCode::VolInvertedTet,
        DefectCode::VolLowQuality,
        DefectCode::MeshBudgetExceeded,
        DefectCode::MeshRecoveryFailed,
    ] {
        let expected = match code {
            DefectCode::VolLowQuality => sliver_input,
            DefectCode::MeshRecoveryFailed => c.conformal_fails,
            _ => false,
        };
        let found = run.records.iter().any(|r| r.code == code);
        if code == DefectCode::MeshRecoveryFailed && c.conformal_fails && !found {
            failures.push(format!("{label}: expected {}", code.as_str()));
        }
        if expected {
            continue;
        }
        if found {
            failures.push(format!("{label}: {}", code.as_str()));
        }
    }
    tally.meshed = run.repaired.len() as u32;
    tally
}

/// The exact constrained kernel's totals over the suite.
#[derive(Default)]
struct ExactTally {
    parts: u32,
    steiner_segments: usize,
    steiner_faces: usize,
    wrapped: usize,
    seconds: f64,
    /// Parts with a tet that rounding to f64 turned over or flattened.
    rounding_invalid: u32,
    untangles: usize,
}

/// Milestone 12: the exact kernel (`mesher_tet::cdt::mesh_region`) on every
/// repaired part of a closed case that passes the surface gate (and whose
/// shells do not touch at duplicated vertices, which no conformal kernel
/// takes): every input triangle recovered, and the inside region exactly
/// the enclosed volume; then (Milestone 13) rounded to f64 and untangled.
fn check_exact(c: &Case, run: &Run, label: &str, failures: &mut Vec<String>, t: &mut ExactTally) {
    if c.open || !c.gate_fails.is_empty() || c.conformal_fails {
        return;
    }
    for r in &run.repaired {
        let started = std::time::Instant::now();
        let region = match mesher_tet::cdt::mesh_region(&r.points, &r.faces, 0, 5) {
            Ok(region) => region,
            Err(e) => {
                failures.push(format!("{label}: part {} exact kernel: {e}", r.name));
                continue;
            }
        };
        t.seconds += started.elapsed().as_secs_f64();
        let det = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
            a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
        };
        // Relative to a corner, so large offsets do not swamp the volume.
        let o = r.points.first().copied().unwrap_or([0.0; 3]);
        let rel = |p: [f64; 3]| [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
        let enclosed: f64 = r
            .faces
            .iter()
            .map(|f| {
                det(
                    rel(r.points[f[0] as usize]),
                    rel(r.points[f[1] as usize]),
                    rel(r.points[f[2] as usize]),
                ) / 6.0
            })
            .sum::<f64>()
            .abs();
        let at = |v: u32| rel(region.dt().points()[v as usize].approx());
        let volume: f64 = region
            .tets
            .iter()
            .map(|v| mesher_tet::quality::float_volume(v.map(at)))
            .sum();
        let error = (volume - enclosed).abs() / enclosed.max(f64::MIN_POSITIVE);
        if error > 1e-9 {
            failures.push(format!(
                "{label}: part {} exact kernel volume off by {error:.3e}",
                r.name
            ));
        }
        t.parts += 1;
        t.steiner_segments += region.recovered.steiner;
        t.steiner_faces += region.recovered.steiner_on_faces;
        t.wrapped += region.recovered.wrapped;
        // Milestone 13: rounded to f64, every tet positive exactly and in
        // floating point, the volume still the enclosed one.
        match mesher_tet::cdt::round_region(region, &r.points, &r.faces) {
            Ok(m) => {
                if m.inverted_before + m.flat_before > 0 {
                    t.rounding_invalid += 1;
                }
                t.untangles += m.untangles.len();
                let volume: f64 = m
                    .tets
                    .iter()
                    .map(|v| {
                        mesher_tet::quality::float_volume(v.map(|x| rel(m.points[x as usize])))
                    })
                    .sum();
                let error = (volume - enclosed).abs() / enclosed.max(f64::MIN_POSITIVE);
                if error > 1e-9 {
                    failures.push(format!(
                        "{label}: part {} rounded volume off by {error:.3e}",
                        r.name
                    ));
                }
            }
            Err(e) => failures.push(format!("{label}: part {} rounding: {e}", r.name)),
        }
    }
}

/// How the suite compares locations for each repairable code.
fn repair_key(code: DefectCode, r: &Record) -> Option<Key> {
    use DefectCode::*;
    match code {
        GeoDegenerateFace | GeoDuplicateFace | TopInconsistentOrientation | TopInvertedShell => {
            Some(Key::Faces(r.location.faces.clone()))
        }
        GeoDuplicateVertex | TopBoundaryEdge | TopNonmanifoldEdge | TopNonmanifoldVertex => {
            Some(Key::Vertices(r.location.vertices.clone()))
        }
        _ => None,
    }
}

fn check_repair(c: &Case, run: &Run, label: &str, failures: &mut Vec<String>) {
    let repaired: BTreeSet<&str> = run.repaired.iter().map(|r| r.name.as_str()).collect();
    let stage = |s: Stage| run.records.iter().filter(move |r| r.stage == s);

    // Every repairable finding on a repaired part has its repair record.
    for finding in stage(Stage::Diagnose) {
        let Some(key) = repair_key(finding.code, finding) else {
            continue;
        };
        if !finding
            .part
            .as_deref()
            .is_some_and(|p| repaired.contains(p))
        {
            continue;
        }
        let want = if c.repair_fails.contains(&finding.code) {
            Outcome::Failed
        } else {
            Outcome::Applied
        };
        let found = stage(Stage::Repair).any(|r| {
            r.code == finding.code
                && r.part == finding.part
                && r.outcome == Some(want)
                && repair_key(r.code, r) == Some(key.clone())
        });
        if !found {
            failures.push(format!(
                "{label}: no {want:?} repair record for {} at {key:?}",
                finding.code.as_str()
            ));
        }
    }

    // A part whose intersections rung 2 resolved keeps no intersection
    // error from diagnosis, whichever kind each became after rung 1.
    for r in stage(Stage::Repair) {
        if r.action != Some(Action::Resolve) || r.outcome != Some(Outcome::Applied) {
            continue;
        }
        for finding in stage(Stage::Diagnose).filter(|f| {
            f.part == r.part
                && matches!(
                    f.code,
                    DefectCode::IntSelfIntersection | DefectCode::IntShellOverlap
                )
        }) {
            if finding.severity == Severity::Error
                || finding.action != Some(Action::Resolve)
                || finding.outcome != Some(Outcome::Applied)
            {
                failures.push(format!(
                    "{label}: {} on a resolved part left as {:?} {:?} {:?}",
                    finding.code.as_str(),
                    finding.severity,
                    finding.action,
                    finding.outcome
                ));
            }
        }
    }

    // Gate verdict and guards.
    let gate: BTreeSet<DefectCode> = stage(Stage::SurfaceGate).map(|r| r.code).collect();
    let want: BTreeSet<DefectCode> = c.gate_fails.iter().copied().collect();
    if gate != want {
        failures.push(format!(
            "{label}: gate failed on {gate:?}, expected {want:?}"
        ));
    }
    for p in &run.repaired {
        if p.summary.gate_passed != want.is_empty() {
            failures.push(format!(
                "{label}: part {} gate_passed = {}",
                p.name, p.summary.gate_passed
            ));
        }
    }
    if run
        .records
        .iter()
        .any(|r| r.code == DefectCode::RepairGuardExceeded)
    {
        failures.push(format!("{label}: a repair guard tripped"));
    }

    // Accounting: counts before and after are explained by the repair
    // records (for a part wrapped as the last resort, its wrap's counts).
    let finals: Vec<&RepairedPart> = run
        .repaired
        .iter()
        .map(|r| {
            run.meshes
                .wrapped
                .iter()
                .find(|w| w.part == r.part)
                .unwrap_or(r)
        })
        .collect();
    for p in finals {
        let applied: Vec<&Record> = stage(Stage::Repair)
            .filter(|r| {
                r.part.as_deref() == Some(p.name.as_str()) && r.outcome == Some(Outcome::Applied)
            })
            .collect();
        let sum = |k: &str| -> i64 {
            applied
                .iter()
                .map(|r| {
                    r.measured
                        .get(k)
                        .and_then(serde_json::Value::as_i64)
                        .unwrap_or(0)
                })
                .sum()
        };
        let s = p.summary;
        let faces = s.faces_after as i64 - s.faces_before as i64;
        if faces != sum("faces_added") - sum("faces_removed") {
            failures.push(format!(
                "{label}: part {} face count change {faces} not explained by records",
                p.name
            ));
        }
        let verts = s.vertices_after as i64 - s.vertices_before as i64;
        if verts != sum("vertices_added") {
            failures.push(format!(
                "{label}: part {} vertex count change {verts} not explained by records",
                p.name
            ));
        }
    }
}

#[derive(Default, Debug, Clone, Copy)]
struct Tally {
    tp: u32,
    fp: u32,
    fn_: u32,
}

/// Golden reports (Milestone 10): for each case and variant, every record by
/// code, severity, action and outcome with its count, and each part's mesh
/// (kernel, acceptance, size and quality). Codes, counts, kernels and
/// acceptance must match exactly; mesh metrics within `GOLDEN_BANDS`.
/// `UPDATE_GOLDEN=1` rewrites the file.
const GOLDEN_FILE: &str = "golden.txt";

/// Tolerances for mesh metrics: relative tet count, and degrees for the
/// minimum dihedral and mean non-orthogonality.
const GOLDEN_BANDS: (f64, f64, f64) = (0.10, 2.0, 2.0);

#[derive(Debug, Clone, PartialEq)]
struct GoldenPart {
    part: usize,
    kernel: String,
    accepted: bool,
    tets: u64,
    min_dihedral: f64,
    mean_non_orthogonality: f64,
}

#[derive(Debug, Clone, Default, PartialEq)]
struct Golden {
    records: BTreeMap<String, u64>,
    parts: Vec<GoldenPart>,
}

fn serde_name<T: serde::Serialize>(v: &T) -> String {
    serde_json::to_value(v)
        .ok()
        .and_then(|v| v.as_str().map(str::to_owned))
        .unwrap_or_else(|| "?".to_owned())
}

fn golden_of(run: &Run) -> Golden {
    let mut g = Golden::default();
    for r in &run.records {
        let action = match (&r.action, &r.outcome) {
            (Some(a), Some(o)) => format!(" {}/{}", serde_name(a), serde_name(o)),
            _ => String::new(),
        };
        let key = format!("{} {}{action}", r.code.as_str(), serde_name(&r.severity));
        *g.records.entry(key).or_default() += 1;
    }
    for (part, m) in &run.meshes.summaries {
        g.parts.push(GoldenPart {
            part: *part,
            kernel: serde_name(&m.kernel),
            accepted: m.accepted,
            tets: m.tets,
            min_dihedral: (m.min_dihedral_degrees * 10.0).round() / 10.0,
            mean_non_orthogonality: (m.mean_non_orthogonality_degrees * 10.0).round() / 10.0,
        });
    }
    g.parts.sort_by_key(|p| p.part);
    g
}

/// One line per case: `label <TAB> records <TAB> parts`.
fn golden_line(label: &str, g: &Golden) -> String {
    let records: Vec<String> = g.records.iter().map(|(k, n)| format!("{k} x{n}")).collect();
    let parts: Vec<String> = g
        .parts
        .iter()
        .map(|p| {
            format!(
                "{} {} {} {} {:.1} {:.1}",
                p.part,
                p.kernel,
                if p.accepted { "accepted" } else { "rejected" },
                p.tets,
                p.min_dihedral,
                p.mean_non_orthogonality
            )
        })
        .collect();
    format!("{label}\t{}\t{}", records.join("; "), parts.join("; "))
}

fn parse_golden(text: &str) -> BTreeMap<String, Golden> {
    let mut out = BTreeMap::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
    {
        let f: Vec<&str> = line.split('\t').collect();
        let mut g = Golden::default();
        for r in f
            .get(1)
            .unwrap_or(&"")
            .split("; ")
            .filter(|s| !s.is_empty())
        {
            let (k, n) = r.rsplit_once(" x").unwrap();
            g.records.insert(k.to_owned(), n.parse().unwrap());
        }
        for p in f
            .get(2)
            .unwrap_or(&"")
            .split("; ")
            .filter(|s| !s.is_empty())
        {
            let w: Vec<&str> = p.split(' ').collect();
            g.parts.push(GoldenPart {
                part: w[0].parse().unwrap(),
                kernel: w[1].to_owned(),
                accepted: w[2] == "accepted",
                tets: w[3].parse().unwrap(),
                min_dihedral: w[4].parse().unwrap(),
                mean_non_orthogonality: w[5].parse().unwrap(),
            });
        }
        out.insert(f[0].to_owned(), g);
    }
    out
}

/// Differences between a run and its golden entry, beyond the bands.
fn golden_diff(label: &str, want: &Golden, got: &Golden) -> Vec<String> {
    let mut out = Vec::new();
    if want.records != got.records {
        out.push(format!(
            "{label}: records changed\n    golden: {:?}\n    now:    {:?}",
            want.records, got.records
        ));
    }
    if want.parts.len() != got.parts.len() {
        out.push(format!(
            "{label}: {} meshed parts, golden {}",
            got.parts.len(),
            want.parts.len()
        ));
        return out;
    }
    let (tets, dihedral, non_orth) = GOLDEN_BANDS;
    for (w, g) in want.parts.iter().zip(&got.parts) {
        if (w.part, &w.kernel, w.accepted) != (g.part, &g.kernel, g.accepted) {
            out.push(format!(
                "{label}: part {} path changed: {w:?} -> {g:?}",
                w.part
            ));
            continue;
        }
        let rel = (g.tets as f64 - w.tets as f64).abs() / (w.tets as f64).max(1.0);
        if rel > tets
            || (g.min_dihedral - w.min_dihedral).abs() > dihedral
            || (g.mean_non_orthogonality - w.mean_non_orthogonality).abs() > non_orth
        {
            out.push(format!(
                "{label}: part {} metrics out of band: {w:?} -> {g:?}",
                w.part
            ));
        }
    }
    out
}

/// Compare every run with the golden file, or rewrite it.
fn check_golden(runs: &[(String, Golden)], failures: &mut Vec<String>) {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN_FILE);
    if std::env::var_os("UPDATE_GOLDEN").is_some() {
        let mut text = String::from(
            "# Golden reports for the synthetic suite (tests/synthetic/tests/suite.rs).\n\
             # label <TAB> code severity [action/outcome] xcount; ... <TAB> part kernel accepted tets min_dihedral mean_non_orthogonality; ...\n\
             # Regenerate with UPDATE_GOLDEN=1 cargo test -p mesher-synthetic\n",
        );
        for (label, g) in runs {
            text.push_str(&golden_line(label, g));
            text.push('\n');
        }
        std::fs::write(&path, text).unwrap();
        return;
    }
    let Ok(text) = std::fs::read_to_string(&path) else {
        failures.push(format!(
            "{} missing; run UPDATE_GOLDEN=1 cargo test -p mesher-synthetic",
            path.display()
        ));
        return;
    };
    let golden = parse_golden(&text);
    for (label, got) in runs {
        match golden.get(label) {
            Some(want) => failures.extend(golden_diff(label, want, got)),
            None => failures.push(format!("{label}: not in the golden file")),
        }
    }
    for label in golden.keys() {
        if !runs.iter().any(|(l, _)| l == label) {
            failures.push(format!("{label}: in the golden file but no longer run"));
        }
    }
}

#[test]
fn synthetic_suite() {
    let mut tally: BTreeMap<DefectCode, Tally> = BTreeMap::new();
    let mut goldens: Vec<(String, Golden)> = Vec::new();
    let mut failures = Vec::new();
    let (mut meshed_parts, mut accepted_parts, mut delaunay_parts) = (0u32, 0u32, 0u32);
    let mut exact = ExactTally::default();
    let variants = [(1.0, 0.0), (1e-3, 1e3), (1e4, -5e4)];

    for c in cases() {
        let runs: &[(f64, f64)] = if c.variants {
            &variants
        } else {
            &variants[..1]
        };
        for &(scale, offset) in runs {
            let label = format!("{} [x{scale} +{offset}]", c.name);
            let started = std::time::Instant::now();
            let outcome = run(&c.scene.transformed(scale, offset), &c.parts);
            if std::env::var("SUITE_TIMING").is_ok() {
                eprintln!("{:>8.2?} {label}", started.elapsed());
            }
            if let Some(n) = c.parts_count {
                if outcome.parts != n {
                    failures.push(format!("{label}: {} parts, expected {n}", outcome.parts));
                }
            }
            goldens.push((label.clone(), golden_of(&outcome)));
            check_repair(&c, &outcome, &label, &mut failures);
            let t = check_meshes(&c, &outcome, &label, &mut failures);
            meshed_parts += t.meshed;
            accepted_parts += t.accepted;
            delaunay_parts += t.delaunay;
            check_exact(&c, &outcome, &label, &mut failures, &mut exact);
            let mut detected: BTreeMap<DefectCode, Vec<&Record>> = BTreeMap::new();
            for r in outcome
                .records
                .iter()
                .filter(|r| matches!(r.stage, Stage::Diagnose | Stage::Select))
            {
                detected.entry(r.code).or_default().push(r);
            }
            let codes: BTreeSet<DefectCode> = detected
                .keys()
                .copied()
                .chain(c.expect.iter().map(|e| e.0))
                .collect();
            for code in codes {
                let got = detected.remove(&code).unwrap_or_default();
                let want: Vec<&Key> = c
                    .expect
                    .iter()
                    .filter(|e| e.0 == code)
                    .map(|e| &e.1)
                    .collect();
                let t = tally.entry(code).or_default();
                if want.contains(&&Key::Any) {
                    if got.is_empty() {
                        t.fn_ += 1;
                        failures.push(format!("{label}: missed {}", code.as_str()));
                    } else {
                        t.tp += 1;
                    }
                    continue;
                }
                let mut keys: Vec<Key> = got
                    .iter()
                    .map(|r| match want.first() {
                        Some(Key::Vertices(_)) => Key::Vertices(r.location.vertices.clone()),
                        _ => Key::Faces(r.location.faces.clone()),
                    })
                    .collect();
                for w in want {
                    let w = &capped(w);
                    match keys.iter().position(|k| k == w) {
                        Some(i) => {
                            keys.remove(i);
                            t.tp += 1;
                        }
                        None => {
                            t.fn_ += 1;
                            failures.push(format!(
                                "{label}: missed {} at {w:?}; got {keys:?}",
                                code.as_str()
                            ));
                        }
                    }
                }
                for k in keys {
                    t.fp += 1;
                    failures.push(format!("{label}: unexpected {} at {k:?}", code.as_str()));
                }
            }
        }
    }

    eprintln!(
        "\n{:<30} {:>4} {:>4} {:>4} {:>10} {:>8}",
        "code", "tp", "fp", "fn", "precision", "recall"
    );
    for (code, t) in &tally {
        let p = f64::from(t.tp) / f64::from((t.tp + t.fp).max(1));
        let r = f64::from(t.tp) / f64::from((t.tp + t.fn_).max(1));
        eprintln!(
            "{:<30} {:>4} {:>4} {:>4} {:>10.3} {:>8.3}",
            code.as_str(),
            t.tp,
            t.fp,
            t.fn_,
            p,
            r
        );
    }

    eprintln!(
        "\nvolume: {meshed_parts} parts meshed ({delaunay_parts} conformal Delaunay, the rest stuffing), {accepted_parts} accepted by the volume gate"
    );
    eprintln!(
        "exact kernel: {} parts with every constraint recovered, {} Steiner points on segments, {} on faces, {} facets gift-wrapped, {:.2} s",
        exact.parts, exact.steiner_segments, exact.steiner_faces, exact.wrapped, exact.seconds
    );
    eprintln!(
        "rounding to f64: {} of {} parts valid as rounded, all after {} untangles",
        exact.parts - exact.rounding_invalid,
        exact.parts,
        exact.untangles
    );

    check_golden(&goldens, &mut failures);

    // Every diagnosable code is exercised (capacity overflow cannot be synthesized).
    let covered: BTreeSet<DefectCode> = tally
        .iter()
        .filter(|(_, t)| t.tp > 0)
        .map(|(c, _)| *c)
        .collect();
    for &code in DefectCode::ALL {
        let family = code.as_str().split('_').next().unwrap_or("");
        if ["GEO", "TOP", "INT", "PART", "QUAL"].contains(&family)
            && code != DefectCode::GeoCapacityExceeded
            && !covered.contains(&code)
        {
            failures.push(format!("no case exercises {}", code.as_str()));
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
