//! Milestone 3 exit criterion: every injected defect is detected with its exact
//! location, and clean shapes produce no findings, at several scales and
//! offsets. Prints precision and recall per defect code.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::vec_init_then_push
)]

use std::collections::{BTreeMap, BTreeSet};

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::report::MAX_LOCATION_ITEMS;
use mesher_core::{DefectCode, JobConfig, Record, RecordLog};
use mesher_io::parse_obj;
use mesher_synthetic::{Scene, Solid};

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
}

fn case(name: &'static str, scene: Scene, expect: Vec<(DefectCode, Key)>) -> Case {
    Case {
        name,
        scene,
        parts: Vec::new(),
        expect,
        parts_count: None,
        variants: true,
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
        out.push(case(
            "hole, open bowl",
            Scene::single(bowl),
            vec![(TopBoundaryEdge, Key::Any)],
        ));
    }
    {
        let s = cube().union(&Solid::cuboid([1.0, 1.0, 0.0], [2.0, 2.0, 1.0]));
        let e = [[1.0, 1.0, 0.0], [1.0, 1.0, 1.0]].map(|p| s.find(p).unwrap());
        out.push(case(
            "non-manifold edge (cubes sharing an edge)",
            Scene::single(s),
            vec![(TopNonmanifoldEdge, verts(e))],
        ));
    }
    {
        let s = cube().union(&Solid::cuboid([1.0; 3], [2.0; 3]));
        let v = s.find([1.0; 3]).unwrap();
        out.push(Case {
            parts_count: Some(1),
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

struct Outcome {
    records: Vec<Record>,
    parts: usize,
}

fn run(scene: &Scene, parts: &[&str]) -> Outcome {
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
    if let Some(selected) = select_parts(&diagnosis, &config, &mut log) {
        let mut overlay = diagnosis.overlay.clone();
        check_part_pairs(&model, &diagnosis, &selected, &mut log, &mut overlay).unwrap();
    }
    Outcome {
        records: log.records().to_vec(),
        parts: diagnosis.parts.len(),
    }
}

#[derive(Default, Debug, Clone, Copy)]
struct Tally {
    tp: u32,
    fp: u32,
    fn_: u32,
}

#[test]
fn synthetic_suite() {
    let mut tally: BTreeMap<DefectCode, Tally> = BTreeMap::new();
    let mut failures = Vec::new();
    let variants = [(1.0, 0.0), (1e-3, 1e3), (1e4, -5e4)];

    for c in cases() {
        let runs: &[(f64, f64)] = if c.variants {
            &variants
        } else {
            &variants[..1]
        };
        for &(scale, offset) in runs {
            let label = format!("{} [x{scale} +{offset}]", c.name);
            let outcome = run(&c.scene.transformed(scale, offset), &c.parts);
            if let Some(n) = c.parts_count {
                if outcome.parts != n {
                    failures.push(format!("{label}: {} parts, expected {n}", outcome.parts));
                }
            }
            let mut detected: BTreeMap<DefectCode, Vec<&Record>> = BTreeMap::new();
            for r in &outcome.records {
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
