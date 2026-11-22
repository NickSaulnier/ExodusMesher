//! Milestone 18, stage 3: size controls per part and per patch
//! (`JobConfig::sizes`), through the whole pipeline. A part's control sets
//! its element size; a patch's control refines the mesh at that patch and
//! grades away from it; controls that take no effect are
//! `SIZE_CONTROL_IGNORED` warnings.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::repair::repair;
use mesher_core::volume::{MeshOutcome, mesh_parts};
use mesher_core::{DefectCode, JobConfig, MeshSummary, Record, RecordLog, SizeControl};
use mesher_io::parse_obj;
use mesher_synthetic::{Scene, Solid};

fn mesh(scene: &Scene, config: &JobConfig) -> (MeshOutcome, Vec<Record>) {
    let obj = scene.to_obj();
    let mut log = RecordLog::new();
    let model = parse_obj(obj.as_bytes(), &mut log, &mut |_| {});
    let diagnosis = diagnose(&model, &Tolerances::default(), &mut log).unwrap();
    let selected = select_parts(&diagnosis, config, &mut log).unwrap();
    let mut overlay = diagnosis.overlay.clone();
    check_part_pairs(&model, &diagnosis, &selected, &mut log, &mut overlay).unwrap();
    let repaired = repair(
        &model,
        &diagnosis,
        &selected,
        &Tolerances::default(),
        &mut log,
    )
    .unwrap();
    let outcome = mesh_parts(&model, &repaired, config, &mut log).unwrap();
    (outcome, log.records().to_vec())
}

fn summary(outcome: &MeshOutcome, part: usize) -> &MeshSummary {
    &outcome
        .summaries
        .iter()
        .find(|(i, _)| *i == part)
        .unwrap()
        .1
}

/// Two unit cubes, `a` and `b`; `b`'s top face (z = 1) is the group `lid`,
/// its other faces `side`.
fn scene() -> Scene {
    let b = Solid::cuboid([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
    // `cuboid` writes z = min first, then z = max, two triangles each.
    let groups = (0..b.tris.len())
        .map(|t| if (2..4).contains(&t) { "lid" } else { "side" }.to_owned())
        .collect();
    Scene::default()
        .with(Some("a"), Solid::unit_cube())
        .with_groups(Some("b"), b, groups)
}

fn control(part: Option<&str>, group: Option<&str>, element_size: f64) -> SizeControl {
    SizeControl {
        part: part.map(str::to_owned),
        group: group.map(str::to_owned),
        element_size,
    }
}

fn base() -> JobConfig {
    JobConfig {
        element_size: 0.2,
        ..JobConfig::default()
    }
}

#[test]
fn a_part_control_sets_its_element_size() {
    let (plain, _) = mesh(&scene(), &base());
    let config = JobConfig {
        sizes: vec![
            control(Some("b"), None, 0.3),
            // The last control for a part wins; by 1-based index too.
            control(Some("2"), None, 0.1),
        ],
        ..base()
    };
    let (sized, records) = mesh(&scene(), &config);
    let diagonal = 3f64.sqrt();
    // Part a keeps the job's size; b takes 0.1 of its diagonal.
    assert_eq!(
        summary(&sized, 0).element_size,
        summary(&plain, 0).element_size
    );
    assert!((summary(&plain, 1).element_size - 0.2 * diagonal).abs() < 1e-12);
    assert!((summary(&sized, 1).element_size - 0.1 * diagonal).abs() < 1e-12);
    assert!(summary(&sized, 1).tets > 2 * summary(&plain, 1).tets);
    assert!(
        !records
            .iter()
            .any(|r| r.code == DefectCode::SizeControlIgnored)
    );
}

#[test]
fn a_patch_control_refines_its_patch() {
    let (plain, _) = mesh(&scene(), &base());
    let config = JobConfig {
        sizes: vec![control(Some("b"), Some("lid"), 0.03)],
        ..base()
    };
    let (sized, records) = mesh(&scene(), &config);
    let (p, s) = (summary(&plain, 1), summary(&sized, 1));
    // The part's size stays; the lid is meshed far finer.
    assert_eq!(p.element_size, s.element_size);
    let min = s.min_element_size.unwrap();
    assert!((min - 0.03 * 3f64.sqrt()).abs() < 1e-9, "{min}");
    assert!(s.tets > 3 * p.tets, "{} against {}", s.tets, p.tets);
    // Part a has no lid and is unchanged.
    assert_eq!(summary(&sized, 0).tets, summary(&plain, 0).tets);
    assert!(
        !records
            .iter()
            .any(|r| r.code == DefectCode::SizeControlIgnored)
    );
}

#[test]
fn controls_without_effect_are_recorded() {
    let config = JobConfig {
        sizes: vec![
            control(Some("nope"), None, 0.1),
            control(Some("a"), Some("lid"), 0.05),
            control(Some("b"), Some("lid"), 0.5),
            control(None, None, 0.1),
        ],
        ..base()
    };
    let (outcome, records) = mesh(&scene(), &config);
    let reasons: Vec<String> = records
        .iter()
        .filter(|r| r.code == DefectCode::SizeControlIgnored)
        .map(|r| r.measured["reason"].as_str().unwrap().to_owned())
        .collect();
    assert_eq!(
        reasons,
        [
            "no selected part has this name or index",
            "no face of the selected parts is in this group",
            "a patch's size cannot exceed its part's element size",
            "it names neither a part nor a group (use `element_size` for every part)",
        ],
        "{reasons:?}"
    );
    // Nothing changed the meshes.
    let (plain, _) = mesh(&scene(), &base());
    for part in 0..2 {
        assert_eq!(summary(&outcome, part).tets, summary(&plain, part).tets);
    }
}
