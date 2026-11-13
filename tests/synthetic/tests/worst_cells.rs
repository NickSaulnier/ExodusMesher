//! Milestone 15: the worst cells are reported and located. A cube with a
//! sliver face, meshed on its own triangulation (`--no-remesh`), keeps poor
//! cells at the sliver; each measure beyond its limit gets a
//! `VOL_LOW_QUALITY` warning located at the worst cell, and the summary
//! locates the worst of every measure. Far beyond the limits, stuffing is
//! kept instead if its cells are better.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::repair::repair;
use mesher_core::volume::{MeshOutcome, mesh_parts};
use mesher_core::{
    Action, AttemptOutcome, DefectCode, JobConfig, Kernel, Outcome, QualityMeasure, Record,
    RecordLog, Severity,
};
use mesher_io::parse_obj;
use mesher_synthetic::{Scene, Solid};
use mesher_tet::fv_quality;

/// The cube with a sliver face, and the sliver's apex.
fn sliver_cube() -> (Scene, [f64; 3]) {
    sliver_cube_at(0.002)
}

/// The cube with a sliver face whose apex is `t` of the way from the
/// middle of its long edge to the opposite corner.
fn sliver_cube_at(t: f64) -> (Scene, [f64; 3]) {
    let mut s = Solid::unit_cube();
    let [a, b, c] = s.tris[0];
    let (pa, pb, pc) = (
        s.points[a as usize],
        s.points[b as usize],
        s.points[c as usize],
    );
    let mid = [0, 1, 2].map(|k| (pa[k] + pb[k]) / 2.0);
    s.points
        .push([0, 1, 2].map(|k| mid[k] + t * (pc[k] - mid[k])));
    let m = s.points.len() as u32 - 1;
    s.tris[0] = [a, b, m];
    s.tris.push([b, c, m]);
    s.tris.push([c, a, m]);
    let apex = s.points[m as usize];
    (Scene::single(s), apex)
}

fn mesh(scene: &Scene, no_remesh: bool) -> (Vec<Record>, MeshOutcome) {
    mesh_at(scene, no_remesh, 0.06)
}

fn mesh_at(scene: &Scene, no_remesh: bool, element_size: f64) -> (Vec<Record>, MeshOutcome) {
    let obj = scene.to_obj();
    let mut log = RecordLog::new();
    let model = parse_obj(obj.as_bytes(), &mut log, &mut |_| {});
    let diagnosis = diagnose(&model, &Tolerances::default(), &mut log).unwrap();
    let config = JobConfig {
        element_size,
        no_remesh,
        ..JobConfig::default()
    };
    let selected = select_parts(&diagnosis, &config, &mut log).unwrap();
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
    let outcome = mesh_parts(&model, &repaired, &config, &mut log).unwrap();
    (log.records().to_vec(), outcome)
}

#[test]
fn worst_cells_are_located() {
    let (scene, apex) = sliver_cube();
    let (records, outcome) = mesh(&scene, true);
    let part = outcome.meshed.first().expect("the part is meshed");
    let (_, summary) = outcome.summaries.first().expect("a summary");
    assert!(summary.accepted && !summary.remeshed);
    let q = fv_quality(&part.points, &part.tets);
    let centroid = |e: u64| {
        let t = part.tets[e as usize - 1];
        let p = t.map(|v| part.points[v as usize]);
        [0, 1, 2].map(|k| (p[0][k] + p[1][k] + p[2][k] + p[3][k]) / 4.0)
    };

    // The summary: every measure, overall and at the walls, in range.
    let measures: Vec<(QualityMeasure, bool)> = summary
        .worst_cells
        .iter()
        .map(|w| (w.measure, w.wall))
        .collect();
    assert_eq!(
        measures,
        [
            (QualityMeasure::MinDihedral, false),
            (QualityMeasure::NonOrthogonality, false),
            (QualityMeasure::Skewness, false),
            (QualityMeasure::NonOrthogonality, true),
            (QualityMeasure::Skewness, true),
        ]
    );
    let tets = part.tets.len() as u64;
    for w in &summary.worst_cells {
        assert!((1..=tets).contains(&w.element), "{w:?}");
        assert!(w.neighbour.is_none_or(|n| (1..=tets).contains(&n)), "{w:?}");
    }
    let dihedral = &summary.worst_cells[0];
    assert_eq!(dihedral.value, summary.min_dihedral_degrees);
    assert_eq!(dihedral.point, centroid(dihedral.element));
    assert_eq!(
        summary.worst_cells[1].value,
        summary.max_non_orthogonality_degrees
    );
    assert_eq!(summary.worst_cells[2].value, summary.max_skewness);
    assert_eq!(
        summary.worst_cells[3].value,
        q.wall_max_non_orthogonality_degrees
    );

    // The warnings: at least one measure is beyond its limit at the sliver,
    // and every quality warning is located at the worst cell.
    let warnings: Vec<&Record> = records
        .iter()
        .filter(|r| r.code == DefectCode::VolLowQuality && r.severity == Severity::Warning)
        .collect();
    assert!(!warnings.is_empty(), "{records:#?}");
    for r in &warnings {
        let point = r.location.point.expect("located");
        let element = r.measured["element"].as_u64().expect("an element");
        let w = summary
            .worst_cells
            .iter()
            .find(|w| !w.wall && w.element == element && w.point == point)
            .expect("the record's cell is a summary's worst");
        // The worst is at the sliver, within an element size of its apex.
        let d = (0..3)
            .map(|k| (point[k] - apex[k]).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(d < 0.06, "{w:?} is {d} from the sliver's apex {apex:?}");
    }
}

/// Remeshing (the default) lifts the same cube within the limits: no
/// warnings, and the summary still locates its worst cells.
#[test]
fn a_remeshed_sliver_cube_has_no_quality_warning() {
    let (records, outcome) = mesh(&sliver_cube().0, false);
    let (_, summary) = outcome.summaries.first().expect("a summary");
    assert!(summary.accepted && summary.remeshed);
    assert_eq!(summary.worst_cells.len(), 5);
    assert!(
        !records
            .iter()
            .any(|r| r.code == DefectCode::VolLowQuality && r.severity == Severity::Warning),
        "{records:#?}"
    );
}

/// A sliver far thinner, meshed on its own triangulation: the conformal
/// mesh's worst cell is far beyond the limits, so stuffing is tried and
/// kept, with a `VOL_LOW_QUALITY` warning (`fallback` / `applied`) giving
/// both worst cells, and attempts saying the conformal mesh was not kept.
#[test]
fn far_beyond_the_limits_stuffing_is_kept() {
    // Fine enough that stuffing keeps the cube's volume to 1 % despite
    // rounding its sharp edges.
    let (records, outcome) = mesh_at(&sliver_cube_at(1e-5).0, true, 0.03);
    let (_, summary) = outcome.summaries.first().expect("a summary");
    assert_eq!(summary.kernel, Kernel::Stuffing);
    assert!(summary.accepted);
    let part = outcome.meshed.first().expect("meshed");
    assert_eq!(part.kernel, Kernel::Stuffing);
    let kinds: Vec<(Kernel, AttemptOutcome)> = outcome
        .attempts
        .iter()
        .map(|a| (a.kernel, a.outcome))
        .collect();
    assert_eq!(
        kinds,
        [
            (Kernel::Delaunay, AttemptOutcome::Failed),
            (Kernel::Stuffing, AttemptOutcome::Succeeded)
        ]
    );
    let switched: Vec<&Record> = records
        .iter()
        .filter(|r| {
            r.code == DefectCode::VolLowQuality
                && r.action == Some(Action::Fallback)
                && r.outcome == Some(Outcome::Applied)
        })
        .collect();
    assert_eq!(switched.len(), 1, "{records:#?}");
    let m = &switched[0].measured;
    let (conformal, stuffing) = (
        m["conformal_worst_cell"].as_f64().unwrap(),
        m["stuffing_worst_cell"].as_f64().unwrap(),
    );
    assert!(conformal < 0.25 && stuffing > conformal, "{m:?}");
    assert!(m["volume_error"].as_f64().unwrap() <= 0.01, "{m:?}");
    // No quality warning is left about the conformal mesh, which was not
    // written; the stuffing mesh is within the limits.
    assert!(
        !records.iter().any(|r| r.code == DefectCode::VolLowQuality
            && r.severity == Severity::Warning
            && r.action.is_none()),
        "{records:#?}"
    );
}
