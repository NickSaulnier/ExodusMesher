//! Part selection (`--parts`) and checks between selected parts.

use std::collections::BTreeSet;

use mesher_geom::{Bvh, Contact, FastWindingNumber, TriangleSoup, classify_contact};

use super::parts::centroid;
use super::{Diagnosis, InternalError};
use crate::model::SurfaceModel;
use crate::overlay::Overlay;
use crate::report::{Action, DefectCode, JobConfig, Outcome, Record, RecordLog, Stage};

/// Resolve `config.parts` (names, or 1-based indices) to part positions. An
/// empty list selects every part. Records `PART_EMPTY_SELECTION` and returns
/// `None` if a requested part does not exist or nothing is selected.
pub fn select_parts(
    diagnosis: &Diagnosis,
    config: &JobConfig,
    log: &mut RecordLog,
) -> Option<Vec<usize>> {
    let parts = &diagnosis.parts;
    let available = || -> Vec<String> { parts.iter().map(|p| p.name.clone()).collect() };
    if config.parts.is_empty() {
        if parts.is_empty() {
            log.push(
                Record::new(DefectCode::PartEmptySelection, Stage::Select)
                    .measured("available", available())
                    .action(Action::Abort, Outcome::Applied)
                    .hint("the model has no parts"),
            );
            return None;
        }
        return Some((0..parts.len()).collect());
    }
    let mut chosen = BTreeSet::new();
    let mut missing = Vec::new();
    for wanted in &config.parts {
        let by_name = parts.iter().position(|p| &p.name == wanted);
        let by_index = || {
            wanted
                .parse::<usize>()
                .ok()
                .filter(|&i| (1..=parts.len()).contains(&i))
                .map(|i| i - 1)
        };
        match by_name.or_else(by_index) {
            Some(i) => {
                chosen.insert(i);
            }
            None => missing.push(wanted.clone()),
        }
    }
    if !missing.is_empty() {
        log.push(
            Record::new(DefectCode::PartEmptySelection, Stage::Select)
                .measured("missing", missing)
                .measured("available", available())
                .action(Action::Abort, Outcome::Applied)
                .hint("select parts by name or by 1-based index"),
        );
        return None;
    }
    Some(chosen.into_iter().collect())
}

/// Check every pair of selected parts for overlapping volumes
/// (`PART_OVERLAP`) and shared interfaces (`PART_SHARED_INTERFACE`).
pub fn check_part_pairs(
    model: &SurfaceModel,
    diagnosis: &Diagnosis,
    selected: &[usize],
    log: &mut RecordLog,
    overlay: &mut Overlay,
) -> Result<(), InternalError> {
    let points = &model.vertices;
    let parts = &diagnosis.parts;
    let soups: Vec<Option<TriangleSoup<'_>>> = parts
        .iter()
        .map(|p| {
            if p.faces.is_empty() {
                Ok(None)
            } else {
                TriangleSoup::new(points, &p.faces)
                    .map(Some)
                    .map_err(|e| InternalError(format!("part {}: {e}", p.name)))
            }
        })
        .collect::<Result<_, _>>()?;
    let bvhs: Vec<Option<Bvh<'_>>> = soups.iter().map(|s| s.map(Bvh::new)).collect();

    for (k, &i) in selected.iter().enumerate() {
        for &j in &selected[k + 1..] {
            let (Some(bi), Some(bj)) = (&bvhs[i], &bvhs[j]) else {
                continue;
            };
            if !bi.bounds().overlaps(&bj.bounds()) {
                continue;
            }
            let (pi, pj) = (&parts[i], &parts[j]);
            let mut overlapping = Vec::new();
            let mut interface = Vec::new();
            let mut any_contact = false;
            for (f, g) in bi.overlaps_with(bj) {
                let (tf, tg) = (pi.faces[f as usize], pj.faces[g as usize]);
                let pos = |t: [u32; 3]| t.map(|v| points[v as usize]);
                match classify_contact(tf, pos(tf), tg, pos(tg)) {
                    Contact::None => {}
                    Contact::Touching => any_contact = true,
                    Contact::Crossing | Contact::CoplanarAligned => {
                        any_contact = true;
                        overlapping.push((f, g));
                    }
                    Contact::CoplanarOpposed => {
                        any_contact = true;
                        interface.push((f, g));
                    }
                }
            }
            // No contact at all: one part may still lie wholly inside the other.
            let contained = !any_contact && (inside(points, pi, pj) || inside(points, pj, pi));
            let faces = |pairs: &[(u32, u32)]| -> Vec<u32> {
                pairs
                    .iter()
                    .flat_map(|&(f, g)| [pi.face_ids[f as usize], pj.face_ids[g as usize]])
                    .collect()
            };
            let line = |t: &u32| model.triangle_line(*t);
            if !overlapping.is_empty() || contained {
                let fs = faces(&overlapping);
                log.push(
                    Record::new(DefectCode::PartOverlap, Stage::Select)
                        .part(&pi.name)
                        .faces(fs.iter().copied())
                        .lines(fs.iter().map(line))
                        .measured("other_part", pj.name.clone())
                        .measured("intersecting_pairs", overlapping.len() as u64)
                        .measured("contained", contained)
                        .hint("the volumes of these parts overlap"),
                );
                overlay.faces(fs, DefectCode::PartOverlap);
            } else if !interface.is_empty() {
                let fs = faces(&interface);
                log.push(
                    Record::new(DefectCode::PartSharedInterface, Stage::Select)
                        .part(&pi.name)
                        .faces(fs.iter().copied())
                        .lines(fs.iter().map(line))
                        .measured("other_part", pj.name.clone())
                        .measured("interface_pairs", interface.len() as u64)
                        .hint("meshed independently in v1; the meshes will not be conformal across the interface"),
                );
                overlay.faces(fs, DefectCode::PartSharedInterface);
            }
        }
    }
    Ok(())
}

/// Is a point of part `a` inside part `b` (by `b`'s winding number)?
fn inside(points: &[[f64; 3]], a: &super::Part, b: &super::Part) -> bool {
    let (Some(&first), Ok(soup)) = (a.faces.first(), TriangleSoup::new(points, &b.faces)) else {
        return false;
    };
    FastWindingNumber::new(soup)
        .winding_number(centroid(points, first))
        .abs()
        > 0.5
}
