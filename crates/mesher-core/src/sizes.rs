//! Size controls per part and per patch (Milestone 18): `JobConfig::sizes`.
//!
//! A control naming a part and no group sets that part's element size (the
//! last one listed wins over `element_size`). One naming a group sets the
//! size on that group's faces, in the named part or in every part, and the
//! size field grades away from them; it only refines. A control that
//! matches no selected part, names a group with no faces there, or asks a
//! patch for more than its part's size is a `SIZE_CONTROL_IGNORED` warning.

use crate::features::face_groups;
use crate::model::SurfaceModel;
use crate::repair::RepairedPart;
use crate::report::{DefectCode, JobConfig, Record, RecordLog, SizeControl, Stage};

/// Whether a control's part matches `r` (by name or 1-based index).
fn names(control: &SizeControl, r: &RepairedPart) -> bool {
    control
        .part
        .as_deref()
        .is_none_or(|p| p == r.name || p.parse::<usize>().is_ok_and(|i| i == r.part + 1))
}

/// The part's element size, relative to its bounding-box diagonal.
#[must_use]
pub fn part_size(r: &RepairedPart, config: &JobConfig) -> f64 {
    config
        .sizes
        .iter()
        .rev()
        .find(|c| c.group.is_none() && c.part.is_some() && names(c, r))
        .map_or(config.element_size, |c| c.element_size)
}

/// The size each face's patch asks for (absolute; infinite where none, or
/// where it would not be below the part's size `h`).
#[must_use]
pub fn face_sizes(
    model: &SurfaceModel,
    r: &RepairedPart,
    config: &JobConfig,
    diagonal: f64,
    h: f64,
) -> Vec<f64> {
    let controls: Vec<(u32, f64)> = config
        .sizes
        .iter()
        .filter(|c| names(c, r))
        .filter_map(|c| {
            let g = c.group.as_deref()?;
            let id = model.groups.iter().position(|n| n == g)? as u32;
            Some((id, c.element_size * diagonal))
        })
        .filter(|&(_, s)| s < h)
        .collect();
    if controls.is_empty() {
        return Vec::new();
    }
    face_groups(model, r)
        .iter()
        .map(|&g| {
            controls
                .iter()
                .filter(|c| c.0 == g)
                .map(|c| c.1)
                .fold(f64::INFINITY, f64::min)
        })
        .collect()
}

/// Record the controls that took no effect.
pub fn record_ignored_controls(
    model: &SurfaceModel,
    repaired: &[RepairedPart],
    config: &JobConfig,
    log: &mut RecordLog,
) {
    for c in &config.sizes {
        let parts: Vec<&RepairedPart> = repaired.iter().filter(|r| names(c, r)).collect();
        let reason = if parts.is_empty() || (c.part.is_none() && c.group.is_none()) {
            Some(if c.part.is_none() && c.group.is_none() {
                "it names neither a part nor a group (use `element_size` for every part)"
            } else {
                "no selected part has this name or index"
            })
        } else if let Some(g) = c.group.as_deref() {
            let id = model.groups.iter().position(|n| n == g).map(|i| i as u32);
            let mut found = false;
            let mut refines = false;
            for r in &parts {
                if id.is_some_and(|id| face_groups(model, r).contains(&id)) {
                    found = true;
                    let diagonal = mesher_geom::Aabb::of_points(&r.points).diagonal();
                    refines |= c.element_size < part_size(r, config) || diagonal <= 0.0;
                }
            }
            if !found {
                Some("no face of the selected parts is in this group")
            } else if !refines {
                Some("a patch's size cannot exceed its part's element size")
            } else {
                None
            }
        } else {
            None
        };
        if let Some(reason) = reason {
            let mut rec = Record::new(DefectCode::SizeControlIgnored, Stage::Sizing)
                .measured("element_size", c.element_size)
                .measured("reason", reason)
                .hint("the size control took no effect");
            if let Some(p) = &c.part {
                rec = rec.measured("part", p.clone());
            }
            if let Some(g) = &c.group {
                rec = rec.measured("group", g.clone());
            }
            log.push(rec);
        }
    }
}
