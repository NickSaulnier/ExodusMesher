//! Diagnosis: a read-only pass that splits the model into parts and records
//! every defect per part. Nothing is changed here; each record's action is
//! `Deferred` to the repair stage.
//!
//! Topology is judged on the surface as repair would leave it after its least
//! invasive steps: duplicate vertices merged (by tolerance), degenerate and
//! duplicate faces set aside. That way one root cause gives one record (a
//! duplicated vertex is reported as such, not as a ring of boundary edges).

mod part;
mod parts;
mod select;
mod uf;
mod weld;

use mesher_geom::Aabb;

use crate::model::SurfaceModel;
use crate::overlay::Overlay;
use crate::report::{DefectCode, PartSource, PartSummary, Record, RecordLog, Stage};

pub use select::{check_part_pairs, select_parts};

/// Diagnosis thresholds. Lengths are relative to the model's bounding-box
/// diagonal; angles and ratios are scale-free.
#[derive(Debug, Clone, PartialEq)]
pub struct Tolerances {
    /// Vertices closer than this × diagonal are duplicates.
    pub merge_relative: f64,
    /// Faces with a smaller angle are slivers.
    pub sliver_angle_degrees: f64,
    /// A shell whose closed-up volume is at most this × area^1.5 encloses nothing.
    pub open_sheet_volume_ratio: f64,
    /// Bounding-box diagonals outside this range suggest wrong units.
    pub suspect_scale: (f64, f64),
}

impl Default for Tolerances {
    fn default() -> Tolerances {
        Tolerances {
            merge_relative: 1e-7,
            sliver_angle_degrees: 1.0,
            open_sheet_volume_ratio: 1e-6,
            suspect_scale: (1e-4, 1e6),
        }
    }
}

/// One part after diagnosis.
#[derive(Debug, Clone)]
pub struct Part {
    pub name: String,
    pub source: PartSource,
    /// Every triangle of the part (indices into the model), ascending.
    pub triangles: Vec<u32>,
    /// Faces kept for meshing (not degenerate, not duplicates), consistently
    /// oriented by majority, as merged-vertex ids.
    pub faces: Vec<[u32; 3]>,
    /// Model triangle index of each kept face.
    pub face_ids: Vec<u32>,
    /// Distinct merged vertices.
    pub vertices: usize,
    pub shells: usize,
    pub bounds: Aabb,
    /// A fatal defect means the part cannot be meshed.
    pub fatal: bool,
}

#[derive(Debug, Clone)]
pub struct Diagnosis {
    pub parts: Vec<Part>,
    /// Merged-vertex id of every model vertex (the smallest index in its cluster).
    pub rep: Vec<u32>,
    /// Absolute merge tolerance used.
    pub merge_tolerance: f64,
    /// Bounding-box diagonal of the referenced vertices.
    pub diagonal: f64,
    pub overlay: Overlay,
}

impl Diagnosis {
    /// Part summaries for the report; `selected` holds part positions.
    #[must_use]
    pub fn summaries(&self, selected: &[usize]) -> Vec<PartSummary> {
        self.parts
            .iter()
            .enumerate()
            .map(|(i, p)| PartSummary {
                index: i as u32 + 1,
                name: p.name.clone(),
                source: p.source,
                triangles: p.triangles.len() as u64,
                vertices: p.vertices as u64,
                shells: p.shells as u64,
                bounds: if p.bounds.is_empty() {
                    [[0.0; 3]; 2]
                } else {
                    [p.bounds.min, p.bounds.max]
                },
                selected: selected.contains(&i),
            })
            .collect()
    }
}

/// A failure that is a bug, not bad input: the worker reports it as an
/// internal error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InternalError(pub String);

/// Diagnose every part of the model.
pub fn diagnose(
    model: &SurfaceModel,
    tol: &Tolerances,
    log: &mut RecordLog,
) -> Result<Diagnosis, InternalError> {
    let mut referenced = vec![false; model.vertices.len()];
    for t in &model.triangles {
        for &v in &t.vertices {
            if let Some(r) = referenced.get_mut(v as usize) {
                *r = true;
            }
        }
    }
    let bounds = model
        .vertices
        .iter()
        .zip(&referenced)
        .filter(|(_, r)| **r)
        .fold(Aabb::EMPTY, |b, (p, _)| b.with_point(*p));
    let diagonal = bounds.diagonal();
    let merge_tolerance = tol.merge_relative * diagonal;
    let rep = weld::weld(&model.vertices, &referenced, merge_tolerance);

    if !bounds.is_empty() && (diagonal < tol.suspect_scale.0 || diagonal > tol.suspect_scale.1) {
        log.push(
            Record::new(DefectCode::QualSuspectScale, Stage::Diagnose)
                .measured("bounding_box_diagonal", diagonal)
                .threshold("min_diagonal", tol.suspect_scale.0)
                .threshold("max_diagonal", tol.suspect_scale.1)
                .hint("the model's size suggests unexpected units; tolerances are relative, so this is informational"),
        );
    }

    let groups = parts::split(model, &rep)?;
    let ctx = part::Context {
        model,
        rep: &rep,
        tol,
        merge_tolerance,
        diagonal,
    };
    let mut overlay = Overlay::default();
    let mut out = Vec::with_capacity(groups.len());
    for group in groups {
        out.push(part::analyze(&ctx, group, log, &mut overlay)?);
    }
    Ok(Diagnosis {
        parts: out,
        rep,
        merge_tolerance,
        diagonal,
        overlay,
    })
}
