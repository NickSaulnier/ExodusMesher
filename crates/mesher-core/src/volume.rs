//! Volume meshing and the volume gate.
//!
//! The fallback strategy of `PLAN.md`: a part whose repaired surface passed
//! the surface gate is meshed first by the conforming Delaunay kernel
//! (boundary-exact), or, with `--kernel exact`, by the exact constrained
//! kernel rounded to f64, the conforming kernel following if that fails. If that fails (budget, recovery) or the surface did not
//! pass the gate, isosurface stuffing on the winding-number field (`0.5 − w`)
//! takes over, unless `--strict` forbids it. A failed conformal attempt is a
//! warning when the fallback then succeeds, an error otherwise. Last, a part
//! whose surface failed the gate and whose stuffing was rejected (or not
//! allowed) is wrapped (repair rung 3) and meshed again from the start.
//!
//! The volume gate checks orientation, quality (the proven bounds for
//! stuffing; the sliver threshold for Delaunay) and the Hausdorff deviation
//! from the input surface, and records finite-volume quality metrics. Accepted
//! parts are assembled into one mesh: an element block and a `<part>_boundary`
//! side set per part, plus a side set per OBJ group on the boundary of a
//! conformal mesh.

use std::collections::{BTreeMap, HashMap};

use mesher_geom::bvh::closest_point_on_triangle;
use mesher_geom::predicates::orient3d;
use mesher_geom::vec3::{cross, dot, norm, scale, sub};
use mesher_geom::{Aabb, Bvh, FastWindingNumber, Point3, Sign, TriangleSoup};
use mesher_tet::cdt::{
    FaceError, RegionError, RoundingError, SegmentError, Untangle, mesh_region, round_region,
};
use mesher_tet::stuffing::{DIHEDRAL_MAX_DEGREES, DIHEDRAL_MIN_DEGREES};
use mesher_tet::{
    ConformingError, ConformingParams, FvQuality, GradeParams, ImproveParams, RemeshParams,
    SizeField, SizingParams, Stall, StuffingError, StuffingParams, WorstAt, conforming_mesh,
    fv_quality, grade, improve, remesh, stuff, stuff_graded,
};

use crate::diagnose::{InternalError, Tolerances};
use crate::mesh::{Block, SideSet, VolumeMesh, tet_side_nodes};
use crate::model::SurfaceModel;
use crate::repair::{RepairedPart, gate_failures, wrap_part};
use crate::report::{
    Action, Attempt, AttemptOutcome, DefectCode, JobConfig, Kernel, KernelChoice, MeshSummary,
    Outcome, QualityMeasure, Record, RecordLog, Severity, Stage, WorstCell,
};

/// Slivers below this dihedral angle in a conformal mesh earn a quality warning.
pub const SLIVER_DEGREES: f64 = 10.0;

/// Faces beyond this non-orthogonality earn a quality warning (the target
/// of `ImproveParams`).
pub const MAX_NON_ORTHOGONALITY_DEGREES: f64 = 70.0;

/// Faces beyond this skewness earn a quality warning (the target of
/// `ImproveParams`).
pub const MAX_SKEWNESS: f64 = 2.0;

/// The largest ratio between face neighbours' sizes (longest edges) before
/// `VOL_GRADATION` warns. Measured Nov 22, 2026: stuffing, whose lattice
/// halves its size at each level, never exceeds 2.85 (78 meshes); the
/// conformal meshes of the subsets exceed 3 on 93 of 197, 4 on 45 and 5 on
/// 20 (up to 31.8), where small cells crowd against large ones. Since then
/// conformal meshes are graded to it (`improve_cells`), so it warns only
/// where grading could not split a cell or ran out of budget.
pub const MAX_SIZE_RATIO: f64 = 4.0;

/// A conformal mesh whose worst cell (`worst_cell`) is below this is far
/// beyond the quality limits (a dihedral angle under 2.5°, a
/// non-orthogonality over 85° or a skewness over 8): stuffing is tried too,
/// and kept if its worst cell is better.
pub const FAR_WORST_CELL: f64 = 0.25;

/// Stuffing kept for its cells must hold the surface's volume to this
/// fraction: its boundary approximates the surface, and thin features it
/// misses would be lost (on the Thingi10k subset, 135674 by 12.5 %).
pub const FAR_MAX_VOLUME_ERROR: f64 = 0.01;

/// A part's accepted volume mesh.
#[derive(Debug, Clone)]
pub struct MeshedPart {
    /// Position of the part in `Diagnosis::parts`.
    pub part: usize,
    pub name: String,
    pub kernel: Kernel,
    pub points: Vec<Point3>,
    pub tets: Vec<[u32; 4]>,
    /// Boundary faces per OBJ group (conformal meshes only), by group name.
    pub face_groups: Vec<(String, Vec<[u32; 3]>)>,
}

#[derive(Debug, Clone, Default)]
pub struct MeshOutcome {
    /// Parts whose mesh passed the volume gate.
    pub meshed: Vec<MeshedPart>,
    /// The summary of each part's final attempt, by part position.
    pub summaries: Vec<(usize, MeshSummary)>,
    pub attempts: Vec<Attempt>,
    /// The parts wrapped by repair rung 3 (the last resort): their wrap
    /// and repair summary replace what repair gave.
    pub wrapped: Vec<RepairedPart>,
}

/// The size field of a part at element size `h`: graded from its curvature
/// and thickness, the feature angle (when set) for its curvature term too,
/// and its patches' size controls; `None` under `uniform_size` with no patch
/// sizes (the size is `h` everywhere).
fn size_field(
    model: &SurfaceModel,
    r: &RepairedPart,
    h: f64,
    config: &JobConfig,
) -> Option<SizeField> {
    let diagonal = Aabb::of_points(&r.points).diagonal();
    let face_sizes = crate::sizes::face_sizes(model, r, config, diagonal, h);
    let smallest = face_sizes.iter().copied().fold(f64::INFINITY, f64::min);
    if config.uniform_size && !smallest.is_finite() {
        return None;
    }
    let mut p = SizingParams::for_size(h);
    if let Some(a) = config.feature_angle {
        p.feature_degrees = a;
    }
    // A patch may ask for less than the usual floor.
    p.min_size = p.min_size.min(smallest);
    p.surface_terms = !config.uniform_size;
    if smallest.is_finite() {
        p.face_sizes = face_sizes;
    }
    Some(SizeField::from_surface(&r.points, &r.faces, &p))
}

/// Mesh every repaired part.
pub fn mesh_parts(
    model: &SurfaceModel,
    repaired: &[RepairedPart],
    config: &JobConfig,
    log: &mut RecordLog,
) -> Result<MeshOutcome, InternalError> {
    crate::features::record_ignored_overrides(model, repaired, config, log);
    crate::sizes::record_ignored_controls(model, repaired, config, log);
    let mut outcome = MeshOutcome::default();
    for r in repaired {
        mesh_part(model, r, config, log, &mut outcome)?;
    }
    Ok(outcome)
}

struct Candidate {
    kernel: Kernel,
    points: Vec<Point3>,
    tets: Vec<[u32; 4]>,
    /// Zero-volume tets the kernel dropped.
    dropped: u64,
    steiner: u64,
    face_groups: Vec<(String, Vec<[u32; 3]>)>,
    /// Smallest target size, when the size was graded.
    min_size: Option<f64>,
    /// Meshed from the remeshed surface.
    remeshed: bool,
}

fn mesh_part(
    model: &SurfaceModel,
    r: &RepairedPart,
    config: &JobConfig,
    log: &mut RecordLog,
    outcome: &mut MeshOutcome,
) -> Result<(), InternalError> {
    let _prof = mesher_tet::prof::scope("volume::mesh_part");
    let name = r.name.as_str();
    let internal = |e: String| InternalError(format!("part {name}: {e}"));
    let soup = TriangleSoup::new(&r.points, &r.faces).map_err(|e| internal(e.to_string()))?;
    let bounds = Aabb::of_points(&r.points);
    let diagonal = bounds.diagonal();
    let h = crate::sizes::part_size(r, config) * diagonal;
    let winding = FastWindingNumber::new(soup);
    let attempt = |kernel: Kernel, ok: bool| Attempt {
        part: Some(name.to_owned()),
        kernel,
        outcome: if ok {
            AttemptOutcome::Succeeded
        } else {
            AttemptOutcome::Failed
        },
        remeshed: false,
        stdout: String::new(),
        stderr: String::new(),
    };

    // 0. With `--kernel exact`: the exact constrained kernel, rounded to f64.
    let mut exact_failure: Option<Record> = None;
    if config.kernel == KernelChoice::Exact && r.summary.gate_passed && !config.stuffing_only {
        match exact(model, r, name) {
            Ok((c, fixes)) => {
                let (accepted, summary) =
                    volume_gate(name, &c, r, &winding, diagonal, h, config, log)?;
                outcome.attempts.push(attempt(Kernel::Exact, accepted));
                outcome.summaries.push((r.part, summary));
                // Every untangle is a record, whether or not the mesh is kept.
                for fix in fixes {
                    log.push(fix);
                }
                if accepted {
                    outcome.meshed.push(MeshedPart {
                        part: r.part,
                        name: name.to_owned(),
                        kernel: Kernel::Exact,
                        points: c.points,
                        tets: c.tets,
                        face_groups: c.face_groups,
                    });
                    return Ok(());
                }
                exact_failure = Some(
                    Record::new(DefectCode::MeshRecoveryFailed, Stage::VolumeMesh)
                        .part(name)
                        .measured("kernel", "exact")
                        .hint("the exact mesh failed the volume gate"),
                );
            }
            Err(record) => {
                outcome.attempts.push(attempt(Kernel::Exact, false));
                exact_failure = Some((*record).part(name));
            }
        }
    }
    // An exact failure is settled when the part is: a warning if a later
    // kernel meshed it, an error if none did.
    let settle = |log: &mut RecordLog, rec: Option<Record>, meshed: bool| {
        if let Some(rec) = rec {
            log.push(if meshed {
                rec.severity(Severity::Warning)
                    .action(Action::Fallback, Outcome::Applied)
            } else {
                rec.severity(Severity::Error)
                    .action(Action::Fallback, Outcome::Failed)
            });
        }
    };

    // 1. The repaired surface as it is, conforming Delaunay, for a surface
    // that passed the gate; and 2. the surface remeshed to the size field,
    // when the first fails or its worst cell is beyond the quality limits
    // (`worst_cell`). A remeshed surface has well-shaped, nearly Delaunay
    // triangles, which the constrained kernel recovers far more easily than
    // many inputs, and the cells on it are far better than on long, thin
    // input triangles (Milestone 15: on the Thingi10k subset, wall-adjacent
    // non-orthogonality 70° against 82°). But it follows the input only to
    // the element size (a curved part loses some volume), so a clean input
    // keeps its own triangulation. Of two meshes, the one with the better
    // worst cell goes to the volume gate first; the remeshed one must be
    // strictly better. `--no-remesh` keeps the input's triangulation.
    let mut conformal_failure: Option<Record> = None;
    if r.summary.gate_passed && !config.stuffing_only {
        let mut candidates: Vec<(Candidate, Option<Remeshing>)> = Vec::new();
        let mut exact_worst = None;
        match conformal(model, r, h, config, false) {
            Ok(c) => {
                exact_worst = Some(worst_cell(&c));
                candidates.push((c, None));
            }
            Err(record) => {
                outcome.attempts.push(attempt(Kernel::Delaunay, false));
                conformal_failure = Some((*record).part(name));
            }
        }
        let mut remesh_failure: Option<String> = None;
        let mut remeshed_worst = None;
        if !config.no_remesh && exact_worst.is_none_or(|w| w < 1.0) {
            match remeshed(model, r, h, config) {
                Ok((remeshing, c)) => {
                    let w = worst_cell(&c);
                    remeshed_worst = Some(w);
                    if exact_worst.is_some_and(|e| w > e) {
                        candidates.insert(0, (c, Some(remeshing)));
                    } else {
                        candidates.push((c, Some(remeshing)));
                    }
                }
                Err(reason) => {
                    outcome.attempts.push(Attempt {
                        remeshed: true,
                        ..attempt(Kernel::Delaunay, false)
                    });
                    remesh_failure = Some(reason);
                }
            }
        }
        // 3. Far beyond the limits: stuffing too, if its cells are better
        // and it keeps the volume. Its boundary only approximates the
        // surface, so it is kept only for a mesh no conformal candidate
        // comes near the limits on. Gated into
        // a scratch log: its records count only if it is kept.
        let best = candidates.first().map(|(c, _)| worst_cell(c));
        if let Some(conformal_worst) = best.filter(|&w| w < FAR_WORST_CELL && !config.strict) {
            let mut trial = RecordLog::new();
            if let Ok(s) = stuffing_candidate(model, r, h, config, &winding, &bounds)? {
                let stuffing_worst = worst_cell(&s);
                if stuffing_worst > conformal_worst {
                    let (accepted, summary) =
                        volume_gate(name, &s, r, &winding, diagonal, h, config, &mut trial)?;
                    let volume_error = if summary.surface_volume > 0.0 {
                        (summary.volume / summary.surface_volume - 1.0).abs()
                    } else {
                        f64::INFINITY
                    };
                    if accepted && volume_error <= FAR_MAX_VOLUME_ERROR {
                        // The conformal meshes are not kept.
                        for (_, remeshing) in &candidates {
                            outcome.attempts.push(Attempt {
                                remeshed: remeshing.is_some(),
                                ..attempt(Kernel::Delaunay, false)
                            });
                        }
                        outcome.attempts.push(attempt(Kernel::Stuffing, true));
                        outcome.summaries.retain(|(i, _)| *i != r.part);
                        outcome.summaries.push((r.part, summary));
                        for rec in trial.records() {
                            log.push(rec.clone());
                        }
                        log.worked_around(name, Action::Fallback);
                        if let Some(rec) = conformal_failure.take() {
                            log.push(
                                rec.severity(Severity::Warning)
                                    .action(Action::Fallback, Outcome::Applied),
                            );
                        }
                        // What remeshing gave, though not kept.
                        let remeshed = candidates.iter().find_map(|(_, m)| m.clone());
                        let mut rec = Record::new(DefectCode::VolLowQuality, Stage::VolumeMesh);
                        if let Some(m) = remeshed {
                            rec = m
                                .measure(rec)
                                .measured("remeshed_worst_cell", remeshed_worst.unwrap_or(0.0));
                        } else if let Some(reason) = remesh_failure.take() {
                            rec = rec.measured("remesh_retry", reason);
                        }
                        log.push(
                            rec.part(name)
                                .severity(Severity::Warning)
                                .measured("conformal_worst_cell", conformal_worst)
                                .measured("stuffing_worst_cell", stuffing_worst)
                                .measured("volume_error", volume_error)
                                .threshold("worst_cell", FAR_WORST_CELL)
                                .threshold("volume_error", FAR_MAX_VOLUME_ERROR)
                                .action(Action::Fallback, Outcome::Applied)
                                .hint("every conformal mesh had cells far beyond the quality limits; the stuffing fallback gave better ones, with a boundary that only approximates the surface (see the mesh's hausdorff)"),
                        );
                        settle(log, exact_failure, true);
                        outcome.meshed.push(MeshedPart {
                            part: r.part,
                            name: name.to_owned(),
                            kernel: Kernel::Stuffing,
                            points: s.points,
                            tets: s.tets,
                            face_groups: s.face_groups,
                        });
                        return Ok(());
                    }
                }
                // Tried, not kept.
                outcome.attempts.push(attempt(Kernel::Stuffing, false));
            }
        }
        for (i, (c, remeshing)) in candidates.into_iter().enumerate() {
            let (accepted, summary) = volume_gate(name, &c, r, &winding, diagonal, h, config, log)?;
            outcome.attempts.push(Attempt {
                remeshed: remeshing.is_some(),
                ..attempt(Kernel::Delaunay, accepted)
            });
            outcome.summaries.retain(|(i, _)| *i != r.part);
            outcome.summaries.push((r.part, summary));
            if accepted {
                if i > 0 {
                    // The first candidate's gate rejections are worked
                    // around by this one.
                    log.worked_around(
                        name,
                        if remeshing.is_some() {
                            Action::Remesh
                        } else {
                            Action::Fallback
                        },
                    );
                }
                if let Some(remeshing) = remeshing {
                    // Findings deferred to remeshing are now dealt with.
                    log.deferred_applied(name, Action::Remesh);
                    let rec = match conformal_failure.take() {
                        // The repaired surface did not mesh: remeshing
                        // worked around it.
                        Some(rec) => remeshing
                            .measure(rec.severity(Severity::Warning))
                            .action(Action::Remesh, Outcome::Applied),
                        // It meshed, with poor cells: the boundary
                        // triangulation is no longer the input's, so say so.
                        None => remeshing
                            .measure(
                                Record::new(DefectCode::VolLowQuality, Stage::VolumeMesh)
                                    .part(name)
                                    .severity(Severity::Info)
                                    .measured("input_worst_cell", exact_worst.unwrap_or(0.0))
                                    .measured("remeshed_worst_cell", remeshed_worst.unwrap_or(0.0)),
                            )
                            .threshold("worst_cell", 1.0)
                                .action(Action::Remesh, Outcome::Applied)
                            .hint("the input's triangulation gave cells beyond the quality limits; the surface was remeshed to the element size, which gave better ones"),
                    };
                    log.push(rec);
                }
                settle(log, exact_failure, true);
                outcome.meshed.push(MeshedPart {
                    part: r.part,
                    name: name.to_owned(),
                    kernel: Kernel::Delaunay,
                    points: c.points,
                    tets: c.tets,
                    face_groups: c.face_groups,
                });
                return Ok(());
            }
            if remeshing.is_some() {
                remesh_failure = Some("rejected by the volume gate".to_owned());
            } else {
                conformal_failure = Some(
                    Record::new(DefectCode::MeshRecoveryFailed, Stage::VolumeMesh)
                        .part(name)
                        .hint("the conformal mesh failed the volume gate"),
                );
            }
        }
        if let Some(reason) = remesh_failure {
            conformal_failure = conformal_failure.map(|rec| rec.measured("remesh_retry", reason));
        }
    }

    if config.strict {
        // Stuffing is not allowed: the wrap is the last resort.
        if wrapped_and_meshed(model, r, config, log, outcome)? {
            return Ok(());
        }
        if let Some(rec) = exact_failure {
            log.push(rec.severity(Severity::Error));
        }
        if let Some(rec) = conformal_failure {
            log.push(rec.severity(Severity::Error));
        }
        log.push(
            Record::new(DefectCode::MeshStrictUnmet, Stage::VolumeMesh)
                .part(name)
                .hint(if r.summary.gate_passed {
                    "the conformal kernel failed and --strict forbids the stuffing fallback"
                } else {
                    "the repaired surface failed the surface gate, so only the stuffing fallback could mesh it"
                }),
        );
        return Ok(());
    }

    // 4. Fallback: isosurface stuffing.
    let stuffed = match stuffing_candidate(model, r, h, config, &winding, &bounds)? {
        Ok(c) => Some(c),
        Err(rec) => {
            log.push(rec.part(name));
            outcome.attempts.push(attempt(Kernel::Stuffing, false));
            None
        }
    };
    let mut fallback_ok = false;
    if let Some(c) = stuffed {
        let (accepted, summary) = volume_gate(name, &c, r, &winding, diagonal, h, config, log)?;
        outcome.attempts.push(attempt(Kernel::Stuffing, accepted));
        outcome.summaries.retain(|(i, _)| *i != r.part);
        outcome.summaries.push((r.part, summary));
        if accepted {
            fallback_ok = true;
            // What sent the part here is now worked around.
            log.worked_around(name, Action::Fallback);
            outcome.meshed.push(MeshedPart {
                part: r.part,
                name: name.to_owned(),
                kernel: Kernel::Stuffing,
                points: c.points,
                tets: c.tets,
                face_groups: c.face_groups,
            });
        }
    }
    // 5. The last resort: the wrap.
    if !fallback_ok && wrapped_and_meshed(model, r, config, log, outcome)? {
        fallback_ok = true;
    }
    settle(log, exact_failure, fallback_ok);
    if let Some(rec) = conformal_failure {
        log.push(if fallback_ok {
            rec.severity(Severity::Warning)
                .action(Action::Fallback, Outcome::Applied)
        } else {
            rec.severity(Severity::Error)
                .action(Action::Fallback, Outcome::Failed)
        });
    }
    Ok(())
}

/// Repair rung 3 for a part that found no mesh (its surface failed the
/// gate, or passed it and the conformal kernel and stuffing both failed):
/// wrap it and mesh the wrap through the whole pipeline. Whether the wrap
/// meshed; its errors that sent the part here are then worked around.
/// Under `--strict` a surface that passed the gate is not wrapped: it is
/// meshed conformally as it is, or the part fails.
fn wrapped_and_meshed(
    model: &SurfaceModel,
    r: &RepairedPart,
    config: &JobConfig,
    log: &mut RecordLog,
    outcome: &mut MeshOutcome,
) -> Result<bool, InternalError> {
    let _prof = mesher_tet::prof::scope("volume::wrapped_and_meshed");
    if r.summary.rung >= 3 || (r.summary.gate_passed && config.strict) {
        return Ok(false);
    }
    let Some(w) = wrap_part(model, r, &Tolerances::for_job(config), log)? else {
        return Ok(false);
    };
    mesh_part(model, &w, config, log, outcome)?;
    outcome.wrapped.retain(|x| x.part != r.part);
    outcome.wrapped.push(w);
    let meshed = outcome.meshed.iter().any(|m| m.part == r.part);
    if meshed {
        log.worked_around(&r.name, Action::Wrap);
    }
    Ok(meshed)
}

/// The part meshed by isosurface stuffing of its winding number, or the
/// budget record (without its part) if the lattice is too large. Where the
/// size field asks for less than the element size, the lattice is graded:
/// uniformly fine at the surface (h halved until it reaches the field's
/// smallest size, at most twice), coarsening to h inside, which keeps the
/// proven angle bounds.
fn stuffing_candidate(
    model: &SurfaceModel,
    r: &RepairedPart,
    h: f64,
    config: &JobConfig,
    winding: &FastWindingNumber<'_>,
    bounds: &Aabb,
) -> Result<Result<Candidate, Record>, InternalError> {
    let _prof = mesher_tet::prof::scope("volume::stuffing_candidate");
    let levels = match size_field(model, r, h, config) {
        Some(field) if !field.is_uniform() => {
            libm::log2(h / field.min()).ceil().clamp(0.0, 2.0) as u32
        }
        _ => 0,
    };
    let leaf = h / f64::from(1u32 << levels);
    let level_fn = |p: Point3| 0.5 - winding.winding_number(p);
    let stuffed = if levels == 0 {
        stuff(level_fn, bounds, &StuffingParams::with_cell_size(h))
    } else {
        // Every surface vertex and face centroid seeds the search for the
        // cubes the surface passes through.
        let seeds: Vec<Point3> = r
            .points
            .iter()
            .copied()
            .chain(r.faces.iter().map(|t| {
                let [a, b, c] = t.map(|v| r.points[v as usize]);
                [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0)
            }))
            .collect();
        stuff_graded(
            level_fn,
            bounds,
            &seeds,
            &StuffingParams::with_cell_size(leaf),
            levels,
        )
    };
    match stuffed {
        Ok(s) => Ok(Ok(Candidate {
            kernel: Kernel::Stuffing,
            face_groups: stuffed_groups(model, r, &s.points, &s.tets),
            points: s.points,
            tets: s.tets,
            dropped: s.dropped,
            steiner: 0,
            min_size: (levels > 0).then_some(leaf),
            remeshed: false,
        })),
        Err(StuffingError::Budget {
            lattice_vertices,
            limit,
        }) => Ok(Err(Record::new(
            DefectCode::MeshBudgetExceeded,
            Stage::VolumeMesh,
        )
        .measured("lattice_vertices", lattice_vertices)
        .threshold("max_lattice_vertices", limit)
        .hint("raise --element-size"))),
        Err(e) => Err(InternalError(format!("part {}: {e}", r.name))),
    }
}

/// The exact constrained mesh of a part's repaired surface, rounded to f64,
/// with a `MESH_ROUNDING_REPAIRED` record for each tet rounding turned over
/// or flattened and untangling fixed; or the record explaining why not.
/// No refinement: the tets join the surface's vertices and the Steiner
/// points recovery needed (refinement is Milestone 14).
fn exact(
    model: &SurfaceModel,
    r: &RepairedPart,
    name: &str,
) -> Result<(Candidate, Vec<Record>), Box<Record>> {
    let _prof = mesher_tet::prof::scope("volume::exact");
    let recovery = |e: &dyn std::fmt::Display| {
        Record::new(DefectCode::MeshRecoveryFailed, Stage::VolumeMesh)
            .measured("kernel", "exact")
            .hint(format!("exact meshing failed: {e}"))
    };
    let region = mesh_region(&r.points, &r.faces, 0, 5).map_err(|e| {
        Box::new(match &e {
            RegionError::Faces(FaceError::Segments(SegmentError::Budget { points, limit })) => {
                Record::new(DefectCode::MeshBudgetExceeded, Stage::VolumeMesh)
                    .measured("kernel", "exact")
                    .measured("points", *points as u64)
                    .threshold("max_points", *limit as u64)
            }
            RegionError::Missing { triangles } => {
                recovery(&e).measured("missing_faces", triangles.len() as u64)
            }
            _ => recovery(&e),
        })
    })?;
    let steiner = (region.recovered.steiner + region.recovered.steiner_on_faces) as u64;
    let m = round_region(region, &r.points, &r.faces).map_err(|e| {
        let rec = Record::new(DefectCode::MeshRoundingFailed, Stage::VolumeMesh)
            .hint(format!("{e}; the next kernel takes over"));
        Box::new(match &e {
            RoundingError::Inverted {
                tets,
                inverted_before,
                untangles,
            } => {
                let first = tets
                    .first()
                    .map(|t| [0, 1, 2].map(|k| (t[0][k] + t[1][k] + t[2][k] + t[3][k]) / 4.0));
                let rec = rec
                    .measured("tets_left", tets.len() as u64)
                    .measured("turned_over_before", *inverted_before as u64)
                    .measured("fixes", *untangles as u64);
                match first {
                    Some(p) => rec.point(p),
                    None => rec,
                }
            }
            _ => rec,
        })
    })?;
    // Where each fix was, in the rounded mesh.
    let local: BTreeMap<u32, u32> = m
        .origin
        .iter()
        .enumerate()
        .map(|(i, &o)| (o, i as u32))
        .collect();
    let centroid = |t: &[u32; 4]| -> Option<Point3> {
        let p: Vec<Point3> = t
            .iter()
            .map(|v| local.get(v).map(|&i| m.points[i as usize]))
            .collect::<Option<_>>()?;
        Some([0, 1, 2].map(|k| (p[0][k] + p[1][k] + p[2][k] + p[3][k]) / 4.0))
    };
    let records = m
        .untangles
        .iter()
        .map(|u| {
            let rec = Record::new(DefectCode::MeshRoundingRepaired, Stage::VolumeMesh)
                .part(name)
                .action(Action::Untangle, Outcome::Applied);
            match u {
                Untangle::Move { from, to, .. } => rec
                    .measured("fix", "vertex_move")
                    .measured(
                        "moved",
                        libm::sqrt(
                            (0..3)
                                .map(|k| (to[k] - from[k]) * (to[k] - from[k]))
                                .sum::<f64>(),
                        ),
                    )
                    .point(*to),
                Untangle::Flip23 { tet, .. } => {
                    let rec = rec.measured("fix", "flip_2_3");
                    match centroid(tet) {
                        Some(p) => rec.point(p),
                        None => rec,
                    }
                }
                Untangle::RemoveEdge { tet, .. } => {
                    let rec = rec.measured("fix", "edge_removal");
                    match centroid(tet) {
                        Some(p) => rec.point(p),
                        None => rec,
                    }
                }
            }
        })
        .collect();
    let group_of = |f: usize| -> Option<u32> {
        r.face_origin
            .get(f)
            .copied()
            .flatten()
            .and_then(|t| model.polygon_of(t))
            .and_then(|p| p.group)
    };
    let mut groups: BTreeMap<String, Vec<[u32; 3]>> = BTreeMap::new();
    for &(face, o) in &m.boundary {
        if let Some(g) = group_of(o as usize).and_then(|g| model.groups.get(g as usize)) {
            groups.entry(g.clone()).or_default().push(face);
        }
    }
    let mut c = Candidate {
        kernel: Kernel::Exact,
        points: m.points,
        tets: m.tets,
        dropped: 0,
        steiner,
        face_groups: groups.into_iter().collect(),
        min_size: None,
        remeshed: false,
    };
    improve_cells(&mut c);
    Ok((c, records))
}

/// Improve a conformal candidate's poor cells by local operations
/// (`mesher_tet::improve`). Boundary faces and their vertices never change,
/// so the face groups only need renumbering; if one could not be, the
/// candidate is left as it was.
/// A candidate's worst cell against the volume gate's quality limits:
/// the smallest of its smallest dihedral angle over `SLIVER_DEGREES`, (90°
/// − its worst non-orthogonality) over (90° −
/// `MAX_NON_ORTHOGONALITY_DEGREES`), and `MAX_SKEWNESS` over its worst
/// skewness. So 1 is at the limits and smaller is worse; 0 for a candidate
/// with no cells.
fn worst_cell(c: &Candidate) -> f64 {
    if c.tets.is_empty() {
        return 0.0;
    }
    let q = fv_quality(&c.points, &c.tets);
    let dihedral = q.min_dihedral_degrees / SLIVER_DEGREES;
    let non_orthogonality =
        (90.0 - q.max_non_orthogonality_degrees) / (90.0 - MAX_NON_ORTHOGONALITY_DEGREES);
    let skewness = if q.max_skewness > 0.0 {
        MAX_SKEWNESS / q.max_skewness
    } else {
        f64::INFINITY
    };
    let w = dihedral.min(non_orthogonality).min(skewness);
    if w.is_nan() { 0.0 } else { w }
}

/// Grade a candidate's cells to `MAX_SIZE_RATIO` (`mesher_tet::grade`),
/// its side sets following split boundary faces; whether it split any.
fn grade_cells(c: &mut Candidate) -> bool {
    let _prof = mesher_tet::prof::scope("volume::grade_cells");
    let tagged: Vec<([u32; 3], u32)> = c
        .face_groups
        .iter()
        .enumerate()
        .flat_map(|(g, (_, faces))| faces.iter().map(move |&f| (f, g as u32)))
        .collect();
    let g = grade(
        &c.points,
        &c.tets,
        &tagged,
        &GradeParams {
            limit: MAX_SIZE_RATIO,
            ..GradeParams::default()
        },
    );
    if g.stats.splits == 0 {
        return false;
    }
    for (_, faces) in &mut c.face_groups {
        faces.clear();
    }
    for (f, tag) in g.faces {
        if let Some((_, faces)) = c.face_groups.get_mut(tag as usize) {
            faces.push(f);
        }
    }
    c.points = g.points;
    c.tets = g.tets;
    true
}

/// Improve a candidate's cells (`mesher_tet::improve`), never making a jump
/// in sizes beyond `MAX_SIZE_RATIO`; grade the jumps the kernel left
/// (`grade_cells`), and if that split any cells, improve again, still
/// without making one. Grading after improvement, not before: improvement
/// works on the kernel's cells, and the second pass repairs the quality
/// bisection costs around the graded cells.
fn improve_cells(c: &mut Candidate) {
    let _prof = mesher_tet::prof::scope("volume::improve_cells");
    improve_once(c);
    if grade_cells(c) {
        improve_once(c);
    }
}

fn improve_once(c: &mut Candidate) {
    let params = ImproveParams {
        max_size_ratio: MAX_SIZE_RATIO,
        ..ImproveParams::default()
    };
    let r = improve(&c.points, &c.tets, &params);
    let renumbered: Option<Vec<(String, Vec<[u32; 3]>)>> = c
        .face_groups
        .iter()
        .map(|(name, faces)| {
            let faces = faces
                .iter()
                .map(|f| {
                    let g = f.map(|v| r.map.get(v as usize).copied().unwrap_or(u32::MAX));
                    (!g.contains(&u32::MAX)).then_some(g)
                })
                .collect::<Option<Vec<_>>>()?;
            Some((name.clone(), faces))
        })
        .collect();
    if let Some(groups) = renumbered {
        c.face_groups = groups;
        c.points = r.points;
        c.tets = r.tets;
    }
}

/// The conforming Delaunay mesh of a repaired part, or the record explaining
/// why not.
/// The repaired surface remeshed to the size field, and its conformal mesh:
/// the number of remeshed triangles and the candidate, or why not.
/// What remeshing did to a part's surface, for its records.
#[derive(Debug, Clone)]
struct Remeshing {
    triangles: usize,
    /// Features too small to mesh that were merged away (each moves the
    /// surface by at most a tenth of the local size).
    merged_small_features: usize,
    /// The ladder rung whose surface this is (`remeshed`).
    rung: String,
    /// The other remeshings tried: each rung's name and why it was not
    /// kept (it did not mesh, or its worst cell was worse).
    retried: Vec<String>,
}

impl Remeshing {
    fn measure(self, rec: Record) -> Record {
        let mut rec = rec
            .measured("remeshed_triangles", self.triangles as u64)
            .measured("remesh_rung", self.rung);
        if self.merged_small_features > 0 {
            rec = rec.measured("merged_small_features", self.merged_small_features as u64);
        }
        if !self.retried.is_empty() {
            rec = rec.measured("other_remeshings", self.retried);
        }
        rec
    }
}

fn remeshed(
    model: &SurfaceModel,
    r: &RepairedPart,
    h: f64,
    config: &JobConfig,
) -> Result<(Remeshing, Candidate), String> {
    let _prof = mesher_tet::prof::scope("volume::remeshed");
    let field = size_field(model, r, h, config).unwrap_or_else(|| SizeField::uniform(h));
    // Group boundaries (side sets) stay feature lines.
    let groups = crate::features::face_groups(model, r);
    let base = crate::features::remesh_params(r, config);
    // Snapping and free small collapses merge the narrow strips along
    // feature edges, which otherwise leave needles. A remeshed surface that
    // does not mesh (most often one that crosses itself) is remeshed again
    // with every operation checked against crossing, then with each of
    // those merges turned off in turn, back to the remesher as it was. Each
    // rung only runs when the one before did not mesh, so a model that
    // meshed on an earlier remesher still does. The checks refuse moves
    // that later passes would undo, so they are never the first attempt:
    // on the subsets they lost as many models as they won.
    let ladder: [(&str, bool, bool, bool); 5] = [
        // (name, snap strips, free small collapses, check crossings)
        ("snapped", true, true, false),
        ("snapped, checked", true, true, true),
        ("free collapses", false, true, false),
        ("free collapses, checked", false, true, true),
        ("plain", false, false, false),
    ];
    // A surface that meshes with cells far beyond the limits
    // (`FAR_WORST_CELL`) does not end the ladder: on some models snapping or
    // free collapses leave worse cells than they remove (72329, 1356636: the
    // plain remesher's cells are within the limits). The rungs below are
    // tried until one meshes within them, and the best kept. A checked rung
    // only runs when the same rung unchecked did not mesh.
    let mut tried: Vec<String> = Vec::new();
    let mut best: Option<(f64, usize, Remeshing, Candidate)> = None;
    let mut unchecked_meshed = false;
    for (name, snap_strips, free_small_collapses, check_crossings) in ladder {
        if check_crossings && unchecked_meshed {
            continue;
        }
        let params = RemeshParams {
            snap_strips,
            free_small_collapses,
            check_crossings,
            ..base.clone()
        };
        match remeshed_with(model, r, h, config, &field, &groups, &params) {
            Ok((remeshing, c)) => {
                unchecked_meshed = !check_crossings;
                let w = worst_cell(&c);
                tried.push(format!("{name}: worst cell {w:.3}"));
                if best.as_ref().is_none_or(|b| w > b.0) {
                    best = Some((w, tried.len() - 1, remeshing, c));
                }
                if w >= FAR_WORST_CELL {
                    break;
                }
            }
            Err(reason) => {
                unchecked_meshed = false;
                tried.push(format!("{name}: {reason}"));
            }
        }
    }
    match best {
        Some((_, kept, remeshing, c)) => {
            let rung = ladder_rung(&tried[kept]);
            tried.remove(kept);
            Ok((
                Remeshing {
                    rung,
                    retried: tried,
                    ..remeshing
                },
                c,
            ))
        }
        None => Err(tried.join("; ")),
    }
}

/// The rung's name in one of `remeshed`'s entries, `name: why`.
fn ladder_rung(entry: &str) -> String {
    entry.split(": ").next().unwrap_or(entry).to_owned()
}

/// `remeshed` with one rung's remesher.
fn remeshed_with(
    model: &SurfaceModel,
    r: &RepairedPart,
    h: f64,
    config: &JobConfig,
    field: &SizeField,
    groups: &[u32],
    params: &RemeshParams,
) -> Result<(Remeshing, Candidate), String> {
    let prof_remesh = mesher_tet::prof::scope("volume::remeshed/remesh");
    let s =
        remesh(&r.points, &r.faces, groups, &|p| field.at(p), params).map_err(|e| e.to_string())?;
    let remeshing = Remeshing {
        triangles: s.triangles.len(),
        merged_small_features: s.merged_small_features,
        rung: String::new(),
        retried: Vec::new(),
    };
    let surface = RepairedPart {
        points: s.points,
        faces: s.triangles,
        vertex_origin: Vec::new(),
        face_origin: s
            .origin
            .iter()
            .map(|&o| r.face_origin.get(o as usize).copied().flatten())
            .collect(),
        ..r.clone()
    };
    drop(prof_remesh);
    let prof_gate = mesher_tet::prof::scope("volume::remeshed/gate");
    let failures = gate_failures(&surface.points, &surface.faces, &Tolerances::default());
    drop(prof_gate);
    if !failures.is_empty() {
        return Err(format!(
            "the remeshed surface failed the surface gate: {}",
            failures.join(", ")
        ));
    }
    let mut c = conformal(model, &surface, h, config, true).map_err(|rec| {
        rec.hint
            .clone()
            .unwrap_or_else(|| "conformal meshing failed".to_owned())
    })?;
    c.remeshed = true;
    Ok((remeshing, c))
}

/// The conformal mesh of a surface. With `own_facets`, every triangle is
/// its own facet (a remeshed surface: its triangles are nearly Delaunay, and
/// merging those that lie on one input facet makes recovery cascade).
fn conformal(
    model: &SurfaceModel,
    r: &RepairedPart,
    h: f64,
    config: &JobConfig,
    own_facets: bool,
) -> Result<Candidate, Box<Record>> {
    let _prof = mesher_tet::prof::scope("volume::conformal");
    // OBJ group of each input face (fill faces have none).
    let group_of = |f: usize| -> Option<u32> {
        r.face_origin
            .get(f)
            .copied()
            .flatten()
            .and_then(|t| model.polygon_of(t))
            .and_then(|p| p.group)
    };
    // The sizing field grades the size down where the surface curves or the
    // part is thin.
    let prof_size = mesher_tet::prof::scope("volume::conformal/size_field");
    let size = size_field(model, r, h, config).filter(|f| !f.is_uniform());
    drop(prof_size);
    let min_size = size.as_ref().map(SizeField::min);
    let side_sets: Vec<u32> = (0..r.faces.len())
        .map(|f| group_of(f).unwrap_or(u32::MAX))
        .collect();
    let params = ConformingParams {
        groups: if own_facets {
            (0..r.faces.len() as u32).collect()
        } else {
            side_sets.clone()
        },
        side_sets,
        size: size.map(std::sync::Arc::new),
        ..ConformingParams::with_size(h)
    };
    let m = conforming_mesh(&r.points, &r.faces, &params).map_err(|e| {
        let rec = match &e {
            ConformingError::Budget { points, limit } => {
                Record::new(DefectCode::MeshBudgetExceeded, Stage::VolumeMesh)
                    .measured("points", *points as u64)
                    .threshold("max_points", *limit as u64)
            }
            ConformingError::Work { steps, limit } => {
                Record::new(DefectCode::MeshBudgetExceeded, Stage::VolumeMesh)
                    .measured("work_steps", *steps)
                    .threshold("max_work_steps", *limit)
            }
            ConformingError::Recovery { missing_faces } => {
                Record::new(DefectCode::MeshRecoveryFailed, Stage::VolumeMesh)
                    .measured("missing_faces", *missing_faces as u64)
            }
            ConformingError::Stalled {
                cause,
                value,
                limit,
            } => {
                let (measure, max) = match cause {
                    Stall::Points => ("new_points", "max_new_points"),
                    Stall::Missing | Stall::Cascade => ("missing_pieces", "max_missing_pieces"),
                };
                Record::new(DefectCode::MeshRecoveryFailed, Stage::VolumeMesh)
                    .measured("stalled", cause.as_str())
                    .measured(measure, *value as u64)
                    .threshold(max, *limit as u64)
            }
            ConformingError::Surface(_) | ConformingError::Delaunay(_) => {
                Record::new(DefectCode::MeshRecoveryFailed, Stage::VolumeMesh)
            }
        };
        Box::new(rec.hint(format!("conformal meshing failed: {e}")))
    })?;
    let mut groups: BTreeMap<String, Vec<[u32; 3]>> = BTreeMap::new();
    for &(face, o) in &m.boundary {
        if let Some(name) = group_of(o as usize).and_then(|g| model.groups.get(g as usize)) {
            groups.entry(name.clone()).or_default().push(face);
        }
    }
    let mut c = Candidate {
        kernel: Kernel::Delaunay,
        points: m.points,
        tets: m.tets,
        dropped: 0,
        steiner: (m.boundary_steiner + m.interior_steiner) as u64,
        face_groups: groups.into_iter().collect(),
        min_size,
        remeshed: false,
    };
    improve_cells(&mut c);
    Ok(c)
}

/// Faces of a tet mesh held by more than two tets (`crowded`), and faces held
/// by two tets on the same side of them (`same_side`). Both are zero for a
/// mesh that is a manifold with boundary.
fn face_manifold(tets: &[[u32; 4]]) -> (u64, u64) {
    // Each face with its vertices sorted, and whether sorting them was an
    // even permutation of the face as the tet lists it (all four listed the
    // same way round, so two tets on opposite sides list it with opposite
    // parities).
    let mut faces: Vec<([u32; 3], bool)> = Vec::with_capacity(4 * tets.len());
    for t in tets {
        for f in [
            [t[1], t[2], t[3]],
            [t[0], t[3], t[2]],
            [t[0], t[1], t[3]],
            [t[0], t[2], t[1]],
        ] {
            let inversions =
                u32::from(f[0] > f[1]) + u32::from(f[0] > f[2]) + u32::from(f[1] > f[2]);
            let mut key = f;
            key.sort_unstable();
            faces.push((key, inversions % 2 == 0));
        }
    }
    faces.sort_unstable();
    let (mut crowded, mut same_side) = (0, 0);
    for group in faces.chunk_by(|a, b| a.0 == b.0) {
        match group {
            [_] => {}
            [a, b] if a.1 == b.1 => same_side += 1,
            [_, _] => {}
            _ => crowded += 1,
        }
    }
    (crowded, same_side)
}

/// Check a candidate mesh; returns whether it is accepted and its summary.
#[allow(clippy::too_many_arguments)]
fn volume_gate(
    name: &str,
    c: &Candidate,
    r: &RepairedPart,
    winding: &FastWindingNumber<'_>,
    diagonal: f64,
    h: f64,
    config: &JobConfig,
    log: &mut RecordLog,
) -> Result<(bool, MeshSummary), InternalError> {
    let _prof = mesher_tet::prof::scope("volume::volume_gate");
    let gate = |code| Record::new(code, Stage::VolumeGate).part(name);
    let mut accepted = true;
    let (points, tets) = (&c.points, &c.tets);
    let not_positive = tets
        .iter()
        .filter(|t| {
            let [a, b, d, e] = t.map(|v| points[v as usize]);
            orient3d(a, b, d, e) != Sign::Positive
        })
        .count() as u64
        + c.dropped;
    if not_positive > 0 {
        accepted = false;
        log.push(
            gate(DefectCode::VolInvertedTet)
                .measured("tets", not_positive)
                .measured("kernel", format!("{:?}", c.kernel))
                .hint("the volume kernels should never do this; please report it"),
        );
    }
    // With every tet positive and the boundary on the input surface, a mesh
    // whose faces each have at most one tet on each side covers every point
    // of the part once: tets overlapping one another would show here.
    let (crowded, same_side) = face_manifold(tets);
    if crowded + same_side > 0 {
        accepted = false;
        log.push(
            gate(DefectCode::VolNonmanifold)
                .measured("crowded_faces", crowded)
                .measured("same_side_faces", same_side)
                .measured("kernel", format!("{:?}", c.kernel))
                .hint("the volume kernels should never do this; please report it"),
        );
    }
    let q = fv_quality(points, tets);
    // Positive by the exact predicate but flat to rounding: a solver
    // computing the volume in floating point gets zero or less.
    if q.flat_in_float > 0 {
        accepted = false;
        log.push(
            gate(DefectCode::VolInvertedTet)
                .measured("flat_in_float", q.flat_in_float as u64)
                .measured("kernel", format!("{:?}", c.kernel))
                .hint("tets that are positive exactly but have no positive volume in floating point; a finite-volume solver cannot use them"),
        );
    }
    let (lo, hi) = match c.kernel {
        Kernel::Stuffing => (DIHEDRAL_MIN_DEGREES, DIHEDRAL_MAX_DEGREES),
        Kernel::Delaunay | Kernel::Exact => (SLIVER_DEGREES, 180.0 - SLIVER_DEGREES),
    };
    // The worst cell of each measure beyond its limit, located: the cell's
    // 1-based element number in the part's block, and its centroid (the
    // face's, for a face measure).
    let located = |rec: Record, w: Option<WorstAt>| match w {
        Some(w) => {
            let rec = rec.measured("element", w.tet as u64 + 1).point(w.point);
            match w.neighbour {
                Some(n) => rec.measured("neighbour", n as u64 + 1),
                None => rec,
            }
        }
        None => rec,
    };
    if !tets.is_empty() && (q.min_dihedral_degrees < lo || q.max_dihedral_degrees > hi) {
        log.push(located(
            gate(DefectCode::VolLowQuality)
                .measured("min_dihedral_degrees", q.min_dihedral_degrees)
                .measured("max_dihedral_degrees", q.max_dihedral_degrees)
                .threshold("min_dihedral_degrees", lo)
                .threshold("max_dihedral_degrees", hi)
                .hint("flat cells remain (located: the one with the smallest dihedral angle); a finite-volume solver will lose some accuracy there"),
            q.worst_dihedral,
        ));
    }
    if q.max_non_orthogonality_degrees > MAX_NON_ORTHOGONALITY_DEGREES {
        log.push(located(
            gate(DefectCode::VolLowQuality)
                .measured("max_non_orthogonality_degrees", q.max_non_orthogonality_degrees)
                .threshold("max_non_orthogonality_degrees", MAX_NON_ORTHOGONALITY_DEGREES)
                .hint("faces whose normal is far from the line joining their cells' centroids (located: the worst); a solver needs strong non-orthogonal correction there"),
            q.worst_non_orthogonality,
        ));
    }
    if q.max_skewness > MAX_SKEWNESS {
        log.push(located(
            gate(DefectCode::VolLowQuality)
                .measured("max_skewness", q.max_skewness)
                .threshold("max_skewness", MAX_SKEWNESS)
                .hint("faces whose centroid is far from where the line joining their cells' centroids crosses them (located: the worst); interpolation to the face loses accuracy there"),
            q.worst_skewness,
        ));
    }
    if q.max_size_ratio > MAX_SIZE_RATIO {
        log.push(located(
            gate(DefectCode::VolGradation)
                .measured("max_size_ratio", q.max_size_ratio)
                .threshold("max_size_ratio", MAX_SIZE_RATIO)
                .hint("a small cell meets a much larger one across a face (located: the worst); the size changes more abruptly than the gradation allows"),
            q.worst_size_ratio,
        ));
    }
    let single = VolumeMesh {
        points: points.clone(),
        blocks: vec![Block {
            name: name.to_owned(),
            tets: tets.clone(),
        }],
        side_sets: Vec::new(),
    };
    let boundary: Vec<[u32; 3]> = single
        .boundary_sides(0)
        .into_iter()
        .filter_map(|(e, s)| tets.get(e as usize).and_then(|&t| tet_side_nodes(t, s)))
        .collect();
    let distance = if tets.is_empty() {
        f64::INFINITY
    } else {
        hausdorff(&r.points, &r.faces, winding, points, &boundary)
            .map_err(|e| InternalError(format!("part {name}: {e}")))?
    };
    let relative = if diagonal > 0.0 {
        distance / diagonal
    } else {
        0.0
    };
    if tets.is_empty() || relative > config.max_deviation {
        accepted = false;
        let rec = gate(DefectCode::VolBoundaryDeviation)
            .measured("tets", tets.len() as u64)
            .threshold("max_deviation_relative", config.max_deviation);
        log.push(if tets.is_empty() {
            rec.hint("nothing was meshed: the part is thinner than the element size; lower --element-size")
        } else {
            rec.measured("hausdorff", distance)
                .measured("hausdorff_relative", relative)
                .hint("lower --element-size for a closer fit, or raise --max-deviation")
        });
    }
    // JSON has no infinities or NaN (they would be written as null, which no
    // reader accepts as a number), so a value that cannot be measured is -1.
    let finite = |x: f64| if x.is_finite() { x } else { -1.0 };
    let summary = MeshSummary {
        kernel: c.kernel,
        accepted,
        tets: tets.len() as u64,
        nodes: points.len() as u64,
        steiner_points: c.steiner,
        element_size: finite(h),
        min_element_size: c.min_size.map(finite),
        remeshed: c.remeshed,
        volume: finite(single.volume()),
        surface_volume: finite(enclosed_volume(&r.points, &r.faces)),
        min_dihedral_degrees: finite(q.min_dihedral_degrees),
        max_dihedral_degrees: finite(q.max_dihedral_degrees),
        dihedral_histogram: q.dihedral_histogram.clone(),
        max_non_orthogonality_degrees: finite(q.max_non_orthogonality_degrees),
        mean_non_orthogonality_degrees: finite(q.mean_non_orthogonality_degrees),
        max_skewness: finite(q.max_skewness),
        max_volume_ratio: finite(q.max_volume_ratio),
        max_size_ratio: finite(q.max_size_ratio),
        hausdorff: finite(distance),
        hausdorff_relative: finite(relative),
        worst_cells: worst_cells(&q),
    };
    Ok((accepted, summary))
}

/// The located worst values of a mesh's quality, for its summary.
fn worst_cells(q: &FvQuality) -> Vec<WorstCell> {
    [
        (QualityMeasure::MinDihedral, false, q.worst_dihedral),
        (
            QualityMeasure::NonOrthogonality,
            false,
            q.worst_non_orthogonality,
        ),
        (QualityMeasure::Skewness, false, q.worst_skewness),
        (
            QualityMeasure::NonOrthogonality,
            true,
            q.wall_worst_non_orthogonality,
        ),
        (QualityMeasure::Skewness, true, q.wall_worst_skewness),
        (QualityMeasure::SizeRatio, false, q.worst_size_ratio),
    ]
    .into_iter()
    .filter_map(|(measure, wall, w)| {
        let w = w?;
        w.value.is_finite().then(|| WorstCell {
            measure,
            wall,
            value: w.value,
            element: w.tet as u64 + 1,
            neighbour: w.neighbour.map(|n| n as u64 + 1),
            point: w.point,
        })
    })
    .collect()
}

/// The OBJ groups of a stuffed mesh's boundary: each boundary face takes the
/// group of a repaired triangle near its centroid (none for a fill face):
/// of the nearest one and those around its corners, the one facing most its
/// way, weighed against how far it is. Stuffing's boundary only approximates
/// the surface and cuts sharp edges at an angle, so faces near a group
/// boundary may go to either side.
fn stuffed_groups(
    model: &SurfaceModel,
    r: &RepairedPart,
    points: &[Point3],
    tets: &[[u32; 4]],
) -> Vec<(String, Vec<[u32; 3]>)> {
    let groups = crate::features::face_groups(model, r);
    if groups.iter().all(|&g| g == u32::MAX) {
        return Vec::new();
    }
    let Ok(soup) = TriangleSoup::new(&r.points, &r.faces) else {
        return Vec::new();
    };
    let bvh = Bvh::new(soup);
    let unit = |a: Point3, b: Point3, c: Point3| {
        let n = cross(sub(b, a), sub(c, a));
        let l = norm(n);
        if l > 0.0 { scale(n, 1.0 / l) } else { [0.0; 3] }
    };
    let normals: Vec<Point3> = r
        .faces
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| r.points[v as usize]);
            unit(a, b, c)
        })
        .collect();
    let mut around: Vec<Vec<u32>> = vec![Vec::new(); r.points.len()];
    for (t, f) in r.faces.iter().enumerate() {
        for &v in f {
            around[v as usize].push(t as u32);
        }
    }
    // Each boundary face, with the tet's fourth vertex (to orient it out).
    let mut count: BTreeMap<[u32; 3], ([u32; 3], u32, u32)> = BTreeMap::new();
    for &t in tets {
        for side in 1..=4u8 {
            if let Some(f) = tet_side_nodes(t, side) {
                let mut k = f;
                k.sort_unstable();
                let opposite = t.iter().copied().find(|v| !f.contains(v)).unwrap_or(f[0]);
                count.entry(k).or_insert((f, 0, opposite)).1 += 1;
            }
        }
    }
    let mut out: BTreeMap<String, Vec<[u32; 3]>> = BTreeMap::new();
    for (f, n, opposite) in count.into_values() {
        if n != 1 {
            continue;
        }
        let [a, b, c] = f.map(|v| points[v as usize]);
        let m = [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0);
        let mut nf = unit(a, b, c);
        if dot(nf, sub(points[opposite as usize], a)) > 0.0 {
            nf = scale(nf, -1.0);
        }
        let Some(near) = bvh.nearest(m) else {
            continue;
        };
        // Score: distance, doubled for a triangle facing across.
        let score = |t: u32| {
            let tri = r.faces[t as usize].map(|v| r.points[v as usize]);
            let d = norm(sub(m, closest_point_on_triangle(m, tri)));
            (d + 1e-300) * (2.0 - dot(nf, normals[t as usize]))
        };
        let mut best = (score(near.triangle), near.triangle);
        for &v in &r.faces[near.triangle as usize] {
            for &t in &around[v as usize] {
                let s = score(t);
                if s < best.0 || (s == best.0 && t < best.1) {
                    best = (s, t);
                }
            }
        }
        let name = groups
            .get(best.1 as usize)
            .copied()
            .and_then(|g| model.groups.get(g as usize));
        if let Some(name) = name {
            out.entry(name.clone()).or_default().push(f);
        }
    }
    out.into_iter().collect()
}

/// One volume mesh from the accepted parts: an element block and a
/// `<part>_boundary` side set each, plus a side set per OBJ group (one for
/// a group in several parts, spanning their blocks).
#[must_use]
pub fn assemble(parts: &[MeshedPart]) -> VolumeMesh {
    let mut mesh = VolumeMesh::default();
    let mut point_offsets = Vec::with_capacity(parts.len());
    for p in parts {
        let offset = mesh.points.len() as u32;
        point_offsets.push(offset);
        mesh.points.extend(&p.points);
        mesh.blocks.push(Block {
            name: p.name.clone(),
            tets: p.tets.iter().map(|t| t.map(|v| v + offset)).collect(),
        });
    }
    mesh.add_boundary_side_sets();
    let offsets = mesh.block_offsets();
    // The group side sets, by name, in order of first appearance.
    let mut named: BTreeMap<String, usize> = BTreeMap::new();
    let mut group_sets: Vec<SideSet> = Vec::new();
    for (b, p) in parts.iter().enumerate() {
        if p.face_groups.is_empty() {
            continue;
        }
        let elements = &mesh.blocks[b].tets;
        let mut by_face: HashMap<[u32; 3], (u32, u8)> = HashMap::new();
        for (e, side) in mesh.boundary_sides(b) {
            let local = (e - offsets[b]) as usize;
            if let Some(mut f) = elements.get(local).and_then(|&t| tet_side_nodes(t, side)) {
                f.sort_unstable();
                by_face.insert(f, (e, side));
            }
        }
        let offset = point_offsets[b];
        for (group, faces) in &p.face_groups {
            let sides: Vec<(u32, u8)> = faces
                .iter()
                .filter_map(|f| {
                    let mut k = f.map(|v| v + offset);
                    k.sort_unstable();
                    by_face.get(&k).copied()
                })
                .collect();
            if sides.is_empty() {
                continue;
            }
            let at = *named.entry(group.clone()).or_insert_with(|| {
                group_sets.push(SideSet {
                    name: group.clone(),
                    sides: Vec::new(),
                });
                group_sets.len() - 1
            });
            group_sets[at].sides.extend(sides);
        }
    }
    for mut set in group_sets {
        set.sides.sort_unstable();
        mesh.side_sets.push(set);
    }
    mesh
}

/// Signed volume enclosed by an outward-oriented closed surface.
fn enclosed_volume(points: &[Point3], faces: &[[u32; 3]]) -> f64 {
    let Some(r) = points.first().copied() else {
        return 0.0;
    };
    faces
        .iter()
        .map(|f| {
            let [a, b, c] = f.map(|v| sub(points[v as usize], r));
            dot(a, cross(b, c))
        })
        .sum::<f64>()
        / 6.0
}

/// Two-sided Hausdorff distance between the input surface `a` and the mesh
/// boundary `b`, sampled at vertices and face centroids. Input samples buried
/// inside the region (winding number of 1 or more, as on the inner faces of
/// overlapping shells) are skipped: the mesh bounds the region, not every
/// input face.
fn hausdorff(
    a_points: &[Point3],
    a_faces: &[[u32; 3]],
    winding: &FastWindingNumber<'_>,
    b_points: &[Point3],
    b_faces: &[[u32; 3]],
) -> Result<f64, String> {
    let a = TriangleSoup::new(a_points, a_faces).map_err(|e| e.to_string())?;
    let b = TriangleSoup::new(b_points, b_faces).map_err(|e| e.to_string())?;
    let (ba, bb) = (Bvh::new(a), Bvh::new(b));
    let samples = |points: &[Point3], faces: &[[u32; 3]]| -> Vec<Point3> {
        let mut s: Vec<Point3> = faces
            .iter()
            .flatten()
            .map(|&v| points[v as usize])
            .collect();
        s.extend(faces.iter().map(|f| {
            let [p, q, r] = f.map(|v| points[v as usize]);
            [0, 1, 2].map(|k| (p[k] + q[k] + r[k]) / 3.0)
        }));
        s
    };
    let one_sided = |from: &[Point3], to: &Bvh<'_>| {
        from.iter()
            .filter_map(|&p| to.nearest(p))
            .map(|n| n.distance_squared.sqrt())
            .fold(0.0, f64::max)
    };
    let bounding: Vec<Point3> = samples(a_points, a_faces)
        .into_iter()
        .filter(|&p| winding.winding_number(p).abs() < 1.0)
        .collect();
    Ok(one_sided(&bounding, &bb).max(one_sided(&samples(b_points, b_faces), &ba)))
}

#[cfg(test)]
mod tests {
    use super::face_manifold;

    #[test]
    fn face_manifold_counts_crowded_and_same_side_faces() {
        // One tet, and two on opposite sides of their shared face 0-1-2.
        assert_eq!(face_manifold(&[[0, 1, 2, 3]]), (0, 0));
        assert_eq!(face_manifold(&[[0, 1, 2, 3], [0, 2, 1, 4]]), (0, 0));
        // The same tet twice (listed from another vertex): all four faces
        // have both copies on the same side.
        assert_eq!(face_manifold(&[[0, 1, 2, 3], [1, 2, 0, 3]]), (0, 4));
        // Three tets on one face.
        assert_eq!(
            face_manifold(&[[0, 1, 2, 3], [0, 2, 1, 4], [0, 2, 1, 5]]),
            (1, 0)
        );
    }
}
