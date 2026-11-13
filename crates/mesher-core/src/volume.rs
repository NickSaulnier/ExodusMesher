//! Volume meshing and the volume gate.
//!
//! The fallback strategy of `PLAN.md`: a part whose repaired surface passed
//! the surface gate is meshed first by the conforming Delaunay kernel
//! (boundary-exact), or, with `--kernel exact`, by the exact constrained
//! kernel rounded to f64, the conforming kernel following if that fails. If that fails (budget, recovery) or the surface did not
//! pass the gate, isosurface stuffing on the winding-number field (`0.5 − w`)
//! takes over, unless `--strict` forbids it. A failed conformal attempt is a
//! warning when the fallback then succeeds, an error otherwise.
//!
//! The volume gate checks orientation, quality (the proven bounds for
//! stuffing; the sliver threshold for Delaunay) and the Hausdorff deviation
//! from the input surface, and records finite-volume quality metrics. Accepted
//! parts are assembled into one mesh: an element block and a `<part>_boundary`
//! side set per part, plus a side set per OBJ group on the boundary of a
//! conformal mesh.

use std::collections::{BTreeMap, HashMap};

use mesher_geom::predicates::orient3d;
use mesher_geom::vec3::{cross, dot, sub};
use mesher_geom::{Aabb, Bvh, FastWindingNumber, Point3, Sign, TriangleSoup};
use mesher_tet::cdt::{
    FaceError, RegionError, RoundingError, SegmentError, Untangle, mesh_region, round_region,
};
use mesher_tet::stuffing::{DIHEDRAL_MAX_DEGREES, DIHEDRAL_MIN_DEGREES};
use mesher_tet::{
    ConformingError, ConformingParams, ImproveParams, RemeshParams, SizeField, SizingParams,
    StuffingError, StuffingParams, conforming_mesh, fv_quality, improve, remesh, stuff,
    stuff_graded,
};

use crate::diagnose::{InternalError, Tolerances};
use crate::mesh::{Block, SideSet, VolumeMesh, tet_side_nodes};
use crate::model::SurfaceModel;
use crate::repair::{RepairedPart, passes_gate};
use crate::report::{
    Action, Attempt, AttemptOutcome, DefectCode, JobConfig, Kernel, KernelChoice, MeshSummary,
    Outcome, Record, RecordLog, Severity, Stage,
};

/// Slivers below this dihedral angle in a conformal mesh earn a quality warning.
pub const SLIVER_DEGREES: f64 = 10.0;

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
}

/// Mesh every repaired part.
pub fn mesh_parts(
    model: &SurfaceModel,
    repaired: &[RepairedPart],
    config: &JobConfig,
    log: &mut RecordLog,
) -> Result<MeshOutcome, InternalError> {
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
    let name = r.name.as_str();
    let internal = |e: String| InternalError(format!("part {name}: {e}"));
    let soup = TriangleSoup::new(&r.points, &r.faces).map_err(|e| internal(e.to_string()))?;
    let bounds = Aabb::of_points(&r.points);
    let diagonal = bounds.diagonal();
    let h = config.element_size * diagonal;
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
        let mut candidates: Vec<(Candidate, Option<usize>)> = Vec::new();
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
                Ok((triangles, c)) => {
                    let w = worst_cell(&c);
                    remeshed_worst = Some(w);
                    if exact_worst.is_some_and(|e| w > e) {
                        candidates.insert(0, (c, Some(triangles)));
                    } else {
                        candidates.push((c, Some(triangles)));
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
        for (i, (c, remeshed_triangles)) in candidates.into_iter().enumerate() {
            let (accepted, summary) = volume_gate(name, &c, r, &winding, diagonal, h, config, log)?;
            outcome.attempts.push(Attempt {
                remeshed: remeshed_triangles.is_some(),
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
                        if remeshed_triangles.is_some() {
                            Action::Remesh
                        } else {
                            Action::Fallback
                        },
                    );
                }
                if let Some(triangles) = remeshed_triangles {
                    // Findings deferred to remeshing are now dealt with.
                    log.deferred_applied(name, Action::Remesh);
                    let rec = match conformal_failure.take() {
                        // The repaired surface did not mesh: remeshing
                        // worked around it.
                        Some(rec) => {
                            rec.severity(Severity::Warning)
                                .measured("remeshed_triangles", triangles as u64)
                                .action(Action::Remesh, Outcome::Applied)
                        }
                        // It meshed, with poor cells: the boundary
                        // triangulation is no longer the input's, so say so.
                        None => {
                            Record::new(DefectCode::VolLowQuality, Stage::VolumeMesh)
                                .part(name)
                                .severity(Severity::Info)
                                .measured("input_worst_cell", exact_worst.unwrap_or(0.0))
                                .measured("remeshed_worst_cell", remeshed_worst.unwrap_or(0.0))
                                .measured("remeshed_triangles", triangles as u64)
                                .threshold("worst_cell", 1.0)
                                .action(Action::Remesh, Outcome::Applied)
                                .hint("the input's triangulation gave cells beyond the quality limits; the surface was remeshed to the element size, which gave better ones")
                        }
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
            if remeshed_triangles.is_some() {
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

    // 2. Fallback: isosurface stuffing. Where the size field asks for less
    // than the element size, the lattice is graded: uniformly fine at the
    // surface (h halved until it reaches the field's smallest size, at most
    // twice), coarsening to h inside, which keeps the proven angle bounds.
    let levels = if config.uniform_size {
        0
    } else {
        let field = SizeField::from_surface(&r.points, &r.faces, &SizingParams::for_size(h));
        if field.is_uniform() {
            0
        } else {
            libm::log2(h / field.min()).ceil().clamp(0.0, 2.0) as u32
        }
    };
    let leaf = h / f64::from(1u32 << levels);
    let level_fn = |p: Point3| 0.5 - winding.winding_number(p);
    let stuffed = if levels == 0 {
        stuff(level_fn, &bounds, &StuffingParams::with_cell_size(h))
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
            &bounds,
            &seeds,
            &StuffingParams::with_cell_size(leaf),
            levels,
        )
    };
    let stuffed = match stuffed {
        Ok(s) => Some(s),
        Err(StuffingError::Budget {
            lattice_vertices,
            limit,
        }) => {
            log.push(
                Record::new(DefectCode::MeshBudgetExceeded, Stage::VolumeMesh)
                    .part(name)
                    .measured("lattice_vertices", lattice_vertices)
                    .threshold("max_lattice_vertices", limit)
                    .hint("raise --element-size"),
            );
            outcome.attempts.push(attempt(Kernel::Stuffing, false));
            None
        }
        Err(e) => return Err(internal(e.to_string())),
    };
    let mut fallback_ok = false;
    if let Some(s) = stuffed {
        let c = Candidate {
            kernel: Kernel::Stuffing,
            points: s.points,
            tets: s.tets,
            dropped: s.dropped,
            steiner: 0,
            face_groups: Vec::new(),
            min_size: (levels > 0).then_some(leaf),
            remeshed: false,
        };
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
                face_groups: Vec::new(),
            });
        }
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
/// − its worst non-orthogonality) over (90° − the target of
/// `ImproveParams`), and the skewness target over its worst skewness. So 1
/// is at the limits and smaller is worse; 0 for a candidate with no cells.
fn worst_cell(c: &Candidate) -> f64 {
    if c.tets.is_empty() {
        return 0.0;
    }
    let targets = ImproveParams::default();
    let q = fv_quality(&c.points, &c.tets);
    let dihedral = q.min_dihedral_degrees / SLIVER_DEGREES;
    let non_orthogonality =
        (90.0 - q.max_non_orthogonality_degrees) / (90.0 - targets.max_non_orthogonality_degrees);
    let skewness = if q.max_skewness > 0.0 {
        targets.max_skewness / q.max_skewness
    } else {
        f64::INFINITY
    };
    let w = dihedral.min(non_orthogonality).min(skewness);
    if w.is_nan() { 0.0 } else { w }
}

fn improve_cells(c: &mut Candidate) {
    let r = improve(&c.points, &c.tets, &ImproveParams::default());
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
fn remeshed(
    model: &SurfaceModel,
    r: &RepairedPart,
    h: f64,
    config: &JobConfig,
) -> Result<(usize, Candidate), String> {
    let field = if config.uniform_size {
        SizeField::uniform(h)
    } else {
        SizeField::from_surface(&r.points, &r.faces, &SizingParams::for_size(h))
    };
    // Group boundaries (side sets) stay feature lines.
    let groups: Vec<u32> = r
        .face_origin
        .iter()
        .map(|o| {
            o.and_then(|t| model.polygon_of(t))
                .and_then(|p| p.group)
                .unwrap_or(u32::MAX)
        })
        .collect();
    let s = remesh(
        &r.points,
        &r.faces,
        &groups,
        &|p| field.at(p),
        &RemeshParams::default(),
    )
    .map_err(|e| e.to_string())?;
    let faces = s.triangles.len();
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
    if !passes_gate(&surface.points, &surface.faces, &Tolerances::default()) {
        return Err("the remeshed surface failed the surface gate".to_owned());
    }
    let mut c = conformal(model, &surface, h, config, true).map_err(|rec| {
        rec.hint
            .clone()
            .unwrap_or_else(|| "conformal meshing failed".to_owned())
    })?;
    c.remeshed = true;
    Ok((faces, c))
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
    let size = (!config.uniform_size)
        .then(|| SizeField::from_surface(&r.points, &r.faces, &SizingParams::for_size(h)))
        .filter(|f| !f.is_uniform());
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
    if !tets.is_empty() && (q.min_dihedral_degrees < lo || q.max_dihedral_degrees > hi) {
        log.push(
            gate(DefectCode::VolLowQuality)
                .measured("min_dihedral_degrees", q.min_dihedral_degrees)
                .measured("max_dihedral_degrees", q.max_dihedral_degrees)
                .measured("max_non_orthogonality_degrees", q.max_non_orthogonality_degrees)
                .threshold("min_dihedral_degrees", lo)
                .threshold("max_dihedral_degrees", hi)
                .hint("a few low-quality cells remain; a finite-volume solver will lose some accuracy there"),
        );
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
        hausdorff: finite(distance),
        hausdorff_relative: finite(relative),
    };
    Ok((accepted, summary))
}

/// One volume mesh from the accepted parts: an element block and a
/// `<part>_boundary` side set each, plus a side set per OBJ group.
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
            let mut sides: Vec<(u32, u8)> = faces
                .iter()
                .filter_map(|f| {
                    let mut k = f.map(|v| v + offset);
                    k.sort_unstable();
                    by_face.get(&k).copied()
                })
                .collect();
            sides.sort_unstable();
            if !sides.is_empty() {
                mesh.side_sets.push(SideSet {
                    name: group.clone(),
                    sides,
                });
            }
        }
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
