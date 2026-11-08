//! Volume meshing and the volume gate.
//!
//! The fallback strategy of `PLAN.md`: a part whose repaired surface passed
//! the surface gate is meshed first by the conforming Delaunay kernel
//! (boundary-exact). If that fails (budget, recovery) or the surface did not
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
use mesher_tet::stuffing::{DIHEDRAL_MAX_DEGREES, DIHEDRAL_MIN_DEGREES};
use mesher_tet::{
    ConformingError, ConformingParams, RemeshParams, SizeField, SizingParams, StuffingError,
    StuffingParams, conforming_mesh, fv_quality, remesh, stuff,
};

use crate::diagnose::{InternalError, Tolerances};
use crate::mesh::{Block, SideSet, VolumeMesh, tet_side_nodes};
use crate::model::SurfaceModel;
use crate::repair::{RepairedPart, passes_gate};
use crate::report::{
    Action, Attempt, AttemptOutcome, DefectCode, JobConfig, Kernel, MeshSummary, Outcome, Record,
    RecordLog, Severity, Stage,
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

    // 1. Direct: conforming Delaunay, for a surface that passed the gate.
    let mut conformal_failure: Option<Record> = None;
    if r.summary.gate_passed && !config.stuffing_only {
        match conformal(model, r, h, config, false) {
            Ok(c) => {
                let (accepted, summary) =
                    volume_gate(name, &c, r, &winding, diagonal, h, config, log)?;
                outcome.attempts.push(attempt(Kernel::Delaunay, accepted));
                outcome.summaries.push((r.part, summary));
                if accepted {
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
                conformal_failure = Some(
                    Record::new(DefectCode::MeshRecoveryFailed, Stage::VolumeMesh)
                        .part(name)
                        .hint("the conformal mesh failed the volume gate"),
                );
            }
            Err(record) => {
                outcome.attempts.push(attempt(Kernel::Delaunay, false));
                conformal_failure = Some((*record).part(name));
            }
        }
    }

    // 1b. Remesh and retry once: a surface remeshed to the size field has
    // well-shaped, nearly Delaunay triangles, which the constrained kernel
    // recovers far more easily than many inputs. The mesh then conforms to
    // the remeshed surface; the volume gate measures its deviation from the
    // repaired one.
    if conformal_failure.is_some() && !config.no_remesh {
        match remeshed(model, r, h, config) {
            Ok((triangles, c)) => {
                let (accepted, summary) =
                    volume_gate(name, &c, r, &winding, diagonal, h, config, log)?;
                outcome.attempts.push(Attempt {
                    remeshed: true,
                    ..attempt(Kernel::Delaunay, accepted)
                });
                outcome.summaries.retain(|(i, _)| *i != r.part);
                outcome.summaries.push((r.part, summary));
                if accepted {
                    outcome.meshed.push(MeshedPart {
                        part: r.part,
                        name: name.to_owned(),
                        kernel: Kernel::Delaunay,
                        points: c.points,
                        tets: c.tets,
                        face_groups: c.face_groups,
                    });
                    if let Some(rec) = conformal_failure {
                        log.push(
                            rec.severity(Severity::Warning)
                                .measured("remeshed_triangles", triangles as u64)
                                .action(Action::Remesh, Outcome::Applied),
                        );
                    }
                    return Ok(());
                }
                conformal_failure = conformal_failure
                    .map(|rec| rec.measured("remesh_retry", "rejected by the volume gate"));
            }
            Err(reason) => {
                outcome.attempts.push(Attempt {
                    remeshed: true,
                    ..attempt(Kernel::Delaunay, false)
                });
                conformal_failure =
                    conformal_failure.map(|rec| rec.measured("remesh_retry", reason));
            }
        }
    }

    if config.strict {
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

    // 2. Fallback: isosurface stuffing.
    let stuffed = match stuff(
        |p| 0.5 - winding.winding_number(p),
        &bounds,
        &StuffingParams::with_cell_size(h),
    ) {
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
            min_size: None,
            remeshed: false,
        };
        let (accepted, summary) = volume_gate(name, &c, r, &winding, diagonal, h, config, log)?;
        outcome.attempts.push(attempt(Kernel::Stuffing, accepted));
        outcome.summaries.retain(|(i, _)| *i != r.part);
        outcome.summaries.push((r.part, summary));
        if accepted {
            fallback_ok = true;
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
    let params = ConformingParams {
        groups: if own_facets {
            (0..r.faces.len() as u32).collect()
        } else {
            (0..r.faces.len())
                .map(|f| group_of(f).unwrap_or(u32::MAX))
                .collect()
        },
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
    Ok(Candidate {
        kernel: Kernel::Delaunay,
        points: m.points,
        tets: m.tets,
        dropped: 0,
        steiner: (m.boundary_steiner + m.interior_steiner) as u64,
        face_groups: groups.into_iter().collect(),
        min_size,
        remeshed: false,
    })
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
    let (lo, hi) = match c.kernel {
        Kernel::Stuffing => (DIHEDRAL_MIN_DEGREES, DIHEDRAL_MAX_DEGREES),
        Kernel::Delaunay => (SLIVER_DEGREES, 180.0 - SLIVER_DEGREES),
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
    let summary = MeshSummary {
        kernel: c.kernel,
        accepted,
        tets: tets.len() as u64,
        nodes: points.len() as u64,
        steiner_points: c.steiner,
        element_size: h,
        min_element_size: c.min_size,
        remeshed: c.remeshed,
        volume: single.volume(),
        surface_volume: enclosed_volume(&r.points, &r.faces),
        min_dihedral_degrees: q.min_dihedral_degrees,
        max_dihedral_degrees: q.max_dihedral_degrees,
        dihedral_histogram: q.dihedral_histogram.clone(),
        max_non_orthogonality_degrees: q.max_non_orthogonality_degrees,
        mean_non_orthogonality_degrees: q.mean_non_orthogonality_degrees,
        max_skewness: q.max_skewness,
        max_volume_ratio: q.max_volume_ratio,
        hausdorff: if distance.is_finite() { distance } else { -1.0 },
        hausdorff_relative: if relative.is_finite() { relative } else { -1.0 },
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
