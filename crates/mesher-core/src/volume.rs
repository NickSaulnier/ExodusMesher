//! Volume meshing and the volume gate.
//!
//! Each repaired part is meshed by isosurface stuffing on its winding-number
//! field (`0.5 − w`), the robust fallback; the conformal Delaunay kernel
//! arrives in Milestones 7–8. The volume gate then checks orientation, the
//! dihedral-angle bounds, and the Hausdorff deviation from the input surface.
//! Parts that pass are assembled into one volume mesh: one element block and
//! one `<part>_boundary` side set per part.

use mesher_geom::predicates::orient3d;
use mesher_geom::vec3::{cross, dot, sub};
use mesher_geom::{Aabb, Bvh, FastWindingNumber, Point3, Sign, TriangleSoup};
use mesher_tet::stuffing::{DIHEDRAL_MAX_DEGREES, DIHEDRAL_MIN_DEGREES};
use mesher_tet::{StuffingError, StuffingParams, dihedral_range, stuff};

use crate::diagnose::InternalError;
use crate::mesh::{Block, VolumeMesh, tet_side_nodes};
use crate::repair::RepairedPart;
use crate::report::{
    Attempt, AttemptOutcome, DefectCode, JobConfig, Kernel, MeshSummary, Record, RecordLog, Stage,
};

/// A part's accepted volume mesh.
#[derive(Debug, Clone)]
pub struct MeshedPart {
    /// Position of the part in `Diagnosis::parts`.
    pub part: usize,
    pub name: String,
    pub points: Vec<Point3>,
    pub tets: Vec<[u32; 4]>,
}

#[derive(Debug, Clone, Default)]
pub struct MeshOutcome {
    /// Parts whose mesh passed the volume gate.
    pub meshed: Vec<MeshedPart>,
    /// Every attempted part's summary, by part position.
    pub summaries: Vec<(usize, MeshSummary)>,
    pub attempts: Vec<Attempt>,
}

/// Mesh every repaired part.
pub fn mesh_parts(
    repaired: &[RepairedPart],
    config: &JobConfig,
    log: &mut RecordLog,
) -> Result<MeshOutcome, InternalError> {
    let mut outcome = MeshOutcome::default();
    for r in repaired {
        mesh_part(r, config, log, &mut outcome)?;
    }
    Ok(outcome)
}

fn mesh_part(
    r: &RepairedPart,
    config: &JobConfig,
    log: &mut RecordLog,
    outcome: &mut MeshOutcome,
) -> Result<(), InternalError> {
    let name = r.name.as_str();
    if config.strict {
        log.push(
            Record::new(DefectCode::MeshStrictUnmet, Stage::VolumeMesh)
                .part(name)
                .hint("the conformal Delaunay kernel arrives in Milestone 8; without --strict the stuffing fallback is used"),
        );
        return Ok(());
    }
    let internal = |e: String| InternalError(format!("part {name}: {e}"));
    let soup = TriangleSoup::new(&r.points, &r.faces).map_err(|e| internal(e.to_string()))?;
    let bounds = Aabb::of_points(&r.points);
    let diagonal = bounds.diagonal();
    let h = config.element_size * diagonal;
    let winding = FastWindingNumber::new(soup);
    let attempt = |outcome: AttemptOutcome| Attempt {
        part: Some(name.to_owned()),
        kernel: Kernel::Stuffing,
        outcome,
        stdout: String::new(),
        stderr: String::new(),
    };

    let stuffed = match stuff(
        |p| 0.5 - winding.winding_number(p),
        &bounds,
        &StuffingParams::with_cell_size(h),
    ) {
        Ok(s) => s,
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
            outcome.attempts.push(attempt(AttemptOutcome::Failed));
            return Ok(());
        }
        Err(e) => return Err(internal(e.to_string())),
    };

    // Volume gate.
    let mut accepted = true;
    let gate = |code| Record::new(code, Stage::VolumeGate).part(name);
    let points = stuffed.points;
    let tets = stuffed.tets;
    let not_positive = tets
        .iter()
        .filter(|t| {
            let [a, b, c, d] = t.map(|v| points[v as usize]);
            orient3d(a, b, c, d) != Sign::Positive
        })
        .count() as u64
        + stuffed.dropped;
    if not_positive > 0 {
        accepted = false;
        log.push(
            gate(DefectCode::VolInvertedTet)
                .measured("tets", not_positive)
                .hint("isosurface stuffing should never do this; please report it"),
        );
    }
    let (min_angle, max_angle) = dihedral_range(&points, &tets).unwrap_or((0.0, 0.0));
    if !tets.is_empty() && (min_angle < DIHEDRAL_MIN_DEGREES || max_angle > DIHEDRAL_MAX_DEGREES) {
        log.push(
            gate(DefectCode::VolLowQuality)
                .measured("min_dihedral_degrees", min_angle)
                .measured("max_dihedral_degrees", max_angle)
                .threshold("min_dihedral_degrees", DIHEDRAL_MIN_DEGREES)
                .threshold("max_dihedral_degrees", DIHEDRAL_MAX_DEGREES),
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
    let hausdorff = if tets.is_empty() {
        f64::INFINITY
    } else {
        hausdorff(&r.points, &r.faces, &winding, &points, &boundary).map_err(internal)?
    };
    let relative = if diagonal > 0.0 {
        hausdorff / diagonal
    } else {
        0.0
    };
    if tets.is_empty() || relative > config.max_deviation {
        accepted = false;
        let mut rec = gate(DefectCode::VolBoundaryDeviation)
            .measured("tets", tets.len() as u64)
            .threshold("max_deviation_relative", config.max_deviation);
        rec = if tets.is_empty() {
            rec.hint("nothing was meshed: the part is thinner than the element size; lower --element-size")
        } else {
            rec.measured("hausdorff", hausdorff)
                .measured("hausdorff_relative", relative)
                .hint("lower --element-size for a closer fit, or raise --max-deviation")
        };
        log.push(rec);
    }

    let summary = MeshSummary {
        kernel: Kernel::Stuffing,
        accepted,
        tets: tets.len() as u64,
        nodes: points.len() as u64,
        element_size: h,
        volume: single.volume(),
        surface_volume: enclosed_volume(&r.points, &r.faces),
        min_dihedral_degrees: min_angle,
        max_dihedral_degrees: max_angle,
        hausdorff: if hausdorff.is_finite() {
            hausdorff
        } else {
            -1.0
        },
        hausdorff_relative: if relative.is_finite() { relative } else { -1.0 },
    };
    outcome.summaries.push((r.part, summary));
    outcome.attempts.push(attempt(if accepted {
        AttemptOutcome::Succeeded
    } else {
        AttemptOutcome::Failed
    }));
    if accepted {
        outcome.meshed.push(MeshedPart {
            part: r.part,
            name: name.to_owned(),
            points,
            tets,
        });
    }
    Ok(())
}

/// One volume mesh from the accepted parts: a block and a boundary side set each.
#[must_use]
pub fn assemble(parts: &[MeshedPart]) -> VolumeMesh {
    let mut mesh = VolumeMesh::default();
    for p in parts {
        let offset = mesh.points.len() as u32;
        mesh.points.extend(&p.points);
        mesh.blocks.push(Block {
            name: p.name.clone(),
            tets: p.tets.iter().map(|t| t.map(|v| v + offset)).collect(),
        });
    }
    mesh.add_boundary_side_sets();
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
