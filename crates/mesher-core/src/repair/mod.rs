//! Repair: conservative, ordered fixes to each selected part, every one logged,
//! followed by guards and the surface gate.
//!
//! Order, least invasive first:
//! 1. merge duplicate vertices, remove degenerate faces (splitting T-junctions
//!    they leave) and duplicate faces, reorient faces against their
//!    neighbours: diagnosis already computed this surface, so repair logs it;
//! 2. flip inverted shells;
//! 3. split non-manifold edges and vertices;
//! 4. fill small holes.
//!
//! Then two guards (the shell count and the enclosed volume must not change
//! beyond a limit) and the surface gate: watertight, edge- and vertex-manifold,
//! consistently and outwardly oriented, enclosing a volume, free of
//! intersections. Repair records have stage `repair`: severity `info` for a
//! fix applied (the finding itself already carries the warning) and `error`
//! for a fix that could not be made. Gate failures have stage `surface_gate`.

mod fill;
mod nonmanifold;

use std::collections::BTreeMap;

use mesher_geom::{GeomError, Point3};
use serde_json::Value;

use crate::diagnose::analysis::{Analysis, analyze_surface};
use crate::diagnose::uf::DisjointSets;
use crate::diagnose::{Diagnosis, InternalError, Part, Tolerances};
use crate::model::SurfaceModel;
use crate::report::{
    Action, DefectCode, Outcome, Record, RecordLog, RepairSummary, Severity, Stage,
};

/// A part after repair: its own vertex list, with provenance.
#[derive(Debug, Clone)]
pub struct RepairedPart {
    /// Position of the part in `Diagnosis::parts`.
    pub part: usize,
    pub name: String,
    pub points: Vec<Point3>,
    pub faces: Vec<[u32; 3]>,
    /// Model vertex (merged id) each vertex came from; `None` for a vertex
    /// added by hole filling.
    pub vertex_origin: Vec<Option<u32>>,
    /// Model triangle each face came from; `None` for a fill face.
    pub face_origin: Vec<Option<u32>>,
    pub summary: RepairSummary,
}

/// Repair every selected part that is not fatally defective.
pub fn repair(
    model: &SurfaceModel,
    diagnosis: &Diagnosis,
    selected: &[usize],
    tol: &Tolerances,
    log: &mut RecordLog,
) -> Result<Vec<RepairedPart>, InternalError> {
    selected
        .iter()
        .filter_map(|&i| diagnosis.parts.get(i).map(|p| (i, p)))
        .filter(|(_, p)| !p.fatal && !p.faces.is_empty())
        .map(|(i, p)| repair_part(model, p, i, tol, log))
        .collect()
}

fn repair_part(
    model: &SurfaceModel,
    part: &Part,
    index: usize,
    tol: &Tolerances,
    log: &mut RecordLog,
) -> Result<RepairedPart, InternalError> {
    let name = part.name.as_str();
    let line = |t: u32| model.triangle_line(t);
    let applied = |code, action| {
        Record::new(code, Stage::Repair)
            .part(name)
            .severity(Severity::Info)
            .action(action, Outcome::Applied)
    };
    let failed = |code, action| {
        Record::new(code, Stage::Repair)
            .part(name)
            .severity(Severity::Error)
            .action(action, Outcome::Failed)
    };
    let err = |e: GeomError| InternalError(format!("part {name}: {e}"));

    // 1. Clean-up, already applied to `part.faces` by diagnosis.
    for members in &part.merged {
        log.push(
            applied(DefectCode::GeoDuplicateVertex, Action::Merge)
                .vertex_set(members.iter().copied())
                .lines(members.iter().map(|&v| model.vertex_line(v)))
                .measured("vertices_before", members.len() as u64)
                .measured("vertices_after", 1u64),
        );
    }
    for &(t, split) in &part.degenerate {
        log.push(
            applied(DefectCode::GeoDegenerateFace, Action::Remove)
                .faces([t])
                .lines([line(t)])
                .measured("faces_removed", 1 + u64::from(split))
                .measured("faces_added", 2 * u64::from(split))
                .measured("t_junction_split", split),
        );
    }
    for copies in &part.duplicates {
        log.push(
            applied(DefectCode::GeoDuplicateFace, Action::Remove)
                .faces(copies.iter().copied())
                .lines(copies.iter().map(|&t| line(t)))
                .measured("faces_removed", copies.len() as u64 - 1),
        );
    }
    for flipped in part.reoriented.iter().filter(|f| !f.is_empty()) {
        log.push(
            applied(DefectCode::TopInconsistentOrientation, Action::Reorient)
                .faces(flipped.iter().copied())
                .lines(flipped.iter().map(|&t| line(t)))
                .measured("faces_flipped", flipped.len() as u64),
        );
    }

    // Local copy of the cleaned surface.
    let mut verts: Vec<u32> = part.faces.iter().flatten().copied().collect();
    verts.sort_unstable();
    verts.dedup();
    let local = |v: u32| verts.binary_search(&v).unwrap_or(0) as u32;
    let mut points: Vec<Point3> = verts.iter().map(|&v| model.vertices[v as usize]).collect();
    let mut origin: Vec<Option<u32>> = verts.iter().map(|&v| Some(v)).collect();
    let mut faces: Vec<[u32; 3]> = part.faces.iter().map(|t| t.map(local)).collect();
    let mut face_origin: Vec<Option<u32>> = part.face_ids.iter().map(|&t| Some(t)).collect();
    let vertices_before = points.len() as u64;
    let n_model = model.vertices.len() as u32;
    let contact = |faces: &[[u32; 3]], origin: &[Option<u32>]| -> Vec<[u32; 3]> {
        faces
            .iter()
            .map(|t| {
                t.map(|v| {
                    origin
                        .get(v as usize)
                        .copied()
                        .flatten()
                        .unwrap_or(n_model.saturating_add(v))
                })
            })
            .collect()
    };
    let origin_faces = |fs: &[u32], face_origin: &[Option<u32>]| -> Vec<u32> {
        fs.iter()
            .filter_map(|&f| face_origin.get(f as usize).copied().flatten())
            .collect()
    };

    let a0 = analyze_surface(&points, &faces, &contact(&faces, &origin), tol).map_err(err)?;
    let volume_before: f64 = a0.volume.iter().map(|v| v.abs()).sum();
    let shells_before = components(&contact(&faces, &origin));
    faces.clone_from(&a0.oriented);

    // 2. Flip inverted shells.
    for s in 0..a0.shell_faces.len() {
        if a0.sheet[s] || !a0.inverted[s] {
            continue;
        }
        for &f in &a0.shell_faces[s] {
            faces[f as usize].swap(1, 2);
        }
        let fs = origin_faces(&a0.shell_faces[s], &face_origin);
        log.push(
            applied(DefectCode::TopInvertedShell, Action::Flip)
                .faces(fs.iter().copied())
                .lines(fs.iter().map(|&t| line(t)))
                .measured("faces_flipped", a0.shell_faces[s].len() as u64)
                .measured("signed_volume_before", a0.volume[s])
                .measured("signed_volume_after", -a0.volume[s]),
        );
    }

    // 3. Split non-manifold edges and vertices.
    let (edge_splits, vertex_splits) =
        nonmanifold::split(&mut points, &mut faces, &mut origin).map_err(err)?;
    let vorigin = |v: u32, origin: &[Option<u32>]| origin.get(v as usize).copied().flatten();
    for s in &edge_splits {
        let fs = origin_faces(&s.faces, &face_origin);
        let vs: Vec<u32> = s.edge.iter().filter_map(|&v| vorigin(v, &origin)).collect();
        let record = if s.paired {
            applied(DefectCode::TopNonmanifoldEdge, Action::SplitEdge)
                .measured("sheets", s.faces.len() as u64 / 2)
                .measured("vertices_added", s.vertices_added)
        } else {
            failed(DefectCode::TopNonmanifoldEdge, Action::SplitEdge)
                .hint("the faces around this edge do not pair up into separate solids")
        };
        log.push(
            record
                .vertex_set(vs)
                .faces(fs.iter().copied())
                .lines(fs.iter().map(|&t| line(t))),
        );
    }
    for s in &vertex_splits {
        let fs = origin_faces(&s.faces, &face_origin);
        log.push(
            applied(DefectCode::TopNonmanifoldVertex, Action::DuplicateVertex)
                .vertex_set(vorigin(s.vertex, &origin))
                .faces(fs)
                .measured("fans", s.fans as u64)
                .measured("vertices_added", s.fans as u64 - 1),
        );
    }

    // 4. Fill small holes.
    let a1 = analyze_surface(&points, &faces, &contact(&faces, &origin), tol).map_err(err)?;
    for hole in &a1.loops {
        if a1.sheet[hole.shell] {
            continue;
        }
        let vs: Vec<u32> = hole
            .edges
            .iter()
            .flat_map(|e| e.0)
            .filter_map(|v| vorigin(v, &origin))
            .collect();
        let fs = origin_faces(
            &hole.edges.iter().map(|e| e.1).collect::<Vec<_>>(),
            &face_origin,
        );
        let base = |r: Record| {
            r.vertex_set(vs.iter().copied())
                .faces(fs.iter().copied())
                .lines(fs.iter().map(|&t| line(t)))
                .measured("edges", hole.edges.len() as u64)
                .threshold("max_hole_edges", tol.max_hole_edges as u64)
                .threshold("max_hole_fraction", tol.max_hole_fraction)
        };
        match fill::fill_loop(
            &mut points,
            &mut faces,
            &mut origin,
            hole,
            a1.area[hole.shell],
            tol,
        ) {
            fill::Fill::Filled {
                faces_added,
                vertex_added,
                area,
                fraction,
            } => {
                face_origin.extend(std::iter::repeat_n(None, faces_added as usize));
                log.push(
                    base(applied(DefectCode::TopBoundaryEdge, Action::FillHole))
                        .measured("faces_added", faces_added)
                        .measured("vertices_added", u64::from(vertex_added))
                        .measured("area", area)
                        .measured("hole_fraction", fraction),
                );
            }
            fill::Fill::Failed { reason, fraction } => {
                let mut r = base(failed(DefectCode::TopBoundaryEdge, Action::FillHole))
                    .measured("reason", reason);
                if let Some(fraction) = fraction {
                    r = r.measured("hole_fraction", fraction);
                }
                log.push(r.hint("left open; a hole this size may be an intentional opening"));
            }
        }
    }

    // Guards, then the gate.
    let ids = contact(&faces, &origin);
    let a2 = analyze_surface(&points, &faces, &ids, tol).map_err(err)?;
    let shells_after = components(&ids);
    if shells_after != shells_before {
        log.push(
            Record::new(DefectCode::RepairGuardExceeded, Stage::Repair)
                .part(name)
                .measured("guard", "shell_count")
                .measured("shells_before", shells_before as u64)
                .measured("shells_after", shells_after as u64),
        );
    }
    let volume_after: f64 = a2.volume.iter().map(|v| v.abs()).sum();
    let change = if volume_before > 0.0 {
        (volume_after - volume_before).abs() / volume_before
    } else {
        0.0
    };
    if change > tol.max_volume_change {
        log.push(
            Record::new(DefectCode::RepairGuardExceeded, Stage::Repair)
                .part(name)
                .measured("guard", "enclosed_volume")
                .measured("volume_before", volume_before)
                .measured("volume_after", volume_after)
                .measured("relative_change", change)
                .threshold("max_relative_change", tol.max_volume_change),
        );
    }
    let gate_passed = gate(name, &a2, &origin, &face_origin, model, log);

    let summary = RepairSummary {
        faces_before: part.triangles.len() as u64,
        faces_after: faces.len() as u64,
        vertices_before,
        vertices_after: points.len() as u64,
        gate_passed,
    };
    Ok(RepairedPart {
        part: index,
        name: part.name.clone(),
        points,
        faces,
        vertex_origin: origin,
        face_origin,
        summary,
    })
}

/// Whether a surface passes every surface-gate check, without recording
/// anything: for a surface built after repair, such as a remeshed one.
pub(crate) fn passes_gate(points: &[Point3], faces: &[[u32; 3]], tol: &Tolerances) -> bool {
    let Ok(a) = analyze_surface(points, faces, faces, tol) else {
        return false;
    };
    a.loops.iter().all(|l| a.sheet[l.shell])
        && a.nonmanifold_edges.is_empty()
        && a.nonmanifold_vertices.is_empty()
        && a.disagreements.is_empty()
        && !a.sheet.iter().any(|&s| s)
        && !a.inverted.iter().any(|&i| i)
        && a.self_pairs.is_empty()
        && a.shell_overlaps.is_empty()
}

/// Connected pieces of a surface by its vertex ids.
fn components(faces: &[[u32; 3]]) -> usize {
    let mut dense: BTreeMap<u32, u32> = BTreeMap::new();
    for &v in faces.iter().flatten() {
        let next = dense.len() as u32;
        dense.entry(v).or_insert(next);
    }
    let mut sets = DisjointSets::new(dense.len());
    for t in faces {
        let [a, b, c] = t.map(|v| dense.get(&v).copied().unwrap_or(0));
        sets.union(a, b);
        sets.union(b, c);
    }
    (0..dense.len() as u32)
        .filter(|&i| sets.find(i) == i)
        .count()
}

/// The surface gate. Records one `surface_gate` error per failed check and
/// returns whether every check passed.
fn gate(
    name: &str,
    a: &Analysis,
    origin: &[Option<u32>],
    face_origin: &[Option<u32>],
    model: &SurfaceModel,
    log: &mut RecordLog,
) -> bool {
    let mut passed = true;
    let mut fail =
        |code: DefectCode, faces: Vec<u32>, vertices: Vec<u32>, measured: Vec<(&str, Value)>| {
            passed = false;
            let fs: Vec<u32> = faces
                .iter()
                .filter_map(|&f| face_origin.get(f as usize).copied().flatten())
                .collect();
            let vs: Vec<u32> = vertices
                .iter()
                .filter_map(|&v| origin.get(v as usize).copied().flatten())
                .collect();
            let mut r = Record::new(code, Stage::SurfaceGate)
                .part(name)
                .severity(Severity::Error)
                .faces(fs.iter().copied())
                .vertex_set(vs)
                .lines(fs.iter().map(|&t| model.triangle_line(t)))
                .hint("the repaired surface is not yet a valid closed solid");
            for (k, v) in measured {
                r = r.measured(k, v);
            }
            log.push(r);
        };

    let open: Vec<_> = a.loops.iter().filter(|l| !a.sheet[l.shell]).collect();
    if !open.is_empty() {
        fail(
            DefectCode::TopBoundaryEdge,
            open.iter()
                .flat_map(|l| l.edges.iter().map(|e| e.1))
                .collect(),
            open.iter()
                .flat_map(|l| l.edges.iter().flat_map(|e| e.0))
                .collect(),
            vec![
                ("holes", open.len().into()),
                (
                    "edges",
                    open.iter().map(|l| l.edges.len()).sum::<usize>().into(),
                ),
            ],
        );
    }
    if !a.nonmanifold_edges.is_empty() {
        fail(
            DefectCode::TopNonmanifoldEdge,
            a.nonmanifold_edges
                .iter()
                .flat_map(|&e| a.mesh.edge_half_edges(e).iter().map(|h| h / 3))
                .collect(),
            a.nonmanifold_edges
                .iter()
                .filter_map(|&e| a.mesh.edge_vertices(e))
                .flatten()
                .collect(),
            vec![("edges", a.nonmanifold_edges.len().into())],
        );
    }
    if !a.nonmanifold_vertices.is_empty() {
        fail(
            DefectCode::TopNonmanifoldVertex,
            Vec::new(),
            a.nonmanifold_vertices.iter().map(|v| v.0).collect(),
            vec![("vertices", a.nonmanifold_vertices.len().into())],
        );
    }
    if !a.disagreements.is_empty() {
        fail(
            DefectCode::TopInconsistentOrientation,
            a.disagreements
                .iter()
                .flat_map(|d| d.minority.iter().copied())
                .collect(),
            Vec::new(),
            vec![("components", a.disagreements.len().into())],
        );
    }
    let sheets: Vec<usize> = (0..a.shell_faces.len()).filter(|&s| a.sheet[s]).collect();
    if !sheets.is_empty() {
        fail(
            DefectCode::TopOpenSheet,
            sheets
                .iter()
                .flat_map(|&s| a.shell_faces[s].iter().copied())
                .collect(),
            Vec::new(),
            vec![("shells", sheets.len().into())],
        );
    }
    let inverted: Vec<usize> = (0..a.shell_faces.len())
        .filter(|&s| !a.sheet[s] && a.inverted[s])
        .collect();
    if !inverted.is_empty() {
        fail(
            DefectCode::TopInvertedShell,
            inverted
                .iter()
                .flat_map(|&s| a.shell_faces[s].iter().copied())
                .collect(),
            Vec::new(),
            vec![("shells", inverted.len().into())],
        );
    }
    if !a.self_pairs.is_empty() {
        fail(
            DefectCode::IntSelfIntersection,
            a.self_pairs.iter().flat_map(|&(f, g)| [f, g]).collect(),
            Vec::new(),
            vec![("intersecting_pairs", a.self_pairs.len().into())],
        );
    }
    if !a.shell_overlaps.is_empty() {
        fail(
            DefectCode::IntShellOverlap,
            a.shell_overlaps
                .values()
                .flatten()
                .flat_map(|&(f, g)| [f, g])
                .collect(),
            Vec::new(),
            vec![("shell_pairs", a.shell_overlaps.len().into())],
        );
    }
    passed
}
