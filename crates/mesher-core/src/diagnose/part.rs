//! Per-part diagnosis: clean-up predictions (merged vertices, degenerate and
//! duplicate faces, T-junctions), then the shared surface analysis turned into
//! records.

use std::collections::BTreeMap;

use mesher_geom::vec3::{norm, sub};
use mesher_geom::{Aabb, GeomError, Point3};

use super::analysis::{Analysis, analyze_surface};
use super::parts::{Group, centroid, degenerate};
use super::{InternalError, Part, Tolerances};
use crate::model::SurfaceModel;
use crate::overlay::Overlay;
use crate::report::{Action, DefectCode, Outcome, Record, RecordLog, Severity, Stage};

pub(crate) struct Context<'a> {
    pub model: &'a SurfaceModel,
    pub rep: &'a [u32],
    pub tol: &'a Tolerances,
    pub merge_tolerance: f64,
    pub diagonal: f64,
}

pub(crate) fn analyze(
    ctx: &Context<'_>,
    group: Group,
    log: &mut RecordLog,
    overlay: &mut Overlay,
) -> Result<Part, InternalError> {
    let model = ctx.model;
    let points = &model.vertices;
    let name = group.name.as_str();
    let rec = |code| Record::new(code, Stage::Diagnose).part(name);
    let line = |t: u32| model.triangle_line(t);

    let bounds = group
        .triangles
        .iter()
        .flat_map(|&t| model.triangles[t as usize].vertices)
        .fold(Aabb::EMPTY, |b, v| b.with_point(points[v as usize]));

    // Duplicate vertices: clusters with several members inside this part.
    let mut clusters: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for &t in &group.triangles {
        for v in model.triangles[t as usize].vertices {
            clusters.entry(ctx.rep[v as usize]).or_default().push(v);
        }
    }
    let mut merged = Vec::new();
    for (&r, members) in clusters.iter_mut() {
        members.sort_unstable();
        members.dedup();
        if members.len() < 2 {
            continue;
        }
        let rp = points[r as usize];
        let spread = members
            .iter()
            .map(|&v| norm(sub(points[v as usize], rp)))
            .fold(0.0, f64::max);
        log.push(
            rec(DefectCode::GeoDuplicateVertex)
                .vertex_set(members.iter().copied())
                .lines(members.iter().map(|&v| model.vertex_line(v)))
                .point(rp)
                .measured("count", members.len() as u64)
                .measured("max_distance", spread)
                .threshold("merge_tolerance", ctx.merge_tolerance)
                .threshold("merge_tolerance_relative", ctx.tol.merge_relative)
                .action(Action::Merge, Outcome::Deferred),
        );
        overlay.vertices(members.iter().copied(), DefectCode::GeoDuplicateVertex);
        merged.push(members.clone());
    }

    // Degenerate faces; the rest are candidates.
    let mut candidates: Vec<([u32; 3], u32)> = Vec::new();
    let mut collinear: Vec<(u32, [u32; 3])> = Vec::new();
    let mut degenerate_faces: Vec<(u32, bool)> = Vec::new();
    for &t in &group.triangles {
        let f = model.triangles[t as usize]
            .vertices
            .map(|v| ctx.rep[v as usize]);
        if degenerate(points, f) {
            let repeated = f[0] == f[1] || f[1] == f[2] || f[0] == f[2];
            if !repeated {
                collinear.push((t, f));
            }
            degenerate_faces.push((t, false));
        } else {
            candidates.push((f, t));
        }
    }

    // Duplicates among the candidates: keep the first copy.
    let mut keyed: Vec<([u32; 3], u32, [u32; 3])> = candidates
        .iter()
        .map(|&(f, t)| {
            let mut k = f;
            k.sort_unstable();
            (k, t, f)
        })
        .collect();
    keyed.sort_unstable();
    let mut kept: Vec<([u32; 3], u32)> = Vec::with_capacity(keyed.len());
    let mut duplicates = Vec::new();
    let mut i = 0;
    while i < keyed.len() {
        let j = (i..keyed.len())
            .find(|&j| keyed[j].0 != keyed[i].0)
            .unwrap_or(keyed.len());
        kept.push((keyed[i].2, keyed[i].1));
        if j - i > 1 {
            let copies: Vec<u32> = keyed[i..j].iter().map(|k| k.1).collect();
            let opposed = keyed[i + 1..j]
                .iter()
                .any(|k| !same_winding(k.2, keyed[i].2));
            log.push(
                rec(DefectCode::GeoDuplicateFace)
                    .faces(copies.iter().copied())
                    .lines(copies.iter().map(|&t| line(t)))
                    .measured("copies", copies.len() as u64)
                    .measured("opposite_winding", opposed)
                    .action(Action::Remove, Outcome::Deferred)
                    .hint("the first copy is kept"),
            );
            overlay.faces(copies[1..].iter().copied(), DefectCode::GeoDuplicateFace);
            duplicates.push(copies);
        }
        i = j;
    }
    kept.sort_unstable_by_key(|k| k.1);

    // A collinear face is a collapsed sliver: removing it leaves a zero-area
    // slit unless the face across its long edge is split at the middle vertex.
    for &(t, f) in &collinear {
        let split = split_t_junction(points, f, &mut kept);
        if let Some(d) = degenerate_faces.iter_mut().find(|d| d.0 == t) {
            d.1 = split;
        }
    }
    for &(t, split) in &degenerate_faces {
        let orig = model.triangles[t as usize].vertices;
        let f = orig.map(|v| ctx.rep[v as usize]);
        let repeated = f[0] == f[1] || f[1] == f[2] || f[0] == f[2];
        log.push(
            rec(DefectCode::GeoDegenerateFace)
                .faces([t])
                .vertex_set(orig)
                .lines([line(t)])
                .measured(
                    "reason",
                    if repeated {
                        "repeated_vertex"
                    } else {
                        "collinear"
                    },
                )
                .measured("t_junction_split", split)
                .action(Action::Remove, Outcome::Deferred),
        );
        overlay.faces([t], DefectCode::GeoDegenerateFace);
    }

    let mut part = Part {
        name: group.name.clone(),
        source: group.source,
        triangles: group.triangles.clone(),
        faces: Vec::new(),
        face_ids: Vec::new(),
        vertices: clusters.len(),
        shells: 0,
        bounds,
        fatal: false,
        merged,
        degenerate: degenerate_faces,
        duplicates,
        reoriented: Vec::new(),
    };
    if kept.is_empty() {
        log.push(
            rec(DefectCode::TopOpenSheet)
                .faces(group.triangles.iter().copied())
                .lines(group.triangles.iter().map(|&t| line(t)))
                .measured("kept_faces", 0u64)
                .hint("every face is degenerate or a duplicate; nothing encloses a volume"),
        );
        part.fatal = true;
        return Ok(part);
    }

    // Compact to local indices and analyse.
    let mut verts: Vec<u32> = kept.iter().flat_map(|k| k.0).collect();
    verts.sort_unstable();
    verts.dedup();
    let local_of = |v: u32| verts.binary_search(&v).unwrap_or(0) as u32;
    let local_points: Vec<Point3> = verts.iter().map(|&v| points[v as usize]).collect();
    let local_tris: Vec<[u32; 3]> = kept.iter().map(|k| k.0.map(local_of)).collect();
    let ids: Vec<u32> = kept.iter().map(|k| k.1).collect();
    let a = match analyze_surface(&local_points, &local_tris, &local_tris, ctx.tol) {
        Ok(a) => a,
        Err(GeomError::TooLarge {
            vertices,
            triangles,
        }) => {
            log.push(
                rec(DefectCode::GeoCapacityExceeded)
                    .measured("vertices", vertices as u64)
                    .measured("triangles", triangles as u64),
            );
            part.fatal = true;
            return Ok(part);
        }
        Err(e) => return Err(InternalError(format!("part {name}: {e}"))),
    };
    let orig_v = |l: u32| verts[l as usize];
    let orig_f = |f: u32| ids[f as usize];
    report(
        ctx,
        name,
        &a,
        &local_points,
        &orig_v,
        &orig_f,
        log,
        overlay,
        &mut part,
    );

    part.shells = a.shell_faces.len();
    part.reoriented = a
        .disagreements
        .iter()
        .map(|d| d.minority.iter().map(|&f| orig_f(f)).collect())
        .collect();
    part.faces = a.oriented.iter().map(|t| t.map(orig_v)).collect();
    part.face_ids = ids.clone();
    Ok(part)
}

/// Turn an analysis into diagnosis records.
#[allow(clippy::too_many_arguments)]
fn report(
    ctx: &Context<'_>,
    name: &str,
    a: &Analysis,
    points: &[Point3],
    orig_v: &dyn Fn(u32) -> u32,
    orig_f: &dyn Fn(u32) -> u32,
    log: &mut RecordLog,
    overlay: &mut Overlay,
    part: &mut Part,
) {
    let model = ctx.model;
    let rec = |code| Record::new(code, Stage::Diagnose).part(name);
    let line = |t: u32| model.triangle_line(t);
    let faces_of = |fs: &[u32]| -> Vec<u32> { fs.iter().map(|&f| orig_f(f)).collect() };

    for &e in &a.nonmanifold_edges {
        let [u, v] = a.mesh.edge_vertices(e).unwrap_or([0, 0]);
        let faces: Vec<u32> = a
            .mesh
            .edge_half_edges(e)
            .iter()
            .map(|&h| orig_f(h / 3))
            .collect();
        log.push(
            rec(DefectCode::TopNonmanifoldEdge)
                .vertex_set([orig_v(u), orig_v(v)])
                .faces(faces.iter().copied())
                .lines(faces.iter().map(|&t| line(t)))
                .point(points[u as usize])
                .measured("incident_faces", faces.len() as u64)
                .threshold("max_incident_faces", 2u64)
                .action(Action::SplitEdge, Outcome::Deferred),
        );
        overlay.edge(orig_v(u), orig_v(v), DefectCode::TopNonmanifoldEdge);
    }

    for &(v, fans) in &a.nonmanifold_vertices {
        let faces = faces_of(&a.mesh.vertex_triangles(v));
        log.push(
            rec(DefectCode::TopNonmanifoldVertex)
                .vertex_set([orig_v(v)])
                .faces(faces)
                .lines([model.vertex_line(orig_v(v))])
                .point(points[v as usize])
                .measured("fans", fans as u64)
                .threshold("max_fans", 1u64)
                .action(Action::DuplicateVertex, Outcome::Deferred),
        );
        overlay.vertices([orig_v(v)], DefectCode::TopNonmanifoldVertex);
    }

    for d in &a.disagreements {
        let faces = faces_of(&d.minority);
        log.push(
            rec(DefectCode::TopInconsistentOrientation)
                .faces(faces.iter().copied())
                .lines(faces.iter().map(|&t| line(t)))
                .measured("faces_to_flip", faces.len() as u64)
                .measured("inconsistent_edges", d.inconsistent_edges.len() as u64)
                .measured("non_orientable", d.non_orientable)
                .action(Action::Reorient, Outcome::Deferred),
        );
        overlay.faces(faces, DefectCode::TopInconsistentOrientation);
        for &e in &d.inconsistent_edges {
            let [u, v] = a.mesh.edge_vertices(e).unwrap_or([0, 0]);
            overlay.edge(orig_v(u), orig_v(v), DefectCode::TopInconsistentOrientation);
        }
    }

    for s in 0..a.shell_faces.len() {
        if !a.sheet[s] {
            continue;
        }
        part.fatal = true;
        let faces = faces_of(&a.shell_faces[s]);
        log.push(
            rec(DefectCode::TopOpenSheet)
                .severity(Severity::Fatal)
                .faces(faces.iter().copied())
                .lines(faces.iter().map(|&t| line(t)))
                .measured("enclosed_volume_ratio", a.ratio[s])
                .threshold("min_enclosed_volume_ratio", ctx.tol.open_sheet_volume_ratio)
                .hint("this surface encloses no volume even with its holes closed"),
        );
        overlay.faces(faces, DefectCode::TopOpenSheet);
    }

    for &(f, g) in &a.self_pairs {
        let (x, y) = (orig_f(f), orig_f(g));
        let tf = a.mesh.triangle(f).unwrap_or([0, 0, 0]);
        log.push(
            rec(DefectCode::IntSelfIntersection)
                .faces([x, y])
                .lines([line(x), line(y)])
                .point(centroid(points, tf))
                .action(Action::Fallback, Outcome::Deferred),
        );
        overlay.faces([x, y], DefectCode::IntSelfIntersection);
    }
    for pairs in a.shell_overlaps.values() {
        let faces: Vec<u32> = pairs
            .iter()
            .flat_map(|&(f, g)| [orig_f(f), orig_f(g)])
            .collect();
        let tf = a.mesh.triangle(pairs[0].0).unwrap_or([0, 0, 0]);
        log.push(
            rec(DefectCode::IntShellOverlap)
                .faces(faces.iter().copied())
                .lines(faces.iter().map(|&t| line(t)))
                .point(centroid(points, tf))
                .measured("intersecting_pairs", pairs.len() as u64)
                .action(Action::Fallback, Outcome::Deferred),
        );
        overlay.faces(faces, DefectCode::IntShellOverlap);
    }

    for s in 0..a.shell_faces.len() {
        if a.sheet[s] || !a.inverted[s] {
            continue;
        }
        let faces = faces_of(&a.shell_faces[s]);
        log.push(
            rec(DefectCode::TopInvertedShell)
                .faces(faces.iter().copied())
                .lines(faces.iter().map(|&t| line(t)))
                .measured("signed_volume", a.volume[s])
                .measured("nesting_depth", a.depth[s])
                .action(Action::Flip, Outcome::Deferred),
        );
        overlay.faces(faces, DefectCode::TopInvertedShell);
    }

    for l in &a.loops {
        if a.sheet[l.shell] {
            continue;
        }
        let verts: Vec<u32> = l.edges.iter().flat_map(|e| e.0).map(orig_v).collect();
        let faces: Vec<u32> = l.edges.iter().map(|e| orig_f(e.1)).collect();
        let perimeter: f64 = l
            .edges
            .iter()
            .map(|&([i, j], _)| norm(sub(points[i as usize], points[j as usize])))
            .sum();
        let first = l.edges[0].0[0];
        log.push(
            rec(DefectCode::TopBoundaryEdge)
                .vertex_set(verts)
                .faces(faces.iter().copied())
                .lines(faces.iter().map(|&t| line(t)))
                .point(points[first as usize])
                .measured("edges", l.edges.len() as u64)
                .measured("perimeter", perimeter)
                .measured(
                    "perimeter_relative",
                    if ctx.diagonal > 0.0 {
                        perimeter / ctx.diagonal
                    } else {
                        0.0
                    },
                )
                .action(Action::FillHole, Outcome::Deferred),
        );
        for &([i, j], _) in &l.edges {
            overlay.edge(orig_v(i), orig_v(j), DefectCode::TopBoundaryEdge);
        }
    }

    for &(f, angle) in &a.slivers {
        let id = orig_f(f);
        log.push(
            rec(DefectCode::QualSliverFace)
                .faces([id])
                .lines([line(id)])
                .measured("min_angle_degrees", angle)
                .threshold("min_angle_degrees", ctx.tol.sliver_angle_degrees)
                .action(Action::Remesh, Outcome::Deferred),
        );
        overlay.faces([id], DefectCode::QualSliverFace);
    }
}

/// Split the single face across the long edge of collinear face `f` at its
/// middle vertex. Returns whether a split happened (not when the long edge has
/// no face, or several, across it).
fn split_t_junction(points: &[Point3], f: [u32; 3], kept: &mut Vec<([u32; 3], u32)>) -> bool {
    let pos = |v: u32| points[v as usize];
    let Some(k) = (0..3).find(|&k| f.iter().any(|&v| pos(v)[k] != pos(f[0])[k])) else {
        return false;
    };
    let mut order = f;
    order.sort_by(|&a, &b| pos(a)[k].total_cmp(&pos(b)[k]));
    let [p, m, q] = order;
    if pos(m)[k] == pos(p)[k] || pos(m)[k] == pos(q)[k] {
        return false;
    }
    let across: Vec<usize> = (0..kept.len())
        .filter(|&i| kept[i].0.contains(&p) && kept[i].0.contains(&q) && !kept[i].0.contains(&m))
        .collect();
    let [i] = across[..] else {
        return false;
    };
    let (t, id) = kept[i];
    let r = (0..3).find(|&c| t[c] != p && t[c] != q).unwrap_or(0);
    let [x, y] = [t[(r + 1) % 3], t[(r + 2) % 3]];
    let rv = t[r];
    kept[i] = ([x, m, rv], id);
    kept.insert(i + 1, ([m, y, rv], id));
    true
}

fn same_winding(a: [u32; 3], b: [u32; 3]) -> bool {
    (0..3).any(|r| [b[r], b[(r + 1) % 3], b[(r + 2) % 3]] == a)
}
