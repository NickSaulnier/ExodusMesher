//! Per-part checks: duplicate vertices, degenerate and duplicate faces,
//! non-manifold edges and vertices, orientation, shells (open sheets, inverted
//! shells, holes), self-intersections, shell overlaps and slivers.

use std::collections::BTreeMap;

use mesher_geom::vec3::{cross, dot, norm, sub};
use mesher_geom::{
    Aabb, Bvh, EdgeKind, FastWindingNumber, GeomError, HalfEdge, Point3, SurfaceMesh, TriangleSoup,
    faces_intersect,
};

use super::parts::{Group, centroid, degenerate};
use super::uf::DisjointSets;
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

/// Local (compacted) view of a part's kept faces.
struct Local {
    /// Merged-vertex id of each local vertex, ascending.
    verts: Vec<u32>,
    points: Vec<Point3>,
    tris: Vec<[u32; 3]>,
    /// Model triangle of each local face.
    ids: Vec<u32>,
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
    }
    let vertices = clusters.len();

    // Degenerate faces, then duplicates among the rest.
    let mut candidates: Vec<([u32; 3], u32)> = Vec::new();
    for &t in &group.triangles {
        let orig = model.triangles[t as usize].vertices;
        let f = orig.map(|v| ctx.rep[v as usize]);
        if degenerate(points, f) {
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
                    .action(Action::Remove, Outcome::Deferred),
            );
            overlay.faces([t], DefectCode::GeoDegenerateFace);
        } else {
            candidates.push((f, t));
        }
    }
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
        }
        i = j;
    }
    kept.sort_unstable_by_key(|k| k.1);

    let mut part = Part {
        name: group.name.clone(),
        source: group.source,
        triangles: group.triangles.clone(),
        faces: Vec::new(),
        face_ids: Vec::new(),
        vertices,
        shells: 0,
        bounds,
        fatal: false,
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

    // Compact to local indices.
    let mut verts: Vec<u32> = kept.iter().flat_map(|k| k.0).collect();
    verts.sort_unstable();
    verts.dedup();
    let local_of = |v: u32| verts.binary_search(&v).unwrap_or(0) as u32;
    let local = Local {
        points: verts.iter().map(|&v| points[v as usize]).collect(),
        tris: kept.iter().map(|k| k.0.map(local_of)).collect(),
        ids: kept.iter().map(|k| k.1).collect(),
        verts: verts.clone(),
    };
    let soup = match TriangleSoup::new(&local.points, &local.tris) {
        Ok(soup) => soup,
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
    let mesh = SurfaceMesh::new(&soup);
    let orig_v = |l: u32| local.verts[l as usize];
    let orig_f = |f: u32| local.ids[f as usize];
    let nf = local.tris.len();

    // Non-manifold edges.
    for e in mesh.edges() {
        if mesh.edge_kind(e) != Some(EdgeKind::NonManifold) {
            continue;
        }
        let [a, b] = mesh.edge_vertices(e).unwrap_or([0, 0]);
        let faces: Vec<u32> = mesh
            .edge_half_edges(e)
            .iter()
            .map(|&h| orig_f(h / 3))
            .collect();
        log.push(
            rec(DefectCode::TopNonmanifoldEdge)
                .vertex_set([orig_v(a), orig_v(b)])
                .faces(faces.iter().copied())
                .lines(faces.iter().map(|&t| line(t)))
                .point(local.points[a as usize])
                .measured("incident_faces", faces.len() as u64)
                .threshold("max_incident_faces", 2u64)
                .action(Action::SplitEdge, Outcome::Deferred),
        );
        overlay.edge(orig_v(a), orig_v(b), DefectCode::TopNonmanifoldEdge);
    }

    // Non-manifold vertices.
    for v in 0..mesh.vertex_count() {
        let fans = mesh.vertex_fan_count(v);
        if fans > 1 {
            let faces: Vec<u32> = mesh.vertex_triangles(v).into_iter().map(orig_f).collect();
            log.push(
                rec(DefectCode::TopNonmanifoldVertex)
                    .vertex_set([orig_v(v)])
                    .faces(faces.iter().copied())
                    .lines([model.vertex_line(orig_v(v))])
                    .point(local.points[v as usize])
                    .measured("fans", fans as u64)
                    .threshold("max_fans", 1u64)
                    .action(Action::DuplicateVertex, Outcome::Deferred),
            );
            overlay.vertices([orig_v(v)], DefectCode::TopNonmanifoldVertex);
        }
    }

    // Shells: faces connected through edges.
    let mut sets = DisjointSets::new(nf);
    for e in mesh.edges() {
        let hs = mesh.edge_half_edges(e);
        for w in hs.windows(2) {
            sets.union(w[0] / 3, w[1] / 3);
        }
    }
    let mut shell_index: BTreeMap<u32, usize> = BTreeMap::new();
    let mut shell_of = vec![0usize; nf];
    for (f, s) in shell_of.iter_mut().enumerate() {
        let root = sets.find(f as u32);
        let next = shell_index.len();
        *s = *shell_index.entry(root).or_insert(next);
    }
    let shells = shell_index.len();
    part.shells = shells;
    let mut shell_faces: Vec<Vec<u32>> = vec![Vec::new(); shells];
    for f in 0..nf {
        shell_faces[shell_of[f]].push(f as u32);
    }

    // Orientation: propagate across edges with exactly two faces; flip the
    // minority of each orientation component.
    let origin = |h: u32| mesh.endpoints(HalfEdge(h)).map_or(0, |e| e.0);
    let mut flag: Vec<Option<bool>> = vec![None; nf];
    let mut flip = vec![false; nf];
    for seed in 0..nf {
        if flag[seed].is_some() {
            continue;
        }
        flag[seed] = Some(false);
        let mut members = vec![seed as u32];
        let mut stack = vec![seed as u32];
        let mut conflict = false;
        let mut inconsistent_edges = Vec::new();
        while let Some(f) = stack.pop() {
            let ff = flag[f as usize].unwrap_or(false);
            for c in 0..3 {
                let h = f * 3 + c;
                let Some(e) = mesh.edge(HalfEdge(h)) else {
                    continue;
                };
                let [x, y] = match *mesh.edge_half_edges(e) {
                    [x, y] => [x, y],
                    _ => continue,
                };
                let other = if x == h { y } else { x };
                let same_direction = origin(other) == origin(h);
                if same_direction && h < other {
                    inconsistent_edges.push(e);
                }
                let want = ff ^ same_direction;
                let g = (other / 3) as usize;
                match flag[g] {
                    None => {
                        flag[g] = Some(want);
                        members.push(g as u32);
                        stack.push(g as u32);
                    }
                    Some(have) => conflict |= have != want,
                }
            }
        }
        let ones: Vec<u32> = members
            .iter()
            .copied()
            .filter(|&f| flag[f as usize] == Some(true))
            .collect();
        let minority_is_ones = ones.len() * 2 <= members.len();
        let minority: Vec<u32> = members
            .iter()
            .copied()
            .filter(|&f| flag[f as usize] == Some(minority_is_ones))
            .collect();
        for &f in &minority {
            flip[f as usize] = true;
        }
        if !minority.is_empty() {
            let faces: Vec<u32> = minority.iter().map(|&f| orig_f(f)).collect();
            log.push(
                rec(DefectCode::TopInconsistentOrientation)
                    .faces(faces.iter().copied())
                    .lines(faces.iter().map(|&t| line(t)))
                    .measured("faces_to_flip", faces.len() as u64)
                    .measured("inconsistent_edges", inconsistent_edges.len() as u64)
                    .measured("non_orientable", conflict)
                    .action(Action::Reorient, Outcome::Deferred),
            );
            overlay.faces(faces, DefectCode::TopInconsistentOrientation);
            for e in inconsistent_edges {
                let [a, b] = mesh.edge_vertices(e).unwrap_or([0, 0]);
                overlay.edge(orig_v(a), orig_v(b), DefectCode::TopInconsistentOrientation);
            }
        }
    }
    let oriented: Vec<[u32; 3]> = local
        .tris
        .iter()
        .zip(&flip)
        .map(|(&t, &f)| if f { [t[0], t[2], t[1]] } else { t })
        .collect();

    // Boundary loops: boundary edges of one shell joined through shared vertices.
    let boundary: Vec<(usize, [u32; 2], u32)> = mesh
        .edges()
        .filter(|&e| mesh.edge_kind(e) == Some(EdgeKind::Boundary))
        .filter_map(|e| {
            let h = *mesh.edge_half_edges(e).first()?;
            let (i, j) = mesh.endpoints(HalfEdge(h))?;
            let f = h / 3;
            let dir = if flip[f as usize] { [j, i] } else { [i, j] };
            Some((shell_of[f as usize], dir, f))
        })
        .collect();
    let mut loops = DisjointSets::new(boundary.len());
    let mut first_at: BTreeMap<(usize, u32), u32> = BTreeMap::new();
    for (k, &(s, [i, j], _)) in boundary.iter().enumerate() {
        for v in [i, j] {
            match first_at.get(&(s, v)) {
                Some(&other) => loops.union(other, k as u32),
                None => {
                    first_at.insert((s, v), k as u32);
                }
            }
        }
    }
    let mut loop_edges: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for k in 0..boundary.len() {
        loop_edges.entry(loops.find(k as u32)).or_default().push(k);
    }

    // Per shell: closed-up signed volume and area.
    let mut volume = vec![0.0f64; shells];
    let mut area = vec![0.0f64; shells];
    let reference: Vec<Point3> = shell_faces
        .iter()
        .map(|fs| local.points[local.tris[fs[0] as usize][0] as usize])
        .collect();
    for (f, t) in oriented.iter().enumerate() {
        let s = shell_of[f];
        let [a, b, c] = t.map(|v| local.points[v as usize]);
        volume[s] += tet(reference[s], a, b, c);
        area[s] += norm(cross(sub(b, a), sub(c, a))) / 2.0;
    }
    for edges in loop_edges.values() {
        let s = boundary[edges[0]].0;
        let mut loop_verts: Vec<u32> = edges.iter().flat_map(|&k| boundary[k].1).collect();
        loop_verts.sort_unstable();
        loop_verts.dedup();
        let n = loop_verts.len() as f64;
        let c = [0, 1, 2].map(|k| {
            loop_verts
                .iter()
                .map(|&v| local.points[v as usize][k])
                .sum::<f64>()
                / n
        });
        for &k in edges {
            let [i, j] = boundary[k].1;
            volume[s] += tet(
                reference[s],
                c,
                local.points[j as usize],
                local.points[i as usize],
            );
        }
    }

    // Open sheets.
    let mut sheet = vec![false; shells];
    for s in 0..shells {
        let ratio = if area[s] > 0.0 {
            volume[s].abs() / area[s].powf(1.5)
        } else {
            0.0
        };
        if ratio <= ctx.tol.open_sheet_volume_ratio {
            sheet[s] = true;
            part.fatal = true;
            let faces: Vec<u32> = shell_faces[s].iter().map(|&f| orig_f(f)).collect();
            log.push(
                rec(DefectCode::TopOpenSheet)
                    .severity(Severity::Fatal)
                    .faces(faces.iter().copied())
                    .lines(faces.iter().map(|&t| line(t)))
                    .measured("enclosed_volume_ratio", ratio)
                    .threshold("min_enclosed_volume_ratio", ctx.tol.open_sheet_volume_ratio)
                    .hint("this surface encloses no volume even with its holes closed"),
            );
            overlay.faces(faces, DefectCode::TopOpenSheet);
        }
    }

    // Intersections between faces beyond their shared vertices and edges.
    let bvh = Bvh::new(soup);
    let mut overlaps: BTreeMap<(usize, usize), Vec<(u32, u32)>> = BTreeMap::new();
    for (f, g) in bvh.self_overlaps() {
        let (tf, tg) = (local.tris[f as usize], local.tris[g as usize]);
        let pos = |t: [u32; 3]| t.map(|v| local.points[v as usize]);
        if !faces_intersect(tf, pos(tf), tg, pos(tg)) {
            continue;
        }
        let (sf, sg) = (shell_of[f as usize], shell_of[g as usize]);
        if sf == sg {
            let (a, b) = (orig_f(f), orig_f(g));
            log.push(
                rec(DefectCode::IntSelfIntersection)
                    .faces([a, b])
                    .lines([line(a), line(b)])
                    .point(centroid(&local.points, tf))
                    .action(Action::Fallback, Outcome::Deferred),
            );
            overlay.faces([a, b], DefectCode::IntSelfIntersection);
        } else {
            overlaps
                .entry((sf.min(sg), sf.max(sg)))
                .or_default()
                .push((f, g));
        }
    }
    for pairs in overlaps.values() {
        let faces: Vec<u32> = pairs
            .iter()
            .flat_map(|&(f, g)| [orig_f(f), orig_f(g)])
            .collect();
        log.push(
            rec(DefectCode::IntShellOverlap)
                .faces(faces.iter().copied())
                .lines(faces.iter().map(|&t| line(t)))
                .point(centroid(&local.points, local.tris[pairs[0].0 as usize]))
                .measured("intersecting_pairs", pairs.len() as u64)
                .action(Action::Fallback, Outcome::Deferred),
        );
        overlay.faces(faces, DefectCode::IntShellOverlap);
    }

    // Inverted shells: outward shells need positive volume, cavities (inside an
    // odd number of other shells) negative.
    let solid: Vec<usize> = (0..shells).filter(|&s| !sheet[s]).collect();
    let shell_tris: Vec<Vec<[u32; 3]>> = shell_faces
        .iter()
        .map(|fs| fs.iter().map(|&f| oriented[f as usize]).collect())
        .collect();
    for &s in &solid {
        let probe = centroid(&local.points, shell_tris[s][0]);
        let mut depth = 0u64;
        if solid.len() > 1 {
            for &t in &solid {
                // Nesting is only meaningful between disjoint shells; shells that
                // intersect are reported as overlapping instead.
                if t == s || overlaps.contains_key(&(s.min(t), s.max(t))) {
                    continue;
                }
                let soup = TriangleSoup::new(&local.points, &shell_tris[t])
                    .map_err(|e| InternalError(e.to_string()))?;
                if FastWindingNumber::new(soup).winding_number(probe).abs() > 0.5 {
                    depth += 1;
                }
            }
        }
        let want_positive = depth % 2 == 0;
        if (volume[s] > 0.0) != want_positive {
            let faces: Vec<u32> = shell_faces[s].iter().map(|&f| orig_f(f)).collect();
            log.push(
                rec(DefectCode::TopInvertedShell)
                    .faces(faces.iter().copied())
                    .lines(faces.iter().map(|&t| line(t)))
                    .measured("signed_volume", volume[s])
                    .measured("nesting_depth", depth)
                    .action(Action::Flip, Outcome::Deferred),
            );
            overlay.faces(faces, DefectCode::TopInvertedShell);
        }
    }

    // Holes, on shells that do enclose a volume.
    for edges in loop_edges.values() {
        let s = boundary[edges[0]].0;
        if sheet[s] {
            continue;
        }
        let verts: Vec<u32> = edges
            .iter()
            .flat_map(|&k| boundary[k].1)
            .map(orig_v)
            .collect();
        let faces: Vec<u32> = edges.iter().map(|&k| orig_f(boundary[k].2)).collect();
        let perimeter: f64 = edges
            .iter()
            .map(|&k| {
                let [i, j] = boundary[k].1;
                norm(sub(local.points[i as usize], local.points[j as usize]))
            })
            .sum();
        let first = boundary[edges[0]].1[0];
        log.push(
            rec(DefectCode::TopBoundaryEdge)
                .vertex_set(verts)
                .faces(faces.iter().copied())
                .lines(faces.iter().map(|&t| line(t)))
                .point(local.points[first as usize])
                .measured("edges", edges.len() as u64)
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
        for &k in edges {
            let [i, j] = boundary[k].1;
            overlay.edge(orig_v(i), orig_v(j), DefectCode::TopBoundaryEdge);
        }
    }

    // Slivers.
    for (f, t) in local.tris.iter().enumerate() {
        let angle = min_angle_degrees(t.map(|v| local.points[v as usize]));
        if angle < ctx.tol.sliver_angle_degrees {
            let id = orig_f(f as u32);
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

    part.faces = oriented.iter().map(|t| t.map(orig_v)).collect();
    part.face_ids = local.ids.clone();
    Ok(part)
}

fn same_winding(a: [u32; 3], b: [u32; 3]) -> bool {
    (0..3).any(|r| [b[r], b[(r + 1) % 3], b[(r + 2) % 3]] == a)
}

fn tet(r: Point3, a: Point3, b: Point3, c: Point3) -> f64 {
    dot(sub(a, r), cross(sub(b, r), sub(c, r))) / 6.0
}

fn min_angle_degrees(p: [Point3; 3]) -> f64 {
    (0..3)
        .map(|i| {
            let u = sub(p[(i + 1) % 3], p[i]);
            let v = sub(p[(i + 2) % 3], p[i]);
            norm(cross(u, v)).atan2(dot(u, v))
        })
        .fold(f64::INFINITY, f64::min)
        .to_degrees()
}
