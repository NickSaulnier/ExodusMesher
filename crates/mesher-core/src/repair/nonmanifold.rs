//! Splitting non-manifold edges and vertices.
//!
//! At an edge with more than two faces, the faces are sorted by angle around
//! the edge (exactly, with `orient3d`) and paired so that each pair bounds the
//! material of one solid: a face traversing the edge `v→u` followed
//! counter-clockwise (around `u→v`) by one traversing it `u→v`. Faces then
//! connect only within their pair. At every vertex, the faces connected
//! through edges at the vertex form fans; each fan beyond the first gets its
//! own copy of the vertex. Vertex copies keep their origin, so later stages
//! still know which faces met there.

use std::cmp::Ordering;

use mesher_geom::predicates::{orient2d_projected, orient3d, plane_axis};
use mesher_geom::{EdgeId, GeomError, HalfEdge, Point3, Sign, SurfaceMesh, TriangleSoup};

use crate::diagnose::uf::DisjointSets;

/// What happened at one non-manifold edge.
#[derive(Debug, Clone)]
pub(crate) struct EdgeSplit {
    /// Local vertex ids of the edge.
    pub edge: [u32; 2],
    pub faces: Vec<u32>,
    /// Whether the faces could be paired into separate solids.
    pub paired: bool,
    /// Vertex copies made because of this edge's split.
    pub vertices_added: u64,
}

/// What happened at one vertex joining several fans (diagnosed as non-manifold).
#[derive(Debug, Clone)]
pub(crate) struct VertexSplit {
    pub vertex: u32,
    pub fans: usize,
    pub faces: Vec<u32>,
}

pub(crate) fn split(
    points: &mut Vec<Point3>,
    faces: &mut [[u32; 3]],
    origin: &mut Vec<Option<u32>>,
) -> Result<(Vec<EdgeSplit>, Vec<VertexSplit>), GeomError> {
    let snapshot = points.clone();
    let corners = faces.to_vec();
    let soup = TriangleSoup::new(&snapshot, &corners)?;
    let mesh = SurfaceMesh::new(&soup);

    // Face links across each edge.
    let mut links: Vec<Vec<(u32, u32)>> = vec![Vec::new(); mesh.edge_count() as usize];
    let mut edge_splits = Vec::new();
    for e in mesh.edges() {
        let hs = mesh.edge_half_edges(e);
        let link = &mut links[e.0 as usize];
        match hs.len() {
            0 | 1 => {}
            2 => link.push((hs[0] / 3, hs[1] / 3)),
            _ => {
                let pairs = pair_around_edge(&snapshot, &mesh, e);
                let paired = pairs.is_some();
                match pairs {
                    Some(pairs) => link.extend(pairs),
                    None => link.extend(hs.windows(2).map(|w| (w[0] / 3, w[1] / 3))),
                }
                edge_splits.push(EdgeSplit {
                    edge: mesh.edge_vertices(e).unwrap_or([0, 0]),
                    faces: hs.iter().map(|&h| h / 3).collect(),
                    paired,
                    vertices_added: 0,
                });
            }
        }
    }

    let mut vertex_splits = Vec::new();
    for v in 0..mesh.vertex_count() {
        let tris = mesh.vertex_triangles(v);
        if tris.len() < 2 {
            continue;
        }
        let local = |f: u32| tris.binary_search(&f).ok().map(|i| i as u32);
        let mut sets = DisjointSets::new(tris.len());
        let mut incident_nonmanifold = None;
        for &h in mesh.outgoing(v) {
            for half in [HalfEdge(h), HalfEdge(h).prev()] {
                let Some(e) = mesh.edge(half) else {
                    continue;
                };
                if mesh.edge_half_edges(e).len() > 2 && incident_nonmanifold.is_none() {
                    incident_nonmanifold = mesh.edge_vertices(e);
                }
                for &(f, g) in &links[e.0 as usize] {
                    if let (Some(a), Some(b)) = (local(f), local(g)) {
                        sets.union(a, b);
                    }
                }
            }
        }
        let mut fans: Vec<Vec<u32>> = Vec::new();
        let mut root_of: Vec<(u32, usize)> = Vec::new();
        for (i, &f) in tris.iter().enumerate() {
            let root = sets.find(i as u32);
            match root_of.iter().find(|r| r.0 == root) {
                Some(&(_, k)) => fans[k].push(f),
                None => {
                    root_of.push((root, fans.len()));
                    fans.push(vec![f]);
                }
            }
        }
        if fans.len() < 2 {
            continue;
        }
        for fan in &fans[1..] {
            let copy = points.len() as u32;
            points.push(snapshot[v as usize]);
            origin.push(origin.get(v as usize).copied().flatten());
            for &f in fan {
                for c in &mut faces[f as usize] {
                    if *c == v {
                        *c = copy;
                    }
                }
            }
        }
        let added = fans.len() as u64 - 1;
        if mesh.vertex_fan_count(v) > 1 {
            vertex_splits.push(VertexSplit {
                vertex: v,
                fans: fans.len(),
                faces: tris.clone(),
            });
        } else if let Some(edge) = incident_nonmanifold {
            if let Some(s) = edge_splits.iter_mut().find(|s| s.edge == edge) {
                s.vertices_added += added;
            }
        }
    }
    Ok((edge_splits, vertex_splits))
}

/// Pair the faces around non-manifold edge `e` into solids, or `None` if they
/// do not alternate in direction the way separate solids would.
fn pair_around_edge(points: &[Point3], mesh: &SurfaceMesh, e: EdgeId) -> Option<Vec<(u32, u32)>> {
    let [u, v] = mesh.edge_vertices(e)?;
    let (pu, pv) = (points[u as usize], points[v as usize]);
    let mut items: Vec<(u32, bool, Point3)> = Vec::new();
    for &h in mesh.edge_half_edges(e) {
        let f = h / 3;
        let (from, _) = mesh.endpoints(HalfEdge(h))?;
        let tri = mesh.triangle(f)?;
        let w = tri.into_iter().find(|&x| x != u && x != v)?;
        items.push((f, from == u, points[w as usize]));
    }
    let n = items.len();
    if n % 2 == 1 {
        return None;
    }
    let pr = items[0].2;
    let axis = plane_axis(pu, pv, pr)?;
    let half = |pw: Point3| match orient3d(pu, pv, pr, pw) {
        Sign::Positive => 0,
        Sign::Negative => 1,
        Sign::Zero => {
            if orient2d_projected(pu, pv, pr, axis) == orient2d_projected(pu, pv, pw, axis) {
                0
            } else {
                1
            }
        }
    };
    items.sort_by(|a, b| {
        half(a.2)
            .cmp(&half(b.2))
            .then_with(|| match orient3d(pu, pv, a.2, b.2) {
                Sign::Positive => Ordering::Less,
                Sign::Negative => Ordering::Greater,
                Sign::Zero => Ordering::Equal,
            })
            .then(a.0.cmp(&b.0))
    });
    for offset in 0..2 {
        let ok =
            (0..n / 2).all(|i| !items[(offset + 2 * i) % n].1 && items[(offset + 2 * i + 1) % n].1);
        if ok {
            return Some(
                (0..n / 2)
                    .map(|i| {
                        (
                            items[(offset + 2 * i) % n].0,
                            items[(offset + 2 * i + 1) % n].0,
                        )
                    })
                    .collect(),
            );
        }
    }
    None
}
