//! Half-edge connectivity for triangle soups that may be non-manifold,
//! inconsistently oriented or degenerate.
//!
//! Half-edge `3t + i` runs from corner `i` to corner `i+1` of triangle `t`, so
//! `next`, `prev` and the owning triangle are arithmetic. Undirected edges are
//! stored separately with every half-edge on them, which is what makes
//! non-manifold input representable: an edge simply has one, two, or more
//! half-edges. A `twin` exists only when an edge has exactly two.
//!
//! Everything is built by sorting, never by hash-map iteration, so all
//! numbering is deterministic.

use crate::soup::TriangleSoup;

/// A half-edge: `3 × triangle + corner`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HalfEdge(pub u32);

/// An undirected edge, numbered in order of its sorted vertex pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeId(pub u32);

impl HalfEdge {
    #[must_use]
    pub fn triangle(self) -> u32 {
        self.0 / 3
    }

    #[must_use]
    pub fn corner(self) -> u32 {
        self.0 % 3
    }

    #[must_use]
    pub fn next(self) -> HalfEdge {
        HalfEdge(self.0 - self.corner() + (self.corner() + 1) % 3)
    }

    #[must_use]
    pub fn prev(self) -> HalfEdge {
        HalfEdge(self.0 - self.corner() + (self.corner() + 2) % 3)
    }
}

/// How the faces around an edge fit together.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EdgeKind {
    /// One face: a hole's rim.
    Boundary,
    /// Two faces traversing it in opposite directions.
    Manifold,
    /// Two faces traversing it in the same direction.
    Inconsistent,
    /// Three or more faces.
    NonManifold,
    /// Both ends are the same vertex (from a face with a repeated index).
    Degenerate,
}

#[derive(Debug, Clone)]
pub struct SurfaceMesh {
    triangles: Vec<[u32; 3]>,
    vertex_count: u32,
    /// Sorted, unique `[min, max]` vertex pairs.
    edges: Vec<[u32; 2]>,
    /// Edge of each half-edge.
    edge_of: Vec<u32>,
    /// CSR: half-edges of edge `e` are `edge_half_edges[edge_start[e]..edge_start[e+1]]`, ascending.
    edge_start: Vec<u32>,
    edge_half_edges: Vec<u32>,
    /// CSR: outgoing half-edges of vertex `v`, ascending.
    vertex_start: Vec<u32>,
    vertex_half_edges: Vec<u32>,
}

impl SurfaceMesh {
    #[must_use]
    pub fn new(soup: &TriangleSoup<'_>) -> SurfaceMesh {
        let triangles = soup.triangles().to_vec();
        // TriangleSoup guarantees both counts fit in u32 (and 3 × triangles too).
        let vertex_count = soup.points().len() as u32;
        let half_edges = triangles.len() * 3;

        let mut keyed: Vec<([u32; 2], u32)> = (0..half_edges as u32)
            .map(|h| {
                let (a, b) = endpoints(&triangles, HalfEdge(h));
                ([a.min(b), a.max(b)], h)
            })
            .collect();
        keyed.sort_unstable();

        let mut edges = Vec::new();
        let mut edge_of = vec![0u32; half_edges];
        let mut edge_start = vec![0u32];
        let mut edge_half_edges = Vec::with_capacity(half_edges);
        for (key, h) in keyed {
            if edges.last() != Some(&key) {
                if !edges.is_empty() {
                    edge_start.push(edge_half_edges.len() as u32);
                }
                edges.push(key);
            }
            edge_of[h as usize] = (edges.len() - 1) as u32;
            edge_half_edges.push(h);
        }
        edge_start.push(edge_half_edges.len() as u32);
        if edges.is_empty() {
            edge_start = vec![0];
        }

        let mut by_origin: Vec<(u32, u32)> = (0..half_edges as u32)
            .map(|h| (endpoints(&triangles, HalfEdge(h)).0, h))
            .collect();
        by_origin.sort_unstable();
        let mut vertex_start = vec![0u32; vertex_count as usize + 1];
        for &(v, _) in &by_origin {
            vertex_start[v as usize + 1] += 1;
        }
        for v in 0..vertex_count as usize {
            vertex_start[v + 1] += vertex_start[v];
        }
        let vertex_half_edges = by_origin.into_iter().map(|(_, h)| h).collect();

        SurfaceMesh {
            triangles,
            vertex_count,
            edges,
            edge_of,
            edge_start,
            edge_half_edges,
            vertex_start,
            vertex_half_edges,
        }
    }

    #[must_use]
    pub fn vertex_count(&self) -> u32 {
        self.vertex_count
    }

    #[must_use]
    pub fn triangle_count(&self) -> u32 {
        self.triangles.len() as u32
    }

    #[must_use]
    pub fn edge_count(&self) -> u32 {
        self.edges.len() as u32
    }

    /// Vertices of triangle `t`, or `None` if out of range.
    #[must_use]
    pub fn triangle(&self, t: u32) -> Option<[u32; 3]> {
        self.triangles.get(t as usize).copied()
    }

    /// Every half-edge, in order.
    pub fn half_edges(&self) -> impl Iterator<Item = HalfEdge> + '_ {
        (0..self.triangles.len() as u32 * 3).map(HalfEdge)
    }

    /// Every edge, in order.
    pub fn edges(&self) -> impl Iterator<Item = EdgeId> + '_ {
        (0..self.edges.len() as u32).map(EdgeId)
    }

    /// `(origin, target)` of a half-edge, or `None` if out of range.
    #[must_use]
    pub fn endpoints(&self, h: HalfEdge) -> Option<(u32, u32)> {
        ((h.0 as usize) < self.edge_of.len()).then(|| endpoints(&self.triangles, h))
    }

    /// The undirected edge a half-edge lies on.
    #[must_use]
    pub fn edge(&self, h: HalfEdge) -> Option<EdgeId> {
        self.edge_of.get(h.0 as usize).map(|&e| EdgeId(e))
    }

    /// `[min, max]` vertices of an edge.
    #[must_use]
    pub fn edge_vertices(&self, e: EdgeId) -> Option<[u32; 2]> {
        self.edges.get(e.0 as usize).copied()
    }

    /// The half-edges on an edge, ascending; empty if `e` is out of range.
    #[must_use]
    pub fn edge_half_edges(&self, e: EdgeId) -> &[u32] {
        let i = e.0 as usize;
        match (self.edge_start.get(i), self.edge_start.get(i + 1)) {
            (Some(&s), Some(&t)) => &self.edge_half_edges[s as usize..t as usize],
            _ => &[],
        }
    }

    /// The other half-edge on this half-edge's edge, if the edge has exactly two.
    #[must_use]
    pub fn twin(&self, h: HalfEdge) -> Option<HalfEdge> {
        match *self.edge_half_edges(self.edge(h)?) {
            [a, b] => Some(HalfEdge(if a == h.0 { b } else { a })),
            _ => None,
        }
    }

    #[must_use]
    pub fn edge_kind(&self, e: EdgeId) -> Option<EdgeKind> {
        let [a, b] = self.edge_vertices(e)?;
        if a == b {
            return Some(EdgeKind::Degenerate);
        }
        Some(match *self.edge_half_edges(e) {
            [_] => EdgeKind::Boundary,
            [h, k] => {
                let origin = |x: u32| endpoints(&self.triangles, HalfEdge(x)).0;
                if origin(h) == origin(k) {
                    EdgeKind::Inconsistent
                } else {
                    EdgeKind::Manifold
                }
            }
            _ => EdgeKind::NonManifold,
        })
    }

    /// Half-edges leaving vertex `v`, ascending; empty if out of range.
    #[must_use]
    pub fn outgoing(&self, v: u32) -> &[u32] {
        let i = v as usize;
        match (self.vertex_start.get(i), self.vertex_start.get(i + 1)) {
            (Some(&s), Some(&t)) => &self.vertex_half_edges[s as usize..t as usize],
            _ => &[],
        }
    }

    /// Triangles containing vertex `v`, ascending and without repeats.
    #[must_use]
    pub fn vertex_triangles(&self, v: u32) -> Vec<u32> {
        let mut t: Vec<u32> = self.outgoing(v).iter().map(|&h| h / 3).collect();
        t.dedup();
        t
    }

    /// Number of fans at `v`: groups of incident triangles connected through
    /// edges at `v`. A manifold vertex has one (zero if isolated); more than
    /// one means separate surface sheets meet at a single point.
    #[must_use]
    pub fn vertex_fan_count(&self, v: u32) -> usize {
        let tris = self.vertex_triangles(v);
        let mut parent: Vec<usize> = (0..tris.len()).collect();
        let local = |t: u32| tris.binary_search(&t).ok();
        for &h in self.outgoing(v) {
            // Both edges of the corner at v: the outgoing one and the incoming one.
            for edge_half in [HalfEdge(h), HalfEdge(h).prev()] {
                let Some(e) = self.edge(edge_half) else {
                    continue;
                };
                if self.edge_kind(e) == Some(EdgeKind::Degenerate) {
                    continue;
                }
                let mut members = self.edge_half_edges(e).iter().filter_map(|&k| local(k / 3));
                if let Some(first) = members.next() {
                    for m in members {
                        union(&mut parent, first, m);
                    }
                }
            }
        }
        (0..tris.len())
            .filter(|&i| find(&mut parent, i) == i)
            .count()
    }
}

fn endpoints(triangles: &[[u32; 3]], h: HalfEdge) -> (u32, u32) {
    let tri = triangles[h.triangle() as usize];
    let i = h.corner() as usize;
    (tri[i], tri[(i + 1) % 3])
}

fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

fn union(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (find(parent, a), find(parent, b));
    if ra != rb {
        parent[ra.max(rb)] = ra.min(rb);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mesh(points: usize, tris: &[[u32; 3]]) -> SurfaceMesh {
        let pts = vec![[0.0; 3]; points];
        SurfaceMesh::new(&TriangleSoup::new(&pts, tris).unwrap())
    }

    #[test]
    fn tetrahedron_is_closed_and_manifold() {
        let m = mesh(4, &[[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]);
        assert_eq!(m.edge_count(), 6);
        assert!(
            m.edges()
                .all(|e| m.edge_kind(e) == Some(EdgeKind::Manifold))
        );
        for h in m.half_edges() {
            let t = m.twin(h).unwrap();
            assert_eq!(m.twin(t), Some(h));
            let (a, b) = m.endpoints(h).unwrap();
            assert_eq!(m.endpoints(t), Some((b, a)));
            assert_eq!(h.next().next().next(), h);
            assert_eq!(h.next().prev(), h);
        }
        assert!((0..4).all(|v| m.vertex_fan_count(v) == 1));
    }

    #[test]
    fn edge_kinds() {
        // 0: boundary pair, 1: inconsistent pair, 2: three faces on edge (4,5).
        let m = mesh(
            10,
            &[
                [0, 1, 2],
                [3, 4, 2],
                [3, 4, 1],
                [4, 5, 6],
                [5, 4, 7],
                [4, 5, 8],
                [9, 9, 0],
            ],
        );
        let kind = |a, b| {
            m.edges()
                .find(|&e| m.edge_vertices(e) == Some([a, b]))
                .and_then(|e| m.edge_kind(e))
        };
        assert_eq!(kind(0, 1), Some(EdgeKind::Boundary));
        assert_eq!(kind(3, 4), Some(EdgeKind::Inconsistent));
        assert_eq!(kind(4, 5), Some(EdgeKind::NonManifold));
        assert_eq!(kind(9, 9), Some(EdgeKind::Degenerate));
    }

    #[test]
    fn bowtie_vertex_has_two_fans() {
        let m = mesh(5, &[[0, 1, 2], [0, 3, 4]]);
        assert_eq!(m.vertex_fan_count(0), 2);
        assert_eq!(m.vertex_fan_count(1), 1);
    }

    #[test]
    fn empty_and_out_of_range_queries() {
        let m = mesh(3, &[]);
        assert_eq!(m.edge_count(), 0);
        assert!(m.outgoing(7).is_empty());
        assert!(m.edge_half_edges(EdgeId(9)).is_empty());
        assert_eq!(m.twin(HalfEdge(0)), None);
        assert_eq!(m.vertex_fan_count(1), 0);
    }
}
