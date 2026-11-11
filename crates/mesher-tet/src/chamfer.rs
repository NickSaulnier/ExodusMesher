//! Chamfering a closed triangle surface so Delaunay refinement is sure to
//! finish (Milestone 14, after Diazzi, Dai, Panozzo and Attene, "Surface
//! chamfering for robust tetrahedral meshing", SIGGRAPH 2026, from the
//! paper's text).
//!
//! Delaunay refinement terminates when the input satisfies the projection
//! condition: no two constraints sharing a vertex meet at an acute angle.
//! This module finds where the condition fails, exactly:
//!
//! - an edge is *acute* when its two faces meet at an acute dihedral angle
//!   (measured either side: refinement meshes both), and *flat* when they
//!   are coplanar, folding out to 180°;
//! - a vertex is *acute* when an incident edge is acute, two incident edges
//!   meet at an acute angle, or an incident edge leans over an incident face
//!   (its projection meets the face's interior). Flat edges are ignored
//!   here: they only triangulate a planar region.
//!
//! The local feature size of a vertex or an edge is its distance to its
//! link (the opposite edges of its triangles, or their opposite vertices):
//! approximate, as it only chooses how far chamfering cuts.

use std::collections::BTreeMap;

use mesher_geom::Point3;
use mesher_geom::Sign;
use mesher_geom::implicit::{self, Point};

/// Why a surface cannot be classified.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ChamferError {
    /// An edge with other than two triangles: the surface is not a closed
    /// 2-manifold.
    NotManifold { edge: (u32, u32) },
    /// A non-finite coordinate.
    NonFinite { vertex: u32 },
    /// A triangle with a repeated vertex.
    Degenerate { triangle: u32 },
}

impl std::fmt::Display for ChamferError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ChamferError::NotManifold { edge } => {
                write!(f, "edge {}-{} does not have two triangles", edge.0, edge.1)
            }
            ChamferError::NonFinite { vertex } => write!(f, "vertex {vertex} is not finite"),
            ChamferError::Degenerate { triangle } => {
                write!(f, "triangle {triangle} repeats a vertex")
            }
        }
    }
}

/// Where a closed surface breaks the projection condition.
#[derive(Debug, Clone, PartialEq)]
pub struct Acute {
    /// Every edge (lower vertex first), with its two triangles' opposite
    /// vertices, whether it is acute, and whether it is flat.
    pub edges: BTreeMap<(u32, u32), EdgeInfo>,
    /// Whether each vertex is acute.
    pub vertices: Vec<bool>,
    /// Local feature size of each vertex: its smallest triangle height at it.
    pub vertex_lfs: Vec<f64>,
}

/// One edge of the surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeInfo {
    /// The opposite vertices of its two triangles.
    pub opposite: [u32; 2],
    /// Its two triangles.
    pub triangles: [u32; 2],
    pub acute: bool,
    pub flat: bool,
    /// Local feature size: the distance from the edge to the nearer opposite
    /// vertex.
    pub lfs: f64,
}

impl Acute {
    /// Acute edges, lower vertex first.
    pub fn acute_edges(&self) -> impl Iterator<Item = (u32, u32)> + '_ {
        self.edges.iter().filter(|(_, e)| e.acute).map(|(&k, _)| k)
    }

    /// Acute vertices.
    pub fn acute_vertices(&self) -> impl Iterator<Item = u32> + '_ {
        self.vertices
            .iter()
            .enumerate()
            .filter(|(_, a)| **a)
            .map(|(v, _)| v as u32)
    }
}

fn distance_to_segment(p: Point3, a: Point3, b: Point3) -> f64 {
    let d = [0, 1, 2].map(|k| b[k] - a[k]);
    let w = [0, 1, 2].map(|k| p[k] - a[k]);
    let dd: f64 = d.iter().map(|x| x * x).sum();
    let t = if dd > 0.0 {
        ((0..3).map(|k| w[k] * d[k]).sum::<f64>() / dd).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (0..3)
        .map(|k| (p[k] - a[k] - t * d[k]).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Classify a closed, oriented, 2-manifold triangle surface. Every decision
/// is an exact predicate on the input coordinates.
///
/// # Errors
/// See [`ChamferError`].
pub fn classify(points: &[Point3], triangles: &[[u32; 3]]) -> Result<Acute, ChamferError> {
    let pts: Vec<Point> = points
        .iter()
        .enumerate()
        .map(|(i, &p)| Point::explicit(p).map_err(|_| ChamferError::NonFinite { vertex: i as u32 }))
        .collect::<Result<_, _>>()?;
    let mut sides: BTreeMap<(u32, u32), Vec<(u32, u32)>> = BTreeMap::new();
    for (ti, t) in triangles.iter().enumerate() {
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            return Err(ChamferError::Degenerate {
                triangle: ti as u32,
            });
        }
        for k in 0..3 {
            let (u, v, w) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            sides
                .entry((u.min(v), u.max(v)))
                .or_default()
                .push((ti as u32, w));
        }
    }
    let mut edges = BTreeMap::new();
    for (&(c, d), around) in &sides {
        let &[(t0, a), (t1, b)] = around.as_slice() else {
            return Err(ChamferError::NotManifold { edge: (c, d) });
        };
        let p = |v: u32| &pts[v as usize];
        let dihedral = implicit::dihedral_sign(p(a), p(b), p(c), p(d));
        let coplanar = implicit::orient3d(p(a), p(b), p(c), p(d)) == Sign::Zero;
        let at = |v: u32| points[v as usize];
        edges.insert(
            (c, d),
            EdgeInfo {
                opposite: [a, b],
                triangles: [t0, t1],
                acute: dihedral == Sign::Positive,
                flat: coplanar && dihedral == Sign::Negative,
                lfs: distance_to_segment(at(a), at(c), at(d)).min(distance_to_segment(
                    at(b),
                    at(c),
                    at(d),
                )),
            },
        );
    }
    // The non-flat edges and the triangles at each vertex.
    let n = points.len();
    let mut feature: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut at_vertex: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut vertices = vec![false; n];
    for (&(c, d), e) in &edges {
        if e.flat {
            continue;
        }
        feature[c as usize].push(d);
        feature[d as usize].push(c);
        if e.acute {
            vertices[c as usize] = true;
            vertices[d as usize] = true;
        }
    }
    for (ti, t) in triangles.iter().enumerate() {
        for &v in t {
            at_vertex[v as usize].push(ti as u32);
        }
    }
    let mut vertex_lfs = vec![f64::INFINITY; n];
    for (v, tris) in at_vertex.iter().enumerate() {
        for &ti in tris {
            let t = triangles[ti as usize];
            let k = t.iter().position(|&x| x as usize == v).unwrap_or(0);
            let (a, b) = (t[(k + 1) % 3], t[(k + 2) % 3]);
            let h = distance_to_segment(points[v], points[a as usize], points[b as usize]);
            vertex_lfs[v] = vertex_lfs[v].min(h);
        }
    }
    for v in 0..n {
        if vertices[v] {
            continue;
        }
        let p = |x: u32| &pts[x as usize];
        let here = &pts[v];
        let near = &feature[v];
        // Two incident edges at an acute angle.
        'pairs: for (i, &q1) in near.iter().enumerate() {
            for &q2 in &near[i + 1..] {
                if implicit::dot_sign(p(q1), here, p(q2)) == Sign::Positive {
                    vertices[v] = true;
                    break 'pairs;
                }
            }
        }
        if vertices[v] {
            continue;
        }
        // An incident edge leaning over an incident face.
        'lean: for &q in near {
            for &ti in &at_vertex[v] {
                let t = triangles[ti as usize];
                if t.contains(&q) {
                    continue;
                }
                let k = t.iter().position(|&x| x as usize == v).unwrap_or(0);
                let (b, c) = (t[(k + 1) % 3], t[(k + 2) % 3]);
                if implicit::edge_projects_into_triangle(here, p(q), p(b), p(c)) {
                    vertices[v] = true;
                    break 'lean;
                }
            }
        }
    }
    Ok(Acute {
        edges,
        vertices,
        vertex_lfs,
    })
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;

    fn cube() -> (Vec<Point3>, Vec<[u32; 3]>) {
        let p: Vec<Point3> = (0..8)
            .map(|i| {
                [
                    f64::from(i & 1),
                    f64::from((i >> 1) & 1),
                    f64::from((i >> 2) & 1),
                ]
            })
            .collect();
        let quads = [
            [0, 2, 3, 1],
            [4, 5, 7, 6],
            [0, 1, 5, 4],
            [2, 6, 7, 3],
            [0, 4, 6, 2],
            [1, 3, 7, 5],
        ];
        let f = quads
            .iter()
            .flat_map(|&[a, b, c, d]| [[a, b, c], [a, c, d]])
            .collect();
        (p, f)
    }

    /// A cube's edges are right angles and its diagonals flat: nothing is
    /// acute, and the local feature size is the shortest height.
    #[test]
    fn a_cube_has_nothing_acute() {
        let (p, f) = cube();
        let a = classify(&p, &f).unwrap();
        assert_eq!(a.acute_edges().count(), 0);
        assert_eq!(a.acute_vertices().count(), 0);
        assert_eq!(a.edges.values().filter(|e| e.flat).count(), 6);
        for &l in &a.vertex_lfs {
            assert!((l - libm::sqrt(0.5)).abs() < 1e-12 || (l - 1.0).abs() < 1e-12);
        }
    }

    /// A regular tetrahedron: every dihedral is about 70.5°, so every edge
    /// and vertex is acute.
    #[test]
    fn a_tetrahedron_is_acute_everywhere() {
        let p = vec![
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        let f = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        let a = classify(&p, &f).unwrap();
        assert_eq!(a.acute_edges().count(), 6);
        assert_eq!(a.acute_vertices().count(), 4);
    }

    /// A triangular prism whose cross-section has angles of 80.5°, 80.5° and
    /// 19°: the three lengthwise edges are acute (their faces meet at those
    /// angles), the end caps meet the sides at right angles, and the quads'
    /// diagonals are flat.
    #[test]
    fn a_wedge_has_an_acute_ridge() {
        // Cross-section: base from (0,0) to (1,0), ridge at (0.5, 3).
        let mut p = Vec::new();
        for z in [0.0, 1.0] {
            p.extend([[0.0, 0.0, z], [1.0, 0.0, z], [0.5, 3.0, z]]);
        }
        let f = vec![
            [0, 2, 1],
            [3, 4, 5],
            [0, 1, 4],
            [0, 4, 3],
            [1, 2, 5],
            [1, 5, 4],
            [2, 0, 3],
            [2, 3, 5],
        ];
        let a = classify(&p, &f).unwrap();
        let acute: Vec<(u32, u32)> = a.acute_edges().collect();
        assert_eq!(acute, vec![(0, 3), (1, 4), (2, 5)]);
        assert_eq!(a.edges.values().filter(|e| e.flat).count(), 3);
        assert_eq!(a.acute_vertices().count(), 6);
    }

    /// A square pyramid whose apex sits right above a base corner: the
    /// apex's edges meet at acute angles, and at the far corner the edge up
    /// to the apex makes an acute angle with both base edges there.
    #[test]
    fn a_leaning_apex_is_acute() {
        let p = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
        ];
        let f = vec![
            [0, 2, 1],
            [0, 3, 2],
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
        ];
        let a = classify(&p, &f).unwrap();
        assert!(a.vertices[4], "the apex");
        assert!(a.vertices[2], "the far corner");
    }

    #[test]
    fn open_surfaces_are_refused() {
        let p = vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(matches!(
            classify(&p, &[[0, 1, 2]]),
            Err(ChamferError::NotManifold { .. })
        ));
    }
}
