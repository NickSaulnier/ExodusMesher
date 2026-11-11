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

/// A vertex of the chamfered surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChamferVertex {
    /// An input vertex, kept.
    Input(u32),
    /// On input edge `edge` (lower vertex first), cut at distance `epsilon`
    /// from its end `from`, an acute vertex.
    OnEdge {
        edge: (u32, u32),
        from: u32,
        epsilon: f64,
    },
    /// Inside input triangle `triangle`: a bridge point.
    InFace { triangle: u32 },
}

/// One face of the chamfered surface: what is left of an input triangle,
/// a polygon in the triangle's orientation.
#[derive(Debug, Clone, PartialEq)]
pub struct ChamferedFace {
    pub triangle: u32,
    /// Vertex indices, in order around the face.
    pub polygon: Vec<u32>,
}

/// The chamfered surface: every acute vertex and edge cut away, so no two
/// of its constraints sharing a vertex meet at an acute angle. It is open:
/// the cuts leave holes, bounded by the bridges.
#[derive(Debug, Clone)]
pub struct Chamfered {
    pub vertices: Vec<ChamferVertex>,
    /// Each vertex's exact position: input vertices explicit, cuts on edges
    /// LNC points, bridge points BPT points of their triangle.
    pub points: Vec<Point>,
    pub faces: Vec<ChamferedFace>,
    /// How far each input vertex was chamfered (zero if it was not).
    pub epsilon: Vec<f64>,
}

fn sub3(a: Point3, b: Point3) -> Point3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot3(a: Point3, b: Point3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm3(a: Point3) -> f64 {
    dot3(a, a).sqrt()
}

fn scale3(a: Point3, s: f64) -> Point3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn add3(a: Point3, b: Point3) -> Point3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn unit(a: Point3) -> Option<Point3> {
    let n = norm3(a);
    (n > 0.0 && n.is_finite()).then(|| scale3(a, 1.0 / n))
}

fn cross3(a: Point3, b: Point3) -> Point3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The unit vector in the plane of `u` and `w`, perpendicular to `u`, on
/// `w`'s side: `(u × w) × u`, by cross products, which stay accurate for a
/// needle (nearly parallel `u` and `w`) where subtracting projections
/// cancels.
fn perpendicular(u: Point3, w: Point3) -> Option<Point3> {
    unit(cross3(cross3(u, w), u))
}

/// The length `L` of a bridge's three edges: the legs from the cuts `U`,
/// `V` (at distance `e` along unit `a`, `b` from the corner) along the unit
/// perpendiculars `na`, `nb`, and the edge between their ends, all equal:
/// `|e(a − b) + L(na − nb)| = L`, the smaller positive root.
fn bridge_length(e: f64, a: Point3, b: Point3, na: Point3, nb: Point3) -> Option<f64> {
    let d = scale3(sub3(a, b), e);
    let f = sub3(na, nb);
    let (qa, qb, qc) = (dot3(f, f) - 1.0, 2.0 * dot3(d, f), dot3(d, d));
    let roots: Vec<f64> = if qa.abs() < 1e-12 {
        if qb == 0.0 {
            Vec::new()
        } else {
            vec![-qc / qb]
        }
    } else {
        let disc = qb * qb - 4.0 * qa * qc;
        if disc < 0.0 {
            Vec::new()
        } else {
            let s = disc.sqrt();
            vec![(-qb - s) / (2.0 * qa), (-qb + s) / (2.0 * qa)]
        }
    };
    roots
        .into_iter()
        .filter(|l| *l > 0.0 && l.is_finite())
        .min_by(f64::total_cmp)
}

/// Barycentric weights of `p` in triangle `t` (assumed in its plane), kept
/// inside it: negative weights are clamped and the rest renormalized.
fn barycentric(p: Point3, t: [Point3; 3]) -> [f64; 3] {
    let n = cross3(sub3(t[1], t[0]), sub3(t[2], t[0]));
    let area = |a: Point3, b: Point3, c: Point3| dot3(n, cross3(sub3(b, a), sub3(c, a)));
    let whole = area(t[0], t[1], t[2]);
    let w = [
        area(p, t[1], t[2]) / whole,
        area(t[0], p, t[2]) / whole,
        area(t[0], t[1], p) / whole,
    ]
    .map(|x| if x.is_finite() { x.max(0.0) } else { 0.0 });
    let sum: f64 = w.iter().sum();
    if sum > 0.0 {
        w.map(|x| x / sum)
    } else {
        [1.0 / 3.0; 3]
    }
}

/// A bridge leg: from a cut on an input edge (at `base`), the unit direction
/// into the triangle, its length, and the bridge point at its end (a vertex
/// index).
#[derive(Debug, Clone, Copy)]
struct Leg {
    base: Point3,
    direction: Point3,
    length: f64,
    end: u32,
}

/// Chamfer a classified surface (see [`classify`]): every acute vertex `P`
/// is cut at `ε = lfs(P) / 3` on its edges, each triangle at `P` taking a
/// bridge of three equal edges whose legs are perpendicular to its edges
/// (so no new acute angle appears), and every acute edge is then cut along
/// its length: in each triangle at it, the bridge legs are shortened to the
/// shortest one and the rectangle between the edge and the legs' ends
/// removed.
///
/// # Errors
/// [`ChamferError::Degenerate`] where a bridge cannot be built (a triangle
/// with collinear corners).
pub fn chamfer(
    points: &[Point3],
    triangles: &[[u32; 3]],
    acute: &Acute,
) -> Result<Chamfered, ChamferError> {
    let n = points.len();
    let epsilon: Vec<f64> = (0..n)
        .map(|v| {
            if acute.vertices[v] {
                acute.vertex_lfs[v] / 3.0
            } else {
                0.0
            }
        })
        .collect();
    let mut vertices: Vec<ChamferVertex> = (0..n as u32).map(ChamferVertex::Input).collect();
    let mut exact: Vec<Point> = points
        .iter()
        .enumerate()
        .map(|(i, &p)| Point::explicit(p).map_err(|_| ChamferError::NonFinite { vertex: i as u32 }))
        .collect::<Result<_, _>>()?;
    let mut approx: Vec<Point3> = points.to_vec();
    // Cuts on edges, shared by the two triangles at the edge.
    let mut cuts: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    let mut cut = |from: u32,
                   to: u32,
                   vertices: &mut Vec<ChamferVertex>,
                   exact: &mut Vec<Point>,
                   approx: &mut Vec<Point3>|
     -> Result<u32, ChamferError> {
        if let Some(&c) = cuts.get(&(from, to)) {
            return Ok(c);
        }
        let (p, x) = (points[from as usize], points[to as usize]);
        let length = norm3(sub3(x, p));
        let e = epsilon[from as usize];
        let t = 1.0 - e / length;
        let point = Point::lnc(p, x, t).map_err(|_| ChamferError::NonFinite { vertex: from })?;
        let id = vertices.len() as u32;
        vertices.push(ChamferVertex::OnEdge {
            edge: (from.min(to), from.max(to)),
            from,
            epsilon: e,
        });
        approx.push(point.approx());
        exact.push(point);
        cuts.insert((from, to), id);
        Ok(id)
    };
    let mut faces = Vec::with_capacity(triangles.len());
    for (ti, t) in triangles.iter().enumerate() {
        let tri = t.map(|v| points[v as usize]);
        let mut polygon: Vec<u32> = Vec::new();
        let mut legs: Vec<Leg> = Vec::new();
        for k in 0..3 {
            let (p, next, prev) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            if !acute.vertices[p as usize] {
                polygon.push(p);
                continue;
            }
            // Corner p, between prev and next: replaced by the cut on p-prev,
            // its bridge point, the other bridge point, the cut on p-next.
            let at = points[p as usize];
            let (a, b) = (
                sub3(points[next as usize], at),
                sub3(points[prev as usize], at),
            );
            let degenerate = ChamferError::Degenerate {
                triangle: ti as u32,
            };
            let (ua, ub) = (
                unit(a).ok_or(degenerate.clone())?,
                unit(b).ok_or(degenerate.clone())?,
            );
            let na = perpendicular(a, b).ok_or(degenerate.clone())?;
            let nb = perpendicular(b, a).ok_or(degenerate.clone())?;
            let e = epsilon[p as usize];
            let length = bridge_length(e, ua, ub, na, nb).ok_or(degenerate)?;
            let u = cut(p, next, &mut vertices, &mut exact, &mut approx)?;
            let v = cut(p, prev, &mut vertices, &mut exact, &mut approx)?;
            let (u_base, v_base) = (approx[u as usize], approx[v as usize]);
            let mut leg = |base: Point3, direction: Point3| {
                let id = vertices.len() as u32;
                vertices.push(ChamferVertex::InFace {
                    triangle: ti as u32,
                });
                exact.push(Point::from(base));
                approx.push(base);
                legs.push(Leg {
                    base,
                    direction,
                    length,
                    end: id,
                });
                id
            };
            let v_end = leg(v_base, nb);
            let u_end = leg(u_base, na);
            polygon.extend([v, v_end, u_end, u]);
        }
        // Edge chamfers: legs shortened to the shortest, then each acute
        // edge's cuts dropped from the face (the rectangle between them and
        // their legs' ends goes).
        let has_acute_edge = (0..3).any(|k| {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            acute
                .edges
                .get(&(u.min(v), u.max(v)))
                .is_some_and(|e| e.acute)
        });
        if has_acute_edge {
            let shortest = legs.iter().map(|l| l.length).fold(f64::INFINITY, f64::min);
            for l in &mut legs {
                l.length = shortest;
            }
            let on_acute_edge = |c: u32| match vertices[c as usize] {
                ChamferVertex::OnEdge { edge, .. } => {
                    acute.edges.get(&edge).is_some_and(|e| e.acute)
                }
                _ => false,
            };
            polygon.retain(|&c| !on_acute_edge(c));
        }
        // Bridge points: the legs' ends, exactly on the triangle.
        for l in &legs {
            let end = add3(l.base, scale3(l.direction, l.length));
            let w = barycentric(end, tri);
            let p = Point::bpt(tri, w).map_err(|_| ChamferError::Degenerate {
                triangle: ti as u32,
            })?;
            approx[l.end as usize] = p.approx();
            exact[l.end as usize] = p;
        }
        faces.push(ChamferedFace {
            triangle: ti as u32,
            polygon,
        });
    }
    Ok(Chamfered {
        vertices,
        points: exact,
        faces,
        epsilon,
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

    /// Chamfer and check: every face vertex on its own triangle, every face
    /// polygon a proper polygon covering less than its triangle, and at every
    /// new vertex the two face edges there not acute (to rounding: bridge
    /// points are rounded onto their triangle). Returns the faces' vertex
    /// counts.
    fn check_chamfer(name: &str, p: &[Point3], f: &[[u32; 3]]) -> Vec<usize> {
        let a = classify(p, f).unwrap();
        let c = chamfer(p, f, &a).unwrap();
        let at = |v: u32| c.points[v as usize].approx();
        for face in &c.faces {
            let t = f[face.triangle as usize];
            assert!(face.polygon.len() >= 3, "{name}: {face:?}");
            for &v in &face.polygon {
                match c.vertices[v as usize] {
                    ChamferVertex::Input(x) => assert!(t.contains(&x), "{name}"),
                    ChamferVertex::OnEdge { edge, .. } => {
                        assert!(t.contains(&edge.0) && t.contains(&edge.1), "{name}");
                    }
                    ChamferVertex::InFace { triangle } => {
                        assert_eq!(triangle, face.triangle, "{name}");
                    }
                }
            }
            let k = face.polygon.len();
            for i in 0..k {
                let v = face.polygon[i];
                if matches!(c.vertices[v as usize], ChamferVertex::Input(_)) {
                    continue;
                }
                let prev = at(face.polygon[(i + k - 1) % k]);
                let next = at(face.polygon[(i + 1) % k]);
                let here = at(v);
                let (u, w) = (sub3(prev, here), sub3(next, here));
                let cos = dot3(u, w) / (norm3(u) * norm3(w));
                assert!(cos <= 1e-6, "{name}: acute at {v} ({cos}) in {face:?}");
            }
            // Area: positive, and under the triangle's.
            let tri = t.map(|v| p[v as usize]);
            let normal = {
                let u = sub3(tri[1], tri[0]);
                let w = sub3(tri[2], tri[0]);
                [
                    u[1] * w[2] - u[2] * w[1],
                    u[2] * w[0] - u[0] * w[2],
                    u[0] * w[1] - u[1] * w[0],
                ]
            };
            let mut area = 0.0;
            for i in 0..k {
                let (q, r) = (at(face.polygon[i]), at(face.polygon[(i + 1) % k]));
                let (q, r) = (sub3(q, tri[0]), sub3(r, tri[0]));
                area += dot3(
                    normal,
                    [
                        q[1] * r[2] - q[2] * r[1],
                        q[2] * r[0] - q[0] * r[2],
                        q[0] * r[1] - q[1] * r[0],
                    ],
                );
            }
            let whole = dot3(normal, normal);
            assert!(
                area > 0.0 && area <= whole * (1.0 + 1e-12),
                "{name}: area {area} of {whole}"
            );
        }
        c.faces.iter().map(|f| f.polygon.len()).collect()
    }

    #[test]
    fn chamfering_cuts_acute_corners_and_edges() {
        // A cube: nothing to cut.
        let (p, f) = cube();
        assert!(check_chamfer("cube", &p, &f).iter().all(|&n| n == 3));
        // A regular tetrahedron: every corner and edge cut, leaving each face
        // a hexagon of bridge points.
        let p = vec![
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        let f = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        assert!(check_chamfer("tetrahedron", &p, &f).iter().all(|&n| n == 6));
        // A wedge prism: its lengthwise edges cut, its end caps' corners too.
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
        check_chamfer("wedge", &p, &f);
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
