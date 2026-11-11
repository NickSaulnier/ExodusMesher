//! Conforming Delaunay refinement of a chamfered surface (Milestone 14,
//! after Diazzi, Dai, Panozzo and Attene, "Surface chamfering for robust
//! tetrahedral meshing", SIGGRAPH 2026, from the paper's text; the
//! refinement itself is Ruppert's and Shewchuk's).
//!
//! The chamfered surface (see [`crate::chamfer`]) has no acute angle, so
//! Delaunay refinement of it is sure to finish. Its vertices and the corners
//! of an outer box are tetrahedralized, and the tetrahedralization stays
//! Delaunay: constraints are recovered by splitting them, never by flips.
//! This first stage recovers the segments: every face edge of the chamfered
//! surface that is a segment, and the box's edges. A subsegment is split at
//! its middle while it is not an edge or a vertex lies strictly inside its
//! diametral sphere (it is *encroached*: exactly, the angle the subsegment
//! subtends there is obtuse). Split points stay exact: on a segment along an
//! input edge (or a box edge) an LNC point of the edge; on a bridge segment
//! inside an input triangle a BPT point of the triangle, at the average of
//! its ends' barycentric weights (exactly on the triangle, not necessarily
//! on the segment, which the paper accepts: the chamfered surface is only an
//! intermediate, and the enriched CDT conforms to the input exactly).

use std::collections::BTreeMap;

use mesher_geom::implicit::{self, Kind, Point};
use mesher_geom::{Point3, Sign};

use crate::cdt::{DelaunayError, ExactDelaunay, Inserted};
use crate::chamfer::{ChamferVertex, FilteredChamfer};

/// Why refinement failed.
#[derive(Debug, Clone, PartialEq)]
pub enum RefineError {
    Delaunay(DelaunayError),
    /// A split point equals an existing vertex.
    Duplicate {
        point: Point3,
    },
    /// A split point could not be constructed (a degenerate carrier).
    Construction,
    /// More points than the budget.
    Budget {
        points: usize,
        limit: usize,
    },
}

impl std::fmt::Display for RefineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RefineError::Delaunay(e) => write!(f, "{e}"),
            RefineError::Duplicate { point } => {
                write!(f, "a split point equals a vertex near {point:?}")
            }
            RefineError::Construction => write!(f, "a split point could not be constructed"),
            RefineError::Budget { points, limit } => {
                write!(f, "{points} points exceed the budget of {limit}")
            }
        }
    }
}

impl From<DelaunayError> for RefineError {
    fn from(e: DelaunayError) -> RefineError {
        RefineError::Delaunay(e)
    }
}

/// What a segment lies on, to split it exactly.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Carrier {
    /// The line `a + s(b − a)`, `a` and `b` explicit.
    Line { a: Point3, b: Point3 },
    /// The plane of triangle `tri`, by barycentric weights.
    Triangle { tri: [Point3; 3] },
}

/// Where a segment end sits on its carrier.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Param {
    S(f64),
    W([f64; 3]),
}

/// A subsegment: its ends (vertices of the tetrahedralization), carrier, and
/// the ends' parameters on it.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Seg {
    ends: [u32; 2],
    carrier: Carrier,
    params: [Param; 2],
}

/// The chamfered surface's vertices and the outer box tetrahedralized, every
/// segment a chain of Delaunay edges.
#[derive(Debug, Clone)]
pub struct RefinedSegments {
    pub dt: ExactDelaunay,
    /// Every subsegment, by its end vertices.
    pub subsegments: Vec<[u32; 2]>,
    /// The tetrahedralization's vertex of each chamfered-surface vertex in
    /// use (`None` for one no face uses).
    pub vertex_of: Vec<Option<u32>>,
    /// The box corners' vertices.
    pub box_corners: [u32; 8],
    /// Points inserted on segments.
    pub steiner: usize,
}

/// The parameter of point `p` on line `a b` (`p = a + s(b − a)`), from its
/// construction: an end, or an LNC point of the line's ends in either
/// order. `None` if `p` is not built on that line.
fn line_param(p: &Point, a: Point3, b: Point3) -> Option<f64> {
    match *p.kind() {
        Kind::Explicit(q) if q == a => Some(0.0),
        Kind::Explicit(q) if q == b => Some(1.0),
        // t·x + (1 − t)·y
        Kind::Lnc { a: x, b: y, t } if x == b && y == a => Some(t),
        Kind::Lnc { a: x, b: y, t } if x == a && y == b => Some(1.0 - t),
        _ => None,
    }
}

/// The normalized barycentric weights of `p` in triangle `tri`, from its
/// construction: a corner, an LNC point of two corners, or a BPT point of
/// the triangle.
fn tri_param(p: &Point, tri: [Point3; 3]) -> Option<[f64; 3]> {
    let corner = |q: Point3| tri.iter().position(|&c| c == q);
    match *p.kind() {
        Kind::Explicit(q) => {
            let k = corner(q)?;
            let mut w = [0.0; 3];
            w[k] = 1.0;
            Some(w)
        }
        Kind::Lnc { a, b, t } => {
            let (i, j) = (corner(a)?, corner(b)?);
            let mut w = [0.0; 3];
            w[i] = t;
            w[j] = 1.0 - t;
            Some(w)
        }
        Kind::Bpt { tri: t, weights } if t == tri => {
            let sum: f64 = weights.iter().sum();
            (sum != 0.0 && sum.is_finite()).then(|| weights.map(|x| x / sum))
        }
        _ => None,
    }
}

impl Seg {
    /// The middle, exactly on the carrier.
    fn middle(&self) -> Option<(Point, Param)> {
        match (self.carrier, self.params) {
            (Carrier::Line { a, b }, [Param::S(s0), Param::S(s1)]) => {
                let s = 0.5 * (s0 + s1);
                // a + s(b − a) = s·b + (1 − s)·a
                Point::lnc(b, a, s).ok().map(|p| (p, Param::S(s)))
            }
            (Carrier::Triangle { tri }, [Param::W(w0), Param::W(w1)]) => {
                let w = [0, 1, 2].map(|k| 0.5 * (w0[k] + w1[k]));
                Point::bpt(tri, w).ok().map(|p| (p, Param::W(w)))
            }
            _ => None,
        }
    }
}

/// Is `p` strictly inside the diametral sphere of `a b` (is the angle
/// `a p b` obtuse)? Exact.
fn encroaches(a: &Point, p: &Point, b: &Point) -> bool {
    implicit::dot_sign(a, p, b) == Sign::Negative
}

/// Is subsegment `a b` a Delaunay edge that no vertex of the tets around it
/// encroaches?
fn segment_ok(dt: &ExactDelaunay, [a, b]: [u32; 2]) -> bool {
    let Some(ring) = dt.edge_ring(a, b) else {
        return false;
    };
    let (pa, pb) = (&dt.points()[a as usize], &dt.points()[b as usize]);
    !ring.iter().any(|&t| {
        dt.tet(t).is_some_and(|v| {
            v.iter()
                .any(|&w| w != a && w != b && encroaches(pa, &dt.points()[w as usize], pb))
        })
    })
}

/// Tetrahedralize a filtered chamfered surface (see
/// [`crate::chamfer::chamfer_filtered`]) and the corners of a box around it
/// (its bounding box grown by `margin` of its diagonal on every side), and
/// recover every segment and box edge as Delaunay edges by splitting what is
/// missing or encroached. `points` and `triangles` are the input the
/// surface was chamfered from; `max_points` bounds the vertex count.
///
/// # Errors
/// See [`RefineError`].
pub fn recover_segments_conforming(
    surface: &FilteredChamfer,
    points: &[Point3],
    triangles: &[[u32; 3]],
    margin: f64,
    max_points: usize,
) -> Result<RefinedSegments, RefineError> {
    let c = &surface.chamfered;
    // The vertices in use, then the box corners.
    let mut used = vec![false; c.vertices.len()];
    for f in &c.faces {
        for &v in &f.polygon {
            used[v as usize] = true;
        }
    }
    let mut vertex_of: Vec<Option<u32>> = vec![None; c.vertices.len()];
    let mut pts: Vec<Point> = Vec::new();
    for (v, &u) in used.iter().enumerate() {
        if u {
            vertex_of[v] = Some(pts.len() as u32);
            pts.push(c.points[v]);
        }
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &pts {
        let a = p.approx();
        for k in 0..3 {
            lo[k] = lo[k].min(a[k]);
            hi[k] = hi[k].max(a[k]);
        }
    }
    if pts.is_empty() {
        lo = [0.0; 3];
        hi = [1.0; 3];
    }
    let diagonal = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
    let grow = margin * diagonal.max(f64::MIN_POSITIVE);
    let (lo, hi) = (lo.map(|x| x - grow), hi.map(|x| x + grow));
    let corner = |i: usize| -> Point3 {
        [
            if i & 1 == 0 { lo[0] } else { hi[0] },
            if i & 2 == 0 { lo[1] } else { hi[1] },
            if i & 4 == 0 { lo[2] } else { hi[2] },
        ]
    };
    let first_corner = pts.len() as u32;
    for i in 0..8 {
        pts.push(Point::explicit(corner(i)).map_err(|_| RefineError::Construction)?);
    }
    let box_corners: [u32; 8] = [0, 1, 2, 3, 4, 5, 6, 7].map(|i| first_corner + i);
    let mut dt = ExactDelaunay::new(pts, 5)?;
    if let Some(&(v, _)) = dt.duplicates().first() {
        return Err(RefineError::Duplicate {
            point: dt.points()[v as usize].approx(),
        });
    }
    // Segments: face edges that are segments, by carrier.
    let mut segs: Vec<Seg> = Vec::new();
    let mut seen: BTreeMap<(u32, u32), ()> = BTreeMap::new();
    for f in &c.faces {
        let k = f.polygon.len();
        let tri_index = surface.kept_ids[f.triangle as usize];
        let corners = triangles[tri_index as usize].map(|v| points[v as usize]);
        for i in 0..k {
            let (u, v) = (f.polygon[i], f.polygon[(i + 1) % k]);
            if !c.is_segment(u, v, &surface.acute)
                || seen.insert((u.min(v), u.max(v)), ()).is_some()
            {
                continue;
            }
            let (Some(du), Some(dv)) = (vertex_of[u as usize], vertex_of[v as usize]) else {
                continue;
            };
            let (pu, pv) = (&c.points[u as usize], &c.points[v as usize]);
            // Along an input edge, or across the triangle.
            let edge = match (c.vertices[u as usize], c.vertices[v as usize]) {
                (ChamferVertex::OnEdge { edge, .. }, _)
                | (_, ChamferVertex::OnEdge { edge, .. })
                    if [u, v].iter().all(|&x| match c.vertices[x as usize] {
                        ChamferVertex::OnEdge { edge: e, .. } => e == edge,
                        ChamferVertex::Input(i) => edge.0 == i || edge.1 == i,
                        ChamferVertex::InFace { .. } => false,
                    }) =>
                {
                    Some(edge)
                }
                (ChamferVertex::Input(i), ChamferVertex::Input(j)) => Some((i.min(j), i.max(j))),
                _ => None,
            };
            let seg = match edge {
                Some((e0, e1)) => {
                    let (a, b) = (points[e0 as usize], points[e1 as usize]);
                    let (Some(s0), Some(s1)) = (line_param(pu, a, b), line_param(pv, a, b)) else {
                        return Err(RefineError::Construction);
                    };
                    Seg {
                        ends: [du, dv],
                        carrier: Carrier::Line { a, b },
                        params: [Param::S(s0), Param::S(s1)],
                    }
                }
                None => {
                    let (Some(w0), Some(w1)) = (tri_param(pu, corners), tri_param(pv, corners))
                    else {
                        return Err(RefineError::Construction);
                    };
                    Seg {
                        ends: [du, dv],
                        carrier: Carrier::Triangle { tri: corners },
                        params: [Param::W(w0), Param::W(w1)],
                    }
                }
            };
            segs.push(seg);
        }
    }
    for (i, j) in [
        (0, 1),
        (2, 3),
        (4, 5),
        (6, 7),
        (0, 2),
        (1, 3),
        (4, 6),
        (5, 7),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ] {
        segs.push(Seg {
            ends: [box_corners[i], box_corners[j]],
            carrier: Carrier::Line {
                a: corner(i),
                b: corner(j),
            },
            params: [Param::S(0.0), Param::S(1.0)],
        });
    }
    // Split what is missing or encroached, rechecking the subsegments at a
    // new vertex's neighbours.
    let mut alive = vec![true; segs.len()];
    let mut at_vertex: Vec<Vec<usize>> = vec![Vec::new(); dt.points().len()];
    for (i, s) in segs.iter().enumerate() {
        for &e in &s.ends {
            at_vertex[e as usize].push(i);
        }
    }
    let mut queue: Vec<usize> = (0..segs.len()).rev().collect();
    let mut steiner = 0;
    loop {
        while let Some(i) = queue.pop() {
            if !alive[i] {
                continue;
            }
            let s = segs[i];
            let [a, b] = s.ends;
            let check = crate::prof::scope("refine segments: check");
            let ok = segment_ok(&dt, [a, b]);
            drop(check);
            if ok {
                continue;
            }
            if dt.points().len() >= max_points {
                return Err(RefineError::Budget {
                    points: dt.points().len(),
                    limit: max_points,
                });
            }
            let (p, param) = s.middle().ok_or(RefineError::Construction)?;
            let insert = crate::prof::scope("refine segments: insert");
            let m = match dt.insert(p)? {
                Inserted::New(m) => m,
                Inserted::Duplicate(v) => {
                    return Err(RefineError::Duplicate {
                        point: dt.points()[v as usize].approx(),
                    });
                }
            };
            drop(insert);
            steiner += 1;
            alive[i] = false;
            at_vertex.resize(dt.points().len(), Vec::new());
            for (ends, params) in [
                ([a, m], [s.params[0], param]),
                ([m, b], [param, s.params[1]]),
            ] {
                let id = segs.len();
                segs.push(Seg {
                    ends,
                    carrier: s.carrier,
                    params,
                });
                alive.push(true);
                for &e in &ends {
                    at_vertex[e as usize].push(id);
                }
                queue.push(id);
            }
            // A subsegment the new vertex is not inside the diametral sphere of
            // stays a Delaunay edge; check the others again. (A final pass
            // catches what this misses: an edge lost in a cospherical tie, or a
            // sphere reached from beyond the neighbours.)
            let _s = crate::prof::scope("refine segments: requeue");
            let pm = dt.points()[m as usize];
            for w in dt.neighbours(m) {
                if let Some(list) = at_vertex.get(w as usize) {
                    for &j in list {
                        let [x, y] = segs[j].ends;
                        if alive[j]
                            && x != m
                            && y != m
                            && implicit::dot_sign(
                                &dt.points()[x as usize],
                                &pm,
                                &dt.points()[y as usize],
                            ) != Sign::Positive
                        {
                            queue.push(j);
                        }
                    }
                }
            }
        }
        // The final pass: every subsegment checked in full.
        let _s = crate::prof::scope("refine segments: final pass");
        queue = (0..segs.len())
            .filter(|&i| alive[i] && !segment_ok(&dt, segs[i].ends))
            .collect();
        if queue.is_empty() {
            break;
        }
    }
    let subsegments = segs
        .iter()
        .zip(&alive)
        .filter(|(_, a)| **a)
        .map(|(s, _)| s.ends)
        .collect();
    Ok(RefinedSegments {
        dt,
        subsegments,
        vertex_of,
        box_corners,
        steiner,
    })
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;
    use crate::chamfer::chamfer_filtered;

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

    /// Every final subsegment is an unencroached Delaunay edge, and the
    /// tetrahedralization is valid and Delaunay.
    fn check(name: &str, p: &[Point3], f: &[[u32; 3]]) -> RefinedSegments {
        let surface = chamfer_filtered(p, f, 1e-6).unwrap();
        let r = recover_segments_conforming(&surface, p, f, 0.2, 100_000).unwrap();
        r.dt.check().unwrap();
        r.dt.check_delaunay().unwrap();
        for s in &r.subsegments {
            assert!(segment_ok(&r.dt, *s), "{name}: {s:?}");
        }
        r
    }

    #[test]
    fn segments_of_a_cube_and_a_chamfered_tetrahedron_are_recovered() {
        let (p, f) = cube();
        let r = check("cube", &p, &f);
        // The cube's 12 edges (its diagonals are flat) and the box's 12.
        assert!(r.subsegments.len() >= 24);
        let p = vec![
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        let f = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        // Every corner and edge chamfered: four hexagons of bridge points.
        let r = check("tetrahedron", &p, &f);
        assert!(r.subsegments.len() >= 24 + 12);
    }
}
