//! Exact intersections of a triangle soup (Milestone 16, stage 1): where its
//! triangles meet, as implicit points and the segments between them, for
//! subdividing each triangle into an arrangement (after Cherchi et al.,
//! 2022).
//!
//! For every pair of triangles that meet beyond the vertices and edges they
//! share by index ([`crate::faces_intersect`]):
//!
//! - *In general position* the two planes meet in a line, and the triangles
//!   in a segment or a point of it: its ends are among the points where each
//!   triangle's boundary meets the other's plane (a vertex on the plane, or
//!   an edge crossing it: the LPI of the edge's line and the plane), those
//!   that lie in the other triangle (closed). By convexity at most two of
//!   them are distinct.
//! - *Coplanar* triangles meet in a convex polygon; each gets the parts of
//!   the other's edges inside it. Two edges crossing in the common plane meet
//!   at the LPI of one edge's line and the plane through the other edge
//!   parallel to the projection axis: that plane's third point differs from
//!   an edge end only in the axis coordinate, so it is exact.
//!
//! Every decision is an exact predicate on explicit or implicit points;
//! nothing is rounded. Points are indexed globally (the input points first,
//! explicit, then the new ones) and deduplicated exactly
//! ([`implicit::coincident`]), so a point found from several pairs is one
//! point. A point on an input edge is given to both triangles of the edge,
//! so their subdivisions agree along it.

mod cdt;
mod label;
mod round;
mod subdivide;

use std::collections::BTreeMap;

use thiserror::Error;

pub use label::{Resolved, resolve, resolve_with};
pub use round::{Rounded, resolve_rounded};
pub use subdivide::{Arrangement, arrangement};

use crate::implicit::{self, ImplicitError, Point};
use crate::predicates::{Sign, orient3d, plane_axis};
use crate::{Bvh, GeomError, Point3, TriangleSoup, faces_intersect};

/// Why the intersections could not be computed.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum ArrangementError {
    #[error("{0}")]
    Surface(#[from] GeomError),
    #[error("triangle {triangle} has zero area")]
    Degenerate { triangle: u32 },
    #[error("an intersection point could not be built: {0}")]
    Implicit(#[from] ImplicitError),
    #[error("a triangle could not be subdivided: {0}")]
    Triangulation(String),
    #[error("the surface is not closed: {edges} edges are not used as often in each direction")]
    Open { edges: usize },
    #[error("no segment could be cast clear of every edge and vertex")]
    Cast,
    #[error("rounded to f64, the surface still has {pairs} intersecting triangle pairs")]
    Rounding { pairs: usize },
}

/// What a constraint segment lies on besides its triangle's plane, so its
/// crossings with other segments can be built exactly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Support {
    /// The plane through three explicit points (the other triangle's).
    Plane([Point3; 3]),
    /// The line through two explicit points (an edge of a coplanar
    /// triangle).
    Line([Point3; 2]),
}

/// A constraint segment: its ends by global point index (lower first) and
/// what it lies on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    pub ends: [u32; 2],
    pub support: Support,
}

/// What a triangle must be subdivided by: points on it other than its
/// corners (inside it or on its edges), and segments across it, by global
/// point index. Sorted and without repeats (one segment per pair of ends).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Constraints {
    pub points: Vec<u32>,
    pub segments: Vec<Segment>,
}

impl Constraints {
    /// Nothing to subdivide by.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.points.is_empty() && self.segments.is_empty()
    }
}

/// The intersections of a soup.
#[derive(Debug, Clone)]
pub struct Intersections {
    /// The input points (explicit, same indices) and then the points the
    /// intersections made; no two coincide, except input points that do.
    pub points: Vec<Point>,
    /// Each input triangle's constraints.
    pub constraints: Vec<Constraints>,
    /// The pairs of triangles that meet beyond their shared vertices and
    /// edges.
    pub pairs: usize,
    /// Of those, the coplanar pairs, sorted.
    pub coplanar: Vec<(u32, u32)>,
}

/// A point of a triangle's boundary, and where on the triangle it lies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum On {
    Corner,
    /// On the edge from corner `i` to corner `i + 1`.
    Edge(usize),
    Inside,
}

fn edge_key(u: u32, v: u32) -> (u32, u32) {
    (u.min(v), u.max(v))
}

/// The intersections of `triangles` over `points`.
///
/// # Errors
/// [`ArrangementError::Surface`] for an invalid soup,
/// [`ArrangementError::Degenerate`] for a triangle of zero area, and
/// [`ArrangementError::Implicit`] if a point that must exist could not be
/// built (a bug).
pub fn intersections(
    points: &[Point3],
    triangles: &[[u32; 3]],
) -> Result<Intersections, ArrangementError> {
    let soup = TriangleSoup::new(points, triangles)?;
    let pos = |t: [u32; 3]| t.map(|v| points[v as usize]);
    let mut axes = Vec::with_capacity(triangles.len());
    for (i, &t) in triangles.iter().enumerate() {
        let [a, b, c] = pos(t);
        let Some(k) = plane_axis(a, b, c) else {
            return Err(ArrangementError::Degenerate { triangle: i as u32 });
        };
        axes.push(k);
    }
    let mut b = Builder {
        inputs: points.len(),
        points: points.iter().map(|&p| Point::from(p)).collect(),
        keys: BTreeMap::new(),
        constraints: vec![Constraints::default(); triangles.len()],
        edge_points: BTreeMap::new(),
    };
    let bvh = Bvh::new(soup);
    let mut pairs = 0;
    let mut coplanar = Vec::new();
    for (f, g) in bvh.self_overlaps() {
        let (tf, tg) = (triangles[f as usize], triangles[g as usize]);
        if !faces_intersect(tf, pos(tf), tg, pos(tg)) {
            continue;
        }
        pairs += 1;
        let ta = Tri {
            index: f,
            v: tf,
            p: pos(tf),
            axis: axes[f as usize],
        };
        let tb = Tri {
            index: g,
            v: tg,
            p: pos(tg),
            axis: axes[g as usize],
        };
        let sides = tb.p.map(|q| orient3d(ta.p[0], ta.p[1], ta.p[2], q));
        if sides.iter().all(|&s| s == Sign::Zero) {
            coplanar.push((f, g));
            b.coplanar(&ta, &tb)?;
            b.coplanar(&tb, &ta)?;
        } else {
            b.crossing(&ta, &tb)?;
        }
    }
    let mut out = b.finish(triangles);
    out.pairs = pairs;
    out.coplanar = coplanar;
    Ok(out)
}

/// A triangle with what the tests need.
struct Tri {
    index: u32,
    v: [u32; 3],
    p: [Point3; 3],
    /// An axis its projection along keeps its area.
    axis: usize,
}

impl Tri {
    /// Where `q`, in this triangle's plane, lies on it; `None` outside.
    fn locate(&self, q: &Point) -> Option<On> {
        let k = self.axis;
        let c = self.p.map(Point::from);
        let s = implicit::orient2d_projected(&c[0], &c[1], &c[2], k);
        let o = [0, 1, 2].map(|i| implicit::orient2d_projected(&c[i], &c[(i + 1) % 3], q, k));
        if o.iter().any(|&x| x != s && x != Sign::Zero) {
            return None;
        }
        let zeros: Vec<usize> = (0..3).filter(|&i| o[i] == Sign::Zero).collect();
        Some(match zeros.as_slice() {
            [] => On::Inside,
            [i] => On::Edge(*i),
            _ => On::Corner,
        })
    }
}

/// Key of a constructed point, so the same construction is made once.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Key {
    /// An edge's line through a triangle's plane.
    EdgePlane((u32, u32), u32),
    /// Two coplanar edges crossing.
    EdgeEdge((u32, u32), (u32, u32)),
}

struct Builder {
    /// How many of `points` are the input's.
    inputs: usize,
    points: Vec<Point>,
    keys: BTreeMap<Key, u32>,
    constraints: Vec<Constraints>,
    /// Points on each input edge, other than its ends.
    edge_points: BTreeMap<(u32, u32), Vec<u32>>,
}

impl Builder {
    fn make(
        &mut self,
        key: Key,
        build: impl FnOnce() -> Result<Point, ImplicitError>,
    ) -> Result<u32, ArrangementError> {
        if let Some(&i) = self.keys.get(&key) {
            return Ok(i);
        }
        let i = self.points.len() as u32;
        self.points.push(build()?);
        self.keys.insert(key, i);
        Ok(i)
    }

    /// Record where point `i` lies on triangle `t`.
    fn place(&mut self, t: &Tri, i: u32, on: On) {
        match on {
            On::Corner => {}
            On::Edge(e) => self
                .edge_points
                .entry(edge_key(t.v[e], t.v[(e + 1) % 3]))
                .or_default()
                .push(i),
            On::Inside => self.constraints[t.index as usize].points.push(i),
        }
    }

    /// `a` and `b` in planes that meet in a line.
    fn crossing(&mut self, a: &Tri, b: &Tri) -> Result<(), ArrangementError> {
        // Each triangle's boundary against the other's plane.
        let mut candidates: Vec<(u32, Option<On>, Option<On>)> = Vec::new();
        for (this, other, first) in [(a, b, true), (b, a, false)] {
            let side = this
                .p
                .map(|q| orient3d(other.p[0], other.p[1], other.p[2], q));
            for i in 0..3 {
                let j = (i + 1) % 3;
                if side[i] == Sign::Zero {
                    let q = Point::from(this.p[i]);
                    let on_other = other.locate(&q);
                    candidates.push(if first {
                        (this.v[i], Some(On::Corner), on_other)
                    } else {
                        (this.v[i], on_other, Some(On::Corner))
                    });
                }
                if side[i].opposes(side[j]) {
                    let line = [this.p[i], this.p[j]];
                    let plane = other.p;
                    let key = Key::EdgePlane(edge_key(this.v[i], this.v[j]), other.index);
                    let id = self.make(key, || Point::lpi(line, plane))?;
                    let on_other = other.locate(&self.points[id as usize]);
                    candidates.push(if first {
                        (id, Some(On::Edge(i)), on_other)
                    } else {
                        (id, on_other, Some(On::Edge(i)))
                    });
                }
            }
        }
        // The ends of the intersection: the candidates in both triangles.
        let mut ends: Vec<(u32, On, On)> = Vec::new();
        for (id, on_a, on_b) in candidates {
            let (Some(on_a), Some(on_b)) = (on_a, on_b) else {
                continue;
            };
            let q = self.points[id as usize];
            if ends
                .iter()
                .any(|&(e, _, _)| implicit::coincident(&self.points[e as usize], &q))
            {
                continue;
            }
            ends.push((id, on_a, on_b));
        }
        for &(id, on_a, on_b) in &ends {
            self.place(a, id, on_a);
            self.place(b, id, on_b);
        }
        if let [(p, ..), (q, ..)] = ends.as_slice() {
            for (t, other) in [(a, b), (b, a)] {
                self.constraints[t.index as usize].segments.push(Segment {
                    ends: [*p, *q],
                    support: Support::Plane(other.p),
                });
            }
        }
        Ok(())
    }

    /// `a` and `b` coplanar: give `a` the parts of `b`'s edges inside it.
    fn coplanar(&mut self, a: &Tri, b: &Tri) -> Result<(), ArrangementError> {
        let k = a.axis;
        let ca = a.p.map(Point::from);
        for e in 0..3 {
            let (u, w) = (b.v[e], b.v[(e + 1) % 3]);
            let (pu, pw) = (b.p[e], b.p[(e + 1) % 3]);
            let mut ends: Vec<(u32, On)> = Vec::new();
            let add = |builder: &Builder, id: u32, on: On, ends: &mut Vec<(u32, On)>| {
                let q = builder.points[id as usize];
                if !ends
                    .iter()
                    .any(|&(x, _)| implicit::coincident(&builder.points[x as usize], &q))
                {
                    ends.push((id, on));
                }
            };
            for (id, p) in [(u, pu), (w, pw)] {
                if let Some(on) = a.locate(&Point::from(p)) {
                    add(self, id, on, &mut ends);
                }
            }
            // Proper crossings of the edge with a's edges.
            let (qu, qw) = (Point::from(pu), Point::from(pw));
            for i in 0..3 {
                let j = (i + 1) % 3;
                let s1 = implicit::orient2d_projected(&ca[i], &ca[j], &qu, k);
                let s2 = implicit::orient2d_projected(&ca[i], &ca[j], &qw, k);
                let s3 = implicit::orient2d_projected(&qu, &qw, &ca[i], k);
                let s4 = implicit::orient2d_projected(&qu, &qw, &ca[j], k);
                if !(s1.opposes(s2) && s3.opposes(s4)) {
                    continue;
                }
                let ea = edge_key(a.v[i], a.v[j]);
                let eb = edge_key(u, w);
                // One construction per pair of edges: the lower edge's line
                // through the plane along the higher edge.
                let (line, along) = if ea < eb {
                    ([a.p[i], a.p[j]], [pu, pw])
                } else {
                    ([pu, pw], [a.p[i], a.p[j]])
                };
                let key = Key::EdgeEdge(ea.min(eb), ea.max(eb));
                let id = self.make(key, || Point::lpi(line, vertical_plane(along, k)))?;
                add(self, id, On::Edge(i), &mut ends);
            }
            for &(id, on) in &ends {
                self.place(a, id, on);
            }
            if let [(p, _), (q, _)] = ends.as_slice() {
                self.constraints[a.index as usize].segments.push(Segment {
                    ends: [*p, *q],
                    support: Support::Line([pu, pw]),
                });
            }
        }
        Ok(())
    }

    /// Merge coincident points, give each triangle the points on its edges,
    /// and tidy the lists.
    fn finish(mut self, triangles: &[[u32; 3]]) -> Intersections {
        let map = dedupe(&self.points);
        let r = |i: u32| map[i as usize];
        let mut edge_points: BTreeMap<(u32, u32), Vec<u32>> = BTreeMap::new();
        for (e, ps) in std::mem::take(&mut self.edge_points) {
            let entry = edge_points.entry(e).or_default();
            entry.extend(ps.into_iter().map(r).filter(|&i| i != e.0 && i != e.1));
        }
        let mut constraints = std::mem::take(&mut self.constraints);
        for (t, c) in triangles.iter().zip(constraints.iter_mut()) {
            let corners = t.map(r);
            let mut ps: Vec<u32> = c.points.iter().map(|&i| r(i)).collect();
            for e in 0..3 {
                if let Some(on) = edge_points.get(&edge_key(t[e], t[(e + 1) % 3])) {
                    ps.extend(on.iter().copied());
                }
            }
            ps.retain(|i| !corners.contains(i));
            ps.sort_unstable();
            ps.dedup();
            let mut ss: Vec<Segment> = c
                .segments
                .iter()
                .map(|s| {
                    let [p, q] = s.ends.map(r);
                    Segment {
                        ends: [p.min(q), p.max(q)],
                        support: s.support,
                    }
                })
                .filter(|s| s.ends[0] != s.ends[1])
                .collect();
            // One segment per pair of ends (any support will do).
            ss.sort_by_key(|s| s.ends);
            ss.dedup_by_key(|s| s.ends);
            *c = Constraints {
                points: ps,
                segments: ss,
            };
        }
        // Keep the input points (with their indices: a point merges into
        // the lowest index) and the new points a constraint uses,
        // renumbering them.
        let mut used = vec![false; self.points.len()];
        for c in &constraints {
            for &i in c
                .points
                .iter()
                .chain(c.segments.iter().flat_map(|s| &s.ends))
            {
                used[i as usize] = true;
            }
        }
        let mut keep = vec![u32::MAX; self.points.len()];
        let mut points = Vec::new();
        for (i, p) in self.points.iter().enumerate() {
            if i < self.inputs || (map[i] as usize == i && used[i]) {
                keep[i] = points.len() as u32;
                points.push(*p);
            }
        }
        let renumber = |i: u32| keep[i as usize];
        for c in &mut constraints {
            for i in &mut c.points {
                *i = renumber(*i);
            }
            for s in &mut c.segments {
                s.ends = s.ends.map(renumber);
            }
        }
        Intersections {
            points,
            constraints,
            pairs: 0,
            coplanar: Vec::new(),
        }
    }
}

/// The plane through edge `ab` parallel to axis `k`, by three explicit
/// points: `a`, `b`, and `a` moved along `k` only.
pub(crate) fn vertical_plane([a, b]: [Point3; 2], k: usize) -> [Point3; 3] {
    let mut c = a;
    c[k] = a[k] + a[k].abs().max(1.0);
    [a, b, c]
}

/// Each point's representative: the lowest index of the points it
/// coincides with. Candidates are pairs whose certain bounds overlap.
pub(crate) fn dedupe(points: &[Point]) -> Vec<u32> {
    let n = points.len();
    let mut rep: Vec<u32> = (0..n as u32).collect();
    let bounds: Vec<[[f64; 2]; 3]> = points.iter().map(Point::bounds).collect();
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by(|&i, &j| {
        bounds[i as usize][0][0]
            .total_cmp(&bounds[j as usize][0][0])
            .then(i.cmp(&j))
    });
    fn find(rep: &mut [u32], i: u32) -> u32 {
        let mut r = i;
        while rep[r as usize] != r {
            r = rep[r as usize];
        }
        let mut x = i;
        while rep[x as usize] != r {
            let next = rep[x as usize];
            rep[x as usize] = r;
            x = next;
        }
        r
    }
    let overlaps = |i: usize, j: usize| {
        (0..3).all(|k| bounds[i][k][0] <= bounds[j][k][1] && bounds[j][k][0] <= bounds[i][k][1])
    };
    // A sweep along x over the boxes' low ends.
    let mut active: Vec<u32> = Vec::new();
    for &i in &order {
        let lo = bounds[i as usize][0][0];
        active.retain(|&j| bounds[j as usize][0][1] >= lo);
        for &j in &active {
            if overlaps(i as usize, j as usize)
                && implicit::coincident(&points[i as usize], &points[j as usize])
            {
                let (ri, rj) = (find(&mut rep, i), find(&mut rep, j));
                let (lo, hi) = (ri.min(rj), ri.max(rj));
                rep[hi as usize] = lo;
            }
        }
        active.push(i);
    }
    (0..n as u32).map(|i| find(&mut rep, i)).collect()
}
