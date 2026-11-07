//! Exact intersection tests for closed segments and triangles, built only on
//! the exact predicates and exact coordinate comparisons.
//!
//! Every test handles degenerate input (zero-length segments, collinear or
//! coincident triangle corners) without panicking and without an epsilon.
//!
//! [`triangles_intersect`] rests on one fact: two closed convex sets that meet
//! have an extreme point of their intersection on the boundary of one of them,
//! so two triangles meet exactly when an edge of one meets the other.

use crate::Point3;
use crate::predicates::{
    Sign, collinear, orient2d, orient2d_projected, orient3d, plane_axis, project,
};

fn minmax(a: f64, b: f64) -> (f64, f64) {
    if a <= b { (a, b) } else { (b, a) }
}

/// Closed box containment, used on points already known to be collinear.
fn within_box(p: Point3, a: Point3, b: Point3) -> bool {
    (0..3).all(|k| {
        let (lo, hi) = minmax(a[k], b[k]);
        lo <= p[k] && p[k] <= hi
    })
}

/// `p` lies on the closed segment `ab` (a point if `a == b`).
#[must_use]
pub fn point_on_segment(p: Point3, a: Point3, b: Point3) -> bool {
    collinear(a, b, p) && within_box(p, a, b)
}

fn on_segment_2d(r: [f64; 2], a: [f64; 2], b: [f64; 2]) -> bool {
    (0..2).all(|k| {
        let (lo, hi) = minmax(a[k], b[k]);
        lo <= r[k] && r[k] <= hi
    })
}

/// Closed 2D segment intersection (CLRS), exact.
fn segments_intersect_2d(p: [f64; 2], q: [f64; 2], x: [f64; 2], y: [f64; 2]) -> bool {
    let o1 = orient2d(p, q, x);
    let o2 = orient2d(p, q, y);
    let o3 = orient2d(x, y, p);
    let o4 = orient2d(x, y, q);
    (o1.opposes(o2) && o3.opposes(o4))
        || (o1 == Sign::Zero && on_segment_2d(x, p, q))
        || (o2 == Sign::Zero && on_segment_2d(y, p, q))
        || (o3 == Sign::Zero && on_segment_2d(p, x, y))
        || (o4 == Sign::Zero && on_segment_2d(q, x, y))
}

/// Closed 2D point-in-triangle for a triangle with nonzero area.
fn point_in_triangle_2d(p: [f64; 2], [a, b, c]: [[f64; 2]; 3]) -> bool {
    let s = orient2d(a, b, c);
    [orient2d(a, b, p), orient2d(b, c, p), orient2d(c, a, p)]
        .into_iter()
        .all(|o| o == s || o == Sign::Zero)
}

/// Closed 3D segment intersection, exact.
#[must_use]
pub fn segments_intersect(p: Point3, q: Point3, x: Point3, y: Point3) -> bool {
    if p == q {
        return point_on_segment(p, x, y);
    }
    if x == y {
        return point_on_segment(x, p, q);
    }
    if orient3d(p, q, x, y) != Sign::Zero {
        return false;
    }
    match plane_axis(p, q, x).or_else(|| plane_axis(p, q, y)) {
        Some(k) => {
            segments_intersect_2d(project(p, k), project(q, k), project(x, k), project(y, k))
        }
        None => {
            // All four on the line pq: compare along an axis where p and q
            // differ; coordinates along a line are monotone in its parameter.
            let k = (0..3).find(|&k| p[k] != q[k]).unwrap_or(0);
            let (lo1, hi1) = minmax(p[k], q[k]);
            let (lo2, hi2) = minmax(x[k], y[k]);
            lo1 <= hi2 && lo2 <= hi1
        }
    }
}

/// The convex hull of three collinear points as a segment (equal ends for a point).
fn collinear_hull(pts: [Point3; 3]) -> (Point3, Point3) {
    let Some(k) = (0..3).find(|&k| pts.iter().any(|p| p[k] != pts[0][k])) else {
        return (pts[0], pts[0]);
    };
    let lo = pts.iter().min_by(|a, b| a[k].total_cmp(&b[k])).copied();
    let hi = pts.iter().max_by(|a, b| a[k].total_cmp(&b[k])).copied();
    (lo.unwrap_or(pts[0]), hi.unwrap_or(pts[0]))
}

/// Closed segment `pq` against closed triangle `tri` (which may be degenerate), exact.
#[must_use]
pub fn segment_intersects_triangle(p: Point3, q: Point3, tri: [Point3; 3]) -> bool {
    let [a, b, c] = tri;
    let Some(k) = plane_axis(a, b, c) else {
        let (x, y) = collinear_hull(tri);
        return segments_intersect(p, q, x, y);
    };
    let flat = |r: Point3| project(r, k);
    let tri2 = [flat(a), flat(b), flat(c)];
    let sp = orient3d(a, b, c, p);
    let sq = orient3d(a, b, c, q);
    match (sp, sq) {
        (Sign::Zero, Sign::Zero) => {
            let (p2, q2) = (flat(p), flat(q));
            point_in_triangle_2d(p2, tri2)
                || point_in_triangle_2d(q2, tri2)
                || (0..3).any(|i| segments_intersect_2d(p2, q2, tri2[i], tri2[(i + 1) % 3]))
        }
        (Sign::Zero, _) => point_in_triangle_2d(flat(p), tri2),
        (_, Sign::Zero) => point_in_triangle_2d(flat(q), tri2),
        (s, t) if s == t => false,
        _ => {
            // pq crosses the plane; the crossing is inside the triangle iff the
            // line pq passes on the same side of all three edges.
            let o = [
                orient3d(p, q, a, b),
                orient3d(p, q, b, c),
                orient3d(p, q, c, a),
            ];
            !(o.contains(&Sign::Positive) && o.contains(&Sign::Negative))
        }
    }
}

fn boxes_disjoint(a: &[Point3; 3], b: &[Point3; 3]) -> bool {
    (0..3).any(|k| {
        let amin = a.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min);
        let amax = a.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max);
        let bmin = b.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min);
        let bmax = b.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max);
        amax < bmin || bmax < amin
    })
}

/// Do two closed triangles share at least one point? Exact; degenerate
/// triangles are treated as the segment or point they cover.
#[must_use]
pub fn triangles_intersect(a: [Point3; 3], b: [Point3; 3]) -> bool {
    if boxes_disjoint(&a, &b) {
        return false;
    }
    (0..3).any(|i| segment_intersects_triangle(a[i], a[(i + 1) % 3], b))
        || (0..3).any(|i| segment_intersects_triangle(b[i], b[(i + 1) % 3], a))
}

/// Do two faces of one mesh meet anywhere other than in the vertices and edge
/// they share by index? This is the self-intersection test: neighbours that
/// merely touch along their shared vertex or edge do not count.
///
/// The shared-vertex and shared-edge cases need both faces to have nonzero
/// area and distinct indices; otherwise this falls back to the plain closed
/// test, which reports any contact (degenerate faces are diagnosed separately).
#[must_use]
pub fn faces_intersect(ia: [u32; 3], a: [Point3; 3], ib: [u32; 3], b: [Point3; 3]) -> bool {
    let distinct = |t: [u32; 3]| t[0] != t[1] && t[1] != t[2] && t[0] != t[2];
    let proper = |p: [Point3; 3]| plane_axis(p[0], p[1], p[2]).is_some();
    if !(distinct(ia) && distinct(ib) && proper(a) && proper(b)) {
        return triangles_intersect(a, b);
    }
    let shared: Vec<usize> = (0..3).filter(|&i| ib.contains(&ia[i])).collect();
    match shared.len() {
        0 => triangles_intersect(a, b),
        1 => {
            let v = ia[shared[0]];
            let ra = rotate_to(a, ia, v);
            let rb = rotate_to(b, ib, v);
            beyond_shared_vertex(ra, rb)
        }
        2 => {
            let third_a = (0..3).find(|i| !shared.contains(i)).unwrap_or(0);
            let third_b = (0..3).find(|&i| !ia.contains(&ib[i])).unwrap_or(0);
            let (u, w) = (a[shared[0]], a[shared[1]]);
            beyond_shared_edge(u, w, a[third_a], b[third_b])
        }
        _ => true,
    }
}

/// Rotate corners so the one with index `v` comes first (keeps cyclic order).
fn rotate_to(p: [Point3; 3], ids: [u32; 3], v: u32) -> [Point3; 3] {
    let i = ids.iter().position(|&x| x == v).unwrap_or(0);
    [p[i], p[(i + 1) % 3], p[(i + 2) % 3]]
}

/// `a = (v, a1, a2)` and `b = (v, b1, b2)` share only `v`. The intersection
/// is convex and contains `v`; its point farthest from `v`, if not `v`, lies on
/// an edge of one face, so checking each edge suffices.
fn beyond_shared_vertex(a: [Point3; 3], b: [Point3; 3]) -> bool {
    segment_intersects_triangle(a[1], a[2], b)
        || segment_intersects_triangle(b[1], b[2], a)
        || edge_enters(a[1], b)
        || edge_enters(a[2], b)
        || edge_enters(b[1], a)
        || edge_enters(b[2], a)
}

/// Does the segment from `t[0]` to `x` meet triangle `t = (v, t1, t2)` in more
/// than `v`? Only if `x` is coplanar and the direction lies in the closed
/// angle of `t` at `v`.
fn edge_enters(x: Point3, t: [Point3; 3]) -> bool {
    let [v, t1, t2] = t;
    if orient3d(v, t1, t2, x) != Sign::Zero {
        return false;
    }
    let Some(k) = plane_axis(v, t1, t2) else {
        return false;
    };
    let s = orient2d_projected(v, t1, t2, k);
    let o1 = orient2d_projected(v, t1, x, k);
    let o2 = orient2d_projected(v, x, t2, k);
    (o1 == s || o1 == Sign::Zero) && (o2 == s || o2 == Sign::Zero)
}

/// Faces `(u, w, a)` and `(u, w, b)` share edge `uw`: they overlap beyond it
/// only when coplanar with `a` and `b` on the same side of the edge.
fn beyond_shared_edge(u: Point3, w: Point3, a: Point3, b: Point3) -> bool {
    if orient3d(u, w, a, b) != Sign::Zero {
        return false;
    }
    let Some(k) = plane_axis(u, w, a) else {
        return true;
    };
    orient2d_projected(u, w, a, k) == orient2d_projected(u, w, b, k)
}

/// How two faces from different solids meet.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Contact {
    /// No common point beyond shared vertices or edges.
    None,
    /// They touch along a vertex, an edge or a line without either passing
    /// through the other, or coplanar faces meet with no common area.
    Touching,
    /// Each face passes through the other's plane where they meet.
    Crossing,
    /// Coplanar with a common area and normals pointing the same way: the
    /// solids behind them are on the same side, so their volumes overlap.
    CoplanarAligned,
    /// Coplanar with a common area and opposite normals: a shared interface.
    CoplanarOpposed,
}

/// Classify the contact between faces of two different closed surfaces, for
/// telling overlapping solids from solids that share an interface. Faces of
/// zero area only ever touch.
#[must_use]
pub fn classify_contact(ia: [u32; 3], a: [Point3; 3], ib: [u32; 3], b: [Point3; 3]) -> Contact {
    if !faces_intersect(ia, a, ib, b) {
        return Contact::None;
    }
    let (Some(ka), Some(_)) = (plane_axis(a[0], a[1], a[2]), plane_axis(b[0], b[1], b[2])) else {
        return Contact::Touching;
    };
    let side_b = b.map(|q| orient3d(a[0], a[1], a[2], q));
    if side_b.iter().all(|&s| s == Sign::Zero) {
        if !coplanar_interiors_overlap(a, b) {
            return Contact::Touching;
        }
        let na = orient2d_projected(a[0], a[1], a[2], ka);
        let nb = orient2d_projected(b[0], b[1], b[2], ka);
        return if na == nb {
            Contact::CoplanarAligned
        } else {
            Contact::CoplanarOpposed
        };
    }
    let side_a = a.map(|q| orient3d(b[0], b[1], b[2], q));
    let straddles = |s: [Sign; 3]| s.contains(&Sign::Positive) && s.contains(&Sign::Negative);
    if straddles(side_a) && straddles(side_b) {
        Contact::Crossing
    } else {
        Contact::Touching
    }
}

/// Do two coplanar triangles of nonzero area share interior points (a common
/// area, not just an edge or a corner)? Separating-axis test with exact
/// orientations: interiors are disjoint exactly when some edge line has the
/// other triangle entirely on its closed outer side.
#[must_use]
pub fn coplanar_interiors_overlap(a: [Point3; 3], b: [Point3; 3]) -> bool {
    let Some(k) = plane_axis(a[0], a[1], a[2]) else {
        return false;
    };
    let a2 = a.map(|p| project(p, k));
    let b2 = b.map(|p| project(p, k));
    let separated_by_edges_of = |t: [[f64; 2]; 3], other: [[f64; 2]; 3]| {
        let s = orient2d(t[0], t[1], t[2]);
        s != Sign::Zero
            && (0..3).any(|i| {
                let (p, q) = (t[i], t[(i + 1) % 3]);
                other.iter().all(|&r| orient2d(p, q, r) != s)
            })
    };
    !(separated_by_edges_of(a2, b2) || separated_by_edges_of(b2, a2))
}

#[cfg(test)]
mod tests {
    use super::*;

    const O: Point3 = [0.0, 0.0, 0.0];
    const X: Point3 = [1.0, 0.0, 0.0];
    const Y: Point3 = [0.0, 1.0, 0.0];
    const Z: Point3 = [0.0, 0.0, 1.0];

    #[test]
    fn piercing_and_touching() {
        let base = [O, [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
        let pierce = [[1.0, 1.0, -1.0], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0]];
        let touch = [[1.0, 1.0, 0.0], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0]];
        let above = [[1.0, 1.0, 0.5], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0]];
        assert!(triangles_intersect(base, pierce));
        assert!(triangles_intersect(base, touch));
        assert!(!triangles_intersect(base, above));
    }

    #[test]
    fn coplanar_cases() {
        let big = [O, [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
        let inside = [[1.0, 1.0, 0.0], [2.0, 1.0, 0.0], [1.0, 2.0, 0.0]];
        let apart = [[5.0, 5.0, 0.0], [6.0, 5.0, 0.0], [5.0, 6.0, 0.0]];
        assert!(triangles_intersect(big, inside));
        assert!(!triangles_intersect(big, apart));
    }

    #[test]
    fn neighbours_sharing_an_edge() {
        let ids_a = [0, 1, 2];
        let ids_b = [1, 0, 3];
        let a = [O, X, Y];
        // Folded neighbour: only the edge in common.
        assert!(!faces_intersect(ids_a, a, ids_b, [X, O, Z]));
        // Coplanar on the other side of the edge: fine.
        assert!(!faces_intersect(ids_a, a, ids_b, [X, O, [0.0, -1.0, 0.0]]));
        // Coplanar on the same side: overlap.
        assert!(faces_intersect(ids_a, a, ids_b, [X, O, [0.5, 0.5, 0.0]]));
    }

    #[test]
    fn neighbours_sharing_a_vertex() {
        let a = [O, [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let ids_a = [0, 1, 2];
        let ids_b = [0, 3, 4];
        // Fan neighbour meeting only at O.
        assert!(!faces_intersect(
            ids_a,
            a,
            ids_b,
            [O, [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]
        ));
        // Passes through a: intersects beyond O.
        assert!(faces_intersect(
            ids_a,
            a,
            ids_b,
            [O, [1.0, 1.0, 1.0], [1.0, 0.5, -1.0]]
        ));
        // Edge from O running inside a, coplanar.
        assert!(faces_intersect(
            ids_a,
            a,
            ids_b,
            [O, [0.5, 0.5, 0.0], [0.0, 0.0, 1.0]]
        ));
    }

    #[test]
    fn degenerate_inputs() {
        assert!(segments_intersect(O, X, [0.5, 0.0, 0.0], [0.5, 0.0, 0.0]));
        assert!(segments_intersect(O, X, [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]));
        assert!(!segments_intersect(O, X, [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]));
        let point = [Y, Y, Y];
        assert!(segment_intersects_triangle(O, [0.0, 2.0, 0.0], point));
        let line = [O, [2.0, 0.0, 0.0], X];
        assert!(triangles_intersect(
            line,
            [[1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 1.0]]
        ));
    }

    #[test]
    fn contact_kinds() {
        let ids_a = [0, 1, 2];
        let ids_b = [3, 4, 5];
        let a = [O, [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let opposed = [[0.5, 0.5, 0.0], [0.5, 1.0, 0.0], [1.0, 0.5, 0.0]];
        let aligned = [[0.5, 0.5, 0.0], [1.0, 0.5, 0.0], [0.5, 1.0, 0.0]];
        let crossing = [[0.5, 0.5, -1.0], [0.5, 0.5, 1.0], [1.0, 0.2, 0.0]];
        let resting = [[0.5, 0.5, 0.0], [0.5, 0.5, 1.0], [1.0, 0.2, 1.0]];
        let edge_to_edge = [[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [2.0, 0.0, 0.0]];
        assert_eq!(
            classify_contact(ids_a, a, ids_b, opposed),
            Contact::CoplanarOpposed
        );
        assert_eq!(
            classify_contact(ids_a, a, ids_b, aligned),
            Contact::CoplanarAligned
        );
        assert_eq!(
            classify_contact(ids_a, a, ids_b, crossing),
            Contact::Crossing
        );
        assert_eq!(
            classify_contact(ids_a, a, ids_b, resting),
            Contact::Touching
        );
        assert_eq!(
            classify_contact(ids_a, a, ids_b, edge_to_edge),
            Contact::Touching
        );
        let far = [[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0, 6.0, 5.0]];
        assert_eq!(classify_contact(ids_a, a, ids_b, far), Contact::None);
    }

    #[test]
    fn coplanar_overlap_needs_area() {
        let a = [O, X, Y];
        assert!(coplanar_interiors_overlap(
            a,
            [[0.2, 0.2, 0.0], [0.4, 0.2, 0.0], [0.2, 0.4, 0.0]]
        ));
        assert!(coplanar_interiors_overlap(a, a));
        assert!(!coplanar_interiors_overlap(a, [X, Y, [1.0, 1.0, 0.0]]));
        assert!(!coplanar_interiors_overlap(
            a,
            [O, [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]
        ));
    }
}
