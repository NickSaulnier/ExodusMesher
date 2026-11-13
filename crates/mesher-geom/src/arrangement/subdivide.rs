//! Stage 2 of the arrangement (Milestone 16): each intersected triangle
//! subdivided along its constraints.
//!
//! Coplanar triangles that overlap are subdivided together, as one cluster,
//! so a region they share becomes the same pieces in each (by construction,
//! not by agreeing tie-breaks in separate triangulations). Per cluster, in
//! the projection along an axis that keeps its area:
//!
//! 1. the points: every member's corners and constraint points;
//! 2. the segments: every member's edges (on their explicit lines) and
//!    constraint segments;
//! 3. where two segments cross properly, a new point, built from what they
//!    lie on: two planes give the TPI with the cluster's plane, a line and a
//!    plane their LPI, two lines the LPI of one with the plane along the
//!    other parallel to the axis;
//! 4. each segment split at the points on it (exactly collinear, strictly
//!    between its ends);
//! 5. a constrained Delaunay triangulation of the points and pieces
//!    ([`super::cdt`]); each triangle goes to every member it lies in, with
//!    that member's orientation.
//!
//! Points made in different clusters (a triple point is found in each of
//! its three triangles) are merged exactly at the end.

use std::collections::BTreeMap;

use super::cdt::{Cdt, CdtError};
use super::{ArrangementError, Intersections, Support, dedupe, intersections, vertical_plane};
use crate::Point3;
use crate::implicit::{self, Point};
use crate::predicates::{Sign, plane_axis};

/// A triangle soup cut along its intersections.
#[derive(Debug, Clone)]
pub struct Arrangement {
    /// The input points (same indices), then the points the intersections
    /// made; no two coincide (unless input points do).
    pub points: Vec<Point>,
    /// Pieces, each oriented like the input triangle it lies on. A piece of
    /// a region shared by coplanar triangles appears once per triangle.
    pub triangles: Vec<[u32; 3]>,
    /// The input triangle each piece lies on.
    pub origin: Vec<u32>,
    /// Intersecting pairs of input triangles, and of them coplanar ones.
    pub pairs: usize,
    pub coplanar_pairs: usize,
}

/// The arrangement of `triangles` over `points`.
///
/// # Errors
/// As [`intersections`], and [`ArrangementError::Triangulation`] if a
/// cluster could not be triangulated (a bug).
pub fn arrangement(
    points: &[Point3],
    triangles: &[[u32; 3]],
) -> Result<Arrangement, ArrangementError> {
    let x = intersections(points, triangles)?;
    subdivide(points, triangles, x)
}

fn subdivide(
    input: &[Point3],
    triangles: &[[u32; 3]],
    x: Intersections,
) -> Result<Arrangement, ArrangementError> {
    let Intersections {
        mut points,
        constraints,
        pairs,
        coplanar,
    } = x;
    // Clusters of coplanar overlapping triangles.
    let n = triangles.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], i: usize) -> usize {
        let mut r = i;
        while parent[r] != r {
            r = parent[r];
        }
        let mut x = i;
        while parent[x] != r {
            let next = parent[x];
            parent[x] = r;
            x = next;
        }
        r
    }
    for &(f, g) in &coplanar {
        let (a, b) = (find(&mut parent, f as usize), find(&mut parent, g as usize));
        if a != b {
            parent[a.max(b)] = a.min(b);
        }
    }
    let mut clusters: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for t in 0..n {
        let r = find(&mut parent, t);
        clusters.entry(r).or_default().push(t);
    }

    let mut out_tris: Vec<[u32; 3]> = Vec::new();
    let mut origin: Vec<u32> = Vec::new();
    for members in clusters.values() {
        if let [t] = members.as_slice()
            && constraints[*t].is_empty()
        {
            out_tris.push(triangles[*t]);
            origin.push(*t as u32);
            continue;
        }
        cluster(
            input,
            triangles,
            &constraints,
            members,
            &mut points,
            &mut out_tris,
            &mut origin,
        )?;
    }

    // Merge points made more than once, and drop those no piece uses.
    let map = dedupe(&points);
    let mut used = vec![false; points.len()];
    for t in &mut out_tris {
        *t = t.map(|v| map[v as usize]);
        for &v in t.iter() {
            used[v as usize] = true;
        }
    }
    let mut keep = vec![u32::MAX; points.len()];
    let mut kept = Vec::new();
    for (i, p) in points.iter().enumerate() {
        if i < input.len() || (map[i] as usize == i && used[i]) {
            keep[i] = kept.len() as u32;
            kept.push(*p);
        }
    }
    let mut tris = Vec::with_capacity(out_tris.len());
    let mut orig = Vec::with_capacity(out_tris.len());
    for (t, o) in out_tris.into_iter().zip(origin) {
        let t = t.map(|v| keep[v as usize]);
        if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
            return Err(ArrangementError::Triangulation(
                "a piece collapsed when its points merged".into(),
            ));
        }
        tris.push(t);
        orig.push(o);
    }
    Ok(Arrangement {
        points: kept,
        triangles: tris,
        origin: orig,
        pairs,
        coplanar_pairs: coplanar.len(),
    })
}

/// A segment in a cluster: its ends (local) and what it lies on.
#[derive(Debug, Clone, Copy)]
struct Seg {
    ends: [usize; 2],
    support: Support,
}

#[allow(clippy::too_many_arguments)]
fn cluster(
    input: &[Point3],
    triangles: &[[u32; 3]],
    constraints: &[super::Constraints],
    members: &[usize],
    points: &mut Vec<Point>,
    out_tris: &mut Vec<[u32; 3]>,
    origin: &mut Vec<u32>,
) -> Result<(), ArrangementError> {
    let corners = |t: usize| triangles[t].map(|v| input[v as usize]);
    let [a, b, c] = corners(members[0]);
    let axis = plane_axis(a, b, c).ok_or(ArrangementError::Degenerate {
        triangle: members[0] as u32,
    })?;
    let plane = [a, b, c];

    // 1-2. Points (global ids, local by position) and segments.
    let mut ids: Vec<u32> = Vec::new();
    for &m in members {
        ids.extend(triangles[m]);
        ids.extend(constraints[m].points.iter().copied());
    }
    ids.sort_unstable();
    ids.dedup();
    let local = |ids: &[u32], g: u32| ids.binary_search(&g).unwrap_or(usize::MAX);
    let mut segs: Vec<Seg> = Vec::new();
    for &m in members {
        let t = triangles[m];
        for k in 0..3 {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            segs.push(Seg {
                ends: [local(&ids, u), local(&ids, v)],
                support: Support::Line([input[u as usize], input[v as usize]]),
            });
        }
        for s in &constraints[m].segments {
            segs.push(Seg {
                ends: s.ends.map(|g| local(&ids, g)),
                support: s.support,
            });
        }
    }
    if segs.iter().flat_map(|s| s.ends).any(|i| i == usize::MAX) {
        return Err(ArrangementError::Triangulation(
            "a segment end is not a point of its cluster".into(),
        ));
    }
    let mut pts: Vec<Point> = ids.iter().map(|&g| points[g as usize]).collect();
    let mut gid: Vec<Option<u32>> = ids.iter().map(|&g| Some(g)).collect();

    // 3. Proper crossings.
    let orient = |pts: &[Point], i: usize, j: usize, k: usize| {
        implicit::orient2d_projected(&pts[i], &pts[j], &pts[k], axis)
    };
    let seg_boxes = |pts: &[Point], segs: &[Seg]| -> Vec<[[f64; 2]; 3]> {
        segs.iter()
            .map(|s| {
                let (p, q) = (pts[s.ends[0]].bounds(), pts[s.ends[1]].bounds());
                [0, 1, 2].map(|k| [p[k][0].min(q[k][0]), p[k][1].max(q[k][1])])
            })
            .collect()
    };
    let boxes = seg_boxes(&pts, &segs);
    for (i, j) in overlapping_pairs(&boxes, &boxes, true) {
        let (s, t) = (segs[i], segs[j]);
        let [p, q] = s.ends;
        let [r, u] = t.ends;
        if p == r || p == u || q == r || q == u {
            continue;
        }
        let proper = orient(&pts, p, q, r).opposes(orient(&pts, p, q, u))
            && orient(&pts, r, u, p).opposes(orient(&pts, r, u, q));
        if !proper {
            continue;
        }
        let x = crossing(plane, axis, s.support, t.support)?;
        pts.push(x);
        gid.push(None);
    }
    // Crossings of three segments at one point, or at an existing point,
    // are one point.
    let map = dedupe(&pts);
    let mut remap = vec![usize::MAX; pts.len()];
    let mut merged_pts = Vec::new();
    let mut merged_gid = Vec::new();
    for i in 0..pts.len() {
        if map[i] as usize == i {
            remap[i] = merged_pts.len();
            merged_pts.push(pts[i]);
            merged_gid.push(gid[i]);
        }
    }
    for i in 0..pts.len() {
        remap[i] = remap[map[i] as usize];
    }
    let pts = merged_pts;
    let mut gid = merged_gid;
    for s in &mut segs {
        s.ends = s.ends.map(|e| remap[e]);
    }

    // 4. Split segments at the points on them.
    let pt_boxes: Vec<[[f64; 2]; 3]> = pts.iter().map(Point::bounds).collect();
    let boxes = seg_boxes(&pts, &segs);
    let mut on: Vec<Vec<usize>> = vec![Vec::new(); segs.len()];
    for (si, pi) in overlapping_pairs(&boxes, &pt_boxes, false) {
        let [p, q] = segs[si].ends;
        if pi == p || pi == q {
            continue;
        }
        if orient(&pts, p, q, pi) == Sign::Zero
            && implicit::dot_sign(&pts[p], &pts[pi], &pts[q]) == Sign::Negative
        {
            on[si].push(pi);
        }
    }
    let mut pieces: Vec<[usize; 2]> = Vec::new();
    for (s, inner) in segs.iter().zip(on.iter_mut()) {
        let [p, q] = s.ends;
        // Order along the segment: r before w when the angle at r between
        // p and w is straight (r lies between p and w).
        let mut ordered: Vec<usize> = Vec::with_capacity(inner.len());
        for &r in inner.iter() {
            let at = ordered
                .iter()
                .position(|&w| implicit::dot_sign(&pts[p], &pts[r], &pts[w]) == Sign::Negative)
                .unwrap_or(ordered.len());
            ordered.insert(at, r);
        }
        let mut chain = vec![p];
        chain.extend(ordered);
        chain.push(q);
        for w in chain.windows(2) {
            pieces.push([w[0].min(w[1]), w[0].max(w[1])]);
        }
    }
    pieces.sort_unstable();
    pieces.dedup();

    // 5. Triangulate, and hand out the pieces.
    let fail = |e: CdtError| ArrangementError::Triangulation(format!("{e:?}"));
    // Every point lies in a member, so the members' corners frame them.
    let frame: Vec<Point3> = members.iter().flat_map(|&m| corners(m)).collect();
    let mut cdt = Cdt::new(&pts, axis, &frame).map_err(fail)?;
    for &[p, q] in &pieces {
        cdt.constrain(p, q).map_err(fail)?;
    }
    // New points get global indices.
    for (i, g) in gid.iter_mut().enumerate() {
        if g.is_none() {
            *g = Some(points.len() as u32);
            points.push(pts[i]);
        }
    }
    let global = |i: usize| gid[i].unwrap_or(u32::MAX);
    let member_corners: Vec<[Point; 3]> = members
        .iter()
        .map(|&m| corners(m).map(Point::from))
        .collect();
    let member_sign: Vec<Sign> = member_corners
        .iter()
        .map(|c| implicit::orient2d_projected(&c[0], &c[1], &c[2], axis))
        .collect();
    for t in cdt.triangles() {
        for (mi, &m) in members.iter().enumerate() {
            let c = &member_corners[mi];
            let s = member_sign[mi];
            let inside = t.iter().all(|&v| {
                (0..3).all(|k| {
                    let o = implicit::orient2d_projected(&c[k], &c[(k + 1) % 3], &pts[v], axis);
                    o == s || o == Sign::Zero
                })
            });
            if !inside {
                continue;
            }
            let g = t.map(global);
            out_tris.push(if s == Sign::Positive {
                g
            } else {
                [g[0], g[2], g[1]]
            });
            origin.push(m as u32);
        }
    }
    Ok(())
}

/// Where two segments in the plane `plane` cross, from what they lie on.
fn crossing(
    plane: [Point3; 3],
    axis: usize,
    s: Support,
    t: Support,
) -> Result<Point, ArrangementError> {
    Ok(match (s, t) {
        (Support::Line(l), Support::Plane(p)) | (Support::Plane(p), Support::Line(l)) => {
            Point::lpi(l, p)?
        }
        (Support::Line(l), Support::Line(m)) => Point::lpi(l, vertical_plane(m, axis))?,
        (Support::Plane(p), Support::Plane(q)) => Point::tpi([plane, p, q])?,
    })
}

/// Pairs `(i, j)` of boxes from `a` and `b` that overlap (closed), by a
/// sweep along x; with `same`, `a` and `b` are one list and only `i < j`.
fn overlapping_pairs(a: &[[[f64; 2]; 3]], b: &[[[f64; 2]; 3]], same: bool) -> Vec<(usize, usize)> {
    let overlap = |p: &[[f64; 2]; 3], q: &[[f64; 2]; 3]| {
        (0..3).all(|k| p[k][0] <= q[k][1] && q[k][0] <= p[k][1])
    };
    // Events: every box by its low x, from a (tag 0) and b (tag 1).
    let mut events: Vec<(f64, u8, usize)> =
        a.iter().enumerate().map(|(i, x)| (x[0][0], 0, i)).collect();
    if !same {
        events.extend(b.iter().enumerate().map(|(i, x)| (x[0][0], 1, i)));
    }
    events.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
    let (mut active_a, mut active_b): (Vec<usize>, Vec<usize>) = (Vec::new(), Vec::new());
    let mut out = Vec::new();
    for (lo, tag, i) in events {
        active_a.retain(|&j| a[j][0][1] >= lo);
        active_b.retain(|&j| b[j][0][1] >= lo);
        if tag == 0 {
            let others = if same { &active_a } else { &active_b };
            for &j in others {
                if overlap(&a[i], &b[j]) {
                    out.push(if same { (i.min(j), i.max(j)) } else { (i, j) });
                }
            }
            active_a.push(i);
        } else {
            for &j in &active_a {
                if overlap(&a[j], &b[i]) {
                    out.push((j, i));
                }
            }
            active_b.push(i);
        }
    }
    out.sort_unstable();
    out
}
