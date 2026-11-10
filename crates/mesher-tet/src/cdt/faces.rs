//! Face recovery: make every input triangle (a facet) a union of faces of
//! the tetrahedralization, after segment recovery.
//!
//! Segment recovery leaves Steiner points on triangle edges, so each facet
//! is a polygon: its corners and the points along its sides. A facet with
//! no points inside is triangulated into subfaces by a dynamic program over
//! the polygon's triangulations that keeps the most faces the
//! tetrahedralization already has, then the most edges, then the best worst
//! shape. No triangle has three vertices on one side; the points lie
//! exactly on the sides (LNC points), so that is decided by which sides each
//! vertex is on, without predicates.
//!
//! Recovery goes in three stages, each only where the last left a facet
//! crossed by an edge:
//!
//! 1. Flips: the subfaces' sides (the diagonals of a facet), then each
//!    subface, by removing the edges that cross it
//!    ([`ExactDelaunay::recover_face`]). Passes repeat while they make
//!    progress; before each, a facet with a missing subface is triangulated
//!    again to fit the tetrahedralization as it is then.
//! 2. The cavity of the tets crossing the facet, retetrahedralized above and
//!    below it ([`ExactDelaunay::recover_facet_by_cavity`]); the cavity
//!    chooses how the uncovered part is triangulated. Where growing the
//!    cavity fails, it is filled by gift-wrapping instead.
//! 3. Steiner points on the facet where edges cross it: the exact
//!    line-plane intersection (an LPI point) when the edge's ends are
//!    explicit, else a BPT point of the input triangle near it, exactly on
//!    its plane. They go in so that no protected edge or face is lost; then
//!    the cavity again. Rounds repeat within the point budget, a facet
//!    taking one point in its first round and twice as many each round it
//!    is still crossed (up to eight), at edges sharing no vertex, so the
//!    points spread over it.
//!
//! Every subsegment, and every side of a subface or of a face already on a
//! facet, is protected, so nothing recovered is lost. A facet is recovered
//! when no edge crosses it; its subfaces are then read off the
//! tetrahedralization. What is left is reported, not repaired.

use std::collections::BTreeSet;

use mesher_geom::Point3;
use mesher_geom::implicit::{self, Kind, Point};

use super::delaunay::{ExactDelaunay, Inserted, Refill};
use super::segments::{SegmentError, Subsegment, recover_segments};

/// A face of the tetrahedralization on an input triangle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Subface {
    /// The input triangle.
    pub triangle: u32,
    pub v: [u32; 3],
}

/// Why face recovery failed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FaceError {
    Segments(SegmentError),
    /// An input triangle's corners are collinear.
    Degenerate {
        triangle: u32,
    },
}

impl std::fmt::Display for FaceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FaceError::Segments(e) => write!(f, "{e}"),
            FaceError::Degenerate { triangle } => {
                write!(f, "triangle {triangle} is degenerate")
            }
        }
    }
}

impl From<SegmentError> for FaceError {
    fn from(e: SegmentError) -> FaceError {
        FaceError::Segments(e)
    }
}

/// A tetrahedralization with every segment recovered and, unless listed in
/// `missing`, every input triangle a union of faces.
#[derive(Debug, Clone)]
pub struct RecoveredFaces {
    pub dt: ExactDelaunay,
    pub subsegments: Vec<Subsegment>,
    /// The faces on each recovered input triangle, by triangle.
    pub subfaces: Vec<Subface>,
    /// Input triangles still crossed by an edge.
    pub missing: Vec<u32>,
    /// Steiner points on segments.
    pub steiner: usize,
    /// Steiner points inside input triangles.
    pub steiner_on_faces: usize,
    /// Subfaces missing after segment recovery.
    pub initially_missing: usize,
    /// Input triangles crossed after segment recovery, after flips, and
    /// after the first cavity retetrahedralizations.
    pub crossed: [usize; 3],
    /// Facet cavities filled by gift-wrapping.
    pub wrapped: usize,
    /// Facet cavity expansions refused for reaching across the facet's
    /// plane (see [`ExactDelaunay::expansions_across_plane`]).
    pub across_plane: usize,
}

/// Flip steps per missing subface or diagonal in a pass.
const FLIP_STEPS: usize = 64;

/// Passes over the missing subfaces while they make progress.
const MAX_PASSES: usize = 8;

/// Rounds of Steiner points on facets.
const MAX_ROUNDS: usize = 64;

/// Steiner points inserted per facet in a round, at most: one at first,
/// doubling each round the facet is still crossed, up to this; and the
/// crossing edges tried for them.
const PER_ROUND: usize = 8;
const TRIES: usize = 4 * PER_ROUND;

/// A facet: its input triangle, boundary polygon (with the sides each
/// vertex is on) and the Steiner points inside it.
struct Facet {
    corners: [u32; 3],
    polygon: Vec<(u32, [bool; 3])>,
    inside: Vec<u32>,
    /// The subfaces flips work towards, while there are no points inside.
    subfaces: Vec<[u32; 3]>,
}

impl Facet {
    fn vertices(&self) -> Vec<u32> {
        self.polygon
            .iter()
            .map(|p| p.0)
            .chain(self.inside.iter().copied())
            .collect()
    }
}

/// The faces of `dt` with every vertex in `on` (the vertices of one input
/// triangle, all on it), so on the triangle.
fn faces_on(dt: &ExactDelaunay, on: &[u32]) -> Vec<[u32; 3]> {
    let set: BTreeSet<u32> = on.iter().copied().collect();
    let mut out: Vec<[u32; 3]> = on
        .iter()
        .flat_map(|&v| dt.faces_at(v))
        .filter(|f| f.iter().all(|w| set.contains(w)))
        .collect();
    out.sort_unstable();
    out.dedup();
    out
}

fn key(u: u32, v: u32) -> (u32, u32) {
    (u.min(v), u.max(v))
}

/// Recover every segment of `triangles` (see [`recover_segments`]), then
/// make every input triangle a union of faces, by flips, cavity
/// retetrahedralization and Steiner points on it. `max_points` bounds the
/// vertex count (0: 20 × the inputs + 10,000).
///
/// # Errors
/// See [`FaceError`].
pub fn recover_faces(
    points: &[Point3],
    triangles: &[[u32; 3]],
    max_points: usize,
    seed: u64,
) -> Result<RecoveredFaces, FaceError> {
    let segments = recover_segments(points, triangles, max_points, seed)?;
    let limit = if max_points == 0 {
        20 * points.len() + 10_000
    } else {
        max_points
    };
    let mut dt = segments.dt;
    let subs = segments.subsegments;
    // Each segment's vertices from its first end to its second.
    let mut chains: std::collections::BTreeMap<(u32, u32), Vec<u32>> =
        std::collections::BTreeMap::new();
    for s in &subs {
        let chain = chains.entry(s.segment).or_default();
        if chain.is_empty() {
            chain.push(s.ends[0]);
        }
        chain.push(s.ends[1]);
    }
    let mut facets: Vec<Facet> = Vec::with_capacity(triangles.len());
    for (ti, t) in triangles.iter().enumerate() {
        let ti = ti as u32;
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            return Err(FaceError::Degenerate { triangle: ti });
        }
        let [p0, p1, p2] = t.map(|v| dt.points()[v as usize]);
        if implicit::collinear(&p0, &p1, &p2) {
            return Err(FaceError::Degenerate { triangle: ti });
        }
        // The polygon, and the sides (0: t0 t1, 1: t1 t2, 2: t2 t0) each
        // vertex lies on.
        let mut polygon: Vec<(u32, [bool; 3])> = Vec::new();
        for k in 0..3 {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            let mut chain = chains
                .get(&key(u, v))
                .cloned()
                .unwrap_or_else(|| vec![u, v]);
            if chain.first() != Some(&u) {
                chain.reverse();
            }
            let mut on = [false; 3];
            on[k] = true;
            for &w in &chain[..chain.len() - 1] {
                polygon.push((w, on));
            }
        }
        // Each corner is also on the side before it.
        for (k, &corner) in t.iter().enumerate() {
            if let Some(c) = polygon.iter_mut().find(|(w, _)| *w == corner) {
                c.1[(k + 2) % 3] = true;
            }
        }
        let subfaces = triangulate(&dt, &polygon).ok_or(FaceError::Degenerate { triangle: ti })?;
        facets.push(Facet {
            corners: *t,
            polygon,
            inside: Vec::new(),
            subfaces,
        });
    }
    let initially_missing = facets
        .iter()
        .flat_map(|f| f.subfaces.iter())
        .filter(|&&f| !dt.has_face(f))
        .count();
    let segment_edges: BTreeSet<(u32, u32)> =
        subs.iter().map(|s| key(s.ends[0], s.ends[1])).collect();
    // Subsegments, the subfaces' sides, and the sides of faces already on
    // a facet with points inside.
    let protect = |dt: &ExactDelaunay, facets: &[Facet]| {
        let mut p = segment_edges.clone();
        for f in facets {
            let faces = if f.inside.is_empty() {
                f.subfaces.clone()
            } else {
                faces_on(dt, &f.vertices())
            };
            for s in faces {
                for k in 0..3 {
                    p.insert(key(s[k], s[(k + 1) % 3]));
                }
            }
        }
        p
    };
    // A facet with no points inside whose subfaces are all faces is covered;
    // otherwise look for an edge meeting it.
    let crossed = |dt: &ExactDelaunay, f: &Facet| {
        !(f.inside.is_empty() && f.subfaces.iter().all(|&s| dt.has_face(s)))
            && !dt.facet_crossings(f.corners, &f.vertices()).is_empty()
    };
    let crossed_before = facets.iter().filter(|f| crossed(&dt, f)).count();
    // 1. Flips.
    for pass in 0..MAX_PASSES {
        if pass > 0 {
            for f in &mut facets {
                if f.subfaces.iter().any(|&s| !dt.has_face(s)) {
                    if let Some(s) = triangulate(&dt, &f.polygon) {
                        f.subfaces = s;
                    }
                }
            }
        }
        let protected = protect(&dt, &facets);
        let is_protected = |x: u32, y: u32| protected.contains(&key(x, y));
        let mut progress = false;
        for f in facets.iter().flat_map(|f| f.subfaces.iter()) {
            if dt.has_face(*f) {
                continue;
            }
            for k in 0..3 {
                let (u, v) = (f[k], f[(k + 1) % 3]);
                if !dt.has_edge(u, v) {
                    let _s = crate::prof::scope("cdt faces: diagonals");
                    progress |= dt.recover_edge(u, v, &is_protected, FLIP_STEPS);
                }
            }
        }
        for f in facets.iter().flat_map(|f| f.subfaces.iter()) {
            if !dt.has_face(*f) {
                let _s = crate::prof::scope("cdt faces: flips");
                progress |= dt.recover_face(*f, &is_protected, FLIP_STEPS);
            }
        }
        if !progress {
            break;
        }
    }
    // A recovered facet stays recovered, its faces protected, so only these
    // are checked again (every facet is, at the end).
    let mut active: Vec<usize> = (0..facets.len())
        .filter(|&ti| crossed(&dt, &facets[ti]))
        .collect();
    let crossed_after_flips = active.len();
    // 2. Cavities, then 3. Steiner points on what is left.
    let mut steiner_on_faces = 0;
    let mut wrapped = 0;
    // Rounds each facet has had Steiner points in.
    let mut steiner_rounds = vec![0u32; facets.len()];
    let mut crossed_after_cavities = None;
    for _ in 0..MAX_ROUNDS {
        let mut left = Vec::new();
        for &ti in &active {
            if !crossed(&dt, &facets[ti]) {
                continue;
            }
            let protected = protect(&dt, &facets);
            let is_protected = |x: u32, y: u32| protected.contains(&key(x, y));
            let f = &facets[ti];
            let _s = crate::prof::scope("cdt faces: cavity");
            let refill = dt.recover_facet_by_cavity(f.corners, &f.vertices(), &is_protected);
            if refill == Refill::GiftWrapped {
                wrapped += 1;
            }
            if refill.recovered() {
                // The cavity chose how the uncovered part is triangulated:
                // read it back, so the faces it made are protected.
                if facets[ti].inside.is_empty() {
                    if let Some(s) = triangulate(&dt, &facets[ti].polygon) {
                        facets[ti].subfaces = s;
                    }
                }
            } else {
                left.push(ti);
            }
        }
        crossed_after_cavities.get_or_insert(left.len());
        if left.is_empty() {
            break;
        }
        let mut inserted = false;
        active.clone_from(&left);
        for ti in left {
            if dt.points().len() >= limit {
                break;
            }
            // At or near the first crossings that give new vertices.
            let crossings = dt.facet_crossings(facets[ti].corners, &facets[ti].vertices());
            let protected = protect(&dt, &facets);
            let is_protected = |x: u32, y: u32| protected.contains(&key(x, y));
            let _s = crate::prof::scope("cdt faces: steiner points");
            // Edges sharing no vertex first, so the points spread over the
            // facet rather than along one fan of edges; the rest after.
            let mut used: BTreeSet<u32> = BTreeSet::new();
            let (spread, rest): (Vec<_>, Vec<_>) = crossings.iter().partition(|&&(x, y)| {
                let fresh = !used.contains(&x) && !used.contains(&y);
                if fresh {
                    used.insert(x);
                    used.insert(y);
                }
                fresh
            });
            let candidates: Vec<Point> = spread
                .iter()
                .chain(&rest)
                .flat_map(|&(x, y)| crossing_point(&dt, points, facets[ti].corners, x, y))
                .take(TRIES)
                .collect();
            let mut added = 0;
            let most = PER_ROUND.min(1 << steiner_rounds[ti].min(8));
            steiner_rounds[ti] += 1;
            for p in candidates {
                if added == most || dt.points().len() >= limit {
                    break;
                }
                if let Ok(Inserted::New(m)) = dt.insert_protected(p, &is_protected) {
                    facets[ti].inside.push(m);
                    steiner_on_faces += 1;
                    inserted = true;
                    added += 1;
                }
            }
        }
        if !inserted {
            break;
        }
    }
    let mut missing = Vec::new();
    let mut subfaces = Vec::new();
    for (ti, f) in facets.iter().enumerate() {
        if crossed(&dt, f) {
            missing.push(ti as u32);
        } else {
            subfaces.extend(faces_on(&dt, &f.vertices()).into_iter().map(|v| Subface {
                triangle: ti as u32,
                v,
            }));
        }
    }
    let across_plane = dt.expansions_across_plane();
    Ok(RecoveredFaces {
        dt,
        subsegments: subs,
        subfaces,
        missing,
        steiner: segments.steiner,
        steiner_on_faces,
        initially_missing,
        crossed: [
            crossed_before,
            crossed_after_flips,
            crossed_after_cavities.unwrap_or(0),
        ],
        wrapped,
        across_plane,
    })
}

/// A point exactly on input triangle `corners` where the edge `xy` crosses
/// it, to try as a Steiner point: exactly the crossing (an LPI point) when
/// both ends are explicit, which lies inside the part of the triangle not
/// covered by faces. Otherwise near it, from approximate positions: a BPT
/// point of the triangle, exactly on its plane, with weights kept a little
/// inside.
fn crossing_point(
    dt: &ExactDelaunay,
    points: &[Point3],
    corners: [u32; 3],
    x: u32,
    y: u32,
) -> Option<Point> {
    let c = corners.map(|v| points[v as usize]);
    if let (Kind::Explicit(a), Kind::Explicit(b)) = (
        dt.points()[x as usize].kind(),
        dt.points()[y as usize].kind(),
    ) {
        if let Ok(p) = Point::lpi([*a, *b], c) {
            return Some(p);
        }
    }
    let (px, py) = (
        dt.points()[x as usize].approx(),
        dt.points()[y as usize].approx(),
    );
    let sub = |a: Point3, b: Point3| [0, 1, 2].map(|k| a[k] - b[k]);
    let cross = |a: Point3, b: Point3| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let dot = |a: Point3, b: Point3| (0..3).map(|k| a[k] * b[k]).sum::<f64>();
    let n = cross(sub(c[1], c[0]), sub(c[2], c[0]));
    let (dx, dy) = (dot(n, sub(px, c[0])), dot(n, sub(py, c[0])));
    let s = dx / (dx - dy);
    if !s.is_finite() {
        return None;
    }
    let q = [0, 1, 2].map(|k| px[k] + s * (py[k] - px[k]));
    let area = dot(n, n);
    let mut w =
        [0, 1, 2].map(|k| dot(n, cross(sub(c[(k + 1) % 3], q), sub(c[(k + 2) % 3], q))) / area);
    // Keep it inside, away from the sides.
    for x in &mut w {
        *x = x.max(1e-3);
    }
    let sum: f64 = w.iter().sum();
    if !(sum.is_finite() && sum > 0.0) {
        return None;
    }
    Point::bpt(c, w.map(|x| x / sum)).ok()
}

/// Triangulate a facet polygon (vertices in order, with the sides each is
/// on): of the triangulations with no triangle on one side, the one with the
/// most existing faces, then the best worst shape. `None` if there is none
/// (the corners are collinear).
fn triangulate(dt: &ExactDelaunay, polygon: &[(u32, [bool; 3])]) -> Option<Vec<[u32; 3]>> {
    let k = polygon.len();
    if k == 3 {
        return Some(vec![[polygon[0].0, polygon[1].0, polygon[2].0]]);
    }
    let mut existing: BTreeSet<[u32; 3]> = BTreeSet::new();
    let ids: BTreeSet<u32> = polygon.iter().map(|p| p.0).collect();
    for &(v, _) in polygon {
        for f in dt.faces_at(v) {
            if f.iter().all(|w| ids.contains(w)) {
                existing.insert(f);
            }
        }
    }
    let mut edges: BTreeSet<(u32, u32)> = BTreeSet::new();
    for &(v, _) in polygon {
        for w in dt.neighbours(v) {
            if ids.contains(&w) {
                edges.insert((v.min(w), v.max(w)));
            }
        }
    }
    let at: Vec<Point3> = polygon
        .iter()
        .map(|p| dt.points()[p.0 as usize].approx())
        .collect();
    let flat = |i: usize, m: usize, j: usize| {
        (0..3).any(|s| polygon[i].1[s] && polygon[m].1[s] && polygon[j].1[s])
    };
    let shape = |i: usize, m: usize, j: usize| {
        let (a, b, c) = (at[i], at[m], at[j]);
        let u = [0, 1, 2].map(|q| b[q] - a[q]);
        let w = [0, 1, 2].map(|q| c[q] - a[q]);
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let area2 = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let d2 = |p: Point3, q: Point3| (0..3).map(|r| (p[r] - q[r]).powi(2)).sum::<f64>();
        let sum = d2(a, b) + d2(b, c) + d2(c, a);
        let s = area2 / sum;
        if s.is_finite() { s } else { 0.0 }
    };
    // best[i][j]: (existing faces, existing edges, worst shape) of the chain
    // i..=j closed by the edge j i, and its middle vertex.
    const NONE: usize = usize::MAX;
    let mut best = vec![vec![(0usize, 0usize, f64::INFINITY, NONE); k]; k];
    for len in 2..k {
        for i in 0..k - len {
            let j = i + len;
            let mut cell = (0usize, 0usize, f64::NEG_INFINITY, NONE);
            for m in i + 1..j {
                let (l, r) = (best[i][m], best[m][j]);
                if (m > i + 1 && l.3 == NONE) || (j > m + 1 && r.3 == NONE) || flat(i, m, j) {
                    continue;
                }
                let mut tri = [polygon[i].0, polygon[m].0, polygon[j].0];
                tri.sort_unstable();
                let count = l.0 + r.0 + usize::from(existing.contains(&tri));
                let (u, w) = (polygon[i].0, polygon[j].0);
                let edge = l.1 + r.1 + usize::from(edges.contains(&(u.min(w), u.max(w))));
                let worst = l.2.min(r.2).min(shape(i, m, j));
                if cell.3 == NONE || (count, edge, worst) > (cell.0, cell.1, cell.2) {
                    cell = (count, edge, worst, m);
                }
            }
            best[i][j] = cell;
        }
    }
    let mut out = Vec::with_capacity(k - 2);
    let mut stack = vec![(0, k - 1)];
    while let Some((i, j)) = stack.pop() {
        if j < i + 2 {
            continue;
        }
        let m = best[i][j].3;
        if m == NONE {
            return None;
        }
        out.push([polygon[i].0, polygon[m].0, polygon[j].0]);
        stack.push((i, m));
        stack.push((m, j));
    }
    Some(out)
}
