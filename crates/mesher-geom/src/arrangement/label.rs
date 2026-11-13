//! Stage 3 of the arrangement (Milestone 16): keep the pieces between the
//! solid and the outside, by winding number.
//!
//! The input is closed and oriented (every edge used as often in each
//! direction), so its winding number is an integer off the surface, and
//! crossing a piece along its normal lowers it by one. The solid is where
//! it is at least 1: the union of everything the surface wraps.
//!
//! 1. Identical pieces (coplanar triangles' shared regions) become one face
//!    with a multiplicity: +1 per copy facing its way, −1 per copy facing
//!    the other. A face of multiplicity 0 lies between solids that touch
//!    face to face, with the same winding on both sides, and goes.
//! 2. Faces joined across an edge they alone share, in opposite directions,
//!    form a patch: no other surface meets the patch, so the winding on
//!    each side is the same all over it.
//! 3. Per patch, one exact count: from a point strictly inside one face (a
//!    barycentric point of its input triangle, checked exactly) to an
//!    explicit point far outside, where the winding is 0, each face the
//!    segment crosses adds or removes its multiplicity. Every test is an
//!    exact `orient3d`; a segment that touches an edge or a vertex, or lies
//!    in a face's plane, is cast again in another direction.
//! 4. A face is kept when exactly one side is solid, oriented so the solid
//!    is behind it.

use std::collections::BTreeMap;

use super::{Arrangement, ArrangementError, arrangement};
use crate::Point3;
use crate::implicit::{self, Point};
use crate::predicates::{Sign, plane_axis};

/// The boundary of the solid a closed, possibly self-intersecting surface
/// wraps, as triangles over exact points.
#[derive(Debug, Clone)]
pub struct Resolved {
    /// The input points (same indices), then the points the intersections
    /// made.
    pub points: Vec<Point>,
    /// Triangles of the solid's boundary, facing out.
    pub triangles: Vec<[u32; 3]>,
    /// The input triangle each lies on.
    pub origin: Vec<u32>,
    /// Intersecting pairs of input triangles, and of them coplanar ones.
    pub pairs: usize,
    pub coplanar_pairs: usize,
    /// Faces of the arrangement (identical pieces merged), patches, and
    /// faces dropped: inside the solid, outside it, or between solids that
    /// touch face to face.
    pub faces: usize,
    pub patches: usize,
    pub dropped: usize,
}

/// The boundary of the solid that the closed surface `triangles` wraps.
///
/// # Errors
/// As [`arrangement`]; [`ArrangementError::Open`] when the input is not
/// closed and oriented; [`ArrangementError::Cast`] if no segment could be
/// cast clear of every edge and vertex for a patch (a bug).
pub fn resolve(points: &[Point3], triangles: &[[u32; 3]]) -> Result<Resolved, ArrangementError> {
    resolve_with(points, triangles, true)
}

/// [`resolve`], casting for every patch (`cast_all`) or only for one patch
/// per connected piece of the arrangement and finding the rest around
/// edges, as for patches too thin to cast from (for tests: both must
/// agree).
///
/// # Errors
/// As [`resolve`].
#[doc(hidden)]
pub fn resolve_with(
    points: &[Point3],
    triangles: &[[u32; 3]],
    cast_all: bool,
) -> Result<Resolved, ArrangementError> {
    let unbalanced = unbalanced_edges(triangles);
    if unbalanced > 0 {
        return Err(ArrangementError::Open { edges: unbalanced });
    }
    let a = arrangement(points, triangles)?;
    label(points, triangles, a, cast_all)
}

/// Edges not used as often in each direction.
fn unbalanced_edges(triangles: &[[u32; 3]]) -> usize {
    let mut net: BTreeMap<(u32, u32), i64> = BTreeMap::new();
    for t in triangles {
        for k in 0..3 {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            *net.entry((u.min(v), u.max(v))).or_insert(0) += if u < v { 1 } else { -1 };
        }
    }
    net.values().filter(|&&n| n != 0).count()
}

/// A face: identical pieces merged.
struct Face {
    v: [u32; 3],
    m: i64,
    origin: u32,
}

fn label(
    input: &[Point3],
    triangles: &[[u32; 3]],
    a: Arrangement,
    cast_all: bool,
) -> Result<Resolved, ArrangementError> {
    // 1. Merge identical pieces.
    let mut by_key: BTreeMap<[u32; 3], usize> = BTreeMap::new();
    let mut faces: Vec<Face> = Vec::new();
    for (t, &o) in a.triangles.iter().zip(&a.origin) {
        let mut key = *t;
        key.sort_unstable();
        match by_key.get(&key) {
            Some(&f) => {
                let same = same_turn(faces[f].v, *t);
                faces[f].m += if same { 1 } else { -1 };
            }
            None => {
                by_key.insert(key, faces.len());
                faces.push(Face {
                    v: *t,
                    m: 1,
                    origin: o,
                });
            }
        }
    }
    let total_faces = faces.len();
    let mut dropped = faces.iter().filter(|f| f.m == 0).count();
    faces.retain(|f| f.m != 0);

    // 2. Patches.
    let mut edges: BTreeMap<(u32, u32), Vec<(usize, bool)>> = BTreeMap::new();
    for (i, f) in faces.iter().enumerate() {
        for k in 0..3 {
            let (u, v) = (f.v[k], f.v[(k + 1) % 3]);
            edges
                .entry((u.min(v), u.max(v)))
                .or_default()
                .push((i, u < v));
        }
    }
    let mut parent: Vec<usize> = (0..faces.len()).collect();
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
    for around in edges.values() {
        if let [(f, d), (g, e)] = around.as_slice()
            && d != e
            && faces[*f].m == faces[*g].m
        {
            let (x, y) = (find(&mut parent, *f), find(&mut parent, *g));
            if x != y {
                parent[x.max(y)] = x.min(y);
            }
        }
    }
    let mut patches: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..faces.len() {
        let r = find(&mut parent, i);
        patches.entry(r).or_default().push(i);
    }

    // 3. One cast per patch.
    let caster = Caster::new(input, triangles, &a.points, &faces);
    let mut patch_of = vec![0; faces.len()];
    for (&p, members) in &patches {
        for &f in members {
            patch_of[f] = p;
        }
    }
    // With `cast_all` off, one patch per connected piece (faces joined by
    // any edge).
    let mut piece: Vec<usize> = (0..faces.len()).collect();
    for around in edges.values() {
        for w in around.windows(2) {
            let (x, y) = (find(&mut piece, w[0].0), find(&mut piece, w[1].0));
            if x != y {
                piece[x.max(y)] = x.min(y);
            }
        }
    }
    let mut cast_in: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let mut wind: BTreeMap<usize, (i64, i64)> = BTreeMap::new();
    for (&p, members) in &patches {
        if !cast_all && !cast_in.insert(find(&mut piece, members[0])) {
            continue;
        }
        // The patch's largest face first; others if it is too thin to hold
        // a point.
        let mut order = members.clone();
        let area = |f: usize| caster.area(f);
        order.sort_by(|&x, &y| area(y).total_cmp(&area(x)).then(x.cmp(&y)));
        for &f in order.iter().take(8) {
            if let Some(s) = caster.sides(f)? {
                wind.insert(p, s);
                break;
            }
        }
    }
    // 3b. Patches no point fits in (needles an ulp wide): from a patch
    // around one of their edges, by the faces' order around it.
    while wind.len() < patches.len() {
        let mut progress = false;
        for (&(u, v), around) in &edges {
            let known = around.iter().any(|&(f, _)| wind.contains_key(&patch_of[f]));
            let unknown = around
                .iter()
                .any(|&(f, _)| !wind.contains_key(&patch_of[f]));
            if !(known && unknown) {
                continue;
            }
            let ring = radial(&a.points, &faces, u, v, around);
            // Wedge k lies after face ring[k], before face ring[k + 1].
            let n = ring.len();
            let start = (0..n)
                .find(|&k| wind.contains_key(&patch_of[ring[k].0]))
                .unwrap_or(0);
            let (f, forward) = ring[start];
            let (behind, front) = wind[&patch_of[f]];
            // A face holding u → v faces the way the ring turns.
            let mut after = if forward { front } else { behind };
            for step in 1..=n {
                let (g, forward) = ring[(start + step) % n];
                let m = faces[g].m;
                let before = after;
                after = if forward { before - m } else { before + m };
                let sides = if forward {
                    (before, after)
                } else {
                    (after, before)
                };
                match wind.get(&patch_of[g]) {
                    Some(&known) if known != sides => {
                        return Err(ArrangementError::Triangulation(
                            "winding numbers around an edge disagree".into(),
                        ));
                    }
                    Some(_) => {}
                    None => {
                        wind.insert(patch_of[g], sides);
                        progress = true;
                    }
                }
            }
        }
        if !progress {
            return Err(ArrangementError::Cast);
        }
    }
    let mut keep: Vec<Option<bool>> = vec![None; faces.len()];
    for (f, k) in keep.iter_mut().enumerate() {
        let (behind, front) = wind[&patch_of[f]];
        let (solid_behind, solid_front) = (behind >= 1, front >= 1);
        *k = (solid_behind != solid_front).then_some(solid_behind);
    }

    // 4. Keep the boundary, facing out.
    let mut out = Vec::new();
    let mut origin = Vec::new();
    for (f, k) in faces.iter().zip(&keep) {
        match k {
            Some(true) => out.push(f.v),
            Some(false) => out.push([f.v[0], f.v[2], f.v[1]]),
            None => {
                dropped += 1;
                continue;
            }
        }
        origin.push(f.origin);
    }
    Ok(Resolved {
        points: a.points,
        triangles: out,
        origin,
        pairs: a.pairs,
        coplanar_pairs: a.coplanar_pairs,
        faces: total_faces,
        patches: patches.len(),
        dropped,
    })
}

/// The faces around edge `uv` (`u < v`), each with whether it holds the
/// edge as `u → v`, in the order they turn about the axis `u → v`
/// (right-handed), from the first. Exact: each face's half-plane is placed
/// against the first's plane, then ordered within its half-turn.
fn radial(
    points: &[Point],
    faces: &[Face],
    u: u32,
    v: u32,
    around: &[(usize, bool)],
) -> Vec<(usize, bool)> {
    let (pu, pv) = (&points[u as usize], &points[v as usize]);
    let third = |f: usize| {
        let w = faces[f]
            .v
            .iter()
            .copied()
            .find(|&x| x != u && x != v)
            .unwrap_or(u);
        &points[w as usize]
    };
    let c0 = third(around[0].0);
    let axis = (0..3)
        .find(|&k| implicit::orient2d_projected(pu, pv, c0, k) != Sign::Zero)
        .unwrap_or(0);
    let side0 = implicit::orient2d_projected(pu, pv, c0, axis);
    // Half-turn classes: 0 the first's half-plane, 1 turning towards the
    // positive side of its plane, 2 the opposite half-plane, 3 the rest.
    let class = |f: usize| -> u8 {
        let c = third(f);
        match implicit::orient3d(pu, pv, c0, c) {
            Sign::Positive => 1,
            Sign::Negative => 3,
            Sign::Zero => {
                if implicit::orient2d_projected(pu, pv, c, axis) == side0 {
                    0
                } else {
                    2
                }
            }
        }
    };
    let mut ring: Vec<((usize, bool), u8)> = around.iter().map(|&x| (x, class(x.0))).collect();
    // Within a class, f before g when g lies further along the turn.
    let mut sorted: Vec<((usize, bool), u8)> = Vec::with_capacity(ring.len());
    ring.sort_by_key(|x| x.1);
    for item in ring {
        let at = sorted
            .iter()
            .position(|&(other, c)| {
                c > item.1
                    || (c == item.1
                        && implicit::orient3d(pu, pv, third(item.0.0), third(other.0))
                            == Sign::Positive)
            })
            .unwrap_or(sorted.len());
        sorted.insert(at, item);
    }
    sorted.into_iter().map(|(x, _)| x).collect()
}

/// Do `a` and `b`, the same three vertices, turn the same way?
fn same_turn(a: [u32; 3], b: [u32; 3]) -> bool {
    let k = b.iter().position(|&v| v == a[0]).unwrap_or(0);
    b[(k + 1) % 3] == a[1]
}

/// Exact segment casts against the faces.
struct Caster<'a> {
    input: &'a [Point3],
    triangles: &'a [[u32; 3]],
    points: &'a [Point],
    faces: &'a [Face],
    /// Each face's corners, and a box certainly holding it.
    corners: Vec<[Point; 3]>,
    boxes: Vec<[[f64; 2]; 3]>,
    /// The centre and radius of a sphere around everything.
    centre: Point3,
    radius: f64,
}

/// Directions to cast in, spread and away from the axes.
const DIRECTIONS: [Point3; 12] = [
    [0.8017, 0.5345, 0.2673],
    [-0.4364, 0.8729, 0.2182],
    [0.3015, -0.3015, 0.9045],
    [-0.6667, -0.3333, -0.6667],
    [0.2357, 0.9428, -0.2357],
    [0.9113, -0.2735, -0.3077],
    [-0.1826, 0.3651, -0.9129],
    [-0.8452, -0.5071, 0.1690],
    [0.4851, 0.1213, -0.8660],
    [-0.2408, -0.9631, 0.1204],
    [0.6247, -0.7063, 0.3313],
    [-0.5774, 0.5774, 0.5774],
];

impl<'a> Caster<'a> {
    fn new(
        input: &'a [Point3],
        triangles: &'a [[u32; 3]],
        points: &'a [Point],
        faces: &'a [Face],
    ) -> Caster<'a> {
        let corners: Vec<[Point; 3]> = faces
            .iter()
            .map(|f| f.v.map(|v| points[v as usize]))
            .collect();
        let boxes = corners
            .iter()
            .map(|c| {
                let b = c.map(|p| p.bounds());
                [0, 1, 2].map(|k| {
                    [
                        b.iter().map(|x| x[k][0]).fold(f64::INFINITY, f64::min),
                        b.iter().map(|x| x[k][1]).fold(f64::NEG_INFINITY, f64::max),
                    ]
                })
            })
            .collect();
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in input {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let centre = [0, 1, 2].map(|k| 0.5 * (lo[k] + hi[k]));
        let radius = (0..3)
            .map(|k| (hi[k] - lo[k]).powi(2))
            .sum::<f64>()
            .sqrt()
            .max(1.0);
        Caster {
            input,
            triangles,
            points,
            faces,
            corners,
            boxes,
            centre,
            radius,
        }
    }

    fn area(&self, f: usize) -> f64 {
        let [a, b, c] = self.faces[f].v.map(|v| self.points[v as usize].approx());
        let u = [0, 1, 2].map(|k| b[k] - a[k]);
        let w = [0, 1, 2].map(|k| c[k] - a[k]);
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
    }

    /// A point strictly inside face `f`: a barycentric point of its input
    /// triangle near the face's centroid, checked exactly.
    fn inside(&self, f: usize) -> Option<Point> {
        let face = &self.faces[f];
        let tri = self.triangles[face.origin as usize].map(|v| self.input[v as usize]);
        let axis = plane_axis(tri[0], tri[1], tri[2])?;
        let c = &self.corners[f];
        let s = implicit::orient2d_projected(&c[0], &c[1], &c[2], axis);
        let approx = c.map(|p| p.approx());
        // Barycentric weights in the projection.
        let [i, j] = [(axis + 1) % 3, (axis + 2) % 3];
        let area2 = |p: Point3, q: Point3, r: Point3| {
            (q[i] - p[i]) * (r[j] - p[j]) - (q[j] - p[j]) * (r[i] - p[i])
        };
        // The centroid, then points nearer each corner's opposite edge.
        for mix in [
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 1.0],
            [2.0, 1.0, 2.0],
            [1.0, 2.0, 2.0],
        ] {
            let total: f64 = mix.iter().sum();
            let target = [0, 1, 2].map(|k| {
                (mix[0] * approx[0][k] + mix[1] * approx[1][k] + mix[2] * approx[2][k]) / total
            });
            let weights = [
                area2(target, tri[1], tri[2]),
                area2(tri[0], target, tri[2]),
                area2(tri[0], tri[1], target),
            ];
            let Ok(q) = Point::bpt(tri, weights) else {
                continue;
            };
            let strictly =
                (0..3).all(|k| implicit::orient2d_projected(&c[k], &c[(k + 1) % 3], &q, axis) == s);
            if strictly {
                return Some(q);
            }
        }
        None
    }

    /// The winding behind and in front of face `f`, or `None` if no point
    /// strictly inside it could be found.
    fn sides(&self, f: usize) -> Result<Option<(i64, i64)>, ArrangementError> {
        let Some(q) = self.inside(f) else {
            return Ok(None);
        };
        let qa = q.approx();
        for d in DIRECTIONS {
            let o = [0, 1, 2].map(|k| self.centre[k] + 4.0 * self.radius * d[k]);
            let Some(count) = self.count(f, o, &q, qa) else {
                continue;
            };
            let m = self.faces[f].m;
            let c = &self.corners[f];
            let op = Point::from(o);
            return Ok(Some(match implicit::orient3d(&c[0], &c[1], &c[2], &op) {
                // o in front: the count is the winding in front.
                Sign::Positive => (count + m, count),
                Sign::Negative => (count, count - m),
                Sign::Zero => continue,
            }));
        }
        Err(ArrangementError::Cast)
    }

    /// The winding at the end `q` of the segment from `o` (winding 0) on
    /// `o`'s side of face `skip`: the multiplicities of the faces crossed,
    /// entering (from front to back) counting up. `None` if the segment
    /// touches an edge or a vertex, or lies in a face's plane.
    fn count(&self, skip: usize, o: Point3, q: &Point, qa: Point3) -> Option<i64> {
        let op = Point::from(o);
        let mut w = 0;
        for (g, c) in self.corners.iter().enumerate() {
            if g == skip || !self.might_cross(g, o, qa) {
                continue;
            }
            let so = implicit::orient3d(&c[0], &c[1], &c[2], &op);
            let sq = implicit::orient3d(&c[0], &c[1], &c[2], q);
            if so == Sign::Zero {
                return None;
            }
            // q is inside face `skip` only, so a plane through it meets the
            // segment only there, outside that face.
            if sq == Sign::Zero || so == sq {
                continue;
            }
            let e = [0, 1, 2].map(|k| implicit::orient3d(&op, q, &c[k], &c[(k + 1) % 3]));
            if e.contains(&Sign::Zero) {
                return None;
            }
            if e[0] == e[1] && e[1] == e[2] {
                let m = self.faces[g].m;
                w += if so == Sign::Positive { m } else { -m };
            }
        }
        Some(w)
    }

    /// Might the segment from `o` to about `q` meet face `g`'s box? A slab
    /// test on the box grown a little, so rounding never hides a crossing.
    fn might_cross(&self, g: usize, o: Point3, q: Point3) -> bool {
        let b = &self.boxes[g];
        if b.iter().flatten().any(|x| !x.is_finite()) {
            return true;
        }
        let (mut t0, mut t1) = (0.0f64, 1.0f64);
        for k in 0..3 {
            let grow = 1e-9 * (self.radius + b[k][0].abs().max(b[k][1].abs()));
            let (lo, hi) = (b[k][0] - grow, b[k][1] + grow);
            let d = q[k] - o[k];
            if d == 0.0 {
                if o[k] < lo || o[k] > hi {
                    return false;
                }
                continue;
            }
            let (mut a, mut z) = ((lo - o[k]) / d, (hi - o[k]) / d);
            if a > z {
                std::mem::swap(&mut a, &mut z);
            }
            t0 = t0.max(a - 1e-9);
            t1 = t1.min(z + 1e-9);
            if t0 > t1 {
                return false;
            }
        }
        true
    }
}
