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

use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, BTreeSet};

use mesher_geom::implicit::{self, Kind, Point};
use mesher_geom::{Point3, Sign};

use super::delaunay::{ExactDelaunay, Inserted, Refill};
use super::segments::{SegmentError, Subsegment, recover_segments_with};

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
    /// Each input triangle's facet, by its first triangle (itself unless
    /// flat triangles were merged).
    pub facet_of: Vec<u32>,
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

/// With `merge_flat`, the longest polygon a lone triangle is triangulated
/// into subfaces for; a longer one is a planar region instead.
const SUBFACE_POLYGON: usize = 512;

/// A facet: its input triangle, boundary polygon (with the sides each
/// vertex is on) and the Steiner points inside it.
struct Facet {
    /// Its input triangle (for a merged facet, its first).
    triangle: u32,
    corners: [u32; 3],
    polygon: Vec<(u32, [bool; 3])>,
    inside: Vec<u32>,
    /// The subfaces flips work towards, while there are no points inside.
    subfaces: Vec<[u32; 3]>,
    /// A merged facet's triangles, with their flat sides (empty for one
    /// input triangle): `polygon` then holds every vertex on it, boundary
    /// and inside.
    tris: Vec<([u32; 3], [bool; 3])>,
    /// A merged facet's boundary pieces, each in the direction that has the
    /// facet on its left.
    boundary: Vec<(u32, u32)>,
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

/// The edges between two coplanar input triangles that fold out flat (each
/// side of the edge has one triangle, in opposite directions, and their far
/// corners lie on opposite sides of it in their plane). Exact.
fn flat_edges(points: &[Point], triangles: &[[u32; 3]]) -> BTreeSet<(u32, u32)> {
    // Each edge's sides: (from, to, far corner) in its triangles.
    type Sides = Vec<(u32, u32, u32)>;
    let mut by_edge: std::collections::BTreeMap<(u32, u32), Sides> =
        std::collections::BTreeMap::new();
    for t in triangles {
        for k in 0..3 {
            let (u, v, w) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            by_edge.entry(key(u, v)).or_default().push((u, v, w));
        }
    }
    let pt = |v: u32| &points[v as usize];
    let mut out = BTreeSet::new();
    for (e, sides) in by_edge {
        let &[(u, v, a), (u2, v2, b)] = sides.as_slice() else {
            continue;
        };
        if (u2, v2) != (v, u) || implicit::orient3d(pt(u), pt(v), pt(a), pt(b)) != Sign::Zero {
            continue;
        }
        let Some(axis) =
            (0..3).find(|&k| implicit::orient2d_projected(pt(u), pt(v), pt(a), k) != Sign::Zero)
        else {
            continue;
        };
        let sa = implicit::orient2d_projected(pt(u), pt(v), pt(a), axis);
        let sb = implicit::orient2d_projected(pt(u), pt(v), pt(b), axis);
        if sb != Sign::Zero && sb != sa {
            out.insert(e);
        }
    }
    out
}

/// The faces of a merged facet: those with every vertex on it, flooded from
/// the inner side of its boundary pieces without crossing one (faces with
/// their vertices on it can span a notch outside it).
fn region_faces(dt: &ExactDelaunay, f: &Facet) -> Vec<[u32; 3]> {
    region_flood(dt, f).0
}

/// [`region_faces`], and whether they cover the facet: every boundary piece
/// has one of them on its inner side and every other edge of theirs two of
/// them. Then they form a patch whose boundary is exactly the facet's, so
/// they cover it (faces of a tetrahedralization do not overlap); where a
/// part is uncovered, the faces around it have an edge with one of them.
/// Exact and combinatorial.
fn region_flood(dt: &ExactDelaunay, f: &Facet) -> (Vec<[u32; 3]>, bool) {
    let _s = crate::prof::scope("cdt faces: region faces");
    let candidates = faces_on(dt, &f.vertices());
    let Some(&(t0, _)) = f.tris.first() else {
        return (candidates, false);
    };
    let pt = |v: u32| &dt.points()[v as usize];
    let Some(axis) = (0..3)
        .find(|&k| implicit::orient2d_projected(pt(t0[0]), pt(t0[1]), pt(t0[2]), k) != Sign::Zero)
    else {
        return (Vec::new(), false);
    };
    let up = implicit::orient2d_projected(pt(t0[0]), pt(t0[1]), pt(t0[2]), axis);
    let mut by_edge: std::collections::BTreeMap<(u32, u32), Vec<usize>> =
        std::collections::BTreeMap::new();
    for (i, c) in candidates.iter().enumerate() {
        for k in 0..3 {
            by_edge
                .entry(key(c[k], c[(k + 1) % 3]))
                .or_default()
                .push(i);
        }
    }
    let walls: BTreeSet<(u32, u32)> = f.boundary.iter().map(|&(a, b)| key(a, b)).collect();
    let mut inside = vec![false; candidates.len()];
    let mut stack = Vec::new();
    let mut covered = true;
    for &(a, b) in &f.boundary {
        let mut inner = 0;
        for &i in by_edge.get(&key(a, b)).map_or(&[][..], |v| v.as_slice()) {
            let c = candidates[i];
            let Some(&w) = c.iter().find(|&&x| x != a && x != b) else {
                continue;
            };
            if implicit::orient2d_projected(pt(a), pt(b), pt(w), axis) == up {
                inner += 1;
                if !inside[i] {
                    inside[i] = true;
                    stack.push(i);
                }
            }
        }
        covered &= inner == 1;
    }
    while let Some(i) = stack.pop() {
        let c = candidates[i];
        for k in 0..3 {
            let e = key(c[k], c[(k + 1) % 3]);
            if walls.contains(&e) {
                continue;
            }
            for &j in by_edge.get(&e).map_or(&[][..], |v| v.as_slice()) {
                if !inside[j] {
                    inside[j] = true;
                    stack.push(j);
                }
            }
        }
    }
    let faces: Vec<[u32; 3]> = candidates
        .into_iter()
        .zip(inside)
        .filter(|&(_, i)| i)
        .map(|(c, _)| c)
        .collect();
    if covered {
        let mut count: std::collections::BTreeMap<(u32, u32), u32> =
            std::collections::BTreeMap::new();
        for c in &faces {
            for k in 0..3 {
                *count.entry(key(c[k], c[(k + 1) % 3])).or_default() += 1;
            }
        }
        covered = count
            .iter()
            .all(|(e, &n)| n == if walls.contains(e) { 1 } else { 2 });
    }
    (faces, covered)
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

/// Is facet `f` not yet a union of faces? A facet with no points inside
/// whose subfaces are all faces is covered; otherwise look for an edge
/// meeting it (for a merged facet, see [`region_flood`]).
fn crossed_now(dt: &ExactDelaunay, f: &Facet) -> bool {
    let _s = crate::prof::scope("cdt faces: crossed");
    if !f.tris.is_empty() {
        return !region_flood(dt, f).1;
    }
    !(f.inside.is_empty() && f.subfaces.iter().all(|&s| dt.has_face(s)))
        && !dt.facet_crossings(f.corners, &f.vertices()).is_empty()
}

/// The edges face recovery must not remove: the subsegments, and the sides
/// of each facet's faces (its subfaces while it has no points inside, the
/// faces on it otherwise). Kept up to date facet by facet: a facet's faces
/// change only where the tetrahedralization made tets at one of its
/// vertices (see [`ExactDelaunay::log_touched`]), or where its own subfaces
/// or points change (marked dirty).
struct Protected {
    segments: BTreeSet<(u32, u32)>,
    /// The edges each facet gives.
    by_facet: Vec<Vec<(u32, u32)>>,
    /// How many facets give each edge.
    count: BTreeMap<(u32, u32), u32>,
    /// The facets at each vertex.
    at: BTreeMap<u32, Vec<usize>>,
    /// Facets to read again at the next refresh.
    dirty: BTreeSet<usize>,
    /// Whether each merged facet's faces cover it (see [`region_flood`]).
    covered: Vec<bool>,
}

impl Protected {
    /// Every facet's edges; `dt` logs from here on.
    fn new(dt: &mut ExactDelaunay, facets: &[Facet], segments: BTreeSet<(u32, u32)>) -> Self {
        dt.log_touched(true);
        let mut p = Protected {
            segments,
            by_facet: vec![Vec::new(); facets.len()],
            count: BTreeMap::new(),
            at: BTreeMap::new(),
            dirty: (0..facets.len()).collect(),
            covered: vec![false; facets.len()],
        };
        for (i, f) in facets.iter().enumerate() {
            for v in f.vertices() {
                p.at.entry(v).or_default().push(i);
            }
        }
        p.refresh(dt, facets);
        p
    }

    /// A facet's edges, and for a merged facet whether its faces cover it.
    fn facet_edges(dt: &ExactDelaunay, f: &Facet) -> (Vec<(u32, u32)>, bool) {
        let (faces, covered) = if !f.tris.is_empty() {
            region_flood(dt, f)
        } else if f.inside.is_empty() {
            (f.subfaces.clone(), false)
        } else {
            (faces_on(dt, &f.vertices()), false)
        };
        let mut edges: Vec<(u32, u32)> = faces
            .iter()
            .flat_map(|s| (0..3).map(move |k| key(s[k], s[(k + 1) % 3])))
            .collect();
        edges.sort_unstable();
        edges.dedup();
        (edges, covered)
    }

    /// Facet `i` gained vertex `v` (a Steiner point on it).
    fn add_vertex(&mut self, i: usize, v: u32) {
        self.at.entry(v).or_default().push(i);
        self.dirty.insert(i);
    }

    /// Read again the facets at the vertices `dt` logged and the dirty ones.
    fn refresh(&mut self, dt: &mut ExactDelaunay, facets: &[Facet]) {
        let _s = crate::prof::scope("cdt faces: protect");
        for v in dt.take_touched() {
            if let Some(fs) = self.at.get(&v) {
                self.dirty.extend(fs.iter().copied());
            }
        }
        for i in std::mem::take(&mut self.dirty) {
            let (new, covered) = Self::facet_edges(dt, &facets[i]);
            self.covered[i] = covered;
            for e in std::mem::replace(&mut self.by_facet[i], new) {
                if let Entry::Occupied(mut c) = self.count.entry(e) {
                    *c.get_mut() -= 1;
                    if *c.get() == 0 {
                        c.remove();
                    }
                }
            }
            for &e in &self.by_facet[i] {
                *self.count.entry(e).or_default() += 1;
            }
        }
        // Debug builds check the facets left alone against reading them all.
        #[cfg(debug_assertions)]
        for (i, f) in facets.iter().enumerate() {
            debug_assert_eq!(
                (self.by_facet[i].clone(), self.covered[i]),
                Self::facet_edges(dt, f),
                "facet {i}"
            );
        }
    }

    fn contains(&self, x: u32, y: u32) -> bool {
        let e = key(x, y);
        self.segments.contains(&e) || self.count.contains_key(&e)
    }
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
    recover_faces_with(points, &[], triangles, false, max_points, seed)
}

/// [`recover_faces`] with `extra` points too (see
/// [`super::segments::recover_segments_with`]): those exactly inside an
/// input triangle start as vertices inside its facet, so its subfaces have
/// them as corners. With `merge_flat`, coplanar triangles across flat edges
/// are one facet, and the flat edges are not segments; its subfaces carry
/// its first triangle. A triangle with points inside, or whose polygon is
/// longer than 512 vertices, is such a planar-region facet too.
///
/// # Errors
/// See [`FaceError`].
pub fn recover_faces_with(
    points: &[Point3],
    extra: &[Point],
    triangles: &[[u32; 3]],
    merge_flat: bool,
    max_points: usize,
    seed: u64,
) -> Result<RecoveredFaces, FaceError> {
    let setup = crate::prof::scope("cdt faces: flat edges");
    let flat = if merge_flat {
        let explicit: Vec<Point> = points.iter().map(|&p| Point::from(p)).collect();
        flat_edges(&explicit, triangles)
    } else {
        BTreeSet::new()
    };
    drop(setup);
    let mut segments = recover_segments_with(points, extra, triangles, &flat, max_points, seed)?;
    let setup = crate::prof::scope("cdt faces: setup");
    let mut inside_of = std::mem::take(&mut segments.inside);
    let limit = if max_points == 0 {
        20 * points.len() + extra.len() + 10_000
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
        let inside = inside_of
            .get_mut(ti as usize)
            .map(std::mem::take)
            .unwrap_or_default();
        // Flips work towards subfaces only where nothing is inside. With
        // `merge_flat`, only on a triangle with no flat side and a short
        // polygon (triangulating one is cubic in its length); the rest
        // become planar regions.
        let alone = !(0..3).any(|k| flat.contains(&key(t[k], t[(k + 1) % 3])));
        let subfaces =
            if inside.is_empty() && (!merge_flat || (alone && polygon.len() <= SUBFACE_POLYGON)) {
                triangulate(&dt, &polygon).ok_or(FaceError::Degenerate { triangle: ti })?
            } else {
                Vec::new()
            };
        facets.push(Facet {
            triangle: ti,
            corners: *t,
            polygon,
            subfaces,
            inside,
            tris: Vec::new(),
            boundary: Vec::new(),
        });
    }
    let mut facet_of: Vec<u32> = (0..triangles.len() as u32).collect();
    drop(setup);
    if merge_flat {
        let _s = crate::prof::scope("cdt faces: merge");
        facets = merge_facets(facets, triangles, &flat, &chains, &mut facet_of);
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
    let mut protected = Protected::new(&mut dt, &facets, segment_edges);
    // A facet with no points inside whose subfaces are all faces is covered;
    // otherwise look for an edge meeting it.
    // A merged facet's answer comes with its protected edges, kept up to
    // date.
    let crossed =
        |dt: &mut ExactDelaunay, protected: &mut Protected, facets: &[Facet], ti: usize| {
            if !facets[ti].tris.is_empty() {
                protected.refresh(dt, facets);
                return !protected.covered[ti];
            }
            crossed_now(dt, &facets[ti])
        };
    let crossed_before = (0..facets.len())
        .filter(|&ti| crossed(&mut dt, &mut protected, &facets, ti))
        .count();
    // 1. Flips.
    for pass in 0..MAX_PASSES {
        if pass > 0 {
            for (i, f) in facets.iter_mut().enumerate() {
                if f.tris.is_empty() && f.subfaces.iter().any(|&s| !dt.has_face(s)) {
                    if let Some(s) = triangulate(&dt, &f.polygon) {
                        f.subfaces = s;
                        protected.dirty.insert(i);
                    }
                }
            }
        }
        protected.refresh(&mut dt, &facets);
        let is_protected = |x: u32, y: u32| protected.contains(x, y);
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
        .filter(|&ti| crossed(&mut dt, &mut protected, &facets, ti))
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
            if !crossed(&mut dt, &mut protected, &facets, ti) {
                continue;
            }
            protected.refresh(&mut dt, &facets);
            let is_protected = |x: u32, y: u32| protected.contains(x, y);
            let f = &facets[ti];
            let _s = crate::prof::scope("cdt faces: cavity");
            let refill = if f.tris.is_empty() {
                dt.recover_facet_by_cavity(f.corners, &f.vertices(), &is_protected)
            } else {
                dt.recover_region_by_cavity(&f.tris, &f.vertices(), &is_protected)
            };
            if refill == Refill::GiftWrapped {
                wrapped += 1;
            }
            if refill.recovered() {
                // The cavity chose how the uncovered part is triangulated:
                // read it back, so the faces it made are protected.
                if facets[ti].inside.is_empty() && facets[ti].tris.is_empty() {
                    if let Some(s) = triangulate(&dt, &facets[ti].polygon) {
                        facets[ti].subfaces = s;
                        protected.dirty.insert(ti);
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
            let crossings = if facets[ti].tris.is_empty() {
                dt.facet_crossings(facets[ti].corners, &facets[ti].vertices())
            } else {
                dt.region_crossings(&facets[ti].tris, &facets[ti].vertices())
            };
            protected.refresh(&mut dt, &facets);
            let is_protected = |x: u32, y: u32| protected.contains(x, y);
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
                .flat_map(|&(x, y)| {
                    // On the triangle of a merged facet the edge crosses.
                    let f = &facets[ti];
                    let corners = f
                        .tris
                        .iter()
                        .find(|t| dt.crosses_region(x, y, std::slice::from_ref(*t)))
                        .map_or(f.corners, |t| t.0);
                    crossing_point(&dt, points, corners, x, y)
                })
                .take(TRIES)
                .collect();
            let mut added = 0;
            let mut gained = Vec::new();
            let most = PER_ROUND.min(1 << steiner_rounds[ti].min(8));
            steiner_rounds[ti] += 1;
            for p in candidates {
                if added == most || dt.points().len() >= limit {
                    break;
                }
                if let Ok(Inserted::New(m)) = dt.insert_protected(p, &is_protected) {
                    facets[ti].inside.push(m);
                    gained.push(m);
                    steiner_on_faces += 1;
                    inserted = true;
                    added += 1;
                }
            }
            for m in gained {
                protected.add_vertex(ti, m);
            }
        }
        if !inserted {
            break;
        }
    }
    dt.log_touched(false);
    let _s = crate::prof::scope("cdt faces: final");
    let mut missing = Vec::new();
    let mut subfaces = Vec::new();
    for f in &facets {
        if crossed_now(&dt, f) {
            missing.push(f.triangle);
        } else {
            let faces = if f.tris.is_empty() {
                faces_on(&dt, &f.vertices())
            } else {
                region_faces(&dt, f)
            };
            subfaces.extend(faces.into_iter().map(|v| Subface {
                triangle: f.triangle,
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
        facet_of,
    })
}

/// Merge the facets of coplanar triangles across `flat` edges: one facet
/// each, holding its triangles, every vertex on it and its boundary pieces;
/// a triangle with no flat side keeps its own facet. `facet_of` gets each
/// triangle's facet's first triangle.
fn merge_facets(
    facets: Vec<Facet>,
    triangles: &[[u32; 3]],
    flat: &BTreeSet<(u32, u32)>,
    chains: &std::collections::BTreeMap<(u32, u32), Vec<u32>>,
    facet_of: &mut [u32],
) -> Vec<Facet> {
    fn find(root: &mut [u32], x: u32) -> u32 {
        let mut r = x;
        while root[r as usize] != r {
            r = root[r as usize];
        }
        let mut y = x;
        while root[y as usize] != r {
            let next = root[y as usize];
            root[y as usize] = r;
            y = next;
        }
        r
    }
    let n = triangles.len();
    let mut root: Vec<u32> = (0..n as u32).collect();
    let mut by_edge: std::collections::BTreeMap<(u32, u32), Vec<u32>> =
        std::collections::BTreeMap::new();
    for (ti, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            by_edge
                .entry(key(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(ti as u32);
        }
    }
    for e in flat {
        if let Some(&[a, b]) = by_edge.get(e).map(Vec::as_slice) {
            let (ra, rb) = (find(&mut root, a), find(&mut root, b));
            if ra != rb {
                root[ra.max(rb) as usize] = ra.min(rb);
            }
        }
    }
    let mut groups: std::collections::BTreeMap<u32, Vec<u32>> = std::collections::BTreeMap::new();
    for ti in 0..n as u32 {
        let r = find(&mut root, ti);
        groups.entry(r).or_default().push(ti);
    }
    let mut old: Vec<Option<Facet>> = facets.into_iter().map(Some).collect();
    let mut out = Vec::with_capacity(groups.len());
    for (_, group) in groups {
        let first = group[0];
        for &ti in &group {
            facet_of[ti as usize] = first;
        }
        let sides = |t: [u32; 3]| [0, 1, 2].map(|k| flat.contains(&key(t[k], t[(k + 1) % 3])));
        // A lone triangle with subfaces keeps its facet.
        if group.len() == 1
            && old[first as usize]
                .as_ref()
                .is_some_and(|f| !f.subfaces.is_empty())
        {
            if let Some(f) = old[first as usize].take() {
                out.push(f);
            }
            continue;
        }
        let mut vertices: BTreeSet<u32> = BTreeSet::new();
        let mut tris = Vec::new();
        let mut boundary = Vec::new();
        for &ti in &group {
            let t = triangles[ti as usize];
            let flats = sides(t);
            tris.push((t, flats));
            vertices.extend(t);
            if let Some(f) = old[ti as usize].take() {
                vertices.extend(f.inside);
            }
            for k in 0..3 {
                if flats[k] {
                    continue;
                }
                let (u, v) = (t[k], t[(k + 1) % 3]);
                let mut chain = chains
                    .get(&key(u, v))
                    .cloned()
                    .unwrap_or_else(|| vec![u, v]);
                if chain.first() != Some(&u) {
                    chain.reverse();
                }
                vertices.extend(&chain);
                boundary.extend(chain.windows(2).map(|w| (w[0], w[1])));
            }
        }
        out.push(Facet {
            triangle: first,
            corners: triangles[first as usize],
            polygon: vertices.into_iter().map(|v| (v, [false; 3])).collect(),
            inside: Vec::new(),
            subfaces: Vec::new(),
            tris,
            boundary,
        });
    }
    out
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
