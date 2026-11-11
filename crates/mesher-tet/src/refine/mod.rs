//! Conforming Delaunay refinement of a chamfered surface (Milestone 14,
//! after Diazzi, Dai, Panozzo and Attene, "Surface chamfering for robust
//! tetrahedral meshing", SIGGRAPH 2026, from the paper's text; the
//! refinement itself is Ruppert's and Shewchuk's).
//!
//! The chamfered surface (see [`crate::chamfer`]) has no acute angle, so
//! Delaunay refinement of it is sure to finish. Its vertices and the corners
//! of an outer box are tetrahedralized, and the tetrahedralization stays
//! Delaunay: constraints are recovered by splitting them, never by flips.
//!
//! - **Segments**: every face edge of the chamfered surface that is a
//!   segment, and the box's edges, each kept as a chain of pieces. A piece
//!   is split at its middle while it is not an edge or is *encroached* (a
//!   vertex of the tets around it strictly inside its diametral sphere:
//!   exactly, the angle it subtends there is obtuse). Split points stay
//!   exact: an LNC point of the input or box edge a segment lies on; a BPT
//!   point of the input triangle, at the average of the ends' barycentric
//!   weights, on a bridge segment (exactly on the triangle, not necessarily
//!   on the segment, which the paper accepts: the chamfered surface is only
//!   an intermediate, and the enriched CDT conforms to the input exactly).
//! - **Faces**: each planar region (chamfered faces merged across flat
//!   edges, and the box's sides) keeps a 2D Delaunay triangulation of its
//!   vertices ([`plane::PlaneDelaunay`]). A triangle inside the region that
//!   is not a face of the tetrahedralization, or whose equatorial sphere
//!   holds a vertex of the two tets at it, is split at its circumcentre (a
//!   BPT point of the region's plane), unless that point encroaches a piece
//!   of the region's boundary, which is split instead. The equatorial-sphere
//!   test is f64: it only chooses where to put points, while every
//!   topological decision stays exact.

use std::collections::{BTreeMap, BTreeSet};

use mesher_geom::implicit::{self, Kind, Point};
use mesher_geom::{Point3, Sign};

use crate::cdt::{DelaunayError, ExactDelaunay, GHOST, Inserted};
use crate::chamfer::{ChamferVertex, FilteredChamfer};

mod plane;

use plane::PlaneDelaunay;

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

/// Where a vertex sits on a carrier.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Param {
    S(f64),
    W([f64; 3]),
}

/// The point halfway between two parameters on a carrier, exactly on it.
fn middle(carrier: Carrier, p0: Param, p1: Param) -> Option<(Point, Param)> {
    match (carrier, p0, p1) {
        (Carrier::Line { a, b }, Param::S(s0), Param::S(s1)) => {
            let s = 0.5 * (s0 + s1);
            // a + s(b − a) = s·b + (1 − s)·a
            Point::lnc(b, a, s).ok().map(|p| (p, Param::S(s)))
        }
        (Carrier::Triangle { tri }, Param::W(w0), Param::W(w1)) => {
            let w = [0, 1, 2].map(|k| 0.5 * (w0[k] + w1[k]));
            Point::bpt(tri, w).ok().map(|p| (p, Param::W(w)))
        }
        _ => None,
    }
}

/// The parameter of point `p` on line `a b` (`p = a + s(b − a)`), from its
/// construction: an end, or an LNC point of the line's ends in either
/// order.
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

/// Is `p` strictly inside the diametral sphere of `a b` (is the angle
/// `a p b` obtuse)? Exact.
fn encroaches(a: &Point, p: &Point, b: &Point) -> bool {
    implicit::dot_sign(a, p, b) == Sign::Negative
}

fn sub(a: Point3, b: Point3) -> Point3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: Point3, b: Point3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Point3, b: Point3) -> Point3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The circumcentre of triangle `a b c` in its plane (f64).
fn circumcentre(a: Point3, b: Point3, c: Point3) -> Option<Point3> {
    let (u, v) = (sub(b, a), sub(c, a));
    let n = cross(u, v);
    let nn = dot(n, n);
    if !(nn > 0.0 && nn.is_finite()) {
        return None;
    }
    let (uu, vv) = (dot(u, u), dot(v, v));
    let x = cross(n, u);
    let y = cross(v, n);
    let o = [0, 1, 2].map(|k| (uu * y[k] + vv * x[k]) / (2.0 * nn));
    Some([a[0] + o[0], a[1] + o[1], a[2] + o[2]])
}

/// Barycentric weights of `p` (in the plane) in triangle `t`, unclamped:
/// the point may lie anywhere on the plane.
fn plane_weights(p: Point3, t: [Point3; 3]) -> [f64; 3] {
    let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
    let area = |a: Point3, b: Point3, c: Point3| dot(n, cross(sub(b, a), sub(c, a)));
    let whole = area(t[0], t[1], t[2]);
    [
        area(p, t[1], t[2]) / whole,
        area(t[0], p, t[2]) / whole,
        area(t[0], t[1], p) / whole,
    ]
}

/// A segment as a chain of pieces: its carrier, its vertices in order with
/// their parameters, and the regions it bounds (with whether the region's
/// faces run along it backwards).
#[derive(Debug, Clone)]
struct Chain {
    carrier: Carrier,
    verts: Vec<u32>,
    params: Vec<Param>,
    regions: Vec<(usize, bool)>,
}

/// A planar region: a plane (three explicit points on it), the projection
/// axis, whether that projection keeps the faces' orientation, its 2D
/// Delaunay triangulation, the chains bounding it, and whether every
/// triangle of it is inside (a box side).
#[derive(Debug, Clone)]
struct Region {
    plane: [Point3; 3],
    positive: bool,
    pd: PlaneDelaunay,
    chains: Vec<usize>,
    all_inside: bool,
    dirty: bool,
}

/// The refined tetrahedralization: the chamfered surface's vertices and the
/// outer box, every segment a chain of Delaunay edges and (after face
/// refinement) every planar region a union of Delaunay faces.
#[derive(Debug, Clone)]
pub struct Refined {
    pub dt: ExactDelaunay,
    /// Every subsegment, by its end vertices.
    pub subsegments: Vec<[u32; 2]>,
    /// Every subface (after face refinement), with its region.
    pub subfaces: Vec<([u32; 3], u32)>,
    /// The tetrahedralization's vertex of each chamfered-surface vertex in
    /// use (`None` for one no face uses).
    pub vertex_of: Vec<Option<u32>>,
    /// The box corners' vertices.
    pub box_corners: [u32; 8],
    /// Points inserted on segments, and on faces.
    pub steiner_segments: usize,
    pub steiner_faces: usize,
    /// Subfaces found encroached or missing whose fix added no point (a
    /// circumcentre already a vertex), left as they are.
    pub stuck_subfaces: usize,
}

/// Refinement under way.
struct Refiner {
    dt: ExactDelaunay,
    chains: Vec<Chain>,
    regions: Vec<Region>,
    /// The chains through each vertex.
    at_vertex: Vec<Vec<u32>>,
    /// The regions each vertex lies on.
    region_of: Vec<Vec<u32>>,
    queue: Vec<(u32, u32, u32)>,
    /// Subfaces (sorted vertices) whose fix added no point: set aside, not
    /// retried.
    stuck: BTreeSet<[u32; 3]>,
    max_points: usize,
    steiner_segments: usize,
    steiner_faces: usize,
}

impl Refiner {
    fn pt(&self, v: u32) -> &Point {
        &self.dt.points()[v as usize]
    }

    /// Insert `p`; check again the pieces at its neighbours it lies inside
    /// or on the diametral sphere of, and mark their regions for checking.
    /// `None` if it equals a vertex.
    fn insert(&mut self, p: Point) -> Result<Option<u32>, RefineError> {
        if self.dt.points().len() >= self.max_points {
            return Err(RefineError::Budget {
                points: self.dt.points().len(),
                limit: self.max_points,
            });
        }
        let m = match self.dt.insert(p)? {
            Inserted::New(m) => m,
            Inserted::Duplicate(_) => return Ok(None),
        };
        let n = self.dt.points().len();
        self.at_vertex.resize(n, Vec::new());
        self.region_of.resize(n, Vec::new());
        for w in self.dt.neighbours(m) {
            for &r in &self.region_of[w as usize] {
                self.regions[r as usize].dirty = true;
            }
            for &c in &self.at_vertex[w as usize] {
                let chain = &self.chains[c as usize];
                let Some(i) = chain.verts.iter().position(|&x| x == w) else {
                    continue;
                };
                for (a, b) in [
                    (i.checked_sub(1).map(|j| chain.verts[j]), Some(w)),
                    (Some(w), chain.verts.get(i + 1).copied()),
                ] {
                    if let (Some(a), Some(b)) = (a, b) {
                        if implicit::dot_sign(self.pt(a), &p, self.pt(b)) != Sign::Positive {
                            self.queue.push((c, a, b));
                        }
                    }
                }
            }
        }
        Ok(Some(m))
    }

    /// Is piece `a b` a Delaunay edge no vertex of the tets around it
    /// encroaches?
    fn piece_ok(&self, a: u32, b: u32) -> bool {
        let Some(ring) = self.dt.edge_ring(a, b) else {
            return false;
        };
        let (pa, pb) = (self.pt(a), self.pt(b));
        !ring.iter().any(|&t| {
            self.dt.tet(t).is_some_and(|v| {
                v.iter()
                    .any(|&w| w != a && w != b && encroaches(pa, self.pt(w), pb))
            })
        })
    }

    /// Split piece `a b` of chain `c` at its middle.
    fn split(&mut self, c: u32, a: u32, b: u32) -> Result<(), RefineError> {
        let chain = &self.chains[c as usize];
        let Some(i) = chain.verts.iter().position(|&x| x == a) else {
            return Ok(());
        };
        if chain.verts.get(i + 1) != Some(&b) {
            return Ok(());
        }
        let (p, param) = middle(chain.carrier, chain.params[i], chain.params[i + 1])
            .ok_or(RefineError::Construction)?;
        let Some(m) = self.insert(p)? else {
            return Err(RefineError::Duplicate { point: p.approx() });
        };
        self.steiner_segments += 1;
        let chain = &mut self.chains[c as usize];
        chain.verts.insert(i + 1, m);
        chain.params.insert(i + 1, param);
        let regions = chain.regions.clone();
        self.at_vertex[m as usize].push(c);
        for (r, _) in regions {
            self.region_of[m as usize].push(r as u32);
            let region = &mut self.regions[r];
            region.pd.insert(m, p);
            region.dirty = true;
        }
        self.queue.push((c, a, m));
        self.queue.push((c, m, b));
        Ok(())
    }

    /// Split the queued pieces that are missing or encroached, and what
    /// that queues in turn.
    fn drain(&mut self) -> Result<(), RefineError> {
        while let Some((c, a, b)) = self.queue.pop() {
            let chain = &self.chains[c as usize];
            let adjacent = chain.verts.windows(2).any(|w| w[0] == a && w[1] == b);
            if !adjacent {
                continue;
            }
            let _s = crate::prof::scope("refine segments: check");
            if self.piece_ok(a, b) {
                continue;
            }
            drop(_s);
            self.split(c, a, b)?;
        }
        Ok(())
    }

    /// Split pieces until every one is an unencroached Delaunay edge,
    /// finishing with a full check.
    fn segments(&mut self) -> Result<(), RefineError> {
        loop {
            self.drain()?;
            let _s = crate::prof::scope("refine segments: final pass");
            for (c, chain) in self.chains.iter().enumerate() {
                for w in chain.verts.windows(2) {
                    if !self.piece_ok(w[0], w[1]) {
                        self.queue.push((c as u32, w[0], w[1]));
                    }
                }
            }
            if self.queue.is_empty() {
                return Ok(());
            }
        }
    }

    /// The triangles of region `r` inside it: all of them for a box side;
    /// otherwise a flood from the inner side of each boundary piece that
    /// never crosses one.
    fn inside(&self, r: usize) -> Vec<[u32; 3]> {
        let region = &self.regions[r];
        let tris = region.pd.triangles();
        if region.all_inside {
            return tris.into_iter().map(|(v, _)| v).collect();
        }
        let mut by_edge: BTreeMap<(u32, u32), usize> = BTreeMap::new();
        for (i, (v, _)) in tris.iter().enumerate() {
            for k in 0..3 {
                by_edge.insert((v[(k + 1) % 3], v[(k + 2) % 3]), i);
            }
        }
        let mut walls: BTreeSet<(u32, u32)> = BTreeSet::new();
        let mut inside = vec![false; tris.len()];
        let mut stack = Vec::new();
        for &c in &region.chains {
            let chain = &self.chains[c];
            let backwards = chain
                .regions
                .iter()
                .find(|(x, _)| *x == r)
                .is_some_and(|(_, b)| *b);
            for w in chain.verts.windows(2) {
                let (a, b) = if backwards {
                    (w[1], w[0])
                } else {
                    (w[0], w[1])
                };
                walls.insert((a.min(b), a.max(b)));
                // The faces run a → b with the inside on their left; the
                // triangulation is counterclockwise in the projection.
                let key = if region.positive { (a, b) } else { (b, a) };
                if let Some(&t) = by_edge.get(&key) {
                    if !inside[t] {
                        inside[t] = true;
                        stack.push(t);
                    }
                }
            }
        }
        while let Some(t) = stack.pop() {
            let (v, n) = tris[t];
            for k in 0..3 {
                let (a, b) = (v[(k + 1) % 3], v[(k + 2) % 3]);
                if walls.contains(&(a.min(b), a.max(b))) {
                    continue;
                }
                if let Some(u) = n[k] {
                    if !inside[u] {
                        inside[u] = true;
                        stack.push(u);
                    }
                }
            }
        }
        tris.iter()
            .zip(&inside)
            .filter(|(_, i)| **i)
            .map(|(t, _)| t.0)
            .collect()
    }

    /// Is triangle `v` a face whose equatorial sphere holds neither apex of
    /// the tets at it? (The sphere test is f64, with an apex on the sphere to
    /// rounding not counted: it only chooses where to refine.)
    fn subface_ok(&self, v: [u32; 3]) -> bool {
        let Some(apexes) = self.dt.face_apexes(v) else {
            return false;
        };
        let [a, b, c] = v.map(|x| self.pt(x).approx());
        let Some(o) = circumcentre(a, b, c) else {
            return true;
        };
        let r2 = dot(sub(a, o), sub(a, o));
        !apexes.iter().any(|&w| {
            w != GHOST && {
                let p = self.pt(w).approx();
                dot(sub(p, o), sub(p, o)) < r2 * (1.0 - 1e-9)
            }
        })
    }

    /// One pass over region `r`: each triangle inside it that is missing or
    /// encroached is split at its circumcentre, or the boundary pieces that
    /// point encroaches are split instead.
    fn region(&mut self, r: usize) -> Result<(), RefineError> {
        self.regions[r].dirty = false;
        let bad: Vec<[u32; 3]> = self
            .inside(r)
            .into_iter()
            .filter(|&v| !self.is_stuck(v) && !self.subface_ok(v))
            .collect();
        for v in bad {
            let before = self.dt.points().len();
            let [a, b, c] = v.map(|x| self.pt(x).approx());
            let Some(o) = circumcentre(a, b, c) else {
                continue;
            };
            let plane = self.regions[r].plane;
            let Ok(p) = Point::bpt(plane, plane_weights(o, plane)) else {
                continue;
            };
            // Boundary pieces it encroaches are split instead.
            let mut pieces = Vec::new();
            for &ch in &self.regions[r].chains {
                for w in self.chains[ch].verts.windows(2) {
                    if encroaches(self.pt(w[0]), &p, self.pt(w[1])) {
                        pieces.push((ch, w[0], w[1]));
                    }
                }
            }
            if pieces.is_empty() {
                if let Some(m) = self.insert(p)? {
                    self.steiner_faces += 1;
                    self.region_of[m as usize].push(r as u32);
                    let region = &mut self.regions[r];
                    region.pd.insert(m, p);
                    region.dirty = true;
                }
            } else {
                for (ch, x, y) in pieces {
                    self.split(ch as u32, x, y)?;
                }
            }
            self.drain()?;
            // No point added: set it aside if it is still there (a triangle
            // an earlier fix in this pass already removed is not stuck).
            if self.dt.points().len() == before && self.regions[r].pd.has_triangle(v) {
                let mut key = v;
                key.sort_unstable();
                self.stuck.insert(key);
            }
        }
        Ok(())
    }

    fn is_stuck(&self, v: [u32; 3]) -> bool {
        let mut key = v;
        key.sort_unstable();
        self.stuck.contains(&key)
    }

    /// Faces: regions checked while any is marked, then a full pass.
    fn faces(&mut self) -> Result<(), RefineError> {
        self.segments()?;
        loop {
            self.drain()?;
            if let Some(r) = self.regions.iter().position(|x| x.dirty) {
                let _s = crate::prof::scope("refine faces: region");
                self.region(r)?;
                continue;
            }
            self.segments()?;
            let _s = crate::prof::scope("refine faces: final pass");
            let mut any = false;
            for r in 0..self.regions.len() {
                if self
                    .inside(r)
                    .into_iter()
                    .any(|v| !self.is_stuck(v) && !self.subface_ok(v))
                {
                    self.regions[r].dirty = true;
                    any = true;
                }
            }
            if !any {
                return Ok(());
            }
        }
    }

    fn finish(self, vertex_of: Vec<Option<u32>>, box_corners: [u32; 8]) -> Refined {
        let mut subsegments = Vec::new();
        for c in &self.chains {
            for w in c.verts.windows(2) {
                subsegments.push([w[0], w[1]]);
            }
        }
        let mut subfaces = Vec::new();
        for r in 0..self.regions.len() {
            for v in self.inside(r) {
                subfaces.push((v, r as u32));
            }
        }
        Refined {
            dt: self.dt,
            subsegments,
            subfaces,
            vertex_of,
            box_corners,
            steiner_segments: self.steiner_segments,
            steiner_faces: self.steiner_faces,
            stuck_subfaces: self.stuck.len(),
        }
    }
}

/// The refiner, each chamfered-surface vertex's vertex, and the box corners.
type SetUp = (Refiner, Vec<Option<u32>>, [u32; 8]);

/// The projection axis for the plane of `tri`: the one along which its
/// projected orientation is not zero (the largest normal component first).
fn axis_of(tri: [Point3; 3]) -> Option<(usize, bool)> {
    let n = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]));
    let mut axes = [0usize, 1, 2];
    axes.sort_by(|&i, &j| n[j].abs().total_cmp(&n[i].abs()));
    let p = tri.map(|q| Point::explicit(q).ok());
    let [Some(a), Some(b), Some(c)] = p else {
        return None;
    };
    axes.into_iter()
        .find_map(|k| match implicit::orient2d_projected(&a, &b, &c, k) {
            Sign::Zero => None,
            s => Some((k, s == Sign::Positive)),
        })
}

/// Set up refinement of a filtered chamfered surface: its vertices in use
/// and a box's corners (its bounding box grown by `margin` of its diagonal)
/// tetrahedralized, its segments and the box's edges as chains, its planar
/// regions and the box's sides with their 2D triangulations.
fn set_up(
    surface: &FilteredChamfer,
    points: &[Point3],
    triangles: &[[u32; 3]],
    margin: f64,
    max_points: usize,
) -> Result<SetUp, RefineError> {
    let c = &surface.chamfered;
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
    let dt = ExactDelaunay::new(pts.clone(), 5)?;
    if let Some(&(v, _)) = dt.duplicates().first() {
        return Err(RefineError::Duplicate {
            point: dt.points()[v as usize].approx(),
        });
    }
    let n = dt.points().len();
    // Planar regions: chamfered faces merged across flat edges.
    let faces = &c.faces;
    let mut root: Vec<usize> = (0..faces.len()).collect();
    fn find(root: &mut [usize], x: usize) -> usize {
        let mut y = x;
        while root[y] != y {
            root[y] = root[root[y]];
            y = root[y];
        }
        y
    }
    for e in surface.acute.edges.values() {
        if e.flat {
            let (a, b) = (e.triangles[0] as usize, e.triangles[1] as usize);
            if a < faces.len() && b < faces.len() {
                let (ra, rb) = (find(&mut root, a), find(&mut root, b));
                if ra != rb {
                    root[ra.max(rb)] = ra.min(rb);
                }
            }
        }
    }
    let mut region_index: BTreeMap<usize, usize> = BTreeMap::new();
    let mut regions: Vec<Region> = Vec::new();
    let mut region_of: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut face_region = vec![0usize; faces.len()];
    for (fi, f) in faces.iter().enumerate() {
        let rt = find(&mut root, fi);
        let r = *region_index.entry(rt).or_insert_with(|| {
            let t = triangles[surface.kept_ids[f.triangle as usize] as usize];
            let plane = t.map(|v| points[v as usize]);
            let (axis, positive) = axis_of(plane).unwrap_or((2, true));
            regions.push(Region {
                plane,
                positive,
                pd: PlaneDelaunay::new(axis),
                chains: Vec::new(),
                all_inside: false,
                dirty: true,
            });
            regions.len() - 1
        });
        face_region[fi] = r;
        for &v in &f.polygon {
            if let Some(d) = vertex_of[v as usize] {
                if !region_of[d as usize].contains(&(r as u32)) {
                    region_of[d as usize].push(r as u32);
                    regions[r].pd.insert(d, dt.points()[d as usize]);
                }
            }
        }
    }
    // Chains: the segments, then the box's edges.
    let mut chains: Vec<Chain> = Vec::new();
    let mut chain_of: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for (fi, f) in faces.iter().enumerate() {
        let k = f.polygon.len();
        let tri_index = surface.kept_ids[f.triangle as usize];
        let corners = triangles[tri_index as usize].map(|v| points[v as usize]);
        for i in 0..k {
            let (u, v) = (f.polygon[i], f.polygon[(i + 1) % k]);
            if !c.is_segment(u, v, &surface.acute) {
                continue;
            }
            let (Some(du), Some(dv)) = (vertex_of[u as usize], vertex_of[v as usize]) else {
                continue;
            };
            let r = face_region[fi];
            if let Some(&id) = chain_of.get(&(du.min(dv), du.max(dv))) {
                let backwards = chains[id].verts[0] != du;
                if !chains[id].regions.iter().any(|(x, _)| *x == r) {
                    chains[id].regions.push((r, backwards));
                    regions[r].chains.push(id);
                }
                continue;
            }
            let (pu, pv) = (&c.points[u as usize], &c.points[v as usize]);
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
            let (carrier, params) = match edge {
                Some((e0, e1)) => {
                    let (a, b) = (points[e0 as usize], points[e1 as usize]);
                    let (Some(s0), Some(s1)) = (line_param(pu, a, b), line_param(pv, a, b)) else {
                        return Err(RefineError::Construction);
                    };
                    (Carrier::Line { a, b }, vec![Param::S(s0), Param::S(s1)])
                }
                None => {
                    let (Some(w0), Some(w1)) = (tri_param(pu, corners), tri_param(pv, corners))
                    else {
                        return Err(RefineError::Construction);
                    };
                    (
                        Carrier::Triangle { tri: corners },
                        vec![Param::W(w0), Param::W(w1)],
                    )
                }
            };
            let id = chains.len();
            chain_of.insert((du.min(dv), du.max(dv)), id);
            chains.push(Chain {
                carrier,
                verts: vec![du, dv],
                params,
                regions: vec![(r, false)],
            });
            regions[r].chains.push(id);
        }
    }
    // The box: its edges as chains, its sides as regions.
    let sides: [[usize; 4]; 6] = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let mut box_chain: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for side in sides {
        let plane = [corner(side[0]), corner(side[1]), corner(side[2])];
        let (axis, positive) = axis_of(plane).unwrap_or((2, true));
        let r = regions.len();
        let mut pd = PlaneDelaunay::new(axis);
        for &i in &side {
            let d = box_corners[i];
            pd.insert(d, dt.points()[d as usize]);
            region_of[d as usize].push(r as u32);
        }
        let mut region_chains = Vec::new();
        for k in 0..4 {
            let (i, j) = (side[k], side[(k + 1) % 4]);
            let key = (i.min(j), i.max(j));
            let id = *box_chain.entry(key).or_insert_with(|| {
                chains.push(Chain {
                    carrier: Carrier::Line {
                        a: corner(key.0),
                        b: corner(key.1),
                    },
                    verts: vec![box_corners[key.0], box_corners[key.1]],
                    params: vec![Param::S(0.0), Param::S(1.0)],
                    regions: Vec::new(),
                });
                chains.len() - 1
            });
            chains[id].regions.push((r, false));
            region_chains.push(id);
        }
        regions.push(Region {
            plane,
            positive,
            pd,
            chains: region_chains,
            all_inside: true,
            dirty: true,
        });
    }
    let mut at_vertex: Vec<Vec<u32>> = vec![Vec::new(); n];
    for (id, ch) in chains.iter().enumerate() {
        for &v in &ch.verts {
            at_vertex[v as usize].push(id as u32);
        }
    }
    let queue: Vec<(u32, u32, u32)> = chains
        .iter()
        .enumerate()
        .rev()
        .map(|(id, ch)| (id as u32, ch.verts[0], ch.verts[1]))
        .collect();
    let refiner = Refiner {
        dt,
        chains,
        regions,
        at_vertex,
        region_of,
        queue,
        stuck: BTreeSet::new(),
        max_points,
        steiner_segments: 0,
        steiner_faces: 0,
    };
    Ok((refiner, vertex_of, box_corners))
}

/// Tetrahedralize a filtered chamfered surface (see
/// [`crate::chamfer::chamfer_filtered`]) and the corners of a box around it
/// (its bounding box grown by `margin` of its diagonal on every side), and
/// recover every segment and box edge as Delaunay edges. `points` and
/// `triangles` are the input the surface was chamfered from; `max_points`
/// bounds the vertex count.
///
/// # Errors
/// See [`RefineError`].
pub fn recover_segments_conforming(
    surface: &FilteredChamfer,
    points: &[Point3],
    triangles: &[[u32; 3]],
    margin: f64,
    max_points: usize,
) -> Result<Refined, RefineError> {
    let (mut r, vertex_of, box_corners) = set_up(surface, points, triangles, margin, max_points)?;
    r.segments()?;
    for region in &mut r.regions {
        region.dirty = false;
    }
    let mut out = r.finish(vertex_of, box_corners);
    out.subfaces.clear();
    Ok(out)
}

/// [`recover_segments_conforming`], then every planar region of the
/// chamfered surface and every side of the box recovered as a union of
/// Delaunay faces, none encroached.
///
/// # Errors
/// See [`RefineError`].
pub fn recover_faces_conforming(
    surface: &FilteredChamfer,
    points: &[Point3],
    triangles: &[[u32; 3]],
    margin: f64,
    max_points: usize,
) -> Result<Refined, RefineError> {
    let (mut r, vertex_of, box_corners) = set_up(surface, points, triangles, margin, max_points)?;
    r.faces()?;
    Ok(r.finish(vertex_of, box_corners))
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

    fn tetrahedron() -> (Vec<Point3>, Vec<[u32; 3]>) {
        let p = vec![
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        let f = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        (p, f)
    }

    fn area(t: [Point3; 3]) -> f64 {
        let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
        0.5 * dot(n, n).sqrt()
    }

    /// Every subsegment an unencroached Delaunay edge, every subface a face,
    /// the tetrahedralization valid and Delaunay; with faces, each region's
    /// subfaces cover its chamfered faces' area, and the box sides theirs.
    fn check(name: &str, p: &[Point3], f: &[[u32; 3]], faces: bool) -> Refined {
        let surface = chamfer_filtered(p, f, 1e-6).unwrap();
        let r = if faces {
            recover_faces_conforming(&surface, p, f, 0.2, 200_000).unwrap()
        } else {
            recover_segments_conforming(&surface, p, f, 0.2, 200_000).unwrap()
        };
        r.dt.check().unwrap();
        r.dt.check_delaunay().unwrap();
        let at = |v: u32| r.dt.points()[v as usize].approx();
        for s in &r.subsegments {
            assert!(r.dt.has_edge(s[0], s[1]), "{name}: {s:?}");
        }
        if faces {
            assert_eq!(r.stuck_subfaces, 0, "{name}");
            for (v, _) in &r.subfaces {
                assert!(r.dt.has_face(*v), "{name}: {v:?}");
            }
            // The chamfered faces' area, against the subfaces off the box.
            let c = &surface.chamfered;
            let chamfered: f64 = c
                .faces
                .iter()
                .map(|f| {
                    let q: Vec<Point3> = f
                        .polygon
                        .iter()
                        .map(|&v| c.points[v as usize].approx())
                        .collect();
                    (1..q.len() - 1)
                        .map(|i| area([q[0], q[i], q[i + 1]]))
                        .sum::<f64>()
                })
                .sum();
            let on_box = |v: &[u32; 3]| {
                let b = r.box_corners.map(at);
                let (lo, hi) = (b[0], b[7]);
                (0..3).any(|k| {
                    v.iter().all(|&x| at(x)[k] == lo[k]) || v.iter().all(|&x| at(x)[k] == hi[k])
                })
            };
            let surface_area: f64 = r
                .subfaces
                .iter()
                .filter(|(v, _)| !on_box(v))
                .map(|(v, _)| area(v.map(at)))
                .sum();
            assert!(
                (surface_area - chamfered).abs() <= 1e-9 * chamfered,
                "{name}: {surface_area} against {chamfered}"
            );
        }
        r
    }

    #[test]
    fn segments_of_a_cube_and_a_chamfered_tetrahedron_are_recovered() {
        let (p, f) = cube();
        let r = check("cube", &p, &f, false);
        assert!(r.subsegments.len() >= 24);
        let (p, f) = tetrahedron();
        let r = check("tetrahedron", &p, &f, false);
        assert!(r.subsegments.len() >= 24 + 12);
    }

    #[test]
    fn faces_of_a_cube_and_a_chamfered_tetrahedron_are_recovered() {
        let (p, f) = cube();
        check("cube", &p, &f, true);
        let (p, f) = tetrahedron();
        check("tetrahedron", &p, &f, true);
    }
}
