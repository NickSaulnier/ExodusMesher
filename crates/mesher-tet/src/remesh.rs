//! Isotropic surface remeshing to a size field (after Botsch and Kobbelt,
//! 2004): repeated passes of splitting long edges, collapsing short ones,
//! Delaunay edge flips (instead of their valence flips; see `try_flip`), and
//! tangential smoothing, with every
//! vertex projected back onto the input surface.
//!
//! Features are kept: edges where the surface turns by more than the feature
//! angle, or where the group (OBJ group, hence side set) changes, are
//! feature lines. Vertices on them only move along them (projected onto the
//! input's feature segments), and corners (feature vertices that are not on
//! a single smooth line) never move. Smooth vertices are projected onto the
//! nearest input triangle facing their way.
//!
//! The result is a closed, consistently oriented, manifold surface whose
//! triangles each record the input triangle they lie on (for groups). It is
//! an approximation of the input: the caller measures its deviation.

use std::collections::{HashMap, HashSet, VecDeque};

use mesher_geom::bvh::{Bvh, closest_point_on_segment, closest_point_on_triangle};
use mesher_geom::vec3::{add, cross, distance_squared, dot, norm, scale, sub};
use mesher_geom::{Point3, TriangleSoup};

/// How the remesher works.
#[derive(Debug, Clone, PartialEq)]
pub struct RemeshParams {
    /// Edges turning by more than this are features, degrees.
    pub feature_degrees: f64,
    /// Passes of split, collapse, flip and smooth.
    pub iterations: usize,
    /// Refuse to grow beyond this many triangles.
    pub max_triangles: usize,
    /// A feature edge between two corners shorter than this fraction of the
    /// target size is a feature too small to mesh: one corner merges into
    /// the other (moving the surface by at most the edge's length).
    pub small_feature: f64,
}

impl Default for RemeshParams {
    fn default() -> RemeshParams {
        RemeshParams {
            feature_degrees: 40.0,
            iterations: 6,
            max_triangles: 2_000_000,
            small_feature: 0.1,
        }
    }
}

/// A remeshed surface.
#[derive(Debug, Clone, PartialEq)]
pub struct Remeshed {
    pub points: Vec<Point3>,
    pub triangles: Vec<[u32; 3]>,
    /// The input triangle each triangle lies on (for its group).
    pub origin: Vec<u32>,
    /// Features too small to mesh that were merged away
    /// (`RemeshParams::small_feature`): corners merged into corners, and
    /// feature vertices merged across a strip into another line.
    pub merged_small_features: usize,
}

/// Why a surface could not be remeshed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RemeshError {
    /// The input is not a closed, consistently oriented manifold surface.
    Surface(String),
    /// The result would exceed `max_triangles`.
    TooLarge { triangles: usize, limit: usize },
    /// The passes broke the surface (a bug; reported, not hidden).
    Broken(String),
}

impl std::fmt::Display for RemeshError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RemeshError::Surface(e) => write!(f, "surface cannot be remeshed: {e}"),
            RemeshError::TooLarge { triangles, limit } => {
                write!(
                    f,
                    "remeshing would make {triangles} triangles (limit {limit})"
                )
            }
            RemeshError::Broken(e) => write!(f, "remeshing broke the surface: {e}"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Smooth,
    Feature,
    Corner,
}

fn key(u: u32, v: u32) -> (u32, u32) {
    (u.min(v), u.max(v))
}

fn normal(a: Point3, b: Point3, c: Point3) -> Point3 {
    cross(sub(b, a), sub(c, a))
}

fn unit(v: Point3) -> Point3 {
    let l = norm(v);
    if l > 0.0 { scale(v, 1.0 / l) } else { [0.0; 3] }
}

/// The input surface, for projection.
struct Input<'a> {
    points: &'a [Point3],
    triangles: &'a [[u32; 3]],
    normals: Vec<Point3>,
    bvh: Bvh<'a>,
    /// Triangles around each input vertex.
    vertex_tris: Vec<Vec<u32>>,
    /// Feature segments (input vertex pairs) and those at each input vertex.
    segments: Vec<(u32, u32)>,
    vertex_segments: Vec<Vec<u32>>,
    /// Whether an input vertex is on one smooth feature line (two segments,
    /// not turning sharply), so a feature vertex may slide past it.
    line: Vec<bool>,
}

impl Input<'_> {
    /// The closest point on the input surface facing along `n`, starting
    /// from triangle `home`: the global nearest point if it faces the same
    /// way, else the nearest among the triangles around `home`.
    fn project(&self, p: Point3, n: Point3, home: u32) -> (Point3, u32) {
        if let Some(near) = self.bvh.nearest(p) {
            if dot(self.normals[near.triangle as usize], n) > 0.0 {
                return (near.point, near.triangle);
            }
        }
        let mut best = (closest(self, p, home), home);
        for &v in &self.triangles[home as usize] {
            for &t in &self.vertex_tris[v as usize] {
                if dot(self.normals[t as usize], n) <= 0.0 {
                    continue;
                }
                let q = closest(self, p, t);
                if distance_squared(p, q) < distance_squared(p, best.0) {
                    best = (q, t);
                }
            }
        }
        best
    }

    /// The closest point on the feature lines of segments `starts` (at
    /// least one), searching outward along them but never past a corner: a
    /// vertex stays on its own line, never jumping to a nearby one (the
    /// other side of a thin slot).
    fn project_feature(&self, p: Point3, starts: &[u32]) -> (Point3, u32) {
        let seg = |s: u32| {
            let (a, b) = self.segments[s as usize];
            closest_point_on_segment(p, self.points[a as usize], self.points[b as usize])
        };
        let home = starts.first().copied().unwrap_or(0);
        let mut best = (seg(home), home);
        let mut seen: HashSet<u32> = HashSet::new();
        let mut queue: VecDeque<u32> = VecDeque::new();
        for &s in starts {
            if seen.insert(s) {
                let q = seg(s);
                if distance_squared(p, q) < distance_squared(p, best.0) {
                    best = (q, s);
                }
                queue.push_back(s);
            }
        }
        while let Some(s) = queue.pop_front() {
            if seen.len() > 512 {
                break;
            }
            let (a, b) = self.segments[s as usize];
            for &v in &[a, b] {
                if !self.line[v as usize] {
                    continue;
                }
                for &t in &self.vertex_segments[v as usize] {
                    if seen.insert(t) {
                        let q = seg(t);
                        if distance_squared(p, q) < distance_squared(p, best.0) {
                            best = (q, t);
                        }
                        queue.push_back(t);
                    }
                }
            }
        }
        best
    }
}

fn closest(input: &Input<'_>, p: Point3, t: u32) -> Point3 {
    let [a, b, c] = input.triangles[t as usize].map(|v| input.points[v as usize]);
    closest_point_on_triangle(p, [a, b, c])
}

/// The surface being remeshed.
struct Mesh {
    pos: Vec<Point3>,
    kind: Vec<Kind>,
    /// Input triangle (smooth) or feature segment (feature, corner).
    home: Vec<u32>,
    tris: Vec<[u32; 3]>,
    alive: Vec<bool>,
    origin: Vec<u32>,
    /// Live triangles at each vertex.
    vt: Vec<Vec<u32>>,
    feature: HashSet<(u32, u32)>,
    /// Collapses that merged a feature too small to mesh.
    merged: usize,
}

impl Mesh {
    fn edge_tris(&self, u: u32, v: u32) -> Vec<u32> {
        self.vt[u as usize]
            .iter()
            .copied()
            .filter(|&t| self.tris[t as usize].contains(&v))
            .collect()
    }

    fn neighbours(&self, u: u32) -> Vec<u32> {
        let mut n: Vec<u32> = self.vt[u as usize]
            .iter()
            .flat_map(|&t| self.tris[t as usize])
            .filter(|&x| x != u)
            .collect();
        n.sort_unstable();
        n.dedup();
        n
    }

    fn tri_normal(&self, t: u32) -> Point3 {
        let [a, b, c] = self.tris[t as usize].map(|v| self.pos[v as usize]);
        normal(a, b, c)
    }

    fn vertex_normal(&self, u: u32) -> Point3 {
        let mut n = [0.0; 3];
        for &t in &self.vt[u as usize] {
            n = add(n, self.tri_normal(t));
        }
        unit(n)
    }

    /// The live edges, sorted.
    fn edges(&self) -> Vec<(u32, u32)> {
        let mut e: Vec<(u32, u32)> = (0..self.tris.len())
            .filter(|&t| self.alive[t])
            .flat_map(|t| {
                let [a, b, c] = self.tris[t];
                [key(a, b), key(b, c), key(c, a)]
            })
            .collect();
        e.sort_unstable();
        e.dedup();
        e
    }

    fn push_tri(&mut self, t: [u32; 3], origin: u32) -> u32 {
        let i = self.tris.len() as u32;
        self.tris.push(t);
        self.alive.push(true);
        self.origin.push(origin);
        for &v in &t {
            self.vt[v as usize].push(i);
        }
        i
    }

    fn unlink(&mut self, v: u32, t: u32) {
        self.vt[v as usize].retain(|&x| x != t);
    }

    /// Split edge `uv` at `m` (position already set).
    fn split(&mut self, u: u32, v: u32, m: u32) {
        for t in self.edge_tris(u, v) {
            let tri = self.tris[t as usize];
            // Rotate to [p, q, w] with pq the edge in this triangle's order.
            let r = (0..3)
                .find(|&k| key(tri[k], tri[(k + 1) % 3]) == key(u, v))
                .unwrap_or(0);
            let (p, q, w) = (tri[r], tri[(r + 1) % 3], tri[(r + 2) % 3]);
            self.tris[t as usize] = [p, m, w];
            self.unlink(q, t);
            self.vt[m as usize].push(t);
            let origin = self.origin[t as usize];
            self.push_tri([m, q, w], origin);
        }
        if self.feature.remove(&key(u, v)) {
            self.feature.insert(key(u, m));
            self.feature.insert(key(m, v));
        }
    }
}

/// Remesh a closed, outward-oriented, manifold surface so its edges follow
/// `size` (an edge of length L at a point of size h is split above 4h/3 and
/// collapsed below 4h/5). `groups` gives each triangle's group (empty: one
/// group); group boundaries are kept as feature lines.
///
/// # Errors
/// `Surface` for an input that is not a closed oriented manifold, `TooLarge`
/// past `max_triangles`, `Broken` if a pass broke the surface.
pub fn remesh(
    points: &[Point3],
    triangles: &[[u32; 3]],
    groups: &[u32],
    size: &dyn Fn(Point3) -> f64,
    params: &RemeshParams,
) -> Result<Remeshed, RemeshError> {
    let soup =
        TriangleSoup::new(points, triangles).map_err(|e| RemeshError::Surface(e.to_string()))?;
    let group = |t: usize| groups.get(t).copied().unwrap_or(0);
    let normals: Vec<Point3> = triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| points[v as usize]);
            unit(normal(a, b, c))
        })
        .collect();
    let mut vertex_tris = vec![Vec::new(); points.len()];
    for (t, tri) in triangles.iter().enumerate() {
        for &v in tri {
            vertex_tris[v as usize].push(t as u32);
        }
    }
    // Every edge must have exactly two triangles, in opposite directions.
    let mut directed: HashMap<(u32, u32), u32> = HashMap::new();
    for (t, tri) in triangles.iter().enumerate() {
        for k in 0..3 {
            if directed
                .insert((tri[k], tri[(k + 1) % 3]), t as u32)
                .is_some()
            {
                return Err(RemeshError::Surface(
                    "an edge is used twice in one direction".into(),
                ));
            }
        }
    }
    let cos_feature = libm::cos(params.feature_degrees.to_radians());
    let mut segments: Vec<(u32, u32)> = Vec::new();
    let mut sorted: Vec<(&(u32, u32), &u32)> = directed.iter().collect();
    sorted.sort_unstable();
    for (&(u, v), &t) in sorted {
        let Some(&s) = directed.get(&(v, u)) else {
            return Err(RemeshError::Surface("the surface is open".into()));
        };
        if u < v {
            let sharp = dot(normals[t as usize], normals[s as usize]) < cos_feature
                || normals[t as usize] == [0.0; 3]
                || normals[s as usize] == [0.0; 3];
            if sharp || group(t as usize) != group(s as usize) {
                segments.push((u, v));
            }
        }
    }
    let mut vertex_segments = vec![Vec::new(); points.len()];
    for (s, &(u, v)) in segments.iter().enumerate() {
        vertex_segments[u as usize].push(s as u32);
        vertex_segments[v as usize].push(s as u32);
    }
    // A vertex on one line through it, unless the line turns sharply there.
    let line: Vec<bool> = (0..points.len())
        .map(|v| {
            let segs = &vertex_segments[v];
            if segs.len() != 2 {
                return false;
            }
            let other = |s: u32| {
                let (a, b) = segments[s as usize];
                if a as usize == v { b } else { a }
            };
            let p = points[v];
            let d0 = unit(sub(points[other(segs[0]) as usize], p));
            let d1 = unit(sub(points[other(segs[1]) as usize], p));
            dot(d0, d1) < -cos_feature
        })
        .collect();
    let input = Input {
        points,
        triangles,
        normals,
        bvh: Bvh::new(soup),
        vertex_tris,
        segments,
        vertex_segments,
        line,
    };

    // The starting mesh: the input itself.
    let mut m = Mesh {
        pos: points.to_vec(),
        kind: Vec::with_capacity(points.len()),
        home: Vec::with_capacity(points.len()),
        tris: Vec::new(),
        alive: Vec::new(),
        origin: Vec::new(),
        vt: vec![Vec::new(); points.len()],
        feature: input.segments.iter().map(|&(u, v)| key(u, v)).collect(),
        merged: 0,
    };
    for v in 0..points.len() {
        let segs = &input.vertex_segments[v];
        let kind = if segs.is_empty() {
            Kind::Smooth
        } else if input.line[v] {
            Kind::Feature
        } else {
            Kind::Corner
        };
        m.kind.push(kind);
        m.home.push(match kind {
            Kind::Smooth => input.vertex_tris[v].first().copied().unwrap_or(0),
            _ => segs[0],
        });
    }
    for (t, tri) in triangles.iter().enumerate() {
        m.push_tri(*tri, t as u32);
    }

    let length = |m: &Mesh, u: u32, v: u32| norm(sub(m.pos[u as usize], m.pos[v as usize]));
    let target =
        |m: &Mesh, u: u32, v: u32| size(scale(add(m.pos[u as usize], m.pos[v as usize]), 0.5));

    for _ in 0..params.iterations {
        // 1. Split long edges, longest first.
        let mut long: Vec<(f64, (u32, u32))> = m
            .edges()
            .into_iter()
            .filter_map(|(u, v)| {
                let l = length(&m, u, v);
                (l > 4.0 / 3.0 * target(&m, u, v)).then_some((l, (u, v)))
            })
            .collect();
        long.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        for (_, (u, v)) in long {
            if m.edge_tris(u, v).len() != 2 {
                continue;
            }
            let live = m.alive.iter().filter(|&&a| a).count();
            if live + 2 > params.max_triangles {
                return Err(RemeshError::TooLarge {
                    triangles: live + 2,
                    limit: params.max_triangles,
                });
            }
            let mid = scale(add(m.pos[u as usize], m.pos[v as usize]), 0.5);
            let i = m.pos.len() as u32;
            let feature = m.feature.contains(&key(u, v));
            let (p, kind, home) = if feature {
                // From a feature end's segment, or from every segment at
                // a corner (corners are input vertices).
                let starts = if m.kind[u as usize] == Kind::Feature {
                    vec![m.home[u as usize]]
                } else if m.kind[v as usize] == Kind::Feature {
                    vec![m.home[v as usize]]
                } else {
                    input.vertex_segments[u as usize].clone()
                };
                let (q, s) = input.project_feature(mid, &starts);
                (q, Kind::Feature, s)
            } else {
                let tris = m.edge_tris(u, v);
                let n = unit(add(m.tri_normal(tris[0]), m.tri_normal(tris[1])));
                let (q, t) = input.project(mid, n, m.origin[tris[0] as usize]);
                (q, Kind::Smooth, t)
            };
            m.pos.push(p);
            m.kind.push(kind);
            m.home.push(home);
            m.vt.push(Vec::new());
            m.split(u, v, i);
        }

        // 2. Collapse short edges, shortest first.
        let mut short: Vec<(f64, (u32, u32))> = m
            .edges()
            .into_iter()
            .filter_map(|(u, v)| {
                let l = length(&m, u, v);
                (l < 0.8 * target(&m, u, v)).then_some((l, (u, v)))
            })
            .collect();
        short.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        for (_, (u, v)) in short {
            if m.edge_tris(u, v).len() != 2 {
                continue;
            }
            let h = target(&m, u, v);
            let small = params.small_feature * h;
            if !(try_collapse(&mut m, u, v, h, small) || try_collapse(&mut m, v, u, h, small)) {
                continue;
            }
        }

        // 3. Flip edges to Delaunay.
        for (u, v) in m.edges() {
            try_flip(&mut m, u, v);
        }

        // 4. Tangential smoothing, projected back onto the input.
        let order: Vec<u32> = (0..m.pos.len() as u32)
            .filter(|&v| !m.vt[v as usize].is_empty() && m.kind[v as usize] != Kind::Corner)
            .collect();
        for v in order {
            smooth(&mut m, &input, v);
        }
    }

    finish(&m)
}

/// Collapse `x` into `y` if that keeps the surface valid and its edges no
/// longer than 4h/3. A corner collapses only into another corner along a
/// feature edge shorter than `small`; a feature vertex along its line, or
/// into a feature vertex or corner of another line closer than `small`.
fn try_collapse(m: &mut Mesh, x: u32, y: u32, h: f64, small: f64) -> bool {
    let (kx, ky) = (m.kind[x as usize], m.kind[y as usize]);
    match kx {
        Kind::Corner => {
            if ky != Kind::Corner
                || !m.feature.contains(&key(x, y))
                || norm(sub(m.pos[x as usize], m.pos[y as usize])) >= small
            {
                return false;
            }
        }
        Kind::Feature => {
            // Along its line, or across to another line closer than
            // `small` (the two sides of a strip too narrow to mesh).
            if ky == Kind::Smooth
                || !(m.feature.contains(&key(x, y))
                    || norm(sub(m.pos[x as usize], m.pos[y as usize])) < small)
            {
                return false;
            }
        }
        Kind::Smooth => {}
    }
    let edge = m.edge_tris(x, y);
    if edge.len() != 2 {
        return false;
    }
    let opposite: Vec<u32> = edge
        .iter()
        .map(|&t| {
            m.tris[t as usize]
                .iter()
                .copied()
                .find(|&w| w != x && w != y)
                .unwrap_or(x)
        })
        .collect();
    // Link condition: x and y share exactly the two opposite vertices.
    let (nx, ny) = (m.neighbours(x), m.neighbours(y));
    let common: Vec<u32> = nx.iter().copied().filter(|w| ny.contains(w)).collect();
    if common.len() != 2 || !opposite.iter().all(|w| common.contains(w)) {
        return false;
    }
    // Opposite vertices keep a valence of at least 3.
    if opposite.iter().any(|&w| m.neighbours(w).len() <= 3) {
        return false;
    }
    let py = m.pos[y as usize];
    for &t in &m.vt[x as usize] {
        if edge.contains(&t) {
            continue;
        }
        let tri = m.tris[t as usize];
        let old = m.tri_normal(t);
        let moved = tri.map(|w| if w == x { py } else { m.pos[w as usize] });
        let new = normal(moved[0], moved[1], moved[2]);
        // No flipped or collapsed triangle, no edge longer than allowed.
        if dot(unit(old), unit(new)) < 0.5 || norm(new) <= 0.0 {
            return false;
        }
        for &w in &tri {
            if w != x && norm(sub(m.pos[w as usize], py)) > 4.0 / 3.0 * h {
                return false;
            }
        }
    }
    // Do it.
    if kx == Kind::Corner || (kx == Kind::Feature && !m.feature.contains(&key(x, y))) {
        m.merged += 1;
    }
    for &t in &edge {
        m.alive[t as usize] = false;
        for w in m.tris[t as usize] {
            m.unlink(w, t);
        }
    }
    let around = std::mem::take(&mut m.vt[x as usize]);
    for t in around {
        for w in &mut m.tris[t as usize] {
            if *w == x {
                *w = y;
            }
        }
        m.vt[y as usize].push(t);
    }
    let rekey: Vec<(u32, u32)> = m
        .feature
        .iter()
        .copied()
        .filter(|&(a, b)| a == x || b == x)
        .collect();
    for (a, b) in rekey {
        m.feature.remove(&(a, b));
        let other = if a == x { b } else { a };
        if other != y {
            m.feature.insert(key(other, y));
        }
    }
    true
}

/// Flip edge `uv` if it is not locally Delaunay (its opposite angles sum
/// past 180°) and the two triangles are nearly coplanar. Delaunay flips
/// rather than the usual valence flips: the volume kernel recovers a surface
/// whose triangles are (close to) Delaunay with little work, while
/// valence-optimal triangles are often missing from the tetrahedralization
/// of their vertices and cost a cascade of Steiner points.
fn try_flip(m: &mut Mesh, u: u32, v: u32) -> bool {
    if m.feature.contains(&key(u, v)) {
        return false;
    }
    let tris = m.edge_tris(u, v);
    let &[t1, t2] = tris.as_slice() else {
        return false;
    };
    // Orient: t1 holds u → v.
    let has = |t: u32, a: u32, b: u32| {
        let x = m.tris[t as usize];
        (0..3).any(|k| x[k] == a && x[(k + 1) % 3] == b)
    };
    let (t1, t2) = if has(t1, u, v) { (t1, t2) } else { (t2, t1) };
    if !has(t1, u, v) || !has(t2, v, u) {
        return false;
    }
    let third = |t: u32| {
        m.tris[t as usize]
            .iter()
            .copied()
            .find(|&w| w != u && w != v)
            .unwrap_or(u)
    };
    let (a, b) = (third(t1), third(t2));
    if a == b || m.neighbours(a).contains(&b) {
        return false;
    }
    if m.neighbours(u).len() <= 3 || m.neighbours(v).len() <= 3 {
        return false;
    }
    let angle = |at: u32, x: u32, y: u32| {
        let p = m.pos[at as usize];
        let (d0, d1) = (
            unit(sub(m.pos[x as usize], p)),
            unit(sub(m.pos[y as usize], p)),
        );
        libm::acos(dot(d0, d1).clamp(-1.0, 1.0))
    };
    // A small margin keeps cocircular quads from flipping back and forth.
    if angle(a, u, v) + angle(b, u, v) <= std::f64::consts::PI + 1e-9 {
        return false;
    }
    let (n1, n2) = (unit(m.tri_normal(t1)), unit(m.tri_normal(t2)));
    if dot(n1, n2) < 0.94 {
        return false;
    }
    let p = |w: u32| m.pos[w as usize];
    let (na, nb) = (
        unit(normal(p(a), p(u), p(b))),
        unit(normal(p(b), p(v), p(a))),
    );
    if [na, nb]
        .iter()
        .any(|&n| dot(n, n1) < 0.5 || dot(n, n2) < 0.5)
    {
        return false;
    }
    m.tris[t1 as usize] = [a, u, b];
    m.tris[t2 as usize] = [b, v, a];
    m.unlink(u, t2);
    m.unlink(v, t1);
    m.vt[a as usize].push(t2);
    m.vt[b as usize].push(t1);
    true
}

/// Move a vertex halfway to the centroid of its neighbours, in its tangent
/// plane (or along its feature line), then back onto the input; undone if
/// any triangle around it would turn over.
fn smooth(m: &mut Mesh, input: &Input<'_>, v: u32) {
    let p = m.pos[v as usize];
    let (target, home) = match m.kind[v as usize] {
        Kind::Corner => return,
        Kind::Feature => {
            let along: Vec<u32> = m
                .neighbours(v)
                .into_iter()
                .filter(|&w| m.feature.contains(&key(v, w)))
                .collect();
            if along.len() != 2 {
                return;
            }
            let c = scale(add(m.pos[along[0] as usize], m.pos[along[1] as usize]), 0.5);
            input.project_feature(add(p, scale(sub(c, p), 0.5)), &[m.home[v as usize]])
        }
        Kind::Smooth => {
            let ring = m.neighbours(v);
            if ring.is_empty() {
                return;
            }
            let mut c = [0.0; 3];
            for &w in &ring {
                c = add(c, m.pos[w as usize]);
            }
            c = scale(c, 1.0 / ring.len() as f64);
            let n = m.vertex_normal(v);
            let d = sub(c, p);
            let tangential = sub(d, scale(n, dot(d, n)));
            input.project(add(p, scale(tangential, 0.5)), n, m.home[v as usize])
        }
    };
    let before: Vec<Point3> = m.vt[v as usize]
        .iter()
        .map(|&t| unit(m.tri_normal(t)))
        .collect();
    m.pos[v as usize] = target;
    let ok = m.vt[v as usize].iter().zip(&before).all(|(&t, &n)| {
        let after = m.tri_normal(t);
        norm(after) > 0.0 && dot(unit(after), n) > 0.5
    });
    if ok {
        m.home[v as usize] = home;
    } else {
        m.pos[v as usize] = p;
    }
}

/// Compact the live mesh and check it is closed and consistently oriented.
fn finish(m: &Mesh) -> Result<Remeshed, RemeshError> {
    let mut index = vec![u32::MAX; m.pos.len()];
    let mut points = Vec::new();
    for (v, around) in m.vt.iter().enumerate() {
        if !around.is_empty() {
            index[v] = points.len() as u32;
            points.push(m.pos[v]);
        }
    }
    let mut triangles = Vec::new();
    let mut origin = Vec::new();
    for t in 0..m.tris.len() {
        if m.alive[t] {
            triangles.push(m.tris[t].map(|v| index[v as usize]));
            origin.push(m.origin[t]);
        }
    }
    let mut directed: HashSet<(u32, u32)> = HashSet::new();
    for t in &triangles {
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            return Err(RemeshError::Broken("a degenerate triangle".into()));
        }
        for k in 0..3 {
            if !directed.insert((t[k], t[(k + 1) % 3])) {
                return Err(RemeshError::Broken(
                    "an edge used twice in one direction".into(),
                ));
            }
        }
    }
    if directed.iter().any(|&(a, b)| !directed.contains(&(b, a))) {
        return Err(RemeshError::Broken("an open edge".into()));
    }
    Ok(Remeshed {
        points,
        triangles,
        origin,
        merged_small_features: m.merged,
    })
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;

    fn boxed(lo: Point3, hi: Point3) -> (Vec<Point3>, Vec<[u32; 3]>) {
        let p: Vec<Point3> = (0..8)
            .map(|i| {
                [
                    if i & 1 == 0 { lo[0] } else { hi[0] },
                    if i & 2 == 0 { lo[1] } else { hi[1] },
                    if i & 4 == 0 { lo[2] } else { hi[2] },
                ]
            })
            .collect();
        let f = vec![
            [0, 2, 3],
            [0, 3, 1],
            [4, 5, 7],
            [4, 7, 6],
            [0, 1, 5],
            [0, 5, 4],
            [2, 6, 7],
            [2, 7, 3],
            [0, 4, 6],
            [0, 6, 2],
            [1, 3, 7],
            [1, 7, 5],
        ];
        (p, f)
    }

    fn volume(p: &[Point3], f: &[[u32; 3]]) -> f64 {
        f.iter()
            .map(|t| {
                let [a, b, c] = t.map(|v| p[v as usize]);
                dot(a, cross(b, c)) / 6.0
            })
            .sum()
    }

    fn edge_lengths(r: &Remeshed) -> Vec<f64> {
        let mut e: Vec<(u32, u32)> = r
            .triangles
            .iter()
            .flat_map(|t| [key(t[0], t[1]), key(t[1], t[2]), key(t[2], t[0])])
            .collect();
        e.sort_unstable();
        e.dedup();
        e.iter()
            .map(|&(a, b)| norm(sub(r.points[a as usize], r.points[b as usize])))
            .collect()
    }

    #[test]
    fn a_box_keeps_its_shape_and_corners() {
        let (p, f) = boxed([0.0; 3], [1.0, 2.0, 0.5]);
        let r = remesh(&p, &f, &[], &|_| 0.2, &RemeshParams::default()).unwrap();
        assert!(r.triangles.len() > 100, "{}", r.triangles.len());
        // Flat faces and sharp edges: the volume is exact to rounding.
        assert!((volume(&r.points, &r.triangles) - 1.0).abs() < 1e-9);
        for corner in &p {
            assert!(r.points.contains(corner), "lost corner {corner:?}");
        }
        let lengths = edge_lengths(&r);
        // Smoothing comes last, so a few edges near the sharp edges end
        // somewhat long.
        let max = lengths.iter().copied().fold(0.0, f64::max);
        assert!(max < 1.8 * 0.2, "{max}");
        let within = lengths
            .iter()
            .filter(|&&l| (0.4 * 0.2..=1.4 * 0.2).contains(&l))
            .count();
        assert!(
            within as f64 > 0.95 * lengths.len() as f64,
            "{within} of {}",
            lengths.len()
        );
        let mean = lengths.iter().sum::<f64>() / lengths.len() as f64;
        assert!((0.12..0.3).contains(&mean), "{mean}");
    }

    /// The unit cube with the edge at x = 1, z = 1 cut by a 45° chamfer `w`
    /// wide, each convex face triangulated as a fan and oriented outwards.
    fn chamfered_cube(w: f64) -> (Vec<Point3>, Vec<[u32; 3]>) {
        let mut p = Vec::new();
        for y in [0.0, 1.0] {
            p.extend([
                [0.0, y, 0.0],
                [1.0, y, 0.0],
                [1.0, y, 1.0 - w],
                [1.0 - w, y, 1.0],
                [0.0, y, 1.0],
            ]);
        }
        let faces: [&[u32]; 7] = [
            &[0, 1, 2, 3, 4],
            &[5, 6, 7, 8, 9],
            &[0, 4, 9, 5],
            &[0, 1, 6, 5],
            &[1, 2, 7, 6],
            &[2, 3, 8, 7],
            &[3, 4, 9, 8],
        ];
        let centre = [0.5; 3];
        let mut f = Vec::new();
        for poly in faces {
            for k in 1..poly.len() - 1 {
                let t = [poly[0], poly[k], poly[k + 1]];
                let [a, b, c] = t.map(|v| p[v as usize]);
                f.push(if dot(normal(a, b, c), sub(a, centre)) > 0.0 {
                    t
                } else {
                    [t[0], t[2], t[1]]
                });
            }
        }
        (p, f)
    }

    fn smallest_angle(r: &Remeshed) -> f64 {
        let angle = |a: Point3, b: Point3, c: Point3| {
            let (u, w) = (sub(b, a), sub(c, a));
            libm::acos((dot(u, w) / (norm(u) * norm(w))).clamp(-1.0, 1.0)).to_degrees()
        };
        r.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|v| r.points[v as usize]);
                angle(a, b, c).min(angle(b, c, a)).min(angle(c, a, b))
            })
            .fold(180.0, f64::min)
    }

    /// A chamfer far narrower than the size is merged away (its two lines
    /// into one, its end corners into each other), so no needle triangles
    /// run along it; the surface moves by at most the chamfer's width.
    #[test]
    fn a_feature_too_small_to_mesh_is_merged() {
        let w = 0.005;
        let (p, f) = chamfered_cube(w);
        assert!((volume(&p, &f) - (1.0 - w * w / 2.0)).abs() < 1e-12);
        let kept = remesh(
            &p,
            &f,
            &[],
            &|_| 0.2,
            &RemeshParams {
                small_feature: 0.0,
                ..RemeshParams::default()
            },
        )
        .unwrap();
        assert_eq!(kept.merged_small_features, 0);
        assert!(smallest_angle(&kept) < 5.0, "{}", smallest_angle(&kept));
        let r = remesh(&p, &f, &[], &|_| 0.2, &RemeshParams::default()).unwrap();
        assert!(r.merged_small_features > 0);
        assert!(smallest_angle(&r) > 20.0, "{}", smallest_angle(&r));
        // A merged corner tilts the faces at it by up to `w`: the volume
        // moves by at most w times their area (here one unit face).
        let moved = (volume(&r.points, &r.triangles) - volume(&p, &f)).abs();
        assert!(moved <= w, "{moved}");
        // The six corners off the chamfer are kept.
        for corner in [0, 1, 4, 5, 6, 9] {
            assert!(r.points.contains(&p[corner]), "lost corner {corner}");
        }
    }

    #[test]
    fn group_boundaries_are_kept() {
        // A box split into two groups at its top face: the top's boundary
        // edges must stay edges of the result, and each triangle keeps the
        // group of its origin.
        let (p, f) = boxed([0.0; 3], [1.0; 3]);
        let groups: Vec<u32> = (0..12).map(|t| u32::from(t == 2 || t == 3)).collect();
        let r = remesh(&p, &f, &groups, &|_| 0.15, &RemeshParams::default()).unwrap();
        for (t, &o) in r.triangles.iter().zip(&r.origin) {
            let z: f64 = t.iter().map(|&v| r.points[v as usize][2]).sum::<f64>() / 3.0;
            let top = (z - 1.0).abs() < 1e-12;
            assert_eq!(groups[o as usize] == 1, top, "triangle {t:?} from {o}");
        }
    }

    /// A UV sphere, coordinates rounded to six decimals.
    fn sphere(segments: u32, rings: u32) -> (Vec<Point3>, Vec<[u32; 3]>) {
        let round = |x: f64| (x * 1e6).round() / 1e6;
        let mut p = vec![[0.0, 0.0, 1.0]];
        for i in 1..rings {
            let t = f64::from(i) / f64::from(rings) * std::f64::consts::PI;
            for j in 0..segments {
                let a = f64::from(j) / f64::from(segments) * std::f64::consts::TAU;
                p.push([
                    round(libm::sin(t) * libm::cos(a)),
                    round(libm::sin(t) * libm::sin(a)),
                    round(libm::cos(t)),
                ]);
            }
        }
        let south = p.len() as u32;
        p.push([0.0, 0.0, -1.0]);
        let ring = |i: u32, j: u32| 1 + i * segments + j % segments;
        let mut f = Vec::new();
        for j in 0..segments {
            f.push([0, ring(0, j), ring(0, j + 1)]);
            f.push([south, ring(rings - 2, j + 1), ring(rings - 2, j)]);
        }
        for i in 0..rings - 2 {
            for j in 0..segments {
                f.push([ring(i, j), ring(i + 1, j), ring(i + 1, j + 1)]);
                f.push([ring(i, j), ring(i + 1, j + 1), ring(i, j + 1)]);
            }
        }
        (p, f)
    }

    #[test]
    fn a_sphere_stays_on_its_surface() {
        let (p, f) = sphere(16, 10);
        let r = remesh(&p, &f, &[], &|_| 0.15, &RemeshParams::default()).unwrap();
        let v0 = volume(&p, &f);
        let v1 = volume(&r.points, &r.triangles);
        assert!((v1 - v0).abs() / v0 < 0.02, "{v0} -> {v1}");
        for q in &r.points {
            let d = norm(*q);
            assert!((0.95..=1.0 + 1e-6).contains(&d), "{q:?} off the sphere");
        }
        // Thin pole triangles are gone: the smallest angle is reasonable.
        let worst = r
            .triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|v| r.points[v as usize]);
                let angle = |p: Point3, q: Point3, s: Point3| {
                    libm::acos(dot(unit(sub(q, p)), unit(sub(s, p))).clamp(-1.0, 1.0)).to_degrees()
                };
                angle(a, b, c).min(angle(b, c, a)).min(angle(c, a, b))
            })
            .fold(180.0, f64::min);
        assert!(worst > 15.0, "{worst}");
    }

    #[test]
    fn open_surfaces_are_refused() {
        let (p, mut f) = boxed([0.0; 3], [1.0; 3]);
        f.pop();
        assert!(matches!(
            remesh(&p, &f, &[], &|_| 0.2, &RemeshParams::default()),
            Err(RemeshError::Surface(_))
        ));
    }

    /// Two feature lines leave one corner 3° apart (the two walls of a thin
    /// slot): a vertex on the first, projected from near the second, stays
    /// on the first; it never jumps across past the corner.
    #[test]
    fn a_feature_vertex_stays_on_its_own_line() {
        let points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 0.05, 0.0],
            [2.0, 0.1, 0.0],
        ];
        let triangles = vec![[0, 1, 3]];
        let segments = vec![(0, 1), (1, 2), (0, 3), (3, 4)];
        let mut vertex_segments = vec![Vec::new(); points.len()];
        for (s, &(a, b)) in segments.iter().enumerate() {
            vertex_segments[a as usize].push(s as u32);
            vertex_segments[b as usize].push(s as u32);
        }
        let input = Input {
            points: &points,
            triangles: &triangles,
            normals: vec![[0.0, 0.0, 1.0]],
            bvh: Bvh::new(TriangleSoup::new(&points, &triangles).unwrap()),
            vertex_tris: vec![vec![0], vec![0], Vec::new(), vec![0], Vec::new()],
            segments,
            vertex_segments,
            // Vertices 1 and 3 are on straight lines; 0 is a corner.
            line: vec![false, true, false, true, false],
        };
        let (q, s) = input.project_feature([1.5, 0.08, 0.0], &[1]);
        assert!(s <= 1, "jumped to segment {s}");
        assert!(distance_squared(q, [1.5, 0.0, 0.0]) < 1e-24, "{q:?}");
        // From the corner both lines are open.
        let (_, s) = input.project_feature([1.5, 0.08, 0.0], &[0, 2]);
        assert_eq!(s, 3);
    }
}
