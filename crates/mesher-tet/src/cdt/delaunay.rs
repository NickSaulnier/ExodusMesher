//! An exact incremental Delaunay tetrahedralization over implicit points
//! (`mesher_geom::implicit`): the base of the exact constrained kernel.
//!
//! The same design as `crate::delaunay`: a flat array of tets, `n[k]` across
//! the face opposite `v[k]`, ghost tets sharing a vertex at infinity
//! ([`GHOST`], always last) outside the convex hull, and real tets positively
//! oriented. Inserting a point walks to it, grows the cavity of tets whose
//! circumsphere strictly contains it (for a ghost: the point is beyond the
//! hull face, or on its plane inside its circumcircle), grows it further
//! across any boundary face the point does not strictly see, which settles
//! cospherical and coplanar points, and fills it with a cone from the point.
//!
//! Every decision is an exact implicit predicate, so points may be explicit
//! or constructed (LNC, BPT, LPI, TPI) and are never rounded. f64
//! approximations are used only to order insertions.

use mesher_geom::Sign;
use mesher_geom::implicit::{self, Point};

/// The vertex at infinity.
pub const GHOST: u32 = u32::MAX;
const NONE: u32 = u32::MAX;

/// Local face opposite each vertex, ordered outward for a positive tet.
const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

/// The six edges of a tet, as local vertex pairs.
const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// What a segment from a vertex first meets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Crossing {
    /// The interior of face `j` of tet `t` (the face opposite the start).
    Face { t: u32, j: usize },
    /// The interior of edge `xy`.
    Edge { x: u32, y: u32 },
    /// A vertex on the segment, or nothing found.
    None,
}

/// What is wrong with a constrained cavity.
enum Violation {
    None,
    /// Take this tet out.
    Remove(u32),
    /// Only a tet containing the point could be taken out.
    Stuck,
}

/// A face's vertices rotated so the smallest is first, keeping orientation.
fn rotated(f: [u32; 3]) -> [u32; 3] {
    let m = (0..3).min_by_key(|&k| f[k]).unwrap_or(0);
    [f[m], f[(m + 1) % 3], f[(m + 2) % 3]]
}

#[derive(Debug, Clone, Copy)]
struct Tet {
    v: [u32; 4],
    n: [u32; 4],
    alive: bool,
}

/// Why a tetrahedralization failed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DelaunayError {
    /// Fewer than four points, or all of them coplanar.
    Flat,
    /// A cavity could not be filled (a bug or a pathological input).
    Cavity { point: u32 },
    /// A protected insertion outside the convex hull.
    Outside { point: u32 },
}

impl std::fmt::Display for DelaunayError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DelaunayError::Flat => write!(f, "the points do not span a volume"),
            DelaunayError::Cavity { point } => {
                write!(f, "could not fill the cavity of point {point}")
            }
            DelaunayError::Outside { point } => {
                write!(f, "point {point} is outside the convex hull")
            }
        }
    }
}

/// Where an inserted point went.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Inserted {
    /// A new vertex with this index.
    New(u32),
    /// Exactly equal to this existing vertex; nothing changed.
    Duplicate(u32),
}

/// The exact Delaunay tetrahedralization.
#[derive(Debug, Clone)]
pub struct ExactDelaunay {
    points: Vec<Point>,
    tets: Vec<Tet>,
    free: Vec<u32>,
    alive: usize,
    /// Insertion stamp per tet (cavity membership).
    stamp: Vec<u32>,
    epoch: u32,
    /// A tet to start the next walk from.
    last: u32,
    rng: u64,
    /// `(point, earlier equal vertex)` for exact duplicates among the inputs.
    duplicates: Vec<(u32, u32)>,
    /// A live tet holding each vertex (`NONE` for a point that is not one).
    vertex_tet: Vec<u32>,
}

/// The face opposite local vertex `k`, outward for a positive tet.
fn face(v: &[u32; 4], k: usize) -> [u32; 3] {
    FACES[k].map(|i| v[i])
}

fn sorted3(mut f: [u32; 3]) -> [u32; 3] {
    f.sort_unstable();
    f
}

impl ExactDelaunay {
    /// Tetrahedralize `points`, inserted in biased randomized order from
    /// `seed` (sorted along a Morton curve by their approximate positions).
    /// Vertex `i` is `points[i]`; exact duplicates are not vertices (see
    /// [`ExactDelaunay::duplicates`]).
    ///
    /// # Errors
    /// [`DelaunayError::Flat`] when the points do not span a volume, or
    /// [`DelaunayError::Cavity`] if an insertion fails.
    pub fn new(points: Vec<Point>, seed: u64) -> Result<ExactDelaunay, DelaunayError> {
        let approx: Vec<_> = points.iter().map(Point::approx).collect();
        let order = crate::delaunay::brio_order(&approx, seed);
        let mut d = ExactDelaunay {
            points,
            tets: Vec::new(),
            free: Vec::new(),
            alive: 0,
            stamp: Vec::new(),
            epoch: 0,
            last: 0,
            rng: seed | 1,
            duplicates: Vec::new(),
            vertex_tet: Vec::new(),
        };
        d.vertex_tet = vec![NONE; d.points.len()];
        let first = d.initial_tet(&order)?;
        for &i in &order {
            if first.contains(&i) {
                continue;
            }
            if let Inserted::Duplicate(v) = d.insert_index(i)? {
                d.duplicates.push((i, v));
            }
        }
        d.duplicates.sort_unstable();
        Ok(d)
    }

    /// Insert a further point.
    ///
    /// # Errors
    /// [`DelaunayError::Cavity`] if the cavity could not be filled; nothing
    /// changes then.
    pub fn insert(&mut self, p: Point) -> Result<Inserted, DelaunayError> {
        let i = self.points.len() as u32;
        self.points.push(p);
        self.vertex_tet.push(NONE);
        let result = self.insert_index(i);
        if !matches!(result, Ok(Inserted::New(_))) {
            self.points.pop();
            self.vertex_tet.pop();
        }
        result
    }

    #[must_use]
    pub fn points(&self) -> &[Point] {
        &self.points
    }

    /// `(point, equal earlier vertex)` for exact duplicates among the inputs.
    #[must_use]
    pub fn duplicates(&self) -> &[(u32, u32)] {
        &self.duplicates
    }

    /// The real (finite) tets, in slot order.
    #[must_use]
    pub fn tets(&self) -> Vec<[u32; 4]> {
        self.tets
            .iter()
            .filter(|t| t.alive && t.v[3] != GHOST)
            .map(|t| t.v)
            .collect()
    }

    /// The live tets (ghosts included) around vertex `v`, found from its
    /// stored tet across the faces that contain it.
    #[must_use]
    pub fn star(&self, v: u32) -> Vec<u32> {
        let Some(&start) = self.vertex_tet.get(v as usize) else {
            return Vec::new();
        };
        if start == NONE || !self.tets.get(start as usize).is_some_and(|t| t.alive) {
            return Vec::new();
        }
        let mut out = vec![start];
        let mut k = 0;
        while k < out.len() {
            let t = self.tets[out[k] as usize];
            k += 1;
            for j in 0..4 {
                // Across face j, which contains v unless v is opposite it.
                if t.v[j] != v {
                    let u = t.n[j];
                    if !out.contains(&u) {
                        out.push(u);
                    }
                }
            }
        }
        out
    }

    /// The vertices joined to `v` by an edge (not the vertex at infinity),
    /// sorted.
    #[must_use]
    pub fn neighbours(&self, v: u32) -> Vec<u32> {
        let mut out: Vec<u32> = self
            .star(v)
            .iter()
            .flat_map(|&t| self.tets[t as usize].v)
            .filter(|&x| x != v && x != GHOST)
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Is `uv` an edge of the tetrahedralization?
    #[must_use]
    pub fn has_edge(&self, u: u32, v: u32) -> bool {
        self.star(u)
            .iter()
            .any(|&t| self.tets[t as usize].v.contains(&v))
    }

    /// Is `f` (any vertex order) a face?
    #[must_use]
    pub fn has_face(&self, f: [u32; 3]) -> bool {
        self.star(f[0])
            .iter()
            .any(|&t| f.iter().all(|x| self.tets[t as usize].v.contains(x)))
    }

    fn pt(&self, v: u32) -> &Point {
        &self.points[v as usize]
    }

    fn orient(&self, f: [u32; 3], p: &Point) -> Sign {
        implicit::orient3d(self.pt(f[0]), self.pt(f[1]), self.pt(f[2]), p)
    }

    fn alloc(&mut self, v: [u32; 4]) -> u32 {
        let t = Tet {
            v,
            n: [NONE; 4],
            alive: true,
        };
        self.alive += 1;
        let id = if let Some(id) = self.free.pop() {
            self.tets[id as usize] = t;
            id
        } else {
            self.tets.push(t);
            self.stamp.push(0);
            (self.tets.len() - 1) as u32
        };
        for x in v {
            if let Some(slot) = self.vertex_tet.get_mut(x as usize) {
                *slot = id;
            }
        }
        id
    }

    fn kill(&mut self, t: u32) {
        self.tets[t as usize].alive = false;
        self.alive -= 1;
        self.free.push(t);
    }

    fn is_ghost(&self, t: u32) -> bool {
        self.tets[t as usize].v[3] == GHOST
    }

    /// The first four affinely independent points (in insertion order) form
    /// the first tet, closed off by four ghosts.
    fn initial_tet(&mut self, order: &[u32]) -> Result<[u32; 4], DelaunayError> {
        let a = *order.first().ok_or(DelaunayError::Flat)?;
        let b = *order
            .iter()
            .find(|&&i| !implicit::coincident(self.pt(a), self.pt(i)))
            .ok_or(DelaunayError::Flat)?;
        let c = *order
            .iter()
            .find(|&&i| !implicit::collinear(self.pt(a), self.pt(b), self.pt(i)))
            .ok_or(DelaunayError::Flat)?;
        let (d, s) = order
            .iter()
            .map(|&i| (i, self.orient([a, b, c], self.pt(i))))
            .find(|&(_, s)| s != Sign::Zero)
            .ok_or(DelaunayError::Flat)?;
        let v = if s == Sign::Positive {
            [a, b, c, d]
        } else {
            [b, a, c, d]
        };
        let t0 = self.alloc(v);
        let mut ghosts = Vec::with_capacity(4);
        for k in 0..4 {
            let f = face(&v, k);
            let g = self.alloc([f[0], f[1], f[2], GHOST]);
            self.tets[t0 as usize].n[k] = g;
            self.tets[g as usize].n[3] = t0;
            ghosts.push(g);
        }
        self.link_around(&ghosts, GHOST)
            .map_err(|()| DelaunayError::Cavity { point: d })?;
        self.last = t0;
        Ok([a, b, c, d])
    }

    /// Pair up the faces containing vertex `s` among `ids`.
    fn link_around(&mut self, ids: &[u32], s: u32) -> Result<(), ()> {
        let mut keys: Vec<([u32; 3], u32, usize)> = Vec::with_capacity(ids.len() * 3);
        for &t in ids {
            let v = self.tets[t as usize].v;
            for k in 0..4 {
                if v[k] != s {
                    keys.push((sorted3(face(&v, k)), t, k));
                }
            }
        }
        keys.sort_unstable();
        if keys.len() % 2 == 1 {
            return Err(());
        }
        for pair in keys.chunks(2) {
            let [(fa, ta, ka), (fb, tb, kb)] = [pair[0], pair[1]];
            if fa != fb {
                return Err(());
            }
            self.tets[ta as usize].n[ka] = tb;
            self.tets[tb as usize].n[kb] = ta;
        }
        Ok(())
    }

    /// Does tet `t` conflict with `p`? A real tet when `p` is strictly inside
    /// its circumsphere; a ghost when `p` is beyond its hull face, or on the
    /// face's plane strictly inside its circumcircle.
    fn conflict(&self, t: u32, p: &Point) -> bool {
        let v = self.tets[t as usize].v;
        let [a, b, c] = [v[0], v[1], v[2]].map(|x| self.pt(x));
        if v[3] == GHOST {
            match implicit::orient3d(a, b, c, p) {
                Sign::Positive => true,
                Sign::Negative => false,
                Sign::Zero => implicit::incircle_coplanar(a, b, c, p) == Sign::Positive,
            }
        } else {
            implicit::insphere(a, b, c, self.pt(v[3]), p) == Sign::Positive
        }
    }

    /// Does `p` strictly see face `j` of cavity tet `t` from inside, so that
    /// the cone tet on it is positive? A ghost face (two real vertices) is
    /// seen when `p` is off their line.
    fn sees(&self, t: u32, j: usize, p: &Point) -> bool {
        let f = face(&self.tets[t as usize].v, j);
        if f.contains(&GHOST) {
            let r: Vec<u32> = f.iter().copied().filter(|&x| x != GHOST).collect();
            !implicit::collinear(self.pt(r[0]), self.pt(r[1]), p)
        } else {
            self.orient(f, p) == Sign::Negative
        }
    }

    fn next_random(&mut self) -> u64 {
        self.rng ^= self.rng << 13;
        self.rng ^= self.rng >> 7;
        self.rng ^= self.rng << 17;
        self.rng
    }

    /// A real tet whose closure contains `p`, or a ghost whose hull face `p`
    /// is beyond.
    fn locate(&mut self, p: &Point) -> u32 {
        let mut t = self.last;
        if !self.tets.get(t as usize).is_some_and(|x| x.alive) {
            t = self.tets.iter().position(|x| x.alive).unwrap_or(0) as u32;
        }
        if self.is_ghost(t) {
            t = self.tets[t as usize].n[3];
        }
        let limit = 4 * self.tets.len() + 64;
        for _ in 0..limit {
            if self.is_ghost(t) {
                return t;
            }
            let v = self.tets[t as usize].v;
            let start = (self.next_random() % 4) as usize;
            let next = (0..4)
                .map(|k| (start + k) % 4)
                .find(|&i| self.orient(face(&v, i), p) == Sign::Positive);
            match next {
                Some(i) => t = self.tets[t as usize].n[i],
                None => return t,
            }
        }
        // The walk always ends; scan if it did not.
        (0..self.tets.len() as u32)
            .find(|&t| {
                let x = self.tets[t as usize];
                x.alive
                    && (if x.v[3] == GHOST {
                        self.conflict(t, p)
                    } else {
                        (0..4).all(|i| self.orient(face(&x.v, i), p) != Sign::Positive)
                    })
            })
            .unwrap_or(t)
    }

    fn insert_index(&mut self, i: u32) -> Result<Inserted, DelaunayError> {
        let fail = DelaunayError::Cavity { point: i };
        let p = *self.pt(i);
        let mut t0 = self.locate(&p);
        if !self.is_ghost(t0) {
            for &v in &self.tets[t0 as usize].v {
                if implicit::coincident(self.pt(v), &p) {
                    return Ok(Inserted::Duplicate(v));
                }
            }
        }
        if !self.conflict(t0, &p) {
            // Defensive: the located tet always conflicts.
            t0 = (0..self.tets.len() as u32)
                .find(|&t| self.tets[t as usize].alive && self.conflict(t, &p))
                .ok_or(fail.clone())?;
        }
        let epoch = self.next_epoch();
        let mut cavity = vec![t0];
        self.stamp[t0 as usize] = epoch;
        let mut k = 0;
        while k < cavity.len() {
            let t = cavity[k];
            k += 1;
            for j in 0..4 {
                let u = self.tets[t as usize].n[j];
                if self.stamp[u as usize] != epoch && self.conflict(u, &p) {
                    self.stamp[u as usize] = epoch;
                    cavity.push(u);
                }
            }
        }
        // Grow across every boundary face p does not strictly see: then each
        // cone tet is positive (or a valid ghost).
        loop {
            let mut grow = Vec::new();
            for &t in &cavity {
                for j in 0..4 {
                    let u = self.tets[t as usize].n[j];
                    if self.stamp[u as usize] != epoch && !self.sees(t, j, &p) {
                        self.stamp[u as usize] = epoch;
                        grow.push(u);
                    }
                }
            }
            if grow.is_empty() {
                break;
            }
            cavity.extend(grow);
            if cavity.len() >= self.alive {
                return Err(fail);
            }
        }
        self.fill_cavity(i, &cavity)
    }

    /// A fresh cavity stamp (never 0, which marks a tet taken back out).
    fn next_epoch(&mut self) -> u32 {
        self.epoch = self.epoch.wrapping_add(1).max(1);
        self.epoch
    }

    /// Replace the cavity (its tets stamped with the current epoch) by a cone
    /// from vertex `i` over its boundary.
    fn fill_cavity(&mut self, i: u32, cavity: &[u32]) -> Result<Inserted, DelaunayError> {
        let fail = DelaunayError::Cavity { point: i };
        let epoch = self.epoch;
        let boundary: Vec<(u32, usize, u32)> = cavity
            .iter()
            .flat_map(|&t| (0..4).map(move |j| (t, j)))
            .map(|(t, j)| (t, j, self.tets[t as usize].n[j]))
            .filter(|&(_, _, u)| self.stamp[u as usize] != epoch)
            .collect();
        // No vertex may be left inside the cavity.
        let mut on_boundary: Vec<u32> = boundary
            .iter()
            .flat_map(|&(t, j, _)| face(&self.tets[t as usize].v, j))
            .collect();
        on_boundary.sort_unstable();
        on_boundary.dedup();
        let lost = cavity
            .iter()
            .flat_map(|&t| self.tets[t as usize].v)
            .any(|x| on_boundary.binary_search(&x).is_err());
        if lost || boundary.is_empty() {
            return Err(fail);
        }
        // Cone from p over the boundary.
        let mut created = Vec::with_capacity(boundary.len());
        for &(t, j, u) in &boundary {
            let f = face(&self.tets[t as usize].v, j);
            let mut v = [f[1], f[0], f[2], i];
            if let Some(g) = v.iter().position(|&x| x == GHOST) {
                if g != 3 {
                    v.swap(g, 3);
                    v.swap(0, 1);
                }
            }
            let id = self.alloc(v);
            let at_p = v.iter().position(|&x| x == i).unwrap_or(3);
            self.tets[id as usize].n[at_p] = u;
            if let Some(back) = (0..4).find(|&m| self.tets[u as usize].n[m] == t) {
                self.tets[u as usize].n[back] = id;
            }
            created.push(id);
        }
        for &t in cavity {
            self.kill(t);
        }
        self.link_around(&created, i).map_err(|()| fail)?;
        self.last = created[0];
        Ok(Inserted::New(i))
    }

    /// Insert a point without removing any edge `protected` accepts: the
    /// cavity starts as the Delaunay conflict region and shrinks until every
    /// boundary face is strictly visible from the point, no protected edge
    /// and no vertex would be lost, keeping every tet whose closure contains
    /// the point (after Si's constrained Bowyer-Watson in TetGen). The result
    /// need not be Delaunay.
    ///
    /// # Errors
    /// [`DelaunayError::Outside`] when the point is outside the convex hull
    /// (segment points never are), [`DelaunayError::Cavity`] when it lies on
    /// a protected edge or the cavity cannot be filled; nothing changes then.
    pub fn insert_protected(
        &mut self,
        p: Point,
        protected: &dyn Fn(u32, u32) -> bool,
    ) -> Result<Inserted, DelaunayError> {
        let i = self.points.len() as u32;
        self.points.push(p);
        self.vertex_tet.push(NONE);
        let result = self.insert_index_protected(i, protected);
        if !matches!(result, Ok(Inserted::New(_))) {
            self.points.pop();
            self.vertex_tet.pop();
        }
        result
    }

    fn insert_index_protected(
        &mut self,
        i: u32,
        protected: &dyn Fn(u32, u32) -> bool,
    ) -> Result<Inserted, DelaunayError> {
        let fail = DelaunayError::Cavity { point: i };
        let p = *self.pt(i);
        let s = crate::prof::scope("cdt insert: locate");
        let t0 = self.locate(&p);
        drop(s);
        if self.is_ghost(t0) {
            return Err(DelaunayError::Outside { point: i });
        }
        for &v in &self.tets[t0 as usize].v {
            if implicit::coincident(self.pt(v), &p) {
                return Ok(Inserted::Duplicate(v));
            }
        }
        // Every tet whose closure contains p: connected across the faces p
        // lies on.
        let s = crate::prof::scope("cdt insert: must and conflicts");
        let mut must = vec![t0];
        let mut k = 0;
        while k < must.len() {
            let t = must[k];
            k += 1;
            let v = self.tets[t as usize].v;
            for j in 0..4 {
                let u = self.tets[t as usize].n[j];
                if must.contains(&u) {
                    continue;
                }
                let f = face(&v, j);
                let on = if f.contains(&GHOST) {
                    let r: Vec<u32> = f.iter().copied().filter(|&x| x != GHOST).collect();
                    implicit::collinear(self.pt(r[0]), self.pt(r[1]), &p)
                } else {
                    self.orient(f, &p) == Sign::Zero
                };
                if on {
                    must.push(u);
                }
            }
        }
        // Start from the conflict region around them.
        let epoch = self.next_epoch();
        for &t in &must {
            self.stamp[t as usize] = epoch;
        }
        let mut cavity = must.clone();
        let mut k = 0;
        while k < cavity.len() {
            let t = cavity[k];
            k += 1;
            for j in 0..4 {
                let u = self.tets[t as usize].n[j];
                if self.stamp[u as usize] != epoch && self.conflict(u, &p) {
                    self.stamp[u as usize] = epoch;
                    cavity.push(u);
                }
            }
        }
        drop(s);
        // Shrink, one tet at a time, until the cavity is valid.
        let _s = crate::prof::scope("cdt insert: shrink cavity");
        loop {
            // Keep what is connected to the tets containing p.
            let mut reach = must.clone();
            let mut seen = must.clone();
            let mut k = 0;
            while k < reach.len() {
                let t = reach[k];
                k += 1;
                for j in 0..4 {
                    let u = self.tets[t as usize].n[j];
                    if self.stamp[u as usize] == epoch && !seen.contains(&u) {
                        seen.push(u);
                        reach.push(u);
                    }
                }
            }
            for &t in &cavity {
                if !seen.contains(&t) {
                    self.stamp[t as usize] = 0;
                }
            }
            cavity = reach;
            match self.cavity_violation(&cavity, &must, &p, protected) {
                Violation::None => break,
                Violation::Remove(t) => {
                    self.stamp[t as usize] = 0;
                    cavity.retain(|&x| x != t);
                }
                Violation::Stuck => {
                    for &t in &cavity {
                        self.stamp[t as usize] = 0;
                    }
                    return Err(fail);
                }
            }
        }
        drop(_s);
        let _s = crate::prof::scope("cdt insert: fill");
        self.fill_cavity(i, &cavity)
    }

    /// The first thing wrong with a constrained cavity, and a tet to take out
    /// for it (never one of `must`).
    fn cavity_violation(
        &self,
        cavity: &[u32],
        must: &[u32],
        p: &Point,
        protected: &dyn Fn(u32, u32) -> bool,
    ) -> Violation {
        let epoch = self.epoch;
        let inside = |t: u32| self.stamp[t as usize] == epoch;
        // A boundary face p does not strictly see.
        for &t in cavity {
            if must.contains(&t) {
                continue;
            }
            for j in 0..4 {
                let u = self.tets[t as usize].n[j];
                if !inside(u) && !self.sees(t, j, p) {
                    return Violation::Remove(t);
                }
            }
        }
        // A protected edge with every tet around it in the cavity.
        for &t in cavity {
            let v = self.tets[t as usize].v;
            for (a, b) in EDGES {
                let (x, y) = (v[a], v[b]);
                if x == GHOST || y == GHOST || !protected(x, y) {
                    continue;
                }
                let Some(ring) = self.edge_ring(x, y) else {
                    continue;
                };
                if ring.iter().all(|&r| inside(r)) {
                    return match ring.iter().find(|r| !must.contains(r)) {
                        Some(&r) => Violation::Remove(r),
                        None => Violation::Stuck,
                    };
                }
            }
        }
        // A vertex inside the cavity, on no boundary face.
        let mut on_boundary: Vec<u32> = cavity
            .iter()
            .flat_map(|&t| (0..4).map(move |j| (t, j)))
            .filter(|&(t, j)| !inside(self.tets[t as usize].n[j]))
            .flat_map(|(t, j)| face(&self.tets[t as usize].v, j))
            .collect();
        on_boundary.sort_unstable();
        on_boundary.dedup();
        for &t in cavity {
            for x in self.tets[t as usize].v {
                if on_boundary.binary_search(&x).is_err() {
                    return match cavity
                        .iter()
                        .find(|&&r| !must.contains(&r) && self.tets[r as usize].v.contains(&x))
                    {
                        Some(&r) => Violation::Remove(r),
                        None => Violation::Stuck,
                    };
                }
            }
        }
        Violation::None
    }

    fn orient4(&self, v: [u32; 4]) -> Sign {
        let [a, b, c, d] = v.map(|x| self.pt(x));
        implicit::orient3d(a, b, c, d)
    }

    /// The tets around edge `xy` in cyclic order (ghosts included), or `None`
    /// if `xy` is not an edge.
    #[must_use]
    pub fn edge_ring(&self, x: u32, y: u32) -> Option<Vec<u32>> {
        let start = *self
            .star(x)
            .iter()
            .find(|&&t| self.tets[t as usize].v.contains(&y))?;
        let mut ring = vec![start];
        let (mut prev, mut t) = (NONE, start);
        loop {
            let tet = self.tets[t as usize];
            // The two faces holding x and y are opposite the other vertices.
            let next = (0..4)
                .filter(|&k| tet.v[k] != x && tet.v[k] != y)
                .map(|k| tet.n[k])
                .find(|&u| u != prev)?;
            if next == start {
                return Some(ring);
            }
            if ring.len() > self.tets.len() {
                return None;
            }
            ring.push(next);
            (prev, t) = (t, next);
        }
    }

    /// Replace the real tets `old` by `new` if that is a valid local
    /// retriangulation: every new tet positive, and the faces of `new` that
    /// are not shared between new tets exactly the boundary of `old`, with the
    /// same orientation. Then both tile the same region. Returns whether it
    /// was done; nothing changes otherwise.
    fn replace(&mut self, old: &[u32], new: &[[u32; 4]]) -> bool {
        // Ghosts (the vertex at infinity last) only take part in hull flips,
        // whose callers check the hull stays the same.
        if new.iter().any(|v| v[..3].contains(&GHOST)) {
            return false;
        }
        if new
            .iter()
            .any(|&v| v[3] != GHOST && self.orient4(v) != Sign::Positive)
        {
            return false;
        }
        // The boundary of `old`, oriented outward, with the tet beyond.
        let mut boundary: Vec<([u32; 3], u32)> = Vec::new();
        for &t in old {
            let tet = self.tets[t as usize];
            for j in 0..4 {
                if !old.contains(&tet.n[j]) {
                    boundary.push((rotated(face(&tet.v, j)), tet.n[j]));
                }
            }
        }
        boundary.sort_unstable();
        let faces: Vec<[u32; 3]> = new
            .iter()
            .flat_map(|v| (0..4).map(move |j| rotated(face(v, j))))
            .collect();
        let mut matched = 0;
        for f in &faces {
            let reverse = rotated([f[0], f[2], f[1]]);
            if faces.contains(&reverse) {
                continue;
            }
            if boundary.binary_search_by(|(b, _)| b.cmp(f)).is_err() {
                return false;
            }
            matched += 1;
        }
        if matched != boundary.len() {
            return false;
        }
        for &t in old {
            self.kill(t);
        }
        let ids: Vec<u32> = new.iter().map(|&v| self.alloc(v)).collect();
        let mut internal: Vec<([u32; 3], u32, usize)> = Vec::new();
        for &id in &ids {
            let v = self.tets[id as usize].v;
            for j in 0..4 {
                let f = rotated(face(&v, j));
                match boundary.binary_search_by(|(b, _)| b.cmp(&f)) {
                    Ok(at) => {
                        let u = boundary[at].1;
                        self.tets[id as usize].n[j] = u;
                        let key = sorted3(f);
                        if let Some(m) =
                            (0..4).find(|&m| sorted3(face(&self.tets[u as usize].v, m)) == key)
                        {
                            self.tets[u as usize].n[m] = id;
                        }
                    }
                    Err(_) => internal.push((sorted3(f), id, j)),
                }
            }
        }
        internal.sort_unstable();
        for pair in internal.chunks(2) {
            if let [(_, ta, ka), (_, tb, kb)] = *pair {
                self.tets[ta as usize].n[ka] = tb;
                self.tets[tb as usize].n[kb] = ta;
            }
        }
        self.last = ids[0];
        true
    }

    /// The 2-3 flip of face `j` of tet `t`: the two tets sharing it become
    /// three around the edge joining their apices, if that is valid.
    pub fn flip23(&mut self, t: u32, j: usize) -> bool {
        let tet = self.tets[t as usize];
        let u = tet.n[j];
        if !tet.alive || self.is_ghost(t) || self.is_ghost(u) {
            return false;
        }
        let f = face(&tet.v, j);
        let a = tet.v[j];
        let Some(&w) = self.tets[u as usize].v.iter().find(|x| !f.contains(x)) else {
            return false;
        };
        let new = [[f[0], f[1], a, w], [f[1], f[2], a, w], [f[2], f[0], a, w]];
        self.replace(&[t, u], &new)
    }

    /// The 3-2 flip of edge `xy` with exactly three tets around it, if valid.
    fn flip32(&mut self, ring: &[u32], x: u32, y: u32) -> bool {
        let mut p: Vec<u32> = ring
            .iter()
            .flat_map(|&t| self.tets[t as usize].v)
            .filter(|&v| v != x && v != y)
            .collect();
        p.sort_unstable();
        p.dedup();
        let [p0, p1, p2] = p[..] else {
            return false;
        };
        self.replace(ring, &[[p0, p1, p2, x], [p0, p2, p1, y]])
            || self.replace(ring, &[[p0, p2, p1, x], [p0, p1, p2, y]])
    }

    /// The 2-2 flip of hull edge `xy` whose ring is `[ghost, real, real,
    /// ghost]`: the hull faces `xya` and `xyc` are coplanar, and the edge
    /// becomes `ac`, if that is valid. The hull does not change.
    fn flip22_hull(&mut self, ring: &[u32], x: u32, y: u32) -> bool {
        let apex = |t: u32, not: &[u32]| {
            self.tets[t as usize]
                .v
                .into_iter()
                .find(|v| !not.contains(v) && *v != GHOST)
        };
        let (Some(a), Some(c)) = (apex(ring[0], &[x, y]), apex(ring[3], &[x, y])) else {
            return false;
        };
        let Some(w) = apex(ring[1], &[x, y, a]) else {
            return false;
        };
        if self.orient4([x, y, a, c]) != Sign::Zero {
            return false;
        }
        let mut new = Vec::with_capacity(4);
        for e in [x, y] {
            let mut v = [a, c, e, w];
            if self.orient4(v) == Sign::Negative {
                v.swap(0, 1);
            }
            // The hull face, outward from the real tet, is opposite w.
            let f = face(&v, 3);
            new.push(v);
            new.push([f[0], f[1], f[2], GHOST]);
        }
        self.replace(ring, &new)
    }

    /// Remove edge `xy` by flips: 2-3 flips of faces around it shrink its
    /// ring to three tets, then a 3-2 flip removes it. Greedy, so it can fail
    /// part way (the tetrahedralization stays valid). Never removes an edge
    /// `protected` accepts.
    fn remove_edge(&mut self, x: u32, y: u32, protected: &dyn Fn(u32, u32) -> bool) -> bool {
        if protected(x, y) {
            return false;
        }
        loop {
            let Some(mut ring) = self.edge_ring(x, y) else {
                return true;
            };
            if ring.iter().any(|&t| self.is_ghost(t)) {
                // A hull edge: its ring is two adjacent ghosts and real tets.
                // Start at the ghost a real tet follows, shrink the real part
                // to two, then 2-2 flip.
                let n = ring.len();
                let Some(g) =
                    (0..n).find(|&k| self.is_ghost(ring[k]) && !self.is_ghost(ring[(k + 1) % n]))
                else {
                    return false;
                };
                ring.rotate_left(g);
                if n < 4 || !self.is_ghost(ring[n - 1]) {
                    return false;
                }
                if n == 4 {
                    return self.flip22_hull(&ring, x, y);
                }
                let flipped = (1..n - 2).any(|i| {
                    let (t, u) = (ring[i], ring[i + 1]);
                    (0..4)
                        .find(|&j| self.tets[t as usize].n[j] == u)
                        .is_some_and(|j| self.flip23(t, j))
                });
                if !flipped {
                    return false;
                }
                continue;
            }
            if ring.len() <= 3 {
                return ring.len() == 3 && self.flip32(&ring, x, y);
            }
            let n = ring.len();
            let flipped = (0..n).any(|i| {
                let (t, u) = (ring[i], ring[(i + 1) % n]);
                (0..4)
                    .find(|&j| self.tets[t as usize].n[j] == u)
                    .is_some_and(|j| self.flip23(t, j))
            });
            if !flipped {
                return false;
            }
        }
    }

    /// What the segment from vertex `a` to vertex `b` first meets after
    /// leaving `a`.
    fn crossing(&self, a: u32, b: u32) -> Crossing {
        let pb = *self.pt(b);
        for t in self.star(a) {
            if self.is_ghost(t) {
                continue;
            }
            let v = self.tets[t as usize].v;
            let Some(k) = v.iter().position(|&x| x == a) else {
                continue;
            };
            let mut zeros = Vec::new();
            let mut outside = false;
            for j in (0..4).filter(|&j| j != k) {
                match self.orient(face(&v, j), &pb) {
                    Sign::Positive => {
                        outside = true;
                        break;
                    }
                    Sign::Zero => zeros.push(j),
                    Sign::Negative => {}
                }
            }
            if outside {
                continue;
            }
            return match zeros[..] {
                [] if self.orient(face(&v, k), &pb) == Sign::Positive => Crossing::Face { t, j: k },
                [] => Crossing::None,
                [z] => {
                    let e: Vec<u32> = face(&v, z).into_iter().filter(|&x| x != a).collect();
                    Crossing::Edge { x: e[0], y: e[1] }
                }
                _ => Crossing::None,
            };
        }
        Crossing::None
    }

    /// One flip, or edge removal, toward making `ab` an edge, from `a`'s end.
    fn flip_toward(&mut self, a: u32, b: u32, protected: &dyn Fn(u32, u32) -> bool) -> bool {
        match self.crossing(a, b) {
            Crossing::Face { t, j } => {
                if self.flip23(t, j) {
                    return true;
                }
                // Not convex across the face: remove the edge in the way.
                let tet = self.tets[t as usize];
                let f = face(&tet.v, j);
                let Some(&w) = self.tets[tet.n[j] as usize]
                    .v
                    .iter()
                    .find(|x| !f.contains(x))
                else {
                    return false;
                };
                (0..3).any(|m| {
                    let (e0, e1) = (f[m], f[(m + 1) % 3]);
                    self.orient4([e0, e1, a, w]) != Sign::Positive
                        && self.remove_edge(e0, e1, protected)
                })
            }
            Crossing::Edge { x, y } => self.remove_edge(x, y, protected),
            Crossing::None => false,
        }
    }

    /// Try to make `ab` an edge by flips alone, at most `max_steps` steps,
    /// never removing an edge `protected` accepts. The result need not be
    /// Delaunay. Returns whether `ab` is an edge.
    pub fn recover_edge(
        &mut self,
        a: u32,
        b: u32,
        protected: &dyn Fn(u32, u32) -> bool,
        max_steps: usize,
    ) -> bool {
        for _ in 0..max_steps {
            if self.has_edge(a, b) {
                return true;
            }
            if !(self.flip_toward(a, b, protected) || self.flip_toward(b, a, protected)) {
                break;
            }
        }
        // A failed attempt keeps its flips: undoing them was measured to
        // recover fewer segments overall.
        self.has_edge(a, b)
    }

    /// Structural and geometric validity: every real tet positive (exactly),
    /// ghosts with the vertex at infinity last, neighbours mutual and sharing
    /// the face.
    ///
    /// # Errors
    /// What is wrong.
    pub fn check(&self) -> Result<(), String> {
        for (ti, t) in self.tets.iter().enumerate() {
            if !t.alive {
                continue;
            }
            if t.v[..3].contains(&GHOST) {
                return Err(format!("tet {ti}: ghost vertex not last"));
            }
            if t.v[3] != GHOST {
                let [a, b, c, d] = t.v.map(|v| self.pt(v));
                if implicit::orient3d(a, b, c, d) != Sign::Positive {
                    return Err(format!("tet {ti} {:?} is not positive", t.v));
                }
            }
            for k in 0..4 {
                let u = t.n[k];
                let Some(nb) = self.tets.get(u as usize).filter(|n| n.alive) else {
                    return Err(format!("tet {ti}: neighbour {k} missing"));
                };
                let Some(back) = (0..4).find(|&j| nb.n[j] == ti as u32) else {
                    return Err(format!("tet {ti}: neighbour {u} does not point back"));
                };
                if sorted3(face(&t.v, k)) != sorted3(face(&nb.v, back)) {
                    return Err(format!("tets {ti} and {u} disagree on their shared face"));
                }
            }
        }
        Ok(())
    }

    /// The Delaunay property, by brute force (for tests): no vertex strictly
    /// inside any real tet's circumsphere.
    ///
    /// # Errors
    /// The first violation.
    pub fn check_delaunay(&self) -> Result<(), String> {
        let mut vertices: Vec<u32> = self
            .tets
            .iter()
            .filter(|t| t.alive)
            .flat_map(|t| t.v)
            .filter(|&v| v != GHOST)
            .collect();
        vertices.sort_unstable();
        vertices.dedup();
        for t in self.tets.iter().filter(|t| t.alive && t.v[3] != GHOST) {
            let [a, b, c, d] = t.v.map(|v| self.pt(v));
            for &v in &vertices {
                if !t.v.contains(&v) && implicit::insphere(a, b, c, d, self.pt(v)) == Sign::Positive
                {
                    return Err(format!("vertex {v} inside the sphere of {:?}", t.v));
                }
            }
        }
        Ok(())
    }
}
