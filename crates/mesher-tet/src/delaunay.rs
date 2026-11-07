//! Incremental 3D Delaunay tetrahedralization (Bowyer–Watson) with exact
//! predicates.
//!
//! Data structure: a flat array of tetrahedra, each with four vertices and
//! four neighbours (`n[i]` is across the face opposite `v[i]`), and a free
//! list for deleted slots. The convex hull is closed off by *ghost* tets that
//! share a vertex at infinity ([`GHOST`], always stored last): ghost
//! `[a, b, c, ∞]` sits outside hull face `abc`, and `orient3d(a, b, c, q)` is
//! positive for points `q` beyond that face. Real tets are positively
//! oriented.
//!
//! Inserting a point: walk from the last new tet towards it (visibility walk,
//! deterministic pseudo-random face order); grow the cavity of tets whose
//! circumsphere strictly contains it (for a ghost: the point is beyond its hull
//! face, or on the face's plane inside its circumcircle); re-triangulate the
//! cavity as a cone from the point. If a cavity boundary face cannot see the
//! point (possible only with cospherical or coplanar points) the cavity is
//! grown across it, keeping the result a valid triangulation.
//!
//! Points are inserted in biased randomized insertion order (BRIO): rounds of
//! doubling size, each sorted along a Morton curve, so walks stay short.

use mesher_geom::predicates::{collinear, in_circumcircle_coplanar, insphere, orient3d};
use mesher_geom::{Point3, Sign};

/// The vertex at infinity.
pub const GHOST: u32 = u32::MAX;
const NONE: u32 = u32::MAX;

/// Local face opposite each vertex, ordered outward for a positive tet.
const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

#[derive(Debug, Clone, Copy)]
struct Tet {
    v: [u32; 4],
    n: [u32; 4],
    alive: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DelaunayError {
    NonFinite {
        point: u32,
    },
    /// Fewer than four points, or all of them coplanar.
    Flat,
    /// A cavity could not be re-triangulated (a bug or a pathological input).
    Cavity {
        point: u32,
    },
}

impl std::fmt::Display for DelaunayError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DelaunayError::NonFinite { point } => write!(f, "point {point} is not finite"),
            DelaunayError::Flat => write!(f, "the points do not span a volume"),
            DelaunayError::Cavity { point } => {
                write!(f, "could not re-triangulate the cavity of point {point}")
            }
        }
    }
}

/// Result of inserting a point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Inserted {
    New,
    /// Exactly equal to an existing vertex, which is returned.
    Duplicate(u32),
}

/// A Delaunay tetrahedralization that accepts further points.
#[derive(Debug, Clone)]
pub struct Delaunay {
    points: Vec<Point3>,
    tets: Vec<Tet>,
    free: Vec<u32>,
    stamp: Vec<u32>,
    epoch: u32,
    alive: usize,
    last: u32,
    rng: SplitMix,
    duplicates: Vec<(u32, u32)>,
}

impl Delaunay {
    /// Tetrahedralize `points` (inserted in BRIO order from `seed`).
    pub fn new(points: Vec<Point3>, seed: u64) -> Result<Delaunay, DelaunayError> {
        if let Some(i) = points.iter().position(|p| !p.iter().all(|c| c.is_finite())) {
            return Err(DelaunayError::NonFinite { point: i as u32 });
        }
        let order = brio_order(&points, seed);
        let mut d = Delaunay {
            points,
            tets: Vec::new(),
            free: Vec::new(),
            stamp: Vec::new(),
            epoch: 0,
            alive: 0,
            last: 0,
            rng: SplitMix(seed ^ 0x9e37_79b9_7f4a_7c15),
            duplicates: Vec::new(),
        };
        let first = d.initial_tet(&order)?;
        for &i in &order {
            if first.contains(&i) {
                continue;
            }
            if let Inserted::Duplicate(v) = d.insert_index(i)? {
                d.duplicates.push((i, v));
            }
        }
        Ok(d)
    }

    #[must_use]
    pub fn points(&self) -> &[Point3] {
        &self.points
    }

    /// `(point, existing vertex)` for every point equal to an earlier one.
    #[must_use]
    pub fn duplicates(&self) -> &[(u32, u32)] {
        &self.duplicates
    }

    /// Add a point and insert it; returns its index and the outcome.
    pub fn insert(&mut self, p: Point3) -> Result<(u32, Inserted), DelaunayError> {
        let i = self.points.len() as u32;
        if !p.iter().all(|c| c.is_finite()) {
            return Err(DelaunayError::NonFinite { point: i });
        }
        self.points.push(p);
        let outcome = self.insert_index(i)?;
        if let Inserted::Duplicate(v) = outcome {
            self.duplicates.push((i, v));
        }
        Ok((i, outcome))
    }

    /// The real tetrahedra, positively oriented.
    #[must_use]
    pub fn tets(&self) -> Vec<[u32; 4]> {
        self.tets
            .iter()
            .filter(|t| t.alive && t.v[3] != GHOST)
            .map(|t| t.v)
            .collect()
    }

    /// Convex hull faces, oriented outward.
    #[must_use]
    pub fn hull(&self) -> Vec<[u32; 3]> {
        self.tets
            .iter()
            .filter(|t| t.alive && t.v[3] == GHOST)
            .map(|t| [t.v[0], t.v[1], t.v[2]])
            .collect()
    }

    /// Check the structure's invariants: neighbours alive and reciprocal,
    /// sharing exactly the face between them; real tets positive; ghost
    /// vertex last.
    pub fn check(&self) -> Result<(), String> {
        for (ti, t) in self.tets.iter().enumerate() {
            if !t.alive {
                continue;
            }
            if t.v[..3].contains(&GHOST) {
                return Err(format!("tet {ti}: ghost vertex not last"));
            }
            if t.v[3] != GHOST {
                let [a, b, c, d] = t.v.map(|v| self.pos(v));
                if orient3d(a, b, c, d) != Sign::Positive {
                    return Err(format!("tet {ti} is not positive"));
                }
            }
            for k in 0..4 {
                let u = t.n[k];
                let Some(nb) = self.tets.get(u as usize).filter(|n| n.alive) else {
                    return Err(format!("tet {ti}: neighbour {k} missing"));
                };
                let mut mine = face(&t.v, k);
                mine.sort_unstable();
                let Some(back) = (0..4).find(|&j| nb.n[j] == ti as u32) else {
                    return Err(format!("tet {ti}: neighbour {u} does not point back"));
                };
                let mut theirs = face(&nb.v, back);
                theirs.sort_unstable();
                if mine != theirs {
                    return Err(format!("tets {ti} and {u} disagree on their shared face"));
                }
            }
        }
        Ok(())
    }

    fn pos(&self, v: u32) -> Point3 {
        self.points[v as usize]
    }

    fn alloc(&mut self, v: [u32; 4]) -> u32 {
        let t = Tet {
            v,
            n: [NONE; 4],
            alive: true,
        };
        self.alive += 1;
        if let Some(id) = self.free.pop() {
            self.tets[id as usize] = t;
            id
        } else {
            self.tets.push(t);
            self.stamp.push(0);
            (self.tets.len() - 1) as u32
        }
    }

    fn kill(&mut self, t: u32) {
        self.tets[t as usize].alive = false;
        self.alive -= 1;
        self.free.push(t);
    }

    /// The first four affinely independent points (in insertion order) form
    /// the first tet, closed off by four ghosts.
    fn initial_tet(&mut self, order: &[u32]) -> Result<[u32; 4], DelaunayError> {
        let pos = |i: u32| self.points[i as usize];
        let a = *order.first().ok_or(DelaunayError::Flat)?;
        let b = *order
            .iter()
            .find(|&&i| pos(i) != pos(a))
            .ok_or(DelaunayError::Flat)?;
        let c = *order
            .iter()
            .find(|&&i| !collinear(pos(a), pos(b), pos(i)))
            .ok_or(DelaunayError::Flat)?;
        let (d, s) = order
            .iter()
            .map(|&i| (i, orient3d(pos(a), pos(b), pos(c), pos(i))))
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
                    let mut f = face(&v, k);
                    f.sort_unstable();
                    keys.push((f, t, k));
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

    fn is_ghost(&self, t: u32) -> bool {
        self.tets[t as usize].v[3] == GHOST
    }

    fn conflict(&self, t: u32, p: Point3) -> bool {
        let v = self.tets[t as usize].v;
        let [a, b, c] = [v[0], v[1], v[2]].map(|x| self.pos(x));
        if v[3] == GHOST {
            match orient3d(a, b, c, p) {
                Sign::Positive => true,
                Sign::Negative => false,
                Sign::Zero => in_circumcircle_coplanar(a, b, c, p),
            }
        } else {
            insphere(a, b, c, self.pos(v[3]), p) == Sign::Positive
        }
    }

    /// A tet whose closure contains `p`, or a ghost whose hull face `p` is beyond.
    fn locate(&mut self, p: Point3) -> u32 {
        let mut t = self.last;
        if !self.tets.get(t as usize).is_some_and(|x| x.alive) {
            t = self.tets.iter().position(|x| x.alive).unwrap_or(0) as u32;
        }
        // Start inside the hull: a ghost's real neighbour is across its hull face.
        if self.is_ghost(t) {
            t = self.tets[t as usize].n[3];
        }
        let limit = 4 * self.tets.len() + 64;
        for _ in 0..limit {
            // A ghost is reached only by crossing a hull face that p is beyond.
            if self.is_ghost(t) {
                return t;
            }
            let v = self.tets[t as usize].v;
            let start = (self.rng.next() % 4) as usize;
            let mut next = None;
            for k in 0..4 {
                let i = (start + k) % 4;
                let f = face(&v, i).map(|x| self.pos(x));
                if orient3d(f[0], f[1], f[2], p) == Sign::Positive {
                    next = Some(self.tets[t as usize].n[i]);
                    break;
                }
            }
            match next {
                Some(u) => t = u,
                None => return t,
            }
        }
        // The walk should always terminate; fall back to a scan if it did not.
        let inside = (0..self.tets.len() as u32).find(|&t| {
            let x = self.tets[t as usize];
            x.alive
                && x.v[3] != GHOST
                && (0..4).all(|i| {
                    let f = face(&x.v, i).map(|v| self.pos(v));
                    orient3d(f[0], f[1], f[2], p) != Sign::Positive
                })
        });
        inside.unwrap_or_else(|| {
            (0..self.tets.len() as u32)
                .find(|&t| self.tets[t as usize].alive && self.is_ghost(t) && self.conflict(t, p))
                .unwrap_or(t)
        })
    }

    fn insert_index(&mut self, i: u32) -> Result<Inserted, DelaunayError> {
        let p = self.pos(i);
        let mut t0 = self.locate(p);
        if self.is_ghost(t0) && !self.conflict(t0, p) {
            // Defensive: never seed a cavity from a ghost p does not conflict with.
            match (0..self.tets.len() as u32)
                .find(|&t| self.tets[t as usize].alive && self.conflict(t, p))
            {
                Some(t) => t0 = t,
                None => return Err(DelaunayError::Cavity { point: i }),
            }
        }
        if !self.is_ghost(t0) {
            for &v in &self.tets[t0 as usize].v {
                if self.pos(v) == p {
                    return Ok(Inserted::Duplicate(v));
                }
            }
        }
        let fail = DelaunayError::Cavity { point: i };

        // Cavity: tets whose circumsphere strictly contains p.
        self.epoch = self.epoch.wrapping_add(1);
        let epoch = self.epoch;
        let mut cavity = vec![t0];
        self.stamp[t0 as usize] = epoch;
        let mut k = 0;
        while k < cavity.len() {
            let t = cavity[k];
            k += 1;
            for j in 0..4 {
                let u = self.tets[t as usize].n[j];
                if self.stamp[u as usize] != epoch && self.conflict(u, p) {
                    self.stamp[u as usize] = epoch;
                    cavity.push(u);
                }
            }
        }

        // Boundary faces must see p; grow the cavity across any that do not.
        let boundary = loop {
            let mut boundary: Vec<(u32, usize, u32)> = Vec::new();
            for &t in &cavity {
                for j in 0..4 {
                    let u = self.tets[t as usize].n[j];
                    if self.stamp[u as usize] != epoch {
                        boundary.push((t, j, u));
                    }
                }
            }
            let blocked = boundary.iter().find(|&&(t, j, _)| {
                let f = face(&self.tets[t as usize].v, j);
                if f.contains(&GHOST) {
                    let r: Vec<Point3> = f
                        .iter()
                        .filter(|&&x| x != GHOST)
                        .map(|&x| self.pos(x))
                        .collect();
                    collinear(r[0], r[1], p)
                } else {
                    let [a, b, c] = f.map(|x| self.pos(x));
                    orient3d(a, b, c, p) != Sign::Negative
                }
            });
            match blocked {
                Some(&(_, _, u)) => {
                    if cavity.len() + 1 >= self.alive {
                        return Err(fail);
                    }
                    self.stamp[u as usize] = epoch;
                    cavity.push(u);
                }
                None => break boundary,
            }
        };
        if boundary.is_empty() {
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
        for &t in &cavity {
            self.kill(t);
        }
        self.link_around(&created, i).map_err(|()| fail.clone())?;
        self.last = created[0];
        Ok(Inserted::New)
    }
}

/// The face opposite local vertex `k`, outward for a positive tet.
fn face(v: &[u32; 4], k: usize) -> [u32; 3] {
    FACES[k].map(|i| v[i])
}

/// A tetrahedralization of points, with the inputs that were not used.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tetrahedralization {
    pub tets: Vec<[u32; 4]>,
    /// `(point, equal earlier point)` for exact duplicates, which are not vertices.
    pub duplicates: Vec<(u32, u32)>,
}

/// The Delaunay tetrahedralization of `points`, deterministic for a seed.
pub fn delaunay(points: &[Point3]) -> Result<Tetrahedralization, DelaunayError> {
    let d = Delaunay::new(points.to_vec(), 0x5eed)?;
    Ok(Tetrahedralization {
        tets: d.tets(),
        duplicates: d.duplicates().to_vec(),
    })
}

/// SplitMix64: a small, seedable generator for deterministic shuffles.
#[derive(Debug, Clone)]
struct SplitMix(u64);

impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}

/// Biased randomized insertion order: shuffle, split into rounds of doubling
/// size, and sort each round along a Morton curve.
fn brio_order(points: &[Point3], seed: u64) -> Vec<u32> {
    let n = points.len();
    let mut order: Vec<u32> = (0..n as u32).collect();
    let mut rng = SplitMix(seed);
    for i in (1..n).rev() {
        let j = (rng.next() % (i as u64 + 1)) as usize;
        order.swap(i, j);
    }
    let lo = [0, 1, 2].map(|k| points.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min));
    let hi = [0, 1, 2].map(|k| {
        points
            .iter()
            .map(|p| p[k])
            .fold(f64::NEG_INFINITY, f64::max)
    });
    let morton = |p: Point3| -> u64 {
        let q = [0, 1, 2].map(|k| {
            let span = hi[k] - lo[k];
            let t = if span > 0.0 {
                (p[k] - lo[k]) / span
            } else {
                0.0
            };
            ((t * f64::from((1u32 << 21) - 1)) as u64).min((1 << 21) - 1)
        });
        let spread = |mut x: u64| {
            x &= 0x1f_ffff;
            x = (x | (x << 32)) & 0x1f_0000_0000_ffff;
            x = (x | (x << 16)) & 0x1f_0000_ff00_00ff;
            x = (x | (x << 8)) & 0x100f_00f0_0f00_f00f;
            x = (x | (x << 4)) & 0x10c3_0c30_c30c_30c3;
            x = (x | (x << 2)) & 0x1249_2492_4924_9249;
            x
        };
        spread(q[0]) | (spread(q[1]) << 1) | (spread(q[2]) << 2)
    };
    let mut bounds = vec![n];
    let mut m = n;
    while m > 64 {
        m /= 2;
        bounds.push(m);
    }
    bounds.push(0);
    bounds.reverse();
    for w in bounds.windows(2) {
        order[w[0]..w[1]].sort_by_key(|&i| (morton(points[i as usize]), i));
    }
    order
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_tet_and_its_ghosts() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let d = Delaunay::new(pts, 1).unwrap();
        d.check().unwrap();
        assert_eq!(d.tets().len(), 1);
        assert_eq!(d.hull().len(), 4);
    }

    #[test]
    fn cube_corners_and_centre() {
        let mut pts: Vec<Point3> = (0..8)
            .map(|i| {
                [
                    f64::from(i & 1),
                    f64::from((i >> 1) & 1),
                    f64::from((i >> 2) & 1),
                ]
            })
            .collect();
        pts.push([0.5, 0.5, 0.5]);
        let d = Delaunay::new(pts, 7).unwrap();
        d.check().unwrap();
        // Twelve tets: the centre coned to the 12 triangles of the cube's faces.
        assert_eq!(d.tets().len(), 12);
        assert_eq!(d.hull().len(), 12);
    }

    #[test]
    fn duplicates_flat_and_non_finite() {
        let mut pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        pts.push([1.0, 0.0, 0.0]);
        let d = Delaunay::new(pts, 3).unwrap();
        assert_eq!(d.duplicates(), &[(4, 1)]);
        let flat = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        assert_eq!(Delaunay::new(flat, 1).unwrap_err(), DelaunayError::Flat);
        let bad = vec![[f64::NAN, 0.0, 0.0]];
        assert_eq!(
            Delaunay::new(bad, 1).unwrap_err(),
            DelaunayError::NonFinite { point: 0 }
        );
    }

    #[test]
    fn incremental_insertion() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [0.0, 0.0, 4.0],
        ];
        let mut d = Delaunay::new(pts, 1).unwrap();
        let (i, out) = d.insert([0.5, 0.5, 0.5]).unwrap();
        assert_eq!((i, out), (4, Inserted::New));
        // Outside the hull.
        d.insert([5.0, 5.0, 5.0]).unwrap();
        d.check().unwrap();
        assert_eq!(d.insert([0.0, 0.0, 0.0]).unwrap().1, Inserted::Duplicate(0));
    }
}
