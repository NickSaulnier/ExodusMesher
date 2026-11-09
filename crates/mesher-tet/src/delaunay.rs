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
//!
//! *Constraints.* Faces and edges can be marked constrained, which makes the
//! tetrahedralization a constrained one: a cavity never grows across a
//! constrained face, and it is shrunk (tets taken back out, never the one
//! holding the point) until every boundary face sees the point strictly and no
//! constrained face, constrained edge or vertex would be deleted. An insertion
//! that would have to delete a constraint is refused ([`DelaunayError::Blocked`])
//! and leaves the structure unchanged.

use crate::fxhash::FxHashSet;

use mesher_geom::predicates::{collinear, in_circumcircle_coplanar, insphere, orient3d};
use mesher_geom::{Point3, Sign};

/// Tets a cavity may grow to in a constrained tetrahedralization.
const MAX_CAVITY: usize = 4096;

/// Cavity rank of a tet that must not be taken out.
const KEEP: u8 = 3;

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
    /// The point could not be inserted without deleting a constraint; nothing
    /// changed.
    Blocked {
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
            DelaunayError::Blocked { point } => {
                write!(f, "point {point} would delete a constrained face or edge")
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

/// Where a straight walk between two points ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Walk {
    /// The end point is reached without crossing a constrained face; it is
    /// in (or on) the real tet at this slot.
    Clear(u32),
    /// The first constrained face crossed (sorted vertices).
    Constraint([u32; 3]),
    /// The walk left the convex hull or did not terminate.
    Lost,
}

/// A Delaunay tetrahedralization that accepts further points.
#[derive(Debug, Clone)]
pub struct Delaunay {
    points: Vec<Point3>,
    tets: Vec<Tet>,
    free: Vec<u32>,
    stamp: Vec<u32>,
    /// Per tet face (4 per slot): the insertion (epoch) that last worked out
    /// whether its point sees the face, and the answer (bit 0: seen, bit 1:
    /// on one of the point's planes). Flat arrays rather than maps: the
    /// cavity loop asks about every boundary face on every pass.
    sight_stamp: Vec<u32>,
    sight_bits: Vec<u8>,
    /// Per tet: the insertion (epoch) that took it out of the cavity.
    excluded_stamp: Vec<u32>,
    /// Per tet: its position in the current insertion's cavity order.
    cavity_rank: Vec<u32>,
    /// Per tet: the insertion (epoch) that last tested it for conflict, and
    /// the answer.
    conflict_stamp: Vec<u32>,
    conflict_bit: Vec<bool>,
    epoch: u32,
    alive: usize,
    last: u32,
    rng: SplitMix,
    duplicates: Vec<(u32, u32)>,
    /// Slots created by the last successful insertion.
    created: Vec<u32>,
    /// A live tet holding each vertex (any one), to find its star.
    vertex_tet: Vec<u32>,
    /// Constrained faces and edges, with sorted vertices.
    faces: FxHashSet<[u32; 3]>,
    edges: FxHashSet<(u32, u32)>,
    lost: Vec<[u32; 3]>,
    /// Plane labels per vertex (see `set_planes`).
    planes: Vec<Vec<u32>>,
}

impl Delaunay {
    /// Tetrahedralize `points` (inserted in BRIO order from `seed`).
    pub fn new(points: Vec<Point3>, seed: u64) -> Result<Delaunay, DelaunayError> {
        Delaunay::with_planes(points, seed, Vec::new())
    }

    /// Tetrahedralize `points` with plane labels (see `set_planes`) from the
    /// start, so no tet is flat to rounding within a plane.
    pub fn with_planes(
        points: Vec<Point3>,
        seed: u64,
        planes: Vec<Vec<u32>>,
    ) -> Result<Delaunay, DelaunayError> {
        let _prof = crate::prof::scope("delaunay::with_planes");
        if let Some(i) = points.iter().position(|p| !p.iter().all(|c| c.is_finite())) {
            return Err(DelaunayError::NonFinite { point: i as u32 });
        }
        let order = brio_order(&points, seed);
        let n = points.len();
        let mut d = Delaunay {
            points,
            tets: Vec::new(),
            free: Vec::new(),
            stamp: Vec::new(),
            sight_stamp: Vec::new(),
            sight_bits: Vec::new(),
            excluded_stamp: Vec::new(),
            cavity_rank: Vec::new(),
            conflict_stamp: Vec::new(),
            conflict_bit: Vec::new(),
            epoch: 0,
            alive: 0,
            last: 0,
            rng: SplitMix(seed ^ 0x9e37_79b9_7f4a_7c15),
            duplicates: Vec::new(),
            created: Vec::new(),
            vertex_tet: vec![NONE; n],
            faces: FxHashSet::default(),
            edges: FxHashSet::default(),
            lost: Vec::new(),
            planes,
        };
        let first = d.initial_tet(&order)?;
        for &i in &order {
            if first.contains(&i) {
                continue;
            }
            let outcome = match d.insert_index(i, &[], &[], false) {
                // Its plane labels can leave a point no cavity it sees
                // whole; insert it unlabelled (flat tets it makes are
                // removed afterwards), then label it.
                Err(DelaunayError::Blocked { .. }) if !d.planes(i).is_empty() => {
                    let planes = std::mem::take(&mut d.planes[i as usize]);
                    let outcome = d.insert_index(i, &[], &[], false);
                    d.planes[i as usize] = planes;
                    outcome?
                }
                other => other?,
            };
            if let Inserted::Duplicate(v) = outcome {
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
        self.insert_with(p, &[], &[], false)
    }

    /// Insert a point whose cavity must also hold the tets in `seeds` (taken
    /// out again if they cannot be kept) and `keep` (refused if they cannot).
    /// A point on a face or edge that is to be split seeds the cavity with the
    /// tets around it, after the face or edge is unconstrained.
    ///
    /// `force` ignores the face constraints (not the edge ones), as plain
    /// Bowyer–Watson does: the constrained faces the insertion deletes are
    /// unconstrained and listed in [`Delaunay::lost`].
    pub fn insert_with(
        &mut self,
        p: Point3,
        seeds: &[u32],
        keep: &[u32],
        force: bool,
    ) -> Result<(u32, Inserted), DelaunayError> {
        let _prof = crate::prof::scope("delaunay::insert_with");
        self.lost.clear();
        let i = self.points.len() as u32;
        if !p.iter().all(|c| c.is_finite()) {
            return Err(DelaunayError::NonFinite { point: i });
        }
        self.points.push(p);
        self.vertex_tet.push(NONE);
        let outcome = match self.insert_index(i, seeds, keep, force) {
            Ok(o) => o,
            Err(e) => {
                if matches!(e, DelaunayError::Blocked { .. }) {
                    self.points.pop();
                    self.vertex_tet.pop();
                }
                return Err(e);
            }
        };
        if let Inserted::Duplicate(v) = outcome {
            self.duplicates.push((i, v));
        }
        Ok((i, outcome))
    }

    /// The constrained faces (sorted) deleted by the last forced insertion.
    #[must_use]
    pub fn lost(&self) -> &[[u32; 3]] {
        &self.lost
    }

    /// Label vertex `v` with the planes (facet ids, ascending) it lies on.
    /// Points sharing a plane label are coplanar by construction, whatever
    /// rounding says: a point never sees a face on one of its planes, so no
    /// tet flat to rounding is made in a plane.
    pub fn set_planes(&mut self, v: u32, planes: Vec<u32>) {
        if self.planes.len() <= v as usize {
            self.planes.resize(v as usize + 1, Vec::new());
        }
        self.planes[v as usize] = planes;
    }

    /// The plane labels of vertex `v`.
    #[must_use]
    pub fn planes(&self, v: u32) -> &[u32] {
        self.planes.get(v as usize).map_or(&[], Vec::as_slice)
    }

    /// Do the vertices share a plane label?
    #[must_use]
    pub fn on_one_plane(&self, vs: &[u32]) -> bool {
        let Some((&first, rest)) = vs.split_first() else {
            return false;
        };
        self.planes(first).iter().any(|f| {
            rest.iter()
                .all(|&v| self.planes(v).binary_search(f).is_ok())
        })
    }

    /// Remove the tets whose four vertices share a plane label (flat by
    /// construction; rounding made them) by edge flips, a few passes.
    /// Returns how many remain.
    pub fn remove_flat_tets(&mut self) -> usize {
        let _prof = crate::prof::scope("delaunay::remove_flat_tets");
        for _ in 0..8 {
            let flat: Vec<u32> = self
                .real_slots()
                .into_iter()
                .filter(|&t| self.on_one_plane(&self.tets[t as usize].v))
                .collect();
            let mut changed = false;
            for t in &flat {
                let Some(v) = self.tet(*t) else {
                    continue;
                };
                if !self.on_one_plane(&v) {
                    continue;
                }
                for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
                    if self.remove_edge(v[a], v[b]) {
                        changed = true;
                        break;
                    }
                }
            }
            if !changed {
                break;
            }
        }
        self.real_slots()
            .into_iter()
            .filter(|&t| self.on_one_plane(&self.tets[t as usize].v))
            .count()
    }

    /// Mark a face constrained (vertex order does not matter).
    pub fn constrain_face(&mut self, f: [u32; 3]) {
        self.faces.insert(sorted3(f));
    }

    pub fn unconstrain_face(&mut self, f: [u32; 3]) {
        self.faces.remove(&sorted3(f));
    }

    #[must_use]
    pub fn is_constrained_face(&self, f: [u32; 3]) -> bool {
        self.faces.contains(&sorted3(f))
    }

    /// Mark an edge constrained.
    pub fn constrain_edge(&mut self, u: u32, v: u32) {
        self.edges.insert((u.min(v), u.max(v)));
    }

    pub fn unconstrain_edge(&mut self, u: u32, v: u32) {
        self.edges.remove(&(u.min(v), u.max(v)));
    }

    /// Slots of the live tets (ghosts included) holding vertex `v`.
    #[must_use]
    pub fn star(&self, v: u32) -> Vec<u32> {
        let _prof = crate::prof::scope("delaunay::star");
        let holds = |t: u32| {
            self.tets
                .get(t as usize)
                .is_some_and(|x| x.alive && x.v.contains(&v))
        };
        let mut start = self.vertex_tet.get(v as usize).copied().unwrap_or(NONE);
        if !holds(start) {
            match (0..self.tets.len() as u32).find(|&t| holds(t)) {
                Some(t) => start = t,
                None => return Vec::new(),
            }
        }
        let mut out = vec![start];
        let mut seen: FxHashSet<u32> = FxHashSet::default();
        seen.insert(start);
        let mut k = 0;
        while k < out.len() {
            let t = self.tets[out[k] as usize];
            k += 1;
            for j in 0..4 {
                // Faces holding v are those opposite the other vertices.
                if t.v[j] != v && seen.insert(t.n[j]) {
                    out.push(t.n[j]);
                }
            }
        }
        out
    }

    /// Slots of the tets (ghosts included) with face `f`.
    #[must_use]
    pub fn face_tets(&self, f: [u32; 3]) -> Vec<u32> {
        let _prof = crate::prof::scope("delaunay::face_tets");
        self.star(f[0])
            .into_iter()
            .filter(|&t| {
                let x = self.tets[t as usize];
                x.v.contains(&f[1]) && x.v.contains(&f[2])
            })
            .collect()
    }

    /// Slots of the tets (ghosts included) with edge `uv`.
    #[must_use]
    pub fn edge_tets(&self, u: u32, v: u32) -> Vec<u32> {
        let _prof = crate::prof::scope("delaunay::edge_tets");
        self.star(u)
            .into_iter()
            .filter(|&t| self.tets[t as usize].v.contains(&v))
            .collect()
    }

    /// Whether `slot` holds a live ghost.
    #[must_use]
    pub fn is_ghost_slot(&self, slot: u32) -> bool {
        self.tets
            .get(slot as usize)
            .is_some_and(|x| x.alive && x.v[3] == GHOST)
    }

    /// The edges (sorted pairs) crossing the open triangle `f` (its vertices
    /// in the tetrahedralization): found from the tets around its corners,
    /// spreading through tets that hold a crossing edge.
    #[must_use]
    pub fn crossing_edges(&self, f: [u32; 3]) -> Vec<(u32, u32)> {
        let _prof = crate::prof::scope("delaunay::crossing_edges");
        let [a, b, c] = f.map(|x| self.pos(x));
        let crosses = |x: u32, y: u32| {
            if x == GHOST || y == GHOST || f.contains(&x) || f.contains(&y) {
                return false;
            }
            let (px, py) = (self.pos(x), self.pos(y));
            // A point on the triangle's plane by label is on it.
            let side = |q: u32, pq: Point3| {
                if self.on_one_plane(&[f[0], f[1], f[2], q]) {
                    Sign::Zero
                } else {
                    orient3d(a, b, c, pq)
                }
            };
            let (sx, sy) = (side(x, px), side(y, py));
            if sx == Sign::Zero || sy == Sign::Zero || sx == sy {
                return false;
            }
            let s = [
                orient3d(px, py, a, b),
                orient3d(px, py, b, c),
                orient3d(px, py, c, a),
            ];
            s[0] != Sign::Zero && s[0] == s[1] && s[1] == s[2]
        };
        let mut seen: Vec<u32> = Vec::new();
        let mut queue: Vec<u32> = f.iter().flat_map(|&x| self.star(x)).collect();
        let mut out = Vec::new();
        while let Some(t) = queue.pop() {
            if seen.contains(&t) {
                continue;
            }
            seen.push(t);
            if seen.len() > 4096 {
                break;
            }
            let v = self.tets[t as usize].v;
            let mut any = false;
            for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
                if crosses(v[i], v[j]) {
                    any = true;
                    let e = (v[i].min(v[j]), v[i].max(v[j]));
                    if !out.contains(&e) {
                        out.push(e);
                    }
                }
            }
            if any {
                queue.extend(self.tets[t as usize].n);
            }
        }
        out.sort_unstable();
        out
    }

    /// Whether face `f` is in the tetrahedralization.
    #[must_use]
    pub fn has_face(&self, f: [u32; 3]) -> bool {
        let _prof = crate::prof::scope("delaunay::has_face");
        self.star(f[0]).into_iter().any(|t| {
            let x = self.tets[t as usize];
            x.v.contains(&f[1]) && x.v.contains(&f[2])
        })
    }

    /// Whether edge `uv` is in the tetrahedralization.
    #[must_use]
    pub fn has_edge(&self, u: u32, v: u32) -> bool {
        let _prof = crate::prof::scope("delaunay::has_edge");
        self.star(u)
            .into_iter()
            .any(|t| self.tets[t as usize].v.contains(&v))
    }

    /// Remove edge `uv` by a flip: the n tets around it (3 to 7, all real)
    /// are replaced by 2(n − 2) tets over a fan triangulation of the ring of
    /// their other vertices; n = 3 is the 3–2 flip. Every fan apex is tried,
    /// the one whose worst tet is fattest first. Refused (`false`,
    /// nothing changed) when the edge or a face around it is constrained or
    /// no fan is valid.
    pub fn remove_edge(&mut self, u: u32, v: u32) -> bool {
        self.remove_edge_for(u, v, None, false)
    }

    /// `remove_edge`, preferring fans that make face `want`; with `improve`,
    /// only if the worst new tet is fatter than the worst old one.
    pub fn remove_edge_for(
        &mut self,
        u: u32,
        v: u32,
        want: Option<[u32; 3]>,
        improve: bool,
    ) -> bool {
        let _prof = crate::prof::scope("delaunay::remove_edge_for");
        if u == v || self.edges.contains(&(u.min(v), u.max(v))) {
            return false;
        }
        let around = self.edge_tets(u, v);
        if !(3..=7).contains(&around.len()) || around.iter().any(|&t| self.is_ghost(t)) {
            return false;
        }
        // The other two vertices of each tet, chained into the ring.
        let mut pairs: Vec<[u32; 2]> = around
            .iter()
            .map(|&t| {
                let mut o = self.tets[t as usize]
                    .v
                    .iter()
                    .copied()
                    .filter(|&x| x != u && x != v);
                [o.next().unwrap_or(NONE), o.next().unwrap_or(NONE)]
            })
            .collect();
        let mut ring = pairs.swap_remove(0).to_vec();
        while !pairs.is_empty() {
            let last = ring[ring.len() - 1];
            let Some(k) = pairs.iter().position(|p| p.contains(&last)) else {
                return false;
            };
            let p = pairs.swap_remove(k);
            let next = if p[0] == last { p[1] } else { p[0] };
            if pairs.is_empty() {
                if next != ring[0] {
                    return false;
                }
            } else {
                ring.push(next);
            }
        }
        let n = ring.len();
        if n != around.len()
            || ring
                .iter()
                .any(|&r| self.faces.contains(&sorted3([u, v, r])))
        {
            return false;
        }
        let (pu, pv) = (self.pos(u), self.pos(v));
        let mut fans: Vec<(f64, Vec<[u32; 4]>)> = Vec::new();
        'apex: for k in 0..n {
            let mut new = Vec::with_capacity(2 * (n - 2));
            let mut worst = f64::INFINITY;
            for i in 1..n - 1 {
                let tri = [ring[k], ring[(k + i) % n], ring[(k + i + 1) % n]];
                let [a, b, c] = tri.map(|x| self.pos(x));
                let (su, sv) = (orient3d(a, b, c, pu), orient3d(a, b, c, pv));
                let pair = match (su, sv) {
                    (Sign::Positive, Sign::Negative) => {
                        [[tri[0], tri[1], tri[2], u], [tri[1], tri[0], tri[2], v]]
                    }
                    (Sign::Negative, Sign::Positive) => {
                        [[tri[1], tri[0], tri[2], u], [tri[0], tri[1], tri[2], v]]
                    }
                    _ => continue 'apex,
                };
                for t in pair {
                    if self.on_one_plane(&t) {
                        continue 'apex;
                    }
                    worst = worst.min(self.shape(t));
                }
                new.extend(pair);
            }
            fans.push((worst, new));
        }
        // Fans making `want` first, then the fattest worst tet; a fan that
        // folds over itself is refused.
        let makes = |new: &[[u32; 4]]| {
            want.is_some_and(|w| {
                let w = sorted3(w);
                new.iter().any(|t| (0..4).any(|j| sorted3(face(t, j)) == w))
            })
        };
        fans.sort_by(|a, b| makes(&b.1).cmp(&makes(&a.1)).then(b.0.total_cmp(&a.0)));
        if improve {
            let old = around
                .iter()
                .map(|&t| self.shape(self.tets[t as usize].v))
                .fold(f64::INFINITY, f64::min);
            fans.retain(|f| f.0 > old);
        }
        fans.into_iter()
            .any(|(_, new)| self.replace(&around, &new).is_ok())
    }

    /// Scale-free fatness of a real tet: volume over the cube of its rms edge
    /// length (zero when flat; about 0.118 for the regular tet).
    #[must_use]
    pub fn shape(&self, t: [u32; 4]) -> f64 {
        let p = t.map(|x| self.pos(x));
        let e = |i: usize, c: usize| p[i][c] - p[0][c];
        let vol = (e(1, 0) * (e(2, 1) * e(3, 2) - e(2, 2) * e(3, 1))
            - e(1, 1) * (e(2, 0) * e(3, 2) - e(2, 2) * e(3, 0))
            + e(1, 2) * (e(2, 0) * e(3, 1) - e(2, 1) * e(3, 0)))
        .abs()
            / 6.0;
        let sq: f64 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            .iter()
            .map(|&(i, j)| (0..3).map(|c| (p[i][c] - p[j][c]).powi(2)).sum::<f64>())
            .sum();
        let rms = (sq / 6.0).sqrt();
        if rms > 0.0 { vol / rms.powi(3) } else { 0.0 }
    }

    /// The 2–3 flip: replace the two real tets sharing unconstrained face `f`
    /// by three around the edge joining their apexes, if that edge crosses
    /// `f`; with `improve`, only if the worst new tet is fatter than both.
    pub fn flip_face(&mut self, f: [u32; 3], improve: bool) -> bool {
        let _prof = crate::prof::scope("delaunay::flip_face");
        if self.faces.contains(&sorted3(f)) {
            return false;
        }
        let pair = self.face_tets(f);
        let [t1, t2] = pair[..] else {
            return false;
        };
        if self.is_ghost(t1) || self.is_ghost(t2) {
            return false;
        }
        let apex = |t: u32| self.tets[t as usize].v.into_iter().find(|x| !f.contains(x));
        let (Some(d), Some(e)) = (apex(t1), apex(t2)) else {
            return false;
        };
        let mut new = Vec::with_capacity(3);
        for (a, b) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
            let [pa, pb, pd, pe] = [a, b, d, e].map(|x| self.pos(x));
            let t = match orient3d(pa, pb, pd, pe) {
                Sign::Positive => [a, b, d, e],
                Sign::Negative => [b, a, d, e],
                Sign::Zero => return false,
            };
            if self.on_one_plane(&t) {
                return false;
            }
            new.push(t);
        }
        if improve {
            let old = self
                .shape(self.tets[t1 as usize].v)
                .min(self.shape(self.tets[t2 as usize].v));
            if new
                .iter()
                .map(|&t| self.shape(t))
                .fold(f64::INFINITY, f64::min)
                <= old
            {
                return false;
            }
        }
        self.replace(&[t1, t2], &new).is_ok()
    }

    /// Replace the tets `old` by `new` (positive, filling the same region),
    /// linking the new tets to each other and to the outside.
    fn replace(&mut self, old: &[u32], new: &[[u32; 4]]) -> Result<(), ()> {
        let _prof = crate::prof::scope("delaunay::replace");
        // Same cyclic order (one face's outward orientation)?
        let same =
            |a: [u32; 3], b: [u32; 3]| (0..3).any(|r| [b[r], b[(r + 1) % 3], b[(r + 2) % 3]] == a);
        // Each outside face (outward from the old tets) and neighbour, with
        // the index of the neighbour's pointer back, found now: freed slots
        // are reused by the new tets.
        let mut outside: Vec<([u32; 3], [u32; 3], u32, usize)> = Vec::new();
        for &t in old {
            for j in 0..4 {
                let nb = self.tets[t as usize].n[j];
                if !old.contains(&nb) {
                    let back = (0..4)
                        .find(|&m| self.tets[nb as usize].n[m] == t)
                        .ok_or(())?;
                    let f = face(&self.tets[t as usize].v, j);
                    outside.push((sorted3(f), f, nb, back));
                }
            }
        }
        outside.sort_unstable();
        // Match the new tets' faces before changing anything: a shared face
        // must be oriented oppositely in its two tets (else they overlap), an
        // outside face as it was.
        let mut keys: Vec<([u32; 3], [u32; 3], usize, usize)> = Vec::with_capacity(4 * new.len());
        for (n, v) in new.iter().enumerate() {
            for j in 0..4 {
                let f = face(v, j);
                keys.push((sorted3(f), f, n, j));
            }
        }
        keys.sort_unstable();
        enum Link {
            Inner(usize, usize, usize, usize),
            Outer(usize, usize, u32, usize),
        }
        let mut links = Vec::with_capacity(keys.len());
        let mut k = 0;
        while k < keys.len() {
            let (key, f, n, j) = keys[k];
            if k + 1 < keys.len() && keys[k + 1].0 == key {
                let (_, f2, n2, j2) = keys[k + 1];
                if same(f, f2) || (k + 2 < keys.len() && keys[k + 2].0 == key) {
                    return Err(());
                }
                links.push(Link::Inner(n, j, n2, j2));
                k += 2;
                continue;
            }
            let at = outside
                .binary_search_by(|o| o.0.cmp(&key))
                .map_err(|_| ())?;
            let (_, was, nb, back) = outside[at];
            if !same(f, was) {
                return Err(());
            }
            links.push(Link::Outer(n, j, nb, back));
            k += 1;
        }
        if links
            .iter()
            .filter(|l| matches!(l, Link::Outer(..)))
            .count()
            != outside.len()
        {
            return Err(());
        }
        for &t in old {
            self.kill(t);
        }
        let ids: Vec<u32> = new.iter().map(|&v| self.alloc(v)).collect();
        for l in links {
            match l {
                Link::Inner(n, j, n2, j2) => {
                    self.tets[ids[n] as usize].n[j] = ids[n2];
                    self.tets[ids[n2] as usize].n[j2] = ids[n];
                }
                Link::Outer(n, j, nb, back) => {
                    self.tets[ids[n] as usize].n[j] = nb;
                    self.tets[nb as usize].n[back] = ids[n];
                }
            }
        }
        self.last = ids[0];
        self.created = ids;
        Ok(())
    }

    /// Restore the (constrained) Delaunay property by Lawson flips: a face
    /// that is not constrained, between two real tets, whose second apex lies
    /// inside the first tet's circumsphere is flipped by a 2–3 flip, or by a
    /// 3–2 flip of one of its edges when that edge has three tets around it.
    /// Faces that cannot flip stay. At most `max_flips`; returns the count.
    pub fn delaunize(&mut self, max_flips: usize) -> usize {
        let _prof = crate::prof::scope("delaunay::delaunize");
        let mut queue: Vec<[u32; 3]> = Vec::new();
        for t in self.real_slots() {
            let v = self.tets[t as usize].v;
            for j in 0..4 {
                queue.push(sorted3(face(&v, j)));
            }
        }
        queue.sort_unstable();
        queue.dedup();
        let mut flips = 0;
        while let Some(f) = queue.pop() {
            if flips >= max_flips {
                break;
            }
            if self.faces.contains(&f) {
                continue;
            }
            let pair = self.face_tets(f);
            let [t1, t2] = pair[..] else {
                continue;
            };
            if self.is_ghost(t1) || self.is_ghost(t2) {
                continue;
            }
            let v1 = self.tets[t1 as usize].v;
            let Some(e) = self.tets[t2 as usize]
                .v
                .into_iter()
                .find(|x| !f.contains(x))
            else {
                continue;
            };
            let [a, b, c, d] = v1.map(|x| self.pos(x));
            if insphere(a, b, c, d, self.pos(e)) != Sign::Positive {
                continue;
            }
            let done = self.flip_face(f, false)
                || [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])]
                    .into_iter()
                    .any(|(x, y)| self.edge_tets(x, y).len() == 3 && self.remove_edge(x, y));
            if done {
                flips += 1;
                for t in self.created.clone() {
                    let v = self.tets[t as usize].v;
                    for j in 0..4 {
                        queue.push(sorted3(face(&v, j)));
                    }
                }
            }
        }
        flips
    }

    /// Walk the straight segment from `from` (inside the real tet `slot`) to
    /// `to`, reporting the first constrained face it crosses.
    #[must_use]
    pub fn walk(&self, slot: u32, from: Point3, to: Point3) -> Walk {
        let _prof = crate::prof::scope("delaunay::walk");
        if !self
            .tets
            .get(slot as usize)
            .is_some_and(|x| x.alive && x.v[3] != GHOST)
        {
            return Walk::Lost;
        }
        let mut t = slot;
        for _ in 0..self.tets.len() + 16 {
            let v = self.tets[t as usize].v;
            if v[3] == GHOST {
                return Walk::Lost;
            }
            // Exit through a face that `to` is beyond and the segment's line
            // passes through (its boundary counts, for degenerate lines).
            let mut exit = None;
            for j in 0..4 {
                let f = face(&v, j);
                let [a, b, c] = f.map(|x| self.pos(x));
                if orient3d(a, b, c, to) != Sign::Positive {
                    continue;
                }
                let s = [
                    orient3d(from, to, a, b),
                    orient3d(from, to, b, c),
                    orient3d(from, to, c, a),
                ];
                let through = !(s.contains(&Sign::Positive) && s.contains(&Sign::Negative));
                if through || exit.is_none() {
                    exit = Some(j);
                }
                if through {
                    break;
                }
            }
            let Some(j) = exit else {
                return Walk::Clear(t);
            };
            if self.face_constrained(t, j) {
                return Walk::Constraint(sorted3(face(&v, j)));
            }
            t = self.tets[t as usize].n[j];
        }
        Walk::Lost
    }

    /// Slots of the real tets created by the last insertion that added a point.
    #[must_use]
    pub fn last_created(&self) -> Vec<u32> {
        self.created
            .iter()
            .copied()
            .filter(|&t| self.tet(t).is_some())
            .collect()
    }

    /// The vertices of the real tet in `slot`, if it is alive.
    #[must_use]
    pub fn tet(&self, slot: u32) -> Option<[u32; 4]> {
        self.tets
            .get(slot as usize)
            .filter(|t| t.alive && t.v[3] != GHOST)
            .map(|t| t.v)
    }

    /// Slots of every alive real tet, ascending.
    #[must_use]
    pub fn real_slots(&self) -> Vec<u32> {
        let _prof = crate::prof::scope("delaunay::real_slots");
        (0..self.tets.len() as u32)
            .filter(|&t| self.tet(t).is_some())
            .collect()
    }

    /// The real tetrahedra, positively oriented.
    #[must_use]
    pub fn tets(&self) -> Vec<[u32; 4]> {
        let _prof = crate::prof::scope("delaunay::tets");
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
        let id = if let Some(id) = self.free.pop() {
            self.tets[id as usize] = t;
            id
        } else {
            self.tets.push(t);
            self.stamp.push(0);
            self.sight_stamp.extend([0; 4]);
            self.sight_bits.extend([0; 4]);
            self.excluded_stamp.push(0);
            self.cavity_rank.push(0);
            self.conflict_stamp.push(0);
            self.conflict_bit.push(false);
            (self.tets.len() - 1) as u32
        };
        for x in v {
            if let Some(slot) = self.vertex_tet.get_mut(x as usize) {
                *slot = id;
            }
        }
        id
    }

    /// Whether point `i` (at `p`) sees face `j` of tet `t` from outside, and
    /// whether the face is on one of `i`'s planes; worked out once per
    /// insertion (`epoch`) and remembered.
    fn sight(&mut self, t: u32, j: usize, p: Point3, i: u32, epoch: u32) -> (bool, bool) {
        let slot = 4 * t as usize + j;
        if self.sight_stamp[slot] == epoch {
            let b = self.sight_bits[slot];
            return (b & 1 != 0, b & 2 != 0);
        }
        let f = face(&self.tets[t as usize].v, j);
        let (seen, planar) = if f.contains(&GHOST) {
            let r: Vec<Point3> = f
                .iter()
                .filter(|&&x| x != GHOST)
                .map(|&x| self.pos(x))
                .collect();
            (!collinear(r[0], r[1], p), false)
        } else {
            let [a, b, c] = f.map(|x| self.pos(x));
            // A face on one of p's planes is seen from neither side: growing
            // across it would only meet it again further on.
            let planar = self.on_one_plane(&[f[0], f[1], f[2], i]);
            (!planar && orient3d(a, b, c, p) == Sign::Negative, planar)
        };
        self.sight_stamp[slot] = epoch;
        self.sight_bits[slot] = u8::from(seen) | (u8::from(planar) << 1);
        (seen, planar)
    }

    /// Is face `j` of tet `t` constrained?
    fn face_constrained(&self, t: u32, j: usize) -> bool {
        if self.faces.is_empty() {
            return false;
        }
        let f = face(&self.tets[t as usize].v, j);
        !f.contains(&GHOST) && self.faces.contains(&sorted3(f))
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
            .find(|&(i, s)| s != Sign::Zero && !self.on_one_plane(&[a, b, c, i]))
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

    /// Does tet `t` conflict with point `i` (at `p`)? A ghost conflicts when
    /// the point is beyond its hull face, or on the face's plane (exactly, or
    /// by plane label) inside its circumcircle.
    fn conflict(&self, t: u32, p: Point3, i: u32) -> bool {
        let v = self.tets[t as usize].v;
        let [a, b, c] = [v[0], v[1], v[2]].map(|x| self.pos(x));
        if v[3] == GHOST {
            let side = if self.on_one_plane(&[v[0], v[1], v[2], i]) {
                Sign::Zero
            } else {
                orient3d(a, b, c, p)
            };
            match side {
                Sign::Positive => true,
                Sign::Negative => false,
                Sign::Zero => in_circumcircle_coplanar(a, b, c, p),
            }
        } else {
            insphere(a, b, c, self.pos(v[3]), p) == Sign::Positive
        }
    }

    /// `conflict`, remembered for the insertion (`epoch`): the cavity loop
    /// asks about the same tets on every pass.
    fn conflict_cached(&mut self, t: u32, p: Point3, i: u32, epoch: u32) -> bool {
        let k = t as usize;
        if self.conflict_stamp[k] != epoch {
            self.conflict_bit[k] = self.conflict(t, p, i);
            self.conflict_stamp[k] = epoch;
        }
        self.conflict_bit[k]
    }

    /// A tet whose closure contains `p`, or a ghost whose hull face `p` is beyond.
    fn locate(&mut self, p: Point3) -> u32 {
        let _prof = crate::prof::scope("delaunay::locate");
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
                .find(|&t| {
                    self.tets[t as usize].alive && self.is_ghost(t) && self.conflict(t, p, NONE)
                })
                .unwrap_or(t)
        })
    }

    /// A cavity tet to take out because the cavity would delete a constrained
    /// face (one between two cavity tets), a vertex (on no boundary face) or a
    /// constrained edge (on no boundary face). `Err` if only fixed tets could
    /// be taken out; `None` if nothing is deleted. The tet of lowest `rank`
    /// goes; rank `KEEP` never does.
    #[allow(clippy::too_many_arguments)]
    fn violation(
        &self,
        cavity: &[u32],
        boundary: &[(u32, usize, u32)],
        epoch: u32,
        rank: &dyn Fn(u32) -> u8,
        force: bool,
        p: Point3,
    ) -> Result<Option<u32>, ()> {
        let _prof = crate::prof::scope("delaunay::violation");
        let pick = |a: u32, b: Option<u32>| {
            let best = std::iter::once(a)
                .chain(b)
                .filter(|&t| rank(t) < KEEP)
                .min_by_key(|&t| rank(t));
            best.map(Some).ok_or(())
        };
        if !force && !self.faces.is_empty() {
            for &t in cavity {
                for j in 0..4 {
                    let u = self.tets[t as usize].n[j];
                    if self.stamp[u as usize] == epoch && self.face_constrained(t, j) {
                        return pick(t, Some(u));
                    }
                }
            }
        }
        let mut on_boundary: Vec<u32> = Vec::with_capacity(boundary.len() * 3);
        let mut edges_on: Vec<(u32, u32)> = Vec::new();
        for &(t, j, _) in boundary {
            let f = face(&self.tets[t as usize].v, j);
            on_boundary.extend(f.iter().copied().filter(|&x| x != GHOST));
            if !self.edges.is_empty() {
                for (a, b) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
                    if a != GHOST && b != GHOST {
                        edges_on.push((a.min(b), a.max(b)));
                    }
                }
            }
        }
        on_boundary.sort_unstable();
        edges_on.sort_unstable();
        // The cavity tet holding all of `of` that is best taken out: the one
        // behind them as seen from p (p beyond most of its faces through
        // them), so that the faces it exposes face p; then the lowest rank.
        let other = |of: &[u32]| {
            cavity
                .iter()
                .copied()
                .filter(|&c| of.iter().all(|x| self.tets[c as usize].v.contains(x)))
                .filter(|&c| rank(c) < KEEP)
                .min_by_key(|&c| {
                    let v = self.tets[c as usize].v;
                    let behind = (0..4)
                        .filter(|&j| {
                            let f = face(&v, j);
                            !f.contains(&GHOST) && of.iter().all(|x| f.contains(x)) && {
                                let [a, b, c] = f.map(|x| self.pos(x));
                                orient3d(a, b, c, p) == Sign::Positive
                            }
                        })
                        .count();
                    (std::cmp::Reverse(behind), rank(c))
                })
        };
        for &t in cavity {
            let v = self.tets[t as usize].v;
            for &x in &v {
                if x != GHOST && on_boundary.binary_search(&x).is_err() {
                    return pick(other(&[x]).unwrap_or(t), None);
                }
            }
            if self.edges.is_empty() {
                continue;
            }
            for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
                let (x, y) = (v[a], v[b]);
                if x == GHOST || y == GHOST {
                    continue;
                }
                let e = (x.min(y), x.max(y));
                if self.edges.contains(&e) && edges_on.binary_search(&e).is_err() {
                    return pick(other(&[x, y]).unwrap_or(t), None);
                }
            }
        }
        // The cone needs a closed 2-manifold boundary: every boundary edge on
        // exactly two boundary faces. Shrinking can pinch it at an edge.
        let mut all: Vec<(u32, u32)> = boundary
            .iter()
            .flat_map(|&(t, j, _)| {
                let f = face(&self.tets[t as usize].v, j);
                [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])].map(|(a, b)| (a.min(b), a.max(b)))
            })
            .collect();
        all.sort_unstable();
        let mut k = 0;
        while k < all.len() {
            let n = all[k..].iter().take_while(|&&e| e == all[k]).count();
            if n != 2 {
                let (x, y) = all[k];
                let of: Vec<u32> = [x, y].into_iter().filter(|&q| q != GHOST).collect();
                return match other(&of) {
                    Some(t) => pick(t, None),
                    None => Err(()),
                };
            }
            k += n;
        }
        Ok(None)
    }

    fn insert_index(
        &mut self,
        i: u32,
        seeds: &[u32],
        keep: &[u32],
        force: bool,
    ) -> Result<Inserted, DelaunayError> {
        let _prof = crate::prof::scope("delaunay::insert_index");
        let p = self.pos(i);
        let mut t0 = self.locate(p);
        if self.is_ghost(t0) && !self.conflict(t0, p, i) {
            // Defensive: never seed a cavity from a ghost p does not conflict with.
            match (0..self.tets.len() as u32)
                .find(|&t| self.tets[t as usize].alive && self.conflict(t, p, i))
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
        let blocked = DelaunayError::Blocked { point: i };

        // Cavity: the seeds, and tets whose circumsphere strictly contains p
        // reached from them without crossing a constrained face.
        self.epoch = self.epoch.wrapping_add(1);
        let epoch = self.epoch;
        let mut cavity = Vec::new();
        for &t in std::iter::once(&t0).chain(seeds).chain(keep) {
            if self.tets.get(t as usize).is_some_and(|x| x.alive) && self.stamp[t as usize] != epoch
            {
                self.stamp[t as usize] = epoch;
                cavity.push(t);
            }
        }
        let grow_scope = crate::prof::scope("delaunay::insert/grow");
        let mut k = 0;
        while k < cavity.len() {
            let t = cavity[k];
            k += 1;
            for j in 0..4 {
                let u = self.tets[t as usize].n[j];
                if self.stamp[u as usize] != epoch
                    && (force || !self.face_constrained(t, j))
                    && self.conflict(u, p, i)
                {
                    self.stamp[u as usize] = epoch;
                    cavity.push(u);
                }
            }
        }
        // Tets that must stay in the cavity, and tets taken out (never re-added).
        // The tet holding p may go too, but only as a last resort: a cavity
        // that does not hold p ends up empty, and the insertion is refused.
        let rank = |t: u32| {
            if keep.contains(&t) {
                KEEP
            } else if t == t0 {
                2
            } else if seeds.contains(&t) {
                1
            } else {
                0
            }
        };
        drop(grow_scope);
        // Position of each cavity tet in the cavity's order (ranks only grow,
        // so taking tets out keeps the order of the rest).
        let mut next_rank = 0u32;
        for &t in &cavity {
            self.cavity_rank[t as usize] = next_rank;
            next_rank += 1;
        }
        let boundary_scope = crate::prof::scope("delaunay::insert/boundary");
        let constrained = !self.faces.is_empty() || !self.edges.is_empty();
        let grow_limit = (4 * cavity.len()).max(64);
        let interior = seeds.is_empty() && keep.is_empty() && self.planes(i).is_empty();
        // Whether p sees each face, and whether the face is on one of its
        // planes: fixed for the insertion, so worked out once.

        // Boundary faces must see p: grow the cavity across any that do not,
        // or, at a constraint, shrink it. Then no constrained face, constrained
        // edge or vertex may be left inside.
        // Boundary faces p does not see, in boundary order (cavity order,
        // then face index): the only faces a pass looks at. Rebuilt in full
        // after growth; after a tet is taken out only its own faces leave and
        // its cavity neighbours' faces towards it join, so the list is
        // updated in place.
        let mut unseen: Vec<(u32, usize, u32)> = Vec::new();
        let mut rebuild = true;
        let boundary = loop {
            let _iter = crate::prof::scope("delaunay::insert/boundary-pass");
            if rebuild {
                unseen.clear();
                for &t in &cavity {
                    for j in 0..4 {
                        let u = self.tets[t as usize].n[j];
                        if self.stamp[u as usize] != epoch && !self.sight(t, j, p, i, epoch).0 {
                            unseen.push((t, j, u));
                        }
                    }
                }
            }
            // Every boundary face p does not see, at once: grow across it,
            // or (at a constraint, a face on one of p's planes, or a tet taken
            // out before) take its cavity tet out.
            // Only the first tet to take out matters (one goes per pass),
            // but every face is still checked for growth and for a kept tet.
            let mut grow = Vec::new();
            let mut first_shrink: Option<u32> = None;
            for &(t, j, u) in &unseen {
                let planar = self.sight_bits[4 * t as usize + j] & 2 != 0;
                // Growing across a face p does not see is how plain
                // Bowyer–Watson handles cospherical points. With constraints
                // (unless forced) the mesh is not Delaunay and such faces are
                // common: growth is allowed locally (it links the point to
                // the tet it refines), then the cavity only shrinks.
                let may_grow = force
                    || !constrained
                    || (interior && cavity.len() < grow_limit)
                    || self.conflict_cached(u, p, i, epoch);
                if may_grow
                    && !planar
                    && (force || !self.face_constrained(t, j))
                    && self.excluded_stamp[u as usize] != epoch
                {
                    if self.stamp[u as usize] != epoch {
                        self.stamp[u as usize] = epoch;
                        grow.push(u);
                    }
                } else if rank(t) == KEEP {
                    return Err(blocked);
                } else if first_shrink.is_none() {
                    first_shrink = Some(t);
                }
            }
            // Grow all at once, but take tets out one at a time: one removal
            // can expose faces that make others unnecessary.
            let mut shrink: Vec<u32> = Vec::new();
            if grow.is_empty() {
                shrink.extend(first_shrink);
            }
            if grow.is_empty() && shrink.is_empty() {
                let mut boundary: Vec<(u32, usize, u32)> = Vec::new();
                for &t in &cavity {
                    for j in 0..4 {
                        let u = self.tets[t as usize].n[j];
                        if self.stamp[u as usize] != epoch {
                            boundary.push((t, j, u));
                        }
                    }
                }
                match self.violation(&cavity, &boundary, epoch, &rank, force, p) {
                    Ok(Some(t)) => shrink.push(t),
                    Ok(None) => break boundary,
                    Err(()) => return Err(blocked),
                }
            }
            if cavity.len() + grow.len() >= self.alive {
                return Err(fail);
            }
            // A cavity this big means growth ran away through a part of the
            // mesh that constraints keep from being Delaunay.
            if constrained && cavity.len() + grow.len() > MAX_CAVITY {
                return Err(blocked);
            }
            rebuild = !grow.is_empty();
            for &g in &grow {
                self.cavity_rank[g as usize] = next_rank;
                next_rank += 1;
            }
            cavity.extend(grow);
            for victim in shrink {
                cavity.retain(|&x| x != victim);
                self.stamp[victim as usize] = 0;
                self.excluded_stamp[victim as usize] = epoch;
                if rebuild {
                    continue;
                }
                unseen.retain(|&(t, _, _)| t != victim);
                let around = self.tets[victim as usize].n;
                for w in around {
                    if self.stamp[w as usize] != epoch {
                        continue;
                    }
                    for j in 0..4 {
                        if self.tets[w as usize].n[j] != victim || self.sight(w, j, p, i, epoch).0 {
                            continue;
                        }
                        let key = (self.cavity_rank[w as usize], j);
                        let at = unseen
                            .partition_point(|&(t, k, _)| (self.cavity_rank[t as usize], k) < key);
                        unseen.insert(at, (w, j, victim));
                    }
                }
            }
            if cavity.is_empty() {
                return Err(blocked);
            }
        };
        drop(boundary_scope);
        let _cone_scope = crate::prof::scope("delaunay::insert/cone");
        if boundary.is_empty() {
            return Err(fail);
        }
        if force && !self.faces.is_empty() {
            for &t in &cavity {
                for j in 0..4 {
                    let u = self.tets[t as usize].n[j];
                    if self.stamp[u as usize] == epoch && self.face_constrained(t, j) {
                        let f = sorted3(face(&self.tets[t as usize].v, j));
                        if !self.lost.contains(&f) {
                            self.lost.push(f);
                        }
                    }
                }
            }
            for f in &self.lost {
                self.faces.remove(f);
            }
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
        self.created = created;
        Ok(Inserted::New)
    }
}

/// The face opposite local vertex `k`, outward for a positive tet.
fn face(v: &[u32; 4], k: usize) -> [u32; 3] {
    FACES[k].map(|i| v[i])
}

fn sorted3(mut f: [u32; 3]) -> [u32; 3] {
    f.sort_unstable();
    f
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
    fn edge_removal_keeps_a_valid_mesh() {
        let mut rng = SplitMix(42);
        let mut unit = || (rng.next() >> 11) as f64 / (1u64 << 53) as f64;
        let pts: Vec<Point3> = (0..300).map(|_| [unit(), unit(), unit()]).collect();
        let mut d = Delaunay::new(pts, 5).unwrap();
        let mut removed = 0;
        for round in 0..2000u32 {
            let tets = d.tets();
            let t = tets[(round as usize * 7919) % tets.len()];
            let (u, v) = (t[(round % 4) as usize], t[((round / 4) % 4) as usize]);
            if u == v {
                continue;
            }
            let before = d.tets().len();
            let vol = |d: &Delaunay| -> f64 {
                d.tets()
                    .iter()
                    .map(|t| {
                        let p = t.map(|x| d.points()[x as usize]);
                        let e = |i: usize, c: usize| p[i][c] - p[0][c];
                        e(1, 0) * (e(2, 1) * e(3, 2) - e(2, 2) * e(3, 1))
                            - e(1, 1) * (e(2, 0) * e(3, 2) - e(2, 2) * e(3, 0))
                            + e(1, 2) * (e(2, 0) * e(3, 1) - e(2, 1) * e(3, 0))
                    })
                    .sum()
            };
            let v0 = vol(&d);
            let star = d.edge_tets(u, v);
            if d.remove_edge(u, v) {
                assert!((3..=7).contains(&star.len()));
                removed += 1;
                d.check().unwrap();
                assert!(d.tets().len() >= before - 1);
                assert!((vol(&d) - v0).abs() < 1e-9);
                assert!(!d.has_edge(u, v));
            }
        }
        assert!(removed > 10, "{removed} edges removed");
    }

    fn random_points(n: usize, seed: u64) -> Vec<Point3> {
        let mut rng = SplitMix(seed);
        let mut unit = || (rng.next() >> 11) as f64 / (1u64 << 53) as f64;
        (0..n).map(|_| [unit(), unit(), unit()]).collect()
    }

    #[test]
    fn constrained_faces_survive_insertion() {
        let mut d = Delaunay::new(random_points(200, 3), 9).unwrap();
        // Constrain every face of a few tets, then insert points around them.
        let held: Vec<[u32; 3]> = d
            .tets()
            .iter()
            .step_by(37)
            .flat_map(|t| (0..4).map(|j| face(t, j)).collect::<Vec<_>>())
            .collect();
        for &f in &held {
            d.constrain_face(f);
        }
        let mut blocked = 0;
        for p in random_points(300, 11) {
            let before = (d.points().len(), d.tets().len());
            match d.insert(p) {
                Ok(_) => {}
                Err(DelaunayError::Blocked { .. }) => {
                    blocked += 1;
                    assert_eq!((d.points().len(), d.tets().len()), before);
                }
                Err(e) => panic!("{e}"),
            }
        }
        d.check().unwrap();
        for f in held {
            assert!(d.has_face(f), "constrained face {f:?} was deleted");
        }
        assert!(blocked < 300);
    }

    #[test]
    fn face_flips_keep_a_valid_mesh() {
        let mut d = Delaunay::new(random_points(300, 21), 2).unwrap();
        let mut flips = 0;
        for round in 0..3000usize {
            let tets = d.tets();
            let t = tets[(round * 7919) % tets.len()];
            let f = face(&t, round % 4);
            let before = d.tets().len();
            if d.flip_face(f, false) {
                flips += 1;
                d.check().unwrap();
                assert_eq!(d.tets().len(), before + 1);
                assert!(!d.has_face(f));
            }
        }
        assert!(flips > 10, "{flips} faces flipped");
    }

    #[test]
    fn plane_labels_prevent_flat_tets() {
        // A 6×6 grid on a tilted plane (coplanar only to rounding) and points
        // on both sides.
        let (s, c) = (0.37f64.sin(), 0.37f64.cos());
        let mut pts = Vec::new();
        let mut planes = Vec::new();
        for i in 0..6 {
            for j in 0..6 {
                let (x, y) = (f64::from(i) / 5.0, f64::from(j) / 5.0);
                pts.push([x, c * y, s * y]);
                planes.push(vec![7]);
            }
        }
        for p in random_points(60, 5) {
            pts.push([
                p[0],
                c * p[1] - s * (p[2] - 0.5),
                s * p[1] + c * (p[2] - 0.5),
            ]);
            planes.push(Vec::new());
        }
        let plain = Delaunay::new(pts.clone(), 4).unwrap();
        let on = |t: &[u32; 4]| t.iter().all(|&v| v < 36);
        let labelled = Delaunay::with_planes(pts, 4, planes).unwrap();
        labelled.check().unwrap();
        assert_eq!(labelled.tets().iter().filter(|t| on(t)).count(), 0);
        // Without labels, rounding makes some (this is what labels are for).
        let _ = plain.tets().iter().filter(|t| on(t)).count();
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
