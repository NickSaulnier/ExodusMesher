//! An exact incremental Delaunay triangulation of points on one plane: the
//! 2D triangulation of a planar region that face refinement works from.
//!
//! The points are implicit (`mesher_geom::implicit::Point`) and all lie on
//! one plane exactly. Orientation is the exact orientation projected along
//! the axis the plane is not parallel to; the Delaunay test is the true
//! in-circle test in the plane ([`implicit::incircle_coplanar`]), not a
//! projected one (projection distorts circles). The structure follows the 3D
//! one ([`crate::cdt::ExactDelaunay`]): ghost triangles outside the hull, a
//! visibility walk, a cavity of strict in-circle conflicts grown across
//! every edge the point does not strictly see, and new triangles made by
//! putting the point in place of the vertex opposite each boundary edge,
//! which keeps the orientation.

use mesher_geom::Sign;
use mesher_geom::implicit::{self, Point};

/// The vertex at infinity, last in a ghost triangle.
const GHOST: u32 = u32::MAX;
const NONE: u32 = u32::MAX;

#[derive(Debug, Clone, Copy)]
struct Tri {
    /// Local vertices, counterclockwise in the projection; a ghost is
    /// `[a, b, GHOST]` with the hull edge `b a` of the real triangle beside
    /// it (the outside on the left of `a → b`).
    v: [u32; 3],
    /// Across the edge opposite `v[k]`.
    n: [u32; 3],
    alive: bool,
}

/// The Delaunay triangulation of points on a plane.
#[derive(Debug, Clone)]
pub(crate) struct PlaneDelaunay {
    axis: usize,
    pts: Vec<Point>,
    /// Each local point's id outside (a vertex of the 3D tetrahedralization).
    ids: Vec<u32>,
    /// The local point of each id outside.
    local_of: std::collections::BTreeMap<u32, u32>,
    /// A live triangle at each local point (`NONE` until it has one).
    vertex_tri: Vec<u32>,
    tris: Vec<Tri>,
    free: Vec<u32>,
    last: u32,
    /// Points waiting for three that are not collinear.
    pending: Vec<u32>,
    stamp: Vec<u32>,
    epoch: u32,
    /// Triangle slots filled since the last insertion began.
    created: Vec<u32>,
}

impl PlaneDelaunay {
    /// An empty triangulation of points on a plane not parallel to `axis`
    /// (the projection drops coordinate `axis`).
    pub(crate) fn new(axis: usize) -> PlaneDelaunay {
        PlaneDelaunay {
            axis,
            pts: Vec::new(),
            ids: Vec::new(),
            local_of: std::collections::BTreeMap::new(),
            vertex_tri: Vec::new(),
            tris: Vec::new(),
            free: Vec::new(),
            last: NONE,
            pending: Vec::new(),
            stamp: Vec::new(),
            epoch: 0,
            created: Vec::new(),
        }
    }

    fn orient(&self, a: u32, b: u32, c: u32) -> Sign {
        implicit::orient2d_projected(
            &self.pts[a as usize],
            &self.pts[b as usize],
            &self.pts[c as usize],
            self.axis,
        )
    }

    fn orient_p(&self, a: u32, b: u32, p: &Point) -> Sign {
        implicit::orient2d_projected(&self.pts[a as usize], &self.pts[b as usize], p, self.axis)
    }

    fn is_ghost(&self, t: u32) -> bool {
        self.tris[t as usize].v[2] == GHOST
    }

    fn alloc(&mut self, v: [u32; 3]) -> u32 {
        let t = Tri {
            v,
            n: [NONE; 3],
            alive: true,
        };
        let id = if let Some(id) = self.free.pop() {
            self.tris[id as usize] = t;
            id
        } else {
            self.tris.push(t);
            self.stamp.push(0);
            (self.tris.len() - 1) as u32
        };
        self.created.push(id);
        for x in v {
            if let Some(slot) = self.vertex_tri.get_mut(x as usize) {
                *slot = id;
            }
        }
        id
    }

    /// Pair up the edges the triangles `new` share (the first triangle and
    /// its three ghosts).
    fn link(&mut self, new: &[u32]) {
        let mut keys: Vec<((u32, u32), u32, usize)> = Vec::new();
        for &t in new {
            let v = self.tris[t as usize].v;
            for k in 0..3 {
                let (a, b) = (v[(k + 1) % 3], v[(k + 2) % 3]);
                keys.push(((a.min(b), a.max(b)), t, k));
            }
        }
        keys.sort_unstable();
        let mut i = 0;
        while i < keys.len() {
            if i + 1 < keys.len() && keys[i].0 == keys[i + 1].0 {
                let ((_, ta, ka), (_, tb, kb)) = (
                    (keys[i].0, keys[i].1, keys[i].2),
                    (keys[i + 1].0, keys[i + 1].1, keys[i + 1].2),
                );
                self.tris[ta as usize].n[ka] = tb;
                self.tris[tb as usize].n[kb] = ta;
                i += 2;
            } else {
                i += 1;
            }
        }
    }

    /// Start the next point location from a triangle at the vertex with id
    /// `id` outside, if it has one: a hint when the next point is near it.
    pub(crate) fn hint(&mut self, id: u32) {
        let Some(&i) = self.local_of.get(&id) else {
            return;
        };
        let t = self.vertex_tri[i as usize];
        if self.tris.get(t as usize).is_some_and(|x| x.alive) {
            self.last = t;
        }
    }

    /// Add point `p` (id `id` outside). Returns `false` if it equals a point
    /// already there (nothing changes).
    pub(crate) fn insert(&mut self, id: u32, p: Point) -> bool {
        self.created.clear();
        let _s = crate::prof::scope("plane insert");
        // Before the first triangle a duplicate is found here; after it,
        // by the insertion.
        if self.last == NONE && self.pts.iter().any(|q| implicit::coincident(q, &p)) {
            return false;
        }
        let i = self.pts.len() as u32;
        self.pts.push(p);
        self.ids.push(id);
        self.vertex_tri.push(NONE);
        if self.last == NONE {
            self.local_of.insert(id, i);
            self.pending.push(i);
            self.start();
            return true;
        }
        if self.insert_index(i) {
            self.local_of.insert(id, i);
            true
        } else {
            self.pts.pop();
            self.ids.pop();
            self.vertex_tri.pop();
            false
        }
    }

    /// The first triangle, once three pending points are not collinear;
    /// then the rest of the pending points.
    fn start(&mut self) {
        let n = self.pending.len();
        if n < 3 {
            return;
        }
        let (a, b) = (self.pending[0], self.pending[1]);
        let Some(&c) = self.pending[2..]
            .iter()
            .find(|&&c| self.orient(a, b, c) != Sign::Zero)
        else {
            return;
        };
        let v = if self.orient(a, b, c) == Sign::Positive {
            [a, b, c]
        } else {
            [b, a, c]
        };
        let t0 = self.alloc(v);
        // Ghosts across each edge: edge opposite v[k] is v[k+1] → v[k+2]
        // in the real triangle; the ghost has it reversed.
        let mut ids = vec![t0];
        for k in 0..3 {
            let (x, y) = (v[(k + 1) % 3], v[(k + 2) % 3]);
            ids.push(self.alloc([y, x, GHOST]));
        }
        self.link(&ids);
        self.last = t0;
        let rest: Vec<u32> = self
            .pending
            .iter()
            .copied()
            .filter(|&x| x != a && x != b && x != c)
            .collect();
        self.pending.clear();
        // A duplicate among the first points stays an unused local point, in
        // no triangle.
        for x in rest {
            self.insert_index(x);
        }
    }

    /// Does triangle `t` conflict with `p`? A real one when `p` is strictly
    /// inside its circumcircle; a ghost when `p` is strictly outside its hull
    /// edge, or on the edge's line strictly between its ends.
    fn conflict(&self, t: u32, p: &Point) -> bool {
        let v = self.tris[t as usize].v;
        if v[2] == GHOST {
            match self.orient_p(v[0], v[1], p) {
                Sign::Positive => true,
                Sign::Negative => false,
                Sign::Zero => {
                    implicit::dot_sign(&self.pts[v[0] as usize], p, &self.pts[v[1] as usize])
                        == Sign::Negative
                }
            }
        } else {
            let [a, b, c] = v.map(|x| &self.pts[x as usize]);
            implicit::incircle_coplanar(a, b, c, p) == Sign::Positive
        }
    }

    /// Does `p` strictly see edge `k` of cavity triangle `t`, so that the
    /// triangle made by putting `p` in place of `v[k]` is valid?
    fn sees(&self, t: u32, k: usize, p: &Point) -> bool {
        let v = self.tris[t as usize].v;
        let (a, b) = (v[(k + 1) % 3], v[(k + 2) % 3]);
        if a == GHOST || b == GHOST {
            // A new ghost: valid unless p is its real vertex.
            let r = if a == GHOST { b } else { a };
            !implicit::coincident(&self.pts[r as usize], p)
        } else {
            self.orient_p(a, b, p) == Sign::Positive
        }
    }

    /// A triangle whose closure holds `p` (a ghost if `p` is outside the
    /// hull), by walking from the last one.
    fn locate(&self, p: &Point) -> u32 {
        let mut t = self.last;
        if !self.tris.get(t as usize).is_some_and(|x| x.alive) {
            t = self.tris.iter().position(|x| x.alive).unwrap_or(0) as u32;
        }
        let limit = 4 * self.tris.len() + 16;
        let mut rot = (self.pts.len() as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        for _ in 0..limit {
            if self.is_ghost(t) {
                return t;
            }
            let v = self.tris[t as usize].v;
            rot = rot.rotate_left(7) ^ 0x2545_f491_4f6c_dd1d;
            let start = (rot % 3) as usize;
            let mut moved = false;
            for d in 0..3 {
                let k = (start + d) % 3;
                let (a, b) = (v[(k + 1) % 3], v[(k + 2) % 3]);
                if self.orient_p(a, b, p) == Sign::Negative {
                    t = self.tris[t as usize].n[k];
                    moved = true;
                    break;
                }
            }
            if !moved {
                return t;
            }
        }
        // Fall back to a scan (never expected).
        (0..self.tris.len() as u32)
            .find(|&t| self.tris[t as usize].alive && self.conflict(t, p))
            .unwrap_or(t)
    }

    fn insert_index(&mut self, i: u32) -> bool {
        let p = self.pts[i as usize];
        let s = crate::prof::scope("plane locate");
        let t0 = self.locate(&p);
        drop(s);
        let v0 = self.tris[t0 as usize].v;
        if v0
            .iter()
            .any(|&x| x != GHOST && implicit::coincident(&self.pts[x as usize], &p))
        {
            return false;
        }
        let start = if self.conflict(t0, &p) {
            t0
        } else {
            match (0..self.tris.len() as u32)
                .find(|&t| self.tris[t as usize].alive && self.conflict(t, &p))
            {
                Some(t) => t,
                None => return false,
            }
        };
        self.epoch = self.epoch.wrapping_add(1).max(1);
        let epoch = self.epoch;
        let mut cavity = vec![start];
        self.stamp[start as usize] = epoch;
        let mut k = 0;
        while k < cavity.len() {
            let t = cavity[k];
            k += 1;
            for j in 0..3 {
                let u = self.tris[t as usize].n[j];
                if self.stamp[u as usize] != epoch && self.conflict(u, &p) {
                    self.stamp[u as usize] = epoch;
                    cavity.push(u);
                }
            }
        }
        loop {
            let mut grow = Vec::new();
            for &t in &cavity {
                for j in 0..3 {
                    let u = self.tris[t as usize].n[j];
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
            if cavity.len() > self.tris.len() {
                return false;
            }
        }
        // New triangles: p in place of the vertex opposite each boundary
        // edge; then link them with each other and with the outside.
        let mut new = Vec::new();
        let mut outside: Vec<(u32, u32)> = Vec::new();
        for &t in &cavity {
            let tri = self.tris[t as usize];
            for k in 0..3 {
                let u = tri.n[k];
                if self.stamp[u as usize] == epoch {
                    continue;
                }
                let mut v = tri.v;
                v[k] = i;
                // Keep the vertex at infinity last.
                if v[0] == GHOST {
                    v = [v[1], v[2], v[0]];
                } else if v[1] == GHOST {
                    v = [v[2], v[0], v[1]];
                }
                let id = self.alloc(v);
                new.push(id);
                outside.push((id, u));
            }
        }
        for &t in &cavity {
            self.tris[t as usize].alive = false;
            self.free.push(t);
        }
        // Each new triangle's edge without `i` faces the outside triangle.
        for &(id, u) in &outside {
            let v = self.tris[id as usize].v;
            let k = v.iter().position(|&x| x == i).unwrap_or(0);
            self.tris[id as usize].n[k] = u;
            let (a, b) = (v[(k + 1) % 3], v[(k + 2) % 3]);
            let key = (a.min(b), a.max(b));
            let ut = self.tris[u as usize];
            for m in 0..3 {
                let (x, y) = (ut.v[(m + 1) % 3], ut.v[(m + 2) % 3]);
                if (x.min(y), x.max(y)) == key {
                    self.tris[u as usize].n[m] = id;
                }
            }
        }
        // Edges between new triangles: those through `i`.
        let mut keys: Vec<((u32, u32), u32, usize)> = Vec::new();
        for &id in &new {
            let v = self.tris[id as usize].v;
            for k in 0..3 {
                if v[k] == i {
                    continue;
                }
                let (a, b) = (v[(k + 1) % 3], v[(k + 2) % 3]);
                keys.push(((a.min(b), a.max(b)), id, k));
            }
        }
        keys.sort_unstable();
        for pair in keys.chunks(2) {
            if let [(ka, ta, ja), (kb, tb, jb)] = *pair {
                if ka == kb {
                    self.tris[ta as usize].n[ja] = tb;
                    self.tris[tb as usize].n[jb] = ta;
                }
            }
        }
        self.last = new
            .iter()
            .copied()
            .find(|&t| !self.is_ghost(t))
            .unwrap_or(new[0]);
        true
    }

    /// Is there a real triangle with these outside vertex ids (any order)?
    /// A walk around the first vertex.
    pub(crate) fn has_triangle(&self, v: [u32; 3]) -> bool {
        let locals: Vec<u32> = v
            .iter()
            .filter_map(|id| self.local_of.get(id).copied())
            .collect();
        let &[a, b, c] = locals.as_slice() else {
            return false;
        };
        let Some(&start) = self.vertex_tri.get(a as usize) else {
            return false;
        };
        if start == NONE
            || !self.tris[start as usize].alive
            || !self.tris[start as usize].v.contains(&a)
        {
            return false;
        }
        // The triangles around a: across each edge through a.
        let mut seen = vec![start];
        let mut k = 0;
        while k < seen.len() {
            let t = self.tris[seen[k] as usize];
            k += 1;
            if t.v[2] != GHOST && t.v.contains(&b) && t.v.contains(&c) {
                return true;
            }
            for j in 0..3 {
                if t.v[j] != a {
                    let u = t.n[j];
                    if u != NONE && !seen.contains(&u) {
                        seen.push(u);
                    }
                }
            }
        }
        false
    }

    /// Is there a first triangle yet? After that, every triangle an
    /// insertion makes has the new point as a vertex.
    pub(crate) fn started(&self) -> bool {
        self.last != NONE
    }

    /// The triangle slots the last insertion filled.
    pub(crate) fn created(&self) -> &[u32] {
        &self.created
    }

    /// Triangle slots so far (live or not).
    pub(crate) fn slots(&self) -> usize {
        self.tris.len()
    }

    /// Slot `t`, if it holds a live real triangle: its vertices' outside ids
    /// (counterclockwise in the projection) and the slots across the edges
    /// opposite them.
    pub(crate) fn real(&self, t: u32) -> Option<([u32; 3], [u32; 3])> {
        let tri = self.tris.get(t as usize)?;
        (tri.alive && tri.v[2] != GHOST).then(|| (tri.v.map(|x| self.ids[x as usize]), tri.n))
    }

    /// The real triangles at the vertex with outside id `id`: each one's slot
    /// (see [`PlaneDelaunay::real`]) and its
    /// vertices' outside ids. Empty for a vertex not in the triangulation.
    pub(crate) fn triangles_at(&self, id: u32) -> Vec<(u32, [u32; 3])> {
        let Some(&a) = self.local_of.get(&id) else {
            return Vec::new();
        };
        let Some(&start) = self.vertex_tri.get(a as usize) else {
            return Vec::new();
        };
        if start == NONE
            || !self.tris[start as usize].alive
            || !self.tris[start as usize].v.contains(&a)
        {
            return Vec::new();
        }
        let mut seen = vec![start];
        let mut k = 0;
        while k < seen.len() {
            let t = self.tris[seen[k] as usize];
            k += 1;
            for j in 0..3 {
                if t.v[j] != a {
                    let u = t.n[j];
                    if u != NONE && !seen.contains(&u) {
                        seen.push(u);
                    }
                }
            }
        }
        seen.into_iter()
            .filter(|&t| self.tris[t as usize].v[2] != GHOST)
            .map(|t| (t, self.tris[t as usize].v.map(|x| self.ids[x as usize])))
            .collect()
    }

    /// The real triangles, counterclockwise in the projection, with their
    /// vertices' outside ids, and each one's neighbours (by position in this
    /// list; `None` across the hull).
    #[cfg(test)]
    pub(crate) fn triangles(&self) -> Vec<([u32; 3], [Option<usize>; 3])> {
        let mut index = vec![usize::MAX; self.tris.len()];
        let mut out = Vec::new();
        for (t, tri) in self.tris.iter().enumerate() {
            if tri.alive && tri.v[2] != GHOST {
                index[t] = out.len();
                out.push((tri.v.map(|x| self.ids[x as usize]), tri.n));
            }
        }
        out.into_iter()
            .map(|(v, n)| {
                let n = n.map(|u| index.get(u as usize).copied().filter(|&x| x != usize::MAX));
                (v, n)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;

    struct Rng(u64);

    impl Rng {
        fn unit(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// Every triangle positive and empty of other points' circumcircles,
    /// neighbours mutual, and the triangle count right for the hull
    /// (2n − h − 2).
    fn check(pd: &PlaneDelaunay, points: &[Point]) {
        let tris = pd.triangles();
        let at = |id: u32| &points[id as usize];
        for (i, (v, n)) in tris.iter().enumerate() {
            let [a, b, c] = v.map(at);
            assert_eq!(
                implicit::orient2d_projected(a, b, c, pd.axis),
                Sign::Positive
            );
            for (j, q) in points.iter().enumerate() {
                if v.contains(&(j as u32)) {
                    continue;
                }
                assert_ne!(
                    implicit::incircle_coplanar(a, b, c, q),
                    Sign::Positive,
                    "{v:?} {j}"
                );
            }
            for u in n.iter().flatten() {
                assert!(tris[*u].1.contains(&Some(i)));
            }
        }
        let hull = tris
            .iter()
            .flat_map(|(_, n)| n.iter())
            .filter(|x| x.is_none())
            .count();
        assert_eq!(tris.len(), 2 * points.len() - hull - 2);
    }

    #[test]
    fn random_points_on_a_tilted_plane() {
        // A plane through three explicit points; BPT points on it.
        let tri = [[0.0, 0.0, 0.0], [3.0, 0.5, 1.0], [0.5, 2.0, -1.5]];
        for seed in 1..=5u64 {
            let mut rng = Rng(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1);
            let points: Vec<Point> = (0..200)
                .map(|_| {
                    let (u, v) = (rng.unit(), rng.unit());
                    Point::bpt(tri, [1.0 - u - v, u, v]).unwrap()
                })
                .collect();
            let mut pd = PlaneDelaunay::new(2);
            for (i, p) in points.iter().enumerate() {
                assert!(pd.insert(i as u32, *p));
            }
            check(&pd, &points);
        }
    }

    /// A square's corners and edge midpoints, then its centre (inside the
    /// circumcircles of the triangles through three midpoints).
    #[test]
    fn a_square_with_midpoints_then_its_centre() {
        let xy = [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.0],
            [1.0, 0.5],
            [0.5, 1.0],
            [0.0, 0.5],
            [0.5, 0.5],
        ];
        let points: Vec<Point> = xy
            .iter()
            .map(|q| Point::explicit([q[0], q[1], 0.0]).unwrap())
            .collect();
        let mut pd = PlaneDelaunay::new(2);
        for (i, p) in points.iter().enumerate() {
            assert!(pd.insert(i as u32, *p), "point {i}");
        }
        check(&pd, &points);
    }

    #[test]
    fn a_cocircular_grid() {
        // A 6 × 6 grid in the plane z = 0: every square's corners cocircular.
        let points: Vec<Point> = (0..36)
            .map(|i| Point::explicit([f64::from(i % 6), f64::from(i / 6), 0.0]).unwrap())
            .collect();
        let mut pd = PlaneDelaunay::new(2);
        for (i, p) in points.iter().enumerate() {
            assert!(pd.insert(i as u32, *p));
        }
        check(&pd, &points);
        assert_eq!(pd.triangles().len(), 50);
    }
}
