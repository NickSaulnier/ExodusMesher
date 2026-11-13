//! A constrained Delaunay triangulation of implicit points in a plane, in
//! its projection along an axis (Milestone 16, stage 2): every decision is
//! an exact [`implicit::orient2d_projected`] or
//! [`implicit::incircle_projected`].
//!
//! Points go in one at a time inside an explicit super-triangle (located by
//! a walk, then Lawson flips), then the segments (if not already edges, the
//! triangles they cross are removed and the two sides retriangulated, each
//! by choosing the vertex whose circle holds no other, after Anglada,
//! 1997). Segments must not cross each other or pass through a point:
//! callers split them first. Triangles are counter-clockwise in the
//! projection.

use std::collections::{HashMap, HashSet};

use crate::Point3;
use crate::implicit::{self, Point};
use crate::predicates::Sign;

/// Why a triangulation failed (a broken precondition: a bug upstream).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum CdtError {
    /// A segment crosses a constrained edge or passes through a point.
    Crossing,
    /// The walk lost its way, or a triangle went missing.
    Lost,
}

pub(crate) struct Cdt {
    /// Local points: the three super-triangle corners, then the input.
    pts: Vec<Point>,
    axis: usize,
    tris: Vec<[usize; 3]>,
    live: Vec<bool>,
    /// Directed edge to the live triangle holding it.
    edges: HashMap<(usize, usize), usize>,
    fixed: HashSet<(usize, usize)>,
    /// A recent triangle at each point.
    at: Vec<usize>,
    last: usize,
}

impl Cdt {
    /// A triangulation of `points` (all in one plane, whose projection along
    /// `axis` keeps area, and all within the box of the explicit `frame`)
    /// with no segments yet. Local index `i + 3` is `points[i]`.
    pub(crate) fn new(points: &[Point], axis: usize, frame: &[Point3]) -> Result<Cdt, CdtError> {
        let [i, j] = [(axis + 1) % 3, (axis + 2) % 3];
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for p in frame {
            for (m, k) in [i, j].into_iter().enumerate() {
                lo[m] = lo[m].min(p[k]);
                hi[m] = hi[m].max(p[k]);
            }
        }
        if !(lo.iter().chain(&hi).all(|x| x.is_finite())) {
            return Err(CdtError::Lost);
        }
        let size = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-300);
        let span = 16.0 * size + lo.iter().chain(&hi).map(|x| x.abs()).fold(0.0, f64::max);
        let centre = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])];
        let corner = |u: f64, v: f64| {
            let mut p: Point3 = [0.0; 3];
            p[i] = centre[0] + u * span;
            p[j] = centre[1] + v * span;
            Point::from(p)
        };
        let mut pts = vec![corner(-2.0, -2.0), corner(2.0, -2.0), corner(0.0, 2.0)];
        pts.extend_from_slice(points);
        let mut cdt = Cdt {
            at: vec![0; pts.len()],
            pts,
            axis,
            tris: Vec::new(),
            live: Vec::new(),
            edges: HashMap::new(),
            fixed: HashSet::new(),
            last: 0,
        };
        let first = if cdt.orient(0, 1, 2) == Sign::Positive {
            [0, 1, 2]
        } else {
            [0, 2, 1]
        };
        cdt.add(first);
        for v in 3..cdt.pts.len() {
            cdt.insert(v)?;
        }
        Ok(cdt)
    }

    fn orient(&self, a: usize, b: usize, c: usize) -> Sign {
        implicit::orient2d_projected(&self.pts[a], &self.pts[b], &self.pts[c], self.axis)
    }

    fn incircle(&self, a: usize, b: usize, c: usize, d: usize) -> Sign {
        implicit::incircle_projected(
            &self.pts[a],
            &self.pts[b],
            &self.pts[c],
            &self.pts[d],
            self.axis,
        )
    }

    fn add(&mut self, t: [usize; 3]) -> usize {
        let id = self.tris.len();
        self.tris.push(t);
        self.live.push(true);
        for k in 0..3 {
            self.edges.insert((t[k], t[(k + 1) % 3]), id);
            self.at[t[k]] = id;
        }
        self.last = id;
        id
    }

    fn kill(&mut self, id: usize) {
        self.live[id] = false;
        let t = self.tris[id];
        for k in 0..3 {
            let e = (t[k], t[(k + 1) % 3]);
            if self.edges.get(&e) == Some(&id) {
                self.edges.remove(&e);
            }
        }
    }

    /// The live triangle holding the directed edge `ab`.
    fn holding(&self, a: usize, b: usize) -> Option<usize> {
        self.edges.get(&(a, b)).copied()
    }

    fn is_fixed(&self, a: usize, b: usize) -> bool {
        self.fixed.contains(&(a.min(b), a.max(b)))
    }

    /// Where `p` lies: in triangle `t`, on its edge from corner `k + 1` to
    /// `k + 2` (`Some(k)`), or (`None` edge) strictly inside.
    fn locate(&self, p: usize) -> Result<(usize, Option<usize>), CdtError> {
        let mut t = if self.live[self.last] {
            self.last
        } else {
            self.live.iter().position(|&l| l).ok_or(CdtError::Lost)?
        };
        let limit = 4 * self.tris.len() + 64;
        for step in 0..limit {
            let tri = self.tris[t];
            let mut moved = false;
            let mut zero = None;
            for r in 0..3 {
                let k = (r + step) % 3;
                let (a, b) = (tri[(k + 1) % 3], tri[(k + 2) % 3]);
                match self.orient(a, b, p) {
                    Sign::Negative => {
                        t = self.holding(b, a).ok_or(CdtError::Lost)?;
                        moved = true;
                        break;
                    }
                    Sign::Zero => zero = Some(k),
                    Sign::Positive => {}
                }
            }
            if !moved {
                return Ok((t, zero));
            }
        }
        // A walk that cycles: look everywhere.
        for (t, &live) in self.live.iter().enumerate() {
            if !live {
                continue;
            }
            let tri = self.tris[t];
            let o: Vec<Sign> = (0..3)
                .map(|k| self.orient(tri[(k + 1) % 3], tri[(k + 2) % 3], p))
                .collect();
            if !o.contains(&Sign::Negative) {
                return Ok((t, o.iter().position(|&s| s == Sign::Zero)));
            }
        }
        Err(CdtError::Lost)
    }

    fn insert(&mut self, p: usize) -> Result<(), CdtError> {
        let (t, edge) = self.locate(p)?;
        let tri = self.tris[t];
        match edge {
            None => {
                self.kill(t);
                let [a, b, c] = tri;
                self.add([a, b, p]);
                self.add([b, c, p]);
                self.add([c, a, p]);
                for (x, y) in [(a, b), (b, c), (c, a)] {
                    self.legalize(p, x, y);
                }
            }
            Some(k) => {
                // On the edge from a to b, opposite c.
                let (a, b, c) = (tri[(k + 1) % 3], tri[(k + 2) % 3], tri[k]);
                if self.is_fixed(a, b) {
                    return Err(CdtError::Crossing);
                }
                let other = self.holding(b, a);
                self.kill(t);
                self.add([a, p, c]);
                self.add([p, b, c]);
                let mut outer = vec![(b, c), (c, a)];
                if let Some(u) = other {
                    let ut = self.tris[u];
                    let d = ut
                        .iter()
                        .copied()
                        .find(|&v| v != a && v != b)
                        .ok_or(CdtError::Lost)?;
                    self.kill(u);
                    self.add([b, p, d]);
                    self.add([p, a, d]);
                    outer.extend([(a, d), (d, b)]);
                }
                for (x, y) in outer {
                    self.legalize(p, x, y);
                }
            }
        }
        Ok(())
    }

    /// Restore the Delaunay property across edge `xy` of triangle `x y p`.
    fn legalize(&mut self, p: usize, x: usize, y: usize) {
        let mut stack = vec![(x, y)];
        while let Some((x, y)) = stack.pop() {
            if self.is_fixed(x, y) {
                continue;
            }
            let (Some(t), Some(u)) = (self.holding(x, y), self.holding(y, x)) else {
                continue;
            };
            if !self.tris[t].contains(&p) {
                continue;
            }
            let Some(d) = self.tris[u].iter().copied().find(|&v| v != x && v != y) else {
                continue;
            };
            if self.incircle(x, y, p, d) != Sign::Positive {
                continue;
            }
            self.kill(t);
            self.kill(u);
            self.add([x, d, p]);
            self.add([d, y, p]);
            stack.push((x, d));
            stack.push((d, y));
        }
    }

    /// The live triangles at `a`.
    fn around(&self, a: usize) -> Vec<usize> {
        let start = self.at[a];
        if self.live[start] && self.tris[start].contains(&a) {
            let mut out = vec![start];
            // Counter-clockwise from start, then clockwise if a boundary
            // stops the turn.
            let mut t = start;
            loop {
                let tri = self.tris[t];
                let k = tri.iter().position(|&v| v == a).unwrap_or(0);
                let prev = tri[(k + 2) % 3];
                match self.holding(a, prev) {
                    Some(u) if u == start => return out,
                    Some(u) => {
                        out.push(u);
                        t = u;
                    }
                    None => break,
                }
            }
            let mut t = start;
            loop {
                let tri = self.tris[t];
                let k = tri.iter().position(|&v| v == a).unwrap_or(0);
                let next = tri[(k + 1) % 3];
                match self.holding(next, a) {
                    Some(u) if out.contains(&u) => return out,
                    Some(u) => {
                        out.push(u);
                        t = u;
                    }
                    None => return out,
                }
            }
        }
        (0..self.tris.len())
            .filter(|&t| self.live[t] && self.tris[t].contains(&a))
            .collect()
    }

    /// Make local points `a` and `b` (indices into the input, without the
    /// super-triangle's three) an edge that no flip removes.
    pub(crate) fn constrain(&mut self, a: usize, b: usize) -> Result<(), CdtError> {
        let (a, b) = (a + 3, b + 3);
        if a == b {
            return Ok(());
        }
        if self.holding(a, b).is_some() || self.holding(b, a).is_some() {
            self.fixed.insert((a.min(b), a.max(b)));
            return Ok(());
        }
        // The triangle at a whose angle holds b.
        let mut first = None;
        for t in self.around(a) {
            let tri = self.tris[t];
            let k = tri.iter().position(|&v| v == a).unwrap_or(0);
            let (x, y) = (tri[(k + 1) % 3], tri[(k + 2) % 3]);
            if self.orient(a, x, b) == Sign::Positive && self.orient(a, b, y) == Sign::Positive {
                first = Some((t, x, y));
                break;
            }
        }
        let (mut t, mut x, mut y) = first.ok_or(CdtError::Crossing)?;
        // Walk across the triangles the segment crosses, x on its right
        // and y on its left.
        let (mut right, mut left) = (vec![x], vec![y]);
        let mut dead = vec![t];
        loop {
            if self.is_fixed(x, y) {
                return Err(CdtError::Crossing);
            }
            let u = self.holding(y, x).ok_or(CdtError::Lost)?;
            let z = self.tris[u]
                .iter()
                .copied()
                .find(|&v| v != x && v != y)
                .ok_or(CdtError::Lost)?;
            dead.push(u);
            t = u;
            if z == b {
                break;
            }
            match self.orient(a, b, z) {
                Sign::Negative => {
                    right.push(z);
                    x = z;
                }
                Sign::Positive => {
                    left.push(z);
                    y = z;
                }
                Sign::Zero => return Err(CdtError::Crossing),
            }
        }
        let _ = t;
        for &d in &dead {
            self.kill(d);
        }
        // Right of a→b: polygon b, a, right...; left: a, b, left reversed.
        let mut new = Vec::new();
        self.fill(b, a, &right, &mut new);
        left.reverse();
        self.fill(a, b, &left, &mut new);
        for t in new {
            self.add(t);
        }
        self.fixed.insert((a.min(b), a.max(b)));
        Ok(())
    }

    /// Triangulate the polygon `p, q, chain...` (counter-clockwise, the
    /// chain on the left of `p→q`) by Delaunay choice.
    fn fill(&self, p: usize, q: usize, chain: &[usize], out: &mut Vec<[usize; 3]>) {
        if chain.is_empty() {
            return;
        }
        let mut c = 0;
        for j in 1..chain.len() {
            if self.incircle(p, q, chain[c], chain[j]) == Sign::Positive {
                c = j;
            }
        }
        out.push([p, q, chain[c]]);
        self.fill(chain[c], q, &chain[..c], out);
        self.fill(p, chain[c], &chain[c + 1..], out);
    }

    /// The triangles, by input index (the super-triangle's left out).
    pub(crate) fn triangles(&self) -> Vec<[usize; 3]> {
        let mut out: Vec<[usize; 3]> = (0..self.tris.len())
            .filter(|&t| self.live[t] && self.tris[t].iter().all(|&v| v >= 3))
            .map(|t| self.tris[t].map(|v| v - 3))
            .collect();
        out.sort_unstable();
        out
    }
}
