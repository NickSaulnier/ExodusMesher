//! Worst-cell improvement by invariant-guarded local operations
//! (Milestone 15), in the style of the Wildmeshing Toolkit.
//!
//! A finished tetrahedral mesh is improved one poor tet at a time, worst
//! first, by local operations:
//!
//! - **Flips:** a 2-3 flip of an interior face, and edge removal of an
//!   interior edge (3-2 and 4-4 flips and their larger cousins), the ring of
//!   tets around the edge replaced by the best triangulation of its ring
//!   polygon (Shewchuk's maximin dynamic program).
//! - **Vertex moves:** a vertex on no boundary face moves to the best of the
//!   sliver smoothing's candidate positions (`smooth::best_position`).
//! - **Collapses:** a vertex on no boundary face merges into a neighbour.
//! - **Splits:** an interior edge splits at a new vertex, which then moves
//!   to the best candidate position.
//!
//! Each is accepted only if every tet it makes is positive exactly (and in
//! floating point, as a solver computes it) and the worst quality of the
//! region it replaces strictly improves. No operation touches a boundary
//! face or moves a boundary vertex: flips and splits act only on faces and
//! edges inside the mesh, and moves and collapses only on vertices off the
//! boundary. So the boundary stays exactly the input's, the volume is
//! unchanged, and a mesh positive before stays positive after.
//!
//! The target is the smallest dihedral angle for now. Deterministic: the
//! poor tets go in order of quality and index, and every choice is by a
//! fixed order.

use mesher_geom::predicates::orient3d;
use mesher_geom::{Point3, Sign};

use crate::fxhash::FxHashMap;
use crate::quality::float_volume;
use crate::smooth::{best_position, min_dihedral};

/// An operation must raise the worst quality of its region by more than
/// this (degrees): a strict gain, so operations cannot undo each other.
const GAIN: f64 = 1e-6;

/// Largest ring of tets an edge is removed from.
const MAX_RING: usize = 10;

/// The six edges of a tet, by local index.
const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// The face opposite each vertex of a positive tet, ordered so that the
/// face followed by that vertex is positive.
const INWARD: [[usize; 3]; 4] = [[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]];

/// How improvement works.
#[derive(Debug, Clone, PartialEq)]
pub struct ImproveParams {
    /// Tets whose smallest dihedral angle is below this are improved,
    /// degrees.
    pub target_degrees: f64,
    /// Passes over the poor tets; stops early when nothing changes.
    pub passes: usize,
    /// Whether vertices may collapse and edges split.
    pub collapses: bool,
    pub splits: bool,
}

impl Default for ImproveParams {
    fn default() -> ImproveParams {
        ImproveParams {
            target_degrees: 15.0,
            passes: 8,
            collapses: true,
            splits: true,
        }
    }
}

/// What improvement did.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ImproveStats {
    /// Passes run.
    pub passes: usize,
    /// 2-3 flips of an interior face.
    pub flips: usize,
    /// Interior edges removed (3-2, 4-4 and larger flips).
    pub edge_removals: usize,
    /// Vertex moves.
    pub moves: usize,
    /// Vertices collapsed into a neighbour.
    pub collapses: usize,
    /// Edges split at a new vertex.
    pub splits: usize,
    /// Smallest dihedral angle before and after, degrees; -1 where a tet is
    /// not positive in floating point.
    pub min_before: f64,
    pub min_after: f64,
    /// Tets below the target before and after.
    pub poor_before: usize,
    pub poor_after: usize,
}

/// An improved mesh.
#[derive(Debug, Clone)]
pub struct Improved {
    pub points: Vec<Point3>,
    /// Positively oriented.
    pub tets: Vec<[u32; 4]>,
    /// Each input vertex's index in `points`, `u32::MAX` for one collapsed
    /// away. Vertices added by splits come after the input's.
    pub map: Vec<u32>,
    pub stats: ImproveStats,
}

/// A tet's quality: its smallest dihedral angle in degrees, or -1 if it is
/// not positive exactly or its floating-point volume is not positive.
fn quality(p: [Point3; 4]) -> f64 {
    let volume = float_volume(p);
    if orient3d(p[0], p[1], p[2], p[3]) != Sign::Positive || volume <= 0.0 || volume.is_nan() {
        return -1.0;
    }
    min_dihedral(p)
}

/// `quality` for a candidate tet: `None` if it is not usable.
fn usable(p: [Point3; 4]) -> Option<f64> {
    let q = quality(p);
    (q >= 0.0).then_some(q)
}

/// Whether `order` is an even permutation of tet `t`'s vertices.
fn even(t: [u32; 4], order: [u32; 4]) -> bool {
    let pos = order.map(|v| t.iter().position(|&w| w == v).unwrap_or(0));
    let mut inversions = 0;
    for i in 0..4 {
        for j in i + 1..4 {
            if pos[i] > pos[j] {
                inversions += 1;
            }
        }
    }
    inversions % 2 == 0
}

/// A local operation: tets out, tets in, and a vertex added (a split's).
struct Op {
    old: Vec<u32>,
    new: Vec<[u32; 4]>,
    /// The worst quality of `new`.
    worst: f64,
    point: Option<Point3>,
}

/// Which operation an `Op` is, for the statistics.
#[derive(Clone, Copy)]
enum Kind {
    Flip,
    EdgeRemoval,
    Collapse,
    Split,
}

struct Mesh {
    points: Vec<Point3>,
    tets: Vec<[u32; 4]>,
    q: Vec<f64>,
    alive: Vec<bool>,
    /// The live tets around each vertex.
    star: Vec<Vec<u32>>,
    /// On a boundary face.
    fixed: Vec<bool>,
}

impl Mesh {
    fn corners(&self, t: [u32; 4]) -> [Point3; 4] {
        t.map(|v| self.points[v as usize])
    }

    fn add(&mut self, t: [u32; 4]) {
        let id = self.tets.len() as u32;
        self.q.push(quality(self.corners(t)));
        self.tets.push(t);
        self.alive.push(true);
        for v in t {
            self.star[v as usize].push(id);
        }
    }

    fn remove(&mut self, id: u32) {
        self.alive[id as usize] = false;
        for v in self.tets[id as usize] {
            self.star[v as usize].retain(|&t| t != id);
        }
    }

    fn apply(&mut self, op: Op) {
        if let Some(x) = op.point {
            self.points.push(x);
            self.star.push(Vec::new());
            self.fixed.push(false);
        }
        for &t in &op.old {
            self.remove(t);
        }
        for t in op.new {
            self.add(t);
        }
    }

    fn worst(&self, ids: &[u32]) -> f64 {
        ids.iter()
            .map(|&t| self.q[t as usize])
            .fold(f64::INFINITY, f64::min)
    }

    /// The worst quality of new tets, -1 as soon as one is not usable.
    fn new_worst(&self, tets: &[[u32; 4]]) -> f64 {
        let mut worst = f64::INFINITY;
        for &t in tets {
            let q = quality(self.corners(t));
            if q < 0.0 {
                return -1.0;
            }
            worst = worst.min(q);
        }
        worst
    }

    fn edge_tets(&self, a: u32, b: u32) -> Vec<u32> {
        self.star[a as usize]
            .iter()
            .copied()
            .filter(|&t| self.tets[t as usize].contains(&b))
            .collect()
    }

    /// The tets around interior edge `ab` and its ring of vertices, ordered
    /// so that `(a, b, r[i], r[i + 1])` is positive; `None` if the ring is
    /// not closed (a boundary edge) or not a simple cycle.
    fn ring(&self, a: u32, b: u32) -> Option<(Vec<u32>, Vec<u32>)> {
        let around = self.edge_tets(a, b);
        if around.len() < 3 {
            return None;
        }
        let mut pairs = Vec::with_capacity(around.len());
        for &t in &around {
            let tet = self.tets[t as usize];
            let mut others = tet.iter().copied().filter(|&x| x != a && x != b);
            let (Some(x), Some(y), None) = (others.next(), others.next(), others.next()) else {
                return None;
            };
            pairs.push(if even(tet, [a, b, x, y]) {
                (x, y)
            } else {
                (y, x)
            });
        }
        let first = pairs[0];
        let mut ring = vec![first.0];
        let mut cur = first.1;
        while cur != first.0 {
            if ring.len() >= pairs.len() {
                return None;
            }
            let mut next = pairs.iter().filter(|p| p.0 == cur);
            let (Some(&(_, to)), None) = (next.next(), next.next()) else {
                return None;
            };
            ring.push(cur);
            cur = to;
        }
        if ring.len() != pairs.len() {
            return None;
        }
        let mut sorted = ring.clone();
        sorted.sort_unstable();
        sorted.dedup();
        (sorted.len() == ring.len()).then_some((around, ring))
    }

    /// Remove interior edge `ab`: its ring of tets replaced by the
    /// triangulation of the ring polygon whose worst tet is best, each
    /// triangle coned to `a` and to `b`.
    fn edge_removal(&self, a: u32, b: u32) -> Option<Op> {
        let (old, r) = self.ring(a, b)?;
        let k = r.len();
        if k > MAX_RING {
            return None;
        }
        let before = self.worst(&old);
        let cones = |i: usize, m: usize, j: usize| [[r[i], r[j], r[m], a], [r[i], r[m], r[j], b]];
        // best[i][j]: the best worst quality over triangulations of the
        // polygon r[i..=j]; choice[i][j]: the apex of its triangle on
        // (r[i], r[j]).
        let mut best = vec![vec![f64::INFINITY; k]; k];
        let mut choice = vec![vec![usize::MAX; k]; k];
        #[allow(clippy::needless_range_loop)] // m indexes two tables and the ring
        for len in 2..k {
            for i in 0..k - len {
                let j = i + len;
                let mut top = -1.0;
                for m in i + 1..j {
                    let sub = best[i][m].min(best[m][j]);
                    if sub <= top {
                        continue;
                    }
                    let value = self.new_worst(&cones(i, m, j)).min(sub);
                    if value > top {
                        top = value;
                        choice[i][j] = m;
                    }
                }
                best[i][j] = top;
            }
        }
        let after = best[0][k - 1];
        if !(after >= 0.0 && after > before + GAIN) {
            return None;
        }
        let mut new = Vec::with_capacity(2 * (k - 2));
        let mut stack = vec![(0, k - 1)];
        while let Some((i, j)) = stack.pop() {
            if j < i + 2 {
                continue;
            }
            let m = choice[i][j];
            if m == usize::MAX {
                return None;
            }
            new.extend(cones(i, m, j));
            stack.push((i, m));
            stack.push((m, j));
        }
        Some(Op {
            old,
            new,
            worst: after,
            point: None,
        })
    }

    /// The 2-3 flip of the face of tet `t` opposite its vertex `i`.
    fn flip(&self, t: u32, i: usize) -> Option<Op> {
        let tet = self.tets[t as usize];
        let d = tet[i];
        let [a, b, c] = INWARD[i].map(|k| tet[k]);
        let mut others = self.star[a as usize].iter().copied().filter(|&u| {
            u != t && self.tets[u as usize].contains(&b) && self.tets[u as usize].contains(&c)
        });
        let (Some(u), None) = (others.next(), others.next()) else {
            return None;
        };
        let e = self.tets[u as usize]
            .into_iter()
            .find(|&x| x != a && x != b && x != c)?;
        if e == d {
            return None;
        }
        let new = vec![[a, b, e, d], [b, c, e, d], [c, a, e, d]];
        let before = self.q[t as usize].min(self.q[u as usize]);
        let after = self.new_worst(&new);
        (after >= 0.0 && after > before + GAIN).then_some(Op {
            old: vec![t, u],
            new,
            worst: after,
            point: None,
        })
    }

    /// Collapse free vertex `v` into its neighbour `w`.
    fn collapse(&self, v: u32, w: u32) -> Option<Op> {
        let old = self.star[v as usize].clone();
        let before = self.worst(&old);
        let new: Vec<[u32; 4]> = old
            .iter()
            .map(|&t| self.tets[t as usize])
            .filter(|t| !t.contains(&w))
            .map(|t| t.map(|x| if x == v { w } else { x }))
            .collect();
        if new.is_empty() {
            return None;
        }
        let after = self.new_worst(&new);
        (after >= 0.0 && after > before + GAIN).then_some(Op {
            old,
            new,
            worst: after,
            point: None,
        })
    }

    /// Split interior edge `ab` at a new vertex, placed at the best of the
    /// candidate positions around the midpoint.
    fn split(&mut self, a: u32, b: u32) -> Option<Op> {
        let (old, _) = self.ring(a, b)?;
        let before = self.worst(&old);
        let m = self.points.len() as u32;
        let (pa, pb) = (self.points[a as usize], self.points[b as usize]);
        let mid = [0, 1, 2].map(|k| 0.5 * (pa[k] + pb[k]));
        let new: Vec<[u32; 4]> = old
            .iter()
            .flat_map(|&t| {
                let tet = self.tets[t as usize];
                [
                    tet.map(|x| if x == a { m } else { x }),
                    tet.map(|x| if x == b { m } else { x }),
                ]
            })
            .collect();
        let star: Vec<u32> = (0..new.len() as u32).collect();
        self.points.push(mid);
        let found = best_position(&self.points, &new, &star, m, &usable);
        self.points.pop();
        let (_, after, x) = found?;
        (after >= 0.0 && after > before + GAIN).then_some(Op {
            old,
            new,
            worst: after,
            point: Some(x),
        })
    }

    /// Move free vertex `v` to its best candidate position, if better.
    fn try_move(&mut self, v: u32) -> bool {
        if self.fixed[v as usize] {
            return false;
        }
        let star = self.star[v as usize].clone();
        let Some((q0, q, x)) = best_position(&self.points, &self.tets, &star, v, &usable) else {
            return false;
        };
        if q <= q0 + GAIN || q.is_nan() {
            return false;
        }
        self.points[v as usize] = x;
        for t in star {
            self.q[t as usize] = quality(self.corners(self.tets[t as usize]));
        }
        true
    }

    /// Improve poor tet `t` by the first kind of operation that can: flips,
    /// then vertex moves, collapses and splits. Among operations of one kind
    /// the one whose new tets' worst is best wins.
    fn fix(&mut self, t: u32, params: &ImproveParams, stats: &mut ImproveStats) -> bool {
        let tet = self.tets[t as usize];
        let mut best: Option<(Op, Kind)> = None;
        let offer = |op: Option<Op>, kind: Kind, best: &mut Option<(Op, Kind)>| {
            if let Some(op) = op {
                if best.as_ref().is_none_or(|(b, _)| op.worst > b.worst) {
                    *best = Some((op, kind));
                }
            }
        };
        for (i, j) in EDGES {
            offer(
                self.edge_removal(tet[i], tet[j]),
                Kind::EdgeRemoval,
                &mut best,
            );
        }
        for i in 0..4 {
            offer(self.flip(t, i), Kind::Flip, &mut best);
        }
        if best.is_none() {
            let mut verts = tet;
            verts.sort_unstable();
            for v in verts {
                if self.try_move(v) {
                    stats.moves += 1;
                    return true;
                }
            }
            for v in verts {
                if self.fixed[v as usize] || !params.collapses {
                    continue;
                }
                let mut neighbours: Vec<u32> = self.star[v as usize]
                    .iter()
                    .flat_map(|&s| self.tets[s as usize])
                    .filter(|&w| w != v)
                    .collect();
                neighbours.sort_unstable();
                neighbours.dedup();
                for w in neighbours {
                    offer(self.collapse(v, w), Kind::Collapse, &mut best);
                }
            }
        }
        if best.is_none() && params.splits {
            for (i, j) in EDGES {
                let op = self.split(tet[i], tet[j]);
                offer(op, Kind::Split, &mut best);
            }
        }
        let Some((op, kind)) = best else {
            return false;
        };
        match kind {
            Kind::Flip => stats.flips += 1,
            Kind::EdgeRemoval => stats.edge_removals += 1,
            Kind::Collapse => stats.collapses += 1,
            Kind::Split => stats.splits += 1,
        }
        self.apply(op);
        true
    }

    fn summary(&self, target: f64) -> (f64, usize) {
        let mut worst = f64::INFINITY;
        let mut poor = 0;
        for (t, &q) in self.q.iter().enumerate() {
            if self.alive[t] {
                worst = worst.min(q);
                poor += usize::from(q < target);
            }
        }
        (if worst.is_finite() { worst } else { 0.0 }, poor)
    }
}

/// Improve the poor tets of a positively oriented mesh whose boundary is
/// the faces on one tet only. The boundary faces and their vertices never
/// change; see the module documentation.
#[must_use]
pub fn improve(points: &[Point3], tets: &[[u32; 4]], params: &ImproveParams) -> Improved {
    let _prof = crate::prof::scope("improve::improve");
    let unchanged = || Improved {
        points: points.to_vec(),
        tets: tets.to_vec(),
        map: (0..points.len() as u32).collect(),
        stats: ImproveStats::default(),
    };
    if tets.iter().flatten().any(|&v| v as usize >= points.len()) {
        return unchanged();
    }
    let mut mesh = Mesh {
        points: points.to_vec(),
        tets: Vec::with_capacity(tets.len()),
        q: Vec::with_capacity(tets.len()),
        alive: Vec::with_capacity(tets.len()),
        star: vec![Vec::new(); points.len()],
        fixed: vec![false; points.len()],
    };
    for &t in tets {
        mesh.add(t);
    }
    let face = |t: [u32; 4], i: usize| {
        let mut f = INWARD[i].map(|k| t[k]);
        f.sort_unstable();
        f
    };
    let mut uses: FxHashMap<[u32; 3], u32> = FxHashMap::default();
    for &t in tets {
        for i in 0..4 {
            *uses.entry(face(t, i)).or_default() += 1;
        }
    }
    for &t in tets {
        for i in 0..4 {
            let f = face(t, i);
            if uses.get(&f).copied() != Some(2) {
                for v in f {
                    mesh.fixed[v as usize] = true;
                }
            }
        }
    }
    let mut stats = ImproveStats::default();
    (stats.min_before, stats.poor_before) = mesh.summary(params.target_degrees);
    for _ in 0..params.passes {
        stats.passes += 1;
        let mut poor: Vec<u32> = (0..mesh.tets.len() as u32)
            .filter(|&t| mesh.alive[t as usize] && mesh.q[t as usize] < params.target_degrees)
            .collect();
        poor.sort_by(|&a, &b| {
            mesh.q[a as usize]
                .total_cmp(&mesh.q[b as usize])
                .then(a.cmp(&b))
        });
        let mut changed = 0;
        for t in poor {
            if mesh.alive[t as usize]
                && mesh.q[t as usize] < params.target_degrees
                && mesh.fix(t, params, &mut stats)
            {
                changed += 1;
            }
        }
        if changed == 0 {
            break;
        }
    }
    (stats.min_after, stats.poor_after) = mesh.summary(params.target_degrees);

    // Compact to the vertices in use, in their order.
    let mut used = vec![false; mesh.points.len()];
    for (t, tet) in mesh.tets.iter().enumerate() {
        if mesh.alive[t] {
            for &v in tet {
                used[v as usize] = true;
            }
        }
    }
    let mut map = vec![u32::MAX; mesh.points.len()];
    let mut out_points = Vec::with_capacity(mesh.points.len());
    for (v, &p) in mesh.points.iter().enumerate() {
        if used[v] {
            map[v] = out_points.len() as u32;
            out_points.push(p);
        }
    }
    let out_tets = mesh
        .tets
        .iter()
        .enumerate()
        .filter(|&(t, _)| mesh.alive[t])
        .map(|(_, tet)| tet.map(|v| map[v as usize]))
        .collect();
    map.truncate(points.len());
    Improved {
        points: out_points,
        tets: out_tets,
        map,
        stats,
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::panic)]
    use super::*;
    use std::collections::BTreeMap;

    fn volume(p: &[Point3], t: &[[u32; 4]]) -> f64 {
        t.iter()
            .map(|t| float_volume(t.map(|v| p[v as usize])))
            .sum()
    }

    /// Boundary faces, oriented outward, as coordinates (indices change).
    fn boundary(p: &[Point3], t: &[[u32; 4]]) -> Vec<[[u64; 3]; 3]> {
        let mut faces: BTreeMap<[u32; 3], (usize, [u32; 3])> = BTreeMap::new();
        for tet in t {
            for f in INWARD {
                // Outward: reverse the inward order.
                let o = [tet[f[0]], tet[f[2]], tet[f[1]]];
                let mut key = o;
                key.sort_unstable();
                faces.entry(key).or_insert((0, o)).0 += 1;
            }
        }
        let mut out: Vec<[[u64; 3]; 3]> = faces
            .values()
            .filter(|(n, _)| *n == 1)
            .map(|(_, o)| {
                // Rotate so the smallest coordinate triple leads.
                let c = o.map(|v| p[v as usize].map(f64::to_bits));
                let k = (0..3).min_by_key(|&k| c[k]).unwrap();
                [c[k], c[(k + 1) % 3], c[(k + 2) % 3]]
            })
            .collect();
        out.sort_unstable();
        out
    }

    /// Every tet positive, every face on at most two tets, and none twice.
    fn check_valid(p: &[Point3], t: &[[u32; 4]]) {
        let mut faces: BTreeMap<[u32; 3], usize> = BTreeMap::new();
        let mut seen = std::collections::BTreeSet::new();
        for tet in t {
            assert!(quality(tet.map(|v| p[v as usize])) >= 0.0, "{tet:?}");
            let mut key = *tet;
            key.sort_unstable();
            assert!(seen.insert(key), "tet twice: {tet:?}");
            for f in INWARD {
                let mut k = f.map(|i| tet[i]);
                k.sort_unstable();
                *faces.entry(k).or_default() += 1;
            }
        }
        assert!(faces.values().all(|&n| n <= 2));
    }

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*seed >> 11) as f64 / (1u64 << 53) as f64
    }

    /// The Delaunay mesh of random points in a cube, with its corners.
    fn random_mesh(n: usize, seed: u64) -> (Vec<Point3>, Vec<[u32; 4]>) {
        let mut s = seed;
        let mut p: Vec<Point3> = (0..8)
            .map(|i| {
                [
                    f64::from(i & 1),
                    f64::from((i >> 1) & 1),
                    f64::from((i >> 2) & 1),
                ]
            })
            .collect();
        for _ in 0..n {
            p.push([lcg(&mut s), lcg(&mut s), lcg(&mut s)]);
        }
        let t = crate::delaunay(&p).unwrap().tets;
        (p, t)
    }

    #[test]
    fn inward_faces_are_positive_with_their_vertex() {
        for (i, f) in INWARD.iter().enumerate() {
            assert!(even(
                [0, 1, 2, 3],
                [f[0] as u32, f[1] as u32, f[2] as u32, i as u32]
            ));
        }
    }

    /// Three tets around a long edge: removing it (a 3-2 flip) leaves two
    /// positive tets.
    #[test]
    fn a_three_ring_edge_is_removed() {
        let mut p = vec![[0.0, 0.0, 1.5], [0.0, 0.0, -1.5]];
        for k in 0..3 {
            let a = 2.0 * std::f64::consts::PI * f64::from(k) / 3.0;
            p.push([libm::cos(a), libm::sin(a), 0.0]);
        }
        let mut t = Vec::new();
        for k in 0..3u32 {
            let mut tet = [0, 1, 2 + k, 2 + (k + 1) % 3];
            if orient3d(p[0], p[1], p[tet[2] as usize], p[tet[3] as usize]) != Sign::Positive {
                tet.swap(2, 3);
            }
            t.push(tet);
        }
        let before = boundary(&p, &t);
        let r = improve(
            &p,
            &t,
            &ImproveParams {
                target_degrees: 60.0,
                passes: 1,
                ..ImproveParams::default()
            },
        );
        assert_eq!(r.stats.edge_removals, 1, "{:?}", r.stats);
        assert_eq!(r.tets.len(), 2);
        check_valid(&r.points, &r.tets);
        assert_eq!(boundary(&r.points, &r.tets), before);
    }

    /// Two flat tets on a shared face, the segment between their apexes
    /// through it: the 2-3 flip makes three better ones.
    #[test]
    fn a_face_is_flipped() {
        let p = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.9, 0.0],
            [0.5, 0.3, 0.1],
            [0.5, 0.3, -0.1],
        ];
        let mut t = vec![[0, 1, 2, 3], [0, 2, 1, 4]];
        for tet in &mut t {
            let q = tet.map(|v| p[v as usize]);
            if orient3d(q[0], q[1], q[2], q[3]) != Sign::Positive {
                tet.swap(0, 1);
            }
        }
        let before = boundary(&p, &t);
        let r = improve(
            &p,
            &t,
            &ImproveParams {
                target_degrees: 60.0,
                passes: 1,
                ..ImproveParams::default()
            },
        );
        assert_eq!(r.stats.flips, 1, "{:?}", r.stats);
        assert_eq!(r.tets.len(), 3);
        check_valid(&r.points, &r.tets);
        assert_eq!(boundary(&r.points, &r.tets), before);
        assert!(r.stats.min_after > r.stats.min_before);
    }

    /// Random points make slivers. Improvement keeps the mesh valid, its
    /// boundary and volume, never worsens the worst tet, and is
    /// deterministic.
    #[test]
    fn improvement_keeps_the_invariants() {
        for seed in [1, 2, 3] {
            let (p, t) = random_mesh(300, seed);
            let params = ImproveParams {
                target_degrees: 25.0,
                ..ImproveParams::default()
            };
            let r = improve(&p, &t, &params);
            check_valid(&r.points, &r.tets);
            assert_eq!(boundary(&r.points, &r.tets), boundary(&p, &t));
            assert!((volume(&r.points, &r.tets) - 1.0).abs() < 1e-12);
            let s = r.stats;
            assert!(s.min_after >= s.min_before, "{s:?}");
            assert!(s.poor_after < s.poor_before, "{s:?}");
            assert!(s.flips + s.edge_removals > 0 && s.moves > 0, "{s:?}");
            // Corners and every boundary vertex stay where they were.
            for (v, &to) in r.map.iter().enumerate().take(8) {
                assert_eq!(r.points[to as usize], p[v]);
            }
            let again = improve(&p, &t, &params);
            assert_eq!(again.points, r.points);
            assert_eq!(again.tets, r.tets);
        }
    }

    #[test]
    fn bad_indices_leave_the_mesh_unchanged() {
        let p = vec![[0.0; 3]; 3];
        let t = vec![[0, 1, 2, 7]];
        let r = improve(&p, &t, &ImproveParams::default());
        assert_eq!(r.tets, t);
        assert_eq!(r.stats, ImproveStats::default());
    }
}
