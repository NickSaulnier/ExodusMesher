//! Worst-cell improvement by invariant-guarded local operations
//! (Milestone 15), in the style of the Wildmeshing Toolkit.
//!
//! A finished tetrahedral mesh is improved one poor cell at a time, cells
//! with a boundary face first, worst first, by local operations:
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
//! **Quality** is what a cell-centred finite-volume solver sees, against
//! targets ([`ImproveParams`]), each measure divided by its target so that 1
//! is on target and smaller is worse:
//!
//! - each tet's smallest dihedral angle;
//! - each face's non-orthogonality: the angle between its normal and the
//!   line joining the centroids of its two cells, or, on the boundary, the
//!   line from its cell's centroid to its own centroid;
//! - each interior face's skewness: the distance from its centroid to where
//!   the centroid line crosses its plane, over the centroid distance.
//!
//! A cell's quality is the worst of its own angle and its four faces',
//! capped at 1: a cell on target in every measure is not improved further.
//!
//! Each operation is accepted only if every tet it makes is positive
//! exactly (and in floating point, as a solver computes it) and the worst
//! quality of the region it changes (the new tets and every face of them,
//! against the old tets and theirs) strictly improves. No operation touches
//! a boundary face or moves a boundary vertex: flips and splits act only on
//! faces and edges inside the mesh, and moves and collapses only on
//! vertices off the boundary. So the boundary stays exactly the input's,
//! the volume is unchanged, and a mesh positive before stays positive
//! after. Deterministic: the poor cells go in order of quality and index,
//! and every choice is by a fixed order.
//!
//! **In parallel, with the serial result.** The poor cells of a pass are
//! taken in batches. Each cell of a batch is evaluated (its best operation
//! found) in parallel against the mesh as it is at the batch's start, and
//! the evaluation records what it read: the tets (their vertices, quality
//! and wall flag) and the vertex stars. The operations are then applied in
//! the serial order; one whose evaluation read a tet or a star that an
//! operation applied before it in the batch changed is evaluated again
//! first. So every operation is the one the serial pass would choose, and
//! the mesh is the same at any number of threads. A new vertex (a split's
//! or an insertion's) is a placeholder until its operation is applied, when
//! it gets the next index, as it would serially.

use std::cell::RefCell;

use mesher_geom::predicates::{insphere, orient3d};
use mesher_geom::vec3::{add, cross, dot, norm, scale, sub};
use mesher_geom::{Point3, Sign};
use rayon::prelude::*;

use crate::fxhash::FxHashMap;
use crate::quality::float_volume;
use crate::smooth::{best_position, min_dihedral};

/// An operation must raise the worst quality of its region by more than
/// this (a fraction of a target): a strict gain, so operations cannot undo
/// each other.
const GAIN: f64 = 1e-6;

/// Largest ring of tets an edge is removed from.
const MAX_RING: usize = 10;

/// Most tets an insertion replaces.
const MAX_CAVITY: usize = 128;

/// Poor cells evaluated together, per thread. Larger batches leave threads
/// less idle at the batch's end but go stale more often (on 894834 with 32
/// threads: 1, 2, 4, 8 and 16 a thread took 28, 24, 24, 27 and 33 s).
const BATCH_PER_THREAD: usize = 2;

/// An operation's new vertex (a split's or an insertion's) until the
/// operation is applied, when it gets the next index.
const NEW: u32 = u32::MAX;

/// What an evaluation read: tets (their vertices, quality and wall flag)
/// and vertex stars. Points need no record: a point is read only for a
/// vertex of a tet read, or one whose star was read, and moving a vertex
/// changes every tet of its star.
#[derive(Default)]
struct Reads {
    tets: Vec<u32>,
    stars: Vec<u32>,
}

thread_local! {
    /// The reads of the evaluation running on this thread, while one is
    /// recorded.
    static READS: RefCell<Option<Reads>> = const { RefCell::new(None) };
}

fn read_tet(t: u32) {
    READS.with(|r| {
        if let Some(r) = r.borrow_mut().as_mut() {
            r.tets.push(t);
        }
    });
}

fn read_star(v: u32) {
    READS.with(|r| {
        if let Some(r) = r.borrow_mut().as_mut() {
            r.stars.push(v);
        }
    });
}

/// `f`'s result and what it read.
fn recording<T>(f: impl FnOnce() -> T) -> (T, Reads) {
    READS.with(|r| *r.borrow_mut() = Some(Reads::default()));
    let out = f();
    let reads = READS.with(|r| r.borrow_mut().take()).unwrap_or_default();
    (out, reads)
}

/// The six edges of a tet, by local index.
const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// The face opposite each vertex of a positive tet, ordered so that the
/// face followed by that vertex is positive.
const INWARD: [[usize; 3]; 4] = [[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]];

/// How improvement works. A target at or beyond its measure's range (a
/// non-orthogonality of 90° or more, an infinite skewness, a dihedral angle
/// of 0°) turns that measure off.
#[derive(Debug, Clone, PartialEq)]
pub struct ImproveParams {
    /// Smallest dihedral angle of a tet, degrees.
    pub target_degrees: f64,
    /// Largest non-orthogonality of a face, degrees.
    pub max_non_orthogonality_degrees: f64,
    /// Largest skewness of an interior face.
    pub max_skewness: f64,
    /// Passes over the poor cells; stops early when nothing changes.
    pub passes: usize,
    /// Whether vertices may collapse, edges split, and vertices be inserted
    /// above boundary faces.
    pub collapses: bool,
    pub splits: bool,
    pub insertions: bool,
}

impl Default for ImproveParams {
    fn default() -> ImproveParams {
        ImproveParams {
            target_degrees: 15.0,
            max_non_orthogonality_degrees: 70.0,
            max_skewness: 2.0,
            passes: 8,
            collapses: true,
            splits: true,
            insertions: true,
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
    /// Vertices inserted above a boundary face.
    pub insertions: usize,
    /// Smallest dihedral angle before and after, degrees; -1 where a tet is
    /// not positive in floating point.
    pub min_before: f64,
    pub min_after: f64,
    /// The worst cell's quality before and after (1 on target; negative
    /// infinity for a tet not positive in floating point).
    pub worst_before: f64,
    pub worst_after: f64,
    /// Cells below target before and after.
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

/// Whether a tet is positive exactly and its floating-point volume positive.
fn usable(p: [Point3; 4]) -> bool {
    let volume = float_volume(p);
    orient3d(p[0], p[1], p[2], p[3]) == Sign::Positive && volume > 0.0
}

fn centroid(p: [Point3; 4]) -> Point3 {
    [0, 1, 2].map(|k| (p[0][k] + p[1][k] + p[2][k] + p[3][k]) / 4.0)
}

/// The measures against their targets.
#[derive(Debug, Clone, Copy)]
struct Targets {
    dihedral: f64,
    non_orthogonality: f64,
    skewness: f64,
}

impl Targets {
    fn new(params: &ImproveParams) -> Targets {
        Targets {
            dihedral: params.target_degrees,
            non_orthogonality: params.max_non_orthogonality_degrees,
            skewness: params.max_skewness,
        }
    }

    /// A usable tet's smallest dihedral angle against the target.
    fn tet(&self, p: [Point3; 4]) -> f64 {
        if self.dihedral > 0.0 {
            min_dihedral(p) / self.dihedral
        } else {
            f64::INFINITY
        }
    }

    /// Face `local` of usable tet `p` (centroid `c`) against the targets:
    /// its non-orthogonality, and its skewness if `other` (the centroid of
    /// the cell across) is given; on the boundary, `other` is `None`.
    fn face(&self, p: [Point3; 4], local: usize, c: Point3, other: Option<Point3>) -> f64 {
        let [a, b, d] = INWARD[local].map(|k| p[k]);
        let inward = cross(sub(b, a), sub(d, a));
        let fc = [0, 1, 2].map(|k| (a[k] + b[k] + d[k]) / 3.0);
        let line = sub(other.unwrap_or(fc), c);
        let (nn, length) = (norm(inward), norm(line));
        if !(nn > 0.0 && length > 0.0) {
            return f64::NEG_INFINITY;
        }
        let n = scale(inward, -1.0 / nn);
        let along = dot(line, n);
        let mut s = f64::INFINITY;
        if self.non_orthogonality < 90.0 {
            let angle = libm::acos((along / length).clamp(-1.0, 1.0)).to_degrees();
            s = (90.0 - angle) / (90.0 - self.non_orthogonality);
        }
        if other.is_some() && self.skewness.is_finite() {
            let skew = if along > 0.0 {
                let x = add(c, scale(line, dot(sub(fc, c), n) / along));
                norm(sub(fc, x)) / length
            } else {
                norm(sub(fc, c)) / length
            };
            s = s.min(self.skewness / skew);
        }
        if s.is_nan() { f64::NEG_INFINITY } else { s }
    }
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

/// What is across a face of a new tet.
#[derive(Debug, Clone, Copy)]
enum Across {
    /// Another new tet, by index.
    New(usize),
    /// A tet the operation keeps, by its centroid.
    Cell(Point3),
    /// The boundary.
    Wall,
}

/// The faces of an operation's new tets, each once: the tet (by index), its
/// local face, and what is across. The same for every position of a moved
/// vertex, since the tets across never contain it.
type Frame = Vec<(usize, usize, Across)>;

/// A local operation: tets out, tets in, and a vertex added (a split's).
struct Op {
    old: Vec<u32>,
    new: Vec<[u32; 4]>,
    /// The region's score after, as a `key`.
    worst: f64,
    /// Whether the region's worst rises (not only the other measure or the
    /// shortfall).
    raises: bool,
    point: Option<Point3>,
    /// A vertex move: the vertex and where; `old` is its star, `new` empty.
    moved: Option<(u32, Point3)>,
}

/// Which operation an `Op` is, for the statistics.
#[derive(Clone, Copy)]
enum Kind {
    Flip,
    EdgeRemoval,
    Move,
    Collapse,
    Split,
    Insertion,
}

/// A region's quality: by measure, each against its target and capped at
/// 1, the smallest dihedral angle and the worst face (non-orthogonality and
/// skewness), negative infinity where a tet is not usable; then the
/// region's shortfall, the sum over its tets and faces of the square of
/// how far each is below 1 (zero in a cached cell score).
type Score = [f64; 3];

const UNUSABLE: Score = [f64::NEG_INFINITY, f64::NEG_INFINITY, f64::INFINITY];

/// A measure's contribution to the shortfall.
fn shortfall(s: f64) -> f64 {
    let d = 1.0 - s.min(1.0);
    d * d
}

/// The worst of a score's measures.
fn low(s: Score) -> f64 {
    s[0].min(s[1])
}

/// The least a measure's worst may fall to when `before` is replaced:
/// its worst before, or halfway from the worst of both to the target, if
/// lower. So one measure is traded for the other only well above the
/// region's worst.
fn floor(before: Score, i: usize) -> f64 {
    before[i].min(0.5 * (1.0 + low(before)))
}

/// The better of a score's measures.
fn high(s: Score) -> f64 {
    s[0].max(s[1])
}

/// Whether `after` may replace `before`: neither measure falls below its
/// floor, the worst of both does not fall, and the worst strictly improves,
/// or the other measure does, or the shortfall falls. So a measure that
/// cannot improve (a face fixed by the boundary) does not stop the other,
/// and a cluster of poor cells, whose worst face no single operation can
/// reach, is worn down one operation at a time.
fn improves(before: Score, after: Score) -> bool {
    after[0] >= floor(before, 0)
        && after[1] >= floor(before, 1)
        && low(after) >= low(before)
        && (low(after) > low(before) + GAIN
            || high(after) > high(before) + GAIN
            || after[2] < before[2] * (1.0 - 1e-6) - GAIN)
}

/// A score as one number, for picking the best of several candidates: by
/// the worst measure, then the other and the shortfall.
fn key(s: Score) -> f64 {
    low(s) + 1e-3 * (high(s) - s[2] / (1.0 + s[2]))
}

struct Mesh {
    points: Vec<Point3>,
    tets: Vec<[u32; 4]>,
    /// Each tet's quality: its own angle and the worst of its faces.
    q: Vec<Score>,
    /// Whether a tet has a boundary face.
    wall: Vec<bool>,
    alive: Vec<bool>,
    /// The live tets around each vertex.
    star: Vec<Vec<u32>>,
    /// On a boundary face.
    fixed: Vec<bool>,
    targets: Targets,
    /// The batch that last changed each tet (removed it, or changed its
    /// quality or wall flag) and each vertex (moved it or changed its star).
    tet_batch: Vec<u32>,
    vertex_batch: Vec<u32>,
    batch: u32,
}

impl Mesh {
    /// Tet `t`'s vertices, recorded as read.
    fn tet(&self, t: u32) -> [u32; 4] {
        read_tet(t);
        self.tets[t as usize]
    }

    fn q_of(&self, t: u32) -> Score {
        read_tet(t);
        self.q[t as usize]
    }

    fn wall_of(&self, t: u32) -> bool {
        read_tet(t);
        self.wall[t as usize]
    }

    /// The live tets around vertex `v`, recorded as read.
    fn star_of(&self, v: u32) -> &[u32] {
        read_star(v);
        &self.star[v as usize]
    }

    /// Whether an evaluation that read `reads` may have read something an
    /// operation of this batch changed since.
    fn stale(&self, reads: &Reads) -> bool {
        reads
            .tets
            .iter()
            .any(|&t| self.tet_batch[t as usize] == self.batch)
            || reads
                .stars
                .iter()
                .any(|&v| self.vertex_batch[v as usize] == self.batch)
    }

    fn corners(&self, t: [u32; 4]) -> [Point3; 4] {
        t.map(|v| self.points[v as usize])
    }

    /// The live tet with face `f`, other than those in `skip`.
    fn across(&self, f: [u32; 3], skip: &[u32]) -> Option<u32> {
        self.star_of(f[0]).iter().copied().find(|&u| {
            let t = self.tet(u);
            !skip.contains(&u) && t.contains(&f[1]) && t.contains(&f[2])
        })
    }

    /// Tet `t`'s quality and whether it has a boundary face.
    fn cell(&self, t: u32) -> (Score, bool) {
        let _prof = crate::prof::scope("improve::cell");
        let tet = self.tets[t as usize];
        let p = self.corners(tet);
        let ok = usable(p);
        let c = centroid(p);
        let mut q = if ok {
            [self.targets.tet(p), f64::INFINITY, 0.0]
        } else {
            UNUSABLE
        };
        let mut wall = false;
        for (i, f) in INWARD.iter().enumerate() {
            let other = self.across(f.map(|k| tet[k]), &[t]);
            wall |= other.is_none();
            if ok {
                let other = other.map(|u| centroid(self.corners(self.tets[u as usize])));
                q[1] = q[1].min(self.targets.face(p, i, c, other));
            }
        }
        (
            [
                q[0].min(1.0),
                q[1].min(1.0),
                if ok { 0.0 } else { f64::INFINITY },
            ],
            wall,
        )
    }

    /// Recompute the quality of tets `ids` and of their neighbours.
    fn refresh(&mut self, ids: &[u32]) {
        let _prof = crate::prof::scope("improve::refresh");
        let mut all = Vec::with_capacity(ids.len() * 5);
        for &t in ids {
            all.push(t);
            let tet = self.tets[t as usize];
            for f in INWARD {
                all.extend(self.across(f.map(|k| tet[k]), &[t]));
            }
        }
        all.sort_unstable();
        all.dedup();
        for t in all {
            self.tet_batch[t as usize] = self.batch;
            (self.q[t as usize], self.wall[t as usize]) = self.cell(t);
        }
    }

    fn add(&mut self, t: [u32; 4]) -> u32 {
        let id = self.tets.len() as u32;
        self.q.push(UNUSABLE);
        self.wall.push(false);
        self.tet_batch.push(0);
        self.tets.push(t);
        self.alive.push(true);
        for v in t {
            self.star[v as usize].push(id);
        }
        id
    }

    fn remove(&mut self, id: u32) {
        self.alive[id as usize] = false;
        for v in self.tets[id as usize] {
            self.star[v as usize].retain(|&t| t != id);
        }
    }

    fn apply(&mut self, mut op: Op) {
        if let Some((v, x)) = op.moved {
            self.points[v as usize] = x;
            self.vertex_batch[v as usize] = self.batch;
            self.refresh(&op.old);
            return;
        }
        if let Some(x) = op.point {
            let m = self.points.len() as u32;
            self.points.push(x);
            self.star.push(Vec::new());
            self.fixed.push(false);
            self.vertex_batch.push(0);
            for v in op.new.iter_mut().flatten() {
                if *v == NEW {
                    *v = m;
                }
            }
        }
        for &t in &op.old {
            self.tet_batch[t as usize] = self.batch;
            for v in self.tets[t as usize] {
                self.vertex_batch[v as usize] = self.batch;
            }
            self.remove(t);
        }
        for &v in op.new.iter().flatten() {
            self.vertex_batch[v as usize] = self.batch;
        }
        let ids: Vec<u32> = op.new.into_iter().map(|t| self.add(t)).collect();
        self.refresh(&ids);
    }

    /// The worst of the cached scores of tets `ids` (no shortfall).
    fn worst(&self, ids: &[u32]) -> Score {
        ids.iter()
            .fold([f64::INFINITY, f64::INFINITY, 0.0], |w, &t| {
                let q = self.q_of(t);
                [w[0].min(q[0]), w[1].min(q[1]), 0.0]
            })
    }

    /// The score of tets `ids` and their faces as they are, with the
    /// shortfall.
    fn region(&self, ids: &[u32]) -> Option<Score> {
        let tets: Vec<[u32; 4]> = ids.iter().map(|&t| self.tet(t)).collect();
        let frame = self.frame(ids, &tets)?;
        Some(self.eval(&tets, &frame, None))
    }

    /// The faces of `new`, replacing `old`; `None` if they do not fit (a
    /// face on three new tets, or a face on one that is neither on a kept
    /// tet nor on the boundary).
    fn frame(&self, old: &[u32], new: &[[u32; 4]]) -> Option<Frame> {
        let _prof = crate::prof::scope("improve::frame");
        let mut keys: Vec<([u32; 3], usize, usize)> = Vec::with_capacity(4 * new.len());
        for (n, t) in new.iter().enumerate() {
            for (i, f) in INWARD.iter().enumerate() {
                let mut k = f.map(|j| t[j]);
                k.sort_unstable();
                keys.push((k, n, i));
            }
        }
        keys.sort_unstable();
        let mut old = old.to_vec();
        old.sort_unstable();
        let mut frame = Vec::with_capacity(keys.len());
        let mut i = 0;
        while i < keys.len() {
            let (k, n, l) = keys[i];
            let same = |j: usize| keys.get(j).is_some_and(|x| x.0 == k);
            if same(i + 1) {
                if same(i + 2) {
                    return None;
                }
                frame.push((n, l, Across::New(keys[i + 1].1)));
                i += 2;
                continue;
            }
            let kept = self.star_of(k[0]).iter().copied().find(|&u| {
                let t = self.tet(u);
                t.contains(&k[1]) && t.contains(&k[2]) && old.binary_search(&u).is_err()
            });
            let across = match kept {
                Some(u) => Across::Cell(centroid(self.corners(self.tet(u)))),
                None if k
                    .iter()
                    .all(|&v| self.fixed.get(v as usize).copied().unwrap_or(false)) =>
                {
                    Across::Wall
                }
                None => return None,
            };
            frame.push((n, l, across));
            i += 1;
        }
        Some(frame)
    }

    /// The worst quality of `new` and its faces (`frame`), with vertex
    /// `moved.0` at `moved.1`, capped at 1; `UNUSABLE` if a tet is not
    /// usable.
    fn eval(&self, new: &[[u32; 4]], frame: &Frame, moved: Option<(u32, Point3)>) -> Score {
        let _prof = crate::prof::scope("improve::eval");
        let pos = |v: u32| match moved {
            Some((m, x)) if m == v => x,
            _ => self.points[v as usize],
        };
        let mut worst = [f64::INFINITY, f64::INFINITY, 0.0];
        let mut corners = Vec::with_capacity(new.len());
        let mut centroids = Vec::with_capacity(new.len());
        for t in new {
            let p = t.map(pos);
            if !usable(p) {
                return UNUSABLE;
            }
            let s = self.targets.tet(p);
            worst[0] = worst[0].min(s);
            worst[2] += shortfall(s);
            corners.push(p);
            centroids.push(centroid(p));
        }
        for &(n, l, across) in frame {
            let other = match across {
                Across::New(j) => Some(centroids[j]),
                Across::Cell(c) => Some(c),
                Across::Wall => None,
            };
            let s = self.targets.face(corners[n], l, centroids[n], other);
            worst[1] = worst[1].min(s);
            worst[2] += shortfall(s);
        }
        [worst[0].min(1.0), worst[1].min(1.0), worst[2]]
    }

    /// `old` replaced by `new` if that improves the region's worst.
    fn offer(&self, old: Vec<u32>, new: Vec<[u32; 4]>) -> Option<Op> {
        self.offer_at(old, new, None)
    }

    /// `offer`, with a new vertex `added.0` at `added.1`.
    fn offer_at(
        &self,
        old: Vec<u32>,
        new: Vec<[u32; 4]>,
        added: Option<(u32, Point3)>,
    ) -> Option<Op> {
        let cached = self.worst(&old);
        let pos = |v: u32| match added {
            Some((m, x)) if m == v => x,
            _ => self.points[v as usize],
        };
        // The tets alone first: most candidates fail here, before their
        // faces are looked up.
        let least = floor(cached, 0).max(low(cached));
        for &t in &new {
            let p = t.map(pos);
            if !usable(p) || self.targets.tet(p).min(1.0) < least {
                return None;
            }
        }
        let frame = self.frame(&old, &new)?;
        let before = self.region(&old)?;
        let after = self.eval(&new, &frame, added);
        improves(before, after).then_some(Op {
            old,
            new,
            worst: key(after),
            raises: low(after) > low(before) + GAIN,
            point: added.map(|(_, x)| x),
            moved: None,
        })
    }

    /// The tet containing `p`, walking from tet `t`; `None` if the walk
    /// leaves the mesh or does not end.
    fn locate(&self, mut t: u32, p: Point3) -> Option<u32> {
        'walk: for _ in 0..4096 {
            let tet = self.tet(t);
            for f in INWARD {
                let f = f.map(|k| tet[k]);
                let [a, b, c] = f.map(|v| self.points[v as usize]);
                if orient3d(a, b, c, p) == Sign::Negative {
                    t = self.across(f, &[t])?;
                    continue 'walk;
                }
            }
            return Some(t);
        }
        None
    }

    /// Insert a vertex at `p`, which tet `t` or a tet near it contains: the
    /// tets whose circumspheres contain `p` (not across the boundary, at
    /// most `MAX_CAVITY`), less those that would hide a face of the cavity
    /// from `p`, replaced by the cone of the cavity's faces to `p`.
    fn insert(&self, t: u32, p: Point3) -> Option<Op> {
        let _prof = crate::prof::scope("improve::insert");
        let start = self.locate(t, p)?;
        let mut cavity = vec![start];
        let mut next = 0;
        while next < cavity.len() && cavity.len() < MAX_CAVITY {
            let c = cavity[next];
            next += 1;
            let tet = self.tet(c);
            for f in INWARD {
                let Some(u) = self.across(f.map(|k| tet[k]), &[c]) else {
                    continue;
                };
                if cavity.contains(&u) || cavity.len() >= MAX_CAVITY {
                    continue;
                }
                let [a, b, d, e] = self.corners(self.tet(u));
                if insphere(a, b, d, e, p) == Sign::Positive {
                    cavity.push(u);
                }
            }
        }
        // Shrink to a cavity every face of which `p` sees.
        let new = loop {
            let mut hidden = Vec::new();
            let mut new = Vec::new();
            for &c in &cavity {
                let tet = self.tet(c);
                for f in INWARD {
                    let f = f.map(|k| tet[k]);
                    if self.across(f, &[c]).is_some_and(|u| cavity.contains(&u)) {
                        continue;
                    }
                    let [a, b, d] = f.map(|v| self.points[v as usize]);
                    if orient3d(a, b, d, p) != Sign::Positive {
                        hidden.push(c);
                    }
                    new.push(f);
                }
            }
            if hidden.is_empty() {
                break new;
            }
            if hidden.contains(&start) {
                return None;
            }
            cavity.retain(|c| !hidden.contains(c));
        };
        let m = NEW;
        let new: Vec<[u32; 4]> = new.into_iter().map(|[a, b, c]| [a, b, c, m]).collect();
        // The cone must fill the cavity exactly: a cavity that is not
        // star-shaped from `p` would overlap itself.
        let before: f64 = cavity
            .iter()
            .map(|&c| float_volume(self.corners(self.tet(c))))
            .sum();
        let after: f64 = new
            .iter()
            .map(|t| float_volume(t.map(|v| if v == m { p } else { self.points[v as usize] })))
            .sum();
        if (after - before).abs() > 1e-9 * before.abs() || after.is_nan() {
            return None;
        }
        cavity.sort_unstable();
        self.offer_at(cavity, new, Some((m, p)))
    }

    /// Insertions for wall cell `t`: above each of its boundary faces, on
    /// the line through the face's centroid along its inward normal, at a
    /// few fractions of a regular tet's height on that face.
    fn insertions(&self, t: u32) -> Vec<Op> {
        let tet = self.tet(t);
        let mut ops = Vec::new();
        for f in INWARD {
            let f = f.map(|k| tet[k]);
            if self.across(f, &[t]).is_some() {
                continue;
            }
            let [a, b, c] = f.map(|v| self.points[v as usize]);
            let n = cross(sub(b, a), sub(c, a));
            let nn = norm(n);
            let edge = (norm(sub(b, a)) + norm(sub(c, b)) + norm(sub(a, c))) / 3.0;
            if !(nn > 0.0 && edge > 0.0) {
                continue;
            }
            let fc = [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0);
            // A regular tet's height is sqrt(2/3) of its edge.
            let height = 0.816_496_580_927_726 * edge;
            for s in [1.0, 0.7, 0.5, 0.35] {
                let p = add(fc, scale(n, s * height / nn));
                ops.extend(self.insert(t, p));
            }
        }
        ops
    }

    fn edge_tets(&self, a: u32, b: u32) -> Vec<u32> {
        self.star_of(a)
            .iter()
            .copied()
            .filter(|&t| self.tet(t).contains(&b))
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
            let tet = self.tet(t);
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

    /// What is across the ring faces of edge `ab`'s tets: for ring edge
    /// `(r[i], r[i + 1])`, the centroid of the kept tet across face
    /// `(a, r[i], r[i + 1])` and across `(b, r[i], r[i + 1])`, `None` for the
    /// boundary. `None` if a face is on neither.
    fn ring_faces(
        &self,
        a: u32,
        b: u32,
        old: &[u32],
        r: &[u32],
    ) -> Option<Vec<[Option<Point3>; 2]>> {
        let k = r.len();
        let side = |apex: u32, i: usize| -> Option<Option<Point3>> {
            let mut f = [apex, r[i], r[(i + 1) % k]];
            f.sort_unstable();
            match self.across(f, old) {
                Some(u) => Some(Some(centroid(self.corners(self.tet(u))))),
                None => f.iter().all(|&v| self.fixed[v as usize]).then_some(None),
            }
        };
        (0..k).map(|i| Some([side(a, i)?, side(b, i)?])).collect()
    }

    /// Remove interior edge `ab`: its ring of tets replaced by the
    /// triangulation of the ring polygon whose worst is best, each triangle
    /// coned to `a` and to `b`. The dynamic program scores each triangle by
    /// its two tets, the face between them and their faces on the ring; the
    /// triangulation it picks is then scored whole.
    fn edge_removal(&self, a: u32, b: u32) -> Option<Op> {
        let _prof = crate::prof::scope("improve::edge_removal");
        let (old, r) = self.ring(a, b)?;
        let k = r.len();
        if k > MAX_RING {
            return None;
        }
        let outer = self.ring_faces(a, b, &old, &r)?;
        let cones = |i: usize, m: usize, j: usize| [[r[i], r[j], r[m], a], [r[i], r[m], r[j], b]];
        // The ring edge of polygon side (p, q), if it is one.
        let edge = |p: usize, q: usize| {
            if q == p + 1 {
                Some(p)
            } else if p == 0 && q == k - 1 {
                Some(k - 1)
            } else {
                None
            }
        };
        let term = |i: usize, m: usize, j: usize| {
            let [ta, tb] = cones(i, m, j);
            let (pa, pb) = (self.corners(ta), self.corners(tb));
            if !(usable(pa) && usable(pb)) {
                return f64::NEG_INFINITY;
            }
            let (ca, cb) = (centroid(pa), centroid(pb));
            let t = &self.targets;
            let mut s = t.tet(pa).min(t.tet(pb)).min(t.face(pa, 3, ca, Some(cb)));
            // Cone faces (apex, r[p], r[q]) by local index in each cone.
            for (p, q, la, lb) in [(i, m, 1, 2), (m, j, 0, 0), (i, j, 2, 1)] {
                if let Some(e) = edge(p, q) {
                    s = s
                        .min(t.face(pa, la, ca, outer[e][0]))
                        .min(t.face(pb, lb, cb, outer[e][1]));
                }
            }
            s.min(1.0)
        };
        // best[i][j]: the best worst over triangulations of the polygon
        // r[i..=j]; choice[i][j]: the apex of its triangle on (r[i], r[j]).
        let mut best = vec![vec![f64::INFINITY; k]; k];
        let mut choice = vec![vec![usize::MAX; k]; k];
        #[allow(clippy::needless_range_loop)] // m indexes two tables and the ring
        for len in 2..k {
            for i in 0..k - len {
                let j = i + len;
                let mut top = f64::NEG_INFINITY;
                for m in i + 1..j {
                    let sub = best[i][m].min(best[m][j]);
                    if sub <= top {
                        continue;
                    }
                    let value = term(i, m, j).min(sub);
                    if value > top {
                        top = value;
                        choice[i][j] = m;
                    }
                }
                best[i][j] = top;
            }
        }
        if best[0][k - 1] < low(self.worst(&old)) {
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
        self.offer(old, new)
    }

    /// The 2-3 flip of the face of tet `t` opposite its vertex `i`.
    fn flip(&self, t: u32, i: usize) -> Option<Op> {
        let _prof = crate::prof::scope("improve::flip");
        let tet = self.tet(t);
        let d = tet[i];
        let [a, b, c] = INWARD[i].map(|k| tet[k]);
        let u = self.across([a, b, c], &[t])?;
        let e = self
            .tet(u)
            .into_iter()
            .find(|&x| x != a && x != b && x != c)?;
        if e == d {
            return None;
        }
        self.offer(vec![t, u], vec![[a, b, e, d], [b, c, e, d], [c, a, e, d]])
    }

    /// Collapse free vertex `v` into its neighbour `w`.
    fn collapse(&self, v: u32, w: u32) -> Option<Op> {
        let _prof = crate::prof::scope("improve::collapse");
        let old = self.star_of(v).to_vec();
        let new: Vec<[u32; 4]> = old
            .iter()
            .map(|&t| self.tet(t))
            .filter(|t| !t.contains(&w))
            .map(|t| t.map(|x| if x == v { w } else { x }))
            .collect();
        if new.is_empty() {
            return None;
        }
        self.offer(old, new)
    }

    /// Split interior edge `ab` at a new vertex, placed at the best of the
    /// candidate positions around the midpoint.
    fn split(&self, a: u32, b: u32) -> Option<Op> {
        let _prof = crate::prof::scope("improve::split");
        let (old, _) = self.ring(a, b)?;
        let before = self.region(&old)?;
        let m = NEW;
        let (pa, pb) = (self.points[a as usize], self.points[b as usize]);
        let mid = [0, 1, 2].map(|k| 0.5 * (pa[k] + pb[k]));
        let new: Vec<[u32; 4]> = old
            .iter()
            .flat_map(|&t| {
                let tet = self.tet(t);
                [
                    tet.map(|x| if x == a { m } else { x }),
                    tet.map(|x| if x == b { m } else { x }),
                ]
            })
            .collect();
        let frame = self.frame(&old, &new)?;
        let star: Vec<u32> = (0..new.len() as u32).collect();
        let pos = |v: u32| if v == m { mid } else { self.points[v as usize] };
        let at = |x: Point3| self.admissible(before, self.eval(&new, &frame, Some((m, x))));
        let (_, _, x) = best_position(&pos, &new, &star, m, &at)?;
        let after = self.eval(&new, &frame, Some((m, x)));
        improves(before, after).then_some(Op {
            old,
            new,
            worst: key(after),
            raises: low(after) > low(before) + GAIN,
            point: Some(x),
            moved: None,
        })
    }

    /// A vertex position's `key` for `best_position`: `None` if a tet is
    /// not usable or a measure falls below its floor (`floor`).
    fn admissible(&self, before: Score, after: Score) -> Option<f64> {
        (low(after) > f64::NEG_INFINITY
            && after[0] >= floor(before, 0)
            && after[1] >= floor(before, 1))
        .then_some(key(after))
    }

    /// Move free vertex `v` to its best candidate position, if better.
    fn move_op(&self, v: u32) -> Option<Op> {
        let _prof = crate::prof::scope("improve::move_op");
        if self.fixed[v as usize] {
            return None;
        }
        let star = self.star_of(v).to_vec();
        let tets: Vec<[u32; 4]> = star.iter().map(|&t| self.tet(t)).collect();
        let frame = self.frame(&star, &tets)?;
        let before = self.eval(&tets, &frame, None);
        let at = |x: Point3| self.admissible(before, self.eval(&tets, &frame, Some((v, x))));
        let local: Vec<u32> = (0..tets.len() as u32).collect();
        let pos = |w: u32| self.points[w as usize];
        let (_, _, x) = best_position(&pos, &tets, &local, v, &at)?;
        let after = self.eval(&tets, &frame, Some((v, x)));
        improves(before, after).then_some(Op {
            old: star,
            new: Vec::new(),
            worst: key(after),
            raises: low(after) > low(before) + GAIN,
            point: None,
            moved: Some((v, x)),
        })
    }

    /// The operation that improves poor tet `t`: flips first; if none
    /// raises the worst of its region, every kind of operation (flips,
    /// vertex moves, collapses, splits, insertions above its boundary
    /// faces), the best by `key`.
    fn evaluate(&self, t: u32, params: &ImproveParams) -> Option<(Op, Kind)> {
        let _prof = crate::prof::scope("improve::evaluate");
        let tet = self.tet(t);
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
        if !best.as_ref().is_some_and(|(b, _)| b.raises) {
            let mut verts = tet;
            verts.sort_unstable();
            for v in verts {
                offer(self.move_op(v), Kind::Move, &mut best);
            }
            for v in verts {
                if self.fixed[v as usize] || !params.collapses {
                    continue;
                }
                let mut neighbours: Vec<u32> = self
                    .star_of(v)
                    .iter()
                    .flat_map(|&s| self.tet(s))
                    .filter(|&w| w != v)
                    .collect();
                neighbours.sort_unstable();
                neighbours.dedup();
                for w in neighbours {
                    offer(self.collapse(v, w), Kind::Collapse, &mut best);
                }
            }
            if params.splits {
                for (i, j) in EDGES {
                    let op = self.split(tet[i], tet[j]);
                    offer(op, Kind::Split, &mut best);
                }
            }
            if params.insertions && self.wall_of(t) {
                for op in self.insertions(t) {
                    offer(Some(op), Kind::Insertion, &mut best);
                }
            }
        }
        best
    }

    /// Apply an operation `evaluate` chose, counting it.
    fn commit(&mut self, op: Op, kind: Kind, stats: &mut ImproveStats) {
        match kind {
            Kind::Flip => stats.flips += 1,
            Kind::EdgeRemoval => stats.edge_removals += 1,
            Kind::Move => stats.moves += 1,
            Kind::Collapse => stats.collapses += 1,
            Kind::Split => stats.splits += 1,
            Kind::Insertion => stats.insertions += 1,
        }
        self.apply(op);
    }

    /// The smallest dihedral angle (-1 if a tet is not usable), the worst
    /// quality, and the cells below target.
    fn summary(&self) -> (f64, f64, usize) {
        let (mut angle, mut worst, mut poor) = (f64::INFINITY, f64::INFINITY, 0);
        for (t, &q) in self.q.iter().enumerate() {
            let q = low(q);
            if self.alive[t] {
                let p = self.corners(self.tets[t]);
                angle = angle.min(if usable(p) { min_dihedral(p) } else { -1.0 });
                worst = worst.min(q);
                poor += usize::from(q < 1.0);
            }
        }
        let finite = |x: f64| if x == f64::INFINITY { 0.0 } else { x };
        (finite(angle), finite(worst), poor)
    }
}

/// Improve the poor cells of a positively oriented mesh whose boundary is
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
        wall: Vec::with_capacity(tets.len()),
        alive: Vec::with_capacity(tets.len()),
        star: vec![Vec::new(); points.len()],
        fixed: vec![false; points.len()],
        targets: Targets::new(params),
        tet_batch: Vec::with_capacity(tets.len()),
        vertex_batch: vec![0; points.len()],
        batch: 0,
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
    for t in 0..mesh.tets.len() as u32 {
        (mesh.q[t as usize], mesh.wall[t as usize]) = mesh.cell(t);
    }
    let mut stats = ImproveStats::default();
    (stats.min_before, stats.worst_before, stats.poor_before) = mesh.summary();
    for _ in 0..params.passes {
        stats.passes += 1;
        // Cells with a boundary face first, then by quality and index.
        let mut poor: Vec<u32> = (0..mesh.tets.len() as u32)
            .filter(|&t| mesh.alive[t as usize] && low(mesh.q[t as usize]) < 1.0)
            .collect();
        poor.sort_by(|&a, &b| {
            let (a, b) = (a as usize, b as usize);
            mesh.wall[b]
                .cmp(&mesh.wall[a])
                .then(low(mesh.q[a]).total_cmp(&low(mesh.q[b])))
                .then(a.cmp(&b))
        });
        let live = |mesh: &Mesh, t: u32| mesh.alive[t as usize] && low(mesh.q[t as usize]) < 1.0;
        let mut changed = 0;
        let threads = rayon::current_num_threads();
        if threads <= 1 {
            for t in poor {
                if live(&mesh, t) {
                    if let Some((op, kind)) = mesh.evaluate(t, params) {
                        mesh.commit(op, kind, &mut stats);
                        changed += 1;
                    }
                }
            }
        } else {
            for chunk in poor.chunks(BATCH_PER_THREAD * threads) {
                mesh.batch += 1;
                let m = &mesh;
                let speculated: Vec<_> = chunk
                    .par_iter()
                    .map(|&t| live(m, t).then(|| recording(|| m.evaluate(t, params))))
                    .collect();
                for (&t, speculated) in chunk.iter().zip(speculated) {
                    if !live(&mesh, t) {
                        continue;
                    }
                    let found = match speculated {
                        Some((found, reads)) if !mesh.stale(&reads) => found,
                        _ => {
                            let _prof = crate::prof::scope("improve::reevaluate");
                            mesh.evaluate(t, params)
                        }
                    };
                    if let Some((op, kind)) = found {
                        mesh.commit(op, kind, &mut stats);
                        changed += 1;
                    }
                }
            }
        }
        if changed == 0 {
            break;
        }
    }
    (stats.min_after, stats.worst_after, stats.poor_after) = mesh.summary();

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
            assert!(usable(tet.map(|v| p[v as usize])), "{tet:?}");
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
            assert!(s.worst_after >= s.worst_before, "{s:?}");
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

    /// The parallel pass gives the serial pass's mesh, bit for bit, at any
    /// number of threads; the meshes are dense in poor cells, so many
    /// speculative evaluations go stale and are redone.
    #[test]
    fn the_same_mesh_at_any_thread_count() {
        for seed in [4, 5] {
            let (p, t) = random_mesh(600, seed);
            let params = ImproveParams {
                target_degrees: 25.0,
                ..ImproveParams::default()
            };
            let at = |threads: usize| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap()
                    .install(|| improve(&p, &t, &params))
            };
            let serial = at(1);
            let s = serial.stats;
            assert!(s.splits + s.insertions + s.collapses > 0, "{s:?}");
            for threads in [2, 4, 16] {
                let r = at(threads);
                assert_eq!(r.points, serial.points, "{threads} threads");
                assert_eq!(r.tets, serial.tets, "{threads} threads");
                assert_eq!(r.map, serial.map, "{threads} threads");
                assert_eq!(r.stats, serial.stats, "{threads} threads");
            }
        }
    }

    /// Two tets mirrored through their shared face: the centroid line is
    /// along the face normal, so the face is orthogonal and unskewed.
    #[test]
    fn a_mirrored_face_is_on_target() {
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let t = Targets::new(&ImproveParams::default());
        // Face opposite vertex 0 is x + y + z = 1. This tet's centroid is
        // (1/4, 1/4, 1/4); its mirror through the face is (5/12, 5/12, 5/12).
        let c = centroid(p);
        let mirror = [5.0 / 12.0; 3];
        let s = t.face(p, 0, c, Some(mirror));
        assert!((s - 90.0 / 20.0).abs() < 1e-9, "{s}");
        // A boundary face: the line from the centroid to the face centroid
        // of the face x = 0 makes acos(0.25 / |(0.25, 1/12, 1/12)|).
        let wall = t.face(p, 1, c, None);
        let line: f64 = (0.0625f64 + 2.0 / 144.0).sqrt();
        let angle = libm::acos(0.25 / line).to_degrees();
        assert!((wall - (90.0 - angle) / 20.0).abs() < 1e-9, "{wall}");
        // Turned off, a measure never limits.
        let off = Targets::new(&ImproveParams {
            max_non_orthogonality_degrees: 90.0,
            max_skewness: f64::INFINITY,
            ..ImproveParams::default()
        });
        assert_eq!(off.face(p, 0, c, Some(mirror)), f64::INFINITY);
    }

    /// The face targets lower the worst non-orthogonality and skewness of
    /// a mesh with slivers, the boundary and volume kept, inserting points
    /// above boundary faces among other operations.
    #[test]
    fn face_targets_improve_the_worst_faces() {
        for seed in [4, 5] {
            let (p, t) = random_mesh(300, seed);
            let before = crate::fv_quality(&p, &t);
            let r = improve(&p, &t, &ImproveParams::default());
            check_valid(&r.points, &r.tets);
            assert_eq!(boundary(&r.points, &r.tets), boundary(&p, &t));
            assert!((volume(&r.points, &r.tets) - 1.0).abs() < 1e-12);
            let after = crate::fv_quality(&r.points, &r.tets);
            assert!(
                after.max_non_orthogonality_degrees < before.max_non_orthogonality_degrees,
                "{before:?} {after:?}"
            );
            assert!(
                after.max_skewness < before.max_skewness,
                "{before:?} {after:?}"
            );
            assert!(
                after.wall_max_non_orthogonality_degrees
                    <= before.wall_max_non_orthogonality_degrees,
                "{before:?} {after:?}"
            );
            assert!(r.stats.worst_after > r.stats.worst_before, "{:?}", r.stats);
            // Points go in above boundary faces, and the boundary stays.
            assert!(r.stats.insertions > 0, "{:?}", r.stats);
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
