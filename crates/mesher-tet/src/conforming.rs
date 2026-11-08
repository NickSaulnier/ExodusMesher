//! Boundary-exact tetrahedral meshing of a closed surface: a constrained
//! Delaunay tetrahedralization of the surface, refined.
//!
//! *Facets* are the edge-connected, same-group triangles that are coplanar
//! (exactly, or to rounding); only edges between facets are *subsegments*.
//! Every vertex carries the facets it lies on, and points sharing a facet are
//! coplanar by construction (see `Delaunay::set_planes`): no tet is made flat
//! within a facet, and a point never "sees" a face on its own facet, whatever
//! rounding says. All other decisions use exact predicates.
//!
//! 1. *Delaunay* of the surface vertices, with those plane labels.
//! 2. *Boundary recovery.* Every subsegment and surface triangle (*subfacet*)
//!    must appear in the tetrahedralization, and what appears is constrained:
//!    no later insertion may delete it. A missing subsegment is split (on
//!    concentric shells around input vertices). A missing subfacet is
//!    recovered, in order of preference, by adopting the tetrahedralization's
//!    own triangulation of its facet, by swapping a facet diagonal to the
//!    mesh's, by flipping away the edges crossing it, or by a Steiner point.
//!    Boundary Steiner points are inserted by 2D Bowyer–Watson on the facet
//!    with the tets on both sides seeding the 3D cavity, so the new subfacets
//!    appear with the point. Steiner points lie on the input surface, so the
//!    geometry is unchanged.
//! 3. *Constrained Delaunay refinement.* A tet inside the region whose
//!    radius-edge ratio exceeds the bound (and whose edges are not below the
//!    size floor), or whose circumradius exceeds the size bound, gets its
//!    circumcentre inserted, unless that point encroaches a subsegment
//!    (inside its diametral sphere) or a subfacet (inside its diametral
//!    sphere and projecting into it), which is split instead; or unless a
//!    subfacet lies between the tet and the point, which is split. A boundary
//!    split the constraints refuse is forced (face constraints ignored), and
//!    the subfacets it deletes are recovered at once.
//! 4. *Sliver removal.* Flips that make the worst tet fatter (edge removal
//!    and 2–3 flips, never across a constraint), a couple of rounds of
//!    circumcentre insertion, and flips again.
//! 5. *Region marking.* A flood from the hull that toggles inside and outside
//!    at every subfacet; the kept tets' boundary must be exactly the
//!    subfacets, or the result is refused.
//!
//! A point budget bounds every stage; exceeding it, or failing to recover the
//! boundary, is an error the caller turns into a fallback.

use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::Arc;

use mesher_geom::predicates::{incircle, insphere, orient2d, orient3d};
use mesher_geom::vec3::{add, cross, distance_squared, dot, norm, scale, sub};
use mesher_geom::{FastWindingNumber, Point3, Sign, TriangleSoup};

use crate::delaunay::{Delaunay, DelaunayError, Inserted, Walk};
use crate::fxhash::FxHashMap;
use crate::quality::dihedral_angles;
use crate::sizing::SizeField;
use crate::smooth::{SmoothParams, smooth_slivers};

#[derive(Debug, Clone, PartialEq)]
pub struct ConformingParams {
    /// Tets with a larger circumradius are refined. With a size field, the
    /// largest size anywhere (the field is capped at it).
    pub max_circumradius: f64,
    /// Target size by position: a tet whose circumradius exceeds the size at
    /// its centroid is refined, and the point budget follows the field.
    /// `None`: `max_circumradius` everywhere.
    pub size: Option<Arc<SizeField>>,
    /// Tets with a larger circumradius-to-shortest-edge ratio are refined.
    pub radius_edge_bound: f64,
    /// Refuse to grow beyond this many points; `0` sizes the budget from the
    /// surface (20 × (volume / size³ + area / size²), at least 10,000).
    pub max_points: usize,
    /// Refuse to take more than this many kernel steps (point insertions,
    /// including refused ones, refinement and repair iterations, recovery
    /// passes); `0` sizes it from the point budget (64 steps a point). A
    /// deterministic bound on work, so a cascade or a loop ends as a
    /// recorded failure instead of a worker timeout.
    pub max_work: u64,
    /// Refinement and sliver removal never split a subsegment or subfacet
    /// smaller than this (recovery always may). Without a floor, input
    /// segments meeting at small angles encroach on each other forever.
    pub min_boundary_size: f64,
    /// Tets with a smaller dihedral angle are slivers.
    pub sliver_degrees: f64,
    pub sliver_rounds: usize,
    /// Passes of vertex smoothing on the finished mesh (free vertices of
    /// tets below 15° move to raise their worst angle; see
    /// `smooth::smooth_slivers`); `0` turns it off.
    pub smooth_passes: usize,
    pub seed: u64,
    /// Group of each input triangle (e.g. its OBJ group). Rounding-flat
    /// slivers between two input triangles are only flattened when the two
    /// share a group, so group boundaries (side sets) stay exact. Empty: one
    /// group.
    pub groups: Vec<u32>,
}

impl ConformingParams {
    #[must_use]
    pub fn with_size(max_circumradius: f64) -> ConformingParams {
        ConformingParams {
            max_circumradius,
            size: None,
            radius_edge_bound: 2.0,
            max_points: 0,
            max_work: 0,
            min_boundary_size: max_circumradius / 16.0,
            sliver_degrees: 10.0,
            sliver_rounds: 2,
            smooth_passes: 4,
            seed: 0x5eed,
            groups: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConformingError {
    /// The surface cannot be meshed conformally as given (coincident
    /// vertices, invalid triangles).
    Surface(String),
    Delaunay(DelaunayError),
    /// The point budget was exceeded.
    Budget {
        points: usize,
        limit: usize,
    },
    /// The work budget was exceeded.
    Work {
        steps: u64,
        limit: u64,
    },
    /// Some boundary triangles are missing from the result.
    Recovery {
        missing_faces: usize,
    },
}

impl std::fmt::Display for ConformingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConformingError::Surface(e) => write!(f, "surface not usable: {e}"),
            ConformingError::Delaunay(e) => write!(f, "{e}"),
            ConformingError::Budget { points, limit } => {
                write!(f, "{points} points exceed the budget of {limit}")
            }
            ConformingError::Work { steps, limit } => {
                write!(f, "{steps} kernel steps exceed the work budget of {limit}")
            }
            ConformingError::Recovery { missing_faces } => {
                write!(
                    f,
                    "{missing_faces} boundary triangles could not be recovered"
                )
            }
        }
    }
}

/// A conforming tetrahedral mesh of the region a surface bounds.
#[derive(Debug, Clone, PartialEq)]
pub struct ConformingMesh {
    pub points: Vec<Point3>,
    /// Positively oriented.
    pub tets: Vec<[u32; 4]>,
    /// Boundary triangles, outward, each with the input triangle it lies on.
    pub boundary: Vec<([u32; 3], u32)>,
    pub boundary_steiner: usize,
    pub interior_steiner: usize,
    /// Tets still below the sliver angle.
    pub slivers_left: usize,
    /// Vertex moves made by smoothing.
    pub smoothed: usize,
    /// Rounding-flat slivers dropped by flipping a planar quad's diagonal.
    pub flattened: usize,
}

/// Local faces of a positive tet, outward.
const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

#[derive(Debug, Clone, Copy)]
struct Sub {
    v: [u32; 3],
    origin: u32,
    alive: bool,
}

/// A uniform hash grid of item ids by bounding box (probed, never iterated).
#[derive(Debug)]
struct Grid {
    cell: f64,
    cells: HashMap<[i64; 3], Vec<u32>>,
    /// Items spanning too many cells to register (the diametral sphere of a
    /// sliver triangle); every probe sees them.
    big: Vec<u32>,
}

/// Items spanning more cells than this along an axis go to `Grid::big`.
const BIG_SPAN: i64 = 16;

impl Grid {
    fn key(&self, p: Point3) -> [i64; 3] {
        p.map(|c| (c / self.cell).floor() as i64)
    }

    fn add(&mut self, id: u32, lo: Point3, hi: Point3) {
        let (a, b) = (self.key(lo), self.key(hi));
        if (0..3).any(|k| b[k] - a[k] > BIG_SPAN) {
            self.big.push(id);
            return;
        }
        let span = |k: usize| (b[k] - a[k]).max(0);
        for i in 0..=span(0) {
            for j in 0..=span(1) {
                for k in 0..=span(2) {
                    self.cells
                        .entry([a[0] + i, a[1] + j, a[2] + k])
                        .or_default()
                        .push(id);
                }
            }
        }
    }

    /// Items whose box may contain `p` (plus every big item).
    fn at(&self, p: Point3) -> impl Iterator<Item = u32> + '_ {
        self.cells
            .get(&self.key(p))
            .map_or(&[][..], Vec::as_slice)
            .iter()
            .chain(&self.big)
            .copied()
    }

    /// Forget a big item (others stay registered; probes filter them).
    fn remove_big(&mut self, id: u32) {
        if let Some(i) = self.big.iter().position(|&x| x == id) {
            self.big.swap_remove(i);
        }
    }
}

enum Split {
    Edge(u32, u32),
    Sub(u32),
}

/// Where a point went.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Placed {
    New(u32),
    /// It equals this existing vertex; nothing changed.
    Existing(u32),
    /// Inserting it would delete a constraint; nothing changed.
    Blocked,
    /// Subfacets in the way were rearranged; compute the split again.
    Retry,
}

/// Is a tet flat to rounding? Rounding moves points by a few ulps of their
/// coordinates, which far from the origin is much more than of the tet's size.
fn rounding_flat(p: [Point3; 4]) -> bool {
    let longest = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        .iter()
        .map(|&(i, j)| distance_squared(p[i], p[j]))
        .fold(0.0, f64::max)
        .sqrt();
    let volume = dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))).abs() / 6.0;
    let magnitude = p.iter().flatten().fold(0.0f64, |m, c| m.max(c.abs()));
    let flat = (1e-12 * longest).max(64.0 * f64::EPSILON * magnitude);
    volume <= flat * longest * longest
}

struct Builder<'a> {
    dt: Delaunay,
    subs: Vec<Sub>,
    edge_subs: HashMap<(u32, u32), Vec<u32>>,
    sub_grid: Grid,
    /// The live subfacet on each face (sorted vertices).
    face_sub: HashMap<[u32; 3], u32>,
    winding: FastWindingNumber<'a>,
    params: &'a ConformingParams,
    inputs: usize,
    boundary_steiner: usize,
    tet_queue: VecDeque<(u32, [u32; 4])>,
    /// Refinement attempts per tet, so no tet is retried forever.
    attempts: HashMap<[u32; 4], u8>,
    /// An orthonormal frame in the plane of each input triangle: origin and
    /// two axes. Flips use these fixed 2D coordinates, so that rounding
    /// off the plane cannot make the in-circle tests disagree.
    frames: Vec<[Point3; 3]>,
    /// The facet of each input triangle: the smallest triangle index among the
    /// edge-connected, coplanar (exactly or to rounding), same-group triangles
    /// facing the same way. A facet may be triangulated any way, and only
    /// edges between facets are subsegments.
    facet: Vec<u32>,
    /// Whether splits respect `min_boundary_size` (off during recovery).
    floored: bool,
    /// New subfacets that are not (yet) in the tetrahedralization, because a
    /// forced insertion deleted them; repaired before refining further.
    missing: Vec<u32>,
    /// Whether subfacets in the tetrahedralization are constrained (off
    /// while recovering).
    face_constraints: bool,
    /// Kernel steps taken, against `params.max_work`.
    work: u64,
    /// Whether each tet (by sorted vertices) is inside; see `tet_inside`.
    inside_cache: FxHashMap<[u32; 4], bool>,
}

/// Missing subfacets past which a growing repair cascade is left to a full
/// recovery pass.
const CASCADE: usize = 256;
/// Attempts allowed per tet before it is left alone.
const MAX_ATTEMPTS: u8 = 8;
/// Flips one legalization may make before giving up (a few per insertion are typical).
const MAX_FLIPS: usize = 10_000;

/// The facet of each triangle (see `Builder::facet`), by union-find over
/// manifold edges whose two triangles are coplanar (exactly, or to rounding:
/// their four points make a tet flat to rounding), face the same way and
/// share a group.
fn facets(points: &[Point3], triangles: &[[u32; 3]], groups: &[u32]) -> Vec<u32> {
    let _prof = crate::prof::scope("conforming::facets");
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            let up = parent[parent[x as usize] as usize];
            parent[x as usize] = up;
            x = up;
        }
        x
    }
    let mut parent: Vec<u32> = (0..triangles.len() as u32).collect();
    let mut by_edge: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
    for (i, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            by_edge
                .entry(key2(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(i as u32);
        }
    }
    let mut pairs: Vec<(u32, u32)> = by_edge
        .values()
        .filter_map(|l| match l[..] {
            [a, b] => Some((a.min(b), a.max(b))),
            _ => None,
        })
        .collect();
    pairs.sort_unstable();
    let normal = |t: [u32; 3]| {
        let [a, b, c] = t.map(|v| points[v as usize]);
        cross(sub(b, a), sub(c, a))
    };
    for (i, j) in pairs {
        let (ti, tj) = (triangles[i as usize], triangles[j as usize]);
        if groups.get(i as usize) != groups.get(j as usize) {
            continue;
        }
        let Some(&d) = tj.iter().find(|v| !ti.contains(v)) else {
            continue;
        };
        let [a, b, c] = ti.map(|v| points[v as usize]);
        let pd = points[d as usize];
        let coplanar = orient3d(a, b, c, pd) == Sign::Zero || rounding_flat([a, b, c, pd]);
        if !coplanar || dot(normal(ti), normal(tj)) <= 0.0 {
            continue;
        }
        let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
        let (lo, hi) = (ri.min(rj), ri.max(rj));
        parent[hi as usize] = lo;
    }
    (0..triangles.len() as u32)
        .map(|i| find(&mut parent, i))
        .collect()
}

/// Origin and in-plane orthonormal axes of a triangle (any frame if degenerate).
fn plane_frame(a: Point3, b: Point3, c: Point3) -> [Point3; 3] {
    let n = cross(sub(b, a), sub(c, a));
    let e = sub(b, a);
    // An axis-aligned plane drops a coordinate: exact, so ties stay ties.
    if let Some(k) = (0..3).find(|&k| n[(k + 1) % 3] == 0.0 && n[(k + 2) % 3] == 0.0 && n[k] != 0.0)
    {
        let mut e1 = [0.0; 3];
        let mut e2 = [0.0; 3];
        e1[(k + 1) % 3] = 1.0;
        e2[(k + 2) % 3] = 1.0;
        return [[0.0; 3], e1, e2];
    }
    let (ln, le) = (norm(n), norm(e));
    if ln == 0.0 || le == 0.0 {
        return [a, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    }
    let e1 = scale(e, 1.0 / le);
    let e2 = cross(scale(n, 1.0 / ln), e1);
    [a, e1, e2]
}

fn key2(u: u32, v: u32) -> (u32, u32) {
    (u.min(v), u.max(v))
}

fn key3(mut f: [u32; 3]) -> [u32; 3] {
    f.sort_unstable();
    f
}

/// Circumcentre and circumradius of a triangle in 3D (`None` if degenerate).
fn tri_circle(a: Point3, b: Point3, c: Point3) -> Option<(Point3, f64)> {
    let (ab, ac) = (sub(b, a), sub(c, a));
    let n = cross(ab, ac);
    let nn = dot(n, n);
    if nn <= 0.0 || !nn.is_finite() {
        return None;
    }
    let offset = scale(
        add(
            scale(cross(n, ab), dot(ac, ac)),
            scale(cross(ac, n), dot(ab, ab)),
        ),
        0.5 / nn,
    );
    Some((add(a, offset), norm(offset)))
}

/// Circumcentre and circumradius of a tet (`None` if flat).
fn tet_sphere(p: [Point3; 4]) -> Option<(Point3, f64)> {
    let (b, c, d) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
    let det = dot(b, cross(c, d));
    if det == 0.0 || !det.is_finite() {
        return None;
    }
    let num = add(
        add(scale(cross(c, d), dot(b, b)), scale(cross(d, b), dot(c, c))),
        scale(cross(b, c), dot(d, d)),
    );
    let offset = scale(num, 0.5 / det);
    Some((add(p[0], offset), norm(offset)))
}

/// Is `p` inside triangle `abc` (coplanar, by barycentric coordinates)?
fn inside_triangle(p: Point3, a: Point3, b: Point3, c: Point3) -> bool {
    let n = cross(sub(b, a), sub(c, a));
    let nn = dot(n, n);
    if nn <= 0.0 {
        return false;
    }
    let w = |u: Point3, v: Point3| dot(cross(sub(v, u), sub(p, u)), n) / nn;
    w(a, b) > 1e-9 && w(b, c) > 1e-9 && w(c, a) > 1e-9
}

/// Does `p` project (orthogonally) onto the closed triangle `abc`? Only then
/// does it encroach the subfacet: the diametral sphere of a thin subfacet is
/// huge, and constraints keep the subfacet whatever lies in it.
fn projects_into(p: Point3, a: Point3, b: Point3, c: Point3) -> bool {
    let n = cross(sub(b, a), sub(c, a));
    let w = |u: Point3, v: Point3| dot(cross(sub(v, u), sub(p, u)), n);
    w(a, b) >= 0.0 && w(b, c) >= 0.0 && w(c, a) >= 0.0
}

impl Builder<'_> {
    fn pos(&self, v: u32) -> Point3 {
        self.dt.points()[v as usize]
    }

    fn register_sub(&mut self, id: u32) {
        let s = self.subs[id as usize];
        let [a, b, c] = s.v.map(|v| self.pos(v));
        let longest = [(a, b), (b, c), (c, a)]
            .iter()
            .map(|&(x, y)| norm(sub(x, y)))
            .fold(0.0, f64::max);
        let r = tri_circle(a, b, c).map_or(longest, |(_, r)| r.max(longest));
        let centre = scale(add(add(a, b), c), 1.0 / 3.0);
        let lo = centre.map(|x| x - 2.0 * r);
        let hi = centre.map(|x| x + 2.0 * r);
        self.sub_grid.add(id, lo, hi);
        for (u, v) in [(s.v[0], s.v[1]), (s.v[1], s.v[2]), (s.v[2], s.v[0])] {
            self.edge_subs.entry(key2(u, v)).or_default().push(id);
        }
        self.face_sub.insert(key3(s.v), id);
    }

    fn add_sub(&mut self, v: [u32; 3], origin: u32) {
        self.subs.push(Sub {
            v,
            origin,
            alive: true,
        });
        let id = (self.subs.len() - 1) as u32;
        self.register_sub(id);
    }

    fn kill_sub(&mut self, id: u32) {
        let s = self.subs[id as usize];
        self.subs[id as usize].alive = false;
        self.sub_grid.remove_big(id);
        if self.face_sub.get(&key3(s.v)) == Some(&id) {
            self.face_sub.remove(&key3(s.v));
            self.dt.unconstrain_face(s.v);
        }
        for (u, v) in [(s.v[0], s.v[1]), (s.v[1], s.v[2]), (s.v[2], s.v[0])] {
            if let Some(list) = self.edge_subs.get_mut(&key2(u, v)) {
                list.retain(|&x| x != id);
            }
        }
    }

    /// Insert a point into the tetrahedralization (see `Delaunay::insert_with`)
    /// and queue the new tets.
    fn insert_point(
        &mut self,
        p: Point3,
        seeds: &[u32],
        keep: &[u32],
        force: bool,
    ) -> Result<Placed, ConformingError> {
        let _prof = crate::prof::scope("conforming::insert_point");
        self.step()?;
        let limit = self.params.max_points;
        if self.dt.points().len() >= limit {
            return Err(ConformingError::Budget {
                points: self.dt.points().len(),
                limit,
            });
        }
        let r = self.dt.insert_with(p, seeds, keep, force);
        match r {
            Ok((i, Inserted::New)) => {
                // Subfacets a forced insertion deleted are recovered next.
                for f in self.dt.lost().to_vec() {
                    if let Some(&s) = self.face_sub.get(&f) {
                        self.missing.push(s);
                    }
                }
                for slot in self.dt.last_created() {
                    if let Some(v) = self.dt.tet(slot) {
                        self.tet_queue.push_back((slot, v));
                    }
                }
                Ok(Placed::New(i))
            }
            Ok((_, Inserted::Duplicate(v))) => Ok(Placed::Existing(v)),
            Err(DelaunayError::Blocked { .. }) => Ok(Placed::Blocked),
            Err(e) => Err(ConformingError::Delaunay(e)),
        }
    }

    /// Remove a tet flat to rounding that touches the boundary (typically
    /// lying in a facet whose points are coplanar only to rounding) by
    /// flipping away one of its edges that is not on the boundary. Such a
    /// tet blocks splits next to it: rounding decides on which side of its
    /// faces a new point on the facet lies.
    fn unflatten(&mut self, slot: u32) -> bool {
        let _prof = crate::prof::scope("conforming::unflatten");
        let Some(t) = self.dt.tet(slot) else {
            return false;
        };
        let touches = FACES
            .iter()
            .any(|f| self.face_sub.contains_key(&key3(f.map(|k| t[k]))));
        if !touches || !(self.dt.on_one_plane(&t) || rounding_flat(t.map(|v| self.pos(v)))) {
            return false;
        }
        // The edge between the two vertices not on a subfacet face first
        // (the other diagonal of a quad planar to rounding).
        let mut edges: Vec<(u32, u32)> = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            .into_iter()
            .map(|(a, b)| key2(t[a], t[b]))
            .filter(|e| self.edge_subs.get(e).is_none_or(Vec::is_empty))
            .collect();
        edges.sort_by_key(|&(a, b)| {
            let on = |x: u32| {
                FACES.iter().any(|f| {
                    let g = f.map(|k| t[k]);
                    g.contains(&x) && self.face_sub.contains_key(&key3(g))
                })
            };
            (on(a) as u8 + on(b) as u8, a, b)
        });
        edges.into_iter().any(|(a, b)| self.dt.remove_edge(a, b))
    }

    /// Remove the tets flat to rounding at the boundary (see `unflatten`),
    /// a few passes.
    fn unflatten_all(&mut self) {
        let _prof = crate::prof::scope("conforming::unflatten_all");
        for _ in 0..4 {
            let mut changed = false;
            for slot in self.dt.real_slots() {
                changed |= self.unflatten(slot);
            }
            if !changed {
                break;
            }
        }
    }

    /// Constrain a subfacet, and its subsegments, that are in the
    /// tetrahedralization; then no insertion can delete them.
    fn protect(&mut self, s: u32) {
        let v = self.subs[s as usize].v;
        if !self.face_constraints {
            // Recovery: subfacets are recovered the conforming way.
        } else if self.dt.has_face(v) {
            self.dt.constrain_face(v);
        } else {
            self.missing.push(s);
        }
        for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            if self.is_segment(a, b) && self.dt.has_edge(a, b) {
                self.dt.constrain_edge(a, b);
            }
        }
    }

    /// The facet's subfacets whose circumcircle strictly contains `p`,
    /// reached from `start` without crossing a subsegment (2D Bowyer–Watson
    /// on the facet, in its exact frame), shrunk until `p` sees every
    /// boundary edge.
    fn facet_cavity(&self, start: u32, p: Point3, split: Option<(u32, u32)>) -> Vec<u32> {
        let _prof = crate::prof::scope("conforming::facet_cavity");
        let f = self.facet[self.subs[start as usize].origin as usize];
        let [o, e1, e2] = self.frames[f as usize];
        let flat = |q: Point3| {
            let d = sub(q, o);
            [dot(d, e1), dot(d, e2)]
        };
        let pp = flat(p);
        let tri = |s: u32| self.subs[s as usize].v.map(|x| flat(self.pos(x)));
        let mut cavity = vec![start];
        let mut k = 0;
        while k < cavity.len() {
            let s = cavity[k];
            k += 1;
            let v = self.subs[s as usize].v;
            for (x, y) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
                if self.is_segment(x, y) {
                    continue;
                }
                for &t in self
                    .edge_subs
                    .get(&key2(x, y))
                    .map_or(&[][..], Vec::as_slice)
                {
                    if t == s || cavity.contains(&t) {
                        continue;
                    }
                    let [a, b, c] = tri(t);
                    let side = orient2d(a, b, c);
                    if side != Sign::Zero && incircle(a, b, c, pp) == side {
                        cavity.push(t);
                    }
                }
            }
        }
        // Keep the fan valid: every boundary edge strictly in front of p.
        loop {
            let unseen = cavity.iter().copied().skip(1).find(|&s| {
                let v = self.subs[s as usize].v;
                let [a, b, c] = tri(s);
                let side = orient2d(a, b, c);
                (0..3).any(|r| {
                    let (x, y) = (v[r], v[(r + 1) % 3]);
                    let shared = cavity.iter().any(|&t| {
                        t != s
                            && self.subs[t as usize].v.contains(&x)
                            && self.subs[t as usize].v.contains(&y)
                    });
                    let q = [a, b, c];
                    !shared
                        && Some(key2(x, y)) != split
                        && orient2d(q[r], q[(r + 1) % 3], pp) != side
                })
            });
            match unseen {
                Some(s) => cavity.retain(|&t| t != s),
                None => return cavity,
            }
        }
    }

    /// Put `p` on the boundary: the subfacets in `cavity` (star-shaped from
    /// `p` in their facets, holding the `primary` ones that hold `p`), and
    /// the subsegment `split` that `p` lies on, are replaced by fans from
    /// `p`. The tets on both sides of the replaced subfacets and around
    /// `split` must all go, or an old face or edge would survive as a flat
    /// tet beside the new ones. If the constraints refuse that, the primary
    /// subfacets alone are tried, and then a forced insertion, whose deleted
    /// subfacets are recovered next.
    fn split_boundary(
        &mut self,
        p: Point3,
        cavity: &[u32],
        primary: &[u32],
        split: Option<(u32, u32)>,
    ) -> Result<Placed, ConformingError> {
        let _prof = crate::prof::scope("conforming::split_boundary");
        // A tet flat to rounding next to what is split would block it.
        let mut near: Vec<u32> = primary
            .iter()
            .flat_map(|&s| self.dt.face_tets(self.subs[s as usize].v))
            .collect();
        if let Some((u, v)) = split {
            near.extend(self.dt.edge_tets(u, v));
        }
        let mut changed = false;
        for t in near {
            changed |= self.unflatten(t);
        }
        if changed {
            return Ok(Placed::Retry);
        }
        for &s in cavity {
            self.dt.unconstrain_face(self.subs[s as usize].v);
        }
        if let Some((u, v)) = split {
            self.dt.unconstrain_edge(u, v);
        }
        let attempts: [(&[u32], bool); 3] = [(cavity, false), (primary, false), (primary, true)];
        let mut outcome = Ok(Placed::Blocked);
        let mut used: &[u32] = cavity;
        for (subs, force) in attempts {
            let (mut seeds, mut keep) = (Vec::new(), Vec::new());
            let mut around: Vec<u32> = subs
                .iter()
                .flat_map(|&s| self.dt.face_tets(self.subs[s as usize].v))
                .collect();
            if let Some((u, v)) = split {
                around.extend(self.dt.edge_tets(u, v));
            }
            for t in around {
                if self.dt.is_ghost_slot(t) {
                    seeds.push(t);
                } else {
                    keep.push(t);
                }
            }
            used = subs;
            let mut planes: Vec<u32> = primary
                .iter()
                .map(|&s| self.facet[self.subs[s as usize].origin as usize])
                .collect();
            planes.sort_unstable();
            planes.dedup();
            self.dt.set_planes(self.dt.points().len() as u32, planes);
            outcome = self.insert_point(p, &seeds, &keep, force);
            if !matches!(outcome, Ok(Placed::Blocked)) {
                break;
            }
        }
        let Ok(Placed::New(i)) = outcome else {
            for &s in cavity {
                self.protect(s);
            }
            self.missing.retain(|s| !cavity.contains(s));
            if let Some((u, v)) = split {
                if self.dt.has_edge(u, v) {
                    self.dt.constrain_edge(u, v);
                }
            }
            return outcome;
        };
        self.boundary_steiner += 1;
        let used = used.to_vec();
        // Subfacets of the larger cavity left in place keep their constraint.
        for &s in cavity {
            if !used.contains(&s) {
                self.protect(s);
            }
        }
        let mut fan = Vec::new();
        for &s in &used {
            let Sub { v, origin, .. } = self.subs[s as usize];
            for r in 0..3 {
                let (x, y) = (v[r], v[(r + 1) % 3]);
                let inner = used.iter().any(|&t| {
                    t != s
                        && self.facet[self.subs[t as usize].origin as usize]
                            == self.facet[origin as usize]
                        && self.subs[t as usize].v.contains(&x)
                        && self.subs[t as usize].v.contains(&y)
                });
                if !inner && Some(key2(x, y)) != split {
                    fan.push(([x, y, i], origin));
                }
            }
        }
        for &s in &used {
            self.kill_sub(s);
        }
        let first = self.subs.len() as u32;
        for (v, origin) in fan {
            self.add_sub(v, origin);
        }
        for s in first..self.subs.len() as u32 {
            self.protect(s);
        }
        Ok(Placed::New(i))
    }

    /// Where to split subsegment `uv`: its midpoint, unless exactly one end
    /// is an input vertex. Then it is split on a concentric shell around that
    /// vertex (a distance of size × 2^k, near the midpoint), so subsegments
    /// meeting there at a small angle are split at the same radii and stop
    /// encroaching each other (Ruppert's concentric shells).
    fn split_point(&self, u: u32, v: u32) -> Point3 {
        let inputs = self.inputs as u32;
        let (apex, other) = match (u < inputs, v < inputs) {
            (true, false) => (u, v),
            (false, true) => (v, u),
            _ => return scale(add(self.pos(u), self.pos(v)), 0.5),
        };
        let (a, b) = (self.pos(apex), self.pos(other));
        let length = norm(sub(b, a));
        let unit = self.params.max_circumradius;
        let d = unit * (length / (2.0 * unit)).log2().round().exp2();
        let t = if d.is_finite() && d > 0.0 {
            (d / length).clamp(0.25, 0.75)
        } else {
            0.5
        };
        add(a, scale(sub(b, a), t))
    }

    /// Split a subsegment. `false` when the floor forbids it.
    fn split_edge(&mut self, u: u32, v: u32) -> Result<bool, ConformingError> {
        let _prof = crate::prof::scope("conforming::split_edge");
        if self.floored
            && distance_squared(self.pos(u), self.pos(v)) < self.params.min_boundary_size.powi(2)
        {
            return Ok(false);
        }
        let mut m = self.split_point(u, v);
        // Other points along the edge, for when the constraints refuse one.
        let (pu, pv) = (self.pos(u), self.pos(v));
        let mut others = [0.5, 0.3, 0.7]
            .into_iter()
            .map(|t| add(pu, scale(sub(pv, pu), t)))
            .filter(|&q| q != m)
            .collect::<Vec<_>>()
            .into_iter();
        let mut attempts = 0;
        let (i, on_edge) = loop {
            let on_edge = self.edge_subs.get(&key2(u, v)).cloned().unwrap_or_default();
            if on_edge.is_empty() {
                // Swapped away: the edge is no longer on the boundary.
                return Ok(true);
            }
            let mut cavity: Vec<u32> = Vec::new();
            for &s in &on_edge {
                for t in self.facet_cavity(s, m, Some(key2(u, v))) {
                    if !cavity.contains(&t) {
                        cavity.push(t);
                    }
                }
            }
            match self.split_boundary(m, &cavity, &on_edge, Some(key2(u, v)))? {
                Placed::New(_) => return Ok(true),
                Placed::Blocked => match others.next() {
                    Some(q) => m = q,
                    None => {
                        return Ok(false);
                    }
                },
                Placed::Retry if attempts < 4 => attempts += 1,
                Placed::Retry => return Ok(false),
                Placed::Existing(i) => break (i, on_edge),
            }
        };
        // The split point is an existing vertex, so it lies on the edge: the
        // apex of a subfacet flat to rounding, or a vertex of a degenerate
        // input triangle along the same line. The subfacets on the edge split
        // at it (one whose apex it is collapses); their halves are recovered
        // like any missing subfacet. Only an end of the edge itself (an edge
        // too short to split) is a failure.
        if i == u || i == v {
            return Err(ConformingError::Recovery { missing_faces: 1 });
        }
        self.dt.unconstrain_edge(u, v);
        let first = self.subs.len() as u32;
        for s in on_edge {
            let Sub { v: t, origin, .. } = self.subs[s as usize];
            let r = (0..3)
                .find(|&k| key2(t[k], t[(k + 1) % 3]) == key2(u, v))
                .unwrap_or(0);
            let (x, y, w) = (t[r], t[(r + 1) % 3], t[(r + 2) % 3]);
            self.kill_sub(s);
            if w != i {
                self.add_sub([x, i, w], origin);
                self.add_sub([i, y, w], origin);
            }
        }
        for s in first..self.subs.len() as u32 {
            self.protect(s);
        }
        Ok(true)
    }

    /// Lawson flips: keep each input triangle's subfacets Delaunay in its
    /// plane. Only edges inside one input triangle (both subfacets with the
    /// same origin) may flip; subsegments never do.
    fn legalize(&mut self, mut stack: Vec<(u32, u32)>) {
        let _prof = crate::prof::scope("conforming::legalize");
        let mut flips = 0usize;
        while let Some((u, v)) = stack.pop() {
            // Exact in-circle tests on fixed 2D coordinates make flipping
            // terminate; the cap is a guard only.
            if flips > MAX_FLIPS {
                return;
            }
            let Some(list) = self.edge_subs.get(&key2(u, v)) else {
                continue;
            };
            let [s1, s2] = match list[..] {
                [p, q] => [p, q],
                _ => continue,
            };
            let (t1, t2) = (self.subs[s1 as usize], self.subs[s2 as usize]);
            if self.facet[t1.origin as usize] != self.facet[t2.origin as usize] {
                continue;
            }
            // Orient s1 as (a, b, w) and expect s2 = (b, a, x).
            let Some(r) = (0..3).find(|&k| key2(t1.v[k], t1.v[(k + 1) % 3]) == key2(u, v)) else {
                continue;
            };
            let (a, b, w) = (t1.v[r], t1.v[(r + 1) % 3], t1.v[(r + 2) % 3]);
            let Some(x) = t2.v.iter().copied().find(|&q| q != a && q != b) else {
                continue;
            };
            let [o, e1, e2] = self.frames[self.facet[t1.origin as usize] as usize];
            let flat = |q: u32| {
                let d = sub(self.pos(q), o);
                [dot(d, e1), dot(d, e2)]
            };
            let (pa, pb, pw, px) = (flat(a), flat(b), flat(w), flat(x));
            // Flip only a strictly convex quad whose fourth point is strictly
            // inside the circle (orientation-independent).
            let side = orient2d(pa, pb, pw);
            if side == Sign::Zero
                || orient2d(pb, pa, px) != side
                || orient2d(pw, px, pa) == orient2d(pw, px, pb)
                || orient2d(pw, px, pa) == Sign::Zero
                || orient2d(pw, px, pb) == Sign::Zero
                || incircle(pa, pb, pw, px) != side
            {
                continue;
            }
            flips += 1;
            self.kill_sub(s1);
            self.kill_sub(s2);
            self.add_sub([a, x, w], t1.origin);
            self.add_sub([x, b, w], t1.origin);
            stack.extend([key2(a, x), key2(x, b), key2(b, w), key2(w, a)]);
        }
    }

    /// Split a subfacet at its circumcentre, or its longest edge if the
    /// circumcentre lies outside it or encroaches a subsegment.
    /// `false` when the floor forbids the split.
    fn split_sub(&mut self, s: u32) -> Result<bool, ConformingError> {
        let _prof = crate::prof::scope("conforming::split_sub");
        let t = self.subs[s as usize].v;
        let [a, b, c] = t.map(|x| self.pos(x));
        let floor = self.params.min_boundary_size.powi(2);
        if self.floored
            && [(a, b), (b, c), (c, a)]
                .iter()
                .all(|&(x, y)| distance_squared(x, y) < floor)
        {
            return Ok(false);
        }
        let centre = tri_circle(a, b, c)
            .map(|(cc, _)| cc)
            .filter(|&cc| inside_triangle(cc, a, b, c));
        let Some(cc) = centre else {
            let (u, v) = [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]
                .into_iter()
                .max_by(|&(p, q), &(r, w)| {
                    distance_squared(self.pos(p), self.pos(q))
                        .total_cmp(&distance_squared(self.pos(r), self.pos(w)))
                })
                .unwrap_or((t[0], t[1]));
            return self.split_edge(u, v);
        };
        if let Some((u, v)) = self.encroached_edge(cc) {
            return self.split_edge(u, v);
        }
        let cavity = self.facet_cavity(s, cc, None);
        match self.split_boundary(cc, &cavity, &[s], None)? {
            // A retry means subfacets changed here, which is progress too.
            Placed::New(_) | Placed::Retry => Ok(true),
            Placed::Blocked => Ok(false),
            Placed::Existing(_) => Err(ConformingError::Recovery { missing_faces: 1 }),
        }
    }

    /// Is `uv` a subsegment: an edge between two facets (or not shared by
    /// exactly two subfacets)?
    fn is_segment(&self, u: u32, v: u32) -> bool {
        match self.edge_subs.get(&key2(u, v)).map(Vec::as_slice) {
            Some(&[s1, s2]) => {
                let f = |s: u32| self.facet[self.subs[s as usize].origin as usize];
                f(s1) != f(s2)
            }
            _ => true,
        }
    }

    /// A subsegment whose diametral sphere strictly contains `p`.
    fn encroached_edge(&self, p: Point3) -> Option<(u32, u32)> {
        let mut best: Option<(u32, u32)> = None;
        for s in self.sub_grid.at(p) {
            let sb = self.subs[s as usize];
            if !sb.alive {
                continue;
            }
            for (u, v) in [(sb.v[0], sb.v[1]), (sb.v[1], sb.v[2]), (sb.v[2], sb.v[0])] {
                let (a, b) = (self.pos(u), self.pos(v));
                let mid = scale(add(a, b), 0.5);
                if distance_squared(p, mid) < 0.25 * distance_squared(a, b)
                    && p != a
                    && p != b
                    && self.is_segment(u, v)
                {
                    let k = key2(u, v);
                    if best.is_none_or(|bk| k < bk) {
                        best = Some(k);
                    }
                }
            }
        }
        best
    }

    /// A subfacet whose diametral sphere strictly contains `p`.
    fn encroached_sub(&self, p: Point3) -> Option<u32> {
        self.sub_grid
            .at(p)
            .filter(|&s| {
                let sb = self.subs[s as usize];
                let [a, b, c] = sb.v.map(|v| self.pos(v));
                sb.alive
                    && ![a, b, c].contains(&p)
                    && tri_circle(a, b, c).is_some_and(|(cc, r)| distance_squared(p, cc) < r * r)
                    && projects_into(p, a, b, c)
            })
            .min()
    }

    /// What to split for a candidate point: a subsegment, a subfacet, or nothing.
    fn encroachment(&self, p: Point3) -> Option<Split> {
        let _prof = crate::prof::scope("conforming::encroachment");
        if let Some((u, v)) = self.encroached_edge(p) {
            return Some(Split::Edge(u, v));
        }
        self.encroached_sub(p).map(Split::Sub)
    }

    /// Whether tet `v` (with centroid `centroid`) is inside the region,
    /// remembered per tet: refinement and sliver removal ask about the same
    /// tets many times, and the winding number is the costliest query.
    fn tet_inside(&mut self, v: [u32; 4], centroid: Point3) -> bool {
        let mut k = v;
        k.sort_unstable();
        if let Some(&known) = self.inside_cache.get(&k) {
            return known;
        }
        let inside = self.inside(centroid);
        self.inside_cache.insert(k, inside);
        inside
    }

    fn inside(&self, p: Point3) -> bool {
        let _prof = crate::prof::scope("conforming::inside");
        self.winding.winding_number(p) >= 0.5
    }

    /// A facet may be triangulated any way: when the tetrahedralization has
    /// the other two triangles of the quad two subfacets of one facet form
    /// across their edge `uv` (missing), swap the subfacets to them. Needs a
    /// strictly convex quad in the facet's frame.
    fn swap_to_mesh(&mut self, u: u32, v: u32) -> bool {
        let _prof = crate::prof::scope("conforming::swap_to_mesh");
        if self.is_segment(u, v) || self.dt.has_edge(u, v) {
            return false;
        }
        let Some(&[s1, s2]) = self.edge_subs.get(&key2(u, v)).map(Vec::as_slice) else {
            return false;
        };
        let (t1, t2) = (self.subs[s1 as usize], self.subs[s2 as usize]);
        let Some(r) = (0..3).find(|&k| key2(t1.v[k], t1.v[(k + 1) % 3]) == key2(u, v)) else {
            return false;
        };
        let (a, b, w) = (t1.v[r], t1.v[(r + 1) % 3], t1.v[(r + 2) % 3]);
        let Some(x) = t2.v.iter().copied().find(|&q| q != a && q != b) else {
            return false;
        };
        if !self.dt.has_face([a, x, w]) || !self.dt.has_face([x, b, w]) {
            return false;
        }
        let [o, e1, e2] = self.frames[self.facet[t1.origin as usize] as usize];
        let flat = |q: u32| {
            let d = sub(self.pos(q), o);
            [dot(d, e1), dot(d, e2)]
        };
        let (pa, pb, pw, px) = (flat(a), flat(b), flat(w), flat(x));
        let side = orient2d(pa, pb, pw);
        if side == Sign::Zero
            || orient2d(pb, pa, px) != side
            || orient2d(pw, px, pa) == orient2d(pw, px, pb)
            || orient2d(pw, px, pa) == Sign::Zero
            || orient2d(pw, px, pb) == Sign::Zero
        {
            return false;
        }
        self.kill_sub(s1);
        self.kill_sub(s2);
        self.add_sub([a, x, w], t1.origin);
        self.add_sub([x, b, w], t1.origin);
        let n = self.subs.len() as u32;
        self.protect(n - 2);
        self.protect(n - 1);
        true
    }

    /// Adopt the tetrahedralization's own triangulation of facet `f` when it
    /// has a complete one: the faces whose vertices all lie on the facet (by
    /// label), inside its subfacets, covering exactly their area, with each
    /// subsegment of the facet an edge of exactly one of them. They replace
    /// the facet's subfacets (each keeping the input triangle holding its
    /// centroid).
    fn adopt_facet(&mut self, f: u32, subs: &[u32]) -> bool {
        let _prof = crate::prof::scope("conforming::adopt_facet");
        let [o, e1, e2] = self.frames[f as usize];
        let flat = |b: &Self, q: u32| {
            let d = sub(b.pos(q), o);
            [dot(d, e1), dot(d, e2)]
        };
        let area = |b: &Self, t: [u32; 3]| {
            let [x, y, z] = t.map(|q| b.pos(q));
            norm(cross(sub(y, x), sub(z, x))) / 2.0
        };
        let mut verts: Vec<u32> = subs.iter().flat_map(|&s| self.subs[s as usize].v).collect();
        verts.sort_unstable();
        verts.dedup();
        // Candidate faces: in the stars of the facet's vertices, all three
        // vertices on the facet's plane.
        let mut cand: Vec<[u32; 3]> = Vec::new();
        for &v in &verts {
            for t in self.dt.star(v) {
                let Some(tv) = self.dt.tet(t) else {
                    continue;
                };
                for g in FACES {
                    let face = key3(g.map(|k| tv[k]));
                    if face.iter().all(|x| verts.binary_search(x).is_ok()) {
                        cand.push(face);
                    }
                }
            }
        }
        cand.sort_unstable();
        cand.dedup();
        // Keep those whose centroid lies in a subfacet; orient them like it.
        let mut new: Vec<([u32; 3], u32)> = Vec::new();
        for c in cand {
            let pc = c.map(|q| flat(self, q));
            let centre = [
                (pc[0][0] + pc[1][0] + pc[2][0]) / 3.0,
                (pc[0][1] + pc[1][1] + pc[2][1]) / 3.0,
            ];
            let home = subs.iter().copied().find(|&s| {
                let [a, b, d] = self.subs[s as usize].v.map(|q| flat(self, q));
                let side = orient2d(a, b, d);
                side != Sign::Zero
                    && orient2d(a, b, centre) == side
                    && orient2d(b, d, centre) == side
                    && orient2d(d, a, centre) == side
            });
            let Some(home) = home else {
                continue;
            };
            let [a, b, d] = self.subs[home as usize].v.map(|q| flat(self, q));
            let side = orient2d(a, b, d);
            let oriented = if orient2d(pc[0], pc[1], pc[2]) == side {
                c
            } else {
                [c[0], c[2], c[1]]
            };
            new.push((oriented, self.subs[home as usize].origin));
        }
        let want: f64 = subs
            .iter()
            .map(|&s| area(self, self.subs[s as usize].v))
            .sum();
        let have: f64 = new.iter().map(|&(t, _)| area(self, t)).sum();
        if new.is_empty() || (want - have).abs() > 1e-9 * want {
            return false;
        }
        // Each edge of the new triangulation: on two of them, or a
        // subsegment of the facet on one.
        let mut edges: HashMap<(u32, u32), usize> = HashMap::new();
        for &(t, _) in &new {
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                *edges.entry(key2(a, b)).or_default() += 1;
            }
        }
        let mut boundary: Vec<(u32, u32)> = Vec::new();
        for &s in subs {
            let v = self.subs[s as usize].v;
            for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
                if self.is_segment(a, b) {
                    boundary.push(key2(a, b));
                }
            }
        }
        for (e, n) in &edges {
            let on = boundary.contains(e);
            if (on && *n != 1) || (!on && *n != 2) {
                return false;
            }
        }
        if boundary.iter().any(|e| !edges.contains_key(e)) {
            return false;
        }
        for &s in subs {
            self.kill_sub(s);
        }
        for (t, origin) in new {
            self.add_sub(t, origin);
        }
        true
    }

    /// Recover a missing subfacet whose edges are all in by flips: remove the
    /// edges crossing it, one edge-removal flip at a time.
    fn flip_recover(&mut self, s: u32) -> bool {
        let _prof = crate::prof::scope("conforming::flip_recover");
        let v = self.subs[s as usize].v;
        for _ in 0..32 {
            if self.dt.has_face(v) {
                return true;
            }
            let crossing = self.dt.crossing_edges(v);
            if !crossing
                .into_iter()
                .any(|(x, y)| self.dt.remove_edge_for(x, y, Some(v), false))
            {
                return false;
            }
        }
        self.dt.has_face(v)
    }

    /// Recover every subsegment and subfacet: what is in the
    /// tetrahedralization is constrained at each pass; missing subsegments are
    /// split; missing subfacets are recovered by facet diagonal swaps, then by
    /// flipping away crossing edges, then by splitting. Fails when the number
    /// missing stops going down.
    fn recover(&mut self) -> Result<(), ConformingError> {
        let _prof = crate::prof::scope("conforming::recover");
        // Passes since the number missing last went down.
        let (mut best, mut stalled) = (usize::MAX, 0usize);
        loop {
            self.step()?;
            let tets = self.dt.tets();
            let mut edges: HashSet<(u32, u32)> = HashSet::with_capacity(tets.len() * 2);
            let mut faces: HashSet<[u32; 3]> = HashSet::with_capacity(tets.len() * 3);
            for t in &tets {
                for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
                    edges.insert(key2(t[i], t[j]));
                }
                for f in FACES {
                    faces.insert(key3(f.map(|k| t[k])));
                }
            }
            let live: Vec<u32> = (0..self.subs.len() as u32)
                .filter(|&s| self.subs[s as usize].alive)
                .collect();
            // A recovered subsegment stays recovered.
            for &s in &live {
                let v = self.subs[s as usize].v;
                if faces.contains(&key3(v)) {
                    self.dt.constrain_face(v);
                }
                for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
                    if edges.contains(&key2(a, b)) && self.is_segment(a, b) {
                        self.dt.constrain_edge(a, b);
                    }
                }
            }
            let mut missing_edges: Vec<(u32, u32)> = live
                .iter()
                .flat_map(|&s| {
                    let v = self.subs[s as usize].v;
                    [key2(v[0], v[1]), key2(v[1], v[2]), key2(v[2], v[0])]
                })
                .filter(|e| !edges.contains(e))
                .filter(|&(u, v)| self.is_segment(u, v))
                .collect();
            missing_edges.sort_unstable();
            missing_edges.dedup();
            if !missing_edges.is_empty() {
                let mut progress = false;
                let count = missing_edges.len();
                for (u, v) in missing_edges {
                    if self
                        .edge_subs
                        .get(&key2(u, v))
                        .is_some_and(|l| !l.is_empty())
                    {
                        progress |= self.split_edge(u, v)?;
                    }
                }
                if !progress {
                    return Err(ConformingError::Recovery {
                        missing_faces: count,
                    });
                }
                continue;
            }
            let missing: Vec<u32> = live
                .iter()
                .copied()
                .filter(|&s| !faces.contains(&key3(self.subs[s as usize].v)))
                .collect();
            if missing.is_empty() {
                // Everything is in: from now on no insertion may delete a
                // subfacet.
                for s in live {
                    self.dt.constrain_face(self.subs[s as usize].v);
                }
                self.face_constraints = true;
                self.missing.clear();
                return Ok(());
            }
            let mut progress = false;
            let count = missing.len();
            if count < best {
                (best, stalled) = (count, 0);
            } else {
                stalled += 1;
                if stalled > 8 {
                    return Err(ConformingError::Recovery {
                        missing_faces: count,
                    });
                }
            }
            // Flips first (no new points): swap facet diagonals to the ones
            // the tetrahedralization has, then flip away crossing edges;
            // split what they cannot recover.
            let mut left = Vec::new();
            // Facets with a missing subfacet that the mesh triangulates
            // completely in its own way take that triangulation.
            let mut facets: Vec<u32> = missing
                .iter()
                .map(|&s| self.facet[self.subs[s as usize].origin as usize])
                .collect();
            facets.sort_unstable();
            facets.dedup();
            let mut by_facet: HashMap<u32, Vec<u32>> = HashMap::new();
            for &s in &live {
                let f = self.facet[self.subs[s as usize].origin as usize];
                if facets.binary_search(&f).is_ok() {
                    by_facet.entry(f).or_default().push(s);
                }
            }
            for f in facets {
                let subs = by_facet.remove(&f).unwrap_or_default();
                if subs.len() > 1 && self.adopt_facet(f, &subs) {
                    progress = true;
                }
            }
            if progress {
                continue;
            }
            for &s in &missing {
                if !self.subs[s as usize].alive {
                    continue;
                }
                let v = self.subs[s as usize].v;
                for (x, y) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
                    if self.swap_to_mesh(x, y) {
                        progress = true;
                        break;
                    }
                }
            }
            if progress {
                continue;
            }
            for s in missing {
                if self.flip_recover(s) {
                    progress = true;
                } else {
                    left.push(s);
                }
            }
            if progress {
                continue;
            }
            for s in left {
                if self.subs[s as usize].alive && !self.dt.has_face(self.subs[s as usize].v) {
                    progress |= self.split_sub(s)?;
                }
            }
            if !progress {
                return Err(ConformingError::Recovery {
                    missing_faces: count,
                });
            }
        }
    }

    /// The tet's circumcentre, if it is trustworthy: exactly inside the tet's
    /// circumsphere. For a near-flat tet the computed centre can be far off.
    fn reliable_centre(p: [Point3; 4]) -> Option<(Point3, f64)> {
        tet_sphere(p).filter(|&(cc, _)| insphere(p[0], p[1], p[2], p[3], cc) == Sign::Positive)
    }

    /// Kill a tet whose circumcentre is unreliable by splitting its longest
    /// boundary edge (an edge's midpoint is always inside the circumsphere).
    /// Returns whether anything was split.
    fn split_flat(&mut self, v: [u32; 4]) -> Result<bool, ConformingError> {
        let _prof = crate::prof::scope("conforming::split_flat");
        let best = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            .into_iter()
            .map(|(i, j)| key2(v[i], v[j]))
            .filter(|e| self.edge_subs.get(e).is_some_and(|l| !l.is_empty()))
            .max_by(|&(a, b), &(c, d)| {
                distance_squared(self.pos(a), self.pos(b))
                    .total_cmp(&distance_squared(self.pos(c), self.pos(d)))
                    .then((a, b).cmp(&(c, d)))
            });
        match best {
            Some((u, w)) => self.split_edge(u, w),
            None => Ok(false),
        }
    }

    /// The target size at `p`: the size field's, capped at the maximum.
    fn size_at(&self, p: Point3) -> f64 {
        let h = self.params.max_circumradius;
        self.params.size.as_ref().map_or(h, |f| f.at(p).min(h))
    }

    /// Count one kernel step; an error once the work budget is spent.
    fn step(&mut self) -> Result<(), ConformingError> {
        self.work += 1;
        if self.work > self.params.max_work {
            return Err(ConformingError::Work {
                steps: self.work,
                limit: self.params.max_work,
            });
        }
        Ok(())
    }

    /// Count an attempt on a tet; `false` once it has had its share.
    fn attempt(&mut self, v: [u32; 4]) -> bool {
        let mut k = v;
        k.sort_unstable();
        let n = self.attempts.entry(k).or_insert(0);
        *n += 1;
        *n <= MAX_ATTEMPTS
    }

    /// Refine tet `slot` (with centroid `from`) at the point `cc`: split the
    /// subsegment or subfacet `cc` encroaches, or the subfacet between the
    /// tet and `cc`; otherwise insert `cc`. `Some(true)` for a boundary
    /// split (the tet may survive it), `Some(false)` for an inserted point,
    /// `None` when nothing could be done.
    fn refine_at(
        &mut self,
        slot: u32,
        from: Point3,
        cc: Point3,
        fallback: bool,
    ) -> Result<Option<bool>, ConformingError> {
        let _prof = crate::prof::scope("conforming::refine_at");
        let encroached = self.encroachment(cc);
        let any = encroached.is_some();
        let split = match encroached {
            Some(Split::Edge(u, w)) => self.split_edge(u, w)?,
            Some(Split::Sub(s)) => self.split_sub(s)?,
            None => false,
        };
        if split {
            return Ok(Some(true));
        }
        if any && !fallback {
            return Ok(None);
        }
        // No boundary split (none needed, or the floor or the constraints
        // refused it): insert the circumcentre if the tet sees it.
        let w = self.dt.walk(slot, from, cc);
        let split = match w {
            Walk::Clear => {
                self.dt
                    .set_planes(self.dt.points().len() as u32, Vec::new());
                let placed = self.insert_point(cc, &[], &[], false)?;
                return Ok(matches!(placed, Placed::New(_)).then_some(false));
            }
            Walk::Constraint(f) => match self.face_sub.get(&f).copied() {
                Some(s) => self.split_sub(s)?,
                None => false,
            },
            Walk::Lost => false,
        };
        Ok(split.then_some(true))
    }

    /// Recover the subfacets a forced split deleted, locally: a missing
    /// subsegment is split; then a facet diagonal is swapped to the mesh's,
    /// or crossing edges are flipped away, or the subfacet is split. What is
    /// still missing after a few rounds goes to a full recovery pass.
    fn repair(&mut self) -> Result<(), ConformingError> {
        let _prof = crate::prof::scope("conforming::repair");
        let floored = self.floored;
        self.floored = false;
        let mut queue = std::mem::take(&mut self.missing);
        queue.sort_unstable();
        queue.dedup();
        // A cascade: local splits that delete more subfacets than they
        // recover, each round slower than the last. Local repair cannot win
        // it, but a full recovery pass settles it quickly (Thingi10k), so hand
        // over once the queue has grown three rounds running past `CASCADE`.
        let mut last = queue.len();
        let mut growing = 0usize;
        for _ in 0..16 {
            if queue.is_empty() {
                break;
            }
            queue.sort_unstable();
            queue.dedup();
            growing = if queue.len() > last { growing + 1 } else { 0 };
            last = queue.len();
            if growing >= 3 && queue.len() > CASCADE {
                break;
            }
            let mut next = Vec::new();
            for s in queue {
                self.step()?;
                if !self.subs[s as usize].alive {
                    continue;
                }
                let v = self.subs[s as usize].v;
                if self.dt.has_face(v) {
                    self.protect(s);
                    continue;
                }
                let edges = [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])];
                let fixed = if let Some(&(a, b)) = edges
                    .iter()
                    .find(|&&(a, b)| self.is_segment(a, b) && !self.dt.has_edge(a, b))
                {
                    self.split_edge(a, b)?
                } else if edges.iter().any(|&(a, b)| self.swap_to_mesh(a, b)) {
                    true
                } else if self.flip_recover(s) {
                    self.protect(s);
                    true
                } else {
                    self.split_sub(s)?
                };
                if !fixed || self.subs[s as usize].alive && !self.dt.has_face(v) {
                    next.push(s);
                }
            }
            // New subfacets the splits and swaps left missing.
            next.append(&mut self.missing);
            // Swapped-in subfacets that are present need their constraint.
            queue = next;
        }
        let result = if queue
            .iter()
            .any(|&s| self.subs[s as usize].alive && !self.dt.has_face(self.subs[s as usize].v))
        {
            self.recover()
        } else {
            Ok(())
        };
        self.floored = floored;
        result
    }

    /// Delaunay refinement on radius-edge ratio and size, inside the region.
    fn refine(&mut self) -> Result<(), ConformingError> {
        let _prof = crate::prof::scope("conforming::refine");
        loop {
            self.step()?;
            if !self.missing.is_empty() {
                self.repair()?;
            }
            let Some((slot, v)) = self.tet_queue.pop_front() else {
                return Ok(());
            };
            if self.dt.tet(slot) != Some(v) {
                continue;
            }
            let p = v.map(|x| self.pos(x));
            let centroid = scale(add(add(p[0], p[1]), add(p[2], p[3])), 0.25);
            if !self.tet_inside(v, centroid) {
                continue;
            }
            let Some((cc, r)) = tet_sphere(p) else {
                continue;
            };
            let reliable = Self::reliable_centre(p).is_some();
            let shortest = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
                .iter()
                .map(|&(i, j)| distance_squared(p[i], p[j]))
                .fold(f64::INFINITY, f64::min)
                .sqrt();
            // A tet with an edge below the floor owes its shape to a small
            // input feature; refining it would chase the feature forever.
            let bad = r > self.size_at(centroid)
                || (r > self.params.radius_edge_bound * shortest
                    && shortest >= self.params.min_boundary_size);
            if !bad || !self.attempt(v) {
                continue;
            }
            let split = if reliable {
                self.refine_at(slot, centroid, cc, true)? == Some(true)
            } else {
                self.split_flat(v)?
            };
            // The tet may survive a boundary split; look at it again. If the
            // floor refused the split, the tet stays as it is.
            if split {
                self.tet_queue.push_back((slot, v));
            }
        }
    }

    /// Remove slivers by flips (no new points): for each sliver inside the
    /// region, an edge removal or 2–3 flip that makes the worst tet fatter.
    /// Constrained faces and edges never flip, so the boundary stays.
    fn flip_slivers(&mut self) -> usize {
        let _prof = crate::prof::scope("conforming::flip_slivers");
        let mut flips = 0;
        for _ in 0..4 {
            let mut changed = 0;
            for slot in self.dt.real_slots() {
                let Some(v) = self.dt.tet(slot) else {
                    continue;
                };
                let p = v.map(|x| self.pos(x));
                let min = dihedral_angles(p)
                    .iter()
                    .copied()
                    .fold(f64::INFINITY, f64::min);
                if min >= self.params.sliver_degrees {
                    continue;
                }
                let centroid = scale(add(add(p[0], p[1]), add(p[2], p[3])), 0.25);
                if !self.tet_inside(v, centroid) {
                    continue;
                }
                let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
                let done = edges
                    .iter()
                    .any(|&(i, j)| self.dt.remove_edge_for(v[i], v[j], None, true))
                    || FACES
                        .iter()
                        .any(|f| self.dt.flip_face(f.map(|k| v[k]), true));
                if done {
                    changed += 1;
                }
            }
            flips += changed;
            if changed == 0 {
                break;
            }
        }
        flips
    }

    /// Refine at the circumcentres of slivers, a few rounds.
    fn remove_slivers(&mut self) -> Result<(), ConformingError> {
        let _prof = crate::prof::scope("conforming::remove_slivers");
        for _ in 0..self.params.sliver_rounds {
            let mut inserted = 0;
            for slot in self.dt.real_slots() {
                let Some(v) = self.dt.tet(slot) else {
                    continue;
                };
                let p = v.map(|x| self.pos(x));
                let centroid = scale(add(add(p[0], p[1]), add(p[2], p[3])), 0.25);
                let min = dihedral_angles(p)
                    .iter()
                    .copied()
                    .fold(f64::INFINITY, f64::min);
                if min >= self.params.sliver_degrees || !self.tet_inside(v, centroid) {
                    continue;
                }
                if !self.attempt(v) {
                    continue;
                }
                let done = match Self::reliable_centre(p) {
                    Some((cc, _)) => self.refine_at(slot, centroid, cc, false)?.is_some(),
                    None => self.split_flat(v)?,
                };
                if done {
                    inserted += 1;
                }
                if !self.missing.is_empty() {
                    self.repair()?;
                }
            }
            self.tet_queue.clear();
            if inserted == 0 {
                break;
            }
        }
        Ok(())
    }
}

/// Mesh the region bounded by a closed, outward, non-self-intersecting
/// surface conformally.
pub fn conforming_mesh(
    points: &[Point3],
    triangles: &[[u32; 3]],
    params: &ConformingParams,
) -> Result<ConformingMesh, ConformingError> {
    let _prof = crate::prof::scope("conforming::conforming_mesh");
    let soup = TriangleSoup::new(points, triangles)
        .map_err(|e| ConformingError::Surface(e.to_string()))?;
    let mut sorted = points.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if sorted.windows(2).any(|w| w[0] == w[1]) {
        return Err(ConformingError::Surface("coincident vertices".into()));
    }
    if !(params.max_circumradius > 0.0 && params.max_circumradius.is_finite()) {
        return Err(ConformingError::Surface(
            "size bound must be positive".into(),
        ));
    }
    let budget = if params.max_points > 0 {
        params.max_points
    } else {
        // Elements fill the volume at the maximum size, plus the surface at
        // its local size; where the size is smaller, its graded layer (about
        // h / gradation thick) adds ∫ 1/(h + g d)³ dd = 1/(2 g h²) per unit
        // area over what the maximum size would put there.
        let h = params.max_circumradius;
        let field = params.size.as_deref();
        let g = field.map_or(0.0, SizeField::gradation);
        let (mut surface, mut volume) = (0.0f64, 0.0f64);
        for t in triangles {
            let [a, b, c] = t.map(|v| points[v as usize]);
            let area = norm(cross(sub(b, a), sub(c, a))) / 2.0;
            let local = field.map_or(h, |f| {
                f.at(scale(add(add(a, b), c), 1.0 / 3.0)).clamp(f.min(), h)
            });
            surface += area / (local * local);
            if g > 0.0 && local < h {
                surface += area * (1.0 / (local * local) - 1.0 / (h * h)) / (2.0 * g);
            }
            // Relative to a vertex: far from the origin the products cancel.
            let o = points[0];
            volume += dot(sub(a, o), cross(sub(b, o), sub(c, o))) / 6.0;
        }
        let estimate = volume.abs() / h.powi(3) + surface + points.len() as f64;
        ((20.0 * estimate) as usize).clamp(10_000, 4_000_000)
    };
    let work = if params.max_work > 0 {
        params.max_work
    } else {
        64 * budget as u64
    };
    let local_params = ConformingParams {
        max_points: budget,
        max_work: work,
        ..params.clone()
    };
    let params = &local_params;
    // Each input vertex lies on the facets of its triangles: points sharing a
    // facet are coplanar by construction, whatever rounding says.
    let facet = facets(points, triangles, &params.groups);
    let mut planes: Vec<Vec<u32>> = vec![Vec::new(); points.len()];
    for (t, tri) in triangles.iter().enumerate() {
        for &v in tri {
            planes[v as usize].push(facet[t]);
        }
    }
    for l in &mut planes {
        l.sort_unstable();
        l.dedup();
    }
    let dt = Delaunay::with_planes(points.to_vec(), params.seed, planes)
        .map_err(ConformingError::Delaunay)?;
    if !dt.duplicates().is_empty() {
        return Err(ConformingError::Surface("coincident vertices".into()));
    }
    let mut b = Builder {
        dt,
        subs: Vec::new(),
        edge_subs: HashMap::new(),
        sub_grid: Grid {
            cell: params.max_circumradius,
            cells: HashMap::new(),
            big: Vec::new(),
        },
        face_sub: HashMap::new(),
        winding: FastWindingNumber::new(soup),
        params,
        inputs: points.len(),
        boundary_steiner: 0,
        tet_queue: VecDeque::new(),
        attempts: HashMap::new(),
        frames: triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|v| points[v as usize]);
                plane_frame(a, b, c)
            })
            .collect(),
        floored: false,
        missing: Vec::new(),
        face_constraints: false,
        facet,
        work: 0,
        inside_cache: FxHashMap::default(),
    };
    for (t, tri) in triangles.iter().enumerate() {
        b.add_sub(*tri, t as u32);
    }
    // Plane labels keep the tetrahedralization free of tets flat within a
    // facet; this only mops up any that slipped through.
    b.dt.remove_flat_tets();
    // Start each facet from its Delaunay triangulation.
    let mut interior: Vec<(u32, u32)> = b
        .edge_subs
        .keys()
        .copied()
        .filter(|&(u, v)| !b.is_segment(u, v))
        .collect();
    interior.sort_unstable();
    b.legalize(interior);

    b.recover()?;
    b.unflatten_all();
    // Recovery's flips and shrunk cavities leave the mesh short of
    // constrained Delaunay; restore it, so refinement's cavities are regular.
    let n = b.dt.points().len();
    b.dt.delaunize(64 * n);
    b.tet_queue =
        b.dt.real_slots()
            .into_iter()
            .filter_map(|s| b.dt.tet(s).map(|v| (s, v)))
            .collect();
    b.floored = true;
    b.refine()?;
    b.flip_slivers();
    b.remove_slivers()?;
    b.flip_slivers();
    // Constraints keep every recovered subfacet; only a split that landed on
    // an existing vertex leaves halves to recover.
    b.floored = false;
    b.recover()?;

    // Region marking: flood from the hull, toggling inside/outside at every
    // subfacet crossed. Exact and topological, so tets thinner than rounding
    // near the surface are classified as reliably as any other.
    let all = b.dt.points().to_vec();
    let tets = b.dt.tets();
    let subs: HashSet<[u32; 3]> = b
        .subs
        .iter()
        .filter(|s| s.alive)
        .map(|s| key3(s.v))
        .collect();
    let mut by_face: HashMap<[u32; 3], Vec<u32>> = HashMap::with_capacity(tets.len() * 2);
    for (i, t) in tets.iter().enumerate() {
        for f in FACES {
            by_face
                .entry(key3(f.map(|k| t[k])))
                .or_default()
                .push(i as u32);
        }
    }
    let mut label: Vec<Option<bool>> = vec![None; tets.len()];
    let mut queue: VecDeque<u32> = VecDeque::new();
    for (i, t) in tets.iter().enumerate() {
        for f in FACES {
            let k = key3(f.map(|j| t[j]));
            if by_face.get(&k).is_some_and(|l| l.len() == 1) && label[i].is_none() {
                label[i] = Some(subs.contains(&k));
                queue.push_back(i as u32);
            }
        }
    }
    while let Some(i) = queue.pop_front() {
        let here = label[i as usize].unwrap_or(false);
        let t = tets[i as usize];
        for f in FACES {
            let k = key3(f.map(|j| t[j]));
            let crossing = subs.contains(&k);
            for &n in by_face.get(&k).map_or(&[][..], Vec::as_slice) {
                if n == i {
                    continue;
                }
                let want = here ^ crossing;
                match label[n as usize] {
                    None => {
                        label[n as usize] = Some(want);
                        queue.push_back(n);
                    }
                    Some(have) if have != want => {
                        return Err(ConformingError::Recovery { missing_faces: 1 });
                    }
                    Some(_) => {}
                }
            }
        }
    }
    let mut kept: Vec<[u32; 4]> = tets
        .iter()
        .zip(&label)
        .filter(|(_, l)| **l == Some(true))
        .map(|(t, _)| *t)
        .collect();

    // A tet flat to rounding with two subfacets as faces fills the gap
    // between two coplanar surface triangles and the other diagonal of their
    // quad. Drop it and let its other two faces be the boundary: the surface
    // moves by rounding only, and each new face keeps the input triangle that
    // holds its centroid.
    let mut origin: HashMap<[u32; 3], u32> = b
        .subs
        .iter()
        .filter(|s| s.alive)
        .map(|s| (key3(s.v), s.origin))
        .collect();
    let mut flattened = 0usize;
    loop {
        let mut drop = None;
        for (i, t) in kept.iter().enumerate() {
            let faces: Vec<[u32; 3]> = FACES.iter().map(|f| f.map(|k| t[k])).collect();
            let on: Vec<usize> = (0..4)
                .filter(|&k| origin.contains_key(&key3(faces[k])))
                .collect();
            if on.len() != 2 {
                continue;
            }
            let group = |f: [u32; 3]| {
                params
                    .groups
                    .get(origin[&key3(f)] as usize)
                    .copied()
                    .unwrap_or(0)
            };
            if group(faces[on[0]]) != group(faces[on[1]]) {
                continue;
            }
            let p = t.map(|v| all[v as usize]);
            let longest = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
                .iter()
                .map(|&(i, j)| distance_squared(p[i], p[j]))
                .fold(0.0, f64::max)
                .sqrt();
            let volume = dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))).abs() / 6.0;
            // Rounding moves points by a few ulps of their coordinates, which
            // far from the origin is much more than of the tet's size.
            let magnitude = p.iter().flatten().fold(0.0f64, |m, c| m.max(c.abs()));
            let flat = (1e-12 * longest).max(64.0 * f64::EPSILON * magnitude);
            // Only a convex quad's two diagonals cover the same area; a flat
            // tet with a vertex inside the other three's triangle does not.
            let area = |f: [u32; 3]| {
                let [a, bb, c] = f.map(|v| all[v as usize]);
                norm(cross(sub(bb, a), sub(c, a))) / 2.0
            };
            let subs_area: f64 = on.iter().map(|&k| area(faces[k])).sum();
            let other_area: f64 = (0..4)
                .filter(|k| !on.contains(k))
                .map(|k| area(faces[k]))
                .sum();
            if volume <= flat * longest * longest
                && (subs_area - other_area).abs() <= 1e-9 * subs_area.max(other_area)
            {
                drop = Some((i, faces, on));
                break;
            }
        }
        let Some((i, faces, on)) = drop else {
            break;
        };
        let old: Vec<([u32; 3], u32)> = on
            .iter()
            .map(|&k| (faces[k], origin[&key3(faces[k])]))
            .collect();
        for &k in &on {
            origin.remove(&key3(faces[k]));
        }
        for k in (0..4).filter(|k| !on.contains(k)) {
            let f = faces[k];
            let c = scale(
                add(
                    add(all[f[0] as usize], all[f[1] as usize]),
                    all[f[2] as usize],
                ),
                1.0 / 3.0,
            );
            let o = old
                .iter()
                .find(|(g, _)| {
                    let [a, bb, cc] = g.map(|v| all[v as usize]);
                    inside_triangle(c, a, bb, cc)
                })
                .map_or(old[0].1, |&(_, o)| o);
            origin.insert(key3(f), o);
        }
        kept.swap_remove(i);
        flattened += 1;
    }
    kept.sort_unstable();
    let mut uses: HashMap<[u32; 3], ([u32; 3], u32)> = HashMap::new();
    let mut order: Vec<[u32; 3]> = Vec::new();
    for t in &kept {
        for f in FACES {
            let face = f.map(|k| t[k]);
            let k = key3(face);
            match uses.get_mut(&k) {
                Some(e) => e.1 += 1,
                None => {
                    uses.insert(k, (face, 1));
                    order.push(k);
                }
            }
        }
    }
    let mut boundary = Vec::new();
    let mut stray = 0usize;
    for k in &order {
        let (face, count) = uses[k];
        if count == 1 {
            match origin.get(k) {
                Some(&o) => boundary.push((face, o)),
                None => stray += 1,
            }
        }
    }
    let missing = origin.len().saturating_sub(boundary.len()) + stray;
    if missing > 0 {
        return Err(ConformingError::Recovery {
            missing_faces: missing,
        });
    }

    // Compact to the vertices in use.
    let mut map = vec![u32::MAX; all.len()];
    let mut out_points = Vec::new();
    let mut remap = |v: u32| {
        if map[v as usize] == u32::MAX {
            map[v as usize] = out_points.len() as u32;
            out_points.push(all[v as usize]);
        }
        map[v as usize]
    };
    let tets: Vec<[u32; 4]> = kept.iter().map(|t| t.map(&mut remap)).collect();
    let boundary: Vec<([u32; 3], u32)> = boundary
        .into_iter()
        .map(|(f, o)| (f.map(&mut remap), o))
        .collect();
    // Smooth the free vertices (those on no boundary face) of poor tets.
    let mut fixed = vec![false; out_points.len()];
    for &(f, _) in &boundary {
        for v in f {
            fixed[v as usize] = true;
        }
    }
    let smoothed = smooth_slivers(
        &mut out_points,
        &tets,
        &fixed,
        &SmoothParams {
            passes: params.smooth_passes,
            ..SmoothParams::default()
        },
    )
    .moves;
    let slivers_left = tets
        .iter()
        .filter(|t| {
            dihedral_angles(t.map(|v| out_points[v as usize]))
                .iter()
                .any(|&a| a < params.sliver_degrees)
        })
        .count();
    let total = b.dt.points().len();
    Ok(ConformingMesh {
        points: out_points,
        tets,
        boundary,
        boundary_steiner: b.boundary_steiner,
        interior_steiner: total - b.inputs - b.boundary_steiner,
        slivers_left,
        flattened,
        smoothed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let t = quads
            .iter()
            .flat_map(|&[a, b, c, d]| [[a, b, c], [a, c, d]])
            .collect();
        (p, t)
    }

    fn volume(m: &ConformingMesh) -> f64 {
        m.tets
            .iter()
            .map(|t| {
                let [a, b, c, d] = t.map(|v| m.points[v as usize]);
                dot(sub(b, a), cross(sub(c, a), sub(d, a))) / 6.0
            })
            .sum()
    }

    #[test]
    fn cube_is_meshed_conformally() {
        let (p, t) = cube();
        let m = conforming_mesh(&p, &t, &ConformingParams::with_size(0.3)).unwrap();
        assert!((volume(&m) - 1.0).abs() < 1e-12, "volume {}", volume(&m));
        // Every boundary triangle lies on its input triangle's plane.
        for &(f, o) in &m.boundary {
            let [a, b, c] = t[o as usize].map(|v| p[v as usize]);
            for q in f.map(|v| m.points[v as usize]) {
                let d = dot(sub(q, a), cross(sub(b, a), sub(c, a)));
                assert!(d.abs() < 1e-12);
            }
        }
        let area: f64 = m
            .boundary
            .iter()
            .map(|(f, _)| {
                let [a, b, c] = f.map(|v| m.points[v as usize]);
                norm(cross(sub(b, a), sub(c, a))) / 2.0
            })
            .sum();
        assert!((area - 6.0).abs() < 1e-12);
    }

    #[test]
    fn spheres_are_equidistant() {
        let tets = [
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            [
                [0.3, -0.2, 0.9],
                [1.7, 0.4, -0.3],
                [-0.5, 1.1, 0.2],
                [0.9, 0.8, 1.6],
            ],
        ];
        for p in tets {
            let (cc, r) = tet_sphere(p).unwrap();
            for q in p {
                assert!((norm(sub(q, cc)) - r).abs() < 1e-12);
            }
            let (c2, r2) = tri_circle(p[0], p[1], p[2]).unwrap();
            for q in &p[..3] {
                assert!((norm(sub(*q, c2)) - r2).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn coincident_vertices_are_refused() {
        let (mut p, t) = cube();
        p.push(p[0]);
        let mut t2 = t;
        t2.push([8, 1, 2]);
        assert!(matches!(
            conforming_mesh(&p, &t2, &ConformingParams::with_size(0.3)),
            Err(ConformingError::Surface(_))
        ));
    }
}
