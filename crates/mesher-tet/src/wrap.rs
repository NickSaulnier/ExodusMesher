//! Wrapping a triangle soup (repair rung 3), after Portaneri, Rouxel-Labbé,
//! Hemmer, Cohen-Steiner and Alliez, "Alpha Wrapping with an Offset"
//! (SIGGRAPH 2022), implemented from the paper's description only.
//!
//! A Delaunay tetrahedralization of a box around the input starts with
//! every real cell *inside* and the ghosts *outside*. A *gate* is a facet
//! between an outside and an inside cell. It is α-traversable when the
//! smallest Delaunay ball through its corners has radius at least α: those
//! balls have their centres on the facet's dual Voronoi edge, from the
//! outside cell's circumcentre (far along the hull normal for a ghost) to the
//! inside cell's. Traversable gates are taken largest ball first. Traversing
//! one into cell `c`:
//!
//! 1. if the Voronoi edge, followed from its outside end, meets the offset
//!    surface (the points at distance δ from the input), the first such
//!    point is inserted;
//! 2. else if the input enters `c`'s circumball, the projection of its
//!    circumcentre onto the offset surface is inserted (the paper's rule
//!    takes the case of a circumcentre within δ of the input; without the
//!    rest, a large cell whose centre is clear of the input could be carved
//!    with the input inside it);
//! 3. else `c` is carved: labelled outside, its other facets become gates.
//!
//! Either inserted point is strictly inside the circumball of a cell of the
//! gate (a ball through the facet's corners centred on the edge lies in the
//! union of the two end balls; rule 2's point is on the segment from the
//! centre to the nearest input point, nearer the centre than the input is
//! or within δ < α of it), so the insertion destroys the gate; the cells it
//! creates are inside. Rule 2's point is more than δ from every vertex (the
//! ball is empty), which bounds the insertions. No carved cell holds input. When no traversable gate is left, the boundary of the outside
//! region is the wrap: its vertices are on the offset surface, it encloses
//! the input, and as a subcomplex of a Delaunay tetrahedralization no two of
//! its triangles cross.
//!
//! *Winding guard* (`WrapParams::winding_guard`, for open inputs). The plain
//! wrap of an open surface whose holes are wider than about 2α goes through
//! them and wraps both sides of its walls. With the guard, points where the
//! input's generalized winding number is at least 1/2 in magnitude count as
//! input: where a Voronoi edge enters that solid (bisection on the winding
//! number between two sphere-tracing samples, both in balls clear of the
//! input) the entry point is inserted (on the cap across the hole, so not at
//! δ from the input), and a cell whose circumcentre is in it stays inside.
//! The entry point is the centre of an empty ball through the facet's corners
//! of radius at least α, so these insertions are bounded too.
//!
//! Distances come from the input's BVH; the offset surface is met by sphere
//! tracing (distance is 1-Lipschitz, so a step of `distance − δ` never
//! passes it). Every combinatorial decision is the Delaunay's, with exact
//! predicates.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use mesher_geom::vec3::{add, cross, dot, norm, scale, sub};
use mesher_geom::{Aabb, Bvh, FastWindingNumber, GeomError, Point3, TriangleSoup};

use crate::delaunay::{Delaunay, DelaunayError, FACE_ORDER, GHOST, Inserted};

/// Wrapping parameters, in model units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WrapParams {
    /// Radius of the ball that must fit through a gate: openings narrower
    /// than about 2α are closed.
    pub alpha: f64,
    /// Distance δ of the wrap from the input; must be below α.
    pub offset: f64,
    /// At most this many inserted points (the wrap fails beyond).
    pub max_points: usize,
    /// Treat the points where the input's generalized winding number is at
    /// least 1/2 in magnitude as part of the input: the wrap then closes an
    /// open surface's holes across the winding number's level set instead
    /// of wrapping both sides of its walls.
    pub winding_guard: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum WrapError {
    /// The input has no triangle, or a non-finite coordinate.
    Input(String),
    /// α and δ must be positive and finite, with δ < α.
    Params { alpha: f64, offset: f64 },
    /// The tetrahedralization failed.
    Delaunay(DelaunayError),
    /// More than `max_points` points were needed.
    Budget { points: usize },
}

impl std::fmt::Display for WrapError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WrapError::Input(e) => write!(f, "the input cannot be wrapped: {e}"),
            WrapError::Params { alpha, offset } => write!(
                f,
                "alpha {alpha} and offset {offset} must be positive, with the offset below alpha"
            ),
            WrapError::Delaunay(e) => write!(f, "tetrahedralization failed: {e}"),
            WrapError::Budget { points } => {
                write!(f, "the wrap needed more than {points} points")
            }
        }
    }
}

impl std::error::Error for WrapError {}

impl From<GeomError> for WrapError {
    fn from(e: GeomError) -> Self {
        WrapError::Input(e.to_string())
    }
}

impl From<DelaunayError> for WrapError {
    fn from(e: DelaunayError) -> Self {
        WrapError::Delaunay(e)
    }
}

/// A wrap: a closed, oriented (outward), intersection-free surface around
/// the input.
#[derive(Debug, Clone, PartialEq)]
pub struct Wrapped {
    pub points: Vec<Point3>,
    pub triangles: Vec<[u32; 3]>,
    /// Points inserted by rule 1 (on Voronoi edges) and rule 2 (projected
    /// circumcentres).
    pub edge_points: usize,
    pub projected_points: usize,
    /// Points inserted where a Voronoi edge enters the winding guard's
    /// solid, and cells left inside because their circumcentre is in it.
    pub cap_points: usize,
    pub blocked: usize,
    /// Cells carved.
    pub carved: usize,
    /// Gates whose insertion left them in place (floating error at a
    /// tangency); the cell behind stays inside.
    pub stuck: usize,
    /// Non-manifold edges and vertices the carving left (each counted once
    /// per time it was found), and the outside cells relabelled inside to
    /// remove them.
    pub nonmanifold_edges: usize,
    pub nonmanifold_vertices: usize,
    pub relabelled: usize,
}

/// A gate in the queue: the inside cell and its face towards the outside.
#[derive(Debug, Clone, Copy)]
struct Gate {
    /// Radius of the smallest Delaunay ball through the facet.
    radius: f64,
    /// Sorted facet vertices, for deterministic ties.
    key: [u32; 3],
    cell: u32,
    generation: u32,
    face: usize,
}

impl PartialEq for Gate {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Gate {}

impl PartialOrd for Gate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Gate {
    /// Largest ball first; then the smaller facet vertices.
    fn cmp(&self, other: &Self) -> Ordering {
        self.radius
            .total_cmp(&other.radius)
            .then_with(|| other.key.cmp(&self.key))
            .then_with(|| other.cell.cmp(&self.cell))
            .then_with(|| other.face.cmp(&self.face))
    }
}

/// Wrap the triangle soup `points`/`triangles` (any soup: open,
/// self-intersecting, non-manifold).
pub fn wrap(
    points: &[Point3],
    triangles: &[[u32; 3]],
    params: WrapParams,
) -> Result<Wrapped, WrapError> {
    let WrapParams {
        alpha,
        offset,
        max_points,
        winding_guard,
    } = params;
    if !(alpha.is_finite() && offset.is_finite() && offset > 0.0 && offset < alpha) {
        return Err(WrapError::Params { alpha, offset });
    }
    if triangles.is_empty() {
        return Err(WrapError::Input("no triangles".into()));
    }
    let soup = TriangleSoup::new(points, triangles)?;
    let used: Vec<Point3> = triangles
        .iter()
        .flatten()
        .filter_map(|&v| points.get(v as usize).copied())
        .collect();
    if used.iter().flatten().any(|c| !c.is_finite()) {
        return Err(WrapError::Input("a coordinate is not finite".into()));
    }
    let bvh = Bvh::new(soup);
    let fwn = winding_guard.then(|| FastWindingNumber::new(soup));
    let bounds = Aabb::of_points(&used);
    let diagonal = bounds.diagonal();
    // The box: far enough out that no ball of radius α around a point of the
    // offset surface reaches it.
    let margin = 2.0 * (alpha + offset) + 0.1 * diagonal;
    let lo = bounds.min.map(|x| x - margin);
    let hi = bounds.max.map(|x| x + margin);
    let corners: Vec<Point3> = (0..8)
        .map(|k| {
            [
                if k & 1 == 0 { lo[0] } else { hi[0] },
                if k & 2 == 0 { lo[1] } else { hi[1] },
                if k & 4 == 0 { lo[2] } else { hi[2] },
            ]
        })
        .collect();
    let far = 4.0 * (norm(sub(hi, lo)) + margin);
    let mut wrapper = Wrapper {
        dt: Delaunay::new(corners, 0x5eed_0a1f)?,
        outside: Vec::new(),
        generation: Vec::new(),
        queue: BinaryHeap::new(),
        bvh: &bvh,
        fwn: fwn.as_ref(),
        alpha,
        offset,
        far,
        tolerance: 1e-4 * offset,
    };
    let n = wrapper.dt.slot_count();
    wrapper.outside = (0..n as u32).map(|s| wrapper.dt.is_ghost_slot(s)).collect();
    wrapper.generation = vec![0; n];
    for s in 0..n as u32 {
        if wrapper.dt.tet(s).is_some() {
            wrapper.push_gates(s);
        }
    }

    let mut out = Wrapped {
        points: Vec::new(),
        triangles: Vec::new(),
        edge_points: 0,
        projected_points: 0,
        cap_points: 0,
        blocked: 0,
        carved: 0,
        stuck: 0,
        nonmanifold_edges: 0,
        nonmanifold_vertices: 0,
        relabelled: 0,
    };
    while let Some(gate) = wrapper.queue.pop() {
        if !wrapper.valid(&gate) {
            continue;
        }
        let Some(step) = wrapper.traverse(&gate) else {
            continue;
        };
        match step {
            Step::Block => out.blocked += 1,
            Step::Carve => {
                wrapper.outside[gate.cell as usize] = true;
                out.carved += 1;
                wrapper.push_gates_around(gate.cell);
            }
            Step::Insert(p, rule) => {
                if out.edge_points + out.projected_points + out.cap_points >= max_points {
                    return Err(WrapError::Budget { points: max_points });
                }
                let (_, inserted) = wrapper.dt.insert(p)?;
                if inserted != Inserted::New {
                    out.stuck += 1;
                    continue;
                }
                match rule {
                    Rule::Edge => out.edge_points += 1,
                    Rule::Projection => out.projected_points += 1,
                    Rule::Cap => out.cap_points += 1,
                }
                wrapper.after_insertion();
                // Floating error at a tangency can leave the gate in place:
                // the cell behind it then stays inside.
                if wrapper.valid(&gate) {
                    out.stuck += 1;
                }
            }
        }
    }

    let (edges, vertices, relabelled) = wrapper.make_manifold();
    out.nonmanifold_edges = edges;
    out.nonmanifold_vertices = vertices;
    out.relabelled = relabelled;

    // The boundary of the outside region, oriented into it.
    let mut index = vec![u32::MAX; wrapper.dt.points().len()];
    for s in 0..wrapper.dt.slot_count() as u32 {
        let Some(v) = wrapper.dt.tet(s) else { continue };
        if wrapper.outside[s as usize] {
            continue;
        }
        let Some((_, nb)) = wrapper.dt.cell(s) else {
            continue;
        };
        for (j, &m) in nb.iter().enumerate() {
            if !wrapper.is_outside(m) {
                continue;
            }
            let f = FACE_ORDER[j].map(|k| v[k]);
            let t = f.map(|x| {
                if index[x as usize] == u32::MAX {
                    index[x as usize] = out.points.len() as u32;
                    out.points.push(wrapper.dt.points()[x as usize]);
                }
                index[x as usize]
            });
            out.triangles.push(t);
        }
    }
    Ok(out)
}

#[derive(Debug, Clone, Copy)]
enum Rule {
    Edge,
    Projection,
    Cap,
}

enum Step {
    Carve,
    /// The cell stays inside (its circumcentre is in the winding guard's
    /// solid).
    Block,
    Insert(Point3, Rule),
}

struct Wrapper<'a> {
    dt: Delaunay,
    /// Per slot: labelled outside (ghosts always).
    outside: Vec<bool>,
    /// Per slot: bumped whenever the slot holds a new tet.
    generation: Vec<u32>,
    queue: BinaryHeap<Gate>,
    bvh: &'a Bvh<'a>,
    fwn: Option<&'a FastWindingNumber<'a>>,
    alpha: f64,
    offset: f64,
    /// How far out a ghost's end of a Voronoi edge is.
    far: f64,
    /// Sphere tracing stops this close to the offset surface.
    tolerance: f64,
}

impl Wrapper<'_> {
    /// The facets between inside and outside cells, by their vertices.
    fn boundary(&self) -> Vec<[u32; 3]> {
        let mut out = Vec::new();
        for s in 0..self.dt.slot_count() as u32 {
            let Some(v) = self.dt.tet(s) else { continue };
            if self.outside[s as usize] {
                continue;
            }
            let Some((_, nb)) = self.dt.cell(s) else {
                continue;
            };
            for (j, &m) in nb.iter().enumerate() {
                if self.is_outside(m) {
                    out.push(FACE_ORDER[j].map(|k| v[k]));
                }
            }
        }
        out
    }

    /// Relabel outside cells inside until the boundary is a 2-manifold:
    /// around each edge and vertex of it the outside cells, and the inside
    /// ones, must form one group joined across faces through that edge or
    /// vertex. Of several outside groups the largest stays (one holding a
    /// ghost, if any: ghosts are never relabelled); a single outside group
    /// around a vertex whose inside cells are split is relabelled whole.
    /// The inside only grows, so the input stays enclosed; every pass
    /// relabels at least one cell, so this ends. Returns the non-manifold
    /// edges and vertices found and the cells relabelled.
    fn make_manifold(&mut self) -> (usize, usize, usize) {
        let (mut edges_found, mut vertices_found, mut relabelled) = (0, 0, 0);
        loop {
            let facets = self.boundary();
            let mut uses: std::collections::BTreeMap<(u32, u32), u32> =
                std::collections::BTreeMap::new();
            let mut vertices: Vec<u32> = Vec::new();
            for f in &facets {
                for k in 0..3 {
                    let (a, b) = (f[k], f[(k + 1) % 3]);
                    *uses.entry((a.min(b), a.max(b))).or_default() += 1;
                }
                vertices.extend(f);
            }
            vertices.sort_unstable();
            vertices.dedup();
            let mut changed = 0;
            for (&(u, v), &n) in &uses {
                if n > 2 {
                    edges_found += 1;
                    changed += self.fix(&self.dt.edge_tets(u, v), &[u, v]);
                }
            }
            if changed == 0 {
                for &v in &vertices {
                    let star = self.dt.star(v);
                    let r = self.fix(&star, &[v]);
                    if r > 0 {
                        vertices_found += 1;
                        changed += r;
                    }
                }
            }
            if changed == 0 {
                return (edges_found, vertices_found, relabelled);
            }
            relabelled += changed;
        }
    }

    /// Around the cells `cells` (those holding all of `through`), group the
    /// outside cells and the inside ones by shared faces holding `through`;
    /// relabel as described in `make_manifold`. Returns the cells
    /// relabelled.
    fn fix(&mut self, cells: &[u32], through: &[u32]) -> usize {
        let groups = |want_outside: bool| -> Vec<Vec<u32>> {
            let mut groups: Vec<Vec<u32>> = Vec::new();
            let mut seen: Vec<u32> = Vec::new();
            for &c in cells {
                if self.is_outside(c) != want_outside || seen.contains(&c) {
                    continue;
                }
                let mut group = vec![c];
                seen.push(c);
                let mut k = 0;
                while k < group.len() {
                    let t = group[k];
                    k += 1;
                    let Some((v, nb)) = self.dt.cell(t) else {
                        continue;
                    };
                    for j in 0..4 {
                        if through.contains(&v[j]) {
                            continue;
                        }
                        let m = nb[j];
                        if cells.contains(&m)
                            && self.is_outside(m) == want_outside
                            && !seen.contains(&m)
                        {
                            seen.push(m);
                            group.push(m);
                        }
                    }
                }
                groups.push(group);
            }
            groups
        };
        let outside = groups(true);
        let inside = groups(false);
        if outside.len() <= 1 && inside.len() <= 1 {
            return 0;
        }
        let ghostly = |g: &Vec<u32>| g.iter().any(|&c| self.dt.is_ghost_slot(c));
        let relabel: Vec<u32> = if outside.len() >= 2 {
            let keep = (0..outside.len())
                .max_by_key(|&i| (ghostly(&outside[i]), outside[i].len(), std::cmp::Reverse(i)))
                .unwrap_or(0);
            outside
                .iter()
                .enumerate()
                .filter(|&(i, g)| i != keep && !ghostly(g))
                .flat_map(|(_, g)| g.iter().copied())
                .collect()
        } else {
            outside
                .iter()
                .filter(|g| !ghostly(g))
                .flatten()
                .copied()
                .collect()
        };
        for &c in &relabel {
            self.outside[c as usize] = false;
        }
        relabel.len()
    }

    fn created(&self) -> Vec<u32> {
        self.dt.last_created_slots().to_vec()
    }

    fn is_outside(&self, slot: u32) -> bool {
        self.outside.get(slot as usize).copied().unwrap_or(true) || self.dt.is_ghost_slot(slot)
    }

    /// Whether the winding guard puts `p` in the input's solid.
    fn solid(&self, p: Point3) -> bool {
        self.fwn
            .is_some_and(|f| libm::fabs(f.winding_number(p)) >= 0.5)
    }

    /// The point where the segment from `a` (outside the guard's solid) to
    /// `b` (in it) enters the solid, by bisection on the winding number
    /// (continuous there: the segment is clear of the input).
    fn entry(&self, a: Point3, b: Point3) -> Point3 {
        let (mut lo, mut hi) = (a, b);
        for _ in 0..40 {
            let mid = scale(add(lo, hi), 0.5);
            if self.solid(mid) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        lo
    }

    fn distance(&self, p: Point3) -> f64 {
        self.bvh
            .nearest(p)
            .map_or(f64::INFINITY, |n| n.distance_squared.sqrt())
    }

    /// Labels and generations for the slots the last insertion created:
    /// real ones inside, ghosts outside; then their gates.
    fn after_insertion(&mut self) {
        let created = self.created();
        let n = self.dt.slot_count();
        self.outside.resize(n, false);
        self.generation.resize(n, 0);
        for &s in &created {
            self.outside[s as usize] = self.dt.is_ghost_slot(s);
            self.generation[s as usize] = self.generation[s as usize].wrapping_add(1);
        }
        for &s in &created {
            if self.dt.tet(s).is_some() {
                self.push_gates(s);
            }
        }
    }

    /// Whether a queued gate still is one: its inside cell unchanged and
    /// inside, its neighbour outside.
    fn valid(&self, g: &Gate) -> bool {
        let s = g.cell as usize;
        if self.dt.tet(g.cell).is_none()
            || self.generation.get(s) != Some(&g.generation)
            || self.outside.get(s).copied().unwrap_or(true)
        {
            return false;
        }
        self.dt
            .cell(g.cell)
            .is_some_and(|(_, nb)| self.is_outside(nb[g.face]))
    }

    /// Queue the gates of the inside cell in `slot`.
    fn push_gates(&mut self, slot: u32) {
        let Some((v, nb)) = self.dt.cell(slot) else {
            return;
        };
        if v[3] == GHOST || self.outside[slot as usize] {
            return;
        }
        for (j, &m) in nb.iter().enumerate() {
            if !self.is_outside(m) {
                continue;
            }
            if let Some(radius) = self.gate_radius(slot, j) {
                let mut key = FACE_ORDER[j].map(|k| v[k]);
                key.sort_unstable();
                self.queue.push(Gate {
                    radius,
                    key,
                    cell: slot,
                    generation: self.generation[slot as usize],
                    face: j,
                });
            }
        }
    }

    /// After carving `slot`, its inside neighbours get a gate towards it.
    fn push_gates_around(&mut self, slot: u32) {
        let Some((_, nb)) = self.dt.cell(slot) else {
            return;
        };
        for m in nb {
            if !self.is_outside(m) {
                self.push_gates(m);
            }
        }
    }

    /// The Voronoi edge dual to face `j` of the inside cell `slot`: its
    /// line (the facet's circumcentre and unit normal, pointing out of the
    /// cell), the facet's circumradius, and the parameters of the outside
    /// and inside ends along the normal.
    fn voronoi_edge(&self, slot: u32, j: usize) -> Option<VoronoiEdge> {
        let (v, nb) = self.dt.cell(slot)?;
        let pts = self.dt.points();
        let f = FACE_ORDER[j].map(|k| pts[v[k] as usize]);
        let (centre, radius) = triangle_circle(f)?;
        let normal = cross(sub(f[1], f[0]), sub(f[2], f[0]));
        let normal = scale(normal, 1.0 / norm(normal));
        let inside_centre = circumcentre(v.map(|x| pts[x as usize]))?;
        let inside_radius = norm(sub(pts[v[0] as usize], inside_centre));
        let t_in = dot(sub(inside_centre, centre), normal);
        let other = nb[j];
        let t_out = match self.dt.tet(other) {
            Some(w) => {
                let c = circumcentre(w.map(|x| pts[x as usize]))?;
                dot(sub(c, centre), normal)
            }
            // A ghost: the edge runs out to infinity along the outward normal.
            None => t_in.max(0.0) + self.far,
        };
        Some(VoronoiEdge {
            centre,
            normal,
            radius,
            t_out,
            t_in,
            inside_centre,
            inside_radius,
        })
    }

    /// The radius of the smallest Delaunay ball through the facet, if the
    /// facet is α-traversable.
    fn gate_radius(&self, slot: u32, j: usize) -> Option<f64> {
        let e = self.voronoi_edge(slot, j)?;
        let (lo, hi) = (e.t_out.min(e.t_in), e.t_out.max(e.t_in));
        let t = if lo <= 0.0 && 0.0 <= hi {
            0.0
        } else {
            lo.abs().min(hi.abs())
        };
        let r = (e.radius * e.radius + t * t).sqrt();
        (r >= self.alpha).then_some(r)
    }

    fn traverse(&self, g: &Gate) -> Option<Step> {
        let e = self.voronoi_edge(g.cell, g.face)?;
        let at = |t: f64| add(e.centre, scale(e.normal, t));
        // Rule 1: sphere-trace from the outside end towards the inside one.
        let length = e.t_out - e.t_in;
        let direction = if length >= 0.0 { -1.0 } else { 1.0 };
        let length = length.abs();
        let mut s = 0.0;
        let mut previous = at(e.t_out);
        for _ in 0..10_000 {
            let p = at(e.t_out + direction * s);
            let d = self.distance(p);
            if d <= self.offset + self.tolerance {
                return Some(Step::Insert(p, Rule::Edge));
            }
            // The guard: entering the solid between two samples (each step
            // stays in a ball clear of the input).
            if self.solid(p) {
                return Some(Step::Insert(self.entry(previous, p), Rule::Cap));
            }
            previous = p;
            s += d - self.offset;
            if s > length {
                let end = at(e.t_in);
                if self.solid(end) {
                    return Some(Step::Insert(self.entry(previous, end), Rule::Cap));
                }
                break;
            }
        }
        // Rule 2: the input enters the inside cell's circumball (its centre
        // within the offset is the paper's case).
        let c = e.inside_centre;
        let nearest = self.bvh.nearest(c)?;
        let d = nearest.distance_squared.sqrt();
        if self.solid(c) {
            return Some(Step::Block);
        }
        if d < e.inside_radius {
            let away = if d > 0.0 {
                scale(sub(c, nearest.point), 1.0 / d)
            } else {
                scale(e.normal, -direction)
            };
            let p = add(nearest.point, scale(away, self.offset));
            return Some(Step::Insert(p, Rule::Projection));
        }
        Some(Step::Carve)
    }
}

struct VoronoiEdge {
    centre: Point3,
    normal: Point3,
    radius: f64,
    t_out: f64,
    t_in: f64,
    inside_centre: Point3,
    inside_radius: f64,
}

/// The circumcentre and circumradius of a triangle, in its plane.
fn triangle_circle([a, b, c]: [Point3; 3]) -> Option<(Point3, f64)> {
    let (u, v) = (sub(b, a), sub(c, a));
    let n = cross(u, v);
    let nn = dot(n, n);
    if !(nn > 0.0 && nn.is_finite()) {
        return None;
    }
    let x = cross(n, u);
    let y = cross(v, n);
    let (uu, vv) = (dot(u, u), dot(v, v));
    let o = [0, 1, 2].map(|k| (uu * y[k] + vv * x[k]) / (2.0 * nn));
    Some((add(a, o), norm(o)))
}

/// The circumcentre of a tetrahedron (f64).
fn circumcentre([a, b, c, d]: [Point3; 4]) -> Option<Point3> {
    let (u, v, w) = (sub(b, a), sub(c, a), sub(d, a));
    let det = 2.0 * dot(u, cross(v, w));
    if !(det != 0.0 && det.is_finite()) {
        return None;
    }
    let num = add(
        add(scale(cross(v, w), dot(u, u)), scale(cross(w, u), dot(v, v))),
        scale(cross(u, v), dot(w, w)),
    );
    let o = scale(num, 1.0 / det);
    o.iter().all(|x| x.is_finite()).then(|| add(a, o))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circumcentres_are_equidistant() {
        let t = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ];
        let c = circumcentre(t).unwrap();
        let r: Vec<f64> = t.iter().map(|&p| norm(sub(p, c))).collect();
        for x in &r {
            assert!((x - r[0]).abs() < 1e-12);
        }
        let (o, r) = triangle_circle([t[0], t[1], t[2]]).unwrap();
        assert!((norm(sub(t[1], o)) - r).abs() < 1e-12);
        assert!((norm(sub(t[2], o)) - r).abs() < 1e-12);
    }
}
