//! The graded background grid of isosurface stuffing (Labelle and Shewchuk
//! 2007, Section 5): uniformly fine BCC tetrahedra wherever the surface
//! passes, grading to coarser tetrahedra in the interior, so the proven
//! dihedral-angle bounds still hold (only BCC tetrahedra are ever cut).
//!
//! 1. Leaf cubes: the finest cubes (width `w`) whose nine probe points (eight
//!    corners and the centre) have a nonpositive and a nonnegative sign,
//!    found by continuation from cubes around seed points on the surface.
//! 2. The Continuation Condition adds cubes across faces with a sign change
//!    and around corners whose sign differs from the centre's; a further
//!    step adds cubes where half-pyramids could be over-warped (a violated
//!    centre and two opposite violated corners on one face). Repeated until
//!    nothing changes.
//! 3. An octree of the leaf cubes and their ancestors, up to root cubes of
//!    width `w · 2^levels` tiling the box; an octant may have any subset of
//!    its eight children.
//! 4. Weak Balance: an octant that meets an edge shorter than half its width
//!    gets every child that meets the edge's interior.
//! 5. The paper's Figure 10 turns octants into BCC, bisected BCC,
//!    quadrisected BCC tetrahedra and half-pyramids.
//!
//! Coordinates are integers in units of `w / 2`, so octant corners and
//! centres at every level are exact.

use std::collections::{BTreeMap, BTreeSet};

use mesher_geom::predicates::orient3d;
use mesher_geom::vec3::{add, distance_squared, scale, sub};
use mesher_geom::{Aabb, Point3, Sign};
use rayon::prelude::*;

use crate::fxhash::{FxHashMap, FxHashSet};
use crate::stuffing::{Background, StuffingError, StuffingParams, sign_of};

type C = [i64; 3];

/// The graded background: vertices, their exact coordinates and signs, and
/// tetrahedra (positively oriented).
pub(crate) struct Graded {
    origin: Point3,
    unit: f64,
    coords: Vec<C>,
    signs: Vec<i8>,
    tets: Vec<[u32; 4]>,
}

impl Background for Graded {
    fn vertex_count(&self) -> usize {
        self.coords.len()
    }

    fn point(&self, id: u32) -> Point3 {
        point(self.origin, self.unit, self.coords[id as usize])
    }

    fn for_each_tet(&self, f: &mut dyn FnMut([u32; 4])) {
        for &t in &self.tets {
            f(t);
        }
    }

    fn edge(&self, u: u32, v: u32) -> (bool, f64) {
        let (a, b) = (self.coords[u as usize], self.coords[v as usize]);
        let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let axes = d.iter().filter(|&&x| x != 0).count();
        let length = (d.iter().map(|&x| (x * x) as f64).sum::<f64>()).sqrt() * self.unit;
        (axes == 1, length)
    }

    fn known_signs(&self) -> Option<Vec<i8>> {
        Some(self.signs.clone())
    }
}

fn point(origin: Point3, unit: f64, c: C) -> Point3 {
    [0, 1, 2].map(|k| origin[k] + unit * c[k] as f64)
}

/// Width of an octant at `level`, in units.
fn width(level: u32) -> i64 {
    2i64 << level
}

fn corner_of(level: u32, idx: C) -> C {
    idx.map(|i| i * width(level))
}

fn centre_of(level: u32, idx: C) -> C {
    let h = width(level) / 2;
    corner_of(level, idx).map(|x| x + h)
}

/// The eight corners, by bits (x = bit 0, y = bit 1, z = bit 2).
fn corners(level: u32, idx: C) -> [C; 8] {
    let c = corner_of(level, idx);
    let w = width(level);
    std::array::from_fn(|b| {
        [
            c[0] + w * (b as i64 & 1),
            c[1] + w * ((b as i64 >> 1) & 1),
            c[2] + w * ((b as i64 >> 2) & 1),
        ]
    })
}

fn shifted(c: C, axis: usize, by: i64) -> C {
    let mut c = c;
    c[axis] += by;
    c
}

/// Signs of the level function at lattice points, evaluated on demand in
/// parallel batches.
struct Signs<'a, F> {
    level: &'a F,
    origin: Point3,
    unit: f64,
    known: FxHashMap<C, i8>,
}

impl<F: Fn(Point3) -> f64 + Sync> Signs<'_, F> {
    fn ensure(&mut self, coords: impl IntoIterator<Item = C>) {
        let mut new: Vec<C> = coords
            .into_iter()
            .filter(|c| !self.known.contains_key(c))
            .collect();
        new.sort_unstable();
        new.dedup();
        let (origin, unit, level) = (self.origin, self.unit, self.level);
        let values: Vec<i8> = new
            .par_iter()
            .map(|&c| sign_of(level(point(origin, unit, c))))
            .collect();
        self.known.extend(new.into_iter().zip(values));
    }

    fn get(&self, c: C) -> i8 {
        self.known.get(&c).copied().unwrap_or(1)
    }
}

/// Probe points of a leaf cube: its eight corners and its centre.
fn probes(idx: C) -> impl Iterator<Item = C> {
    corners(0, idx)
        .into_iter()
        .chain(std::iter::once(centre_of(0, idx)))
}

/// Build the graded background for the region `level < 0` in `bounds`.
/// `params.cell_size` is the leaf width `w`; roots are `levels` levels up.
pub(crate) fn graded_background<F>(
    level: &F,
    bounds: &Aabb,
    seeds: &[Point3],
    params: &StuffingParams,
    levels: u32,
) -> Result<Graded, StuffingError>
where
    F: Fn(Point3) -> f64 + Sync,
{
    let w = params.cell_size;
    if !(w > 0.0 && w.is_finite()) {
        return Err(StuffingError::BadCellSize);
    }
    if bounds.is_empty() {
        return Err(StuffingError::EmptyBounds);
    }
    let unit = w / 2.0;
    let root = w * f64::from(1u32 << levels);
    let origin = bounds.min.map(|c| c - 1.5 * root);
    let roots: [i64; 3] = [0, 1, 2].map(|k| {
        let span = bounds.max[k] - bounds.min[k] + 3.0 * root;
        ((span / root).ceil() as i64).max(2)
    });
    let leaves_per_axis = roots.map(|n| n << levels);
    let in_range = |c: C| (0..3).all(|a| c[a] >= 0 && c[a] < leaves_per_axis[a]);
    let mut signs = Signs {
        level,
        origin,
        unit,
        known: FxHashMap::default(),
    };
    let budget = |n: usize| -> Result<(), StuffingError> {
        if n as u64 > params.max_lattice_vertices {
            Err(StuffingError::Budget {
                lattice_vertices: n as u64,
                limit: params.max_lattice_vertices,
            })
        } else {
            Ok(())
        }
    };

    // 1. Leaf cubes the surface passes through, by continuation from the
    //    cubes around the seeds.
    let leaf_of = |p: Point3| -> C { [0, 1, 2].map(|k| ((p[k] - origin[k]) / w).floor() as i64) };
    let mut visited: FxHashSet<C> = FxHashSet::default();
    let mut frontier: Vec<C> = Vec::new();
    for &p in seeds {
        let c = leaf_of(p);
        for d in 0..27 {
            let n = [c[0] + d % 3 - 1, c[1] + (d / 3) % 3 - 1, c[2] + d / 9 - 1];
            if in_range(n) && visited.insert(n) {
                frontier.push(n);
            }
        }
    }
    frontier.sort_unstable();
    let straddles = |signs: &Signs<'_, F>, idx: C| {
        let s: Vec<i8> = probes(idx).map(|c| signs.get(c)).collect();
        s.iter().any(|&x| x <= 0) && s.iter().any(|&x| x >= 0)
    };
    let mut leaves: BTreeSet<C> = BTreeSet::new();
    while !frontier.is_empty() {
        signs.ensure(frontier.iter().flat_map(|&c| probes(c)));
        budget(signs.known.len())?;
        let mut next = Vec::new();
        for &c in &frontier {
            if straddles(&signs, c) {
                leaves.insert(c);
                for a in 0..3 {
                    for by in [-1, 1] {
                        let n = shifted(c, a, by);
                        if in_range(n) && visited.insert(n) {
                            next.push(n);
                        }
                    }
                }
            }
        }
        next.sort_unstable();
        frontier = next;
    }

    // 2. The Continuation Condition and the half-pyramid step, until stable.
    loop {
        let before = leaves.len();
        continuation(&mut leaves, &mut signs, &in_range, params)?;
        half_pyramid_guard(&mut leaves, &mut signs, &in_range, params)?;
        if leaves.len() == before {
            break;
        }
    }

    // 3. The octree: leaves, their ancestors, and every root.
    let mut octants: Vec<FxHashSet<C>> = vec![FxHashSet::default(); levels as usize + 1];
    for &c in &leaves {
        let mut idx = c;
        for l in 0..=levels {
            octants[l as usize].insert(idx);
            idx = idx.map(|i| i >> 1);
        }
    }
    for k in 0..roots[2] {
        for j in 0..roots[1] {
            for i in 0..roots[0] {
                octants[levels as usize].insert([i, j, k]);
            }
        }
    }

    // 4. Weak Balance.
    weak_balance(&mut octants, levels, leaves_per_axis);

    // 5. Background tetrahedra (Figure 10).
    let mut vertices: FxHashSet<C> = FxHashSet::default();
    for (l, set) in octants.iter().enumerate() {
        for &idx in set {
            vertices.extend(corners(l as u32, idx));
            vertices.insert(centre_of(l as u32, idx));
        }
    }
    budget(vertices.len())?;
    let mut centres_needed: Vec<C> = Vec::new();
    for (l, set) in octants.iter().enumerate() {
        for &idx in set {
            centres_needed.push(centre_of(l as u32, idx));
        }
    }
    signs.ensure(centres_needed);
    let mut quads: BTreeSet<[C; 4]> = BTreeSet::new();
    for l in 0..=levels {
        let mut list: Vec<C> = octants[l as usize].iter().copied().collect();
        list.sort_unstable();
        for idx in list {
            figure_10(&octants, &vertices, &signs, levels, l, idx, &mut quads);
        }
    }

    // Number the vertices in coordinate order, orient, deduplicate.
    let mut used: Vec<C> = quads.iter().flatten().copied().collect();
    used.sort_unstable();
    used.dedup();
    signs.ensure(used.iter().copied());
    let id: BTreeMap<C, u32> = used
        .iter()
        .enumerate()
        .map(|(i, &c)| (c, i as u32))
        .collect();
    let mut tets = Vec::with_capacity(quads.len());
    for q in &quads {
        let mut t = q.map(|c| id[&c]);
        let p = q.map(|c| point(origin, unit, c));
        match orient3d(p[0], p[1], p[2], p[3]) {
            Sign::Positive => {}
            Sign::Negative => t.swap(0, 1),
            Sign::Zero => continue,
        }
        tets.push(t);
    }
    let signs_of = used.iter().map(|&c| signs.get(c)).collect();
    Ok(Graded {
        origin,
        unit,
        coords: used,
        signs: signs_of,
        tets,
    })
}

/// The Continuation Condition (Section 5), to a fixed point.
fn continuation<F: Fn(Point3) -> f64 + Sync>(
    leaves: &mut BTreeSet<C>,
    signs: &mut Signs<'_, F>,
    in_range: &dyn Fn(C) -> bool,
    params: &StuffingParams,
) -> Result<(), StuffingError> {
    let mut work: Vec<C> = leaves.iter().copied().collect();
    while !work.is_empty() {
        signs.ensure(work.iter().flat_map(|&c| probes(c)));
        if signs.known.len() as u64 > params.max_lattice_vertices {
            return Err(StuffingError::Budget {
                lattice_vertices: signs.known.len() as u64,
                limit: params.max_lattice_vertices,
            });
        }
        let mut add: BTreeSet<C> = BTreeSet::new();
        for &o in &work {
            let cs = corners(0, o);
            let s: [i8; 8] = cs.map(|c| signs.get(c));
            let centre = signs.get(centre_of(0, o));
            // A face with a nonpositive and a nonnegative corner.
            for a in 0..3 {
                for side in 0..2 {
                    let face: Vec<i8> = (0..8)
                        .filter(|&b| (b >> a) & 1 == side)
                        .map(|b| s[b])
                        .collect();
                    if face.iter().any(|&x| x <= 0) && face.iter().any(|&x| x >= 0) {
                        let n = shifted(o, a, if side == 1 { 1 } else { -1 });
                        if in_range(n) && !leaves.contains(&n) {
                            add.insert(n);
                        }
                    }
                }
            }
            // A corner whose sign differs from the centre's, or is zero.
            for (b, &sb) in s.iter().enumerate() {
                if sb * centre < 0 || sb == 0 || centre == 0 {
                    for a in 0..3 {
                        let n = shifted(o, a, if (b >> a) & 1 == 1 { 1 } else { -1 });
                        if in_range(n) && !leaves.contains(&n) {
                            add.insert(n);
                        }
                    }
                }
            }
        }
        work = add.iter().copied().collect();
        leaves.extend(add);
    }
    Ok(())
}

/// Add a leaf across a face where a leaf has a violated centre and two
/// diagonally opposite violated corners on that face, so the triangle two
/// half-pyramids share never has all three vertices warped.
fn half_pyramid_guard<F: Fn(Point3) -> f64 + Sync>(
    leaves: &mut BTreeSet<C>,
    signs: &mut Signs<'_, F>,
    in_range: &dyn Fn(C) -> bool,
    params: &StuffingParams,
) -> Result<(), StuffingError> {
    // The BCC edges among the leaves: centre to corner (short), cube edges
    // (long), centre to a neighbouring leaf's centre (long).
    let mut edges: BTreeSet<(C, C, bool)> = BTreeSet::new();
    for &o in leaves.iter() {
        let c = centre_of(0, o);
        let cs = corners(0, o);
        for &q in &cs {
            edges.insert((c.min(q), c.max(q), false));
        }
        for b in 0..8usize {
            for a in 0..3 {
                if (b >> a) & 1 == 0 {
                    let (p, q) = (cs[b], cs[b | (1 << a)]);
                    edges.insert((p, q, true));
                }
            }
        }
        for a in 0..3 {
            let n = shifted(o, a, 1);
            if leaves.contains(&n) {
                edges.insert((c, centre_of(0, n), true));
            }
        }
    }
    let cut: Vec<(C, C, bool)> = edges
        .into_iter()
        .filter(|&(p, q, _)| signs.get(p) * signs.get(q) < 0)
        .collect();
    let (origin, unit, level) = (signs.origin, signs.unit, signs.level);
    let steps = params.bisection_steps;
    let violated: Vec<Vec<C>> = cut
        .par_iter()
        .map(|&(p, q, long)| {
            let (pp, pq) = (point(origin, unit, p), point(origin, unit, q));
            let inside_p = level(pp) < 0.0;
            let (a, b) = if inside_p { (pp, pq) } else { (pq, pp) };
            let (mut lo, mut hi) = (0.0f64, 1.0f64);
            for _ in 0..steps {
                let mid = 0.5 * (lo + hi);
                if level(add(a, scale(sub(b, a), mid))) < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let x = add(a, scale(sub(b, a), 0.5 * (lo + hi)));
            let length = distance_squared(pp, pq).sqrt();
            let alpha = if long {
                params.alpha_long
            } else {
                params.alpha_short
            };
            let limit = (alpha * length).powi(2);
            let mut out = Vec::new();
            if distance_squared(x, pp) < limit {
                out.push(p);
            }
            if distance_squared(x, pq) < limit {
                out.push(q);
            }
            out
        })
        .collect();
    let violated: FxHashSet<C> = violated.into_iter().flatten().collect();
    let mut add: BTreeSet<C> = BTreeSet::new();
    for &o in leaves.iter() {
        if !violated.contains(&centre_of(0, o)) {
            continue;
        }
        let cs = corners(0, o);
        for a in 0..3 {
            for side in 0..2usize {
                // The face's corners in order around it.
                let (b1, b2) = ((a + 1) % 3, (a + 2) % 3);
                let at = |u: usize, v: usize| cs[(side << a) | (u << b1) | (v << b2)];
                let ring = [at(0, 0), at(1, 0), at(1, 1), at(0, 1)];
                let opposite = (violated.contains(&ring[0]) && violated.contains(&ring[2]))
                    || (violated.contains(&ring[1]) && violated.contains(&ring[3]));
                if opposite {
                    let n = shifted(o, a, if side == 1 { 1 } else { -1 });
                    if in_range(n) && !leaves.contains(&n) {
                        add.insert(n);
                    }
                }
            }
        }
    }
    leaves.extend(add);
    Ok(())
}

/// Weak Balance: an octant that meets an edge of the octree shorter than
/// half its width gets every child that meets the edge's interior.
fn weak_balance(octants: &mut [FxHashSet<C>], levels: u32, leaves_per_axis: [i64; 3]) {
    let range = |l: u32, idx: C| (0..3).all(|a| idx[a] >= 0 && idx[a] < (leaves_per_axis[a] >> l));
    let mut work: Vec<(u32, C)> = Vec::new();
    for l in 0..=levels {
        let mut list: Vec<C> = octants[l as usize].iter().copied().collect();
        list.sort_unstable();
        work.extend(list.into_iter().map(|c| (l, c)));
    }
    while let Some((l, q)) = work.pop() {
        let cs = corners(l, q);
        // The 12 edges: corner b to corner b + axis a.
        for (b, &p0) in cs.iter().enumerate() {
            for a in 0..3 {
                if (b >> a) & 1 == 1 {
                    continue;
                }
                for m in l + 2..=levels {
                    let wm = width(m);
                    // Level-m cubes meeting the closed edge.
                    let along = p0[a].div_euclid(wm);
                    let options = |x: i64, w: i64| -> Vec<i64> {
                        if x.rem_euclid(w) == 0 {
                            vec![x / w - 1, x / w]
                        } else {
                            vec![x.div_euclid(w)]
                        }
                    };
                    let (b1, b2) = ((a + 1) % 3, (a + 2) % 3);
                    for i1 in options(p0[b1], wm) {
                        for i2 in options(p0[b2], wm) {
                            let mut o = [0i64; 3];
                            o[a] = along;
                            o[b1] = i1;
                            o[b2] = i2;
                            if !range(m, o) || !octants[m as usize].contains(&o) {
                                continue;
                            }
                            // Children meeting the edge's interior.
                            let wc = width(m - 1);
                            let child_along = p0[a].div_euclid(wc);
                            for j1 in options(p0[b1], wc) {
                                for j2 in options(p0[b2], wc) {
                                    let mut ch = [0i64; 3];
                                    ch[a] = child_along;
                                    ch[b1] = j1;
                                    ch[b2] = j2;
                                    if ch.map(|x| x >> 1) != o {
                                        continue;
                                    }
                                    if octants[(m - 1) as usize].insert(ch) {
                                        work.push((m - 1, ch));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Figure 10 of the paper, for one octant: the background tetrahedra it is
/// responsible for, as coordinate quadruples (sorted, so a tetrahedron two
/// octants both create is kept once).
fn figure_10<F: Fn(Point3) -> f64 + Sync>(
    octants: &[FxHashSet<C>],
    vertices: &FxHashSet<C>,
    signs: &Signs<'_, F>,
    levels: u32,
    l: u32,
    idx: C,
    out: &mut BTreeSet<[C; 4]>,
) {
    let has_child = |b: usize| -> bool {
        l > 0 && {
            let ch = [0, 1, 2].map(|a| 2 * idx[a] + ((b >> a) & 1) as i64);
            octants[(l - 1) as usize].contains(&ch)
        }
    };
    let leaf = (0..8).all(|b| !has_child(b));
    let c = centre_of(l, idx);
    // Inside or on (the paper's "nonnegative": its sign convention is ours
    // reversed).
    if !(leaf || signs.get(c) <= 0) {
        return;
    }
    let mut push = |mut q: [C; 4]| {
        q.sort_unstable();
        out.insert(q);
    };
    let half = width(l) / 2;
    let cs = corners(l, idx);
    for a in 0..3 {
        for side in 0..2usize {
            let (b1, b2) = ((a + 1) % 3, (a + 2) % 3);
            let bit = |u: usize, v: usize| (side << a) | (u << b1) | (v << b2);
            let ring_bits = [bit(0, 0), bit(1, 0), bit(1, 1), bit(0, 1)];
            let ring = ring_bits.map(|b| cs[b]);
            let face_centre = shifted(c, a, if side == 1 { half } else { -half });
            let mid = |p: C, q: C| [0, 1, 2].map(|k| (p[k] + q[k]) / 2);
            if !vertices.contains(&face_centre) {
                let n = shifted(idx, a, if side == 1 { 1 } else { -1 });
                if octants[l as usize].contains(&n) {
                    let c2 = centre_of(l, n);
                    for e in 0..4 {
                        let (p, q) = (ring[e], ring[(e + 1) % 4]);
                        let m = mid(p, q);
                        if vertices.contains(&m) {
                            push([p, m, c, c2]);
                            push([q, m, c, c2]);
                        } else {
                            push([p, q, c, c2]);
                        }
                    }
                } else {
                    // Shared with a larger octant or the boundary: two
                    // half-pyramids, split by the diagonal through a corner
                    // or the centre of the parent.
                    let parent: Option<(C, [C; 8])> = (l < levels).then(|| {
                        let pi = idx.map(|x| x >> 1);
                        (centre_of(l + 1, pi), corners(l + 1, pi))
                    });
                    let on_parent =
                        |v: C| parent.is_some_and(|(pc, pcs)| v == pc || pcs.contains(&v));
                    let k = (0..4).find(|&k| on_parent(ring[k])).unwrap_or(0);
                    let (d1, d2) = (ring[k], ring[(k + 2) % 4]);
                    push([d1, d2, ring[(k + 1) % 4], c]);
                    push([d1, d2, ring[(k + 3) % 4], c]);
                }
            } else {
                let d = face_centre;
                for e in 0..4 {
                    let (p, q) = (ring[e], ring[(e + 1) % 4]);
                    let m = mid(p, q);
                    if vertices.contains(&m) {
                        for (end, b) in [(p, ring_bits[e]), (q, ring_bits[(e + 1) % 4])] {
                            if !has_child(b) {
                                push([end, m, d, c]);
                            }
                        }
                    } else {
                        push([p, q, d, c]);
                    }
                }
            }
        }
    }
}
