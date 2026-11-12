//! Untangling a tetrahedralization after its implicit points are rounded to
//! f64 (Milestone 13).
//!
//! Rounding moves every constructed point by up to an ulp or so, which can
//! turn a thin tet over. Every orientation here is exact on the rounded
//! coordinates. A tet inside the region that is not usable (not positive
//! exactly, or with no positive volume in floating point, which a
//! finite-volume solver computes) is fixed, in
//! order of preference, by moving one of its constructed vertices to a
//! nearby f64 point (within a few ulps of the exact position, then up to
//! [`FAR_ULPS`], so the boundary moves by rounding only), by a flip that
//! removes it (a 2-3 flip across a face, or edge removal), never across a
//! boundary face or through a protected edge, or by moving two of its
//! constructed vertices together. Each fix is returned, for the report.
//!
//! The far and paired moves are for slivers between near-coincident walls
//! (two facets of the input less than an ulp apart, with flat tets
//! between them): no single vertex can move within a few ulps without
//! turning over another of its tets.

use mesher_geom::Sign;
use mesher_geom::implicit::Point;

use super::{EDGES, ExactDelaunay, GHOST, face, sorted3};

/// Ulps either way a coordinate may move from its rounded value, trying
/// every offset.
const MOVE_ULPS: i64 = 2;

/// The farthest a coordinate may move, in ulps: beyond [`MOVE_ULPS`], each
/// coordinate moves by `-r`, 0 or `r` for `r` doubling up to this. At 64
/// ulps a point moves by under 1e-14 of its distance from the origin.
const FAR_ULPS: i64 = 64;

/// Unusable tets per untangling that may try paired moves, which weigh
/// tens of thousands of positions each.
const PAIR_TETS: usize = 32;

/// Passes over the tets still inverted.
const PASSES: usize = 8;

/// What untangling came to.
#[derive(Debug, Clone, Default)]
pub(crate) struct Untangled {
    /// Inside tets turned over (not positive exactly) before.
    pub inverted: usize,
    /// Inside tets positive exactly but with no positive volume in floating
    /// point, before.
    pub flat: usize,
    pub fixes: Vec<Untangle>,
    /// Inside tets still not usable, by vertices.
    pub left: Vec<[u32; 4]>,
}

/// One fix of a tet rounding turned over or flattened.
#[derive(Debug, Clone, PartialEq)]
pub enum Untangle {
    /// A constructed vertex moved, within a few ulps of its exact position.
    Move {
        vertex: u32,
        from: [f64; 3],
        to: [f64; 3],
    },
    /// The tet was flipped away: a 2-3 flip across face `face`.
    Flip23 { tet: [u32; 4], face: [u32; 3] },
    /// The tet was flipped away by removing edge `edge`.
    RemoveEdge { tet: [u32; 4], edge: [u32; 2] },
}

/// The next f64 above `x` (finite `x`; through zero; the MSRV predates
/// `f64::next_up`).
fn next_up(x: f64) -> f64 {
    if !x.is_finite() || x == f64::MAX {
        return x;
    }
    if x == 0.0 {
        return f64::from_bits(1);
    }
    let bits = x.to_bits();
    f64::from_bits(if x > 0.0 { bits + 1 } else { bits - 1 })
}

/// Offsets in ulps per coordinate, nearest first: every one within `near`,
/// then each coordinate `-r`, 0 or `r` for `r` doubling from `2 * near` up
/// to [`FAR_ULPS`]. A fixed order: deterministic.
fn offsets(near: i64) -> Vec<[i64; 3]> {
    let norm = |o: &[i64; 3]| o.iter().map(|x| x * x).sum::<i64>();
    let mut all: Vec<[i64; 3]> = Vec::new();
    for a in -near..=near {
        for b in -near..=near {
            for c in -near..=near {
                all.push([a, b, c]);
            }
        }
    }
    all.sort_by_key(|o| (norm(o), *o));
    let mut r = 2 * near;
    while r <= FAR_ULPS {
        let mut ring: Vec<[i64; 3]> = Vec::new();
        for a in [-r, 0, r] {
            for b in [-r, 0, r] {
                for c in [-r, 0, r] {
                    if [a, b, c] != [0, 0, 0] {
                        ring.push([a, b, c]);
                    }
                }
            }
        }
        ring.sort_by_key(|o| (norm(o), *o));
        all.extend(ring);
        r *= 2;
    }
    all
}

/// The f64 `k` steps from `x` (`k` may be negative), through zero.
fn step(x: f64, k: i64) -> f64 {
    let mut y = x;
    for _ in 0..k.unsigned_abs() {
        y = if k > 0 { next_up(y) } else { -next_up(-y) };
    }
    y
}

impl ExactDelaunay {
    /// Replace point `i` (a vertex) by `p`. The structure is unchanged; tets
    /// may then be turned over.
    pub(crate) fn set_point(&mut self, i: u32, p: Point) {
        if let Some(slot) = self.points.get_mut(i as usize) {
            *slot = p;
        }
    }

    /// Which live tets (by index) are inside: a parity flood from the
    /// ghosts that toggles at every wall face; `None` where labels disagree.
    pub(crate) fn inside_tets(&self, wall: &dyn Fn([u32; 3]) -> bool) -> Option<Vec<bool>> {
        let mut label: Vec<Option<bool>> = vec![None; self.tets.len()];
        let mut queue: std::collections::VecDeque<u32> = std::collections::VecDeque::new();
        for (t, tet) in self.tets.iter().enumerate() {
            if tet.alive && tet.v[3] == GHOST {
                label[t] = Some(false);
                queue.push_back(t as u32);
            }
        }
        while let Some(t) = queue.pop_front() {
            let tet = self.tets[t as usize];
            let here = label[t as usize].unwrap_or(false);
            for j in 0..4 {
                let u = tet.n[j];
                let f = sorted3(face(&tet.v, j));
                let want = here ^ (!f.contains(&GHOST) && wall(f));
                match label[u as usize] {
                    None => {
                        label[u as usize] = Some(want);
                        queue.push_back(u);
                    }
                    Some(have) if have != want => return None,
                    Some(_) => {}
                }
            }
        }
        Some(
            self.tets
                .iter()
                .zip(&label)
                .map(|(t, l)| t.alive && t.v[3] != GHOST && *l == Some(true))
                .collect(),
        )
    }

    /// Is the tet positive exactly and in floating point (the volume a
    /// finite-volume solver computes, from its vertices in the best of the
    /// orders that keep its orientation; that order is the one written)?
    fn usable(&self, v: [u32; 4]) -> bool {
        self.orient4(v) == Sign::Positive
            && crate::quality::best_float_order(v, &|x| self.points[x as usize].approx()).1 > 0.0
    }

    /// The inside tets (by index) that are not usable.
    fn unusable_inside(&self, inside: &[bool]) -> Vec<u32> {
        (0..self.tets.len() as u32)
            .filter(|&t| inside[t as usize] && !self.usable(self.tets[t as usize].v))
            .collect()
    }

    /// Unusable real tets around `v`: inside ones, and the rest.
    fn bad_around(&self, v: u32, inside: &[bool]) -> (usize, usize) {
        let mut bad = (0, 0);
        for t in self.star(v) {
            if self.is_ghost(t) || self.usable(self.tets[t as usize].v) {
                continue;
            }
            if inside[t as usize] {
                bad.0 += 1;
            } else {
                bad.1 += 1;
            }
        }
        bad
    }

    /// Move constructed vertex `v` to the f64 point near `exact` (by
    /// [`offsets`] of its rounding) that leaves no inside tet around it
    /// turned over and no more outside ones than now, nearest first.
    /// Returns where it went.
    fn try_move(&mut self, v: u32, exact: &Point, inside: &[bool]) -> Option<([f64; 3], [f64; 3])> {
        let from = self.points[v as usize].approx();
        let base = exact.approx();
        let (_, outside_now) = self.bad_around(v, inside);
        for o in offsets(MOVE_ULPS) {
            let to = [0, 1, 2].map(|k| step(base[k], o[k]));
            if to == from || !to.iter().all(|x| x.is_finite()) {
                continue;
            }
            let Ok(p) = Point::explicit(to) else {
                continue;
            };
            self.set_point(v, p);
            let (ins, outs) = self.bad_around(v, inside);
            if ins == 0 && outs <= outside_now {
                return Some((from, to));
            }
        }
        if let Ok(p) = Point::explicit(from) {
            self.set_point(v, p);
        }
        None
    }

    /// Move constructed vertices `a` and `b` together (each by
    /// [`offsets`] within one ulp, then doubling) so that no inside tet
    /// around either is unusable and no more outside ones are than now.
    /// Returns the two moves.
    fn try_pair(
        &mut self,
        [a, b]: [u32; 2],
        exact: &[Point],
        inside: &[bool],
    ) -> Option<[Untangle; 2]> {
        let (from_a, from_b) = (self.points[a as usize], self.points[b as usize]);
        let (base_a, base_b) = (exact[a as usize].approx(), exact[b as usize].approx());
        let (_, outside_a) = self.bad_around(a, inside);
        let (_, outside_b) = self.bad_around(b, inside);
        let offsets = offsets(1);
        for oa in &offsets {
            let to_a = [0, 1, 2].map(|k| step(base_a[k], oa[k]));
            let Ok(pa) = Point::explicit(to_a) else {
                continue;
            };
            self.set_point(a, pa);
            for ob in &offsets {
                let to_b = [0, 1, 2].map(|k| step(base_b[k], ob[k]));
                let Ok(pb) = Point::explicit(to_b) else {
                    continue;
                };
                self.set_point(b, pb);
                let (in_a, out_a) = self.bad_around(a, inside);
                let (in_b, out_b) = self.bad_around(b, inside);
                if in_a == 0 && in_b == 0 && out_a <= outside_a && out_b <= outside_b {
                    return Some([
                        Untangle::Move {
                            vertex: a,
                            from: from_a.approx(),
                            to: to_a,
                        },
                        Untangle::Move {
                            vertex: b,
                            from: from_b.approx(),
                            to: to_b,
                        },
                    ]);
                }
            }
        }
        self.set_point(a, from_a);
        self.set_point(b, from_b);
        None
    }

    /// Fix the inside tets that rounding left unusable, by moving
    /// constructed vertices (`movable`, with their exact positions in
    /// `exact`) and by flips that never remove a `wall` face or a `protected`
    /// edge. `None` if the walls stop bounding a region.
    pub(crate) fn untangle(
        &mut self,
        movable: &[bool],
        exact: &[Point],
        wall: &dyn Fn([u32; 3]) -> bool,
        protected: &dyn Fn(u32, u32) -> bool,
    ) -> Option<Untangled> {
        let mut inside = self.inside_tets(wall)?;
        let before = self.unusable_inside(&inside);
        let inverted = before
            .iter()
            .filter(|&&t| self.orient4(self.tets[t as usize].v) != Sign::Positive)
            .count();
        let flat = before.len() - inverted;
        let mut fixes = Vec::new();
        let mut paired = 0;
        // Flips must give tets a solver can use, not just positive ones.
        self.float_positive = true;
        for _ in 0..PASSES {
            let bad = self.unusable_inside(&inside);
            if bad.is_empty() {
                break;
            }
            let mut progress = false;
            for t in bad {
                let tet = self.tets[t as usize];
                if !tet.alive || !inside[t as usize] || self.usable(tet.v) {
                    continue;
                }
                // 1. Move a constructed vertex within its rounding.
                let mut fixed = false;
                for &v in &tet.v {
                    if !movable.get(v as usize).copied().unwrap_or(false) {
                        continue;
                    }
                    if let Some((from, to)) = self.try_move(v, &exact[v as usize], &inside) {
                        fixes.push(Untangle::Move {
                            vertex: v,
                            from,
                            to,
                        });
                        fixed = true;
                        break;
                    }
                }
                // 2. Flip it away: across a face, then by removing an edge.
                if !fixed {
                    for j in 0..4 {
                        let f = face(&tet.v, j);
                        if self.is_ghost(tet.n[j]) || wall(sorted3(f)) {
                            continue;
                        }
                        if self.flip23(t, j) {
                            fixes.push(Untangle::Flip23 {
                                tet: tet.v,
                                face: f,
                            });
                            fixed = true;
                            break;
                        }
                    }
                }
                if !fixed {
                    for (a, b) in EDGES {
                        let (x, y) = (tet.v[a], tet.v[b]);
                        if protected(x, y) {
                            continue;
                        }
                        if self.remove_edge(x, y, protected) {
                            fixes.push(Untangle::RemoveEdge {
                                tet: tet.v,
                                edge: [x, y],
                            });
                            fixed = true;
                            break;
                        }
                    }
                }
                // 3. Move two constructed vertices together.
                if !fixed && paired < PAIR_TETS {
                    paired += 1;
                    let mine: Vec<u32> = tet
                        .v
                        .into_iter()
                        .filter(|&v| movable.get(v as usize).copied().unwrap_or(false))
                        .collect();
                    'pairs: for (i, &x) in mine.iter().enumerate() {
                        for &y in &mine[i + 1..] {
                            if let Some(moves) = self.try_pair([x, y], exact, &inside) {
                                fixes.extend(moves);
                                fixed = true;
                                break 'pairs;
                            }
                        }
                    }
                }
                if fixed {
                    progress = true;
                    // Flips replace tets on one side of the walls; label anew.
                    let Some(now) = self.inside_tets(wall) else {
                        self.float_positive = false;
                        return None;
                    };
                    inside = now;
                }
            }
            if !progress {
                break;
            }
        }
        self.float_positive = false;
        let left = self
            .unusable_inside(&inside)
            .into_iter()
            .map(|t| self.tets[t as usize].v)
            .collect();
        Some(Untangled {
            inverted,
            flat,
            fixes,
            left,
        })
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use mesher_geom::implicit::Point;

    use super::super::ExactDelaunay;
    use super::{FAR_ULPS, MOVE_ULPS, Untangle, offsets, step};

    /// A cube's corners and a constructed point at its middle (an LNC point
    /// on the diagonal): every tet joins the middle to a hull triangle.
    /// Knocked far off, the middle turns tets over; untangling moves it back
    /// within its rounding.
    #[test]
    fn a_vertex_knocked_off_is_moved_back() {
        let mut p: Vec<Point> = (0..8)
            .map(|i| {
                Point::explicit([
                    f64::from(i & 1),
                    f64::from((i >> 1) & 1),
                    f64::from((i >> 2) & 1),
                ])
                .unwrap()
            })
            .collect();
        let middle = Point::lnc([0.0; 3], [1.0; 3], 0.5).unwrap();
        p.push(middle);
        let mut d = ExactDelaunay::new(p.clone(), 3).unwrap();
        let m = 8u32;
        d.set_point(m, Point::explicit([1.5, 0.5, 0.5]).unwrap());
        let wall = |f: [u32; 3]| !f.contains(&m);
        let movable: Vec<bool> = (0..9).map(|i| i == 8).collect();
        let u = d.untangle(&movable, &p, &wall, &|_, _| false).unwrap();
        assert!(u.inverted > 0, "the knock should turn tets over");
        assert!(u.left.is_empty(), "{:?}", u.left);
        assert!(matches!(u.fixes[0], Untangle::Move { vertex: 8, .. }));
        d.check().unwrap();
    }

    /// Every offset within the near cube comes first, nearest first, then
    /// rings doubling out to the farthest; none repeats.
    #[test]
    fn offsets_go_from_near_to_far() {
        let o = offsets(MOVE_ULPS);
        let near = (2 * MOVE_ULPS + 1).pow(3) as usize;
        assert_eq!(o[0], [0, 0, 0]);
        assert!(
            o[..near]
                .iter()
                .all(|v| v.iter().all(|x| x.abs() <= MOVE_ULPS))
        );
        assert!(
            o[near..]
                .iter()
                .all(|v| v.iter().any(|x| x.abs() > MOVE_ULPS))
        );
        let norm = |v: &[i64; 3]| v.iter().map(|x| x * x).sum::<i64>();
        assert!(o[..near].windows(2).all(|w| norm(&w[0]) <= norm(&w[1])));
        assert_eq!(
            o.last().map(|v| v.iter().map(|x| x.abs()).max()),
            Some(Some(FAR_ULPS))
        );
        let mut sorted = o.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), o.len());
    }

    #[test]
    fn steps_walk_the_f64_lattice() {
        assert_eq!(step(1.0, 0), 1.0);
        assert_eq!(step(1.0, 1), 1.0 + f64::EPSILON);
        assert_eq!(step(1.0, -1), 1.0 - f64::EPSILON / 2.0);
        assert_eq!(step(step(0.0, 1), -2), -f64::from_bits(1));
    }
}
