//! Gift-wrapping a facet cavity, where filling it by Delaunay
//! tetrahedralization and expansion fails (Diazzi et al. 2023, Section 4.4,
//! after Shewchuk's gift-wrapping).
//!
//! Each half of the cavity is filled one tet at a time from a front of
//! faces, the boundary of what is still empty: a front face takes the apex,
//! among the half's vertices on its inner side, whose tet overlaps nothing
//! already there and whose circumsphere holds no other such vertex. The
//! half filled first leaves the plane open, and the faces its tets put on
//! the plane close the other half. Ties between cospherical points are
//! broken by symbolic perturbation of the lifting (after Edelsbrunner and
//! Mücke), by vertex index, so both halves see the same order. No points
//! are added; where a front face has no valid apex the fill fails.

use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};

use mesher_geom::Sign;
use mesher_geom::implicit;

use super::{Bounds, ExactDelaunay, apart_boxes, face, join, rotated, sorted3};
use crate::fxhash::FxHashMap;

/// Tets per cavity tet and vertex a fill may make before it is abandoned.
const WRAP_STEPS: usize = 16;

/// Orientations of vertex quadruples, remembered: a fill asks the same ones
/// of many candidate tets. Also each cavity vertex's box
/// ([`mesher_geom::implicit::Point::bounds`]), to pass over vertices and
/// front faces that certainly miss a candidate tet.
struct Orients<'a> {
    dt: &'a ExactDelaunay,
    seen: RefCell<FxHashMap<[u32; 4], Sign>>,
    /// The cavity's vertices, sorted, and their boxes.
    verts: Vec<u32>,
    boxes: Vec<Bounds>,
}

impl Orients<'_> {
    fn bounds(&self, v: u32) -> Bounds {
        match self.verts.binary_search(&v) {
            Ok(k) => self.boxes[k],
            Err(_) => self.dt.pt(v).bounds(),
        }
    }

    fn bounds_of(&self, vs: &[u32]) -> Bounds {
        vs.iter()
            .map(|&v| self.bounds(v))
            .reduce(|a, b| join(&a, &b))
            .unwrap_or([[f64::NEG_INFINITY, f64::INFINITY]; 3])
    }

    fn orient(&self, v: [u32; 4]) -> Sign {
        // Sort, counting swaps: an odd permutation flips the sign.
        let mut k = v;
        let mut odd = false;
        for i in 0..4 {
            for j in 0..3 - i {
                if k[j] > k[j + 1] {
                    k.swap(j, j + 1);
                    odd = !odd;
                }
            }
        }
        if k[0] == k[1] || k[1] == k[2] || k[2] == k[3] {
            return Sign::Zero;
        }
        let known = self.seen.borrow().get(&k).copied();
        let s = known.unwrap_or_else(|| {
            let s = self.dt.orient4(k);
            self.seen.borrow_mut().insert(k, s);
            s
        });
        if odd { s.flipped() } else { s }
    }
}

impl ExactDelaunay {
    /// Is `v[4]` inside the circumsphere of tet `v[0..4]`, with every lift
    /// lowered by an infinitesimal that shrinks with vertex index? Positive
    /// inside, for a positive tet; never zero for five distinct points not
    /// all coplanar.
    pub(super) fn insphere_perturbed(&self, v: [u32; 5]) -> Sign {
        let [a, b, c, d, e] = v.map(|x| self.pt(x));
        let s = implicit::insphere(a, b, c, d, e);
        if s != Sign::Zero {
            return s;
        }
        // The determinant alternates with the order of its rows, so lowering
        // the lift of the point in slot k counts as that point being the
        // query: +orient3d of the other four, after swapping slot k into the
        // last place (a sign change unless k is last).
        let mut by_index: Vec<usize> = (0..5).collect();
        by_index.sort_unstable_by_key(|&k| v[k]);
        for k in by_index {
            let mut w = v;
            w.swap(k, 4);
            let o = self.orient4([w[0], w[1], w[2], w[3]]);
            if o != Sign::Zero {
                return if k == 4 { o } else { o.flipped() };
            }
        }
        Sign::Zero
    }

    /// Do the open segments `ab` and `cd`, in one plane, cross at a point
    /// inside both? Exact.
    fn segments_cross(&self, o4: &Orients, a: u32, b: u32, c: u32, d: u32) -> bool {
        if o4.orient([a, b, c, d]) != Sign::Zero {
            return false;
        }
        let [pa, pb, pc, pd] = [a, b, c, d].map(|x| self.pt(x));
        let Some(axis) =
            (0..3).find(|&k| implicit::orient2d_projected(pa, pb, pc, k) != Sign::Zero)
        else {
            return false;
        };
        let o = |p, q, r| implicit::orient2d_projected(p, q, r, axis);
        let opposite = |x: Sign, y: Sign| x != Sign::Zero && y != Sign::Zero && x != y;
        opposite(o(pa, pb, pc), o(pa, pb, pd)) && opposite(o(pc, pd, pa), o(pc, pd, pb))
    }

    /// Can positive tet `t`, on front face `face(t, 3)`, go into the empty
    /// part bounded by `front` (faces outward, keyed by sorted vertices)?
    /// No vertex of `verts` other than its own lies in it or on it, a face
    /// it shares with the front is matched from the inside, and no edge or
    /// face of it meets a front face or edge other than in a common vertex,
    /// edge or face.
    fn wrap_fits(
        &self,
        o4: &Orients,
        t: [u32; 4],
        front: &BTreeMap<[u32; 3], [u32; 3]>,
        verts: &[u32],
        below: &[Sign],
    ) -> bool {
        let _p = crate::prof::scope("cdt wrap: fits");
        let faces = [0, 1, 2, 3].map(|k| face(&t, k));
        let tet_box = o4.bounds_of(&t);
        // Each vertex's side of the tet's four face planes (outward); the
        // last, the front face's, is `below` (reversed: the tet is on its
        // negative side). Only vertices beyond the front face can be inside,
        // and only those whose box meets the tet's.
        let beyond = |x: u32, s: Sign| {
            let [a, b, c] = [0, 1, 2].map(|k| {
                let f = faces[k];
                o4.orient([f[0], f[1], f[2], x])
            });
            [a, b, c, s]
        };
        for (k, &x) in verts.iter().enumerate() {
            if below[k] == Sign::Negative
                && !t.contains(&x)
                && !apart_boxes(&o4.bounds(x), &tet_box)
                && beyond(x, below[k]).iter().all(|&s| s != Sign::Positive)
            {
                return false;
            }
        }
        let sign = |x: u32| match verts.binary_search(&x) {
            Ok(k) => match below[k] {
                Sign::Negative => beyond(x, below[k]),
                // Outside the front face's plane: nothing else matters.
                Sign::Positive => [Sign::Zero, Sign::Zero, Sign::Zero, Sign::Positive],
                Sign::Zero => faces.map(|f| o4.orient([f[0], f[1], f[2], x])),
            },
            Err(_) => faces.map(|f| o4.orient([f[0], f[1], f[2], x])),
        };
        // Points all strictly outside one face plane miss the tet.
        let apart = |s: &[[Sign; 4]]| (0..4).any(|k| s.iter().all(|x| x[k] == Sign::Positive));
        let w = t[3];
        let new_edges = [(w, t[0]), (w, t[1]), (w, t[2])];
        let mut front_edges: BTreeSet<(u32, u32)> = BTreeSet::new();
        for (key, &tau) in front {
            if let Some(k) = faces.iter().position(|&f| sorted3(f) == *key) {
                // Matched from inside: the tet's outward face is the front's.
                if k != 3 && rotated(faces[k]) != rotated(tau) {
                    return false;
                }
                continue;
            }
            // A face whose box misses the tet's meets nothing of it.
            if apart_boxes(&o4.bounds_of(&tau), &tet_box) {
                continue;
            }
            let st = tau.map(sign);
            if apart(&st) {
                continue;
            }
            for (x, y) in new_edges {
                if !(tau.contains(&x) && tau.contains(&y))
                    && self.meets_triangle_by(&|v| o4.orient(v), x, y, tau)
                {
                    return false;
                }
            }
            for k in 0..3 {
                let (x, y) = (tau[k], tau[(k + 1) % 3]);
                if !apart(&[st[k], st[(k + 1) % 3]]) {
                    front_edges.insert((x.min(y), x.max(y)));
                }
            }
        }
        for &(x, y) in &front_edges {
            if t.contains(&x) && t.contains(&y) {
                continue;
            }
            for f in &faces[..3] {
                if self.meets_triangle_by(&|v| o4.orient(v), x, y, *f) {
                    return false;
                }
            }
            for (u, v) in new_edges {
                if ![x, y].contains(&u)
                    && ![x, y].contains(&v)
                    && self.segments_cross(o4, u, v, x, y)
                {
                    return false;
                }
            }
        }
        true
    }

    /// Fill the region bounded by `boundary` (faces outward) with tets on
    /// `verts`, by gift-wrapping. With `open`, faces whose vertices all lie on
    /// the plane `open` accepts are left open, and those of the tets made come
    /// back (outward from them). `None` where a front face has no apex.
    #[allow(clippy::type_complexity)]
    fn wrap_half(
        &self,
        o4: &Orients,
        boundary: &[[u32; 3]],
        verts: &[u32],
        open: Option<&dyn Fn(u32) -> bool>,
    ) -> Option<(Vec<[u32; 4]>, Vec<[u32; 3]>)> {
        let mut front: BTreeMap<[u32; 3], [u32; 3]> = BTreeMap::new();
        for &f in boundary {
            if front.insert(sorted3(f), f).is_some() {
                return None;
            }
        }
        let on_plane = |f: [u32; 3]| open.is_some_and(|on| f.iter().all(|&v| on(v)));
        let mut tets = Vec::new();
        let mut region = Vec::new();
        let budget = WRAP_STEPS * (boundary.len() + verts.len());
        while let Some((&key, &sigma)) = front.iter().next() {
            if tets.len() > budget {
                return None;
            }
            let below: Vec<Sign> = verts
                .iter()
                .map(|&w| o4.orient([sigma[0], sigma[1], sigma[2], w]))
                .collect();
            let mut candidates: Vec<u32> = verts
                .iter()
                .zip(&below)
                .filter(|&(w, &s)| !sigma.contains(w) && s == Sign::Negative)
                .map(|(&w, _)| w)
                .collect();
            // The candidate no other is inside the sphere through the face and
            // it (spheres through the face are ordered along their pencil);
            // if it does not fit, the rest in that order until one does.
            let _s = crate::prof::scope("cdt wrap: apex");
            let nearer = |a: u32, b: u32| {
                self.insphere_perturbed([sigma[0], sigma[2], sigma[1], b, a]) == Sign::Positive
            };
            let apex = |w: u32| [sigma[0], sigma[2], sigma[1], w];
            let &first = candidates.first()?;
            let best = candidates[1..]
                .iter()
                .fold(first, |best, &c| if nearer(c, best) { c } else { best });
            let tet = if self.wrap_fits(o4, apex(best), &front, verts, &below) {
                apex(best)
            } else {
                candidates.retain(|&c| c != best);
                // Binary insertion: a few comparisons each, and no panic
                // whatever the comparisons say.
                let mut order: Vec<u32> = Vec::with_capacity(candidates.len());
                for &c in &candidates {
                    let at = order.partition_point(|&o| nearer(o, c));
                    order.insert(at, c);
                }
                order
                    .iter()
                    .map(|&w| apex(w))
                    .find(|&t| self.wrap_fits(o4, t, &front, verts, &below))?
            };
            drop(_s);
            front.remove(&key);
            for k in 0..3 {
                let f = face(&tet, k);
                let fk = sorted3(f);
                if front.remove(&fk).is_some() {
                    continue;
                }
                if on_plane(f) {
                    region.push(f);
                } else {
                    front.insert(fk, [f[0], f[2], f[1]]);
                }
            }
            tets.push(tet);
        }
        Some((tets, region))
    }

    /// Fill both halves of a facet cavity by gift-wrapping (see
    /// [`ExactDelaunay::wrap_half`]): one leaving the plane open, then the
    /// other closed by the faces it put there; then the other way round.
    #[allow(clippy::type_complexity)]
    pub(super) fn wrap_facet_cavity(
        &self,
        halves: &[Vec<([u32; 3], Option<u32>)>; 2],
        cavity: &[u32],
        side: &dyn Fn(&Self, u32) -> Sign,
    ) -> Option<Vec<[u32; 4]>> {
        let mut verts: Vec<u32> = cavity
            .iter()
            .flat_map(|&t| self.tets[t as usize].v)
            .collect();
        verts.sort_unstable();
        verts.dedup();
        let signs: Vec<Sign> = verts.iter().map(|&v| side(self, v)).collect();
        let upper: Vec<u32> = (0..verts.len())
            .filter(|&k| signs[k] != Sign::Negative)
            .map(|k| verts[k])
            .collect();
        let lower: Vec<u32> = (0..verts.len())
            .filter(|&k| signs[k] != Sign::Positive)
            .map(|k| verts[k])
            .collect();
        let on_plane = |v: u32| {
            verts
                .binary_search(&v)
                .is_ok_and(|k| signs[k] == Sign::Zero)
        };
        let sides = [&upper, &lower];
        let boxes = verts.iter().map(|&v| self.pt(v).bounds()).collect();
        let o4 = Orients {
            dt: self,
            seen: RefCell::new(FxHashMap::default()),
            verts: verts.clone(),
            boxes,
        };
        for (open, closed) in [(0, 1), (1, 0)] {
            let faces: Vec<[u32; 3]> = halves[open].iter().map(|&(f, _)| f).collect();
            let Some((mut tets, region)) =
                self.wrap_half(&o4, &faces, sides[open], Some(&on_plane))
            else {
                continue;
            };
            let mut faces: Vec<[u32; 3]> = halves[closed].iter().map(|&(f, _)| f).collect();
            faces.extend(region.iter().map(|f| [f[0], f[2], f[1]]));
            if let Some((more, _)) = self.wrap_half(&o4, &faces, sides[closed], None) {
                tets.extend(more);
                return Some(tets);
            }
        }
        None
    }

    /// Does every vertex of the tets `old` appear in `new`?
    pub(super) fn keeps_vertices(&self, old: &[u32], new: &[[u32; 4]]) -> bool {
        let kept: BTreeSet<u32> = new.iter().flatten().copied().collect();
        old.iter()
            .flat_map(|&t| self.tets[t as usize].v)
            .all(|v| kept.contains(&v))
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use mesher_geom::Sign;
    use mesher_geom::implicit::Point;

    use super::super::ExactDelaunay;

    /// The perturbed in-sphere test is a consistent tie-break: on five
    /// cospherical points it is never zero, flips with the tet's orientation
    /// and with which point is the query, as the determinant does.
    #[test]
    fn perturbed_insphere_breaks_ties_consistently() {
        // Corners of a cube: every five are cospherical.
        let p: Vec<Point> = (0..8)
            .map(|i| {
                Point::explicit([
                    f64::from(i & 1),
                    f64::from((i >> 1) & 1),
                    f64::from((i >> 2) & 1),
                ])
                .unwrap()
            })
            .collect();
        let d = ExactDelaunay::new(p, 1).unwrap();
        let tet = [0u32, 1, 2, 4];
        let t = if d.orient4(tet) == Sign::Positive {
            tet
        } else {
            [1, 0, 2, 4]
        };
        for e in [3u32, 5, 6, 7] {
            let s = d.insphere_perturbed([t[0], t[1], t[2], t[3], e]);
            assert_ne!(s, Sign::Zero);
            // A negative tet gives the opposite sign.
            let r = d.insphere_perturbed([t[1], t[0], t[2], t[3], e]);
            assert_eq!(r, s.flipped());
            // Swapping the query with a tet vertex flips it.
            let q = d.insphere_perturbed([e, t[1], t[2], t[3], t[0]]);
            assert_eq!(q, s.flipped());
        }
    }
}
