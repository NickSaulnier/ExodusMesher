//! Generalized winding number (Jacobson, Kavan and Sorkine-Hornung, 2013).
//!
//! `w(p) = Σ Ω_t(p) / 4π` over triangles, with signed solid angles. For a closed,
//! outward-oriented surface it is 1 inside and 0 outside; on surfaces with holes
//! or self-intersections it degrades smoothly, which is what makes it a usable
//! inside/outside field for the stuffing kernel and region marking.
//!
//! The fast path is exact in exact arithmetic, not an approximation. For a
//! patch `P` of triangles (a BVH node) with boundary cycle `∂P`, the cone `C`
//! from any point `o` of the patch over `∂P` has the same boundary, so `P − C`
//! is closed and its winding number vanishes outside its convex hull. Hence for
//! `p` outside the node's box, `w_P(p) = w_C(p)`, which costs one solid angle
//! per boundary edge instead of one per triangle. Closed regions have empty
//! boundaries and cost nothing.

use std::f64::consts::PI;

use crate::Point3;
use crate::bvh::Bvh;
use crate::soup::TriangleSoup;
use crate::vec3::{cross, dot, norm, sub};

/// Signed solid angle of triangle `tri` seen from `p` (Van Oosterom and
/// Strackee, 1983), in `[-2π, 2π]`. Positive when `p` is behind the triangle,
/// i.e. on the side opposite its normal `(b-a)×(c-a)`.
#[must_use]
pub fn solid_angle(p: Point3, [a, b, c]: [Point3; 3]) -> f64 {
    let (a, b, c) = (sub(a, p), sub(b, p), sub(c, p));
    let (la, lb, lc) = (norm(a), norm(b), norm(c));
    let numerator = dot(a, cross(b, c));
    let denominator = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
    2.0 * libm::atan2(numerator, denominator)
}

/// Winding number by direct summation over every triangle: the reference.
#[must_use]
pub fn winding_number_brute(soup: &TriangleSoup<'_>, p: Point3) -> f64 {
    soup.triangles()
        .iter()
        .map(|&t| solid_angle(p, soup.corners(t)))
        .sum::<f64>()
        / (4.0 * PI)
}

/// Boundary of a patch: `(i, j, k)` means edge `i→j` with net multiplicity `k`
/// (negative for `j→i`), keyed by `i < j`, sorted.
type Boundary = Vec<(u32, u32, i64)>;

/// Winding number with the hierarchical fast path.
#[derive(Debug, Clone)]
pub struct FastWindingNumber<'a> {
    bvh: Bvh<'a>,
    /// Per node: the boundary, kept only where it is shorter than the node's
    /// triangle list (otherwise the fast path would not be faster).
    caps: Vec<Option<Boundary>>,
}

impl<'a> FastWindingNumber<'a> {
    #[must_use]
    pub fn new(soup: TriangleSoup<'a>) -> FastWindingNumber<'a> {
        let bvh = Bvh::new(soup);
        let mut caps = vec![None; bvh.node_count()];
        if bvh.node_count() > 0 {
            boundary_of(&bvh, 0, &mut caps);
        }
        FastWindingNumber { bvh, caps }
    }

    #[must_use]
    pub fn bvh(&self) -> &Bvh<'a> {
        &self.bvh
    }

    #[must_use]
    pub fn winding_number(&self, p: Point3) -> f64 {
        if self.bvh.node_count() == 0 {
            return 0.0;
        }
        let soup = self.bvh.soup();
        let mut total = 0.0;
        let mut stack = vec![0u32];
        while let Some(n) = stack.pop() {
            if !self.bvh.node_box(n).contains(p) {
                if let Some(cap) = &self.caps[n as usize] {
                    total += cap_solid_angle(soup.points(), cap, p);
                    continue;
                }
            }
            match self.bvh.node_children(n) {
                Ok(ids) => {
                    for &t in ids {
                        total += solid_angle(p, soup.corners(soup.triangles()[t as usize]));
                    }
                }
                Err((left, right)) => stack.extend([right, left]),
            }
        }
        total / (4.0 * PI)
    }
}

/// Sum of solid angles of the cone from the first boundary vertex over the
/// boundary edges. Cone triangles touching the apex are degenerate and
/// contribute zero.
fn cap_solid_angle(points: &[Point3], cap: &Boundary, p: Point3) -> f64 {
    let Some(&(apex, _, _)) = cap.first() else {
        return 0.0;
    };
    let o = points[apex as usize];
    cap.iter()
        .map(|&(i, j, k)| k as f64 * solid_angle(p, [o, points[i as usize], points[j as usize]]))
        .sum()
}

/// Post-order boundary computation; stores useful caps and returns the node's
/// boundary for its parent.
fn boundary_of(bvh: &Bvh<'_>, n: u32, caps: &mut [Option<Boundary>]) -> Boundary {
    let boundary = match bvh.node_children(n) {
        Ok(ids) => {
            let mut edges: Boundary = Vec::with_capacity(ids.len() * 3);
            for &t in ids {
                let tri = bvh.soup().triangles()[t as usize];
                for c in 0..3 {
                    let (i, j) = (tri[c], tri[(c + 1) % 3]);
                    if i < j {
                        edges.push((i, j, 1));
                    } else if j < i {
                        edges.push((j, i, -1));
                    }
                }
            }
            edges.sort_unstable();
            merge(Vec::new(), edges)
        }
        Err((left, right)) => {
            let l = boundary_of(bvh, left, caps);
            let r = boundary_of(bvh, right, caps);
            merge(l, r)
        }
    };
    if boundary.len() < bvh.node_triangle_count(n) {
        caps[n as usize] = Some(boundary.clone());
    }
    boundary
}

/// Merge two sorted boundaries, summing multiplicities and dropping zeros.
fn merge(a: Boundary, b: Boundary) -> Boundary {
    let mut out: Boundary = Vec::with_capacity(a.len() + b.len());
    let mut all = a.into_iter().chain(b).collect::<Vec<_>>();
    all.sort_unstable_by_key(|&(i, j, _)| (i, j));
    for (i, j, k) in all {
        match out.last_mut() {
            Some(last) if (last.0, last.1) == (i, j) => last.2 += k,
            _ => out.push((i, j, k)),
        }
    }
    out.retain(|&(_, _, k)| k != 0);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tetra() -> (Vec<Point3>, Vec<[u32; 3]>) {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        // Outward-facing.
        let tris = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
        (pts, tris)
    }

    #[test]
    fn closed_tetra_inside_outside() {
        let (pts, tris) = tetra();
        let soup = TriangleSoup::new(&pts, &tris).unwrap();
        let inside = [0.1, 0.1, 0.1];
        let outside = [2.0, 2.0, 2.0];
        assert!((winding_number_brute(&soup, inside) - 1.0).abs() < 1e-12);
        assert!(winding_number_brute(&soup, outside).abs() < 1e-12);
        let fast = FastWindingNumber::new(soup);
        assert!((fast.winding_number(inside) - 1.0).abs() < 1e-12);
        assert!(fast.winding_number(outside).abs() < 1e-12);
    }

    #[test]
    fn single_triangle_hemisphere() {
        // A huge triangle seen from just behind its centre covers half the sphere.
        let tri = [[-1e6, -1e6, 0.0], [1e6, -1e6, 0.0], [0.0, 1e6, 0.0]];
        let w = solid_angle([0.0, 0.0, -1.0], tri) / (4.0 * PI);
        assert!((w - 0.5).abs() < 1e-5);
    }

    #[test]
    fn closed_region_has_empty_root_boundary() {
        let (pts, tris) = tetra();
        let fast = FastWindingNumber::new(TriangleSoup::new(&pts, &tris).unwrap());
        assert_eq!(fast.caps[0].as_deref(), Some(&[][..]));
    }

    #[test]
    fn fast_path_is_used_on_large_patches() {
        // A flat 16×16 grid of quads: BVH nodes are compact patches whose
        // boundaries are much shorter than their triangle lists.
        let n = 16u32;
        let pts: Vec<Point3> = (0..=n)
            .flat_map(|i| (0..=n).map(move |j| [f64::from(i), f64::from(j), 0.0]))
            .collect();
        let id = |i: u32, j: u32| i * (n + 1) + j;
        let tris: Vec<[u32; 3]> = (0..n)
            .flat_map(|i| (0..n).map(move |j| (i, j)))
            .flat_map(|(i, j)| {
                [
                    [id(i, j), id(i + 1, j), id(i + 1, j + 1)],
                    [id(i, j), id(i + 1, j + 1), id(i, j + 1)],
                ]
            })
            .collect();
        let soup = TriangleSoup::new(&pts, &tris).unwrap();
        let fast = FastWindingNumber::new(soup);
        let cached = fast.caps.iter().filter(|c| c.is_some()).count();
        assert!(cached >= 3, "only {cached} nodes cache a boundary");
        let p = [3.3, 7.1, 0.25];
        let (f, b) = (fast.winding_number(p), winding_number_brute(&soup, p));
        assert!((f - b).abs() < 1e-12, "{f} vs {b}");
    }

    #[test]
    fn empty_soup_is_zero() {
        let fast = FastWindingNumber::new(TriangleSoup::new(&[], &[]).unwrap());
        assert_eq!(fast.winding_number([0.0; 3]), 0.0);
    }
}
