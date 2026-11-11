//! The enriched CDT (Milestone 14, the paper's Section 6): the refined
//! mesh's points with the unchamfered input, through the exact constrained
//! kernel, so the mesh conforms to the input exactly, with the refined
//! quality everywhere but near the former acute angles and filtered
//! features.
//!
//! The paper's input is the refined points outside the chamfering balls
//! and strips, the input's vertices, the refined faces partitioning the
//! chamfered surface, and the cut-away pieces of the input. Here the
//! constraints are the input triangles themselves: every refined point on
//! the chamfered surface lies exactly on an input edge or triangle (cuts and
//! split points on edges are LNC points of the input edge, bridge and face
//! points BPT points of an input triangle, snapped points on axis-aligned
//! planes exactly on them), and the kernel places each one there exactly
//! (`cdt::recover_faces_with`). The constraints still partition the input
//! surface, and every point the kernel adds is built from the input's own
//! explicit vertices, so no implicit point is built on another (the paper's
//! Appendix C, cascaded implicit points, is not needed).
//!
//! Refined points off the surface are left out where the refinement could
//! not see the input: inside a chamfered vertex's ball, an acute edge's
//! strip, or within τ of a filtered triangle. So are the split points inside
//! bridge segments: they lie exactly on their input triangle but only near
//! the bridge (BPT points at averaged weights), so a chain of them zig-zags
//! at the rounding level, and with the bridge no constraint here the
//! kernel's facet can take a face across three of them, at a face angle of
//! 0° (128915, 74494 and 73154 had 61 to 136 tets under 1° with them, none
//! without). Coplanar input triangles across flat edges are one
//! facet, as refinement's planar regions are, so flat edges need no points.
//! Which tets are inside the input is the kernel's parity flood
//! (`cdt::mesh_region_with`).

use std::collections::BTreeSet;

use mesher_geom::Point3;
use mesher_geom::implicit::Point;

use super::Refined;
use crate::cdt::{Region, RegionError, mesh_region_with};
use crate::chamfer::{ChamferVertex, FilteredChamfer};

/// The enriched CDT's region and what went into it.
#[derive(Debug, Clone)]
pub struct Enriched {
    /// The tets inside the input, with exact points: the input's vertices
    /// first, then the refined points kept, then the kernel's Steiner
    /// points.
    pub region: Region,
    /// Refined points kept (copies of input vertices aside).
    pub kept: usize,
    /// Refined points left out: in a chamfered vertex's ball, in an acute
    /// edge's strip, near a filtered triangle, and split points on bridges.
    pub cleared: [usize; 4],
}

fn sub(a: Point3, b: Point3) -> Point3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: Point3, b: Point3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn add_scaled(a: Point3, d: Point3, t: f64) -> Point3 {
    [a[0] + t * d[0], a[1] + t * d[1], a[2] + t * d[2]]
}

fn distance2(a: Point3, b: Point3) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}

/// Squared distance from `p` to segment `a b`.
fn segment_distance2(p: Point3, a: Point3, b: Point3) -> f64 {
    let d = sub(b, a);
    let l2 = dot(d, d);
    let t = if l2 > 0.0 {
        (dot(sub(p, a), d) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    distance2(p, add_scaled(a, d, t))
}

/// Squared distance from `p` to triangle `a b c` (the nearest point by
/// regions, after Ericson, Real-Time Collision Detection, 5.1.5).
fn triangle_distance2(p: Point3, a: Point3, b: Point3, c: Point3) -> f64 {
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return distance2(p, a);
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return distance2(p, b);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return distance2(p, add_scaled(a, ab, d1 / (d1 - d3)));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return distance2(p, c);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return distance2(p, add_scaled(a, ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return distance2(p, add_scaled(b, sub(c, b), w));
    }
    let denom = 1.0 / (va + vb + vc);
    let (v, w) = (vb * denom, vc * denom);
    distance2(p, add_scaled(add_scaled(a, ab, v), ac, w))
}

/// Boxes in a uniform grid, for the boxes a point may be in.
struct BoxGrid {
    lo: Point3,
    size: Point3,
    n: usize,
    cells: Vec<Vec<u32>>,
}

impl BoxGrid {
    fn new(boxes: &[(Point3, Point3)]) -> BoxGrid {
        let (lo, hi) = boxes.iter().fold(
            ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
            |(lo, hi), (a, b)| {
                (
                    [0, 1, 2].map(|k| lo[k].min(a[k])),
                    [0, 1, 2].map(|k| hi[k].max(b[k])),
                )
            },
        );
        let mut n = 1usize;
        while n * n * n < 2 * boxes.len() && n < 64 {
            n += 1;
        }
        let size = [0, 1, 2].map(|k| ((hi[k] - lo[k]) / n as f64).max(f64::MIN_POSITIVE));
        let mut grid = BoxGrid {
            lo,
            size,
            n,
            cells: vec![Vec::new(); n * n * n],
        };
        for (i, (a, b)) in boxes.iter().enumerate() {
            let (c0, c1) = (grid.cell(*a), grid.cell(*b));
            for x in c0[0]..=c1[0] {
                for y in c0[1]..=c1[1] {
                    for z in c0[2]..=c1[2] {
                        grid.cells[(x * n + y) * n + z].push(i as u32);
                    }
                }
            }
        }
        grid
    }

    fn cell(&self, p: Point3) -> [usize; 3] {
        [0, 1, 2].map(|k| {
            (((p[k] - self.lo[k]) / self.size[k]).floor().max(0.0) as usize).min(self.n - 1)
        })
    }

    /// The boxes in `p`'s cell (a superset of those holding `p`, if `p` is
    /// in the grid's bounds).
    fn near(&self, p: Point3) -> &[u32] {
        if self.cells.is_empty() || (0..3).any(|k| !p[k].is_finite()) {
            return &[];
        }
        let n = self.n;
        let c = self.cell(p);
        &self.cells[(c[0] * n + c[1]) * n + c[2]]
    }
}

/// A region a refined point off the surface is left out of.
enum Zone {
    Ball(Point3, f64),
    Strip(Point3, Point3, f64),
    Near([Point3; 3], f64),
}

impl Zone {
    fn bounds(&self) -> (Point3, Point3) {
        let grow = |pts: &[Point3], r: f64| {
            pts.iter().fold(
                ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
                |(lo, hi), p| {
                    (
                        [0, 1, 2].map(|k| lo[k].min(p[k] - r)),
                        [0, 1, 2].map(|k| hi[k].max(p[k] + r)),
                    )
                },
            )
        };
        match self {
            Zone::Ball(c, r) => grow(&[*c], *r),
            Zone::Strip(a, b, r) => grow(&[*a, *b], *r),
            Zone::Near(t, r) => grow(t, *r),
        }
    }

    /// Strictly inside? (f64: it only chooses points to leave out.)
    fn holds(&self, p: Point3) -> bool {
        match self {
            Zone::Ball(c, r) => distance2(p, *c) < r * r,
            Zone::Strip(a, b, r) => segment_distance2(p, *a, *b) < r * r,
            Zone::Near(t, r) => triangle_distance2(p, t[0], t[1], t[2]) < r * r,
        }
    }

    fn kind(&self) -> usize {
        match self {
            Zone::Ball(..) => 0,
            Zone::Strip(..) => 1,
            Zone::Near(..) => 2,
        }
    }
}

/// Make the enriched CDT of `refined` (from [`super::refine_conforming`]
/// or the face stage), `surface` (its filtered chamfer) and the input
/// `points` and `triangles` it was chamfered from. `tau` is the filtering
/// tolerance (absolute: τ of the diagonal); `max_points` and `seed` go to
/// the kernel.
///
/// # Errors
/// See [`RegionError`].
pub fn enriched_cdt(
    refined: &Refined,
    surface: &FilteredChamfer,
    points: &[Point3],
    triangles: &[[u32; 3]],
    tau: f64,
    max_points: usize,
    seed: u64,
) -> Result<Enriched, RegionError> {
    let chamfered = &surface.chamfered;
    let mut zones: Vec<Zone> = Vec::new();
    for (v, &r) in chamfered.radius.iter().enumerate() {
        if r > 0.0 {
            zones.push(Zone::Ball(points[v], r));
        }
    }
    for (&(a, b), &w) in &chamfered.strips {
        if w > 0.0 {
            zones.push(Zone::Strip(points[a as usize], points[b as usize], w));
        }
    }
    if tau > 0.0 {
        for (t, tri) in triangles.iter().enumerate() {
            if !surface.kept.get(t).copied().unwrap_or(true) {
                zones.push(Zone::Near(tri.map(|v| points[v as usize]), tau));
            }
        }
    }
    let grid = BoxGrid::new(&zones.iter().map(Zone::bounds).collect::<Vec<_>>());
    // Points on the chamfered surface stay; copies of input vertices are the
    // input's own.
    let mut on_surface: BTreeSet<u32> = BTreeSet::new();
    for s in &refined.subsegments {
        on_surface.extend(s);
    }
    for (f, _) in &refined.subfaces {
        on_surface.extend(f);
    }
    let mut copies: BTreeSet<u32> = BTreeSet::new();
    for (c, v) in refined.vertex_of.iter().enumerate() {
        if let Some(v) = v {
            on_surface.insert(*v);
            if matches!(chamfered.vertices.get(c), Some(ChamferVertex::Input(_))) {
                copies.insert(*v);
            }
        }
    }
    let bridges: BTreeSet<u32> = refined.on_bridges.iter().copied().collect();
    let mut extra: Vec<Point> = Vec::new();
    let mut cleared = [0usize; 4];
    for (v, p) in refined.dt.points().iter().enumerate() {
        let v = v as u32;
        if copies.contains(&v) {
            continue;
        }
        if bridges.contains(&v) {
            cleared[3] += 1;
            continue;
        }
        if !on_surface.contains(&v) {
            let q = p.approx();
            if let Some(z) = grid
                .near(q)
                .iter()
                .map(|&z| &zones[z as usize])
                .find(|z| z.holds(q))
            {
                cleared[z.kind()] += 1;
                continue;
            }
        }
        extra.push(*p);
    }
    let kept = extra.len();
    let region = mesh_region_with(points, &extra, triangles, true, max_points, seed)?;
    Ok(Enriched {
        region,
        kept,
        cleared,
    })
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;

    #[test]
    fn distances_to_a_triangle_by_region() {
        let (a, b, c) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let d = |p: Point3| triangle_distance2(p, a, b, c).sqrt();
        assert!((d([0.25, 0.25, 2.0]) - 2.0).abs() < 1e-12);
        assert!((d([-1.0, -1.0, 0.0]) - 2f64.sqrt()).abs() < 1e-12);
        assert!((d([0.5, -1.0, 0.0]) - 1.0).abs() < 1e-12);
        assert!((d([1.0, 1.0, 0.0]) - 0.5f64.sqrt()).abs() < 1e-12);
        assert!((d([2.0, 0.0, 0.0]) - 1.0).abs() < 1e-12);
        assert!((segment_distance2([0.5, 1.0, 0.0], a, b).sqrt() - 1.0).abs() < 1e-12);
    }
}
