//! Where extra points lie on the input surface, for the enriched CDT
//! (Milestone 14): a point exactly on an input segment becomes a vertex of
//! its chain, and a point exactly inside an input triangle a vertex inside
//! its facet, so the constraints still partition the surface.
//!
//! Every decision is exact: `orient3d` for being on a triangle's plane, then
//! the projected orientation against each side. f64 boxes and plane
//! distances (to 1e-9 of the bounding-box diagonal, far above the points'
//! rounding) only skip triangles a point is clearly off.

use std::collections::{BTreeMap, BTreeSet};

use mesher_geom::implicit::{self, Point};
use mesher_geom::{Point3, Sign};

/// The extra points found on the input surface.
#[derive(Debug, Clone, Default)]
pub(crate) struct OnSurface {
    /// Each segment's (lower vertex first) extra points strictly inside it,
    /// in order from its lower vertex.
    pub on_segment: BTreeMap<(u32, u32), Vec<u32>>,
    /// Each input triangle's extra points strictly inside it.
    pub inside: Vec<Vec<u32>>,
    /// Extra points equal to an input vertex.
    pub coincident: Vec<u32>,
}

fn sub(a: Point3, b: Point3) -> Point3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: Point3, b: Point3) -> Point3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: Point3, b: Point3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Where `p` is on triangle `abc`: `None` off it; `Some(None)` strictly
/// inside; `Some(Some(k))` strictly inside side `k` (0: a b, 1: b c, 2:
/// c a); a corner is reported as `Some(Some(3 + corner))`.
fn on_triangle(
    a: &Point,
    b: &Point,
    c: &Point,
    p: &Point,
    normal: Point3,
) -> Option<Option<usize>> {
    if implicit::orient3d(a, b, c, p) != Sign::Zero {
        return None;
    }
    // Project along the normal's largest component (any axis the triangle
    // does not project flat along works).
    let mut axes = [0usize, 1, 2];
    axes.sort_by(|&x, &y| normal[y].abs().total_cmp(&normal[x].abs()));
    let axis = axes
        .into_iter()
        .find(|&k| implicit::orient2d_projected(a, b, c, k) != Sign::Zero)?;
    let up = implicit::orient2d_projected(a, b, c, axis);
    let side = |x: &Point, y: &Point| {
        let s = implicit::orient2d_projected(x, y, p, axis);
        if s == Sign::Zero {
            0
        } else if s == up {
            1
        } else {
            -1
        }
    };
    let s = [side(a, b), side(b, c), side(c, a)];
    if s.contains(&-1) {
        return None;
    }
    let zeros: Vec<usize> = (0..3).filter(|&k| s[k] == 0).collect();
    match zeros.as_slice() {
        [] => Some(None),
        [k] => Some(Some(*k)),
        // On two sides: at the corner they share (sides k and k + 1 meet at
        // corner k + 1).
        [0, 1] => Some(Some(4)),
        [1, 2] => Some(Some(5)),
        [0, 2] => Some(Some(3)),
        _ => None,
    }
}

/// Locate the points `all[first_extra..]` on the input triangles (whose
/// corners are `all[..first_extra]`, at `points`).
pub(crate) fn locate(
    points: &[Point3],
    triangles: &[[u32; 3]],
    all: &[Point],
    first_extra: u32,
) -> OnSurface {
    let _s = crate::prof::scope("cdt locate");
    let mut out = OnSurface {
        inside: vec![Vec::new(); triangles.len()],
        ..OnSurface::default()
    };
    if all.len() <= first_extra as usize || triangles.is_empty() {
        return out;
    }
    let (lo, hi) = points.iter().fold(
        ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
        |(lo, hi), p| {
            (
                [0, 1, 2].map(|k| lo[k].min(p[k])),
                [0, 1, 2].map(|k| hi[k].max(p[k])),
            )
        },
    );
    let diagonal = dot(sub(hi, lo), sub(hi, lo)).sqrt();
    let tol = 1e-9 * diagonal;
    // A uniform grid of about two triangles a cell, by their boxes.
    let mut n = 1usize;
    while n * n * n < 2 * triangles.len() && n < 64 {
        n += 1;
    }
    let size = [0, 1, 2].map(|k| ((hi[k] - lo[k]) / n as f64).max(f64::MIN_POSITIVE));
    let cell = |x: f64, k: usize| (((x - lo[k]) / size[k]).floor().max(0.0) as usize).min(n - 1);
    let mut grid: Vec<Vec<u32>> = vec![Vec::new(); n * n * n];
    let mut boxes: Vec<(Point3, Point3, Point3)> = Vec::with_capacity(triangles.len());
    for (ti, t) in triangles.iter().enumerate() {
        let [a, b, c] = t.map(|v| points[v as usize]);
        let tlo = [0, 1, 2].map(|k| a[k].min(b[k]).min(c[k]) - tol);
        let thi = [0, 1, 2].map(|k| a[k].max(b[k]).max(c[k]) + tol);
        boxes.push((tlo, thi, cross(sub(b, a), sub(c, a))));
        let (c0, c1) = (
            [0, 1, 2].map(|k| cell(tlo[k], k)),
            [0, 1, 2].map(|k| cell(thi[k], k)),
        );
        for i in c0[0]..=c1[0] {
            for j in c0[1]..=c1[1] {
                for k in c0[2]..=c1[2] {
                    grid[(i * n + j) * n + k].push(ti as u32);
                }
            }
        }
    }
    let mut on_segment: BTreeMap<(u32, u32), BTreeSet<u32>> = BTreeMap::new();
    let mut coincident: BTreeSet<u32> = BTreeSet::new();
    for x in first_extra..all.len() as u32 {
        let p = &all[x as usize];
        let q = p.approx();
        if !q.iter().all(|c| c.is_finite()) {
            continue;
        }
        let at = (0..3).map(|k| cell(q[k], k)).collect::<Vec<_>>();
        for &ti in &grid[(at[0] * n + at[1]) * n + at[2]] {
            let (tlo, thi, normal) = boxes[ti as usize];
            if (0..3).any(|k| q[k] < tlo[k] || q[k] > thi[k]) {
                continue;
            }
            let t = triangles[ti as usize];
            let a = points[t[0] as usize];
            let length = dot(normal, normal).sqrt();
            if dot(normal, sub(q, a)).abs() > tol * length {
                continue;
            }
            let [pa, pb, pc] = t.map(|v| &all[v as usize]);
            match on_triangle(pa, pb, pc, p, normal) {
                None => {}
                Some(None) => out.inside[ti as usize].push(x),
                Some(Some(k)) if k < 3 => {
                    let (u, v) = (t[k], t[(k + 1) % 3]);
                    on_segment
                        .entry((u.min(v), u.max(v)))
                        .or_default()
                        .insert(x);
                }
                Some(Some(_)) => {
                    coincident.insert(x);
                }
            }
        }
    }
    for list in &mut out.inside {
        list.sort_unstable();
        list.dedup();
        list.retain(|x| !coincident.contains(x));
    }
    for ((u, v), set) in on_segment {
        let mut list: Vec<u32> = set
            .into_iter()
            .filter(|x| !coincident.contains(x))
            .collect();
        let start = &all[u as usize];
        // p comes before q when it lies between the segment's start and q.
        list.sort_by(|&p, &q| {
            if p == q {
                std::cmp::Ordering::Equal
            } else if implicit::dot_sign(start, &all[p as usize], &all[q as usize])
                == Sign::Negative
            {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            }
        });
        if !list.is_empty() {
            out.on_segment.insert((u, v), list);
        }
    }
    // A point on a segment is not inside the triangles beside it.
    let on_any: BTreeSet<u32> = out.on_segment.values().flatten().copied().collect();
    for list in &mut out.inside {
        list.retain(|x| !on_any.contains(x));
    }
    out.coincident = coincident.into_iter().collect();
    out
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;

    #[test]
    fn points_on_a_tetrahedron_are_found_on_its_edges_and_faces() {
        let points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let triangles = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
        let mut all: Vec<Point> = points.iter().map(|&p| Point::from(p)).collect();
        // On edge 0-1 (twice, out of order), inside face 0 2 1, on the slanted
        // face 1 2 3 as an LNC of its corners, off the surface, at a corner.
        all.push(Point::from([0.75, 0.0, 0.0]));
        all.push(Point::from([0.25, 0.0, 0.0]));
        all.push(Point::from([0.25, 0.25, 0.0]));
        all.push(Point::lnc([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.3).unwrap());
        all.push(Point::bpt(points_of(&points, [1, 2, 3]), [1.0, 2.0, 3.0]).unwrap());
        all.push(Point::from([0.2, 0.2, 0.2]));
        all.push(Point::from([0.0, 0.0, 1.0]));
        let found = locate(&points, &triangles, &all, 4);
        assert_eq!(found.on_segment.get(&(0, 1)), Some(&vec![5, 4]));
        assert_eq!(found.on_segment.get(&(1, 2)), Some(&vec![7]));
        assert_eq!(found.inside[0], vec![6]);
        assert_eq!(found.inside[2], vec![8]);
        assert!(found.inside[1].is_empty() && found.inside[3].is_empty());
        assert_eq!(found.coincident, vec![10]);
    }

    fn points_of(points: &[Point3], t: [usize; 3]) -> [Point3; 3] {
        t.map(|v| points[v])
    }
}
