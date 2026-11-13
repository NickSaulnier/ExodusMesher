//! Intersections of triangle soups as implicit points (Milestone 16, stage 1).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeMap;

use mesher_geom::Point3;
use mesher_geom::arrangement::{ArrangementError, Intersections, intersections};
use mesher_geom::implicit::{self, Point};

/// An axis-aligned box as 12 outward triangles.
fn cuboid(lo: Point3, hi: Point3, offset: u32) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let p: Vec<Point3> = (0..8)
        .map(|i| {
            [
                if i & 1 == 0 { lo[0] } else { hi[0] },
                if i & 2 == 0 { lo[1] } else { hi[1] },
                if i & 4 == 0 { lo[2] } else { hi[2] },
            ]
        })
        .collect();
    let f = [
        [0, 2, 3],
        [0, 3, 1],
        [4, 5, 7],
        [4, 7, 6],
        [0, 1, 5],
        [0, 5, 4],
        [2, 6, 7],
        [2, 7, 3],
        [0, 4, 6],
        [0, 6, 2],
        [1, 3, 7],
        [1, 7, 5],
    ];
    (p, f.iter().map(|t| t.map(|v| v + offset)).collect())
}

fn union(parts: &[(Vec<Point3>, Vec<[u32; 3]>)]) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let mut p = Vec::new();
    let mut f = Vec::new();
    for (pp, ff) in parts {
        let o = p.len() as u32;
        p.extend(pp.iter().copied());
        f.extend(ff.iter().map(|t| t.map(|v| v + o)));
    }
    (p, f)
}

fn approx_eq(a: Point3, b: Point3) -> bool {
    (0..3).all(|k| (a[k] - b[k]).abs() < 1e-12)
}

/// Every segment of every triangle, once, with how many triangles hold it.
fn segments(x: &Intersections) -> BTreeMap<[u32; 2], usize> {
    let mut out = BTreeMap::new();
    for c in &x.constraints {
        for &s in &c.segments {
            *out.entry(s).or_insert(0) += 1;
        }
    }
    out
}

/// The intersection curves' points by how many distinct segments meet there.
fn degrees(x: &Intersections) -> BTreeMap<u32, usize> {
    let mut d = BTreeMap::new();
    for s in segments(x).keys() {
        for &v in s {
            *d.entry(v).or_insert(0) += 1;
        }
    }
    d
}

/// No two points coincide.
fn distinct(x: &Intersections) {
    for i in 0..x.points.len() {
        for j in i + 1..x.points.len() {
            assert!(
                !implicit::coincident(&x.points[i], &x.points[j]),
                "points {i} and {j} coincide"
            );
        }
    }
}

#[test]
fn two_crossing_triangles_meet_in_a_segment() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [1.0, 1.0, -1.0],
        [1.0, 1.0, 1.0],
        [3.0, 1.0, 1.0],
    ];
    let f = vec![[0, 1, 2], [3, 4, 5]];
    let x = intersections(&p, &f).unwrap();
    assert_eq!(x.pairs, 1);
    assert_eq!(x.points.len(), 8);
    // The second triangle's edges 3-4 and 5-3 pierce the plane z = 0 at
    // (1, 1, 0) and (2, 1, 0), both inside the first.
    let new: Vec<Point3> = x.points[6..].iter().map(Point::approx).collect();
    assert!(
        new.iter().any(|&q| approx_eq(q, [1.0, 1.0, 0.0])),
        "{new:?}"
    );
    assert!(
        new.iter().any(|&q| approx_eq(q, [2.0, 1.0, 0.0])),
        "{new:?}"
    );
    for c in &x.constraints {
        assert_eq!(c.segments, vec![[6, 7]]);
    }
    // Inside the first triangle; on the second's edges.
    assert_eq!(x.constraints[0].points, vec![6, 7]);
    assert_eq!(x.constraints[1].points, vec![6, 7]);
}

#[test]
fn a_triangle_resting_on_another_splits_it_at_the_contact() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 1.0, 1.0],
        [2.0, 1.0, 1.0],
    ];
    let f = vec![[0, 1, 2], [3, 4, 5]];
    let x = intersections(&p, &f).unwrap();
    assert_eq!(x.points.len(), 6);
    assert_eq!(x.constraints[0].points, vec![3]);
    assert!(x.constraints[0].segments.is_empty());
    assert!(x.constraints[1].is_empty());
}

#[test]
fn coplanar_triangles_take_each_others_edges() {
    // Two triangles in z = 0 overlapping in a quadrilateral.
    let p = vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [1.0, 1.0, 0.0],
        [5.0, 1.0, 0.0],
        [1.0, 5.0, 0.0],
    ];
    let f = vec![[0, 1, 2], [3, 4, 5]];
    let x = intersections(&p, &f).unwrap();
    assert_eq!(x.coplanar_pairs, 1);
    distinct(&x);
    // The second's corner 3 is inside the first; the edges 3-4 and 3-5
    // cross the first's hypotenuse at (3, 1) and (1, 3).
    let new: Vec<Point3> = x.points[6..].iter().map(Point::approx).collect();
    assert_eq!(new.len(), 2, "{new:?}");
    assert!(new.iter().any(|&q| approx_eq(q, [3.0, 1.0, 0.0])));
    assert!(new.iter().any(|&q| approx_eq(q, [1.0, 3.0, 0.0])));
    // The first is cut along the second's edges 3-4 and 3-5 up to its
    // hypotenuse; the second along the first's hypotenuse, the only edge of
    // the first that enters it.
    assert_eq!(x.constraints[0].segments.len(), 2);
    assert_eq!(x.constraints[0].points, vec![3, 6, 7]);
    assert_eq!(x.constraints[1].segments, vec![[6, 7]]);
    assert_eq!(x.constraints[1].points, vec![6, 7]);
}

#[test]
fn overlapping_boxes_meet_in_one_closed_curve() {
    // The suite's overlap: the unit cube and a shifted copy. The shifted
    // box's edges meet the cube's faces on their diagonals, so many points
    // land exactly on input edges.
    let (p, f) = union(&[
        cuboid([0.0; 3], [1.0; 3], 0),
        cuboid([0.5, 0.25, 0.25], [1.5, 1.25, 1.25], 0),
    ]);
    let x = intersections(&p, &f).unwrap();
    distinct(&x);
    assert!(x.pairs > 0);
    let d = degrees(&x);
    assert!(!d.is_empty());
    for (&v, &n) in &d {
        assert_eq!(n, 2, "point {v} at {:?}", x.points[v as usize].approx());
    }
    // One loop.
    let segs: Vec<[u32; 2]> = segments(&x).into_keys().collect();
    let mut seen = vec![false; segs.len()];
    let mut loops = 0;
    for start in 0..segs.len() {
        if seen[start] {
            continue;
        }
        loops += 1;
        let mut stack = vec![start];
        while let Some(s) = stack.pop() {
            if seen[s] {
                continue;
            }
            seen[s] = true;
            for (t, other) in segs.iter().enumerate() {
                if !seen[t] && other.iter().any(|v| segs[s].contains(v)) {
                    stack.push(t);
                }
            }
        }
    }
    assert_eq!(loops, 1);
    // Every curve point lies on the boundary of both boxes:
    // on a face of each (x, y or z at one of its bounds).
    for &v in d.keys() {
        let q = x.points[v as usize].approx();
        let on = |lo: Point3, hi: Point3| {
            (0..3).any(|k| (q[k] - lo[k]).abs() < 1e-12 || (q[k] - hi[k]).abs() < 1e-12)
        };
        assert!(on([0.0; 3], [1.0; 3]) && on([0.5, 0.25, 0.25], [1.5, 1.25, 1.25]));
    }
}

#[test]
fn points_on_an_edge_reach_both_of_its_triangles() {
    let (p, f) = union(&[
        cuboid([0.0; 3], [1.0; 3], 0),
        cuboid([0.3, 0.2, 0.1], [1.4, 0.7, 0.6], 0),
    ]);
    let x = intersections(&p, &f).unwrap();
    // Each triangle's point list holds every curve point on its boundary
    // that the neighbour across the edge has on that edge: check by the
    // shared edges of the first box.
    let mut by_edge: BTreeMap<(u32, u32), Vec<usize>> = BTreeMap::new();
    for (t, tri) in f.iter().enumerate() {
        for e in 0..3 {
            let (u, v) = (tri[e], tri[(e + 1) % 3]);
            by_edge.entry((u.min(v), u.max(v))).or_default().push(t);
        }
    }
    for ((u, v), ts) in by_edge {
        assert_eq!(ts.len(), 2);
        let on_edge = |t: usize| -> Vec<u32> {
            x.constraints[t]
                .points
                .iter()
                .copied()
                .filter(|&i| {
                    implicit::collinear(
                        &Point::from(p[u as usize]),
                        &Point::from(p[v as usize]),
                        &x.points[i as usize],
                    )
                })
                .collect()
        };
        assert_eq!(on_edge(ts[0]), on_edge(ts[1]), "edge {u}-{v}");
    }
}

#[test]
fn a_degenerate_triangle_is_refused() {
    let p = vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    assert!(matches!(
        intersections(&p, &[[0, 1, 2]]),
        Err(ArrangementError::Degenerate { triangle: 0 })
    ));
}
