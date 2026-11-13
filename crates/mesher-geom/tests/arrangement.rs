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
        for s in &c.segments {
            *out.entry(s.ends).or_insert(0) += 1;
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
        let ends: Vec<[u32; 2]> = c.segments.iter().map(|s| s.ends).collect();
        assert_eq!(ends, vec![[6, 7]]);
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
    assert_eq!(x.coplanar, vec![(0, 1)]);
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
    assert_eq!(x.constraints[1].segments.len(), 1);
    assert_eq!(x.constraints[1].segments[0].ends, [6, 7]);
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

// Stage 2: the arrangement.

use mesher_geom::arrangement::{Arrangement, arrangement};

fn area(p: Point3, q: Point3, r: Point3) -> f64 {
    let u = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
    let v = [r[0] - p[0], r[1] - p[1], r[2] - p[2]];
    let c = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
}

fn normal(a: Point3, b: Point3, c: Point3) -> Point3 {
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let e = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    [
        d[1] * e[2] - d[2] * e[1],
        d[2] * e[0] - d[0] * e[2],
        d[0] * e[1] - d[1] * e[0],
    ]
}

/// Each input triangle's pieces cover its area and face its way; for a
/// closed input, every directed edge has its reverse as often.
fn check(p: &[Point3], f: &[[u32; 3]], a: &Arrangement, closed: bool) {
    let pos = |i: u32| a.points[i as usize].approx();
    let mut covered = vec![0.0; f.len()];
    for (t, &o) in a.triangles.iter().zip(&a.origin) {
        let [x, y, z] = t.map(pos);
        covered[o as usize] += area(x, y, z);
        let [u, v, w] = f[o as usize].map(|i| p[i as usize]);
        let (n1, n2) = (normal(x, y, z), normal(u, v, w));
        assert!(
            n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2] > 0.0,
            "piece {t:?} of {o} is turned over"
        );
    }
    for (o, t) in f.iter().enumerate() {
        let [u, v, w] = t.map(|i| p[i as usize]);
        let want = area(u, v, w);
        assert!(
            (covered[o] - want).abs() <= 1e-9 * want.max(1.0),
            "triangle {o}: pieces cover {} of {want}",
            covered[o]
        );
    }
    if closed {
        let mut directed: BTreeMap<(u32, u32), i64> = BTreeMap::new();
        for t in &a.triangles {
            for k in 0..3 {
                let (u, v) = (t[k], t[(k + 1) % 3]);
                *directed.entry((u.min(v), u.max(v))).or_insert(0) += if u < v { 1 } else { -1 };
            }
        }
        for (e, n) in directed {
            assert_eq!(n, 0, "edge {e:?} is not balanced");
        }
    }
}

#[test]
fn crossing_triangles_are_cut_along_their_segment() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [1.0, 1.0, -1.0],
        [1.0, 1.0, 1.0],
        [3.0, 1.0, 1.0],
    ];
    let f = vec![[0, 1, 2], [3, 4, 5]];
    let a = arrangement(&p, &f).unwrap();
    check(&p, &f, &a, false);
    // The segment is an edge of pieces on both triangles.
    for o in 0..2 {
        assert!(
            a.triangles
                .iter()
                .zip(&a.origin)
                .any(|(t, &g)| g == o && t.contains(&6) && t.contains(&7))
        );
    }
}

#[test]
fn overlapping_boxes_make_a_closed_arrangement() {
    for second in [
        ([0.5, 0.25, 0.25], [1.5, 1.25, 1.25]),
        ([0.3, 0.2, 0.1], [1.4, 0.7, 0.6]),
        // Sharing the planes z = 0 and x = 0: coplanar overlaps.
        ([0.0, 0.5, 0.0], [0.5, 1.5, 1.0]),
    ] {
        let (p, f) = union(&[cuboid([0.0; 3], [1.0; 3], 0), cuboid(second.0, second.1, 0)]);
        let a = arrangement(&p, &f).unwrap();
        check(&p, &f, &a, true);
    }
}

#[test]
fn coplanar_overlap_gives_the_same_pieces_to_both() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [1.0, 1.0, 0.0],
        [5.0, 1.0, 0.0],
        [1.0, 5.0, 0.0],
    ];
    // The second faces the other way.
    let f = vec![[0, 1, 2], [3, 5, 4]];
    let a = arrangement(&p, &f).unwrap();
    check(&p, &f, &a, false);
    let key = |t: [u32; 3]| {
        let mut k = t;
        k.sort_unstable();
        k
    };
    let of = |o: u32| -> Vec<[u32; 3]> {
        a.triangles
            .iter()
            .zip(&a.origin)
            .filter(|(_, g)| **g == o)
            .map(|(t, _)| key(*t))
            .collect()
    };
    let (first, second) = (of(0), of(1));
    let overlap: f64 = first
        .iter()
        .filter(|k| second.contains(k))
        .map(|k| {
            let [x, y, z] = k.map(|i| a.points[i as usize].approx());
            area(x, y, z)
        })
        .sum();
    // The overlap is the triangle (1, 1), (3, 1), (1, 3): area 2.
    assert!((overlap - 2.0).abs() < 1e-12, "{overlap}");
}

#[test]
fn a_pushed_through_vertex_makes_a_closed_arrangement() {
    // An octahedron with its top vertex pushed through to other places.
    let mut p = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let f = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    let mut crossed = 0;
    for top in [
        [0.1, 0.05, -1.7],
        [1.5, 0.0, -0.2],
        [0.6, 0.7, -0.4],
        [-0.3, 1.4, -0.6],
        [0.2, 0.1, -0.5],
    ] {
        p[4] = top;
        let a = arrangement(&p, &f).unwrap();
        crossed += usize::from(a.pairs > 0);
        check(&p, &f, &a, true);
    }
    assert!(crossed >= 2, "{crossed}");
}

/// A deterministic generator for the random cases.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn tetrahedron(c: [Point3; 4]) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let n = normal(c[0], c[1], c[2]);
    let d = [c[3][0] - c[0][0], c[3][1] - c[0][1], c[3][2] - c[0][2]];
    // Corner 3 above the base 0 1 2 (counter-clockwise from above): the base
    // faces down.
    let f = if n[0] * d[0] + n[1] * d[1] + n[2] * d[2] > 0.0 {
        vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]
    } else {
        vec![[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]]
    };
    (c.to_vec(), f)
}

#[test]
fn random_overlapping_solids_make_closed_arrangements() {
    let mut rng = Lcg(16);
    let mut checked = 0;
    for case in 0..80 {
        let grid = case % 2 == 0;
        let mut parts = Vec::new();
        for _ in 0..3 {
            // On a coarse grid in half the cases, so coplanar and collinear
            // contacts happen.
            let mut q = || {
                let v = rng.next() * 4.0;
                if grid { (2.0 * v).round() / 2.0 } else { v }
            };
            if case % 3 == 0 {
                let lo = [q(), q(), q()];
                let hi = [lo[0] + 0.5 + q(), lo[1] + 0.5 + q(), lo[2] + 0.5 + q()];
                parts.push(cuboid(lo, hi, 0));
            } else {
                let c = [
                    [q(), q(), q()],
                    [q(), q(), q()],
                    [q(), q(), q()],
                    [q(), q(), q()],
                ];
                parts.push(tetrahedron(c));
            }
        }
        let (p, f) = union(&parts);
        let proper = f.iter().all(|t| {
            let [a, b, c] = t.map(|i| p[i as usize]);
            mesher_geom::predicates::plane_axis(a, b, c).is_some()
        });
        if !proper {
            continue;
        }
        // Weld coincident corners, so touching solids share vertices.
        let mut welded: Vec<Point3> = Vec::new();
        let ids: Vec<u32> = p
            .iter()
            .map(|q| match welded.iter().position(|w| w == q) {
                Some(i) => i as u32,
                None => {
                    welded.push(*q);
                    welded.len() as u32 - 1
                }
            })
            .collect();
        let f: Vec<[u32; 3]> = f.iter().map(|t| t.map(|v| ids[v as usize])).collect();
        let a = arrangement(&welded, &f).unwrap_or_else(|e| panic!("case {case}: {e}"));
        check(&welded, &f, &a, true);
        checked += 1;
    }
    assert!(checked > 40, "{checked}");
}

// Stage 3: the solid's boundary.

use mesher_geom::TriangleSoup;
use mesher_geom::arrangement::{Resolved, resolve, resolve_with};
use mesher_geom::bvh::closest_point_on_triangle;
use mesher_geom::winding::winding_number_brute;

fn volume(r: &Resolved) -> f64 {
    r.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| r.points[v as usize].approx());
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

fn balanced(triangles: &[[u32; 3]]) -> bool {
    let mut net: BTreeMap<(u32, u32), i64> = BTreeMap::new();
    for t in triangles {
        for k in 0..3 {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            *net.entry((u.min(v), u.max(v))).or_insert(0) += if u < v { 1 } else { -1 };
        }
    }
    net.values().all(|&n| n == 0)
}

#[test]
fn overlapping_boxes_resolve_to_their_union() {
    let (p, f) = union(&[
        cuboid([0.0; 3], [1.0; 3], 0),
        cuboid([0.5, 0.25, 0.25], [1.5, 1.25, 1.25], 0),
    ]);
    let r = resolve(&p, &f).unwrap();
    assert!(balanced(&r.triangles));
    // 1 + 1 − 0.5 × 0.75 × 0.75.
    assert!((volume(&r) - 1.71875).abs() < 1e-12, "{}", volume(&r));
}

#[test]
fn a_nested_box_is_dropped() {
    let (p, f) = union(&[cuboid([0.0; 3], [3.0; 3], 0), cuboid([1.0; 3], [2.0; 3], 0)]);
    let r = resolve(&p, &f).unwrap();
    assert!(balanced(&r.triangles));
    assert_eq!(r.triangles.len(), 12);
    assert!(r.origin.iter().all(|&o| o < 12));
    assert!((volume(&r) - 27.0).abs() < 1e-12);
}

#[test]
fn boxes_touching_face_to_face_become_one() {
    // Side by side, sharing the plane x = 1 over part of a face.
    let (p, f) = union(&[
        cuboid([0.0; 3], [1.0; 3], 0),
        cuboid([1.0, 0.25, 0.0], [2.0, 0.75, 1.0], 0),
    ]);
    let r = resolve(&p, &f).unwrap();
    assert!(balanced(&r.triangles));
    assert!((volume(&r) - 1.5).abs() < 1e-12, "{}", volume(&r));
    // Nothing is left in the shared region of x = 1.
    for t in &r.triangles {
        let c = t.map(|v| r.points[v as usize].approx());
        let on_interface = c
            .iter()
            .all(|q| q[0] == 1.0 && q[1] >= 0.25 && q[1] <= 0.75);
        assert!(!on_interface, "{c:?}");
    }
    assert!(r.dropped > 0);
}

#[test]
fn a_box_inside_out_within_another_is_removed() {
    // An inverted box inside a box: winding 0 inside it, so it is a cavity.
    let (p, mut f) = union(&[cuboid([0.0; 3], [3.0; 3], 0), cuboid([1.0; 3], [2.0; 3], 0)]);
    for t in &mut f[12..] {
        t.swap(1, 2);
    }
    let r = resolve(&p, &f).unwrap();
    assert!(balanced(&r.triangles));
    assert_eq!(r.triangles.len(), 24);
    assert!((volume(&r) - 26.0).abs() < 1e-12);
}

#[test]
fn an_open_surface_is_refused() {
    let (p, f) = cuboid([0.0; 3], [1.0; 3], 0);
    let f = f[1..].to_vec();
    assert!(matches!(
        resolve(&p, &f),
        Err(ArrangementError::Open { .. })
    ));
}

/// Away from both surfaces, the input wraps a point (winding number at
/// least a half) exactly when the resolved surface encloses it.
fn same_solid(p: &[Point3], f: &[[u32; 3]], r: &Resolved, rng: &mut Lcg, samples: usize) {
    let rp: Vec<Point3> = r.points.iter().map(|q| q.approx()).collect();
    let soup_in = TriangleSoup::new(p, f).unwrap();
    let soup_out = TriangleSoup::new(&rp, &r.triangles).unwrap();
    let (lo, hi) = p.iter().fold(
        ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
        |(lo, hi), q| {
            (
                [0, 1, 2].map(|k| lo[k].min(q[k])),
                [0, 1, 2].map(|k| hi[k].max(q[k])),
            )
        },
    );
    let near = |x: Point3, pts: &[Point3], tris: &[[u32; 3]]| {
        tris.iter().any(|t| {
            let c = closest_point_on_triangle(x, t.map(|v| pts[v as usize]));
            (0..3).map(|k| (c[k] - x[k]).powi(2)).sum::<f64>() < 1e-12
        })
    };
    for _ in 0..samples {
        let x = [0, 1, 2].map(|k| lo[k] - 0.1 + (hi[k] - lo[k] + 0.2) * rng.next());
        if near(x, p, f) || near(x, &rp, &r.triangles) {
            continue;
        }
        let w_in = winding_number_brute(&soup_in, x);
        let w_out = winding_number_brute(&soup_out, x);
        assert_eq!(
            w_in >= 0.5,
            w_out >= 0.5,
            "at {x:?}: {w_in} against {w_out}"
        );
        assert!(
            (w_out - w_out.round()).abs() < 1e-6 && w_out.round() <= 1.0,
            "{w_out}"
        );
    }
}

#[test]
fn random_overlapping_solids_resolve_to_their_union() {
    let mut rng = Lcg(61);
    let mut checked = 0;
    for case in 0..60 {
        let grid = case % 2 == 0;
        let mut parts = Vec::new();
        for _ in 0..3 {
            let mut q = || {
                let v = rng.next() * 4.0;
                if grid { (2.0 * v).round() / 2.0 } else { v }
            };
            if case % 3 == 0 {
                let lo = [q(), q(), q()];
                let hi = [lo[0] + 0.5 + q(), lo[1] + 0.5 + q(), lo[2] + 0.5 + q()];
                parts.push(cuboid(lo, hi, 0));
            } else {
                let c = [
                    [q(), q(), q()],
                    [q(), q(), q()],
                    [q(), q(), q()],
                    [q(), q(), q()],
                ];
                parts.push(tetrahedron(c));
            }
        }
        let (p, f) = union(&parts);
        let proper = f.iter().all(|t| {
            let [a, b, c] = t.map(|i| p[i as usize]);
            mesher_geom::predicates::plane_axis(a, b, c).is_some()
        });
        if !proper {
            continue;
        }
        let mut welded: Vec<Point3> = Vec::new();
        let ids: Vec<u32> = p
            .iter()
            .map(|q| match welded.iter().position(|w| w == q) {
                Some(i) => i as u32,
                None => {
                    welded.push(*q);
                    welded.len() as u32 - 1
                }
            })
            .collect();
        let f: Vec<[u32; 3]> = f.iter().map(|t| t.map(|v| ids[v as usize])).collect();
        let r = resolve(&welded, &f).unwrap_or_else(|e| panic!("case {case}: {e}"));
        assert!(balanced(&r.triangles), "case {case}");
        same_solid(&welded, &f, &r, &mut rng, 200);
        // Found around edges from one cast per piece, the same surface.
        let by_edges = resolve_with(&welded, &f, false)
            .unwrap_or_else(|e| panic!("case {case}, around edges: {e}"));
        assert_eq!(by_edges.triangles, r.triangles, "case {case}");
        checked += 1;
    }
    assert!(checked > 40, "{checked}");
}
