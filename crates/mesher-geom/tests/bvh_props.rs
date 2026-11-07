//! BVH queries against linear scans over every triangle.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

mod common;

use common::{float_soup, grid_cube, grid_point, grid_soup, to_f64};
use mesher_geom::bvh::{closest_point_on_triangle, ray_triangle};
use mesher_geom::vec3::distance_squared;
use mesher_geom::{Aabb, Bvh, Point3, Ray, TriangleSoup};
use proptest::prelude::*;

fn brute_overlaps(soup: &TriangleSoup<'_>) -> Vec<(u32, u32)> {
    let boxes: Vec<Aabb> = (0..soup.len())
        .map(|t| Aabb::of_points(&soup.triangle(t).unwrap()))
        .collect();
    let mut out = Vec::new();
    for i in 0..boxes.len() {
        for j in i + 1..boxes.len() {
            if boxes[i].overlaps(&boxes[j]) {
                out.push((i as u32, j as u32));
            }
        }
    }
    out
}

fn brute_ray(soup: &TriangleSoup<'_>, ray: &Ray) -> Option<(f64, u32)> {
    (0..soup.len())
        .filter_map(|t| ray_triangle(ray, soup.triangle(t).unwrap(), t as u32))
        .map(|h| (h.t, h.triangle))
        .min_by(|a, b| a.partial_cmp(b).unwrap())
}

fn brute_nearest(soup: &TriangleSoup<'_>, p: Point3) -> Option<f64> {
    (0..soup.len())
        .map(|t| distance_squared(p, closest_point_on_triangle(p, soup.triangle(t).unwrap())))
        .min_by(f64::total_cmp)
}

fn check(
    points: &[Point3],
    tris: &[[u32; 3]],
    probes: &[(Point3, Point3)],
) -> Result<(), TestCaseError> {
    let soup = TriangleSoup::new(points, tris).unwrap();
    let bvh = Bvh::new(soup);
    prop_assert_eq!(bvh.self_overlaps(), brute_overlaps(&soup));
    let scale = bvh.bounds().diagonal() + 1.0;
    for &(origin, dir) in probes {
        let ray = Ray { origin, dir };
        let fast = bvh.ray_cast(&ray).map(|h| (h.t, h.triangle));
        let slow = brute_ray(&soup, &ray);
        match (fast, slow) {
            (Some(f), Some(s)) => prop_assert!(
                f == s || (f.0 - s.0).abs() <= 1e-9 * (1.0 + s.0.abs()),
                "{f:?} vs {s:?}"
            ),
            (f, s) => prop_assert_eq!(f, s),
        }
        let fast = bvh.nearest(origin).map(|n| n.distance_squared);
        let slow = brute_nearest(&soup, origin);
        match (fast, slow) {
            (Some(f), Some(s)) => {
                prop_assert!((f.sqrt() - s.sqrt()).abs() <= 1e-9 * scale, "{f} vs {s}")
            }
            (f, s) => prop_assert_eq!(f, s),
        }
    }
    Ok(())
}

fn probes(range: f64) -> impl Strategy<Value = Vec<(Point3, Point3)>> {
    prop::collection::vec(
        (
            prop::array::uniform3(-range..range),
            prop::array::uniform3(-1.0f64..1.0),
        ),
        1..8,
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(400))]

    #[test]
    fn random_float_soups((pts, tris) in float_soup(40, 120), probes in probes(15.0)) {
        check(&pts, &tris, &probes)?;
    }

    #[test]
    fn degenerate_grid_soups(
        (pts, tris) in grid_soup(30, 80),
        grid_probes in prop::collection::vec((grid_point(), grid_point()), 1..8),
    ) {
        // Rays from grid points along grid directions hit vertices and edges exactly.
        let probes: Vec<_> = grid_probes.into_iter().map(|(o, d)| (to_f64(o), to_f64(d))).collect();
        check(&pts, &tris, &probes)?;
    }

    #[test]
    fn grid_cube_surfaces(n in 1i64..5, probes in probes(6.0)) {
        let (pts, tris) = grid_cube(n);
        check(&pts, &tris, &probes)?;
    }

    #[test]
    fn overlaps_between_two_soups((pa, ta) in grid_soup(20, 40), (pb, tb) in float_soup(20, 40)) {
        let pb: Vec<Point3> = pb.iter().map(|p| p.map(|c| c / 3.0)).collect();
        let (sa, sb) = (TriangleSoup::new(&pa, &ta).unwrap(), TriangleSoup::new(&pb, &tb).unwrap());
        let mut expected = Vec::new();
        for i in 0..sa.len() {
            for j in 0..sb.len() {
                let (bi, bj) = (Aabb::of_points(&sa.triangle(i).unwrap()), Aabb::of_points(&sb.triangle(j).unwrap()));
                if bi.overlaps(&bj) {
                    expected.push((i as u32, j as u32));
                }
            }
        }
        prop_assert_eq!(Bvh::new(sa).overlaps_with(&Bvh::new(sb)), expected);
    }
}
