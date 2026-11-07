//! The hierarchical winding number against the direct sum, and both against
//! known inside/outside answers on closed surfaces.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

mod common;

use common::{float_soup, grid_cube, grid_point, icosphere, to_f64};
use mesher_geom::{FastWindingNumber, Point3, TriangleSoup, winding_number_brute};
use proptest::prelude::*;

fn agree(points: &[Point3], tris: &[[u32; 3]], queries: &[Point3]) -> Result<(), TestCaseError> {
    let soup = TriangleSoup::new(points, tris).unwrap();
    let fast = FastWindingNumber::new(soup);
    let tol = 1e-10 * (tris.len() as f64 + 1.0);
    for &p in queries {
        let (f, b) = (fast.winding_number(p), winding_number_brute(&soup, p));
        prop_assert!((f - b).abs() <= tol, "at {p:?}: fast {f} vs brute {b}");
    }
    Ok(())
}

#[test]
fn closed_surfaces_give_integers() {
    for (pts, tris, inside, outside) in [
        {
            let (p, t) = icosphere(2);
            (p, t, [0.1, -0.2, 0.3], [1.5, 0.2, -0.1])
        },
        {
            let (p, t) = grid_cube(3);
            (p, t, [1.5, 1.2, 0.7], [-0.5, 1.0, 1.0])
        },
    ] {
        let soup = TriangleSoup::new(&pts, &tris).unwrap();
        let fast = FastWindingNumber::new(soup);
        for (p, expected) in [(inside, 1.0), (outside, 0.0)] {
            assert!((winding_number_brute(&soup, p) - expected).abs() < 1e-10);
            assert!((fast.winding_number(p) - expected).abs() < 1e-10);
        }
    }
}

#[test]
fn reversed_surface_gives_minus_one() {
    let (pts, tris) = icosphere(1);
    let flipped: Vec<[u32; 3]> = tris.iter().map(|&[a, b, c]| [a, c, b]).collect();
    let fast = FastWindingNumber::new(TriangleSoup::new(&pts, &flipped).unwrap());
    assert!((fast.winding_number([0.0; 3]) + 1.0).abs() < 1e-10);
}

#[test]
fn open_hemisphere_is_fractional() {
    let (pts, tris) = icosphere(2);
    // Keep triangles with all corners at z >= 0: an open cap.
    let cap: Vec<[u32; 3]> = tris
        .into_iter()
        .filter(|t| t.iter().all(|&v| pts[v as usize][2] >= 0.0))
        .collect();
    let soup = TriangleSoup::new(&pts, &cap).unwrap();
    let w = FastWindingNumber::new(soup).winding_number([0.0, 0.0, 0.0]);
    assert!(w > 0.3 && w < 0.7, "{w}");
    assert!((w - winding_number_brute(&soup, [0.0; 3])).abs() < 1e-10);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(300))]

    #[test]
    fn icosphere_queries(sub in 0u32..3, qs in prop::collection::vec(prop::array::uniform3(-2.0f64..2.0), 1..20)) {
        let (pts, tris) = icosphere(sub);
        agree(&pts, &tris, &qs)?;
        // Away from the surface the answer is the indicator of the polyhedron.
        let soup = TriangleSoup::new(&pts, &tris).unwrap();
        let fast = FastWindingNumber::new(soup);
        for p in qs {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            if r < 0.7 {
                prop_assert!((fast.winding_number(p) - 1.0).abs() < 1e-9);
            } else if r > 1.05 {
                prop_assert!(fast.winding_number(p).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn grid_cube_queries_on_and_off_the_surface(n in 1i64..5, qs in prop::collection::vec(grid_point(), 1..30)) {
        // Integer queries land on vertices, edges and faces of the surface.
        let (pts, tris) = grid_cube(n);
        let queries: Vec<Point3> = qs.into_iter().map(|q| to_f64(q).map(|c| c + 1.0)).collect();
        agree(&pts, &tris, &queries)?;
    }

    #[test]
    fn random_soups((pts, tris) in float_soup(30, 80), qs in prop::collection::vec(prop::array::uniform3(-30.0f64..30.0), 1..20)) {
        agree(&pts, &tris, &qs)?;
    }
}
