//! Exact intersection tests against the linear-programming reference in
//! `common::basic_feasible_solutions`, on integer grids where coincident,
//! collinear and coplanar configurations are common.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

mod common;

use common::{basic_feasible_solutions, cross_i, grid_point, to_f64};
use mesher_geom::{
    faces_intersect, segment_intersects_triangle, segments_intersect, triangles_intersect,
};
use proptest::prelude::*;

fn tri_f(t: [[i64; 3]; 3]) -> [[f64; 3]; 3] {
    t.map(to_f64)
}

fn proper(t: [[i64; 3]; 3]) -> bool {
    cross_i(t[0], t[1], t[2]) != [0, 0, 0]
}

/// Three distinct indices into a pool of 5.
fn face_ids() -> impl Strategy<Value = [u32; 3]> {
    Just(()).prop_perturb(|(), mut rng| {
        let mut ids = [0u32, 1, 2, 3, 4];
        for i in (1..5).rev() {
            let j = rng.random_range(0..=i);
            ids.swap(i, j);
        }
        [ids[0], ids[1], ids[2]]
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(3000))]

    #[test]
    fn triangles_match_oracle(a in prop::array::uniform3(grid_point()), b in prop::array::uniform3(grid_point())) {
        let expected = !basic_feasible_solutions(a, b).is_empty();
        prop_assert_eq!(triangles_intersect(tri_f(a), tri_f(b)), expected);
        prop_assert_eq!(triangles_intersect(tri_f(b), tri_f(a)), expected);
        let t = |tri: [[i64; 3]; 3]| tri_f(tri).map(|p| p.map(|x| x * 2f64.powi(-40) + 1024.0));
        prop_assert_eq!(triangles_intersect(t(a), t(b)), expected);
    }

    #[test]
    fn segment_triangle_matches_oracle(p in grid_point(), q in grid_point(), b in prop::array::uniform3(grid_point())) {
        let expected = !basic_feasible_solutions([p, q, q], b).is_empty();
        prop_assert_eq!(segment_intersects_triangle(to_f64(p), to_f64(q), tri_f(b)), expected);
    }

    #[test]
    fn segments_match_oracle(p in grid_point(), q in grid_point(), x in grid_point(), y in grid_point()) {
        let expected = !basic_feasible_solutions([p, q, q], [x, y, y]).is_empty();
        prop_assert_eq!(segments_intersect(to_f64(p), to_f64(q), to_f64(x), to_f64(y)), expected);
    }

    #[test]
    fn mesh_faces_match_oracle(
        pool in prop::array::uniform5(grid_point()),
        ia in face_ids(),
        ib in face_ids(),
    ) {
        let a = ia.map(|i| pool[i as usize]);
        let b = ib.map(|i| pool[i as usize]);
        let bfs = basic_feasible_solutions(a, b);
        let shared: Vec<usize> = (0..3).filter(|&i| ib.contains(&ia[i])).collect();
        let expected = if !(proper(a) && proper(b)) || shared.is_empty() {
            // Documented fallback: plain closed intersection.
            !bfs.is_empty()
        } else if shared.len() == 1 {
            // Some common point other than the shared vertex: λ_v < 1.
            bfs.iter().any(|x| x[shared[0]].n < x[shared[0]].d)
        } else if shared.len() == 2 {
            // Some common point off the shared edge: λ_third > 0.
            let third = (0..3).find(|i| !shared.contains(i)).unwrap();
            bfs.iter().any(|x| x[third].n > 0)
        } else {
            true
        };
        prop_assert_eq!(faces_intersect(ia, tri_f(a), ib, tri_f(b)), expected);
        prop_assert_eq!(faces_intersect(ib, tri_f(b), ia, tri_f(a)), expected);
    }
}
