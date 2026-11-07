//! Milestone 7 exit criterion: on random point sets in general position, the
//! Delaunay tetrahedralization matches TetGen's exactly (it is unique there):
//! same element count, same tetrahedra.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeSet;
use std::time::Instant;

use mesher_geom::Point3;
use mesher_tet::delaunay;
use tritet::Tetgen;

/// Deterministic uniform points in the unit cube (SplitMix64).
fn random_points(n: usize, seed: u64) -> Vec<Point3> {
    let mut s = seed;
    let mut next = move || {
        s = s.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = s;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n).map(|_| [next(), next(), next()]).collect()
}

fn canonical(tets: impl IntoIterator<Item = [u32; 4]>) -> BTreeSet<[u32; 4]> {
    tets.into_iter()
        .map(|mut t| {
            t.sort_unstable();
            t
        })
        .collect()
}

fn tetgen(points: &[Point3]) -> BTreeSet<[u32; 4]> {
    let mut tg = Tetgen::new(points.len(), None, None, None).unwrap();
    for (i, p) in points.iter().enumerate() {
        tg.set_point(i, 0, p[0], p[1], p[2]).unwrap();
    }
    tg.generate_delaunay(false).unwrap();
    assert_eq!(
        tg.out_npoint(),
        points.len(),
        "TetGen renumbered or added points"
    );
    canonical((0..tg.out_ncell()).map(|c| [0, 1, 2, 3].map(|m| tg.out_cell_point(c, m) as u32)))
}

#[test]
fn matches_tetgen_on_random_points() {
    for (n, seed) in [
        (5, 1),
        (10, 2),
        (50, 3),
        (200, 4),
        (1000, 5),
        (5000, 6),
        (20_000, 7),
    ] {
        let points = random_points(n, seed);
        let t0 = Instant::now();
        let ours = canonical(delaunay(&points).unwrap().tets);
        let ours_time = t0.elapsed();
        let t1 = Instant::now();
        let theirs = tetgen(&points);
        let tetgen_time = t1.elapsed();
        eprintln!(
            "n = {n:>6}: {} tets ours ({ours_time:?}), {} TetGen ({tetgen_time:?})",
            ours.len(),
            theirs.len()
        );
        assert_eq!(ours.len(), theirs.len(), "element count differs at n = {n}");
        assert!(ours == theirs, "connectivity differs at n = {n}");
    }
}

#[test]
fn matches_tetgen_on_clustered_points() {
    // Points packed near a few centres, with widely varying scale.
    for seed in 10..20 {
        let base = random_points(400, seed);
        let points: Vec<Point3> = base
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let c = [(i % 4) as f64 * 10.0, 0.0, 0.0];
                let s = 10f64.powi(-((i % 3) as i32));
                [c[0] + s * p[0], c[1] + s * p[1], c[2] + s * p[2]]
            })
            .collect();
        assert_eq!(
            canonical(delaunay(&points).unwrap().tets),
            tetgen(&points),
            "seed {seed}"
        );
    }
}
