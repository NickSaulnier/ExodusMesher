//! Delaunay tetrahedralization against brute-force checks: valid structure,
//! every distinct point used, empty circumspheres, conforming faces, and a
//! closed convex hull. Inputs range from general position to exactly
//! cospherical and coplanar-heavy integer sets.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::{BTreeMap, BTreeSet};

use mesher_geom::predicates::{insphere, orient3d};
use mesher_geom::{Point3, Sign};
use mesher_tet::{Delaunay, DelaunayError};
use proptest::prelude::*;

/// Exactly coplanar (or fewer than four distinct points)?
fn flat(points: &[Point3]) -> bool {
    let distinct: Vec<Point3> = {
        let mut v = points.to_vec();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v.dedup();
        v
    };
    if distinct.len() < 4 {
        return true;
    }
    for a in 0..distinct.len() {
        for b in a + 1..distinct.len() {
            for c in b + 1..distinct.len() {
                for d in c + 1..distinct.len() {
                    if orient3d(distinct[a], distinct[b], distinct[c], distinct[d]) != Sign::Zero {
                        return false;
                    }
                }
            }
        }
    }
    true
}

fn verify(points: &[Point3], seed: u64) -> Result<(), TestCaseError> {
    let d = match Delaunay::new(points.to_vec(), seed) {
        Ok(d) => d,
        Err(DelaunayError::Flat) => {
            prop_assert!(flat(points), "reported flat but the points span a volume");
            return Ok(());
        }
        Err(e) => return Err(TestCaseError::fail(format!("{e}"))),
    };
    if let Err(e) = d.check() {
        return Err(TestCaseError::fail(e));
    }
    let tets = d.tets();

    // Every point is a vertex or an exact duplicate of one.
    let used: BTreeSet<u32> = tets.iter().flatten().copied().collect();
    let dups: BTreeMap<u32, u32> = d.duplicates().iter().copied().collect();
    for i in 0..points.len() as u32 {
        match dups.get(&i) {
            Some(&v) => prop_assert_eq!(points[i as usize], points[v as usize]),
            None => prop_assert!(used.contains(&i), "point {} is not a vertex", i),
        }
    }

    // Empty circumspheres.
    for t in &tets {
        let [a, b, c, e] = t.map(|v| points[v as usize]);
        prop_assert_eq!(orient3d(a, b, c, e), Sign::Positive);
        for &q in &used {
            if !t.contains(&q) {
                prop_assert!(
                    insphere(a, b, c, e, points[q as usize]) != Sign::Positive,
                    "point {} inside the circumsphere of {:?}",
                    q,
                    t
                );
            }
        }
    }

    // Conforming faces; the faces used once are exactly the hull, outward.
    let mut faces: BTreeMap<[u32; 3], Vec<[u32; 3]>> = BTreeMap::new();
    for t in &tets {
        for f in [
            [t[1], t[2], t[3]],
            [t[0], t[3], t[2]],
            [t[0], t[1], t[3]],
            [t[0], t[2], t[1]],
        ] {
            let mut k = f;
            k.sort_unstable();
            faces.entry(k).or_default().push(f);
        }
    }
    let mut boundary: BTreeSet<[u32; 3]> = BTreeSet::new();
    for uses in faces.values() {
        prop_assert!(uses.len() <= 2);
        if uses.len() == 1 {
            let f = uses[0];
            let r = (0..3).min_by_key(|&k| f[k]).unwrap();
            boundary.insert([f[r], f[(r + 1) % 3], f[(r + 2) % 3]]);
        }
    }
    let hull: BTreeSet<[u32; 3]> = d
        .hull()
        .into_iter()
        .map(|f| {
            let r = (0..3).min_by_key(|&k| f[k]).unwrap();
            [f[r], f[(r + 1) % 3], f[(r + 2) % 3]]
        })
        .collect();
    prop_assert_eq!(&boundary, &hull);
    for f in &hull {
        let [a, b, c] = f.map(|v| points[v as usize]);
        for &q in &used {
            prop_assert!(
                orient3d(a, b, c, points[q as usize]) != Sign::Positive,
                "point {} beyond hull face {:?}",
                q,
                f
            );
        }
    }
    Ok(())
}

/// Integer points with x² + y² + z² = 25: exactly cospherical.
fn sphere25() -> Vec<Point3> {
    let mut out = Vec::new();
    for x in -5i32..=5 {
        for y in -5i32..=5 {
            for z in -5i32..=5 {
                if x * x + y * y + z * z == 25 {
                    out.push([x, y, z].map(f64::from));
                }
            }
        }
    }
    out
}

#[test]
fn exactly_cospherical_points() {
    let pts = sphere25();
    assert_eq!(pts.len(), 30);
    for seed in 0..20 {
        verify(&pts, seed).unwrap();
    }
    // With the centre too.
    let mut with_centre = pts;
    with_centre.push([0.0; 3]);
    for seed in 0..20 {
        verify(&with_centre, seed).unwrap();
    }
}

#[test]
fn grid_on_a_cube_surface() {
    let mut pts = Vec::new();
    for x in 0..=3 {
        for y in 0..=3 {
            for z in 0..=3 {
                if [x, y, z].iter().any(|&c| c == 0 || c == 3) {
                    pts.push([x, y, z].map(f64::from));
                }
            }
        }
    }
    for seed in 0..20 {
        verify(&pts, seed).unwrap();
    }
}

#[test]
fn incremental_insertion_stays_valid() {
    let mut d = Delaunay::new(sphere25(), 3).unwrap();
    for p in [
        [0.0, 0.0, 0.0],
        [1.0, 2.0, 2.0],
        [9.0, 0.0, 0.0],
        [-9.0, -9.0, 9.0],
        [0.5, 0.25, 0.125],
    ] {
        d.insert(p).unwrap();
        d.check().unwrap();
    }
    let all = d.points().to_vec();
    // The incremental result and a fresh build are both valid Delaunay
    // tetrahedralizations of the same points.
    assert!(!Delaunay::new(all.clone(), 3).unwrap().tets().is_empty());
    verify(&all, 3).unwrap();
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(300))]

    #[test]
    fn random_points(pts in prop::collection::vec(prop::array::uniform3(-10.0f64..10.0), 4..150), seed in any::<u64>()) {
        verify(&pts, seed)?;
    }

    #[test]
    fn integer_grid_points(pts in prop::collection::vec(prop::array::uniform3((-2i32..=2).prop_map(f64::from)), 4..70), seed in any::<u64>()) {
        verify(&pts, seed)?;
    }

    #[test]
    fn nearly_flat_points(pts in prop::collection::vec((-5.0f64..5.0, -5.0f64..5.0, prop_oneof![Just(0.0), Just(1e-12), -1e-6f64..1e-6]), 4..60), seed in any::<u64>()) {
        let pts: Vec<Point3> = pts.into_iter().map(|(x, y, z)| [x, y, z]).collect();
        verify(&pts, seed)?;
    }
}
