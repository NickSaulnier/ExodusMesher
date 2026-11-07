//! Tolerance-based vertex clustering. Every referenced vertex maps to the
//! smallest vertex index within merge tolerance of it (transitively).
//!
//! Neighbour search uses a hash grid with cells of the tolerance size; the map
//! is only ever probed, never iterated, so results do not depend on hash order.

use std::collections::HashMap;

use mesher_geom::Point3;
use mesher_geom::vec3::distance_squared;

use super::uf::DisjointSets;

pub(crate) fn weld(points: &[Point3], referenced: &[bool], tolerance: f64) -> Vec<u32> {
    let n = points.len();
    let mut sets = DisjointSets::new(n);
    let used = (0..n as u32).filter(|&v| referenced[v as usize]);

    if tolerance > 0.0 && tolerance.is_finite() {
        let cell = |p: Point3| p.map(|c| (c / tolerance).floor() as i64);
        let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
        let tol2 = tolerance * tolerance;
        for v in used {
            let p = points[v as usize];
            let c = cell(p);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let key = [
                            c[0].saturating_add(dx),
                            c[1].saturating_add(dy),
                            c[2].saturating_add(dz),
                        ];
                        if let Some(members) = grid.get(&key) {
                            for &u in members {
                                if distance_squared(p, points[u as usize]) <= tol2 {
                                    sets.union(u, v);
                                }
                            }
                        }
                    }
                }
            }
            grid.entry(c).or_default().push(v);
        }
    } else {
        // Zero tolerance (all referenced vertices coincide, or none): merge
        // exactly equal positions only.
        let mut sorted: Vec<u32> = used.collect();
        let key = |v: u32| points[v as usize].map(|c| c + 0.0);
        sorted.sort_by(|&a, &b| {
            let (ka, kb) = (key(a), key(b));
            ka.iter()
                .zip(&kb)
                .map(|(x, y)| x.total_cmp(y))
                .find(|o| o.is_ne())
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.cmp(&b))
        });
        for w in sorted.windows(2) {
            if key(w[0]) == key(w[1]) {
                sets.union(w[0], w[1]);
            }
        }
    }
    (0..n as u32).map(|v| sets.find(v)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clusters_within_tolerance() {
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1e-9],
            [1.0, 0.0, 0.0],
            [5.0, 5.0, 5.0],
        ];
        let rep = weld(&pts, &[true, true, true, true, false], 1e-6);
        assert_eq!(rep, vec![0, 1, 0, 1, 4]);
    }

    #[test]
    fn zero_tolerance_merges_exact_copies() {
        let pts = [
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [-0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ];
        assert_eq!(weld(&pts, &[true; 4], 0.0), vec![0, 0, 2, 2]);
    }
}
