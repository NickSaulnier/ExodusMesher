//! Element quality measures. Milestone 6 needs dihedral angles; the
//! finite-volume metrics (non-orthogonality, skewness, volume ratio) follow
//! in Milestone 8.

use mesher_geom::Point3;
use mesher_geom::vec3::{cross, dot, norm, scale, sub};

/// The six dihedral angles of a tetrahedron, in degrees, at edges
/// (0,1), (0,2), (0,3), (1,2), (1,3), (2,3).
#[must_use]
pub fn dihedral_angles(p: [Point3; 4]) -> [f64; 6] {
    const EDGES: [(usize, usize, usize, usize); 6] = [
        (0, 1, 2, 3),
        (0, 2, 1, 3),
        (0, 3, 1, 2),
        (1, 2, 0, 3),
        (1, 3, 0, 2),
        (2, 3, 0, 1),
    ];
    EDGES.map(|(i, j, k, l)| {
        let e = sub(p[j], p[i]);
        let ee = dot(e, e);
        let perp = |q: Point3| {
            let d = sub(q, p[i]);
            if ee > 0.0 {
                sub(d, scale(e, dot(d, e) / ee))
            } else {
                d
            }
        };
        let (u, v) = (perp(p[k]), perp(p[l]));
        norm(cross(u, v)).atan2(dot(u, v)).to_degrees()
    })
}

/// Smallest and largest dihedral angle over a mesh, in degrees; `None` for no tets.
#[must_use]
pub fn dihedral_range(points: &[Point3], tets: &[[u32; 4]]) -> Option<(f64, f64)> {
    tets.iter()
        .filter_map(|t| {
            let p = t.map(|v| points.get(v as usize).copied());
            let [Some(a), Some(b), Some(c), Some(d)] = p else {
                return None;
            };
            let angles = dihedral_angles([a, b, c, d]);
            Some((
                angles.iter().copied().fold(f64::INFINITY, f64::min),
                angles.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            ))
        })
        .reduce(|x, y| (x.0.min(y.0), x.1.max(y.1)))
}

/// Finite-volume quality of a tetrahedral mesh: what a cell-centred solver
/// sees. Non-orthogonality and skewness are measured on interior faces.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FvQuality {
    pub cells: usize,
    pub min_dihedral_degrees: f64,
    pub max_dihedral_degrees: f64,
    /// Every dihedral angle (six per tet) in 10° bins, [0°, 10°) to [170°, 180°].
    pub dihedral_histogram: Vec<u64>,
    /// Angle between a face normal and the line joining the two cell centroids.
    pub max_non_orthogonality_degrees: f64,
    pub mean_non_orthogonality_degrees: f64,
    /// Distance from the face centroid to where the centroid line crosses the
    /// face, over the centroid distance.
    pub max_skewness: f64,
    /// Largest volume ratio between face neighbours whose floating-point
    /// volumes are both positive; always finite.
    pub max_volume_ratio: f64,
    /// Tets whose volume, computed in floating point as a solver would, is
    /// not positive. Exact predicates can find such a tet positive (it is
    /// flat to rounding), but a finite-volume solver cannot use it.
    pub flat_in_float: usize,
}

/// The signed volume of a tet in floating point, the way a solver computes it.
#[must_use]
pub fn float_volume(p: [Point3; 4]) -> f64 {
    dot(sub(p[1], p[0]), cross(sub(p[2], p[0]), sub(p[3], p[0]))) / 6.0
}

/// Faces of a positive tet, outward.
const OUTWARD: [[usize; 3]; 4] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

#[must_use]
pub fn fv_quality(points: &[Point3], tets: &[[u32; 4]]) -> FvQuality {
    let corners = |t: &[u32; 4]| t.map(|v| points.get(v as usize).copied().unwrap_or([0.0; 3]));
    let centroid =
        |p: [Point3; 4]| [0, 1, 2].map(|k| (p[0][k] + p[1][k] + p[2][k] + p[3][k]) / 4.0);
    let volume = float_volume;

    let mut q = FvQuality {
        cells: tets.len(),
        min_dihedral_degrees: f64::INFINITY,
        max_dihedral_degrees: f64::NEG_INFINITY,
        dihedral_histogram: vec![0; 18],
        ..FvQuality::default()
    };
    for t in tets {
        let v = volume(corners(t));
        if v <= 0.0 || v.is_nan() {
            q.flat_in_float += 1;
        }
        for a in dihedral_angles(corners(t)) {
            q.min_dihedral_degrees = q.min_dihedral_degrees.min(a);
            q.max_dihedral_degrees = q.max_dihedral_degrees.max(a);
            let bin = ((a / 10.0).floor().max(0.0) as usize).min(17);
            q.dihedral_histogram[bin] += 1;
        }
    }
    if tets.is_empty() {
        q.min_dihedral_degrees = 0.0;
        q.max_dihedral_degrees = 0.0;
        return q;
    }

    let mut faces: std::collections::HashMap<[u32; 3], (usize, usize)> =
        std::collections::HashMap::with_capacity(tets.len() * 2);
    let (mut sum, mut count) = (0.0f64, 0usize);
    for (i, t) in tets.iter().enumerate() {
        for (k, f) in OUTWARD.iter().enumerate() {
            let mut key = f.map(|j| t[j]);
            key.sort_unstable();
            let Some((j, kj)) = faces.remove(&key) else {
                faces.insert(key, (i, k));
                continue;
            };
            // Interior face between cell j (first seen) and cell i.
            let (pj, pi) = (corners(&tets[j]), corners(t));
            let (cj, ci) = (centroid(pj), centroid(pi));
            let face = OUTWARD[kj].map(|m| pj[m]);
            let fc = [0, 1, 2].map(|m| (face[0][m] + face[1][m] + face[2][m]) / 3.0);
            let n = cross(sub(face[1], face[0]), sub(face[2], face[0]));
            let nn = norm(n);
            let d = sub(ci, cj);
            let dd = norm(d);
            if nn <= 0.0 || dd <= 0.0 {
                continue;
            }
            let n = scale(n, 1.0 / nn);
            let cos = (dot(d, n) / dd).clamp(-1.0, 1.0);
            let angle = cos.acos().to_degrees();
            q.max_non_orthogonality_degrees = q.max_non_orthogonality_degrees.max(angle);
            sum += angle;
            count += 1;
            let skew = if dot(d, n) > 0.0 {
                let t = dot(sub(fc, cj), n) / dot(d, n);
                let x = [0, 1, 2].map(|m| cj[m] + t * d[m]);
                norm(sub(fc, x)) / dd
            } else {
                norm(sub(fc, cj)) / dd
            };
            q.max_skewness = q.max_skewness.max(skew);
            // Flat cells are counted in `flat_in_float`; a ratio against one
            // would be unbounded.
            let (vj, vi) = (volume(pj), volume(pi));
            if vj > 0.0 && vi > 0.0 {
                let ratio = (vj.max(vi) / vj.min(vi)).min(f64::MAX);
                q.max_volume_ratio = q.max_volume_ratio.max(ratio);
            }
        }
    }
    q.mean_non_orthogonality_degrees = if count > 0 { sum / count as f64 } else { 0.0 };
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regular_tetrahedron() {
        let p = [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        let expected = (1.0f64 / 3.0).acos().to_degrees();
        for a in dihedral_angles(p) {
            assert!((a - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn corner_tetrahedron() {
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let a = dihedral_angles(p);
        // Edges at the right-angled corner are 90°.
        assert!((a[0] - 90.0).abs() < 1e-12);
        assert!((a[1] - 90.0).abs() < 1e-12);
        assert!((a[2] - 90.0).abs() < 1e-12);
    }

    #[test]
    fn fv_quality_of_two_mirror_tets() {
        // Two tets mirrored through the face x + y + z = 1 opposite the
        // corner: the centroid line is perpendicular to the shared face.
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
        ];
        let tets = [[0, 1, 2, 3], [4, 1, 3, 2]];
        let q = fv_quality(&p, &tets);
        assert_eq!(q.cells, 2);
        assert!(
            q.max_non_orthogonality_degrees < 1e-6,
            "{}",
            q.max_non_orthogonality_degrees
        );
        assert!(q.max_skewness < 1e-9);
        assert!((q.max_volume_ratio - 1.0).abs() < 1e-12);
        assert_eq!(q.dihedral_histogram.iter().sum::<u64>(), 12);
        assert_eq!(q.flat_in_float, 0);
    }

    /// A neighbour flat in floating point is counted, not divided by: the
    /// volume ratio stays finite (an infinite one cannot be written as JSON).
    #[test]
    fn flat_neighbour_is_counted_and_ratio_stays_finite() {
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            // Exactly on the shared face's plane x + y + z = 1: a flat tet.
            [0.5, 0.5, 0.0],
        ];
        let tets = [[0, 1, 2, 3], [4, 1, 3, 2]];
        let q = fv_quality(&p, &tets);
        assert_eq!(q.flat_in_float, 1);
        assert!(q.max_volume_ratio.is_finite());
    }
}
