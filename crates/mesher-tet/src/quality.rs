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
}
