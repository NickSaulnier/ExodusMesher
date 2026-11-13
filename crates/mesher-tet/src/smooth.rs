//! Sliver repair by moving free vertices of a finished tetrahedral mesh.
//!
//! A vertex not on the boundary may move anywhere its surrounding tets stay
//! positive, so each free vertex of a poor tet is moved to the candidate
//! position that most raises the smallest dihedral angle among its tets:
//! the centroid of its neighbours (and points towards it), steps along the
//! direction that improves the worst angle (a finite-difference gradient),
//! or steps away from the face opposite it in its worst tet. A move is kept only if every surrounding tet stays positively
//! oriented (exact `orient3d`) and the worst angle improves. Connectivity
//! and the boundary never change, so a conformal mesh stays conformal and
//! its volume is unchanged. Deterministic: fixed order, no randomness.

use mesher_geom::predicates::orient3d;
use mesher_geom::vec3::{add, cross, dot, norm, scale, sub};
use mesher_geom::{Point3, Sign};

/// How smoothing works.
#[derive(Debug, Clone, PartialEq)]
pub struct SmoothParams {
    /// Free vertices of tets whose smallest dihedral angle is below this are
    /// moved, degrees.
    pub target_degrees: f64,
    /// Passes over the poor tets; stops early when nothing moves.
    pub passes: usize,
}

impl Default for SmoothParams {
    fn default() -> SmoothParams {
        SmoothParams {
            target_degrees: 15.0,
            passes: 4,
        }
    }
}

/// What smoothing did.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SmoothStats {
    /// Vertex moves kept.
    pub moves: usize,
    /// Passes run.
    pub passes: usize,
}

/// The smallest dihedral angle of a tet, degrees: the angle at an edge is
/// π minus the angle between the outward normals of its two faces, so the
/// smallest angle comes from the largest cosine between normals. Cheaper
/// than `dihedral_angles` (one `acos`, not six), which matters here: every
/// candidate position re-evaluates every tet around the vertex.
pub(crate) fn min_dihedral(p: [Point3; 4]) -> f64 {
    // Outward normals of a positive tet (all inward for a negative one,
    // which leaves every angle between them the same). They are oriented by
    // the vertex order, not by which side of its face the opposite vertex
    // lies: for a tet flat in floating point that side is noise, and a flat
    // tet would come out at 180°, not 0°.
    let normals = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]].map(|[a, b, c]| {
        let n = cross(sub(p[b], p[a]), sub(p[c], p[a]));
        let l = norm(n);
        (l > 0.0).then(|| scale(n, 1.0 / l))
    });
    let mut max_cos = -1.0f64;
    for i in 0..4 {
        for j in i + 1..4 {
            let (Some(a), Some(b)) = (normals[i], normals[j]) else {
                return 0.0;
            };
            max_cos = max_cos.max(-dot(a, b));
        }
    }
    libm::acos(max_cos.clamp(-1.0, 1.0)).to_degrees()
}

/// Move free vertices (`fixed[v]` false) of poor tets to raise the smallest
/// dihedral angle around them. `tets` must be positively oriented.
pub fn smooth_slivers(
    points: &mut [Point3],
    tets: &[[u32; 4]],
    fixed: &[bool],
    params: &SmoothParams,
) -> SmoothStats {
    let _prof = crate::prof::scope("smooth::smooth_slivers");
    let mut star: Vec<Vec<u32>> = vec![Vec::new(); points.len()];
    for (t, tet) in tets.iter().enumerate() {
        for &v in tet {
            star[v as usize].push(t as u32);
        }
    }
    let free = |v: u32| !fixed.get(v as usize).copied().unwrap_or(true);
    let mut stats = SmoothStats::default();
    for _ in 0..params.passes {
        stats.passes += 1;
        let mut candidates: Vec<u32> = tets
            .iter()
            .filter(|t| min_dihedral(t.map(|v| points[v as usize])) < params.target_degrees)
            .flat_map(|t| t.iter().copied().filter(|&v| free(v)))
            .collect();
        candidates.sort_unstable();
        candidates.dedup();
        let mut moved = 0;
        for v in candidates {
            if improve(points, tets, &star[v as usize], v) {
                moved += 1;
            }
        }
        stats.moves += moved;
        if moved == 0 {
            break;
        }
    }
    stats
}

/// A tet's quality for a vertex move, or `None` if the tet is not usable.
type Quality<'a> = &'a dyn Fn([Point3; 4]) -> Option<f64>;

/// The worst quality around `v` with `v` at `x`, or `None` if a tet would
/// not be usable.
fn quality_at(
    points: &[Point3],
    tets: &[[u32; 4]],
    star: &[u32],
    v: u32,
    x: Point3,
    quality: Quality<'_>,
) -> Option<f64> {
    let mut worst = 180.0f64;
    for &t in star {
        let p = tets[t as usize].map(|w| if w == v { x } else { points[w as usize] });
        worst = worst.min(quality(p)?);
    }
    Some(worst)
}

/// Positive tets, by their smallest dihedral angle.
fn positive(p: [Point3; 4]) -> Option<f64> {
    (orient3d(p[0], p[1], p[2], p[3]) == Sign::Positive).then(|| min_dihedral(p))
}

/// Try candidate positions for `v`; keep the best if it improves.
fn improve(points: &mut [Point3], tets: &[[u32; 4]], star: &[u32], v: u32) -> bool {
    let at = |x: Point3| quality_at(points, tets, star, v, x, &positive);
    let pos = |w: u32| points[w as usize];
    match best_position(&pos, tets, star, v, &at) {
        Some((q0, q, x)) if q > q0 + 1e-9 => {
            points[v as usize] = x;
            true
        }
        _ => false,
    }
}

/// The quality around a vertex with it at a position, or `None` if a tet
/// there would not be usable.
pub(crate) type At<'a> = &'a dyn Fn(Point3) -> Option<f64>;

/// The quality around `v` where it is (negative infinity if a tet there is
/// not usable), and the best of the candidate positions with its quality
/// (`v`'s own position if none is better), by `at`. `None` if `v` has no
/// neighbours. `star` indexes `tets`; every tet in it has `v`; `pos` gives
/// each vertex's position.
pub(crate) fn best_position(
    pos: &dyn Fn(u32) -> Point3,
    tets: &[[u32; 4]],
    star: &[u32],
    v: u32,
    at: At<'_>,
) -> Option<(f64, f64, Point3)> {
    let p = pos(v);
    let q0 = at(p).unwrap_or(f64::NEG_INFINITY);
    // Neighbours and the local edge length.
    let mut neighbours: Vec<u32> = star
        .iter()
        .flat_map(|&t| tets[t as usize])
        .filter(|&w| w != v)
        .collect();
    neighbours.sort_unstable();
    neighbours.dedup();
    if neighbours.is_empty() {
        return None;
    }
    let length = neighbours
        .iter()
        .map(|&w| norm(sub(pos(w), p)))
        .sum::<f64>()
        / neighbours.len() as f64;
    if !(length > 0.0 && length.is_finite()) {
        return None;
    }
    let mut centroid = [0.0; 3];
    for &w in &neighbours {
        centroid = add(centroid, pos(w));
    }
    centroid = scale(centroid, 1.0 / neighbours.len() as f64);

    let mut candidates: Vec<Point3> = [1.0, 0.5, 0.25]
        .iter()
        .map(|&s| add(p, scale(sub(centroid, p), s)))
        .collect();
    // Finite-difference gradient of the worst angle.
    let h = 1e-3 * length;
    let mut grad = [0.0; 3];
    for k in 0..3 {
        let (mut a, mut b) = (p, p);
        a[k] += h;
        b[k] -= h;
        if let (Some(qa), Some(qb)) = (at(a), at(b)) {
            grad[k] = (qa - qb) / (2.0 * h);
        }
    }
    let g = norm(grad);
    if g > 0.0 && g.is_finite() {
        let dir = scale(grad, 1.0 / g);
        for s in [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625] {
            candidates.push(add(p, scale(dir, s * length)));
        }
    }
    // Away from the face opposite `v` in its worst tet: a sliver standing on
    // a face with its apex nearly in the face's plane.
    if let Some(&worst) = star.iter().min_by(|&&a, &&b| {
        let q = |t: u32| min_dihedral(tets[t as usize].map(pos));
        q(a).total_cmp(&q(b)).then(a.cmp(&b))
    }) {
        let others: Vec<Point3> = tets[worst as usize]
            .iter()
            .filter(|&&w| w != v)
            .map(|&w| pos(w))
            .collect();
        let n = cross(sub(others[1], others[0]), sub(others[2], others[0]));
        let l = norm(n);
        if l > 0.0 {
            let side = if dot(n, sub(p, others[0])) >= 0.0 {
                1.0
            } else {
                -1.0
            };
            let away = scale(n, side / l);
            for s in [0.5, 0.25, 0.125, 0.0625, 0.03125] {
                candidates.push(add(p, scale(away, s * length)));
            }
        }
    }
    let mut best = (q0, p);
    for c in candidates {
        if let Some(q) = at(c) {
            if q > best.0 {
                best = (q, c);
            }
        }
    }
    Some((q0, best.0, best.1))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;

    /// A unit cube split into tets around a centre vertex, with the centre
    /// pushed almost onto one face: the tets on that face become slivers.
    fn cube_with_centre(centre: Point3) -> (Vec<Point3>, Vec<[u32; 4]>) {
        let mut p: Vec<Point3> = (0..8)
            .map(|i| {
                [
                    f64::from(i & 1),
                    f64::from((i >> 1) & 1),
                    f64::from((i >> 2) & 1),
                ]
            })
            .collect();
        p.push(centre);
        let c = 8;
        let quads = [
            [0, 2, 3, 1],
            [4, 5, 7, 6],
            [0, 1, 5, 4],
            [2, 6, 7, 3],
            [0, 4, 6, 2],
            [1, 3, 7, 5],
        ];
        let mut t = Vec::new();
        for [a, b, cc, d] in quads {
            for tri in [[a, b, cc], [a, cc, d]] {
                // Outward triangle; the tet with the centre is positive when
                // the centre is behind it.
                let mut tet = [tri[0], tri[1], tri[2], c];
                let q = tet.map(|v| p[v as usize]);
                if orient3d(q[0], q[1], q[2], q[3]) == Sign::Negative {
                    tet.swap(0, 1);
                }
                t.push(tet);
            }
        }
        (p, t)
    }

    fn worst(p: &[Point3], t: &[[u32; 4]]) -> f64 {
        t.iter()
            .map(|t| min_dihedral(t.map(|v| p[v as usize])))
            .fold(180.0, f64::min)
    }

    fn volume(p: &[Point3], t: &[[u32; 4]]) -> f64 {
        t.iter()
            .map(|t| {
                let [a, b, c, d] = t.map(|v| p[v as usize]);
                let (u, v, w) = (sub(b, a), sub(c, a), sub(d, a));
                (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                    + u[2] * (v[0] * w[1] - v[1] * w[0]))
                    / 6.0
            })
            .sum()
    }

    #[test]
    fn a_free_vertex_moves_out_of_a_sliver() {
        let (mut p, t) = cube_with_centre([0.5, 0.5, 0.002]);
        let before = worst(&p, &t);
        assert!(before < 2.0, "{before}");
        let fixed: Vec<bool> = (0..9).map(|v| v < 8).collect();
        let stats = smooth_slivers(&mut p, &t, &fixed, &SmoothParams::default());
        let after = worst(&p, &t);
        assert!(stats.moves > 0);
        assert!(after > 20.0, "{before} -> {after}");
        assert!((volume(&p, &t) - 1.0).abs() < 1e-12);
        let original = cube_with_centre([0.5; 3]).0;
        assert_eq!(p[..8], original[..8], "a fixed vertex moved");
    }

    #[test]
    fn the_fast_minimum_matches_the_dihedral_angles() {
        let mut seed = 7u64;
        let mut next = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..1000 {
            let p: [Point3; 4] = std::array::from_fn(|_| [next(), next(), next()]);
            let slow = crate::quality::dihedral_angles(p)
                .into_iter()
                .fold(180.0, f64::min);
            assert!((min_dihedral(p) - slow).abs() < 1e-6, "{p:?}");
        }
    }

    /// A tet flat in floating point (four points in a vertical plane,
    /// positive exactly) has a dihedral angle of 0°, not 180°.
    #[test]
    fn a_flat_tet_has_no_angle() {
        let p = [
            [-0.0180615, -5.61284, 0.41910145931375087],
            [-0.010918675000000001, -5.5563, 0.3841763377042716],
            [-0.00377585, -5.49976, 0.4540265809232301],
            [-0.00377585, -5.49976, 0.4365640201184905],
        ];
        assert_eq!(orient3d(p[0], p[1], p[2], p[3]), Sign::Positive);
        assert!(min_dihedral(p) < 1e-3, "{}", min_dihedral(p));
    }

    #[test]
    fn fixed_vertices_never_move() {
        let (mut p, t) = cube_with_centre([0.5, 0.5, 0.002]);
        let original = p.clone();
        let stats = smooth_slivers(&mut p, &t, &[true; 9], &SmoothParams::default());
        assert_eq!(stats.moves, 0);
        assert_eq!(p, original);
    }

    #[test]
    fn smoothing_is_deterministic() {
        let (mut a, t) = cube_with_centre([0.31, 0.77, 0.004]);
        let mut b = a.clone();
        let fixed: Vec<bool> = (0..9).map(|v| v < 8).collect();
        smooth_slivers(&mut a, &t, &fixed, &SmoothParams::default());
        smooth_slivers(&mut b, &t, &fixed, &SmoothParams::default());
        assert_eq!(a, b);
    }
}
