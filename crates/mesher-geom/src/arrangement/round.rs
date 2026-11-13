//! Stage 4 of the arrangement (Milestone 16): the resolved surface in f64.
//!
//! Rounding the implicit points can break what the exact surface had, so
//! each round is repaired and checked, and the whole repeated:
//!
//! 1. every point to its f64 position ([`Point::approx`], within ulps);
//!    points that land on one position merge, triangles that collapse
//!    go, and two copies of one triangle facing opposite ways cancel;
//! 2. a triangle left flat (three corners exactly collinear) goes, and its
//!    neighbour across the longest edge is split at the middle corner;
//! 3. the result checked for triangles meeting beyond their shared
//!    vertices and edges ([`super::intersections`], exact on f64); if any
//!    do, the rounded surface is resolved again (stages 1 to 3, now on
//!    explicit points) and rounded again, as autorefine-and-snap does.
//!
//! Every step keeps each edge used as often in each direction, so the
//! surface stays closed and oriented. Each triangle keeps the input
//! triangle it lies on, through every round.

use std::collections::BTreeMap;

use super::{ArrangementError, Resolved, intersections, resolve};
use crate::Point3;
use crate::implicit::{self, Point};
use crate::predicates::{Sign, plane_axis};

/// Rounds of resolving again before giving up.
const MAX_ROUNDS: usize = 6;
/// Flat triangles one round may repair.
const MAX_FLAT_FIXES: usize = 100_000;

/// The resolved surface in f64: closed, oriented, with no zero-area
/// triangle and no two triangles meeting beyond their shared vertices and
/// edges.
#[derive(Debug, Clone)]
pub struct Rounded {
    pub points: Vec<Point3>,
    pub triangles: Vec<[u32; 3]>,
    /// The input triangle each lies on.
    pub origin: Vec<u32>,
    /// The first resolution's counts (see [`Resolved`]).
    pub pairs: usize,
    pub coplanar_pairs: usize,
    pub faces: usize,
    pub patches: usize,
    pub dropped: usize,
    /// Times the rounded surface had to be resolved again.
    pub rounds: usize,
    /// Over all rounds: points merged by rounding, triangles that
    /// collapsed or cancelled, and flat triangles repaired.
    pub merged_points: usize,
    pub removed_triangles: usize,
    pub flat_repaired: usize,
}

/// [`resolve`], then the result in f64.
///
/// # Errors
/// As [`resolve`] (also for a rounded surface resolved again), and
/// [`ArrangementError::Rounding`] if rounding still breaks the surface
/// after the last round.
pub fn resolve_rounded(
    points: &[Point3],
    triangles: &[[u32; 3]],
) -> Result<Rounded, ArrangementError> {
    let first = resolve(points, triangles)?;
    let mut out = Rounded {
        points: Vec::new(),
        triangles: Vec::new(),
        origin: Vec::new(),
        pairs: first.pairs,
        coplanar_pairs: first.coplanar_pairs,
        faces: first.faces,
        patches: first.patches,
        dropped: first.dropped,
        rounds: 0,
        merged_points: 0,
        removed_triangles: 0,
        flat_repaired: 0,
    };
    let mut current = first;
    // Origins of the current resolution's triangles in the input.
    let mut to_input: Option<Vec<u32>> = None;
    loop {
        let mesh = round(&current, &mut out)?;
        let origin: Vec<u32> = match &to_input {
            None => mesh.origin.clone(),
            Some(map) => mesh.origin.iter().map(|&o| map[o as usize]).collect(),
        };
        let crossing = intersections(&mesh.points, &mesh.triangles)?.pairs;
        if crossing == 0 {
            out.points = mesh.points;
            out.triangles = mesh.triangles;
            out.origin = origin;
            return Ok(out);
        }
        if out.rounds == MAX_ROUNDS {
            return Err(ArrangementError::Rounding { pairs: crossing });
        }
        out.rounds += 1;
        current = resolve(&mesh.points, &mesh.triangles)?;
        to_input = Some(origin);
    }
}

/// One rounded mesh: triangles over f64 points, with the triangle of the
/// resolution each came from.
struct Mesh {
    points: Vec<Point3>,
    triangles: Vec<[u32; 3]>,
    origin: Vec<u32>,
}

fn round(r: &Resolved, stats: &mut Rounded) -> Result<Mesh, ArrangementError> {
    // 1. Positions, merged where they coincide.
    let mut index: BTreeMap<[u64; 3], u32> = BTreeMap::new();
    let mut points: Vec<Point3> = Vec::new();
    let mut map: Vec<u32> = vec![u32::MAX; r.points.len()];
    let mut used = vec![false; r.points.len()];
    for t in &r.triangles {
        for &v in t {
            used[v as usize] = true;
        }
    }
    for (i, p) in r.points.iter().enumerate() {
        if !used[i] {
            continue;
        }
        let q = p.approx();
        if !q.iter().all(|x| x.is_finite()) {
            return Err(ArrangementError::Triangulation(
                "a point has no finite position".into(),
            ));
        }
        // +0 and -0 are one position.
        let key = q.map(|x| (x + 0.0).to_bits());
        let next = points.len() as u32;
        let id = *index.entry(key).or_insert(next);
        if id == next {
            points.push(q);
        } else {
            stats.merged_points += 1;
        }
        map[i] = id;
    }
    let mut tris: Vec<([u32; 3], u32)> = Vec::new();
    for (t, &o) in r.triangles.iter().zip(&r.origin) {
        let t = t.map(|v| map[v as usize]);
        if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
            stats.removed_triangles += 1;
            continue;
        }
        tris.push((t, o));
    }
    // Copies facing opposite ways cancel.
    let mut by_key: BTreeMap<[u32; 3], Vec<usize>> = BTreeMap::new();
    for (i, (t, _)) in tris.iter().enumerate() {
        let mut k = *t;
        k.sort_unstable();
        by_key.entry(k).or_default().push(i);
    }
    let mut gone = vec![false; tris.len()];
    for copies in by_key.values() {
        let (mut ups, mut downs): (Vec<usize>, Vec<usize>) = (Vec::new(), Vec::new());
        let first = tris[copies[0]].0;
        for &i in copies {
            if same_turn(first, tris[i].0) {
                ups.push(i);
            } else {
                downs.push(i);
            }
        }
        for (&u, &d) in ups.iter().zip(&downs) {
            gone[u] = true;
            gone[d] = true;
            stats.removed_triangles += 2;
        }
    }
    let mut triangles = Vec::new();
    let mut origin = Vec::new();
    for (i, (t, o)) in tris.into_iter().enumerate() {
        if !gone[i] {
            triangles.push(t);
            origin.push(o);
        }
    }
    // 2. Flat triangles.
    repair_flat(&points, &mut triangles, &mut origin, stats)?;
    Ok(Mesh {
        points,
        triangles,
        origin,
    })
}

fn same_turn(a: [u32; 3], b: [u32; 3]) -> bool {
    let k = b.iter().position(|&v| v == a[0]).unwrap_or(0);
    b[(k + 1) % 3] == a[1]
}

/// Remove each flat triangle `a b c` (`c` exactly between `a` and `b`) and
/// split its neighbour `b a d` across `ab` into `b c d` and `c a d`.
fn repair_flat(
    points: &[Point3],
    triangles: &mut Vec<[u32; 3]>,
    origin: &mut Vec<u32>,
    stats: &mut Rounded,
) -> Result<(), ArrangementError> {
    let flat = |t: [u32; 3]| {
        let [a, b, c] = t.map(|v| points[v as usize]);
        plane_axis(a, b, c).is_none()
    };
    let mut queue: Vec<usize> = (0..triangles.len())
        .filter(|&i| flat(triangles[i]))
        .collect();
    if queue.is_empty() {
        return Ok(());
    }
    let mut live = vec![true; triangles.len()];
    let mut holding: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for (i, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            holding.insert((t[k], t[(k + 1) % 3]), i);
        }
    }
    let mut fixes = 0;
    while let Some(i) = queue.pop() {
        if !live[i] || !flat(triangles[i]) {
            continue;
        }
        fixes += 1;
        if fixes > MAX_FLAT_FIXES {
            return Err(ArrangementError::Rounding { pairs: 0 });
        }
        let t = triangles[i];
        // The corner between the other two: the one where they lie on
        // opposite sides (the angle there is straight).
        let p = t.map(|v| Point::from(points[v as usize]));
        let Some(k) = (0..3).find(|&k| {
            implicit::dot_sign(&p[(k + 1) % 3], &p[k], &p[(k + 2) % 3]) == Sign::Negative
        }) else {
            return Err(ArrangementError::Triangulation(
                "a flat triangle with no corner between the others".into(),
            ));
        };
        let (c, a, b) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
        // The neighbour holding b → a.
        let Some(&n) = holding.get(&(b, a)) else {
            return Err(ArrangementError::Triangulation(
                "a flat triangle with no neighbour across its long edge".into(),
            ));
        };
        let nt = triangles[n];
        let Some(&d) = nt.iter().find(|&&v| v != a && v != b) else {
            return Err(ArrangementError::Triangulation("a broken neighbour".into()));
        };
        live[i] = false;
        for x in [i, n] {
            let tx = triangles[x];
            for e in 0..3 {
                let key = (tx[e], tx[(e + 1) % 3]);
                if holding.get(&key) == Some(&x) {
                    holding.remove(&key);
                }
            }
        }
        // b c d replaces the neighbour; c a d is new.
        triangles[n] = [b, c, d];
        let m = triangles.len();
        triangles.push([c, a, d]);
        origin.push(origin[n]);
        live.push(true);
        for (x, tx) in [(n, triangles[n]), (m, triangles[m])] {
            for e in 0..3 {
                holding.insert((tx[e], tx[(e + 1) % 3]), x);
            }
        }
        queue.push(n);
        queue.push(m);
        stats.flat_repaired += 1;
    }
    let mut kept_t = Vec::new();
    let mut kept_o = Vec::new();
    for (i, (t, o)) in triangles.iter().zip(origin.iter()).enumerate() {
        if live[i] {
            kept_t.push(*t);
            kept_o.push(*o);
        }
    }
    *triangles = kept_t;
    *origin = kept_o;
    Ok(())
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;

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

    /// A tetrahedron whose face over the edge a b is split at its midpoint
    /// c, with the flat triangle a b c closing the surface: the repair
    /// splits the other face at c too.
    #[test]
    fn a_flat_triangle_is_repaired() {
        let points: Vec<Point3> = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [1.0, 0.0, 0.0],
        ];
        // a = 0, b = 1, c = 4; faces of the tetrahedron 0 1 2 3, outward.
        let mut triangles = vec![
            [0, 4, 3],
            [4, 1, 3],
            [0, 1, 4],
            [1, 0, 2],
            [0, 3, 2],
            [1, 2, 3],
        ];
        assert!(balanced(&triangles));
        let mut origin = vec![0, 0, 0, 1, 2, 3];
        let mut stats = Rounded {
            points: Vec::new(),
            triangles: Vec::new(),
            origin: Vec::new(),
            pairs: 0,
            coplanar_pairs: 0,
            faces: 0,
            patches: 0,
            dropped: 0,
            rounds: 0,
            merged_points: 0,
            removed_triangles: 0,
            flat_repaired: 0,
        };
        repair_flat(&points, &mut triangles, &mut origin, &mut stats).unwrap();
        assert_eq!(stats.flat_repaired, 1);
        assert_eq!(triangles.len(), 6);
        assert!(balanced(&triangles));
        for t in &triangles {
            let [a, b, c] = t.map(|v| points[v as usize]);
            assert!(plane_axis(a, b, c).is_some(), "{t:?}");
        }
        // The split face keeps its origin.
        assert_eq!(origin.iter().filter(|&&o| o == 1).count(), 2);
    }
}
