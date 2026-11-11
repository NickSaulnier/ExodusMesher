//! Rounding the exact mesh to f64 coordinates (Milestone 13).
//!
//! Exodus stores doubles. Every constructed (implicit) point is rounded to
//! the f64 point its construction approximates; input points are f64
//! already. Rounding can turn thin tets over, so every inside tet's
//! orientation is then checked exactly on the rounded coordinates, and those
//! turned over, or with no positive volume in floating point from any vertex
//! order that keeps their orientation, are untangled ([`Untangle`]): a constructed vertex moves
//! within a few ulps, or the tet is flipped away without touching the
//! boundary. Each fix is returned for the report. Afterwards the mesh is
//! checked again: every inside tet strictly positive, the boundary exactly
//! the recovered subfaces, and each input triangle's area covered by its
//! subfaces (the boundary moves by rounding only). What fails is an error,
//! never a mesh with a non-positive tet.

use std::collections::{BTreeMap, BTreeSet};

use mesher_geom::implicit::Point;
use mesher_geom::{Point3, Sign};

pub use super::delaunay::Untangle;
use super::regions::Region;

/// Local faces of a positive tet, outward.
const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

/// The mesh in f64 coordinates.
#[derive(Debug, Clone)]
pub struct RoundedMesh {
    pub points: Vec<Point3>,
    /// Every tet strictly positive on `points`, exactly and in floating
    /// point.
    pub tets: Vec<[u32; 4]>,
    /// Boundary faces, outward, each with its input triangle.
    pub boundary: Vec<([u32; 3], u32)>,
    /// Each point's vertex in the exact tetrahedralization; below the input
    /// count it is the input vertex.
    pub origin: Vec<u32>,
    /// Constructed points rounded (Steiner points of segment and face
    /// recovery).
    pub rounded: usize,
    /// Inside tets that rounding turned over, before untangling.
    pub inverted_before: usize,
    /// Inside tets positive exactly but with no positive volume in floating
    /// point (as a finite-volume solver computes it), before untangling.
    pub flat_before: usize,
    /// The fixes, in order.
    pub untangles: Vec<Untangle>,
}

/// Why rounding failed.
#[derive(Debug, Clone, PartialEq)]
pub enum RoundingError {
    /// A constructed point is not finite in f64.
    NonFinite { vertex: u32 },
    /// The boundary no longer bounds a region (a bug: rounding and flips
    /// keep the subfaces).
    Region,
    /// Tets still turned over or flat in floating point after untangling,
    /// by vertex positions.
    Inverted {
        tets: Vec<[Point3; 4]>,
        inverted_before: usize,
        untangles: usize,
    },
    /// An input triangle's subfaces are no longer its boundary faces, or no
    /// longer cover it.
    Boundary { triangle: u32 },
}

impl std::fmt::Display for RoundingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RoundingError::NonFinite { vertex } => {
                write!(f, "constructed point {vertex} is not finite in f64")
            }
            RoundingError::Region => write!(f, "the boundary no longer bounds a region"),
            RoundingError::Inverted {
                tets,
                inverted_before,
                untangles,
            } => write!(
                f,
                "{} tets stay turned over after rounding ({inverted_before} before, {untangles} fixes)",
                tets.len()
            ),
            RoundingError::Boundary { triangle } => {
                write!(
                    f,
                    "input triangle {triangle} is not its subfaces after rounding"
                )
            }
        }
    }
}

fn key(u: u32, v: u32) -> (u32, u32) {
    (u.min(v), u.max(v))
}

fn area(t: [Point3; 3]) -> f64 {
    let u = [0, 1, 2].map(|k| t[1][k] - t[0][k]);
    let w = [0, 1, 2].map(|k| t[2][k] - t[0][k]);
    let n = [
        u[1] * w[2] - u[2] * w[1],
        u[2] * w[0] - u[0] * w[2],
        u[0] * w[1] - u[1] * w[0],
    ];
    0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
}

/// Round the region of `mesh_region` to f64 and untangle what turned over.
/// `points` and `triangles` are the input it was made from.
///
/// # Errors
/// See [`RoundingError`].
pub fn round_region(
    region: Region,
    points: &[Point3],
    triangles: &[[u32; 3]],
) -> Result<RoundedMesh, RoundingError> {
    let Region { recovered, .. } = region;
    let mut dt = recovered.dt;
    let exact: Vec<Point> = dt.points().to_vec();
    let movable: Vec<bool> = exact.iter().map(|p| !p.is_explicit()).collect();
    let mut rounded = 0;
    for (i, p) in exact.iter().enumerate() {
        if !p.is_explicit() {
            let q = Point::explicit(p.approx())
                .map_err(|_| RoundingError::NonFinite { vertex: i as u32 })?;
            dt.set_point(i as u32, q);
            rounded += 1;
        }
    }
    let walls: BTreeMap<[u32; 3], u32> = recovered
        .subfaces
        .iter()
        .map(|s| {
            let mut v = s.v;
            v.sort_unstable();
            (v, s.triangle)
        })
        .collect();
    // Subsegments and every subface side stay.
    let mut protected: BTreeSet<(u32, u32)> = recovered
        .subsegments
        .iter()
        .map(|s| key(s.ends[0], s.ends[1]))
        .collect();
    for s in &recovered.subfaces {
        for k in 0..3 {
            protected.insert(key(s.v[k], s.v[(k + 1) % 3]));
        }
    }
    let wall = |f: [u32; 3]| walls.contains_key(&f);
    let is_protected = |x: u32, y: u32| protected.contains(&key(x, y));
    let u = dt
        .untangle(&movable, &exact, &wall, &is_protected)
        .ok_or(RoundingError::Region)?;
    let (inverted_before, flat_before, untangles, left) = (u.inverted, u.flat, u.fixes, u.left);
    let at = |v: u32| dt.points()[v as usize].approx();
    if !left.is_empty() {
        return Err(RoundingError::Inverted {
            tets: left.iter().map(|t| t.map(at)).collect(),
            inverted_before,
            untangles: untangles.len(),
        });
    }
    // The inside, as it is now.
    let labelled = dt.flood_regions(&wall).map_err(|_| RoundingError::Region)?;
    // Each in the vertex order (keeping its orientation) whose volume in
    // floating point is largest.
    let inside: Vec<[u32; 4]> = labelled
        .into_iter()
        .filter(|&(_, i)| i)
        .map(|(v, _)| crate::quality::best_float_order(v, &at).0)
        .collect();
    // Every inside tet strictly positive, exactly, on the rounded points.
    let rounded_points: Vec<Point> = dt.points().to_vec();
    for v in &inside {
        let [a, b, c, d] = v.map(|x| &rounded_points[x as usize]);
        if mesher_geom::implicit::orient3d(a, b, c, d) != Sign::Positive
            || crate::quality::float_volume(v.map(at)) <= 0.0
        {
            return Err(RoundingError::Inverted {
                tets: vec![v.map(at)],
                inverted_before,
                untangles: untangles.len(),
            });
        }
    }
    // The boundary is exactly the subfaces, each once.
    let mut boundary = Vec::new();
    let mut seen: BTreeSet<[u32; 3]> = BTreeSet::new();
    for v in &inside {
        for f in FACES {
            let face = f.map(|k| v[k]);
            let mut s = face;
            s.sort_unstable();
            if let Some(&t) = walls.get(&s) {
                if !seen.insert(s) {
                    return Err(RoundingError::Boundary { triangle: t });
                }
                boundary.push((face, t));
            }
        }
    }
    if let Some((_, &t)) = walls.iter().find(|(s, _)| !seen.contains(*s)) {
        return Err(RoundingError::Boundary { triangle: t });
    }
    // Each input triangle covered by its subfaces, to rounding (relative to
    // the bounding box, as needles have next to no area).
    let (lo, hi) = points.iter().fold(
        ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
        |(lo, hi), p| {
            (
                [0, 1, 2].map(|k| lo[k].min(p[k])),
                [0, 1, 2].map(|k| hi[k].max(p[k])),
            )
        },
    );
    let diagonal2: f64 = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum();
    // By facet (merged flat triangles together).
    let facet_of = &recovered.facet_of;
    let facet = |t: usize| facet_of.get(t).copied().unwrap_or(t as u32) as usize;
    let mut covered = vec![0.0f64; triangles.len()];
    let mut whole = vec![0.0f64; triangles.len()];
    for &(f, t) in &boundary {
        if let Some(c) = covered.get_mut(facet(t as usize)) {
            *c += area(f.map(at));
        }
    }
    for (t, tri) in triangles.iter().enumerate() {
        if let Some(w) = whole.get_mut(facet(t)) {
            *w += area(tri.map(|v| points[v as usize]));
        }
    }
    for t in 0..triangles.len() {
        if facet(t) == t && (covered[t] - whole[t]).abs() > 1e-9 * whole[t] + 1e-12 * diagonal2 {
            return Err(RoundingError::Boundary { triangle: t as u32 });
        }
    }
    // Compact to the vertices used.
    let mut index: BTreeMap<u32, u32> = BTreeMap::new();
    for v in &inside {
        for &x in v {
            let n = index.len() as u32;
            index.entry(x).or_insert(n);
        }
    }
    let mut origin = vec![0u32; index.len()];
    for (&x, &i) in &index {
        origin[i as usize] = x;
    }
    let local = |x: u32| index.get(&x).copied().unwrap_or(0);
    Ok(RoundedMesh {
        points: origin.iter().map(|&x| at(x)).collect(),
        tets: inside.iter().map(|v| v.map(local)).collect(),
        boundary: boundary.iter().map(|&(f, t)| (f.map(local), t)).collect(),
        origin,
        rounded,
        inverted_before,
        flat_before,
        untangles,
    })
}
