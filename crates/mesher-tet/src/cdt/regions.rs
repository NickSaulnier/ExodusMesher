//! Region marking: which tets of the recovered tetrahedralization are inside
//! the surface.
//!
//! A flood from the hull ([`ExactDelaunay::flood_regions`]) starts outside
//! and toggles at every recovered subface it crosses, as the current kernel
//! does. It is topological, so a tet thinner than rounding is classified as
//! reliably as any other, and nested shells (holes) come out right. Labels
//! that disagree mean the subfaces do not bound a region, which face
//! recovery rules out, so that is reported rather than guessed around.
//!
//! The points are still exact (implicit): rounding them to f64 is
//! Milestone 13.

use std::collections::BTreeMap;

use mesher_geom::Point3;

use super::delaunay::ExactDelaunay;
use super::faces::{FaceError, RecoveredFaces, recover_faces};

/// Local faces of a positive tet, outward.
const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

/// Why the region could not be marked.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RegionError {
    Faces(FaceError),
    /// Input triangles face recovery left crossed by an edge.
    Missing {
        triangles: Vec<u32>,
    },
    /// The flood's labels disagree across this face (vertices sorted).
    Parity {
        face: [u32; 3],
    },
}

impl std::fmt::Display for RegionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RegionError::Faces(e) => write!(f, "{e}"),
            RegionError::Missing { triangles } => write!(
                f,
                "{} input triangles are not recovered (first {})",
                triangles.len(),
                triangles.first().copied().unwrap_or(0)
            ),
            RegionError::Parity { face } => write!(
                f,
                "inside and outside meet across face {}-{}-{}",
                face[0], face[1], face[2]
            ),
        }
    }
}

impl From<FaceError> for RegionError {
    fn from(e: FaceError) -> RegionError {
        RegionError::Faces(e)
    }
}

/// The tets inside the surface, with exact points.
#[derive(Debug, Clone)]
pub struct Region {
    /// Segment and face recovery, whose tetrahedralization holds the points.
    pub recovered: RecoveredFaces,
    /// The tets inside, positively oriented.
    pub tets: Vec<[u32; 4]>,
    /// The boundary faces, outward from the inside, each with its input
    /// triangle.
    pub boundary: Vec<([u32; 3], u32)>,
}

impl Region {
    #[must_use]
    pub fn dt(&self) -> &ExactDelaunay {
        &self.recovered.dt
    }
}

/// Recover every segment and input triangle (see [`recover_faces`]), then
/// keep the tets inside the surface.
///
/// # Errors
/// See [`RegionError`].
pub fn mesh_region(
    points: &[Point3],
    triangles: &[[u32; 3]],
    max_points: usize,
    seed: u64,
) -> Result<Region, RegionError> {
    let recovered = recover_faces(points, triangles, max_points, seed)?;
    if !recovered.missing.is_empty() {
        return Err(RegionError::Missing {
            triangles: recovered.missing.clone(),
        });
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
    let labelled = recovered
        .dt
        .flood_regions(&|f| walls.contains_key(&f))
        .map_err(|face| RegionError::Parity { face })?;
    let tets: Vec<[u32; 4]> = labelled
        .into_iter()
        .filter(|&(_, inside)| inside)
        .map(|(v, _)| v)
        .collect();
    // Every wall face of an inside tet has the outside beyond it.
    let mut boundary = Vec::new();
    for v in &tets {
        for f in FACES {
            let face = f.map(|k| v[k]);
            let mut key = face;
            key.sort_unstable();
            if let Some(&t) = walls.get(&key) {
                boundary.push((face, t));
            }
        }
    }
    Ok(Region {
        recovered,
        tets,
        boundary,
    })
}
