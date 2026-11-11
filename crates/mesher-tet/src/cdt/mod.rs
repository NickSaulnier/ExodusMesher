//! The exact constrained Delaunay kernel (Milestone 12), following Diazzi et
//! al. (2023), built beside the current kernel (`crate::conforming`).
//!
//! Every point is a `mesher_geom::implicit::Point` and every decision an exact
//! implicit predicate: input vertices are explicit, and Steiner points are
//! constructed (a point on an input segment is an LNC of the segment's
//! endpoints), so nothing is rounded until the mesh is written. There is no
//! facet merging by tolerance: every input triangle is a constraint.
//!
//! - [`delaunay`]: the exact incremental Delaunay tetrahedralization, with
//!   flips, edge removal by retriangulating an edge's ring, and insertion
//!   that protects chosen edges.
//! - [`segments`]: segment recovery by flips, then splitting with LNC
//!   Steiner points.
//! - [`faces`]: face recovery by flips, then cavity retetrahedralization
//!   (gift-wrapping where expansion fails), then Steiner points on the
//!   faces.
//! - [`regions`]: the tets inside the surface, by a parity flood from the
//!   hull.
//! - `locate`: for the enriched CDT (Milestone 14), extra points (implicit or
//!   not) placed exactly on the input segments and triangles they lie on;
//!   with `merge_flat`, coplanar triangles across flat edges are one facet
//!   (`recover_faces_with`), as is a triangle with points inside or a long
//!   polygon, recovered by the same cavity and Steiner steps
//!   and read off by a flood from its boundary.
//! - [`rounding`]: the mesh in f64 coordinates, with what rounding turned
//!   over untangled (Milestone 13).

pub mod delaunay;
pub mod faces;
mod locate;
pub mod regions;
pub mod rounding;
pub mod segments;

pub use delaunay::{DelaunayError, ExactDelaunay, GHOST, Inserted, Refill};
pub use faces::{FaceError, RecoveredFaces, Subface, recover_faces, recover_faces_with};
pub use regions::{Region, RegionError, mesh_region, mesh_region_with};
pub use rounding::{RoundedMesh, RoundingError, Untangle, round_region};
pub use segments::{RecoveredSegments, SegmentError, recover_segments, recover_segments_with};
