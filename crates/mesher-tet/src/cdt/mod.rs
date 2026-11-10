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
//! - [`faces`]: face recovery by flips, then cavity retetrahedralization,
//!   then Steiner points on the faces.
//! - [`regions`]: the tets inside the surface, by a parity flood from the
//!   hull.
//!
//! Rounding to f64 (Milestone 13) and the pipeline integration follow.

pub mod delaunay;
pub mod faces;
pub mod regions;
pub mod segments;

pub use delaunay::{DelaunayError, ExactDelaunay, GHOST, Inserted, Refill};
pub use faces::{FaceError, RecoveredFaces, Subface, recover_faces};
pub use regions::{Region, RegionError, mesh_region};
pub use segments::{RecoveredSegments, SegmentError, recover_segments};
