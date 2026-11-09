//! The exact constrained Delaunay kernel (Milestone 12), following Diazzi et
//! al. (2023), built beside the current kernel (`crate::conforming`).
//!
//! Every point is a `mesher_geom::implicit::Point` and every decision an exact
//! implicit predicate: input vertices are explicit, and Steiner points are
//! constructed (a point on an input segment is an LNC of the segment's
//! endpoints), so nothing is rounded until the mesh is written. There is no
//! facet merging by tolerance: every input triangle is a constraint.
//!
//! - [`delaunay`]: the exact incremental Delaunay tetrahedralization.
//! - [`segments`]: segment recovery by splitting with LNC Steiner points.
//!
//! Face recovery, region marking and the pipeline integration follow.

pub mod delaunay;
pub mod segments;

pub use delaunay::{DelaunayError, ExactDelaunay, GHOST, Inserted};
pub use segments::{RecoveredSegments, SegmentError, recover_segments};
