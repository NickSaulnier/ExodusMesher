//! Volume meshing kernels. Milestone 6: isosurface stuffing (Labelle and
//! Shewchuk, 2007). The constrained Delaunay kernel follows in Milestones 7–8.

pub mod quality;
pub mod stuffing;

pub use quality::{dihedral_angles, dihedral_range};
pub use stuffing::{Stuffed, StuffingError, StuffingParams, stuff};
