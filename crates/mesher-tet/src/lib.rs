//! Volume meshing kernels: isosurface stuffing (Labelle and Shewchuk, 2007)
//! and incremental Delaunay tetrahedralization (Bowyer-Watson), the base of
//! the constrained Delaunay kernel (Milestone 8).

pub mod delaunay;
pub mod quality;
pub mod stuffing;

pub use delaunay::{Delaunay, DelaunayError, Inserted, Tetrahedralization, delaunay};
pub use quality::{dihedral_angles, dihedral_range};
pub use stuffing::{Stuffed, StuffingError, StuffingParams, stuff};
