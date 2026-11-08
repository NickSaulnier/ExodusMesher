//! Volume meshing kernels: isosurface stuffing (Labelle and Shewchuk, 2007)
//! and incremental Delaunay tetrahedralization (Bowyer-Watson), the base of
//! the constrained Delaunay kernel (Milestone 8).

pub mod conforming;
pub mod delaunay;
pub mod quality;
pub mod stuffing;

pub use conforming::{ConformingError, ConformingMesh, ConformingParams, conforming_mesh};
pub use delaunay::{Delaunay, DelaunayError, Inserted, Tetrahedralization, delaunay};
pub use quality::{FvQuality, dihedral_angles, dihedral_range, fv_quality};
pub use stuffing::{Stuffed, StuffingError, StuffingParams, stuff};
