//! Volume meshing kernels: isosurface stuffing (Labelle and Shewchuk, 2007)
//! and incremental Delaunay tetrahedralization (Bowyer-Watson), the base of
//! the constrained Delaunay kernel (Milestone 8).

pub mod cdt;
pub mod chamfer;
pub mod conforming;
pub mod delaunay;
mod fxhash;
mod graded;
pub mod prof;
pub mod quality;
pub mod refine;
pub mod remesh;
pub mod sizing;
pub mod smooth;
pub mod stuffing;

pub use conforming::{ConformingError, ConformingMesh, ConformingParams, conforming_mesh};
pub use delaunay::{Delaunay, DelaunayError, Inserted, Tetrahedralization, delaunay};
pub use quality::{FvQuality, dihedral_angles, dihedral_range, fv_quality};
pub use remesh::{RemeshError, RemeshParams, Remeshed, remesh};
pub use sizing::{SizeField, SizingParams};
pub use smooth::{SmoothParams, SmoothStats, smooth_slivers};
pub use stuffing::{Stuffed, StuffingError, StuffingParams, stuff, stuff_graded};
