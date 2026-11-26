//! Readers and writers: the OBJ reader and writer, the defect overlay (.vtp),
//! and the volume-mesh writers: Exodus II, Gmsh .msh 4.1, VTK .vtu, an
//! OpenFOAM polyMesh case and a Fluent mesh file.

pub mod exodus;
pub mod fluent;
pub mod foam;
pub mod msh;
pub mod obj;
mod shared;
pub mod vtp;
pub mod vtu;

pub use mesher_core::model::{Polygon, SurfaceModel, Triangle};
pub use obj::{NamedSurface, ObjModel, parse_obj, read_obj, write_obj};
