//! Readers and writers: the OBJ reader and writer, the defect overlay (.vtp),
//! and the volume-mesh writers: Exodus II, Gmsh .msh 4.1 and VTK .vtu.

pub mod exodus;
pub mod msh;
pub mod obj;
pub mod vtp;
pub mod vtu;

pub use mesher_core::model::{Polygon, SurfaceModel, Triangle};
pub use obj::{NamedSurface, ObjModel, parse_obj, read_obj, write_obj};
