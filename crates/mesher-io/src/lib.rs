//! Readers and writers: the OBJ reader and the defect overlay (.vtp). The
//! Exodus II, .msh and .vtu writers arrive in Milestone 5.

pub mod obj;
pub mod vtp;

pub use mesher_core::model::{Polygon, SurfaceModel, Triangle};
pub use obj::{ObjModel, parse_obj, read_obj};
