//! Readers and writers. Milestone 1 has the OBJ reader; the Exodus II, .msh,
//! .vtu and .vtp writers arrive in later milestones.

pub mod obj;

pub use obj::{ObjModel, Polygon, Triangle, parse_obj, read_obj};
