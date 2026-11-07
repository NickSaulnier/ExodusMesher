//! Exact predicates and geometry kernels.
//!
//! Every orientation decision in ExodusMesher goes through the adaptive-precision
//! predicates in [`predicates`]; nothing here compares floats against an epsilon to
//! decide topology.

pub mod polygon;
pub mod predicates;

pub use polygon::{FanReason, MAX_EAR_CLIP_VERTICES, Triangulation, TriangulationMethod};
pub use predicates::Sign;
