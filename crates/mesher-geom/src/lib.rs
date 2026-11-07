//! Exact predicates and geometry kernels.
//!
//! Every orientation decision in ExodusMesher goes through the adaptive-precision
//! predicates in [`predicates`]; nothing here compares floats against an epsilon to
//! decide topology. Floating-point arithmetic ([`vec3`]) is used only for
//! measurements: distances, ray parameters, solid angles.

pub mod bvh;
pub mod halfedge;
pub mod intersect;
pub mod polygon;
pub mod predicates;
pub mod soup;
pub mod vec3;
pub mod winding;

/// A point or vector in 3D.
pub type Point3 = [f64; 3];

pub use bvh::{Aabb, Bvh, Nearest, Ray, RayHit};
pub use halfedge::{EdgeId, EdgeKind, HalfEdge, SurfaceMesh};
pub use intersect::{
    faces_intersect, segment_intersects_triangle, segments_intersect, triangles_intersect,
};
pub use polygon::{FanReason, MAX_EAR_CLIP_VERTICES, Triangulation, TriangulationMethod};
pub use predicates::Sign;
pub use soup::{GeomError, TriangleSoup};
pub use winding::{FastWindingNumber, solid_angle, winding_number_brute};
