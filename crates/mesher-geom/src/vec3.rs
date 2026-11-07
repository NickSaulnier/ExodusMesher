//! Small floating-point vector helpers on `[f64; 3]`. These are for measurements
//! (distances, solid angles, ray parameters), never for topological decisions.

use crate::Point3;

#[must_use]
pub fn add(a: Point3, b: Point3) -> Point3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[must_use]
pub fn sub(a: Point3, b: Point3) -> Point3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[must_use]
pub fn scale(a: Point3, s: f64) -> Point3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[must_use]
pub fn dot(a: Point3, b: Point3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[must_use]
pub fn cross(a: Point3, b: Point3) -> Point3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[must_use]
pub fn norm(a: Point3) -> f64 {
    dot(a, a).sqrt()
}

#[must_use]
pub fn distance_squared(a: Point3, b: Point3) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}
