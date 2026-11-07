//! A validated, borrowed triangle soup: the input to every surface algorithm.

use thiserror::Error;

use crate::Point3;

/// Why a soup cannot be built. Callers map these to report records.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum GeomError {
    #[error("vertex {vertex} has a non-finite coordinate")]
    NonFiniteCoordinate { vertex: u32 },
    #[error("triangle {triangle} references vertex {vertex}, past the last vertex")]
    VertexOutOfRange { triangle: u32, vertex: u32 },
    #[error("{vertices} vertices or {triangles} triangles exceed 32-bit half-edge indexing")]
    TooLarge { vertices: usize, triangles: usize },
}

/// Points and triangles, checked once so algorithms can index without failing:
/// every coordinate is finite, every index is in range, and `3 × triangles`
/// fits in `u32`. Degenerate and repeated triangles are allowed.
#[derive(Debug, Clone, Copy)]
pub struct TriangleSoup<'a> {
    points: &'a [Point3],
    triangles: &'a [[u32; 3]],
}

impl<'a> TriangleSoup<'a> {
    pub fn new(points: &'a [Point3], triangles: &'a [[u32; 3]]) -> Result<Self, GeomError> {
        let too_large = GeomError::TooLarge {
            vertices: points.len(),
            triangles: triangles.len(),
        };
        let n = u32::try_from(points.len()).map_err(|_| too_large.clone())?;
        if triangles.len() > (u32::MAX / 3) as usize {
            return Err(too_large);
        }
        if let Some(v) = points.iter().position(|p| !p.iter().all(|c| c.is_finite())) {
            return Err(GeomError::NonFiniteCoordinate { vertex: v as u32 });
        }
        for (t, tri) in triangles.iter().enumerate() {
            if let Some(&vertex) = tri.iter().find(|&&v| v >= n) {
                return Err(GeomError::VertexOutOfRange {
                    triangle: t as u32,
                    vertex,
                });
            }
        }
        Ok(TriangleSoup { points, triangles })
    }

    #[must_use]
    pub fn points(&self) -> &'a [Point3] {
        self.points
    }

    #[must_use]
    pub fn triangles(&self) -> &'a [[u32; 3]] {
        self.triangles
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.triangles.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    /// Corner positions of triangle `t`, or `None` if `t` is out of range.
    #[must_use]
    pub fn triangle(&self, t: usize) -> Option<[Point3; 3]> {
        self.triangles.get(t).map(|&tri| self.corners(tri))
    }

    /// Corner positions of a triangle taken from this soup (indices are valid).
    pub(crate) fn corners(&self, tri: [u32; 3]) -> [Point3; 3] {
        tri.map(|v| self.points[v as usize])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validation() {
        let pts = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(TriangleSoup::new(&pts, &[[0, 1, 2], [2, 2, 2]]).is_ok());
        assert_eq!(
            TriangleSoup::new(&pts, &[[0, 1, 3]]).unwrap_err(),
            GeomError::VertexOutOfRange {
                triangle: 0,
                vertex: 3
            }
        );
        let bad = [[0.0, f64::INFINITY, 0.0]];
        assert_eq!(
            TriangleSoup::new(&bad, &[]).unwrap_err(),
            GeomError::NonFiniteCoordinate { vertex: 0 }
        );
    }
}
