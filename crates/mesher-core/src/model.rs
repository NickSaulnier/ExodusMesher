//! The parsed surface: triangles with provenance back to the input file. The
//! OBJ reader produces it; diagnosis and later stages consume it.

use crate::report::InputStats;

/// A source face as written in the file, with vertex indices into
/// [`SurfaceModel::vertices`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polygon {
    pub vertices: Vec<u32>,
    /// 1-based input line (the first line of a continued statement).
    pub line: u32,
    /// Index into [`SurfaceModel::objects`] of the active `o` statement.
    pub object: Option<u32>,
    /// Index into [`SurfaceModel::groups`] of the active `g` statement.
    pub group: Option<u32>,
}

/// A triangle from triangulating a polygon; winding follows the polygon.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Triangle {
    pub vertices: [u32; 3],
    /// Index into [`SurfaceModel::polygons`].
    pub polygon: u32,
}

/// The parsed geometry plus provenance. Coordinates are finite and every index
/// is in range (the reader guarantees both).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SurfaceModel {
    /// Well-formed vertices, in file order.
    pub vertices: Vec<[f64; 3]>,
    /// 1-based input line of each vertex.
    pub vertex_lines: Vec<u32>,
    pub polygons: Vec<Polygon>,
    pub triangles: Vec<Triangle>,
    /// `o` names in order of first appearance.
    pub objects: Vec<String>,
    /// `g` names in order of first appearance.
    pub groups: Vec<String>,
    pub bytes: u64,
    pub lines: u64,
}

impl SurfaceModel {
    #[must_use]
    pub fn stats(&self) -> InputStats {
        InputStats {
            bytes: self.bytes,
            lines: self.lines,
            vertices: self.vertices.len() as u64,
            polygons: self.polygons.len() as u64,
            triangles: self.triangles.len() as u64,
            objects: self.objects.clone(),
            groups: self.groups.clone(),
        }
    }

    /// The polygon a triangle came from.
    #[must_use]
    pub fn polygon_of(&self, triangle: u32) -> Option<&Polygon> {
        let t = self.triangles.get(triangle as usize)?;
        self.polygons.get(t.polygon as usize)
    }

    /// Input line of a triangle's polygon (0 if unknown).
    #[must_use]
    pub fn triangle_line(&self, triangle: u32) -> u32 {
        self.polygon_of(triangle).map_or(0, |p| p.line)
    }

    /// Input line of a vertex (0 if unknown).
    #[must_use]
    pub fn vertex_line(&self, vertex: u32) -> u32 {
        self.vertex_lines.get(vertex as usize).copied().unwrap_or(0)
    }
}
