//! Defect geometry for the overlay: which faces, edges and vertices carry
//! which codes. Filled by diagnosis, written as .vtp by `mesher-io` and drawn
//! live by the GUI. Vertex and face indices are model indices.

use crate::report::DefectCode;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Overlay {
    pub faces: Vec<(u32, DefectCode)>,
    pub edges: Vec<([u32; 2], DefectCode)>,
    pub vertices: Vec<(u32, DefectCode)>,
}

impl Overlay {
    pub fn faces(&mut self, faces: impl IntoIterator<Item = u32>, code: DefectCode) {
        self.faces.extend(faces.into_iter().map(|f| (f, code)));
    }

    pub fn edge(&mut self, a: u32, b: u32, code: DefectCode) {
        self.edges.push(([a.min(b), a.max(b)], code));
    }

    pub fn vertices(&mut self, vertices: impl IntoIterator<Item = u32>, code: DefectCode) {
        self.vertices
            .extend(vertices.into_iter().map(|v| (v, code)));
    }
}
