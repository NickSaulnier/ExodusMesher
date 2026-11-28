//! VTK XML UnstructuredGrid (.vtu, ASCII) writer: points, tetrahedra
//! (VTK type 10, same node order as ours) and a `block` cell array holding
//! each tet's 1-based block number. A polyhedral mesh ([`vtu_poly_text`])
//! has `VTK_POLYHEDRON` cells (type 42): each cell's points in
//! `connectivity`, and its faces, out of the cell, in `faces` and
//! `faceoffsets`.

use std::fmt::Write as _;

use mesher_core::mesh::{PolyMesh, VolumeMesh};

const VTK_TETRA: u8 = 10;
const VTK_POLYHEDRON: u8 = 42;

/// Render `mesh` as .vtu text. Deterministic.
#[must_use]
pub fn vtu_text(mesh: &VolumeMesh) -> String {
    let n_cells = mesh.element_count();
    let mut s = String::new();
    let _ = writeln!(s, r#"<?xml version="1.0"?>"#);
    let _ = writeln!(
        s,
        r#"<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">"#
    );
    let _ = writeln!(s, "  <UnstructuredGrid>");
    let _ = writeln!(
        s,
        r#"    <Piece NumberOfPoints="{}" NumberOfCells="{n_cells}">"#,
        mesh.points.len()
    );

    let _ = writeln!(s, "      <Points>");
    let _ = write!(
        s,
        r#"        <DataArray type="Float64" NumberOfComponents="3" format="ascii">"#
    );
    for p in &mesh.points {
        let _ = write!(s, " {} {} {}", p[0], p[1], p[2]);
    }
    let _ = writeln!(s, " </DataArray>");
    let _ = writeln!(s, "      </Points>");

    let _ = writeln!(s, "      <Cells>");
    let _ = write!(
        s,
        r#"        <DataArray type="Int64" Name="connectivity" format="ascii">"#
    );
    for t in mesh.elements() {
        let _ = write!(s, " {} {} {} {}", t[0], t[1], t[2], t[3]);
    }
    let _ = writeln!(s, " </DataArray>");
    let _ = write!(
        s,
        r#"        <DataArray type="Int64" Name="offsets" format="ascii">"#
    );
    for k in 1..=n_cells {
        let _ = write!(s, " {}", 4 * k);
    }
    let _ = writeln!(s, " </DataArray>");
    let _ = write!(
        s,
        r#"        <DataArray type="UInt8" Name="types" format="ascii">"#
    );
    for _ in 0..n_cells {
        let _ = write!(s, " {VTK_TETRA}");
    }
    let _ = writeln!(s, " </DataArray>");
    let _ = writeln!(s, "      </Cells>");

    let _ = writeln!(s, r#"      <CellData Scalars="block">"#);
    let _ = write!(
        s,
        r#"        <DataArray type="Int32" Name="block" format="ascii">"#
    );
    for (i, block) in mesh.blocks.iter().enumerate() {
        for _ in &block.tets {
            let _ = write!(s, " {}", i + 1);
        }
    }
    let _ = writeln!(s, " </DataArray>");
    let _ = writeln!(s, "      </CellData>");
    let _ = writeln!(s, "    </Piece>");
    let _ = writeln!(s, "  </UnstructuredGrid>");
    let _ = writeln!(s, "</VTKFile>");
    s
}

/// Render a polyhedral mesh (`cells` `dual`) as .vtu text: every cell a
/// `VTK_POLYHEDRON`, its points (each once, in the order its faces meet
/// them) in `connectivity`, its faces in `faces` (per cell: the face count,
/// then each face's point count and points, ordered out of the cell) with
/// `faceoffsets` (where each cell's list ends), and its 1-based block in
/// `block`. Deterministic.
#[must_use]
pub fn vtu_poly_text(mesh: &PolyMesh) -> String {
    let f = &mesh.faces;
    let cell_faces = f.cell_faces();
    let n_cells = f.cells;
    let mut connectivity: Vec<u32> = Vec::new();
    let mut offsets: Vec<usize> = Vec::with_capacity(n_cells);
    let mut faces: Vec<u32> = Vec::new();
    let mut face_offsets: Vec<usize> = Vec::with_capacity(n_cells);
    let mut seen: Vec<u32> = Vec::new();
    for c in 0..n_cells {
        let list = cell_faces.get(c);
        faces.push(list.len() as u32);
        seen.clear();
        for &k in list {
            let nodes = f.faces.get(k as usize);
            faces.push(nodes.len() as u32);
            let outward = f.owner.get(k as usize) == Some(&(c as u32));
            let ordered: Vec<u32> = if outward {
                nodes.to_vec()
            } else {
                nodes.iter().rev().copied().collect()
            };
            for &v in &ordered {
                if !seen.contains(&v) {
                    seen.push(v);
                }
            }
            faces.extend(ordered);
        }
        connectivity.extend_from_slice(&seen);
        offsets.push(connectivity.len());
        face_offsets.push(faces.len());
    }

    let mut s = String::new();
    let _ = writeln!(s, r#"<?xml version="1.0"?>"#);
    let _ = writeln!(
        s,
        r#"<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">"#
    );
    let _ = writeln!(s, "  <UnstructuredGrid>");
    let _ = writeln!(
        s,
        r#"    <Piece NumberOfPoints="{}" NumberOfCells="{n_cells}">"#,
        mesh.points.len()
    );
    let _ = writeln!(s, "      <Points>");
    let _ = write!(
        s,
        r#"        <DataArray type="Float64" NumberOfComponents="3" format="ascii">"#
    );
    for p in &mesh.points {
        let _ = write!(s, " {} {} {}", p[0], p[1], p[2]);
    }
    let _ = writeln!(s, " </DataArray>");
    let _ = writeln!(s, "      </Points>");
    let _ = writeln!(s, "      <Cells>");
    let array =
        |s: &mut String, kind: &str, name: &str, values: &mut dyn Iterator<Item = String>| {
            let _ = write!(
                s,
                r#"        <DataArray type="{kind}" Name="{name}" format="ascii">"#
            );
            for v in values {
                let _ = write!(s, " {v}");
            }
            let _ = writeln!(s, " </DataArray>");
        };
    array(
        &mut s,
        "Int64",
        "connectivity",
        &mut connectivity.iter().map(u32::to_string),
    );
    array(
        &mut s,
        "Int64",
        "offsets",
        &mut offsets.iter().map(usize::to_string),
    );
    array(
        &mut s,
        "UInt8",
        "types",
        &mut (0..n_cells).map(|_| VTK_POLYHEDRON.to_string()),
    );
    array(
        &mut s,
        "Int64",
        "faces",
        &mut faces.iter().map(u32::to_string),
    );
    array(
        &mut s,
        "Int64",
        "faceoffsets",
        &mut face_offsets.iter().map(usize::to_string),
    );
    let _ = writeln!(s, "      </Cells>");
    let _ = writeln!(s, r#"      <CellData Scalars="block">"#);
    let mut blocks = f
        .zones
        .iter()
        .enumerate()
        .flat_map(|(i, z)| std::iter::repeat_n((i + 1).to_string(), z.count));
    array(&mut s, "Int32", "block", &mut blocks);
    let _ = writeln!(s, "      </CellData>");
    let _ = writeln!(s, "    </Piece>");
    let _ = writeln!(s, "  </UnstructuredGrid>");
    let _ = writeln!(s, "</VTKFile>");
    s
}
