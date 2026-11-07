//! VTK XML UnstructuredGrid (.vtu, ASCII) writer: points, tetrahedra
//! (VTK type 10, same node order as ours) and a `block` cell array holding
//! each tet's 1-based block number.

use std::fmt::Write as _;

use mesher_core::mesh::VolumeMesh;

const VTK_TETRA: u8 = 10;

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
