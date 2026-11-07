//! Defect overlay as VTK XML PolyData (.vtp), for ParaView.
//!
//! Cells are, in VTK's order: one vertex cell per vertex defect, one line per
//! edge defect, then every model triangle. Cell arrays:
//! - `defect`: index into the code list in the file's header comment (the
//!   order of `DefectCode::ALL`), `-1` for none; a face with several codes
//!   shows its most severe.
//! - `severity`: 0 info, 1 warning, 2 error, 3 fatal, `-1` none.
//! - `part`: 1-based part index of a triangle, 0 for vertex and line cells.

use std::fmt::Write as _;

use mesher_core::model::SurfaceModel;
use mesher_core::overlay::Overlay;
use mesher_core::{DefectCode, Severity};

fn code_index(code: DefectCode) -> i32 {
    DefectCode::ALL
        .iter()
        .position(|&c| c == code)
        .map_or(-1, |i| i as i32)
}

fn severity_index(code: DefectCode) -> i32 {
    match code.default_severity() {
        Severity::Info => 0,
        Severity::Warning => 1,
        Severity::Error => 2,
        Severity::Fatal => 3,
    }
}

/// Render the overlay. `part_of_triangle` holds each triangle's 1-based part
/// index (0 if none). Deterministic: same inputs, same bytes.
#[must_use]
pub fn overlay_vtp(model: &SurfaceModel, part_of_triangle: &[u32], overlay: &Overlay) -> String {
    let n_points = model.vertices.len();
    let n_tris = model.triangles.len();

    // Most severe code per face (ties: earlier in the taxonomy).
    let mut face_code: Vec<Option<DefectCode>> = vec![None; n_tris];
    for &(f, code) in &overlay.faces {
        if let Some(slot) = face_code.get_mut(f as usize) {
            let better = slot.is_none_or(|c| {
                (severity_index(code), -code_index(code)) > (severity_index(c), -code_index(c))
            });
            if better {
                *slot = Some(code);
            }
        }
    }
    let mut verts: Vec<(u32, DefectCode)> = overlay
        .vertices
        .iter()
        .copied()
        .filter(|&(v, _)| (v as usize) < n_points)
        .collect();
    verts.sort_unstable();
    verts.dedup();
    let mut lines: Vec<([u32; 2], DefectCode)> = overlay
        .edges
        .iter()
        .copied()
        .filter(|&([a, b], _)| (a as usize) < n_points && (b as usize) < n_points)
        .collect();
    lines.sort_unstable();
    lines.dedup();

    let mut s = String::new();
    let _ = writeln!(s, r#"<?xml version="1.0"?>"#);
    let _ = writeln!(
        s,
        "<!-- ExodusMesher defect overlay. defect codes: {} -->",
        DefectCode::ALL
            .iter()
            .enumerate()
            .map(|(i, c)| format!("{i}={}", c.as_str()))
            .collect::<Vec<_>>()
            .join(" ")
    );
    let _ = writeln!(
        s,
        r#"<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">"#
    );
    let _ = writeln!(s, "  <PolyData>");
    let _ = writeln!(
        s,
        r#"    <Piece NumberOfPoints="{n_points}" NumberOfVerts="{}" NumberOfLines="{}" NumberOfStrips="0" NumberOfPolys="{n_tris}">"#,
        verts.len(),
        lines.len()
    );

    let _ = writeln!(s, "      <Points>");
    let _ = write!(
        s,
        r#"        <DataArray type="Float64" NumberOfComponents="3" format="ascii">"#
    );
    for p in &model.vertices {
        let _ = write!(s, " {} {} {}", p[0], p[1], p[2]);
    }
    let _ = writeln!(s, " </DataArray>");
    let _ = writeln!(s, "      </Points>");

    let mut section = |tag: &str, conn: Vec<u32>, per_cell: usize| {
        if conn.is_empty() {
            return;
        }
        let _ = writeln!(s, "      <{tag}>");
        let _ = write!(
            s,
            r#"        <DataArray type="Int64" Name="connectivity" format="ascii">"#
        );
        for c in &conn {
            let _ = write!(s, " {c}");
        }
        let _ = writeln!(s, " </DataArray>");
        let _ = write!(
            s,
            r#"        <DataArray type="Int64" Name="offsets" format="ascii">"#
        );
        for k in 1..=conn.len() / per_cell {
            let _ = write!(s, " {}", k * per_cell);
        }
        let _ = writeln!(s, " </DataArray>");
        let _ = writeln!(s, "      </{tag}>");
    };
    section("Verts", verts.iter().map(|v| v.0).collect(), 1);
    section("Lines", lines.iter().flat_map(|l| l.0).collect(), 2);
    section(
        "Polys",
        model.triangles.iter().flat_map(|t| t.vertices).collect(),
        3,
    );

    let defect: Vec<i32> = verts
        .iter()
        .map(|v| code_index(v.1))
        .chain(lines.iter().map(|l| code_index(l.1)))
        .chain(face_code.iter().map(|c| c.map_or(-1, code_index)))
        .collect();
    let severity: Vec<i32> = verts
        .iter()
        .map(|v| severity_index(v.1))
        .chain(lines.iter().map(|l| severity_index(l.1)))
        .chain(face_code.iter().map(|c| c.map_or(-1, severity_index)))
        .collect();
    let part: Vec<u32> = std::iter::repeat_n(0, verts.len() + lines.len())
        .chain((0..n_tris).map(|t| part_of_triangle.get(t).copied().unwrap_or(0)))
        .collect();
    let _ = writeln!(s, r#"      <CellData Scalars="defect">"#);
    for (name, values) in [
        (
            "defect",
            defect.iter().map(|&v| i64::from(v)).collect::<Vec<_>>(),
        ),
        ("severity", severity.iter().map(|&v| i64::from(v)).collect()),
        ("part", part.iter().map(|&p| i64::from(p)).collect()),
    ] {
        let _ = write!(
            s,
            r#"        <DataArray type="Int32" Name="{name}" format="ascii">"#
        );
        for v in values {
            let _ = write!(s, " {v}");
        }
        let _ = writeln!(s, " </DataArray>");
    }
    let _ = writeln!(s, "      </CellData>");
    let _ = writeln!(s, "    </Piece>");
    let _ = writeln!(s, "  </PolyData>");
    let _ = writeln!(s, "</VTKFile>");
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use mesher_core::model::Triangle;

    #[test]
    fn counts_and_arrays_line_up() {
        let model = SurfaceModel {
            vertices: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            triangles: vec![Triangle {
                vertices: [0, 1, 2],
                polygon: 0,
            }],
            ..SurfaceModel::default()
        };
        let mut overlay = Overlay::default();
        overlay.edge(1, 0, DefectCode::TopBoundaryEdge);
        overlay.faces([0], DefectCode::QualSliverFace);
        overlay.faces([0], DefectCode::IntSelfIntersection);
        let text = overlay_vtp(&model, &[1], &overlay);
        assert!(text.contains(r#"NumberOfVerts="0" NumberOfLines="1""#));
        assert!(!text.contains("<Verts>"));
        // The line cell, then the face showing its most severe code (an error).
        let sev = code_index(DefectCode::IntSelfIntersection);
        let line = code_index(DefectCode::TopBoundaryEdge);
        assert!(text.contains(&format!(
            r#"Name="defect" format="ascii"> {line} {sev} </DataArray>"#
        )));
        assert!(text.contains(r#"Name="part" format="ascii"> 0 1 </DataArray>"#));
        assert_eq!(text, overlay_vtp(&model, &[1], &overlay));
    }
}
