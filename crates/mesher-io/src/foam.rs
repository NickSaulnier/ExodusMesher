//! OpenFOAM polyMesh writer (ASCII), from the format's public description.
//!
//! A case directory: `constant/polyMesh/` holds `points`, `faces` (each a
//! triangle, listed out of its owner cell), `owner` (per face), `neighbour`
//! (per internal face), `boundary` (a patch of type `patch` per
//! [`FaceMesh`] patch: its face count and first face) and `cellZones` (a zone
//! per element block), and `system/` holds the least a utility such as
//! `checkMesh` needs to run on the case as written: `controlDict`, and
//! `fvSchemes` and `fvSolution` with generic defaults (a solver's case
//! replaces them). The face order is
//! [`VolumeMesh::faces`]': internal faces first, by owner and then
//! neighbour (upper-triangular order), then each patch's faces together.
//! Points no cell uses are left out, the rest keep their order.
//!
//! Patch and zone names must be OpenFOAM words; a name that is not one (a
//! space, a quote, `/`, `;`, braces, brackets or parentheses, a leading digit
//! or sign) or that would repeat another is changed, and recorded as
//! `OUT_NAME_CHANGED`.

use std::fmt::Write as _;

use mesher_core::RecordLog;
use mesher_core::mesh::{FaceMesh, MeshProblem, VolumeMesh};

use crate::shared::{NameRule, used_points};

/// The files of a case, relative to its directory, in the order written.
pub const POINTS: &str = "constant/polyMesh/points";
pub const FACES: &str = "constant/polyMesh/faces";
pub const OWNER: &str = "constant/polyMesh/owner";
pub const NEIGHBOUR: &str = "constant/polyMesh/neighbour";
pub const BOUNDARY: &str = "constant/polyMesh/boundary";
pub const CELL_ZONES: &str = "constant/polyMesh/cellZones";
pub const CONTROL_DICT: &str = "system/controlDict";
pub const FV_SCHEMES: &str = "system/fvSchemes";
pub const FV_SOLUTION: &str = "system/fvSolution";

/// Words for `names`, distinct, each recorded where it changed.
fn words(entity: &str, names: &[&str], log: &mut RecordLog) -> Vec<String> {
    let items: Vec<(&str, &str)> = names.iter().map(|&n| (entity, n)).collect();
    crate::shared::names(NameRule::Foam, &items, log)
}

fn header(s: &mut String, class: &str, location: &str, object: &str, note: Option<&str>) {
    let _ = writeln!(s, "// Written by ExodusMesher.");
    let _ = writeln!(s, "FoamFile");
    let _ = writeln!(s, "{{");
    let _ = writeln!(s, "    version     2.0;");
    let _ = writeln!(s, "    format      ascii;");
    let _ = writeln!(s, "    class       {class};");
    if let Some(note) = note {
        let _ = writeln!(s, "    note        \"{note}\";");
    }
    let _ = writeln!(s, "    location    \"{location}\";");
    let _ = writeln!(s, "    object      {object};");
    let _ = writeln!(s, "}}");
    let _ = writeln!(s);
}

/// A list of labels, one per line.
fn label_list(s: &mut String, labels: impl ExactSizeIterator<Item = u32>) {
    let _ = writeln!(s, "{}", labels.len());
    let _ = writeln!(s, "(");
    for l in labels {
        let _ = writeln!(s, "{l}");
    }
    let _ = writeln!(s, ")");
}

/// Render `mesh` as the files of an OpenFOAM case, paths relative to the
/// case directory. Deterministic.
///
/// # Errors
/// The mesh's [`MeshProblem`] when it has no face view (a face on more than
/// two cells: a bug upstream, as the volume gate rejects such meshes).
pub fn foam_files(
    mesh: &VolumeMesh,
    log: &mut RecordLog,
) -> Result<Vec<(&'static str, String)>, MeshProblem> {
    let f = mesh.faces()?;
    Ok(render(mesh, &f, log))
}

fn render(mesh: &VolumeMesh, f: &FaceMesh, log: &mut RecordLog) -> Vec<(&'static str, String)> {
    // Points in use, renumbered in order.
    let (new, kept) = used_points(mesh.points.len(), f);
    let note = format!(
        "nPoints:{} nCells:{} nFaces:{} nInternalFaces:{}",
        kept.len(),
        f.cells,
        f.faces.len(),
        f.internal()
    );
    let location = "constant/polyMesh";

    let mut points = String::new();
    header(&mut points, "vectorField", location, "points", None);
    let _ = writeln!(points, "{}", kept.len());
    let _ = writeln!(points, "(");
    for &v in &kept {
        let [x, y, z] = mesh.points[v];
        let _ = writeln!(points, "({x:?} {y:?} {z:?})");
    }
    let _ = writeln!(points, ")");

    let mut faces = String::new();
    header(&mut faces, "faceList", location, "faces", None);
    let _ = writeln!(faces, "{}", f.faces.len());
    let _ = writeln!(faces, "(");
    for face in &f.faces {
        let [a, b, c] = face.map(|v| new.get(v as usize).copied().unwrap_or(0));
        let _ = writeln!(faces, "3({a} {b} {c})");
    }
    let _ = writeln!(faces, ")");

    let mut owner = String::new();
    header(&mut owner, "labelList", location, "owner", Some(&note));
    label_list(&mut owner, f.owner.iter().copied());

    let mut neighbour = String::new();
    header(
        &mut neighbour,
        "labelList",
        location,
        "neighbour",
        Some(&note),
    );
    label_list(&mut neighbour, f.neighbour.iter().copied());

    let patch_names: Vec<&str> = f.patches.iter().map(|p| p.name.as_str()).collect();
    let patch_words = words("patch", &patch_names, log);
    let mut boundary = String::new();
    header(
        &mut boundary,
        "polyBoundaryMesh",
        location,
        "boundary",
        None,
    );
    let _ = writeln!(boundary, "{}", f.patches.len());
    let _ = writeln!(boundary, "(");
    for (p, name) in f.patches.iter().zip(&patch_words) {
        let _ = writeln!(boundary, "    {name}");
        let _ = writeln!(boundary, "    {{");
        let _ = writeln!(boundary, "        type            patch;");
        let _ = writeln!(boundary, "        nFaces          {};", p.count);
        let _ = writeln!(boundary, "        startFace       {};", p.start);
        let _ = writeln!(boundary, "    }}");
    }
    let _ = writeln!(boundary, ")");

    let zone_names: Vec<&str> = f.zones.iter().map(|z| z.name.as_str()).collect();
    let zone_words = words("cell_zone", &zone_names, log);
    let mut zones = String::new();
    header(&mut zones, "regIOobject", location, "cellZones", None);
    let _ = writeln!(zones, "{}", f.zones.len());
    let _ = writeln!(zones, "(");
    for (z, name) in f.zones.iter().zip(&zone_words) {
        let _ = writeln!(zones, "{name}");
        let _ = writeln!(zones, "{{");
        let _ = writeln!(zones, "    type cellZone;");
        let _ = write!(zones, "    cellLabels List<label> ");
        label_list(&mut zones, (z.start..z.start + z.count).map(|c| c as u32));
        let _ = writeln!(zones, ";");
        let _ = writeln!(zones, "}}");
    }
    let _ = writeln!(zones, ")");

    let mut control = String::new();
    header(&mut control, "dictionary", "system", "controlDict", None);
    for line in [
        "application     none;",
        "startFrom       startTime;",
        "startTime       0;",
        "stopAt          endTime;",
        "endTime         0;",
        "deltaT          1;",
        "writeControl    timeStep;",
        "writeInterval   1;",
        "writeFormat     ascii;",
        "writePrecision  17;",
    ] {
        let _ = writeln!(control, "{line}");
    }

    let mut schemes = String::new();
    header(&mut schemes, "dictionary", "system", "fvSchemes", None);
    for line in [
        "ddtSchemes           { default steadyState; }",
        "gradSchemes          { default Gauss linear; }",
        "divSchemes           { default none; }",
        "laplacianSchemes     { default Gauss linear corrected; }",
        "interpolationSchemes { default linear; }",
        "snGradSchemes        { default corrected; }",
    ] {
        let _ = writeln!(schemes, "{line}");
    }
    let mut solution = String::new();
    header(&mut solution, "dictionary", "system", "fvSolution", None);
    let _ = writeln!(solution, "solvers {{}}");

    vec![
        (POINTS, points),
        (FACES, faces),
        (OWNER, owner),
        (NEIGHBOUR, neighbour),
        (BOUNDARY, boundary),
        (CELL_ZONES, zones),
        (CONTROL_DICT, control),
        (FV_SCHEMES, schemes),
        (FV_SOLUTION, solution),
    ]
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;
    use mesher_core::mesh::{Block, SideSet};

    fn file<'a>(files: &'a [(&'static str, String)], name: &str) -> &'a str {
        let text = &files.iter().find(|(n, _)| *n == name).unwrap().1;
        // The data after the header.
        text.split_once("}\n\n").unwrap().1
    }

    /// Two tets sharing the face (1, 2, 3); tet 1's other faces are the
    /// group `top`, and the second block's name is not a word.
    fn two_tets() -> VolumeMesh {
        let mut mesh = VolumeMesh {
            points: vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
                // Unused.
                [9.0, 9.0, 9.0],
            ],
            blocks: vec![Block {
                name: "b".into(),
                tets: vec![[0, 1, 2, 3], [1, 4, 2, 3]],
            }],
            side_sets: Vec::new(),
        };
        mesh.add_boundary_side_sets();
        mesh.side_sets.push(SideSet {
            name: "top".into(),
            sides: vec![(1, 1), (1, 2), (1, 4)],
        });
        mesh
    }

    #[test]
    fn two_tets_write_exactly() {
        let mut log = RecordLog::new();
        let files = foam_files(&two_tets(), &mut log).unwrap();
        assert!(log.records().is_empty());
        let names: Vec<&str> = files.iter().map(|(n, _)| *n).collect();
        assert_eq!(
            names,
            [
                POINTS,
                FACES,
                OWNER,
                NEIGHBOUR,
                BOUNDARY,
                CELL_ZONES,
                CONTROL_DICT,
                FV_SCHEMES,
                FV_SOLUTION
            ]
        );
        assert_eq!(
            file(&files, POINTS),
            "5\n(\n(0.0 0.0 0.0)\n(1.0 0.0 0.0)\n(0.0 1.0 0.0)\n(0.0 0.0 1.0)\n(1.0 1.0 1.0)\n)\n"
        );
        let f = two_tets().faces().unwrap();
        let mut faces = format!("{}\n(\n", f.faces.len());
        for [a, b, c] in &f.faces {
            faces += &format!("3({a} {b} {c})\n");
        }
        faces += ")\n";
        assert_eq!(file(&files, FACES), faces);
        assert_eq!(f.faces.len(), 7);
        assert_eq!(
            file(&files, OWNER),
            "7\n(\n0\n0\n0\n0\n1\n1\n1\n)\n",
            "the shared face, tet 0's three, tet 1's three (`top`)"
        );
        assert_eq!(file(&files, NEIGHBOUR), "1\n(\n1\n)\n");
        assert_eq!(
            file(&files, BOUNDARY),
            "2\n(\n    b_boundary\n    {\n        type            patch;\n        nFaces          3;\n        startFace       1;\n    }\n    top\n    {\n        type            patch;\n        nFaces          3;\n        startFace       4;\n    }\n)\n"
        );
        assert_eq!(
            file(&files, CELL_ZONES),
            "1\n(\nb\n{\n    type cellZone;\n    cellLabels List<label> 2\n(\n0\n1\n)\n;\n}\n)\n"
        );
        let owner = &files.iter().find(|(n, _)| *n == OWNER).unwrap().1;
        assert!(owner.contains("note        \"nPoints:5 nCells:2 nFaces:7 nInternalFaces:1\";"));
    }

    #[test]
    fn names_that_are_not_words_are_changed_and_recorded() {
        let mut log = RecordLog::new();
        let w = words(
            "patch",
            &["a b", "a_b", "top", "2nd", "x/y;{z}", "", "top", "é"],
            &mut log,
        );
        assert_eq!(
            w,
            [
                "a_b_2", "a_b", "top", "_2nd", "x_y__z_", "_", "top_2", "__2"
            ]
        );
        let changed: Vec<(String, String)> = log
            .records()
            .iter()
            .map(|r| {
                assert_eq!(r.code, mesher_core::DefectCode::OutNameChanged);
                let get = |k: &str| r.measured[k].as_str().unwrap().to_owned();
                assert_eq!(get("format"), "openfoam");
                (get("name"), get("written_as"))
            })
            .collect();
        assert_eq!(changed.len(), 6, "{changed:?}");
        assert_eq!(changed[0], ("a b".to_owned(), "a_b_2".to_owned()));
    }

    #[test]
    fn a_block_name_that_is_not_a_word_is_a_changed_zone_and_patch() {
        let mut mesh = two_tets();
        mesh.blocks[0].name = "my part".into();
        mesh.side_sets.retain(|s| s.name == "top");
        mesh.add_boundary_side_sets();
        let mut log = RecordLog::new();
        let files = foam_files(&mesh, &mut log).unwrap();
        assert!(file(&files, CELL_ZONES).contains("\nmy_part\n"));
        assert!(file(&files, BOUNDARY).contains("my_part_boundary"));
        let entities: Vec<&str> = log
            .records()
            .iter()
            .map(|r| r.measured["entity"].as_str().unwrap())
            .collect();
        assert_eq!(entities, ["patch", "cell_zone"]);
    }
}
