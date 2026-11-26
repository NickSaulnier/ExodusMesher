//! Fluent mesh file writer (ASCII `.msh`), from the format's public
//! description (the Fluent User's Guide, appendix "Case and data file
//! formats": grid and zone sections).
//!
//! Sections, in order: a comment (0), the dimensions (2: 3), the node, cell
//! and face declarations (10, 12 and 13 with zone 0: the totals), then one
//! node zone (type 1, every point a cell uses, in order), one cell zone of
//! tetrahedra (element type 2, type 1: fluid) per element block, one
//! interior face zone (bc-type 2) with every internal face, a wall face zone
//! (bc-type 3) per [`FaceMesh`] patch, and a zone section (45) naming each
//! zone. Indices are 1-based. Counts, indices and the grid sections' zone
//! ids are hexadecimal, the zone sections' ids decimal, coordinates decimal.
//!
//! Each face line is its three nodes and its cells `c0 c1`, the nodes
//! ordered so that the right-hand rule points to `c0`: `c0` is the face's
//! owner in [`VolumeMesh::faces`] (whose faces are listed out of the owner,
//! so here reversed), `c1` the neighbour, or 0 on the boundary.
//!
//! Zone names must be Scheme symbols, here lowercase; a name that is not
//! one, or that would repeat another zone's, is changed and recorded as
//! `OUT_NAME_CHANGED`. The interior zone is `interior`.

use std::fmt::Write as _;

use mesher_core::RecordLog;
use mesher_core::mesh::{FaceMesh, MeshProblem, VolumeMesh};

use crate::shared::{NameRule, names, used_points};

/// Fluent's codes, from its file format description.
const ELEMENT_TETRAHEDRON: u32 = 2;
const CELL_FLUID: u32 = 1;
const NODE_ANY: u32 = 1;
const FACE_TRIANGLE: u32 = 3;
const BC_INTERIOR: u32 = 2;
const BC_WALL: u32 = 3;

/// Render `mesh` as a Fluent mesh file. Deterministic.
///
/// # Errors
/// The mesh's [`MeshProblem`] when it has no face view (a face on more than
/// two cells: a bug upstream, as the volume gate rejects such meshes).
pub fn fluent_text(mesh: &VolumeMesh, log: &mut RecordLog) -> Result<String, MeshProblem> {
    let f = mesh.faces()?;
    Ok(render(mesh, &f, log))
}

fn faces(s: &mut String, f: &FaceMesh, new: &[u32], range: std::ops::Range<usize>) {
    let node = |v: u32| new.get(v as usize).copied().unwrap_or(0) + 1;
    for k in range {
        let [a, b, c] = f.faces[k];
        let c0 = f.owner[k] + 1;
        let c1 = f.neighbour.get(k).map_or(0, |n| n + 1);
        // Reversed, so the right-hand rule points into the owner, c0.
        let _ = writeln!(s, "{:x} {:x} {:x} {c0:x} {c1:x}", node(a), node(c), node(b));
    }
}

fn render(mesh: &VolumeMesh, f: &FaceMesh, log: &mut RecordLog) -> String {
    let (new, kept) = used_points(mesh.points.len(), f);
    let mut s = String::new();
    let _ = writeln!(s, "(0 \"Written by ExodusMesher\")");
    let _ = writeln!(s, "(2 3)");
    let _ = writeln!(s, "(10 (0 1 {:x} 0 3))", kept.len());
    let _ = writeln!(s, "(12 (0 1 {:x} 0))", f.cells);
    let _ = writeln!(s, "(13 (0 1 {:x} 0))", f.faces.len());

    // Zone ids: the nodes 1, then the cell zones, the interior, the patches.
    let mut next = 1u32;
    let mut id = || {
        next += 1;
        next
    };
    let node_zone = 1u32;
    let _ = writeln!(s, "(10 ({node_zone:x} 1 {:x} {NODE_ANY} 3)(", kept.len());
    for &v in &kept {
        let [x, y, z] = mesh.points[v];
        let _ = writeln!(s, "{x:?} {y:?} {z:?}");
    }
    let _ = writeln!(s, "))");

    let mut zones: Vec<(u32, &str, &str, &str)> = Vec::new();
    for z in &f.zones {
        if z.count == 0 {
            continue;
        }
        let zid = id();
        let _ = writeln!(
            s,
            "(12 ({zid:x} {:x} {:x} {CELL_FLUID} {ELEMENT_TETRAHEDRON}))",
            z.start + 1,
            z.start + z.count
        );
        zones.push((zid, "fluid", "cell_zone", &z.name));
    }
    if f.internal() > 0 {
        let zid = id();
        let _ = writeln!(
            s,
            "(13 ({zid:x} 1 {:x} {BC_INTERIOR} {FACE_TRIANGLE})(",
            f.internal()
        );
        faces(&mut s, f, &new, 0..f.internal());
        let _ = writeln!(s, "))");
        zones.push((zid, "interior", "interior", "interior"));
    }
    for p in &f.patches {
        let zid = id();
        let _ = writeln!(
            s,
            "(13 ({zid:x} {:x} {:x} {BC_WALL} {FACE_TRIANGLE})(",
            p.start + 1,
            p.start + p.count
        );
        faces(&mut s, f, &new, p.start..p.start + p.count);
        let _ = writeln!(s, "))");
        zones.push((zid, "wall", "patch", &p.name));
    }

    // Names: the interior's first, so it keeps `interior`.
    let mut order: Vec<usize> = (0..zones.len()).collect();
    order.sort_by_key(|&i| zones[i].1 != "interior");
    let items: Vec<(&str, &str)> = order.iter().map(|&i| (zones[i].2, zones[i].3)).collect();
    let written = names(NameRule::Fluent, &items, log);
    let mut name = vec![String::new(); zones.len()];
    for (&i, w) in order.iter().zip(written) {
        name[i] = w;
    }
    for ((zid, kind, _, _), n) in zones.iter().zip(&name) {
        let _ = writeln!(s, "(45 ({zid} {kind} {n})())");
    }
    s
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;
    use mesher_core::mesh::{Block, SideSet};
    use mesher_geom::vec3::{cross, dot, sub};

    /// Two tets sharing the face (1, 2, 3); tet 1's other faces are the
    /// group `Top`.
    fn two_tets() -> VolumeMesh {
        let mut mesh = VolumeMesh {
            points: vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
            ],
            blocks: vec![Block {
                name: "b".into(),
                tets: vec![[0, 1, 2, 3], [1, 4, 2, 3]],
            }],
            side_sets: Vec::new(),
        };
        mesh.add_boundary_side_sets();
        mesh.side_sets.push(SideSet {
            name: "Top".into(),
            sides: vec![(1, 1), (1, 2), (1, 4)],
        });
        mesh
    }

    #[test]
    fn two_tets_write_exactly() {
        let mesh = two_tets();
        let mut log = RecordLog::new();
        let text = fluent_text(&mesh, &mut log).unwrap();
        let f = mesh.faces().unwrap();
        let line = |k: usize| {
            let [a, b, c] = f.faces[k];
            let n = f.neighbour.get(k).map_or(0, |n| n + 1);
            format!(
                "{:x} {:x} {:x} {:x} {n:x}\n",
                a + 1,
                c + 1,
                b + 1,
                f.owner[k] + 1
            )
        };
        let mut want = String::from(
            "(0 \"Written by ExodusMesher\")\n(2 3)\n(10 (0 1 5 0 3))\n(12 (0 1 2 0))\n(13 (0 1 7 0))\n\
             (10 (1 1 5 1 3)(\n0.0 0.0 0.0\n1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1.0 1.0 1.0\n))\n\
             (12 (2 1 2 1 2))\n(13 (3 1 1 2 3)(\n",
        );
        want += &line(0);
        want += "))\n(13 (4 2 4 3 3)(\n";
        for k in 1..4 {
            want += &line(k);
        }
        want += "))\n(13 (5 5 7 3 3)(\n";
        for k in 4..7 {
            want += &line(k);
        }
        want += "))\n(45 (2 fluid b)())\n(45 (3 interior interior)())\n(45 (4 wall b_boundary)())\n(45 (5 wall top)())\n";
        assert_eq!(text, want);
        // `Top` is lowered, and recorded.
        assert_eq!(log.records().len(), 1);
        assert_eq!(
            log.records()[0].measured["written_as"].as_str(),
            Some("top")
        );
    }

    /// The right-hand rule of each written face points to its c0.
    #[test]
    fn each_face_points_to_its_c0() {
        let mesh = two_tets();
        let text = fluent_text(&mesh, &mut RecordLog::new()).unwrap();
        let tets: Vec<[u32; 4]> = mesh.elements().collect();
        let centroid = |c: usize| {
            let p = tets[c].map(|v| mesh.points[v as usize]);
            [0, 1, 2].map(|m| (p[0][m] + p[1][m] + p[2][m] + p[3][m]) / 4.0)
        };
        let mut checked = 0;
        for l in text.lines() {
            let v: Vec<u32> = l
                .split(' ')
                .map(|t| u32::from_str_radix(t, 16))
                .collect::<Result<_, _>>()
                .unwrap_or_default();
            if v.len() != 5 || l.contains('.') {
                continue;
            }
            let [a, b, c] = [v[0], v[1], v[2]].map(|n| mesh.points[n as usize - 1]);
            let normal = cross(sub(b, a), sub(c, a));
            let c0 = v[3] as usize - 1;
            assert!(dot(normal, sub(centroid(c0), a)) > 0.0, "{l}");
            if v[4] > 0 {
                assert!(
                    dot(normal, sub(centroid(v[4] as usize - 1), a)) < 0.0,
                    "{l}"
                );
            }
            checked += 1;
        }
        assert_eq!(checked, 7);
    }

    /// Zone ids over 9 are hexadecimal in the grid sections and decimal in
    /// the zone sections.
    #[test]
    fn zone_ids_are_hex_in_grid_sections_and_decimal_in_zone_sections() {
        let mut mesh = two_tets();
        // Each of block b's boundary faces its own patch, p0 to p5.
        mesh.side_sets.truncate(1);
        let sides = mesh.side_sets[0].sides.clone();
        for (i, &side) in sides.iter().enumerate() {
            mesh.side_sets.push(SideSet {
                name: format!("p{i}"),
                sides: vec![side],
            });
        }
        // A second block, a tet of its own.
        let base = mesh.points.len() as u32;
        mesh.points.extend([
            [5.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [5.0, 1.0, 0.0],
            [5.0, 0.0, 1.0],
        ]);
        mesh.blocks.push(Block {
            name: "c".into(),
            tets: vec![[base, base + 1, base + 2, base + 3]],
        });
        let c = mesh.boundary_sides(1);
        mesh.side_sets.push(SideSet {
            name: "c_boundary".into(),
            sides: c,
        });
        let text = fluent_text(&mesh, &mut RecordLog::new()).unwrap();
        // Zones: nodes 1, cells 2 and 3, the interior 4, p0 to p5 5 to 10,
        // c_boundary 11 (b_boundary keeps no face, so is left out).
        assert!(text.contains("(12 (3 3 3 1 2))"), "{text}");
        assert!(text.contains("(13 (a 7 7 3 3)("), "{text}");
        assert!(text.contains("(13 (b 8 b 3 3)("), "{text}");
        assert!(text.contains("(45 (10 wall p5)())"), "{text}");
        assert!(text.contains("(45 (11 wall c_boundary)())"), "{text}");
        assert!(!text.contains("b_boundary"), "{text}");
    }
}
