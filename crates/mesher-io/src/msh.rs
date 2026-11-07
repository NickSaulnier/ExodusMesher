//! Gmsh MSH 4.1 (ASCII) writer.
//!
//! Each element block becomes a volume entity carrying a physical volume of
//! the same name (tag = block number); each side set becomes a surface entity
//! carrying a physical surface (tag = set number) with one outward triangle
//! per side. All nodes are listed in one block on the first volume entity.
//! Gmsh's tetrahedron node order matches ours (positive orientation).

use std::fmt::Write as _;

use mesher_core::mesh::{VolumeMesh, tet_side_nodes};

/// Gmsh element type numbers.
const GMSH_TRIANGLE: u32 = 2;
const GMSH_TETRAHEDRON: u32 = 4;

fn bbox(points: impl Iterator<Item = [f64; 3]>) -> [f64; 6] {
    let mut b = [
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    let mut any = false;
    for p in points {
        any = true;
        for k in 0..3 {
            b[k] = b[k].min(p[k]);
            b[k + 3] = b[k + 3].max(p[k]);
        }
    }
    if any { b } else { [0.0; 6] }
}

fn quoted(name: &str) -> String {
    format!("\"{}\"", name.replace('"', "'"))
}

/// Render `mesh` as MSH 4.1 text. Deterministic.
#[must_use]
pub fn msh_text(mesh: &VolumeMesh) -> String {
    let elements: Vec<[u32; 4]> = mesh.elements().collect();
    let point = |v: u32| mesh.points.get(v as usize).copied().unwrap_or([0.0; 3]);
    let set_faces: Vec<Vec<[u32; 3]>> = mesh
        .side_sets
        .iter()
        .map(|set| {
            set.sides
                .iter()
                .filter_map(|&(e, s)| elements.get(e as usize).and_then(|&t| tet_side_nodes(t, s)))
                .collect()
        })
        .collect();

    let mut s = String::new();
    let _ = writeln!(s, "$MeshFormat\n4.1 0 8\n$EndMeshFormat");

    let _ = writeln!(s, "$PhysicalNames");
    let _ = writeln!(s, "{}", mesh.blocks.len() + mesh.side_sets.len());
    for (i, set) in mesh.side_sets.iter().enumerate() {
        let _ = writeln!(s, "2 {} {}", i + 1, quoted(&set.name));
    }
    for (i, block) in mesh.blocks.iter().enumerate() {
        let _ = writeln!(s, "3 {} {}", i + 1, quoted(&block.name));
    }
    let _ = writeln!(s, "$EndPhysicalNames");

    let _ = writeln!(s, "$Entities");
    let _ = writeln!(s, "0 0 {} {}", mesh.side_sets.len(), mesh.blocks.len());
    for (i, faces) in set_faces.iter().enumerate() {
        let b = bbox(faces.iter().flatten().map(|&v| point(v)));
        let _ = writeln!(
            s,
            "{} {} {} {} {} {} {} 1 {} 0",
            i + 1,
            b[0],
            b[1],
            b[2],
            b[3],
            b[4],
            b[5],
            i + 1
        );
    }
    for (i, block) in mesh.blocks.iter().enumerate() {
        let b = bbox(block.tets.iter().flatten().map(|&v| point(v)));
        let _ = writeln!(
            s,
            "{} {} {} {} {} {} {} 1 {} 0",
            i + 1,
            b[0],
            b[1],
            b[2],
            b[3],
            b[4],
            b[5],
            i + 1
        );
    }
    let _ = writeln!(s, "$EndEntities");

    let n = mesh.points.len();
    let _ = writeln!(s, "$Nodes");
    if n == 0 {
        let _ = writeln!(s, "0 0 0 0");
    } else {
        let _ = writeln!(s, "1 {n} 1 {n}");
        let _ = writeln!(s, "3 1 0 {n}");
        for i in 1..=n {
            let _ = writeln!(s, "{i}");
        }
        for p in &mesh.points {
            let _ = writeln!(s, "{} {} {}", p[0], p[1], p[2]);
        }
    }
    let _ = writeln!(s, "$EndNodes");

    let blocks = mesh.blocks.iter().filter(|b| !b.tets.is_empty()).count()
        + set_faces.iter().filter(|f| !f.is_empty()).count();
    let total = elements.len() + set_faces.iter().map(Vec::len).sum::<usize>();
    let _ = writeln!(s, "$Elements");
    let _ = writeln!(s, "{blocks} {total} {} {total}", u32::from(total > 0));
    let mut tag = 1usize;
    for (i, block) in mesh.blocks.iter().enumerate() {
        if block.tets.is_empty() {
            continue;
        }
        let _ = writeln!(s, "3 {} {GMSH_TETRAHEDRON} {}", i + 1, block.tets.len());
        for t in &block.tets {
            let _ = writeln!(
                s,
                "{tag} {} {} {} {}",
                t[0] + 1,
                t[1] + 1,
                t[2] + 1,
                t[3] + 1
            );
            tag += 1;
        }
    }
    for (i, faces) in set_faces.iter().enumerate() {
        if faces.is_empty() {
            continue;
        }
        let _ = writeln!(s, "2 {} {GMSH_TRIANGLE} {}", i + 1, faces.len());
        for f in faces {
            let _ = writeln!(s, "{tag} {} {} {}", f[0] + 1, f[1] + 1, f[2] + 1);
            tag += 1;
        }
    }
    let _ = writeln!(s, "$EndElements");
    s
}
