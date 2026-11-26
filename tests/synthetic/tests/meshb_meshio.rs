//! Milestone 19, stage 5: our reader on meshio's `.meshb` files. With
//! `MESHB_MESHIO_DIR=<root>`, every case directory there holding both our
//! `mesh.meshb` and meshio's `meshio.meshb` (written by
//! `tests/meshb/meshio_roundtrip.py` from what it read of ours) must read to
//! the same vertices, triangles, tetrahedra and references. Without the
//! variable the test does nothing (CI sets it after the script).
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_io::meshb::{MeshbFile, read_meshb};

fn read(path: &std::path::Path) -> MeshbFile {
    let bytes = std::fs::read(path).unwrap();
    let (file, error) = read_meshb(&bytes);
    assert_eq!(error, None, "{}", path.display());
    file
}

#[test]
fn meshio_files_read_as_ours() {
    let Some(root) = std::env::var_os("MESHB_MESHIO_DIR").filter(|d| !d.is_empty()) else {
        return;
    };
    let mut cases: Vec<std::path::PathBuf> = std::fs::read_dir(&root)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.join("mesh.meshb").exists() && p.join("meshio.meshb").exists())
        .collect();
    cases.sort();
    assert!(!cases.is_empty(), "no meshio.meshb files under {root:?}");
    let mut versions = std::collections::BTreeSet::new();
    for case in &cases {
        let ours = read(&case.join("mesh.meshb"));
        let theirs = read(&case.join("meshio.meshb"));
        versions.insert(theirs.version);
        let name = case.display();
        assert_eq!(theirs.dimension, 3, "{name}");
        assert_eq!(theirs.vertices, ours.vertices, "{name}: vertices");
        assert_eq!(
            theirs.vertex_refs, ours.vertex_refs,
            "{name}: vertex references"
        );
        assert_eq!(theirs.triangles, ours.triangles, "{name}: triangles");
        assert_eq!(theirs.tetrahedra, ours.tetrahedra, "{name}: tetrahedra");
    }
    eprintln!(
        "{} meshio files read as ours (meshio's versions {versions:?})",
        cases.len()
    );
}
