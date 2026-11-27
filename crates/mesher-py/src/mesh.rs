//! Reading a run's mesh back (Milestone 20 stage 4): the `.meshb` file the
//! worker writes (`mesh.meshb`, with `mesh_refs.txt` naming its
//! references), as flat arrays for Python.
//!
//! Tetrahedra carry their part, boundary triangles their patch, each as an
//! index into the names; vertex numbers become 0-based. Names come from the
//! references file beside the mesh (`<stem>_refs.txt`, so `mesh_refs.txt`
//! for `mesh.meshb`); a reference it does not name (or a file without one)
//! is `ref<k>`. Parts and patches are in order of their reference numbers.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use mesher_core::protocol::{MESH_MESHB_FILE, MESH_MESHB_REFS_FILE};
use mesher_io::meshb::{MeshbFile, read_meshb};

/// A volume mesh as read: 0-based, flat.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MeshData {
    pub points: Vec<[f64; 3]>,
    pub tets: Vec<[u32; 4]>,
    /// Per tetrahedron: its index in `parts`.
    pub tet_parts: Vec<u32>,
    pub parts: Vec<String>,
    pub triangles: Vec<[u32; 3]>,
    /// Per boundary triangle: its index in `patches`.
    pub triangle_patches: Vec<u32>,
    pub patches: Vec<String>,
}

/// Why a mesh could not be read.
#[derive(Debug)]
pub enum ReadError {
    /// The file could not be read (missing, unreadable).
    Io(PathBuf, std::io::Error),
    /// The file is not a mesh this reader can follow.
    Invalid(PathBuf, String),
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadError::Io(p, e) => write!(f, "{}: {e}", p.display()),
            ReadError::Invalid(p, e) => write!(f, "{}: {e}", p.display()),
        }
    }
}

impl std::error::Error for ReadError {}

/// The mesh file and its references file for `path`: a run's output
/// directory (its `mesh.meshb`) or a `.meshb` file.
#[must_use]
pub fn mesh_files(path: &Path) -> (PathBuf, PathBuf) {
    if path.is_dir() {
        return (path.join(MESH_MESHB_FILE), path.join(MESH_MESHB_REFS_FILE));
    }
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("mesh");
    let refs = path.with_file_name(format!("{stem}_refs.txt"));
    (path.to_path_buf(), refs)
}

/// Read the mesh at `path` (a run's output directory or a `.meshb` file).
///
/// # Errors
/// The file cannot be read, or is not a complete volume mesh (a reading
/// problem, a vertex number out of range).
pub fn read_mesh(path: &Path) -> Result<MeshData, ReadError> {
    let (file, refs) = mesh_files(path);
    let bytes = std::fs::read(&file).map_err(|e| ReadError::Io(file.clone(), e))?;
    // A missing references file leaves the references unnamed.
    let names = match std::fs::read_to_string(&refs) {
        Ok(text) => text,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => String::new(),
        Err(e) => return Err(ReadError::Io(refs, e)),
    };
    mesh_from_bytes(&bytes, &names).map_err(|e| ReadError::Invalid(file, e))
}

/// The mesh in `.meshb` bytes, named by the references file's text.
///
/// # Errors
/// As [`read_mesh`], as a message.
pub fn mesh_from_bytes(bytes: &[u8], refs: &str) -> Result<MeshData, String> {
    let (file, error) = read_meshb(bytes);
    if let Some(e) = error {
        return Err(format!("at byte {}: {}", e.offset, e.reason));
    }
    convert(&file, refs)
}

/// The names `refs` gives each reference of `kind` (`tetrahedra` or
/// `triangles`): lines `<kind> <ref> <name>`, the name the rest of the line.
fn names(refs: &str, kind: &str) -> BTreeMap<i64, String> {
    let mut out = BTreeMap::new();
    for line in refs.lines() {
        let mut it = line.splitn(3, ' ');
        if let (Some(k), Some(r), Some(name)) = (it.next(), it.next(), it.next())
            && k == kind
            && let Ok(r) = r.parse::<i64>()
        {
            out.entry(r).or_insert_with(|| name.to_owned());
        }
    }
    out
}

/// Each element's group index, and the groups' names: every reference used
/// or named (a named one may have no elements), in order of number.
fn groups(references: &[i64], named: &BTreeMap<i64, String>) -> (Vec<u32>, Vec<String>) {
    let used: BTreeMap<i64, u32> = {
        let mut refs: Vec<i64> = references.to_vec();
        refs.extend(named.keys());
        refs.sort_unstable();
        refs.dedup();
        refs.into_iter().zip(0u32..).collect()
    };
    let index = references
        .iter()
        .map(|r| used.get(r).copied().unwrap_or(0))
        .collect();
    let names = used
        .keys()
        .map(|r| named.get(r).cloned().unwrap_or_else(|| format!("ref{r}")))
        .collect();
    (index, names)
}

fn convert(file: &MeshbFile, refs: &str) -> Result<MeshData, String> {
    let n = file.vertices.len() as u64;
    let vertex = |v: u64| -> Result<u32, String> {
        if v == 0 || v > n {
            return Err(format!("vertex number {v} out of range 1..={n}"));
        }
        u32::try_from(v - 1).map_err(|_| format!("vertex number {v} beyond 32 bits"))
    };
    let mut tets = Vec::with_capacity(file.tetrahedra.len());
    for (t, _) in &file.tetrahedra {
        tets.push([vertex(t[0])?, vertex(t[1])?, vertex(t[2])?, vertex(t[3])?]);
    }
    let mut triangles = Vec::with_capacity(file.triangles.len());
    for (t, _) in &file.triangles {
        triangles.push([vertex(t[0])?, vertex(t[1])?, vertex(t[2])?]);
    }
    let tet_refs: Vec<i64> = file.tetrahedra.iter().map(|(_, r)| *r).collect();
    let tri_refs: Vec<i64> = file.triangles.iter().map(|(_, r)| *r).collect();
    let (tet_parts, parts) = groups(&tet_refs, &names(refs, "tetrahedra"));
    let (triangle_patches, patches) = groups(&tri_refs, &names(refs, "triangles"));
    Ok(MeshData {
        points: file.vertices.clone(),
        tets,
        tet_parts,
        parts,
        triangles,
        triangle_patches,
        patches,
    })
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;
    use mesher_core::mesh::{Block, SideSet, VolumeMesh};

    /// Two parts of one tet each (sharing a face), a named side set on the
    /// first, and a point no cell uses.
    fn two_parts() -> VolumeMesh {
        let mut mesh = VolumeMesh {
            points: vec![
                [9.0, 9.0, 9.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
            ],
            blocks: vec![
                Block {
                    name: "base part".into(),
                    tets: vec![[1, 2, 3, 4]],
                },
                Block {
                    name: "cap".into(),
                    tets: vec![[2, 5, 3, 4]],
                },
            ],
            side_sets: Vec::new(),
        };
        mesh.add_boundary_side_sets();
        mesh.side_sets.push(SideSet {
            name: "floor".into(),
            sides: vec![(0, 1)],
        });
        mesh
    }

    #[test]
    fn a_written_mesh_reads_back_zero_based_and_named() {
        let mesh = two_parts();
        let written = mesher_io::meshb::meshb(&mesh).unwrap();
        let m = mesh_from_bytes(&written.bytes, &written.refs).unwrap();
        // The unused point is left out, so every index moves down one.
        assert_eq!(m.points, mesh.points[1..].to_vec());
        assert_eq!(m.tets, vec![[0, 1, 2, 3], [1, 4, 2, 3]]);
        assert_eq!(m.tet_parts, vec![0, 1]);
        assert_eq!(m.parts, vec!["base part", "cap"]);
        // Each boundary face once, in its patch; the shared face is internal.
        let f = mesh.faces().unwrap();
        assert_eq!(m.triangles.len(), f.faces.len() - f.internal());
        assert_eq!(m.triangles.len(), 6);
        let names: Vec<&str> = m
            .triangle_patches
            .iter()
            .map(|&p| m.patches[p as usize].as_str())
            .collect();
        assert_eq!(names.iter().filter(|n| **n == "floor").count(), 1);
        assert_eq!(
            m.patches.iter().map(String::as_str).collect::<Vec<_>>(),
            f.patches
                .iter()
                .map(|p| p.name.as_str())
                .collect::<Vec<_>>()
        );
        for t in &m.triangles {
            assert!(t.iter().all(|&v| (v as usize) < m.points.len()));
        }
    }

    #[test]
    fn unnamed_references_are_numbered() {
        let written = mesher_io::meshb::meshb(&two_parts()).unwrap();
        let m = mesh_from_bytes(&written.bytes, "tetrahedra 2 cap\n").unwrap();
        assert_eq!(m.parts, vec!["ref1", "cap"]);
        assert!(m.patches.iter().all(|p| p.starts_with("ref")));
    }

    #[test]
    fn a_broken_file_is_an_error() {
        let written = mesher_io::meshb::meshb(&two_parts()).unwrap();
        let cut = &written.bytes[..written.bytes.len() - 20];
        assert!(mesh_from_bytes(cut, "").unwrap_err().contains("truncated"));
        assert!(mesh_from_bytes(b"not a mesh", "").is_err());
        // A vertex number past the vertices: the last tet's last vertex.
        let mut bad = written.bytes.clone();
        let at = bad.len() - 8 - 8;
        bad[at..at + 4].copy_from_slice(&99u32.to_le_bytes());
        assert!(
            mesh_from_bytes(&bad, "")
                .unwrap_err()
                .contains("out of range")
        );
    }

    #[test]
    fn a_directory_means_its_mesh_meshb() {
        let dir = std::env::temp_dir().join(format!("mesher-py-mesh-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let written = mesher_io::meshb::meshb(&two_parts()).unwrap();
        std::fs::write(dir.join(MESH_MESHB_FILE), &written.bytes).unwrap();
        std::fs::write(dir.join(MESH_MESHB_REFS_FILE), &written.refs).unwrap();
        let from_dir = read_mesh(&dir).unwrap();
        assert_eq!(read_mesh(&dir.join(MESH_MESHB_FILE)).unwrap(), from_dir);
        assert_eq!(from_dir.parts, vec!["base part", "cap"]);
        let missing = read_mesh(&dir.join("none.meshb")).unwrap_err();
        assert!(matches!(missing, ReadError::Io(..)), "{missing}");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
