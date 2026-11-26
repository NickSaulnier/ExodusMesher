//! libMeshb binary files (`.meshb`, the Gamma Mesh Format), from the format's
//! public description: write a volume mesh, read a surface.
//!
//! Layout: a header of two 32-bit integers, `1` (so a reader can tell the
//! byte order) and the version; then keywords, each its 32-bit code, the
//! byte position of the next keyword, and its data; the last is `End` (54),
//! whose next position is 0. The version sets the field sizes:
//!
//! | version | reals | integers (counts, indices, references) | positions |
//! | --- | --- | --- | --- |
//! | 1 | 32-bit | 32-bit | 32-bit |
//! | 2 | 64-bit | 32-bit | 32-bit |
//! | 3 | 64-bit | 32-bit | 64-bit |
//! | 4 | 64-bit | 64-bit | 64-bit |
//!
//! `Dimension` (3) holds one integer; `Vertices` (4) a count, then per
//! vertex its coordinates and a reference; `Triangles` (6) and `Tetrahedra`
//! (8) a count, then per element its 1-based vertex numbers and a reference.
//!
//! The writer writes little-endian, version 2 when every position fits 32
//! bits, else 3 (4 if a count did not fit 32 bits): `Dimension` 3, the
//! points cells use (reference 0), the tetrahedra (reference: the block's
//! 1-based number) and the boundary triangles out of the mesh (reference:
//! the [`FaceMesh`] patch's 1-based number). The format has no names, so
//! [`Meshb::refs`] lists what each reference is (`mesh_refs.txt`).
//!
//! The reader takes any version and byte order. As input, a file's
//! triangles are the surface, their references its groups (`ref<k>`), which
//! are patches only: parts are the connected components. Problems become
//! records (`PARSE_FORMAT_INVALID` for a file it cannot follow, keeping what
//! was read before; `PARSE_BAD_LINE` for a vertex that is not finite;
//! `PARSE_BAD_INDEX` for a triangle naming a missing vertex;
//! `PARSE_UNSUPPORTED_ELEMENT` for other keywords, tetrahedra included).

use std::fmt::Write as _;

use mesher_core::mesh::{FaceMesh, MeshProblem, VolumeMesh};
use mesher_core::model::{Polygon, SurfaceModel, Triangle};
use mesher_core::{Action, DefectCode, Outcome, Record, RecordLog, Severity, Stage};

use crate::shared::used_points;

/// Keyword codes.
pub const KW_DIMENSION: u32 = 3;
pub const KW_VERTICES: u32 = 4;
pub const KW_TRIANGLES: u32 = 6;
pub const KW_TETRAHEDRA: u32 = 8;
pub const KW_END: u32 = 54;

/// What the writer produces: the file, and the meaning of its references.
#[derive(Debug, Clone, PartialEq)]
pub struct Meshb {
    pub bytes: Vec<u8>,
    /// One line per reference: `tetrahedra <ref> <block>` and
    /// `triangles <ref> <patch>`.
    pub refs: String,
    pub version: u32,
}

/// Field sizes of a version.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Sizes {
    real: usize,
    int: usize,
    pos: usize,
}

const V1: Sizes = Sizes {
    real: 4,
    int: 4,
    pos: 4,
};
const V2: Sizes = Sizes {
    real: 8,
    int: 4,
    pos: 4,
};
const V3: Sizes = Sizes {
    real: 8,
    int: 4,
    pos: 8,
};
const V4: Sizes = Sizes {
    real: 8,
    int: 8,
    pos: 8,
};

fn sizes(version: u32) -> Option<Sizes> {
    match version {
        1 => Some(V1),
        2 => Some(V2),
        3 => Some(V3),
        4 => Some(V4),
        _ => None,
    }
}

struct Writer {
    out: Vec<u8>,
    s: Sizes,
}

impl Writer {
    fn int(&mut self, v: u64) {
        match self.s.int {
            4 => self.out.extend_from_slice(&(v as u32).to_le_bytes()),
            _ => self.out.extend_from_slice(&v.to_le_bytes()),
        }
    }
    fn real(&mut self, v: f64) {
        self.out.extend_from_slice(&v.to_le_bytes());
    }
    /// A keyword's code and the position of the next, `size` bytes of data
    /// after this header.
    fn keyword(&mut self, code: u32, size: usize) {
        self.out.extend_from_slice(&code.to_le_bytes());
        let next = if code == KW_END {
            0
        } else {
            (self.out.len() + self.s.pos + size) as u64
        };
        match self.s.pos {
            4 => self.out.extend_from_slice(&(next as u32).to_le_bytes()),
            _ => self.out.extend_from_slice(&next.to_le_bytes()),
        }
    }
}

/// Render `mesh` as a `.meshb` file. Deterministic.
///
/// # Errors
/// The mesh's [`MeshProblem`] when it has no face view (a face on more than
/// two cells: a bug upstream, as the volume gate rejects such meshes).
pub fn meshb(mesh: &VolumeMesh) -> Result<Meshb, MeshProblem> {
    let f = mesh.faces()?;
    Ok(render(mesh, &f))
}

fn render(mesh: &VolumeMesh, f: &FaceMesh) -> Meshb {
    let (new, kept) = used_points(mesh.points.len(), f);
    let boundary = f.faces.len() - f.internal();
    let tets = mesh.element_count();
    let largest = kept.len().max(boundary).max(tets) as u64;
    let size = |s: Sizes| {
        let header = 4 + s.pos;
        8 + (header + s.int)
            + (header + s.int + kept.len() * (3 * s.real + s.int))
            + (header + s.int + boundary * 4 * s.int)
            + (header + s.int + tets * 5 * s.int)
            + header
    };
    let (version, s): (u32, Sizes) = if largest > i32::MAX as u64 {
        (4, V4)
    } else if size(V2) as u64 <= i32::MAX as u64 {
        (2, V2)
    } else {
        (3, V3)
    };
    let mut w = Writer {
        out: Vec::with_capacity(size(s)),
        s,
    };
    w.out.extend_from_slice(&1u32.to_le_bytes());
    w.out.extend_from_slice(&version.to_le_bytes());
    w.keyword(KW_DIMENSION, s.int);
    w.int(3);

    w.keyword(KW_VERTICES, s.int + kept.len() * (3 * s.real + s.int));
    w.int(kept.len() as u64);
    for &v in &kept {
        for c in mesh.points[v] {
            w.real(c);
        }
        w.int(0);
    }
    let node = |v: u32| u64::from(new.get(v as usize).copied().unwrap_or(0)) + 1;

    w.keyword(KW_TRIANGLES, s.int + boundary * 4 * s.int);
    w.int(boundary as u64);
    for (r, p) in f.patches.iter().enumerate() {
        for face in &f.faces[p.start..p.start + p.count] {
            for &v in face {
                w.int(node(v));
            }
            w.int(r as u64 + 1);
        }
    }

    w.keyword(KW_TETRAHEDRA, s.int + tets * 5 * s.int);
    w.int(tets as u64);
    for (r, b) in mesh.blocks.iter().enumerate() {
        for t in &b.tets {
            for &v in t {
                w.int(node(v));
            }
            w.int(r as u64 + 1);
        }
    }
    w.keyword(KW_END, 0);

    let mut refs = String::new();
    for (r, b) in mesh.blocks.iter().enumerate() {
        let _ = writeln!(refs, "tetrahedra {} {}", r + 1, b.name);
    }
    for (r, p) in f.patches.iter().enumerate() {
        let _ = writeln!(refs, "triangles {} {}", r + 1, p.name);
    }
    Meshb {
        bytes: w.out,
        refs,
        version,
    }
}

/// A `.meshb` file's contents, as read: vertex numbers 1-based as written.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MeshbFile {
    pub version: u32,
    pub dimension: u32,
    pub vertices: Vec<[f64; 3]>,
    pub vertex_refs: Vec<i64>,
    pub triangles: Vec<([u64; 3], i64)>,
    pub tetrahedra: Vec<([u64; 4], i64)>,
    /// Keywords read past: code and, where it has one, count.
    pub skipped: Vec<(u32, u64)>,
}

/// Why a file could not be followed to its end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MeshbError {
    /// Byte offset where reading stopped.
    pub offset: u64,
    pub reason: String,
}

struct Reader<'a> {
    bytes: &'a [u8],
    at: usize,
    big_endian: bool,
}

impl Reader<'_> {
    fn take<const N: usize>(&mut self) -> Option<[u8; N]> {
        let end = self.at.checked_add(N)?;
        let b: [u8; N] = self.bytes.get(self.at..end)?.try_into().ok()?;
        self.at = end;
        Some(b)
    }
    fn u32(&mut self) -> Option<u32> {
        let b = self.take::<4>()?;
        Some(if self.big_endian {
            u32::from_be_bytes(b)
        } else {
            u32::from_le_bytes(b)
        })
    }
    fn u64(&mut self) -> Option<u64> {
        let b = self.take::<8>()?;
        Some(if self.big_endian {
            u64::from_be_bytes(b)
        } else {
            u64::from_le_bytes(b)
        })
    }
    /// A signed integer of `size` bytes.
    fn int(&mut self, size: usize) -> Option<i64> {
        if size == 4 {
            self.u32().map(|v| i64::from(v as i32))
        } else {
            self.u64().map(|v| v as i64)
        }
    }
    fn real(&mut self, size: usize) -> Option<f64> {
        if size == 4 {
            self.u32().map(|v| f64::from(f32::from_bits(v)))
        } else {
            self.u64().map(f64::from_bits)
        }
    }
}

/// Read a `.meshb` file: everything read before any problem, and the
/// problem, if there was one.
#[must_use]
pub fn read_meshb(bytes: &[u8]) -> (MeshbFile, Option<MeshbError>) {
    let mut file = MeshbFile::default();
    let error = read_into(bytes, &mut file).err();
    (file, error)
}

fn read_into(bytes: &[u8], file: &mut MeshbFile) -> Result<(), MeshbError> {
    let fail = |offset: usize, reason: &str| MeshbError {
        offset: offset as u64,
        reason: reason.to_owned(),
    };
    let mut r = Reader {
        bytes,
        at: 0,
        big_endian: false,
    };
    match bytes.get(0..4) {
        Some([1, 0, 0, 0]) => {}
        Some([0, 0, 0, 1]) => r.big_endian = true,
        _ => return Err(fail(0, "not a libMeshb binary file (no byte-order code 1)")),
    }
    r.at = 4;
    let version = r.u32().ok_or_else(|| fail(4, "truncated header"))?;
    let s = sizes(version).ok_or_else(|| fail(4, &format!("unknown version {version}")))?;
    file.version = version;
    let truncated = |r: &Reader<'_>| fail(r.at, "truncated");
    loop {
        let start = r.at;
        let code = r.u32().ok_or_else(|| fail(start, "missing End keyword"))?;
        let next = match s.pos {
            4 => r.u32().map(u64::from),
            _ => r.u64(),
        }
        .ok_or_else(|| truncated(&r))?;
        if code == KW_END {
            return Ok(());
        }
        let count = |r: &mut Reader<'_>| -> Result<usize, MeshbError> {
            let n = r.int(s.int).ok_or_else(|| truncated(r))?;
            usize::try_from(n).map_err(|_| fail(r.at, "negative count"))
        };
        match code {
            KW_DIMENSION => {
                let d = r.int(s.int).ok_or_else(|| truncated(&r))?;
                if d != 3 {
                    return Err(fail(start, &format!("dimension {d}: only 3 is read")));
                }
                file.dimension = 3;
            }
            KW_VERTICES => {
                if file.dimension != 3 {
                    return Err(fail(start, "vertices before a 3D dimension"));
                }
                let n = count(&mut r)?;
                if n.saturating_mul(3 * s.real + s.int) > bytes.len() - r.at.min(bytes.len()) {
                    return Err(truncated(&r));
                }
                for _ in 0..n {
                    let mut p = [0.0; 3];
                    for c in &mut p {
                        *c = r.real(s.real).ok_or_else(|| truncated(&r))?;
                    }
                    file.vertices.push(p);
                    file.vertex_refs
                        .push(r.int(s.int).ok_or_else(|| truncated(&r))?);
                }
            }
            KW_TRIANGLES | KW_TETRAHEDRA => {
                let k = if code == KW_TRIANGLES { 3 } else { 4 };
                let n = count(&mut r)?;
                if n.saturating_mul((k + 1) * s.int) > bytes.len() - r.at.min(bytes.len()) {
                    return Err(truncated(&r));
                }
                for _ in 0..n {
                    let mut v = [0u64; 4];
                    for x in v.iter_mut().take(k) {
                        // Negative numbers become 0, never a vertex.
                        *x = u64::try_from(r.int(s.int).ok_or_else(|| truncated(&r))?).unwrap_or(0);
                    }
                    let reference = r.int(s.int).ok_or_else(|| truncated(&r))?;
                    if k == 3 {
                        file.triangles.push(([v[0], v[1], v[2]], reference));
                    } else {
                        file.tetrahedra.push((v, reference));
                    }
                }
            }
            _ => {
                // Most keywords start with a count; it is only reported.
                let n = r.int(s.int).unwrap_or(0).max(0) as u64;
                file.skipped.push((code, n));
            }
        }
        let next = usize::try_from(next).map_err(|_| fail(start, "bad next position"))?;
        if next <= start || next > bytes.len() {
            return Err(fail(
                start,
                &format!("keyword {code}: bad next position {next}"),
            ));
        }
        r.at = next;
    }
}

/// The name of a keyword the input reader does not mesh.
fn keyword_name(code: u32) -> String {
    match code {
        5 => "edges".to_owned(),
        7 => "quadrilaterals".to_owned(),
        KW_TETRAHEDRA => "tetrahedra".to_owned(),
        9 => "prisms".to_owned(),
        10 => "hexahedra".to_owned(),
        c => format!("keyword {c}"),
    }
}

/// Parse `.meshb` bytes as an input surface. Never panics; every problem
/// becomes a record in `log`.
pub fn parse_meshb(bytes: &[u8], log: &mut RecordLog) -> SurfaceModel {
    let (file, error) = read_meshb(bytes);
    if let Some(e) = &error {
        log.push(
            Record::new(DefectCode::ParseFormatInvalid, Stage::Parse)
                .measured("offset", e.offset)
                .measured("reason", e.reason.clone())
                .action(Action::Abort, Outcome::Applied)
                .hint("reading stopped here; what came before is kept"),
        );
    }
    let mut model = SurfaceModel {
        bytes: bytes.len() as u64,
        groups_are_patches: true,
        ..SurfaceModel::default()
    };
    let mut usable = Vec::with_capacity(file.vertices.len());
    for (i, p) in file.vertices.iter().enumerate() {
        let ok = p.iter().all(|c| c.is_finite());
        if !ok {
            log.push(
                Record::new(DefectCode::ParseBadLine, Stage::Parse)
                    .measured("vertex", i as u64 + 1)
                    .action(Action::SkipLine, Outcome::Applied)
                    .hint("a vertex coordinate is not finite"),
            );
        }
        usable.push(ok);
        model.vertices.push(if ok { *p } else { [0.0; 3] });
        model.vertex_lines.push(0);
    }
    let mut group_of: std::collections::BTreeMap<i64, u32> = std::collections::BTreeMap::new();
    for (k, (v, reference)) in file.triangles.iter().enumerate() {
        let index = |n: u64| {
            usize::try_from(n)
                .ok()
                .and_then(|n| n.checked_sub(1))
                .filter(|&i| usable.get(i).copied().unwrap_or(false))
        };
        let (Some(a), Some(b), Some(c)) = (index(v[0]), index(v[1]), index(v[2])) else {
            let bad = v.iter().copied().find(|&n| index(n).is_none()).unwrap_or(0);
            log.push(
                Record::new(DefectCode::ParseBadIndex, Stage::Parse)
                    .measured("triangle", k as u64 + 1)
                    .measured("index", bad)
                    .measured("vertices_defined", file.vertices.len() as u64)
                    .action(Action::DropFace, Outcome::Applied)
                    .hint("triangle names a vertex that is missing or not finite"),
            );
            continue;
        };
        let next = group_of.len() as u32;
        let group = *group_of.entry(*reference).or_insert_with(|| {
            model.groups.push(format!("ref{reference}"));
            next
        });
        let polygon = model.polygons.len() as u32;
        model.polygons.push(Polygon {
            vertices: vec![a as u32, b as u32, c as u32],
            line: 0,
            object: None,
            group: Some(group),
        });
        model.triangles.push(Triangle {
            vertices: [a as u32, b as u32, c as u32],
            polygon,
        });
    }
    let mut skipped = file.skipped.clone();
    if !file.tetrahedra.is_empty() {
        skipped.push((KW_TETRAHEDRA, file.tetrahedra.len() as u64));
    }
    for (code, count) in skipped {
        log.push(
            Record::new(DefectCode::ParseUnsupportedElement, Stage::Parse)
                .measured("statement", keyword_name(code))
                .measured("count", count)
                .action(Action::SkipLine, Outcome::Applied)
                .hint("only triangles are read from a .meshb file"),
        );
    }
    if model.triangles.is_empty() {
        log.push(
            Record::new(DefectCode::ParseNoGeometry, Stage::Parse)
                .severity(Severity::Fatal)
                .measured("vertices", model.vertices.len() as u64)
                .measured("polygons", 0u64)
                .action(Action::Abort, Outcome::Applied)
                .hint("the file contains no usable triangles"),
        );
    }
    model
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;
    use mesher_core::mesh::{Block, SideSet};

    fn two_tets() -> VolumeMesh {
        let mut mesh = VolumeMesh {
            points: vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
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
    fn two_tets_round_trip() {
        let mesh = two_tets();
        let m = meshb(&mesh).unwrap();
        assert_eq!(m.version, 2);
        // Header 8; dimension 12; vertices 12 + 5 x 28; triangles 12 + 6 x
        // 16; tetrahedra 12 + 2 x 20; End 8.
        assert_eq!(m.bytes.len(), 8 + 12 + 152 + 108 + 52 + 8);
        assert_eq!(&m.bytes[..8], &[1, 0, 0, 0, 2, 0, 0, 0]);
        assert_eq!(
            m.refs,
            "tetrahedra 1 b\ntriangles 1 b_boundary\ntriangles 2 top\n"
        );
        let (file, error) = read_meshb(&m.bytes);
        assert_eq!(error, None);
        assert_eq!(file.dimension, 3);
        // The unused point is left out.
        assert_eq!(file.vertices, mesh.points[..5].to_vec());
        assert_eq!(file.tetrahedra, vec![([1, 2, 3, 4], 1), ([2, 5, 3, 4], 1)]);
        let f = mesh.faces().unwrap();
        let want: Vec<([u64; 3], i64)> = (f.internal()..f.faces.len())
            .map(|k| {
                let r = if k < 4 { 1 } else { 2 };
                (f.faces[k].map(|v| u64::from(v) + 1), r)
            })
            .collect();
        assert_eq!(file.triangles, want);
    }

    /// Hand-written files of the other versions and byte order.
    fn file(version: u32, big: bool) -> Vec<u8> {
        let s = sizes(version).unwrap();
        let mut out: Vec<u8> = Vec::new();
        let put = |out: &mut Vec<u8>, v: u64, n: usize| {
            let b = v.to_be_bytes();
            let l = v.to_le_bytes();
            if big {
                out.extend_from_slice(&b[8 - n..]);
            } else {
                out.extend_from_slice(&l[..n]);
            }
        };
        put(&mut out, 1, 4);
        put(&mut out, u64::from(version), 4);
        let kw = |out: &mut Vec<u8>, code: u32, size: usize| {
            put(out, u64::from(code), 4);
            let next = if code == KW_END {
                0
            } else {
                (out.len() + s.pos + size) as u64
            };
            put(out, next, s.pos);
        };
        kw(&mut out, KW_DIMENSION, s.int);
        put(&mut out, 3, s.int);
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0f64],
        ];
        kw(&mut out, KW_VERTICES, s.int + 4 * (3 * s.real + s.int));
        put(&mut out, 4, s.int);
        for p in points {
            for c in p {
                if s.real == 4 {
                    put(&mut out, u64::from((c as f32).to_bits()), 4);
                } else {
                    put(&mut out, c.to_bits(), 8);
                }
            }
            put(&mut out, 0, s.int);
        }
        // An unknown keyword (Corners, 13) to skip.
        kw(&mut out, 13, 2 * s.int);
        put(&mut out, 1, s.int);
        put(&mut out, 1, s.int);
        let tris = [[1, 3, 2], [1, 2, 4], [1, 4, 3], [2, 3, 4u64]];
        kw(&mut out, KW_TRIANGLES, s.int + 4 * 4 * s.int);
        put(&mut out, 4, s.int);
        for (i, t) in tris.iter().enumerate() {
            for &v in t {
                put(&mut out, v, s.int);
            }
            put(&mut out, if i == 0 { 7 } else { 3 }, s.int);
        }
        kw(&mut out, KW_END, 0);
        out
    }

    #[test]
    fn every_version_and_byte_order_reads_as_a_surface() {
        for version in 1..=4 {
            for big in [false, true] {
                let mut log = RecordLog::new();
                let m = parse_meshb(&file(version, big), &mut log);
                assert_eq!(m.vertices.len(), 4, "v{version} big {big}");
                assert_eq!(m.vertices[3], [0.0, 0.0, 1.0]);
                assert_eq!(m.triangles.len(), 4);
                assert_eq!(m.triangles[0].vertices, [0, 2, 1]);
                assert_eq!(m.groups, ["ref7", "ref3"]);
                assert!(m.groups_are_patches);
                let codes: Vec<DefectCode> = log.records().iter().map(|r| r.code).collect();
                assert_eq!(codes, [DefectCode::ParseUnsupportedElement], "{codes:?}");
            }
        }
    }

    #[test]
    fn a_broken_file_keeps_what_came_before() {
        let good = file(2, false);
        // Cut inside the triangles: the vertices survive, no triangles do.
        let mut log = RecordLog::new();
        let m = parse_meshb(&good[..good.len() - 20], &mut log);
        assert_eq!(m.vertices.len(), 4);
        let codes: Vec<DefectCode> = log.records().iter().map(|r| r.code).collect();
        assert!(codes.contains(&DefectCode::ParseFormatInvalid), "{codes:?}");
        assert!(codes.contains(&DefectCode::ParseNoGeometry), "{codes:?}");
        // Not a meshb file at all, an unknown version, garbage: records, no panic.
        for bytes in [
            &b"v 0 0 0\n"[..],
            &[1, 0, 0, 0, 9, 0, 0, 0][..],
            &[1, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0xff, 0xff, 0xff, 0x7f][..],
            &[][..],
        ] {
            let mut log = RecordLog::new();
            let m = parse_meshb(bytes, &mut log);
            assert!(m.triangles.is_empty());
            assert!(log.has_fatal(), "{bytes:?}");
        }
        // A triangle naming a missing vertex is dropped.
        let mut bad = good.clone();
        let at = bad.len() - 8 - 16 * 4 + 4;
        bad[at..at + 4].copy_from_slice(&99u32.to_le_bytes());
        let mut log = RecordLog::new();
        let m = parse_meshb(&bad, &mut log);
        assert_eq!(m.triangles.len(), 3);
        assert!(
            log.records()
                .iter()
                .any(|r| r.code == DefectCode::ParseBadIndex)
        );
    }
}
