//! Exodus II writer on the netCDF-3 64-bit-offset format (pure Rust, via the
//! `netcdf3` crate).
//!
//! Layout, following the Exodus II specification for a file with 32-bit
//! integers and double-precision reals:
//! - global attributes `api_version`, `version`, `floating_point_word_size`
//!   (8), `file_size` (1), `maximum_name_length` (32), `int64_status` (0), `title`;
//! - `coordx`/`coordy`/`coordz` and `coor_names`;
//! - one element block per part: `connect<i>` (1-based node numbers, attribute
//!   `elem_type = "TETRA4"`), with `eb_status`, `eb_prop1` (IDs) and `eb_names`;
//!   for a polyhedral mesh ([`write_exodus_poly`]), an `NFACED` block per
//!   part (`facconn<i>`, 1-based face numbers, and `ebepecnt<i>`, faces per
//!   element) with an `NSIDED` face block of its faces (`fbconn<i>`, 1-based
//!   node numbers, and `fbepecnt<i>`, nodes per face; `num_face`,
//!   `fa_status`, `fa_prop1`, `fa_names`), as the Exodus library defines
//!   them;
//! - side sets: `elem_ss<i>` (1-based global element numbers) and `side_ss<i>`
//!   (Exodus TETRA4 side ordinals; for polyhedra a face's 1-based place in
//!   its element's list), with `ss_status`, `ss_prop1`, `ss_names`;
//! - a node set mirroring each side set: `node_ns<i>`, `ns_status`,
//!   `ns_prop1`, `ns_names`;
//! - QA and info records; an unlimited `time_step` dimension with no steps.
//!
//! Empty blocks and sets are written as null entities (status 0, no data), as
//! the Exodus library does. Names longer than 32 characters are truncated and
//! recorded as `EXO_NAME_TRUNCATED`.

use std::path::Path;

use mesher_core::mesh::{PolyMesh, VolumeMesh, tet_side_nodes};
use mesher_core::{Action, DefectCode, Outcome, Record, RecordLog, Stage};
use netcdf3::{DataSet, FileWriter, Version};

/// Longest entity name Exodus readers accept by default.
pub const MAX_NAME_LENGTH: usize = 32;
const LEN_NAME: usize = MAX_NAME_LENGTH + 1;
const LEN_STRING: usize = 33;
const LEN_LINE: usize = 81;
/// Largest variable the 64-bit-offset format can hold.
pub const MAX_VARIABLE_BYTES: u64 = (1 << 32) - 4;

/// A QA record: program, version, date, time (each at most 32 characters).
/// Leave date and time empty to keep output byte-identical across runs.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct QaRecord {
    pub program: String,
    pub version: String,
    pub date: String,
    pub time: String,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ExodusOptions {
    pub title: String,
    pub qa: Vec<QaRecord>,
    /// Free-text lines (at most 80 characters each; longer ones are wrapped).
    pub info: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExodusError {
    /// A variable would exceed the format's limit (`EXO_CAPACITY_EXCEEDED`).
    Capacity { variable: String, bytes: u64 },
    /// The mesh has no elements; there is nothing to write.
    Empty,
    /// The file could not be defined or written (`OUT_WRITE_FAILED`).
    Write(String),
}

impl std::fmt::Display for ExodusError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExodusError::Capacity { variable, bytes } => write!(
                f,
                "variable {variable} would need {bytes} bytes, over the 64-bit-offset limit"
            ),
            ExodusError::Empty => write!(f, "the mesh has no elements"),
            ExodusError::Write(e) => write!(f, "{e}"),
        }
    }
}

fn write_err(e: impl std::fmt::Debug) -> ExodusError {
    ExodusError::Write(format!("{e:?}"))
}

/// Sizes that decide whether a mesh fits; checked before anything is written.
#[derive(Debug, Clone, Copy)]
pub struct Sizes {
    pub nodes: u64,
    /// Elements per block.
    pub largest_block: u64,
    pub largest_side_set: u64,
    pub largest_node_set: u64,
    /// Polyhedral meshes: the most nodes of one `NSIDED` face block's faces.
    pub largest_face_block_nodes: u64,
    /// Polyhedral meshes: the most faces listed by one `NFACED` block.
    pub largest_element_faces: u64,
}

/// Check every variable against the 64-bit-offset limit and every count
/// against 32-bit integers.
pub fn check_capacity(s: Sizes) -> Result<(), ExodusError> {
    let checks = [
        ("coordx", s.nodes, 8),
        ("connect", s.largest_block, 16),
        ("elem_ss", s.largest_side_set, 4),
        ("node_ns", s.largest_node_set, 4),
        ("fbconn", s.largest_face_block_nodes, 4),
        ("facconn", s.largest_element_faces, 4),
    ];
    for (variable, count, bytes_each) in checks {
        let bytes = count.saturating_mul(bytes_each);
        if bytes > MAX_VARIABLE_BYTES || count > i32::MAX as u64 {
            return Err(ExodusError::Capacity {
                variable: variable.to_owned(),
                bytes,
            });
        }
    }
    Ok(())
}

/// Truncate a name to [`MAX_NAME_LENGTH`] bytes at a character boundary,
/// recording a warning when it is cut.
fn exodus_name(kind: &str, name: &str, log: &mut RecordLog) -> String {
    if name.len() <= MAX_NAME_LENGTH {
        return name.to_owned();
    }
    let mut cut = MAX_NAME_LENGTH;
    while !name.is_char_boundary(cut) {
        cut -= 1;
    }
    let short = name[..cut].to_owned();
    log.push(
        Record::new(DefectCode::ExoNameTruncated, Stage::Emit)
            .measured("entity", kind)
            .measured("name", name)
            .measured("written_as", short.clone())
            .threshold("max_name_length", MAX_NAME_LENGTH as u64)
            .action(Action::Truncate, Outcome::Applied),
    );
    short
}

/// Fixed-width, NUL-padded character rows.
fn char_rows(rows: &[String], width: usize) -> Vec<u8> {
    let mut out = vec![0u8; rows.len() * width];
    for (i, row) in rows.iter().enumerate() {
        let bytes = row.as_bytes();
        let n = bytes.len().min(width - 1);
        out[i * width..i * width + n].copy_from_slice(&bytes[..n]);
    }
    out
}

fn i32s(values: impl IntoIterator<Item = u64>) -> Vec<i32> {
    values
        .into_iter()
        .map(|v| i32::try_from(v).unwrap_or(i32::MAX))
        .collect()
}

/// An element block for [`write_model`].
enum BlockData<'a> {
    /// `TETRA4` connectivity, 0-based nodes.
    Tets(&'a [[u32; 4]]),
    /// An `NFACED` block with its `NSIDED` face block: the faces' nodes
    /// (0-based) and counts, each element's faces (1-based global face
    /// numbers) and counts.
    Polyhedra {
        face_nodes: Vec<u32>,
        face_counts: Vec<u32>,
        element_faces: Vec<u32>,
        element_counts: Vec<u32>,
    },
}

impl BlockData<'_> {
    fn len(&self) -> usize {
        match self {
            BlockData::Tets(t) => t.len(),
            BlockData::Polyhedra { element_counts, .. } => element_counts.len(),
        }
    }
}

/// What [`write_model`] writes: points, blocks by name, side sets as
/// `(0-based global element, 1-based side)` by name, and a node set per side
/// set.
struct Model<'a> {
    points: &'a [mesher_geom::Point3],
    blocks: Vec<(&'a str, BlockData<'a>)>,
    side_sets: Vec<(&'a str, Vec<(u32, u32)>)>,
    node_sets: Vec<Vec<u32>>,
}

/// Write `mesh` as Exodus II to `path`. Name truncations are recorded in `log`.
pub fn write_exodus(
    path: &Path,
    mesh: &VolumeMesh,
    options: &ExodusOptions,
    log: &mut RecordLog,
) -> Result<(), ExodusError> {
    // Node sets mirror side sets: the distinct nodes of each set's faces.
    let elements: Vec<[u32; 4]> = mesh.elements().collect();
    let node_sets: Vec<Vec<u32>> = mesh
        .side_sets
        .iter()
        .map(|set| {
            let mut nodes: Vec<u32> = set
                .sides
                .iter()
                .filter_map(|&(e, s)| elements.get(e as usize).and_then(|&t| tet_side_nodes(t, s)))
                .flatten()
                .collect();
            nodes.sort_unstable();
            nodes.dedup();
            nodes
        })
        .collect();
    let model = Model {
        points: &mesh.points,
        blocks: mesh
            .blocks
            .iter()
            .map(|b| (b.name.as_str(), BlockData::Tets(&b.tets)))
            .collect(),
        side_sets: mesh
            .side_sets
            .iter()
            .map(|s| {
                (
                    s.name.as_str(),
                    s.sides.iter().map(|&(e, k)| (e, u32::from(k))).collect(),
                )
            })
            .collect(),
        node_sets,
    };
    write_model(path, &model, options, log)
}

/// Write a polyhedral mesh (`cells` `dual`) as Exodus II to `path`: per
/// zone, an `NSIDED` face block (the faces its cells own, in face order;
/// faces are numbered globally, block by block) and an `NFACED` element
/// block (each cell's faces in face order); side sets as `(element, side)`,
/// the side a face's 1-based place in its owner's list. Faces are written as
/// the face view orients them, out of their owner, so a reader orients each
/// for its neighbour. Name truncations are recorded in `log`.
pub fn write_exodus_poly(
    path: &Path,
    mesh: &PolyMesh,
    options: &ExodusOptions,
    log: &mut RecordLog,
) -> Result<(), ExodusError> {
    let f = &mesh.faces;
    let cell_faces = f.cell_faces();
    // Each zone's faces (those its cells own) and every face's global number.
    let mut block_faces: Vec<Vec<u32>> = vec![Vec::new(); f.zones.len()];
    let mut zone = vec![0usize; f.cells];
    for (z, zn) in f.zones.iter().enumerate() {
        for c in zn.start..zn.start + zn.count {
            if let Some(slot) = zone.get_mut(c) {
                *slot = z;
            }
        }
    }
    for (k, &o) in f.owner.iter().enumerate() {
        let z = zone.get(o as usize).copied().unwrap_or(0);
        if let Some(list) = block_faces.get_mut(z) {
            list.push(k as u32);
        }
    }
    let mut global = vec![0u32; f.faces.len()];
    let mut next = 1u32;
    for list in &block_faces {
        for &k in list {
            global[k as usize] = next;
            next += 1;
        }
    }
    let blocks: Vec<(&str, BlockData<'_>)> = f
        .zones
        .iter()
        .zip(&block_faces)
        .map(|(zn, faces)| {
            let mut face_nodes = Vec::new();
            let mut face_counts = Vec::with_capacity(faces.len());
            for &k in faces {
                let nodes = f.faces.get(k as usize);
                face_nodes.extend_from_slice(nodes);
                face_counts.push(nodes.len() as u32);
            }
            let mut element_faces = Vec::new();
            let mut element_counts = Vec::with_capacity(zn.count);
            for c in zn.start..zn.start + zn.count {
                let list = cell_faces.get(c);
                element_faces.extend(list.iter().map(|&k| global[k as usize]));
                element_counts.push(list.len() as u32);
            }
            (
                zn.name.as_str(),
                BlockData::Polyhedra {
                    face_nodes,
                    face_counts,
                    element_faces,
                    element_counts,
                },
            )
        })
        .collect();
    let side_sets: Vec<(&str, Vec<(u32, u32)>)> = mesh
        .side_sets
        .iter()
        .map(|set| {
            let sides = set
                .faces
                .iter()
                .filter_map(|&k| {
                    let owner = *f.owner.get(k as usize)?;
                    let place = cell_faces
                        .get(owner as usize)
                        .iter()
                        .position(|&g| g == k)?;
                    Some((owner, place as u32 + 1))
                })
                .collect();
            (set.name.as_str(), sides)
        })
        .collect();
    let node_sets: Vec<Vec<u32>> = mesh
        .side_sets
        .iter()
        .map(|set| {
            let mut nodes: Vec<u32> = set
                .faces
                .iter()
                .filter(|&&k| (k as usize) < f.faces.len())
                .flat_map(|&k| f.faces.get(k as usize).iter().copied())
                .collect();
            nodes.sort_unstable();
            nodes.dedup();
            nodes
        })
        .collect();
    let model = Model {
        points: &mesh.points,
        blocks,
        side_sets,
        node_sets,
    };
    write_model(path, &model, options, log)
}

fn write_model(
    path: &Path,
    mesh: &Model<'_>,
    options: &ExodusOptions,
    log: &mut RecordLog,
) -> Result<(), ExodusError> {
    let n_elem: usize = mesh.blocks.iter().map(|(_, b)| b.len()).sum();
    if n_elem == 0 {
        return Err(ExodusError::Empty);
    }
    let polyhedra = mesh
        .blocks
        .iter()
        .any(|(_, b)| matches!(b, BlockData::Polyhedra { .. }));
    let (mut face_block_nodes, mut element_face_list, mut n_faces) = (0u64, 0u64, 0usize);
    for (_, b) in &mesh.blocks {
        if let BlockData::Polyhedra {
            face_nodes,
            face_counts,
            element_faces,
            ..
        } = b
        {
            face_block_nodes = face_block_nodes.max(face_nodes.len() as u64);
            element_face_list = element_face_list.max(element_faces.len() as u64);
            n_faces += face_counts.len();
        }
    }
    check_capacity(Sizes {
        nodes: mesh.points.len() as u64,
        largest_block: mesh
            .blocks
            .iter()
            .filter(|(_, b)| matches!(b, BlockData::Tets(_)))
            .map(|(_, b)| b.len() as u64)
            .max()
            .unwrap_or(0),
        largest_side_set: mesh
            .side_sets
            .iter()
            .map(|(_, s)| s.len() as u64)
            .max()
            .unwrap_or(0),
        largest_node_set: mesh
            .node_sets
            .iter()
            .map(|s| s.len() as u64)
            .max()
            .unwrap_or(0),
        largest_face_block_nodes: face_block_nodes,
        largest_element_faces: element_face_list,
    })?;

    let block_names: Vec<String> = mesh
        .blocks
        .iter()
        .map(|(name, _)| exodus_name("element block", name, log))
        .collect();
    let set_names: Vec<String> = mesh
        .side_sets
        .iter()
        .map(|(name, _)| exodus_name("side set", name, log))
        .collect();
    let info: Vec<String> = options
        .info
        .iter()
        .flat_map(|line| {
            let chars: Vec<char> = line.chars().collect();
            chars
                .chunks(LEN_LINE - 1)
                .map(|c| c.iter().collect::<String>())
                .collect::<Vec<_>>()
        })
        .collect();

    let mut ds = DataSet::new();
    let def = |r: Result<(), netcdf3::InvalidDataSet>| r.map_err(write_err);
    def(ds.add_global_attr_f32("api_version", vec![8.03]))?;
    def(ds.add_global_attr_f32("version", vec![8.03]))?;
    def(ds.add_global_attr_i32("floating_point_word_size", vec![8]))?;
    def(ds.add_global_attr_i32("file_size", vec![1]))?;
    def(ds.add_global_attr_i32("maximum_name_length", vec![MAX_NAME_LENGTH as i32]))?;
    def(ds.add_global_attr_i32("int64_status", vec![0]))?;
    def(ds.add_global_attr_string("title", &options.title))?;

    def(ds.add_fixed_dim("len_string", LEN_STRING))?;
    def(ds.add_fixed_dim("len_line", LEN_LINE))?;
    def(ds.add_fixed_dim("four", 4))?;
    def(ds.add_fixed_dim("len_name", LEN_NAME))?;
    def(ds.set_unlimited_dim("time_step", 0))?;
    def(ds.add_fixed_dim("num_dim", 3))?;
    def(ds.add_fixed_dim("num_nodes", mesh.points.len()))?;
    def(ds.add_fixed_dim("num_elem", n_elem))?;
    def(ds.add_fixed_dim("num_el_blk", mesh.blocks.len()))?;
    if polyhedra {
        def(ds.add_fixed_dim("num_face", n_faces))?;
        def(ds.add_fixed_dim("num_fa_blk", mesh.blocks.len()))?;
    }
    let n_sets = mesh.side_sets.len();
    if n_sets > 0 {
        def(ds.add_fixed_dim("num_side_sets", n_sets))?;
        def(ds.add_fixed_dim("num_node_sets", n_sets))?;
    }
    if !options.qa.is_empty() {
        def(ds.add_fixed_dim("num_qa_rec", options.qa.len()))?;
    }
    if !info.is_empty() {
        def(ds.add_fixed_dim("num_info", info.len()))?;
    }

    def(ds.add_var_f64("time_whole", &["time_step"]))?;
    def(ds.add_var_i32("eb_status", &["num_el_blk"]))?;
    def(ds.add_var_i32("eb_prop1", &["num_el_blk"]))?;
    def(ds.add_var_attr_string("eb_prop1", "name", "ID"))?;
    def(ds.add_var_u8("eb_names", &["num_el_blk", "len_name"]))?;
    if polyhedra {
        def(ds.add_var_i32("fa_status", &["num_fa_blk"]))?;
        def(ds.add_var_i32("fa_prop1", &["num_fa_blk"]))?;
        def(ds.add_var_attr_string("fa_prop1", "name", "ID"))?;
        def(ds.add_var_u8("fa_names", &["num_fa_blk", "len_name"]))?;
    }
    def(ds.add_var_f64("coordx", &["num_nodes"]))?;
    def(ds.add_var_f64("coordy", &["num_nodes"]))?;
    def(ds.add_var_f64("coordz", &["num_nodes"]))?;
    def(ds.add_var_u8("coor_names", &["num_dim", "len_name"]))?;
    for (i, (_, block)) in mesh.blocks.iter().enumerate() {
        if block.len() == 0 {
            continue;
        }
        let k = i + 1;
        match block {
            BlockData::Tets(tets) => {
                def(ds.add_fixed_dim(format!("num_el_in_blk{k}"), tets.len()))?;
                def(ds.add_fixed_dim(format!("num_nod_per_el{k}"), 4))?;
                let connect = format!("connect{k}");
                def(ds.add_var_i32(
                    &connect,
                    &[format!("num_el_in_blk{k}"), format!("num_nod_per_el{k}")],
                ))?;
                def(ds.add_var_attr_string(&connect, "elem_type", "TETRA4"))?;
            }
            BlockData::Polyhedra {
                face_nodes,
                face_counts,
                element_faces,
                element_counts,
            } => {
                // The face block, as the Exodus library defines an `NSIDED`
                // block: all its nodes in one list, a count per face.
                def(ds.add_fixed_dim(format!("num_fa_in_blk{k}"), face_counts.len()))?;
                def(ds.add_fixed_dim(format!("num_nod_per_fa{k}"), face_nodes.len()))?;
                let fbconn = format!("fbconn{k}");
                def(ds.add_var_i32(&fbconn, &[format!("num_nod_per_fa{k}")]))?;
                def(ds.add_var_attr_string(&fbconn, "elem_type", "NSIDED"))?;
                let fbepecnt = format!("fbepecnt{k}");
                def(ds.add_var_i32(&fbepecnt, &[format!("num_fa_in_blk{k}")]))?;
                def(ds.add_var_attr_string(&fbepecnt, "entity_type1", "NODE"))?;
                def(ds.add_var_attr_string(&fbepecnt, "entity_type2", "FACE"))?;
                // The element block: all its faces in one list, a count per
                // element.
                def(ds.add_fixed_dim(format!("num_el_in_blk{k}"), element_counts.len()))?;
                def(ds.add_fixed_dim(format!("num_fac_per_el{k}"), element_faces.len()))?;
                let facconn = format!("facconn{k}");
                def(ds.add_var_i32(&facconn, &[format!("num_fac_per_el{k}")]))?;
                def(ds.add_var_attr_string(&facconn, "elem_type", "NFACED"))?;
                let ebepecnt = format!("ebepecnt{k}");
                def(ds.add_var_i32(&ebepecnt, &[format!("num_el_in_blk{k}")]))?;
                def(ds.add_var_attr_string(&ebepecnt, "entity_type1", "FACE"))?;
                def(ds.add_var_attr_string(&ebepecnt, "entity_type2", "ELEM"))?;
            }
        }
    }
    if n_sets > 0 {
        for (prefix, count_dim) in [("ss", "num_side_sets"), ("ns", "num_node_sets")] {
            let (status, prop, names) = (
                format!("{prefix}_status"),
                format!("{prefix}_prop1"),
                format!("{prefix}_names"),
            );
            def(ds.add_var_i32(&status, &[count_dim]))?;
            def(ds.add_var_i32(&prop, &[count_dim]))?;
            def(ds.add_var_attr_string(&prop, "name", "ID"))?;
            def(ds.add_var_u8(&names, &[count_dim, "len_name"]))?;
        }
        for (i, (_, sides)) in mesh.side_sets.iter().enumerate() {
            let k = i + 1;
            if !sides.is_empty() {
                let dim = format!("num_side_ss{k}");
                def(ds.add_fixed_dim(&dim, sides.len()))?;
                def(ds.add_var_i32(&format!("elem_ss{k}"), &[&dim]))?;
                def(ds.add_var_i32(&format!("side_ss{k}"), &[&dim]))?;
            }
            if !mesh.node_sets[i].is_empty() {
                let dim = format!("num_nod_ns{k}");
                def(ds.add_fixed_dim(&dim, mesh.node_sets[i].len()))?;
                def(ds.add_var_i32(&format!("node_ns{k}"), &[&dim]))?;
            }
        }
    }
    if !options.qa.is_empty() {
        def(ds.add_var_u8("qa_records", &["num_qa_rec", "four", "len_string"]))?;
    }
    if !info.is_empty() {
        def(ds.add_var_u8("info_records", &["num_info", "len_line"]))?;
    }

    let mut w = FileWriter::open(path).map_err(write_err)?;
    w.set_def(&ds, Version::Offset64Bit, 0).map_err(write_err)?;
    let put = |r: Result<(), netcdf3::WriteError>| r.map_err(write_err);

    let status = |empty: bool| i32::from(!empty);
    let block_status: Vec<i32> = mesh
        .blocks
        .iter()
        .map(|(_, b)| status(b.len() == 0))
        .collect();
    put(w.write_var_i32("eb_status", &block_status))?;
    put(w.write_var_i32("eb_prop1", &i32s(1..=mesh.blocks.len() as u64)))?;
    put(w.write_var_u8("eb_names", &char_rows(&block_names, LEN_NAME)))?;
    if polyhedra {
        put(w.write_var_i32("fa_status", &block_status))?;
        put(w.write_var_i32("fa_prop1", &i32s(1..=mesh.blocks.len() as u64)))?;
        put(w.write_var_u8("fa_names", &char_rows(&block_names, LEN_NAME)))?;
    }
    for (k, var) in ["coordx", "coordy", "coordz"].iter().enumerate() {
        let c: Vec<f64> = mesh.points.iter().map(|p| p[k]).collect();
        put(w.write_var_f64(var, &c))?;
    }
    put(w.write_var_u8(
        "coor_names",
        &char_rows(&["x".into(), "y".into(), "z".into()], LEN_NAME),
    ))?;
    for (i, (_, block)) in mesh.blocks.iter().enumerate() {
        if block.len() == 0 {
            continue;
        }
        let k = i + 1;
        match block {
            BlockData::Tets(tets) => {
                let conn = i32s(tets.iter().flatten().map(|&v| u64::from(v) + 1));
                put(w.write_var_i32(&format!("connect{k}"), &conn))?;
            }
            BlockData::Polyhedra {
                face_nodes,
                face_counts,
                element_faces,
                element_counts,
            } => {
                let nodes = i32s(face_nodes.iter().map(|&v| u64::from(v) + 1));
                put(w.write_var_i32(&format!("fbconn{k}"), &nodes))?;
                let counts = i32s(face_counts.iter().map(|&n| u64::from(n)));
                put(w.write_var_i32(&format!("fbepecnt{k}"), &counts))?;
                let faces = i32s(element_faces.iter().map(|&n| u64::from(n)));
                put(w.write_var_i32(&format!("facconn{k}"), &faces))?;
                let counts = i32s(element_counts.iter().map(|&n| u64::from(n)));
                put(w.write_var_i32(&format!("ebepecnt{k}"), &counts))?;
            }
        }
    }
    if n_sets > 0 {
        let ss_status: Vec<i32> = mesh
            .side_sets
            .iter()
            .map(|(_, s)| status(s.is_empty()))
            .collect();
        let ns_status: Vec<i32> = mesh
            .node_sets
            .iter()
            .map(|s| status(s.is_empty()))
            .collect();
        put(w.write_var_i32("ss_status", &ss_status))?;
        put(w.write_var_i32("ns_status", &ns_status))?;
        put(w.write_var_i32("ss_prop1", &i32s(1..=n_sets as u64)))?;
        put(w.write_var_i32("ns_prop1", &i32s(1..=n_sets as u64)))?;
        put(w.write_var_u8("ss_names", &char_rows(&set_names, LEN_NAME)))?;
        put(w.write_var_u8("ns_names", &char_rows(&set_names, LEN_NAME)))?;
        for (i, (_, sides)) in mesh.side_sets.iter().enumerate() {
            let k = i + 1;
            if !sides.is_empty() {
                put(w.write_var_i32(
                    &format!("elem_ss{k}"),
                    &i32s(sides.iter().map(|s| u64::from(s.0) + 1)),
                ))?;
                put(w.write_var_i32(
                    &format!("side_ss{k}"),
                    &i32s(sides.iter().map(|s| u64::from(s.1))),
                ))?;
            }
            if !mesh.node_sets[i].is_empty() {
                put(w.write_var_i32(
                    &format!("node_ns{k}"),
                    &i32s(mesh.node_sets[i].iter().map(|&v| u64::from(v) + 1)),
                ))?;
            }
        }
    }
    if !options.qa.is_empty() {
        let fields: Vec<String> = options
            .qa
            .iter()
            .flat_map(|q| {
                [
                    q.program.clone(),
                    q.version.clone(),
                    q.date.clone(),
                    q.time.clone(),
                ]
            })
            .collect();
        put(w.write_var_u8("qa_records", &char_rows(&fields, LEN_STRING)))?;
    }
    if !info.is_empty() {
        put(w.write_var_u8("info_records", &char_rows(&info, LEN_LINE)))?;
    }
    w.close().map_err(write_err)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capacity_limits() {
        let ok = Sizes {
            nodes: 1_000_000,
            largest_block: 5_000_000,
            largest_side_set: 100_000,
            largest_node_set: 50_000,
            largest_face_block_nodes: 0,
            largest_element_faces: 0,
        };
        assert!(check_capacity(ok).is_ok());
        // 300 million tets need 4.8 GB of connectivity.
        let big = Sizes {
            largest_block: 300_000_000,
            ..ok
        };
        assert_eq!(
            check_capacity(big),
            Err(ExodusError::Capacity {
                variable: "connect".into(),
                bytes: 4_800_000_000
            })
        );
        let many_nodes = Sizes {
            nodes: 600_000_000,
            ..ok
        };
        assert!(matches!(
            check_capacity(many_nodes),
            Err(ExodusError::Capacity { .. })
        ));
    }

    #[test]
    fn names_are_truncated_at_character_boundaries() {
        let mut log = RecordLog::new();
        assert_eq!(exodus_name("block", "short", &mut log), "short");
        assert!(log.records().is_empty());
        let long = "é".repeat(20); // 40 bytes
        let short = exodus_name("block", &long, &mut log);
        assert_eq!(short.len(), 32);
        assert_eq!(log.records()[0].code, DefectCode::ExoNameTruncated);
    }

    #[test]
    fn char_rows_pad_and_cut() {
        let rows = char_rows(&["ab".into(), "x".repeat(40)], 5);
        assert_eq!(rows, b"ab\0\0\0xxxx\0".to_vec());
    }
}
