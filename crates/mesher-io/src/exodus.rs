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
//! - side sets: `elem_ss<i>` (1-based global element numbers) and `side_ss<i>`
//!   (Exodus TETRA4 side ordinals), with `ss_status`, `ss_prop1`, `ss_names`;
//! - a node set mirroring each side set: `node_ns<i>`, `ns_status`,
//!   `ns_prop1`, `ns_names`;
//! - QA and info records; an unlimited `time_step` dimension with no steps.
//!
//! Empty blocks and sets are written as null entities (status 0, no data), as
//! the Exodus library does. Names longer than 32 characters are truncated and
//! recorded as `EXO_NAME_TRUNCATED`.

use std::path::Path;

use mesher_core::mesh::{VolumeMesh, tet_side_nodes};
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
}

/// Check every variable against the 64-bit-offset limit and every count
/// against 32-bit integers.
pub fn check_capacity(s: Sizes) -> Result<(), ExodusError> {
    let checks = [
        ("coordx", s.nodes, 8),
        ("connect", s.largest_block, 16),
        ("elem_ss", s.largest_side_set, 4),
        ("node_ns", s.largest_node_set, 4),
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

/// Write `mesh` as Exodus II to `path`. Name truncations are recorded in `log`.
pub fn write_exodus(
    path: &Path,
    mesh: &VolumeMesh,
    options: &ExodusOptions,
    log: &mut RecordLog,
) -> Result<(), ExodusError> {
    let n_elem = mesh.element_count() as u64;
    if n_elem == 0 {
        return Err(ExodusError::Empty);
    }
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
    check_capacity(Sizes {
        nodes: mesh.points.len() as u64,
        largest_block: mesh
            .blocks
            .iter()
            .map(|b| b.tets.len() as u64)
            .max()
            .unwrap_or(0),
        largest_side_set: mesh
            .side_sets
            .iter()
            .map(|s| s.sides.len() as u64)
            .max()
            .unwrap_or(0),
        largest_node_set: node_sets.iter().map(|s| s.len() as u64).max().unwrap_or(0),
    })?;

    let block_names: Vec<String> = mesh
        .blocks
        .iter()
        .map(|b| exodus_name("element block", &b.name, log))
        .collect();
    let set_names: Vec<String> = mesh
        .side_sets
        .iter()
        .map(|s| exodus_name("side set", &s.name, log))
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
    def(ds.add_fixed_dim("num_elem", n_elem as usize))?;
    def(ds.add_fixed_dim("num_el_blk", mesh.blocks.len()))?;
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
    def(ds.add_var_f64("coordx", &["num_nodes"]))?;
    def(ds.add_var_f64("coordy", &["num_nodes"]))?;
    def(ds.add_var_f64("coordz", &["num_nodes"]))?;
    def(ds.add_var_u8("coor_names", &["num_dim", "len_name"]))?;
    for (i, block) in mesh.blocks.iter().enumerate() {
        if block.tets.is_empty() {
            continue;
        }
        let k = i + 1;
        def(ds.add_fixed_dim(format!("num_el_in_blk{k}"), block.tets.len()))?;
        def(ds.add_fixed_dim(format!("num_nod_per_el{k}"), 4))?;
        let connect = format!("connect{k}");
        def(ds.add_var_i32(
            &connect,
            &[format!("num_el_in_blk{k}"), format!("num_nod_per_el{k}")],
        ))?;
        def(ds.add_var_attr_string(&connect, "elem_type", "TETRA4"))?;
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
        for (i, set) in mesh.side_sets.iter().enumerate() {
            let k = i + 1;
            if !set.sides.is_empty() {
                let dim = format!("num_side_ss{k}");
                def(ds.add_fixed_dim(&dim, set.sides.len()))?;
                def(ds.add_var_i32(&format!("elem_ss{k}"), &[&dim]))?;
                def(ds.add_var_i32(&format!("side_ss{k}"), &[&dim]))?;
            }
            if !node_sets[i].is_empty() {
                let dim = format!("num_nod_ns{k}");
                def(ds.add_fixed_dim(&dim, node_sets[i].len()))?;
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
    put(w.write_var_i32(
        "eb_status",
        &mesh
            .blocks
            .iter()
            .map(|b| status(b.tets.is_empty()))
            .collect::<Vec<_>>(),
    ))?;
    put(w.write_var_i32("eb_prop1", &i32s(1..=mesh.blocks.len() as u64)))?;
    put(w.write_var_u8("eb_names", &char_rows(&block_names, LEN_NAME)))?;
    for (k, var) in ["coordx", "coordy", "coordz"].iter().enumerate() {
        let c: Vec<f64> = mesh.points.iter().map(|p| p[k]).collect();
        put(w.write_var_f64(var, &c))?;
    }
    put(w.write_var_u8(
        "coor_names",
        &char_rows(&["x".into(), "y".into(), "z".into()], LEN_NAME),
    ))?;
    for (i, block) in mesh.blocks.iter().enumerate() {
        if block.tets.is_empty() {
            continue;
        }
        let conn = i32s(block.tets.iter().flatten().map(|&v| u64::from(v) + 1));
        put(w.write_var_i32(&format!("connect{}", i + 1), &conn))?;
    }
    if n_sets > 0 {
        let ss_status: Vec<i32> = mesh
            .side_sets
            .iter()
            .map(|s| status(s.sides.is_empty()))
            .collect();
        let ns_status: Vec<i32> = node_sets.iter().map(|s| status(s.is_empty())).collect();
        put(w.write_var_i32("ss_status", &ss_status))?;
        put(w.write_var_i32("ns_status", &ns_status))?;
        put(w.write_var_i32("ss_prop1", &i32s(1..=n_sets as u64)))?;
        put(w.write_var_i32("ns_prop1", &i32s(1..=n_sets as u64)))?;
        put(w.write_var_u8("ss_names", &char_rows(&set_names, LEN_NAME)))?;
        put(w.write_var_u8("ns_names", &char_rows(&set_names, LEN_NAME)))?;
        for (i, set) in mesh.side_sets.iter().enumerate() {
            let k = i + 1;
            if !set.sides.is_empty() {
                put(w.write_var_i32(
                    &format!("elem_ss{k}"),
                    &i32s(set.sides.iter().map(|s| u64::from(s.0) + 1)),
                ))?;
                put(w.write_var_i32(
                    &format!("side_ss{k}"),
                    &i32s(set.sides.iter().map(|s| u64::from(s.1))),
                ))?;
            }
            if !node_sets[i].is_empty() {
                put(w.write_var_i32(
                    &format!("node_ns{k}"),
                    &i32s(node_sets[i].iter().map(|&v| u64::from(v) + 1)),
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
