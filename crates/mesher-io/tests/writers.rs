//! Volume-mesh writers, checked by reading every file back with an
//! independent reader: `netcdf3` for Exodus II, `mshio` for MSH 4.1, `vtkio`
//! for VTU. Side-set ordinals are verified per face against geometry.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::path::Path;

use mesher_core::mesh::{Block, SideSet, VolumeMesh, tet_side_nodes};
use mesher_core::{DefectCode, RecordLog};
use mesher_geom::predicates::orient3d;
use mesher_geom::{Point3, Sign};
use mesher_io::exodus::{ExodusOptions, QaRecord, write_exodus};
use mesher_io::{msh, vtu};
use netcdf3::FileReader;

/// Kuhn triangulation of the box `[x0, x0 + 1] × [0, 1]²`: six positive tets.
fn kuhn_box(points: &mut Vec<Point3>, x0: f64) -> Vec<[u32; 4]> {
    let base = points.len() as u32;
    for i in 0..8 {
        points.push([
            x0 + f64::from(i & 1),
            f64::from((i >> 1) & 1),
            f64::from((i >> 2) & 1),
        ]);
    }
    let mut tets = Vec::new();
    for [a, b] in [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]] {
        let mut t = [base, base + (1 << a), base + (1 << a) + (1 << b), base + 7];
        let p = t.map(|v| points[v as usize]);
        if orient3d(p[0], p[1], p[2], p[3]) == Sign::Negative {
            t.swap(0, 1);
        }
        tets.push(t);
    }
    tets
}

/// Two unit boxes as two blocks (separate nodes), boundary side sets, and an
/// empty third block.
fn sample() -> VolumeMesh {
    let mut points = Vec::new();
    let left = kuhn_box(&mut points, 0.0);
    let right = kuhn_box(&mut points, 3.0);
    let mut mesh = VolumeMesh {
        points,
        blocks: vec![
            Block {
                name: "left".into(),
                tets: left,
            },
            Block {
                name: "a_part_name_that_is_far_too_long_for_exodus".into(),
                tets: right,
            },
            Block {
                name: "empty".into(),
                tets: Vec::new(),
            },
        ],
        side_sets: Vec::new(),
    };
    mesh.add_boundary_side_sets();
    assert!(mesh.validate(10).is_empty());
    mesh
}

fn chars(bytes: &[u8], width: usize) -> Vec<String> {
    bytes
        .chunks(width)
        .map(|row| {
            let end = row.iter().position(|&b| b == 0).unwrap_or(row.len());
            String::from_utf8_lossy(&row[..end]).into_owned()
        })
        .collect()
}

fn write(mesh: &VolumeMesh, path: &Path) -> RecordLog {
    let mut log = RecordLog::new();
    let options = ExodusOptions {
        title: "sample".into(),
        qa: vec![QaRecord {
            program: "ExodusMesher".into(),
            version: "0.1.0".into(),
            date: String::new(),
            time: String::new(),
        }],
        info: vec!["first line".into(), "x".repeat(100)],
    };
    write_exodus(path, mesh, &options, &mut log).unwrap();
    log
}

#[test]
fn exodus_round_trip() {
    let mesh = sample();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("sample.exo");
    let log = write(&mesh, &path);

    // The long block name was truncated, and recorded.
    assert_eq!(log.records().len(), 2, "{:?}", log.records());
    assert!(
        log.records()
            .iter()
            .all(|r| r.code == DefectCode::ExoNameTruncated)
    );

    let (ds, version) = FileReader::open(&path).unwrap().close();
    assert_eq!(version, netcdf3::Version::Offset64Bit);
    let mut r = FileReader::open(&path).unwrap();
    let gi = |name: &str| {
        ds.get_global_attr(name)
            .and_then(|a| a.get_i32())
            .map(<[i32]>::to_vec)
    };
    assert_eq!(gi("floating_point_word_size"), Some(vec![8]));
    assert_eq!(gi("file_size"), Some(vec![1]));
    assert_eq!(gi("maximum_name_length"), Some(vec![32]));
    assert_eq!(gi("int64_status"), Some(vec![0]));
    assert_eq!(ds.dim_size("num_nodes"), Some(16));
    assert_eq!(ds.dim_size("num_elem"), Some(12));
    assert_eq!(ds.dim_size("num_el_blk"), Some(3));
    assert_eq!(ds.dim_size("num_side_sets"), Some(2));
    assert!(ds.has_unlimited_dim());

    // Coordinates exactly.
    for (k, var) in ["coordx", "coordy", "coordz"].iter().enumerate() {
        let c = r.read_var_f64(var).unwrap();
        let expected: Vec<f64> = mesh.points.iter().map(|p| p[k]).collect();
        assert_eq!(c, expected);
    }
    // Blocks: status, ids, names, 1-based connectivity, element type.
    assert_eq!(r.read_var_i32("eb_status").unwrap(), vec![1, 1, 0]);
    assert_eq!(r.read_var_i32("eb_prop1").unwrap(), vec![1, 2, 3]);
    let names = chars(&r.read_var_u8("eb_names").unwrap(), 33);
    assert_eq!(names[0], "left");
    assert_eq!(names[1], "a_part_name_that_is_far_too_long");
    assert_eq!(names[2], "empty");
    assert_eq!(
        ds.get_var_attr_as_string("connect1", "elem_type")
            .as_deref(),
        Some("TETRA4")
    );
    assert!(!ds.has_var("connect3"));
    let mut read_elements: Vec<[u32; 4]> = Vec::new();
    for k in 1..=2 {
        let conn = r.read_var_i32(&format!("connect{k}")).unwrap();
        assert!(conn.iter().all(|v| (1..=16).contains(v)));
        read_elements.extend(
            conn.chunks(4)
                .map(|c| [0, 1, 2, 3].map(|i| (c[i] - 1) as u32)),
        );
    }
    assert_eq!(read_elements, mesh.elements().collect::<Vec<_>>());

    // Side sets: every side, rebuilt from the read-back connectivity, is a
    // face on the block's boundary whose normal points out of the block.
    let block_of = |e: usize| if e < 6 { 0.0 } else { 3.0 };
    for k in 1..=2 {
        let elems = r.read_var_i32(&format!("elem_ss{k}")).unwrap();
        let sides = r.read_var_i32(&format!("side_ss{k}")).unwrap();
        assert_eq!(elems.len(), 12, "a unit box has 12 boundary triangles");
        let mut nodes_seen = Vec::new();
        for (&e, &s) in elems.iter().zip(&sides) {
            assert!((1..=4).contains(&s));
            let e = (e - 1) as usize;
            let face = tet_side_nodes(read_elements[e], s as u8).unwrap();
            nodes_seen.extend(face);
            let [a, b, c] = face.map(|v| mesh.points[v as usize]);
            let centre = [block_of(e) + 0.5, 0.5, 0.5];
            // Outward: the box centre is behind the face.
            assert_eq!(
                orient3d(a, b, c, centre),
                Sign::Negative,
                "set {k} element {e} side {s}"
            );
            // On the boundary: all three nodes share an extreme coordinate.
            let on_wall = (0..3).any(|axis| {
                let lo = if axis == 0 { block_of(e) } else { 0.0 };
                [a, b, c].iter().all(|p| p[axis] == lo)
                    || [a, b, c].iter().all(|p| p[axis] == lo + 1.0)
            });
            assert!(
                on_wall,
                "set {k} element {e} side {s} is not on the box surface"
            );
        }
        nodes_seen.sort_unstable();
        nodes_seen.dedup();
        let ns: Vec<u32> = r
            .read_var_i32(&format!("node_ns{k}"))
            .unwrap()
            .iter()
            .map(|&v| (v - 1) as u32)
            .collect();
        assert_eq!(ns, nodes_seen, "node set {k} mirrors side set {k}");
    }
    let set_names = chars(&r.read_var_u8("ss_names").unwrap(), 33);
    assert_eq!(set_names[0], "left_boundary");
    assert_eq!(set_names[1].len(), 32);

    // QA and info records.
    let qa = chars(&r.read_var_u8("qa_records").unwrap(), 33);
    assert_eq!(qa, vec!["ExodusMesher", "0.1.0", "", ""]);
    let info = chars(&r.read_var_u8("info_records").unwrap(), 81);
    assert_eq!(info.len(), 3);
    assert_eq!(info[1].len(), 80);
    assert_eq!(r.read_var_f64("time_whole").unwrap(), Vec::<f64>::new());
}

#[test]
fn exodus_output_is_byte_identical() {
    let mesh = sample();
    let dir = tempfile::tempdir().unwrap();
    let (a, b) = (dir.path().join("a.exo"), dir.path().join("b.exo"));
    write(&mesh, &a);
    write(&mesh, &b);
    assert_eq!(std::fs::read(a).unwrap(), std::fs::read(b).unwrap());
}

#[test]
fn exodus_refuses_an_empty_mesh() {
    let dir = tempfile::tempdir().unwrap();
    let mesh = VolumeMesh {
        points: vec![[0.0; 3]],
        blocks: vec![Block {
            name: "e".into(),
            tets: Vec::new(),
        }],
        side_sets: vec![SideSet {
            name: "s".into(),
            sides: Vec::new(),
        }],
    };
    let err = write_exodus(
        &dir.path().join("e.exo"),
        &mesh,
        &ExodusOptions::default(),
        &mut RecordLog::new(),
    );
    assert_eq!(err, Err(mesher_io::exodus::ExodusError::Empty));
}

#[test]
fn msh_round_trip() {
    let mesh = sample();
    let text = msh::msh_text(&mesh);
    assert_eq!(text, msh::msh_text(&mesh));
    let parsed = mshio::parse_msh_bytes(text.as_bytes()).unwrap();
    let data = parsed.data;
    let nodes = data.nodes.unwrap();
    let points: Vec<Point3> = nodes.node_blocks[0]
        .nodes
        .iter()
        .map(|n| [n.x, n.y, n.z])
        .collect();
    assert_eq!(points, mesh.points);
    let entities = data.entities.unwrap();
    assert_eq!(entities.volumes.len(), 3);
    assert_eq!(entities.surfaces.len(), 2);
    let elements = data.elements.unwrap();
    let mut tets = Vec::new();
    let mut triangles = 0;
    for block in &elements.element_blocks {
        for e in &block.elements {
            match block.entity_dim {
                3 => tets.push([0, 1, 2, 3].map(|i| (e.nodes[i] - 1) as u32)),
                2 => {
                    triangles += 1;
                    assert_eq!(e.nodes.len(), 3);
                }
                d => panic!("unexpected entity dimension {d}"),
            }
        }
    }
    assert_eq!(tets, mesh.elements().collect::<Vec<_>>());
    assert_eq!(triangles, 24);
    assert!(text.contains("3 1 \"left\""));
    assert!(text.contains("2 1 \"left_boundary\""));
}

#[test]
fn vtu_round_trip() {
    let mesh = sample();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("sample.vtu");
    std::fs::write(&path, vtu::vtu_text(&mesh)).unwrap();
    let vtk = vtkio::Vtk::import(&path).unwrap();
    let vtkio::model::DataSet::UnstructuredGrid { pieces, .. } = vtk.data else {
        panic!("not an unstructured grid");
    };
    let piece = pieces[0].load_piece_data(None).unwrap();
    let points: Vec<f64> = piece.points.into_vec::<f64>().unwrap();
    let expected: Vec<f64> = mesh.points.iter().flatten().copied().collect();
    assert_eq!(points, expected);
    let (n, legacy) = piece.cells.cell_verts.into_legacy();
    assert_eq!(n, 12);
    let tets: Vec<[u32; 4]> = legacy
        .chunks(5)
        .map(|c| {
            assert_eq!(c[0], 4);
            [c[1], c[2], c[3], c[4]]
        })
        .collect();
    assert_eq!(tets, mesh.elements().collect::<Vec<_>>());
    assert!(
        piece
            .cells
            .types
            .iter()
            .all(|t| *t == vtkio::model::CellType::Tetra)
    );
}
