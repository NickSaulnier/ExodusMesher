//! Milestone 21 stage 4: the polyhedral writers (`cells` `dual`), checked by
//! reading each file back: Exodus II's `NSIDED` face blocks and `NFACED`
//! element blocks with `netcdf3`, VTU's `VTK_POLYHEDRON` cells from the text.
//! Every cell read back is closed, the VTU's faces point out of their cell
//! (the cells' volumes are the mesh's), and each side set holds its patch's
//! faces.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeMap;

use mesher_core::RecordLog;
use mesher_core::dual;
use mesher_core::mesh::{Block, PolyMesh, VolumeMesh};
use mesher_geom::predicates::orient3d;
use mesher_geom::{Point3, Sign};
use mesher_io::exodus::{ExodusOptions, write_exodus_poly};
use mesher_io::vtu;
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

/// Two boxes as two parts, their duals joined as the worker joins them, with
/// the tet mesh's side sets (each block's boundary).
fn sample() -> (PolyMesh, VolumeMesh) {
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
                name: "right".into(),
                tets: right,
            },
        ],
        side_sets: Vec::new(),
    };
    mesh.add_boundary_side_sets();
    let d = dual::dual(&mesh).unwrap();
    let names: Vec<String> = mesh.side_sets.iter().map(|s| s.name.clone()).collect();
    let poly = PolyMesh::join(&[(&d.points, &d.faces)], &names);
    poly.faces.check(poly.points.len()).unwrap();
    (poly, mesh)
}

/// Whether the faces close a cell: every edge on two of them.
fn closed(faces: &[Vec<u32>]) -> bool {
    let mut edges: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    for f in faces {
        for i in 0..f.len() {
            let (a, b) = (f[i], f[(i + 1) % f.len()]);
            *edges.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    edges.values().all(|&n| n == 2)
}

/// The faces' nodes as a sorted list, for comparing faces whatever their
/// orientation and first node.
fn key(face: &[u32]) -> Vec<u32> {
    let mut k = face.to_vec();
    k.sort_unstable();
    k
}

#[test]
fn join_keeps_faces_cells_and_side_sets() {
    let (poly, mesh) = sample();
    let f = &poly.faces;
    assert_eq!(f.zones.len(), 2);
    assert_eq!(f.zones[0].start + f.zones[0].count, f.zones[1].start);
    assert_eq!(poly.side_sets.len(), 2);
    let boundary: usize = f.patches.iter().map(|p| p.count).sum();
    assert_eq!(
        poly.side_sets.iter().map(|s| s.faces.len()).sum::<usize>(),
        boundary
    );
    // Each cell's faces, in face order.
    let cf = f.cell_faces();
    for c in 0..f.cells {
        let list = cf.get(c);
        assert!(list.windows(2).all(|w| w[0] < w[1]));
        let faces: Vec<Vec<u32>> = list
            .iter()
            .map(|&k| f.faces.get(k as usize).to_vec())
            .collect();
        assert!(closed(&faces), "cell {c}");
    }
    assert_eq!(cf.nodes().len(), 2 * f.internal() + boundary);
    // Two parts joined as two: the same mesh.
    let mut parts = Vec::new();
    for b in &mesh.blocks {
        let mut one = VolumeMesh {
            points: mesh.points.clone(),
            blocks: vec![b.clone()],
            side_sets: Vec::new(),
        };
        one.add_boundary_side_sets();
        parts.push(dual::dual(&one).unwrap());
    }
    let names: Vec<String> = mesh.side_sets.iter().map(|s| s.name.clone()).collect();
    let views: Vec<(&[Point3], &mesher_core::mesh::FaceMesh)> = parts
        .iter()
        .map(|d| (d.points.as_slice(), &d.faces))
        .collect();
    let joined = PolyMesh::join(&views, &names);
    joined.faces.check(joined.points.len()).unwrap();
    assert_eq!(joined.faces.cells, f.cells);
    assert_eq!(joined.faces.faces.len(), f.faces.len());
    assert_eq!(joined.faces.internal(), f.internal());
    let counts = |p: &PolyMesh| {
        p.side_sets
            .iter()
            .map(|s| s.faces.len())
            .collect::<Vec<_>>()
    };
    assert_eq!(counts(&joined), counts(&poly));
}

#[test]
fn exodus_nfaced_round_trip() {
    let (poly, _) = sample();
    let f = &poly.faces;
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("dual.exo");
    let mut log = RecordLog::new();
    let options = ExodusOptions {
        title: "dual".into(),
        ..ExodusOptions::default()
    };
    write_exodus_poly(&path, &poly, &options, &mut log).unwrap();
    assert!(log.records().is_empty(), "{:?}", log.records());

    let (ds, version) = FileReader::open(&path).unwrap().close();
    assert_eq!(version, netcdf3::Version::Offset64Bit);
    let mut r = FileReader::open(&path).unwrap();
    assert_eq!(ds.dim_size("num_nodes"), Some(poly.points.len()));
    assert_eq!(ds.dim_size("num_elem"), Some(f.cells));
    assert_eq!(ds.dim_size("num_el_blk"), Some(2));
    assert_eq!(ds.dim_size("num_fa_blk"), Some(2));
    assert_eq!(ds.dim_size("num_face"), Some(f.faces.len()));
    assert_eq!(r.read_var_i32("fa_status").unwrap(), vec![1, 1]);
    assert_eq!(r.read_var_i32("fa_prop1").unwrap(), vec![1, 2]);

    // Faces, numbered globally block by block, and each element's faces.
    let mut faces: Vec<Vec<u32>> = Vec::new();
    let mut elements: Vec<Vec<u32>> = Vec::new();
    for k in 1..=2 {
        assert_eq!(
            ds.get_var_attr_as_string(&format!("fbconn{k}"), "elem_type"),
            Some("NSIDED".into())
        );
        assert_eq!(
            ds.get_var_attr_as_string(&format!("facconn{k}"), "elem_type"),
            Some("NFACED".into())
        );
        assert_eq!(
            ds.get_var_attr_as_string(&format!("ebepecnt{k}"), "entity_type1"),
            Some("FACE".into())
        );
        assert!(!ds.has_var(&format!("connect{k}")));
        let nodes = r.read_var_i32(&format!("fbconn{k}")).unwrap();
        let counts = r.read_var_i32(&format!("fbepecnt{k}")).unwrap();
        assert_eq!(counts.iter().sum::<i32>() as usize, nodes.len());
        let mut at = 0;
        for &n in &counts {
            let face: Vec<u32> = nodes[at..at + n as usize]
                .iter()
                .map(|&v| {
                    assert!(v >= 1 && v as usize <= poly.points.len());
                    v as u32 - 1
                })
                .collect();
            faces.push(face);
            at += n as usize;
        }
        let list = r.read_var_i32(&format!("facconn{k}")).unwrap();
        let counts = r.read_var_i32(&format!("ebepecnt{k}")).unwrap();
        assert_eq!(counts.len(), f.zones[k - 1].count);
        assert_eq!(counts.iter().sum::<i32>() as usize, list.len());
        let mut at = 0;
        for &n in &counts {
            elements.push(
                list[at..at + n as usize]
                    .iter()
                    .map(|&g| g as u32)
                    .collect(),
            );
            at += n as usize;
        }
    }
    assert_eq!(faces.len(), f.faces.len());
    let written: Vec<Vec<u32>> = faces.iter().map(|g| key(g)).collect();
    let mut ours: Vec<Vec<u32>> = f.faces.iter().map(key).collect();
    let mut sorted = written.clone();
    sorted.sort();
    ours.sort();
    assert_eq!(sorted, ours, "the faces, once each");

    // Every element closed; each face on one element (boundary) or two.
    let mut uses = vec![0u32; faces.len()];
    for (e, list) in elements.iter().enumerate() {
        let cell: Vec<Vec<u32>> = list
            .iter()
            .map(|&g| {
                uses[g as usize - 1] += 1;
                faces[g as usize - 1].clone()
            })
            .collect();
        assert!(closed(&cell), "element {}", e + 1);
    }
    assert!(uses.iter().all(|&n| n == 1 || n == 2));
    let boundary = uses.iter().filter(|&&n| n == 1).count();
    assert_eq!(boundary, f.faces.len() - f.internal());

    // Side sets: (element, place in its face list), each a boundary face of
    // the set's patch; node sets their nodes.
    assert_eq!(ds.dim_size("num_side_sets"), Some(2));
    for (i, set) in poly.side_sets.iter().enumerate() {
        let k = i + 1;
        let elems = r.read_var_i32(&format!("elem_ss{k}")).unwrap();
        let sides = r.read_var_i32(&format!("side_ss{k}")).unwrap();
        let mut got: Vec<Vec<u32>> = elems
            .iter()
            .zip(&sides)
            .map(|(&e, &s)| {
                let g = elements[e as usize - 1][s as usize - 1];
                assert_eq!(uses[g as usize - 1], 1, "a side set face on two cells");
                key(&faces[g as usize - 1])
            })
            .collect();
        let mut want: Vec<Vec<u32>> = set
            .faces
            .iter()
            .map(|&k| key(f.faces.get(k as usize)))
            .collect();
        got.sort();
        want.sort();
        assert_eq!(got, want, "side set {}", set.name);
        let mut nodes: Vec<u32> = want.into_iter().flatten().collect();
        nodes.sort_unstable();
        nodes.dedup();
        let ns: Vec<u32> = r
            .read_var_i32(&format!("node_ns{k}"))
            .unwrap()
            .iter()
            .map(|&v| v as u32 - 1)
            .collect();
        assert_eq!(ns, nodes);
    }
}

/// The integers of a VTU data array by name.
fn array(text: &str, name: &str) -> Vec<i64> {
    let tag = format!("Name=\"{name}\" format=\"ascii\">");
    let start = text.find(&tag).unwrap() + tag.len();
    let end = start + text[start..].find("</DataArray>").unwrap();
    text[start..end]
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect()
}

#[test]
fn vtu_polyhedra_point_out() {
    let (poly, mesh) = sample();
    let text = vtu::vtu_poly_text(&poly);
    let n = poly.faces.cells;
    let types = array(&text, "types");
    assert_eq!(types.len(), n);
    assert!(types.iter().all(|&t| t == 42));
    let connectivity = array(&text, "connectivity");
    let offsets = array(&text, "offsets");
    let faces = array(&text, "faces");
    let face_offsets = array(&text, "faceoffsets");
    let blocks = array(&text, "block");
    assert_eq!((offsets.len(), face_offsets.len(), blocks.len()), (n, n, n));
    let (mut start, mut face_start) = (0usize, 0usize);
    let mut volume = 0.0;
    for c in 0..n {
        let points: Vec<u32> = connectivity[start..offsets[c] as usize]
            .iter()
            .map(|&v| v as u32)
            .collect();
        let stream = &faces[face_start..face_offsets[c] as usize];
        let mut cell: Vec<Vec<u32>> = Vec::new();
        let mut at = 1;
        for _ in 0..stream[0] {
            let m = stream[at] as usize;
            cell.push(
                stream[at + 1..at + 1 + m]
                    .iter()
                    .map(|&v| v as u32)
                    .collect(),
            );
            at += 1 + m;
        }
        assert_eq!(at, stream.len());
        assert!(closed(&cell), "cell {c}");
        // Out of the cell: every edge once each way.
        let mut directed: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        for f in &cell {
            for i in 0..f.len() {
                *directed.entry((f[i], f[(i + 1) % f.len()])).or_default() += 1;
            }
        }
        assert!(
            directed
                .keys()
                .all(|&(a, b)| directed.get(&(b, a)) == Some(&1))
        );
        // The cell's points, each once, are its faces'.
        let mut used: Vec<u32> = cell.iter().flatten().copied().collect();
        used.sort_unstable();
        used.dedup();
        let mut listed = points.clone();
        listed.sort_unstable();
        assert_eq!(listed, used);
        // Its volume, by the divergence theorem about the origin, each face
        // a fan about its nodes' average (the same for either orientation).
        for f in &cell {
            let p = |v: u32| poly.points[v as usize];
            let mut a = [0.0; 3];
            for &v in f {
                for (ak, pk) in a.iter_mut().zip(p(v)) {
                    *ak += pk / f.len() as f64;
                }
            }
            for i in 0..f.len() {
                let (b, d) = (p(f[i]), p(f[(i + 1) % f.len()]));
                volume += (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0])
                    + a[2] * (b[0] * d[1] - b[1] * d[0]))
                    / 6.0;
            }
        }
        start = offsets[c] as usize;
        face_start = face_offsets[c] as usize;
    }
    assert!((volume - mesh.volume()).abs() < 1e-12, "{volume}");
    assert_eq!(
        blocks.iter().filter(|&&b| b == 1).count(),
        poly.faces.zones[0].count
    );
}
