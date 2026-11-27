//! Milestone 21 stage 1: polyhedral face views through the face-based
//! writers. Hand-made polyhedral meshes (a 2 × 2 × 2 grid of hexahedra in
//! two zones, and two stacked hexagonal prisms) are built as `FaceMesh`es
//! of `CellShape::Polyhedron`, checked (layout, every cell closed with its
//! volume), and written as an OpenFOAM case and a Fluent file. With
//! `SUITE_FOAM_DIR=<dir>` the files go there, for `tests/foam/check.sh`
//! (`checkMesh`, and the Fluent file round-tripped through OpenFOAM).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeMap;

use mesher_core::RecordLog;
use mesher_core::mesh::{CellShape, FaceMesh, Patch, Polygons, Zone};
use mesher_geom::Point3;
use mesher_geom::vec3::{add, cross, dot, scale, sub};

/// A face's area vector by Newell's method (right-hand rule).
fn area_vector(points: &[Point3], face: &[u32]) -> Point3 {
    let mut s = [0.0; 3];
    for (i, &a) in face.iter().enumerate() {
        let b = face[(i + 1) % face.len()];
        s = add(s, cross(points[a as usize], points[b as usize]));
    }
    scale(s, 0.5)
}

fn average(points: &[Point3], nodes: &[u32]) -> Point3 {
    let sum = nodes
        .iter()
        .fold([0.0; 3], |s, &v| add(s, points[v as usize]));
    scale(sum, 1.0 / nodes.len() as f64)
}

/// A polyhedral mesh: each cell a list of faces (any orientation; turned
/// out of the cell here, the cells being convex). Faces two cells share
/// are internal, owned by the lower cell; the rest go to the patch
/// `patch(face centre)` names, patches in the order of `patches`.
fn face_mesh(
    points: &[Point3],
    cells: &[Vec<Vec<u32>>],
    zones: &[(&str, usize)],
    patches: &[&str],
    patch: impl Fn(Point3) -> &'static str,
) -> FaceMesh {
    let mut by_key: BTreeMap<Vec<u32>, Vec<(u32, Vec<u32>)>> = BTreeMap::new();
    for (c, faces) in cells.iter().enumerate() {
        let all: Vec<u32> = faces.iter().flatten().copied().collect();
        let centre = average(points, &all);
        for face in faces {
            let mut face = face.clone();
            if dot(
                area_vector(points, &face),
                sub(average(points, &face), centre),
            ) < 0.0
            {
                face.reverse();
            }
            let mut key = face.clone();
            key.sort_unstable();
            by_key.entry(key).or_default().push((c as u32, face));
        }
    }
    let mut internal = Vec::new();
    let mut boundary: Vec<Vec<(u32, Vec<u32>)>> = vec![Vec::new(); patches.len()];
    for (_, mut on) in by_key {
        on.sort_unstable();
        match on.len() {
            1 => {
                let (c, face) = on.pop().unwrap();
                let name = patch(average(points, &face));
                let p = patches.iter().position(|&n| n == name).unwrap();
                boundary[p].push((c, face));
            }
            2 => internal.push((on[0].0, on[1].0, on[0].1.clone())),
            n => panic!("a face on {n} cells"),
        }
    }
    internal.sort_unstable();
    let mut f = FaceMesh {
        faces: Polygons::default(),
        owner: Vec::new(),
        neighbour: Vec::new(),
        patches: Vec::new(),
        zones: Vec::new(),
        cells: cells.len(),
        shape: CellShape::Polyhedron,
    };
    for (o, n, face) in internal {
        f.faces.push(&face);
        f.owner.push(o);
        f.neighbour.push(n);
    }
    for (name, mut faces) in patches.iter().zip(boundary) {
        faces.sort_unstable();
        f.patches.push(Patch {
            name: (*name).to_owned(),
            start: f.faces.len(),
            count: faces.len(),
        });
        for (c, face) in faces {
            f.faces.push(&face);
            f.owner.push(c);
        }
    }
    let mut start = 0;
    for &(name, count) in zones {
        f.zones.push(Zone {
            name: name.to_owned(),
            start,
            count,
        });
        start += count;
    }
    f
}

/// Each cell's volume by the divergence theorem over its faces (each face
/// a fan of triangles about its average point, as a solver computes it),
/// and the largest closure error: the norm of the sum of a cell's outward
/// area vectors, over the largest face area.
fn volumes(points: &[Point3], f: &FaceMesh) -> (Vec<f64>, f64) {
    let mut v = vec![0.0; f.cells];
    let mut closure = vec![[0.0; 3]; f.cells];
    let mut largest: f64 = 0.0;
    for (k, face) in f.faces.iter().enumerate() {
        let c = average(points, face);
        let mut flux = 0.0;
        let mut s = [0.0; 3];
        for (i, &a) in face.iter().enumerate() {
            let b = face[(i + 1) % face.len()];
            let t = scale(
                cross(sub(points[a as usize], c), sub(points[b as usize], c)),
                0.5,
            );
            s = add(s, t);
            let centroid = scale(
                add(add(c, points[a as usize]), points[b as usize]),
                1.0 / 3.0,
            );
            flux += dot(centroid, t) / 3.0;
        }
        largest = largest.max(dot(s, s).sqrt());
        let owner = f.owner[k] as usize;
        v[owner] += flux;
        closure[owner] = add(closure[owner], s);
        if let Some(&n) = f.neighbour.get(k) {
            v[n as usize] -= flux;
            closure[n as usize] = sub(closure[n as usize], s);
        }
    }
    let worst = closure
        .iter()
        .map(|s| dot(*s, *s).sqrt())
        .fold(0.0, f64::max);
    (v, worst / largest)
}

/// A 2 × 2 × 2 grid of unit cubes on [0, 2]³, each a polyhedron of six
/// quadrilaterals; the cells below z = 1 in zone `lower`, the rest in
/// `upper`; patches `bottom` (z = 0) and `walls` (the rest).
fn hex_grid() -> (Vec<Point3>, FaceMesh) {
    let p = |i: u32, j: u32, k: u32| i + 3 * j + 9 * k;
    let mut points = Vec::new();
    for k in 0..3 {
        for j in 0..3 {
            for i in 0..3 {
                points.push([f64::from(i), f64::from(j), f64::from(k)]);
            }
        }
    }
    let mut cells = Vec::new();
    for c in 0..2 {
        for b in 0..2 {
            for a in 0..2 {
                cells.push(vec![
                    vec![
                        p(a, b, c),
                        p(a + 1, b, c),
                        p(a + 1, b + 1, c),
                        p(a, b + 1, c),
                    ],
                    vec![
                        p(a, b, c + 1),
                        p(a + 1, b, c + 1),
                        p(a + 1, b + 1, c + 1),
                        p(a, b + 1, c + 1),
                    ],
                    vec![
                        p(a, b, c),
                        p(a + 1, b, c),
                        p(a + 1, b, c + 1),
                        p(a, b, c + 1),
                    ],
                    vec![
                        p(a, b + 1, c),
                        p(a + 1, b + 1, c),
                        p(a + 1, b + 1, c + 1),
                        p(a, b + 1, c + 1),
                    ],
                    vec![
                        p(a, b, c),
                        p(a, b + 1, c),
                        p(a, b + 1, c + 1),
                        p(a, b, c + 1),
                    ],
                    vec![
                        p(a + 1, b, c),
                        p(a + 1, b + 1, c),
                        p(a + 1, b + 1, c + 1),
                        p(a + 1, b, c + 1),
                    ],
                ]);
            }
        }
    }
    let f = face_mesh(
        &points,
        &cells,
        &[("lower", 4), ("upper", 4)],
        &["bottom", "walls"],
        |c| if c[2] == 0.0 { "bottom" } else { "walls" },
    );
    (points, f)
}

/// Two hexagonal prisms stacked on z ∈ [0, 1] and [1, 2], sharing a
/// hexagon: polygons of six nodes, internal and on the boundary; patches
/// `ends` (z = 0 and 2) and `sides`.
fn hex_prisms() -> (Vec<Point3>, FaceMesh) {
    let mut points = Vec::new();
    for k in 0..3 {
        for i in 0..6 {
            let a = std::f64::consts::PI / 3.0 * f64::from(i);
            points.push([libm::cos(a), libm::sin(a), f64::from(k)]);
        }
    }
    let ring = |k: u32| -> Vec<u32> { (0..6).map(|i| 6 * k + i).collect() };
    let mut cells = Vec::new();
    for k in 0..2 {
        let mut faces = vec![ring(k), ring(k + 1)];
        for i in 0..6 {
            let j = (i + 1) % 6;
            faces.push(vec![6 * k + i, 6 * k + j, 6 * (k + 1) + j, 6 * (k + 1) + i]);
        }
        cells.push(faces);
    }
    let f = face_mesh(&points, &cells, &[("stack", 2)], &["ends", "sides"], |c| {
        if c[2] == 0.0 || c[2] == 2.0 {
            "ends"
        } else {
            "sides"
        }
    });
    (points, f)
}

/// The case's files, written under `SUITE_FOAM_DIR` when it is set.
fn write(name: &str, points: &[Point3], f: &FaceMesh) -> Vec<(&'static str, String)> {
    let mut log = RecordLog::new();
    let files = mesher_io::foam::case_files(points, f, &mut log).unwrap();
    let fluent = mesher_io::fluent::text(points, f, &mut log).unwrap();
    assert!(log.records().is_empty(), "{:?}", log.records());
    if let Some(root) = std::env::var_os("SUITE_FOAM_DIR").filter(|d| !d.is_empty()) {
        let case = std::path::Path::new(&root).join(name);
        for (file, text) in &files {
            let path = case.join(file);
            std::fs::create_dir_all(path.parent().unwrap()).unwrap();
            std::fs::write(&path, text).unwrap();
        }
        std::fs::write(case.join("mesh_fluent.msh"), &fluent).unwrap();
    }
    files
}

fn file<'a>(files: &'a [(&str, String)], name: &str) -> &'a str {
    &files.iter().find(|(n, _)| *n == name).unwrap().1
}

#[test]
fn a_grid_of_hexahedra_as_polyhedra() {
    let (points, f) = hex_grid();
    f.check(points.len()).unwrap();
    // 12 internal faces, 24 on the boundary (4 at the bottom).
    assert_eq!((f.internal(), f.faces.len(), f.cells), (12, 36, 8));
    let counts: Vec<(&str, usize)> = f
        .patches
        .iter()
        .map(|p| (p.name.as_str(), p.count))
        .collect();
    assert_eq!(counts, [("bottom", 4), ("walls", 20)]);
    let (v, closure) = volumes(&points, &f);
    assert!(v.iter().all(|&v| (v - 1.0).abs() < 1e-12), "{v:?}");
    assert!(closure < 1e-12, "{closure}");
    let files = write("polyhedra_hex_grid", &points, &f);
    let faces = file(&files, mesher_io::foam::FACES);
    assert_eq!(
        faces.matches("\n4(").count(),
        36,
        "every face a quadrilateral"
    );
    assert!(file(&files, mesher_io::foam::CELL_ZONES).contains("lower"));
}

#[test]
fn stacked_hexagonal_prisms() {
    let (points, f) = hex_prisms();
    f.check(points.len()).unwrap();
    assert_eq!((f.internal(), f.faces.len(), f.cells), (1, 15, 2));
    let (v, closure) = volumes(&points, &f);
    // A regular hexagon of circumradius 1: area 3√3/2.
    let area = 1.5 * libm::sqrt(3.0);
    assert!(v.iter().all(|&v| (v - area).abs() < 1e-12), "{v:?}");
    assert!(closure < 1e-12, "{closure}");
    let files = write("polyhedra_hex_prisms", &points, &f);
    let faces = file(&files, mesher_io::foam::FACES);
    assert_eq!(faces.matches("\n6(").count(), 3, "the hexagons");
    assert_eq!(faces.matches("\n4(").count(), 12, "the sides");
    // The shared hexagon is internal: out of cell 0, so upwards.
    let shared = f.faces.get(0);
    assert_eq!(shared.len(), 6);
    assert!(area_vector(&points, shared)[2] > 0.0);
}

/// The Fluent file: polyhedral cell zones (element type 7), and face zones
/// of polygons (face type 5) with each line's node count first.
#[test]
fn polyhedra_in_fluent_are_element_type_7_and_polygon_faces() {
    let (points, f) = hex_prisms();
    let mut log = RecordLog::new();
    let text = mesher_io::fluent::text(&points, &f, &mut log).unwrap();
    assert!(text.contains("(12 (2 1 2 1 7))"), "{text}");
    // The interior: one hexagon, its nodes reversed after the first, so the
    // right-hand rule points into its owner, cell 1 (1-based).
    let shared: Vec<String> = f
        .faces
        .get(0)
        .iter()
        .map(|v| format!("{:x}", v + 1))
        .collect();
    let mut reversed = vec![shared[0].clone()];
    reversed.extend(shared[1..].iter().rev().cloned());
    let line = format!("6 {} 1 2\n", reversed.join(" "));
    assert!(
        text.contains(&format!("(13 (3 1 1 2 5)(\n{line}))")),
        "{text}"
    );
    // A face view that breaks the layout is refused before anything is written.
    let mut broken = f;
    broken.zones.clear();
    assert!(mesher_io::fluent::text(&points, &broken, &mut log).is_err());
    assert!(mesher_io::foam::case_files(&points, &broken, &mut log).is_err());
}
