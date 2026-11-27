//! Milestone 18, stage 2: OBJ groups become side sets on every kernel. Shapes
//! with named faces are meshed conformally, by the exact kernel and by
//! stuffing; every group must come out as one side set of its name whose
//! area is the group's (exactly for a conformal mesh, which keeps the
//! surface; within a tolerance for stuffing, whose boundary approximates it
//! and takes each face's group from the nearest input triangle). A group in
//! two parts is one side set spanning both blocks. (The remeshed path keeps
//! groups too: `group_boundaries_are_kept` in `mesher_tet::remesh` and the
//! remeshed box in `crates/mesher-tet/tests/conforming_shapes.rs`; these flat
//! shapes are never remeshed, as each face is one facet.)
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeMap;

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::mesh::{VolumeMesh, tet_side_nodes};
use mesher_core::repair::repair;
use mesher_core::volume::{MeshOutcome, assemble, mesh_parts};
use mesher_core::{JobConfig, KernelChoice, RecordLog};
use mesher_geom::vec3::{cross, norm, sub};
use mesher_io::parse_obj;
use mesher_synthetic::{Point, Scene, Solid};

/// As in the suite: coarser than the default, to keep it fast.
const ELEMENT_SIZE: f64 = 0.06;

fn mesh(scene: &Scene, config: &JobConfig) -> MeshOutcome {
    let obj = scene.to_obj();
    let mut log = RecordLog::new();
    let model = parse_obj(obj.as_bytes(), &mut log, &mut |_| {});
    let diagnosis = diagnose(&model, &Tolerances::default(), &mut log).unwrap();
    let selected = select_parts(&diagnosis, config, &mut log).unwrap();
    let mut overlay = diagnosis.overlay.clone();
    check_part_pairs(&model, &diagnosis, &selected, &mut log, &mut overlay).unwrap();
    let repaired = repair(
        &model,
        &diagnosis,
        &selected,
        &Tolerances::default(),
        &mut log,
    )
    .unwrap();
    mesh_parts(&model, &repaired, config, &mut log).unwrap()
}

fn area(a: Point, b: Point, c: Point) -> f64 {
    0.5 * norm(cross(sub(b, a), sub(c, a)))
}

/// Each group's area on the input scene.
fn group_areas(scene: &Scene) -> BTreeMap<String, f64> {
    let mut out = BTreeMap::new();
    for (s, (_, solid)) in scene.solids.iter().enumerate() {
        for (t, tri) in solid.tris.iter().enumerate() {
            let [a, b, c] = tri.map(|v| solid.points[v as usize]);
            let name = scene.groups[&s][t].clone();
            *out.entry(name).or_insert(0.0) += area(a, b, c);
        }
    }
    out
}

/// Each side set's area (and how many there are of each name).
fn side_set_areas(mesh: &VolumeMesh) -> BTreeMap<String, (f64, usize)> {
    let tets: Vec<[u32; 4]> = mesh
        .blocks
        .iter()
        .flat_map(|b| b.tets.iter().copied())
        .collect();
    let mut out: BTreeMap<String, (f64, usize)> = BTreeMap::new();
    for set in &mesh.side_sets {
        let a: f64 = set
            .sides
            .iter()
            .map(|&(e, side)| {
                let f = tet_side_nodes(tets[e as usize], side).unwrap();
                let [p, q, r] = f.map(|v| mesh.points[v as usize]);
                area(p, q, r)
            })
            .sum();
        let entry = out.entry(set.name.clone()).or_insert((0.0, 0));
        entry.0 += a;
        entry.1 += 1;
    }
    out
}

/// The unit cube's six faces as groups `z0`, `z1`, `y0`, `y1`, `x0`, `x1`.
fn named_box(lo: Point, hi: Point, names: [&str; 6]) -> (Solid, Vec<String>) {
    let solid = Solid::cuboid(lo, hi);
    // `cuboid` writes its faces z = min, z = max, y = min, y = max, x = min,
    // x = max, two triangles each.
    let groups = names
        .iter()
        .flat_map(|n| [n.to_string(), n.to_string()])
        .collect();
    (solid, groups)
}

/// An `n`-sided can: `bottom` and `top` caps (fans from one vertex) and its
/// `wall`.
fn can(n: u32) -> (Solid, Vec<String>) {
    let polygon: Vec<[f64; 2]> = (0..n)
        .map(|i| {
            let a = 2.0 * std::f64::consts::PI * f64::from(i) / f64::from(n);
            [libm::cos(a), libm::sin(a)]
        })
        .collect();
    let solid = Solid::prism(&polygon, 0.0, 2.0, 0);
    // `prism` writes the caps' fans first, bottom and top alternating, then
    // the wall.
    let caps = 2 * (n as usize - 2);
    let groups = (0..solid.tris.len())
        .map(|t| {
            if t >= caps {
                "wall"
            } else if t % 2 == 0 {
                "bottom"
            } else {
                "top"
            }
            .to_owned()
        })
        .collect();
    (solid, groups)
}

fn configs() -> Vec<(&'static str, JobConfig, f64)> {
    let base = JobConfig {
        element_size: ELEMENT_SIZE,
        ..JobConfig::default()
    };
    vec![
        ("conformal", base.clone(), 1e-9),
        (
            "exact",
            JobConfig {
                kernel: KernelChoice::Exact,
                ..base.clone()
            },
            1e-9,
        ),
        (
            "stuffing",
            JobConfig {
                stuffing_only: true,
                ..base
            },
            0.15,
        ),
    ]
}

/// With `SUITE_FOAM_DIR`, the mesh as an OpenFOAM case, a Fluent file and a
/// `.meshb` file there, as the suite writes them (for `tests/foam/check.sh`
/// and `tests/meshb/meshio_roundtrip.py`).
fn write_cases(name: &str, mesh: &VolumeMesh) {
    let Some(root) = std::env::var_os("SUITE_FOAM_DIR").filter(|d| !d.is_empty()) else {
        return;
    };
    let name: String = name
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    let case = std::path::Path::new(&root).join(name);
    let mut log = RecordLog::new();
    for (file, text) in mesher_io::foam::foam_files(mesh, &mut log).expect("faces") {
        let path = case.join(file);
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, text).unwrap();
    }
    let fluent = mesher_io::fluent::fluent_text(mesh, &mut log).expect("faces");
    std::fs::write(case.join("mesh_fluent.msh"), fluent).unwrap();
    let m = mesher_io::meshb::meshb(mesh).expect("faces");
    std::fs::write(case.join("mesh.meshb"), m.bytes).unwrap();
    std::fs::write(case.join("mesh_refs.txt"), m.refs).unwrap();
}

/// Every group is one side set of its name, of the group's area.
fn check(label: &str, scene: &Scene) {
    let want = group_areas(scene);
    for (kernel, config, tolerance) in configs() {
        let outcome = mesh(scene, &config);
        assert_eq!(
            outcome.meshed.len(),
            scene.solids.len(),
            "{label}, {kernel}: every part meshed"
        );
        let assembled = assemble(&outcome.meshed);
        write_cases(&format!("side_sets_{label}_{kernel}"), &assembled);
        let got = side_set_areas(&assembled);
        // As patches (the face-based formats), each group has its faces
        // exactly, and `<part>_boundary` keeps only what no group has.
        let faces = assembled.faces().expect("faces");
        let patch_area = |name: &str| -> f64 {
            faces
                .patches
                .iter()
                .filter(|p| p.name == name)
                .flat_map(|p| p.start..p.start + p.count)
                .map(|k| {
                    let [a, b, c] = <[u32; 3]>::try_from(faces.faces.get(k))
                        .unwrap()
                        .map(|v| assembled.points[v as usize]);
                    0.5 * norm(cross(sub(b, a), sub(c, a)))
                })
                .sum()
        };
        for (name, &(a, _)) in &got {
            if !name.ends_with("_boundary") {
                let p = patch_area(name);
                assert!(
                    (p - a).abs() <= 1e-12 * a,
                    "{label}, {kernel}: patch {name} {p} against its side set {a}"
                );
            }
        }
        let total: f64 = faces.patches.iter().map(|p| patch_area(&p.name)).sum();
        let boundary: f64 = got
            .iter()
            .filter(|(n, _)| n.ends_with("_boundary"))
            .map(|(_, (a, _))| a)
            .sum();
        assert!(
            (total - boundary).abs() <= 1e-9 * boundary,
            "{label}, {kernel}: patches {total} against the boundary {boundary}"
        );
        for (name, &area) in &want {
            let Some(&(a, sets)) = got.get(name) else {
                panic!("{label}, {kernel}: no side set {name}: {got:?}");
            };
            assert_eq!(sets, 1, "{label}, {kernel}: {sets} side sets named {name}");
            let error = (a - area).abs() / area;
            assert!(
                error <= tolerance,
                "{label}, {kernel}: side set {name} has area {a}, the group {area} ({error:.3e})"
            );
            eprintln!("{label}, {kernel}: {name} area error {error:.2e}");
        }
        for (name, _) in &scene.solids {
            let boundary = format!("{}_boundary", name.as_deref().unwrap());
            assert!(got.contains_key(&boundary), "{label}, {kernel}: {boundary}");
        }
    }
}

#[test]
fn a_box_keeps_its_six_named_faces() {
    let (solid, groups) = named_box(
        [0.0; 3],
        [1.0, 2.0, 0.5],
        ["z0", "z1", "y0", "y1", "x0", "x1"],
    );
    check(
        "box",
        &Scene::default().with_groups(Some("box"), solid, groups),
    );
}

#[test]
fn a_can_keeps_its_caps_and_wall() {
    let (solid, groups) = can(24);
    check(
        "can",
        &Scene::default().with_groups(Some("can"), solid, groups),
    );
}

/// Two boxes, each with its own top and a shared `skin` group: one `skin`
/// side set spans both blocks.
#[test]
fn a_group_in_two_parts_is_one_side_set() {
    let (a, ga) = named_box(
        [0.0; 3],
        [1.0; 3],
        ["skin", "a_top", "skin", "skin", "skin", "skin"],
    );
    let (b, gb) = named_box(
        [2.0, 0.0, 0.0],
        [3.0, 1.0, 1.0],
        ["skin", "b_top", "skin", "skin", "skin", "skin"],
    );
    let scene = Scene::default()
        .with_groups(Some("a"), a, ga)
        .with_groups(Some("b"), b, gb);
    check("two boxes", &scene);
}
