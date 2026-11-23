//! Milestone 18, stage 5: sharp edges and named patches survive meshing
//! (`mesher_core::survival`). On shapes with named faces, meshed
//! conformally, by the exact kernel, and remeshed (a cube with a sliver
//! face, whose surface the pipeline remeshes): every input feature edge is
//! on an output feature edge to rounding, and every group's side set has
//! its area. Stuffing rounds sharp edges off; it is measured, not held to
//! that.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::repair::repair;
use mesher_core::survival::{Survival, survival};
use mesher_core::volume::mesh_parts;
use mesher_core::{JobConfig, KernelChoice, RecordLog};
use mesher_io::parse_obj;
use mesher_synthetic::{Scene, Solid};

/// Mesh a one-part scene and measure its survival.
fn measure(scene: &Scene, config: &JobConfig) -> (Survival, bool) {
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
    let outcome = mesh_parts(&model, &repaired, config, &mut log).unwrap();
    let m = outcome.meshed.first().expect("meshed");
    let (_, summary) = outcome
        .summaries
        .iter()
        .find(|(i, _)| *i == m.part)
        .unwrap();
    let r = outcome
        .wrapped
        .iter()
        .find(|w| w.part == m.part)
        .unwrap_or(&repaired[m.part]);
    (
        survival(&model, r, m, config, summary.element_size),
        summary.remeshed,
    )
}

fn named_box() -> Scene {
    let s = Solid::cuboid([0.0; 3], [1.0, 2.0, 0.5]);
    let names = ["z0", "z1", "y0", "y1", "x0", "x1"];
    let groups = names
        .iter()
        .flat_map(|n| [n.to_string(), n.to_string()])
        .collect();
    Scene::default().with_groups(Some("box"), s, groups)
}

/// An L-shaped prism: a reflex edge as well as convex ones.
fn l_prism() -> Scene {
    let s = Solid::prism(
        &[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ],
        0.0,
        0.6,
        0,
    );
    let n = s.tris.len();
    let groups = (0..n)
        .map(|t| if t < 8 { "caps" } else { "walls" }.to_owned())
        .collect();
    Scene::default().with_groups(Some("ell"), s, groups)
}

/// The unit cube with a sliver face (its z = 0 face split at a point
/// 0.002 off its diagonal): the input's cells are far beyond the limits, so
/// the pipeline remeshes the surface.
fn sliver_cube() -> Scene {
    let mut s = Solid::unit_cube();
    let [a, b, c] = s.tris[0];
    let (pa, pb, pc) = (
        s.points[a as usize],
        s.points[b as usize],
        s.points[c as usize],
    );
    let mid = [0, 1, 2].map(|k| (pa[k] + pb[k]) / 2.0);
    s.points
        .push([0, 1, 2].map(|k| mid[k] + 0.002 * (pc[k] - mid[k])));
    let m = s.points.len() as u32 - 1;
    s.tris[0] = [a, b, m];
    s.tris.push([b, c, m]);
    s.tris.push([c, a, m]);
    // `cuboid`'s faces in order, two triangles each; the two added
    // triangles are on z = 0.
    let names = ["z0", "z1", "y0", "y1", "x0", "x1"];
    let mut groups: Vec<String> = names
        .iter()
        .flat_map(|n| [n.to_string(), n.to_string()])
        .collect();
    groups.push("z0".into());
    groups.push("z0".into());
    Scene::default().with_groups(Some("cube"), s, groups)
}

fn config(kernel: KernelChoice, stuffing_only: bool) -> JobConfig {
    JobConfig {
        element_size: 0.06,
        kernel,
        stuffing_only,
        ..JobConfig::default()
    }
}

/// Every feature edge kept to rounding, every patch's area to 1e-9.
fn assert_kept(label: &str, s: &Survival) {
    assert!(s.input_length > 0.0, "{label}: no features");
    assert!(
        s.max_distance < 1e-9,
        "{label}: a feature edge {} sizes off",
        s.max_distance
    );
    assert!((s.kept - 1.0).abs() < 1e-12, "{label}: kept {}", s.kept);
    for (name, a, b) in &s.patches {
        assert!(
            (a - b).abs() <= 1e-9 * a,
            "{label}: patch {name} {a} against {b}"
        );
    }
}

#[test]
fn conformal_and_exact_meshes_keep_features_and_patches() {
    for (label, scene) in [("box", named_box()), ("ell", l_prism())] {
        for (kernel, name) in [
            (KernelChoice::Conformal, "conformal"),
            (KernelChoice::Exact, "exact"),
        ] {
            let (s, remeshed) = measure(&scene, &config(kernel, false));
            assert!(!remeshed);
            assert_kept(&format!("{label}, {name}"), &s);
            assert_eq!(s.patches.len(), if label == "box" { 6 } else { 2 });
        }
    }
}

#[test]
fn a_remeshed_surface_keeps_features_and_patches() {
    let (s, remeshed) = measure(&sliver_cube(), &config(KernelChoice::Conformal, false));
    assert!(remeshed, "the sliver cube is remeshed");
    assert_kept("sliver cube, remeshed", &s);
    assert_eq!(s.patches.len(), 6);
}

/// Stuffing rounds sharp edges off, within about an element size, and keeps
/// patches to its boundary's approximation.
#[test]
fn stuffing_is_measured() {
    let (s, _) = measure(&named_box(), &config(KernelChoice::Conformal, true));
    eprintln!(
        "box, stuffing: max distance {:.2} sizes, kept {:.2}",
        s.max_distance, s.kept
    );
    assert!(s.max_distance <= 1.0, "{}", s.max_distance);
    for (name, a, b) in &s.patches {
        assert!((a - b).abs() <= 0.15 * a, "patch {name} {a} against {b}");
    }
}
