//! Milestone 9 exit criterion, headless: select parts of a multi-object OBJ
//! the way the GUI does (part tree and viewport picking), export through the
//! worker, and get an Exodus file with one block per selected part; then
//! replay the saved project with `mesher run project.toml` and get the same.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::AtomicBool;
use std::time::Duration;

use mesher_core::project::Project;
use mesher_core::protocol::MESH_EXODUS_FILE;
use mesher_gui::camera::Camera;
use mesher_gui::scene::Scene;
use mesher_gui::session::Session;
use netcdf3::FileReader;

/// Three unit tetrahedra side by side, as named objects.
fn three_tets() -> String {
    let mut s = String::new();
    for (k, name) in ["left", "middle", "right"].iter().enumerate() {
        let x = 3.0 * k as f64;
        s.push_str(&format!("o {name}\n"));
        for p in [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ] {
            s.push_str(&format!("v {} {} {}\n", p[0] + x, p[1], p[2]));
        }
        let b = 4 * k + 1;
        for f in [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]] {
            s.push_str(&format!("f {} {} {}\n", b + f[0], b + f[1], b + f[2]));
        }
    }
    s
}

fn binary(name: &str) -> PathBuf {
    // Tests run from target/<profile>/deps; binaries sit one level up.
    let exe = std::env::current_exe().unwrap();
    let dir = exe.parent().and_then(Path::parent).unwrap();
    let path = dir.join(format!("{name}{}", std::env::consts::EXE_SUFFIX));
    assert!(
        path.exists(),
        "{} not built; run `cargo test --workspace`",
        path.display()
    );
    path
}

/// Element block names and element counts of an Exodus file.
fn blocks(path: &Path) -> Vec<(String, usize)> {
    let (ds, _) = FileReader::open(path).unwrap().close();
    let mut r = FileReader::open(path).unwrap();
    let n = ds.dim_size("num_el_blk").unwrap();
    let raw = r.read_var_u8("eb_names").unwrap();
    let names: Vec<String> = raw
        .chunks(33)
        .map(|c| String::from_utf8_lossy(c).trim_end_matches('\0').to_owned())
        .collect();
    (0..n)
        .map(|i| {
            let count = ds.dim_size(&format!("num_el_in_blk{}", i + 1)).unwrap();
            (names[i].clone(), count)
        })
        .collect()
}

#[test]
fn selected_parts_export_one_block_each_and_replay_from_the_cli() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("three.obj");
    std::fs::write(&input, three_tets()).unwrap();

    // Open: the parts the worker will see, all selected.
    let scene = Scene::load(&input);
    let names: Vec<&str> = scene.parts.iter().map(|p| p.name.as_str()).collect();
    assert_eq!(names, ["left", "middle", "right"]);
    let mut session = Session::new(scene);
    assert!(session.selected.iter().all(|&s| s));

    // Pick the middle tet in the viewport: a ray from the camera through its
    // centroid toggles it off; then toggle "right" off and on in the tree.
    let cam = Camera::fit(session.scene.bounds);
    let centroid = [3.25f32, 0.25, 0.25];
    let eye = cam.eye();
    let d = [
        f64::from(centroid[0] - eye[0]),
        f64::from(centroid[1] - eye[1]),
        f64::from(centroid[2] - eye[2]),
    ];
    let picked = session.scene.pick(eye.map(f64::from), d).unwrap();
    assert_eq!(session.scene.parts[picked].name, "middle");
    session.toggle(picked);
    session.toggle(2);
    session.toggle(2);
    assert_eq!(session.part_args(), ["left", "right"]);

    // Export: save the project beside the outputs and run the worker.
    session.output = dir.path().join("gui_out");
    let project = session.project();
    std::fs::create_dir_all(&session.output).unwrap();
    project.write(&session.project_path()).unwrap();
    let limits = mesher_client::Limits {
        timeout: Duration::from_secs(300),
        memory: None,
    };
    let never = AtomicBool::new(false);
    let mut progress = 0usize;
    let (report, _) = mesher_client::run_job(
        &binary("mesher-worker"),
        &project.job(),
        limits,
        &mut |_| progress += 1,
        &never,
    );
    assert!(progress > 0, "messages stream to the caller");
    assert!(report.summary.exit_code <= 10, "{:?}", report.summary);
    let exported = blocks(&session.output.join(MESH_EXODUS_FILE));
    let block_names: Vec<&str> = exported.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(block_names, ["left", "right"]);
    assert!(exported.iter().all(|(_, n)| *n > 0));

    // Reproducible from the CLI: the saved project, into another directory.
    let saved = Project::read(&session.project_path()).unwrap();
    assert_eq!(saved, project);
    let cli_out = dir.path().join("cli_out");
    let run = Command::new(binary("mesher"))
        .arg("run")
        .arg(session.project_path())
        .arg("--out")
        .arg(&cli_out)
        .arg("--worker")
        .arg(binary("mesher-worker"))
        .output()
        .unwrap();
    assert!(run.status.code().is_some_and(|c| c <= 10), "{run:?}");
    assert_eq!(blocks(&cli_out.join(MESH_EXODUS_FILE)), exported);
}
