//! Milestone 1 exit criterion: any file, including garbage, yields a valid
//! `report.json` and never a panic. Drives the real `mesher` and
//! `mesher-worker` binaries.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use mesher_core::protocol::TEST_FAULT_ENV;
use mesher_core::{DefectCode, Kernel, Outcome, Report, Severity};
use proptest::prelude::*;

fn worker_path() -> PathBuf {
    let mesher = PathBuf::from(env!("CARGO_BIN_EXE_mesher"));
    let worker = mesher.with_file_name(format!("mesher-worker{}", std::env::consts::EXE_SUFFIX));
    assert!(
        worker.exists(),
        "{} not built; run `cargo test --workspace` or `cargo build -p mesher-worker`",
        worker.display()
    );
    worker
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/obj")
        .join(name)
}

fn mesher(input: &Path, out: &Path, extra: &[&str], fault: Option<&str>) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_mesher"));
    cmd.arg("run")
        .arg(input)
        .arg("--out")
        .arg(out)
        .arg("--worker")
        .arg(worker_path())
        .args(extra);
    match fault {
        Some(f) => cmd.env(TEST_FAULT_ENV, f),
        None => cmd.env_remove(TEST_FAULT_ENV),
    };
    cmd.output().unwrap()
}

/// Read and schema-validate `report.json`.
fn report(out: &Path) -> Report {
    let text = std::fs::read_to_string(out.join("report.json")).expect("report.json written");
    let schema: serde_json::Value = serde_json::from_str(&Report::json_schema()).unwrap();
    let validator = jsonschema::validator_for(&schema).unwrap();
    let instance: serde_json::Value = serde_json::from_str(&text).unwrap();
    let errors: Vec<String> = validator
        .iter_errors(&instance)
        .map(|e| e.to_string())
        .collect();
    assert!(errors.is_empty(), "schema errors: {errors:?}");
    serde_json::from_str(&text).unwrap()
}

/// Run on `bytes` and assert the exit criterion; returns the exit code.
fn assert_handles(bytes: &[u8]) -> i32 {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.obj");
    std::fs::write(&input, bytes).unwrap();
    let out = dir.path().join("out");
    let output = mesher(&input, &out, &[], None);
    let code = output.status.code().unwrap();
    assert!(
        [0, 10, 20, 30].contains(&code),
        "exit {code}; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report = report(&out);
    assert_eq!(report.summary.exit_code, code);
    assert!(
        !report
            .records
            .iter()
            .any(|r| r.code == DefectCode::MeshWorkerFailed)
    );
    code
}

#[test]
fn clean_and_defective_fixtures() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("cube");
    let output = mesher(&fixture("cube.obj"), &out, &[], None);
    assert_eq!(output.status.code(), Some(0));
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.contains("8 vertices, 6 polygons, 12 triangles"),
        "{stdout}"
    );
    assert!(stdout.contains("(Delaunay)"), "{stdout}");
    assert!(stdout.contains("exit 0: clean"), "{stdout}");

    let out = dir.path().join("defects");
    let output = mesher(&fixture("cube_parse_defects.obj"), &out, &[], None);
    assert_eq!(output.status.code(), Some(10));
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.contains("PARSE_BAD_LINE x2, first at line 11"),
        "{stdout}"
    );
    assert_eq!(report(&out).records.len(), 4);
}

#[test]
fn hostile_inputs_yield_reports() {
    let mut huge_line = b"f ".to_vec();
    huge_line.extend(std::iter::repeat_n(b'9', 1 << 20));
    let mut many_bad_lines = Vec::new();
    for _ in 0..5000 {
        many_bad_lines.extend_from_slice(b"v 1 2\n");
    }
    let mut huge_polygon = String::new();
    for i in 0..2000 {
        let a = f64::from(i) / 2000.0 * std::f64::consts::TAU;
        huge_polygon.push_str(&format!("v {} {} 0\n", libm::cos(a), libm::sin(a)));
    }
    huge_polygon.push('f');
    for i in 1..=2000 {
        huge_polygon.push_str(&format!(" {i}"));
    }
    let cases: Vec<Vec<u8>> = vec![
        Vec::new(),
        b"\x00\x01\x02\xff\xfe\xfd".repeat(1000),
        b"\n\n\r\r\n\\\\\\".to_vec(),
        huge_line,
        b"v 1 2 3\nf 1 1 1\n".to_vec(),
        b"f -1 -2 -3\n".to_vec(),
        b"v 1e308 1e308 1e308\nv -1e308 0 0\nv 0 -1e308 0\nf 1 2 3\n".to_vec(),
        many_bad_lines,
        huge_polygon.into_bytes(),
    ];
    for case in cases {
        assert_handles(&case);
    }
}

#[test]
fn record_detail_is_capped_but_counted() {
    let mut text = Vec::new();
    for _ in 0..5000 {
        text.extend_from_slice(b"v 1 2\n");
    }
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("bad.obj");
    std::fs::write(&input, text).unwrap();
    let out = dir.path().join("out");
    mesher(&input, &out, &[], None);
    let r = report(&out);
    let bad = r.summary.by_code[&DefectCode::ParseBadLine];
    assert_eq!((bad.total, bad.reported), (5000, 100));
}

#[test]
fn directory_as_input_is_a_read_failure() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    let output = mesher(dir.path(), &out, &[], None);
    assert_eq!(output.status.code(), Some(30));
    assert_eq!(report(&out).records[0].code, DefectCode::ParseReadFailed);
}

#[test]
fn worker_panic_exits_70_with_report() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    let output = mesher(&fixture("cube.obj"), &out, &[], Some("panic"));
    assert_eq!(output.status.code(), Some(70));
    let r = report(&out);
    let failure = r.records.last().unwrap();
    assert_eq!(failure.code, DefectCode::MeshWorkerFailed);
    assert_eq!(failure.measured["reason"], "died");
    assert!(
        failure.measured["stderr_tail"]
            .as_str()
            .unwrap()
            .contains("deliberate worker panic")
    );
    assert_eq!(r.summary.exit_code, 70);
}

#[test]
fn worker_hang_times_out_with_report() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    let started = std::time::Instant::now();
    let output = mesher(
        &fixture("cube.obj"),
        &out,
        &["--timeout", "1"],
        Some("hang"),
    );
    assert!(started.elapsed() < std::time::Duration::from_secs(30));
    assert_eq!(output.status.code(), Some(30));
    let r = report(&out);
    assert_eq!(r.records.last().unwrap().measured["reason"], "timeout");
}

/// Two overlapping cubes in one object: repair rung 2 replaces them by the
/// boundary of their union (Milestone 16; they went to the stuffing
/// fallback before), so the part meshes conformally, the overlap findings
/// are warnings with `resolve` / `applied`, and the run exits 10, also under
/// `--strict` (rung 2 is a repair, not a fallback).
#[test]
fn overlapping_shells_resolved_by_rung_2_exit_10() {
    let mut obj = String::from("o overlap\n");
    for (k, offset) in [0.0, 0.5].into_iter().enumerate() {
        for v in 0..8 {
            let c = |bit: usize| f64::from(u8::from(v & bit != 0)) + offset;
            obj.push_str(&format!("v {} {} {}\n", c(1), c(2), c(4)));
        }
        let b = k * 8;
        for f in [
            [1, 3, 4, 2],
            [5, 6, 8, 7],
            [1, 2, 6, 5],
            [2, 4, 8, 6],
            [4, 3, 7, 8],
            [3, 1, 5, 7],
        ] {
            let f = f.map(|i| (i + b).to_string());
            obj.push_str(&format!("f {}\n", f.join(" ")));
        }
    }
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("overlap.obj");
    std::fs::write(&input, obj).unwrap();
    for (dir_name, args) in [("out", &[][..]), ("strict", &["--strict"][..])] {
        let out = dir.path().join(dir_name);
        let output = mesher(&input, &out, args, None);
        assert_eq!(output.status.code(), Some(10), "{dir_name}: {output:?}");
        let r = report(&out);
        assert_eq!(r.summary.errors, 0);
        let overlap = r
            .records
            .iter()
            .find(|r| {
                r.code == DefectCode::IntShellOverlap && r.stage == mesher_core::Stage::Diagnose
            })
            .expect("the overlap is diagnosed");
        assert_eq!(overlap.severity, Severity::Warning);
        assert_eq!(overlap.action, Some(mesher_core::Action::Resolve));
        assert_eq!(overlap.outcome, Some(Outcome::Applied));
        let part = &r.parts[0];
        assert_eq!(part.repair.as_ref().map(|s| s.rung), Some(2));
        assert!(
            part.mesh
                .as_ref()
                .is_some_and(|m| m.accepted && m.kernel == Kernel::Delaunay)
        );
    }
}

/// A worker that times out on the conformal path is retried once with the
/// stuffing fallback only; the report keeps the first failure as a warning.
#[test]
fn worker_timeout_falls_back_to_stuffing() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    let output = mesher(
        &fixture("cube.obj"),
        &out,
        &["--timeout", "2"],
        Some("hang_conformal"),
    );
    assert_eq!(output.status.code(), Some(20), "{output:?}");
    let r = report(&out);
    let first = &r.records[0];
    assert_eq!(first.code, DefectCode::MeshWorkerFailed);
    assert_eq!(first.measured["reason"], "timeout");
    assert_eq!(first.severity, Severity::Warning);
    assert_eq!(first.outcome, Some(Outcome::Applied));
    assert!(r.config.stuffing_only);
    assert_eq!(r.attempts[0].kernel, Kernel::Delaunay);
    assert!(
        r.parts[0]
            .mesh
            .as_ref()
            .is_some_and(|m| m.accepted && m.kernel == Kernel::Stuffing)
    );
    assert!(out.join("mesh.exo").exists());

    // --strict forbids the fallback: the timeout is the result.
    let strict = dir.path().join("strict");
    let output = mesher(
        &fixture("cube.obj"),
        &strict,
        &["--timeout", "2", "--strict"],
        Some("hang_conformal"),
    );
    assert_eq!(output.status.code(), Some(30));
    assert_eq!(
        report(&strict).records.last().unwrap().measured["reason"],
        "timeout"
    );
}

#[test]
fn missing_worker_exits_70_with_report() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    let output = Command::new(env!("CARGO_BIN_EXE_mesher"))
        .args(["run", "x.obj", "--worker"])
        .arg(dir.path().join("no-such-worker"))
        .arg("--out")
        .arg(&out)
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(70));
    assert_eq!(report(&out).records[0].measured["reason"], "spawn_failed");
}

#[test]
fn schema_subcommand_prints_the_schema() {
    let output = Command::new(env!("CARGO_BIN_EXE_mesher"))
        .arg("schema")
        .output()
        .unwrap();
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        Report::json_schema()
    );
}

#[test]
fn bad_timeout_is_a_usage_error() {
    let output = Command::new(env!("CARGO_BIN_EXE_mesher"))
        .args(["run", "x.obj", "--timeout", "-1"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn reruns_are_byte_identical() {
    let dir = tempfile::tempdir().unwrap();
    let (a, b) = (dir.path().join("a"), dir.path().join("b"));
    mesher(&fixture("cube_parse_defects.obj"), &a, &[], None);
    mesher(&fixture("cube_parse_defects.obj"), &b, &[], None);
    let read = |d: &Path| std::fs::read(d.join("report.json")).unwrap();
    assert_eq!(read(&a), read(&b));
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn random_bytes_yield_reports(bytes in prop::collection::vec(any::<u8>(), 0..4096)) {
        assert_handles(&bytes);
    }
}

#[test]
fn memory_cap_stops_a_runaway_worker() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    let output = mesher(
        &fixture("cube.obj"),
        &out,
        &["--memory-limit", "64"],
        Some("alloc"),
    );
    let r = report(&out);
    let failure = r.records.last().unwrap();
    assert_eq!(
        failure.code,
        DefectCode::MeshWorkerFailed,
        "{:?}",
        r.records
    );
    assert_eq!(failure.measured["reason"], "out_of_memory", "{failure:?}");
    // Running out of memory is a failure, not a bug.
    assert_eq!(output.status.code(), Some(30));
}

#[test]
fn project_file_replays_a_selection() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::copy(fixture("two_objects.obj"), dir.path().join("two.obj")).unwrap();
    let project = dir.path().join("job.toml");
    std::fs::write(
        &project,
        "input = \"two.obj\"\noutput = \"out\"\n\n[settings]\nparts = [\"right\"]\n",
    )
    .unwrap();
    let run = |extra: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_mesher"))
            .arg("run")
            .arg(&project)
            .arg("--worker")
            .arg(worker_path())
            .args(extra)
            .env_remove(TEST_FAULT_ENV)
            .output()
            .unwrap()
    };
    // Paths resolve against the project file, not the current directory.
    let out = run(&[]);
    assert_eq!(out.status.code(), Some(0), "{out:?}");
    let r = report(&dir.path().join("out"));
    let selected: Vec<&str> = r
        .parts
        .iter()
        .filter(|p| p.selected)
        .map(|p| p.name.as_str())
        .collect();
    assert_eq!(selected, ["right"]);
    assert!(dir.path().join("out/mesh.exo").exists());
    // Options given on the command line override the project.
    let elsewhere = dir.path().join("elsewhere");
    let out = run(&["--out", elsewhere.to_str().unwrap(), "--parts", "left"]);
    assert_eq!(out.status.code(), Some(0), "{out:?}");
    let r = report(&elsewhere);
    assert!(r.parts.iter().all(|p| p.selected == (p.name == "left")));
}

#[test]
fn invalid_project_file_is_reported() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("bad.toml");
    std::fs::write(&project, "input = \"x.obj\"\nunknown = 3\n").unwrap();
    let out = dir.path().join("out");
    let run = mesher(&project, &out, &[], None);
    assert_eq!(run.status.code(), Some(30), "{run:?}");
    let r = report(&out);
    assert!(
        r.records
            .iter()
            .any(|x| x.code == DefectCode::ParseProjectInvalid)
    );
}
