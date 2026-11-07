//! End-to-end tests of the worker binary: job on stdin, JSON lines on stdout,
//! `report.json` on disk.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use mesher_core::protocol::{Job, PROTOCOL_VERSION, WorkerMessage};
use mesher_core::{DefectCode, JobConfig, Report, Stage};

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/obj")
        .join(name)
        .display()
        .to_string()
}

struct Run {
    exit_code: i32,
    messages: Vec<WorkerMessage>,
    report_file: Option<Report>,
}

fn job_for(input: &str, out: &Path) -> Job {
    Job {
        input: input.to_owned(),
        output_dir: out.display().to_string(),
        config: JobConfig::default(),
    }
}

fn run_worker(input: &str, out: &Path) -> Run {
    run_job(&job_for(input, out), out)
}

fn run_job(job: &Job, out: &Path) -> Run {
    let mut child = Command::new(env!("CARGO_BIN_EXE_mesher-worker"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(serde_json::to_string(job).unwrap().as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();
    let messages = stdout
        .lines()
        .map(|l| WorkerMessage::from_line(l).unwrap_or_else(|e| panic!("bad line {l:?}: {e}")))
        .collect();
    let report_file = std::fs::read_to_string(out.join("report.json"))
        .ok()
        .map(|text| {
            validate(&text);
            serde_json::from_str(&text).unwrap()
        });
    Run {
        exit_code: output.status.code().unwrap(),
        messages,
        report_file,
    }
}

fn validate(report_json: &str) {
    let schema: serde_json::Value = serde_json::from_str(&Report::json_schema()).unwrap();
    let validator = jsonschema::validator_for(&schema).unwrap();
    let instance: serde_json::Value = serde_json::from_str(report_json).unwrap();
    let errors: Vec<String> = validator
        .iter_errors(&instance)
        .map(|e| e.to_string())
        .collect();
    assert!(errors.is_empty(), "report fails schema: {errors:?}");
}

fn result_of(run: &Run) -> (&Report, &Option<String>) {
    match run.messages.last() {
        Some(WorkerMessage::Result {
            report,
            report_path,
            exit_code,
        }) => {
            assert_eq!(*exit_code, run.exit_code);
            (report, report_path)
        }
        other => panic!("last message is not a result: {other:?}"),
    }
}

#[test]
fn clean_cube_exits_zero() {
    let dir = tempfile::tempdir().unwrap();
    let run = run_worker(&fixture("cube.obj"), dir.path());
    assert_eq!(run.exit_code, 0);
    assert!(matches!(
        run.messages.first(),
        Some(WorkerMessage::Started { protocol, .. }) if *protocol == PROTOCOL_VERSION
    ));
    let (report, path) = result_of(&run);
    assert!(path.is_some());
    assert_eq!(Some(report), run.report_file.as_ref());
    let stats = report.input.stats.as_ref().unwrap();
    assert_eq!(
        (stats.vertices, stats.polygons, stats.triangles),
        (8, 6, 12)
    );
}

#[test]
fn parse_defects_are_streamed_and_reported() {
    let dir = tempfile::tempdir().unwrap();
    let run = run_worker(&fixture("cube_parse_defects.obj"), dir.path());
    assert_eq!(run.exit_code, 10);
    let streamed: Vec<(DefectCode, Vec<u32>)> = run
        .messages
        .iter()
        .filter_map(|m| match m {
            WorkerMessage::Record { record } => {
                Some((record.code, record.location.obj_lines.clone()))
            }
            _ => None,
        })
        .collect();
    assert_eq!(
        streamed,
        vec![
            (DefectCode::ParseBadLine, vec![11]),
            (DefectCode::ParseUnsupportedElement, vec![19]),
            (DefectCode::ParseBadLine, vec![20]),
            (DefectCode::ParseBadIndex, vec![18]),
        ]
    );
    let report = run.report_file.unwrap();
    assert_eq!(report.records.len(), 4);
    assert_eq!(report.summary.warnings, 3);
    assert_eq!(report.summary.info, 1);
}

#[test]
fn multi_object_and_concave_crlf_inputs() {
    let dir = tempfile::tempdir().unwrap();
    let run = run_worker(&fixture("two_objects.obj"), dir.path());
    assert_eq!(run.exit_code, 0);
    let report = run.report_file.unwrap();
    let stats = report.input.stats.as_ref().unwrap();
    assert_eq!(stats.objects, vec!["left", "right"]);
    assert_eq!(stats.groups, vec!["left_shell"]);
    let parts: Vec<(&str, u64, bool)> = report
        .parts
        .iter()
        .map(|p| (p.name.as_str(), p.triangles, p.selected))
        .collect();
    assert_eq!(parts, vec![("left", 4, true), ("right", 4, true)]);
    assert_eq!(
        report.stages_run,
        vec![Stage::Parse, Stage::Diagnose, Stage::Select]
    );

    // A single flat pentagon: parsed into 3 triangles, then diagnosed as an
    // open sheet, which is fatal for its part.
    let dir = tempfile::tempdir().unwrap();
    let run = run_worker(&fixture("concave_pentagon_crlf.obj"), dir.path());
    assert_eq!(run.exit_code, 30);
    let report = run.report_file.unwrap();
    assert_eq!(report.input.stats.unwrap().triangles, 3);
    let codes: Vec<DefectCode> = report.records.iter().map(|r| r.code).collect();
    assert_eq!(codes, vec![DefectCode::TopOpenSheet]);
    assert_eq!(report.records[0].location.faces, vec![0, 1, 2]);
}

#[test]
fn diagnosis_writes_an_overlay_and_honours_selection() {
    let dir = tempfile::tempdir().unwrap();
    let run = run_worker(&fixture("cube.obj"), dir.path());
    assert_eq!(run.exit_code, 0);
    let vtp = std::fs::read_to_string(dir.path().join("defects.vtp")).unwrap();
    assert!(vtp.contains(r#"NumberOfPoints="8""#) && vtp.contains(r#"NumberOfPolys="12""#));

    let dir = tempfile::tempdir().unwrap();
    let mut job = job_for(&fixture("two_objects.obj"), dir.path());
    job.config.parts = vec!["right".into()];
    let run = run_job(&job, dir.path());
    assert_eq!(run.exit_code, 0);
    let selected: Vec<bool> = run
        .report_file
        .unwrap()
        .parts
        .iter()
        .map(|p| p.selected)
        .collect();
    assert_eq!(selected, vec![false, true]);

    let dir = tempfile::tempdir().unwrap();
    let mut job = job_for(&fixture("two_objects.obj"), dir.path());
    job.config.parts = vec!["middle".into()];
    let run = run_job(&job, dir.path());
    assert_eq!(run.exit_code, 30);
    let report = run.report_file.unwrap();
    assert_eq!(report.records[0].code, DefectCode::PartEmptySelection);
}

#[test]
fn missing_input_is_a_fatal_report_not_a_crash() {
    let dir = tempfile::tempdir().unwrap();
    let run = run_worker(&fixture("does_not_exist.obj"), dir.path());
    assert_eq!(run.exit_code, 30);
    let report = run.report_file.unwrap();
    assert_eq!(report.records[0].code, DefectCode::ParseReadFailed);
    assert!(report.input.stats.is_none());
}

#[test]
fn unwritable_output_dir_is_reported_on_stdout() {
    let dir = tempfile::tempdir().unwrap();
    // A file where the output directory should be.
    let blocker = dir.path().join("blocker");
    std::fs::write(&blocker, b"").unwrap();
    let run = run_worker(&fixture("cube.obj"), &blocker.join("out"));
    assert_eq!(run.exit_code, 30);
    let (report, path) = result_of(&run);
    assert!(path.is_none());
    assert!(
        report
            .records
            .iter()
            .any(|r| r.code == DefectCode::OutWriteFailed)
    );
}

#[test]
fn reruns_are_byte_identical() {
    let a = tempfile::tempdir().unwrap();
    let b = tempfile::tempdir().unwrap();
    run_worker(&fixture("cube_parse_defects.obj"), a.path());
    run_worker(&fixture("cube_parse_defects.obj"), b.path());
    let read = |d: &Path| std::fs::read(d.join("report.json")).unwrap();
    assert_eq!(read(a.path()), read(b.path()));
}

#[test]
fn malformed_job_exits_70() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_mesher-worker"))
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(b"{not json").unwrap();
    assert_eq!(child.wait().unwrap().code(), Some(70));
}
