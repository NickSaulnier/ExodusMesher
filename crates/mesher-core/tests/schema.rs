//! Keeps `schema/report.schema.json` and the `PLAN.md` taxonomy table in step
//! with the code.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::collections::BTreeSet;
use std::path::PathBuf;

use mesher_core::DefectCode;
use mesher_core::report::Report;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

#[test]
fn committed_schema_is_current() {
    let path = repo_root().join("schema/report.schema.json");
    let generated = Report::json_schema();
    if std::env::var_os("UPDATE_SCHEMA").is_some() {
        std::fs::write(&path, &generated).unwrap();
        return;
    }
    let committed = std::fs::read_to_string(&path)
        .unwrap_or_default()
        .replace("\r\n", "\n");
    assert!(
        committed == generated,
        "schema/report.schema.json is stale; run `UPDATE_SCHEMA=1 cargo test -p mesher-core --test schema`"
    );
}

#[test]
fn schema_accepts_a_real_report() {
    let schema: serde_json::Value = serde_json::from_str(&Report::json_schema()).unwrap();
    let validator = jsonschema::validator_for(&schema).unwrap();
    let mut report = Report::new("in.obj", Default::default());
    let mut log = mesher_core::RecordLog::new();
    log.push(mesher_core::Record::new(DefectCode::ParseBadLine, mesher_core::Stage::Parse).line(1));
    report.finish(log);
    let instance: serde_json::Value = serde_json::from_str(&report.to_json()).unwrap();
    assert!(validator.is_valid(&instance));
}

/// The defect taxonomy table in PLAN.md and `DefectCode` must list the same codes.
#[test]
fn plan_taxonomy_matches_enum() {
    let plan = std::fs::read_to_string(repo_root().join("PLAN.md")).unwrap();
    let section = plan
        .split("**Defect taxonomy (v1).**")
        .nth(1)
        .expect("PLAN.md has a defect taxonomy section");
    let in_plan: BTreeSet<&str> = section
        .lines()
        .skip_while(|l| !l.starts_with('|'))
        .take_while(|l| l.starts_with('|'))
        .filter_map(|l| l.strip_prefix("| `"))
        .filter_map(|l| l.split('`').next())
        .collect();
    let in_code: BTreeSet<&str> = DefectCode::ALL.iter().map(|c| c.as_str()).collect();
    assert_eq!(in_plan, in_code);
}

/// TetGen (AGPL) may only appear in test-only crates under `tests/`.
#[test]
fn shipped_crates_do_not_depend_on_tetgen() {
    let crates = repo_root().join("crates");
    for entry in std::fs::read_dir(crates).unwrap() {
        let manifest = entry.unwrap().path().join("Cargo.toml");
        if let Ok(text) = std::fs::read_to_string(&manifest) {
            assert!(
                !text.contains("tritet"),
                "{} depends on tritet",
                manifest.display()
            );
        }
    }
}
