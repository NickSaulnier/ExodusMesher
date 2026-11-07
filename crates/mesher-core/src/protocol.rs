//! The worker protocol: the job goes to the worker as one JSON document on stdin,
//! and the worker answers with JSON lines on stdout, one [`WorkerMessage`] each.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::report::{JobConfig, Record, Report, Stage};

pub const PROTOCOL_VERSION: u32 = 1;

/// Name of the report file the worker writes into the job's output directory.
pub const REPORT_FILE: &str = "report.json";

/// Environment variable that makes the worker fail on purpose (`panic` or
/// `hang`), so tests can exercise the front ends' failure handling.
pub const TEST_FAULT_ENV: &str = "MESHER_TEST_FAULT";

/// One meshing job.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Job {
    /// Path of the input .obj file, as the user gave it.
    pub input: String,
    /// Directory the worker writes `report.json` and meshes into; created if missing.
    pub output_dir: String,
    #[serde(default)]
    pub config: JobConfig,
}

/// One line of worker stdout.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum WorkerMessage {
    Started {
        protocol: u32,
        worker_version: String,
    },
    Progress {
        stage: Stage,
        /// Fraction of the stage done, in `[0, 1]`.
        fraction: f64,
    },
    Record {
        record: Box<Record>,
    },
    Result {
        exit_code: i32,
        /// Path of the written report, or `None` if writing it failed (the
        /// report then carries an `OUT_WRITE_FAILED` record).
        report_path: Option<String>,
        report: Box<Report>,
    },
}

impl WorkerMessage {
    /// Serialize as a single protocol line (no trailing newline).
    #[must_use]
    pub fn to_line(&self) -> String {
        // Serializing these plain data types cannot fail: no maps with non-string
        // keys, no fallible Serialize impls.
        serde_json::to_string(self).unwrap_or_default()
    }

    /// Parse one protocol line.
    pub fn from_line(line: &str) -> Result<WorkerMessage, serde_json::Error> {
        serde_json::from_str(line)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::{DefectCode, Record};

    #[test]
    fn messages_round_trip_as_tagged_lines() {
        let messages = [
            WorkerMessage::Started {
                protocol: PROTOCOL_VERSION,
                worker_version: "0.1.0".into(),
            },
            WorkerMessage::Progress {
                stage: Stage::Parse,
                fraction: 0.5,
            },
            WorkerMessage::Record {
                record: Box::new(Record::new(DefectCode::ParseBadLine, Stage::Parse).line(3)),
            },
        ];
        for m in messages {
            let line = m.to_line();
            assert!(!line.contains('\n'));
            assert_eq!(WorkerMessage::from_line(&line).ok(), Some(m));
        }
        let progress = WorkerMessage::Progress {
            stage: Stage::Parse,
            fraction: 1.0,
        };
        assert!(progress.to_line().starts_with(r#"{"type":"progress""#));
    }

    #[test]
    fn job_rejects_unknown_fields() {
        let ok = r#"{"input":"a.obj","output_dir":"out"}"#;
        assert!(serde_json::from_str::<Job>(ok).is_ok());
        let bad = r#"{"input":"a.obj","output_dir":"out","bogus":1}"#;
        assert!(serde_json::from_str::<Job>(bad).is_err());
    }
}
