//! `mesher-py`: the Rust side of ExodusMesher's Python API (Milestone 20).
//!
//! It is deliberately thin. A job arrives as the JSON of a
//! [`Job`](mesher_core::protocol::Job), runs in a `mesher-worker` process
//! through [`mesher_client::run_job`] (so a panic or a runaway job never
//! reaches the caller's interpreter), and its report goes back as JSON text
//! with the path of `report.json`. Reading a run's mesh back is [`mesh`]
//! (its `.meshb`). Everything else, the typed report included, is Python
//! (`python/mesher/`).
//!
//! The Python module itself (`mesher._mesher`) is behind the `python`
//! feature, which maturin turns on; without it this crate is just
//! [`run`], so the workspace builds and tests without Python.

use std::path::Path;
use std::sync::atomic::AtomicBool;

use mesher_client::{JobRequest, Limits, RequestError};
use mesher_core::protocol::{Job, WorkerMessage};

pub mod mesh;

#[cfg(feature = "python")]
#[allow(unsafe_code)]
mod python;

/// A finished run: the report as JSON, and where `report.json` was written
/// (`None` if it could not be).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Run {
    pub report_json: String,
    pub report_path: Option<String>,
}

/// Run the job given as JSON in the worker at `worker`.
///
/// # Errors
/// The JSON is not a valid job (nothing runs). Everything after that,
/// failures of the worker included, is in the report, as for the CLI.
pub fn run(
    job_json: &str,
    worker: &Path,
    limits: Limits,
    on_message: &mut dyn FnMut(&WorkerMessage),
    cancel: &AtomicBool,
) -> Result<Run, String> {
    let job: Job = serde_json::from_str(job_json).map_err(|e| format!("not a valid job: {e}"))?;
    let (report, report_path) = mesher_client::run_job(worker, &job, limits, on_message, cancel);
    Ok(Run {
        report_json: report.to_json(),
        report_path,
    })
}

/// Run a request (a [`JobRequest`] as JSON: an input or project file and
/// the CLI's options), exactly as `mesher run` does: its job in the worker,
/// or, for a project file that cannot be used, the report that says so.
///
/// # Errors
/// The JSON is not a valid request, or an option is out of range (nothing
/// runs, no report is written).
pub fn run_request(
    request_json: &str,
    worker: &Path,
    limits: Limits,
    on_message: &mut dyn FnMut(&WorkerMessage),
    cancel: &AtomicBool,
) -> Result<Run, String> {
    let request = parse_request(request_json)?;
    let (report, report_path) =
        mesher_client::run_request(worker, &request, limits, on_message, cancel)?;
    Ok(Run {
        report_json: report.to_json(),
        report_path,
    })
}

/// The job a request builds, as JSON, or `None` for a project file that
/// cannot be used (running it gives a `PARSE_PROJECT_INVALID` report).
///
/// # Errors
/// The JSON is not a valid request, or an option is out of range.
pub fn build_job(request_json: &str) -> Result<Option<String>, String> {
    match parse_request(request_json)?.job() {
        Ok(job) => serde_json::to_string(&job)
            .map(Some)
            .map_err(|e| e.to_string()),
        Err(RequestError::Project(_)) => Ok(None),
        Err(RequestError::Invalid(e)) => Err(e),
    }
}

fn parse_request(request_json: &str) -> Result<JobRequest, String> {
    serde_json::from_str(request_json).map_err(|e| format!("not a valid request: {e}"))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;
    use std::time::Duration;

    #[test]
    fn a_job_that_is_not_json_is_refused_before_anything_runs() {
        let limits = Limits {
            timeout: Duration::from_secs(1),
            memory: None,
        };
        let never = AtomicBool::new(false);
        for bad in ["", "{", "{\"input\": 3}", "[]"] {
            let r = run(
                bad,
                Path::new("no-such-worker"),
                limits,
                &mut |_| {},
                &never,
            );
            assert!(r.unwrap_err().starts_with("not a valid job"), "{bad}");
        }
    }

    #[test]
    fn a_request_builds_the_clis_job() {
        let job = build_job(r#"{"input": "m.obj", "element_size": 0.05}"#)
            .unwrap()
            .unwrap();
        let job: Job = serde_json::from_str(&job).unwrap();
        assert_eq!(job.output_dir, "m_out");
        assert_eq!(job.config.element_size, 0.05);
        assert!(build_job(r#"{"input": "m.obj", "element_size": -1}"#).is_err());
        assert!(build_job(r#"{"input": "m.obj", "colour": 1}"#).is_err());
        assert_eq!(build_job(r#"{"input": "missing.toml"}"#), Ok(None));
    }
}
