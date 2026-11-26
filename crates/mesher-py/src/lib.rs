//! `mesher-py`: the Rust side of ExodusMesher's Python API (Milestone 20).
//!
//! It is deliberately thin. A job arrives as the JSON of a
//! [`Job`](mesher_core::protocol::Job), runs in a `mesher-worker` process
//! through [`mesher_client::run_job`] (so a panic or a runaway job never
//! reaches the caller's interpreter), and its report goes back as JSON text
//! with the path of `report.json`. Everything else, the typed report
//! included, is Python (`python/mesher/`).
//!
//! The Python module itself (`mesher._mesher`) is behind the `python`
//! feature, which maturin turns on; without it this crate is just
//! [`run`], so the workspace builds and tests without Python.

use std::path::Path;
use std::sync::atomic::AtomicBool;

use mesher_client::Limits;
use mesher_core::protocol::{Job, WorkerMessage};

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
}
