//! The Python module `mesher._mesher`: one function, [`_run`].
//!
//! The job runs on a thread of its own while this thread, holding the
//! interpreter only between waits, passes progress to the Python callback
//! and checks Python's signals, so Ctrl-C (or an exception from the
//! callback) cancels the worker. Cancelling still ends in a report; the
//! exception is raised after it is written.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, mpsc};
use std::time::Duration;

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use mesher_client::Limits;
use mesher_core::protocol::WorkerMessage;

/// How often the waiting thread wakes to check Python's signals.
const POLL: Duration = Duration::from_millis(100);

/// Run a job (its JSON) in the worker at `worker`; returns the report's JSON
/// and the path of `report.json` (or `None`). `progress(stage, fraction)`
/// is called as the worker reports progress.
#[pyfunction]
#[pyo3(signature = (job_json, worker, timeout, memory=None, progress=None))]
fn _run(
    py: Python<'_>,
    job_json: &str,
    worker: PathBuf,
    timeout: f64,
    memory: Option<usize>,
    progress: Option<Py<PyAny>>,
) -> PyResult<(String, Option<String>)> {
    let timeout = Duration::try_from_secs_f64(timeout)
        .map_err(|_| PyValueError::new_err("timeout must be a positive number of seconds"))?;
    let limits = Limits { timeout, memory };
    let cancel = AtomicBool::new(false);
    let (tx, rx) = mpsc::channel::<(String, f64)>();
    // Waited on with the interpreter released, so it must be `Sync`.
    let rx = Mutex::new(rx);
    let mut raised: Option<PyErr> = None;

    let result = std::thread::scope(|s| {
        let handle = s.spawn(|| {
            let tx = tx;
            crate::run(
                job_json,
                &worker,
                limits,
                &mut |m| {
                    if let WorkerMessage::Progress { stage, fraction } = m {
                        let stage = serde_json::to_value(stage)
                            .ok()
                            .and_then(|v| v.as_str().map(str::to_owned))
                            .unwrap_or_default();
                        let _ = tx.send((stage, *fraction));
                    }
                },
                &cancel,
            )
        });
        loop {
            let next = py.detach(|| match rx.lock() {
                Ok(rx) => rx.recv_timeout(POLL),
                Err(_) => Err(mpsc::RecvTimeoutError::Disconnected),
            });
            match next {
                Ok((stage, fraction)) => {
                    if let (Some(cb), None) = (&progress, &raised)
                        && let Err(e) = cb.call1(py, (stage, fraction))
                    {
                        raised = Some(e);
                        cancel.store(true, Ordering::SeqCst);
                    }
                }
                Err(mpsc::RecvTimeoutError::Timeout) => {}
                Err(mpsc::RecvTimeoutError::Disconnected) => break,
            }
            if raised.is_none()
                && let Err(e) = py.check_signals()
            {
                raised = Some(e);
                cancel.store(true, Ordering::SeqCst);
            }
        }
        py.detach(|| handle.join())
    });
    if let Some(e) = raised {
        return Err(e);
    }
    match result {
        Ok(Ok(run)) => Ok((run.report_json, run.report_path)),
        Ok(Err(message)) => Err(PyValueError::new_err(message)),
        Err(_) => Err(PyRuntimeError::new_err(
            "the job's thread panicked (a bug in mesher-client)",
        )),
    }
}

#[pymodule]
fn _mesher(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_run, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
