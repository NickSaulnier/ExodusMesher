//! The Python module `mesher._mesher`: [`_run`] (a job), [`_run_request`]
//! (a request, as `mesher run` takes it), [`_build_job`],
//! [`_read_mesh`] (a run's mesh, as bytes for numpy or `array`) and
//! [`_next_job`] (a report's first applicable remedy, as a job).
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
use pyo3::types::{PyBytes, PyDict};

use mesher_client::Limits;
use mesher_core::protocol::WorkerMessage;

/// How often the waiting thread wakes to check Python's signals.
const POLL: Duration = Duration::from_millis(100);

/// Run `run` (one of the crate's runners) on a thread of its own, passing
/// progress to `progress(stage, fraction)` and cancelling it on Ctrl-C or an
/// exception from the callback, which is raised once it returns.
fn supervise<F>(
    py: Python<'_>,
    timeout: f64,
    memory: Option<usize>,
    progress: Option<Py<PyAny>>,
    run: F,
) -> PyResult<(String, Option<String>)>
where
    F: FnOnce(Limits, &mut dyn FnMut(&WorkerMessage), &AtomicBool) -> Result<crate::Run, String>
        + Send,
{
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
            run(
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
    supervise(py, timeout, memory, progress, |limits, on, cancel| {
        crate::run(job_json, &worker, limits, on, cancel)
    })
}

/// Run a request (its JSON: an input or project file and the CLI's options)
/// as `mesher run` does; returns as [`_run`].
#[pyfunction]
#[pyo3(signature = (request_json, worker, timeout, memory=None, progress=None))]
fn _run_request(
    py: Python<'_>,
    request_json: &str,
    worker: PathBuf,
    timeout: f64,
    memory: Option<usize>,
    progress: Option<Py<PyAny>>,
) -> PyResult<(String, Option<String>)> {
    supervise(py, timeout, memory, progress, |limits, on, cancel| {
        crate::run_request(request_json, &worker, limits, on, cancel)
    })
}

/// The job a request builds, as JSON, or `None` for a project file that
/// cannot be used.
#[pyfunction]
fn _build_job(request_json: &str) -> PyResult<Option<String>> {
    crate::build_job(request_json).map_err(PyValueError::new_err)
}

/// The job to run next after a report (its JSON), writing to
/// `output_dir`: the job's JSON with the record's and the remedy's
/// indices, or `None` when no record at `min` severity or above has a
/// remedy that applies.
#[pyfunction]
fn _next_job(
    report_json: &str,
    output_dir: &str,
    min: &str,
) -> PyResult<Option<(String, usize, usize)>> {
    crate::next_job(report_json, output_dir, min).map_err(PyValueError::new_err)
}

/// Native-endian bytes of 64-bit values.
fn ints<'py, const N: usize>(py: Python<'py>, rows: &[[u32; N]]) -> Bound<'py, PyBytes> {
    let mut out = Vec::with_capacity(rows.len() * N * 8);
    for row in rows {
        for &v in row {
            out.extend_from_slice(&i64::from(v).to_ne_bytes());
        }
    }
    PyBytes::new(py, &out)
}

/// The mesh at `path` (a run's output directory or a `.meshb` file): a dict
/// of `points` (3 doubles each), `tets` (4 int64 each, 0-based),
/// `tet_parts`, `triangles` (3 each) and `triangle_patches` (one int64
/// each) as native-endian bytes, and the names `parts` and `patches`.
/// A missing or unreadable file raises `OSError`; one that is not a
/// complete mesh, `ValueError`.
#[pyfunction]
fn _read_mesh(py: Python<'_>, path: PathBuf) -> PyResult<Bound<'_, PyDict>> {
    let m = py
        .detach(|| crate::mesh::read_mesh(&path))
        .map_err(|e| match e {
            crate::mesh::ReadError::Io(p, e) => {
                let message = format!("{}: {e}", p.display());
                PyErr::from(std::io::Error::new(e.kind(), message))
            }
            e @ crate::mesh::ReadError::Invalid(..) => PyValueError::new_err(e.to_string()),
        })?;
    let mut points = Vec::with_capacity(m.points.len() * 24);
    for p in &m.points {
        for c in p {
            points.extend_from_slice(&c.to_ne_bytes());
        }
    }
    let single = |v: &[u32]| v.iter().map(|&x| [x]).collect::<Vec<[u32; 1]>>();
    let d = PyDict::new(py);
    d.set_item("points", PyBytes::new(py, &points))?;
    d.set_item("tets", ints(py, &m.tets))?;
    d.set_item("tet_parts", ints(py, &single(&m.tet_parts)))?;
    d.set_item("parts", m.parts)?;
    d.set_item("triangles", ints(py, &m.triangles))?;
    d.set_item("triangle_patches", ints(py, &single(&m.triangle_patches)))?;
    d.set_item("patches", m.patches)?;
    Ok(d)
}

#[pymodule]
fn _mesher(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_run, m)?)?;
    m.add_function(wrap_pyfunction!(_run_request, m)?)?;
    m.add_function(wrap_pyfunction!(_build_job, m)?)?;
    m.add_function(wrap_pyfunction!(_read_mesh, m)?)?;
    m.add_function(wrap_pyfunction!(_next_job, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
