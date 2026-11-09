//! `mesher-client`: runs a job in a `mesher-worker` process for a front end
//! (the CLI or the GUI). It sends the job, streams protocol messages to the
//! caller as they arrive, enforces the timeout and memory cap
//! (`process_control`: a job object on Windows, `setrlimit` on Linux), honours
//! cancellation, and builds a `MESH_WORKER_FAILED` report when the worker
//! does not deliver one. A worker that timed out or ran out of memory gets
//! one retry with the stuffing fallback only (unless the job is strict), so
//! a slow conformal kernel still ends in a mesh. No meshing logic lives here.

use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;
use std::time::{Duration, Instant};

use process_control::{ChildExt, Control};

use mesher_core::protocol::{Job, REPORT_FILE, WorkerMessage};
use mesher_core::{
    Action, Attempt, AttemptOutcome, DefectCode, ExitStatus, Kernel, Outcome, Record, RecordLog,
    Report, Severity, Stage,
};
use serde_json::Value;

/// Keep at most this much of the worker's stderr for the report.
const STDERR_TAIL_BYTES: usize = 4096;

/// How often the supervisor wakes to deliver messages and check for cancel.
const POLL: Duration = Duration::from_millis(50);

enum Event {
    Line(String),
    Eof,
}

/// How the worker run ended, when it did not deliver a result.
enum Failure {
    Spawn(std::io::Error),
    Timeout,
    Cancelled,
    OutOfMemory,
    Died {
        status: Option<i64>,
    },
    /// The worker sent a result line that does not parse (a protocol bug).
    UnreadableResult {
        error: String,
    },
}

/// Limits the worker runs under.
#[derive(Debug, Clone, Copy)]
pub struct Limits {
    pub timeout: Duration,
    /// Memory cap in bytes (enforced on Windows and Linux).
    pub memory: Option<usize>,
}

/// Where `report.json` goes for a job.
#[must_use]
pub fn report_file(job: &Job) -> PathBuf {
    Path::new(&job.output_dir).join(REPORT_FILE)
}

/// The worker binary: `$MESHER_WORKER`, else `mesher-worker` next to the
/// running executable.
#[must_use]
pub fn default_worker_path() -> PathBuf {
    if let Some(path) = std::env::var_os("MESHER_WORKER") {
        return PathBuf::from(path);
    }
    let name = format!("mesher-worker{}", std::env::consts::EXE_SUFFIX);
    std::env::current_exe()
        .ok()
        .and_then(|exe| exe.parent().map(|dir| dir.join(&name)))
        .unwrap_or_else(|| PathBuf::from(name))
}

/// Run the job in a worker process. `on_message` sees every protocol message
/// as it arrives; setting `cancel` kills the worker. Always returns a report,
/// plus the path it was written to (`None` if it could not be written).
///
/// If the worker times out or runs out of memory while the conformal kernel
/// may run, it is run once more with `stuffing_only` (same limits). The
/// retry's report then leads with the `MESH_WORKER_FAILED` record of the
/// first run: a warning when the fallback produced a mesh, an error when not.
pub fn run_job(
    worker: &Path,
    job: &Job,
    limits: Limits,
    on_message: &mut dyn FnMut(&WorkerMessage),
    cancel: &AtomicBool,
) -> (Report, Option<String>) {
    let mut streamed = Vec::new();
    let (failure, stderr) =
        match spawn_and_read(worker, job, limits, on_message, cancel, &mut streamed) {
            Ok(result) => return result,
            Err(e) => e,
        };
    let retry = matches!(failure, Failure::Timeout | Failure::OutOfMemory)
        && !job.config.strict
        && !job.config.no_repair
        && !job.config.stuffing_only;
    if !retry {
        return failure_report(job, streamed, failure, &stderr, limits);
    }
    let mut fallback = job.clone();
    fallback.config.stuffing_only = true;
    let mut again = Vec::new();
    match spawn_and_read(worker, &fallback, limits, on_message, cancel, &mut again) {
        Ok((mut report, path)) => {
            let (record, _) = failure_record(&failure, &stderr, limits);
            let meshed = report.summary.exit_code <= ExitStatus::Fallback.code();
            report.attempts.insert(
                0,
                Attempt {
                    part: None,
                    kernel: Kernel::Delaunay,
                    outcome: AttemptOutcome::Failed,
                    remeshed: false,
                    stdout: String::new(),
                    stderr: stderr.clone(),
                },
            );
            report.prepend(if meshed {
                record
                    .severity(Severity::Warning)
                    .action(Action::Fallback, Outcome::Applied)
            } else {
                record.action(Action::Fallback, Outcome::Failed)
            });
            let written = path
                .as_ref()
                .is_some_and(|p| std::fs::write(p, report.to_json()).is_ok());
            (report, path.filter(|_| written))
        }
        // The retry failed too: report the first failure, as without a retry.
        Err(_) => failure_report(job, streamed, failure, &stderr, limits),
    }
}

type Delivered = (Report, Option<String>);

fn spawn_and_read(
    worker: &Path,
    job: &Job,
    limits: Limits,
    on_message: &mut dyn FnMut(&WorkerMessage),
    cancel: &AtomicBool,
    streamed: &mut Vec<Record>,
) -> Result<Delivered, (Failure, String)> {
    let mut child = Command::new(worker)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| (Failure::Spawn(e), String::new()))?;

    if let Some(mut stdin) = child.stdin.take() {
        // A worker that died before reading its job shows up below as a death.
        let _ = stdin.write_all(serde_json::to_string(job).unwrap_or_default().as_bytes());
    }

    let (tx, rx) = mpsc::channel();
    if let Some(stdout) = child.stdout.take() {
        std::thread::spawn(move || {
            for line in BufReader::new(stdout).lines() {
                match line {
                    Ok(line) => {
                        if tx.send(Event::Line(line)).is_err() {
                            return;
                        }
                    }
                    Err(_) => break,
                }
            }
            let _ = tx.send(Event::Eof);
        });
    }
    let stderr_thread = child.stderr.take().map(|mut stderr| {
        std::thread::spawn(move || {
            let mut buf = Vec::new();
            let _ = stderr.read_to_end(&mut buf);
            buf
        })
    });

    let mut delivered = None;
    let mut eof = false;
    let mut unreadable: Option<String> = None;
    let mut handle = |event: Event, delivered: &mut Option<Delivered>, eof: &mut bool| match event {
        // A stray non-protocol line surfaces only if no result follows.
        Event::Line(line) => match WorkerMessage::from_line(&line) {
            Ok(message) => {
                on_message(&message);
                match message {
                    WorkerMessage::Record { record } => streamed.push(*record),
                    WorkerMessage::Result {
                        report,
                        report_path,
                        ..
                    } => *delivered = Some((*report, report_path)),
                    _ => {}
                }
            }
            Err(e) if line.starts_with(r#"{"type":"result""#) => {
                unreadable = Some(e.to_string());
            }
            Err(_) => {}
        },
        Event::Eof => *eof = true,
    };

    // Wait in short slices so messages reach the caller as they arrive and
    // cancellation is noticed; the memory cap is set on the first slice.
    let start = Instant::now();
    let mut first = true;
    let (status, ended) = loop {
        let slice = POLL.min(limits.timeout.saturating_sub(start.elapsed()));
        let waited = wait_slice(&mut child, slice, if first { limits.memory } else { None });
        first = false;
        while let Ok(event) = rx.try_recv() {
            handle(event, &mut delivered, &mut eof);
        }
        match waited {
            Ok(Some(code)) => break (code, None),
            Ok(None) => {
                if cancel.load(Ordering::Relaxed) {
                    let _ = child.kill();
                    let _ = child.wait();
                    break (None, Some(Failure::Cancelled));
                }
                if start.elapsed() >= limits.timeout {
                    let _ = child.kill();
                    let _ = child.wait();
                    break (None, Some(Failure::Timeout));
                }
            }
            Err(()) => break (None, None),
        }
    };
    while !eof {
        match rx.recv() {
            Ok(event) => handle(event, &mut delivered, &mut eof),
            Err(_) => break,
        }
    }
    let stderr = stderr_thread
        .and_then(|t| t.join().ok())
        .map(|bytes| tail(&String::from_utf8_lossy(&bytes)))
        .unwrap_or_default();

    match (delivered, ended) {
        (_, Some(failure)) => Err((failure, stderr)),
        (Some(result), None) => Ok(result),
        (None, None) if unreadable.is_some() => Err((
            Failure::UnreadableResult {
                error: unreadable.unwrap_or_default(),
            },
            stderr,
        )),
        (None, None) if stderr.contains("memory allocation of") => {
            Err((Failure::OutOfMemory, stderr))
        }
        (None, None) => Err((Failure::Died { status }, stderr)),
    }
}

/// Wait up to `slice` for the child: `Ok(Some(code))` when it exited,
/// `Ok(None)` when it is still running, `Err` when waiting failed.
#[allow(clippy::option_option)]
fn wait_slice(
    child: &mut Child,
    slice: Duration,
    memory: Option<usize>,
) -> Result<Option<Option<i64>>, ()> {
    let control = child.controlled().time_limit(slice);
    #[cfg(any(
        windows,
        target_os = "android",
        all(target_os = "linux", any(target_env = "gnu", target_env = "musl"))
    ))]
    let control = match memory {
        Some(bytes) => control.memory_limit(bytes),
        None => control,
    };
    #[cfg(not(any(
        windows,
        target_os = "android",
        all(target_os = "linux", any(target_env = "gnu", target_env = "musl"))
    )))]
    let _ = memory;
    match control.wait() {
        Ok(Some(status)) => Ok(Some(status.code())),
        Ok(None) => Ok(None),
        Err(_) => Err(()),
    }
}

fn tail(text: &str) -> String {
    let start = text.len().saturating_sub(STDERR_TAIL_BYTES);
    let start = (start..=text.len())
        .find(|&i| text.is_char_boundary(i))
        .unwrap_or(text.len());
    text[start..].trim_end().to_owned()
}

fn failure_report(
    job: &Job,
    streamed: Vec<Record>,
    failure: Failure,
    stderr: &str,
    limits: Limits,
) -> (Report, Option<String>) {
    let mut log = RecordLog::new();
    for record in streamed {
        log.push(record);
    }
    let (record, internal_bug) = failure_record(&failure, stderr, limits);
    log.push(record.action(Action::Abort, Outcome::Applied));

    let mut report = Report::new(&job.input, job.config.clone());
    report.finish(log);
    if internal_bug {
        report.summary.exit_code = ExitStatus::InternalBug.code();
    }

    let path = report_file(job);
    let written = std::fs::create_dir_all(&job.output_dir)
        .and_then(|()| std::fs::write(&path, report.to_json()))
        .is_ok();
    (report, written.then(|| path.display().to_string()))
}

/// The `MESH_WORKER_FAILED` record for a failed run, and whether the failure
/// is a bug (the worker could not start, or died without a result).
fn failure_record(failure: &Failure, stderr: &str, limits: Limits) -> (Record, bool) {
    let mut record = Record::new(DefectCode::MeshWorkerFailed, Stage::Worker);
    let internal_bug = match failure {
        Failure::Spawn(err) => {
            record = record
                .measured("reason", "spawn_failed")
                .measured("io_error", err.to_string())
                .hint("could not start mesher-worker; pass --worker or set MESHER_WORKER");
            true
        }
        Failure::Timeout => {
            record = record
                .measured("reason", "timeout")
                .threshold("timeout_seconds", limits.timeout.as_secs_f64())
                .hint("the worker was killed after the timeout");
            false
        }
        Failure::Cancelled => {
            record = record
                .measured("reason", "cancelled")
                .hint("the run was cancelled");
            false
        }
        Failure::OutOfMemory => {
            record = record
                .measured("reason", "out_of_memory")
                .hint("the worker hit its memory cap; raise --memory-limit or --element-size");
            if let Some(bytes) = limits.memory {
                record = record.threshold("memory_limit_bytes", bytes as u64);
            }
            false
        }
        Failure::Died { status } => {
            record = record
                .measured("reason", "died")
                .measured("exit_status", status.map_or(Value::Null, Value::from))
                .hint("the worker exited without a result; this is a bug");
            true
        }
        Failure::UnreadableResult { error } => {
            record = record
                .measured("reason", "unreadable_result")
                .measured("parse_error", error.as_str())
                .hint("the worker finished, but its result could not be read; this is a bug");
            true
        }
    };
    if !stderr.is_empty() {
        record = record.measured("stderr_tail", stderr);
    }
    (record, internal_bug)
}

#[cfg(test)]
mod tests {
    use super::tail;

    #[test]
    fn tail_respects_char_boundaries() {
        let text = "é".repeat(super::STDERR_TAIL_BYTES);
        let t = tail(&text);
        assert!(t.len() <= super::STDERR_TAIL_BYTES);
        assert!(t.chars().all(|c| c == 'é'));
    }
}
