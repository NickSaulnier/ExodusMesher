//! Runs the worker process: sends the job, reads protocol lines, enforces the
//! timeout, and builds a `MESH_WORKER_FAILED` report when the worker dies.

use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::time::{Duration, Instant};

use mesher_core::protocol::{Job, WorkerMessage};
use mesher_core::{Action, DefectCode, ExitStatus, Outcome, Record, RecordLog, Report, Stage};
use serde_json::Value;

/// Keep at most this much of the worker's stderr for the report.
const STDERR_TAIL_BYTES: usize = 4096;

enum Event {
    Line(String),
    Eof,
}

/// How the worker run ended, when it did not deliver a result.
enum Failure {
    Spawn(std::io::Error),
    Timeout,
    Died { status: Option<i32> },
}

/// Run the job in a worker process. Always returns a report, plus the path it
/// was written to (`None` if it could not be written).
pub fn run_job(worker: &Path, job: &Job, timeout: Duration) -> (Report, Option<String>) {
    let mut streamed = Vec::new();
    let (failure, stderr) = match spawn_and_read(worker, job, timeout, &mut streamed) {
        Ok(result) => return result,
        Err(failure) => failure,
    };
    failure_report(job, streamed, failure, &stderr, timeout)
}

type Delivered = (Report, Option<String>);

fn spawn_and_read(
    worker: &Path,
    job: &Job,
    timeout: Duration,
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

    let deadline = Instant::now().checked_add(timeout);
    let mut delivered = None;
    let mut timed_out = false;
    loop {
        let event = match deadline {
            Some(deadline) => rx.recv_timeout(deadline.saturating_duration_since(Instant::now())),
            None => rx.recv().map_err(|_| mpsc::RecvTimeoutError::Disconnected),
        };
        match event {
            Ok(Event::Line(line)) => match WorkerMessage::from_line(&line) {
                Ok(WorkerMessage::Record { record }) => streamed.push(*record),
                Ok(WorkerMessage::Result {
                    report,
                    report_path,
                    ..
                }) => delivered = Some((*report, report_path)),
                // Started and progress are informational; a stray non-protocol
                // line is ignored here and surfaces only if no result follows.
                Ok(_) | Err(_) => {}
            },
            Ok(Event::Eof) | Err(mpsc::RecvTimeoutError::Disconnected) => break,
            Err(mpsc::RecvTimeoutError::Timeout) => {
                timed_out = true;
                let _ = child.kill();
                break;
            }
        }
    }
    let status = child.wait().ok().and_then(|s| s.code());
    let stderr = stderr_thread
        .and_then(|t| t.join().ok())
        .map(|bytes| tail(&String::from_utf8_lossy(&bytes)))
        .unwrap_or_default();

    match delivered {
        Some(result) if !timed_out => Ok(result),
        _ if timed_out => Err((Failure::Timeout, stderr)),
        _ => Err((Failure::Died { status }, stderr)),
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
    timeout: Duration,
) -> (Report, Option<String>) {
    let mut log = RecordLog::new();
    for record in streamed {
        log.push(record);
    }
    let mut record = Record::new(DefectCode::MeshWorkerFailed, Stage::Worker);
    let internal_bug = match &failure {
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
                .threshold("timeout_seconds", timeout.as_secs_f64())
                .hint("the worker was killed after the timeout");
            false
        }
        Failure::Died { status } => {
            record = record
                .measured("reason", "died")
                .measured("exit_status", status.map_or(Value::Null, Value::from))
                .hint("the worker exited without a result; this is a bug");
            true
        }
    };
    if !stderr.is_empty() {
        record = record.measured("stderr_tail", stderr);
    }
    log.push(record.action(Action::Abort, Outcome::Applied));

    let mut report = Report::new(&job.input, job.config.clone());
    report.finish(log);
    if internal_bug {
        report.summary.exit_code = ExitStatus::InternalBug.code();
    }

    let path = crate::report_file(job);
    let written = std::fs::create_dir_all(&job.output_dir)
        .and_then(|()| std::fs::write(&path, report.to_json()))
        .is_ok();
    (report, written.then(|| path.display().to_string()))
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
