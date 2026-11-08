//! `mesher-worker`: reads one job as JSON from stdin (or `--job <path>`), runs
//! it, and writes JSON-lines protocol messages to stdout. Its exit code is the
//! run's exit code; a malformed job exits 70, since only a buggy front end sends one.

use std::io::{Read, Write};
use std::process::ExitCode;

use mesher_core::ExitStatus;
use mesher_core::protocol::{Job, PROTOCOL_VERSION, TEST_FAULT_ENV, WorkerMessage};

const USAGE: &str = "usage: mesher-worker [--job <path|->] [--version]";

fn main() -> ExitCode {
    let mut job_source = "-".to_owned();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--job" => match args.next() {
                Some(path) => job_source = path,
                None => return usage_error(),
            },
            "--version" => {
                println!("mesher-worker {}", env!("CARGO_PKG_VERSION"));
                return ExitCode::SUCCESS;
            }
            _ => return usage_error(),
        }
    }

    let job = match read_job(&job_source) {
        Ok(job) => job,
        Err(message) => {
            eprintln!("mesher-worker: invalid job: {message}");
            return exit(ExitStatus::InternalBug.code());
        }
    };

    let mut stdout = std::io::stdout().lock();
    let mut emit = |message: &WorkerMessage| {
        // A closed pipe means the front end has gone; nothing useful to do.
        let _ = writeln!(stdout, "{}", message.to_line());
        let _ = stdout.flush();
    };
    emit(&WorkerMessage::Started {
        protocol: PROTOCOL_VERSION,
        worker_version: env!("CARGO_PKG_VERSION").to_owned(),
    });
    inject_test_fault(job.config.stuffing_only);
    exit(mesher_worker::run(&job, &mut emit))
}

fn read_job(source: &str) -> Result<Job, String> {
    let text = if source == "-" {
        let mut text = String::new();
        std::io::stdin()
            .read_to_string(&mut text)
            .map_err(|e| format!("reading stdin: {e}"))?;
        text
    } else {
        std::fs::read_to_string(source).map_err(|e| format!("reading {source}: {e}"))?
    };
    serde_json::from_str(&text).map_err(|e| e.to_string())
}

/// Deliberate failures for testing the front ends' crash and timeout handling.
/// `hang_conformal` hangs only when the conformal kernel may run, so the
/// client's stuffing-only retry gets through.
#[allow(clippy::panic)]
fn inject_test_fault(stuffing_only: bool) {
    match std::env::var(TEST_FAULT_ENV).as_deref() {
        Ok("hang_conformal") if !stuffing_only => loop {
            std::thread::sleep(std::time::Duration::from_secs(1));
        },
        Ok("panic") => panic!("{TEST_FAULT_ENV}=panic: deliberate worker panic"),
        Ok("hang") => loop {
            std::thread::sleep(std::time::Duration::from_secs(1));
        },
        Ok("alloc") => {
            // Touch 16 MiB at a time, up to 1 GiB, so a memory cap trips.
            let mut held: Vec<Vec<u8>> = Vec::new();
            for _ in 0..64 {
                held.push(vec![1u8; 16 << 20]);
            }
            eprintln!("allocated {} MiB without hitting a cap", held.len() * 16);
        }
        _ => {}
    }
}

fn usage_error() -> ExitCode {
    eprintln!("{USAGE}");
    ExitCode::from(2)
}

fn exit(code: i32) -> ExitCode {
    ExitCode::from(u8::try_from(code).unwrap_or(u8::MAX))
}
