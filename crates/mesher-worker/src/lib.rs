//! The worker: runs one job and reports through [`WorkerMessage`]s. The binary
//! wraps [`run`] with the stdin/stdout protocol; keeping the logic here lets
//! tests drive it in-process too.

use std::path::Path;

use mesher_core::protocol::{Job, REPORT_FILE, WorkerMessage};
use mesher_core::{Action, DefectCode, Outcome, Record, RecordLog, Report, Stage};

/// Run a job to completion, sending progress, records and the final result to
/// `emit`. Returns the exit code. Never panics on bad input.
pub fn run(job: &Job, emit: &mut dyn FnMut(&WorkerMessage)) -> i32 {
    let mut log = RecordLog::new();
    let mut report = Report::new(&job.input, job.config.clone());

    emit(&WorkerMessage::Progress {
        stage: Stage::Parse,
        fraction: 0.0,
    });
    let model = mesher_io::read_obj(Path::new(&job.input), &mut log, &mut |fraction| {
        emit(&WorkerMessage::Progress {
            stage: Stage::Parse,
            fraction,
        });
    });
    report.stages_run.push(Stage::Parse);
    report.input.stats = model.as_ref().map(mesher_io::ObjModel::stats);

    // Milestone 1 stops after parsing; diagnosis and meshing follow in later
    // milestones.

    for record in log.records() {
        emit(&WorkerMessage::Record {
            record: Box::new(record.clone()),
        });
    }

    let report_path = write_report(job, &mut report, log, emit);
    let exit_code = report.summary.exit_code;
    emit(&WorkerMessage::Result {
        exit_code,
        report_path,
        report: Box::new(report),
    });
    exit_code
}

/// Finish the report and write it. If writing fails, the report gains an
/// `OUT_WRITE_FAILED` record (streamed too) and the path is `None`.
fn write_report(
    job: &Job,
    report: &mut Report,
    mut log: RecordLog,
    emit: &mut dyn FnMut(&WorkerMessage),
) -> Option<String> {
    let dir = Path::new(&job.output_dir);
    let path = dir.join(REPORT_FILE);
    let attempt = std::fs::create_dir_all(dir).and_then(|()| {
        let mut finished = report.clone();
        finished.finish(log.clone());
        std::fs::write(&path, finished.to_json()).map(|()| finished)
    });
    match attempt {
        Ok(finished) => {
            *report = finished;
            Some(path.display().to_string())
        }
        Err(err) => {
            let record = Record::new(DefectCode::OutWriteFailed, Stage::Emit)
                .measured("path", path.display().to_string())
                .measured("io_error", err.kind().to_string())
                .action(Action::Abort, Outcome::Applied)
                .hint(format!("could not write the report: {err}"));
            if let Some(kept) = log.push(record) {
                emit(&WorkerMessage::Record {
                    record: Box::new(kept.clone()),
                });
            }
            report.finish(log);
            None
        }
    }
}
