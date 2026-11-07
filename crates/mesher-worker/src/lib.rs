//! The worker: runs one job and reports through [`WorkerMessage`]s. The binary
//! wraps [`run`] with the stdin/stdout protocol; keeping the logic here lets
//! tests drive it in-process too.
//!
//! Pipeline so far: parse, diagnose, select, then (unless `no_repair`) repair
//! and the surface gate. Diagnosis runs if parsing produced geometry; repair
//! runs on every selected part without a fatal defect of its own.

use std::path::Path;

use mesher_core::diagnose::{Tolerances, check_part_pairs, diagnose, select_parts};
use mesher_core::model::SurfaceModel;
use mesher_core::protocol::{Job, OVERLAY_FILE, REPAIRED_FILE, REPORT_FILE, WorkerMessage};
use mesher_core::repair::repair;
use mesher_core::{Action, DefectCode, ExitStatus, Outcome, Record, RecordLog, Report, Stage};

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
    report.input.stats = model.as_ref().map(SurfaceModel::stats);

    let mut files: Vec<(&str, String)> = Vec::new();
    let mut internal_error = None;
    if let Some(model) = model.as_ref().filter(|_| !log.has_fatal()) {
        match process(job, model, &mut log, &mut report, emit) {
            Ok(outputs) => files = outputs,
            Err(message) => internal_error = Some(message),
        }
    }
    if let Some(message) = &internal_error {
        log.push(
            Record::new(DefectCode::MeshWorkerFailed, Stage::Worker)
                .measured("reason", "internal_error")
                .action(Action::Abort, Outcome::Applied)
                .hint(format!("internal error (a bug): {message}")),
        );
    }

    for record in log.records() {
        emit(&WorkerMessage::Record {
            record: Box::new(record.clone()),
        });
    }

    let report_path = write_outputs(
        job,
        &mut report,
        log,
        &files,
        internal_error.is_some(),
        emit,
    );
    let exit_code = report.summary.exit_code;
    emit(&WorkerMessage::Result {
        exit_code,
        report_path,
        report: Box::new(report),
    });
    exit_code
}

/// Diagnose, select parts, check selected pairs, then repair and gate.
/// Returns the files to write, or an internal-error message.
fn process(
    job: &Job,
    model: &SurfaceModel,
    log: &mut RecordLog,
    report: &mut Report,
    emit: &mut dyn FnMut(&WorkerMessage),
) -> Result<Vec<(&'static str, String)>, String> {
    let tol = Tolerances::default();
    emit(&WorkerMessage::Progress {
        stage: Stage::Diagnose,
        fraction: 0.0,
    });
    let diagnosis = diagnose(model, &tol, log).map_err(|e| e.0)?;
    report.stages_run.push(Stage::Diagnose);

    emit(&WorkerMessage::Progress {
        stage: Stage::Select,
        fraction: 0.0,
    });
    let selected = select_parts(&diagnosis, &job.config, log);
    report.stages_run.push(Stage::Select);
    let mut overlay = diagnosis.overlay.clone();
    if let Some(selected) = &selected {
        check_part_pairs(model, &diagnosis, selected, log, &mut overlay).map_err(|e| e.0)?;
    }
    report.parts = diagnosis.summaries(selected.as_deref().unwrap_or(&[]));

    let mut part_of = vec![0u32; model.triangles.len()];
    for (i, part) in diagnosis.parts.iter().enumerate() {
        for &t in &part.triangles {
            if let Some(slot) = part_of.get_mut(t as usize) {
                *slot = i as u32 + 1;
            }
        }
    }
    let mut files = vec![(
        OVERLAY_FILE,
        mesher_io::vtp::overlay_vtp(model, &part_of, &overlay),
    )];

    let Some(selected) = selected.filter(|_| !job.config.no_repair) else {
        return Ok(files);
    };
    emit(&WorkerMessage::Progress {
        stage: Stage::Repair,
        fraction: 0.0,
    });
    let repaired = repair(model, &diagnosis, &selected, &tol, log).map_err(|e| e.0)?;
    report.stages_run.push(Stage::Repair);
    report.stages_run.push(Stage::SurfaceGate);
    for r in &repaired {
        if let Some(p) = report.parts.get_mut(r.part) {
            p.repair = Some(r.summary);
        }
    }
    if !repaired.is_empty() {
        let surfaces: Vec<mesher_io::NamedSurface<'_>> = repaired
            .iter()
            .map(|r| (r.name.as_str(), &r.points[..], &r.faces[..]))
            .collect();
        files.push((REPAIRED_FILE, mesher_io::write_obj(&surfaces)));
    }
    Ok(files)
}

/// Write auxiliary files, then finish and write the report. Any write failure
/// becomes an `OUT_WRITE_FAILED` record (streamed too); the report path is
/// `None` if the report itself could not be written.
fn write_outputs(
    job: &Job,
    report: &mut Report,
    mut log: RecordLog,
    files: &[(&str, String)],
    internal_error: bool,
    emit: &mut dyn FnMut(&WorkerMessage),
) -> Option<String> {
    let dir = Path::new(&job.output_dir);
    let mut fail = |log: &mut RecordLog, path: &Path, err: &std::io::Error| {
        let record = Record::new(DefectCode::OutWriteFailed, Stage::Emit)
            .measured("path", path.display().to_string())
            .measured("io_error", err.kind().to_string())
            .action(Action::Abort, Outcome::Applied)
            .hint(format!("could not write {}: {err}", path.display()));
        if let Some(kept) = log.push(record) {
            emit(&WorkerMessage::Record {
                record: Box::new(kept.clone()),
            });
        }
    };
    let finish = |report: &mut Report, log: RecordLog| {
        report.finish(log);
        if internal_error {
            report.summary.exit_code = ExitStatus::InternalBug.code();
        }
    };

    if let Err(err) = std::fs::create_dir_all(dir) {
        fail(&mut log, dir, &err);
        finish(report, log);
        return None;
    }
    for (name, text) in files {
        let path = dir.join(name);
        if let Err(err) = std::fs::write(&path, text) {
            fail(&mut log, &path, &err);
        }
    }
    let path = dir.join(REPORT_FILE);
    let mut finished = report.clone();
    finish(&mut finished, log.clone());
    match std::fs::write(&path, finished.to_json()) {
        Ok(()) => {
            *report = finished;
            Some(path.display().to_string())
        }
        Err(err) => {
            fail(&mut log, &path, &err);
            finish(report, log);
            None
        }
    }
}
