//! What a front end asks for, turned into a worker [`Job`]: the one place a
//! job is assembled, so the CLI and the Python API given the same options
//! run the same job and get the same report.
//!
//! A request names an input (an OBJ or `.meshb` file, or a project file,
//! `.toml`) and options. A project supplies the input, output directory and
//! settings; options given in the request override it. Without a project or
//! an `out`, the output directory is `<input stem>_out` in the current
//! directory.

use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;

use serde::{Deserialize, Serialize};

use mesher_core::project::{Project, ProjectError};
use mesher_core::protocol::{Job, REPORT_FILE, WorkerMessage};
use mesher_core::{
    Action, CellChoice, DefectCode, JobConfig, KernelChoice, Outcome, OutputFormat, Record,
    RecordLog, Report, Stage,
};

use crate::Limits;

/// A front end's request: an input and the options that override the
/// defaults (or a project's settings). Unset options keep them.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JobRequest {
    /// An OBJ or `.meshb` file, or a project file (`.toml`).
    pub input: PathBuf,
    /// The output directory [default: the project's, else `<input stem>_out`
    /// in the current directory].
    #[serde(default)]
    pub out: Option<PathBuf>,
    /// Names or 1-based indices of the parts to mesh [default: all].
    #[serde(default)]
    pub parts: Vec<String>,
    #[serde(default)]
    pub no_repair: bool,
    #[serde(default)]
    pub strict: bool,
    #[serde(default)]
    pub uniform_size: bool,
    #[serde(default)]
    pub no_remesh: bool,
    #[serde(default)]
    pub kernel: Option<KernelChoice>,
    #[serde(default)]
    pub max_deviation: Option<f64>,
    #[serde(default)]
    pub element_size: Option<f64>,
    #[serde(default)]
    pub wrap_alpha: Option<f64>,
    #[serde(default)]
    pub feature_angle: Option<f64>,
    /// The formats written [default: the project's, else Exodus, Gmsh and
    /// VTU]; written sorted, each once.
    #[serde(default)]
    pub formats: Vec<OutputFormat>,
    /// The cells written: tets, or the polyhedral dual [default: the
    /// project's, else tets].
    #[serde(default)]
    pub cells: Option<CellChoice>,
}

/// Why a request could not become a job.
#[derive(Debug, Clone, PartialEq)]
pub enum RequestError {
    /// An option is out of range (the CLI rejects these as usage errors);
    /// nothing runs and no report is written.
    Invalid(String),
    /// The project file cannot be used: a report with `PARSE_PROJECT_INVALID`.
    Project(ProjectError),
}

impl std::fmt::Display for RequestError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RequestError::Invalid(e) => write!(f, "{e}"),
            RequestError::Project(e) => write!(f, "{e}"),
        }
    }
}

impl JobRequest {
    /// Whether the input is a project file (extension `.toml`, any case).
    #[must_use]
    pub fn is_project(&self) -> bool {
        self.input
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("toml"))
    }

    /// Check the options' ranges, as the CLI's argument parser does.
    ///
    /// # Errors
    /// `Invalid` naming the first option out of range.
    pub fn validate(&self) -> Result<(), RequestError> {
        let positive = |name: &str, v: Option<f64>| match v {
            Some(x) if !(x > 0.0 && x.is_finite()) => Err(RequestError::Invalid(format!(
                "`{name}` must be a positive number, not {x}"
            ))),
            _ => Ok(()),
        };
        positive("max_deviation", self.max_deviation)?;
        positive("element_size", self.element_size)?;
        positive("wrap_alpha", self.wrap_alpha)?;
        if let Some(a) = self.feature_angle
            && !(a > 0.0 && a < 180.0)
        {
            return Err(RequestError::Invalid(format!(
                "`feature_angle` must be between 0 and 180 degrees, not {a}"
            )));
        }
        Ok(())
    }

    /// The job this request runs.
    ///
    /// # Errors
    /// `Invalid` for an option out of range; `Project` for a project file
    /// that cannot be read or is not a project.
    pub fn job(&self) -> Result<Job, RequestError> {
        self.validate()?;
        let mut job = if self.is_project() {
            Project::read(&self.input)
                .map_err(RequestError::Project)?
                .job()
        } else {
            Job {
                input: self.input.display().to_string(),
                output_dir: default_out_dir(&self.input).display().to_string(),
                config: JobConfig::default(),
            }
        };
        if let Some(out) = &self.out {
            job.output_dir = out.display().to_string();
        }
        let c = &mut job.config;
        if !self.parts.is_empty() {
            c.parts.clone_from(&self.parts);
        }
        c.no_repair |= self.no_repair;
        c.strict |= self.strict;
        c.uniform_size |= self.uniform_size;
        c.no_remesh |= self.no_remesh;
        if let Some(k) = self.kernel {
            c.kernel = k;
        }
        if let Some(d) = self.max_deviation {
            c.max_deviation = d;
        }
        if let Some(e) = self.element_size {
            c.element_size = e;
        }
        if let Some(a) = self.wrap_alpha {
            c.wrap_alpha = a;
        }
        if let Some(a) = self.feature_angle {
            c.feature_angle = Some(a);
        }
        if !self.formats.is_empty() {
            let mut formats = self.formats.clone();
            formats.sort_unstable();
            formats.dedup();
            c.formats = formats;
        }
        if let Some(cells) = self.cells {
            c.cells = cells;
        }
        Ok(job)
    }

    /// The output directory a failed project still reports into: `out`, else
    /// `<input stem>_out`.
    fn failure_out_dir(&self) -> PathBuf {
        self.out
            .clone()
            .unwrap_or_else(|| default_out_dir(&self.input))
    }
}

/// `<input stem>_out`, relative to the current directory.
#[must_use]
pub fn default_out_dir(input: &Path) -> PathBuf {
    let stem = input
        .file_stem()
        .map_or_else(|| "mesh".to_owned(), |s| s.to_string_lossy().into_owned());
    PathBuf::from(format!("{stem}_out"))
}

/// The report for a project file that cannot be used, written as
/// `report.json` in the request's output directory; its path is `None` if it
/// could not be written.
#[must_use]
pub fn project_failure(request: &JobRequest, error: &ProjectError) -> (Report, Option<String>) {
    let mut log = RecordLog::new();
    log.push(
        Record::new(DefectCode::ParseProjectInvalid, Stage::Parse)
            .hint(error.to_string())
            .action(Action::Abort, Outcome::Applied),
    );
    let mut report = Report::new(&request.input.display().to_string(), JobConfig::default());
    report.finish(log);
    let out = request.failure_out_dir();
    let path = out.join(REPORT_FILE);
    let written = std::fs::create_dir_all(&out)
        .and_then(|()| std::fs::write(&path, report.to_json()))
        .is_ok();
    (report, written.then(|| path.display().to_string()))
}

/// Run a request: its job in the worker, or, for a project file that cannot
/// be used, the report that says so.
///
/// # Errors
/// `Invalid` for an option out of range (nothing runs, no report).
pub fn run_request(
    worker: &Path,
    request: &JobRequest,
    limits: Limits,
    on_message: &mut dyn FnMut(&WorkerMessage),
    cancel: &AtomicBool,
) -> Result<(Report, Option<String>), String> {
    match request.job() {
        Ok(job) => Ok(crate::run_job(worker, &job, limits, on_message, cancel)),
        Err(RequestError::Project(e)) => Ok(project_failure(request, &e)),
        Err(RequestError::Invalid(e)) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;

    fn request(input: &str) -> JobRequest {
        JobRequest {
            input: input.into(),
            ..JobRequest::default()
        }
    }

    #[test]
    fn a_plain_input_meshes_into_its_stem_out_with_defaults() {
        let job = request("models/part.obj").job().unwrap();
        assert_eq!(
            job.input,
            Path::new("models/part.obj").display().to_string()
        );
        assert_eq!(job.output_dir, "part_out");
        assert_eq!(job.config, JobConfig::default());
    }

    #[test]
    fn options_override_the_defaults() {
        let r = JobRequest {
            out: Some("o".into()),
            parts: vec!["lid".into(), "2".into()],
            strict: true,
            kernel: Some(KernelChoice::Exact),
            element_size: Some(0.05),
            feature_angle: Some(30.0),
            formats: vec![OutputFormat::Vtu, OutputFormat::Exodus, OutputFormat::Vtu],
            cells: Some(CellChoice::Dual),
            ..request("m.obj")
        };
        let job = r.job().unwrap();
        assert_eq!(job.output_dir, "o");
        let c = &job.config;
        assert_eq!(c.parts, ["lid", "2"]);
        assert!(c.strict && !c.no_repair);
        assert_eq!(c.kernel, KernelChoice::Exact);
        assert_eq!(c.element_size, 0.05);
        assert_eq!(c.feature_angle, Some(30.0));
        assert_eq!(c.formats, [OutputFormat::Exodus, OutputFormat::Vtu]);
        assert_eq!(c.cells, CellChoice::Dual);
    }

    #[test]
    fn options_out_of_range_are_refused() {
        for r in [
            JobRequest {
                element_size: Some(0.0),
                ..request("m.obj")
            },
            JobRequest {
                max_deviation: Some(f64::NAN),
                ..request("m.obj")
            },
            JobRequest {
                wrap_alpha: Some(-1.0),
                ..request("m.obj")
            },
            JobRequest {
                feature_angle: Some(180.0),
                ..request("m.obj")
            },
        ] {
            assert!(matches!(r.job(), Err(RequestError::Invalid(_))), "{r:?}");
        }
    }

    #[test]
    fn a_project_supplies_the_job_and_options_override_it() {
        let dir = tempfile::tempdir().unwrap();
        let project = dir.path().join("p.toml");
        std::fs::write(
            &project,
            "input = \"m.obj\"\noutput = \"o\"\n[settings]\nparts = [\"a\"]\nelement_size = 0.1\n",
        )
        .unwrap();
        let job = request(project.to_str().unwrap()).job().unwrap();
        assert_eq!(job.config.parts, ["a"]);
        assert_eq!(job.config.element_size, 0.1);
        assert_eq!(Path::new(&job.output_dir), dir.path().join("o"));
        let job = JobRequest {
            element_size: Some(0.2),
            out: Some("elsewhere".into()),
            ..request(project.to_str().unwrap())
        }
        .job()
        .unwrap();
        assert_eq!(job.config.element_size, 0.2);
        assert_eq!(job.config.parts, ["a"]);
        assert_eq!(job.output_dir, "elsewhere");
    }

    #[test]
    fn a_project_that_cannot_be_used_is_a_report() {
        let dir = tempfile::tempdir().unwrap();
        let project = dir.path().join("bad.toml");
        std::fs::write(&project, "not a project").unwrap();
        let r = JobRequest {
            out: Some(dir.path().join("out")),
            ..request(project.to_str().unwrap())
        };
        let Err(RequestError::Project(e)) = r.job() else {
            panic!("expected a project error");
        };
        let (report, path) = project_failure(&r, &e);
        assert_eq!(report.records[0].code, DefectCode::ParseProjectInvalid);
        let written = std::fs::read_to_string(path.unwrap()).unwrap();
        assert_eq!(written, report.to_json());
    }

    #[test]
    fn a_request_reads_from_json_with_defaults_and_refuses_unknown_options() {
        let r: JobRequest =
            serde_json::from_str(r#"{"input": "m.obj", "kernel": "exact", "formats": ["vtu"]}"#)
                .unwrap();
        assert_eq!(r.kernel, Some(KernelChoice::Exact));
        assert_eq!(r.formats, [OutputFormat::Vtu]);
        assert!(serde_json::from_str::<JobRequest>(r#"{"input": "m.obj", "size": 1}"#).is_err());
        assert!(serde_json::from_str::<JobRequest>(r#"{"out": "o"}"#).is_err());
    }
}
