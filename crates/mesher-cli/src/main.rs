//! `mesher`: the command-line front end. It holds no meshing logic: it builds a
//! job, runs `mesher-worker` under a timeout, prints a summary, and turns worker
//! crashes and timeouts into `MESH_WORKER_FAILED` reports.

mod supervise;

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Duration;

use clap::{Parser, Subcommand};
use mesher_core::protocol::{Job, REPORT_FILE};
use mesher_core::{ExitStatus, JobConfig, Outcome, Report, Severity, Stage};

#[derive(Parser)]
#[command(
    name = "mesher",
    version,
    about = "Automatic surface and volume mesher for OBJ files"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Mesh an OBJ file (Milestone 1: parse and report only).
    Run(RunArgs),
    /// Print the JSON schema of report.json.
    Schema,
}

#[derive(clap::Args)]
struct RunArgs {
    /// Input .obj file.
    input: PathBuf,
    /// Output directory [default: <input stem>_out in the current directory].
    #[arg(short, long)]
    out: Option<PathBuf>,
    /// Names or indices of the parts to mesh, comma-separated [default: all].
    #[arg(long, value_delimiter = ',')]
    parts: Vec<String>,
    /// Diagnose and report only; never repair.
    #[arg(long)]
    no_repair: bool,
    /// Conformal paths only; never use the stuffing fallback.
    #[arg(long)]
    strict: bool,
    /// Largest acceptable fallback deviation, relative to the bounding-box diagonal.
    #[arg(long, default_value_t = JobConfig::default().max_deviation)]
    max_deviation: f64,
    /// Kill the worker after this many seconds.
    #[arg(long, default_value_t = 3600.0, value_parser = positive_seconds)]
    timeout: f64,
    /// Path to the mesher-worker binary [default: $MESHER_WORKER, else next to this executable].
    #[arg(long)]
    worker: Option<PathBuf>,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let code = match cli.command {
        Command::Schema => {
            print!("{}", Report::json_schema());
            0
        }
        Command::Run(args) => run(args),
    };
    ExitCode::from(u8::try_from(code).unwrap_or(u8::MAX))
}

fn run(args: RunArgs) -> i32 {
    let out = args
        .out
        .clone()
        .unwrap_or_else(|| default_out_dir(&args.input));
    let job = Job {
        input: args.input.display().to_string(),
        output_dir: out.display().to_string(),
        config: JobConfig {
            parts: args.parts,
            no_repair: args.no_repair,
            strict: args.strict,
            max_deviation: args.max_deviation,
        },
    };
    let worker = args.worker.unwrap_or_else(default_worker_path);
    let timeout = Duration::try_from_secs_f64(args.timeout).unwrap_or(Duration::MAX);

    let (report, report_path) = supervise::run_job(&worker, &job, timeout);
    print_summary(&report, report_path.as_deref());
    report.summary.exit_code
}

fn positive_seconds(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(s) if s > 0.0 && s.is_finite() => Ok(s),
        _ => Err("expected a positive number of seconds".to_owned()),
    }
}

fn default_out_dir(input: &Path) -> PathBuf {
    let stem = input
        .file_stem()
        .map_or_else(|| "mesh".to_owned(), |s| s.to_string_lossy().into_owned());
    PathBuf::from(format!("{stem}_out"))
}

fn default_worker_path() -> PathBuf {
    if let Some(path) = std::env::var_os("MESHER_WORKER") {
        return PathBuf::from(path);
    }
    let name = format!("mesher-worker{}", std::env::consts::EXE_SUFFIX);
    std::env::current_exe()
        .ok()
        .and_then(|exe| exe.parent().map(|dir| dir.join(&name)))
        .unwrap_or_else(|| PathBuf::from(name))
}

fn print_summary(report: &Report, report_path: Option<&str>) {
    let s = &report.summary;
    println!("mesher: {}", report.input.path);
    if let Some(stats) = &report.input.stats {
        let mut line = format!(
            "  input: {} vertices, {} polygons, {} triangles",
            stats.vertices, stats.polygons, stats.triangles
        );
        if !stats.objects.is_empty() {
            line.push_str(&format!("; objects: {}", stats.objects.join(", ")));
        }
        println!("{line}");
    }
    if !report.parts.is_empty() {
        println!(
            "  parts ({} selected of {}):",
            report.parts.iter().filter(|p| p.selected).count(),
            report.parts.len()
        );
        for p in &report.parts {
            let repair = p.repair.map_or_else(String::new, |r| {
                format!(
                    "; repaired to {} faces, gate {}",
                    r.faces_after,
                    if r.gate_passed { "passed" } else { "failed" }
                )
            });
            println!(
                "    {}{} {}: {} triangles, {} shell{}{repair}",
                if p.selected { "*" } else { " " },
                p.index,
                p.name,
                p.triangles,
                p.shells,
                if p.shells == 1 { "" } else { "s" }
            );
        }
    }
    println!(
        "  findings: {} info, {} warnings, {} errors, {} fatal",
        s.info, s.warnings, s.errors, s.fatal
    );
    // Findings per code, then what repair and the gate made of them. Counts come
    // from the detailed records (the summary's per-code totals mix stages).
    for (code, count) in &s.by_code {
        let of_code = || report.records.iter().filter(|r| r.code == *code);
        let found: Vec<_> = of_code()
            .filter(|r| !matches!(r.stage, Stage::Repair | Stage::SurfaceGate))
            .collect();
        let tally = |stage: Stage, outcome: Option<Outcome>| {
            of_code()
                .filter(|r| r.stage == stage && (outcome.is_none() || r.outcome == outcome))
                .count()
        };
        let mut after = Vec::new();
        for (n, label) in [
            (tally(Stage::Repair, Some(Outcome::Applied)), "repaired"),
            (tally(Stage::Repair, Some(Outcome::Failed)), "repair failed"),
            (tally(Stage::SurfaceGate, None), "still failing the gate"),
        ] {
            if n > 0 {
                after.push(format!("{label} {n}"));
            }
        }
        let first = found.first().copied().or_else(|| of_code().next());
        let severity = first.map_or(code.default_severity(), |r| r.severity);
        let at = first
            .and_then(|r| r.location.obj_lines.first())
            .map(|l| format!(", first at line {l}"))
            .unwrap_or_default();
        let truncated = if count.reported < count.total {
            "+"
        } else {
            ""
        };
        let after = if after.is_empty() {
            String::new()
        } else {
            format!(" ({})", after.join(", "))
        };
        println!(
            "    {:<7} {} x{}{truncated}{after}{at}",
            severity_label(severity),
            code.as_str(),
            found.len()
        );
    }
    let stages: Vec<String> = report
        .stages_run
        .iter()
        .filter_map(|st| serde_json::to_value(st).ok())
        .filter_map(|v| v.as_str().map(str::to_owned))
        .collect();
    if !stages.is_empty() {
        println!("  stages run: {}", stages.join(", "));
    }
    match report_path {
        Some(path) => println!("  report: {path}"),
        None => println!("  report: not written"),
    }
    println!("  exit {}: {}", s.exit_code, exit_label(s.exit_code));
}

fn severity_label(severity: Severity) -> &'static str {
    match severity {
        Severity::Info => "info",
        Severity::Warning => "warning",
        Severity::Error => "error",
        Severity::Fatal => "fatal",
    }
}

fn exit_label(code: i32) -> &'static str {
    match ExitStatus::from_code(code) {
        Some(ExitStatus::Clean) => "clean",
        Some(ExitStatus::Repaired) => "succeeded after repairs",
        Some(ExitStatus::Fallback) => "succeeded via fallback",
        Some(ExitStatus::Failed) => "failed; see report",
        Some(ExitStatus::InternalBug) => "internal error; please file a bug with the report",
        None => "unknown",
    }
}

/// Where the CLI writes a report the worker could not.
pub(crate) fn report_file(job: &Job) -> PathBuf {
    Path::new(&job.output_dir).join(REPORT_FILE)
}
