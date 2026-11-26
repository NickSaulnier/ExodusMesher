//! `mesher`: the command-line front end. It holds no meshing logic: it builds a
//! job, runs `mesher-worker` under a timeout, prints a summary, and turns worker
//! crashes and timeouts into `MESH_WORKER_FAILED` reports.

use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::atomic::AtomicBool;
use std::time::Duration;

use clap::{Parser, Subcommand};
use mesher_client::JobRequest;
use mesher_core::{ExitStatus, KernelChoice, Outcome, OutputFormat, Report, Severity, Stage};

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
    /// Mesh an OBJ file, or replay a project file (.toml) saved by the GUI.
    Run(Box<RunArgs>),
    /// Print the JSON schema of report.json.
    Schema,
}

#[derive(clap::Args)]
struct RunArgs {
    /// Input .obj or .meshb file, or a project file (.toml). A project supplies the
    /// input, output directory and settings; options given here override it.
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
    /// Largest acceptable fallback deviation, relative to each part's
    /// bounding-box diagonal [default: 0.05].
    #[arg(long, value_parser = positive)]
    max_deviation: Option<f64>,
    /// Target element size, relative to each part's bounding-box diagonal
    /// [default: 0.03].
    #[arg(long, value_parser = positive)]
    element_size: Option<f64>,
    /// Mesh at the element size everywhere, without grading it down where
    /// the surface curves or the part is thin.
    #[arg(long)]
    uniform_size: bool,
    /// Keep each part's own surface triangulation: never remesh it (by
    /// default a part's surface is remeshed to the element size when its
    /// conformal mesh fails or its cells are beyond the quality limits).
    #[arg(long)]
    no_remesh: bool,
    /// The volume kernel tried first: `conformal` (refined to the element
    /// size), or `exact` (the input surface tetrahedralized exactly and
    /// rounded to f64, without refinement; the conformal kernel and
    /// stuffing follow if it fails) [default: conformal].
    #[arg(long, value_enum)]
    kernel: Option<KernelArg>,
    /// The wrap's α for repair rung 3, relative to the element size: a part
    /// the surface gate still rejects after local repairs and resolving its
    /// intersections, and whose stuffing fallback is rejected, is replaced
    /// as the last resort by a wrap that closes openings narrower than about
    /// 2α [default: 1].
    #[arg(long, value_parser = positive)]
    wrap_alpha: Option<f64>,
    /// The feature angle in degrees: edges where the surface turns by more
    /// than this are feature lines, kept by remeshing and listed in the
    /// report and features.vtp [default: 40 for remeshing, 45 for the size
    /// field's curvature term].
    #[arg(long, value_parser = angle)]
    feature_angle: Option<f64>,
    /// The mesh files to write, comma-separated: `exodus` (mesh.exo),
    /// `gmsh` (mesh.msh), `vtu` (mesh.vtu), `openfoam` (an OpenFOAM case in
    /// openfoam/), `fluent` (mesh_fluent.msh), `meshb` (mesh.meshb and
    /// mesh_refs.txt) [default: exodus, gmsh and vtu].
    #[arg(long, value_enum, value_delimiter = ',', num_args = 1..)]
    formats: Vec<FormatArg>,
    /// Kill the worker after this many seconds.
    #[arg(long, default_value_t = 3600.0, value_parser = positive)]
    timeout: f64,
    /// Cap the worker's memory at this many MiB (Windows and Linux).
    #[arg(long, default_value_t = 8192)]
    memory_limit: u64,
    /// Path to the mesher-worker binary [default: $MESHER_WORKER, else next to this executable].
    #[arg(long)]
    worker: Option<PathBuf>,
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum FormatArg {
    Exodus,
    Gmsh,
    Vtu,
    Openfoam,
    Fluent,
    Meshb,
}

impl FormatArg {
    fn format(self) -> OutputFormat {
        match self {
            FormatArg::Exodus => OutputFormat::Exodus,
            FormatArg::Gmsh => OutputFormat::Gmsh,
            FormatArg::Vtu => OutputFormat::Vtu,
            FormatArg::Openfoam => OutputFormat::Openfoam,
            FormatArg::Fluent => OutputFormat::Fluent,
            FormatArg::Meshb => OutputFormat::Meshb,
        }
    }
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum KernelArg {
    Conformal,
    Exact,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let code = match cli.command {
        Command::Schema => {
            print!("{}", Report::json_schema());
            0
        }
        Command::Run(args) => run(*args),
    };
    ExitCode::from(u8::try_from(code).unwrap_or(u8::MAX))
}

impl RunArgs {
    /// The request these arguments make (everything but the run's limits).
    fn request(&self) -> JobRequest {
        JobRequest {
            input: self.input.clone(),
            out: self.out.clone(),
            parts: self.parts.clone(),
            no_repair: self.no_repair,
            strict: self.strict,
            uniform_size: self.uniform_size,
            no_remesh: self.no_remesh,
            kernel: self.kernel.map(|k| match k {
                KernelArg::Conformal => KernelChoice::Conformal,
                KernelArg::Exact => KernelChoice::Exact,
            }),
            max_deviation: self.max_deviation,
            element_size: self.element_size,
            wrap_alpha: self.wrap_alpha,
            feature_angle: self.feature_angle,
            formats: self.formats.iter().map(|f| f.format()).collect(),
        }
    }
}

fn run(args: RunArgs) -> i32 {
    let worker = args
        .worker
        .clone()
        .unwrap_or_else(mesher_client::default_worker_path);
    let limits = mesher_client::Limits {
        timeout: Duration::try_from_secs_f64(args.timeout).unwrap_or(Duration::from_secs(3600)),
        memory: usize::try_from(args.memory_limit.saturating_mul(1 << 20)).ok(),
    };
    let never = AtomicBool::new(false);
    match mesher_client::run_request(&worker, &args.request(), limits, &mut |_| {}, &never) {
        Ok((report, report_path)) => {
            print_summary(&report, report_path.as_deref());
            report.summary.exit_code
        }
        // The argument parser already checks the ranges; this is a backstop.
        Err(e) => {
            eprintln!("mesher: {e}");
            2
        }
    }
}

fn positive(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(s) if s > 0.0 && s.is_finite() => Ok(s),
        _ => Err("expected a positive number".to_owned()),
    }
}

fn angle(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(a) if a > 0.0 && a < 180.0 => Ok(a),
        _ => Err("expected an angle between 0 and 180 degrees".to_owned()),
    }
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
            if let Some(m) = &p.mesh {
                println!(
                    "       {} {} tets, dihedral {:.1}°–{:.1}°, deviation {:.2}% of diagonal ({:?})",
                    if m.accepted { "meshed:" } else { "rejected:" },
                    m.tets,
                    m.min_dihedral_degrees,
                    m.max_dihedral_degrees,
                    100.0 * m.hausdorff_relative,
                    m.kernel
                );
            }
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

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::panic)]
    use super::*;

    /// The command line and the Python API (which sends a request as JSON)
    /// given the same options build the same request, and so the same job.
    #[test]
    fn arguments_and_json_build_the_same_request() {
        let cases: &[(&[&str], &str)] = &[
            (&["m.obj"], r#"{"input": "m.obj"}"#),
            (
                &["m.obj", "-o", "out", "--parts", "lid,2", "--strict"],
                r#"{"input": "m.obj", "out": "out", "parts": ["lid", "2"], "strict": true}"#,
            ),
            (
                &[
                    "m.obj",
                    "--element-size",
                    "0.05",
                    "--max-deviation",
                    "0.01",
                    "--kernel",
                    "exact",
                    "--feature-angle",
                    "30",
                    "--wrap-alpha",
                    "2",
                ],
                r#"{"input": "m.obj", "element_size": 0.05, "max_deviation": 0.01, "kernel": "exact", "feature_angle": 30.0, "wrap_alpha": 2.0}"#,
            ),
            (
                &["m.obj", "--no-repair", "--uniform-size", "--no-remesh"],
                r#"{"input": "m.obj", "no_repair": true, "uniform_size": true, "no_remesh": true}"#,
            ),
            (
                &["m.obj", "--formats", "vtu,openfoam,vtu,meshb"],
                r#"{"input": "m.obj", "formats": ["vtu", "openfoam", "vtu", "meshb"]}"#,
            ),
        ];
        for (args, json) in cases {
            let cli = Cli::try_parse_from(["mesher", "run"].iter().chain(args.iter())).unwrap();
            let Command::Run(run) = cli.command else {
                panic!("not a run");
            };
            let from_args = run.request();
            let from_json: JobRequest = serde_json::from_str(json).unwrap();
            assert_eq!(from_args, from_json, "{args:?}");
            assert_eq!(
                from_args.job().unwrap(),
                from_json.job().unwrap(),
                "{args:?}"
            );
        }
    }
}
