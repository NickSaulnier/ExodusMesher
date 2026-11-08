//! Milestone 8 exit criterion: the clean Thingi10k subset (ids in
//! `tests/thingi10k/subset.txt`), meshed end to end by the `mesher` CLI and,
//! for comparison, by TetGen on the same surface with the same targets
//! (radius-edge 2 and the volume of a regular tet of the element size).
//!
//! The meshes are not in the repository (Thingiverse licences); download
//! them to `target/thingi10k/raw/<id>.stl` (or `$THINGI10K_DIR/raw`), build
//! the binaries, then run:
//!
//! ```text
//! cargo build --release -p mesher-cli -p mesher-worker
//! cargo test --release -p mesher-oracle --test thingi10k -- --ignored --nocapture
//! ```
//!
//! The table goes to `target/tmp/oracle/thingi10k.txt`.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use mesher_core::{Kernel, Report, Severity};
use mesher_geom::predicates::orient3d;
use mesher_geom::{Point3, Sign};
use mesher_oracle::{Surface, read_stl, to_obj};
use mesher_tet::{FvQuality, fv_quality};
use tritet::Tetgen;

fn data_dir() -> PathBuf {
    std::env::var_os("THINGI10K_DIR").map_or_else(
        || PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/thingi10k"),
        PathBuf::from,
    )
}

/// `mesher` and `mesher-worker` from the same profile as this test.
fn binaries() -> (PathBuf, PathBuf) {
    let exe = std::env::current_exe().unwrap();
    let dir = exe.parent().and_then(Path::parent).unwrap().to_path_buf();
    let name = |n: &str| dir.join(format!("{n}{}", std::env::consts::EXE_SUFFIX));
    let (cli, worker) = (name("mesher"), name("mesher-worker"));
    assert!(
        cli.exists() && worker.exists(),
        "build them first: cargo build [--release] -p mesher-cli -p mesher-worker (looked in {})",
        dir.display()
    );
    (cli, worker)
}

fn subset() -> Vec<u64> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../thingi10k/subset.txt");
    let text = std::fs::read_to_string(path).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.split_whitespace().next().unwrap().parse().unwrap())
        .collect()
}

fn tetgen(s: &Surface, max_volume: f64) -> Result<(Vec<Point3>, Vec<[u32; 4]>), String> {
    let mut tg = Tetgen::new(s.points.len(), Some(vec![3; s.triangles.len()]), None, None)
        .map_err(str::to_owned)?;
    for (i, q) in s.points.iter().enumerate() {
        tg.set_point(i, 0, q[0], q[1], q[2])
            .map_err(str::to_owned)?;
    }
    for (i, t) in s.triangles.iter().enumerate() {
        for (m, &v) in t.iter().enumerate() {
            tg.set_facet_point(i, m, v as usize)
                .map_err(str::to_owned)?;
        }
    }
    tg.generate_mesh(false, false, Some(max_volume), None)
        .map_err(str::to_owned)?;
    let points: Vec<Point3> = (0..tg.out_npoint())
        .map(|i| [tg.out_point(i, 0), tg.out_point(i, 1), tg.out_point(i, 2)])
        .collect();
    let tets = (0..tg.out_ncell())
        .map(|c| {
            let mut t = [0, 1, 2, 3].map(|m| tg.out_cell_point(c, m) as u32);
            let q = t.map(|v| points[v as usize]);
            if orient3d(q[0], q[1], q[2], q[3]) == Sign::Negative {
                t.swap(0, 1);
            }
            t
        })
        .collect();
    Ok((points, tets))
}

/// Aggregates for one mesher.
#[derive(Default)]
struct Totals {
    histogram: Vec<u64>,
    min_dihedral: Vec<f64>,
    mean_non_orthogonality: Vec<f64>,
    seconds: f64,
}

impl Totals {
    fn add(&mut self, histogram: &[u64], min: f64, non_orthogonality: f64, seconds: f64) {
        self.histogram
            .resize(histogram.len().max(self.histogram.len()), 0);
        for (a, b) in self.histogram.iter_mut().zip(histogram) {
            *a += b;
        }
        self.min_dihedral.push(min);
        self.mean_non_orthogonality.push(non_orthogonality);
        self.seconds += seconds;
    }

    fn describe(&self, who: &str, out: &mut String) {
        writeln!(
            out,
            "{who}: {} meshes; per-model minimum dihedral: worst {:.2}°, 5th percentile {:.2}°, median {:.2}°; median mean non-orthogonality {:.1}°; total {:.1} s",
            self.min_dihedral.len(),
            percentile(&self.min_dihedral, 0.0),
            percentile(&self.min_dihedral, 0.05),
            percentile(&self.min_dihedral, 0.5),
            percentile(&self.mean_non_orthogonality, 0.5),
            self.seconds
        )
        .unwrap();
        let total: u64 = self.histogram.iter().sum();
        let bins: Vec<String> = self
            .histogram
            .iter()
            .map(|&n| format!("{:4.1}", 100.0 * n as f64 / total.max(1) as f64))
            .collect();
        writeln!(
            out,
            "{who}: all dihedral angles, % per 10° bin: {}",
            bins.join(" ")
        )
        .unwrap();
    }
}

fn percentile(v: &[f64], p: f64) -> f64 {
    let mut v = v.to_vec();
    v.sort_by(f64::total_cmp);
    let i = (v.len().saturating_sub(1) as f64 * p).round() as usize;
    v.get(i).copied().unwrap_or(f64::NAN)
}

#[test]
#[ignore = "needs the downloaded Thingi10k subset and built binaries"]
fn clean_subset_meshes_conformally() {
    let (cli, worker) = binaries();
    let dir = data_dir();
    let (obj_dir, out_dir) = (dir.join("obj"), dir.join("out"));
    std::fs::create_dir_all(&obj_dir).unwrap();
    let ids = subset();
    let (mut ours, mut theirs) = (Totals::default(), Totals::default());
    let mut failures = Vec::new();
    let mut conformal = 0;
    let mut table = format!(
        "{:>8} {:>6} {:>4} {:>8} {:>9} | {:>8} {:>12} {:>8} {:>8} | {:>8} {:>12} {:>8} {:>8}\n",
        "id",
        "faces",
        "exit",
        "kernel",
        "size",
        "tets",
        "dihedral",
        "non-orth",
        "time",
        "TetGen",
        "dihedral",
        "non-orth",
        "time"
    );
    for &id in &ids {
        let stl = dir.join("raw").join(format!("{id}.stl"));
        let Ok(bytes) = std::fs::read(&stl) else {
            failures.push(format!("{id}: {} missing", stl.display()));
            continue;
        };
        let surface = match read_stl(&bytes) {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!("{id}: unreadable STL: {e}"));
                continue;
            }
        };
        let obj = obj_dir.join(format!("{id}.obj"));
        std::fs::write(&obj, to_obj(&surface, &format!("thing_{id}"))).unwrap();
        let out = out_dir.join(id.to_string());
        let _ = std::fs::remove_dir_all(&out);

        let started = Instant::now();
        let status = Command::new(&cli)
            .arg("run")
            .arg(&obj)
            .arg("--out")
            .arg(&out)
            .arg("--worker")
            .arg(&worker)
            .args(["--timeout", "600"])
            .output()
            .unwrap();
        let ours_seconds = started.elapsed().as_secs_f64();
        let code = status.status.code().unwrap_or(-1);
        let report: Option<Report> = std::fs::read(out.join("report.json"))
            .ok()
            .and_then(|b| serde_json::from_slice(&b).ok());
        let Some(report) = report else {
            failures.push(format!("{id}: exit {code}, no report"));
            continue;
        };
        let meshes: Vec<_> = report.parts.iter().filter_map(|p| p.mesh.clone()).collect();
        if !meshes.is_empty()
            && meshes
                .iter()
                .all(|m| m.accepted && m.kernel == Kernel::Delaunay)
        {
            conformal += 1;
        } else {
            let codes: Vec<&str> = report
                .records
                .iter()
                .filter(|r| r.severity != Severity::Info)
                .map(|r| r.code.as_str())
                .collect();
            failures.push(format!(
                "{id}: exit {code}, not conformal; records {codes:?}"
            ));
        }
        let Some(m) = meshes.first() else {
            writeln!(
                table,
                "{id:>8} {:>6} {code:>4} (no mesh)",
                surface.triangles.len()
            )
            .unwrap();
            continue;
        };
        let h = m.element_size;
        ours.add(
            &m.dihedral_histogram,
            m.min_dihedral_degrees,
            m.mean_non_orthogonality_degrees,
            ours_seconds,
        );

        let started = Instant::now();
        let reference = tetgen_in_child(id, h);
        let theirs_seconds = started.elapsed().as_secs_f64();
        let them = match reference {
            Ok(q) => {
                theirs.add(&q.histogram, q.min, q.non_orthogonality, theirs_seconds);
                format!(
                    "{:>8} {:>5.1}..{:>5.1} {:>8.1} {theirs_seconds:>7.2}s",
                    q.cells, q.min, q.max, q.non_orthogonality
                )
            }
            Err(e) => format!("failed: {e}"),
        };
        writeln!(
            table,
            "{id:>8} {:>6} {code:>4} {:>8} {h:>9.3e} | {:>8} {:>5.1}..{:>5.1} {:>8.1} {ours_seconds:>7.2}s | {them}",
            surface.triangles.len(),
            format!("{:?}", m.kernel),
            m.tets,
            m.min_dihedral_degrees,
            m.max_dihedral_degrees,
            m.mean_non_orthogonality_degrees,
        )
        .unwrap();
        eprintln!("{}", table.lines().last().unwrap_or(""));
    }

    let mut summary = format!(
        "Thingi10k clean subset: {} models, {conformal} meshed conformally by the mesher CLI at default settings.\n\
         TetGen: pzAfq on the same welded surface, with -a the regular-tet volume of the same element size.\n\
         Our times include parsing, diagnosis, repair and writing every output file.\n",
        ids.len()
    );
    ours.describe("ours  ", &mut summary);
    theirs.describe("TetGen", &mut summary);
    if !failures.is_empty() {
        writeln!(
            summary,
            "\nNot meshed conformally:\n{}",
            failures.join("\n")
        )
        .unwrap();
    }
    eprintln!("\n{summary}");
    let tmp = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("oracle");
    std::fs::create_dir_all(&tmp).unwrap();
    std::fs::write(tmp.join("thingi10k.txt"), format!("{summary}\n{table}")).unwrap();
    assert!(
        failures.is_empty(),
        "{} of {} models failed",
        failures.len(),
        ids.len()
    );
}

/// TetGen's quality summary, as the child process reports it.
struct Reference {
    cells: usize,
    min: f64,
    max: f64,
    non_orthogonality: f64,
    histogram: Vec<u64>,
}

/// Seconds TetGen gets per model before it is killed.
const TETGEN_SECONDS: u64 = 120;

/// Run TetGen on one model in a child copy of this test binary, so that a
/// model TetGen loops on cannot stall the benchmark.
fn tetgen_in_child(id: u64, h: f64) -> Result<Reference, String> {
    let mut child = Command::new(std::env::current_exe().unwrap())
        .args([
            "tetgen_one",
            "--exact",
            "--ignored",
            "--nocapture",
            "--test-threads=1",
        ])
        .env("THINGI10K_TETGEN_ID", id.to_string())
        .env("THINGI10K_TETGEN_SIZE", format!("{h:e}"))
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null())
        .spawn()
        .unwrap();
    let started = Instant::now();
    loop {
        if child.try_wait().unwrap().is_some() {
            break;
        }
        if started.elapsed().as_secs() >= TETGEN_SECONDS {
            let _ = child.kill();
            let _ = child.wait();
            return Err(format!("timed out after {TETGEN_SECONDS} s"));
        }
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
    let output = child.wait_with_output().unwrap();
    let text = String::from_utf8_lossy(&output.stdout);
    let line = text
        .lines()
        .find_map(|l| l.find("TETGEN ").map(|i| &l[i + 7..]))
        .ok_or_else(|| format!("no result (exit {:?})", output.status.code()))?;
    if let Some(e) = line.strip_prefix("error ") {
        return Err(e.to_owned());
    }
    let f: Vec<&str> = line.split_whitespace().collect();
    Ok(Reference {
        cells: f[0].parse().unwrap(),
        min: f[1].parse().unwrap(),
        max: f[2].parse().unwrap(),
        non_orthogonality: f[3].parse().unwrap(),
        histogram: f[4].split(',').map(|n| n.parse().unwrap()).collect(),
    })
}

/// The child side of `tetgen_in_child`; does nothing when run directly.
#[test]
#[ignore = "run by clean_subset_meshes_conformally"]
fn tetgen_one() {
    let (Ok(id), Ok(size)) = (
        std::env::var("THINGI10K_TETGEN_ID"),
        std::env::var("THINGI10K_TETGEN_SIZE"),
    ) else {
        return;
    };
    let h: f64 = size.parse().unwrap();
    let (surface, _) = load_surface(&data_dir(), id.parse().unwrap()).unwrap();
    let max_volume = 8.0 / (9.0 * 3f64.sqrt()) * h.powi(3);
    match tetgen(&surface, max_volume) {
        Ok((p, t)) => {
            let q: FvQuality = fv_quality(&p, &t);
            let hist: Vec<String> = q.dihedral_histogram.iter().map(u64::to_string).collect();
            println!(
                "TETGEN {} {} {} {} {}",
                q.cells,
                q.min_dihedral_degrees,
                q.max_dihedral_degrees,
                q.mean_non_orthogonality_degrees,
                hist.join(",")
            );
        }
        Err(e) => println!("TETGEN error {e}"),
    }
}

/// A downloaded model as a welded triangle surface: `raw/<id>.stl`, or for
/// the few models the dataset has only as OBJ, `raw/<id>.obj` (also returned,
/// so the CLI can read the original).
fn load_surface(dir: &Path, id: u64) -> Result<(Surface, Option<PathBuf>), String> {
    let raw = dir.join("raw");
    let stl = raw.join(format!("{id}.stl"));
    if stl.exists() {
        let bytes = std::fs::read(&stl).map_err(|e| e.to_string())?;
        return read_stl(&bytes).map(|s| (s, None));
    }
    let obj = raw.join(format!("{id}.obj"));
    let mut log = mesher_core::RecordLog::new();
    let model = mesher_io::read_obj(&obj, &mut log, &mut |_| {})
        .ok_or_else(|| format!("{}: unreadable", obj.display()))?;
    let surface = Surface {
        points: model.vertices.clone(),
        triangles: model.triangles.iter().map(|t| t.vertices).collect(),
    };
    Ok((surface, Some(obj)))
}

/// One model of the stratified sample: what each mesher did with it.
struct SampleRow {
    id: u64,
    class: String,
    faces: usize,
    exit: i32,
    path: &'static str,
    kernels: String,
    tets: u64,
    min_dihedral: f64,
    non_orthogonality: f64,
    seconds: f64,
    tetgen: Result<Reference, String>,
    tetgen_seconds: f64,
    histogram: Vec<u64>,
    records: String,
}

/// The pipeline path a run took, from its exit code and kernels.
fn path_of(exit: i32, kernels: &[Kernel]) -> &'static str {
    let conformal = !kernels.is_empty() && kernels.iter().all(|k| *k == Kernel::Delaunay);
    match exit {
        0 if conformal => "direct",
        10 if conformal => "repaired",
        0 | 10 | 20 => "fallback",
        30 => "failed",
        70 => "bug",
        _ => "no report",
    }
}

fn sample_ids() -> Vec<(u64, String, usize)> {
    let name = std::env::var("THINGI10K_SAMPLE").unwrap_or_else(|_| "sample1000.txt".into());
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../thingi10k")
        .join(name);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            (
                f[0].parse().unwrap(),
                f.get(1).copied().unwrap_or("unknown").to_owned(),
                f.get(2).and_then(|n| n.parse().ok()).unwrap_or(0),
            )
        })
        .collect()
}

fn env_number(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

/// Mesh one model of the sample with the CLI, then TetGen at the same size.
fn sample_one(
    cli: &Path,
    worker: &Path,
    dir: &Path,
    (id, class, faces): &(u64, String, usize),
    timeout: u64,
) -> SampleRow {
    let mut row = SampleRow {
        id: *id,
        class: class.clone(),
        faces: *faces,
        exit: -1,
        path: "no report",
        kernels: String::new(),
        tets: 0,
        min_dihedral: f64::NAN,
        non_orthogonality: f64::NAN,
        seconds: 0.0,
        tetgen: Err("not run".into()),
        tetgen_seconds: 0.0,
        histogram: Vec::new(),
        records: String::new(),
    };
    let (surface, native) = match load_surface(dir, *id) {
        Ok(s) => s,
        Err(e) => {
            row.path = "unreadable";
            row.records = e;
            return row;
        }
    };
    row.faces = surface.triangles.len();
    // A model that exists only as OBJ is meshed from that file as it is.
    let obj = native.unwrap_or_else(|| {
        let obj = dir.join("obj").join(format!("{id}.obj"));
        std::fs::write(&obj, to_obj(&surface, &format!("thing_{id}"))).unwrap();
        obj
    });
    let out = dir.join("sample").join(id.to_string());
    let _ = std::fs::remove_dir_all(&out);
    let started = Instant::now();
    let status = Command::new(cli)
        .arg("run")
        .arg(&obj)
        .arg("--out")
        .arg(&out)
        .arg("--worker")
        .arg(worker)
        .args(["--timeout", &timeout.to_string(), "--memory-limit", "3072"])
        .output()
        .unwrap();
    row.seconds = started.elapsed().as_secs_f64();
    row.exit = status.status.code().unwrap_or(-1);
    let report: Option<Report> = std::fs::read(out.join("report.json"))
        .ok()
        .and_then(|b| serde_json::from_slice(&b).ok());
    let mut size = None;
    if let Some(report) = &report {
        let meshes: Vec<_> = report
            .parts
            .iter()
            .filter_map(|p| p.mesh.clone())
            .filter(|m| m.accepted)
            .collect();
        let kernels: Vec<Kernel> = meshes.iter().map(|m| m.kernel).collect();
        row.path = path_of(row.exit, &kernels);
        row.kernels = kernels
            .iter()
            .map(|k| format!("{k:?}").to_lowercase())
            .collect::<Vec<_>>()
            .join("+");
        let mut codes: Vec<&str> = report
            .records
            .iter()
            .filter(|r| r.severity != Severity::Info)
            .map(|r| r.code.as_str())
            .collect();
        codes.sort_unstable();
        codes.dedup();
        row.records = codes.join(" ");
        if !meshes.is_empty() {
            row.tets = meshes.iter().map(|m| m.tets).sum();
            row.min_dihedral = meshes
                .iter()
                .map(|m| m.min_dihedral_degrees)
                .fold(f64::INFINITY, f64::min);
            let weight: f64 = meshes.iter().map(|m| m.tets as f64).sum();
            row.non_orthogonality = meshes
                .iter()
                .map(|m| m.mean_non_orthogonality_degrees * m.tets as f64)
                .sum::<f64>()
                / weight.max(1.0);
            for m in &meshes {
                row.histogram
                    .resize(m.dihedral_histogram.len().max(row.histogram.len()), 0);
                for (a, b) in row.histogram.iter_mut().zip(&m.dihedral_histogram) {
                    *a += b;
                }
            }
            size = meshes.first().map(|m| m.element_size);
        }
    }
    // TetGen at the same size: ours, else 0.03 of the model's diagonal.
    let h = size.unwrap_or_else(|| {
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in &surface.points {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        0.03 * (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt()
    });
    if h.is_finite() && h > 0.0 {
        let started = Instant::now();
        row.tetgen = tetgen_in_child(*id, h);
        row.tetgen_seconds = started.elapsed().as_secs_f64();
    }
    row
}

/// Milestone 10 benchmark: the stratified sample (`sample1000.txt`, or
/// `$THINGI10K_SAMPLE`), meshed in parallel (`$THINGI10K_JOBS`, default 12)
/// by the CLI and by TetGen. Writes `thingi10k_sample.txt` (outcomes by path
/// and class, quality against TetGen) and `thingi10k_sample.csv` (one row a
/// model) to `target/tmp/oracle/`. Download the meshes first with
/// `python tests/thingi10k/fetch.py tests/thingi10k/sample1000.txt`.
#[test]
#[ignore = "needs the downloaded Thingi10k sample and built binaries; takes hours"]
fn stratified_sample_benchmark() {
    let (cli, worker) = binaries();
    let dir = data_dir();
    std::fs::create_dir_all(dir.join("obj")).unwrap();
    let ids = sample_ids();
    let jobs = env_number("THINGI10K_JOBS", 12).max(1) as usize;
    let timeout = env_number("THINGI10K_TIMEOUT", 600);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let rows = std::sync::Mutex::new(Vec::new());
    let started = Instant::now();
    std::thread::scope(|scope| {
        for _ in 0..jobs {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some(model) = ids.get(i) else {
                        return;
                    };
                    let row = sample_one(&cli, &worker, &dir, model, timeout);
                    eprintln!(
                        "[{:>4}/{}] {:>8} {:<17} {:>6} faces: {:<9} {:>6.1}s | TetGen {}",
                        i + 1,
                        ids.len(),
                        row.id,
                        row.class,
                        row.faces,
                        row.path,
                        row.seconds,
                        match &row.tetgen {
                            Ok(q) => format!("{} cells", q.cells),
                            Err(e) => e.clone(),
                        }
                    );
                    rows.lock().unwrap().push(row);
                }
            });
        }
    });
    let mut rows = rows.into_inner().unwrap();
    rows.sort_by_key(|r| r.id);
    let wall = started.elapsed().as_secs_f64();

    // Outcomes by path, per class.
    let paths = [
        "direct",
        "repaired",
        "fallback",
        "failed",
        "bug",
        "no report",
        "unreadable",
    ];
    let mut classes: Vec<String> = rows.iter().map(|r| r.class.clone()).collect();
    classes.sort();
    classes.dedup();
    let mut text = format!(
        "Thingi10k stratified sample: {} models, {jobs} at a time, {wall:.0} s wall.\n\
         Paths: direct = conformal Delaunay, no warnings; repaired = conformal after recorded repairs;\n\
         fallback = isosurface stuffing (recorded); failed = report with errors, no usable mesh.\n\
         TetGen: pzAfq on the same welded surface, -a the regular-tet volume of the same size, 120 s limit.\n\n",
        rows.len()
    );
    write!(text, "{:<18}", "class").unwrap();
    for p in paths {
        write!(text, "{p:>11}").unwrap();
    }
    writeln!(text, "{:>8} {:>13}", "total", "TetGen meshed").unwrap();
    let mut by_class = |name: &str, sel: &dyn Fn(&SampleRow) -> bool| {
        let chosen: Vec<&SampleRow> = rows.iter().filter(|r| sel(r)).collect();
        write!(text, "{name:<18}").unwrap();
        for p in paths {
            write!(
                text,
                "{:>11}",
                chosen.iter().filter(|r| r.path == p).count()
            )
            .unwrap();
        }
        writeln!(
            text,
            "{:>8} {:>13}",
            chosen.len(),
            chosen.iter().filter(|r| r.tetgen.is_ok()).count()
        )
        .unwrap();
    };
    for c in &classes {
        by_class(c, &|r: &SampleRow| &r.class == c);
    }
    by_class("all", &|_| true);

    // Quality on the models both meshed.
    let (mut ours, mut theirs) = (Totals::default(), Totals::default());
    let mut ours_conformal = Totals::default();
    for r in &rows {
        if let (Ok(q), true) = (&r.tetgen, r.tets > 0) {
            ours.add(&r.histogram, r.min_dihedral, r.non_orthogonality, r.seconds);
            theirs.add(&q.histogram, q.min, q.non_orthogonality, r.tetgen_seconds);
            if r.path == "direct" || r.path == "repaired" {
                ours_conformal.add(&r.histogram, r.min_dihedral, r.non_orthogonality, r.seconds);
            }
        }
    }
    writeln!(
        text,
        "\nQuality on the {} models both meshed:",
        ours.min_dihedral.len()
    )
    .unwrap();
    ours.describe("ours  ", &mut text);
    ours_conformal.describe("ours (conformal only)", &mut text);
    theirs.describe("TetGen", &mut text);
    let seconds: Vec<f64> = rows.iter().map(|r| r.seconds).collect();
    writeln!(
        text,
        "\nOur time per model: median {:.1} s, 90th percentile {:.1} s, max {:.1} s",
        percentile(&seconds, 0.5),
        percentile(&seconds, 0.9),
        percentile(&seconds, 1.0)
    )
    .unwrap();

    let mut csv = String::from(
        "id,class,faces,exit,path,kernels,tets,min_dihedral,mean_non_orthogonality,seconds,tetgen,tetgen_cells,tetgen_min_dihedral,tetgen_non_orthogonality,tetgen_seconds,records\n",
    );
    for r in &rows {
        let (tg, cells, min, non) = match &r.tetgen {
            Ok(q) => (
                "meshed".to_owned(),
                q.cells.to_string(),
                format!("{:.3}", q.min),
                format!("{:.3}", q.non_orthogonality),
            ),
            Err(e) => (
                e.replace(',', ";"),
                String::new(),
                String::new(),
                String::new(),
            ),
        };
        writeln!(
            csv,
            "{},{},{},{},{},{},{},{:.3},{:.3},{:.2},{tg},{cells},{min},{non},{:.2},{}",
            r.id,
            r.class,
            r.faces,
            r.exit,
            r.path,
            r.kernels,
            r.tets,
            r.min_dihedral,
            r.non_orthogonality,
            r.seconds,
            r.tetgen_seconds,
            r.records
        )
        .unwrap();
    }
    eprintln!("\n{text}");
    let tmp = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("oracle");
    std::fs::create_dir_all(&tmp).unwrap();
    std::fs::write(tmp.join("thingi10k_sample.txt"), &text).unwrap();
    std::fs::write(tmp.join("thingi10k_sample.csv"), csv).unwrap();
    let bugs: Vec<u64> = rows
        .iter()
        .filter(|r| r.path == "bug" || r.path == "no report")
        .map(|r| r.id)
        .collect();
    assert!(
        bugs.is_empty(),
        "internal bugs or missing reports: {bugs:?}"
    );
}
