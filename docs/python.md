# The Python API (`mesher`)

The distribution is `exodusmesher`; the import is `mesher`. It runs the same worker as the command line and the GUI, so a call with the same options as `mesher run` builds the same job and writes the same files, `report.json` among them, byte for byte (CI checks this on every synthetic input on Linux and Windows).

```python
import mesher

report = mesher.run("part.obj", element_size=0.02, formats=["exodus", "meshb"])
print(report.exit_code)                     # 0, 10, 20, 30 or 70, as `mesher run`
for part in report.parts:
    if part.mesh is not None:
        print(part.name, part.mesh.tets, part.mesh.min_dihedral_degrees)
mesh = report.mesh()                        # points, tets, boundary triangles
print(mesh.points.shape, mesh.parts, mesh.patches)
```

## Install

There is no release on PyPI. Build the wheel from the repository (Rust 1.85 or newer, and maturin 1.16 on `PATH`):

```bash
python -I crates/mesher-py/build_wheel.py        # writes target/wheels/exodusmesher-*.whl
python -m pip install target/wheels/exodusmesher-*.whl
```

`build_wheel.py` builds the release `mesher-worker` and copies it into the package (`mesher/bin/`), then builds the wheel: one per platform, for Python 3.9 and later (abi3). CI uploads the Linux and Windows wheels as the artifacts `wheel-ubuntu-latest` and `wheel-windows-latest`. numpy is optional: with it, meshes are numpy arrays, without it tuples.

The worker runs as a separate process, under a timeout and a memory cap, so a crash or a runaway job ends in a report and never takes the interpreter down. `$MESHER_WORKER` (or a function's `worker=` argument) runs another worker binary instead of the packaged one.

## Running jobs

### `mesher.run(input, out=None, *, ...) -> Report`

Meshes `input` (an OBJ or `.meshb` file, or a project file) as `mesher run` does with the same options, and returns the typed report.

| Argument | Command line | Meaning |
| --- | --- | --- |
| `out` | `--out` | Output directory (default: the project's, else `<input stem>_out` in the current directory) |
| `parts` | `--parts` | Names or 1-based indices of the parts to mesh (default: all) |
| `element_size` | `--element-size` | Relative to each part's bounding-box diagonal (default 0.03) |
| `max_deviation` | `--max-deviation` | Largest acceptable fallback deviation, relative to the diagonal (default 0.05) |
| `kernel` | `--kernel` | `"conformal"` (default) or `"exact"`, or a `KernelChoice` |
| `formats` | `--formats` | Any of `"exodus"`, `"gmsh"`, `"vtu"`, `"openfoam"`, `"fluent"`, `"meshb"`, or `OutputFormat`s (default the first three) |
| `cells` | `--cells` | `"tet"` (default) or `"dual"` (experimental: the polyhedral dual, not written to Gmsh or `.meshb`), or a `CellChoice` |
| `strict` | `--strict` | Conformal meshes only; never the stuffing fallback |
| `no_repair` | `--no-repair` | Diagnose and report only |
| `no_remesh` | `--no-remesh` | Never remesh a part's surface |
| `uniform_size` | `--uniform-size` | No grading for curvature or thin walls |
| `wrap_alpha` | `--wrap-alpha` | The wrap's α, relative to the element size (default 1) |
| `feature_angle` | `--feature-angle` | Degrees |
| `timeout` | `--timeout` | Seconds (default 3600) |
| `memory_limit` | `--memory-limit` | MiB (default 8192; `None`: no cap) |
| `worker` | | The worker binary (default `worker_path()`) |
| `progress` | | `progress(stage, fraction)`, called as the worker reports progress |

Unset options keep the project's value or the default; given with a project file, they override its settings, as on the command line. The options are checked before anything runs: one out of range raises `ValueError`. Everything after that, including a file that cannot be read or parsed, is a record in the report, never an exception.

### `mesher.run_project(path, *, timeout, memory_limit, worker, progress) -> Report`

Runs a project file (`.toml`, as the GUI saves with every export): its input, parts, output directory and settings. A project that cannot be used gives a report with `PARSE_PROJECT_INVALID`; a path that is not a `.toml` file raises `ValueError`. To override a project's settings, pass it to `run` with the options.

### Cancelling

Ctrl-C, or an exception raised by the `progress` callback, cancels the worker. The exception is raised once the worker has written its report (`MESH_WORKER_FAILED` with reason `cancelled`), so the output directory always holds a report.

```python
def progress(stage, fraction):
    print(f"{stage}: {fraction:.0%}")

report = mesher.run("part.obj", progress=progress, timeout=600)
```

### Lower levels

| Function | Returns | |
| --- | --- | --- |
| `build_job(request)` | `dict` or `None` | The job `mesher run` would build: `request` has `input` and the options above by their names (`out`, `parts`, ...). `None` for a project that cannot be used; `ValueError` for an unknown option or one out of range. |
| `run_request(request, *, timeout, memory_limit, worker, progress)` | `RawRun` | Runs a request; what `run` calls. |
| `run_job(job, *, timeout, memory_limit, worker, progress)` | `RawRun` | Runs a job in the worker protocol's form (`input`, `output_dir`, `config`, the last a `JobConfig` as JSON). |
| `worker_path()` | `str` | `$MESHER_WORKER`, else the packaged worker. |

`RawRun` is a named tuple of `report_json` (the worker's text) and `report_path` (where `report.json` was written, or `None`); `RawRun.report()` types it.

## The report

`Report` is `report.json` as frozen dataclasses and `str` enums, generated from the [report schema](../schema/report.schema.json); every field is described in the [report reference](report.md). Lists are tuples, an absent optional field is `None`, and enum members are named by their values in upper case (`Severity.FATAL`, `DefectCode.TOP_NONMANIFOLD_EDGE`, `Kernel.STUFFING`). Being `str` enums, they compare equal to the JSON's strings.

| Member | |
| --- | --- |
| `Report.from_json(text, path=None)`, `mesher.read_report(path)` | A report from the worker's JSON, or from a `report.json` file |
| `to_json()` | The worker's text, exactly as `report.json` holds it |
| `path` | Where `report.json` was written (`None` if not read from a run) |
| `exit_code` | The exit code `mesher run` gives: 0 clean, 10 succeeded after repairs, 20 succeeded by the stuffing fallback, 30 failed, 70 an internal bug |
| `records_by_code` | The records grouped by `DefectCode`, in report order |
| `part(name)` | The `PartSummary` named `name` (`KeyError` if none) |
| `mesh(*, numpy=None)` | The run's mesh (below); needs `"meshb"` among the run's formats |
| `next_job(out=None, *, severity=Severity.ERROR)` | The job to run next (below) |
| `input`, `config`, `stages_run`, `parts`, `records`, `attempts`, `summary`, `tool`, `schema_version` | The report's fields |

```python
blocking = {mesher.Severity.ERROR, mesher.Severity.FATAL}
for record in report.records:
    if record.severity in blocking:
        print(record.code, record.part, record.hint)
```

Enum members are strings, so `<` and `>=` compare their text, not their severity: test membership as above, or count with `report.summary.errors` and `report.summary.fatal`.

A part's `mesh` (`MeshSummary`) holds its tet and node counts, the kernel that meshed it, whether the volume gate accepted it, and its quality: `min_dihedral_degrees`, `max_dihedral_degrees`, `max_non_orthogonality_degrees`, `mean_non_orthogonality_degrees`, `max_skewness`, `max_volume_ratio`, `max_size_ratio`, `hausdorff_relative`, and its `worst_cells`. Its `repair` (`RepairSummary`) holds the face and vertex counts before and after repair and whether the surface gate passed.

## Reading meshes

`report.mesh()`, or `mesher.read_mesh(path)` for an output directory or a `.meshb` file, reads the mesh back from the run's `mesh.meshb` (with `mesh_refs.txt`, which names its references). A run writes it only with `"meshb"` among its formats; `report.mesh()` raises `ValueError` for a run that did not, and `FileNotFoundError` when no part meshed.

| `Mesh` field | |
| --- | --- |
| `points` | n × 3 coordinates (float64) |
| `tets` | m × 4 indices into `points`, 0-based (int64) |
| `tet_parts` | Per tetrahedron, its index in `parts` |
| `parts` | The names of the parts that meshed |
| `triangles` | k × 3: each boundary face once, its normal pointing out of its cell |
| `triangle_patches` | Per triangle, its index in `patches` |
| `patches` | The boundary's side sets that hold faces: each OBJ group, and `<part>_boundary` for faces in none (faces between two parts are internal, in no patch) |
| `path` | The `.meshb` file read |
| `part_tets(name)`, `patch_triangles(name)` | One part's tetrahedra, one patch's triangles (`KeyError` for an unknown name) |

With numpy (`numpy=None` uses it when installed, `True` requires it) the arrays are read-only numpy arrays; with `numpy=False`, or without numpy, they are tuples of tuples.

```python
import numpy as np

mesh = mesher.run("part.obj", formats=["exodus", "meshb"]).mesh()
p = mesh.points[mesh.tets]                                  # m × 4 × 3
volumes = np.einsum("ij,ij->i", p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])) / 6
print(volumes.sum(), volumes.min() > 0)
```

## Remedies and the next job

A record that left a part without a mesh, or a conformal failure that a fallback covered, carries `remedies` when a settings change is known to help: each a `Remedy` with a `setting` (a `Setting`, naming a field of `config`), either a new `value` or a `scale` to multiply it by, and a `reason`, best first. The table of remedies per code is in the [report reference](report.md#remedies).

`report.next_job(out=None, *, severity=Severity.ERROR)` returns a `NextJob` (`job`, `record`, `remedy`) for the first remedy that applies and changes the settings, of the worst record at `severity` or above, or `None` when the run needs nothing or nothing known would help. The job writes to `out`, by default the run's own output directory, replacing its files. `severity=Severity.WARNING` also acts on failures a fallback covered, to reach a conformal mesh.

```python
report = mesher.run("part.obj", "out/try0")
for attempt in range(1, 4):
    following = report.next_job(f"out/try{attempt}")
    if following is None:
        break
    print(following.record.code, "->", following.remedy.setting, following.remedy.reason)
    report = mesher.run_job(following.job).report()
print(report.exit_code)
```

## Types

`mesher` exports every report type: `Report`, `PartSummary`, `MeshSummary`, `RepairSummary`, `FeatureSummary`, `Record`, `Location`, `Remedy`, `Attempt`, `Summary`, `CodeCount`, `JobConfig`, `SizeControl`, `InputInfo`, `InputStats`, `ToolInfo`, `WorstCell`, and the enums `DefectCode`, `Severity`, `Stage`, `Action`, `Outcome`, `Kernel`, `KernelChoice`, `OutputFormat`, `PartSource`, `AttemptOutcome`, `QualityMeasure` and `Setting`. The package ships type information (`py.typed`).

The types are generated from the schema by a Rust test, which CI runs: after changing report types, `UPDATE_SCHEMA=1 cargo test -p mesher-core --test schema`, then `UPDATE_REPORT_MODULE=1 cargo test -p mesher-py --test report_module` (it writes `crates/mesher-py/python/mesher/_schema.py`).

## Checks

- `python -I crates/mesher-py/tests/python/smoke.py <dir>`, against an installed wheel: every function here on the cube fixture, with numpy and without.
- `python -I crates/mesher-py/tests/python/compare_cli.py <inputs> target/release/mesher`: every input the synthetic suite writes with `SUITE_INPUT_DIR=<inputs>` meshed by the CLI and by `mesher.run`, every file the same byte for byte.
- `cargo clippy -p mesher-py --features python --all-targets -- -D warnings` lints the module (PyO3 is behind the `python` feature, off by default).
