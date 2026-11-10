# ExodusMesher

A pure-Rust automatic mesher that turns OBJ files, as they occur in the wild, into tetrahedral volume meshes for a finite-volume solver. Every run ends in a valid mesh or a report that says exactly why not: structured defect codes, locations back to OBJ lines, every repair and every meshing attempt.

- **Inputs:** OBJ with triangles, quads and n-gons, several objects and groups, and the usual defects (duplicate vertices, degenerate faces, holes, non-manifold edges, flipped or inconsistent orientation, self-intersections, touching or nested shells, unknown units).
- **Outputs:** Exodus II (`mesh.exo`, one `TETRA4` block per part, side sets per part boundary and OBJ group), Gmsh `.msh` v4, VTK `.vtu`, the repaired surface (`repaired.obj`), a defect overlay for ParaView (`defects.vtp`) and `report.json`.
- **Meshing:** a constrained Delaunay kernel that conforms exactly to the surface, with a size field graded for curvature and thin walls. When it cannot mesh a part, the part's surface is remeshed and the kernel retried; isosurface stuffing is the last, recorded fallback. Both are written here; TetGen is used only as a test oracle.
- **Front ends:** a CLI and a desktop GUI (egui with a wgpu viewport). Both run the meshing in a separate `mesher-worker` process under a timeout and memory cap, so a crash or a runaway refinement becomes a report, not a hang.

The design, decisions and benchmark results are in [PLAN.md](PLAN.md); the phases after it are in [ROADMAP.md](ROADMAP.md).

## Build

Rust 1.85 or newer:

```bash
cargo build --release -p mesher-cli -p mesher-worker -p mesher-gui
```

The binaries land in `target/release/`; `mesher` and `mesher-gui` look for `mesher-worker` next to themselves (or in `$MESHER_WORKER`). On Linux the GUI needs the GTK 3, xkbcommon and Wayland development packages.

## Command line

```bash
mesher run model.obj
```

meshes every part of `model.obj` into `model_out/` and prints a summary:

```text
mesher: model.obj
  parts (3 selected of 3):
    *1 base: 12 triangles, 1 shell; repaired to 12 faces, gate passed
       meshed: 7923 tets, dihedral 10.4°–163.0°, deviation 0.00% of diagonal (Delaunay)
  ...
  exit 0: clean
```

| Option | Meaning |
| --- | --- |
| `-o, --out DIR` | Output directory (default `<input stem>_out`) |
| `--parts a,b,3` | Parts to mesh, by name or index (default: all) |
| `--element-size F` | Target element size, relative to each part's bounding-box diagonal (default 0.03) |
| `--max-deviation F` | Largest acceptable fallback deviation, relative to the diagonal (default 0.05) |
| `--strict` | Conformal meshes only; never use the stuffing fallback |
| `--uniform-size` | Mesh at the element size everywhere (by default the size is graded down to a quarter of it where the surface curves or the part is thin) |
| `--no-remesh` | When a part's conformal mesh fails, skip remeshing its surface and retrying before the stuffing fallback |
| `--kernel exact` | Try the exact constrained Delaunay kernel first: the surface tetrahedralized exactly and rounded to f64, without refinement (so few, often thin, tets); the conformal kernel and stuffing follow if it fails. Default `conformal` |
| `--no-repair` | Diagnose and report only (a linter) |
| `--timeout S`, `--memory-limit MiB` | Limits for the worker (defaults 3600 s, 8192 MiB) |

`mesher schema` prints the JSON schema of `report.json` ([schema/report.schema.json](schema/report.schema.json)).

**Parts.** Each `o` object is a part; in a file without objects, each `g` group; in a file with neither, each connected component (with enclosed cavities kept in the component around them).

**Exit codes**, for scripts and CI:

| Code | Meaning |
| --- | --- |
| 0 | Clean: conformal mesh, no warnings |
| 10 | Succeeded after recorded repairs |
| 20 | Succeeded via the stuffing fallback (the deviation is reported) |
| 30 | Failed; the report says why |
| 70 | Internal bug (a worker panic) |

## Project files

The GUI saves a project with every export, and the CLI replays it:

```toml
input = "model.obj"
output = "model_out"

[settings]
parts = ["body", "lid"]
element_size = 0.03
```

```bash
mesher run model_out/project.toml
```

Paths are relative to the project file. Options on the command line override the project's settings.

## GUI

```bash
mesher-gui model.obj
```

Open an OBJ or a project (or drop one on the window); pick parts in the list or by clicking them in the viewport; inspect diagnosis records, which highlight their faces and point to their OBJ lines; set the element size and policy; Run and export. The viewport orbits (drag), pans (right-drag), zooms (scroll) and fits (double-click); wireframe, defect overlay and a clipping plane are toggles in the toolbar.

## The report

`report.json` holds the run's config and tool versions, every record, each meshing attempt and the final mesh quality per part. A record looks like this:

```json
{
  "code": "TOP_NONMANIFOLD_EDGE",
  "severity": "error",
  "stage": "diagnose",
  "location": { "vertices": [1042, 1043], "faces": [2201, 2209, 3310], "obj_lines": [1088, 1096, 2197] },
  "measured": { "incident_faces": 3 },
  "threshold": { "max_incident_faces": 2 },
  "action": "split_edge",
  "outcome": "repaired",
  "hint": "Three faces share this edge; likely two touching parts exported as one object."
}
```

Severities are `info`, `warning` (repaired; output still valid), `error` (blocks a stage; a fallback may still succeed) and `fatal`. Mesh quality is reported with finite-volume metrics: face non-orthogonality, skewness, neighbour volume ratio and the dihedral angle range and histogram.

**[docs/report.md](docs/report.md)** is the full reference: every field, every defect code with its default severity and action, what the meshing records measure, the quality metrics, the exit codes, and the other output files.

## Development

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

- `tests/synthetic/`: the defect-injection suite. Every injected defect must be found at its exact location and none reported on clean shapes; `cargo test -p mesher-synthetic -- --nocapture` prints precision and recall per code. Each run is compared with `tests/synthetic/golden.txt`; regenerate it with `UPDATE_GOLDEN=1 cargo test -p mesher-synthetic` after an intended change.
- `tests/oracle/`: comparisons against TetGen (needs a C++ compiler; MSVC on Windows). Quality samples go to `target/tmp/oracle/`.
- After changing report types, regenerate the schema: `UPDATE_SCHEMA=1 cargo test -p mesher-core --test schema`.
- Thingi10k benchmarks (local only; the meshes carry their authors' licences and are not in the repository):

  ```bash
  python tests/thingi10k/fetch.py tests/thingi10k/sample1000.txt
  cargo build --release -p mesher-cli -p mesher-worker
  cargo test --release -p mesher-oracle --test thingi10k stratified_sample_benchmark -- --ignored --nocapture
  ```

## Benchmark

The stratified 1,000-model Thingi10k sample was run on Nov 8, 2026, with default settings, against TetGen (`pzAfq`) at the same element size. Full details are in [PLAN.md](PLAN.md) (Milestone 10).

| Class | Models | Conformal (exit 0 or 10) | Stuffing fallback (20) | Failed (30) | TetGen meshed |
| --- | ---: | ---: | ---: | ---: | ---: |
| Clean | 350 | 304 | 40 | 6 | 323 |
| Multi-component | 150 | 114 | 32 | 4 | 126 |
| Non-manifold | 120 | 59 | 57 | 4 | 72 |
| Open | 80 | 65 | 11 | 4 | 19 |
| Self-intersecting | 300 | 26 | 235 | 39 | 3 |
| All | 1,000 | 568 | 375 | 57 | 543 |

On the 534 models both meshed, the median per-model minimum dihedral is 2.98° for ExodusMesher (0.90° for its conformal meshes alone) and 0.72° for TetGen. The median mean non-orthogonality is 23.2° against 27.2°. ExodusMesher takes about 15 times TetGen's total time (median 38 s a model).

## Licensing

MIT (see [LICENSE](LICENSE)). The shipped crates depend only on permissively licensed Rust crates. TetGen (AGPL) is a dev-dependency of the oracle tests and is never linked into the binaries. The Thingi10k models used for benchmarking are not redistributed.
