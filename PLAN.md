# OBJ Automatic Mesher — Project Plan

Nov 7, 2026 · @Nicholas Saulnier

## Goals and scope

Build, entirely in Rust, a mesher that turns any .obj file into a validated surface mesh and tetrahedral volume mesh, with its own volume mesher as an experiment and a desktop GUI for choosing which parts to mesh. Every input ends in a valid mesh or a report explaining exactly why not.

**Input class.** OBJ files as they occur in the wild: triangles, quads and n-gons; multiple objects and groups; duplicate vertices; degenerate faces; holes; non-manifold edges; self-intersections; inconsistent orientation; several or nested shells; unknown units. Normals, texture coordinates and materials are read and ignored.

**Parts.** Each `o` object is a part (or, in a file without `o`, each `g` group); a file with neither is split into connected components, with enclosed cavities kept in the component around them. The GUI and the CLI's `--parts` flag select which parts to mesh.

**Outputs.**

- Exodus II (.exo): one element block per selected part, side sets for part boundaries
- Gmsh .msh v4 and VTK .vtu
- Repaired surface mesh (.obj or .ply)
- `report.json`: every defect, repair, attempt and quality metric
- A defect overlay (.vtp) for ParaView, also drawn live in the GUI
- A short human-readable summary on stdout

**Non-goals for v1.** Boundary-layer prisms, hex meshing, CAD (STEP/IGES) input, TETRA10 elements, conformal interfaces between touching parts, distributed parallel meshing, and physics-driven sizing. Each is a candidate for later phases.

**Target solver.** A finite-volume solver. Cell quality is judged by finite-volume metrics (face non-orthogonality, skewness, volume ratio between neighbouring cells), the downstream check is a finite-volume solve, and the solver reads Exodus II directly, so its side sets become the boundary patches and no other solver format is needed. Tetrahedra are the v1 cell type; a polyhedral dual mesh, which finite-volume codes often prefer for fewer cells and better gradients, is a later phase.

**Success criteria (proposed).**

- Every input produces either a valid mesh or a report with specific defect codes and locations: zero crashes, zero silent repairs.
- On the synthetic defect suite, every injected defect is detected and none are reported on clean shapes.
- The in-house mesher is measured against TetGen on the same inputs: validity, element quality and run time.
- Exodus files open in ParaView and read cleanly with SEACAS tools.
- Same input and config give byte-identical outputs.

## Key decisions

Everything ships as pure Rust in one Cargo workspace; the only C++ in the project is TetGen, used as a test oracle and never linked into the shipped binaries.

| Concern | Choice | Why |
| --- | --- | --- |
| Structure | Cargo workspace: core libraries, a worker binary, CLI, GUI | One engine, two front ends |
| Exact predicates | `robust` crate (orient3d, insphere) | Adaptive-precision Shewchuk predicates; every topological decision rests on them |
| Linear algebra | `[f64; 3]` helpers now; `nalgebra` once matrices are needed | Points need nothing more; avoids a dependency until smoothing and sizing |
| Spatial queries | Own BVH | Intersection, nearest-point and ray queries; exact boxes, so overlap candidates are never missed |
| Inside/outside | Generalized winding number, own implementation, with an exact boundary-cap fast path | Stays meaningful on surfaces with holes and self-intersections |
| Parallelism | `rayon` | Diagnosis and quality checks run per face or per tet |
| Volume meshing | Own kernels: isosurface stuffing and constrained Delaunay | The experiment; see Volume mesher |
| Test oracle | TetGen via [`tritet`](https://docs.rs/crate/tritet/latest), dev-dependency only | Benchmark against a mature mesher without shipping AGPL code |
| Exodus II | Own writer on the pure-Rust [`netcdf3`](https://docs.rs/netcdf3) crate | No libnetcdf or libhdf5 required; see Exodus II output |
| Other formats | Own .msh v4, .vtu and .vtp writers | Small, well-documented formats |
| Report | `serde` plus `schemars` | Versioned JSON with a generated schema |
| Errors | `thiserror` in libraries, enumerated defect codes | Every failure path maps to a code |
| CLI | `clap` | Standard |
| GUI | `egui`/`eframe` with a `wgpu` viewport | Immediate-mode UI, native 3D, single binary |
| Downstream check | Own cell-centred finite-volume Poisson solver with a sparse solver from `faer` or `sprs` | Proves output meshes are solver-ready |

**Licensing.** Shipped binaries carry only permissive Rust dependencies. TetGen (AGPL) appears only in the test harness; confirm crate licences before reusing anything at work.

**Boundary with work.** Develop on your own machine with public data only, so nothing touches ITAR-controlled geometry or your employer's AFLR pipeline.

## Pipeline architecture

```text
Parse OBJ -> Diagnose -> Select parts -> Repair -> [Surface valid?]
   yes -> Sizing + remesh -> Delaunay mesher -> [Volume valid?]
            yes -> Write mesh + report (Exodus, .msh, .vtu, report.json)
            no  -> Stuffing fallback
   no  -> Stuffing fallback -> [Deviation OK?]
            yes -> Write mesh + report
            no  -> Fail with report (defects, attempts, fixes)

Runs per selected part; meshing runs in mesher-worker; every stage writes to report.json.
```

The pipeline runs once per selected part, top to bottom through two validation gates; a failure at either gate routes to the isosurface-stuffing fallback instead of stopping.

1. **Parse.** Own tolerant parser: triangulates quads and n-gons, handles negative indices and `o`/`g` groups, and keeps a map from every vertex and face to its OBJ line.
2. **Diagnose.** A read-only pass that splits the model into parts and inventories every defect per part. `--no-repair` simply stops here.
3. **Select parts.** Picked in the GUI or named with `--parts`; the choice is saved in a project file.
4. **Repair.** Ordered, conservative operations (see Fallback and repair strategy), each logged with before and after.
5. **Surface gate.** Watertight, edge-manifold, consistently and outwardly oriented, no self-intersections.
6. **Sizing and remesh.** Sizing field from curvature and local feature size, isotropic remeshing, then the gate checks again.
7. **Delaunay mesher.** The in-house constrained Delaunay kernel, run inside the worker process.
8. **Volume gate.** Positive orientation, quality bounds, volume match, boundary conformity.
9. **Emit.** Exodus II, .msh and .vtu, plus `report.json`, the defect overlay and the stdout summary.

**Process model.** The CLI and GUI are thin front ends. Both launch `mesher-worker`, which runs one job and streams JSON lines (progress, report records, result) on stdout. Panics, runaway refinement and hangs stay inside the worker, and cancelling is just killing it.

**Workspace layout.**

```text
crates/
  mesher-geom/    predicates, half-edge mesh, BVH, winding number
  mesher-core/    parts, diagnose, repair, gates, sizing, report types
  mesher-tet/     volume kernels: stuffing, Delaunay, recovery, refinement
  mesher-io/      OBJ reader; Exodus, .msh, .vtu, .vtp writers
  mesher-worker/  binary: one job, JSON-lines protocol on stdout
  mesher-cli/     clap front end, exit codes, project files
  mesher-gui/     egui/eframe app with wgpu viewport
  mesher-fvm/     finite-volume Poisson solver for downstream checks
tests/
  synthetic/      defect injection suite
  oracle/         comparisons against TetGen via tritet
  regression/     Thingi10k subset and golden reports
```

## Volume mesher

Two pure-Rust kernels: an isosurface-stuffing mesher that almost always succeeds, and a constrained Delaunay mesher that conforms exactly. Build the robust one first, so the pipeline produces meshes end to end while the hard kernel matures.

**Kernel A: isosurface stuffing (fallback, built first).** Following Labelle and Shewchuk (2007), fill space with an octree-graded body-centred cubic lattice, cut it against an inside/outside function, then snap and warp lattice vertices so every tetrahedron stays within the method's proven dihedral-angle bounds. Inside/outside comes from the generalized winding number (Jacobson, Kavan and Sorkine-Hornung, 2013), which stays meaningful on surfaces with holes and self-intersections. The boundary is approximated rather than matched, so every run reports its Hausdorff deviation.

**Kernel B: constrained Delaunay with refinement (primary).**

1. **Delaunay tetrahedralization.** Incremental Bowyer–Watson over the surface vertices, with exact orient3d and insphere predicates, randomized insertion order with spatial locality (BRIO), and point location by walking.
2. **Boundary recovery.** Recover missing surface edges, then facets, inserting Steiner points where flips cannot, in the spirit of Si's work behind TetGen. A Steiner-point budget bounds it; an unrecoverable facet becomes `MESH_RECOVERY_FAILED` with its location, and the job falls back to Kernel A.
3. **Region marking.** Remove exterior tets and cavities using the winding number.
4. **Quality refinement.** Delaunay refinement on radius-edge ratio and the size field, then sliver removal with flips and vertex smoothing, with refinement capped by an element budget.

**Data structure.** A flat array of tetrahedra with four vertex and four neighbour indices, a ghost vertex at infinity to simplify hull handling, and free-list reuse of deleted slots. Use `u32` newtype indices rather than reference-counted graphs, which keeps it cache-friendly and borrow-checker friendly.

**Testing against TetGen.** For the same inputs, run TetGen through `tritet` and compare validity, element count, dihedral-angle histograms and wall time. Property tests (`proptest`) check the Delaunay kernel on random point sets, including deliberately degenerate ones: co-spherical, co-planar and duplicated points.

## GUI

An `egui`/`eframe` desktop app that drives the same worker as the CLI and holds no meshing logic, so anything done in the GUI can be replayed from its saved project file with `mesher run project.toml`.

- **Viewport.** A `wgpu` render pass inside eframe: orbit, pan and zoom; shaded plus wireframe; one colour per part; a clipping plane for looking inside volume meshes.
- **Part tree.** A side panel listing each part by its OBJ name (or component number) with triangle count, bounding-box size and defect badges, plus a checkbox to select it.
- **Picking.** Click a part in the viewport to toggle it, via a ray cast against the BVH; box and lasso selection later.
- **Defect overlay.** Report records drawn as coloured faces, edges and markers. Clicking one shows its record, including the OBJ line it came from.
- **Run panel.** Target size per part and the policy flags. Run launches the worker; progress streams into a log and progress bar; Cancel kills the worker.
- **Export.** Writes Exodus, .msh and .vtu plus the report, and saves the project file (input path, selection, settings).

The GUI can start in parallel once diagnosis works, since the viewport, part tree and overlay need nothing from the volume mesher. If you would rather reuse your frontend experience, Tauri with a three.js viewport is a reasonable alternative, still with a Rust core.

## Exodus II output

Write Exodus II directly on the netCDF-3 64-bit-offset format using the pure-Rust [`netcdf3`](https://docs.rs/netcdf3) crate, which reads and writes both classic and 64-bit-offset files. No C libraries are involved, and SEACAS tools and ParaView read the result like any other Exodus file.

| Mesher concept | Exodus entity |
| --- | --- |
| Selected part | Element block, type `TETRA4`, named after the part |
| Mesh vertices | `coordx`, `coordy`, `coordz` |
| A part's outer boundary | Side set named `<part>_boundary` |
| OBJ sub-groups inside a part | Side sets named after the group |
| Boundary vertices | Node sets mirroring each side set (optional) |
| Run provenance | QA records (program, version, date) and info records (report summary) |
| Shared interface between touching parts (later phase) | Side set named `<a>__<b>` |

**Details that bite.**

- Connectivity is 1-based, and side sets use Exodus's own face ordinals for tetrahedra; a wrong ordinal silently puts boundary conditions on the wrong face. Unit-test every ordinal on a single tet.
- Entity names default to a 32-character limit; longer part names are truncated and reported as `EXO_NAME_TRUNCATED`.
- The 64-bit-offset format caps each variable at roughly 4 GiB, enough for meshes into the hundreds of millions of tets. Past that, report `EXO_CAPACITY_EXCEEDED`; CDF-5 or netCDF-4 would be the next step.
- `netcdf3` writes whole variables, not slices, so the full mesh is held in memory while writing.

**Validation.** Round-trip every written file through the `netcdf3` reader in tests, open samples in ParaView, and run SEACAS `io_info` and `exodiff` where installed. The [`mesh-sieve`](https://docs.rs/crate/mesh-sieve/latest/source/docs/exodus.md) crate has an Exodus II writer, but it supports only the HDF5-backed layout and needs libhdf5, so it is useful as a reference rather than a dependency.

## Error reporting design

Every finding is a structured record with a stable code, a severity, a location you can click to, the measured value against its threshold, and what was done about it. Nothing is ever repaired without a record.

**Severity levels.**

- `info`: noteworthy, no action needed (for example, multiple shells)
- `warning`: repaired automatically; output still valid
- `error`: blocks a stage; a fallback may still succeed
- `fatal`: no usable output possible

**Defect taxonomy (v1).**

| Code | Meaning | Default action |
| --- | --- | --- |
| `PARSE_BAD_LINE` | Malformed or non-numeric OBJ line | Skip line, warning |
| `PARSE_BAD_INDEX` | Face references a missing vertex | Drop face, warning |
| `PARSE_READ_FAILED` | Input file cannot be opened or read | Fatal |
| `PARSE_NO_GEOMETRY` | No faces survive parsing | Fatal |
| `PARSE_UNSUPPORTED_ELEMENT` | Lines, points, curves or surfaces (`l`, `p`, `curv`, `surf`, ...) | Skip line, info |
| `PARSE_POLYGON_FALLBACK` | Polygon could not be ear-clipped (degenerate, self-overlapping or over 256 vertices) | Fan-triangulate, warning |
| `PART_EMPTY_SELECTION` | No parts selected, or a named part does not exist | Fatal |
| `PART_OVERLAP` | Two selected parts intersect | Error |
| `PART_SHARED_INTERFACE` | Two selected parts touch along a surface | Info; meshed independently in v1 |
| `GEO_DUPLICATE_VERTEX` | Vertices closer than tolerance | Merge, warning |
| `GEO_DEGENERATE_FACE` | Zero-area or collinear face | Remove, warning |
| `GEO_DUPLICATE_FACE` | Same face listed twice | Remove, warning |
| `TOP_BOUNDARY_EDGE` | Hole: edge with one adjacent face | Fill if small, else error |
| `TOP_NONMANIFOLD_EDGE` | Edge shared by three or more faces | Split, else error |
| `TOP_NONMANIFOLD_VERTEX` | Vertex joining separate fans | Duplicate vertex, warning |
| `TOP_INCONSISTENT_ORIENTATION` | Neighbouring faces disagree on winding | Reorient, warning |
| `TOP_INVERTED_SHELL` | Shell with negative signed volume | Flip, warning |
| `TOP_OPEN_SHEET` | Surface encloses no volume | Fatal for that part |
| `INT_SELF_INTERSECTION` | Faces of one shell intersect | Error, triggers fallback |
| `INT_SHELL_OVERLAP` | Two shells of one part intersect | Error, triggers fallback |
| `QUAL_SLIVER_FACE` | Minimum angle below threshold | Remesh, warning |
| `QUAL_SUSPECT_SCALE` | Bounding box suggests wrong units | Info; tolerances are relative |
| `MESH_RECOVERY_FAILED` | Constrained Delaunay could not recover a facet | Error, triggers fallback |
| `MESH_BUDGET_EXCEEDED` | Steiner-point or element budget hit | Error, triggers fallback |
| `GEO_CAPACITY_EXCEEDED` | Part too large for 32-bit vertex and face indices | Fatal |
| `REPAIR_GUARD_EXCEEDED` | A repair changed the shell count or enclosed volume beyond its limit | Error |
| `MESH_STRICT_UNMET` | `--strict` was given but no conformal mesh could be produced | Error |
| `MESH_WORKER_FAILED` | Worker panicked, timed out or ran out of memory | Error, triggers fallback |
| `VOL_INVERTED_TET` | Tetrahedron with non-positive volume | Error |
| `VOL_LOW_QUALITY` | Non-orthogonality, skewness or volume ratio out of bounds | Warning or error by threshold |
| `VOL_BOUNDARY_DEVIATION` | Hausdorff distance to input over tolerance | Error |
| `EXO_NAME_TRUNCATED` | Part or set name over the Exodus limit | Truncate, warning |
| `EXO_CAPACITY_EXCEEDED` | Mesh too large for the 64-bit-offset format | Error |
| `OUT_WRITE_FAILED` | An output file could not be written | Fatal |

**Record shape.**

```json
{
  "code": "TOP_NONMANIFOLD_EDGE",
  "severity": "error",
  "stage": "diagnose",
  "location": {
    "vertices": [1042, 1043],
    "faces": [2201, 2209, 3310],
    "obj_lines": [1088, 1096, 2197],
    "point": [0.412, -1.07, 3.3]
  },
  "measured": {"incident_faces": 3},
  "threshold": {"max_incident_faces": 2},
  "action": "split_edge",
  "outcome": "repaired",
  "hint": "Three faces share this edge; likely two touching parts exported as one object."
}
```

The report also carries the run's config, tool versions, an ordered `attempts` list (one entry per meshing attempt, with captured stdout and stderr) and final quality metrics.

**Exit codes.** `0` clean, `10` succeeded after repairs, `20` succeeded via fallback, `30` failed with a full report, `70` internal bug. Scripts and CI can branch on these without parsing JSON.

**Isolation.** Pure Rust removes segfaults, but not panics, infinite loops or runaway refinement in an experimental mesher. Every meshing job runs in `mesher-worker` under a timeout and memory cap, so any of those becomes a `MESH_WORKER_FAILED` record rather than a frozen GUI or a dead CLI.

**Provenance.** The parser keeps a map from each vertex and face back to its OBJ line, and every repair updates it, so a defect found late in the pipeline still points at the line that caused it.

## Fallback and repair strategy

Try the most faithful path first and degrade one step at a time, recording each attempt, so the report always says how far from the input the result had to move.

1. **Direct.** Run the constrained Delaunay kernel on the repaired surface. Exact boundary conformity.
2. **Local repair, then retry once.** If facets cannot be recovered, remesh the surface locally around them and retry. Bounded to one retry so it cannot loop.
3. **Stuffing fallback.** Run isosurface stuffing on the winding-number field. The result is valid with bounded dihedral angles but approximates the boundary, so measure and report the Hausdorff distance to the original surface.
4. **Fail with a report.** If stuffing fails or the deviation exceeds `--max-deviation`, stop and emit the full report and defect overlay for that part; other selected parts continue.

**Repair rules.**

- Least invasive first: merge, remove degenerates and reorient before anything that adds geometry.
- Hole filling is automatic only below a size threshold (edge count and area relative to the shell); larger holes may be intentional openings and become errors.
- No repair may change the number of shells or the enclosed volume by more than a set fraction without becoming an error.

**Policy flags.**

- `--parts`: names or indices of the parts to mesh (the GUI writes these into the project file).
- `--strict`: conformal paths only; never use the stuffing fallback.
- `--no-repair`: diagnose and report only, a pure linting mode.
- `--max-deviation`: the largest acceptable fallback deviation, relative to the bounding-box diagonal.

## Validation and benchmarking

Two test sets answer two different questions: a synthetic defect suite proves the error reporting is correct, and Thingi10k measures how robust the mesher is on real-world files.

**Synthetic defect suite.** Start from clean shapes (cube, sphere, torus, a cube with an internal cavity, two touching boxes) and inject one known defect at a time: punch a hole, flip a patch, duplicate vertices, add a sliver, push one shell into another, append a malformed line. Each case asserts the exact codes and locations reported. Track detection precision and recall per code; the target is 100% on both.

**Thingi10k.** About 10,000 real models from Thingiverse, many of them defective, and the benchmark TetWild was evaluated on. Run the full set and publish a table of outcomes by path taken. Keep a fixed subset of around 200 models in CI as a regression gate.

**Metrics per run.**

- Outcome and path: direct, repaired, fallback, failed
- Cell quality: maximum face non-orthogonality, skewness, neighbour volume ratio and minimum dihedral angle, with distribution histograms
- Fidelity: Hausdorff distance to the input surface, relative enclosed-volume error
- Cost: element count, wall time, peak memory

**Volume mesh checks.** All tetrahedra positively oriented; total tet volume matches the surface's enclosed volume; boundary faces of the volume mesh reproduce the surface mesh (conformal paths) or stay within tolerance (fallback).

**Oracle comparison.** For every input, also run TetGen through `tritet` and record both results side by side: validity, element count, dihedral-angle histograms, Hausdorff deviation and wall time. The gap between the two is the experiment's headline result.

**Exodus conformance.** Round-trip every .exo through the `netcdf3` reader and check blocks, side sets and ordinals; spot-check samples in ParaView and with SEACAS tools.

**Downstream usability.** Solve a Poisson problem with an in-house cell-centred finite-volume solver on every output and check convergence against a manufactured solution; non-orthogonal cells show up directly as lost accuracy. This catches meshes that pass geometric checks but are useless to a solver.

**Regression.** Store golden reports for the synthetic suite and CI subset; any change in codes, paths or metrics beyond tolerance fails the build.

## Milestones

Build the reporting spine before any meshing, so every later feature lands with its diagnostics already wired in. Each phase ends at a testable exit criterion.

1. **Foundations.** Workspace, report schema and exit codes, tolerant OBJ parser with line provenance, worker protocol, CLI skeleton.
   - Exit: any file, including garbage, yields a valid `report.json` and never a panic.
2. **Geometry core.** Predicate wrappers, half-edge mesh, BVH, generalized winding number.
   - Exit: property tests agree with brute-force reference implementations on random and degenerate inputs.
3. **Diagnosis and parts.** Part splitting, the full defect taxonomy, defect overlay, synthetic defect suite.
   - Exit: every injected defect detected, zero false positives on clean shapes.
4. **Repair.** Conservative repairs with logging and provenance updates; the surface gate.
   - Exit: repaired synthetic cases pass the gate, and every repair appears in the report.
5. **Writers.** Exodus II, .msh and .vtu, tested first on TetGen output from the oracle.
   - Exit: files round-trip, open in ParaView, and pass SEACAS checks; side-set ordinals verified per face.
6. **Stuffing kernel.** Isosurface stuffing on the winding-number field.
   - Exit: valid tets on every synthetic case with a closed or nearly closed surface; dihedral bounds hold; first end-to-end OBJ to Exodus runs.
7. **Delaunay kernel.** Incremental 3D Delaunay with exact predicates.
   - Exit: on random point sets in general position, element count and connectivity match TetGen's Delaunay output.
8. **Recovery and refinement.** Boundary recovery, region marking, quality refinement, sliver removal.
   - Exit: conformal meshes on the clean part of the Thingi10k subset, with quality histograms reported against TetGen.
9. **GUI.** Viewport, part tree, picking, overlay, run panel with progress and cancel, export, project files. Can start after milestone 3.
   - Exit: select parts of a multi-object OBJ and export an Exodus file with one block per part, reproducible from the CLI.
10. **Sizing, benchmark and hardening.** Sizing field and surface remeshing, full Thingi10k run, profiling, CI regression subset, documentation.
    - Exit: published outcome table by path, against TetGen, with CI green.

## Phases and task breakdown

The milestones group into six phases. Phases run in order, except the GUI phase, which can start once Phase 2 is done. A milestone's tasks are ticked when they land; the milestone is done only when its exit criterion passes under `cargo fmt`, `cargo clippy -D warnings` and `cargo test --workspace`.

| Phase | Milestones | Delivers |
| --- | --- | --- |
| 1. Reporting spine | 1 | Workspace, report schema, parser, worker, CLI: any file in, a valid report out |
| 2. Geometry and diagnosis | 2, 3 | Exact geometry core, part splitting, the full defect taxonomy, synthetic suite |
| 3. Repair and output | 4, 5 | Logged repairs and surface gate; Exodus, .msh and .vtu writers |
| 4. Volume meshing | 6, 7, 8 | Stuffing fallback end to end, then the constrained Delaunay kernel |
| 5. GUI | 9 | egui front end over the worker (can start after Phase 2) |
| 6. Hardening | 10 | Sizing, Thingi10k benchmark, CI regression gate |

Crates are created when their first milestone starts, so early builds do not pull in GUI or solver dependencies.

### Phase 1: Reporting spine

**Milestone 1: Foundations.** Done (Nov 7, 2026): the exit criterion is enforced by `crates/mesher-cli/tests/cli.rs`, which runs the real binaries on hostile and random inputs.

- [x] Workspace: `mesher-geom`, `mesher-core`, `mesher-io`, `mesher-worker`, `mesher-cli`; `.gitattributes` forcing LF; `tests/fixtures/obj/`
- [x] Report types in `mesher-core::report`: `DefectCode` (one variant per taxonomy row, with default severity), `Severity`, `Stage`, `Record`, `Location`, `Attempt`, `Report`, `Summary`
- [x] JSON schema generated with `schemars`, committed as `schema/report.schema.json`; a test fails when it is stale (regenerate with `UPDATE_SCHEMA=1 cargo test -p mesher-core`)
- [x] A test that the `DefectCode` enum and the taxonomy table in this file list the same codes
- [x] Record capping: at most 100 detailed records per code; the summary keeps exact totals per code, so garbage input cannot produce an unbounded report
- [x] Exit-code rule (below) in `mesher-core::exit`
- [x] Polygon triangulation in `mesher-geom`: ear clipping in the plane of the Newell normal with `robust::orient2d`; fan fallback, recorded as `PARSE_POLYGON_FALLBACK`, for degenerate, self-overlapping or very large (over 256 vertices) polygons
- [x] Tolerant OBJ parser in `mesher-io::obj`: byte-level line splitting (LF, CRLF, CR), `\` continuations, comments, invalid UTF-8, negative and forward indices, `v/vt/vn` face tokens, `o` and `g` names; vertex and face OBJ-line provenance
- [x] Worker protocol in `mesher-core::protocol`, worker binary reading the job from stdin and writing `report.json`
- [x] CLI skeleton: `mesher run`, `mesher schema`; worker timeout; panics and timeouts become `MESH_WORKER_FAILED`
- [x] Tests: parser unit tests, property tests on arbitrary bytes and OBJ-like token streams, worker and CLI end-to-end tests on garbage, schema validation of every report, byte-identical reruns
- Deferred: the worker memory cap (Windows job object, `setrlimit` elsewhere) moves to Milestone 6, when meshing can first allocate without bound; project files (`mesher run project.toml`) move to Milestone 9

**Decisions made in Milestone 1.**

- *Exit codes.* Computed from the records: any `error` or `fatal` gives `30`; otherwise any attempt that used the stuffing fallback gives `20`; otherwise any `warning` gives `10`; otherwise `0`. A worker panic gives `70`; a timeout is a `MESH_WORKER_FAILED` error and gives `30`. Command-line usage errors exit with clap's `2`.
- *Index conventions.* Record `vertices` and `faces` are 0-based indices into the parsed arrays; `obj_lines` are 1-based line numbers in the input file.
- *Determinism.* Reports carry no timestamps or absolute paths beyond what the user passed; maps are `BTreeMap`.
- *Worker protocol (version 1).* The job is one JSON document on the worker's stdin. Stdout is JSON lines tagged by `type`: `started` (protocol and worker version), `progress` (stage, fraction), `record` (a report record), `result` (exit code, report path). The worker writes `report.json` itself; if the worker dies, the CLI writes it from the records already streamed plus a `MESH_WORKER_FAILED` record carrying the tail of the worker's stderr. The environment variable `MESHER_TEST_FAULT=panic|hang` makes the worker fail on purpose, for testing that path.
- *Pipeline extent.* The report lists `stages_run`; in Milestone 1 the pipeline stops after `parse`.

### Phase 2: Geometry and diagnosis

**Milestone 2: Geometry core.** Done (Nov 7, 2026): every property test agrees with its brute-force reference, and the oracle tests were also run at 100,000 cases each.

- [x] Predicate wrappers (`orient2d`, `orient3d`, `insphere`) returning a sign enum, plus exact helpers (`collinear`, `plane_axis`, cyclic `project`)
- [x] `TriangleSoup`: validated borrowed input (finite coordinates, indices in range, `3 × triangles` fits `u32`), failing with `GeomError`
- [x] Half-edge mesh built from triangle soup, tolerating non-manifold input (edge-to-faces map alongside): edge kinds `Boundary`, `Manifold`, `Inconsistent`, `NonManifold`, `Degenerate`; vertex fan counts
- [x] BVH over triangles: ray cast, nearest point, triangle-pair overlap candidates (within one soup and between two)
- [x] Exact triangle-triangle intersection test on the predicates, plus `faces_intersect`, which ignores contact along the vertices or edge two mesh faces share
- [x] Generalized winding number, with a hierarchical fast path checked against the brute-force sum
- [x] Property tests against brute-force references, including degenerate inputs (`crates/mesher-geom/tests/`)

**Decisions made in Milestone 2.**

- *Orientation convention.* `orient3d(a, b, c, d)` is positive for a positive-volume tetrahedron (the opposite of Shewchuk's raw sign, which the wrapper flips); `insphere` is positive inside for a positively oriented tetrahedron. Unit tests pin both.
- *Exact intersection.* Two closed convex sets that meet have an extreme point of their intersection on one set's boundary, so triangles meet exactly when an edge of one meets the other. Coplanar work happens in a coordinate projection chosen by an exact nonzero normal component, so no coordinates are ever computed. The reference is an exact rational linear program (enumerating basic feasible solutions of the barycentric system) on integer grids, sharing no code with the implementation.
- *Winding-number fast path.* Not the dipole approximation of Barill et al. (2018) but an exact identity: outside a node's box, a patch has the same winding number as the cone over its boundary loop (Jacobson et al., 2013), so a node costs one solid angle per boundary edge. Caps are stored only where shorter than the node's triangle list; closed regions cost nothing.
- *Vectors.* Points are plain `[f64; 3]` with small helpers in `mesher-geom::vec3`; `nalgebra` waits until a milestone needs matrices (smoothing, sizing).
- *Floating-point queries.* Ray casts and nearest points are measurements, not topology. Their BVH culling is padded so it never drops a triangle a linear scan would accept; tests compare results within `1e-9` relative.

**Milestone 3: Diagnosis and parts.** Done (Nov 7, 2026): the synthetic suite (`tests/synthetic/`: 33 cases, 95 runs across three scales and offsets) shows precision and recall of 1.000 for every code it covers, with zero findings on the clean shapes.

- [x] Map `GeomError` to records: `TooLarge` becomes the new `GEO_CAPACITY_EXCEEDED` (fatal for the part); the other variants cannot come from the parser, so they surface as `MESH_WORKER_FAILED` with exit code 70 (a bug)
- [x] Exclude degenerate faces from intersection tests (they get `GEO_DEGENERATE_FACE`; `faces_intersect` falls back to plain contact for them)
- [x] Part splitting: `o`/`g` names, else connected components (cavities stay with the component that encloses them)
- [x] Detectors for every `GEO_`, `TOP_`, `INT_`, `PART_` and `QUAL_` code (all but `GEO_CAPACITY_EXCEEDED` exercised by the suite)
- [x] Relative tolerances from the bounding-box diagonal; `QUAL_SUSPECT_SCALE`
- [x] Defect overlay writer (.vtp), written by the worker as `defects.vtp`
- [x] Synthetic suite in `tests/synthetic/`: clean shapes plus one-defect injectors, asserting exact codes and locations; precision and recall per code
- [x] Worker runs parse, diagnose and select; the report gains `parts` (index, name, source, triangle, vertex and shell counts, bounds, selection) for the GUI part tree; the CLI prints them

**Decisions made in Milestone 3.**

- *Diagnose on the repaired-to-be surface.* Topology is judged after merging duplicate vertices and setting aside degenerate and duplicate faces, the first steps repair will take. One root cause gives one record: a duplicated vertex is reported once, not as a ring of boundary edges. Every diagnosis record has outcome `deferred`.
- *Locations.* `faces` are triangle indices (a polygon's triangles are consecutive); `vertices` are model indices, with merged vertices named by the smallest index in their cluster. Location lists are capped at 256 entries, with full counts in `measured` (`faces_total` and so on).
- *Records are per root cause, not per element:* one `TOP_BOUNDARY_EDGE` per hole (boundary loop), one `TOP_INCONSISTENT_ORIENTATION` per orientation component (listing the minority faces to flip), one `INT_SHELL_OVERLAP` per pair of shells, but one `INT_SELF_INTERSECTION` per face pair.
- *Thresholds* (`diagnose::Tolerances`, not yet user-configurable): merge tolerance 1e-7 × diagonal; sliver below 1°; open sheet when the shell's volume with holes closed is at most 1e-6 × area^1.5; suspect scale when the diagonal is outside [1e-4, 1e6].
- *Inverted shells use nesting.* A shell's expected volume sign alternates with how many other disjoint shells of its part contain it (by winding number), so cavities are legitimately inward-facing. Shells that intersect each other are excluded from nesting and reported as `INT_SHELL_OVERLAP` instead.
- *Open sheets suppress their holes.* A shell that encloses nothing gets `TOP_OPEN_SHEET` (fatal for the part) and no `TOP_BOUNDARY_EDGE` records, since filling holes cannot help it.
- *Parts.* `o` objects take precedence over `g` groups; faces before the first statement form the part `default`. Without either, parts are connected components through shared merged vertices, and an inward-facing component inside another joins it as a cavity. Parts are ordered by their first triangle; `--parts` takes names or 1-based indices.
- *Part pairs.* Exact contact classification (`mesher_geom::classify_contact`): faces that cross each other, or coplanar faces sharing area with normals aligned, mean overlapping volumes (`PART_OVERLAP`); coplanar faces sharing area with opposed normals form an interface (`PART_SHARED_INTERFACE`); contact along edges or points only is ignored. With no contact at all, a winding-number test catches one part inside another.
- *Known limitation.* Two faces that cross each other's planes but meet at a single point are classified as crossing, so parts touching edge-to-edge at a point can be reported as overlapping. Revisit if real inputs hit it.
- *Exit code during Milestone 3.* Diagnosis warnings gave exit 10 ("succeeded after repairs") before repair existed; since Milestone 4 that is accurate, except in `--no-repair` mode, where 10 means "repairable defects found".
- *Deferred to Milestone 10:* a per-file time budget for the intersection tests, which are quadratic in the worst case (every face box overlapping every other).

### Phase 3: Repair and output

**Milestone 4: Repair.** Done (Nov 7, 2026): in the synthetic suite every repairable finding has a matching applied repair record, every repaired part passes the surface gate except the three cases built to fail it (a hole too large to fill, a self-intersection, overlapping shells), no guard trips, and each part's face and vertex counts before and after are fully explained by its repair records. A mutation check (disabling hole filling, then edge pairing) makes the suite fail exactly on the affected cases.

- [x] Ordered repairs: merge, remove degenerates and duplicates, reorient, flip inverted shells, split non-manifold edges and vertices, fill small holes
- [x] Provenance maps updated by every repair; before-and-after logged per record
- [x] Shell-count and enclosed-volume guards (`REPAIR_GUARD_EXCEEDED`)
- [x] Surface gate; `--no-repair` stops after diagnosis
- [x] T-junctions: a collinear (collapsed) face is removed by splitting the face across its long edge at the middle vertex, in the shared clean-up, so diagnosis and repair agree
- [x] Worker writes `repaired.obj` (one object per repaired part); the report's part summaries gain face and vertex counts and the gate verdict; the CLI shows both
- [x] Property test: diagnosis and repair are total and deterministic on random degenerate-heavy soups (20,000 cases run once)

**Decisions made in Milestone 4.**

- *One analysis, three stages.* `diagnose::analysis::analyze_surface` computes every topological and geometric fact once; diagnosis, repair and the gate each turn it into their own records, so the gate judges exactly what diagnosis judges.
- *Diagnosis predicts the clean-up.* Merging, removing degenerate and duplicate faces, T-junction splits and reorientation happen in diagnosis's clean-up (they are needed to diagnose root causes); repair records them as applied and continues from that surface.
- *Repair record severity.* A fix applied is `info` (its finding already carries the warning, so warnings are not counted twice); a fix that could not be made is `error`. Records carry `faces_added`, `faces_removed` and `vertices_added`, so counts before and after can be audited.
- *Non-manifold splitting is exact.* Faces around an edge are sorted by angle with `orient3d` and paired as solids (a face traversing the edge `v→u` followed counter-clockwise by one traversing `u→v`); fans at each vertex then get their own vertex copies. Copies keep their origin, and the gate's intersection test uses origins, so solids that met at a split edge or vertex are not flagged for touching there.
- *Hole filling.* Only simple loops (each vertex with one boundary edge in and one out) of at most 256 edges whose fill is at most 20% of the closed shell's area; ear clipping in the loop's plane, else a fan from a new centroid vertex. Larger holes are left open as possible intentional openings (`TOP_BOUNDARY_EDGE`, outcome `failed`, severity `error`).
- *Guards.* Vertex-connected pieces (by origin) must not change, and the total enclosed volume (sum of absolute shell volumes, holes closed) must change by at most 1%.
- *Gate.* One `surface_gate` error per failed check per part: holes, non-manifold edges or vertices, inconsistent orientation, inverted shells, open sheets, self-intersections, shell overlaps. Slivers are not gated (remeshing is Milestone 10).
- *Not repaired here:* self-intersections and shell overlaps (the stuffing fallback, Milestone 6), slivers (remeshing, Milestone 10), and bare T-junction slits with no collapsed face to remove (they show as zero-area holes that filling refuses).

**Milestone 5: Writers.** Code, automated checks and the ParaView spot check done (Nov 7, 2026); the SEACAS check is still open because SEACAS is not installed on the development machine (it is easiest to run under Linux or WSL). Files to check are written by `cargo test -p mesher-oracle` to `target/tmp/oracle/` (`l_prism.exo`, `.msh`, `.vtu`).

- [x] Volume mesh type (`mesher-core::mesh::VolumeMesh`): points, one block of positive tets per part, side sets as (global element, Exodus side), boundary side sets, validation
- [x] Exodus II on `netcdf3`: blocks, side sets with tet face ordinals, node sets, QA and info records, name truncation, capacity check
- [x] Gmsh .msh v4 and VTK .vtu writers
- [x] TetGen oracle harness in `tests/oracle/` (dev-dependency `tritet`, built with MSVC); a guard test fails if any crate under `crates/` mentions `tritet`
- [x] Round-trip tests with independent readers (`netcdf3` for Exodus, `mshio` for MSH, `vtkio` for VTU); per-ordinal unit tests on a single tet; byte-identical reruns
- [x] Side-set ordinals verified per face: on the TetGen L-prism, every side read back from the file lies exactly on one facet plane, points out of the solid (exact `orient3d`), and the sides cover each facet's area
- [x] ParaView spot check of `l_prism.exo` (Nov 7, 2026, by eye): the element block, side set and node set load; the L-prism is filled with tets under a clip; Mesh Quality shows every tet volume positive; normal glyphs on the side set point out of the solid on every face, including the notch
- [ ] SEACAS (`io_info`, `exodiff`) check of the oracle files; deferred by decision (Nov 7, 2026) so Milestone 6 could start with ParaView verification in hand. Do it before the Thingi10k benchmark in Milestone 10, under Linux or WSL

**Decisions made in Milestone 5.**

- *Exodus layout.* 32-bit integers (`int64_status` 0), doubles, `file_size` 1 (64-bit offset), `maximum_name_length` 32 with `len_name` 33, an unlimited `time_step` with no steps. Element numbering is global, block by block. Each side set gets a mirroring node set. Empty blocks and sets are written as null entities (status 0, no data), as the Exodus library does; a mesh with no elements at all is refused (`ExodusError::Empty`) rather than written.
- *Side ordinals.* `mesher_core::mesh::TET_SIDES` is the single source: (1,2,4), (2,3,4), (1,4,3), (1,3,2) in Exodus's 1-based numbering, each ordered outward for a positive tet. The .msh writer reuses it for boundary triangles.
- *Determinism.* QA date and time are left empty by default so reruns are byte-identical; callers may fill them.
- *Capacity.* Checked before writing: any variable over 2^32 − 4 bytes, or any count over `i32::MAX`, returns `ExodusError::Capacity`, which the worker will record as `EXO_CAPACITY_EXCEEDED`.
- *Formats.* .msh and .vtu are ASCII for now (simple and diffable); binary variants can follow if file size matters. Side sets appear in .msh as physical surfaces and are not represented in .vtu.
- *Wiring.* Nothing produces tetrahedra yet, so the worker does not call the writers; Milestone 6 connects them (blocks per selected part, `<part>_boundary` side sets, `EXO_NAME_TRUNCATED`, `EXO_CAPACITY_EXCEEDED` and `OUT_WRITE_FAILED` records).
- *Dev-only dependencies.* `vtkio` 0.6.3 and `mshio` pull in `nom` 3 and `quick-xml` 0.22, which trigger future-incompatibility notices; they are test-only and can be replaced if the notices become errors.

### Phase 4: Volume meshing

**Milestone 6: Stuffing kernel.** Done (Nov 7, 2026): every repaired part of every closed or nearly closed synthetic case (98 part meshes across three scales and offsets) is stuffed with positive tets inside the proven dihedral bounds; accepted meshes validate and match the enclosed volume within 15% at the suite's coarse element size; the cube fixture runs end to end from OBJ to Exodus, .msh and .vtu.

- [x] `mesher-tet` crate: uniform BCC lattice, level-function signs (parallel, `rayon`), cut points by bisection, warping with the paper's α (0.24999 long, 0.41189 short), stencils
- [x] Dihedral angles (`mesher_tet::quality`); property test on random spheres, ellipsoids, boxes and tori at random cell sizes: angles within [10.7°, 164.8°], no degenerate tets, every face shared by at most two tets, closed and consistently oriented boundary (1,000 shapes run once)
- [x] Volume stage (`mesher-core::volume`): stuffing on each repaired part's winding-number field, volume gate (`VOL_INVERTED_TET`, `VOL_LOW_QUALITY`, `VOL_BOUNDARY_DEVIATION`), `MESH_BUDGET_EXCEEDED`, attempts recorded; accepted parts assembled into one mesh, a block and a `<part>_boundary` side set each
- [x] Worker writes `mesh.exo`, `mesh.msh`, `mesh.vtu` (`EXO_NAME_TRUNCATED`, `EXO_CAPACITY_EXCEEDED`, `OUT_WRITE_FAILED` recorded); the report's part summaries gain mesh statistics; the CLI prints them
- [x] Worker memory cap: the CLI runs the worker under `process_control` (a job object on Windows, `setrlimit` on Linux), `--memory-limit` MiB (default 8192); allocation failure becomes `MESH_WORKER_FAILED` with reason `out_of_memory` (exit 30, not a bug); tested with a deliberate runaway allocation
- [x] `--element-size` (relative to each part's diagonal, default 0.03); `--strict` records `MESH_STRICT_UNMET` until the conformal kernel exists
- [ ] Octree grading (adaptive lattice): deferred to Milestone 10, where the sizing field gives it something to grade to

**Decisions made in Milestone 6.**

- *All-zero tets are dropped.* A lattice tet whose four vertices were all warped onto the surface is a sliver lying in it; keeping it (as an early version did when its centroid was inside) produced dihedral angles down to 2.8°. Dropping it leaves its faces, which lie on the surface, as boundary. With that, the empirical range on all tests stays inside the paper's [10.7°, 164.8°].
- *Stencil filling by coning.* Each lattice tet's inside part is filled by a cone from its lowest-numbered inside vertex over its boundary triangles; a quadrilateral (only on faces with two inside vertices and one outside) is split from the lower-numbered inside vertex. Both rules depend only on the shared face, so neighbours conform. This replaces the paper's colour-based parity rule; the property tests above are the evidence that the angle bounds still hold.
- *Deviation is measured against the region's boundary.* The Hausdorff distance is two-sided between the mesh boundary and the input surface, sampled at vertices and face centroids, but input samples with winding number ≥ 1 (faces buried inside the region, as with overlapping shells) are skipped.
- *Defaults.* Element size 0.03 × part diagonal; maximum deviation 0.05 × part diagonal (raised from 0.01: stuffing rounds sharp corners, e.g. 1.0% on a lattice-aligned cube and 3.5% on a tetrahedron at the default size).
- *Exit code 20 for every meshed run* until the Delaunay kernel lands: stuffing is by definition the fallback.
- *Known gap:* a feature thinner than the element size is lost (the synthetic spike case is rejected at 15% deviation); the fix is a smaller element size or, later, grading.

**Milestone 7: Delaunay kernel.** Flat tet array with ghost vertex, BRIO insertion, walking point location, Bowyer–Watson; property tests and TetGen comparison.

**Milestone 8: Recovery and refinement.** Edge and facet recovery with a Steiner budget, region marking, Delaunay refinement, sliver removal, finite-volume quality metrics, `mesher-fvm` Poisson check.

### Phase 5: GUI

**Milestone 9: GUI.** `mesher-gui` with wgpu viewport, part tree, picking, overlay, run panel, export, project files shared with the CLI.

### Phase 6: Hardening

**Milestone 10: Sizing, benchmark and hardening.** Sizing field and isotropic remeshing, full Thingi10k run, profiling, CI regression subset with golden reports, documentation.

**After v1.** Conformal interfaces between touching parts (shared side sets in Exodus), a polyhedral dual mesh for the finite-volume solver, TETRA10 elements, boundary-layer prisms, a model that predicts failure from surface features before meshing, and an agent loop that reads the report and chooses repairs.

## Risks and open questions

The biggest technical risk is boundary recovery in the in-house constrained Delaunay kernel: it is where TetGen spent years on degenerate cases, and it may never match TetGen's robustness. The stuffing fallback keeps the pipeline useful regardless.

| Risk | Mitigation |
| --- | --- |
| Boundary recovery fails or loops on hard inputs | Steiner and element budgets, one bounded retry, stuffing fallback, oracle comparison to measure the gap |
| Exact predicates slow the Delaunay kernel | Adaptive predicates take the fast float path in most cases; profile before optimizing |
| Self-intersection tests slow on large meshes | BVH culling, parallel tests, a per-file time budget that reports partial results |
| Experimental code panics or hangs | Worker process with timeouts and memory caps |
| Exodus side-set ordinals or names subtly wrong | Per-face unit tests, round-trip reads, SEACAS checks |
| Unknown units make absolute tolerances meaningless | All tolerances relative to the bounding-box diagonal; optional `--units` |
| Repairs override intent (an opening filled shut) | Size thresholds on hole filling; large changes become errors |
| GUI absorbs the project | GUI drives the worker only; it can wait until diagnosis is solid |

**Open questions.**

- [x] Which downstream solver will consume the meshes: finite element or finite volume? It decides output formats and which quality metrics matter most.
- [ ] Single region only, or multi-region meshes where nested or touching shells become separate material zones?
- [ ] Should open sheets be rejected, or optionally closed by thickening?
- [ ] What mesh-size controls should users get in v1: a global target size, or per-region overrides too?
- [x] Does the finite-volume solver read Exodus II directly, or does it need OpenFOAM polyMesh, CGNS or Fluent .msh as well?
