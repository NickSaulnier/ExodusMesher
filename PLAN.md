# OBJ Automatic Mesher — Project Plan

Nov 7, 2026 · @Nicholas Saulnier

## Goals and scope

Build, entirely in Rust, a mesher that turns any .obj file into a validated surface mesh and tetrahedral volume mesh, with its own volume mesher as an experiment and a desktop GUI for choosing which parts to mesh. Every input ends in a valid mesh or a report explaining exactly why not.

**Input class.** OBJ files as they occur in the wild: triangles, quads and n-gons; multiple objects and groups; duplicate vertices; degenerate faces; holes; non-manifold edges; self-intersections; inconsistent orientation; several or nested shells; unknown units. Normals, texture coordinates and materials are read and ignored.

**Parts.** Each `o` object or `g` group is a part; a file with neither is split into connected components. The GUI and the CLI's `--parts` flag select which parts to mesh.

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
| Linear algebra | `nalgebra` | Mature, well-typed small vectors and matrices |
| Spatial queries | Own BVH, or `parry3d` | Intersection, nearest-point and ray queries |
| Inside/outside | Generalized winding number, own implementation | Stays meaningful on surfaces with holes and self-intersections |
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
| `MESH_WORKER_FAILED` | Worker panicked, timed out or ran out of memory | Error, triggers fallback |
| `VOL_INVERTED_TET` | Tetrahedron with non-positive volume | Error |
| `VOL_LOW_QUALITY` | Non-orthogonality, skewness or volume ratio out of bounds | Warning or error by threshold |
| `VOL_BOUNDARY_DEVIATION` | Hausdorff distance to input over tolerance | Error |
| `EXO_NAME_TRUNCATED` | Part or set name over the Exodus limit | Truncate, warning |
| `EXO_CAPACITY_EXCEEDED` | Mesh too large for the 64-bit-offset format | Error |

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
