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

**Non-goals for v1.** Boundary-layer prisms, hex meshing, CAD (STEP/IGES) input, TETRA10 elements, conformal interfaces between touching parts, distributed parallel meshing, and physics-driven sizing. Each is a candidate for later phases; [ROADMAP.md](ROADMAP.md) places them.

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
   yes -> Sizing field -> Delaunay mesher -> [Volume valid?]
            yes -> Write mesh + report (Exodus, .msh, .vtu, report.json)
            no  -> Remesh surface -> Delaunay mesher -> [Volume valid?]
                     yes -> Write mesh + report
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
6. **Sizing.** A sizing field from curvature and local thickness grades the target size. If the conformal mesh fails, the surface is remeshed isotropically to the field (features and group boundaries kept), the surface gate checks it again, and the Delaunay mesher retries once.
7. **Delaunay mesher.** The in-house constrained Delaunay kernel, run inside the worker process.
8. **Volume gate.** Positive orientation, quality bounds, volume match, boundary conformity.
9. **Emit.** Exodus II, .msh and .vtu, plus `report.json`, the defect overlay and the stdout summary.

**Process model.** The CLI and GUI are thin front ends. Both launch `mesher-worker`, which runs one job and streams JSON lines (progress, report records, result) on stdout. Panics, runaway refinement and hangs stay inside the worker, and cancelling is just killing it.

**Workspace layout.**

```text
crates/
  mesher-geom/    predicates, half-edge mesh, BVH, winding number
  mesher-core/    parts, diagnose, repair, gates, sizing, report types, project files
  mesher-tet/     volume kernels: stuffing, Delaunay, recovery, refinement
  mesher-io/      OBJ reader; Exodus, .msh, .vtu, .vtp writers
  mesher-worker/  binary: one job, JSON-lines protocol on stdout
  mesher-client/  worker supervision shared by the front ends
  mesher-cli/     clap front end, exit codes, `mesher run project.toml`
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
2. **Boundary recovery.** Recover missing subsegments by Steiner points, then missing facet triangles by adopting the mesh's own facet triangulation, swapping facet diagonals, flipping away crossing edges, and Steiner points as a last resort. Whatever is recovered is constrained, so no later insertion deletes it (a constrained Delaunay tetrahedralization; see the Milestone 8 decisions). Lawson flips then restore the constrained Delaunay property. A point budget bounds it; an unrecoverable facet becomes `MESH_RECOVERY_FAILED`, and the job falls back to Kernel A.
3. **Region marking.** Remove exterior tets and cavities by a parity flood fill from the hull across the recovered surface triangles.
4. **Quality refinement.** Constrained Delaunay refinement on radius-edge ratio and size, then sliver removal by flips (edge removal and 2–3 flips that improve the worst tet) and two rounds of circumcentre insertion, with refinement capped by the point budget. Vertex smoothing is deferred.

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
| `PARSE_PROJECT_INVALID` | Project file cannot be read or is not a valid project | Fatal |
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
| `MESH_ROUNDING_REPAIRED` | A tet that rounding to f64 turned over or flattened was untangled | Untangle, warning |
| `MESH_ROUNDING_FAILED` | Tets that rounding to f64 turned over or flattened could not be untangled | Error, triggers fallback |
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
2. **Remesh, then retry once.** If the conformal mesh fails, remesh the whole part's surface to the size field (keeping sharp edges, corners and group boundaries) and run the Delaunay mesher on it once more. The mesh then conforms to the remeshed surface, and the volume gate measures its deviation from the repaired one. Bounded to one retry so it cannot loop; `--no-remesh` skips it.
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
11. **Indirect predicates.** Points defined by construction (on a segment, on a triangle, line–plane and three-plane intersections), with exact predicates on any mix of them.
    - Exit: every predicate agrees with an exact rational reference on random, degenerate and nearly degenerate inputs; the filter settles at least 95% of random queries; property tests at 100,000 cases each.
12. **Exact constrained Delaunay kernel.** `mesher_tet::cdt` beside the current kernel, after Diazzi et al. (2023).
    - Exit: in exact form, every constraint recovered on the synthetic suite, the 100-model clean subset and the rounded UV sphere; a disjoint-cavity case reproduced.
13. **Rounding to f64.** Round the exact mesh, check every tet exactly, untangle what flipped.
    - Exit: rounding validity measured before and after untangling; no non-positive tet ever written; every untangle a record.
14. **Chamfer-and-refine.** Quality refinement that terminates, after the surface-chamfering paper (2026).
    - Exit: refinement ends within its budget on every clean-subset model or records why; angles and cost tabulated.
15. **Worst-cell quality for finite volume.** Invariant-guarded local operations aimed at non-orthogonality and skewness at the wall.
    - Exit: worst wall-adjacent non-orthogonality and skewness below both the current kernel's and TetGen's; byte-identical at 1, 4 and 16 threads.
16. **Exact arrangement (repair rung 2).** Self-intersections and overlapping shells resolved exactly, cells kept by winding number.
    - Exit: every part of the self-intersecting and multi-component strata resolves to a surface that passes the gate or fails with a located record.
17. **Wrapping (repair rung 3).** Alpha wrapping for what the arrangement cannot fix.
    - Exit (and the phase's): the full benchmark on all 10,000 models and the valid subset; outputs checked by `mesher-fvm` and by the user's solver.

## Phases and task breakdown

The milestones group into seven phases. Phases run in order, except the GUI phase, which can start once Phase 2 is done. A milestone's tasks are ticked when they land; the milestone is done only when its exit criterion passes under `cargo fmt`, `cargo clippy -D warnings` and `cargo test --workspace`.

| Phase | Milestones | Delivers |
| --- | --- | --- |
| 1. Reporting spine | 1 | Workspace, report schema, parser, worker, CLI: any file in, a valid report out |
| 2. Geometry and diagnosis | 2, 3 | Exact geometry core, part splitting, the full defect taxonomy, synthetic suite |
| 3. Repair and output | 4, 5 | Logged repairs and surface gate; Exodus, .msh and .vtu writers |
| 4. Volume meshing | 6, 7, 8 | Stuffing fallback end to end, then the constrained Delaunay kernel |
| 5. GUI | 9 | egui front end over the worker (can start after Phase 2) |
| 6. Hardening | 10 | Sizing, Thingi10k benchmark, CI regression gate |
| 7. Exact kernel and repair ladder | 11–17 | Indirect predicates, exact constrained Delaunay kernel, f64 rounding, terminating refinement, finite-volume quality, exact arrangement, wrapping (Phase 1 of `ROADMAP.md`) |

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
- [x] Octree grading (adaptive lattice): done in Milestone 10 (graded interior, uniform surface, as the paper proves)

**Decisions made in Milestone 6.**

- *All-zero tets are dropped.* A lattice tet whose four vertices were all warped onto the surface is a sliver lying in it; keeping it (as an early version did when its centroid was inside) produced dihedral angles down to 2.8°. Dropping it leaves its faces, which lie on the surface, as boundary. With that, the empirical range on all tests stays inside the paper's [10.7°, 164.8°].
- *Stencil filling by coning.* Each lattice tet's inside part is filled by a cone from its lowest-numbered inside vertex over its boundary triangles; a quadrilateral (only on faces with two inside vertices and one outside) is split from the lower-numbered inside vertex. Both rules depend only on the shared face, so neighbours conform. This replaces the paper's colour-based parity rule; the property tests above are the evidence that the angle bounds still hold.
- *Deviation is measured against the region's boundary.* The Hausdorff distance is two-sided between the mesh boundary and the input surface, sampled at vertices and face centroids, but input samples with winding number ≥ 1 (faces buried inside the region, as with overlapping shells) are skipped.
- *Defaults.* Element size 0.03 × part diagonal; maximum deviation 0.05 × part diagonal (raised from 0.01: stuffing rounds sharp corners, e.g. 1.0% on a lattice-aligned cube and 3.5% on a tetrahedron at the default size).
- *Exit code 20 for every meshed run* until the Delaunay kernel lands: stuffing is by definition the fallback. (Superseded in Milestone 8: a conformal mesh exits 0 or 10.)
- *Known gap:* a feature thinner than the element size is lost (the synthetic spike case is rejected at 15% deviation); the fix is a smaller element size or, later, grading.

**Milestone 7: Delaunay kernel.** Done (Nov 7, 2026): on random point sets of 5 to 20,000 points, and on clustered point sets spanning three orders of magnitude in scale, the tetrahedralization matches TetGen's exactly (same count, same tets).

- [x] `mesher_tet::delaunay`: flat tet array with four vertices and four neighbours each, free list, ghost vertex at infinity (always stored last), positively oriented real tets
- [x] Bowyer–Watson insertion with exact predicates; ghost conflict: beyond the hull face, or coplanar and strictly inside its circumcircle; visibility walk with deterministic pseudo-random face order (scan fallback)
- [x] BRIO: seeded shuffle, rounds of doubling size, Morton order within each round
- [x] Incremental `insert` for later Steiner points; `check()` for structural invariants; exact duplicates reported, not inserted; all-coplanar input refused (`DelaunayError::Flat`)
- [x] Property tests against brute force: structure valid, every distinct point a vertex, no point strictly inside any circumsphere, faces conform, hull closed, outward and convex. Inputs: random points, integer grids, exactly cospherical points (all 30 integer points with x² + y² + z² = 25, with and without the centre), cube-surface grids, nearly flat sets (6,000 cases each run once)
- [x] TetGen comparison in `tests/oracle/tetgen_delaunay.rs`

**Decisions made in Milestone 7.**

- *Exact coplanar in-circle without projection.* `in_circumcircle_coplanar` uses `insphere` with an off-plane point: every sphere through three points meets their plane in their circumcircle. A first version projected onto a coordinate plane, which keeps orientations but turns circles into ellipses; the cube-grid property test caught it. A regression test on the tilted plane x + y + z = 3 guards it.
- *Walks never start on a ghost.* The last-created tet can be a ghost; starting the walk there located points in the wrong tet (found by inserting a cube's centre). The walk now steps into the hull first, and a cavity is never seeded from a ghost the point does not conflict with.
- *Cavity growth for degeneracies.* If a cavity boundary face cannot see the new point (only possible with cospherical or coplanar points), the cavity grows across it; the result stays a valid triangulation, and on the degenerate test sets it stays Delaunay.
- *Performance.* About 2.3× TetGen's time in release builds (146 ms against 63 ms for 20,000 points). Profiling is deferred to Milestone 10, as the risk table says.

**Milestone 8: Recovery and refinement.** Landed (Nov 7, 2026) with its exit criterion open by decision: the pipeline meshes the clean Thingi10k models conformally or by the recorded stuffing fallback, and Milestone 9 starts. The gaps are listed under "Things to consider improving" at the end of this plan. On synthetic shapes every mesh is conformal, with dihedral angles in [10°, 165°] on the shape tests and a minimum of 9.8° on the oracle's finer torus, well above TetGen's 5.5°. On the Thingi10k subset 54 of 100 clean models mesh conformally (67 with the earlier conforming kernel), with much better quality than before (see the benchmark).

- [x] `mesher_tet::conforming`: constrained Delaunay mesh of a closed surface (the module keeps its name). Recovered subsegments and subfacets are constraints. Boundary Steiner points go in by 2D Bowyer–Watson on the facet, with the tets on both sides of the replaced subfacets seeding the 3D cavity, so the new subfacets appear with the point. Subsegments touching an input vertex split on concentric shells
- [x] Kernel constraints (`mesher_tet::delaunay`): constrained faces and edges. Cavities grow by conflict without crossing a constraint, then shrink (the tet holding the point goes last) until every boundary face sees the point and no constrained face, constrained edge or vertex would be deleted; otherwise the insertion is refused and nothing changes (`DelaunayError::Blocked`). A forced insertion ignores face constraints and reports the faces it deleted, which are recovered at once
- [x] Plane labels: every vertex records the facets it lies on, and points sharing a facet are coplanar by construction. No tet is ever made flat within a facet, and a point never "sees" a face on its own facet. The initial Delaunay tetrahedralization is built with them
- [x] Flips: edge removal (3 to 7 tets around an edge, every fan tried, validated by face orientations before anything changes), the 2–3 flip, and Lawson flips back to constrained Delaunay after recovery. Unit tests check validity and volume
- [x] Facet recovery by flips first: adopt the mesh's own triangulation of a facet when it is complete, swap a facet diagonal to the mesh's, flip away edges crossing a missing subfacet, and use Steiner points only where these fail. A recovery that stops making progress is `MESH_RECOVERY_FAILED`
- [x] Constrained refinement on radius-edge ratio 2 and circumradius ≤ size. A point encroaches a subfacet only inside its diametral sphere and projecting into it; a subfacet between a tet and its circumcentre is split; tets whose shortest edge is below the size floor are not refined for shape. Sliver removal: improving flips, two rounds of circumcentre insertion, flips again
- [x] Region marking by parity flood fill from the hull; a conflict is `MESH_RECOVERY_FAILED`. Tets flat to rounding between two subfacets of one OBJ group are flattened into the boundary when both pairs of faces cover the same area
- [x] Point budget sized from the surface (20 × (volume/size³ + area/size²), at least 10,000) as `MESH_BUDGET_EXCEEDED`; `ConformingParams::min_boundary_size` (size/16) floors refinement-driven boundary splits
- [x] Pipeline: conformal first for every part that passed the surface gate, stuffing as the recorded fallback (`MESH_RECOVERY_FAILED` or `MESH_BUDGET_EXCEEDED` as a warning with `fallback/applied`, or an error with `fallback/failed`); `--strict` meshes conformally or fails with `MESH_STRICT_UNMET`; OBJ groups become side sets
- [x] Finite-volume quality (`mesher_tet::fv_quality`): non-orthogonality, skewness, neighbour volume ratio, minimum dihedral, 18-bin dihedral histogram; all in the report's mesh summary and schema; Delaunay meshes gate at [10°, 170°] (`VOL_LOW_QUALITY`, a warning)
- [x] `mesher-fvm`: cell-centred Poisson solver (over-relaxed non-orthogonal correction, least-squares gradients, Jacobi-preconditioned CG) on a manufactured solution; L2 error on conformal cube meshes 1.7e-2, 3.2e-3, 7.9e-4 at sizes 0.25, 0.125, 0.0625
- [x] Synthetic suite: 86 of 98 part meshes conformal, every one exact in volume; touching shells fall back to stuffing as expected. The synthetic shapes (`crates/mesher-tet/tests/conforming_shapes.rs`) include a tilted cube and a torus whose quads are planar only to rounding
- [x] TetGen oracle (`tests/oracle/tests/tetgen_quality.rs`): cube, L-prism, icosphere and torus at the same size targets. Ours keeps every dihedral in [9.8°, 165.6°] (10.0° except one boundary sliver on the torus), against TetGen's minimum of about 5.5°, with lower mean non-orthogonality and a similar cell count (`target/tmp/oracle/quality_vs_tetgen.txt`)
- [x] Facets: edge-connected, same-group triangles facing the same way that are coplanar, exactly or to rounding (their four points make a tet flat to rounding), form one facet. Any triangulation of a facet is allowed; only edges between facets are subsegments
- [x] Thingi10k subset: 100 clean models (closed, manifold, oriented, one component, no self-intersections, solid, at most 20,000 faces; evenly spaced by id from 3,701), ids in `tests/thingi10k/subset.txt`; meshes downloaded to `target/thingi10k/` (not committed: Thingiverse licences). Harness `tests/oracle/tests/thingi10k.rs` (ignored; STL welded to OBJ, the `mesher` CLI end to end, TetGen per model in a child process with a 120 s limit)
- [ ] Exit criterion: conformal meshes on the clean part of a Thingi10k subset, with quality histograms against TetGen. Deferred by decision (Nov 7, 2026); see "Things to consider improving"

**Thingi10k benchmark (Nov 7, 2026, default settings, release build).**

| Measure | Constrained kernel (now) | Conforming kernel (earlier the same day) | TetGen (`pzAfq`, same size) |
| --- | --- | --- | --- |
| Meshed | 98 (54 conformal, 44 stuffing fallback; 2 worker timeouts with no mesh) | 100 (67 conformal, 33 stuffing fallback) | 93 (2 timed out at 120 s, 3 crashed) |
| Exit codes | 0: 6, 10: 48, 20: 31, 30: 13 | 0: 9, 10: 58, 20: 25, 30: 8 | |
| Per-model minimum dihedral, median | 10.1° | 5.6° | 3.4° |
| Dihedral histogram (% per 10° bin) | 0.5 2.6 5.1 7.3 9.4 15.1 14.8 9.6 10.7 10.5 5.0 3.6 2.4 1.6 1.0 0.6 0.1 0.1 | 0.6 2.2 5.3 8.2 10.6 13.7 13.4 10.7 10.4 9.4 5.5 3.9 2.6 1.6 1.0 0.6 0.1 0.1 | 0.6 2.4 6.1 8.2 10.5 12.2 12.3 11.6 10.9 9.0 6.0 4.1 2.6 1.6 1.0 0.7 0.2 0.0 |
| Median mean non-orthogonality | 20.7° | 25.0° | 26.6° |
| Total time | 2,680 s | 3,214 s | 140 s |

On the first 48 models the constrained kernel newly meshes 3 conformally (34785, 45009, 112922; 45009 at 277k tets against TetGen's 285k) and loses 10 that the conforming kernel meshed (40179, 42936, 43853, 57657, 64992, 69403, 73154, 87721, 108330, 118032). Most are recovery stalls on facets whose neighbours are coplanar only to float32 noise; two are point-budget overruns. 91607 and 677341 (about 420 faces each) loop until the worker timeout.

**Decisions made in Milestone 8.**

- *Constrained, not conforming (revised the same day).* The first kernel recovered the boundary by Steiner points alone and kept subfacets Delaunay (Lawson flips, and splitting every subfacet a vertex encroached). On Thingi10k every fallback was a point-budget overrun: recovery cascades, where a split destroyed what an earlier one had recovered, and refinement splitting thin subfacets with huge diametral spheres. Constraints fix both: what is recovered stays, and subfacets no longer have to stay Delaunay.
- *Symbolic coplanarity.* Steiner points on a facet are off its plane by rounding, and facets coplanar only to rounding are common in synthetic f64 data. Exact orientation then decides visibility at random, and flat tets block every split near them. Plane labels make "on the same facet" a topological fact, so orientation stays exact everywhere else (no epsilon).
- *Facets coplanar to rounding merge.* Like TetGen's coplanarity tolerance, but with the rounding-relative test that flattening already used (64 ulps of the coordinates, or 1e-12 of the size), so only rounding is absorbed. Float32 noise on CAD tessellations is far larger, and those stay separate facets.
- *Flips before points.* Facet recovery adopts the mesh's own facet triangulation, swaps diagonals and removes crossing edges before any Steiner point. Lawson flips restore constrained Delaunay afterwards; without them, refinement's cavities had to grow and shrink unpredictably.
- *Cavity growth across unseen faces.* Plain Bowyer–Watson grows across a face the point cannot see (cospherical points). In a constrained mesh such faces are common, and unbounded growth ran away (cavities of 8,000 tets). Growth is now allowed only without constraints, for forced insertions, into conflicting tets, or locally (up to 4× the initial cavity) for interior points; a cavity beyond 4,096 tets is refused.
- *Encroachment by projection.* A subfacet is encroached only by a point inside its diametral sphere that projects into it. With the full sphere, almost every circumcentre near a thin input triangle was rejected and no interior point was ever inserted.
- *Parity region marking.* Winding-number classification of tet centroids misclassified rounding-thin tets on the icosphere. The flood fill is exact and topological.
- *Unreliable circumcentres.* A near-flat tet's float circumcentre can be far off, so a centre counts only if the exact `insphere` agrees; otherwise the tet's longest subsegment is split. Each tet gets at most eight attempts.
- *Small input angles.* Segments meeting at small angles encroach each other forever (classic Ruppert). Concentric-shell splitting fixes moderate angles (the L-prism's 18° fan); the size floor bounds the rest.
- *Flatness is relative to coordinates.* Tets flat to rounding are judged against 64 ulps of their coordinates as well as their size, since far from the origin rounding is much larger than the tet.
- *TetGen runs in a child process.* It looped on two benchmark models and crashed on three; in-process, that stalled the harness.
- *Touching shells fall back.* Repair separates shells that share an edge or a vertex by duplicating it; the conformal kernel refuses coincident vertices, so these parts use stuffing (recorded). Merging them is deferred.
- *Exit codes.* A clean cube now exits 0; parse defects exit 10; the stuffing fallback exits 20.
- *Exit criterion deferred.* By the user's decision (Nov 7, 2026), Milestone 9 starts with 54 of 100 conformal; the improvements are listed at the end of the plan.

### Phase 5: GUI

**Milestone 9: GUI.** `mesher-gui` with wgpu viewport, part tree, picking, overlay, run panel, export, project files shared with the CLI. Exit criterion passes headless (Nov 7, 2026); the interactive check is open. Run it with `cargo run -p mesher-gui --release -- model.obj` (build `mesher-worker` too: the GUI runs the worker next to it).

- [x] Project files: `mesher_core::project::Project` (input, output directory, job settings) as TOML, paths relative to the project file; `mesher run project.toml` runs it, with command-line flags overriding; an unusable project file is `PARSE_PROJECT_INVALID`
- [x] `mesher-client`: worker supervision (spawn, stream messages, timeout, memory cap, cancel, crash and timeout reports) shared by the CLI and the GUI
- [x] `mesher-gui` shell: eframe 0.37 app; open an OBJ or a project (file dialog, drag and drop, or a command-line argument); the model is read and diagnosed with the worker's own code on a background thread, so parts, records and the overlay match the worker exactly (no meshing in the GUI)
- [x] Viewport: wgpu render pass into an offscreen colour and depth target shown as an egui image; orbit, pan, zoom, double-click to fit; flat shading with a headlight plus wireframe; one colour per part, grey when unselected; a clipping plane along x
- [x] Part tree: name, triangle count, bounding-box size, defect badges by severity, selection checkbox, all and none
- [x] Picking: click a part to toggle it (a ray from the camera against every display triangle; box and lasso selection later)
- [x] Defect overlay: diagnosis overlay faces in red, edges as red lines, vertices as markers; a record list where clicking one highlights its faces and vertices in yellow and shows its code, hint, measurements and OBJ lines
- [x] Run panel: element size, maximum deviation, strict and no-repair; Run streams progress into a bar and records into a log; Cancel kills the worker (a `MESH_WORKER_FAILED` report with reason `cancelled`)
- [x] Export: output directory chooser; the worker writes Exodus, .msh, .vtu and the report; `project.toml` is saved beside them (and "Save project..." saves one anywhere)
- [x] Exit criterion: select parts of a multi-object OBJ and export an Exodus file with one block per part, reproducible from the CLI. `crates/mesher-gui/tests/exit_criterion.rs` picks a part with a camera ray, toggles another in the session, exports through the worker (two blocks of three parts), and replays the saved project with `mesher run project.toml` to the same blocks
- [ ] Interactive check of the window by the user (the automated tests cover the logic, not the drawing)

### Phase 6: Hardening

**Milestone 10: Sizing, benchmark and hardening.** Sizing field and isotropic remeshing, Thingi10k run, profiling, CI regression subset with golden reports, documentation. Started Nov 8, 2026.

Scope decisions (Nov 8, 2026, by the user):

- *Benchmark set.* A stratified sample of 1,000 Thingi10k models rather than all 9,997 (438 M faces, tens of GB): sampled across face counts and defect classes (clean, holes, non-manifold, self-intersecting, several components), capped at 100,000 faces, so every pipeline path is exercised in a few hours on one machine.
- *CI.* GitHub Actions runs `cargo fmt --check`, `cargo clippy -D warnings`, `cargo test --workspace` and golden reports for the synthetic suite. The Thingi10k regression subset stays a local, ignored test: the meshes cannot be committed (Thingiverse licences).
- *SEACAS.* Installed in WSL Ubuntu from conda-forge; the deferred Milestone 5 check runs first.

- [ ] SEACAS `io_info` and `exodiff` on the oracle and pipeline Exodus files (closes the Milestone 5 item). Deferred by the user (Nov 8, 2026). SEACAS is not on conda-forge or apt, so it has to be built from Sandia's source in WSL (`fmt` from conda-forge, system netCDF and HDF5)
- [x] Work budget: a deterministic step count in the constrained kernel (insertion attempts, refinement and repair iterations, recovery passes), 64 steps per budgeted point, `ConformingError::Work` → `MESH_BUDGET_EXCEEDED`. A local repair cascade that has grown three rounds running past 256 missing subfacets is handed to a full recovery pass at once (91607 now falls back in about 100 s instead of timing out)
- [x] Worker failure triggers the fallback, as the taxonomy says: a worker that times out or runs out of memory is rerun once by `mesher-client` with `stuffing_only` (not under `--strict` or `--no-repair`); the report leads with the first run's `MESH_WORKER_FAILED` record, a warning with `fallback/applied` when the retry meshed (677341: no mesh before, stuffing mesh now)
- [x] Sizing field (`mesher_tet::sizing`): curvature (about 25° of turning per element, i.e. 0.44 × the radius of curvature, on edges turning by less than the 45° feature angle) and thickness (a ray cast inward from each face; no element larger than the local thickness), clamped to [h/4, h], graded at 0.4 through a background grid (Dijkstra over the 26-neighbourhood) and interpolated trilinearly. The constrained kernel refines each tet against the size at its centroid, and the point budget is estimated from the field. On by default; `--uniform-size` (and a GUI switch) meshes at h everywhere; the report's mesh summary gains `min_element_size`
- [x] Isotropic surface remeshing (`mesher_tet::remesh`): split above 4h/3, collapse below 4h/5 (link condition, no flipped triangles), Delaunay edge flips, tangential smoothing; every vertex projected back onto the input (smooth vertices onto the nearest input triangle facing their way, feature vertices along the input's feature lines, corners fixed); feature lines are edges turning by more than 40° and group boundaries, so side sets survive. Used as the retry when the conformal mesh fails (fallback ladder step 2): the remeshed surface is gated again, meshed with each triangle its own facet, and judged by the volume gate against the repaired surface. `MESH_RECOVERY_FAILED` becomes a warning with `remesh/applied`; the attempt and the mesh summary carry `remeshed`; `--no-remesh` (and a GUI switch) skips it
- [x] Graded stuffing lattice (`mesher_tet::graded`, `stuff_graded`), following Labelle and Shewchuk's Section 5 so the proven bounds hold: leaf cubes found by continuation from surface seeds, the Continuation Condition and the half-pyramid guard, an octree with partial children up to element-size roots, Weak Balance, and the Figure 10 tetrahedra (BCC, bisected, quadrisected, half-pyramids). The cutting, warping and stencil code is shared with the uniform lattice through a `Background` trait. In the pipeline, a part whose size field asks for less than h is stuffed with surface cells of h/2 or h/4 grading to h inside; uniform fields keep the plain lattice. Property test: 500 random graded shapes all within [10.7°, 164.8°], conforming and closed; a graded sphere has the uniform one's volume to 1e-9 with less than half the tets. 472131's fallback: deviation 4.24% → 0.61%, 2,012 → 82,958 tets
- [x] Stratified 1,000-model sample (`tests/thingi10k/sample1000.txt`), download script, parallel harness with per-model timeouts, outcome table by path against TetGen (the tooling; the run itself is part of the exit criterion)
- [x] Profiling of the constrained kernel (`mesher_tet::prof`: timing scopes behind the `profile` feature, compiled out otherwise; `cargo run --release -p mesher-tet --features profile --example profile -- model.obj [--max-points N]`), and fixes for the hot spots, all leaving results bit-identical (golden reports unchanged; 43853 gives the same 92,009 tets): inside tests cached per tet; FxHash for the Delaunay core's constraint sets; the cavity loop keeps only the first tet to take out, keeps visibility, exclusion and conflict answers in flat per-slot arrays stamped by insertion, and updates its list of unseen boundary faces in place after a removal instead of rebuilding it. Kernel time on 43853 8.4 s to 5.5 s; 677341 to 8,000 points 26.4 s to 3.8 s (it now ends at its point budget in 97 s instead of timing out); end to end, 472020 and 135674 (failure, then remesh) about 50 s instead of 3 minutes
- [x] Sliver repair by vertex smoothing (`mesher_tet::smooth`): after the constrained kernel, every vertex on no boundary face that belongs to a tet below 15° moves to the best of a few candidates (towards its neighbours' centroid, along a finite-difference gradient of the worst angle, away from the opposite face of its worst tet), kept only if every tet around it stays positive (exact `orient3d`) and the worst angle improves; four passes. Connectivity and boundary never change, so conformity and volume are exact. Cells below 1°: 43853 131 → 50, 45009 524 → 42, 113890 132 → 64; synthetic suite minimum dihedral up 2–5° on every conformal part (golden regenerated; mean non-orthogonality up by at most 0.6°). About 8% of kernel time (a one-`acos` minimum dihedral keeps it there)
- [x] Golden reports for the synthetic suite (`tests/synthetic/golden.txt`, 98 runs): records by code, severity and action with counts, and each part's kernel and acceptance, must match exactly; tet count within 10%, minimum dihedral and mean non-orthogonality within 2°. `UPDATE_GOLDEN=1` regenerates it
- [x] GitHub Actions workflow (`.github/workflows/ci.yml`: fmt and clippy on Linux; `cargo test --workspace --release` on Linux and Windows); green on GitHub (confirmed by the user, Nov 10, 2026)
- [x] Documentation: README (install, CLI, GUI, project files, exit codes, report, benchmark table) and the report reference (`docs/report.md`: every field, defect code, meshing measurement, quality metric and output file; a test in `mesher-core/tests/schema.rs` keeps its code table in step with `DefectCode`)
- [x] Stratified benchmark run (Nov 8, 2026, at the user's request): the outcome table by path against TetGen is below. The run found three bugs, fixed before 315 affected models were rerun: errors that the fallback worked around still exited 30; a tet flat in floating point made the volume ratio infinite, which JSON writes as `null`, so the client could not read the result and reported a dead worker; and the harness counted TetGen's empty meshes as meshed
- [x] Exit criterion: published outcome table by path, against TetGen (done, below), with CI green (confirmed Nov 10, 2026)

**Stratified Thingi10k benchmark (Nov 8, 2026, 1,000 models, default settings, release build).** 685 models from the first run and 315 rerun after the fixes above. 12 meshing jobs at a time on the first run and 8 on the rerun. Each worker run is limited to 600 s (a timeout is retried once with stuffing only) and 3 GiB of memory; TetGen (`pzAfq`, same size) runs in a child process limited to 120 s.

| Class | Models | Conformal, no repair (exit 0) | Conformal after repair (10) | Stuffing fallback (20) | Failed (30) | TetGen meshed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Clean | 350 | 22 | 282 | 40 | 6 | 323 |
| Multi-component | 150 | 6 | 108 | 32 | 4 | 126 |
| Non-manifold | 120 | 0 | 59 | 57 | 4 | 72 |
| Open | 80 | 0 | 65 | 11 | 4 | 19 |
| Self-intersecting | 300 | 0 | 26 | 235 | 39 | 3 |
| All | 1,000 | 28 | 540 | 375 | 57 | 543 |

- *Meshed:* 943 of 1,000 (568 conformal), against TetGen's 543. No internal errors (exit 70) and no missing reports. TetGen returned no tetrahedra on 417 models (most of the self-intersecting and open ones), crashed on 21, timed out on 17 and failed with an internal error on 2. 48 models fail in both; 9 mesh only in TetGen. In 5 of those, conformal recovery failed and the fallback was more than 5% off the surface (52134, 130986, 593082, 1436274, 1489676). The other 4 meshed conformally in the first run, but with tets flat in floating point. The new flat-tet gate rejects those meshes, so the models go down the ladder, and they timed out under load (76541, 500100, 500110, 500120). Alone, 500100 took 130 s and got a stuffing mesh of 1.09 M tets where the rejected conformal one had 107 k. Since the kernel fixes of Nov 9 all four mesh conformally again (500100: 106,886 tets in 22 s); see the findings below.
- *Quality on the 534 models both meshed* (per-model minimum dihedral, median): ours 2.98° (the stuffing fallback's 12.3° lifts it), ours on the 452 conformal meshes 0.90°, TetGen 0.72°; median mean non-orthogonality 23.2° against 27.2°. The 5th percentile is 0° for both: some cells are flat to a fraction of a degree. Dihedral histogram over every accepted mesh (% per 10° bin): 0.1 0.7 2.1 3.4 11.8 19.1 19.8 6.1 9.7 18.8 2.8 2.4 1.7 0.7 0.4 0.2 0.0 0.0 (the stuffing lattice's angles dominate it).
- *Failures* (models carrying each code): `QUAL_SLIVER_FACE` 43, `VOL_BOUNDARY_DEVIATION` 38 (the fallback mesh is further than 5% of the diagonal from the surface), `INT_SELF_INTERSECTION` 30, `TOP_BOUNDARY_EDGE` 20, `MESH_WORKER_FAILED` 12 (all timeouts), `TOP_OPEN_SHEET` 11, `MESH_RECOVERY_FAILED` 9, `VOL_LOW_QUALITY` 7.
- *Time:* median 37.8 s a model, 90th percentile 255 s, about 15× TetGen in total. The rerun's models took a median 2.2× their first-run time: the rerun ran the heaviest models (mostly self-intersecting fallbacks) side by side.

**Findings and decisions in Milestone 10.**

- *The two "hangs" were not loops.* 91607 was a repair cascade: a forced boundary split deleted 3 subfacets, and each round of local repair deleted more than it recovered (3, 9, 16, 21, 60, 100, 216, 327, 478, 658, ...), each round slower. 677341 refines normally but slowly (about 150 points a second against a budget of 88,720 points; boundary splits and constrained insertions dominate), so it ran into the worker timeout. Cascades are common and usually harmless: on the converging models local repair grows its queue every round for all 16 rounds (to 300–1,300 subfacets), and the full recovery pass that follows settles it. 91607 never got there because its rounds became too slow (refused insertions). Failing a cascade cost conformal models (13 with a tight threshold, 8 at 400), so a growing cascade is instead handed to full recovery early. The handover point is a trade-off and not monotonic: at 64, 91607 falls back in 2 s but 34785 and 148292 lose their conformal meshes; at 256, only 34785 does. 256 was kept; the 100-model subset has not been rerun with it.
- *Determinism versus timeouts.* The kernel's own limits count steps, not seconds, so a given input and config always end the same way. The worker timeout is the only wall-clock limit, and the stuffing retry after it is recorded, so a slower machine can change the path but never silently.
- *A small reproducer for improvement 3.* A 16 × 10 UV sphere meshes conformally, but the same sphere off-centre with coordinates rounded to six decimals (as an OBJ file stores it) fails recovery (47 subfacets missing then; with later kernel changes it still fails). It meshes conformally once remeshed (`rounded_uv_sphere_meshes_after_remeshing`). It is the ignored test `rounded_uv_sphere_meshes_conformally` in `crates/mesher-tet/tests/conforming_shapes.rs`.
- *Sizing defaults are a trade-off.* With 15° and a floor of h/8, every graded part hit the floor and three conformal models doubled in cells and time; with 25° and h/4, 43853 drops from 315k to 92k tets (and from 31 s to 12 s), while 45009, 113890 and 148292 grow by 36%, 13% and 9% where they curve or thin out. On the synthetic suite only the tori change (about 70% more cells, slightly better non-orthogonality). The field does not fix the point-budget overruns (57657, 73154): those come from shape refinement chasing small input features, not from the target size, and a size-relative shape floor did not help either (it cost 148292 its conformal mesh), so it was dropped.
- *Remeshing is the retry, not a first step.* A part whose conformal mesh succeeds keeps its input facets exactly; only failures are remeshed. On models that fell back before, the retry meshes conformally where stuffing did not do well: 472020 (stuffing at 2.1% deviation, now conformal at 0.17%), 472131 (4.2%, now 0.12%), 135674 (a failure at 30%, now conformal at 0.57%), and the demo's UV sphere. The cost is a second kernel run, after the first has failed.
- *Delaunay flips, not valence flips.* Botsch and Kobbelt flip towards valence 6. Those triangles were often missing from the tetrahedralization of their own vertices, and recovery spent 18,000 Steiner points on a 3,000-triangle sphere. With Delaunay flips every remeshed triangle of that sphere is a face of the plain Delaunay tetrahedralization (worst angle 37°).
- *Remeshed triangles are their own facets.* Even then the remeshed sphere ran out of points, because the kernel merged remeshed triangles lying on one input facet (coplanar to rounding) into one facet. Each remeshed triangle as its own facet meshes it in 0.4 s with no boundary Steiner points. The merge does not explain the rounded input sphere's failure (it fails either way), so improvement 3 remains.
- *What smoothing can and cannot reach.* Classifying the cells below 1° (`--features profile` example) showed most have a vertex off the boundary (43853: about 120 of 130; 45009: about 490 of 550), which smoothing can move. The rest stand on degenerate boundary triangles from the input, or have all four vertices on the boundary; see improvement 8.
- *Stuffing grades the interior only (user's decision, Nov 8, 2026).* The paper proves its bounds only for meshes uniformly fine at the surface: grading on the surface has no guarantee better than 1.66°–174.72°. So the fallback's surface takes the field's smallest size everywhere (at most h/4), and only the interior coarsens. Fidelity rises sharply for thin or curved parts at the cost of more cells near the surface.
- *Worked-around errors are warnings.* The exit code is derived from severity counts, so a part that the remesh retry or the fallback meshed still exited 30 for the errors that sent it there (245 of the first run's 300 "failures"). Once a later step's mesh passes the volume gate, the part's errors from diagnosis, repair and the gates become warnings with that step's action (`RecordLog::worked_around`). `REPAIR_GUARD_EXCEEDED` stays an error: it means a repair was refused, not that one failed.
- *Flat in floating point is rejected.* The exact `orient3d` keeps every tet positive, but 16 accepted conformal meshes had tets (1 to 60 each) whose volume computed in floating point is zero or less. A finite-volume solver would divide by that volume. The volume gate now rejects them with `VOL_INVERTED_TET` (`flat_in_float`), the volume ratio skips such neighbours, and every number in the mesh summary is finite (−1 when it is not), because JSON cannot carry inf or NaN. After the rerun no written mesh has a flat tet.
- *Where the flat tets came from (Nov 9, 2026).* The `profile` example now lists every tet flat in floating point. There were three causes.
  - Clusters of interior Steiner points a few ulps apart (69548). Circumcentres of cospherical tets agree to rounding, and each was inserted. A refinement point that is a vertex of the tet containing it to rounding (64 ulps of the coordinates, or 1e-12 of the refined tet's circumradius) is now refused; `Walk::Clear` carries the tet the walk ends in.
  - Rounding-flat tets between two remeshed triangles (904472, and the 500100 family). Flattening only merged triangles of one group, and a remeshed surface makes every triangle its own facet group. `ConformingParams::side_sets` now keeps the OBJ groups for flattening apart from the facet groups (`remeshed_box_flattens_rounding_flat_tets_within_a_side_set`).
  - Needle subfacets along finely split input edges (50659, 45136). Still open: improvement 11.

  Of the 16 models, 13 now mesh conformally with no flat tets, all 4 timeouts are gone, and 116870, 50659 and 45136 use the fallback. On 40 random conformal models from the benchmark, every path is unchanged, but 23 change their tet count. Most changes are under 2%. One, 80657, changes by 64%: the refusal changes the kernel's path, its first attempt now misses one boundary triangle, and the remesh retry meshes it conformally in 15,851 tets instead of 43,873. The synthetic golden reports are unchanged.
- *An unreadable result is a bug, not a dead worker.* A result line that does not parse now becomes `MESH_WORKER_FAILED` with reason `unreadable_result` (exit 70).
- *The same mesh on every platform (Nov 9, 2026).* CI's first Linux run failed the golden reports on the synthetic torus: 1,897 tets against 1,810, and a minimum dihedral of 15.0° against 12.3°. The torus generator used `sin` and `cos` from the platform's maths library, which differs from Windows' in the last bits, and the kernel is sensitive to those bits. The mesher also used the platform's `acos`, `atan2`, `cbrt`, `log2` and `pow` in decisions: dihedral angles for sliver flips and removal, smoothing, the sizing field, remeshing's feature lines, winding numbers and the shell splits. All of them now go through `libm`, a pure-Rust port of musl's maths library (MIT/Apache) that uses only correctly rounded IEEE operations and so gives the same bits everywhere. `clippy.toml` disallows the `f64` methods outside the GUI, whose camera maths never reaches a mesh. Only the torus lines of the golden reports changed. Not yet confirmed on Linux: the next CI run checks it.
- *Benchmark tooling.* `tests/thingi10k/sample.py` writes `sample1000.txt` (350 clean, 300 self-intersecting, 150 multi-component, 120 non-manifold, 80 open; 12.4 M faces); `fetch.py` downloads any id list (two models exist only as OBJ and are meshed from the original file); `run_subset.py` runs the CLI alone in parallel for quick checks; the ignored test `stratified_sample_benchmark` in `tests/oracle/tests/thingi10k.rs` runs the CLI and TetGen in parallel and writes the outcome table and a CSV.

**After v1.** Conformal interfaces between touching parts (shared side sets in Exodus), a polyhedral dual mesh for the finite-volume solver, TETRA10 elements, boundary-layer prisms, a model that predicts failure from surface features before meshing, and an agent loop that reads the report and chooses repairs. The phases beyond Milestone 10 are laid out as Milestones 11 to 28 in [ROADMAP.md](ROADMAP.md); each phase's milestones move into this plan when it starts (Milestones 11 to 17 are Phase 7 below).

### Phase 7: Exact kernel and repair ladder

This is Phase 1 of [ROADMAP.md](ROADMAP.md), moved here when it started on Nov 9, 2026.

**Goal.** A robust, fully automatic tet mesher for messy OBJ input that never fails silently, with the kernel rebuilt on exact predicates.

**Phase exit.**

- On the valid Thingi10k subset, conformal success comparable to the published exact CDT results, reported on both the full and the valid denominators.
- Every failure carries a defect code and a location. Zero crashes, zero silent repairs, and no invalid tet ever written.
- The Exodus files load in the user's solver and pass a manufactured-solution check there and in `mesher-fvm`.

Decisions (Nov 9, 2026, by the user):

- *The new kernel goes beside the old one.* `mesher_tet::cdt` is built next to `mesher_tet::conforming` and chosen per job by a setting. It becomes the default only when it beats the current kernel on the 1,000-model benchmark; the Nov 8 table is the baseline. The golden reports keep guarding the old path until it retires.
- *Milestone 11 starts before Milestone 10's CI run is confirmed green.* The push of Nov 9 starts that run.

**Milestone 11: Indirect predicates.** In `mesher-geom`, module `implicit`. Started and done Nov 9, 2026.

- [x] `implicit::Point`, stored by its defining explicit points and f64 parameters and never rounded:
  - `Explicit`.
  - `Lnc`: t·a + (1−t)·b.
  - `Bpt`: barycentric weights u, v, w on a triangle.
  - `Lpi`: line pq meeting plane rst.
  - `Tpi`: the meeting point of three planes, each through three points.

  The checked constructors reject non-finite input. An exact test also rejects an intersection that does not exist (a line parallel to the plane or lying in it, three planes without a single common point, a degenerate plane) and barycentric weights summing to zero. Each point keeps the exact sign of its `W`.
- [x] Homogeneous coordinates (X, Y, Z, W) for each kind, as polynomials in the defining coordinates. They are computed relative to an explicit point of the query inside the arithmetic, so the filter's error scales with the configuration's size, not its distance from the origin
- [x] Three arithmetic stages behind one trait (`arith::Ring`), with each predicate written once:
  - `Filtered`: f64 with a rigorous running error bound, which also proves exact zeros of exact inputs.
  - `Interval`: outward rounding.
  - `Exact`: dyadic big integers.

  Queries on explicit points only go straight to Shewchuk's predicates.
- [x] Predicates for every combination of point kinds, with the conventions of `predicates`, each with a `_staged` form that also reports the stage that settled it:
  - `orient3d` and `insphere`.
  - `orient2d_projected` and `incircle_projected`: the projection along an axis, the consistent tests for triangulating a facet in its projection.
  - `incircle_coplanar`: the true circumcircle in 3D, via `insphere` with an explicit point exactly off the plane.
- [x] `Point::approx`: the f64 position, for spatial queries and output only
- [x] Stage statistics through the `_staged` predicates
- [x] Tests in `crates/mesher-geom/tests/implicit_oracle.rs` against an exact `BigRational` reference. The reference divides out each point's coordinates, checks them against the definition (an LPI on its line and plane, a TPI on its three planes), and evaluates the textbook formulas; for insphere it also checks the side of the exact circumcentre. The cases:
  - random: five points of random kinds per case, near the origin or 10⁴ away, of size 1 or 10⁻³
  - degenerate: coplanar, cospherical and coincident points in mixed implicit forms, and points on a segment
  - nearly degenerate: intersections of nearly parallel planes and lines (tilts from 10⁻² to 10⁻¹⁵)

  The cases run on all cores in 64 seeded chunks, so they do not depend on the core count. Release builds run 100,000 cases each (about 46 s), debug builds 2,000; `IMPLICIT_CASES` overrides.
- [x] Exit criterion: all 100,000 cases of each test agree with the reference, and the filter settles 100.00% of random queries for all four predicates

**Decisions and findings in Milestone 11.**

- *The exact stage uses big integers, not floating-point expansions.* The roadmap named Shewchuk-style expansions. Insphere on TPI points is a polynomial of degree about 30 in the input coordinates, where f64 expansions can overflow or underflow and lose exactness without notice. Dyadic numbers (`m·2^e` with a `BigInt` mantissa, `num-bigint`) are exact at any degree. Expansions could come back later as a faster stage if profiling shows the exact stage is hot.
- *The filter tracks its error bound at run time, not with the paper's constants.* Diazzi et al. give static error constants per predicate (Appendix A). Rather than transcribe and verify them, the filter carries a rigorous bound through every operation: the textbook first-order bound, plus a term for underflow, inflated to cover the bound's own rounding. It also proves zeros: sums and products of exact values that round without error stay exact. It costs about twice a bare f64 evaluation. Static constants are the optimisation if Milestone 12 profiles it as hot.
- *Intervals rarely help.* Over 190,000 nearly degenerate queries, the interval stage settled 37; the filter and the exact stage settled the rest. Keep it for now (it is cheap), and drop it if profiling says so.
- *Two in-circle tests.* A projection keeps orientation but turns circles into ellipses. `incircle_projected` is the consistent test for triangulating a facet in its projected coordinates, as the current kernel does. `incircle_coplanar` is the true circumcircle, for when the metric matters.

**Milestone 12: Exact constrained Delaunay kernel.** `mesher_tet::cdt`, following Diazzi et al. (2023), beside the current kernel. Started Nov 9, 2026.

- [x] Exact incremental Delaunay core over implicit points (`cdt::ExactDelaunay`): ghost tets, a visibility walk, a cavity of strict in-sphere conflicts grown across every face the point does not strictly see, and BRIO order by approximate positions. It also provides exact duplicate detection (`implicit::coincident`, new), vertex stars, and `check` and `check_delaunay`. Tests in `crates/mesher-tet/tests/cdt.rs`:
  - on random points it matches the float kernel tet for tet (the Delaunay tetrahedralization is unique there)
  - a fully cospherical 5×5×5 grid is valid, Delaunay, and fills the cube
  - LNC, BPT and TPI points insert exactly
- [x] Segment recovery (`cdt::recover_segments`). Every triangle edge is a segment. Each round tries flips on every missing subsegment, then splits what is still missing with an LNC point (concentric power-of-two shells at input vertices, the middle elsewhere) and tries flips on its two halves.
  - Baseline, splitting only (Nov 9, 2026): 77 of the 100-model clean subset recover; the other 23 hit the point budget (20 × inputs + 10,000). Median 3,996 Steiner points, maximum 47,897. The cause was cascades: 43853's first tetrahedralization lacks 684 of its 2,196 segments, and recovery took 36 rounds of splitting near needle triangles.
  - Placing splits by a reference point, as TetGen does, was tried two ways in the baseline (the encroacher nearest an end, and the most encroaching) and made it worse: cascades to f64-short pieces. It was reverted.
  - Flips (Nov 9, 2026). `ExactDelaunay` gained a validated local `replace` (every new tet positive, the outer faces matching the old boundary with the same orientation), 2-3 and 3-2 flips, a hull 2-2 flip, and greedy edge removal (2-3 flips shrink an edge's ring to three tets, then a 3-2 flip). `recover_edge` walks from one end of a missing segment and flips away the face or edge it crosses first, from either end, at most 64 steps.
  - Protected insertion (`insert_protected`, after Si's constrained Bowyer-Watson in TetGen). The cavity starts as the conflict region and shrinks, one tet at a time, until every boundary face is strictly visible, no protected edge and no vertex would be lost, keeping every tet that contains the point. Recovered subsegments are protected, so a split never undoes them. The tetrahedralization is then no longer Delaunay.
  - Result on the clean subset: 92 of 100 recover, at a median of 422 Steiner points (maximum 28,040). On the 76 models the baseline also recovered, the median falls from 3,704 to 322 and no model needs more points. The synthetic spheres now need none; the cube needs 2 (cospherical corners where greedy edge removal finds no valid flip), the wedge 1.
  - Tried and dropped: undoing a failed flip sequence (as TetGen does) recovered 62 of 100, not 92. Flipping the halves right after each split was kept: the same 92 recover, the median falls from 518 to 422.
  - Still failing: 4 models hit the budget, 3 split a piece until it is too short in f64, and 101548 takes over 300 s. Total time is about 1,070 s against the baseline's 525 s, mostly in the conflict region of protected insertion: about 270 exact in-sphere tests per insertion, because the cavities of a non-Delaunay tetrahedralization are large. On 87721, 2,815 of 4,144 splits leave neither half an edge: long segments pass through densely packed vertices, and halving only succeeds once the pieces are as short as the local spacing.
  - Edge removal by one flip (Nov 9, 2026). An edge with up to 20 tets around it is removed by triangulating the polygon of its ring vertices and joining each triangle to both ends (Shewchuk's edge removal, an n-to-2(n−2) flip). A dynamic program keeps the triangulation whose worst tet is best among those whose every tet is exactly positive, so it finds one whenever one exists, where greedy 2-3 flips could dead-end. A hull edge is removed where its two hull faces are coplanar, and the ring's closing edge becomes the hull edge. Where no triangulation is valid, it first removes an edge from either end to a ring vertex, one level deep (after TetGen's n-to-m flips), then tries again. Greedy flips remain for rings over 20.
  - Protected insertion's conflict search now crosses only faces the point strictly sees: the final cavity is star-shaped from the point, so nothing else could stay in it. Same meshes, two to three times faster.
  - Result on the clean subset: 97 of 100 recover, at a median of 130 Steiner points (maximum 8,306). On the 91 models both recover, the median falls from 416 to 118, and one model needs more points. The cube now needs no points. Total time is about 1,280 s, measured with six runs sharing the machine (the previous 1,065 s was not), so roughly the same.
  - Recursion depth: none recovered 96 at a median of 304 on the models both recover (against 416); depth 2 times out on 7 of the first 79 models. Depth 1 is kept.
  - Reference-point placement, tried again with flips and edge removal in: the projection of the vertex most in the way (largest angle subtended), kept a quarter away from the ends, with or without shells at input ends. The median fell 7 to 13 % against the middle, but it was worse on about 30 % of models and lost a model to a too-short piece. Not kept.
  - Still failing: 87721 and 101548 take over 300 s; 462529 splits a piece until it is too short in f64. The too-short failures seen (196190, 95496 and 39245 before recursion) are segments meeting at well under a degree (0.58° at both ends of 196190's segment 24-85), where shells at input vertices do not stop the cascade.
  - Cavity growth in protected insertion (Nov 10, 2026). 462529's cascade was not a sharp angle: a wall of the model lies in a plane that meets segment 320-346 almost exactly at a point, and nearly flat tets in that wall held an edge across the segment there. No flip removes it (three of its ring vertices are coplanar with it), and every split point, a hair off the wall's plane, could not strictly see the flat tets' faces, so the cavity shrank to keep them and the pieces converged on the point. The cavity now first grows across boundary faces the point does not strictly see, as unconstrained insertion does (never across a protected face or into a ghost, and at most to 4 × its size + 64 tets), then shrinks as before. Result on the clean subset: every model's segments recover (462529 with 6,382 points); on the 97 models both runs finish, Steiner points on segments fall from 92,677 to 72,518 in all (fewer on 57 models, more on 5; median 130 to 125), and Steiner points on faces from 316 to 38, in about the same time.
  - A round now checks only the subsegments not yet known to be edges (recovered ones are protected, so they stay); same results; the check had taken 113 s of 101548's profiled run. Tried and dropped: skipping a round's flips on pieces tried just after their split, 30 to 40 % faster but 2 to 20 % more Steiner points.
  - 101548 is the slow one: needle triangles with angles down to 0.006° (7,864 angles under a degree, edges up to 110 long), so segments meeting at such angles need pieces roughly in proportion to 1/angle: 70,156 Steiner points, in 1,500 to 2,000 s alone (two runs). Most of the time is protected insertion (about 7 ms each, mostly shrinking the cavity) and flips (about 3 to 5 ms a try).
  - The example `cdt_segments` measures it on OBJ files (`CDT_MAX_POINTS` caps the vertex count; build with `--features profile` for phase timings).


- [x] Input a valid PLC (the surface gate guarantees it); every input triangle is a constraint, with no facet merging by tolerance. Invalid inputs come back as errors, never panics: on the synthetic suite's open, non-manifold, self-intersecting and overlapping cases the kernel reports a parity conflict, a duplicated vertex, unrecovered triangles or a cavity it cannot fill (mapping these to defect codes is pipeline work, Milestone 13).
- [x] Segment recovery by Si's encroachment-based splitting, with LNC Steiner points (above: splitting with flips and protected insertion first)
- [x] Face recovery (`cdt::recover_faces`, Nov 9, 2026). Each input triangle (a facet) is a polygon: its corners and the segment Steiner points on its sides. Three stages, each only where the last left a facet crossed by an edge:
  - Flips. A dynamic program triangulates the polygon, keeping the most faces the tetrahedralization already has, then the most edges, then the best worst shape; no triangle has three vertices on one side, decided by which sides each vertex is on. Missing diagonals are recovered by `recover_edge`, then missing subfaces by `recover_face`, which removes the edges crossing them. Facets with a missing subface are triangulated again before each pass.
  - Cavity retetrahedralization (`ExactDelaunay::recover_facet_by_cavity`, after Si and the paper). The cavity is every tet whose interior meets the facet. Because the facet's sides are edges, each of those meets the plane only inside the facet, so no boundary face crosses it and the cavity splits into the part above and the part below. One half is filled with the Delaunay tetrahedralization of its vertices, its boundary faces recovered there by flips and its tets found by flood fill. The fill may leave only through faces on the plane, and those faces, which choose how the uncovered part is triangulated, bound the other half; then the other way round. The cavity grows across old boundary faces a fill lacks (TetGen's expansion, at most 4 times; every success on the subset needed 4 or fewer, most none), never across a protected face. The validated `replace` checks the result: every new tet positive and the outer faces matching the old boundary means the two fills tile the cavity. An expansion whose cavity reaches across the plane (the flaw in TetGen's theory the paper's Section 3.5 found, where the halves would overlap) is refused and counted (`ExactDelaunay::expansions_across_plane`), and the cavity also checks that no vertex is lost.
  - Gift-wrapping (Nov 10, 2026, the paper's Section 4.4 after Shewchuk; `cdt/delaunay/wrap.rs`). Where expansion fails, the cavity at each stage, largest first, is filled by gift-wrapping where its halves still lie either side of the plane: each front face takes the apex whose sphere through it holds no other candidate, among those whose tet overlaps nothing (no vertex in or on it, no edge or face meeting the front other than in a shared simplex, exactly), one half leaving the plane open and the other closed by the faces it put there, then the other way round. Cospherical ties are broken by symbolic perturbation of the lifting by vertex index (Edelsbrunner and Mücke), so both halves agree on the plane. No points are added. It can still dead-end, since our segment recovery leaves the mesh non-Delaunay and the paper's guarantee assumes a Delaunay one; on 1120771 most attempts on cavities of 300 to 700 tets do. Orientations are memoized per cavity and candidates ordered lazily: about 0.5 s for a 200-tet cavity.
  - Steiner points on the facet where edges cross it: the exact line-plane intersection (an LPI point) when both ends are explicit, else a BPT point of the input triangle near the crossing. Protected insertion now also keeps every face whose three sides are protected. Then the cavity again, up to 64 rounds. A facet takes one point in its first round and twice as many each round it is still crossed, up to 8 (Nov 10, 2026; it was one a round, and 1120771's facets, crossed by 70 to 100 edges, ran out of rounds), at crossing edges sharing no vertex first, so the points spread rather than follow one fan of edges. A facet with points inside is recovered when no edge crosses it, and its subfaces are read off the tetrahedralization.
  - Every subsegment and the sides of every subface, or of every face already on a facet with points inside, are protected.
  - Result on the clean subset (the 97 models whose segments recover): input triangles crossed fall from 4,177 after segment recovery to 378 after flips, 47 after cavities, and 2 after Steiner points; 95 models have every triangle recovered. 280 Steiner points on faces in all (median 0, maximum 103). The segment Steiner median is unchanged at 130. Total time is about 2,310 s against about 1,280 s for segments alone, mostly in the cavities of five models (77922, 472020, 135674, 64992, 82377, 150 to 250 s each); 87721 now passes 300 s.
  - Tried and dropped: one candidate point per tet around a crossing edge (the middle of its section by the plane) for edges with an implicit end; points went in without removing the crossing, and 1120771 and 64992 ended worse.
  - Tests: Schönhardt's twisted prism (recovered with one segment Steiner point) and a large tetrahedron with small ones just above and below its base, whose crossed base flips recover. Cavities (Nov 10, 2026): a large tetrahedron with twisted prisms hovering just above and below its base reaches the cavity stage (five seeds), and one with spiky star-shaped shells there (seed 127) needs gift-wrapping; a star-shaped sphere (24 × 14, seed 161, after cavity growth made the earlier seed's case go away) has an expansion reaching across its plane, refused, the facet recovered another way (`a_cavity_reaching_across_its_plane_is_refused`). Forcing gift-wrapping on every cavity (a temporary switch), 1,350 random shapes of these kinds meshed with the volume exact.
  - Fixes found by region marking (Nov 9, 2026): an edge lying in a facet's plane across it, with an end off the facet, now counts as crossing it (an exact 2D separating-line test); before, 500113 reported every triangle recovered while one was not. A facet the cavity recovers now has its triangulation read back, so the faces the cavity chose are the ones protected; a later Steiner point had removed one. Every facet is checked at the end, cheaply where its subfaces are all faces.
  - Result with gift-wrapping and the Steiner schedule (Nov 10, 2026, `cdt_region`): 97 of 100 mesh, and every model that finishes recovers every triangle (472078 and 1120771 now do). 21 facets gift-wrapped; 316 Steiner points on faces in all (maximum 110; with gift-wrapping alone it was 74, but 1120771 still failed and 64992 timed out). Total time about 1,920 s with eight runs sharing the machine. Failing: 101548 and 87721 over 300 s (87721 takes 225 s alone, all in segment recovery), and 462529 in segment recovery.
  - Result after cavity growth in protected insertion (Nov 10, 2026): 99 of 100 mesh within 300 s, eight runs sharing the machine, and 101548 does too, in 1,500 to 2,000 s alone; every model recovers every segment and triangle. 16 facets gift-wrapped, 38 Steiner points on faces in all (maximum 11), segment Steiner median 129. Total time about 2,120 s for the 99; slowest 135674 (294 s shared), 472020 (269 s), 1620057 (196 s).
  - Next: the exit criterion run (the synthetic suite through the exact kernel, and the comparison with the paper), and the time of 101548.
- [x] Region marking by the parity flood fill (`cdt::mesh_region`, Nov 9, 2026); deterministic insertion order (BRIO) is in the core. `ExactDelaunay::flood_regions` starts from the ghosts, outside, and toggles at every recovered subface; labels that disagree are an error (`RegionError::Parity`), as are triangles face recovery left (`RegionError::Missing`). The region keeps the exact points (rounding is Milestone 13) and returns the inside tets and the boundary faces, outward, each with its input triangle.
  - Tests: on the synthetic shapes (cube, a hollow cube for nested shells, the spheres, the wedge, Schönhardt's prism, the sandwich with a void), the inside tets' volume equals the enclosed volume, the divergence theorem over the boundary gives it again (so the boundary faces outward), and the boundary is one face per subface with its input triangle.
  - Result on the clean subset (`cdt_region` example): 96 of 100 models mesh, every inside volume within 1.4e-10 of the enclosed volume (from approximate positions); 834,000 tets inside in all. Segment Steiner median 128.5, 226 Steiner points on faces in all (maximum 103). Total time about 1,740 s; slowest 87721 (221 s), 77922 (201 s), 64992 (177 s), 472020 (149 s). Failing: 101548 (over 300 s), 462529 (segment too short to split), 472078 and 1120771 (one triangle each not recovered).
- [x] Exit criterion (met Nov 10, 2026):
  - in exact form, every constraint recovered on the synthetic suite, the 100-model clean subset and the rounded UV sphere. The synthetic suite (`tests/synthetic`, `check_exact`) runs the exact kernel on every repaired part of a closed case that passes the surface gate, at all three scales and offsets: 86 parts, every triangle recovered and the region's volume within 1e-9 of the enclosed volume, 62 Steiner points on segments and none on faces, 0.2 s in all. The subset: all 100 models (above). The rounded UV sphere (16 × 10, rounded to 6 decimals, as `rounded_uv_sphere_meshes_conformally` builds it): recovered with region marking in `crates/mesher-tet/tests/cdt.rs`.
  - a test reproduces a disjoint-cavity case (`a_cavity_reaching_across_its_plane_is_refused`)
  - runtime and Steiner counts against the paper's (4.3 s average over 4,408 Thingi10k models; most runs under 50,000 Steiner points). On the clean subset (eight runs sharing the machine; 101548 alone): median 1.4 s, mean 36 s (21 s without 101548), 72 models under 10 s and 90 under a minute, slowest 101548 at 1,500 to 2,000 s. Steiner points: median 130, 3 models need none, 99 under 10,000, maximum 70,156 (101548, needle angles down to 0.006°), median 0.12 per input triangle. So the counts are in the paper's range but the time is not: the cost is protected insertion and flips on a mesh that is no longer Delaunay (see 101548 above), the first thing to profile if the exact kernel becomes the default.

**Milestone 13: Rounding to f64.**

Started Nov 10, 2026.

- [x] Round every implicit point; check the exact `orient3d` of every tet on the rounded coordinates (`cdt::round_region`, `cdt/rounding.rs`). Each constructed point becomes the f64 point its construction approximates; input points are f64 already. A tet counts as usable when it is positive exactly and has positive volume in floating point (as a finite-volume solver computes it, and as the volume gate's `flat_in_float` checks), from the best of its 12 orientation-preserving vertex orders, which is the order written: a tet nearly flat to rounding can come out zero from one corner and positive from another (the synthetic torus's quads, coplanar in exact arithmetic but not in f64, made the only flat tets; reordering fixed every one).
- [x] Untangle what flipped (`cdt/delaunay/untangle.rs`): a constructed vertex of the bad tet moves to the nearest f64 point within 2 ulps per coordinate of its rounding that leaves no inside tet around it unusable and no more outside ones; else a 2-3 flip across a face that is not a boundary face, or edge removal of an edge that is not protected (subsegments and subface sides); while untangling, flips must give tets usable in floating point too. Up to 8 passes. Each fix is returned. Extra explicit points were not needed on the subset or the suite, so they are not implemented; a mesh still tangled is `MESH_ROUNDING_FAILED`.
- [x] Re-check boundary conformity after rounding: every inside tet usable, the boundary faces exactly the recovered subfaces, each once, and each input triangle's area covered by its subfaces to 1e-9 of its area plus 1e-12 of the squared bounding-box diagonal (the boundary moves by rounding only).
- [x] Codes `MESH_ROUNDING_REPAIRED` (warning, one record per fix: `fix` is `vertex_move` with `moved`, `flip_2_3` or `edge_removal`, at `location.point`; action `untangle`, a new `Action`) and `MESH_ROUNDING_FAILED` (error with `tets_left`, `turned_over_before` and `fixes`, which triggers the next fallback). Report schema, `docs/report.md` and the taxonomy table updated.
- [x] `--kernel exact` in the pipeline (`JobConfig::kernel`, `KernelChoice::{Conformal, Exact}`, `Kernel::Exact` in attempts and summaries): for a part that passed the surface gate, the exact kernel and rounding run first, without refinement (the tets join the surface's vertices and the Steiner points recovery needed; refinement is Milestone 14), judged by the volume gate with the conforming kernel's dihedral bounds; then the conforming kernel, the remesh retry and stuffing as before. An exact failure (`MESH_RECOVERY_FAILED` or `MESH_BUDGET_EXCEEDED` with `kernel: exact`, or `MESH_ROUNDING_FAILED`) becomes a warning with `fallback` / `applied` when a later kernel meshes the part, an error otherwise. Boundary faces carry their input triangle, so OBJ groups become side sets as for the conforming kernel. Test: `exact_kernel_meshes_a_cube_and_falls_back_where_it_cannot` (worker protocol). The GUI has an "Exact kernel" switch in the Run panel (Nov 10, 2026), saved in project files as `kernel = "exact"`.
- [x] Exit criterion (met Nov 10, 2026):
  - rounding validity before and after untangling reported on the clean subset and the synthetic suite. Subset (`cdt_round` example): of the 99 models whose exact mesh finishes within 300 s (101548, 1,500 to 2,000 s, not measured; 135674 alone), 91 round validly as they are (92 %; the review's figure is 93.2 %) and all 99 after untangling. The other 8 had 47 inside tets turned over (36 in 87721), fixed by 24 vertex moves and one edge removal; none was flat in floating point once tets are written in their best vertex order (before that, 27 were, in 4 models). Rounding takes under 0.3 s a model. Synthetic suite (`check_exact`): all 86 parts valid as rounded, the volume still within 1e-9.
  - no non-positive tet written: `round_region` returns a mesh only when every tet is usable, and the volume gate checks again
  - every untangle a record (`MESH_ROUNDING_REPAIRED`)

**Milestone 14: Chamfer-and-refine.** Started Nov 10, 2026, from Diazzi, Dai, Panozzo and Attene, "Surface chamfering for robust tetrahedral meshing" (ACM TOG 45(4), SIGGRAPH 2026), clean-room from the paper's text (it is CC BY-NC-ND and its code is not read). The method, as the paper describes it: filter out input triangles responsible for a local feature size below τ (1e-4 of the diagonal; they come back at the end); chamfer every acute vertex and edge, so the surface meets the projection condition; Delaunay-refine the chamfered surface in an outer box (split encroached segments at LNC midpoints, then missing or encroached Delaunay triangles of each face at their circumcentres as BPT points unless they encroach a segment, then tets whose circumradius over shortest edge exceeds 2, about 14.47°, at f64 circumcentres unless they encroach a subsegment or subface; optionally slivers by shortest edge over shortest opposite-edge distance, 7.2); then the enriched CDT: the refined vertices outside the chamfered regions with the input, through the exact constrained kernel, for exact conformity with bad tets only near the former acute angles. The paper reports 3.6° average minimum face angle (1.3° for the CDT alone), about 120 s a model and one 36-hour outlier.

Tasks, in order:

- [x] Exact angle predicates on implicit points (Nov 10, 2026): `implicit::dot_sign` for the sign of (p − q)·(r − q) (is the angle at q acute?), `implicit::dihedral_sign` for [(b − d)×(c − d)]·[(a − d)×(c − d)] (do faces a c d and b c d meet at an acute dihedral?), and `implicit::edge_projects_into_triangle` (both dihedrals at the shared vertex's sides acute), in homogeneous form like `orient3d`, filtered and staged. Oracle (`angle_queries_agree`, `edge_projection_into_a_triangle`): 100,000 random queries each agree with the rational reference, the filter settles all of them, and about 50,000 each of exactly right angles and right dihedrals, built on an integer lattice in implicit forms, come out zero
- [x] Acute classification of the input surface (`mesher_tet::chamfer::classify`, Nov 10, 2026): acute edges (dihedral, either side, since refinement meshes the whole box), flat edges (coplanar, folding out to 180°), acute vertices (an acute incident edge, two incident non-flat edges at an acute angle, or a non-flat edge whose projection meets an incident face); every decision exact. Local feature size per vertex (its smallest triangle height there) and per edge (distance to the nearer opposite vertex), approximate: it only chooses ε. Tests: a cube has nothing acute, a regular tetrahedron everything, a prism exactly its three lengthwise edges, a leaning pyramid its apex and far corner. On the clean subset (`acute` example, 0.4 s for all 100): a median of 56 % of vertices are acute (curved surfaces bend a little at every edge, so almost every vertex there is; 12 models have none) and 2 % of edges; 54 models have a local feature size under 1e-4 of the diagonal, the smallest 8e-19, so input filtering is needed as in the paper
- [x] Orthogonal chamfering (`mesher_tet::chamfer::chamfer`, Nov 10, 2026). At each acute vertex P, every incident edge is cut at ε = local lfs(P) / 3 (an LNC point, shared by the edge's two triangles), and in each triangle at P the corner is replaced by a bridge of three equal edges whose legs are perpendicular to the triangle's edges (the leg length is the smaller root of a quadratic); the bridge points are BPT points of the triangle (f64 position, barycentric weights kept inside it). Perpendiculars come from cross products, which stay accurate for needles. Then in each triangle at an acute edge the legs are shortened to the triangle's shortest and the cuts on the acute edge dropped, removing the rectangle between them. Each input triangle leaves one polygon (in its orientation); flat edges are cut like any other but are not segments. Tests: a cube is untouched, a regular tetrahedron's faces become hexagons of bridge points, a wedge prism; every face vertex lies on its triangle, every face keeps positive area under its triangle's, and no new vertex has an acute angle between its face edges. On the clean subset (`chamfer` example, 1.1 s for all 100): every model chamfers; in 98, no two segments (face edges off flat input edges) meet under 89.9° (the rest is rounding the bridge points); in 360073 and 39245, local feature sizes near 1e-14 shrink the legs to a few ulps, so their rounded positions show 0° angles (input filtering removes features that small). The chamfered surface has a median of 11 times the input's vertices (19 at most): curved surfaces are acute at almost every vertex.
- [x] Bridge simplification (`chamfer::simplify_bridges`, the paper's Section 4.4, Nov 10, 2026): bridges in order of their shortest edge, each replaced by one edge between its cuts, else by two edges dropping one bridge point (the longer new edge first), wherever the new edge makes no acute angle with the segments at its ends and leans over no other face there (exact `dot_sign` and `edge_projects_into_triangle`). `Chamfered::is_segment` tells segments from face edges along flat input edges. Test: a UV sphere, every vertex acute, loses bridge points and keeps every face angle non-acute. On the clean subset it drops a median of 42 % of the chamfered surface's vertices, from 10.6 to 5.5 times the input's (18 to 10.3 at most), and no two segments meet under 89.98° except in the two models with features near 1e-14; 2.8 s for all 100
- [ ] Input filtering by τ, alternated with chamfering until no small feature remains
- [ ] Conforming Delaunay refinement of the chamfered surface in an outer box (on `cdt::ExactDelaunay`, which keeps the Delaunay property here: no protected insertion), with a priority queue by radius-edge ratio and the work budget
- [ ] Marking the tets inside the input, then the enriched CDT through `cdt::recover_faces` (which then has to take implicit input vertices, as the paper's Appendix C does) and `cdt::round_region`
- [ ] Refinement against the sizing field. Opt-in and under the work budget; it hands over to Milestone 15 on `MESH_BUDGET_EXCEEDED`
- [ ] Exit criterion, on the clean subset:
  - refinement ends within budget on every model, or records why
  - minimum face angle reported against the paper's (1.3° to 3.6°)
  - cost per model tabulated

**Milestone 15: Worst-cell quality for finite volume.**

- [ ] Invariant-guarded local operations (flips, splits, collapses, vertex moves). Each is accepted only if every tet stays positive exactly, the boundary stays exact and the target improves. Sliver smoothing becomes one of them
- [ ] A WSVM-style energy as the vertex-move objective, measured against the candidate search
- [ ] Targets: per-face non-orthogonality and skewness, cells with a boundary face first
- [ ] Report the worst cells and their locations, wall-adjacent cells separately. `VOL_LOW_QUALITY` is gated on non-orthogonality and skewness too
- [ ] Absorb improvements 5, 6 and 8; parallel by independent sets with a fixed colouring
- [ ] Exit criterion:
  - on the clean subset, worst wall-adjacent non-orthogonality and skewness below the current kernel's and TetGen's
  - worst cells located
  - byte-identical results at 1, 4 and 16 threads

**Milestone 16: Exact arrangement (repair rung 2).**

- [ ] Arrangement of each part's triangle soup with LPI and TPI points (Cherchi et al. 2022); Guo and Fu's indirect offset predicates (2024) evaluated as the faster variant
- [ ] Cells labelled by generalized winding number, keeping the solid region. This resolves `INT_SELF_INTERSECTION`, `INT_SHELL_OVERLAP`, `PART_OVERLAP` and touching shells (improvement 9)
- [ ] Winding numbers by the Antipodal method's ray-crossing form, where faster
- [ ] The resolved surface goes to the exact kernel, or is snap-rounded and re-checked
- [ ] Report which rung succeeded, and the components and Euler characteristic before and after; proposed `TOP_AMBIGUOUS_INSIDE`
- [ ] Exit criterion:
  - every part of the self-intersecting and multi-component strata resolves to a surface that passes the gate, or fails with a located record
  - component counts and Euler characteristic match a brute-force check on the synthetic intersection cases

**Milestone 17: Wrapping (repair rung 3).**

- [ ] Alpha wrapping after Portaneri et al. (2022), clean-room. Alpha and offset come from the element size and the diagonal; proposed `REPAIR_WRAPPED`
- [ ] The wrapped surface goes through the exact kernel; stuffing stays last
- [ ] Exit criterion (and the phase's): the full benchmark on all 10,000 models and the valid subset; outputs checked by `mesher-fvm` and by the user's solver

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

## Things to consider improving

Collected at the end of Milestone 8 (Nov 7, 2026), roughly in order of value. Each is a candidate, not a commitment; the Thingi10k benchmark (`tests/oracle/tests/thingi10k.rs`) measures them.

1. *Hangs with no mesh.* Done in Milestone 10: a deterministic work budget, a cascade handover to full recovery, a stuffing retry after a worker timeout, and the profiling speed-ups; 91607 and 677341 now end in a recorded fallback.
2. *A per-model fallback to conforming recovery.* The earlier conforming kernel meshed 13 models that the constrained one does not. When constrained recovery fails, retrying the old way (subfacets unconstrained, Gabriel splits) before stuffing would recover them, at the cost of keeping both paths.
3. *Facets coplanar only to float32 noise.* Neighbouring facets at a dihedral of about 1e-7 leave thin wedge tets whose faces overlap a facet's region, so adopting the mesh's facet triangulation fails and recovery stalls (most of the 10 regressions). Likely fixes: choose a non-overlapping subset of the candidate faces, or flip the wedges away first.
4. *Near-degenerate input triangles.* Nearly collinear triangles (reported as `QUAL_SLIVER_FACE`) produce coincident subfacets of opposite orientation on one edge. Detect and handle them before meshing, or merge them with a recorded repair.
5. *Early stop.* Recovery and refinement should give up when progress relative to the input stalls, rather than running to the budget or the timeout; some models take minutes before falling back.
6. *Point-budget overruns* (57657, 73154; 43853 is fixed by the sizing field). Shape refinement near small input features still adds too many points, and the size field does not bound it (see the Milestone 10 findings). Candidates: local feature size from the input's own small features (short edges, small angles) to stop shape refinement there, or protecting small input angles as TetGen does.
7. *Speed.* About 20× TetGen's total time before Milestone 10's profiling; boundary insertions were the hot spot (a cavity shrunk one tet per pass, about 150 passes per insertion on some models). Those passes are now cheap, but their number is not: shrinking several tets at once, or smaller cavities for boundary splits, would change results and needs a benchmark run to judge. Winding-number inside tests are the next largest cost (about a third on 43853); inheriting inside-ness across new tets would remove most of them.
8. *Boundary slivers.* Smoothing of free vertices (Milestone 10) removed most cells below 1° that have a vertex off the boundary. What remains is of two kinds. Tets on input triangles that are themselves needles (113890 has 426 input triangles with an angle under 1°): nothing that keeps the boundary can fix them, but remeshing can. And tets whose four vertices are all on the boundary (thin regions, and most of the remeshed 472020's 0.4°): sliding boundary Steiner points within their facet, or 4–4 flips and multi-face removal, would be next.
9. *Touching shells conformally.* Shells sharing an edge or a vertex fall back to stuffing; meshing them as one constrained region would keep them exact.
10. *The SEACAS check* of the Exodus output (deferred from Milestone 5; do it before Milestone 10).
11. *Tets flat in floating point.* Mostly done (Nov 9, 2026; see the Milestone 10 findings): of the 16 benchmark models whose conformal meshes had them, 13 now mesh conformally with none. Two still have them (50659: 57, 45136: 14). Recovery splits an input edge near a corner into Steiner points about 1e-7 apart and a few ulps off the line, so the facet's triangulation needs needle triangles through three nearly collinear points, and any tet standing on one is flat. Collapsing such a Steiner point into its neighbour along the edge, or not splitting an edge that finely, would remove them; until then those parts use the fallback.
