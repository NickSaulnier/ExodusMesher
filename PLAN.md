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
| `INT_SELF_INTERSECTION` | Faces of one shell intersect | Resolve (repair rung 2, Milestone 16) for a closed part, warning; else error, triggers fallback |
| `INT_SHELL_OVERLAP` | Two shells of one part intersect | Resolve (repair rung 2, Milestone 16) for a closed part, warning; else error, triggers fallback |
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
2. **Remesh.** If the conformal mesh fails, or its worst cell is beyond the volume gate's quality limits (smallest dihedral angle under 10°, non-orthogonality over 70°, or skewness over 2: `worst_cell`), remesh the whole part's surface to the size field (keeping sharp edges, corners and group boundaries, but merging features under a tenth of the size: `RemeshParams::small_feature`; and along sharp edges choosing the walls' diagonals so the cell on both walls is not twisted: `align_sharp_edges`) and run the Delaunay mesher on it once more. The mesh then conforms to the remeshed surface, and the volume gate measures its deviation from the repaired one. Of two meshes, the one with the better worst cell goes to the volume gate first (the remeshed one only if strictly better), the other if it is rejected; a remeshed mesh kept for quality gets a `VOL_LOW_QUALITY` info record (`remesh/applied`), and findings deferred to remeshing (`QUAL_SLIVER_FACE`) become applied. A clean input keeps its own triangulation and so its exact volume (a remeshed curved surface follows the input only to the element size: the synthetic torus loses 2 %). Bounded to one remesh so it cannot loop; `--no-remesh` skips it (Milestone 15, Nov 12, 2026; it was a retry on failure only).
3. **Stuffing fallback.** Run isosurface stuffing on the winding-number field. The result is valid with bounded dihedral angles but approximates the boundary, so measure and report the Hausdorff distance to the original surface. Stuffing is also tried when every conformal mesh is far beyond the quality limits (the best one's `worst_cell` under `FAR_WORST_CELL`, 0.25: a dihedral angle under 2.5°, non-orthogonality over 85° or skewness over 8), and kept if it passes the volume gate with a better worst cell and holds the surface's volume to 1 % (`FAR_MAX_VOLUME_ERROR`: thin features it misses would be lost); it is gated into a scratch log first, so its records count only if it is kept. Kept, it gets a `VOL_LOW_QUALITY` warning (`fallback/applied`) with both worst cells, and the conformal attempts are listed as failed; `--strict` never does this (Milestone 15, Nov 12, 2026).
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
- [x] Record capping: at most 100 detailed records per code and stage (per code until Nov 13, 2026, when diagnosis filling a code's share dropped the rung-2 and surface-gate records of it); the summary keeps exact totals per code, so garbage input cannot produce an unbounded report
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
- [x] Dihedral angles (`mesher_tet::quality`); property test on random spheres, ellipsoids, boxes and tori at random cell sizes: angles within [9.3171°, 161.6432°] (corrected Nov 12, 2026; see the decisions below), no degenerate tets, every face shared by at most two tets, closed and consistently oriented boundary (1,000 shapes run once)
- [x] Volume stage (`mesher-core::volume`): stuffing on each repaired part's winding-number field, volume gate (`VOL_INVERTED_TET`, `VOL_LOW_QUALITY`, `VOL_BOUNDARY_DEVIATION`), `MESH_BUDGET_EXCEEDED`, attempts recorded; accepted parts assembled into one mesh, a block and a `<part>_boundary` side set each
- [x] Worker writes `mesh.exo`, `mesh.msh`, `mesh.vtu` (`EXO_NAME_TRUNCATED`, `EXO_CAPACITY_EXCEEDED`, `OUT_WRITE_FAILED` recorded); the report's part summaries gain mesh statistics; the CLI prints them
- [x] Worker memory cap: the CLI runs the worker under `process_control` (a job object on Windows, `setrlimit` on Linux), `--memory-limit` MiB (default 8192); allocation failure becomes `MESH_WORKER_FAILED` with reason `out_of_memory` (exit 30, not a bug); tested with a deliberate runaway allocation
- [x] `--element-size` (relative to each part's diagonal, default 0.03); `--strict` records `MESH_STRICT_UNMET` until the conformal kernel exists
- [x] Octree grading (adaptive lattice): done in Milestone 10 (graded interior, uniform surface, as the paper proves)

**Decisions made in Milestone 6.**

- *All-zero tets are dropped.* A lattice tet whose four vertices were all warped onto the surface is a sliver lying in it; keeping it (as an early version did when its centroid was inside) produced dihedral angles down to 2.8°. Dropping it leaves its faces, which lie on the surface, as boundary. With that, the empirical range on all tests stays inside the paper's bounds.
- *Stencil filling by coning.* Each lattice tet's inside part is filled by a cone from its lowest-numbered inside vertex over its boundary triangles; a quadrilateral (only on faces with two inside vertices and one outside) is split from the lower-numbered inside vertex. Both rules depend only on the shared face, so neighbours conform. This replaces the paper's colour-based parity rule; the property tests above are the evidence that the angle bounds still hold.
- *The proven bounds are [9.3171°, 161.6432°], not [10.7°, 164.8°] (corrected Nov 12, 2026).* The warping parameters (α 0.24999 long, 0.41189 short) are the paper's Table 1 row "minimum dihedral angle, safe", whose computer-verified bounds are 9.3171° and 161.6432°. The 10.7° and 164.8° quoted in the abstract belong to the "unsafe" row (α 0.28511 and 0.39882), under which warping can invert a BCC tetrahedron. The code had paired the safe parameters with the unsafe row's bounds, and a property-test case found the gap: an ellipsoid (axes 0.675, 1.530, 1.311 at h = 0.1929) has a whole BCC tetrahedron (two red and one black vertex warped about 0.40 of a short edge onto the surface, one inside black vertex) with a 10.23° dihedral angle. No stencil is involved; each warp is within α, so the mesh is correct and the bound was wrong. We keep the safe parameters, because an inverted BCC tetrahedron would be re-oriented by the stencil step and overlap its neighbours, and use the safe row's bounds in `stuffing::DIHEDRAL_MIN_DEGREES` and `DIHEDRAL_MAX_DEGREES`, which the volume gate and the synthetic suite read. The tighter maximum also holds with our coning stencils (property tests raised to 2,000 uniform and 1,000 graded cases, run three times: 9,015 meshes, extremes 10.18° and 159.17°; 10.18° is a second shape below the old 10.7°).
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
- [x] Graded stuffing lattice (`mesher_tet::graded`, `stuff_graded`), following Labelle and Shewchuk's Section 5 so the proven bounds hold: leaf cubes found by continuation from surface seeds, the Continuation Condition and the half-pyramid guard, an octree with partial children up to element-size roots, Weak Balance, and the Figure 10 tetrahedra (BCC, bisected, quadrisected, half-pyramids). The cutting, warping and stencil code is shared with the uniform lattice through a `Background` trait. In the pipeline, a part whose size field asks for less than h is stuffed with surface cells of h/2 or h/4 grading to h inside; uniform fields keep the plain lattice. Property test: 500 random graded shapes all within the proven bounds, conforming and closed; a graded sphere has the uniform one's volume to 1e-9 with less than half the tets. 472131's fallback: deviation 4.24% → 0.61%, 2,012 → 82,958 tets
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
- [x] Input filtering by τ (`chamfer::small_features`, `chamfer::chamfer_filtered`, Nov 10, 2026), τ = 1e-4 of the bounding-box diagonal as in the paper. First pass: a triangle with an edge or a height under τ goes, and so do both triangles of any pair with non-incident features closer than τ (a vertex and a triangle not containing it, or two edges sharing no vertex; candidate pairs by sweeping τ-grown boxes along x; f64 distances, as they only choose what is filtered). Then chamfer (classification now takes open surfaces: an edge with one triangle is a segment, never acute or flat) and simplify bridges, and while a chamfered segment is shorter than τ/3, filter out its triangle and chamfer again. A third, not τ: after the first pass every local feature size is at least τ, so ε = lfs/3 cuts as little as τ/3 by construction, and holding chamfered segments to τ filtered out every triangle a chamfer touched, cascading (472020: 347 triangles after the first pass, 8,170 in the end); the paper does not give its threshold here. Test: two cubes 1e-3 apart lose exactly the triangles near the gap at a tolerance of 1e-2 and none at 1e-4. On the clean subset (`chamfer` example, 11.7 s for all 100): 73 models filter something, a median of 6.7 % of triangles (15.7 % with τ for the second check); a few CAD-like models lose most of theirs (101548 84 % over 30 rounds, 135674 80 %, 1620057 64 %: a quarter to two fifths of their triangles already fail the first pass), which come back with CDT quality in the enriched CDT. Every chamfered surface then has no two segments meeting under 89.9°.
- [x] Conforming Delaunay refinement of the chamfered surface in an outer box (on `cdt::ExactDelaunay`, which keeps the Delaunay property here: no protected insertion), with a priority queue by radius-edge ratio and the work budget
  - [x] Segments (`mesher_tet::refine::recover_segments_conforming`, Nov 10, 2026). The chamfered surface's vertices in use and the corners of a box (its bounding box grown by 20 % of its diagonal) are tetrahedralized; every segment (face edges off flat input edges) and box edge is split at its middle while it is not an edge or is encroached (a vertex of the tets around it strictly inside its diametral sphere, exactly: `dot_sign` negative). Split points: LNC points of the input or box edge a segment lies on; BPT points of the input triangle, at the average of the ends' barycentric weights, for bridge segments (exactly on the triangle, not necessarily on the segment). After an insertion only the subsegments at the new vertex's neighbours that it lies inside or on the diametral sphere of are checked again, and a final full pass catches the rest (83 s to 10 s on 55277, the same points). `ExactDelaunay::star` switches from a linear scan to a set past 32 tets (box corners have thousands); the exact kernel's meshes are unchanged. Test: a cube's and a fully chamfered tetrahedron's segments come out unencroached Delaunay edges, the tetrahedralization valid and Delaunay. On the clean subset (`refine` example): all 100 recover, a median of 2 Steiner points per chamfered vertex (472175 needs 41,020 for 403: features near τ force pieces that short along long segments), median 7 s, slowest 102 s with eight runs sharing the machine. Most of the time is the conflict search of exact insertion (about 330 µs a point). The shortest subsegment is about τ/3 of the diagonal, except 43853 at 2.8e-7: the second filtering step checks only chamfered segment lengths, not distances between non-incident chamfered features.
  - [x] Faces (`refine::recover_faces_conforming`, Nov 10, 2026). The module became `refine/` with one engine for both stages: segments are chains of pieces (carrier, vertices, parameters, the regions they bound), and each planar region (chamfered faces merged across flat edges by union-find, and the six box sides) keeps a 2D Delaunay triangulation of its vertices (`refine/plane.rs`: exact, implicit points on one plane, projected orientation and the true in-plane circle test `incircle_coplanar`; ghost triangles and in-place cavity filling as in 3D; tests on random BPT points on a tilted plane, a fully cocircular grid and a square with midpoints then centre). The triangles inside a region come from a flood from the inner side of its oriented boundary pieces; each missing from the tetrahedralization, or with an apex of its two tets inside its equatorial sphere (f64, an apex on the sphere to 1e-9 not counted: it only chooses where to refine), is split at its circumcentre (a BPT point of the region's plane), unless that point encroaches a boundary piece, which is split instead. Only regions a new point's neighbours lie on are checked again, then a full pass. A bad triangle whose fix adds no point and that is still there is set aside and counted (`Refined::stuck_subfaces`); without that the loop could spin forever (first run: 19 of 100 timed out, several after 3 s of segments). Test: a cube's and a chamfered tetrahedron's faces and segments recovered with nothing set aside, the subfaces off the box covering the chamfered faces' area to 1e-9, the tetrahedralization valid and Delaunay. On the clean subset (`refine` example): 97 of 100 finish within 300 s with eight runs sharing the machine; a median of 194 Steiner points on faces (49,410 at most) against 6,797 on segments; 29 subfaces set aside in 16 models (4 at most); median 12 s, slowest 190 s (500113). 95496, 104738 and 225958 take longer (95496 is still adding points after 20 minutes): each pass rescans a whole region (one of 11,000 triangles there is marked again by every nearby insertion), and the box sides need refining too (one with 897 of 922 triangles bad). Then, with the tet stage (the same day): a triangle no longer in its plane's triangulation when its turn comes is skipped (the triangles of a cocircular quad share a circumcentre, and fixing the second put a point an ulp from the first: 49423 ran away to 390,000 points), split points on axis-aligned lines and planes are explicit, and passes check only the triangles near a change (see Tets). On the clean subset: 98 of 100 finish within 300 s, median 2.1 s, slowest 39 s apart from 95496 (169 s); a median of 76 Steiner points on faces (61,383 at most); nothing set aside. 104738 and 225958 still time out.
  - [x] Tets (`refine::refine_conforming`, Nov 10, 2026). Every tet whose circumradius over shortest edge exceeds the bound (2) is queued, and they are split worst first (a heap by ratio, then vertices) at their f64 circumcentre, unless it lies outside the box or encroaches a subsegment or subface; those are split instead. Segments and faces are repaired before each tet. The encroachment candidates are the edges and faces of the tets the point conflicts with, checked on the conflict region the insertion itself finds (`ExactDelaunay::insert_unless`), so it is searched once. Each insertion queues the bad tets it makes. A bad tet whose fix adds no point is set aside and counted (`Refined::stuck_tets`); one whose fix split something else goes back in once the queue is empty. A full check of segments and faces ends it, and it starts again if that added points.
    - Speed-ups, each leaving results identical (checked on 40 models, segments, faces and tets):
      - The implicit determinants expand each minor once (Laplace expansion over column subsets: 75 products and one table for a 5×5, against 120 products and 85 allocations). This takes the exact kernel's conflict search from 436 to 40 µs a point, and the constrained kernel shares it.
      - A region pass checks only the triangles at the vertices touched since the last pass, plus any a fix left in place, in slot order. An insertion changes only faces whose vertices all become its neighbours, and a plane's new triangles all contain the new point.
      - Inside labels are kept per triangle slot: a new triangle's label is set from the edge opposite the new point (a boundary piece's side, else its neighbour's label, else outside at the hull), instead of a flood each pass.
      - The subface fix's scan of the region boundary has an f64 distance pre-filter.
      - The plane's duplicate scan runs only before its first triangle.
      - Survivors held back until the queue empties replace the full rescans between rounds (55277: 15 rounds to 1).
    - Timings: 40179 went from 9.9 s to 1.6 s, and 55277 from 139 s to 49 s.
    - Test: a cube and a chamfered tetrahedron end with no tet over the bound bar any set aside, every subsegment and subface present, and the tetrahedralization valid and Delaunay (0.7 s; the cube's faces leave nothing to split).
    - On the clean subset (`refine` example, eight runs sharing the machine, 600 s each):
      - 98 of 100 finish, median 14 s, slowest 312 s (98939, 5.5 million tets).
      - A median of 352,000 tets and 29,000 Steiner points in tets (476,209 at most), and 8,400 more on faces, since tet circumcentres encroach subfaces.
      - 3,003 tets set aside in 65 models (284 at most); no subface set aside.
      - 104738 and 225958 still time out in the face stage.
    - Locating plane insertions from a nearby vertex would halve the walk but changes slot order, and so results, and gained little overall; not done.
    - Slivers are not removed yet (optional in the paper).
- [x] Marking the tets inside the input, then the enriched CDT through the exact kernel and `cdt::round_region` (`refine::enriched_cdt`, Nov 11, 2026; `--example enriched` runs the whole chain).
  - **Constraints:**
    - The constraints are the input triangles themselves, not the paper's refined faces plus cut-away pieces.
    - Every refined point on the chamfered surface lies exactly on an input edge or triangle:
      - cuts and edge splits are LNC points of the input edge;
      - bridge and face points are BPT points of an input triangle;
      - snapped points lie exactly on their axis-aligned plane.
    - `cdt::recover_faces_with` takes them as extra points. `cdt::locate` places each one exactly on the segment or triangle it lies on: a grid of triangle boxes and f64 plane distances pick candidates, then `orient3d` and projected orientations decide, and `dot_sign` orders points along a segment. A copy of an input vertex is left out.
    - So every point the kernel adds is built from the input's explicit vertices, and no implicit point is built on another: the paper's Appendix C (cascaded implicit points) is not needed.
  - **Merged facets:**
    - Coplanar triangles across flat edges are one facet (`merge_flat`), as refinement's planar regions are, and flat edges are not segments.
    - Without this, every flat diagonal had to be recovered with Steiner points: the cube's smallest face angle was 6°, against 23° with merging, and 40179's was 0.36°, against 5.4°.
    - A merged facet's crossings, cavity and Steiner points work on its triangles (`ExactDelaunay::crosses_region`, `recover_region_by_cavity`). An edge through a flat edge's inside counts as crossing.
    - Its subfaces come from a flood from the inner side of its boundary pieces, since faces with vertices on it can span a notch outside it.
    - It is covered exactly when that flood reaches every boundary piece once and pairs every other edge. This exact, combinatorial check replaced listing crossings: on 40179 the kernel went from 11.1 s to 0.74 s.
    - Rounding checks coverage by facet (`RecoveredFaces::facet_of`).
    - The existing exact-kernel path is unchanged; the synthetic golden reports match.
  - **Refined points left out:**
    - Points off the surface are left out where refinement could not see the input:
      - in a chamfered vertex's ball, of radius its farthest bridge point (`Chamfered::radius`);
      - in an acute edge's strip (`Chamfered::strips`);
      - within τ of a filtered triangle.
    - So are the split points inside bridge segments (`Refined::on_bridges`). They lie exactly on their triangle but only near the bridge (BPT points at averaged weights), so a chain of them zig-zags at the rounding level. With the bridge no constraint here, the kernel's facet takes faces across three of them, at 0° face angles.
    - Leaving them out took 128915, 74494 and 73154 from 61 to 136 tets under 1° to none, at about 1° minimum, and rounding there from up to 115 tets turned over to at most 1.
  - **Untangling fix:** while untangling, `replace` now refuses a new inner face or edge that already exists outside the tets it replaces. Rounded positions can be tangled, so positive new tets need not tile the old region alone; a duplicated face had broken the inside labels in 3 models.
  - **Test:** `enriched_cdt_of_a_cube_and_a_chamfered_tetrahedron`. The cube keeps radius-edge 2, and both round to their enclosed volume. Also `extra_points_and_merged_flat_facets` in `tests/cdt.rs`.
  - **On the clean subset** (`enriched` example, eight runs sharing the machine, 900 s each): 82 of 100 mesh, conform exactly and round.
    - Failures and timeouts:
      - 16 time out: 104738 and 225958 in refinement; most others in the kernel, whose median is 38 s against refinement's 10 s.
      - 500113 leaves one input triangle unrecovered.
      - 1501546 keeps one tet turned over after rounding.
    - Kernel Steiner points: segment Steiner points in 70 models (median 186, at most 16,306), for input edges in the cleared balls and strips and filtered features.
    - Rounding: 18 models had tets to untangle (27 at most).
    - Face angles: only 18 of the 82 inputs have no angle under 1° (median smallest input angle 0.06°), and a conforming mesh cannot beat its input's needles.
      - For those 18, the smallest face angle has a median of 1.0°, against the paper's 1.3° to 3.6°.
      - 591216 still has a tet at 0°.
      - 33 of the 82 meshes have no face angle below their input's smallest.
    - Median 87,500 tets (1.1 million at most).
  - **Kernel speed** (Nov 11, 2026). Profiled on 332810, 42936, 112922 and 472175, whose kernels took 140 to 855 s. The first five changes leave every mesh identical (checked on those four models; the synthetic golden reports match):
    - `crosses_region` tests the plane once per edge, not once per triangle. On 332810 the kernel went from 140 s to 14 s.
    - The protected edges are kept up to date facet by facet (`Protected` in `cdt::faces`), instead of being rebuilt in full before every cavity and Steiner round. `ExactDelaunay::log_touched` logs the vertices of every tet made, and only facets at those vertices, or whose own subfaces or points changed, are read again. Debug builds check the result against a full rebuild. A merged facet's covered flag comes from the same read.
    - `Point::bounds`: a box certainly holding an implicit point, by interval arithmetic rounded outwards (checked against exact rationals in `implicit_oracle.rs`, nearly parallel line-plane intersections included). Gift-wrapping passes over vertices and front faces whose boxes miss the candidate tet's, and the cavity collection over triangles whose boxes miss the edge's. Each vertex's side of the facet plane is remembered within a cavity.
    - A triangle's polygon is triangulated into subfaces only where they are used: the dynamic program is cubic in the polygon's length, and its result was thrown away for every merged triangle (on 472175, 478 s of the kernel's 532 s).
    - With `merge_flat`, a lone triangle whose polygon is longer than 512 vertices is a planar-region facet too. This changes results. Making every triangle a region was faster still (median kernel 7.7 s, against 10.9 s) but lost 472175, 1536355 and 1224385; at 512, 472175's 231-vertex triangle is recovered by flips, as before. Turning a triangle that flips leave crossed into a region afterwards was tried and recovered nothing more: its Steiner points then fall on protected faces.
    - **On the clean subset** (same settings): 91 of 100 mesh, conform and round, against 82.
      - On the 82 meshed in both runs, the kernel's median went from 37.6 s to 8.5 s, and its total from 11,042 s to 2,275 s. 76 meshes are the same; 6 differ by the long-polygon regions (42936, 112922, 124043, 472175, 780017, 1536355), by at most 0.3% of their tets. 112922 has 125 more tets under 1° (5,511) and 56 fewer face Steiner points; 780017 has 2 fewer tets under 1°.
      - Newly meshed: 37883, 45009, 55277, 98939, 129988, 135674, 472020, 822698 and 894834. The kernel is slowest on these (135674 at 638 s, 472020 at 505 s, 37883 at 493 s).
      - Over all 91, the kernel's median is 11.0 s against refinement's 10.9 s; 74 take under 60 s.
      - Six still time out: 104738 and 225958 in refinement, and 77922, 101548, 200965 and 1620057.
      - 39245 now leaves one input triangle unrecovered (it timed out before), and 500113 still does; 1501546 still keeps one tet turned over after rounding.
  - **Steiner points on facets** (Nov 11, 2026). Two causes left input triangles unrecovered:
    - Where a crossing edge has an implicit end, the point near the crossing (a BPT point from approximate positions) lies at a corner of the tets' sections by the plane, so on a facet that flips have partly covered it can fall on a protected face, and every insertion fails (39245, 500113). `ExactDelaunay::point_in_crossing_tet` puts it instead at the approximate centroid of the section of a tet around the edge, as a BPT point of the triangle, kept only if exactly inside that tet: then it is in the uncovered part and on no edge or face. The point near the crossing is still tried where that one does not go in; using only the new one lost 106779, 472175 and 1146192. Test: `points_in_crossing_tets_lie_on_the_plane_inside_a_tet`.
    - A large facet ran out of rounds, not points: 106779's, crossed by 505 edges, took every point but ended at 264 crossings after 64 rounds of at most 8. The points a round now double up to 32 (`PER_ROUND`), from 128 crossing edges tried.
    - The synthetic golden reports are unchanged.
    - **On the clean subset:** 92 of 100 mesh, conform and round (91 before). 39245, 500113 and 77922 are newly meshed; no model now fails face recovery.
      - The machine ran slower for this run: refinement, which did not change, took 28% longer at the median. Run back to back, the old and new kernels take the same time (332810, 42936, 119918). So 135674 and 472020 (638 s and 505 s before) went past 900 s, as did 200965, which meshes in 765 s alone.
      - 72 of the 89 meshed in both runs have the same Steiner counts; the other 17 differ in face Steiner points (24 fewer to 136 more), by at most 0.44% of their tets. The smallest face angle fell to 0.000° in 64992, 106779 and 822698 (from 0.021°, 0.080° and 0.027°) and rose in 132060.
      - Seven time out (104738 and 225958 in refinement; 101548, 135674, 200965, 472020 and 1620057), and 1501546 still keeps one tet turned over after rounding.
  - **Timeouts** (Nov 11, 2026). Profiled with `PROF_EVERY=<seconds>`, which prints the totals so far to stderr at that interval, for runs that never finish:
    - 1620057 spent its kernel inserting segment Steiner points: shrinking each cavity to protect the recovered faces took 139 ms, re-asking the same face's side and the same edge's ring of tets. Both are now remembered within one insertion (`Known`), at 5 ms a point, and the meshes are identical (`cdt_segments` on 40179, 76579, 1120771 and 42936).
    - 104738 and 225958 spent refinement locating points, at 2.5 to 4.3 ms a walk from the last tet made. Each split or subface point now starts its walk, in the volume and in the plane, from a vertex of what it splits (`ExactDelaunay::hint_vertex`, `PlaneDelaunay::hint`), and approximate positions are kept rather than recomputed. 225958's refinement fell from 1,264 s to under 500 s, 104738's to 1,076 s. The walk's start changes its random steps and so where new tets are stored, which changes the order refinement visits them: the refined points differ (332810: 32,742 against 32,811), and only 2 of the 90 models meshed in both runs keep their tet count (most within 1%, 70380 by 15%).
    - 101548 spent minutes per gift-wrapping call, on cavities of 1,400 to 9,400 vertices; the apex search weighs every vertex. A cavity over 1,000 vertices (`WRAP_VERTICES`) is now left to the Steiner points, and a facet crossed by more than 128 edges gets up to a quarter as many points a round, not 32, so it does not run out of rounds.
    - **On the clean subset:** 95 of 100 mesh, conform and round (92 before). 101548, 135674, 200965, 472020 and 1620057 are newly meshed; the kernel is slowest on 135674 (598 s), 37883 (491 s) and 472020 (367 s).
      - On the 90 meshed in both runs, the kernel's median went from 17.4 s to 11.4 s and its total from 5,574 s to 3,671 s; refinement's median from 15.7 s to 12.7 s.
      - 104738 and 225958 still time out. Alone, 225958 fails at 1,730 s: two single-triangle facets, each carrying about 74,000 refined points, are still crossed by about 660 edges after 64 Steiner rounds, and the kernel spends its time on their cavities (546 s), protected edges (225 s) and faces (217 s).
      - 77922 and 1224385 now keep tets turned over after rounding (2 and 1), as 1501546 still does (3): the different refinement gives rounding tets it cannot untangle. With the earlier build 1224385 meshes, and with every build from the new walk start on it fails.
  - **Rounding between near-coincident walls** (Nov 11, 2026). The tets left turned over in 77922, 1224385 and 1501546 were slivers, flat to within an ulp, between input facets less than an ulp apart (in 1224385 and 1501546 a layer of such tets, each with all its vertices within a few ulps of one plane). No single constructed vertex could move within 2 ulps without turning over another of its tets. Untangling (`cdt::delaunay::untangle`) now also:
    - moves one constructed vertex further, each coordinate by -r, 0 or r for r doubling from 4 to 64 ulps (`FAR_ULPS`), still under 1e-14 of the point's distance from the origin; this fixed 1224385 (one vertex by 16 ulps) and two of 1501546's three tets;
    - after the flips, moves two of the tet's constructed vertices together (every offset within one ulp, then doubling to 64), for at most 32 tets an untangling (`PAIR_TETS`); this fixed the rest (1501546 by 2 ulps, 77922 with one vertex by 32).
    - Two other repairs were tried and dropped. Replacing the tet and the inside tets around it by the star of a new f64 point fails: the whole neighbourhood is in the sub-ulp layer, so every candidate point is on its plane. Moves that only lessen the unusable tets around a vertex added nothing once the far moves were in.
    - Test: `offsets_go_from_near_to_far`. The synthetic golden reports are unchanged; `docs/report.md` says a move may reach 64 ulps and that a paired move gives two records.
    - **On the clean subset:** 98 of 100 mesh, conform and round (95 before); the other 95 meshes are identical. Rounding's median is 0.6 s; the three models take 0.5 to 7 s to round.
  - **Open** (carried to Milestone 15 and the kernel's speed work):
    - the two timeouts, both now in the enriched CDT: 104738's refinement ends within its budget (923,154 points of 1,201,640, 603 s alone) but the kernel takes over 900 s on that many points, and 225958's two large facets (above);
    - the kernel's time on the slowest models (135674 at 640 s);
    - the 0° tets left (591216);
    - slivers (dihedral angles reach 0 in 60 meshes; the paper's sliver criterion is optional).
- [x] Refinement against the sizing field, opt-in and under a deterministic budget (`refine::refine_with`, `RefineParams`, Nov 11, 2026).
  - **Size:** with a size field, a tet whose centroid is inside the input (the input's fast winding number over 0.5, asked only where the size wants a split) and whose circumradius exceeds the field at its centroid is bad too. It is queued by the larger of its radius-edge ratio and the bound times its circumradius over the size, so without a field the queue and the meshes are exactly the quality-only ones. Tets in the outer box are refined for shape only. `Refined::steiner_size` counts the points split for size alone.
  - **Budget:** a vertex budget and a step budget (insertion attempts, 64 a budgeted point), both deterministic. The vertex budget is `BUDGET_FACTOR` (8) times the meshing estimate at the element size (`sizing::point_budget`, now shared with the constrained kernel, which is unchanged). The estimate alone was too tight: shape refinement puts points where the input's features are small, which the size does not see, and 46 of 98 models stopped at it (tets under 1° up from 677,000 to 962,000); unbudgeted, models needed up to 6.5 times the estimate (98939) and 90 % of them 3.4 times or less.
  - **Hand-over:** at either budget refinement stops and returns the points it has, with `Refined::stopped` saying why (`RefineError::Budget` or `RefineError::Work`). The enriched CDT recovers the input from any set of points, so the mesh is only coarser or worse shaped where refinement had not reached; Milestone 15's local operations take it from there. In the pipeline this is to be recorded as `MESH_BUDGET_EXCEEDED`, a warning, when the chain is wired in.
  - Tests: `a_size_field_refines_inside_the_input` (a cube at size 0.25: no tet centred inside has a larger circumradius, the box keeps its coarse tets, and without a field the mesh is `refine_conforming`'s), `a_budget_stops_refinement_and_hands_over` (point and step budgets on a chamfered tetrahedron; the stopped mesh still makes an enriched CDT that rounds).
  - `--example enriched`: `REFINE_SIZE` (element size over the diagonal, 0.03 as in the pipeline) sizes the budget, `REFINE_SIZING=1` refines to the sizing field for it, and the output gives the points split for size, the steps, the budget and why refinement stopped.
  - **On the clean subset** (900 s each, eight runs sharing the machine):
    - Without sizing: the 98 meshes are identical to the unbudgeted run, and no model stops at its budget.
    - With sizing: 96 of 100. On the 96 meshed in both runs, sizing adds 1.6 % of tets in all and 4 % more tets under 1°; the smallest face angle rises on 11 models and falls on 8. 101548 and 500113 newly leave one input triangle unrecovered in the enriched CDT. No model stops at its budget.
- [x] Exit criterion, on the clean subset (met Nov 11, 2026; `enriched` example, without sizing, 900 s a model, eight runs sharing the machine):
  - **Refinement ends within its budget on every model.** None stops at it, with or without sizing; 104738's is the largest at 923,154 of 1,201,640 points. The two models that do not mesh, 104738 and 225958, time out in the enriched CDT, after refinement.
  - **Minimum face angle** (98 meshed): on the 18 inputs with no angle under 1°, the smallest face angle in the mesh has a median of 1.01° (1.08° with sizing), from 0.00° to 19.5°, and 7 reach the paper's lowest average, 1.3° (the paper reports 3.6° average minimum face angle, 1.3° for the CDT alone). The other 80 inputs have needles under 1°, which a conforming mesh cannot beat; 35 of the 98 meshes (38 of 96 with sizing) have no face angle below their input's smallest.
  - **Cost per model** (seconds; chamfering, refinement, enriched CDT, rounding):

    | | Median | 90 % | Largest |
    |---|---|---|---|
    | Whole chain | 31.1 | 268.5 | 716.4 |
    | Chamfering | 0.03 | | 2.0 |
    | Refinement | 13.7 | | 339.8 |
    | Enriched CDT | 13.3 | | 643.0 |
    | Rounding | 0.66 | | 116.9 |

    59 models take under a minute and 90 under five; the median mesh has 133,000 tets (2.4 million at most). With sizing the whole chain's median is 29.1 s (90 %: 227.7 s). The paper reports about 120 s a model and one 36-hour outlier.

**Milestone 15: Worst-cell quality for finite volume.** Done (Nov 13, 2026).

- [x] Invariant-guarded local operations (`mesher_tet::improve`, Nov 12, 2026). Each poor tet (smallest dihedral angle under 15° for now), worst first, gets the first kind of operation that improves it, the best of that kind:
  - **flips:** 2-3 flips of interior faces, and removal of interior edges with rings of up to 10 tets (3-2, 4-4 and larger), the ring polygon retriangulated by Shewchuk's maximin dynamic program;
  - **vertex moves:** a vertex on no boundary face moves to the best of the sliver smoothing's candidates (`smooth::best_position`, now shared, so smoothing is one of these operations);
  - **collapses:** a vertex on no boundary face merges into a neighbour;
  - **splits:** an interior edge splits at a new vertex, placed at the best candidate around its midpoint.
  - **Guards:** every new tet must be positive exactly and in floating point, and the worst tet of the region replaced must improve by more than 10⁻⁶°. Flips and splits act only on faces and edges inside the mesh, moves and collapses only on vertices off the boundary, so boundary faces and their vertices never change and the volume is kept. Eight passes at most; deterministic (poor tets by quality, then index).
  - **Pipeline:** runs on every conformal candidate (conforming and exact kernels) before the volume gate; face groups are renumbered. The synthetic golden reports were regenerated: minimum dihedral angles rise on 45 runs (most to 15°, the sliver-face cases from 2.1° to 8.8–10°), with tet counts within 1 % and mean non-orthogonality within 0.4°.
  - **Fixed in passing:** the fast smallest dihedral angle (`smooth::min_dihedral`) oriented each face normal by which side of its face the opposite vertex lay, which is noise for a tet flat in floating point, so such a tet came out at 180°, not 0°. Normals are now oriented by the vertex order (`a_flat_tet_has_no_angle`). The conforming kernel's meshes on 43853, 45009 and 113890 are unchanged.
  - Tests: `a_three_ring_edge_is_removed`, `a_face_is_flipped`, `improvement_keeps_the_invariants` (Delaunay meshes of random points: valid, boundary and volume kept, worst tet never worse, deterministic), `bad_indices_leave_the_mesh_unchanged`.
  - `--example improve` meshes as the pipeline does (`IMPROVE_KERNEL=exact` for the exact kernel), improves (`IMPROVE_TARGET`, `IMPROVE_SPLITS=0`, `IMPROVE_COLLAPSES=0`) and prints quality before and after.
  - **On the clean subset** (conforming kernel, the 53 models it meshes): the smallest dihedral angle rises on all 53 (median 0.82° to 7.18°; 30 reach 5°, from 17), tets under 1° fall from 6,805 to 1,082 (models with any: 27 to 10), tets under 10° from 65,796 to 55,683. Tets grow by 3.5 % at the median and 43 % at most (41106; 894834 by 34 %), almost all from splits: without them the smallest angle barely moves (43853: 0.0003° to 0.0004°, against 2.28° with them). The largest non-orthogonality falls on 30 models and rises on 11 (median 87.7° to 86.6°); the largest skewness rises on 19 (median 3.05 to 4.74), as the target is not yet a finite-volume one. Improvement takes 0.13 s at the median (12 % of meshing), 29 s at most (894834).
- [ ] A WSVM-style energy as the vertex-move objective, measured against the candidate search
- [x] Targets: per-face non-orthogonality and skewness, cells with a boundary face first (Nov 12, 2026).
  - **Quality** is what a cell-centred solver sees, each measure over its target (`ImproveParams`: smallest dihedral angle 15°, non-orthogonality 70°, skewness 2), so 1 is on target: a tet's smallest dihedral angle; each face's non-orthogonality (normal against the centroid line; on the boundary, against the line from the cell's centroid to the face's); each interior face's skewness. A cell's quality is its angle and its four faces', each capped at 1, so a cell on target is left alone. A target of 90°, infinity or 0° turns its measure off.
  - **Regions:** an operation is scored on its new tets and every face of them, against the old tets and theirs, so the cells across the cavity count. A vertex move's face structure is built once and scored at each candidate position.
  - **Acceptance:** the region's worst of the two measures (angle, faces) must rise, or stay while the other rises; neither may fall below its floor, the lower of its old worst and halfway from the region's worst to the target. So a measure the boundary fixes does not stop the other improving, and one measure is traded for the other only well above the worst. Tried and set aside: the worst alone (dihedral angles traded down for faces: tets under 10° 55,788 to 81,922), and a strict guard that no measure's worst may fall (no better on the subset, and 655 s on 894834 before the cheap rejection below).
  - **Order:** the poor cells with a boundary face first, then by quality and index.
  - **Speed:** a candidate's tets are checked (positive, angle above the floor) before its faces are looked up; that took 97659 from 50 s to 4.8 s.
  - **Measured:** `FvQuality` gains `wall_cells`, `wall_max_non_orthogonality_degrees` (every face of a cell with a boundary face, boundary faces by the centroid-to-face line) and `wall_max_skewness`. `--example improve` prints them, takes `IMPROVE_NONORTH`, `IMPROVE_SKEW` and `IMPROVE_FV=0`, and with `IMPROVE_WORST=1` describes the worst face left. `tests/oracle/examples/tetgen_quality.rs` meshes the same files with TetGen (`pzAfq`, the same size) and prints the same line.
  - Tests: `a_mirrored_face_is_on_target`, `face_targets_improve_the_worst_faces`, and the wall measures in `fv_quality_of_two_mirror_tets`. Goldens regenerated: records, kernels and acceptance unchanged, tet counts within 1 %, mean non-orthogonality within 0.2°; the three sliver-face runs' smallest angle falls from 8.8–10° to 3.5–6.9° (the angle-only target's moves there worsen a boundary face that cannot be fixed).
  - **On the clean subset** (the 53 models the conforming kernel meshes; medians of the per-model worst):

    | | Kernel | Angle only | Angle and faces | TetGen |
    |---|---|---|---|---|
    | Smallest dihedral | 0.82° | 6.00° | 3.87° | 4.67° |
    | Tets under 1° (sum) | 6,805 | 1,137 | 2,974 | 4,418 |
    | Tets under 10° (sum) | 65,796 | 55,788 | 49,499 | 80,192 |
    | Largest non-orthogonality | 87.7° | 86.1° | 81.9° | 76.0° |
    | Largest skewness | 3.05 | 4.05 | 2.43 | 1.99 |
    | Wall-adjacent non-orthogonality | 88.3° | 86.1° | 81.9° | 71.4° |
    | Wall-adjacent skewness | 2.21 | 2.81 | 2.01 | 1.16 |

    Against the kernel, the wall-adjacent non-orthogonality falls on 36 models and the skewness on 29. Against TetGen it is lower on only 12 models (skewness 16). Tets grow 5 % at the median (46 % at most); improvement takes 0.43 s at the median, 115 s at most. The sum of the largest skewness falls from 42,784 to 2,634.
  - **What is left at the walls:** the worst faces sit between two cells with three boundary vertices each and one free vertex low over a small boundary triangle (height 3–12 % of its edge), in parts 8 to 140 edges thick, so not thin features. No operation here may change the boundary, and TetGen splits boundary faces and segments. Boundary operations are the next task.
- [x] Reach the walls (Nov 12, 2026). Planned as boundary operations (coplanar 2-2 flips, boundary edge and face splits, in-plane vertex moves); what the worst faces turned out to need was different.
  - **Diagnosis** (`IMPROVE_WORST=1`): the worst wall faces' free vertex is the apex of a fan over a patch of small boundary triangles, much wider than they are (68370: one vertex in 240 tets, 120 of them standing on boundary triangles). The conforming kernel exempts tets whose shortest edge is under its size floor (h/16) from radius-edge refinement, so such patches keep one apex. And on skinny boundary triangles no tet is good. The vertex is boxed in: of 729 trial positions around it, 727 invert a tet; every collapse makes an unusable tet; no split or insertion improves its region.
  - **Built in `improve`:**
    - **insertion** above a wall cell's boundary faces (a point on the face's inward normal at 1, 0.7, 0.5 and 0.35 of a regular tet's height; its Bowyer–Watson cavity, at most 128 tets and never across the boundary, shrunk until the point sees every cavity face; the cone checked to fill the cavity by volume);
    - a **shortfall** criterion (the sum of each measure's squared shortfall below 1 over the region's tets and faces may fall while the worst holds, so a cluster whose worst face no single operation can reach is worn down);
    - **every kind scored together** when no flip raises the worst (vertex moves are now candidates like the rest, not applied on first gain).

    On the original surfaces, against the previous task's results: smallest angle median 3.87° to 4.52°, tets under 1° 2,974 to 2,574, wall-adjacent non-orthogonality 81.9° to 81.5°, skewness 2.43 to 2.39; tets 6.4 % more at the median (was 5.1 %). The sliver-face goldens recover (5.2–6.9° to 10.0–10.1°; one `VOL_LOW_QUALITY` warning gone).
  - **Tried and dropped:** off-centres in the conforming kernel, a point above each small subfacet of a poor tet under the floor. On 68370 it added 48,000 tets without moving the worst face (89.9°), and 97659 then failed boundary recovery.
  - **What works: a remeshed surface.** `IMPROVE_REMESH=1` remeshes the surface to the size field first (the pipeline's fallback path). On the subset, 49 of 53 mesh (45009, 148292, 225958 and 894834 fail recovery on the remeshed surface), and the medians beat TetGen's:

    | Median of per-model worst | Original surface, improved | Remeshed, improved | TetGen |
    |---|---|---|---|
    | Wall-adjacent non-orthogonality | 81.5° | 69.3° | 71.4° |
    | Wall-adjacent skewness | 2.03 | 0.99 | 1.16 |
    | Non-orthogonality | 81.5° | 68.2° | 76.0° |
    | Smallest dihedral | 4.52° | 15.0° | 4.67° |
    | Tets | 15,163 | 11,930 | 12,485 |

    Per model, remeshed and improved beats TetGen on wall-adjacent non-orthogonality on 28 of 49 and on skewness on 35 of 49 (the largest gaps: 93554, 59546, 233187). Improvement there takes 0.01 s at the median.
- [x] The pipeline's surface policy for quality (decided Nov 12, 2026: remesh when the input is poor).
  - **Considered:** remeshing every part first costs volume on curved parts (the synthetic icosphere 1 %, torus 2 %, failing the suite's volume check) and, on clean flat shapes, 37 % more tets for a slightly worse mean non-orthogonality. No measure of the input triangulation predicts the gain (skinny triangles are everywhere, in 12-triangle boxes that gain nothing too; almost no model has triangles far finer than the size field).
  - **Rule:** the repaired surface is meshed first; if that fails, or its worst cell is beyond the quality limits (`worst_cell` in `volume.rs`: smallest dihedral angle under `SLIVER_DEGREES`, non-orthogonality over 70°, skewness over 2, the `ImproveParams` targets), the remeshed surface is meshed too, and the better worst cell goes to the volume gate first (the remeshed one only if strictly better). See fallback ladder step 2. Records: a remeshed mesh kept for quality gets `VOL_LOW_QUALITY` info `remesh/applied` with both worst cells; findings deferred to remeshing become `applied` (`RecordLog::deferred_applied`, test `a_deferred_remesh_becomes_applied`). `--no-remesh` (and the GUI switch) now means "keep the input's triangulation".
  - **Goldens:** clean shapes keep their exact surfaces; the three sliver-face runs now remesh (smallest angle 3.5–6.9° to 15.0°, `QUAL_SLIVER_FACE` `remesh/applied`, the dihedral warning an info record).
  - **On the clean subset with the CLI** (`run_subset.py`, 8 jobs): all 100 meshed, 94 conformal (67 of them remeshed) and 6 by stuffing; exit codes 0: 13, 10: 81, 20: 6; 3,308 s in all. On the 94 conformal meshes the largest non-orthogonality has median 69.4° and skewness 1.47, the smallest dihedral 15.0°. Volume: exact on the 27 kept surfaces; on the 67 remeshed ones 0.07 % off at the median, 2.6 % at most. Against TetGen on the 53 models both mesh (per model, whole-mesh worst): non-orthogonality lower on 40, skewness lower on 39; medians 69.3° against 76.0° and 1.46 against 1.99.
- [x] Report the worst cells and their locations, wall-adjacent cells separately; `VOL_LOW_QUALITY` gated on non-orthogonality and skewness too (Nov 12, 2026).
  - **Located:** `fv_quality` keeps where each worst value is (`WorstAt`: the cell, the cell across a face, and the cell's centroid or the face's): the smallest dihedral angle, the largest non-orthogonality and skewness over interior faces, and both over the faces of wall-adjacent cells. Ties go to the lower cell indices, so the boundary faces' hash-map order does not show.
  - **Summary:** `MeshSummary::worst_cells` lists them (`WorstCell`: `measure`, `wall`, `value`, `element` and `neighbour` as 1-based numbers in the part's element block, which the writers keep in order, and `point`).
  - **Gate:** one `VOL_LOW_QUALITY` warning per measure beyond its limit (dihedral angles against the kernel's bounds as before; non-orthogonality over `MAX_NON_ORTHOGONALITY_DEGREES`, 70°; skewness over `MAX_SKEWNESS`, 2, the improvement targets and `worst_cell`'s limits), each located at the worst cell (`element`, `neighbour`, the location's `point`, which the GUI shows).
  - Tests: `worst_values_are_located_independently_of_tet_order` and the located worst values in `fv_quality_of_two_mirror_tets`; `tests/synthetic/tests/worst_cells.rs` meshes the sliver-face cube with `--no-remesh` (the dihedral and non-orthogonality warnings are located within an element of the sliver's apex, and are the summary's worst cells) and without (remeshed, no warning). The suite's goldens are unchanged: no synthetic part ends beyond the limits.
  - **On the clean subset** (from the last CLI run's summaries; the meshes are unchanged): 35 of 100 parts get a quality warning, 27 before (dihedral angles only): 33 for non-orthogonality, 15 for skewness, none of them stuffing meshes. Exit codes are unchanged (the 8 newly warned parts had warnings already).
- [x] In parallel, with the serial result (Nov 12, 2026). Planned as independent sets with a fixed colouring; that would change the result with the colouring, so instead the serial pass is kept and sped up.
  - **Scheme:** a pass's poor cells go in batches (2 a thread). Each cell of a batch is evaluated in parallel against the mesh at the batch's start (`Mesh::evaluate`, now `&self`: splits no longer push their vertex), recording what it reads: tets (vertices, quality, wall flag) and vertex stars, through accessors and a thread-local recorder that is off otherwise. The operations are applied in the serial order; one whose reads include a tet or star changed by an earlier operation of the batch (stamped per batch on apply and refresh) is evaluated again first. Points need no record: a point is read only for a vertex of a tet read or a star read, and a move refreshes the vertex's whole star. A new vertex is a placeholder (`NEW`) until applied, then gets the next index, as serially. One thread runs the plain serial loop. `smooth::best_position` takes a position lookup instead of the point array.
  - **Identical:** on the 53 conforming subset models, every line of `--example improve` (tets, quality, operation counts) matches the serial code of 67e3603; the CLI's outputs (Exodus, msh, vtu, defects, repaired OBJ, and the report less its timings) are byte-identical at 1, 4 and 16 threads on 41106, 43853, 68370 and 97659. Test: `the_same_mesh_at_any_thread_count` (random Delaunay meshes, 1, 2, 4 and 16 threads, points, tets, map and statistics equal).
  - **Speed** (32 threads, i9-13900K): improvement on the subset 353 s in all to 100 s; 2.0× at the median, 5.6× on 894834 (137 s to 25 s). With one thread it is as fast as before (894834: 129 s), so recording costs nothing when off. Batch size per thread on 894834: 1, 2, 4, 8, 16 took 28, 24, 24, 27, 33 s. What is left serial there: re-evaluating stale cells (6 s), applying and refreshing (3 s), and each batch waits for its slowest evaluation.
- [x] Absorb improvement 8 (remeshing; below). Improvements 5 and 6 (early stop, point-budget overruns: the worker timeouts on 39245, 101548 and 472078 and the budget hits on 132060 and 677341) are deferred to the future by decision (Nov 13, 2026); they stay under "Things to consider improving".
  - [x] Remeshing keeps feature vertices on their own line (Nov 12, 2026; improvement 8's "remeshing can"). Of the 10 subset models with the worst cells or a stuffing fallback, the remeshed surface failed on all: on 4 it crossed itself (the surface gate), on the rest a few triangles were not recovered. The crossings were at thin slots: a feature vertex was projected onto the nearest of up to 512 feature segments reachable from its own, and the search passed through corners, so it could jump onto the parallel edge of the slot's other wall, and its triangles spanned the gap. `Input::project_feature` now never passes a corner (an input vertex not on one smooth line: `Input::line`), and a feature edge between two corners searches from every segment at the corner. Test: `a_feature_vertex_stays_on_its_own_line`. A remeshed surface that fails the surface gate now says which checks and how often (`gate_failures`; `remesh_retry` reads, for one, `the remeshed surface failed the surface gate: 3 self-intersecting triangle pairs`).
    - **Tried and dropped:** splitting the edges of crossing triangles at midpoints projected onto the input, round after round. The crossings grew each round on every failing model (4 to 20,380 on 1146192), since they were not chords sagging through a thin wall.
    - **On the subset with the CLI:** 95 models unchanged; 290506 now meshes conformally (remeshed; it was stuffing); 148292 and 894834 now take the remeshed surface (skewness 12.1 to 6.1 and 317 to 47; tets 191,082 to 37,788 and 679,911 to 108,816); 95496's remeshed surface now passes the gate but one triangle is not recovered. 101548 and 472078 went to stuffing on the worker timeout (600 s), with their meshes unchanged when they finish: alone they take 665 s and 597 s, and timings on this machine vary by about 40 % (132060 alternated with the build before this change: 205 and 278 s against 289 and 292 s, the same mesh). So the subset is 93 conformal as run, 95 without the two timeouts.
    - **Recovery failures on remeshed surfaces** (investigated Nov 12, 2026; no change). On 95496, 45009 and 39245 the conforming kernel cannot split a missing subsegment between two of its own Steiner points: every split point's forced insertion is refused because the cavity shrinks to nothing. The point lies on a remeshed triangle the kernel has already put Steiner points on, several of them packed along the edge between two nearly coplanar remeshed triangles (each its own facet), so the tet holding the point has a face in the point's plane, and the tets around it are held by the recovered pieces of that edge (constrained edges, which forcing does not override). On 39245 the cascade reaches an edge 1.4e-14 long. Growing the cavity across in-plane faces (for the held tets, then for any tet of a forced insertion) did not help and was reverted. Remeshing at 0.8 of the size lets 5 of the 6 failing models mesh, but with non-orthogonality 88–90° and skewness 9 to 84 (in the example, which does not check for self-crossing), while the stuffing they get now has 61–66° and 0.6–0.8 at 0.35–1.6 % boundary deviation. So recovering these surfaces would not improve their worst cells. 225958 and 45009 keep their own surfaces with far worse cells (skewness 1,310 and 9.5) than stuffing would give: whether to prefer stuffing when every conformal mesh is far beyond the limits is a policy question.
  - [x] Stuffing for conformal meshes far beyond the limits (decided Nov 12, 2026, after the recovery investigation above). When the best conformal candidate's `worst_cell` is under 0.25, stuffing is meshed and gated into a scratch log; it is kept if it passes the volume gate, its worst cell is better and its volume is within 1 % of the surface's (see fallback ladder step 3). Test: `far_beyond_the_limits_stuffing_is_kept` in `tests/synthetic/tests/worst_cells.rs` (the sliver cube, its apex 1e-5 of the way across, with `--no-remesh`). Goldens unchanged.
    - **On the subset with the CLI:** 15 models switch from conformal to stuffing (all 16 under 0.25 but 135674, whose stuffing mesh lost 12.5 % of the volume, hence the volume rule). Their worst non-orthogonality falls from 82–90° to 61–67°, skewness from 1.8–1,310 to 0.58–0.79; the boundary deviation is 0.27–2.9 % and the volume 0.01–0.65 % short. No mesh on the subset is left with a worst cell under 0.25 but 135674. Medians over all 100: worst non-orthogonality 69.3° to 68.3°, skewness 1.44 to 1.40. Exit codes 0: 13, 10: 65, 20: 22 (20 was 7). 225958's graded stuffing mesh has a smallest dihedral angle of 8.35°, under the proven 10.7° (a `VOL_LOW_QUALITY` warning; the graded lattice's bound is open).
    - **Left:** 1146192, 129988 and 98939 still cross: a smoothed vertex beside a cluster of tiny sharp corners is projected onto a nearby but different patch, and its triangles cut the input triangles there. That needs each remeshing operation checked against the triangles near it.
  - [x] Remeshing merges features too small to mesh (Nov 13, 2026, after the exit-criterion measurement below found 11 models beyond a limit and worse than TetGen).
    - **Diagnosis** (`IMPROVE_WORST=1` now also describes a worst cell with every vertex on the wall: its boundary faces, the angle between them, their smallest angles, its corners). On most of the 11 the worst face is on such a cell, which no vertex move reaches, standing on needle triangles of the remeshed surface (smallest angles 1.2–9°) beside an input feature far smaller than the element size: on 233187 a jog of 0.02 between two corners (h = 1.85), on 93554 a 45° strip 0.019 wide between two feature lines, on 59546 a wall thinner than the size floor. The remesher never moved a corner or collapsed across two feature lines, so these survived at any size. 132060 and 822698 are genuinely sharp wedges (9.3° and 14.9°).
    - **Tried and dropped:** lower size floors (the size field's h/4 to h/16, the kernel's `min_boundary_size` h/16 to h/64; `IMPROVE_SIZE_FLOOR` and `IMPROVE_BOUNDARY_FLOOR` in the example): on the 11 as many models got worse as better (77922 70.8° to 83.6°), and 233187 did not move.
    - **Rule** (`RemeshParams::small_feature`, 0.1): a corner collapses into another corner along a feature edge shorter than a tenth of the target size, and a feature vertex into a vertex of another feature line closer than that (the two sides of a strip that narrow merge into one line). Each moves the surface by at most the edge's length; the volume gate measures the remeshed surface's deviation as before. `Remeshed::merged_small_features` counts them, and the remesh records carry it (`merged_small_features`). Test: `a_feature_too_small_to_mesh_is_merged` (a cube with a 0.005 chamfer at size 0.2: smallest triangle angle under 5° kept, over 20° merged, volume within the chamfer's width). The synthetic goldens are unchanged.
    - **On the subset with the CLI** (against the measurement below): 18 models change. Corners alone took the list of models beyond a limit and worse than TetGen from 11 to 7; with strips too, 8 models that fell back to stuffing now mesh conformally on the remeshed surface (exit codes 0: 13, 10: 73, 20: 14), and the list is 45009, 77922, 97659, 113890, 132060, 263198, 822698. 233187 goes from 80.4° to 63.3°, 59546 from 81.5° to 68.5°, 1501546 from 87.0° to 73.5°, 119918 from 79.1° to 69.0°, 43853 from 80.3° to 71.0°. The volume error of the changed models stays at 0.35 % or less (135674's 2.6 % is unchanged), the boundary deviation at 1 % or less. Medians over the 92 TetGen meshes: wall-adjacent non-orthogonality 69.7° against 77.3°, skewness 1.00 against 1.65; the CLI is lower on 68 and 75 of 92 (65 and 73 before).
    - **Left:** 45009 and 97659 are now conformal with cells within the far-off limit but worse than their stuffing meshes were (79.6° and skewness 2.28, against 70.0° and 0.73), as the far-off rule intends; 132060 and 822698 are sharp wedges; 113890, 263198 and 77922 are not yet diagnosed (113890's worst cell stands beside a 0.008 offset that was not merged).
  - **The remaining 7 diagnosed** (Nov 13, 2026; no change), from the CLI's own meshes (each worst cell read from `mesh.msh` by its `worst_cells` element):
    - *Sharp wedges:* 132060 (9.3°), 822698 (15°) and 113890 (17°, not the offset). The worst cell stands on both walls at the wedge's edge, all four vertices on the wall, flat across the wedge. Every cell at such an edge has a dihedral angle no larger than the wedge's, so improvement cannot split it without lowering its angle further, and good faces there need cells graded down towards the edge. TetGen leads by 1.3–4.3° (79.5° against 83.8°, 75.9° against 77.2°, 76.0° against 79.3°).
    - *Features between a tenth and half the local size:* 263198 (an input edge of 0.079 beside a 90° edge; the local target is about 0.44, so the merge threshold is 0.044), 97659 (a strip 0.25 wide) and 45009 (strips 0.1 wide, with needle triangles of 1.5° beside a 0.0045 tilted strip). Raising `small_feature` to 0.2–0.5 (`IMPROVE_SMALL_FEATURE` in the example) helped 263198 (73.8° to 69.1°) and 1501546 but made 43853 worse (72.6° to 75.7°), 97659 moved either way and 45009 not at all, so 0.1 stays. TetGen meshes 45009 with 285,405 tets against 29,223 for its 75.2°. On 97659 the CLI's non-orthogonality is far better than TetGen's (75.7° against 88.2°); only its skewness is worse (2.28 against 2.02).
    - *A tie:* 77922 (76.73° against 76.7°), an ordinary interior face over a free vertex.
    - The stuffing meshes 45009 and 97659 had before strip merging were better than both (70.0° and 71.6°): preferring stuffing whenever a conformal mesh is beyond the limits and stuffing within them would be a policy change from the far-off rule.
  - [x] Sharp wedges: untwisted cells along sharp edges (Nov 13, 2026).
    - **Why grading cannot help:** the cell standing on both walls of a sharp edge (the edge's segment and the two triangles' far vertices) is flat at any size, so its faces' non-orthogonality does not change as cells shrink towards the edge. What decides it is where the far vertices sit along the segment. At a 15° wedge: 77.5° when they sit at opposite ends (twisted), 72.1° over the same end, 67.7° both over the middle (82.2°, 78.8° and 76.1° at 9.3°; 75.9°, 69.8° and 64.9° at 17°). 822698's worst cell was twisted at 77.2°. On the walls beside a sharp edge the remesher's vertices often line up across the wall, and the two walls' diagonals then go opposite ways.
    - **Tried and dropped:** insertion candidates inside a wall cell (a quarter and half of the way from a boundary face's centroid to the opposite corner, between two boundary faces' centroids, the centroid), since the candidates above a face lie outside a thin wedge: at most 0.5° on the wedge models.
    - **Rule** (`align_sharp_edges`, the remesher's last step): along each feature edge whose solid angle is under 60°, each wall's triangle on the segment may flip the edge from an end of the segment to its far vertex, making the next triangle's far vertex the segment's, if the flipped triangles keep at least half their smallest angle; the combination of flips that brings the two far vertices nearest each other and the middle is applied. Twisted segments: 132060 16 to 6 of 68, 822698 20 to 3 of 52, 113890 8 to 0 of 32. `try_flip`'s checks are now `flippable`, shared. Test: `sharp_edges_are_not_twisted` (prisms with 15°, 20° and 30° edges, 11, 5 and 3 of 16 segments twisted without the rule, none with it).
    - **On the subset with the CLI:** 12 models change; medians the same (wall-adjacent non-orthogonality 69.75°, skewness 1.00). Better: 822698 77.2° to 75.8° (now below TetGen's 75.9°), 77922 76.7° to 72.3°, 263198 75.8° to 73.8°; worse: 97659 75.7° to 79.5° (its skewness 2.28 to 1.90, below TetGen's), 472131 68.0° to 70.5°, 1224385 68.4° to 69.8°, 113890 79.3° to 80.4°. Lower than TetGen on wall-adjacent non-orthogonality on 70 of 92 (68). Beyond a limit and worse than TetGen: 45009, 77922 (now on its skewness, 2.19 against 1.28), 113890, 132060, 263198.
    - **Left:** on 132060 and 113890 the worst cell is still at the wedge, where the edge's segments are shorter than the wall triangles (0.126 against 0.28 on 132060), so the far vertices sit beyond the segment's ends and no flip brings them over it. That needs the feature edge's vertices spaced like the walls', or the walls' like the edge's.
- [x] Exit criterion (met Nov 13, 2026, on the medians by decision; the per-model gaps left are recorded above and below):
  - on the clean subset, worst wall-adjacent non-orthogonality and skewness below the current kernel's and TetGen's (the medians of each model's worst: 69.75° and 1.00 for the CLI at aa89043, against TetGen's 77.3° and 1.65 over the 92 models it meshes, and the kernel's 88.3° and 2.21; the CLI is lower than TetGen on 70 and 75 of 92; beyond a limit and worse than TetGen on 5: 45009, 77922, 113890, 132060, 263198)
  - worst cells located (done: `worst_cells` and located `VOL_LOW_QUALITY` warnings)
  - byte-identical results at 1, 4 and 16 threads (done: see the parallel task)
  - **Measured Nov 13, 2026** (at dd98e8a). The CLI (`run_subset.py`, 8 jobs: all 100 meshed, 78 conformal, 22 by stuffing; exit codes 0: 13, 10: 65, 20: 22; 472078 timed out at 600 s this run), the conforming kernel without improvement (`--example improve`'s first measure, 53 of 100 meshed) and TetGen (`tetgen_quality`, 92 of 100), each model's worst over its wall-adjacent cells:

    | Median of per-model worst | Kernel | TetGen | CLI |
    |---|---|---|---|
    | Wall-adjacent non-orthogonality (53 all three mesh) | 88.3° | 71.4° | 69.7° |
    | Wall-adjacent skewness (53) | 2.21 | 1.16 | 1.03 |
    | Wall-adjacent non-orthogonality (92 CLI and TetGen mesh) | | 77.3° | 69.8° |
    | Wall-adjacent skewness (92) | | 1.65 | 0.95 |
    | Non-orthogonality (53) | 87.7° | 76.0° | 68.4° |
    | Skewness (53) | 3.05 | 1.99 | 1.38 |
    | Smallest dihedral (53) | 0.82° | 4.67° | 15.0° |
    | Tets (53) | 14,863 | 12,485 | 12,507 |

    The medians are below both on both measures, the non-orthogonality only 1.7° below TetGen's. Per model, the CLI is lower than the kernel on 46 of 53 (non-orthogonality) and 42 of 53 (skewness), than TetGen on 65 of 92 and 73 of 92, than both on 26 and 34 of 53. Most models where TetGen is lower are within the limits on both meshes, since improvement stops at its 70° and 2 targets: beyond a limit, 43 of 53 kernel meshes, 64 of 92 TetGen meshes and 33 of 92 CLI meshes. The CLI is beyond a limit and worse than TetGen there on 11: 43853, 59546, 77922, 93554, 113890, 119918, 132060, 233187, 263198, 822698, 1501546 (wall non-orthogonality 76–87° against TetGen's 70–90°; 93554 is stuffing at 71.1° against 70.1°). TetGen fails on 8 the CLI meshes (39245, 129988, 135674, 196190, 360073, 472131, 780017, 1620057).

**Milestone 16: Exact arrangement (repair rung 2).** Done (Nov 13, 2026): exit criterion met on 448 of the strata's 450 models; the two worker timeouts are deferred by decision with improvements 5 and 6.

Staged plan (decided at the start; the outline tasks below are ticked as the stages deliver them):

1. *Intersections as implicit points* (`mesher_geom::arrangement`). For each pair of triangles a BVH finds and the exact tests confirm, the points and segments where they meet: a vertex on the other triangle (explicit), an edge through the other triangle (LPI: the edge's line and the triangle's plane), and, for coplanar pairs, edges crossing (LPI with an auxiliary plane through the other edge and parallel to the projection axis, built from explicit points so it is exact). Points are deduplicated exactly (`implicit::coincident`) and indexed globally; each triangle gets its constraint points and segments.
2. *Triangle subdivision.* Each intersected triangle is retriangulated in its projection (exact `orient2d_projected`) with its constraints; where two constraint segments cross inside it, the crossing is the TPI of the three planes. Coplanar triangles' shared sub-triangles are matched, so each piece of surface appears once with its multiplicity.
3. *Labelling.* The winding number of the input on each side of each piece, from patches (pieces joined across uncut edges) and one exact ray count per patch, the generalized winding number where the input is open; keep the pieces between a side inside (at least 1/2) and a side outside. This resolves self-intersections, overlapping shells and nested shells (the inner one is dropped), and touching shells mesh as one region (improvement 9). Values near 1/2 are reported (`TOP_AMBIGUOUS_INSIDE`, added with this stage). *Changed when the stage landed (Nov 13, 2026):* closed inputs are labelled exactly (one exact segment cast per patch, integer windings), so nothing is ambiguous there; open inputs are refused by this stage and `TOP_AMBIGUOUS_INSIDE` waits for the stage that decides what to do with them (stage 5).
4. *Output.* The kept pieces rounded to f64 (snap rounding: implicit points to their nearest f64, coincident ones merged), then re-checked by the surface gate, since rounding can bring intersections back. The exact kernel taking implicit vertices directly is evaluated after.
5. *Pipeline.* Repair rung 2 runs on a part the gate rejects for `INT_SELF_INTERSECTION` or `INT_SHELL_OVERLAP`; records give the rung that succeeded, and the components and Euler characteristic before and after (a new action for the resolution); the resolved surface then goes through the volume pipeline like any other.
6. *Checks.* Synthetic intersection cases (the pushed-through icosphere, overlapping cubes, a cube through a cube, nested and touching boxes, coplanar overlaps) against a brute-force count of components and Euler characteristic and the exact union volume; a routine 100-model subset of the self-intersecting and multi-component strata (`tests/thingi10k/subset_dirty.txt`), with the strata's full 450 models for the exit criterion (a long run, on the user's go-ahead).

- [x] Arrangement of each part's triangle soup with LPI and TPI points (Cherchi et al. 2022); Guo and Fu's indirect offset predicates (2024) as the faster variant not evaluated, deferred at the close (Nov 13, 2026): the arrangement is not the bottleneck (median 0.05 s on the routine subset)
  - [x] Stage 1, intersections as implicit points (`mesher_geom::arrangement::intersections`, Nov 13, 2026). Per intersecting pair (the BVH's box overlaps, then `faces_intersect`): in general position, each triangle's boundary against the other's plane (vertices on it; edges through it as LPI points, one construction per edge and plane), keeping the candidates `locate` finds in the other triangle (exact `orient2d_projected` in its projection; corner, edge or inside); coplanar pairs give each triangle the parts of the other's edges inside it, edge crossings as the LPI of one edge's line and the plane along the other edge parallel to the projection axis (third point the edge's end moved along that axis only, so the plane is exact; one construction per pair of edges). At most two distinct ends per piece by convexity, deduplicated exactly. Points on an input edge go to both of its triangles. All points are then merged exactly (interval bounds swept along x, `coincident` on overlapping boxes; the lowest index wins, so input points keep theirs), and points no constraint uses are dropped. Output: the points (inputs first), each triangle's `Constraints` (points other than its corners, segments), and the pair counts. Tests in `crates/mesher-geom/tests/arrangement.rs`: a crossing pair, a triangle resting on another, coplanar overlap, the suite's overlapping boxes (points on both boxes' faces; one closed curve, every point of degree 2, though many points land exactly on the cubes' face diagonals), points on an edge reaching both its triangles, a degenerate triangle refused.
    - **Measured** (`cargo run --release -p mesher-geom --example arrange -- <obj files>`, `ARRANGE_DUMP=1` lists curve ends) on the new routine subset `tests/thingi10k/subset_dirty.txt` (100 of the 1,000-model sample's self-intersecting and multi-component strata, `random.Random(16)`: 67 and 33; median 6,504 faces, at most 99,980): all 100 run, median 0.03 s, at most 8.6 s (424211); 67 have intersecting pairs, 47 coplanar ones. On closed models without coplanar pairs every curve end is an input vertex resting on another face (touching shells: 81262's one segment is an edge of one shell lying flat on a face of the other), none a constructed point.
  - [x] Stage 2, triangle subdivision (`mesher_geom::arrangement::arrangement`, Nov 13, 2026). Stage 1 now records what each segment lies on besides its triangle (`Support`: the other triangle's plane, or the line of a coplanar triangle's edge) and the coplanar pairs. Coplanar overlapping triangles form clusters (union-find over the coplanar pairs) subdivided together, so a shared region becomes the same pieces in each by construction, not by matching tie-breaks between separate triangulations (axis-aligned faces make cocircular rectangles everywhere). Per cluster, in the projection along the first member's `plane_axis`: the members' corners and constraint points; their edges (on explicit lines) and constraint segments; each proper crossing of two segments built from their supports (two planes: the TPI with the cluster's plane; a line and a plane: their LPI; two lines: the LPI of one with the plane along the other parallel to the axis), merged exactly; each segment split at the points exactly on it (`orient2d_projected` zero, `dot_sign` negative), ordered by `dot_sign`; then a constrained Delaunay triangulation (`arrangement::cdt`: an explicit super-triangle framed by the members' corners, a walk with a full search as fallback, Lawson flips by `incircle_projected`, segments inserted by removing the triangles they cross and filling each side by Delaunay choice after Anglada). Each triangle goes to every member it lies in (all three corners in the closed member), oriented and attributed like it. Points made in several clusters (a triple point in each of its three triangles) are merged exactly at the end. Output: the points (inputs first), the pieces, and each piece's input triangle.
    - **`Point::approx` made accurate** (in passing): the f64 evaluation of an LPI whose line is nearly parallel to its plane, or a TPI of planes nearly through one line, could land far from the point (on 424211, rounded positions of the arrangement's points were far enough off to misplace whole pieces in the area check). It now keeps the f64 evaluation only when its error bound puts it within 2^-42 of the configuration's size, and otherwise divides the exact coordinates (`Exact::quotient`, within an ulp). Test: `approximate_positions_are_close_when_nearly_degenerate` in `implicit_oracle.rs` (fails without the fallback). `num-traits` became a normal dependency of `mesher-geom` (it was a dev-dependency; `num-bigint` already depends on it).
    - Tests in `crates/mesher-geom/tests/arrangement.rs`: each input triangle's pieces cover its area and face its way, and on closed inputs every edge of the arrangement is used as often in each direction (balanced): a crossing pair, overlapping boxes (in general position, with points on the cubes' diagonals, and sharing two face planes), a coplanar overlap whose shared region is the same pieces in both triangles (area 2 exactly), an octahedron with a vertex pushed to five places, and 80 random sets of three overlapping boxes or tetrahedra (half on a coarse grid, so coplanar and collinear contacts occur, welded where corners coincide).
    - **Measured** (`--example arrange`, now the whole arrangement: pieces, unbalanced edges in and out, the largest error in a triangle's covered area over the squared diagonal, `ARRANGE_DUMP=1` for the triangles that miss) on `subset_dirty.txt`: all 100 run, median 0.05 s, 26.5 s in all, at most 11.5 s (424211, 1,238 pairs, 635 coplanar); every one of the 68 closed inputs gives a balanced arrangement; area errors at most 1.7e-16 of the squared diagonal; median 7,351 pieces from 6,504 triangles. Before the fix above the frame of the super-triangle came from the points' interval bounds, which are unbounded for such ill-conditioned points, and 424211 failed.
- [x] Cells labelled by winding number, keeping the solid region (the outline task below; stage 3 does it for closed parts, the pipeline wiring is stage 5)
  - [x] Stage 3, labelling (`mesher_geom::arrangement::resolve`, Nov 13, 2026). For a closed, oriented input (every edge used as often in each direction; otherwise `ArrangementError::Open`): identical pieces merge into faces with a multiplicity (+1 per copy facing the first's way, −1 per opposite copy; multiplicity 0, between solids touching face to face, is dropped); faces joined across an edge they alone share in opposite directions form patches, with one winding number behind and one in front each; per patch, one exact cast from a point strictly inside one of its faces (a barycentric point of the face's input triangle, checked with `orient2d_projected`) to an explicit point far outside, counting each face the segment crosses by its multiplicity (exact `orient3d`; a segment touching an edge or a vertex, or with its far end in a face's plane, is cast again in another of 12 directions; a float slab test on rigorous boxes prunes the faces). A face is kept when exactly one side has winding at least 1, facing out. Patches no point fits in (needles an ulp wide: 40923 has two faces of area 3e-18) get their windings around an edge from a known patch: the faces around it ordered exactly (`radial`: each half-plane against the first face's plane, then by `orient3d` within its half-turn), each crossing changing the winding by its multiplicity, the way set by its orientation; a disagreement around an edge is an error. `resolve_with(.., false)` casts once per connected piece and finds every other patch that way.
    - Tests: overlapping boxes resolve to their union (volume 1.71875 exactly), a nested box is dropped, boxes touching face to face become one (the shared region gone), an inside-out box inside another becomes a cavity, an open surface is refused, and 60 random sets of overlapping boxes and tetrahedra (half on a grid) give balanced surfaces whose inside agrees with the input's winding number at 200 sample points each (brute-force solid angles, points near either surface skipped), and the same surface when found around edges.
    - **Measured** (`ARRANGE_RESOLVE=1` in `--example arrange`; `ARRANGE_PROPAGATE=1` compares with the around-edge labelling) on `subset_dirty.txt`: 68 closed models resolve, 32 open ones are refused; median 0.06 s, at most 2.0 s; every result balanced; the around-edge labelling gives the same surface on all 68. Rounded to f64 (a preview of stage 4), one model gets flat triangles and crossings back (904481: 4 and 21). Volumes against the input's signed volume range from 0.54 (70558: heavily overlapping shells, counted twice by the input) to 1.59 (96066: an inside-out shell sticking out of the solid; in the pipeline rung 1 flips such shells first).
- [x] The resolved surface goes to the exact kernel, or is snap-rounded and re-checked (the outline task below; stage 4 rounds; the exact kernel taking implicit vertices directly is deferred at the close, Nov 13, 2026)
  - [x] Stage 4, rounding (`mesher_geom::arrangement::resolve_rounded`, Nov 13, 2026). Each round: every point to its f64 position (`Point::approx`, within ulps), points landing on one position merged, collapsed triangles removed and opposite copies of one triangle cancelled; each flat triangle (corners exactly collinear) removed and its neighbour across the long edge split at the middle corner (the middle one found by an exact `dot_sign`); then the f64 surface checked exactly for triangles meeting beyond their shared vertices and edges (`intersections`), and if any do, resolved again (stages 1 to 3 on explicit points) and rounded again, at most 6 times (`ArrangementError::Rounding` after). Every step keeps each edge balanced; each triangle keeps its input triangle through all rounds. `Rounded` reports the rounds, points merged, triangles removed and flat triangles repaired.
    - Tests: overlapping boxes round to their union with no further rounds (volume exact); 60 random sets of overlapping solids (a quarter on a grid) round to closed surfaces with no flat triangle and no intersecting pair, volume within 1e-9 of the exact resolution's, origins valid; a flat triangle closing a tetrahedron whose face was split at an edge's midpoint is repaired (unit test in `round.rs`).
    - **Measured** (`ARRANGE_ROUND=1` in `--example arrange`) on `subset_dirty.txt`: all 68 closed models round to clean surfaces, median 0.09 s (resolution included), at most 1.8 s. Only 904481 changes in rounding: 2 points round onto one, 4 triangles collapse, and nothing crosses after (the preview's 4 flat triangles and 21 crossing pairs came from leaving those unmerged); no model needs a flat repair or another round. The volume of every model is unchanged by rounding to the printed precision.
- [x] Cells labelled by generalized winding number, keeping the solid region. This resolves `INT_SELF_INTERSECTION`, `INT_SHELL_OVERLAP`, `PART_OVERLAP` and touching shells (improvement 9). Done for closed parts with exact integer winding numbers (self-intersections, overlapping, nested and touching shells); the generalized winding number of open parts moves to Milestone 17 with them, and `PART_OVERLAP` (between parts, which rung 2 resolves one at a time) stays open (decided at the close, Nov 13, 2026)
- [ ] Winding numbers by the Antipodal method's ray-crossing form, where faster. Deferred at the close (Nov 13, 2026): one exact cast per patch is already fast (labelling at most 2.0 s on the routine subset)
- [x] The resolved surface goes to the exact kernel, or is snap-rounded and re-checked (snap-rounded, stage 4; as above)
- [x] Report which rung succeeded, and the components and Euler characteristic before and after (stage 5, below); proposed `TOP_AMBIGUOUS_INSIDE` deferred by decision (Nov 13, 2026): rung 2 takes closed parts only and labels them with exact integer winding numbers, so nothing is ambiguous there; ambiguity belongs to the winding-number field of open parts, which go to stuffing now and to wrapping in Milestone 17
  - [x] Stage 5, the pipeline (Nov 13, 2026). In `repair_part`, after rung 1 and its guards, a part with self-intersecting pairs or overlapping shells goes through `resolve_rounded`. On success its surface is replaced: triangles keep their model triangle (`origin`), vertices their model vertex where they did not move; non-manifold edges and vertices the union leaves where shells touch are split again as in rung 1 (with rung 1's records); one record per intersection code found (stage `repair`, info, new action `resolve` / `applied`: `rung`, `intersecting_pairs`, `coplanar_pairs`, `faces_before`/`after`, `pieces_dropped`, `components_before`/`after`, `euler_before`/`after`, `rounding_rounds`; the first also `faces_removed`, `faces_added`, `vertices_added` for the accounting); diagnosis's findings of those codes, and non-manifold errors the resolved surface no longer has (a failed split keeps its action), become warnings with `resolve` / `applied` (`RecordLog::repaired`). On failure (an open part, rounding that cannot be made clean) the records have `resolve` / `failed` and `reason`, and the part goes to the fallback as before. `RepairSummary` gains `rung` and `components_before`/`after`, `euler_before`/`after` for every part (the per-code cap can leave a rung-2 record out: on a square frame of four bars diagnosis had already reported ten `INT_SELF_INTERSECTION` pairs). The schema, `docs/report.md` and the defect table are updated.
    - **Fixed on the way:** rung 1 leaves split vertices at one position with different indices; the arrangement merged them as segment ends but not as triangle corners ("a segment end is not a point of its cluster", 8 models). `arrangement` now welds exactly coincident input points first (`coincident_input_points_are_one_point` reproduces it).
    - **Tests:** synthetic cases `bar through a cube`, `boxes overlapping on a common floor`, `boxes touching face to face`, and the existing self-intersection and shell-overlap cases now resolve and mesh conformally (goldens regenerated: those cases from `stuffing` (rejected for the pushed-through icosphere) to `delaunay accepted`, records `resolve/applied`; nothing else changed). `tests/synthetic/tests/resolve.rs` checks the summary's components and Euler characteristic against a brute-force count on 48³ voxels classified by the input's winding number (the surface's Euler characteristic is twice the solid's) for six scenes, among them a square frame of four overlapping bars (one component, Euler characteristic 0) and the pushed-through icosphere. The CLI test of overlapping cubes now exits 10 with the part resolved and conformal, also under `--strict` (rung 2 is a repair, not a fallback); the fallback's handling of surface-gate errors keeps its unit tests in `report.rs`.
    - **On the dirty subset with the CLI** (`run_subset.py tests/thingi10k/subset_dirty.txt`, 8 jobs, against the same run at 25aaaa1): rung 2 resolves 45 of 100 models; conformal 31 to 46 (14 of the resolved ones); exit codes 0: 2, 10: 29 to 44, 20: 55 to 44, 30: 14 to 10; no model gets a worse exit code. Of the resolved models still meshed by stuffing, most pass the gate but fail conformal recovery or remeshing (34 `MESH_RECOVERY_FAILED`; remeshed surfaces that cross themselves), or time out: the resolved surfaces carry sliver pieces along the intersection curves. Open parts (5 of the self-intersecting models here) are refused by rung 2 and go to stuffing as before. About 3 h of CPU per run (8 jobs).
- [x] Exit criterion (met Nov 13, 2026, on 448 of 450 models; the two worker timeouts are deferred by decision with improvements 5 and 6):
  - every part of the self-intersecting and multi-component strata resolves to a surface that passes the gate, or fails with a located record
  - component counts and Euler characteristic match a brute-force check on the synthetic intersection cases
  - **Exit run** (Nov 13, 2026; `run_subset.py` over the 450 models of the 1,000-model sample's two strata, 300 self-intersecting and 150 multi-component, 8 jobs, about 18 h of CPU): 448 of 450 models meet the first part; the second is met by `tests/synthetic/tests/resolve.rs` (stage 5). Of the 448 parts of the 448 models with parts: 379 pass the gate (rung 2 resolves 196, 193 of them passing the gate, 55 meshed conformally, 138 by stuffing; rung 1 alone passes 186); 60 fail it with located `surface_gate` records (40 open parts that rung 2 refuses, with `resolve` / `failed` and the reason, one whose arrangement cannot subdivide a flat triangle, 3 resolved parts with non-manifold edges or an open sheet, and parts rung 1 leaves with non-manifold edges or holes); 9 are not repaired (fatal or empty parts, with located errors). Exit codes: self-intersecting 10: 79, 20: 193, 30: 28; multi-component 0: 12, 10: 85, 20: 50, 30: 3; conformal 177 of 450.
    - **The two that do not:** 466329 (67,888 faces) and 498459 (40,542) time out in the worker (600 s) before any part is reported, so their only record is `MESH_WORKER_FAILED` without a location. On 466329 rung 2 takes 3 s (16 intersecting pairs); the time is elsewhere. On 498459 `resolve_rounded` takes 87 s (`--example arrange`: 76,189 intersecting pairs, 368,368 pieces, 229,997 dropped, rounded clean to 125,916 triangles), so most of the time is also elsewhere (diagnosis of the pairs, meshing the resolved surface). Timeouts and point budgets are deferred by decision (improvements 5 and 6).
    - **Fixed by the run** (and rerun on the 13 models they touched): (1) vertices rung 2 makes have no model vertex, and the copies a later split makes of one were different points of contact, so the gate counted the triangles meeting there as intersecting (9 resolved parts failed the gate; now 6 pass, the other 3 fail on non-manifold edges and an open sheet, located); contact ids of such vertices now go by position. (2) The record cap was per code, so when diagnosis had filled its 100 records of `INT_SELF_INTERSECTION`, rung 2's outcome and the gate's located records of it were dropped (7 parts); the cap is now per code and stage (`log_caps_per_code_and_stage`). The synthetic goldens are unchanged by both.

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
5. *Early stop* (deferred by decision, Nov 13, 2026: look into it in the future). Recovery and refinement should give up when progress relative to the input stalls, rather than running to the budget or the timeout; some models take minutes before falling back.
6. *Point-budget overruns* (deferred with improvement 5, Nov 13, 2026; on the Milestone 15 subset runs 39245, 101548 and 472078 hit the 600 s worker timeout on some runs and fall back to stuffing, and 132060 and 677341 hit the point budget; 57657, 73154; 43853 is fixed by the sizing field; in the Milestone 16 exit run 466329 and 498459 time out before any part is reported, so the timeout record has no location, and 30 models in all end with a worker timeout). Shape refinement near small input features still adds too many points, and the size field does not bound it (see the Milestone 10 findings). Candidates: local feature size from the input's own small features (short edges, small angles) to stop shape refinement there, or protecting small input angles as TetGen does.
7. *Speed.* About 20× TetGen's total time before Milestone 10's profiling; boundary insertions were the hot spot (a cavity shrunk one tet per pass, about 150 passes per insertion on some models). Those passes are now cheap, but their number is not: shrinking several tets at once, or smaller cavities for boundary splits, would change results and needs a benchmark run to judge. Winding-number inside tests are the next largest cost (about a third on 43853); inheriting inside-ness across new tets would remove most of them.
8. *Boundary slivers.* Smoothing of free vertices (Milestone 10) removed most cells below 1° that have a vertex off the boundary. What remains is of two kinds. Tets on input triangles that are themselves needles (113890 has 426 input triangles with an angle under 1°): nothing that keeps the boundary can fix them, but remeshing can. And tets whose four vertices are all on the boundary (thin regions, and most of the remeshed 472020's 0.4°): sliding boundary Steiner points within their facet, or 4–4 flips and multi-face removal, would be next.
9. *Touching shells conformally.* Shells sharing an edge or a vertex fall back to stuffing; meshing them as one constrained region would keep them exact.
10. *The SEACAS check* of the Exodus output (deferred from Milestone 5; do it before Milestone 10).
11. *Tets flat in floating point.* Mostly done (Nov 9, 2026; see the Milestone 10 findings): of the 16 benchmark models whose conformal meshes had them, 13 now mesh conformally with none. Two still have them (50659: 57, 45136: 14). Recovery splits an input edge near a corner into Steiner points about 1e-7 apart and a few ulps off the line, so the facet's triangulation needs needle triangles through three nearly collinear points, and any tet standing on one is flat. Collapsing such a Steiner point into its neighbour along the edge, or not splitting an edge that finely, would remove them; until then those parts use the fallback.
