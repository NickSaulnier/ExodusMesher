# ExodusMesher Roadmap

Nov 8, 2026 · @Nicholas Saulnier

Turns the open-source mesher roadmap (Nov 7, 2026) and the 2022–2026 literature review into milestones for this repository. `PLAN.md` stays the source of truth for the milestone in progress. When a phase here starts, its milestones move into `PLAN.md` with full task lists, and this file keeps only the outline.

## Where this aims

An open mesher cannot match Ansys's breadth. It can beat it on five things, and every milestone below serves at least one:

| Advantage | What it means here |
| --- | --- |
| Robustness you can see | Every failure is a typed record with a location and a hint (the defect taxonomy and `report.json`) |
| Determinism | Same input and config, same bytes out, whatever the core count |
| No per-core licensing | MIT; GPL and AGPL code never in a shipped crate |
| Automation first | Exit codes, a JSON schema, a Python API, diagnostics an agent can act on |
| Openness for research | Published benchmarks, an outcome dataset, every algorithm inspectable |

The review's main conclusion shapes the order of everything. Exact predicates on implicit points have made constrained Delaunay tetrahedrization provably robust on valid input: Diazzi et al. (2023) mesh all 4,408 valid Thingi10k models, and the 2026 surface-chamfering follow-up adds refinement that is guaranteed to converge. Robustness now has two weak points: rounding the exact result to f64 coordinates, and repairing input that is not a valid piecewise-linear complex. Our kernel's 54 of 100 conformal models (`PLAN.md`, Milestone 8) is the gap this closes.

## Where we stand

The roadmap's Phase 1 is the current ExodusMesher plan plus the literature changes. Against its deliverables:

| Phase 1 deliverable | State | Where |
| --- | --- | --- |
| Tolerant OBJ parser, part selection, defect taxonomy, `report.json` | Done | Milestones 1–4, `docs/report.md` |
| Repair ladder, rung 1: local fixes | Done | Milestone 4 |
| Repair ladder, rung 2: exact self-intersection resolution | Not started; self-intersections go to stuffing | Milestone 16 below |
| Repair ladder, rung 3: wrapping | Partly: stuffing on the winding-number field is a lattice wrap, with no feature fidelity | Milestone 17 below |
| Constrained Delaunay kernel with a stuffing fallback | Done, but on f64 Steiner points and plane labels, not indirect predicates; 54% conformal on the 100-model clean subset (Nov 7); on the 1,000-model sample (Nov 8), 304 of 350 clean models and 568 of 1,000 overall | Milestone 8; Milestones 11–13 below |
| Exodus II, CLI, worker, GUI | Done (SEACAS check deferred; GUI interactive check open) | Milestones 5, 9 |
| Generalized winding numbers | Done for inside/outside; ambiguity near 0.5 not reported | Milestone 16 below |
| Validity check after rounding to f64 | Not needed yet (no implicit points); required once the exact kernel lands | Milestone 13 below |
| Manufactured-solution check | Done in-house (`mesher-fvm`); not yet on your own solver | Phase 1 exit |

Phase 2 is partly done already: the sizing field (curvature, thickness, gradation 0.4), feature-preserving isotropic remeshing (40° features, group boundaries), named side sets per OBJ group, golden reports and CI all landed in Milestone 10.

## Decisions needed before Milestone 11

1. **Milestone 10's exit criterion.** The 1,000-model benchmark ran on Nov 8, 2026 (943 meshed, 568 conformal; TetGen meshed 543; table in `PLAN.md`). It is the baseline the new kernel must beat. Only CI green on GitHub remains.
2. **New kernel beside the old one.** Recommended: build the exact kernel as `mesher_tet::cdt` next to `mesher_tet::conforming`, choose it per job by a setting, and make it the default only when it beats the current kernel on the benchmark. Results stay comparable, and the golden reports keep guarding the old path until it retires.
3. **Your solver's Exodus reader.** Whether it reads polyhedral blocks (`NFACED`), and the classic and 64-bit-offset netCDF formats, decides Milestone 21 and Phase 4's output format. Check before Phase 2.
4. **Going public.** Check the employment agreement's terms on outside projects before publishing benchmarks or the dataset (yours to do).

## Ground rules the roadmap adds

- **Clean-room from papers.** Implement from the published descriptions only. Attene's CDT and DelOptim code and the Antipodal code have unverified licences, the chamfering paper is CC BY-NC-ND, and CGAL (alpha wrapping, autorefine) is GPL. None of their code is read, ported or linked.
- **Exact references in tests.** Exact predicates are checked against pure-Rust rational arithmetic (`num-rational`, MIT/Apache), not GMP-based crates (LGPL). The `delaunay` crate (BSD-3) can serve as an extra oracle for the unconstrained case.
- **External tools stay external.** OpenFOAM (GPL), SU2 (LGPL), Gmsh (GPL) and fTetWild run as separate processes in test harnesses, as TetGen does, never linked.
- **New defect codes land with their milestone.** Each code proposed below is added to `DefectCode`, the `PLAN.md` table and `docs/report.md` together when its milestone starts; the names are proposals.
- **Benchmark conventions.** Every published result reports success on all of Thingi10k (fTetWild's convention) and on the valid subset (Diazzi's), a failure taxonomy (algorithmic limit, budget, timeout, out of memory), whether the boundary is exact or within an envelope, and the rate of valid f64 rounding.

## Phase 1 · Finish ExodusMesher

**Goal.** A robust, fully automatic tet mesher for messy OBJ input that never fails silently, with the kernel rebuilt on exact predicates.

**Phase exit criteria.**

- On the valid Thingi10k subset, conformal success comparable to the published exact CDT results (100% of valid models), reported on both the full and valid denominators.
- Every failure carries a defect code and a location; zero crashes, zero silent repairs, no invalid tet ever written.
- Exodus files load in your finite-volume solver and pass a manufactured-solution convergence check there, as well as in `mesher-fvm`.

**Order.** Indirect predicates come first because both the exact kernel (Milestones 12–14) and exact intersection resolution (Milestone 16) are built on them. The kernel comes before the repair ladder because it closes the largest measured gap (46% of clean models fall back today).

### Milestone 10 · Sizing, benchmark and hardening (in progress)

See `PLAN.md`. Open: CI green on GitHub, and the SEACAS check (deferred by you). The benchmark table is in `PLAN.md`.

### Milestone 11 · Indirect predicates

The exact-geometry foundation, in `mesher-geom` (module `implicit`).

- `enum Point { Explicit, Lnc, Bpt, Lpi, Tpi }`: a point on a segment (t·a + (1−t)·b), a barycentric point on a triangle, a line–plane intersection and a three-plane intersection, each stored by its defining explicit points and parameters, never rounded during construction.
- `orient3d`, `insphere` and the coplanar in-circle test for every combination of point kinds, each as a floating-point filter (the constants from Diazzi et al. 2023, Appendix A), then interval arithmetic, then exact expansion arithmetic.
- An approximate f64 position per implicit point, for spatial queries and output only; never used for a topological decision.

*Exit.* Every predicate agrees with an exact rational reference on random, degenerate and nearly degenerate inputs (cospherical, coplanar, points on segments and faces, intersections of nearly parallel planes); the filter settles at least 95% of random queries; property tests at 100,000 cases each.

### Milestone 12 · Exact constrained Delaunay kernel

`mesher_tet::cdt`, following Diazzi et al. (2023), beside the current kernel.

- Input: a valid PLC (the surface gate guarantees it). Every input triangle is a constraint; no facet merging by tolerance, which removes improvement 3 (float32-noise facets) at its root.
- Segment recovery by Si's encroachment-based splitting, with LNC Steiner points.
- Face recovery by splitting the cavity of tets the missing face pierces and refilling each half. Before refilling, check that the two half-cavities have disjoint interiors: the flaw in TetGen's theory that the paper's Section 3.5 found (2 of 4,408 models), which a port of TetGen would inherit. When cavity expansion fails, fall back to gift-wrapping with symbolic perturbation; build it from the start, not as a patch.
- Region marking by the existing parity flood fill; deterministic insertion order (the existing BRIO).
- No refinement yet: the output is exact and may have poor tets.

*Exit.* In exact form (before rounding), the kernel recovers every constraint on the synthetic suite, the 100-model clean subset and the rounded UV sphere (`rounded_uv_sphere_meshes_conformally`, ignored today); a test reproduces a disjoint-cavity case. Runtime and Steiner counts recorded against the paper's (4.3 s average, most runs under 50k Steiner points).

### Milestone 13 · Rounding to f64

Exodus stores doubles, so the exact mesh must round. The review's warning: 93.2% of exact meshes round validly as they are, 99.8% after post-processing, and some inputs (the Schönhardt polyhedron) cannot be represented in floats without extra points. This milestone gets as much care as the kernel.

- Round every implicit point, then check the exact `orient3d` sign of every tet on the rounded coordinates.
- Local untangling of what flipped: vertex moves within the point's rounding box, flips, and re-insertion of extra explicit points where nothing else works, each recorded.
- Boundary conformity re-checked after rounding (the boundary may move by rounding only).
- Proposed codes: `MESH_ROUNDING_REPAIRED` (warning, with the tets fixed and how) and `MESH_ROUNDING_FAILED` (error, triggers the next fallback).
- Pipeline integration: the exact kernel becomes selectable (`--kernel exact`), gated by the existing volume gate, with the current kernel and stuffing behind it.

*Exit.* The rounding validity rate, before and after untangling, measured and reported on the clean subset and the synthetic suite; no tet with non-positive volume is ever written; every untangle is a record.

### Milestone 14 · Chamfer-and-refine

Quality refinement that is guaranteed to finish, following the surface-chamfering paper (SIGGRAPH 2026).

- Chamfer acute vertices and edges within an ε from local feature size; Ruppert-style refinement with BPT points for face splits; restore the chamfered regions with the enriched CDT so the output still conforms exactly.
- Refinement against the existing sizing field.
- Opt-in at first: the paper reports about 80 s per model and one 36-hour outlier, so it runs under the work budget, and `MESH_BUDGET_EXCEEDED` hands over to Milestone 15's optimizer on the unrefined mesh.

*Exit.* On the clean subset, refinement terminates within its budget on every model or records why not; minimum face angle reported against the paper's (1.3° to 3.6°); cost per model tabulated, to settle the review's open question of how much optimization can replace it.

### Milestone 15 · Worst-cell quality for finite volume

Finite-volume accuracy hinges on a handful of skewed cells at the wall, not on averages (Cortellessa, Giacomini and Huerta 2025: errors stagnated near 10⁻¹ under refinement with average non-orthogonality at only 2–5°). Neither exact-CDT paper bounds dihedral angles.

- Invariant-guarded local operations in the Wildmeshing Toolkit's style: flips, splits, collapses and vertex moves, each accepted only if every tet stays positive (exact), the boundary stays exact, and the target metric improves. The existing sliver smoothing becomes one operation among them.
- A WSVM-style energy (a weighted sum of squared tet volumes, convex per vertex, 3×3 Newton steps) as the vertex-move objective, measured against the existing candidate search.
- Aimed at finite-volume metrics: per-face non-orthogonality and |s|/|d| skewness, prioritising cells with a boundary face.
- Report: the worst cells with their locations, separately for wall-adjacent cells; `VOL_LOW_QUALITY` gated on non-orthogonality and skewness as well as dihedral angles; thresholds with headroom for run-to-run variation (Sandia's re-meshing study saw near-surface minimum quality move by about 0.15).
- Absorbs `PLAN.md` improvements 5 (early stop), 6 (point-budget overruns) and 8 (boundary slivers).
- Parallel by independent sets with a fixed colouring, so output does not depend on the core count.

*Exit.* On the clean subset, the worst wall-adjacent non-orthogonality and skewness are lower than both the current kernel's and TetGen's; the report locates the worst cells; results are byte-identical at 1, 4 and 16 threads.

### Milestone 16 · Repair ladder, rung 2: exact arrangement

About 45% of Thingi10k models self-intersect and 29% have several components; these never reach the exact kernel today.

- Mesh arrangement of the whole part's triangle soup with indirect predicates (LPI and TPI points), following Cherchi et al. (2022), with Guo and Fu's indirect offset predicates (2024) as the faster variant to evaluate.
- Cells labelled by generalized winding number; keep the region the winding number says is solid, which resolves `INT_SELF_INTERSECTION`, `INT_SHELL_OVERLAP` and `PART_OVERLAP`, and lets touching shells mesh as one region (improvement 9).
- Winding numbers by the Antipodal method's ray-crossing form (a BVH crossing count plus a boundary term that vanishes on closed meshes), replacing the current evaluation where it is faster.
- The resolved surface carries implicit points; it goes straight to the exact kernel, or is snap-rounded to f64 (as CGAL 6.1's autorefine-and-snap does) and re-checked, since rounding can reintroduce intersections.
- Report: the rung that succeeded per part; components and Euler characteristic before and after; winding-number ambiguity (proposed `TOP_AMBIGUOUS_INSIDE`: sample points with values near 0.5, located, marking holes and flipped patches).

*Exit.* On the self-intersecting and multi-component strata of the 1,000-model sample, every part either resolves to a surface that passes the gate or fails with a located record; component counts and Euler characteristic match a brute-force check on the synthetic intersection cases.

### Milestone 17 · Repair ladder, rung 3: wrapping

- Alpha wrapping from Portaneri et al. (2022), clean-room: guaranteed to terminate, output oriented, 2-manifold, intersection-free and strictly enclosing the input.
- Alpha and offset from the element size and the part's diagonal; both, and the volume change, in the report (proposed `REPAIR_WRAPPED`, a warning naming lost thin features where the deviation shows them).
- The wrapped surface goes through the exact kernel, so dirty input keeps sharp-ish features; stuffing stays as the last resort behind it.

*Exit, and Phase 1's.* The full benchmark with the conventions above, on all 10,000 models and the valid subset; outputs checked by `mesher-fvm` and by your solver.

**Research hooks.** Guaranteed-valid f64 rounding (Milestone 13); a public dataset of meshing outcomes from the reports (see Hardening).

**Main risk.** Milestones 11–13 are the hardest code in the project. The current kernel and stuffing keep the tool useful while they mature, and nothing replaces them until the benchmark says so.

## Phase 2 · Surfaces and formats

**Goal.** Meshes that respect sharp features and named boundaries, size themselves sensibly, and load in the major open solvers.

**Phase exit criteria.** Sharp edges and named patches survive remeshing on the synthetic suite and a set of public CAD-derived OBJ files; neighbouring element sizes stay within the gradation limit; output passes OpenFOAM's `checkMesh` with no errors.

### Milestone 18 · Features, patches and size controls

- Feature curves as first-class output: dihedral threshold and OBJ group boundaries (both exist in `mesher_tet::remesh`), every detected feature reported, thresholds and per-edge overrides in the project file.
- Patches carried through every path, including stuffing (today a stuffed part gets only `<part>_boundary`).
- User size controls per part and per patch in the project file (closes `PLAN.md`'s open question on size controls); proximity between separate surfaces added to the thickness term.
- A gradation check in the volume gate: neighbouring cell size ratio against the limit.

*Exit.* As the phase's first two criteria, on the synthetic suite plus the public CAD-derived set.

### Milestone 19 · Solver formats

- OpenFOAM polyMesh and Fluent `.msh` writers (both text), side sets as patches.
- `.meshb` (libMeshb) read and write, the format of the UGAWG adaptation benchmarks.
- CGNS later: it needs HDF5, and pure-Rust HDF5 writing is immature.
- `checkMesh` in a Linux CI job on the synthetic outputs (OpenFOAM installed in the runner, run as a process).

*Exit.* `checkMesh` passes on every synthetic output; a Fluent `.msh` and a `.meshb` file round-trip through independent readers.

### Milestone 20 · Python API

- `mesher-py` with PyO3 and maturin: run a job, get the report as typed Python objects, list parts and records, read back the mesh.
- Jobs still run in `mesher-worker` through `mesher-client`, so a kernel panic or a runaway job cannot take down the caller's interpreter.
- Agent-facing diagnostics: every record already has a code, location and hint; add machine-readable suggested remedies (the setting or repair to try next), since LLM simulation agents fail mostly at meshing and work through generate-validate-repair loops.

*Exit.* A Python script meshes the synthetic suite and reproduces the CLI's reports byte for byte; wheels build in CI on Linux and Windows.

### Milestone 21 · Experimental polyhedral dual

- A vertex-centred dual of the tet mesh (about one cell per tet vertex, against roughly five tets per vertex today), with boundary cells clipped to the surface.
- `mesher-fvm` extended to polyhedral cells.
- Polyhedral output, if Decision 3 confirms your solver reads it.
- Generated after any adaptation, never before.

*Exit and research hook.* A controlled tet-versus-dual comparison on the same geometries and solver: error per degree of freedom and iterations to converge. Current evidence rests on one weak 2021 study, so the result is publishable either way.

**Main risk.** False sharp edges on noisy triangle soups; expose the thresholds and report every detected feature.

## Phase 3 · Boundary layers

**Goal.** Prism layers that resolve viscous walls automatically and stay valid in concave corners and narrow gaps, where most meshers give up or silently thin the layers.

**Phase exit criteria.** Full layers, or explicitly reported reductions, on concave-corner and narrow-gap test cases; on a public aircraft model, solver-reported wall y+ within the requested target.

### Milestone 22 · Advancing-layer baseline

- Layer controls per patch: first-cell height or a target y+, growth ratio, layer count.
- Advancing layers along smoothed normals, with several normals at convex corners and height compression in concave and narrow regions.
- Collision detection, with the "air mesh" idea from Ye and Chen's group used locally: a thin shell of tets ahead of the front turns a collision into an inverted element the advance refuses.
- Every reduction a record with its location (proposed `BL_LAYERS_REDUCED`, `BL_COLLISION`).
- Synthetic cases: a concave corner, a narrow gap, a U-shape.

### Milestone 23 · Prisms joined to the core

- `WEDGE6` blocks in Exodus, `.msh`, `.vtu` and polyMesh; side sets on the walls and the layer top.
- The layer top becomes the constrained surface for the tet core, so the join is conformal by construction.
- Boundary-layer quality: non-orthogonality at the layer top, the volume jump into the core, inverted prisms (proposed `BL_INVERTED_PRISM`).

*Exit.* As the phase criteria; the y+ check needs a solver run (your solver, OpenFOAM or SU2).

### Milestone 24 · Research track: fast guaranteed-valid layers

The mapping-based method (IMR 2023, JCP 2025) keeps full layers where others stop, but took about 15 hours for one layer on DLR-F6. Making it fast is the roadmap's most valuable research contribution: the Milestone 22 result as the initial guess, GPU optimisation, or a learned initialisation (learn the decision, keep the guarantee).

**Main risk.** Effort. Ship the reported-reduction baseline first; treat Milestone 24 as research.

## Phase 4 · Hex core

**Goal.** Fewer, better-aligned cells away from walls: an octree hex core joined to the layers through polyhedral transition cells, the open equivalent of poly-hexcore.

**Phase exit criteria.** Materially fewer cells than all-tet meshes at equal wall resolution; solver accuracy and iterations to converge equal or better.

### Milestone 25 · Octree core

- An octree refined by the size field with 2:1 balance; the graded stuffing octree (`mesher_tet::graded`, already weakly balanced) is the starting point.
- The core trimmed back from the layer top (or the surface) by a buffer zone.
- `HEX8` blocks.

### Milestone 26 · Transition cells and mixed output

- The buffer filled with tets from the exact kernel, then merged into polyhedra that meet the hex faces conformally; hanging nodes at level changes as polyhedra or templates.
- Mixed-element output: `HEX8`, `WEDGE6`, `TETRA4` and polyhedral blocks (Exodus `NFACED`, subject to Decision 3).
- Worst-cell reporting (Milestone 15) extended to transition cells, where poly-hexcore meshes hide their worst skewness.

*Exit and research hook.* The tet versus polyhedral dual versus hex-core comparison on the same geometries and solver; the review found no solid 2022–2026 study, so a reproducible one has value on its own.

## Phase 5 · CAD and scale

**Goal.** Take engineering CAD directly, and mesh aircraft- and vehicle-scale models on HPC systems.

**Phase exit criteria.** A public aircraft model meshed from STEP with named patches intact; a 100M+ cell mesh on one HPC node, with measured scaling and output identical at every core count.

### Milestone 27 · STEP import

- `mesher-cad` over OpenCASCADE (LGPL), behind a feature flag and dynamically linked. It breaks the pure-Rust rule, so that rule gains a written exception for this crate first.
- CAD faces become named patches and CAD edges feature curves, so Phase 2 and 3 controls apply per face.
- Defeaturing (small holes, fillets, slivers) as reported, reversible operations, never silent.

### Milestone 28 · Parallel and distributed

- Parallel local operations in the Wildmeshing Toolkit's model, deterministic regardless of thread count (Milestone 15's colouring generalised).
- Domain decomposition with MPI for distributed meshing.
- GPU passes for winding numbers and quality (the Antipodal method reports about 10⁹ queries a second on GPU).
- Large output: `netcdf3` writes whole variables only and the 64-bit-offset format caps each variable near 4 GiB (`EXO_CAPACITY_EXCEEDED`), so CDF-5, netCDF-4 or a streaming writer of our own.

**Research hook.** Learned defeaturing: predicting which CAD features a simulation can drop (Sandia's random forests on B-rep features already work for related decisions).

**Main risk.** CAD import is a long tail of kernel quirks; budget it as ongoing maintenance.

## Hardening, throughout

What separates a research mesher from a commercial one is mostly this.

- **Corpus.** Thingi10k in full (as a local run; meshes never committed), then ABC, then public aerospace and automotive geometries (NASA CRM, DLR-F6, DrivAer), each after checking its licence. Success rate, worst-cell quality and runtime tracked per run, with the benchmark conventions above.
- **Regression gate.** The synthetic golden reports in CI on every change (done); a small public-domain subset in CI once one exists that can be committed.
- **Fuzzing and property tests.** `cargo-fuzz` on the OBJ parser, the predicates (Milestone 11) and the repair ladder, in a scheduled Linux CI job; the existing property tests stay.
- **Solver in the loop.** `mesher-fvm` on every output now; OpenFOAM or SU2 after Milestone 19; AIAA workshop cases once boundary layers exist. A good-looking mesh is not enough.
- **Head-to-head benchmarks.** TetGen (done), then Gmsh and fTetWild, run as processes. Commercial tools only where their licence terms allow publishing results.
- **Outcome dataset.** The reports over Thingi10k, published as a dataset pairing geometry ids with meshing outcomes, failure codes and locations. The IMR survey notes that no public corpus does this; it supports failure-prediction models (`PLAN.md`'s "after v1" idea) and agent research.
- **Venue.** The SIAM IMR meshing contest and its geometry-repair and AI-for-meshing tracks.

## People, licence and funding

Not milestones; recorded so the decisions stay visible.

- **Effort.** The roadmap's ranges assume a team of three to five: Phase 1 12–18 months, Phase 2 6–9, Phases 3 and 4 12–18 each, Phase 5 18–24. Phases 2 and 3 can overlap once the Phase 1 kernel is stable.
- **Licence.** MIT already. Keep GPL and AGPL out of shipped crates (the schema test already enforces it for TetGen); OpenCASCADE's LGPL stays behind its feature flag.
- **Funding and governance.** SBIR/STTR, DOE and NASA programmes (the CFD Vision 2030 study names meshing as a bottleneck), an industry consortium or a foundation, on the OpenFOAM, Gmsh and SU2 models. Pairing with an open solver community brings users and validation data.
- **Separation from work.** No employer code, geometry or knowledge gained under work licences; public data only (already a project rule).

## References

From the literature review (Nov 2026); entries after June 2026 are as the review gives them.

- Diazzi, Panozzo, Vaxman, Attene. Constrained Delaunay tetrahedrization: a robust and practical approach. 2023. https://cims.nyu.edu/gcl/papers/2023-CDT.pdf
- Surface chamfering for robust tetrahedral meshing. ACM TOG 45(4), SIGGRAPH 2026. https://cims.nyu.edu/gcl/papers/2026-chamfering.pdf
- WSVM tetrahedral mesh optimisation. 2024. https://arxiv.org/html/2409.05525
- Wildmeshing Toolkit (MIT). https://github.com/wildmeshing/wildmeshing-toolkit
- Zhou, Jacobson. Thingi10K. https://arxiv.org/pdf/1605.04797
- Hu et al. fTetWild. https://ar5iv.arxiv.org/html/1908.03581
- Cherchi, Pellacini, Attene, Livesu. Interactive and robust mesh booleans. TOG 2022. https://arxiv.org/pdf/2205.14151.pdf
- Guo, Fu. Indirect offset predicates for mesh arrangements. TOG 2024. https://mangoleaves.github.io/projects/mesh-arrangements
- The Antipodal method for generalized winding numbers. 2026. https://arxiv.org/abs/2605.01536
- Portaneri et al. Alpha wrapping with an offset. SIGGRAPH 2022. https://doc.cgal.org/latest/Alpha_wrap_3/
- Cortellessa, Giacomini, Huerta. Finite-volume accuracy on distorted simplicial meshes. 2025. https://arxiv.org/pdf/2501.00450
- Owen et al. AI and machine learning for meshing (survey). https://arxiv.org/html/2512.23719v2
- Ye et al. Mapping-based boundary-layer prisms. IMR 2023, JCP 2025. https://internationalmeshingroundtable.com/assets/papers/2023/06-Ye-compressed.pdf
- Foam-Agent 2.0. https://arxiv.org/pdf/2509.18178
- UGAWG adaptation benchmarks. https://ugawg.github.io/
