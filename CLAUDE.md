# ExodusMesher

Pure-Rust automatic surface + volume mesher for .obj files, with Exodus II output for a finite-volume solver, structured error reporting, and an egui GUI for selecting parts. The full design lives in `PLAN.md`; treat it as the source of truth and update it whenever a decision changes.

## How to work

- Build milestone by milestone (`PLAN.md` > Milestones). Do not start the next milestone until the current one's exit criterion passes. Track progress in the task checklists under `PLAN.md` > Phases and task breakdown, ticking items as they land. Milestones 1–4 are done. Milestone 5 (Writers) is done and checked in ParaView; only its SEACAS check is open, deferred by decision (do it before Milestone 10). Milestones 6 (Stuffing) and 7 (Delaunay) are done. Milestone 8 (Recovery and refinement, a constrained Delaunay kernel) has landed with its exit criterion deferred by decision: 54 of 100 Thingi10k models mesh conformally (see `PLAN.md` > Thingi10k benchmark, and `PLAN.md` > Things to consider improving for the open work). Milestone 9 (GUI, `crates/mesher-gui`) passes its exit criterion headless; the user's interactive check is open. Run the GUI with `cargo run -p mesher-gui --release -- model.obj` after building `mesher-worker`. Milestone 10 (Sizing, benchmark and hardening) is in progress; its checklist and decisions are in `PLAN.md`. The 1,000-model benchmark ran on Nov 8, 2026 (table in `PLAN.md` and `README.md`); its exit criterion now waits only on CI green. The SEACAS check and long benchmark runs are deferred by the user: don't start them unasked. The benchmark is `tests/oracle/tests/thingi10k.rs` (ignored; needs the meshes in `target/thingi10k/raw/` and release binaries). Milestones 11 to 17 (Phase 7: the exact kernel and repair ladder) are in `PLAN.md`. Milestone 11 (indirect predicates, `mesher_geom::implicit`) is done; its oracle is `crates/mesher-geom/tests/implicit_oracle.rs` (100,000 cases a test in release, 2,000 in debug, `IMPLICIT_CASES` overrides). Milestone 12 (the exact constrained Delaunay kernel, `mesher_tet::cdt`, beside the current kernel) is in progress: the exact Delaunay core and segment recovery by flips, edge removal and protected splitting are in (97 of the 100-model subset; `crates/mesher-tet/tests/cdt.rs`; measure with `cargo run --release -p mesher-tet --example cdt_segments -- <obj files>`, results in `PLAN.md`); face recovery is next. Later milestones are outlined in `ROADMAP.md` (phases, exit criteria, clean-room and licence rules); move a phase's milestones into `PLAN.md` with full task lists when it starts.
- The synthetic defect suite lives in `tests/synthetic/` (workspace crate `mesher-synthetic`); run `cargo test -p mesher-synthetic -- --nocapture` to see precision and recall per code (set `SUITE_TIMING=1` for per-case times). Any new diagnosis or repair behaviour needs a case there.
- TetGen oracle tests live in `tests/oracle/` (crate `mesher-oracle`, needs MSVC); they write sample files to `target/tmp/oracle/`, including `quality_vs_tetgen.txt`.
- Integration tests in `crates/*/tests/` allow `unwrap`/`expect`/`panic` at file level; library and binary code must not use them.
- After changing report types, regenerate the committed schema: `UPDATE_SCHEMA=1 cargo test -p mesher-core --test schema`, and update the report reference `docs/report.md` (a test checks its defect code table against `DefectCode`).
- The synthetic suite compares every run with `tests/synthetic/golden.txt` (records by code, severity and action; each part's kernel and acceptance exactly; tet count, minimum dihedral and non-orthogonality within bands). After an intended change, regenerate it with `UPDATE_GOLDEN=1 cargo test -p mesher-synthetic` and review the diff.
- Profile the constrained kernel with `cargo run --release -p mesher-tet --features profile --example profile -- target/thingi10k/obj/<id>.obj [--max-points N]`: timing scopes (`mesher_tet::prof::scope`) compile to nothing without the feature. Performance changes must leave results identical; check the golden reports and a model's tet count.
- CI (`.github/workflows/ci.yml`) runs fmt, clippy and `cargo test --workspace --release` on Linux and Windows. The Thingi10k benchmarks stay local: `python tests/thingi10k/fetch.py tests/thingi10k/sample1000.txt`, then the ignored `stratified_sample_benchmark` test (or `tests/thingi10k/run_subset.py` for a quick CLI-only check).
- Follow the workspace layout in `PLAN.md` (`crates/mesher-*`, `tests/`).
- Before calling a task done: `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace`.

## Non-negotiables

- Every failure path maps to a defect code from the taxonomy in `PLAN.md`. Never repair silently. When adding a code, update the report enum and the `PLAN.md` table together.
- Library code never panics on bad input; return errors that become report records. Panics are bugs (exit code 70).
- All orientation and in-sphere decisions go through exact predicates (`robust` crate); no epsilon comparisons for topology.
- Tolerances are relative to the bounding-box diagonal, never absolute.
- Deterministic output: no hash-map iteration order in anything written to disk; seeded RNGs only; the same mesh on every platform, so transcendental functions (`sin`, `acos`, `atan2`, `cbrt`, `log2`, `pow`, ...) go through `libm`, never the platform's `f64` methods (`clippy.toml` enforces it outside the GUI).
- Meshing runs in `mesher-worker`. The CLI and GUI contain no meshing logic; both run the worker through `mesher-client`.
- TetGen (`tritet`) is a dev-dependency for oracle tests only; never a dependency of a shipped crate.
- Exodus: 1-based connectivity, Exodus tet face ordinals for side sets, `TETRA4` blocks per part. Unit-test ordinals on a single tet.
- Quality metrics target finite-volume use: non-orthogonality, skewness, neighbour volume ratio, plus minimum dihedral angle.

## Environment

- Windows dev machine; keep paths and line endings cross-platform. `tritet` builds C/C++ code, so oracle tests need the MSVC toolchain.
- Personal project: public data only (Thingi10k, synthetic shapes). No work or ITAR-controlled geometry.
