# ExodusMesher

Pure-Rust automatic surface + volume mesher for .obj files, with Exodus II output for a finite-volume solver, structured error reporting, and an egui GUI for selecting parts. The full design lives in `PLAN.md`; treat it as the source of truth and update it whenever a decision changes.

## How to work

- Build milestone by milestone (`PLAN.md` > Milestones). Do not start the next milestone until the current one's exit criterion passes. Track progress in the task checklists under `PLAN.md` > Phases and task breakdown, ticking items as they land. Milestones 1 and 2 are done; next is Milestone 3: Diagnosis and parts.
- Integration tests in `crates/*/tests/` allow `unwrap`/`expect`/`panic` at file level; library and binary code must not use them.
- After changing report types, regenerate the committed schema: `UPDATE_SCHEMA=1 cargo test -p mesher-core --test schema`.
- Follow the workspace layout in `PLAN.md` (`crates/mesher-*`, `tests/`).
- Before calling a task done: `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace`.

## Non-negotiables

- Every failure path maps to a defect code from the taxonomy in `PLAN.md`. Never repair silently. When adding a code, update the report enum and the `PLAN.md` table together.
- Library code never panics on bad input; return errors that become report records. Panics are bugs (exit code 70).
- All orientation and in-sphere decisions go through exact predicates (`robust` crate); no epsilon comparisons for topology.
- Tolerances are relative to the bounding-box diagonal, never absolute.
- Deterministic output: no hash-map iteration order in anything written to disk; seeded RNGs only.
- Meshing runs in `mesher-worker`. The CLI and GUI contain no meshing logic.
- TetGen (`tritet`) is a dev-dependency for oracle tests only; never a dependency of a shipped crate.
- Exodus: 1-based connectivity, Exodus tet face ordinals for side sets, `TETRA4` blocks per part. Unit-test ordinals on a single tet.
- Quality metrics target finite-volume use: non-orthogonality, skewness, neighbour volume ratio, plus minimum dihedral angle.

## Environment

- Windows dev machine; keep paths and line endings cross-platform. `tritet` builds C/C++ code, so oracle tests need the MSVC toolchain.
- Personal project: public data only (Thingi10k, synthetic shapes). No work or ITAR-controlled geometry.
