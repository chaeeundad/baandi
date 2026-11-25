TASK_ID: R7B-CELL-RELAX
GOAL: |
  Stage R7b of docs/roadmap/refactor-r-plan.md: variable-cell relaxation on the reference path. R7a (merged) gives frozen-basis strained cells in the reference start and delivered stress; R3 (merged) gives refscf::api::solve.
READ FIRST: docs/roadmap/refactor-r-plan.md (section 2 rows L4, L5; section 4 "R7"; section 7 decision D5); crates/dft-cli/src/relax.rs (cell-mode rejection ~:761, extrapolation ~:767) and relax/cell.rs; crates/dft-cli/src/driver_select.rs (L4 predicate); crates/dft-engine/src/refscf/{api.rs,delivered.rs}; tests/reference/relaxation/cell/validation.json; docs/performance/attempts-ledger.md rows on RELAX_EXTRAPOLATION / C3 predictor.
TASK: |
  1. Cell relaxation modes in relax.rs/relax/cell.rs call api::solve with the strained cell and an ionic+cell seed: density and orbitals transferred between cells with the frozen Miller map (reuse ionic_seed machinery). Remove the cell-mode rejection and the L4 predicate in driver_select.rs.
  2. D5: read the ledger rows for RELAX_EXTRAPOLATION/C3; if it was measured as beneficial, port it into the seed construction; otherwise retire the switch and the start_extrapolation key (input error as an unknown key; no migration message, user decision D4). State which and why.
  3. Tests: the stored cell-relax protocol cases (tests/reference/relaxation/cell/validation.json) reference vs legacy pin: final lattice within 1e-5 bohr, pressure within 1e-6 Ha/bohr^3, per-step dE within 5e-9 Ha at tightened thresholds; ion-only relax results unchanged.
OWNED_FILES:
  - crates/dft-cli/src/relax.rs
  - crates/dft-cli/src/relax/cell.rs
  - crates/dft-cli/src/driver_select.rs   # L4/L5 predicates only
  - crates/dft-cli/tests/all/cell_relaxation.rs
  - crates/dft-cli/tests/fixtures/routing/   # L4/L5 expectations only
  - crates/dft-engine/src/refscf/delivered.rs   # cell seed only
  - crates/dft-numerics/src/switches.rs   # RELAX_EXTRAPOLATION entry only, if retired
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Forces, stress and relaxation" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-cli relax
  - cargo test -p dft-cli --test all cell_relaxation
  - cargo test -p dft-cli --test driver_routing
  - cargo fmt --check && cargo clippy -p dft-cli -p dft-engine --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: variable-cell relaxation on the reference path agrees with legacy on the stored protocol; D5 resolved with reasons.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
