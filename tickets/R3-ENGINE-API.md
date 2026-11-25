TASK_ID: R3-ENGINE-API
GOAL: |
  Stage R3 of docs/roadmap/refactor-r-plan.md: one engine facade refscf::api::solve for the reference path, used by the CLI run and relax reference branches (later by Python, MCP, bands).
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 0, 1a, 3, 4 "R3"); crates/dft-cli/src/run.rs reference branches (k_driver::run/run_prepared, electronic::run/run_prepared, delivered::force_state[_prepared], symmetry preparation) and relax.rs reference branch (delivered::run, force_state, ionic_seed); crates/dft-engine/src/refscf/{mod.rs,driver.rs,k_driver.rs,electronic.rs,delivered.rs,symmetry.rs}.
TASK: |
  1. crates/dft-engine/src/refscf/api.rs: solve(input: StartInput, request: Request) -> Result<Delivered, ScfError>, where Request holds profile, points, Option<symmetry::Prepared>, Option<electronic::Policy>, Option<IonicSeed>, and which derivatives are wanted (forces, stress); it dispatches exactly as run.rs and relax.rs do today (same calls, same order, same collectives), returning one Delivered struct with the run, records and optional force/stress results.
  2. Replace the CLI run and relax reference branches with calls to api::solve. No behaviour change: results.json byte-identical except wall-clock fields; add a JSON-diff helper test with timing keys masked covering Gamma, k222s, FD, fixed spin, symmetry and forces/stress fixtures.
  3. Bit-identity: reference_identity, bench_scf --require-identical against baseline-reference.json (default, no switches).
OWNED_FILES:
  - crates/dft-engine/src/refscf/api.rs   # new
  - crates/dft-engine/src/refscf/mod.rs   # declaration and re-exports only
  - crates/dft-cli/src/run.rs   # reference branches only
  - crates/dft-cli/src/relax.rs   # reference branch only
  - crates/dft-cli/tests/api_identity.rs   # new
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-cli
  - cargo test -p dft-engine --test refscf
  - scripts/check.sh (BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (release, default baseline)
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: facade in place, CLI uses it, outputs byte-identical modulo timing, baselines identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
