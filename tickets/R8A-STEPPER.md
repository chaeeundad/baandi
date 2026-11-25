TASK_ID: R8A-STEPPER
GOAL: |
  Stage R8a of docs/roadmap/refactor-r-plan.md (user decision D1: port runtime control): split the reference SCF loops into a stepper so checkpoint/restart, observe/agent runtime, MCP and Python stepping can drive them (R8b-R8d). This is the highest-risk refactor ticket: it touches the hot loop of the fast path.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 2 rows L6-L8, 4 "R8", 5 risk 1); docs/performance/AGENT-BRIEF.md (hot-path rules: no allocation, logging or JSON in the loop); crates/dft-engine/src/refscf/{driver.rs,k_driver.rs,electronic.rs,api.rs,record.rs,kept_hx.rs,map.rs}; the legacy ScfSession step contract in crates/dft-engine/src/session.rs (what a step returns, safe points) and runtime.rs, for the interface shape only.
TASK: |
  1. refscf/stepper.rs: Stepper { initialize(), step() -> StepReport, is_done(), finish() -> delivered result } for the Gamma (driver.rs), k (k_driver.rs) and spin/FD (electronic.rs) drivers. One step = one SCF map including its mixing decision; the verification map(s) are steps too, flagged as such. StepReport carries iteration index, estimator value, dr2, energy if computed, ethr, H columns, kind (ordinary / verification / retry), and whether the state is at a safe point (after mixing, before the next eigensolve).
  2. run() of each driver becomes Stepper::run_to_end() with exactly the same arithmetic order, collectives and allocations: no new allocation per map (preallocate in initialize), no extra clones of orbitals/HX/density. api::solve uses run_to_end.
  3. Identity: reference_identity, bench_scf --require-identical (baseline-reference.json), and a stepwise test that drives step() to completion and compares with run() bit-for-bit (energy bits, density hash, orbitals hash) on Gamma, k222s, FD and fixed spin fixtures, ranks 1/3/8.
  4. Write docs/performance/lanes/RSCF/r8a-dev08-request.md: ONE run per arm, si128 MPI8 verify off and on, main vs R8a: wall, maps, H columns; identity of energy bits and density hash.
OWNED_FILES:
  - crates/dft-engine/src/refscf/stepper.rs   # new
  - crates/dft-engine/src/refscf/driver.rs
  - crates/dft-engine/src/refscf/k_driver.rs
  - crates/dft-engine/src/refscf/electronic.rs
  - crates/dft-engine/src/refscf/api.rs
  - crates/dft-engine/src/refscf/record.rs
  - crates/dft-engine/src/refscf/mod.rs   # declaration only
  - crates/dft-engine/tests/refscf_stepper.rs   # new
  - docs/design/reference-scf.md   # a "Stepper" section only
  - docs/performance/lanes/RSCF/r8a-dev08-request.md   # new
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf_stepper
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli --test reference_identity
  - scripts/check.sh (BAANDI_TEST_TMP outside the checkout, e.g. /private/tmp/r8a-check)
  - python3 -I scripts/bench_scf.py --require-identical (release, default baseline)
  - cargo fmt --check && cargo clippy -p dft-engine --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: stepper in place; run() and stepwise runs bit-identical to main; no per-map allocation added; dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
