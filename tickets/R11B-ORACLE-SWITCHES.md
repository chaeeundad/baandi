TASK_ID: R11B-ORACLE-SWITCHES
GOAL: |
  Stage R11b of docs/roadmap/refactor-r-plan.md (criterion b, user decision D3): the strict oracle in crates/dft-engine/src/oracle is frozen at its default arithmetic. Remove every switch that is read only by oracle code (or by legacy solvers/physics used only by the oracle) by hard-coding the default arm and deleting the non-default arms. A tunable oracle is not an oracle; only REFERENCE_STRICT selects it.
READ FIRST: docs/roadmap/refactor-r-plan.md (section 4 "R11", classification table rows "read only in legacy engine files" and "legacy solver and physics code"); docs/performance/switch-retirement.md (R11a table: add rows); crates/dft-numerics/src/switches.rs; git grep "Id::NAME" for each candidate; crates/dft-engine/AGENTS.md "Frozen legacy oracle".
TASK: |
  1. For each switch: list readers; if all readers are under crates/dft-engine/src/oracle or legacy-only code in dft-solvers/dft-physics/dft-numerics used only by the oracle, hard-code the default and delete the switch and the non-default branches (and code only reachable from them). Keep switches also read by the reference path or shared kernels.
  2. Delete legacy-only solver/physics code that becomes unreachable (legacy Davidson variants, lean/distributed legacy paths) only if no non-oracle code and no oracle default path uses it.
  3. Extend docs/performance/switch-retirement.md and the ledger (one row per removed switch).
  4. Identity: RG-ORC locally (oracle_identity, bench_scf --require-identical with BAANDI_REFERENCE_STRICT=on against baseline-oracle.json) and the reference baseline; report before/after switches, lines, tests.
OWNED_FILES:
  - crates/dft-engine/src/oracle/
  - crates/dft-solvers/src/   # legacy-only paths
  - crates/dft-physics/src/   # legacy-only paths
  - crates/dft-numerics/src/
  - crates/dft-numerics/tests/
  - crates/dft-engine/tests/   # tests of deleted branches only
  - crates/dft-solvers/tests/
  - crates/dft-physics/tests/
  - crates/dft-cli/tests/   # switch counts/pins only
  - docs/performance/switch-retirement.md
  - docs/performance/attempts-ledger.md   # appended rows
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo check -p dft-engine --no-default-features
  - cargo check -p dft-solvers -p dft-numerics -p dft-mkl --tests --features dft-numerics/mkl
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (both baselines)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: no oracle-only switches remain; oracle and reference bit-identical; counts reported.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
