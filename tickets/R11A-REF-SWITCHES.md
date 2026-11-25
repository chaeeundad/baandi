TASK_ID: R11A-REF-SWITCHES
GOAL: |
  Stage R11 of docs/roadmap/refactor-r-plan.md for the reference path, with user decision D3 (section 7): remove the default-on rollback switches of the reference path together with their pre-adoption code ("they can be recreated quickly if needed"), and remove default-off experiment switches whose ledger status is rejected/superseded, with their code. Keep: diagnostics/profiling switches, infrastructure, REFERENCE_STRICT and REFERENCE_SCF (R10c handles routing), XC_SCREEN_NEGATIVE (ADR 0046 off-switch kept for the oracle comparison), isolated algorithm modules protected by the Generality rule (REFSCF_RMM*, REFSCF_LADDER*, Kerker/TF preconditioning if present, MKL_* backend selectors), and anything the strict oracle reads (R11b handles those after R9c).
READ FIRST: docs/roadmap/refactor-r-plan.md (section 4 "R11" criteria a-e and the classification table, section 7 D3); crates/dft-numerics/src/switches.rs; docs/performance/attempts-ledger.md; docs/performance/AGENT-BRIEF.md (hot-path rules); every reader of each candidate switch (git grep "Id::NAME").
TASK: |
  1. Build a table (docs/performance/switch-retirement.md, new): every switch read by the reference path or shared kernels, with class (rollback default-on / experiment default-off / diagnostic / infra / algorithm module / oracle), its readers, its ledger row, and the decision (remove / keep + reason).
  2. Remove: for each rollback switch, delete the switch and the non-default code branch, leaving the default arithmetic exactly as is; for each rejected default-off experiment, delete the switch and its code. Add one ledger row per removed switch ("retired in R11a, recreate from git <commit> if needed").
  3. Tests that pinned removed switches: delete the tests of deleted branches, keep identity tests by removing the now-meaningless pins; report every deleted test.
  4. Identity: the default reference path must be bit-identical: reference_identity, bench_scf --require-identical against baseline-reference.json (all cases), and the oracle baseline unchanged; also run the existing hot-path identity tests (ref_davidson, kept HX, fused paths). Report before/after counts: switches, lines, tests.
  Do not touch files R9C-ORACLE-FEATURE is moving (legacy driver files: scf.rs, session.rs, convergence.rs, collinear.rs, finite_temperature.rs, gamma_real.rs, m2.rs, mix_dist.rs, cutoff_ladder.rs, band_density.rs, ...); if a candidate switch is read there too, leave it for R11b.
OWNED_FILES:
  - crates/dft-numerics/src/   # switches.rs and kernels whose rollback branches are removed
  - crates/dft-numerics/tests/
  - crates/dft-solvers/src/   # reference Davidson and shared kernels only (not legacy solvers used only by the oracle)
  - crates/dft-solvers/tests/
  - crates/dft-physics/src/   # distributed kernels' rollback branches only
  - crates/dft-physics/tests/
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/src/refstart.rs
  - crates/dft-engine/tests/refscf*.rs
  - crates/dft-mkl/src/   # only if a removed switch selected an adapter path
  - crates/dft-cli/src/run.rs   # switch report fields only
  - crates/dft-cli/tests/   # switch counts/pins only
  - docs/performance/switch-retirement.md   # new
  - docs/performance/attempts-ledger.md   # appended rows
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo check -p dft-solvers -p dft-numerics -p dft-mkl --tests --features dft-numerics/mkl
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (reference baseline; and BAANDI_REFERENCE_STRICT=on with baseline-oracle.json)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: retirement table; rollback and rejected switches removed with their code; default path bit-identical; counts reported.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
