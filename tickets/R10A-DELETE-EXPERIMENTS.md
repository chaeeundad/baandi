TASK_ID: R10A-DELETE-EXPERIMENTS
GOAL: |
  Stage R10a of docs/roadmap/refactor-r-plan.md, pulled ahead of R9c (user decision D3: rejected experiment switches and dead code are removed): delete the legacy-only experimental SCF paths and their tests and switches. They are not used by the reference path (R9a moved the prolongation helpers to refscf/prolong.rs) nor by the strict oracle's default arithmetic.
READ FIRST: docs/roadmap/refactor-r-plan.md (section 1c group "Legacy standard-path experiments", section 2 row L12, section 3 line counts, section 4 "R10", section 7 decisions); crates/dft-engine/src/lib.rs and scf.rs/session.rs (where these modules are wired); crates/dft-numerics/src/switches.rs; docs/performance/attempts-ledger.md (for each switch you delete, its row).
TASK: |
  1. Delete: cutoff_ladder.rs (everything except what refscf still uses via re-exports: keep only the re-export shims or redirect callers to refscf::prolong), direct_min.rs, advance.rs, fixed_map.rs, spectral_scf.rs, replay_capture.rs, first_block.rs, verification_profile.rs, floor_tests.rs, and their wiring in scf.rs/session.rs/lib.rs (the strict oracle path must keep its exact default arithmetic: only remove branches gated by the deleted switches, at their non-default values).
  2. Delete their tests (engine tests: cutoff_ladder, direct_min, scf_advance, first_block_diag, and snapshot/fp32/instrument only where they test deleted modules) and the CLI ladder JSON in run.rs. Never delete files under tests/reference (stored evidence stays even when its generator test is deleted).
  3. Remove the switches read only by the deleted code: CUTOFF_LADDER*, LADDER_* (legacy ones, not REFSCF_LADDER*), SCF_ADVANCE*, FIRST_BLOCK, SNAPSHOT_DIR if orphaned, VERIFY_PROFILE if orphaned; update switch counts. Add one row per removed switch to docs/performance/attempts-ledger.md ("retired in R10a, code deleted; recreate from git history <commit> if needed").
  4. Report exact before/after: lines in src and tests, number of tests (passed/ignored), test binaries, switches; list every deleted test with "tests deleted code".
  Identity: bench_scf --require-identical for both baselines (reference and BAANDI_REFERENCE_STRICT=on oracle), reference_identity.
OWNED_FILES:
  - crates/dft-engine/src/   # the listed modules and their wiring only; not refscf/ or refstart.rs
  - crates/dft-engine/tests/   # deletions of tests of deleted code only
  - crates/dft-engine/Cargo.toml   # test target entries only
  - crates/dft-cli/src/run.rs   # legacy ladder JSON only
  - crates/dft-cli/tests/   # deletions/pins of tests of deleted code only
  - crates/dft-cli/examples/   # only if an example uses deleted code
  - crates/dft-numerics/src/switches.rs   # removals only
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/attempts-ledger.md   # appended rows only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (reference baseline; and BAANDI_REFERENCE_STRICT=on with baseline-oracle.json)
  - python3 -I scripts/make_readme_figures.py && python3 -I scripts/tests/test_readme_figures.py
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: experiments deleted with their tests and switches; baselines identical; README block unchanged; exact before/after counts reported.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
