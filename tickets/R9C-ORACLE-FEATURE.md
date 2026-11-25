TASK_ID: R9C-ORACLE-FEATURE
GOAL: |
  Stage R9c of docs/roadmap/refactor-r-plan.md (user decision D2: the strict oracle is moved whole): put the legacy driver (the strict oracle and the remaining legacy standard path) behind a Cargo feature `strict-oracle` of dft-engine, default on, so that `cargo check -p dft-engine --no-default-features` builds the production engine without any legacy driver code. Pure moves: bit-identical by construction.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 3 layout and "Where the oracle lives", 4 "R9", 5 risks, 7 decisions); crates/dft-engine/src/lib.rs (after R9a/R9b: crate::model, top-level refscf, the independence test crates/dft-engine/tests/refscf_independence.rs); crates/dft-engine/Cargo.toml; crates/dft-cli/Cargo.toml and how dft-cli uses legacy entry points (driver_select fallbacks, strict oracle, legacy pins in tests); scripts/check.sh; crates/dft-engine/AGENTS.md.
TASK: |
  1. git mv the legacy driver files into crates/dft-engine/src/oracle/ (scf.rs, session.rs, convergence.rs, collinear.rs, finite_temperature.rs, band_density.rs, atomic_start.rs if not shared, start_density.rs if not shared, spectral.rs, diagnostics.rs, nscf.rs (legacy), events.rs, runtime.rs/policy.rs if only legacy uses them, gamma_real.rs, m2.rs, mix_dist.rs, cutoff_ladder.rs remains, fixed_map.rs if legacy-only, verification_tests.rs, checkpoint.rs v1, snapshot.rs...): decide per file from actual dependencies (anything the reference path or shared physics needs stays outside oracle/). Keep `pub use` re-exports at the old paths under the feature so dft-cli/baandi-py/baandi-mcp compile unchanged with the feature on.
  2. Cargo: `[features] default = ["strict-oracle"]` in crates/dft-engine/Cargo.toml; `#[cfg(feature = "strict-oracle")]` on the oracle module and re-exports; tests that use legacy entry points get `required-features = ["strict-oracle"]` (or a cfg at the top); dft-cli declares the feature explicitly.
  3. scripts/check.sh gains a gate: `cargo check -p dft-engine --no-default-features` and `cargo clippy -p dft-engine --no-default-features --all-targets -- -D warnings` (fast; report the time added).
  4. crates/dft-engine/AGENTS.md: oracle section (where it lives, that it is frozen, how tests use it).
  5. Identity: identical test counts (report), reference_identity, oracle_identity, bench_scf --require-identical for both baselines.
OWNED_FILES:
  - crates/dft-engine/src/   # moves, cfg attributes, re-exports only
  - crates/dft-engine/Cargo.toml
  - crates/dft-engine/tests/   # required-features / cfg only
  - crates/dft-engine/AGENTS.md
  - crates/dft-cli/Cargo.toml   # feature declaration only
  - scripts/check.sh   # the new gate only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo check -p dft-engine --no-default-features
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (reference baseline; and BAANDI_REFERENCE_STRICT=on with baseline-oracle.json)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: production engine builds without legacy code; default build unchanged and bit-identical; new check.sh gate.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
