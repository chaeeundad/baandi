TASK_ID: R8B-CHECKPOINT
GOAL: |
  Stage R8b of docs/roadmap/refactor-r-plan.md (user decision D1): checkpoint and restart for the reference path on top of the R8a stepper (merged; crates/dft-engine/src/refscf/stepper.rs). Restart must be bit-identical to the uninterrupted run.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 2 row L6, 4 "R8", 5 risk 5, 7 decisions); docs/design/reference-scf.md "Stepper"; crates/dft-engine/src/refscf/{stepper.rs,driver.rs,k_driver.rs,electronic.rs,broyden.rs,kept_hx.rs}; the legacy checkpoint v1: crates/dft-engine/src/checkpoint.rs, crates/dft-cli/src/checkpoint_io.rs, checkpoint_inspect.rs, docs/checkpoint-format.md, tests/reference/checkpoint/results.json (the identity protocol to reproduce).
TASK: |
  1. crates/dft-engine/src/refscf/checkpoint.rs: serialize the reference state at a stepper safe point: density (owned planes gathered or per-rank files, choose and document), Broyden history, kept orbitals/HX (or a documented rule that they are rebuilt bit-identically), ethr ratchet, iteration counters, verification state, start metadata, input/UPF hashes and binary identity. Restore it into a stepper so that continuing gives the same bits as the uninterrupted run.
  2. CLI: checkpoint/restart commands on the reference driver (format v2 with a version field; v1 files keep loading on the legacy path); checkpoint_inspect reads v2. Remove the L6 routing predicate for the reference path.
  3. Tests: memory and file restart identity (energy bits, density hash, orbitals hash) at several safe points on Gamma, k222s and a skew two-species fixture, ranks 1/3/8; corrupted/mismatched-input files are rejected with clear errors. Write new evidence to tests/reference/checkpoint/reference-v2-results.json (keep the v1 file).
  4. docs/checkpoint-format.md: v2 section.
OWNED_FILES:
  - crates/dft-engine/src/refscf/checkpoint.rs   # new
  - crates/dft-engine/src/refscf/stepper.rs   # restore hooks only, no arithmetic change
  - crates/dft-engine/src/refscf/mod.rs   # declaration only
  - crates/dft-engine/src/refscf/broyden.rs   # serialization only
  - crates/dft-engine/src/refscf/kept_hx.rs   # serialization only
  - crates/dft-cli/src/checkpoint_io.rs
  - crates/dft-cli/src/checkpoint_inspect.rs
  - crates/dft-cli/src/driver_select.rs   # L6 only
  - crates/dft-cli/tests/fixtures/routing/   # L6 expectations only
  - crates/dft-cli/tests/driver_routing.rs   # L6 expectations only
  - crates/dft-cli/tests/checkpoint_v2.rs   # new
  - crates/dft-engine/tests/refscf_checkpoint.rs   # new
  - tests/reference/checkpoint/reference-v2-results.json   # new
  - docs/checkpoint-format.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf_checkpoint
  - cargo test -p dft-cli --test checkpoint_v2
  - cargo test -p dft-cli --test checkpoint_io
  - scripts/check.sh (BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (release, default baseline)
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: reference checkpoint/restart bit-identical to uninterrupted runs; v1 still loads on legacy; baselines identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
