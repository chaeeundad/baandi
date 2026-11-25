TASK_ID: R8C-RUNTIME
GOAL: |
  Stage R8c (and the CLI part of R8d) of docs/roadmap/refactor-r-plan.md (user decision D1): the observe/agent runtime, the control channel and the self-healing demo drive the reference stepper (R8a, merged) instead of the legacy ScfSession; then the routing predicate L7 (runtime modes) is removed for the reference path.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 1a telemetry/control_channel/heal rows, 2 rows L7-L8, 4 "R8", 5 risk 5); docs/design/reference-scf.md "Stepper"; crates/dft-engine/src/refscf/stepper.rs and api.rs; crates/dft-cli/src/telemetry.rs, control_channel.rs, heal.rs; crates/dft-engine/src/runtime.rs, policy.rs, events.rs (legacy runtime contracts: events, safe points, pause/resume/cancel, policy changes); tests/reference/runtime/generated/results.json (README-cited evidence).
TASK: |
  1. A reference runtime adapter: run --observe/--agent on the reference driver emits the same event stream schema per stepper step (iteration, estimator, dr2, energy when computed, ethr, H columns, kind), supports pause/resume/cancel at safe points, and accepts only reference-meaningful policy changes (mixing beta, history, threshold, verification on/off); other legacy policy knobs are rejected with a clear message on the reference driver.
  2. control_channel and heal use the adapter; the self-healing demo's scenario is reproduced on the reference path (or documented as not applicable with the reason). New evidence goes to a new file next to tests/reference/runtime/generated/results.json (do not overwrite the old one).
  3. Remove the L7 predicate for the reference path in driver_select.rs (legacy still reachable with driver="legacy").
  4. Tests: an observe run equals a plain run bit-for-bit (energy bits, density hash) with a recorded event stream; pause/resume and a beta change at a safe point behave deterministically; cancel leaves a consistent partial record; ranks 1/3 where the runtime supports it.
OWNED_FILES:
  - crates/dft-cli/src/telemetry.rs
  - crates/dft-cli/src/control_channel.rs
  - crates/dft-cli/src/heal.rs
  - crates/dft-cli/src/run.rs   # observe/agent dispatch only
  - crates/dft-cli/src/driver_select.rs   # L7 only
  - crates/dft-cli/tests/fixtures/routing/   # L7 expectations only
  - crates/dft-cli/tests/driver_routing.rs   # L7 expectations only
  - crates/dft-cli/tests/reference_runtime.rs   # new
  - crates/dft-engine/src/refscf/runtime.rs   # new adapter, if engine-side is cleaner
  - crates/dft-engine/src/refscf/mod.rs   # declaration only
  - tests/reference/runtime/generated/reference-results.json   # new
  - docs/design/reference-scf.md   # a "Runtime control" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-cli --test reference_runtime
  - cargo test -p dft-cli
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - cargo fmt --check && cargo clippy -p dft-cli -p dft-engine --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: observe/agent/control/heal work on the reference stepper with bit-identical results to plain runs; L7 removed for reference; new evidence file.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
