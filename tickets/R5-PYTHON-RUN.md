TASK_ID: R5-PYTHON-RUN
GOAL: |
  Stage R5 of docs/roadmap/refactor-r-plan.md: the Python binding's Session.run(), result(), forces() and stress() use the reference engine facade refscf::api::solve (R3, merged). Stepwise control (initialize/step/status/iteration) stays on an explicit legacy driver until the reference stepper (R8) exists.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 1b, 4 "R5", 7 decisions); crates/baandi-py/src/lib.rs, VALIDATION.md, AGENTS.md; crates/dft-engine/src/refscf/api.rs; crates/dft-cli/src/run.rs (how the CLI builds StartInput/Request and the driver selection in crates/dft-cli/src/driver_select.rs).
TASK: |
  1. Options gains driver = "reference" (default) | "legacy". run()/result()/forces()/stress() with reference go through api::solve (same option mapping as the CLI: mixing beta 0.7 / history 8 defaults of the reference path, verification on). With driver="legacy" the current behaviour is kept. initialize()/step()/status()/iteration() on the reference driver raise NotImplementedError("stepwise control requires driver='legacy' until the reference stepper (R8)").
  2. Native test (cargo test, no Python needed if the crate has Rust-side tests; otherwise a pytest that is skipped when the extension is not built): energy, forces and stress from the binding equal the CLI on the same input (bit-identical, same engine call) for Si2 Gamma and k222.
  3. Document the default change (mixing defaults, verification) in VALIDATION.md; keep old evidence, add new evidence as a new file if regenerated.
OWNED_FILES:
  - crates/baandi-py/src/lib.rs
  - crates/baandi-py/VALIDATION.md
  - crates/baandi-py/AGENTS.md
  - crates/baandi-py/tests/   # new tests only
  - crates/baandi-py/python/   # tests/stubs only if present
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p baandi-py
  - cargo fmt --check && cargo clippy -p baandi-py --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: Python run/forces/stress on the reference facade, bit-identical to the CLI; stepwise calls clearly gated.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
