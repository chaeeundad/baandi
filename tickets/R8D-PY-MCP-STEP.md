TASK_ID: R8D-PY-MCP-STEP
GOAL: |
  Stage R8d of docs/roadmap/refactor-r-plan.md (user decision D1): Python stepwise control (Session.initialize/step/status/iteration) and the MCP tools (create_calculation, run, step, pause, resume, cancel, propose_policy_change) drive the reference stepper (R8a) and the reference runtime adapter (R8c, crates/dft-engine/src/refscf/runtime.rs) instead of the legacy ScfSession; remove the L8 routing predicate for the reference path.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 1b, 2 rows L7-L8, 4 "R8"); crates/dft-engine/src/refscf/{stepper.rs,runtime.rs,api.rs}; crates/dft-cli/src/telemetry.rs (how R8c drives the runtime); crates/baandi-py/src/lib.rs (R5 gated stepwise calls with NotImplementedError); crates/baandi-mcp/src/lib.rs, tests/server.rs, examples/tool_schemas.rs; crates/dft-cli/src/driver_select.rs (L8).
TASK: |
  1. baandi-py: initialize/step/status/iteration work on the reference driver through the stepper/runtime adapter; step returns the StepReport fields; a stepwise run to completion is bit-identical to run().
  2. baandi-mcp: tools act on the reference runtime; propose_policy_change accepts only the reference-meaningful policies (mixing beta, history, threshold, verification) and rejects others with a clear message; tool schemas regenerated (examples/tool_schemas.rs) if they change.
  3. Remove L8 for the reference path; legacy remains reachable with driver="legacy".
  4. Tests: Python native test (skipped if the extension is not built) stepwise == run bit-identical; MCP server tests updated for the reference runtime (create, step, pause/resume, cancel, policy change at a safe point).
OWNED_FILES:
  - crates/baandi-py/
  - crates/baandi-mcp/
  - crates/dft-cli/src/driver_select.rs   # L8 only
  - crates/dft-cli/tests/fixtures/routing/   # L8 expectations only
  - crates/dft-cli/tests/driver_routing.rs   # L8 expectations only
  - crates/dft-engine/src/refscf/runtime.rs   # API additions only if needed
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p baandi-py
  - cargo test -p baandi-mcp
  - cargo test -p dft-cli --test driver_routing
  - cargo fmt --check && cargo clippy -p baandi-py -p baandi-mcp --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: Python stepping and MCP tools on the reference runtime, bit-identical to run(); L8 removed for reference.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
