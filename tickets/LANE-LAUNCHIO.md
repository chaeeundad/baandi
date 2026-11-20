PHASE 0 ONLY (user decision 2026-11-20): profiling is NOT done on the Mac. See docs/performance/lanes/CHARTER.md top section and tickets/LANE-PHASE0-msg.md (read both first; they are authoritative over the text below). This run = Phase 0: add default-off or PHASE_PROFILE-gated instrumentation for your lane, write docs/performance/lanes/LAUNCHIO/profile-request.md with exact dev08 commands and scripts/lane_launchio_parse.py, verify bit-identical results on the Mac, return status PROFILE_REQUEST and stop. Known dev08 context for your lane (si64 np8 A profile, instrumented): density 0.746 s (vs QE sum_band 0.34), energy 0.291 s, potential 0.134 s, mixing 0.032 s, io/materialize 0.296 s, between 0.364 s (session_setup 0.235, ladder_handoff 0.123), outside-process launch/exit/output about 0.84 s vs QE 0.12 s.

TASK_ID: LANE-LAUNCHIO
GOAL: |
  Lane LAUNCHIO: process launch, input parsing, pseudopotential loading, output writing, results.json, MPI init/finalize, between-map glue in the engine driver. Follow docs/performance/lanes/CHARTER.md exactly (diagnose, propose >=2 candidates, implement your choice behind a default-off switch, validate on the Mac, write the dev08 request, stop). You decide what to improve inside your lane; the orchestrator approves from dev08 numbers. Current budget context (dev08, MPI8 si64): total about 10 s; outside the ordinary fine maps about 3.9 s (init 0.88, coarse 1.2, verification 0.98, launcher/IO 0.54, between 0.29); one-core si64: coarse 5.07 s, verification 4.84 s, init 2.56 s. Target: BAANDI wall <= 0.90 x fastest QE (same parallel model).
SCOPE: Launcher/IO about 0.54 s at si64 np8 (F-LAUNCH-IO closed pieces). Measure launch-to-exit minus SCF; find per-rank redundant work, file reads/writes, serialization, allocator/page-fault costs.
OWNED_FILES:
  - crates/dft-cli/src/run.rs
  - crates/dft-cli/src/checkpoint_io.rs
  - crates/dft-engine/src/runtime.rs
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/lanes/LAUNCHIO/
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/performance/lanes/CHARTER.md, docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md, docs/performance/backlog.md
TEST_COMMANDS:
  - scripts/check.sh fast <changed crates>   # 3 CLI cache tests reject in-repo TMPDIR on the Mac
ACCEPTANCE_CRITERIA: diagnosis with numbers; >=2 candidates with expected gains; chosen candidate behind a default-off switch with Mac equivalence evidence; dev08 request with pre-registered counter and threshold; defaults bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
