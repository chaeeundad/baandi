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
