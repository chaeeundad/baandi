TASK_ID: LANE-LADDERMAPS
GOAL: |
  Lane LADDERMAPS: cutoff-ladder structure costs that are NOT the fraction/stop rules (closed, see docs/design/ladder-fraction.md and ladder-2x2.md): per-map cost of the coarse stage, handoff construction, coarse operator setup, coarse-stage tolerance schedule. Follow docs/performance/lanes/CHARTER.md exactly (diagnose, propose >=2 candidates, implement your choice behind a default-off switch, validate on the Mac, write the dev08 request, stop). You decide what to improve inside your lane; the orchestrator approves from dev08 numbers. Current budget context (dev08, MPI8 si64): total about 10 s; outside the ordinary fine maps about 3.9 s (init 0.88, coarse 1.2, verification 0.98, launcher/IO 0.54, between 0.29); one-core si64: coarse 5.07 s, verification 4.84 s, init 2.56 s. Target: BAANDI wall <= 0.90 x fastest QE (same parallel model).
SCOPE: Coarse stage about 1.2 s at si64 np8 and 5.07 s one core. Fraction rules, coarse-stop-only, owned handoff, minimum grid, margin indicator are closed negative. Look for per-map overhead inside the coarse stage and handoff build, not for a new stop rule.
OWNED_FILES:
  - crates/dft-engine/src/cutoff_ladder.rs
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/lanes/LADDERMAPS/
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/performance/lanes/CHARTER.md, docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md, docs/performance/backlog.md, docs/design/ladder-fraction.md, docs/design/ladder-2x2.md, docs/design/ladder-handoff.md
TEST_COMMANDS:
  - scripts/check.sh fast <changed crates>   # 3 CLI cache tests reject in-repo TMPDIR on the Mac
ACCEPTANCE_CRITERIA: diagnosis with numbers; >=2 candidates with expected gains; chosen candidate behind a default-off switch with Mac equivalence evidence; dev08 request with pre-registered counter and threshold; defaults bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
