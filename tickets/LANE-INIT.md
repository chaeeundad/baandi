TASK_ID: LANE-INIT
GOAL: |
  Lane INIT: startup before the first SCF map (projector tables, structure factors, atomic start, Ewald, preflight, session construction). Follow docs/performance/lanes/CHARTER.md exactly (diagnose, propose >=2 candidates, implement your choice behind a default-off switch, validate on the Mac, write the dev08 request, stop). You decide what to improve inside your lane; the orchestrator approves from dev08 numbers. Current budget context (dev08, MPI8 si64): total about 10 s; outside the ordinary fine maps about 3.9 s (init 0.88, coarse 1.2, verification 0.98, launcher/IO 0.54, between 0.29); one-core si64: coarse 5.07 s, verification 4.84 s, init 2.56 s. Target: BAANDI wall <= 0.90 x fastest QE (same parallel model).
SCOPE: Time from process start to the first SCF map on one core and MPI ranks. Closed/adopted already: owned projector init, Ewald cache/split/fast, ionic fast, atomic block, init owned states, phase validation fast path (default off).
OWNED_FILES:
  - crates/dft-physics/src/hamiltonian.rs   # construction paths only
  - crates/dft-engine/src/scf.rs   # init section only
  - crates/dft-cli/src/preflight.rs
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/lanes/INIT/
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/performance/lanes/CHARTER.md, docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md, docs/performance/backlog.md
TEST_COMMANDS:
  - scripts/check.sh fast <changed crates>   # 3 CLI cache tests reject in-repo TMPDIR on the Mac
ACCEPTANCE_CRITERIA: diagnosis with numbers; >=2 candidates with expected gains; chosen candidate behind a default-off switch with Mac equivalence evidence; dev08 request with pre-registered counter and threshold; defaults bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
