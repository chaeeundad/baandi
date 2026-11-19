TASK_ID: LANE-DAVIDSON
GOAL: |
  Lane DAVIDSON: Davidson eigensolver dense algebra: subspace Gram/Rayleigh-Ritz, reduced eigenproblem, orthogonalization, residual blocks, restarts, cluster handling, MPI reductions. Follow docs/performance/lanes/CHARTER.md exactly (diagnose, propose >=2 candidates, implement your choice behind a default-off switch, validate on the Mac, write the dev08 request, stop). You decide what to improve inside your lane; the orchestrator approves from dev08 numbers. Current budget context (dev08, MPI8 si64): total about 10 s; outside the ordinary fine maps about 3.9 s (init 0.88, coarse 1.2, verification 0.98, launcher/IO 0.54, between 0.29); one-core si64: coarse 5.07 s, verification 4.84 s, init 2.56 s. Target: BAANDI wall <= 0.90 x fastest QE (same parallel model).
SCOPE: Per-iteration time outside H application: GEMM/SYRK shapes, reduced eigen cost (REL-SLICES closed; Advance/bounded-subspace and restart-keep rejected), allreduce count and sizes, work per map; also tolerance schedule per map. Read davidson-algebra-gap-analysis and diagonalization-plan.
OWNED_FILES:
  - crates/dft-solvers/src/lean_davidson.rs
  - crates/dft-solvers/src/distributed_davidson.rs
  - crates/dft-solvers/src/real_davidson.rs
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/lanes/DAVIDSON/
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/performance/lanes/CHARTER.md, docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md, docs/performance/backlog.md, docs/performance/davidson-algebra-gap-analysis.md, docs/performance/diagonalization-plan.md
TEST_COMMANDS:
  - scripts/check.sh fast <changed crates>   # 3 CLI cache tests reject in-repo TMPDIR on the Mac
ACCEPTANCE_CRITERIA: diagnosis with numbers; >=2 candidates with expected gains; chosen candidate behind a default-off switch with Mac equivalence evidence; dev08 request with pre-registered counter and threshold; defaults bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
