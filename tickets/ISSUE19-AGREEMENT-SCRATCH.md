TASK_ID: ISSUE19-AGREEMENT-SCRATCH
GOAL: |
  Implement GitHub issue #19 (tickets/issues/ISSUE-19.md). dev08 measured (master, 2026-11-17) that build 405b729 is +6.24 % slower per step than 0fb09c0 with identical switch state and bit-identical counters (si32 np=8), i.e. non-arithmetic code cost. Prime suspect: crates/dft-solvers/src/distributed_davidson.rs ~1605: a zero-filled agreement scratch vec![0.0; 2*cap*cap+2*cap+1] is allocated on every solve right before the STAGE_SETUP timer; also ~1776 vec![0.0; 1 + bands + capacity*bands]; ApplyAgreement::tree_sum_f64 copies the whole input into scratch and the result back to append one failure-status value.
  Do: (1) give the agreement scratch a persistent owner with the lifetime of the distributed solver workspace (grow on capacity change, reuse across SCF maps, no unkeyed global cache, no stale failed-status state); (2) where the reduction buffer is already owned, reserve a status slot and reduce directly in it, removing full in/out copies, deterministic reduction order and failure propagation unchanged; (3) look at the other per-solve allocations on the same path (the 1776 line and any other per-solve vec! in distributed_davidson.rs / workspace/*.rs) and make them workspace-owned too if the same reasoning holds; (4) add counters for allocation bytes/count per solve (profile-off cost zero) so the orchestrator can verify 'after warm-up no new heap allocation' on dev08. Keep all checks (returned-vector X^T X check, collective error agreement, pre-apply clear/nonfinite scan).
OWNED_FILES:
  - crates/dft-solvers/src/distributed_davidson.rs
  - crates/dft-solvers/src/distributed_davidson/workspace.rs
  - crates/dft-solvers/src/distributed_davidson/workspace/lean.rs
  - crates/dft-solvers/tests/distributed_lean_algebra.rs
  - crates/dft-solvers/tests/agreement_scratch.rs   # new
  - docs/performance/attempts-ledger.md   # one row (no timing claim; timing comes from dev08)
READ_ONLY_CONTEXT: AGENTS.md, crates/dft-solvers/AGENTS.md, tickets/issues/ISSUE-19.md, docs/design/mpi-overlap-davidson.md
IMPLEMENTATION_STEPS:
  1. Record energies/eigenvalue bits/collective counts of the existing distributed fixtures on current main BEFORE the change; the same must hold after (bit-identical, same collective count and order).
  2. Add tests: adversarial rank-local operator failure, empty ranks, nonfinite output, unstable overlap factor, returned-vector defect all still detected without hang (reuse existing harness); scratch capacity growth and reuse across two consecutive solves of different shape (no stale status after a failed solve followed by a healthy one).
  3. Implement; add the allocation counter with a test that after one warm-up solve of fixed shape the second solve allocates zero agreement scratch bytes.
TEST_COMMANDS:
  - cargo test -p dft-solvers --offline
  - scripts/check.sh dft-solvers
ACCEPTANCE_CRITERIA: bit-identical healthy-path results and collective ordering; all adversarial tests pass; zero new allocation after warm-up; clippy/fmt clean.
NON_GOALS: removing any hardening check; algorithm changes.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
