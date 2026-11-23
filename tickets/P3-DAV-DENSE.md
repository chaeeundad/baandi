TASK_ID: P3-DAV-DENSE
GOAL: |
  Cut the non-H dense-algebra cost of the reference Davidson solver (crates/dft-solvers/src/ref_davidson.rs and the distributed workspace it calls). dev08 MPI8 profile at main 2aee083 (median of 3, rank max, SCF + verification solves; QE 7.6 Gamma regterg lines for the same cells, approximate pairing):
    si128: restart 8.500 s (QE: no separate cost) | border 5.495 vs QE over 1.655 | corrections 5.953 vs QE upda+g_psi 1.615 | return 3.874 vs QE last 1.780 | factor 0.384 | reduced_eigen 2.647 vs rdiaghg 2.330 | apply_h 16.850 vs h_psi 11.140. Work: 20 solver calls, 83 Davidson iterations, 14203 H columns (QE 16 c_bands, 82 h_psi calls).
    si64: restart 0.696 | border 0.707 vs 0.240 | corrections 0.906 vs 0.19 | return 0.378 vs 0.220.
  So about 10 s of si128's 30 s gap to QE is in restart + border + corrections + return.
READ FIRST: docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md (do not repeat rejected ideas without new evidence), docs/performance/backlog.md, docs/research/qe-study/davidson-iteration-structure-gamma.md (incl. its preconditioner addendum), docs/research/qe-study/davidson-gemm-accounting.md, docs/research/qe-study/qe-solvers-and-subspace.md, docs/design/reference-scf.md ("Inner solver"). Never read or request QE source.
TASK: |
  1. Account the flops and data movement of each stage per iteration (restart, border, corrections, return) for n bands, basis m <= 2n, local rows r, as the code does them now; compare with the minimal work the algorithm needs (the QE note's accounting). Write it into docs/design/reference-scf.md section "Inner solver cost accounting".
  2. Remove the excess: e.g. restart that rebuilds or re-orthonormalizes what can be rotated, rotations done as several GEMMs that can be one, recomputation of overlaps already known, per-iteration allocation/copies, replicated reductions that can be fused into one collective, corrections computed for converged bands. Keep the algorithm and its stopping/selection decisions identical: the band-convergence decisions, iteration counts and H-column counts must not change except by rounding (record any rounding branch as the rank-rule does).
  3. Each independent change behind its own switch (crates/dft-numerics/src/switches.rs, default ON only if it is an overhead-only change proven bit-identical in a test; otherwise default off), so dev08 can A/B it.
  4. Mac checks only (user rule 2026-11-24; never check.sh full or release matrices on the Mac): cargo check --workspace --tests; cargo test -p dft-solvers; cargo test -p dft-engine --test refscf; cargo test -p dft-numerics switches; fmt; clippy -D warnings on touched crates. Add unit tests proving identical decisions/columns on small fixtures, and a flop/byte count test where practical.
  5. Write docs/performance/lanes/RSCF/dav-dense-dev08-request.md: same si64/si128 MPI8 inputs as the verify-cost A/B, default vs each switch (and all together), A/A then crossed pairs, plus a profiled run with BAANDI_DAVIDSON_STAGES=on per arm for the stage table. Add a ledger row draft in the report (the orchestrator commits the ledger).
OWNED_FILES:
  - crates/dft-solvers/src/ref_davidson.rs
  - crates/dft-solvers/src/distributed_davidson/   # only code reached from ref_davidson; the old path's behaviour must stay bit-identical
  - crates/dft-solvers/src/tests.rs   # new tests only
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Inner solver cost accounting" section only
  - docs/performance/lanes/RSCF/dav-dense-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-numerics --all-targets -- -D warnings
  - (orchestrator, on dev08) scripts/check.sh full, rank matrix, A/B
ACCEPTANCE_CRITERIA: accounting written; changes switchable; identical decisions/H columns on fixtures; old path untouched; dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
