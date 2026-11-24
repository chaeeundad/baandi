TASK_ID: P3-RMM
GOAL: |
  User decision 2026-11-25 (A): add QE-style rmm-davidson to the reference SCF driver as a switchable module (REFSCF_RMM, default off). Implement it from docs/research/qe-study/qe-rmm-davidson.md ONLY (an own-words note; you never get QE source and must not look for it). On dev08 si128 MPI8 the non-H Davidson work is ~6.5 s per run (border, corrections, restart, return, reduced eigensolve); QE's rmm-davidson is 2-3% faster than its own Davidson on our cells, and its gain comes from removing most non-H work at the price of more H columns.
READ FIRST: the note above (sections 2, 4-8 are the contract, section 9 explains why our earlier hybrid A22/A23 failed and what must be replicated), docs/design/reference-scf.md "Inner solver" and "Inner solver cost accounting", crates/dft-solvers/src/ref_davidson.rs and ref_davidson/ (reuse its distributed operator trait, GEMM workspace routes, preconditioner pieces and stage counters), docs/performance/AGENT-BRIEF.md, attempts-ledger.md A22/A23.
TASK: |
  1. New solver module crates/dft-solvers/src/ref_rmm.rs (distributed, Gamma real, same operator trait and collective interface as ref_davidson): per-call Rayleigh-Ritz rotation in the new H, one sweep of at most 4 DIIS steps per unconverged band with the note's constrained DIIS (canonical orthogonalization, fallback to newest entry), the published Teter-Payne-Allan / Kresse-Furthmuller preconditioner scaled by band kinetic energy, the parabola step length clamped to [0.05, 1], band freezing by eigenvalue change (occupied ethr, empty max(5 ethr, 1e-5 Ry) converted to Ha consistently), block Gram-Schmidt in band order carrying HX, diagonal Rayleigh quotients and sorting. Record the smallest Gram-Schmidt norm before normalization. Stage counters like ref_davidson (seconds, H columns, GEMM FLOPs) at the full diagnostic level.
  2. Driver: with REFSCF_RMM on, choose per SCF iteration by the note's schedule (Davidson at iterations 1, 2 and every 5th; RMM otherwise; first-step redo keeps index 1); verification maps always use Davidson (note section 8 proposal). Everything else (ethr rule, stopping estimator, verification criterion, mixing) unchanged. Keep driver edits minimal and local (P3-ENERGY edits the per-map energy code in refscf/ in parallel).
  3. Tests on the Mac (only cargo check and touched crates' tests): fixed-H tests against the dense oracle (eigenvalues of a converged sweep sequence), SCF fixtures (cubic, skew, elongated, synthetic nonlocal) with RMM on: converge and verify, delivered energy within 1e-10 Ha of the Davidson run, report maps, H columns and dense FLOPs per arm; simulated ranks 1/3/8 repeatability.
  4. docs/performance/lanes/RSCF/rmm-dev08-request.md: ONE profiled run per arm (D default, M = REFSCF_RMM on) on si64/si128/elongated/rskew MPI8, verify on and off for si128, with stage table, H columns, maps, E - E_strict; plus QE rmm-davidson vs QE Davidson on the same cells for reference.
OWNED_FILES:
  - crates/dft-solvers/src/ref_rmm.rs   # new
  - crates/dft-solvers/src/lib.rs   # module declaration only
  - crates/dft-solvers/src/ref_davidson.rs   # pub(crate) visibility of reusable pieces only
  - crates/dft-solvers/src/ref_davidson/   # same
  - crates/dft-engine/src/refscf/driver.rs   # solver selection per iteration only
  - crates/dft-engine/src/refscf/map.rs   # solver call dispatch only
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-cli/src/run.rs   # rmm stage JSON only
  - docs/design/reference-scf.md   # a "RMM-Davidson" section only
  - docs/performance/lanes/RSCF/rmm-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-engine -p dft-numerics -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: REFSCF_RMM default off; fixtures converge and verify with delivered energy within 1e-10 Ha of Davidson; per-arm counts reported; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
