TASK_ID: P4-DAVREM
GOAL: |
  The reference Davidson's work outside H and the reduced eigensolve now exceeds QE's regterg remainder (dev08 MPI8 profiled snapshot at fa4cfad): si216 davidson_remainder 25.69 s vs QE Davidson regterg - h_psi - rdiaghg = 19.2 s (+6.5 s); si128 6.46 vs 5.21 s (+1.25 s). Stage table si216 (rank max): restart 5.12 s (259 GFLOP, 50.6 GF/s), border_s 4.94 + border_h 4.60 (254 GFLOP each, 53-55 GF/s), return 3.89 (201 GFLOP), corrections 3.41 (40 GF/s), setup 2.11 (5.5 GB streamed), transform 1.23, coefficients 0.83, factor 0.78, normalize 0.76, precondition 0.55, apply_check 0.52. The GEMMs already run near the per-core peak, so the lever is fewer FLOPs and fewer, larger calls.
READ FIRST: docs/design/reference-scf.md "Inner solver cost accounting"; crates/dft-solvers/src/ref_davidson.rs and ref_davidson/workspace.rs; crates/dft-mkl/src/adapter.rs (dgemm entries); attempts-ledger rows P3 (WIDE_ROTATE, SINGLE_GEMM, IDEAL_RESET, STREAM), the earlier rejected "upper-border split" experiment and why; docs/research/qe-study/davidson-gemm-accounting.md. Never read QE source.
TASK: |
  Behind separate switches (default off), each with H-column/decision identity tests on the fixtures:
  1. Symmetric borders with MKL dgemmt (only the needed triangle of V^T E and V^T HE where the block is square/symmetric), or one GEMM V^T [E | HE] (both borders in one call, twice the width) if dgemmt is not usable for the rectangular border: report FLOPs saved.
  2. Restart and return as one GEMM each over the row-stacked [V; HV] (2r x m) times C instead of two GEMMs, when memory layout allows (or an equivalent single call with a strided batch).
  3. Davidson setup/append/normalize/apply_check: remaining streamed GB; fuse what is left.
  Add a dgemmt entry to the MKL adapter if needed (unsafe only in the adapter, preconditions documented and tested). Rounding class allowed; decisions and H columns must stay equal on fixtures (rank-rule branch recording if not).
  Mac checks: cargo check (also --features dft-numerics/mkl), cargo test -p dft-solvers, -p dft-engine --test refscf, -p dft-mkl, fmt, clippy. docs/performance/lanes/RSCF/p4-davrem-dev08-request.md: ONE profiled run per arm, si128 and si216 MPI8 verify off, stage table with GFLOP and GF/s.
OWNED_FILES:
  - crates/dft-solvers/src/ref_davidson.rs
  - crates/dft-solvers/src/ref_davidson/
  - crates/dft-numerics/src/real.rs   # new dispatch entries only
  - crates/dft-mkl/src/adapter.rs   # dgemmt entry only
  - crates/dft-mkl/tests/   # adapter contract tests only
  - crates/dft-solvers/src/tests.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Inner solver cost accounting" section only
  - docs/performance/lanes/RSCF/p4-davrem-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo check -p dft-solvers -p dft-numerics -p dft-mkl --tests --features dft-numerics/mkl
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-mkl
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-numerics -p dft-mkl --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: FLOP/call reductions implemented behind switches with identity tests; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
