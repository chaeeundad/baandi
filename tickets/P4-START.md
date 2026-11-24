TASK_ID: P4-START
GOAL: |
  Setup is now the largest single excess over QE (dev08 MPI8 profiled snapshot at fa4cfad, verify off): si128 engine setup 2.125 s vs QE init_run 0.76 s; si216 5.834 s vs 2.90 s. Rows (rank 0): si216 atomic_block_total 3.604 (atomic_h_apply 1.764, atomic_h_apply_gram 2.056 incl. the apply), local_tables_structure_factors 0.802, overlap_factor_inverse 0.552, generalized_diagonalize_rotate 0.391, orbital_fill 0.292, field_fft_plan 0.286, operator plans 0.248, overlap_gram 0.239, nonlocal_diagonal 0.234, preflight 0.233; si128 atomic_block_total 1.043, local tables 0.405, field_fft_plan 0.164, overlap_factor_inverse 0.128, preflight 0.126. Target: setup at or below QE's init_run.
READ FIRST: docs/design/reference-scf.md "Start" and "Setup and output cost"; attempts-ledger rows REFSTART_MKL, P3-SETUP..P3-SETUP4, P3t; docs/performance/AGENT-BRIEF.md; docs/research/qe-study/qe-scf-cycle.md. Never read QE source.
TASK: |
  1. atomic_h_apply: H applied to all atomic orbitals (si216: 1152 columns). Check it goes through the same batched pair FFT, NONLOCAL_SINGLE_GEMM/BACK_ACCUM and kept plan paths as the SCF H (per-column cost vs the SCF's ~0.38 ms/column at si216); fix any slower route.
  2. overlap_factor_inverse and generalized_diagonalize_rotate: replicated dense work on 1152x1152 (si216); use MKL dpotrf/dtrsm or a Cholesky-based generalized eigensolve (dsygvd-style: one factor, one standard dsyevd) instead of an explicit inverse; check it is not done redundantly.
  3. local_tables_structure_factors and orbital_fill: owned G rows only, vectorized, one pass per species; nonlocal_diagonal and operator plans reuse shared tables.
  4. field_fft_plan / operator plans: share with the plans already created (REFSTART_FFT_SHARE exists for the density plan; extend to field/operator plans).
  5. Everything behind switches (default off; bit-identical ones may be proposed default on), Mac checks only (cargo check, touched crates' tests, fmt, clippy); grep tests before any default flip. docs/performance/lanes/RSCF/p4-start-dev08-request.md: ONE profiled run per arm, si64/si128/si216 MPI8, setup_s split next to QE init_run.
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs
  - crates/dft-engine/src/refscf/   # plan sharing only
  - crates/dft-solvers/src/distributed_fft.rs   # plan API only
  - crates/dft-numerics/src/linalg.rs   # only if a Cholesky-based generalized eigensolve helper is added
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Setup and output cost" section only
  - docs/performance/lanes/RSCF/p4-start-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-solvers
  - cargo test -p dft-numerics
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-solvers -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: each setup item addressed or explained; start state within rounding; switches default off; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
