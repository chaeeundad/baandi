TASK_ID: P3-ADDV
GOAL: |
  Per H call on dev08 si128 MPI8 (one run at 31b23e11, rank 0, 93 calls): nonlocal back-application add_vuspsi 0.0246 s vs QE 0.0185 s (+33%), projection calbec 0.0193 vs 0.0172 (+12%), local part vloc_psi 0.1046 vs 0.097 (+8%). The projection now uses one MKL call for half the H blocks (NONLOCAL_SINGLE_GEMM); the back-application P W (rows x projectors x columns) and the coupling step were not changed.
READ FIRST: docs/performance/AGENT-BRIEF.md, attempts-ledger.md rows P3l, P3-NL, M2_DGEMM, M2-OWNED-FIELDS; docs/research/qe-study/qe-hpsi-and-fft.md. Never read QE source.
TASK: |
  1. Account the back-application and coupling per H call: GEMM shapes, routes (gemm_dispatch), whether the panelled 256-row path is used, copies of W, zero fills of the output before accumulation, and the add into HX. Why does NONLOCAL_SINGLE_GEMM cover only 103 of 206 projection calls on si128 (which blocks fall back and why)?
  2. Behind switches (default off): back-application as one dgemm over all local rows accumulating directly into HX (beta = 1, no temporary), coupling fused into the projection result, single-GEMM path for the remaining projection calls when the shape allows. Bit-identical where the operation order is unchanged, otherwise rounding class.
  3. Mac checks only: cargo check --workspace --tests; cargo test -p dft-physics; cargo test -p dft-engine --test refscf; cargo test -p dft-numerics; fmt; clippy -D warnings.
  4. docs/performance/lanes/RSCF/addv-dev08-request.md: ONE profiled run per arm on si64/si128/elongated/rskew MPI8 with calbec/add_vuspsi per call and the nonlocal GEMM routes.
OWNED_FILES:
  - crates/dft-physics/src/distributed_real_gamma.rs   # nonlocal projection/back-application only
  - crates/dft-physics/src/   # nonlocal operator files reached from it only
  - crates/dft-numerics/src/real.rs   # only if a beta=1 accumulate entry is needed for the single-call route
  - crates/dft-numerics/src/switches.rs   # append at the end (P3-RESET also appends)
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/lanes/RSCF/addv-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-physics
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics
  - cargo fmt --check && cargo clippy -p dft-physics -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: per-call accounting; switchable back-application/projection improvements with H equality tests; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
