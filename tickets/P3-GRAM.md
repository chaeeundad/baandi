TASK_ID: P3-GRAM
GOAL: |
  Speed up the dense GEMMs of the reference Davidson that still run as fixed 256-row/256-term MKL panels. dev08 si128 MPI8 profile at 66d45e69 (W/R/S on): border 4.43 s vs QE regterg:over 1.66 s; corrections 2.73 vs upda+g_psi 1.62; gemm_dispatch "other/mkl" 2353 calls, 328 GFLOP per rank (shape up to 5981 x 5981 x 516: transposed Gram/border products V^T E, V^T HE with the reduction over the ~5981 local rows done as 24 panel GEMMs plus a fixed-order sum). The rotations were fixed by REF_DAV_WIDE_ROTATE (pure fallback -> wide entry); this ticket is about the remaining panelled products.
READ FIRST: docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md (search "panel", "dgemm", "M2_DGEMM", "wide"; do not repeat a rejected idea without new evidence), docs/research/qe-study/davidson-gemm-accounting.md, docs/design/reference-scf.md "Inner solver cost accounting". Never read or request QE source.
TASK: |
  1. Find why the panels exist (determinism across thread counts, scratch bounds) and what the contract requires for the reference driver at one thread per rank (results must stay thread-count independent and rank-layout repeatable; see the rank-count acceptance section of docs/design/reference-scf.md).
  2. Add a single-call route for the reference Davidson's transposed Gram/border products and the correction projections (one dgemm per product over all local rows, MKL in mkl builds), resolved once per solve from a new switch REF_DAV_SINGLE_GEMM (default off). Only the reference Davidson's call sites use it; every other caller and the old path keep their current route bit-identically.
  3. Tests: same decisions and H columns on the existing fixtures with the switch on (rounding class), fixed-layout repeatability (run twice, identical bits), and a thread-count independence test if the pure build can express it. Count the GEMM calls/FLOPs in gemm_dispatch under a distinct route name (e.g. "mkl_single").
  4. Mac checks only (user rule 2026-11-24): cargo check --workspace --tests; cargo test -p dft-numerics; cargo test -p dft-solvers; cargo test -p dft-engine --test refscf; fmt; clippy -D warnings on touched crates.
  5. Write docs/performance/lanes/RSCF/gram-dev08-request.md: untimed screen (default vs on) on si64/si128/elongated/rskew MPI8 with the stage table and gemm_dispatch rows, then A/A + crossed A/B if not slower.
OWNED_FILES:
  - crates/dft-numerics/src/real.rs   # new single-call entry and its route only
  - crates/dft-numerics/src/gemm_dispatch.rs   # new route name only
  - crates/dft-numerics/tests/gemm_dispatch.rs
  - crates/dft-solvers/src/ref_davidson/workspace.rs
  - crates/dft-solvers/src/ref_davidson.rs   # switch resolution only
  - crates/dft-numerics/src/switches.rs   # append the new switch at the end of the list
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/lanes/RSCF/gram-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-numerics
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo fmt --check && cargo clippy -p dft-numerics -p dft-solvers --all-targets -- -D warnings
  - (orchestrator, on dev08) scripts/check.sh full, rank matrix, screen, A/B
ACCEPTANCE_CRITERIA: switch default off; other callers bit-identical; decisions/H columns unchanged on fixtures; repeatability test; dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
