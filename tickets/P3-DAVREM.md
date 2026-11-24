TASK_ID: P3-DAVREM
GOAL: |
  The reference Davidson's work outside H application and the reduced eigensolve is 9.25 s on dev08 si128 MPI8 (one run, QE-matched clock, 18 solves) against QE's regterg minus h_psi minus rdiaghg of 5.29 s (16 solves); per solve 0.51 s vs 0.33 s. This "remainder" holds the border/Gram products, corrections and preconditioning, the overlap factor, restart rotations, return rotations, the residual norms and every collective. Its GEMMs already run on MKL (REF_DAV_WIDE_ROTATE, REF_DAV_SINGLE_GEMM on). Find where the remaining 1.5x per solve goes and remove it.
READ FIRST: docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md (P3 rows, P3c, P3d, P3f/P3g, M2_* rows), docs/design/reference-scf.md sections "Inner solver cost accounting" and "QE-matched clock table", docs/research/qe-study/davidson-gemm-accounting.md and davidson-iteration-structure-gamma.md. Never read QE source.
TASK: |
  1. Per-stage counters in the reference Davidson at the full diagnostic level (BAANDI_DAVIDSON_STAGES=on), zero cost when off: for border, corrections, factor, restart, return and every other sub-step: nominal FLOPs, bytes copied, number of GEMM calls with their shapes, number and payload of collectives, and seconds. Serialize as davidson_stages.work (rank 0 and max) so achieved GFLOP/s and GB/s per stage can be computed from one run.
  2. From the code and these counters on the Mac fixtures, list non-GEMM work that scales with rows x bands (copies, packing, zero fills, elementwise passes, separate norm passes, per-iteration allocations, replicated small-matrix work that every rank repeats, collectives that could be fused) and remove what is removable behind switches (default off; bit-identical ones may be proposed default on). Keep band decisions and H-column counts unchanged.
  3. Mac checks only (user rule): cargo check --workspace --tests; cargo test -p dft-solvers; cargo test -p dft-engine --test refscf; cargo test -p dft-numerics switches; fmt; clippy -D warnings on touched crates.
  4. docs/performance/lanes/RSCF/davrem-dev08-request.md: ONE profiled run per arm on si64/si128/elongated/rskew MPI8 (development protocol), with the per-stage GFLOP/s and GB/s table.
OWNED_FILES:
  - crates/dft-solvers/src/ref_davidson.rs
  - crates/dft-solvers/src/ref_davidson/
  - crates/dft-solvers/src/tests.rs   # new tests only
  - crates/dft-cli/src/run.rs   # davidson_stages.work JSON only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Inner solver cost accounting" section only
  - docs/performance/lanes/RSCF/davrem-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-numerics -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: per-stage FLOP/byte/collective/seconds counters at full diagnostics, zero cost off; removable non-GEMM work behind switches with unchanged decisions; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
