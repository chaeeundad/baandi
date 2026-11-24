TASK_ID: P3-KEPTHX2
GOAL: |
  Work-order screen (dev08 Si64 MPI8, bb416c73, ledger "Work-order screen"): REFSCF_KEPT_HX replaced 1452 starting full-H columns (fine full-H columns 6094 -> 4642, 11 cache hits) and solve fell 0.15 s, but density rose 0.371 -> 0.850 s (the retained-block density pass that keeps real-space pair slabs) and peak RSS +53 MB/rank; wall 6.41 -> 6.60 s. Make the warm start pay: the density step must cost no more than without the switch.
READ FIRST: docs/design/reference-scf.md "Kept HX warm start"; crates/dft-engine/src/refscf/kept_hx.rs; the density path of the reference map; band_density_sum_keeping and apply_local_change in crates/dft-physics/src/distributed_real_gamma.rs; docs/performance/AGENT-BRIEF.md; attempts-ledger rows M2_REUSE, P3-KEPTHX, Work-order screen. Never read QE source.
TASK: |
  1. Account the density pass with and without keeping (FFTs, copies into slabs, extra passes, allocations, slab layout and stride). The inverse FFT of each band pair is already done for the density; keeping should only add a store of the real-space values that are already in cache.
  2. Fix: write the retained slab directly from the density pass's inverse transform buffer (no second transform, no extra full pass, preallocated slabs reused across maps, contiguous layout that apply_local_change reads), or, if storing is inherently costly, a variant that keeps only what apply_local_change needs. Behind the existing switch; results unchanged (bit-identical to the current kept-HX path).
  3. Tests: fast_order_ kept-HX test still passes; add a check that the keeping density pass performs the same number of FFTs as the plain density pass (counters).
  4. Mac checks only. docs/performance/lanes/RSCF/kepthx2-dev08-request.md: ONE Si64 MPI8 run D vs H as in the work order, plus si128 D vs H.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-physics/src/distributed_real_gamma.rs   # keeping density pass and apply_local_change only
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - docs/design/reference-scf.md   # "Kept HX warm start" section only
  - docs/performance/lanes/RSCF/kepthx2-dev08-request.md
TEST_COMMANDS:
  - cargo test -p dft-engine --test refscf fast_order_
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-physics
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-physics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: keeping density pass has the plain pass's FFT count and no extra full pass; kept-HX results bit-identical to before; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
