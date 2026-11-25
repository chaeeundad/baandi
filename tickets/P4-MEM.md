TASK_ID: P4-MEM
GOAL: |
  BAANDI uses 2-4x QE's memory (precise campaign 2026-11-26, peak RSS per rank): MPI8 si128 834 vs 269 MB, si216 2083 vs 597 MB; one core si128 3998 vs 1589 MB. Kept-HX slabs are ~25% (si216 ~500 MB/rank); the rest (~1.4-2.9 GB at si216) has no breakdown. Find where it goes and cut it without slowing the run.
READ FIRST: docs/performance/AGENT-BRIEF.md; docs/design/reference-scf.md (Kept HX, Start, workspaces); results.json "peak_payload_bytes_per_rank" / memory_estimate fields; ledger rows KEPTHX2, P4 (KEPT_RECOMPUTE). Never read QE source.
TASK: |
  1. Account every large allocation of a reference-driver run per rank (start arrays and atomic block, Davidson workspace at capacity 2n: V, HV, scratch; GEMM tiles/scratch; FFT plans and buffers per lane; kept-HX slabs and HX; Broyden history; field arrays; verification copies; output/gather buffers) with sizes for si128/si216 MPI8, and record them in results.json (a memory breakdown object, measured from the actual allocations, not estimates).
  2. Remove what is not needed: start arrays freed after the first map, duplicate copies of V/HV or density, history stored at full grid when the owned half sphere suffices, oversize scratch, FFT buffers per lane that can be shared. Behind switches if it changes arithmetic order; bit-identical changes may be proposed default on.
  3. Mac checks only. docs/performance/lanes/RSCF/p4-mem-dev08-request.md: ONE run per arm si128/si216 MPI8 with peak RSS and the new breakdown, plus wall.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/src/refstart.rs   # freeing start arrays only
  - crates/dft-solvers/src/ref_davidson/   # workspace sizing only
  - crates/dft-cli/src/run.rs   # memory breakdown JSON only
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/lanes/RSCF/p4-mem-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-solvers
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-solvers -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: measured memory breakdown in results.json; reductions with unchanged results; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
