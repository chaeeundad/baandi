TASK_ID: P4-DENSKEEP
GOAL: |
  With REFSCF_KEPT_HX on (default), the density step costs more than QE's sum_band: dev08 si128 MPI8 profiled snapshot at fa4cfad baandi_density 2.327 s / 16 vs QE sum_band 1.78 s; before KEPT_HX it was ~1.83-2.12 s. The kept real-space slabs (185 MB/rank si128, ~500 MB si216, about 25% of RSS) are written during the density pass. Remove the remaining density overhead and the memory cost while keeping the kept-HX gain (solve -2.1 s on si128).
READ FIRST: docs/design/reference-scf.md "Kept HX warm start"; crates/dft-engine/src/refscf/kept_hx.rs; ledger rows P3-KEPTHX, KEPTHX2, Work-order screen; docs/performance/AGENT-BRIEF.md. Never read QE source.
TASK: |
  1. Measure on the Mac fixtures (counters) where the extra density time goes with keeping on: slab store bandwidth (bytes written per map), layout (strided vs contiguous), first-touch page faults of large slabs, cache pollution of the density accumulation.
  2. Options behind switches (default off): store slabs in f32 (rounding: dV X computed in f64 from f32 values; quantify the starting-HX error and its effect on Davidson counts and delivered energy; must stay within the strict class), store only the slab part that apply_local_change reads, non-temporal/streaming stores, pre-faulted reusable pages, or recompute instead of store where cheaper.
  3. Tests (Mac): fast_order_ kept-HX test passes (with tolerance documented if f32); fixtures' delivered energies within 1e-10 Ha of the f64 path; counts unchanged or reported.
  4. docs/performance/lanes/RSCF/p4-denskeep-dev08-request.md: ONE run per arm si128 and si216 MPI8 verify off: density s, solve s, wall, peak RSS, E - E_strict.
OWNED_FILES:
  - crates/dft-engine/src/refscf/kept_hx.rs
  - crates/dft-engine/src/refscf/map.rs   # density/keeping call only
  - crates/dft-physics/src/distributed_real_gamma.rs   # keeping density pass and apply_local_change only
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Kept HX warm start" section only
  - docs/performance/lanes/RSCF/p4-denskeep-dev08-request.md
TEST_COMMANDS:
  - cargo test -p dft-engine --test refscf fast_order_
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-physics
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-physics -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: overhead source measured; at least one switchable option that cuts density overhead or slab memory with accuracy within the strict class; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
