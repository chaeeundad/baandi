TASK_ID: P3-PERMAP
GOAL: |
  Two per-map items on dev08 si128 MPI8 (verify off, 15 maps, P3s at e5ec8ae): (a) the 1D FFT kernels are ~4.9 s; MKL DFTI (BAANDI_M2_FFT=mkl) was 20-25% faster on the long axes of si64/si128 but slower on short lengths (elongated cell +8%), so the global switch was not adopted (ledger M2_FFT, P3e). (b) per-map work outside solve and density is ~0.66 s vs QE ~0.27 s: energy_correction 0.25, mix 0.10, density_transform 0.07, energy 0.20, potential 0.04.
READ FIRST: docs/performance/AGENT-BRIEF.md; attempts-ledger.md rows M2_FFT, P3e screen (M2_FFT=mkl), P3q (LEAN_ENERGY), P3-FFT, P3-FFT2; docs/design/reference-scf.md "H-apply FFT cost accounting" and "Per-map energy cost". Never read QE source.
TASK: |
  1. FFT kernel selection by length (switch FFT_KERNEL_BY_LENGTH, default off): per 1D transform length (and batch size) choose rustfft or MKL DFTI by a rule that is measured, not tuned to a material: a one-time per-run calibration (time both kernels on the actual lengths/batches at plan creation, pick the faster, record the choice) or a documented length threshold derived from a micro-benchmark over many lengths. Pure builds keep rustfft. Record choices in results.json. Results are rounding class; on/off transform tests over cubic, skewed, elongated, odd and prime grids. Do not touch plan creation/sharing (P3-SETUP4 owns that part of distributed_fft.rs; keep edits to the kernel dispatch).
  2. Per-map work: account energy_correction (the descf integral and fields of the mixed density), mix and density_transform; reuse fields already computed (LEAN_ENERGY caches) and fuse the integrals into existing passes and collectives; behind a switch if not bit-identical.
  3. Mac checks only: cargo check --workspace --tests; cargo test -p dft-solvers; cargo test -p dft-engine --test refscf; cargo test -p dft-numerics; fmt; clippy -D warnings. Grep tests for any switch before proposing a default.
  4. docs/performance/lanes/RSCF/permap-dev08-request.md: ONE profiled run per arm, four cases MPI8, with FFT kernel split, kernel choices and per-map phase seconds.
OWNED_FILES:
  - crates/dft-solvers/src/distributed_fft.rs   # 1D kernel dispatch only
  - crates/dft-numerics/src/fused_pair_plan.rs   # kernel dispatch only
  - crates/dft-engine/src/refscf/driver.rs   # energy correction / mixing per-map code only
  - crates/dft-engine/src/refscf/broyden.rs
  - crates/dft-engine/src/refscf/metric.rs
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-cli/src/run.rs   # kernel-choice JSON only
  - docs/performance/lanes/RSCF/permap-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-numerics -p dft-engine --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: measured length-based kernel choice behind a switch with transform tests; per-map correction/mix work reduced; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
