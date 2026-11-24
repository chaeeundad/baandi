TASK_ID: P3-FFT
GOAL: |
  Cut the FFT part of the reference driver's H application. dev08 si128 MPI8 at 66d45e69 (rank max): FFT stages without the potential multiply 11.16 s vs QE fftw 6.70 s (QE vloc_psi 7.98 total) for a comparable number of transforms (BAANDI 14204 H columns, two real bands per complex FFT; QE 13592 fftw calls). BAANDI detail: inverse 5.80 / forward 4.95; local 2D FFTs 2.12 + 2.14; alltoall 1.16 + 0.88; transpose pack/unpack 1.37 + 1.10; grid pack/unpack (ORC rows to the FFT grid and back, mirror) 1.49; multiply 0.43. si64 has the same shape (FFT 2.29 vs 1.33). So about 4 s of si128's gap is data movement around the transforms (pack/unpack/transposes) and about 1-2 s may be the 1D/2D transform speed itself.
READ FIRST: docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md (search "FFT", "M2_FFT", "alltoall", "pack", "fused"; several FFT ideas were tried on the old path, do not repeat a rejected one without new evidence), docs/performance/backlog.md, docs/research/qe-study/fft-call-structure.md, qe-hpsi-and-fft.md, qe-gamma-real-path.md. Never read or request QE source.
TASK: |
  1. Account the bytes moved per band pair in each stage (grid pack/unpack, transpose pack/unpack, alltoall payload, 2D/1D transforms) as the code does it now, against the minimum (one scatter into the sticks/columns that hold nonzero G, transforms only on nonzero columns/planes, one transpose each way). Write it into docs/design/reference-scf.md section "H-apply FFT cost accounting".
  2. Remove the excess behind switches (default off), e.g.: fuse grid pack with the first 1D transform (transform only the columns that hold G vectors, as the sphere occupies a fraction of the grid), fuse transpose pack/unpack with the neighbouring transform, avoid zero-filling the full grid when only sticks are non-zero, process several band pairs per alltoall to amortize latency, reuse plans/buffers. Keep the result identical within rounding; prefer bit-identical where the operation order is unchanged.
  3. Do not touch crates/dft-numerics/src/real.rs (owned by P3-GRAM in parallel).
  4. Mac checks only (user rule 2026-11-24): cargo check --workspace --tests; cargo test -p dft-solvers; cargo test -p dft-physics; cargo test -p dft-engine --test refscf; fmt; clippy -D warnings on touched crates. Add tests comparing transforms with the switch on/off on several grids (cubic, skewed, elongated, odd sizes) and rank counts (simulated).
  5. Write docs/performance/lanes/RSCF/fft-dev08-request.md: untimed screen with the apply_h split (BAANDI_DAVIDSON_STAGES=on) on si64/si128/elongated/rskew MPI8, then A/A + crossed A/B per switch if not slower.
OWNED_FILES:
  - crates/dft-solvers/src/distributed_fft.rs
  - crates/dft-solvers/src/distributed_fft/   # if it exists
  - crates/dft-physics/src/distributed_real_gamma.rs   # FFT call sites and pack/unpack only
  - crates/dft-numerics/src/switches.rs   # append new switches at the end of the list (P3-GRAM also appends one)
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "H-apply FFT cost accounting" section only
  - docs/performance/lanes/RSCF/fft-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-physics
  - cargo test -p dft-engine --test refscf
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-physics --all-targets -- -D warnings
  - (orchestrator, on dev08) scripts/check.sh full, screen, A/B
ACCEPTANCE_CRITERIA: accounting written; switches default off with on/off transform tests on several grids and simulated rank counts; old path bit-identical; dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
