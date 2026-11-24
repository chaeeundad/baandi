TASK_ID: P3-FFT2
GOAL: |
  Cut the data movement of the distributed FFT used by the reference driver's H application. dev08 si128 MPI8 at 9e99bd3, one profiled run (rank max): FFT forward+inverse 10.2 s vs QE fftw 6.7 s. The 1D kernels are only ~4.9 s (MKL DFTI would save ~1 s of that but loses on small lengths; not adopted). The rest is movement: alltoall fwd 0.86 + inv 0.98, middle pack 0.69 + unpack 0.91, slab gather 0.53 + scatter 0.70, transpose pack 0.47 + unpack 0.40, potential multiply 0.56. P3-FFT's direct ORC pack and sparse clear were bit-identical but showed no measurable gain.
READ FIRST: docs/performance/AGENT-BRIEF.md; docs/performance/attempts-ledger.md entries A26, A27, A29, A11, PFFT-3, PFFT-4, P39A-1..4 (the mirror x-slab pair plan was measured 0.87-0.91 on the Mac and blocked only by an accuracy gate that P39A-4 later resolved), M2_FFT, P3-FFT; docs/research/qe-study/fft-call-structure.md, qe-hpsi-and-fft.md, qe-gamma-real-path.md; docs/design/reference-scf.md "H-apply FFT cost accounting". Never read QE source.
TASK: |
  1. From the accounting, list for one band pair and one direction every pass over memory (reads/writes of the stick, middle, slab and plane buffers, the alltoall send/recv buffers) and how many of them a minimal one-transpose scheme needs (sticks -> one alltoall -> planes, transforms fused with the pack/unpack of the neighbouring stage, potential multiply fused into the plane pass). Record the comparison in the accounting section.
  2. Implement the largest removable part behind a switch (default off), e.g. fuse the middle pack/unpack into the 1D transform of that axis (transform in place in the receive buffer layout), fuse slab gather/scatter with the plane FFT and the potential multiply (as PFFT-3 did serially), or batch several band pairs per alltoall to cut latency and pack passes. Prefer a change whose operation order keeps results bit-identical; otherwise it is rounding class.
  3. Tests (Mac only: cargo check, touched crates' tests): on/off transforms on cubic, skewed, elongated, odd and prime grids and simulated rank counts 1/2/3/8; H apply bits or rounding-class agreement; reference SCF fixtures unchanged in counts.
  4. docs/performance/lanes/RSCF/fft2-dev08-request.md: ONE profiled run per arm on si64/si128/elongated/rskew MPI8 with the FFT split (development protocol, no repeated runs).
OWNED_FILES:
  - crates/dft-solvers/src/distributed_fft.rs
  - crates/dft-physics/src/distributed_real_gamma.rs   # FFT call sites only (P3-QECLOCK adds counters in this file in parallel; keep edits local)
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "H-apply FFT cost accounting" section only
  - docs/performance/lanes/RSCF/fft2-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-physics
  - cargo test -p dft-engine --test refscf
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-physics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: pass accounting written; one switchable movement reduction with on/off tests on several grids and rank counts; old path untouched; dev08 single-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
