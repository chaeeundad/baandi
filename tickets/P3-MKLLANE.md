TASK_ID: P3-MKLLANE
GOAL: |
  Twice now (P3e M2_FFT=mkl and the work-order si128 FFT_KERNEL_BY_LENGTH run, dev08 MPI8) the MKL 1D kernels were 16-40% faster than rustfft (si128 kernels -0.9 s), but the stages AROUND them got slower by ~0.5 s: slab_gather 0.520 -> 0.625, slab_scatter 0.690 -> 0.836, middle pack/unpack 0.613/0.760 -> 0.685/0.844, potential multiply 0.455 -> 0.547 s. Those stages do not call MKL. Find out why and remove the penalty so the kernel gain shows in wall.
READ FIRST: docs/performance/AGENT-BRIEF.md; attempts-ledger rows M2_FFT, P3e, Work-order screen; crates/dft-solvers/src/distributed_fft.rs and crates/dft-numerics/src/fused_pair_plan.rs (lane dispatch, MKL descriptors); crates/dft-mkl/src/ (DFTI adapter: threading layer, descriptor settings, in/out of place, any internal copies or allocations). Never read QE source.
TASK: |
  Hypotheses to check in code and with counters, in this order: (a) MKL threading: is the sequential layer or mkl_set_num_threads(1)/DFTI_THREAD_LIMIT 1 enforced per descriptor? With 8 ranks x several MKL threads the surrounding single-threaded stages lose their cores/cache; (b) the MKL lane path copies lanes into a separate contiguous buffer or uses a different layout than rustfft, leaving the next stage with cold or strided data; (c) descriptor commit or workspace allocation inside the timed region of the neighbours; (d) denormals/FTZ or AVX frequency effects. Fix what is found behind the existing switches; add counters that make the cause visible (bytes copied, threads used by MKL per call).
  Mac: cargo check including `--features dft-numerics/mkl`, touched crates' tests, fmt, clippy (with and without mkl feature). docs/performance/lanes/RSCF/mkllane-dev08-request.md: ONE profiled si128 MPI8 run off vs on, plus the MKL thread count actually used.
OWNED_FILES:
  - crates/dft-mkl/src/
  - crates/dft-numerics/src/fused_pair_plan.rs
  - crates/dft-solvers/src/distributed_fft.rs   # MKL lane dispatch only
  - docs/performance/lanes/RSCF/mkllane-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo check -p dft-solvers -p dft-numerics -p dft-mkl --tests --features dft-numerics/mkl
  - cargo test -p dft-solvers
  - cargo test -p dft-numerics
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-numerics -p dft-mkl --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: cause identified with evidence; fix behind the existing switches; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
