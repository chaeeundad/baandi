TASK_ID: P3-A2A
GOAL: |
  H application is now the dominant reference-Davidson stage: dev08 si128 MPI8 15.9 s of a 28 s stage sum (one profiled run at ec553c93). Its collectives: alltoall 206 calls, 4134 MB per rank (about 20 MB per call), plus 103 tree sums of 58 MB; collective time about 2.0 s; FFT movement stages (middle pack/unpack, slab gather/scatter, transposes) about 3.5 s; 1D kernels about 4.9 s. QE's h_psi on the same cell is 11.25 s for 82 calls. Cut the bytes moved and the time spent waiting.
READ FIRST: docs/performance/AGENT-BRIEF.md; docs/performance/attempts-ledger.md rows MPI-M2D, MPI-SCALING-GAP-STICKS (mirror sticks, adopted), MPI-REMAINING-ASYNC-ALLTOALL (async overlap rejected by analysis at si32 sizes: re-check the premise at si128 message sizes before reopening), M2-OWNED-FIELDS, P3-FFT, P3-FFT2, P3-QECLOCK, P3i; docs/design/reference-scf.md "H-apply FFT cost accounting"; docs/research/qe-study/fft-call-structure.md and qe-hpsi-and-fft.md. Never read QE source.
TASK: |
  1. Account the alltoall payload per band pair and direction as sent now (including padding, zero columns, both halves of mirror pairs, complex vs real packing) against the minimum: only the G-sphere sticks needed after the Gamma real-pair packing, in the narrowest type. Also count calls per H application and the per-call latency share. Write it into the accounting section.
  2. Implement the largest removable part behind a switch (default off): e.g. send only nonzero sticks (sphere) instead of full columns; larger batches of band pairs per alltoall when the H block is wide; or, with new evidence at si128 sizes, overlap the alltoall of batch k with the transforms of batch k+1 (nonblocking collectives through the dft-mpi wrapper only if it supports them; otherwise document what is missing).
  3. Tests (Mac, simulated ranks 1/2/3/8 and several grids incl. odd and prime sizes): transforms and H application equal on/off (bit-identical where the arithmetic order is unchanged), reference SCF fixtures unchanged in counts.
  4. docs/performance/lanes/RSCF/a2a-dev08-request.md: ONE profiled run per arm on si64/si128/elongated/rskew MPI8 with alltoall bytes/calls/seconds and the FFT split.
OWNED_FILES:
  - crates/dft-solvers/src/distributed_fft.rs
  - crates/dft-physics/src/distributed_real_gamma.rs   # FFT call sites and batching only
  - crates/dft-mpi/src/   # only if a nonblocking alltoall wrapper is added; MPI behind the feature, single-rank path unchanged
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "H-apply FFT cost accounting" section only
  - docs/performance/lanes/RSCF/a2a-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-physics
  - cargo test -p dft-mpi
  - cargo test -p dft-engine --test refscf
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-physics -p dft-mpi --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: payload accounting written; one switchable reduction with on/off tests over grids and simulated ranks; old path untouched; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
