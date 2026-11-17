# dev08 request: replicated reduced eigenproblem and tree sum on the lean MPI path (REDUCED-EIGEN-LEAN, 2026-11-18)

Ledger section REDUCED-EIGEN-LEAN (rows REL-FINDING, REL-CHECKS, REL-SLICES, REL-TREE-RHD). Every job below is sized under 10 minutes. Load below 0.30, `--bind-to core`, one socket (CPUs 0-7). Report binary sha256, raw output, load before and after.

## Stage 0a: MKL microbenchmark of a sliced reduced eigensolve (REL-SLICES)
```
MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-mkl --features mkl --example evd_slices
sha256sum target/release/examples/evd_slices
taskset -c 2 target/release/examples/evd_slices 0.34 20            # one busy core
for c in 0 1 2 3 4 5 6 7; do taskset -c $c target/release/examples/evd_slices 0.34 20 > ~/fastbench/evd_slices-8x-$c.txt & done; wait   # 8 busy cores (shared node, as at np=8)
```
Decision (pre-registered): implement `BAANDI_M2_EIGEN=slices` only if `sliced-mrrr` or `sliced-stein` at P=8, plus the m x k allgatherv time at np=8 (relayed collective timers), is at most 0.80 x `dsyevd` at m = 204 and 408 in both matrix models, with orthogonality <= 1e-13; otherwise the distributed-EVD lever is closed for orders below 1000.

## Stage 0b: halving/doubling tree sum (REL-TREE-RHD)
```
cargo build --release -p dft-mpi --features mpi --example mpi_tree_bench
sha256sum target/release/examples/mpi_tree_bench
for P in 2 4 8; do mpirun -np $P --bind-to core --cpu-set 0-7 target/release/examples/mpi_tree_bench; done
```
The program aborts if binomial, rd, rhd and rhd with threshold 0 ever differ in one bit on a power-of-two rank count. Report the per-payload lines (columns binomial / rd / rhd / rhd-forced, worst rank, us per call). The crossover payload where `rhd-forced` beats `rd` calibrates `dft_mpi::HALVING_MIN_DOUBLES` (now 2048).

## Stage 1: 10-step fine-solve screens at np 4 and 8 (the eecd57f finding protocol)
Same binary (`cargo build --release -p dft-cli --features mpi,mkl`), si32 and si64 standard f1e-3 seed 601, `scf.max_iterations = 10`, `BAANDI_PHASE_PROFILE=on`, three alternating pairs per cell:
- arm A: main at this branch's merge (contains REL-CHECKS, always on; compare its `solve.reduced_eigen` against the eecd57f rows si32 0.0952 / 0.1155 s at np 4 / 8: expected -2 to -6 %, E);
- arm B: A plus `BAANDI_MPI_TREE=rhd` (expected: fine `coll.tree_sum` seconds down by up to 20 to 35 % of their latency/bandwidth share at np=8, energy and density bits identical to A).
```
for P in 4 8; do for ARM in A B A B A B; do
  env BAANDI_PHASE_PROFILE=on $( [ $ARM = B ] && echo BAANDI_MPI_TREE=rhd ) \
    mpirun -np $P --bind-to core --cpu-set 0-7 target/release/dft-cli run <si32-10step.toml> --output-dir ~/fastbench/rel-si32-np$P-$ARM-$RANDOM
done; done          # then the same for si64
```
Rows to report per cell: fine `solve.reduced_eigen`, `coll.tree_sum` (seconds and calls), `solve.finish`, `solve.real_state`, wall, energy bits, density SHA-256.

## Optional diagnostic (not an adoption arm)
The np4 -> np8 per-unit-work slowdown hypothesis (AVX-512 / all-core frequency): one pair of the stage-1 si32 np=8 cell with `MKL_ENABLE_INSTRUCTIONS=AVX2` (rounding class); compare `solve.reduced_eigen` and `solve.rotate` against arm A.
