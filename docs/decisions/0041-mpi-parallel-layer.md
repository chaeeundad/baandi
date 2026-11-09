# ADR 0041: MPI parallel layer (staged)

Status: Proposed (scope approved by the user on 2026-11-09 as U9 in AGENTS.md; the decomposition and gates below are the design contract to confirm)

## Context and measured need
dev08 comparison against QE 7.6 at QE's best start (docs/performance/qe-vs-baandi-current.md): BAANDI is 2.2x to 19.4x slower; si32 strict is 11.6x at one core and 19.4x at eight cores, so about 1.7x of the 19.4x is parallel efficiency and the rest is single-core cost. QE scales 4.26x on 8 MPI ranks (identical work per rank count after its atomic start), BAANDI 2.55x on 8 threads: 2.77 of 8.5 thread-stacks work and 67.6 percent of the thread time waits. The thread layer cannot leave one node; multi-node runs need MPI. Single-core cost (algebra, step count) remains the larger part of the gap and keeps its own track; MPI is not expected to close it.

## Options
A. Band-parallel H application, basis replicated: each rank applies H to its share of the columns of the Davidson block (FFT, nonlocal, local multiply are independent per column pair) and allgathers the products. Traffic per Davidson step is the block's H columns (about 2 MB at si32), 0.4 GB per run; parallelises the H application (FFT 20 to 37 percent plus nonlocal 8 to 22 percent of the wall) but leaves the dense algebra (about 50 percent in strict) replicated.
B. Plane-wave-row distributed algebra: the basis and products are stored as row blocks; rotate, Gram, projection and the reduced matrix are local GEMMs on the local rows plus an allreduce of small m x m or m x b matrices (about 1300 reductions, 50 to 185 KB each, latency dominated); needs transposes between row blocks and whole columns for the FFT stage (all-to-all of the block, about 17 GB per run at si32).
C. QE-style full G-vector distribution including a distributed 3D FFT (slab or pencil all-to-all): the most scalable, the largest change.
D. k-point parallelism (each rank owns k points): only for k meshes; trivially deterministic; complements A to C.

## Decision (recommended, staged)
M0 (feasibility, no numerical change): `dft-mpi` adapter crate (or rsmpi after the library checks of ADR 0003: license, maintenance, build with Open MPI and Intel MPI), build feature `mpi` off by default, single-rank and mock paths, a smoke test (rank/size, allgatherv, ordered reduction) on dev08.
M1: option A inside the real Davidson (threads inside each rank, MPI across ranks), plus D for k meshes. Gate: pre-registered after M0 from a Mac and dev08 microbenchmark of the allgather; target to be fixed then (for example si32 strict at 8 ranks at least 1.5x faster than 8 threads).
M2: option B for the dense algebra, only if M1 leaves the algebra as the limit (Amdahl check from the dev08 phase timers).
M3: option C only if a measured need (large cells beyond the memory of one rank).

## Determinism contract
- Rank count 1: bit-identical to the threads-only path (a hard gate).
- Fixed rank count and thread count: bit-reproducible (reductions in a fixed order: allgather then a local ordered sum, never an order-undefined MPI_Allreduce on floating-point data; row blocks and band shares are fixed functions of rank count).
- Different rank counts: equal within the accuracy classes (strict 1e-10 Ha, standard 5e-9 Ha), densities not bit-identical; documented in the results (`parallel.mpi_ranks`).

## Scope and non-goals
Real Gamma path, unpolarized fixed occupation first; k meshes through D; spin, finite temperature and forces follow only after M1 passes. No MPI in numerical hot loops other than the planned collective calls; no I/O from non-root ranks; `unsafe` only in the adapter crate with documented preconditions.

## Validation
Single rank versus threads (bit-identical); two to eight ranks versus single rank (accuracy classes, counts); perturbed-start distributions as for MKL-R1 and MGS-B; scaling measured on dev08 only (Mac observations excluded).

## Open items
MPI environment on dev08 (implementation, version, launcher, interconnect), MPI on the Mac for development (not installed), crate or own adapter, hybrid defaults (ranks times threads) and pinning rules.
