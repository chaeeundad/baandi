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

## Revision 2026-11-09 after the dev08 collective measurements and the Amdahl check
Measured (Open MPI 5.0.8, shared memory, one socket, us per call, P = 2 / 4 / 8): Allgatherv 11877 x 19 columns 82 / 254 / 400, 11877 x 68 columns 303 / 952 / 1982; Allreduce 408 x 68 44 / 80 / 129; allgather-plus-rank-ordered-sum (`ordered_sum`) 49 / 204 / 497 (1.1x / 2.6x / 3.8x the Allreduce; 6.3x across two sockets at P = 8); Alltoall of a 2.99 MB si32 FFT grid 115 / 72 / 51. Consequences:
1. Option A alone (band-parallel H, replicated algebra) is bounded by Amdahl: the H application is about 40 to 55 percent of the strict wall, so at most about 1.8x at 8 ranks, less than the 2.55x that threads already give (threads also parallelise the algebra). A is therefore not a stand-alone milestone.
2. M1 is redefined as the layout B-lite: the basis and products are stored as plane-wave row blocks, the algebra (rotation, Gram, projection, reduced matrix) is local GEMM on the local rows plus small reductions, and the H application is distributed by column pairs with an all-to-all transpose of the correction block between the two layouts per Davidson step (volume about 2 x 1.8 MB per step, about 0.1 ms per transpose, a few percent of the run); the FFT stays local per column pair (no distributed 3D FFT).
3. Reductions: `ordered_sum` costs 3.8x an Allreduce at P = 8 and would be about 19 percent of an ideally scaled si32 strict wall; the deterministic replacement is `World::tree_sum_f64` (fixed binomial-tree association, 2 log2(P) small messages), which keeps bit-reproducibility for a fixed rank count; its cost is to be measured against Allreduce; the reduced-matrix eigenproblem stays replicated per rank (about 5 percent of the T=1 wall, serial: it becomes the Amdahl floor and is addressed by dsyevd and by the subspace size, not by MPI).
4. Keep ranks inside one socket while P <= 8 (cross-socket reductions cost 1.15 to 1.9x).
5. Multi-node evidence is not obtainable on dev08 (shared memory and TCP only); M3 stays a design item.

## Reduction rule (decided 2026-11-09 from the dev08 Rust measurements)
Rust `tree_sum_f64` against `ordered_sum_f64` (us per call, 408 x 68 block): P=2 94 against 78, P=4 154 against 247, P=8 228 against 454; 68 x 68: 17 / 13, 29 / 38, 45 / 87; 16 x 68: 6.7 / 4.5, 11.8 / 12.3, 14.1 / 28.3; MPI_Allreduce (C, not deterministic by contract) 43.6 / 90 / 127 for 408 x 68. The tree is the single deterministic reduction for every rank count (its P=2 penalty is 16 us per 408 x 68 call, about 20 ms per run); `ordered_sum_f64` stays as a reference and test oracle. The two differ by 1 ulp at P=8 (identical at P=2 and 4), so the rule never switches between them by size or rank count: for a fixed rank count the result is one fixed function. Measured at the new default (si32 strict, MKL, unperturbed): T=8 8.179 s (QE 8 ranks 0.827 s), SCF 10, Davidson 91, H columns 3397; the Davidson count against QE falls from 12.1x to 3.6x, so the algebra and per-step cost remain the larger gap.

## Determinism contract
- Rank count 1: bit-identical to the threads-only path (a hard gate).
- Fixed rank count and thread count: bit-reproducible (reductions in a fixed order: allgather then a local ordered sum, never an order-undefined MPI_Allreduce on floating-point data; row blocks and band shares are fixed functions of rank count).
- Different rank counts: equal within the accuracy classes (strict 1e-10 Ha, standard 5e-9 Ha), densities not bit-identical; documented in the results (`parallel.mpi_ranks`).

## Scope and non-goals
Real Gamma path, unpolarized fixed occupation first; k meshes through D; spin, finite temperature and forces follow only after M1 passes. No MPI in numerical hot loops other than the planned collective calls; no I/O from non-root ranks; `unsafe` only in the adapter crate with documented preconditions.

## Validation
Single rank versus threads (bit-identical); two to eight ranks versus single rank (accuracy classes, counts); perturbed-start distributions as for MKL-R1 and MGS-B; scaling measured on dev08 only (Mac observations excluded).

## Decisions of 2026-11-09 (user)
- No MPI installation on the development Mac; all MPI measurements and runs happen on dev08 through the master session (neutral machine). The Mac only type-checks.
- Therefore the adapter is an own thin FFI crate `dft-mpi` (declarations only, linked at build time from `MPI_ROOT`/`mpicc -show` like `dft-mkl` with `MKLROOT`), not rsmpi, whose build script needs a working `mpicc` even for `cargo check`.
- Multi-core comparisons against QE so far are MPI (QE) against threads (BAANDI); they become like for like only after M1.

## dev08 MPI facts (relayed 2026-11-09)
Open MPI 5.0.8, user-built at ~/fq/deps (2026-10-19) with the ~/fq/tc gcc; `mpicc -show` = gcc with -I/-L/-rpath of ~/fq/deps and -lmpi; no Intel MPI or MPICH. BTL components: self, sm, tcp only; PML ob1 (plus cm, monitoring, v); no UCX, no InfiniBand component. Topology: 2 x Xeon Gold 6230, 20 cores per socket, 2 threads per core, 80 logical CPUs, two NUMA nodes (node0 = CPUs 0 to 19 and 40 to 59, node1 = 20 to 39 and 60 to 79); numactl not installed. Consequence: M0 to M2 (one node, shared-memory transport) are measurable on dev08; M3 and any multi-node scaling are not (a TCP number over Ethernet is not representative of a cluster interconnect) and must not be reported as scaling evidence from this machine.

## Open items
MPI environment on dev08 (implementation, version, launcher, interconnect), hybrid defaults (ranks times threads) and pinning rules.
