# MPI reduced eigenproblem: the next ceiling of the distributed Davidson

Work order: docs/performance/requests/work-order-beat-qe-2026-11-14.md, phase E2 (the distributed counterpart of B4). MPI against MPI only; no threads number appears here.
Evidence tags: (R) relayed from dev08, (O) Mac observation (pure build, faer EVD, simulated ranks = threads of one shared, loaded Mac; single runs), (C) counted from code or from exact work counters, (E) estimate.
Status 2026-11-15: analysis; two switches implemented, both OFF by default (`BAANDI_M2_EIGEN=root`, `BAANDI_M2_B5=on`), plus the order distribution in results.json `parallel.reduced_solves`; dev08 measurement pending (request in section 9).

## 1. Starting numbers (R, ledger M2-ENTRY-COST-AB-DEV08 and the relay of 2026-11-15)
MKL, one thread per rank, f1e-3 seed 601, entry-cost switches ORTHO/DGEMM/LEAN on. `parallel.stages` reduced_eigen:

| case | np | reduced_eigen s | wall s | share |
|---|---:|---:|---:|---:|
| si32 strict | 8 | 1.37 | 10.0 | 13.7 % |
| si64 strict | 2 | 5.45 | 73.6 | 7.4 % |
| si64 strict | 8 | 7.20 | 47.7 | 15.1 % |

No entry-cost switch changes it (it is not touched by them), and it does not shrink with ranks: the same work takes 32 % longer per rank at np=8 than at np=2.

## 2. What the distributed solver solves (C)
`solve_core` (crates/dft-solvers/src/distributed_davidson.rs) does exactly one dense reduced solve per Davidson iteration, on every rank, all m pairs (`RealEighWorkspace::solve_into`: MKL `dsyevd` on MKL builds, faer otherwise), then uses only the lowest `bands` coefficient columns. The order m starts at `bands` after a (contraction) restart and grows by the number of unconverged roots c per iteration up to the capacity 6 x bands (si32 408, si64 792); a contraction restart returns to `bands`. Solves per SCF step = Davidson iterations of the step (the step's final iteration included).

Distribution, Mac analogues of the dev08 cells (same structure and cutoff as dev08 si32/si64, the harness's own f1e-3 displacements; counts are trajectory properties, not machine properties):

| case (O/C) | Davidson / solves per run | order mean | p90 | max | sum m^3 | single-column appends |
|---|---:|---:|---:|---:|---:|---:|
| si32 strict, np=2 and 4 | 279 / 293 | 224.6 | 393 | 408 | 6.51e9 | 74 (25 %) |
| si64 strict, np=4 | 214 / 228 | 460.8 | 743 | 787 | 3.87e10 | 42 (18 %) |

The cubic term is dominated by the large orders: at si64 the solves with m >= 500 carry most of the sum (the serial solver: 82 percent, davidson-algebra-gap-analysis.md 3.3).

## 3. Cost model and where each rank spends its time
- `dsyevd` with eigenvectors: tridiagonalisation 4/3 m^3 (half of it BLAS2), divide and conquer on the tridiagonal (up to about 4/3 m^3 before deflation), back-transformation 2 m^3; measured, not counted: MKL-EVD microbenchmark 11.7 ms at m = 408 on one dev08 core = 0.17 ns per m^3 (R).
- dev08 si64 strict (E): 161 Davidson iterations, about 172 solves; at the Mac analogue's sum m^3 scaled to 172 solves (2.9e10) the relayed 5.45 s is 0.19 ns per m^3 at np=2 (consistent with the microbenchmark) and 7.20 s is 0.25 ns per m^3 at np=8. Per solve: about 32 ms at np=2, 42 ms at np=8.
- The stage contains no communication (timer around the solve only). Every rank computes identical bits (rank energies equal rank 0, R), so the time per rank is the time of the critical path.
- Why it grows with ranks (E): the same single-thread solve runs concurrently on every rank. Its working set (matrix plus `dsyevd` workspace, about 3 m^2 doubles = 15 MB at m = 792) times 8 ranks exceeds one socket's L3, so the concurrent copies compete for L3 and memory bandwidth; the replicated non-Davidson stages show the same doubling mechanism beyond 8 ranks (ledger MPI-SCALING-BEYOND-8). Mac O: si32 np=4 reduced_eigen 1.400 s replicated against 1.194 s when only one rank solves (section 6).

## 4. Serial structure that the distributed solver lacks
| serial item | distributed solver | effect on this stage |
|---|---|---|
| B5 bordered secular update when one column was appended (default ON since B5-PART-C) | missing | 18 to 25 % of the solves qualify (table 2); bordered solve about 28 % cheaper at m = 408 (`bordered_timing`, ALGEBRA-STREAM-B5-PAR: 14.3 against 10.3 ms) -> implemented as `BAANDI_M2_B5=on` |
| B1 selected roots / narrow iterations | missing | changes the rotation and residual stages, not the solve (the serial solve is also all pairs); separate E2 item |
| selected-root range solve (`dsyevr`) | not used | rejected by MKL-EVD: slower than `dsyevd` below order/b about 9; here order/b <= 6. Not repeated |
| capacity 6 x bands | same | see section 5 (smaller capacity does not pay) |

Only the lowest `bands` pairs are used, but with dense LAPACK drivers the measured fastest way to get them at these ratios is still all pairs (MKL-EVD).

## 5. Alternatives to the replicated dense solve
- Solve on one rank and broadcast (`BAANDI_M2_EIGEN=root`): the critical path is still one solve, but only one copy runs, so the contention of section 3 disappears. Broadcast per iteration: 1 + bands + m x bands doubles (si64 mean 61 k doubles = 0.49 MB). Bit-identical to the replicated solve (same matrix, same single-thread solver).
- Disjoint column blocks of the eigenvector problem: impossible with dense `dsyevd`. Possible with a distributed dense eigensolver (ScaLAPACK `pdsyevd`/`pdsyevr` from MKL with the Open MPI BLACS, or ELPA), which splits tridiagonalisation and back-transformation over a process grid; bits are fixed for a fixed grid (ADR 0041 holds). This is the only option that shortens the critical path itself; it needs a library adoption (ADR 0003 checks, docs/dependencies.md) and cannot run on the Mac, so the first step is a dev08 microbenchmark (orders 400, 600, 792; P = 2, 4, 8) against `dsyevd`.
- Smaller capacity (fewer and smaller solves, more iterations), Mac counts (O, si32 strict np=2, entry-cost stack):

| capacity per band | SCF / Davidson / H | sum m^3 | reduced_eigen s | other Davidson stages s | wall s |
|---|---|---:|---:|---:|---:|
| 6 | 16 / 279 / 6732 | 6.51e9 | 1.20 | 10.63 | 13.93 |
| 4 | 16 / 355 / 7725 | 2.85e9 | 0.74 | 12.36 | 15.31 |
| 3 | 17 / 464 / 9012 | 1.64e9 | 0.48 | 12.27 | 15.02 |

  At si64 strict np=8 (E): capacity 4 saves about 0.39 x 7.2 = 2.8 s of eigen but adds about 16 percent to the other distributed Davidson stages (about 35 s): +5.6 s. Net loss; not proposed (the premise of A13 holds on the MPI path too).
- Locking converged roots out of the reduced problem (smaller m): an algorithm change with its own history (A56 gated, P24A not adopted); not proposed here.

## 6. Collectives: where the tree_sum calls come from (C)
Per Davidson iteration with `BAANDI_M2_ORTHO=serial`, panel width p (16 on MKL builds), c appended corrections:

| site | tree_sum calls |
|---|---:|
| reduced Gram (new border) | 1 |
| residual norms | 1 (+1 per certificate attempt) |
| orthonormalize: normalisation, basis projection | 2 |
| panel Gram-Schmidt: cross-panel projections | 2 (ceil(c/p) - 1) |
| panel Gram-Schmidt: per column two dot reductions and one norm | 3c - 2 ceil(c/p) |
| second pass when DGKS fails | the three rows above again |
| H apply (nonlocal P^T X) | 1 per apply |
| restart defect | 1 |
| per SCF step: orthonormalization of the warm start (bands columns) | about 3 x bands |

Measured totals (O, exact counts): si32 23507 per rank (about 84 per Davidson iteration), si64 46722 per rank (about 205 per solve); by the table about 85 percent are the per-column in-panel reductions and about 14 percent the warm-start orthonormalization (si64: 17 steps x 3 x 132 = 6.7 k).

Fusion: no material exact fusion exists. Each column's second dot reduction depends on the first, its norm on the second, and the next column's first reduction on the kept column (modified Gram-Schmidt order); the remaining independent pairs (residual norms with the correction norms; normalisation with the projection Gram) save one or two calls per iteration (under 1 percent). One Allreduce per panel would be CholeskyQR on the panel Gram, which drops near-dependent corrections (A52, MPI-M1-GATE1): excluded. What can be reduced:
- Latency per call, exact: `tree_sum_f64` is a binomial reduce to rank 0 followed by a broadcast (2 log2 P message rounds). A recursive-doubling exchange in the same pairing computes, on every rank, the same block partial sums in the same association (IEEE addition is commutative, so both partners obtain identical bits), in log2 P rounds: bit-identical to the current tree sum for P a power of two (other P keep the current code). dev08 microbenchmark (wiki/mpi.md): tree sum 228 us against Allreduce 129 us at P = 8 for 408 x 68 values, i.e. the reduce-plus-broadcast pattern costs about twice one exchange pass.
- Count, rounding class: the second-pass norm from the Pythagorean identity, reduced together with the second-pass dots (3 -> 2 reductions per column, about -28 percent of the calls), with an explicit norm reduction when the second-pass coefficients are not small against the first-pass norm.
- Count, rounding class: skip re-orthonormalizing a warm start whose Gram defect is at most 1e-12 (the contraction-restart rule), one Gram reduction instead of 3 x bands calls per SCF step (about -14 percent).

Collectives do not shrink with ranks because the call count is fixed by the trajectory (identical at every P) while the latency per call grows with log2 P (Mac O si32: 20 us per call at np=2, 38 us at np=4).

## 7. Ranked reductions (si64 strict, per rank, E unless marked)
| rank | change (switch) | rounding | np=2 | np=8 | status |
|---|---|---|---:|---:|---|
| 1 | rank 0 solves, broadcast of the lowest pairs (`BAANDI_M2_EIGEN=root`) | exact (bit-identical) | about 0 | 0 to -1.7 s (removes the 32 % contention excess; less if the slowdown is clock frequency, which spinning ranks keep) | implemented, off |
| 2 | bordered secular update on single-column iterations (`BAANDI_M2_B5=on`) | rounding class | -0.2 s | -0.2 to -0.4 s | implemented, off |
| 3 | recursive-doubling tree sum for power-of-two P | exact | -0.1 s | -0.3 to -0.5 s | not done (FFI; no MPI on the Mac) |
| 4 | Pythagorean second-pass norm | rounding class | -0.1 s | -0.2 to -0.3 s | not done |
| 5 | skip re-orthonormalizing an orthonormal warm start | rounding class | -0.1 s | -0.1 to -0.3 s | not done |
| 6 | distributed dense EVD (ScaLAPACK / ELPA) | rounding class, fixed per P | | -3 to -5 s if 2-3x at P = 8 | needs microbenchmark and library adoption |
| - | capacity 4 or 3 per band | trajectory | | net loss (section 5) | not proposed |
| - | range solve `dsyevr` | | | slower (MKL-EVD) | not repeated |

Mac observations of 1 and 2 (O, single runs, load 4.7 to 6.2, work counts identical in every row, energies identical to 12 printed digits):

| case | default | EIGEN=root | B5=on |
|---|---:|---:|---:|
| si32 np=2 reduced_eigen / wall | 1.202 / 13.93 | 1.176 / 13.76 | 1.130 / 13.73 |
| si32 np=4 reduced_eigen / wall | 1.400 / 10.35 | 1.194 / 10.16 | (run disturbed by load: 1.546 / 12.35) |
| si64 np=4 reduced_eigen / wall | 6.737 / 55.09 | root and B5 together: 5.481 / 57.01 (load 2.9 -> 8.7 during the run; every other stage 5 to 12 % slower) | |

B5 fired on 70 of the 74 single-column solves per rank at si32 and on 42 of 42 at si64 (the others: streak bound or a base broken by a restart). The Mac walls are not evidence either way (shared machine, single runs); the stage drops (-15 % root at si32 np=4, -19 % root plus B5 at si64 np=4) are the observation.

## 8. Implementation and tests
- `crates/dft-solvers/src/distributed_davidson.rs`: `DistributedTuning.eigen_root` and `.bordered`; `solve_reduced` (B5 with its own capacity x capacity GEMM workspace, the serial streak bound `MAX_STREAK` and invalidation at every restart); root mode packs status, `values[..bands]` and the first `bands` coefficient columns into one broadcast, and a failure on rank 0 is broadcast as a status so no rank waits forever; `reduced_solve_stats()` (order histogram, sum of m^3, single-column, bordered, broadcast). Non-root ranks in root mode allocate no eigen workspace.
- `crates/dft-numerics/src/switches.rs`: `BAANDI_M2_EIGEN` (replicated, root; default replicated) and `BAANDI_M2_B5` (on, off; default off); `crates/dft-cli/tests/run.rs` switch count 56 -> 58; results.json `parallel.reduced_solves`.
- Tests (Mac, pure build, simulated ranks): `crates/dft-solvers/tests/distributed_davidson.rs` `root_eigen_broadcast_is_bit_identical_to_replicated` (eigenvectors, eigenvalues and work counts bit-identical to the replicated solve at P = 2, 3, 4, 5; entry-cost stack on scattered and clustered spectra, classic order on clustered) and `bordered_reduced_update_stays_in_the_rounding_class` (entry-cost stack, scattered and clustered: eigenvalues within 1e-9 of the dense solve and the serial solver, residuals below 1e-8, ranks bit-equal, reproducible per P, iterations within 10 percent plus 2, B5 plus root bit-identical to B5 replicated, B5 fires); debug times 8.8 and 8.7 s; `crates/dft-engine/tests/distributed_scf.rs` `m2_reduced_eigen_switches_keep_classes` (SCF, cubic and skewed cells, strict and standard, 2 and 3 ranks: root bit-identical energy and density to replicated; B5 plus root against the threads-only run |dE| 2.2e-16 to 4.4e-16 Ha strict, at most 6.6e-11 Ha standard (limits 1e-10 and 5e-9), Davidson 44 vs 42 / 67 vs 67 strict and 15 to 17 vs 15 / 24 vs 24 standard, ranks bit-equal, reproducible); switch default tests in `switches.rs` and `tests/switches.rs`.
- Profile harness: `crates/dft-engine/tests/m2_entry_cost.rs` now has si64, `M2_PROFILE_CAPACITY` and `M2_PROFILE_SERIAL=off`, and prints the reduced-solve distribution.

## 9. Ledger cross-check
- MKL-EVD (range solve slower): respected, no `dsyevr`. Threaded `dsyevd` (thread-dependent bits) stays excluded.
- B5-PART-C / ALGEBRA-STREAM-B5-PAR: the distributed B5 is the adopted serial code path (`bordered.rs`), not a new algorithm; `B5=par` is not used (one thread per rank).
- A13 / A08 (smaller capacity): re-examined on the MPI path because the replicated solve does not divide by P (new premise); the Mac counts confirm the rejection, so no capacity change is proposed.
- A52 / MPI-M1-GATE1 (Cholesky-based block orthonormalisation): not used; the collective fusions above keep per-column dependence decisions.
- A56 / P24A (locking): not proposed.
- MPI-SCALING-BEYOND-8: the contention mechanism of the replicated stages is the same one the root solve removes for this stage.

## 10. dev08 request (MPI against MPI only)
Build: this branch with `--features mkl,mpi` (MKLROOT=/opt/intel/oneapi/mkl/2026.1), one binary, sha256 recorded; MKL_NUM_THREADS=1, one thread per rank, `--map-by core --bind-to core` on socket 0, load gated as in BAANDI-MPI-RANK-SWEEP. Every arm sets the adopted entry-cost switches explicitly (`BAANDI_M2_ORTHO=serial BAANDI_M2_DGEMM=on BAANDI_M2_LEAN=on`) so the result does not depend on whether the build contains the 2026-11-15 default change.
Cells: si32 and si64, strict and standard, f1e-3 seed 601 (si32 also 602, 603), ranks 2, 4, 8 (np=1 is the serial path, not affected).
Arms (same binary, environment A/B, interleaved per seed): 1. reference (new switches unset); 2. `BAANDI_M2_EIGEN=root`; 3. `BAANDI_M2_B5=on`; 4. both.
Report per row: wall; `parallel.stages` (reduced_eigen, reduced_gram, rotate, residual, orthonormalize, apply_total); `parallel.reduced_solves` (solves, mean/p90/max order, sum of m^3, single_column, bordered, broadcast); `parallel.collectives` (bcast calls and seconds in arms 2 and 4); SCF / Davidson / H counts; total energy. Correctness: arm 2 against arm 1 at the same rank count must be bit-identical (energy bits and density SHA); arms 3 and 4: |dE| against arm 1 within 1e-10 Ha strict, 5e-9 Ha standard, same convergence; every rank's energy equal to rank 0; density SHA of two repeats of arm 4 at np=4 identical.
Optional (one cell, to check the section 5 estimate): si64 strict np=8 arm 1 with `scf.davidson_capacity_per_band = 4`.
Adoption (AGENTS.md adoption rule): each arm that improves the median wall at np=2, 4 and 8 without a correctness failure becomes the default, its variable kept for rollback.
Next (separately requested): a dev08 microbenchmark of MKL ScaLAPACK `pdsyevd` against `dsyevd` at orders 400, 600, 792 on P = 2, 4, 8 (item 6), and the recursive-doubling tree sum (item 3).
