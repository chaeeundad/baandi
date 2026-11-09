# MPI M1 plan: row-block Davidson with column-pair H application (2026-11-09)

Contract: [ADR 0041](../decisions/0041-mpi-parallel-layer.md) (staged, deterministic, behind feature `mpi`). Target code: real Gamma Davidson of `dft-solvers` and the SCF phases around it in `dft-engine`. The threads-only default path stays bit-identical; no MPI is touched without the feature.

## Layout
- Plane-wave rows (real half-sphere coordinates, n = 5915 si16, 11877 si32) are split into fixed contiguous row blocks per rank (`block_counts`). The Davidson basis V, the products HV, the rotated Ritz blocks and the residuals are stored as row blocks; the small matrices (reduced matrix, coefficients, overlaps, m x b) are replicated on every rank.
- Algebra (rotate, Gram, projection, reduced-matrix update): local GEMM on the local rows; where a row-summed matrix is needed (V^T HV, Q^T X) the local partial is reduced with the deterministic tree sum, so every rank gets identical bits.
- H application: the correction block (n x w, w about 19 in strict, 41 in standard) is transposed from row blocks to column ownership by one all-to-all (about 0.03 ms per call, measured), each rank applies H to its own column pairs with the existing pair FFT and nonlocal code (full vectors, local FFT, no distributed 3D FFT), and the products are transposed back by a second all-to-all.
- Replicated serial parts: the reduced eigenproblem (dsyevd on MKL builds), mixing, the potential; the density and the audit nonlocal term are computed from the bands a rank owns and combined with the tree sum.

## Sub-steps and gates (each ends with check.sh and a ledger row)
- M1a (done): `Collective` trait, `World` (MPI or single rank), in-process `simulate` (threads as ranks) so the algorithms are unit-tested on the Mac without MPI.
- M1b: distributed algebra primitives (row-block Gram, rotate, projection, reduced update) generic over `Collective`; gate: simulated 1..8 ranks equal the serial result within 4 k eps relative, bit-identical to serial at one rank, bit-reproducible at fixed rank count, independent of the thread budget inside a rank.
- M1c: `DistributedRealDavidson` on the primitives plus the two transposes; gate: one rank bit-identical to `RealDavidson` (the same GEMM calls in the same order); simulated ranks within the accuracy classes on the existing Davidson tests.
- M1d: engine integration (density, audit, certificate, atomic start Ritz step) with rank-0 output; gate: single rank bit-identical to the threads-only SCF (`bench_scf --require-identical` with `[parallel]` absent); simulated multi-rank SCF within the classes.
- M1e: CLI and results (`[parallel] mpi = "auto"`, launcher under mpirun, rank-0-only output, `parallel` block in results.json with rank count and reduction rule).
- Scaling gate (pre-registered now, dev08 only, MKL native, one socket, ranks x 1 thread, si32 strict and standard at the atomic default, unperturbed and perturbed seeds as in MGS-B R2): speed-up of P ranks over 1 rank at least 2.2x for P = 4 and 3.5x for P = 8 on si32 strict (threads reach 2.43x at T = 8); all converged, |dE| in class, SHA identical run to run at fixed P; no case at P = 4 slower than threads at T = 4 by more than 3 percent. If the gate fails report which part and stop; no tuning.

## Not in M1
Distributed 3D FFT (M3), spin, finite temperature, forces and stress, k-point parallelism (option D, a separate small step), multi-node (not measurable on dev08).

## Status 2026-11-09
M1a, M1b done. M1c (`distributed_davidson.rs`), M1d (`dft_engine::parallel`, branch in `gamma_real.rs`) and M1e (`dft-cli --features mpi`) implemented and tested on simulated ranks; the dev08 smoke run and the pre-registered scaling gate are the open items (ledger row MPI-M1CDE).
