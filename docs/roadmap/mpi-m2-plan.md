# MPI M2: G-vector distributed FFT and distributed SCF stages (plan, not implemented)

Status 2026-11-10. Approved by the user as scope decision U10 (AGENTS.md); builds on ADR 0041 and `docs/roadmap/mpi-m1-plan.md`. Evidence that motivates it: ledger rows MPI-M1-GATE5 and MPI-M1-CEILING.

## Why M1 stops scaling
M1 divides the H application over column pairs. On si32 strict 66 percent of apply calls carry one column, so at most one rank works in those calls; a perfectly balanced pair split cannot exceed 3.9x at 8 ranks. The replicated stages (potential, density, total_energy, mixing) are 7.4 s of 41 s at one rank and grow with the rank count (10.3 s at 8 ranks). Only splitting the work of ONE column, i.e. the 3D FFT, removes both limits.

## Design
- Row partition of the ORC basis by FFT sticks (columns of G vectors along z), balanced by stick length. The M1 algebra (`distributed::{gram, rotate_local, project_out, ...}`, the MGS panel, the reduced matrix) is independent of how rows are partitioned, so it is reused unchanged.
- Distributed paired real FFT: pack the pair, 1D transforms along z on the local sticks, one all-to-all transpose to xy planes owned by plane, 2D transforms on planes, local multiplication by the potential on the owned planes, and the inverse path. The transpose uses `Collective::alltoall_f64` (already in the trait).
- Nonlocal operator: projector columns are row-partitioned like the basis; P^T x is a local GEMM plus the deterministic tree sum of a (projectors x columns) matrix, the coupling is replicated, P W is local.
- Kinetic term is diagonal and local.
- Replicated stages: density accumulation uses the same distributed transforms (planes reduced into the owner), potential build and total_energy work on the plane slabs with scalar reductions; mixing keeps the density in the plane distribution (the Pulay history is distributed, dot products by tree sum).
- Determinism contract (ADR 0041): one rank reproduces the threads-only path bit for bit (the distributed code is only taken when ranks >= 2); a fixed rank count is bit-reproducible; different rank counts agree within the accuracy classes.

## Stages and gates (fixed before coding)
- M2a: stick partition and distributed paired FFT with simulated ranks. Gate: forward then inverse reproduces the input to 4 n eps relative, equals the serial transform within 8 n eps relative for 1..8 simulated ranks, bit-reproducible at fixed rank count.
- M2b: distributed H application (kinetic, local potential, nonlocal). Gate: H applied to random pairs equals the serial `RealHamiltonian::apply_block` within 1e-12 relative for 1..8 simulated ranks, on the existing synthetic pseudopotential fixtures.
- M2c: Davidson on the stick partition (reuse of the M1 solver with the distributed operator). Gate: eigenvalues within the accuracy classes of the serial solver on the existing Davidson tests and the clustered fixtures; total Davidson work (recorded plus certificate) in the perturbed-start distribution of the serial solver.
- M2d: distributed density, potential, energy and mixing. Gate: SCF total energy within 5e-9 Ha (standard) and 1e-10 Ha (strict) of the threads-only run on the engine tests, all converged, density SHA reproducible at fixed rank count.
- M2e: dev08 scaling gate through the master session, perturbed seeds 601-610, MKL build, load-gated, np=1 of the same binary as baseline, with the proof that the distributed branch ran (rank directories, nonzero dE at 1e-13, `parallel.collectives` call counts). Pass requires all of: (1) np=4 speedup >= 2.2x and np=8 >= 3.5x over np=1 (median over seeds) on si32 strict; (2) at least one MPI configuration on the same 8 physical cores is as fast as or faster than the threads-only path (np=1 x T=8); (3) communication seconds (`parallel.collectives`) below 25 percent of the wall at np=8; (4) accuracy and determinism as in M2d. If (2) fails the result is recorded as "MPI does not beat threads on one node" and M2 is not extended further without a multi-node measurement.

## Risks and open points
- dev08 is one node (shared memory and TCP only): the scaling gate measures the algorithm, not network behaviour.
- The paired real FFT pairs two columns per transform; the distributed version keeps pairing (a single-column call still pays a pair), unless the second-vector work (docs/roadmap/second-vector-plan.md) changes the call mix.
- The plan fixes interfaces only; the grid slab partition, the stick balancing rule and the FFT backend per rank (MKL DFTI, rustfft) are decided in M2a with measurements on the transform sizes of si16 and si32.
