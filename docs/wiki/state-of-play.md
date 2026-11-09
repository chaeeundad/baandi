# State of play (2026-11-09)

## Speed against QE 7.6 (dev08, Xeon Gold 6230; full tables in docs/performance/qe-vs-baandi-current.md)
- Like for like, 1 core, QE at its BEST start (`startingwfc = atomic`), BAANDI MKL at current defaults, unperturbed: BAANDI over QE wall, si8 standard 1.4x / strict 2.7x, si16 2.2x / 3.0x, si32 3.9x / 5.7x (median 2.9x, worst 5.7x; it was 4.0x and 11.6x before the atomic start).
- 8 cores, mixed models (QE 8 MPI ranks against BAANDI 8 threads, same 8 physical cores, measured): si32 strict 9.9x, standard 7.3x, si16 strict 5.0x, standard 3.5x. The like-for-like multi-core comparison needs BAANDI's MPI layer.
- Start sensitivity: QE's SCF count does not move under 1e-3 Bohr displacements; BAANDI with the free-electron start moved 40 to 59 percent on 1e-7 Bohr; the atomic start removes the spread (SCF spread 7 to 0, 2 to 0).
- Hard cells (si-skew, si-elongated, si-slab, fixed occupations): QE needs 51 to 200 iterations and fails on two of eight runs; BAANDI numbers on those cells were requested (cell C of the dev08 package).

## Defaults that changed in this campaign (each has an ADR or ledger row)
gamma_real auto (ADR 0040), start_wavefunctions atomic (ADR 0042), MKL builds: dgemm for the real Davidson products (MKL-R1/R2), blocked Gram-Schmidt (MGS-B), dsyevd for the reduced eigenproblem (MKL-EVD), band-parallel audit (AUD1), no redundant finiteness scans (FIN1), results.json always records gamma_real and start_wavefunctions.

## Where the gap lives now (stage decomposition, 1 core)
Davidson algebra was 64 to 86 percent of the gap before the atomic start; FFT is at or below QE per call (BAANDI's mirror pair transform is faster than MKL DFTI and QE's fftw); the nonlocal operator is 6 to 12x QE's calbec plus add_vuspsi (candidate: GEMM path); serial pieces limit threads (reduced eigensolve, energy audit before AUD1); QE's thread scaling (MPI, 4.26x on 8 ranks at si32 strict after its atomic start) is better than BAANDI's (2.43 to 2.55x on 8 threads).

## Open work
1. MPI M1 (ADR 0041, docs/roadmap/mpi-m1-plan.md): row-block basis, local algebra, column-pair H, deterministic tree reduction; pre-registered scaling gate (si32 strict 2.2x at 4 ranks, 3.5x at 8 ranks).
2. Nonlocal operator through GEMM (MKL feature path) and batching of single-column projection products.
3. Remaining gap in per-step cost after the atomic start (Davidson count ratio to QE 3.6x on si32 strict).
4. QE-equal-target table with `local_mesh = "qe_like_g0"` (closes the Si.pz-vbc energy offset, ADR 0037 addendum).
5. Parked by the user (2026-11-09): the idea of starting from the density of a similar system (no design exists; the Hanyang session never discussed it).

## Update 2026-11-10: hard cells
si-elongated and si-slab (the non-converging cells) are fixed-occupation problems whose integer cut lies inside a near or exactly degenerate cluster at the Fermi level; QE limit-cycles on the first (137 iterations) and fails on the second (200). BAANDI with Fermi-Dirac smearing (strict profile) converges both in 20 to 35 iterations. Do not tune mixing for them (ledger ELONGATED-GAP, ELONGATED-SMEARING, SLAB-SMEARING). MPI M1 is correct but slower than threads on one node (ledger MPI-M1-GATE5, MPI-M1-CEILING); decision on M2 pending with the user.
