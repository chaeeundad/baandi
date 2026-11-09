# Davidson interface review and MPI sequencing (design note, 2026-11-09)

Status: proposal. No MPI code is written; MPI stays out of scope (AGENTS.md) until the user approves a ticket. Source: the user's four Davidson design points relayed by the dev08 master session, checked against the code at f75a817 and the dev08 measurements.

## The four points against the code
1. **Do not materialise H (apply it as FFT + kinetic + nonlocal): already the design.** `dft-solvers` `LinearOperator` and `RealLinearOperator` (`apply_block`); no dense Hamiltonian is built anywhere; H psi runs through `real_gamma::apply_block` (pair FFT, kinetic, `NonlocalOperator::add`). Nothing to do.
2. **Block Davidson for BLAS-3: the interface is blocked, the blocks are as wide as the number of unconverged bands.** Measured si32 strict: 5787 H columns over 303 Davidson iterations, 19.1 columns per iteration (standard 40.7); the dgemm census has 25 percent of the calls with n = 1 at 9.1 GFLOP/s against 58 GFLOP/s for wide calls. Careful: the width is the unconverged-band mask (the same mask QE uses with `notcnv`), so widening it to all bands adds H applications for converged bands; it does not simply remove the n = 1 calls. Candidates that do not add work: pad an odd column count to even for the pair FFT (the second band of a pair is free); batch the single-column projection products into the next correction block. To be measured, not assumed.
3. **Orthogonalisation and reductions may outgrow the FFT under MPI: there is no MPI yet.** Measured on threads: T=8 si32 strict 2.77 of 8.5 thread-stacks working, wall speedup 2.55x; orthogonalisation share at T=1: projection 5.3 percent plus dgemm panels 3.6 percent; blocked Gram-Schmidt (MGS-B) bought -9 to -13 percent per step at T=8.
4. **Preconditioner as an independent interface reusable by LOBPCG: already there.** `Preconditioner` with `apply` and `apply_with_vectors(eigenvalues, residuals, vectors, selected, output)`.

## MPI sequencing (estimates, assumptions stated)
Single node, shared-memory all-to-all, about 20 GB/s effective, no communication/compute overlap, current 19-column blocks, si32 strict (grid 72x72x36, complex 16 B, one 3D FFT per H column):
| decomposition | dominant traffic | estimate |
|---|---|---|
| G-vector distributed | FFT all-to-all | 5787 columns x 2.99 MB, about 17 GB, about 0.9 s at 8 ranks on a socket |
| band distributed | orthogonalisation allreduce | about 1300 reductions x 50 to 185 KB, about 0.2 GB, latency dominated |
So the user's concern (orthogonalisation and reduction outgrowing the FFT) holds for band distribution and is reversed for G distribution. The choice of layout depends on how wide the Davidson blocks become, which argues for settling the block width before an MPI layout. MPI is eventually needed (QE's 4 and 8 core figures are MPI; multi-node needs it) but the 1-core gap (si32 strict 11.6x against QE's best start) and the thread scaling (2.55x) come first.

## Suggested order
1. Block-width work that does not add H applications (see 2).
2. Early-iteration eigen tolerance: note that `initial_eigen_tolerance` and `eigen_tolerance_decay` only apply to the legacy schedule (`residual_eigen_schedule = false`); the default schedule is tau = min(1e-2, 0.1 x density residual), floor 0.1 x eigen_tolerance, so the knobs to screen are the 0.1 factor, the 1e-2 cap and the stopping predicate of the first three SCF steps (strict uses a residual-norm criterion, QE a Ritz-value change). Screen in a clone first.
3. Atomic start wavefunctions (R3 stage 2 in flight).
4. MPI design ticket after 1 to 3, with the layout chosen on the block width and the measured thread scaling.
