# QE 7.6 study note 1: the iterative eigensolver and the SCF threshold rule

Status: analysis note written by Claude Code in its own words after reading the QE 7.6 sources named below (policy: AGENTS.md, "QE source study policy"). No code or excerpts. File and routine names are references only. Everything in the section "Tuning hypotheses" is a hypothesis to be tested on rust-dft by measurement; "QE does it" is never the justification.
References read: `PW/src/electrons.f90` (SCF loop and threshold update), `PW/src/c_bands.f90` (band driver), `KS_Solvers/Davidson/cegterg.f90` (block Davidson for k points), `PW/src/g_psi.f90` and `PW/src/usnldiag.f90` (preconditioner), `Modules/read_namelists.f90` (defaults), `PW/src/sum_band.f90` (band types).

## 1. How the block Davidson is organized
Let `m` be the number of requested bands and `nvecx = d*m` the maximum subspace dimension, with the default multiplier `d = 2` (input `diago_david_ndim`), so a restart is needed whenever more than `m` new vectors would not fit next to the current `m`.
1. The subspace basis is never orthonormalized. It starts from the `m` current vectors. New basis vectors are only normalized to unit length. The reduced problem is therefore the generalized Hermitian problem `H_red c = e S_red c` with `S_red` the overlap of the basis; it is built and solved by a dense generalized solver on matrices of dimension at most `nvecx`.
2. The products of H with all basis vectors are stored. Each expansion applies H only to the newly added vectors, and only the new rows and columns of `H_red` and `S_red` are computed (one block product of the new vectors with the whole basis).
3. The correction for an unconverged root `i` is the residual `r_i = (H - e_i S) psi_i`, where `psi_i` is the current Ritz vector. It is assembled from the stored basis and the stored H products by one block product with the Ritz coefficients (no new application of H), then multiplied by an approximate inverse (preconditioner), then normalized and appended. Only the unconverged roots get a correction, so the number of H applications in an iteration equals the number of unconverged roots.
4. Convergence of a root is judged by the CHANGE of its Ritz value between two consecutive iterations, `|e_new - e_old| < threshold`, not by the residual norm. The threshold for bands whose occupation is negligible (below 1 percent of the maximum; they are tagged by the band-sum step) is looser: `max(5*ethr, 1e-5)`.
5. A restart ("refresh") happens when all roots converge, or when `nbase + (number unconverged) > nvecx`, or at the iteration limit. It contracts the stored basis AND the stored H products with the Ritz coefficients to form the new first `m` basis vectors and their H products. H is not applied again. The reduced matrices are reset to the diagonal of Ritz values and the identity.
6. Energies and thresholds are in Ry. The count of H applications is tracked per call and reported as the average number of iterations per k point.

## 2. The preconditioner
The diagonal approximation of `H` for plane wave `G` is the sum of the kinetic energy of `|k+G|^2`, the G=0 component of the local potential (the mean value of the local part), and the diagonal of the nonlocal pseudopotential in the plane-wave basis; the overlap diagonal is 1 for norm-conserving pseudopotentials. The correction is divided by a smooth, strictly positive regularization of `x = h_diag - e`: `denominator = 0.5*(1 + x + sqrt(1 + (x-1)^2))`, which behaves like `x` for large positive `x` and never crosses zero for negative `x` (the older variant clips `|x|` to 1e-4).

## 3. The SCF threshold rule
`dr2` is the estimated SCF error in energy units (a Hartree-metric norm of the density residual, not a plain L2 norm), `nelec` the electron count.
- Iteration 1 uses the initial value (the namelist default 0 means: set internally); from iteration 2 on the cap is `1e-2` (`1e-5` for a special mode).
- From iteration 2 on: `ethr = min(ethr, 0.1 * dr2 / max(1, nelec))` (0.05 for the orbital-resolved Hubbard case), floored at `1e-13`.
- Safeguard after the first diagonalization: if the measured `dr2` turns out smaller than the estimated diagonalization error `ethr * max(1, nelec)`, the threshold was too loose: set `ethr = max(0.1*dr2/max(1,nelec), 1e-13)` and diagonalize that iteration again.
- For insulators QE by default computes only the occupied bands (`nbnd = nelec/2`); extra empty bands exist only on request.
- Davidson is initialized with a Rayleigh-Ritz rotation in the span of atomic plus random starting functions at the first SCF iteration; later iterations start from the previous eigenvectors without rotation.

## 4. Differences from rust-dft today (verified in our code and logs)
| Item | rust-dft | QE 7.6 |
|---|---|---|
| Basis | orthonormalized against the whole basis (two Gram-Schmidt passes per correction) | not orthonormalized; generalized reduced problem |
| Subspace capacity | 48 for 8 bands (6 times) | 2 times the band count |
| Restart | orthonormalize the Ritz vectors, re-apply H to all of them | contract stored vectors and stored H products, no new application |
| Convergence test per band | residual norm below the tolerance (strict; H applications per iteration about 8.1 of 8 bands) | Ritz-value change below the threshold; empty bands looser |
| Empty bands | same tolerance as occupied bands; default 4 extra bands | looser threshold; none by default for insulators |
| Preconditioner | kinetic diagonal | kinetic plus mean local potential plus nonlocal diagonal, smooth positive denominator |
| SCF threshold | `min(0.1 * L2 density residual, 1e-2)`, floor `0.1*eigen_tolerance` | `min(ethr, 0.1*dr2/nelec)`, floor 1e-13, with a first-iteration redo safeguard; `dr2` is an energy-weighted estimate |
| Davidson iterations per k point over the whole SCF (8 k points) | about 53 (after O4c) | about 34 |

## 5. Tuning hypotheses (to test one by one on rust-dft; keep only what helps and passes the gates)
H1. Judge in-loop Davidson convergence by the Ritz-value change (an eigenvalue-error measure, quadratically smaller than the residual): fewer iterations, no residual evaluation for the test. The strict residual criterion remains for the final fixed-point certificate, so the guarantees do not change.
H2. Use a looser in-loop threshold for empty bands (their eigenvalues do not enter the density); keep the lowest empty band accurate enough for the sampled-gap criterion and solve all bands tightly in the certificate pass.
H3. Contract stored basis and H products at a restart instead of orthonormalizing and re-applying H (saves about `m` applications per restart), and try the subspace multiplier 2 and 3 instead of 6 (smaller reduced problems, more frequent restarts).
H4. Replace explicit orthonormalization of corrections by normalization plus a generalized reduced eigenproblem. Expected to remove most of the orthonormalization and `hypot` cost (about 20 percent of the profile after O4c). Risk: loss of accuracy when the basis becomes nearly dependent; needs a conditioning guard and a fallback.
H5. Improve the preconditioner with the mean local potential and the nonlocal diagonal and the smooth positive denominator: fewer iterations.
H6. Drive the SCF threshold by an energy-weighted estimate of the SCF error (a Hartree-metric density residual divided by the electron count) with a first-iteration safeguard, instead of the plain L2 residual.
H7. Measure the effect of the number of extra bands (8 versus 5 versus 4 bands) on time; the default is not changed by this note (validation inputs use 8 bands to match the QE comparison), but the cost is reported.
H8 (FFT, from the stick layout): QE transforms whole sets of 1D columns contiguously (many 1D FFTs per call along the contiguous axis, with the other directions done on columns that exist) which corresponds to the pruning already done in O3; the remaining difference is the memory layout: contiguous batches of lanes instead of one strided lane at a time. Test lane batching and a transposed layout for the strided axes (bit-identical).
Order by risk: H2, H3 (contraction, then multiplier), H5, H1, H6, H8, H4 last.

## 6. Measured outcomes on rust-dft (2026-11-03, ticket P7P7 and P7P9; numbers in docs/performance/optimization-log.md)
Hypotheses are not adopted because the reference code does it; this is what measurement on our code showed.
- H3a (restart by contraction of the stored products): kept. Zero change for the default benchmark cases (they converge before the first restart: capacity 48 for 8 bands), bit-identical density, and on test problems that restart the operator applications fall from 165 to 114 and from 92 to 64. Expected to matter for larger systems with more bands.
- H3b (smaller subspace capacity, 2, 3, 4 times the band count instead of 6): rejected. Capacity 2 times: more applications (452 against 373 on Gamma) and slower; 3 and 4 times: no gain or a loss. With corrections appended after orthonormalization in our solver, the larger subspace converges in fewer total applications.
- H5 (fuller diagonal preconditioner): rejected. The Gamma case needed 11 SCF steps instead of 10 and the iteration counts did not drop.
- H1 (Ritz-value-change convergence with strict refinement): rejected. More solver iterations and slower on the benchmark cases.
- H4 (generalized non-orthonormal reduced problem): not attempted, because the stages that it builds on were negative and it carries the highest numerical risk.
- H8 (FFT lane batching with four adjacent lanes, bit-identical): measured, small gain: 1 to 4 percent for the transform alone on production-size grids, no measurable gain end to end (1.00 to 1.02 times in an alternating A/B comparison), losses on tiny grids. Not merged.
- H2 (looser threshold for empty bands), H6 (energy-weighted SCF estimate with the first-iteration safeguard), H7 (number of extra bands): not yet tested; they act at the SCF level and are the next candidates.
Reading: where our solver differs from the reference by design (orthonormalized basis, capacity 6 times the band count) the reference's parameter choices do not transfer; the transferable idea so far is the restart contraction.
