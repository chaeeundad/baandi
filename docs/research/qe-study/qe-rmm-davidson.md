# QE 7.6 study note: the `rmm-davidson` band solver (Davidson warm-up, then single-sweep RMM-DIIS), real Gamma path

Status: analysis note written by Claude Code in its own words after reading QE 7.6 sources (policy: AGENTS.md, "QE source study policy"). No code, no excerpts, no copied tables; file and routine names are references only. Purpose: a Codex implementer builds a switchable `rmm-davidson` mode of the reference SCF driver from this note alone.
References read in full: `KS_Solvers/RMM/rrmmdiagg.f90` (real Gamma RMM-DIIS), `KS_Solvers/DENSE/gram_schmidt_gamma.f90`, `KS_Solvers/DENSE/rotate_xpsi_gamma.f90` (product structure), `KS_Solvers/DENSE/rotate_driver.f90`, `PW/src/c_bands.f90` (gamma branch of the band driver, the exit test, the two schedule functions at its end). Read in the relevant parts: `PW/src/electrons.f90` (iteration counter, threshold schedule, first-step redo), `PW/src/input.f90` (solver names), `Modules/read_namelists.f90` (input defaults), `PW/src/setup.f90` (Davidson workspace in this mode), `PW/src/memory_report.f90`. Compared, not re-read line by line: `KS_Solvers/RMM/crmmdiagg.f90` (complex k-point twin; same constants and structure, complex Hermitian small problems).
Not repeated here: the Davidson call itself (`davidson-iteration-structure-gamma.md`), the half-sphere real layout (`qe-gamma-real-path.md`), the short survey of the same solver in `qe-solvers-and-subspace.md` section 3.3 (this note supersedes it where they differ; they do not contradict).

Notation: `n` bands requested, `r` real coordinates per vector on this rank (BAANDI ORC rows; QE's `2*npw` real numbers), `u` unconverged bands at a step (`u <= n`), `m` the DIIS depth (input `diago_rmm_ndim`, default 4), `k` the DIIS step index `1..m`, `b` the Gram-Schmidt block size (input `diago_gs_nblock`, default 16), `P` plane-wave ranks. S is the overlap operator; for BAANDI (norm-conserving) S = 1 and every "S-product" below is the vector itself. Energies: QE works in Ry; every rule below is stated so it is unit-free or the unit is named.

## 1. Purpose

`rmm-davidson` replaces the full iterative diagonalization in most SCF iterations by a cheap, deliberately unconverged refinement: one Rayleigh-Ritz rotation of the current orbitals in the new Hamiltonian, then at most `m` band-by-band residual-minimization steps, then an orthonormalization. It never builds a growing subspace, never solves a dense problem larger than `n` (once per SCF iteration) or `m` (per band), and spends at most `n + m*n` H columns per SCF iteration with no other block-subspace machinery. The SCF loop itself is relied on to converge the orbitals over iterations; a full Davidson solve is inserted periodically to repair what band-by-band refinement cannot (band crossings, lost states, ordering).

## 2. When QE uses Davidson and when RMM-DIIS (the switching rule)

2.1 Names. The input values `rmm`, `rmm-diis` and `rmm-davidson` are synonyms: all select the RMM solver alternated with Davidson. `rmm-paro` alternates with ParO instead (not in scope here).

2.2 Rule. Let `i` be the SCF iteration index as printed by pw.x (1 for the first iteration, incremented once per pass of the outer SCF loop). The band solver for iteration `i` is
- Davidson if `i < 3` or `i` is a multiple of 5,
- RMM otherwise.
So the sequence is D D R R D R R R R D R R R R D ... (Davidson at 1, 2, 5, 10, 15, ...; RMM at 3, 4, 6-9, 11-14, ...). The rule depends only on the iteration index; it does not look at dr2, residuals, eigenvalue changes or the threshold. The first-iteration "redo with lowered threshold" (when dr2 at iteration 1 is below the diagonalization error estimate) reuses the same index and is therefore always Davidson.

2.3 Why (principle, partly INFERRED; the source states only "alternate to stabilize the SCF loop"):
- Iterations 1 and 2 start from atomic orbitals in a potential that is still far from self-consistent; a band-by-band local minimizer started there converges to whatever state is nearest and can skip or duplicate states. A full subspace solve gives RMM a start that already spans the right low spectrum. (BAANDI A23 is a direct measurement of what happens otherwise: section 9.)
- Every fifth iteration a full Davidson solve re-establishes a converged, correctly ordered, complete set at the current density. RMM between these points only polishes vectors already close to the right ones.
- The Davidson calls in this mode are ordinary Davidson calls (same threshold rule, same exit test, workspace capacity 2n as in the default Davidson mode because the default Davidson depth is 2), with the GS-orthonormalized RMM output as their start.

2.4 Davidson-or-RMM is decided per SCF iteration, globally (all k points the same). There is no fallback from RMM to Davidson inside an iteration and no test that triggers an extra Davidson call when RMM does badly; the only safeguard is the fixed schedule.

## 3. One RMM iteration of the SCF loop (the driver level)

For one k point (Gamma: once per SCF iteration), in order:

3.1 Subspace rotation (Rayleigh-Ritz in the span of the current orbitals). Input: the `n` orbitals delivered by the previous SCF iteration (orthonormal). Apply H to all `n` (one block of n columns); build the `n x n` reduced matrices `H_red = X^T H X` and `S_red = X^T S X` (both computed, the overlap is not assumed to be identity); solve the dense generalized symmetric problem for all `n` pairs; rotate both X and HX by the eigenvector matrix (two block products). Output: rotated orbitals, their H products (kept, not recomputed), and the Ritz values as the starting eigenvalues `e_i`. This step is mandatory: H has changed since the previous iteration and the H products of the previous iteration are not kept.

3.2 RMM-DIIS sweep (section 4), using the rotated orbitals, their H products and the Ritz values. It does NOT apply H to the starting orbitals again (the H products of 3.1 are passed in).

3.3 The sweep result is accepted whatever its convergence: in an SCF run with the default `diago_rmm_conv = false` the driver sets the unconverged count to zero after the sweep, so the retry loop of the band driver ends after one sweep and the "too many bands not converged" error can never fire in RMM iterations. With `diago_rmm_conv = true` the driver repeats rotation + sweep up to the usual retry limit while more than 5 bands remain unconverged (not the default; not recommended as the first implementation).

3.4 Orthonormalization (block classical Gram-Schmidt), with eigenvalue recomputation and sorting:
- Bands are processed in their current index order (ascending Ritz value from 3.1), in blocks of `b` consecutive bands.
- Inside a block, each band in turn: subtract its projections on the already finished bands of the same block, then normalize. Projections onto earlier blocks were already removed (next point).
- After a block is finished, the projections of ALL later bands onto this block are removed in one block step: overlap `(b x r)(r x n_later)` and update `(r x b)(b x n_later)`.
- The overlaps are taken against the ORIGINAL (input) vectors, i.e. classical (not modified) Gram-Schmidt, one pass, no reorthogonalization. This is adequate only because the input is already nearly orthonormal (it moved little from the orthonormal rotated set).
- Every linear operation applied to a vector is applied identically to its H product, so H X stays consistent without any new H application.
- Then each eigenvalue is set to the diagonal Rayleigh quotient `e_i = x_i^T H x_i` (vectors now normalized), and the bands are sorted by ascending `e_i` (swapping vectors and H products together).
- Equivalent statement: the result is the Q factor (positive diagonal R) of the QR factorization of the band matrix in that column order, followed by a sort. Any QR method gives the same vectors up to rounding (relevant for section 8).
- There is NO Rayleigh-Ritz after the orthonormalization: the delivered eigenvalues are diagonal Rayleigh quotients of a non-rotated set, and the H products are then discarded (the next iteration's rotation recomputes them for the new H).

3.5 What the SCF loop receives: orthonormal orbitals sorted by Rayleigh quotient and those quotients as eigenvalues. Density, band energy and occupations are built from them exactly as after a Davidson call.

3.6 QE's "avg # of iterations" in an RMM iteration is 1 (rotation) + (sum over sweep steps of u)/n + 0.5 (Gram-Schmidt, a nominal charge, not H work). Useful to compare per-iteration work with a dev08 QE output.

## 4. One RMM-DIIS sweep in mathematics

Everything in 4.2 to 4.6 is done band by band (each band has its own history and small problems), but every H application is batched over all currently unconverged bands. Bands do not interact during the sweep: no orthogonalization between bands, no shared subspace.

4.1 State per band `i`: current vector `x`, `Hx` (and `Sx`), the scalars `h = x^T H x` and `s = x^T S x`, the current eigenvalue estimate `e_i`, a converged flag, and a history of up to `m` triples `(phi_j, H phi_j, S phi_j)` with the scalar `eps_j = h` at the time phi_j was stored. Initialization from 3.1: `x` = rotated vector, `Hx` = rotated H product, `h = e_i`, `s = 1`, all flags false, empty history.

4.2 Step `k` (k = 1..m), for each unconverged band:

(a) Store. Append the current `(x, Hx, Sx)` as history entry `k` with `eps_k = h`. (Because every update below ends normalized, `s = 1` whenever something is stored, so `eps_k` is the Rayleigh quotient of `phi_k`.)

(b) Residual history. For each stored entry `j <= k`: `R_j = H phi_j - eps_j S phi_j`. Each history vector carries its OWN Rayleigh quotient (no common shift).

(c) Small matrices (grow by one row and column per step; only the new row is computed, the matrices are symmetric):
- `A_jl = R_j^T R_l` (residual overlap),
- `B_jl = phi_j^T S phi_l` (state overlap),
for `j, l <= k`.

(d) Extrapolation (DIIS).
- If `k = 1`: no extrapolation; `x` is normalized by `1/sqrt(s)` (a no-op in practice), and the working residual is `R = Hx - h Sx`.
- If `k > 1`: find coefficients `v` (length k) minimizing `v^T A v` subject to `v^T B v = 1`, i.e. the lowest eigenpair of the generalized problem `A v = lambda B v` (A is positive semidefinite, so "lowest" equals "smallest magnitude", which is how QE picks it). It is solved by canonical orthogonalization: eigendecompose `B = U diag(beta) U^T`, keep the directions with `beta > 1e-14` (absolute threshold; B has unit diagonal so this is effectively relative), form `W = U_kept diag(beta_kept^-1/2)`, eigendecompose `W^T A W`, take its eigenvector `y` with smallest-magnitude eigenvalue, `v = W y`, rescale `v` so that `v^T B v = 1`. If at most one direction is kept, use `v = (0, ..., 0, 1)` (the newest entry alone).
  Then `x = sum_j v_j phi_j`, `Hx = sum_j v_j H phi_j`, `Sx = sum_j v_j S phi_j`, and the working residual is the SAME combination of the history residuals, `R = sum_j v_j R_j`. This is not the true residual of the new `x` (that would need the new Rayleigh quotient); it is the standard linearized DIIS residual. No H application. `h` and `s` are NOT updated here (4.2(f) recomputes them).
  Difference from textbook Pulay DIIS: the constraint is the normalization `v^T B v = 1`, not `sum v_j = 1`; the problem is a 1-root generalized eigenproblem, not a bordered linear system (this is the Kresse-Furthmüller 1996 "RMM-DIIS" form, PRB 54, 11169, section on residual minimization).

(e) Preconditioning. Band kinetic energy `T = x^T (-1/2 nabla^2) x` (in Ha, with the vector normalized; in QE's Ry units the same quantity is `sum |G|^2 |c_G|^2`). For every plane-wave coordinate with plane-wave kinetic energy `t_G = |G|^2/2` (same unit as T), with `y = t_G / (1.5 T)`:
  `K_G = -(4/3) / T * (27 + 18y + 12y^2 + 8y^3) / (27 + 18y + 12y^2 + 8y^3 + 16y^4)`.
  This is the published Teter-Payne-Allan polynomial preconditioner (PRB 40, 12255, 1989) in the scaling of Kresse-Furthmüller 1996; it is not a QE invention. Properties: diagonal in G; depends on the band only through its kinetic energy T (NOT through its eigenvalue, unlike the Davidson preconditioner); `K -> -(4/3)/T` for small `t_G`; `K -> -1/t_G` for large `t_G` (an approximate inverse of `-(H - e)` at high G). It is negative by construction, so the descent direction is `+K R`. Unit check: `K R` is dimensionless in any consistent unit system (K scales as 1/energy, R as energy), so the step lengths below are unit-free. Direction: `p = K R` (elementwise).

(f) Trial step and step length (one H application per band, batched).
- Apply H (and S) to the block of all `p` of unconverged bands: `Hp`, `Sp`.
- Compute six scalars per band: `a11 = x^T H x`, `a12 = p^T H x`, `a22 = p^T H p`, `b11 = x^T S x`, `b12 = p^T S x`, `b22 = p^T S p`.
- The Rayleigh quotient along `x + t p` is `E(t) = (a11 + 2t a12 + t^2 a22) / (b11 + 2t b12 + t^2 b22)`.
- QE does NOT minimize E(t) exactly. It fits a parabola `E(t) ~ E0 + g t + c t^2` from three facts: `E0 = E(0) = a11/b11`, the exact slope `g = E'(0) = 2 (a12 b11 - a11 b12) / b11^2`, and the value `E1 = E(t_ref)` at the fixed trial length `t_ref = 0.5`; `c = (E1 - E0 - g t_ref) / t_ref^2`.
- Step: if `|c| > 1e-16`, `t = -g / (2c)`; otherwise `t = t_max` when `g < 0` and `t = t_min` when `g >= 0`. Then clamp `t` to `[t_min, t_max] = [0.05, 1.0]`. (A concave fit, c < 0, with g < 0 gives a negative t and therefore the minimum step 0.05: QE accepts a short step rather than refusing the update. A step is always taken.)
- Update with normalization: `nrm = sqrt(b11 + 2t b12 + t^2 b22)`, `x <- (x + t p)/nrm`, `Hx <- (Hx + t Hp)/nrm`, `Sx <- (Sx + t Sp)/nrm`, and the scalars `h <- (a11 + 2t a12 + t^2 a22)/nrm^2`, `s <- 1` (exactly; computed by the same formula). No further H application: the new `Hx` is the linear combination.
- Why a clamped parabola and not the exact 2x2 minimum (INFERRED): the fit with a fixed trial length and the clamp act as a trust region; with an unconverged, history-extrapolated residual an exact minimization along p can take long steps that pull the band toward a neighbouring state. The exact minimizer is a legitimate variant to test, but it is not what QE does.

(g) Convergence test (section 5) with the new eigenvalue estimate `h/s`.

4.3 The sweep ends after step `m` or as soon as every band is converged. The output vector of a band is its last updated `x` (after the last step's line search). Note what the history holds: entry 1 is the rotated start, entry `j > 1` is the line-search OUTPUT of step `j-1` (not the DIIS-extrapolated vector). With `m = 4`: four steps, four H applications per still-unconverged band, history of at most four entries, histories discarded at the end of the sweep (nothing is carried to the next SCF iteration).

4.4 Converged bands are frozen for the rest of the sweep: no new history entry, no extrapolation, no preconditioning, no H application. The unconverged bands are packed contiguously for the batched H application of each step, so the H block shrinks as bands converge.

## 5. Convergence test per band (sticky within the sweep)

After each step, with `e_new = h/s` and `e_old` the value from the previous step (from the rotation, for step 1): band i becomes converged if `|e_new - e_old| < thr_i`, with `thr_i = ethr` for an occupied band and `max(5 ethr, 1e-5 Ry)` for an empty band (the same occupied/empty flags and the same `ethr` schedule as the Davidson calls; 1e-5 Ry = 5e-6 Ha). Once converged it stays converged for the sweep. Then `e_old <- e_new` for all bands. There is no residual-norm test. INFERRED (not measured): at iteration 3 and later `ethr` is already small (it follows the dr2-driven schedule), so few bands should meet the test and most bands do all `m` steps: plan with `n + m n` H columns per RMM iteration and measure the actual per-step u.

The test only freezes work; it does not decide acceptance (section 3.3). The per-band eigenvalues delivered to the SCF loop are the Gram-Schmidt Rayleigh quotients of 3.4, not the sweep's estimates.

## 6. Gamma-real specifics

- QE stores half the G sphere as complex numbers and treats each vector as a real array of length `2 npw`; an inner product is twice the real dot product minus the product of the G = 0 real parts (counted once), and the imaginary part of the G = 0 coefficient is reset to zero after every update of x, Hx, Sx and of the preconditioned direction "for numerical stability". The band kinetic energy for the preconditioner uses the same weighting (pairs twice, G = 0 once).
- In BAANDI's orthonormal real coordinates (ORC) these weights are already encoded: every inner product above is a plain real dot product over the local rows, the G = 0 coordinate has no imaginary partner, and no reset is needed. Do not add the half-sphere weight again (same rule as `ref_davidson`, docs/design/reference-scf.md "Inner solver"). The kinetic energy is `T = sum_rows t_G x_row^2`, with both coordinates of a pair carrying the same `t_G`; the preconditioner is a per-row diagonal built from `t_G`.
- All small problems (n x n rotation, k x k DIIS, 2 x 2 line search) are real symmetric. Coefficients are real; combinations keep the real structure.
- The complex k-point twin (`crmmdiagg`) has identical constants and steps with Hermitian small problems; not needed for the Gamma reference driver.

## 7. Cost of one RMM SCF iteration

Real FLOPs (one multiply-add = 2), per rank, NC (S = 1, so S-products are free and the S overlaps are Gram products). Totals assume all bands do all `m` steps (the realistic case, section 5).

| Stage | H columns | Dense block products (shapes) | Vector work (BLAS-1/2) | Reductions over plane-wave ranks |
|---|---|---|---|---|
| Rotation | n | `(n x r)(r x n)` twice (H_red, S_red; symmetric, upper triangle suffices); `(r x n)(n x n)` twice (X and HX): about 8 r n^2 | none | one fused n x n pair (QE: two) |
| Dense solve | 0 | n x n generalized symmetric EVD, all pairs (replicated or root+broadcast) | | (broadcast if root) |
| Sweep step k | u | none: every product is per band | per band: rebuild k residuals (or store them, see 8.3), new row of A and B (two length-k matrix-vector products over r), combine x, Hx, R with v (3 k axpys), kinetic energy, preconditioner, six dots, four axpy/scale: about (14 k + 33) r FLOPs per band | A and B new rows (2 k u values), kinetic energies (u), the six line-search scalars (6 u); QE also broadcasts per-band DIIS coefficients and step coefficients from a root (avoidable, 8.4) |
| Sweep total (m = 4) | 4 n | none | about 270 r FLOPs per band, i.e. about 270 r n | about 3 reductions per step, 12 per sweep, payloads O(m u) |
| Gram-Schmidt | 0 | block overlap `(b x r)(r x n_later)` and two updates `(r x b)(b x n_later)` (X and HX), summed over blocks: about 3 r n^2 | in-block: two matrix-vector products per band against up to b-1 vectors, norms | QE: one reduction per block for the block overlap plus one per band in-block and one per band norm (about 2n + n/b tiny reductions) |
| Eigenvalues, sort | 0 | | n dots | one of n values |

Per RMM iteration: about `(1 + m) n = 5n` H columns, about `11 r n^2` FLOPs in block products (rotation 8, Gram-Schmidt 3), about `270 r n` FLOPs of band-local vector work, one dense `n x n` generalized EVD, `n` tiny `k x k` (k <= 4) EVD pairs per step. No growing subspace, no `(2n)^3` solve, no `r n^2`-sized work inside the sweep.

Memory: history of `m` vectors and `m` H products per band (`2 m r n` doubles; QE's memory report counts this term), plus direction and its H product (`2 r n`), plus X and HX (`2 r n`): about `12 r n` doubles at `m = 4`, versus `4 r n` for the V and HV blocks of a capacity-2n Davidson. RMM mode uses MORE wavefunction memory than Davidson.

Comparison with a warm Davidson call (from `davidson-iteration-structure-gamma.md`): Davidson's floor is `2n` H columns per call (init n + one full correction block) and typical warm calls take more; it does `O(r q n)` subspace GEMMs per step and refresh rotations. RMM trades more H columns at small n (up to 5n against Davidson's 2n to 3n late in the SCF) for the removal of almost all non-H work. It pays only where non-H work dominates the Davidson call (large n, where `r n^2` terms and the `(2n)^3` dense solves grow) or where its band-local structure helps parallel scaling (fewer large collectives). That is the regime to look for in measurements; at small n a loss is expected.

## 8. Implementation contract for BAANDI (what must match, what is free)

8.1 Must match (defines "QE's scheme"):
- Schedule of section 2.2 on the reference driver's SCF iteration index (same index QE prints; the first-step redo keeps index 1). Davidson iterations call the existing `ref_davidson` unchanged.
- Per RMM iteration: rotation with H products kept, ONE sweep of at most `m` steps, unconverged count ignored for acceptance, no retry, no fallback, then orthonormalization in band order with H products carried, Rayleigh-quotient eigenvalues, ascending sort.
- The sweep math of section 4: own-Rayleigh-quotient residuals per history entry, normalization-constrained DIIS solved by canonical orthogonalization (threshold 1e-14) with newest-only fallback, kinetic-energy polynomial preconditioner, trial length 0.5, parabola step clamped to [0.05, 1], normalized update, eigenvalue-change freeze test with the existing occupied/empty thresholds.
- Defaults `m = 4`, `b = 16` as switchable parameters, not tuned constants (generality rule: any later change needs the multi-system evidence).

8.2 Proposed for BAANDI because QE has no counterpart:
- Verification maps (ADR 0043, after the estimator stop) always use Davidson, independent of the iteration index: the certificate needs a converged solve. An estimator stop may occur right after an RMM iteration; that is allowed (QE does the same), the verification then decides.
- Switch: an input key plus an environment variable (for A/B), default off; resolved once at driver entry; recorded in the run record together with the per-iteration solver kind.
- Record per SCF iteration: solver kind, H columns per stage (rotation, each sweep step, Davidson), per-step u, the number of DIIS steps that fell back to newest-only, the number of clamped steps (low, high) and the Gram-Schmidt minimum norm before normalization (a lost or duplicated state shows up there).

8.3 Free choices (same mathematics, implementer may choose):
- Store the residual vectors `R_j` instead of rebuilding them each step (one more `m r n` of memory, removes the rebuild term).
- Orthonormalization by any QR in band order (block classical GS as in QE, two-pass GS, or CholeskyQR on the nearly orthonormal block with a pivot guard and fallback). Result equals QE's up to rounding. If a Cholesky route is used, the H products must be transformed by the same triangular inverse.
- The rotation may reuse the existing reduced-solve adapters (guarded Cholesky + symmetric EVD). The overlap must be measured, not assumed identity.
- Batched per-band small products may be written as block-diagonal batched dot products; no cross-band GEMM is implied.

8.4 MPI (G-distributed rows): every small problem should be solved REPLICATED on all ranks from the reduced data (deterministic, no broadcast), which removes QE's per-band coefficient broadcasts. Fuse the per-step reductions where the data dependence allows: the A/B rows of step k in one reduction; the kinetic energies after the extrapolation; the six line-search scalars in one. Gram-Schmidt per-band reductions (QE's) should be avoided: compute the in-block overlaps of a whole block in one reduction (or use CholeskyQR, one reduction per pass). The threads-only path must remain bit-identical when the switch is off.

## 9. Differences from BAANDI's earlier hybrid (A22/A23)

The P10A hybrid (ledger A22, `docs/research/rmm-diis-design.md`, `docs/performance/rmm-diis-benchmark.md`) was a different algorithm in a different regime. It needed 2 to 10 times the H columns of Davidson (k444s: hybrid 136022 against Davidson 17628) and lost on every complete SCF case. The differences that explain this, and what must be replicated to avoid repeating it:

| Aspect | A22 hybrid (old path) | QE `rmm-davidson` |
|---|---|---|
| Regime | Every SCF step's eigensolve converged to a strict residual (1e-9 class), up to 80 sweeps, stagnation test | One sweep of at most 4 steps per SCF iteration, result accepted unconverged; the SCF loop converges the orbitals |
| When Davidson runs | Warm-up to a residual threshold (1e-3) when the start was screened as poor; fallback on stagnation | Fixed schedule: SCF iterations 1, 2 and every fifth; never as a fallback |
| Extra H work inside a solve | Residual screening of the start, re-orthonormalization + Rayleigh-Ritz + fresh H X for certification after each inner block, fallback Davidson work | Exactly n (rotation, H products kept) + at most m n (one trial direction per step); Gram-Schmidt carries H X, no certification |
| DIIS form | Affine (coefficients sum to 1), regularized residual Gram inverse, fixed Ritz shifts within a sweep | Normalization-constrained generalized eigenproblem, each history entry with its own Rayleigh quotient, canonical orthogonalization |
| Step length | Least-squares residual step at fixed shift | Clamped parabola fit of the Rayleigh quotient (trial 0.5, range [0.05, 1]) |
| Preconditioner | Davidson's diagonal with the Ritz value | Kinetic-energy polynomial preconditioner using the band kinetic energy |
| Convergence test | True residual norm <= tolerance | Eigenvalue change < ethr, only freezes work |
| Inter-band orthogonalization | After every inner block | Once per SCF iteration, after the sweep |
| Path | Complex, old legacy driver with its tolerance schedule | Real Gamma, QE threshold schedule (the reference driver already has it) |

A23 (warm-up skipped for an orthonormal coordinate start: Si Gamma 2750 vs 470, four-cell Si 46480 vs 2679) is the failure that QE's schedule rules out structurally: RMM never sees the atomic start, because iterations 1 and 2 are always Davidson.

What must be replicated to avoid the A22/A23 failure: (1) the fixed iteration schedule with Davidson at 1 and 2; (2) one bounded sweep per SCF iteration with no convergence requirement and no internal certification or fallback; (3) reuse of the rotation's H products and carrying H X through Gram-Schmidt, so that H columns per RMM iteration are exactly `n + sum_k u_k`; (4) the certificate only at the end (verification map with Davidson). A22's negative result says nothing about this regime (the ledger's "premise changed" remark on A22 already records this).

## 10. What to measure to decide adoption

Under the adoption rule (any measured dev08 wall gain, reject only for hard-coding or correctness failure):

1. Correctness (must hold on every case): converged and verified; energy difference to the Davidson-only reference path within the standard class 5e-9 Ha (strict 1e-10 Ha where the case is run strict); verification pass on the first attempt as often as with Davidson; thread-count independence of results; no lost or duplicated state (Gram-Schmidt minimum pre-normalization norm, occupied eigenvalues against the Davidson run's).
2. Work: SCF iteration count (the main risk: unconverged maps can add iterations and verification attempts), total H columns, H columns per RMM iteration against Davidson iterations at the same index, and per-stage times (rotation, sweep, Gram-Schmidt, dense solve) next to Davidson's restart/border/correction clocks.
3. Wall: median over alternating A/B runs on dev08, one core and MPI lanes separately (comparison rule), at the same cells as the QE comparison.
4. Generality set (AGENTS.md generality rule): cubic, skewed, elongated, slab; small (si8) to large (si128 and above) band counts, where the crossover is expected at large n (section 7); small and large gap; cases with empty bands; at least two projector structures.
5. QE side, same conditions: run pw.x with `diagonalization='rmm-davidson'` (defaults) on the benchmark cells. If it is faster than QE Davidson for a cell, it becomes QE's best configuration for that cell and the P0.3 target moves (work order 2026-11-14 item 2). Its per-iteration "avg # of iterations" lines (section 3.6) give the H work per band per iteration to compare with BAANDI's record.
6. Diagnostic sweeps (information, not tuning): m in {2, 3, 4, 6}, exact 2x2 minimization instead of the parabola, and the Davidson period 5 versus 3 and 10, each recorded on the full generality set; a default changes only with evidence on all classes.

## 11. Uncertain or not determined from the source

- The reason for the schedule (iterations 1, 2, every fifth) and for the constants (trial 0.5, clamp [0.05, 1], depth 4, block 16) is not documented in the source beyond "stabilize the SCF loop"; the explanations above are INFERRED from the algorithm.
- In the first DIIS step with S = 1 the source forms the working residual from the vector's stored S product, which for norm-conserving data is the vector itself; the formula `R = Hx - h x` is therefore exact for BAANDI. (For ultrasoft data the choice of product there would matter; irrelevant here.)
- The convergence test compares the band's estimate across steps by band INDEX; bands are not reordered during the sweep, so this is consistent, but two bands drifting toward the same state are not detected until Gram-Schmidt (where the later one loses norm) or the next Davidson iteration. QE has no explicit guard; the only failure signal is a Gram-Schmidt "linearly dependent" error (norm below 1e-16).
- Whether QE's RMM mode is faster than its Davidson on our benchmark cells is unknown (no dev08 run yet); section 10 item 5.
- `diago_rmm_conv = true` (converged RMM with retries) was read but is not described in detail because it is not the default and resembles the A22 regime.

## 12. Hypotheses supported by this note

- R1. On large-band cells (n of order 100 and more) the RMM iterations remove most of Davidson's non-H work (restart rotations, borders, `(2n)^3` dense solves) at the price of about 5n H columns per iteration; net wall gain if the SCF iteration count rises by at most about one. Measure per section 10.
- R2. On small-band cells RMM loses (more H columns, little non-H work to remove); expected outcome, keep switch off there unless measured otherwise.
- R3. In the MPI lane the band-local sweep needs only small reductions (O(m n) values per step) and one `n x n` reduction per iteration, so it may scale better than Davidson's border reductions; measure collective time per iteration on dev08.
