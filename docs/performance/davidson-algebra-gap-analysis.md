# Davidson algebra gap to QE at si64 strict (analysis, 2026-11-14)

Question: at si64 strict, one core, MKL, BAANDI's Davidson algebra (rotate, orthonormalize, eigen, reduced) is 5.8x QE's (regterg:last, overlap, update, rdiaghg, about 8.5 s, [qe-stage-decomposition-2026-11-14.md](qe-stage-decomposition-2026-11-14.md)), while h_apply is 1.6x. This note says what each stage does per Davidson iteration, what it costs, what is avoidable, and which change is worth trying first. It is an analysis, not a measurement of a change. QE facts come only from `docs/research/qe-study/` (no QE source read).

## Evidence and method
- Code read: `crates/dft-solvers/src/real_davidson.rs` (solve_research, reduced_matrix_into, ritz_and_residuals, Projection::run_into), `crates/dft-numerics/src/real.rs` (dispatch of the GEMM sites).
- Counts and shapes: a temporary per-iteration trace (removed; source restored) on the Mac, pure build of 40b6b84, T=8, current defaults (B1 on with count sketch, B5 on, B6 on, DGKS 0.1). Inputs from `scripts/bench_suite.py` helpers, seed 601, mixing beta removed: (a) si64 f1e-3 (fractional 0.001, the QE-table condition) and (b) si64 at fractional 4.87e-9 (about 1e-7 Bohr, the P48A condition), plus si32 f1e-3. Counts are exact for these trajectories; dev08 MKL trajectories differ a little (dev08 f1e-3: 176 fine Davidson iterations and 10068 H columns, here 168 and 10140).
- Seconds are a dev08 model, not a measurement: four effective rates fitted to the relayed P48A si64 strict T=1 MKL stage times (b1e-7: rotate 20.6, orthonormalize 18.3, reduced 4.0, eigen 5.6 s) using trajectory (b): rotate (NN dgemm) 48.8 GFLOP/s, reduced (TN dgemm with tiles) 38.9 GFLOP/s, orthonormalize 31 GFLOP/s for the dgemm part plus 8.4 GFLOP/s for the part that runs on the nano kernel (finding 1; 8.4 GFLOP/s is the pure-build rotate rate of P48A on the same flops), eigen 0.281 ns per m^3 of each full EVD. Cross-check: the model gives algebra 55.8 s plus h_apply about 38.6 s for the f1e-3 trajectory, about 95 s against the measured 101.1 s (diagonalization 80.42 + ladder coarse 20.63) of the QE stage table.
- Flops are 2mnk. n = 23871 real rows (fine, 72^3 grid), 8409 (ladder coarse); 132 bands (128 strict, 4 buffer); capacity 6 x 132 = 792.

## 1. What one Davidson iteration does (si64 strict, default path)
| step | operation | shape (rows x inner x cols) | needed for |
|---|---|---|---|
| reduced | V^T (HV new columns), new border only (prefix kept; full m x m after a restart or solve start) | m x n x a | reduced matrix |
| eigen | dense EVD of the m x m reduced matrix, all pairs (MKL dsyevd); B5 secular update when exactly one column was appended | m^3 | Ritz values, coefficients |
| sketch | B1 count-sketch estimates of all 132 residuals (s = 128 rows) | small | decides which roots get explicit residuals |
| rotate | X = V C and H X = HV C for the explicit set: all 132 bands on the full path, the B1 set on narrow iterations; then R = HX - X diag(theta) and norms | 2 x (n x m x k) | residuals, corrections, restart, return |
| precondition | kinetic diagonal on the uncorrected residuals of unconverged roots | n x c | corrections |
| orthonormalize | normalize; Q^T X and X -= Q S against the whole basis; second pass only if DGKS 0.1 asks; in-block MGS (blocked 16 on MKL) | 2 x (m x n x c) per pass, plus 4 n c^2 | orthonormal basis |
| h_apply | H on the c kept columns | n x c | basis products |

Run totals, si64 strict f1e-3 (fine / coarse ladder stage):
| item | fine | coarse |
|---|---|---|
| Davidson iterations (all solves incl. certificate solves) | 168 | 68 |
| full EVDs / B5 bordered | 150 / 18 | 60 / 8 |
| narrow (B1) / full-path / returning iterations | 83 / 73 / 12 | 24 / 38 / 6 |
| mean basis m (max) | 375 (791) | 559 (763) |
| mean rotated columns: narrow / full | 47.6 / 132 | 55.2 / 132 |
| mean corrections c: narrow / full | 5.0 / 104.3 | 3.8 / 86.5 |
| appends / appended columns / second projection passes (columns) | 156 / 8028 / 63 (653) | 62 / 3379 / 4 (16) |
| contraction restarts / rebuilds | 5 / 0 | 0 / 0 |
| H columns: initial + appended + certificate | 1584 + 8028 + 528 = 10140 | 4171 |

So 8028 of the 10140 fine H columns (79 percent) were projected against the whole current basis (mean 392 columns on full-path appends, 298 on narrow ones), 653 of them twice; initial and certificate columns are not projected.

## 2. Flops and seconds per stage (si64 strict, f1e-3, dev08 model)
| stage | per call (typical) | calls fine/coarse | GFLOP fine + coarse | s fine + coarse | QE counterpart (s) |
|---|---|---|---|---|---|
| rotate | full: 2 x (23871 x 433 x 132) = 5.46 GFLOP; narrow: 2 x (23871 x 298 x 48) = 1.35; return: 2 x (23871 x 546 x 132) = 6.9 | 168 / 68 | 593 + 127 = 720 | 12.2 + 2.6 = 14.8 | update 2.20 + last 1.72 |
| reduced | 2 n m a: full-path 23871 x 433 x 104 = 2.2 GFLOP; narrow 0.07 | 168 / 68 | 195 + 29 = 224 | 5.0 + 0.7 = 5.8 | (H border, inside overlap) |
| orthonormalize | full-path append, one pass: 4 x 23871 x 392 x 104 = 3.9 GFLOP + MGS 1.0; narrow: 0.14 GFLOP but streams Q twice (2 x 57 MB) | 156 / 62 appends | 438 (of which 154 on nano, see 3.1) | 24.2 + 3.3 = 27.5 | overlap 2.51 (H and S border) |
| eigen | dsyevd order m: sum of m^3 14.8e9 + 13.0e9 | 150 / 60 | (BLAS2-bound) | 4.2 + 3.6 = 7.8 | rdiaghg 2.04 |
| total algebra | | | 1382 | 45.5 + 10.3 = 55.8 | 8.5 |

Mapping note: by the step structure in the QE notes (qe-stage-by-stage.md section 4, davidson-and-scf-thresholds.md section 1), `update` forms the residuals of the unconverged roots from stored V and HV, `overlap` builds the new H and S rectangles, `last` is the restart/final contraction. The stage table's pairing (rotate with last, orthonormalize with overlap, reduced with update) is therefore not like for like; only the totals compare. QE has no ladder coarse stage: 10.3 s (18 percent) of BAANDI's algebra is there, 3.6 s of it in eigen (coarse m is larger than fine m).

Rotate split (GFLOP, fine + coarse): residuals of roots that get a correction 308 (6.3 s); Ritz vectors and H X for 5 contractions and 18 returns 148 (3.0 s); explicit residual checks of roots that turn out converged plus all-band rotation on full-path iterations 265 (5.4 s). For trajectory (b), the P48A regime, the last item is 688 of 1006 GFLOP (14 s): 43 full-path iterations append one column only, because the B1 explicit set held on average 98 of 132 roots on full-path iterations (at least half the bands, so B1 kept the all-band path).

## 3. Findings
### 3.1 The projection apply runs on the slow pure kernel on MKL builds
`Projection::run_into` calls `multiply_slices_into` (site ProjectionApply) for X -= Q S. Its dgemm branch requires rows x cols <= capacity^2 for a non-transposed product (`real.rs`, multiply_slices_impl); otherwise it falls back to `multiply_into` -> `rotate_legacy` with the MKL-build default `BAANDI_PURE_GEMM=nano`. The bound is c <= 792^2/23871 = 26 columns at si64, 14 at si32 (408^2/11877), 7 at si16. At si64 f1e-3, 7605 of 8028 appended columns (95 percent) and 154 GFLOP take the nano path. The same bound sends the MGS-B cross-panel apply (16 columns, site Other) to nano at si16 and si32 (not at si64). Existing dev08 evidence: DGEMM-SITES (si32 strict) shows projection_apply 8.72 GFLOP through dgemm against 53.7 GFLOP for projection_overlap, although both products have identical flops per pass: about 45 GFLOP of apply ran outside dgemm. The model attributes 18.3 s of the 27.5 s orthonormalize stage at f1e-3 (12.2 of 18.3 s at b1e-7) to this fallback. This is a dispatch gap left by MKL-R2 (which meant to move the projection products to dgemm), not an algorithm property.

### 3.2 Most rotation is spent on roots that receive no correction
Necessary per iteration: residual vectors of the corrected roots (4 n m c, same as QE's update), plus Ritz vectors and H X at restarts and returns. BAANDI also forms residuals for every root whose sketch estimate is at least half its threshold (strict mode), and on full-path iterations for all 132 roots. Distribution of those explicit roots that turned out converged (fine, f1e-3, narrow iterations): residual/threshold below 0.6: 1639, 0.6-0.7: 1046, 0.7-0.8: 487, 0.8-0.9: 278, 0.9-1.0: 82; corrected: 415. Sketch accuracy over all fine iterations of si32 and si64 (both amplitudes): estimate/true residual between 0.76 and 1.26 (s = 128). The half-threshold margin is about 6 standard deviations; most explicit checks are of roots well inside the threshold.

### 3.3 The basis is twice QE's
Every stage except h_apply scales with m (rotate, reduced, projection) or m^3 (eigen). BAANDI's mean m is 375 fine / 559 coarse (capacity 6 x bands = 792); QE restarts at 2 x bands = 264. 82 percent (fine) and 96 percent (coarse) of the EVD m^3 comes from m >= 500.

### 3.4 Already at the minimum
Reduced matrix: one rectangle V^T HV[:, prefix..m] per iteration (the "two rectangles" in older notes is outdated); the only excess is the full recompute after the 5 contractions (6.6 GFLOP, 0.2 s), where the leading block is diag(theta) by construction. Contraction restarts reuse stored products (no H re-application, 0 rebuilds). DGKS skips the second pass on 93 of 156 appends.

## 4. Minimum flops for an equivalent algorithm (same trajectory: same steps, H columns, thresholds, certificate)
| stage | BAANDI now (GFLOP / s) | minimum, orthonormal basis kept | QE step structure on the same trajectory |
|---|---|---|---|
| rotate | 720 / 14.8 | about 520 / 10.7 (corrected roots 308, contractions 47, returns X 51 + H X for uncertified returns about 34, irreducible residual checks about 80) | 406 / 8.3 (update 308, last 98) |
| reduced | 224 / 5.8 | 218 / 5.6 | 218 / 5.6 (H border) |
| orthonormalize | 438 / 27.5 | 412 / 13.3 (one pass 313, MGS 99, all on dgemm) | 206 / 5.3 (S border instead) |
| eigen | m^3 sum 27.8e9 / 7.8 | 7.8 | 7.8 (generalized, same m) |
| total | 1382 / 55.8 | about 1150 / 37 | about 830 / 27 |
QE's measured 8.5 s is a further factor of about 3 below the last column: its basis is at most 264 columns, it has no coarse stage, and it never checks residuals (Ritz-change stopping). Of BAANDI's 47 s gap at this model, about 19 s is kernel dispatch and avoidable residual checks (3.1, 3.2), about 10 s is the QE step structure (non-orthonormal basis), and the rest is basis size, the coarse stage and the stopping rule.

## 5. Ranked reductions (si64 strict T=1, savings from the model; f1e-3 first, P48A-regime b1e-7 second)
| rank | change | expected saving | accuracy / risk | touches (do not repeat without this evidence) |
|---|---|---|---|---|
| 1 | Projection apply (and MGS-B cross-panel apply) through dgemm for any width: call `multiply_slices_wide_into` at site ProjectionApply and in the MGS-B panel apply, as ORTHO_BLOCK already does; new switch for A/B, default on | 15.3 s (orthonormalize 27.5 -> 12.2); 10.2 s at b1e-7; about 1.9 s at si32 strict | rounding class (dgemm panels instead of nano), thread-independent (fixed 256-row panels, per-panel tile bound unchanged); no algorithm change; low | MKL-R2 (intended coverage), DGEMM-SITES (8.72 vs 53.7 GFLOP), P48B (nano kept on MKL builds on the premise that dgemm covers these sites), MGS-B, P48D; NARROW_PROJ is the separate c <= 4 kernel and stays independent |
| 2 | B1 rotation only where needed: (a) narrow path whenever the explicit set is smaller than the band count (drop the half-band gate; restart, return and limit iterations keep the all-band path); (b) sketch margin from 0.5 to about 0.7 with 256 to 512 count-sketch rows (keeps the margin at about 6 standard deviations) | (a) 62 GFLOP, 1.3 s; (b) about 130 GFLOP, 2.7 s; together about 4 s at f1e-3, about 10 s at b1e-7 (rotate 20.6 -> about 10.5) | (a) rounding class; (b) more sketch disagreements are possible (a root above threshold skipped for one iteration changes the trajectory), accuracy unchanged because returns stay decided on the all-band explicit path plus certificate; medium-low | B1-SELECTED-ROOTS, B1-PART-B, ALGEBRA-STREAM-B1-MAC; PRT1 (stale flags) is not repeated: flags stay fresh each iteration; RITZ_FUSED is complementary (traffic only) |
| 3 | Smaller basis: retest the capacity multiplier (3x, 4x) after 1 and 2, coarse ladder stage first (it only supplies the start) | gross about 11 s if steps and H columns were unchanged (m-linear stages after 1 and 2 about 29 s x 25 percent, eigen 7.8 s x 55 percent); each 1 percent more H columns costs about 0.39 s, so break-even near +28 percent H columns | trajectory change, iteration counts may rise; accuracy unaffected; medium; generality set required | A13, A08, P31A, backlog B17. New evidence: B1 makes rotation proportional to the explicit set, the nano fallback inflated the per-column cost of m, and si64 eigen (7.8 s) is cubic in m with 82 to 96 percent from m >= 500; A13 was measured with 8 bands and all-root rotation |
| 4 | Non-orthonormal basis with generalized small problem (QE structure) | about 7 s after rank 1 (orthonormalize 12.2 -> S border 5.3) | high (A52 rank/rounding mechanism undetermined); gated | A16, B13 |
| 5 | Small items: diag(theta) after a contraction instead of the full reduced recompute (0.2 s); no H X rotation on returns that are certified (0.3 s); NARROW_PROJ and NARROW_REDUCED (69 of 156 fine appends have c <= 4, about 13 GB streamed; at most about 1 s); RITZ_FUSED (about 6 GB of H X traffic, at most 0.5 s); B5 for appends of 2 to 4 columns (narrow iterations) | about 1 to 2 s together | rounding class | NARROW_PROJ, NARROW_REDUCED, RITZ_FUSED, B5, ALGEBRA-STREAM rows |

Not proposed: range-selected dsyevr (MKL-EVD: slower at these orders), CholeskyQR2 (A52), residual norms from a Gram matrix (HV)^T HV (cancellation floor about sqrt(eps x m) x ||H||, a few 1e-6 Ha, too coarse for strict thresholds), ENERGY_IDENTITY (not algebra; rejected on dev08), SECOND_VECTOR (inactive while B1 is on; Mac signal unfavourable).

## 6. What dev08 should confirm first
1. DGEMM_TIMERS build (`--features mkl,dgemm-timers`), si64 strict f1e-3 seed 601 T=1: prediction for the current main, projection_apply gflops_total about 14 GFLOP against projection_overlap about 168 GFLOP. After rank 1 both should be equal and orthonormalize should fall by 13 to 16 s.
2. `BAANDI_DAVIDSON_STAGES=on` stage times at f1e-3 (the P48A times are b1e-7; the two trajectories differ: 92 width-1 iterations at b1e-7, 27 at f1e-3).
3. Rank 1 and rank 2 as separate same-binary switch A/Bs (seeds 601-610, si16/si32/si64 strict and standard, T=1 and T=8, k444s and mgo-40), under the adoption rule.

Limits: counts come from Mac pure trajectories; seconds are a four-rate model fitted to one relayed dev08 run, so individual stage seconds carry about 10 to 20 percent uncertainty; QE flops are not measured; the residual-distribution and sketch-ratio numbers are from seed 601 only.
