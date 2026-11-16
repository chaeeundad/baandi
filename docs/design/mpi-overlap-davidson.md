# Overlap-matrix (generalized reduced problem) Davidson on the M2 stick rows (ticket MPI-OVERLAP-DAVIDSON, 2026-11-17)

Switch `BAANDI_M2_LEAN_ALGEBRA=on|off` (default off; `DistributedTuning.lean_algebra`). Code: `crates/dft-solvers/src/distributed_davidson/workspace/lean.rs` (the iteration), `workspace.rs` (dispatch, start), `lean_davidson.rs` (the serial routines it reuses). Tests: `crates/dft-solvers/tests/distributed_lean_algebra.rs`, `crates/dft-engine/tests/m2_lean_algebra.rs`. Screens: section 8; dev08 request: docs/performance/requests/mpi-overlap-davidson.md. MPI against MPI only (comparison-rules.md). Labels: **O** Mac observation (simulated ranks), **C** exact count or read from the code, **E** estimate.

## 0. Plain summary
- What changes: on standard (Ritz-change) solves, the MPI solver stops orthonormalizing every new correction column against the basis by Gram-Schmidt. It appends the corrections as they are (unit length), keeps the overlap matrix S = V^T V of the basis beside the reduced Hamiltonian, and solves the generalized reduced problem through the inverse Cholesky factor of S (the serial lean solver's algebra, QE's structure). Strict solves are unchanged.
- Reductions: the orthonormal path spends about 25 to 60 small tree sums per Davidson iteration (normalization, basis projection, and a chain of about one reduction per appended column in the panel Gram-Schmidt); the new path spends two: one for the overlap border (with the residual norms riding along) and one for the Hamiltonian border. Measured on the warm si32 solve (O, C): 664 -> 257 tree sums per solve at np=8, 580 -> 120 at np=2 (where no fallback fires); the remaining count is mostly the start orthonormalization (M2_START) and the fallback stage 1 (section 4).
- What it could give: the panel Gram-Schmidt chain is the stage that does not shrink from np 4 to 8 (about a quarter of the np=8 warm si32 solve, mpi-remaining.md section 1). Mac screens (O): stage 0 warm solve -12 to -25 percent (si32 np 2/4/8, si64 np 4/8: -14 / -15 percent), cold -6 to -14 percent; stage 1 production SCF -6 to -10 percent of the wall (si32 np 4/8, si64 np 4), tree sums per run -66 percent, -89 percent together with the pending stack (mirror sticks, Gram start, kept start; si32 np=8 S -> S + L -10.3 percent); H columns, Davidson iterations, SCF steps unchanged everywhere. dev08 decides (request).

## 1. Background
FINDING (MPI-REMAINING-MAP): from np 4 to 8 the panel Gram-Schmidt of the correction block (49.8 -> 48.0 ms, warm si32), the replicated reduced eigenproblem and the wait after the rotation do not shrink: 58 percent of the np=8 warm si32 solve. QE never orthonormalizes the expanded basis: it forms the overlap border with one product per iteration and solves a generalized reduced problem (study notes qe-solvers-and-subspace.md, davidson-and-scf-thresholds.md; no QE source read). The serial lean solver (B-LEAN-DAVIDSON) already does this on one rank; the M2 solver needed an orthonormal basis. Review item F05 step (2): port the overlap S, bordered Cholesky and generalized reduced solve to the distributed solver, separately from the selected-root rotation (step 1, adopted as M2_SELECTED_ROOTS).

## 2. Distributed versus replicated
| quantity | where | size | how it is formed |
|---|---|---|---|
| basis V, H images HV | distributed (local stick rows) | n_loc x m | corrections appended without projection, scaled to unit norm after the S border |
| Ritz vectors / residuals of the selected roots, corrections | distributed | n_loc x k | local GEMM V C_s, (HV) C_s; preconditioner row-local |
| S border V^T E (m x k), with the k squared residual norms in an extra column | partial local, then ONE tree sum | m x (k + 1) | local GEMM V^T E, tree sum (site `lean_border`) |
| H border V^T (H E) | partial local, ONE tree sum | m x k | as the orthonormal path's `reduced_gram` |
| S, H_r | replicated | capacity x capacity, leading dimension capacity | written from the reduced borders (upper triangle defines) |
| T = L^{-1} (inverse Cholesky factor of S) | replicated | order m <= capacity (3 or 6 x bands) | serial `extend_factor` (panels of 64, pivot rule), small GEMM workspace (capacity x capacity, fixed panels) |
| A = T H_r T^T, eigenpairs, C = T^T Y | replicated | m x m | serial `reduced_border` (border only), dense solver or B5, serial `ritz_coefficients` |

Every replicated quantity is computed from tree-summed data by the same code on every rank (same GEMM workspace shape, thread-count independent), so every rank holds the same bits and takes the same branches (pivot decisions, fallbacks, convergence flags).

## 3. Reductions per iteration (C)
M2 stack (serial-order orthonormalization, fused panel Gram-Schmidt, selected roots), one Davidson iteration that appends k corrections:
| step | orthonormal path | lean algebra |
|---|---|---|
| H border | 1 (`reduced_gram`) | 1 (`reduced_gram`) |
| residual norms of the selected roots | 1 (`residual_norms`) | 0 (ride with the S border) |
| normalize corrections | 1 (`ortho_normalize`) | 0 (S diagonal) |
| project against the basis | 1, a second pass when DGKS fails (`ortho_projection`) | 0 (S border) |
| panel Gram-Schmidt | 2 per panel after the first (`panel_cross`) + about 1 per column with fused, up to 3 (`panel_dot`, `panel_norm`) | 0 |
| S border | | 1 (`lean_border`) |
| total | about 4 + k (+ panels, + second passes) | 2 |
| fallback stage 1 (pivot test failed) | | + 1 norm reduction + the in-block Gram-Schmidt of the block (1 + about 1 per column) + 1 S border |
| contraction | + 1 (`restart_defect`) | 0 (the S border covers the fresh basis X) |
| return | 1 (all-band residual norms) | 1 |

Payload: the S border carries m x (k + 1) doubles per call (bucket 65-512 to 513-4096), the H border m x k; total doubles per solve are about equal to the orthonormal path (warm si32: 1.19 -> 1.13 MB per rank at np 4 and 8, O).

## 4. Dependence and rank decisions (unchanged meaning)
Same constants and order as the serial lean solver:
- Pivot rule: every appended column's squared Cholesky pivot over its overlap diagonal must be at least `PIVOT_RATIO_MIN = eps / 1e-12` (the serial `extend_factor`, called unchanged). The pivot test runs on the reduced S before H is applied.
- Stage 1 (pivot test failed): E' = E - V T^T T (V^T E) by a local GEMM (V^T E is the S border already reduced; no new Gram), one norm reduction, columns with |E'_j| <= `DEPENDENCE_TOLERANCE` dropped (the classic rank floor on a unit input), the rest through the serial-order in-block Gram-Schmidt (its rank floor), then a new S border and pivot test.
- Stage 2 (stage 1 failed again, or nothing survived): V and HV are rotated once to the full Ritz basis V C (span unchanged, no H application) unless the basis is the fresh contraction basis; the surviving columns are projected with the classic projection against the basis (random block if nothing survives, as the orthonormal path); a failure of that block is an error.
- A fresh basis (start, contraction X, stage-2 rotation) that fails the pivot test is rebuilt by Gram-Schmidt and H application (counted `rebuilds`; expected zero, zero in every test and screen).
- Contraction restart: X = V C with C^T S C = I, taken without the 1e-12 orthogonality-defect check of the orthonormal path (the generalized path needs a factorable S, not an orthonormal X); its S comes with the next S border.

Tests (C): on 1 to 4 ranks the distributed solver takes exactly the serial lean solver's decisions (appends, stage-1 re-projections, stage-2 fallbacks), iterations and H columns, on the cubic, degenerate and skewed fixtures and on the serial lean test's ill-conditioned fixture (coupling 2e-6: stage 1 on both appends; coupling 0: stage 1 and stage 2 down to the random block).

## 5. Reproducibility and accuracy
- Fixed rank count: bit-reproducible (fixed-order tree sums; replicated small GEMMs on fixed panels). Ranks agree bit for bit. Different rank counts differ by rounding (the tree sums group the rows differently), as on the orthonormal path.
- Accuracy class unchanged: standard solves stop on the Ritz-value change as before; returned residual norms are formed explicitly from the returned X and H X (all bands); the returned X is S-orthonormal to about eps / pivot ratio (<= 1e-12 by the pivot rule; observed <= 1e-11 defect in every test). Strict solves (residual criterion, certificate) never take this path. No tolerance changed.

## 6. Rejected variant: one tree sum per iteration (C, O)
First version: S border, H border and residual norms in ONE tree sum after H was applied to the raw corrections (the pivot test then runs after the application; H E scaled with E). On the warm si32 solve (capacity 3) 15 of 22 appends fail the pivot test (stage 1), so 42 H columns per solve were applied to rejected blocks: H columns 272 -> 314 (+15 percent), alltoall calls +64 percent. The H application is the stage that scales; trading it for one reduction is a loss on dev08. Rejected before measurement on dev08; the two-reduction form applies H only to accepted blocks (H columns equal to the orthonormal path in every screen).

## 7. Interaction with the other switches
- `M2_WORKSPACE`: the lean algebra runs on the workspace path (implied, like M2_SELECTED_ROOTS); its replicated state (S, T, A, C, small GEMM workspace, border buffers) lives in the per-rank workspace, allocated once per (capacity, bands, pool), counted in `workspace_stats`. No allocation in the iteration except in the stage-2 random block (as the orthonormal path) and the allocating classic orthonormalization order (`M2_ORTHO=classic`, not the default).
- `M2_SELECTED_ROOTS`: the lean path always uses the selected-root rule (residuals of unconverged roots only; all bands at return and contraction), independent of that switch.
- `M2_PANEL_GS`: unused on the lean path except inside the fallback stages (in-block Gram-Schmidt of a rejected block, rebuild), where the switch applies as before.
- `M2_EIGEN=root`: ignored on the lean path (stage 2 needs all reduced eigenvectors on every rank); the reduced problem is solved replicated. `M2_B5=on` applies to A (its leading block is kept bit for bit by the border update).
- `M2_START=gram` (default since 0d87d6e): on the lean path the start block is accepted when its Gram matrix passes the pivot rule (instead of the 1e-12 defect), and the factor is kept (no second S reduction). `M2_START=ortho`: the start is orthonormalized as before (its panel Gram-Schmidt is then most of the remaining reductions of a solve).
- `M2_REUSE` (kept start, default on since 0d87d6e): the generalized path needs a factorable S of the kept X, not an orthonormal X: the kept start is taken when initial equals X (unchanged flag reduction) and the Gram of X (unchanged reduction) passes the pivot rule; the factor is kept. Counter `factored_starts`.
- `M2_LAZY`: unchanged (the return keeps the local rows; no gather).
- M1 path (`solve_distributed`): takes the same code when the switch is on (workspace path).

## 8. Screens (O, Mac pure build, simulated ranks, shared host at load 6 to 12; counters exact C)
Measured on 2963c86 defaults (greedy sticks, ortho start, no kept start); main 0d87d6e has since adopted mirror sticks, the Gram start and the kept start, i.e. the stack arm S below is the new default.
### Stage 0 (`scripts/fast_screen_dav.py --ranks P --b BAANDI_M2_LEAN_ALGEBRA=on`, f1e-3 seed 601, converged base density, 3 alternating pairs x 5 repeats, davidson_bench sha256 9f3833e8c557cea6; capacity rule 6 per band at np=2, 3 at np >= 4)
| cell | cold solve B/A (spread) | warm solve B/A (spread) | tree sums cold / warm A -> B | H columns, Davidson iterations | warm stage-1 re-projections | GEMM flop warm |
|---|---|---|---|---|---|---|
| si32 np=2 | 0.865 (5.0 %) | 0.747 (6.1 %) | 252 -> 98 / 580 -> 120 | equal (190 / 3, 244 / 10) | 0 of 9 appends | -25.3 % |
| si32 np=4 | 0.936 (8.6 %, not beyond) | 0.876 (1.9 %) | 224 -> 88 / 683 -> 261 | equal (190 / 3, 272 / 23) | 15 of 22 | -16.0 % |
| si32 np=8 | 0.897 (3.4 %) | 0.815 (2.6 %) | 214 -> 84 / 664 -> 257 | equal | 15 of 22 | -10.8 % |
| si64 np=4 | 0.926 (1.9 %) | 0.862 (1.6 %) | 482 -> 178 / 1489 -> 395 | equal (374 / 3, 603 / 24) | 14 of 23 | -19.9 % |
| si64 np=8 | 0.936 (4.4 %) | 0.853 (11.4 %, not beyond) | 426 -> 158 / 1354 -> 375 | equal | 14 of 23 | -16.1 % |
Occupied eigenvalues moved by at most 3.3e-14 Ha; eigenvalues reproducible within each arm; no rebuild, no stage 2, no wasted H column. Stage clocks (si32 np=8, `BAANDI_DAVIDSON_STAGES=on` on both arms, medians of 3, ms per rank, never used for timing), warm A -> B: orthonormalize 67.1 -> 7.3 (stage-1 fallbacks only), reduced_gram 16.6 -> 43.9 (both borders and the wait for the slowest rank that `residual` held before), residual 11.3 -> 1.6, reduced_eigen 25.2 -> 31.2 (factor extension and border of A, replicated), rotate 12.7 -> 15.1, apply 82.7 -> 86.9, tree_sum waits 75.9 -> 47.9; the four algebra stages 120.2 -> 84.0 ms.

### Stage 1 (`scripts/fast_screen_mpi.py`, production standard with verification and ladder, si32/si64 f1e-3 seed 601, 3 pairs, m2_entry_cost sha256 1b73aa48f6e33ebc)
| cell | arm A -> B | wall A -> B (spread) | tree sums per rank | per site (A -> B) | dE |
|---|---|---|---|---|---|
| si32 np=4 | defaults -> L | 3.684 -> 3.409 s (-7.5 %, 0.034 s) | 6120 -> 2103 (-65.6 %) | panel_dot 4632 -> 1426, panel_norm 454 -> 70, ortho_projection 113 -> 0, residual_norms 127 -> 22, lean_border 0 -> 127, reduced_gram 127 -> 105 | 1.4e-14 Ha |
| si32 np=8 | defaults -> L | 3.523 -> 3.154 s (-10.5 %, spread 0.47 s: within) | 6016 -> 2033 (-66.2 %) | as np=4 | 0 (energy bits equal) |
| si64 np=4 | defaults -> L | 13.98 -> 13.12 s (-6.1 %, 0.58 s) | 9118 -> 3014 (-66.9 %) | panel_dot 6543 -> 2027, panel_cross 1302 -> 400 | 0 (energy bits equal) |
| si32 np=8 | S -> S + L (S = mirror sticks, Gram start, kept start) | 3.162 -> 2.835 s (-10.3 %, 0.23 s) | 4494 -> 511 (-88.6 %) | every panel and ortho site -> 0, lean_border 105, restart_defect 31 -> 20 | 0 (energy bits equal) |
SCF, Davidson and H counts unchanged in every cell; density reproducible within each arm; 22 lean solves per rank per si32 run carry every standard iteration (127 iterations: 127 S borders + 105 H borders), no fallback, no rebuild in any SCF screen. With `M2_START=ortho` (default) the start orthonormalization is the remaining chain (panel_dot 1426 per rank); with the stack's Gram start (factored by the pivot rule, `factored_starts` 22) it disappears: 511 tree sums per rank for the whole si32 run.
