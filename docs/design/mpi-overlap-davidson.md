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
- Pivot rule: every appended column's squared Cholesky pivot over its overlap diagonal must be at least `PIVOT_RATIO_MIN = eps / 1e-12`. This is a local dependence test, **not** a bound on the conditioning of the whole overlap. The shared serial `extend_factor` also bounds weighted inverse-factor growth (section 9). Both tests run on the reduced S before H is applied.
- Stage 1 (pivot test failed): E' = E - V T^T T (V^T E) by a local GEMM (V^T E is the S border already reduced; no new Gram), one norm reduction, columns with |E'_j| <= `DEPENDENCE_TOLERANCE` dropped (the classic rank floor on a unit input), the rest through the serial-order in-block Gram-Schmidt (its rank floor), then a new S border and pivot test.
- Stage 2 (stage 1 failed again, or nothing survived): V and HV are rotated once to the full Ritz basis V C (span unchanged, no H application) unless the basis is the fresh contraction basis; the surviving columns are projected with the classic projection against the basis (random block if nothing survives, as the orthonormal path); a failure of that block is an error.
- A fresh basis (start, contraction X, stage-2 rotation) that fails the pivot test is rebuilt by Gram-Schmidt and H application (counted `rebuilds`; expected zero, zero in every test and screen).
- Contraction restart: X = V C with C^T S C = I in exact arithmetic. The actual Gram defect of X is checked against 1e-12 using the leading block of the next S border. A failed check rebuilds X before any classical projection. No new collective is needed.

Tests (C): on 1 to 4 ranks the distributed solver takes exactly the serial lean solver's decisions (appends, stage-1 re-projections, stage-2 fallbacks), iterations and H columns, on the cubic, degenerate and skewed fixtures and on the serial lean test's ill-conditioned fixture (coupling 2e-6: stage 1 on both appends; coupling 0: stage 1 and stage 2 down to the random block).

## 5. Reproducibility and accuracy
- Fixed rank count: bit-reproducible (fixed-order tree sums; replicated small GEMMs on fixed panels). Ranks agree bit for bit. Different rank counts differ by rounding (the tree sums group the rows differently), as on the orthonormal path.
- Accuracy class unchanged: standard solves stop on the Ritz-value change as before; returned residual norms are formed explicitly from the returned X and H X (all bands); the returned block must pass a maximum-entry Gram-defect check of 1e-12. The pivot ratio alone does not imply this bound. Section 10 evaluates the check from the current measured S and C, retaining a direct Gram for cancellation-sensitive or borderline cases. Failure rebuilds and reconfirms Ritz changes, or returns an error at the iteration limit. Strict solves (residual criterion, certificate) never take this path. No tolerance changed.

## 6. Rejected variant: one tree sum per iteration (C, O)
First version: S border, H border and residual norms in ONE tree sum after H was applied to the raw corrections (the pivot test then runs after the application; H E scaled with E). On the warm si32 solve (capacity 3) 15 of 22 appends fail the pivot test (stage 1), so 42 H columns per solve were applied to rejected blocks: H columns 272 -> 314 (+15 percent), alltoall calls +64 percent. The H application is the stage that scales; trading it for one reduction is a loss on dev08. Rejected before measurement on dev08; the two-reduction form applies H only to accepted blocks (H columns equal to the orthonormal path in every screen).

## 7. Interaction with the other switches
- `M2_WORKSPACE`: the lean algebra runs on the workspace path (implied, like M2_SELECTED_ROOTS); its replicated state (S, T, A, C, small GEMM workspace, border buffers) lives in the per-rank workspace, allocated once per (capacity, bands, pool), counted in `workspace_stats`. No allocation in the iteration except in the stage-2 random block (as the orthonormal path) and the allocating classic orthonormalization order (`M2_ORTHO=classic`, not the default).
- `M2_SELECTED_ROOTS`: the lean path always uses the selected-root rule (residuals of unconverged roots only; all bands at return and contraction), independent of that switch.
- `M2_PANEL_GS`: unused on the lean path except inside the fallback stages (in-block Gram-Schmidt of a rejected block, rebuild), where the switch applies as before.
- `M2_EIGEN=root`: ignored on the lean path (stage 2 needs all reduced eigenvectors on every rank); the reduced problem is solved replicated. `M2_B5=on` applies to A (its leading block is kept bit for bit by the border update).
- `M2_START=gram` (default since 0d87d6e): on the lean path a complete start block is accepted when its Gram matrix passes the pivot and inverse-growth rules, and the factor is kept (no second S reduction). A partial start is accepted unchanged only when its Gram defect is at most 1e-12; otherwise it is orthonormalized before random completion. `M2_START=ortho`: the start is orthonormalized as before (its panel Gram-Schmidt is then most of the remaining reductions of a solve).
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

## 9. LEAN-ALGEBRA-HARDEN (2026-11-17; implementation ad7a5c0)

This section records the first hardening implementation. Section 10 supersedes
its unconditional direct return Gram and additional status-double payload;
the thresholds, recovery rules and regression evidence remain in force.

The review's bidiagonal chain (`v0 = e0`, `vj = sqrt(0.99) e[j-1] +
0.1 e[j]`) passes every original pivot test but amplifies overlap rounding.
For each newly computed row of T, the shared factor extension computes
`g_i = sum_j |T_ij| sqrt(S_jj)`. A componentwise Gram perturbation bounded by
`eps sqrt(S_ii S_jj)` is amplified by at most `eps g_i g_j`. The guard rejects
`eps g_i^2 > sqrt(eps)`, i.e. loss of half the working digits
(`g_i > eps^(-1/4) = 8192` for binary64). This is a conservative early warning,
not a forward-error theorem for dot products of arbitrary length. The actual
1e-12 return/contraction checks enforce the orthogonality contract. No existing
pivot, dependence, Ritz-change or residual tolerance is relaxed. Rejected
extensions use the existing stage-1/stage-2 rules; rejected fresh bases rebuild.
A rebuild with fewer columns than required bands returns a replicated error
before the reduced solve can index missing columns.

On the row-operator interface, operator errors and nonfinite outputs are recorded
locally, with the failed output replaced by zeros. A solver-only communicator
packs one status double with each existing solver tree sum, using buffers
reserved at solve entry. The operator retains the underlying communicator, so
its own collectives are unaffected. The next solver sum makes every rank return
the same error before using its numerical payload. This covers both workspace
and allocating row solvers. It requires the operator to complete its own
collectives before returning, as in the operator contract; it cannot recover a
rank that exits inside an unfinished operator collective. Initial local row and
diagonal errors ride with the existing row-count exchange (or the layout-cache
sum). There is no new collective on the ordinary path. Payload sizes increase;
reduction histogram buckets can therefore change even though call counts do not.

Empty row partitions retain the replicated panel column count explicitly;
deriving it from `local_elements / local_rows` loses that information on an
empty rank and previously produced mismatched reduction lengths during repair.

Validation sources: `tests/lean_harden.rs`, the factor-chain unit test, the
forced-stage-2-after-contraction unit test, and
`tests/distributed_lean_algebra.rs`. The independent oracle diagonalizes the
original dense operator, never the solver's reduced overlap problem. The
original loose buffer policy is retained: dense eigenvalue errors are bounded
by 1e-10 Ha on required bands and 1e-9 Ha on buffer bands; independently computed
true residuals of **all** bands are bounded by 1e-4 Ha. This is deliberately
separate from a strict residual certificate: the unchanged skew fixture has
true residuals up to 1.67e-5 Ha under its 1e-9 Ha buffer Ritz-change policy.
The exact-spectrum diagonal and forced-contraction tridiagonal regressions
bound every required true residual by 1e-6 Ha and eigenvalue error by 1e-10 Ha.
No convergence flag is used as a substitute for either check.

Baseline for bit comparisons: `0fb09c0edf189f6ced15ae2d4fca4091f5a364aa`,
Apple Silicon pure Rust, offline debug. The test prints FNV-1a-style hashes of
all eigenvalue, vector and reported-residual bits, plus iteration, H-column,
append/fallback and tree-sum counts. Hashes are evidence for this backend only,
not portable golden values across BLAS/EVD implementations. Before/after records
are below; no SCF energy or dev08 timing claim is made by this correctness task.

```text
BASELINE cubic cap 3 ranks 1: 0aaf10ee411555c8 801b4171736bfe14 (5, 3, 27, 18, 6, 1, 0) sums 31
BASELINE cubic cap 3 ranks 2: 0f4eb50c65f441b9 1a0c1819d4f78ac2 (5, 3, 27, 18, 12, 2, 0) sums 62
BASELINE cubic cap 3 ranks 3: d88d10327b6f429b db336239dddd6c7a (5, 3, 27, 18, 18, 3, 0) sums 93
BASELINE cubic cap 3 ranks 4: da5c184ac2dd0e89 4c53dc228c465e08 (5, 3, 27, 18, 24, 4, 0) sums 124
SERIAL-BASELINE cubic cap 3: 31fec50e8f18d4f1 (5, 27, 4, 1, 0)
BASELINE cubic cap 6 ranks 1: 8e2966f46bb91b9e 3cdf69b7c7bafea7 (5, 3, 27, 18, 6, 1, 0) sums 31
BASELINE cubic cap 6 ranks 2: ff33298b31d74997 a09f0903c9a1d5f7 (5, 3, 27, 18, 12, 2, 0) sums 62
BASELINE cubic cap 6 ranks 3: 949432f178454588 5a179a192d3425cc (5, 3, 27, 18, 18, 3, 0) sums 93
BASELINE cubic cap 6 ranks 4: 38af73844d42e4c3 1a38378436d16202 (5, 3, 27, 18, 24, 4, 0) sums 124
SERIAL-BASELINE cubic cap 6: 07f170aa53ffa99b (5, 27, 4, 1, 0)
BASELINE degenerate cap 3 ranks 1: 975eeb5ba2cbce5a 2791c32b75f1c48b (4, 4, 24, 18, 6, 0, 0) sums 18
BASELINE degenerate cap 3 ranks 2: 3d13b470b3461a91 f52b2b67e737cc2c (4, 4, 24, 18, 12, 0, 0) sums 36
BASELINE degenerate cap 3 ranks 3: 30e26e96cca0ea2c a5d945be1ab291de (4, 4, 24, 18, 18, 0, 0) sums 54
BASELINE degenerate cap 3 ranks 4: 82cc6e775155625b 066106ab6f356615 (4, 4, 24, 18, 24, 0, 0) sums 72
SERIAL-BASELINE degenerate cap 3: a3424a575f933baf (4, 24, 3, 0, 0)
BASELINE degenerate cap 6 ranks 1: 87c35157bcdf0474 dfbba373d65f2e5e (4, 4, 24, 18, 6, 0, 0) sums 18
BASELINE degenerate cap 6 ranks 2: fc829e758052d8f3 11d03e7855ae3ee6 (4, 4, 24, 18, 12, 0, 0) sums 36
BASELINE degenerate cap 6 ranks 3: cab70bfcdd48b751 6ce20f7c326df057 (4, 4, 24, 18, 18, 0, 0) sums 54
BASELINE degenerate cap 6 ranks 4: 63b2831caa63b833 51537d86d7514d4b (4, 4, 24, 18, 24, 0, 0) sums 72
SERIAL-BASELINE degenerate cap 6: 6c11734df57325ce (4, 24, 3, 0, 0)
BASELINE skew cap 3 ranks 1: 8a4b1624217198af caa45072667de0b0 (22, 8, 103, 47, 28, 0, 0) sums 62
BASELINE skew cap 3 ranks 2: 85185de0cbe8e53b d3bbb00510fd6dd9 (22, 8, 103, 47, 56, 0, 0) sums 124
BASELINE skew cap 3 ranks 3: c0d38a23c1266943 fef7df4cbad8076a (22, 8, 103, 47, 84, 0, 0) sums 186
BASELINE skew cap 3 ranks 4: 843a390731e7185d 6a2860cc9f742b00 (22, 8, 103, 47, 112, 0, 0) sums 248
SERIAL-BASELINE skew cap 3: b48272857239eb90 (22, 103, 21, 0, 0)
BASELINE skew cap 6 ranks 1: 3aa7148e65b12529 2efe667fe63281a4 (9, 6, 64, 44, 13, 0, 0) sums 32
BASELINE skew cap 6 ranks 2: 6528e7b4838a7f8e 6941d1b2325b083f (9, 6, 64, 44, 26, 0, 0) sums 64
BASELINE skew cap 6 ranks 3: c8d03164a56b4485 16e0f14fe577ac7d (9, 6, 64, 44, 39, 0, 0) sums 96
BASELINE skew cap 6 ranks 4: c1a78814a5e71296 b6d3a78df2488573 (9, 6, 64, 44, 52, 0, 0) sums 128
SERIAL-BASELINE skew cap 6: 838277921965fd57 (9, 64, 8, 0, 0)
```

Verification on the patched tree: `cargo test -p dft-solvers --offline` exited
0 (133 passed, 25 ignored across 27 test binaries/doc-test groups).
`cargo clippy -p dft-solvers --offline --all-targets -- -D warnings` exited 0.
The original-source adversarial run exited 101 (8 failed in 3.01 s); the
patched integration run passed 9/9 in 0.01 s. The factor-growth, forced
contraction fallback, and nonempty rank-loss rebuild unit regressions pass.
The dense-oracle fixture test completes within the approximately 10 s debug
budget (whole integration binary 6.43 s in the recorded run).

`scripts/check.sh dft-solvers` exited 0: build, test, fmt and clippy all PASS.

## 10. Hardening overhead follow-up (2026-11-17)

Motivation, **relayed by the orchestrator**, not measured here: on dev08 the
hardening added 2.15% / 2.36% GEMM FLOP (si32 / si64, np=8), 1.1–1.4% GEMM
calls, approximately 5% reduction payload, and 0.65% si64 wall. Safety checks
remain enabled; this change removes their unnecessary large-array work.

### Return check and its exceptional direct path

Both lean solvers now compute `G = C^T S C`, with `C` the actual current Ritz
coefficients and `S` the measured overlap of their basis. Small fixed-order
compensated dot products reuse `LeanScratch.p` and `.z`, without allocation or
GEMM dispatch. They cost O(m^2 B + m B^2) scalar operations instead of an
N-by-B Gram. These scalar operations are **not included** in the GEMM FLOP
counter: identical GEMM counters do not mean the safety check does zero work.
No eigenvalue, vector, residual, factor or convergence operand is changed.

The defect and the decision to use a direct Gram use replicated matrices only:
all ranks take the same branch. The acceptance threshold remains **1e-12**.
Let `d = max |G-I|`, `a = max (|C|^T |S| |C|)`, and
`g = max_j sum_i |C_ij| sqrt(S_ii)`. The direct check is retained exactly when
an input/intermediate is nonfinite or

```
d + 4 eps (a + g^2) >= 1e-12.
```

Compensated summation removes the first-order m factor from reduced dot-product
summation error. The positive `a` and `g^2` terms screen product rounding and
amplification of overlap error; they are a conservative sensitivity screen,
not a theorem about arbitrary corrupt/stale S. A naive m-scaled screen was
tried first: it sent ordinary warm solves back to the full Gram and did not
remove the cost, so it is not the delivered implementation. The final screen
is based on the same general rounding quantities for every system, without a
material, band-count, rank-count or benchmark-specific branch.

**Overlap freshness is a required invariant.** A contraction or stage-2 basis
rotation resets `factored` to zero; its next border measures the new basis
before solving/returning. A rebuild does the same. Kept/Gram starts use their
measured start Gram. Thus no path substitutes I or reuses the old basis's S
for a contracted basis. S still contains floating-point Gram/normalization
error; large coefficient cancellation or a defect near the threshold is the
condition for retaining the direct check. In that case MPI folds the original
full Gram into the existing residual reduction, with exactly the old 1e-12
check and rebuild. Ordinary returns reduce residual norms alone, with their
original arithmetic. Contraction checks and the inverse-factor growth guard
are unchanged. The latter scans new triangular entries in O(to^2-from^2),
at most O(w^2) for width w, and makes no additional GEMM call.

### Error agreement without an extra scalar

A rank that has completed the operator's collectives but received an operator
error/nonfinite output retains the local failure flag and substitutes a quiet
NaN for the first scalar of the next nonempty solver reduction. IEEE addition
propagates the sentinel to every rank; the wrapper rejects the reduced payload
before it is used. The failed numerical payload was already unusable. On the
normal path the wrapper forwards the original slice directly, without packing
or copying. Empty reductions are forwarded and the flag is retained until a
nonempty reduction; operator collectives still use the underlying communicator.
The startup row-count exchange and layout-cache sum similarly use NaN only on
invalid input. Counts, bytes and numerical data are unchanged on valid input.
The existing error/NaN/shape timeout tests still cover first and later applies.

### Measurements: pure Rust, simulated ranks, counters only

Built with `cargo build --release -p dft-cli --offline --example davidson_bench`.
Three binaries: original numerical sources from `0fb09c0`, hardening `ad7a5c0`,
and this uncommitted follow-up. The original binary was built by temporarily
restoring the four solver source files from `0fb09c0`; all other numerical
sources are unchanged between those commits. Sources were restored afterwards.
SHA-256, including the actual binary build metadata:

| Artifact | SHA-256 |
|---|---|
| Original binary `/tmp/lean-follow-base-bench` | `03073477bb1a64fd712d205625a1f8423ff6d69daf088ae9e07754895043c84e` |
| Hardening binary `/tmp/lean-follow-hard-bench` | `7cbb740270b874b731936a785ecc13513d449996fcc13bf3e900f4a5aa878469` |
| Follow-up binary `/tmp/lean-follow-final-bench` | `fc7f393235168013243760ec03d206a3b03e9fd694a652c02a3bc274cd929fe5` |
| Input `/tmp/lean-follow-si32.toml` | `a34cba3f80788eb18acd238db099e4c47808dddf6f304a136d1bb77a8bfd2e2d` |
| Density `/tmp/lean-follow-si32-density.bin` | `9491950d470503bf3e5c0efbb2f0991912fbcfea1904bb7eeb9f3a44dccb24af` |
| External Si UPF | `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217` |

Input generated with `python3 -I scripts/fast_screen_dav.py --upf
../pseudo/Si.pz-vbc.UPF --case si32 --emit-input /tmp/lean-follow-si32.toml`:
seed 601, fractional displacement 0.001, one thread per simulated rank,
11,877 plane waves, 68 bands, 72 x 72 x 36 FFT. The supplied density is the
existing converged `fbdav-work/si32-after/density.bin` (its adjacent results
record convergence); its input `si32-gamma-s601-a0.001-T1.toml` is byte-identical
to this input. Both cold and warm use this fixed base density. Command for
each binary and P=4,8:

```
BINARY --input /tmp/lean-follow-si32.toml --density /tmp/lean-follow-si32-density.bin --sim-ranks P --repeats 2
```

Counters are per rank, per timed map containing one Davidson solve and its
unchanged surrounding map work (two repeats identical). `red_payload`
is doubles (`dist.reductions.doubles`); `coll_bytes` sums every `coll.*.bytes`.

| np | kind | version | gemm_calls | gemm_gflop | red_payload | coll_bytes |
|---|---|---|---:|---:|---:|---:|
| 4 | cold | original = follow-up | 44 | 0.949565804 | 73373 | 37450848 |
| 4 | cold | hardening | 45 | 0.977025428 | 78005 | 37487936 |
| 4 | warm | original = follow-up | 396 | 1.441010662 | 119740 | 47111520 |
| 4 | warm | hardening | 397 | 1.468470286 | 124524 | 47149824 |
| 8 | cold | original = follow-up | 44 | 0.489712118 | 73373 | 22106240 |
| 8 | cold | hardening | 45 | 0.503441930 | 78005 | 22143360 |
| 8 | warm | original = follow-up | 396 | 0.745362895 | 119740 | 27343232 |
| 8 | warm | hardening | 397 | 0.759092707 | 124524 | 27381568 |

**Every** benchmark counter, all eigenvalues, residuals and occupied band-energy
bits match the original in all four ordinary cases, not only the four counters
in the table. Eigenvalue FNV: np4 cold `0xf58b72257aba54b9`, warm
`0xdf5548b326d935c4`; np8 cold `0x772dbd8d468fd617`, warm
`0x4a3e106bfd5aadc9`. The 24 synthetic fixture `WORK` records likewise match
original GEMM calls/FLOP/payload exactly, and all 30 prior result-bit/decision
records match. Local raw logs: `/tmp/lean-follow-{base,hard,final}-rho-np{4,8}.log`
and `/tmp/lean-follow-{base,final}-fixtures.log`.

### Deliberately retained cost outside the ordinary cases

Also ran the same input **without** `--density` (uniform density, two follow-up
repeats). Both cold solves match original counts/bits. Both warm solves trigger
the cancellation/borderline screen and retain the full Gram: removing this cost
would weaken the check. At np4 the original hardening already rebuilds once
(32 iterations / 753 H columns versus original 30 / 617); the follow-up keeps
that repair and its eigenvalue FNV `0x85440bf8b00d238f`. This is not introduced by
the follow-up. At np8 no rebuild is needed and original eigenvalue bits remain,
but the direct check is retained because the reduced check is sensitive.

| uniform warm | original | hardening | follow-up |
|---|---:|---:|---:|
| np4 GEMM calls | 435 | 479 | 478 |
| np4 GEMM FLOP | 4119666493 | 4762850445 | 4735390821 |
| np4 reduction doubles | 361898 | 420949 | 416176 |
| np4 collective bytes | 107690192 | 127903304 | 127865088 |
| np8 GEMM calls | 435 | 436 | 436 |
| np8 GEMM FLOP | 2142353840 | 2156083652 | 2156083652 |
| np8 reduction doubles | 361898 | 366592 | 366522 |
| np8 collective bytes | 59997424 | 60035040 | 60034416 |

No wall-time speedup is claimed from these shared-host simulations. The
orchestrator should repeat its dev08 measurement on this change; no remote
server was accessed or modified.

Follow-up verification: `cargo test -p dft-solvers --offline` exit 0
(134 passed, 25 ignored across 27 groups); all nine `lean_harden` integration
tests pass unchanged, including their direct 1e-12 Gram checks. The added
unit test covers nonidentity S, the near-threshold direct branch, and
cancellation-sensitive coefficients. `scripts/check.sh dft-solvers` exit 0:
build, test, fmt, clippy all PASS. Only owned files changed; uncommitted.
