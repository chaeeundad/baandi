# History-preserving restart and same-H continuation (ticket BREAKTHROUGH-1-RESTART, 2026-11-17)

Source: review "BAANDI breakthrough performance review 2026-11-17", sections 3.2, 5 (R1), 13.2, 14 (E1, E2), 15; backlog B05. Ledger rows R1A-SAME-H-CONTINUE and RESTART-KEEP-R1B (section "Breakthrough restart"). Labels: **C** exact count or read from the code, **O** Mac observation (simulated ranks are threads of one process), **E** estimate. MPI is compared with MPI only, threads with threads only (comparison-rules.md).

## 0. Plain summary
- R1-a (same-H continuation of the first-step "standard redo"): **not implemented, ceiling zero.** The redo fires in 0 of 9 default runs (si8, si16, si32, seeds 601 to 603, cutoff ladder on) and 0 of 6 runs with the ladder off (si8, si32): the first map's stop estimate is never below its tolerance times the electron count. The kept-start counter `unchanged_potential` (the only way a solve sees the same H twice) is 0 in all 15 runs. Nothing to save in these workloads; reopen only for a workload that shows redos (a restart from a converged density, a relaxation step with density transfer).
- R1-b (keep previous search directions at a contraction restart, `BAANDI_RESTART_KEEP`): **implemented, REJECT.** Keeping directions does not remove work. On the only case with contractions (the MPI capacity rule, 3 columns per band) the H columns never fall: stage 0 (si32 np=8 simulated, warm solve from the start density, 4 contractions) +0.3 % (r = 1) to +6.5 % (r = B/4); stage 1 (whole production SCF, si32 np=8) -0.6 % (r = 1) and +0.6 % (r = B/4), Davidson iterations +4.5 % / +1.5 %, rotation FLOPs -1.2 % / +8.4 %, sum m^3 +2.7 % / +5.9 %. The kept directions take room from the next correction blocks, so the basis fills sooner and contractions become more frequent (11 -> 13 per SCF at r = B/4; 4 -> 6 or 7 per warm solve). The threads path (serial lean, capacity 6 per band) has **no contraction at all** in the default si8/si16/si32 SCF: the switch is a no-op there (energy bit-identical).

## 1. R1-a: where the redo is and why it does not fire
`ScfSession::advance` (session.rs) redoes the first SCF step at the same input density with the ratchet tolerance when `iteration == 1` and the standard stop estimate E_H[rho_out - rho_in] is already below tolerance x N_e. With the cutoff ladder (default) the fine first step starts from the prolonged coarse density with a seeded tolerance; without the ladder it starts from the atomic (or uniform) density at the 5e-3 Ha cap. In neither case is the first map near its fixed point, so the condition is false. Counted from results.json (`convergence.history[].standard_first_iteration_redo`, `convergence.work.reuse_psi.unchanged_potential`):

| runs | redo | unchanged-potential kept starts |
|---|---|---|
| si8/si16/si32, seeds 601-603, f1e-3, T=1, defaults (ladder on) | 0 of 9 | 0 |
| si8/si32, seeds 601-603, `BAANDI_CUTOFF_LADDER=off` | 0 of 6 | 0 |

The continuation contract (keep basis, H V, S, T, reduced matrices and the last Ritz-value change; re-evaluate the flags at the tighter tolerance; arm only when the potential is bit-equal and the operator was not rebuilt) is written down here for the case a workload with redos appears; it was not coded because its saving is zero in every measured run (no gold plating).

## 2. R1-b: the algorithm (C)
Lean Davidson with basis V (m columns, not orthonormal), H V, overlap S = V^T V (replicated, from the existing border), Ritz coefficients C (X = V C, C^T S C = I).
- History: at every iteration without a contraction, the coefficients D of the selected (unconverged) roots' Ritz vectors (`ritz_coefficients` with the selected roots, already formed: a copy of m x k numbers); after a contraction or a stage-2 Ritz-basis rotation the history is the basis columns of those roots (unit vectors); a rebuild or a new solve invalidates it.
- At a contraction (m + k > capacity): for the roots selected now that the history holds, in ascending order: D (zero below the history prefix), A = C^T S D, G = D - C A, a second pass G -= C (C^T S G), then |P|^2 = G^T S G. Kept when |P| > DEPENDENCE_TOLERANCE |Y| (the classic rank floor on a unit input) and the squared pivot against the directions kept before, over |P|^2, is at least PIVOT_RATIO_MIN = eps / 1e-12 (the lean pivot rule); the kept columns are S-orthonormalized in coefficient space (G L^{-T}). At most min(r, capacity - bands - k) directions: room for this iteration's correction block is always left.
- X and P come from ONE rotation with [C | G] (P = V G, H P = (H V) G): no H application, no FFT, no projector product. Extra rotation cost 2 x 2 n m r' flops.
- The new basis [X | P] is factored with the same pivot rule (X first, then P); any failure falls back to the plain restart from X alone (serial: truncate; distributed: truncate, re-append the raw corrections and repeat the border, one extra reduction, counted `keep_fallbacks`; zero in every test and screen).
- r: `BAANDI_RESTART_KEEP=off|1|2|4|8|b/8|b/4|b/2` (a count or a fraction of the band count, at least one). Capacity rule unchanged (3 per band at np >= 4, 6 otherwise). No element or material rule.

Distributed (C): P and H P are local rows from the same local GEMM as X; D, A, G, the Cholesky of G^T S G are small replicated matrices computed from the reduced S in a fixed order, so every rank takes the same decisions. **Extra reductions per restart: zero** in the normal case (the S border after a contraction already covers the whole fresh basis; it is m x (m + 1) doubles wider by r' columns); one per R1-b fallback.

Paths: serial lean (`lean_davidson.rs`, threads path) and the distributed lean algebra (`distributed_davidson/workspace/lean.rs`, `BAANDI_M2_LEAN_ALGEBRA=on`). The default orthonormal M2 workspace path is not covered (its contraction needs an orthonormal X with the 1e-12 defect rule; P would need an explicit orthonormalization against X); not pursued after the lean result. `BAANDI_LEAN_PROJ_CACHE=on` (default off): a contraction that keeps directions invalidates the projection cache (the return then computes P^T X directly).

## 3. Tests (C, `crates/dft-solvers/tests/restart_keep.rs`)
- Synthetic fixtures cubic, exactly degenerate, clustered (1e-6 Ha splittings), skewed; capacities 3 and 6; Ritz-change tolerances 1e-11 and 1e-9; r = 1, B/4, B/2; cold and warm solves. Returned vectors orthonormal to 1e-11, reported residuals equal the true residuals to 1e-11, all roots converged, no rebuild.
- Without a contraction (capacity 12) the switch on is bit-identical to off (serial and distributed, 1 and 3 ranks).
- Distributed (1 to 4 simulated ranks) takes exactly the serial lean decisions with the switch on: iterations, H columns, contractions, kept directions, fallbacks; every rank the same bits; a fixed rank count bit-reproducible; border reductions at most 2 per iteration plus the counted fallback stages.
- Accuracy finding: the Ritz-change stop is not an error bound for either path. On the skew fixture at capacity 3 the plain restart itself misses the strict and standard classes against the dense eigenvalues (1.8e-10 Ha at tolerance 1e-11, 1.1e-8 Ha at 1e-9), and the kept directions move the iteration where the stop fires (r = B/4, 1e-11: 4.7e-10 Ha; r = B/2, 1e-9: 9.0e-8 Ha; other settings equal or better). The test bounds the error by max(class, 10 x plain error) and prints the table. Synthetic H columns (skew, capacity 3, cold + warm): plain 103+47 (1e-11) / 80+35 (1e-9); r = 1: 119+47 / 95+34; r = B/4: 111+46 / 88+35; r = B/2: 102+46 / 78+35. Every other fixture: one contraction or none, identical counts.
- Debug run 4 s (smoke subsets); full grids `#[ignore]`: `cargo test --release -p dft-solvers --test restart_keep -- --include-ignored --nocapture --test-threads=1`.

## 4. Stage 0 (O, Mac pure build, simulated ranks np=8, si32-gamma f1e-3 seed 601, `BAANDI_M2_LEAN_ALGEBRA=on` on both arms, capacity 3 per band)
Finding the restarts: with the converged base density (the documented late-SCF warm solve) the warm solve fills the basis exactly to capacity and returns (H 204 = 3 x 68, 0 contractions) and the cold solve has none either: R1-b is bit-identical there. With the start density (step-2-like warm solve, ratchet tolerance 1e-8 Ha) the warm solve has 4 contractions in the baseline; that is the stage-0 case.

| r | contractions | kept directions | H columns | Davidson iterations | rotation FLOP | GEMM FLOP | sum m^3 | tree sums | occupied eigenvalue moved |
|---|---|---|---|---|---|---|---|---|---|
| off (A) | 4 | 0 | 617 | 30 | 6.68e8 | 2.14e9 | 9.94e7 | 104 | |
| 1 | 6 | 6 | 619 (+0.3 %) | 25 (-17 %) | -2.2 % | +3.5 % | -2.5 % | -22 % | 1.8e-6 Ha |
| B/8 | 6 | 44 | 618 (+0.2 %) | 30 | +6.5 % | +8.9 % | -4.2 % | -8 % | 1.6e-6 Ha |
| B/4 | 7 | 106 | 657 (+6.5 %) | 30 | +26 % | +26 % | +19 % | +12 % | 1.3e-6 Ha |
| B/2 | 7 | 208 | 623 (+1.0 %) | 24 (-20 %) | +45 % | +41 % | -5.5 % | +61 % | 2.0e-6 Ha |

Three-pair screen r = 1 (`fast_screen_dav.py --ranks 8 --expect h_columns:down`, davidson_bench sha256 02952a546d25388e, load 5 to 6): cold 1.010 (spread 1.5 %, 1 kept direction, H equal), warm 1.087 (spread 6.5 %); **VERDICT REJECT** (H columns did not move down, 824 -> 826; occupied eigenvalue moved 1.8e-6 Ha against the 1e-7 gate). The moved eigenvalues are lower in B (band energy -4.6e-5 Ha at r = 1): at this rough potential the plain restart stalls at higher Ritz values (the same 1.8e-6 Ha difference persists at warm tolerance 1e-12, where the baseline needs 1036 H columns and 7 contractions and r = 1 1075, +3.8 %), so the kept directions change which subspace the solve ends in, not how much work it does.

## 5. Stage 1 (O, `fast_screen_mpi.py`, production standard + verification + ladder, si32 np=8 simulated, one pair, m2_entry_cost sha256 7e5ffd710a69315d)
| arm B | SCF | Davidson | H columns | contractions | kept | rotation FLOP | sum m^3 | |dE| |
|---|---|---|---|---|---|---|---|---|
| A (lean algebra) | 11 | 67 | 1888 | 11 | | 3.12e9 | 4.21e8 | |
| r = 1 | 11 | 70 (+4.5 %) | 1877 (-0.6 %) | 11 | 11 | -1.2 % | +2.7 % | 2.0e-11 Ha |
| r = B/4 | 11 | 68 (+1.5 %) | 1899 (+0.6 %) | 13 | 206 | +8.4 % | +5.9 % | 7.5e-13 Ha |
Wall within noise in both directions (one pair; simulated ranks). Ceiling (C): 11 contractions in a whole si32 np=8 SCF of 22 solves, so even a perfect restart could only touch a small part of the 1888 H columns.

Threads path (serial lean, T=1, `dft-cli run`, capacity 6 per band): 0 contractions in si8/si16/si32 seed 601 (H 430 / 916 / 1717); the switch is a no-op (energy bits equal).

## 6. Capacity 2 per band: where kept directions do remove work (O, C; stage 1 as section 5, `M2_PROFILE_CAPACITY=2`, host load 24 to 29 during these runs, so walls are void, counters exact)
The idea (GD+k) pays when restarts are frequent relative to the iterations. At capacity 2 per band (the setting of MPI-SCALING-GAP-CAP2: sum m^3 -47 %, H +23 % on one warm solve, rejected) the whole si32 np=8 SCF has 51 contractions:

| arm | SCF | Davidson | H columns | contractions | kept | rotation FLOP | sum m^3 | tree sums | |dE| vs cap 2 |
|---|---|---|---|---|---|---|---|---|---|
| capacity 3, lean (current MPI rule) | 11 | 67 | 1888 | 11 | | 3.12e9 | 4.21e8 | 511 | |
| capacity 2, lean | 13 | 116 | 2720 | 51 | | 3.97e9 | 2.70e8 | 739 | |
| capacity 2, r = 1 | 13 | 115 | 2689 | 52 | 45 | 4.04e9 | 2.56e8 | 703 | 1.4e-11 Ha |
| capacity 2, r = B/4 | 11 | 91 | 2206 | 50 | 612 | 3.91e9 | 2.28e8 | 547 | 1.3e-10 Ha |
| capacity 2, r = B/2 | 11 | 91 | 2189 | 56 | 911 | 4.35e9 | 2.26e8 | 544 | 9.9e-11 Ha |
| capacity 3, r = B/2 | 11 | 65 | 1857 | 13 | 368 | 3.63e9 | 4.83e8 | 502 | 3.2e-11 Ha (vs cap 3) |

At capacity 2 the kept directions remove 62 percent of the H-column excess over capacity 3 (2720 -> 2206 at r = B/4, two SCF steps fewer: the plain capacity-2 solves are inexact enough to cost two extra SCF steps). Against the current rule (capacity 3) the pair "capacity 2 + r = B/4" still applies 17 percent more H columns and 25 percent more rotation FLOPs, for 46 percent less sum m^3 (the replicated reduced eigenproblem) and 7 percent more tree sums. One seed (601), one cell, simulated ranks.

## 7. Verdict and what would change it
REJECT at the current capacity rules (3 per band at np >= 4: H -0.6 to +0.6 %, rotation and GEMM up; 6 per band on the threads path: no contraction at all). Open lead, not a candidate yet: capacity 2 + r = B/4 against capacity 3 at np=8, where the replicated reduced eigenproblem is the term that does not scale. On the Mac counters it trades +17 % H columns for -46 % sum m^3; MPI-SCALING-GAP-CAP2 showed that on dev08 a +23 % H excess was not paid back by halving the eigenproblem, so the expected outcome is negative; one dev08 stage-1 job decides it (requests/breakthrough-restart.md). A capacity change is a default change under the generality rule (needs stage 2/3 on non-Si cells). The switch stays, default off.
