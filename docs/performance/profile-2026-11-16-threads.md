# Stage profile of main 62c821f, one core (threads path, T=1), 2026-11-16

Labels: **O** = Mac observation (Apple Silicon, pure build without MKL, shared host, load 2 to 20 during the session); **R** = relayed dev08 number; **E** = estimate. Nothing here is acceptance evidence. The Mac pure build runs every GEMM on the pure-Rust kernel, so the GEMM share is higher than on the MKL build that dev08 compares with QE.

## Case and method
- Case: `benchmarks/cases/si32-gamma.toml` (32 Si, 30 Ry, 68 bands, 11877 plane waves, grid 72 x 72 x 36), fractional displacement +-1e-3 seed 601 (generator of `scripts/fast_screen.py`), single-path defaults (standard 5e-9 Ha, verification on, lean Davidson, cutoff ladder, REUSE-PSI on), T=1, release build of 62c821f, binary sha256 924af148...
- Stage seconds: `results.json` `timings` (closure, phases, ladder_coarse, verification_s) and `BAANDI_DAVIDSON_STAGES=on` stage timers (O).
- Function shares: macOS `sample` at 1 ms, one run, 3841 main-thread samples (inclusive counts, first occurrence on each stack; scratch parsers, not committed) (O).
- Work counters (exact, machine-independent): fine SCF 10 steps, coarse 7 steps; H vectors 1717 fine + 201 verification + 1018 coarse (coarse grid); Davidson iterations 53 + 34 coarse; real GEMM 64.9 GFLOP (gemm_dispatch rows); verification 1 attempt, 5 Davidson iterations.

## Wall breakdown (si32 seed 601, T=1, total 5.77 s O)
| stage | seconds (timers) | percent of wall | sample share | note |
|---|---|---|---|---|
| process + setup | 0.07 | 1.2 | 0.1 | process start not visible (CLI wall 5.79 vs total_s 5.77) |
| ladder coarse stage (7 SCF steps at 0.5 cutoff, includes its atomic start 0.17 s) | 1.15 | 19.8 | 18.6 | no QE counterpart |
| fine SCF steps without verification | 4.05 | 70.2 | 72.8 | |
| - Davidson solves (fine) | | | 55.8 | |
| - kept-start potential update (REUSE-PSI forward pair FFTs, all stages) | | | 6.1 | |
| - density build (fine) | 0.31 | 5.4 | 6.0 | |
| - nonlocal energy for the audit (`expectation_values`) | | | 2.1 | C3a cache removes it |
| - potential, energy, mixing, Coulomb metric | 0.13 | 2.2 | 3.4 | |
| verification (1 attempt) | 0.47 | 8.2 | 7.9 | 82 % of it is the Davidson solve |
| step overhead and remainder | 0.06 | 1.0 | | |

By kernel (all stages, sample shares, O):
| kernel | percent of wall |
|---|---|
| H apply (all) | 39.2 |
| - local part (pair FFTs and potential) | 29.6 |
| - nonlocal (P^T X and back product, GEMM) | 9.7 |
| Davidson algebra | 34.0 |
| - Ritz vectors and residuals of selected roots, return rotation | 12.7 |
| - reduced border V^T (H V_new) | 9.0 |
| - overlap border V^T V_new (lean basis is not orthonormal) | 8.8 |
| - reduced EVD (faer) and factor | 4.6 |
| density build | 8.3 |
| kept-start potential update (REUSE-PSI) | 6.1 |
| all FFT (rustfft) | 34.3 |
| all pure-Rust GEMM (`bounded_matmul`) | 42.0 |

Davidson stage timers (fine stage only, O): h_apply 1.55 s, rotate 0.69, orthonormalize (lean: overlap border) 0.39, reduced 0.38, eigen 0.15, precondition + residual 0.03; corrections by width: 15 iterations with 65+ columns, 16 with 9-64, 7 with 1-8, 10 with none.

## Verification and seeds (the per-run variance driver)
si32 f1e-3 seeds 601-610 (O, main): verification attempts 1, 1, **4**, 1, **3**, 1, 1, **3**, **2**, **5** (half of the seeds fail at least once); a failed attempt costs one verification solve (0.45 to 0.5 s, about 190 H vectors at the ratchet tolerance) plus one more fine SCF step (about 200 H vectors). Seed 603: 13 fine steps and 4 attempts, 7.8 to 9.0 s against 5.3 to 5.8 s for seed 601; verification 1.9 s (21 %). The stop estimator E_H[rho_out - rho_in] was 6.5e-10, 5.5e-10, 4.6e-10 Ha at the three failing stops while the verification measured E_H[F(rho_out) - rho_out] = 1.05e-8, 8.7e-9, 1.33e-8 Ha: on slowly contracting trajectories the stop rule fires 10 to 20 times too early, and every failure discards the verification map. Relayed (coordinator, 2026-11-16, 8f4d921 MKL): si64 standard T=1 verification 6.68 of 42.1 s (15.9 %), si32 0.56 of 8.3 s (R); master's mixer rows show attempts 1 to 5 per run under the default mixer (R, MIXER-BETA07-DEV08-CANDIDATE).

## si8 fixed overhead (seed 601, T=1, O; dev08 R: BAANDI 0.73 s vs QE 0.28 s)
Total 0.50 s (CLI wall equal within 5 ms, so process start is negligible on the Mac): setup 0.014 s (2.8 %), ladder coarse 0.058 s (11.6 %, 5 coarse steps, its atomic start 0.008 s), fine phases 0.41 s (11 fine steps; diagonalization 0.30 s), verification 0.10 s (20 %, **2 attempts**, 101 of 541 fine H vectors), step overhead 0.014 s. So at si8 the overhead that QE does not pay is the verification (20 %, one failed attempt here) and the coarse stage (12 %); setup and process start are small on the Mac. The 11 fine steps against QE's 7 to 9 (atomic start) are the other half.

## Against QE (stage decomposition 2026-11-14, R, strict-equivalent si32 one core)
QE si32: h_psi 58 % of PWSCF (vloc_psi FFT 48 %), regterg algebra (overlap, update, last, rdiaghg) 24 %, sum_band 9 %, init 3 %; QE needs no coarse stage and no separate verification and runs 7 to 10 SCF steps from the atomic start.
Biggest gap stages for the threads/one-core path, in order:
1. **Extra maps**: coarse stage (19.8 %) plus 10 fine steps plus 1 to 5 verification attempts against QE's 7 to 10 maps in total. A failed verification costs a solve and a step; on the seeds with failures this is 15 to 25 % of the run (O). Mixer beta (master's MIXER-BETA07 candidate) and VERIFY-REUSE (this session) act here.
2. **Verification cost when it passes** (8 % si32, 16 % si64 R): one warm-started solve at tau_v = ratchet(tau, E, N_e) (ADR 0043 V1), 5 Davidson iterations; the Ritz-change rule needs at least two iterations and a full correction block in the first. tau_v is fixed by the ADR (never looser than tau); reducing it is a user decision, not a performance change (see the dev08 request section "question for the user").
3. **Davidson algebra** 34 % (Mac pure; dev08 MKL is lower): the reduced and overlap borders are at QE's structure (QE also forms the H and S borders); the remaining excess is the basis size (capacity 6 per band; serial capacity 3 and 4 screened: REJECT, more verification attempts and Davidson iterations, ledger SERIAL-CAPACITY-SCREEN).
4. H apply local FFT: per-column cost near QE (backlog note H), count higher because of item 1.
