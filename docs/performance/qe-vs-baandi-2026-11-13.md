# QE 7.6 against BAANDI after the P48 round (si32 and si64, dev08)

Status 2026-11-13. Measured by the master session on dev08. BAANDI: main 3fea9e8, `cli-main-mkl` 3e6a9ea476d77548 (MKL build, all adoption-round defaults including S1 and the standard residual guard; the P48B/C/D/F switches are not in this build). QE: 7.6, `K_POINTS gamma`, `orig-mpi-dfti`, same cell, cutoff, FFT grid and band count as BAANDI (si32 68 bands, si64 132 bands). Unperturbed coordinates for both codes, T=1 (BAANDI) against np=1 (QE), T=8 against np=8 (threads against MPI, mixed parallel models). Accuracy class: BAANDI standard against QE `conv_thr` 1e-8, BAANDI strict against `conv_thr` 1e-10.

## Wall seconds and ratio (BAANDI divided by QE; above 1 means BAANDI is slower)
| cell | class | BAANDI 1 core | QE np=1 | ratio | BAANDI 8 | QE np=8 | ratio |
|---|---|---|---|---|---|---|---|
| si32 | standard | 6.32 | 5.00 | 1.26x | 4.36 | 1.70 | 2.56x |
| si32 | strict | 12.60 | 6.59 | 1.91x | 7.00 | 1.97 | 3.55x |
| si64 | standard | 30.06 | 21.39 | 1.41x | 17.86 | 4.32 | 4.13x |
| si64 | strict | 96.78 | 27.68 | 3.50x | 39.27 | 5.81 | 6.76x |

Work counts, BAANDI (SCF / Davidson / H vectors): si32 standard 5/24/1119, si32 strict 8/54/2900, si64 standard 5/18/2086, si64 strict 8/147/6143. QE SCF iterations: si32 10 (1e-8), 13 (1e-10); si64 11 (1e-8), 15 (1e-10).

## Growth with size (log2 of the wall ratio, si32 to si64, one core)
| | QE | BAANDI |
|---|---|---|
| standard | 2.10 | 2.25 (perturbed-median table) / 2.25 here |
| strict | 2.07 | 2.68 (perturbed medians) / 2.94 here |
Standard has almost caught QE's growth rate (S1 and the residual guard); strict is the open gap. Thread scaling T=1 to T=8: QE 2.9 to 5.0 times, BAANDI 1.8 to 2.4 times; the 8-core ratio is the product of the single-core gap and the scaling gap (si64 strict 3.50 x 2.13 = 7.45 on the perturbed-median table).

## Energy
BAANDI minus QE is +4.07e-5 Ha (si32) and +8.15e-5 Ha (si64), +1.27e-6 Ha per atom at both sizes: the known discretisation offset, not a convergence error.

## Caveats
1. The si32 QE wall here (5.00 s at 1e-8) is larger than the 2.993 s of `qe-vs-baandi-2026-11-11.md`; the earlier runs may have used QE's default band count (64) against BAANDI's 68. The new numbers use equal band counts, so ratios in the two tables must not be mixed. To be confirmed by master (nbnd in the old input).
2. Perturbed ratios are not available for si64 on the QE side; the ratios above are unperturbed. BAANDI perturbed (+-1e-7 Bohr, seeds 601-610) medians at 3fea9e8 are slower: si32 standard 8.17 s, strict 15.14 s, si64 standard 38.94 s, strict 96.97 s (T=1); see ledger P48A-SCALING. The realistic-amplitude headline is pending the displacement scan (work order M2).
3. Since this build the following were merged but not yet measured against QE: P48B (pure GEMM, off), P48C (auto, helps k-point cells only), P48D (blocked orthonormalization, off), P48F (phase threading, off); dev08 A/B for all four is queued.

Sources: ledger rows QE-SI64-REFERENCE, QE-SI64-CORRECTION, P48A-SCALING, P48A-OTHER (master relays 2026-11-13); raw data on dev08 under /scr/mino/rustdft.
