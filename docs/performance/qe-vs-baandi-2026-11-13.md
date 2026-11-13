# QE 7.6 against BAANDI after the P48 round (si8 to si64, dev08, equal conditions)

Status 2026-11-14. This replaces the first 2026-11-13 table, whose QE input lacked `startingwfc = 'atomic'` (QE then runs a worse start and is too slow; see ledger QE-STARTINGWFC-CORRECTION). Measured by the master session on dev08.

QE 7.6, `K_POINTS gamma`, nosym/noinv, `startingwfc = 'atomic'`, `diagonalization = david`, `diago_full_acc`, same cell, cutoff and FFT grid as BAANDI, `nbnd` = BAANDI's band count (si32 68, si64 132), median of 3 runs, UPF SHA-256 51b1... (full hash in the raw log), 569 rows, raw data /scr/mino/rustdft/qe-remeasure-2026-11-13. BAANDI: main 3fea9e8 (MKL build, adoption-round defaults incl. S1 and the residual guard; P48B/C/D/F not in this build), P48A seeds 601-610. Same displaced coordinates for both codes (b1e-7 = +-1e-7 Bohr per Cartesian component, f1e-3 = +-0.001 fractional). Accuracy: BAANDI standard against QE `conv_thr` 1e-8, BAANDI strict against 1e-10.

## QE medians, wall seconds (nbnd = BAANDI value)
| cell | conv_thr | unperturbed np=1 | b1e-7 np=1 | f1e-3 np=1 | unperturbed np=8 | b1e-7 np=8 | f1e-3 np=8 |
|---|---|---|---|---|---|---|---|
| si8 | 1e-8 | 0.393 | 0.331 | 0.441 | 0.258 | 0.225 | 0.272 |
| si8 | 1e-10 | 0.468 | 0.414 | 0.531 | 0.270 | 0.244 | 0.306 |
| si16 | 1e-8 | 0.857 | 0.984 | 1.300 | 0.371 | 0.403 | 0.504 |
| si16 | 1e-10 | 1.227 | 1.415 | 1.669 | 0.484 | 0.525 | 0.631 |
| si32 | 1e-8 | 2.969 | 2.988 | 5.948 | 0.741 | 0.750 | 1.737 |
| si32 | 1e-10 | 3.484 | 3.512 | 7.413 | 0.833 | 0.860 | 2.138 |
| si64 | 1e-8 | 16.076 | 15.966 | 26.548 | 3.122 | 3.169 | 6.449 |
| si64 | 1e-10 | 21.510 | 21.896 | 30.467 | 4.115 | 4.564 | 7.497 |

## BAANDI divided by QE, geometry b1e-7 (BAANDI P48A medians, MKL; above 1 means BAANDI is slower)
| cell | class | T=1 / np=1 | T=8 / np=8 (threads against MPI) |
|---|---|---|---|
| si8 | standard | 0.459 / 0.331 = 1.39 | 0.449 / 0.225 = 2.00 |
| si8 | strict | 0.939 / 0.414 = 2.27 | 0.816 / 0.244 = 3.34 |
| si16 | standard | 1.964 / 0.984 = 2.00 | 1.471 / 0.403 = 3.65 |
| si16 | strict | 3.302 / 1.415 = 2.33 | 2.136 / 0.525 = 4.07 |
| si32 | standard | 8.168 / 2.988 = 2.73 | 4.768 / 0.750 = 6.36 |
| si32 | strict | 15.140 / 3.512 = 4.31 | 7.519 / 0.860 = 8.74 |
| si64 | standard | 38.935 / 15.966 = 2.44 | 20.283 / 3.169 = 6.40 |
| si64 | strict | 96.969 / 21.896 = 4.43 | 43.309 / 4.564 = 9.49 |

Unperturbed BAANDI (3fea9e8, same build) against the QE unperturbed medians, one core: si32 standard 6.32 / 2.969 = 2.13, strict 12.60 / 3.484 = 3.62; si64 standard 30.06 / 16.076 = 1.87, strict 96.78 / 21.510 = 4.50. Eight cores (threads against MPI): si32 standard 4.36 / 0.741 = 5.88, strict 7.00 / 0.833 = 8.40; si64 standard 17.86 / 3.122 = 5.72, strict 39.27 / 4.115 = 9.54.

## Threads against threads (the T=8 column above is threads against MPI)
QE si32 strict unperturbed: np=1 3.484 s, np=1 with OMP 8 threads 1.957 s (1.78x), np=8 MPI 0.833 s (4.18x). QE's OpenMP scaling is 43 percent of its MPI scaling, so BAANDI T=8 against QE np=8 compares BAANDI threads with QE's best mode. Thread against thread (si32 strict): 7.519 / 1.957 = 3.84 instead of 8.74.

## Displacement sensitivity (si32 strict, T=1 / np=1, seeds 601-605, increase over unperturbed; M2 partial)
| amplitude | BAANDI | QE | sensitivity ratio (BAANDI / QE) |
|---|---|---|---|
| b1e-7 | +23.3 % | +0.8 % | 29.1 |
| b1e-5 | +38.9 % | +43.6 % | 0.89 |
| b1e-3 | +49.8 % | +72.7 % | 0.68 |
| f1e-4 | +59.1 % | +102.7 % | 0.58 |
BAANDI is much more sensitive than QE only in the narrow regime where the symmetry has just been broken (b1e-7); at amplitudes closer to a real structure QE becomes the more sensitive code (QE si32 1e-10 f1e-3: +112.8 %). BAANDI's cost at +-1e-7 Bohr fell from 63 percent (before the adoption round) to 23.3 percent. The remaining M2 rows (f1e-3, f1e-2, standard, ladder on/off, T=8/np=8) are pending.

## Band-count effect on QE (si32, atomic start): QE default nbnd is 3.9 to 7.3 percent faster than nbnd = 68 for unperturbed and b1e-7 geometries. For displaced geometries it can reverse: si8 1e-10 f1e-3 with nbnd 20 needs 13.0 h_psi calls against 26.5 with the default nbnd 16 (wall 0.531 against 0.559 s), because QE without buffer bands runs Davidson twice as long. The fair comparison uses BAANDI's band count.

## Energy
BAANDI minus QE is +1.27e-6 Ha per atom at every size (si32 +4.07e-5 Ha, si64 +8.15e-5 Ha): the known discretisation offset, not a convergence error.

## Reading
1. One core, b1e-7: BAANDI is 1.4x (si8 standard) to 4.4x (si64 strict) slower. Strict is the larger gap and grows with size (2.3x at si8/si16, 4.3x at si32, 4.4x at si64).
2. Threads against MPI the gap is 2.0x to 9.5x; threads against threads (one cell measured) it is 3.84x for si32 strict.
3. Since this build the following were merged and are not yet in the table: P48B (pure GEMM, opt-in), P48C (auto), P48D (pure builds default on), P48F (phase threading, opt-in). Their dev08 A/B is queued.
