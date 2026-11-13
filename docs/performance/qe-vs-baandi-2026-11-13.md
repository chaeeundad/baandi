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

## Displacement sensitivity, M2 (si32 Gamma, BAANDI T=1 / QE np=1, same geometries, seeds 601-605, wall seconds and increase over unperturbed)
| class | amplitude | BAANDI | increase | SCF / H | QE | increase | it / h_psi | sensitivity ratio (BAANDI / QE) |
|---|---|---|---|---|---|---|---|---|
| strict | unperturbed | 12.326 | - | 8 / 2900 | 3.444 | - | 8 / 35 | - |
| strict | b1e-7 | 15.197 | +23.3 % | 9 / 3556 | 3.470 | +0.8 % | 8 / 35 | 29.1 |
| strict | b1e-5 | 17.119 | +38.9 % | 10 / 4082 | 4.946 | +43.6 % | 10 / 109 | 0.89 |
| strict | b1e-3 | 18.467 | +49.8 % | 10 / 4467 | 5.948 | +72.7 % | 11 / 163 | 0.68 |
| strict | f1e-4 | 19.608 | +59.1 % | 11 / 4795 | 6.981 | +102.7 % | 13 / 193 | 0.58 |
| strict | f1e-3 | 19.935 | +61.7 % | 11 / 4656 | 7.347 | +113.3 % | 14 / 222 | 0.54 |
| strict | f1e-2 | 26.103 | +111.8 % | 15 / 6416 | 8.142 | +136.4 % | 17 / 221 | 0.82 |
| standard | unperturbed | 6.270 | - | 5 / 1119 | 2.937 | - | 7 / 28 | - |
| standard | b1e-7 | 8.178 | +30.4 % | 7 / 1704 | 2.953 | +0.6 % | 7 / 28 | 50.7 |
| standard | b1e-5 | 8.938 | +42.5 % | 9 / 1954 | 3.764 | +28.2 % | 8 / 67 | 1.51 |
| standard | b1e-3 | 9.074 | +44.7 % | 9 / 1963 | 5.222 | +77.8 % | 10 / 135 | 0.57 |
| standard | f1e-4 | 9.213 | +46.9 % | 9 / 1986 | 5.436 | +85.1 % | 10 / 142 | 0.55 |
| standard | f1e-3 | 10.443 | +66.6 % | 10 / 2305 | 6.124 | +108.6 % | 11 / 180 | 0.61 |
| standard | f1e-2 | 11.875 | +89.4 % | 11 / 2449 | 6.346 | +116.1 % | 13 / 160 | 0.77 |
No nonconvergence; BAANDI standard error against the strict reference of the same coordinates is at most 3.879e-9 Ha at every amplitude (limit 5e-9). BAANDI is far more sensitive than QE only in the narrow regime where the symmetry has just been broken (b1e-7, b1e-5 for standard); from b1e-3 on QE is 1.2 to 1.9 times more sensitive. BAANDI's cost at +-1e-7 Bohr fell from 63 percent (before the adoption round) to 23.3 percent. f1e-2 (about 0.2 Bohr in this cell) is a physically large displacement; f1e-4 and f1e-3 are the representative amplitudes.

Ladder effect (si32 standard T=1, ladder on against off, each arm against its own unperturbed): the ladder lowers the unperturbed baseline strongly (6.29 against 7.30 s, SCF 5 against 10), so its percentage increase is larger at every amplitude (b1e-7 +30.1 against +13.6 %, f1e-3 +65.6 against +56.2 %), yet the absolute wall is lower with the ladder at every amplitude except b1e-5 (8.933 against 8.707 s). Quote absolute times, not only increases. The growth of BAANDI's wall with amplitude is work growth (si32 standard SCF 5 to 9, Davidson iterations 24 to 59, H columns 1119 to 1954 at b1e-5), mostly narrower correction blocks (columns per Davidson iteration 46.6 to 33.1); the cost per H column does not rise (-6 % standard, +9.5 % strict).

Eight cores (M2c, si32, BAANDI T=8 against QE np=8): strict unperturbed 6.334 / 0.845 s; b1e-7 +19.1 % / +1.7 % (sensitivity ratio 11.2); f1e-3 +51.4 % / +152.6 % (ratio 0.34). Standard unperturbed 3.968 / 0.740 s; b1e-7 +20.0 % / +1.1 % (ratio 18.1); f1e-3 +44.3 % / +141.9 % (ratio 0.31). QE's MPI absorbs the extra work of displaced geometries worse (np=1 to np=8 speed-up 3.44 at f1e-3 against 4.08 unperturbed) while BAANDI's thread speed-up barely changes (2.08 against 1.95). Ladder at f1e-2 (si32 standard T=1): on 11.894 s (SCF 11, H 2449) against off 15.233 s (SCF 22, H 4393): the ladder gains most on the hardest geometry (-21.9 %). Threads against MPI, f1e-3, strict: 9.589 / 2.134 = 4.49x; standard 5.726 / 1.792 = 3.20x.

## Headline (decision of 2026-11-11, realistic amplitude f1e-3), si32 Gamma, one core
| class | BAANDI | QE (atomic start, nbnd 68, conv_thr 1e-8 / 1e-10) | BAANDI / QE |
|---|---|---|---|
| standard | 10.443 | 6.124 | 1.71x |
| strict | 19.935 | 7.347 | 2.71x |
Other cells need BAANDI at f1e-3 (queued); the QE f1e-3 medians are in the table above.

## Band-count effect on QE (si32, atomic start): QE default nbnd is 3.9 to 7.3 percent faster than nbnd = 68 for unperturbed and b1e-7 geometries. For displaced geometries it can reverse: si8 1e-10 f1e-3 with nbnd 20 needs 13.0 h_psi calls against 26.5 with the default nbnd 16 (wall 0.531 against 0.559 s), because QE without buffer bands runs Davidson twice as long. The fair comparison uses BAANDI's band count.

## Energy
BAANDI minus QE is +1.27e-6 Ha per atom at every size (si32 +4.07e-5 Ha, si64 +8.15e-5 Ha): the known discretisation offset, not a convergence error.

## Reading
1. One core, b1e-7: BAANDI is 1.4x (si8 standard) to 4.4x (si64 strict) slower. Strict is the larger gap and grows with size (2.3x at si8/si16, 4.3x at si32, 4.4x at si64).
2. Threads against MPI the gap is 2.0x to 9.5x; threads against threads (one cell measured) it is 3.84x for si32 strict.
3. Since this build the following were merged and are not yet in the table: P48B (pure GEMM, opt-in), P48C (auto), P48D (pure builds default on), P48F (phase threading, opt-in). Their dev08 A/B is queued.
