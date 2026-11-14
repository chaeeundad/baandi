# QE 7.6 and BAANDI, stage decomposition (dev08, f1e-3, strict-equivalent, one core)

Source: QE's own timers from the stored re-measurement outputs (/scr/mino/rustdft/qe-remeasure-2026-11-13, extraction script /scr/mino/rustdft/qestage.py, full output qestage_out.txt, relayed by the master session 2026-11-14); BAANDI: f3 run, main 38048b9 MKL, `timings` of results.json, median. The QE harness kept only the last of the 3 repetitions, so these QE numbers can differ from the 3-run medians in the headline table at the percent level. Gamma cells use QE's real-wavefunction Davidson (regterg); s_psi, newd, sum_bec and addusdens do not occur (norm-conserving).

## si32, conv_thr 1e-10 equivalent, f1e-3, np=1 (n=10) against BAANDI T=1
| QE stage | seconds | calls | BAANDI counterpart | BAANDI seconds |
|---|---|---|---|---|
| init_run | 0.24 | 1 | setup | 0.28 |
| c_bands | 6.26 | 15 | diagonalization | 13.61 |
| regterg (Davidson) | 6.14 | 15 | | |
| h_psi | 4.31 | 220 | h_apply | |
| vloc_psi (local FFT) | 3.51 | 220 | | |
| h_psi:calbec | 0.45 | 220 | nonlocal projection | |
| add_vuspsi | 0.32 | 220 | nonlocal back product | |
| rdiaghg | 0.61 | 218 | dense eigen | |
| regterg:overlap | 0.53 | 204 | orthonormalize | |
| regterg:update | 0.38 | 204 | reduced / update | |
| regterg:last | 0.27 | 40 | rotate | |
| sum_band | 0.64 | 15 | density | 0.74 |
| v_of_rho | 0.13 | 15 | potential | 0.23 |
| mix_rho | 0.11 | 15 | mixing | 0.10 |
| fftw | 2.96 | 3914 | (all FFT) | not separately measured |
| PWSCF | 7.37 | 1 | total | 19.84 |

## si64, same conditions (n=3)
| QE stage | seconds | calls | BAANDI counterpart | BAANDI seconds |
|---|---|---|---|---|
| init_run | 0.91 | 1 | setup | 0.85 |
| c_bands | 26.55 | 15 | diagonalization | 80.42 |
| h_psi | 17.22 | 186 | h_apply | |
| vloc_psi | 13.16 | 186 | | |
| calbec / add_vuspsi | 2.00 / 1.89 | 186 | nonlocal | |
| rdiaghg | 2.04 | 184 | dense eigen | |
| regterg:overlap | 2.51 | 170 | orthonormalize | |
| regterg:update | 2.20 | 170 | reduced / update | |
| regterg:last | 1.72 | 35 | rotate | |
| sum_band | 2.33 | 15 | density | 2.91 |
| v_of_rho | 0.27 | 15 | potential | 0.46 |
| mix_rho | 0.28 | 15 | mixing | 0.21 |
| fftw | 11.15 | 7114 | (all FFT) | not separately measured |
| PWSCF | 30.42 | 1 | total | 110.11 |

## BAANDI stage times, f1e-3, MKL 38048b9 (seconds)
| case | class | T | potential | density | mixing | energy | diagonalization | ladder coarse | setup | total |
|---|---|---|---|---|---|---|---|---|---|---|
| si32 | standard | 1 | 0.21 | 0.62 | 0.09 | 0.23 | 6.12 | 2.50 | 0.29 | 10.25 |
| si32 | standard | 8 | 0.20 | 0.20 | 0.09 | 0.23 | 2.79 | 1.57 | 0.35 | 5.59 |
| si32 | strict | 1 | 0.23 | 0.74 | 0.10 | 0.28 | 13.61 | 3.98 | 0.28 | 19.84 |
| si32 | strict | 8 | 0.22 | 0.24 | 0.10 | 0.27 | 5.52 | 2.40 | 0.34 | 9.42 |
| si64 | standard | 1 | 0.42 | 2.44 | 0.18 | 0.54 | 29.43 | 10.43 | 0.88 | 45.98 |
| si64 | standard | 8 | 0.42 | 1.07 | 0.18 | 0.58 | 13.77 | 5.99 | 0.96 | 23.65 |
| si64 | strict | 1 | 0.46 | 2.91 | 0.21 | 0.65 | 80.42 | 20.63 | 0.85 | 110.11 |
| si64 | strict | 8 | 0.46 | 1.43 | 0.21 | 0.70 | 32.45 | 10.58 | 0.93 | 48.70 |

## Stage ratio BAANDI / QE (f1e-3, strict-equivalent, T=1 against np=1)
| stage | si32 | si64 |
|---|---|---|
| potential | 1.77 | 1.70 |
| density | 1.16 | 1.25 |
| mixing | 0.91 | 0.75 |
| setup | 1.17 | 0.93 |
| diagonalization | 2.17 | 3.03 |
| total | 2.69 | 3.62 |
The gap is almost entirely diagonalization (and the ladder coarse stage, which QE has no counterpart for): at si64 BAANDI is 79.7 s slower, 53.9 s of it in diagonalization and 20.6 s in the ladder coarse stage; the other four stages differ by 1.5 s in total.

## Inside diagonalization (QE si64 strict f1e-3 np=1, c_bands 26.55 s)
h_psi 17.22 s (64.9 %), dense eigen 2.04, overlap 2.51, update 2.20, last 1.72, rest about 2.9 s. BAANDI MKL si64 strict (P48A, b1e-7, different amplitude, size comparison only): h_apply 28.0 s, rotate 20.6, orthonormalize 18.8, eigen 5.9, reduced 4.0. Ratios BAANDI / QE: h_apply 1.6x, rotate 12x (against regterg:last), orthonormalize 7.5x (regterg:overlap), eigen 2.9x (rdiaghg), reduced 1.8x (regterg:update). BAANDI's Davidson algebra (rotate + orthonormalize + eigen + reduced, 49.3 s) is 5.8x QE's 8.5 s. The Davidson iteration counts are similar (QE 170 overlap calls, BAANDI 176 Davidson iterations) but BAANDI computes 10068 H columns at much larger block width, which is an inference, not verified; the capacity study A13 (smaller basis, rejected) and the algebra cost model (V, about 5600N vs 1100N multiply-adds per step) are earlier evidence.

## Not comparable yet
QE's fftw timer separates the FFT (11.15 s, 37 percent of PWSCF at si64); BAANDI has no FFT timer ("included in enclosing phases"). A BAANDI FFT timer is needed for an FFT-only comparison.
