# dev08 stage decomposition, relayed by the master session (2026-11-08)
Relayed data, not produced by the orchestrator. Plan: docs/performance/stage-decomposition.md. Raw logs on dev08: /scr/mino/rustdft/stageA/.

## (A) QE 7.6 complete timer tables
Program PWSCF v.7.6, pw.x sha256 07c7cd7935f42b15157500ab719a5c3948825d8ce8a4487d992d7501c659658f (MKL DFTI build), `mpirun -np 1 pw.x -in <case>.in`, OMP/MKL threads 1, UPF Si.pz-vbc.UPF sha256 e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217, disk_io nowf, verbosity high, diagonalization david, diago_full_acc true, mixing_beta 0.7, FFT grid forced to BAANDI's. Gamma cases K_POINTS gamma (regterg/rdiaghg, nosym noinv), k-mesh cases symmetry on (cegterg/cdiaghg). n=3, load < 0.30. Input SHA-256: si8_1e-8 7850981919d50de372ba2e321c6dda8b64783aee58f986e03a0a888b0d4c70e9; si8_1e-10 6e8fc9c467825936968e785ee26e6a47847fc3734e92576aec79ef758ac9df30; si16_1e-8 cc5da991cba197d6e5231a95f5a05b221a6e33b483022f2c334a3ebd27c9488f; si16_1e-10 87450aac3661b7e2c0a910df77069b0dd1d098cd270f651087bc662030d15921; si32_1e-8 ea280746f1db62e8c117babc0901b969dc322aad16a74e89a148fae62b158b60; si32_1e-10 ba1b27028c6ee29604ed5676d9db9d86b1957372811736b22dc9e05808301545; k222s_1e-8 9da3a69c3a2b403892b1c42d9f5b75fabe91eeb927bbb88377acf2a00579b545; k222s_1e-10 312ee599a278836ad36b31124c34435c6c07acdec8c4382c43207baafd61387a; k444s_1e-8 00ee2b6e3cfecfb8cc6be38359973f1050c6f2aed52b4b37d1283c67db91e0f4; k444s_1e-10 249ee68a900c812bcfabb91b337ea1e977fb76a85583e8cf2f62c5fda5355ddf. The k444s 1e-8 input only changes conv_thr of the existing qe_k444s_sym.in; the k222s inputs were generated from BAANDI's sym_out_k222s results with make_qe_input.py --fft-grid ours and symmetry switched on (newly built inputs, not earlier QE references).

### Case metadata and work counters (QE side)
| case | wall (s), n=3 | SCF | k pts | Sym ops | FFT grid | G-vectors |
|---|---|---|---|---|---|---|
| si8 1e-8 | 0.4694 / 0.4726 / 0.4688 | 12 | 1 | none | 36^3 | 11936 |
| si8 1e-10 | 0.5518 / 0.5523 / 0.5504 | 15 | 1 | none | 36^3 | 11936 |
| si16 1e-8 | 1.4127 / 1.4110 / 1.4142 | 10 | 1 | none | 72x36x36 | 23926 |
| si16 1e-10 | 1.8275 / 1.9005 / 1.8313 | 13 | 1 | none | 72x36x36 | 23926 |
| si32 1e-8 | 5.0417 / 5.0283 / 5.0256 | 10 | 1 | none | 72x72x36 | 47906 |
| si32 1e-10 | 6.7000 / 6.7596 / 6.6973 | 13 | 1 | none | 72x72x36 | 47906 |
| k222s 1e-8 | 0.2554 / 0.2489 / 0.2476 (below 0.3 s: timer-resolution rule, not valid) | 6 | 2 | 48 | 25^3 | 5985 |
| k222s 1e-10 | 0.2587 / 0.2567 / 0.2554 (not valid) | 7 | 2 | 48 | 25^3 | 5985 |
| k444s 1e-8 | 0.5159 / 0.5143 / 0.5133 | 5 | 10 | 48 | 25^3 | 5985 |
| k444s 1e-10 | 0.6228 / 0.6268 / 0.6207 | 7 | 10 | 48 | 25^3 | 5985 |
`avg # of iterations` per SCF step: si8 3,2,4,10,2,7,2,16,12,2,17,4 (1e-10 adds 8,2); si16 3,2,2,4,16,20,9,20,20,20,20 (1e-10 adds 20,20,20); si32 2,2,1,2,15,8,20,20,20,13,8 (1e-10 adds 20,20,20); k222s 4,3,1.5,3,3,10.5,2.5; k444s 4.5,3.5,1.3,3.6,4.0,7.0. Many si16 and si32 steps sit at 20.0 (iteration cap).

### Timer table (WALL seconds / calls, run 1); columns: si8 1e-8 | si8 1e-10 | si16 1e-8 | si16 1e-10 | si32 1e-8 | si32 1e-10 | k222s 1e-8 | k222s 1e-10 | k444s 1e-8 | k444s 1e-10
- init_run: 0.04 | 0.04 | 0.09 | 0.09 | 0.26 | 0.26 | 0.03 | 0.03 | 0.05 | 0.05
- electrons: 0.31 | 0.39 | 1.21 | 1.62 | 4.66 | 6.31 | 0.10 | 0.11 | 0.35 | 0.46
- wfcinit: 0.01 | 0.01 | 0.05 | 0.05 | 0.20 | 0.20 | 0.00 | 0.00 | 0.02 | 0.02
- potinit: 0.01 | 0.01 | 0.01 | 0.01 | 0.02 | 0.02 | 0.01 | 0.01 | 0.01 | 0.01
- hinit0: 0.02 | 0.02 | 0.02 | 0.02 | 0.03 | 0.03 | 0.02 | 0.02 | 0.02 | 0.02
- c_bands (calls): 0.23 (12) | 0.30 (15) | 1.02 (11) | 1.38 (14) | 4.02 (11) | 5.48 (14) | 0.07 (7) | 0.08 (8) | 0.29 (6) | 0.38 (8)
- sum_band: 0.03 | 0.04 | 0.11 | 0.14 | 0.47 | 0.60 | 0.02 | 0.02 | 0.05 | 0.07 (sum_band:loo 0.03 | 0.03 | 0.10 | 0.13 | 0.44 | 0.56 | 0.01 | 0.01 | 0.04 | 0.06)
- v_of_rho: 0.03 | 0.03 | 0.04 | 0.06 | 0.09 | 0.12 | 0.00 | 0.01 | 0.00 | 0.01 (v_h 0.00 | 0.01 | 0.01 | 0.01 | 0.02 | 0.03; v_xc 0.02 | 0.03 | 0.04 | 0.05 | 0.07 | 0.09)
- mix_rho: 0.02 | 0.02 | 0.03 | 0.04 | 0.07 | 0.10 | 0 | 0 | 0 | 0
- regterg (Gamma): 0.23 | 0.29 | 1.00 | 1.35 | 3.93 | 5.37; cegterg (k): k222s 0.07 (14 calls) | 0.08 (16); k444s 0.29 (60) | 0.37 (80)
- rdiaghg (calls): 0.02 (93) | 0.03 (125) | 0.11 (146) | 0.16 (209) | 0.32 (121) | 0.47 (184); cdiaghg: 0.00 (67) | 0.01 (74) | 0.02 (289) | 0.03 (377)
- regterg:over / :upda / :last: si8 0.01/0.00/0.01 | 0.02/0.01/0.01; si16 0.06/0.04/0.03 | 0.09/0.05/0.04; si32 0.29/0.21/0.17 | 0.41/0.31/0.25 (calls over/upda 81, 110, 136, 196, 111, 171; last 34, 43, 32, 43, 27, 37); cegterg:over/:upda/:last k444s 0.01/0.00/0.00 | 0.01/0.01/0.01
- h_psi (calls): 0.18 (94) | 0.23 (126) | 0.76 (148) | 1.01 (211) | 2.95 (123) | 3.91 (186) | 0.06 (71) | 0.07 (78) | 0.25 (309) | 0.33 (397)
- h_psi:calbec / calbec: 0.01 | 0.01 | 0.06 | 0.09 | 0.26 | 0.36 | 0 | 0 | 0 | 0.01
- vloc_psi: 0.16 | 0.21 | 0.66 | 0.87 | 2.48 | 3.25 | 0.06 | 0.06 | 0.24 | 0.31
- add_vuspsi: 0.00 | 0.00 | 0.02 | 0.03 | 0.19 | 0.27 | 0 | 0 | 0 | 0
- g_psi: 0.00 | 0.00 | 0.01 | 0.01 | 0.03 | 0.04 | 0 | 0 | 0 | 0
- fft: 0.01 (38) | 0.01 (47) | 0.01 (33) | 0.02 (42) | 0.03 | 0.03; fftw (calls): 0.13 (964) | 0.17 (1208) | 0.55 (1612) | 0.72 (2120) | 2.11 (2768) | 2.75 (3612) | 0.06 (902) | 0.06 (988) | 0.24 (3812) | 0.31 (4932)
- PWSCF: 0.36 | 0.45 | 1.31 | 1.72 | 4.94 | 6.59 | 0.14 | 0.16 | 0.41 | 0.52
- newd, interpolate absent (norm-conserving PP; ecutrho = 4 ecutwfc so smooth and dense grids coincide); init_us_1, init_vloc and setup are not separate lines (inside init_run and hinit0).

### Master's S0 to S13 attribution for QE (1e-8, seconds)
| stage | QE timers | si8 | si16 | si32 |
|---|---|---|---|---|
| S0 | init_run - wfcinit - potinit | 0.020 | 0.000 | 0.040 |
| S1 | wfcinit + potinit | 0.020 | 0.060 | 0.220 |
| S2 | vloc_psi | 0.160 | 0.660 | 2.480 |
| S3 | h_psi:calbec + add_vuspsi | 0.010 | 0.080 | 0.450 |
| S4 | regterg:over + regterg:upda | 0.010 | 0.100 | 0.500 |
| S5 | rdiaghg | 0.020 | 0.110 | 0.320 |
| S6 | regterg:last | 0.010 | 0.030 | 0.170 |
| S7 | g_psi | 0.000 | 0.010 | 0.030 |
| S9 | sum_band | 0.030 | 0.110 | 0.470 |
| S10 | v_of_rho | 0.030 | 0.040 | 0.090 |
| S11 | mix_rho | 0.020 | 0.030 | 0.070 |
Rules: fftw is inside vloc_psi and sum_band, so it is a reference value and not added; regterg:over and :upda are S4, regterg:last is S6. QE timers have 0.01 s quantization: S0, S7, S10, S11 of si8 and si16 are not usable as ratios; si32 values (>= 0.03 s) are safe.

## (D) Ceilings, single core (relayed 2026-11-08)
Xeon Gold 6230 2.10 GHz, taskset -c 0, oneMKL 2026.1 sequential layer, best of 7, load < 0.30; benchmark crates outside the repository (/scr/mino/rustdft/{dgemmbench,fftbench}), RUSTFLAGS target-cpu=native. Shapes: real basis d = 1480 / 2958 / 5939 (si8 / si16 / si32), bands 20 / 36 / 68, subspace s = 6b = 96 / 216 / 408, projector channels n_proj = 4 per atom = 32 / 64 / 128.

D1 single-core MKL dgemm (GFLOP/s; square 2000^3 reference 98.47): si8 rotate 74.3, gram 76.9, reduced 100.6, nonlocal project 88.9, apply 90.0; si16 rotate 76.6, gram 87.7, reduced 83.9, project 80.5, apply 81.7; si32 rotate 83.8, gram 101.2, reduced 89.1, project 93.2, apply 82.2. All dominant Davidson and nonlocal shapes reach 75 to 103 percent of the square ceiling.

D2 one 3D FFT at production grids (microseconds): 25^3 MKL in place 76.09, Plan3d::transform 75.39, FFTW3 60.38; 36^3 147.81 / 193.98 (copy share 23.8%) / 128.58; 72x36x36 327.59 / 494.93 (33.8%) / 292.05; 72x72x36 727.19 / 1098.86 (33.8%) / 667.32. FFTW3 (GPL, reference only) is 8 to 21 percent faster than MKL.

D3 BAANDI pair transform (microseconds per band; library plan / mirror plan / single): 25^3 105.4 / 79.2 / 216.4; 36^3 239.1 / 127.8 / 497.2; 48^3 709.8 / 425.9 / 1365.7; 72x36x36 556.7 / 276.2 / 1064.4. 72x72x36 is not in the test list (added to `pair_microbenchmark` afterwards). Library-plan numbers vary 5 to 12 percent between sessions (cache interference); use ratios within one run.

D4 mirror vs ceilings (us per transform): 36^3 mirror 127.8, MKL 147.8, FFTW3 128.6, QE fftw timer 135 (0.130 s / 964 calls); 72x36x36 mirror 276.2, MKL 327.6, FFTW3 292.1, QE 341; 72x72x36 MKL 727.2, FFTW3 667.3, QE 762, mirror not measured. BAANDI's mirror pair transform is 14 to 16 percent faster than the MKL in-place call and 19 to 23 percent faster than QE's per-call FFT time; the FFT transform itself is no longer a cause of the gap (limits: different processes for MKL/FFTW3 and the pair microbenchmark).

## (D) addendum: si32 grid pair microbenchmark at ee86f07 (relayed 2026-11-09)
taskset -c 0, load 0.25, five grids in one run (us per band; library / mirror / single): 25^3 105.505 / 79.302 / 216.584; 36^3 238.139 / 127.305 / 492.675; 48^3 708.510 / 429.509 / 1370.457; 72x36x36 555.045 / 276.666 / 1065.582; 72x72x36 1165.697 / 647.459 / 2351.510 (mirror/library 0.555, pair/single 0.496). Completed table (us per transform): mirror vs MKL in place vs FFTW3 vs QE fftw timer: 36^3 127.31 / 147.81 / 128.58 / 135; 72x36x36 276.67 / 327.59 / 292.05 / 341; 72x72x36 647.46 / 727.19 / 667.32 / 762. The mirror pair transform is 11 to 16 percent faster than the MKL in-place call and 6 to 19 percent faster than QE at all three grids. Use ratios within one run; mirror/MKL and mirror/QE compare different processes.
