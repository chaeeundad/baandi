# QE fastest configuration per cell and model, and the 0.90 targets (dev08, 2026-11-15)

Source: master session, stored outputs of the P0.3 sweep (qeomp_out_omp_partial, qempi_out, qe64p03_out, qe32omp_out_part1/2/3); 280 usable rows (converged, `diago_full_acc = false`), 112 (cell, threshold, geometry, model, n) combinations, 0 missing. QE 7.6, `K_POINTS gamma`, `startingwfc = 'atomic'`, nbnd = BAANDI value, best of david and rmm-davidson per entry (ppcg does not exist in QE 7.6), f1e-3 = median over seeds of the per-seed minimum (si64 seed 601 only for some), unp = one unperturbed run. OpenMP: np=1 with OMP_NUM_THREADS = T on cores 0..T-1; MPI: `--map-by core --bind-to core`, one thread per rank, socket 0. Models are never mixed (AGENTS.md comparison rule).

## QE fastest wall seconds
| cell | conv_thr | geom | OMP T=1 | T=2 | T=4 | T=8 | MPI np=2 | np=4 | np=8 |
|---|---|---|---|---|---|---|---|---|---|
| si8 | 1e-8 | unp | 0.274 | 0.264 | 0.255 | 0.271 | 0.238 | 0.220 | 0.240 |
| si8 | 1e-8 | f1e-3 | 0.324 | 0.290 | 0.281 | 0.297 | 0.267 | 0.225 | 0.216 |
| si8 | 1e-10 | unp | 0.407 | 0.340 | 0.307 | 0.324 | 0.308 | 0.224 | 0.220 |
| si8 | 1e-10 | f1e-3 | 0.409 | 0.364 | 0.346 | 0.371 | 0.314 | 0.296 | 0.244 |
| si16 | 1e-8 | unp | 0.722 | 0.657 | 0.599 | 0.564 | 0.507 | 0.361 | 0.324 |
| si16 | 1e-8 | f1e-3 | 1.063 | 0.958 | 0.871 | 0.816 | 0.644 | 0.454 | 0.397 |
| si16 | 1e-10 | unp | 0.972 | 0.824 | 0.737 | 0.682 | 0.620 | 0.416 | 0.432 |
| si16 | 1e-10 | f1e-3 | 1.293 | 1.164 | 1.047 | 0.989 | 0.765 | 0.551 | 0.463 |
| si32 | 1e-8 | unp | 2.834 | 2.051 | 1.749 | 1.675 | 1.521 | 0.911 | 0.696 |
| si32 | 1e-8 | f1e-3 | 4.491 | 3.190 | 2.784 | 2.702 | 2.407 | 1.476 | 1.169 |
| si32 | 1e-10 | unp | 3.408 | 2.431 | 2.055 | 1.950 | 1.811 | 1.081 | 0.825 |
| si32 | 1e-10 | f1e-3 | 5.454 | 3.901 | 3.307 | 3.138 | 2.910 | 1.782 | 1.391 |
| si64 | 1e-8 | unp | 15.000 | 10.463 | 8.396 | 7.798 | 7.061 | 3.892 | 2.901 |
| si64 | 1e-8 | f1e-3 | 20.170 | 13.563 | 11.008 | 10.340 | 10.400 | 6.120 | 4.249 |
| si64 | 1e-10 | unp | 19.355 | 13.708 | 11.108 | 9.507 | 9.640 | 5.334 | 3.707 |
| si64 | 1e-10 | f1e-3 | 24.888 | 16.707 | 13.558 | 12.732 | 12.791 | 7.489 | 5.307 |
rmm-davidson wins only on si32 and si64 unperturbed or si32 1e-8 f1e-3; si8 and si16 are all david. OpenMP speed-up at T=8: si8 1.09x, si16 1.28x, si32 1.69x, si64 1.92x (saturates; si8 slower at 8 than at 4 threads). MPI speed-up np=2 to np=8: si8 1.14x, si16 1.57x, si32 2.19x, si64 2.60x. The displaced geometry is consistently slower (si64 1e-10 T=1: 19.355 -> 24.888 s, +28.6 %): the f1e-3 columns are the main comparison axis.

## Targets (0.90 x QE, f1e-3, seconds) and BAANDI today (main 38048b9 MKL, before the P48/P49/MPI work; BAANDI standard against QE 1e-8, strict against 1e-10)
One core (BAANDI T=1 against QE OMP T=1):
| cell | class | QE | target 0.90 x QE | BAANDI | BAANDI / QE |
|---|---|---|---|---|---|
| si8 | standard | 0.324 | 0.292 | 0.607 | 1.87x |
| si8 | strict | 0.409 | 0.368 | 0.917 | 2.24x |
| si16 | standard | 1.063 | 0.957 | 2.580 | 2.43x |
| si16 | strict | 1.293 | 1.164 | 4.305 | 3.33x |
| si32 | standard | 4.491 | 4.042 | 10.372 | 2.31x |
| si32 | strict | 5.454 | 4.909 | 19.955 | 3.66x |
| si64 | standard | 20.170 | 18.153 | 46.272 | 2.29x |
| si64 | strict | 24.888 | 22.399 | 110.388 | 4.44x |
Threads (BAANDI T=8 against QE OpenMP T=8):
| cell | class | QE | target | BAANDI | ratio |
|---|---|---|---|---|---|
| si8 | standard | 0.297 | 0.267 | 0.569 | 1.92x |
| si8 | strict | 0.371 | 0.334 | 0.788 | 2.12x |
| si16 | standard | 0.816 | 0.734 | 1.841 | 2.26x |
| si16 | strict | 0.989 | 0.890 | 2.672 | 2.70x |
| si32 | standard | 2.702 | 2.432 | 5.749 | 2.13x |
| si32 | strict | 3.138 | 2.824 | 9.586 | 3.05x |
| si64 | standard | 10.340 | 9.306 | 24.042 | 2.33x |
| si64 | strict | 12.732 | 11.459 | 49.086 | 3.86x |
MPI (BAANDI MPI np=8 before the 2026-11-15 entry-cost switches, against QE MPI np=8):
| cell | class | QE | target | BAANDI | ratio |
|---|---|---|---|---|---|
| si8 | standard | 0.216 | 0.194 | 0.626 | 2.90x |
| si8 | strict | 0.244 | 0.220 | 0.805 | 3.30x |
| si16 | standard | 0.397 | 0.357 | 1.973 | 4.97x |
| si16 | strict | 0.463 | 0.417 | 2.868 | 6.19x |
| si32 | standard | 1.169 | 1.052 | 7.699 | 6.59x |
| si32 | strict | 1.391 | 1.252 | 13.677 | 9.83x |
| si64 | standard | 4.249 | 3.824 | 35.608 | 8.38x |
| si64 | strict | 5.307 | 4.776 | 77.314 | 14.57x |
Notes: (1) The BAANDI columns do not yet include the merged P48B/C/D/F work (pure builds), PROJECT_WIDE and B1 (pending), the MPI entry-cost switches (si64 strict np=8 -22.6 %, np=2 -53.7 %) or the single path (verification, off by default). (2) "Strict" is compared with QE conv_thr 1e-10 only as the nearest threshold; under the single-path decision the comparison becomes standard threshold against QE conv_thr at each value. (3) The MPI axis is the furthest from the target.
