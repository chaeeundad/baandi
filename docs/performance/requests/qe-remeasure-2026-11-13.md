# dev08 request: QE 7.6 re-measurement at equal conditions (2026-11-13)

Why: the si32 QE standard wall of the 11-11 table (2.993 s) differs from the new equal-band-count run (5.00 s), and the QE perturbed columns exist only up to si32. One consistent QE data set is needed for the headline table.

Common: QE 7.6, `K_POINTS gamma`, nosym/noinv, best start (`startingwfc = atomic`), `diagonalization = david`, same cell, cutoff, FFT grid and UPF as BAANDI (record the UPF SHA-256), idle node (load below 0.30, record uptime before and after), np=1 on core 0, np=8 on socket 0 (`mpirun -np N`, OMP and MKL threads 1), 3 repetitions per cell, report the median and spread. Record command, version, input and band count for every run.

Matrix:
1. Cells: si8, si16, si32, si64 Gamma.
2. Classes: `conv_thr` 1e-8 (compare with BAANDI standard) and 1e-10 (strict).
3. Band count: BAANDI's count (si8 16+buffer as BAANDI reports in results `bands`, si16, si32 68, si64 132) AND QE's default `nbnd`, so the effect of nbnd is visible (si32 first, then the rest if it matters by more than 5 percent).
4. Geometry: unperturbed, and the displaced coordinates BAANDI uses: seeds 601-610 at +-1e-7 Bohr (Cartesian per component, `random.Random(seed)` generator of the MPI-M1-GATE3 protocol) and fractional +-0.001 (same stream). For si64 use seeds 601-603.
5. Report per run: wall (mpirun wall and PWSCF timer), SCF iterations, h_psi calls, final estimated accuracy, total energy (to 1e-9 Ry).
6. Also keep the threading comparison honest: QE with np=8 (MPI) is the reference; if cheap, add QE with OMP_NUM_THREADS=8 and np=1 for one cell (si32 strict) to see whether QE threads match its MPI scaling.

Deliver: one table per (cell, class, nbnd, geometry) with medians, energies, and a short statement of anything unexpected. Raw data stays under /scr/mino/rustdft/qe-remeasure-2026-11-13.
