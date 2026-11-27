# R3 BAANDI observations on the Mac (2026-11-27; not oracle evidence)

Uncommitted worktree RPA-R3-ACFDT on main `33731e6d`, `cargo build -p
dft-cli` (dev profile, optimized), Apple Silicon shared with other jobs
(walls are observations only), `[scf] threads = 4`, one rank.
Si.pz-vbc.UPF SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`
(not committed). Base input `../si2-444.toml` (Si2 fcc a = 10.2612 Bohr,
LDA PZ, 30 Ry, FFT 25^3, SCF -7.919059559 Ha), variants differ only in the
named key (the base then had `ecut_chi_ha_series = [2.0, 2.5, 3.0, 4.0]`, `points = 16`; the committed input now carries the dev08 series). Command: `baandi rpa <input> --output-dir <name>`. Records:
`si2-444-<name>-rpa.json` (+ CSVs of `main`).

Cutoff series (`main`, k 4x4x4, crystal q symmetry: 8 irreducible q,
bands = N_G, 16 points, omega0 auto 3.902 Ha):

| E_cut^chi (Ha) | N_G = bands | E_c (meV/atom) |
| --- | --- | --- |
| 4 | 113 | -5782.99 |
| 5 | 137 | -5909.43 |
| 6 | 181 | -6012.57 |
| 8 | 283 | -6122.50 |

Fit `E_inf + A E^-1.5`: -6301.07 meV/atom, residuals -8.7 / 14.7 / 1.7 /
-7.7 meV/atom: the plan-2.5 linearity check FAILS (needs <= 0.5): cutoffs
4-8 Ha are below the asymptotic regime (shell structure of N_G). The 1/N_G
diagnostic fit gives -6353.8 (residuals up to 15 meV/atom). dev08 must run
higher cutoffs (e.g. 8, 10, 12, 16 Ha; bands up to ~1000) before any
extrapolated number is quoted. Max NSCF residual 1.5e-13 Ha; min Cholesky
pivot 1.0000002; Hermiticity defect 1.0e-16; q = 0 direction sets (Lebedev
6-50 vs 24x48 product): <= 2.5e-12 meV/atom (cubic: eps_M isotropic).

Frequency grid (E_cut^chi 4 Ha alone, omega0 2.073 Ha): 8 points
-5784.177, 16 points -5783.185, 32 points -5783.179 meV/atom: 16 -> 32
changes 0.006 meV/atom (gate 0.05: passed at 16).
Series dependence of the auto omega0: in `main` the largest transition
comes from 283 bands (omega0 3.902 Ha) and the 4 Ha point at 16 points is
-5782.988 (0.20 meV/atom off); `series32` ([4, 8] Ha, 32 points) gives
-5783.178 and -6122.690 (8 Ha moved by 0.19 meV/atom). Rule for series
runs: 32 points (or a convergence pair 16/32 recorded) whenever the series
spans a large band range; the per-system frequency check of plan 2.4 stays
mandatory.

Band count at E_cut^chi 4 Ha (fixed N_G 113): 113 bands -5783.19, 200
bands -5883.44, 400 bands -5894.00 meV/atom (monotone, converging).

q reduction at 3 Ha: crystal (8 q) / time reversal (36 q) / none (64 q):
-5397.5325 meV/atom each; differences to none -4.5e-11 and -5.4e-12 Ha
(SCF-density symmetry level). Walls 38 / 87 / 148 s.

k (= q) mesh at 4 Ha: 2x2x2 -7503.20, 3x3x3 -5943.88, 4x4x4 -5783.19
meV/atom: the q = 0 direction-average term carries the expected O(1/N_q)
error; 6x6x6 and 8x8x8 belong to dev08.

Absolute values are NC (valence-only) RPA correlation energies; they are
not comparable with all-electron or PAW literature numbers (plan 5).
