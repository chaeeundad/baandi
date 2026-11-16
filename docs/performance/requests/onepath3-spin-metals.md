# dev08 request: ONEPATH-3 collinear spin and Fermi-Dirac on the single standard path

Status: pre-registered (2026-11-16), not yet run. Ticket `tickets/ONEPATH-3-SPIN-METALS.md`,
design `docs/design/onepath3-spin-metals.md`, ADR 0043 amendment ONEPATH-3. BAANDI-only:
no QE run and no QE ratio is requested, so the parallel-model rule
(`comparison-rules.md`) is satisfied trivially; every table is one model at one thread
count. Accuracy classes unchanged (strict 1e-10 Ha, standard 5e-9 Ha).

## What is being checked
Collinear spin (fixed sector and Fermi-Dirac) and Fermi-Dirac occupations moved from the
strict drivers to the standard path with the fixed-point verification. Collinear stop
estimate E_H[d_n] + E_H[d_m] + E_0[d_m]; residual guard with the actual occupations; the
verification repeats the whole map including the Fermi-level solve. The magnetization term
is an estimator without an energy-bound proof, so delivered accuracy on real
pseudopotentials and larger cells is the question. Mac evidence (synthetic, ledger
ONEPATH-3): every pre-registered bound met with margins of 20x (energy) and 15x
(first-order quantities).

## Build
Branch HEAD of ONEPATH-3 after merge with main. MKL build (`--features mkl`,
`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`). Binary sha256, build command
and exit code recorded. Same binary for every arm.

## Arms
- REF: `BAANDI_REFERENCE_STRICT=on`, plain input (strict drivers with their certificates).
- TIGHT: REF with `[scf] energy_tolerance_ha = 1e-12`, `density_residual_tolerance = 1e-10`,
  `eigen_tolerance = 1e-11`, `max_iterations = 400` (truth for the error columns).
- S9: plain input (default 5e-9 Ha, `verify = "fixed_point"`).
- S11: `[scf] standard_threshold_ha = 5e-11`.
T=1 (taskset core 0) for every run; T=8 (cores 0-7) for the S9 arm of F2 and M2 only
(thread independence). Node idle (load below 0.30), uptime before and after, arms
interleaved per case.

## Fixtures
Real NC files from `DFT_NC_DIR` (PseudoDojo NC, SHA-256 recorded; never committed).
- F1 ferromagnetic fixed sector: bcc Fe-like 2-atom cell with the Fe NC file
  (`[spin] mode = "collinear"`, `total_moment` set to the strict REF's integer-nearest
  sector that the REF certifies; if REF reports NonInsulating, record and move to F2).
- F2 ferromagnetic FD, free moment: the same Fe cell, `[occupations] mode = "fermi_dirac"`,
  `smearing_ha = 0.005`, `band_cap` from the REF band study, `initial_moment_seed = 4.0`,
  k mesh 4x4x4 (`ktime_reversal = false`).
- F3 fixed-moment FD: F2 with `total_moment = 2.0`.
- M1 nonmagnetic metal, scalar FD: fcc Al primitive cell (Al NC), `smearing_ha = 0.01`,
  k mesh 6x6x6 shifted.
- M2 larger metal: Al 4-atom conventional cell, 4x4x4 mesh, `smearing_ha = 0.01`.
- U1 unpolarized limit: si8-gamma (`benchmarks/cases/si8-gamma.toml`) with
  `[spin] mode = "collinear"`, `total_moment = 0.0`, `xc.functional = "lda-pw"` (and an
  LDA-PW Si file), against the same input without `[spin]`, both S9 and both with
  `gamma_real = "off"` (spin runs are complex).
- G1 generality: si-skew and si-slab cells with one atom replaced by Al (odd electron count),
  scalar FD `smearing_ha = 0.01` and collinear FD free moment.

## Report per run
exit code, wall (s), `timings.verification_s`, SCF iterations, verification attempts and
the last attempt's `coulomb_error_ha`, `residual_bound_ha`, `chemical_potential_change_ha`,
`moment_change_mu_b`; energy E or free energy F (17 digits), `spin.moment_mu_b`, every
chemical potential, the occupation block; `single_path` and `execution` blocks; density
SHA-256.

## Accuracy table (pre-registered bounds, not widened)
Per fixture: |F(S) - F(TIGHT)| and |F(REF) - F(TIGHT)| (Ha), |M(S) - M(TIGHT)| (mu_B),
max |mu(S) - mu(TIGHT)| (Ha), sum_k w_k |f(S) - f(TIGHT)| (electrons), SCF iterations and
wall of S9, S11, REF.
- Pass: S9 converges with `verified = true` wherever REF converges; |dF| <= 5e-9 Ha
  (S9) and <= 5e-11 Ha + |REF - TIGHT| (S11, information); |dM| <= 1e-4 mu_B,
  |dmu| <= 1e-4 Ha, sum w|df| <= 1e-4 electrons (S9). U1: |E(spin) - E(no spin)| <=
  1e-10 Ha and SCF counts within 1.
- Thread independence (F2, M2, S9): energy bits and density SHA-256 equal at T=1 and T=8.
- A failure is a finding: report it with the run; the bound is not widened, the default
  threshold for spin/FD inputs is revisited instead (as for stress in ONEPATH-2).

## Mac observations already made (not acceptance)
Synthetic fixtures, Mac release, ledger ONEPATH-3: zero-moment spin vs closed shell 1e-15 Ha
(same SCF count); ferro fixed sector |dE| 2.3e-10 Ha; scalar FD metal |dF| 2.0e-12 Ha,
|dmu| 8.9e-7 Ha, sum w|df| 3.8e-6 e; collinear FD free moment |dF| 5.5e-11 Ha, |dM|
6.7e-6 mu_B, |dmu| 6.6e-6 Ha; SCF 15/5/12 standard against 45/9/39 strict.

## What to send back
Per run: arm, fixture, T, exit code, wall, `timings`, `convergence` block, `verification`,
`single_path`, `execution`, `occupations`, `spin`, `thermodynamics`, density SHA-256; NC
file SHA-256s, binary sha256, commit, build command, uptime before and after, interleave
order. Original result files stored on dev08 with sha256.
