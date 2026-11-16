# dev08 request: ONEPATH-1 forces and fixed-cell relaxation on the single standard path

Status: pre-registered (2026-11-16), not yet run. Ticket `tickets/ONEPATH-1-FORCES-RELAX.md`,
ADR 0043 amendment ONEPATH-1. BAANDI-only: no QE run and no QE ratio is requested, so the
parallel-model rule (`comparison-rules.md`) is satisfied trivially; every table is one model
at one thread count. Accuracy classes unchanged (strict 1e-10 Ha, standard 5e-9 Ha).

## What is being checked
Forces and fixed-cell relaxation moved from the strict path to the standard path with the
fixed-point verification (forces need a passed verification; every ionic step verifies; a
failed verification continues that step's SCF). The user chooses the relaxation threshold
from the force differences against strict (work order P0.1, last bullet). Variable-cell
relaxation and stress stay strict and are not part of this list.

## Build
Branch HEAD of ONEPATH-1 after merge with main. MKL build (`--features mkl`,
`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`). Binary sha256, build command
and exit code recorded. Same binary for every arm.

## Arms
- REF: `BAANDI_REFERENCE_STRICT=on`, plain input (strict path with its certificate).
- S8: plain input (default standard path, `standard_threshold_ha = 5e-9`, verification on).
- S10: plain input with `[scf] standard_threshold_ha = 5e-11`.
Defaults otherwise (atomic start, gamma_real auto, cutoff ladder as on main for `run`; the
relaxation driver has no ladder). T=1 (taskset core 0) for every run; T=8 (cores 0-7) for the
si16 force arms only (thread independence). Node idle (load below 0.30), uptime before and
after, arms interleaved per geometry.

## Part F: forces (`dft-cli run`, `output.forces = true`)
- Cells: si8-gamma and si16-gamma (`benchmarks/cases/si8-gamma.toml`,
  `si16-gamma.toml`), fractional displacement +-1e-3 from the MPI-M1-GATE3 generator
  (`random.Random(seed)`, same stream as `qe-remeasure-2026-11-13.md`), seeds 601-610.
- Arms REF, S8, S10 on every geometry (60 runs).
- Report per geometry: max |F(S) - F(REF)| over components (Ha/Bohr), RMS over atoms, the
  largest REF component for scale and the ratio, |E(S) - E(REF)| (Ha), and the S-arm
  `forces.preconditions` block (threshold, attempts, Coulomb residual) and `verification`
  block (attempts, `timings.verification_s`). Medians and maxima over seeds per cell and arm.
- Thread independence (si16, seeds 601-603, S8): forces and density SHA-256 equal at T=1
  and T=8.
- Not pass/fail except: every S8 and S10 run converges with `verified = true` and
  `forces.status = "available"` where REF converges (a failure is a finding, reported with
  the run).

## Part R: fixed-cell relaxation (`dft-cli relax`)
- R1, displaced Si8: si8-gamma with fractional +-1e-2 displacement (seeds 601-603, same
  generator, amplitude chosen so that the relaxation takes several ionic steps), atom 0
  fixed, `[relax] force_tol_ha_bohr = 1e-4`, `max_steps = 60`, defaults otherwise.
- R2, synthetic slab: si-slab (`benchmarks/cases/si-slab.toml`, cutoff 8 Ry), fractional
  +-1e-2 displacement (seed 601), atom 0 fixed, same `[relax]` settings.
- R3, skewed cell: si-skew (`benchmarks/cases/si-skew.toml`), fractional +-1e-2 (seed 601),
  atom 0 fixed, same `[relax]` settings.
- Arms REF, S8, S10 on every geometry (also `start_extrapolation = "second_order"` with S8 on
  R1 seed 601 as an information row).
- Report per run: exit code, `relaxation.status`, ionic steps, total SCF iterations (sum of
  `frames[].scf_iterations`), total H vectors (sum of `frames[].scf_h_vectors`), verification
  attempts per step (S arms: `frames[].verification.attempts`), wall (s), final max force
  (Ha/Bohr), final energy (17 digits), and max |x(S) - x(REF)| of the final Cartesian
  positions (Bohr). Wall ratio S8/REF and S10/REF per geometry.
- Not pass/fail except: every S8 and S10 run that REF converges also converges (exit 0) and
  every S frame has `verification.passed = true`; a failure is a finding with the run.

## Mac observations already made (not acceptance)
Synthetic fixtures, Mac, see ledger row ONEPATH-1: forces max |dF| relative to the largest
component 5.8e-6 (cubic, 5e-9 Ha) and 3.3e-5 to 3.7e-5 (skewed four-atom nonlocal fixture,
5e-9 Ha); relaxation of a compact two-atom fixture: final positions 1.6e-6 Bohr apart, SCF
iterations 41 (standard) against 122 (strict).

## What to send back
Per run: arm, cell, seed, T, exit code, wall, `timings` (including `verification_s`),
`convergence` summary or relaxation summary as listed, `single_path` block,
`forces.preconditions`, total energy (17 digits), density SHA-256 where written; binary
sha256, commit, build command, uptime before and after, interleave order. Original result
files stored on dev08 with sha256.
