# dev08 request: ONEPATH-2 stress and variable-cell relaxation on the single standard path

Status: pre-registered (2026-11-16), not yet run. Ticket `tickets/ONEPATH-2-STRESS.md`,
ADR 0043 amendment ONEPATH-2. BAANDI-only: no QE run and no QE ratio is requested, so the
parallel-model rule (`comparison-rules.md`) is satisfied trivially; every table is one model
at one thread count. Accuracy classes unchanged (strict 1e-10 Ha, standard 5e-9 Ha); stress
class 1e-7 Ha/Bohr^3 (synthetic stress oracle cap, `docs/validation/stress.md`).

## What is being checked
Stress (`output.stress`) and variable-cell relaxation (`relax.mode = "cell"`, `"ions_cell"`)
moved from the strict path to the standard path with the fixed-point verification; the
three cutoff-scan SCFs of the Pulay estimate run standard with verification at
`min(1e-9 Ha, scf.standard_threshold_ha)`. On the Mac the default 5e-9 Ha failed the
pre-registered stress class on one synthetic nonlocal fixture (1.305e-7 Ha/Bohr^3) and
5e-11 Ha met it on both (margin 1.6x). Orchestrator decision 2026-11-16: stress inputs
without an explicit threshold now default to 5e-11 Ha, PROVISIONAL until this list
confirms that 5e-11 meets the stress class on every geometry (if not, the default is
revisited, never the bound).

## Build
Branch HEAD of ONEPATH-2 after merge with main. MKL build (`--features mkl`,
`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`). Binary sha256, build command
and exit code recorded. Same binary for every arm.

## Arms
- REF: `BAANDI_REFERENCE_STRICT=on`, plain input (strict path with its certificate).
- TIGHT: REF with `[scf] energy_tolerance_ha = 1e-12`, `density_residual_tolerance = 1e-10`,
  `eigen_tolerance = 1e-11`, `max_iterations = 300` (stress truth for the error columns).
- S9: `[scf] standard_threshold_ha = 5e-9` written explicitly (stress inputs no
  longer default to 5e-9; `threshold_source = "user"`).
- S10: `[scf] standard_threshold_ha = 5e-10`.
- S11: plain input (the provisional stress default 5e-11 Ha,
  `threshold_source = "stress_default"`); this arm is the one that confirms or
  refutes the default.
Defaults otherwise. T=1 (taskset core 0) for every run; T=8 (cores 0-7) for the si16
arms S9 and S11 only (thread independence). Node idle (load below 0.30), uptime before and
after, arms interleaved per geometry.

## Part S: stress (`dft-cli run`, `output.stress = true`)
- Cells: si8-gamma, si16-gamma, si-skew, si-slab (`benchmarks/cases/*.toml`), each at the
  case geometry and at a fractional displacement +-1e-3 from the MPI-M1-GATE3 generator
  (`random.Random(seed)`, seeds 601-603), plus si8-gamma isotropically strained by +1 %
  (lattice vectors x 1.01) so the stress is far from zero.
- Arms REF, TIGHT, S9, S10, S11 on every geometry.
- Report per geometry and arm: max |sigma(arm) - sigma(TIGHT)| over components and the same
  against REF (Ha/Bohr^3), |dP|, |d Pulay correction|, max |sigma(TIGHT)| for scale, |dE| (Ha),
  SCF iterations, verification attempts, `pulay_assessment.scan_thresholds_ha`, wall of the
  run and of the Pulay scans (`pulay_assessment.elapsed_seconds`), `timings.verification_s`.
- Thread independence (si16, seeds 601-603, S9 and S11): stress tensor bits and density
  SHA-256 equal at T=1 and T=8.
- Pass/fail (pre-registered, not widened): every S run converges with `verified = true` and
  `stress.status = "available"` where REF converges; the stress class (1e-7 Ha/Bohr^3,
  against TIGHT) is evaluated per arm and reported as a count of geometries inside/outside,
  which is the evidence for the default-threshold decision.

## Part C: variable-cell relaxation (`dft-cli relax`)
- C1: si8-gamma compressed by 2 % (lattice x 0.98), `mode = "cell"`,
  `cell_dofs = "isotropic"`, `stress_tol_ha_bohr3 = 1e-6`, `max_steps = 60`.
- C2: si8-gamma, `mode = "ions_cell"`, `cell_dofs = "full"`, fractional +-1e-2 displacement
  (seed 601) and lattice x 1.01, `force_tol_ha_bohr = 1e-4`, `stress_tol_ha_bohr3 = 1e-6`.
- C3: si-skew, `mode = "ions_cell"`, `cell_dofs = "full"`, seed 601 +-1e-2, same tolerances.
- Arms REF, S9, S11 on every case (C1 also with `stress_tol_ha_bohr3 = 1e-7` as an
  information row).
- Report per run: exit code, `relaxation.status`, steps, total SCF iterations, total H
  vectors, verification attempts per frame, wall (s), final max force, final max |sigma|,
  final enthalpy (17 digits), rebuilt pressure check, and max |a(S) - a(REF)| of the final
  lattice columns and max |x(S) - x(REF)| of the final Cartesian positions (Bohr).
- Pass/fail: every S run that REF converges also converges (exit 0) with every frame
  `verification.passed = true`; final cell and positions within 1e-3 Bohr of REF (the Mac
  bound). A failure is a finding, reported with the run.

## Mac observations already made (not acceptance)
Synthetic fixtures, Mac release, ledger row ONEPATH-2: stress standard vs strict, cubic
two-H 2.1e-9 (5e-9) and 6.5e-10 (5e-11) Ha/Bohr^3; skewed nonlocal four-atom 1.305e-7
(5e-9, above the class) and 6.1e-8 (5e-11); default strict within 8.6e-11 of a tightened
strict reference. Variable-cell relaxations reach the strict cell within 1.6e-8 Bohr (BCC
two-H, cell isotropic) and 5.4e-6 Bohr (compact skewed two-H, ions_cell full).

## What to send back
Per run: arm, cell, seed, T, exit code, wall, `timings` (including `verification_s`),
`stress` block (values, preconditions, pulay_assessment), relaxation summary as listed,
`single_path` block, total energy (17 digits), density SHA-256 where written; binary
sha256, commit, build command, uptime before and after, interleave order. Original result
files stored on dev08 with sha256.
