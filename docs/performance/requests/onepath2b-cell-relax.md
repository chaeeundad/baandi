# dev08 request: ONEPATH-2B variable-cell relaxation on the standard path after the metric fix

Status: pre-registered (2026-11-16), not yet run. Ticket `tickets/ONEPATH-2B-CELL-RELAX-GRID.md`.
BAANDI-only (no QE run or ratio; the parallel-model rule is satisfied trivially). Replaces the
failed Part C rows of `onepath2-stress.md` measured on binary 7a58d83.

## What is being checked
On 7a58d83 every standard-path variable-cell relaxation of real Si8 failed at the first step
that expanded the cell ("Miller index [0, 0, 18] is outside the canonical frequencies for
[36, 36, 36]"): the standard stop rule's Coulomb metric was built from the current cell's
density sphere instead of the frozen reference Miller set. Fixed in ONEPATH-2B. This list
checks that standard relaxations now converge to the strict reference's cell and positions.

## Build
Branch HEAD of ONEPATH-2B after merge with main. MKL build (`--features mkl`,
`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`). Binary sha256, build command
and exit code recorded. Same binary for every arm.

## Arms
- REF: `BAANDI_REFERENCE_STRICT=on`.
- S9: `[scf] standard_threshold_ha = 5e-9` written explicitly.
- S11: no threshold key (stress default 5e-11 Ha, `scf.threshold_source = "stress_default"`).
T=1 (taskset core 0). Node idle (load below 0.30), uptime before and after, arms interleaved.

## Cases (si8-gamma, `benchmarks/cases/si8-gamma.toml`, real PseudoDojo Si PBE, 30 Ry, FFT 36^3)
- C1: lattice x 0.98, `mode = "cell"`, `cell_dofs = "isotropic"`, `method = "bfgs"`,
  `force_tol_ha_bohr = 1e-3`, pressure 0 (the failing dev08 case).
- C1x: lattice x 1.02, same settings (expansion is the direction that failed; compression
  checks the other).
- C2: `mode = "ions_cell"`, `cell_dofs = "full"`, lattice x 1.01, fractional +-1e-2
  displacement (seed 601, MPI-M1-GATE3 generator), `force_tol_ha_bohr = 1e-4`,
  `stress_tol_ha_bohr3 = 1e-6`.
- C3: si-skew (`benchmarks/cases/si-skew.toml`), `mode = "ions_cell"`, `cell_dofs = "full"`,
  seed 601 +-1e-2, same tolerances.

## Report per run
Exit code, `relaxation.status`, steps, failures (step, attempt, reason), total SCF
iterations, total H vectors, verification attempts per frame, wall (s), final max force,
final max |sigma|, final enthalpy (17 digits), `rebuilt_pressure_check`, and max
|a(S) - a(REF)| of the final lattice columns and max |x(S) - x(REF)| of the final Cartesian
positions (Bohr).

## Pass/fail (pre-registered, not widened)
Every S9 and S11 run that REF converges converges too (exit 0, no failures entry with a
Miller-index reason), every frame `verification.passed = true`, final cell and positions
within 1e-3 Bohr of REF. A failure is a finding, reported with the run.

## Mac observations already made (not acceptance)
Real PseudoDojo Si PBE, Si8, 30 Ry, Mac release: before the fix the standard arm failed at
step 3 exactly as on dev08; after it C1-like (x 0.98, isotropic, force_tol 1e-3) converged
in 8 steps like the strict reference; numbers in ledger row ONEPATH-2B. Separate finding:
Si8 ions_cell with one displaced atom (lattice x 1.01, displacement 0.005 to 0.01) failed on
the Mac in BOTH arms, strict included (SCF divergence near step 20 to 25, same trajectory);
report C2/C3 REF failures as such, they are not regressions of the standard path.

## What to send back
Per run: arm, case, exit code, wall, relaxation summary as listed, `single_path` and `scf`
blocks, binary sha256, commit, build command, uptime before and after, interleave order.
Original result files stored on dev08 with sha256.
