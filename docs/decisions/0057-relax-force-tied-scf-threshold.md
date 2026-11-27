# ADR 0057: Force-tied default SCF threshold for relaxations

- Status: **Proposed** (ticket RELAX-FORCE-NOISE, 2026-11-28; orchestrator
  review pending; dev08 confirmation not yet requested).
- Scope: `crates/dft-cli/src/config.rs` (default threshold resolution),
  `crates/dft-cli/src/relax.rs` (coefficient and helper), `run.rs`
  (reference profile uses the same rule). No change to SCF, forces or
  optimizer numerics. Inputs without `[relax]` and inputs that set
  `scf.standard_threshold_ha` are unchanged bit for bit.
- Relates to: ADR 0043 (single standard path; ONEPATH-1 item 3 left the
  relaxation threshold to the user, ONEPATH-2 item 5 introduced the
  observable-based stress default), ADR 0056 (default force tolerance
  0.02 eV/A = 3.889e-4 Ha/Bohr on the atomic norm).

## Finding

At the default threshold 5e-9 Ha, warm-started ionic SCFs delivered forces with
errors up to 2.3e-4 Ha/Bohr on MgO PBE (Gamma, 120 Ry) and 2.8e-4 Ha/Bohr on
wurtzite ZnO PBE (k 2x2x2, 80 Ry), i.e. 0.6 to 0.7 of the default force
tolerance. A relaxation could therefore stop at a noise-influenced point, and
BFGS/FIRE cycled at a noise floor (MgO with `force_tol_ha_bohr = 1e-5`: no
convergence in 80 steps).

## Measurement (Mac, release, 2026-11-28)

Protocol: a default `baandi relax` per system; for every accepted frame a cold
`baandi run` with forces at the same geometry at 5e-12 Ha (reference) and at
5e-9 Ha; error = max over free atoms of |F - F_ref| (Ha/Bohr). The 5e-12 Ha
reference agrees with the strict reference (`BAANDI_REFERENCE_STRICT=on`) to
about 5e-7 Ha/Bohr on MgO. Scripts and raw data:
`/Volumes/990-2TB/Develop/rust-dft-scratch/relax-force-noise` (study.py,
cfit.py, runs/).

| system | warm frames, 5e-9 (max) | cold, 5e-9 (max) | warm frames, force-tied (max) | evals / SCF iter / H vectors, 5e-9 -> force-tied |
|---|---|---|---|---|
| Si8 LDA k222, PRELAX seed 1, 0.2 Bohr | 2.65e-5 | 1.80e-5 | 1.93e-6 | 18/119/44932 -> 17/156/54223 |
| Si8-AlP PBE k222, seed 1 | 4.31e-5 | 2.00e-5 | 2.80e-6 | 16/129/52082 -> 16/171/61104 |
| NaCl 8 atoms PBE Gamma 70 Ry, 0.15 Bohr | 5.88e-5 | 1.88e-5 | 6.18e-6 | 8/64/5067 -> 8/76/5533 |
| ZnO wurtzite PBE k222 80 Ry, 0.1 Bohr | 2.80e-4 | 7.67e-5 | 2.62e-5 | 14/111/51182 -> 13/131/57409 |
| MgO PBE Gamma 120 Ry (real_si_and_mgo geometry), last 5 frames | 1.25e-4 | 2.45e-5 | 8.69e-6 | 55/368/8318 -> 55/443/9473 |
| MgO, force_tol 1e-5 probe, frames near the minimum | 2.33e-4 (29 frames, median 7.9e-5) | 3.66e-5 | 1.15e-5 | 81 (not converged, final 6.5e-5) -> 58 (converged, final 3.4e-6) |

The force-tied threshold here is 4.727e-11 Ha (default tolerance). Cost:
SCF iterations +18 to +33 percent, H vectors +9 to +21 percent per relaxation
on these cases; ionic evaluations equal or fewer.

## Root cause

1. The force error of a verified standard state is first order in the
   residual of its delivered density: |dF| ~ C sqrt(2 eps), eps the
   verification's Coulomb residual E_H[F(y) - y] (Ha). Measured C (warm
   frames, max): Si8 0.29, Si8-AlP 0.45, NaCl 1.7, MgO 3.9, ZnO 5.5
   Ha^1/2/Bohr. The Coulomb metric (4 pi / G^2) down-weights the short-
   wavelength density error around hard ions (O 2p, semicore), which is what
   forces see; C therefore varies by 20x between systems.
2. Cold SCFs converge geometrically and stop with residuals 10 to 100 times
   below the threshold (MgO frame 60: stop dr2 7e-10 Ry, verification
   2.5e-11 Ha). A warm SCF starts with a residual already near the threshold
   (the previous SCF's leftover error plus the small move), stops after one
   or two maps with no overshoot, and delivers eps close to the threshold
   (MgO frame 60: one map, stop dr2 9.1e-9 Ry, verification 3.98e-9 Ha against
   5e-9). The estimate is not inherited: the warm first map is redone at
   ethr = 0.1 dr2/N (QE-style first-step redo) and dr2 is measured at the new
   geometry; it is correct, the threshold simply allows that much error.
   Warm energies are as accurate as cold ones (NaCl |dE| <= 1.5e-9 Ha); only
   the first-order force sees the difference.
3. Not the cause: the iteration count (MgO frames with 2 to 7 iterations had
   errors up to 1.7e-4, so a minimum iteration count does not help), the
   start extrapolation (`first_order`: MgO probe warm max 1.9e-4, median
   7.2e-5, unchanged), or a stale warm-start residual estimate.

## Decision

A relaxation that moves atoms (`relax.mode = "ions"` or `"ions_cell"`) and does
not set `scf.standard_threshold_ha` runs at

    threshold_ha = max(FLOOR, min(base, force_tol^2 / 3200))

with `force_tol` = `relax.force_tol_ha_bohr` (Ha/Bohr), base = the class default
5e-9 Ha (ions) or the stress default 5e-11 Ha (ions_cell), FLOOR = 5e-14 Ha.
The coefficient 1/3200 Bohr^2/Ha = (0.1 / C_ref)^2 / 2 with C_ref = 4: the
worst-case force error C_ref sqrt(2 threshold) is 0.1 force_tol. At the default
tolerance the threshold is 4.73e-11 Ha; 1e-3 Ha/Bohr gives 3.1e-10; 1e-4 gives
3.1e-12; tolerances above 4e-3 keep 5e-9. ZnO's warm C (5.5) exceeds C_ref,
but its measured errors under the rule stay at 0.067 force_tol (margin 1.5x);
the constant is a rule by observable (force tolerance), not by material, and
is set by the hardest case of the set.

- An explicit `scf.standard_threshold_ha` always wins (looser or tighter).
- `results.json` records `single_path.threshold_source = "relax_force_default"`
  (also in the relaxation `scf` block) when the rule set the threshold.
- `relax.mode = "cell"` keeps the stress default.
- The trust-radius energy-rise tolerance (twice the threshold when
  `energy_tol_ha` is unset, ADR 0056) follows the tighter threshold.
- Like the stress default, the rule follows the input: a `run` of an input
  that contains a `[relax]` table resolves the same threshold. Inputs without
  `[relax]` are bit-identical (checked: MgO cold run, energy and forces equal).

## Alternatives considered

- Minimum SCF iterations after a geometry change: rejected by the data
  (2-iteration stops 1.2e-4 to 1.6e-4 on MgO).
- Fix the warm residual estimate: there is no stale estimate (item 2).
- QE-style progressive threshold (tight only near convergence, documented in
  QE's `upscale` input) would recover most of the +9 to 21 percent H-vector
  cost; needs a per-step threshold and a confirmation SCF before a stop.
  Follow-up candidate.
- SCF correction force term (Chan, Bohnen and Ho, PRB 47, 4771 (1993);
  described in `docs/research/qe-study/qe-stage-by-stage.md` section 7): makes
  the force error second order in the density error for every run, not only
  relaxations. A forces-numerics change with its own ADR and validation;
  follow-up candidate, after which the coefficient could be loosened.

## Open

dev08 confirmation on the generality set (cost and errors), and the two
follow-ups above.
