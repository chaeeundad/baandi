# ADR 0056: Ionic stop rule on the atomic force norm, optional energy test, FIRE in ASE units

- Status: **Proposed** (ticket RELAX-STOP-RULE, 2026-11-27; user decision
  option A of 2026-11-27 on the PRELAX-P0 audit; orchestrator review pending).
- Scope: `crates/dft-cli/src/relax.rs`, `relax/optimizer.rs`, `relax/cell.rs`
  (force measure and logs only). No change to SCF, forces or stress numerics.
- Evidence: [precon-relax-audit.md](../roadmap/precon-relax-audit.md)
  sections 2.2, 2.4, 3 and 5 (Mac), PRELAX-P1 stages 1-2 on dev08 (orchestrator
  record): native BFGS reached 0.02 eV/A (maximum atomic force norm) at
  evaluation 12/17/16 (medians) but stopped only after 17/22/22-40; native FIRE
  needed 3.6-4.2x the evaluations of native BFGS.

## Context

Before this ADR the fixed-cell stop rule was: the largest absolute free
Cartesian force COMPONENT below `force_tol_ha_bohr` AND |E_k - E_{k-1}| below
`energy_tol_ha` (default 1e-8 Ha, twice the default SCF threshold), both on two
consecutive accepted steps. The component measure is looser than the atomic
norm used by ASE (`fmax`) and VASP (`EDIFFG < 0`) by up to sqrt(3), and the
energy test kept runs going 4 to 23 evaluations after the force target was met.
Native FIRE applied the ASE dimensionless parameters (dt 0.1, dtmax 1.0) in
Hartree atomic units with unit masses, so its steps were about 97 times shorter
than ASE FIRE's for the same physical force.

## Decision

1. **Force measure.** `[relax] force_criterion = "norm"` (default): converged
   when max over atoms of |F_i| (Cartesian norm of the force after the
   constraint projection, Ha/Bohr) < `force_tol_ha_bohr`. `"component"`
   restores the old measure (largest free Cartesian component). Since
   component <= norm <= sqrt(3) component, the norm rule at equal tolerance is
   stricter by up to a factor sqrt(3). The inequality stays strict (as ASE).
   Fully frozen atoms contribute zero; partially frozen atoms contribute the
   norm of their free components.
   **Default tolerance:** `force_tol_ha_bohr` defaults to 0.02 eV/A on the
   norm, computed from `dft_units` (0.02 eV / 27.2114 eV/Ha x 0.529177 A/Bohr
   = 3.889e-4 Ha/Bohr). Previous default: 1e-3 Ha/Bohr on the component plus
   the 1e-8 Ha energy test. With the energy test gone the force threshold alone
   sets the accuracy, so the old 1e-3 would have stopped earlier and less
   accurately (Si8 smoke: 9 evaluations, final norm 9.1e-4 Ha/Bohr, 6.9e-5 Ha
   above the minimum); 0.02 eV/A is the PRELAX benchmark target and a common
   production setting. The default applies to `ions` and to the force part of
   `ions_cell`.
2. **Energy test optional.** `energy_tol_ha` is unset by default. Unset, in ions
   mode, the first frame that passes the force test stops the run, including
   the initial frame (exit 0 with `steps = 0`). Set, the pre-ADR rule applies:
   force test and |dE| < `energy_tol_ha` on two consecutive accepted steps.
   Trust-radius halving after an energy rise uses `energy_tol_ha` when set,
   otherwise twice the SCF energy accuracy (`scf.standard_threshold_ha`, or the
   strict energy tolerance on the oracle path); at the default 5e-9 Ha this is
   the former 1e-8 Ha, so BFGS/FIRE trajectories are unchanged until the stop.
3. **Cell modes** (`cell`, `ions_cell`): the stress test and the energy rule
   (two consecutive steps; enthalpy backtracking tolerance) are unchanged, with
   `energy_tol_ha` defaulting to 1e-8 Ha there as before. The force part of
   `ions_cell` follows `force_criterion` (default norm).
4. **FIRE units.** `fire_dt` (default 0.1) and `fire_dt_max` (default 1.0) are
   given in ASE FIRE units, Angstrom/sqrt(eV) (unit masses, forces in eV/A,
   positions in A), and internally scaled by s = sqrt(E)/L with E = 27.211386
   eV/Ha and L = 0.52917721 A/Bohr (s = 9.8576, s^2 = 97.17). Derivation: an
   Euler step moves x = dt^2 F. ASE: x_A = dt^2 F_eV/A. With F_eV/A = F_Ha/Bohr
   E/L and x_Bohr = x_A/L, x_Bohr = dt^2 (E/L^2) F_Ha/Bohr, i.e. the native time
   step dt s. The FIRE velocity mixing v <- (1-a) v + a |v| F/|F| and the power
   sign P = v.F are homogeneous in v and F, so this one factor makes the whole
   trajectory identical to ASE FIRE in eV/A (unit test
   `fire_trajectory_equals_ase_units_fire`, 60 steps, 1e-10 relative). The
   other parameters are the ASE ones (N_min 5, f_inc 1.1, f_dec 0.5, alpha 0.1,
   f_alpha 0.99), and two orderings now follow ASE: the first step is a pure
   kick (no power test, so dt is not halved on step 1), and dt grows when the
   count of positive-power steps before this one exceeds N_min. Nothing is
   tuned on a system: the factor is a unit conversion and the parameters are
   ASE's defaults. The old Hartree-unit behaviour is not kept; its time
   steps correspond to `fire_dt = 0.01014`, `fire_dt_max = 0.1014` (0.1/s and
   1/s), apart from the two ordering changes above.
   Remaining differences from ASE FIRE: the displacement cap is per atom
   (`max_displacement_bohr`, default 0.2 Bohr = 0.106 A, and the trust radius)
   instead of ASE's cap on the norm of the whole 3N step vector (`maxstep`
   0.2 A); a failed-SCF retry zeroes the velocity (the next step halves dt).
5. **Records.** Every frame gains `max_force_norm_ha_bohr` next to
   `max_force_ha_bohr` (component; meaning unchanged). results.json
   `relaxation` gains `final_max_force_norm_ha_bohr` and a `stop_rule` block
   (criterion, tolerances, consecutive steps required). `relax.log` appends a
   seventh column `max_force_norm_Ha_Bohr` (the first six are unchanged; the
   cell-mode log appends the same column). The per-step wall is printed on the
   stdout step line (` wall=...s`) and is derivable from `progress.jsonl`
   `elapsed_s`; it is deliberately not in frames, restart.json or relax.log,
   which stay deterministic so that a strict-oracle resume reproduces them
   bit for bit.

## Input compatibility

No key is removed or renamed. Old inputs parse unchanged. An input that sets
`energy_tol_ha` keeps the energy test (now opt-in); an input that omitted it
previously got the 1e-8 Ha energy test and now gets none. `force_tol_ha_bohr`
keeps its name; its default changes from 1e-3 to 3.889e-4 (0.02 eV/A) and its
measure to the norm unless
`force_criterion = "component"`. New keys: `force_criterion`, `fire_dt`,
`fire_dt_max` (validated: positive, finite, `fire_dt <= fire_dt_max`;
`energy_tol_ha`, when set, positive and finite). Unknown keys still reject
(`deny_unknown_fields`). Restart files of an older build were never accepted
(the restart identity includes the executable hash). The stored P44A
byte-identity gate (`reference_mode_relaxation`, needs `DFT_REFERENCE_BASELINE`)
compares against a pre-change executable and now differs by the added fields;
it is a historical gate for P44A, not for this change.

## Consequences

- Fewer evaluations with the same `force_tol_ha_bohr` in the common case
  (energy test gone), more where the norm rule binds (up to sqrt(3) on
  diagonal forces). Final structures are converged to a stricter force measure.
- FIRE is now comparable to ASE FIRE; it remains slower than BFGS on small
  bulk cells (audit section 5: ASE FIRE 1.6-2.4x ASE BFGS evaluations).
- Preconditioning (PRELAX Phase 2) is a separate ticket.
