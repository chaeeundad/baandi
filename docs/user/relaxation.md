# Fixed-cell atomic relaxation

Current routing: [ADR 0047](../decisions/0047-reference-driver-single-path.md).
Production ionic and variable-cell steps use the reference driver; strict is an
explicit oracle. The optimizer and its restart are separate from checkpoint v2.

`baandi relax input.toml --output-dir relaxed` uses complete P20A
forces from verified reference states (or strict-certified oracle states) to move atoms while retaining the input cell. This initial driver supports
unpolarized, fixed insulating occupations, Gamma or the existing k meshes,
LDA/PBE, NC nonlocal projectors, optional NLCC and auto symmetry.
Collinear and finite-temperature drivers reject explicitly. For cell modes, see
[variable-cell relaxation](cell-relaxation.md). The ordinary `run` path is unchanged.

```toml
[relax]
method = "bfgs"                 # alternatively "fire"
force_tol_ha_bohr = 3.889e-4    # default (= 0.02 eV/A), on the max atomic force norm |F_i| (Ha/Bohr)
force_criterion = "norm"        # default; "component" = largest free Cartesian component
# energy_tol_ha = 1e-8          # optional (unset by default): also require |dE| < tol
max_steps = 200
max_displacement_bohr = 0.2
trust_radius_bohr = 0.2
warm_start = true
fixed = [{atom = 0}, {atom = 2, components = [false, true, true]}]
```

Atom indices are zero-based in input order. Omitted `components` fixes all three
Cartesian components. `true` fixes that component; `false` leaves it free.
The same constraints can be written per atom in the compact `[structure]` form
(`Si 0 0 0 fixed=xyz`, `fixed=z`, ...) or come from POSCAR selective dynamics or
an extended XYZ `move_mask`; give them in one place only
([Structure input](../cli.md#structure-input)).
Constraints project both optimizer forces and steps; the original physical force
arrays remain unchanged in the trajectory/results.

Stop rule ([ADR 0056](../decisions/0056-relax-stop-rule-and-fire-units.md),
2026-11-27, the ASE `fmax` / VASP `EDIFFG < 0` convention): the run is converged
when the **largest atomic force norm** |F_i| of the projected force is below
`force_tol_ha_bohr`, default 0.02 eV/A = 3.889e-4 Ha/Bohr (before ADR 0056:
1e-3 Ha/Bohr on the component plus a 1e-8 Ha energy test). Fully frozen atoms count
zero, partially frozen atoms the norm of their free components.
`force_criterion = "component"` selects the rule used before ADR 0056 (largest
free Cartesian component); the norm is larger by up to sqrt(3), so at equal
tolerance the norm rule is stricter. There is no energy requirement by
default: the first frame that passes the force test stops, the initial one
included (exit 0 with zero steps). Setting `energy_tol_ha` adds the former
energy test: the force test and |dE| < `energy_tol_ha` must then hold on two
consecutive accepted steps. Inputs that set `energy_tol_ha` explicitly keep
that behaviour. The initial SCF is step 0;
`max_steps` counts subsequent ionic moves. All force and energy tolerances are
strict inequalities. Variable-cell modes keep their stress and energy rule
(energy default 1e-8 Ha there), see [cell relaxation](cell-relaxation.md). Exit 0 means converged, 3 means ionic step budget exhausted,
4 means SCF/derivative failure after the retry, 1 invalid input, 5 I/O failure.
An exhausted ionic trajectory still contains certified SCFs and usable forces.

Since ONEPATH-1 (user decision 2026-11-16) every ionic step of a fixed-cell
relaxation runs the single standard SCF path with the fixed-point verification
(ADR 0043); forces are computed only from a step whose verification passed, and
a failed verification continues that step's SCF. The SCF accuracy is
`scf.standard_threshold_ha` (default 5e-9 Ha); `scf.verify = "off"` is an input
error. The strict path with its certificate remains as a reference tool
(`BAANDI_REFERENCE_STRICT=on`). Since ONEPATH-2 the variable-cell modes (`cell`,
`ions_cell`) run the same standard path: forces and stress of every step come
from its verified state, and the stress cutoff scans are verified standard runs
(see [cell relaxation](cell-relaxation.md)); without an explicit
`scf.standard_threshold_ha` they default to 5e-11 Ha (stress default,
provisional), while the ions mode starts from the class default 5e-9 Ha. **Force-tied default
([ADR 0057](../decisions/0057-relax-force-tied-scf-threshold.md), 2026-11-28):**
when `scf.standard_threshold_ha` is not set, `ions` and `ions_cell` relaxations
run at min(5e-9 Ha or the stress default, `force_tol_ha_bohr`^2 / 3200), floor
5e-14 Ha: 4.73e-11 Ha at the default 0.02 eV/A. Reason: warm-started ionic SCFs
stop with a delivered density residual close to the threshold, and at 5e-9 Ha
that gave force errors up to 2.3e-4 (MgO) and 2.8e-4 Ha/Bohr (ZnO), more than
half the default tolerance; under the rule the measured errors were at most
2.6e-5 Ha/Bohr (0.07 of the tolerance) on Si8, Si8-AlP, NaCl, ZnO and MgO, for
9 to 21 percent more H applications per relaxation.
`single_path.threshold_source = "relax_force_default"` records it; an explicit
threshold always wins. Tighten `standard_threshold_ha` to resolve a
small ionic energy threshold or tight force tolerance (for example 5e-11 Ha, or
5e-13 Ha for `energy_tol_ha = 1e-12`); nothing is loosened automatically. A warning (stderr and `progress.jsonl`) is printed when `energy_tol_ha` is below twice the SCF threshold, where the energy test and the forces near convergence are limited by SCF noise (RELAX-GATE-FIXES: MgO PBE at 5e-9 Ha scattered forces by about 2e-4 Ha/Bohr in 1-iteration warm SCFs).
Legacy `spectral_guard` and `xc.density_floor` input keys were removed in R10c.
No tolerances are automatically loosened.

## Optimizers and symmetry

BFGS stores a dense inverse Hessian. Its initial independent-spring model has
stiffness 1 Ha/Bohr^2, hence inverse Hessian I Bohr^2/Ha. Positive-curvature
secant pairs update it; unresolved/nonpositive curvature pairs are skipped.
A nonfinite or uphill search direction resets the spring model. Uniform scaling
limits the largest **per-atom Euclidean displacement** by both the current trust
radius and `max_displacement_bohr`. The trust radius halves after an energy rise
larger than `energy_tol_ha` (when unset: twice the effective SCF threshold,
9.5e-11 Ha with the force-tied default of ADR 0057, 1e-8 Ha at 5e-9 Ha) or a failed-SCF retry; otherwise it
grows by 1.1 up to its input value. There is no energy line search or global
minimum guarantee. This is a capped, safeguarded BFGS step, not a solution of an
exact quadratic trust-region subproblem. Energy monotonicity is measured on the
Si oracle, not guaranteed for every system.

FIRE is a selectable alternative, not an automatic switch. Since ADR 0056 it
moves atoms as ASE FIRE does in eV/Angstrom: unit fictitious masses, an Euler
velocity kick, `fire_dt` (initial, default 0.1) and `fire_dt_max` (default 1.0)
in ASE units Angstrom/sqrt(eV), internally multiplied by
sqrt(27.2114)/0.529177 = 9.858 for Hartree atomic units (displacement
dt^2 F, so the same dt gives the same step as ASE for the same physical force),
alpha=0.1. The first step is a pure kick. Positive power mixes velocity toward
the force; once more than five consecutive positive-power steps preceded, dt
grows by 1.1 up to `fire_dt_max` and alpha shrinks by 0.99. Nonpositive power
clears velocity, halves dt and restores alpha. Its displacements use the same
per-atom cap (ASE caps the norm of the whole step at 0.2 A instead). Before ADR
0056 the same numbers were applied in Hartree units, which made the steps about
97 times shorter. These fictitious integration parameters do not represent
physical dynamics.

With `[symmetry] mode="auto"`, P19A/P20A symmetrize density and raw forces.
The isotropic starting Hessian, invariant updates and scalar cap preserve the
selected group up to arithmetic roundoff. Cartesian constraints must commute
with every selected spatial rotation and atom permutation. A constraint breaking
the group rejects before SCF; use symmetry off or an invariant set of constraints.

## Density transfer, failure and restart

Every accepted structure starts a new reference SCF through `refscf::api`.
Within a live calculation, `warm_start=true` retains an `IonicSeed` (density and
compatible orbitals) from the preceding delivery, rebuilding potentials and
clearing SCF history; it is not exact electronic replay. Model, basis, grid,
charge and symmetry compatibility are checked. Cell transfer uses
`IonicSeed::in_cell` with the frozen reference support. `warm_start=false` uses
fresh atomic starts. StrictOracle retains its legacy density-transfer API.

The optimizer restart currently does **not** persist the reference orbital seed.
After process restart the reference route records
`reference restart has no orbital seed; fresh atomic start`. Historical bytewise
optimizer-resume evidence applies to its original legacy path; do not infer
bit-identical reference continuation from it. Exact electronic continuation is
separately provided by [checkpoint v2](../checkpoint-format.md), not `relax --resume`.

Optional start extrapolation (C3, ions mode with `warm_start=true`):
`start_extrapolation = "none"` (default, the transfer above), `"first_order"` or
`"second_order"`, overridden by `BAANDI_RELAX_EXTRAPOLATION`. The transferred
density becomes n_at(R_new) + d, where n_at is the superposition of the UPF
atomic charges (`PP_RHOATOM`) and d = n_scf - n_at of the previous step (first
order) or its least-squares extrapolation along the last three ionic steps
(second order). Negative samples of the prediction are clipped and the charge
restored (`start_clipped_fraction` per frame); without atomic charges the plain
density is used (`start_fallback`). Provenance is `extrapolated_first_order` or
`extrapolated_second_order`; certificates and forces are unchanged. The predictor history
is stored in restart.json; this does not persist the reference orbital seed or
guarantee identical subsequent reference SCF trajectories. Every frame records
`scf_h_vectors`. Evidence and dev08 request:
`docs/performance/requests/c3-ionic-extrapolation.md`.

If a candidate SCF fails or exhausts its budget, the driver retries exactly once
from a fresh density at half the proposed displacement from the last good state.
If that also fails, the result status is `failed`, the failure records retain
attempted coordinates/reasons, and all final scientific outputs retain the last
good structure. No forces are queried for an uncertified state. A failure of the
initial SCF is retried fresh at the same structure; if neither succeeds, there
is no last good structure and the result has no numeric force payload.

After every accepted step, `restart.json` is published by same-directory
exclusive temporary creation, file fsync, atomic rename and directory fsync.
The checksum-protected payload stores all f64 values as integer IEEE-754 bits,
including coordinates, optimizer memory, last density and accepted frames.
Density storage is necessary for identical future warm-start arithmetic, so the
file scales with the density grid and history; it is not a tiny geometry-only
record and is not a BAANDICP electronic checkpoint. The 512 MiB import bound
rejects oversized artifacts. Checksums establish integrity, not authenticity.

```sh
baandi relax input.toml --output-dir relaxed --resume
```

Use the same input bytes, UPFs, executable, architecture and OS. Resume rejects
identity/shape/checksum changes and does not silently change optimization limits.
It regenerates derived outputs from the last committed optimizer boundary,
avoiding duplicate trajectory frames after interrupted output publication.
The historical legacy same-build ionic continuation test is bit-identical;
reference resume has the missing-orbital-seed limitation above. Each later SCF
is a new solve, not exact electronic continuation. Initial-SCF
failure produces only `results.json`, since no resumable good boundary exists.
A finished failed/unconverged run is not silently restarted with a new budget.

## Files and units

- `trajectory.extxyz`: input atom order, periodic cell, positions in Angstrom,
  energies in eV, raw forces in eV/Angstrom; standard ASE extended XYZ fields.
- `final.xyz`: elements and Angstrom coordinates; loses cell, constraints,
  energy, forces, PP identity and species-label distinctions.
- `POSCAR`: scale 1, lattice rows and Cartesian coordinates in Angstrom,
  species grouped in first-appearance order. Loses original interleaved atom
  order, PP identity, constraints and energies/forces; selective-dynamics flags
  are deliberately absent because the requested constraints are Cartesian.
- `final.toml`: replacement `[[system.atoms]]` entries in the original cell,
  with original species labels and unwrapped fractional coordinates.
- `results.json`: additive `relaxation` block, settings, constraints, status,
  all accepted energies/raw forces/positions, projected maxima per frame
  (`max_force_ha_bohr` = largest free component, `max_force_norm_ha_bohr` =
  largest atomic norm), `final_max_force_norm_ha_bohr`, the applied `stop_rule`,
  failure records, input/executable/PP SHA-256 provenance. Internal record
  units are Ha and Bohr.
- `relax.log`: step, energy (Ha), projected maximum force component (Ha/Bohr),
  maximum atom displacement (Bohr), SCF iterations, fresh/warm provenance and
  (seventh column, ADR 0056) projected maximum atomic force norm (Ha/Bohr).
- stdout: one line per step with both force measures and the step's wall
  (`wall=...s`); the wall is not stored in the deterministic files above
  (`progress.jsonl` carries `elapsed_s` per step).

Existing output files require `--resume`; symlink, nonregular, hard-linked and
input-alias destinations reject. Ordinary XYZ and POSCAR are lossy views; retain
the original input and JSON records for complete calculation provenance.
ASE independently reads the extxyz and POSCAR outputs in the acceptance test.

ASE optimizers continue to use `Baandi.get_forces()`. Unsupported properties now
raise ASE's `PropertyNotImplementedError` (a `NotImplementedError` subclass),
allowing BFGS to detect that `free_energy` is unavailable and use ordinary energy.
No fabricated free-energy or stress alias is introduced.

See [validation evidence](../validation/relaxation.md) for measured scope and limits.
