# Fixed-cell atomic relaxation

`dft-cli relax input.toml --output-dir relaxed` uses complete, certified P20A
forces to move atoms while retaining the input cell. This initial driver supports
unpolarized, fixed insulating occupations, Gamma or the existing k meshes,
LDA/PBE, NC nonlocal projectors, optional NLCC/density floor and auto symmetry.
Collinear and finite-temperature drivers reject explicitly; stress/cell relaxation
are not implemented. The ordinary `run` path is unchanged.

```toml
[relax]
method = "bfgs"                 # alternatively "fire"
force_tol_ha_bohr = 1e-3
energy_tol_ha = 1e-8
max_steps = 200
max_displacement_bohr = 0.2
trust_radius_bohr = 0.2
warm_start = true
fixed = [{atom = 0}, {atom = 2, components = [false, true, true]}]
```

Atom indices are zero-based in input order. Omitted `components` fixes all three
Cartesian components. `true` fixes that component; `false` leaves it free.
Constraints project both optimizer forces and steps; the original physical force
arrays remain unchanged in the trajectory/results. The maximum force used for
stopping is the maximum absolute **free Cartesian component**, not the maximum
atom-vector norm. Both this force threshold and the absolute energy-change
threshold must pass on two consecutive accepted steps. The initial SCF is step 0;
`max_steps` counts subsequent ionic moves. All force and energy tolerances are
strict inequalities. Exit 0 means converged, 3 means ionic step budget exhausted,
4 means SCF/derivative failure after the retry, 1 invalid input, 5 I/O failure.
An exhausted ionic trajectory still contains certified SCFs and usable forces.

The standard SCF profile is rejected. The P20A force requirements still apply:
SCF density <=1e-7, energy/eigen tolerances <=1e-8 Ha and a fresh strict certificate.
Tighter SCF tolerances may be necessary to resolve a small ionic energy threshold.
The spectral guard remains available. No tolerances are automatically loosened.

## Optimizers and symmetry

BFGS stores a dense inverse Hessian. Its initial independent-spring model has
stiffness 1 Ha/Bohr^2, hence inverse Hessian I Bohr^2/Ha. Positive-curvature
secant pairs update it; unresolved/nonpositive curvature pairs are skipped.
A nonfinite or uphill search direction resets the spring model. Uniform scaling
limits the largest **per-atom Euclidean displacement** by both the current trust
radius and `max_displacement_bohr`. The trust radius halves after an energy rise
larger than the configured energy tolerance or a failed-SCF retry; otherwise it
grows by 1.1 up to its input value. There is no energy line search or global
minimum guarantee. This is a capped, safeguarded BFGS step, not a solution of an
exact quadratic trust-region subproblem. Energy monotonicity is measured on the
Si oracle, not guaranteed for every system.

FIRE is a selectable alternative, not an automatic switch. It uses unit
fictitious masses, an Euler velocity kick, initial dt=0.1 and alpha=0.1. Positive
power mixes velocity toward the force; after five positive-power steps dt grows
by 1.1 up to 1 and alpha shrinks by 0.99. Nonpositive power clears velocity,
halves dt and restores alpha. Its displacements use the same per-atom cap.
These fictitious integration parameters do not represent physical dynamics.

With `[symmetry] mode="auto"`, P19A/P20A symmetrize density and raw forces.
The isotropic starting Hessian, invariant updates and scalar cap preserve the
selected group up to arithmetic roundoff. Cartesian constraints must commute
with every selected spatial rotation and atom permutation. A constraint breaking
the group rejects before SCF; use symmetry off or an invariant set of constraints.

## Density transfer, failure and restart

Every accepted structure starts a new SCF. `warm_start=true` passes the previous
evaluated density, with no orbitals, mixer, iteration history or certificate.
The engine compares independently retained model descriptions, allowing only
ordered atom positions to differ. Cell, ordered species, full parsed PP content,
cutoff, XC, density floor, grid, sampling/weights, occupations and symmetry must
match. Shape, C-contiguous array order, finiteness, XC negative-roundoff allowance and relative charge
error <=1e-8 are checked before initialization. Density is not renormalized or
silently projected. A change of detected symmetry may require a fresh retry.
Frame provenance is `warm_transfer`, never exact SCF replay; `warm_start=false`
uses fresh SCFs throughout. The additive engine entry point is
`ScfSession::initialize_density_transfer(source_problem, initial_density)`.

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
dft-cli relax input.toml --output-dir relaxed --resume
```

Use the same input bytes, UPFs, executable, architecture and OS. Resume rejects
identity/shape/checksum changes and does not silently change optimization limits.
It regenerates derived outputs from the last committed optimizer boundary,
avoiding duplicate trajectory frames after interrupted output publication.
Same-build ionic continuation is bit-identical in the tested case; each later
SCF still uses warm transfer, not exact electronic continuation. Initial-SCF
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
  all accepted energies/raw forces/positions, projected maxima, failure records,
  input/executable/PP SHA-256 provenance. Internal record units are Ha and Bohr.
- `relax.log`: step, energy (Ha), projected maximum force (Ha/Bohr), maximum
  atom displacement (Bohr), SCF iterations and fresh/warm provenance.

Existing output files require `--resume`; symlink, nonregular, hard-linked and
input-alias destinations reject. Ordinary XYZ and POSCAR are lossy views; retain
the original input and JSON records for complete calculation provenance.
ASE independently reads the extxyz and POSCAR outputs in the acceptance test.

ASE optimizers continue to use `Baandi.get_forces()`. Unsupported properties now
raise ASE's `PropertyNotImplementedError` (a `NotImplementedError` subclass),
allowing BFGS to detect that `free_energy` is unavailable and use ordinary energy.
No fabricated free-energy or stress alias is introduced.

See [validation evidence](../validation/relaxation.md) for measured scope and limits.
