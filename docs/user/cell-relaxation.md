# Variable-cell relaxation

`dft-cli relax input.toml --output-dir relaxed` accepts these additional settings:

```toml
[relax]
mode = "ions_cell"              # ions (default), cell, ions_cell
method = "bfgs"                 # or fire
cell_dofs = "full"              # full, isotropic, xyz_diagonal
fixed_axes = [false, false, false]
pressure_ha_bohr3 = 0.0
stress_tol_ha_bohr3 = 1e-6
pulay_warning_ha_bohr3 = 5e-6
rebuild_every = 0
```

`ions` retains the existing fixed-cell driver. `cell` retains fractional ionic
positions; `ions_cell` varies fractional positions and strain together. The
objective is enthalpy **H = E + P V**, with constant external pressure P in
Ha/Bohr^3. Positive P compresses. Internal tensile-positive stress is sigma;
internal pressure is -trace(sigma)/3. Equilibrium requires projected sigma+P I
and free ionic forces to vanish. Both derivative thresholds and the enthalpy
change tolerance must pass for two consecutive accepted moves.

`full` uses all six symmetric Cartesian strain components; rotation is excluded.
`isotropic` preserves shape and angles. `xyz_diagonal` permits independent
Cartesian x/y/z stretches: angles remain fixed for an orthogonal Cartesian cell,
but not for arbitrary skew lattice vectors. `fixed_axes` removes all strain
components involving a marked Cartesian axis. Isotropic strain combined with
any fixed axis rejects. Cartesian atom constraints (`relax.fixed`) currently
reject in the variable-cell modes; use the fixed-cell driver for those constraints.

The supported electronic models are those of the fixed-cell driver: unpolarized
fixed occupations with complete strict force/stress certificates. No standard
profile, spin/FD driver, ESLE continuation or uncertified derivative is admitted.
Cell runs use the requested strict eigen tolerance from the initial solve,
disabling residual-based eigen scheduling in this path. Ordinary SCF and ionic
relaxation retain their existing schedules.

## Coordinates, gradients and symmetry

Let A0 be the initial lattice matrix, L=V0^(1/3), D=I+q_cell/L, and A=D A0.
The generalized ionic coordinates are q_I=L s_I, where s_I are fractional
coordinates. The optimizer stores nine tensor components and projects onto
symmetric tensors; off-diagonal components are both present, with no implicit
Voigt factor. The exact differential at any current D is

    dH = sum_I [-(A^T F_I)/L] . dq_I
         + [V/L (sigma + P I) D^(-T)] : dq_cell.

Thus the generalized negative gradients supplied to BFGS/FIRE are A^T F/L and
-V/L (sigma+P I) D^(-T). V is the current volume. L is fixed throughout a run,
including basis rebuilds; optimizer coordinates and secant pairs use one chart.
The stress stopping measure is the maximum component of the projected cell
negative gradient times L/V. It is a constrained generalized residual in
Ha/Bohr^3; raw physical stress remains separately available.

With symmetry auto, P19A symmetrizes the SCF/derivatives. Cell directions use the
initial group's Cartesian rotation projector. Fractional displacement rotations
need not be orthogonal, so the driver orthonormalizes the invariant displacement
subspace and projects both gradients and steps there. Cell constraints must
commute with the symmetry projector or they reject. The initial group is
preserved; enhanced symmetry can be detected at later structures. The final
operation count is reported without assigning a space-group symbol.

BFGS/FIRE reuse the P21A implementation. A scalar cap additionally bounds actual
Cartesian atomic displacement, including affine cell motion. Non-positive-definite
cell deformations reject. Up to 16 halvings backtrack an enthalpy rise, allowing
only the configured energy tolerance as SCF noise. An SCF failure receives one
fresh half-step retry. Failure preserves the last accepted structure and records
the attempted generalized coordinates/reason. Initial failure writes a result
with no numeric last-good structure; periodic/final rebuild failures are explicit.

## Frozen basis, rebuild and warm transfer

Every optimization segment freezes the reference cell's ordered Miller support
for both wavefunction and density and freezes FFT dimensions. Reciprocal lengths
change with the cell. Under isotropic contraction the effective maximum kinetic
energy increases; expansion decreases it. Anisotropic changes have no single
scalar effective cutoff. This is the smooth finite-basis enthalpy differentiated
by frozen stress, not a fixed-cutoff or cutoff-converged energy curve.

At completion, including step-budget exhaustion, a fresh single point rebuilds
at the final cell with the nominal cutoff and an appropriate FFT grid. The JSON
contains both frames, their energy/force/stress/pressure differences, and a
warning if the largest stress difference exceeds `pulay_warning_ha_bohr3`.
The rebuilt frame also includes P22A's independent isotropic cutoff-derivative
Pulay estimate and uncertainty. The rebuild difference and that estimate are
distinct diagnostics; neither proves cutoff convergence.

`rebuild_every=N` optionally begins a new smooth segment after N accepted moves.
It resets optimizer memory and consecutive-convergence counting, starts fresh,
and records a `basis_rebuild` frame. Those frames do not count as moves. Energy
monotonicity applies within a segment; changing finite-basis functionals can
introduce a jump. Default zero performs only the final rebuild.

Warm starts retain samples on the identical fractional grid and multiply density
by V_old/V_new. This preserves quadrature charge without interpolation. Ordered
species, full parsed potentials, XC, sampling, occupations, frozen reference,
symmetry and grid must still match. Changed support rejects this transfer.
Shape, finite/XC-feasible density and relative charge <=1e-8 are checked before
initialization. This transfers no mixer, orbitals or certificate. Provenance is
`warm_transfer`, never an exact electronic continuation.

`--resume` restores a checksum-protected, float-bit optimizer restart containing
q, initial/reference cells, grid, optimizer memory, accepted frames, next density,
move count and any final rebuilt result. Same-build/input/PP identity is required.
The independent continuation test compares output files byte for byte. This is
cell optimizer continuation, not an ESLE or BAANDICP electronic checkpoint.

## Outputs and limits

- `trajectory.extxyz` carries each frame's own cell, Cartesian positions, energy
  and physical forces in Angstrom/eV units.
- `POSCAR`, `final.xyz` and `final.cif` describe the last accepted geometry.
  CIF uses core CIF lengths/angles and explicit fractional sites in P1, preserving
  all atoms without claiming a detected space-group symbol. XYZ loses the cell.
  These formats lose pseudopotential identity and optimizer state.
- `final.toml` replaces the original `system` and atom sections, including lattice.
- `relax.log` includes E, H, volume, pressure, force, displacement and provenance.
- `results.json` includes frozen/rebuilt frames, pressure and basis diagnostics,
  detected spatial operation count and failed trials. Retain it and the original
  input for full provenance. Restarts store the optimization boundary before
  derived output publication.

No agreement with experiment or another electronic-structure code is claimed.
The real-file oracle uses the P13D shifted meshes (Si 4x4x4 at 30/40 Ha,
MgO 3x3x3 at 50/60 Ha), PBE and the same external UPFs. It constrains the cell
to isotropic strain for a cubic lattice-constant comparison and reports cutoff
variation; it is not bulk qualification. See [validation](../validation/cell-relaxation.md).

### Final rebuilt-pressure warning

`relaxation.rebuilt_pressure_check` separately compares the final nominal-cutoff
pressure with `pressure_ha_bohr3`, using the unchanged `stress_tol_ha_bohr3`.
At or above that tolerance, the JSON contains a warning and actionable message,
and the CLI prints the warning to stderr. This is distinct from the tensor
rebuild-difference warning controlled by `pulay_warning_ha_bohr3` (default 5e-6).
A converged frozen optimization does not imply a nominal-cutoff equilibrium;
`rebuild_every=0` remains the default. This additional check does not silently
reoptimize a frozen-only result or loosen its tolerance.

### Recommended material-equilibrium protocol

For an equilibrium at the nominal cutoff, use `rebuild_every=3` (the N=5
variant is also covered by the material validation), then check the **rebuilt**
pressure rather than only the frozen pressure. Keep the cutoff, k mesh, XC,
pseudopotentials and tolerances identical between starting structures.

A fresh-basis polish consists of a new `relax` calculation whose starting
geometry is the preceding accepted final structure. Replace the `system`/atom
sections from `final.toml`, keep the electronic model unchanged, and use a new
output directory. Use `rebuild_every=0` within each short polish cycle; each new
cycle independently selects the nominal-cutoff basis at its starting cell.
This is not `--resume`: changing the basis requires reset optimizer/electronic
history, rather than exact continuation of the earlier finite-basis function.

The material validation executes two polish cycles, then a third only if needed.
Acceptance requires the final rebuilt pressure to meet the unchanged stress
tolerance and the two starts to agree within the previously estimated lattice
tolerance. If three cycles do not suffice, report the result as inconclusive;
do not increase tolerances. For nonzero applied pressure, compare Pr with the
specified external pressure. Atomic degrees of freedom also require their force
criterion when enabled. The current material oracle holds fractional sites fixed
and uses isotropic cell freedom.

This protocol is implemented reproducibly in
`tests/reference/relaxation/cell/rebuild_protocol.py`, invoked by the ignored
`real_files_rebuild_start_agreement` test. Its frozen-only input/evidence is kept
separate. The original `rebuild_every=0` default is unchanged: it minimizes the
initial Miller set's smooth function, which depends on the starting cell. Such
runs remain useful derivative oracles, but two different frozen functions need
not have identical minima at a finite cutoff.
