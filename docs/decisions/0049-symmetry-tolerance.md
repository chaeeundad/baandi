# ADR 0049: Symmetry tolerance default and projection of the input onto the found group

- Status: **Accepted** (orchestrator, 2026-11-27, under the user's autonomy grant for general-path speed work; reported to the user) (ticket GS1-SYM-TOL, 2026-11-27).
- Scope: `[symmetry] mode = "auto"` only. `mode = "off"` (the default) is
  unchanged and bit-identical.

## Context

The dev08 Materials Project campaign ran hcp metals (Ru, Hf, Sr) whose
fractional coordinates are rounded to six decimals (0.333333 / 0.666667, about
1.7e-6 Bohr from the exact 2c/2d sites; the Angstrom lattice is rounded the same
way). With the P19A default `tolerance_bohr = 1e-7` the detector found 4
operations on the dev08 Ru input, while QE found 24 (nsym 24) on the identical
input, so BAANDI solved 3.7 to 3.8 times more irreducible k points (Ru: 244
against 64). The 1e-7 Bohr default was a "no tolerance" choice for synthetic,
exactly written inputs; it is below the rounding error of every published
structure database.

Common practice (public documentation, recalled, not re-fetched for this
ticket; no QE source was read): VASP `SYMPREC` default 1e-5 Angstrom
(1.9e-5 Bohr); spglib takes `symprec` in Cartesian length units, with 1e-5
Angstrom (phonopy) up to 1e-2 Angstrom (pymatgen `SpacegroupAnalyzer`) as usual
defaults; QE accepts an atom match at about 1e-5 in crystal coordinates (from
the task statement and QE output behaviour; for a 5 to 10 Bohr cell this is
5e-5 to 1e-4 Bohr). The repository's QE notes (`docs/research/qe-study/`) do not
cover symmetry.

## Decision

1. **Default `symmetry.tolerance_bohr = 1e-4` Bohr** (about 5.3e-5 Angstrom),
   constant `dft_geometry::symmetry::DEFAULT_TOLERANCE_BOHR`. Unit: Cartesian
   length, as before (the P19A detector contract in `docs/validation/symmetry.md` is a
   Cartesian distance with metric-derived fractional bounds). A Cartesian length
   is material- and cell-shape-independent: a fractional tolerance would grow
   with the cell (a slab with 60 Bohr of vacuum would accept 6e-4 Bohr along the
   normal) and shrink in small cells. 1e-4 Bohr lies between VASP's default and
   a 1e-5 crystal-coordinate criterion on a 10 Bohr cell, covers six-decimal
   rounding (5e-7 per fractional component) in cells up to roughly 50 Bohr, and
   is two to three orders of magnitude below any distortion that changes
   physics at the accuracy classes we gate (a 1e-4 Bohr displacement against a
   0.1 Ha/Bohr force is 1e-5 Ha at first order; but see 2: we remove it). The
   key stays user-settable; a larger value is the user's explicit choice and is
   recorded.
2. **Project the input onto the found group before the calculation**
   (preferred alternative). At parse time (`Config::parse` ->
   `Config::symmetrize_input`, `dft_geometry::symmetry::symmetrize`):
   - detect the group at the tolerance;
   - lattice: average the metric over the rotations, M' = mean R^T M R, and
     replace A by the symmetric stretch A' = A M^-1/2 M'^1/2 (no rigid rotation;
     the replaced lattice is written back in Bohr);
   - translations: least-squares estimate from all atoms; if every component
     snaps to a rational p/q (q <= 24) within the tolerance and the snapped set
     closes exactly as a group, the snapped set is used (this removes a small
     origin offset and keeps the translations commensurate with FFT grids);
     otherwise the estimate (an origin-shift-consistent set) is used;
   - positions: Reynolds average x_i' = mean_g (R_g x_j + tau_g) over all
     operations g mapping atom j to atom i;
   - refuse (error, exit 1 with the `symmetry.operations` diagnostic) if any
     atom or lattice vector would move by more than the tolerance, or if the
     result does not carry the group at a tight 1e-9 Bohr re-detection;
   - changes at or below 1e-12 Bohr are roundoff: the input bits are kept, so
     exactly written inputs give bit-identical results;
   - record in `results.json.symmetry`: `tolerance_bohr`,
     `crystal_operations_count`, and `input_symmetrization` (positions_changed,
     lattice_changed, translations_snapped, max_atom_displacement_bohr,
     max_lattice_change_bohr).
3. **Relaxation.** Fixed-cell and cell relaxations start from the projected,
   exactly symmetric geometry; forces and stress are symmetrized with the run's
   operations (existing P19A/P20A/P23A paths; the relaxation test asserts a
   constant operation count over the trajectory), so the trajectory keeps the start group to roundoff. The group is
   still re-detected per geometry at the configured tolerance (unchanged code);
   a trajectory that comes within the tolerance of a higher symmetry can gain
   operations. That is a documented limitation, not pinned in this ticket.

## Alternatives considered

- **Keep positions, symmetrize density, forces and stress only (QE-like).** The
  Hamiltonian is then not exactly invariant: k points of one orbit have
  eigenvalues that differ at first order in the residual asymmetry, and the
  reduced k sum weights only the representative, so the reduced energy differs
  from the full-mesh energy at first order (lattice rounding of 1e-7 relative
  can give 1e-8 to 1e-7 Ha, above the strict 1e-10 Ha class). Translations
  derived from rounded positions also miss FFT commensurability by ~1e-6, which
  our grid check (1e-8) rejects. Rejected for consistency.
- **Fractional (crystal-coordinate) tolerance.** Matches how databases round,
  but makes the physical tolerance depend on the cell size and anisotropy.
  Rejected; a user with a very large cell from a six-decimal source can raise
  `tolerance_bohr`.
- **Keep 1e-7 Bohr and ask users to clean inputs.** Leaves the 3.8x k-point
  penalty on every database structure. Rejected.
- **Symmetrize positions but not the lattice.** Leaves the first-order lattice
  term above. Rejected.

## Consequences and validation

- `crates/dft-geometry/tests/symmetry.rs`: six-decimal Ru-like and Mg-like hcp
  find 24 operations (inversion included) at the default and fewer at 1e-7;
  the projection lands on the exact sites to 1e-14 (fractional), carries the
  group at 1e-10 Bohr, is FFT-commensurate, and gives the same irreducible k
  count as the exact structure (Gamma-centred 6x6x4: 39 -> 21; 12x12x8:
  215 -> 95; 9x9x6: 100 -> 48 on this synthetic Ru-like record, where 1e-7 finds
  8 operations; dev08's Ru record found 4). Exact diamond, exact hcp and a
  low-symmetry skew cell are left bit-identical with the same group as at 1e-7.
  A 8.1e-3 Bohr (and a 2e-4 Bohr) distortion keeps the lower group; a uniform
  3e-6 origin offset is snapped away; an arbitrary origin (0.0123) is left
  untouched.
- `crates/dft-cli/tests/all/symmetry.rs`: six-decimal hcp synthetic run with
  symmetry auto against symmetry off on exact coordinates, strict path, 3x3x2
  mesh, Fermi-Dirac 0.01 Ha: energy difference 7.8e-12 Ha (strict class 1e-10).
- A dev08 rerun of the MP hcp cases measures the real k count and wall effect.
