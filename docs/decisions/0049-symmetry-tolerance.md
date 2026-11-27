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
   constant operation count over the trajectory), so the trajectory keeps the start group to roundoff.
   **Amendment GS1B (2026-11-27):** the loose tolerance is applied once, to
   the input, by the projection. Every later detection (run setup, k
   reduction, reference symmetry, every relaxation step, checkpoint and bands)
   uses `min(tolerance_bohr, 1e-9 Bohr)` (`Config::detection_tolerance`,
   recorded as `results.json.symmetry.detection_tolerance_bohr`). The
   projected input carries its group at 1e-9 Bohr (checked by `symmetrize`),
   so the start group is unchanged and relaxations keep it. Reason (dev08
   gen59, reproduced on the Mac in release): in
   `cell_relaxation::tetragonal_shape_returns_to_cubic` the tetragonal cell
   approached cubic; the cell optimizer writes each new lattice into the
   config, whose group was re-detected at 1e-4 Bohr, so once |c - a| fell
   below 1e-4 Bohr the cubic group was found, the SCF stress was symmetrized
   with cubic operations, the anisotropic stress that drives c -> a was
   averaged away, and the cell froze at |c - a| = 9.6e-5 Bohr (relative
   2.48e-5 against the 1e-5 test bound; with `tolerance_bohr = 1e-7` the same
   test reaches 6.4e-8). The tetragonal case of `reference_cell_protocol`
   stalled the same way (rebuilt pressure residual 4.4e-8 against
   1e-9 Ha/Bohr^3). Remaining limitation: a trajectory that comes within
   1e-9 Bohr of a higher symmetry can still gain operations; the driving force
   left at that distance is below every stress and force tolerance we gate.
   The cell optimizer's ionic and strain projectors are built once from the
   start group (unchanged).

4. **Amendment GS1C (2026-11-27): explicit FFT grids keep a subgroup.** With
   an explicit `basis.fft_grid`, the detected operations that do not permute
   the grid (rotation or fractional translation not commensurate with it) are
   dropped and the largest grid-preserving subgroup is used everywhere
   (`dft_geometry::symmetry::grid_subgroup`; `PreparedSymmetry::new_for_grid`;
   relaxation steps). The grid-preserving affine maps are closed under
   composition and inversion, so the kept subset of a group is a group; closure
   is re-checked and the identity-only group, always grid-preserving, is the
   fallback, so this never refuses. The input projection uses the kept subgroup
   only (translations chosen as in 2, then filtered), so atoms are not moved
   toward a symmetry that is then dropped. Preflight prints a
   `symmetry.fft_grid` warning and `results.json.symmetry` records
   `detected_operations_count`, `crystal_operations_count` (kept) and
   `dropped_operations` {count, reason}. Reason: dev08 MP campaign entry 18
   (Cl2, QE nsym 8, explicit grid copied from QE) was refused at 1e-4 Bohr
   while QE discards incommensurate operations and continues. Without an
   explicit grid, or with a compatible one, nothing changes (bit-identical).
   Tests: six-decimal hcp on [24,24,25] keeps 12 of 24 operations (no 6_3
   screw) and is a closed group; [24,25,25] keeps fewer; compatible grids give
   the identical projection; CLI run of the six-decimal hcp on an odd explicit
   grid against symmetry off on exact coordinates, same grid, strict path:
   2.2e-16 Ha. Engine-level `enable_symmetry` without a grid argument (oracle
   API) still rejects an incompatible grid.

5. **Amendment GS1D (2026-11-27): what "preserves the grid" means, and the
   Cl2 case.** An operation x -> R x + tau maps the grid {m_j / n_j} onto
   itself iff (a) n_i R_ij / n_j is an integer for all i, j (the image of one
   grid step along axis j must be a whole number of steps along every axis i)
   and (b) n_i tau_i is an integer (to 1e-8). Both are necessary and together
   sufficient; this is exactly the rule in use (`grid_rotation_compatible`,
   `grid_translation_compatible`), and `results.json.symmetry.dropped_operations`
   now splits the count into `rotation_count` and `translation_count`.
   Cl2 (MP campaign entry 18, Cmca) on the explicit grid [40, 36, 80]: the
   primitive cell has a1 = b_conv and a2 = (a_conv + b_conv)/2 with unequal
   lengths, hence n1 = 40 != n2 = 36. The mirrors and two-fold axes normal to
   x_conv or y_conv send a2 to +-(a1 - a2), i.e. R_12 = +-1, so a step of 1/36
   along a2 becomes a step of 1/36 along a1, which is not a multiple of 1/40.
   These four operations (two of them with a fractional translation) truly do
   not map this grid onto itself and are dropped by rotation, not by
   translation; every Cmca translation (halves) is commensurate. Fixture
   `cl2_like_explicit_grid_drops_only_operations_that_do_not_map_the_grid`
   reproduces 8 detected, 4 with translations, 4 kept, 4 dropped (rotation 4,
   translation 0), and shows that any grid with n2 dividing n1 ([40,40,80],
   [72,36,80]) keeps all eight. QE reports all 8 on [40,36,80], so it does not
   require the rotation part (a); presumably its density symmetrization works in
   reciprocal space, where any lattice rotation maps the G sphere onto itself.
   Our engine keeps (a) because the XC energy and potential are evaluated on the
   real-space grid: if a rotation does not permute the grid points, the discrete
   XC quadrature is not invariant under it (aliasing of the nonlinear f(rho)),
   so a reduced-k run would no longer reproduce the full-mesh run exactly.
   Accepting such operations (QE-like) is a separate, measurable change: a
   reciprocal-space density symmetrizer for non-grid-preserving rotations plus
   a measured reduced-versus-full energy difference against the accuracy
   classes. It is not made here.

6. **Amendment GS1D-B (2026-11-27, ticket GT2B-GRID-OPS): rotations that do
   not map an explicit grid are kept (switch `GRID_INCOMPATIBLE_OPS`, default
   `reciprocal`; `drop` restores item 5).** Cause: dev08 gen64 Ta6Zn7
   (explicit [60,60,108]) kept 4 of 12 operations and Hg ([60,60,54]) 2 of 8,
   while QE keeps 12 and 8 on the same grids. A hand check (every grid point
   mapped through x -> R x + tau, also with R^T) confirmed that the dropped
   rotations really move grid points off the grid: Ta6Zn7 rotations with
   R_13 = -1 send (0,0,1/108) to (-1/108,0,1/108), -0.56 steps along a1; the
   Hg body-centred rotations mix all three axes with 54 not dividing 60. All
   their translations are exact after the input projection (Ta6Zn7: 0; Hg:
   (3/4,1/4,1/2), N tau = (45,15,27)); the 1e-8 translation tolerance is not
   involved. Decision with `reciprocal`:
   - an operation is kept when its fractional translation is commensurate
     with the grid (translation-incompatible operations are dropped as
     before, as in QE); its rotation need not permute the grid;
   - the density is projected onto the whole group in G space on the density
     sphere (the existing projector; it never needed the grid to be permuted,
     only the guard required it);
   - the local potential (ionic + Hartree + XC, each collinear channel) is
     projected onto the group in G space after every evaluation
     (reference maps: `refscf::map::fields` and the spin fields; strict
     oracle: `input_fields`, collinear channels), so the Kohn-Sham Hamiltonian
     is exactly invariant although the real-space XC quadrature is not.
     Components outside the density sphere are zeroed; on an alias-free grid
     they do not enter Hamiltonian matrix elements, and the double-counting
     integral is unchanged because the density is invariant and confined to
     the sphere. Only runs where some kept operation does not permute the grid
     project the potential (`needs_potential_projection`); all other runs are
     bit-identical;
   - k reduction (point-wise, ADR 0051), forces, stress and relaxation use
     the whole group (they read the same operation list);
   - record: `results.json.symmetry.grid_incompatible_ops` {policy,
     explicit_fft_grid, rotations_not_permuting_grid_kept,
     potential_projection, rule}; preflight names the kept operations.
   Validation (`crates/dft-cli/tests/spin.rs
   grid_incompatible_rotations_reciprocal_matches_projected_full_mesh`,
   release; raw output `tests/reference/symmetry/gs1db-grid-incompatible-ops.txt`):
   synthetic Z = 1 nonlocal potentials, 4 Ry, Fermi-Dirac 0.02 Ha, forces and
   stress; arms full (symmetry off, raw full mesh), unreduced (full mesh with
   density and potential projected: the exact oracle), reciprocal, drop.

   | cell, grid, k | ops kept drop / reciprocal | nk drop / reciprocal / full | (b) reciprocal - unreduced: dE, dF, dsigma | (a) reciprocal - full: dE, dF, dsigma |
   |---|---|---|---|---|
   | Ta6Zn7 [60,60,108], 2x2x2 | 4 / 12 | 6 / 4 / 8 | 2.6e-13 Ha, 5.2e-12 Ha/Bohr, 2.0e-14 Ha/Bohr^3 | 2.7e-11, 1.4e-6, 6.8e-10 |
   | Ta6Zn7 [20,20,36], 2x2x2 | 4 / 12 | 6 / 4 / 8 | 7.6e-14, 5.6e-12, 1.8e-14 | 3.1e-11, 1.4e-6, 7.0e-10 |
   | Hg [60,60,54], 2x2x2 | 2 / 8 | 6 / 4 / 8 | 1.1e-14, 2.4e-16, 3.4e-17 | 3.1e-14, 1.2e-11, 6.9e-14 |
   | Hg [20,20,18], 2x2x2 | 2 / 8 | 6 / 4 / 8 | 2.2e-14, 8.4e-16, 6.6e-17 | 0, 9.9e-10, 4.8e-11 |
   | Cl2-like [40,36,80], 2x2x1 | 4 / 8 | 4 / 3 / 4 | 7.0e-14, 1.3e-14, 2.3e-16 | 1.1e-13, 5.2e-11, 5.2e-16 |
   | Cl2-like [20,18,40], 2x2x1 | 4 / 8 | 4 / 3 / 4 | 1.1e-13, 1.3e-13, 7.7e-16 | 1.4e-13, 4.2e-11, 1.0e-13 |

   All meshes here are invariant under the group, so (a) isolates the cost of
   the non-invariant XC quadrature (largest on Ta6Zn7: 3e-11 Ha, 1.4e-6
   Ha/Bohr). The standard reference driver (threshold 1e-12 Ha) on the scaled
   grids gives (b) <= 7.1e-13 Ha, 1.5e-8 Ha/Bohr, 2.6e-11 Ha/Bohr^3. Default
   `reciprocal` because (b) holds within the strict class and the force and
   stress bounds on every fixture and path. With a non-invariant mesh the
   (a) difference also contains the ADR 0051 quadrature change (Ta6Zn7 at
   2x2x1, first Mac run, not in the stored log: 1.4e-5 Ha, 4.4e-3 Ha/Bohr,
   while (b) stayed at 1.6e-13 Ha).

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
