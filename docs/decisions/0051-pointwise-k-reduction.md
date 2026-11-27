# ADR 0051: Point-wise k reduction with the full crystal group

- Status: **Proposed** (ticket GT2-FULL-GROUP-K, 2026-11-27; orchestrator review pending).
- Scope: `[symmetry] mode = "auto"` with a k grid, every model (unpolarized,
  collinear, fixed and Fermi-Dirac occupations, forces, stress, relaxation).
  Switch `BAANDI_K_REDUCTION` (`[performance] k_reduction`): `pointwise`
  (default), `mesh` (former behaviour), `unreduced` (oracle). Table 94 -> 95.
- Follows the open item of ADR 0050 ("Follow-up").

## Context

The dev08 Materials Project campaign (gen61) ran Ta6Zn7 (13 atoms, 6x6x3
Gamma mesh) with 34 irreducible k in BAANDI against 26 in QE 7.6, both with
12 operations. In the Materials Project cell setting only 4 of the 12
operations map the whole 6x6x3 mesh onto itself, and `reduce_mesh` reduced
with those 4 only. QE (behaviour from its output, no QE source read) merges a
mesh point with every image that lands on the mesh, using all operations, and
symmetrizes the density with all of them. The same situation occurs in
ordinary inputs: a shifted 4x4x4 mesh in the fcc primitive cell is preserved
by 12 of 48 operations, a 4x4x2 mesh in a cubic cell by 8, any mesh with
unequal divisions along axes that a rotation mixes.

## Theory

Let M be the mesh (N points), G the crystal group (with fractional
translations), P the density projector onto G-invariant fields, and V a
G-invariant potential. Then eps_n(S k) = eps_n(k) and
rho_{S k}(r) = rho_k(S^-1 r), so P rho_{S k} = P rho_k.

1. Point-wise orbits. k ~ k' iff k' = S k (mod G) for some S in G (or -S k
   with time reversal). Because G is a group this is an equivalence on M, and
   the orbit of k is star(k) intersected with M; every member is the image of
   the representative under one operation (the orbit map stays
   `k_full = (+/- R^-T) k_rep + shift`). Weight = orbit size / N.
2. Identity. P(sum_M rho_k / N) = sum_orbits (|orbit| / N) P rho_rep. The
   same holds for band energies, Fermi level and entropy (functions of the
   eigenvalues), and for forces and stress after symmetrization with G
   (F(S k) is the rotated and permuted F(k)). Hence the point-wise reduced run
   with the G projection and G symmetrization is **exactly** the full mesh with
   the G projection and G symmetrization in every SCF map, at every iterate.
3. What it is not. sum_M rho_k is G-invariant only if M is. For a
   non-invariant mesh the unsymmetrized full-mesh density breaks the crystal
   symmetry (Si shifted 4^3: 2.1e-3 e/Bohr^3 under each of the 36
   non-preserving operations), so no reduction with the full group can
   reproduce the unsymmetrized full mesh. The point-wise result is the full
   mesh averaged over the group: the k quadrature (1/|G|) sum_S sum_{k in M}
   delta(k - S k_i) / N, a valid quadrature (positive weights summing to one)
   that converges to the same Brillouin-zone integral, differing from the
   mesh quadrature at the k-sampling error level. The mesh-subgroup reduction
   (H = mesh-preserving subgroup) is exact for the unsymmetrized full mesh and
   keeps only H symmetry in density, forces and stress.
4. Invariant meshes. If every operation (and time reversal, when used)
   preserves M, both policies give the same orbits, representatives,
   operations and weights; the code then returns the former `reduce_mesh`
   result bit for bit.

## Decision

1. Default `K_REDUCTION = pointwise`: k reduction by star(k) on the mesh with
   every crystal operation kept for the run (after ADR 0049 GS1C, the
   operations that permute an explicit FFT grid), plus k ~ -k partners when
   time reversal is on (ADR 0050); images off the mesh are ignored. The
   operation list of the orbit map is the whole group, so the "mesh" density
   group, the force and stress symmetrization and the relaxation projectors
   (all read from that list) use the crystal group. No other code path
   changes; `density_group = "crystal"` now selects the same group.
2. Accepted approximation: on a non-invariant mesh the k quadrature is the
   group average of the mesh (QE-identical sampling; results differ from the
   unsymmetrized full mesh at k-sampling level, measured below). Equality with
   the unsymmetrized full mesh is not a criterion for this policy; its exact
   oracle is `unreduced` (full mesh, crystal projection).
3. `K_REDUCTION = mesh` restores the mesh-preserving subgroup reduction (exact
   for the unsymmetrized full mesh); `unreduced` runs the full mesh with the
   crystal projection and symmetrization (oracle, no k saving).
4. Record: `results.json.symmetry.k_reduction` = {requested, applied
   (`pointwise_crystal`, `mesh_subgroup` for an invariant mesh or the `mesh`
   switch, `unreduced_crystal`), crystal_operations_count,
   mesh_preserving_operations_count, irreducible_k_count, full_k_count,
   policy}. The preflight symmetry line names the applied reduction.
5. Speed of the reduction itself: the mesh lookups use a hash index (wrapped
   coordinates rounded to 1e-6, 27 neighbouring buckets, the same 1e-12
   acceptance test), replacing all-pairs searches of cost N^2 |G|; results are
   unchanged (all former geometry tests pass, mesh-subgroup runs bit-identical).

Why default on although results change on non-invariant meshes: (a) the
policy is exact against its oracle; (b) the former default delivered forces,
stress and densities that break the crystal symmetry on such meshes, and
relaxations then follow only the subgroup; (c) the work-order rule of
identical conditions for QE comparisons (AGENTS.md, 2026-11-14 item 4)
requires the same k quadrature, which QE defines point-wise; (d) the
measured k-sampling error against a dense invariant mesh is equal or smaller
with the point-wise quadrature on every fixture below (forces 2 to 3 times
smaller on the tetrahedron). The orchestrator's stated gate "reduced against
the unsymmetrized full mesh within 1e-10 Ha" cannot hold by item 3 of the
theory for any non-invariant mesh; it holds for invariant meshes (bit
identity) and for `mesh`.

## Validation (Mac, strict path unless noted)

Exact identities, `crates/dft-cli/tests/spin.rs
pointwise_k_reduction_matches_crystal_projected_full_mesh` (release, about
27 s; raw output `tests/reference/symmetry/gt2-pointwise-k-reduction.txt`).
Synthetic nonlocal H potential (l = 1 projector), 6 Ry, explicit grid, Gamma
4x4x2 (32 points), Fermi-Dirac 0.02 Ha, forces (stress for unpolarized FD):

| fixture | ops / mesh-preserving | k: full / pointwise / mesh | pointwise - unreduced: dE, max dF, max dsigma | mesh - full: dE, max dF |
|---|---|---|---|---|
| tetrahedron -43m, cubic 5 Bohr, FD | 24 / 8 | 32 / 9 / 12 | 1.4e-13 Ha, 3.8e-13 Ha/Bohr, 3.1e-14 Ha/Bohr^3 | 1.3e-13, 4.4e-12 |
| same, collinear FD (M = 0.5) | 24 / 8 | 32 / 9 / 12 | 5.4e-14, 3.7e-11 | 1.1e-13, 5.9e-11 |
| rhombohedral -3m, hexagonal axes (Ta6Zn7-like), FD | 12 / 4 | 32 / 11 / 14 | 6.4e-14, 1.4e-13, 6.8e-15 | 1.4e-14, 9.3e-13 |

Size of the quadrature change and k-sampling error (not bounded):

| fixture | pointwise - full: dE, max dF, max dsigma | full - dense 8^3: dE, max dF | pointwise - dense 8^3: dE, max dF |
|---|---|---|---|
| tetrahedron FD | 1.8e-4 Ha, 1.5e-2 Ha/Bohr, 1.2e-3 Ha/Bohr^3 | -2.27e-2, 2.27e-2 | -2.25e-2, 7.3e-3 |
| tetrahedron collinear FD | 2.3e-4, 1.6e-2 | -2.17e-2, 2.28e-2 | -2.15e-2, 6.9e-3 |
| rhombohedral FD | 1.3e-4, 1.3e-2, 6.4e-4 | -9.35e-3, 3.37e-2 | -9.22e-3, 3.34e-2 |

(These toy metals at 4x4x2 are far from k convergence; the numbers show the
direction, not a production error.)

Real pseudopotential, `crates/dft-engine/tests/all/symmetry_full_group.rs
silicon_full_group_oracles` (Si.pz-vbc, 15 Ha, 32^3, tightened tolerances;
raw output `tests/reference/symmetry/gt2-silicon-k444s.txt`): unshifted 4^3
(invariant) pointwise = mesh = crystal, dE = 0, density identical; shifted 4^3
(12 of 48 preserve; 10 k either way) pointwise = mesh reduction with crystal
density group: dE = 0, max density difference 0 (bit-identical), and the
change against the former default (mesh group, = unsymmetrized full mesh):
-1.55e-5 Ha, 1.6e-3 e/Bohr^3.

Geometry (`crates/dft-geometry/tests/symmetry.rs`): Ta6Zn7 6x6x3 point-wise
26 k with and without time reversal (QE 26; union-find orbits identical),
4 mesh-preserving operations; `pointwise_reduction_contract` (diamond meshes:
invariant meshes bit-identical to `reduce_mesh`, orbit maps reconstructed
independently, star(k) closure); `pointwise_smaller_than_mesh_subgroup_fixtures`
(4x4x2: tetrahedron 12 -> 9, diamond 12 -> 8, rhombohedral 14 -> 11).
Existing tests that compare against the unsymmetrized full mesh on a
non-invariant mesh (GT1 `time_reversal_*` on 3x3x1, the G6B diamond
fingerprint, the Si mesh-versus-crystal oracle) now pin `K_REDUCTION=mesh`.

## Consequences

- Ta6Zn7 6x6x3: 34 -> 26 irreducible k (QE count); any cell whose mesh is
  not invariant under the whole group solves fewer k points.
- Results change on non-invariant meshes even when the k count does not (Si
  shifted 4^3: same 10 k, -1.55e-5 Ha) because the density is now projected
  with all operations. Invariant meshes are bit-identical.
- dev08 MP campaign entries with non-invariant meshes change their energies
  at k-sampling level and become QE-identical in k sampling; earlier
  BAANDI-QE differences on such entries partly reflected the different
  quadrature.
