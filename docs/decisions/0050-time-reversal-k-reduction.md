# ADR 0050: Time reversal in the k reduction (unpolarized and collinear)

- Status: **Proposed** (ticket GT1-TIME-REVERSAL, 2026-11-27; orchestrator review pending).
- Scope: k-point reduction by k ~ -k for spin-unpolarized and collinear
  (non-SOC) runs, with and without `[symmetry] mode = "auto"`. Switch
  `BAANDI_COLLINEAR_TIME_REVERSAL` (`[performance] collinear_time_reversal`),
  default `on`; `off` restores the former collinear rejection.

## Context

The dev08 Materials Project campaign (gen61) reported Ta6Zn7 (mp-aaackwro,
13 atoms, 6x6x3 Gamma mesh, PBE, Fermi-Dirac, QE nspin 1) with 34 irreducible
k points in BAANDI against 26 in QE 7.6, both with 12 operations. The preflight
line "magnetic spin swap/time reversal disabled without proof" suggested that
time reversal (TR) was off. It was not:

- For unpolarized runs TR has always been on for `kgrid` input
  (`basis.ktime_reversal` defaults to true for a grid, false for an explicit
  list; ADR 0013). With `symmetry.mode = "auto"`, `reduce_mesh` adds one
  antiunitary partner (k -> -R^-T k) per mesh-preserving operation. The dev08
  results.json lists these partners (`antiunitary: true`). The preflight text
  was a fixed string, printed whatever the setting.
- Ta6Zn7 is centrosymmetric (QE: "12 Sym. Ops., with inversion"), so TR adds
  nothing there.
- The real cause is the mesh-preserving subgroup. In the Materials Project
  cell setting only 4 of the 12 operations map the 6x6x3 mesh onto itself,
  and `reduce_mesh` reduces with those 4 only (34 k). If all 12 operations are
  applied point by point, keeping every image that lands on a mesh point (the
  usual practice, and what QE's count matches), the result is 26 orbits. The
  geometry test `ta6zn7_mesh_subgroup_explains_k_count` reproduces 12 / 4 /
  34 / 26 on the campaign input. Changing this is outside this ADR (see
  "Follow-up").

What was really missing was TR for **collinear** spin. Collinear input had to
set `basis.ktime_reversal = false` (CLI rejection, plus engine guards in the
collinear oracle and the reference spin/FD driver), so every magnetic run used
twice the k points of a non-centrosymmetric unpolarized run.

## Decision

1. **Unpolarized (fixed and Fermi-Dirac):** unchanged. TR is on by default for
   `kgrid`, explicit `basis.ktime_reversal` wins.
2. **Collinear (fixed sector, fixed moment, free-moment FD):** TR within each
   spin channel is accepted and is the default for `kgrid`, as for unpolarized
   runs. Each channel's Kohn-Sham Hamiltonian H_s = T + V_loc + V_H + V_xc,s +
   V_NL has a real local potential and real, spin-independent nonlocal
   coefficients (norm-conserving, scalar-relativistic, no SOC, no field acting on
   orbitals, no orbital-dependent Hubbard term, which the CLI rejects). Hence
   H_s(-k) = conj(H_s(k)) for each s separately: eps_s,n(-k) = eps_s,n(k),
   psi_s,n,-k = conj(psi_s,n,k), and rho_s,-k(r) = rho_s,k(r). The k ~ -k
   operation acts inside each channel and never exchanges up and down; it is
   not the magnetic time reversal of a spin-full model (that one flips the spin
   and stays disabled, as do spin swaps). Uniform fixed total moment and channel
   counts are invariant under it. Spatial operations stay identity-spin and the
   start magnetisation is uniform, so crystal operations remain valid.
3. **Density:** no extra projection is needed. The density from a weighted
   irreducible set equals the full-mesh density because rho_-k = rho_k per
   channel; with crystal operations the unitary projector is applied as before
   (antiunitary partners are excluded from the projector, as for unpolarized
   runs). Forces and stress are TR-even and need no change.
4. **Record:** results carry `time_reversal: {used, source, scope, policy}`
   (fixed, unpolarized FD, collinear and reference spin/FD writers); the
   preflight line now states the crystal/mesh operation counts and whether TR
   is on, and within each channel for collinear runs.
5. **Switch:** `COLLINEAR_TIME_REVERSAL` (default `on`; table 87 -> 88).
   `off` restores the former rejection in the CLI and both engine guards.
6. **Inversion:** when inversion is in the mesh-preserving group, the
   irreducible points and weights with and without TR are identical (geometry
   test `time_reversal_folds_only_without_inversion`), so such runs are
   unchanged bit for bit.

## Validation (Mac, debug-optimized test profile)

Fixture: four synthetic Z=1/2 atoms on a regular tetrahedron in a 5 Bohr cubic
cell (point group -43m, no inversion), nonlocal synthetic potential with an
l = 1 projector, 6 Ry, explicit 18^3 grid, 3x3x1 Gamma mesh. Oracle: full mesh
(9 k), symmetry off, no TR. Arms: TR only (5 k), crystal + TR (3 k).
Strict path (`BAANDI_REFERENCE_STRICT=on`), |dE| (free energy for FD) and
max |dF|:

| model | TR only | crystal + TR |
|---|---|---|
| unpolarized fixed | 0, 1.6e-8 Ha/Bohr | 1.8e-15, 1.5e-9 |
| unpolarized FD | 1.7e-13, 2.6e-12 | 9.3e-14, 1.3e-12 |
| collinear fixed (1,1) | 5.8e-15, 2.8e-9 | 2.2e-15, 1.9e-9 |
| collinear FD (M = 0.5) | 8.1e-14, 7.5e-11 | 3.5e-14, 3.8e-11 |

Standard reference spin/FD driver (threshold tightened to 1e-12 Ha): collinear
fixed 6.4e-14 / 1.3e-12 Ha, collinear FD 4.3e-13 / 5.0e-14 Ha, forces within
2.1e-7 Ha/Bohr. Energies are inside the strict class 1e-10 Ha; forces inside
the 1e-6 Ha/Bohr bound of the existing path comparisons (the strict SCF leaves
about 1e-8 Ha/Bohr of force noise). Tests: `dft-cli --test spin
time_reversal_matches_full_mesh`, `time_reversal_standard_spin_driver`;
`dft-geometry --test symmetry time_reversal_folds_only_without_inversion`
(3x3x3 tetrahedron: 27 full, 14 TR, 5 crystal, 4 crystal + TR).

## Consequences

- Collinear inputs with a k grid no longer need `basis.ktime_reversal = false`;
  existing inputs that set it keep the full mesh and stay bit-identical.
- Non-centrosymmetric magnetic cells solve up to half the k points.
- Ta6Zn7-type gaps against QE are not closed by this ADR.

## Follow-up (not decided here)

Point-wise reduction with the full crystal group when the mesh is not
invariant under every operation (Ta6Zn7: 34 -> 26 k). This is exact only if
the density is projected with the full crystal group (`symmetry.density_group
= "crystal"`), because the full-mesh sum of a non-invariant mesh is not
symmetric; it therefore changes the converged density at the k-sampling error
level against the present default (`mesh`), and forces/stress must be
symmetrized with the same group. It needs its own ticket and ADR.

**Resolved by ADR 0051 (GT2, 2026-11-27):** point-wise reduction with the full
crystal group is the default (`K_REDUCTION=pointwise`); Ta6Zn7 6x6x3 now
reduces to 26 k.
