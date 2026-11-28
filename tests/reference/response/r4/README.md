# R4 oracle protocol: atomic projection (R4a) and Wannier90 subspace (R4b)

Reference plan: `docs/roadmap/rpa-crpa-plan.md` Section 4 (R4) and ADR 0059 (R4a).
No pseudopotentials are stored in this repository. External oracle runs execute on
`dev08` at identical conditions (cell, UPF file, cutoff, FFT grid, k mesh, band
window and orbital selectors).

## R4a: Atomic-orbital projection with Loewdin orthonormalization

### Cases (generality set, plan Section 4)
1. **SrVO3 (cubic perovskite, V $t_{2g}$ isolated manifold)**: `species = "V"`, `orbitals = ["t2g"]`, `k_mesh = [4, 4, 4]` and `[6, 6, 6]`.
2. **Si (diamond, $sp^3$ entangled manifold)**: `species = "Si"`, `orbitals = ["sp3"]`, `window_ha` spanning valence + lower conduction bands, `k_mesh = [4, 4, 4]`.
3. **Diamond C ($sp^3$)**: `species = "C"`, `orbitals = ["sp3"]`, `k_mesh = [4, 4, 4]`.
4. **Wurtzite AlN (hexagonal, N $p$ / Al $s,p$)**: non-cubic WS degeneracy and anisotropic hoppings.

### Oracle comparison (QE `projwfc.x` 7.6 on dev08)
1. Run `pw.x` SCF and NSCF on the identical Gamma-centered full k-mesh (`nosym = .true.`, `noinv = .true.`).
2. Run `projwfc.x` on the same prefix to obtain Loewdin-orthonormalized orbital populations and spilling parameter.
3. Run `baandi subspace <case>.toml` with the same UPF, cutoff, k-mesh and orbital selection.

### Acceptance gates (plan Section 4)
- Orthonormality defect `max |U^H U - I| <= 1e-12` at every k-point (internal).
- Real-space hopping Hermiticity `max |H_ij(R) - H_ji^*(-R)| <= 1e-12` Ha (internal).
- Isolated-manifold eigenvalue recovery on every mesh k-point `<= 1e-12` Ha (internal).
- Projected orbital populations versus QE `projwfc.x` within `1e-3` electrons.
- Real-space hoppings `subspace_hr.dat` vs Wannier90 (with `num_iter = 0`, projection-only) within `1 meV`.

## Status
- **2026-11-28 (R4a)**: Numerical core (`dft_response::subspace`), engine driver (`dft_engine::response::subspace`), CLI (`baandi subspace`, `[task] type = "subspace"`) and synthetic/analytic gates implemented and verified locally on Mac. External `dev08` runs against QE `projwfc.x` and Wannier90 3.x pending.
- **R4b** (`.win/.amn/.mmn/.eig` export and `_u.mat`/`_u_dis.mat` import): next ticket.
