# ADR 0059: Atomic-orbital subspace projection, Loewdin orthonormalization and real-space hoppings (R4a)

- Status: **Proposed** (ticket `RPA-R4A-ATOMIC-PROJECTION`, 2026-11-28; internal numerical, engine and CLI gates pass locally; dev08 oracle checks against QE `projwfc.x` and Wannier90 3.x pending).
- Scope: `dft_response::subspace`, `dft_engine::response::subspace`, `baandi subspace` / `[task] type = "subspace"` / `[subspace]`.
- Relates to: `docs/roadmap/rpa-crpa-plan.md` Sections 2.8, 2.11, 3.1, 3.6 and Stage R4 (user decisions D1-D10 of 2026-11-26, option A each); ADR 0017 (writers in CLI only); ADR 0052 (NSCF bands); ADR 0054 (response conventions).

## Context

Stages R1-R3 established many-band NSCF solves, dielectric matrices and ACFDT RPA correlation energies. Stage R4 constructs a localized correlated subspace for cRPA (R5) and future DFT+U consumption (R6), as well as tight-binding models and Wannier interpolation. Plan Section 2.8 splits R4 into:

- **R4a**: atomic-orbital projection using UPF `PSWFC` radial wavefunctions, projected onto a band window and Loewdin-orthonormalized ("projected Wannier" functions, no gauge minimization), together with real-space hoppings `H(R)`. This is always defined for isolated and entangled manifolds and matches the ortho-atomic projector convention needed by R6 (DFT+U).
- **R4b**: external Wannier90 3.x file export (`.win`, `.amn`, `.mmn`, `.eig`) and gauge import (`_u.mat`, `_u_dis.mat`).

## Decision

### 1. Trial atomic orbitals in the plane-wave basis

For an atom at Cartesian position `tau` (Bohr) with UPF radial pseudo-wavefunction `chi(r) = r R(r)` (angular momentum `l <= 3`) normalized on the UPF radial mesh (`int |chi(r)|^2 dr = 1`), the cell-normalized plane-wave representation at wavevector `k + G` is

```text
g_{tau,l,m}(k + G) = (-i)^l exp(-i (k + G) . tau) * Omega^{-1/2}
                     * [4 pi int_0^inf r chi(r) j_l(|k + G| r) dr] * Y_{lm}(hat{k + G})
```

where `Y_{lm}` is the real spherical harmonic of `dft_upf::radial::real_spherical_harmonic` (Condon-Shortley `P_l^m` with `(-1)^m` real combination; at `|k + G| = 0`, `Y_00 = 1/sqrt(4 pi)` for `l = 0` and `0` for `l > 0`). This matches `dft-engine::refstart` bit-for-bit in Fourier and phase conventions.

Standard real-harmonic names map to `(l, m)` as:
- `l = 0`: `s` (`m = 0`)
- `l = 1`: `py` (`m = -1`), `pz` (`m = 0`), `px` (`m = +1`); shorthand `p` expands to `["pz", "px", "py"]` (Wannier90 standard order)
- `l = 2`: `dxy` (`m = -2`), `dyz` (`m = -1`), `dz2` (`m = 0`), `dxz` (`m = +1`), `dx2-y2` (`m = +2`); shorthand `d` expands to `["dz2", "dxz", "dyz", "dx2-y2", "dxy"]`, `t2g` to `["dxy", "dyz", "dxz"]`, `eg` to `["dz2", "dx2-y2"]`
- `l = 3`: `f-3`..`f+3` (`m = -3..=+3`); shorthand `f` expands to `m = 0, +1, -1, +2, -2, +3, -3` (Wannier90 standard order)
- Tetrahedral hybrids: `sp3` expands on each selected atom (requiring both `l = 0` and `l = 1` PSWFC channels) to the four normalized hybrids `(s + px + py + pz)/2`, `(s + px - py - pz)/2`, `(s - px + py - pz)/2`, `(s - px - py + pz)/2`.

When a pseudopotential has multiple `PSWFC` channels of the same `l`, an explicit `channel` label selects one; if omitted, the channel with the highest occupation (breaking ties by the last channel, i.e. valence over semicore) is chosen and recorded in `subspace.json`. Trial orbitals are NOT renormalized per k-point in the finite plane-wave sphere, so raw overlaps and spilling measure plane-wave and window completeness honestly.

### 2. Band window, projection and Loewdin orthonormalization

At each k-point on a Gamma-centered uniform `N_1 x N_2 x N_3` mesh, a band window `W(k)` is selected either by:
- `window_ha = [e_min, e_max]` relative to the reference energy (insulator VBM across the mesh), or
- `band_range = [first, last]` (1-based inclusive global band range).

Every k-point must satisfy `|W(k)| >= N_orb`, otherwise the run rejects with an explicit error.

With the Bloch eigenvectors `C_W(k)` (`N_pw x |W(k)|`) and trial orbitals `G(k)` (`N_pw x N_orb`):
1. Raw overlap: `A(k) = C_W(k)^H G(k)` (`|W(k)| x N_orb`).
2. Overlap metric: `S(k) = A(k)^H A(k)` (Hermitized explicitly).
3. Hermitian EVD `S(k) = V diag(lambda) V^H`. If `min(lambda) <= singular_tol` (default `1e-10`), reject; never clip or pseudo-invert a rank-deficient projection silently.
4. Loewdin gauge: `U(k) = A(k) V diag(lambda^{-1/2}) V^H` (`|W(k)| x N_orb`), satisfying `U(k)^H U(k) = I` to machine precision (gate: `<= 1e-12`).
5. Per-band subspace weights (for R5 cRPA projector exclusion, plan 2.9):
   `p_{mk} = sum_{i=1}^{N_orb} |U_{mi}(k)|^2` for `m in W(k)` and `0` outside, obeying `sum_m p_{mk} = N_orb`.
6. Projected Hamiltonian: `H^(W)(k) = U(k)^H diag(eps_W(k)) U(k)` (Hermitized explicitly). For an isolated band manifold (`|W(k)| = N_orb`), `U(k)` is unitary and the eigenvalues of `H^(W)(k)` reproduce the Kohn-Sham eigenvalues in `W(k)` to roundoff.
7. Spilling parameter: `s(k) = 1 - (1 / N_orb) Tr S(k)`.
8. Projected populations per orbital (occupied bands, spin factor 2 for unpolarized): both raw (all solved bands and in-window) and Loewdin-orthonormalized in-window populations are recorded.

### 3. Wigner-Seitz supercell and real-space hoppings `H(R)`

For the `N_1 x N_2 x N_3` Gamma-centered mesh and primitive lattice vectors `a_0, a_1, a_2` (columns of `A`), the supercell vectors are `A_i = N_i a_i`. The Wigner-Seitz lattice points `R in Z^3` satisfy `|A R|^2 <= min_T |A (R - T)|^2 + tol` over all supercell lattice shifts `T = (t_0 N_0, t_1 N_1, t_2 N_2)`, with search bounds derived analytically from the longest supercell half-diagonal and the reciprocal basis norms so that arbitrarily skewed cells are exact. The degeneracy `g(R)` is the multiplicity of achieving the minimum; the exact partition-of-unity sum rule `sum_R 1/g(R) = N_1 N_2 N_3` is verified on every call.

Real-space hoppings (Hartree internally):
```text
H_{ij}(R) = (1 / N_k) sum_k exp(-2 pi i k . R) H^(W)_{ij}(k)
```
and Wannier interpolation at arbitrary `k'`:
```text
H^(W)_{ij}(k') = sum_{R in WS} (1 / g(R)) exp(2 pi i k' . R) H_{ij}(R)
```
which satisfies `H_{ij}(R) = H_{ji}^*(-R)` and reproduces `H^(W)(k)` identically on every mesh k-point.

### 4. Projector fingerprint and CLI outputs

`baandi subspace <input.toml>` (and `[task] type = "subspace"`) writes:
- `subspace.json` (and `results.json`): metadata, units, orbital table, deterministic SHA-256 `projector_fingerprint` (UPF hashes, species/atom/channel/orbital definitions, real-harmonic convention, orthonormalization, window, k-mesh), per-k overlap eigenvalues/orthogonality defects/spilling, populations, and real-space hopping summary.
- `subspace_hr.dat`: standard Wannier90 `seedname_hr.dat` layout (1-based orbital indices, WS degeneracies 15 per line, hoppings converted to eV via `dft_units`).
- `subspace_populations.csv`: per-orbital raw and Loewdin populations.
- `subspace_weights.csv`: per-k, per-band energies (Ha and eV), window flag and subspace weights `p_{mk}`.
- `subspace_bands.csv` (optional, when `[bands]` path is supplied): Wannier-interpolated bands along the path.

## Consequences

- `dft-response` gains `subspace` with no new crate dependencies (`ownership.toml` unchanged).
- R4b (Wannier90 `.win/.amn/.mmn/.eig` export and `_u.mat`/`_u_dis.mat` import) can reuse the exact trial-orbital builder, window selector, WS grid and `H(R)` representation introduced here.
- R5 (cRPA) directly consumes `KProjection::band_weights` (`p_{mk}`) and `KProjection::gauge` (`U(k)`).
- Collinear spin and metallic (Fermi-Dirac) ground states remain deferred until R5b / post-R5 (decision D6); spin-polarized inputs reject explicitly in R4a.
