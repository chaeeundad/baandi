# ADR 0060: Wannier90 file export (.win/.amn/.mmn/.eig) and gauge import (R4b)

- Status: **Proposed** (ticket `RPA-R4B-WANNIER90`, 2026-11-28; internal numerical, engine and CLI gates verified locally/dev08; external Wannier90 3.x vs QE `pw2wannier90.x` oracle runs on dev08 pending).
- Scope: `dft_response::wannier`, `dft_engine::response::wannier`, `baandi wannier-export` / `[task] type = "wannier-export"`, `[subspace] method = "wannier90"`.
- Relates to: `docs/roadmap/rpa-crpa-plan.md` Sections 2.8, 2.11, 3.6 and Stage R4b (user decision D2 of 2026-11-26: Wannier90 only as an external program through files, respecting the MIT/GPL boundary); ADR 0017 (writers in CLI only); ADR 0059 (R4a atomic-orbital projection).

## Context

Decision D2 establishes that maximally localized Wannier functions (MLWF) and Souza-Marzari-Vanderbilt disentanglement are performed by running Wannier90 3.x externally through its standard text files, never by linking GPL code into BAANDI. R4a (ADR 0059) already constructs the atomic trial orbitals, Loewdin projections, subspace weights `p_mk` and Wigner-Seitz real-space Hamiltonian `H(R)`. R4b adds:

1. Export of `<seedname>.win`, `<seedname>.eig`, `<seedname>.amn` and `<seedname>.mmn`.
2. Support for matching an external `wannier90.x -pp` neighbor list via `--nnkp <file>` as well as standalone Marzari-Vanderbilt B1 shell generation.
3. Import of Wannier90's `<seedname>_u.mat` (and optional `<seedname>_u_dis.mat`) into `baandi subspace` (`method = "wannier90"`), producing the same typed `p_mk` weights and `H(R)` representation consumed by R5 (cRPA).

## Decision

### 1. Bloch overlap convention `M_mn^(k,b)` (`.mmn`)

For cell-normalized Bloch orbitals `psi_{nk}(r) = Omega^{-1/2} sum_G c_{nk}(G) exp(i(k_basis + G).r)` (where `nscf` wraps `k_mesh in [0, 1)^3` to `k_basis in [-1/2, 1/2)^3`) and a neighbor relation `k_mesh + b = k'_mesh + G0` (`G0 in Z^3`), the periodic cell parts satisfy `u_{n,k+b}(r) = exp(-i G0.r) u_{n,k'}(r)`, giving

```text
M_mn^(k,b) = <u_mk | u_{n,k+b}> = sum_G c_{mk}^*(G) c_{n,k'}(G + dG_basis)
```

with exact integer basis shift `dG_basis = round(k_basis + b - k'_basis)`. This satisfies the exact adjoint identity `M_nm^(k', -b) = (M_mn^(k,b))^*` (verified on every export and reported as `max_adjoint_defect`).

### 2. Neighbor b-vector shells (Marzari-Vanderbilt B1) and `.nnkp`

- **Standalone mode** (no `--nnkp`): `dft_response::wannier::find_bvector_shells` enumerates reciprocal sub-mesh lattice vectors `b = n_0 (b_0/N_0) + n_1 (b_1/N_1) + n_2 (b_2/N_2)` in ascending `|b|^2` shells and solves the B1 completeness linear system `sum_s w_s (sum_{b in s} b_a b_c) = delta_ac` (6 symmetric equations) by adding linearly independent shells until residual `< 1e-11`. Works for cubic, fcc, bcc, hexagonal, orthorhombic and arbitrary triclinic cells.
- **External `.nnkp` mode** (`--nnkp <seedname.nnkp>`): parses the `begin nnkpts` / `end nnkpts` block written by `wannier90.x -pp` so `<seedname>.mmn` matches Wannier90's exact neighbor ordering and shifts `G0`.

### 3. Marzari-Vanderbilt spread diagnostic

Given neighbor weights `w_b`, overlaps `M^(k,b)` and semi-unitary gauge matrices `U^(tot)(k)` (`N_bands x N_wann`):
- Transformed overlaps `M^(W)(k,b) = U^(tot)(k)^H M^(k,b) U^(tot)(k')`.
- Wannier centers: `r_n = -(1 / N_k) sum_{k,b} w_b b Im ln M^(W)(k,b)_{nn}` (Bohr).
- Gauge-invariant spread: `Omega_I = (1 / N_k) sum_{k,b} w_b (N_wann - sum_{m,n} |M^(W)(k,b)_{mn}|^2)` (Bohr^2).
- Total spread: `Omega = (1 / N_k) sum_{k,b} w_b sum_n [(1 - |M^(W)_{nn}|^2) + (Im ln M^(W)_{nn} + b . r_n)^2]` (Bohr^2).

Recorded in `wannier_export.json` and `subspace.json` in both Bohr^2 and Angstrom^2 (via `dft_units`).

### 4. File formats and gauge import

Writers and parsers live exclusively in `dft-cli` (ADR 0017) and follow the public Wannier90 file format specification:
- `<seedname>.win`: lattice in Angstrom (`dft_units`), fractional atom coordinates, `mp_grid`, `num_bands`, `num_wann`, `kpoints` block, optional `dis_win_min`/`dis_win_max` (eV), projections block.
- `<seedname>.eig`: `m k E_mk` (1-based indices, eV via `dft_units`).
- `<seedname>.amn`: `m n k Re(A_mn(k)) Im(A_mn(k))` (`m` fastest, then `n`, then `k`), using the exact R4a `build_trial_matrix_at_k` overlaps.
- `<seedname>.mmn`: per `(ik, ik_nn, G0)` block, `Re(M_mn) Im(M_mn)` with `m` fastest then `n` (column-major).
- Import (`[subspace] method = "wannier90"`, `wannier_dir`, `seedname`): parses `<seedname>_u.mat` and optional `<seedname>_u_dis.mat`, checks k-point coordinate agreement (`<= 1e-6`) and semi-unitarity `max |U^H U - I| <= 1e-8`, forms `U^(tot)(k) = U_dis(k) U(k)` (or `U(k)` when square), and produces the same `subspace.json`, `subspace_hr.dat`, `subspace_weights.csv` (`p_mk`), `subspace_populations.csv` and optional `subspace_bands.csv` outputs as R4a.

## Consequences

- Complete MIT-clean interoperability with Wannier90 3.x without linking GPL libraries.
- R5 (cRPA) can seamlessly accept either R4a projected subspaces or R4b imported MLWF subspaces through one unified gauge/weight representation.
- External `dev08` validation comparing `.amn`/`.mmn`, spreads and `hr.dat` against QE `pw2wannier90.x` + Wannier90 3.x is tracked in `tests/reference/response/r4/README.md`.
