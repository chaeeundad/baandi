TASK_ID: RPA-R4B-WANNIER90
GOAL: Stage R4b of `docs/roadmap/rpa-crpa-plan.md` (decisions D2, D3, D10 of 2026-11-26): Wannier90 3.x external file interface — export of `.win`, `.eig`, `.amn`, `.mmn` (with automatic Marzari-Vanderbilt B1 neighbor shells and optional `.nnkp` neighbor-list matching), Marzari-Vanderbilt spread diagnostics, and import of Wannier90 `_u.mat` / `_u_dis.mat` gauge matrices into `[subspace] method = "wannier90"`.

SCOPE:
1. `dft-response::wannier`:
   - Automatic Marzari-Vanderbilt B1 completeness shell finder `find_bvector_shells(reciprocal_vectors, mesh)` solving `sum_b w_b b_a b_c = delta_ac` to `< 1e-11` for cubic, fcc, bcc, hexagonal, orthorhombic and arbitrary skewed cells.
   - Exact Bloch overlap matrix `bloch_overlap(left, right, g_shift)` computing `M_mn^(k,b) = <u_mk | u_n,k+b> = sum_G c_mk^*(G) c_n,k'(G + dG)` where `dG = round(k_basis + b - k'_basis)` handles arbitrary k-point wrapping cleanly.
   - Marzari-Vanderbilt spread evaluator (`evaluate_spreads`): gauge-invariant spread Omega_I, total spread Omega, and Wannier centers r_bar_n from `M^(k,b)` and gauge matrices `U^(tot)(k)` (Bohr^2 and Bohr).
   - Imported-gauge projection builder `project_with_gauge(eigenvalues, u_total)` verifying semi-unitarity `U^H U = I`, computing band weights `p_mk = sum_i |U_mi|^2` (for R5 cRPA) and `H^(W)(k) = U^H diag(eps) U`.
   - Pure numerical core: no I/O, serde or logging; Hartree atomic units throughout.
2. `dft-engine::response::wannier`:
   - `wannier_export_with` / `wannier_export`: fixed-density NSCF on a full Gamma-centered mesh, raw projection overlaps `A_mn(k)` reusing R4a `build_trial_matrix_at_k`, neighbor list from automatic B1 shells or caller-supplied `.nnkp` neighbor entries, Bloch overlaps `M_mn^(k,b)` for every `(k, b)` pair, adjointness check `max |M^(k',-b) - (M^(k,b))^H|`, and optional initial spread diagnostics.
   - `subspace_from_wannier_with` / `subspace_from_wannier`: load imported `_u.mat` / `_u_dis.mat` gauges on the matching mesh, verify k-coordinates and semi-unitarity, compute subspace weights `p_mk`, populations, Wigner-Seitz `H(R)`, and spreads.
3. `dft-cli`:
   - `baandi wannier-export <input.toml> [--output-dir DIR] [--seedname NAME] [--nnkp FILE]` and `[task] type = "wannier-export"`.
   - Writers (public Wannier90 format spec, no GPL code copied): `<seedname>.win` (Angstrom via `dft_units`), `<seedname>.eig` (eV via `dft_units`), `<seedname>.amn`, `<seedname>.mmn`, `wannier_export.json`, `results.json`.
   - Parsers: `<seedname>.nnkp` (`nnkpts` block), `<seedname>_u.mat`, `<seedname>_u_dis.mat`.
   - Wire `[subspace] method = "wannier90"` (with `wannier_dir` and `seedname` keys) into `baandi subspace`.
4. Tests and records:
   - `crates/dft-response/tests/wannier.rs`: B1 completeness across 5 Bravais lattice families, plane-wave `M_mn` exactness & adjoint symmetry, gauge invariance of Omega_I under random k-dependent unitary rotations, exact center shift under real-space translation, imported gauge checks.
   - `crates/dft-engine/tests/response_r4b.rs`: end-to-end synthetic skewed cell export (A_mn, M_mn adjointness, initial spreads) + gauge import round-trip + MPI rank independence.
   - `crates/dft-cli/tests/wannier_r4b.rs`: end-to-end `baandi wannier-export` + independent parsers of `.win`, `.eig`, `.amn`, `.mmn` + `.nnkp` matching + round-trip `[subspace] method = "wannier90"` gauge import + single entry + guards.
   - ADR `docs/decisions/0060-wannier90-file-interface.md` (Proposed until dev08 Wannier90 3.x oracle runs).
   - Update `tests/reference/response/r4/README.md`, `docs/roadmap/rpa-crpa-plan.md`, `crates/dft-response/AGENTS.md`, `artifacts/task_titles.json`.

OWNED_FILES:
- crates/dft-response/src/lib.rs
- crates/dft-response/src/wannier.rs
- crates/dft-response/tests/wannier.rs
- crates/dft-response/AGENTS.md
- crates/dft-engine/Cargo.toml
- crates/dft-engine/src/response.rs
- crates/dft-engine/src/response/wannier.rs
- crates/dft-engine/tests/response_r4b.rs
- crates/dft-cli/src/config.rs
- crates/dft-cli/src/task.rs
- crates/dft-cli/src/main.rs
- crates/dft-cli/src/subspace.rs
- crates/dft-cli/src/wannier.rs
- crates/dft-cli/tests/wannier_r4b.rs
- docs/decisions/0060-wannier90-file-interface.md
- docs/roadmap/rpa-crpa-plan.md
- tests/reference/response/r4/README.md
- artifacts/task_titles.json
- tickets/RPA-R4B-WANNIER90.md

RULES:
- Follow AGENTS.md, ownership.toml (no new crate dependencies), dft-units (no hardcoded conversion literals).
- Never copy GPL code (Wannier90/QE); file formats follow public Wannier90 documentation only.
- Verify with dev08 official `/scr/mino/gate/gate.sh` (full workspace + mpi,mkl) and `test_no_hardcoded_units.py`.
