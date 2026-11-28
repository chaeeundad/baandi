TASK_ID: RPA-R4A-ATOMIC-PROJECTION
GOAL: Stage R4a of `docs/roadmap/rpa-crpa-plan.md` (decisions D1-D10 of 2026-11-26): atomic-orbital projection with Loewdin orthonormalization, band window selection, subspace weights p_nk (for R5 cRPA), projected populations, real-space tight-binding hoppings H(R) with Wigner-Seitz degeneracies, and Wannier band interpolation.

SCOPE:
1. `dft-response::subspace`:
   - Projection matrix A_mn(k) = <psi_mk | g_nk> on a band window W(k), overlap S(k) = A^H A, Loewdin orthonormalization U(k) = A(k) S(k)^{-1/2} via Hermitian EVD, strict singularity rejection when min eigenvalue <= tol.
   - Orthonormality defect max |U^H U - I|, spilling parameter, per-band subspace weights p_mk = sum_i |U_mi(k)|^2 (sum rule: sum_{m in W} p_mk = N_orb), projected Hamiltonian H^(W)(k) = U^H diag(eps) U.
   - Exact Wigner-Seitz supercell lattice points and degeneracies for arbitrary (including skewed) cells, verified by sum_R 1/g(R) = N_k.
   - Real-space hoppings H(R) = (1/N_k) sum_k exp(-2 pi i k.R) H^(W)(k), Hermiticity check H_ij(R) = H_ji^*(-R), and Wannier interpolation H^(W)(k') and eigenvalues at arbitrary k'.
   - Pure numerical core: no I/O, no serde, no logging; Hartree atomic units throughout.
2. `dft-engine::response::subspace`:
   - Build atomic trial orbitals in the plane-wave basis from UPF PSWFC radial wavefunctions (normalized on the radial mesh) and `dft_upf::radial::real_spherical_harmonic` (matching `refstart.rs` Fourier and phase conventions).
   - Support orbital selectors by species, optional atom subset, angular momentum l, optional radial channel label, orbital names (s, p, px/py/pz, d, dz2/dxz/dyz/dx2-y2/dxy, t2g, eg, f, and tetrahedral sp3 hybrids).
   - Band window selection by energy window `window_ha` (relative to VBM/reference) or explicit 1-based inclusive `band_range`.
   - Deterministic execution under serial and MPI communicators (np = 1 bit-identical to serial; multi-rank agrees to roundoff).
   - Return raw (all-band and in-window) and Loewdin projected populations, per-k projections, subspace weights p_nk, and real-space Hamiltonian H(R).
3. `dft-cli`:
   - `[subspace]` TOML table (`method = "atomic-projection"`, `k_mesh`, `bands`, `band_solver`, `eigen_tolerance_ha`, `window_ha` / `band_range`, shorthand `species`/`l`/`orbitals`/`atoms`/`channel` or multi-entry `[[subspace.projections]]`).
   - `baandi subspace <input.toml> [--output-dir DIR]` subcommand and `[task] type = "subspace"` single-entry support.
   - Writers: `subspace.json` (with deterministic SHA-256 `projector_fingerprint`), `subspace_hr.dat` (Wannier90 `seedname_hr.dat` layout in eV via `dft_units`), `subspace_populations.csv`, `subspace_weights.csv`, optional `subspace_bands.csv` (when `[bands]` path is present), `results.json`.
4. Tests and records:
   - `crates/dft-response/tests/subspace.rs`: orthonormality <= 1e-12, isolated-manifold exact eigenvalue recovery, entangled-window sum rule & Loewdin optimality, WS partition-of-unity on cubic/fcc/hexagonal/skewed cells, H(R) Hermiticity, exact interpolation on mesh k-points, singular-overlap rejection.
   - `crates/dft-engine/tests/response_r4a.rs`: synthetic skewed cell with atomic UPF, isolated & entangled windows, sp3 hybrids, multi-rank determinism, population sum rules.
   - `crates/dft-cli/tests/subspace_r4a.rs`: end-to-end CLI & single-entry test, independent parser check of `subspace_hr.dat` against the public Wannier90 format specification, CSV/JSON consistency, bad-input guards.
   - ADR `docs/decisions/0059-subspace-atomic-projection.md` (Proposed until dev08 `projwfc.x`/Wannier90 oracle runs).
   - Oracle protocol `tests/reference/response/r4/README.md`, `docs/roadmap/rpa-crpa-plan.md` status update, `crates/dft-response/AGENTS.md` update, `artifacts/task_titles.json`.

OWNED_FILES:
- crates/dft-response/src/lib.rs
- crates/dft-response/src/subspace.rs
- crates/dft-response/tests/subspace.rs
- crates/dft-response/AGENTS.md
- crates/dft-engine/src/response.rs
- crates/dft-engine/src/response/subspace.rs
- crates/dft-engine/tests/response_r4a.rs
- crates/dft-cli/src/config.rs
- crates/dft-cli/src/task.rs
- crates/dft-cli/src/main.rs
- crates/dft-cli/src/subspace.rs
- crates/dft-cli/tests/subspace_r4a.rs
- docs/decisions/0059-subspace-atomic-projection.md
- docs/roadmap/rpa-crpa-plan.md
- tests/reference/response/r4/README.md
- artifacts/task_titles.json
- tickets/RPA-R4A-ATOMIC-PROJECTION.md

RULES:
- Follow AGENTS.md, ownership.toml (no new crate dependencies), dft-units (no hardcoded conversion literals).
- Never copy GPL code (Wannier90/QE); format writers follow public documentation only.
- Verify with `scripts/check.sh full` and `python3 -I scripts/tests/test_no_hardcoded_units.py`.
