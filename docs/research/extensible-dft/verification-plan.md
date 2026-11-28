# Track B Gate 3: End-to-End Verification Plan (Levels 1–5)

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Track B — Phase B1)
- **Phase**: Gate 3 — Pre-Declared Verification Plan for Gate 4 SCF Integration
- **Date**: 2026-11-28
- **Policy**: All tolerances below are pre-declared before Gate 4 SCF integration begins. Post-hoc tolerance relaxation is forbidden (`decision-log.md` Standing Rule 5). Every Level-1 algebraic/variational invariant must be backed by a Lean 4 + Mathlib proof under `formal/BaandiMath/` (`Standing Rule 6`).

---

## 1. Level-1: Formal Mathematical Verification (Lean 4 + Mathlib)

The core Level-1 mathematical guarantees in `formal/BaandiMath/BaandiMath/Solvation.lean` (verified on `dev08` under Lean `v4.36.0` via `bash scripts/check_lean.sh`, mapped hypothesis-by-hypothesis in [formal-verification-map.md](formal-verification-map.md)) remain mandatory CI checks. For Gate 4 SCF integration, `formal/BaandiMath/BaandiMath/Solvation.lean` also formalizes:

1. `spin_polarized_solvation_chain_rule`: $E_{\mathrm{solv}}[n_\uparrow + n_\downarrow]$ yields identical spin-channel potential $\frac{\delta E_{\mathrm{solv}}}{\delta n_\uparrow} = \frac{\delta E_{\mathrm{solv}}}{\delta n_\downarrow} = V_{\mathrm{solv}}[n_\uparrow + n_\downarrow]$ and total double-counting $\langle n_\uparrow, V_{\mathrm{solv}}\rangle + \langle n_\downarrow, V_{\mathrm{solv}}\rangle = \langle n_{\mathrm{tot}}, V_{\mathrm{solv}}\rangle$.
2. `electrostatic_sign_involution`: $U_{\mathrm{el}}(\phi) = -\phi$ is an involution preserving the variational energy sign $\int \rho_{\mathrm{el}}\phi = -\int n\phi = +\int n U_{\mathrm{el}}$.
3. `bregman_bound_nonneg_triangle`: triangle-inequality upper bound combining Hartree, XC, and solvation Bregman remainders in `EnergyAudit::functional_difference_bound`.

- **Pass criterion**: `lake build` exits `0`, zero `sorry`, `#print axioms` restricted strictly to `[propext, Classical.choice, Quot.sound]`.

---

## 2. Level-2: Standalone Solver & Kernel Regression (Already Passing in Gate 2)

Retain all 3 unit tests (`solvation::cavity::tests::*`) and 8 independent integration tests (`crates/dft-physics/tests/solvation_gate2.rs`) plus add:

- `evaluate_with_density_g_matches_evaluate_bit_for_bit`: pre-computed `n_g` matches `evaluate()` to `< 1e-15` Ha across `P = 1, 2, 4` ranks.
- `typed_electrostatic_potential_sign_contract`: unit test in `dft-units` verifying `ElectrostaticPotentialHaPerE(phi).to_electron_energy().0 == -phi` and round-trip involution.

---

## 3. Level-3: Full Self-Consistent Field (SCF) Invariants (Gate 4 Core)

New integration test suite `crates/dft-engine/tests/solvation_scf.rs` (synthetic local/NC pseudopotentials, zero external UPF commits):

| # | Test Case | Driver(s) | Pre-Declared Tolerance / Requirement |
| --- | --- | --- | --- |
| L3-1 | **Off-Path Bit Identity** (`solvation: None`) | Gamma, k-point, spin/FD | **Bit-identical** (`0.0` diff) total energy, eigenvalues, forces, stress, and `results.json` vs `main` baseline |
| L3-2 | **Vacuum-Limit SCF Reduction** ($\epsilon_b = 1, \kappa = 0, \tau = 0, \beta = 0$) | Gamma, k-point, spin/FD | `|E_solv_vac - E_pure_pw| < 1e-12 Ha`, `|F_solv_vac - F_pure_pw| < 1e-10 Ha/Bohr`, `audit.solvation.total_ha == 0.0` |
| L3-3 | **Harris–Foulkes & Band-Form Identity** (polar dimer in water, $\epsilon_b = 78.4, \kappa \in \{0, 0.15\}, \tau > 0$) | Gamma, k-point, collinear spin, Fermi–Dirac | At converged SCF (`dr2_ry < 1e-12`): `|total_hartree - band_form| < 1e-9 Ha`, `|total_hartree - harris_foulkes| < 1e-9 Ha`, `|total_hartree - harris_foulkes| <= functional_difference_bound` |
| L3-4 | **Self-Consistent Hellmann–Feynman Force vs Total-Energy Central FD** | Gamma + k-point (skewed cell, polar heteronuclear dimer, $\epsilon_b = 78.4, \kappa = 0.15, \tau = 2\times 10^{-5}$) | Converge full SCF at $R_I \pm h$ ($h = 5\times 10^{-4}$ Bohr, `dr2_ry < 1e-13`, `solvation_tol = 1e-12`); `max |F_ana - F_fd| < 1e-6 Ha/Bohr` on all 6 axes, exact rigid-translation FD `|(sum_I F_I)_x - F_rigid_fd_x| < 1e-6 Ha/Bohr`, `>10×` grid egg-box suppression (`32^3 -> 48^3`; verified `< 1e-6` on `80^3`), plus `|sum_I F_I| < 1e-6 Ha/Bohr` on the symmetric solvated dimer (`L3-8`) |
| L3-5 | **Ionic Pseudo-Charge Width $r_c$ Invariance in Full SCF** | Neutral & charged ($q = +1, -1$) systems, $r_c \in \{0.45, 0.60, 0.80\}$ Bohr (or well-inside-cavity widths) | `|E_tot(r_c1) - E_tot(r_c2)| < 1e-8 Ha` (neutral & charged with $\kappa > 0$) |
| L3-6 | **Physical Solvation Stabilization & Dielectric Monotonicity** | Polar dimer across $\epsilon_b \in \{1, 4, 20, 78.4\}$ and $\kappa \in \{0, 0.15, 0.30\}$ | Strictly monotonic stabilization $E_{\mathrm{tot}}(\epsilon_b = 78.4) < E_{\mathrm{tot}}(20) < E_{\mathrm{tot}}(4) < E_{\mathrm{tot}}(1)$ and $E_{\mathrm{es}}(\kappa = 0.30) < E_{\mathrm{es}}(0.15) < E_{\mathrm{es}}(0)$ |
| L3-7 | **Collinear Spin-Polarized & Metallic (Fermi–Dirac) SCF** | Open-shell doublet ($N_\uparrow \neq N_\downarrow$) & small smearing cell | Converges cleanly; spin channels receive identical $V_{\mathrm{solv}}[n_\uparrow + n_\downarrow]$; `band_form` & `harris_foulkes` agree `< 1e-9 Ha` |
| L3-8 | **Crystal Symmetry Invariance** | Symmetric dimer/slab with `symmetry` enabled vs disabled | `|E_sym - E_nosym| < 1e-9 Ha`, `max |F_sym - F_nosym| < 1e-8 Ha/Bohr` |
| L3-9 | **MPI Rank Independence in Full SCF** (`dft_mpi::simulate` $P \in \{1, 2, 4\}$) | Gamma & k-point solvated SCF | `|E(P) - E(1)| < 1e-10 Ha`, `max |F(P) - F(1)| < 1e-9 Ha/Bohr` |
| L3-10 | **Fail-Closed Contracts** | `stress = true`, `oracle` path, mismatched `ModelDescription`, invalid TOML `[solvation]` | Explicit structured error (`observables.stress: incomplete derivative coverage: solvation`, etc.); never silently returns zero or vacuum stress |

---

## 4. Level-4: Real-Pseudopotential Physical Benchmark Suite (`dev08`, External UPF Only)

Executed on `dev08` using existing external PseudoDojo NC-SR v0.4 PBE UPFs (never committed to repo):

1. **Isolated $\text{H}_2\text{O}$ molecule in periodic box (`L = 16–20 Bohr`, $E_{\mathrm{cut}} = 60\text{ Ry}$, VASPsol-erfc $\epsilon_b = 78.4$)**:
   - Measure total solvation energy $\Delta E_{\mathrm{solv}} = E_{\mathrm{tot}}(\epsilon_b = 78.4) - E_{\mathrm{tot}}(\text{vac})$ (expected physically in the published VASPsol range `-0.008 to -0.012 Ha` / `-5 to -7.5 kcal/mol` for water in water, Mathew et al. 2014).
   - Measure SCF iteration count (solvated vs vacuum) and verify smooth Broyden convergence without oscillation.
   - Verify full-SCF finite-difference force on O/H atoms matches analytic `ForceResult::values` to `< 5e-6 Ha/Bohr`.
2. **Charged ion / polar molecule with electrolyte screening ($\kappa = 0.10\text{–}0.15\text{ Bohr}^{-1}$)**:
   - Verify $\phi_{\mathrm{bulk}} \to 0$ cleanly away from the solute and box-size convergence vs $L = 16, 20, 24\text{ Bohr}$.
3. **Solver parity on real DFT density (`DirectPcg` vs `FisicaroPcg` at $E_{\mathrm{cut}} = 60\text{ Ry}$)**:
   - Measure exact `|E_tot(FisicaroPcg) - E_tot(DirectPcg)|`, `max |F(FisicaroPcg) - F(DirectPcg)|`, and total SCF wall-time speedup.

---

## 5. Level-5: Workspace Gate & Regression Sign-Off

Before any Gate 5 production merge proposal:

1. `/scr/mino/gate/gate.sh <commit>` on `dev08` passes with `check_rc=0 hard=0 known_platform=0 new=0` across `build`, `test`, `fmt`, `clippy`, `production-check` (`--no-default-features`), `production-clippy`, and `mpi+mkl` `feat_build` + `feat_test`.
2. `lake build` on `formal/BaandiMath/` passes with exit code `0`, zero `sorry`, standard axioms only.
3. Zero-regression benchmark ([benchmark-plan.md](benchmark-plan.md)) confirms pure-PW wall time delta is within measurement noise (`<= 1%`) and bit-identical in output.
