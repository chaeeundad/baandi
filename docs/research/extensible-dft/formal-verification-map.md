# Track B Formal-to-Rust Assumption Mapping (`BaandiMath.Solvation`)

- **Scope**: Explicit hypothesis-by-hypothesis traceability between the Lean 4 + Mathlib formalization in [`formal/BaandiMath/BaandiMath/Solvation.lean`](../../../formal/BaandiMath/BaandiMath/Solvation.lean) and the Rust implementation across `dft-units`, `dft-physics::solvation`, and `dft-engine`.
- **Policy**: Per `AGENTS.md` (Formal mathematical verification policy) and [`decision-log.md`](decision-log.md) (`D-014`), formal verification is judged strictly by the physical and algorithmic guarantees it establishes—never by theorem counts.
- **Reproducible CI audit**: `bash scripts/check_lean.sh` (verifies zero `sorry`, zero custom `axiom` declarations, clean `lake build`, and `#print axioms ⊆ {propext, Classical.choice, Quot.sound}`).

---

## 1. Guarantee-by-Guarantee Assumption Mapping

### 1.1 Nyquist-Zeroed Spectral Derivative & Realistic PCG-Space SPD

| Theorem (`BaandiMath.Solvation`) | Formal Hypotheses | Rust Enforcement Site (File, Function, Lines) | Verification Notes |
| --- | --- | --- | --- |
| `nyquist_zeroed_multiplier_skew_adjoint`, `nyquist_nonzero_breaks_reality` | `mirror` is `-G mod N`; `hg_odd: g (mirror k) = -g k` (forcing `g k = 0` when `mirror k = k`); `(ur, ui)` conjugate-symmetric | `crates/dft-physics/src/solvation/cavity.rs`: `is_even_nyquist` (`L21-31`), `skew_adjoint_g_cart` (`L36-59`); `crates/dft-physics/src/solvation/poisson.rs`: `FftEngine::new` (`L88-104`, `L121-147` mirror index table), `FftEngine::apply_operator` (`L266-439`, `-i G_a` adjoint pass at `L344-349`, `L423-436`) | Even-axis Nyquist Miller components (`2 * ijk[a] == dims[a]`) are zeroed *before* multiplying by `recip_cols`, guaranteeing real-field closure and exact `D_a^T = -D_a` on skewed and cubic cells. |
| `gpOp_add`, `gpOp_quad_form`, `gpOp_symmetric`, `gpOp_psd` | `dv > 0`, `hskew: ⟨u, D_a v⟩ = -⟨D_a u, v⟩`, `∀ i, 1 ≤ eps i`, `∀ i, 0 ≤ k2 i` | `poisson.rs`: `FftEngine::new` (`volume > 0`, `dv = volume / total_points`, `L60-65`, `L156`); `PdeSolver::solve_pde` collective validation rejecting `eps < 1.0` and `kappa2 < 0.0` (`L796-841`); `cavity.rs`: `evaluate_cavity_fields` (`L90-154`); `mod.rs`: `SolvationParams::validate` (`L175-231`) | `FftEngine::apply_operator` (`L266-439`) computes `∑_a (-i G_a) F[eps * F^{-1}(i G_a x)] + F[kappa2 * x]`, matching `gpOp` algebraically under Parseval's identity. |
| `gpOp_quad_eq_zero_iff` | Same as `gpOp_psd` (`dv > 0`, `eps ≥ 1`, `k2 ≥ 0`) | Same enforcement as above | Exact characterization of `ker(A)`: `⟨u, A u⟩ = 0` iff `u ∈ ker(D)` and `u = 0` on `{i \| k2 i > 0}`. |
| `gpOp_spd_on_grad_complement` | `dv > 0`, `eps ≥ 1`, `k2 ≥ 0` (allows `k2 ≡ 0`), `∃ a i, D a u i ≠ 0` | `poisson.rs`: `FftEngine::new` collects `zero_g2_indices` (`norm2 == 0.0`, covering `G = 0` and all even-Nyquist corner modes where `D_a = 0` for all `a`, `L92-104`); `PdeSolver::solve_pde` projects `rhs_eff_g` onto `zero_g2_indices = 0` when `!has_kappa` (`L850-855`); `fill_direct_precond` sets `precond = 0` at `g2 == 0` when `!has_kappa` (`L939-950`); `solve_direct_pcg` & `solve_direct_pcg_from_guess` zero `out_phi_g` and `w_g` on `zero_g2_indices` (`L992-996`, `L1040-1044`, `L1139-1143`) | Proves strict SPD of `DirectPcg` on the exact pure-GP (`κ = 0`) computation subspace. Every nonzero PCG search direction `p_g` has nonzero energy outside `zero_g2_indices`, hence `∃ a i, D a p i ≠ 0`. |
| `gpOp_spd_lpb_cavity_vanishing` | `dv > 0`, `eps ≥ 1`, `k2 ≥ 0` (explicitly allowing `k2 = 0` inside solute cavity), `(∃ a i, D a u i ≠ 0) ∨ (∃ i, 0 < k2 i ∧ u i ≠ 0)` | `poisson.rs`: `PdeSolver::solve_pde` sets `has_kappa = sum(kappa2) > 0.0` (`L843-844`), keeps full space unprojected, and uses strictly positive diagonal preconditioner `1 / (eps_bar * g2 + kappa2_bar)` (`fill_direct_precond`, `L943-949`); `solve_direct_pcg_from_guess` checks `pap <= 0.0` (`L1101-1106`) | On even-axis FFT grids, `ker(D)` is spanned by `G = 0` plus up to 7 Nyquist corner parity modes `(-1)^(p · ijk)`. Full-grid SPD holds iff `{i \| k2 i > 0}` hits every parity coset (satisfied by any bulk solvent containing a contiguous `2 × 2 × 2` grid block). |
| `gpOp_spd_const_mode_lpb` | `dv > 0`, `eps ≥ 1`, `k2 ≥ 0`, `∃ i0, 0 < k2 i0`, `c ≠ 0` | `poisson.rs`: `PdeSolver::solve_pde` (`has_kappa = checks[4] > 0.0`, `L843`), `SolvationWorkspace::evaluate_inner` (`has_kappa = global_kappa2_sum[0] > 0.0`, `L2175-2178`) | Guarantees a single bulk electrolyte point lifts the `G = 0` constant-potential null mode and pins `φ_bulk = 0` for charged cells and Phase B2 GC-DFT. |
| `gpOp_spd_uniform_screening` | `dv > 0`, `eps ≥ 1`, `∀ i, 0 < k2 i`, `∃ i, u i ≠ 0` | Uniform-screening MMS fixtures in `crates/dft-physics/tests/solvation_gate2.rs` (`mms_variable_dielectric_and_lpb_cubic_and_skewed_all_solvers`, `L329-470`) | Restricted corollary of `gpOp_spd_lpb_cavity_vanishing`; physical cavities with `k2 = 0` inside the solute rely on `gpOp_spd_lpb_cavity_vanishing` instead. |

---

### 1.2 Dual Variational Stationarity, Envelope Theorem & Exact Vacuum Reduction

| Theorem (`BaandiMath.Solvation`) | Formal Hypotheses | Rust Enforcement Site (File, Function, Lines) | Verification Notes |
| --- | --- | --- | --- |
| `dualEnergy_on_shell`, `dualEnergy_perturb_exact`, `dualEnergy_global_max` | `pi_val > 0`, `A` additive, symmetric, PSD (`gpOp_add`, `gpOp_symmetric`, `gpOp_psd`), `A phi_star = 4 * pi * rho` | `poisson.rs`: `SolvationWorkspace::evaluate_inner` solves `A phi_pol = 4*pi*rho_sol - A phi_vac` (`L2261-2305`), forms `phi_solv = phi_vac + phi_pol` (`L2307-2315`), evaluates `electrostatic_ha = 0.5 * volume * Re⟨rho_sol, phi_pol⟩` (`L2327-2334`), and uses stationary Hellmann–Feynman derivatives `v_es`, `v_diel`, `v_scr` (`L2340-2363`) without implicit `dphi/dn` terms | Stationarity makes the dual electrostatic functional globally maximal at the exact PDE solution and eliminates implicit `δφ/δn` and `∂φ/∂τ_I` response terms from KS potentials and ionic forces. |
| `gpOp_vacuum_eq`, `vacuumOp_symmetric`, `gpOp_apply_zero`, `polarization_rhs_vanishes_in_vacuum` | `eps = 1`, `k2 = 0`, `L0 phi_vac = 4 * pi * rho` using the identical skew-adjoint `D_a` as `gpOp` | `poisson.rs`: `SolvationWorkspace::evaluate_inner` constructs `phi_vac_g` with the identical `self.solver.engine.g2` (`L2208-2220`), executes exact vacuum fast path when `params.is_vacuum() && |q_ion - n_el| < 1e-12` (`L2222-2255`), and forms `b_pol_g` via `eps - 1.0` (`L2261-2277`); `mod.rs`: `SolvationParams::is_vacuum` (`L167-172`) | Using the exact same skew-adjoint `g_cart` / `g2` for both `A` and `L0` ensures `b_pol_g = 0` identically in vacuum (`Δφ_pol = 0`, `ΔE_es = 0`). |

---

### 1.3 Zero Double Counting, Background-Offset Cancellation & PDE Width Invariance

| Theorem (`BaandiMath.Solvation`) | Formal Hypotheses | Rust Enforcement Site (File, Function, Lines) | Verification Notes |
| --- | --- | --- | --- |
| `zero_double_counting_cancellation`, `bg_offset_algebraic_cancellation` | Real scalar decomposition of `E_H + E_loc + E_ewald`, `Delta E_es`, and `E_bg,corr = (Q_ion - N_e) * sum_I alpha_z,I^sm` | `cavity.rs`: `prepare_ions` (`alpha_z_sm = PI * q * r_c^2 / omega`, `L201`), `build_ionic_pseudo_charge_g` (`sum_alpha_z_sm`, `L222-262`); `poisson.rs`: `evaluate_inner` (`background_correction_ha = (q_ion - n_el) * sum_alpha_z_sm`, `L2204-2205`, `total_ha = electrostatic_ha + background_correction_ha + cavitation_ha`, `L2423`) | Proves algebraic cancellation of bare `E_es_vac` and the `(N_e - Q_ion) * alpha_sm` background offset when `state.ionic`, `hartree.rs`, and `ewald.rs` remain untouched. Does **not** alone prove PDE width invariance. |
| `full_electrostatic_width_invariance` | `D_a` additive, subtractive, annihilating constants, skew-adjoint; `L0 phi_vac_k = 4*pi*(rho_k - rho_bar)` with zero mean; Newton-shell non-overlap (`D_a phi_tilde = 0` where `eps ≠ 1` and `phi_tilde = 0` where `k2 ≠ 0`); `A phi_pol1 = 4*pi*rho1 - A phi_vac1` | `cavity.rs`: `build_ionic_pseudo_charge_g` (`L222-262`, Gaussian width `width_bohr`, default `0.60 Bohr` in `crates/dft-engine/src/model/options.rs:28`); `poisson.rs`: `evaluate_inner` (`phi_vac_g` zero-mean at `L2215-2216`, `b_pol_g` at `L2261-2277`, `electrostatic_ha` at `L2327-2334`, `background_correction_ha` at `L2205`, `bg_pot_shift` at `L2340`) | Genuine PDE-level invariance proof: shows `phi_pol2 = phi_pol1 - (alpha_sm2 - alpha_sm1)` solves the shifted polarization PDE and `(1/2)⟨rho, phi_pol⟩ + Q_net * alpha_sm` is strictly width-invariant for neutral and charged cells via Green's reciprocity. |

---

### 1.4 Harris–Foulkes & Grand Canonical Bregman Remainders

| Theorem (`BaandiMath.Solvation`) | Formal Hypotheses | Rust Enforcement Site (File, Function, Lines) | Verification Notes |
| --- | --- | --- | --- |
| `harris_foulkes_bregman_exact` | Hartree homogeneity `⟨n_in, V_H[n_in]⟩ = 2 * E_H[n_in]` and consistent `hxc = V_H + V_xc + V_solv` double counting | `crates/dft-engine/src/refscf/map.rs`: `fields` (`eh = 0.5 * dv * ⟨n, V_H⟩` at `L141`, `hxc += v_total_planes` & `v_solv_integral` at `L194-221`), `audit` (`band_form` `L592-595`, `harris_foulkes` `L596-600`, `functional_difference_bound` `L606-626`); `crates/dft-engine/src/refscf/electronic_map.rs` (`L1026-1040`, `L1122-1141`); `crates/dft-engine/src/energy.rs`: `EnergyAudit` (`L8-44`), `total_hartree` (`L87-94`) | Guarantees `E_direct - E_HF = Q_H + Q_xc + Q_solv + delta_eig`, vanishing to second order in `n_out - n_in` at SCF convergence. |
| `grand_canonical_bregman_exact` | Same as `harris_foulkes_bregman_exact`, with `Omega = E - TS - mu * N_out` | Phase B2 GC-DFT readiness target (`refscf/electronic.rs` & `electronic_map.rs` grand-potential accounting) | Establishes that subtracting `TS + mu * N_out` preserves the exact Bregman remainder decomposition for fixed-`mu` SCF audits. |
| `bregman_bound_nonneg_triangle` | Real remainders `R_H, R_xc, R_solv` | `crates/dft-engine/src/refscf/map.rs:618-621`, `crates/dft-engine/src/refscf/electronic_map.rs:1133-1136` | Justifies the non-negative triangle-inequality upper bound reported in `EnergyAudit::functional_difference_bound`. |

---

### 1.5 Analytic Cavity Switching Derivatives, Monotonicity & Repulsive Potentials

| Theorem (`BaandiMath.Solvation`) | Formal Hypotheses | Rust Enforcement Site (File, Function, Lines) | Verification Notes |
| --- | --- | --- | --- |
| `vaspsolSwitch_hasDerivAt`, `vaspsolSwitch_hasDerivAt_nonpos`, `vaspsol_deriv_nonpos` | `0 < n_c`, `0 < sigma`, `0 < n`, analytic derivative of `erfc` | `crates/dft-physics/src/solvation/mod.rs`: `SolvationParams::validate` (`nc_bohr3 > 0`, `sigma > 0`, `density_floor_bohr3 > 0`, `L201-213`); `crates/dft-physics/src/solvation/cavity.rs`: `evaluate_cavity_fields` (`VaspSol` branch, `L99-118`, guarding `n <= floor` with `(1.0, 0.0)`) | Machine-checked `HasDerivAt.comp` chain rule matching `cavity.rs:106-108` (`x = ln(n / nc) / (sigma * sqrt(2))`, `dsv = -exp(-x^2) / (n * sigma * sqrt(2*pi))`). |
| `sccsSwitch_hasDerivAt`, `sccsSwitch_hasDerivAt_nonpos`, `sccs_deriv_nonpos` | `0 < rho_min < n < rho_max` | `mod.rs`: `SolvationParams::validate` (`0 < rho_min_bohr3 < rho_max_bohr3`, `L214-228`); `cavity.rs`: `evaluate_cavity_fields` (`Sccs` branch, `L119-153`) | Note equivalence of parameterization: `cavity.rs:135-137` uses `t_rust = ln(rho_max / n) / ln(rho_max / rho_min) = 1 - t_lean`; since `sin(2π(1 - t)) = -sin(2π t)` and `cos(2π(1 - t)) = cos(2π t)`, `t_rust - sin(2π t_rust)/(2π) = 1 - (t_lean - sin(2π t_lean)/(2π))` identically. |
| `diel_and_screening_potentials_nonneg` | `1 ≤ eps_b`, `0 ≤ eps`, `0 ≤ kappa`, `ds_dn ≤ 0`, `0 < pi_val`, `0 ≤ grad_phi_sq` | `mod.rs`: `SolvationParams::validate` (`epsilon_bulk >= 1.0`, `kappa_inv_bohr >= 0.0`, `L176-185`); `cavity.rs`: `evaluate_cavity_fields` (`L113-116`, `L140-151`); `poisson.rs`: `reconstruct_potentials_and_grad_sq` (`grad_phi_sq = sum_a (D_a phi)^2 >= 0`, `L476-488`), `evaluate_inner` (`v_diel = -INV_EIGHT_PI * deps * grad2`, `v_scr = -INV_EIGHT_PI * dkap * phi_solv^2`, `L2361-2362`) | Proves repulsive (`≥ 0`) dielectric and ionic screening potentials for both linear-`eps` and SCCS log-`eps` models. |

---

### 1.6 Gate 3–4 SCF Invariants (Sign Safety & Spin-Polarized Chain Rule)

| Theorem (`BaandiMath.Solvation`) | Formal Hypotheses | Rust Enforcement Site (File, Function, Lines) | Verification Notes |
| --- | --- | --- | --- |
| `electrostatic_sign_involution` | Arbitrary `dv`, `rho_e`, `phi` | `crates/dft-units/src/lib.rs`: `ElectrostaticPotentialHaPerE` (`L60-61`), `ElectronPotentialEnergyHa` (`L65-66`), `to_electron_energy` / `to_electrostatic_potential` (`L68-82`); `poisson.rs:2358-2360` (`ElectrostaticPotentialHaPerE(phi_pol + bg_pot_shift).to_electron_energy().0`) | Compile-time type boundary preventing `phi` vs `U_el = -phi` sign inversion (ADR 0017 / ADR 0061). |
| `spin_polarized_solvation_chain_rule` | Total-density functional `E[n_up + n_dn]` Gateaux-differentiable with derivative `V_solv[n_tot]` | `crates/dft-engine/src/refscf/electronic_map.rs`: evaluates scalar `map::fields` on `total = n_up + n_dn` (`L371-385`), adds identical `s.v_solv` to both spin channels before `potential_symmetry` (`L453-459` for PBE, `L527-531` for LDA), and contracts `(n_out_tot - n_in_tot) * v_solv` in `l_solv` (`L1026-1040`) | Proves spin-independent `V_solv[n_tot]` is the exact functional derivative for both `n_up` and `n_dn` and that spin-summed double counting reduces to `⟨n_tot, V_solv⟩`. |

---

## 2. Formal-to-Numerical Verification Boundary

Lean 4 proofs operate over exact mathematical reals (`ℝ`) and exact algebraic identities. The following gaps between the formal model and the Rust numerical implementation are **outside** the scope of Level-1 Lean proofs and are explicitly verified by Level-2/Level-3 Rust test suites:

| Formal Assumption (Lean 4 over `ℝ`) | Numerical Reality (Rust over IEEE-754 `f64`) | Bridging Rust Verification Test |
| --- | --- | --- |
| Exact `ℝ` arithmetic and exact Parseval `⟨u, v⟩_r = ⟨û, v̂⟩_G` | IEEE-754 `f64` rounding, `ndrustfft`/`rustfft` transform rounding, MPI `tree_sum_f64` reduction ordering | `crates/dft-physics/tests/solvation_gate2.rs`: `mpi_rank_independence_and_mirror_vs_greedy_parity` (`< 1e-12` across `P ∈ {1, 2, 4}` and `mirror` vs `greedy`) |
| Exact PDE solution `A phi_star = 4*pi*rho` | Iterative `DirectPcg` / `FisicaroPcg` stopping at finite residual `‖r‖_2 ≤ tol * max(1, ‖b‖_2)` | `solvation_gate2.rs`: `mms_variable_dielectric_and_lpb_cubic_and_skewed_all_solvers` (`< 1e-12`), `variational_potential_matches_central_finite_differences_vaspsol_and_sccs` (`< 1e-9` rel / `< 1e-10` Ha abs) |
| Exact non-overlap (`D_a phi_tilde = 0` on `eps ≠ 1`, `phi_tilde = 0` on `k2 ≠ 0`) | Exponentially small Gaussian tails `exp(-r^2 / r_c^2)` leaking into the dielectric boundary when `r_c` is large | `solvation_gate2.rs`: `ionic_forces_match_central_finite_differences_and_width_invariance` (sweeping `r_c` inside the cavity, verifying sensitivity `< 1e-6 Ha` and force FD `< 1e-9 Ha/Bohr`) |
| Exact vacuum limit `eps = 1, k2 = 0` | Neutral fast path triggers at `|q_ion - n_el| < 1e-12`; otherwise general path solves `b_pol_g = 0` | `solvation_gate2.rs`: `vacuum_limit_is_identically_zero_on_skewed_cell` (`< 1e-14`); `crates/dft-engine/tests/solvation_scf.rs` (vacuum SCF identity vs bare PW) |
| Analytic `erfc` derivative hypothesis `herfc` | `libm::erfc(x)` IEEE-754 approximation and `n <= floor` cutoff clamp | `solvation_gate2.rs`: `variational_potential_matches_central_finite_differences_vaspsol_and_sccs` (VASPsol & SCCS central FD check) |
| Discrete parity-coset coverage for full-grid LPB SPD on even-axis FFT grids | Requires bulk solvent `{i \| k2 i > 0}` to intersect all `2^(# even axes)` parity cosets; `DirectPcg` guards runtime breakdown (`pap <= 0.0`) | `solvation_gate2.rs`: `born_ion_and_debye_huckel_analytic_convergence`, `invalid_solvation_parameters_are_rejected`; `crates/dft-engine/tests/solvation_scf.rs` (LPB SCF convergence & Bregman audit) |
| Exact Gateaux differentiability in SCF (`harris_foulkes_bregman_exact`, `spin_polarized_solvation_chain_rule`) | Finite eigensolver residual `delta_eig`, adaptive inner tolerance ratchet (`next_solvation_tol`), symmetry projection (`ADR 0049`), and `FisicaroPcg` product-rule aliasing (`max_dx <= 0.25 Bohr` guard) | `crates/dft-engine/tests/solvation_scf.rs` (L3 end-to-end SCF energy, Bregman remainder, spin-polarized SCF, finite-difference SCF forces, and `DirectPcg` vs `FisicaroPcg` parity) |

---

## 3. Scope of Optional Rust-to-Lean Extraction (`Aeneas`)

Translating small pure-Rust kernels into Lean via **Aeneas** (for example, `dft-units` scalar conversions or the pointwise `evaluate_cavity_fields` switching math) may be explored in the future as an **optional, self-contained research study** to tighten the gap between Rust source and Lean definitions. Per `AGENTS.md` and `D-014`:
- It is **not** a prerequisite for engine development, Gate 4/5 milestones, or production integration.
- It must never introduce restrictions, allocations, or abstractions that degrade BAANDI's production hot-path performance.
