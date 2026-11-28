# `BaandiMath` — Formal Mathematical Verification for BAANDI Track B

Self-contained Lean 4 (`v4.36.0`) + Mathlib (`v4.36.0`) formalization of the Level-1 physical and algorithmic guarantees for BAANDI's implicit solvation and Grand Canonical DFT formulation ([`docs/research/extensible-dft/implicit-solvation-theory.md`](../../docs/research/extensible-dft/implicit-solvation-theory.md)).

Every theorem is proved with **zero `sorry`** and **zero custom `axiom`** declarations. For the complete hypothesis-by-hypothesis mapping to Rust code and the formal-to-numerical verification boundary, see [`docs/research/extensible-dft/formal-verification-map.md`](../../docs/research/extensible-dft/formal-verification-map.md).

## Reproducible Verification (`dev08` / CI)

Run the automated build, `sorry`/custom-`axiom` scan, and `#print axioms` audit from the repository root:

```bash
bash scripts/check_lean.sh
# Or point at an explicit synced build directory on dev08:
bash scripts/check_lean.sh /scr/mino/lean/BaandiMath
```

On a fresh checkout on `dev08`, fetch Mathlib precompiled caches once before running `scripts/check_lean.sh`:

```bash
(cd formal/BaandiMath && lake exe cache get)
```

## Guarantee Mapping to `implicit-solvation-theory.md`

| Guarantee Group in `Solvation.lean` | Core Theorems | Specification & Rust Mapping |
| --- | --- | --- |
| **1. Nyquist-Zeroed Skew-Adjointness & Realistic PCG-Space SPD** | `nyquist_zeroed_multiplier_skew_adjoint`, `nyquist_nonzero_breaks_reality`, `gpOp_quad_form`, `gpOp_symmetric`, `gpOp_psd`, `gpOp_quad_eq_zero_iff`, `gpOp_spd_on_grad_complement`, `gpOp_spd_lpb_cavity_vanishing`, `gpOp_spd_const_mode_lpb`, `gpOp_spd_uniform_screening` | §4.1 (`D_a^T = -D_a`, real-field closure via Nyquist zeroing, pure-GP `zero_g2_indices` subspace SPD, cavity-vanishing LPB SPD, single-point bulk `G = 0` gauge pinning) |
| **2. Dual Variational Stationarity & Exact Vacuum Limit** | `dualEnergy_on_shell`, `dualEnergy_perturb_exact`, `dualEnergy_global_max`, `gpOp_vacuum_eq`, `vacuumOp_symmetric`, `gpOp_apply_zero`, `polarization_rhs_vanishes_in_vacuum` | §4.2 (dual functional $E_{\mathrm{es}}$, envelope/stationarity, $\Delta\phi_{\mathrm{pol}} = 0$ vacuum reduction) |
| **3. Zero Double Counting & PDE Width Invariance** | `zero_double_counting_cancellation`, `bg_offset_algebraic_cancellation`, `full_electrostatic_width_invariance` | §4.2 (exact reconciliation with `state.ionic`/`ewald`, algebraic `G = 0` cancellation, and PDE-level Green's reciprocity proof of $\frac{1}{2}\langle\rho, \Delta\phi_{\mathrm{pol}}\rangle + Q_{\mathrm{net}}\alpha_{\mathrm{sm}}$ width invariance for neutral & charged cells) |
| **4. Harris–Foulkes & GC-DFT Bregman Identities** | `harris_foulkes_bregman_exact`, `grand_canonical_bregman_exact` | §4.5 ($E_{\mathrm{direct}} - E_{\mathrm{HF}}$ and $\Omega_{\mathrm{direct}} - \Omega_{\mathrm{HF}}$ with $Q_{\mathrm{solv}}$) |
| **5. Analytic Cavity Switching & Repulsive Potentials** | `vaspsol_deriv_nonpos`, `sccs_deriv_nonpos`, `sccsSwitch_hasDerivAt`, `sccsSwitch_hasDerivAt_nonpos`, `vaspsolSwitch_hasDerivAt`, `vaspsolSwitch_hasDerivAt_nonpos`, `diel_and_screening_potentials_nonneg` | §§3, 4.3 (machine-checked `HasDerivAt` for SCCS log-sine and VASPsol `erfc` chain rule, $S' \le 0$, $V_{\mathrm{diel}} \ge 0$ for linear-$\epsilon$ and log-$\epsilon$, $V_{\mathrm{scr}} \ge 0$) |
| **6. Gate 3–4 SCF Invariants** | `electrostatic_sign_involution`, `spin_polarized_solvation_chain_rule`, `bregman_bound_nonneg_triangle` | ADR 0017 / ADR 0061 (`U_el = -phi`, spin-polarized Gateaux chain rule & double counting, `EnergyAudit::functional_difference_bound` triangle bound) |
