# `BaandiMath` — Formal Mathematical Verification for BAANDI Track B

Self-contained Lean 4 (`v4.36.0`) + Mathlib (`v4.36.0`) formalization of the Level-1 theoretical invariants for BAANDI's implicit solvation and Grand Canonical DFT formulation ([`docs/research/extensible-dft/implicit-solvation-theory.md`](../../docs/research/extensible-dft/implicit-solvation-theory.md)).

Every theorem is proved with **zero `sorry`** and **zero custom `axiom`** declarations.

## Building on `dev08`

Run from `formal/BaandiMath/`:

```bash
lake exe cache get && lake build
```

`BaandiMath/Solvation.lean` concludes with `#print axioms` checks for all main theorems so `lake build` prints an axiom audit confirming only standard Lean kernel axioms (`propext`, `Classical.choice`, `Quot.sound`).

## Theorem Mapping to `implicit-solvation-theory.md`

| Section in `Solvation.lean` | Theorems | Specification Reference |
| --- | --- | --- |
| **1. Weighted Inner Product & Operator SPD** | `gpOp_quad_form`, `gpOp_symmetric`, `gpOp_psd`, `gpOp_spd` | Section 4.1 (skew-adjoint $D_a^\top = -D_a$, SPD GP/LPB operator $A$) |
| **2. Dual Variational Stationarity & Vacuum Limit** | `dualEnergy_on_shell`, `dualEnergy_perturb_exact`, `dualEnergy_global_max`, `gpOp_vacuum_eq`, `gpOp_apply_zero`, `polarization_rhs_vanishes_in_vacuum` | Section 4.2 (dual functional $E_{\mathrm{es}}$, envelope/stationarity, $\Delta\phi_{\mathrm{pol}} = 0$ vacuum limit) |
| **3. Zero Double Counting & Width Invariance** | `zero_double_counting_cancellation`, `width_invariance_arbitrary_charge` | Section 4.2 ($E_{\mathrm{bg,corr}} = (Q_{\mathrm{ion}} - N_e)\sum_I \alpha_{Z,I}^{\mathrm{sm}}$, $r_{c,I}$-invariance for neutral & charged cells) |
| **4. Harris–Foulkes & GC-DFT Bregman Identities** | `harris_foulkes_bregman_exact`, `grand_canonical_bregman_exact` | Section 4.5 ($E_{\mathrm{direct}} - E_{\mathrm{HF}}$ and $\Omega_{\mathrm{direct}} - \Omega_{\mathrm{HF}}$ with $Q_{\mathrm{solv}}$) |
| **5. Cavity Monotonicity & Repulsive Potentials** | `vaspsol_deriv_nonpos`, `sccs_deriv_nonpos`, `diel_and_screening_potentials_nonneg` | Sections 3 & 4.3 ($S'_{\mathrm{VASPsol}} \le 0$, $S'_{\mathrm{SCCS}} \le 0$, $V_{\mathrm{diel}} \ge 0$ for linear-$\epsilon$ and log-$\epsilon$, $V_{\mathrm{scr}} \ge 0$) |
