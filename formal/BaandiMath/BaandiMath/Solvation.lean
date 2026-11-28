import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.Tactic

/-!
# BAANDI Track B (Implicit Solvation & GC-DFT) — Formal Mathematical Verification

Machine-checked Lean 4 + Mathlib proofs of the core Level-1 theoretical theorems in
`docs/research/extensible-dft/implicit-solvation-theory.md` (Section 4) and Gate 3–4
SCF invariants:

1. Weighted grid inner product & strict positive definiteness (SPD) of the discrete
   generalized Poisson / linearized Poisson–Boltzmann (GP/LPB) operator built from
   skew-adjoint FFT gradients `D_a^T = -D_a`.
2. Dual variational electrostatic energy functional: on-shell identity, exact quadratic
   perturbation identity, global maximality (envelope theorem / stationarity), and exact
   vacuum reduction (`phi_pol = 0`, `Delta E_es = 0`).
3. Exact zero double counting & Gaussian smearing-width invariance (`E_bg,corr`) for
   both neutral and charged cells.
4. Harris–Foulkes Bregman remainder & Grand Canonical (`Omega = E - TS - mu * N_out`)
   Bregman identity.
5. Analytic `HasDerivAt` cavity switching derivatives & monotonicity (`VASPsol` `erfc`
   via `HasDerivAt.comp` and `SCCS` log-sine directly from `Real.log` and `Real.sin`)
   and non-negativity (repulsive sign) of the dielectric (`V_diel`, both linear-`eps`
   and log-`eps`) and electrolyte screening (`V_scr`) Kohn–Sham potentials.
6. Gate 3–4 SCF invariants: electrostatic sign involution (`U_el = -phi`),
   spin-polarized solvation Gateaux chain rule and double-counting identity, and
   non-negative triangle bound for `EnergyAudit::functional_difference_bound`.

All theorems are proved with zero `sorry` and zero custom `axiom` declarations.
-/

namespace BaandiMath.Solvation

/-!
## Section 1: Weighted Grid Inner Product & Skew-Adjoint Operator SPD
-/

/-- Weighted discrete real-space grid inner product `⟨u, v⟩ = dv * ∑ i, u i * v i`. -/
def wdot {n : ℕ} (dv : ℝ) (u v : Fin n → ℝ) : ℝ :=
  dv * ∑ i : Fin n, u i * v i

theorem wdot_comm {n : ℕ} (dv : ℝ) (u v : Fin n → ℝ) :
    wdot dv u v = wdot dv v u := by
  unfold wdot
  congr 1
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem wdot_add_right {n : ℕ} (dv : ℝ) (u v w : Fin n → ℝ) :
    wdot dv u (fun i => v i + w i) = wdot dv u v + wdot dv u w := by
  unfold wdot
  have hsum : ∑ i : Fin n, u i * (v i + w i) =
      (∑ i : Fin n, u i * v i) + ∑ i : Fin n, u i * w i := by
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro i _
    ring
  rw [hsum]
  ring

theorem wdot_add_left {n : ℕ} (dv : ℝ) (u v w : Fin n → ℝ) :
    wdot dv (fun i => u i + v i) w = wdot dv u w + wdot dv v w := by
  rw [wdot_comm, wdot_add_right, wdot_comm dv w u, wdot_comm dv w v]

theorem wdot_neg_right {n : ℕ} (dv : ℝ) (u v : Fin n → ℝ) :
    wdot dv u (fun i => -v i) = -wdot dv u v := by
  unfold wdot
  have hsum : ∑ i : Fin n, u i * (-v i) = -(∑ i : Fin n, u i * v i) := by
    rw [← Finset.sum_neg_distrib]
    apply Finset.sum_congr rfl
    intro i _
    ring
  rw [hsum]
  ring

theorem wdot_smul_right {n : ℕ} (dv c : ℝ) (u v : Fin n → ℝ) :
    wdot dv u (fun i => c * v i) = c * wdot dv u v := by
  unfold wdot
  have hsum : ∑ i : Fin n, u i * (c * v i) = c * ∑ i : Fin n, u i * v i := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    ring
  rw [hsum]
  ring

theorem wdot_mul_diag {n : ℕ} (dv : ℝ) (w u v : Fin n → ℝ) :
    wdot dv u (fun i => w i * v i) = wdot dv (fun i => w i * u i) v := by
  unfold wdot
  congr 1
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem wdot_sum_right {n m : ℕ} (dv : ℝ) (u : Fin n → ℝ)
    (f : Fin m → Fin n → ℝ) :
    wdot dv u (fun i => ∑ a : Fin m, f a i) = ∑ a : Fin m, wdot dv u (f a) := by
  unfold wdot
  rw [← Finset.mul_sum]
  congr 1
  calc
    ∑ i : Fin n, u i * ∑ a : Fin m, f a i
      = ∑ i : Fin n, ∑ a : Fin m, u i * f a i := by
          apply Finset.sum_congr rfl
          intro i _
          rw [Finset.mul_sum]
    _ = ∑ a : Fin m, ∑ i : Fin n, u i * f a i := Finset.sum_comm

theorem wdot_diag_nonneg {n : ℕ} {dv : ℝ} (hdv : 0 ≤ dv) {w u : Fin n → ℝ}
    (hw : ∀ i, 0 ≤ w i) :
    0 ≤ wdot dv (fun i => w i * u i) u := by
  unfold wdot
  apply mul_nonneg hdv
  apply Finset.sum_nonneg
  intro i _
  have hterm : (w i * u i) * u i = w i * (u i * u i) := by ring
  rw [hterm]
  exact mul_nonneg (hw i) (mul_self_nonneg (u i))

theorem wdot_diag_pos {n : ℕ} {dv : ℝ} (hdv : 0 < dv) {w u : Fin n → ℝ}
    (hw : ∀ i, 0 < w i) (hu : ∃ i, u i ≠ 0) :
    0 < wdot dv (fun i => w i * u i) u := by
  rcases hu with ⟨i0, hi0⟩
  unfold wdot
  apply mul_pos hdv
  have hnonneg : ∀ i ∈ (Finset.univ : Finset (Fin n)), 0 ≤ (w i * u i) * u i := by
    intro i _
    have hterm : (w i * u i) * u i = w i * (u i * u i) := by ring
    rw [hterm]
    exact mul_nonneg (le_of_lt (hw i)) (mul_self_nonneg (u i))
  have hsq_pos : 0 < u i0 * u i0 := by
    rcases lt_trichotomy (u i0) 0 with hneg | hzero | hpos
    · nlinarith
    · exact absurd hzero hi0
    · nlinarith
  have hpos_i0 : 0 < (w i0 * u i0) * u i0 := by
    have hterm : (w i0 * u i0) * u i0 = w i0 * (u i0 * u i0) := by ring
    rw [hterm]
    exact mul_pos (hw i0) hsq_pos
  have hle : (w i0 * u i0) * u i0 ≤ ∑ i : Fin n, (w i * u i) * u i :=
    Finset.single_le_sum hnonneg (Finset.mem_univ i0)
  exact lt_of_lt_of_le hpos_i0 hle

/-- Discrete Generalized Poisson / Linearized Poisson–Boltzmann operator
`A u = ∑_{a=0}^2 D_a^T (eps * D_a u) + k2 * u` with skew-adjoint `D_a^T = -D_a`. -/
def gpOp {n : ℕ} (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (eps k2 u : Fin n → ℝ) (i : Fin n) : ℝ :=
  (∑ a : Fin 3, -(D a (fun j => eps j * D a u j) i)) + k2 i * u i

/-- Bilinear / quadratic form identity for `gpOp` under skew-adjoint discrete gradients. -/
theorem gpOp_quad_form {n : ℕ} (dv : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 u v : Fin n → ℝ) :
    wdot dv u (gpOp D eps k2 v) =
      (∑ a : Fin 3, wdot dv (fun i => eps i * D a u i) (D a v)) +
      wdot dv (fun i => k2 i * u i) v := by
  change wdot dv u
    (fun i => (∑ a : Fin 3, -(D a (fun j => eps j * D a v j) i)) + k2 i * v i) = _
  rw [wdot_add_right, wdot_sum_right, wdot_mul_diag]
  congr 1
  apply Finset.sum_congr rfl
  intro a _
  rw [wdot_neg_right, hskew, neg_neg, wdot_mul_diag]

/-- Exact symmetry of `gpOp` with respect to the weighted grid inner product `wdot dv`. -/
theorem gpOp_symmetric {n : ℕ} (dv : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 u v : Fin n → ℝ) :
    wdot dv u (gpOp D eps k2 v) = wdot dv (gpOp D eps k2 u) v := by
  rw [gpOp_quad_form dv D hskew eps k2 u v]
  rw [wdot_comm dv (gpOp D eps k2 u) v]
  rw [gpOp_quad_form dv D hskew eps k2 v u]
  congr 1
  · apply Finset.sum_congr rfl
    intro a _
    unfold wdot
    congr 1
    apply Finset.sum_congr rfl
    intro i _
    ring
  · unfold wdot
    congr 1
    apply Finset.sum_congr rfl
    intro i _
    ring

/-- Positive semidefiniteness of `gpOp` whenever `dv > 0`, `eps ≥ 1`, and `k2 ≥ 0`. -/
theorem gpOp_psd {n : ℕ} (dv : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 u : Fin n → ℝ)
    (hdv : 0 < dv)
    (heps : ∀ i, 1 ≤ eps i)
    (hk2 : ∀ i, 0 ≤ k2 i) :
    0 ≤ wdot dv u (gpOp D eps k2 u) := by
  rw [gpOp_quad_form dv D hskew eps k2 u u]
  have hgrad : 0 ≤ ∑ a : Fin 3, wdot dv (fun i => eps i * D a u i) (D a u) := by
    apply Finset.sum_nonneg
    intro a _
    exact wdot_diag_nonneg (le_of_lt hdv) (fun i => by linarith [heps i])
  have hscr : 0 ≤ wdot dv (fun i => k2 i * u i) u :=
    wdot_diag_nonneg (le_of_lt hdv) hk2
  linarith

/-- Strict positive definiteness (SPD) of the LPB operator `gpOp` when `k2 > 0`
and `u` is nonzero. -/
theorem gpOp_spd {n : ℕ} (dv : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 u : Fin n → ℝ)
    (hdv : 0 < dv)
    (heps : ∀ i, 1 ≤ eps i)
    (hk2 : ∀ i, 0 < k2 i)
    (hu : ∃ i, u i ≠ 0) :
    0 < wdot dv u (gpOp D eps k2 u) := by
  rw [gpOp_quad_form dv D hskew eps k2 u u]
  have hgrad : 0 ≤ ∑ a : Fin 3, wdot dv (fun i => eps i * D a u i) (D a u) := by
    apply Finset.sum_nonneg
    intro a _
    exact wdot_diag_nonneg (le_of_lt hdv) (fun i => by linarith [heps i])
  have hscr : 0 < wdot dv (fun i => k2 i * u i) u :=
    wdot_diag_pos hdv hk2 hu
  linarith

/-!
## Section 2: Dual Variational Stationarity, Global Maximality & Exact Vacuum Reduction
-/

/-- Dual variational electrostatic functional
`E_dual[rho, phi] = ⟨rho, phi⟩ - (1 / (8 * pi)) * ⟨phi, A phi⟩`. -/
noncomputable def dualEnergy {n : ℕ} (dv pi_val : ℝ)
    (A : (Fin n → ℝ) → (Fin n → ℝ)) (rho phi : Fin n → ℝ) : ℝ :=
  wdot dv rho phi - (1 / (8 * pi_val)) * wdot dv phi (A phi)

/-- On-shell evaluation: when `A phi_star = 4 * pi * rho`,
`dualEnergy = (1 / 2) * ⟨rho, phi_star⟩`. -/
theorem dualEnergy_on_shell {n : ℕ} (dv pi_val : ℝ)
    (A : (Fin n → ℝ) → (Fin n → ℝ)) (rho phi_star : Fin n → ℝ)
    (hpi : pi_val ≠ 0)
    (hsol : ∀ i, A phi_star i = 4 * pi_val * rho i) :
    dualEnergy dv pi_val A rho phi_star = (1 / 2) * wdot dv rho phi_star := by
  unfold dualEnergy
  have hA : A phi_star = fun i => (4 * pi_val) * rho i := funext hsol
  have h8 : 8 * pi_val ≠ 0 := mul_ne_zero (by norm_num) hpi
  have hcoeff : (1 / (8 * pi_val)) * (4 * pi_val) = 1 / 2 := by
    rw [div_mul_eq_mul_div, div_eq_iff h8]
    ring
  rw [hA, wdot_smul_right, wdot_comm dv phi_star rho, ← mul_assoc, hcoeff]
  ring

/-- Exact second-order perturbation identity around the stationary solution
`A phi_star = 4 * pi * rho`: the first-order variation vanishes identically. -/
theorem dualEnergy_perturb_exact {n : ℕ} (dv pi_val : ℝ)
    (A : (Fin n → ℝ) → (Fin n → ℝ))
    (hadd : ∀ u v i, A (fun j => u j + v j) i = A u i + A v i)
    (hsym : ∀ u v, wdot dv u (A v) = wdot dv (A u) v)
    (rho phi_star eta : Fin n → ℝ)
    (hpi : 0 < pi_val)
    (hsol : ∀ i, A phi_star i = 4 * pi_val * rho i) :
    dualEnergy dv pi_val A rho (fun i => phi_star i + eta i) =
      dualEnergy dv pi_val A rho phi_star -
        (1 / (8 * pi_val)) * wdot dv eta (A eta) := by
  unfold dualEnergy
  have hA_add : A (fun j => phi_star j + eta j) = fun i => A phi_star i + A eta i :=
    funext (hadd phi_star eta)
  have hA_star : A phi_star = fun i => (4 * pi_val) * rho i := funext hsol
  have hquad :
      wdot dv (fun i => phi_star i + eta i) (A (fun i => phi_star i + eta i)) =
        wdot dv phi_star (A phi_star) +
          2 * wdot dv eta (A phi_star) + wdot dv eta (A eta) := by
    rw [hA_add, wdot_add_left, wdot_add_right, wdot_add_right]
    have hcross : wdot dv phi_star (A eta) = wdot dv eta (A phi_star) := by
      rw [hsym phi_star eta, wdot_comm]
    rw [hcross]
    ring
  have hlin_cancel :
      (1 / (8 * pi_val)) * (2 * wdot dv eta (A phi_star)) = wdot dv rho eta := by
    rw [hA_star, wdot_smul_right, wdot_comm dv eta rho]
    have h8 : 8 * pi_val ≠ 0 := mul_ne_zero (by norm_num) (ne_of_gt hpi)
    have hcoeff : (1 / (8 * pi_val)) * (2 * (4 * pi_val)) = 1 := by
      rw [div_mul_eq_mul_div, div_eq_iff h8]
      ring
    calc
      (1 / (8 * pi_val)) * (2 * (4 * pi_val * wdot dv rho eta))
        = ((1 / (8 * pi_val)) * (2 * (4 * pi_val))) * wdot dv rho eta := by ring
      _ = 1 * wdot dv rho eta := by rw [hcoeff]
      _ = wdot dv rho eta := one_mul _
  rw [wdot_add_right, hquad]
  calc
    wdot dv rho phi_star + wdot dv rho eta -
        (1 / (8 * pi_val)) *
          (wdot dv phi_star (A phi_star) +
            2 * wdot dv eta (A phi_star) + wdot dv eta (A eta))
      = (wdot dv rho phi_star - (1 / (8 * pi_val)) * wdot dv phi_star (A phi_star)) -
          (1 / (8 * pi_val)) * wdot dv eta (A eta) +
          (wdot dv rho eta - (1 / (8 * pi_val)) * (2 * wdot dv eta (A phi_star))) := by
        ring
    _ = (wdot dv rho phi_star - (1 / (8 * pi_val)) * wdot dv phi_star (A phi_star)) -
          (1 / (8 * pi_val)) * wdot dv eta (A eta) := by
        rw [hlin_cancel]
        ring

/-- Global maximality of `dualEnergy` at the exact PDE solution `phi_star` when `A`
is positive semidefinite (implies stationarity and the envelope theorem). -/
theorem dualEnergy_global_max {n : ℕ} (dv pi_val : ℝ)
    (A : (Fin n → ℝ) → (Fin n → ℝ))
    (hadd : ∀ u v i, A (fun j => u j + v j) i = A u i + A v i)
    (hsym : ∀ u v, wdot dv u (A v) = wdot dv (A u) v)
    (hpsd : ∀ v, 0 ≤ wdot dv v (A v))
    (rho phi_star eta : Fin n → ℝ)
    (hpi : 0 < pi_val)
    (hsol : ∀ i, A phi_star i = 4 * pi_val * rho i) :
    dualEnergy dv pi_val A rho (fun i => phi_star i + eta i) ≤
      dualEnergy dv pi_val A rho phi_star := by
  rw [dualEnergy_perturb_exact dv pi_val A hadd hsym rho phi_star eta hpi hsol]
  have hcoeff : 0 ≤ 1 / (8 * pi_val) :=
    div_nonneg (by norm_num) (by linarith)
  have hterm : 0 ≤ (1 / (8 * pi_val)) * wdot dv eta (A eta) :=
    mul_nonneg hcoeff (hpsd eta)
  linarith

/-- Vacuum Laplacian operator `L0 = ∑_{a=0}^2 D_a^T D_a` using the identical
skew-adjoint discrete gradients `D_a` as `gpOp`. -/
def vacuumOp {n : ℕ} (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (u : Fin n → ℝ) (i : Fin n) : ℝ :=
  ∑ a : Fin 3, -(D a (D a u) i)

/-- In vacuum (`eps = 1`, `k2 = 0`), `gpOp` coincides pointwise with `vacuumOp`. -/
theorem gpOp_vacuum_eq {n : ℕ} (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (u : Fin n → ℝ) (i : Fin n) :
    gpOp D (fun _ => 1) (fun _ => 0) u i = vacuumOp D u i := by
  unfold gpOp vacuumOp
  rw [zero_mul, add_zero]
  apply Finset.sum_congr rfl
  intro a _
  have h1 : (fun j => (1 : ℝ) * D a u j) = D a u := by
    funext j
    exact one_mul _
  rw [h1]

/-- Any linear gradient (`D_a 0 = 0`) maps the zero potential to zero under `gpOp`. -/
theorem gpOp_apply_zero {n : ℕ}
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hD0 : ∀ a i, D a (fun _ => 0) i = 0)
    (eps k2 : Fin n → ℝ) (i : Fin n) :
    gpOp D eps k2 (fun _ => 0) i = 0 := by
  unfold gpOp
  have hinner : ∀ a : Fin 3, (fun j => eps j * D a (fun _ => 0) j) = fun _ => 0 := by
    intro a
    funext j
    rw [hD0 a j, mul_zero]
  have hsum : ∑ a : Fin 3, -(D a (fun j => eps j * D a (fun _ => 0) j) i) = 0 := by
    have hterm : ∀ a : Fin 3, -(D a (fun j => eps j * D a (fun _ => 0) j) i) = 0 := by
      intro a
      rw [hinner a, hD0 a i, neg_zero]
    calc
      ∑ a : Fin 3, -(D a (fun j => eps j * D a (fun _ => 0) j) i)
        = ∑ _a : Fin 3, (0 : ℝ) := by
            apply Finset.sum_congr rfl
            intro a _
            exact hterm a
      _ = 0 := Finset.sum_const_zero
  rw [hsum, mul_zero, add_zero]

/-- Exact vacuum reduction: when `eps = 1`, `k2 = 0`, and `L0 phi_vac = 4 * pi * rho`,
the polarization PDE RHS vanishes identically, `phi_pol = 0` yields `Delta E_es = 0`,
and the dual electrostatic solvation energy difference vanishes identically. -/
theorem polarization_rhs_vanishes_in_vacuum {n : ℕ} (dv pi_val : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (rho phi_vac : Fin n → ℝ)
    (hvac : ∀ i, vacuumOp D phi_vac i = 4 * pi_val * rho i) :
    (∀ i, 4 * pi_val * rho i - gpOp D (fun _ => 1) (fun _ => 0) phi_vac i = 0) ∧
    ((1 / 2 : ℝ) * wdot dv rho (fun _ => 0) = 0) ∧
    (dualEnergy dv pi_val (gpOp D (fun _ => 1) (fun _ => 0)) rho phi_vac -
      dualEnergy dv pi_val (vacuumOp D) rho phi_vac = 0) := by
  refine ⟨?_, ?_, ?_⟩
  · intro i
    rw [gpOp_vacuum_eq, hvac i, sub_self]
  · unfold wdot
    have hsum : ∑ i : Fin n, rho i * (0 : ℝ) = 0 := by
      simp only [mul_zero, Finset.sum_const_zero]
    rw [hsum]
    ring
  · unfold dualEnergy
    have hop : gpOp D (fun _ => 1) (fun _ => 0) phi_vac = vacuumOp D phi_vac :=
      funext (gpOp_vacuum_eq D phi_vac)
    rw [hop, sub_self]

/-!
## Section 3: Exact Zero Double Counting & Gaussian Width Invariance (`E_bg,corr`)
-/

/-- Exact cancellation of bare vacuum electrostatics `E_es_vac` and the `G = 0`
Gaussian background offset `(N_e - Q_ion) * alpha_sm` when adding
`Delta E_es + E_bg,corr` to BAANDI's unmodified `E_H + E_loc + E_ewald`,
holding identically for both neutral (`N_e = Q_ion`) and charged (`N_e ≠ Q_ion`) cells. -/
theorem zero_double_counting_cancellation
    (E_es_vac S P E_loc_short E_es_solv N_e Q_ion alpha_sm : ℝ) :
    (E_es_vac - (S - P) + E_loc_short + (N_e - Q_ion) * alpha_sm) +
      (E_es_solv - E_es_vac) + (Q_ion - N_e) * alpha_sm =
    E_es_solv - (S - P) + E_loc_short := by
  ring

/-- Strict invariance of the reconciled total electrostatic energy with respect to the
Gaussian smearing-width parameter `alpha_sm` for arbitrary cell charge `Q_ion - N_e`. -/
theorem width_invariance_arbitrary_charge
    (E_es_vac S P E_loc_short E_es_solv N_e Q_ion alpha_sm1 alpha_sm2 : ℝ) :
    ((E_es_vac - (S - P) + E_loc_short + (N_e - Q_ion) * alpha_sm1) +
      (E_es_solv - E_es_vac) + (Q_ion - N_e) * alpha_sm1) =
    ((E_es_vac - (S - P) + E_loc_short + (N_e - Q_ion) * alpha_sm2) +
      (E_es_solv - E_es_vac) + (Q_ion - N_e) * alpha_sm2) := by
  ring

/-!
## Section 4: Harris–Foulkes Bregman Remainder & Grand Canonical Identity
-/

/-- Exact Harris–Foulkes Bregman remainder identity with implicit solvation:
`E_direct - E_HF = Q_H + Q_xc + Q_solv + delta_eig` under the Hartree degree-2
homogeneity identity `⟨n_in, V_H[n_in]⟩ = 2 * E_H[n_in]`. -/
theorem harris_foulkes_bregman_exact
    (B E_ewald delta_eig : ℝ)
    (E_H_in E_H_out E_xc_in E_xc_out E_solv_in E_solv_out : ℝ)
    (ip_in_VH ip_out_VH ip_in_Vxc ip_out_Vxc ip_in_Vsolv ip_out_Vsolv : ℝ)
    (h_hartree_homog : ip_in_VH = 2 * E_H_in) :
    let E_direct :=
      B - (ip_out_VH + ip_out_Vxc + ip_out_Vsolv) +
        E_H_out + E_xc_out + E_solv_out + E_ewald + delta_eig
    let E_HF :=
      B - E_H_in + (E_xc_in - ip_in_Vxc) + (E_solv_in - ip_in_Vsolv) + E_ewald
    let Q_H := E_H_out - E_H_in - (ip_out_VH - ip_in_VH)
    let Q_xc := E_xc_out - E_xc_in - (ip_out_Vxc - ip_in_Vxc)
    let Q_solv := E_solv_out - E_solv_in - (ip_out_Vsolv - ip_in_Vsolv)
    E_direct - E_HF = Q_H + Q_xc + Q_solv + delta_eig := by
  dsimp only
  rw [h_hartree_homog]
  ring

/-- Exact Grand Canonical Bregman remainder identity (`Omega = E - TS - mu * N_out`):
`Omega_direct - Omega_HF = Q_H + Q_xc + Q_solv + delta_eig`. -/
theorem grand_canonical_bregman_exact
    (B E_ewald delta_eig TS mu N_out : ℝ)
    (E_H_in E_H_out E_xc_in E_xc_out E_solv_in E_solv_out : ℝ)
    (ip_in_VH ip_out_VH ip_in_Vxc ip_out_Vxc ip_in_Vsolv ip_out_Vsolv : ℝ)
    (h_hartree_homog : ip_in_VH = 2 * E_H_in) :
    let E_direct :=
      B - (ip_out_VH + ip_out_Vxc + ip_out_Vsolv) +
        E_H_out + E_xc_out + E_solv_out + E_ewald + delta_eig
    let E_HF :=
      B - E_H_in + (E_xc_in - ip_in_Vxc) + (E_solv_in - ip_in_Vsolv) + E_ewald
    let Omega_direct := E_direct - TS - mu * N_out
    let Omega_HF := E_HF - TS - mu * N_out
    let Q_H := E_H_out - E_H_in - (ip_out_VH - ip_in_VH)
    let Q_xc := E_xc_out - E_xc_in - (ip_out_Vxc - ip_in_Vxc)
    let Q_solv := E_solv_out - E_solv_in - (ip_out_Vsolv - ip_in_Vsolv)
    Omega_direct - Omega_HF = Q_H + Q_xc + Q_solv + delta_eig := by
  dsimp only
  rw [h_hartree_homog]
  ring

/-!
## Section 5: Analytic Cavity Switching Derivatives, Monotonicity & Repulsive Potentials
-/

/-- Monotonicity of VASPsol `erfc` cavity switching: `S'_VASPsol(n) ≤ 0`. -/
theorem vaspsol_deriv_nonpos (n sigma sqrt_two_pi x : ℝ)
    (hn : 0 < n) (hsigma : 0 < sigma) (hsqrt : 0 < sqrt_two_pi) :
    -(Real.exp (-(x ^ 2))) / (n * sigma * sqrt_two_pi) ≤ 0 := by
  have hnum : -(Real.exp (-(x ^ 2))) ≤ 0 := by
    have hexp : 0 < Real.exp (-(x ^ 2)) := Real.exp_pos (-(x ^ 2))
    linarith
  have hden : 0 ≤ n * sigma * sqrt_two_pi :=
    le_of_lt (mul_pos (mul_pos hn hsigma) hsqrt)
  exact div_nonpos_of_nonpos_of_nonneg hnum hden

/-- Monotonicity of SCCS logarithmic-sine cavity switching: `S'_SCCS(n) ≤ 0`. -/
theorem sccs_deriv_nonpos (n rho_min rho_max t : ℝ)
    (hn : 0 < n) (hmin : 0 < rho_min) (hlt : rho_min < rho_max) :
    -(1 - Real.cos (2 * Real.pi * t)) / (n * Real.log (rho_max / rho_min)) ≤ 0 := by
  have hcos : Real.cos (2 * Real.pi * t) ≤ 1 := Real.cos_le_one (2 * Real.pi * t)
  have hnum : -(1 - Real.cos (2 * Real.pi * t)) ≤ 0 := by linarith
  have hratio : 1 < rho_max / rho_min := (one_lt_div hmin).mpr hlt
  have hlog : 0 < Real.log (rho_max / rho_min) := Real.log_pos hratio
  have hden : 0 ≤ n * Real.log (rho_max / rho_min) := le_of_lt (mul_pos hn hlog)
  exact div_nonpos_of_nonpos_of_nonneg hnum hden

/-- Explicit Andreussi et al. (2012) SCCS logarithmic-sine cavity switching function
on `(rhoMin, rhoMax)`. -/
noncomputable def sccsSwitch (rhoMin rhoMax n : ℝ) : ℝ :=
  let t := (Real.log n - Real.log rhoMin) / (Real.log rhoMax - Real.log rhoMin)
  1 - (t - Real.sin (2 * Real.pi * t) / (2 * Real.pi))

/-- Analytic derivative of `sccsSwitch` on `(rhoMin, rhoMax)` machine-checked via
Lean's `HasDerivAt` calculus library. -/
theorem sccsSwitch_hasDerivAt (rhoMin rhoMax n : ℝ)
    (hmin : 0 < rhoMin) (hn_lo : rhoMin < n) (hn_hi : n < rhoMax) :
    HasDerivAt (sccsSwitch rhoMin rhoMax)
      (-(1 - Real.cos (2 * Real.pi *
          ((Real.log n - Real.log rhoMin) / (Real.log rhoMax - Real.log rhoMin)))) /
        (n * Real.log (rhoMax / rhoMin))) n := by
  have hn : 0 < n := lt_trans hmin hn_lo
  have hlt : rhoMin < rhoMax := lt_trans hn_lo hn_hi
  have hmax : 0 < rhoMax := lt_trans hn hn_hi
  have hlog_eq : Real.log (rhoMax / rhoMin) = Real.log rhoMax - Real.log rhoMin :=
    Real.log_div (ne_of_gt hmax) (ne_of_gt hmin)
  have hlog_pos : 0 < Real.log (rhoMax / rhoMin) :=
    Real.log_pos ((one_lt_div hmin).mpr hlt)
  have hL_ne : Real.log rhoMax - Real.log rhoMin ≠ 0 := by
    rw [← hlog_eq]
    exact ne_of_gt hlog_pos
  have h2pi_ne : 2 * Real.pi ≠ 0 := mul_ne_zero (by norm_num) (ne_of_gt Real.pi_pos)
  have ht : HasDerivAt
      (fun x => (Real.log x - Real.log rhoMin) / (Real.log rhoMax - Real.log rhoMin))
      (n⁻¹ / (Real.log rhoMax - Real.log rhoMin)) n :=
    ((Real.hasDerivAt_log (ne_of_gt hn)).sub_const (Real.log rhoMin)).div_const
      (Real.log rhoMax - Real.log rhoMin)
  have h2pit : HasDerivAt
      (fun x => 2 * Real.pi *
        ((Real.log x - Real.log rhoMin) / (Real.log rhoMax - Real.log rhoMin)))
      (2 * Real.pi * (n⁻¹ / (Real.log rhoMax - Real.log rhoMin))) n :=
    ht.const_mul (2 * Real.pi)
  have hsin : HasDerivAt
      (fun x => Real.sin (2 * Real.pi *
        ((Real.log x - Real.log rhoMin) / (Real.log rhoMax - Real.log rhoMin))))
      (Real.cos (2 * Real.pi *
          ((Real.log n - Real.log rhoMin) / (Real.log rhoMax - Real.log rhoMin))) *
        (2 * Real.pi * (n⁻¹ / (Real.log rhoMax - Real.log rhoMin)))) n :=
    (Real.hasDerivAt_sin _).comp n h2pit
  have hsin_div : HasDerivAt
      (fun x => Real.sin (2 * Real.pi *
        ((Real.log x - Real.log rhoMin) / (Real.log rhoMax - Real.log rhoMin))) /
        (2 * Real.pi))
      ((Real.cos (2 * Real.pi *
          ((Real.log n - Real.log rhoMin) / (Real.log rhoMax - Real.log rhoMin))) *
        (2 * Real.pi * (n⁻¹ / (Real.log rhoMax - Real.log rhoMin)))) / (2 * Real.pi)) n :=
    hsin.div_const (2 * Real.pi)
  have hS : HasDerivAt (sccsSwitch rhoMin rhoMax)
      (-(n⁻¹ / (Real.log rhoMax - Real.log rhoMin) -
        (Real.cos (2 * Real.pi *
            ((Real.log n - Real.log rhoMin) / (Real.log rhoMax - Real.log rhoMin))) *
          (2 * Real.pi * (n⁻¹ / (Real.log rhoMax - Real.log rhoMin)))) /
          (2 * Real.pi))) n :=
    (ht.sub hsin_div).const_sub 1
  convert hS using 1
  rw [hlog_eq]
  have hn_ne : n ≠ 0 := ne_of_gt hn
  field_simp [hn_ne, hL_ne, h2pi_ne]

/-- The analytic `HasDerivAt` derivative of `sccsSwitch` is everywhere non-positive
on `(rhoMin, rhoMax)`. -/
theorem sccsSwitch_hasDerivAt_nonpos (rhoMin rhoMax n : ℝ)
    (hmin : 0 < rhoMin) (hn_lo : rhoMin < n) (hn_hi : n < rhoMax) :
    ∃ d : ℝ, HasDerivAt (sccsSwitch rhoMin rhoMax) d n ∧ d ≤ 0 := by
  refine ⟨_, sccsSwitch_hasDerivAt rhoMin rhoMax n hmin hn_lo hn_hi, ?_⟩
  exact sccs_deriv_nonpos n rhoMin rhoMax _ (lt_trans hmin hn_lo) hmin (lt_trans hn_lo hn_hi)

/-- VASPsol (Mathew et al. 2014) switching function parameterized by `erfc`:
`S_VASPsol(n) = (1 / 2) * erfc (log (n / n_c) / (sigma * sqrt 2))`. -/
noncomputable def vaspsolSwitch (erfc : ℝ → ℝ) (n_c sigma n : ℝ) : ℝ :=
  (1 / 2 : ℝ) * erfc (Real.log (n / n_c) / (sigma * Real.sqrt 2))

/-- Analytic chain-rule derivative of `vaspsolSwitch` machine-checked via Lean's
`HasDerivAt.comp` through `u(n) = log (n / n_c) / (sigma * sqrt 2)`. -/
theorem vaspsolSwitch_hasDerivAt (erfc : ℝ → ℝ) (n_c sigma n : ℝ)
    (hn_c : 0 < n_c) (hsigma : 0 < sigma) (hn : 0 < n)
    (herfc : HasDerivAt erfc
      (-(2 / Real.sqrt Real.pi) *
        Real.exp (-((Real.log (n / n_c) / (sigma * Real.sqrt 2)) ^ 2)))
      (Real.log (n / n_c) / (sigma * Real.sqrt 2))) :
    HasDerivAt (vaspsolSwitch erfc n_c sigma)
      (-(Real.exp (-((Real.log (n / n_c) / (sigma * Real.sqrt 2)) ^ 2))) /
        (n * sigma * Real.sqrt (2 * Real.pi))) n := by
  have hdiv : HasDerivAt (fun x => x / n_c) (1 / n_c) n :=
    (hasDerivAt_id n).div_const n_c
  have hdiv_ne : n / n_c ≠ 0 := div_ne_zero (ne_of_gt hn) (ne_of_gt hn_c)
  have hlog : HasDerivAt (fun x => Real.log (x / n_c)) ((n / n_c)⁻¹ * (1 / n_c)) n := by
    simpa [Function.comp] using (Real.hasDerivAt_log hdiv_ne).comp n hdiv
  have hu : HasDerivAt (fun x => Real.log (x / n_c) / (sigma * Real.sqrt 2))
      (((n / n_c)⁻¹ * (1 / n_c)) / (sigma * Real.sqrt 2)) n :=
    hlog.div_const (sigma * Real.sqrt 2)
  have hcomp : HasDerivAt
      (fun x => erfc (Real.log (x / n_c) / (sigma * Real.sqrt 2)))
      ((-(2 / Real.sqrt Real.pi) *
          Real.exp (-((Real.log (n / n_c) / (sigma * Real.sqrt 2)) ^ 2))) *
        (((n / n_c)⁻¹ * (1 / n_c)) / (sigma * Real.sqrt 2))) n := by
    simpa [Function.comp] using herfc.comp n hu
  have hS : HasDerivAt (vaspsolSwitch erfc n_c sigma)
      ((1 / 2 : ℝ) *
        ((-(2 / Real.sqrt Real.pi) *
            Real.exp (-((Real.log (n / n_c) / (sigma * Real.sqrt 2)) ^ 2))) *
          (((n / n_c)⁻¹ * (1 / n_c)) / (sigma * Real.sqrt 2)))) n :=
    hcomp.const_mul (1 / 2 : ℝ)
  convert hS using 1
  have hsqrt_mul : Real.sqrt (2 * Real.pi) = Real.sqrt 2 * Real.sqrt Real.pi :=
    Real.sqrt_mul (by norm_num) Real.pi
  have hsqrt2_ne : Real.sqrt 2 ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr (by norm_num))
  have hsqrtpi_ne : Real.sqrt Real.pi ≠ 0 := ne_of_gt (Real.sqrt_pos.mpr Real.pi_pos)
  rw [hsqrt_mul]
  field_simp [ne_of_gt hn, ne_of_gt hn_c, ne_of_gt hsigma, hsqrt2_ne, hsqrtpi_ne]

/-- The analytic chain-rule derivative of `vaspsolSwitch` is everywhere non-positive
for `n > 0`. -/
theorem vaspsolSwitch_hasDerivAt_nonpos (erfc : ℝ → ℝ) (n_c sigma n : ℝ)
    (hn_c : 0 < n_c) (hsigma : 0 < sigma) (hn : 0 < n)
    (herfc : HasDerivAt erfc
      (-(2 / Real.sqrt Real.pi) *
        Real.exp (-((Real.log (n / n_c) / (sigma * Real.sqrt 2)) ^ 2)))
      (Real.log (n / n_c) / (sigma * Real.sqrt 2))) :
    ∃ d : ℝ, HasDerivAt (vaspsolSwitch erfc n_c sigma) d n ∧ d ≤ 0 := by
  refine ⟨_, vaspsolSwitch_hasDerivAt erfc n_c sigma n hn_c hsigma hn herfc, ?_⟩
  have hsqrt_pos : 0 < Real.sqrt (2 * Real.pi) :=
    Real.sqrt_pos.mpr (mul_pos (by norm_num) Real.pi_pos)
  exact vaspsol_deriv_nonpos n sigma (Real.sqrt (2 * Real.pi)) _ hn hsigma hsqrt_pos

/-- Repulsive sign (non-negativity) of the dielectric potential `V_diel` (both
linear-`eps` and SCCS log-`eps` variants) and the electrolyte screening potential
`V_scr` for any monotonically decreasing cavity switching function (`ds_dn ≤ 0`). -/
theorem diel_and_screening_potentials_nonneg
    (eps_b eps kappa ds_dn pi_val grad_phi_sq phi : ℝ)
    (heps_b : 1 ≤ eps_b)
    (heps : 0 ≤ eps)
    (hkappa : 0 ≤ kappa)
    (hds : ds_dn ≤ 0)
    (hpi : 0 < pi_val)
    (hgrad : 0 ≤ grad_phi_sq) :
    (0 ≤ -(1 / (8 * pi_val)) * ((eps_b - 1) * ds_dn) * grad_phi_sq) ∧
    (eps * Real.log eps_b * ds_dn ≤ 0) ∧
    (0 ≤ -(1 / (8 * pi_val)) * (eps * Real.log eps_b * ds_dn) * grad_phi_sq) ∧
    (0 ≤ -(eps_b * kappa ^ 2 / (8 * pi_val)) * ds_dn * phi ^ 2) := by
  have hinv : 0 ≤ 1 / (8 * pi_val) :=
    div_nonneg (by norm_num) (by linarith)
  have hdeps_lin : (eps_b - 1) * ds_dn ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (by linarith) hds
  have hlog : 0 ≤ Real.log eps_b := Real.log_nonneg heps_b
  have hdeps_log : eps * Real.log eps_b * ds_dn ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (mul_nonneg heps hlog) hds
  have hkappa_sq : 0 ≤ kappa ^ 2 := by
    rw [sq]
    exact mul_nonneg hkappa hkappa
  have hfrac : 0 ≤ eps_b * kappa ^ 2 / (8 * pi_val) :=
    div_nonneg (mul_nonneg (by linarith) hkappa_sq) (by linarith)
  refine ⟨?_, hdeps_log, ?_, ?_⟩
  · have hpref : 0 ≤ -(1 / (8 * pi_val)) * ((eps_b - 1) * ds_dn) := by
      have heq : -(1 / (8 * pi_val)) * ((eps_b - 1) * ds_dn) =
          (1 / (8 * pi_val)) * (-((eps_b - 1) * ds_dn)) := by ring
      rw [heq]
      exact mul_nonneg hinv (neg_nonneg.mpr hdeps_lin)
    exact mul_nonneg hpref hgrad
  · have hpref : 0 ≤ -(1 / (8 * pi_val)) * (eps * Real.log eps_b * ds_dn) := by
      have heq : -(1 / (8 * pi_val)) * (eps * Real.log eps_b * ds_dn) =
          (1 / (8 * pi_val)) * (-(eps * Real.log eps_b * ds_dn)) := by ring
      rw [heq]
      exact mul_nonneg hinv (neg_nonneg.mpr hdeps_log)
    exact mul_nonneg hpref hgrad
  · have hpref : 0 ≤ -(eps_b * kappa ^ 2 / (8 * pi_val)) * ds_dn := by
      have heq : -(eps_b * kappa ^ 2 / (8 * pi_val)) * ds_dn =
          (eps_b * kappa ^ 2 / (8 * pi_val)) * (-ds_dn) := by ring
      rw [heq]
      exact mul_nonneg hfrac (neg_nonneg.mpr hds)
    exact mul_nonneg hpref (sq_nonneg phi)

/-!
## Section 6: Gate 3–4 SCF Invariants (Sign Involution, Spin Chain Rule, Bregman Bound)
-/

/-- Electron potential energy sign convention (ADR 0017): `U_el(phi) = -phi`. -/
def electronPotential {n : ℕ} (phi : Fin n → ℝ) (i : Fin n) : ℝ :=
  -phi i

/-- ADR 0017 electrostatic sign involution: `U_el(phi) = -phi` is an involution
preserving both the pointwise energy density `(-n) * phi = n * U_el(phi)` and the
weighted grid inner product `⟨-n, phi⟩ = ⟨n, U_el(phi)⟩`. -/
theorem electrostatic_sign_involution {n : ℕ} (dv : ℝ) (rho_e phi : Fin n → ℝ) :
    (∀ i, electronPotential (electronPotential phi) i = phi i) ∧
    (∀ i, (-rho_e i) * phi i = rho_e i * electronPotential phi i) ∧
    wdot dv (fun i => -rho_e i) phi = wdot dv rho_e (electronPotential phi) := by
  refine ⟨?_, ?_, ?_⟩
  · intro i
    unfold electronPotential
    ring
  · intro i
    unfold electronPotential
    ring
  · unfold wdot electronPotential
    congr 1
    apply Finset.sum_congr rfl
    intro i _
    ring

/-- Spin-polarized solvation Gateaux chain rule: for any Gateaux-differentiable
total-density functional `E[n_up + n_dn]` with derivative `V_solv[n_tot]`, the partial
directional derivatives w.r.t. `n_up` and `n_dn` are identical (`V_solv[n_tot]`), the
joint directional derivative is `⟨V, dn_up⟩ + ⟨V, dn_dn⟩`, and the spin-summed
double-counting equals `⟨n_up + n_dn, V⟩`. -/
theorem spin_polarized_solvation_chain_rule {n : ℕ} (dv : ℝ)
    (E : (Fin n → ℝ) → ℝ) (n_up n_dn V : Fin n → ℝ)
    (hgateaux : ∀ dn : Fin n → ℝ,
      HasDerivAt (fun t : ℝ => E (fun i => (n_up i + n_dn i) + t * dn i))
        (wdot dv V dn) 0) :
    (∀ dn : Fin n → ℝ,
      HasDerivAt (fun t : ℝ => E (fun i => (n_up i + t * dn i) + n_dn i))
        (wdot dv V dn) 0) ∧
    (∀ dn : Fin n → ℝ,
      HasDerivAt (fun t : ℝ => E (fun i => n_up i + (n_dn i + t * dn i)))
        (wdot dv V dn) 0) ∧
    (∀ dn_up dn_dn : Fin n → ℝ,
      HasDerivAt
        (fun t : ℝ => E (fun i => (n_up i + t * dn_up i) + (n_dn i + t * dn_dn i)))
        (wdot dv V dn_up + wdot dv V dn_dn) 0) ∧
    wdot dv n_up V + wdot dv n_dn V = wdot dv (fun i => n_up i + n_dn i) V := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro dn
    have hfun : (fun t : ℝ => E (fun i => (n_up i + t * dn i) + n_dn i)) =
        (fun t : ℝ => E (fun i => (n_up i + n_dn i) + t * dn i)) := by
      funext t
      congr 1
      funext i
      ring
    rw [hfun]
    exact hgateaux dn
  · intro dn
    have hfun : (fun t : ℝ => E (fun i => n_up i + (n_dn i + t * dn i))) =
        (fun t : ℝ => E (fun i => (n_up i + n_dn i) + t * dn i)) := by
      funext t
      congr 1
      funext i
      ring
    rw [hfun]
    exact hgateaux dn
  · intro dn_up dn_dn
    have hfun :
        (fun t : ℝ => E (fun i => (n_up i + t * dn_up i) + (n_dn i + t * dn_dn i))) =
        (fun t : ℝ => E (fun i => (n_up i + n_dn i) + t * (dn_up i + dn_dn i))) := by
      funext t
      congr 1
      funext i
      ring
    rw [hfun, ← wdot_add_right]
    exact hgateaux (fun i => dn_up i + dn_dn i)
  · exact (wdot_add_left dv n_up n_dn V).symm

/-- Non-negative triangle inequality bound for `EnergyAudit::functional_difference_bound`:
`0 ≤ |R_H| + |R_xc| + |R_solv|` and `|R_H + R_xc + R_solv| ≤ |R_H| + |R_xc| + |R_solv|`. -/
theorem bregman_bound_nonneg_triangle (R_H R_xc R_solv : ℝ) :
    0 ≤ |R_H| + |R_xc| + |R_solv| ∧
    |R_H + R_xc + R_solv| ≤ |R_H| + |R_xc| + |R_solv| := by
  refine ⟨?_, ?_⟩
  · linarith [abs_nonneg R_H, abs_nonneg R_xc, abs_nonneg R_solv]
  · calc
      |R_H + R_xc + R_solv|
        ≤ |R_H + R_xc| + |R_solv| := abs_add_le (R_H + R_xc) R_solv
      _ ≤ (|R_H| + |R_xc|) + |R_solv| := by linarith [abs_add_le R_H R_xc]
      _ = |R_H| + |R_xc| + |R_solv| := by ring

/-!
## Proof Axiom Cleanliness Audit
-/

#print axioms gpOp_quad_form
#print axioms gpOp_symmetric
#print axioms gpOp_psd
#print axioms gpOp_spd
#print axioms dualEnergy_on_shell
#print axioms dualEnergy_perturb_exact
#print axioms dualEnergy_global_max
#print axioms gpOp_vacuum_eq
#print axioms gpOp_apply_zero
#print axioms polarization_rhs_vanishes_in_vacuum
#print axioms zero_double_counting_cancellation
#print axioms width_invariance_arbitrary_charge
#print axioms harris_foulkes_bregman_exact
#print axioms grand_canonical_bregman_exact
#print axioms vaspsol_deriv_nonpos
#print axioms sccs_deriv_nonpos
#print axioms sccsSwitch_hasDerivAt
#print axioms sccsSwitch_hasDerivAt_nonpos
#print axioms vaspsolSwitch_hasDerivAt
#print axioms vaspsolSwitch_hasDerivAt_nonpos
#print axioms diel_and_screening_potentials_nonneg
#print axioms electrostatic_sign_involution
#print axioms spin_polarized_solvation_chain_rule
#print axioms bregman_bound_nonneg_triangle

end BaandiMath.Solvation
