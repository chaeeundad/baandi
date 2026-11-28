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

Machine-checked Lean 4 + Mathlib proofs of the core Level-1 theoretical guarantees in
`docs/research/extensible-dft/implicit-solvation-theory.md` (Section 4) and Gate 3–4
SCF invariants (mapped hypothesis-by-hypothesis to Rust in
`docs/research/extensible-dft/formal-verification-map.md`):

1. Weighted grid inner product, Nyquist-zeroed Fourier-multiplier skew-adjointness
   (`D_a^T = -D_a` and real-field closure), exact kernel characterization of the
   discrete GP/LPB operator `A = ∑_a D_a^T (eps * D_a) + k2`, and strict positive
   definiteness (SPD) on the actual solver computation spaces:
   - Pure GP (`k2 = 0`) on the `ker(D)`-projected subspace (`zero_g2_indices`),
   - Realistic cavity-vanishing LPB (`0 ≤ k2`, allowing `k2 = 0` inside the solute),
   - Single-bulk-point LPB positivity on constant modes (`G = 0` gauge pinning).
2. Dual variational electrostatic energy functional: on-shell identity, exact quadratic
   perturbation identity, global maximality (envelope theorem / stationarity), and exact
   vacuum reduction (`phi_pol = 0`, `Delta E_es = 0`).
3. Exact zero double counting, algebraic `G = 0` background-offset cancellation, and
   genuine PDE-level full electrostatic functional Gaussian smearing-width invariance
   (`(1/2) ⟨rho, phi_pol⟩ + Q_net * alpha_sm`) for both neutral and charged cells.
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
## Section 1: Weighted Grid Inner Product, Nyquist-Zeroed Skew-Adjointness & Operator SPD
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

theorem wdot_sub_right {n : ℕ} (dv : ℝ) (u v w : Fin n → ℝ) :
    wdot dv u (fun i => v i - w i) = wdot dv u v - wdot dv u w := by
  have h : (fun i => v i - w i) = fun i => v i + (-w i) := by
    funext i
    ring
  rw [h, wdot_add_right, wdot_neg_right]
  ring

theorem wdot_sub_left {n : ℕ} (dv : ℝ) (u v w : Fin n → ℝ) :
    wdot dv (fun i => u i - v i) w = wdot dv u w - wdot dv v w := by
  rw [wdot_comm, wdot_sub_right, wdot_comm dv w u, wdot_comm dv w v]

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

theorem wdot_smul_left {n : ℕ} (dv c : ℝ) (u v : Fin n → ℝ) :
    wdot dv (fun i => c * u i) v = c * wdot dv u v := by
  rw [wdot_comm, wdot_smul_right, wdot_comm]

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

theorem wdot_diag_eq_zero_iff {n : ℕ} {dv : ℝ} (hdv : 0 < dv) {w u : Fin n → ℝ}
    (hw : ∀ i, 0 ≤ w i) :
    wdot dv (fun i => w i * u i) u = 0 ↔ ∀ i, w i = 0 ∨ u i = 0 := by
  unfold wdot
  constructor
  · intro h i
    have hsum : ∑ j : Fin n, (w j * u j) * u j = 0 :=
      (mul_eq_zero.mp h).resolve_left (ne_of_gt hdv)
    have hnonneg : ∀ j ∈ (Finset.univ : Finset (Fin n)), 0 ≤ (w j * u j) * u j := by
      intro j _
      have hterm : (w j * u j) * u j = w j * (u j * u j) := by ring
      rw [hterm]
      exact mul_nonneg (hw j) (mul_self_nonneg (u j))
    have hle : (w i * u i) * u i ≤ ∑ j : Fin n, (w j * u j) * u j :=
      Finset.single_le_sum hnonneg (Finset.mem_univ i)
    have hge : 0 ≤ (w i * u i) * u i := hnonneg i (Finset.mem_univ i)
    have hi0 : (w i * u i) * u i = 0 := by linarith
    have hprod : w i * (u i * u i) = 0 := by
      have hterm : (w i * u i) * u i = w i * (u i * u i) := by ring
      rw [← hterm]
      exact hi0
    rcases mul_eq_zero.mp hprod with hw0 | huu0
    · exact Or.inl hw0
    · rcases mul_eq_zero.mp huu0 with hu0 | hu0
      · exact Or.inr hu0
      · exact Or.inr hu0
  · intro h
    have hsum : ∑ i : Fin n, (w i * u i) * u i = 0 := by
      calc
        ∑ i : Fin n, (w i * u i) * u i
          = ∑ _i : Fin n, (0 : ℝ) := by
              apply Finset.sum_congr rfl
              intro i _
              rcases h i with hw0 | hu0
              · rw [hw0]
                ring
              · rw [hu0]
                ring
        _ = 0 := Finset.sum_const_zero
    rw [hsum, mul_zero]

theorem wdot_diag_pos_of_exists {n : ℕ} {dv : ℝ} (hdv : 0 < dv) {w u : Fin n → ℝ}
    (hw : ∀ i, 0 ≤ w i) (hpos : ∃ i, 0 < w i ∧ u i ≠ 0) :
    0 < wdot dv (fun i => w i * u i) u := by
  rcases hpos with ⟨i0, hw0, hu0⟩
  unfold wdot
  apply mul_pos hdv
  have hnonneg : ∀ i ∈ (Finset.univ : Finset (Fin n)), 0 ≤ (w i * u i) * u i := by
    intro i _
    have hterm : (w i * u i) * u i = w i * (u i * u i) := by ring
    rw [hterm]
    exact mul_nonneg (hw i) (mul_self_nonneg (u i))
  have hsq_pos : 0 < u i0 * u i0 := by
    rcases lt_trichotomy (u i0) 0 with hneg | hzero | hpos_u
    · nlinarith
    · exact absurd hzero hu0
    · nlinarith
  have hpos_i0 : 0 < (w i0 * u i0) * u i0 := by
    have hterm : (w i0 * u i0) * u i0 = w i0 * (u i0 * u i0) := by ring
    rw [hterm]
    exact mul_pos hw0 hsq_pos
  have hle : (w i0 * u i0) * u i0 ≤ ∑ i : Fin n, (w i * u i) * u i :=
    Finset.single_le_sum hnonneg (Finset.mem_univ i0)
  exact lt_of_lt_of_le hpos_i0 hle

theorem wdot_diag_pos {n : ℕ} {dv : ℝ} (hdv : 0 < dv) {w u : Fin n → ℝ}
    (hw : ∀ i, 0 < w i) (hu : ∃ i, u i ≠ 0) :
    0 < wdot dv (fun i => w i * u i) u := by
  rcases hu with ⟨i0, hi0⟩
  exact wdot_diag_pos_of_exists hdv (fun i => le_of_lt (hw i)) ⟨i0, hw i0, hi0⟩

/-- Parseval real inner product of complex Fourier-space fields `(ur + i ui)` and
`(vr + i vi)` representing real-space fields (`Re ⟨û, v̂⟩_Ω = Ω ∑_k (u_r v_r + u_i v_i)`). -/
def fourierDot {n : ℕ} (vol : ℝ) (ur ui vr vi : Fin n → ℝ) : ℝ :=
  vol * ∑ k : Fin n, (ur k * vr k + ui k * vi k)

/-- Mathematical characterization of BAANDI's Nyquist-zeroed spectral derivative
(`cavity::skew_adjoint_g_cart`): for the reciprocal-lattice mirror map `mirror : Fin n → Fin n`
(`-G mod N`), multiplying Fourier coefficients by `i * g_k` (`(ur, ui) ↦ (-g * ui, g * ur)`)
satisfies:
1. **Self-conjugate (Nyquist & `G = 0`) zeroing**: oddness `g (mirror k) = -g k` forces
   `g k = 0` on every self-conjugate mode `mirror k = k` (i.e., `G = 0` and all even-axis
   Nyquist planes `2 * ijk[a] == N[a]`).
2. **Real-field closure**: whenever `g` is odd under `mirror` and `(ur, ui)` is
   conjugate-symmetric (`ur (mirror k) = ur k`, `ui (mirror k) = -ui k`), the derivative
   `(-g * ui, g * ur)` is again conjugate-symmetric (maps real grids to real grids).
3. **Exact Parseval skew-adjointness**: `⟨u, D v⟩_Ω = -⟨D u, v⟩_Ω`. -/
theorem nyquist_zeroed_multiplier_skew_adjoint {n : ℕ} (vol : ℝ)
    (mirror : Fin n → Fin n)
    (g ur ui vr vi : Fin n → ℝ)
    (hg_odd : ∀ k, g (mirror k) = -g k)
    (hur : ∀ k, ur (mirror k) = ur k)
    (hui : ∀ k, ui (mirror k) = -ui k) :
    (∀ k, mirror k = k → g k = 0) ∧
    (∀ k, -g (mirror k) * ui (mirror k) = -g k * ui k) ∧
    (∀ k, g (mirror k) * ur (mirror k) = -(g k * ur k)) ∧
    fourierDot vol ur ui (fun k => -g k * vi k) (fun k => g k * vr k) =
      -fourierDot vol (fun k => -g k * ui k) (fun k => g k * ur k) vr vi := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro k hfix
    have h := hg_odd k
    rw [hfix] at h
    linarith
  · intro k
    rw [hg_odd k, hui k]
    ring
  · intro k
    rw [hg_odd k, hur k]
    ring
  · unfold fourierDot
    have hsum :
        ∑ k : Fin n, (ur k * (-g k * vi k) + ui k * (g k * vr k)) =
          -(∑ k : Fin n, ((-g k * ui k) * vr k + (g k * ur k) * vi k)) := by
      rw [← Finset.sum_neg_distrib]
      apply Finset.sum_congr rfl
      intro k _
      ring
    rw [hsum]
    ring

/-- Converse necessity of Nyquist zeroing: if a self-conjugate mode (`mirror k = k`,
such as an even-axis Nyquist plane) keeps a nonzero multiplier `g_k ≠ 0`, then a nonzero
real amplitude `ur_k ≠ 0` produces a nonzero imaginary component `g_k * ur_k ≠ -(g_k * ur_k)`,
violating Hermitian conjugate symmetry and breaking real-field closure of `D_a`. -/
theorem nyquist_nonzero_breaks_reality (g_k ur_k : ℝ)
    (hg : g_k ≠ 0) (hu : ur_k ≠ 0) :
    g_k * ur_k ≠ -(g_k * ur_k) := by
  intro h
  have hzero : g_k * ur_k = 0 := by linarith
  rcases mul_eq_zero.mp hzero with hg0 | hu0
  · exact hg hg0
  · exact hu hu0

/-- Discrete Generalized Poisson / Linearized Poisson–Boltzmann operator
`A u = ∑_{a=0}^2 D_a^T (eps * D_a u) + k2 * u` with skew-adjoint `D_a^T = -D_a`. -/
def gpOp {n : ℕ} (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (eps k2 u : Fin n → ℝ) (i : Fin n) : ℝ :=
  (∑ a : Fin 3, -(D a (fun j => eps j * D a u j) i)) + k2 i * u i

/-- Pointwise additivity of `gpOp` when each discrete gradient `D_a` is additive. -/
theorem gpOp_add {n : ℕ}
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hD_add : ∀ a u v i, D a (fun j => u j + v j) i = D a u i + D a v i)
    (eps k2 u v : Fin n → ℝ) (i : Fin n) :
    gpOp D eps k2 (fun j => u j + v j) i = gpOp D eps k2 u i + gpOp D eps k2 v i := by
  unfold gpOp
  have hsum :
      ∑ a : Fin 3, -(D a (fun j => eps j * D a (fun l => u l + v l) j) i) =
        (∑ a : Fin 3, -(D a (fun j => eps j * D a u j) i)) +
        (∑ a : Fin 3, -(D a (fun j => eps j * D a v j) i)) := by
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro a _
    have hinner :
        (fun j => eps j * D a (fun l => u l + v l) j) =
        (fun j => eps j * D a u j + eps j * D a v j) := by
      funext j
      rw [hD_add a u v j]
      ring
    rw [hinner, hD_add]
    ring
  rw [hsum]
  ring

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

/-- Exact kernel characterization of the GP/LPB quadratic form under realistic physical
bounds (`dv > 0`, `eps ≥ 1`, `k2 ≥ 0` allowing `k2 = 0` inside the solute cavity):
`⟨u, A u⟩ = 0` iff `u` lies in `ker(D)` AND vanishes on the electrolyte support `{i | k2 i > 0}`. -/
theorem gpOp_quad_eq_zero_iff {n : ℕ} (dv : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 u : Fin n → ℝ)
    (hdv : 0 < dv)
    (heps : ∀ i, 1 ≤ eps i)
    (hk2 : ∀ i, 0 ≤ k2 i) :
    wdot dv u (gpOp D eps k2 u) = 0 ↔
      (∀ a i, D a u i = 0) ∧ (∀ i, k2 i = 0 ∨ u i = 0) := by
  rw [gpOp_quad_form dv D hskew eps k2 u u]
  have heps_nonneg : ∀ i, 0 ≤ eps i := fun i => by linarith [heps i]
  have hgrad_term_nonneg :
      ∀ a ∈ (Finset.univ : Finset (Fin 3)),
        0 ≤ wdot dv (fun i => eps i * D a u i) (D a u) := by
    intro a _
    exact wdot_diag_nonneg (le_of_lt hdv) heps_nonneg
  have hgrad_nonneg :
      0 ≤ ∑ a : Fin 3, wdot dv (fun i => eps i * D a u i) (D a u) :=
    Finset.sum_nonneg hgrad_term_nonneg
  have hscr_nonneg : 0 ≤ wdot dv (fun i => k2 i * u i) u :=
    wdot_diag_nonneg (le_of_lt hdv) hk2
  constructor
  · intro hzero
    have hgrad_zero :
        ∑ a : Fin 3, wdot dv (fun i => eps i * D a u i) (D a u) = 0 := by linarith
    have hscr_zero : wdot dv (fun i => k2 i * u i) u = 0 := by linarith
    refine ⟨?_, (wdot_diag_eq_zero_iff hdv hk2).mp hscr_zero⟩
    intro a i
    have hle :
        wdot dv (fun j => eps j * D a u j) (D a u) ≤
          ∑ b : Fin 3, wdot dv (fun j => eps j * D b u j) (D b u) :=
      Finset.single_le_sum hgrad_term_nonneg (Finset.mem_univ a)
    have hge : 0 ≤ wdot dv (fun j => eps j * D a u j) (D a u) :=
      hgrad_term_nonneg a (Finset.mem_univ a)
    have ha_zero : wdot dv (fun j => eps j * D a u j) (D a u) = 0 := by linarith
    rcases (wdot_diag_eq_zero_iff hdv heps_nonneg).mp ha_zero i with heps0 | hdu0
    · linarith [heps i]
    · exact hdu0
  · rintro ⟨hDu, hk2u⟩
    have hgrad_zero :
        ∑ a : Fin 3, wdot dv (fun i => eps i * D a u i) (D a u) = 0 := by
      calc
        ∑ a : Fin 3, wdot dv (fun i => eps i * D a u i) (D a u)
          = ∑ _a : Fin 3, (0 : ℝ) := by
              apply Finset.sum_congr rfl
              intro a _
              exact (wdot_diag_eq_zero_iff hdv heps_nonneg).mpr
                (fun i => Or.inr (hDu a i))
        _ = 0 := Finset.sum_const_zero
    have hscr_zero : wdot dv (fun i => k2 i * u i) u = 0 :=
      (wdot_diag_eq_zero_iff hdv hk2).mpr hk2u
    linarith

/-- Pure-GP (`k2 = 0` everywhere) PCG gauge-projected subspace SPD:
for any `u` outside `ker(D)` (`∃ a i, D a u i ≠ 0`), `0 < ⟨u, A u⟩` holds under
`dv > 0`, `eps ≥ 1`, and `k2 ≥ 0` (including `k2 ≡ 0`). This matches BAANDI's
`zero_g2_indices` null-space projection in `poisson.rs`, which zeroes `G = 0` and all
even-axis Nyquist corner modes where every component of `skew_adjoint_g_cart` vanishes. -/
theorem gpOp_spd_on_grad_complement {n : ℕ} (dv : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 u : Fin n → ℝ)
    (hdv : 0 < dv)
    (heps : ∀ i, 1 ≤ eps i)
    (hk2 : ∀ i, 0 ≤ k2 i)
    (hu_grad : ∃ a i, D a u i ≠ 0) :
    0 < wdot dv u (gpOp D eps k2 u) := by
  have hpsd := gpOp_psd dv D hskew eps k2 u hdv heps hk2
  rcases lt_trichotomy (wdot dv u (gpOp D eps k2 u)) 0 with hneg | hzero | hpos
  · linarith
  · exfalso
    have hker := (gpOp_quad_eq_zero_iff dv D hskew eps k2 u hdv heps hk2).mp hzero
    rcases hu_grad with ⟨a, i, hne⟩
    exact hne (hker.1 a i)
  · exact hpos

/-- Realistic LPB (`κ > 0`) cavity-vanishing SPD:
with `∀ i, 0 ≤ k2 i` (explicitly allowing `k2 i = 0` inside the solute cavity!),
`0 < ⟨u, A u⟩` holds whenever `u` either has a nonzero discrete gradient
(`∃ a i, D a u i ≠ 0`) or is nonzero at some electrolyte grid point
(`∃ i, 0 < k2 i ∧ u i ≠ 0`).

Discrete Nyquist parity-coset note: on an FFT grid with even axes, `ker(D)` is spanned
by the constant mode `G = 0` plus the up-to-7 even-Nyquist corner parity modes
`(-1)^(p · ijk)` (`p ∈ {0, 1}^3 \ {0}`). Consequently, `gpOp` is strictly SPD on the
entire real-space grid iff the bulk electrolyte support `{i | 0 < k2 i}` intersects all
`2^(# even axes)` parity cosets (automatically satisfied whenever the bulk solvent
contains any contiguous `2 × 2 × 2` block of grid points). -/
theorem gpOp_spd_lpb_cavity_vanishing {n : ℕ} (dv : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 u : Fin n → ℝ)
    (hdv : 0 < dv)
    (heps : ∀ i, 1 ≤ eps i)
    (hk2 : ∀ i, 0 ≤ k2 i)
    (hu : (∃ a i, D a u i ≠ 0) ∨ (∃ i, 0 < k2 i ∧ u i ≠ 0)) :
    0 < wdot dv u (gpOp D eps k2 u) := by
  have hpsd := gpOp_psd dv D hskew eps k2 u hdv heps hk2
  rcases lt_trichotomy (wdot dv u (gpOp D eps k2 u)) 0 with hneg | hzero | hpos
  · linarith
  · exfalso
    have hker := (gpOp_quad_eq_zero_iff dv D hskew eps k2 u hdv heps hk2).mp hzero
    rcases hu with ⟨a, i, hne⟩ | ⟨i, hkpos, hne⟩
    · exact hne (hker.1 a i)
    · rcases hker.2 i with hk0 | hu0
      · linarith
      · exact hne hu0
  · exact hpos

/-- Constant-mode LPB positivity: for any nonzero constant potential `u = fun _ => c`
(`c ≠ 0`), strict positivity `0 < ⟨u, A u⟩` requires ONLY a single bulk electrolyte
point `∃ i0, 0 < k2 i0` (not pointwise positivity everywhere). This pins the `G = 0`
gauge (`φ_bulk = 0`) for charged cells and Phase B2 GC-DFT. -/
theorem gpOp_spd_const_mode_lpb {n : ℕ} (dv c : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 : Fin n → ℝ)
    (hdv : 0 < dv)
    (heps : ∀ i, 1 ≤ eps i)
    (hk2 : ∀ i, 0 ≤ k2 i)
    (hbulk : ∃ i0, 0 < k2 i0)
    (hc : c ≠ 0) :
    0 < wdot dv (fun _ => c) (gpOp D eps k2 (fun _ => c)) := by
  rcases hbulk with ⟨i0, hi0⟩
  exact gpOp_spd_lpb_cavity_vanishing dv D hskew eps k2 (fun _ => c) hdv heps hk2
    (Or.inr ⟨i0, hi0, hc⟩)

/-- Restricted corollary: strict positive definiteness of `gpOp` under strictly positive
uniform/everywhere screening (`∀ i, 0 < k2 i`, e.g., analytic MMS benchmarks).
For physical cavities where `k2 = 0` inside the solute, use `gpOp_spd_lpb_cavity_vanishing`
(LPB) or `gpOp_spd_on_grad_complement` (pure GP). -/
theorem gpOp_spd_uniform_screening {n : ℕ} (dv : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 u : Fin n → ℝ)
    (hdv : 0 < dv)
    (heps : ∀ i, 1 ≤ eps i)
    (hk2 : ∀ i, 0 < k2 i)
    (hu : ∃ i, u i ≠ 0) :
    0 < wdot dv u (gpOp D eps k2 u) := by
  rcases hu with ⟨i0, hi0⟩
  exact gpOp_spd_lpb_cavity_vanishing dv D hskew eps k2 u hdv heps
    (fun i => le_of_lt (hk2 i)) (Or.inr ⟨i0, hk2 i0, hi0⟩)

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

/-- Exact symmetry of the discrete vacuum Laplacian `vacuumOp` under skew-adjoint `D_a`. -/
theorem vacuumOp_symmetric {n : ℕ} (dv : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (u v : Fin n → ℝ) :
    wdot dv u (vacuumOp D v) = wdot dv (vacuumOp D u) v := by
  have hu : vacuumOp D u = gpOp D (fun _ => 1) (fun _ => 0) u :=
    funext (fun i => (gpOp_vacuum_eq D u i).symm)
  have hv : vacuumOp D v = gpOp D (fun _ => 1) (fun _ => 0) v :=
    funext (fun i => (gpOp_vacuum_eq D v i).symm)
  rw [hu, hv]
  exact gpOp_symmetric dv D hskew (fun _ => 1) (fun _ => 0) u v

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
## Section 3: Zero Double Counting, Background-Offset Cancellation & PDE Width Invariance
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

/-- Algebraic `G = 0` background-offset cancellation ONLY: for a fixed `E_es_solv`,
the Gaussian background term `(N_e - Q_ion) * alpha_sm` in the bare decomposition cancels
against `E_bg,corr = (Q_ion - N_e) * alpha_sm` identically for any `alpha_sm1, alpha_sm2`.

Scope note: this algebraic identity alone does NOT prove width invariance of the coupled
solvation PDE solution `phi_pol(r_c)` or `Delta E_es(r_c)`. For the genuine PDE-level
electrostatic functional width-invariance theorem under non-overlapping cavity support
(Newton shell theorem), see `full_electrostatic_width_invariance` below. -/
theorem bg_offset_algebraic_cancellation
    (E_es_vac S P E_loc_short E_es_solv N_e Q_ion alpha_sm1 alpha_sm2 : ℝ) :
    ((E_es_vac - (S - P) + E_loc_short + (N_e - Q_ion) * alpha_sm1) +
      (E_es_solv - E_es_vac) + (Q_ion - N_e) * alpha_sm1) =
    ((E_es_vac - (S - P) + E_loc_short + (N_e - Q_ion) * alpha_sm2) +
      (E_es_solv - E_es_vac) + (Q_ion - N_e) * alpha_sm2) := by
  ring

/-- Genuine PDE-level full electrostatic functional Gaussian-width invariance theorem
(arbitrary net charge `Q_net = ⟨rho1, 1⟩ = Q_ion - N_e`, neutral or charged):

Let `rho1, rho2` be two solute charge densities differing only by the Gaussian ionic
pseudo-charge smearing width (`alpha_sm1` vs `alpha_sm2`), with equal mean charge density
`rho_bar` and zero-mean vacuum potentials `phi_vac1, phi_vac2` solving
`L0 phi_vac_k = 4 * pi * (rho_k - rho_bar)`.
Define `d_alpha = alpha_sm2 - alpha_sm1` and the shifted vacuum difference
`phi_tilde = (phi_vac2 - phi_vac1) - d_alpha`.

Under the physical Newton-shell non-overlap hypotheses:
- `h_supp_eps`: `D_a phi_tilde` vanishes wherever `eps ≠ 1` (Gaussian difference field
  stays strictly inside the `eps = 1` solute cavity),
- `h_supp_k2`: `phi_tilde` vanishes wherever `k2 ≠ 0` (shifted Gaussian potential
  difference vanishes everywhere in the mobile electrolyte),

and given `phi_pol1` solving the polarization PDE `A phi_pol1 = 4 * pi * rho1 - A phi_vac1`:
1. Charge conservation holds: `⟨rho2, 1⟩ = ⟨rho1, 1⟩ = Q_net`.
2. The uniformly shifted polarization potential `phi_pol2 = phi_pol1 - d_alpha` solves
   the exact polarization PDE `A phi_pol2 = 4 * pi * rho2 - A phi_vac2`.
3. By Green's reciprocity for `A` and `L0`, the full reconciled solvation electrostatic
   energy `(1/2) ⟨rho, phi_pol⟩ + Q_net * alpha_sm` is **strictly invariant** under the
   smearing-width change `alpha_sm1 ↦ alpha_sm2` (with half of `-Q_net * d_alpha` arising
   from the `-d_alpha` shift of `phi_pol2` and half from `⟨rho2 - rho1, phi_pol1⟩`,
   cancelling `+Q_net * d_alpha` from `E_bg,corr` identically). -/
theorem full_electrostatic_width_invariance {n : ℕ} (dv pi_val : ℝ)
    (D : Fin 3 → (Fin n → ℝ) → (Fin n → ℝ))
    (hD_add : ∀ a u v i, D a (fun j => u j + v j) i = D a u i + D a v i)
    (hD_sub : ∀ a u v i, D a (fun j => u j - v j) i = D a u i - D a v i)
    (hD_const : ∀ a c i, D a (fun _ => c) i = 0)
    (hskew : ∀ a u v, wdot dv u (D a v) = -wdot dv (D a u) v)
    (eps k2 rho1 rho2 phi_vac1 phi_vac2 phi_pol1 : Fin n → ℝ)
    (rho_bar alpha_sm1 alpha_sm2 : ℝ)
    (hpi : pi_val ≠ 0)
    (h_vac1 : ∀ i, vacuumOp D phi_vac1 i = 4 * pi_val * (rho1 i - rho_bar))
    (h_vac2 : ∀ i, vacuumOp D phi_vac2 i = 4 * pi_val * (rho2 i - rho_bar))
    (h_mean1 : wdot dv phi_vac1 (fun _ => 1) = 0)
    (h_mean2 : wdot dv phi_vac2 (fun _ => 1) = 0)
    (h_supp_eps : ∀ a i, eps i ≠ 1 →
      D a (fun j => (phi_vac2 j - phi_vac1 j) - (alpha_sm2 - alpha_sm1)) i = 0)
    (h_supp_k2 : ∀ i, k2 i ≠ 0 →
      (phi_vac2 i - phi_vac1 i) - (alpha_sm2 - alpha_sm1) = 0)
    (h_pol1 : ∀ i, gpOp D eps k2 phi_pol1 i = 4 * pi_val * rho1 i - gpOp D eps k2 phi_vac1 i) :
    let d_alpha := alpha_sm2 - alpha_sm1
    let phi_pol2 : Fin n → ℝ := fun i => phi_pol1 i - d_alpha
    let Q_net := wdot dv rho1 (fun _ => 1)
    (wdot dv rho2 (fun _ => 1) = Q_net) ∧
    (∀ i, gpOp D eps k2 phi_pol2 i = 4 * pi_val * rho2 i - gpOp D eps k2 phi_vac2 i) ∧
    ((1 / 2 : ℝ) * wdot dv rho2 phi_pol2 + Q_net * alpha_sm2 =
      (1 / 2 : ℝ) * wdot dv rho1 phi_pol1 + Q_net * alpha_sm1) := by
  dsimp only
  let phi_tilde : Fin n → ℝ :=
    fun i => (phi_vac2 i - phi_vac1 i) - (alpha_sm2 - alpha_sm1)
  have h4pi : 4 * pi_val ≠ 0 := mul_ne_zero (by norm_num) hpi
  have h_vac_tilde : ∀ i, vacuumOp D phi_tilde i = 4 * pi_val * (rho2 i - rho1 i) := by
    intro i
    have hsum :
        vacuumOp D phi_tilde i = vacuumOp D phi_vac2 i - vacuumOp D phi_vac1 i := by
      unfold vacuumOp
      rw [← Finset.sum_sub_distrib]
      apply Finset.sum_congr rfl
      intro a _
      have h1 : D a phi_tilde = fun j => D a phi_vac2 j - D a phi_vac1 j := by
        funext j
        have hstep1 :
            D a phi_tilde j =
              D a (fun l => phi_vac2 l - phi_vac1 l) j -
                D a (fun _ => alpha_sm2 - alpha_sm1) j :=
          hD_sub a (fun l => phi_vac2 l - phi_vac1 l) (fun _ => alpha_sm2 - alpha_sm1) j
        rw [hstep1, hD_sub a phi_vac2 phi_vac1 j, hD_const a (alpha_sm2 - alpha_sm1) j]
        ring
      rw [h1, hD_sub a (D a phi_vac2) (D a phi_vac1) i]
      ring
    rw [hsum, h_vac2 i, h_vac1 i]
    ring
  have h_gp_tilde : ∀ i, gpOp D eps k2 phi_tilde i = 4 * pi_val * (rho2 i - rho1 i) := by
    intro i
    have heq : gpOp D eps k2 phi_tilde i = vacuumOp D phi_tilde i := by
      unfold gpOp vacuumOp
      have hk0 : k2 i * phi_tilde i = 0 := by
        by_cases hk : k2 i = 0
        · rw [hk, zero_mul]
        · have hpt : phi_tilde i = 0 := h_supp_k2 i hk
          rw [hpt, mul_zero]
      rw [hk0, add_zero]
      apply Finset.sum_congr rfl
      intro a _
      have heps_eq : (fun j => eps j * D a phi_tilde j) = D a phi_tilde := by
        funext j
        by_cases he : eps j = 1
        · rw [he, one_mul]
        · have hdt : D a phi_tilde j = 0 := h_supp_eps a j he
          rw [hdt, mul_zero]
      rw [heps_eq]
    rw [heq, h_vac_tilde i]
  have h_vac_one : ∀ i, vacuumOp D (fun _ => (1 : ℝ)) i = 0 := by
    intro i
    unfold vacuumOp
    have hterm : ∀ a : Fin 3, -(D a (D a (fun _ => (1 : ℝ))) i) = 0 := by
      intro a
      have h1 : D a (fun _ => (1 : ℝ)) = fun _ => (0 : ℝ) := by
        funext j
        exact hD_const a 1 j
      rw [h1, hD_const a 0 i, neg_zero]
    calc
      ∑ a : Fin 3, -(D a (D a (fun _ => (1 : ℝ))) i)
        = ∑ _a : Fin 3, (0 : ℝ) := by
            apply Finset.sum_congr rfl
            intro a _
            exact hterm a
      _ = 0 := Finset.sum_const_zero
  have hQ : wdot dv rho2 (fun _ => 1) = wdot dv rho1 (fun _ => 1) := by
    have hdiff_zero : 4 * pi_val * wdot dv (fun i => rho2 i - rho1 i) (fun _ => 1) = 0 := by
      calc
        4 * pi_val * wdot dv (fun i => rho2 i - rho1 i) (fun _ => 1)
          = wdot dv (fun i => 4 * pi_val * (rho2 i - rho1 i)) (fun _ => 1) :=
              (wdot_smul_left dv (4 * pi_val) (fun i => rho2 i - rho1 i) (fun _ => 1)).symm
        _ = wdot dv (vacuumOp D phi_tilde) (fun _ => 1) := by
              have hfun : (fun i => 4 * pi_val * (rho2 i - rho1 i)) = vacuumOp D phi_tilde :=
                funext (fun i => (h_vac_tilde i).symm)
              rw [hfun]
        _ = wdot dv phi_tilde (vacuumOp D (fun _ => 1)) :=
              (vacuumOp_symmetric dv D hskew phi_tilde (fun _ => 1)).symm
        _ = wdot dv phi_tilde (fun _ => 0) := by
              have hfun : vacuumOp D (fun _ => (1 : ℝ)) = fun _ => (0 : ℝ) :=
                funext h_vac_one
              rw [hfun]
        _ = 0 := by
              unfold wdot
              have h0 : ∑ i : Fin n, phi_tilde i * (0 : ℝ) = 0 := by
                simp only [mul_zero, Finset.sum_const_zero]
              rw [h0, mul_zero]
    have hsub0 : wdot dv (fun i => rho2 i - rho1 i) (fun _ => 1) = 0 :=
      (mul_eq_zero.mp hdiff_zero).resolve_left h4pi
    rw [wdot_sub_left] at hsub0
    linarith
  have h_pol2 : ∀ i,
      gpOp D eps k2 (fun j => phi_pol1 j - (alpha_sm2 - alpha_sm1)) i =
        4 * pi_val * rho2 i - gpOp D eps k2 phi_vac2 i := by
    intro i
    have hfun_eq :
        (fun j => phi_vac2 j + (phi_pol1 j - (alpha_sm2 - alpha_sm1))) =
        (fun j => (phi_vac1 j + phi_pol1 j) + phi_tilde j) := by
      funext j
      dsimp only [phi_tilde]
      ring
    have hsum :
        gpOp D eps k2 phi_vac2 i +
          gpOp D eps k2 (fun j => phi_pol1 j - (alpha_sm2 - alpha_sm1)) i =
          4 * pi_val * rho2 i := by
      calc
        gpOp D eps k2 phi_vac2 i +
            gpOp D eps k2 (fun j => phi_pol1 j - (alpha_sm2 - alpha_sm1)) i
          = gpOp D eps k2
              (fun j => phi_vac2 j + (phi_pol1 j - (alpha_sm2 - alpha_sm1))) i :=
              (gpOp_add D hD_add eps k2 phi_vac2
                (fun j => phi_pol1 j - (alpha_sm2 - alpha_sm1)) i).symm
        _ = gpOp D eps k2 (fun j => (phi_vac1 j + phi_pol1 j) + phi_tilde j) i := by
              rw [hfun_eq]
        _ = gpOp D eps k2 (fun j => phi_vac1 j + phi_pol1 j) i +
              gpOp D eps k2 phi_tilde i :=
              gpOp_add D hD_add eps k2 (fun j => phi_vac1 j + phi_pol1 j) phi_tilde i
        _ = (gpOp D eps k2 phi_vac1 i + gpOp D eps k2 phi_pol1 i) +
              gpOp D eps k2 phi_tilde i := by
              rw [gpOp_add D hD_add eps k2 phi_vac1 phi_pol1 i]
        _ = (gpOp D eps k2 phi_vac1 i +
              (4 * pi_val * rho1 i - gpOp D eps k2 phi_vac1 i)) +
              4 * pi_val * (rho2 i - rho1 i) := by
              rw [h_pol1 i, h_gp_tilde i]
        _ = 4 * pi_val * rho2 i := by ring
    linarith
  have h_vac_recip : wdot dv phi_vac2 rho1 = wdot dv rho2 phi_vac1 := by
    have h1 : wdot dv phi_vac2 (fun i => rho1 i - rho_bar) = wdot dv phi_vac2 rho1 := by
      rw [wdot_sub_right]
      have hbar : wdot dv phi_vac2 (fun _ => rho_bar) = 0 := by
        have heq : (fun _ : Fin n => rho_bar) = fun _ => rho_bar * (1 : ℝ) := by
          funext _
          ring
        rw [heq, wdot_smul_right, h_mean2, mul_zero]
      rw [hbar, sub_zero]
    have h2 : wdot dv (fun i => rho2 i - rho_bar) phi_vac1 = wdot dv rho2 phi_vac1 := by
      rw [wdot_sub_left]
      have hbar : wdot dv (fun _ => rho_bar) phi_vac1 = 0 := by
        have heq : (fun _ : Fin n => rho_bar) = fun _ => rho_bar * (1 : ℝ) := by
          funext _
          ring
        rw [heq, wdot_smul_left, wdot_comm, h_mean1, mul_zero]
      rw [hbar, sub_zero]
    have h4 : 4 * pi_val * wdot dv phi_vac2 rho1 = 4 * pi_val * wdot dv rho2 phi_vac1 := by
      calc
        4 * pi_val * wdot dv phi_vac2 rho1
          = 4 * pi_val * wdot dv phi_vac2 (fun i => rho1 i - rho_bar) := by rw [h1]
        _ = wdot dv phi_vac2 (fun i => 4 * pi_val * (rho1 i - rho_bar)) :=
              (wdot_smul_right dv (4 * pi_val) phi_vac2 (fun i => rho1 i - rho_bar)).symm
        _ = wdot dv phi_vac2 (vacuumOp D phi_vac1) := by
              have hfun : (fun i => 4 * pi_val * (rho1 i - rho_bar)) = vacuumOp D phi_vac1 :=
                funext (fun i => (h_vac1 i).symm)
              rw [hfun]
        _ = wdot dv (vacuumOp D phi_vac2) phi_vac1 :=
              vacuumOp_symmetric dv D hskew phi_vac2 phi_vac1
        _ = wdot dv (fun i => 4 * pi_val * (rho2 i - rho_bar)) phi_vac1 := by
              have hfun : vacuumOp D phi_vac2 = fun i => 4 * pi_val * (rho2 i - rho_bar) :=
                funext h_vac2
              rw [hfun]
        _ = 4 * pi_val * wdot dv (fun i => rho2 i - rho_bar) phi_vac1 :=
              wdot_smul_left dv (4 * pi_val) (fun i => rho2 i - rho_bar) phi_vac1
        _ = 4 * pi_val * wdot dv rho2 phi_vac1 := by rw [h2]
    have hsub : 4 * pi_val * (wdot dv phi_vac2 rho1 - wdot dv rho2 phi_vac1) = 0 := by
      linarith
    have h0 : wdot dv phi_vac2 rho1 - wdot dv rho2 phi_vac1 = 0 :=
      (mul_eq_zero.mp hsub).resolve_left h4pi
    linarith
  have h_pol_recip :
      wdot dv (fun i => rho2 i - rho1 i) phi_pol1 =
        -(wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1)) := by
    have h_rhs :
        wdot dv phi_tilde (gpOp D eps k2 phi_pol1) =
          4 * pi_val * wdot dv phi_tilde rho1 -
            wdot dv phi_tilde (gpOp D eps k2 phi_vac1) := by
      have hfun : gpOp D eps k2 phi_pol1 =
          fun i => 4 * pi_val * rho1 i - gpOp D eps k2 phi_vac1 i :=
        funext h_pol1
      rw [hfun, wdot_sub_right, wdot_smul_right]
    have h_vac_term :
        wdot dv phi_tilde (gpOp D eps k2 phi_vac1) =
          4 * pi_val * wdot dv (fun i => rho2 i - rho1 i) phi_vac1 := by
      rw [gpOp_symmetric dv D hskew eps k2 phi_tilde phi_vac1]
      have hfun : gpOp D eps k2 phi_tilde =
          fun i => 4 * pi_val * (rho2 i - rho1 i) :=
        funext h_gp_tilde
      rw [hfun, wdot_smul_left]
    have h_tilde_rho1 :
        wdot dv phi_tilde rho1 =
          wdot dv phi_vac2 rho1 - wdot dv phi_vac1 rho1 -
            wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1) := by
      change wdot dv (fun i => (phi_vac2 i - phi_vac1 i) - (alpha_sm2 - alpha_sm1)) rho1 = _
      rw [wdot_sub_left, wdot_sub_left]
      have hconst :
          wdot dv (fun _ => alpha_sm2 - alpha_sm1) rho1 =
            wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1) := by
        have heq : (fun _ : Fin n => alpha_sm2 - alpha_sm1) =
            fun _ => (alpha_sm2 - alpha_sm1) * (1 : ℝ) := by
          funext _
          ring
        rw [heq, wdot_smul_left, wdot_comm]
        ring
      rw [hconst]
    have h_drho_vac1 :
        wdot dv (fun i => rho2 i - rho1 i) phi_vac1 =
          wdot dv rho2 phi_vac1 - wdot dv phi_vac1 rho1 := by
      rw [wdot_sub_left, wdot_comm dv rho1 phi_vac1]
    have h4 :
        4 * pi_val * wdot dv (fun i => rho2 i - rho1 i) phi_pol1 =
          4 * pi_val * (-(wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1))) := by
      calc
        4 * pi_val * wdot dv (fun i => rho2 i - rho1 i) phi_pol1
          = wdot dv (fun i => 4 * pi_val * (rho2 i - rho1 i)) phi_pol1 :=
              (wdot_smul_left dv (4 * pi_val) (fun i => rho2 i - rho1 i) phi_pol1).symm
        _ = wdot dv (gpOp D eps k2 phi_tilde) phi_pol1 := by
              have hfun : (fun i => 4 * pi_val * (rho2 i - rho1 i)) = gpOp D eps k2 phi_tilde :=
                funext (fun i => (h_gp_tilde i).symm)
              rw [hfun]
        _ = wdot dv phi_tilde (gpOp D eps k2 phi_pol1) :=
              (gpOp_symmetric dv D hskew eps k2 phi_tilde phi_pol1).symm
        _ = 4 * pi_val * wdot dv phi_tilde rho1 -
              wdot dv phi_tilde (gpOp D eps k2 phi_vac1) := h_rhs
        _ = 4 * pi_val *
              (wdot dv phi_vac2 rho1 - wdot dv phi_vac1 rho1 -
                wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1)) -
              4 * pi_val * (wdot dv rho2 phi_vac1 - wdot dv phi_vac1 rho1) := by
              rw [h_tilde_rho1, h_vac_term, h_drho_vac1]
        _ = 4 * pi_val * (-(wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1))) := by
              rw [h_vac_recip]
              ring
    have hsub :
        4 * pi_val *
          (wdot dv (fun i => rho2 i - rho1 i) phi_pol1 -
            (-(wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1)))) = 0 := by
      linarith
    have h0 :
        wdot dv (fun i => rho2 i - rho1 i) phi_pol1 -
          (-(wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1))) = 0 :=
      (mul_eq_zero.mp hsub).resolve_left h4pi
    linarith
  refine ⟨hQ, h_pol2, ?_⟩
  have h_e2 :
      wdot dv rho2 (fun i => phi_pol1 i - (alpha_sm2 - alpha_sm1)) =
        wdot dv rho2 phi_pol1 - wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1) := by
    rw [wdot_sub_right]
    have hconst :
        wdot dv rho2 (fun _ => alpha_sm2 - alpha_sm1) =
          wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1) := by
      have heq : (fun _ : Fin n => alpha_sm2 - alpha_sm1) =
          fun _ => (alpha_sm2 - alpha_sm1) * (1 : ℝ) := by
        funext _
        ring
      rw [heq, wdot_smul_right, hQ]
      ring
    rw [hconst]
  have h_rho2_pol1 :
      wdot dv rho2 phi_pol1 =
        wdot dv rho1 phi_pol1 - wdot dv rho1 (fun _ => 1) * (alpha_sm2 - alpha_sm1) := by
    have hsub := h_pol_recip
    rw [wdot_sub_left] at hsub
    linarith
  rw [h_e2, h_rho2_pol1]
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
  have hlog : HasDerivAt (fun x => Real.log (x / n_c)) ((1 / n_c) / (n / n_c)) n :=
    hdiv.log hdiv_ne
  have hu : HasDerivAt (fun x => Real.log (x / n_c) / (sigma * Real.sqrt 2))
      (((1 / n_c) / (n / n_c)) / (sigma * Real.sqrt 2)) n :=
    hlog.div_const (sigma * Real.sqrt 2)
  have hcomp := herfc.comp n hu
  have hS : HasDerivAt (vaspsolSwitch erfc n_c sigma)
      ((1 / 2 : ℝ) *
        ((-(2 / Real.sqrt Real.pi) *
            Real.exp (-((Real.log (n / n_c) / (sigma * Real.sqrt 2)) ^ 2))) *
          (((1 / n_c) / (n / n_c)) / (sigma * Real.sqrt 2)))) n :=
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

#print axioms wdot_comm
#print axioms wdot_add_right
#print axioms wdot_add_left
#print axioms wdot_neg_right
#print axioms wdot_sub_right
#print axioms wdot_sub_left
#print axioms wdot_smul_right
#print axioms wdot_smul_left
#print axioms wdot_mul_diag
#print axioms wdot_sum_right
#print axioms wdot_diag_nonneg
#print axioms wdot_diag_eq_zero_iff
#print axioms wdot_diag_pos_of_exists
#print axioms wdot_diag_pos
#print axioms nyquist_zeroed_multiplier_skew_adjoint
#print axioms nyquist_nonzero_breaks_reality
#print axioms gpOp_add
#print axioms gpOp_quad_form
#print axioms gpOp_symmetric
#print axioms gpOp_psd
#print axioms gpOp_quad_eq_zero_iff
#print axioms gpOp_spd_on_grad_complement
#print axioms gpOp_spd_lpb_cavity_vanishing
#print axioms gpOp_spd_const_mode_lpb
#print axioms gpOp_spd_uniform_screening
#print axioms dualEnergy_on_shell
#print axioms dualEnergy_perturb_exact
#print axioms dualEnergy_global_max
#print axioms gpOp_vacuum_eq
#print axioms vacuumOp_symmetric
#print axioms gpOp_apply_zero
#print axioms polarization_rhs_vanishes_in_vacuum
#print axioms zero_double_counting_cancellation
#print axioms bg_offset_algebraic_cancellation
#print axioms full_electrostatic_width_invariance
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
