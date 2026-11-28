# Track B Gate 3: Architecture Proposal — Typed Electrostatics & Implicit Solvation Seam

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Track B — Phase B1)
- **Phase**: Gate 3 — Architecture Proposal (with Phase B2 GC-DFT & Track C Electrostatic Embedding Readiness)
- **Date**: 2026-11-28
- **Prerequisites**: [Gate 0 Assessment](repository-assessment.md), [Gate 1 Theory](implicit-solvation-theory.md), [Gate 1 Algorithms](algorithm-comparison.md), [Gate 2 Feasibility Report](feasibility-report.md), [Decision Log](decision-log.md) (`D-001`–`D-012`)
- **Proposed ADR**: [ADR 0061: Typed Electrostatics and Implicit Solvation Seam](../../decisions/0061-typed-electrostatics-and-solvation-seam.md) (Draft — pending Gate 5 acceptance)

---

## 1. Architectural Objectives & Invariants

1. **Zero Pure-PW Regression (Bit-Identical Off-Path)**: When `ScfOptions::solvation` is `None`, zero solvation allocations, zero extra FFT plans, zero extra MPI collectives, zero branches inside `DistributedRealHamiltonian::apply_block`, `DistributedComplexHamiltonian::apply_block`, `ref_davidson`, or `ref_davidson_complex`, and zero new JSON keys are emitted.
2. **ADR 0047 Single-Path Compliance**: Integrate strictly into the production `refscf/` driver family (`driver`, `k_driver`, `electronic`) and `refstart.rs`. Frozen `oracle/` (`strict-oracle` feature) has zero behavioral changes (only two mechanical `solvation: None` struct-literal completions in `oracle/collinear.rs` and `oracle/scf.rs`), and `cargo check/clippy -p dft-engine --no-default-features` stays green.
3. **Zero Crate-Graph Churn**: No changes to `ownership.toml` or workspace `Cargo.toml` dependencies (`dft-units` → `dft-physics` → `dft-engine` → `dft-cli` edges already exist).
4. **Compile-Time Electrostatic Sign Safety (ADR 0017 / `D-003`)**: Eliminate $\phi \leftrightarrow U_{\mathrm{el}} = -\phi$ sign confusion across Track B (solvation, GC-DFT) and Track C (QM/MM, polarizable MLIP) via explicit scalar newtypes in `dft-units`.
5. **Fail-Closed Derivative Coverage**: Phase B1 provides analytic Gaussian ionic solvation forces (plus NLCC core-cavity forces when core density enters the cavity) and explicitly marks solvation stress `stress_available = false` so stress requests fail closed until Phase B1b/B2 implements analytic solvation stress.
6. **Future Track B2 / Track C Compatibility Without Premature Abstraction**: Keep the Phase B1 boundary minimal and concrete (`SolvationState` / `SolvationFieldSummary` / `SolvationEnergyLedger`), while keeping the $\phi_{\mathrm{bulk}} = 0$ gauge, Bregman grand-potential hooks (`grand_canonical_bregman_exact`), and `ElectrostaticPotentialHaPerE` types directly reusable by Phase B2 GC-DFT and Track C electrostatic embedding.

---

## 2. Module-by-Module Design & Exact Seams

### 2.1 Typed Electrostatic Semantics (`crates/dft-units/src/lib.rs`)

Add explicit scalar newtypes immediately after existing `Bohr`/`Hartree`/`Rydberg`/`ElectronVolt` definitions (`crates/dft-units/src/lib.rs:53`):

```rust
/// Electrostatic potential phi in Hartree per +1e (positive near nuclei).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct ElectrostaticPotentialHaPerE(pub f64);

/// Potential energy U_el = -phi of an electron (charge -1e) in Hartree.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct ElectronPotentialEnergyHa(pub f64);

impl ElectrostaticPotentialHaPerE {
    pub const fn to_electron_energy(self) -> ElectronPotentialEnergyHa {
        ElectronPotentialEnergyHa(-self.0)
    }
}

impl ElectronPotentialEnergyHa {
    pub const fn to_electrostatic_potential(self) -> ElectrostaticPotentialHaPerE {
        ElectrostaticPotentialHaPerE(-self.0)
    }
}
```

- **No implicit `From<f64>`** and no direct addition between `ElectrostaticPotentialHaPerE` and `ElectronPotentialEnergyHa`: every boundary between `dft-physics::solvation` ($\phi_{\mathrm{pol}}, \phi_{\mathrm{solv}}$) and `dft-engine` ($V_{\mathrm{eff}} = U_{\mathrm{el}}$) converts explicitly via `.to_electron_energy()`.

### 2.2 Physics Kernels (`crates/dft-physics/src/solvation/`)

Carry the Gate 2 verified module (`mod.rs`, `cavity.rs`, `poisson.rs`) forward with three targeted Gate 3–4 upgrades:

1. **Zero-copy `n_g` reuse overload (`evaluate_with_density_g`)**:
   - Accept caller-supplied `density_g: &[Complex64]` computed by `refscf::map::fields()` (`map.rs:75`) right before Hartree overwrites `coeff`, eliminating 1 full 3D forward FFT per SCF map.
2. **Dual-solver production policy (`D-012`)**:
   - Default solver: `SolverKind::DirectPcg` (strictly SPD, exact `< 1e-9` discrete variational consistency on all grids, `1.5e-16` Ha MPI parity, 0 setup FFTs).
   - Opt-in fast solver: `SolverKind::FisicaroPcg` (3.4×–12.7× faster, 2 FFTs/iter), with:
     - Grid-resolution guard emitting a diagnostic preflight warning if $\max_a \Delta x_a > 0.25\text{ Bohr}$ (where Gate 2 proved coarse-grid operator aliasing can exceed `1e-4 Ha`; consistent with `SolvationWorkspace::max_grid_spacing_bohr()` and ADR 0061),
     - Automatic fallback to `DirectPcg` on `SolvationError::Breakdown(_)`,
     - Gate 4 R&D experiment (`solvation_suite_c.rs`): evaluated consistent $A_{\mathrm{Fis}}\phi_{\mathrm{vac}}$ polarization RHS and inexact Fisicaro preconditioning inside Flexible PCG (documented as a Gate 4 negative result for replacing `DirectPcg`/`FisicaroPcg`; see ADR 0061 and `D-014`).
   - Remove `HybridFisicaroDirect` from the production CLI enum (retained only in Gate 2 benchmarks/docs as a documented negative result).
3. **Ionic pseudo-charge construction & cavitation regularization**:
   - Construct `Vec<IonicPseudoCharge>` directly from `(atoms, species_z_val, width_bohr)` in `refstart.rs` and `forces.rs` with default `width_bohr = 0.60 Bohr` (verified width-invariant in Lean 4 via `bg_offset_algebraic_cancellation` and PDE-level `full_electrostatic_width_invariance`, plus Rust FD tests).
   - Default `cavity_grad_reg = 1.0e-4 Bohr^-1` in `SolvationParams` (`vacuum()`, `vaspsol_water()`, `sccs_water()`): avoids amplifying IEEE-754 FP64 FFT rounding noise (`~1e-15`) in $|\mathbf{D}S|$ into `O(1e-3 Ha)` asymmetric noise in $V_{\mathrm{cav}}$, while perturbing $E_{\mathrm{cav}}$ by `< 2e-8 Ha` and preserving exact discrete variational consistency.

### 2.3 Engine Options, State & Identity (`crates/dft-engine/src/`)

| File & Seam | Addition | Off-Path Behavior (`None`) |
| --- | --- | --- |
| `model/options.rs:150` (`ScfOptions`) | `pub solvation: Option<Box<SolvationConfig>>` (boxed to keep `ScfOptions` compact; wraps `SolvationParams`, `SolverConfig`, `ionic_width_bohr`, `adaptive_tolerance`) | `None` in `default()` and `strict_reference()` |
| `force_identity.rs:52` (`ModelDescription`) | `solvation: Option<SolvationIdentity>` (bit-encoded params, cavity model, solver kind, ionic width) | `None`; rejects mismatched solvation models on derivative calls |
| `refstart.rs:128` (`StartState`) | `pub solvation: Option<Box<SolvationState>>` (reused by `KStartState.common` via `DerefMut` at L1811) | `None` (0 bytes heap, 0 FFT plans) |
| `refscf/profile.rs:55` | `pub fn next_solvation_tol(previous, dr2_ry, electrons, target) -> f64` | Unused when `solvation.is_none()` |

### 2.4 Density-to-Potential Maps & `EnergyAudit` (`refscf/{map,k_map,electronic_map}.rs`, `energy.rs`)

1. **`refscf/map.rs`**:
   - Between L101 (end of PBE borrow of `&coeff`) and L102 (Hartree overwrite of `coeff`), if `let Some(solv) = &mut state.solvation`, evaluate `SolvationWorkspace::evaluate_with_density_g(...)` using the current ratcheted inner tolerance.
   - Add `v_total_planes` into `hxc` at L170–172 **before** `potential_symmetry` (L175, ADR 0049), ensuring full crystal-symmetry invariance of $V_H + V_{\mathrm{xc}} + V_{\mathrm{solv}}$.
   - Store `SolvationFieldSummary { ledger, v_solv, v_solv_integral }` in `Fields::solvation: Option<...>`.
   - `cached_fields()` (L388–410) caches `Fields` by bit-identical density, automatically avoiding redundant PDE solves on repeated evaluations.
2. **`refscf/k_map.rs`**:
   - Delegates directly to `map::cached_fields(&mut state.common, ...)` (L60–69); zero duplicate solver code.
3. **`refscf/electronic_map.rs`**:
   - Unpolarized (`density.len() == 1`, L386): inherits `scalar.hxc` (already containing `V_solv`) and forwards `scalar.solvation`.
   - **Collinear spin-polarized (`density.len() == 2`, L446–522)**: explicitly add spin-independent `V_solv[n_tot]` to **both** `potentials[0]` and `potentials[1]` before `potential_symmetry`, and subtract `v_solv` alongside `vh` when computing the `dvxc` response diagnostic at L1033–1034.
4. **`energy.rs` (`EnergyAudit`) & Harris–Foulkes / Bregman Accounting**:
   - Add `pub solvation: Option<SolvationEnergyLedger>` to `EnergyAudit` (`electrostatic_ha`, `background_correction_ha`, `cavitation_ha`, `total_ha`, `potential_integral_ha`, `iterations`, `relative_residual`, `fft_count`).
   - `total_hartree()` (`energy.rs:72`): only adds `s.total_ha` inside `if let Some(s) = &self.solvation` (strict pure-PW bit identity).
   - `band_form`: adds `+ E_solv[n_out]` (`double_counting` already includes $\langle n_{\mathrm{out}}, V_{\mathrm{solv}}[n_{\mathrm{in}}]\rangle$ via `hxc`).
   - `harris_foulkes`: adds $+E_{\mathrm{solv}}[n_{\mathrm{in}}] - \langle n_{\mathrm{in}}, V_{\mathrm{solv}}[n_{\mathrm{in}}]\rangle$, matching the Lean 4 theorem `harris_foulkes_bregman_exact` (`HF = E[n_out]` at convergence).
   - `functional_difference_bound`: separates $L_{\mathrm{solv}} = \langle n_{\mathrm{out}} - n_{\mathrm{in}}, V_{\mathrm{solv}}[n_{\mathrm{in}}]\rangle$ from the XC remainder and adds $|E_{\mathrm{solv}}[n_{\mathrm{out}}] - E_{\mathrm{solv}}[n_{\mathrm{in}}] - L_{\mathrm{solv}}|$.
5. **Adaptive Inner Tolerance Ratchet (`driver.rs`, `k_driver.rs`, `electronic.rs`)**:
   - Early SCF iterations use ratcheted $\tau_{\mathrm{inner}}^{(t)} = \operatorname{clamp}(\min(\tau_{\mathrm{prev}}, 0.05\,\max(\mathrm{dr2\_ry}, 0)), \tau_{\mathrm{target}}, \max(\tau_{\mathrm{target}}, 10^{-4}))$ (`refscf::profile::next_solvation_tol`, linear in `dr2_ry` mirroring `next_ethr_ry` so PDE noise never creates a Broyden stall floor) with warm-start enabled.
   - Final verification map (`driver.rs:736`, `k_driver.rs:552`, `electronic.rs:1305`) always solves to exact `config.tolerance` (default `1e-10`).

### 2.5 Forces & Fail-Closed Stress (`forces.rs`, `stress.rs`)

- **`forces.rs`**:
  - Keep base `derivative_ledger()` unchanged (8 base NC entries for pure-PW runs) and append `"solvation"` (`owner: "physics::solvation"`, `force_available: true`, `stress_available: false`, `included_in: "E"`) inside `validate()` only when `audit.solvation.is_some()` (preserving bit-identical off-path `ForceResult::ledger` and `results.json`).
  - Add `pub solvation: Option<Vec<[f64; 3]>>` to `ForceResult`.
  - In `forces()`, when `state.options.solvation.is_some()`, evaluate exact ionic solvation forces on the converged density at tight tolerance (`<= 1e-11`), symmetrize alongside `[ewald, local, nonlocal, core]` (L476), and accumulate into `ForceResult::values`.
  - In `validate()`, verify `require_complete(&ledger)` after appending active `"solvation"`.
- **`stress.rs` (Phase B1 Fail-Closed)**:
  - Reject in `preconditions()` (`stress.rs:114`) when `state.options.solvation.is_some()` with `"observables.stress: incomplete derivative coverage: solvation"`.
  - Guard `frozen_pass()` (`stress.rs:454`) so `entry.stress_available = entry.term_id != "solvation"` and enforce `require_stress_complete(&ledger)`.
  - Because `refscf/api.rs:326–340` (`forces_maybe_joint`) falls back gracefully to forces-only when `forces_and_frozen_stress` declines, `forces = true` succeeds while `stress = true` fails closed with an explicit diagnostic.

### 2.6 CLI TOML Input & `results.json` Schema (`crates/dft-cli/src/`)

Additive `[solvation]` TOML table in `crates/dft-cli/src/config.rs` (`#[serde(deny_unknown_fields)]`):

```toml
[solvation]
model = "vaspsol"           # "vaspsol" | "sccs"
epsilon_bulk = 78.4         # >= 1.0 (1.0 with kappa=0, tau=0, beta=0 => exact vacuum)
nc_bohr3 = 0.0003704618     # VASPsol default (0.0025 Angstrom^-3 via dft_units::BOHR_IN_ANGSTROM)
sigma = 0.6                 # VASPsol default
# rho_min_bohr3 = 1e-4      # SCCS optional
# rho_max_bohr3 = 1.5e-3    # SCCS optional
# log_epsilon = true        # SCCS optional
kappa_inv_bohr = 0.0        # Debye screening wavevector (Bohr^-1), 0.0 = neutral solvent
tau_ha_bohr2 = 0.0          # Cavitation surface tension (Ha/Bohr^2)
beta_ha_bohr3 = 0.0         # Cavitation volume pressure (Ha/Bohr^3)
solver = "direct_pcg"       # "direct_pcg" (default) | "fisicaro_pcg"
tolerance = 1e-10           # Final PDE relative residual target
max_iterations = 300
warm_start = true
adaptive_tolerance = true
ionic_width_bohr = 0.60     # Gaussian ionic pseudo-charge width (width-invariant)
```

- **Fail-closed legacy paths**: Reject `[solvation]` cleanly in `oracle_options` (`run.rs:348`) and legacy `spin_run.rs`/`finite_run.rs` with an explicit message directing users to the production reference path.
- **Additive `results.json`**: Emit `"solvation"` inside `"energy"` (`energy_json`, L1236), top-level `"solvation"` telemetry block (solver, iterations, FFTs, residual, cavity model), and `"per_term"."solvation"` in `force_json` (L3604) **only** when `[solvation]` is active. When omitted, `results.json` is byte-for-byte unchanged.

---

## 3. Phase B2 (GC-DFT) & Track C (Electrostatic Embedding) Readiness

- **Phase B2 (Grand-Canonical DFT)**:
  - $\kappa > 0$ already pins $\phi_{\mathrm{bulk}} = 0$ uniquely, providing the exact absolute electrostatic reference needed to map $\mu_e$ to electrode potential $U_{\mathrm{SHE}}$.
  - `SolvationState` is stateless with respect to $N_e$ (charge neutrality is NOT assumed; `E_bg,corr` algebraic offset cancellation and full PDE-level electrostatic width invariance for arbitrary $Q_{\mathrm{ion}} - N_e$ are formally proven in Lean 4 via `bg_offset_algebraic_cancellation` and `full_electrostatic_width_invariance`).
  - `grand_canonical_bregman_exact` is already formalized in `BaandiMath.Solvation` and plugs straight into `refscf/electronic.rs` when adding a fixed-$\mu$ outer root-finder.
- **Track C (Electrostatic QM/MM & Polarizable MLIP)**:
  - Reuses `ElectrostaticPotentialHaPerE` / `ElectronPotentialEnergyHa` in `dft-units`, the exact `hxc` injection seam in `map.rs:170` + `electronic_map.rs:446–522`, the `EnergyAudit` additive ledger pattern, and the `derivative_ledger()` force/stress coverage contract.

---

## 4. Gate 4 Implementation Scope & Rollout Plan (Research Branch Only)

Gate 4 executes on an isolated research branch (`research/track-b-gate4-scf`) without modifying `main` until Gate 5 approval:

1. **Step 4.1**: Add `ElectrostaticPotentialHaPerE` / `ElectronPotentialEnergyHa` to `dft-units` + `evaluate_with_density_g` to `dft-physics::solvation`.
2. **Step 4.2**: Wire `ScfOptions::solvation`, `StartState::solvation`, `map.rs` / `k_map.rs` / `electronic_map.rs`, `EnergyAudit`, and adaptive inner tolerance ratchet across all three reference drivers.
3. **Step 4.3**: Wire `forces.rs`, fail-closed `stress.rs`, `force_identity.rs`, and `[solvation]` CLI TOML + `results.json` output.
4. **Step 4.4**: Execute full [Gate 3 Verification Plan](verification-plan.md) (L1–L4) and [Gate 3 Benchmark Plan](benchmark-plan.md) (L5) on `dev08`, including real-molecule/slab SCF convergence, finite-difference SCF total-energy force validation, zero-regression off-path audit, and Gate 4 Flexible-PCG / consistent-Fisicaro speed exploration.
5. **Step 4.5**: Deliver Gate 5 Production Review Package to Research Lead for merge sign-off.
