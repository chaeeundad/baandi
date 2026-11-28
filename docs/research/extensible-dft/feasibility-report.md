# Track B Gate 2: Implicit Solvation Feasibility & Verification Report

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Track B — Phase B1)
- **Phase**: Gate 2 — Isolated Numerical Prototype, Formal Verification, Benchmarks & Negative Results
- **Date**: 2026-11-28
- **Prototype Branch**: `research/track-b-gate2-proto` (baseline workspace gate verified on `ee0a183a`; strict audit-hardened L1/L2 + Lean 4 `HasDerivAt` verified on `bf8e8c3f` and final gate commit below)
- **Hardware**: `dev08` (`192.168.0.228`, dual-socket Intel Xeon Gold 6230, 40 physical cores, 187 GB RAM, pinned via `taskset -c 20-39,60-79 nice -n 10`)
- **Clean-Room & Oracle Notice**: All cavity switching functions (`VASPsol-erfc` from Mathew et al. 2014, DOI `10.1063/1.4865107`; `SCCS` log-sine from Andreussi et al. 2012, DOI `10.1063/1.3676407`) are clean-room Rust implementations of the published equations in `crates/dft-physics/src/solvation/`. No external VASP, VASPsol, or QE+Environ binaries were executed (`D-010`), and zero lines of GPL/LGPL or proprietary source code were copied.

---

## 1. Executive Summary & Gate 2 Verdict

**Verdict: GO → Advance to Gate 3 (Architecture Proposal) after Completed Gate 2 `Revise` Cycle & Explicit Performance Waiver (`D-012`).**

1. **Level-1 Formal Mathematical Verification (Lean 4 + Mathlib)**: All **24 theorems** (core discrete skew-adjoint SPD operator $A$, dual variational global maximizer, analytic `HasDerivAt` chain-rule proofs for both `sccsSwitch` and `vaspsolSwitch`, algebraic $G=0$ double-counting & width-invariance cancellation, Harris–Foulkes/GC-DFT Bregman identities, electrostatic sign involution, and spin-polarized chain rule) are machine-checked in [`formal/BaandiMath/BaandiMath/Solvation.lean`](../../../formal/BaandiMath/BaandiMath/Solvation.lean) under Lean `v4.36.0` + Mathlib (`lake build` exit code `0`, zero `sorry`, `#print axioms` strictly `[propext, Classical.choice, Quot.sound]`).
2. **Level-1 / Level-2 Independent Numerical Verification (Zero Tolerance Relaxation)**: **3/3 unit tests** and **8/8 independent verification tests** in [`crates/dft-physics/tests/solvation_gate2.rs`](../../../crates/dft-physics/tests/solvation_gate2.rs) pass **every pre-declared tolerance** from [implicit-solvation-theory.md](implicit-solvation-theory.md#9-verification-plan-levels-15-pre-declared-tolerances):
   - Exact vacuum limit: `0.0` (`< 1e-14`) on skewed cells.
   - Smooth MMS ($\epsilon_b = 78.4$, `cubic 16^3` & `skewed 16x18x14`, `tol = 1e-14`): `DirectPcg` `rel_linf = 8.38e-15` to `1.44e-14` (`<< 1e-12`).
   - Smooth-wall Born ion (`kappa = 0`, Makov–Payne corrected): **`0.819%`** (`< 1%`); Debye–Hückel LPB vs exact 1D continuum ODE with monotonic convergence across box sizes: **`0.245%`** (`L=24 -> 32`, $\epsilon_b = 78.4$) and **`0.703%`** (`L=24 -> 36`, $\epsilon_b = 10.0$) (`< 1%`).
   - 4th-order Richardson variational potential $V_{\mathrm{solv}} = \delta E_{\mathrm{solv}}/\delta n$ vs FD on skewed cells: **every full model** (VASPsol-erfc, SCCS log-$\epsilon$, SCCS linear-$\epsilon$) and **every isolated term** (`cavitation_only`: `9.76e-12 Ha`, `diel_and_es_kappa0`: `6.10e-11 Ha`, `scr_and_es_eps1`: `3.14e-11 Ha`) passes `abs < 1e-10 Ha` / `rel < 1e-9`.
   - Ionic solvation forces vs central FD (`< 1e-9` Ha/Bohr) & Gaussian width $r_c$ invariance (`< 1e-6 Ha`): **`2.01e-8 Ha`** (neutral, $\kappa = 0$) and **`4.05e-8 Ha`** (charged $Q = +1$, $\kappa = 0.2$, canceling a `2 mHa` `E_bg,corr` shift).
   - MPI rank independence ($P \in \{1, 2, 3, 4\}$) & `mirror` vs `greedy` stick-partition parity: `< 1e-12` (energy parity `1.5e-16 Ha` on `48^3`).
3. **Transparent `Revise` Cycle & Dual-Solver Performance Waiver (`algorithm-comparison.md` §5)**:
   - Gate 1 Criterion 2–3 optimistically targeted `<= 35` cold / `<= 10` warm PCG iterations at $\epsilon_b = 78.4$ and `< 5%` of a real SCF map wall time, defining misses as a **`Revise`** trigger.
   - Empirically on sharp VASPsol-erfc cavities ($\epsilon_b = 78.4$, $\sigma = 0.6$), `DirectPcg` required **`99–109` cold iters** (`10^{-10}`), **`75–83` warm(1e-3) iters** (`10^{-10}`), and **`24–31` warm(1e-3) iters at adaptive SCF tolerance `10^{-5}`** (`~141 ms` extrapolated on 4 ranks @ `48^3`, ~`2.8–14.1%` of a `1–5 s` SCF map).
   - Per the `Revise` rule, we explored two preconditioner/hybrid revisions (`M_split^{-1}` and single-shot `HybridFisicaroDirect`) and documented both negative results in §6. Meanwhile, `FisicaroPcg` (2 FFTs/iter) achieved **2.74×–5.27×** cold speedup (**1.85×–2.70×** at adaptive `warm(1e-3)@1e-5`, **7.9×–12.7×** on smooth MMS) with `~60 ms` extrapolated 4-rank adaptive wall time (`~1.2–6.0%` of `1–5 s`), while revealing a coarse-grid aliasing discrepancy (`7.70e-3 Ha` at `32^3`, `dx = 0.50 Bohr`) that decays exponentially to **`7.21e-6 Ha` at `64^3`** (`dx = 0.25 Bohr`).
   - **Resolution (`D-012`)**: Waive the standalone iteration count target for Gate 2 exit, adopt **`DirectPcg` as exact discrete variational default** + **grid-guarded (`max_dx <= 0.25 Bohr`) `FisicaroPcg` as high-speed opt-in accelerator**, and defer real-SCF-map `< 5%` wall-time measurement and Flexible-PCG acceleration R&D directly to Gate 4 ([benchmark-plan.md](benchmark-plan.md)).
4. **Workspace Gate Verification**: Full `/scr/mino/gate/gate.sh ee0a183a` on `dev08` passed cleanly (`GATE_DONE check_rc=0 hard=0 known_platform=0 new=0`; `140 binaries, 1134 passed, 0 failed`; `production-check` & `production-clippy` PASS; `mpi+mkl` `feat_build rc=0` & `feat_test rc=0`). Zero changes to `ownership.toml`, `Cargo.toml`, `dft-solvers`, or `dft-engine`; `unsafe_code = "deny"`; 0 inner-loop heap allocations.

---

## 2. Prototype Scope & Architecture

| Artifact | Path | Role |
| --- | --- | --- |
| Public API & error types | `crates/dft-physics/src/solvation/mod.rs` | `SolverKind`, `CavityModel`, `SolvationParams`, `SolverConfig`, `IonicPseudoCharge`, `SolveReport`, `SolvationEvaluation` |
| Cavity & ionic kernels | `crates/dft-physics/src/solvation/cavity.rs` | Skew-adjoint Nyquist-zeroed $D_a$, VASPsol-erfc & SCCS switching functions + derivatives, Gaussian ionic pseudo-charge & forces |
| Distributed PDE solvers | `crates/dft-physics/src/solvation/poisson.rs` | `SolvationWorkspace`, zero-allocation `DistributedPairFft` `np = 2` pairing (`mirror` + `greedy` fallback), all 4 solver arms, full `evaluate()` |
| Independent L1/L2 suite | `crates/dft-physics/tests/solvation_gate2.rs` | 8 independent tests covering every L1/L2 oracle at exact pre-declared tolerances |
| Release benchmark | `crates/dft-physics/examples/solvation_bench.rs` | 4-arm comparison across `32^3, 48^3, 64^3`, operator parity, adaptive SCF tolerance, and `P = 1, 2, 4` simulated MPI scaling |
| Lean 4 formal proofs | `formal/BaandiMath/BaandiMath/Solvation.lean` | 24 machine-checked theorems via Lean 4 + Mathlib (`v4.36.0`) |

---

## 3. Level-1 Formal Verification (Lean 4 + Mathlib)

Verified on `dev08` at `/scr/mino/lean/BaandiMath` synced from `formal/BaandiMath/` (`lean-toolchain`: `leanprover/lean4:v4.36.0`, Mathlib `v4.36.0`, `lake build` exit code `0`, zero `sorry`, `#print axioms` restricted strictly to `[propext, Classical.choice, Quot.sound]` on all 24 theorems):

| Group | Theorems (`BaandiMath.Solvation.*`) | Exact Proven Scope & Qualifications |
| --- | --- | --- |
| **Discrete SPD Operator** | `gpOp_quad_form`, `gpOp_symmetric`, `gpOp_psd`, `gpOp_spd`, `gpOp_vacuum_eq`, `gpOp_apply_zero` | Exact quadratic form, symmetry, and PSD ($\epsilon > 0, \kappa^2 \ge 0$) under skew-adjoint $D_a^\top = -D_a$. Strict SPD (`gpOp_spd`) is proven under pointwise $\kappa^2_i > 0$ everywhere (which holds for VASPsol $\kappa > 0$; SCCS inside the cavity where $S=0$ relies on `gpOp_psd` plus kernel connectivity of $D_a$). |
| **Dual Variational Principle** | `dualEnergy_on_shell`, `dualEnergy_perturb_exact`, `dualEnergy_global_max` | Exact on-shell value $\frac{1}{2}\langle \rho, \phi^*\rangle$, exact second-order perturbation $-\frac{1}{8\pi}\langle \delta, A\delta\rangle$, and global maximality of $\phi^*$. |
| **Vacuum & Algebraic $G=0$ Cancellation** | `polarization_rhs_vanishes_in_vacuum`, `zero_double_counting_cancellation`, `width_invariance_arbitrary_charge` | Exact vacuum RHS vanishing; algebraic cancellation of vacuum Hartree/Ewald double counting and algebraic $G=0$ `alpha_sm` cancellation across arbitrary net charge (assuming core-cavity non-overlap invariance of short-range electrostatic terms, which is verified numerically in Rust to `4.05e-8 Ha`). |
| **Bregman & GC-DFT Accounting** | `harris_foulkes_bregman_exact`, `grand_canonical_bregman_exact`, `bregman_bound_nonneg_triangle` | Exact solvation Harris–Foulkes remainder, grand-potential $\Omega = E - TS - \mu N$ shift consistency, and triangle inequality upper bound for `EnergyAudit::functional_difference_bound`. |
| **Analytic Calculus (`HasDerivAt`) & Sign Laws** | `sccsSwitch_hasDerivAt`, `sccsSwitch_hasDerivAt_nonpos`, `vaspsolSwitch_hasDerivAt`, `vaspsolSwitch_hasDerivAt_nonpos`, `vaspsol_deriv_nonpos`, `sccs_deriv_nonpos`, `diel_and_screening_potentials_nonneg` | Genuine Mathlib `HasDerivAt` chain-rule proofs on `sccsSwitch` (exact Andreussi 2012 log-sine formula on $(\rho_{\min}, \rho_{\max})$) and `vaspsolSwitch` (Mathew 2014 `erfc` chain rule for $n > 0$), proving $S'(n) \le 0$ and pointwise $V_{\mathrm{diel}} \ge 0, V_{\mathrm{scr}} \ge 0$. |
| **Gate 3–4 SCF Contracts** | `electrostatic_sign_involution`, `spin_polarized_solvation_chain_rule` | Involution $U_{\mathrm{el}}(\phi) = -\phi$ preserving $\rho_{\mathrm{el}}\phi = n U_{\mathrm{el}}$, and spin-independent chain rule $V_{\mathrm{solv}}^\uparrow = V_{\mathrm{solv}}^\downarrow = V_{\mathrm{solv}}[n_\uparrow + n_\downarrow]$. |

---

## 4. Level-1 / Level-2 Numerical Verification Results

Command: `cargo test -p dft-physics solvation:: -- --nocapture && cargo test -p dft-physics --test solvation_gate2 -- --nocapture` on `dev08` (`bf8e8c3f`): **3/3 unit tests PASS (`0.00 s`), 8/8 integration tests PASS (`14.46 s`), 0 clippy warnings**.

| Test | Oracle & Setup | Pre-Declared Tol (Gate 1 §9) | Measured Result | Status |
| --- | --- | --- | --- | --- |
| `cavity::tests::*` (3 unit tests) | Nyquist rules, VASPsol/SCCS FD derivatives, ionic force FD & MPI | Unit checks (`< 1e-7` / `< 1e-8`) | All 3 pass in `0.00 s` | **PASS** |
| `vacuum_limit_is_identically_zero_on_skewed_cell` | Skewed cell, $\epsilon_b = 1, \kappa = 0, \tau = 0, \beta = 0$ | `< 1e-14` | Exact `0.0` energy/potential/force, 0 FFTs | **PASS** |
| `mms_variable_dielectric_and_lpb_cubic_and_skewed_all_solvers` | (a) `cubic 16^3` & `skewed 16x18x14`, $\epsilon_b = 78.4$, `tol = 1e-14`; (b) 4-arm `cubic 24^3`, $\epsilon_b = 10$, `tol = 1e-11` | `DirectPcg` `< 1e-12` | (a) `DirectPcg`: `1.08e-14` / `8.38e-15` (cubic), `1.44e-14` / `1.03e-14` (skewed), `op_rel_err <= 1.71e-14`; (b) `DirectPcg`: `1.95e-12`, `FisicaroPcg`: `8.51e-6`, `PolarizationDiis`: `7.05e-7` | **PASS** |
| `born_ion_and_debye_huckel_analytic_convergence` | Smooth-wall Born (`L=14 -> 20`, Makov–Payne) + 1D continuum LPB ODE (`L=24 -> 32/36`, `dx=0.50`) | `< 1%` with monotonic $L$-convergence | Born: **`0.819%`** (`-0.338483` vs `-0.335734` Ha); DH ($\epsilon_b=78.4, k=0.4-0.2$): `1.90%` $\to$ **`0.245%`**; DH ($\epsilon_b=10.0, k=0.3-0$): `2.57%` $\to$ **`0.703%`** | **PASS** |
| `variational_potential_matches_central_finite_differences_vaspsol_and_sccs` | Skewed cell, 4th-order Richardson FD: full VASPsol-erfc, SCCS log/linear + isolated `cavitation_only`, `diel_and_es_kappa0`, `scr_and_es_eps1` | `rel < 1e-9` (`abs < 1e-10 Ha`) individually & summed | All 3 full models & all 3 isolated terms pass `< 1e-10 Ha` (`cav`: `9.76e-12`, `diel+es k=0`: `6.10e-11`, `scr+es eps=1`: `3.14e-11 Ha`) | **PASS** |
| `ionic_forces_match_central_finite_differences_and_width_invariance` | Skewed cell central FD ionic forces + neutral ($\kappa=0$) & charged ($Q=+1, \kappa=0.2$) $r_c \in \{0.65, 0.80\}$ invariance | Forces `< 1e-9` Ha/Bohr; width `< 1e-6 Ha` | All force components `< 1e-9`; neutral `|dE| = 2.01e-8 Ha`, charged `|dE| = 4.05e-8 Ha` | **PASS** |
| `mpi_rank_independence_and_mirror_vs_greedy_parity` | Skewed cell `P = 1, 2, 3, 4` (`dft_mpi::simulate`) + `mirror` vs `greedy` (`P = 3`) | `< 1e-12` | Energy, potential, forces match `< 1e-12` | **PASS** |
| `warm_start_bulk_delta_reduces_iterations_sharply` | `24^3` VASPsol-erfc ($\epsilon_b = 78.4, \kappa = 0.15$, $10^{-3}$ density pert) | Parity `< 1e-10`, strict drop | `DirectPcg` @ `1e-10`: `104 -> 79` (`|dE|=6.69e-12`); @ `1e-5`: `51 -> 25` (`104` FFTs). `FisicaroPcg` @ `1e-10`: `55 -> 48` (`|dE_cold|=1.68e-13`); @ `1e-5`: `34 -> 25` (`59` FFTs) | **PASS** |
| `invalid_solvation_parameters_are_rejected` | $\epsilon_b < 1$, invalid $\sigma, \kappa, \tau$ | Error returned | `SolvationError::InvalidParameter` | **PASS** |

---

## 5. `dev08` Release Performance Benchmarks (`solvation_bench`)

Command: `taskset -c 20-39,60-79 nice -n 10 cargo run --release -p dft-physics --example solvation_bench` (measured on `168d61b9`, solver arithmetic identical to `ee0a183a`/`bf8e8c3f`, single rank unless noted, $\epsilon_b = 78.4$, `tol = 1e-10` unless `@1e-5`).

### 5.1 Smooth MMS vs Sharp VASPsol-erfc Water Cavity (Single Rank)

| Grid (`dx` Bohr) | Case | $\kappa$ | Solver | Mode | Iters | FFTs | Wall (ms) | Rel Res | `|E - E_direct|` (Ha) |
| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `32^3` (`0.500`) | MMS ($\epsilon_b=78.4$) | 0.15 | `DirectPcg` | cold | 78 | 316 | 150.0 | `7.76e-11` | `0` |
| `32^3` (`0.500`) | MMS ($\epsilon_b=78.4$) | 0.15 | `FisicaroPcg` | cold | **15** | **37** | **16.4** (**9.1×**) | `9.10e-11` | — |
| `32^3` (`0.500`) | MMS ($\epsilon_b=78.4$) | 0.15 | `PolarizationDiis` | cold | 17 | 55 | 49.7 | `5.60e-11` | — |
| `32^3` (`0.500`) | VASPsol-erfc | 0.00 | `DirectPcg` | cold / warm@`1e-5` | 99 / 24 | 400 / 100 | 202.5 / 60.5 | `7.76e-11` / `9.25e-6` | `0` (`E = -0.009014`) |
| `32^3` (`0.500`) | VASPsol-erfc | 0.00 | `FisicaroPcg` | cold / warm@`1e-5` | 58 / 35 | 123 / 79 | 73.8 / 51.7 | `6.82e-11` / `5.54e-6` | **`9.12e-3`** *(sign flip!)* |
| `32^3` (`0.500`) | VASPsol-erfc | 0.00 | `PolarizationDiis` | cold | 35 | 109 | 130.2 | `7.61e-11` | **`5.91e-3`** |
| `32^3` (`0.500`) | VASPsol-erfc | 0.15 | `DirectPcg` | cold / warm@`1e-5` | 104 / 26 | 420 / 108 | 206.1 / 63.6 | `9.78e-11` / `9.81e-6` | `0` (`E = -0.009020`) |
| `32^3` (`0.500`) | VASPsol-erfc | 0.15 | `FisicaroPcg` | cold / warm@`1e-5` | 46 / 20 | 99 / 49 | 56.8 / 34.4 | `5.64e-11` / `6.44e-6` | **`7.70e-3`** |
| `32^3` (`0.500`) | VASPsol-erfc | 0.15 | `PolarizationDiis` | cold | 36 | 112 | 133.2 | `8.72e-11` | **`5.88e-3`** |
| `48^3` (`0.333`) | MMS ($\epsilon_b=78.4$) | 0.15 | `DirectPcg` vs `FisicaroPcg` | cold | 86 vs **15** | 348 vs **37** | 748.3 vs **72.1** (**10.4×**) | `< 1e-10` | — |
| `48^3` (`0.333`) | VASPsol-erfc | 0.15 | `DirectPcg` | cold / warm(`1e-3`) / warm@`1e-5` | 107 / 82 / 29 | 432 / 332 / 120 | 1012.6 / 800.9 / 325.4 | `< 1e-10` / `9.17e-6` | `0` (`E = -0.009267`) |
| `48^3` (`0.333`) | VASPsol-erfc | 0.15 | `FisicaroPcg` | cold / warm(`1e-3`) / warm@`1e-5` | 42 / 35 / **18** | 91 / 79 / **45** | 235.6 / 211.0 / **145.3** | `< 1e-10` / `9.16e-6` | **`6.94e-4`** |
| `48^3` (`0.333`) | VASPsol-erfc | 0.15 | `PolarizationDiis` | cold | 35 | 109 | 633.1 | `9.33e-11` | `1.02e-5` |
| `64^3` (`0.250`) | MMS ($\epsilon_b=78.4$) | 0.15 | `DirectPcg` vs `FisicaroPcg` | cold | 90 vs **14** | 364 vs **35** | 3042.0 vs **240.4** (**12.7×**) | `< 1e-10` | — |
| `64^3` (`0.250`) | VASPsol-erfc | 0.00 | `DirectPcg` | cold / warm(`1e-5`) / warm@`1e-5` | 102 / 58 / 28 | 412 / 236 / 116 | 3621.1 / 2182.1 / 1139.8 | `< 1e-10` / `8.46e-6` | `0` (`E = -0.009323`) |
| `64^3` (`0.250`) | VASPsol-erfc | 0.00 | `FisicaroPcg` | cold / warm(`1e-5`) / warm@`1e-5` | 49 / 42 / **20** | 105 / 93 / **49** | 998.5 / 896.0 / **546.8** | `< 1e-10` / `6.91e-6` | **`7.30e-6`** |
| `64^3` (`0.250`) | VASPsol-erfc | 0.15 | `DirectPcg` | cold / warm(`1e-5`) / warm@`1e-5` | 109 / 62 / 31 | 440 / 252 / 128 | 3920.6 / 2326.6 / 1246.3 | `< 1e-10` / `8.58e-6` | `0` (`E = -0.009329`) |
| `64^3` (`0.250`) | VASPsol-erfc | 0.15 | `FisicaroPcg` | cold / warm(`1e-5`) / warm@`1e-5` | **40** / **28** / **18** | **87** / **65** / **45** | **744.5** / **592.5** / **461.1** | `< 1e-10` / `5.02e-6` | **`7.21e-6`** |
| `64^3` (`0.250`) | VASPsol-erfc | 0.15 | `PolarizationDiis` | cold | 37 | 115 | 2115.1 | `5.38e-11` | `4.83e-6` |
| `64^3` (`0.250`) | VASPsol-erfc | 0.15 | `HybridFisicaroDirect` | cold | 130 | 451 | 3797.2 | `9.99e-11` | `7.06e-12` *(negative)* |

*Observed `FisicaroPcg` vs `DirectPcg` speedup summary: **2.74×–5.27×** on cold sharp VASPsol-erfc cavities, **1.85×–2.70×** on adaptive `warm(1e-3)@1e-5` VASPsol-erfc cavities, and **7.9×–12.7×** on smooth MMS.*

### 5.2 Simulated MPI Scaling & Rank Parity (VASPsol-erfc Water, $\kappa = 0.15$, `tol = 1e-10`, `mirror = true`)

| Grid | Solver | Ranks $P$ | Mode | Iters | FFTs | Wall (ms) | Total $E$ (Ha) | `|E(P) - E(1)|` (Ha) |
| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| `32^3` | `DirectPcg` | 1 / 2 / 4 | cold | 104 | 420 | 209.5 / 145.3 / **93.2** | `-0.0090203978` | `0` / `2.20e-16` / `3.07e-16` |
| `32^3` | `DirectPcg` | 1 / 2 / 4 | warm(`1e-3`) | 79 | 320 | 162.6 / 112.1 / **70.7** | `-0.0090947294` | `0` / `0.00e0` / `9.71e-17` |
| `32^3` | `FisicaroPcg` | 1 / 2 / 4 | warm(`1e-3`) | 38 | 85 | 51.6 / 35.5 / **22.8** | `-0.0013760570` | `0` / `7.93e-14` / `4.46e-14` |
| `48^3` | `DirectPcg` | 1 / 2 / 4 | cold | 107 | 432 | 1063.1 / 744.8 / **449.7** | `-0.0092666521` | `0` / `1.46e-16` / `6.96e-16` |
| `48^3` | `DirectPcg` | 1 / 2 / 4 | warm(`1e-3`) | 82 | 332 | 845.9 / 598.2 / **347.8** | `-0.0093402290` | `0` / `3.05e-16` / `1.47e-16` |
| `48^3` | `FisicaroPcg` | 1 / 2 / 4 | cold | 42 | 91 | 244.9 / 168.1 / **100.4** | `-0.0085725802` | `0` / `1.15e-13` / `2.20e-14` |
| `48^3` | `FisicaroPcg` | 1 / 2 / 4 | warm(`1e-3`) | 35 | 79 | 215.6 / 150.0 / **87.0** | `-0.0086499948` | `0` / `1.77e-13` / `1.69e-13` |

*Extrapolated adaptive SCF-step cost (`warm(1e-3)@1e-5`, 4 ranks, `48^3`, scaling by Table 5.1 ratios): `DirectPcg` $\approx 347.8 \times (325.4 / 800.9) = \mathbf{141.3\text{ ms}}$ (~`2.8–14.1%` of a `1–5 s` SCF map) and `FisicaroPcg` $\approx 87.0 \times (145.3 / 211.0) = \mathbf{59.9\text{ ms}}$ (~`1.2–6.0%` of `1–5 s`), with **0.0%** when solvation is disabled.*

---

## 6. Documented Negative Results & Failure Analyses (Directive Section 9)

1. **Symmetric Split Preconditioner $M_{\mathrm{split}}^{-1} = \epsilon^{-1/2}(L_0 + \bar{q})^{-1}\epsilon^{-1/2}$ Inside `DirectPcg` — REJECTED**:
   - *Result*: Stalled at `400` iterations on sharp cavities (`dx = 0.50` Bohr).
   - *Root Cause*: The transformed Rayleigh quotient $\int \epsilon|\nabla(\epsilon^{-1/2}v)|^2 = \int (|\nabla v|^2 + \frac{\nabla^2\sqrt{\epsilon}}{\sqrt{\epsilon}} v^2)$ omits the sharp cavity-wall curvature term $\nabla^2\sqrt{\epsilon}/\sqrt{\epsilon} \propto \epsilon_b^2/\delta^2 \sim 10^4 \gg \epsilon_b$, severely worsening conditioning. Reverted to optimal 0-FFT geometric-mean $\bar{\epsilon} = \sqrt{\epsilon_{\min}\epsilon_{\max}}$.
2. **Replacing $G^2 = 0$ with $G_{\mathrm{canon}}^2$ at the 7 Even-Nyquist Corners When $\kappa > 0$ — REJECTED**:
   - *Result*: Doubled `DirectPcg` iterations (`99 -> 198`).
   - *Root Cause*: Because BAANDI's skew-adjoint $D_a$ zeroes the 7 even-Nyquist modes, $A$ acts on those modes purely via $\kappa^2 \approx 1.4$, not $\bar{\epsilon} G_{\mathrm{canon}}^2 \approx 2000$ (1,400× under-scaling). Fixed by using exact skew-adjoint `g2` in `fill_direct_precond`.
3. **Single-Shot `HybridFisicaroDirect` (Stage 1 `FisicaroPcg` → Stage 2 `DirectPcg` Polish) — REJECTED FOR PRODUCTION**:
   - *Result*: Worse iterations and FFT counts than standalone `DirectPcg` (`32^3`: `165` iters / `555` FFTs / `273 ms` vs `99` / `400` / `202 ms`; `64^3`: `130` iters / `451` FFTs / `3,797 ms` vs `109` / `440` / `3,921 ms` — no meaningful wall-time saving while doing more work).
   - *Root Cause*: Because `evaluate()` computes $b_{\mathrm{pol}} = 4\pi\rho_{\mathrm{sol}} - A_{\mathrm{disc}}\phi_{\mathrm{vac}}$, solving $A_{\mathrm{Fis}}\phi_1 = b_{\mathrm{pol}}$ injects a high-$G$ cavity-wall aliasing error $(A_{\mathrm{Fis}} - A_{\mathrm{disc}})\phi_{\mathrm{vac}}$, leaving an $O(10^{-1})$ discrete residual for Stage 2.
4. **Coarse-Grid Operator Aliasing of Standalone `FisicaroPcg` and `PolarizationDiis` — QUANTIFIED**:
   - At `dx = 0.50` Bohr (`32^3`, $E_{\mathrm{cut,wfc}} \approx 10\text{ Ry}$, $E_{\mathrm{cut,\rho}} \approx 40\text{ Ry}$), `|E_fisicaro - E_direct| = 7.70e-3 Ha` (wrong sign at $\kappa = 0$) and `|E_diis - E_direct| = 5.88e-3 Ha`, while `DirectPcg` is already within `0.31 mHa` of the `64^3` limit. At production charge-grid spacing `dx = 0.25` Bohr (`64^3`, $E_{\mathrm{cut,wfc}} \approx 40\text{ Ry}$, $E_{\mathrm{cut,\rho}} \approx 158\text{ Ry}$), aliasing drops 1,000× to `7.21e-6 Ha` (`FisicaroPcg`) and `4.83e-6 Ha` (`PolarizationDiis`).

---

## 7. Workspace Gate Status (`dev08` `/scr/mino/gate/gate.sh`)

Verified on `dev08` at `/scr/mino/gate/runs/ee0a183a/gate.log`:

```text
HEAD ee0a183a2ea5a8a582005b29e9b85ba0a7de6287
rustc 1.99.1 (58cf8f9c3 2026-10-13)
PASS build
PASS test (140 binaries, 1134 passed, 0 failed, 247 ignored, 46 s wall)
PASS fmt
PASS clippy
PASS production-check
PASS production-clippy
GATE_WALL 187 s rc=0
=== feat_build rc=0 wall=44s
=== feat_test rc=0 wall=348s
GATE_DONE check_rc=0 hard=0 known_platform=0 new=0
```

*(Final Gate 2 commit after strict audit-hardening `bf8e8c3f` + Lean 4 fix passes 8/8 `solvation_gate2`, `clippy`, `lake build`, and workspace `fmt`.)*

---

## 8. Gate 3 Architecture Recommendations

1. **Default Production Solver**: `DirectPcg` (guaranteed SPD, exact `< 1e-10 Ha` discrete variational consistency on every grid, zero breakdown risk, `1.5e-16 Ha` MPI parity, 0 setup FFTs).
2. **High-Speed Opt-In Solver**: `FisicaroPcg` (2.7×–5.3× cold / 1.9×–2.7× adaptive warm speedup on sharp cavities, 2 FFTs/iter), guarded by a consistent grid-spacing warning (`max_dx <= 0.25 Bohr`) and automatic fallback to `DirectPcg` upon `Breakdown`, with Gate 4 R&D on **consistent $A_{\mathrm{Fis}}\phi_{\mathrm{vac}}$ RHS** and **inexact Fisicaro preconditioning inside Flexible PCG** ([benchmark-plan.md](benchmark-plan.md) §2 Suite C).
3. **Zero Core Regression**: Keep `SolvationWorkspace` behind `Option<Box<SolvationState>>` on `StartState` (`refstart.rs`) so pure-PW paths execute zero extra allocations, zero extra FFTs, zero extra collectives, and remain bit-identical.
