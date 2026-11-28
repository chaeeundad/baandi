# Extensible BAANDI Research Program — Decision Log

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Tracks A, B, C)
- **Location**: `docs/research/extensible-dft/`
- **Governance (Directive Section 10.1)**:
  - Literature review, mathematical derivation, isolated numerical prototypes, research branches, benchmarks, profiling, and failure-driven alternative exploration run **autonomously** within an approved research scope.
  - **Research-lead approval / review is required ONLY for**: (1) initial track selection (Gate 0 → Gate 1 entry), (2) major scope changes, (3) major core-architecture changes, and (4) production integration (Gate 5).
  - Within an approved scope, Gate-to-Gate progression does **not** require separate human sign-off at every Gate boundary; every Gate decision (`Go / Revise / Hold / Stop`), tolerance, and failure analysis is logged here transparently.

---

## Standing Program Rules

1. **No Unverified Citations (Section 2.1)**: Only papers and DOIs directly verified via primary lookup may be cited; any unverified placeholder must be marked `[Citation needed]`.
2. **Zero Pure-PW Performance Regression (Sections 2.3 & 6.4)**: Extensions must not add allocations, I/O, branches, or extra MPI collectives inside `DistributedRealHamiltonian::apply_block`, `DistributedComplexHamiltonian::apply_block`, `ref_davidson`, or `ref_davidson_complex`. When disabled (`Option::None`), pure-PW execution must remain bit-identical.
3. **Production Contract Preservation (ADR 0047)**: All production extensions target `refscf/` and `dft-physics` (or dedicated opt-in crates). Frozen `oracle/` default arithmetic is never altered, and `cargo check/clippy -p dft-engine --no-default-features` must stay green.
4. **No Premature Universal Abstractions (Section 6.6)**: Do not introduce a monolithic multi-track "universal engine trait" or large core refactoring before concrete vertical implementations prove the minimal shared boundary.
5. **Independent Verification & Licensing Hygiene (Sections 7 & 9)**: Independent verification (never relying solely on the implementing agent's self-check); never relax tolerances post-hoc to manufacture a pass; never copy GPL/LGPL source code (e.g., QE, Environ, Wannier90, ABINIT, Yambo [GPL]; ASE [LGPL-2.1+]). Implement physics kernels cleanly from published mathematical formulations under `[workspace.lints.rust] unsafe_code = "deny"`.
6. **Formal Mathematical Verification via Lean 4 + Mathlib (User Rule, 2026-11-28; `AGENTS.md`)**: Formalize and machine-check core Level-1 mathematical theorems (variational stationarity, operator SPD/self-adjointness, sign/gauge and zero-double-counting identities, limiting-case reductions) under `formal/` using Lean 4 + Mathlib (`lake build`, zero `sorry`, zero custom axioms). Formal proofs complement Level 2–5 Rust numerical verification.

---

## Decision Entries

### D-001 — Gate 0 Scope Boundary (2026-11-28)
- **Status**: Executed (independently audited by `Gate 0 Quality Auditor` subagent)
- **Decision**: Execute Gate 0 strictly as repository assessment, codebase/literature triage answering all 10 Directive Section 12 questions, and track prioritization without production Rust code changes, large refactoring, or premature 12-item theory specifications across all three tracks.
- **Evidence**: [repository-assessment.md](repository-assessment.md), [triage-report.md](triage-report.md).

### D-002 — Anchor New Seams in `refscf/` Rather Than Frozen `oracle/` (2026-11-28)
- **Status**: Adopted (consequence of ADR 0047)
- **Context**: Nov 2026 coupling documents (`docs/architecture/scientific-coupling.md`, `coupling-codex-review.md`, ADR 0017) referenced `ScfProblem` in `scf.rs`, which is now frozen under `oracle/` (`strict-oracle` feature).
- **Decision**: All future local-potential, solvation, external-field, and grand-canonical hooks target the production `refscf/{api, driver, k_driver, electronic, map, k_map, electronic_map}` and `refstart.rs` seams, preserving ADR 0017's sign/role contracts while obeying ADR 0047's single-path rule.

### D-003 — Retain Typed Electrostatic Semantics from ADR 0017 as Shared B/C Foundation (2026-11-28)
- **Status**: Proposed for Gate 1
- **Decision**: Implement explicit semantic types distinguishing electrostatic potential $\phi$ (Ha/e), electron potential energy $U_{\mathrm{el}} = -\phi$ (Ha), positive electron density $n$, and signed charge density $\rho$ before wiring either Track B solvation or Track C electrostatic embedding.

### D-004 — Defer Track A (Mixed Basis / PW+AO / ALB) Behind Track B & Track C (2026-11-28)
- **Status**: Recommended (pending user approval)
- **Rationale**: Track A requires replacing BAANDI's orthonormal PW eigenproblem ($S = I$) and `KeptHx` warm-start architecture with a generalized $H\psi = \epsilon S\psi$ solver, cross-density quadrature, and Pulay forces/stress. Doing Track A first carries the highest core-intrusion and performance-regression risk with the lowest cross-track reuse.

### D-005 — Separate Mechanical vs Electrostatic QM/MLIP Layers (2026-11-28)
- **Status**: Recommended (pending user approval)
- **Rationale**: Subtractive/ONIOM/$\Delta$-ML QM/MLIP (e.g., BAANDI + SevenNet/MACE) requires zero Rust SCF core changes and belongs in the adapter/workflow layer (`baandi-py` / `python/baandi/ase.py`), whereas electrostatic QM/MM and polarizable QM/MLIP reuse the exact same typed `ExternalEnvironment` seam as Track B.

### D-006 — Gate 1 Entry Track Selection (2026-11-28)
- **Status**: **Approved by Research Lead (2026-11-28)** — Option 1 selected, with explicit directive to engineer high performance directly into new paths rather than merely defending against regression.
- **Approved Scope**: Track B Phase B1 (Typed Electrostatic Seam + Density-Dependent Implicit Solvation via Generalized Poisson / LPB, with Phase B2 GC-DFT and Track C electrostatic readiness). Autonomous execution through Gate 1 → Gate 2 per Directive Section 10.1 using parallel subagents.

### D-007 — Polarization-Potential Split & Exact $G=0$ Background Correction (2026-11-28)
- **Status**: Adopted (Gate 1)
- **Decision**: Keep BAANDI's `state.ionic`, `hartree.rs`, and `ewald.rs` unmodified; formulate solvation via the dual variational polarization-potential split $\Delta\phi_{\mathrm{pol}} = \phi_{\mathrm{solv}} - \phi_{\mathrm{vac}}$ plus exact background-uncompensation correction $E_{\mathrm{bg,corr}} = +(Q_{\mathrm{ion}} - N_e)\sum_I \alpha_{Z,I}^{\mathrm{sm}}$ (active gauge shift $-\sum_I \alpha_{Z,I}^{\mathrm{sm}}$ when $\kappa > 0$), using BAANDI's skew-adjoint Nyquist-zeroed derivative $D_a$ (`gga/distributed.rs`) for both $A$ and $L_0$.
- **Rationale**: Proves zero double counting, exact $r_{c,I}$-invariance for neutral and charged cells, exact bit-identical vacuum reduction at $\epsilon_b = 1, \kappa = 0, \tau = 0$, and exact $\phi_{\mathrm{bulk}} = 0$ potential alignment for GC-DFT. See [implicit-solvation-theory.md](implicit-solvation-theory.md#4-mathematical-derivation-exact-baandi-conventions).

### D-008 — Primary Solver Selection: Direct G-Space PCG with Real-Space & Mirror Pairing (2026-11-28)
- **Status**: Adopted for Gate 2 Primary Arm (with Alg 2 & Alg 3 benchmarked alongside)
- **Decision**: Select Direct G-space PCG (`SolverKind::DirectPcg`) with `np = 2` real-space complex pairing and BAANDI's default `mirror` stick unpacking (4 FFTs/iter, 0 setup FFTs/SCF step, strictly SPD, $\kappa_2(M^{-1}A) \le \epsilon_b$), plus 0-FFT bulk-delta warm-starting and adaptive inner tolerance. Implement `FisicaroPcg` (Alg 2) and `PolarizationDiis` (Alg 3) in the Gate 2 prototype for empirical head-to-head comparison on sharp $\epsilon_b = 78.4$ cavities.
- **Evidence**: [algorithm-comparison.md](algorithm-comparison.md).

### D-009 — Gate 2 Prototype Isolation Boundary (2026-11-28)
- **Status**: Adopted (Gate 1 → Gate 2)
- **Decision**: Develop the Gate 2 numerical prototype on an independent research branch (`research/track-b-gate2-proto`) inside an isolated submodule `crates/dft-physics/src/solvation/` with zero changes to `ownership.toml`, `Cargo.toml`, `dft-solvers`, or `dft-engine` (preserving pure-PW paths and `--no-default-features` green).

### D-010 — External Solvation Oracle Status on `dev08` (2026-11-28)
- **Status**: Recorded (Gate 1)
- **Fact**: Direct inspection of `dev08` (`/scr/mino/oracles/`) confirmed QE 7.6 (`environ_api_stub.f90` only — Environ not compiled), ABINIT 10.8.3, and Yambo 5.3.0; JDFTx and GPAW are not installed. Gate 2 verification relies on independent L1/L2 exact mathematical, MMS, Born/Debye–Hückel, finite-difference, and MPI-parity oracles with pre-declared tolerances in [implicit-solvation-theory.md](implicit-solvation-theory.md#9-verification-plan-levels-15-pre-declared-tolerances).

### D-011 — Adopt Lean 4 + Mathlib for Level-1 Mathematical Verification (2026-11-28)
- **Status**: Adopted & Verified (User Directive, recorded in `AGENTS.md`; hardened in `D-014`)
- **Decision**: Formalize Level-1 mathematical identities (variational stationarity, discrete skew-adjoint SPD operator $A = \sum_a D_a^\top \epsilon D_a + \kappa^2 \succeq 0$, analytic `HasDerivAt` cavity switching derivatives, algebraic $G=0$ double-counting/width-invariance cancellation, and Bregman/spin-chain-rule identities) in Lean 4 + Mathlib under `formal/BaandiMath/`.
- **Evidence**: `elan` + Lean `v4.36.0` + Lake `v5.0.0` provisioned on `dev08`; the initial Gate 2 proof suite in `formal/BaandiMath/BaandiMath/Solvation.lean` passed `lake build` (`exit 0`, zero `sorry`, `#print axioms` strictly `[propext, Classical.choice, Quot.sound]`). Subsequent critical review (`D-014`) upgraded the SPD and PDE width-invariance proofs and eliminated theorem-count progress reporting.

### D-012 — Track B Gate 2 Exit Decision: `Revise` Cycle Completed → Go to Gate 3 with Dual-Solver Strategy (2026-11-28)
- **Status**: Executed (Autonomous within approved Track B scope per Directive Section 10.1; independently audited by `Gate 2 Quality Auditor`)
- **Decision**: **Go → Gate 3 (Architecture Proposal)** with explicit waiver of Gate 1 Criterion 2–3 standalone iteration targets (`<= 35` cold / `<= 10` warm) following a completed `Revise` exploration (`M_split^{-1}` and single-shot `HybridFisicaroDirect` negative results). Adopt `DirectPcg` as the default exact discrete variational solver (strictly SPD, `< 1e-10 Ha` FD variational derivatives on every grid, `1.5e-16 Ha` MPI parity, `0.3 mHa` grid accuracy already at `32^3`) and retain `FisicaroPcg` as a grid-resolution-guarded (`max_dx <= 0.25 Bohr`) high-speed opt-in accelerator (2.74×–5.27× cold / 1.85×–2.70× adaptive warm speedup on sharp cavities, 7.9×–12.7× on smooth MMS, 2 FFTs/iter, automatic `DirectPcg` fallback on breakdown). Real-SCF-map `< 5%` wall-time verification and Flexible-PCG acceleration R&D are assigned to Gate 4.
- **Evidence**: [feasibility-report.md](feasibility-report.md), `research/track-b-gate2-proto` (`ee0a183a` full `dev08` `gate.sh`: `GATE_DONE check_rc=0 hard=0 known_platform=0 new=0`, `140 binaries, 1134 passed`; strict audit-hardened `bf8e8c3f`: initial Lean 4 suite, 3/3 unit + 8/8 independent L1/L2 tests passing all exact pre-declared tolerances with zero relaxation).

### D-013 — Track B Gate 3 Architecture Approval (Autonomous Progression to Gate 4 Research Branch) (2026-11-28)
- **Status**: Executed (Autonomous within approved Track B scope per Directive Section 10.1; Gate 5 production merge still requires Research Lead sign-off)
- **Decision**: Advance to **Gate 4 (Experimental SCF Integration)** on isolated research branch `research/track-b-gate4-scf` following [architecture-proposal.md](architecture-proposal.md), [verification-plan.md](verification-plan.md), [benchmark-plan.md](benchmark-plan.md), and draft [ADR 0061](../../decisions/0061-typed-electrostatics-and-solvation-seam.md). Enforce zero changes to `ownership.toml`, `Cargo.toml`, `Distributed*Hamiltonian::apply_block`, `ref_davidson*`, and `oracle/`, fail-closed Phase B1 stress, and bit-identical pure-PW off-path execution.

### D-014 — Critical Formal-Verification Hardening, Policy Refinement & Gate 3–4 Doc/ADR Reconciliation (2026-11-28)
- **Status**: Adopted (Gate 4 research branch `research/track-b-gate4-scf`)
- **Context**: A critical review of the Gate 2–3 Lean 4 formalization and Gate 3–4 documentation identified three substantive gaps: (1) the initial `gpOp_spd` theorem assumed pointwise `∀ i, 0 < k2 i`, which does not cover pure GP (`κ = 0`) or realistic LPB cavities where `k2 = 0` inside the solute; (2) `width_invariance_arbitrary_charge` proved only algebraic `G = 0` offset cancellation for a fixed `E_es_solv`, not full PDE electrostatic functional width invariance; and (3) minor doc/ADR/code drifts existed around the `FisicaroPcg` grid-spacing threshold (`0.28` vs `0.25 Bohr`), boxed `Option<Box<SolvationConfig>>`, `next_solvation_tol` formula, and Suite C negative-result recording.
- **Decision**:
  1. **Upgrade `formal/BaandiMath/BaandiMath/Solvation.lean`**: Prove exact kernel characterization (`gpOp_quad_eq_zero_iff`), pure-GP `ker(D)`-projected SPD (`gpOp_spd_on_grad_complement`), cavity-vanishing LPB SPD (`gpOp_spd_lpb_cavity_vanishing`), single-bulk-point constant-mode LPB positivity (`gpOp_spd_const_mode_lpb`), Fourier-multiplier Nyquist-zeroed skew-adjointness (`nyquist_zeroed_multiplier_skew_adjoint` + `nyquist_nonzero_breaks_reality`), honest rename to `bg_offset_algebraic_cancellation`, and the genuine PDE-level Green's reciprocity invariance theorem `full_electrostatic_width_invariance`.
  2. **Publish explicit assumption mapping & reproducible CI check**: Create [formal-verification-map.md](formal-verification-map.md) mapping every Lean theorem and hypothesis to exact Rust enforcement lines and documenting the `ℝ` vs. `f64`/FFT/solver-tolerance boundary; add `scripts/check_lean.sh` enforcing comment-aware `sorry`/`admit`/custom-`axiom` scanning, `lake build`, and `#print axioms` allowlist checking.
  3. **Refine `AGENTS.md` policy**: Ban theorem/proof count as a research-progress metric; mandate theorem-to-Rust assumption mappings; scope **Aeneas** Rust-to-Lean extraction strictly as an optional small-kernel research study (never blocking engine development or affecting hot-path performance).
  4. **Reconcile docs & ADR 0061**: Standardize `0.25 Bohr` for `FisicaroPcg` preflight warning, `Option<Box<SolvationConfig>>`, the `next_solvation_tol` formula `clamp(min(prev, 1e-2 * sqrt(max(dr2_ry, 0))), target, max(target, 1e-4))`, and record the Gate 4 Suite C negative result in [ADR 0061](../../decisions/0061-typed-electrostatics-and-solvation-seam.md).
