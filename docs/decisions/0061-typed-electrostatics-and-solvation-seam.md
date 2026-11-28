# ADR 0061: Typed Electrostatics and Implicit Solvation Seam (Track B Phase B1)

- Status: **Proposed** (Gate 3 Architecture Proposal, 2026-11-28; pending Gate 4 experimental SCF verification on research branch and Gate 5 production merge review).
- Scope: `dft-units`, `dft-physics::solvation`, `dft-engine` (`ScfOptions`, `refstart`, `refscf/{map,k_map,electronic_map,driver,k_driver,electronic,profile}`, `energy`, `forces`, `stress`, `force_identity`), `dft-cli` (`[solvation]` TOML input and additive `results.json` output).
- Relates to: `docs/research/extensible-dft/architecture-proposal.md`, `implicit-solvation-theory.md`, `feasibility-report.md`, `formal-verification-map.md`, `decision-log.md` (`D-001`–`D-014`); ADR 0017 (coupling sign/role contracts); ADR 0033 (derivative ledger and force/stress signs); ADR 0047 (single production reference driver); ADR 0049 (potential symmetry projection).

## Context

Track B Phase B1 introduces density-dependent continuum implicit solvation (Generalized Poisson and Linearized Poisson–Boltzmann) into BAANDI while establishing the shared electrostatic foundation for Phase B2 (grand-canonical DFT) and Track C (electrostatic QM/MM and polarizable MLIP embedding).

Two recurring failure modes in plane-wave electronic structure extensions must be prevented structurally:
1. **Electrostatic sign confusion**: Mixing electrostatic potential $\phi(\mathbf{r})$ (energy per $+1e$, Ha/$e$) with Kohn–Sham electron potential energy $U_{\mathrm{el}}(\mathbf{r}) = -\phi(\mathbf{r})$ (energy of $-1e$, Ha), which breaks variational consistency and double-counting cancellation.
2. **Pure-PW hot-path intrusion**: Modifying `DistributedRealHamiltonian::apply_block`, `DistributedComplexHamiltonian::apply_block`, `ref_davidson*`, frozen `oracle/`, or adding unconditional allocations/collectives when solvation is disabled.

## Decision

### 1. Compile-Time Electrostatic Potential Types (`dft-units`)
Add explicit scalar newtypes `ElectrostaticPotentialHaPerE(pub f64)` ($\phi$) and `ElectronPotentialEnergyHa(pub f64)` ($U_{\mathrm{el}} = -\phi$) with explicit `const fn` conversions (`to_electron_energy`, `to_electrostatic_potential`) and no implicit `From<f64>` or cross-type addition. All physics-to-engine boundaries convert explicitly.

### 2. Polarization-Potential Split and Exact $G=0$ Background Correction (`dft-physics::solvation`)
Keep `state.ionic`, `hartree.rs`, and `ewald.rs` unchanged. Solve for the polarization potential $\Delta\phi_{\mathrm{pol}} = \phi_{\mathrm{solv}} - \phi_{\mathrm{vac}}$ using BAANDI's skew-adjoint Nyquist-zeroed spectral derivative $D_a$ (`gga/distributed.rs`) for both $A = \sum_a D_a^\top \epsilon D_a + \epsilon_b\kappa^2\lambda$ and $L_0 = \sum_a D_a^\top D_a$, with exact background correction $E_{\mathrm{bg,corr}} = (Q_{\mathrm{ion}} - N_e)\sum_I \alpha_{Z,I}^{\mathrm{sm}}$ and consistent potential shift $-\sum_I \alpha_{Z,I}^{\mathrm{sm}}$.
- **Default solver**: `SolverKind::DirectPcg` (strictly SPD on every grid, exact `< 1e-9` discrete variational consistency, `1.5e-16` Ha MPI parity, 0 setup FFTs, `np = 2` real-space complex pairing + `mirror` unpacking with automatic `greedy` fallback).
- **Opt-in fast solver**: `SolverKind::FisicaroPcg` (2 FFTs/iter, 3.4×–12.7× faster on well-resolved grids), guarded by a preflight grid-spacing diagnostic warning when $\max_a \Delta x_a > 0.25\text{ Bohr}$ (`SolvationWorkspace::max_grid_spacing_bohr`) with automatic fallback to `DirectPcg` upon `Breakdown`.
- **Rejected arms**: Single-shot `HybridFisicaroDirect` and `M_split^{-1}` (Gate 2 negative results in `feasibility-report.md` §6), as well as Gate 4 Suite C (`crates/dft-physics/examples/solvation_suite_c.rs`: consistent Fisicaro RHS and inexact Fisicaro preconditioning inside Flexible PCG, which did not outperform standalone `DirectPcg`/`FisicaroPcg` enough to warrant production adoption), are not carried into production.

### 3. Production `refscf/` Seams and Zero Off-Path Regression (`dft-engine`)
- Attach `pub solvation: Option<Box<SolvationConfig>>` to `ScfOptions` (default `None`, boxed to keep `ScfOptions` compact) and `pub solvation: Option<Box<SolvationState>>` to `StartState` (shared by `KStartState` via `DerefMut`).
- Evaluate solvation inside `refscf::map::fields()` reusing `n_g` (`coeff`) before Hartree overwrites it (saving 1 3D FFT per map), fold `V_solv` into `hxc` before `potential_symmetry` (ADR 0049), and in spin-polarized `electronic_map::fields()` add spin-independent `V_solv[n_tot]` to both spin channels.
- Extend `EnergyAudit` with `pub solvation: Option<SolvationEnergyLedger>` and exact Harris–Foulkes / Bregman remainder accounting; when `None`, `total_hartree()` executes the identical floating-point summation as pure PW.
- Ratchet inner PDE tolerance adaptively via `refscf::profile::next_solvation_tol`, tightening to the exact target tolerance on the final verification map and derivative evaluation.
- Record solvation parameters in `force_identity::ModelDescription`.

### 4. Forces and Fail-Closed Stress (`forces.rs`, `stress.rs`)
- Keep base `derivative_ledger()` unchanged (exact 8 NC terms for pure-PW runs) and append `"solvation"` (`force_available = true`, `stress_available = false`) inside `forces::validate()` only when solvation is active.
- Compute analytic Gaussian ionic solvation forces (plus NLCC core-cavity forces when core density enters the cavity) on the converged density and include them in `ForceResult::solvation` and `ForceResult::values`.
- Fail closed in `stress::preconditions()` and `stress::frozen_pass()` whenever `options.solvation.is_some()` with `"observables.stress: incomplete derivative coverage: solvation"` until analytic solvation stress is implemented in a subsequent phase.

### 5. Additive CLI Input and Output (`dft-cli`)
- Add optional `[solvation]` TOML section (`#[serde(deny_unknown_fields)]`), rejected cleanly on frozen `oracle` and legacy spin/finite-T paths.
- Emit `"solvation"` blocks in `results.json` (`energy`, top-level telemetry, `forces.per_term`, and `forces.derivative_ledger`) only when active.

## Consequences

- Zero changes to `ownership.toml`, workspace `Cargo.toml`, `Distributed*Hamiltonian::apply_block`, or `ref_davidson*`, and zero behavioral changes to frozen `oracle/` (only two mechanical `solvation: None` struct-literal completions for the new `EnergyAudit` field).
- Bit-identical pure-PW behavior and `results.json` output when `[solvation]` is absent.
- Directly reusable electrostatic types, $\phi_{\mathrm{bulk}} = 0$ gauge, and Bregman audit seams for Track B Phase B2 (GC-DFT) and Track C (electrostatic QM/MM and polarizable MLIP).
- Core Level-1 mathematical guarantees are formally machine-checked under `formal/BaandiMath/` via Lean 4 + Mathlib (`bash scripts/check_lean.sh`, zero `sorry`, standard axioms only), with hypothesis-by-hypothesis Rust traceability in `docs/research/extensible-dft/formal-verification-map.md`.
