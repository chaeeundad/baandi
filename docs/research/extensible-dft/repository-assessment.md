# BAANDI Repository Assessment (Gate 0)

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Tracks A, B, C)
- **Phase**: Gate 0 — Repository Assessment & Triage (analysis only; no production code or refactoring)
- **Date**: 2026-11-28
- **Assessed Baseline**: `main` at `24d48b65` (includes `RPA-R4A-ATOMIC-PROJECTION` and `RPA-R4B-WANNIER90`, verified via `dev08` gate `a9e50e58`)
- **Evidence Convention**: All file paths, structs, functions, and invariants below were directly inspected in this checkout. Pre-Stage-R coupling documents (ADR 0017, ADR 0027, `docs/architecture/scientific-coupling.md`) are explicitly reconciled against the post-Stage-R codebase.

---

## 1. Workspace Topology and Ownership Boundaries

BAANDI is a 15-crate Rust workspace (`Cargo.toml`, Rust 2024 edition, pinned toolchain `1.99.1`) enforcing `[workspace.lints.rust] unsafe_code = "deny"` (`Cargo.toml:29`). Unsafe FFI code is strictly confined to three adapter crates:
- `dft-mpi` (`crates/dft-mpi/Cargo.toml` sets `unsafe_code = "allow"`; `crates/dft-mpi/src/ffi.rs` calls C MPI via `unsafe extern "C"`).
- `dft-mkl` (link-time `unsafe extern "C"` FFI configured via `crates/dft-mkl/build.rs` under Cargo feature `mkl`).
- `dft-dispersion` (runtime `dlopen` via `libloading` in `crates/dft-dispersion/src/ffi.rs` under Cargo feature `dftd3`, ADR 0048).

`dft-numerics` inherits `unsafe_code = "deny"` and contains zero `unsafe` in `src/`. Dependency direction is strictly acyclic and governed by `ownership.toml`:

```text
dft-units
  ├── dft-upf
  └── dft-geometry
        └── dft-basis  (+ dft-numerics [+ optional dft-mkl])
              └── dft-solvers  (+ dft-mpi, dft-numerics)
                    └── dft-physics  (+ dft-upf, dft-basis, dft-geometry, dft-units, dft-numerics, dft-solvers, dft-mpi)
                          ├── dft-response  (+ dft-geometry, dft-basis, dft-numerics, dft-physics, dft-upf, dft-mpi)
                          └── dft-engine    (+ dft-response, dft-physics, dft-solvers, dft-basis, dft-geometry, dft-upf, dft-units, dft-numerics, dft-mpi)
                                ├── dft-cli        (+ optional dft-dispersion via dlopen)
                                ├── baandi-py      (PyO3 + python/baandi/ase.py)
                                └── baandi-mcp     (agent-native MCP adapter)
```

Key architectural properties verified in code:

1. **Single production SCF contract (ADR 0047)**: `crates/dft-engine/src/refscf/api.rs` (`solve`, `stepper`, `Observables`) is the sole production entry point across CLI, Python/ASE, and MCP. The frozen strict oracle lives under `crates/dft-engine/src/oracle/` behind Cargo feature `strict-oracle` and is isolated by `cargo check/clippy -p dft-engine --no-default-features` in `scripts/check.sh`.
2. **Strict unit and sign boundaries**: All internal physics operates in Hartree atomic units (Bohr, Ha, Ha/Bohr, Ha/Bohr^3, positive electron number density $n(\mathbf{r}) \ge 0$). Unit conversions live exclusively in `dft-units` at CLI/Python/UPF boundaries, enforced by `scripts/tests/test_no_hardcoded_units.py`.
3. **Optional external-library pattern (ADR 0048)**: `crates/dft-dispersion/src/lib.rs` establishes BAANDI's boundary for optional LGPL/third-party scientific libraries: runtime `dlopen`, zero overhead when absent, typed availability/error handling, and additive composition with energy, forces, stress, and relaxation.
4. **No `dft-core` crate**: Pure mathematical kernels belong in `dft-physics` (with FFT/linear primitives in `dft-numerics`/`dft-solvers`), SCF orchestration and energy/derivative accounting in `dft-engine`, and I/O, serialization, and external workflow drivers in `dft-cli`, `baandi-py`, or dedicated adapter crates.

---

## 2. Production SCF, Hamiltonian, and Derivative Architecture

### 2.1 Three Specialized Production SCF Paths

`crates/dft-engine/src/refscf/api.rs` dispatches into three rank-distributed SCF drivers sharing `refstart.rs`, `refscf/broyden.rs`, `refscf/metric.rs`, and `refscf/map.rs::cached_fields`:

| Path | Driver & Map Modules | Hamiltonian Struct | Basis / State Representation | Scope |
| --- | --- | --- | --- | --- |
| **Gamma real** | `refscf/driver.rs`, `refscf/map.rs` | `dft_physics::distributed_real_gamma::DistributedRealHamiltonian` | Half-sphere real ORC (`2*nh - 1` real rows), paired band FFT (`DistributedPairFft`) | Spin-unpolarized, fixed even $N$, $\Gamma$-only |
| **Complex $k$-point** | `refscf/k_driver.rs`, `refscf/k_map.rs` | `dft_physics::distributed_complex::DistributedComplexHamiltonian` | Full-sphere complex Bloch basis per $k$, optional $k$-pools (`refscf/k_pool.rs`) | Spin-unpolarized, fixed even $N$, multi-$k$ + symmetry |
| **Spin / Fermi–Dirac** | `refscf/electronic.rs`, `refscf/electronic_map.rs` | `DistributedComplexHamiltonian` (per spin channel & $k$) | 1 or 2 collinear spin channels, FD occupations (`dft_physics::occupations::FermiDirac`), free or fixed moment, automatic band growth (`BandGrowth`) | Metals (FD smearing), collinear spin, NLCC, PBE/LDA |

*(Note: Free-moment collinear FD with PBE currently falls back to `StrictOracle` per ADR 0047 `G4D_FALLBACK`.)*

### 2.2 Exact Local Potential, Poisson, and Energy Accounting Seams

Across all three production paths, the effective local potential is assembled on rank-owned real-space planes of the density FFT grid (`dims`), distributed through `StartFields` (`crates/dft-engine/src/refstart.rs:217–220`, which wraps `fft: DistributedPairFft` and implements `dft_physics::gga::distributed::Fields<K>`):

1. **Ionic local potential (`state.ionic`)**: Precomputed once per geometry in `crates/dft-engine/src/refstart.rs` from norm-conserving UPF local transforms (`dft_upf::radial`), retaining the finite Coulomb-subtracted $G=0$ mode `alpha_z()`.
2. **Hartree + XC assembly (`refscf/map.rs::fields` and `refscf/electronic_map.rs::fields`)**:
   - Transforms rank-owned `density` to reciprocal space via `state.fields.forward`.
   - Computes Hartree directly on owned G-sticks: $V_H(\mathbf{G}) = 4\pi n(\mathbf{G}) / |\mathbf{G}|^2$ for $\mathbf{G} \ne \mathbf{0}$, with $V_H(\mathbf{0}) = 0$ hard-coded (`crates/dft-physics/src/hartree.rs` and `refscf/map.rs`).
   - Evaluates LDA/PBE XC (with optional NLCC `state.core_density` / `state.pbe_core`, PBE gradient via `StartFields`, and optional `potential_symmetry` G-space projection per ADR 0049).
   - Forms `potential = state.ionic + input_fields.hxc`, computes global mean `v0` for the Davidson diagonal preconditioner (`diagonal_static + v0`), redistributes planes from field layout to operator layout via `crate::grid_dist::PlaneMap` (or `k_pool::FullField`), and installs it via `operator.set_owned_potential(&op_potential, epoch)`.
3. **Energy accounting (`crates/dft-engine/src/energy.rs::EnergyAudit`)**:
   - Tracks `kinetic`, `local` (including `local_zero_mode`), `nonlocal`, `hartree`, `xc`, `ewald`, `band_sum`, `double_counting`, `band_form`, `harris_foulkes`, and `functional_difference_bound`.
   - `EnergyAudit::total_hartree()` sums `kinetic + local + nonlocal + hartree + xc + ewald`.
   - For finite-temperature FD runs, `refscf/electronic.rs` additionally maintains `dft_engine::model::ThermodynamicLedger` (defined in `crates/dft-engine/src/model/spin.rs`: $E$, $-TS$, $F = E - TS$) and `occupations::Evaluation` ($\mu$, $N$, entropy, edge coverage).
   - **Critical Gap**: Neither `Fields`, `StartState`, `KStartState`, `EnergyAudit`, nor `ScfOptions` currently carries an `ExternalFieldTerm`, generalized dielectric/solvation potential $V_{\mathrm{solv}}[n]$, or charged-cell / grand-canonical energy ledger entry.

### 2.3 Hard-Coded Neutrality and Electron-Count Assumptions

Direct inspection of `crates/dft-engine/src/refstart.rs` reveals explicit charge-neutrality enforcement:

- `refstart.rs:509–515` (SCF startup): `let electrons: f64 = ions.iter().map(|i| i.charge).sum();` followed by strict integer check `(electrons - electrons.round()).abs() > 1e-10` and even-count enforcement when `!general_occupations`.
- `refstart.rs:2746–2813` (`start_kpoints_from_density`, fixed-density NSCF startup): Re-enforces the same neutrality and integer-count invariant (`"reference NSCF: invalid electron count"`). (SCF ionic warm starts separately use `IonicSeed` in `refscf/delivered.rs`.)
- `crates/dft-physics/src/ewald.rs`: Computes point-charge Ewald energy/forces/stress with a uniform neutralizing background $-\pi Q_{\mathrm{ion}}^2 / (2 \Omega \eta^2)$ tied solely to $\sum Z_I$.
- `crates/dft-physics/src/occupations.rs::FermiDirac::evaluate`: Solves $N(\mu) = N_{\mathrm{target}}$ for fixed `target: f64` by monotone bracketing + bisection (`electron_tolerance <= 1e-13`), rejecting endpoints (`wanted <= 0` or `wanted >= capacity`). No fixed-$\mu$ evaluator exists yet, though `fraction_entropy((eps - mu)/theta)` (`occupations.rs:59–66`) is already cleanly factored and overflow-free.

### 2.4 Forces, Stress, Relaxation, Subspace, and External Adapters

- **Forces (`crates/dft-physics/src/forces.rs`, `crates/dft-engine/src/forces.rs`, `refscf::api::Observables::{evaluate_forces, evaluate_forces_electronic}`)**: Post-convergence Hellmann–Feynman evaluator combining local structure-factor derivative (`local_or_core_cached`), NLCC core derivative, nonlocal Kleinman–Bylander derivative, Ewald derivative, and optional post-SCF `dft-dispersion` (D3). Requires verified convergence (`verify = "on"`).
- **Stress (`crates/dft-physics/src/stress.rs`, `crates/dft-engine/src/stress.rs`, `refscf::api::Observables::{evaluate_stress, evaluate_stress_electronic}`)**: Analytic frozen-Miller strain derivatives (`hartree`, `local_or_core_cached`, kinetic, nonlocal, XC, Ewald, optional D3) with opt-in Pulay cutoff scan (`PulayAssessmentMode`) and joint force+stress pass (`forces_and_frozen_stress`).
- **Relaxation**: Fixed-cell and variable-cell BFGS / preconditioned L-BFGS (ADR 0058) / FIRE in `crates/dft-cli/src/relax.rs`, with checked density warm-transfer and symmetry preservation.
- **Subspace & Atomic Orbitals (`dft-upf/src/radial.rs`, `refstart.rs`, `dft-response/src/{subspace.rs, wannier.rs}`)**: UPF `PSWFC` radial wavefunctions are parsed, Bessel-transformed via `RadialTransformTable`, projected onto PWs for atomic start (`refstart.rs`), Löwdin-orthonormalized (`RPA-R4A`, ADR 0059), and interfaced with Wannier90 `.win/.eig/.amn/.mmn/_u.mat` (`RPA-R4B`, ADR 0060).
- **Python / ASE (`crates/baandi-py/src/lib.rs`, `python/baandi/ase.py`)**: Wraps `Session` for **3D-periodic, neutral, spin-unpolarized, $\Gamma$-only** energy, forces, and stress (`ase.py:52–58` explicitly rejects non-periodic cells, initial magnetic moments, and initial charges; `use_cache = False`).

---

## 3. Performance-Critical Hot Paths and Invariants to Preserve

BAANDI's performance engineering (P10–P48, P3–P4, GK1–GV1) targets $\le 0.90\times$ Quantum ESPRESSO wall time at equal accuracy on multi-rank MPI + MKL (`dev08`, 40-core Xeon Gold 6230). Extensions must preserve:

1. **Zero-allocation, zero-I/O operator apply**:
   - `DistributedRealHamiltonian::apply_block` (`distributed_real_gamma.rs`) pre-allocates `hbuf`, `xsend`, `xrecv`, `hsend`, `projections`, and `gemm`.
   - `DistributedComplexHamiltonian::apply_block` (`distributed_complex.rs`) pre-allocates `gemm`, `projection`, `coupled`, and `back`.
   - Neither path permits heap allocation, logging, external I/O, or rank-divergent NaN checks (which would deadlock MPI collectives) inside `apply_block` or Davidson iterations.
2. **Fixed MPI collective budget**: `DistributedRealHamiltonian::apply_block` executes a fixed, rank-independent collective count (`2 + 2` batched `alltoall` + `1` `tree_sum_f64`). Extensions must not add collectives inside `apply_block`.
3. **Orthonormal PW Davidson & `KeptHx` warm-start invariant**: `ref_davidson` and `ref_davidson_complex` assume $S = I$ on PW rows and exploit `KeptHx` ($H_{\mathrm{new}} X = H_{\mathrm{old}} X + \Delta V X$ via `owned_epoch` in `refscf/{kept_hx.rs, k_kept.rs}`). Additive scalar local potentials ($V_{\mathrm{solv}}(\mathbf{r})$ or $U_{\mathrm{ext}}(\mathbf{r})$) preserve both invariants automatically with **zero overhead inside Davidson**; non-orthogonal bases (Track A) break both unless isolated.
4. **Bit-identical, zero-overhead off-path**: When an extension is disabled (`Option::None`), pure-PW arithmetic, collective sequence, and memory footprint remain bit-identical (following `dft-dispersion` ADR 0048 and `potential_symmetry` ADR 0049).
5. **Independent verification and oracle isolation**: All production extensions preserve `VerificationReport`, leave frozen `oracle/` default arithmetic untouched, and keep `cargo check/clippy -p dft-engine --no-default-features` green.

---

## 4. Reconciliation of Pre-Stage-R Coupling Design Documents

Prior to Stage R, forces, stress, spin, FD smearing, symmetry, NLCC, D3, and `refscf`, BAANDI authored `docs/architecture/scientific-coupling.md`, `docs/research/coupling-physics.md`, `docs/research/coupling-codex-review.md`, **ADR 0017**, and **ADR 0027**. Reconciling them against the current codebase:

| Design Commitment (Nov 2026) | Status in Current Codebase (2026-11-28) | Action Needed for Extensible-DFT |
| --- | --- | --- |
| Typed physical roles (`ElectrostaticPotentialHaPerE` $\phi$ vs `ElectronEnergyHa` $U_{\mathrm{el}} = -\phi$, `SignedChargeDensityEPerBohr3` vs positive $n$) (ADR 0017) | **Unimplemented** (`dft-units` lacks electric potential and charge types) | Implement as shared foundation for Tracks B and C |
| Local potential seam in `ScfProblem` (`scf.rs`) | **Superseded by ADR 0047**: production assembly lives in `refscf/{map.rs, k_map.rs, electronic_map.rs}` on rank-owned FFT planes | Anchor all new seams in `refscf/` + `dft-physics`, never in frozen `oracle/` |
| Typed `Availability` for unsupported forces/stress (no zero placeholders) | **Resolved for pure PW**: analytic forces & stress exist; `Availability` pattern established | Enforce strict `Unsupported`/`Unavailable` when an extension lacks verified derivatives (e.g., cavity stress or QM/MM box stress) |
| Zero-mode & charged-cell policy (never silently zero $G=0$ of arbitrary $U_{\mathrm{ext}}$; explicit background) | **Still enforced**: `refstart.rs` hard-codes $N = \sum Z_I$; `hartree.rs` hard-zeros $G=0$ | Introduce explicit `NetChargePolicy` / `PotentialReference` before any charged-cell or grand-canonical solve |
| Workflow API vs Low-level Coupling API split (ADR 0027) | **Partially realized**: `refscf::Stepper` + `Runtime` + `baandi-py` + `baandi-mcp` exist, without field injection or electrostatic grid exchange | Extend `refscf::api` and `baandi-py` cleanly along ADR 0027 boundaries |

---

## 5. Concrete Extension Points and Bottlenecks by Track

### 5.1 Track A — Mixed Basis / Plane Wave + Atomic Orbitals / Adaptive Local Basis

- **Reusable assets**: `dft-upf/src/radial.rs` (`RadialTransformTable`, `PSWFC`), `refstart.rs` (atomic-orbital PW projection), `dft-response/src/{subspace.rs, wannier.rs}` (Löwdin orthonormalization and localized subspace transforms).
- **Bottlenecks & high-intrusion points**:
  1. **Generalized eigenproblem $H\psi = \epsilon S\psi$**: `ref_davidson` and `ref_davidson_complex` assume $S = I$; PW+AO bases have non-trivial overlap blocks $S$ with severe ill-conditioning as $E_{\mathrm{cut}}$ rises (requiring canonical orthogonalization or $(1-\hat{P}_{\mathrm{PW}})\phi_\mu$ projection).
  2. **Cross-density & quadrature**: Assembling $n(\mathbf{r}) = n_{\mathrm{PW-PW}} + n_{\mathrm{PW-AO}} + n_{\mathrm{AO-AO}}$ and matrix elements $\langle \chi_\mu | V_{\mathrm{loc}} | \chi_\nu \rangle$, $\langle \mathbf{G} | V_{\mathrm{loc}} | \chi_\nu \rangle$ without FFT aliasing requires dual-grid/multigrid real-space collocation (GPW-style) or DG element quadrature (ALB-style), neither present in `dft-basis`/`dft-physics`.
  3. **Pulay forces & stress + BSSE**: Atom-centered orbitals add $\partial S/\partial \mathbf{R}_I$ and $\partial H/\partial \mathbf{R}_I$ Pulay terms to `forces.rs`/`stress.rs` and introduce Basis Set Superposition Error.
  4. **MPI layout**: `DistributedRealHamiltonian` and `DistributedComplexHamiltonian` partition G-sticks across ranks; hybrid G-stick + atom-centered partitioning requires an isolated operator/solver backend to avoid regressing pure-PW paths.

### 5.2 Track B — Implicit Solvation & Electrochemical Environment (Grand Canonical DFT)

- **Reusable assets**:
  - **Distributed FFT & density gradient**: `StartFields` (`refstart.rs:217–220`) wraps `DistributedPairFft` and already evaluates distributed FFTs and $\nabla n(\mathbf{r})$ for PBE in `refscf/{map.rs, electronic_map.rs}`.
  - **Local potential assembly seam**: `refscf/{map.rs, k_map.rs, electronic_map.rs}` assemble `ionic + vh + vxc` per SCF step; injecting $V_{\mathrm{solv}}[n]$ touches **zero lines** inside `DistributedRealHamiltonian`, `DistributedComplexHamiltonian`, `ref_davidson`, or `KeptHx`.
  - **Fermi–Dirac thermodynamics**: `dft-physics/src/occupations.rs` (`FermiDirac`, `fraction_entropy`) and `refscf/electronic.rs` (`dft_engine::model::ThermodynamicLedger`, `BandGrowth`) already implement verified finite-temperature occupations and $F = E - TS$.
  - **Forces**: For valence-density-only cavities ($\epsilon[n_{\mathrm{val}}]$), $\delta E_{\mathrm{solv}}/\delta n_{\mathrm{val}}$ enters $V_{\mathrm{eff}}$ directly, avoiding cavity Pulay terms (ionic Gaussian pseudo-charge electrostatic forces are still required). **NLCC qualification**: if $\epsilon$ instead depends on $n_{\mathrm{val}} + n_{\mathrm{core}}(\mathbf{R}_I)$, an additional explicit core-cavity force term $\int V_{\mathrm{solv}}(\mathbf{r}) \frac{\partial n_{\mathrm{core}}}{\partial \mathbf{R}_I}\,d\mathbf{r}$ arises alongside the existing NLCC XC core force in `forces.rs::local_or_core_cached`.
- **Extension points & bottlenecks**:
  1. **Generalized Poisson / LPB solver in `dft-physics`**: Solve $\nabla\cdot[\epsilon(\mathbf{r})\nabla\phi] - \kappa^2(\mathbf{r})\phi = -4\pi\rho_{\mathrm{tot}}$ on the existing distributed `StartFields` / `DistributedPairFft` grid, plus dielectric functional $\epsilon[n]$, dielectric derivative $-\frac{1}{8\pi}\frac{d\epsilon}{dn}|\nabla\phi|^2$, and non-electrostatic cavitation/surface-tension terms.
  2. **Ionic pseudo-charge reconciliation**: Smooth Gaussian ionic pseudo-charges $\rho_{\mathrm{ion}}^{\mathrm{sm}}(\mathbf{r})$ for the generalized Poisson RHS must be reconciled cleanly with BAANDI's existing local-pseudopotential + Ewald decomposition without double counting.
  3. **Neutrality relaxation & Grand Canonical mode**: Relax `refstart.rs:509, 2802` under explicit compensating-charge policy, add `FermiDirac::evaluate_at_mu` in `occupations.rs`, and track $\Omega = F - \mu N$ in `ThermodynamicLedger` and `EnergyAudit`.
  4. **Electrostatic potential reference ($\phi_{\mathrm{ref}}$)**: Define bulk-electrolyte ($\kappa > 0$) or vacuum-slab alignment to map $\mu$ to physical electrode potential $U_{\mathrm{SHE}}$.

### 5.3 Track C — Multiscale Coupling (QM/MM and QM/MLIP)

- **Reusable assets**: `refscf::api::{solve, stepper}`, `dft-cli`, `baandi-py`, `python/baandi/ase.py`, post-SCF additive pattern (`dft-dispersion`, ADR 0048), and the **exact same typed electrostatic potential seam as Track B**.
- **Extension points & bottlenecks**:
  1. **`baandi-py` / `python/baandi/ase.py` coverage gap**: Currently restricted to $\Gamma$-only, neutral, spin-unpolarized, 3D-periodic calculations (`ase.py:52–58`). Must expose $k$-points, spin, FD smearing, and typed field exchange.
  2. **Electrostatic embedding core seam**: Requires validated `ExternalPotentialGrid` ($U_{\mathrm{el}}(\mathbf{r})$ in Ha), regularized Gaussian-smeared MM point charges (GEEP-style) to prevent electron spill-out, QM reaction forces on MM charges $\mathbf{F}_a^{\mathrm{MM}} = q_a \mathbf{E}_{\mathrm{QM}}(\mathbf{R}_a)$, QM ion forces $+Z_I \mathbf{E}_{\mathrm{MM}}(\mathbf{R}_I)$, and explicit periodic vs isolated-box Coulomb policy.
  3. **QM/MLIP layering**: Subtractive ONIOM / $\Delta$-ML requires zero Rust SCF changes (adapter/workflow layer with SevenNet/MACE), while electrostatic/polarizable QM/MLIP reuses the Track B/C typed electrostatic seam.

---

## 6. Summary of Architectural Readiness

| Capability Dimension | Current Readiness | Primary Missing Building Block |
| --- | --- | --- |
| Pure-PW SCF (LDA/PBE, NLCC, $k$, Spin, FD, Symmetry, Forces, Stress, D3, MPI+MKL) | **Production-ready** (ADR 0047) | Preserve bit-identical off-path & hot-path speed |
| Typed electrostatic / external field seam (`U_ext`, $\phi$, $\rho_{\mathrm{signed}}$, `EnergyAudit` extension) | **Designed (ADR 0017), unimplemented** | Shared `dft-units` + `dft-physics` + `refscf` local-potential hook (unlocks Tracks B & C) |
| Generalized Poisson / LPB & implicit dielectric cavity ($\epsilon[n]$, smooth ionic charge) | **Not implemented** (`StartFields`/`DistributedPairFft` & $\nabla n$ ready) | Isolated `dft-physics` generalized Poisson solver + `refscf` map hook |
| Charged cells & Grand Canonical DFT (fixed $\mu$, $\Omega = F - \mu N$) | **Partially ready** (`FermiDirac` & `ThermodynamicLedger` exist) | Relax neutrality in `refstart.rs`, add `evaluate_at_mu`, add compensating electrolyte/background |
| Non-orthogonal mixed PW+AO / ALB basis (Track A) | **Low readiness** (only UPF `PSWFC` startup & R4a/R4b subspace exist) | Generalized $H\psi = \epsilon S\psi$ solver, cross-density quadrature, Pulay forces/stress |
| External QM/MLIP & QM/MM driver layer (Track C) | **Partial** (`baandi-py`/`ase.py` exists, Gamma/neutral/unpolarized-only) | Broaden `baandi-py`/`ase.py` + add typed `U_ext`/MM-force interface |
