# Extensible BAANDI Research Program — Gate 0 Triage Report

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Tracks A, B, C)
- **Phase**: Gate 0 — Repository Assessment & Triage (analysis only; no production code, no large refactoring, no premature full-track theory specs)
- **Date**: 2026-11-28
- **Assessed Baseline**: `main` at `24d48b65` (includes `RPA-R4A-ATOMIC-PROJECTION` and `RPA-R4B-WANNIER90`, verified via `dev08` gate `a9e50e58`)
- **Companion Documents**: [repository-assessment.md](repository-assessment.md), [decision-log.md](decision-log.md)

---

## 1. Grounded Answers to the 10 Mandatory Gate 0 Questions (Directive Section 12)

### Q1. Under what conditions is Mixed-basis DFT substantially advantageous over pure PW-DFT?

From Louie, Ho, and Cohen (1979, DOI `10.1103/PhysRevB.19.1774`), Lippert, Hutter, and Parrinello (1997, DOI `10.1080/002689797170220`), and Lin et al. (2012, DOI `10.1016/j.jcp.2011.11.032`) evaluated against BAANDI's norm-conserving PW cost structure:

- **Physical regimes where mixed PW+AO can win**:
  1. **Severe species cutoff mismatch in hybrid systems** (e.g., $\text{CO}/\text{Al}(111)$ or $\text{H}_2\text{O}/\text{MgO}(001)$): A soft metal substrate (Al, Na) converges at $E_{\mathrm{cut}} \approx 30\text{–}40$ Ry, whereas a localized $2p$ adsorbate (C, O) or $3d$ transition-metal active site (Cu, Pt, Fe) forces the *entire* supercell PW sphere to $80\text{–}120$ Ry ($3\text{–}5\times$ more PWs everywhere, including vacuum). Localizing the high-momentum curvature on a handful of atom-centered orbitals lets the global PW cutoff stay at the soft substrate value.
  2. **Large vacuum fraction / low-dimensional systems** (slabs, isolated molecules, nanoparticles, electrocatalytic interfaces): Pure PW basis size scales with total cell volume $\Omega E_{\mathrm{cut}}^{3/2}$ regardless of empty vacuum, whereas atom-centered enrichment scales with the number of "hard" atoms $N_{\mathrm{hard}}$.
- **Why PW count reduction alone is NOT sufficient (full Section 3.5 cost model)**:
  - Even if $N_{\mathrm{PW}}$ drops by $3\text{–}5\times$, a localized AO with sharp radial structure still introduces high-$G$ Fourier components into $n_{\mathrm{AO-AO}}(\mathbf{r})$ and the cross density $n_{\mathrm{PW-AO}}(\mathbf{r})$. Unless the **real-space integration grid is decoupled from the wavefunction PW cutoff** (dual-grid / multigrid collocation as in GPW, or separate radial quadrature for AO–AO blocks), the density FFT and XC/Hartree cost remain high.
  - Furthermore, mixed-basis introduces four costs absent from BAANDI's pure-PW path: (i) generalized eigenproblem $H\psi = \epsilon S\psi$ and $O(N_{\mathrm{AO}}^3)$ / rank-revealing orthogonalization to control ill-conditioning of $S$, (ii) PW–AO transformation and overlap/Hamiltonian block assembly, (iii) potential degradation in Davidson preconditioning and breakdown of `KeptHx` warm starts (`refscf/{kept_hx.rs, k_kept.rs}`), and (iv) Pulay force/stress derivatives ($\partial S/\partial \mathbf{R}_I$, $\partial H/\partial \mathbf{R}_I$).
  - **Net assessment**: On homogeneous dense solids with smooth PseudoDojo NC pseudopotentials (e.g., Si bulk), BAANDI's highly optimized `DistributedRealHamiltonian` / `DistributedComplexHamiltonian` will almost certainly beat mixed-basis. Mixed-basis becomes substantially advantageous only in **vacuum-rich or strongly cutoff-mismatched adsorbate/surface/3d-oxide systems** *and only if* the cross-density quadrature and Pulay-force overhead stay well below the $O(N_b^2 N_{\mathrm{PW}})$ subspace-algebra savings. *(Quantifying the exact crossover is a working hypothesis deferred to Track A Gate 1–2.)*

### Q2. What research value does adaptive local basis enrichment have compared to existing mixed-basis methods?

Comparing Directive Section 3.4 Approaches A–E (Nonorthogonal Mixed Basis, Orthogonalized Enrichment $(1-\hat{P}_{\mathrm{PW}})\phi_\mu$, Optimized Local Basis, Adaptive Basis Enrichment, Independent AO Backend) and Section 3.6 variational properties:

- **What classical mixed-basis lacks**: Fixed atom-centered AO recipes (Louie et al. 1979) require hand-tuned orbital sets per species/oxidation state, suffer from **near-linear dependence** against PWs as $E_{\mathrm{cut}}$ increases, and produce uncontrollable **Basis Set Superposition Error (BSSE)** when geometries or adsorption distances vary.
- **Where adaptive local basis enrichment adds genuine research value**:
  1. **Variational error-targeted enrichment**: For a fixed Hermitian $H$, nested spaces $\mathcal{V}_{\mathrm{PW}} \subseteq \mathcal{V}_{\mathrm{PW+AO}}$ guarantee $\epsilon_i^{\mathrm{PW+AO}} \le \epsilon_i^{\mathrm{PW}}$ by Rayleigh–Ritz. Unlike static AOs, an adaptive scheme adds orthogonalized complement functions $(1-\hat{P}_{\mathrm{PW}})\phi_\mu$ (or domain-adapted functions as in DG-ALB, Lin et al. 2012) *only* on atoms/regions where a local error indicator exceeds tolerance—directly targeting **observable accuracy** (adsorption energy, active-site forces) rather than global total energy.
  2. **Critical mathematical subtlety (outside-subspace residual)**: Evaluating the Kohn–Sham residual $R_n = (\hat{H} - \epsilon_n)\psi_n$ *inside* the current truncated PW+AO subspace vanishes at convergence and **cannot measure basis incompleteness**. A valid a posteriori error estimator must project $R_n$ onto a larger test/enrichment space (e.g., a higher-$E_{\mathrm{cut}}$ shell or candidate local-orbital pool $\{(1-\hat{P}_{\mathrm{PW}})\phi_\alpha\}$) via perturbation estimates $\Delta E_\alpha \approx -\sum_n f_n |\langle \tilde{\phi}_\alpha | R_n \rangle|^2 / (\langle \tilde{\phi}_\alpha | \hat{H} | \tilde{\phi}_\alpha \rangle - \epsilon_n)$—a genuinely non-trivial research problem closely aligned with BAANDI's existing `dft-response/src/subspace.rs` projection infrastructure.
  3. **Key unresolved risks**: Discrete enrichment steps cause **energy and force discontinuities** along relaxation/MD trajectories unless the enrichment basis is frozen per geometry segment (analogous to BAANDI's frozen-Miller variable-cell segments in ADR 0047/P23A) or smoothly weighted, and BSSE still requires explicit counterpoise or environment-independent local-space control.

### Q3. What changes does implicit solvation require in BAANDI's current computational structure?

Anchored directly in BAANDI's post-Stage-R codebase ([repository-assessment.md](repository-assessment.md)):

1. **Zero changes to pure-PW Hamiltonian apply and eigensolvers**: `DistributedRealHamiltonian::apply_block`, `DistributedComplexHamiltonian::apply_block`, `ref_davidson`, `ref_davidson_complex`, and `KeptHx` consume only the rank-owned real-space local potential `op_potential` via `set_owned_potential(&op_potential, epoch)`. Adding a scalar solvation potential $V_{\mathrm{solv}}[n](\mathbf{r})$ preserves all hot-path invariants verbatim.
2. **New isolated generalized Poisson / LPB kernel in `dft-physics`**:
   - Solves $\nabla\cdot[\epsilon(\mathbf{r})\nabla\phi(\mathbf{r})] - \kappa^2(\mathbf{r})\phi(\mathbf{r}) = -4\pi\rho_{\mathrm{tot}}(\mathbf{r})$ on the existing distributed density FFT grid (`StartFields` in `refstart.rs:217–220`, wrapping `DistributedPairFft` and implementing `dft_physics::gga::distributed::Fields<K>`), using FFT-preconditioned CG (Fisicaro et al. 2016, DOI `10.1063/1.4939125`) for SCCS (Andreussi, Dabo, Marzari 2012, DOI `10.1063/1.3676407`) or VASPsol/LPB (Mathew et al. 2014, DOI `10.1063/1.4865107`).
   - Evaluates smooth dielectric switching $\epsilon[n]$, dielectric derivative potential $V_\epsilon = -\frac{1}{8\pi}\frac{d\epsilon}{dn}|\nabla\phi|^2$, and non-electrostatic cavitation/surface-tension terms.
3. **Ionic pseudo-charge & local-potential seam**:
   - `refscf/{map.rs, k_map.rs, electronic_map.rs}` assemble `ionic + vh + vxc` per SCF step; an optional `SolvationState` adds $V_{\mathrm{solv}}[n]$ (or replaces vacuum $V_H$ with the polarization difference $\Delta\phi_{\mathrm{pol}} = \phi_{\mathrm{solv}} - \phi_{\mathrm{vac}}$) while reconciling smooth Gaussian ionic pseudo-charges $\rho_{\mathrm{ion}}^{\mathrm{sm}}$ with `state.ionic` and `ewald.rs` without double counting.
4. **Energy, force, and stress ledgers**:
   - Extend `crates/dft-engine/src/energy.rs::EnergyAudit` with explicit solvation electrostatic + non-electrostatic energy and band double-counting subtraction ($-\int n_{\mathrm{out}} V_{\mathrm{solv}}[n_{\mathrm{in}}]\,d\mathbf{r} + E_{\mathrm{solv}}[n_{\mathrm{out}}]$).
   - Forces (`crates/dft-physics/src/forces.rs`, `refscf::api::Observables::{evaluate_forces, evaluate_forces_electronic}`): add ionic Gaussian pseudo-charge electrostatic forces in $\phi_{\mathrm{pol}}$. For valence-density-only cavities $\epsilon[n_{\mathrm{val}}]$, no cavity Pulay term arises; **if NLCC is enabled and $\epsilon$ depends on $n_{\mathrm{val}} + n_{\mathrm{core}}(\mathbf{R}_I)$**, an explicit core-cavity force term $\int V_{\mathrm{solv}}(\mathbf{r})\frac{\partial n_{\mathrm{core}}}{\partial \mathbf{R}_I}\,d\mathbf{r}$ must also be added (reusing the `local_or_core_cached` structure-factor derivative). Stress initially returns typed `Availability::Unsupported` until strain derivatives of $\epsilon[n]$ and cavitation are independently verified.
5. **Grand Canonical extension (Phase B2, Sundararaman, Goddard, Arias 2017, DOI `10.1063/1.4978411`)**:
   - Relax hard-coded $N = \sum Z_I$ in `refstart.rs:509–515` (SCF) and `refstart.rs:2746–2813` (NSCF) ONLY under an explicit compensating-charge / LPB electrolyte policy.
   - Add `FermiDirac::evaluate_at_mu` in `crates/dft-physics/src/occupations.rs` (reusing `fraction_entropy` at lines 59–66), track grand potential $\Omega = F - \mu N$ in `dft_engine::model::ThermodynamicLedger` (`crates/dft-engine/src/model/spin.rs`), and fix the electrostatic potential reference $\phi_{\mathrm{ref}}$ via bulk electrolyte ($\kappa > 0$) or slab vacuum alignment (Otani & Sugino 2006, DOI `10.1103/PhysRevB.73.115407`).

### Q4. What are the minimal conditions for physically consistent QM/MM and QM/MLIP coupling?

Grounded in Senn & Thiel (2009, DOI `10.1002/anie.200802019`), Laino et al. (2005, DOI `10.1021/ct050123f`), Chung et al. (2015, DOI `10.1021/cr5004419`), Batatia et al. (2022, MACE, arXiv:2206.07697), Park et al. (2024, SevenNet, DOI `10.1021/acs.jctc.4c00190`), ADR 0017, and Directive Sections 5.3–5.5:

1. **Exact energy decomposition, reference alignment, and zero double-counting**:
   - In subtractive coupling ($E_{\mathrm{total}} = E_{\mathrm{low}}(\text{Full}) + E_{\mathrm{QM}}(\text{QM}) - E_{\mathrm{low}}(\text{QM})$), every energy term must carry explicit ownership; constant per-species atomic energy offsets cancel in $E_{\mathrm{low}}(\text{Full}) - E_{\mathrm{low}}(\text{QM})$ ONLY when the QM region composition/termination matches exactly (link atoms require consistent capping subtraction).
   - In electrostatic embedding, adding $U_{\mathrm{ext}}(\mathbf{r}) = -\phi_{\mathrm{MM}}(\mathbf{r})$ to $H_{\mathrm{KS}}$ is **not** a complete energy: one must consistently include $E_{\mathrm{el-MM}} = \int n(\mathbf{r}) U_{\mathrm{ext}}(\mathbf{r})\,d\mathbf{r}$ (without a $1/2$ factor), $E_{\mathrm{ion-MM}} = \sum_I Z_I \phi_{\mathrm{MM}}(\mathbf{R}_I)$, MM–MM self-energy (owned by the MM driver, never duplicated), and zero double-counting in `EnergyAudit::band_form`.
2. **Sign, unit, and zero-mode gauge rigor (ADR 0017)**: Enforce typed distinction between $\phi$ (Ha/e) and $U_{\mathrm{el}} = -\phi$ (Ha), preserve non-zero mean $U_0$ explicitly (neutral QM total shift $(\sum Z_I - N)\phi_0 = 0$ requires both electron and ion terms), and forbid silently applying Hartree's $G=0 \mapsto 0$ rule to external potentials.
3. **Regularized short-range electrostatics (preventing electron spill-out)**: Bare $1/r$ attractive classical point charges on a periodic PW grid cause unphysical electron density collapse ("spill-out") and overpolarization. Regularized Gaussian-smeared charges (GEEP, Laino et al. 2005) or repulsive core pseudization are mandatory for physical electrostatic embedding.
4. **Covalent boundary / link-atom consistency**: When QM/MM or QM/MLIP cuts a covalent bond, link-atom placement must be a deterministic function of host QM/MM positions so that chain-rule force projection onto real atoms preserves exact $\mathbf{F} = -\nabla E_{\mathrm{total}}$ and Newton's third law.
5. **Variational force consistency & MD energy conservation**:
   - QM ions must receive $+Z_I \mathbf{E}_{\mathrm{MM}}(\mathbf{R}_I)$ (and any grid-moving/cavity derivatives), while MM charges must receive exact reaction forces $\mathbf{F}_a^{\mathrm{MM}} = q_a \mathbf{E}_{\mathrm{QM}}(\mathbf{R}_a)$ evaluated on the converged QM density.
   - Tight inner SCF convergence (verified via `VerificationReport`) is mandatory to prevent systematic energy drift in NVE MD.
6. **Boundary-condition consistency & MLIP validity bounds**:
   - Periodic 3D PW cells cannot silently masquerade as isolated QM clusters; periodic-image interactions vs vacuum-padded / truncated Coulomb boundaries must be declared explicitly.
   - MLIPs (SevenNet, MACE) are typically short-/medium-range message-passing models trained on neutral bulk/slab DFT; coupling must never assume an MLIP captures long-range electrostatics, charge transfer, or explicit polarization unless specifically parameterized and validated for that regime, and must log model provenance/version/uncertainty alongside BAANDI metadata.

### Q5. What must be commonized vs. kept independent when combining different physical models?

Mapping Directive Section 6.2 (the 8 independent modeling axes) and Section 6.6 (basis-independent physics vs backend-specific numerics) onto BAANDI's crate structure:

| Layer | **Commonize (Shared Contracts)** | **Keep Independent (Specialized Backends)** |
| --- | --- | --- |
| **Units, Signs & Physical Roles** | `dft-units` + ADR 0017 semantic newtypes ($\phi$, $U_{\mathrm{el}}$, $n$, $\rho_{\mathrm{signed}}$, $\mu$, $\Omega$) | Internal raw `f64` / `Complex64` slice loops inside hot kernels |
| **Modeling Axes & Compatibility (Sec. 6.2)** | Explicit typed configuration & preflight compatibility matrix across the 8 axes (**Basis**, **Boundary**, **Environment**, **Embedding**, **Ensemble**, **Solver**, **Execution**, **Accuracy Control**); fail-closed rejection of unsupported combinations before SCF starts | Each supported combination dispatches statically/via enum to its own tuned driver (no dynamic dispatch inside inner loops) |
| **Thermodynamic & Energy Ledger (Sec. 6.3)** | Unified `EnergyAudit` + `ThermodynamicLedger` for generalized functional $\mathcal{F}[\{\psi_n\},\{f_n\},\rho,\{\mathbf{R}_I\},\boldsymbol{\xi}]$ with explicit component ownership, double-counting checks, and typed `Availability` for forces/stress | Basis-specific kinetic/nonlocal/Pulay evaluations and environment-specific internal auxiliary solvers (e.g., PCG for $\phi_{\mathrm{solv}}$, Ewald, external MLIP) |
| **Continuous Grid Physics (Sec. 6.6)** | Real-space XC functionals, density-dependent dielectric $\epsilon[n]$, smooth charge deposition, Fermi–Dirac occupation math (`occupations.rs`) | Discretization & operator application: half-sphere real Gamma FFT (`DistributedRealHamiltonian`), full-sphere complex Bloch (`DistributedComplexHamiltonian`), or any future AO/mixed-basis integral & $H\psi = \epsilon S\psi$ solver |
| **Runtime, Provenance & Agent Interface (Sec. 6.5)** | `refscf::Stepper` + `Runtime` + structured diagnostics, provenance hashes, verification reports, and `baandi-py`/`baandi-mcp` adapters | External MM/MLIP engines (SevenNet, MACE, OpenMM) remain decoupled external processes/calculators, never linked into `dft-engine` |

### Q6. Can BAANDI be extended while minimizing impact on existing high-performance PW paths?

**Yes—with near-zero impact for Tracks B and C, and only via strict backend isolation for Track A.**

- **Verified pure-PW hot paths**: `DistributedRealHamiltonian::apply_block`, `DistributedComplexHamiltonian::apply_block`, `ref_davidson`, `ref_davidson_complex`, `KeptHx`, and `forces_and_frozen_stress`.
- **Track B (Implicit Solvation & GC-DFT)**: Hooks strictly into the outer per-SCF-step potential assembly (`refscf/{map,k_map,electronic_map}`), `occupations.rs`, and post-SCF `forces.rs` behind `Option<SolvationConfig>` / `Option<GrandCanonicalConfig>`. When `None`, zero extra code executes in hot loops, preserving bit-identical pure-PW arithmetic, zero extra MPI collectives in `apply_block`, and unchanged memory envelope. Even when *enabled*, generalized Poisson runs once per SCF iteration (~1–5% of total SCF time on medium/large cells once warm-started from the previous SCF step's $\phi$!).
- **Track C (QM/MM & QM/MLIP)**: Mechanical/subtractive QM/MLIP touches **zero** Rust SCF code (only `baandi-py`/`ase.py` coverage expansion). Electrostatic QM/MM injects a precomputed (or outer-updated) real-space `Option<&[f64]>` grid at the exact same `refscf` local-potential seam—again zero `apply_block`/Davidson intrusion.
- **Track A (Mixed Basis)**: Cannot share `DistributedRealHamiltonian` or `ref_davidson` ($S = I$) without polluting hot paths. Minimizing impact requires keeping pure-PW drivers untouched and isolating any Track A solver in a separate module/crate.

### Q7. What are the expected implementation difficulty (code intrusion) and verification difficulty of each research track?

| Criterion | Track A: Mixed Basis (PW+AO / ALB) | Track B: Implicit Solvation & GC-DFT | Track C: QM/MM & QM/MLIP |
| --- | --- | --- | --- |
| **Code Intrusion** | **Very High** (new generalized $H\psi = \epsilon S\psi$ eigensolver, overlap/cross-density discretization, Pulay forces/stress, new MPI row partition) | **Low–Moderate** (isolated generalized Poisson/LPB solver in `dft-physics` using existing `StartFields`/`DistributedPairFft` + `refscf` local-potential/ledger hook; `occupations.rs` fixed-$\mu$ for Phase B2) | **Low** (Layer 1 mechanical QM/MLIP = zero Rust SCF intrusion; Layer 2 electrostatic QM/MM reuses Track B's typed `U_ext` & signed-charge seam) |
| **Pure-PW Perf Risk** | **High** (unless completely duplicated as a separate backend) | **Very Low** (`Option::None` bit-identical off-path; zero `apply_block`/Davidson edits) | **Very Low** (external driver + `Option::None` off-path) |
| **Verification Difficulty (Levels 1–5)** | **Hard** (Pulay forces, BSSE separation, ill-conditioning of $S$, cross-code comparison requires matching AO basis definitions) | **Excellent / High-Confidence** (L1: exact $\epsilon=1$ vacuum limit, MMS PDE residual $< 10^{-12}$, analytic Born & Debye–Hückel limits; L2: central-FD force & $\partial\Omega/\partial\mu = -N$ checks; L3–L4: published SCCS/VASPsol/JDFTx hydration energies & slab capacitance) | **Moderate** (L1–L2: exact $U_0$ gauge shift, subtractive energy/force conservation, QM/MM FD forces; L4: electrostatic spill-out, link-atom & MLIP validity bounds need careful physical benchmarking) |
| **Accessible Oracles (`dev08` & Local)** | Pure-PW high-$E_{\mathrm{cut}}$ BAANDI & QE 7.6 (`/scr/mino/oracles/qe-7.6/bin/pw.x`) | Internal MMS/Born/FD/vacuum oracles + published SCCS/VASPsol/JDFTx tables *(QE+Environ / JDFTx installability on `dev08` to check in Gate 1)* | Internal $U_0$/FD/subtractive oracles + Python/ASE *(SevenNet/MACE availability on `dev08` to check in Gate 1)* |

### Q8. Is there substantial differentiation and research originality compared to existing DFT codes?

Comparing honestly against mature incumbents (QE+Environ, JDFTx, VASPsol, CP2K, GPAW, SIESTA/OpenMX):

- **Where BAANDI would merely duplicate incumbents (low differentiation)**:
  - Re-implementing a static GPW engine (competing head-on with 25 years of CP2K) or standard LCAO (SIESTA/OpenMX) would offer little research originality.
  - Implementing a bare-bones VASPsol clone without electrochemical/multiscale composability or verification rigor would add utility but limited scientific novelty.
- **Where BAANDI has genuine, high-leverage differentiation**:
  1. **Memory-safe, MPI+MKL high-performance Rust PW+Electrochemistry core**: Existing electrochemical PW-DFT stacks are Fortran/C++ monoliths (QE+Environ, VASPsol, JDFTx) with fragile global state, ad-hoc unit/sign conventions, and brittle restart/control loops. BAANDI combines QE-class speed, `#![deny(unsafe_code)]` physics kernels, **compile-time typed electrostatic semantics (ADR 0017)**, and **verified fixed-point SCF certificates (`VerificationReport`)**.
  2. **Stepwise, agent-native multiscale & electrochemical control (`refscf::Stepper` + `Runtime` + MCP)**: No existing PW-DFT code exposes a deterministic, checkpoint-v2-backed, step-by-step runtime where an AI agent or outer multiscale/electrochemical controller can inspect residuals, diagnose dielectric/GC-SCF charge sloshing, adjust mixing/potential schedules safely, and compose QM + Implicit Electrolyte + External MLIP (SevenNet/MACE) with full machine-readable provenance.
  3. **Composable path to the Target Demonstration (Section 8 — Electrochemical Interface)**: Combining verified implicit solvation + LPB Grand Canonical DFT + typed electrostatic embedding (allowing explicit MLIP/MM solvent shells inside a continuum bulk electrolyte without double counting) addresses one of the hardest open engineering/scientific problems in computational catalysis and batteries.
  4. *(Longer-term Track A originality)*: Observable-targeted **adaptive orthogonalized local enrichment** driven by outside-subspace residual estimators is a genuinely novel research direction, best explored *after* the shared electrostatic/multiscale seams are proven.

### Q9. Which of the three research tracks is most rational as the first vertical breakthrough (Gate 1 → Gate 2)?

**Recommendation: Track B — Phase B1 (Typed Electrostatic Seam + Density-Dependent Implicit Solvation via Generalized Poisson / Linearized Poisson–Boltzmann), designed with explicit hooks for Phase B2 (Grand Canonical DFT) and Track C (Electrostatic QM/MM & QM/MLIP).**

Rationale:
1. **Maximum cross-track leverage per line of core code**: Phase B1 builds the exact ADR 0017 typed electrostatic foundation (`dft-units` / `dft-physics`), distributed generalized Poisson/LPB solver (reusing `StartFields` / `DistributedPairFft`), and `refscf` local-potential + `EnergyAudit` seam that **unlocks both Track B Phase B2 (GC-DFT, where LPB $\kappa > 0$ fixes the $G=0$ potential reference)** AND **Track C electrostatic QM/MM & polarizable QM/MLIP**.
2. **Lowest regression risk to BAANDI's pure-PW asset**: Zero edits to `DistributedRealHamiltonian`, `DistributedComplexHamiltonian`, `ref_davidson`, or `KeptHx`.
3. **Strongest Gate 2 numerical verification**: Exact vacuum reduction ($\epsilon_b=1, \kappa=0, \tau=0$), Level-1 MMS PDE verification ($< 10^{-12}$), analytic Born/Debye–Hückel limits, Level-2 force finite-difference tests, and MPI rank-independence tests—all executable directly on `dev08` without external-code lock-in.

*(Alternative options preserved in [decision-log.md](decision-log.md): Option 2 starts with Track C `baandi-py`/`ase.py` + mechanical QM/MLIP + fixed $U_{\mathrm{ext}}$; Option 3 starts with an isolated Track A mixed-basis prototype.)*

### Q10. Are there more promising or fundamentally different approaches beyond those proposed?

Yes—four genuinely distinct paradigms merit awareness and selective evaluation during Gate 1–2 rather than assuming the classical formulations are optimal:

1. **Direct Free-Energy / Grand-Potential Minimization instead of Broyden Mixing for Solvation & GC-DFT** (Freysoldt, Boeck, Neugebauer 2009, DOI `10.1103/PhysRevB.79.241103`; Sundararaman, Goddard, Arias 2017, DOI `10.1063/1.4978411`):
   - Standard charge-density Broyden mixing frequently oscillates or stalls in implicit-solvent charged slabs and constant-potential GC-DFT because the dielectric derivative $\frac{d\epsilon}{dn}|\nabla\phi|^2$ and $G=0$ charge exchange amplify low-$q$ charge sloshing. Preconditioned **direct minimization** of the coupled free energy $\mathcal{F}[\{\psi_n\}, \{f_n\}, \phi]$ or grand potential $\Omega(\mu)$ guarantees monotonic descent and naturally couples orbital, occupation, and dielectric degrees of freedom. *Action*: Evaluate Broyden stability vs. inner-outer / direct-minimization alternatives explicitly in Track B Gate 1–2.
2. **Explicit-MLIP First-Solvation-Shell + Implicit Continuum Bulk Hybrid (Skipping Pure Continuum Parameter Fitting)**:
   - Pure implicit solvation struggles with directional hydrogen bonding and specific ion adsorption at electrochemical interfaces. Because modern universal MLIPs (MACE, SevenNet) evaluate thousands of explicit water molecules in milliseconds, coupling a small QM active site + explicit MLIP water shell + bulk continuum dielectric (via Track B + Track C) may deliver far higher physical accuracy than tuning increasingly complex pure-continuum cavity models—provided explicit–implicit boundary double-counting is rigorously avoided.
3. **Compressed Subspace / Interpolative Separable Density Fitting (ISDF) or Dual-Grid Enrichment instead of Classical Non-Orthogonal PW+AO** (Lu & Ying 2015, DOI `10.1016/j.jcp.2015.09.014`):
   - Rather than replacing BAANDI's orthonormal PW Davidson with an ill-conditioned non-orthogonal PW+AO eigenproblem, real-space **interpolative separable density fitting (ISDF)** or **PAW/ultrasoft dual-grid augmentation** directly addresses the high-$E_{\mathrm{cut}}$ density/exchange bottleneck while retaining well-conditioned matrix-free operators (and also accelerates BAANDI's existing RPA/cRPA pipeline in `dft-response`!).
4. **Agent-Guided Adaptive Numerical Control on Top of `refscf::Stepper`**:
   - Instead of hard-coding brittle heuristics for dielectric-cavity SCF damping, GC-SCF $\mu$-stepping, or adaptive QM-region expansion, BAANDI's existing `Stepper` + `Runtime` + `baandi-mcp` architecture allows structured, verifiable policy control (e.g., adaptive solver tolerance coupling between inner generalized-Poisson PCG and outer KS Davidson, or MLIP uncertainty-triggered QM checks) while keeping physics validation fail-closed inside Rust.

---

## 2. Explicit Working Hypotheses & Open Questions (Deferred to Gate 1)

1. **Ionic pseudo-charge vs. Ewald decomposition**: Whether solving for the *polarization potential difference* $\Delta\phi_{\mathrm{pol}} = \phi_{\mathrm{solv}} - \phi_{\mathrm{vac}}$ (avoiding modification of short-range vacuum Ewald/local-PP terms) or solving directly for total $\phi$ with Gaussian ionic pseudo-charges $\rho_{\mathrm{ion}}^{\mathrm{sm}}$ gives cleaner variational forces and lower cutoff sensitivity in BAANDI.
2. **Dielectric parametrization & NLCC density argument**: SCCS (Andreussi et al. 2012) vs. VASPsol/LPB (Mathew et al. 2014) on PseudoDojo NC pseudopotentials—specifically whether $\epsilon$ should depend on $n_{\mathrm{val}}$ alone or $n_{\mathrm{val}} + n_{\mathrm{core}}(\mathbf{R}_I)$ (which introduces an explicit core-cavity force term).
3. **Generalized Poisson preconditioning & MPI scaling**: Fisicaro et al. (2016) $\tilde{\phi} = \sqrt{\epsilon}\phi$ PCG vs. direct FFT-preconditioned CG iteration counts and collective overhead on `dev08` at $\epsilon_b = 78.4$ (water) for molecule and slab cells.
4. **Broyden SCF stability under nonlinear $V_\epsilon[n]$**: Whether `refscf/broyden.rs` converges reliably with $V_\epsilon = -\frac{1}{8\pi}\frac{d\epsilon}{dn}|\nabla\phi|^2$, or whether inner-outer damping / direct minimization is needed.
5. **External oracle availability on `dev08`**: Checking in Gate 1 whether QE+Environ, JDFTx, or SevenNet/MACE can be run in an isolated environment on `dev08` alongside existing QE 7.6, Yambo 5.3.0, and ABINIT 10.8.3.

---

## 3. Verified Literature References

All 17 references below were verified via primary bibliographic lookup during Gate 0 (zero unverified DOIs):

- **Track A (Mixed Basis / GPW / ALB / Compression)**:
  - S. G. Louie, K.-M. Ho, and M. L. Cohen (1979), *Self-consistent mixed-basis approach to the electronic structure of solids*, *Phys. Rev. B* 19, 1774–1782. DOI: `10.1103/PhysRevB.19.1774`
  - G. Lippert, J. Hutter, and M. Parrinello (1997), *A hybrid Gaussian and plane wave density functional scheme*, *Mol. Phys.* 92(3), 477–488. DOI: `10.1080/002689797170220`
  - L. Lin, J. Lu, L. Ying, and W. E (2012), *Adaptive local basis set for Kohn–Sham density functional theory in a discontinuous Galerkin framework I: Total energy calculation*, *J. Comput. Phys.* 231(4), 2140–2154. DOI: `10.1016/j.jcp.2011.11.032`
  - J. Lu and L. Ying (2015), *Compression of the electron repulsion integral tensor in tensor hypercontraction format with cubic scaling cost*, *J. Comput. Phys.* 302, 329–335. DOI: `10.1016/j.jcp.2015.09.014`
- **Track B (Implicit Solvation, Electrochemistry & Direct Minimization)**:
  - M. Otani and O. Sugino (2006), *First-principles calculations of charged surfaces and interfaces: A plane-wave nonrepeated slab approach*, *Phys. Rev. B* 73, 115407. DOI: `10.1103/PhysRevB.73.115407`
  - C. Freysoldt, S. Boeck, and J. Neugebauer (2009), *Direct minimization technique for metals in density functional theory*, *Phys. Rev. B* 79, 241103. DOI: `10.1103/PhysRevB.79.241103`
  - O. Andreussi, I. Dabo, and N. Marzari (2012), *Revised self-consistent continuum solvation in electronic-structure calculations*, *J. Chem. Phys.* 136, 064102. DOI: `10.1063/1.3676407`
  - K. Mathew, R. Sundararaman, K. Letchworth-Weaver, T. A. Arias, and R. G. Hennig (2014), *Implicit solvation model for density-functional study of nanocrystal surfaces and reaction pathways*, *J. Chem. Phys.* 140, 084106. DOI: `10.1063/1.4865107`
  - G. Fisicaro, L. Genovese, O. Andreussi, N. Marzari, and S. Goedecker (2016), *A generalized Poisson and Poisson-Boltzmann solver for electrostatic environments*, *J. Chem. Phys.* 144, 014103. DOI: `10.1063/1.4939125`
  - P. Giannozzi, O. Andreussi, et al. (2017), *Advanced capabilities for materials modelling with QUANTUM ESPRESSO*, *J. Phys.: Condens. Matter* 29, 465901. DOI: `10.1088/1361-648X/aa8f79`
  - R. Sundararaman, K. Letchworth-Weaver, K. A. Schwarz, D. Gunceler, Y. Ozhabes, and T. A. Arias (2017), *JDFTx: software for joint density-functional theory*, *SoftwareX* 6, 278–284. DOI: `10.1016/j.softx.2017.10.006`
  - R. Sundararaman, W. A. Goddard III, and T. A. Arias (2017), *Grand canonical electronic density-functional theory: Algorithms and applications to electrochemistry*, *J. Chem. Phys.* 146, 114104. DOI: `10.1063/1.4978411`
- **Track C (QM/MM & QM/MLIP)**:
  - T. Laino, F. Mohamed, A. Laio, and M. Parrinello (2005), *An Efficient Real Space Multigrid QM/MM Electrostatic Coupling*, *J. Chem. Theory Comput.* 1(6), 1176–1184. DOI: `10.1021/ct050123f`
  - H. M. Senn and W. Thiel (2009), *QM/MM Methods for Biomolecular Systems*, *Angew. Chem. Int. Ed.* 48(7), 1198–1229. DOI: `10.1002/anie.200802019`
  - L. W. Chung, W. M. C. Sameera, R. Ramozzi, K. Morokuma, et al. (2015), *The ONIOM Method and Its Applications*, *Chem. Rev.* (2015). DOI: `10.1021/cr5004419`
  - I. Batatia, D. P. Kovács, G. N. C. Simm, C. Ortner, and G. Csányi (2022), *MACE: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields*, *NeurIPS* 35 (arXiv:2206.07697, DOI: `10.52202/068431-0830`)
  - Y. Park, J. Kim, S. Hwang, and S. Han (2024), *Scalable Parallel Algorithm for Graph Neural Network Interatomic Potentials in Molecular Dynamics Simulations* (SevenNet), *J. Chem. Theory Comput.* 20. DOI: `10.1021/acs.jctc.4c00190`
