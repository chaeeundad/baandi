# Track B Gate 1: Implicit Solvation Theory & Variational Formulation for BAANDI

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Track B — Phase B1 Implicit Solvation + Phase B2 GC-DFT Readiness)
- **Phase**: Gate 1 — Theoretical Review (complete 12-item specification per Directive Section 11; zero production Rust changes)
- **Date**: 2026-11-28
- **Baseline**: `main` at `939092d8` (includes `RPA-R4A`, `RPA-R4B`, `EXTENSIBLE-DFT-GATE0`)
- **Companion Documents**: [algorithm-comparison.md](algorithm-comparison.md), [repository-assessment.md](repository-assessment.md), [triage-report.md](triage-report.md), [decision-log.md](decision-log.md)

---

## 1. Problem Definition

Standard BAANDI calculations solve periodic Kohn–Sham DFT in vacuum under fixed integer neutrality ($N_e = \sum_I Z_I$). Real catalytic, molecular, and electrochemical systems occur inside a polarizable liquid dielectric and mobile electrolyte ions, often at controlled electrode potential $\mu_e$.

This document establishes a **variationally exact, zero-double-counting, MPI-native continuum solvation theory** tailored to BAANDI's norm-conserving plane-wave conventions (`docs/physics-conventions.md`, ADR 0047), satisfying four strict requirements:
1. **Exact physical consistency**: Total energy $E_{\mathrm{solv}}$, Kohn–Sham effective potential $V_{\mathrm{solv}} = \delta E_{\mathrm{solv}}/\delta n$, and ionic forces $\mathbf{F}_I^{\mathrm{solv}} = -\partial E_{\mathrm{solv}}/\partial\boldsymbol{\tau}_I$ derive from a single stationary Lagrangian.
2. **Zero double counting & width invariance**: Smooth Gaussian ionic pseudo-charges $\rho_{\mathrm{ion}}^{\mathrm{sm}}$ reconcile algebraically with BAANDI's unmodified `state.ionic` + `ewald.rs` so results are independent of smearing width $r_{c,I}$ (provided each Gaussian core stays inside the solute cavity).
3. **Exact vacuum reduction**: At $\epsilon_b = 1, \kappa = 0, \tau = 0$ (neutral $Q_{\mathrm{net}} = 0$), every solvation energy, potential, and force vanishes identically, reducing bit-identically to vacuum BAANDI.
4. **GC-DFT and Track C readiness**: The electrostatic reference ($\phi_{\mathrm{bulk}} = 0$ at $\kappa > 0$) and typed electrostatic conventions ($\phi$ vs $U_{\mathrm{el}} = -\phi$, ADR 0017) directly support Phase B2 Grand Canonical DFT ($\Omega = F - \mu_e N_e$) and Track C electrostatic embedding.

---

## 2. Prior Work & Verified References

All references below were verified via primary bibliographic lookup (zero unverified DOIs):

- **SCCS / Quantum ESPRESSO Environ**:
  - O. Andreussi, I. Dabo, and N. Marzari (2012), *Revised self-consistent continuum solvation in electronic-structure calculations*, *J. Chem. Phys.* 136, 064102. DOI: `10.1063/1.3676407`
  - G. Fisicaro, L. Genovese, O. Andreussi, N. Marzari, and S. Goedecker (2016), *A generalized Poisson and Poisson-Boltzmann solver for electrostatic environments*, *J. Chem. Phys.* 144, 014103. DOI: `10.1063/1.4939125`
  - P. Giannozzi, O. Andreussi, et al. (2017), *Advanced capabilities for materials modelling with QUANTUM ESPRESSO*, *J. Phys.: Condens. Matter* 29, 465901. DOI: `10.1088/1361-648X/aa8f79`
- **VASPsol / Density-Dependent Linearized Poisson–Boltzmann (LPB)**:
  - K. Mathew, R. Sundararaman, K. Letchworth-Weaver, T. A. Arias, and R. G. Hennig (2014), *Implicit solvation model for density-functional study of nanocrystal surfaces and reaction pathways*, *J. Chem. Phys.* 140, 084106. DOI: `10.1063/1.4865107`
- **JDFTx / CANDLE / Grand Canonical DFT**:
  - R. Sundararaman and W. A. Goddard III (2015), *The charge-asymmetric nonlocally determined local-electric (CANDLE) solvation model*, *J. Chem. Phys.* 142, 064107. DOI: `10.1063/1.4907731`
  - R. Sundararaman, K. Letchworth-Weaver, K. A. Schwarz, D. Gunceler, Y. Ozhabes, and T. A. Arias (2017), *JDFTx: software for joint density-functional theory*, *SoftwareX* 6, 278–284. DOI: `10.1016/j.softx.2017.10.006`
  - R. Sundararaman, W. A. Goddard III, and T. A. Arias (2017), *Grand canonical electronic density-functional theory: Algorithms and applications to electrochemistry*, *J. Chem. Phys.* 146, 114104. DOI: `10.1063/1.4978411`
- **Charged Surfaces & Direct Minimization**:
  - M. Otani and O. Sugino (2006), *First-principles calculations of charged surfaces and interfaces: A plane-wave nonrepeated slab approach*, *Phys. Rev. B* 73, 115407. DOI: `10.1103/PhysRevB.73.115407`
  - C. Freysoldt, S. Boeck, and J. Neugebauer (2009), *Direct minimization technique for metals in density functional theory*, *Phys. Rev. B* 79, 241103. DOI: `10.1103/PhysRevB.79.241103`

---

## 3. Physical Background

In continuum implicit solvation, the solvent is represented by:
1. A **spatially varying local relative permittivity** $\epsilon(\mathbf{r}) = 1 + (\epsilon_b - 1) S(n_{\mathrm{cav}}(\mathbf{r}))$, transitioning smoothly from $\epsilon = 1$ inside the solute ($S = 0$) to bulk solvent $\epsilon_b$ (e.g., $\epsilon_b \approx 78.4$ for water) outside ($S = 1$).
2. An optional **mobile electrolyte screening** in the linearized Poisson–Boltzmann (LPB) regime, characterized by inverse bulk Debye length $\kappa = \sqrt{8\pi I / (\epsilon_b k_B T_{\mathrm{solv}})}$ ($\mathrm{Bohr}^{-1}$ for 1:1 electrolyte of ionic strength $I$ in $\mathrm{Bohr}^{-3}$) modulated by an ion-accessibility function $\lambda(n_{\mathrm{cav}}(\mathbf{r})) \in [0,1]$ (typically $\lambda = S$). The counterion charge density is $\rho_{\mathrm{mob}}(\mathbf{r}) = -\frac{\epsilon_b\kappa^2}{4\pi}\lambda(\mathbf{r})\phi_{\mathrm{solv}}(\mathbf{r})$.
3. A **non-electrostatic cavitation (plus dispersion/repulsion) energy** proportional to the effective cavity surface area $\int |\nabla S|\,d\mathbf{r}$ (VASPsol $\tau$, SCCS $\alpha + \gamma$) and optionally cavity volume $\int (1 - S)\,d\mathbf{r}$ (SCCS $\beta$).

Two canonical smooth cavity switching functions $S(n_{\mathrm{cav}})$ are considered:
- **VASPsol `erfc` switching** (Mathew et al. 2014):
  $$S_{\mathrm{VASPsol}}(n) = \frac{1}{2}\operatorname{erfc}\!\left(\frac{\ln(n / n_c)}{\sigma\sqrt{2}}\right), \qquad S'_{\mathrm{VASPsol}}(n) = -\frac{1}{n\,\sigma\sqrt{2\pi}}\exp\!\left(-\frac{[\ln(n/n_c)]^2}{2\sigma^2}\right) \le 0,$$
  with water defaults $n_c = 0.0025\,\mathrm{Bohr}^{-3}$, $\sigma = 0.6$, $\epsilon_b = 78.4$ (for $n \le n_{\mathrm{floor}}$, $S = 1, S' = 0$ smoothly to machine precision).
- **SCCS logarithmic-sine switching** (Andreussi et al. 2012):
  For $t(n) = \frac{\ln \rho_{\max} - \ln n}{\ln \rho_{\max} - \ln \rho_{\min}}$:
  $$S_{\mathrm{SCCS}}(n) = \begin{cases} 1, & n \le \rho_{\min}, \\ t(n) - \frac{1}{2\pi}\sin(2\pi t(n)), & \rho_{\min} < n < \rho_{\max}, \\ 0, & n \ge \rho_{\max}, \end{cases} \qquad S'_{\mathrm{SCCS}}(n) = -\frac{1 - \cos(2\pi t(n))}{n\ln(\rho_{\max}/\rho_{\min})} \le 0$$
  (or the log-$\epsilon$ variant $\ln\epsilon = (\ln\epsilon_b) S_{\mathrm{SCCS}}(n)$, which keeps $\epsilon(n)$ $C^2$-smooth across $[\rho_{\min}, \rho_{\max}]$; water defaults $\rho_{\min} = 10^{-4}\,\mathrm{Bohr}^{-3}$, $\rho_{\max} = 1\text{–}5 \times 10^{-3}\,\mathrm{Bohr}^{-3}$).

---

## 4. Mathematical Derivation (Exact BAANDI Conventions)

### 4.1 BAANDI Units, FFT Normalization, and Skew-Adjoint Derivative
All equations use Hartree atomic units ($e = m_e = \hbar = 4\pi\varepsilon_0 = 1$; `docs/physics-conventions.md:5-7`):
- Grid $(n_1,n_2,n_3)$, $N_{\mathrm{grid}} = n_1 n_2 n_3$, cell volume $\Omega = |\det A|$, quadrature weight $\Delta V = \Omega / N_{\mathrm{grid}}$ (`refscf/map.rs:49`).
- Forward FFT (`StartFields::forward`): $f(\mathbf{G}) = \frac{1}{N_{\mathrm{grid}}}\sum_{\mathbf{r}} f(\mathbf{r}) e^{-i\mathbf{G}\cdot\mathbf{r}}$; inverse: $f(\mathbf{r}) = \sum_{\mathbf{G}} f(\mathbf{G}) e^{+i\mathbf{G}\cdot\mathbf{r}}$; Parseval: $\Delta V \sum_{\mathbf{r}} f^* g = \Omega \sum_{\mathbf{G}} f^*(\mathbf{G}) g(\mathbf{G})$.
- **Discrete skew-adjoint gradient $D_a$** (`crates/dft-physics/src/gga.rs:15-47`, `gga/distributed.rs:109-163`, `docs/physics-conventions.md:1112-1117`): Fourier multiplier $+i G_a$ with each even-axis Nyquist Miller index zeroed *before* multiplying by the reciprocal lattice matrix $B$. This guarantees $D_a$ maps real grids to real grids and satisfies $D_a^\top = -D_a$ strictly under $\langle u, v\rangle = \Delta V \sum_{\mathbf{r}} u(\mathbf{r}) v(\mathbf{r})$.
- **Solute charge density & potential signs**:
  - Valence electron number density $n(\mathbf{r}) \ge 0$, $\Delta V\sum_{\mathbf{r}} n = N_e$, electric charge density $\rho_{\mathrm{el}} = -n$.
  - Smooth Gaussian ionic pseudo-charge density ($Z_I = Z_{\mathrm{val},I} > 0$, $\int_\Omega \rho_{\mathrm{ion}}^{\mathrm{sm}} d\mathbf{r} = Q_{\mathrm{ion}} = \sum_I Z_I$):
    $$\rho_{\mathrm{ion}}^{\mathrm{sm}}(\mathbf{G}) = \frac{1}{\Omega}\sum_I Z_I e^{-G^2 r_{c,I}^2 / 4} e^{-i\mathbf{G}\cdot\boldsymbol{\tau}_I}, \qquad \rho_{\mathrm{sol}}(\mathbf{r}) = \rho_{\mathrm{ion}}^{\mathrm{sm}}(\mathbf{r}) - n(\mathbf{r}), \qquad \bar{\rho}_{\mathrm{sol}} = \frac{Q_{\mathrm{ion}} - N_e}{\Omega}.$$
  - Electrostatic potential $\phi(\mathbf{r})$ (energy per $+1e$, in Ha) satisfies the discrete SPD GP/LPB equation:
    $$A\,\phi_{\mathrm{solv}} \equiv \left(\sum_{a=1}^3 D_a^\top \epsilon(\mathbf{r}) D_a + \epsilon_b\kappa^2\lambda(\mathbf{r})\right)\phi_{\mathrm{solv}}(\mathbf{r}) = 4\pi \rho_{\mathrm{sol}}(\mathbf{r}).$$
  - Electron potential energy is $U_{\mathrm{el}}(\mathbf{r}) = -\phi(\mathbf{r})$ (ADR 0017). In vacuum ($\epsilon = 1, \kappa = 0$), $-\phi_{\mathrm{vac}}(\mathbf{G}\ne\mathbf{0}) = +\frac{4\pi}{G^2}n(\mathbf{G}) - \frac{4\pi}{G^2}\rho_{\mathrm{ion}}^{\mathrm{sm}}(\mathbf{G}) = V_H(\mathbf{G}) + V_{\mathrm{ion}}^{\mathrm{sm}}(\mathbf{G})$, matching `hartree.rs:97-103` and `dft-upf/src/radial.rs:285-288` identically.

### 4.2 Exact Reconciliation with `state.ionic` and `ewald.rs` (Zero Double Counting Proof)
Decompose BAANDI's existing local UPF potential around $\rho_{\mathrm{ion},I}^{\mathrm{sm}}$:
$$V_{\mathrm{ion}}(\mathbf{G}\ne\mathbf{0}) = V_{\mathrm{ion}}^{\mathrm{sm}}(\mathbf{G}) + \Delta V_{\mathrm{loc}}^{\mathrm{short}}(\mathbf{G}), \qquad \alpha_{Z,I} = \underbrace{\frac{\pi Z_I r_{c,I}^2}{\Omega}}_{\alpha_{Z,I}^{\mathrm{sm}}} + \Delta\alpha_{Z,I}^{\mathrm{short}}.$$
By BAANDI's Section 24 Gaussian Ewald identity (`docs/physics-conventions.md:938-968`, `ewald.rs:645-672`), for **arbitrary $N_e$** (neutral or charged):
$$E_H[n] + E_{\mathrm{loc}}[n] + E_{\mathrm{ewald}} = E_{\mathrm{es}}^{\mathrm{vac}}[\tilde{\rho}_{\mathrm{sol}}] - (S - P) + E_{\mathrm{loc}}^{\mathrm{short}}[n] + (N_e - Q_{\mathrm{ion}})\sum_I \alpha_{Z,I}^{\mathrm{sm}},$$
where $\tilde{\rho}_{\mathrm{sol}} = \rho_{\mathrm{sol}} - \bar{\rho}_{\mathrm{sol}}$, $S = \sum_I Z_I^2/(\sqrt{2\pi}r_{c,I})$ is the Gaussian self-energy, and $P$ is the short-range Gaussian pair overlap.

We keep BAANDI's `state.ionic`, `hartree`, and `ewald` **100% unmodified** and define the dual variational solvation functional $E_{\mathrm{solv}} = \Delta E_{\mathrm{es}} + E_{\mathrm{bg,corr}} + E_{\mathrm{cav}}$:
1. **Dual electrostatic solvation functional** (stationary w.r.t. $\phi_{\mathrm{solv}}$ and $\phi_{\mathrm{vac}}$):
   $$\Delta E_{\mathrm{es}} = \underbrace{\left(\langle \rho_{\mathrm{sol}}, \phi_{\mathrm{solv}}\rangle - \frac{1}{8\pi}\langle \phi_{\mathrm{solv}}, A\,\phi_{\mathrm{solv}}\rangle\right)}_{E_{\mathrm{es}}^{\mathrm{solv}}} - \underbrace{\left(\langle \tilde{\rho}_{\mathrm{sol}}, \phi_{\mathrm{vac}}\rangle - \frac{1}{8\pi}\langle \phi_{\mathrm{vac}}, L_0\,\phi_{\mathrm{vac}}\rangle\right)}_{E_{\mathrm{es}}^{\mathrm{vac}}},$$
   where $L_0 = \sum_a D_a^\top D_a$ uses the **same skew-adjoint $D_a$** as $A$, and $\langle\phi_{\mathrm{vac}}\rangle = 0$. On-shell ($\phi_{\mathrm{solv}} = A^{-1}(4\pi\rho_{\mathrm{sol}})$, $\phi_{\mathrm{vac}} = L_0^+(4\pi\tilde{\rho}_{\mathrm{sol}})$), with $\Delta\phi_{\mathrm{pol}} \equiv \phi_{\mathrm{solv}} - \phi_{\mathrm{vac}}$:
   $$\Delta E_{\mathrm{es}} = \frac{1}{2}\langle \rho_{\mathrm{sol}}, \Delta\phi_{\mathrm{pol}}\rangle.$$
2. **Exact $G=0$ background-uncompensation correction**:
   $$E_{\mathrm{bg,corr}} = +(Q_{\mathrm{ion}} - N_e)\sum_I \alpha_{Z,I}^{\mathrm{sm}} = +(Q_{\mathrm{ion}} - N_e)\sum_I \frac{\pi Z_I r_{c,I}^2}{\Omega}.$$
   Summing $E_H + E_{\mathrm{loc}} + E_{\mathrm{ewald}} + \Delta E_{\mathrm{es}} + E_{\mathrm{bg,corr}} = E_{\mathrm{es}}^{\mathrm{solv}}[\rho_{\mathrm{sol}}] - (S-P) + E_{\mathrm{loc}}^{\mathrm{short}}[n]$ cancels all bare Coulomb terms and $(N_e - Q_{\mathrm{ion}})\sum_I \alpha_{Z,I}^{\mathrm{sm}}$ **identically**, guaranteeing strict $r_{c,I}$-invariance for both neutral and charged cells. Moreover, at $\kappa > 0$, $\delta E_{\mathrm{bg,corr}}/\delta n = -\sum_I \alpha_{Z,I}^{\mathrm{sm}}$ cancels BAANDI's $V_{\mathrm{ion}}(\mathbf{0})$ offset so $V_{\mathrm{eff}}(\mathbf{r}_{\mathrm{bulk}}) = 0$ exactly, aligning $\mu_e$ to the bulk-electrolyte reference.
3. **Polarization-potential PDE & exact vacuum limit**:
   Substituting $\phi_{\mathrm{solv}} = \phi_{\mathrm{vac}} + \Delta\phi_{\mathrm{pol}}$ yields:
   $$A\,\Delta\phi_{\mathrm{pol}} = 4\pi\bar{\rho}_{\mathrm{sol}} - \sum_a D_a^\top\big[(\epsilon - 1)D_a\phi_{\mathrm{vac}}\big] - \epsilon_b\kappa^2\lambda\phi_{\mathrm{vac}}.$$
   Solving directly for $\Delta\phi_{\mathrm{pol}}$ avoids $O(Z/r_c)$ core-peak subtraction cancellation in $\Delta E_{\mathrm{es}}$, and at $\epsilon_b = 1, \kappa = 0, \tau = 0, Q_{\mathrm{net}} = 0$ the RHS is identically zero ($\Delta\phi_{\mathrm{pol}} \equiv 0$, exact vacuum reduction). Total $\phi_{\mathrm{solv}} = \phi_{\mathrm{vac}} + \Delta\phi_{\mathrm{pol}}$ is then reconstructed for the nonlinear cavity derivatives below.
4. **Regularized discrete cavitation energy**:
   $$E_{\mathrm{cav}} = \tau\,\Delta V \sum_{\mathbf{r}}\left(\sqrt{\sum_{a=1}^3 (D_a S(\mathbf{r}))^2 + \delta_S^2} - \delta_S\right) + \beta\,\Delta V\sum_{\mathbf{r}}\big(1 - S(\mathbf{r})\big),$$
   with tiny regularization $\delta_S > 0$ (e.g., $10^{-14}\,\mathrm{Bohr}^{-1}$) ensuring smooth differentiability everywhere ($E_{\mathrm{cav}} = 0$ at $\tau = \beta = 0$).

### 4.3 Variational Potential $V_{\mathrm{solv}} = \delta E_{\mathrm{solv}}/\delta n$
By stationarity of the dual Lagrangian w.r.t. $\phi_{\mathrm{solv}}$ and $\phi_{\mathrm{vac}}$, implicit $\delta\phi/\delta n$ terms vanish identically:
$$V_{\mathrm{solv}}(\mathbf{r}) = V_{\mathrm{es}}^{\mathrm{solv}}(\mathbf{r}) + V_{\mathrm{diel}}(\mathbf{r}) + V_{\mathrm{scr}}(\mathbf{r}) + V_{\mathrm{cav}}(\mathbf{r}),$$
where:
- $V_{\mathrm{es}}^{\mathrm{solv}}(\mathbf{r}) = -\Delta\phi_{\mathrm{pol}}(\mathbf{r}) - \delta_{\kappa>0}\sum_I \alpha_{Z,I}^{\mathrm{sm}}$
- $V_{\mathrm{diel}}(\mathbf{r}) = -\frac{1}{8\pi}\frac{d\epsilon}{dn_{\mathrm{cav}}}\sum_a (D_a \phi_{\mathrm{solv}})^2 \ge 0$
- $V_{\mathrm{scr}}(\mathbf{r}) = -\frac{\epsilon_b\kappa^2}{8\pi}\frac{d\lambda}{dn_{\mathrm{cav}}}\phi_{\mathrm{solv}}^2 \ge 0$
- $V_{\mathrm{cav}}(\mathbf{r}) = -\tau\,S'(n_{\mathrm{cav}})\sum_a D_a\!\left(\frac{D_a S}{\sqrt{|\mathbf{D}S|^2 + \delta_S^2}}\right) - \beta\,S'(n_{\mathrm{cav}})$

**Critical discrete identity**: Differentiate $D_a S$ directly (never substitute $S'(n) D_a n$) so $V_{\mathrm{cav}}$ is the *exact* discrete gradient of $E_{\mathrm{cav}}$ on the finite FFT grid (matching BAANDI's GGA convention at `gga.rs:226-231`).

### 4.4 Variational Ionic Forces $\mathbf{F}_I^{\mathrm{solv}} = -\partial E_{\mathrm{solv}}/\partial\boldsymbol{\tau}_I$
At SCF stationarity (fixed-$N$ or fixed-$\mu$), all implicit orbital, occupation, and $\phi$ responses vanish:
1. **Gaussian ionic pseudo-charge force** (always active):
   $$F_{I,\alpha}^{\mathrm{solv,ion}} = -\Omega \sum_{\mathbf{G}\in\mathcal{B}_{\mathrm{dens}}} G_\alpha\,\mathrm{Im}\!\left[\Delta\phi_{\mathrm{pol}}^*(\mathbf{G})\,\frac{Z_I}{\Omega} e^{-G^2 r_{c,I}^2 / 4}\,e^{-i\mathbf{G}\cdot\boldsymbol{\tau}_I}\right],$$
   identical in form to `dft_physics::forces::local_or_core_cached` (`forces.rs:90-95`).
2. **NLCC core-cavity force** (active ONLY when $n_{\mathrm{cav}} = n_{\mathrm{val}} + n_{\mathrm{core}}$ with NLCC):
   $$\mathbf{F}_I^{\mathrm{solv,core}} = -\int_\Omega \big(V_{\mathrm{diel}} + V_{\mathrm{scr}} + V_{\mathrm{cav}}\big)\,\frac{\partial n_{\mathrm{core},I}}{\partial\boldsymbol{\tau}_I}\,d\mathbf{r},$$
   reusing `forces::local_or_core_cached(..., /* core = */ true, ...)` verbatim (note $V_{\mathrm{es}}^{\mathrm{solv}}$ is excluded since $n_{\mathrm{core}}$ is not part of $\rho_{\mathrm{sol}}$).

### 4.5 `EnergyAudit`, Harris–Foulkes Remainder, and GC-DFT ($\Omega = F - \mu_e N_e$)
- Add `solvation: Option<SolvationEnergyAudit>` to `EnergyAudit` (`energy.rs:8-29`), folding `vsolv` into `hxc` (and subtracting `vsolv` from the `dvxc` PBE diagnostic at `electronic_map.rs:1033-1035`).
- Double counting: $D = \sum_s \int n_{\mathrm{out},s}(V_H + V_{\mathrm{xc},s} + V_{\mathrm{solv}})\,d\mathbf{r}$; `band_form` $= B - D + E_H + E_{\mathrm{xc}} + E_{\mathrm{solv}} + E_{\mathrm{ewald}} \equiv E_{\mathrm{direct}}$.
- Harris–Foulkes & Bregman remainder (`energy.rs:339-380`):
  $$E_{\mathrm{HF}}[n_{\mathrm{in}}] = B - E_H[n_{\mathrm{in}}] + \big(E_{\mathrm{xc}} - \langle n_{\mathrm{in}}, V_{\mathrm{xc}}\rangle\big) + \big(E_{\mathrm{solv}} - \langle n_{\mathrm{in}}, V_{\mathrm{solv}}\rangle\big) + E_{\mathrm{ewald}},$$
  $$Q_{\mathrm{solv}} = E_{\mathrm{solv}}[n_{\mathrm{out}}] - E_{\mathrm{solv}}[n_{\mathrm{in}}] - \langle n_{\mathrm{out}} - n_{\mathrm{in}}, V_{\mathrm{solv}}[n_{\mathrm{in}}]\rangle.$$
- **Phase B2 GC-DFT readiness**: When $\kappa > 0$ fixes $\phi_{\mathrm{bulk}} = 0$, adding `FermiDirac::evaluate_fixed_mu` in `occupations.rs` (reusing `fraction_entropy` at lines 59–66), tracking `grand_potential_ha = F - mu * N_e` in `dft_engine::model::ThermodynamicLedger`, using $N_{\mathrm{out}}$ (not fixed $Q_{\mathrm{ion}}$) in `local_zero_mode` (`map.rs:531`, `electronic_map.rs:1089`), and letting the $G=0$ density mode mix freely satisfies $\partial\Omega/\partial\mu = -N_e$ and $\mathbf{F}_I = -(\partial\Omega/\partial\boldsymbol{\tau}_I)_{\mu,T}$ identically.

---

## 5. Numerical Algorithms

Three distributed G-space solvers for $A\,\phi = b$ on BAANDI's `DistributedPairFft` are analyzed in detail in [algorithm-comparison.md](algorithm-comparison.md):
1. **Algorithm 1 (Primary) — Direct G-space PCG with Real-Space & Mirror Complex Pairing**: Strictly SPD on any grid, mesh-independent $\kappa_2(M^{-1}A) \le \epsilon_b$, 0 setup FFTs, **4 FFTs/iter** under BAANDI's default `mirror` stick partition (`np = 2` batched `G -> r` and `r -> G`).
2. **Algorithm 2 (Benchmark Arm) — Fisicaro $\sqrt{\epsilon}$-PCG** (Fisicaro et al. 2016): Solves $(-\nabla^2 + q)\tilde{\phi} = 4\pi\rho/\sqrt{\epsilon}$ at 2 FFTs/iter, but pays 3 setup FFTs/SCF step and risks indefiniteness/CG breakdown at sharp cavity interfaces where $q = \nabla^2\sqrt{\epsilon}/\sqrt{\epsilon}$ has large negative lobes.
3. **Algorithm 3 (Benchmark Arm) — Polarization-Charge Fixed-Point + Anderson/DIIS**: 3 FFTs/iter (paired `G -> r`), 4 setup FFTs/SCF step, non-symmetric operator with potential convergence degradation at large $\|\nabla\ln\epsilon\|_\infty$.

---

## 6. Computational Complexity

Let $N = N_0 N_1 N_2$ be the density FFT grid size, $P$ the MPI rank count, and $K_{\mathrm{PCG}}$ the inner PCG iterations per SCF step:
- **Per PCG iteration**: $O\!\left(\frac{N}{P}\log N\right)$ FFT arithmetic + 2 batched `alltoall_f64` (payload $4U$ for mirror-paired Alg 1) + 2 tiny `tree_sum_f64` (1 and 2 `f64` scalars).
- **Per SCF step with warm-start** (0-FFT bulk-delta predictor + adaptive inner tolerance $\tau_{\mathrm{inner}} \propto \sqrt{\mathrm{dr2\_ry}}$): $K_{\mathrm{PCG}} \approx 2\text{–}6$ iterations $\implies$ **8–24 full-grid FFTs/SCF step**, versus 200–500 full-grid FFT equivalents in Davidson + GEMMs for typical 64–128 band / multi-$k$ runs (**< 3–5% SCF wall-time overhead**, and **0% when disabled**).
- **Memory footprint**: $\sim 128\text{–}176\,(N/P)$ bytes/rank for Alg 1 workspace (e.g., $\sim 20\text{–}25$ MB/rank on a $100^3$ grid at $P=8$), allocated once per run; zero heap allocations in inner iterations.

---

## 7. Key Approximations & Applicability Limits

1. **Local linear dielectric response**: Ignores dielectric saturation at extreme interfacial electric fields ($|\mathbf{E}| \gtrsim 1\,\mathrm{V/\AA}$) and nonlocal hydrogen-bonding structure (first solvation shell).
2. **Linearized Poisson–Boltzmann ($\kappa > 0$)**: Valid for moderate interfacial potentials ($|e\phi| \lesssim k_B T \approx 26$ mV, with reasonable qualitative behavior up to $\sim 0.3\text{–}0.5$ V); high-potential double layers require nonlinear PB or modified steric PB.
3. **Cavity-charge non-overlap**: Exact $r_{c,I}$-invariance requires $\rho_{\mathrm{ion},I}^{\mathrm{sm}}(\mathbf{r})$ to decay before reaching the cavity boundary ($S > 0$). Typical $r_{c,I} \approx 0.5\text{–}0.8\,\mathrm{Bohr}$ satisfies this for standard atoms; pathological tiny-cutoff grids must be rejected at preflight.
4. **Periodic image interactions**: At $\kappa = 0$, a periodic cell still has dipolar/quadrupolar solvent-screened images (attenuated by $1/\epsilon_b$, not eliminated); at $\kappa > 0$, images decay exponentially as $e^{-\kappa L}/L$.
5. **Stress deferred**: Analytic solvation stress ($\partial E_{\mathrm{solv}}/\partial\varepsilon_{\alpha\beta}$) is explicitly marked `Availability::Unsupported` in Phase B1/Gate 2 rather than returning placeholder zeros.

---

## 8. Implementation Candidates in BAANDI

- **Gate 2 Prototype (isolated in `dft-physics`, zero `dft-engine` edits)**:
  `crates/dft-physics/src/solvation/{mod.rs, cavity.rs, poisson.rs, tests.rs}` implementing all three solver arms (`DirectPcg`, `FisicaroPcg`, `PolarizationDiis`), both cavity models (`VaspSol`, `Sccs`), discrete cavitation energy/potential, Gaussian ionic pseudo-charge energy/force, and `dft_mpi::simulate` multi-rank verification.
- **Gate 3–4 Integration (after Gate 2 Go decision)**:
  - Add typed semantic wrappers (`ElectrostaticPotentialHaPerE`, `ElectronEnergyHa`, `SignedChargeDensityEPerBohr3`, ADR 0017) in `dft-units`/`dft-physics`.
  - Wire `Option<Box<SolvationWorkspace>>` into `refstart::StartFields` and `refscf/{map.rs, k_map.rs, electronic_map.rs}`, `energy.rs::EnergyAudit`, and `forces.rs` (`DerivativeEntry`).
  - Expose `[solvation]` TOML table in `dft-cli/src/config.rs` (`deny_unknown_fields`).

---

## 9. Verification Plan (Levels 1–5, Pre-Declared Tolerances)

All tolerances below are fixed **before** Gate 2 coding and must not be relaxed post-hoc:

| Level | Oracle / Test | Pre-Declared Pass Criterion |
| --- | --- | --- |
| **L1 — Exact Vacuum Limit** | $\epsilon_b = 1, \kappa = 0, \tau = 0, Q_{\mathrm{net}} = 0$ on arbitrary skewed cell | $\|\Delta\phi_{\mathrm{pol}}\|_\infty < 10^{-14}$, $|\Delta E_{\mathrm{es}}| < 10^{-14}$ Ha, $\|V_{\mathrm{solv}}\|_\infty < 10^{-14}$ Ha, $\|\mathbf{F}^{\mathrm{solv}}\|_\infty < 10^{-14}$ Ha/Bohr |
| **L1 — Method of Manufactured Solutions (MMS)** | Analytic 3D periodic $\epsilon(\mathbf{r})$ ($\epsilon_b = 78.4$), $\kappa^2(\mathbf{r}) > 0$, $\phi_{\mathrm{exact}}(\mathbf{r})$ on cubic & skewed cells | $\|\phi_{\mathrm{num}} - \phi_{\mathrm{exact}}\|_\infty / \|\phi_{\mathrm{exact}}\|_\infty < 10^{-12}$ (both $\kappa = 0$ and $\kappa > 0$) |
| **L1 — Analytic Born & Debye–Hückel Limits** | Spherical Gaussian charge in spherical dielectric cavity ($\epsilon_b = 78.4$, $\kappa = 0$ and $\kappa > 0$) across box sizes $L$ | Monotonic $O(1/L)$ (or $e^{-\kappa L}/L$) convergence to analytic $\Delta G_{\mathrm{Born}} = -\frac{q^2}{2R}(1 - 1/\epsilon_b)$ and Debye–Hückel limit within $< 1\%$ at converged $L$ |
| **L2 — Variational Derivative FD Check** | Central FD of $E_{\mathrm{solv}}[n \pm \delta n]$ vs $V_{\mathrm{es}}^{\mathrm{solv}} + V_{\mathrm{diel}} + V_{\mathrm{scr}} + V_{\mathrm{cav}}$ on random grid points (skewed cell, VASPsol & SCCS) | Relative error $< 10^{-9}$ (absolute $< 10^{-10}$ Ha) for every term individually and summed |
| **L2 — Ionic Force & Width-Invariance Check** | Central FD of $E_{\mathrm{solv}}(\boldsymbol{\tau}_I \pm \delta)$ vs $\mathbf{F}_I^{\mathrm{solv,ion}}$; $r_c$ sweep inside cavity | Force FD error $< 10^{-9}$ Ha/Bohr; $r_c$ sensitivity $< 10^{-6}$ Ha when Gaussian tail $< 10^{-8}$ at cavity edge |
| **L2 — MPI & Stick-Partition Determinism** | `dft_mpi::simulate` for $P \in \{1, 2, 4\}$, comparing `mirror` vs `greedy` stick partitions | Cross-rank agreement $< 10^{-12}$ (energy & max field diff); `mirror` vs `greedy` agreement $< 10^{-12}$ |
| **L3–L4 — Electronic Structure & Physical Validation (Gate 4)** | Neutral molecules ($\text{H}_2\text{O}, \text{CO}, \text{CH}_4, \text{NH}_3$) & slab ($\text{Pt}(111)$ / $\text{Al}(111)$) SCF hydration energies, forces, and work-function shift | SCF FD forces $< 10^{-6}$ Ha/Bohr; hydration energies consistent with published VASPsol/SCCS benchmarks (Mathew 2014 / Andreussi 2012) |
| **L5 — Performance Evaluation (Gate 2 & Gate 4 on `dev08`)** | PCG wall time, iteration counts (cold vs warm-start), MPI scaling, vacuum off-path regression | Warm-started solvation $< 5\%$ of SCF wall time; off-path pure-PW wall time & bits **0% regression** |

*(External oracle status verified on `dev08` in Gate 1: QE 7.6 is built with `environ_api_stub.f90` only—Environ is not yet compiled; JDFTx and GPAW are not installed. Therefore L1–L2 mathematical/FD/Born oracles + published Mathew 2014 / Andreussi 2012 reference values serve as Gate 2–4 primary gates, with optional QE+Environ compilation on `dev08` tracked for Gate 4 cross-code checks.)*

---

## 10. Expected Failure Modes & Mitigations

1. **CG breakdown in Fisicaro $\sqrt{\epsilon}$-PCG (Alg 2) on sharp cavities**: High-$G$ aliasing of $\nabla^2\sqrt{\epsilon}$ creates negative eigenvalues $\implies$ *Mitigation*: Use strictly SPD Direct G-space PCG (Alg 1) as primary solver.
2. **Non-physical cavity penetration by Gaussian ionic core**: If $r_{c,I}$ is too large relative to bond/cavity radius, $\rho_{\mathrm{ion}}^{\mathrm{sm}}$ leaks into $\epsilon > 1$ and breaks $r_c$-invariance $\implies$ *Mitigation*: Automatic preflight diagnostic checking $\max_{\mathbf{r}} \rho_{\mathrm{ion}}^{\mathrm{sm}}(\mathbf{r}) S(\mathbf{r})$ and default conservative $r_{c,I} \approx 0.5\,\mathrm{Bohr}$ (subject to grid Nyquist check $r_{c,I} \ge 2.5\,\Delta x$).
3. **Broyden SCF oscillation from steep $V_{\mathrm{diel}} \propto \frac{d\epsilon}{dn}|\nabla\phi|^2$**: On grids with small negative density ringing in vacuum (`G4E-VACUUM-DENSITY`), $\ln(n/n_c)$ is undefined if $n \le 0$ $\implies$ *Mitigation*: Smooth floor ($n \le n_{\mathrm{floor}} \implies S = 1, S' = 0$) preventing NaN/div-by-zero; test Broyden stability and evaluate damped dielectric updates / direct minimization in Gate 4 if stalling occurs.
4. **Even-Nyquist mismatch between `hartree.rs` and `gga.rs`**: Using unzeroed Nyquist modes in $\phi_{\mathrm{vac}}$ while using skew-adjoint $D_a$ in $A$ breaks exact vacuum cancellation $\implies$ *Mitigation*: Always construct $\phi_{\mathrm{vac}}$ with the identical skew-adjoint $L_0 = \sum_a D_a^\top D_a$ operator.

---

## 11. BAANDI Applicability & High-Performance Design

Implicit solvation is an ideal architectural fit for BAANDI:
- Reuses `DistributedPairFft`, `StartFields`, `gga/distributed.rs` skew-adjoint derivative convention, `Collective::tree_sum_f64`, and `forces::local_or_core_cached`.
- Zero changes to `DistributedRealHamiltonian`, `DistributedComplexHamiltonian`, `ref_davidson`, `ref_davidson_complex`, or `KeptHx`.
- Complex pairing (`np = 2`) + BAANDI's default `mirror` stick layout + 0-FFT bulk-delta warm-start predictor make BAANDI's generalized Poisson solver faster per iteration than textbook 6–8 FFT implementations while keeping memory flat and allocations strictly zero.

---

## 12. Researcher Judgment & Open Questions for Gate 2

- **Judgment**: The theory is mathematically closed, variationally exact, and directly implementable in safe Rust (`unsafe_code = "deny"`) inside `dft-physics::solvation` without touching `ownership.toml` or `dft-engine` in Gate 2.
- **Open empirical questions to resolve in Gate 2 (prototype & `dev08` benchmarks)**:
  1. Actual iteration counts and wall times of **Alg 1 (Direct mirror-paired PCG)** vs **Alg 2 (Fisicaro $\sqrt{\epsilon}$-PCG)** vs **Alg 3 (Polarization-Charge DIIS)** on realistic sharp molecular/slab density cavities at $\epsilon_b = 78.4$ and $\kappa \in \{0, 0.1, 0.3\}\,\mathrm{Bohr}^{-1}$.
  2. Exact speedup of the **0-FFT bulk-delta warm-start predictor** across realistic density perturbations ($\|\delta n\| / \|n\| \sim 10^{-1} \to 10^{-6}$).
  3. Whether `VaspSol` (`erfc`) vs `Sccs` (log-sine) shows better discrete cavitation smoothness and lower PCG iteration counts on typical PW density grids ($E_{\mathrm{cut}} = 40\text{–}80$ Ry).
