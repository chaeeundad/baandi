# Track B Gate 1: Generalized Poisson / LPB Algorithm & Performance Comparison

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Track B — Phase B1)
- **Phase**: Gate 1 — Algorithm Comparison & Performance Engineering Plan
- **Date**: 2026-11-28
- **Target PDE (SPD form, Hartree a.u.)**:
  $$A\,\phi \equiv \left(\sum_{a=1}^3 D_a^\top \epsilon(\mathbf{r}) D_a + \epsilon_b\kappa^2\lambda(\mathbf{r})\right)\phi(\mathbf{r}) = 4\pi\rho_{\mathrm{tot}}(\mathbf{r})$$
  where $D_a$ is BAANDI's skew-adjoint Nyquist-zeroed spectral derivative (`dft_physics::gga::distributed::multiplier`).

---

## 1. BAANDI Distributed Primitives & Zero-Allocation Strategy

All candidate solvers run strictly inside `dft-physics` (`unsafe_code = "deny"`, zero new crate dependencies, zero changes to frozen `ownership.toml` or `Cargo.toml`):

- **FFT engine**: Directly borrow `&mut DistributedPairFft` (`crates/dft-solvers/src/distributed_fft.rs`), pre-warmed once via `reserve_batch(np)`. (`StartFields::{forward, inverse}` allocate a `Vec` per call, so hot solver iterations operate directly on `DistributedPairFft`.)
- **Real-space complex pairing (key BAANDI optimization)**:
  - Because $\phi(\mathbf{r})$ is real, packing $\hat{u}_0(\mathbf{G}) = (i G_x - G_y)\hat{\phi}(\mathbf{G})$ and $\hat{u}_1(\mathbf{G}) = i(G_z + 1)\hat{\phi}(\mathbf{G})$ (`np = 2`) recovers $\partial_x\phi = \operatorname{Re} u_0$, $\partial_y\phi = \operatorname{Im} u_0$, $\partial_z\phi = \operatorname{Re} u_1$, and $\phi = \operatorname{Im} u_1$ in **2 inverse FFTs + 1 batched `alltoall_f64`** on **any** stick partition (`mirror` or `greedy`).
  - On the reverse leg (`r -> G`, Alg 1), packing $(\epsilon\partial_x p + i\epsilon\partial_y p)$ and $(\epsilon\partial_z p + i\epsilon_b\kappa^2\lambda p)$ (`np = 2`) and unpacking locally via $\hat{a}(\mathbf{G}) = \frac{1}{2}(\hat{u}(\mathbf{G}) + \hat{u}(-\mathbf{G})^*)$ uses **2 forward FFTs + 1 batched `alltoall_f64`** whenever $-\mathbf{G}$ is rank-local (`BAANDI_M2_STICKS=mirror`, BAANDI default), with automatic exact fallback (`np = 3` for $\kappa=0$, `np = 4` for $\kappa > 0$) on `greedy` partitions.
- **G-space formulation**: Running PCG directly on owned reciprocal sticks (`local_positions()`) makes the diagonal preconditioner $M^{-1}(\mathbf{G}) = (\bar{\epsilon} G^2 + \bar{\kappa}^2)^{-1}$ cost **0 FFTs and 0 MPI collectives** per iteration.
- **Deterministic MPI reductions**: All inner products use `Collective::tree_sum_f64` (ADR 0041) on fixed caller-owned stack arrays (zero allocations).

---

## 2. Head-to-Head Algorithm Comparison (Per Inner Iteration & Per SCF Step)

Let $N_{\mathrm{loc}} \approx N/P$ be owned grid points/sticks per rank, and $U = 2\cdot\text{block}\cdot P$ `f64`s be the `alltoall_f64` payload per complex field.

| Metric | **Alg 1: Direct G-Space PCG** (`mirror`-paired / `greedy`-paired) | **Alg 2: Fisicaro $\sqrt{\epsilon}$-PCG** (Fisicaro et al. 2016) | **Alg 3: Polarization-Charge + Anderson/DIIS** ($m=8$) |
| --- | --- | --- | --- |
| **3D FFTs / iter ($\kappa = 0$)** | **`4` (2i+2f)** / `5` (2i+3f) | **`2`** (1i+1f via `apply_local_batch_planes`) | **`3`** (2i+1f, paired `G -> r`) |
| **3D FFTs / iter ($\kappa > 0$)** | **`4` (2i+2f)** / `6` (2i+4f) | **`2`** (1i+1f) | **`3`** (2i+1f) |
| **`alltoall_f64` calls / iter** | **2** (0 on $P=1$ direct) | **2** (0 on $P=1$ direct) | **2** (0 on $P=1$ direct) |
| **`alltoall_f64` payload / iter** | **`4U`** (`mirror`) / `5U–6U` (`greedy`) | **`2U`** | **`3U`** |
| **`tree_sum_f64` calls / iter** | **2** (1 `f64` + fused 2 `f64`) | **2** (1 `f64` + fused 2 `f64`) | **1** ($m+1 \approx 9$ `f64`s) |
| **Setup FFTs per SCF step** (when $\epsilon[n]$ updates) | **`0`** ($\epsilon, \lambda$ pointwise; reuse $\hat{n}$ from `map.rs`) | **`3`** (2 for $\nabla^2\sqrt{\epsilon}$ + 1 for $4\pi\rho/\sqrt{\epsilon}$) | **`4`** (1 fwd $\ln\epsilon$ + 2 inv $\nabla\ln\epsilon$ + 1 fwd $4\pi\rho/\epsilon$) |
| **Post-solve FFTs** (for $\phi(\mathbf{r})$ & $|\nabla\phi|^2$) | **2** paired inv (1 `alltoall`) | **2** paired inv + $\sqrt{\epsilon}$ chain rule | **0–2** |
| **Workspace memory / rank** | $\sim 128\text{–}176\,N_{\mathrm{loc}}$ B | $\sim 104\,N_{\mathrm{loc}}$ B | $\sim 360\text{–}420\,N_{\mathrm{loc}}$ B |
| **Heap allocations / iter** | **0** | **0** | **0** (`RealEighWorkspace` preallocated) |
| **Discrete SPD & Stability** | **Strictly SPD on every grid**; mesh-independent $\kappa_2(M^{-1}A) \le \epsilon_b$; exact discrete variational match | **Indefinite risk**: $q = \frac{\nabla^2\sqrt{\epsilon}}{\sqrt{\epsilon}} + \frac{\epsilon_b\kappa^2\lambda}{\epsilon}$ has negative lobes ($\propto -\delta^{-2}$) at sharp cavities; CG can break down; $\kappa=0$ nullspace is $\sqrt{\epsilon}$, not $G=0$ | **Non-symmetric**: convergence degrades with $\|\nabla\ln\epsilon\|_\infty$; DIIS Gram can become singular; aliasing mismatch vs variational $\nabla\cdot(\epsilon\nabla\phi)$ |
| **Total FFTs for warm-started SCF step ($K=4$ iters)** | **$0 + 4\times 4 + 2 = 18$ FFTs** | $3 + 4\times 2 + 2 = 13$ FFTs *(if stable & same $K$)* | $4 + 4\times 3 + 0 = 16$ FFTs *(if same $K$)* |

### Analytical Selection
- **Primary candidate: Algorithm 1 (Direct G-Space PCG with Real-Space & Mirror Pairing)**.
  1. Guaranteed strict SPD ($A = \sum_a D_a^\top \epsilon D_a + \epsilon_b\kappa^2\lambda \succeq 0$) with mesh- and cavity-sharpness-independent condition number $\kappa_2(M^{-1}A) \le \epsilon_b \approx 78.4$ ($\sqrt{\kappa_2} \le 8.85$).
  2. Zero setup FFTs per SCF step — critical in warm-started SCF where $K \approx 2\text{–}6$, where Alg 2/3's 3–4 setup FFTs erase much of their per-iteration advantage.
  3. Exact discrete consistency with `EnergyAudit` and cavity forces ($V_{\mathrm{diel}}, V_{\mathrm{scr}}, V_{\mathrm{cav}}$).
- **Gate 2 empirical arms**: Implement **all three** (`DirectPcg`, `FisicaroPcg`, `PolarizationDiis`) in `dft-physics::solvation` behind `pub enum SolverKind` and benchmark them head-to-head on identical manufactured solutions and sharp $\epsilon_b = 78.4$ molecular cavities on `dev08`.

---

## 3. Cavity Parametrization & SCF Coupling Comparison

| Choice | Option A | Option B | Gate 1–2 Plan |
| --- | --- | --- | --- |
| **Switching function $S(n)$** | **VASPsol `erfc`** ($n_c = 2.5\times 10^{-3}$, $\sigma = 0.6$; $C^\infty$ for $n > 0$) | **SCCS log-sine** ($\rho_{\min}=10^{-4}$, $\rho_{\max}=1.5\times 10^{-3}$; compact support) | Implement both in `solvation/cavity.rs`; compare PCG iteration counts & cavitation smoothness in Gate 2 |
| **NLCC density argument** | Valence-only $n_{\mathrm{val}}$ (zero core-cavity force) | Total $n_{\mathrm{val}} + n_{\mathrm{core}}$ (adds `local_or_core_cached` core-cavity force) | Support via explicit `CavityDensityMode` enum; test valence-only first in Gate 2, audit NLCC in Gate 4 |
| **SCF outer coupling** | BAANDI Broyden (`refscf/broyden.rs`) with $V_{\mathrm{solv}}$ folded into `hxc` | Damped dielectric / direct free-energy minimization (Freysoldt 2009 / Sundararaman 2017) | Target Broyden first in Gate 3–4; reserve direct minimization as fallback if charged-slab / GC-DFT Broyden stalls |

---

## 4. High-Performance Engineering Plan (< 3–5% SCF Overhead, 0% Off-Path)

1. **Zero-FFT bulk-delta warm-start predictor across SCF steps**:
   Retain $\hat{\phi}^{(t-1)}(\mathbf{G})$ and $\hat{\rho}_{\mathrm{tot}}^{(t-1)}(\mathbf{G})$ across SCF steps; start step $t$ from:
   $$\hat{\phi}^{(t,0)}(\mathbf{G}) = \hat{\phi}^{(t-1)}(\mathbf{G}) + \frac{4\pi\big(\hat{\rho}_{\mathrm{tot}}^{(t)}(\mathbf{G}) - \hat{\rho}_{\mathrm{tot}}^{(t-1)}(\mathbf{G})\big)}{\bar{\epsilon} G^2 + \bar{\kappa}^2},$$
   with geometric-mean $\bar{\epsilon} = \sqrt{\epsilon_{\min}\epsilon_{\max}}$ (optimal for $\kappa_2(M^{-1}A)$) and mean $\bar{\kappa}^2 = \frac{1}{N}\sum_{\mathbf{r}}\epsilon_b\kappa^2\lambda(\mathbf{r})$ (exact $G=0$ scaling).
2. **Adaptive inner tolerance ratchet** (mirroring `refscf::profile::next_ethr_ry`):
   $$\tau_{\mathrm{inner}}^{(t)} = \operatorname{clamp}\!\left(10^{-2}\cdot\min(\sqrt{\mathrm{dr2\_ry}},\, \tau_{\mathrm{prev}}),\; \tau_{\min},\; 10^{-4}\right),$$
   tightening to $\tau_{\mathrm{final}} = 10^{-11}\text{–}10^{-12}$ on the verification map and post-SCF force pass.
3. **Fused G-space updates & single post-solve batched transform**:
   Fuse $\phi \mathrel{+}= \alpha p$, $r \mathrel{-}= \alpha w$, $z = M^{-1}r$, and both dot products $\langle r, z\rangle, \|r\|_2^2$ into one loop + one 2-scalar `tree_sum_f64`. After convergence, a single `np = 2` batched inverse FFT yields both $\phi(\mathbf{r})$ and $|\nabla\phi(\mathbf{r})|^2$ for $V_{\mathrm{diel}}$ simultaneously.
4. **Zero off-path cost**: `Option<Box<SolvationWorkspace>>` on SCF state ensures 0 bytes, 0 plans, and bit-identical arithmetic whenever solvation is disabled.

---

## 5. Gate 2 Go / Revise / Hold / Stop Decision Criteria (Pre-Declared)

Before writing Gate 2 prototype code, the exit gate is fixed as:
- **Go (advance to Gate 3 Architecture Proposal)**:
  1. Primary solver (`DirectPcg`) passes **all** L1–L2 oracles in [implicit-solvation-theory.md](implicit-solvation-theory.md#9-verification-plan-levels-15-pre-declared-tolerances) (vacuum limit $< 10^{-14}$, MMS $< 10^{-12}$, Born/Debye–Hückel analytic convergence, variational potential & ionic force FD $< 10^{-9}$, MPI $P \in \{1,2,4\}$ & `mirror`/`greedy` parity $< 10^{-12}$).
  2. Cold-start PCG at $\epsilon_b = 78.4$ converges in $\le 35$ iterations; warm-started ($10^{-3}$ density perturbation) converges in $\le 10$ iterations with **0 heap allocations** in the inner loop.
  3. Wall-time benchmark on `dev08` confirms warm-started solve is $< 5\%$ of a comparable BAANDI SCF map, and `scripts/check.sh full` + MPI+MKL gate pass cleanly with zero regressions.
- **Revise**: L1/L2 accuracy passes, but PCG iteration count or wall time misses target (e.g., evaluate improved preconditioner or hybrid Alg 1/Alg 2 switching).
- **Hold / Stop**: Fundamental numerical instability, irreconcilable variational defect, or prohibitive overhead ($> 25\%$ of SCF wall time) with no viable algorithmic remedy.
