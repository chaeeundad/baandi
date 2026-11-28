# Track B Gate 3: Performance & Scaling Benchmark Plan (Level 5)

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Track B — Phase B1)
- **Phase**: Gate 3 — Pre-Declared Performance Targets, Profiling & Regression Protocol
- **Date**: 2026-11-28
- **Hardware**: `dev08` (`192.168.0.228`, dual-socket Intel Xeon Gold 6230, 40 physical cores, 187 GB RAM; pinned CPU sets, release + `mpi+mkl` builds)

---

## 1. Performance Objectives

Per Research Lead directive (`D-006`): do not merely guard against pure-PW regression — actively engineer high performance into the solvated SCF path.

| Benchmark Class | Metric | Pre-Declared Target |
| --- | --- | --- |
| **Pure-PW Off-Path (`solvation = None`)** | Wall time vs `main` baseline (Si8, MgO, k444s) | `<= 1.01×` (within noise) & **bit-identical** energies/forces |
| **Pure-PW Hot Kernels** | Allocations / branches / collectives added to `apply_block`, `ref_davidson*` | **Strictly `0`** (audited by diff inspection) |
| **Solvated SCF Overhead (`DirectPcg`, adaptive + warm-start)** | `(T_solv - T_vac) / T_vac` on 64³–96³ molecular/slab cells (`P = 4–8` ranks) | **`<= 15%`** (target `< 10%` with `n_g` reuse + adaptive tolerance ratchet) |
| **Solvated SCF Overhead (`FisicaroPcg`, `dx <= 0.25 Bohr`, adaptive + warm-start)** | `(T_solv - T_vac) / T_vac` on 64³–96³ cells (`P = 4–8` ranks) | **`<= 5%`** |
| **Per-Step Heap Allocations in Inner PDE Loop** | Allocations during PCG iterations | **`0`** |
| **MPI Strong Scaling (Solvation Kernel)** | Wall time at $P = 1, 2, 4, 8$ ranks (`48^3`, `64^3`, `96^3`) | Monotonic speedup with `|E(P) - E(1)| < 1e-12 Ha` (`DirectPcg`) |

---

## 2. Benchmark Matrix on `dev08`

### Suite A — Pure-PW Zero-Regression Check (`solvation = None`)
Run existing BAANDI reference cases (`si8` Gamma, `mgo` k-point, `k444s`) on `main` vs `research/track-b-gate4-scf` (3 repetitions each, median wall time, `mpi+mkl` release binary):
- Assert bit-identical `total_energy_ha`, `scf_iterations`, and force components.
- Record wall-time ratio `T_branch / T_main` and stage breakdown (`PhaseTimings`).

### Suite B — Solvated SCF Wall-Time & Iteration Decomposition (`solvation = Some(...)`)
Three representative DFT workloads (real external PseudoDojo NC UPFs on `dev08`, $\epsilon_b = 78.4$, VASPsol-erfc & SCCS):
1. **Small polar molecule**: $\text{H}_2\text{O}$ in `16 Bohr` box ($E_{\mathrm{cut}} = 60\text{ Ry}$, ~`64^3` charge grid), $P \in \{1, 4\}$ ranks.
2. **Larger polar/ionic system**: $\text{H}_2\text{O}$ / $\text{OH}^-$ in `20 Bohr` box ($E_{\mathrm{cut}} = 60\text{ Ry}$, ~`80^3` charge grid), $\kappa \in \{0, 0.15\}\text{ Bohr}^{-1}$, $P \in \{1, 4, 8\}$ ranks.
3. **Slab with vacuum/electrolyte**: Metallic or semiconductor slab (`z = 24–30 Bohr`, Fermi–Dirac smearing, $\kappa = 0.15$), $P \in \{4, 8\}$ ranks, verifying Broyden stability and adaptive tolerance behavior.

For every run, record from `PhaseTimings` + `SolvationEnergyLedger`:
- Total SCF iterations (vacuum vs solvated)
- Total wall time (vacuum vs solvated) and exact solvation fraction (`T_solv_kernel / T_total`)
- Cumulative PDE iterations & 3D FFTs across the entire SCF cycle:
  - Cold-start fixed `1e-10` vs warm-start fixed `1e-10` vs **warm-start + adaptive tolerance ratchet** (`next_solvation_tol`)
- `DirectPcg` vs `FisicaroPcg` wall time, energy difference `|E_fis - E_direct|`, and force difference `max |F_fis - F_direct|`.

### Suite C — Gate 4 Solver Acceleration R&D (Closing the 5× Exact-Variational Gap)
Gate 2 revealed a clean engineering opportunity: `DirectPcg` is exact-variational (`< 1e-9`) and strictly SPD on every grid, while `FisicaroPcg` is 3.4×–5.3× faster on sharp `64^3` cavities (`7.2 µHa` aliasing diff). In Gate 4, benchmark two targeted upgrades on `solvation_bench` (`32^3, 48^3, 64^3, 96^3`):
1. **Consistent Fisicaro RHS** ($b_{\mathrm{pol,Fis}} = 4\pi\rho_{\mathrm{sol}} - A_{\mathrm{Fis}}\phi_{\mathrm{vac}}$): measure whether eliminating the spurious $(A_{\mathrm{Fis}} - A_{\mathrm{disc}})\phi_{\mathrm{vac}}$ source term cuts `FisicaroPcg`'s coarse-grid (`32^3`/`48^3`) energy discrepancy by 10×–100×.
2. **Inexact Fisicaro Preconditioner inside Flexible PCG (FPCG / Polak–Ribière)** targeting the **exact discrete operator $A_{\mathrm{disc}}$**: run 4–8 inner Fisicaro steps as a preconditioner to outer $A_{\mathrm{disc}}$ iterations and test whether total wall time beats standalone `DirectPcg` by `>= 2×` while keeping exact `< 1e-9` discrete variational consistency!
   - *Adopt-what-improves rule*: retain only variants that measurably beat `DirectPcg` on `dev08` release benchmarks without compromising variational accuracy or stability; document rejected variants in `decision-log.md`.

---

## 3. Reproducibility & Reporting Standard

Every recorded benchmark table in the Gate 4 report must include:
- Exact git commit SHA, Rust compiler & MKL version (`dev08`), CPU pinning (`taskset`), rank/thread counts (`P`, `BAANDI_THREADS`), and exact CLI/example command.
- Raw medians across `>= 3` runs (never cherry-picked single runs) plus negative results if any solver/preconditioner arm fails to beat baseline.
