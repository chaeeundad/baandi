# Decision table: ideas and their status

Status words: ADOPTED (default), OPT-IN (in main, off by default), REJECTED, OPEN. Evidence ids refer to rows of docs/performance/attempts-ledger.md or ADRs.

| Idea | Status | Why | Evidence |
|---|---|---|---|
| Real Gamma path (half sphere, paired FFT, real Davidson) | ADOPTED default auto | 1.1 to 3.2x over complex path | ADR 0040, GR-AUTO |
| Atomic start wavefunctions | ADOPTED default | si32 strict 40.6 to 19.9 s at T=1, SCF spread to 0 | ADR 0042, ATOMIC-DEFAULT, R3 |
| Atomic charge superposition start density | OPT-IN | no gain in SCF count on si8/si16 | QE-START |
| MKL dgemm for real Davidson products (R1, R2) | ADOPTED on MKL builds | -8 to -27 percent combined | MKL-R1, MKL-R2 |
| Blocked Gram-Schmidt | ADOPTED on MKL builds only | serial MGS is a T=8 bottleneck; pure builds lose | MGS-B, R2 |
| MKL dsyevd reduced eigensolver | ADOPTED on MKL builds | 1.6 to 1.9x faster than faer, -5.7 percent at T=8 si32 | MKL-EVD |
| Remove redundant dgemm input finiteness scans | ADOPTED | bit-identical, -2 to -16 percent | FIN1 |
| Band-parallel nonlocal energy in the audit | ADOPTED | total_energy -74 percent at T=8 si32, bit-identical | AUD1 |
| Radial transform caches in stress and forces | ADOPTED | forces 50x, stress frozen 54x | STR1 |
| Pulay scan in standard profile (stress) | ADOPTED | 2.4 to 2.8x on the scan, estimate moves 1e-11 | STR2 |
| Optional blocked rotation kernel (pure Rust) | REJECTED | auto-vectorised plain Rust loses to faer/gemm micro-kernels | ROT1 |
| MKL-style panels of other sizes, adaptive TN panel | REJECTED | trajectory changes dominate; dgemm GFLOP/s flat | PANEL-INSITU |
| Partial rotation (rotate unconverged roots only) | REJECTED | converged roots do not stay converged | PRT1 |
| Davidson preconditioners TPA, Ritz shift | REJECTED | Si8 standard SCF count rises | P47A |
| fp32 H application early in the SCF | REJECTED | upper bound 6 to 7 percent | P47B |
| Locking, generalized Davidson, band-count rules | REJECTED | accuracy gates or extra H work | P24A, P37A, P38A, P42A, P43A |
| Mixing alternatives | REJECTED as default | case dependent | P40A |
| Cutoff-scan warm start for stress | REJECTED | -4 percent only | STR2 |
| Panel-size change for threads | REJECTED | panel count was never the thread limit | PANEL-INSITU |
| Reference compat mode (PBE) | ADOPTED option | closes PseudoDojo PBE to 1e-8 Ry | ADR 0038 |
| `local_mesh = qe_like_g0` | OPT-IN | closes the Si.pz-vbc LDA offset (head interval of the log mesh) | ADR 0037 addendum, OFFSET-G0 |
| MPI layer | OPEN (approved, staged) | gap partly parallel efficiency, multi-node needs it | ADR 0041 |
| Nonlocal operator as GEMM (MKL path) | OPEN | 6 to 12x QE's calbec/add_vuspsi; earlier traversal rewrites lost | stage-decomposition-results |
| Batching single-column projection products | OPEN | 25 percent of dgemm calls have n = 1 | DGEMM-CENSUS |
| Start from the density of a similar system | PARKED by the user | no design | wiki state-of-play |

- 2026-11-10 BAANDI_MKL_NONLOCAL (MKL dgemm for nonlocal projector products): -6.4 percent strict and -4.2 percent standard wall at T=1, only -2.0 and -1.5 percent at T=8; NOT adopted as default (pre-registered 3 percent in both thread budgets); opt-in. Ledger NONLOCAL-DGEMM-AB.
- 2026-11-10 BAANDI_CERT_SCREEN=on (skip fixed-point certificate re-solves that a first-order eigenvalue-shift estimate shows cannot pass; trajectory-neutral): every si32 seed 601-610 fails its first two certificates on the eigenvalue criterion (26 percent of H vectors); Mac -17 to -21 percent wall, density SHA identical; OPT-IN until the pre-registered dev08 A/B (docs/performance/requests/single-node-speed.md). Ledger CERT-SCREEN.
- 2026-11-10 BAANDI_NL_ENERGY=fast (nonlocal energy of the audit from projections only, occupied bands only): energy phase -55 to -60 percent on the Mac, energy within 3e-13 Ha; OPT-IN until the dev08 A/B. Ledger NL-ENERGY. Density pairing (DENSITY-PAIR) and potential reuse (POTENTIAL-REUSE) are recorded as below-2-percent ideas, not implemented.
- 2026-11-10 (proposal B3) POTENTIAL-REUSE implemented with the screen (bit-identical; `BAANDI_POTENTIAL_REUSE=off` for A/B) and the screen's cluster transforms made pool-parallel (bit-identical for any thread count): Mac si16 T=8 screen 34 -> 16 ms per certificate; the potential build was only 2.5 ms of it (the proposal's premise was wrong, the band transforms were the cost). Second registration with fresh seeds 701-710 for CERT_SCREEN + NL_ENERGY=fast as defaults written: docs/performance/requests/default-candidates.md. USER DECISION after dev08 data.
- 2026-11-10 (proposal B2) `BAANDI_ENERGY_IDENTITY=on`: non-final iterations take the nonlocal energy from the band-sum identity (history energy = band-structure form), explicit audit at the certificate and the final state, final energy bit-identical to off on the Mac. Changes the meaning of history energies: USER DECISION before any default. Ledger ENERGY-IDENTITY.
- 2026-11-10 (proposal B6) physics-side cached-projection API only (docs/performance/requests/b6.md); solver half deferred to Track A. C3 and C1 are design notes (docs/roadmap/c3-ionic-extrapolation.md, c1-cutoff-ladder.md).
- 2026-11-10 Davidson stage timers `BAANDI_DAVIDSON_STAGES=on` (measurement only): on the Mac a width-1 iteration of si32 strict costs about 20 H columns; rotate 40, orthonormalize 27, H 24, eigen 6 percent of Davidson time at T=1 (eigen 14 at T=8). Ledger DAVIDSON-STAGES.
- 2026-11-10 B1 selected-root products `BAANDI_B1=standard|strict|on` (+ `BAANDI_B1_SKETCH=count`): decisions and H work unchanged, energies within 3.4e-13 Ha; Mac: dense sketch roughly neutral (its maintenance cancels the rotation saving), count sketch positive; OPT-IN until the dev08 A/B (requests/b1-b4-b5-davidson.md). Ledger B1-SELECTED-ROOTS (premise change against PRT1).
- 2026-11-10 B4 inner iterations for narrow blocks: REJECTED on the Mac harness (no fixture with a single-column tail; -3 to -9 percent steps where it acts). Ledger B4-INNER.
- 2026-11-10 B5 bordered reduced eigen-update `BAANDI_B5=on`: gate passed (eigen above 5 percent on strict cells), update 0.5 to 0.7 of a full solve, Mac effect inside noise; OPT-IN until the dev08 A/B. Ledger B5-BORDERED.
