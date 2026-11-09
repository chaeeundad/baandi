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
