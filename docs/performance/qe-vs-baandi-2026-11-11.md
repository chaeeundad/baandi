# QE 7.6 against BAANDI after the 2026-11 adoption round (like-for-like, dev08)

Status 2026-11-11. Measured by the master session on dev08 at origin/main 208b226 (binaries cli-mkl-208b226, cli-mpi-208b226, pw.x orig-mpi-dfti): same day, same node, load gated, every run of both codes converged. QE: `K_POINTS gamma`, nosym/noinv, `startingwfc = atomic`, `diago_full_acc`, same FFT grid and `nbnd`, OMP and MKL threads 1, conv_thr 1e-8 (standard) and 1e-10 (strict). BAANDI: current defaults (MKL build). "unp" = unperturbed single start (median of 3); "pert" = the same displaced coordinates (+-1e-7 Bohr, seeds 601-605, median) given to BOTH codes.

## One core (QE one process, BAANDI T=1), wall seconds
| cell | QE unp | BAANDI unp | ratio | QE pert | BAANDI pert | ratio |
|---|---|---|---|---|---|---|
| si8 standard | 0.464 | 0.520 | 1.12x | 0.377 | 0.522 | 1.38x |
| si8 strict | 0.526 | 0.986 | 1.87x | 0.464 | 1.180 | 2.54x |
| si16 standard | 1.079 | 1.737 | 1.61x | 1.112 | 1.899 | 1.71x |
| si16 strict | 1.324 | 3.124 | 2.36x | 1.562 | 3.772 | 2.41x |
| si32 standard | 2.993 | 7.252 | 2.42x | 3.050 | 8.210 | 2.69x |
| si32 strict | 3.512 | 13.337 | 3.80x | 3.581 | 21.759 | 6.08x |

## Eight cores, mixed models (QE 8 MPI ranks; BAANDI 8 threads; BAANDI MPI M2 np=8), unperturbed
| cell | QE np=8 | BAANDI threads | ratio | BAANDI MPI np=8 | ratio |
|---|---|---|---|---|---|
| si8 standard | 0.218 | 0.447 | 2.05x | 0.577 | 2.64x |
| si8 strict | 0.257 | 0.871 | 3.39x | 0.878 | 3.42x |
| si16 standard | 0.400 | 1.226 | 3.06x | 1.536 | 3.84x |
| si16 strict | 0.460 | 2.053 | 4.47x | 2.204 | 4.79x |
| si32 standard | 0.717 | 4.115 | 5.74x | 5.955 | 8.30x |
| si32 strict | 0.825 | 6.720 | 8.15x | 8.399 | 10.19x |

## Work counts (unperturbed)
BAANDI SCF / Davidson / H vectors: si8 standard 9/28/495, si8 strict 9/174/1054, si16 standard 9/28/873, si16 strict 8/72/1553, si32 standard 10/31/1843, si32 strict 8/64/3060. QE SCF / h_psi calls (standard, strict): si8 13/82 and 13/82, si16 8/86 and 12/146, si32 7/28 and 8/35. QE needs more SCF steps on si8 and far fewer H applications on the large cells.

## Energy
The known offset reproduces exactly: BAANDI minus QE is +1.273e-6 Ha per atom in every cell and class (+1.018e-5 Ha si8, +2.036e-5 si16, +4.07e-5 si32).

## Findings
1. One core: BAANDI is 1.1x to 3.8x slower on the unperturbed starts and 1.4x to 6.1x slower on the perturbed starts. Before the round (qe-vs-baandi-current.md) the unperturbed gap was 2.1x to 6.1x.
2. BAANDI is far more sensitive to breaking the exact symmetry of the cell than QE: a +-1e-7 Bohr displacement costs BAANDI +0.4, +19.7, +9.3, +20.7, +13.2 and +63.1 percent (si8 standard to si32 strict) and costs QE -18.8, -11.8, +3.1, +18.0, +1.9 and +2.0 percent on identical geometry. A perfectly symmetric cell is the easiest case; real structures are never exactly symmetric, so the PERTURBED column is the honest one (si32 strict on one core: 6.1x, not 3.8x).
3. At eight cores the gap widens (up to 8.2x threads, 10.2x MPI on si32 strict, unperturbed) because QE scales better with MPI at these sizes; BAANDI's MPI M2 is at parity with its own threads on one socket (ledger MPI-M2D-AB) but starts from a slower one-core base on the large strict case.
4. Lead: the symmetry-breaking sensitivity of the strict class is a larger single effect than anything in the micro-stack; hypothesis under test (switch BAANDI_EMPTY_LOOSE): the lone unconverged first-empty band (strict threshold, inside a near-degenerate cluster once the symmetry is broken) drives the single-column tail.
