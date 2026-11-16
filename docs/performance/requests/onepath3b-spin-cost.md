# dev08 request: ONEPATH-3B zero-moment spin cost (atomic start for collinear fixed sectors)

Status: pre-registered (2026-11-16), not yet run. Ticket `tickets/ONEPATH-3B-SPIN-COST.md`,
ADR 0042 amendment ONEPATH-3B. BAANDI-only: no QE run, no QE ratio; every table is one
model at one thread count (`comparison-rules.md` satisfied trivially). Accuracy classes
unchanged (strict 1e-10 Ha, standard 5e-9 Ha); no tolerance changed.

## Finding being fixed
dev08, binary 006baaf, 3 repeats identical: Si8 Gamma PBE (PseudoDojo NC Si.upf),
standard path, plain input. Unpolarized 9 SCF / 1 verification attempt / 2.67 s,
E = -33.472972013925983 Ha. Zero-moment collinear (`total_moment = 0.0`) 22 / 2 / 6.33 s,
E = -33.47297201386214 Ha. Cause (ticket): the spin driver did not use the atomic start
(ADR 0042 exclusion), the unpolarized run did. The fix gives cold fixed-sector solves the
closed-shell atomic start (span once per k, Rayleigh-Ritz per channel). Fermi-Dirac
problems keep the free-electron start.

## Build
Branch HEAD of ONEPATH-3B after merge with main. MKL build (`--features mkl`,
`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`) and the pure build. Binary
sha256, build command and exit code recorded. Same binary for every arm of a build.

## Arms (T=1, taskset core 0; node idle, load below 0.30; uptime before and after; arms interleaved, 3 repeats)
- U0: si8-gamma (`benchmarks/cases/si8-gamma.toml`), `xc.functional = "pbe"`, PseudoDojo
  NC PBE Si, 30 Ry, plain input (unpolarized; record whether the real Gamma path and the
  cutoff ladder ran).
- U0c: U0 with `[scf] gamma_real = "off"` (the complex unpolarized path the spin run uses).
- S0: U0 with `[spin] mode = "collinear"`, `total_moment = 0.0`.
- S0fe: S0 with `[scf] start_wavefunctions = "free_electron"` (the pre-fix start).
- F1: ONEPATH-3 fixture F1 (Fe fixed sector) on the new binary and on 006baaf (changes:
  it now uses the atomic start when the Fe file has PP_PSWFC).
- F2, F3, M1, M2: ONEPATH-3 fixtures (Fe FD free / fixed moment, Al FD) on the new binary
  and on 006baaf (expected unchanged: Fermi-Dirac keeps the free-electron start).

## Pre-registered checks (bounds fixed now, not widened)
- S0 vs U0c: SCF counts equal or within 1, equal verification attempt counts,
  |dE| <= 1e-10 Ha. U0 and U0c energies bit-identical to 006baaf (unpolarized code path
  untouched).
- S0 wall at most 2.5x U0c (two channel solves per step; Mac observation 2.1x).
- S0fe reproduces the 006baaf S0 row (same SCF count and attempts, |dE| <= 1e-12 Ha).
- F2, F3, M1, M2: energy bits, SCF counts and density SHA-256 equal to 006baaf.
- F1: converged and verified, |dE| <= 5e-9 Ha against the ONEPATH-3 TIGHT row; SCF count
  and wall reported against 006baaf (information).

## Mac observations already made (not acceptance)
Mac release, pure, T=1, same Si.upf (sha256 39822757...5282): before S0 18 SCF / 3
attempts / 3.20 s; after S0 9 / 1 / 1.58 s, |dE vs U0c| 5.4e-13 Ha (U0c 9 / 1 / 0.74 s);
S0fe 18 / 3 / 3.20 s; U0 6 / 2 / 0.47 s (real Gamma + ladder). Synthetic dimer: SCF counts
equal for both starts; atomic-start |dE| 2.6e-10 to 5.8e-10 Ha inside the closed-shell
atomic-start rounding envelope (up to 1.5e-9 Ha), see the ticket.

## What to send back
Per run: arm, build, exit code, wall, `timings`, `convergence` (iterations, history),
`verification`, `atomic_start`, `execution`, energy (17 digits), density SHA-256; binary
sha256, commit, NC file SHA-256s, uptime before and after, interleave order.
