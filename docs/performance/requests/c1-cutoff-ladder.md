# dev08 request: C1 cutoff ladder (`BAANDI_CUTOFF_LADDER`), 2026-11-10

Status: implemented on branch `worktree-agent-a719eab3f9d2f2305` (commit in the hand-off report), default OFF. Design: `docs/roadmap/c1-cutoff-ladder.md`. Per the adoption rule (AGENTS.md, 2026-11-10) there is no percentage threshold: the ladder is adopted as the default when the median wall improves and nothing regresses beyond run-to-run noise, and it stays switchable by the environment variable.

## What the switch does
`BAANDI_CUTOFF_LADDER=on` (or `[scf] cutoff_ladder = "on"`; `BAANDI_CUTOFF_LADDER=off` forces it off; a number in (0, 1) sets the fraction) runs `dft-cli run` as:
1. a coarse SCF at E1 = 0.5 x ecutwfc on the FFT grid E1 implies (a fixed `fft_grid` is dropped for this stage only), ordinary session, stopped when its density residual is at or below 1e-4 Bohr^-3/2 or after 30 coarse iterations (no certificate at E1);
2. exact prolongation: orbital coefficients copied by Miller index into the full basis (zero for new G; an injection, orthonormality and Gamma Hermitian symmetry preserved bit for bit), density by Fourier zero padding (exact for the band-limited coarse field, charge unchanged);
3. the unchanged strict session at the full cutoff, started from that density and those warm orbitals, with an empty mixer. Convergence, the energy audit, eigenpairs and the fixed-point certificate are all at the full cutoff.
The rule is general: one stated fraction of the cutoff (the design's pre-registered arm 0.5), one fixed hand-off residual (pre-registered 1e-4), one budget bound. Nothing is tuned per system.
Eligible: one exact Gamma point on the real path (`gamma_real` not off), fixed uniform integer occupations, strict profile, no symmetry projector, no strain reference, Davidson without locking/guard. Everything else (k meshes, standard profile, Fermi-Dirac, spin, symmetry, observe/agent modes) runs the ordinary SCF; the reason is in `results.json` `cutoff_ladder.fallback_reason` (the block is present whenever the switch is on).

## Mac observations (release, T=1, pure Rust backend, one run each, not acceptance)
Perturbed atoms (fractional +-0.01 uniform, seeds as listed), strict (energy tolerance 1e-10 Ha), Si.pz-vbc 30 Ry unless noted; synthetic = Gaussian ion with s,s,p projectors (non-physical).

| case | plain SCF / H cols / wall s | ladder coarse SCF (H full-grid eq.) + fine SCF / H cols / wall s | dE (Ha) |
| --- | --- | --- | --- |
| si8 unperturbed | 11 / 1319 / 0.688 | 6 (184) + 10 / 1092 / 0.626 | 1e-13 |
| si8 seed 601 | 19 / 2163 / 0.966 | 9 (315) + 12 / 1316 / 0.734 | 9e-14 |
| si8 seed 602 | 17 / 1897 / 0.857 | 8 (277) + 11 / 1207 / 0.662 | 4e-13 |
| si16 seed 601 | 20 / 4434 / 5.056 | 10 (720) + 12 / 2607 / 3.768 | 8e-13 |
| si16 seed 602 | 23 / 5047 / 5.435 | 10 (727) + 16 / 3421 / 4.473 | 1.4e-13 |
| si16 seed 603 | 22 / 4742 / 5.345 | 12 (842) + 14 / 3056 / 4.375 | 4e-13 |
| si32 seed 601 | 24 / 11122 / 39.73 | 9 (1388) + 15 / 6747 / 28.69 | 2.7e-13 |
| si primitive (Gamma) | 12 / 357 / 0.057 | 7 (74) + 8 / 231 / 0.053 | 8e-14 |
| si-skew (sheared) | 21 / 614 / 0.194 | 12 (110) + 11 / 318 / 0.109 | <1e-13 |
| syn fcc 8 Ry | 20 / 615 / 0.020 | 18 (66) + 12 / 378 / 0.016 | 2e-15 |
| syn fcc seed 601 | 16 / 522 / 0.016 | 11 (47) + 12 / 381 / 0.016 | 9e-15 |
| syn skewed 10 Ry | 12 / 433 / 0.024 | 6 (56) + 12 / 411 / 0.034 | <1e-13 |
| syn slab-like 10 Ry | 18 / 651 / 0.046 | 17 (129) + 12 / 431 / 0.031 | 1e-13 |

H cols = recorded plus certificate H vectors at the full cutoff; the coarse H count is scaled by the FFT grid-point ratio (e.g. si32 54x54x25 / 72x72x36 = 0.42). Every run converged with a passed full-cutoff certificate. The only worse wall is a 24 ms synthetic cell, where the coarse setup dominates. Fraction sweep (reported only): 0.6 and 0.7 were not better than 0.5 on si16 (4.09 s and 5.68 s against 3.77 s) and mixed on si8.

## Protocol
One MKL build of the branch HEAD (`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`), the same binary for both arms; the A/B is the environment only (`BAANDI_CUTOFF_LADDER=off` against `on`). Cells `benchmarks/cases/si16-gamma` and `si32-gamma`, perturbed seeds 601-610 (as in MPI-M1-GATE3), strict and standard profiles, T=1 (taskset core 0) and T=8 (cores 0-7). Interleave off/on per seed; node idle (load below 0.30); uptime before and after.
Standard profile: the ladder is ineligible (the real path excludes it), so the standard cells measure only that the fallback costs nothing (expected identical results, `cutoff_ladder.used = false`, reason `standard profile`).

## Correctness clauses (any violation rejects)
1. Every ladder run converges with a passed full-cutoff certificate (`convergence.converged = true`).
2. Accuracy class: |E(on) - E(off)| below 1e-10 Ha (strict) and 5e-9 Ha (standard) for every seed and T; occupied eigenvalues within the certificate's 1e-7 Ha. Density SHA is NOT required to match between arms (different start), but must be identical between T=1 and T=8 within each arm (thread-count independence; the ladder's own reductions are serial and fixed order).
3. Standard cells: energy and density SHA bit-identical between arms (fallback path).

## Adoption
Adopted as the default (then `BAANDI_CUTOFF_LADDER=off` is the rollback) when the median wall over seeds 601-610 improves on si16 and si32 strict at T=1 and T=8 and no cell regresses beyond the run-to-run noise of the off arm.

## What to send back
Per run: wall, `convergence.iterations`, `convergence.work` (H vectors recorded and certificate), the `cutoff_ladder` block (coarse iterations, coarse H vectors, full-grid equivalents, coarse stage seconds), energy (17 digits), density SHA, exit code; binary sha256, commit, build command, uptime before and after, interleave order.
