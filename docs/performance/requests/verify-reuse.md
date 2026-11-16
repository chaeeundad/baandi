# dev08 request: fast-loop stage 2 for VERIFY-REUSE (and C3a), threads / one core

Status: pre-registered 2026-11-16, not yet run. Stage 2 of `docs/performance/fast-loop-protocol.md` (five to six cells, paired env A/B, same binary). Screening numbers (Mac, O) are in the ledger rows VERIFY-REUSE-FEED, VERIFY-REUSE-WARM, C3a-SCREEN and in `docs/performance/profile-2026-11-16-threads.md`. BAANDI-only: no QE column, no QE ratio; every table is one model at one thread count (`comparison-rules.md`). Accuracy classes unchanged (strict 1e-10 Ha, standard 5e-9 Ha); no tolerance, threshold or stop rule changed; acceptance still needs a passed verification.

## What changes
`BAANDI_VERIFY_REUSE=off|warm|feed` (default off; `crates/dft-engine/src/session.rs`). Acts only after a FAILED fixed-point verification (ADR 0043 V1/V2); with no failure the run is bit-identical to off.
- warm: the next SCF step starts from the verification solve's orbitals (eigenpairs at V(rho_out)); the REUSE-PSI kept start is then taken instead of declined.
- feed: warm, and the pair (rho_out, F(rho_out)) enters the mixer after the step's own pair (rho_in, rho_out); the next input density is the mixer's proposal from that history. The verification map becomes an SCF map instead of being discarded.
Rounding class (the trajectory after the first failed attempt changes). Thread-count independence: no new parallel code; density SHA must stay equal across T within each arm.

## Build
Branch HEAD (commit given by the orchestrator) after merge with main; MKL build (`--features mkl`, `MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`) and pure build, separate target directories, binary sha256 and build exit codes recorded. One binary per build for every arm.

## Arms (same binary, environment only, interleaved per geometry)
- A: default (switch unset).
- F: `BAANDI_VERIFY_REUSE=feed`.
- W: `BAANDI_VERIFY_REUSE=warm` (si32 only; cost attribution).
- P: `BAANDI_LEAN_PROJ_CACHE=on` (C3a, independent code path; stage-2 data for the existing C3a request).
- FP: F + P (stacking check).

## Cells and geometries (stage 2: not the full generality set)
f1e-3 fractional displacement, p48lib generator, plain single-path inputs:
- si16, si32: seeds 601-610 (the Mac shows failures on about half of the si32 seeds and a near-threshold loss on si16 seeds 601/602/605, so ten seeds per cell are needed).
- si64: seeds 601-603.
- si-k444s: seeds 601-603.
- mgo-40 (default mixer of the cell): seeds 601-603.
- si-skew: seeds 601-603.
T=1 (taskset core 0) for all; T=8 (cores 0-7) for si32 and si64 arms A and F. Node idle (load below 0.30 before each run), uptime before and after, 3 repeats per arm and geometry for si16/si32 T=1 (1 repeat elsewhere; counters are exact), INVALID guard (exit code, empty output, wall below 0.1 s).

## Decisive clauses (fixed now)
1. Identity: on every geometry whose arm-A run has exactly one verification attempt, arm F has the same energy bits, density SHA and counts as arm A (the switch must not act without a failure).
2. Convergence: every F run ends `converged = true`, `verified = true` within the default budget wherever A does.
3. Accuracy: |E(F) - E(REF)| <= 5e-9 Ha where REF = `BAANDI_REFERENCE_STRICT=on` (one REF run per geometry, si16/si32/si64/mgo-40).
4. Thread independence: density SHA of F at T=1 equals T=8 (si32, si64).
5. Speed, by the adoption rule, med_by_seed per cell: F adopted as default only if the median over seeds of the per-seed median wall improves in the cells where failures occur and no cell regresses beyond its run-to-run spread. Report per seed: verification attempts A/F, fine SCF steps, H vectors (recorded + certificate), wall.
A failure of 1 rejects the change; 2 to 4 keep it off and go to the ledger with numbers; 5 decides default on/off (stage 3 then runs the full set once).

## Question for the user (not a switch; reported, no code)
The cost of a PASSED verification (si32 about 8 %, si64 16 % relayed) is a warm solve at tau_v = ratchet(tau, E, N_e), fixed by ADR 0043 V1 (never looser than the last step's tau). For si32 seed 601 tau_v = 8.4e-14 Ha against the last step's 5.9e-12 Ha. A looser tau_v (for example tau_v = tau) would shorten every verification but is an ADR change; it is not proposed here as a performance switch. Please record, per run, `verification.attempts[].eigen_tolerance_ha`, `davidson_iterations`, `h_vectors`, `wall_s`, so the user can decide with numbers.

## What to send back
Per run: arm, cell, seed, T, exit code, wall, `timings` (with `verification_s`), `converged`, `verified`, verification attempts block, fine and coarse SCF steps, `convergence.work` (h_vectors_recorded, certificate_h_vectors, davidson iterations, gemm_dispatch rows, davidson_lean.proj_cache), total energy (17 digits), density SHA-256; binary sha256, commit, build command, uptime before and after, interleave order; result files stored on dev08 with sha256.
