# dev08 request: P0.1 single SCF path with fixed-point verification (ADR 0043)

Status: pre-registered, not yet run. Nothing here changes a default. `BAANDI_SINGLE_PATH` stays `off` until this check passes and the user approves the flip (ADR 0043 amendment). BAANDI-only: no QE run and no QE ratio is requested, so the parallel-model rule (`comparison-rules.md`) is satisfied trivially; every table below is one model at one thread count.

## What is being checked
The work order (P0.1) makes the standard path the only SCF path and runs the certificate once after convergence as a verification (ADR 0043: pass when the Coulomb fixed-point residual of the delivered density, E_H[F(rho_out) - rho_out], is below the threshold; a failure continues SCF and verifies again). The standard path has never been measured at the tight threshold, and forces have never been computed from a standard state.

## Build
- Candidate commit (branch HEAD after merge) and current `main` for the identity clause. MKL build (`--features mkl`, `MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`) for everything; pure build for the identity clause only. Separate target directories, binary sha256, build exit codes.

## Thresholds
QE `conv_thr` is in Ry; BAANDI `standard_threshold_ha` in Ha (same estimator definition, ADR 0035). The two classes of the work order:
- "1e-8": `standard_threshold_ha = 5e-9` (today's standard default; matches QE `conv_thr = 1e-8` Ry).
- "1e-10": `standard_threshold_ha = 5e-11` (matches QE `conv_thr = 1e-10` Ry; never measured on the standard path).

## Arms (same binary, same coordinates per geometry)
- REF: explicit `convergence_profile = "strict"` (switch unset). The strict reference of each geometry.
- S8: `BAANDI_SINGLE_PATH=on`, plain input (no profile), `standard_threshold_ha = 5e-9`, verification default on.
- S10: as S8 with `standard_threshold_ha = 5e-11`.
- S8-NV (cost attribution, si32 and si64 only): S8 with `scf.verify = "off"`.
Defaults as on main otherwise (cutoff ladder on, atomic start, gamma_real auto).

## Cells and geometries
- si16, si32, si64 Gamma (`benchmarks/cases/si{16,32,64}-gamma.toml`), fractional displacement +-1e-3 (headline amplitude) from the MPI-M1-GATE3 generator (`random.Random(seed)`, same stream as `qe-remeasure-2026-11-13.md`), fresh seeds 901-905 (si64: 901-903). Also unperturbed (one geometry each).
- Generality set (unperturbed and seeds 901-903 at fractional +-1e-3): k444s, mgo-40, syn-p3, si-skew, syn-slab, si-elongated. si-elongated and syn-slab hit the iteration limit on older strict defaults (P10G): record whatever happens on REF and on S8/S10, a stall on both arms is a finding, not a failure of the single path.
- Threads: T=1 (taskset core 0) for every run; T=8 (cores 0-7) for si32 and si64 S8, S10 and REF. Node idle (load below 0.30), uptime before and after, arms interleaved per geometry.

## Decisive clauses (fixed now)
1. Identity of the switch: `BAANDI_SINGLE_PATH` unset and `=off` on the candidate give density SHA-256, total energy bits and SCF counts identical to `main` on gamma, k222s, k444s, si8, si16 (default strict input, and an explicit `convergence_profile = "standard"` input), pure and MKL, T=1 and T=8. Any difference rejects.
2. Convergence: every S8 and S10 run ends `convergence.converged = true` and `verified = true` within the default iteration budget, on every geometry where REF converges.
3. Delivered error S8: |E(S8) - E(REF)| <= 5e-9 Ha (standard class) on every geometry where both converge.
4. No stall: S8 and S10 SCF iteration count (full-cutoff history) at most 2 x REF's plus 3 on every geometry.
5. Thread independence: S8 density SHA-256 equal at T=1 and T=8 (si32, si64).
A failure of 1 rejects the change; a failure of 2 to 5 keeps the switch default off and goes to the ledger with the numbers.

## Reported, not decisive (for the user's decisions)
- S10 delivered error |E(S10) - E(REF)| against the strict class 1e-10 Ha, and the S10 iteration count and wall against REF (does the tight standard path replace strict in speed and accuracy?).
- Delivered density error: L2 norm of rho(S) - rho(REF) from `density.bin` (Bohr^-3/2, same normalisation as `density_residual`), and the largest occupied-eigenvalue difference (Ha).
- Verification: pass rate on the first attempt, number of attempts per run, `verification.attempts[].coulomb_error_ha / threshold_ha`, `strict_criteria_met` rate, `timings.verification_s` and its share of `timings.total_s`, its H vectors (`verification.attempts[].h_vectors`). Cost attribution: wall(S8) - wall(S8-NV) against `timings.verification_s` (si32, si64, medians).
- Work: total H vectors (recorded plus certificate) and wall for REF, S8, S10 per cell; medians over seeds with spread.

## Forces (reported to the user; no pass/fail)
- Geometries: si16 and si32 Gamma at fractional +-1e-3 (seeds 901-903), and mgo8 unperturbed plus seed 901.
- Arms: REF with `output.forces = true` (strict certificate), S8 and S10 with `output.forces = true` (forces require the passed verification).
- Report per geometry: maximum absolute component difference |F(S) - F(REF)| (Ha/Bohr) and the RMS over atoms, the largest REF component for scale, and the S-arm `verification` block. The user decides from this which threshold relaxations and MD use (work order P0.1, last bullet). Relaxation and stress stay on the strict path in this stage (ADR 0043 item 5).

## Mac observations already made (not acceptance)
Synthetic fixtures only (debug build): engine unit fixture (one Xx atom, skewed 4 Bohr cell) single path 5e-9: one attempt, Coulomb residual 8.9e-11 Ha, 6 H vectors, max |dF| against strict 6.7e-7 Ha/Bohr; CLI fixture (two H atoms, 8 Bohr cube): one attempt, verification 1.9 ms, max |dF| 1.5e-7 Ha/Bohr, dE 2.1e-11 Ha. A forced failed verification continues SCF and passes on the next stop (unit test).

## What to send back
Per run: arm, cell, seed, T, exit code, wall, `timings` (including `verification_s` and `certificate`), `convergence.converged`, iterations, `convergence.work`, `verification` block, `single_path` block, total energy (17 digits), density SHA-256, forces block where requested; binary sha256, commit, build command, uptime before and after, interleave order. Original result files stored on dev08 with sha256.
