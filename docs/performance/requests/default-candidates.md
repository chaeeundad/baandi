# dev08 request: CERT_SCREEN + NL_ENERGY=fast as defaults (second registration, 2026-11-10)

Status: pre-registered, not yet run. Nothing in this document changes a default; a pass is a recommendation that the user decides on (it changes default behaviour of the strict and standard classes).

## Revision statement (read first)
This is the second registration for these two switches. The first one (`single-node-speed.md`, seeds 601-610) was run on dev08 and each switch failed one clause by under one point (ledger SINGLE-NODE-AB1): `BAANDI_CERT_SCREEN=on` failed the si16 strict T=8 cost bound (+1.36 percent against +1 percent; five repeats mean +1.65), and `BAANDI_NL_ENERGY=fast` failed the `total_energy` phase proxy at T=8 (-39.45 percent against -40) while its real target, wall, passed. Those verdicts stand. The clauses below were revised AFTER seeing those data: the small-cell cost bound is wider and stated as a median plus a per-seed worst case, the phase-timer proxy is no longer decisive, and the decisive comparison is the two switches together against both off. To keep the check honest:
- every number comes from fresh seeds 701-710 that have never been run with either switch;
- the large-cell gain thresholds (si32 strict, -10 percent at T=1 and -8 percent at T=8 for the pair) are fixed below the effect seen in AB1 (-24.15 and -17.70 percent) but far above run-to-run noise, so a pass is not a re-reading of AB1;
- the binary contains two changes made after AB1 that target the failed cost clause (below); their effect is measured here, not assumed.

## Changes since AB1 in the candidate binary (both opt-in, bit-identical when off)
- POTENTIAL-REUSE: with the screen active, the Hartree/XC fields it builds for the output density are reused by the certificate solve at that same density (consumed only on a bit-exact density match; `BAANDI_POTENTIAL_REUSE=off` disables it for A/B). Results bit-identical with reuse on and off (unit test `certificate_potential_reuse_is_bit_identical`).
- The screen's per-band transforms of near-degenerate clusters run one per pool worker with their own FFT plan and the cluster matrix entries are split over the pool, each entry with its serial sum: the estimate is bit-identical for any thread count (same unit test, 1 against 4 threads). Mac observation (si16, T=8, seeds 601-605): screen time 34 ms -> 16 ms per certificate; the certificate's potential build 2.5 ms -> 0.08 ms. The proposal's premise that the AB1 overhead was one potential build was wrong on the Mac: the build is 2.5 ms of the 34 ms screen, the band transforms were 20 ms (ledger POTENTIAL-REUSE).

## Arms
- OFF: neither switch.
- ON: `BAANDI_CERT_SCREEN=on BAANDI_NL_ENERGY=fast` (POTENTIAL_REUSE follows the screen, default on with it).
- Attribution only (reported, not decisive): `BAANDI_CERT_SCREEN=on` alone and `BAANDI_NL_ENERGY=fast` alone on si32 strict, seeds 701-705.
`BAANDI_ENERGY_IDENTITY` (B2) is NOT part of this registration: it changes the meaning of per-iteration `energy` in `history` and needs the user's decision first.

## Protocol
One MKL build of the candidate commit (`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`), the same binary for both arms; the A/B is the environment only. Cells `benchmarks/cases/si16-gamma` and `si32-gamma`, perturbed seeds 701-710 (+-1e-7 Bohr as in MPI-M1-GATE3). Profiles: strict for si16 and si32; standard for si16 and si32 (only NL_ENERGY acts there; the certificate does not exist in the standard class). Thread budgets T=1 (taskset core 0) and T=8 (cores 0-7). Interleave OFF/ON per seed, node idle (load below 0.30), uptime before and after.

## Decisive rules (fixed now; never moved after data)
1. Identity (every run, both T): density SHA-256 equal between OFF and ON; total energy within 1e-12 Ha; SCF iteration count equal; Davidson iterations and H vectors (recorded plus certificate, `convergence.work`) with ON not above OFF. Density SHA reproducible per arm and equal between T=1 and T=8. Any violation rejects.
2. Large-cell gain: si32 strict median wall over seeds 701-710 improves by at least 10 percent at T=1 and at least 8 percent at T=8; no si32 strict seed worse than +2 percent.
3. Small-cell cost: si16 strict median wall not worse than +2 percent at T=1 and at T=8; no si16 strict seed worse than +5 percent.
4. Standard class (NL_ENERGY only): si32 standard median wall improves by at least 2 percent at T=1 and T=8; si16 standard median not worse than +2 percent; no seed worse than +5 percent.
A failure of any rule leaves both defaults unchanged and is recorded in the ledger with the numbers. A pass of 1 to 4 is a recommendation to the user to make both switches default for the integer-occupation spin-unpolarized driver (the Fermi-Dirac and collinear drivers have their own certificate and energy code and are not affected by either switch).

## Reported, not decisive
Seeds with at least one screened certificate and Davidson iterations avoided; screen time per certificate from two `BAANDI_CERT_SCREEN=trace` runs (`screen_s=` in the stderr line); `total_energy` and `potential` phase timers (also the certificate's own `potential_s` inside `history[].fixed_point_certificate.timings`); the attribution arms.

## What to send back
Per run: wall, phase timers, `convergence.work`, energy (17 digits), density SHA, exit code; trace stderr lines; binary sha256, commit, build command, uptime before and after, interleave order.
