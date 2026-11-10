# dev08 request: single-node speed, two opt-in switches (pre-registered 2026-11-10)

Branch: worktree-agent-a4d64fc7edb774116 (threads-only; no MPI file touched). Both switches are environment variables read once per process, default off, and the default path is bit-identical to the parent commit (checked on Mac, see the ledger rows CERT-SCREEN and NL-ENERGY). Build exactly like the other A/B binaries (MKL, `MKLROOT=/opt/intel/oneapi/mkl/2026.1`, one binary, the A/B is the environment variable only).

## What the switches do
- `BAANDI_CERT_SCREEN=on` (strict class only). When an SCF iteration meets the density, energy and eigen-residual criteria, the fixed-point certificate re-solves the Kohn-Sham problem at the output density with a 10x tighter tolerance and requires the occupied eigenvalues to move by at most 1e-7 Ha. Before that solve, a first-order estimate of that eigenvalue movement is computed from the stored per-band densities (single bands) and a small cluster matrix (bands within 1e-3 Ha of each other): cost one potential build and a few FFTs. If the estimate exceeds 1.02 x 1e-7 Ha the certificate is declared failed without the solve. The re-solve result of a failed certificate is discarded by the SCF loop in any case, so the SCF trajectory is unchanged; only the wasted work disappears. A screened record in `history[].fixed_point_certificate` has `davidson_iterations = 0`, `operator_applications = 0`, the step's own residuals, and `occupied_eigenvalue_change_ha` = the estimate. `BAANDI_CERT_SCREEN=trace` runs the solve as usual and prints one `CERT-SCREEN iter=.. bound=.. actual=..` line per certificate to stderr (diagnostic).
- `BAANDI_NL_ENERGY=fast`. The nonlocal energy of the energy audit (the `total_energy` phase) is computed from the projector products only (no back-application), for occupied bands only, four bands per pass over a projector column. Rounding differs from the legacy path (observed total-energy difference 1.3e-13 Ha on si16); the per-band values do not depend on grouping or thread count.

## Pre-registered acceptance (fixed before any dev08 data; never moved after seeing data)
Cells: si16 and si32 (`benchmarks/cases/si16-gamma`, `si32-gamma`), perturbed seeds 601-610 (+-1e-7 Bohr as in MPI-M1-GATE3), strict class for both switches and standard class for `BAANDI_NL_ENERGY` only (the certificate does not exist in the standard class). Thread budgets T=1 (taskset core 0) and T=8 (cores 0-7). MKL build, `MKL_NUM_THREADS=1`, interleaved off/on per seed, node idle (load below 0.30), the same binary for both. One switch at a time, then both together as a final check.

`BAANDI_CERT_SCREEN=on` (trajectory-neutral by construction):
1. Identity: density SHA-256 and total energy bit-identical between off and on in every run (every seed, T=1 and T=8). Any difference rejects the switch.
2. No run does more Davidson iterations or H vectors with the switch on (the sum of recorded and certificate work in `convergence.work`).
3. Wall: median over the ten si32 strict seeds improves by at least 3 percent at T=1 and at least 2 percent at T=8; the worst seed is not worse than +2 percent; the si16 strict median is not worse than +1 percent.
4. Reported (not decisive): number of seeds with at least one screened certificate, Davidson iterations avoided, the screen's own cost (estimate-only runs with `trace` on two seeds).

`BAANDI_NL_ENERGY=fast` (rounding change):
1. Energy within 1e-12 Ha of the off run in every run; Davidson totals and SCF iteration counts unchanged (a changed count is reported and judged by the energy tolerance and convergence, as in the P14B rule); density SHA reproducible per setting and equal between T=1 and T=8 within the setting.
2. The `total_energy` phase timer drops by at least 40 percent (median over seeds) at T=1 and T=8 on si32 strict.
3. Wall: median improves by at least 2 percent at T=1 and T=8 on si32 strict and on si32 standard; no seed worse than +1 percent.

A switch that fails a rule is recorded as rejected in the attempts ledger with the numbers; the default does not change. Passing both is a recommendation to the orchestrator (a default change is a separate decision with its own A/B of the default build).

## What to send back
Per run: wall, phase timers, `convergence.work`, energy (17 digits), density SHA, exit code; for `trace` runs the stderr lines. Binary sha256, commit, build commands, uptime before and after, interleave order.
