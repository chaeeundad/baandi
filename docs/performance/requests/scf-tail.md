# dev08 request: SCF tail and the verification estimator (SCF-TAIL), 2026-11-18

Ledger rows SCF-TAIL-ACCOUNTING, SCF-TAIL-PREDICTOR, SCF-TAIL-STALL. MPI against MPI and threads against threads only (comparison-rules.md); no QE column. Load below 0.30, `--map-by core --bind-to core`, one socket (CPUs 0-7), the master's start/done/failed/idle signal files on the shared directory. Report binary sha256, raw stdout AND stderr (the trace lines go to stderr), load before and after.

## Why
dev08 si64 np=8 seed 601: the verification at map 10 failed with coulomb_error 1.0e-8 Ha while the true error of the verified state was 7.7e-10 Ha; map 11 and a second attempt cost 11.3 % of the wall. Mac (counts, 100 runs, 10 systems x 5 seeds x thresholds 5e-9 / 1e-10): every one of 64 failed attempts verified a state already inside its class. In the linear Hartree model the verification criterion reads A (1 + A) times the true error of the delivered state (A = Coulomb-weighted static response; eps (eps - 1) for screened long-wavelength modes, so it grows with the cell), the stop estimate (1 + 1/A) and the new cross estimate -X = -<rho_out - rho_in, F(rho_out) - rho_out>_C reads (1 + A). Measured -X / T ranges 0.654 to 18 on the Mac set (not a bound), so the opt-in `on` passes only below threshold / 4; on the Mac set it changes 10 of 100 runs and saves 0.36 % of the H columns. The open question for dev08 is -X / T on si32 and si64 (larger A): if -X stays well above T there, `on` does nothing on the cells that matter, and the switch is closed.

## Arms (to convergence, production defaults, f1e-3, seeds 601-603)
The cross estimate needs the full stop-map densities; the owned-plane loop (default at np >= 2) keeps them on the planes and the trace prints nothing there. MPI arms therefore pin `BAANDI_M2_LOCAL_DENSITY=off` on all three arms.
- A: `BAANDI_M2_LOCAL_DENSITY=off` (MPI) / defaults (threads).
- T: A plus `BAANDI_STOP_PREDICTOR=trace` (must be bit-identical to A: energy bits, SCF count, density SHA).
- O: A plus `BAANDI_STOP_PREDICTOR=on`.
```
for SEED in 601 602 603; do for ARM in A T O; do
  case $ARM in A) E="";; T) E="BAANDI_STOP_PREDICTOR=trace";; O) E="BAANDI_STOP_PREDICTOR=on";; esac
  IN=<si64-production-seed$SEED.toml>; R=tail-si64-np8-$SEED-$ARM
  env BAANDI_M2_LOCAL_DENSITY=off $E mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run $IN \
      --output-dir ~/fastbench/$R > ~/fastbench/$R.stdout 2> ~/fastbench/$R.stderr
done; done      # then si32 np=8; threads: si32 T=8 (arms A T O without the M2 variable)
# strict reference per geometry (true errors): BAANDI_REFERENCE_STRICT=on, same input, same np
```
Per run: wall, maps, H columns (fine recorded + certificate), verification attempts (passed/failed, coulomb_error_ha, h_vectors). Per map (history): energy error against the strict reference, standard_error_ha, density_residual, requested_eigen_tolerance_ha, operator_applications, davidson_iterations, pulay.kept and pulay.gram_condition. Per attempt from stderr (`STOP-PREDICTOR` lines, printed once per rank with identical values; take one): stop_ha, verify_ha, cross_ha (= -X), verify_energy_change_ha.

## Decision
`on` is a default candidate only if (1) on si32 and si64 the ratio -X / T of every attempt with T >= 1e-11 Ha stays above 1 / 4 x 2 = 0.5 (the Mac floor 0.654 with a factor two margin), (2) the median wall improves and nothing regresses, (3) the delivered true error stays below 5e-9 Ha against the strict reference. Otherwise close the switch (keep `trace` as an instrument). Changing the verification criterion itself (ADR 0043 V2) remains a user decision.
