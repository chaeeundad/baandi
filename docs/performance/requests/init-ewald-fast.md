# dev08 request: Ewald and ionic local construction (INIT-EWALD-FAST, INIT-IONIC-FAST, 2026-11-18)

Ledger section INIT-EWALD-FAST (rows INIT-EWALD-SLOWDOWN, INIT-EWALD-PRUNE, INIT-EWALD-FAST, INIT-IONIC-SHELL, INIT-IONIC-FAST). One binary, environment-only A/B:

- arm R (defaults): reference Ewald scheme with pruned image enumeration (bit-identical to eecd57f and to main), direct structure-factor phases with the last-shell form-factor lookup (bit-identical);
- arm F: `BAANDI_EWALD_FAST=on BAANDI_IONIC_FAST=on` (rounding class; Ewald truncation <= min(1e-13 (sum |Z|)^2 / Omega^(1/3), 1e-13 N) Ha);
- optional base B: main 9d942a0 (expected bits equal to R).

Load below 0.30, uptime before and after, `--map-by core --bind-to core` for MPI, `taskset -c 0` for one core; write the start/done/failed/idle signal files on the shared NFS dir as usual. Report the binary sha256 and the full raw stdout.

## Build
```
cargo build --release -p dft-cli --features mpi,mkl
sha256sum target/release/dft-cli
```

## Runs: si64 (and si32) standard, `scf.max_iterations = 10`, three alternating pairs
```
for CELL in si32 si64; do
  for ARM in R F R F R F; do
    case $ARM in R) ENV="";; F) ENV="BAANDI_EWALD_FAST=on BAANDI_IONIC_FAST=on";; esac
    # one core (input threads = 1)
    env $ENV BAANDI_PHASE_PROFILE=on taskset -c 0 $BIN run <$CELL-10step.toml> \
      --output-dir ~/fastbench/initfast-$CELL-t1-$ARM-$RANDOM > ~/fastbench/initfast-$CELL-t1-$ARM-$RANDOM.stdout 2>&1
    # np = 8
    env $ENV BAANDI_PHASE_PROFILE=on mpirun -np 8 --map-by core --bind-to core $BIN run <$CELL-10step.toml> \
      --output-dir ~/fastbench/initfast-$CELL-np8-$ARM-$RANDOM > ~/fastbench/initfast-$CELL-np8-$ARM-$RANDOM.stdout 2>&1
  done
done
```

## Rows to report (median of 3, max over ranks)
`init` wall and the construction rows `init.ewald` and `init.ionic_local` (plus `init.states`, `init.basis`, `init.rest`), total wall, total energy, Ewald energy, density SHA-256. Arm R energy and density SHA must equal base B (bit-identical changes). Arm F: |E_F - E_R| (expected about 2e-9 Ha at si64: the reference scheme's own summation rounding over 3.9M real-space terms; the balanced scheme agrees with a compensated tight reference to 8e-15 Ha/atom on the Mac).

## Expected (Mac O, release, one thread)
- `init.ewald` si64 one core: R about 0.95x of the eecd57f evaluator (dev08 eecd57f replicated about 0.245 s, main 0.813 s), F about 0.6 ms (116x below R).
- `init.ionic_local` si64 one core: R unchanged to -2 %, F about 4x lower (Mac 0.055 -> 0.013 s; dev08 0.390 s, where sin/cos is the larger share, so the gain may be larger).
- np = 8: the Ewald split stays on (bit-identical either side); with F the Ewald work is below the broadcast chain latency.
