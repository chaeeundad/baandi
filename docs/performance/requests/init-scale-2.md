# dev08 request: replicated init pieces (INIT-SCALE-2, 2026-11-18)

Ledger section INIT-SCALE-2 (rows INIT-EWALD-GATE, INIT-PACKED, INIT-REUSE, INIT-CLI, INIT-DENSITY-SPLIT). Every change is bit-identical (energy and density SHA must equal the base at every rank count) and has no switch, so the A/B is two binaries: base = main e902853, head = the INIT-SCALE-2 branch head. Load below 0.30, `--map-by core --bind-to core`, one socket. Report both binaries' sha256, raw stdout (full), load before and after; write the start/done/failed/idle signal files as usual.

## Build (both commits, separate trees)
```
cargo build --release -p dft-cli --features mpi,mkl
sha256sum target/release/dft-cli
```

## Runs: si32 and si64 standard f1e-3 seed 601, `scf.max_iterations = 10`, three alternating pairs per cell
```
for CELL in si32 si64; do for P in 2 4 8; do for ARM in base head base head base head; do
  env BAANDI_PHASE_PROFILE=on mpirun -np $P --map-by core --bind-to core $BIN_$ARM run <$CELL-10step.toml> \
    --output-dir ~/fastbench/init2-$CELL-np$P-$ARM-$RANDOM > ~/fastbench/init2-$CELL-np$P-$ARM-$RANDOM.stdout 2>&1
done; done; done
```

## Rows to report per cell and rank count (median of 3, max over ranks)
`init` wall and its parts `init.cli`, `init.basis`, `init.ionic_local`, `init.densities`, `init.states`, `init.ewald`, `init.rest`, `init.atomic_prepare`; the atomic-start parts (`atomic.span`, `atomic.orthonormalize`, `atomic.ritz`); total wall; energy bits and density SHA-256 (must equal base); peak RSS per rank if available (INIT-PACKED removes one projector table copy, about 98 MB per rank at si64).

## Expected (E, from Mac O)
- `init.states` -25 to -45 %, `init.rest` about -80 % at every rank count; `init.cli` down by one density sphere and two bases.
- `init.ewald` at np = 2: the gate makes it the serial sum (no collective). The Mac showed the serial sum SLOWER than the split at np 2 (0.089 -> 0.182 s, simulated ranks). Suspicion: the split's chained broadcasts are the first collective after the replicated `states` and absorb its rank skew in `init.ewald`. Decision rule: keep `EWALD_SPLIT_MIN_RANKS = 4` only if head's `init` total at np = 2 is not worse than base beyond the session noise floor; otherwise set it to 2 (one-line change, bits unchanged).
