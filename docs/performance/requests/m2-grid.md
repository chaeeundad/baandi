# dev08 request: distributed full-grid stages of the M2 path (M2-GRID, 2026-11-18)

Ledger section M2-GRID (rows M2-GRID-MAP, M2-GRID-POTENTIAL-AUDIT, M2-GRID-MIX, M2-GRID-PLAN). MPI against MPI only (comparison-rules.md); no QE column. Load below 0.30, `--map-by core --bind-to core`, one socket (CPUs 0-7), the master's start/done/failed/idle signal files on the shared directory. Report binary sha256, raw stdout (full), load before and after.

Trigger (master, 2026-11-18, e902853, per SCF step, np=8): outside the solve `potential` 0.0228 / 0.0522 s, `mixing` 0.0124 / 0.0262 s, `energy` 0.0208 / 0.0482 s (si32 / si64): replicated full-grid work that does not shrink with the rank count.

## Switches (both default off, rounding class, independent)
- `BAANDI_M2_GRID=on`: Hartree FFT pair, LDA XC and effective potential of the input density on the planes of a slab-decomposed full-grid FFT (one allgatherv of the effective potential and partial sums); the energy audit's output-density transform and grid sums on the planes (one allgatherv of partial sums); the standard stop-rule estimate from the kept owned coefficients (one small allgatherv). Expected (E): `potential` and `energy` per step down by 40 to 70 percent at np=8, `between.stop_rule` unchanged or lower.
- `BAANDI_M2_MIX=on`: plain Pulay mixing with the window on the rank's grid block (4096-point chunks), Gram row and charge sums as chunk partials, one allgatherv of the new density. Expected (E): `mixing` per step down by 60 to 80 percent at np=8.

Both keep every rank bit-equal and give the same bits for every rank count (canonical plane/stick/chunk sums); against the default they move the energy by rounding only (Mac si32 np=8: |dE| 2.9e-12 Ha, SCF/Davidson/H counts identical).

## Stage 1: 10-step fixed screens at np 4 and 8 (the eecd57f finding protocol)
Same binary (`cargo build --release -p dft-cli --features mpi,mkl` at the merge of this branch), si32 and si64 standard f1e-3 seed 601, `scf.max_iterations = 10`, `BAANDI_PHASE_PROFILE=on`, three alternating pairs per cell and arm, first pair discarded if the session noise floor demands it:
- arm A: defaults;
- arm G: `BAANDI_M2_GRID=on`;
- arm GM: `BAANDI_M2_GRID=on BAANDI_M2_MIX=on`.
```
for P in 4 8; do for ARM in A G GM A G GM A G GM; do
  case $ARM in A) E="";; G) E="BAANDI_M2_GRID=on";; GM) E="BAANDI_M2_GRID=on BAANDI_M2_MIX=on";; esac
  T=m2grid-si32-np$P-$ARM-$RANDOM
  env BAANDI_PHASE_PROFILE=on $E \
    mpirun -np $P --map-by core --bind-to core target/release/dft-cli run <si32-10step.toml> \
    --output-dir ~/fastbench/$T > ~/fastbench/$T.stdout 2>&1
done; done          # then the same for si64
```
Rows to report per cell and arm (fine bucket, per step and total, max over ranks): `potential`, `energy`, `mixing`, `density`, `solve`, `between.stop_rule`, `coll.allgatherv` and `coll.alltoall` (seconds and calls), wall, SCF count, energy bits, density SHA-256. The switch manifest in results.json (`scf.switches`) must show the arm's values.

Decision (adoption rule of 2026-11-10): adopt each switch as the default when its median wall improves at np 4 and 8 on both cells and nothing regresses beyond run-to-run noise; correctness gate: standard class (|dE| <= 5e-9 Ha against arm A), identical SCF counts or a documented rounding-edge difference, verification passed.

## Stage 2 (after Stage 1): production run to convergence
si32 and si64 production defaults (ladder, verification) at np=8, arms A and GM, three alternating pairs: wall, SCF count, energy, verification result. Note: with `BAANDI_M2_MIX=on` an exact checkpoint capture rejects the mixer (opaque state); report if the production CLI path writes one.
