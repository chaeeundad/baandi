# dev08 request: owned defaults and the distributed Coulomb guard (M2-OWNED-DEFAULT, MIX-COULOMB-GUARD-DIST, 2026-11-18)

Ledger rows M2-OWNED-DECLINE, M2-OWNED-DEFAULT, MIX-COULOMB-GUARD-DIST. MPI against MPI and threads against threads only (comparison-rules.md); no QE column. Load below 0.30, `--map-by core --bind-to core`, one socket (CPUs 0-7), the master's start/done/failed/idle signal files on the shared directory. Report binary sha256, full raw stdout, load before and after.

Binary: `cargo build --release -p dft-cli --features mpi,mkl` at the merge of this branch.

## Stage 1: owned defaults (10 fixed maps, the M2-OWNED protocol)
`BAANDI_M2_LOCAL_FIELDS` and `BAANDI_M2_LOCAL_DENSITY` are now default on. Arm D (defaults, no switch env) must reproduce the earlier FD arm (energy bits, density SHA, counters); arm A pins them off.
```
for CELL in si32 si64; do for ARM in D A D A D A; do
  E=""; [ $ARM = A ] && E="BAANDI_M2_LOCAL_FIELDS=off BAANDI_M2_LOCAL_DENSITY=off"
  T=owneddef-$CELL-np8-$ARM-$RANDOM
  env BAANDI_PHASE_PROFILE=on $E mpirun -np 8 --map-by core --bind-to core \
    target/release/dft-cli run <$CELL-10step.toml> --output-dir ~/fastbench/$T > ~/fastbench/$T.stdout 2>&1
done; done
```
Report wall, energy bits, density SHA-256, `scf.switches.BAANDI_M2_LOCAL_{FIELDS,DENSITY}.{path,fallback_reason}` (expected `owned`, null), `parallel.m2` `owned.*` counters including `owned.density_declined.*`.

## Stage 2: safety checks (one run each, np=4, si32)
- `dft-cli checkpoint <si32.toml> cp.bin --iteration 3 --kind exact` then `dft-cli restart <si32.toml> cp.bin out.bin`: both succeed (before: `Rejection::Mixer`/`Boundary`); the restart result equals an uninterrupted `checkpoint --kind result` run's energy within 5e-9 Ha.
- `dft-cli run <si32.toml> --observe` (owned path, `path = owned`) and `--agent` (declined, `fallback_reason` contains `runtime_control`).
- a PBE input and a spin input: `fallback_reason` `functional` / `spin`, no error.

## Stage 3: distributed Coulomb guard (to convergence, production defaults)
Arms: A defaults; G `BAANDI_MIX_COULOMB_GUARD=on` (now acts on the distributed mixer, owned density on; no `BAANDI_M2_MIX=off` needed); T `BAANDI_MIX_COULOMB_GUARD=trace` (must equal A bit for bit). si32 and si64 at np 4 and 8, seeds 601-603, three alternating pairs.
```
for CELL in si32 si64; do for NP in 4 8; do for ARM in A G A G A G T; do
  E=""; case $ARM in G) E="BAANDI_MIX_COULOMB_GUARD=on";; T) E="BAANDI_MIX_COULOMB_GUARD=trace";; esac
  T=guarddist-$CELL-np$NP-$ARM-$RANDOM
  env $E mpirun -np $NP --map-by core --bind-to core target/release/dft-cli run <$CELL-production.toml> \
    --output-dir ~/fastbench/$T > ~/fastbench/$T.stdout 2>&1
done; done; done
```
Report per run: wall, SCF maps (coarse + fine), H columns, Davidson iterations, verification attempts, |dE| against A, density SHA, and per map `history[].pulay` (`differences`, `coulomb_predicted_ratio`, `coulomb_refit`).

Expected (Mac counts, E): H to convergence si16 -4 to -16 %, si32 -7 to -11 % (production + ladder, simulated P 2..4), counts equal to the replicated guard; |dE| <= 1e-10 Ha.

Decision (adoption rule of 2026-11-10): G becomes the default (`BAANDI_MIX_COULOMB_GUARD=on`, replicated and distributed together) when the median wall to convergence improves on si32 and si64 at np 4 and 8 (and on the threads arms of scf-trajectory.md) and nothing regresses beyond noise; correctness: standard class against A, verification passed.
