# dev08 request: owned potential and density of the M2 path (M2-OWNED, 2026-11-18)

Ledger section M2-OWNED. MPI against MPI only (comparison-rules.md); no QE column. Load below 0.30, `--map-by core --bind-to core`, one socket (CPUs 0-7), the master's start/done/failed/idle signal files on the shared directory. Report binary sha256, full raw stdout, load before and after.

Trigger: with `BAANDI_M2_GRID=on BAANDI_M2_MIX=on` (dev08 6be8f09, np=8, si64 10 fixed maps) every map still makes three full-grid allgathervs of N doubles per rank: the effective potential (for the operator), the output density (M2 density build) and the mixed density (distributed Pulay). The distributed H apply reads the potential only on the rank's own FFT planes.

## Switches (default off, independent; both need `BAANDI_M2_GRID=on`, DENSITY also `BAANDI_M2_MIX=on`)
- `BAANDI_M2_LOCAL_FIELDS=on`: the effective potential stays on the slab planes and goes straight into the distributed operator (epoch-tagged; one alltoall only when the slab and operator plane axes differ, which never happens for the si cells), the kept-start delta V runs on owned planes (one extra 1-double tree sum per kept start). The complex Hamiltonian's full potential is gathered only for a consumer that needs it (cold atomic start, complex or M1 fallback). Bit-identical to GRID alone. Expected (E): one N-double allgatherv fewer per warm map, `potential` per step down a little.
- `BAANDI_M2_LOCAL_DENSITY=on` (with FIELDS): the output density stays on the planes through the audit, the stop rule, the session metrics (one small allgatherv of per-plane partials), the plane-mode Pulay mixer and the next potential; the full density is gathered once at delivery. Rounding class (per-plane canonical sums). Expected (E): two N-double allgathervs fewer per map, `density` and `mixing` per step down.

Counters (results.json `parallel.m2.counters`, `owned.*`; seconds in `parallel.m2.stages`, `owned.*`; also with the switches off, so arm A gives the baseline): `potential_gathers`, `density_gathers`, `mixer_gathers` (+ `_values`), `consumer_potential_gathers`, `density_materializations`, `redistributions`, `full_copies` (+ `_values`), `fields_owned_maps`, `density_owned_maps`, `fields_declined.<reason>`, `epoch_checks`, `owned_reductions`; seconds `potential_gather`, `density_gather`, `mixer_gather`, `consumer_potential_gather`, `density_materialize`, `redistribute`, `owned_reduce` (0 with `BAANDI_STAGE_TIMERS=off`). Values are doubles per rank (8 bytes each).

## Stage 1: 10-step fixed screens, np=8 (the M2-GRID protocol)
Same binary (`cargo build --release -p dft-cli --features mpi,mkl` at the merge of this branch), si32 and si64 standard f1e-3 seed 601, `scf.max_iterations = 10`, `BAANDI_PHASE_PROFILE=on`, three alternating pairs per cell and arm:
- arm A (`off`): `BAANDI_M2_GRID=on BAANDI_M2_MIX=on`;
- arm F: arm A + `BAANDI_M2_LOCAL_FIELDS=on`;
- arm FD: arm F + `BAANDI_M2_LOCAL_DENSITY=on`.
```
for CELL in si32 si64; do for ARM in A F FD A F FD A F FD; do
  E="BAANDI_M2_GRID=on BAANDI_M2_MIX=on"
  case $ARM in F) E="$E BAANDI_M2_LOCAL_FIELDS=on";; FD) E="$E BAANDI_M2_LOCAL_FIELDS=on BAANDI_M2_LOCAL_DENSITY=on";; esac
  T=m2owned-$CELL-np8-$ARM-$RANDOM
  env BAANDI_PHASE_PROFILE=on $E \
    mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run <$CELL-10step.toml> \
    --output-dir ~/fastbench/$T > ~/fastbench/$T.stdout 2>&1
done; done
```
Report per cell and arm (fine bucket, per step and total, max over ranks): `potential`, `density`, `energy`, `mixing`, `solve`, `between.stop_rule`, `coll.allgatherv` and `coll.alltoall` (seconds, calls), the `owned.*` counters and seconds, wall, SCF count, energy bits, density SHA-256, the switch manifest (`scf.switches`). Arm F must have the energy bits and density SHA of arm A; arm FD the standard class (|dE| <= 5e-9 Ha) and identical or documented SCF counts.

Decision (adoption rule of 2026-11-10): adopt a switch as default when its median wall improves on both cells and nothing regresses beyond run-to-run noise; correctness gate as above, verification passed. Note: with DENSITY an exact checkpoint at an active boundary rejects (`Boundary`: the input density is owned; `ScfSession::materialize_density` gathers it), and runtime mixer replacements are refused.
