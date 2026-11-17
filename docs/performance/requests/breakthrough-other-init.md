# dev08 request: BREAKTHROUGH-4-OTHER-AND-INIT (what was in solve.other, init and the remainder)

Ticket BREAKTHROUGH-4-OTHER-AND-INIT. The exclusive phase profile of a real-MPI si32 np=8 production run
(a220a85, relayed 2026-11-17) left 0.2601 s of the fine solve (24.7 %) as `solve.other`, 0.2552 s of the
wall outside every bucket and 0.5606 s in `init`. This ticket adds finer labels (all count-only friendly,
results bit-identical with the profile on or off) and three default-off switches that remove the
replicated work the labels exposed. Analysis and Mac numbers (label O, simulated ranks):
`docs/performance/profile-2026-11-17-other-init.md`; ledger section "BREAKTHROUGH-4".

Hard rule (fast-loop-protocol.md): every job is ONE case (si32, seed 601, f1e-3), np=8, finishes in well
under 10 minutes (a production si32 np=8 run is about 3.5 s; 3 pairs = 6 runs per job). MPI against MPI
only; no threads or one-core number appears here. QE is not involved.

## Build (one commit; record the sha256)
```
git fetch && git checkout <branch HEAD of this ticket>      # the orchestrator fills the hash after merge
cargo build --release -p dft-cli --features mpi,mkl --bin dft-cli
sha256sum target/release/dft-cli
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si32 --emit-input ~/fastbench/si32.toml
```
Load gate as always (below 0.30 before each run, `uptime` before and after), `MKL_NUM_THREADS=1`,
`mpirun -np 8 --map-by core --bind-to core`.

## Job 1 (stage 0/1 labels): profile with the new labels, default stack (about 30 s)
```
env MKL_NUM_THREADS=1 BAANDI_PHASE_PROFILE=on \
  mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si32.toml \
  --output-dir ~/fastbench/b4-profile-np8
env MKL_NUM_THREADS=1 BAANDI_PHASE_PROFILE=on BAANDI_PHASE_WAIT=on \
  mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si32.toml \
  --output-dir ~/fastbench/b4-wait-np8
```
Deliver both console tables (the block that starts with `phase profile (8 ranks; ...)`, including the new
lines `construction:`, `atomic start ... (span ... = radial ... )`, `between`, and `solve labels
(BREAKTHROUGH-4 ...)`) and both results.json `parallel.phases`. Expected: `solve.other` of the fine bucket
falls from 0.26 s to a few milliseconds and the named labels carry it (on the Mac: `real_state` 0.129 s,
`kept_update` 0.075 to 0.099 s, `setup` 0.015 to 0.045 s, `prepare` 0.017 to 0.028 s of the fine solve);
the remainder line falls below 1 % (the coarse construction is now `init`, the ladder hand-off, session
set-up, stop rule and verification glue are the `between` bucket). The second run is a DIAGNOSTIC (extra
barriers): its `solve.wait` / `entry_wait` / barrier seconds say how much of each stage is waiting for
the last rank; never use its wall. Energy bits and density sha256 must equal across the two runs.

## Job 2 (stage 1, SCF level): the three removal switches, stacked, against the default (about 1 min)
Three alternating pairs (A B, B A, A B), profile OFF (timing runs), to convergence:
```
A: env MKL_NUM_THREADS=1 mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si32.toml --output-dir ~/fastbench/b4-A-<i>
B: env MKL_NUM_THREADS=1 BAANDI_M2_LIGHT_STATE=on BAANDI_EWALD_CACHE=on BAANDI_M2_LAYOUT_CACHE=on \
   mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si32.toml --output-dir ~/fastbench/b4-B-<i>
```
Deliver per run: wall (results.json `timings.total_s` or the external wall), SCF steps, verification
attempts, H columns, energy bits and density.bin sha256. Expected: identical energy bits and density
sha256 in all six runs (all three switches are bit-identical), identical counts; Mac observation (O, np=8
simulated, 7 pairs): B/A median 0.936.

## Job 3 (stage 1, one switch at a time; only if Job 2 shows a gain; about 3 min)
Same as Job 2 with B = each switch alone (`BAANDI_M2_LIGHT_STATE=on`, `BAANDI_EWALD_CACHE=on`,
`BAANDI_M2_LAYOUT_CACHE=on`), three pairs each, plus one profile run per switch (`BAANDI_PHASE_PROFILE=on`)
so the moved label is visible (`solve.real_state` of the coarse and fine buckets, `init.ewald`,
`solve.setup`). Mac (O, 5 pairs each): light state B/A 0.918, Ewald cache 0.985, layout cache (see the
note).

## Expansion list (stage 2; NOT part of this request, queued only after Job 2/3 show a candidate)
si64 np=8, si32 np=4 and np=2, mgo-40 or a skewed cell np=8, the threads path (Ewald cache only), MKL and
pure builds.

## What comes back into the repository
The orchestrator files the tables under `docs/performance/profile-2026-11-17-other-init.md` (section
"dev08, label D") and the Job 2/3 rows in the ledger section BREAKTHROUGH-4 (relayed, with binary sha256,
load lines and the stored result files).
