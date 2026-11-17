# dev08 request: BREAKTHROUGH-0-INSTRUMENT (real MPI runs of the new measurement switches)

Ticket BREAKTHROUGH-0-INSTRUMENT (review 2026-11-17 sections 13.1 to 13.3). Everything here is measurement
infrastructure: three default-off switches (`BAANDI_PHASE_PROFILE`, `BAANDI_FIRST_BLOCK`, `BAANDI_SNAPSHOT_DIR`)
and a replay mode of `davidson_bench`. No speed claim is made and nothing is adopted; the request asks for the
real-MPI numbers that the Mac (simulated ranks, label O) cannot give: the exclusive phase profile of a production
si32 np=8 run, the first-block residual table of the same run, and one snapshot replay that proves the capture
reproduces the captured solve under real MPI. Hard rule (fast-loop-protocol.md): every job below finishes in
well under 10 minutes (a si32 np=8 production SCF is about 3 s on dev08; the capture run writes about 30 files
of 16 MB per rank).

Comparison rule: MPI against MPI only; no threads or one-core number appears in this request. QE is not involved.

## Build (one commit, record the sha256)
```
git fetch && git checkout <branch HEAD of this ticket>      # the orchestrator fills the hash after merge
cargo build --release -p dft-cli --features mpi,mkl --bin dft-cli --example davidson_bench
cargo build --release -p dft-engine --features mpi,mkl --test m2_entry_cost   # not needed for MPI; the CLI is the MPI harness
sha256sum target/release/dft-cli target/release/examples/davidson_bench
```
Input: `python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si32 --emit-input ~/fastbench/si32.toml`
(the stage-1 displacement stream, seed 601, amplitude 1e-3; same file for every job below). Load gate as always:
load below 0.30 before each run, `uptime` before and after, `MKL_NUM_THREADS=1`, `--map-by core --bind-to core`.

## Job A: exclusive phase profile, si32 np=8, one production run (about 10 s)
```
env MKL_NUM_THREADS=1 BAANDI_PHASE_PROFILE=on \
  mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si32.toml \
  --output-dir ~/fastbench/b0-profile-np8
```
Deliver: `~/fastbench/b0-profile-np8/results.json` (fields `phase_profile` and `parallel.phases`: per bucket the
maximum and the minimum over ranks of every entry) and the console table (the block that starts with
`phase profile (8 ranks; ...)`). The stage timers are on by default (`BAANDI_STAGE_TIMERS=on`); run the same
command once more with `BAANDI_PHASE_PROFILE=off` and once with `BAANDI_STAGE_TIMERS=off BAANDI_PHASE_PROFILE=on`
(count-only mode: the counts and the map clocks fill, the solve stages read 0) and report the three walls, so the
timer cost of the profile is visible (the review's "collective timer" caveat). Expected: the three runs have the
same energy bits and density sha256 (the switches never touch the arithmetic).

## Job B: first-block residual table, si32 np=8, one production run (about 10 s)
```
env MKL_NUM_THREADS=1 BAANDI_FIRST_BLOCK=on \
  mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si32.toml \
  --output-dir ~/fastbench/b0-firstblock-np8
```
Deliver the console table (`first-block residual diagnostic (...)`, one line per warm solve) and
`results.json` field `convergence.work.first_block`. The run's energy and density must equal Job A's (the
diagnostic is read-only; it applies H to the start block only when no kept H X exists, which on the default
stack with `M2_REUSE=on` happens for the first warm solve of each ladder stage only). Not a timing run.

## Job C: snapshot capture and replay, si32 np=8 (capture about 15 s, replay seconds)
```
mkdir -p ~/fastbench/b0-snap-np8
env MKL_NUM_THREADS=1 BAANDI_SNAPSHOT_DIR=~/fastbench/b0-snap-np8 \
  mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si32.toml \
  --output-dir ~/fastbench/b0-snap-run
ls ~/fastbench/b0-snap-np8/*.json          # one manifest per solve: s<seq>-<stage>-<kind>-it<n>.json
```
Pick three manifests: the first fine kept-start map (`s0xx-fine-ordinary-it02.json`, `"start": "kept"`), a late
fine map (the last `fine-ordinary`), and the verification map (`s0xx-verification-verify-itNN.json`). Replay each
with the same rank count:
```
env MKL_NUM_THREADS=1 mpirun -np 8 --map-by core --bind-to core \
  target/release/examples/davidson_bench --input ~/fastbench/si32.toml \
  --snapshot ~/fastbench/b0-snap-np8/<manifest>.json --repeats 5 | grep -E "replay of|DAVBENCH_JSON"
```
Deliver the `replay of ...` line of each (it prints `work reproduced true/false` and `eigenvalues reproduced
true/false` against the manifest's `result`) and the JSON lines. Expected: work reproduced `true` for every
manifest; eigenvalues reproduced `true` as well on the same binary and switches (a coarse-stage manifest needs
no extra option: the bench rebuilds the coarse problem from the manifest's cutoff). This is the precondition for
the R1/R2/R3 experiments to use replays of real SCF states instead of the synthetic warm map.
A/B on a snapshot: run the replay twice with an environment switch on one arm (for example
`BAANDI_M2_LEAN_ALGEBRA=on`): the solve time and the counters of the two arms are the stage-0 comparison at a
REAL mid-SCF state; three alternating pairs as in the protocol.

## What comes back into the repository
The orchestrator files the tables under `docs/performance/profile-2026-11-17-phases.md` (section "dev08, label D")
and `docs/performance/first-block-residual-2026-11-17.md`, with the binary sha256, the load lines and the
relayed JSON kept on dev08 (`~/fastbench/b0-*`). No ledger row changes its verdict from these runs; they fix the
measurement baseline that the R1 to R3 tickets will be judged against.
