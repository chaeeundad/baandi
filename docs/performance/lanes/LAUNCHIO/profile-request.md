# LANE-LAUNCHIO Phase 0 profile request

Status: **PROFILE_REQUEST**. No optimization, candidate selection, default change,
or Mac performance conclusion. Phase 2 starts only after the master returns the
raw dev08 profiles. Base: `61822bca3b470df1212ee47b47c7eb4e786e22e3`.

The ticket's supplied dev08 Si64 MPI8 profile is context, not remeasured evidence:
materialize/io 0.296 s; between 0.364 s (session setup 0.235 s, ladder handoff
0.123 s); outside-process launch/exit/output approximately 0.84 s. These nested
and approximate numbers must not be added into a new closure. No QE comparison
is made here. Prior F-LAUNCH-IO already instruments MPI lifecycle and output;
G-HANDOFF-OWNED is closed after a whole-run regression. Neither optimization is
reopened in this phase.

## Implementation and counter contract

Only `crates/dft-cli/src/run.rs` changes. Existing `BAANDI_PHASE_PROFILE=on`
(environment or `[profile] phase_profile = true`) enables all additions.
Off adds no clock reads, allocated records, file reads/writes or collectives;
there are only conditional checks at CLI boundaries. Existing four main-entry
bootstrap stamps are unchanged. `BAANDI_STAGE_TIMERS=off` does **not** disable
launch or between-map wall clocks; it disables the existing fine-stage clocks.
No numerical expression, iteration, reduction, delivery or default is changed.

Every existing rank directory retains its own results; parser aggregation is
post-process. There is no diagnostic MPI collective or barrier. Successful
closed-shell `run` is the supported measurement path, including Gamma and k
points. Spin/FD early returns and checkpoint commands are outside this report.
Use classic, verified, full-convergence runs, since output and finalization are
part of the question. Do not use a ten-step nonconverged screen for this lane.

Exact paths in each `results.json`:

| Path | Meaning |
| --- | --- |
| `launchio.schema`, `.rank` | `launchio-1`, local rank |
| `launchio.sizes` | `input_bytes`, `density_bytes`, `grid_points`, `plane_waves`, `bands`, `history_entries`, `atoms` |
| `launchio.upf_records[]` | One record per species in the run's second UPF pass: `species`, `bytes`, `projectors`, `read_s`, `parse_s`, `hash_s` |
| `launchio.buckets[]` | All eight existing bucket labels; rank-local `seconds.wall`, `seconds.io`, `seconds["between.session_setup"]`, `seconds["between.ladder_handoff"]`, `seconds["between.stop_rule"]`, `seconds["between.verify_glue"]`, `seconds["between.ladder_truncation"]`; `counts.maps` |
| `launch.seconds` | Existing lifecycle timers, plus `input_read`, `input_parse`, `scf_driver`, `density_encode`, `density_sha256`, `history_json` |
| `launch.milestones_s.scf_entry`, `.scf_return` | SCF driver boundaries on this process's monotonic clock |
| `launch.seconds.mpi_init`, `.mpi_finalize` | Existing lifecycle spans; finalize includes waiting for other ranks |
| `launch.seconds.preflight_including_prepare_and_upf`, `.upf_read_parse_hash` | First pass inclusive vs second pass inclusive; new per-species records refine only the latter |
| `launch.seconds.materialize_including_outputs`, `.results_serialize`, `.results_write`, `.density_bin_write`, `.density_json_write`, `.volumetric_write` | Existing output boundaries, plus new density/history subspans |
| `launch.seconds.profile_gather`, `.profile_report`, `.run_block`, `.binary_sha256`, `.console_write`, `.run_teardown`, `.communicator_drop` | Existing reporting, hashing, output and drop spans |
| `launch.main_wall_s`, `.covered_s`, `.outside_covered_s`, `.closure` | Existing rank-local closed wall tree, ending before the late launch append |

UPF bytes are bytes delivered by `read`, not physical disk reads. The preflight
pass still includes its own UPF read/parse, basis preparation and checks; it is
outside owned files and is not split here. `density_bytes = 8 * grid_points`
is the encoded payload; encoding reads an equal density payload and writes that
many output bytes, while SHA reads the payload again. This is a lower-bound
traffic input, **not** allocator or DRAM traffic. Parser also reports actual
`results_bytes_including_launch_append`. Hash/compiler allocations, filesystem
cache residency, system calls and serialization allocation counts remain unknown.
The GNU time per-rank sidecars supply peak RSS, page faults, CPU time and I/O
operation counts; they do not supply allocation counts.

`scf_driver` includes the cutoff ladder, session setup, handoff and verification;
never add these child clocks to it. Bucket wall values are exclusive **on one
rank**; `between.*` and `io` are subsets, not extra buckets. Parser exports
bucket means per map (null when maps=0), not a per-map latency distribution.
Between-glue call counts are not available; no engine code is changed to invent
them. Do not divide the between bucket's zero map count or attribute all handoff
time to an arbitrary fine map. The local clocks diagnose rank imbalance without
summing maxima from different ranks.

## dev08 commands (master relay only)

Reserve an exclusive socket-0 window. Use the same production MKL/MPI toolchain,
flags and algorithm switches as the supplied A profile. Record the reviewed
commit and binary digest. The script requires Linux host `dev08`, Python 3.11+
(`tomllib`), GNU time, and OpenMPI. It deliberately refuses to profile on Mac.
No package installation or remote execution is performed by Phase 0.

Build once outside the measurement window, from the reviewed checkout:

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cargo build --offline --release -p dft-cli --features mpi,mkl
sha256sum target/release/dft-cli
rustc --version
git rev-parse HEAD
```

Do not change previously registered release flags or switches for this build.
The master stages the exact production Si16/32/64 A-profile inputs at the paths
below. They must use absolute external UPF paths and `scf.threads = 1`, unchanged
cutoff/grid/start/bands/verification/geometry/seed. The runner records hashes of
the original input and every referenced UPF; its preserved input copy is a
provenance copy, not a relocated runnable input. No UPF is committed.

One command per job; do not overlap jobs or builds:

```sh
python3 -I docs/performance/lanes/LAUNCHIO/profile_dev08.py --case si64 --ranks 8 --input /data/sig-dft/lanes/LAUNCHIO/inputs/si64.toml --binary target/release/dft-cli --output /data/sig-dft/lanes/LAUNCHIO/si64-np8
python3 -I docs/performance/lanes/LAUNCHIO/profile_dev08.py --case si32 --ranks 8 --input /data/sig-dft/lanes/LAUNCHIO/inputs/si32.toml --binary target/release/dft-cli --output /data/sig-dft/lanes/LAUNCHIO/si32-np8
python3 -I docs/performance/lanes/LAUNCHIO/profile_dev08.py --case si16 --ranks 8 --input /data/sig-dft/lanes/LAUNCHIO/inputs/si16.toml --binary target/release/dft-cli --output /data/sig-dft/lanes/LAUNCHIO/si16-np8
python3 -I docs/performance/lanes/LAUNCHIO/profile_dev08.py --case si32 --ranks 1 --input /data/sig-dft/lanes/LAUNCHIO/inputs/si32.toml --binary target/release/dft-cli --output /data/sig-dft/lanes/LAUNCHIO/si32-onecore
python3 -I docs/performance/lanes/LAUNCHIO/profile_dev08.py --case si64 --ranks 1 --input /data/sig-dft/lanes/LAUNCHIO/inputs/si64.toml --binary target/release/dft-cli --output /data/sig-dft/lanes/LAUNCHIO/si64-onecore
```

The input and output paths are explicit staging contracts, not claims these
files already exist. Preserve one thread per rank. One core uses `taskset -c 0`
with the same MPI-capable binary directly; MPI8 uses
`mpirun -np 8 --cpu-set 0-7 --map-by core --bind-to core --report-bindings`.
Confirm these logical IDs are socket-0 physical cores on dev08 using `lscpu -e`;
preserve its output and bindings. Reject a job with incorrect binding instead
of interpreting it as a performance result. Do not compare across parallel models.

Each job discards warm-up pair 0, then keeps **5 crossed off/on pairs**. Off uses
PHASE_PROFILE/STAGE_TIMERS off; on uses both on. This estimates instrumentation
perturbation, not an optimization gain. Whole launches include verification,
output, GNU-time wrappers and MPI teardown. The load gate is <0.30 before each
launch, waiting in ten-second intervals for at most ten minutes. Output must
be fresh; failed runs stop and remain preserved. No reusing truncated data.

The runner writes `manifest.json` (all exact argv, binary/input/UPF hashes and
relevant inherited environment), per-launch stdout/stderr, `launch-wall.json`
(`wall_s`, `exit_code`, loads), `resource-rank-N.json` (GNU time units), results,
density files, and `tables.json`. Relay all raw files, not just medians. Also
preserve `lscpu`, toolchain/build log, commit/source diff, and file SHA-256 list.

Parser replay, independently of the runner:

```sh
python3 -I scripts/lane_launchio_parse.py /data/sig-dft/lanes/LAUNCHIO/si64-np8/on-{1,2,3,4,5}/results.json --require-external --output /data/sig-dft/lanes/LAUNCHIO/si64-np8/tables.json
```

Output `runs[].per_rank[]` has lane sizes, per-species records, bucket means,
local lifecycle clocks and optional resources. `timer_ranges` reports min/max
per timer without adding them. `external.wall_minus_max_main_s` is a residual
including launch/loader/late append/post-main and staggered-rank effects, not
an isolated launcher measurement. `external.wall_minus_max_scf_s` also includes
constructor/input/output; it is not a sum of separately reduced phases.

## Pre-registered data gates and handoff

- Every retained on launch must have all ranks, consistent run/input/binary IDs,
  completed/converged/verified output, valid density hashes, and the lane timers.
- All ranks and off/on runs must have bit-identical energy/density; iteration
  and verification histories must agree. The runner checks root energy/density/
  iterations; parser checks science across ranks in on runs. The master compares
  full scientific histories using the same timing-field exclusion as `verify_cli.py`.
- Each local bucket wall sum must equal local `launch.covered_s`, and launch
  closure must reconstruct `main_wall_s - covered_s` within 1e-9 s. No sum of
  cross-rank maxima is used as closure. Bad data is returned as a failed gate.
- All five retained repetitions and all sidecars are required; report median,
  range and rank imbalance per component. Zero clocks in an active nontrivial
  span are unusable and need investigation, not an inferred zero cost.
- No speed threshold or expected gain is registered in Phase 0. The counter
  target here is **complete attribution data with zero scientific changes**.
  Candidate proposals, expected gains and the Phase-2 A/B pre-registration
  (including resolved gain |mean|/SE >= 3) await these data.

Stop at PROFILE_REQUEST. No request to adopt or optimize is implicit.
