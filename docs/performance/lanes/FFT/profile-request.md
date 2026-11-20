# LANE-FFT Phase 0: dev08 profile request

Status: PROFILE_REQUEST. Phase-0 instrumentation is wired and locally verified.
Scope extension approved by the orchestrator on 2026-11-20.
This supersedes the original optimization charter for this run. No optimization,
Mac performance conclusion, acceptance threshold, or default change is proposed.
The prior exploratory code and Mac timing interpretation have been removed.

## Question for dev08

In **ordinary fine maps**, how much local-potential time is spent in each FFT
axis, local permutation/zeroing, all-to-all, and pointwise multiplication? How
much is fixed per pair/call, and how much follows local grid points or padded
collective volume? Are rank differences computational or collective waiting?
Return the profile before any candidate is designed (Phase 2).

## Instrumentation contract

`BAANDI_FFT_PROFILE=on` AND `BAANDI_PHASE_PROFILE=on` enables this lane.
`BAANDI_STAGE_TIMERS=off` gives counts without clocks. All defaults unchanged.
Switch values are captured at plan construction. Disabled sites do no clock,
TLS-counter, allocation or I/O work (there is a cached branch, not literally
zero machine instructions). Enabled timers perturb small FFT stages: profiles
are diagnostic, never acceptance timing. No numerical operation is reordered.

Stages: scatter, inverse_fft0, transpose_pack, inverse_alltoall, middle_unpack,
inverse_fft1, slab_scatter, inverse_fft2, multiply, forward_fft2, slab_gather,
forward_fft1, middle_pack, forward_alltoall, transpose_unpack, forward_fft0,
gather_scale, kept_copy. FFT axes refer to the plan's chosen first/second/z axes, not
necessarily Cartesian x/y/z. Zero-sized lanes are not transformed.

Per stage: seconds, calls, complex points, logical bytes, FFT operation estimate
`round(5 * points * log2(length))`. This is a conventional complex FFT model,
not an instruction count. Butterfly internal traffic is unknown: FFT bytes
count only a read/write of the lane buffer, a lower bound. Permutation bytes
count reads plus writes, including explicit zero fill; index metadata and
write-allocate traffic are omitted. All-to-all bytes are **outgoing payload per
rank**, including padding and the self block, not bidirectional memory traffic.
Multiply: 40 logical bytes and 2 real multiplies per complex point. Gather/scale:
32 bytes and 2 real divisions per coefficient. Do not sum nested pre-existing
H/FFT timers with these stage times. Stage clocks exclude validation, batch-buffer reservation, and the accounting
between clocks; their sum is not a complete H-apply wall. Retain existing
`convergence.work.davidson_stages.by_run_stage.fine.fft` fields
`pair_unsplit_s` (M2), `pair_inverse_s`, `pair_potential_s`, `pair_forward_s`
and `kept_forward_s` for enclosing context, and
`phase_profile.buckets` entries `solve.apply`, `h_columns`, `h_delta_columns`.
Those existing maxima must not be mixed with one rank's lane times as a
closed decomposition. No hardware bandwidth is assumed. The probe counts only local-potential
applications (including kept-grid potential changes); density and Hartree FFTs
are excluded. Serial kept-grid multiplication includes its grid read and slab
write in `multiply`, while M2 records the separate plane copy as `kept_copy`.

This request supports **one thread per rank only**. TLS counters on worker
threads are not merged into the caller's profile for T>1; such a run cannot
support a complete lane table. One-core runs use threads=1; MPI jobs use 8
ranks with threads=1 in the generated input. No threading experiment is requested.

Verified JSON contract:
`phase_profile.fft_lane` has schema `lane-fft-1`, enabled, stage_timers,
and `per_rank[]`. Each rank has rank and buckets[]; each bucket has bucket,
maps, stages[]; each stage has stage, seconds, calls, points, logical_bytes,
fft_flops_estimate. Exact paths (arrays must be selected by label/rank):

- `phase_profile.fft_lane.schema` = `lane-fft-1`
- `phase_profile.fft_lane.enabled` and `.stage_timers`
- `phase_profile.fft_lane.per_rank[r].rank`
- `phase_profile.fft_lane.per_rank[r].buckets[b].bucket` and `.maps`
- `phase_profile.fft_lane.per_rank[r].buckets[b].stages[s].stage`
- the same stage object's `.seconds`, `.calls`, `.points`, `.logical_bytes`,
  `.fft_flops_estimate`
- `basis.fft_grid`, `basis.n_plane_waves`, `bands`, `threading.threads` give
  full-cutoff sizing/context. Coarse sizes differ: use that bucket's own
  points/work totals; do not attach the full-cutoff grid size to coarse data.

All FFT integer metrics survive the final allgatherv exactly (two 32-bit limbs
per u64); the added transfer is only in enabled profiling runs. FFT records
are a separate section, not added to exclusive wall closure or shared stage
slot lists. Per-map figures are same-rank bucket totals divided by maps, **not
individual map traces**. Keep fine, fine_redo, coarse, verify_pass and
verify_fail separate. Never subtract or sum different ranks' maxima.

## Build and provenance (orchestrator only, no remote action by implementer)

At the reviewed commit, record `git rev-parse HEAD`, `rustc -Vv`, build flags,
MKL version, MPI version, input/UPF/binary SHA-256, CPU binding and host load.
Use separate directories for the serial and MPI builds. Do not change the FFT
backend from the currently approved default. `PAIR_FFT=grid` fallback is not
covered by detailed serial stick-lane instrumentation; flag it as unavailable.

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export FFT_WORK="$PWD/target/lane-fft-profile"
mkdir -p "$FFT_WORK"
CARGO_TARGET_DIR=target/lane-fft-serial cargo build --release -p dft-cli --features mkl
CARGO_TARGET_DIR=target/lane-fft-mpi cargo build --release -p dft-cli --features mpi,mkl
sha256sum target/lane-fft-{serial,mpi}/release/dft-cli "$DFT_SI_UPF"
```

`DFT_SI_UPF` must name the same Si UPF as the orchestrator's matched baseline.
One case/model per queue job, timeout 600 seconds. Prepare inputs using the
existing seeded realistic-displacement generator (seed 601, amplitude 1e-3),
which sets threads=1. Do not compare MPI with one core as a speed ratio.

```sh
# Set CASE=si64, then submit si32 and si16 as separate MPI jobs.
CASE=si64
python3 -I scripts/fast_screen_dav.py --case "$CASE" --upf "$DFT_SI_UPF" \
  --threads 1 --seed 601 --amp 0.001 --emit-input "$FFT_WORK/$CASE.toml"
for REP in 1 2 3; do
  env BAANDI_FFT_PROFILE=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
    BAANDI_DAVIDSON_STAGES=on \
    mpirun -np 8 --map-by ppr:8:socket --bind-to core --report-bindings \
    target/lane-fft-mpi/release/dft-cli run "$FFT_WORK/$CASE.toml" \
    --output-dir "$FFT_WORK/$CASE-np8-$REP"
done
# A count-only control: all deterministic FFT counters must equal timed runs.
env BAANDI_FFT_PROFILE=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=off \
  BAANDI_DAVIDSON_STAGES=off \
  mpirun -np 8 --map-by ppr:8:socket --bind-to core --report-bindings \
  target/lane-fft-mpi/release/dft-cli run "$FFT_WORK/$CASE.toml" \
  --output-dir "$FFT_WORK/$CASE-np8-counts"
```

For each of si32 and si64, a separate one-core job, same input generator:

```sh
CASE=si32
python3 -I scripts/fast_screen_dav.py --case "$CASE" --upf "$DFT_SI_UPF" \
  --threads 1 --seed 601 --amp 0.001 --emit-input "$FFT_WORK/$CASE.toml"
for REP in 1 2 3; do
  env BAANDI_FFT_PROFILE=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
    BAANDI_DAVIDSON_STAGES=on \
    taskset -c 0 target/lane-fft-serial/release/dft-cli run "$FFT_WORK/$CASE.toml" \
    --output-dir "$FFT_WORK/$CASE-onecore-$REP"
done
```

Controls (execute the appropriate model's commands inside its job):

```sh
# MPI OFF equivalence control, same binary/input/binding:
env BAANDI_FFT_PROFILE=off BAANDI_PHASE_PROFILE=off BAANDI_STAGE_TIMERS=off \
  BAANDI_DAVIDSON_STAGES=off \
  mpirun -np 8 --map-by ppr:8:socket --bind-to core --report-bindings \
  target/lane-fft-mpi/release/dft-cli run "$FFT_WORK/$CASE.toml" \
  --output-dir "$FFT_WORK/$CASE-np8-off"
python3 -I docs/performance/lanes/FFT/verify_cli.py \
  "$FFT_WORK/$CASE-np8-off" "$FFT_WORK/$CASE-np8-1" "$FFT_WORK/$CASE-np8-counts"
```

One-core controls (one-core jobs only):

```sh
env BAANDI_FFT_PROFILE=off BAANDI_PHASE_PROFILE=off BAANDI_STAGE_TIMERS=off \
  BAANDI_DAVIDSON_STAGES=off \
  taskset -c 0 target/lane-fft-serial/release/dft-cli run "$FFT_WORK/$CASE.toml" \
  --output-dir "$FFT_WORK/$CASE-onecore-off"
env BAANDI_FFT_PROFILE=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=off \
  BAANDI_DAVIDSON_STAGES=off \
  taskset -c 0 target/lane-fft-serial/release/dft-cli run "$FFT_WORK/$CASE.toml" \
  --output-dir "$FFT_WORK/$CASE-onecore-counts"
python3 -I docs/performance/lanes/FFT/verify_cli.py \
  "$FFT_WORK/$CASE-onecore-off" "$FFT_WORK/$CASE-onecore-1" \
  "$FFT_WORK/$CASE-onecore-counts"
```

Return energy bits, density.bin SHA-256, eigenvalues, convergence trajectory
and work counts; require bit identity, not just an energy tolerance. Timed
values are not expected to match. Report any absent stage as unavailable,
never as zero cost. Profiles include verification but tables separate it.

## Return bundle and parser

Return all results.json files and input/UPF/binary hashes, actual commands and
exit codes, density hashes, convergence outcomes, load and CPU placement.
Keep each repetition separate. No QE run is requested in Phase 0.

```sh
# Inside an MPI case job, after the three runs and count-only control:
python3 -I scripts/lane_fft_parse.py "$FFT_WORK/$CASE-np8-"[123]/results.json \
  "$FFT_WORK/$CASE-np8-counts/results.json" > "$FFT_WORK/$CASE-np8-lane.csv"
# Inside a one-core case job:
python3 -I scripts/lane_fft_parse.py --json "$FFT_WORK/$CASE-onecore-"[123]/results.json \
  > "$FFT_WORK/$CASE-onecore-lane.json"
```

Parse enabled runs only (OFF controls intentionally reject as unavailable).
Review fine-bucket tables by rank and stage, with seconds/call, seconds/map,
counts and floor inputs. Zero calls are labeled unusable. Count-only runs have
null time-derived quantities. Repetition medians/ranges should be computed
within the same case/model/rank/stage, not across models or machines.


## Local verification (correctness only)

[cli-equivalence.json](cli-equivalence.json) records a real Si8 CLI run in OFF,
timed and count-only modes: 11 fine SCF iterations, verified in all arms,
energy difference 0 Ha, identical energy/eigenvalue/occupation/history bits
and density SHA-256. Only `timings` and verification `wall_s` are excluded
from scientific comparisons. Both enabled files parse to 144 rows with equal
integer counters. [mac-counts-results.json](mac-counts-results.json) is the
actual unmodified count-only CLI results.json, retained for parser regression;
its ordinary CLI wall fields must not be used for performance conclusions.

```sh
python3 -I scripts/lane_fft_parse.py docs/performance/lanes/FFT/mac-counts-results.json
python3 -I docs/performance/lanes/FFT/test_parser.py
```

For another correctness check, run the three arms on the same input/binary and:

```sh
python3 -I docs/performance/lanes/FFT/verify_cli.py OFF_DIR TIMED_DIR COUNTS_DIR
```

After the dev08 profile is returned, stop for the orchestrator's Phase-1
analysis. No candidate or optimization is proposed in this Phase-0 request.
