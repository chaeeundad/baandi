# LANE-FFT: collective arrival profile request

Status: **READY FOR DEV08 PROFILE**. Instrumentation-only scope approved;
producer, actual CLI JSON wiring and parser implemented. FFT_KEPT_FUSE remains
OFF/CLOSED NEUTRAL. No optimization or acceptance timing is requested.

## Capture contract

BAANDI_COLLECTIVE_TRACE=on with PHASE_PROFILE=on enables a preallocated
131072-event buffer per rank (18 MiB). STAGE_TIMERS=off records counts/labels/
shapes without trace clocks; on records timestamps. FFT_PROFILE=on enables
FFT preceding-stage detail. PHASE_WAIT=off is mandatory: no inserted barriers.
One thread/rank, single host only. Every SharedComm operation during a pending
map is captured: alltoall, tree_sum, allgatherv, bcast, barrier. NONLOCAL reduce
is op=tree_sum/site=nonlocal.reduce, not a separate MPI operation.

SCF numerical operations and application collective sequence/counts are
unchanged. Two diagnostic allgatherv calls (lengths, limb-encoded events) occur
**after SCF**, only when enabled. Disabled capture has only cached branches,
no trace clocks, allocation, buffer updates or extra communication. Profile
initialization allocates the buffer; no capture-loop allocation, logging or
I/O. Overflow preserves scientific results but sets complete=false: reject the
trace, never use a truncated prefix. Output/gather costs are not lane speed.

Stamps use bracketed SystemTime/Instant anchors, translated after gathering to
nanoseconds relative to the earliest rank's trace initialization (common trace
run start, not launcher start). Single-host clock-domain assumption is required.
Post-run anchor drift is added to uncertainty; >1 ms marks clock invalid.
Clock steps are invalid; two anchors cannot detect a step that later undoes
itself. Save host clock-status information. Parser leaves latest-rank ties
unresolved within twice the reported error. No cross-host clocks supported.

Previous-stage labels are the **last completed instrumented stage** since the
previous collective. Stage durations can include nested computation/communication;
they are explicitly not an exclusive decomposition. A missing stage is
unattributed; semantic sites still label solver reductions and FFT/NONLOCAL
communication. Do not invent a preceding kernel or divide an enclosing
timer by a leaf's work. FFT stage shape is chosen [n0,n1,n2]; nonlocal projector
shape is [local_rows,projector_channels,columns]. Its logical bytes model is
8*(rows*channels + rows*columns + channels*columns), not measured DRAM traffic.
Other/unknown shapes are [0,0,0]. See wait-diagnostic.md for interpretation.

## Build (outside the profile queue jobs)

Use the reviewed tree including this patch. No concurrent builds during
capture. Set DFT_SI_UPF to the exact external file used in the baseline.

~~~sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export FFT_TARGET="$PWD/target/lane-arrival-mpi"
CARGO_TARGET_DIR="$FFT_TARGET" cargo build --release -p dft-cli --features mpi,mkl
export FFT_BIN="$FFT_TARGET/release/dft-cli"
git rev-parse HEAD
sha256sum "$FFT_BIN" "$DFT_SI_UPF"
~~~

## Two separate jobs: si64 np8 and si32 np8

Run the following blocks once with CASE=si64 and once with CASE=si32, separate
queue jobs/windows, each capped at 600 s. Baseline switches remain identical;
record full manifests and keep unadopted NONLOCAL experiments at the baseline
setting. No guessed NONLOCAL variable names. No one-core or cross-model ratios.

~~~sh
set -euo pipefail
export CASE=si64
# For the second job: export CASE=si32
export WAIT_WORK="$PWD/target/arrival-$CASE"
mkdir -p "$WAIT_WORK"
python3 -I scripts/fast_screen_dav.py --case "$CASE" --upf "$DFT_SI_UPF" \
  --threads 1 --seed 601 --amp 0.001 --emit-input "$WAIT_WORK/input.toml"
sha256sum "$FFT_BIN" "$DFT_SI_UPF" "$WAIT_WORK/input.toml" > "$WAIT_WORK/hashes.txt"
hostname > "$WAIT_WORK/host.txt"
lscpu -e=CPU,CORE,SOCKET,NODE > "$WAIT_WORK/topology.txt"
uptime > "$WAIT_WORK/load-before.txt"
~~~

Run one OFF scientific control, then C,T,C,T,C,T (three count-only controls
and three timed repetitions). Same binary/input, CPUs 0–7 socket 0. Single-host
launcher is explicit; inspect binding reports and reject any other placement.

~~~sh
env MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  BAANDI_COLLECTIVE_TRACE=off BAANDI_PHASE_PROFILE=on BAANDI_FFT_PROFILE=on \
  BAANDI_PHASE_WAIT=off BAANDI_STAGE_TIMERS=off BAANDI_DAVIDSON_STAGES=on \
  BAANDI_FFT_KEPT_FUSE=off \
  mpirun --host localhost -np 8 --map-by ppr:8:socket --bind-to core --report-bindings \
  "$FFT_BIN" run "$WAIT_WORK/input.toml" --output-dir "$WAIT_WORK/off" \
  > "$WAIT_WORK/off.log" 2>&1
for REP in 1 2 3; do
  for TIMER in off on; do
    env MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
      BAANDI_COLLECTIVE_TRACE=on BAANDI_PHASE_PROFILE=on BAANDI_FFT_PROFILE=on \
      BAANDI_PHASE_WAIT=off BAANDI_STAGE_TIMERS="$TIMER" BAANDI_DAVIDSON_STAGES=on \
      BAANDI_FFT_KEPT_FUSE=off \
      mpirun --host localhost -np 8 --map-by ppr:8:socket --bind-to core --report-bindings \
      "$FFT_BIN" run "$WAIT_WORK/input.toml" --output-dir "$WAIT_WORK/$REP-$TIMER" \
      > "$WAIT_WORK/$REP-$TIMER.log" 2>&1
    python3 -I docs/performance/lanes/FFT/wait_parse.py \
      "$WAIT_WORK/$REP-$TIMER/results.json" > "$WAIT_WORK/$REP-$TIMER/arrival-table.json"
  done
done
uptime > "$WAIT_WORK/load-after.txt"
~~~

Run under a 600 s queue timeout, recording every exit. All seven SCFs must
converge and verify; no fixed-step non-convergence exception here. If
serialization/transport breaches the cap, retain failure and size a smaller
capture window through an explicit follow-up; do not truncate traces or
silently change the input.

## Real JSON paths (implemented)

results.json.phase_profile.arrival_trace (also emitted by CLI's DAVBENCH_JSON
profile serializer) contains:

- schema="arrival-skew-2", enabled, mode="counts"|"timed", rank_count,
  capacity_per_rank, complete, clock_domain, clock_valid, clock_error_ns.
- ranks[]: rank, clock_error_ns, clock_drift_ns, events[].
- events[]: map_id, bucket, seq (within map), call_index (within rank trace),
  op, site, payload_bytes, enter_ns, exit_ns, previous_exit_ns.
- events[]: previous_stage, previous_stage_ns, previous_stage_end_ns,
  previous_stage_points, previous_stage_logical_bytes, previous_stage_shape,
  stage_duration_is_exclusive=false.

Timed fields are JSON null in count-only mode. Disabled/no-event trace is
null, not zero-cost evidence. Bucket is finalized using the existing profile
map classification; parser selects ordinary fine only. Capture also retains
other maps for audit. Legacy profile bucket counters may classify lazy
materialization separately, so do not demand equality with a different
bucket's aggregate counter (Mac warm solve: trace fine 408 vs legacy fine 391).
Compare matching per-map event signatures and existing work across arms;
wrapper coverage is separately unit-tested for all five operations.

## Gates and return tables

Before interpretation, OFF/count/timed energy/eigenvalues/occupations, density
bytes and scientific histories/work must be bit-identical at fixed rank count.
Return verification outcome, hashes and work counts. Across C/T compare exact
(map_id,bucket,seq,call_index,op,site,payload_bytes,previous_stage,points,bytes,shape)
per rank. Zero differences required; no tolerance relaxation. Parser rejects
missing ranks, overflow, invalid clocks/intervals, mismatched collective
sequences and missing fine events. Count-only gives no invented last rank.

Return input/results.json, logs/exits, input/UPF/binary/source hashes, binding,
host/load/topology and parsed tables. Per event/map/site: arrival_gap_ns with
bounds, last_rank_candidates (unique only when resolved),
collective_residence_ns, pre_last_entry_residence_ns, preceding_interval_ns,
previous-stage context/work, and last-versus-earliest previous-release/local-
interval decomposition. Aggregate latest-rank frequencies and median/p95 skew
by site/map, retaining raw per-event rows. No rank-summed exposure as wall saving.

Arrival gap is not measured idle time: communication can progress before the
last entry, and some ranks can complete early. Pre-last-entry residence is
exposure, not pure wait; the remainder is not pure transfer. Identify repeated
slow-rank/work-shape patterns before a candidate. Test whether skew propagates
from FFT to nonlocal or starts at a new stage. Equal points but different time
does not distinguish NUMA, memory contention and scheduling: request isolated
hardware/page-residency evidence if that ambiguity matters. No static balancing
or additional permutation optimization is authorized by a trace alone.

## Mac validation

See arrival-validation.md and arrival-equivalence.json. Actual CLI profile
serializer on real-UPF Si8 simulated MPI3 produced 408 events/rank, parsed to
408 fine-event rows; count/timed signatures identical. Actual trace fixture:
mac-arrival-trace.json. No Mac timing supports performance conclusions.
