# Collective arrival-skew diagnostic design

Status: design and offline parser delivered; **the trace producer is not
implemented or runnable on the current binary**. No new switches registered.
Cross-lane instrumentation scope is needed before the full request can run.
The existing PHASE_WAIT pilot below is runnable but cannot answer every-call
attribution. FFT_KEPT_FUSE is closed NEUTRAL, default off.

## What PHASE_WAIT actually records

Source inspected: distributed_davidson.rs `Counted::probe`, `wait_probe`,
`ApplyAgreement`; engine gamma_real.rs solve entry; engine parallel.rs
`SharedComm::timed`; numerics stage_clock.rs; engine profile.rs map lifecycle.

- `Counted` inserts a barrier before solver allgatherv, tree_sum, bcast and
  alltoall. Barrier duration accumulates in solver `wait`; not a timestamp
  trace, not per-operation/per-call attribution. Stage timers off suppress
  the clock **but do not suppress the inserted barriers**.
- M2 entry inserts another barrier, recorded as `entry_wait`, after ensure_state
  and before operator update/start preparation. This measures a different boundary.
- The distributed H operator uses the underlying communicator, not the solver
  Counted wrapper (`ApplyAgreement` explicitly documents this). PHASE_WAIT
  therefore does not insert a barrier before every FFT/nonlocal collective.
- `SharedComm` sees all five logical operations: allgatherv, tree_sum, bcast,
  alltoall, barrier. It currently retains counts/bytes and accumulated duration,
  with per-thread phase mirrors. No per-call entry/exit, rank arrival order,
  map sequence, or adjacent-stage trace is stored.
- `reduce` is a NONLOCAL stage label, usually implemented by tree_sum_f64;
  it is not a separate Collective trait method. Record `op=tree_sum`,
  `site=nonlocal.reduce` once. Do not double count it as both operations or
  assume an MPI_Reduce interceptor captures the engine's deterministic tree.

Thus neither current PHASE_WAIT nor sum/max stage totals can identify which
rank arrives last. A rank that spends the longest time inside a collective
may be the **earliest**, not the latest. Added barriers change progress and
release order; their time includes barrier protocol cost. The comments calling
the remaining collective clock "transfer" are an approximation, not a clean
measurement of production transfer time.

## Proposed producer (separate instrumentation-only implementation)

One default-off `COLLECTIVE_TRACE` switch, valid only with PHASE_PROFILE on;
STAGE_TIMERS off selects count-only, on selects timestamps. PHASE_WAIT must
be off for the main trace. One thread/rank, one communicator, MPI8 initially.
Instrument logical entry/exit once in SharedComm; audit every fine-map caller
for bypasses. Compare operation counts against the existing counters. Include
bcast/barrier as well as the four requested labels, so sequence matching and
preceding intervals are complete. Trace maps from map_begin through map_commit;
assign stable ordinal first, finalize bucket after verification classification.
Keep all fine-map phases (potential, solve, density, audit), not only solve.

No extra collective, barrier, logging, allocation or I/O in the measured loop.
Use fixed-capacity per-rank buffers allocated before SCF, capacity sized from
the count pass; no overwrite/ring eviction. Record overflow and mark the trace
incomplete, retaining scientific results. OFF does not clock, allocate, update
TLS or serialize events. Output after the run through CLI, with a separate
results.json section/sidecar to avoid merging into other lanes' slot tables.

Record stable call-site IDs/scopes (FFT inverse/forward, nonlocal reduce,
solver projection/orthogonalization, density, mixing, etc.), a communicator ID,
map ordinal and collective sequence. Insert explicit site scopes; don't collect
allocating backtraces in hot loops. Existing counters' stage enums may label
intervals, but nested apply/FFT/collective spans must be made **exclusive**
offline (interval splitting). The parser's preceding_stages are disjoint leaf
compute spans between previous collective exit and next entry. Missing coverage
is `unattributed`, never silently assigned to the preceding FFT stage.

Snapshot structural metadata outside loops: shape, first FFT axis, local and
max stick counts, local coefficient/projector row counts, owned plane count,
active z support, pair/column width, matrix leading dimensions/strides, kept
count, logical/padded/useful bytes, projector channel count. Layout records
have IDs shared by spans/events; preserve per-call variable widths. Rank
placement metadata includes hostname, CPU/core/socket/NUMA node and affinity.
For NUMA placement inspect /proc/<pid>/numa_maps outside timed sections or
run a separate sampling job; do not infer remote pages from socket binding.
Do not retune partitioning or memory placement in the diagnostic baseline.

### Clocks and interpretation

std::Instant origins are process-local; never compare raw elapsed values from
independently created origins. MPI_Wtime is not automatically synchronized.
For this one-host request, pair SystemTime before/after an Instant anchor
outside the map, with another anchor after the map. Express timestamps on a
common host epoch, bounded by bracket widths and observed drift; reject clock
steps, cross-host use or inconsistent offsets. Use the smallest bounded anchor
bracket from repeated samples; do not pretend barrier exit is simultaneous.
A subsequent native common monotonic clock adapter is an alternative requiring
its own FFI ownership. Target calibration error below 1 us; report actual
error and mark unresolved rank order rather than manufacturing a winner.

For event j, common-clock entry a_r and exit e_r:

- Latest entry a_max, arrival gap g_r = a_max - a_r, with uncertainty bounds
  max(0,g_r-2 epsilon)..g_r+2 epsilon. Candidate latest ranks have entries
  within 2 epsilon of a_max; give a unique last rank only when resolved.
- Collective residence e_r-a_r; pre-last-entry residence
  min(g_r,e_r-a_r). This is **exposure to late arrival**, not proven idle wait:
  communication can progress while other ranks compute; some collective ranks
  can complete early. An exact idle/transfer split needs MPI progress/message
  tracing or a controlled perturbation. Do not label the remainder pure transfer.
- Compare last vs earliest rank:
  entry_difference = previous_exit_difference + local_interval_difference.
  This separates inherited release skew from new computation skew. Attribute
  local intervals using exclusive stage times and per-stage last-rank minus
  other-rank median excess. The longest stage alone is not necessarily causal.
- Report per event/map/rank and aggregate by site/map. Rank-summed exposure is
  not wall savings. Report first-map-entry skew separately from ongoing release.

The NONLOCAL observation supports waiting as the working hypothesis and
weakens a payload-only explanation. Equal work and 66.6% duration spread do
not identify the producer of that wait. PACKED_REDUCE's neutral 28.1% traffic
reduction is scoped to its measured case, not proof that bytes never matter.

## Attribution tests with returned data

1. Does the same rank repeatedly arrive last? Compare si64 to balanced si32.
2. Before inverse FFT transpose, normalize FFT0 by s_r*n0*pairs, scatter by
   local coefficients, pack by useful bytes. Excess explained by work counts
   suggests stick assignment; equal work but unequal ns/point suggests locality,
   contention, scheduling or placement. Preserve mirror-pair ownership constraints.
3. Before forward transpose, compare planes, FFT1/2 points and kept/fresh widths.
   On 72^3/8 each owns 9 planes; no plane imbalance can be claimed from timing
   alone. Test whether previous release skew rather than local work explains it.
4. Before nonlocal tree_sum, compare projector rows, channels, widths and GEMM
   flops/strides. Determine whether FFT arrival skew propagates to this site or
   new skew starts in projector computation. Interleave both timelines by sequence.
5. Separate memory-bandwidth explanations from NUMA: equal volume + variable
   runtime is not a hardware counter. On a trace survivor, request isolated
   perf/memory-controller data with supported dev08 events and affinity/page
   residency, without simultaneous acceptance timing. No hard-coded PMU events.

Only then select static ownership balancing or removal of a measured skew
producer. Current evidence does not support a specific >=0.3 s pass-removal
candidate; no optional optimization was implemented.

## Proposed JSON and offline parser

`arrival-skew-1` is a **proposed new sidecar**, not today's results schema.
Top-level: schema, mode counts|timed, complete, rank_count, clock_valid,
clock_error_ns, ranks[]. Each rank: rank, maps[]. Each map: map_id, bucket,
start_ns (common epoch), layout (structural metadata), events[]. Each event:
seq (0-based contiguous within map), op, site, payload_bytes, enter_ns, exit_ns,
preceding_stages[]. Each span: stage, start_ns, end_ns, points, logical_bytes,
flops, plus optional layout ID/shape/stride details. In counts mode timestamps
may be null and are not read by the parser. Sequence/site/count identity is
required across ranks; payload bytes may differ legitimately by rank. Exclude
all output/calibration collectives from this application sequence.

`python3 -I docs/performance/lanes/FFT/wait_parse.py trace.json` writes JSON
rows with each event's last_rank_candidates, arrival gaps/bounds, collective
residence, exposure, preceding stage durations/excess and release decomposition.
Counts mode returns identity/payload rows and no invented times/last rank.
The caller must compare counts and metadata between count/timed runs before
interpreting time. Parser rejects incomplete trace, uncalibrated clocks,
missing ranks, map/sequence mismatches and overlapping compute spans.

`python3 -I docs/performance/lanes/FFT/test_wait_parse.py`: three synthetic
contract tests pass on Mac, including known last rank, uncertainty ties,
early completion, count-only and corrupt traces. **No actual CLI-generated
arrival trace or producer equivalence validation is claimed.**
