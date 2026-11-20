# Collective arrival-skew diagnostic: implemented contract

Status: producer complete after the orchestrator's instrumentation-only scope
approval. Follow wait-profile-request.md for current runnable dev08 commands.
The previous arrival-skew-1 proposal is superseded by arrival-skew-2; the parser
retains its synthetic v1 tests, not a claim that v1 was deployed.

## Why PHASE_WAIT is insufficient

distributed_davidson::Counted inserts a barrier before solver allgatherv,
tree_sum, bcast and alltoall. The wait stage accumulates duration; it does not
retain call timestamps or last-rank identity. M2 solve entry also inserts a
barrier (solve.entry_wait). H uses underlying SharedComm, bypassing Counted,
so FFT/nonlocal collectives are not all barrier-probed. STAGE_TIMERS=off
suppresses clocks but PHASE_WAIT still inserts barriers. Those probes alter
arrival/release patterns. Neither barrier duration nor the remaining collective
clock cleanly splits unperturbed production wait versus transfer.

## Capture and scope

SharedComm now records all five logical collective kinds during each pending
map. No new communication in SCF, no unsafe, no numerical/reduction changes.
A separate numerics stage_clock::arrival module owns fixed-capacity TLS storage.
Each event records map/call indices, entry/exit, previous boundary, latest
completed stage with work/shape, and optional semantic site. Default off costs
cached branches only; no trace clocks, buffers or communication. Count-only
reads no trace clocks. Allocation is at profile initialization, storage is
bounded, and overflow invalidates the trace. The application collectives
continue even when overflow occurs.

Map lifecycle uses the existing engine profile classifier, including verification
pass/fail/redo. Trace records pending-map intervals; legacy stage aggregates
can separately reattribute lazy materialization. Keep denominators distinct.
Wrapper audit: gamma_real uses SharedComm::active; M2 phase consumers use
m2's SharedComm accessor; solver/H receive that wrapper. Direct simulator/World
calls outside the engine are not claimed covered.

Semantic sites: both FFT transpose directions, nonlocal.reduce and all eleven
existing solver reduction sites. Other calls retain op/sequence and the last
completed stage, with site=unattributed. No allocating backtrace in hot loops.
The last stage resets on collective exit to avoid stale attribution. Parent
timers are inclusive; stage_duration_is_exclusive=false is intentional.
No claim of full exclusive stage-span reconstruction. Unknown work is
zero/unknown, not an assertion of zero work. FFT shape gives chosen axes;
projector shape gives local rows/channels/columns.

After SCF, two diagnostic gathers transport lengths/events losslessly as
u64 low/high limbs through f64. They are absent when off and are not counted
as SCF work. CLI writes phase_profile.arrival_trace. Separate fields preserve
existing slot layouts.

## Clocks and attribution

Single-host SystemTime epoch brackets an Instant anchor. Monotonic elapsed
stamps become relative to the minimum rank anchor at output. Post-run anchor
disagreement increases uncertainty; drift >1 ms invalidates the trace.
Common-host epoch is a premise, not an MPI_Wtime synchronization claim. Two
anchors cannot detect an intervening clock step that later cancels; reject
runs with host clock corrections. Report uncertainty, not forced winners.

For a_r=entry, e_r=exit, a_max=max entry, epsilon=clock error:
arrival gap=a_max-a_r, bounded by +/-2 epsilon (lower clipped to zero).
Candidate last ranks lie within 2 epsilon of a_max. Residence=e_r-a_r;
pre-last-entry residence=min(arrival gap,residence), an exposure metric only:
communication may progress while ranks arrive, and early completion is legal.
Exact idle wait requires deeper MPI progress tracing.

For resolved last and earliest rank, entry skew equals previous-boundary skew
plus the difference of local inter-collective intervals. This separates
inherited release skew from newly generated skew. Compare preceding stage
labels, shapes and counts; inclusive durations are context, not additive cost.
First event uses map entry as previous boundary.

At FFT inverse compare sticks/FFT0 point counts and coefficient volume; at
forward compare planes/middle/FFT2 work (72/8 gives nine planes/rank).
At nonlocal reduce compare projector rows/channels/widths. Look for consistent
last-rank identity and cross-site propagation. Equal work but differing time
does not prove NUMA or bandwidth contention: correlate binding/page placement
and, if needed, separate hardware counters. Static stick reassignment trades
off FFT balance and projector/coefficient balance and must retain mirror pairs.

## Prior optimization closure

FFT_KEPT_FUSE is CLOSED NEUTRAL, default off: exact copy removal did not resolve
a speed gain. NONLOCAL's 66.6% reduce spread and neutral 28.1% payload reduction
support arrival attribution as the next investigation. They do not determine
which rank/stage generated wait or establish volume irrelevance everywhere.
No specific permutation removal supports >=0.3 s critical-path gain; none is
implemented. Diagnostic observations select a later candidate, not adoption.
