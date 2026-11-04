# Runtime equivalence acceptance protocol

Status: acceptance protocol with the executed P11E budget-recovery slice below; the broader mode/performance matrix remains future work. The P11R docs-only check does not establish runtime equivalence. Grounding: `ScfProblem::run/run_outcome/run_impl`, `step_kpoints`, `IterationRecord`, `FixedPointCertificate`; existing `crates/dft-engine/tests/{scf,kpoints,threading,pbe}.rs` and `docs/performance/baseline*.json`.

## Controlled comparisons

Pin source/build/backend/CPU, requested and effective threads, seed, synthetic generator revision, model and input hashes. Keep an independently runnable pre-extraction driver during P11A; comparing two wrappers over the same refactored implementation is insufficient. For each completion mode compare old driver, session with no observer, session with snapshots/events, session paused/resumed at each boundary, and exact-checkpoint resumed session. Match operation ordering and immutable model. Normalize only timings, timestamps and external identity fields; never round numerical values.

Require bit-identical total and component energies, density elements, per-k eigenpairs/occupations and numerical history including requested tolerances, operator counts, mixing fraction/reset and certificates. Preserve error variants, history length and failure stage for invalid input, divergence, eigensolver failure, noninsulating result and exhausted budget. Explicitly cover the last-step `strict` mixing versus `retain_last` early-return distinction and failed-certificate continuation. First iteration has null energy change. Candidate convergence is not enough; certificate and gap gate must be identical.

Matrix: local and multiple-projector synthetic potentials, cubic/skew/elongated/slab-like cells, Gamma and multiple canonical k points, small and large band counts, small and large insulating gaps. Reuse existing fixtures where available and add separately ticketed fixtures for missing cells/gaps; do not describe missing coverage as passed. The real Si case is an additional local regression, not general validation. Pin P10A/P10D/P10C/P7PC merged baseline before extraction. Keep existing cross-thread tests; same-thread replay is required at thread counts 1, 2 and 4 where supported. Cross-platform equality is not promised.

## Interventions and persistence

Replay a recorded accepted/rejected policy log at exact `(iteration_done,revision)` boundaries; compare against the original controlled trajectory bitwise, excluding audit actor/timestamps. Do not demand bit equality between different policies. Verify changed-policy final energies against a separately fully converged same-model reference within 1e-10 Ha, with unchanged or tighter certificate tolerances and density/charge/eigen residuals reported. This is a demo acceptance threshold, not evidence that all recovered states must represent the same physical branch.

Checkpoint each Ready/Paused boundary, resume and compare the complete subsequent trajectory. Corrupt arrays, omit mixer history, change model/backend/thread count or completion mode: exact load must reject without mutation. Warm starts reset history/certificate and get new ids; compare final quality, never claim trajectory identity. Exported unconverged results must remain visibly unconverged. Event polling frequency and replay duplicates must not change numerics.

Record full command, exit code, versions, hashes, input, baseline id, thread topology, comparisons and failures in implementation-owned evidence. Operator total is sum of ordinary `operator_applications` plus nested certificate applications, with no double counting. Include wall time separately and observer overhead; diagnostics must not trigger hidden eigensolves. Run `scripts/check.sh`; >10 s debug cases become ignored release tests with documented commands. Use `scripts/bench_scf.py` against the selected baseline as required by performance policy (energy within 1e-12 Ha and bit-identical density for overhead-only changes). This document creates no new measured baseline.

## Three-mode overhead targets and P11J measurement plan

Owner of future measurements: P11J-RUNTIME-BENCH (ticket on main `98f0b50`). P11R owns this design section; P11J appends actual results and evidence under `tests/reference/runtime/` after P11A/P11B merge. All numbers in the target column are acceptance goals, not measured performance.

| Comparison | Target relative overhead | Gate |
| --- | --- | --- |
| Classic vs pinned pre-refactor path | <0.5% | No observer/queue/telemetry serialization; numerical outputs identical |
| Observed vs same pre-refactor path | <1% | Standard bounded scalar sink, no intervention; outputs identical |
| AgentControlled vs same pre-refactor path | <1% | Control channel configured, no intervention/LLM wait; outputs identical |
| Intervention, checkpoint, slow consumer | Report separately; no invented percentage target | Preserve safety/identity promises; explicitly account for intentional pauses, algorithm changes and I/O |

Fix input, executable build profile/compiler/dependencies, hardware/backend, requested/effective thread counts, CPU affinity where available, initial state/seed, completion mode and storage device per pair. Record a runnable pre-refactor executable/revision and paired post-refactor executable, not merely a stale historical timing value from a different machine. Use `docs/performance/baseline*.json` for traceability and scientific comparison, then rerun the reference path on the same host. Repeat source hashes for P10A/D/P7PC merged features so a changed solver is not mistaken for observer overhead.

Run at least 3 warm-ups and 20 measured paired runs per mode/case on a quiet host; alternate/interleave reference and candidate order to expose drift. Report every sample, sample count, median, p5/p95 and paired overhead distribution; define paired overhead_i=100*(T_mode_i/T_reference_i-1). Use a recorded bootstrap seed and 10000 paired resamples to report a 95% confidence interval for median overhead. Mark target met only when the interval upper bound is below the target and equivalence gates pass; overlapping intervals are inconclusive, not rounded into a pass. Run count may increase under a predeclared noise rule, never cherry-pick quiet samples. Report background load and exclusions with reasons. Claude re-measures accepted performance on quiet dev08 through dft5090, as P11J requires; P11R does not alter remote servers.

Include complete normal converged calculations and deliberately short SCF workloads, not only long runs where setup hides observer cost. Matrix includes the tiny synthetic fixture, representative local/multiple-projector problems, Gamma and k meshes, different band counts, skew/elongated/slab-like cells and small/large gaps from P10G where available. Declare missing fixture coverage. A one-iteration exhausted case measures short-path overhead and is labelled unconverged, not a physical benchmark result. Reuse frozen inputs across all modes and count certificate work identically.

| Metric | Required measurement and units |
| --- | --- |
| Total time | External monotonic wall seconds, process CPU seconds, setup/SCF/teardown split; agent reasoning time separate |
| Per-iteration overhead | Boundary cost nanoseconds and total SCF time/ordinary iteration; separate iterations with certificate from those without |
| Event cost | Typed construction/enqueue ns/event, serialization ns/event and delivery p50/p95 latency ns; event type/bytes/count, sink and subscriber count |
| Memory | Peak RSS bytes and observer allocation counts/bytes; queue high-water marks and journal bytes; subtract matched reference but report absolute values too |
| Checkpoint cost | Capture pause seconds, async persistence seconds, payload/manifest bytes and peak extra memory; exact vs warm vs export separately |
| Control latency | Receipt->validation, approval wait, accepted->applied, and receipt->applied nanoseconds plus boundary count; stale/rejected and applied counts |
| Saturation | Producer rate, consumer rate, drop/coalesce counts, missing sequence intervals, critical retained/replayed count, disconnect count and fault boundary |
| Reproducibility | Bitwise density/energy/history/certificate comparisons across modes and replay; final energy in Ha for intentional intervention/reference comparison |

Instrumentation belongs in the benchmark/driver, not new timers/logging in numerical kernels. First measure timer/instrumentation overhead and timer resolution; report results below resolution as unresolved. A disabled-event allocation trace must show no queue/event formatting in Classic. Per-event microbenchmarks complement, never replace, matched end-to-end runs. Async serialization/journal drain is included through final flush in normal end-to-end wall timing; report SCF-only latency separately so deferred I/O cannot conceal overhead.

## Mode, stale-control and saturation regression matrix

- Compare Classic/Observed/AgentControlled-no-intervention to the pinned path at identical threads; snapshots and subscriber attach/detach must leave numerical output unchanged.
- Exercise every mode-switch table edge, including Classic upgrade rejection, approved Observed<->AgentControlled and terminal/downgrade rejection. For permitted switches compare next-input density/orbitals/mixer state and unchanged numerical policy; verify state-version and audit behavior.
- Inject oscillation-resolved, persists, expired, changed-revision and terminal replies using a scripted agent and virtual clock. Assert exact accepted/rejected/applied event order and no rejected-state mutation. Delay replies across 1, 10 and more-than-total-run iteration counts without introducing an implicit pause.
- Rate-limit agent absence/invalid responses; compare continue_current with no agent to the same trajectory with a silent agent. Replay local rules at the same boundaries independent of transport timestamps.
- Saturate routine queues with a stalled subscriber: verify bounded memory, latest coalesced snapshot, accurate `events_dropped`, unaffected numerics and retained certificates. Disconnect/reconnect slow critical consumers and replay all critical records exactly once logically by seq.
- Exhaust/fail the critical journal: assert no next step starts without credits, explicit telemetry/resource terminal error, reserved final events retained, and no false numerical convergence. Lossless infinite progress is not an expected result under finite storage.
- Crash between audit intent/commit/application and between checkpoint payload/manifest writes; verify exactly-once accepted commits or discarded uncommitted intent, never an unaudited policy mutation. Atomicity fault injection belongs to runtime implementation tests, not this docs task.
- Requests for unmerged Kerker/hybrid or incompatible strict/buffer policies reject. After merge, test compatible selection using P10D/P10A/P7PC's own numerical tests, plus state reset/replay and fixed-point certification.

P11J overhead measurement status: NOT RUN. No target is marked achieved by P11R's repository build/test/fmt/clippy gates, which check existing code rather than an unimplemented runtime.

## P11E measured budget-recovery equivalence (2026-11-04)

Command: `scripts/demo_self_healing.sh`. Evidence:
`tests/reference/runtime/generated/results.json` and each arm's event/audit,
array manifest and normalized complete result. The two-ion synthetic Gamma
fixture uses the same build, seed 2026, one thread and RetainLast in
AgentControlled for all runs. Each of baseline, rules, recorded agent and fresh
1000-step reference is run twice. All initial density hashes agree.

Baseline exhausts after one ordinary solve (43 operator applications) without
a certificate. Both recovery arms converge in 30 ordinary iterations with
three certificate solves and 951 total operator applications. Their energy is
-4.684561184277879 Ha: absolute error 0 Ha versus the fully converged reference,
within 1e-10 Ha. Every final density/orbital IEEE bit, eigenvalue, complete
normalized history and certificate agrees across repeat runs and between
rules, agent and reference. Timers are excluded from identity. The rejected
beta=1.1 request in the agent arm has no effect on this trajectory.

Final certificate density residual: 4.551166513318054e-9 Bohr^-3/2;
occupied eigenvalue change: 2.7803310209684184e-8 Ha; relative charge error:
4.440892098500626e-16. There are zero dropped events. Operator totals include
certificate work once. Retained baseline energy is not a certified reference
comparison, and its certified energy/error fields are null.

This is budget repair at the initial safe point, not a demonstrated cure for
physical instability. No hot path or default is changed. Cross-mode,
pre-refactor overhead, cross-thread/platform, large-band and physical recovery
claims remain outside this result; no P11J performance target is marked met.
