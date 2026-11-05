# Runtime equivalence acceptance protocol

Status update (2026-11-05): P11A/B/C/E/J implementations and local records are
merged. Later dev08 results are relayed in
`tests/reference/perf/dev08_final_relayed.md`, Part 3: the 300-pair campaign reports
all 18 combinations meeting targets, provisionally. Raw remote files are still
missing; the local-only failed/inconclusive measurements below remain valid
historical evidence, not superseded by an independently verified remote run.


Status: acceptance protocol with executed P11E recovery and P11J local benchmark slices below; broader performance acceptance remains pending. The P11R docs-only check does not establish runtime equivalence. Grounding: `ScfProblem::run/run_outcome/run_impl`, `step_kpoints`, `IterationRecord`, `FixedPointCertificate`; existing `crates/dft-engine/tests/{scf,kpoints,threading,pbe}.rs` and `docs/performance/baseline*.json`.

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

P11J overhead measurement status: local matched measurements are recorded below; quiet-host acceptance remains pending. P11R's historical build/test/fmt/clippy gates did not establish any overhead target.

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

## P11J local matched measurements (2026-11-04)

Status: local collection and equivalence pass; full acceptance **BLOCKED**, pending quiet-dev08 remeasurement and the unmeasured metrics below. No numerical code or default changed. Evidence: [benchmark.json](../../tests/reference/runtime/benchmark.json); reproduction and exact measurement boundaries: [runtime README](../../tests/reference/runtime/README.md).

Command: `python3 -I scripts/bench_runtime.py` (exit 0). Candidate `43868d9`, independent pre-extraction reference `9495d9fbaaf0e0bf6bf141ea61675e1d803eeeec` (parent of P11A), Rust 1.99.1, release, macOS 27.0.1 arm64, one requested/effective engine thread. CPU-brand/affinity inspection was unavailable in the sandbox. Backend is the repository default CPU/faer/rustfft stack. Same synthetic fixtures, compiler and machine for every pair. Three warm-ups and 20 measured pairs per case/mode; alternating order, no exclusions; bootstrap seed 20261104, 10,000 resamples. The host was **not quiet**: 1-minute load rose from 11.64 to 14.50. The first pilot overlapped workspace checks; this final collection followed their completion. These shared-host measurements remain provisional and are not accepted speedups.

### Targets and total wall time

External wall time includes process launch, full computation/drain, identity hashing, output and teardown. Internal total time includes setup and the completed driver/drain, excluding fingerprinting. The table reports both to expose launch/instrumentation noise on short workloads. A local per-cell statistical pass does not establish overall acceptance: all required cases, equivalence and quiet-host confirmation are required. Refer to each local statistical status below; no mode may be accepted without the quiet-host confirmation and the outstanding metrics.

| Fixture / budget | Mode (target) | Internal total median (ms) | External wall median (ms) | External paired overhead median % [95% CI]; local status | Internal paired overhead median % [95% CI]; local status |
| --- | --- | ---: | ---: | --- | --- |
| 0 / 1 | classic (<0.5%) | 1.699 | 3.574 | 1.166 [-0.389, 4.480]; inconclusive | 1.333 [0.605, 2.167]; not_met |
| 0 / 1 | observed (<1.0%) | 1.691 | 3.491 | 1.266 [-0.236, 4.513]; inconclusive | 1.878 [0.784, 3.047]; inconclusive |
| 0 / 1 | agent (<1.0%) | 1.701 | 3.535 | 3.718 [1.860, 6.038]; not_met | 2.168 [1.778, 3.454]; not_met |
| 0 / 200 | classic (<0.5%) | 2.102 | 4.017 | 2.844 [1.538, 5.500]; not_met | 2.491 [1.804, 2.848]; not_met |
| 0 / 200 | observed (<1.0%) | 2.109 | 4.048 | 2.563 [-0.251, 4.995]; inconclusive | 2.145 [1.269, 3.336]; not_met |
| 0 / 200 | agent (<1.0%) | 2.100 | 3.994 | 0.934 [-0.179, 4.391]; inconclusive | 2.506 [1.222, 2.959]; not_met |
| 1 / 1 | classic (<0.5%) | 2.214 | 4.137 | 0.714 [-1.077, 2.752]; inconclusive | 1.511 [0.927, 1.797]; not_met |
| 1 / 1 | observed (<1.0%) | 2.222 | 4.157 | 3.400 [2.776, 5.080]; not_met | 1.669 [0.575, 2.362]; inconclusive |
| 1 / 1 | agent (<1.0%) | 2.199 | 4.093 | 2.789 [-1.203, 4.053]; inconclusive | 1.356 [0.477, 2.209]; inconclusive |
| 1 / 200 | classic (<0.5%) | 3.846 | 5.837 | 1.723 [-0.178, 2.859]; inconclusive | 0.738 [-0.033, 1.684]; inconclusive |
| 1 / 200 | observed (<1.0%) | 3.845 | 5.806 | 0.691 [0.067, 1.629]; inconclusive | 0.620 [-0.149, 1.247]; inconclusive |
| 1 / 200 | agent (<1.0%) | 3.868 | 5.858 | 1.228 [0.732, 3.274]; inconclusive | 1.165 [0.494, 1.635]; inconclusive |
| 2 / 1 | classic (<0.5%) | 3.624 | 5.592 | 2.629 [0.359, 3.454]; inconclusive | 2.350 [0.848, 2.999]; not_met |
| 2 / 1 | observed (<1.0%) | 3.599 | 5.605 | 2.174 [0.605, 3.396]; inconclusive | 1.545 [0.547, 2.040]; inconclusive |
| 2 / 1 | agent (<1.0%) | 3.621 | 5.523 | 1.492 [-2.130, 2.000]; inconclusive | 1.177 [0.023, 2.147]; inconclusive |
| 2 / 200 | classic (<0.5%) | 4.746 | 6.719 | -1.501 [-2.337, 0.253]; met | -0.821 [-1.645, -0.174]; met |
| 2 / 200 | observed (<1.0%) | 4.742 | 6.700 | 0.833 [-1.640, 3.172]; inconclusive | -0.073 [-0.567, 1.318]; inconclusive |
| 2 / 200 | agent (<1.0%) | 4.714 | 6.584 | -0.221 [-1.766, 1.307]; inconclusive | 1.032 [0.411, 1.848]; inconclusive |

Fixture 0: cubic Gamma/nonlocal, three bands. Fixture 1: four k points/nonlocal, three bands. Fixture 2: skew/two species, local plus three-projector nonlocal, four bands. Budget 1 is deliberately unconverged; budget 200 converges. The stored JSON contains every reference/candidate sample, absolute RSS and paired deltas, CPU time, p5/p95 and confidence intervals. Historical Si timings are never used as the denominator of a synthetic run.

### Equivalence, SCF iteration cost and memory

All 360 measured pairs (720 processes) and all 240 separate probe runs match the independent reference. Every full density, eigenvalue and complex orbital bit is hashed; complete normalized result/history/certificate hashes also match. Only timing fields are normalized. Warm-ups are excluded from this count. Full runs finish in 11/18/13 ordinary iterations for fixtures 0/1/2. This is same-platform evidence, not a cross-platform promise or a QE comparison.

| Full fixture | Energy (Ha) | Certificates | Classic SCF ns/ordinary iteration | Observed SCF ns/ordinary iteration | Agent SCF ns/ordinary iteration | Paired RSS delta medians C/O/A (bytes) |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | -0.62379772639925168 | 1 | 43136.3 | 43784.1 | 43365.5 | [196608.0, 286720.0, 278528.0] |
| 1 | -0.51412521505579167 | 1 | 101243.1 | 101597.2 | 101905.1 | [163840.0, 237568.0, 245760.0] |
| 2 | -1.5249141148565482 | 1 | 101081.7 | 101444.0 | 101152.2 | [172032.0, 204800.0, 253952.0] |

SCF time per ordinary iteration includes initialization, result extraction, typed drain and nested certificates. It is **not** isolated boundary overhead. Whole-step samples are available for instrumented probes; separate certificate/noncertificate boundary timing and observer allocation counts/bytes remain unmeasured. RSS includes setup, allocator behavior and fingerprint buffers; negative paired differences are noise, not saved observer allocations.

### Event, checkpoint, control and saturation probes

20 repetitions per full fixture, separate from normal-mode timing. One synchronous typed consumer, no agent/LLM wait. Event probe serializes with the production CLI schema; normal timing drains typed events without serialization. Checkpoint exact continuation captures after iteration one, round-trips the production atomic file writer, restores and completes with identical output.

| Fixture | Typed construction ns/event | JSON encoding ns/event | Delivery p50/p95 ns | Exact capture median ns | Exact payload median bytes | Sync persistence median ms | Pause/resume receipt-to-applied median ns | Stalled queue delivered/dropped/high-water |
| --- | ---: | ---: | --- | ---: | ---: | ---: | --- | --- |
| 0 | 16.21 | 1747.3 | 0/76000 | 2250 | 7407 | 7.525 | 2833/875 | 8/10/8 |
| 1 | 16.50 | 1498.1 | 0/6000 | 3125 | 9726 | 7.394 | 2334/958 | 6/17/6 |
| 2 | 16.41 | 1825.5 | 0/161000 | 2562 | 13952 | 6.689 | 2334/917 | 6/12/6 |

Construction is a 10,000-event aggregate of typed construction/Arc clone/SystemTime, not internal enqueue cost. Timer-pair medians are approximately 28–30 ns; a 16–17 ns construction average is batch throughput, not resolution of an individually timed event. Delivery uses event wall timestamps and may include previous event serialization. Critical events remain retained in the tested routine-saturation runs; missing sequence intervals are stored. There is no slow-subscriber reconnect journal in this measurement.

Persistence is synchronous encoding plus production atomic write/fsync on local temporary storage, not asynchronous service latency. Payload size varies slightly with timestamp encoding. Restore/setup times and encoding times are stored. Pause/resume is already approved; zero approval wait and zero crossed boundaries are intentional, not latency claims for external agents.

### Outstanding acceptance work

- Claude must repeat on quiet dev08 through dft5090. The Mac run does not establish the 0.5%/1%/1% targets.
- Allocation traces (including disabled-event Classic), isolated generation/enqueue and boundary cost, separately timed validation/application, stale/rejected request latency, warm/export and incremental checkpoint memory, async persistence, producer/consumer rates, critical-storage exhaustion and reconnect/replay remain unmeasured.
- Expand slab/elongated, larger band counts and gap-range fixtures; repeat same-thread runs at 2/4 threads. No default or optimization is accepted on this limited matrix.
- Dashboard upload/registry, new README-generator result support and remote measurements are outside this ticket ownership. README regeneration and its existing tests pass with no generated diff.

### Historical Si baseline regression and independent extraction check

`bench_scf.py --cases gamma,k222 --repeats 3 --threads 1 --baseline
docs/performance/baseline-after-o4c.json --require-identical` exits **5**:
energy differences are 0 Ha (Gamma) and 3.552713678800501e-15 Ha (k222), both
within 1e-12 Ha, but both historical density-file hashes differ. This failure is
retained in `tests/reference/runtime/si-baseline-comparison.{json,log}`; the gate
was not relaxed and the baseline was not rewritten.

An independently archived/built pre-P11A CLI (`9495d9f`) run on the same Si inputs
matches the current CLI (`43868d9`) **exactly** in energy and density-file hash for
both cases (exit 0). Gamma is -7.295530468210769 Ha; k222 is
-7.831230057121137 Ha; both use 10 ordinary iterations. Thus the density difference
against the older O4c baseline predates runtime extraction, rather than being
introduced by this benchmark or the extraction. No claim is made about which
intervening change first produced it. Evidence: `si-pre-refactor.json`,
`si-provenance.json` (actual UPF/executable hashes and command), and
`si-pre-refactor.log`; reproduction helper: `compare_si_pre_refactor.py`, all in
`tests/reference/runtime/`. These are engine regressions, not newly produced QE
references. Historical/shared-host speedup columns are not acceptance evidence.

Validation: `scripts/check.sh` exits 0 (build/test/fmt/clippy); release checkpoint,
runtime-policy and stepwise suites total **31 passed, 0 failed**, including ignored
release tests; three Python statistics regressions pass. Commands and exits are
recorded in `tests/reference/runtime/p11j-checks.json`. No baseline or numerical
threshold was changed. The standalone Rust driver also passes explicit rustfmt.

The final measurement uses a corrected normal driver without probe-statistics
initialization inside the timed region. The complete first pilot and intermediate
invalid-provenance collection remain archived; see the runtime README's pilot
audit trail. No individual sample was removed. Final source hashes match the
current working-tree sources.
