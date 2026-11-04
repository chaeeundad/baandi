# P11J matched runtime benchmark

Run from the repository root (offline vendor directory required):

```sh
python3 -I scripts/bench_runtime.py
python3 -I tests/reference/runtime/test_bench_runtime.py
scripts/check.sh
```

The default is three warm-ups and 100 measured pairs (20 pairs gave 95 percent intervals up to 3.9 points wide and could not resolve the 0.5 and 1.0 percent targets, dev08 2026-11-04; the sufficiency rule below still requires at least 20) for each of six cases and
three modes, plus 20 repetitions of each separate probe on the three full cases.
`--threads 2` and `--threads 4` select additional same-thread comparisons;
`--output PATH` preserves each independent report. The default output is
`tests/reference/runtime/benchmark.json`. A nonzero exit is a build/run failure
or identity mismatch (5); exit zero means collection succeeded, **not** that an
overhead target passed. Incomplete reports are explicitly marked running/failed.

The reference is an independent `git archive` of
`9495d9fbaaf0e0bf6bf141ea61675e1d803eeeec` (parent of P11A), built offline with the
same pinned compiler and release settings. No checkout, worktree, git metadata,
remote host or numerical kernel is modified. Both drivers use the same fixture
source and synthetic generator. The candidate and reference binaries are hashed;
source commit, lockfile, generator and driver hashes are recorded. Temporary
archives/binaries are removed; compilation artifacts use the worktree's target.
No new dependency is installed. The standalone Rust driver is compiled against
Cargo's reported rlib artifacts, not guessed library names.

Fixture kinds reuse the stepwise equivalence matrix: (0) cubic Gamma, three
bands and three nonlocal projectors; (1) four k points with the same potential;
(2) skew cell, two species (nonlocal plus local), four bands. Each has a retained
one-iteration exhausted case and a max-200 converged case. Defaults are unchanged.
This matrix does not establish slab/elongated, large-band or gap-range coverage.

Each mode is paired with a fresh reference process; order alternates. No sample
is excluded. The sample budget is fixed before running. Medians, interpolated
p5/p95, paired percentage overhead and 10,000 seeded bootstrap resamples of the
median overhead are recorded. A target is met only if the upper 95% bound is
strictly below it, there are at least 20 pairs/three warm-ups and identity passes.
An interval straddling the threshold is inconclusive. These local statistics
are provisional; acceptance requires Claude's quiet-dev08/dft5090 remeasurement.

## Measurement boundaries

- `total_ns`: monotonic fixture setup through completed driver and final typed
  event drain. `setup_ns` and `scf_ns` split it. Fingerprinting/printing are outside
  this interval. Classic calls `run_outcome`; Observed and Agent use Runtime with
  a single synchronous typed scalar consumer. There is no deferred sink work.
- `process_wall_ns`: external monotonic launch-through-exit, including fingerprint,
  JSON output and teardown. `process_cpu_s` and `peak_rss_bytes` come from per-child
  `wait4`, with macOS/Linux RSS units normalized. RSS includes the entire process,
  not just observer allocations. Absolute samples and paired RSS deltas remain
  available; a negative delta is allocator/process noise, not negative allocation.
- Normal runs do not serialize events. Separate `events` probes measure real CLI
  event JSON encoding, bytes, count and event timestamp-to-consumer latency.
  Delivery timestamps use SystemTime, not a monotonic producer stamp; clock
  adjustments and timer granularity limit interpretation. Synthetic typed-event
  construction includes Arc cloning and SystemTime but **not private queue push**.
  Timer-pair cost is reported; sub-resolution costs are unresolved.
- `steps_ns` in probes is whole Runtime::step latency, not isolated boundary
  overhead. Dividing `scf_ns` by ordinary iterations includes certificate work;
  certificate counts are reported. No isolated boundary or allocation claim is made.
- `saturation` stalls the consumer for the entire run with routine capacity one.
  It records delivered/drop counts and missing sequence intervals; critical
  capacity stays 128. No reconnect/replay journal or critical exhaustion timing
  is inferred from routine saturation.
- `control` measures preapproved pause/resume at the initial boundary. Submit and
  receipt-through-application nanoseconds are recorded; approval wait and crossed
  boundaries are zero by construction. Validation and application are synchronous
  inside safe_point, so their separate internal costs are unavailable. This does
  not measure a network client, LLM or stale request latency.
- `checkpoint` captures exact runtime state after step one, encodes it, persists
  with the production synchronous atomic writer, reads/round-trips the encoding,
  restores a new Runtime and finishes. Capture, encoding, persistence, restore
  (including new model setup), and payload bytes are separate. Persistence includes
  encoding again. The fixture identity comes from existing synthetic checkpoint
  tests; the executable hash is real, but placeholder model hashes are **not** UPF
  provenance. Warm/export, async persistence and incremental peak memory are not
  measured by this probe.

Every repeated result must match its reference in convergence, iteration and
certificate counts, energy, full density/eigenvalue/complex-orbital bit hash and
normalized complete result/history/certificate hash. Only duration fields are
normalized. Probe trajectories, including checkpoint continuation and stalled
queues, must match too. This is same-build-platform numerical evidence, not QE
physics validation or a promise across platforms.

Historical `docs/performance/baseline*.json` files and hashes are embedded for
traceability. Their Si timings cannot be used as denominators for these synthetic
fixtures. Use `scripts/bench_scf.py` separately for the historical Si regression;
never label this synthetic pre-extraction engine as QE.

## Handoff

Claude runs the same command on quiet dev08 via dft5090, records CPU/backend,
affinity/load and storage, and repeats at required thread counts. No remote action
is taken here. Full acceptance additionally needs observer allocation traces,
private enqueue/boundary instrumentation outside numerical kernels, broader
fixture coverage, warm/export costs and transport/journal fault measurements.
README result-generator and dashboard updates belong to orchestrator-owned files;
this ticket does not extend their formats or publish datasets.

## Historical Si comparison

After building the **candidate** CLI (`cargo build --offline --release -p dft-cli`),
run the following. The historical comparison can exit 5 for a density hash change;
retain that failure and the evidence rather than relaxing the gate.

```sh
python3 -I scripts/bench_scf.py --binary target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222 --repeats 3 --threads 1 \
  --baseline docs/performance/baseline-after-o4c.json --require-identical \
  --record tests/reference/runtime/si-baseline-comparison.json
python3 -I tests/reference/runtime/compare_si_pre_refactor.py
```

The second command archives/builds the pre-extraction CLI and compares it with the
candidate evidence using exactly the same Si input generator. It records actual
UPF and executable hashes in `si-provenance.json`; it never copies the UPF.
Both build helpers leave `target/release/dft-cli` containing the last built
**reference** CLI. Rebuild the candidate before any subsequent candidate CLI run.
These Si timings are three-run scientific regression checks, not 20-pair overhead
acceptance measurements. The UPF path is local to this Mac; other hosts must use
the same read-only file hash and adapt the path before running this local helper.

## Pilot audit trail

`benchmark-instrumented-pilot.json` retains the first full collection, whose
normal driver unnecessarily created a probe statistics object inside the timer.
It is superseded because this fixed adapter cost biases short SCF timing. No
sample was selectively removed. `benchmark-transition-run.json` retains an
intermediate full collection marked **invalid_provenance**: the driver was edited
while its build was in progress, so the recorded initial source hash does not
identify the compiled source. Do not use it for acceptance. `benchmark.json` is
the fresh, fixed-driver collection with stable source hashes. Workload inputs,
physics, run counts, bootstrap seed and thresholds are unchanged.
