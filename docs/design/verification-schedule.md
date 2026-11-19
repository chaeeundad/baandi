# Verification scheduling experiment (F-VERIFY-SCHEDULE)

Status: default ON since 4c3d09d (user-relayed dev08 acceptance).
G-SCHED-OVERHEAD removes redundant diagnostic storage; its performance effect
awaits the orchestrator's dev08 remeasurement. The original Mac NO-GO below is historical
and is superseded as a universal conclusion by the relayed dev08 Si64 result.
Certificate and stop rule unchanged. Scope: closed-shell `ScfSession`; MPI
measurements below are relayed, not reproduced locally. Spectral-audit sessions
decline scheduling.

## G-SCHED-OVERHEAD: allocation audit and lazy diagnostic storage

Source baseline: `b58a5a35523419a2ae829fd39c36a1a326307e23`; Mac pure-Rust
release, one thread. The user-relayed +0.29% Si32 regression is a **whole-run**
measurement. Dividing its 23 ms by ten maps does not identify where it occurs.
Source inspection and a counting allocator show **no schedule operations or
allocations on ordinary non-eligible maps even before this patch**. There is no
per-map history/spectrum clone, string formatting, JSON, timer or environment
lookup attributable to the switch. The pre-existing clocks, history construction
and numerical work are shared by off/on.

| Boundary | Before | After |
| --- | --- | --- |
| Ordinary, not eligible for verification | No scheduler branch/call/allocation/copy | Unchanged, zero |
| First eligible attempt, before verification | Mode comparison, cold helper call, history-last read, schedule initialization, `Vec::with_capacity(max_iterations)`, zero prediction, scalar decision row push | Failed-attempt check is false; no scheduler call or allocation |
| First successful attempt | Fill scalar error/ratio/pass fields | Use the existing actual attempt; no retry journal |
| First failure | Fill the already allocated row | Lazily allocate retry journal and record that same zero-forecast attempt |
| Eligible retries | Existing scalar predictor, counter and decision rows | Same policy, counters and rows |
| Final result | Clone the journal's live rows (second allocation); drop the original with the session | Move an existing journal; if the first attempt passed, create exactly one diagnostic row at output |
| CLI serialization | Serialize required optional schedule metadata | Unchanged; outside this ticket's ownership |

The first-failure journal still reserves the maximum iteration budget, so retry
pushes have the same capacity behavior. Forecast arithmetic and the complete
`verify_fixed_point` / inner certificate region are byte-identical to HEAD.
Their SHA-256 values are respectively
`82c6e5f53f368c4fc7cdf008db59f8d80a0ebeb89f2147573251e797d26a16a2` and
`749fdf5f373f7584300abea7a019f9100b2ca18f324803c09cb0c1ee073b250b`.
Attempts remain in the session as before; only diagnostic storage is moved to
its terminal result. A restored session emits no synthetic row for an attempt
before the restored boundary. No new public type, dependency or switch is added.

### Allocation evidence

`tests/verify_schedule_alloc.rs` installs a test-only forwarding System allocator
and counts alloc/alloc_zeroed/realloc around each real `ScfSession::step`. One
test function and a one-thread skew-cell, nonlocal synthetic fixture avoid other
test threads. Setup, warmup, assertions and output are outside counting. A known
allocation checks that the probe is live. Off/on/trace scientific state is fully
compared (energy, density, orbitals, history, attempts and work; clocks and the
diagnostic journal excluded).

Before: off allocation counts by map `[440, 59, 60, 45, 58, 57, 96]`;
on/trace `[440, 59, 60, 45, 58, 57, 98]`. After: on/trace
`[440, 59, 60, 45, 58, 57, 97]` in release. Thus **zero added allocations for
each nonterminal ordinary map**, before and after; **two -> one** added
allocation at the terminal/output boundary. The surviving allocation is the
required one-row public `Vec` diagnostic, with capacity one instead of 100.
This is deliberately not a claim of zero total allocations in the SCF or zero
allocation for serialized diagnostics. The regression test asserts this exact
boundary distinction rather than hiding the output allocation.

The no-failure unit test also calls the production helper with empty history
and no attempts, proving it returns before touching history or creating a
journal. The forced-failure trace test checks every real-attempt row, first
zero forecast, measured ratios and prediction count. The synthetic wait test
clears the fixture's default-on journal before replacing its history; otherwise
the real passing row incorrectly precedes the artificial wait row (an existing
test assumption invalidated by default-on).

### Failed-attempt identity and verification gates

A separate old/new release comparison ran Si16 at 1e-10 Ha and Si64 at 5e-9 Ha,
seed 601 / fractional amplitude 1e-3, in off/on/trace: **six same-mode
comparisons (12 SCFs), all bit-identical** in energy, density bytes, eigenvalues,
occupations, energy components, scientific history, attempts, work and the
complete schedule journal (only timing fields excluded).

| Case | Off / trace attempts, before = after | On attempts, before = after | On journal |
| --- | --- | --- | --- |
| Si16, 1e-10 Ha | [10,11,12] | [10,12] | First attempt fails; wait at 11 with predicted n=1; pass at 12 |
| Si64, 5e-9 Ha | [11,12,13] | [11,13] | First attempt fails; wait at 12 with predicted n=1; pass at 13 |

Both on/trace journals have two predictor evaluations. Trace always attempts,
including the counterfactual wait, and its measured errors/ratios/pass values
match the baseline exactly. On Si64, the measured first failure is
5.310026608990284e-8 Ha; the final pass is 1.1717501795885816e-9 Ha, ratio
0.022066747793781642. These are local internal measurements, not the relayed
dev08 energy values. Same-mode before/after identity does not imply on/off
identity when scheduling skips a failed solve.

- Requested `CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli`:
  exit **1**, build/fmt/clippy PASS; 51 binaries, **466 passed, 3 failed,
  134 ignored**, 184 s wall. The three failures are the known CLI pseudo-cache
  rejection of an in-repository TMPDIR. The sandbox also denies `ps` load
  collection; the build/test/fmt/clippy gates still execute.
- Final release `--lib verification_schedule`: **4 passed**, exit 0.
- Final release `--test verify_schedule_alloc`: **1 passed**, exit 0.
- `python3 -I target/scratch/analyze.py` and `identity.py`: exit 0. The analysis
  compares every timing run's full scientific JSON (excluding clocks, timing
  derivatives and journal) and density SHA against the baseline.
- Optional `python3 -I scripts/tests/test_readme_figures.py`: exit **1**, 8 of
  9 pass. The unchanged structural-scope test expects the obsolete phrase
  `standard profile is unsupported`; generated README and figure freshness
  tests pass. Reporting sources were not edited outside ownership.
- An initial schedule-test run exposed the inherited-journal assumption
  described above; after fixing the fixture, all four release schedule tests
  and the engine fast suite pass. No tolerance was changed.

Summary SHA-256:
`d91cef1bdc53fc43205605fc0bb6f4f46a422cfc00cfd7bc05669153b0d246a5`.
Identity summary SHA-256:
`8ad80c4b33b825993f131d1cba60548ca0a789942622df380aaaea11632a102c`.

### Measurement protocol and remaining acceptance boundary

Use the checked-in `fast_screen.py` input stream: seed 601, fractional
amplitude 1e-3, Gamma, one thread, 30 Ry, LDA-PZ, production ladder/defaults,
standard threshold 5e-9 Ha. Both cases finish ten fine maps. Si32 is the
no-failure control (one actual attempt at map 10, zero predictions); Si16 has
attempts at maps 9 and 10, one failure and one predictor evaluation, without a
wait. Do not call the Si16 run a no-failure control.

For each metric X and pair i, compute `d_i = 1000 * (X_on - X_off) / 10`
in ms per fine map. Report `mean(d_i)`, sample `stdev(d_i)/sqrt(N)`, and
`abs(mean)/SE`. No outliers are removed. Eight AB/BA alternating pairs were
collected per case/build. The last baseline Si32 and first baseline Si16 pairs
are conservatively excluded because their launch windows might overlap: the
baseline tables use seven pairs each. Candidate cases run sequentially, with
no local build/test running during either timing series, using eight pairs.
The Mac is shared; load and raw timings are retained. These are screening
observations, not small-gain adoption evidence. Separate before/after windows
must not be subtracted to claim a speedup.

`ordinary` means `timings.steps_s - timings.certificate.wall_s`; it includes
fine ordinary numerical work and session overhead. `bookkeeping` is the existing
`timings.closure.step_overhead_s` bucket, not a scheduler-exclusive timer.
`wall` is process wall time including setup, coarse ladder, verification and
output, divided by the number of fine maps only for normalization. It is not
an actual map timer. No additional clocks are placed in production code.

### Mac paired observations (milliseconds per fine map, on minus off)

| Build / case | Pairs | Ordinary mean ± SE | abs(mean)/SE | Bookkeeping mean ± SE | Whole process mean ± SE | 1-min load range |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| before / si16 | 7 | +0.811029 ± 3.256482 | 0.25 | +0.003192 ± 0.021432 | +0.472154 ± 4.888665 | 4.46–4.95 |
| before / si32 | 7 | -0.566346 ± 3.405864 | 0.17 | +0.018304 ± 0.038778 | -2.626031 ± 9.529771 | 4.25–6.25 |
| after / si16 | 8 | -1.092127 ± 2.031926 | 0.54 | +0.009663 ± 0.007714 | -1.242311 ± 2.475990 | 6.43–8.12 |
| after / si32 | 8 | -5.159820 ± 3.574622 | 1.44 | +0.012437 ± 0.020274 | -12.563878 ± 7.306373 | 8.12–20.04 |

No delta is resolved by this sample. In particular, these numbers do **not**
prove that the dev08 +0.29% disappeared, or support a Mac speedup claim. There
was no dedicated process warmup; first runs remain in the sample except for the
two potentially overlapping pairs explicitly excluded above. Candidate Si32
starts after the gate build, with elevated decaying load. Source/counter evidence
is stronger than these timings for locating work.

All 64 timing runs converge and verify. Every before/after off/on result for a
given case has identical energy, density bytes, eigenvalues, occupations,
energy components, full scientific history, attempts and work counters
(timing fields and the optional schedule journal excluded). Si32 energy is
`-126.70343204998852 Ha`; Si16 is `-63.129095464190456 Ha`. These are internal
same-input identities, not QE or absolute-accuracy references.

Reproduce each timing series (use distinct output directories):

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export TMPDIR_990="$PWD/target/scratch" TMPDIR="$PWD/target/scratch"
cargo build --offline --release -p dft-cli
python3 -I scripts/fast_screen.py --binary target/release/dft-cli \
  --upf ../pseudo/Si.pz-vbc.UPF --case benchmarks/cases/si32-gamma.toml \
  --pairs 8 --env-a BAANDI_VERIFY_SCHEDULE=off \
  --env-b BAANDI_VERIFY_SCHEDULE=on \
  --out target/scratch/after-si32 --json target/scratch/after-si32.json
# Repeat sequentially for si16-gamma.toml and after-si16.
cargo test --offline --release -p dft-engine --test verify_schedule_alloc -- --nocapture
cargo test --offline --release -p dft-engine --lib verification_schedule -- --nocapture
CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli
```

Baseline binary SHA-256:
`919a2ebfc36a6ae1b3aee030342cec0e18b1fbd8359e765a78a0940a52420594`.
Candidate binary SHA-256:
`0e3a4bbfdefdb522aebc20589ccdd6365ce053392768845035a7ac0cb71213b5`.
External Si UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Platform: macOS 27.0.1 arm64, no extra release flags. Raw evidence and the
analysis/identity scripts are in `target/scratch/{before,after}-si{16,32}*`,
`alloc-{before,after}.log`, `identity-*`, `manifest.json`, `source-identity.json`,
`summary.json` and `check.log`; preserve this scratch before worktree cleanup.
No UPF is committed. README, dashboard/session JSON and optimization-log files
are outside the ticket's owned files; the master owns their publication.

**dev08 handoff:** repeat one-core Si32 eight alternating off/on pairs and
one-core Si64 crossed A/B within a single clean window. Retain the existing
input/PP/build/binding/threshold/start conditions, full phase timings, binary
hashes and scientific fingerprints. Compare candidate on with baseline on,
as well as candidate off/on; verify Si64 on attempts `[11,13]` and the exact
trace records. No remote run is claimed here. The source/allocator audit does
not establish the cause of the relayed 23 ms, nor that this patch removes it.
If it persists, profile the verification solve, allocator/cache effects and
CLI output on dev08 rather than attributing it to nonexistent per-map records.

## Follow-up: lazy prediction and durable trace

The earlier version printed trace only to **stderr**, not stdout or results.json.
The follow-up removes that output and writes an optional `verification.schedule`
section in results.json in both trace and on modes; off has no section. The
follow-up request authorizes the necessary additive report types in
`convergence.rs`, CLI serialization in `run.rs`, and CLI output tests in addition
to the original ticket's files.

Scheduling now lives **inside the existing eligible-verification branch**. There
is no schedule condition, call, prediction, allocation, logging or counter update
on ordinary maps before the first failed attempt. At the first real attempt the
journal is allocated lazily (capacity = current maximum iterations); prediction
is bypassed. After a failure, only eligible retry boundaries call the pure,
allocation-free scalar predictor. The helper is cold and not inlined into the
ordinary-map path. The journal is preallocated for a fixed iteration budget;
explicit runtime budget increases may require journal growth at later attempt
boundaries. No numerical kernel knows JSON or writes files.

The counter `prediction_evaluations` counts actual predictor calls, including
fallback forecasts. It is **zero** for a no-failure run. Output allocation and
JSON serialization at the terminal attempt are intentionally additional work;
zero **per-map/predictor** work is not a claim of zero total reporting cost or a
measured zero wall-time delta. Scientific bits and counters must match, while
the requested schedule metadata necessarily differs from off.

Each `verification.schedule.records` row contains:

| Key | Meaning |
| --- | --- |
| `map_index` | One-based ordinary map index at the eligible boundary |
| `predicted_additional_maps` | Forecast n, before the last-budget-map override |
| `action` | Actual `attempt` or `wait`; trace always attempts |
| `would_wait` | Counterfactual on decision including budget/scope guards |
| `coulomb_error_ha` | Measured real-attempt error, null on a wait |
| `ratio_to_previous_attempt` | Current/previous real-attempt Coulomb error; null on waits, first attempt, zero denominator or nonfinite ratio |
| `pass` | Real attempt's verdict; null on a wait |

The section also contains `mode` and `prediction_evaluations`. No schedule
record declares convergence: `VerificationReport::passed` reads only actual
attempts. The caller forces `passed = false` on a wait. Measurement fields are
filled only **after** `verify_fixed_point` returns. A final-budget boundary
attempts even when the forecast would otherwise defer.

The journal is diagnostic, not checkpoint state. Exact continuation restores the
existing attempts/history, so the numerical schedule remains reconstructible;
the emitted journal covers the resumed session's boundaries only. It does not
invent predictions for attempts made before that checkpoint. Spin/FD drivers do
not use this closed-shell scheduler and do not emit this section.

### Relayed dev08 evidence (user, commit 41d2d79; before this follow-up)

These are supplied measurements, not local reruns. Each ratio compares on/off
within the stated parallel model; no cross-model speed comparison is made.

| Model / case | On vs off wall | Crossed-pair evidence | Work and state |
| --- | ---: | --- | --- |
| One core / Si64 | -7.06% | 4/4 pairs negative; abs(mean)/SE = 218 | Attempts [11,12,13] -> [11,13]; failures 2 -> 1; Davidson 71 -> 69; H 4008 -> 3871; abs(dE) = 4.11e-11 Ha |
| One core / Si32 | +0.28% | 4/4 pairs positive; abs(mean)/SE = 4.87 | No failed attempts; identical bits and counts; predictor-machinery overhead identified by the master |
| MPI np=8 / Si64 | No effect | Schedule unchanged | No scheduling gain reported |

**Rounding/path class, not path-identical:** skipping an attempt that would have
failed removes that solve and its state/reuse perturbation from the continuation.
The delivered state may therefore differ at the 1e-11 Ha level. The reported
4.11e-11 Ha difference is inside the strict 1e-10 Ha class and was accepted by the
user; it must not be called path identity. This on/off agreement is not by itself
an independent absolute-error oracle. Certificate, criterion, class thresholds,
and default OFF remain unchanged. Trace preserves the original numerical path.


### Follow-up Mac validation and counters

Source baseline: local HEAD `6ca9dad` (the prior scheduler); uncommitted follow-up.
Release, pure Rust, one thread; same eight systems, seeds 601–603, fractional
amplitude 1e-3 and thresholds 5e-9 / 1e-10 Ha as the historical study. **48 triples
(144 SCFs) completed, all converged and verified.** For every geometry, each
new arm is bit-identical to its corresponding previous arm in energy, density
bytes, eigenvalues, occupations, energy components, scientific history, attempts
and work counts (timings and the new diagnostic section excluded). Off vs trace
also matches in all 48 pairs. Thus the scheduling policy and its path-class
behavior have not changed in this overhead/reporting fix.

**19 no-failure geometries, 120 fine ordinary maps:** on and trace each perform
**zero predictor evaluations**, emit exactly one attempt record at the terminal
map, and match off scientific bits and counts. All modes have the same total
14,017 H columns over those geometries (coarse + coarse start + fine + attempts).
The following small-case rows cover all three seeds, not selected examples.

| Case | Seeds | Threshold (Ha) | Fine maps off/on per seed | Total H off/on per seed | Attempts off/on | Predictor evaluations on/trace |
| --- | --- | ---: | --- | --- | --- | --- |
| syn-p3 | 601,602,603 | 5e-9 | 5/5, 5/5, 5/5 | 13/13, 13/13, 13/13 | 1/1 each | 0/0 each |
| syn-p3 | 601,602,603 | 1e-10 | 6/6, 6/6, 6/6 | 14/14, 14/14, 14/14 | 1/1 each | 0/0 each |
| syn-spd | 601,602,603 | 5e-9 | 3/3, 3/3, 3/3 | 11/11, 11/11, 11/11 | 1/1 each | 0/0 each |
| syn-spd | 601,602,603 | 1e-10 | 5/5, 5/5, 5/5 | 13/13, 13/13, 13/13 | 1/1 each | 0/0 each |
| si32 | 601,602 | 5e-9 | 10/10, 10/10 | 3064/3064, 3077/3077 | 1/1 each | 0/0 each |

Across the whole on study, 37 eligible post-failure boundaries call the predictor;
81 actual attempts and **4 wait rows** are reported. All 4 waits have null measured
error, ratio and pass. Every real-attempt row matches the attempt array exactly,
including its ratio to the previous real attempt. Trace reports all 85 actual
attempts, with `action=attempt` even on its counterfactual waits. Off emits no
`schedule` key. Original aggregate on/off H counts remain 48,676 / 48,944;
there is no additional H or Davidson work from this follow-up.

Example of an actual on journal (si16 seed 601, 1e-10 Ha):

```json
{
  "mode": "on",
  "prediction_evaluations": 2,
  "records": [
    {
      "action": "attempt",
      "coulomb_error_ha": 2.522448026517058e-09,
      "map_index": 10,
      "pass": false,
      "predicted_additional_maps": 0,
      "ratio_to_previous_attempt": null,
      "would_wait": false
    },
    {
      "action": "wait",
      "coulomb_error_ha": null,
      "map_index": 11,
      "pass": null,
      "predicted_additional_maps": 1,
      "ratio_to_previous_attempt": null,
      "would_wait": true
    },
    {
      "action": "attempt",
      "coulomb_error_ha": 1.1032070222183943e-11,
      "map_index": 12,
      "pass": true,
      "predicted_additional_maps": 0,
      "ratio_to_previous_attempt": 0.004373556999474351,
      "would_wait": false
    }
  ]
}
```

Validation:

- `CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli`: exit 1;
  build/fmt/clippy PASS, 50 test binaries: **463 passed, 3 failed, 134 ignored**.
  Only the three known CLI pseudo-cache tests reject the required in-repository
  TMPDIR. Sandbox `ps` load collection is denied; the gates still execute.
- New engine tests prove no-failure on/trace identity and zero predictions,
  wait rows cannot pass, even a corrupted diagnostic pass cannot certify a
  state, trace attempts despite a positive forecast, and budget exhaustion
  overrides a wait. Journal capacity remains unchanged on the tested retry.
- New CLI test reads results.json for all three modes and checks the optional
  section, measured fields, zero-prediction counter, density bytes and counts.
- The complete `verify_fixed_point`/inner region remains byte-identical to HEAD:
  SHA-256 `4d6528a178dcd5c67f232cf92af9eaf861d60bd9ac1ec9072a7360b6cda53acc`.
- Mac timings are not an overhead acceptance measurement: this counter study
  overlapped builds and uses one count sample per arm/geometry. In particular,
  it does not establish that dev08's former +0.28% wall delta is now zero.

Reproduce from the scratch harness printed at the end of this document:

```sh
export TMPDIR_990="$PWD/target/scratch" TMPDIR="$PWD/target/scratch"
cargo build --release -p dft-cli
SCHEDULE_BINARY="$PWD/target/release/dft-cli" \
  python3 -I target/scratch/study.py off,on,trace followup-ab-fresh
# Each */{on,trace}/results.json has verification.schedule; off omits it.
```

Raw current evidence: `target/scratch/followup-ab`, `followup-summary.json`,
`followup-analyze.py`, `followup-check.log`. Analysis verifies every recorded
attempt/ratio, null wait values, all numerical identities, and no-failure counters.
Raw files are uncommitted scratch; preserve them before cleanup.

Follow-up binary SHA-256: `d0f77cc8895a39e57028c4ef5ba9afb0930432d24b4aa7dcff32ffd95e7b9f5b`.

Follow-up summary SHA-256: `967d9d48899e757e89a60bbf878da6d548432790bcce73fe8aebd301589bbebc`.

## Predictor and cost model (revised after initial trace)

An ordinary map costs m; an attempt is approximately 1.3 m (ticket's dev08
observation). Waiting one map costs m; attempting and failing costs 2.3 m
including the forced ordinary map. The proposed decision is to wait when the
estimated failure probability exceeds 1 / 2.3 = 0.4347826. This is a screening
cost model, not a measured probability or a wall-time guarantee. Waiting can
lose when an attempt would have passed, and failed-attempt reuse changes the
subsequent trajectory; actual A/B runs are required.

The first attempt is unconditional. There is no system-independent conversion
from the ordinary stop estimator to the verification error (see ledger
VERIFY-FORESIGHT-ANALYSIS). After failure at map j with Coulomb error V_j, the
prediction at map i is V_hat_i = V_j S_i / S_j. S is the existing Coulomb-metric
density residual energy, not the L2 `density_residual`; no FFT is added. This
re-anchors the forecast to the ordinary maps actually performed, including the
effect of failed-attempt reuse. For future maps, q = S_i / S_(i-1); with two
failed attempts, combine q with their geometric per-map error ratio using a
geometric mean. Predict V_hat_i q^n. The failure vote is the point forecast
indicator V_hat_i >= threshold (0 or 1); it must be empirically calibrated and
is not a probability bound. Choose the first n in 0..3 for which the vote falls
below 1/2.3, capped by three maps since the last actual attempt. Never defer the
last budgeted map. Missing, nonfinite, nonpositive or noncontracting evidence,
or a mixer reset/recovery since the anchor, means attempt now.

The initial trace-only candidate used three pre-failure S ratios and a vote
across their extrapolations. It was rejected before enabling on: the first 14
runs had seven proposed waits, five at an actual passing attempt. Pre-failure
history does not represent the continuation after reuse. The revised candidate
above uses an observed post-failure change, with no fitted multiplier. This is
an explicit revision, not a successful pre-registered test of the initial rule.

Trace records this prediction at eligible attempt boundaries and performs the
original attempt; the CLI writes the journal to `verification.schedule`. All predictor inputs come from checkpointed history and
attempts; no countdown or hidden restart state. Off avoids prediction entirely.
On was implemented only after the revised trace calibration below, then tested in real A/B runs. Passing always requires the
unchanged real verification. Missing history, nonfinite data, reset histories,
and noncontracting evidence fall back to attempting.

## Evidence protocol

Mac: si8, si16, si32, si-skew, MgO8 (PBE), synthetic projector fixtures;
seeds 601–603, fractional displacement 1e-3, thresholds 5e-9 and 1e-10 Ha.
Use production defaults, including cutoff ladder and verification reuse.
Record all ordinary, coarse and verification H columns, attempt/failure counts,
maps after first attempt, predictions vs actual pass distance, and regressions.
Trace must match off in energy bits, density, work and scientific history.
If predictions do not calibrate, keep trace only. If fewer than four systems
benefit, close as negative; no default adoption from Mac times.

Scratch and raw evidence live under `target/scratch`; UPFs remain external.
Dashboard/README/reference artifacts are outside the ticket's owned files.

## Original Mac result: NO-GO (historical)

Completed 2026-11-19, source `dba7d8242fe74494300ba9bc6e2474f8dc3b047d` plus this uncommitted patch; pure Rust release, one thread, macOS-27.0.1-arm64-arm-64bit.
Only **3 of 8 systems** show any H-column reduction (4 of 48 runs). Even the liberal aggregate-per-system test misses the required four systems; every system's median paired reduction is zero. **Original local decision: do not promote. Default remains off pending the follow-up dev08 decision.** On is retained only as the bounded, explicit experiment for review and reproduction. No dev08 timing was performed.

Final A/B: **48/48 converge and pass a real verification in both arms**; total H columns **48,944 -> 48,676 (-0.547565%)**, attempts **85 -> 81**, failed attempts **37 -> 33**. Fine maps after the first attempt total **37 -> 37**. No run adds an ordinary map or increases total H. All 24 runs at 5e-9 Ha are unchanged. At 1e-10 Ha, four runs improve. Mac wall times were deliberately not used: studies overlapped builds and other study processes; these are deterministic work counts.

Pristine HEAD binary vs final off: **48/48 bit-identical** energy, density bytes, eigenvalues, occupations, energy components, scientific history and attempts (wall times excluded). Final trace vs off: the same **48/48** identities. The initial trace-only candidate and revised trace each also ran 48 geometries. A separate strict-reference run was made for each of the 24 geometries. The 1e-10 rows below use the same standard path with the tighter class threshold, not a change to the strict-reference algorithm.

Maximum delivered energy difference from the same-geometry strict reference: **1.4630785472036223e-10 Ha** at the 5e-9 Ha class and **8.185452315956354e-12 Ha** at the 1e-10 Ha class. Maximum on/off difference **3.169020601490047e-12 Ha** (standard rows: exactly zero). These are internal references, not QE measurements or independently established material accuracy.

### Distribution over all seeds and both thresholds

H includes ordinary history, every real verification, coarse history and coarse atomic-start columns. No double counting of verification in history and the attempt array. A coarse column and fine column have different cost: this is work accounting, not a wall-time estimator. The latter distinction does not affect these deltas because the coarse work is identical.

| System | Total H off / on | Attempts off / on | Failed off / on | Changed runs / 6 | Paired H change min / median / max (%) | Tail maps off / on |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| mgo8 | 9285 / 9285 | 7 / 7 | 1 / 1 | 0 | 0.0000 / 0.0000 / 0.0000 | 1 / 1 |
| si-skew | 1655 / 1655 | 15 / 15 | 9 / 9 | 0 | 0.0000 / 0.0000 / 0.0000 | 9 / 9 |
| si16 | 10879 / 10682 | 14 / 12 | 8 / 6 | 2 | -5.5584 / 0.0000 / 0.0000 | 8 / 8 |
| si32 | 19882 / 19882 | 10 / 10 | 4 / 4 | 0 | 0.0000 / 0.0000 / 0.0000 | 4 / 4 |
| si8 | 4880 / 4812 | 13 / 12 | 7 / 6 | 1 | -7.3276 / 0.0000 / 0.0000 | 7 / 7 |
| syn-gaussian-fcc | 2210 / 2207 | 14 / 13 | 8 / 7 | 1 | -0.7194 / 0.0000 / 0.0000 | 8 / 8 |
| syn-p3 | 81 / 81 | 6 / 6 | 0 / 0 | 0 | 0.0000 / 0.0000 / 0.0000 | 0 / 0 |
| syn-spd | 72 / 72 | 6 / 6 | 0 / 0 | 0 | 0.0000 / 0.0000 / 0.0000 | 0 / 0 |

### Complete per-run table

All rows passed verification. S = 5e-9 Ha; T = 1e-10 Ha. H is the complete column count above. A/F = attempts/failed attempts. Tail = ordinary maps after the first attempt. Entries are off -> on.

| System | Seed | Class | H | A/F | Tail | On error against strict (Ha) |
| --- | ---: | --- | ---: | --- | ---: | ---: |
| mgo8 | 601 | T | 1577 -> 1577 | 1/0 -> 1/0 | 0 -> 0 | 8.185452e-12 |
| mgo8 | 601 | S | 1513 -> 1513 | 1/0 -> 1/0 | 0 -> 0 | 3.501555e-11 |
| mgo8 | 602 | T | 1527 -> 1527 | 1/0 -> 1/0 | 0 -> 0 | 2.216893e-12 |
| mgo8 | 602 | S | 1562 -> 1562 | 2/1 -> 2/1 | 1 -> 1 | 2.671641e-11 |
| mgo8 | 603 | T | 1575 -> 1575 | 1/0 -> 1/0 | 0 -> 0 | 2.273737e-13 |
| mgo8 | 603 | S | 1531 -> 1531 | 1/0 -> 1/0 | 0 -> 0 | 4.143885e-11 |
| si-skew | 601 | T | 273 -> 273 | 2/1 -> 2/1 | 1 -> 1 | 3.614886e-13 |
| si-skew | 601 | S | 278 -> 278 | 3/2 -> 3/2 | 2 -> 2 | 2.626344e-12 |
| si-skew | 602 | T | 274 -> 274 | 2/1 -> 2/1 | 1 -> 1 | 4.245493e-13 |
| si-skew | 602 | S | 280 -> 280 | 3/2 -> 3/2 | 2 -> 2 | 2.662759e-12 |
| si-skew | 603 | T | 274 -> 274 | 2/1 -> 2/1 | 1 -> 1 | 2.531308e-13 |
| si-skew | 603 | S | 276 -> 276 | 3/2 -> 3/2 | 2 -> 2 | 1.421974e-12 |
| si16 | 601 | T | 1959 -> 1870 | 3/2 -> 2/1 | 2 -> 2 | 6.110668e-13 |
| si16 | 601 | S | 1823 -> 1823 | 2/1 -> 2/1 | 1 -> 1 | 6.944134e-11 |
| si16 | 602 | T | 1943 -> 1835 | 3/2 -> 2/1 | 2 -> 2 | 1.158185e-12 |
| si16 | 602 | S | 1788 -> 1788 | 2/1 -> 2/1 | 1 -> 1 | 1.089191e-10 |
| si16 | 603 | T | 1698 -> 1698 | 2/1 -> 2/1 | 1 -> 1 | 2.359002e-12 |
| si16 | 603 | S | 1668 -> 1668 | 2/1 -> 2/1 | 1 -> 1 | 3.606715e-11 |
| si32 | 601 | T | 3518 -> 3518 | 2/1 -> 2/1 | 1 -> 1 | 0.000000e+00 |
| si32 | 601 | S | 3064 -> 3064 | 1/0 -> 1/0 | 0 -> 0 | 1.044782e-10 |
| si32 | 602 | T | 3519 -> 3519 | 2/1 -> 2/1 | 1 -> 1 | 1.662670e-12 |
| si32 | 602 | S | 3077 -> 3077 | 1/0 -> 1/0 | 0 -> 0 | 1.347615e-10 |
| si32 | 603 | T | 3487 -> 3487 | 2/1 -> 2/1 | 1 -> 1 | 3.552714e-13 |
| si32 | 603 | S | 3217 -> 3217 | 2/1 -> 2/1 | 1 -> 1 | 6.393464e-11 |
| si8 | 601 | T | 829 -> 829 | 2/1 -> 2/1 | 1 -> 1 | 5.353940e-12 |
| si8 | 601 | S | 769 -> 769 | 2/1 -> 2/1 | 1 -> 1 | 8.691359e-11 |
| si8 | 602 | T | 928 -> 860 | 3/2 -> 2/1 | 2 -> 2 | 4.458656e-12 |
| si8 | 602 | S | 763 -> 763 | 2/1 -> 2/1 | 1 -> 1 | 1.326228e-10 |
| si8 | 603 | T | 830 -> 830 | 2/1 -> 2/1 | 1 -> 1 | 2.057021e-12 |
| si8 | 603 | S | 761 -> 761 | 2/1 -> 2/1 | 1 -> 1 | 1.050218e-10 |
| syn-gaussian-fcc | 601 | T | 360 -> 360 | 2/1 -> 2/1 | 1 -> 1 | 2.145839e-12 |
| syn-gaussian-fcc | 601 | S | 338 -> 338 | 2/1 -> 2/1 | 1 -> 1 | 2.201528e-11 |
| syn-gaussian-fcc | 602 | T | 417 -> 414 | 4/3 -> 3/2 | 3 -> 3 | 1.678657e-13 |
| syn-gaussian-fcc | 602 | S | 367 -> 367 | 2/1 -> 2/1 | 1 -> 1 | 4.676259e-12 |
| syn-gaussian-fcc | 603 | T | 373 -> 373 | 2/1 -> 2/1 | 1 -> 1 | 8.375522e-13 |
| syn-gaussian-fcc | 603 | S | 355 -> 355 | 2/1 -> 2/1 | 1 -> 1 | 3.236522e-11 |
| syn-p3 | 601 | T | 14 -> 14 | 1/0 -> 1/0 | 0 -> 0 | 2.131628e-14 |
| syn-p3 | 601 | S | 13 -> 13 | 1/0 -> 1/0 | 0 -> 0 | 6.052936e-12 |
| syn-p3 | 602 | T | 14 -> 14 | 1/0 -> 1/0 | 0 -> 0 | 2.486900e-14 |
| syn-p3 | 602 | S | 13 -> 13 | 1/0 -> 1/0 | 0 -> 0 | 5.857537e-12 |
| syn-p3 | 603 | T | 14 -> 14 | 1/0 -> 1/0 | 0 -> 0 | 1.776357e-14 |
| syn-p3 | 603 | S | 13 -> 13 | 1/0 -> 1/0 | 0 -> 0 | 5.911716e-12 |
| syn-spd | 601 | T | 13 -> 13 | 1/0 -> 1/0 | 0 -> 0 | 8.171241e-14 |
| syn-spd | 601 | S | 11 -> 11 | 1/0 -> 1/0 | 0 -> 0 | 1.449090e-10 |
| syn-spd | 602 | T | 13 -> 13 | 1/0 -> 1/0 | 0 -> 0 | 7.549517e-14 |
| syn-spd | 602 | S | 11 -> 11 | 1/0 -> 1/0 | 0 -> 0 | 1.463079e-10 |
| syn-spd | 603 | T | 13 -> 13 | 1/0 -> 1/0 | 0 -> 0 | 7.993606e-14 |
| syn-spd | 603 | S | 11 -> 11 | 1/0 -> 1/0 | 0 -> 0 | 1.425056e-10 |

### Calibration and adverse cases

Initial pre-failure-ratio trace: 85 actual attempts, 30 forecastable retry boundaries, 14 proposed waits. Only 6 of those 14 were failures; **8 would delay a passing attempt** (42.86% failure rate, below the 43.48% break-even proxy). This predictor was rejected, never enabled. The revised trace uses the measured post-failure change: 36 forecastable retry boundaries, 5 positive votes, **5/5 actual failures**, and 3 failures among the 31 negative votes. One other retry fell back. There are 48 unconditional first attempts. Positive-vote precision is 100% on just five correlated samples, recall on forecastable retry failures is 5/8 = 62.5%; this is not a general probability calibration. First failures are intentionally not predicted.

The revised formula was motivated by the initial silicon observations, so this is not an untouched held-out trial. Seed 601 has one positive proposal; seeds 602/603 have four, all failures. The synthetic positive proposals were not used to choose the formula. No material-dependent multiplier or fitted threshold was introduced.

| System / seed / threshold (Ha) | Candidate map | Predicted additional maps | Actual maps to pass on off trajectory |
| --- | ---: | ---: | ---: |
| si16-601-1e-10 | 11 | 1 | 1 |
| si16-602-1e-10 | 11 | 1 | 1 |
| si8-602-1e-10 | 12 | 2 | 1 |
| syn-gaussian-fcc-602-1e-10 | 10 | 1 | 2 |
| syn-gaussian-fcc-602-1e-10 | 11 | 1 | 1 |

Prediction error (predicted minus actual) on these five proposals: min -1, median 0, max +1 map; 3 exact, 1 short, 1 long. Si8 seed 602 at 1e-10 predicts two additional maps where off passes after one. On recomputes after each actual map and attempts after one: **no extra-map regression in the 48 real pairs**. For syn-gaussian-fcc seed 602, skipping map 10 removes one failed attempt but also its density reuse. The next ordinary/verification work changes; ordinary H increases **241 -> 253** while attempt H decreases **43 -> 28**; net savings are only **3 H columns (417 -> 414)**. Thus avoided attempt H is not the total benefit.

No claim is made that off and on share the same state after a skipped attempt: the default half-reuse continuation changes. The trace pass distances are counterfactual on the off trajectory only. An exact same-state oracle at every skipped on state was not run. This limitation and the zero extra-map regressions are both reported, rather than inferring universal safety of waiting.

### Provenance and validation

MgO uses the existing `ladder_handoff.fixture` geometry: 7.956 Bohr cubic rocksalt, eight atoms, Gamma, 70 Ry, PBE, Mg/O PseudoDojo files. Si and synthetic cases use checked-in benchmark cells/cutoffs/projector parameters. All atoms use `fast_screen.displaced`, seed 601/602/603, fractional amplitude 1e-3. Each reference uses identical coordinates and PP bytes. No UPF is committed.

- `target/scratch/dft-cli-baseline` SHA-256: `c0d77bcbfbc5029c797fba25f52b6e5dc7ce9625bf52f5219b2feaef19ca6dca`.
- `target/scratch/dft-cli-trace-v1` SHA-256: `01b970fcc91df95a47bdced5350e005e58fc64cbc7ea9947967209f6f461dee0`.
- `target/scratch/dft-cli-trace-v2` SHA-256: `81db8cd965a291582aad077d0e6d1571b174e28b703b03104a80dc5e5bafa17c`.
- `target/scratch/dft-cli-on` SHA-256: `3faea8caf09f890cf8ebe4704e6b8949841a528f3fd1dbbb5cc53a23fe36f319`.
- `target/scratch/rows.json` SHA-256: `c1f49c5f0b1d655c9453320b306a17de81d8afd0ef7f78e5f8264d513e040bf3`.
- `target/scratch/attempts.csv` SHA-256: `cbbfc5f7705063d5d35a8c09be91e66f50b1177dcc91f6dd0fd586f885651d3b`.
- `target/scratch/summary.json` SHA-256: `b2968dba9fbf3364c7aa2c5fe4970615bb85361ed1dd2dde8f837c01c99bdf5d`.
- External `Si.pz-vbc.UPF` SHA-256: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
- External `Mg.upf` SHA-256: `15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14`.
- External `O.upf` SHA-256: `2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11`.

- `scripts/check.sh fast dft-engine dft-numerics dft-cli`: exit **1**. Final build, fmt and clippy PASS; 50 test binaries: **460 passed, 3 failed, 134 ignored**. The only failures are the three known CLI pseudo-cache tests rejecting an in-repository TMPDIR. No tolerance or cache boundary was changed. `ps` load collection is also denied by the sandbox, but the gates still ran.
- Initial gate additionally found a clippy partially-ordered negated comparison in the new predictor; fixed before the final gate above.
- `cargo test --release -p dft-engine --lib verification_schedule -- --nocapture`: exit **0**, 2 tests pass (forecast/fallbacks and trace identity with forced failures). The final gate also covers the added real-pass assertion for on.
- `python3 -I scripts/tests/test_readme_figures.py`: exit **1**, 8 pass / 1 fail. Existing `test_structural_scope_and_sources` expects the removed phrase `standard profile is unsupported`. README, generator, test and reference data are unchanged by this patch. No new README result is generated for this rejected optimization; those files are outside OWNED_FILES.
- `git diff --check`: exit **0**. The complete `verify_fixed_point`/inner function region through `verify_reuse` is byte-identical to HEAD (SHA-256 `4d6528a178dcd5c67f232cf92af9eaf861d60bd9ac1ec9072a7360b6cda53acc`). The stopping-rule block is unchanged.
- Exact run logs, full input TOMLs, raw stdout/stderr, density bytes and per-attempt tables are under `target/scratch/{baseline,trace-study,trace-v2,references,ab}` and `target/scratch/attempts.csv`. Raw scratch is uncommitted and must be retained by the orchestrator before worktree cleanup.
- No MPI, thread-count sweep or dev08 timing claim. Spectral-audit sessions decline scheduling; spin/FD drivers are outside this session implementation. Default/invalid switch value is off. No dashboard or session registry entry for this task exists in this checkout; orchestrator-owned dashboard/README/reference JSON was not modified. No recurring orchestration timer was started by this implementer.
- Changes remain uncommitted for Claude review. No remote server changes.

## dev08 reproduction request (not executed here)

The orchestrator decides adoption after remeasuring this follow-up, including
the residual reporting cost on no-failure Si32. These commands are supplied for the
master, not authorization for this implementer to access dev08. Use the reviewed
commit, the master's quiet-node protocol (load below 0.30), and no overlapping
builds/runs. Preserve stdout/stderr, inputs, binary and UPF hashes. Compare arms
only within each `(model, case, seed, threshold)`; no QE comparison is made.
The cutoff, FFT choice, bands, starting state and binding are identical between
arms. No iteration truncation is used. Inputs contain the full displaced
structure from the checked-in benchmark, generated by the existing harness.

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export TMPDIR_990="$PWD/target/scratch"
export TMPDIR="$TMPDIR_990"
mkdir -p "$TMPDIR_990"
# Use the master's existing MKL environment/toolchain; no installs.
CARGO_TARGET_DIR=target/schedule-serial cargo build --release -p dft-cli --features mkl
CARGO_TARGET_DIR=target/schedule-mpi cargo build --release -p dft-cli --features mkl,mpi
export SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF"
python3 - <<'PY'
import csv, hashlib, importlib.util, json, os, subprocess, time
from pathlib import Path
root = Path.cwd()
spec = importlib.util.spec_from_file_location('screen', root/'scripts/fast_screen.py')
screen = importlib.util.module_from_spec(spec); spec.loader.exec_module(screen)
out = root/'target/scratch/dev08-verify-schedule'
out.mkdir(parents=True, exist_ok=False)
upf = Path(os.environ['SI_UPF']).resolve()
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
models = [
    ('onecore-si64', 'si64-gamma', 'schedule-serial', ['taskset', '-c', '0']),
    ('onecore-si32', 'si32-gamma', 'schedule-serial', ['taskset', '-c', '0']),
    ('mpi8-si64', 'si64-gamma', 'schedule-mpi',
     ['mpirun', '-np', '8', '--map-by', 'core', '--bind-to', 'core', '--cpu-set', '0-7']),
]
meta = {'commit': subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        'upf_sha256': sha(upf), 'binaries': {m: sha(root/'target'/b/'release/dft-cli')
                                           for m,_,b,_ in models}}
(out/'provenance.json').write_text(json.dumps(meta,indent=2))
with (out/'attempts.csv').open('w') as af, (out/'runs.csv').open('w') as rf, (out/'schedule.csv').open('w') as sf:
    attempts = csv.writer(af); runs = csv.writer(rf); schedule_rows = csv.writer(sf)
    schedule_rows.writerow(["model","seed","threshold_ha","pair","mode","map","predicted_n",
                            "action","would_wait","measured_coulomb_ha","ratio_to_previous_attempt","pass"])
    attempts.writerow(['model','seed','threshold_ha','pair','arm','map','H_columns',
                       'dav_iterations','tol_ha','coulomb_error_ha','pass'])
    runs.writerow(['model','seed','threshold_ha','pair','arm','exit','wall_s',
                   'verify_s','total_H_columns','attempts','failed','maps_after_first',
                   'energy_ha','density_sha256'])
    for model, case_name, build, launcher in models:
        binary = root/'target'/build/'release/dft-cli'
        case = screen.load_case(root/'benchmarks/cases'/f'{case_name}.toml')
        case['max_iterations'] = 160
        for seed in (601,602,603):
            for threshold in (5e-9,1e-10):
                src = out/f'{model}-{seed}-{threshold:g}.toml'
                src.write_text(screen.case_input(case,upf,seed,1e-3,1,
                    [f'standard_threshold_ha = {threshold}', 'verify = "fixed_point"']))
                for pair in range(5):
                    # Four timed crossed pairs, then one diagnostic trace.
                    order = ('trace',) if pair == 4 else (('off','on') if pair%2 == 0 else ('on','off'))
                    for arm in order:
                        dest = out/f'{model}-{seed}-{threshold:g}-{pair}-{arm}'
                        dest.mkdir()
                        env = {k:v for k,v in os.environ.items() if not k.startswith('BAANDI_')}
                        env['BAANDI_VERIFY_SCHEDULE'] = arm
                        command = launcher+[str(binary),'run',str(src),'--output-dir',str(dest)]
                        (dest/'command.json').write_text(json.dumps(command))
                        started = time.perf_counter()
                        with (dest/'stdout').open('w') as stdout, (dest/'stderr').open('w') as stderr:
                            rc = subprocess.run(command,env=env,stdout=stdout,stderr=stderr).returncode
                        wall = time.perf_counter()-started
                        assert rc == 0, (dest,rc)
                        r = json.loads((dest/'results.json').read_text())
                        assert r['convergence']['converged'] and r['verified']
                        v = r['verification']; aa = v['attempts']; w = r['convergence']['work']
                        if arm == 'off':
                            assert 'schedule' not in v
                        else:
                            journal = v['schedule']
                            if len(aa) == 1:
                                assert journal['prediction_evaluations'] == 0
                            for q in journal['records']:
                                schedule_rows.writerow([model,seed,threshold,pair,arm,q['map_index'],
                                    q['predicted_additional_maps'],q['action'],q['would_wait'],
                                    q['coulomb_error_ha'],q['ratio_to_previous_attempt'],q['pass']])
                            sf.flush()
                        ladder = w.get('cutoff_ladder') or {}
                        total = (w['h_vectors_recorded']+w['certificate_h_vectors']
                                 +ladder.get('coarse_h_vectors',0)+ladder.get('coarse_start_h_vectors',0))
                        runs.writerow([model,seed,threshold,pair,arm,rc,wall,v['wall_s'],total,
                            len(aa),sum(not a['passed'] for a in aa),
                            r['convergence']['iterations']-aa[0]['iteration'],
                            r['total_energy']['ha'],sha(dest/'density.bin')]); rf.flush()
                        for a in aa:
                            attempts.writerow([model,seed,threshold,pair,arm,a['iteration'],
                                a['h_vectors'],a['davidson_iterations'],a['eigen_tolerance_ha'],
                                a['coulomb_error_ha'],a['passed']])
                        af.flush()
PY
```

Run a strict reference for each identical geometry through its **same launcher**
with `BAANDI_REFERENCE_STRICT=on BAANDI_VERIFY_SCHEDULE=off`, separate from the
timed pairs, and compare delivered energy to 5e-9 / 1e-10 Ha respectively. The fifth pass above is trace (excluded from timing ratios): its energy bits, density SHA, scientific history
and counters must match off. Take median paired wall ratios within each model
and report the entire seed distribution, including regressions. MPI uses the
production owned-density path; no `M2_LOCAL_DENSITY=off` override is needed for
this scalar-history predictor. Verify rank agreement and unchanged thread-count
independence on dev08 before considering any future adoption.

## Mac harness reproduction

Save the following block as `target/scratch/study.py`; it composes the existing
`fast_screen`, `bench_suite` and `ladder_handoff` helpers. Use fresh output
folders (the harness resumes existing per-arm exit records). It records exit
codes, binary SHA-256, full inputs and all ordinary/verification history in
results.json. The complete per-attempt table is obtained with the same extraction
as the dev08 block above. Do not infer speed from this single work-count sample.

```sh
mkdir -p target/scratch
export TMPDIR_990="$PWD/target/scratch" TMPDIR="$PWD/target/scratch"
cargo build --release -p dft-cli
python3 -I target/scratch/study.py off,on,trace reproduced-ab
python3 -I target/scratch/study.py strict reproduced-references
```

```python
import importlib.util,json,os,subprocess,sys,hashlib,time,csv
from pathlib import Path
R=Path.cwd(); S=R/'target/scratch'
def module(name):
 spec=importlib.util.spec_from_file_location(name,R/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
fast=module('fast_screen');bench=module('bench_suite');ladder=module('ladder_handoff')
B=Path(os.environ.get('SCHEDULE_BINARY',str(R/'target/release/dft-cli'))); arms=sys.argv[1].split(','); out=S/sys.argv[2];out.mkdir(exist_ok=True)
cases=['si8','si16','si32','si-skew','mgo8','syn-gaussian-fcc','syn-p3','syn-spd']
if len(sys.argv)>3:cases=sys.argv[3].split(',')
for case in cases:
 for seed in [601,602,603]:
  for threshold in [5e-9,1e-10]:
   d=out/f'{case}-{seed}-{threshold:g}';d.mkdir(exist_ok=True)
   if case=='mgo8':
    fixture,_=ladder.fixture(case,seed,1e-3,R.parent/'pseudo');rows=[];atoms=[];pp=[]
    for l in fixture.splitlines():
     w=l.split()
     if w[0]=='row':rows.append(list(map(float,w[1:])))
     elif w[0]=='atom':atoms.append((w[1],list(map(float,w[2:]))))
     elif w[0]=='pp':pp.append((w[1],w[2]))
    inp='[system]\nname="mgo8"\nlattice_unit="bohr"\nlattice_vectors='+json.dumps(rows)+'\n'
    for label,a in atoms:inp+='[[system.atoms]]\nspecies='+json.dumps(label)+'\nfractional='+json.dumps(a)+'\n'
    inp+='[pseudopotentials]\n'+''.join(k+'='+json.dumps(v)+'\n' for k,v in pp)
    inp+='[basis]\necutwfc_ry=70\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional="pbe"\n[scf]\nthreads=1\nmax_iterations=160\nseed='+str(seed)+'\n'
   else:
    alias={'si8':'si8-gamma','si16':'si16-gamma','si32':'si32-gamma'}.get(case,case)
    c=fast.load_case(R/'benchmarks/cases'/f'{alias}.toml');upf=R.parent/'pseudo/Si.pz-vbc.UPF'
    if case.startswith('syn'):
     upf=S/f'{case}.UPF';upf.write_text(bench.synthetic(c))
    inp=fast.case_input(c,upf,seed,1e-3,1,[])
   inp+=f'standard_threshold_ha={threshold}\n'
   src=d/'input.toml';src.write_text(inp)
   for arm in arms:
    if arm == "strict" and threshold == 1e-10:continue
    dest=d/arm
    if (dest/'exit.json').exists():continue
    dest.mkdir(exist_ok=True)
    env={k:v for k,v in os.environ.items() if not k.startswith('BAANDI_')}
    env.update(TMPDIR=str(S),TMPDIR_990=str(S),BAANDI_VERIFY_SCHEDULE=arm)
    if arm=='strict':env.update(BAANDI_VERIFY_SCHEDULE='off',BAANDI_REFERENCE_STRICT='on')
    cmd=[str(B),'run',str(src),'--output-dir',str(dest)]
    start=time.time()
    with (dest/'stdout').open('w') as stdout,(dest/'stderr').open('w') as stderr:
     try:rc=subprocess.run(cmd,env=env,stdout=stdout,stderr=stderr,timeout=240).returncode
     except subprocess.TimeoutExpired:rc=124
    (dest/'exit.json').write_text(json.dumps({'command':cmd,'exit_code':rc,'wall_s':time.time()-start,'binary_sha256':hashlib.sha256(B.read_bytes()).hexdigest()}))
    print(case,seed,threshold,arm,rc,round(time.time()-start,2),flush=True)

```
