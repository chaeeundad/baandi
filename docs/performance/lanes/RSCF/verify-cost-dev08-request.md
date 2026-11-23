# I25-VERIFY-COST: reference-driver verification tail

Status: experimental, both switches default off; no dev08 result or speed
claim. The master session owns remote execution, promotion, dashboard,
README evidence generation and ledger integration. This worktree is left
uncommitted. Design: [Verification cost variants](../../../design/reference-scf.md#verification-cost-variants).

## Premise and prior attempts

User-relayed issue 25 measurements at main `5920aee`, dev08 MPI8: Si64 first
stop 11, pass 13, three verification errors 5.80e-8 / 6.75e-9 / 2.09e-9 Ha;
verification 19% and subsequent recovery 9.2% of delivered wall. Si128 first
stop 15, pass 17, errors 1.51e-7 / 3.09e-8 / 1.47e-9 Ha; verification 17%
and recovery 7.7%. These are motivation supplied in the ticket, not newly
reproduced measurements; obtain original raw evidence before adoption.

Earlier legacy `VERIFY-REUSE-FEED`/`FEED2`/`REPLACE` trials regressed some
Si16 cases; legacy `half` was later adopted (attempts ledger). That does not
establish a win for reference Broyden. The changed premise is the new
reference driver's discarded full verification maps and its measured
three-attempt tails. V1 consumes exactly those maps as ordinary Broyden
iterates; no claim is inherited from the legacy mixer. V2 uses measured
error/estimator scaling, not the legacy schedule or a fitted silicon factor.
Both are trajectory experiments and require full converged runs.

## Arms and identity

Build reviewed candidate once on dev08 with the same release MPI/MKL features
and toolchain as the recorded main baseline. Record baseline, current-main
and candidate source commit/tree hashes, complete build commands, binary
SHA-256, compiler/MPI/MKL versions and effective run configuration. Perform
`./scripts/check.sh full` and release correctness tests there, never on Mac.
A pre-change main/candidate-default identity check must compare energy,
density.bin, eigenvalues, history excluding clocks, maps and attempts.

| Arm | REFSCF_VERIFY_REUSE | REFSCF_VERIFY_GATE |
| --- | --- | --- |
| A | off | off |
| B1 | on | off |
| B2 | off | on |

Set explicit environment values `BAANDI_REFSCF_VERIFY_REUSE` and
`BAANDI_REFSCF_VERIFY_GATE` for every invocation. Use `[scf] driver='reference'`,
`verify='fixed_point'`, unchanged standard 5e-9 Ha threshold and reference
Broyden controls. `BAANDI_DAVIDSON_STAGES=off` (lightweight); no detailed
profiler in acceptance samples. Combined B1+B2 is a separate correctness
check, not a substitute for either independent A/B.

All arms use identical cell/coordinates, cutoff, FFT dimensions, UPF bytes,
bands, atomic start, iteration limit, thresholds, eight MPI ranks, one thread
per rank, socket-0 core binding and library-thread settings. Capture input
and UPF SHA-256 plus `mpirun --report-bindings` output. Use the same eight
physical cores for A/A and A/B; record the actual CPU list and MPI mapping.
No threads-versus-MPI comparison or cross-model ratio is permitted. This
request compares BAANDI arms only, with strict reference as an accuracy oracle.

## Cases and bounded jobs

Run Si64 and Si128 MPI8 from the issue-25 fixture manifest, plus the established
skewed and slab-like generality cells from the reference SCF/benchmark suite.
Freeze the exact four inputs and hashes before measuring. If a named slab
fixture has no converged strict oracle, establish one first; do not silently
substitute a cubic cell or invent reference energy. Do not use an open-shell
fixture with the fixed closed-shell driver. Include a non-Si synthetic
projector fixture in the separate correctness widening before promotion.

Each job is **one case, MPI8, at most 600 seconds**, including all warmups,
pairs, JSON/density output and analysis. First a disposable warmup and an A/A
pair establish runtime; choose a small pair batch that fits the remaining
budget. For Si128 start with one pair/job. Accumulate samples across bounded
jobs in the same controlled campaign; never fit a matrix into one job or
extend a timed-out job beyond the cap. Failed, censored and exhausted samples
are retained and reported, not discarded as timing outliers.

1. A/A: both invocations A, alternating labelled order, 20 pairs/case total.
   Store paired log differences and the measured noise floor/uncertainty.
2. Crossed A/B1: alternate A-B1 and B1-A; 30 pairs/case total. Repeat A/A
   controls in the same windows. Restore all switches explicitly per process.
3. Crossed A/B2: same independently, never reuse B1 as the baseline.
4. Stop on a correctness failure, preserve evidence and diagnose. A runtime
   effect inside the A/A uncertainty is inconclusive; no minimum percentage
   gain is imposed. Any reproducible gain with unchanged accuracy and no
   regression beyond noise is eligible for the orchestrator's adoption review.

Time the whole `mpirun -np 8 ... dft-cli run INPUT --output-dir UNIQUE_DIR`
process externally using the established runner or `/usr/bin/time`, including
stdout/stderr, all requested files, final JSON and MPI exit. Preserve external
wall per invocation and exit code. Never use only engine or solve time as the
acceptance denominator. Use medians consistently for A and B; report median
ratio, paired log-effect and SE, A/A floor and raw samples. Timing-job setup,
builds and strict-oracle creation are separate untimed jobs.

## Required evidence per run

- Completed result, convergence/pass verdict, energy in Ha and absolute error
  against that input's converged strict reference. Require error <5e-9 Ha;
  retain B-minus-A energy too. Preserve strict input, UPF, binary/commit,
  command, tolerance, certification and energy provenance (strict class
  1e-10 Ha). No accuracy claim based only on two loose arms agreeing.
- Every actual attempt: iteration, next-ratchet tolerance in Ha, predicted
  versus actual Coulomb error in Ha, residual estimate, threshold, pass/fail.
  Every skip and reused source iteration; never count a skip as a failure or
  a pass. Undefined predictions remain null, not zero.
- `refscf-check --stop-to-verified RESULTS_JSON`: first stop, passing map,
  interval seconds, verification seconds, recovery seconds, unassigned seconds,
  closure check, actual attempts and skipped count. Preserve the map timestamps
  so the interval is reconstructible. Report each part as seconds and as a
  fraction of the *same run's* external wall.
- Counts: logical SCF iterations; reused iterations; ordinary map requests =
  logical iterations minus reused plus first-map redo; verification requests =
  actual attempts. Evaluated ordinary/verification maps are the lengths of
  their solver-call records, including bounded retries; total evaluated maps
  equals total solver calls. Also report Davidson corrections and H columns.
  Reused maps have zero new calls: their work belongs to the source attempt.
- Density.bin/hash, eigenvalues and complete lightweight history; final state
  must be the candidate checked by the final passing attempt. On representative
  cases additionally compare rank-count correctness using the existing
  reference-driver determinism protocol; known Si16 rank-3 threshold sensitivity
  remains a separate unresolved baseline issue, not permission to relax gates.

Publish one table per case with A/A floor, A/B medians, uncertainty, tail
parts, maps/attempts/skips/reuse, strict-reference errors and failures. Store
raw per-run outputs externally with a manifest of paths/hashes. Add the final
attempts-ledger/optimization-log entry through the orchestrator-owned files,
including rejected/inconclusive results. No default is changed by this request.

## Local validation

Only the ticket's Mac checks are authorized. The small synthetic skew-cell
fixture tests default/off identity, V1, V2 and both, unchanged pass and ratchet,
energy class, counters and clock closure. The analytic predictor test covers
undefined/nonfinite/noncontracting inputs, threshold equality and mandatory
refresh. CLI tests cover explicitly skipped candidates and invalid evidence.
The displaced-Si8 test is ignored: external UPF and full-SCF comparisons belong
on dev08, using the documented release command in the test.

Local measured results and any remaining integration limitation will be
appended after the checks finish. Mac wall observations are not performance
evidence and no release matrix is run here.

### Mac observations (2026-11-24)

`cargo test -p dft-engine --test refscf -- --nocapture`: 14 passed,
1 ignored; test execution 0.22 s (debug profile with the repository's opt-level).
The new tests use the existing synthetic atomic-start Gaussian PP helper,
2 Ha cutoff and a weakly skew cell. Active fixture: z=4, rc=1.2 Bohr,
11 Bohr diagonals, 3 bands, beta=0.7, history=4. Negative-control fixture:
z=2 local-only, rc=1.8 Bohr, 7 Bohr diagonals, 2 bands, history=8.
These are correctness/control-flow fixtures, not physical material references.

| Active fixture arm | Logical iterations | Evaluated maps / solver calls | Attempts | Reused | Skipped | Delivered energy (Ha) | Delta from A (Ha) | Passing Coulomb error (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Default / explicit off | 27 | 29 | 2 | 0 | 0 | -2.61337340577661736 | 0 | 4.378277166e-11 |
| Reuse | 29 | 30 | 2 | 1 | 0 | -2.61337340649649974 | -7.198823759e-10 | 1.458563922e-11 |
| Gate | 28 | 30 | 2 | 0 | 1 | -2.61337340593422773 | -1.576103692e-10 | 5.896222250e-10 |
| Both | 29 | 30 | 2 | 1 | 0 | -2.61337340649649974 | -7.198823759e-10 | 1.458563922e-11 |

No redo/retry occurs in these fixtures, so evaluated maps equal solver calls.
The active fixture **increases total work** under both experiments; the gate
skips what would have been a passing default attempt. This falsifies any
universal map-saving claim and motivates the crossed multi-system experiment.
All energies stay within 5e-9 Ha of A; no strict-oracle accuracy is established
by this local test. All arms of the first-pass negative control take 6 SCF
iterations, 7 maps, 1 attempt, no reuse/skip, and deliver
-1.10503701485504457 Ha with passing error 1.649390463e-10 Ha.
Default and explicitly disabled variants have bit-identical scientific
snapshots (density/orbitals/eigenvalues/energy/history/work, clocks excluded).
The active test also exhausts the budget exactly on a failed attempt and on
a skipped candidate; neither returns a converged/passing result.

### Adapter integration (approved by orchestrator)

The orchestrator extended OWNED_FILES to `crates/dft-cli/src/run.rs` for
`reference_record_json` and the `reference_scf` summary only. Integration is
complete: the five per-iteration `verification_cost` fields appear in both
history copies when either switch is enabled; both-off rows retain their
previous shape. The summary records attempts, skips, reuse, logical SCF
iterations, grouped map requests, evaluated maps and solver calls. Each
`map::evaluate` performs one solver call; retries count as additional evaluated
maps, while grouped requests exclude bounded re-entry. First-step redo is a
separate request/evaluation; reuse incurs neither a new request nor evaluation.

The CLI roundtrip tests run default, reuse, gate and both on the existing small
first-pass fixture and the active-failure fixture, at lightweight and full
diagnostic levels. They check all five fields, both history copies, summary
counts against raw call records, predicted versus measured error, actual
nonzero skip/reuse, and `refscf-check --stop-to-verified` output/clock closure.
No serializer gap remains. The strict reader accepts missing attempt clocks
only with an explicit valid skipped record; contradictory evidence still fails.

### Displaced Si8 per-variant report (dev08 release only)

The test remains `#[ignore]`: it needs an external Si UPF and runs full SCFs.
No new Si8 count is claimed from Mac. Execute this exact command on dev08 with
the established UPF path substituted; preserve UPF SHA-256, source/binary hashes,
exit code and stdout as the count report:

```sh
DFT_SI_UPF=/absolute/path/to/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf stop_to_verified_si8_displaced -- --ignored --exact --nocapture --test-threads=1
```

One case only, within a <=600 s job. The command prints one `STOP_TO_VERIFIED
si8` row for each arm below (default has lightweight/full rows). Each row
contains `variant`, `passing` (logical maps), `evaluated_maps` (all physical
map evaluations including retries), `map_requests`, `attempts`, `skipped`,
`reused`, energy/delta in Ha and interval parts in seconds. It asserts the
unchanged pass/ratchet criteria and the 5e-9 Ha default-energy comparison.

| Variant | Logical maps | Evaluated maps | Attempts | Skipped | Reused |
| --- | --- | --- | --- | --- | --- |
| default | 12, existing fixture assertion | pending execution | 3, existing fixture assertion | 0 by construction | 0 by construction |
| reuse | pending execution | pending execution | pending execution | 0 by construction | pending execution |
| gate | pending execution | pending execution | pending execution | pending execution | 0 by construction |
| both | pending execution | pending execution | pending execution | pending execution | pending execution |

The baseline assertions are inherited expectations, not newly reproduced
measurements. Fill the pending cells only from the command's actual output;
do not substitute the synthetic fixture counts above or the ticket's MPI8
Si64/Si128 counts. No implementer remote execution is performed.

Local command results (initial implementation; refreshed after adapter integration below):

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace and test compilation passes. Initial invocation found a new test calling nonexistent `RunConfig::default`; corrected to environment-free `from_pairs(&[])` and rerun. |
| `cargo test -p dft-engine --test refscf -- --nocapture` | 0 | 14 passed, 1 ignored (external Si8). |
| `cargo test -p dft-cli stop_to_verified` | 0 | 7 selected interval tests passed; current default CLI roundtrip at both diagnostic levels passes. |
| `cargo test -p dft-numerics switches` | 0 | 8 unit plus 6 integration tests passed; registry count 135. |
| `cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings` | 0 | Formatting and lint gates pass. |
| `git diff --check` | 0 | No whitespace errors. |

No full gate, release test/matrix, remote command, UPF installation, commit,
README regeneration or orchestrator-owned JSON/dashboard edit was performed.
The narrow ticket ownership and Mac-command list take precedence over the
repository's broader default reporting/validation workflow.

Historical-default identity also passed (both invocations exit 0): temporarily
build the pre-change driver from HEAD
`4d3580e55ca83625f1919ae7f2dc9a10af9b5018`, adding only `verification_cost: None`
to construct the extended record and the unused predictor API needed to
compile the tests; compare with the candidate using
`cargo test -p dft-engine --test refscf verification_cost_ -- --nocapture`.
For this comparison only, restrict the fixture loop to default/explicit-off
and print the complete timing-stripped scientific snapshot. Both fixtures
match exactly, including density/orbital bits, eigenvalues, energy, history,
map/call records and verdicts. Temporary sources/tests were restored in a
`finally` block. Local raw snapshots are in
`target/tmp/verify-default-{baseline,candidate}.log`; no baseline algorithm
was modified and no performance ratio is inferred. The full final tests above
retain all five arms, active-branch assertions and budget-exhaustion checks.

### Refreshed checks after approved adapter integration

| Command | Exit |
| --- | ---: |
| `cargo check --workspace --tests` | 0 |
| `cargo test -p dft-engine --test refscf -- --nocapture` | 0 |
| `cargo test -p dft-cli stop_to_verified -- --nocapture` | 0 |
| `cargo test -p dft-numerics switches` | 0 |
| `cargo fmt --check` | 0 |
| `cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings` | 0 |

Engine: 14 passed, 1 ignored (displaced Si8), 0.14 s test execution.
CLI interval: 8 passed, 0.82 s; includes 16 actual CLI runs across two
fixtures, four switch arms and two diagnostics levels, with separate reader
processes. Numerics switches: 8 unit + 6 integration passed. Times only
establish the debug-test budget, not performance benefits.

The active CLI fixture matches the engine counts at both diagnostic levels:
reuse = 30 evaluated maps / 2 attempts / 0 skips / 1 reuse; gate =
30 / 2 / 1 / 0; both = 30 / 2 / 0 / 1. The default engine control remains
29 / 2 / 0 / 0. All gated CLI outputs are accepted by the independent interval
reader with closure checked. The two history copies agree on all five new
fields, actual errors match the verification report and summary counts match
the raw solve records. No CLI metadata fields are added for both-off runs.

This closes the earlier ownership/serialization blocker. Displaced-Si8
per-variant measurements and all dev08 timing/strict-oracle gates remain
explicitly deferred to the command and protocol above; no remote results
are inferred from this local completion. All changes remain uncommitted.
