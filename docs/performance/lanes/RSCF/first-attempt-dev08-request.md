# I25-FIRST-ATTEMPT: first verification prediction screen

Status: experimental, `REFSCF_VERIFY_GATE_FIRST` **default off**. No remote
run or speed claim is made by the implementer. The master session owns dev08
execution, adoption, the dashboard, README evidence generation and ledger
integration. Source changes are left uncommitted for review.

Design and field definitions:
[Verification cost variants](../../../design/reference-scf.md#verification-cost-variants).
The existing empirical gate remains default on; reuse remains default off.
Do not use the historical both-off baseline of
[verify-cost-dev08-request.md](verify-cost-dev08-request.md) as arm A here.

## Stage 1: untimed counts first (required)

Build the reviewed tree once with the established dev08 release `mpi,mkl`
configuration. Preserve commit/tree and binary SHA-256, exact build command,
compiler/MPI/MKL versions and effective switch configuration. Run
`scripts/check.sh full` and the release correctness checks on dev08 only.
For identity, also build the pre-change parent with identical settings and
compare its default against candidate default: density.bin bytes, energy,
eigenvalues, history excluding clocks, evaluated maps and attempt/skip counts.
New effective-V3 audit fields are absent from default cost records; the run
configuration registry gains the default-off switch.

Freeze the four exact issue-25 inputs: **Si64, Si128, elongated (24 atoms),
rskew (2 atoms), MPI8**. Recover their inputs, UPF SHA-256, strict-reference
oracles, setup commands and mapping from `/data/sig-dft/vc/<case>/setup.md`
and the campaign manifest; record the actual resolved paths and hashes in
this campaign's manifest before running. No invented fixture or substituted
geometry. The input coordinates, cutoff, FFT grid, bands, atomic start,
threshold, iteration budget and library configuration are identical per pair.
Keep `[scf] driver='reference'`, `verify='fixed_point'`, 5e-9 Ha standard
threshold, existing Broyden controls and lightweight diagnostics.

| Arm | REFSCF_VERIFY_REUSE | REFSCF_VERIFY_GATE | REFSCF_VERIFY_GATE_FIRST |
| --- | --- | --- | --- |
| A: current defaults | off | on | off |
| B: defaults + first | off | on | on |

Use clean per-process environments. A true-default setup invocation leaves
these three switches absent and confirms the resolved values above. Explicit
values above are permitted in controlled A/A and A/B invocations. All unrelated
switches must be identical, with `BAANDI_DAVIDSON_STAGES=off`. No legacy
`VERIFY_SCHEDULE=first` or reuse experiment is selected.

Execute **one untimed complete run per arm per case**, eight runs total.
Use eight MPI ranks, one thread per rank, the same eight socket-0 physical
cores, identical mapping/binding and single-thread library settings. Preserve
`mpirun --report-bindings`, the exact CPU list and every command/exit code.
Each job has one case and a <=600 s cap; split arms into separate jobs when
needed, never extend a timed-out fast job. No detailed profiler or timing
conclusion from these runs. Example command shape (resolve placeholders from
that frozen manifest, not guessed paths):

```sh
BAANDI_REFSCF_VERIFY_REUSE=off BAANDI_REFSCF_VERIFY_GATE=on \
BAANDI_REFSCF_VERIFY_GATE_FIRST=on BAANDI_DAVIDSON_STAGES=off \
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 8 <established-binding-options> <candidate-binary> \
run <frozen-case-input> --output-dir <unique-B-output>
<candidate-binary> refscf-check --stop-to-verified <unique-B-output>/results.json
```

For **every actual attempt** retain iteration, all three prediction fields,
actual Coulomb error (Ha), occupied residual bound (Ha), eigen tolerance (Ha),
threshold (Ha) and pass/fail. Undefined predictions are null, never zero. A's
first estimate is null by construction; B records the uncalibrated estimate.
Record each skip as well. The first skipped candidate has no actual error in
B: if A evaluated the identical candidate, report A's value as a separately
labelled counterfactual only after confirming identical preceding science.
The reader's `verification_predictions` array includes eligible B candidates;
the raw verification report contains only actual attempts.

For each arm report logical iterations, ordinary/verification requests,
evaluated maps including retries, solver calls, attempts, skips, reuse (zero),
Davidson corrections and H columns. Verify counts against raw call ranges.
Preserve outputs including density.bin, eigenvalues and both history copies;
check interval closure with `refscf-check --stop-to-verified`. A first skip
lies in the existing unassigned part until an actual failure; never fabricate
a failed attempt or move that work outside delivered wall.

**Advance to timing only if BOTH Si64 and Si128 have strictly fewer total
evaluated maps in B than A, all four cases converge through an unchanged
passing attempt, and no correctness gate fails.** An attempt reduction alone
is insufficient. If either target is equal/worse, stop after the untimed
screen, keep the switch off, and return all counts/predictions as a negative
or inconclusive result. Do not increase the skip limit or tune constants to
those cells. Inspect elongated/rskew regressions before any adoption decision.

## Stage 2: conditional A/A and crossed A/B

Only after Stage 1 qualifies, follow the bounded-job, external-wall and
statistics protocol of `verify-cost-dev08-request.md`, using the new A/B arms:

- 20 A/A pairs per case to establish the noise floor, alternating labels.
- 30 crossed A/B pairs per case, alternating A-B and B-A, with A/A controls
  in the same windows. Split into <=600 s, one-case jobs; Si128 starts with
  one pair/job. Warmups and controls count against each job's cap.
- Time the entire `mpirun ... run ... --output-dir ...` externally, through
  stdout/stderr, all files, JSON and MPI exit. Include every verification,
  predictor and recovery cost. Report medians consistently, paired log effects,
  SE and the A/A floor. Use the existing campaign resolution rule
  `abs(mean)/SE >= 3` and `abs(mean) >= 3F`; no minimum percentage gain.
- Retain failed, exhausted and censored runs. Stop on a correctness failure;
  uncertainty inside noise is inconclusive, not a win. Compare only MPI8 with
  MPI8; no cross-model table, ratio or conclusion, and no QE speed claim here.

Every setup/timed run must pass the unchanged independent verification and
have absolute energy error <5e-9 Ha against that geometry's converged strict
oracle (strict accuracy class 1e-10 Ha). Preserve strict-oracle command,
input/UPF/binary hashes, tolerances and certificate. Agreement of A and B
alone is not accuracy evidence. Preserve rank-count determinism checks under
the existing reference-driver protocol and report known baseline issues
separately; never relax gates. Widen to the established non-Si/projector
correctness fixtures before any promotion, since the four named cells are Si.

The ignored displaced-Si8 test now labels explicit both-off as `off` and
includes true `default` arms with no verification-switch override plus `first`.
Run only on dev08, preserving source/binary/input/UPF hashes and stdout:

```sh
DFT_SI_UPF=/absolute/path/to/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf stop_to_verified_si8_displaced -- --ignored --exact --nocapture --test-threads=1
```

## Local evidence and limitations

Mac work is limited to the ticket's compile, small debug correctness tests,
format and lint commands. No full gate, release matrix, remote execution,
UPF installation, default promotion or timing adoption occurred. The master
must integrate this report into the ledger/dashboard/README after review;
those files are outside this ticket's owned paths.

### Synthetic Mac observations (2026-11-24)

Parent source: `f8d14c281eb2468839e6f9c158db8b6d1f4d367a`. The two fixtures
are unchanged from the previous verification-cost tests: weakly skew cells,
2 Ha cutoff, beta=0.7, synthetic atomic-start Gaussian pseudopotentials. The
failed-first fixture has z=4, rc=1.2 Bohr, 11 Bohr diagonals, 3 bands, history=4;
the first-pass control has z=2, rc=1.8 Bohr, 7 Bohr diagonals, 2 bands,
history=8. No physical material or QE reference is claimed.

| Failed-first fixture arm | Logical iterations | Evaluated maps | Attempts | Skips | Reuse | Delivered energy (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| default | 28 | 30 | 2 | 1 | 0 | -2.61337340593422773 |
| explicit FIRST=off | 28 | 30 | 2 | 1 | 0 | -2.61337340593422773 |
| default + FIRST=on | 27 | 28 | 1 | 1 | 0 | -2.61337340577661736 |
| reuse/gate off | 27 | 29 | 2 | 0 | 0 | -2.61337340577661736 |
| FIRST=on, GATE=off | 27 | 29 | 2 | 0 | 0 | -2.61337340577661736 |
| reuse on, gate off | 29 | 30 | 2 | 0 | 1 | -2.61337340649649974 |
| explicit gate on | 28 | 30 | 2 | 1 | 0 | -2.61337340593422773 |
| reuse and gate on | 29 | 30 | 2 | 0 | 1 | -2.61337340649649974 |

All eight arms of the first-pass control use 6 logical iterations, 7 evaluated
maps, 1 attempt, 0 skips/reuse, and energy -1.10503701485504457 Ha.
There are no retries/redos in these fixtures; evaluated maps equal solver
calls. Each delivered state passes the unchanged 5e-9 Ha criterion and
occupied-residual guard. FIRST-on minus default energy is 1.576103692e-10 Ha
on the active fixture and 0 Ha on the control; agreement is not a strict oracle.

| Fixture / candidate | Iteration | Secant prediction (Ha) | Mixed contraction (Ha) | Actual error (Ha) | Action |
| --- | ---: | ---: | ---: | ---: | --- |
| Failed-first, initial candidate | 23 | 5.504054419e-9 | 1.660904272e-9 | null in B; 1.783487664e-8 in A | skip in B; fail in A |
| Failed-first, next eligible candidate | 27 | 1.755352352e-8 | 2.043335126e-9 | 4.378277166e-11 | mandatory attempt in B, passes |
| First-pass control | 6 | 3.951894732e-10 | 5.327520181e-11 | 1.649390463e-10 | attempt, passes |

A's iteration-23 error is a baseline counterfactual for the skipped B
candidate, not a measured B verification. A's iteration-27 candidate is
skipped by its already calibrated retry gate; A passes at iteration 28 with
actual error 5.896222250e-10 Ha (prediction 6.029534067e-9 Ha). These facts
explain the 30 -> 28 evaluated-map reduction, not just an attempt count change.
B's overprediction on its *first actual attempt* is roughly 401 times the
actual error, explicitly illustrating the refresh safety requirement.

Prototype comparison (same first stopping candidates, test logs retained in
`target/tmp/i25-engine{,-inverse,-amplification}.log`): direct forward
Jacobian fit predicted 6.375661735e-10 / 5.666906701e-11 Ha on active/control;
inverted inverse-Broyden predicted 1.504757613e-10 / undefined. Neither skipped
the first failure; both retained 30 maps on the active fixture. The selected
observed amplification rule needs no inverse or fitted safety multiplier and
has the counts above. The mathematical contraction limitation and prototype
failures are documented in the design section; no claim that the selected
estimate is an error bound follows.

The analytic scalar-linear oracle checks the amplification exactly, and checks
that prediction leaves subsequent mixing/gamma/history unchanged. Additional
tests cover ring eviction, raw secants under plain/Thomas-Fermi mixing,
no history, zero secant, nonfinite values, absent/disabled first switch,
GATE-off inertness, forced refresh, unchanged ratchet and pass judgment, and
budget exhaustion on a first skip without any actual attempt. Existing Si8
expectations remain dev08-only and are not newly reproduced here.

CLI roundtrip checks cover both synthetic fixtures, seven switch arms and
lightweight/full diagnostics (28 complete CLI runs), both reference-history
copies, counts against raw call records, predictions versus attempt reports,
interval reconstruction and the reader's prediction rows. Default JSON retains
its five cost fields; effective FIRST adds three fields and the summary flag.


Default identity: a temporary harness ran the parent driver's exact default
algorithm and original `Metric::batch` from the parent commit above against
the candidate, on both synthetic fixtures. Only additive record-construction
fields and unused new helper symbols were supplied to compile the common test.
Complete timing-stripped snapshots (energy, density/orbital bits, eigenvalues,
occupations, history, mixing/work records and verification verdicts) matched
exactly for both fixtures. Both filtered engine test commands exited 0;
`target/tmp/i25-identity-{parent,candidate}.log` and
`target/tmp/i25-identity-results.json` retain the local evidence. Temporary
source/test changes were restored in a `finally` block before the final gates.

### Final authorized Mac gates

| Command | Exit |
| --- | ---: |
| `cargo check --workspace --tests` | 0 |
| `cargo test -p dft-engine --test refscf -- --nocapture` | 0 |
| `cargo test -p dft-cli stop_to_verified -- --nocapture` | 0 |
| `cargo test -p dft-numerics switches` | 0 |
| `cargo fmt --check` | 0 |
| `cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings` | 0 |
| `git diff --check` | 0 |

Engine: 16 passed, 1 ignored (external Si8), 0.21 s test execution.
CLI: 8 selected tests passed, 0.98 s; numerics: 8 unit plus 6 integration
switch tests passed, registry count 136. Times establish the debug-test budget,
not a performance claim. Final logs: `target/tmp/i25-verified-*.log`; command
exit codes: `target/tmp/i25-verified-commands.json`. An intermediate formatting
check flagged a newly added test; `cargo fmt --all` fixed it before these final
checks. No final build, test, formatting or lint gate is failing.

Remote correctness, strict-oracle accuracy, MPI scheduling/determinism and
speed remain unmeasured. No tracked source/document paths outside ticket ownership were changed;
README, JSON registries/dashboard and ledger integration remain master-owned.
