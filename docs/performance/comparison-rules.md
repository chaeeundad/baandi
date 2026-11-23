# Speed comparison rules (user rule, 2026-11-14; binding for every later task)

## Current state — 2026-11-20..24 (issue #26)

The new reference SCF path (issue #25) is the active work; the threads lane is
dropped. Defaults are unchanged by this reconciliation. dev08 is operated by
`baandi-dev08`; the `dft5090` relay is retired. See the
[current operations and protocol map](fast-loop-protocol.md#current-state--2026-11-2024-issues-25-and-26).
Use fixed-H replay for kernels, common-input replay for density-update
discrepancies, and full convergence for tail, verification or trajectory changes.
Timing uses [the same median for both engines](comparison-rules.md#statistics-and-delivered-wall);
fast jobs are limited to 600 s, authorized `final-` jobs to 7200 s.

## Comparable conditions (live)

1. Compare one core with one core, or BAANDI MPI with QE MPI at the same rank
   count, one OpenMP/MKL thread per rank and the same physical cores/binding.
   The threads lane is dropped. If reopened by a later ticket, threads must be
   compared only with QE OpenMP at the same thread count and one MPI rank.
   Never include mixed-model columns, even as information; name each model in
   every table header, ratio, ledger conclusion and report.
2. Match cell/coordinates, cutoff, FFT grid, UPF bytes (SHA-256), band count,
   occupations, XC, start and threshold definition. Record any mismatch as a
   diagnostic limitation, not an equal-conditions speed claim. Use the fastest
   eligible QE configuration for that cell/model under the matched conditions.
3. Keep strict 1e-10 Ha and standard 5e-9 Ha delivered-accuracy classes against
   independent tight references. Record QE's actual `conv_thr` and delivered
   error (the existing comparison pairing is standard with `conv_thr=1e-8` Ry,
   strict with `conv_thr=1e-10` Ry); a printed threshold or rounded energy
   agreement alone is not proof of the class. BAANDI independent fixed-point verification remains required.
4. These rules apply to exploration, promotion, final tables, README claims,
   ledger rows and every dev08 request. Historical tables/requests below or
   elsewhere retain provenance; their superseded instructions do not authorize
   a new mixed-model comparison.

## Diagnostics level: QE comparison versus hypothesis A/B (user rule, 2026-11-23)

- A run timed against QE uses fully identical conditions on both sides, including the
  output level: BAANDI runs with detailed profiling OFF (production output only: required
  results, lightweight counters, the verification record and completion marker). Inner
  clocks, per-band Ritz traces and extra trace JSON are never on in a QE-comparison run.
- Detailed profiling is a tool for testing an improvement hypothesis only. In such an A/B
  both arms (before and after the change) run with it ON, so its cost is the same in both;
  these numbers are never put in a QE comparison table.
- The fixed-point verification is a production step (decision 2026-11-14), not profiling;
  it stays in the timed QE-comparison run and is charged as described below.

## Statistics and delivered wall

Use **the median for both BAANDI and QE**, with the same repetition/geometry
protocol and at least three samples per condition (five for a close call).
Report every raw sample, sample count, spread, failures, input/UPF hashes,
command, source/binary hashes and measurement window. The final ratio is
`median(BAANDI delivered wall) / median(QE delivered wall)` for each matched
case/model. A mean from one engine must never be paired with the other's median.
Do not pool distinct cases or geometry ensembles without labelling the protocol.

For candidate A/B, measure the A/A floor first, use crossed order-alternated
pairs within one window, and report per-arm medians plus paired differences,
spread and uncertainty. Paired mean/SE diagnostics may assess resolution;
they do not replace the median headline statistic. Preserve noisy/failing
batches and report repeats; do not select the fastest samples.

BAANDI delivered wall is the external elapsed time to the delivered verified result,
including setup, SCF, **all verification attempts (passed and failed), failed
attempts and recovery**, and requested output. If a recovery restarts the
process, charge those elapsed attempts too. Report verification as its own
component without subtracting it from total wall. Unrecovered failures and
timeouts remain failures, never omitted to manufacture a successful median.
Use the same external start-to-output wall boundary for QE, including any failed
attempts/recovery; report its internal timers only as
separate diagnostics, not as a different denominator. `verify=off`, fixed-H
and truncated screens cannot establish delivered-result speed acceptance.

See [protocol selection and promotion](fast-loop-protocol.md#protocol-map-live)
and [final tables](final-measurement-protocol.md). No default changes or new
numerical evidence are implied by this rule reconciliation.

<details>
<summary>Historical — 2026-11-14 comparison-rule wording and rationale; superseded, not live instructions</summary>

1. **Compare like with like in the parallel model.** Multithreading is compared with multithreading, MPI with MPI. Never put BAANDI threads against QE MPI ranks (or the reverse) in a headline column, a ratio, a ledger conclusion or an adoption argument.
   - Thread against thread: BAANDI `threads = N` (shared memory, one process) against QE with `OMP_NUM_THREADS = N`, one MPI rank (`mpirun -np 1`), MKL threads as QE's OpenMP build uses them; same physical cores, same socket binding.
   - MPI against MPI: BAANDI built with `--features mpi` (M2 path, `mpirun -np N`, threads 1) against QE `mpirun -np N`, OpenMP and MKL threads 1; same ranks, same cores.
   - One core is both models at once (BAANDI T=1, QE np=1, OMP 1) and may always be compared.
2. **Every table states the model of each column** ("threads vs threads", "MPI vs MPI") in the header. A column that mixes models is not allowed, not even labelled as informational (tightened 2026-11-14). Existing mixed columns are removed from the documents; the raw numbers stay in the ledger only as single-model measurements.
3. **Equal conditions beyond the model:** same cell, cutoff, FFT grid, UPF (record SHA-256), band count, start (QE `startingwfc = 'atomic'`), same displaced coordinates from one generator, accuracy class matched (BAANDI standard against QE `conv_thr` 1e-8, strict against 1e-10), idle node, 3 or more repetitions with the median reported, command and binary hashes recorded.
4. **Consequence for existing tables:** the "T=8 against np=8" columns of `qe-vs-baandi-2026-11-11.md` and `qe-vs-baandi-2026-11-13.md` mix models and are removed from the comparison documents. The valid 8-core comparisons are (a) BAANDI threads against QE OpenMP (needs the QE OMP measurement for every cell; only si32 strict unperturbed exists: QE OMP 8 = 1.957 s) and (b) BAANDI MPI (M2d, np=8) against QE MPI np=8.
5. **Where the rule applies:** QE comparison tables, README and report statements about being faster or slower than QE, ledger rows that quote a QE ratio, the performance target in AGENTS.md ("match and then exceed QE at equal core count"), and every dev08 request that involves QE.

Reason: QE's OpenMP scaling (1.78x at 8 threads for si32 strict) is far below its MPI scaling (4.18x), so comparing across models misstates the gap (si32 strict 8.74x threads-against-MPI against 3.84x threads-against-threads).

</details>
