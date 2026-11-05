# Study note 3: P14A Davidson work accounting (2026-11-05)

Outcome: **no numerical default or optimization adopted**. The fixed-factor
experiment does not preserve SCF iteration count on any of the four systems.
The conditional gate for further H6/H2 implementation is not met. H7 also fails
the generality and same-result requirements. This is a negative experimental
result, not a claim that every possible adaptive schedule must fail.

Only rust-dft, existing study notes and QE outputs were read; no QE source was
read. Local source is `a1c6644d59ed58d1e89b610716343bdf696bfc70`, pure Rust,
Apple Silicon, release, one thread unless specified. Binary and UPF hashes are
stored with every experiment. Timings overlap other validation/build work and
are observations only, not speedup evidence. No new sampling profile or accepted
performance ratio is claimed; the relayed dev08 profile motivated this work.

## Evidence and measured bound

[Complete tables](../../../tests/reference/perf/p14a/work-tables.md) include
all SCF iterations, thresholds, Davidson steps, H block calls, H vector counts,
restarts, early/late convergence, buffer-origin fraction, SCF phase fractions,
QE alignment and every rejected variant. The adjoining JSON files are the
machine-readable evidence. [The summary generator](../../../tests/p14a_summarize.py)
checks the independent accounting identities before generating those tables.

Baseline totals, with the final certificate included in H counts but separate
from ordinary Davidson steps:

| Case | SCF count | Ordinary Davidson steps | H blocks including certificate | H vectors including certificate | Restarts |
| --- | ---: | ---: | ---: | ---: | ---: |
| si8 Gamma | 15 | 340 | 410 | 2416 | 14 |
| k444s | 10 | 1686 | 1891 | 11181 | 0 |
| synthetic 20-atom skewed | 6 | 13 | 17 | 255 | 0 |
| synthetic 4-atom elongated | 9 | 28 | 33 | 171 | 0 |

The two synthetic fixtures use the existing soft potential generator, valence
2, channels `[0,0,1]`, and cutoff 3 Ry. Their lattice and atom coordinates are
in the inputs and runner; they are synthetic numerical cases, not material
predictions or QE references.

For factors 1, 3, 10, 30 and 100, the si8 SCF counts are respectively
15, 16, 16, 21 and 57; ordinary Davidson steps are 340, 242, 264, 278 and 288;
total H vectors are 2416, 2072, 1794, 1722 and 2280. A work reduction obtained
by increasing SCF iterations does not meet this ticket's acceptance rule.
The corresponding k444s SCF counts are 10, 12, 15, 35 and 58. The synthetic
20-atom case diverges at factors 3 and 10 after six completed SCF iterations
(13 Davidson steps and 174/170 H vectors respectively); factors 30 and 100
increase its SCF count. The elongated case increases its count for every
factor above one. Exact energies and all counts are in the complete table.

Thus at unchanged SCF count, the tested family leaves the full baseline work:
340/1686/13/28 ordinary Davidson steps and 2416/11181/255/171 total H vectors.
No factor above one qualifies. Equal iteration count would itself be weaker
than identical SCF trajectories; we do not infer trajectory identity merely
from that count. The factor-one runs match the original baseline energies and
work exactly. No relaxed candidate is accepted as a strict-reference result.

## Experimental definition and certificate protection

The actual schedule is in `crates/dft-engine/src/session.rs`. The original
schedule uses `min(0.1 * previous_density_L2_residual, 1e-2)`, with floor
`0.1 * eigen_tolerance`. The experiment multiplies the residual-dependent term
before those bounds. It preserves the initial 1e-2 Ha residual tolerance and
the 1e-9 Ha floor. Multiplying that floor too would prevent factors 30 and 100
from ever meeting the existing ordinary convergence gate, which requires a
requested tolerance at most 1e-8 Ha.

The exact temporary change is retained as
[tolerance-factor.patch](../../../tests/reference/perf/p14a/tolerance-factor.patch).
It was applied to the actual schedule, so history records the tolerance used.
The fixed-point certificate, its invocation criterion, and all its tolerances
were unchanged. After rejection the patch was removed: the final session file
is byte-identical to the checkout. No hidden compile-time numerical policy is
left enabled or shipped. H6's energy-metric rule and first-step redo are not
implemented because this task explicitly makes further implementation
conditional on the measured bound. This result does not disprove H6.

H2 already exists in the current engine: strict bands are occupied bands plus
the first empty band; higher buffers use `max(tolerance, 1e-5)`. The lowest
empty band therefore retains the strict tolerance for the sampled gap.
The historical difference table in note 1 predates that implementation.
Buffer-origin H fractions including the certificate are approximately 24.92%
(si8), 34.39% (k444s), 10.98% (skewed) and 28.65% (elongated). These establish
where work goes, but do not validate another threshold policy.

## H7: rejected smaller band counts

One extra band instead of four reduces k444s H vectors from 11181 to 8310
without changing its 10 SCF iterations; energy differs by 6.57e-14 Ha. But si8
increases to 2444 vectors at 15 iterations, and the skewed case fails to
converge in 100 iterations (3260 vectors). Reject as a general default.

Two extra bands give k444s 9068 vectors at 10 iterations, but si8 needs 17 SCF
iterations and 3249 vectors. The skewed case needs nine iterations and reaches
an energy 29.3839772182878 Ha above the default result despite passing its
certificate. That different fixed point fails the same-result rule; its cause
has not been established. Reject as a general default. `nbands` already exposes
this experiment, so no new option was needed.

## Instrumentation and attribution limits

`P14A_WORK_ACCOUNTING=1 cargo build --release -p dft-cli` enables the source-only
measurement collector. It is off in an ordinary build. It allocates root arrays
once per solve and reserves a bounded 16384-record collector before iterations;
no timers, logging, file I/O, serialization, or per-iteration allocations were
added. The CLI drains records after the SCF and adds `davidson_work` only to
measurement-build results. Production scientific result schemas are unchanged.

Records contain actual block calls and vector columns, restart counts,
first-below/last-above threshold iteration for every root, and root-origin
counts. Successful reduced directions can mix several roots during projection;
"buffer-origin" means the seed correction came from a buffer root, not that the
resulting vector has only buffer character. Rank loss and random fallback are
unattributed. Four elongated-case vector applications are unattributed, so its
reported buffer fraction is a lower bound by at most 4/171. The first empty
band belongs to the strict prefix and is not included in the buffer fraction.

The table's early/late counts use the last above-threshold iteration, so a root
that first converges early but later loses convergence is not falsely counted
as early. Collector order is solve-completion order; the per-SCF grouping is
validated for these one-thread measurements only. No ordered trace guarantee
is made for concurrent k solves. Instrumentation is not an exact checkpoint API.

Davidson already selects corrections only for unconverged roots. Its 14 si8
restarts reuse cached H products; there is no measured restart-related H
reapplication to remove in these successful runs. Rotations/residual formation
still span all roots and are plausible future kernel work, but no optimization
is adopted without profiling and equivalence evidence. Existing phase timings
are SCF-region wall times, not a decomposition of internal Davidson GEMMs.

## QE calibration and provenance

[The relayed file](../../../tests/reference/perf/p14a/qe_iteration_data_relayed.md)
was read directly from main as instructed; git metadata was not writable, so
main was not merged. The relay gives dev08 paths, QE version, complex-Gamma
input mode, 20 bands and one MPI process. It is not an independently rerun
reference or a substitute for complete input/UPF provenance.

The si8 table aligns all 14 QE SCF iterations with both relayed and local
BAANDI steps/tolerances, and retains BAANDI's fifteenth iteration. QE reports
83 Davidson steps and 98 H block calls. Relayed BAANDI reports 290 steps and
2070 H vector applications; local results above differ, including when the
certificate is excluded. No same-build/node identity is asserted. QE block
widths are unavailable, so its vector count cannot be reconstructed by simply
multiplying calls by 20. QE ethr is a Ritz-value-change threshold in Ry, whereas
our threshold is a residual norm in Ha: unit conversion alone does not equate
the stopping criteria. The relay also acknowledges a total-energy discrepancy.

The k444s relay is ambiguous about five versus six iterations and does not give
one-to-one threshold/error pairings. Its aggregate 1010 H block calls and six
reported average-step values are retained as relayed; missing rows are not
invented. Stored `k222_c/qe.out` and `k222s_c/qe.out` supply the complete older
iteration records and SHA-256 hashes in baseline.json, including both shifted
case first-iteration solves. These are identified separately from dev08 si8.

## Reproduction and verification

```sh
PATH=/opt/homebrew/opt/rustup/bin:$PATH cargo build --release -p dft-cli --offline
python3 -I tests/p14a_work_accounting.py --output /tmp/p14a-default.json
python3 -I tests/p14a_work_accounting.py --extra-bands 1 --output /tmp/p14a-extra1.json
python3 -I tests/p14a_work_accounting.py --extra-bands 2 --output /tmp/p14a-extra2.json
P14A_WORK_ACCOUNTING=1 cargo build --release -p dft-cli --offline
python3 -I tests/p14a_work_accounting.py --output /tmp/p14a-accounting.json
P14A_WORK_ACCOUNTING=1 cargo test --release -p dft-solvers work_accounting_counts -- --ignored
python3 -I tests/p14a_summarize.py
```

For the rejected sweep, on a clean session.rs apply the archived patch with
`patch -p1 < tests/reference/perf/p14a/tolerance-factor.patch`, then build/run
once for each factor with `P14A_TOLERANCE_FACTOR=3 cargo build --release -p
dft-cli --offline` (replace 3 by 1, 10, 30, 100). Give the accounting runner a
fresh output and `--label factor=3`. Reverse the patch with `patch -R -p1` using
the same file, unset measurement variables and rebuild before normal use.

Final checks:

- `scripts/check.sh`: exit 0, build/test/fmt/clippy pass. A preceding invocation
  caught a constant assertion lint in the diagnostic-only test; fixed and the
  complete command rerun. The independent operator wrapper test passes and
  agrees on blocks, vectors, iterations and nonzero restarts.
- `scripts/bench_scf.py`: historical baseline comparison exit 0 (Gamma delta
  1.63e-13 Ha, k222 6.04e-14 Ha; k444s absent from that baseline). Final comparison
  against the stored pre-change run exits 0: Gamma/k222/k444s energies exactly
  equal and density files bit-identical.
- `scripts/bench_suite.py --robustness`: historical comparison exits 5 because
  synthetic inputs/UPFs used Xx then and Si now. Real-Si comparisons show no
  regression. The suite against the pre-change report with identical current
  inputs exits 0 with `--require-identical`. Existing nonconvergence of
  si-elongated, syn-s3 and syn-z14 remains, not hidden by rewriting a baseline.
- All four accounting systems have identical energy, density bytes, SCF count
  and H vector totals at threads 1 and 4 within the pure-Rust backend; see
  thread-comparison.json. MKL was not run locally.
- README generation makes no changes and all three README tests pass. No new
  accepted performance result exists to add to its generated validation block.

The canonical registry outside this worktree records this delegated run; the
local registry snapshot predates it. Dashboard/session files and their upload
remain orchestrator-owned and outside this ticket's scope. No active recurring
report timer or dashboard upload is claimed. All changes remain uncommitted.
