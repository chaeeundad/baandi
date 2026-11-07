# P39A: pure-Rust mirror-slab pair FFT

**DONE under the architect-approved 2026-11-08 gate.** All strict fixtures
pass delivered accuracy against their same-input complex-tight reference; all
six standard fixtures pass the unchanged relative criterion. The original
inter-trajectory gate failure remains reported below as superseded evidence.
Mac measurements are complete; dev08 performance acceptance remains separate.

Reviewable implementation on `task/P39A-PAIR-FUSED-PLAN`, starting HEAD
`f66e3966eceb5cedaca6411910f6738e3c8f4866`, uncommitted. Raw evidence,
executables, commands, inputs and scientific output are in `P39A-evidence/`,
outside the single `target/`. Mac observations only; no dev08 execution,
new QE run, dependency, tolerance or default change. Artifacts are untouched
as explicitly required by the ticket. No dashboard upload or recurring session
timer was established in this implementer session.

## Design and ownership

Implements T1 of [pair-fft-efficiency.md](pair-fft-efficiency.md), section 5.
The premise is PFFT-1/2/3: intermediate traffic and normalization overhead
in the pruned pair path, with a small slab working set. This does not revive
A27's existing plane plan, A29's vendor pruning threshold, A34's FFT threads,
A35's scoped threads or PFFT-4's SIMD/index-table alternatives. Those mechanisms
remain unchanged. No QE source was read. The design's connection to the
own-words [FFT study note](../research/qe-study/qe-hpsi-and-fft.md) is inherited;
implementation uses rustfft's public API and this engine's conventions.

`PairFusedPlan` owns plane-major sticks, a compact middle-axis batch, one
last-axis slab and maximum-sized rustfft scratch. Inverse transforms are
unnormalized; local application mirrors the inverse order and divides by the
full grid size once on gather. This deliberately changes rounding. Density
uses the same inverse half, adding weighted squared real/imaginary components
without a forward transform. Single-band density retains weight-times-norm
ordering. Duplicate indices retain last-write-wins scatter semantics; empty
support, invalid dimensions, invalid slices and invalid density parameters
have explicit tests. All scratch is constructed once, with no internal threads.

The design's axis-choice wording needs a precedence rule when z has the fewest
sticks: z remains last for contiguous potential lines; choose the fewest
occupied sticks **among x and y**, breaking ties toward x. If y is first,
map every slab row back to the original C-order grid. Tests include a y-first
elongated shape. There is no occupancy threshold or machine tuning.

`dft-physics` has no Cargo `mkl` feature: the existing unified feature belongs
to `dft-numerics`. Consequently the cfg-gated `RealPairPlan` alias lives in
the new numerics module. It selects `PairFusedPlan` only without MKL and the
unchanged `fused::FusedPlan` with MKL. Physics changes its local workspace type
and constructor to this alias. An initial physics-local cfg produced compiler
warnings; it was corrected before final verification. No manifest changes
or new feature propagation are needed. The global complex `FusedPlan` alias,
`fused.rs` and `fused_pair.rs` are unchanged byte for byte.

The current engine computes density from completed complex orbitals and has
no `PairDensityWorkspace` caller. This ticket implements and tests the new
pair-density backend without changing engine density dispatch. The new library
plan is usable by the existing generic workspace, including its fixed four-pair
pieces and ordered reduction. No other production plumbing was required.

## Fixed validation and limitations

The release kernel matrix uses the nine design shapes plus a y-first elongated
shape, including shear 0.35, odd lengths, vacuum and empty planes. Relative error
means L2/Frobenius norm of the difference divided by the reference norm; no
entrywise-relative assertion is made at cancellation zeros. Local and pair
(and single-band) density are compared against unchanged `OrderedFusedPlan`.
Repeated calls and T=1,2,4,8,12 compare every output bit through local and density
workspaces, including a nine-band block with an odd tail.

A test-only System allocator adapter counts allocations from the first apply
following construction through repeated local and density calls. It forwards
pointer/layout contracts unchanged, following the existing P33A test adapter;
production has no unsafe or instrumentation. Run this ignored test with one
test thread to exclude unrelated test-runner activity. Six shapes report zero
allocations. A separate unit test checks persistent numerical-buffer pointers,
lengths and capacities. These claims concern the pair FFT, not the solver/EVD.

The existing snapshot harness supplies seven systems (the six R5 fixtures
already include primitive Si and Si8, plus Si16), strict and standard, at all
five budgets. Each saved file contains energy/density/eigenvalue bits, SCF work
and complete solver results. All 14 groups are byte-identical across budgets.
Strict before/after energy deltas pass 1e-12 Ha, maximum
7.815970093361102e-13 Ha on Si16. Before/after trajectories are not identical.
In particular Si8 strict SCF changes 15 to 14 and Si16 strict 13 to 15.
The architect accepts changed trajectory counts with the delivered-accuracy
gate; no unchanged-count claim is made. No arithmetic
variant or tolerance was tuned to recover a preferred iteration count.

The historical material strict test still fails its real-versus-complex
primitive count condition (10 versus 12). The standard synthetic slab still
fails its 5e-9 Ha absolute cap: real error 6.009213793234380e-9 Ha, complex error
4.598657676879725e-9 Ha. Its relative cap 1.379597303063917e-8 Ha passes.
The architect disposition at the end of `real-gamma-r5-results.md` keeps this
absolute failure as an observation; it is not represented as a passing gate.
Both standard material fixtures pass the unchanged delivered-accuracy checks.
The architect explicitly excludes these pre-existing observations from this
ticket on 2026-11-08. They are not new failures caused by the pair plan.

## Reproduction

Use the pinned toolchain, `/opt/homebrew/opt/rustup/bin` on PATH. External UPF:
`../pseudo/Si.pz-vbc.UPF` (never copied/committed), SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Binary and protected-source hashes are in `P39A-evidence/metadata.json`.

```sh
cargo build --release -p dft-cli
# Preserve before-cli before edits, after-cli after rebuilding, in P39A-evidence.
python3 -I P39A-evidence/identity.py before
python3 -I P39A-evidence/identity.py after
DFT_SI_UPF=/absolute/Si.pz-vbc.UPF P33A_SNAPSHOTS=/absolute/snapshot-directory \
  cargo test --release -p dft-engine --test gamma_workspace_snapshot -- --ignored --nocapture
python3 -I P39A-evidence/compare_snapshots.py
cargo test --release -p dft-numerics --test pair_fused_plan -- --include-ignored --nocapture --test-threads=1
cargo test --release -p dft-numerics --test real_gamma_pair pair_microbenchmark -- --ignored --nocapture --test-threads=1
python3 -I P39A-evidence/measure.py
scripts/check.sh
python3 -I scripts/make_readme_figures.py
python3 -I scripts/tests/test_readme_figures.py
```

Identity driver runs `bench_scf --require-identical` on gamma/k222s/k222/k444s/si8
at every requested budget, with same-HEAD before/after executables. All 25
comparisons pass, energy delta 0 Ha and identical density hashes. Correctness
runs may overlap builds; their wall times are not performance evidence.

README regeneration and nine freshness tests pass, with no generated changes.
No accepted reference dataset or speed claim was invented for README. The
snapshot parser's first attempt found an ambiguous suffix resembling a history
count; the corrected parser selects the earliest valid history boundary and
retains raw bytes. Snapshot determinism and energies are independent of that
work-count decoding; CLI counts provide independent material confirmation.

## dev08 request

After review, run the release microbenchmark above on the pinned idle node,
using its normal native flags, and retain all 30 alternating sample vectors,
binary/source/UPF hashes and load. It prints microseconds **per band**, with
library pair, new mirror pair and the preserved unpaired reference rows.
Then alternate parent and candidate for Si8/Si16 strict/standard at T=1/8
(and the full thread matrix for determinism). Report main and certificate
SCF/Davidson/H work, accuracy, density identity and trajectories separately
from time. Profile at least 30 working samples in local FFT; check the L2
hypothesis on the elongated grid. Do not apply Mac ratios to QE as a measured
result. MKL and SIMD kernels remain separate tickets; no remote action occurred.

## Kernel measurements (Mac observations)

30 retained alternating samples after two warmup rounds, each sample 20 pair
applies (40 bands). The unpaired comparison rows remain in the benchmark.
No task-owned build or other measurement overlaps this run; host load at its
end was 8.13/11.75/8.41, so absolute times are not idle-node acceptance.

| Grid | Library pair (us/band) | Mirror pair (us/band) | Mirror/library | Library spread | Mirror spread |
|---|---:|---:|---:|---:|---:|
| 36x36x36 | 223.668225 | 180.353125 | 0.806342 | 32.6% | 33.4% |
| 48x48x48 | 504.942700 | 419.628137 | 0.831041 | 30.5% | 20.3% |
| 72x36x36 | 389.630738 | 334.583325 | 0.858719 | 56.2% | 24.0% |

Spread is (max minus min)/median, including every retained sample. Observed
ratios 0.806–0.859 are favorable on this host; they do not establish the dev08
cache hypothesis or any QE speed ratio. Raw ordered sample arrays are retained
in `microbenchmark.log` and `micro-samples.json`.

## SCF snapshot work and energy

Main-history counts; H counts vectors, not FFTs or block calls. Certificate
state is retained in the binary snapshots. Every row is deterministic at all
five budgets; energy delta is candidate minus parent.

| Fixture/profile | Delta E (Ha) | SCF before/after | Davidson before/after | H before/after |
|---|---:|---:|---:|---:|
| si16-strict | 7.81597009336110204e-13 | 13/15 | 237/290 | 3043/3607 |
| si16-standard | 1.03567288078920683e-09 | 16/14 | 84/66 | 1526/1383 |
| si2-strict | 1.06581410364015028e-14 | 12/12 | 54/55 | 345/346 |
| si2-standard | 2.66453525910037570e-15 | 9/9 | 26/26 | 176/176 |
| si8-strict | 2.27373675443232059e-13 | 15/14 | 253/220 | 2001/1844 |
| si8-standard | 8.11962053148818086e-10 | 17/15 | 89/78 | 877/793 |
| synthetic-0-strict | -2.22044604925031308e-16 | 10/10 | 42/42 | 122/123 |
| synthetic-0-standard | 1.15840226300178983e-10 | 6/6 | 17/15 | 67/64 |
| synthetic-1-strict | 3.33066907387546962e-16 | 11/11 | 50/50 | 161/161 |
| synthetic-1-standard | -7.77156117237609578e-16 | 8/8 | 24/24 | 94/94 |
| synthetic-2-strict | 4.44089209850062616e-16 | 11/11 | 34/34 | 118/118 |
| synthetic-2-standard | 0.00000000000000000e+00 | 6/6 | 14/14 | 65/65 |
| synthetic-3-strict | 6.66133814775093924e-16 | 10/10 | 59/59 | 205/205 |
| synthetic-3-standard | -1.77635683940025046e-15 | 17/17 | 55/55 | 212/212 |

Kernel maximum relative local/density errors: 2.70340812071420126e-16 / 4.82130584299562711e-16 (dimensionless); both below 1e-13.

## Alternating end-to-end Mac observations

One discarded warmup pair and five retained alternating pairs per condition;
no task-owned builds or other benchmarks overlapped. Wall includes CLI
setup/output. Every repeat within a binary retains identical energy, density
and work. These are changed-trajectory timings, not equal-work FFT speedups.
Five pairs on a shared host do not resolve a small-percent performance claim.

| Case/profile/T | Before median (s) | After median (s) | After/before | SCF before/after | Davidson before/after | H before/after | Same trajectory |
|---|---:|---:|---:|---:|---:|---:|---|
| si8-strict-T1 | 1.290882 | 1.170578 | 0.906805 | 15/14 | 253/220 | 2001/1844 | False |
| si8-strict-T8 | 0.861383 | 0.795594 | 0.923624 | 15/14 | 253/220 | 2001/1844 | False |
| si8-standard-T1 | 0.535252 | 0.453623 | 0.847493 | 17/15 | 89/78 | 877/793 | False |
| si8-standard-T8 | 0.378819 | 0.325930 | 0.860382 | 17/15 | 89/78 | 877/793 | False |
| si16-strict-T1 | 5.831980 | 4.775239 | 0.818802 | 14/14 | 264/201 | 3254/3151 | False |
| si16-strict-T8 | 3.160054 | 2.610292 | 0.826028 | 14/14 | 264/201 | 3254/3151 | False |
| si16-standard-T1 | 2.113658 | 1.669617 | 0.789918 | 16/13 | 88/68 | 1521/1310 | False |
| si16-standard-T8 | 1.263471 | 1.020844 | 0.807968 | 16/13 | 88/68 | 1521/1310 | False |

| Condition | Before spread | After spread | 1-minute load range |
|---|---:|---:|---:|
| si8-strict-T1 | 0.4% | 0.6% | 7.33–7.80 |
| si8-strict-T8 | 0.6% | 1.2% | 6.89–7.33 |
| si8-standard-T1 | 2.2% | 0.9% | 6.74–6.89 |
| si8-standard-T8 | 0.9% | 1.7% | 6.74–6.76 |
| si16-strict-T1 | 3.7% | 0.4% | 5.60–6.49 |
| si16-strict-T8 | 1.1% | 1.6% | 5.26–6.07 |
| si16-standard-T1 | 0.7% | 0.8% | 5.58–5.91 |
| si16-standard-T8 | 0.8% | 2.0% | 5.45–5.90 |

The exact benchmark-suite Si16 input differs from the snapshot fixture in
seed and atom ordering. It is retained verbatim except for the explicit
thread and real/standard options. Its strict before/after delta is
-1.438849039914203e-11 Ha at both T=1 and T=8: this failed the original
inter-trajectory gate. That gate was superseded by the architect on 2026-11-08;
the delivered-accuracy table below is now the acceptance criterion. FFT norm errors
are near roundoff, while SCF histories change; no unit/normalization error is
indicated by the kernel oracle. The changed finite-tolerance trajectory is a
plausible explanation, not an independently proved diagnosis. Tolerances,
seeds, solver and defaults were not adjusted to suppress the failed row.

All eight before/after history comparisons differ (timing fields excluded).
Raw results include certificate counts separately from the main-history totals
shown here. The faster timings are changed-work Mac observations, not isolated FFT
speedups or dev08 performance acceptance. The revised accuracy gate passes.

## Completed commands and disposition

| Command (logs under P39A-evidence) | Exit | Result |
|---|---:|---|
| `cargo build --release -p dft-cli` before/after/final | 0 each | Saved before/after binaries |
| `cargo test --release -p dft-numerics --test pair_fused_plan -- --include-ignored --nocapture --test-threads=1` | 0 | Three tests, ten shapes, six zero-allocation probes |
| `cargo test --release -p dft-engine --test gamma_workspace_snapshot -- --ignored --nocapture` before/after, with UPF/snapshot variables | 0 each | 70 snapshots per build |
| `python3 -I P39A-evidence/compare_snapshots.py` corrected parser | 0 | Strict snapshot energy and thread checks pass |
| `python3 -I P39A-evidence/identity.py before` and `after` (final binary rerun) | 0 each | Five bench_scf commands per invocation, all exit 0 |
| `python3 -I P39A-evidence/verify_off.py` | 0 | Exact energy bits, density hashes and main/certificate work in all 25 comparisons |
| `cargo test --release -p dft-engine --test gamma_real_scf standard_real_complex_and_threads -- --ignored --nocapture --test-threads=1` | 101 | Historical slab absolute standard cap fails; relative caps pass |
| `DFT_SI_UPF=... cargo test --release -p dft-engine --test gamma_real_materials -- --include-ignored --nocapture --test-threads=1` | 101 | Standard material test passes; historical primitive strict count test fails |
| `scripts/check.sh` (whole workspace, no crate argument; repeated 2026-11-08) | 0 | PASS build/test/fmt/clippy |
| `cargo test --release -p dft-numerics --test real_gamma_pair pair_microbenchmark -- --ignored --nocapture --test-threads=1` | 0 | All original rows preserved; mirror rows and raw samples added |
| `python3 -I P39A-evidence/measure.py` | 0 | Eight conditions, five retained alternating pairs each |
| `python3 -I P39A-evidence/summarize.py` (original gate, 2026-11-07) | 1 | Historical inter-trajectory gate failure; superseded by the architect decision |
| `python3 -I P39A-evidence/tight_references.py` | 0 | Recomputed seven snapshot and two exact benchmark tight references, same target tree |
| `python3 -I P39A-evidence/summarize.py` (revised gate, 2026-11-08) | 0 | All 11 strict rows and six standard relative rows pass; idempotent report section |
| `python3 -I P39A-evidence/si16_threads.py` | 0 | Exact benchmark Si16 strict/standard at T=2/4/12 matches T=1/8 energy/density/work |
| `python3 -I scripts/make_readme_figures.py` | 0 | Generated output unchanged |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | Nine tests pass |
| `git diff --check` | 0 | Clean whitespace |

**Final status: DONE under the revised gate.** Only the architect-authorized
acceptance criterion changed; production arithmetic, solver tolerances,
component caps and measured real inputs were not changed. Keep this patch
uncommitted for review, with the raw evidence intact. Dev08 performance
acceptance remains pending; Mac observations are not a remote speed claim.

The reference-only harness is generated from the retained snapshot test into
`P39A-evidence/tight_reference_snapshot.rs`, compiled through a temporary
integration-test file, then removed by `tight_references.py`. This is test
plumbing only; no extra target directory or persistent engine file was added.
The two CLI tight inputs and full results are retained. The recomputed
Si8/Si16 reference energies differ from the relayed dev08 values by
1.77635683940025046e-14 / -2.62900812231237069e-13 Ha, respectively; the
local same-input references, not the remote rounding, determine the gate.

The completed Mac timing table above requires no new runs because production
code and measured inputs did not change. The summarizer preserves that table
and replaces only its generated delivered-accuracy section, so repeated
invocations do not duplicate tables. Original-gate script/logs remain archived
as historical evidence. README generation and its nine freshness tests were
repeated successfully; generated output remains unchanged.

## Architect-approved delivered accuracy (2026-11-08)

The architect replaced the inter-trajectory 1e-12 Ha gate with
`abs(E_new - E_tight) <= max(2*abs(E_previous_real - E_tight), 2e-11 Ha)`.
Every reference below was recomputed locally with gamma_real off, energy
tolerance 1e-12 Ha and density residual 1e-10; eigen tolerance and every
other fixture setting remain unchanged. The six R5 fixtures include primitive
Si and Si8; the seventh snapshot fixture is Si16. Benchmark-suite Si16 has
a different seed/atom order and is checked separately, alongside benchmark Si8.
The stored dev08 values are provenance cross-checks, not substituted inputs.

All printed numbers use 17 digits after the decimal in scientific notation
(round-trip precision for f64); exact binary energies are also recorded as
hexadecimal strings in `delivered-accuracy.json`. No abbreviated energy delta
is used in these gate tables.

| Fixture | Previous real E (Ha) | New real E (Ha) | Tight E (Ha) |
|---|---:|---:|---:|
| si16 | -6.31292444216146009e+01 | -6.31292444216138193e+01 | -6.31292444216169102e+01 |
| si2 | -7.29553046821089346e+00 | -7.29553046821088280e+00 | -7.29553046821087925e+00 |
| si8 | -3.13218036010316254e+01 | -3.13218036010313980e+01 | -3.13218036010317533e+01 |
| synthetic-0 | -6.23797726399250574e-01 | -6.23797726399250796e-01 | -6.23797726399251906e-01 |
| synthetic-1 | -8.47124717547730555e-01 | -8.47124717547730222e-01 | -8.47124717547730666e-01 |
| synthetic-2 | -3.26376645234357010e+00 | -3.26376645234356966e+00 | -3.26376645234360829e+00 |
| synthetic-3 | -1.44295754190504733e+00 | -1.44295754190504666e+00 | -1.44295754190504910e+00 |
| benchmark-si8-strict-T1 | -3.13218036010316254e+01 | -3.13218036010313980e+01 | -3.13218036010317533e+01 |
| benchmark-si8-strict-T8 | -3.13218036010316254e+01 | -3.13218036010313980e+01 | -3.13218036010317533e+01 |
| benchmark-si16-strict-T1 | -6.31292444216030262e+01 | -6.31292444216174147e+01 | -6.31292444216179334e+01 |
| benchmark-si16-strict-T8 | -6.31292444216030262e+01 | -6.31292444216174147e+01 | -6.31292444216179334e+01 |

| Fixture | New minus previous (Ha) | New minus tight (Ha) | Previous absolute error (Ha) | New absolute error (Ha) | Allowed error (Ha) | Result |
|---|---:|---:|---:|---:|---:|---|
| si16 | 7.81597009336110204e-13 | 3.09086090055643581e-12 | 2.30926389122032560e-12 | 3.09086090055643581e-12 | 1.99999999999999988e-11 | PASS |
| si2 | 1.06581410364015028e-14 | -3.55271367880050093e-15 | 1.42108547152020037e-14 | 3.55271367880050093e-15 | 1.99999999999999988e-11 | PASS |
| si8 | 2.27373675443232059e-13 | 3.55271367880050093e-13 | 1.27897692436818033e-13 | 3.55271367880050093e-13 | 1.99999999999999988e-11 | PASS |
| synthetic-0 | -2.22044604925031308e-16 | 1.11022302462515654e-15 | 1.33226762955018785e-15 | 1.11022302462515654e-15 | 1.99999999999999988e-11 | PASS |
| synthetic-1 | 3.33066907387546962e-16 | 4.44089209850062616e-16 | 1.11022302462515654e-16 | 4.44089209850062616e-16 | 1.99999999999999988e-11 | PASS |
| synthetic-2 | 4.44089209850062616e-16 | 3.86357612569554476e-14 | 3.81916720471053850e-14 | 3.86357612569554476e-14 | 1.99999999999999988e-11 | PASS |
| synthetic-3 | 6.66133814775093924e-16 | 2.44249065417534439e-15 | 1.77635683940025046e-15 | 2.44249065417534439e-15 | 1.99999999999999988e-11 | PASS |
| benchmark-si8-strict-T1 | 2.27373675443232059e-13 | 3.55271367880050093e-13 | 1.27897692436818033e-13 | 3.55271367880050093e-13 | 1.99999999999999988e-11 | PASS |
| benchmark-si8-strict-T8 | 2.27373675443232059e-13 | 3.55271367880050093e-13 | 1.27897692436818033e-13 | 3.55271367880050093e-13 | 1.99999999999999988e-11 | PASS |
| benchmark-si16-strict-T1 | -1.43884903991420288e-11 | 5.18696197104873136e-13 | 1.49071865962469019e-11 | 5.18696197104873136e-13 | 2.98143731924938038e-11 | PASS |
| benchmark-si16-strict-T8 | -1.43884903991420288e-11 | 5.18696197104873136e-13 | 1.49071865962469019e-11 | 5.18696197104873136e-13 | 2.98143731924938038e-11 | PASS |

| Standard fixture | New absolute error (Ha) | Unchanged relative cap (Ha) | Relative result | Historical absolute cap exceeded |
|---|---:|---:|---|---|
| kind=0 | 2.18075779656601298e-10 | 5.18349363431980237e-10 | PASS | False |
| kind=1 | 2.38222552795264164e-10 | 7.14669989854144205e-10 | PASS | False |
| kind=2 | 2.87425638845206777e-10 | 1.20181331553226300e-09 | PASS | False |
| kind=3 | 6.00921379323438032e-09 | 1.37959730306391748e-08 | PASS | True |
| gamma | 2.33610464306366339e-10 | 3.24204663115779113e-10 | PASS | False |
| si8 | 2.07768380278139375e-09 | 3.00267544162124977e-09 | PASS | False |

The standard slab absolute exceedance and primitive strict 10-versus-12 SCF
count remain visible as pre-existing observations, excluded from this ticket
by the architect decision. Component tolerances, allocation, thread identity
and complex-path identity requirements are unchanged. No production source
or input tolerance of either measured real trajectory changed in this follow-up.
