# P42A band-boundary policy study

**Decision: Phase A is complete; do not add an automatic option from this
study.** A converged gap at the last requested band is neither necessary nor
sufficient for stable standard-profile counts in the measured cases. The frozen
converged-spectrum diagnostic violates held-out work and accuracy bounds,
does not reduce median work or spread, and loses a previously convergent seed.
No production code, solver, tolerance, profile or integer-four default changed.
All results are local pure-Rust observations, not dev08 speed acceptance.

## Verified final coverage

All 25 cases / 44 applicable case-path groups / 676 extra-band conditions were
attempted: **8,112 starts**, comprising 6,760 standard perturbations and 1,352
strict reference attempts. Impossible band counts remain setup rejections;
they are not executed numerical SCFs or successful zero-work results.

| Outcome | Standard starts | Strict reference attempts |
|---|---:|---:|
| Converged | 4044 (uncertified) | 810 (certified) |
| Basis/band setup rejection | 1540 | 308 |
| Iteration-budget exhaustion | 1042 | 212 |
| Density-growth divergence | 134 | 22 |
| Total | 6760 | 1352 |

Both strict levels converge in 404 conditions; **396** also meet the unchanged
1e-12 Ha energy-plateau test. Eight converged reference pairs fail that plateau
and are not accuracy oracles. All 404 available tighter-reference densities
received a wide probe. A usable N+1 boundary is available in 356 conditions,
348 also have a valid reference plateau, and 335 additionally have all ten
standard starts converged for the primary wide-spectrum correlations. Missing
boundaries (basis ceiling or prefix disagreement) remain unknown. The narrow
cross-count audit flags 58 disputed spectra; raw values are retained.

The full evidence audit exits 0: every planned condition has 12 actual terminal
records, no reservation remains, all scientific execution processes exit 0,
work sums and density provenance match, and failures never receive a fabricated
converged energy. Numerical failures above are study outcomes, not hidden
collector successes. Measurement-binary and preregistration hashes are constant
across all conditions; fixture and UPF hashes are constant within each case.
See [complete case/path tables](../../P42A-evidence/tables.md),
[all condition rows](../../P42A-evidence/summary.csv), and
[audit log](../../P42A-evidence/audit-final.log).

## Preregistered protocol (before P42A measurements, 2026-11-08)

Measurement only until Phase A establishes a usable relation. Default integer
four, production SCF, solver, tolerances and profile rules remain unchanged.
A19/P38A are reopened because P41A and the orchestrator's band-cut-stability
scan concern perturbation distributions and the LAST requested band, rather
than first-empty multiplet completion or a single strict trajectory. A08/A13
capacity tuning is not repeated. No QE source is read.

Training: si8-gamma, si16-gamma, si-gamma, si-skew, si-low,
syn-cubic-degenerate, material-si-20 (original PseudoDojo Si input).
Held out: si-k222, si-k222s, si-elongated, si-slab, si-k444s,
material-mgo-40, and every other benchmarks/cases/syn-*.toml fixture.
Held-out descriptive scans must not select parameters. Freeze a candidate on
training evidence before evaluating its held-out recommendation gates.

Extras: all integers 1 through 16, except si16-gamma and syn-100-bands use
1,2,4,6,8,10,12,14,16. Preserve basis, cutoff, mesh, XC, material floor/mixing
and synthetic projector parameters. An impossible band count is a recorded
setup failure, not grounds to enlarge the basis. Complex and forced real for
single exact Gamma; multi-k only complex. T=1 for count measurements, including
materials (P29A used T=4 for materials; arithmetic determinism is a separate gate).
Use the engine default orbital seed as P29A, not P38A's seed 42. Preserve
synthetic fixture FFT grid and iteration budget when explicitly specified;
otherwise use P29A's engine defaults. Ignore explicit fixture nbands in favor
of occupied+extra. No external UPF is copied into the worktree.

Each condition: P29A unperturbed strict references 0 and 1, then ten standard
starts using identical mean-centered LCG seeds 1–10 and relative amplitude
1e-10. Reference tolerances (density, energy Ha, eigen Ha) are
(1e-9,1e-11,1e-10) and (1e-10,1e-12,1e-11). Energy accuracy is available only
when both certify and agree within 1e-12 Ha. Compare every extra with the
extra=4 tight reference as well as its own reference to expose different fixed
points. Export occupied eigenvalues and densities for eigenvalue maximum error
(Ha) and density RMS error (Bohr^-3); do not call standard results certified.
Report SCF min/median/max among successes, all failures and their attempt work,
median H vectors including certificates, separate ordinary/certificate work,
and strict reference iteration counts. Do not substitute censored lengths for
converged SCF counts. Phase B additionally requires ten ordinary strict starts
per condition, not treating two tight references as a strict distribution.

For each converged tight reference, run an independent frozen-density spectrum
probe with one band beyond its own requested boundary N. Probe residual target
is 1e-5 Ha on every root (the existing buffer floor); export actual residuals.
Also retain the reference's original eigenvalues from occupied-2 through N.
Let absolute gap g=e[N+1]-e[N] in Ha, relative gap r=g/(e[N]-e[max(1,N-7)]),
with one-based indices. A zero denominator is unknown, never infinity. At
multiple k points report every gap and use the minimum over k for correlation.
A gap <=2e-5 Ha is unresolved at the probe tolerance, not a certified degeneracy.
Insufficient basis or failed probe is unknown. No boundary inference from a
spectrum that ends at N. Frozen probes are diagnostic work, excluded from SCF
work but explicitly recorded; any deployed predictor must charge its own work.

Descriptive relation: per-case Spearman rank correlations of absolute/relative
gap with SCF spread and median H (ties use average ranks), plus pooled results
labelled confounded by case. Report undefined correlations for constant arrays.
Also report number of distinct iteration counts, distinct per-iteration integer
work sequences, and exact scalar-history sequences; equal final iteration count
alone does not establish equal trajectories. Compare high/low gap groups using
fixed thresholds 2e-5 Ha and relative 0.01,0.05,0.1 descriptively, without tuning
on held-out data. Failed references/SCFs remain explicit missing observations.

## Frozen recommendation rule

Recommend an opt-in rule only if, on held-out cases and perturbations:
- Median SCF falls or stays equal on at least 80 percent of cases.
- Median H work rises by no more than 10 percent on any case and falls on
  the median over cases (equally weighted default-relative ratios).
- SCF spread (max-min) increases by no more than 20 percent on any case
  and shrinks on the median. A default zero spread must remain zero.
- No default-convergent case/start becomes a failure.
- Delivered energy, occupied-eigenvalue and density errors are each within
  1.5 times default error against the same tight reference; strict certifies.
  Zero default error requires zero candidate error, no post-hoc floor.
All conditions apply separately to strict/standard and each applicable path.
Missing accuracy or coverage cannot count as a pass. Never relax these bounds.
No default change is authorized even if the rule passes.

A Phase B rule must use only setup kinetic levels or a first loose solve with
buffer roots, possibly a fixed candidate ladder; deterministic across threads,
one common count across k (all-k boundary criterion), fallback four if unknown,
chosen count in results and exact-restart metadata. A converged-spectrum oracle
is descriptive, not such a rule. No Phase B production option is yet justified.

## Reproduction and gates

Raw evidence and collector scripts live in P42A-evidence/, never target/.
One existing local target directory only. Heavy integration tests are ignored.
Full workspace command is `scripts/check.sh` (literal `all` names a nonexistent
crate). Identity: bench_scf --require-identical, gamma/k222s/k222/k444s/si8,
T=1,2,4,8,12, strict/standard, real off/on; unsupported forced-real multi-k
conditions are recorded rejections, not passing comparisons. Mac counts are
observations; dev08 acceptance requires the orchestrator, no remote login here.
No artifacts/* edits are authorized by this ticket.

### Spectral reliability observation during collection

The primitive real-Gamma N+1 probe at extra=6 omits a member of a higher
multiplet despite residual certification (a wider converged spectrum exposes
it). Residual norms alone cannot certify lowest-root completeness. All raw
gaps remain available, but a row is flagged `spectrum_disputed` and excluded
from the primary correlation if any wider same-path probe at the same fixed
point (tight energies within 1e-9 Ha) disagrees on its common prefix by more
than 1e-4 Ha. This is a diagnostic validity filter, not a tuned gap threshold
or relaxed recommendation gate. Both the original preregistration and raw
spectra remain preserved. Undisputed spectra are not an independent dense
completeness proof; recommendations cannot rely on an unresolved root count.

Collection scheduling: training and held-out descriptive scans run in separate
T=1 processes after the initial primitive training smoke. No candidate has been
selected and held-out values will not tune one. Wall times record loaded-host
execution only; no timing claim or acceptance comparison uses them. Builds and
verification may overlap these count-only observations; each process records
load, source/binary hashes, command and wall duration. This scheduling does not
change the scientific split, perturbations or decision thresholds.

### Frozen descriptive oracle (not a production candidate)

Before inspecting held-out spectrum/trajectory summaries, evaluate the fixed
hypothesis: keep default four if its boundary is resolved and r>=0.05 at every
k; otherwise take the smallest scanned extra >=4 with g>2e-5 Ha and r>=0.05
at every k. Disputed or reference-unqualified spectra are ineligible candidates;
fallback four if no eligible candidate is available.
The 0.05 threshold was listed in the preregistration. This oracle knows the
converged spectrum and pays no probe cost, so even a pass would not authorize
Phase B or a recommendation. Its role is to test whether gap closure alone is
a useful sufficient criterion. Do not retune it after held-out results. Report
training and held-out results separately. It does not exclude more elaborate
future predictors if it fails.

A second ignored test `band_boundary_spectrum` reads the exact saved tight
reference density and requests occupied+17 roots (or the largest constructible
count in a limited basis). It checks the common prefix against that reference
and exports all residuals. This resolves disputed N+1 probes without changing
SCF or transferring spectra between densities. `wide_gap_ha` and
`wide_relative_gap` are separate columns/correlations; the original probes and
preregistered measures remain retained. A wide probe with prefix disagreement
>1e-4 Ha, residual >1e-5 Ha, or no N+1 root remains unknown. This still does not
prove dense spectral completeness; no optimization is accepted from it.

Rank sensitivity is also reported after collapsing absolute gaps <=2e-5 Ha to
one unresolved rank. Ranking roundoff-scale splittings as distinct gaps can
produce meaningless correlations. The absolute/relative raw-rank correlations
remain reported, alongside this resolution-aware sensitivity and the previously
fixed gap-group comparisons. Pooled correlations are confounded by case.

## Findings: boundary closure is not necessary or sufficient

The motivating Si8 SCF distribution is reproduced exactly, but the physical
interpretation in the earlier orchestrator scan does not survive the spectrum
measurement. Si8 extra=8 is **inside** the next three-member empty multiplet,
not at a gap. A wider probe at the saved tight reference density gives a gap
8.48e-12 Ha (complex) / 4.85e-12 Ha (real), unresolved at the 1e-5 Ha residual
target, while all ten standard starts take ten SCFs. The gap at extra=6 is
0.0719762645 Ha, yet its real-path SCF range is 24–30, median 25; the complex
path takes 10–11, median 10. Thus a large converged gap does not guarantee
stable behaviour on both arithmetic paths, and a gap is not necessary for
stable counts. A first-iteration or subspace-size mechanism remains possible;
this study does not attribute the effect to a specific Davidson operation.
Changing the band count also changes the starting block and the default
Davidson capacity (max(48, 6*bands)); the scan does not isolate those effects
from the location of the spectral cut. No capacity experiment is repeated here.

| Si8 path | Extra | Standard SCF min/median/max | Median H vectors | Wide boundary gap (Ha) |
|---|---:|---|---:|---:|
| complex | 4 | 20/25.5/32 | 1167.5 | 6.83e-12 |
| complex | 6 | 10/10/11 | 585.5 | 0.0719762645 |
| complex | 8 | 10/10/10 | 629 | 8.48e-12 |
| complex | 10 | 9/9/9 | 593 | 0.1550245601 |
| complex | 12 | 10/10/10 | 692 | 1.82e-13 |
| real | 4 | 21/27.5/30 | 1281.5 | 1.30e-11 |
| real | 6 | 24/25/30 | 1232 | 0.0719762644 |
| real | 8 | 10/10/10 | 609 | 4.85e-12 |

“Same trajectory” requires care: Si8 complex extra=8 has one integer work
sequence but ten distinct exact scalar histories. The original perturbations
remain numerically distinguishable. Equal iteration counts alone are not
bit-identical trajectories; the raw histories preserve those distinctions.

Primitive Gamma gives another counterexample: extra=1 has an unresolved gap
on both paths, but complex is 11/11/11 SCFs and real 62/65/68. Low-cutoff Si
similarly gives complex 9/9/9 versus real 52/61/69. No threshold on the final
gap can distinguish these pairs by itself. Defaults have not changed.

The frozen descriptive oracle already fails its work bound: training
syn-cubic-degenerate selects five extras and raises median H 98 to 112
(+14.29%) on both paths, with unchanged 7/7/7 SCFs. Held-out si-k222s selects
six extras, raises median H 731 to 1083 (+48.15%), and worsens SCF 9/9/9 to
11/11/11. Its two strict-reference H ratios are 1.1740 and 1.2337 (these are
not ten-start ordinary strict distributions). These optimistic oracle costs
exclude the spectrum acquisition itself. A production first-loose predictor
has neither evidence of reliable selection nor budget to ignore its probe cost.
No auto option is justified from these observations; Phase B remains unentered.

Complete per-condition records, including every failure, original/reference
spectra, wide probes, density bytes, accuracy and exact histories:
[P42A-evidence/summary.csv](../../P42A-evidence/summary.csv),
[summary.json](../../P42A-evidence/summary.json),
[reference-spectra.csv](../../P42A-evidence/reference-spectra.csv),
[spectra.csv](../../P42A-evidence/spectra.csv),
[wide-spectra.csv](../../P42A-evidence/wide-spectra.csv),
[correlations.json](../../P42A-evidence/correlations.json),
[relations.json](../../P42A-evidence/relations.json),
[oracle.json](../../P42A-evidence/oracle.json).
Collection and the full coverage audit are complete; the linked tables include
all attempted conditions, not only successful or reference-backed rows.

## Commands and evidence interpretation

Run from the repository root, with the existing sibling external pseudo files
(`../pseudo/Si.pz-vbc.UPF` and `../pseudo/pseudodojo/ex/{Si,Mg,O}.upf`). The
collector generates synthetic UPFs only in temporary directories, retaining
the generator inputs and hashes, never committing or copying UPFs here.

```sh
cargo test --release -p dft-engine --test band_boundary_scan --no-run
python3 -I P42A-evidence/collect.py train
python3 -I P42A-evidence/large_partition.py
python3 -I P42A-evidence/collect.py heldout
cargo test --release -p dft-engine --test band_boundary_spectrum --no-run
python3 -I P42A-evidence/wide.py
python3 -I P42A-evidence/analyze.py
python3 -I P42A-evidence/relations.py
python3 -I P42A-evidence/oracle.py
python3 -I P42A-evidence/tables.py
python3 -I P42A-evidence/decision.py
python3 -I P42A-evidence/audit.py
cargo build --release -p dft-cli
python3 -I P42A-evidence/identity.py
python3 -I P42A-evidence/audit_identity.py
scripts/check.sh
python3 -I scripts/make_readme_figures.py
python3 -I scripts/tests/test_readme_figures.py
```

Collectors are resumable: completed execution records are not overwritten.
The wide collector is rerun after the SCF collectors finish so every available
tight reference density is audited. The full evidence audit rejects missing
conditions/runs, unexpected process failures, invalid success/certification
labels, missing density bytes, inconsistent work totals, and gap-index errors.
Raw evidence contains 17-digit scalar histories; density arrays are little-endian
f64 in Bohr^-3. Accuracy density columns use RMS sample differences (Bohr^-3),
not the engine's electron/volume-normalized SCF convergence residual. Eigenvalue
errors use the occupied prefix at every sampled k, in Ha; energy errors are Ha.

The identity matrix has 70 supported repeat pairs and 30 expected forced-real
multi-k rejections. All supported pairs have bit-identical energy, density hash,
SCF/Davidson/H counts and certificate counts. `bench_scf --require-identical`
checks are supplemented by explicit floating-point bit comparisons in
`identity/exact-bits.json`. These compare the same unchanged production build,
not two independent implementations. There is no production diff to benchmark.

The first workspace gate passed. After the wide probe was added, a second run
failed clippy's `chunks_exact_to_as_chunks`; the byte reader was corrected to
`as_chunks::<8>()` and the full build/test/fmt/clippy gate passed again (exit 0).
README generation and nine freshness tests pass (exit 0); generated README
and figures have no diff. This is a measurement-only negative result, not a
new production validation claim requiring a README evidence type.
No artifacts/*, dependency, checkpoint, option parser or numerical source changed.

Si16 complex strengthens the counterexample. Extra=8 has a resolved wide gap
0.0543755783 Ha (relative 3.3481), yet only 9/10 starts converge, with successful
SCF 39/44/48 and all-attempt median H 3956.5. Seed 5 diverges after 31 steps;
all ten default-extra=4 starts converge (18/21.5/27, H 1854.5). Extra=10 has
an unresolved wide gap 7.13e-13 Ha, yet all starts take ten SCFs (H 1078).
Both wide spectra agree with their original reference prefixes within 1.1e-10
Ha. This directly contradicts a sufficient monotone gap-stability explanation.

Some strict-reference **buffer** roots themselves omit states: for Si16 complex
extra=12 the reference's band 44 is 0.3477658054 Ha, while the wider frozen-density
probe has a two-member group near 0.3372719371 Ha. Occupied roots and total energy
remain consistent. The strict certificate does not promise completeness of all
empty buffers. Such boundary rows are flagged unknown, not silently relabelled
as physical shell closure or used as reliable predictor training data.

Correlations/group comparisons require all ten standard starts to converge and
a valid two-reference plateau. Failed/censored conditions are reported separately
and are not assigned a zero spread. This successful-condition statistic alone
must not conceal the explicit new-failure counterexample at Si16 extra=8.

Large-fixture scheduling amendment before its first P42A execution: the existing
P38A timing records show syn-100-bands takes 25–131 s per single run even at a
100-step budget. Its P42A default 200-step, 216-run scan would take hours
serially. `large_partition.py` assigns six disjoint three-condition shards,
each still T=1, alongside the training/material collector. It writes explicitly
labelled non-result reservations so the original sequential collector skips
those assigned conditions. `analyze.py` excludes all reservations; the final
audit requires actual execution for every condition. No reserved row is counted
as measured and no wall-time inference uses this loaded execution. Scientific
parameters, split and decision thresholds are unchanged. Commands, shard IDs,
loads and actual exit codes remain in P42A-evidence/partition*.json and logs.

The material training case also separates iteration savings from work/accuracy:
PseudoDojo Si extra=1 keeps 11/11/11 SCFs and reduces H 6974 to 4439, but its
maximum energy error against the default tight reference is 1.476513e-9 Ha
versus 2.593357e-10 Ha at extra=4 (about 5.69x, above the 1.5x limit).
Extra=6 lowers SCFs to 10/10/10 while raising H to 7947 (+13.95%) and energy
error to 5.610090e-10 Ha (about 2.16x). These are descriptive scan rows, not
selected production policies. Small SCF counts alone cannot establish a win.

The saved identity records also pass a cross-budget audit: all 14 supported
(profile, path, case) groups have exactly identical energy bits, density hashes,
SCF/Davidson/H and certificate counts at T=1,2,4,8,12. This uses the existing
records, not additional SCF runs (`identity/cross-thread.json`).

Fixture adapter details: the unchanged P29A occupation policy requires a sampled
insulating gap >1e-6 Ha. Unspecified iteration budgets use the engine default
200 (not bench_suite's CLI default 100); material inputs retain 120 and explicit
synthetic budgets, where present, are honored. These settings explain why this
scan must not be merged numerically with P38A's seed-42, 100-step CLI scan.

### Unchanged-harness fidelity check of new failures

Si16 real extra=8 has 7/10 successful starts (41/46/50 SCFs among successes,
median all-attempt H 4127); extra=10 has 0/10 successes, with density-growth
rejections after 10–17 iterations. Its strict references still converge.
These differ from treating a constant short attempt length as stable convergence.
The unchanged `standard_stall_trajectory` test was built and executed separately
for both conditions using the same P29A flags. Every ROW (all scalar values and
work counts) and every END status, including both strict references, is exactly
identical to the new harness. Thus the failures are not introduced by its added
spectrum/output diagnostics. Evidence: `fidelity/comparison.json`, both raw logs,
and `fidelity.py` (exit 0). These are current-build local results, not replacements
for differently sourced historical/relayed distributions.

### Quantified gap–spread relation

For wide, own-density probes with valid prefix agreement, a valid reference
plateau and ten converged standard starts, Spearman correlations after collapsing
all gaps <=2e-5 Ha to one unresolved rank are:

| Case | Path | Valid conditions | Absolute gap vs spread | Relative gap vs spread |
|---|---|---:|---:|---:|
| Si8 | complex | 13 | -0.487739 | -0.433167 |
| Si8 | real | 13 | -0.556813 | -0.497508 |
| Si16 | complex | 7 | +0.081786 | +0.204465 |
| Si16 | real | 6 | -0.516100 | -0.394665 |

These are descriptive within-case associations over extra-band choices, not
independent population samples or a causal law. The direction is inconsistent
on Si16 complex, and the explicitly failed conditions are outside these
successful-only correlations. Raw-rank coefficients, median-H correlations,
all other cases, pooled/confounded coefficients and the fixed absolute/relative
threshold groups are retained in correlations.json and relations.json. Many
small-case spreads are constant zero, for which Spearman is undefined rather
than zero. The substantial unresolved/omitted-root scope is visible in the
per-condition spectrum validity fields; it is not silently filled from another
case or another fixed point.


## Frozen decision rule: failed conditions and limits

The uncharged, converged-spectrum diagnostic is **not** a deployable predictor.
Its final held-out standard results are in
[decision.json](../../P42A-evidence/decision.json) and
[per-case oracle records](../../P42A-evidence/oracle.json).

- Median H ratio is **1.000** on both paths: the required median reduction fails.
  Per-case H increases breach 10% on k222s (+48.15%), MgO (+37.87%), and syn-z1
  (+13.95% complex, +18.60% real).
- Median spread change is **0 iterations** on both paths, not a reduction.
  No fully converged diagnostic pair exceeds the 20% spread-increase bound,
  but that does not satisfy the separate shrinking-median requirement.
- Of the fully converged paired holdouts, SCF median does not rise on 12/13
  complex and 7/7 real cases. This successful subset is not a claim that the
  full-case 80% gate passes; failed/undefined comparisons remain listed.
- MgO selects ten extras: SCFs stay 13/13/13, H rises 5909 -> 8147, and maximum
  energy / occupied-eigenvalue / density errors are respectively **5.9995x,
  2.5067x, 2.8943x** default. Every paired seed violates all three 1.5x bounds.
  The Gaussian-fcc fixture also has paired-seed accuracy violations, despite
  better aggregate work/error statistics. No roundoff floor is added.
- The large synthetic real case selects ten extras and improves the success
  count 6/10 -> 8/10, but seed 2, which converged at default, now fails.
  More successes in aggregate do not preserve every previously convergent
  perturbation. Its strict-reference H ratios are 1.6556 and 1.5898; they are
  reference diagnostics, not ordinary-strict ten-seed acceptance distributions.
- Some reference or paired-start errors are unavailable. Those rows cannot
  count as accuracy passes. No ordinary-strict candidate ensemble or first-loose
  predictor is claimed, because the Phase A premise did not justify Phase B.

The median-spread requirement is particularly restrictive for this split:
11 of 13 fully convergent complex holdouts and 6 of 7 fully convergent real
holdouts already have zero default spread. Zero must remain zero under the
per-case 20% bound. Their median cannot strictly shrink in those strata;
failed defaults do not supply a defined successful spread. We retain the
original criterion rather than redefining the aggregate after observing this.

Across the entire integer scan (not an accepted policy), the largest measured
standard errors against the default tight references are 5.042674e-3 Ha in
energy (syn-s2 real, extra=1), 3.226893e-3 Ha in occupied eigenvalues (syn-s2
real, extra=5), and 3.123424e-2 Bohr^-3 density RMS (syn-s2 complex, extra=2).
These fixed-point-sensitive toy rows must not be advertised as equal-accuracy
speedups. Per-run errors and reference density/eigenvalue plateaus are retained
in summary.json; no failed reference is assigned zero error.

Final commands in final-analysis-exits.json all exit 0: wide probes, analysis,
relations, oracle, tables, decision reporting, full audit and exact-bit audit.
A preliminary partial audit correctly rejected a stale summary while a wide
collector was finishing; regenerating the summary resolved it without changing
raw measurements. The final audit was run only after all collectors completed.

Limits: spectra are residual/prefix-audited, not independently dense-certified
lowest-root completeness proofs; convergence and accuracy are unavailable where
reported; Mac loaded-host wall times establish no speed claim. The early loose
subspace mechanism, possible kinetic/first-iteration predictors, and dev08
performance remain unestablished. This negative result rules out recommending
the tested closure diagnostic, not every possible future predictor. Defaults
remain unchanged and all changes are left uncommitted for review.
