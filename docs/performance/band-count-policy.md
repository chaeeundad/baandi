# P38A band-count policy study

**Decision: no automatic option is justified by this study.** The default stays
integer four and no production source is changed. Both the converged first-empty
group oracle and the simple free-electron proxy violate the prespecified
per-case work bound. This does not exclude a future, independently validated
first-iteration predictor; it does not authorize fitting one to these held-out
results after inspecting them.

## Preregistered protocol (2026-11-07, before new SCF measurements)

Phase A only until the evidence supports a broadly useful rule. Integer default
four, solver, tolerances and real-path implementation remain unchanged. The new
premise relative to rejected A19 is conditional multiplet completion, not a
universal reduction to one or two buffers. A08/A13/P31A capacity tuning and A21
tolerance schedules are not repeated.

Training cases: si-gamma, si-k222, si-k444s, si8-gamma, si-skew, syn-s1,
syn-p3, syn-cubic-degenerate and PseudoDojo Si (si-20 material input).
Held out: si-k222s, si16-gamma, si-elongated, si-slab, all remaining syn-*
fixtures and PseudoDojo MgO (mgo-40 material input). Held-out scans may be
reported descriptively but must not tune a candidate. Any Phase B evaluation
also includes seeds 1–5 of the P29A mean-centered LCG initial-density noise
(relative amplitude <=1e-10) on si8 and si16. No candidate has been selected.

Adoption requires reduced median H vectors over held-out cases and perturbations,
no individual increase above 10%, no newly failed default-convergent case,
unchanged occupied-band tolerance class, strict certificate success, and
energy error against tight references no larger than 1.5 times default error.
Unknown or failed reference accuracy cannot count as a pass. No default change.

Phase A: T=1, complex path, strict and standard, extras 1–12, three runs each;
record process wall seconds, load before/after, SCF steps, Davidson steps,
main H vectors and certificate H vectors separately, converged energies in Ha,
and per-k eigenvalues/occupations. Use original benchmark geometries and
material cutoffs/meshes. Explicit fixture nbands is replaced by occupied+extra
for this experiment. Retain errors and timeouts, never assign failed runs a
converged energy. All raw evidence and the resumable driver live in
`P38A-evidence/`; UPFs remain external or temporary and are never committed.
The existing `bench_suite.case_input` adapter receives seed 42 for benchmark
and synthetic fixtures; material input files retain their omitted seed (engine
default 2026). The seed is identical across band counts and repetitions within
each case. “Default-band result” means extra=4 under that case's recorded input,
not an assertion that every input option equals an unrelated historical run.
Material fixtures are the established si-20/mgo-40 inputs (40/80 Ry, respectively),
with their PBE, XC floor, meshes and mixing retained and threads changed to one.

Multiplet analysis must distinguish the first-empty multiplet from the multiplet
intersected by the **last requested** band. A gap at the truncation boundary
requires at least one eigenvalue beyond that boundary. A spectrum ending at
occupied+12 alone cannot certify completion there. Strict buffer residual targets
have a 1e-5 Ha floor; standard uses Ritz-change stopping without a strict
certificate. Tiny buffer splittings are not independently certified degeneracy.
Report gap sensitivity and unresolved classifications, and extend the spectral
probe beyond 12 when required. Correlation of converged gaps with work alone
does not prove a reliable first-iteration predictor.
Reported multiplet counts are the **empty members needed above the occupied
prefix**, not the total multiplicity when a grouping threshold also spans the
HOMO–LUMO gap. The spectra explicitly flag that occupied-boundary crossing;
near-degenerate grouping at 1e-4 Ha is not a claim of exact symmetry degeneracy.

Mac wall times are observations only; dev08 acceptance remains through the
orchestrator. Existing scan and QE Q1(a) results are historical/relayed evidence,
not new measurements. A17 identifies a six-member cluster cut at band 20;
A19 already records failure and fixed-point changes on a larger skewed fixture.
Thus neither band-count monotonicity nor causality is assumed.

Scheduling amendment during collection: after the Si/shape scans, the expensive
nonconvergent `syn-100-bands` scan and the remaining synthetic/material scan run
in two independent processes, each still T=1. Exact transition time and case
assignment are in `P38A-evidence/scheduling.json`; each sample records load before
and after. No builds overlap timing collection. These are deliberately labelled
loaded-host observations, not isolated single-core speed acceptance. Reserved
rows are explicitly unmeasured until a process writes its actual result; final
coverage must exclude reservations. No scientific parameter or split changed.

The full default identity gate is gamma/k222s/k222/k444s/si8 ×
T=1/2/4/8/12 × strict/standard × gamma_real off/on. Ineligible forced-real
k-point combinations must be reported as rejected, not successful comparisons.
The workspace gate is `scripts/check.sh` (without an argument: literal `all`
would request a nonexistent crate). No artifacts/* changes are authorized.

## Findings

The simple converged-spectrum oracle “retain the entire first-empty multiplet,
and no further empty bands” already fails the per-case work gate on held-out
Si16. At every exploratory gap threshold 1e-7, 1e-6, 1e-5 and 1e-4 Ha, its
first-empty multiplet has two members. Strict H work at extra=2 is 3205 versus
2657 at default extra=4 (**+20.62%**); SCF steps are 14 versus 12. These are
three repeated local runs, not dev08 acceptance measurements. This is an oracle
diagnostic, not an implemented candidate: it uses information unavailable before
convergence and does not charge for obtaining it. A first-iteration predictor
would additionally need its own evidence and cost accounting.

Si8 behaves differently: the first-empty multiplet has six members across those
four thresholds. Strict main H vectors are 2063 at extra=4, 1253 at extra=6,
and 1024 at extra=7 (the minimum of the 1–12 sweep). SCF steps are respectively
15, 12 and 10. Standard minimum is 588 at extra=8 versus 742 at extra=4.
This supports the multiplet-boundary mechanism on Si8, but completing that
multiplet is neither the exact optimum nor a general sufficient rule.

Primitive Gamma has a three-member first-empty multiplet but its minimum main
H count is at extra=1. The observational scan therefore must not be presented
as proof that all incomplete empty multiplets cause stagnation. The allowed
1e-5 Ha buffer residual floor also limits spectral interpretation.

A second diagnostic uses the pre-SCF free-electron shells of cubic Si8. Integer
Miller squared norms 0, 1 and 2 have multiplicities 1, 6 and 12, ending at bands
1, 7 and 19. With 16 occupied bands, completing the first-empty free-electron
shell selects extra=3, not the physical extra=6. Its measured strict work is
2321 versus default 2063 (**+12.51%**). This is a training-case diagnostic of a
fixed physical hypothesis, not an independently validated candidate; no
parameters were fitted. See `P38A-evidence/free-electron-si8.json`. Neither this
pre-SCF proxy nor the converged first-empty-multiplet oracle passes the stated
per-case work bound. A more elaborate first-iteration predictor remains untested.

### Non-cubic cells

The Si16 benchmark itself is a 1:1:2 orthorhombic supercell. Its wide strict
spectrum separates a two-member first-empty group (bands 33–34) from a four-member
group (35–38) by about 4.12e-4 Ha. Default band 36 cuts the latter, but completing
it at extra=6 **increases** main H work from 2657 to 2754. Extra=10 minimizes
work at 1960 despite cutting the next high doublet (42–43); completing that
doublet at extra=11 costs 1988. Common spectra agree within 4.72e-8 Ha for
these rows, so this is not simply an unmatched high-band list. These observations
contradict a universal rule that completing the last-band multiplet reduces work.
They do not rule out a transient first-iteration cluster effect; no such
predictor was measured here.

The sheared primitive `si-skew` has a nondegenerate first empty state at all four
gap thresholds. Strict minimum H work is 296 at extra=2 versus default 376;
standard minimum is 170 at extra=1 versus default 304. It does not reproduce
Si8's large completion benefit. `si-elongated` and `si-slab` exhaust their
100-step budgets at the default in both profiles, so no converged-state
degeneracy conclusion can be drawn from those trajectories. Failed runs remain
in the tables, not in accuracy or successful-speedup statistics.

### Synthetic scope and accuracy counterexample

The seven-plane-wave fixtures cannot represent occupied+12 bands. For example,
`syn-s1/s2/p3/d3/spd/fft-primes` have two occupied bands and at most five empty
bands in their fixed basis; extras 6–12 are rejected by unchanged preflight.
Their widest available converged spectra are retained, with the missing range
explicitly marked unavailable. Cutoffs were not raised to conceal this limit.

Held-out `syn-s2` strict converges at extras 1 and 2 to energies respectively
4.416513597409022e-3 and 4.416513597401028e-3 Ha above the extra=4 result.
Main H work also increases (165/152 versus 115). Extras 3 and 5 agree with the
default energy to 6.22e-15 Ha. This is fresh evidence consistent with the earlier
[fixed-point multiplicity narrative](../validation/scf-fixed-point-multiplicity.md),
not a waiver of the ticket's accuracy gate. The common-spectrum discrepancy is
about 2.89e-3 Ha at extras 1/2, so their labels were initially left unknown rather
than transferred between different states. The frozen-density probes below
resolve them using their own converged densities. Repeated outcomes are
deterministic within each condition.

## Complete scan coverage

The scan contains **1728 measured runs**: 24 cases × two profiles × twelve extra
band counts × three repetitions. All 576 conditions have three actual runs.
There are 1095 converged runs and 633 unsuccessful attempts: 294 preflight
basis/band incompatibilities, 306 SCF-budget exhaustions and 33 density-divergence
errors. Failed attempts retain their available history and work; their energies
are not reported as converged values. The full condition table is
[summary.csv](../../P38A-evidence/summary.csv), and the compact per-case table is
[tables.md](../../P38A-evidence/tables.md). Timing medians are in seconds and
energy differences are in Ha. Recorded one-minute load spans 2.6733–29.6738;
wall observations cannot establish isolated-core speed improvements.

The large synthetic case converges in standard only at extras 11 and 12
(three of three each): respectively 99/66 SCF steps, 482/358 Davidson steps and
24298/17084 main H vectors. Default extra=4 fails, so no default-relative
converged energy or tight-reference accuracy is available. The widest converged
spectrum needs only one empty member to finish the first-empty group, despite
the need for many more buffers for SCF convergence. All 36 strict attempts fail
at the unchanged 100-step budget. Standard convergence is not certification.

The PseudoDojo strict minima are at extra=1: Si 11223 versus default 15364 H
vectors (energy difference -8.88e-15 Ha); MgO 9085 versus 10930 H vectors
(-2.56e-13 Ha). These useful cases do not erase the held-out work, accuracy and
convergence counterexamples. No default or recommended opt-in policy follows.

### Frozen-density resolution and audit (completed 2026-11-08, local date)

The ignored release test passes. All five replays reproduce the original energy
bits and density SHA-256 exactly. The 36 fixed-density k-point spectra have
maximum recomputed residual 9.9977605473e-6 Ha, within the unchanged diagnostic
target of 1e-5 Ha. Results and density identities are in
[frozen-spectra.json](../../P38A-evidence/frozen-spectra.json).

| Original condition | Empty members needed, maximum over k | First-empty group complete? |
| --- | ---: | --- |
| Si16 standard, extra 3 | 2 | yes |
| Si16 standard, extra 4 | 2 | yes |
| PseudoDojo Si standard, extra 2 | 1 | yes |
| syn-s2 strict, extra 1 | 2 | no |
| syn-s2 strict, extra 2 | 2 | yes |

Every successful condition now has an explicitly sourced grouping flag at
1e-4 Ha. All observed closing boundaries require at most six empty bands, so no
extra value outside the requested 1–12 range was needed. Failed SCFs have no
converged spectrum; the seven-plane-wave cases cannot supply the unavailable
upper range. These limitations are retained, not filled with extrapolated data.

The independent raw-result audit passes: all 1728 records are accounted for and
available histories reconcile exactly; repeated successful densities and
numerical outcomes are identical.
All 543 converged strict runs carry their certificate, with maximum certified
occupied/first-empty residual 9.9881473239e-10 Ha. No production source changed
and no UPF exists in the evidence directory. The 70 supported default comparisons
also match T=1 bit-for-bit across thread budgets. See
[audit.json](../../P38A-evidence/audit.json) and the verification logs.

Raw condition tables: [summary.csv](../../P38A-evidence/summary.csv), with
separate main/certificate work, repeat counts, load ranges, wall medians, energy
and occupied-eigenvalue differences against the same-profile default. Failed
attempt wall times, when present, are censored costs rather than speedups.
[spectra.json](../../P38A-evidence/spectra.json) retains per-k gaps and the
four-threshold multiplet sensitivity; no buffer spectrum is called certified.
The portable [spectral CSV](../../P38A-evidence/spectra.csv) records band indices,
eigenvalues, adjacent gaps and available residuals in Ha. Upper buffer values can
differ substantially between band counts even when occupied energies are close;
these are recorded SCF Ritz spectra, not independent proofs that every low root
was found. The completeness flag concerns only the first-empty group, with its
data source stated, and does not certify all higher multiplets.

## Reproduction and verification

All commands below use the one existing `target/` directory. The evidence
directory is outside target and contains no UPFs. Source/binary/UPF hashes and
platform are in `P38A-evidence/manifest.json`; complete generated inputs and
results.json files accompany each sample. Synthetic UPFs are generated only in
a temporary directory from the existing benchmark generator.

```
cargo build --release -p dft-cli --offline
python3 -I P38A-evidence/scan.py
python3 -I P38A-evidence/summarize.py
python3 -I P38A-evidence/report_tables.py
python3 -I P38A-evidence/identity.py
python3 -I P38A-evidence/frozen_probe.py
scripts/check.sh
python3 -I scripts/tests/test_readme_figures.py
python3 -I P38A-evidence/audit.py
```

The collector is resumable and leaves successful/failed measurements intact;
fresh execution is serial unless `P38A_ONLY` explicitly partitions cases.
The first collector attempt exited 1 on Si16 standard extra=8 because failed
CLI results omit `total_energy`. The parser was repaired to retain failure
history without inventing a converged energy, and raw execution information is
now written before summary parsing. The original failure JSON and attempt log
are preserved. This was an evidence-collector failure, not a solver repair.

Final workspace build/test/fmt/clippy passes with the ignored spectrum harness
included (`scripts/check.sh`, exit 0), and README freshness passes nine tests
(exit 0). The release spectrum test and its independent density/energy identity
checks pass (exit 0), as does the complete raw-evidence audit (exit 0).
The default matrix has 70 supported
same-HEAD baseline/comparison pairs, each with bit-identical energy, density and
work counts. All 30 forced-real multi-k combinations reject the unsupported
input; they are not counted as passing comparisons. `identity-bits.json` checks
exact energy bits in addition to bench_scf's tolerance/hash gate. No production
source or default changed, so these are repeated stability checks of the same
executable, not claims about an independent modified build.

The main CLI study is an explicit external measurement command. The additional
`band_count_spectrum` test is ignored in debug/default gates and run explicitly
in release. It replays the five ambiguous converged rows, exports their density
bytes for independent SHA-256 comparison against the original CLI metadata, and
uses the existing fixed-point API at that frozen density. The diagnostic uses
equal mandatory/buffer residual targets of 1e-5 Ha, checked for every root,
below the 1e-4 Ha grouping threshold. Original SCF tolerances are unchanged;
diagnostic costs are not included in the SCF timing/work table. Build and run
this test only after timing collection to avoid build interference.

Phase B perturbations, tight-reference
error ratios, first-iteration predictors and restart-option tests are conditional
on selecting a justified candidate; none has been selected. No auto option or
recommended band-count policy is claimed.
