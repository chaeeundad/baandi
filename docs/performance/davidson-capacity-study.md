# P31A Davidson capacity study

## Prespecified decision rule (before measurement)

Keep the default at 6 in this ticket. A smaller multiplier is a candidate for
future default review only if median wall time over perturbations and over all
cases falls, H-vector work increases by no more than 10%, strict delivered
accuracy retains the 1e-9 certificate, standard energy error is no worse than
1.5 times the multiplier-6 error, and no case converging at 6 fails at the
candidate. Otherwise classify it as case-dependent. Mac timings are observations,
not dev08 acceptance or a new QE comparison.

A08/A13 tested eight bands. The new premise is the measured growth in per-vector
cost on 20–68 bands in the relayed dev08 evidence, consistent with increasing
subspace algebra cost. This reopens B17 without reversing those earlier results.
Derived from the prose study note `docs/research/qe-study/qe-solvers-and-subspace.md`;
no QE source was accessed.

Capacity is min(plane waves, max(absolute subspace floor, multiplier * bands)).
The absolute floor remains 48, so small-band cases may see no capacity change.
Multipliers must be integers 2 through 12. Minimum 2 reserves a full correction
block beside the current bands; a basis smaller than twice the bands is still
allowed because the full basis can be reached. Default 6 preserves arithmetic.

Raw evidence, build logs, binaries and measurement scripts reside in
`P31A-evidence/`, outside target. No UPF is copied. Only one target directory is
used. Measurements and gate results will be added after execution.

## Implementation and measurement semantics

Research input (omission has the same value):

```toml
[scf]
davidson_capacity_per_band = 6
```

`ScfOptions::davidson_capacity_per_band` defaults to 6. CLI accepts only integer
2..=12, records the chosen multiplier in result options (and the existing
`davidson_subspace_min_per_band` key), and uses it in preflight. All scalar,
collinear, spectral-audit and real-Gamma solve capacity floors use this field.
The real solver's existing public entry points retain multiplier 6; its additive
research entry point accepts the same multiplier. Exact restart compares it;
v1 checkpoint omission decodes to 6 and default serialization omits the field.
The sole additional production plumbing file is `dft-engine/src/checkpoint.rs`.

The ignored `capacity_trajectory` test derives from the P29A harness. Inputs are
read directly from benchmarks/cases and tests/reference/materials, with external
UPFs. Materials retain PBE, their documented density floor, cutoff and mixer;
all runs use T=1. Strict references use multiplier 6 at density/energy/eigen
thresholds 1e-9/1e-11/1e-10 and then 1e-10/1e-12/1e-11 (atomic units).
Measured runs use unchanged strict defaults or unchanged standard defaults.
Runs 2..6 use five independent mean-centered LCG density perturbations <=1e-10;
other conditions use run 2. Runs 0 and 1 are unperturbed reference checks.

Wall seconds are uninstrumented harness elapsed time, including two problem
constructions needed for the density-transfer validation API, initialization and
SCF; output formatting/I/O is excluded. They are not CLI wall timings. Work
counters are collected in a separate `P14A_WORK_ACCOUNTING=1` build; its wall
values are not used for speed comparisons. Raw rows distinguish completed SCF
history work, certificate work, and total solver work including failed solves.
Failures are retained; their last available energy is never a delivered result.
Same SCF count alone does not establish the same trajectory, so both count equality
and recorded scientific-history equality are recorded against the paired multiplier-6
run. Timing comparisons are local observations, not a dev08 acceptance claim.

The absolute floor means 2, 3, 4 and 6 all yield 48 columns for eight bands.
No implementation attempts to force those small cases below the absolute floor.

## Reproduction

From the task worktree, using only its existing target directory:

```sh
scripts/check.sh
cargo test --release -p dft-engine --test capacity_trajectory --no-run
# Copy the printed test executable to P31A-evidence/study-normal.
python3 -I P31A-evidence/study.py P31A-evidence/study-normal normal
P14A_WORK_ACCOUNTING=1 cargo test --release -p dft-engine --test capacity_trajectory --no-run
# Copy the printed executable to P31A-evidence/study-account.
python3 -I P31A-evidence/study.py P31A-evidence/study-account account
python3 -I P31A-evidence/compare_accounting.py
python3 -I P31A-evidence/summarize.py
RUSTFLAGS='--cfg feature="worker-analysis"' cargo build --release -p dft-cli
# Copy the instrumented CLI to P31A-evidence/profile-binary.
python3 -I P31A-evidence/profile.py P31A-evidence/profile-binary
# Combined named-site coverage, including reduced_extension:
bash P31A-evidence/full_profile.sh
python3 -I P31A-evidence/validate.py
```

The runner resolves external pseudopotentials under `../pseudo/` and
`../pseudo/pseudodojo/ex/`; edit only those external paths if reproducing on
another host. Its commands JSON records every condition and exit code. A sweep
test exit 0 means the measurement completed, not that every SCF converged;
`END` records and the row status carry the scientific outcome. `hashes.txt`
records UPF/input/binary identities, and `head.txt` records the source baseline.
Raw evidence should be archived with the review; binaries and UPFs are not
source deliverables. No new QE run or remote timing is claimed.

## Decision

**Keep 6. Smaller capacities are case-dependent and none meets the prespecified
candidate-default rule.** The uninstrumented sweep completed 152 runs (including
18 reference checks): 137 converged, 15 failed explicitly. All failures were at
multiplier 2 in strict Si8/Si16: Si8 complex 0/5 and real 0/5 converged; Si16
complex 1/5 and real 4/5 converged. Every corresponding multiplier-6 run converged.
The failures were Davidson's unchanged 400-step limit, not silent acceptance.

For each case/profile/representation, take its perturbation median, divide by
its multiplier-6 median, then take the equally weighted median over conditions:
wall ratios for 2/3/4 are **1.17475 / 1.05999 / 1.00160**. Failed-attempt elapsed
times are included as censored costs, not successful speedups; failure already
vetoes multiplier 2. The 0.16% difference at 4 is below this experiment's timing
resolution and is not a slowdown claim, but supplies no evidence for promotion.
Unchanged-capacity small cases have identical work yet variable wall times,
illustrating shared-host noise. Results do not establish dev08 performance.

Standard median energy-error ratios exceed the 1.5 limit at multiplier 2 for
Si16 complex (2.382), Si16 real (3.502), and Si32 complex (2.155), and at multiplier
3 for Si8 complex (2.104). Standard remains explicitly uncertified. One baseline
Si8 real/standard perturbation has absolute error 1.092e-8 Ha; the stopping metric
is not an absolute energy-error bound. Across all converged measured strict runs,
the largest delivered energy difference is 1.16245e-11 Ha, below 1e-9 Ha; all carry
the unchanged engine certificate. Failed strict runs have no delivered accuracy.

Si32 strict illustrates the tradeoff directly: multipliers 2/3/6 take
143.9721 / 68.7006 / 64.9406 s and 15/13/13 SCF iterations. Including fresh
certificate work, H vectors are 15441/7567/5825, Davidson steps 1859/510/272,
and restarts 208/41/4.
The whole-harness wall divided by total H vectors is 9.32/9.08/11.15 ms
for 2/3/6 (an amortized cost, not the Hamiltonian kernel time). The lower
amortized cost is outweighed by additional H work. Multiplier 3 has the same
SCF count as 6 but performs 29.9% more H work; multiplier
2 performs 165.1% more. Reduced per-step algebra therefore does not translate to
a whole-solve speedup here. A subspace-capacity reduction alone is not supported
as an explanation/remedy for the relayed size-dependent QE gap. This does not
rule out dense algebra as a bottleneck; the named-site attribution below separates
that question from whether shrinking the capacity is an effective intervention.

## Named-site algebra attribution: Si16, T=1, complex Gamma

Combined `worker-analysis` and `P14B_GEMM_ACCOUNTING=1` measurement uses the first
predefined perturbation (run 2), which converges at both multipliers. All four
instrumented scientific histories match their uninstrumented rows exactly.
The failed-perturbation distribution above is not replaced by this successful
profile. These are instrumented observations, not acceptance wall timings or
30 independent SCF repetitions; every named-site sample count exceeds 30.

Algebra time is the interval union of `projection`, `evd`, `gram`, `gram_assembly`,
`rotation`/`rotation_pair`, and `residual_norms`, plus the disjoint
`reduced_extension` GEMM timings. Taking a union avoids double-counting GEMMs
inside projection. The extension is added separately because the older
worker-analysis trace does not span that incremental reduced-matrix branch.
Other copies, preconditioning, bookkeeping and uninstrumented glue remain outside
this named-algebra subtotal. The denominator is the union of Davidson spans,
including certificate solves. Trace/collector overhead can perturb these shares.

| Profile | Multiplier | Davidson s | Named algebra s | Algebra share | H application s | EVD samples | Reduced-extension calls |
|---|---:|---:|---:|---:|---:|---:|---:|
| strict | 6 | 9.2652 | 5.6340 | 60.81% | 3.2668 | 224 | 410 |
| strict | 2 | 37.2943 | 20.1646 | 54.07% | 14.4196 | 2905 | 5026 |
| standard | 6 | 3.4781 | 1.5450 | 44.42% | 1.8144 | 108 | 176 |
| standard | 2 | 4.0212 | 1.5899 | 39.54% | 2.2610 | 138 | 180 |

Strict named algebra per step falls from about 25.15 to 6.94 ms, but the number
of steps grows from 224 to 2905; absolute algebra time therefore grows 3.58x.
Standard reduces algebra per step too, but absolute algebra time is slightly
larger and H application time grows. Thus algebra is a substantial bottleneck,
yet capacity reduction trades it for much more repeated work. Improving algebra
without throwing away the useful subspace is a separate future hypothesis.

Raw traces, dimensions, call counts and timing sums are in
`P31A-evidence/profile-full-*.log` and `profiles-full.json`. An additional
unperturbed CLI profile is retained in `profile-*.csv`/`profiles.json`; its strict
2x run fails at the unchanged Davidson limit and is not used as a converged timing
comparison. The first CLI profile input omitted explicit Gamma k points; that
configuration error was corrected before any profile numbers were used.

## Verification and scope

- `scripts/check.sh`: build, test, fmt, clippy all PASS, exit 0 on final source.
  New tests cover the 2..=12 range, integer-only CLI input, default 6, old checkpoint
  omission, nondefault wire round trips, and exact-restart capacity mismatch.
- Same-HEAD baseline versus final default: 70 energy/density/SCF-count identities
  across strict/standard, T=1/2/4/8/12, all five requested cases with real off and
  Gamma/Si8 with real on. `bench_scf --require-identical` commands all exit 0.
- Forced real on at k222s/k222/k444s is unsupported by the existing engine:
  all 30 combinations reject explicitly in both binaries (60 rejection checks),
  rather than claiming numerical identity for impossible runs.
- All 70 accounting-condition histories match the uninstrumented histories;
  all 152 rows contain complete solver work and restart counts. Failed solver
  work is larger than the last completed SCF-history subtotal and is retained.
- All nine paired tight references agree within 5.82646e-13 Ha. These are internal
  BAANDI reference refinements following the P29 protocol, not independent QE data.
- `make_readme_figures.py` and its 9 freshness tests pass; generated files are
  unchanged because this ticket stores raw research evidence in P31A-evidence,
  not the frozen README reference schema. Artifacts JSON files remain unchanged.
- Mac/pure-Rust observations only. No remote acceptance, no new default, no
  collinear/finite-T performance claim, and no generality claim from single-run
  multipliers 3/4. Raw timing includes two setup constructions and grouped run
  order; small timing differences are unresolved. The extra checkpoint comparison
  is the only production plumbing outside the named implementation files.

Machine-readable evidence validation: `P31A-evidence/validation.json`; command
logs and input/UPF/binary hashes are retained alongside it. Changes are uncommitted
on `task/P31A-DAVIDSON-CAPACITY` for architect review.

## Profile-separated median wall ratios

This avoids hiding a strict regression behind a standard improvement. Each row
is the median over case-specific perturbation medians; real-on covers Si8/Si16
only. Multiplier 4 excludes Si32 as specified by the ticket.

| Multiplier | Strict complex | Standard complex | Strict real | Standard real |
|---|---:|---:|---:|---:|
| 2 | 1.3827 | 1.0054 | 2.0324 | 1.2119 |
| 3 | 1.0579 | 1.0616 | 1.5689 | 0.9662 |
| 4 | 1.0038 | 0.9968 | 1.2414 | 0.9601 |

Ratios near one are noise-level observations. For multipliers 2/3/4 respectively,
SCF counts match the paired 6x row in 9/50, 11/18, 13/16 measured runs; recorded
history fields match exactly in 6/50, 6/18, 9/16. The CSV includes both flags;
matching counts alone does not establish an identical trajectory.
## Measured rows

All individual rows, including failed attempts and reference runs, are in
[rows.csv](../../P31A-evidence/rows.csv) and [rows.json](../../P31A-evidence/rows.json).
Times below are uninstrumented medians (seconds). H/steps/restarts include
certificate and failed-solve work from the separately verified accounting build.
A failed group has censored work/time, not a successful time-to-solution.
Energy error is the median absolute difference from its tightened strict reference
(Ha), computed only for converged runs. n is the number of perturbations.

| Case | Profile | Real | Multiplier | Converged/n | Wall s | SCF range | Davidson steps | H vectors | Restarts | Energy error Ha |
|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|
| material-mgo-40 | standard | off | 2 | 1/1 | 3.7523 | 13–13 | 613 | 5956 | 14 | 5.694e-11 |
| material-mgo-40 | standard | off | 3 | 1/1 | 3.7011 | 13–13 | 613 | 5956 | 14 | 5.694e-11 |
| material-mgo-40 | standard | off | 4 | 1/1 | 3.4385 | 13–13 | 613 | 5956 | 14 | 5.694e-11 |
| material-mgo-40 | standard | off | 6 | 1/1 | 3.4597 | 13–13 | 609 | 5909 | 0 | 1.509e-10 |
| material-mgo-40 | strict | off | 2 | 1/1 | 8.3776 | 18–18 | 1259 | 11443 | 39 | 5.400e-13 |
| material-mgo-40 | strict | off | 3 | 1/1 | 7.8539 | 18–18 | 1259 | 11443 | 39 | 5.400e-13 |
| material-mgo-40 | strict | off | 4 | 1/1 | 6.8151 | 18–18 | 1259 | 11443 | 39 | 5.400e-13 |
| material-mgo-40 | strict | off | 6 | 1/1 | 7.2314 | 18–18 | 1306 | 12202 | 14 | 3.553e-13 |
| material-si-20 | standard | off | 2 | 1/1 | 3.4034 | 11–11 | 1042 | 6974 | 0 | 2.593e-10 |
| material-si-20 | standard | off | 3 | 1/1 | 3.4758 | 11–11 | 1042 | 6974 | 0 | 2.593e-10 |
| material-si-20 | standard | off | 4 | 1/1 | 3.4513 | 11–11 | 1042 | 6974 | 0 | 2.593e-10 |
| material-si-20 | standard | off | 6 | 1/1 | 3.5749 | 11–11 | 1042 | 6974 | 0 | 2.593e-10 |
| material-si-20 | strict | off | 2 | 1/1 | 6.6682 | 15–15 | 2865 | 16526 | 0 | 3.375e-14 |
| material-si-20 | strict | off | 3 | 1/1 | 6.6970 | 15–15 | 2865 | 16526 | 0 | 3.375e-14 |
| material-si-20 | strict | off | 4 | 1/1 | 6.2906 | 15–15 | 2865 | 16526 | 0 | 3.375e-14 |
| material-si-20 | strict | off | 6 | 1/1 | 7.7998 | 15–15 | 2865 | 16526 | 0 | 3.375e-14 |
| si-gamma | standard | off | 2 | 1/1 | 0.0474 | 9–9 | 27 | 180 | 0 | 1.045e-10 |
| si-gamma | standard | off | 3 | 1/1 | 0.0477 | 9–9 | 27 | 180 | 0 | 1.045e-10 |
| si-gamma | standard | off | 4 | 1/1 | 0.0479 | 9–9 | 27 | 180 | 0 | 1.045e-10 |
| si-gamma | standard | off | 6 | 1/1 | 0.0477 | 9–9 | 27 | 180 | 0 | 1.045e-10 |
| si-gamma | strict | off | 2 | 1/1 | 0.1108 | 12–12 | 56 | 353 | 0 | 2.132e-14 |
| si-gamma | strict | off | 3 | 1/1 | 0.0811 | 12–12 | 56 | 353 | 0 | 2.132e-14 |
| si-gamma | strict | off | 4 | 1/1 | 0.0844 | 12–12 | 56 | 353 | 0 | 2.132e-14 |
| si-gamma | strict | off | 6 | 1/1 | 0.0801 | 12–12 | 56 | 353 | 0 | 2.132e-14 |
| si-k222s | standard | off | 2 | 1/1 | 0.1648 | 9–9 | 106 | 731 | 0 | 4.553e-10 |
| si-k222s | standard | off | 3 | 1/1 | 0.1638 | 9–9 | 106 | 731 | 0 | 4.553e-10 |
| si-k222s | standard | off | 4 | 1/1 | 0.1638 | 9–9 | 106 | 731 | 0 | 4.553e-10 |
| si-k222s | standard | off | 6 | 1/1 | 0.1639 | 9–9 | 106 | 731 | 0 | 4.553e-10 |
| si-k222s | strict | off | 2 | 1/1 | 0.3278 | 12–12 | 261 | 1581 | 0 | 2.665e-14 |
| si-k222s | strict | off | 3 | 1/1 | 0.3271 | 12–12 | 261 | 1581 | 0 | 2.665e-14 |
| si-k222s | strict | off | 4 | 1/1 | 0.3269 | 12–12 | 261 | 1581 | 0 | 2.665e-14 |
| si-k222s | strict | off | 6 | 1/1 | 0.3287 | 12–12 | 261 | 1581 | 0 | 2.665e-14 |
| si16-gamma | standard | off | 2 | 5/5 | 4.6526 | 17–23 | 121 | 1855 | 23 | 1.635e-09 |
| si16-gamma | standard | off | 3 | 1/1 | 5.6087 | 27–27 | 108 | 2273 | 1 | 8.673e-10 |
| si16-gamma | standard | off | 4 | 1/1 | 4.6950 | 21–21 | 99 | 1802 | 1 | 9.198e-10 |
| si16-gamma | standard | off | 6 | 5/5 | 5.0908 | 20–27 | 108 | 1926 | 0 | 6.865e-10 |
| si16-gamma | standard | on | 2 | 5/5 | 3.3460 | 19–23 | 121 | 2034 | 24 | 2.025e-09 |
| si16-gamma | standard | on | 3 | 1/1 | 2.8726 | 20–20 | 83 | 1782 | 1 | 2.825e-10 |
| si16-gamma | standard | on | 4 | 1/1 | 2.4846 | 15–15 | 82 | 1452 | 0 | 8.273e-10 |
| si16-gamma | standard | on | 6 | 5/5 | 2.7143 | 15–25 | 91 | 1510 | 0 | 5.781e-10 |
| si16-gamma | strict | off | 2 | 1/5 | 17.4659 | 7–17 | 1263 | 6757 | 183 | 1.016e-12 |
| si16-gamma | strict | off | 3 | 1/1 | 14.4400 | 15–15 | 727 | 5368 | 61 | 5.684e-13 |
| si16-gamma | strict | off | 4 | 1/1 | 12.0698 | 15–15 | 463 | 4233 | 25 | 9.592e-13 |
| si16-gamma | strict | off | 6 | 5/5 | 11.9127 | 14–14 | 224 | 3198 | 4 | 3.325e-12 |
| si16-gamma | strict | on | 2 | 4/5 | 13.9854 | 6–14 | 1523 | 7272 | 186 | 6.033e-12 |
| si16-gamma | strict | on | 3 | 1/1 | 7.4643 | 14–14 | 555 | 4167 | 43 | 5.826e-13 |
| si16-gamma | strict | on | 4 | 1/1 | 7.2605 | 14–14 | 453 | 3730 | 22 | 4.761e-13 |
| si16-gamma | strict | on | 6 | 5/5 | 5.8624 | 12–13 | 206 | 2768 | 2 | 6.182e-12 |
| si32-gamma | standard | off | 2 | 1/1 | 37.2468 | 33–33 | 166 | 5506 | 28 | 2.259e-09 |
| si32-gamma | standard | off | 3 | 1/1 | 27.0129 | 21–21 | 97 | 3625 | 2 | 2.121e-10 |
| si32-gamma | standard | off | 6 | 1/1 | 19.1297 | 14–14 | 70 | 2569 | 0 | 1.049e-09 |
| si32-gamma | strict | off | 2 | 1/1 | 143.9721 | 15–15 | 1859 | 15441 | 208 | 2.274e-12 |
| si32-gamma | strict | off | 3 | 1/1 | 68.7006 | 13–13 | 510 | 7567 | 41 | 7.063e-12 |
| si32-gamma | strict | off | 6 | 1/1 | 64.9406 | 13–13 | 272 | 5825 | 4 | 4.803e-12 |
| si8-gamma | standard | off | 2 | 5/5 | 1.1685 | 22–33 | 111 | 1344 | 10 | 1.099e-09 |
| si8-gamma | standard | off | 3 | 1/1 | 1.1484 | 28–28 | 112 | 1311 | 5 | 1.682e-09 |
| si8-gamma | standard | off | 4 | 1/1 | 1.0855 | 26–26 | 108 | 1211 | 0 | 7.236e-10 |
| si8-gamma | standard | off | 6 | 5/5 | 1.0817 | 21–27 | 108 | 1211 | 0 | 7.995e-10 |
| si8-gamma | standard | on | 2 | 5/5 | 0.7728 | 23–36 | 109 | 1251 | 11 | 1.312e-09 |
| si8-gamma | standard | on | 3 | 1/1 | 0.5671 | 19–19 | 88 | 937 | 2 | 1.840e-09 |
| si8-gamma | standard | on | 4 | 1/1 | 0.6519 | 22–22 | 97 | 1038 | 0 | 2.741e-09 |
| si8-gamma | standard | on | 6 | 5/5 | 0.6488 | 21–29 | 99 | 1038 | 0 | 2.184e-09 |
| si8-gamma | strict | off | 2 | 0/5 | 3.9148 | 6–7 | 1004 | 4595 | 169 | unavailable |
| si8-gamma | strict | off | 3 | 1/1 | 4.4772 | 14–14 | 1135 | 4764 | 110 | 1.741e-13 |
| si8-gamma | strict | off | 4 | 1/1 | 2.6788 | 14–14 | 541 | 2717 | 34 | 2.203e-13 |
| si8-gamma | strict | off | 6 | 5/5 | 2.2804 | 14–14 | 315 | 2066 | 9 | 1.208e-13 |
| si8-gamma | strict | on | 2 | 0/5 | 2.2857 | 7–9 | 1001 | 4588 | 169 | unavailable |
| si8-gamma | strict | on | 3 | 1/1 | 2.5380 | 14–14 | 1082 | 4523 | 102 | 2.451e-13 |
| si8-gamma | strict | on | 4 | 1/1 | 1.6939 | 14–14 | 552 | 2875 | 37 | 1.350e-13 |
| si8-gamma | strict | on | 6 | 5/5 | 1.3612 | 14–14 | 322 | 2032 | 9 | 2.380e-13 |

## Reference checks

| Case | Real | Tight energy Ha | Difference between two tight levels Ha |
|---|---|---:|---:|
| material-mgo-40 | off | -76.050880851711000 | 7.105e-14 |
| material-si-20 | off | -8.462317836906076 | 5.329e-15 |
| si-gamma | off | -7.295530468210913 | 8.882e-15 |
| si-k222s | off | -7.919423412469110 | 3.642e-14 |
| si16-gamma | off | -63.129244421617933 | 1.847e-13 |
| si16-gamma | on | -63.129244421617713 | 2.132e-13 |
| si32-gamma | off | -126.703998096108791 | 5.826e-13 |
| si8-gamma | off | -31.321803601031732 | 1.776e-14 |
| si8-gamma | on | -31.321803601031853 | 3.553e-14 |

## Distribution and trajectory checks

Five perturbations at multipliers 2 and 6 are paired by seed. Equal SCF counts
are reported separately from exact recorded-history equality in rows.csv.
Even equal iteration counts do not imply identical internal solver trajectories.
The 3/4 sweeps have one perturbation per condition, so do not establish a robust
default on their own. Runs were grouped by condition rather than randomized;
shared-host drift and the extra setup work limit small wall-time claims.

| Case | Profile | Real | Multiplier | Wall min/median/max s | SCF min/median/max | H vectors min/median/max |
|---|---|---|---:|---|---|---|
| si16-gamma | standard | off | 2 | 3.8418 / 4.6526 / 5.4544 | 17 / 20 / 23 | 1635 / 1855 / 2160 |
| si16-gamma | standard | off | 6 | 4.5653 / 5.0908 / 6.2471 | 20 / 23 / 27 | 1739 / 1926 / 2318 |
| si16-gamma | standard | on | 2 | 3.1180 / 3.3460 / 3.6824 | 19 / 22 / 23 | 1763 / 2034 / 2133 |
| si16-gamma | standard | on | 6 | 2.4932 / 2.7143 / 3.6941 | 15 / 17 / 25 | 1430 / 1510 / 2088 |
| si16-gamma | strict | off | 2 | 16.8286 / 17.4659 / 37.8292 | 7 / 7 / 17 | 6568 / 6757 / 14012 |
| si16-gamma | strict | off | 6 | 10.4966 / 11.9127 / 23.7541 | 14 / 14 / 14 | 3175 / 3198 / 3293 |
| si16-gamma | strict | on | 2 | 8.1736 / 13.9854 / 17.5799 | 6 / 14 / 14 | 4952 / 7272 / 9238 |
| si16-gamma | strict | on | 6 | 5.7667 / 5.8624 / 6.2995 | 12 / 12 / 13 | 2724 / 2768 / 2920 |
| si8-gamma | standard | off | 2 | 0.8945 / 1.1685 / 1.2715 | 22 / 29 / 33 | 1035 / 1344 / 1488 |
| si8-gamma | standard | off | 6 | 0.9405 / 1.0817 / 1.1546 | 21 / 26 / 27 | 1025 / 1211 / 1280 |
| si8-gamma | standard | on | 2 | 0.6796 / 0.7728 / 0.9966 | 23 / 28 / 36 | 1103 / 1251 / 1663 |
| si8-gamma | standard | on | 6 | 0.6056 / 0.6488 / 0.8025 | 21 / 22 / 29 | 980 / 1038 / 1320 |
| si8-gamma | strict | off | 2 | 2.7549 / 3.9148 / 4.2516 | 6 / 7 / 7 | 3226 / 4595 / 5011 |
| si8-gamma | strict | off | 6 | 2.2382 / 2.2804 / 2.3297 | 14 / 14 / 14 | 2016 / 2066 / 2097 |
| si8-gamma | strict | on | 2 | 2.1051 / 2.2857 / 3.7535 | 7 / 7 / 9 | 4155 / 4588 / 7437 |
| si8-gamma | strict | on | 6 | 1.2755 / 1.3612 / 1.4121 | 14 / 14 / 14 | 1995 / 2032 / 2136 |
