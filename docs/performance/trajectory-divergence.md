# P41A trajectory divergence: measurement protocol and evidence

Measurement only, source `d8a0168e2117071c4ef71fdaf8a5ad2cbe0dcb4b`.
No default or optimization is selected. Local Apple Silicon/pure-Rust counts
are observations, not dev08 acceptance. QE source was not read.

**Conclusion: the dominant measured channel in the large-Gamma standard runs
is density built from an under-resolved occupied subspace, already at the first
loose solve.** Later band-threshold decisions amplify and reshape the trajectory;
the first difference in total Davidson steps is not the first physical-state
separation. This is observed on both Gamma arithmetic paths and both cell sizes.

At amplitude 1e-10, standard SCF min/median/max becomes:

| Case / path | Default | All ordinary solves tight | Tight density only |
|---|---|---|---|
| Si8 complex | 20/25.5/32 | 9/9/9 | 14/14/15 |
| Si8 real | 21/27.5/30 | 9/9/9 | 13/14/15 |
| Si16 complex | 18/21.5/27 | 7/7/7 | 13/13/13 |
| Si16 real | 15/19/25 | 7/7/7 | 7/7/7 |

Tight solves remove all observed SCF-count variance in these four ensembles.
The density-only control removes 98.4%/97.4% of Si8 complex/real count variance
(13.04 to 0.21; 11.21 to 0.29 iterations squared) and all Si16 count variance.
This does **not** eliminate Davidson-work variation: tighter solves still have
different step counts. It does not establish a universal strict-profile cure:
Si8 default strict already has a constant 14 SCFs, while the density-only
control takes 17–19 complex / 18–19 real SCFs. All added H work is charged below.

The first-map diagnostic finds order-one occupied-subspace defects and output
density differences up to 0.00512 Bohr^-3/2 from input differences below
1.25e-12 Bohr^-3/2; tight outputs differ below 5.92e-13 Bohr^-3/2. No mixer has
acted at this point. That localizes the initiating sensitivity upstream of
mixing and excludes a pure occupied-space gauge rotation as its explanation.
The exact EVD/expansion operation responsible is not isolated.

Generality is deliberately bounded: k222 and PseudoDojo Si standard are already
stable at 10/10/10 and 11/11/11 SCFs, with identical work across seeds. Tight
oracles cost more on every reported standard case. No setting is recommended.
The noise-size experiment retains wide Gamma spreads at all three amplitudes,
but their medians differ; it does not prove exact amplitude independence.
QE's atomic start and solver structure remain ranked hypotheses, not matched
cross-code causal evidence. Full distributions, work counts, threshold events,
protocol limitations and source-only instrumentation follow.

## Protocol

The ignored `trajectory_divergence` integration test copies P29A's fixture and
mean-centered LCG perturbation helpers. Seeds 1–10 multiply the uniform density
by `1 + delta`. The amplitude parameter scales the centered uniform noise;
it is an upper bound, not an RMS perturbation. All arms use the same seeds.
Si8 and Si16 run both complex and forced real Gamma, strict and standard.
Generality uses unshifted `si-k222` and actual PseudoDojo PBE `material-si-20`
(shifted 4x4x4, 40 Ry, original Kerker-Pulay and density-floor settings).

| Arm | Intervention |
|---|---|
| default | Existing schedule, solver, mixing and stopping |
| tight | Every ordinary map uses a separate strict problem, residual tolerance 1e-10 Ha and certification, retaining outer profile stopping and mixing |
| density | Default map and warm orbitals; second strict 1e-10 Ha solve of exactly the same input density/potential, warm-started from loose states; only its density replaces the loose density |
| noise controls | Default at amplitude 1e-11 and 1e-9, versus 1e-10 |

The tight criterion applies to occupied plus first empty bands; higher empty
buffers retain the existing 1e-5 Ha criterion and do not contribute density.
The density arm retains the loose-map energy/eigenvalue diagnostics, including
the strict outer energy-change predicate; it is a deliberately hybrid oracle,
not a consistent new production SCF algorithm. The extra solve's energy is not
substituted. Final fixed-point certificate work is counted separately. In the tight arm
the strict final certificate also uses 1e-10 Ha; other arms retain the
ordinary final certificate tolerance. Three early tight-strict logs using
the ordinary 1e-9 Ha certificate were superseded and rerun; they are preserved
under `P41A-evidence/superseded-certificate/` and excluded from all final tables. No
accuracy or speedup claim is made for these diagnostic interventions.

A copied measurement driver implements the ordinary session sequence, first
standard redo, mixer initialization/feasibility damping, stopping and strict
certificate. Every default run at every amplitude is also run through the real
ScfSession: density-residual, requested tolerance and energy bits, ordinary
Davidson and H counts, trajectory length and convergence must agree. This is
an executable fidelity check rather than assuming that a copied loop is exact.
The copied driver is limited to these fixtures (no symmetry/spectral rescue,
finite temperature or control policies). Native default runs validate gap and
material clipping; the diagnostic driver does not independently certify those
conditions for intervention arms.

The sole production-file edit is making the existing
`ScfProblem::step_kpoints_with_certification` public and doc-hidden. Its body
and every production caller are unchanged. The existing public `step_kpoints`
unconditionally certifies: it cannot reproduce early ordinary solves exactly.
This visibility-only accessor is necessary to use the real kernels in the
measurement driver without introducing a production oracle switch.

## Observables and interpretation

ROW records: case, path, profile, arm, amplitude, seed, SCF iteration, normalized
density residual (Bohr^-3/2), Hartree error estimate (Ha), requested tolerance
(Ha), maximum residual of the strict prefix (Ha, includes first empty), their
ratio, ordinary Davidson steps, ordinary H vectors, density-oracle H vectors,
fixed-point-certificate H vectors, final unconverged-band count, loose/main
energy (Ha), and loose-versus-tight output-density distance (Bohr^-3/2).
The density-distance field is measured only in the density arm; its zero in
other arms is an unmeasured sentinel, not a claim of zero eigenstate error.
SOLVE records give per-k steps/H, strict-prefix length and unconverged count.
First-step redo work is combined into that iteration's per-k counts, matching
native history. It is not represented as two separate solves in histograms.

The active initial schedules are 1e-2 Ha for default strict and 5e-3 Ha for
standard. Strict then uses the density-L2 schedule, with its 1e-9 Ha floor;
standard ratchets with the Hartree error and a 5e-14 Ha floor. All 60 baseline
standard first records retain 5e-3 Ha: the implemented first-step redo does
not trigger in these fixtures. The first approximate occupied subspace is
therefore accepted under that initial loose solve. This is an observation,
not a recommendation to alter either schedule.

The residual/tolerance ratio is dimensionless but compares different stopping
quantities in standard: the threshold is a consecutive Ritz-value change,
not a residual norm. Residuals are the solver-reported recurrence norms unless
that solve requests true-residual certification; they are not all independent
Hamiltonian reapplications. A ratio above one is **not** an unconverged-band count.
Strict buffers use max(tol, 1e-5 Ha); standard buffers actually use
max(5 tol, 5e-6 Ha), not the strict buffer rule quoted in the question.
The first empty root is mandatory in the required prefix but also uses that
looser **empty** threshold in standard; only occupied roots use the requested
Ritz threshold. Required-prefix length and tight-threshold-prefix length must
not be confused. The required prefixes are 17 of 20 bands for Si8, 33 of 36
for Si16, and 5 of 8 for k222/PseudoDojo Si; occupied counts are 16, 32 and 4.
This also explains very large maximum-prefix overshoot ratios.

`pairs.json` contains all 45 seed pairs per group, first per-k step-count
difference, every common iteration's absolute scalar residual difference,
and its ratio to the preceding difference (null for zero denominator).
This is residual-norm separation, not the distance between two density fields.
A missing first divergence means no step-count difference on the common prefix;
it does not imply equal trajectory lengths. Histograms retain per-k counts,
not k-summed counts. All H totals include the second oracle and certificates;
main-only and oracle-only distributions are retained separately in JSON.
The largest standard prefix residual can belong to the unoccupied first-empty
root and therefore need not measure density error; it is not used alone to
attribute the density amplification.

### Per-band observation without production instrumentation

The public result lacks per-step Ritz changes/flags. A disposable source copy
(created by `P41A-evidence/trace_snapshot.py`) adds flat, preallocated stopping
and residual traces to the existing work-accounting collector on both solver
paths. It records the stopping values and required-prefix threshold metadata.
The analyzer applies the actual standard first-empty rule (above) to that
metadata before comparing flags. In raw TRACE rows that first empty root
carries the nominal prefix threshold; `threshold-events.json` carries the
correct effective empty threshold. Other roots use the raw threshold directly. All numerical
expressions and branches are unchanged; storage is drained outside the solve,
and output is written after the trajectory. The source diff is retained as
`trace-instrumentation.patch`. The original production solver files are never
edited. The temporary workspace shares this task's one `target/` and is removed
after collection; no UPFs are copied.

For each default 1e-10 group, every ordinary ROW/SOLVE/END field must match the
uninstrumented measurement. The traced driver's ordinary path also passes the
native ScfSession comparison. `threshold-events.json` records the actual
passing/failing band values at the shorter solve's last step for every first
step-count divergence at Gamma, plus the earliest internal flag difference.
Multiple bands may cross in one pair. This identifies an observed branch,
not a counterfactual proof that that band alone caused all later separation.
Multi-k collector ordering does not identify k indices; if a multi-k pair
branches, it is reported unresolved rather than assigning a band to a guessed k.

## Reproduction

External pseudopotentials remain outside the repository. No new target directory
is used; raw evidence and analysis are in `P41A-evidence/`.

```sh
P41A_SI_UPF=/absolute/Si.pz-vbc.UPF \
DFT_NC_DIR=/absolute/pseudodojo/ex P41A_CASES=si8-gamma \
P41A_GAMMA_REAL=on P41A_PROFILE=standard P41A_ARM=default \
P41A_AMPLITUDE=1e-10 cargo test --release -p dft-engine \
  --test trajectory_divergence sweep -- --ignored --nocapture --test-threads=1
python3 -u P41A-evidence/run.py
python3 P41A-evidence/analyze.py
python3 -u P41A-evidence/trace_snapshot.py
python3 P41A-evidence/analyze_trace.py
```

The batch runner stores individual exit codes and fails on a failed group;
resumption skips only logs with a successful test result. It runs ten seeds in
each of 60 groups (600 measured trajectories, plus 360 native default replays).
The test is ignored in debug. Machine timing is not acceptance evidence.

## QE mechanisms: hypotheses, not implementations

Rank the **missing/default-different reference mechanisms** using
[the own-words stage notes](../research/qe-study/qe-stage-by-stage.md),
sections 2, 4 and 5, checked against current BAANDI source:

1. **Atomic starting density/subspace (section 2): strongest remaining QE-side
   hypothesis.** BAANDI uses uniform density and a free-electron starting space.
   The measured density/subspace separation occurs before any mixing and before
   different stopping flags in the example pair. A better resolved initial
   occupied space is therefore the most directly motivated hypothesis.
   No atomic-start ablation or QE density perturbation was performed, so this
   is not a proven explanation of QE's narrower distribution.
2. **Generalized subspace and root-update/acceptance structure (section 4):
   plausible, not isolated.** BAANDI retains its orthonormal formulation,
   larger working capacity and all-root products, and requires the occupied
   plus first-empty prefix rather than QE's more permissive block acceptance.
   The band traces demonstrate threshold branching, but tight runs still
   branch in Davidson work while their outer iteration counts become constant.
   The generalized formulation or short restarts have not been ablated here;
   cheaper algebra alone does not establish stability. No measured solve
   returned an unconverged root, so the buffer-at-limit exception is not
   exercised in these ensembles.
3. **Default mixer metric/history/formulation (section 5): later amplifier,
   weaker as the initiating explanation.** The reference note describes
   Hartree-weighted difference-form Broyden, beta 0.7/history 8; these benchmarks
   use their existing Pulay/Kerker-Pulay defaults. Hartree mixing already exists
   as an opt-in in BAANDI, so it is absent as the default, not absent as code.
   First-map separation precedes mixing; P40A's case-dependent gains and strict
   regressions further rule out claiming a universal mixer-only remedy.

**Not missing in standard:** section 2's error ratchet and first-step redo,
and section 4's Ritz-change stopping, are already implemented. They are not
valid missing-mechanism explanations for standard. Default strict lacks the
standard ratchet/redo, but the first-map instability appears in both profiles.
Section 5's cheaper energy bookkeeping explains work after the state is built,
not the observed first occupied-subspace separation. All of these are
hypotheses/controls; nothing was implemented or tuned.

QE comparison uses the **relayed** dev08 perturbation distribution in
`tests/reference/perf/dev08_materials_qe_relayed.md`: si8 12/12/15 SCF and
85/93/95 h_psi calls; si16 10/10/10 SCF and 148/148/148 calls. QE perturbs atomic
positions by 1e-7 Bohr, changing the Hamiltonian, whereas this experiment only
perturbs initial density. These are not matched perturbation operators or
units of H work, and no causal claim about QE stability follows from the
cross-code ranges alone. No QE run is performed here.

## Early diagnostic finding: step-count branching is not the first separation

For Si8 complex standard seeds 1/2, the first different ordinary Davidson count
is SCF 6: 2 versus 5 steps. At Davidson step 2, occupied band 15 (zero-based)
has Ritz change 1.6544400750e-4 Ha against threshold 1.6786871818e-4 Ha in
seed 1 (passes), versus 1.7393027493e-4 Ha against 1.5181856252e-4 Ha in seed 2
(fails). The density-residual scalar difference changes from 1.6272342665e-4
at SCF 5 to 8.0448694936e-4 Bohr^-3/2 at SCF 6, a 4.94389x amplification,
then contracts by 0.120856x and grows by 12.0698x. These are local ratios,
not one exponential Lyapunov rate.

However, these runs **already differ at the first output**: their residual
norms differ by 2.5140723594e-5 Bohr^-3/2, despite equal six-step first solves.
At that first solve the per-band residual-norm differences remain at most
1.1353e-13 Ha across the six Davidson steps. The first internal stopping-flag
difference occurs at SCF 5, Davidson step 2, band 16, also later than the
first output separation. Equal eigen-residual norms or equal step counts do
not imply equal densities. Tight and density-only oracle first-output
residual differences for the same seeds are 7.0152e-15 and 4.9544e-15
Bohr^-3/2 respectively. Thus a later discrete stopping event amplifies an
already sensitive loose-state density map; it is not established as the
initiating cause. The separate first-map diagnostic below directly measures occupied-subspace
separation; the particular kernel that first selects it and spectral
completeness remain unmeasured.

## Verification scope

`scripts/check.sh` (no argument means all workspace crates) exits 0: build,
test, format and clippy with warnings denied. The ignored experiment is not
part of the debug-time numerical workload. README figure freshness tests exit 0;
README, its figures, reference datasets and artifacts/* are unchanged under
the ticket's explicit ownership restriction.

The release CLI build exits 0. `bench_scf.py --require-identical` compares a
same-source T=1 recording against T=1,2,4,8,12 on gamma, k222s, k222, k444s,
and si8: all 25 comparisons have zero energy delta and identical density,
with all five comparison commands exiting 0. This is a determinism/unchanged
production-path confirmation, not a historical-build speed comparison.
The production change is visibility only; no arithmetic changed. Benchmark
wall ratios are intentionally not interpreted because experiments/builds ran
concurrently. The full commands, exit codes, logs, source/UPF hashes and base
commit are in `P41A-evidence/`.

## First-map occupied-subspace diagnostic

A separate ignored `first_map` experiment uses the same ten seeds, both Gamma
sizes, both paths and both initial schedules. Each first loose map is paired
with a strict 1e-10 Ha solve from the same cold input. No SCF/mixing step has
occurred yet. For each seed, measure distances to seed 1's input and output,
the gap between the last occupied and first empty Ritz values, and the stable
projection defect `||(I - U_1 U_1^*) U_seed||_F` of the occupied columns.
The defect is dimensionless, invariant under phases or occupied-space unitary
rotations; density distances use Bohr^-3/2. There are ten starts per group,
with nine nontrivial comparisons against seed 1. The self-comparison is
excluded from the reported nontrivial extrema.

For Si8 complex, the loose standard output-density distances are
0.00351788–0.00512185 Bohr^-3/2, whereas tight outputs differ by only
2.55107e-13–3.90072e-13 Bohr^-3/2. Loose occupied-subspace defects span
0.900883–1.31653. The loose standard boundary gap is about 0.00169802 Ha,
versus 0.11072348615 Ha in the tight solve of the same potential. Under the
strict profile's loose first tolerance, the boundary is actually numerically
degenerate: gaps 7.08e-15–1.03e-13 Ha despite that same tight gap.
The loose Ritz eigenvalues themselves change by only roundoff-scale amounts
across seeds. These are different **occupied subspaces**, not just different
representations of one occupied subspace. The initial loose diagonalization
leaves the occupied/empty distinction poorly resolved; integer occupations
then feed a materially different density to mixing. Resolving this first
subspace with the tight oracle removes almost all seed dependence of that
first density. The data do not isolate which small-EVD/expansion operation
first selects the different approximate subspace; no such kernel attribution
or solver redesign is claimed. The tight gap is the gap returned by the tight
iterative oracle, not an independent dense spectral-completeness certificate.

Reproduce with the same environment and release command, replacing `sweep`
with `first_map`, then `python3 P41A-evidence/analyze_first.py`.
Raw field order after case/real/standard/seed is: input density distance,
loose output distance, tight output distance, loose subspace defect, tight
subspace defect, loose gap, tight gap, maximum loose strict-prefix eigenvalue
difference from seed 1, and maximum loose strict-prefix residual-norm difference.

| Case | Path | Profile | Input density distance max | Loose output distance max | Tight output distance max | Loose subspace defect max | Loose gap median (Ha) | Tight gap median (Ha) |
|---|---|---|---|---|---|---|---|---|
| si8-gamma | complex | strict | 1.245977e-12 | 5.118135e-03 | 3.900723e-13 | 1.315061e+00 | 7.358003e-14 | 1.107235e-01 |
| si8-gamma | complex | standard | 1.245977e-12 | 5.121850e-03 | 3.900723e-13 | 1.316529e+00 | 1.698021e-03 | 1.107235e-01 |
| si8-gamma | real | strict | 1.245977e-12 | 5.119444e-03 | 5.918421e-13 | 1.315402e+00 | 7.319145e-14 | 1.107235e-01 |
| si8-gamma | real | standard | 1.245977e-12 | 5.124282e-03 | 5.918421e-13 | 1.317016e+00 | 1.698021e-03 | 1.107235e-01 |
| si16-gamma | complex | strict | 8.800776e-13 | 1.922002e-03 | 1.589771e-13 | 1.355209e+00 | 3.047562e-14 | 1.107235e-01 |
| si16-gamma | complex | standard | 8.800776e-13 | 2.067696e-03 | 1.589771e-13 | 1.389839e+00 | 1.628914e-03 | 1.107235e-01 |
| si16-gamma | real | strict | 8.800776e-13 | 1.924696e-03 | 2.339132e-13 | 1.357665e+00 | 3.089196e-14 | 1.107235e-01 |
| si16-gamma | real | standard | 8.800776e-13 | 1.969828e-03 | 2.339132e-13 | 1.392539e+00 | 1.698021e-03 | 1.107235e-01 |

All density-distance columns are in Bohr^-3/2; subspace defects are dimensionless.

## Complete ensemble tables

All 600 reported trajectories converged. Each row is ten starts; all values
are min/median/max. H means vector applications in ordinary solves, density-only
oracles and every fixed-point certificate attempt; setup is excluded. The
strict-prefix maximum residual includes the first empty state. Full per-k
steps, overshoot, returned flags and per-pair amplification are in the raw logs
and JSON. No returned ordinary solve has an unconverged band (maximum 0).
The first-divergence column compares final per-k step counts on the common
trajectory prefix; it is not the first state difference.

| Case | Path | Profile | Arm | Amplitude | SCF min/median/max | Total H min/median/max | First step divergence min/median/max (pairs/45) |
|---|---|---|---|---|---|---|---|
| material-si-20 | complex | standard | default | 1.0e-10 | 11/11/11 | 6974/6974/6974 | none (0/45) |
| material-si-20 | complex | standard | default | 1.0e-11 | 11/11/11 | 6974/6974/6974 | none (0/45) |
| material-si-20 | complex | standard | default | 1.0e-9 | 11/11/11 | 6974/6974/6974 | none (0/45) |
| material-si-20 | complex | standard | density | 1.0e-10 | 10/10/10 | 40038/40038/40038 | none (0/45) |
| material-si-20 | complex | standard | tight | 1.0e-10 | 10/10/10 | 35715/35715/35715 | none (0/45) |
| material-si-20 | complex | strict | default | 1.0e-10 | 15/15/15 | 16526/16526/16526 | none (0/45) |
| material-si-20 | complex | strict | default | 1.0e-11 | 15/15/15 | 16526/16526/16526 | none (0/45) |
| material-si-20 | complex | strict | default | 1.0e-9 | 15/15/15 | 16526/16526/16526 | none (0/45) |
| material-si-20 | complex | strict | density | 1.0e-10 | 15/15/15 | 50612/50612/50612 | none (0/45) |
| material-si-20 | complex | strict | tight | 1.0e-10 | 14/14/14 | 44252/44252/44252 | none (0/45) |
| si-k222 | complex | standard | default | 1.0e-10 | 10/10/10 | 1600/1600/1600 | none (0/45) |
| si-k222 | complex | standard | default | 1.0e-11 | 10/10/10 | 1600/1600/1600 | none (0/45) |
| si-k222 | complex | standard | default | 1.0e-9 | 10/10/10 | 1600/1600/1600 | none (0/45) |
| si-k222 | complex | standard | density | 1.0e-10 | 9/9/9 | 8637/8648/8658 | none (0/45) |
| si-k222 | complex | standard | tight | 1.0e-10 | 9/9/9 | 7081/7092.5/7106 | 1/1/1 (45/45) |
| si-k222 | complex | strict | default | 1.0e-10 | 11/11/11 | 2986/2986/2986 | none (0/45) |
| si-k222 | complex | strict | default | 1.0e-11 | 11/11/11 | 2986/2986/2986 | none (0/45) |
| si-k222 | complex | strict | default | 1.0e-9 | 11/11/11 | 2986/2986/2986 | none (0/45) |
| si-k222 | complex | strict | density | 1.0e-10 | 10/10/10 | 8672/8681.5/8687 | none (0/45) |
| si-k222 | complex | strict | tight | 1.0e-10 | 10/10/10 | 8048/8060/8083 | 1/1/1 (45/45) |
| si16-gamma | complex | standard | default | 1.0e-10 | 18/21.5/27 | 1580/1854.5/2318 | 5/5/9 (45/45) |
| si16-gamma | complex | standard | default | 1.0e-11 | 20/21/24 | 1718/1849.5/2017 | 4/5/13 (45/45) |
| si16-gamma | complex | standard | default | 1.0e-9 | 18/21/28 | 1576/1840/2364 | 5/5/11 (45/45) |
| si16-gamma | complex | standard | density | 1.0e-10 | 13/13/13 | 9351/9444/9535 | 7/7/12 (29/45) |
| si16-gamma | complex | standard | tight | 1.0e-10 | 7/7/7 | 4561/4654/4671 | 1/1/2 (45/45) |
| si16-gamma | complex | strict | default | 1.0e-10 | 14/14/15 | 3135/3223.5/3966 | 2/3/8 (45/45) |
| si16-gamma | complex | strict | default | 1.0e-11 | 13/14/15 | 2943/3184/3987 | 2/3/9 (45/45) |
| si16-gamma | complex | strict | default | 1.0e-9 | 14/15/15 | 3232/3822/4150 | 2/3/8 (45/45) |
| si16-gamma | complex | strict | density | 1.0e-10 | 15/16/16 | 9984/10666/10836 | 4/4/14 (44/45) |
| si16-gamma | complex | strict | tight | 1.0e-10 | 9/9/9 | 5546/5654/5681 | 1/1/2 (45/45) |
| si16-gamma | real | standard | default | 1.0e-10 | 15/19/25 | 1412/1695.5/2086 | 3/4/5 (45/45) |
| si16-gamma | real | standard | default | 1.0e-11 | 15/20/25 | 1459/1744/2089 | 3/4/6 (45/45) |
| si16-gamma | real | standard | default | 1.0e-9 | 15/18/26 | 1417/1630/2163 | 3/4/5 (45/45) |
| si16-gamma | real | standard | density | 1.0e-10 | 7/7/7 | 5801/5834.5/5873 | none (0/45) |
| si16-gamma | real | standard | tight | 1.0e-10 | 7/7/7 | 4599/4639.5/4668 | 1/1/4 (45/45) |
| si16-gamma | real | strict | default | 1.0e-10 | 12/12/14 | 2659/2765.5/3486 | 2/2/7 (45/45) |
| si16-gamma | real | strict | default | 1.0e-11 | 12/13/13 | 2734/2865/3245 | 2/2/7 (45/45) |
| si16-gamma | real | strict | default | 1.0e-9 | 12/14/14 | 2683/3491/3918 | 2/2/7 (45/45) |
| si16-gamma | real | strict | density | 1.0e-10 | 11/11/12 | 7757/7822/7913 | 8/10/11 (42/45) |
| si16-gamma | real | strict | tight | 1.0e-10 | 9/9/9 | 5587/5655/5680 | 1/1/4 (45/45) |
| si8-gamma | complex | standard | default | 1.0e-10 | 20/25.5/32 | 958/1167.5/1442 | 4/5/9 (45/45) |
| si8-gamma | complex | standard | default | 1.0e-11 | 18/20/30 | 875/985/1344 | 4/5/9 (45/45) |
| si8-gamma | complex | standard | default | 1.0e-9 | 18/22/32 | 864/1051.5/1440 | 4/5/9 (45/45) |
| si8-gamma | complex | standard | density | 1.0e-10 | 14/14/15 | 5460/5854.5/6026 | 5/5/14 (45/45) |
| si8-gamma | complex | standard | tight | 1.0e-10 | 9/9/9 | 2767/2834.5/2846 | 1/1/2 (45/45) |
| si8-gamma | complex | strict | default | 1.0e-10 | 14/14/14 | 1994/2059.5/2113 | 2/3/6 (45/45) |
| si8-gamma | complex | strict | default | 1.0e-11 | 14/14/14 | 2011/2074/2100 | 2/3/6 (45/45) |
| si8-gamma | complex | strict | default | 1.0e-9 | 14/14/14 | 2000/2044.5/2125 | 2/3/6 (45/45) |
| si8-gamma | complex | strict | density | 1.0e-10 | 17/18/19 | 6165/6940/7260 | 4/4/14 (45/45) |
| si8-gamma | complex | strict | tight | 1.0e-10 | 9/9/9 | 2894/2963/2973 | 1/1/2 (45/45) |
| si8-gamma | real | standard | default | 1.0e-10 | 21/27.5/30 | 997/1281.5/1374 | 4/4/12 (45/45) |
| si8-gamma | real | standard | default | 1.0e-11 | 20/24/30 | 949/1121.5/1384 | 4/4/9 (45/45) |
| si8-gamma | real | standard | default | 1.0e-9 | 20/28/33 | 950/1286.5/1550 | 4/4/12 (45/45) |
| si8-gamma | real | standard | density | 1.0e-10 | 13/14/15 | 5456/5849.5/5926 | 5/5/12 (45/45) |
| si8-gamma | real | standard | tight | 1.0e-10 | 9/9/9 | 2744/2824/2847 | 1/1/2 (45/45) |
| si8-gamma | real | strict | default | 1.0e-10 | 14/14/14 | 2007/2058.5/2147 | 2/2/8 (45/45) |
| si8-gamma | real | strict | default | 1.0e-11 | 14/14/14 | 2010/2088/2151 | 2/2/8 (45/45) |
| si8-gamma | real | strict | default | 1.0e-9 | 14/14/14 | 2002/2051.5/2154 | 2/2/8 (45/45) |
| si8-gamma | real | strict | density | 1.0e-10 | 18/18.5/19 | 6907/7068.5/7402 | 4/5/12 (45/45) |
| si8-gamma | real | strict | tight | 1.0e-10 | 9/9/9 | 2872/2953/2974 | 1/1/2 (45/45) |

## Default 1e-10 stopping histograms and branch attribution

These are true per-solve histograms from the instrumented default runs,
including separate first-step redo solves if present, excluding certificates.
The >=10 column groups the tail for readability; `trace-summary.json` retains
every exact step count. Other arms' per-k histograms are in `analysis.json`.
The crossing columns count bands at the shorter solve's final step in each
pair's first differing iteration, so one pair can contribute several bands.
“Strict” here denotes membership in the required occupied-plus-first-empty
prefix; the first empty state still uses the empty threshold in standard.

| Case / path / profile | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | >=10 | Strict / buffer crossings |
|---|---|---|---|---|---|---|---|---|---|---|---|
| material-si-20 / complex / standard | 0 | 1270 | 1350 | 730 | 140 | 0 | 30 | 0 | 0 | 0 | 0 / 0 |
| material-si-20 / complex / strict | 0 | 0 | 10 | 360 | 2270 | 1340 | 670 | 150 | 0 | 0 | 0 / 0 |
| si-k222 / complex / standard | 0 | 330 | 200 | 200 | 40 | 30 | 0 | 0 | 0 | 0 | 0 / 0 |
| si-k222 / complex / strict | 0 | 0 | 0 | 300 | 410 | 70 | 70 | 30 | 0 | 0 | 0 / 0 |
| si16-gamma / complex / standard | 0 | 76 | 36 | 17 | 27 | 24 | 10 | 6 | 5 | 19 | 57 / 9 |
| si16-gamma / complex / strict | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 10 | 3 | 128 | 58 / 7 |
| si16-gamma / real / standard | 0 | 49 | 37 | 34 | 18 | 11 | 7 | 12 | 6 | 19 | 77 / 5 |
| si16-gamma / real / strict | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 10 | 4 | 109 | 59 / 30 |
| si8-gamma / complex / standard | 0 | 99 | 47 | 9 | 14 | 23 | 15 | 19 | 9 | 11 | 70 / 9 |
| si8-gamma / complex / strict | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 17 | 123 | 42 / 31 |
| si8-gamma / real / standard | 0 | 121 | 44 | 11 | 14 | 26 | 10 | 16 | 8 | 13 | 54 / 3 |
| si8-gamma / real / strict | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 19 | 120 | 33 / 29 |

| Case / path / profile | Occupied crossings | First empty crossings | Higher buffer crossings | Unresolved pairs |
|---|---|---|---|---|
| material-si-20 / complex / standard | 0 | 0 | 0 | 0 |
| material-si-20 / complex / strict | 0 | 0 | 0 | 0 |
| si-k222 / complex / standard | 0 | 0 | 0 | 0 |
| si-k222 / complex / strict | 0 | 0 | 0 | 0 |
| si16-gamma / complex / standard | 57 | 0 | 9 | 0 |
| si16-gamma / complex / strict | 42 | 16 | 7 | 0 |
| si16-gamma / real / standard | 76 | 1 | 5 | 0 |
| si16-gamma / real / strict | 24 | 35 | 30 | 0 |
| si8-gamma / complex / standard | 67 | 3 | 9 | 0 |
| si8-gamma / complex / strict | 18 | 24 | 31 | 0 |
| si8-gamma / real / standard | 54 | 0 | 3 | 0 |
| si8-gamma / real / strict | 18 | 15 | 29 | 0 |

The histograms do not support a universal “2 versus 9” two-mode mechanism.
The large-Gamma standard solves have a low-step concentration and a longer
variable tail; strict distributions differ. Many stopping branches are in
occupied roots, not only in the higher buffers. Neither the 400-step limit
nor the allowance for unconverged buffer roots is exercised as a returned
unconverged-band event in this ensemble.

All 540 default seed pairs were checked; the four multi-k/profile groups have
no per-k step-count divergence, so their collector's lack of k identity does
not leave an unresolved branch. The eight Gamma/profile groups have their
first differing-step events and internal flag comparisons in
`threshold-events.json`. Later differences are not monotone amplification;
the complete pair histories retain both growth and contraction.

The all-tight oracle collapses standard SCF-count variance, while its own
Davidson step counts still vary. The density-only oracle isolates the density
channel more closely but retains loose warm orbitals and, in strict mode,
loose energy/eigen residual stopping predicates. Its strict behavior must not
be read as a fully consistent high-accuracy SCF algorithm. This study does not
separate the effect of making only the first solve tight from making all later
solves tight, and it does not test atomic starts. Those are causal questions,
not recommendations to change a default.

## Evidence checks and remaining scope

`P41A-evidence/verify.py` checks the full 60-group matrix, seeds 1–10 and all
2700 seed pairs, all 120 traced default histories against the ordinary logs,
80 first-map pairs, successful command records, workspace gates, unchanged
production function bodies and untouched artifacts. The analyzer independently
checks every traced returned band's effective stopping predicate, including
the standard first-empty rule. The first parser attempt missed libtest's
inline first-ROW prefix and failed; the parser now normalizes that prefix and
asserts contiguous iteration numbers and END counts. No partial group is
included in final results.

Unmeasured: the specific EVD/expansion operation first selecting the different
approximate occupied subspace; independent dense spectral completeness; an
initial-tight-only ablation; a matched QE density perturbation / atomic-start
ablation; dev08 timing or acceptance. Reference behavior is relayed and uses
different perturbations and work units. Nothing here changes a production
numerical default, mixer, tolerance, initialization or hot loop. Leave the
work uncommitted for orchestrator review and archive `P41A-evidence/` before
worktree cleanup.
