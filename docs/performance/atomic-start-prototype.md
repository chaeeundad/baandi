# P45A atomic starting subspace prototype

**Decision: reject promotion under the frozen rule; stop before Phase B.**
The cubic Si8/Si16 distributions collapse on both arithmetic paths, but skewed
Si16 changes from a constant 18 SCFs to 11–12 SCFs: its zero spread becomes one.
It also fails the 1.5x delivered-accuracy limit in energy (2.2566x) and density
(1.5692x). Lower medians do not override these independent failures. No
production option or default change is justified by this experiment.


## Frozen decision rule (before measurement, 2026-11-08)

Phase B is justified only if, on **each** Si8, Si16 and skewed Si16 and each
applicable complex/real path, the standard-profile SCF spread (max minus min)
is at most half the paired free-electron spread, median SCF does not rise,
and median total H work including atomic initialization does not rise by more
than 5%. No previously convergent start may fail. Strict references must still
certify and delivered accuracy must stay within 1.5 times baseline. Stop on
failure; no threshold, seed, band count or tolerance retuning. A zero baseline
spread requires a zero candidate spread. Phase B requires a separate ticket.

Ten paired density-noise seeds 201–210, mean-centered LCG noise of relative
amplitude 1e-10, follow P29A. Two tight strict references per arm/path use the
P29A settings. Compare energy, occupied eigenvalues and density against the
tighter reference, with a reference energy plateau of 1e-12 Ha. Required order:
Si8 and Si16, then externally extracted skewed Si16, PseudoDojo Si and MgO,
then k222 (complex only). Standard results are not certified.

## Premise and scope

Derived from the own-words stage 1 note in
`docs/research/qe-study/qe-stage-by-stage.md` and P41A first-map observations.
No QE source is read. This tests orbitals, not starting density or mixing.
A07 used free-electron waves versus random starts; the new premise is atomic
orbitals versus free-electron waves under the P41A perturbation ensemble.
P43A's rejected buffer-fraction rule is not repeated.

Production source, parser, defaults, artifacts and UPFs remain untouched.
Only the existing target directory is used. Raw evidence belongs in
`P45A-evidence/`, outside target. Counts are local Mac observations, not dev08
performance acceptance. The relayed QE comparison perturbs atomic positions
by 1e-7 Bohr, whereas this experiment perturbs density: these are different
perturbation protocols, not identical inputs.

## Implementation fixed before ensemble runs

Test-only UPF-v2 PP_PSWFC/PP_CHI scan requires finite, full-mesh radial arrays,
angular momentum l <= 3 and nonnegative occupations. Every supplied orbital is
included, with all 2l+1 real harmonics, including zero-occupation orbitals.
Coefficients are `(-i)^l exp(-i q.tau) Y_lm(qhat) 4 pi / sqrt(volume)` times
`integral r chi(r) j_l(q r) dr`, q=G+k. The existing projector radial integration
routine and configured even-mesh rule are reused. At q=0 only l=0 survives.
Two-pass modified Gram-Schmidt uses the library's relative 1e-12 rank test;
a dropped atomic column rejects setup. If the atomic rank is below the band
count, ascending kinetic-energy plane waves are projected and appended until
full rank. No random admixture is used (no orbital RNG seed): atomic orbitals
already cover the requested bands in the initial silicon cases. Density seeds
are independent, the P29A LCG with multiplier 6364136223846793005 and increment 1.

The first-density Hamiltonian is assembled using the engine's public physics
functions and conventions. Rayleigh-Ritz uses the full atomic span; the lowest
requested number of Ritz vectors enter the existing public SCF step API.
The copied P41A loop checks every free-electron history against ScfSession,
including bitwise density residual, tolerance and energy, work and certification.
All density perturbations precede the initial rotation. The first-map pair is
seeds 201/202 before mixing and before a possible standard first-step redo.

`START` records charge the full initial H block (one application per column),
radial sample evaluations, Gram/rotation complex multiply-adds and EVD dimension.
Radial evaluation, orthonormalization and EVD are additional non-H work, not
zero-cost operations. H counts alone are a lower bound on full initialization
cost and cannot justify a positive performance decision without costing these.

Accuracy decision statistic is the maximum absolute error over the ten standard
starts for each of energy (Ha), density (Bohr^-3/2) and occupied eigenvalues
(Ha), divided by the corresponding free-electron maximum. No additive floor
is used to turn an unresolved ratio into a pass. Reference plateau must hold
within and across arms. Missing references make the decision unqualified.

Reproduction (workspace root, same target directory):

```sh
PATH=/opt/homebrew/opt/rustup/bin:$PATH cargo test --release -p dft-engine --test atomic_start_prototype
P45A_SI_UPF=/absolute/Si.pz-vbc.UPF P45A_CASES=si8-gamma,si16-gamma P45A_ARM=atomic P45A_GAMMA_REAL=off cargo test --release -p dft-engine --test atomic_start_prototype sweep -- --ignored --nocapture --test-threads=1
```

Repeat with arm `default` and Gamma path `on`; each invocation runs two tight
references and ten standard seeds per case. `DFT_NC_DIR` supplies external
material UPFs; `P45A_SKEW` supplies the externally extracted skewed input.

Before the atomic ensemble, the Gamma reduced matrix is required to have
imaginary entries below 1e-12 Ha and these roundoff entries are set to zero.
This preserves the real-function subspace under degenerate Ritz rotations and
satisfies the existing real bridge's 1e-14 conjugate-symmetry check. The same
initial rotation is used for complex and real Gamma. This is test-only physics
symmetry enforcement, not a relaxation of a solver tolerance.

## Verification evidence

`P45A-evidence/identity.py` executes five cases (gamma, k222s, k222, k444s,
si8) at threads 1/2/4/8/12, one recorded run and one `--require-identical`
comparison per case/thread. All ten commands exit 0: 25 repeat pairs have
energy delta 0 Ha and bit-identical density. These are same-source repeat
checks, not a historical baseline comparison; production source is unchanged.
Commands, input paths and exit codes are in `identity-commands.json`.

The first atomic invocation failed before numerical work because the narrow
scanner matched the trailing `l` of `label`. The scanner now requires a whole
attribute token, and its regression fixture includes `label` before `l`.
`off-atomic-parser-failure.log` and `commands-initial.json` retain the exit-101
failure; it is not counted as a numerical SCF attempt or a converged sample.
The corrected analytic/scanner tests pass (two tests, one heavy test ignored).
Full `scripts/check.sh` after correction exits 0 (build/test/fmt/clippy).
README freshness passes nine tests; no README or reference dataset was changed.

External PseudoDojo headers have whitespace inside quoted occupation values.
The scanner accepts this while enforcing the attribute-name boundary; a
synthetic fixture tests both features together. This parser-only correction
leaves the measured Si.pz-vbc coefficients unchanged.

## Measured distributions and first map

All primary groups use ten standard seeds 201–210 and two strict references
per arm. 120 primary SCFs completed: 100 standard (uncertified), 20 strict
(certified by the copied fresh-map predicate). All converged. Every free arm
also passed the bitwise native-session history checks. Reference plateaus
within each arm are at most 8.242296e-13 Ha; cross-arm tighter-reference energy
differences are at most 3.552714e-13 Ha. These certify consistency of the
reference energy plateaus, not an independent external physics oracle.

H totals include every ordinary solve, first-step redo, strict certificate
and the initial atomic H block. They exclude non-H initialization work, which
is itemized below: **these totals are lower bounds on full starting cost and
do not establish a positive full-cost gate**. No further costing is needed to
reject the candidate because the spread and accuracy gates already fail.

| Case / path | Free SCF min/median/max | Atomic SCF | Free median H | Atomic median H, including start H | First-map free / atomic density separation (Bohr^-3/2) |
|---|---|---|---:|---:|---|
| si8-gamma / complex | 20/25/29 | 13/13/13 | 1138.5 | 658 | 3.859053e-03 / 4.120564e-14 |
| si16-gamma / complex | 17/21/26 | 9/9/9 | 1863.5 | 948 | 1.775925e-03 / 3.790700e-14 |
| si16-skew / complex | 18/18/18 | 11/12/12 | 1612 | 1200 | 1.147978e-13 / 2.598881e-05 |
| si8-gamma / real | 19/25.5/31 | 13/13/13 | 1197.5 | 658 | 3.859264e-03 / 4.121071e-14 |
| si16-gamma / real | 15/20/21 | 9/9/9 | 1721.5 | 948 | 1.586649e-03 / 3.795426e-14 |

The first-map pair is seeds 201/202, before mixing. In the skewed cell the
candidate **increases** first-map sensitivity from 1.15e-13 to 2.60e-5
Bohr^-3/2. The experiment does not isolate which internal Ritz/stopping event
causes this amplification. It cannot establish universal smoothness merely
from a physically motivated starting span.

## Delivered accuracy

Each number is the maximum absolute error among ten standard starts against
that arm/path's tighter strict reference. Ratios are atomic/free, frozen
before measurement; no error floors or tolerance changes are used.

| Case / path | Free / atomic energy error (Ha) | Free / atomic density error (Bohr^-3/2) | Free / atomic occupied eigenvalue error (Ha) | Ratios E / density / eigenvalue |
|---|---|---|---|---|
| si8-gamma / complex | 2.686441e-09 / 1.348539e-10 | 6.681705e-07 / 1.137594e-07 | 1.612165e-06 / 1.889880e-06 | 0.050198/0.170255/1.172263 |
| si16-gamma / complex | 1.836774e-09 / 3.753371e-10 | 2.030149e-07 / 6.215677e-08 | 1.637955e-06 / 1.217969e-06 | 0.204346/0.306168/0.743591 |
| si16-skew / complex | 9.283383e-10 / 2.094893e-09 | 1.281516e-07 / 2.010910e-07 | 2.268904e-06 / 2.394947e-06 | 2.256605/1.569165/1.055552 |
| si8-gamma / real | 8.285760e-09 / 1.348504e-10 | 1.162133e-06 / 1.137596e-07 | 2.863250e-06 / 1.889856e-06 | 0.016275/0.097889/0.660039 |
| si16-gamma / real | 1.535312e-09 / 3.750671e-10 | 2.019462e-07 / 6.215524e-08 | 1.703355e-06 / 1.217930e-06 | 0.244294/0.307781/0.715018 |

## Initialization work and limits

Every measured atomic start has full rank (32 for Si8, 64 for both Si16 cells),
so neither random noise nor plane-wave completion is exercised. Radial sample
counts include each uncached magnitude once per orbital/species. Small-H
formation and final rotation count complex multiply-adds as m*n*k (not 2mnk
real-flop convention). Two-pass MGS projections contribute exactly
`2 * Npw * n_atomic * (n_atomic - 1)` complex multiply-adds for these full-rank
blocks, in addition to norms, divisions and scalar reduction overhead.

| Case | Npw | Initial H applications | Radial samples | Gram + rotation complex MACs | MGS projection complex MACs | Dense EVD |
|---|---:|---:|---:|---:|---:|---|
| Si8 | 2969 | 32 | 85338 | 4940416 | 5890496 | one 32 x 32 |
| Si16 | 5915 | 64 | 205156 | 37856000 | 47698560 | one 64 x 64 |
| skewed Si16 | 5959 | 64 | 2568760 | 38137600 | 48053376 | one 64 x 64 |

No exact backend EVD operation count, conversion of scalar setup work to H
applications, total-work acceptance or wall-time improvement is claimed.
The H decrease is a component result, not a zero-cost treatment of rotation.
The radial transform reuses magnitude results within each species/orbital;
the skewed cell has many more distinct magnitudes and thus more setup work.
Default atomic-span fallback to free-electron columns is implemented but not
qualified by these full-rank cases. Material/floor/multi-k variants are not
validated by the stopped study. In particular, no real-path claim is made for
plane-wave completion when an atomic span is smaller than the requested band
count.

## Stop, scope and reproducibility

The complex skew distribution was collected in full before analysis established
the violations. The queued real-skew free control had already started; it was
interrupted (runner exit 130) after one completed strict reference, and is
excluded from the 120-primary-SCF matrix. Its partial log is retained, not
reported as a failed numerical calculation or a successful ensemble.
The paired real-skew distribution, PseudoDojo Si, MgO and k222 were **not run**
after the mandatory stop. The skew input was extracted to
`/tmp/P45A-skew-input/P43A-evidence/si16-skew.toml`; no UPF or archive was copied
into this repository. No remote measurements were attempted.

Raw evidence: `P45A-evidence/`. `run.py` and `run-skew.py` reproduce the primary
runs; the latter now stops automatically after the completed paired condition
fails the frozen rule. `analyze.py` generates `analysis.json`; `verify.py`
checks all ten seeds, 120 completions, within/cross-arm reference plateaus,
negative decision and identity-command exits. `commands-initial.json` retains
the pre-fix scanner failure; `commands.json` and `skew-commands.json` record the
successful collectors, `interrupted.json` records the deliberate stop.
`provenance.json` records HEAD, binary and external-input hashes. Parser-only
changes between batches do not change Si.pz-vbc values; source hash is stored
separately. All files remain uncommitted.

Verification commands (latest logs carry the actual results):

- `scripts/check.sh` (no argument means all workspace crates): build/test/fmt/clippy.
- `cargo test --release -p dft-engine --test atomic_start_prototype`: two fast tests, heavy sweep ignored.
- `python3 -I scripts/tests/test_readme_figures.py`: nine passing freshness tests.
- `python3 P45A-evidence/identity.py`: ten successful benchmark commands, 25 identity pairs.
- `python3 P45A-evidence/analyze.py` and `python3 P45A-evidence/verify.py`: recorded evidence checks.

No artifacts/dashboard or README modification is made, consistent with the
explicit ticket ownership. No Phase B implementation or acceptance claim is
made. This negative result does not negate the cubic-cell sensitivity finding;
it rejects this specific universal promotion rule without retuning it.

Final gate log `P45A-evidence/check-final2.log`: all four gates PASS, process
exit 0. Evidence verifier exits 0. `git diff --check` exits 0; tracked production
source, artifacts and README have no diff. Python helper syntax is validated
with `ast.parse` (exit 0); an earlier `py_compile` command failed because this
Mac redirects bytecode into a sandbox-protected Library cache, not because of
a syntax error. No permission escalation or external cache write was used.
