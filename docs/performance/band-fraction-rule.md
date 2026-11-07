# P43A fixed occupied-band fraction: preregistered evaluation

## Frozen protocol (before any P43A SCF measurement)

Rule: extra = max(4, ceil(3 * occupied / 8)), f = 0.375 exactly.
Si8 selects 6, Si16 selects 12, Si32 selects 24. Only unpolarized,
fixed-occupation, single exact Gamma is eligible. Multi-k, spin and FD keep
integer four. Setup knows the occupied count; the proposed count is independent
of threads and arithmetic path, common to all k points. No production option
exists unless every gate passes; default integer four remains unchanged.

This reopens A19/P38A on the P41A loose-density sensitivity and P42A ensemble
premise, not on a converged multiplet-boundary claim. P42A already contradicts
universal stability at Si8 extra=6: real-path SCF 24/25/30, versus complex
10/10/11. Archived evidence is extracted only to /tmp/P43A-archive.

Evaluation order is skewed Si16 first (real then complex), elongated Si12,
PseudoDojo MgO 2x2x1, any available eight-atom synthetic projector fixture,
Si32 real, then Si8/Si16 consistency controls. Stop at the first failed
recommendation gate, retaining complete paired ensembles for that condition.
Unmeasured systems after this stopping rule are not passes. No tuning of f.
Skewed Si16 uses the existing Si16 fractional coordinates and 30 Ry cutoff,
with a fixed shear: lattice row 1 gains y=0.7 Bohr and row 3 gains x=0.3 Bohr.
Its input is retained in P43A-evidence/si16-skew.toml. The geometry was chosen
before measurement, not selected by convergence/work results.

New systems use P29A mean-centered LCG relative density noise 1e-10 with NEW
seeds 101 through 110 (paired between counts); consistency controls use 1–10.
T=1, pure Rust. Compare four extras with the rule count, without changing
mixing, cutoff, orbital seed, solver capacity policy or tolerances. Standard
profile is unchanged and uncertified. Two unperturbed strict references per
count use (density, energy, eigen) = (1e-9,1e-11,1e-10) and
(1e-10,1e-12,1e-11), in Bohr^-3/2, Ha, Ha. They must certify and their energies
must agree within 1e-12 Ha before accuracy is interpreted. Si32 alone permits
one tight reference and five standard starts, real only. Elongated Si12 is
dropped explicitly if its strict reference cannot converge. Optional extra
counts 6, 10 and one above the rule are deferred until primary pairs pass.

Report all failures, SCF min/median/max among successes, median total H-vector
applications including certificate work, and separate strict counts/work.
Accuracy is maximum absolute total-energy error (Ha), occupied eigenvalue
maximum error (Ha), and density sample RMS error (Bohr^-3), across standard
starts against the SAME default tight reference. Retain raw histories,
eigenvalues, density bytes, inputs, commands, exits and SHA-256 provenance.

## Frozen decision rule

On EVERY new system where count changes, for each measured arithmetic path:
- Median SCF must fall or remain equal.
- Median total H must be <=1.10 times default; median of per-system H ratios
  must be below one across new systems.
- SCF max-minus-min must be <=1.20 times default spread; zero must stay zero.
- Every default-convergent start must still converge.
- Each delivered accuracy error must be <=1.5 times default error; zero
  default error requires zero candidate error, with no added floor.
- Strict references must still certify; missing references cannot pass.
Si8 and Si16 must reproduce the claimed stability as consistency checks,
never as new-system evidence. Any failure rejects recommendation and stops
further scientific evaluation. No Phase B implementation after rejection.

## Verification and reproduction

Only one local target directory. No artifacts/* changes. No UPFs copied.
The test is ignored in debug; release command:

```sh
P29A_SI_UPF=../pseudo/Si.pz-vbc.UPF P43A_INPUT=P43A-evidence/si16-skew.toml \
P29A_CASES=si16-skew P29A_EXTRA_BANDS=4 P29A_GAMMA_REAL=on \
P43A_OUT=P43A-evidence/runs/on/4 cargo test --release -p dft-engine \
  --test band_fraction_scan sweep -- --ignored --nocapture --test-threads=1
```

Workspace gate is scripts/check.sh (without literal `all`, which names a
nonexistent package). Default identity matrix: gamma/k222s/k222/k444s/si8,
T=1,2,4,8,12, strict/standard, real off/on, bench_scf --require-identical.
Forced real multi-k must reject explicitly. With no production diff these are
same-build repeat comparisons, not independent implementation comparisons.
Mac counts are observations, not dev08 performance acceptance. No QE run or
QE source reading is part of this experiment.

## Result: reject at the first new system (2026-11-08)

The fixed rule **fails the preregistered per-system H-work bound** on the
complex path of the new skewed Si16 cell. Do not implement or recommend
`extra_bands = "auto"`. Integer four and all numerical sources remain unchanged.
The real path passed the individual bounds but cannot rescue a failure on
another supported path. The fraction was not retuned and the rule was not
restricted after measurement to hide the complex-path regression.

| Path | Extra | Standard SCF min/median/max | Median total H vectors | Tight strict SCF (two levels) | Strict total H vectors |
|---|---:|---|---:|---|---|
| real | 4 | 19/19/19 | 1753 | 17, 21 | 3405, 4494 |
| real | 12 | 17/17/17 | 1814.5 | 19, 22 | 4085, 4641 |
| complex | 4 | 18/18/18 | 1612 | 16, 19 | 3142, 3655 |
| complex | 12 | 17/17/17 | 1815 | 19, 22 | 4056, 4610 |

Complex-path work increases **12.5930521%**, exceeding the fixed 10% upper
bound (1773.2 H vectors). Real-path work increases 3.5082715%. Fewer SCF steps
do not imply less H work when each solve retains more bands. Both paths were
already stable in iteration count at default four; all spreads remain zero.
Equal iteration counts are not a claim of identical floating-point trajectories.

All 40 standard starts converged and remain uncertified; all eight strict
reference runs converged and certified. There were no setup or numerical
failures. The four two-level energy plateaus are 1.27898e-13, 1.13687e-13,
8.24230e-13 and 1.63425e-13 Ha in table order, all below the frozen 1e-12 Ha
limit. The two default tighter references differ across arithmetic paths by
3.12639e-13 Ha. Accuracy below uses each path's own default tight reference,
not a comparison against a different candidate fixed point or against QE.

| Path | Extra | Maximum energy error (Ha) | Maximum occupied eigenvalue error (Ha) | Maximum density RMS error (Bohr^-3) |
|---|---:|---:|---:|---:|
| real | 4 | 8.1382012e-10 | 1.5372408e-6 | 1.9139325e-7 |
| real | 12 | 4.5392312e-10 | 1.3349789e-6 | 1.4195820e-7 |
| complex | 4 | 9.2848040e-10 | 2.2689037e-6 | 1.7643616e-7 |
| complex | 12 | 4.5881876e-10 | 1.4247363e-6 | 1.2487365e-7 |

All accuracy gates pass; accuracy improvement does not override the failed
work gate. Strict work also increases, but it is reported separately: rejection
is already established by standard work, without adding a new strict-work gate.

### Stopping and coverage limits

Per the requested stop-on-regression rule, **no further scientific conditions
were measured**: elongated Si12, MgO 2x2x1, Si32, optional extra counts and
Si8/Si16 consistency reruns are unmeasured, not successful and not dropped for
nonconvergence. Si12 and MgO input files were frozen before either was run;
their geometry descriptions and hashes remain in
`P43A-evidence/additional-fixtures-frozen.json`. No claim that either converges
is made. Existing synthetic benchmark fixtures have 1, 2 or 16 atoms; none
has the requested eight atoms, so no new synthetic potential was invented.

The Si8/Si16 historical comparison in `prior-comparison.json` is an extract
from the archived P42A summary, **not a new consistency reproduction**. In
particular, Si8 real extra=6 historically remains unstable in count. The
new-system complex-path counterexample independently rejects this fixed rule.
No median-over-all-systems improvement, universal stability threshold, or
Si32 result is claimed. This is a stopped negative study, not full coverage
of the originally proposed matrix. A future different rule requires a new
preregistration and independent evaluation; there is no per-case tuning here.

### Evidence and completed checks

Preregistration: `P43A-evidence/preregistered.md`, UTC and SHA-256 in
`freeze.json`; the audit verifies it predates all four condition executions.
Source HEAD, toolchain and numerical source hashes are in
`source-provenance.json`. Every condition retains command, environment,
fixture/UPF/binary hashes, process exit 0, full scalar histories, eigenvalues,
and little-endian f64 density samples. `summary.csv` and `summary.json` are
computed by `analyze.py`. All raw evidence is in the worktree, outside target.
No UPF was copied. The archive was extracted outside the repository.

Reproduce the collected study and checks from the repository root:

```sh
cargo test --release -p dft-engine --test band_fraction_scan --no-run
python3 -I P43A-evidence/collect.py on
python3 -I P43A-evidence/collect.py off
python3 -I P43A-evidence/analyze.py
cargo build --release -p dft-cli
python3 -I P43A-evidence/identity.py
python3 -I P43A-evidence/audit_identity.py
scripts/check.sh
python3 -I scripts/make_readme_figures.py
python3 -I scripts/tests/test_readme_figures.py
python3 -I P43A-evidence/verify.py
```

Collectors refuse to overwrite retained executions; use a fresh worktree
with the retained input/protocol files for a new collection. Analyses and
audits may be rerun directly. The existing target directory is the only one.

Workspace build/test/fmt/clippy all pass (exit 0). The identity matrix passes
70 supported repeat pairs with exact energy bits, density hashes, ordinary
and certificate work; 14 groups agree across T=1,2,4,8,12. All 30 forced-real
multi-k combinations reject as expected (individual exit 3, matrix exit 0).
No production diff exists: these are same-build repeat and thread comparisons,
not a claim of two independent implementations. README generator and nine
freshness tests pass; generated README/figures have no diff. This negative
measurement study adds no production validation claim to README.

The evidence audit passes: 48 terminal runs, eight certified references,
40 successful standard starts, frozen provenance, no UPFs, unchanged numerical
sources and artifacts, and all identity/workspace gates. The only additions
outside the owned test/document/ledger are raw evidence, fixture inputs and
collector/analysis/audit plumbing under `P43A-evidence/`. No dependencies,
config parser, metadata, restart or solver code changed. Changes are uncommitted
on `task/P43A-BAND-FRACTION-RULE` for review.
