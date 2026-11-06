# P26E real Gamma SCF integration (Mac observations)

## Frozen calibration (before any real-path SCF comparison)

On 2026-11-06, the complex path was run on the four synthetic fixtures in
`gamma_real_scf.rs` (cubic, skew, elongated, slab; one/two atoms; local-only or
three radial projectors with s/s/p channels; five bands, an odd count).
Default free-electron initialization gives exactly zero spread for seeds
2026/2027 and T=1/8. To measure independent trajectories rather than identical
starts, calibration disables free-electron initialization ONLY for these complex
calibration runs, and repeats seeds 2026/2027 at T=1/8. All tolerances remain
unchanged. Command (exit 0):

```
cargo test --release -p dft-engine --test gamma_real_scf calibrate_complex_noise -- --ignored --nocapture
```

| Shape | max energy spread (Ha) | max occupied eigenvalue spread (Ha) | density L2 spread (Bohr^-3/2) |
| --- | ---: | ---: | ---: |
| cubic | 9.992007221626409e-16 | 2.612987048955517e-10 | 4.602238199079704e-10 |
| skew | 2.220446049250313e-16 | 3.469446951953614e-16 | 1.013226413057826e-16 |
| elongated | 5.151434834260726e-14 | 1.474633470888165e-9 | 6.999604117172750e-8 |
| slab | 2.220446049250313e-16 | 2.576080737615172e-9 | 1.480297112993362e-8 |

Frozen suite thresholds: occupied eigenvalues **3e-8 Ha**, density L2
**7e-7 Bohr^-3/2** (ten times observed suite maxima, rounded upward).
The factor ten allows finite sampling of SCF trajectories, not a change to
solver/certificate tolerances. Density L2 is sqrt(volume/N * sum(delta_n^2)),
without division by electron count. Energy gate remains **1e-10 Ha** from the
design. These gates are fixed before real comparisons and will not be changed.
Synthetic SCF counts can differ (observed complex slab 13 vs 11).

## Implementation scope

R1–R4 components are connected through one engine module. Completed complex
orbitals feed existing density, audit and derivative paths. R4 always certifies
strict residuals, including intermediate solves; this extra work is counted.
Real Hamiltonian snapshots allocate once per eigensolve, outside operator loops;
R3 has no potential-update API. No performance acceptance is inferred.

The engine requires the existing workspace `dft-numerics` dependency to name
RealMatrix in the RealLinearOperator bridge. This adds no external crate or
version. Cargo.lock changes only the engine's dependency edge.

Checkpoint v1 adds an optional string option, omitted when off and defaulting
to off on read. Exact restart compares the option; path changes use WarmStart.
No format version change or coefficient repair is introduced.

## Measured comparisons and acceptance failure

All numbers below are BAANDI versus BAANDI, not QE references. Portable pure
Rust on this shared Apple Silicon Mac; baseline HEAD
`21b033795c36418e0357431db92c6ccdcd716424`, implementation uncommitted.
The default and all numerical thresholds are unchanged.

| Fixture | real energy (Ha) | absolute delta E (Ha) | max occupied eigenvalue delta (Ha) | density L2 delta (Bohr^-3/2) | SCF complex → real |
| --- | ---: | ---: | ---: | ---: | ---: |
| synthetic cubic, 1 atom | -0.623797726399250241 | 3.330669073875470e-16 | 4.532495490039423e-10 | 9.685741052919343e-10 | 10 → 10 |
| synthetic skew, 1 atom | -0.847124717547730444 | 3.330669073875470e-16 | 1.804112415015879e-16 | 1.490413458162324e-16 | 11 → 11 |
| synthetic elongated, 2 atoms, local only | -3.26376645234356966 | 3.730349362740526e-14 | 1.893469769864709e-9 | 5.571693726425282e-8 | 11 → 11 |
| synthetic slab, 2 atoms | -1.44295754190504688 | 1.332267629550188e-15 | 8.566122589037661e-9 | 1.063228310644686e-8 | 13 → 10 |
| Si primitive Gamma, 2 atoms, 8 bands | -7.29553046821088103 | 1.207922650792170e-13 | 2.255534431450101e-8 | 1.556923501195097e-7 | **10 → 12 (FAIL)** |
| Si conventional Gamma, 8 atoms, 20 bands | -31.3218036010315828 | 3.126388037344441e-13 | 2.221413531594552e-8 | 1.752960372901859e-7 | 15 → 15 |

Energy, occupied eigenvalues and density pass the frozen thresholds in every
fixture. Synthetic occupied gaps are 0.6794520945, 0.3956538578, 0.1721838564,
and 0.3247995163 Ha, respectively. Si adds the smaller-gap real-potential cases.
No coefficient comparison is made between different solver paths.

**The design section 11 realistic SCF-count gate fails:** primitive Si takes two
more iterations, exceeding the allowed one. The external-file release test
checks accuracy and all five thread budgets for BOTH materials before asserting
this gate; it exits **101**. This is an acceptance failure, not a build error.
No solver, initial-space algorithm, certification policy or threshold was changed
to hide it. R5 cannot be reported as passing every design gate. Default stays off.
Further numerical investigation belongs to the architect/solver owner; this
integration does not alter R4's numerics.

### Work accounting and the cost of mandatory certification

The following are main SCF-map counts (as in bench_scf), with separate fresh
SCF-certificate counts where available. One H application means one vector,
never an FFT or a block call.

| Fixture | main Davidson steps, complex → real | main H vectors, complex → real | fresh SCF-certificate steps, complex → real | fresh SCF-certificate H, complex → real |
| --- | ---: | ---: | ---: | ---: |
| cubic | 42 → 40 | 124 → 165 | — | — |
| skew | 50 → 50 | 161 → 211 | — | — |
| elongated | 31 → 34 | 116 → 168 | — | — |
| slab | 70 → 59 | 247 → 250 | — | — |
| Si primitive | 41 → 55 | 281 → 435 | 8 → 5 | 38 → 31 |
| Si conventional | 283 → 320 | 2063 → 2406 | 62 → 68 | 202 → 204 |

Total main steps over these six cases: 517 → 558 (+7.93%). Arithmetic mean of
per-fixture step ratios: +6.07%; both are below the design's suite-mean +10%
criterion, although Si alone regresses. Total main H vectors: 2992 → 3635
(+21.49%). Primitive main H increases by 154; conventional by 343. Including
the fresh SCF certificates: primitive 319 → 466 (+147), conventional
2265 → 2610 (+345).

R4 explicitly applies H to all returned bands for every strict certification
attempt. It does so even when the engine's `certify=false` permits recurrence
acceptance. On the real primitive trajectory, 11 maps request tolerance above
1e-8 Ha; with 8 bands this imposes **at least 88 extra certification H vectors**.
On the real conventional trajectory, 14 such maps × 20 bands imposes **at least
280 extra certification H vectors**. These are measured-history lower bounds,
not exact attribution of the whole H delta: failed certification attempts,
changed starts, changed SCF lengths and recurrence rounding also change work.
Neither a lower-bound estimate nor the whole H delta is labelled an isolated
performance gain/loss. The complex and ORC cold-start spaces differ when a
kinetic shell is cut, as anticipated in the design.

## Fallback, determinism, derivatives and restart

The release synthetic suite has **8 passing tests**. It checks default Off;
Auto versus Off complete solver-result and density equality for nonzero k,
shifted mesh, standard profile, locking Auto, spectral guard, hybrid/RMM,
incompatible warm vectors, collinear spin and FD. On errors for all those
ineligible cases. Scalar FD is labelled `Fermi-Dirac occupations`, not spin;
its internal channel marker otherwise resembles spin setup. Eligibility and
metadata plumbing are the only changes to that module.

Each of the four synthetic fixtures and the two Si fixtures passes the real
thread matrix T=1,2,4,8,12, with identical energy and density bits, SCF counts
and available work counts. CLI raw density SHA-256 checks independently confirm
the two material cases at every budget:

- Si primitive real: `c1582034b4f8363e33390b6c65d300eab1f518bb6a7961c718eb267a54ef4e28`
- Si conventional real: `8dcead01adaaf207e48403eee0b1135ea9f86e2ee7b5c53776dc00966dd4630f`

Derivative comparisons use the existing strict post-SCF APIs, tightened SCF
(density 1e-10, energy 1e-12 Ha, eigen residual 1e-11 Ha), and a 1e-9 identity
cap for each force/stress component, stricter than the existing synthetic finite
difference gates. No derivative implementation changed.

| Shape | max force delta (Ha/Bohr) | max frozen stress delta (Ha/Bohr^3) |
| --- | ---: | ---: |
| cubic | 2.852050843757392e-12 | 1.618843947781556e-14 |
| skew | 3.515165555281665e-16 | 3.469446951953614e-18 |
| elongated | 1.454047993121321e-11 | 3.548376870110559e-14 |
| slab | 1.974927166248364e-11 | 1.677564337443371e-14 |

Exact checkpoints resume bit-identically on the same path and reject option
changes with `Rejection::Options`. Explicit warm checkpoints interoperate:

| Source → target | delta E (Ha) | occupied eigenvalue delta (Ha) | density L2 delta (Bohr^-3/2) |
| --- | ---: | ---: | ---: |
| Off → Auto | 6.661338147750939e-16 | 2.684118649565193e-9 | 4.828828595539734e-9 |
| On → Off | 4.440892098500626e-16 | 5.116955192896455e-10 | 1.978845744107983e-9 |

The Off → Auto warm run correctly stays complex: its stored gauge is not ORC
compatible. This is intentional fallback, not a successful ORC conversion.
On rejects incompatible saved vectors; a compatible real checkpoint can resume
On exactly. CLI tests independently check absent-key decoding to Off, omission
of Off on write, v1 preservation, Auto/On round trips, config defaults and invalid
string/type rejection. No claim of cross-build exact restart is made.

## Reproduction and evidence provenance

External UPF (not committed): `/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF`,
SHA-256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Baseline executable `/tmp/r5-baseline-cli` SHA-256
`4229fde33aefe0eff4b382431dcf197165d546ac2f52134550d5109d50bab2b9`.
Final release CLI SHA-256
`01462f73b78b03a00f440f41e0608d7b36b284f8af83381a122c91d07e1a33e6`.

```
# Synthetic accuracy, calibrated noise, fallback, threads, derivatives and restart:
cargo test --release -p dft-engine --test gamma_real_scf -- --include-ignored --nocapture
# CLI config and optional checkpoint field:
cargo test --release -p dft-cli --test gamma_real_config
# External materials; currently exits 101 for Gamma SCF 10 -> 12:
DFT_SI_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-engine --test gamma_real_materials -- --ignored --nocapture
# Run for T=1,2,4,8,12; baseline is the pre-edit executable on the same HEAD:
python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /path/Si.pz-vbc.UPF --cases gamma,k222s,k222,k444s,si8 --threads T --repeats 1 --baseline /tmp/r5-baseline.json --require-identical
```

The real material test uses the same 30 Ry, beta=0.5 Pulay, 100-iteration inputs
as bench_scf. Synthetic matrices use 2.5 Ha and five bands. Test fixtures are
not fitted to Si. Temporary full CLI inputs/results/density artifacts are under
`/tmp/r5-real-matrix-final/`; release logs under `/tmp/r5-final-release.log` and
`/tmp/r5-material-gate.log`. These are local working artifacts; the numerical
summary and reproducible Rust tests are retained here.

Initial compilation caught import/error-conversion and branch-expression errors
(exit 101); these were corrected. An intermediate engine gate failed fmt/clippy
(exit 1), corrected without numerical changes. During CLI checking and the first
material CLI matrix, ENOSPC caused a linker/temp-file failure and CLI exit 5;
the architect freed disk. Final checks use the same worktree target, no copied
or extra target tree. These environmental failures are not numerical failures.

## Review and dev08 handoff

No dev08 or MKL execution and no new QE run were performed. Request after the
SCF-count investigation: pure/MKL real Gamma at T=1,4,8,16, equal accuracy,
energy/density/thread identity, H-vector and Davidson counts, operator profiles
with at least 30 samples, and same-node QE 7.6 `K_POINTS gamma` bundles with
input/binary/UPF hashes. Mac timings during builds are not accepted performance
numbers and are omitted here.

The architect accepted the engine→numerics dependency edge and owns the later
ownership.toml update. Additive plumbing beyond the original list: engine and
CLI finite-temperature metadata, CLI spin metadata, engine Cargo.toml and its
lockfile edge. The fallback FD check changes no FD numerical operations.
No frozen contracts, default tolerance, complex kernel or checkpoint version
were edited. No feature-status/CLI documentation generator exists for this new
option in the checked scripts; the authorized results note documents the key.
README generation and freshness are checked without inventing accepted evidence
for a gate-failing implementation. The orchestrator can transfer this entry to
the out-of-scope performance ledger/backlog and dashboard artifacts:

`R5: default-off real Gamma bridge implemented; scientific accuracy/thread/
derivative/restart and off identity pass; Si Gamma SCF 10→12 fails ±1 design
criterion; six-case main Davidson +7.93%, H +21.49%; no default adoption. Derived
from real-gamma-design / qe-gamma-real-path study note, no QE source read.`

## Final off-identity and command gates

All 25 case/thread combinations passed `--require-identical` (five cases ×
T=1,2,4,8,12), with delta energy exactly 0 Ha. A separate exact-value check
confirmed SCF counts, Davidson steps, H counts and SCF-certificate work too.
Baseline signatures, identical for every final Off budget:

| Case | Energy (Ha) | SCF count | Density SHA-256 |
| --- | ---: | ---: | --- |
| gamma | -7.2955304682107602 | 10 | `6077d62f5ddcf7d003e61b29f12e83d73e0dc332a208097da0836dfb63ec1b26` |
| k222s | -7.9194234124690865 | 10 | `9e311684370831efd1f4e5a3f3fd41479dcadfc4794b51e2db0183f05e7e9104` |
| k222 | -7.8312300571211351 | 10 | `21238c7026d96860233edfc427807389029229e7ce26162b5025d471e255bdbe` |
| k444s | -7.9259761266480915 | 10 | `767602ec4395d37bf7a9da12cf4a2a399c3fee89d7cdf4129a05b39dee127660` |
| si8 | -31.32180360103127 | 15 | `826e1046923c4674f619f1ce2dbac27906bb9c6b66ae09b02aa1687d7680766e` |

| Final command | Exit | Outcome |
| --- | ---: | --- |
| `scripts/check.sh dft-engine` | 0 | PASS |
| `scripts/check.sh dft-cli` | 0 | PASS |
| `scripts/check.sh` | 0 | PASS |
| `python3 -I scripts/make_readme_figures.py` | 0 | PASS |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | PASS |
| Release `gamma_real_scf -- --include-ignored --nocapture` | 0 | 8 tests passed |
| Release CLI `gamma_real_config` | 0 | 2 tests passed |
| Release external `gamma_real_materials -- --ignored --nocapture` | 101 | Accuracy and all thread budgets pass; Si Gamma SCF-count assertion fails |
| Off `bench_scf --require-identical`, each T=1,2,4,8,12 | 0 each | 25/25 case-budget comparisons pass |

README/figure regeneration produced no tracked changes. Overall R5 acceptance
status: **FAILED**, solely for the explicit realistic SCF-count design gate.
All implementation changes are left uncommitted for architect review.

## P26F certification follow-up — 2026-11-06

Base `cbf604799f5b3b8f728f7748292d93a57c949305`, branch
`task/P26F-REAL-DAVIDSON-CERTIFY`; changes uncommitted. These are portable
pure-Rust Mac observations, BAANDI versus BAANDI, not QE measurements.
The previous sections remain historical R5 evidence.

### Change and exact work attribution

The engine now forwards its existing `certify` flag to the additive
`RealDavidson::solve_with_band_tolerances_and_certification` entry point.
Ordinary strict returns with false use recurrence residuals and report
`residuals_certified=false`. True requests perform a fresh bands-wide H
application before returning; a failed certificate still corrects/restarts.
As in complex Davidson, iteration-limit returns always check true residuals,
even with false. Existing convenience methods still request certification.
No tolerance, capacity, start, correction selection, complex Davidson, standard
profile, or default changes. The new operator-count test independently observes
one initial 3-column block with false (3 H vectors), versus two blocks with true
(6 H vectors), and two blocks at the iteration limit with false (6 H vectors).

All six real energies, density differences and SCF counts reproduce R5 exactly
at the printed precision. Skipping certificates changes H work, not trajectories.
The old material test's `extra_certificate_H_lower_bound` label is retained
unchanged as requested: its 88/280 now counts **removed** intermediate certificate
vectors, not remaining extra applications.

| Fixture | real energy (Ha) | absolute delta E (Ha) | max occupied eigenvalue delta (Ha) | density L2 delta (Bohr^-3/2) | SCF complex → real |
| --- | ---: | ---: | ---: | ---: | ---: |
| synthetic cubic, 1 atom | -0.623797726399250241 | 3.330669073875470e-16 | 4.532495490039423e-10 | 9.685741052919343e-10 | 10 → 10 |
| synthetic skew, 1 atom | -0.847124717547730444 | 3.330669073875470e-16 | 1.804112415015879e-16 | 1.490413458162324e-16 | 11 → 11 |
| synthetic elongated, 2 atoms, local only | -3.26376645234356966 | 3.730349362740526e-14 | 1.893469769864709e-9 | 5.571693726425282e-8 | 11 → 11 |
| synthetic slab, 2 atoms | -1.44295754190504688 | 1.332267629550188e-15 | 8.566122589037661e-9 | 1.063228310644686e-8 | 13 → 10 |
| Si primitive Gamma, 2 atoms, 8 bands | -7.29553046821088103 | 1.207922650792170e-13 | 2.255534431450101e-8 | 1.556923501195097e-7 | **10 → 12 (FAIL)** |
| Si conventional Gamma, 8 atoms, 20 bands | -31.3218036010315828 | 3.126388037344441e-13 | 2.221413531594552e-8 | 1.752960372901859e-7 | 15 → 15 |

| Fixture | main Davidson steps, complex → real | main H vectors, complex → real | fresh SCF-certificate steps, complex → real | fresh SCF-certificate H, complex → real |
| --- | ---: | ---: | ---: | ---: |
| cubic | 42 → 40 | 124 → 120 | — | — |
| skew | 50 → 50 | 161 → 161 | — | — |
| elongated | 31 → 34 | 116 → 118 | — | — |
| slab | 70 → 59 | 247 → 205 | — | — |
| Si primitive | 41 → 55 | 281 → 347 | 8 → 5 | 38 → 31 |
| Si conventional | 283 → 320 | 2063 → 2126 | 62 → 68 | 202 → 204 |

Main H total is 2992 → 3077 (+2.8409%; arithmetic mean of per-fixture ratios
+1.3393%), compared with R5's 3635 real H vectors. Exactly 558 applications
were removed: 45/50/50/45/88/280 by row. Main steps remain 517 → 558
(+7.93%, mean per-fixture +6.07%). Primitive including SCF certificates is
319 → 378 H vectors; conventional 2265 → 2330. These are work counts, not
wall-time speedups. The realistic primitive SCF-count gate still fails.

### Primitive count-gap diagnosis (no start changes applied)

The new ignored external-file diagnostic uses the original production runs
and four separate fixed-H first-map probes. The first-map probes use uniform
input density and the actual first SCF tolerance, 0.01 Ha. The public `step`
API certifies, so those probes count 47 H vectors versus production's 39;
all four require 5 Davidson steps.

The 749-plane-wave primitive basis cuts a kinetic shell at
0.562410736087223029 Ha: one direction below it, eight in the shell, seven
selected. The squared Frobenius distance of the native complex and
Hermitian-completed ORC start projectors is 1.00000000000000178 (dimensionless).
The first output-density difference is 9.787095689877689e-4 Bohr^-3/2,
although the first Ritz spectra differ only at roundoff. Feeding the same
ORC start to both solvers **in this diagnostic only** reduces the first-map
density difference to 2.526622377017739e-15 Bohr^-3/2. This isolates the
cut-shell orientation as the initial trajectory difference, rather than an
H implementation or certification mismatch. It does not prove that a new
start policy would remove exactly two SCF iterations on every fixture.

| SCF iteration | density residual complex → real (engine normalized atomic-unit norm) | steps complex → real | main H complex → real |
| --- | ---: | ---: | ---: |
| 1 | 6.826656818e-2 → 6.826667780e-2 | 5 → 5 | 39 → 39 |
| 4 | 1.876849334e-3 → 1.871749852e-3 | 4 → 4 | 26 → 26 |
| 5 | 1.382313393e-4 → 1.539525457e-4 | 4 → 5 | 25 → 29 |
| 8 | 2.344661869e-6 → 1.735723700e-5 | 3 → 4 | 21 → 28 |
| 9 | 8.465570655e-8 → 1.749445711e-6 | 4 → 4 | 20 → 27 |
| 10 | 4.242642576e-8 → 1.946414001e-7 | 4 → 4 | 30 → 21 |
| 11 | — → 5.273904051e-8 | — → 4 | — → 23 |
| 12 | — → 6.294606848e-9 | — → 7 | — → 32 |

At step 10 the real density residual still exceeds 1e-7. Step 11 requests
1.946414001e-8 Ha, above the final 1e-8 eigen-tolerance gate; step 12 requests
5.273904051e-9 Ha and certifies. The existing density-dependent schedule
therefore explains how the divergent trajectory adds the two final iterations.

Correction-order inspection found **the same rule**, not a different sort:
both select unconverged roots in ascending band order, up to bands columns
(the engine sets complex block_size=bands). Roundoff and different subspaces
can still rotate/order degenerate Ritz vectors. Buffer handling is also equal:
primitive has 5 strict roots (4 occupied plus first empty) and 3 buffers;
strict threshold is the requested tolerance, buffers use max(tolerance,1e-5).
All eight first-map residuals lie between 0.0040440062 and 0.0088339407 Ha,
so all meet the 0.01 Ha threshold in both paths. The existing fixed-H R4 buffer
oracle still gives 8/8 steps and 31/31 H, and the iteration-limit buffer and
failed-certificate tests pass. No evidence supports tuning buffer thresholds
or reordering corrections to fix this gap.

**Proposal only:** separately evaluate a common real-structured, shell-complete
starting subspace (including its added width/work) across the general fixture
suite, with an independent spectral oracle. This requires a new design decision;
neither production start vectors nor band/subspace counts were changed here.
No full SCF with a changed start was used to claim acceptance.

### Verification, provenance and remaining handoff

Commands (PATH includes `/opt/homebrew/opt/rustup/bin`; external UPF path and
SHA-256 are the same as R5):

| Command | Exit | Result |
| --- | ---: | --- |
| `scripts/check.sh` (workspace, no crate argument) | 0 | build/test/fmt/clippy PASS |
| `cargo test --release -p dft-solvers --test real_davidson -- --include-ignored --nocapture` | 0 | 8 passed, including unchanged R4 tests |
| `cargo test --release -p dft-engine --test gamma_real_scf -- --include-ignored --nocapture` | 0 | 8 passed; calibration, derivatives, restart, fallback unchanged |
| `DFT_SI_UPF=... cargo test --release -p dft-engine --test gamma_real_materials --test gamma_real_scf -- --include-ignored --nocapture` | 101 | material accuracy and T=1,2,4,8,12 pass for both; unchanged primitive SCF-count assertion fails, stopping before synthetic binary (run separately above) |
| `DFT_SI_UPF=... cargo test --release -p dft-engine --test gamma_real_diagnosis -- --ignored --nocapture` | 0 | 1 diagnostic passed |
| `cargo build --release -p dft-cli` | 0 | portable build |
| `python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf ... --cases gamma,k222s,k222,k444s,si8 --threads T --repeats 1 --baseline /tmp/r5-baseline.json --require-identical` | 0 each | T=1,2,4,8,12; 25/25 identical, delta E exactly 0 Ha |
| `python3 -I scripts/make_readme_figures.py` | 0 | no tracked generated changes |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | 9 tests passed |

All six real fixtures retain energy/density bit identity and work counts across
T=1,2,4,8,12. No MKL/dev08 run, new QE reference, or timing acceptance.
Off comparisons reuse the preserved R5 baseline, not a newly rebuilt same-HEAD
binary. The unchanged complex code and the previous R5 identity bridge support
this historical comparison; no new baseline was fabricated. Release CLI SHA-256:
`fb5311290e2dc383f1b3837abf0093b1ba05b7bfb164ac9bf3e37ca93db433b0`.
Raw local logs: `/tmp/p26f-{r4,r5,synthetic,diagnosis,check-final}.log`,
`/tmp/p26f-off-t{1,2,4,8,12}.log`. An initial diagnostic compile failed (101)
because Hartree has no FromStr implementation; corrected by adding a test-only
engine dependency on the existing dft-units crate (and its lockfile edge).
This is additive test plumbing outside OWNED_FILES, no new external dependency.

README regeneration is unchanged because no new accepted reference dataset is
claimed. The out-of-scope ledger/backlog/dashboard/session registry remain for
the orchestrator; this worktree has no P26F session entry or delegate report,
and no dashboard upload or recurring session timer was established here.
Ledger handoff: `P26F: certify semantics aligned, 558 redundant H vectors removed;
six-case H +2.84% vs complex, steps +7.93%; primitive 10→12 SCF still fails;
cut-shell first-density difference isolated with common-start diagnostic;
no start/default/tolerance changes, no performance adoption.`
Derived from real-gamma-design (and its own-words study-note inputs); no QE
source read. Existing dev08 handoff remains necessary after the start-policy
investigation: pure/MKL T=1,4,8,16, frozen accuracy and density identity, per-map
H/steps, and caller-resolved profiles with at least 30 samples; equal-core QE
7.6 real-Gamma bundles with input/binary/UPF hashes. Overall design acceptance
remains **FAILED**, while the requested certification fix and diagnosis are
complete and left uncommitted for review.

## Architect wall-time observation after P26F (2026-11-06, Mac, release, pure Rust, 1 thread, median of 3)

Command: `python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf <Si.pz-vbc.UPF> --cases gamma,si8 --repeats 3 --threads 1 --extra-toml 'gamma_real = "off"|"on"'`.

| Case | off wall (s) | on wall (s) | ratio | SCF iterations off / on |
| --- | ---: | ---: | ---: | ---: |
| Si primitive Gamma (2 atoms) | 0.109 | 0.098 | 1.11x | 10 / 12 |
| si8 Gamma (8 atoms, 20 bands) | 2.546 | 1.665 | 1.53x | 15 / 15 |

Mac observation only (shared machine, 3 samples: directional, below the 30-sample rule). The design estimated about 2.1x; the measured 1.53x on si8 is lower because density, mixing, Hartree and XC still use the complex path and the real Davidson algebra is not yet tuned. dev08 acceptance (ticket R0) is still pending.

## P26G standard-profile engine wiring — 2026-11-06

Base `9da560861a50807768844b9d3f800d5cc3c4f9e7`, branch
`task/P26G-REAL-GAMMA-R6`, uncommitted. Pure Rust, Apple Silicon Mac.
**R6 acceptance FAILED:** four of six total-energy differences exceed the
unchanged 1e-10 Ha cap. The implementation remains explicitly opt-in; no
threshold, ratchet, starts, band policy, solver numerics or default was changed.

### Change and scope

The only production edit is in `crates/dft-engine/src/gamma_real.rs`: remove
standard from the exclusion list and select the existing R4 Ritz-change API
when `options.standard` is present. Both serial and pooled callers already
forward the same tolerance and band policy. The CLI key is
`[scf] convergence_profile = "standard"`, combined with
`gamma_real = "on"` or `"auto"`. Strict mode retains P26F certification policy.
No edits outside OWNED_FILES and no extra target directory were needed.
Derived from real-gamma-design sections 8/9/11 and its own-words study inputs;
no QE source was accessed. This connects P17A's separately labelled contract;
it does not revive A15's rejected strict-refinement optimization.

### Six-fixture standard comparison

These are BAANDI real versus BAANDI complex with the **same default standard
profile**, not QE references and not certified results. Fixture geometry,
pseudopotentials, bands and all other settings are the R5 fixtures. Eigenvalue
and density differences are observations; the R5 calibrated strict caps do not
become standard-profile acceptance thresholds. Energy cap remains 1e-10 Ha.

| Fixture | real energy (Ha) | absolute delta E (Ha) | max occupied eigenvalue delta (Ha) | density L2 delta (Bohr^-3/2) | SCF complex → real |
| --- | ---: | ---: | ---: | ---: | ---: |
| synthetic cubic, 1 atom | -0.623797726280578613 | 5.411060488569319e-11 | 9.502641645786625e-8 | 2.282897735016452e-6 | 6 → 6 |
| synthetic skew, 1 atom | -0.847124717309507225 | 1.110223024625157e-16 | 1.415534356397075e-15 | 1.634654385369390e-14 | 8 → 8 |
| synthetic elongated, 2 atoms, local only | -3.26376645205618132 | **1.131787996655476e-10 FAIL** | 1.466575613534582e-6 | 1.636569729828954e-6 | 6 → 6 |
| synthetic slab, 2 atoms | -1.44295753589583353 | **1.410557448622285e-9 FAIL** | 2.144911592671939e-6 | 6.900536771091753e-6 | 17 → 17 |
| Si primitive Gamma, 2 atoms, 8 bands | -7.29553046797730431 | **1.255404669109339e-10 FAIL** | 7.542224024215116e-7 | 3.908261978832460e-6 | 9 → 9 |
| Si conventional Gamma, 8 atoms, 20 bands | -31.3218036004633120 | **4.324363089835970e-10 FAIL** | 1.109661570082576e-6 | 1.420452106095599e-5 | 14 → 15 |

| Fixture | main Davidson steps, complex → real | main H vectors, complex → real | fresh SCF-certificate steps, complex → real | fresh SCF-certificate H, complex → real |
| --- | ---: | ---: | ---: | ---: |
| cubic | 19 → 15 | 71 → 64 | 0 → 0 | 0 → 0 |
| skew | 24 → 24 | 94 → 94 | 0 → 0 | 0 → 0 |
| elongated | 16 → 14 | 66 → 65 | 0 → 0 | 0 → 0 |
| slab | 58 → 55 | 224 → 212 | 0 → 0 | 0 → 0 |
| Si primitive | 26 → 26 | 176 → 176 | 0 → 0 | 0 → 0 |
| Si conventional | 69 → 81 | 742 → 808 | 0 → 0 | 0 → 0 |

Totals: main steps 212 → 215; H vectors 1373 → 1419 (+3.3503%). This is
work accounting, not performance adoption. No tuning was attempted after the
failed energy comparisons. Standard stopping controls a Coulomb error metric,
not an absolute inter-trajectory total-energy bound. Different native starts
and rounding already create distinct trajectories in R5/P26F. These facts
explain why equality is not guaranteed, but do not establish a unique cause
of each R6 failure or justify relaxing its acceptance cap.

| Fixture | final standard error complex (Ha) | final standard error real (Ha) | first redo complex / real |
| --- | ---: | ---: | --- |
| cubic | 1.832934204791843e-9 | 1.839348069588690e-9 | false / false |
| skew | 4.341693744891783e-10 | 4.341693746238191e-10 | false / false |
| elongated | 2.447342378243545e-9 | 2.240860890580877e-9 | false / false |
| slab | 1.510636745840327e-9 | 1.257230960716203e-9 | false / false |
| Si primitive | 1.354298380633612e-9 | 3.509790784764162e-9 | false / false |
| Si conventional | 2.142608186477116e-9 | 8.165916938424420e-10 | false / false |

All final errors meet the frozen 5e-9 Ha standard criterion. Tests assert
finite nonnegative errors throughout, exact existing ratchet transitions,
no later redo, absent SCF certificates, `certified=false`, and
`residuals_certified=false`. A separate cubic test transfers a converged density
through the public density-transfer API to exercise **true** first redo:
complex/real first error 8.342819337594015e-18 / 2.361096376322530e-17 Ha,
7/7 steps, 28/26 H vectors. Both remain uncertified and the real redo retains
real dispatch. Native starts are unchanged.

### Isolation, fallback and derivatives

All six standard real fixtures were run at T=1,2,4,8,12 before the aggregate
energy assertions. Energy bits, density entries, complete returned eigen
solutions, SCF counts, work, ratchet thresholds and redo/error records agree
across budgets. Auto also matches On for the four synthetic fixtures.
CLI tests execute both On and Auto and independently inspect results.json
profile, uncertified flag, real-use flag, null fallback, threshold and absence
of certificate work.

Fallback tests preserve complete complex solutions for nonzero k, shifted
multi-k, locking, spectral guard, hybrid/RMM and incompatible warm vectors;
On rejects them. Nonzero k and spectral guard are additionally checked with
standard enabled. Spin and FD fallback tests remain unchanged and pass.
Existing setup validation still rejects standard+locking and standard+hybrid.
No new RMM enum or eligibility was introduced.

**Derivative requirement conflict:** design section 9, the engine force
preconditions and CLI preflight require strict convergence. Both complex and
real standard states reject force and frozen-stress evaluation with the same
`strict convergence profile required` error. Changing that would alter an
existing scientific contract beyond this plumbing ticket. No standard force or
stress values are fabricated. The existing tightened-strict completion tests
still pass all four synthetic shapes, maximum force difference
1.974927166248364e-11 Ha/Bohr and maximum frozen-stress difference
3.548376870110559e-14 Ha/Bohr^3 (each cap 1e-9 in its units).

The same-HEAD release CLI was built before edits in the only target directory,
then its binary copied to `/tmp/p26g-baseline-cli`. Baseline SHA-256:
`63bfa167e19c80f19b9544b7435ae315f0980a370a68b110c7f55f9de3b3231c`.
Final SHA-256:
`411eda7f5df70c2ef255e13daa34ac497f84a8bd144cef03a21fcd2b57fe960e`.
External UPF SHA-256 remains
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
No UPF was added to the repository.

All five benchmark cases, both profiles, T=1,2,4,8,12 are exactly identical to
their respective pre-edit baseline: **50/50** energy values, density SHA-256,
SCF counts, Davidson steps, H counts and certificate work. Strict benchmark
commands exit 0. Standard commands exit **5**, including the pre-edit baseline:
`bench_scf.py` independently hardcodes a strict si8 energy check against
-31.3218036010316 Ha with cap 1e-9 Ha. Unchanged complex standard si8 is
-31.321803600030876 Ha, which exceeds that separate cap. The comparison lines
all report zero energy delta and identical density; an independent JSON exact
comparison passes. The harness was not edited to suppress this failure.

### Reproduction and command outcomes

Use the pinned toolchain, with `/opt/homebrew/opt/rustup/bin` in PATH. The actual
all-workspace command is `scripts/check.sh` without an argument; `all` would
be interpreted as a crate name by this script.

```sh
cargo build --release -p dft-cli
# Before editing only:
cp target/release/dft-cli /tmp/p26g-baseline-cli
# For PROFILE = strict and standard, record the corresponding baseline;
# add --extra-toml 'convergence_profile = "standard"' for standard only.
python3 -I scripts/bench_scf.py --binary /tmp/p26g-baseline-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --threads 1 --repeats 1 \
  --record /tmp/p26g-base-PROFILE.json
# After editing/building: same cases/UPF/profile, T = 1,2,4,8,12,
# --binary target/release/dft-cli --baseline /tmp/p26g-base-PROFILE.json
# --require-identical --record /tmp/p26g-off-PROFILE-tT.json
cargo test --release -p dft-engine --test gamma_real_scf \
  -- --include-ignored --nocapture --test-threads=1
DFT_SI_UPF=/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  cargo test --release -p dft-engine --test gamma_real_materials standard \
  -- --include-ignored --nocapture --test-threads=1
cargo test --release -p dft-cli --test gamma_real_config
scripts/check.sh
python3 -I scripts/make_readme_figures.py
python3 -I scripts/tests/test_readme_figures.py
```

| Completed command/check | Exit | Outcome |
| --- | ---: | --- |
| Baseline and final `cargo build --release -p dft-cli` | 0 each | same target directory, binary hashes above |
| `scripts/check.sh` | 0 | build, test, fmt, clippy all PASS |
| Release `gamma_real_scf`, all tests including ignored | 101 | 10 passed; energy gate failed for elongated/slab; all thread comparisons completed |
| Release `gamma_real_materials standard`, UPF provided | 101 | both Si energy gates failed after all thread comparisons completed |
| Release `gamma_real_scf standard_first_redo` | 0 | true redo covered for both representations |
| Release CLI `gamma_real_config` | 0 | 3 passed, including On/Auto standard output |
| Strict baseline benchmark and five strict identity invocations | 0 each | 25/25 exact comparisons |
| Standard baseline benchmark and five standard identity invocations | 5 each | existing hardcoded strict si8 check; 25/25 exact comparisons |
| Independent recorded-JSON exact comparison | 0 | 50/50 energy/density/count/work identities |
| README generation and freshness test | 0 each | 9 tests; no tracked generated changes |
| `git diff --check` | 0 | clean |

An initial release material command omitted `DFT_SI_UPF` and exited 101 before
running physics. It was rerun with the external file; the resulting numerical
failures above are separate. The first release synthetic run (before adding
the true-redo/fallback-standard probes) had 9 passes and the same energy-gate
failure; the final 10-pass result supersedes that test count, not its numbers.

Local raw evidence: `/tmp/p26g-{base-build,final-build,check,cli,redo,standard,
synthetic-final,identity,timing}.log`, `/tmp/p26g-base-{strict,standard}.json`,
`/tmp/p26g-off-{strict,standard}-t{1,2,4,8,12}.json`,
`/tmp/p26g-identity-exits.json`, `/tmp/p26g-timing.json`.
The `/tmp` files are not durable repository artifacts; these appended tables
and committed test sources are the reviewable reproduction record.

README generation has no changes because no accepted reference dataset was
added. Ownership excludes shared ledger/backlog/dashboard/session files;
no R6 registry row or `.delegate/P26G-REAL-GAMMA-R6/report.json` exists in this
worktree. Dashboard upload/dataset updates and recurring session timer were
not established here. Orchestrator handoff row:
`P26G: standard real-Gamma dispatch implemented, default off; 4/6 energy gates
fail (max 1.41055745e-9 Ha); all real thread identities and 50 off identities
pass; standard derivatives intentionally reject; no solver/default tuning.`
No MKL/dev08 measurement or new QE reference was performed. Before any default
or speed acceptance, the architect must resolve the energy-gate failures and
the standard-derivative requirement conflict under a separately approved
contract. Existing dev08 equal-core pure/MKL and QE real-Gamma bundle requests
remain applicable; these Mac observations cannot replace them.

### Mac wall-time observation (R6, four combinations)

macOS 27.0.1, arm64, portable release, pure Rust, T=1. Five cyclic repetitions
of strict/off, strict/on, standard/off, standard/on, each with gamma then si8.
No builds or other task benchmarks overlapped this timing batch. Shared-host
load before/after: 21:17  2 users, load averages: 4.66 4.43 2.94 / 21:18  2 users, load averages: 3.88 4.26 2.93.
These are **Mac observations only**, not dev08 acceptance, and the failed
standard accuracy gates prevent an equal-accuracy speed claim. Process wall
includes CLI startup and output. Each condition reproduced identical energy
and density hashes in all five repetitions.

| Case | complex strict median (s) | real strict median (s) | complex standard median (s) | real standard median (s) | SCF counts in same order |
| --- | ---: | ---: | ---: | ---: | --- |
| gamma | 0.109247750 | 0.097314083 | 0.082767292 | 0.071499958 | 10 / 12 / 9 / 9 |
| si8 | 2.549217417 | 1.662199917 | 0.799915458 | 0.626848959 | 15 / 15 / 14 / 15 |

Raw process wall seconds in execution order (five per condition):

| Condition | Samples (s) |
| --- | --- |
| gamma-strict-off | 0.120082125, 0.109751833, 0.109247750, 0.108480084, 0.109208125 |
| si8-strict-off | 2.554553542, 2.553947334, 2.537447042, 2.537356791, 2.549217417 |
| gamma-strict-on | 0.097993291, 0.097314083, 0.097131541, 0.097091416, 0.097432291 |
| si8-strict-on | 1.657777500, 1.662199917, 1.659906875, 1.664081958, 1.668809167 |
| gamma-standard-off | 0.082875708, 0.082767292, 0.082451334, 0.083274417, 0.082392459 |
| si8-standard-off | 0.801972458, 0.797655708, 0.799915458, 0.802367083, 0.797678875 |
| gamma-standard-on | 0.071584167, 0.071208500, 0.071589167, 0.071499958, 0.071028125 |
| si8-standard-on | 0.626848959, 0.627413708, 0.624768833, 0.626713709, 0.627536958 |

Reproducer (run after all builds; this invokes the existing harness's
`run_case` function to retain raw timing and avoids its unrelated strict-si8
summary check; no acceptance check is redefined):

```python
import runpy, statistics
bench = runpy.run_path("scripts/bench_scf.py")
rows = {}
for repeat in range(5):
    for profile, mode in [("strict", "off"), ("strict", "on"),
                          ("standard", "off"), ("standard", "on")]:
        for case in ["gamma", "si8"]:
            extra = [f'gamma_real = "{mode}"']
            if profile == "standard":
                extra += ['convergence_profile = "standard"']
            r = bench["run_case"]("target/release/dft-cli",
                "/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF",
                case, 1, 1, extra)
            rows.setdefault((case, profile, mode), []).extend(r["wall_s"])
for key, samples in rows.items():
    print(key, statistics.median(samples), samples)
```

Executed as `python3 /tmp/p26g-timing.py`, exit 0. Full raw records and load
are in `/tmp/p26g-timing.json`; the numerical samples above remain in this
owned document. No speed ratio against QE or default adoption is claimed.

## P26G architect-approved delivered-accuracy gate — 2026-11-06

This section supersedes the preceding R6 **standard-mode acceptance criterion**
and derivative-conflict disposition, preserving all historical measurements.
The architect explicitly replaced direct real/complex standard agreement at
1e-10 Ha with delivered accuracy against the same certified complex strict
reference. Strict-mode tests and caps are unchanged. Standard force/stress
rejection equality and the pre-existing standard benchmark exit 5 are accepted
as sufficient; neither is now an R6 blocker. No production source, solver,
threshold, ratchet, starts or profile was changed in this follow-up.

**Final R6 status: FAILED, one of six fixtures.** The slab real-standard energy
error is 6.009215125502010e-9 Ha, above the approved 5e-9 Ha absolute cap.
Its complex-standard error is 4.598657676879725e-9 Ha; both satisfy the relative
comparison cap, but only complex satisfies the absolute cap. The other five
real-standard fixtures satisfy both caps. No tuning or further cap change was
attempted. The profile's Coulomb error is not itself a rigorous total-energy
error bound; this is an observed delivered-accuracy failure on this fixture.

### Shared certified reference and exact gate

For each R5 fixture the tests compute complex strict Gamma at P17A's
`(energy_tolerance Ha, density_tolerance Bohr^-3/2, eigen_tolerance Ha)` settings
`(1e-12,1e-9,1e-10)` and `(1e-13,1e-10,1e-11)`. The latter supplies **the same
reference** to both standard paths. All other fixture settings, including
native starts, mixing, cutoff, band count and thread budget T=1, are preserved.
Tests assert strict certification and final eigen-residual certification for
both reference runs, and verify their energy plateau <=1e-12 Ha. Maximum
observed plateau difference is 2.131628207280301e-14 Ha.
These are new calculations on the six R5 fixtures, not borrowed P17A energies;
the R5 synthetic slab is different from P17A's nonconverging stress slab.

With `c = abs(E_complex_standard - E_reference)` and
`r = abs(E_real_standard - E_reference)`, acceptance requires both
`r <= 5e-9 Ha` and `r <= max(3*c, 1e-10 Ha)`. The historical direct difference
is retained only as an observation. Density L2 is the unnormalized
`sqrt(volume/N * sum(delta_density^2))` in Bohr^-3/2; eigenvalue errors are
maximum absolute indexed differences in Ha, reported for occupied roots and
all retained roots separately. Buffer roots do not acquire tighter tolerances.

| Fixture | tight energy (Ha) | tighter/reference energy (Ha) | plateau delta (Ha) | reference SCF |
| --- | ---: | ---: | ---: | ---: |
| cubic | -6.23797726399252017e-1 | -6.23797726399251129e-1 | 8.88178419700125232e-16 | 14 |
| skew | -8.47124717547730777e-1 | -8.47124717547730444e-1 | 3.33066907387546962e-16 | 17 |
| elongated | -3.26376645234360652e0 | -3.26376645234360696e0 | 4.44089209850062616e-16 | 15 |
| slab | -1.44295754190504733e0 | -1.44295754190504866e0 | 1.33226762955018785e-15 | 18 |
| Si primitive | -7.29553046821090234e0 | -7.29553046821091300e0 | 1.06581410364015028e-14 | 16 |
| Si conventional | -3.13218036010317888e1 | -3.13218036010317675e1 | 2.13162820728030056e-14 | 23 |

| Fixture | complex standard energy (Ha) | complex error (Ha) | real error (Ha) | relative cap (Ha) | direct difference, observation (Ha) | both real gates |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| cubic | -6.23797726226468008e-1 | 1.72783121143993412e-10 | 1.18672516258300220e-10 | 5.18349363431980237e-10 | 5.41106048856931920e-11 | PASS |
| skew | -8.47124717309507114e-1 | 2.38223329951381402e-10 | 2.38223218929078939e-10 | 7.14669989854144205e-10 | 1.11022302462515654e-16 | PASS |
| elongated | -3.26376645194300252e0 | 4.00604438510754335e-10 | 2.87425638845206777e-10 | 1.20181331553226300e-9 | 1.13178799665547558e-10 | PASS |
| slab | -1.44295753730639098e0 | 4.59865767687972493e-9 | 6.00921512550200987e-9 | 1.37959730306391748e-8 | 1.41055744862228494e-9 | **FAIL absolute cap** |
| Si primitive | -7.29553046810284478e0 | 1.08068221038593038e-10 | 2.33608687949526939e-10 | 3.24204663115779113e-10 | 1.25540466910933901e-10 | PASS |
| Si conventional | -3.13218036000308757e1 | 1.00089181387374992e-9 | 5.68455504890152952e-10 | 3.00267544162124977e-9 | 4.32436308983596973e-10 | PASS |

Real-standard energies remain exactly those in the preceding R6 six-fixture
table. Every standard run remains uncertified; no tight certificate is attached
to a standard result.

| Fixture | occupied eigen error complex (Ha) | occupied eigen error real (Ha) | all-root eigen error complex (Ha) | all-root eigen error real (Ha) | density L2 error complex (Bohr^-3/2) | density L2 error real (Bohr^-3/2) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cubic | 5.74499871663913808e-6 | 5.84002513309700433e-6 | 1.74958661765867163e-4 | 1.71158483179834064e-4 | 2.61104724359920426e-6 | 2.31996324761337002e-6 |
| skew | 1.46710461865207087e-6 | 1.46710462006760523e-6 | 8.18133368496454239e-6 | 8.18133368296614094e-6 | 2.46201050774313884e-6 | 2.46201050959926628e-6 |
| elongated | 3.52200170783412148e-6 | 2.05542609429953949e-6 | 3.87631912834862646e-6 | 2.67780267726447363e-6 | 5.27812000032777299e-6 | 4.69981793158528814e-6 |
| slab | 9.63924698338392005e-7 | 1.18098689433354664e-6 | 1.38962525308983231e-6 | 1.18098689433354664e-6 | 8.33534301484108849e-6 | 1.09431833569317700e-5 |
| Si primitive | 1.17585140413378753e-6 | 1.93007380655529914e-6 | 1.17585140413378753e-6 | 1.93007380655529914e-6 | 2.21530625951193330e-6 | 4.71932610835889137e-6 |
| Si conventional | 1.65592824571958452e-6 | 1.04310554713271664e-6 | 1.65592824571958452e-6 | 1.04310554713271664e-6 | 1.22481987953357553e-5 | 7.72080897910040592e-6 |

### Verification and unchanged timing evidence

The updated ignored release tests retain all thread, metadata, ratchet and
work comparisons before aggregating accuracy failures. Synthetic release suite:
10 passed, 1 failed solely for slab delivered accuracy. Material standard suite:
1 passed, covering both Si cases. All six real-standard trajectories still
match exactly across T=1,2,4,8,12. The true-first-redo, standard derivative
rejection, fallback, strict synthetic equivalence, strict derivative and restart
tests pass. No numerical production edit follows the original R6 wiring.

The previously recorded 50/50 Off identity comparisons and five-sample cyclic
Mac wall table above remain applicable: the release CLI SHA-256 is unchanged
(`411eda7f5df70c2ef255e13daa34ac497f84a8bd144cef03a21fcd2b57fe960e`).
They were not rerun solely for test/document changes. In particular the si8
standard-to-standard identity still passes even though the old strict absolute
benchmark guard produces exit 5, now explicitly accepted by the architect.
No extra target directory was created. No new QE, MKL or dev08 claim is made.

Reproduction uses the same release commands above with the revised tests.
New raw logs: `/tmp/p26g-revised-synthetic.log`,
`/tmp/p26g-revised-material.log`, `/tmp/p26g-revised-material-all.log`,
`/tmp/p26g-revised-check.log`, and `/tmp/p26g-revised-accuracy.json`.
The final machine-readable report is `.delegate/P26G-REAL-GAMMA-R6/report.json`
(local ignored handoff artifact). Tests and these appended numerical tables
are the durable review record. Ledger handoff replaces the old four-failure
summary with: `P26G revised delivered-accuracy gate: 5/6 pass; slab real error
6.00921513e-9 Ha exceeds 5e-9 Ha, complex error 4.59865768e-9 Ha; all references
certified with plateau <=2.14e-14 Ha; no production/tolerance change.`

Final completed follow-up commands:

| Command | Exit | Outcome |
| --- | ---: | --- |
| `scripts/check.sh` | 0 | build/test/fmt/clippy PASS on revised tests |
| `cargo test --release -p dft-engine --test gamma_real_scf -- --include-ignored --nocapture --test-threads=1` | 101 | 10 passed; revised slab absolute-energy gate failed |
| `DFT_SI_UPF=... cargo test --release -p dft-engine --test gamma_real_materials standard -- --include-ignored --nocapture --test-threads=1` | 0 | both material standard cases pass revised gates |
| `DFT_SI_UPF=... cargo test --release -p dft-engine --test gamma_real_materials -- --include-ignored --nocapture --test-threads=1` | 101 | standard passed; unchanged R5 strict primitive SCF-count failure (10→12); strict energy/eigen/density/thread checks pass |
| `python3 -I scripts/make_readme_figures.py` | 0 | generated output unchanged |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | 9 passed |
| `git diff --check` | 0 | clean |

The strict primitive count failure is pre-existing R5/P26F evidence, separate
from the new standard slab delivered-accuracy failure. Neither is hidden by
the successful default debug gates (these full matrices are release-only).
The standard derivative requirement is resolved by the architect decision;
no outstanding approval is requested in this report. Shared dashboard/registry
ownership and the absence of a recurring reporting timer remain unchanged.

### Architect disposition of the R6 standard-profile gate (2026-11-06)

The absolute cap of 5e-9 Ha was set by the architect for the standard profile and is exceeded on the synthetic slab by the real path (6.01e-9 Ha) while the complex+standard path itself sits at 4.60e-9 Ha, 92 percent of the cap. The relative criterion (real error <= max(3 x complex error, 1e-10 Ha)) passes on all six fixtures. The slab result therefore shows that the fixture/profile combination delivers about 5e-9 Ha on either path, not a solver defect specific to the real path; the real error is 1.31 x the complex one. The absolute cap is kept as a reported observation and is NOT claimed as passed: the slab row stays marked as exceeding it. Acceptance of `gamma_real` with the standard profile as default remains pending dev08 measurements on the benchmark set. The strict-mode gates are unchanged.
