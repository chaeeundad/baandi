# P29B robust gain of the standard convergence profile

2026-11-07. **The Si8 standard-profile median work gain is 1.785x, not the
historical single-run 2.65x.** Its seed-paired 95% bootstrap interval is
1.650–2.004x; worst observed paired gain is 1.404x. This includes the mandatory
strict certificate work. The ratio of marginal medians is 1.764x; the more
conservative min(strict)/max(standard) envelope is 1.383x. These are distinct
statistics, not interchangeable estimates. No perturbed Si8 pair erases the
standard work advantage, but the median is clearly below the historical claim.
The historical dev08 3.0x wall gain is not a measured robust wall gain.

**The one research rule tested does not beat standard across the held-out set.**
It helps Si8 and Si16 but increases all held-out H work by 5.00%, successful-only
work by 0.75%, and whole-case held-out work by 6.31%. Keep it experimental and
OFF. No default, strict tolerance, certificate or mixing change is recommended.

All results here are local Mac/pure-Rust observations at baseline
`1825372d372aced7bd6302a131e8b859841de08f`, not dev08 acceptance or QE runs.
Raw evidence, analysis scripts, source snapshots, hashes and commands are in
[`P29B-evidence/`](../../P29B-evidence/), outside `target/`, for orchestrator
archival. The original P29A harness and artifacts/* are unchanged.

## Protocol and statistical definitions

Ten P29A seeds per case are paired: strict then standard in the same process,
with identical geometry, initialization, mixer, mesh and thread budget. Eight
benchmark geometries plus material Si/MgO make ten cases; six benchmark cases
and both materials converge, while elongated/slab fail. Perturbations are the
same mean-centered LCG relative-density noise <=1e-10. Ordinary benchmark T=1;
material T=4, PBE, Kerker-Pulay, density floor and iteration budget match P29A.
No timings are claimed; some counting runs overlapped local builds/other runs.
Both arms within each Phase A pair shared the process environment.

Runs 0/1 repeat the exact P29A two-level tight reference protocol. **All 120
reference/standard scientific histories match P29A field-for-field**, including
failures; see p29a-identity.json. Strict default trajectories are new. Two
converged references must agree within 1e-12 Ha to support energy error.
The original reference JSON is preserved with provenance; no rounded document
energy is used as an oracle.

Primary Davidson/H counts include all recorded strict certificate attempts,
including rejected attempts, plus ordinary solves/first-step redo. Setup work
is excluded. P29A counted ordinary work only; main-only numbers are retained in
the JSON, and a separate main-only gain is below. Successful SCF distributions
exclude failures; work distributions include every attempted trajectory.
Failed cases have no convergence gain or delivered-accuracy claim.

Tables use min / median / max; IQR. Quartiles interpolate at (n-1)*p, with
IQR=Q3-Q1; JSON also stores Q1 and Q3. Gain is median of ten seed-paired
strict/standard total-H ratios, with 20,000 paired bootstrap resamples,
Random(2902), percentile indices 499/19499. These ten deterministic starts are
a sensitivity ensemble, not a probability model for arbitrary physical starts.
A zero-width bootstrap interval is sample degeneracy, not universal certainty.

## Phase A: paired distributions

SCF is successful convergence count; Davidson steps and H vectors include
certificate work. Each arm has ten trajectories.

| Case | Profile | Converged | SCF min/median/max; IQR | Davidson min/median/max; IQR | H vectors min/median/max; IQR |
|---|---|---:|---|---|---|
| si-gamma | strict | 10/10 | 12 / 12 / 12; 0 | 56 / 56 / 61; 5 | 353 / 353 / 358; 5 |
| si-gamma | standard | 10/10 | 9 / 9 / 9; 0 | 27 / 27 / 27; 0 | 180 / 180 / 180; 0 |
| si-k222s | strict | 10/10 | 12 / 12 / 12; 0 | 261 / 261 / 261; 0 | 1581 / 1581 / 1581; 0 |
| si-k222s | standard | 10/10 | 9 / 9 / 9; 0 | 106 / 106 / 106; 0 | 731 / 731 / 731; 0 |
| si-skew | strict | 10/10 | 21 / 21 / 21; 0 | 97 / 97 / 97; 0 | 609 / 609 / 609; 0 |
| si-skew | standard | 10/10 | 12 / 12 / 12; 0 | 32 / 32 / 32; 0 | 224 / 224 / 224; 0 |
| si-elongated | strict | 0/10 | unavailable | 749 / 957.5 / 1051; 69.5 | 4348 / 5499 / 5736; 453.5 |
| si-elongated | standard | 0/10 | unavailable | 348 / 531.5 / 666; 73.5 | 3216 / 5165.5 / 5878; 657 |
| si-slab | strict | 0/10 | unavailable | 840 / 912 / 994; 75.25 | 6090 / 6486 / 7044; 693.75 |
| si-slab | standard | 0/10 | unavailable | 471 / 506 / 550; 35 | 5170 / 5420 / 5607; 207 |
| si8-gamma | strict | 10/10 | 14 / 14 / 14; 0 | 306 / 321 / 352; 15.5 | 1994 / 2059.5 / 2113; 60.25 |
| si8-gamma | standard | 10/10 | 20 / 25.5 / 32; 5.75 | 86 / 104.5 / 124; 17 | 958 / 1167.5 / 1442; 229.75 |
| si16-gamma | strict | 10/10 | 14 / 14 / 15; 0 | 215 / 227 / 339; 25.75 | 3135 / 3223.5 / 3966; 109.5 |
| si16-gamma | standard | 10/10 | 18 / 21.5 / 27; 3 | 75 / 101.5 / 119; 19 | 1580 / 1854.5 / 2318; 283.25 |
| si-low | strict | 10/10 | 12 / 12 / 12; 0 | 47 / 47 / 47; 0 | 306 / 306 / 306; 0 |
| si-low | standard | 10/10 | 9 / 9 / 9; 0 | 24 / 24 / 24; 0 | 176 / 176 / 176; 0 |
| material-si-20 | strict | 10/10 | 15 / 15 / 15; 0 | 2865 / 2865 / 2865; 0 | 16526 / 16526 / 16526; 0 |
| material-si-20 | standard | 10/10 | 11 / 11 / 11; 0 | 1042 / 1042 / 1042; 0 | 6974 / 6974 / 6974; 0 |
| material-mgo-40 | strict | 10/10 | 18 / 18 / 18; 0 | 1306 / 1306 / 1306; 0 | 12202 / 12202 / 12202; 0 |
| material-mgo-40 | standard | 10/10 | 13 / 13 / 13; 0 | 609 / 609 / 609; 0 | 5909 / 5909 / 5909; 0 |

| Case | Paired median gain [95% bootstrap] | Worst paired gain | Conservative observed envelope | Main-only paired median gain |
|---|---|---:|---:|---:|
| si-gamma | 1.9611 [1.9611, 1.9889] | 1.9611 | 1.9611 | 1.7944 |
| si-k222s | 2.1628 [2.1628, 2.1628] | 2.1628 | 2.1628 | 1.9836 |
| si-skew | 2.7188 [2.7188, 2.7188] | 2.7188 | 2.7188 | 2.4330 |
| si-elongated | unavailable | unavailable | unavailable | unavailable |
| si-slab | unavailable | unavailable | unavailable | unavailable |
| si8-gamma | 1.7849 [1.6503, 2.0044] | 1.4036 | 1.3828 | 1.6155 |
| si16-gamma | 1.7275 [1.5872, 1.9619] | 1.3822 | 1.3525 | 1.5731 |
| si-low | 1.7386 [1.7386, 1.7386] | 1.7386 | 1.7386 | 1.6023 |
| material-si-20 | 2.3697 [2.3697, 2.3697] | 2.3697 | 2.3697 | 2.2030 |
| material-mgo-40 | 2.0650 [2.0650, 2.0650] | 2.0650 | 2.0650 | 1.8497 |

Even excluding certificates, the Si8 median gain is below 2.65x. The dev08
single-run H and wall figures come from a different backend/build/machine
context and are historical comparisons, not same-host causal speed factors.
Standard Si8 SCF is 20/25.5/32 and H is 958/1167.5/1442, versus the supplied
unperturbed 15 SCF/772 H. Strict is also perturbed: its Si8 SCF count is 14 in
all ten runs, but total H varies 1994–2113. Using an unperturbed strict numerator
with perturbed standard denominators would be unfair.

| Case | Tight reference E (Ha) | Reference plateau difference (Ha) | Max strict absolute error (Ha) | Max standard absolute error (Ha) |
|---|---:|---:|---:|---:|
| si-gamma | -7.29553046821091300e+00 | 8.881784e-15 | 2.131628e-14 | 1.045279e-10 |
| si-k222s | -7.91942341246910964e+00 | 3.641532e-14 | 2.842171e-14 | 4.552607e-10 |
| si-skew | -7.25055567296080383e+00 | 8.881784e-16 | 1.598721e-14 | 6.180780e-10 |
| si-elongated | unavailable | unavailable | unavailable | unavailable |
| si-slab | unavailable | unavailable | unavailable | unavailable |
| si8-gamma | -3.13218036010317320e+01 | 1.776357e-14 | 2.451372e-13 | 2.245319e-09 |
| si16-gamma | -6.31292444216179334e+01 | 1.847411e-13 | 3.694112e-11 | 1.880501e-09 |
| si-low | -7.23375645282481461e+00 | 7.105427e-15 | 2.930989e-14 | 2.975309e-11 |
| material-si-20 | -8.46231783690607564e+00 | 5.329071e-15 | 3.730349e-14 | 2.593357e-10 |
| material-mgo-40 | -7.60508808517109998e+01 | 7.105427e-14 | 3.836931e-13 | 1.508766e-10 |

Maximum strict error is 3.694112e-11 Ha (Si16); maximum standard error is
2.245319e-9 Ha (Si8). A strict certificate is not a proof of agreement with the
tighter solution to 1e-12 Ha. Standard's 5e-9 Ha stopping estimate is not an
energy-error certificate. All standard results remain explicitly uncertified.

| Failed case | Profile | Attempt length min/median/max; IQR |
|---|---|---|
| si-elongated | strict | 170 / 200 / 200; 0 |
| si-elongated | standard | 108 / 200 / 200; 31.5 |
| si-slab | strict | 200 / 200 / 200; 0 |
| si-slab | standard | 200 / 200 / 200; 0 |

## Mechanism and one prospective research rule

P29A's monotone ratchet only lowers the next requested Ritz-change tolerance
when 0.1*Hartree_error/electrons beats its previous value. After a rebound,
it can remain flat while density error stalls. For example, baseline Si8 seed
1 holds 1.222751e-4 Ha at iterations 7–14 and seed 9 holds 1.419566e-4 Ha at
iterations 9–12. Repeated two-step eigensolves coexist with fluctuating density
residuals. Ritz-value change can be small while eigenvector residual remains
large; the inconsistent approximate density map can then impair Pulay progress.
This is a mechanism consistent with histories, not proof that it is the sole
cause. Elongated/slab failures with strict solves argue against that stronger
claim. The relayed pure/MKL histories are also consistent with trajectory
sensitivity; they do not establish an MKL defect.

The rule is motivated by relating inner accuracy to outer residual in
[Dembo, Eisenstat and Steihaug, Inexact Newton Methods (1982)](https://epubs.siam.org/doi/10.1137/0719025).
Our SCF is not Newton and its Ritz test is not the theorem's residual test.
For an isolated eigenvector, Rayleigh-quotient error is quadratic in its angular
error: for x=cos(theta)u+sin(theta)v, u orthogonal to v, the difference is
sin(theta)^2*(v*Hv-lambda). This suggests a squared-norm experiment, but supplies
no universal density-to-eigenvector bound. In our own atomic-unit normalization:

`tau_next = max(5e-14 Ha, min(tau_ordinary, c Ha*(rho_L2/(1 Bohr^-3/2))^2))`.

Use the latest completed outer residual, c=1. First-step cap/redo, ordinary
ratchet and floor remain. It never loosens the requested tolerance on the same
history, changes no outer tolerance, allocates no controller state, and does no
hot-loop I/O. It acts before stagnation, instead of repeating P29A's rejected
five-interval tightening guard (+49.61% work). The cap is a heuristic, not a
new certified accuracy claim.

The plan was saved before candidate runs in phase-b-plan.txt. Training is
seeds 1–5 of gamma/skew/Si8/low. Candidate c=1 reduced training aggregate H
8784→7359 (-16.22%), all 20 converged, max error 1.213515e-9 Ha. Skew and low
regressed; selection was aggregate across four cases. Per the prespecified
conditional plan, c=0.1 was not tried after c=1 passed. Only one rule/parameter
was tested; selection.txt freezes it before held-out execution. Training numbers
are selection diagnostics, never the held-out performance estimate.

## Phase B: held-out results only

All rows below use seeds 6–10, never training seeds. H marks whole-case holdouts;
P marks new perturbations of training cases. Every profile has five trajectories.
Strict counts include certificate work; research outcomes remain uncertified.

| Case/group | Profile | Converged | SCF min/median/max; IQR | Davidson min/median/max; IQR | H vectors min/median/max; IQR | Max absolute E error (Ha) |
|---|---|---:|---|---|---|---:|
| si-gamma (P) | strict | 5/5 | 12 / 12 / 12; 0 | 56 / 56 / 61; 5 | 353 / 353 / 358; 5 | 1.953993e-14 |
| si-gamma (P) | standard | 5/5 | 9 / 9 / 9; 0 | 27 / 27 / 27; 0 | 180 / 180 / 180; 0 | 1.045279e-10 |
| si-gamma (P) | forcing | 5/5 | 8 / 8 / 8; 0 | 25 / 25 / 25; 0 | 167 / 167 / 167; 0 | 1.762324e-10 |
| si-skew (P) | strict | 5/5 | 21 / 21 / 21; 0 | 97 / 97 / 97; 0 | 609 / 609 / 609; 0 | 7.105427e-15 |
| si-skew (P) | standard | 5/5 | 12 / 12 / 12; 0 | 32 / 32 / 32; 0 | 224 / 224 / 224; 0 | 6.180709e-10 |
| si-skew (P) | forcing | 5/5 | 12 / 12 / 12; 0 | 34 / 34 / 34; 0 | 229 / 229 / 229; 0 | 3.762368e-10 |
| si8-gamma (P) | strict | 5/5 | 14 / 14 / 14; 0 | 314 / 330 / 352; 19 | 1994 / 2042 / 2113; 65 | 1.918465e-13 |
| si8-gamma (P) | standard | 5/5 | 20 / 21 / 32; 5 | 86 / 97 / 117; 18 | 958 / 1000 / 1442; 214 | 1.555893e-09 |
| si8-gamma (P) | forcing | 5/5 | 11 / 11 / 15; 3 | 67 / 76 / 94; 21 | 769 / 795 / 988; 175 | 8.633592e-10 |
| si-low (P) | strict | 5/5 | 12 / 12 / 12; 0 | 47 / 47 / 47; 0 | 306 / 306 / 306; 0 | 2.930989e-14 |
| si-low (P) | standard | 5/5 | 9 / 9 / 9; 0 | 24 / 24 / 24; 0 | 176 / 176 / 176; 0 | 2.975309e-11 |
| si-low (P) | forcing | 5/5 | 9 / 9 / 9; 0 | 26 / 26 / 26; 0 | 185 / 185 / 185; 0 | 2.573852e-11 |
| si-k222s (H) | strict | 5/5 | 12 / 12 / 12; 0 | 261 / 261 / 261; 0 | 1581 / 1581 / 1581; 0 | 2.842171e-14 |
| si-k222s (H) | standard | 5/5 | 9 / 9 / 9; 0 | 106 / 106 / 106; 0 | 731 / 731 / 731; 0 | 4.552607e-10 |
| si-k222s (H) | forcing | 5/5 | 9 / 9 / 9; 0 | 107 / 107 / 107; 0 | 762 / 762 / 762; 0 | 2.177947e-10 |
| si16-gamma (H) | strict | 5/5 | 14 / 14 / 15; 0 | 221 / 237 / 339; 32 | 3135 / 3282 / 3966; 72 | 3.694112e-11 |
| si16-gamma (H) | standard | 5/5 | 18 / 20 / 23; 3 | 75 / 95 / 112; 20 | 1580 / 1752 / 2050; 235 | 1.880501e-09 |
| si16-gamma (H) | forcing | 5/5 | 11 / 12 / 12; 1 | 67 / 80 / 88; 15 | 1427 / 1645 / 1668; 217 | 4.505196e-10 |
| material-si-20 (H) | strict | 5/5 | 15 / 15 / 15; 0 | 2865 / 2865 / 2865; 0 | 16526 / 16526 / 16526; 0 | 3.730349e-14 |
| material-si-20 (H) | standard | 5/5 | 11 / 11 / 11; 0 | 1042 / 1042 / 1042; 0 | 6974 / 6974 / 6974; 0 | 2.593357e-10 |
| material-si-20 (H) | forcing | 5/5 | 11 / 11 / 11; 0 | 1077 / 1077 / 1077; 0 | 7231 / 7231 / 7231; 0 | 7.472778e-11 |
| material-mgo-40 (H) | strict | 5/5 | 18 / 18 / 18; 0 | 1306 / 1306 / 1306; 0 | 12202 / 12202 / 12202; 0 | 3.694822e-13 |
| material-mgo-40 (H) | standard | 5/5 | 13 / 13 / 13; 0 | 609 / 609 / 609; 0 | 5909 / 5909 / 5909; 0 | 1.508482e-10 |
| material-mgo-40 (H) | forcing | 5/5 | 13 / 13 / 13; 0 | 614 / 614 / 614; 0 | 6252 / 6252 / 6252; 0 | 1.723492e-10 |
| si-elongated (H) | strict | 0/5 | unavailable | 749 / 954 / 1008; 77 | 4348 / 5271 / 5713; 303 | unavailable |
| si-elongated (H) | standard | 0/5 | unavailable | 538 / 553 / 666; 23 | 5227 / 5299 / 5878; 83 | unavailable |
| si-elongated (H) | forcing | 0/5 | unavailable | 546 / 670 / 714; 100 | 4588 / 5878 / 6311; 114 | unavailable |
| si-slab (H) | strict | 0/5 | unavailable | 882 / 927 / 994; 53 | 6293 / 6776 / 7044; 628 | unavailable |
| si-slab (H) | standard | 0/5 | unavailable | 480 / 494 / 514; 7 | 5221 / 5304 / 5454; 101 | unavailable |
| si-slab (H) | forcing | 0/5 | unavailable | 614 / 644 / 661; 29 | 6053 / 6247 / 6495; 218 | unavailable |

| Case | Aggregate H change vs standard | Paired median standard/research H gain | Paired median strict/research H gain |
|---|---:|---:|---:|
| si-gamma | -7.22% | 1.0778 | 2.1138 |
| si-skew | +2.23% | 0.9782 | 2.6594 |
| si8-gamma | -23.08% | 1.2594 | 2.5718 |
| si-low | +5.11% | 0.9514 | 1.6541 |
| si-k222s | +4.24% | 0.9593 | 2.0748 |
| si16-gamma | -13.55% | 1.2032 | 2.1969 |
| material-si-20 | +3.69% | 0.9645 | 2.2854 |
| material-mgo-40 | +5.80% | 0.9451 | 1.9517 |
| si-elongated | +5.65% | unavailable | unavailable |
| si-slab | +18.00% | unavailable | unavailable |

Held-out research max energy error is 8.633592e-10 Ha; all 40 formerly successful
trajectories still converge and all ten stress trajectories still fail. Across
the 50 held-out runs, H rises 139332→146296 (+5.00%). Successful-only work rises
85631→86270 (+0.75%). Whole-case holdouts rise 130822→139076 (+6.31%), whereas
new perturbations of training cases improve 8510→7220 (-15.16%). These slices
are deliberately separate: the latter must not conceal failed generalization.
No tested rule beats the current profile on the declared held-out aggregate
criterion. Si8's held-out improvement is real within this experiment, but is
insufficient for adoption. No further parameters were tuned on these results.

## Switch, preservation and verification

```toml
[scf]
convergence_profile = "standard"
standard_density_forcing = "on" # research only; default "off"
```

On requires standard and stall guard off. Unknown values, strict+on and both
research mechanisms together reject. Option/result/checkpoint provenance carries
the setting; off fields are omitted in the wire format. Exact continuation
rejects a changed setting and reproduces uninterrupted results with it enabled.

Additional authorized plumbing: engine session.rs (actual tolerance dispatcher,
result/checkpoint propagation), checkpoint.rs (identity); CLI checkpoint_io.rs
(wire default/omission) and checkpoint_inspect.rs (inspection provenance).
No other numerical path changes. The original P29A test remains byte-identical.

Verification evidence:

- `scripts/check.sh` (all workspace crates; no literal `all` argument) exit 0:
  build, test, fmt, clippy. Final log: check-final.log.
- Release standard_robust_contract and perturbation tests exit 0. Long sweep is
  ignored in debug; Phase A and held-out release sweeps each exit 0.
- Five bench_scf cases gamma/k222s/k222/k444s/si8 × T=1,2,4,8,12 × strict/standard:
  50 comparisons with `--require-identical` and explicit forcing off, all ten
  commands exit 0. Energy delta 0 Ha, identical density and unchanged SCF counts.
  Baseline binary was saved before production edits; no second target directory.
- `python3 -I scripts/make_readme_figures.py` and
  `python3 -I scripts/tests/test_readme_figures.py` exit 0 (9 tests). Generated
  README/figures are unchanged: this provisional Mac study does not create a
  new accepted performance claim or modify tests/reference outside ownership.
- First measurement attempt exit 101 was a harness-only assertion on a retained
  failed strict result. Corrected assertion conditions certification on successful
  status; failure log retained, full Phase A rerun completed. No tolerance changed.

## Reproduction and dev08 request

Run from this worktree with the external UPFs at the paths recorded in provenance.
`python3 -I P29B-evidence/run_phase_a.py` produces all 220 reference/paired runs;
`run_candidate.py train-c1` and `run_candidate.py heldout-c1` reproduce the
20 training and 50 held-out runs. The harness runs via
`cargo test --release -p dft-engine --test standard_robust_trajectory sweep --
--ignored --nocapture --test-threads=1`. Environment selectors and exact commands
are recorded by the runners. Analysis: analyze.py, analyze_training.py,
analyze_heldout.py, verify_reference.py; identity.py captures before/final gates.

The orchestrator should archive P29B-evidence before cleanup and request the
same paired sweep from dev08 for pure/native and MKL/native, retaining command,
commit/binary/UPF hashes, thread affinity, all histories, failures and certificates.
For a robust wall claim, interleave both profiles on a quiet node, with the same
seeds and case/thread conditions, and bootstrap paired wall ratios separately.
Repeat the frozen research holdout without retuning. Do not infer robust 3x wall
speed from these Mac work ratios. No remote login/server changes or QE source
access occurred here. Evidence is uncommitted on task/P29B-STANDARD-ROBUST-GAIN;
artifacts/dashboard are unchanged under the ticket's explicit prohibition.

