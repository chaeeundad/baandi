# P40A mixing survey (Q3 step 1)

**Completed: all five alternatives are case-dependent; none is a candidate
universal default under the prespecified rule. No defaults or production code
changed.** The complete survey contains 2,244 SCFs: six configurations × 17
case/path groups × (two tight references + ten strict + ten standard starts).
All 102 collector processes exit 0. There are 1,716 successful SCFs and 528
nonconverged attempts, all from the already failing elongated/slab cases (264
each); 48 of these are references and 480 are perturbed starts. No configuration
introduces a new failure. All 936 successful strict/reference results certify;
all 780 successful standard results remain uncertified.

F is C plus the existing Hartree metric. The frozen B–D selection score (median
of case medians over paired successful strict/standard starts and applicable
paths) is C -17.735%, B -4.873%, D -3.709%; all have zero new failures. This
selects an exploratory combination, not a default. C improves standard SCF
medians on 8/10 complex and 5/7 real cases; F does so on 8/10 and 4/7. Both
improve only 6/10 complex and 2/7 real cases in strict, and both have >10%
strict regressions. Every alternative also violates at least one available
case's 1.5x accuracy requirement; the decision does not rely solely on the two
unqualified stress cases.

F's standard paired SCF change is -28.205% [95% hierarchical bootstrap
-35.282%, -22.222%] on complex (8 cases, 80 successful pairs) and -22.222%
[-40.000%, 0.000%] on real (5 cases, 50 successful pairs). Standard total H-work
change, including failed attempts, is -20.166% [-26.786%, -10.673%] on complex
(10 cases, 100 pairs) and -12.429% [-21.455%, -0.803%] on real (7 cases, 70 pairs).
These intervals are exploratory and conditional on selecting F; they are not
selection-adjusted confidence or a wall-time speed claim.

The Si8 example explains why one attractive distribution cannot set a universal
default: F changes standard SCF medians 25.5→16 (complex) and 27.5→16 (real),
with H medians 1167.5→868.5 and 1281.5→864.5 vectors. Strict SCF rises 14→16
on both paths (+14.29%), with H 2059.5→2516.5 and 2058.5→2474.5. F's Si16
standard medians fall 21.5→13 and 19→13, but its strict real median rises 12→15.
E alone raises Si8 complex standard maximum energy error from 2.245319e-9 to
1.666601e-8 Ha. The metric and path interact with the trajectory; P17A's neutral
single-trajectory counts did not establish distributional neutrality.

All converged two-level reference pairs meet the 1e-12 Ha plateau condition;
the largest qualified tight cross-arm difference is 7.247536e-13 Ha (E MgO).
Reference-backed strict/default errors can nevertheless exceed 1e-12 Ha:
baseline Si16 complex reaches 3.694112e-11 Ha. Relative accuracy gates near
roundoff are intentionally evaluated without changing the prespecified floor.
See the full per-case tables, failed-attempt lengths in `analysis.json`, and
[raw evidence](../../P40A-evidence/). This is a Mac/pure-Rust observation at
source HEAD `49877d03a0ef8476cffc2fafd45021c6fe1d3d43`, not dev08 acceptance.

## Frozen protocol, declared before measurements (2026-11-08)

Measurement only. No production code, mixer, default, initialization algorithm,
certificate or tolerance is changed. Derived from the own-words stage 2 and 5
notes in [qe-stage-by-stage.md](../research/qe-study/qe-stage-by-stage.md), not
from QE source. A21's fixed-factor tolerance loosening and A22's converged RMM
are not repeated. P17A-5 judged the existing Hartree metric on single trajectories;
P29A/P29B demonstrated that those are an inadequate sampling unit (20–34 SCF
iterations). The changed premise is paired distributions, not a new mixer.
The Q1 a3 relay reports beta 0.5 versus 0.7 helping QE si16 (9 versus 10 SCF,
-9.5% wall), si32 (-5.8%) and k444s (-4.7%), while si8 favors 0.7.
Those are historical dev08 QE observations, not new reference calculations here.

**Decision rule, frozen before running:** a configuration is a candidate default
only if its median successful SCF count falls on at least 70% of cases, no case
rises more than 10%, there are no new failures, and delivered accuracy is not
worse than 1.5 times baseline, separately in strict and standard and separately
on complex and real paths where applicable. Otherwise report case-dependent.
Missing convergence medians or unqualified references cannot pass a gate; failed
cases remain in the case denominator. No default will change in this ticket.
Accuracy means maximum absolute total-energy error over the ten starts in each
case/profile/path, measured against the common A tight reference. Compare the
errors directly, without an added numerical floor (zero baseline requires zero
candidate error). Every successful strict result must be certified; standard is
explicitly uncertified. Reference qualification and cross-config fixed points
are reported separately, never silently treated as equal.

Ten cases: the eight P29A benchmark geometries si-gamma, si-k222s, si-skew,
si-elongated, si-slab, si8-gamma, si16-gamma, si-low, plus material-si-20 and
material-mgo-40. They cover cubic, skewed, elongated and slab-like cells,
2/4/8/16 atoms, Gamma and shifted meshes, LDA-PZ and PBE, and two distinct
pseudopotential families. All actual external UPFs are available; no synthetic
substitution is necessary. Benchmark defaults: Pulay beta 0.3, history 6, T=1,
200 iterations. Material fixtures: Kerker-Pulay beta 0.3, history 6, q0=1 Bohr^-1,
T=4, 120 iterations, density floor 1e-10, their original 40/80 Ry cutoffs and
4x4x4/3x3x3 shifted meshes. These are P29A defaults, not the beta 0.5 CLI
benchmark setting used in some historical dev08 tables.

| Arm | Change relative to its fixture defaults |
|---|---|
| A | None |
| B | beta 0.5 |
| C | beta 0.7, history 8 |
| D | beta 0.5, history 8 |
| E | `ScfOptions.hartree_mixing=true`, all other defaults |
| F | E plus a globally selected B/C/D, using the selection below |

Each arm/case/path has 22 runs: unperturbed strict references 0/1 at density
1e-9/1e-10, energy 1e-11/1e-12 Ha, eigen 1e-10/1e-11 Ha, followed by ten
seed-paired default-strict/standard runs (even/odd IDs 2–21). Perturbations are
the exact P29A mean-centered LCG multiplicative uniform-density noise <=1e-10.
Two references must converge, certify and agree within 1e-12 Ha; the tighter
A energy is the common oracle, and each candidate's own two-reference plateau
must also qualify before reporting its delivered accuracy. Cross-arm tight
energy differences are retained. Failed reference energies never become oracles.
Seven Gamma fixtures use both off and `P29A_GAMMA_REAL=on`; the three mesh
fixtures use complex only. `on` never silently falls back. `auto` is supported
for diagnostics but is not substituted for the mandatory on observations.

SCF distributions exclude failures. H-vector work includes all attempts and all
recorded fixed-point certificate attempts (including rejected certificates),
excluding setup. Failure counts and attempted lengths are retained independently.
Raw full histories, certified flags and actual real-path flags are stored outside
`target/` in `P40A-evidence/`. A test process succeeding means collection completed,
not that all SCFs converged. Missing/duplicate rows fail the analysis.

F selection is frozen as follows: after all A–E runs, rank B/C/D first by number
of new failed paired starts, then by median paired SCF percentage change across
all mutually successful case/path/profile/seed pairs (equal weight per case,
with within-case median first), then B/C/D lexical order. Run the winning
combination with E on the complete set. This adaptive F selection is exploratory;
its bootstrap does not correct for selection and cannot support an independent
confirmation claim.

Report median paired percentage change `100*(candidate/default - 1)` for SCF and
H work. Use 20,000 hierarchical paired bootstrap replicates, deterministic seed
4001: resample cases with replacement, then paired seeds within each sampled
case, then median of case medians. Report percentile 2.5/97.5% intervals for
each profile/path and arm. SCF requires success in both members; H includes
failures. Missing SCF pairs are explicitly counted; never use exhausted iteration
budgets as successful counts. The sensitivity ensemble of ten deterministic
starts is not a random population of physical starting densities.

## Reproduction

Use only the worktree's ordinary target directory. All raw files go in
`P40A-evidence/`; external UPFs must not be copied or committed.

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export P29A_SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF"
export DFT_NC_DIR="$PWD/../pseudo/pseudodojo/ex"
P40A_CONFIG=A P29A_GAMMA_REAL=off cargo test --release -p dft-engine \
  --test mixing_survey sweep -- --ignored --nocapture --test-threads=1 \
  > P40A-evidence/A-off.log 2>&1
# Repeat A–E with off and on; on explicitly skips the three mesh cases.
# After frozen selection, repeat F with P40A_F_BASE=B (or C/D as selected).
```

`P40A_CASES` selects a comma-separated subset for resumable collection;
`P40A_START_RUN`/`P40A_END_RUN` select a half-open run-ID interval within 0..22.
A full analysis requires all 22 runs in all requested groups. Builds, command
exit codes, source/input/binary hashes and analysis programs are archived with
the raw evidence. Local measurements are Mac observations; dev08 acceptance
must be requested through the orchestrator, using the same matrix and seeds.

## Verification and interpretation boundaries

`scripts/check.sh` (without a crate argument: all workspace crates) exited 0:
build, test, fmt and clippy all PASS. The ordinary release harness test passed
its perturbation check and ignored the heavy sweep. README/figure freshness and
specification tests: 9 passed, exit 0. The generated validation block is unchanged;
this measurement-only study does not introduce an accepted implementation result.

The existing `bench_scf.py --require-identical` confirmation used five cases
(gamma, k222s, k222, k444s, si8), T=1/2/4/8/12 and one recorded/repeated pair
each: 25/25 identical energy, density SHA-256, SCF counts, Davidson steps,
H applications and certificate work, maximum energy difference 0 Ha. This is a
same-source repeat confirmation, not an independent pre/post implementation
comparison: no production source changed. Exact commands and return codes are
in `P40A-evidence/identity-commands.json`; the field-by-field audit is in
`identity-summary.json`.

All measurements use the pure-Rust Mac build and original fixture thread budgets.
Three independent collector processes share the machine; builds and confirmation
checks also overlapped early collection. No wall-time conclusion is drawn.
The harness binary hash, Rust version, environment, source HEAD, fixture and UPF
hashes are in `binary.json` and `provenance.json`. `preregistered.md` preserves
the protocol before the first measurement. OPTIONS records describe the mixer
and path before applying the documented per-run reference/profile controls;
run IDs and ROW tolerances specify the latter. Each log has a separate exit-code
and command/environment record. Output and analysis scripts must be archived
with this note; they are intentionally outside `target/` and contain no UPFs.

A passed strict certificate is not proof of agreement with the tighter reference
to 1e-12 Ha or of spectral completeness. A standard stopping estimate is not a
certified energy bound. Previously nonconvergent elongated/slab cases remain in
the survey, without a manufactured reference or a convergence-gain claim.
The ten tiny deterministic perturbations probe trajectory sensitivity, not a
population of materially different starts. No new QE run, MKL result, remote
speed acceptance, atomic starting density or tolerance adjustment is claimed.

For dev08 replication through the orchestrator: archive these source/evidence
files; build this ignored test with the approved pure native backend and then
MKL separately; record HEAD, binary and UPF hashes, rustc and command exit codes;
run the same A–F matrix and paired seeds, with Gamma explicitly on/off, keeping
strict certification and both tight reference levels. Report each backend's
hierarchical intervals and the frozen decision rule independently. Freeze F
from this survey for a genuinely independent confirmation, or explicitly label
a fresh B–D reselection exploratory. Wall-speed acceptance would additionally
require idle-node, equal-core, alternating comparisons against QE's best setup;
it cannot be inferred from these work counts.

## Measured results (Mac observations)

F selected C by the frozen rule. Full data: `P40A-evidence/analysis.json`.

| Arm | Case | Path | Profile | SCF min/median/max | H median | Fail / 10 | Max energy error (Ha) |
|---|---|---|---|---|---:|---:|---:|
| A | si-gamma | off | strict | 12 / 12 / 12 | 353 | 0 | 2.131628e-14 |
| A | si-gamma | off | standard | 9 / 9 / 9 | 180 | 0 | 1.045279e-10 |
| A | si-gamma | on | strict | 13 / 13 / 13 | 384 | 0 | 3.108624e-14 |
| A | si-gamma | on | standard | 9 / 9 / 9 | 180 | 0 | 1.047358e-10 |
| A | si-k222s | off | strict | 12 / 12 / 12 | 1581 | 0 | 2.842171e-14 |
| A | si-k222s | off | standard | 9 / 9 / 9 | 731 | 0 | 4.552607e-10 |
| A | si-skew | off | strict | 21 / 21 / 21 | 609 | 0 | 1.598721e-14 |
| A | si-skew | off | standard | 12 / 12 / 12 | 224 | 0 | 6.180780e-10 |
| A | si-skew | on | strict | 23 / 23 / 23 | 650 | 0 | 1.065814e-14 |
| A | si-skew | on | standard | 14 / 14 / 14 | 249 | 0 | 1.308475e-09 |
| A | si-elongated | off | strict | unavailable | 5499 | 10 | unavailable |
| A | si-elongated | off | standard | unavailable | 5165.5 | 10 | unavailable |
| A | si-elongated | on | strict | unavailable | 5254.5 | 10 | unavailable |
| A | si-elongated | on | standard | unavailable | 5198 | 10 | unavailable |
| A | si-slab | off | strict | unavailable | 6486 | 10 | unavailable |
| A | si-slab | off | standard | unavailable | 5420 | 10 | unavailable |
| A | si-slab | on | strict | unavailable | 6500 | 10 | unavailable |
| A | si-slab | on | standard | unavailable | 5281 | 10 | unavailable |
| A | si8-gamma | off | strict | 14 / 14 / 14 | 2059.5 | 0 | 2.451372e-13 |
| A | si8-gamma | off | standard | 20 / 25.5 / 32 | 1167.5 | 0 | 2.245319e-09 |
| A | si8-gamma | on | strict | 14 / 14 / 14 | 2058.5 | 0 | 1.207923e-13 |
| A | si8-gamma | on | standard | 21 / 27.5 / 30 | 1281.5 | 0 | 2.915563e-09 |
| A | si16-gamma | off | strict | 14 / 14 / 15 | 3223.5 | 0 | 3.694112e-11 |
| A | si16-gamma | off | standard | 18 / 21.5 / 27 | 1854.5 | 0 | 1.880501e-09 |
| A | si16-gamma | on | strict | 12 / 12 / 14 | 2765.5 | 0 | 1.855938e-11 |
| A | si16-gamma | on | standard | 15 / 19 / 25 | 1695.5 | 0 | 3.007159e-09 |
| A | si-low | off | strict | 12 / 12 / 12 | 306 | 0 | 2.930989e-14 |
| A | si-low | off | standard | 9 / 9 / 9 | 176 | 0 | 2.975309e-11 |
| A | si-low | on | strict | 11 / 11 / 11 | 294 | 0 | 2.842171e-14 |
| A | si-low | on | standard | 9 / 9 / 9 | 177 | 0 | 2.144240e-11 |
| A | material-si-20 | off | strict | 15 / 15 / 15 | 16526 | 0 | 3.730349e-14 |
| A | material-si-20 | off | standard | 11 / 11 / 11 | 6974 | 0 | 2.593357e-10 |
| A | material-mgo-40 | off | strict | 18 / 18 / 18 | 12202 | 0 | 3.836931e-13 |
| A | material-mgo-40 | off | standard | 13 / 13 / 13 | 5909 | 0 | 1.508766e-10 |
| B | si-gamma | off | strict | 10 / 10 / 10 | 319 | 0 | 1.669775e-13 |
| B | si-gamma | off | standard | 9 / 9 / 9 | 176 | 0 | 1.080673e-10 |
| B | si-gamma | on | strict | 12 / 12 / 12 | 377 | 0 | 5.151435e-14 |
| B | si-gamma | on | standard | 9 / 9 / 9 | 176 | 0 | 2.336069e-10 |
| B | si-k222s | off | strict | 10 / 10 / 10 | 1379 | 0 | 2.753353e-14 |
| B | si-k222s | off | standard | 9 / 9 / 9 | 731 | 0 | 4.688960e-11 |
| B | si-skew | off | strict | 13 / 13 / 13 | 410 | 0 | 2.664535e-14 |
| B | si-skew | off | standard | 17 / 17 / 17 | 304 | 0 | 2.806644e-11 |
| B | si-skew | on | strict | 17 / 17 / 17 | 498 | 0 | 1.598721e-14 |
| B | si-skew | on | standard | 16 / 16 / 16 | 286 | 0 | 4.329843e-10 |
| B | si-elongated | off | strict | unavailable | 5902 | 10 | unavailable |
| B | si-elongated | off | standard | unavailable | 5465 | 10 | unavailable |
| B | si-elongated | on | strict | unavailable | 5975 | 10 | unavailable |
| B | si-elongated | on | standard | unavailable | 5561 | 10 | unavailable |
| B | si-slab | off | strict | unavailable | 6749 | 10 | unavailable |
| B | si-slab | off | standard | unavailable | 5595.5 | 10 | unavailable |
| B | si-slab | on | strict | unavailable | 6824 | 10 | unavailable |
| B | si-slab | on | standard | unavailable | 5645.5 | 10 | unavailable |
| B | si8-gamma | off | strict | 15 / 16 / 17 | 2359.5 | 0 | 1.325162e-12 |
| B | si8-gamma | off | standard | 15 / 16 / 19 | 836.5 | 0 | 2.526292e-09 |
| B | si8-gamma | on | strict | 15 / 15 / 16 | 2252.5 | 0 | 1.257661e-12 |
| B | si8-gamma | on | standard | 15 / 16 / 18 | 851 | 0 | 1.951818e-09 |
| B | si16-gamma | off | strict | 12 / 15 / 15 | 3864.5 | 0 | 3.746692e-11 |
| B | si16-gamma | off | standard | 12 / 13 / 15 | 1350 | 0 | 2.245663e-09 |
| B | si16-gamma | on | strict | 14 / 15 / 16 | 3903.5 | 0 | 2.017941e-12 |
| B | si16-gamma | on | standard | 12 / 12.5 / 16 | 1305.5 | 0 | 2.300538e-09 |
| B | si-low | off | strict | 11 / 11 / 11 | 281 | 0 | 2.398082e-14 |
| B | si-low | off | standard | 9 / 9 / 9 | 175 | 0 | 2.692069e-11 |
| B | si-low | on | strict | 11 / 11 / 11 | 287 | 0 | 2.131628e-14 |
| B | si-low | on | standard | 9 / 9 / 9 | 179 | 0 | 2.972751e-10 |
| B | material-si-20 | off | strict | 14 / 14 / 14 | 15319 | 0 | 1.243450e-14 |
| B | material-si-20 | off | standard | 11 / 11 / 11 | 6923 | 0 | 3.702105e-11 |
| B | material-mgo-40 | off | strict | 15 / 15 / 15 | 10519 | 0 | 4.121148e-13 |
| B | material-mgo-40 | off | standard | 11 / 11 / 11 | 5197 | 0 | 1.755751e-10 |
| C | si-gamma | off | strict | 10 / 10 / 10 | 307 | 0 | 4.085621e-14 |
| C | si-gamma | off | standard | 7 / 7 / 7 | 147 | 0 | 1.724327e-10 |
| C | si-gamma | on | strict | 10 / 10 / 10 | 323.5 | 0 | 4.973799e-14 |
| C | si-gamma | on | standard | 7 / 7 / 7 | 140 | 0 | 6.109873e-09 |
| C | si-k222s | off | strict | 9 / 9 / 9 | 1264 | 0 | 6.039613e-14 |
| C | si-k222s | off | standard | 7 / 7 / 7 | 599 | 0 | 6.224177e-11 |
| C | si-skew | off | strict | 13 / 13 / 13 | 469 | 0 | 1.776357e-14 |
| C | si-skew | off | standard | 8 / 8 / 8 | 164 | 0 | 2.622800e-10 |
| C | si-skew | on | strict | 16 / 16 / 16 | 462 | 0 | 2.486900e-14 |
| C | si-skew | on | standard | 13 / 13 / 13 | 233 | 0 | 3.933236e-10 |
| C | si-elongated | off | strict | unavailable | 5232.5 | 10 | unavailable |
| C | si-elongated | off | standard | unavailable | 5082.5 | 10 | unavailable |
| C | si-elongated | on | strict | unavailable | 5061 | 10 | unavailable |
| C | si-elongated | on | standard | unavailable | 5291 | 10 | unavailable |
| C | si-slab | off | strict | unavailable | 6072.5 | 10 | unavailable |
| C | si-slab | off | standard | unavailable | 5813 | 10 | unavailable |
| C | si-slab | on | strict | unavailable | 6715 | 10 | unavailable |
| C | si-slab | on | standard | unavailable | 5832.5 | 10 | unavailable |
| C | si8-gamma | off | strict | 16 / 17 / 17 | 2688 | 0 | 1.644906e-12 |
| C | si8-gamma | off | standard | 16 / 17 / 19 | 936.5 | 0 | 3.752323e-09 |
| C | si8-gamma | on | strict | 16 / 17 / 18 | 2724.5 | 0 | 1.350031e-12 |
| C | si8-gamma | on | standard | 16 / 17 / 21 | 934 | 0 | 3.105914e-09 |
| C | si16-gamma | off | strict | 14 / 15 / 16 | 4080.5 | 0 | 2.941647e-12 |
| C | si16-gamma | off | standard | 13 / 13 / 14 | 1432 | 0 | 1.982613e-09 |
| C | si16-gamma | on | strict | 15 / 16 / 17 | 4423.5 | 0 | 5.329071e-12 |
| C | si16-gamma | on | standard | 13 / 13 / 16 | 1447.5 | 0 | 3.302198e-09 |
| C | si-low | off | strict | 10 / 10 / 10 | 261 | 0 | 2.486900e-14 |
| C | si-low | off | standard | 7 / 7 / 7 | 137 | 0 | 2.173168e-10 |
| C | si-low | on | strict | 11 / 11 / 11 | 291 | 0 | 1.598721e-14 |
| C | si-low | on | standard | 8 / 8 / 8 | 152 | 0 | 9.830980e-11 |
| C | material-si-20 | off | strict | 13 / 13 / 13 | 15784 | 0 | 2.309264e-14 |
| C | material-si-20 | off | standard | 10 / 10 / 10 | 6474 | 0 | 1.075637e-10 |
| C | material-mgo-40 | off | strict | 13 / 13 / 13 | 8733 | 0 | 7.105427e-14 |
| C | material-mgo-40 | off | standard | 11 / 11 / 11 | 5225 | 0 | 2.734168e-11 |
| D | si-gamma | off | strict | 10 / 10 / 10 | 319 | 0 | 1.669775e-13 |
| D | si-gamma | off | standard | 9 / 9 / 9 | 176 | 0 | 1.080673e-10 |
| D | si-gamma | on | strict | 12 / 12 / 12 | 377 | 0 | 5.151435e-14 |
| D | si-gamma | on | standard | 9 / 9 / 9 | 176 | 0 | 2.336069e-10 |
| D | si-k222s | off | strict | 10 / 10 / 10 | 1410 | 0 | 4.796163e-14 |
| D | si-k222s | off | standard | 9 / 9 / 9 | 726 | 0 | 2.389742e-10 |
| D | si-skew | off | strict | 13 / 13 / 13 | 408 | 0 | 2.131628e-14 |
| D | si-skew | off | standard | 17 / 17 / 17 | 305 | 0 | 1.312752e-09 |
| D | si-skew | on | strict | 16 / 16 / 16 | 451 | 0 | 2.753353e-14 |
| D | si-skew | on | standard | 18 / 18 / 18 | 320 | 0 | 1.187157e-10 |
| D | si-elongated | off | strict | unavailable | 4698 | 10 | unavailable |
| D | si-elongated | off | standard | unavailable | 5092.5 | 10 | unavailable |
| D | si-elongated | on | strict | unavailable | 4646.5 | 10 | unavailable |
| D | si-elongated | on | standard | unavailable | 5069.5 | 10 | unavailable |
| D | si-slab | off | strict | unavailable | 6001.5 | 10 | unavailable |
| D | si-slab | off | standard | unavailable | 5438.5 | 10 | unavailable |
| D | si-slab | on | strict | unavailable | 6197 | 10 | unavailable |
| D | si-slab | on | standard | unavailable | 5304.5 | 10 | unavailable |
| D | si8-gamma | off | strict | 15 / 16 / 16 | 2355 | 0 | 9.663381e-13 |
| D | si8-gamma | off | standard | 14 / 16 / 20 | 827.5 | 0 | 2.271356e-09 |
| D | si8-gamma | on | strict | 15 / 16 / 16 | 2401.5 | 0 | 1.300293e-12 |
| D | si8-gamma | on | standard | 15 / 16 / 18 | 831 | 0 | 6.475346e-09 |
| D | si16-gamma | off | strict | 14 / 15 / 15 | 4244 | 0 | 1.129763e-12 |
| D | si16-gamma | off | standard | 12 / 13 / 16 | 1354 | 0 | 2.834845e-09 |
| D | si16-gamma | on | strict | 14 / 14.5 / 15 | 3694 | 0 | 2.763301e-11 |
| D | si16-gamma | on | standard | 12 / 13 / 17 | 1351 | 0 | 2.723070e-09 |
| D | si-low | off | strict | 12 / 12 / 12 | 295 | 0 | 1.598721e-14 |
| D | si-low | off | standard | 9 / 9 / 9 | 175 | 0 | 3.503065e-11 |
| D | si-low | on | strict | 11 / 11 / 11 | 291 | 0 | 1.776357e-14 |
| D | si-low | on | standard | 9 / 9 / 9 | 179 | 0 | 4.591794e-10 |
| D | material-si-20 | off | strict | 14 / 14 / 14 | 15237 | 0 | 3.552714e-15 |
| D | material-si-20 | off | standard | 11 / 11 / 11 | 6941 | 0 | 2.185097e-11 |
| D | material-mgo-40 | off | strict | 14 / 14 / 14 | 9268 | 0 | 6.110668e-13 |
| D | material-mgo-40 | off | standard | 11 / 11 / 11 | 5154 | 0 | 3.827125e-10 |
| E | si-gamma | off | strict | 12 / 12 / 12 | 358 | 0 | 3.907985e-14 |
| E | si-gamma | off | standard | 9 / 9 / 9 | 176 | 0 | 9.746977e-10 |
| E | si-gamma | on | strict | 13 / 13 / 13 | 376 | 0 | 3.375078e-14 |
| E | si-gamma | on | standard | 9 / 9 / 9 | 176 | 0 | 9.431211e-10 |
| E | si-k222s | off | strict | 12 / 12 / 12 | 1549 | 0 | 4.085621e-14 |
| E | si-k222s | off | standard | 9 / 9 / 9 | 745 | 0 | 3.738316e-10 |
| E | si-skew | off | strict | 24 / 24 / 24 | 600 | 0 | 2.575717e-14 |
| E | si-skew | off | standard | 18 / 18 / 18 | 330 | 0 | 1.603411e-10 |
| E | si-skew | on | strict | 27 / 27 / 27 | 687 | 0 | 1.598721e-14 |
| E | si-skew | on | standard | 20 / 20 / 20 | 349 | 0 | 5.165468e-11 |
| E | si-elongated | off | strict | unavailable | 4867 | 10 | unavailable |
| E | si-elongated | off | standard | unavailable | 4964 | 10 | unavailable |
| E | si-elongated | on | strict | unavailable | 4734 | 10 | unavailable |
| E | si-elongated | on | standard | unavailable | 5134.5 | 10 | unavailable |
| E | si-slab | off | strict | unavailable | 6564.5 | 10 | unavailable |
| E | si-slab | off | standard | unavailable | 5239 | 10 | unavailable |
| E | si-slab | on | strict | unavailable | 6298.5 | 10 | unavailable |
| E | si-slab | on | standard | unavailable | 5205 | 10 | unavailable |
| E | si8-gamma | off | strict | 14 / 14 / 15 | 1950.5 | 0 | 1.172396e-13 |
| E | si8-gamma | off | standard | 22 / 27.5 / 32 | 1264.5 | 0 | 1.666601e-08 |
| E | si8-gamma | on | strict | 14 / 14 / 15 | 1986.5 | 0 | 1.989520e-13 |
| E | si8-gamma | on | standard | 19 / 26 / 34 | 1181.5 | 0 | 3.031360e-09 |
| E | si16-gamma | off | strict | 13 / 13 / 14 | 3011.5 | 0 | 8.220979e-12 |
| E | si16-gamma | off | standard | 15 / 16.5 / 22 | 1480 | 0 | 2.561990e-09 |
| E | si16-gamma | on | strict | 13 / 13 / 13 | 2831.5 | 0 | 3.851142e-12 |
| E | si16-gamma | on | standard | 14 / 19 / 22 | 1639 | 0 | 2.653749e-09 |
| E | si-low | off | strict | 12 / 12 / 12 | 305 | 0 | 2.398082e-14 |
| E | si-low | off | standard | 9 / 9 / 9 | 176 | 0 | 5.407497e-11 |
| E | si-low | on | strict | 13 / 13 / 13 | 312 | 0 | 1.421085e-14 |
| E | si-low | on | standard | 9 / 9 / 9 | 179 | 0 | 2.807354e-11 |
| E | material-si-20 | off | strict | 15 / 15 / 15 | 16561 | 0 | 3.375078e-14 |
| E | material-si-20 | off | standard | 11 / 11 / 11 | 7005 | 0 | 1.227018e-10 |
| E | material-mgo-40 | off | strict | 18 / 18 / 18 | 12382 | 0 | 1.122658e-12 |
| E | material-mgo-40 | off | standard | 14 / 14 / 14 | 6613 | 0 | 1.301714e-11 |
| F | si-gamma | off | strict | 10 / 10 / 10 | 307 | 0 | 3.552714e-14 |
| F | si-gamma | off | standard | 6 / 6 / 6 | 130 | 0 | 5.291980e-10 |
| F | si-gamma | on | strict | 10 / 10 / 10 | 328 | 0 | 3.907985e-14 |
| F | si-gamma | on | standard | 7 / 7 / 7 | 143 | 0 | 6.776766e-10 |
| F | si-k222s | off | strict | 9 / 9 / 9 | 1256 | 0 | 3.907985e-14 |
| F | si-k222s | off | standard | 7 / 7 / 7 | 616 | 0 | 1.291482e-10 |
| F | si-skew | off | strict | 13 / 13 / 13 | 369 | 0 | 1.598721e-14 |
| F | si-skew | off | standard | 8 / 8 / 8 | 164 | 0 | 4.842722e-10 |
| F | si-skew | on | strict | 15 / 15 / 15 | 440 | 0 | 5.329071e-15 |
| F | si-skew | on | standard | 14 / 14 / 14 | 247 | 0 | 2.790559e-10 |
| F | si-elongated | off | strict | unavailable | 4464.5 | 10 | unavailable |
| F | si-elongated | off | standard | unavailable | 4963 | 10 | unavailable |
| F | si-elongated | on | strict | unavailable | 4835.5 | 10 | unavailable |
| F | si-elongated | on | standard | unavailable | 5015 | 10 | unavailable |
| F | si-slab | off | strict | unavailable | 6107.5 | 10 | unavailable |
| F | si-slab | off | standard | unavailable | 5452.5 | 10 | unavailable |
| F | si-slab | on | strict | unavailable | 5360.5 | 10 | unavailable |
| F | si-slab | on | standard | unavailable | 5439.5 | 10 | unavailable |
| F | si8-gamma | off | strict | 15 / 16 / 17 | 2516.5 | 0 | 1.712408e-12 |
| F | si8-gamma | off | standard | 14 / 16 / 16 | 868.5 | 0 | 2.746447e-09 |
| F | si8-gamma | on | strict | 15 / 16 / 16 | 2474.5 | 0 | 1.996625e-12 |
| F | si8-gamma | on | standard | 15 / 16 / 17 | 864.5 | 0 | 2.683603e-09 |
| F | si16-gamma | off | strict | 14 / 14 / 15 | 4309 | 0 | 7.467804e-12 |
| F | si16-gamma | off | standard | 12 / 13 / 14 | 1392.5 | 0 | 2.365418e-09 |
| F | si16-gamma | on | strict | 15 / 15 / 15 | 4206 | 0 | 2.756906e-12 |
| F | si16-gamma | on | standard | 12 / 13 / 14 | 1410 | 0 | 3.304883e-09 |
| F | si-low | off | strict | 10 / 10 / 10 | 262 | 0 | 2.398082e-14 |
| F | si-low | off | standard | 7 / 7 / 7 | 139 | 0 | 1.229958e-10 |
| F | si-low | on | strict | 11 / 11 / 11 | 289 | 0 | 2.131628e-14 |
| F | si-low | on | standard | 8 / 8 / 8 | 155 | 0 | 3.454304e-11 |
| F | material-si-20 | off | strict | 12 / 12 / 12 | 13676 | 0 | 3.730349e-14 |
| F | material-si-20 | off | standard | 9 / 9 / 9 | 5888 | 0 | 1.130720e-09 |
| F | material-mgo-40 | off | strict | 13 / 13 / 13 | 9533 | 0 | 5.258016e-13 |
| F | material-mgo-40 | off | standard | 10 / 10 / 10 | 4768 | 0 | 1.114330e-09 |

### Strict reference qualification

Both levels must certify and plateau within 1e-12 Ha. Cross-arm delta uses the A tight result on the same path; unavailable means no qualified pair.

| Arm | Case | Path | Reference failures / 2 | Tight energy (Ha) | Plateau difference (Ha) | Tight minus A (Ha) |
|---|---|---|---:|---:|---:|---:|
| A | si-gamma | off | 0 | -7.29553046821091300e+00 | 8.881784e-15 | +0.000000e+00 |
| A | si-gamma | on | 0 | -7.29553046821090412e+00 | 1.332268e-14 | +0.000000e+00 |
| A | si-k222s | off | 0 | -7.91942341246910964e+00 | 3.641532e-14 | +0.000000e+00 |
| A | si-skew | off | 0 | -7.25055567296080383e+00 | 8.881784e-16 | +0.000000e+00 |
| A | si-skew | on | 0 | -7.25055567296081271e+00 | 1.154632e-14 | +0.000000e+00 |
| A | si-elongated | off | 2 | unavailable | unavailable | unavailable |
| A | si-elongated | on | 2 | unavailable | unavailable | unavailable |
| A | si-slab | off | 2 | unavailable | unavailable | unavailable |
| A | si-slab | on | 2 | unavailable | unavailable | unavailable |
| A | si8-gamma | off | 0 | -3.13218036010317320e+01 | 1.776357e-14 | +0.000000e+00 |
| A | si8-gamma | on | 0 | -3.13218036010317391e+01 | 5.684342e-14 | +0.000000e+00 |
| A | si16-gamma | off | 0 | -6.31292444216179334e+01 | 1.847411e-13 | +0.000000e+00 |
| A | si16-gamma | on | 0 | -6.31292444216179263e+01 | 8.526513e-14 | +0.000000e+00 |
| A | si-low | off | 0 | -7.23375645282481461e+00 | 7.105427e-15 | +0.000000e+00 |
| A | si-low | on | 0 | -7.23375645282481194e+00 | 7.105427e-15 | +0.000000e+00 |
| A | material-si-20 | off | 0 | -8.46231783690607564e+00 | 5.329071e-15 | +0.000000e+00 |
| A | material-mgo-40 | off | 0 | -7.60508808517109998e+01 | 7.105427e-14 | +0.000000e+00 |
| B | si-gamma | off | 0 | -7.29553046821091300e+00 | 1.065814e-14 | +0.000000e+00 |
| B | si-gamma | on | 0 | -7.29553046821093076e+00 | 1.332268e-14 | -2.664535e-14 |
| B | si-k222s | off | 0 | -7.91942341246906878e+00 | 4.085621e-14 | +4.085621e-14 |
| B | si-skew | off | 0 | -7.25055567296081360e+00 | 5.329071e-15 | -9.769963e-15 |
| B | si-skew | on | 0 | -7.25055567296078074e+00 | 1.598721e-14 | +3.197442e-14 |
| B | si-elongated | off | 2 | unavailable | unavailable | unavailable |
| B | si-elongated | on | 2 | unavailable | unavailable | unavailable |
| B | si-slab | off | 2 | unavailable | unavailable | unavailable |
| B | si-slab | on | 2 | unavailable | unavailable | unavailable |
| B | si8-gamma | off | 0 | -3.13218036010318670e+01 | 7.815970e-14 | -1.350031e-13 |
| B | si8-gamma | on | 0 | -3.13218036010317959e+01 | 1.421085e-13 | -5.684342e-14 |
| B | si16-gamma | off | 0 | -6.31292444216178765e+01 | 1.207923e-13 | +5.684342e-14 |
| B | si16-gamma | on | 0 | -6.31292444216178481e+01 | 2.131628e-14 | +7.815970e-14 |
| B | si-low | off | 0 | -7.23375645282481017e+00 | 7.993606e-15 | +4.440892e-15 |
| B | si-low | on | 0 | -7.23375645282481194e+00 | 4.440892e-15 | +0.000000e+00 |
| B | material-si-20 | off | 0 | -8.46231783690607742e+00 | 0.000000e+00 | -1.776357e-15 |
| B | material-mgo-40 | off | 0 | -7.60508808517108150e+01 | 4.405365e-13 | +1.847411e-13 |
| C | si-gamma | off | 0 | -7.29553046821093432e+00 | 3.019807e-14 | -2.131628e-14 |
| C | si-gamma | on | 0 | -7.29553046821089257e+00 | 2.664535e-15 | +1.154632e-14 |
| C | si-k222s | off | 0 | -7.91942341246912296e+00 | 3.286260e-14 | -1.332268e-14 |
| C | si-skew | off | 0 | -7.25055567296078074e+00 | 1.421085e-14 | +2.309264e-14 |
| C | si-skew | on | 0 | -7.25055567296079850e+00 | 1.776357e-15 | +1.421085e-14 |
| C | si-elongated | off | 2 | unavailable | unavailable | unavailable |
| C | si-elongated | on | 2 | unavailable | unavailable | unavailable |
| C | si-slab | off | 2 | unavailable | unavailable | unavailable |
| C | si-slab | on | 2 | unavailable | unavailable | unavailable |
| C | si8-gamma | off | 0 | -3.13218036010317604e+01 | 4.618528e-14 | -2.842171e-14 |
| C | si8-gamma | on | 0 | -3.13218036010317178e+01 | 4.618528e-14 | +2.131628e-14 |
| C | si16-gamma | off | 0 | -6.31292444216177770e+01 | 1.207923e-13 | +1.563194e-13 |
| C | si16-gamma | on | 0 | -6.31292444216174715e+01 | 2.842171e-13 | +4.547474e-13 |
| C | si-low | off | 0 | -7.23375645282480484e+00 | 7.105427e-15 | +9.769963e-15 |
| C | si-low | on | 0 | -7.23375645282479418e+00 | 1.776357e-15 | +1.776357e-14 |
| C | material-si-20 | off | 0 | -8.46231783690606143e+00 | 1.243450e-14 | +1.421085e-14 |
| C | material-mgo-40 | off | 0 | -7.60508808517109571e+01 | 5.826450e-13 | +4.263256e-14 |
| D | si-gamma | off | 0 | -7.29553046821089968e+00 | 9.769963e-15 | +1.332268e-14 |
| D | si-gamma | on | 0 | -7.29553046821092899e+00 | 3.552714e-14 | -2.486900e-14 |
| D | si-k222s | off | 0 | -7.91942341246907144e+00 | 4.174439e-14 | +3.819167e-14 |
| D | si-skew | off | 0 | -7.25055567296080561e+00 | 3.552714e-15 | -1.776357e-15 |
| D | si-skew | on | 0 | -7.25055567296082426e+00 | 3.730349e-14 | -1.154632e-14 |
| D | si-elongated | off | 2 | unavailable | unavailable | unavailable |
| D | si-elongated | on | 2 | unavailable | unavailable | unavailable |
| D | si-slab | off | 2 | unavailable | unavailable | unavailable |
| D | si-slab | on | 2 | unavailable | unavailable | unavailable |
| D | si8-gamma | off | 0 | -3.13218036010318741e+01 | 4.973799e-14 | -1.421085e-13 |
| D | si8-gamma | on | 0 | -3.13218036010317640e+01 | 3.552714e-15 | -2.486900e-14 |
| D | si16-gamma | off | 0 | -6.31292444216178126e+01 | 2.842171e-14 | +1.207923e-13 |
| D | si16-gamma | on | 0 | -6.31292444216176989e+01 | 5.684342e-14 | +2.273737e-13 |
| D | si-low | off | 0 | -7.23375645282480484e+00 | 1.776357e-15 | +9.769963e-15 |
| D | si-low | on | 0 | -7.23375645282479951e+00 | 1.776357e-14 | +1.243450e-14 |
| D | material-si-20 | off | 0 | -8.46231783690607386e+00 | 3.552714e-15 | +1.776357e-15 |
| D | material-mgo-40 | off | 0 | -7.60508808517112982e+01 | 1.563194e-13 | -2.984279e-13 |
| E | si-gamma | off | 0 | -7.29553046821091300e+00 | 2.309264e-14 | +0.000000e+00 |
| E | si-gamma | on | 0 | -7.29553046821090856e+00 | 6.217249e-15 | -4.440892e-15 |
| E | si-k222s | off | 0 | -7.91942341246908477e+00 | 3.197442e-14 | +2.486900e-14 |
| E | si-skew | off | 0 | -7.25055567296079317e+00 | 2.309264e-14 | +1.065814e-14 |
| E | si-skew | on | 0 | -7.25055567296078607e+00 | 9.769963e-15 | +2.664535e-14 |
| E | si-elongated | off | 2 | unavailable | unavailable | unavailable |
| E | si-elongated | on | 2 | unavailable | unavailable | unavailable |
| E | si-slab | off | 2 | unavailable | unavailable | unavailable |
| E | si-slab | on | 2 | unavailable | unavailable | unavailable |
| E | si8-gamma | off | 0 | -3.13218036010318244e+01 | 6.394885e-14 | -9.237056e-14 |
| E | si8-gamma | on | 0 | -3.13218036010319416e+01 | 1.740830e-13 | -2.025047e-13 |
| E | si16-gamma | off | 0 | -6.31292444216177984e+01 | 1.278977e-13 | +1.350031e-13 |
| E | si16-gamma | on | 0 | -6.31292444216177699e+01 | 5.684342e-14 | +1.563194e-13 |
| E | si-low | off | 0 | -7.23375645282480306e+00 | 1.776357e-14 | +1.154632e-14 |
| E | si-low | on | 0 | -7.23375645282479596e+00 | 6.217249e-15 | +1.598721e-14 |
| E | material-si-20 | off | 0 | -8.46231783690607386e+00 | 1.421085e-14 | +1.776357e-15 |
| E | material-mgo-40 | off | 0 | -7.60508808517117245e+01 | 6.110668e-13 | -7.247536e-13 |
| F | si-gamma | off | 0 | -7.29553046821091122e+00 | 0.000000e+00 | +1.776357e-15 |
| F | si-gamma | on | 0 | -7.29553046821091655e+00 | 3.552714e-15 | -1.243450e-14 |
| F | si-k222s | off | 0 | -7.91942341246908832e+00 | 3.907985e-14 | +2.131628e-14 |
| F | si-skew | off | 0 | -7.25055567296080383e+00 | 4.440892e-15 | +0.000000e+00 |
| F | si-skew | on | 0 | -7.25055567296080028e+00 | 5.329071e-15 | +1.243450e-14 |
| F | si-elongated | off | 2 | unavailable | unavailable | unavailable |
| F | si-elongated | on | 2 | unavailable | unavailable | unavailable |
| F | si-slab | off | 2 | unavailable | unavailable | unavailable |
| F | si-slab | on | 2 | unavailable | unavailable | unavailable |
| F | si8-gamma | off | 0 | -3.13218036010317427e+01 | 6.750156e-14 | -1.065814e-14 |
| F | si8-gamma | on | 0 | -3.13218036010316823e+01 | 1.136868e-13 | +5.684342e-14 |
| F | si16-gamma | off | 0 | -6.31292444216178907e+01 | 4.973799e-14 | +4.263256e-14 |
| F | si16-gamma | on | 0 | -6.31292444216175213e+01 | 2.771117e-13 | +4.050094e-13 |
| F | si-low | off | 0 | -7.23375645282481372e+00 | 1.776357e-14 | +8.881784e-16 |
| F | si-low | on | 0 | -7.23375645282481106e+00 | 1.332268e-14 | +8.881784e-16 |
| F | material-si-20 | off | 0 | -8.46231783690607564e+00 | 1.776357e-15 | +0.000000e+00 |
| F | material-mgo-40 | off | 0 | -7.60508808517108292e+01 | 6.394885e-13 | +1.705303e-13 |

### Paired changes and hierarchical bootstrap

Negative percentages mean fewer iterations or H vectors. Successful paired SCF only; all attempts for H. Cases and pairs show coverage.

| Arm | Path | Profile | Metric | Paired median % [95% interval] | Cases / pairs |
|---|---|---|---|---|---|
| B | off | strict | iterations | -12.500 [-16.667, +3.571] | 8 / 80 |
| B | off | strict | h | -7.737 [-12.777, +9.992] | 10 / 100 |
| B | off | standard | iterations | +0.000 [-26.786, +0.000] | 8 / 80 |
| B | off | standard | h | -0.650 [-12.049, +4.128] | 10 / 100 |
| B | on | strict | iterations | +0.000 [-26.087, +16.026] | 5 / 50 |
| B | on | strict | h | +5.960 [-2.381, +16.295] | 7 / 70 |
| B | on | standard | iterations | +0.000 [-36.250, +14.286] | 5 / 50 |
| B | on | standard | h | +1.130 [-18.911, +7.721] | 7 / 70 |
| C | off | strict | iterations | -16.667 [-27.778, +4.048] | 8 / 80 |
| C | off | strict | h | -8.831 [-21.469, +7.549] | 10 / 100 |
| C | off | standard | iterations | -22.222 [-33.333, -15.385] | 8 / 80 |
| C | off | standard | h | -16.634 [-22.159, -7.154] | 10 / 100 |
| C | on | strict | iterations | +0.000 [-30.435, +27.885] | 5 / 50 |
| C | on | strict | h | +3.308 [-15.885, +37.913] | 7 / 70 |
| C | on | standard | iterations | -22.222 [-31.136, -7.143] | 5 / 50 |
| C | on | standard | h | -12.463 [-22.222, +3.382] | 7 / 70 |
| D | off | strict | iterations | -11.667 [-22.222, +7.143] | 8 / 80 |
| D | off | strict | h | -8.749 [-17.431, +2.698] | 10 / 100 |
| D | off | standard | iterations | +0.000 [-25.481, +0.000] | 8 / 80 |
| D | off | standard | h | -1.203 [-12.777, -0.473] | 10 / 100 |
| D | on | strict | iterations | +0.000 [-30.435, +16.026] | 5 / 50 |
| D | on | strict | h | -1.823 [-13.437, +17.967] | 7 / 70 |
| D | on | standard | iterations | +0.000 [-38.000, +28.571] | 5 / 50 |
| D | on | standard | h | -1.421 [-17.274, +2.751] | 7 / 70 |
| E | off | strict | iterations | +0.000 [+0.000, +0.000] | 8 / 80 |
| E | off | strict | h | -0.902 [-4.905, +0.460] | 10 / 100 |
| E | off | standard | iterations | +0.000 [+0.000, +13.462] | 8 / 80 |
| E | off | standard | h | +0.222 [-2.440, +11.914] | 10 / 100 |
| E | on | strict | iterations | +8.333 [+0.000, +18.182] | 5 / 50 |
| E | on | strict | h | -2.083 [-4.387, +5.692] | 7 / 70 |
| E | on | standard | iterations | +0.000 [-5.000, +42.857] | 5 / 50 |
| E | on | standard | h | -1.271 [-3.219, +9.046] | 7 / 70 |
| F | off | strict | iterations | -18.333 [-27.778, +0.000] | 8 / 80 |
| F | off | strict | h | -14.550 [-21.215, +3.731] | 10 / 100 |
| F | off | standard | iterations | -28.205 [-35.282, -22.222] | 8 / 80 |
| F | off | standard | h | -20.166 [-26.786, -10.673] | 10 / 100 |
| F | on | strict | iterations | +0.000 [-34.783, +20.192] | 5 / 50 |
| F | on | strict | h | -3.898 [-16.444, +23.337] | 7 / 70 |
| F | on | standard | iterations | -22.222 [-40.000, +0.000] | 5 / 50 |
| F | on | standard | h | -12.429 [-21.455, -0.803] | 7 / 70 |

### Frozen decision rule evaluation

| Arm | Path | Profile | Cases improved | >10% regression | Unavailable medians | Accuracy failed/unavailable | New failed starts | Pass |
|---|---|---|---|---|---|---|---:|---|
| B | off | strict | 6/10 | si8-gamma | si-elongated, si-slab | si-gamma, si-skew, si-elongated:unavailable, si-slab:unavailable, si8-gamma | 0 | False |
| B | off | standard | 3/10 | si-skew | si-elongated, si-slab | si-elongated:unavailable, si-slab:unavailable | 0 | False |
| B | on | strict | 2/7 | si16-gamma | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable, si8-gamma | 0 | False |
| B | on | standard | 2/7 | si-skew | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable, si-low | 0 | False |
| C | off | strict | 6/10 | si8-gamma | si-elongated, si-slab | si-gamma, si-k222s, si-elongated:unavailable, si-slab:unavailable, si8-gamma | 0 | False |
| C | off | standard | 8/10 | none | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable, si8-gamma, si-low | 0 | False |
| C | on | strict | 2/7 | si8-gamma, si16-gamma | si-elongated, si-slab | si-gamma, si-skew, si-elongated:unavailable, si-slab:unavailable, si8-gamma | 0 | False |
| C | on | standard | 5/7 | none | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable, si-low | 0 | False |
| D | off | strict | 5/10 | si8-gamma | si-elongated, si-slab | si-gamma, si-k222s, si-elongated:unavailable, si-slab:unavailable, si8-gamma, material-mgo-40 | 0 | False |
| D | off | standard | 3/10 | si-skew | si-elongated, si-slab | si-skew, si-elongated:unavailable, si-slab:unavailable, si16-gamma, material-mgo-40 | 0 | False |
| D | on | strict | 2/7 | si8-gamma, si16-gamma | si-elongated, si-slab | si-gamma, si-skew, si-elongated:unavailable, si-slab:unavailable, si8-gamma | 0 | False |
| D | on | standard | 2/7 | si-skew | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable, si8-gamma, si-low | 0 | False |
| E | off | strict | 1/10 | si-skew | si-elongated, si-slab | si-gamma, si-skew, si-elongated:unavailable, si-slab:unavailable, material-mgo-40 | 0 | False |
| E | off | standard | 1/10 | si-skew | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable, si8-gamma, si-low | 0 | False |
| E | on | strict | 0/7 | si-skew, si-low | si-elongated, si-slab | si-elongated:unavailable, si-slab:unavailable, si8-gamma | 0 | False |
| E | on | standard | 1/7 | si-skew | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable | 0 | False |
| F | off | strict | 6/10 | si8-gamma | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable, si8-gamma | 0 | False |
| F | off | standard | 8/10 | none | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable, si-low, material-si-20, material-mgo-40 | 0 | False |
| F | on | strict | 2/7 | si8-gamma, si16-gamma | si-elongated, si-slab | si-elongated:unavailable, si-slab:unavailable, si8-gamma | 0 | False |
| F | on | standard | 4/7 | none | si-elongated, si-slab | si-gamma, si-elongated:unavailable, si-slab:unavailable, si-low | 0 | False |
