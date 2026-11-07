# P29A standard-profile stagnation study

2026-11-07. **Stagnation reproduced; guard is an experimental, default-off
robustness tradeoff, not an accepted speed optimization.** On Si8, the guard
reduces the worst of ten trajectories from 32 to 27 SCF iterations, but median
H work increases 0.99%. It does not cure the two already nonconvergent stress
fixtures and makes them substantially more expensive. Do not enable it generally.

All new numbers below are local Mac/pure-Rust observations at source baseline
`d0756cd30796c367415fdb5bc602dfb4860c9430`, not MKL or QE measurements.
Raw logs, source snapshots, analysis programs, commands and hashes are in
[`P29A-evidence/`](../../P29A-evidence/) outside `target/`, per the architect's
explicit storage decision. The orchestrator archives that directory on merge.
No UPF is copied into the worktree or committed.

## Protocol and scope

The ignored release test `standard_stall_trajectory` uses the existing public
`initialize_density_transfer` API. Ten deterministic seeds multiply the uniform
starting density by `1 + delta`, with mean-centered LCG noise scaled by 1e-10.
Each relative perturbation is <=1e-10; the existing API validates charge and
feasibility. No random seed, Davidson algorithm, mixer or production initialization
hook is changed. No thread perturbations were needed to reproduce the pattern.

Eight exact benchmark geometries/cutoffs/meshes come from `benchmarks/cases`:
si-gamma, si-k222s, si-skew, si-elongated, si-slab, si8-gamma, si16-gamma, si-low.
These use the external Si.pz-vbc.UPF, engine defaults, T=1, standard threshold
5e-9 Ha, complex Gamma, and existing time-reversal mesh reduction. The two
additional material fixtures are `tests/reference/materials/si-20.toml` and
`mgo-40.toml`, using their actual PseudoDojo NC UPFs, PBE, cutoffs 40/80 Ry,
shifted 4x4x4/3x3x3 meshes, density floor 1e-10, Kerker-Pulay beta 0.3,
history 6, T=4 and 120-iteration budget as specified by those fixtures.
Other solver/profile settings remain at defaults. Thus the study is not restricted
to silicon, cubic cells, Gamma or one gap size. No synthetic replacement was
needed because all real files were present. Shapes include cubic, sheared,
elongated and slab-like cells, with 2, 4, 8 and 16 atoms.

Runs 0/1 are unperturbed strict reference candidates at density tolerances
1e-9/1e-10, energy tolerances 1e-11/1e-12 Ha and eigen tolerances 1e-10/1e-11 Ha.
Runs 2–11 are standard seeds 1–10. Both strict candidates must converge and
agree within 1e-12 Ha before an energy difference is reported. Eight cases pass;
elongated and slab exhaust even with strict settings, so their delivered
accuracy is **unavailable**, not zero. Certificate code and default tolerances
were never modified. Every standard outcome remains explicitly uncertified.
A final Hartree error <5e-9 Ha is a stopping condition, not a mathematical
bound on total-energy error; measured energy differences are reported separately.

The prespecified density stagnation event is a six-record window (five SCF
intervals) whose ending residual exceeds half its starting residual. A separate
MKL-like event additionally requires three consecutive bit-identical requested
tolerances inside that window. Counts below are trajectories with at least one
window, not overlapping-window totals. We also report successful trajectories
above 1.5 times their case median in the raw JSON; none occur. Failed trajectories
are excluded from successful SCF distributions but included in work distributions.
The failed elongated trajectories stop at 108/200/200 iterations (off) and
88/200/200 (on): each has three divergence errors and seven budget exhaustions.
Slab exhausts at 200/200/200 in both arms. These are censored attempt lengths,
not successful convergence counts.

Work means Davidson steps and Hamiltonian vector applications in recorded SCF
solves, excluding setup and reference certificate work. The final failed record
is included where supplied by the engine. No wall-time speed claim is made.

## Phase 1: guard off

Ten perturbations per case, 100 trajectories total. Eighty converge; the remaining
20 belong to elongated/slab. 32/100 have a density-stagnation event and 30/100
have an MKL-like event. In reference-backed successful runs those counts are
12/80 and 10/80. The Si8 count distribution and plateau pattern reproduce the
reported weakness without MKL; they do not establish that MKL itself is faulty.

All triples are min / median / max. Work columns include failures.

| Case | Converged / 10 | SCF (successful) | Davidson steps | H vectors | Density / MKL-like events |
|---|---:|---|---|---|---|
| si-gamma | 10 | 9 / 9 / 9 | 27 / 27 / 27 | 180 / 180 / 180 | 0 / 0 |
| si-k222s | 10 | 9 / 9 / 9 | 106 / 106 / 106 | 731 / 731 / 731 | 0 / 0 |
| si-skew | 10 | 12 / 12 / 12 | 32 / 32 / 32 | 224 / 224 / 224 | 0 / 0 |
| si-elongated | 0 | unavailable | 348 / 531.5 / 666 | 3216 / 5165.5 / 5878 | 10 / 10 |
| si-slab | 0 | unavailable | 471 / 506 / 550 | 5170 / 5420 / 5607 | 10 / 10 |
| si8-gamma | 10 | 20 / 25.5 / 32 | 86 / 104.5 / 124 | 958 / 1167.5 / 1442 | 8 / 8 |
| si16-gamma | 10 | 18 / 21.5 / 27 | 75 / 101.5 / 119 | 1580 / 1854.5 / 2318 | 4 / 2 |
| si-low | 10 | 9 / 9 / 9 | 24 / 24 / 24 | 176 / 176 / 176 | 0 / 0 |
| material-si-20 | 10 | 11 / 11 / 11 | 1042 / 1042 / 1042 | 6974 / 6974 / 6974 | 0 / 0 |
| material-mgo-40 | 10 | 13 / 13 / 13 | 609 / 609 / 609 | 5909 / 5909 / 5909 | 0 / 0 |

| Case | Tight reference E (Ha) | Strict plateau difference (Ha) | Max absolute energy difference (Ha) | Sampled gap min–max (Ha) |
|---|---:|---:|---:|---|
| si-gamma | -7.295530e+00 | 8.881784e-15 | 1.045279e-10 | 7.851461e-02 – 7.851461e-02 |
| si-k222s | -7.919423e+00 | 3.641532e-14 | 4.552607e-10 | 7.994936e-02 – 7.994936e-02 |
| si-skew | -7.250556e+00 | 8.881784e-16 | 6.180780e-10 | 3.146852e-03 – 3.146852e-03 |
| si-elongated | unavailable | unavailable | unavailable | unavailable |
| si-slab | unavailable | unavailable | unavailable | unavailable |
| si8-gamma | -3.132180e+01 | 1.776357e-14 | 2.245319e-09 | 1.691838e-02 – 1.691998e-02 |
| si16-gamma | -6.312924e+01 | 1.847411e-13 | 1.880501e-09 | 2.069342e-02 – 2.079582e-02 |
| si-low | -7.233756e+00 | 7.105427e-15 | 2.975309e-11 | 7.743359e-02 – 7.743359e-02 |
| material-si-20 | -8.462318e+00 | 5.329071e-15 | 2.593357e-10 | 4.335866e-02 – 4.335866e-02 |
| material-mgo-40 | -7.605088e+01 | 7.105427e-14 | 1.508766e-10 | 2.345956e-01 – 2.345956e-01 |

The supplied [dev08 histories](../../tests/reference/perf/p29a/si8_standard_histories_dev08.md)
are relayed evidence: 15 pure versus 27 MKL iterations, with repeated requested
tolerances at MKL iterations 9–13 and 15–17. Their first tabulated energies and
density residuals already differ; the pure table lacks the eigen residual needed
to independently verify the narrative's iteration-one 1.3e-16 agreement.
No new QE reference or backend causality claim follows from this study.

## Phase 2: parameter choice and held-out evaluation

The actual loop is `session.rs:advance`. The original ratchet only lowers the
requested Ritz-change tolerance when `0.1 * error / electrons` beats its previous
value. Thus a flat/rising error can preserve a loose eigensolve tolerance.
Tightening changes some plateau trajectories, but the work results below show
that this is not the sole cause of difficult SCF convergence.

The training/held-out split and candidate policy were frozen in
[`guard-plan.txt`](../../P29A-evidence/guard-plan.txt) before any guarded run.
Training: gamma, skew, elongated, Si8 and low-cutoff Si. Held out: k222s, slab,
Si16, material Si and MgO. The latter five never selected parameters.

The five-interval window requires **every** newer error to exceed half the
window-entry error; a transient rebound after genuine progress does not fire.
Baseline training replays of windows 3/5/8 flag 8/7/2 Si8 trajectories and all
10 elongated trajectories, with zero flags on gamma/skew/low for each window.
Five matches the prespecified stagnation interval; eight misses most short
plateaus. Candidates cap the next tolerance at 0.5 or 0.1 times the requested
tolerance at window entry, also retaining the ordinary ratchet and original
5e-14 Ha floor. The window-entry cap avoids multiplying a cut again on every
iteration of the same window. It needs no new mutable state or hot-loop allocation.

The 0.5 candidate does not change any Si8 trajectory; the 0.1 candidate reduces
Si8's maximum SCF count but raises its median work. **Neither passes the frozen
work-saving selection criterion.** The 0.1 candidate was therefore evaluated on
the held-out cases as a diagnostic of tail robustness, not as a selected speed
winner. It is retained only as an explicitly experimental opt-in switch for
review, with this negative result visible. There was no tuning after held-out
results. Ledger A21's rejected strict-schedule loosening is not repeated or
superseded; the strict path is unchanged.

```toml
[scf]
convergence_profile = "standard"
standard_stall_guard = "on" # experimental; default is "off"
```

Guard-on requires standard; unknown strings and strict+on reject. Guard-off
preserves the original arithmetic. Enabled provenance is carried in results,
options, checkpoint export/continuation and inspection. Legacy guard-off wire
fields are omitted and decode to false. Exact restart rejects a changed switch.

### Whole-set comparison (same ten perturbations and same tight references)

| Case | SCF off | SCF on | Davidson on (min/median/max) | H on (min/median/max) | Median H change | Max absolute energy difference on (Ha) |
|---|---|---|---|---|---:|---:|
| si-gamma | 9 / 9 / 9 | 9 / 9 / 9 | 27 / 27 / 27 | 180 / 180 / 180 | +0.00% | 1.045279e-10 |
| si-k222s | 9 / 9 / 9 | 9 / 9 / 9 | 106 / 106 / 106 | 731 / 731 / 731 | +0.00% | 4.552607e-10 |
| si-skew | 12 / 12 / 12 | 12 / 12 / 12 | 32 / 32 / 32 | 224 / 224 / 224 | +0.00% | 6.180780e-10 |
| si-elongated | unavailable | unavailable | 551 / 1519.5 / 1607 | 4626 / 11803.5 / 12956 | +128.51% | unavailable |
| si-slab | unavailable | unavailable | 1499 / 1683 / 1731 | 11394 / 13080 / 13550 | +141.33% | unavailable |
| si8-gamma | 20 / 25.5 / 32 | 20 / 24.5 / 27 | 86 / 101.5 / 113 | 958 / 1179 / 1262 | +0.99% | 3.164040e-09 |
| si16-gamma | 18 / 21.5 / 27 | 18 / 21.5 / 27 | 75 / 101.5 / 130 | 1580 / 1854.5 / 2389 | +0.00% | 1.880501e-09 |
| si-low | 9 / 9 / 9 | 9 / 9 / 9 | 24 / 24 / 24 | 176 / 176 / 176 | +0.00% | 2.975309e-11 |
| material-si-20 | 11 / 11 / 11 | 11 / 11 / 11 | 1042 / 1042 / 1042 | 6974 / 6974 / 6974 | +0.00% | 2.593357e-10 |
| material-mgo-40 | 13 / 13 / 13 | 13 / 13 / 13 | 609 / 609 / 609 | 5909 / 5909 / 5909 | +0.00% | 1.508766e-10 |

All 60 nonstagnating trajectories have identical complete scientific histories
with guard off/on (not just equal median work). Held-out Si16 has unchanged SCF
min/median/max and median H work, but its maximum H work increases 2318→2389.
No converged case becomes nonconvergent. The 20 difficult trajectories still do
not converge; no rescue or delivered-accuracy claim is made for them.

Across all 100 trajectories H work rises **275581→412309 vectors (+49.61%)**.
Across the 80 successful trajectories alone it changes **172409→172309 vectors
(-0.058%)**. This is a negative general performance result despite the Si8 tail
reduction. The largest measured standard energy difference is 3.164040e-9 Ha
(on), versus 2.245319e-9 Ha (off); all successful stopping estimates remain below
5e-9 Ha. These observations do not certify standard results or establish an
energy bound for untested systems. Repeated tightening to the existing floor
can waste substantial work on intrinsically difficult/nonconvergent trajectories.

## Two complete example histories

Both are guard-off Si8 trajectories. Units: energy/error/tolerance/eigen residual
in Ha, density residual in the engine's normalized L2 convention (Bohr^-3/2).
Their initial densities differ only by the documented tiny perturbations.

### Seed 1: 26 SCF iterations

| it | E | density residual | standard error | requested tol | eigen residual | Davidson | H |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | -30.6918914143 | 2.619278e-02 | 3.330154e+00 | 5.000000e-03 | 2.942412e-02 | 6 | 78 |
| 2 | -30.9653053937 | 1.665165e-02 | 1.364513e+00 | 5.000000e-03 | 5.527007e-02 | 2 | 40 |
| 3 | -31.0768218355 | 6.273439e-03 | 3.049387e-01 | 4.264102e-03 | 7.515511e-02 | 3 | 42 |
| 4 | -31.1686766727 | 3.566147e-03 | 7.733181e-02 | 9.529334e-04 | 5.851632e-02 | 3 | 42 |
| 5 | -31.1949597173 | 3.196216e-03 | 5.371799e-02 | 2.416619e-04 | 9.933641e-03 | 5 | 55 |
| 6 | -31.2011255506 | 2.743035e-03 | 3.912804e-02 | 1.678687e-04 | 1.236239e-02 | 2 | 40 |
| 7 | -31.1893779560 | 3.093362e-03 | 5.963296e-02 | 1.222751e-04 | 9.554474e-03 | 2 | 40 |
| 8 | -31.3076081800 | 5.627548e-03 | 1.065667e-01 | 1.222751e-04 | 9.584958e-03 | 8 | 59 |
| 9 | -31.3100017919 | 5.202293e-03 | 8.865902e-02 | 1.222751e-04 | 6.395445e-03 | 2 | 40 |
| 10 | -31.3142915167 | 4.393171e-03 | 6.092355e-02 | 1.222751e-04 | 5.376551e-03 | 2 | 40 |
| 11 | -31.3168067774 | 3.697554e-03 | 4.212841e-02 | 1.222751e-04 | 4.360524e-03 | 2 | 40 |
| 12 | -31.3150639148 | 4.175443e-03 | 5.480720e-02 | 1.222751e-04 | 3.563182e-03 | 2 | 40 |
| 13 | -31.3145372964 | 4.315528e-03 | 5.851321e-02 | 1.222751e-04 | 2.890901e-03 | 2 | 40 |
| 14 | -31.3170987214 | 2.778032e-03 | 2.576557e-02 | 1.222751e-04 | 3.546906e-03 | 2 | 40 |
| 15 | -31.3205580835 | 1.373559e-03 | 6.440458e-03 | 8.051741e-05 | 2.658603e-03 | 2 | 40 |
| 16 | -31.3215719272 | 8.964340e-04 | 2.223345e-03 | 2.012643e-05 | 1.888015e-03 | 2 | 40 |
| 17 | -31.3217548602 | 3.030610e-04 | 2.362114e-04 | 6.947954e-06 | 1.609341e-03 | 2 | 40 |
| 18 | -31.3217616610 | 1.669961e-04 | 1.030091e-04 | 7.381605e-07 | 7.785166e-04 | 8 | 56 |
| 19 | -31.3217973386 | 7.119988e-05 | 1.408916e-05 | 3.219034e-07 | 7.620945e-04 | 5 | 45 |
| 20 | -31.3218029466 | 3.341329e-05 | 2.978505e-06 | 4.402864e-08 | 1.277455e-04 | 11 | 54 |
| 21 | -31.3218027679 | 2.007688e-05 | 1.836173e-06 | 9.307828e-09 | 7.962312e-05 | 6 | 54 |
| 22 | -31.3218034416 | 7.709516e-06 | 2.166193e-07 | 5.738041e-09 | 1.100764e-04 | 4 | 42 |
| 23 | -31.3218035224 | 6.005767e-06 | 1.623130e-07 | 6.769353e-10 | 1.398504e-05 | 7 | 60 |
| 24 | -31.3218035775 | 3.255314e-06 | 3.553472e-08 | 5.072282e-10 | 2.794527e-05 | 3 | 41 |
| 25 | -31.3218035961 | 1.427134e-06 | 8.431216e-09 | 1.110460e-10 | 1.133734e-05 | 8 | 53 |
| 26 | -31.3218036003 | 6.838431e-07 | 1.649232e-09 | 2.634755e-11 | 3.476010e-06 | 7 | 50 |

### Seed 9: 32 SCF iterations

| it | E | density residual | standard error | requested tol | eigen residual | Davidson | H |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | -30.7122504292 | 2.616658e-02 | 3.307132e+00 | 5.000000e-03 | 2.942412e-02 | 6 | 78 |
| 2 | -30.9439222055 | 1.667529e-02 | 1.385094e+00 | 5.000000e-03 | 5.658302e-02 | 2 | 40 |
| 3 | -31.0665067232 | 6.295712e-03 | 3.310677e-01 | 4.328419e-03 | 7.246694e-02 | 3 | 42 |
| 4 | -31.1641736020 | 3.858353e-03 | 9.932647e-02 | 1.034587e-03 | 4.457856e-02 | 3 | 42 |
| 5 | -31.2508003055 | 5.024975e-03 | 9.068737e-02 | 3.103952e-04 | 1.062737e-02 | 11 | 63 |
| 6 | -31.2575749745 | 4.838714e-03 | 7.736278e-02 | 2.833980e-04 | 5.014986e-03 | 2 | 40 |
| 7 | -31.2494082154 | 4.673304e-03 | 8.445866e-02 | 2.417587e-04 | 3.545192e-03 | 2 | 40 |
| 8 | -31.2633572103 | 3.579747e-03 | 4.542612e-02 | 2.417587e-04 | 3.752826e-03 | 2 | 40 |
| 9 | -31.2500265611 | 2.933391e-03 | 5.163400e-02 | 1.419566e-04 | 7.086855e-03 | 2 | 40 |
| 10 | -31.2450428705 | 2.601459e-03 | 5.249274e-02 | 1.419566e-04 | 4.992836e-03 | 2 | 40 |
| 11 | -31.2496097732 | 2.686347e-03 | 4.828849e-02 | 1.419566e-04 | 7.269597e-03 | 2 | 40 |
| 12 | -31.2553880138 | 2.530900e-03 | 3.936548e-02 | 1.419566e-04 | 6.522435e-03 | 2 | 40 |
| 13 | -31.2665328648 | 1.627058e-03 | 1.855329e-02 | 1.230171e-04 | 6.164839e-03 | 2 | 40 |
| 14 | -31.3175210349 | 3.791572e-03 | 3.907214e-02 | 5.797903e-05 | 4.596680e-03 | 9 | 53 |
| 15 | -31.3176217227 | 3.744587e-03 | 3.809990e-02 | 5.797903e-05 | 2.093007e-03 | 2 | 40 |
| 16 | -31.3184463747 | 3.333217e-03 | 3.006720e-02 | 5.797903e-05 | 1.536640e-03 | 2 | 40 |
| 17 | -31.3195325144 | 2.744288e-03 | 2.054765e-02 | 5.797903e-05 | 1.313386e-03 | 2 | 40 |
| 18 | -31.3212852679 | 1.087685e-03 | 3.366620e-03 | 5.797903e-05 | 2.916193e-03 | 2 | 40 |
| 19 | -31.3214943645 | 8.495480e-04 | 2.109303e-03 | 1.052069e-05 | 1.538061e-03 | 2 | 40 |
| 20 | -31.3216829378 | 4.825478e-04 | 7.216230e-04 | 6.591573e-06 | 1.007650e-03 | 2 | 40 |
| 21 | -31.3217902651 | 1.607066e-04 | 7.711754e-05 | 2.255072e-06 | 1.027129e-03 | 3 | 43 |
| 22 | -31.3217995648 | 7.447524e-05 | 1.970725e-05 | 2.409923e-07 | 5.313183e-04 | 3 | 44 |
| 23 | -31.3218026357 | 4.216501e-05 | 5.340645e-06 | 6.158515e-08 | 1.539176e-04 | 8 | 46 |
| 24 | -31.3218033115 | 2.622582e-05 | 2.025883e-06 | 1.668951e-08 | 1.253824e-04 | 2 | 40 |
| 25 | -31.3218034550 | 6.309458e-06 | 1.365508e-07 | 6.330886e-09 | 1.161703e-04 | 3 | 41 |
| 26 | -31.3218034691 | 8.510263e-06 | 3.016725e-07 | 4.267214e-10 | 8.291621e-06 | 8 | 63 |
| 27 | -31.3218035020 | 7.536075e-06 | 2.222327e-07 | 4.267214e-10 | 8.180110e-06 | 2 | 40 |
| 28 | -31.3218035798 | 2.553126e-06 | 1.844808e-08 | 4.267214e-10 | 4.453786e-05 | 4 | 46 |
| 29 | -31.3218035908 | 2.424794e-06 | 2.363453e-08 | 5.765026e-11 | 4.826813e-06 | 9 | 61 |
| 30 | -31.3218035959 | 1.478148e-06 | 8.554647e-09 | 5.765026e-11 | 7.090246e-06 | 2 | 40 |
| 31 | -31.3218035983 | 1.112349e-06 | 5.611437e-09 | 2.673327e-11 | 2.665125e-06 | 9 | 60 |
| 32 | -31.3218035998 | 6.231286e-07 | 1.635502e-09 | 1.753574e-11 | 3.934698e-06 | 2 | 40 |

## Reproduction and verification

Run from the worktree root; only its ordinary `target/` is used. Real UPFs stay
external. `provenance.json` contains their SHA-256 hashes. All experiments use
release builds; the long sweep is ignored by ordinary debug gates.

```sh
P29A_SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF" \
  cargo test --release -p dft-engine --test standard_stall_trajectory sweep \
  -- --ignored --nocapture --test-threads=1 > P29A-evidence/off-si.log 2>&1
P29A_SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF" \
  DFT_NC_DIR="$PWD/../pseudo/pseudodojo/ex" \
  P29A_CASES=material-si-20,material-mgo-40 \
  cargo test --release -p dft-engine --test standard_stall_trajectory sweep \
  -- --ignored --nocapture --test-threads=1 > P29A-evidence/off-material.log 2>&1
# Repeat with P29A_START_RUN=2 P29A_GUARD=on, using the same cases/seeds.
```

The final harness emits complete terminal error histories. Early training logs
omitted the final failed record: elongated off/on were rerun with full histories
and override those incomplete runs in `final_analysis.py`; original logs remain
preserved. Earlier 0.5-candidate failure work is a lower bound and is not used in
the final comparison. The parser also handles Rust's test-name prefix on the
first output row and replaces, rather than concatenates, explicit reruns.

Additional plumbing-only edits: `session.rs` (the actual standard-profile loop,
result/checkpoint provenance), engine `checkpoint.rs` (option identity), CLI
`checkpoint_io.rs` and `checkpoint_inspect.rs` (additive wire/provenance).
No artifacts/dashboard/registry edits are retained. README/figures are regenerated
from existing stored validation; this experimental Mac study adds no performance
claim to that generated block. No remote server, MKL or mixing algorithm changed.

Final gate commands and exit codes are recorded in the completion report and
`P29A-evidence/commands.json`; `before-commands.json` / `final-commands.json`
retain the full five-case × five-thread × two-profile identity commands.
The default-path comparisons require zero energy difference, identical density
hash, unchanged SCF/work and strict certificate work. Standard stays uncertified.

Limitations: one machine/backend; ten deterministic perturbations per case;
no tight reference for elongated/slab; no new spectral-completeness proof;
no timing or QE speed acceptance. A future improvement must address the large
failure-case work penalty and validate independently, not tune on these held-out
cases and call them held out again. The guard is not recommended as a default.
