# Cutoff-ladder handoff: current-orbital indicator experiment

F-LADDER-HANDOFF, 2026-11-19, source checkout `7243075`.

**Scientific result: CLOSED, negative for the registered policy family.**
**Production decision: NO-GO.** The trained indicator policy reduces the
registered total-work proxy on **3 of 7 systems**, below the required four,
and does not win on either skewed fixture. No `BAANDI_HANDOFF_INDICATOR` switch
was added. The adopted ladder, final cutoff/grid, PP, final tolerance,
verification, and all production defaults are unchanged. All Rust additions
are under `#[cfg(test)]`; the long experiment is `#[ignore]`.

There are **126 handoff candidates** (7 systems, 3 seeds, 6 stops) and **105
coarse-configuration candidates** (the same systems/seeds, 5 configurations).
All 231 completed and passed fixed-point verification. Maximum absolute energy
difference from the same-geometry strict fine reference: **1.30466659698e-9 Ha**
for the handoff scan and **4.31099822507e-9 Ha** for the configuration scan,
both below the unchanged **5e-9 Ha total-cell standard class**. The latter
maximum is syn-slab, seed 603, fraction 0.4. These are internal strict-reference
comparisons, not new QE references or cutoff-convergence evidence.

**Delivery gate remains red:** the required three-crate check exits 1 because
three existing CLI pseudopotential-cache tests reject the required in-worktree
scratch directory. Build, fmt and clippy pass; the test runner reports 446
passed, 3 failed, 133 ignored across 48 binaries. The scientific experiment is
complete; a completely green required gate is not claimed. Details below.

## Fixtures and immutable fine problem

Geometry generation uses `scripts/fast_screen.py::displaced`, Python
`random.Random(seed)`, fractional amplitude +/-0.001, seeds 601/602/603.
This is not the old ladder bench's Rust xorshift geometry stream. Si cells
come from `benchmarks/cases`; MgO and synthetic cells reproduce the existing
ladder fixture definitions with the same displacement stream. Threads = 1,
Gamma, fixed closed-shell occupations, pure Rust on the Mac. No MPI or
threads-versus-OpenMP comparison was performed.

| System | Geometry / PP | Final cutoff (Ha) | Fine FFT grid | Bands / occupied |
| --- | --- | ---: | --- | --- |
| si8 | Cubic Si8, Si pz-vbc | 15 | 36x36x36 | 20 / 16 |
| si16 | Elongated Si16, Si pz-vbc | 15 | 36x36x72 | 36 / 32 |
| si32 | Si32 supercell, Si pz-vbc | 15 | 72x72x36 | 68 / 64 |
| si-skew | Skewed primitive Si, Si pz-vbc | 15 | 25x25x25 | 8 / 4 |
| mgo8 | Cubic rocksalt MgO, PseudoDojo Mg/O PBE + NLCC | 35 | 45x45x45 | 36 / 32 |
| syn-skew2 | Skewed two-species synthetic s/s/p + local PP | 4 | 9x8x12 | 4 / 2 |
| syn-slab | Slab-like two-species synthetic s/s/p + local PP | 4 | 8x8x20 | 4 / 2 |

MgO uses the existing NLCC density floor 1e-10; the 35 Ha cutoff is an
algorithmic fixture, not a demonstrated material cutoff limit. Synthetic PPs
are generated in memory by the existing `tests/common` helper. No UPF file
was copied into the repository.

The reference uses the same final problem, with `standard=None`, strict
mixing, density tolerance 1e-9, eigen residual tolerance 1e-10 Ha and energy
tolerance 1e-11 Ha, completed through the strict certificate path. Production
candidates keep `ScfOptions::default()` (standard threshold 5e-9 Ha and
fixed-point verification on). Their fine first tolerance is the existing
coarse-ratchet handoff rule; the fine convergence criteria never change.

## State and indicator validation

Existing `BAANDI_SNAPSHOT_DIR` capture records **solver-entry** state, not
solved output orbitals. Each candidate therefore also saves its actual solved
coarse orbitals and next mixed density, and the exact prolonged handoff state.
All are in task scratch. The run uses `RunConfig::from_pairs` for phase,
stage and Davidson instrumentation; the runner clears inherited `BAANDI_*`
overrides and enables only capture when requested. The fine and coarse run
configuration is shared as in the production ladder.

Orbital prolongation uses the existing exact Miller-index injection; density
uses the existing Fourier interpolation. The probe builds the fine H at the
**actual prolonged handoff density** and applies it to the occupied block
once, using the real Gamma operator. With occupied X and R = HX - X(X†HX),

    eta_high(s) = sum_i f_i sum_(G outside coarse PW space)
                  |R_i(G)|^2 / (|G|^2/2 + s).

The projection is outside the coarse **plane-wave** space, not merely outside
its occupied subspace. Since prolonged X is zero there, its high-G residual
is exactly the high-G part of HX. The denominator and shift are in Ha;
eta has units Ha. Shifts 0.5, 1 and 2 Ha use the same H block; the registered
policy uses 1 Ha. Full-complex +/-G coefficients are summed once each, with
occupations f_i (two electrons per occupied orbital); there is no additional
half-sphere multiplicity factor. The low-G diagnostic uses the occupied
subspace residual. Empty-state defects are not separately probed, although
all empty-band solve work is included in the final recovery counts.

This is a heuristic, **not a certified LDA/PBE error bound**. It is never
compared directly with the density residual in Bohr^-3/2. Calibration is
against completed fine/verification work, with strict energy accuracy checked
on every completed run; no handoff energy-bound claim is made.

Validation actually executed:

- Analytic residual fixture: known weighted high-G value, zero low-G residual,
  zero occupation and empty complement checks pass.
- Each of 21 geometries compares the real packed probe HX against an
  independent complex H block, maximum-entry tolerance 1e-11 Ha; all pass.
- Each baseline candidate compares against ordinary `run_ladder`: final
  energy bits, delivered density, coarse H and fine map count match exactly.
- All 21 coarse tolerance/error trajectory prefixes match across stopping
  arms, validating causal replay of the sampled checkpoints.
- All 21 capture-on/off baseline pairs agree in final/reference energy,
  H/block/delta counts, map counts and coarse tolerance/error history.
- An independent Python reader validates **2,202** snapshot manifests, density
  lengths, rank-file magic/header/layout and payload lengths, plus terminal
  coarse/prolonged state sizes. No snapshot write failures remain in the data.
- `fast_screen.py`, Si8 seed 601, one-core A/A: energy
  -31.321776563202498 Ha, coarse/fine maps 5/11, recorded fine H 430,
  verification H 101, coarse solve H 206 plus atomic-start H 32. Both arms
  match the experiment. Density SHA-256 is
  `43f2a33d8b0d0a5a90129076353d8e4ad26b2305695c3d834c13b7598fb7342b`.
  The harness correctly labels wall unjudged with only one pair.

## Accounting and probe cost

The fixed-stop scan uses fraction 0.5 and residual levels 1e-2, 3e-3, 1e-3,
3e-4, 1e-4 (adopted) and 3e-5. Duplicate stopping maps are retained in the
candidate table but deduplicated for calibration and probe charging.

Hc/Hf/Hv are actual full-H **column counts**, including coarse atomic start,
all ordinary/redo fine solves and **every** failed or passed verification
attempt. C counts full-H block calls. D counts delta-potential block updates
of cached products, separately from full H. The counters reconcile with the
history and start records in every row. Atomic-start H is already in the
real-H counters and is not added twice. Handoff interpolation and import
checks add no H calls; their time is inside the capture-disabled run wall.
The reported handoff timer covers prolongation/checks; subsequent session
import belongs to fine setup in that wall.

A delta update uses full inverse/forward transforms for unkept columns and
only forward transforms for kept columns. Let alpha be coarse grid points /
fine grid points. The raw JSON/CSV records raw total H and full-H-equivalent
Hf + Hv + alpha*Hc. The table's **W** additionally charges unkept delta
columns at one transform-column equivalent and kept delta columns at one
half, weighted by alpha at the coarse level. Delta block counts are checked
against the independent reuse counters. This is a work proxy, not a wall
model: density/potential/audit and import have real costs even when they add
zero H applications. Maps and measured coarse phase times are reported too.
Grid-size weighting also does not model FFT factorization costs.

Every requested probe costs one occupied fine H block (2 to 64 columns in
these fixtures) and one block call, including probes whose criterion fails.
No credit is assumed for reusing its HX in the eventual fine solve. The
policy cost tables add that entire charge. Probe H-only and complete offline
diagnostic times are tabulated below. The diagnostic total includes separate
problem/operator construction, the fine potential, three shifts and low-G
projection work; it is not a timing of an optimized production controller.
Even the optimistic one-block charge fails the four-system work gate, so
there is no evidence here that adding setup costs would justify adoption.

Capture-enabled total wall includes snapshot I/O and is **not** used as a
speed result. The separate configuration batch has capture off. Its three
geometries per arm, fixed arm order and shared-Mac load support observations
only, not a speed-adoption decision. No dev08 claim is inferred from them.

## Registered policy, calibration and decision

Before completion of the main matrix, the policy family was fixed: probe at
distinct states reached by coarse residual triggers 1e-3 and 3e-4, only before
the ordinary 1e-4 handoff. Accept a probed state when its coarse Coulomb
residual energy <= q*eta_high(1 Ha); otherwise continue, falling back to the
adopted handoff. Charge all probes. Train q in {0.001, 0.01, 0.1, 1} on seed
601 using equally weighted mean W/baseline over systems; evaluate on seeds
602 and 603. This does not claim the unmeasured intervening maps are optimal.

The best q is **0.1**, but even its training mean ratio is **1.004606**.
The other q ratios are 1.014994, 1.007892 and 1.059802 (q=0.001, 0.01, 1).
A separate globally fixed threshold selected on the same training seed is
1e-2, with training mean ratio 0.904042. It is compared on the same held-out
seeds below; it is not uniformly better and is not adopted.

Held-out median paired indicator W/baseline: MgO 0.906958, Si8 0.980086,
synthetic slab 0.853172 (wins); Si16 1.074884, Si32 1.074247, synthetic skew
1.000143 (losses); Si-skew 1.000000 (no early probe, tie). Thus **3/7**, with
no skewed-system win: **NO-GO** under the user's four-system requirement.
No extra minimum percentage threshold was imposed.

The following descriptive calibration uses unique (seed, coarse-map) states.
Correlations pool dependent trajectory points and are not statistical tests.
The high-G indicator changes little while downstream work changes materially;
its sign of association with work is not consistent across fixtures. The
result closes this tested kinetic-weighted, sparse-trigger policy family, not
every possible current-state estimator or future richer response model.

| System | Unique states | eta high range (Ha) | Relative span | Remaining fine+verification W range | Pearson high/work | Pearson low/work |
| --- | ---: | --- | ---: | --- | ---: | ---: |
| mgo8 | 15 | 3.53723–3.64853 | 3.108% | 755–1116 | -0.5003 | -0.8077 |
| si-skew | 12 | 0.00840611–0.008465848 | 0.709% | 178–219 | -0.8081 | -0.7180 |
| si16 | 17 | 0.04951565–0.05039066 | 1.739% | 1110–2166 | -0.6419 | 0.5915 |
| si32 | 12 | 0.09548765–0.09764389 | 2.214% | 2278–3968 | -0.4446 | 0.3923 |
| si8 | 13 | 0.02566121–0.02581079 | 0.583% | 577–696 | 0.1885 | 0.1872 |
| syn-skew2 | 15 | 0.01334169–0.01345856 | 0.876% | 77–92 | 0.0130 | 0.0174 |
| syn-slab | 15 | 0.03307941–0.03476324 | 5.079% | 89–111 | -0.2550 | -0.2024 |

## Cheaper coarse configurations

The supplementary scan holds stop=1e-4 and compares the baseline with
fractions 0.35/0.4, the smallest grid that still represents the coarse density
sphere (not an aliased or invalid grid), and a coarse-only loose Ritz schedule.
The loose schedule multiplies each ratchet proposal by four, retaining the
previous-tolerance ceiling and the 0.005 Ha cap. The existing first-step redo
is retained. Its small offline driver first proves factor=1 density/orbital
bit identity and H-count identity against `coarse_stage` for **all 21**
geometries before factor=4 is measured. Fine tolerances and verification stay
on the original session path. This driver has no production hook.

- Fraction 0.35 reduces median paired W on MgO, Si-skew, Si16 and syn-slab,
  but all four require more fine work. Si8, Si32 and syn-skew2 regress.
- Fraction 0.4 reduces W on MgO, Si-skew and Si8 without a median fine-work
  increase; only Si-skew/Si8 avoid additional fine work in every seed. It
  regresses on the other four systems.
- Minimum valid grids reduce W on five systems; four have no median fine
  increase (Si-skew, Si16, Si32, syn-slab). This is not a universal time win:
  MgO's observed wall ratio is 1.361, syn-skew2 1.144 and syn-slab 1.063.
  Smaller nonsmooth FFT dimensions can cost more despite fewer grid points.
  Si8 already has the minimum grid and is unchanged numerically.
- Loose coarse tolerance improves W and fine work on MgO and Si16 only;
  Si8, Si32 and both skew/slab synthetic systems regress in W. Reducing the
  coarse solve alone is insufficient (the BT3-COARSE lesson persists).

No configuration default is changed. Some candidates merit a separately
controlled dev08 screen, but none establishes a generally faster default here.
All per-seed rows, complete tolerance schedules, phase seconds and counts are
retained in JSON/CSV. Summary columns below are independent medians; summing
median components need not equal the median of their totals.

### Coarse configuration screen (medians over seeds 601–603)

W includes all fine work and verification; ratios are medians of paired ratios. Fine ratio includes verification and delta updates. Wall ratios are loaded-Mac observations, not acceptance evidence. No capture I/O is in these walls.

| System | Arm (fraction:stop:mode) | Coarse grid | Coarse/fine/verify maps | Hc/Hf/Hv | Total H | C/D | W | W/base | Fine/base | Wall/base (O) |
| --- | --- | --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: |
| mgo8 | 0.35:1e-4:default | 25x25x25 | 8/8/1 | 579/845/68 | 1501 | 50/15 | 1198.86 | 0.9269 | 1.0056 | 1.038 |
| mgo8 | 0.4:1e-4:default | 27x27x27 | 8/8/1 | 597/830/70 | 1502 | 52/15 | 1215.57 | 0.9387 | 0.9878 | 1.044 |
| mgo8 | 0.5:1e-4:default | 30x30x30 | 8/8/1 | 624/827/68 | 1531 | 50/15 | 1293.37 | 1.0000 | 1.0000 | 1.000 |
| mgo8 | 0.5:1e-4:loose4 | 30x30x30 | 8/8/1 | 547/819/71 | 1436 | 48/15 | 1253.26 | 0.9695 | 0.9659 | 0.953 |
| mgo8 | 0.5:1e-4:minimal-grid | 29x29x29 | 9/8/1 | 683/851/73 | 1608 | 53/16 | 1310.62 | 1.0133 | 1.0169 | 1.361 |
| si-skew | 0.35:1e-4:default | 15x15x15 | 13/8/3 | 136/137/35 | 308 | 51/22 | 276.93 | 0.9732 | 1.0691 | 0.843 |
| si-skew | 0.4:1e-4:default | 18x15x15 | 13/7/2 | 137/120/21 | 278 | 45/20 | 243.17 | 0.8546 | 0.8710 | 0.776 |
| si-skew | 0.5:1e-4:default | 18x18x18 | 11/8/3 | 122/121/36 | 278 | 48/20 | 284.56 | 1.0000 | 1.0000 | 1.000 |
| si-skew | 0.5:1e-4:loose4 | 18x18x18 | 12/8/4 | 127/122/46 | 295 | 47/22 | 306.04 | 1.0703 | 1.0734 | 0.905 |
| si-skew | 0.5:1e-4:minimal-grid | 17x17x17 | 12/8/1 | 134/130/14 | 278 | 47/19 | 254.89 | 0.8957 | 0.8848 | 0.927 |
| si16 | 0.35:1e-4:default | 24x24x45 | 10/9/2 | 697/930/196 | 1774 | 87/19 | 1589.61 | 0.9605 | 1.0098 | 0.829 |
| si16 | 0.4:1e-4:default | 24x24x45 | 10/11/3 | 689/1011/287 | 1989 | 96/22 | 1801.39 | 1.1025 | 1.1818 | 1.017 |
| si16 | 0.5:1e-4:default | 25x25x54 | 10/10/2 | 688/916/195 | 1788 | 88/20 | 1633.95 | 1.0000 | 1.0000 | 1.000 |
| si16 | 0.5:1e-4:loose4 | 25x25x54 | 7/9/1 | 502/880/94 | 1501 | 75/15 | 1403.97 | 0.8592 | 0.8845 | 0.815 |
| si16 | 0.5:1e-4:minimal-grid | 25x25x51 | 10/10/2 | 688/916/195 | 1788 | 88/20 | 1616.51 | 0.9893 | 1.0000 | 0.961 |
| si32 | 0.35:1e-4:default | 45x45x24 | 7/11/2 | 1135/2042/354 | 3531 | 82/18 | 3199.60 | 1.1261 | 1.2339 | 1.033 |
| si32 | 0.4:1e-4:default | 45x45x24 | 6/12/3 | 1027/1843/549 | 3421 | 86/19 | 3216.02 | 1.0529 | 1.1229 | 1.026 |
| si32 | 0.5:1e-4:default | 54x54x25 | 7/10/1 | 1145/1735/201 | 3077 | 79/16 | 2823.64 | 1.0000 | 1.0000 | 1.000 |
| si32 | 0.5:1e-4:loose4 | 54x54x25 | 7/11/2 | 956/1848/390 | 3165 | 89/18 | 3098.81 | 1.0910 | 1.1449 | 1.108 |
| si32 | 0.5:1e-4:minimal-grid | 51x51x25 | 7/10/1 | 1145/1735/201 | 3077 | 79/16 | 2766.21 | 0.9797 | 1.0000 | 0.962 |
| si8 | 0.35:1e-4:default | 24x24x24 | 5/9/2 | 235/442/99 | 780 | 72/14 | 747.44 | 1.0321 | 1.0798 | 0.932 |
| si8 | 0.4:1e-4:default | 24x24x24 | 5/9/2 | 237/408/92 | 738 | 63/14 | 706.44 | 0.9335 | 0.9422 | 0.829 |
| si8 | 0.5:1e-4:default | 25x25x25 | 5/9/2 | 238/423/100 | 763 | 69/15 | 738.78 | 1.0000 | 1.0000 | 1.000 |
| si8 | 0.5:1e-4:loose4 | 25x25x25 | 6/9/2 | 232/434/108 | 802 | 75/16 | 799.13 | 1.0454 | 1.0489 | 1.019 |
| si8 | 0.5:1e-4:minimal-grid | 25x25x25 | 5/9/2 | 238/423/100 | 763 | 69/15 | 738.78 | 1.0000 | 1.0000 | 0.990 |
| syn-skew2 | 0.35:1e-4:default | 5x5x8 | 15/6/1 | 86/79/7 | 161 | 55/20 | 126.69 | 1.1020 | 1.2683 | 0.957 |
| syn-skew2 | 0.4:1e-4:default | 5x5x8 | 15/6/1 | 86/79/7 | 161 | 55/20 | 126.69 | 1.1020 | 1.2683 | 0.974 |
| syn-skew2 | 0.5:1e-4:default | 8x5x8 | 10/5/1 | 62/61/5 | 128 | 44/14 | 113.96 | 1.0000 | 1.0000 | 1.000 |
| syn-skew2 | 0.5:1e-4:loose4 | 8x5x8 | 10/7/1 | 58/61/5 | 124 | 40/16 | 118.48 | 1.0396 | 1.0741 | 0.984 |
| syn-skew2 | 0.5:1e-4:minimal-grid | 7x5x7 | 10/6/1 | 61/64/5 | 132 | 46/15 | 111.95 | 0.9824 | 1.0741 | 1.144 |
| syn-slab | 0.35:1e-4:default | 5x5x12 | 7/6/1 | 35/85/6 | 126 | 42/12 | 121.42 | 0.9865 | 1.0680 | 0.899 |
| syn-slab | 0.4:1e-4:default | 5x5x12 | 9/7/1 | 61/90/8 | 159 | 56/15 | 138.92 | 1.1195 | 1.1553 | 1.057 |
| syn-slab | 0.5:1e-4:default | 5x5x15 | 8/6/1 | 51/78/7 | 136 | 47/13 | 124.09 | 1.0000 | 1.0000 | 1.000 |
| syn-slab | 0.5:1e-4:loose4 | 5x5x15 | 11/6/1 | 60/80/6 | 146 | 48/16 | 130.37 | 1.0506 | 1.0097 | 0.987 |
| syn-slab | 0.5:1e-4:minimal-grid | 5x5x13 | 8/6/1 | 51/78/7 | 136 | 46/13 | 120.28 | 0.9764 | 1.0000 | 1.063 |

### Coarse phase and probe observations (seconds)

Baseline, capture off, median of three different geometries. Probe total includes independent problem construction and all three shifts plus low-G diagnostics; H-only is the single occupied-block apply. Neither is a measured production-controller wall.

| System | Maps (601/602/603) | Solve | Density | Potential | Energy | Mixing | Session | Handoff | Probe H columns | Probe H-only | Probe total |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| mgo8 | 8/8/8 | 0.117671 | 0.017930 | 0.072097 | 0.055243 | 0.002813 | 0.266511 | 0.001305 | 32 | 0.011811 | 0.102477 |
| si-skew | 11/11/11 | 0.005657 | 0.000887 | 0.003848 | 0.002061 | 0.000976 | 0.013607 | 0.000299 | 4 | 0.000306 | 0.010151 |
| si16 | 10/10/6 | 0.216233 | 0.035338 | 0.050215 | 0.011505 | 0.005385 | 0.322131 | 0.001975 | 32 | 0.013292 | 0.078789 |
| si32 | 7/7/6 | 0.721192 | 0.096858 | 0.201439 | 0.015735 | 0.005959 | 1.038415 | 0.003362 | 64 | 0.050688 | 0.402147 |
| si8 | 5/6/5 | 0.030718 | 0.004766 | 0.011976 | 0.002932 | 0.000943 | 0.051116 | 0.000807 | 16 | 0.003375 | 0.016466 |
| syn-skew2 | 10/10/10 | 0.000288 | 0.000023 | 0.000075 | 0.000073 | 0.000062 | 0.000565 | 0.000014 | 2 | 0.000003 | 0.001256 |
| syn-slab | 8/8/8 | 0.000284 | 0.000020 | 0.000074 | 0.000070 | 0.000048 | 0.000543 | 0.000022 | 2 | 0.000007 | 0.000565 |

### Coarse tolerance schedules (Ha)

All baseline schedules start at the standard Ritz cap 0.005 Ha; later entries are the actual requested tolerances, not a fitted schedule.

- si8 seed 601: 0.005, 0.005, 0.00179937, 1.98896e-05, 5.72067e-06
- si8 seed 602: 0.005, 0.005, 0.00179933, 1.9879e-05, 5.72459e-06, 1.39071e-07
- si8 seed 603: 0.005, 0.005, 0.00179932, 1.9881e-05, 5.71697e-06
- si16 seed 601: 0.005, 0.005, 0.00161958, 7.54073e-05, 7.54073e-05, 7.54073e-05, 3.24184e-05, 1.67805e-05, 1.12951e-05, 7.89652e-06
- si16 seed 602: 0.005, 0.005, 0.0016206, 0.000105664, 0.000105664, 0.000105664, 5.07371e-05, 2.51443e-05, 1.52824e-05, 5.57362e-06
- si16 seed 603: 0.005, 0.005, 0.00161642, 3.70442e-05, 3.70442e-05, 3.34617e-05
- si32 seed 601: 0.005, 0.005, 0.00153223, 4.76577e-05, 4.76577e-05, 4.54906e-05, 3.15299e-06
- si32 seed 602: 0.005, 0.005, 0.00153373, 5.25561e-05, 5.25561e-05, 4.98353e-05, 4.02366e-06
- si32 seed 603: 0.005, 0.005, 0.00153228, 4.95775e-05, 4.95775e-05, 4.23151e-05
- si-skew seed 601: 0.005, 0.005, 0.00239069, 0.00198872, 2.01429e-05, 1.46206e-05, 5.91921e-06, 5.28841e-06, 4.35332e-06, 4.35332e-06, 3.69465e-06
- si-skew seed 602: 0.005, 0.005, 0.00239098, 0.00198665, 2.0152e-05, 1.46251e-05, 5.92833e-06, 5.29549e-06, 4.35715e-06, 4.35715e-06, 3.69832e-06
- si-skew seed 603: 0.005, 0.005, 0.00239135, 0.00197721, 2.01354e-05, 1.45698e-05, 5.89118e-06, 5.26243e-06, 4.34906e-06, 4.34906e-06, 3.66747e-06
- mgo8 seed 601: 0.005, 0.005, 0.005, 0.005, 0.000351496, 0.000253043, 2.16055e-05, 3.90679e-07
- mgo8 seed 602: 0.005, 0.005, 0.005, 0.005, 0.000380389, 0.000251043, 1.77411e-05, 2.87386e-07
- mgo8 seed 603: 0.005, 0.005, 0.005, 0.005, 0.000359424, 0.000253674, 2.0522e-05, 3.39631e-07
- syn-skew2 seed 601: 0.005, 0.00382066, 0.000817789, 8.49313e-06, 8.49313e-06, 1.81814e-06, 6.23137e-07, 2.53341e-08, 2.29717e-08, 1.77119e-08
- syn-skew2 seed 602: 0.005, 0.00382121, 0.000817845, 8.36642e-06, 8.36642e-06, 2.01864e-06, 6.64366e-07, 2.19512e-08, 2.19512e-08, 1.8011e-08
- syn-skew2 seed 603: 0.005, 0.00382148, 0.000816913, 8.36522e-06, 8.36522e-06, 2.04145e-06, 6.71239e-07, 2.19242e-08, 2.19242e-08, 1.7727e-08
- syn-slab seed 601: 0.005, 0.005, 0.00282102, 0.00038636, 0.000261742, 1.71617e-06, 1.14436e-07, 1.29359e-08
- syn-slab seed 602: 0.005, 0.005, 0.00282327, 0.000387182, 0.000262533, 1.74777e-06, 1.13035e-07, 1.31819e-08
- syn-slab seed 603: 0.005, 0.005, 0.00282411, 0.000387472, 0.00026281, 1.7975e-06, 1.88289e-07, 2.22165e-08

## Reproduction and evidence

Use only worktree scratch (the user's correction resolves the initial sandbox
blocker):

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export TMPDIR_990="$PWD/target/scratch"
mkdir -p "$TMPDIR_990"
export TMPDIR="$TMPDIR_990"
cargo test --release -p dft-engine --lib --no-run
```

Use the test executable printed by cargo as `--binary` (in this checkout,
`target/release/deps/dft_engine-5a2965bdd204a9df`). The runner refuses a nonempty
output directory to prevent stale evidence. The recorded commands used:

```sh
python3 -I scripts/ladder_handoff.py --binary target/release/deps/dft_engine-5a2965bdd204a9df --out target/scratch/handoff-main --capture
python3 -I scripts/ladder_handoff.py --binary target/release/deps/dft_engine-5a2965bdd204a9df --out target/scratch/handoff-config --arms 0.5:1e-4:default,0.35:1e-4:default,0.4:1e-4:default,0.5:1e-4:minimal-grid,0.5:1e-4:loose4
python3 -I scripts/ladder_handoff.py analyze --rows target/scratch/handoff-main/rows.json --out target/scratch/handoff-analysis
python3 -I scripts/ladder_handoff.py analyze-config --rows target/scratch/handoff-config/rows.json --out target/scratch/handoff-analysis
```

All four exit 0. For reruns, choose fresh output directories and point the
analysis commands at them. Inputs and PP hashes live in each job's `input.json`;
its exact fixture text, stdout, command/exit record and arm measurements are
retained. `rows.json`/`rows.csv` are complete records. Snapshots and terminal
states are internal binary exchange, not new public volumetric formats.
Archived source files in each batch match the source hashes in its provenance.
These scratch artifacts are untracked; the complete handoff count table below
preserves the key evidence in the reviewable document.

| Batch | Binary SHA-256 | Engine measurement source SHA-256 | Rows SHA-256 |
| --- | --- | --- | --- |
| handoff-main | `39d6bfcbc3dbdfdd5a71032951efa249e213c2515cbfff418e09b6a020880ded` | `995fc0e4a022ae6a4398400d9bac9669c1cc59316030c8f93dc11c0f4599bb75` | `47134317c7acba6849b36f165205ced0d9db4394f6dfcfde216919945bd8e9df` |
| handoff-config | `0f5d6888bfc1be73a5b1d8be8bec1bd292f609a800d1bfb3416aaf09bf42d0e4` | `ae8eb05a2f87449a5b221ef87ad2d03f39f8556683d7bfb3e7beb53f5477c761` | `02c956c69bf4cc16b08fb6c125194d35c2a06b1e25f00399b8f6e7bef00b09eb` |

External PP SHA-256 (real files were read, never committed):

- `Si.pz-vbc.UPF`: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
- `Mg.upf`: `15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14`.
- `O.upf`: `2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11`.

The later lint-only change reuses the parent's existing synthetic fixture
module instead of loading it twice; it changes no numerical method. Final
source/binary hashes for that smoke are in `handoff-final-smoke/provenance.json`.

## Verification and limits

- Scope isolation: `scripts/check.sh fast dft-engine dft-numerics` exits **0**;
  all four gates pass, **315 passed, 0 failed, 93 ignored** in 35 binaries.
  This does not replace or waive the required three-crate gate.
- Required `scripts/check.sh fast dft-engine dft-numerics dft-cli`: **exit 1**.
  Build/fmt/clippy **PASS**; tests **446 passed, 3 failed, 133 ignored** across
  48 binaries. `ps` is also denied for the advisory load header; the cargo
  gates still run. Logs: `target/scratch/check-final.log`.
- Existing failures: CLI `pseudo::cache_resolution_and_no_network_side_effects`,
  `pseudo::cached_fetch_checks_hash_and_provenance`, and
  `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`.
  Their temp caches are below `TMPDIR`; `pseudo.rs::outside_repository` rejects
  any path under the compiled repository or a `.git` ancestor. They invoke
  cache operations before any SCF. No CLI source/tests are changed by this task.
  Fixing that fixture/location contract is outside OWNED_FILES (CLI run.rs is
  authorized only for switch counts). No `/tmp`, external scratch, policy
  bypass or wider permission was used to turn this gate green.
- The first check also caught a duplicate-module lint introduced by the
  experiment; reusing `scf::pool_fixture` fixes it, and final clippy passes.
- Early smoke captures with a missing snapshot subdirectory and an initial
  atomic-start/delta counter-accounting error were discarded. The accepted
  batches create/validate capture directories and reconcile every counter.
  The initial Rust build also caught use of nonexistent real-matrix column
  accessors; it was fixed before accepted measurements. No failed smoke is
  used in the tables.
- Final analytic indicator test: `cargo test --release -p dft-engine --lib
  handoff_experiment::high_g_indicator -- --nocapture` exits **0** (one test).
  Final Si8 seed-601 default/loose4 smoke exits **0**; final/reference energies,
  H/block counts, map counts, W and tolerance schedules match the pre-lint-fix
  binary exactly. Logs: `indicator-final.log`, `final-smoke.log` in scratch.
  The required full gate is not relabelled as passing.
- No production implementation or switch was adopted, so no switch-count
  update or new README validation-result type is introduced. Dashboard/session
  publication remains with the master session; those files are outside this
  ticket's owned outputs.

Existing LADDER-REEVAL-DEV08 supports keeping the ladder. D1-LADDER measured
fixed stops and the atomic-density proxy; this experiment measures **current**
solved orbitals instead. BT3-COARSE warned that fine recovery can repay coarse
savings; the present full-cost counts confirm that risk. No old ledger timing
is reused as a measurement of this indicator, and no QE source was read.

## Full handoff scan and held-out policy costs

### Complete candidate costs (per geometry)

Hc includes atomic start; Hf includes redo; Hv includes every verification attempt. C is full H block calls; D is delta-potential block calls. W is the grid-weighted FFT-column proxy including delta updates, before probes. eta is shift=1 Ha. All final energies passed verification.

| System | Seed | Stop | Coarse/fine/verify maps | Hc/Hf/Hv | Total H | C/D | W | eta (Ha) | Error (Ha) |
| --- | ---: | ---: | --- | --- | ---: | --- | ---: | ---: | ---: |
| mgo8 | 601 | 1e-02 | 4/8/1 | 373/522/73 | 968 | 36/11 | 883.30 | 3.64853 | 1.05e-09 |
| mgo8 | 601 | 3e-03 | 6/9/1 | 506/652/69 | 1227 | 41/14 | 1080.56 | 3.53723 | 1.10e-11 |
| mgo8 | 601 | 1e-03 | 7/8/1 | 560/731/67 | 1358 | 46/14 | 1159.48 | 3.57469 | 6.83e-11 |
| mgo8 | 601 | 3e-04 | 7/8/1 | 560/731/67 | 1358 | 46/14 | 1159.48 | 3.57469 | 6.83e-11 |
| mgo8 | 601 | 1e-04 | 8/8/1 | 620/827/66 | 1513 | 50/15 | 1278.19 | 3.58026 | 3.26e-11 |
| mgo8 | 601 | 3e-05 | 9/8/1 | 673/871/73 | 1617 | 56/16 | 1350.81 | 3.58081 | 5.46e-11 |
| mgo8 | 602 | 1e-02 | 4/8/1 | 373/541/72 | 986 | 38/11 | 901.30 | 3.648 | 1.04e-09 |
| mgo8 | 602 | 3e-03 | 6/9/1 | 507/678/53 | 1238 | 42/14 | 1090.85 | 3.5429 | 8.47e-12 |
| mgo8 | 602 | 1e-03 | 7/8/1 | 561/736/67 | 1364 | 46/14 | 1164.78 | 3.57824 | 7.47e-11 |
| mgo8 | 602 | 3e-04 | 7/8/1 | 561/736/67 | 1364 | 46/14 | 1164.78 | 3.57824 | 7.47e-11 |
| mgo8 | 602 | 1e-04 | 8/8/2 | 626/788/148 | 1562 | 51/16 | 1342.96 | 3.58272 | 2.61e-11 |
| mgo8 | 602 | 3e-05 | 9/8/1 | 687/864/68 | 1619 | 52/16 | 1342.96 | 3.58086 | 6.06e-11 |
| mgo8 | 603 | 1e-02 | 4/8/1 | 373/540/72 | 985 | 37/11 | 900.30 | 3.64839 | 1.22e-09 |
| mgo8 | 603 | 3e-03 | 6/9/1 | 506/654/71 | 1231 | 41/14 | 1084.56 | 3.53885 | 8.30e-12 |
| mgo8 | 603 | 1e-03 | 7/8/1 | 560/733/67 | 1360 | 46/14 | 1161.48 | 3.57584 | 7.22e-11 |
| mgo8 | 603 | 3e-04 | 7/8/1 | 560/733/67 | 1360 | 46/14 | 1161.48 | 3.57584 | 7.22e-11 |
| mgo8 | 603 | 1e-04 | 8/8/1 | 624/839/68 | 1531 | 49/15 | 1293.37 | 3.58104 | 4.15e-11 |
| mgo8 | 603 | 3e-05 | 9/8/1 | 676/870/74 | 1620 | 56/16 | 1351.70 | 3.58088 | 7.36e-11 |
| si-skew | 601 | 1e-02 | 4/8/2 | 62/102/22 | 186 | 30/12 | 207.86 | 0.00843896 | 8.12e-11 |
| si-skew | 601 | 3e-03 | 5/8/2 | 71/102/23 | 196 | 32/13 | 214.46 | 0.00845459 | 1.26e-11 |
| si-skew | 601 | 1e-03 | 11/8/3 | 121/121/36 | 278 | 48/20 | 284.56 | 0.00840611 | 2.64e-12 |
| si-skew | 601 | 3e-04 | 11/8/3 | 121/121/36 | 278 | 48/20 | 284.56 | 0.00840611 | 2.64e-12 |
| si-skew | 601 | 1e-04 | 11/8/3 | 121/121/36 | 278 | 48/20 | 284.56 | 0.00840611 | 2.64e-12 |
| si-skew | 601 | 3e-05 | 13/7/1 | 142/135/12 | 289 | 51/19 | 268.88 | 0.00840897 | 7.52e-11 |
| si-skew | 602 | 1e-02 | 4/8/2 | 62/102/22 | 186 | 30/12 | 207.86 | 0.00844284 | 8.14e-11 |
| si-skew | 602 | 3e-03 | 5/8/2 | 71/102/23 | 196 | 32/13 | 214.46 | 0.00845878 | 1.25e-11 |
| si-skew | 602 | 1e-03 | 11/8/3 | 122/122/36 | 280 | 48/20 | 285.93 | 0.00841487 | 2.66e-12 |
| si-skew | 602 | 3e-04 | 11/8/3 | 122/122/36 | 280 | 48/20 | 285.93 | 0.00841487 | 2.66e-12 |
| si-skew | 602 | 1e-04 | 11/8/3 | 122/122/36 | 280 | 48/20 | 285.93 | 0.00841487 | 2.66e-12 |
| si-skew | 602 | 3e-05 | 14/7/3 | 150/121/35 | 306 | 50/22 | 295.10 | 0.00841325 | 1.60e-12 |
| si-skew | 603 | 1e-02 | 4/8/2 | 62/102/22 | 186 | 30/12 | 207.86 | 0.00845005 | 7.65e-11 |
| si-skew | 603 | 3e-03 | 5/8/2 | 71/102/23 | 196 | 32/13 | 214.46 | 0.00846585 | 1.23e-11 |
| si-skew | 603 | 1e-03 | 11/8/3 | 122/120/34 | 276 | 46/20 | 281.93 | 0.0084217 | 1.46e-12 |
| si-skew | 603 | 3e-04 | 11/8/3 | 122/120/34 | 276 | 46/20 | 281.93 | 0.0084217 | 1.46e-12 |
| si-skew | 603 | 1e-04 | 11/8/3 | 122/120/34 | 276 | 46/20 | 281.93 | 0.0084217 | 1.46e-12 |
| si-skew | 603 | 3e-05 | 14/7/3 | 150/132/33 | 315 | 53/22 | 304.10 | 0.00842026 | 9.43e-13 |
| si16 | 601 | 1e-02 | 2/14/2 | 227/1126/194 | 1547 | 85/16 | 1709.34 | 0.0496874 | 2.70e-11 |
| si16 | 601 | 3e-03 | 3/19/3 | 290/1484/262 | 2036 | 98/23 | 2285.36 | 0.0500596 | 4.09e-12 |
| si16 | 601 | 1e-03 | 5/12/1 | 495/823/103 | 1421 | 71/16 | 1373.97 | 0.0503869 | 1.90e-10 |
| si16 | 601 | 3e-04 | 8/12/2 | 616/916/179 | 1711 | 87/20 | 1628.44 | 0.0503641 | 9.31e-12 |
| si16 | 601 | 1e-04 | 10/10/2 | 710/916/197 | 1823 | 88/20 | 1654.90 | 0.0503258 | 6.93e-11 |
| si16 | 601 | 3e-05 | 12/7/2 | 808/787/185 | 1780 | 80/19 | 1503.82 | 0.0503217 | 1.88e-10 |
| si16 | 602 | 1e-02 | 2/20/3 | 225/1250/281 | 1756 | 87/23 | 2059.61 | 0.0495157 | 9.30e-12 |
| si16 | 602 | 3e-03 | 3/18/3 | 293/1456/292 | 2041 | 121/22 | 2268.44 | 0.0500627 | 5.27e-11 |
| si16 | 602 | 1e-03 | 5/12/1 | 486/845/102 | 1433 | 86/16 | 1391.72 | 0.0503907 | 1.42e-10 |
| si16 | 602 | 3e-04 | 9/11/3 | 635/944/290 | 1869 | 100/21 | 1781.54 | 0.0503537 | 5.35e-11 |
| si16 | 602 | 1e-04 | 10/10/2 | 688/909/191 | 1788 | 90/20 | 1633.95 | 0.0503286 | 1.09e-10 |
| si16 | 602 | 3e-05 | 11/9/1 | 738/838/92 | 1668 | 82/19 | 1449.27 | 0.0503277 | 2.42e-10 |
| si16 | 603 | 1e-02 | 2/13/3 | 226/1092/296 | 1614 | 87/16 | 1776.98 | 0.0496812 | 7.88e-11 |
| si16 | 603 | 3e-03 | 3/15/3 | 290/1155/276 | 1721 | 86/19 | 1890.36 | 0.0502257 | 1.21e-11 |
| si16 | 603 | 1e-03 | 5/11/2 | 485/962/177 | 1624 | 78/16 | 1583.35 | 0.05034 | 6.32e-12 |
| si16 | 603 | 3e-04 | 6/10/2 | 541/932/195 | 1668 | 75/16 | 1578.84 | 0.0503107 | 3.59e-11 |
| si16 | 603 | 1e-04 | 6/10/2 | 541/932/195 | 1668 | 75/16 | 1578.84 | 0.0503107 | 3.59e-11 |
| si16 | 603 | 3e-05 | 8/9/3 | 688/995/274 | 1957 | 85/18 | 1788.48 | 0.0503198 | 5.35e-12 |
| si32 | 601 | 1e-02 | 2/14/2 | 435/2216/365 | 3016 | 88/16 | 3304.98 | 0.0954876 | 2.08e-11 |
| si32 | 601 | 3e-03 | 3/18/3 | 557/2662/586 | 3805 | 103/22 | 4213.70 | 0.09712 | 7.74e-11 |
| si32 | 601 | 1e-03 | 3/18/3 | 557/2662/586 | 3805 | 103/22 | 4213.70 | 0.09712 | 7.74e-11 |
| si32 | 601 | 3e-04 | 6/10/2 | 1031/1757/393 | 3181 | 85/16 | 3019.05 | 0.0976384 | 5.60e-11 |
| si32 | 601 | 1e-04 | 7/10/1 | 1146/1717/201 | 3064 | 72/16 | 2810.03 | 0.0976349 | 1.06e-10 |
| si32 | 601 | 3e-05 | 7/10/1 | 1146/1717/201 | 3064 | 72/16 | 2810.03 | 0.0976349 | 1.06e-10 |
| si32 | 602 | 1e-02 | 2/14/2 | 436/2173/375 | 2984 | 78/16 | 3272.38 | 0.0956369 | 3.15e-11 |
| si32 | 602 | 3e-03 | 3/17/2 | 560/2503/392 | 3455 | 96/20 | 3789.88 | 0.0970177 | 1.19e-10 |
| si32 | 602 | 1e-03 | 3/17/2 | 560/2503/392 | 3455 | 96/20 | 3789.88 | 0.0970177 | 1.19e-10 |
| si32 | 602 | 3e-04 | 6/11/2 | 1027/1762/391 | 3180 | 88/17 | 3056.48 | 0.0976439 | 7.26e-11 |
| si32 | 602 | 1e-04 | 7/10/1 | 1145/1735/197 | 3077 | 79/16 | 2823.64 | 0.0976358 | 1.35e-10 |
| si32 | 602 | 3e-05 | 7/10/1 | 1145/1735/197 | 3077 | 79/16 | 2823.64 | 0.0976358 | 1.35e-10 |
| si32 | 603 | 1e-02 | 2/14/2 | 434/2178/362 | 2974 | 89/16 | 3263.59 | 0.0955007 | 2.35e-11 |
| si32 | 603 | 3e-03 | 3/17/2 | 559/2405/390 | 3354 | 88/20 | 3689.48 | 0.0971074 | 4.69e-11 |
| si32 | 603 | 1e-03 | 3/17/2 | 559/2405/390 | 3354 | 88/20 | 3689.48 | 0.0971074 | 4.69e-11 |
| si32 | 603 | 3e-04 | 6/11/2 | 1030/1809/378 | 3217 | 87/17 | 3091.66 | 0.0976366 | 6.39e-11 |
| si32 | 603 | 1e-04 | 6/11/2 | 1030/1809/378 | 3217 | 87/17 | 3091.66 | 0.0976366 | 6.39e-11 |
| si32 | 603 | 3e-05 | 7/10/1 | 1171/1751/189 | 3111 | 78/16 | 2841.80 | 0.0976411 | 7.70e-11 |
| si8 | 601 | 1e-02 | 3/10/2 | 155/396/99 | 650 | 58/13 | 686.95 | 0.0258082 | 1.01e-10 |
| si8 | 601 | 3e-03 | 3/10/2 | 155/396/99 | 650 | 58/13 | 686.95 | 0.0258082 | 1.01e-10 |
| si8 | 601 | 1e-03 | 4/9/2 | 211/396/97 | 704 | 68/13 | 695.72 | 0.0257189 | 7.14e-11 |
| si8 | 601 | 3e-04 | 5/11/2 | 238/430/101 | 769 | 66/16 | 770.78 | 0.0256789 | 8.69e-11 |
| si8 | 601 | 1e-04 | 5/11/2 | 238/430/101 | 769 | 66/16 | 770.78 | 0.0256789 | 8.69e-11 |
| si8 | 601 | 3e-05 | 8/8/2 | 313/410/98 | 821 | 78/16 | 748.95 | 0.0256624 | 9.02e-11 |
| si8 | 602 | 1e-02 | 3/11/2 | 155/436/99 | 690 | 72/14 | 738.95 | 0.0258083 | 3.61e-11 |
| si8 | 602 | 3e-03 | 3/11/2 | 155/436/99 | 690 | 72/14 | 738.95 | 0.0258083 | 3.61e-11 |
| si8 | 602 | 1e-03 | 4/10/2 | 212/400/98 | 710 | 69/14 | 713.05 | 0.0257229 | 6.22e-11 |
| si8 | 602 | 3e-04 | 5/9/3 | 241/407/157 | 805 | 79/15 | 792.79 | 0.025673 | 1.39e-10 |
| si8 | 602 | 1e-04 | 6/9/2 | 269/398/96 | 763 | 69/15 | 724.18 | 0.0256612 | 1.33e-10 |
| si8 | 602 | 3e-05 | 9/7/2 | 353/379/102 | 834 | 69/16 | 727.37 | 0.0256636 | 1.27e-10 |
| si8 | 603 | 1e-02 | 3/10/2 | 155/394/100 | 649 | 61/13 | 685.95 | 0.0258108 | 3.55e-10 |
| si8 | 603 | 3e-03 | 3/10/2 | 155/394/100 | 649 | 61/13 | 685.95 | 0.0258108 | 3.55e-10 |
| si8 | 603 | 1e-03 | 4/9/2 | 210/391/95 | 696 | 77/13 | 688.38 | 0.0257211 | 5.11e-11 |
| si8 | 603 | 3e-04 | 5/9/2 | 238/423/100 | 761 | 79/14 | 738.78 | 0.0256782 | 1.05e-10 |
| si8 | 603 | 1e-04 | 5/9/2 | 238/423/100 | 761 | 79/14 | 738.78 | 0.0256782 | 1.05e-10 |
| si8 | 603 | 3e-05 | 8/8/2 | 341/399/96 | 836 | 76/16 | 745.33 | 0.0256644 | 1.46e-11 |
| syn-skew2 | 601 | 1e-02 | 3/9/1 | 22/61/4 | 87 | 27/11 | 102.37 | 0.0134546 | 7.73e-10 |
| syn-skew2 | 601 | 3e-03 | 3/9/1 | 22/61/4 | 87 | 27/11 | 102.37 | 0.0134546 | 7.73e-10 |
| syn-skew2 | 601 | 1e-03 | 6/8/1 | 43/61/5 | 109 | 36/13 | 111.48 | 0.0133417 | 6.00e-11 |
| syn-skew2 | 601 | 3e-04 | 7/7/1 | 47/62/6 | 115 | 36/13 | 113.07 | 0.0133448 | 2.05e-10 |
| syn-skew2 | 601 | 1e-04 | 10/5/1 | 62/60/5 | 127 | 44/14 | 112.96 | 0.013343 | 8.03e-10 |
| syn-skew2 | 601 | 3e-05 | 12/5/1 | 70/65/5 | 140 | 48/16 | 123.15 | 0.0133423 | 2.86e-10 |
| syn-skew2 | 602 | 1e-02 | 3/9/1 | 22/61/4 | 87 | 27/11 | 102.37 | 0.0134577 | 6.13e-10 |
| syn-skew2 | 602 | 3e-03 | 3/9/1 | 22/61/4 | 87 | 27/11 | 102.37 | 0.0134577 | 6.13e-10 |
| syn-skew2 | 602 | 1e-03 | 6/8/1 | 43/60/6 | 109 | 36/13 | 111.48 | 0.0133424 | 1.02e-10 |
| syn-skew2 | 602 | 3e-04 | 7/7/1 | 47/62/6 | 115 | 36/13 | 113.07 | 0.0133464 | 2.15e-10 |
| syn-skew2 | 602 | 1e-04 | 10/5/1 | 62/61/5 | 128 | 44/14 | 113.96 | 0.0133447 | 8.46e-10 |
| syn-skew2 | 602 | 3e-05 | 12/5/1 | 70/63/8 | 141 | 48/16 | 124.15 | 0.013344 | 1.13e-09 |
| syn-skew2 | 603 | 1e-02 | 3/7/1 | 22/52/4 | 78 | 25/9 | 87.37 | 0.0134586 | 7.42e-10 |
| syn-skew2 | 603 | 3e-03 | 3/7/1 | 22/52/4 | 78 | 25/9 | 87.37 | 0.0134586 | 7.42e-10 |
| syn-skew2 | 603 | 1e-03 | 6/8/1 | 43/61/7 | 111 | 37/13 | 113.48 | 0.0133422 | 1.07e-10 |
| syn-skew2 | 603 | 3e-04 | 7/7/1 | 47/62/8 | 117 | 38/13 | 115.07 | 0.0133462 | 2.52e-10 |
| syn-skew2 | 603 | 1e-04 | 10/5/1 | 62/62/5 | 129 | 45/14 | 114.96 | 0.0133445 | 7.22e-10 |
| syn-skew2 | 603 | 3e-05 | 12/5/1 | 70/63/8 | 141 | 48/16 | 124.15 | 0.0133438 | 1.15e-09 |
| syn-slab | 601 | 1e-02 | 4/10/1 | 31/66/4 | 101 | 34/13 | 111.72 | 0.034691 | 1.30e-10 |
| syn-slab | 601 | 3e-03 | 5/9/1 | 35/76/7 | 118 | 41/13 | 123.77 | 0.0332006 | 1.04e-10 |
| syn-slab | 601 | 1e-03 | 5/9/1 | 35/76/7 | 118 | 41/13 | 123.77 | 0.0332006 | 1.04e-10 |
| syn-slab | 601 | 3e-04 | 6/8/1 | 43/78/5 | 126 | 42/13 | 123.99 | 0.033102 | 1.30e-09 |
| syn-slab | 601 | 1e-04 | 8/6/1 | 51/78/7 | 136 | 47/13 | 124.09 | 0.0330795 | 3.80e-10 |
| syn-slab | 601 | 3e-05 | 9/6/1 | 55/83/6 | 144 | 50/14 | 130.14 | 0.0330794 | 1.39e-10 |
| syn-slab | 602 | 1e-02 | 4/10/1 | 31/65/4 | 100 | 33/13 | 110.72 | 0.034744 | 2.76e-10 |
| syn-slab | 602 | 3e-03 | 5/7/1 | 35/63/6 | 104 | 34/11 | 103.77 | 0.0332503 | 1.23e-09 |
| syn-slab | 602 | 1e-03 | 5/7/1 | 35/63/6 | 104 | 34/11 | 103.77 | 0.0332503 | 1.23e-09 |
| syn-slab | 602 | 3e-04 | 6/9/1 | 43/80/4 | 127 | 41/14 | 127.99 | 0.0331494 | 1.81e-10 |
| syn-slab | 602 | 1e-04 | 8/6/1 | 51/78/7 | 136 | 47/13 | 124.09 | 0.0331272 | 3.96e-10 |
| syn-slab | 602 | 3e-05 | 9/6/1 | 55/83/6 | 144 | 50/14 | 130.14 | 0.0331272 | 1.37e-10 |
| syn-slab | 603 | 1e-02 | 4/10/1 | 31/64/5 | 100 | 34/13 | 110.72 | 0.0347632 | 2.26e-10 |
| syn-slab | 603 | 3e-03 | 5/7/1 | 35/63/5 | 103 | 34/11 | 102.77 | 0.0332891 | 8.76e-10 |
| syn-slab | 603 | 1e-03 | 5/7/1 | 35/63/5 | 103 | 34/11 | 102.77 | 0.0332891 | 8.76e-10 |
| syn-slab | 603 | 3e-04 | 6/7/1 | 45/70/8 | 123 | 43/12 | 116.58 | 0.0331757 | 5.77e-10 |
| syn-slab | 603 | 1e-04 | 8/6/1 | 53/76/7 | 136 | 46/13 | 122.68 | 0.0331433 | 4.85e-10 |
| syn-slab | 603 | 3e-05 | 9/6/1 | 58/77/6 | 141 | 48/14 | 125.02 | 0.0331437 | 3.91e-10 |

### Held-out policy comparison

Training seed 601 selects q=0.1, fixed stop=0.01; evaluation seeds 602/603. Probe columns are included in the policy cost.

| System | Baseline W | Policy W | Fixed W | Paired policy/baseline | Win |
| --- | ---: | ---: | ---: | ---: | --- |
| mgo8 | 1318.17 | 1195.13 | 900.80 | 0.90696 | True |
| si-skew | 283.93 | 283.93 | 207.86 | 1.00000 | False |
| si16 | 1606.39 | 1728.19 | 1918.29 | 1.07488 | False |
| si32 | 2957.65 | 3170.07 | 3267.98 | 1.07425 | False |
| si8 | 731.48 | 716.72 | 712.45 | 0.98009 | True |
| syn-skew2 | 114.46 | 114.48 | 94.87 | 1.00014 | False |
| syn-slab | 123.39 | 105.27 | 110.72 | 0.85317 | True |
