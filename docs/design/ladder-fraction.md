# F-LADDER-FRACTION: fraction and fine-stage recovery

**Scientific result: CLOSED, negative for the tested rule families. No default
change and no adoption recommendation.** Mac pure-Rust release, one core.
This is not dev08 acceptance or a QE run.

The seven-system, three-seed matrix completed **126/126** arms with unchanged
fixed-point verification. The grid-plateau rule wins grid-weighted total H on
**5/7** systems but leaves increased fine work/maps on **5/7**. It wins raw
(unweighted) total H on only **2/7**. The mode-normalized stop wins H-eq on
4/7 and leaves five regression systems. Combining both still has six
regression systems. The watch is numerically inert and cannot remove work.
Thus no rule meets the requested acceptance condition under either H measure.

Worst absolute difference from the same-geometry strict reference:
**4.310998447110137e-9 Ha**, syn-slab seed 603, grid (also grid+modes).
The remaining margin to 5e-9 Ha is **6.89001552889863e-10 Ha (13.78%)**.
That observed margin is narrow; it is not a transferable error guarantee.
The default 0.5 arm's worst difference is 8.463689749049763e-10 Ha.

The focused si32 seed-601 loss is explained by **both** the coarse stop
crossing one map earlier and the resulting roughly 27x looser auto tolerance
seed. A factorial state/seed exchange confirms both effects. Si64 has no
such stop/seed discontinuity. The grid and modes rules keep its large H-eq
win (about 16.7–16.8%, one geometry), but do not remove losses generally.

**Required local delivery gate: red**, exit 1: build/fmt/clippy pass;
468 passed, 3 failed, 134 ignored across 51 binaries. All three failures
are the existing CLI cache tests rejecting in-repository TMPDIR. No gate
or scientific threshold was relaxed; detailed validation is below.

## Protocol fixed before candidate implementation

Reuse `scripts/ladder_handoff.py` and the seven F-LADDER-HANDOFF systems,
seeds 601–603, fractional displacement +/-0.001; add si64 seed 601 for the
mechanism screen. Sweep fractions 0.5, 0.45, 0.4, 0.35 with the existing
1e-4 coarse residual stop and auto first tolerance. Every geometry has its
own accurately converged strict fine reference, unchanged verification on
all standard arms, and full coarse/fine/verification H accounting.

Per-map files retain the standard error, requested eigen tolerance, full-H
columns, Davidson iterations, density residual, mixer resets and the complete
Pulay record. Coulomb guard `trace` measures the predicted Coulomb norm
without refitting. Diagnostics execute after the timed calculation in the
ignored release test; they do not enter production.

Handoff density errors use sqrt(volume * mean(delta_density^2))/electrons,
Bohr^-3/2, the engine convention. Low G means the coarse **wavefunction**
sphere. A separate column uses the larger coarse **density** sphere (4Ec).
Hartree and XC response errors are RMS potential differences (Ha), not
energy-error bounds. Occupied-subspace principal angles are acos of the
singular values of X_handoff^H X_strict (radians, Gamma-real overlaps).
The handoff orbital energy is the full-grid variational functional of the
injected occupied orbitals, using their own density. The mixed handoff input
density does not define a kinetic energy by itself; its error is reported
separately, together with the first fine map's energy error.

Prior evidence: F-LADDER-HANDOFF rejected a trained high-G residual probe
policy (3/7 wins). Its indicator varied only 0.58–5.08% within a system and
was not transferable to remaining fine work. This task does not repeat that
policy. SCF-TRAJECTORY identifies a first-Pulay Coulomb bounce that persists
under tight inner solves; this task measures its dependence on coarse basis.

Candidate families to evaluate independently, without fitting to materials:

1. Fill the FFT plateau: starting from the requested lower fraction, choose
   the largest wavefunction sphere that fits its already selected coarse
   density grid. The discrete grid and Miller vectors define the choice.
   The extra basis may cost algebra; charge every H column and measure wall
   on dev08 before adoption.
2. A mode-normalized coarse residual stop: use Parseval's relation to keep
   the RMS error budget per represented density Fourier mode equal to the
   adopted half-cutoff level. This is a testable heuristic, not an estimate
   of the missing wavefunctions or a certified discretization bound.
3. A diagnostic safeguard: record a failure of the initial fine-map
   contraction forecast, without restarting or discarding any work. It is
   an observation, not a claim to recover work already spent.

Acceptance: better grid-weighted total full H on at least five of seven
systems and no increased fine work on any; also report unweighted total H,
all fine and verification maps, and the worst strict-reference energy error
against 5e-9 Ha. A failed family is a negative result, never a default change.

## First mechanism result, before runtime-rule implementation

Si32 seed 601 reproduces 10 / 10 / 13 / 11 fine maps at fractions
0.5 / 0.45 / 0.4 / 0.35. At 0.4 the coarse stage ends after six maps with
residual 9.9494722265e-5, barely below 1e-4; 0.5 and 0.45 take seven maps
and end at 2.6783792145e-5 and 2.5777905581e-5. The 0.4 handoff density
error is 6.2312274823e-4 versus 3.2041795911e-4 at 0.5, and maximum
occupied angle is 0.042234 versus 0.019960 rad. All three have an early
Pulay bounce: fine map 3 standard errors are 0.00395034 / 0.00402208 /
0.00404397 Ha. Thus the presence of that bounce alone does not identify
the three-extra-map arm. Residual-stop and first-tolerance ablations are
needed before attributing the longer tail to one cause.

## Implemented rule contracts

All three switches are run-owned, default off, accepted in `[performance]`
and overridden by the corresponding `BAANDI_` environment variable. Existing
`cutoff_ladder` input/environment precedence is unchanged. No numerical hot
loop reads a new switch or writes a diagnostic file. `session.rs` is unchanged.

### `ladder_fraction_rule = "grid"`

Let D be the smooth FFT dimensions selected for the requested fraction f.
The density sphere has kinetic cutoff 4 f Ef. Enumerate the final density
sphere and find its lowest-energy vector g for which any
`2*abs(miller_i)+1 > D_i`. The coarse wavefunction cutoff must remain below
`|g|²/8`; divide out the basis inclusion tolerance and round toward the
interior. Rebuild and check that all three dimensions are still exactly D.
Otherwise keep the requested fraction. A second application is idempotent.
The construction works for skewed cells by using actual reciprocal vectors,
not a cubic approximation. The requested fraction (0.4 in this experiment)
is an input, not a hard-coded branch or newly adopted default.

This fills an FFT plateau; it is not a free operation. The additional plane
waves can increase algebra and projector work, and enumeration adds setup.
Those costs are not represented by a grid-volume-only model. Production wall
measurements must include setup. The effective fraction and dimensions are
recorded in the switch path note and the ordinary ladder report.

### `ladder_stop_rule = "modes"`

For M nonzero modes Parseval gives
`R² = volume/electrons² * sum_G |delta_n_G|²`. Equal RMS Fourier coefficient
budget therefore means equal `R/sqrt(M)`. With Mhalf counted at the existing
half-cutoff density sphere, use

    Rstop = Rrequested * min(1, sqrt(Mcoarse / Mhalf)).

The zero mode is excluded because charge conservation fixes it. There is no
new fitted numerical constant; the existing half-cutoff policy defines the
comparison. This rule tightens a smaller-basis stop and cannot loosen the
requested stop. It is **not** a reliable discretization-error estimator:
mode amplitudes need not be equal and the excluded wavefunctions also alter
represented density modes. It is retained as an independently switchable
negative/positive experiment, not justified as an error bound. Its counts
and effective residual level are in the switch path note.

### `ladder_watch = "on"`

With E1 and E2 the first two fine standard errors, let q=E2/E1. When
0<q<1, forecast `1+ceil(log(threshold/E1)/log(q))` total maps. Flag an early
excess if map 3 or 4 exceeds `E1*q^(k-1)`. The path note, filed after the
fine session, contains this forecast, the first excess and the observed map
count. Invalid/noncontracting initial data get an unavailable note. No map
is discarded, repeated or changed. The nonlinear Pulay map need not obey this
forecast, and a flag is not a diagnosis of a failed calculation.

The test-only `first-previous` arm transfers the previous arm's coarse-next
tolerance to the current arm's handoff state. This isolates density/orbital
quality from the tolerance seed in a two-by-two mechanism ablation. It is
not a production switch. The `first-tight` arm uses the existing standard
Ritz floor only as an ablation and is not an optimization candidate.

## Causal isolation on si32, seed 601

These are fresh completed runs, not estimated work saved. `seed` means the
coarse-next tolerance supplied to the fine session; per-map files record the
actual tolerance after the ordinary first-map redo rule.

| Coarse maps / state | Seed source | First fine tolerance (Ha) | Fine maps | Fine + verification H |
| --- | --- | ---: | ---: | ---: |
| 6, stop 1e-4 | same 6-map state | 2.6120908199e-6 | 13 | 2233 |
| 6, stop 1e-4 | 7-map state | 8.1258231119e-8 | 11 | 2471 |
| 7, stop 3e-5 | 6-map state | 2.6120908199e-6 | 11 | 1966 |
| 7, stop 3e-5 | same 7-map state | 8.1258231119e-8 | 10 | 2017 |

Both factors matter. A more accurate coarse state at the old loose seed
removes two maps. The tighter seed at the original state also removes two
maps, but **increases** fine H by 10.66%. Together they remove three maps.
An extreme first seed at the existing 5e-14 Ha Ritz floor reaches ten maps
but costs 5206 fine+verification H. It is not an optimization. A cold-cap
start invokes the existing redo, ends up at 1.39e-6 Ha and takes eleven maps;
using the last coarse tolerance takes twelve. No final tolerance or
verification threshold was changed in these ablations.

The residual stop has a discrete crossing: 0.4 hands off one map earlier,
while the auto seed reads **the last coarse standard error**, not the
representation error on the fine basis. The resulting seed loosens from
9.59e-8 (0.5) or 8.79e-8 Ha (0.45) to 2.61e-6 Ha (0.4), roughly 27–30x.
The fine map then starts from a worse density/orbital state with a looser
inner solve. The matrix tests how well general rules prevent that crossing.

For si64 seed 601 the coarse count stays six for all four fractions. The
first fine tolerance is 4.30/4.16/3.86/3.83e-6 Ha, so there is no analogous
27x discontinuity. At 0.4 the fine count falls from thirteen to twelve,
despite the occupied angle and initial density error becoming larger.
This is why a rule based solely on discretization error cannot predict the
remaining nonlinear SCF work.

All focused si32/si64 arms start with empty fine Pulay history. Their map-2
fit extrapolates and map 3 rises in Coulomb error, including the ten-map
si32 controls. The 0.4 arm does not have a uniquely large map-3 error:
0.003950/0.004022/0.004044 Ha at 0.5/0.45/0.4. The slower recovery is visible
at map 4 (7.56e-6/7.61e-6/1.47e-4 Ha). Thus a bounce flag alone is not a
classifier for this loss. We did not change the handoff history reset or
repeat the rejected secant-transfer policy; an interaction with its shared
history-building transient is not excluded by this study.

The representation error is real, but is not a sufficient explanation of the
extra maps. At fraction 0.4, si32/si64 have almost the same Hartree potential
error (0.001528/0.001529 Ha RMS), XC error (0.005221/0.005236 Ha RMS) and
maximum angle (0.04223/0.04258 rad). Both roughly double their handoff
density error relative to 0.5, yet their map-count changes have opposite
signs. In si32/0.4 the density error outside the coarse **density** sphere is
1.4957e-5 versus total 6.2312e-4 Bohr^-3/2: only about 0.058% of the squared
norm. Missing orbital response mostly changes density modes the coarse
*density grid* can already represent. Enlarging/interpolating that grid alone
cannot repair it. The high-PW density-error component has a substantial XC
response (0.004047 Ha RMS), so it must not be treated as energetically inert.

The safeguard flags 11 of 21 lower-fraction runs. Si32 seed 601 is flagged,
but its initial forecast is thirteen maps and it actually takes thirteen;
si16 seed 601 forecasts 27 and actually takes eleven. Only syn-slab seeds
601/602 exceed the final forecast (seven versus six maps). This is a useful
record of nonlinear departure, not a validated predictor of the si32 loss.
All 21 watch/off pairs have identical energies, work, map counts, handoff
arrays and scalar per-map histories.

## Supplemental si64 rule screen (one geometry, seed 601)

This is outside the seven-system acceptance matrix. It checks whether the
large-system benefit survives, not a three-seed adoption result. No timing
ratio is inferred from these work counts.

| Arm | Effective f | C/F/V maps | Hc/Hf/Hv | Raw total H | H-eq | H-eq / 0.5 | abs dE (Ha) |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: |
| 0.5:1e-4:default | 0.50000000 | 6/13/2 | 2023/3871/750 | 6644 | 5474.453 | 1.000000 | 5.31202e-11 |
| 0.4:1e-4:default | 0.40000000 | 6/12/2 | 2013/3334/727 | 6074 | 4552.455 | 0.831582 | 3.63229e-11 |
| 0.4:1e-4:grid | 0.41321567 | 6/12/2 | 2014/3342/727 | 6083 | 4560.699 | 0.833088 | 3.76588e-11 |
| 0.4:1e-4:modes | 0.40000000 | 6/12/2 | 2013/3334/727 | 6074 | 4552.455 | 0.831582 | 3.63229e-11 |
| 0.4:1e-4:grid-modes | 0.41321567 | 6/12/2 | 2014/3342/727 | 6083 | 4560.699 | 0.833088 | 3.76588e-11 |

## Candidate matrix (Mac one core, three geometries per system)

H includes coarse atomic start, fine redo and every verification. H-eq weights coarse H by the coarse/fine grid volume. Ratios are medians of paired ratios. Fine H includes verification; `increases` counts seeds with more fine H / more fine maps. W also charges delta-potential FFT columns. No wall claim is made.

| System | Arm | Effective f | C/F/V maps | Total H | H-eq | H/base | H-eq/base | W/base | Fine H/base | Increases H/maps | Worst abs dE (Ha) |
| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| mgo8 | 0.5 | 0.5000000 | 8/8/1 | 1531 | 1091.89 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0/0 | 4.1496e-11 |
| mgo8 | 0.4 | 0.4000000 | 8/8/1 | 1502 | 1025.33 | 0.9731 | 0.9390 | 0.9387 | 0.9861 | 1/0 | 6.3812e-10 |
| mgo8 | grid | 0.4365844 | 9/8/1 | 1636 | 1126.40 | 1.0810 | 1.0462 | 1.0335 | 1.1041 | 3/0 | 5.5223e-10 |
| mgo8 | modes | 0.4000000 | 8/8/1 | 1520 | 1051.95 | 0.9731 | 0.9416 | 0.9387 | 1.0045 | 2/0 | 6.4642e-10 |
| mgo8 | watch | 0.4000000 | 8/8/1 | 1502 | 1025.33 | 0.9731 | 0.9390 | 0.9387 | 0.9861 | 1/0 | 6.3812e-10 |
| mgo8 | grid-modes | 0.4365844 | 9/8/1 | 1636 | 1126.40 | 1.0810 | 1.0462 | 1.0335 | 1.1041 | 3/0 | 5.5223e-10 |
| si-skew | 0.5 | 0.5000000 | 11/8/3 | 278 | 202.16 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0/0 | 2.6628e-12 |
| si-skew | 0.4 | 0.4000000 | 13/7/2 | 278 | 176.51 | 1.0000 | 0.8731 | 0.8546 | 0.8981 | 0/0 | 1.5852e-11 |
| si-skew | grid | 0.4007499 | 13/7/2 | 278 | 176.51 | 1.0000 | 0.8731 | 0.8546 | 0.8981 | 0/0 | 1.5863e-11 |
| si-skew | modes | 0.4000000 | 13/7/2 | 278 | 176.51 | 1.0000 | 0.8731 | 0.8546 | 0.8981 | 0/0 | 1.5852e-11 |
| si-skew | watch | 0.4000000 | 13/7/2 | 278 | 176.51 | 1.0000 | 0.8731 | 0.8546 | 0.8981 | 0/0 | 1.5852e-11 |
| si-skew | grid-modes | 0.4007499 | 13/7/2 | 278 | 176.51 | 1.0000 | 0.8731 | 0.8546 | 0.8981 | 0/0 | 1.5863e-11 |
| si16 | 0.5 | 0.5000000 | 10/10/2 | 1788 | 1348.84 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0/0 | 1.0898e-10 |
| si16 | 0.4 | 0.4000000 | 10/11/2 | 1858 | 1360.39 | 1.0391 | 1.0086 | 1.0101 | 1.0627 | 2/2 | 8.7880e-11 |
| si16 | grid | 0.4132157 | 10/10/2 | 1816 | 1322.44 | 1.0258 | 0.9998 | 0.9945 | 1.0408 | 3/1 | 3.9158e-11 |
| si16 | modes | 0.4000000 | 10/11/2 | 1858 | 1360.39 | 1.0391 | 1.0086 | 1.0101 | 1.0627 | 2/2 | 8.7880e-11 |
| si16 | watch | 0.4000000 | 10/11/2 | 1858 | 1360.39 | 1.0391 | 1.0086 | 1.0101 | 1.0627 | 2/2 | 8.7880e-11 |
| si16 | grid-modes | 0.4132157 | 10/10/2 | 1816 | 1322.44 | 1.0258 | 0.9998 | 0.9945 | 1.0408 | 3/1 | 3.9158e-11 |
| si32 | 0.5 | 0.5000000 | 7/10/1 | 3077 | 2379.27 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0/0 | 1.3490e-10 |
| si32 | 0.4 | 0.4000000 | 6/12/2 | 3238 | 2478.45 | 1.0318 | 0.9807 | 0.9738 | 1.0543 | 3/2 | 3.0015e-10 |
| si32 | grid | 0.4132157 | 7/10/1 | 3175 | 2332.61 | 1.0318 | 0.9804 | 0.9735 | 1.0538 | 2/1 | 2.9036e-10 |
| si32 | modes | 0.4000000 | 7/10/1 | 3175 | 2333.35 | 1.0304 | 0.9781 | 0.9716 | 1.0516 | 3/1 | 3.0086e-10 |
| si32 | watch | 0.4000000 | 6/12/2 | 3238 | 2478.45 | 1.0318 | 0.9807 | 0.9738 | 1.0543 | 3/2 | 3.0015e-10 |
| si32 | grid-modes | 0.4132157 | 7/10/1 | 3175 | 2332.61 | 1.0318 | 0.9804 | 0.9735 | 1.0538 | 2/1 | 2.9036e-10 |
| si8 | 0.5 | 0.5000000 | 5/9/2 | 763 | 602.71 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0/0 | 1.3308e-10 |
| si8 | 0.4 | 0.4000000 | 5/9/2 | 738 | 572.22 | 0.9597 | 0.9388 | 0.9335 | 0.9492 | 0/0 | 3.3551e-10 |
| si8 | grid | 0.4499286 | 5/9/2 | 721 | 554.93 | 0.9376 | 0.9087 | 0.8941 | 0.9134 | 0/0 | 2.3528e-10 |
| si8 | modes | 0.4000000 | 6/9/2 | 796 | 597.93 | 1.0351 | 0.9921 | 0.9686 | 0.9790 | 0/0 | 5.2024e-10 |
| si8 | watch | 0.4000000 | 5/9/2 | 738 | 572.22 | 0.9597 | 0.9388 | 0.9335 | 0.9492 | 0/0 | 3.3551e-10 |
| si8 | grid-modes | 0.4499286 | 6/8/2 | 750 | 561.41 | 0.9753 | 0.9193 | 0.8915 | 0.9077 | 0/1 | 2.3528e-10 |
| syn-skew2 | 0.5 | 0.5000000 | 10/5/1 | 128 | 88.96 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0/0 | 8.4637e-10 |
| syn-skew2 | 0.4 | 0.4000000 | 15/6/1 | 161 | 101.74 | 1.2677 | 1.1309 | 1.1020 | 1.2836 | 3/3 | 2.6135e-10 |
| syn-skew2 | grid | 0.4513493 | 9/6/1 | 122 | 80.19 | 0.9531 | 0.8914 | 0.9025 | 1.0000 | 1/3 | 1.2515e-09 |
| syn-skew2 | modes | 0.4000000 | 15/6/1 | 161 | 101.74 | 1.2677 | 1.1309 | 1.1020 | 1.2836 | 3/3 | 2.6135e-10 |
| syn-skew2 | watch | 0.4000000 | 15/6/1 | 161 | 101.74 | 1.2677 | 1.1309 | 1.1020 | 1.2836 | 3/3 | 2.6135e-10 |
| syn-skew2 | grid-modes | 0.4513493 | 9/6/1 | 122 | 80.19 | 0.9531 | 0.8914 | 0.9025 | 1.0000 | 1/3 | 1.2515e-09 |
| syn-slab | 0.5 | 0.5000000 | 8/6/1 | 136 | 99.94 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0/0 | 4.8523e-10 |
| syn-slab | 0.4 | 0.4000000 | 9/7/1 | 159 | 112.30 | 1.1691 | 1.1236 | 1.1195 | 1.1529 | 3/2 | 4.3110e-09 |
| syn-slab | grid | 0.4441322 | 9/7/1 | 159 | 112.30 | 1.1691 | 1.1236 | 1.1195 | 1.1529 | 3/2 | 4.3110e-09 |
| syn-slab | modes | 0.4000000 | 9/7/1 | 159 | 112.30 | 1.1691 | 1.1236 | 1.1195 | 1.1529 | 3/2 | 4.3110e-09 |
| syn-slab | watch | 0.4000000 | 9/7/1 | 159 | 112.30 | 1.1691 | 1.1236 | 1.1195 | 1.1529 | 3/2 | 4.3110e-09 |
| syn-slab | grid-modes | 0.4441322 | 9/7/1 | 159 | 112.30 | 1.1691 | 1.1236 | 1.1195 | 1.1529 | 3/2 | 4.3110e-09 |

### Per-geometry accounting

| System | Seed | Arm | C/F/V maps | Hc/Hf/Hv | Total H | H-eq | abs dE (Ha) |
| --- | ---: | --- | --- | --- | ---: | ---: | ---: |
| si8 | 601 | 0.5:1e-4:default | 5/11/2 | 238/430/101 | 769 | 610.706 | 8.688517e-11 |
| si8 | 601 | 0.4:1e-4:default | 5/10/2 | 234/412/92 | 738 | 573.333 | 4.999379e-11 |
| si8 | 601 | 0.4:1e-4:grid | 5/9/2 | 236/389/96 | 721 | 554.926 | 2.276188e-10 |
| si8 | 601 | 0.4:1e-4:modes | 6/9/2 | 266/430/100 | 796 | 608.815 | 7.988277e-11 |
| si8 | 601 | 0.4:1e-4:watch | 5/10/2 | 234/412/92 | 738 | 573.333 | 4.999379e-11 |
| si8 | 601 | 0.4:1e-4:grid-modes | 6/8/2 | 268/383/99 | 750 | 561.407 | 7.566214e-11 |
| si8 | 602 | 0.5:1e-4:default | 6/9/2 | 269/398/96 | 763 | 584.088 | 1.330775e-10 |
| si8 | 602 | 0.4:1e-4:default | 5/9/1 | 237/391/55 | 683 | 516.222 | 3.355112e-10 |
| si8 | 602 | 0.4:1e-4:grid | 6/7/1 | 290/315/45 | 650 | 445.926 | 2.352820e-10 |
| si8 | 602 | 0.4:1e-4:modes | 6/7/1 | 288/333/47 | 668 | 465.333 | 5.202381e-10 |
| si8 | 602 | 0.4:1e-4:watch | 5/9/1 | 237/391/55 | 683 | 516.222 | 3.355112e-10 |
| si8 | 602 | 0.4:1e-4:grid-modes | 6/7/1 | 290/315/45 | 650 | 445.926 | 2.352820e-10 |
| si8 | 603 | 0.5:1e-4:default | 5/9/2 | 238/423/100 | 761 | 602.706 | 1.048903e-10 |
| si8 | 603 | 0.4:1e-4:default | 5/9/2 | 237/408/94 | 739 | 572.222 | 1.727507e-10 |
| si8 | 603 | 0.4:1e-4:grid | 5/9/2 | 237/413/95 | 745 | 578.222 | 5.675815e-11 |
| si8 | 603 | 0.4:1e-4:modes | 6/9/2 | 290/415/97 | 802 | 597.926 | 4.880718e-11 |
| si8 | 603 | 0.4:1e-4:watch | 5/9/2 | 237/408/94 | 739 | 572.222 | 1.727507e-10 |
| si8 | 603 | 0.4:1e-4:grid-modes | 6/10/2 | 293/421/101 | 815 | 608.815 | 1.100346e-10 |
| si16 | 601 | 0.5:1e-4:default | 10/10/2 | 710/916/197 | 1823 | 1369.800 | 6.934897e-11 |
| si16 | 601 | 0.4:1e-4:default | 10/11/2 | 701/1003/194 | 1898 | 1391.722 | 1.705303e-11 |
| si16 | 601 | 0.4:1e-4:grid | 10/11/2 | 703/999/182 | 1884 | 1376.278 | 1.153211e-11 |
| si16 | 601 | 0.4:1e-4:modes | 10/11/2 | 701/1003/194 | 1898 | 1391.722 | 1.705303e-11 |
| si16 | 601 | 0.4:1e-4:watch | 10/11/2 | 701/1003/194 | 1898 | 1391.722 | 1.705303e-11 |
| si16 | 601 | 0.4:1e-4:grid-modes | 10/11/2 | 703/999/182 | 1884 | 1376.278 | 1.153211e-11 |
| si16 | 602 | 0.5:1e-4:default | 10/10/2 | 688/909/191 | 1788 | 1348.843 | 1.089759e-10 |
| si16 | 602 | 0.4:1e-4:default | 10/11/2 | 689/974/195 | 1858 | 1360.389 | 4.178702e-11 |
| si16 | 602 | 0.4:1e-4:grid | 10/10/2 | 688/942/186 | 1816 | 1319.111 | 3.915801e-11 |
| si16 | 602 | 0.4:1e-4:modes | 10/11/2 | 689/974/195 | 1858 | 1360.389 | 4.178702e-11 |
| si16 | 602 | 0.4:1e-4:watch | 10/11/2 | 689/974/195 | 1858 | 1360.389 | 4.178702e-11 |
| si16 | 602 | 0.4:1e-4:grid-modes | 10/10/2 | 688/942/186 | 1816 | 1319.111 | 3.915801e-11 |
| si16 | 603 | 0.5:1e-4:default | 6/10/2 | 541/932/195 | 1668 | 1322.674 | 3.593925e-11 |
| si16 | 603 | 0.4:1e-4:default | 6/9/2 | 538/900/191 | 1629 | 1240.444 | 8.787993e-11 |
| si16 | 603 | 0.4:1e-4:grid | 6/10/2 | 538/979/194 | 1711 | 1322.444 | 2.884804e-11 |
| si16 | 603 | 0.4:1e-4:modes | 6/9/2 | 538/900/191 | 1629 | 1240.444 | 8.787993e-11 |
| si16 | 603 | 0.4:1e-4:watch | 6/9/2 | 538/900/191 | 1629 | 1240.444 | 8.787993e-11 |
| si16 | 603 | 0.4:1e-4:grid-modes | 6/10/2 | 538/979/194 | 1711 | 1322.444 | 2.884804e-11 |
| si32 | 601 | 0.5:1e-4:default | 7/10/1 | 1146/1717/201 | 3064 | 2365.656 | 1.057998e-10 |
| si32 | 601 | 0.4:1e-4:default | 6/13/2 | 1022/1872/361 | 3255 | 2499.146 | 2.067679e-11 |
| si32 | 601 | 0.4:1e-4:grid | 7/10/1 | 1143/1830/198 | 3171 | 2325.656 | 2.903562e-10 |
| si32 | 601 | 0.4:1e-4:modes | 7/10/1 | 1140/1823/194 | 3157 | 2313.875 | 3.008580e-10 |
| si32 | 601 | 0.4:1e-4:watch | 6/13/2 | 1022/1872/361 | 3255 | 2499.146 | 2.067679e-11 |
| si32 | 601 | 0.4:1e-4:grid-modes | 7/10/1 | 1143/1830/198 | 3171 | 2325.656 | 2.903562e-10 |
| si32 | 602 | 0.5:1e-4:default | 7/10/1 | 1145/1735/197 | 3077 | 2379.266 | 1.349036e-10 |
| si32 | 602 | 0.4:1e-4:default | 7/10/1 | 1138/1843/194 | 3175 | 2333.354 | 3.001475e-10 |
| si32 | 602 | 0.4:1e-4:grid | 7/10/1 | 1139/1841/195 | 3175 | 2332.615 | 2.453930e-10 |
| si32 | 602 | 0.4:1e-4:modes | 7/10/1 | 1138/1843/194 | 3175 | 2333.354 | 3.001475e-10 |
| si32 | 602 | 0.4:1e-4:watch | 7/10/1 | 1138/1843/194 | 3175 | 2333.354 | 3.001475e-10 |
| si32 | 602 | 0.4:1e-4:grid-modes | 7/10/1 | 1139/1841/195 | 3175 | 2332.615 | 2.453930e-10 |
| si32 | 603 | 0.5:1e-4:default | 6/11/2 | 1030/1809/378 | 3217 | 2589.344 | 6.393464e-11 |
| si32 | 603 | 0.4:1e-4:default | 6/12/2 | 1027/1825/386 | 3238 | 2478.448 | 5.407230e-11 |
| si32 | 603 | 0.4:1e-4:grid | 6/12/2 | 1030/1771/390 | 3191 | 2429.229 | 1.140990e-10 |
| si32 | 603 | 0.4:1e-4:modes | 6/12/2 | 1027/1825/386 | 3238 | 2478.448 | 5.407230e-11 |
| si32 | 603 | 0.4:1e-4:watch | 6/12/2 | 1027/1825/386 | 3238 | 2478.448 | 5.407230e-11 |
| si32 | 603 | 0.4:1e-4:grid-modes | 6/12/2 | 1030/1771/390 | 3191 | 2429.229 | 1.140990e-10 |
| si-skew | 601 | 0.5:1e-4:default | 11/8/3 | 121/121/36 | 278 | 202.163 | 2.638778e-12 |
| si-skew | 601 | 0.4:1e-4:default | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.353229e-11 |
| si-skew | 601 | 0.4:1e-4:grid | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.353051e-11 |
| si-skew | 601 | 0.4:1e-4:modes | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.353229e-11 |
| si-skew | 601 | 0.4:1e-4:watch | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.353229e-11 |
| si-skew | 601 | 0.4:1e-4:grid-modes | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.353051e-11 |
| si-skew | 602 | 0.5:1e-4:default | 11/8/3 | 122/122/36 | 280 | 203.536 | 2.662759e-12 |
| si-skew | 602 | 0.4:1e-4:default | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.286260e-11 |
| si-skew | 602 | 0.4:1e-4:grid | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.286438e-11 |
| si-skew | 602 | 0.4:1e-4:modes | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.286260e-11 |
| si-skew | 602 | 0.4:1e-4:watch | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.286260e-11 |
| si-skew | 602 | 0.4:1e-4:grid-modes | 13/7/2 | 137/120/21 | 278 | 176.510 | 1.286438e-11 |
| si-skew | 603 | 0.5:1e-4:default | 11/8/3 | 122/120/34 | 276 | 199.536 | 1.456613e-12 |
| si-skew | 603 | 0.4:1e-4:default | 13/7/2 | 136/127/23 | 286 | 185.251 | 1.585221e-11 |
| si-skew | 603 | 0.4:1e-4:grid | 13/7/2 | 136/127/23 | 286 | 185.251 | 1.586287e-11 |
| si-skew | 603 | 0.4:1e-4:modes | 13/7/2 | 136/127/23 | 286 | 185.251 | 1.585221e-11 |
| si-skew | 603 | 0.4:1e-4:watch | 13/7/2 | 136/127/23 | 286 | 185.251 | 1.585221e-11 |
| si-skew | 603 | 0.4:1e-4:grid-modes | 13/7/2 | 136/127/23 | 286 | 185.251 | 1.586287e-11 |
| mgo8 | 601 | 0.5:1e-4:default | 8/8/1 | 620/827/66 | 1513 | 1076.704 | 3.262812e-11 |
| mgo8 | 601 | 0.4:1e-4:default | 7/8/1 | 541/827/70 | 1438 | 1013.856 | 6.337473e-10 |
| mgo8 | 601 | 0.4:1e-4:grid | 9/8/1 | 650/918/68 | 1636 | 1126.400 | 5.522338e-10 |
| mgo8 | 601 | 0.4:1e-4:modes | 7/8/1 | 541/827/70 | 1438 | 1013.856 | 6.337473e-10 |
| mgo8 | 601 | 0.4:1e-4:watch | 7/8/1 | 541/827/70 | 1438 | 1013.856 | 6.337473e-10 |
| mgo8 | 601 | 0.4:1e-4:grid-modes | 9/8/1 | 650/918/68 | 1636 | 1126.400 | 5.522338e-10 |
| mgo8 | 602 | 0.5:1e-4:default | 8/8/2 | 626/788/148 | 1562 | 1121.481 | 2.609113e-11 |
| mgo8 | 602 | 0.4:1e-4:default | 8/8/1 | 597/846/77 | 1520 | 1051.952 | 5.666720e-10 |
| mgo8 | 602 | 0.4:1e-4:grid | 9/8/1 | 652/915/68 | 1635 | 1123.832 | 5.488801e-10 |
| mgo8 | 602 | 0.4:1e-4:modes | 8/8/1 | 597/846/77 | 1520 | 1051.952 | 5.666720e-10 |
| mgo8 | 602 | 0.4:1e-4:watch | 8/8/1 | 597/846/77 | 1520 | 1051.952 | 5.666720e-10 |
| mgo8 | 602 | 0.4:1e-4:grid-modes | 9/8/1 | 652/915/68 | 1635 | 1123.832 | 5.488801e-10 |
| mgo8 | 603 | 0.5:1e-4:default | 8/8/1 | 624/839/68 | 1531 | 1091.889 | 4.149570e-11 |
| mgo8 | 603 | 0.4:1e-4:default | 8/8/1 | 608/830/64 | 1502 | 1025.328 | 6.381242e-10 |
| mgo8 | 603 | 0.4:1e-4:grid | 9/8/1 | 653/934/68 | 1655 | 1143.048 | 5.521770e-10 |
| mgo8 | 603 | 0.4:1e-4:modes | 9/8/1 | 662/1002/74 | 1738 | 1218.992 | 6.464234e-10 |
| mgo8 | 603 | 0.4:1e-4:watch | 8/8/1 | 608/830/64 | 1502 | 1025.328 | 6.381242e-10 |
| mgo8 | 603 | 0.4:1e-4:grid-modes | 9/8/1 | 653/934/68 | 1655 | 1143.048 | 5.521770e-10 |
| syn-skew2 | 601 | 0.5:1e-4:default | 10/5/1 | 62/60/5 | 127 | 87.963 | 8.025247e-10 |
| syn-skew2 | 601 | 0.4:1e-4:default | 15/6/1 | 86/71/4 | 161 | 94.907 | 2.613521e-10 |
| syn-skew2 | 601 | 0.4:1e-4:grid | 8/6/1 | 51/65/5 | 121 | 81.806 | 7.066370e-11 |
| syn-skew2 | 601 | 0.4:1e-4:modes | 15/6/1 | 86/71/4 | 161 | 94.907 | 2.613521e-10 |
| syn-skew2 | 601 | 0.4:1e-4:watch | 15/6/1 | 86/71/4 | 161 | 94.907 | 2.613521e-10 |
| syn-skew2 | 601 | 0.4:1e-4:grid-modes | 8/6/1 | 51/65/5 | 121 | 81.806 | 7.066370e-11 |
| syn-skew2 | 602 | 0.5:1e-4:default | 10/5/1 | 62/61/5 | 128 | 88.963 | 8.463690e-10 |
| syn-skew2 | 602 | 0.4:1e-4:default | 15/6/1 | 86/81/7 | 174 | 107.907 | 1.576346e-10 |
| syn-skew2 | 602 | 0.4:1e-4:grid | 9/6/1 | 57/61/4 | 122 | 78.194 | 1.251490e-09 |
| syn-skew2 | 602 | 0.4:1e-4:modes | 15/6/1 | 86/81/7 | 174 | 107.907 | 1.576346e-10 |
| syn-skew2 | 602 | 0.4:1e-4:watch | 15/6/1 | 86/81/7 | 174 | 107.907 | 1.576346e-10 |
| syn-skew2 | 602 | 0.4:1e-4:grid-modes | 9/6/1 | 57/61/4 | 122 | 78.194 | 1.251490e-09 |
| syn-skew2 | 603 | 0.5:1e-4:default | 10/5/1 | 62/62/5 | 129 | 89.963 | 7.224288e-10 |
| syn-skew2 | 603 | 0.4:1e-4:default | 11/6/1 | 68/79/7 | 154 | 101.741 | 1.335179e-10 |
| syn-skew2 | 603 | 0.4:1e-4:grid | 9/6/1 | 57/63/4 | 124 | 80.194 | 3.848122e-11 |
| syn-skew2 | 603 | 0.4:1e-4:modes | 11/6/1 | 68/79/7 | 154 | 101.741 | 1.335179e-10 |
| syn-skew2 | 603 | 0.4:1e-4:watch | 11/6/1 | 68/79/7 | 154 | 101.741 | 1.335179e-10 |
| syn-skew2 | 603 | 0.4:1e-4:grid-modes | 9/6/1 | 57/63/4 | 124 | 80.194 | 3.848122e-11 |
| syn-slab | 601 | 0.5:1e-4:default | 8/6/1 | 51/78/7 | 136 | 99.941 | 3.804927e-10 |
| syn-slab | 601 | 0.4:1e-4:default | 9/7/1 | 61/93/5 | 159 | 112.297 | 1.913716e-10 |
| syn-slab | 601 | 0.4:1e-4:grid | 9/7/1 | 61/93/5 | 159 | 112.297 | 1.913700e-10 |
| syn-slab | 601 | 0.4:1e-4:modes | 9/7/1 | 61/93/5 | 159 | 112.297 | 1.913716e-10 |
| syn-slab | 601 | 0.4:1e-4:watch | 9/7/1 | 61/93/5 | 159 | 112.297 | 1.913716e-10 |
| syn-slab | 601 | 0.4:1e-4:grid-modes | 9/7/1 | 61/93/5 | 159 | 112.297 | 1.913700e-10 |
| syn-slab | 602 | 0.5:1e-4:default | 8/6/1 | 51/78/7 | 136 | 99.941 | 3.958160e-10 |
| syn-slab | 602 | 0.4:1e-4:default | 9/7/1 | 61/90/8 | 159 | 112.297 | 1.650890e-09 |
| syn-slab | 602 | 0.4:1e-4:grid | 9/7/1 | 61/90/8 | 159 | 112.297 | 1.650889e-09 |
| syn-slab | 602 | 0.4:1e-4:modes | 9/7/1 | 61/90/8 | 159 | 112.297 | 1.650890e-09 |
| syn-slab | 602 | 0.4:1e-4:watch | 9/7/1 | 61/90/8 | 159 | 112.297 | 1.650890e-09 |
| syn-slab | 602 | 0.4:1e-4:grid-modes | 9/7/1 | 61/90/8 | 159 | 112.297 | 1.650889e-09 |
| syn-slab | 603 | 0.5:1e-4:default | 8/6/1 | 53/76/7 | 136 | 98.527 | 4.852314e-10 |
| syn-slab | 603 | 0.4:1e-4:default | 9/6/1 | 61/84/9 | 154 | 107.297 | 4.310998e-09 |
| syn-slab | 603 | 0.4:1e-4:grid | 9/6/1 | 61/84/9 | 154 | 107.297 | 4.310998e-09 |
| syn-slab | 603 | 0.4:1e-4:modes | 9/6/1 | 61/84/9 | 154 | 107.297 | 4.310998e-09 |
| syn-slab | 603 | 0.4:1e-4:watch | 9/6/1 | 61/84/9 | 154 | 107.297 | 4.310998e-09 |
| syn-slab | 603 | 0.4:1e-4:grid-modes | 9/6/1 | 61/84/9 | 154 | 107.297 | 4.310998e-09 |

## Fraction mechanism screen (seed 601)

All density norms are Bohr^-3/2, energies and tolerances Ha, angles radians. Low G is the coarse orbital sphere; outside-density means beyond the 4Ec coarse density sphere. Raw histories contain every Pulay coefficient, conditioning and predicted norm.

| System | f | C/F/V maps | Coarse residual | Coarse last standard error (Ha) | First fine tolerance (Ha) | H-eq | abs dE (Ha) |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| si32 | 0.5 | 7/10/1 | 2.67838e-05 | 1.22765e-04 | 9.59099e-08 | 2365.66 | 1.058e-10 |
| si32 | 0.45 | 7/10/1 | 2.57779e-05 | 1.12548e-04 | 8.79280e-08 | 2373.13 | 1.702e-10 |
| si32 | 0.4 | 6/13/2 | 9.94947e-05 | 3.34348e-03 | 2.61209e-06 | 2499.15 | 2.068e-11 |
| si32 | 0.35 | 7/11/2 | 2.23164e-05 | 8.58090e-05 | 6.70383e-08 | 2711.35 | 1.448e-11 |
| si64 | 0.5 | 6/13/2 | 8.75574e-05 | 1.09985e-02 | 4.29628e-06 | 5474.45 | 5.312e-11 |
| si64 | 0.45 | 6/13/2 | 8.65231e-05 | 1.06416e-02 | 4.15688e-06 | 5322.15 | 2.754e-11 |
| si64 | 0.4 | 6/12/2 | 8.35163e-05 | 9.87347e-03 | 3.85682e-06 | 4552.46 | 3.632e-11 |
| si64 | 0.35 | 6/13/3 | 8.21680e-05 | 9.79398e-03 | 3.82577e-06 | 5292.19 | 7.645e-12 |
| si8 | 0.5 | 5/11/2 | 9.83239e-05 | 4.88965e-05 | 1.52802e-07 | 610.71 | 8.689e-11 |
| si8 | 0.45 | 5/10/2 | 9.39800e-05 | 4.43859e-05 | 1.38706e-07 | 579.04 | 1.317e-10 |
| si8 | 0.4 | 5/10/2 | 9.19423e-05 | 4.28459e-05 | 1.33893e-07 | 573.33 | 4.999e-11 |
| si8 | 0.35 | 5/9/2 | 8.45393e-05 | 3.76800e-05 | 1.17750e-07 | 573.33 | 8.359e-11 |
| si16 | 0.5 | 10/10/2 | 5.09752e-05 | 3.10411e-04 | 4.85017e-07 | 1369.80 | 6.935e-11 |
| si16 | 0.45 | 10/10/2 | 5.17973e-05 | 3.27070e-04 | 5.11047e-07 | 1365.77 | 1.214e-10 |
| si16 | 0.4 | 10/11/2 | 5.09281e-05 | 3.09325e-04 | 4.83320e-07 | 1391.72 | 1.705e-11 |
| si16 | 0.35 | 10/10/2 | 4.82117e-05 | 3.10813e-04 | 4.85646e-07 | 1319.61 | 1.191e-10 |
| si-skew | 0.5 | 11/8/3 | 8.64832e-05 | 1.45411e-06 | 1.81764e-08 | 202.16 | 2.639e-12 |
| si-skew | 0.45 | 11/9/2 | 8.84925e-05 | 1.28363e-06 | 1.60453e-08 | 202.54 | 8.148e-12 |
| si-skew | 0.4 | 13/7/2 | 7.30132e-05 | 9.13473e-07 | 1.14184e-08 | 176.51 | 1.353e-11 |
| si-skew | 0.35 | 13/8/3 | 5.32306e-05 | 4.88608e-07 | 6.10760e-09 | 201.38 | 6.875e-12 |
| mgo8 | 0.5 | 8/8/1 | 6.79212e-05 | 1.60193e-05 | 2.50301e-08 | 1076.70 | 3.263e-11 |
| mgo8 | 0.45 | 8/8/1 | 4.47140e-05 | 1.15515e-05 | 1.80492e-08 | 1079.15 | 4.798e-10 |
| mgo8 | 0.4 | 7/8/1 | 7.78201e-05 | 3.77335e-05 | 5.89586e-08 | 1013.86 | 6.337e-10 |
| mgo8 | 0.35 | 8/8/1 | 4.64974e-05 | 1.13970e-05 | 1.78078e-08 | 1024.59 | 1.052e-09 |
| syn-skew2 | 0.5 | 10/5/1 | 3.58211e-05 | 4.22314e-08 | 1.05578e-09 | 87.96 | 8.025e-10 |
| syn-skew2 | 0.45 | 8/6/1 | 8.45110e-05 | 2.48726e-07 | 6.21814e-09 | 81.81 | 7.066e-11 |
| syn-skew2 | 0.4 | 15/6/1 | 3.34101e-05 | 4.67739e-08 | 1.16935e-09 | 94.91 | 2.614e-10 |
| syn-skew2 | 0.35 | 15/6/1 | 3.34101e-05 | 4.67739e-08 | 1.16935e-09 | 94.91 | 2.613e-10 |
| syn-slab | 0.5 | 8/6/1 | 3.87128e-05 | 6.48868e-08 | 1.62217e-09 | 99.94 | 3.805e-10 |
| syn-slab | 0.45 | 8/6/1 | 3.87128e-05 | 6.48868e-08 | 1.62217e-09 | 99.94 | 3.805e-10 |
| syn-slab | 0.4 | 9/7/1 | 5.25670e-05 | 2.28339e-07 | 5.70847e-09 | 112.30 | 1.914e-10 |
| syn-slab | 0.35 | 7/6/1 | 2.25624e-05 | 8.88863e-08 | 2.22216e-09 | 98.20 | 2.586e-09 |

### Handoff state compared with the strict solution

| System | f | Density error | Low-G error | Outside-density error | Orbital energy error (Ha) | First-map energy error (Ha) | Max angle (rad) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| si32 | 0.5 | 3.20418e-04 | 2.51250e-04 | 2.54674e-06 | 1.45068e-01 | 1.06579e-04 | 1.99596e-02 |
| si32 | 0.45 | 4.50864e-04 | 3.42300e-04 | 5.95776e-06 | 2.30727e-01 | 1.85289e-04 | 2.67773e-02 |
| si32 | 0.4 | 6.23123e-04 | 4.09451e-04 | 1.49572e-05 | 3.54735e-01 | 3.67129e-04 | 4.22340e-02 |
| si32 | 0.35 | 8.69836e-04 | 5.59379e-04 | 3.07624e-05 | 5.46380e-01 | 5.98528e-04 | 5.15017e-02 |
| si64 | 0.5 | 2.22117e-04 | 1.71467e-04 | 1.79839e-06 | 2.90088e-01 | 4.03177e-04 | 2.02116e-02 |
| si64 | 0.45 | 3.20028e-04 | 2.43939e-04 | 4.18622e-06 | 4.73862e-01 | 5.34228e-04 | 2.69707e-02 |
| si64 | 0.4 | 4.42520e-04 | 2.90330e-04 | 1.05132e-05 | 7.16321e-01 | 7.84839e-04 | 4.25825e-02 |
| si64 | 0.35 | 6.07269e-04 | 3.89760e-04 | 2.16630e-05 | 1.08346e+00 | 1.29453e-03 | 5.20195e-02 |
| si8 | 0.5 | 6.59676e-04 | 5.18355e-04 | 6.31760e-06 | 3.85126e-02 | 2.17072e-05 | 1.99660e-02 |
| si8 | 0.45 | 8.20963e-04 | 5.72805e-04 | 1.30330e-05 | 5.11402e-02 | 3.21167e-05 | 2.45061e-02 |
| si8 | 0.4 | 1.23640e-03 | 8.02037e-04 | 3.01930e-05 | 8.76034e-02 | 6.68504e-05 | 3.52277e-02 |
| si8 | 0.35 | 1.70777e-03 | 1.12028e-03 | 5.99850e-05 | 1.34172e-01 | 1.29737e-04 | 5.11853e-02 |
| si16 | 0.5 | 4.67199e-04 | 3.73949e-04 | 3.87853e-06 | 7.52985e-02 | 1.39971e-04 | 1.99455e-02 |
| si16 | 0.45 | 6.21439e-04 | 4.62353e-04 | 8.69618e-06 | 1.10561e-01 | 1.80669e-04 | 2.66647e-02 |
| si16 | 0.4 | 8.90153e-04 | 5.86089e-04 | 2.12600e-05 | 1.77033e-01 | 2.59916e-04 | 4.16030e-02 |
| si16 | 0.35 | 1.23722e-03 | 8.06511e-04 | 4.32420e-05 | 2.74595e-01 | 4.04701e-04 | 5.13894e-02 |
| si-skew | 0.5 | 1.51306e-03 | 7.05010e-04 | 2.87751e-05 | 1.25918e-02 | 4.48189e-06 | 1.89812e-02 |
| si-skew | 0.45 | 2.33138e-03 | 1.00645e-03 | 4.42832e-05 | 2.06133e-02 | 1.19663e-05 | 2.86725e-02 |
| si-skew | 0.4 | 3.10713e-03 | 1.18296e-03 | 7.86570e-05 | 2.91212e-02 | 2.32587e-05 | 3.54546e-02 |
| si-skew | 0.35 | 3.64208e-03 | 1.23376e-03 | 1.49441e-04 | 3.71661e-02 | 3.13086e-05 | 4.55899e-02 |
| mgo8 | 0.5 | 1.67269e-02 | 1.48516e-02 | 6.83116e-04 | 4.89900e+00 | 1.73041e-01 | 1.22554e-01 |
| mgo8 | 0.45 | 2.27318e-02 | 1.96700e-02 | 1.35335e-03 | 7.79143e+00 | 3.56485e-01 | 1.80680e-01 |
| mgo8 | 0.4 | 2.98624e-02 | 2.64177e-02 | 2.32413e-03 | 1.18365e+01 | 6.75777e-01 | 2.33076e-01 |
| mgo8 | 0.35 | 3.74637e-02 | 3.20341e-02 | 3.98812e-03 | 1.68196e+01 | 1.15282e+00 | 2.72844e-01 |
| syn-skew2 | 0.5 | 8.38205e-03 | 4.42048e-03 | 2.64278e-04 | 2.32271e-02 | 3.27170e-05 | 6.72874e-02 |
| syn-skew2 | 0.45 | 8.45390e-03 | 4.73656e-03 | 1.23718e-03 | 2.51769e-02 | 3.77441e-05 | 7.56731e-02 |
| syn-skew2 | 0.4 | 1.10035e-02 | 6.42538e-03 | 1.37586e-03 | 3.30103e-02 | 7.60052e-05 | 1.13944e-01 |
| syn-skew2 | 0.35 | 1.10035e-02 | 6.42538e-03 | 2.21578e-03 | 3.30103e-02 | 7.60052e-05 | 1.13944e-01 |
| syn-slab | 0.5 | 1.03246e-02 | 5.57452e-03 | 3.83613e-04 | 5.94979e-02 | 1.41784e-03 | 1.24047e-01 |
| syn-slab | 0.45 | 1.03246e-02 | 5.57452e-03 | 6.05910e-04 | 5.94979e-02 | 1.41784e-03 | 1.24047e-01 |
| syn-slab | 0.4 | 1.55962e-02 | 9.35374e-03 | 1.80138e-03 | 1.23129e-01 | 1.75720e-02 | 1.52889e-01 |
| syn-slab | 0.35 | 3.52170e-02 | 2.18868e-02 | 3.46392e-03 | 3.04445e-01 | 4.07541e-02 | 3.04753e-01 |

### Handoff potential response

The high-PW columns restore the high-PW part of the density error to its reference value and measure the actual nonlinear potential change. They are not bounds and are unavailable if that diagnostic density is negative.

| System | f | Total Hartree RMS (Ha) | High-PW Hartree RMS (Ha) | Total XC RMS (Ha) | High-PW XC RMS (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: |
| si32 | 0.5 | 6.71402e-04 | 2.16906e-04 | 2.91032e-03 | 1.71686e-03 |
| si32 | 0.45 | 1.01223e-03 | 3.95231e-04 | 3.97207e-03 | 2.86982e-03 |
| si32 | 0.4 | 1.52789e-03 | 7.64213e-04 | 5.22101e-03 | 4.04661e-03 |
| si32 | 0.35 | 2.36635e-03 | 1.16783e-03 | 6.88258e-03 | 4.23884e-03 |
| si64 | 0.5 | 7.51799e-04 | 2.18935e-04 | 2.85762e-03 | 1.71326e-03 |
| si64 | 0.45 | 1.04960e-03 | 3.94231e-04 | 3.97123e-03 | 2.83698e-03 |
| si64 | 0.4 | 1.52922e-03 | 7.68909e-04 | 5.23579e-03 | 4.04462e-03 |
| si64 | 0.35 | 2.29705e-03 | 1.15417e-03 | 6.80774e-03 | 4.18394e-03 |
| si8 | 0.5 | 6.28648e-04 | 2.14926e-04 | 2.97134e-03 | 1.89549e-03 |
| si8 | 0.45 | 8.27436e-04 | unavailable | 3.67958e-03 | unavailable |
| si8 | 0.4 | 1.48099e-03 | 7.65806e-04 | 5.10744e-03 | 4.55588e-03 |
| si8 | 0.35 | 2.36978e-03 | 1.11258e-03 | 6.71714e-03 | 4.16803e-03 |
| si16 | 0.5 | 9.76882e-04 | 2.11533e-04 | 2.97819e-03 | 1.75134e-03 |
| si16 | 0.45 | 1.23007e-03 | 3.94893e-04 | 3.89730e-03 | 3.03450e-03 |
| si16 | 0.4 | 1.78998e-03 | 7.72526e-04 | 5.24902e-03 | 4.18290e-03 |
| si16 | 0.35 | 2.64408e-03 | 1.15810e-03 | 6.90710e-03 | 4.23061e-03 |
| si-skew | 0.5 | 8.04065e-04 | 4.28624e-04 | 2.88025e-03 | 2.63280e-03 |
| si-skew | 0.45 | 1.42911e-03 | 7.95267e-04 | 4.14545e-03 | 3.65137e-03 |
| si-skew | 0.4 | 2.18071e-03 | 1.21861e-03 | 5.38408e-03 | 4.32454e-03 |
| si-skew | 0.35 | 2.56851e-03 | 1.51556e-03 | 6.22500e-03 | 5.78069e-03 |
| mgo8 | 0.5 | 4.51957e-02 | unavailable | 3.61951e-02 | unavailable |
| mgo8 | 0.45 | 6.65789e-02 | unavailable | 3.90529e-02 | unavailable |
| mgo8 | 0.4 | 9.47532e-02 | unavailable | 4.65743e-02 | unavailable |
| mgo8 | 0.35 | 1.28748e-01 | unavailable | 5.58398e-02 | unavailable |
| syn-skew2 | 0.5 | 9.95667e-03 | unavailable | 1.07405e-02 | unavailable |
| syn-skew2 | 0.45 | 1.09343e-02 | unavailable | 1.20376e-02 | unavailable |
| syn-skew2 | 0.4 | 1.62195e-02 | unavailable | 1.81012e-02 | unavailable |
| syn-skew2 | 0.35 | 1.62195e-02 | unavailable | 1.81012e-02 | unavailable |
| syn-slab | 0.5 | 1.83145e-02 | 6.46851e-03 | 1.60757e-02 | 1.66345e-02 |
| syn-slab | 0.45 | 1.83145e-02 | 6.46851e-03 | 1.60757e-02 | 1.66345e-02 |
| syn-slab | 0.4 | 6.02555e-02 | unavailable | 2.29234e-02 | unavailable |
| syn-slab | 0.35 | 9.40023e-02 | unavailable | 4.52306e-02 | unavailable |

### Fine-map trajectories: si32 and si64

Reset is the cumulative feasibility-reset count; all arms start with empty fine history. Cpred is the predicted Coulomb norm ratio of the current L2 Pulay proposal, measured with the inert trace switch. `none` means no mixed proposal or no difference yet. All other systems have the same fields in their raw `fine-history.json` files.

| System | f | Map | Standard error (Ha) | Eigen tolerance (Ha) | H | Davidson | Density residual | Resets | Pulay differences | Cpred |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| si32 | 0.5 | 1 | 7.36224e-04 | 9.59099e-08 | 403 | 7 | 3.11070e-04 | 0 | 0 | none |
| si32 | 0.5 | 2 | 3.17122e-04 | 9.59099e-08 | 151 | 4 | 1.60489e-04 | 0 | 1 | 1.6083 |
| si32 | 0.5 | 3 | 3.95034e-03 | 9.59099e-08 | 239 | 9 | 7.95217e-05 | 0 | 2 | 0.035975 |
| si32 | 0.5 | 4 | 7.55681e-06 | 9.59099e-08 | 210 | 6 | 4.70393e-06 | 0 | 3 | 0.42674 |
| si32 | 0.5 | 5 | 1.58710e-06 | 5.90376e-09 | 122 | 4 | 2.50257e-06 | 0 | 4 | 0.33827 |
| si32 | 0.5 | 6 | 1.74871e-07 | 1.23992e-09 | 112 | 4 | 7.93498e-07 | 0 | 5 | 0.30869 |
| si32 | 0.5 | 7 | 2.69936e-08 | 1.36618e-10 | 105 | 3 | 2.91930e-07 | 0 | 6 | 0.69976 |
| si32 | 0.5 | 8 | 5.61931e-08 | 2.10887e-11 | 149 | 4 | 3.12935e-07 | 0 | 6 | 0.21016 |
| si32 | 0.5 | 9 | 7.54678e-09 | 2.10887e-11 | 128 | 4 | 1.26530e-07 | 0 | 6 | 0.14403 |
| si32 | 0.5 | 10 | 1.07032e-10 | 5.89592e-12 | 98 | 3 | 2.13896e-08 | 0 | none | none |
| si32 | 0.45 | 1 | 1.59369e-03 | 8.79280e-08 | 412 | 7 | 4.37453e-04 | 0 | 0 | none |
| si32 | 0.45 | 2 | 5.38001e-04 | 8.79280e-08 | 161 | 4 | 2.25169e-04 | 0 | 1 | 1.2442 |
| si32 | 0.45 | 3 | 4.02208e-03 | 8.79280e-08 | 245 | 9 | 8.01169e-05 | 0 | 2 | 0.043599 |
| si32 | 0.45 | 4 | 7.61077e-06 | 8.79280e-08 | 214 | 6 | 5.26187e-06 | 0 | 3 | 0.71503 |
| si32 | 0.45 | 5 | 5.55538e-06 | 5.94592e-09 | 138 | 4 | 4.11728e-06 | 0 | 4 | 0.19148 |
| si32 | 0.45 | 6 | 2.11725e-07 | 4.34014e-09 | 108 | 4 | 8.49722e-07 | 0 | 5 | 0.32102 |
| si32 | 0.45 | 7 | 2.51420e-08 | 1.65410e-10 | 119 | 4 | 3.37250e-07 | 0 | 6 | 0.81234 |
| si32 | 0.45 | 8 | 7.41260e-08 | 1.96422e-11 | 151 | 5 | 3.57759e-07 | 0 | 6 | 0.26066 |
| si32 | 0.45 | 9 | 1.77353e-08 | 1.96422e-11 | 142 | 4 | 1.88788e-07 | 0 | 6 | 0.09475 |
| si32 | 0.45 | 10 | 1.02933e-10 | 1.38557e-11 | 92 | 3 | 2.21004e-08 | 0 | none | none |
| si32 | 0.4 | 1 | 3.34579e-03 | 2.61209e-06 | 339 | 6 | 6.04166e-04 | 0 | 0 | none |
| si32 | 0.4 | 2 | 9.91611e-04 | 2.61209e-06 | 68 | 2 | 3.09662e-04 | 0 | 1 | 1.1195 |
| si32 | 0.4 | 3 | 4.04397e-03 | 7.74696e-07 | 196 | 10 | 8.07929e-05 | 0 | 2 | 0.083556 |
| si32 | 0.4 | 4 | 1.46808e-04 | 7.74696e-07 | 154 | 6 | 1.72811e-05 | 0 | 3 | 0.19632 |
| si32 | 0.4 | 5 | 1.01649e-05 | 1.14693e-07 | 97 | 3 | 6.12087e-06 | 0 | 4 | 0.35314 |
| si32 | 0.4 | 6 | 1.20781e-06 | 7.94133e-09 | 115 | 4 | 2.45877e-06 | 0 | 5 | 0.5032 |
| si32 | 0.4 | 7 | 4.43842e-07 | 9.43602e-10 | 122 | 4 | 1.14640e-06 | 0 | 6 | 0.55064 |
| si32 | 0.4 | 8 | 6.02224e-07 | 3.46752e-10 | 135 | 4 | 1.03681e-06 | 0 | 6 | 0.16691 |
| si32 | 0.4 | 9 | 3.35570e-08 | 3.46752e-10 | 110 | 4 | 3.01289e-07 | 0 | 6 | 0.26759 |
| si32 | 0.4 | 10 | 4.18839e-09 | 2.62164e-11 | 103 | 3 | 1.15280e-07 | 0 | 6 | 0.036057 |
| si32 | 0.4 | 11 | 3.53586e-10 | 3.27218e-12 | 185 | 5 | 3.42910e-08 | 0 | 6 | 0.63086 |
| si32 | 0.4 | 12 | 2.23830e-10 | 2.76239e-13 | 141 | 4 | 2.61970e-08 | 0 | 6 | 0.27321 |
| si32 | 0.4 | 13 | 2.73343e-11 | 1.74867e-13 | 107 | 4 | 9.59367e-09 | 0 | none | none |
| si32 | 0.35 | 1 | 7.19635e-03 | 6.70383e-08 | 438 | 7 | 8.48410e-04 | 0 | 0 | none |
| si32 | 0.35 | 2 | 1.85192e-03 | 6.70383e-08 | 194 | 4 | 4.32638e-04 | 0 | 1 | 0.56424 |
| si32 | 0.35 | 3 | 2.39969e-03 | 6.70383e-08 | 253 | 10 | 6.34517e-05 | 0 | 2 | 0.19814 |
| si32 | 0.35 | 4 | 2.10515e-04 | 6.70383e-08 | 217 | 6 | 2.10323e-05 | 0 | 3 | 0.18758 |
| si32 | 0.35 | 5 | 1.18611e-05 | 6.70383e-08 | 125 | 4 | 5.77592e-06 | 0 | 4 | 0.13093 |
| si32 | 0.35 | 6 | 1.67175e-07 | 9.26646e-09 | 109 | 4 | 9.16667e-07 | 0 | 5 | 0.42823 |
| si32 | 0.35 | 7 | 1.99726e-08 | 1.30605e-10 | 149 | 5 | 3.82736e-07 | 0 | 6 | 0.62339 |
| si32 | 0.35 | 8 | 3.09526e-08 | 1.56036e-11 | 138 | 4 | 2.62036e-07 | 0 | 6 | 0.57835 |
| si32 | 0.35 | 9 | 4.57217e-08 | 1.56036e-11 | 145 | 4 | 2.79610e-07 | 0 | 6 | 0.077433 |
| si32 | 0.35 | 10 | 3.24217e-10 | 1.56036e-11 | 120 | 4 | 3.30610e-08 | 0 | 6 | 0.047966 |
| si32 | 0.35 | 11 | 1.86347e-11 | 2.53295e-13 | 179 | 5 | 8.87935e-09 | 0 | none | none |
| si64 | 0.5 | 1 | 1.57219e-03 | 4.29628e-06 | 593 | 5 | 2.15175e-04 | 0 | 0 | none |
| si64 | 0.5 | 2 | 1.17775e-03 | 6.14138e-07 | 232 | 6 | 1.12586e-04 | 0 | 1 | 2.0326 |
| si64 | 0.5 | 3 | 1.93785e-02 | 4.60060e-07 | 416 | 14 | 9.15383e-05 | 0 | 2 | 0.074597 |
| si64 | 0.5 | 4 | 2.31077e-04 | 4.60060e-07 | 380 | 6 | 1.09976e-05 | 0 | 3 | 0.088229 |
| si64 | 0.5 | 5 | 2.89141e-06 | 9.02646e-08 | 207 | 4 | 2.07740e-06 | 0 | 4 | 0.36151 |
| si64 | 0.5 | 6 | 4.93156e-07 | 1.12946e-09 | 302 | 5 | 1.13890e-06 | 0 | 5 | 0.54936 |
| si64 | 0.5 | 7 | 2.08882e-07 | 1.92639e-10 | 217 | 4 | 7.05082e-07 | 0 | 6 | 0.59305 |
| si64 | 0.5 | 8 | 2.18298e-07 | 8.15947e-11 | 289 | 4 | 3.59918e-07 | 0 | 6 | 0.22113 |
| si64 | 0.5 | 9 | 3.96946e-08 | 8.15947e-11 | 204 | 3 | 1.51097e-07 | 0 | 6 | 0.27665 |
| si64 | 0.5 | 10 | 7.79579e-09 | 1.55057e-11 | 210 | 4 | 8.31359e-08 | 0 | 6 | 0.41719 |
| si64 | 0.5 | 11 | 2.37027e-09 | 3.04523e-12 | 208 | 4 | 4.10857e-08 | 0 | 6 | 0.052705 |
| si64 | 0.5 | 12 | 3.18791e-10 | 9.25887e-13 | 378 | 6 | 1.46577e-08 | 0 | 6 | 0.38996 |
| si64 | 0.5 | 13 | 7.73817e-11 | 1.24528e-13 | 235 | 4 | 7.81596e-09 | 0 | none | none |
| si64 | 0.45 | 1 | 3.31626e-03 | 4.15688e-06 | 630 | 5 | 3.09708e-04 | 0 | 0 | none |
| si64 | 0.45 | 2 | 1.48046e-03 | 1.29541e-06 | 194 | 6 | 1.60410e-04 | 0 | 1 | 1.7459 |
| si64 | 0.45 | 3 | 1.63483e-02 | 5.78305e-07 | 397 | 14 | 8.49348e-05 | 0 | 2 | 0.08059 |
| si64 | 0.45 | 4 | 2.83055e-04 | 5.78305e-07 | 352 | 6 | 1.18800e-05 | 0 | 3 | 0.087045 |
| si64 | 0.45 | 5 | 4.32213e-06 | 1.10568e-07 | 203 | 3 | 2.62002e-06 | 0 | 4 | 0.41571 |
| si64 | 0.45 | 6 | 3.91991e-07 | 1.68833e-09 | 296 | 5 | 1.31448e-06 | 0 | 5 | 1.0261 |
| si64 | 0.45 | 7 | 5.93964e-07 | 1.53121e-10 | 308 | 5 | 7.86983e-07 | 0 | 6 | 0.41585 |
| si64 | 0.45 | 8 | 3.50230e-07 | 1.53121e-10 | 279 | 4 | 4.33699e-07 | 0 | 6 | 0.30207 |
| si64 | 0.45 | 9 | 1.15498e-07 | 1.36809e-10 | 233 | 4 | 2.33676e-07 | 0 | 6 | 0.17644 |
| si64 | 0.45 | 10 | 3.36263e-09 | 4.51166e-11 | 203 | 4 | 6.62418e-08 | 0 | 6 | 0.085692 |
| si64 | 0.45 | 11 | 4.88550e-10 | 1.31353e-12 | 373 | 6 | 2.04456e-08 | 0 | 6 | 0.55274 |
| si64 | 0.45 | 12 | 2.79423e-10 | 1.90840e-13 | 258 | 4 | 1.43294e-08 | 0 | 6 | 0.29982 |
| si64 | 0.45 | 13 | 3.75603e-11 | 1.09150e-13 | 218 | 4 | 5.68772e-09 | 0 | none | none |
| si64 | 0.4 | 1 | 6.76004e-03 | 3.85682e-06 | 652 | 6 | 4.28979e-04 | 0 | 0 | none |
| si64 | 0.4 | 2 | 2.02730e-03 | 2.64064e-06 | 132 | 2 | 2.20048e-04 | 0 | 1 | 1.1969 |
| si64 | 0.4 | 3 | 7.61330e-03 | 7.91915e-07 | 371 | 15 | 5.76776e-05 | 0 | 2 | 0.10954 |
| si64 | 0.4 | 4 | 4.13940e-04 | 7.91915e-07 | 297 | 6 | 1.42635e-05 | 0 | 3 | 0.16282 |
| si64 | 0.4 | 5 | 2.60275e-05 | 1.61695e-07 | 196 | 3 | 5.32227e-06 | 0 | 4 | 0.3019 |
| si64 | 0.4 | 6 | 2.59937e-06 | 1.01670e-08 | 225 | 4 | 2.00276e-06 | 0 | 5 | 0.45319 |
| si64 | 0.4 | 7 | 6.69297e-07 | 1.01538e-09 | 219 | 4 | 8.27819e-07 | 0 | 6 | 0.41212 |
| si64 | 0.4 | 8 | 3.32044e-07 | 2.61444e-10 | 248 | 4 | 4.62571e-07 | 0 | 6 | 0.34733 |
| si64 | 0.4 | 9 | 1.19990e-07 | 1.29705e-10 | 234 | 4 | 2.49072e-07 | 0 | 6 | 0.26005 |
| si64 | 0.4 | 10 | 1.26978e-08 | 4.68711e-11 | 216 | 4 | 1.00498e-07 | 0 | 6 | 0.32327 |
| si64 | 0.4 | 11 | 5.16713e-10 | 4.96008e-12 | 209 | 3 | 2.81582e-08 | 0 | 6 | 0.11706 |
| si64 | 0.4 | 12 | 5.27644e-11 | 2.01841e-13 | 335 | 5 | 8.45376e-09 | 0 | none | none |
| si64 | 0.35 | 1 | 1.39417e-02 | 3.82577e-06 | 685 | 6 | 5.91403e-04 | 0 | 0 | none |
| si64 | 0.35 | 2 | 3.64187e-03 | 3.82577e-06 | 132 | 2 | 3.00596e-04 | 0 | 1 | 0.83498 |
| si64 | 0.35 | 3 | 5.43691e-03 | 1.42261e-06 | 307 | 4 | 4.98394e-05 | 0 | 2 | 0.19923 |
| si64 | 0.35 | 4 | 8.92591e-04 | 1.42261e-06 | 270 | 6 | 2.06878e-05 | 0 | 3 | 0.2264 |
| si64 | 0.35 | 5 | 6.15666e-05 | 3.48668e-07 | 212 | 4 | 7.30821e-06 | 0 | 4 | 0.20888 |
| si64 | 0.35 | 6 | 1.20093e-06 | 2.40494e-08 | 223 | 4 | 1.62762e-06 | 0 | 5 | 0.6297 |
| si64 | 0.35 | 7 | 3.81503e-07 | 4.69114e-10 | 271 | 5 | 6.88391e-07 | 0 | 6 | 0.48533 |
| si64 | 0.35 | 8 | 1.82136e-07 | 1.49025e-10 | 244 | 4 | 3.82273e-07 | 0 | 6 | 0.45855 |
| si64 | 0.35 | 9 | 1.12895e-07 | 7.11470e-11 | 248 | 4 | 2.39672e-07 | 0 | 6 | 0.37173 |
| si64 | 0.35 | 10 | 4.25656e-08 | 4.40998e-11 | 245 | 4 | 1.49515e-07 | 0 | 6 | 0.176 |
| si64 | 0.35 | 11 | 6.87809e-10 | 1.66272e-11 | 200 | 3 | 2.88994e-08 | 0 | 6 | 0.075003 |
| si64 | 0.35 | 12 | 1.97528e-10 | 2.68675e-13 | 353 | 5 | 1.22483e-08 | 0 | 6 | 0.046825 |
| si64 | 0.35 | 13 | 1.16009e-11 | 7.71593e-14 | 371 | 5 | 3.89410e-09 | 0 | none | none |

## Validation and reproduction

All scratch is in this worktree's `target/scratch`. The executable is copied
into each fresh batch directory before any geometry is run; later builds
cannot replace it mid-batch. No UPF file is copied or committed. Each batch
has input files, PP SHA-256, source and binary SHA-256, exit records,
measurement JSON/CSV, per-map histories and binary handoff states. Reference
options and the seven fixtures are inherited unchanged from the documented
F-LADDER-HANDOFF harness. The matrix is pure Rust, Gamma, one core, fixed
occupations: it is not MPI evidence, nor a metallic/k-point generality claim.

Actual successful checks:

- Release unit `fraction_rules_preserve_grid_and_parseval_budget`: three
  cell shapes, five requested fractions; grid preservation, nested basis,
  idempotent grid choice, off identity and Parseval residual-budget relation.
- Release unit `fraction_watch_is_bit_identical_and_forecast_is_diagnostic`:
  exact delivered energy, density and orbital equality; known contraction,
  forecast excess and unavailable cases. Both new tests pass in the fast
  debug gate too (neither is a slow test).
- 21 baseline offline/production controls compare delivered density, energy
  bits, coarse H and fine map counts exactly. Every grid/modes/combined/watch
  arm also checks the same fields against the production `run_ladder` path
  with its requested configuration (84 additional controls).
- 21 watch/off pairs: exact energies, H, map counts, scalar histories and
  handoff binaries; 14 trace/off pairs across all seven systems: exact
  energy, reference energy, stage H, map counts, coarse schedules and all
  four handoff/coarse array files. The pre-rule si32 executable's four
  fractions also match the final off path in those fields and both complete
  scalar histories. This is not a cross-platform bit-identity claim.
- The Python reader checks density Parseval decomposition and angle ranges
  for all 126 matrix arms. The mechanism reader reconciles H sums and map
  counts on 32 arms / **288 fine maps**. All standard arms passed the unchanged
  fixed-point verification. All references completed the strict certificate.
- Required `scripts/check.sh fast dft-engine dft-numerics dft-cli`: **exit 1**;
  build/fmt/clippy pass, **468 passed / 3 failed / 134 ignored**, 51 binaries,
  259 s loaded-Mac wall. Failures are
  `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`,
  `pseudo::cache_resolution_and_no_network_side_effects`, and
  `pseudo::cached_fetch_checks_hash_and_provenance`. Each rejects the required
  in-repository TMPDIR (`pseudopotential cache must be outside repositories`).
  An isolated rerun of the first exits 101 with the same reason. The load
  reporter also cannot execute sandboxed `ps`; the four gates still run.
- Python compile checks and `git diff --check`: exit 0. No test tolerance was
  relaxed. README generators, dashboard JSON, task titles and the master
  session registry are outside this ticket's `OWNED_FILES`; no such files
  were modified. No default/performance result was adopted for README.

The first exploratory `fraction-mechanism` batch preceded executable
snapshotting; builds changed its path during the batch. Only its first
si32 job is used for the pre-rule comparison. All reported mechanism tables
use the separately rerun immutable `fraction-mechanism-final` batch. Timing
is deliberately unjudged: concurrent local builds/experiments and profiling
make these walls unsuitable for performance acceptance. New rule setup costs
are included in each primary calculation wall but are not assigned fictitious
H columns. Offline quality probes and exact-driver controls are outside that
wall and outside the reported candidate H counts.

```sh
export TMPDIR_990="$PWD/target/scratch"
export TMPDIR="$TMPDIR_990"
mkdir -p "$TMPDIR_990"
cargo test --release -p dft-engine --lib --no-run
# Set BIN to the executable printed by cargo (in this run: the path below).
BIN=target/release/deps/dft_engine-5a2965bdd204a9df
python3 -I scripts/ladder_handoff.py --binary "$BIN" \
  --out target/scratch/fraction-matrix \
  --cases si8,si16,si32,si-skew,mgo8,syn-skew2,syn-slab --seeds 601,602,603 \
  --arms '0.5:1e-4:default,0.4:1e-4:default,0.4:1e-4:grid,0.4:1e-4:modes,0.4:1e-4:watch,0.4:1e-4:grid-modes' --timeout 1200
python3 -I scripts/ladder_handoff.py --binary "$BIN" \
  --out target/scratch/fraction-mechanism-final \
  --cases si32,si64,si8,si16,si-skew,mgo8,syn-skew2,syn-slab --seeds 601 \
  --arms '0.5:1e-4:default,0.45:1e-4:default,0.4:1e-4:default,0.35:1e-4:default' \
  --trace-pulay --timeout 1800
python3 -I scripts/ladder_handoff.py --binary "$BIN" \
  --out target/scratch/fraction-ablation --cases si32 --seeds 601 \
  --arms '0.4:1e-4:default,0.4:1e-4:first-cap,0.4:1e-4:first-coarse,0.4:1e-4:first-tight,0.4:3e-5:default,0.45:3e-5:default' --trace-pulay --timeout 1200
python3 -I scripts/ladder_handoff.py --binary "$BIN" \
  --out target/scratch/fraction-factorial --cases si32 --seeds 601 \
  --arms '0.4:3e-5:default,0.4:1e-4:first-previous,0.4:3e-5:first-previous' \
  --trace-pulay --timeout 1200
python3 -I scripts/ladder_handoff.py --binary "$BIN" \
  --out target/scratch/fraction-si64-rules --cases si64 --seeds 601 \
  --arms '0.4:1e-4:grid,0.4:1e-4:modes,0.4:1e-4:grid-modes' --timeout 1800
python3 -I scripts/ladder_fraction_analyze.py \
  --rows target/scratch/fraction-matrix/rows.json --out target/scratch/fraction-analysis
python3 -I scripts/ladder_fraction_mechanism.py \
  --root target/scratch/fraction-mechanism-final --out target/scratch/fraction-analysis
CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli
```

Use fresh output directories on a repeat; the runner rejects nonempty ones.
Preserve the batch directories before eventual worktree cleanup. Reports here
are sufficient to review the negative result; the arrays and raw maps are
needed for any further orbital/spectral investigation.

## dev08 command lines for the orchestrator (not executed here)

No candidate passed the Mac generality condition, so this is a
**characterization request, not an adoption recommendation**. If the master
chooses to measure it, start with the grid rule and its counterexamples.
Compare BAANDI against BAANDI within one parallel model only. Use the same
binary, input bytes, UPF, fine grid, band count, start, thresholds, verification
and binding in A and B. A is the adopted 0.5; B requests 0.4 plus `grid`.
A separate `modes` comparison uses 0.4/off/modes. Do not silently stack them.

Prerequisites: current master-approved production inputs for seeds 601–603,
matching strict-reference energies, binary built with the intended backend,
UPF hashes, a fresh output tree, and the established dev08 signal protocol.
Obtain an A/A floor first and gate node load per the master protocol. At
least four crossed, order-alternated pairs per geometry; retain all runs and
exclude the first pair from timing statistics. Report the mean paired wall
change and its standard error, not a ratio of unrelated medians. Small or
unresolved changes need more pairs. Keep a separately labelled diagnostic
run with per-map histories and profiling; keep profiling/snapshots off for
timed runs. The fixed-point verification stays on and inside wall time.

```bash
# Run only through the master/dev08 session. Substitute approved absolute paths.
BIN=/approved/build/target/release/dft-cli
INPUT_ROOT=/approved/fraction-inputs  # <cell>-<seed>.toml, threads=1
OUT="$PWD/target/scratch/fraction-dev08-NEW"
mkdir -p "$OUT/scratch"
export TMPDIR_990="$OUT/scratch" TMPDIR="$OUT/scratch"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# Start from the manifest's clean, frozen BAANDI environment.
sha256sum "$BIN" > "$OUT/binary.sha256"
run_case() {
  local cell=$1 ranks=$2 model=$3 seed pair arm dest input
  local -a launch order
  if [ "$model" = one_core ]; then
    launch=(taskset -c 0)
  else
    launch=(taskset -c "0-$((ranks-1))" mpirun -np "$ranks"
            --map-by core --bind-to core --report-bindings)
  fi
  for seed in 601 602 603; do
    input="$INPUT_ROOT/$cell-$seed.toml"
    for pair in 0 1 2 3; do
      if (( (pair+seed)%2 )); then order=(B A); else order=(A B); fi
      for arm in "${order[@]}"; do
        dest="$OUT/$cell-$model-$ranks-$seed-$pair-$arm"
        mkdir -p "$dest"
        cp "$input" "$dest/input.toml"
        local fraction=0.5 rule=off
        if [ "$arm" = B ]; then fraction=0.4; rule=grid; fi
        /usr/bin/time -f '%e' -o "$dest/wall.s" env \
          BAANDI_CUTOFF_LADDER="$fraction" BAANDI_LADDER_FRACTION_RULE="$rule" \
          BAANDI_LADDER_STOP_RULE=off BAANDI_LADDER_WATCH=off \
          BAANDI_LADDER_COARSE_STOP=1e-4 BAANDI_LADDER_COARSE_RULE=fixed \
          BAANDI_LADDER_FIRST_TOL=auto BAANDI_LADDER_STANDARD_SEED=coarse \
          BAANDI_PHASE_PROFILE=off BAANDI_STAGE_TIMERS=off BAANDI_DAVIDSON_STAGES=off \
          "${launch[@]}" "$BIN" run "$input" --output-dir "$dest/calc" \
          > "$dest/stdout" 2> "$dest/stderr"
        local rc=$?
        echo "$rc" > "$dest/exit-code"
        [ "$rc" = 0 ] || return "$rc"
      done
    done
  done
}
# Submit each call as its own bounded signal job; do not run jobs concurrently.
run_case si64 8 mpi
run_case si64 1 one_core
run_case si32 8 mpi
run_case si32 1 one_core
run_case si16 8 mpi
run_case si8 4 mpi
run_case si-skew 4 mpi
# Add MgO and synthetic skew/slab counterexamples before any adoption decision.
```

Record load before/after, effective requested/effective fractions, grid and
PW counts, coarse/fine/verification H and maps, Davidson iterations, every
verification attempt, the per-map Pulay record, |E-Estrict| and density hash.
The master decides whether any further rule family is worth pursuing.
This experiment supports no automatic material-specific fallback or default
change. A future first-tolerance rule must account for the measured tradeoff:
fewer fine maps alone can mean **more** fine H.

## Stored provenance

Base commit: `032f6823788518966f8d223d66af8216153021c2`; changes deliberately uncommitted.

| Batch | Arms | Binary SHA-256 | Rows SHA-256 |
| --- | ---: | --- | --- |
| `fraction-matrix` | 126 | `4b4225ceeadca0b0f1ac84225e972224f346f77252ce87e84f446cd8a26347f8` | `6dcaac8a95906a0f622bc9a981e4492fe179026b7f0b798678ad549dd02300c8` |
| `fraction-mechanism-final` | 32 | `0cb215efd8e1990a10a40d9ccc5702a211ddac5497f069874632f50ea251bc89` | `30325ff1cf43a3efabf86ee977885a348071fbb3fcba93d3607415deb0678716` |
| `fraction-ablation` | 6 | `4b4225ceeadca0b0f1ac84225e972224f346f77252ce87e84f446cd8a26347f8` | `d54593059c78787f7f559df38b079800aec47e6d400f249a3b5eeafba0eae7a5` |
| `fraction-factorial` | 3 | `ca574900cd647cd95adc6173518119012c9acec831ed281f66634608f8dbfe7e` | `d8176c073248b85746f231d137478e7d1c4d512d495b588c423d21240124c287` |
| `fraction-si64-rules` | 3 | `0cb215efd8e1990a10a40d9ccc5702a211ddac5497f069874632f50ea251bc89` | `2a6e28c88fc1449dddd525d8c8304676dd76a81c51257f375f5a4c7ffddf2142` |

The binaries differ because test instrumentation and controls were added between batches;
the run-owned production rule implementation is unchanged across the candidate batches.
All PP hashes and exact fixture text are in each geometry directory.
