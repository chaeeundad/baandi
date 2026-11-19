# BAANDI vs QE — follow-up review of the 36 commits to 3effafef

19 November 2026 · CPU · One core and pure MPI · Engineering and numerical work order

**Decision.** Keep owned projector initialization and verification scheduling as adopted. Close the owned orbital transfer and minimum-smooth-grid experiments. Put the next numerical effort into a controlled test of the first fine SCF map, and the next systems effort into a replicated validation loop that survived the owned-projector change. Request a current matched BAANDI/QE pair immediately. This review found both an explanation for the reported 46-fold change in the inner tolerance and a cheap experiment that can distinguish it from a poorer starting density.

The review is pinned to [3effafef][head], 36 commits after [the previous review][delta]. Evidence below distinguishes relayed dev08 observations, independent calculations on published fixtures, and proposed experiments. I inspected the production paths, test driver, solver records and benchmark rules. I ran local arithmetic, reciprocal-lattice enumeration and fixture-loader checks. I did not run either DFT code, a dev08 job, or the Rust test suite. [Previous plan][previous] · [Comparison rules][rules]

## 1. Which changes earned their place

| Pure MPI, dev08; BAANDI candidate vs BAANDI control | Result | Interpretation |
|---|---|---|
| Si64 np8, owned projector initialization | **−11.36%** through convergence, −11.14% at ten steps. Reported complex table per rank 132.2 → 8.3 MB. Same SCF, Davidson and H counts. | Adopted. The measured wall includes any benefit the smaller table provides later. The records do not assign the whole gain to the constructor. |
| Si64 np8, verification scheduling | −0.18%, unresolved; retry schedule unchanged. | Zero demonstrated MPI gain on this case. |
| Si64 np8, direct rank-owned coarse→fine transfer | Original +3.12%. Revised version +1.27%, slower in 9/10 pairs, even as its own stage fell 0.097 → 0.032 s. | **Close.** Tree sums 13 → 0, allgatherv 9 → 1 in the revision. Padded exchange still sent 3.95 MB for 1.11 MB of useful data, and owned peak bytes rose 32%. These are possible costs, not a full attribution of the regression. |
| Si32 np8, revised transfer | −0.27%, unresolved. | No case for enabling it. Even the whole Si64 baseline handoff was only about a tenth of a second. |

[Ledger rows F-INIT-STATES-DEFAULT-ON, G-HANDOFF-OWNED-CLOSED and F-VERIFY-SCHEDULE-DEFAULT-ON][ledger] · [Current switch defaults][defaults]

| One core, dev08; BAANDI candidate vs BAANDI control | Result | Interpretation |
|---|---|---|
| Si64, verification scheduling | **−6.81%**. Attempts [11,12,13] → [11,13], H columns 4,008 → 3,871, energy difference 4.11e−11 Ha. | Adopted. A real certificate still passes. Count the changed continuation as a change in trajectory. |
| Si32, same switch, no failures | +0.29% on the adoption binary; 8/8 pairs slower. | The follow-up source audit found **zero extra allocations on ordinary maps even before its fix**. Final reporting allocations fell from two to one and first-success capacity from 100 to one. A recent dev08 result proving the 23 ms disappeared was not found. |

These are different experiments in different parallel models. None of these percentages is a BAANDI/QE ratio. The recent 51-binary, 466-pass/3-failure/134-ignored figure is a selected test invocation. The earlier inventory covered the workspace. No evidence here establishes a further workspace reduction from 81 to 51 or a passing latest full suite. The three failures are recorded as CLI pseudo-cache/TMPDIR cases. [Allocation audit][sched] · [Test inventory][inventory] · [Retry implementation][retry]

**Status of the previous recommendations.** Direct handoff has a decisive negative result. G-GUARD-GRAM has a ticket and a blocked session entry. The reviewed switch registry and mixer do not show a merged Gram guard or a dev08 A/B/C. Price the existing guard before funding its speedup. The occupied-plus-complementary reduced-density predictor remains an open hypothesis. I found no merged production switch or measured full-run win for it in this change set. The rejected high-G handoff *indicator* was a different test. [Gram ticket][gram] · [Session registry][sessions] · [Existing guard][mixer] · [High-G study][study]

## 2. Where the opportunity moved

Changing only the coarse energy fraction from 0.50 to 0.40 gave the following paired outcomes. Fine cutoff, fine FFT grid and final verification are fixed. This is an experimental arm. Signs below refer to BAANDI elapsed time; a negative sign is a saving.

| Pure MPI, dev08 | Si8 np4 | Si16 np8 | Si-skew np4 | Si32 np4 | Si32 np8 | Si64 np4 | Si64 np8 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Fraction 0.40 vs 0.50 | −8.43% | −7.41% | −2.88% | −0.81% unresolved | −2.33% | −18.55% | **−16.85%** |

| One core, dev08 | Si32 | Si64 |
|---|---:|---:|
| Fraction 0.40 vs 0.50 | **+9.71%** | **−12.83%** |

At Si32 MPI8, 0.45 saved **6.91%** with ten fine maps; 0.40 needed thirteen. At Si64 MPI8, fractions 0.50/0.45/0.40 gave 11/12/11 fine maps and 3,504/3,510/3,044 fine H columns. A global move to 0.40 would knowingly ship a one-core regression. The reported worst candidate-vs-default energy difference of about 2.4e−10 Ha fits the standard 5e−9 Ha budget and exceeds 1e−10 Ha. A difference between two runs cannot by itself establish the absolute error of either. Test strict-class candidates against independently tightened references. [Fraction sweep][fraction] · [Follow-up data][follow]

**Distance from QE.** The earlier diagnostic B = 8.8835 s and 0.90Q ≈ 4.16 s remain an old MPI8 budget. The new protocol shows same-binary Si64 MPI8 medians spanning 8.67–10.71 s in different windows. Subtracting or multiplying those times by later adoption percentages would manufacture a new scoreboard. A useful test independent of any old clock is that a 16.85% saving can meet the target by itself only if its **contemporaneous** baseline B/Q is already ≤0.90/0.8315 = **1.0824**. Get today's B, Q, B/Q and B−0.90Q in one window for each case, with profile off. Also compare SCF counts, H/FFT work definitions and exposed stage seconds separately. An untimed QE routine has unknown exclusive cost. Winning every required row this round has low demonstrated support. [Old absolute budget][launch] · [Si64 instability and ten-step limitation][protocol]

## 3. A causal test for the Si32 cliff

**The source explains the tolerance numerically.** The coarse exit compares a density-residual norm to the coarse stop; other exit reasons include a terminal status and the map limit. The seed for the first standard fine solve uses another quantity, the last coarse standard error in Hartree. With the usual seed, the code computes the smaller of the preceding requested tolerance and 0.1 × error / electron count, then clamps to its allowed range. The final coarse residual and the error in Hartree cannot be substituted for one another. For Si32, 128 valence electrons and errors 6.8e−5 / 3.15e−3 Ha give **5.31e−8 / 2.46e−6 Ha**, reproducing the recorded fine tolerances. The fine session can also redo the first solve and charge it to that map. The ratchet retains a tight request later, so changing the initial seed can affect several maps. [Coarse stop and seed][coarse] · [Fine seed selection][first] · [Ratchet][ratchet] · [Map-one redo][redo]

**The first map is a misleading target.** In the fraction comparison reported for Si32 MPI8, first-map H columns fell from 430 to 359, saving 71. Fine H totals rose from 1,856 to 1,900. Under those rows' accounting, later maps spent **115 extra H columns**, for a net loss of 44. Three extra maps also cost density construction, potential evaluation and mixing. The separate stop sweep in the same note quotes a baseline of 1,717 fine H. Recover the two raw manifests, histories, redo and verification inclusions before joining them. [Published map counts][follow]

**Request F-LADDER-FRACTION to run a 2×2 intervention before fitting a rule.** Capture the two natural handoff states at fractions 0.50 and 0.40. Each state includes its own input density and warm orbitals. Continue each twice, using the tight and loose scalar actually observed in the two natural runs. Keep the same fine basis, Hamiltonian construction, state-import semantics, solver, mixer and certificate. Within each state pair the initial fine Hamiltonian is identical. Between the two states its density-dependent potential can legitimately differ. Keep the original 0.45 trajectory as an external cross-check. Count effective requested tolerance, first redo, H and block calls per map, mixing resets/refits, per-map Coulomb residual, charge, rejected attempts and full verified seconds. Report the fine-only time and the full price of constructing the corresponding starting state. Diagonal arms must reproduce their unsplit production controls. Use an existing test injection point for an arbitrary seed; the literal values are counterfactual controls for this case and must never become production defaults. The current offline driver assigns `c.standard_next` unconditionally, so changing a CLI first-tolerance choice there is insufficient. [Offline continuation][driverseed] · [Production seed setter][seedsetter]

| What the crossed runs show | Decision |
|---|---|
| Tightening the 0.40 state recovers the missing maps and improves full wall, while loosening the 0.50 state reproduces most of the loss | Prioritize a handoff-aware tolerance using only information already available during the candidate run. Validate its additional inner work on Si64 and held-out systems. |
| The 0.40 state stays expensive at both tolerances | Focus on the representation or mixer history. Compare common-low-G density error and occupied-projector principal angles, including degeneracies. Keep the inexpensive natural tolerance. |
| There is a strong interaction or outcomes vary by geometry | Report all four cells and their paired differences. Keep the default. Test the cutoff/grid experiment below before adding another controller. |

The claimed 9/9 separator uses the margin from the **default** trajectory. On a new cold start it is unavailable when choosing the cheaper initial fraction. The adverse observations cited are the same Si32 family in three execution configurations. First try leave-one-geometry-family-out validation, and account for the fact that np4/np8 and one core are not independent structures. A margin from the candidate itself can inform a subsequent decision only after its coarse maps have run. Do not execute a second default SCF to obtain the predictor. [Margin definition and sample][follow]

**How much can a stop-only rule buy?** I re-read the 126 published rows at fixed fraction 0.50. For each held-out seed, let an oracle pick the cheapest of all six tested stopping levels, including the baseline 1e−4, with no cost for selecting it. Using the published rounded FFT-column work proxy W, Si32 seed 602 gives 2,823.64 → 2,823.64; seed 603 gives 3,091.66 → 2,841.80. The median paired upside is **4.04%**. The identical retrospective calculation gives **7.41%** for Si16 and **31.64%** for MgO8. These are upper limits for that six-option menu on those one-core Mac fixtures and seeds, in work units. They do not bound changing the fraction, a new first-step rule or dev08 wall. They explain why more work on a universal coarse stop is a poor first bet for Si32. [Underlying six-stop table][rawstops]

## 4. Two inexpensive opportunities in the code

### 4.1 Change the wavefunction basis while keeping the FFT grid fixed

A fraction changes several things together. It changes the wavefunction sphere, density sphere, quantized smooth FFT dimensions, projector/orthogonalization costs and initial fine state. The old minimum-grid experiment fixed the cutoff and looked for a smaller valid **smooth** grid. That branch reproduced existing production dimensions across 42 conditions and closed. The new experiment varies the cutoff *inside a plateau of the existing smooth-grid routine*. It requires no new grid policy and initially no production code. [Basis inclusion][sphere] · [Density sphere and smooth dimensions][grid] · [Closed minimum-grid result][ledger]

I independently enumerated reciprocal Miller vectors for the shipped 30 Ry, Γ-point orthogonal cells. The counts below are full-sphere wavefunctions, prior to real-Γ packing. The density sphere is at four times the wavefunction cutoff. Each dimension is the first 2/3/5-smooth integer ≥ 2 max|Miller index| + 1. These are calculated basis sizes, **not timings**. The 0.35/0.40/0.50 Si32 dimensions match the published independent study. [Si32 fixture][si32] · [Si64 fixture][si64] · [Published grid table][table]

| Coarse fraction | Si32 G vectors | Si32 FFT | Si64 G vectors | Si64 FFT |
|---:|---:|---|---:|---|
| 0.35 | 2,463 | 45×45×24 | 4,945 | 45×45×45 |
| 0.40 | 3,031 | 45×45×24 | 6,043 | 45×45×45 |
| 0.41 | 3,143 | 45×45×24 | 6,355 | 45×45×45 |
| 0.42 | 3,279 | 48×48×24 | 6,619 | 48×48×48 |
| 0.45 | 3,581 | 50×50×25 | 7,153 | 50×50×50 |
| 0.50 | 4,189 | 54×54×25 | 8,409 | 54×54×54 |

**First screen.** Try 0.40 vs 0.41 on Si32 MPI8 at identical native FFT dimensions, and use 0.45 as the known positive control. The 0.41 arm buys 112 more wavefunctions (+3.7%) for the same FFT dimensions. Its eigensolver, density-sphere loops and projector costs can still increase. Test the hypothesis that a better warm state repays that cost. Si64 is a control for harm. For production, derive candidate shell endpoints from the actual cell and existing grid builder, then validate on unseen materials; 0.41 is an experimental point. Equal-size grids do not imply the largest possible basis is optimal: at identical 24×24×45, the one-core Si16 published work proxy was 0.9605× baseline at 0.35 and 1.1025× at 0.40. If cutoff comparisons remain ambiguous, hold a larger *valid* coarse grid fixed at both cutoffs to separate sphere changes from FFT dimensions. Validate candidate self-consistency and changed XC quadrature. [Si16 counterexample][table]

### 4.2 Remove a remaining per-atom × all-plane-waves validity scan

The newly adopted owned projector constructor retains a replicated loop. For every atom, every rank visits **every** momentum, computes a three-term dot product and checks that it and the position are finite. The actual phase computation for owned rows is separate. For the Si64 fixture at fraction 0.50 this is **64 × (23,871 + 8,409) = 2,065,920 pair checks per rank**, assuming one eligible full and coarse constructor each and the recorded 64 Si positions. This is a source-derived count. Instrument invocations and exclusive time to confirm. [Engine call][owncall] · [Validation loop][validation] · [Numeric phase materialization][recipe]

**Suggested exact fast path.** Scan the components of G and atom positions once for finiteness and maximum absolute values. For example, if both maxima are ≤ 2^500, each product is ≤ 2^1000 and a sum of three remains well within binary64 range. Ordinary finite inputs then provably pass the existing finite-dot test. For a failed sufficient test, execute the original loop with its exact error behavior. Keep the factor/norm checks, the actual phase computation, lazy-table guarantees, and collective ordering. Cover NaN, infinity, extremely imbalanced inputs, signed zero, multiple species and large finite components with cancellation. This gate concerns overflow, so it uses numerical range properties and no physical constants. It reduces the typical validation work to O(G + atoms), with identical stored coefficients as the acceptance condition. Give it a default-off switch for the A/B only if the isolated timer is resolvable. Count pair checks, fast accepts, fallbacks and critical-rank wall. Claim at most the measured fraction of this scan; it cannot by itself explain the full QE gap.

## 5. Larger experiment if state quality dominates

**Correct the missing orbitals once instead of estimating when to switch.** F-LADDER-HANDOFF's scalar high-G indicator changed only 0.58–5.08% within a system; its paid probe won on three of seven held-out systems. Keep that policy closed. A distinct test is available because the first genuine fine Hamiltonian application already computes H X for the prolonged orbitals X. With an insulating, fixed-occupation Γ case, project its residual onto the newly introduced plane waves, construct a kinetic-preconditioned complementary block, and test whether that block improves the **starting density** sufficiently to remove later full maps. The relevant idea in the literature is plane-wave perturbative orbital post-processing; here it would serve only as a starting prediction for the unchanged fine SCF. [Failed policy][study] · [Existing residual operator][residual] · [Cancès et al., Journal of Computational Physics][paper]

In an occupied Ritz basis, a diagnostic candidate on the omitted shell is Z_Gi ≈ −R_Gi/(T_G−ε_i). Use cluster-invariant subspaces and a stable projected solve or fallback for a small denominator. Reorthogonalize and compute the trial density from normalized orbitals. Preserve electron number, Γ pairing, frozen occupations and complete occupied clusters. Price the extra reductions, any H update and density formation. Use the existing first-Davidson correction as the same-work control: if it already contains the same span, there is value only if changing the density update saves *additional* certified time. Start with held-out replay. No extra paid fine-H call is free merely because the production solver will later perform a similar call. Stop if the implementation cannot actually share the operator product or beat the equal-work control. The previous broader reduced-subspace predictor remains the fallback research direction if the measured excess is in mid-SCF rather than at handoff. [Prior research proposal and BT3 limits][previous]

## 6. Measurement and implementation order

1. **Scoreboard now.** In one quiet dev08 window, compare uninstrumented current default against the fastest qualifying QE configuration for Si64 MPI8, Si32 MPI8 and then Si32 one core, each separately. Independently lock input hashes, both binaries, identical final basis/grid, bands, PP, physical start, coordinates, thresholds/units, model, binding and complete converged outputs. Revalidate headline candidates on fresh seeds and the actual winner. Do not pool or rank the different parallel models. Record certificate and independent tightened-reference accuracy in the appropriate class.
2. **Fund the crossed seed experiment and the plateau pair first.** One case/model/seed per job, real MPI for wall, certificate enabled, whole SCF when trajectory changes. An earlier ten-step run at 6c31696 stopped unconverged at 7.57 s; the otherwise matched converging input took 10.86 s. Fixed-length screens miss the very maps under study. Baseline/baseline floor first, alternate A/B and B/A, discard the first pair, require an effect ≥3× floor, separated intervals and a moving counter. Record raw per-pair values and within-window spread. Use separate reservations for diagnosis. A planning-only current scoreboard has no numeric winner until these jobs return. [Fast-loop rules][protocol]
3. **Profile the finite-dot validation independently.** A constructor-only test can run in parallel with the ladder work (`dft-physics/src/hamiltonian.rs` vs `dft-engine/src/cutoff_ladder.rs`). The correction and crossed-seed experiments both touch fine initialization; order them sequentially and share the same controls. If the two leading experiments cannot recover sufficient certified seconds, take a current fine solve apart by actual H columns, FFTs per H, bandwidth and collective wait, and prioritize the dominant measured term.
4. **Generality and failures.** Include cubic/skewed/elongated/slab, small and large bands, distinct real PP species and synthetic projectors, ordinary and hard convergence, separate standard 5e−9 / strict 1e−10 Ha references. No successful exit without an actual certificate. Report declines, nonconvergence and losing cases. Retain switch source/effective value, activation counters and full stdout. Adoption is a resolved same-model dev08 gain over the supported domain without physical hard-coding or correctness regression.

**Two harness findings before replay.** I reproduced a file lookup failure for `ladder_handoff.py --cases si64`: the alias table lists si8/si16/si32 and falls through to `benchmarks/cases/si64.toml`. The committed fixture is `si64-gamma.toml`. Use `--cases si64-gamma` or add/test the alias before the new campaign. Also, in the offline `minimal-grid` arm the timer starts before the original coarse problem is constructed, then that arm constructs a second coarse problem with the override. The default constructs one. The reported end-to-end Mac ratios for that arm therefore include asymmetric setup; attributing a ratio such as 1.361 solely to FFT factorization is unsupported. This correction leaves intact the production finding that the valid smooth minimum was already selected. [Fixture loader][loader] · [Timed constructor sequence][double] · [Clock end][driverseed]

**Requests for the dev08 owner.** Use the existing start/done/failed/idle signals at `/data/sig-dft/signals`. Keep each job within 600 s except `final-` comparisons. Each arm below gets a fresh output and restart directory, with no undeclared BAANDI overrides. Use absolute executable and input paths from the matched production manifest. Verify the effective fraction in the result; an input-specified ladder configuration must agree with the requested arm. For the standard controls, leave the adopted defaults on.

```bash
# One arm of plateau-si32-np8-601. FRACTION is 0.40 or 0.41.
# The owner wraps each entire launch in the same high-resolution external clock.
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    BAANDI_CUTOFF_LADDER="$FRACTION" \
    mpirun -np 8 --map-by core --bind-to core --report-bindings \
    "$BIN" run "$INPUT" --output-dir "$OUT/calc" >"$OUT/stdout" 2>"$OUT/stderr"
# final-si64-np8-refresh: BAANDI uses the same invocation with the ladder override unset.
# QE_FLAGS is the bash array of the winning MPI settings, locked in the manifest.
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    mpirun -np 8 --map-by core --bind-to core --report-bindings \
    "$QE_BIN" "${QE_FLAGS[@]}" -in "$QE_INPUT" >"$OUT/qe.stdout" 2>"$OUT/qe.stderr"
```

**Requested return.** For each system and model send a row with current B, Q and deficit in seconds; coarse/ordinary fine/verification map counts and H work; per-phase and external wall; actual first and post-redo tolerance; certificate and tightened-reference error. For the new interventions also send all four crossed-seed results and the same-FFT pair. This determines whether the next engineer should repair the quality of the first density, reduce the cost of a fine Hamiltonian application, or pursue the later subspace predictor.

[head]: https://github.com/chaeeundad/baandi/tree/3effafef4fe9d722ae75ac34d03abf0a170e37b0
[delta]: https://github.com/chaeeundad/baandi/compare/63aa4adbaa792a429762ff8c71fda32ba73a1e5d...3effafef4fe9d722ae75ac34d03abf0a170e37b0
[previous]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/docs/performance/plan-2026-11-19.md
[rules]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/docs/performance/comparison-rules.md
[ledger]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/docs/performance/attempts-ledger.md
[defaults]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-numerics/src/switches.rs#L307
[sched]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/docs/performance/attempts-ledger.md#L2142
[inventory]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/docs/performance/test-targets.md
[retry]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/session.rs#L1043
[gram]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/tickets/G-GUARD-GRAM.md
[sessions]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/artifacts/sessions.json
[mixer]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/mix_dist.rs#L195
[study]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/docs/design/ladder-handoff.md#L157
[fraction]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/tickets/F-LADDER-FRACTION.md
[follow]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/tickets/F-LADDER-FRACTION-msg1.md
[launch]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/tickets/F-LAUNCH-IO.md
[protocol]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/docs/performance/fast-loop-protocol.md#L154
[coarse]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/cutoff_ladder.rs#L1187
[first]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/cutoff_ladder.rs#L1423
[ratchet]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/convergence.rs#L23
[redo]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/session.rs#L763
[driverseed]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/cutoff_ladder.rs#L2035
[seedsetter]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/session.rs#L1934
[sphere]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-basis/src/lib.rs#L102
[grid]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-basis/src/lib.rs#L237
[si32]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/benchmarks/cases/si32-gamma.toml
[si64]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/benchmarks/cases/si64-gamma.toml
[table]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/docs/design/ladder-handoff.md#L228
[owncall]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/scf.rs#L1568
[validation]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-physics/src/hamiltonian.rs#L2356
[recipe]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-physics/src/hamiltonian.rs#L342
[residual]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/cutoff_ladder.rs#L1678
[paper]: https://doi.org/10.1016/j.jcp.2015.12.012
[loader]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/scripts/ladder_handoff.py#L33
[double]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/crates/dft-engine/src/cutoff_ladder.rs#L1990
[rawstops]: https://github.com/chaeeundad/baandi/blob/3effafef4fe9d722ae75ac34d03abf0a170e37b0/docs/design/ladder-handoff.md#L403
