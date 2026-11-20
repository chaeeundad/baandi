# First-attempt scheduling experiment (reopened)

Final outcome: accuracy gates pass; the intended Si64 attempt-count reduction
is not demonstrated. Defaults stay unchanged. See [final evidence](first-schedule-validation.md).

The orchestrator explicitly reopened attempt timing; the ADR 0043 acceptance
rule remains unchanged. `BAANDI_VERIFY_SCHEDULE=first` is opt-in. Default `on`
retains the adopted retry scheduler, including its arithmetic and decisions.
`first-trace` predicts the first attempt without deferring it, then retains
`on` retry scheduling. This differs intentionally from legacy `trace`, which
also declines to act on retry forecasts. No switch is added; registry count
is unchanged.

## Corrected diagnosis and candidates

The original scheduler never predicts the first attempt: it returns immediately
until an actual failure exists. The relayed np8 `prediction_evaluations=1`
therefore describes the retry, not a first-attempt false pass prediction.

Relayed dev08 `bad02ee`, six runs per model, Si64:

| Model | Median external wall (s) | Verification sum (s) | First failure | Final pass |
|---|---:|---:|---|---|
| MPI8 | 8.5466 | 1.0475 (12.26%) | iteration 10, 1.009e-8 Ha, 380 H | iteration 11, 2.822e-9 Ha, 388 H |
| One core | 41.5704 | 4.8284 (11.62%) | iteration 11, 5.310e-8 Ha, 395 H | iteration 13, 1.172e-9 Ha, 356 H |

One-core phases: solve 3.9565 s, density 0.6604 s, energy 0.1140 s,
potential 0.0889 s. Solve apply 1.6961 s, other 0.7427 s (alias
other_engine, count once), rotate 0.6582 s, orthonormalize 0.3240 s,
reduced_gram 0.2815 s. Independently summarized phases need not sum to the
reported median verification total. These are relayed measurements, not Mac
profiles or measured candidate savings.

1. **Selected: bounded first-attempt geometric-tail envelope.** If it removes
   the initial failure with no other trajectory cost, the opportunity is about
   0.54 s on MPI8 and 2.52 s on one core, separately. Skipped attempts are never
   passes. Unlike the earlier closing scope, timing is now explicitly allowed.
2. **Alternative: fit a universal first-attempt amplification from historical
   verification/stop-error ratios.** Those ratios vary strongly even for equal
   thresholds; a fitted global multiplier could exploit the data but needs
   independent geometry/solver holdouts and would introduce a calibration
   constant. Not selected for this experiment.

Failed HALF verification feeds the mixer and supplies the next warm start.
Skipping it changes that continuation. Thus the ordinary maps on the baseline
are not a proof that the candidate performs the same maps/counts. The Mac
study measures the actual candidate as well as counterfactual first decisions.
The optimistic opportunity is not a forecast of net gain. In a one-extra-map
binary cost model the required correct-skip fraction is Cmap/(Cfail+Cmap):
about 0.491 for MPI8 and 0.430 for one core from the relayed costs. Changed
later attempts and two-map delays invalidate that simplified cost model;
whole-run paired dev08 measurements decide.

## Model, derivation and limits

Let E_j be the ordinary Coulomb stop estimate (Ha), d_j its density L2 residual,
and T the existing verification threshold. At the first crossing j of E<T,
use a geometric residual-amplitude tail anchored at the preceding map, with

    q = max(E_j / E_(j-1), (d_j / d_(j-1))^2).
    L = E_(j-1) / (1 - sqrt(q))^2.
    predicted verification error after k extra maps = L q^k.

The squared L2 ratio and Coulomb-energy ratio have the same residual-squared
scaling. Using the slower rate avoids letting a rapid change in only one norm
drive the forecast. The squared tail follows from summing residual amplitudes:
`sqrt(E_(j-1)) * sum_n sqrt(q)^n`, then squaring. This models the
remaining sequence of corrections rather than only its first term. The
lagged anchor and its use as an unmixed-response surrogate are explicit
conservative heuristics for
a Pulay-cancelled stop residual: accelerated mixed-map contraction does not
measure the response of the unmixed verification map. **It is not a derived
upper bound on F, a probability, or a calibrated guarantee of failure.** In
particular a valid threshold crossing can still pass immediately; those wrong
skips must be counted, not explained away. At a fresh eligible crossing this
model normally delays at least one map; it is a bounded persistence policy,
not a classifier proven to distinguish all successful attempts.

Choose the first k with L q^k<T, capped at two ordinary maps since crossing.
Require finite positive quantities, 0<q<1, and a preceding mixed step with at
least two history entries. Decline after a mixer reset/raw recovery, without
sufficient history, during spectral audit, or with invalid/noncontracting
ratios. Attempt unconditionally on the last iteration budget. Subsequent
eligible calls reconstruct the same crossing and rate from scientific history;
the journal is never control state. No environment reads, array copies,
collectives or clocks are added by this policy. At most a history scan and a
few scalars execute at an otherwise eligible verification boundary.

Threshold, eigen tolerance, electron count and volume are retained in study
features; the latter two are deliberately not model parameters. E/T is already
dimensionless and the norm ratios cancel volume normalization. No atom count,
material identity or silicon constant enters the rule. The rule does not
loosen the inner eigensolve, change gap/residual guards or change delivered
state semantics. Every delivered result still uses the existing verification.

## Evidence protocol

`first_schedule_study.py` reads the existing F-VERIFY-SCHEDULE records for eight
systems (Si8/16/32/skew, MgO8, synthetic Gaussian FCC, p3, spd), seeds 601–603,
thresholds 5e-9 and 1e-10 Ha: 48 first-decision labels. It saves the stop error,
verification error/amplification, both contraction ratios, density residual,
eigen tolerance, electron count and cell volume. Historical labels use the
adopted `on` arm, not the old `off` retry policy.

Fresh parent/default/first-trace/first executions use identical stored inputs.
Trace provides the directly observed counterfactual label at the same initial
state; actual candidate counts quantify the consequences of omitting HALF
products. "Skipped failure" means a predicted wait at a baseline failure,
not a guarantee of one fewer total attempt. "Wrong skip" means a predicted
wait at a baseline pass. Every arm must converge and verify. Strict-reference
parent/first identity is checked separately for each geometry. Both accuracy
classes are checked against that strict energy, not just against each other.
Mac clocks are disabled; no Mac timing supports a performance conclusion.

Reproduction and final tables are in `first-schedule-validation.md`.

## Historical first-attempt feature study

Existing Mac records, correctness/counter evidence only; six observations per row (three seeds, two thresholds). No coefficients were fitted to these labels.

| System | First failures / 6 | Proposed failure skips | Wrong skips | Verification / stop error range | E/T range | Electrons | Volume (Bohr^3) |
|---|---:|---:|---:|---|---|---:|---:|
| si8 | 6 | 6 | 0 | 5.83–8.83 | 0.357–0.911 | 32 | 1080 |
| si16 | 6 | 6 | 0 | 19.3–30.5 | 0.0599–0.934 | 64 | 2161 |
| si32 | 4 | 4 | 2 | 16.2–23 | 0.0214–0.558 | 128 | 4322 |
| si-skew | 6 | 6 | 0 | 12.4–32.4 | 0.12–0.79 | 8 | 243.8 |
| mgo8 | 1 | 1 | 5 | 0.343–2.99 | 0.0184–0.832 | 64 | 503.6 |
| syn-gaussian-fcc | 6 | 5 | 0 | 16.9–51.4 | 0.0504–0.952 | 8 | 270.1 |
| syn-p3 | 0 | 0 | 6 | 0.0362–0.0372 | 0.0048–0.0246 | 4 | 64 |
| syn-spd | 0 | 0 | 3 | 0.00637–0.00665 | 0.0225–0.841 | 4 | 64 |

At seed 601 / 5e-9 Ha: Si16 has E/T=0.103 yet verification/stop=30.5 and fails; Si32 has E/T=0.0214 with amplification=16.2 and passes. The skew cell has E/T=0.587 and amplification=22.2. MgO has E/T=0.0184 and amplification=0.463. The synthetic spd case has E/T=0.833 but amplification=0.00664 and passes. Thus neither size nor proximity to the threshold is itself a failure certificate.

The saved rows also contain the requested eigen tolerance and both residual ratios. These do not bound the unmixed response: a tolerance is a solver stopping parameter, not a density-error theorem, and the geometric mixed-map rate can coexist with amplification greater than one. The initial lag-only hypothesis deliberately accepted false positives; its 28 proposed failure skips and 16 wrong skips (one missed failure) made a dev08 control essential. The table above records that initial screening; final tail-model results are in the validation document. These are baseline-label counts, not measured reductions in total attempts.

## Rejected lag-only variant and revision

The first implementation used L=E_(j-1), without the geometric sum. Across
48 runs it reduced attempts 81->58 and total H 48,676->47,564, with all accuracy
gates passing. However the additional Si64 one-core correctness run moved the
first failure from map 11 to map 12 (395->394 H); two attempts and 13 maps
remained, with total H 6,644->6,622. This did not support removing the relayed
2.52 s failure cost. Raw work/accuracy evidence is retained in
`first-schedule-lag-rejected.json`; this is not a Mac performance conclusion.

The selected revision includes the geometric-amplitude tail above. It adds no
fitted coefficient, material identifier or new acceptance condition. Its
whole 8-system/3-seed matrix is rerun, including known wrong-skip controls.
The same dataset informed development, so these are regression/generalization
checks, not a statistically independent held-out calibration. dev08 must still
measure the actual cost and withdrawals; the model is not a guarantee.

Eigen-tolerance example from the same seed/standard records: N_e*tau/T is
about 0.505 for failing Si16, 1.67 for passing MgO, and 5.64e5 for passing
synthetic spd (whose tiny subspace solves essentially exactly). Thus the
requested eigentolerance is not the actual eigen residual and cannot alone
predict a verification failure. Raw density residuals also depend on cell
normalization; the selected model uses only their dimensionless contraction
ratio and the Coulomb threshold already defined by the contract.
