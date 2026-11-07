# P37A generalized Gram stage 0

## Prespecified decision rule (before measurement)

Recommend a generalized module only if **all** of the following hold over all
analysed solves, cases and capacities (2, 3 and 6 times the band count):

1. Canonical-orthogonalization fallback is needed in **under 5%** of steps at
   the relative overlap eigenvalue threshold **1e-8**.
2. In the strict profile the smallest overlap eigenvalue **never falls below
   1e-12 relative** to the largest overlap eigenvalue.
3. Counted algebra multiply-adds per step are **at least 40% below** the current
   real path on the same solves.

Otherwise identify each failed condition and case, especially the A52 twenty-atom
synthetic skew fixture. Missing coverage is not a pass. These conditions will
not be relaxed after measurement.

This is measurement only. SCF uses its existing solver; diagnostic trajectories
must never feed vectors or densities back into SCF. No production defaults,
tolerances, source files, artifacts or pseudopotential files are changed.

## Premise and provenance

Derived from the own-words notes `qe-stage-by-stage.md` and
`qe-solvers-and-subspace.md`; no QE source is read or copied. Q2 reopens B13/A16
because the relayed real strict D0 algebra share is 48.3% at T=1 and 57.1% at
T=8. Standard-profile D0 remains unmeasured. A12 contraction is retained;
A13 capacity-only, A52 CholeskyQR2 and P24A hard locking are not being proposed
again. The generalized small problem changes the rank-decision location;
that premise needs measurement, especially on A52's synthetic skew geometry.
P32A/P33A timing results do not establish generalized-basis safety.

## Status

**Do not recommend the generalized module under the prespecified rule.**
Both numerical safety conditions fail on strict Si8 (all capacities) and strict
Si16 (capacities 2 and 3). Several measured work comparisons also miss 40%.
Standard-only results do not override a rule that explicitly includes strict.

The final dataset contains 28 SCF trajectories on nine fixtures, 2,410 original
Davidson solves and 67,257 diagnostic steps across three capacities. Si8 has
five perturbed starts per profile; A52 has uniform and seed-1 starts. Twenty
SCFs converge; eight retain the existing elongated/slab failure or the A52
six-step budget exhaustion. Every solve encountered on those censored histories
is included. Forty-five diagnostic capacity/solve trajectories reach the
iteration limit; no failed diagnostic is treated as a successful generalized
solve. The original SCF continues using its original solver.

## Measurement implementation

The ignored engine test writes an isolated source copy under
`P37A-evidence/source`, then injects the solver test's diagnostic payload there.
Both builds share the worktree's one `target/`; the original production files
are never edited. The copy calls the diagnostic after the original solver has
prepared its initial basis. Each diagnostic borrows the exact frozen operator,
initial basis, band policy, kinetic preconditioner and requested tolerance of
that SCF solve, including certificate solves. Its result is discarded before
the original solver continues. The existing P31A fixture/perturbation driver
supplies SCF histories. No UPF is copied; A52's synthetic UPF is generated and
parsed in memory using the existing generator.

Gamma uses real ORC coordinates and real products/eigendecompositions; k points
use complex Hermitian products. Cholesky is faer's library factorization.
Above relative 1e-8 a Cholesky whitening transform is used. Below that threshold
(or on factorization failure), canonical whitening retains eigenvalues above
relative 1e-12. Pre-fallback spectra are recorded, including negative roundoff
eigenvalues; rank uses relative 1e-10. Too few retained directions terminates
that diagnostic, with an explicit rank-exhaustion record. Iteration limits and
nonconverged SCF histories remain visible, not successful safety evidence.
The 1e-6 fallback rate is a counterfactual diagnostic on the same trajectories,
not a separate continuation policy. Every correction is normalized without any
projection against the basis or other corrections. Restarts contract V and HV,
and contract the small H and S with the same coefficients. Only the new
rectangle is charged on the next step. An independent full Gram rebuild checks
the actual floating-point basis at every step, records drift from the cached S,
and supplies the measured overlap to the diagnostic solve. This extra
measurement work is excluded from the candidate work estimate; it is not an
optimized production implementation. No exact-identity Gram reset is assumed.

The strict diagnostic evaluates all root residuals to decide convergence.
Only standard can select roots from small-space Ritz changes before forming
long vectors. Neither profile freezes converged roots. No diagnostic result is
a strict SCF certificate. Analytic three-dimensional real and complex oracles
check the whitening spectrum and residuals to 1e-12 Ha; a duplicate-column
overlap checks Cholesky failure.

### Work units and limits

One real MAC is one multiply-add (two real FLOPs); one complex MAC corresponds
to four real MACs (eight real FLOPs). Baseline real counters wrap every actual
solver GEMM call, including projection, incremental H, rotations and restarts,
and count each executed MGS dot/update pair. Initial basis preparation is
outside both compared trajectories. Hamiltonian/FFT work, norm/division costs,
Gram-panel assembly additions and diagnostic spectrum collection are excluded.
Current real H already uses one rectangle, so there is no fictitious second-H
rectangle saving.

Generalized rectangle and rotation counts use the actual dimensions visited,
including both contraction products. Small H/S contractions are counted in the
small-work column. The independent full Gram reconstruction is diagnostic work.
Small EVD work is separately labelled a model: 2 m^3 MAC-equivalents per full
EVD, plus explicitly shaped whitening products. It is **not** an instrumented
faer EVD operation count and not a range-selected generalized-driver timing.
Actual-call MAC comparisons and model-inclusive comparisons must remain
separate. A modeled saving alone cannot certify the exact 40% acceptance gate.
There is no wall-time speed claim or dev08 performance acceptance.

## Reproduction and coverage

Baseline HEAD: `6c34eed6125c3e05338ddef486edcf1a8bb30c9b`. Ordinary release,
local Mac, pure Rust, T=1. Stored evidence: `P37A-evidence/` (outside target),
with raw logs, command/exit-code JSON, input UPF hashes, generated source copy,
per-case CSV/JSON and table generator. `initial-full-rebuild/` preserves the
superseded pilot that rebuilt retained H/S blocks after restart; its numbers
are excluded from final tables. Fixing that cost-accounting issue did not change
the prespecified decision thresholds.

```sh
P37A_CASES=si8-gamma,si16-gamma,si-skew,si-elongated,si-slab,syn20-skew,si-k222s,material-si-20,material-mgo-40 \
P37A_PROFILE=strict P37A_START_RUN=2 P37A_END_RUN=7 \
P37A_SI_UPF=/absolute/Si.pz-vbc.UPF DFT_NC_DIR=/absolute/pseudodojo/ex \
cargo test --offline --release -p dft-engine --test generalized_gram_stage0 \
histories -- --ignored --nocapture --test-threads=1
```

Repeat with `P37A_PROFILE=standard`. The five runs (2..7) mean seeds 1..5
only for Si8; other cases run seed 1 once per profile. The seed generator and
relative density perturbation amplitude 1e-10 are the existing P29A/P31A
protocol. Additionally run `P37A_CASES=syn20-skew P37A_UNPERTURBED=1` for each
profile to include the original uniform start. After the first wrapper run,
`python3 -I P37A-evidence/run.py` can repeat batches with the same environment
without rebuilding the wrapper. `python3 -I P37A-evidence/analyze.py` regenerates
tables. `P37A_CONTROL=1` runs the identical fixture driver with original,
uninstrumented solver sources for scientific-history identity checks.

The A52 fixture is not `si-skew` (which has two atoms), and has no checked-in
`syn20-skew.toml`. Its source is
`crates/dft-numerics/src/pool/orthonormalize_synthetic.py`: 20 atoms, lattice
rows [14,0,0], [1.5,14,0], [0.5,0.7,14] Bohr, 20 Ry, 24 bands, three synthetic
projectors with angular momenta 0/0/1. The six-iteration budget matches A52's
bounded experiment and is not a converged material result. The real path and
current source are intentionally different from the historical complex
CholeskyQR2 attempt; this study cannot retrospectively identify its exact
233-to-329 iteration mechanism.

Other fixtures use their checked-in geometry, cutoff and k mesh. PseudoDojo
Si and MgO use `tests/reference/materials/si-20.toml` and `mgo-40.toml`, PBE,
density floor 1e-10, their Kerker-Pulay settings and 120-iteration budget.
These and `si-k222s` use complex coordinates; six Gamma geometries use ORC.
Strict and standard are separate throughout. No QE reference run was performed.

UPF SHA-256 (external files only):

- Si.pz-vbc.UPF: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`
- PseudoDojo Si: `39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282`
- Mg: `15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14`
- O: `2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11`

## Decision by condition

1. **Fallback below 5%: fails.** Strict Si8 needs fallback in
   64.60% / 64.57% / 34.83% of steps at capacities 2 / 3 / 6. Strict Si16 needs
   41.39% / 36.98% at 2 / 3. The other case/profile/capacity aggregates are
   below 5%, including A52; this does not rescue the all-case condition.
2. **Strict minimum relative eigenvalue at least 1e-12: fails.** Si8 minima are
   -2.429e-16 / -1.175e-15 / -1.388e-16; Si16 minima at 2 / 3 are
   3.537e-15 / 4.074e-15. Negative measured Gram eigenvalues are roundoff at
   numerical rank loss, not a claim that an exact Gram matrix is indefinite.
3. **At least 40% counted algebra saving: not established generally.** Even the
   counted large-product/MGS comparison misses 40% for strict capacity 6 on
   Si8 (18.80%) and Si16 (18.62%), every si-skew and si-slab row, every strict
   si-elongated row, standard si-elongated capacity 6, and strict A52 capacities
   2 and 6. Small-solve model columns are separate; an exact faer EVD operation
   count was not obtained and is not silently substituted by a cubic estimate.
   The table also reports the worst per-solve saving rather than accepting an
   average as a universal bound. Conditions 1 and 2 independently reject the
   recommendation, irrespective of this small-solve accounting limitation.

**A52 specifically:** strict 2 / 3 / 6 has fallback 0.00% / 1.08% / 2.03%,
minimum relative eigenvalue 3.727e-6 / 5.768e-11 / 4.760e-10, and counted
MAC/step saving 20.97% / 44.99% / 18.99%. Thus its observed safety conditions
pass, but capacities 2 and 6 miss the counted-product work target. At capacity
3, 8 of 12 diagnostics still exhaust their 400-step budget; total counted MACs
are 3.59 times the original solves despite the lower mean cost per step.
Capacity 2 exhausts all 12 diagnostics and costs 5.70 times the original total.
These are current real-path histories on the A52 geometry, not a replay of
historical CholeskyQR2's complex trajectory. No causal claim about A52's old
rank test versus rounding amplification is made.

All standard Gram aggregates pass the observed fallback test
(the 1e-12 minimum rule is strict-only). Standard Si8 MAC/step savings are
45.54% / 50.28% / 50.32%; that restricted result is not general adoption
evidence. The largest observed difference between cached and independently
rebuilt overlap entries is 1.584e-4, dimensionless. This diagnostic deliberately
uses the rebuilt overlap in the small solve; a fully cached production module
would need its own drift analysis and cannot inherit a safety guarantee here.

## Verification and remaining scope

- Original-source controls match all 1,176 recorded scientific lines across
  four paired batches: every SCF history row, final energy/certification flag,
  certificate summary, terminal status and reported gap. Diagnostic results
  never replace SCF vectors. This is exact record identity, not a claim that
  orbital arrays were separately persisted and compared.
- Five standard benchmark cases (gamma, k222s, k222, k444s, si8), each at
  T=1,2,4,8,12, pass `bench_scf.py --require-identical`: all 25 comparisons
  have energy difference **0 Ha** and identical density SHA-256 against the
  same-HEAD T=1 production run. Single repetitions are equivalence checks,
  not timing measurements. Source hashes show no production change.
- The standalone real/complex analytic spectrum and residual oracle passes
  1e-12 Ha, and duplicate-column Cholesky failure is detected in every batch.
- No production accessor, manifest, default, tolerance, README result block,
  artifact JSON or UPF was changed. All changes remain uncommitted.
- Exact small-EVD MAC accounting, cached-only Gram trajectories and dev08
  performance measurements remain outside these results. They are not reasons
  to waive the failed safety criteria. No generalized module is recommended.

The literal `scripts/check.sh all` preflight exits 1 because this repository
interprets `all` as a package name. Its documented workspace form is
`scripts/check.sh` with no argument. Two intermediate full checks found the
same constant-assertion clippy warning in the ignored runner; it was replaced
by an explicit debug-mode guard. Those logs are retained. Final gate results
are recorded below, not inferred from the failed invocations.

## Measured overlap statistics

Eigenvalues and relative ratios are dimensionless. Percentiles use linear interpolation across steps; all five Si8 seeds are pooled; A52 pools the uniform start and seed 1 when both are present. Fallback includes Cholesky failure or a relative eigenvalue below the threshold.

| Case | Profile | Capacity | Solves / steps | min lambda | p1 | p10 | p50 | min relative | Cholesky fail % | fallback 1e-6 % | fallback 1e-8 % | max rank deficit | unfinished solves |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| material-mgo-40 | standard | 2 | 182 / 667 | 1.056e-03 | 2.190e-03 | 1.365e-02 | 1.507e-01 | 2.256e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| material-mgo-40 | standard | 3 | 182 / 623 | 6.133e-04 | 1.290e-03 | 6.885e-03 | 3.432e-02 | 1.607e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| material-mgo-40 | standard | 6 | 182 / 609 | 6.038e-04 | 6.201e-04 | 4.672e-03 | 2.959e-02 | 1.418e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| material-mgo-40 | strict | 2 | 280 / 1965 | 8.445e-09 | 2.512e-06 | 1.079e-03 | 3.468e-02 | 1.449e-09 | 0.00 | 1.12 | 0.31 | 0 | 0 |
| material-mgo-40 | strict | 3 | 280 / 1491 | 1.363e-06 | 3.535e-04 | 4.042e-03 | 4.071e-02 | 2.084e-07 | 0.00 | 0.40 | 0.00 | 0 | 0 |
| material-mgo-40 | strict | 6 | 280 / 1306 | 4.578e-04 | 6.916e-04 | 3.303e-03 | 2.504e-02 | 9.349e-05 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| material-si-20 | standard | 2 | 352 / 1013 | 3.228e-03 | 4.421e-03 | 1.983e-02 | 1.255e-01 | 5.754e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| material-si-20 | standard | 3 | 352 / 1031 | 2.010e-03 | 3.022e-03 | 1.208e-02 | 7.764e-02 | 3.620e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| material-si-20 | standard | 6 | 352 / 1042 | 2.010e-03 | 3.030e-03 | 9.802e-03 | 6.675e-02 | 3.620e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| material-si-20 | strict | 2 | 512 / 4671 | 8.209e-07 | 3.481e-05 | 1.547e-03 | 2.155e-02 | 1.666e-07 | 0.00 | 0.39 | 0.00 | 0 | 0 |
| material-si-20 | strict | 3 | 512 / 3661 | 1.812e-07 | 3.474e-05 | 3.689e-03 | 3.286e-02 | 1.488e-08 | 0.00 | 0.74 | 0.00 | 0 | 0 |
| material-si-20 | strict | 6 | 512 / 2865 | 2.237e-03 | 2.543e-03 | 4.386e-03 | 2.218e-02 | 3.302e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-elongated | standard | 2 | 105 / 305 | 6.517e-04 | 3.517e-03 | 3.795e-02 | 2.094e-01 | 1.325e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-elongated | standard | 3 | 105 / 297 | 3.308e-03 | 1.190e-02 | 4.102e-02 | 1.435e-01 | 1.043e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-elongated | standard | 6 | 105 / 292 | 3.308e-03 | 1.445e-02 | 3.457e-02 | 1.410e-01 | 1.043e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-elongated | strict | 2 | 200 / 1132 | 2.433e-09 | 9.705e-07 | 5.531e-04 | 9.824e-02 | 2.942e-10 | 0.00 | 1.41 | 0.62 | 0 | 0 |
| si-elongated | strict | 3 | 200 / 884 | 1.541e-08 | 4.029e-08 | 1.003e-03 | 8.257e-02 | 3.492e-09 | 0.00 | 3.28 | 1.58 | 0 | 0 |
| si-elongated | strict | 6 | 200 / 780 | 2.189e-08 | 5.579e-06 | 4.424e-03 | 8.557e-02 | 1.877e-09 | 0.00 | 1.54 | 0.38 | 0 | 0 |
| si-k222s | standard | 2 | 36 / 104 | 3.021e-02 | 3.021e-02 | 4.965e-02 | 2.666e-01 | 1.132e-02 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-k222s | standard | 3 | 36 / 106 | 3.021e-02 | 3.021e-02 | 4.032e-02 | 1.282e-01 | 9.912e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-k222s | standard | 6 | 36 / 106 | 3.021e-02 | 3.021e-02 | 3.853e-02 | 1.282e-01 | 7.254e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-k222s | strict | 2 | 52 / 374 | 2.150e-05 | 5.016e-05 | 2.650e-03 | 4.667e-02 | 7.199e-06 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-k222s | strict | 3 | 52 / 300 | 3.404e-04 | 3.410e-04 | 6.761e-03 | 6.740e-02 | 6.346e-05 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-k222s | strict | 6 | 52 / 261 | 2.996e-03 | 3.858e-03 | 1.253e-02 | 6.637e-02 | 6.067e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-skew | standard | 2 | 14 / 33 | 4.828e-03 | 9.364e-03 | 2.586e-02 | 3.563e-01 | 1.349e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-skew | standard | 3 | 14 / 34 | 4.828e-03 | 9.506e-03 | 2.469e-02 | 2.766e-01 | 1.349e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-skew | standard | 6 | 14 / 34 | 4.828e-03 | 9.506e-03 | 2.469e-02 | 2.702e-01 | 1.349e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-skew | strict | 2 | 25 / 121 | 4.156e-03 | 5.943e-03 | 3.445e-02 | 1.709e-01 | 1.219e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-skew | strict | 3 | 25 / 107 | 2.819e-02 | 2.900e-02 | 5.624e-02 | 1.629e-01 | 7.677e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-skew | strict | 6 | 25 / 101 | 1.346e-02 | 2.828e-02 | 4.953e-02 | 1.629e-01 | 3.057e-03 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-slab | standard | 2 | 200 / 546 | 1.051e-03 | 3.981e-03 | 2.525e-02 | 1.862e-01 | 3.991e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-slab | standard | 3 | 200 / 524 | 1.211e-03 | 2.829e-03 | 1.631e-02 | 1.209e-01 | 2.935e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-slab | standard | 6 | 200 / 527 | 5.853e-04 | 1.349e-03 | 1.021e-02 | 1.003e-01 | 1.560e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si-slab | strict | 2 | 200 / 1363 | 3.807e-07 | 3.328e-06 | 2.819e-04 | 4.112e-02 | 3.380e-08 | 0.00 | 1.39 | 0.00 | 0 | 0 |
| si-slab | strict | 3 | 200 / 1048 | 2.834e-08 | 8.899e-06 | 4.921e-04 | 3.506e-02 | 3.748e-09 | 0.00 | 0.76 | 0.19 | 0 | 0 |
| si-slab | strict | 6 | 200 / 880 | 2.834e-08 | 6.474e-06 | 1.379e-03 | 2.661e-02 | 3.748e-09 | 0.00 | 1.14 | 0.23 | 0 | 0 |
| si16-gamma | standard | 2 | 15 / 62 | 6.617e-05 | 3.321e-04 | 1.102e-03 | 4.242e-02 | 1.052e-05 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si16-gamma | standard | 3 | 15 / 65 | 8.183e-05 | 9.825e-05 | 7.145e-04 | 3.643e-02 | 8.505e-06 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si16-gamma | standard | 6 | 15 / 77 | 1.198e-05 | 1.213e-05 | 7.303e-05 | 2.128e-02 | 8.330e-07 | 0.00 | 3.90 | 0.00 | 0 | 0 |
| si16-gamma | strict | 2 | 13 / 894 | 9.689e-14 | 7.424e-13 | 1.220e-10 | 5.418e-07 | 3.537e-15 | 0.00 | 62.98 | 41.39 | 9 | 0 |
| si16-gamma | strict | 3 | 13 / 530 | 1.409e-13 | 1.762e-13 | 4.578e-11 | 3.741e-06 | 4.074e-15 | 0.00 | 58.11 | 36.98 | 14 | 0 |
| si16-gamma | strict | 6 | 13 / 194 | 9.614e-07 | 9.623e-07 | 2.569e-06 | 1.821e-03 | 5.173e-08 | 0.00 | 14.95 | 0.00 | 0 | 0 |
| si8-gamma | standard | 2 | 125 / 463 | 4.270e-06 | 1.922e-05 | 5.815e-04 | 3.217e-02 | 3.977e-07 | 0.00 | 0.65 | 0.00 | 0 | 0 |
| si8-gamma | standard | 3 | 125 / 511 | 1.471e-05 | 1.914e-05 | 4.607e-04 | 2.105e-02 | 1.529e-06 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si8-gamma | standard | 6 | 125 / 515 | 1.471e-05 | 1.599e-05 | 4.027e-04 | 1.865e-02 | 1.267e-06 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| si8-gamma | strict | 2 | 75 / 15002 | -3.381e-15 | -2.053e-17 | 1.363e-13 | 3.055e-09 | -2.429e-16 | 1.00 | 79.75 | 64.60 | 13 | 23 |
| si8-gamma | strict | 3 | 75 / 5970 | -3.727e-15 | -1.964e-15 | -5.968e-17 | 7.779e-10 | -1.175e-15 | 10.23 | 77.89 | 64.57 | 40 | 2 |
| si8-gamma | strict | 6 | 75 / 1777 | -2.571e-15 | -2.192e-15 | 1.033e-11 | 2.890e-06 | -1.388e-16 | 5.51 | 52.84 | 34.83 | 26 | 0 |
| syn20-skew | standard | 2 | 12 / 88 | 9.845e-04 | 1.014e-03 | 1.107e-02 | 9.936e-02 | 3.347e-04 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| syn20-skew | standard | 3 | 12 / 80 | 2.317e-04 | 2.317e-04 | 6.783e-03 | 5.252e-02 | 5.512e-05 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| syn20-skew | standard | 6 | 12 / 94 | 3.695e-05 | 6.517e-05 | 4.145e-03 | 1.866e-02 | 6.755e-06 | 0.00 | 0.00 | 0.00 | 0 | 0 |
| syn20-skew | strict | 2 | 12 / 4800 | 4.627e-05 | 3.241e-04 | 9.825e-03 | 2.758e-02 | 3.727e-06 | 0.00 | 0.00 | 0.00 | 0 | 12 |
| syn20-skew | strict | 3 | 12 / 4342 | 9.440e-10 | 8.895e-08 | 3.162e-05 | 7.540e-04 | 5.768e-11 | 0.00 | 4.88 | 1.08 | 1 | 8 |
| syn20-skew | strict | 6 | 12 / 590 | 1.833e-08 | 3.339e-08 | 1.358e-06 | 2.214e-04 | 4.760e-10 | 0.00 | 26.10 | 2.03 | 0 | 0 |

## Measured algebra dimensions

MAC/step columns are millions of real MACs for Gamma; complex cases are millions of complex MACs and have no real baseline. The small-solve columns are explicitly modeled MAC-equivalents. Total MAC ratio includes the change in diagnostic step count. These are work counts, not timings.

| Case | Profile | Capacity | H+S rectangle | Packed Ritz | Contraction | Small model | Current real MAC/step | Current small model | MAC saving % | Model-inclusive saving % | Worst solve MAC saving % | Total MAC ratio |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| material-mgo-40 | standard | 2 | 0.5241 | 0.3096 | 0.5086 | 0.0581 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-mgo-40 | standard | 3 | 0.6091 | 0.3304 | 0.3763 | 0.0976 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-mgo-40 | standard | 6 | 0.6984 | 0.3722 | 0.3523 | 0.2201 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-mgo-40 | strict | 2 | 0.4174 | 0.7055 | 0.4282 | 0.0621 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-mgo-40 | strict | 3 | 0.5623 | 0.8532 | 0.3616 | 0.1184 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-mgo-40 | strict | 6 | 0.7103 | 1.0570 | 0.3252 | 0.2477 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-si-20 | standard | 2 | 0.1930 | 0.0969 | 0.1639 | 0.0167 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-si-20 | standard | 3 | 0.2107 | 0.1033 | 0.1318 | 0.0273 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-si-20 | standard | 6 | 0.2189 | 0.1031 | 0.1235 | 0.0354 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-si-20 | strict | 2 | 0.1352 | 0.2414 | 0.1462 | 0.0188 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-si-20 | strict | 3 | 0.1756 | 0.3042 | 0.1100 | 0.0389 | n/a | n/a | n/a | n/a | n/a | n/a |
| material-si-20 | strict | 6 | 0.2525 | 0.4115 | 0.1031 | 0.1101 | n/a | n/a | n/a | n/a | n/a | n/a |
| si-elongated | standard | 2 | 0.1903 | 0.0800 | 0.1481 | 0.0518 | 0.7104 | 0.0386 | 41.10 | 37.21 | 24.23 | 0.615 |
| si-elongated | standard | 3 | 0.2035 | 0.0852 | 0.1320 | 0.0813 | 0.7104 | 0.0386 | 40.78 | 32.97 | 28.00 | 0.602 |
| si-elongated | standard | 6 | 0.2141 | 0.0900 | 0.1304 | 0.1233 | 0.7104 | 0.0386 | 38.83 | 25.51 | 28.00 | 0.612 |
| si-elongated | strict | 2 | 0.0860 | 0.2413 | 0.0889 | 0.0459 | 0.5986 | 0.0615 | 30.48 | 30.01 | -100.00 | 1.009 |
| si-elongated | strict | 3 | 0.1117 | 0.2982 | 0.0803 | 0.0963 | 0.5986 | 0.0615 | 18.11 | 11.16 | -100.00 | 0.928 |
| si-elongated | strict | 6 | 0.1343 | 0.3528 | 0.0790 | 0.1946 | 0.5986 | 0.0615 | 5.42 | -15.24 | -100.00 | 0.946 |
| si-k222s | standard | 2 | 0.1290 | 0.0666 | 0.1127 | 0.0169 | n/a | n/a | n/a | n/a | n/a | n/a |
| si-k222s | standard | 3 | 0.1414 | 0.0720 | 0.0899 | 0.0263 | n/a | n/a | n/a | n/a | n/a | n/a |
| si-k222s | standard | 6 | 0.1507 | 0.0725 | 0.0827 | 0.0394 | n/a | n/a | n/a | n/a | n/a | n/a |
| si-k222s | strict | 2 | 0.0939 | 0.1554 | 0.0973 | 0.0186 | n/a | n/a | n/a | n/a | n/a | n/a |
| si-k222s | strict | 3 | 0.1226 | 0.1929 | 0.0756 | 0.0369 | n/a | n/a | n/a | n/a | n/a | n/a |
| si-k222s | strict | 6 | 0.1651 | 0.2466 | 0.0697 | 0.0903 | n/a | n/a | n/a | n/a | n/a | n/a |
| si-skew | standard | 2 | 0.1240 | 0.0527 | 0.0948 | 0.0160 | 0.3940 | 0.0080 | 31.08 | 28.48 | 27.27 | 0.669 |
| si-skew | standard | 3 | 0.1273 | 0.0536 | 0.0821 | 0.0221 | 0.3940 | 0.0080 | 33.24 | 29.07 | 27.27 | 0.668 |
| si-skew | standard | 6 | 0.1305 | 0.0536 | 0.0796 | 0.0259 | 0.3940 | 0.0080 | 33.08 | 27.98 | 27.27 | 0.669 |
| si-skew | strict | 2 | 0.0974 | 0.1407 | 0.0958 | 0.0187 | 0.5014 | 0.0214 | 33.42 | 32.58 | -49.07 | 0.798 |
| si-skew | strict | 3 | 0.1172 | 0.1667 | 0.0745 | 0.0331 | 0.5014 | 0.0214 | 28.52 | 25.12 | -49.07 | 0.757 |
| si-skew | strict | 6 | 0.1426 | 0.1991 | 0.0665 | 0.0675 | 0.5014 | 0.0214 | 18.59 | 9.01 | -49.07 | 0.814 |
| si-slab | standard | 2 | 0.2723 | 0.1233 | 0.2164 | 0.0554 | 0.9887 | 0.0396 | 38.10 | 35.10 | 27.14 | 0.641 |
| si-slab | standard | 3 | 0.2950 | 0.1308 | 0.1925 | 0.0764 | 0.9887 | 0.0396 | 37.46 | 32.44 | 28.00 | 0.622 |
| si-slab | standard | 6 | 0.3088 | 0.1361 | 0.1819 | 0.1264 | 0.9887 | 0.0396 | 36.62 | 26.76 | 28.00 | 0.634 |
| si-slab | strict | 2 | 0.1417 | 0.3223 | 0.1447 | 0.0527 | 0.9751 | 0.0719 | 37.58 | 36.83 | -100.00 | 0.967 |
| si-slab | strict | 3 | 0.1858 | 0.3970 | 0.1241 | 0.1043 | 0.9751 | 0.0719 | 27.51 | 22.52 | -100.00 | 0.863 |
| si-slab | strict | 6 | 0.2391 | 0.4839 | 0.1216 | 0.2274 | 0.9751 | 0.0719 | 13.38 | -2.39 | -100.00 | 0.866 |
| si16-gamma | standard | 2 | 15.0268 | 8.2411 | 13.9716 | 1.4273 | 66.9023 | 1.6624 | 44.34 | 43.61 | 20.88 | 0.448 |
| si16-gamma | standard | 3 | 16.3984 | 8.5931 | 10.5356 | 3.2371 | 66.9023 | 1.6624 | 46.90 | 43.46 | 28.95 | 0.448 |
| si16-gamma | standard | 6 | 14.7419 | 7.6202 | 7.9977 | 5.2800 | 66.9023 | 1.6624 | 54.62 | 48.02 | 28.95 | 0.454 |
| si16-gamma | strict | 2 | 4.1416 | 23.7121 | 4.8171 | 1.3368 | 104.2429 | 7.8719 | 68.66 | 69.67 | 50.29 | 1.444 |
| si16-gamma | strict | 3 | 6.1703 | 30.4303 | 4.1543 | 3.0138 | 104.2429 | 7.8719 | 60.90 | 60.96 | 40.23 | 1.068 |
| si16-gamma | strict | 6 | 18.4341 | 60.3718 | 6.0304 | 24.4763 | 104.2429 | 7.8719 | 18.62 | 2.50 | 12.83 | 0.814 |
| si8-gamma | standard | 2 | 2.2017 | 0.9633 | 1.8802 | 0.2317 | 9.2638 | 0.1920 | 45.54 | 44.19 | 18.93 | 0.490 |
| si8-gamma | standard | 3 | 2.1908 | 1.0058 | 1.4098 | 0.4713 | 9.2638 | 0.1920 | 50.28 | 46.30 | 28.57 | 0.493 |
| si8-gamma | standard | 6 | 2.2411 | 1.0199 | 1.3414 | 0.6140 | 9.2638 | 0.1920 | 50.32 | 44.83 | 28.57 | 0.497 |
| si8-gamma | strict | 2 | 0.6738 | 3.7319 | 0.8469 | 0.2321 | 13.8702 | 1.0970 | 62.13 | 63.36 | 41.18 | 3.524 |
| si8-gamma | strict | 3 | 1.2223 | 4.8868 | 0.8354 | 0.5213 | 13.8702 | 1.0970 | 49.93 | 50.12 | -23.09 | 1.854 |
| si8-gamma | strict | 6 | 2.2068 | 8.3010 | 0.7544 | 3.4705 | 13.8702 | 1.0970 | 18.80 | 1.57 | 0.82 | 0.895 |
| syn20-skew | standard | 2 | 4.1792 | 3.0184 | 4.3560 | 0.4689 | 30.9931 | 1.2068 | 62.72 | 62.66 | 24.55 | 0.349 |
| syn20-skew | standard | 3 | 5.5304 | 4.1391 | 3.7273 | 0.9556 | 30.9931 | 1.2068 | 56.77 | 55.43 | 28.71 | 0.368 |
| syn20-skew | standard | 6 | 7.2129 | 5.7827 | 2.8900 | 3.8143 | 30.9931 | 1.2068 | 48.74 | 38.82 | 28.71 | 0.513 |
| syn20-skew | strict | 2 | 7.8000 | 8.9602 | 8.4850 | 0.7601 | 31.9452 | 2.0388 | 20.97 | 23.48 | -14.08 | 5.696 |
| syn20-skew | strict | 3 | 4.3750 | 10.1462 | 3.0523 | 1.0930 | 31.9452 | 2.0388 | 44.99 | 45.07 | 27.27 | 3.586 |
| syn20-skew | strict | 6 | 6.2700 | 17.5962 | 2.0122 | 6.2592 | 31.9452 | 2.0388 | 18.99 | 5.43 | 0.28 | 0.718 |

## SCF outcomes

| Case | Profile | Seed | Original SCF outcome |
|---|---|---:|---|
| si8-gamma | standard | 1 | Converged |
| si8-gamma | standard | 2 | Converged |
| si8-gamma | standard | 3 | Converged |
| si8-gamma | standard | 4 | Converged |
| si8-gamma | standard | 5 | Converged |
| si16-gamma | standard | 1 | Converged |
| si-skew | standard | 1 | Converged |
| si-elongated | standard | 1 | error:Diverged { reason: "density residual exceeded growth limit" |
| si-slab | standard | 1 | NotConverged |
| syn20-skew | standard | 1 | NotConverged |
| si-k222s | standard | 1 | Converged |
| material-si-20 | standard | 1 | Converged |
| material-mgo-40 | standard | 1 | Converged |
| si8-gamma | strict | 1 | Converged |
| si8-gamma | strict | 2 | Converged |
| si8-gamma | strict | 3 | Converged |
| si8-gamma | strict | 4 | Converged |
| si8-gamma | strict | 5 | Converged |
| si16-gamma | strict | 1 | Converged |
| si-skew | strict | 1 | Converged |
| si-elongated | strict | 1 | NotConverged |
| si-slab | strict | 1 | NotConverged |
| syn20-skew | strict | 1 | NotConverged |
| si-k222s | strict | 1 | Converged |
| material-si-20 | strict | 1 | Converged |
| material-mgo-40 | strict | 1 | Converged |
| syn20-skew | standard | 0 | NotConverged |
| syn20-skew | strict | 0 | NotConverged |

## Final gates

| Command / evidence | Exit | Result |
|---|---:|---|
| `scripts/check.sh` (`P37A-evidence/check-pass.log`) | 0 | build, test, fmt, clippy all PASS |
| `cargo build --offline --release -p dft-cli` | 0 | production executable built; SHA-256 in `cli-sha256.txt` |
| `bench_scf.py --require-identical`, five cases at T=1,2,4,8,12 | 0 each | 25/25: 0 Ha and identical density hashes |
| `python3 -I P37A-evidence/compare_control.py` | 0 | 4 paired batches; 1,176 scientific records identical |
| `python3 -I P37A-evidence/analyze.py` | 0 | 54 case/profile/capacity rows, 67,257 steps |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | existing generated README evidence remains fresh |
| Final ignored release runner, `P37A_CONTROL=1`, uniform A52 standard (`final-wrapper.log`) | 0 | current wrapper, source-identity guard and control path executed |

Expanded argument lists, per-command exit codes and log paths are stored in
`P37A-evidence/final-commands.json` and the individual `*.command.json` records.
All four final measurement batches exit 0; failed SCF or diagnostic convergence
is reported explicitly in the tables rather than hidden by the harness exit.
The two full-profile batches and two uniform-A52 batches constitute the final
measurement dataset; the earlier pilot batches are archived separately.
