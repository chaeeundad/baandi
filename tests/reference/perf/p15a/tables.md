# P15A frozen-H replay measurements

Current native build; ordinary SCF solves only. Candidate counts are prefixes of the baseline trajectory.
QE-like columns use the ratcheted schedule on the **BAANDI** history and occupied/empty thresholds.
Censored solves are excluded from candidate sums and explicitly counted. No new SCF policy was run.

| Case | SCF status / iterations | Native steps / H vectors | QE-like steps / H vectors (0 allowed) | Censored | E_res < 5e-9 Ha first iteration |
| --- | --- | ---: | ---: | ---: | ---: |
| k444s | converged / 10 | 1686 / 10290 | 740 / 5631 | 0 | 7 |
| si8 | converged / 15 | 283 / 2063 | 54 / 1000 | 0 | 11 |
| syn16-gamma | converged / 10 | 65 / 514 | 28 / 350 | 1 | 7 |
| syn16-k211 | exhausted / 100 | 622 / 5973 | 230 / 3389 | 92 | not reached |
| syn20-skew | converged / 6 | 13 / 205 | 8 / 140 | 2 | 4 |
| syn4-elongated | converged / 9 | 28 / 150 | 12 / 76 | 3 | 5 |
| syn4-slab | exhausted / 100 | 539 / 2506 | 192 / 1372 | 21 | not reached |

## Dense spectral-completeness check

Counts below cover ordinary solves. A mismatch is max sorted eigenvalue difference > 1e-8 Ha.
Dense spectra pass an independent true-residual check <= 1e-10 Ha.

| Case | Dense checks | Tight replay mismatches | Maximum tight vs dense difference (Ha) |
| --- | ---: | ---: | ---: |
| k444s | 5 | 0 | 7.27196081e-15 |
| si8 | 3 | 1 | 0.110723486 |
| syn16-gamma | 10 | 1 | 8.22560551 |
| syn16-k211 | 200 | 2 | 8.22560551 |
| syn20-skew | 6 | 1 | 2.95886039 |
| syn4-elongated | 9 | 1 | 0.380087693 |
| syn4-slab | 100 | 1 | 0.523220374 |

## Equal occupied-eigenvalue accuracy on selected k=0 snapshots

Best scanned prefix on the native trajectory, selected in hindsight at error <= 1e-10 Ha
against a dense oracle when available, otherwise the tight all-root replay. Counts include the initial reduced solve.
A dash means no scanned stop delivered that accuracy; this includes missing low roots, not just a short trace.

| Case | SCF iteration | Residual steps / H vectors / allowed | Ritz steps / H vectors / allowed |
| --- | ---: | ---: | ---: |
| k444s | 1 | - | - |
| k444s | 5 | - | - |
| k444s | 10 | 1 / 8 / 0 | 2 / 13 / 0 |
| si8 | 1 | - | - |
| si8 | 8 | 13 / 140 / 1 | 10 / 128 / 5 |
| si8 | 15 | 1 / 20 / 0 | 2 / 37 / 0 |
| syn16-gamma | 1 | - | - |
| syn16-gamma | 5 | - | - |
| syn16-gamma | 10 | 1 / 20 / 0 | 2 / 36 / 0 |
| syn16-k211 | 1 | - | - |
| syn16-k211 | 50 | - | - |
| syn16-k211 | 100 | - | - |
| syn20-skew | 1 | - | - |
| syn20-skew | 3 | - | - |
| syn20-skew | 6 | - | - |
| syn4-elongated | 1 | 2 / 12 / 0 | - |
| syn4-elongated | 5 | 3 / 15 / 0 | 3 / 15 / 2 |
| syn4-elongated | 9 | 1 / 8 / 0 | 2 / 9 / 0 |
| syn4-slab | 1 | - | - |
| syn4-slab | 50 | - | - |
| syn4-slab | 100 | - | - |

## Per-root rate summary (k=0, tight all-root trajectory)

Each entry is the range across roots of the first / sustained step below the residual threshold.
Sustained means no subsequent crossing in the recorded trace, not a locking guarantee.

| Case | SCF iteration | 1e-2 Ha | 1e-3 Ha | 1e-4 Ha | 1e-5 Ha | 1e-6 Ha |
| --- | ---: | --- | --- | --- | --- | --- |
| k444s | 1 | 4-5 / 4-5 | 6-8 / 6-8 | 8-10 / 8-10 | 9-13 / 9-13 | 11-15 / 11-15 |
| k444s | 5 | 1-1 / 1-1 | 2-3 / 2-3 | 4-5 / 4-5 | 5-7 / 5-7 | 7-9 / 7-9 |
| k444s | 10 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-1 | 1-4 / 1-4 |
| si8 | 1 | 5-9 / 5-42 | 7-12 / 7-44 | 8-16 / 8-46 | 10-25 / 10-48 | 12-28 / 12-50 |
| si8 | 8 | 1-1 / 1-7 | 1-1 / 1-30 | 2-32 / 2-32 | 5-34 / 5-34 | 6-36 / 6-36 |
| si8 | 15 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-20 | 1-1 / 1-40 | 1-1 / 1-42 |
| syn16-gamma | 1 | 1-2 / 1-3 | 1-2 / 1-3 | 1-2 / 1-3 | 1-2 / 1-3 | 1-2 / 1-3 |
| syn16-gamma | 5 | 1-1 / 1-1 | 1-1 / 1-1 | 1-3 / 1-3 | 1-10 / 1-12 | 1-16 / 1-16 |
| syn16-gamma | 10 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-1 |
| syn16-k211 | 1 | 1-2 / 1-3 | 1-2 / 1-3 | 1-2 / 1-3 | 1-2 / 1-3 | 1-2 / 1-3 |
| syn16-k211 | 50 | 1-1 / 1-1 | 1-1 / 1-1 | 1-3 / 1-3 | 2-18 / 2-18 | 3-50 / 3-50 |
| syn16-k211 | 100 | 1-1 / 1-1 | 1-1 / 1-1 | 1-2 / 1-2 | 1-10 / 1-10 | 2-28 / 2-28 |
| syn20-skew | 1 | 1-2 / 1-3 | 1-2 / 1-3 | 1-2 / 1-3 | 1-2 / 1-3 | 1-2 / 1-3 |
| syn20-skew | 3 | 1-1 / 1-5 | 1-3 / 1-5 | 1-3 / 1-5 | 1-3 / 1-5 | 1-3 / 1-5 |
| syn20-skew | 6 | 1-1 / 1-5 | 1-1 / 1-5 | 1-1 / 1-5 | 1-1 / 1-5 | 1-3 / 1-6 |
| syn4-elongated | 1 | 1-2 / 1-2 | 1-2 / 1-2 | 1-2 / 1-2 | 1-2 / 1-2 | 1-2 / 1-2 |
| syn4-elongated | 5 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-1 | 1-3 / 1-3 | 2-4 / 2-4 |
| syn4-elongated | 9 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-1 | 1-1 / 1-1 |
| syn4-slab | 1 | 3-5 / 3-8 | 3-5 / 3-8 | 3-5 / 3-8 | 3-5 / 3-8 | 3-5 / 3-8 |
| syn4-slab | 50 | 1-1 / 1-1 | 2-4 / 2-4 | 4-6 / 4-32 | 5-38 / 5-45 | 7-47 / 7-47 |
| syn4-slab | 100 | 1-2 / 1-13 | 3-5 / 3-16 | 5-20 / 5-20 | 6-26 / 6-26 | 9-32 / 9-34 |

Full per-root rate curves, scan hits, residuals at each stop, matched occupied-eigenvalue errors,
and SCF schedules are in `analysis.json`. Raw per-step values, Ritz changes, residuals,
unconverged-root counts and H vector/block/restart counts are in each `*-traces.json.gz`.
Reference: full-space dense spectrum where available; otherwise a tight replay (1e-10 Ha true residual).
Some tight replays miss low roots: residual convergence alone is not spectral completeness. Eigenvalue
differences at 1e-12 Ha are floating-point comparisons, not a certified absolute error bound.
A 5-root early allowance is a separate diagnostic; the `qe_cap` column applies it only
after 20 expansion steps (BAANDI trace row 21). No repeated-call QE implementation is inferred.
First-iteration redo flags use the stored native residual. A redo changes the residual
and starting vectors; its actual work and subsequent SCF trajectory are not inferred.
