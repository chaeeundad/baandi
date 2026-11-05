# P15C: final reported-state spectral audit

**The tested realistic final samples are complete:** Si8 occupied error
1.39e-14 Ha, k444s sampled occupied error 9.55e-15 Ha, and all six PseudoDojo
Si/MgO cutoffs have final ordinary occupied error <=3.60e-14 Ha. Their fresh
certificate samples are also complete (maximum 7.21e-14 Ha). These findings do
not invalidate the corresponding stored QE comparisons.

**Residual-only certification remains unsafe when the guard is off.** The
20-atom skew synthetic fixture still demonstrates a legacy false success:
6 iterations, occupied-root discrepancy 1.515 Ha and fixed-H occupied band-sum
excess 27.958 Ha. `spectral_guard="final"` detects it, repairs the trial space,
and ultimately reports nonconvergence / inability to verify the reported state.
It never labels the legacy incomplete state as spectrally verified.

## Default decision and interface

```toml
[scf]
spectral_guard = "final" # or "off"; default is "off"
```

Default **off** follows the architect's conditional approval. A final-only audit
adds 90 H vector applications on Si8 (3.974% of original SCF work) and 2618 on
k444s (23.415%). Neither meets the <3% target. The previously reported 3.93% and
23.41% already came from certificate-only checking; auditing every iteration
was not their cause. The small Si8 change from 89 to 90 columns reflects checking
the actual reported ordinary Hamiltonian instead of the fresh certificate's H.

No optional first-iteration trigger is installed. Root clustering and small
Ritz changes are not proofs of a missing state; adding a tuned threshold could
miss the invariant-subspace failure or charge benign degeneracies. The bounded
independent search retains its 1e-10 Ha residual target rather than declaring
success from an unconverged short search. Failure may incur extra work.

Off mode preserves the legacy top-level results schema and schema_version.
`scf.spectral_guard="off"` and the console explicitly say spectral completeness
was not verified. Final mode adds `spectral_completeness` without changing
schema_version (1 unpolarized, 2 collinear). CLI tests check the additive block,
unchanged physical data, both spin channels, and rejection of unknown modes.

## Exact audit ordering and failure handling

1. Run the ordinary SCF map and unchanged fixed-point eligibility checks.
2. Run the existing fresh fixed-point certificate. No spectral H work is spent
   if its density/occupied-eigenvalue checks fail.
3. Reconstruct the Hamiltonian that produced the **reported ordinary**
   eigenpairs, using that map's input density. Check every k point (and both spin
   channels). Do not compare ordinary eigenpairs against the different fresh
   certificate Hamiltonian. The existing fixed-point tolerances still relate
   input/output densities; they are not changed to force exact self-consistency.
4. On agreement, retain the original eigenvalues, eigenvectors and density
   bit for bit. A normal complete run has exactly one final audit per block.
5. On a detected missing root, re-solve from an independent seeded initial
   space of 2*nbands, expanding to 3*nbands if needed. Small problems use a dense
   reference to validate the repair. A replacement cannot certify the old
   candidate, including when only empty states changed. Repaired vectors enter
   another ordinary iteration; stale mixer history and growth baselines reset.
   Any later successful candidate needs its own final audit.
6. Failed independent probes never produce a passing flag. The retained outcome
   says `could_not_verify_reported_state` and exits nonzero. An exhausted run
   also has `passed=false`, even when an earlier repair was performed.

For dimension <=128 the check assembles and diagonalizes the full operator.
Larger systems search independently in the orthogonal complement using
`P H P + shift Q Q*`. A detected lower state invokes independent enlarged
full-block solves. This is **an independent seeded numerical check, not an
inertia count or a universal mathematical proof of completeness**; the results
block states that scope. Seeds are independent of thread counts. Scratch in
operator applications is preallocated; output and JSON remain outside kernels.

The cumulative block records method, passed, total retries, total extra H
columns, guard wall time, and ordered events with iteration/spin/k, method,
retries, discrepancy and failure reason. Counters include attempted columns in
failed operator blocks. Timing is segregated under `timings`; it is not part of
deterministic scientific comparisons. A successful check after earlier rescues
keeps the earlier events and cumulative retry count. Unit tests exercise
multiple recoveries followed by a failed probe. The real skew run retains its
iteration-6 rescue after 94 subsequent iterations.

Checkpoint v1 cannot preserve this new audit history and policy. Final-guard
sessions explicitly reject checkpoint export/restore, and the CLI encoder
rejects a manually supplied final-guard continuation. Off-mode checkpoint wire
layout stays unchanged and decodes its policy as off; option identity includes
the guard. This avoids silently losing the record or changing policy on resume.

## Same-H dense validation

Caps remain frozen at **1e-9 Ha** for sorted eigenvalues and occupied band sums.
The independent oracle applies the production H to a complete coordinate basis
then uses the full-space faer Hermitian EVD, not Davidson's reduced space.
Hermiticity is checked; requested dense roots have true residual <=1e-10 Ha.
This shares the existing faer adapter, not a second EVD implementation. The
source-only oracle is bounded at 4096 plane waves.

Original iteration scope and multiplicity matching are retained in
[`summary.json`](../../tests/reference/spectral/summary.json) and the adjacent
per-case traces. Unmatched approximate roots are not called missing roots unless
all returned values can be matched at the frozen cap. For the final skew
ordinary and certificate spectra, missing zero-based indices are **2,10,11,12**:
partial doublets [1,2], [9,10], and the entire doublet [11,12]. A degenerate
member at the requested-band boundary is not itself a missing-root failure.

New ignored CLI tests run the actual P13D inputs with external `DFT_NC_DIR`
files, verify their recorded hashes, compare final/off energy and density bytes,
and compare ordinary and fresh-certificate spectra with full-space EVD at
k indices 0, floor(nk/2), nk-1. All k points receive the final guard; dense
validation is explicitly sampled. No UPF contents are committed.

| Fixture | Final ordinary occupied max (Ha) | Fresh certificate occupied max (Ha) |
| --- | ---: | ---: |
| PseudoDojo Si 20 Ha | 1.09357e-14 | 1.04916e-14 |
| PseudoDojo Si 30 Ha | 2.58960e-14 | 1.45717e-14 |
| PseudoDojo Si 40 Ha | 1.72362e-14 | 2.40086e-14 |
| PseudoDojo MgO 40 Ha | 3.04201e-14 | 4.13003e-14 |
| PseudoDojo MgO 50 Ha | 3.59712e-14 | 5.60663e-14 |
| PseudoDojo MgO 60 Ha | 3.57492e-14 | 7.20535e-14 |

All requested roots, including buffers, pass the 1e-9 Ha cap in these samples;
the largest material buffer error is 2.83826e-11 Ha. The exact P13B Gaussian
spin fixture also passes its stored final dense comparisons; final-guard CLI
tests check both channels at T=1,2,4,8. No new QE run was executed or inferred.

The band sum is `sum_b f_b*(lambda_returned-lambda_dense)` at fixed H. Sampled
k-point values are not invented full-mesh integrals. The skew Gamma value is
27.958043613987726 Ha, **not** a claimed corrected SCF total energy: changing
occupied density changes Hartree/XC and the self-consistent Hamiltonian.

## Why the corrected skew run is hard

A source-only dense-SCF build replaces **every** eigensolve with full-space dense
diagonalization, starting from the same uniform density and using unchanged
mixing, occupations and tolerances. This is feasible: only **41 plane waves**
for 40 valence electrons, 20 occupied roots and 24 returned bands at 3 Ry.
All 100 dense solves have maximum eigen residual <=2.889e-14 Ha, yet the SCF
still exhausts 100 iterations. The final density residual is 0.0020786151
Bohr^-3/2 (minimum 0.0020786070), far above the original 1e-7 target.

| Dense SCF iteration | Density residual (Bohr^-3/2) | Occupied/empty gap (Ha) | Total energy (Ha) |
| ---: | ---: | ---: | ---: |
| 1 | 9.32750e-3 | 1.77636e-15 | 143.0275570526 |
| 6 | 2.16908e-3 | 1.84867e-3 | 143.0275223342 |
| 10 | 2.08197e-3 | 6.87573e-5 | 143.0274467641 |
| 50 | 2.07861e-3 | 3.38494e-8 | 143.0274438971 |
| 100 | 2.07862e-3 | 1.65385e-7 | 143.0274439139 |

The initial occupied boundary splits a doublet. The later gap again becomes
very small while the density residual stalls. This is consistent with an
unstable/nonunique fixed-integer occupation map near a partially occupied
manifold. It is evidence that rescue/mixer history is **not the sole cause**:
dense SCF fails without any rescue. It does not establish a unique converged
DFT ground state or identify a precise limit cycle. The nearly stationary energy
must not substitute for density convergence.

This is a soft, generated multi-projector potential at an extremely low cutoff,
not a validated material or a representative insulating benchmark. Retain it as
a spectral-completeness and fail-closed regression, not as a target for tuning
mixing, occupation or solver defaults. No smearing, fractional occupations,
tolerance relaxation or fixture-specific production parameters were added.

## Work and wall observations

T=1, pure Rust. Off H counts include ordinary and fixed-point-certificate work;
guard H columns are separately reported. Wall samples are overlapping Mac
observations, **not** a statistically resolved overhead/speedup measurement.

| Case | Off H columns | Extra guard H | Extra / off | Off/final iterations | Guard wall (s) |
| --- | ---: | ---: | ---: | --- | ---: |
| Si8 | 2265 | 90 | 3.974% | 15 / 15 | 0.11957 |
| k444s | 11181 | 2618 | 23.415% | 10 / 10 | 0.70251 |
| Skew 20 | 255 | 130 | 50.980% | 6 / 100 exhausted | 0.00153 |
| Elongated 4 | 171 | 29 | 16.959% | 9 / 9 | 0.000114 |
| Gamma 16 | 593 | 79 | 13.322% | 10 / 10 | 0.001046 |
| Slab 4 | 2506 | 0 | 0% | 100 / 100 exhausted | 0 |
| Mesh 16 | 5973 | 0 | 0% | 100 / 100 exhausted | 0 |

The skew final-policy SCF uses 3032 ordinary/certificate H columns plus 130 guard
columns. Its longer trajectory is not equal-state performance evidence versus
the old wrong state. Exhausted slab/mesh cases never pass the fixed-point gate,
so no final audit runs and their block explicitly says not verified.

Si8 and k444s off/final wall observations were 2.905/2.775 s and 2.537/3.275 s;
noise makes the first pair unsuitable for an overhead conclusion. Guard timers
and all raw pairs are retained. Request controlled measurements using
[the dev08 handoff](spectral-dev08-request.md), covering native and default builds.

## Verification and limits

- Build/test/fmt/clippy and release spectral solver regressions pass. The
  unchanged default `bench_scf` passes at 0 Ha energy difference with identical
  density hashes on all five established cases.
- Seven cases x two policies x T=1,2,4,8: **56 combinations** agree bit for bit
  in energy, density hashes, eigenvalues, scientific convergence history,
  spectral events, retries and H counts (wall timers excluded). Four additional
  spin thread settings are covered by the CLI test.
- The six ignored real NC dense tests pass; the existing robustness suite reports
  16 pass and 3 nonconverged, not 19 successful SCFs. README freshness passes.
- `bench_scf.py` remains unchanged: its original counters exclude the new audit.
  Total work is those counters **plus** `spectral_completeness.extra_operator_applications`.
- Full logs, source/binary hashes, failures and commands are in
  [`final/verification.json`](../../tests/reference/spectral/final/verification.json).
  Earlier experimental failures remain in the original verification manifest.
- No default-on or universal mathematical completeness claim is made. Off mode
  cannot prevent the demonstrated wrong-state certificate. Final mode has a
  declared independent-search scope for large operators and can fail closed.
  Physical-manifold initial overlap ranks and controlled dev08 wall times are
  not established here.

Changes remain uncommitted. README generator/dashboard updates remain with the
orchestrator because their files are outside the ticket's ownership; no accepted
README performance claim was introduced.
