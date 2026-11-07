# P32A real Gamma algebra investigation

Baseline: `40118fa62d2924de0500e21ed8e93b40035ea34f`. Worktree-local raw
commands, inputs, binary and measurements: `P32A-evidence/`. No QE source read.
This investigation is derived from real-gamma-design section 8 and its own-words
study notes. Mac timings are observations; dev08 acceptance belongs to the
orchestrator. Capacity, FFT pairing, complex kernels and defaults are untouched.

## Premise and candidates

D0 relayed dev08 evidence reports algebra at 60.2% of working samples at T=1
(495 working samples), 66.4% at T=8 (292). This changes the premise of A49:
we investigate the real backend's 16-cubed microcalls, not complex GEMM tile
sizes or scheduling thresholds. Fixed 256-row outer tiles and ordered partial
sums remain. Unlike the rejected A52, no factorization or orthogonalization
algorithm is changed.

The current real reduced-matrix update already computes one rectangle plus the
new square and fills the transpose. There is no second rectangle to remove.

Source inspection of vendored faer 0.24.4 `linalg/matmul/mod.rs` shows general
large native-f64 GEMM uses gemm/private-gemm-x86 without a caller scratch API.
The gemm-common implementation allocates packing/bookkeeping memory per call;
warming its thread-local slab alone is insufficient. Sequential column-major
and row-major matrix-vector branches use SIMD with no heap allocation. The
first candidate therefore passed each whole 256-row panel and one RHS column
through the safe public faer matmul API. It passed accuracy but lost rotation
throughput and was rejected. Straight slice products also lost. General large
GEMM is timed only as a rejected allocating comparison, not used in production.
No new dependency or unsafe is needed.

The retained experiment uses 4x4x256 rotation calls (at most 4096 MACs), replacing
16x16x16 calls. Gram retains 16x16x16: tested alternatives were slower or too
close to justify changing it. The outer 256-row tiling, ordered Gram sums,
parallel threshold and worker scheduling remain unchanged. This is a bounded
backend-call improvement, **not** delivery of a whole-panel general GEMM with
caller-owned scratch. That larger improvement needs backend API work.

Other retained changes use contiguous column slices for finite checks, partial
sum assembly, corrections, residuals and projection/MGS. They preserve scalar
association. Reduced-prefix copies use column slices; symmetric filling uses
split slices. Ritz coefficients are borrowed prefix views. Two projection passes
swap their equally sized buffers rather than copying the whole intermediate
block. Owned solver results and EVD still allocate at existing orchestration
boundaries; this is not an allocation-free entire solver claim.

Allocation reasoning: every backend call has `m*n*k <= 4096`, selecting faer's
existing nano-GEMM branch (or earlier allocation-free sequential matvec/rank-one
branches). No packing storage is allocated by these calls. Workspace partials
are allocated at construction, and the new slice loops/views/swaps do not
allocate. This source audit is not a runtime global-allocator count.

## Kernel candidate observations

`panel-shapes.log`: ordinary release build, seven batches of 100 calls, median
seconds/call, one test thread. No other task-owned build or benchmark overlapped.
Four shapes span reduction widths 20, 70, 120 and 216, output widths 20 and 36.
For rotations, old vs 4x4x256: 7.395/7.458, 26.466/23.749,
44.070/40.552, 142.330/129.288 microseconds. The smallest shape is neutral;
larger three improve approximately 8–10%. This is not a universal speed claim.
Gram alternatives did not improve consistently; its existing shape stays.

The full-panel column candidate takes 9.163/36.689/65.597/211.821 microseconds
for those rotations and is rejected. Straight slice loops take
8.866/37.465/66.938/215.004 microseconds and are rejected. Large allocating GEMM
is faster (3.857/11.367/19.913/71.955 microseconds), but violates the fixed
no-allocation contract. Other allocation-free call shapes (256x4x4, 64x4x16,
32x8x16, 8x4x128) and Gram results are retained in the raw log.

The old implementation is preserved in `tests/real_algebra.rs`. Relative
Frobenius error is used because entrywise relative errors are undefined for
cancellation zeros. Gates remain 1e-13 against old and 1e-14 against naive;
only summation association inside tiles changes. The original unit tests and
fixed-operator tolerances are unmodified.

## Measurement limitations discovered before editing

`sample` of a live child process exits 255 in this sandbox, with "cannot examine
process ... try running with sudo" (`sample-preflight.log`). Permissions are
not escalatable. Thus no local sampling-share claim or claim of achieving the
30% target is possible here. The orchestrator must obtain before/after samples
with at least 30 observations at each named site. No timing instrumentation is
substituted for a sampling profile.

Completed commands and results are recorded below.

## End-to-end Mac observations

Five repetitions per row; wall includes CLI setup/output. Work counts are main
SCF history totals (certificate work remains separately available in results).
Every run converged. All energies, densities and scientific histories are exact
across T=1,2,4,8 within each build; timings are excluded from history comparison.
Before/after histories differ in every on row. These are trajectory-changing
observations, not equal-work kernel speedups. No task-owned build overlapped
the timing series. Shared host load averaged 11.32/7.66/6.46 at baseline start
and 5.79/6.45/6.27 at its end; other tasks were not controlled. The unchanged
complex-off timings also improved, demonstrating environmental noise.

| Case / profile / threads | Before median (s) | After median (s) | SCF before/after | Davidson before/after | H vectors before/after | Same trajectory |
|---|---:|---:|---:|---:|---:|---|
| si8-strict-T1-on | 1.817735 | 1.270901 | 15/15 | 320/253 | 2126/2001 | False |
| si8-strict-T2-on | 1.656811 | 1.124341 | 15/15 | 320/253 | 2126/2001 | False |
| si8-strict-T4-on | 1.423391 | 1.030011 | 15/15 | 320/253 | 2126/2001 | False |
| si8-strict-T8-on | 1.664327 | 1.051500 | 15/15 | 320/253 | 2126/2001 | False |
| si8-standard-T1-on | 0.607639 | 0.502983 | 15/17 | 81/89 | 808/877 | False |
| si8-standard-T2-on | 0.539966 | 0.455660 | 15/17 | 81/89 | 808/877 | False |
| si8-standard-T4-on | 0.569663 | 0.428309 | 15/17 | 81/89 | 808/877 | False |
| si8-standard-T8-on | 0.573425 | 0.435849 | 15/17 | 81/89 | 808/877 | False |
| si16-strict-T1-on | 9.551941 | 5.336214 | 15/14 | 329/264 | 3545/3254 | False |
| si16-strict-T2-on | 8.614247 | 4.410742 | 15/14 | 329/264 | 3545/3254 | False |
| si16-strict-T4-on | 7.809494 | 3.877212 | 15/14 | 329/264 | 3545/3254 | False |
| si16-strict-T8-on | 7.526996 | 3.692307 | 15/14 | 329/264 | 3545/3254 | False |
| si16-standard-T1-on | 3.202838 | 1.980327 | 18/16 | 80/88 | 1714/1521 | False |
| si16-standard-T2-on | 2.389246 | 1.708113 | 18/16 | 80/88 | 1714/1521 | False |
| si16-standard-T4-on | 2.222637 | 1.554919 | 18/16 | 80/88 | 1714/1521 | False |
| si16-standard-T8-on | 2.161663 | 1.499449 | 18/16 | 80/88 | 1714/1521 | False |

T=1 complex-off comparison (five-run medians, same measurements):

| Case/profile | Before off/on ratio | After off/on ratio | Before off (s) | After off (s) |
|---|---:|---:|---:|---:|
| si8/strict | 1.553 | 2.020 | 2.823060 | 2.567339 |
| si8/standard | 1.283 | 1.382 | 0.779863 | 0.695057 |
| si16/strict | 1.136 | 1.738 | 10.854805 | 9.273093 |
| si16/standard | 1.250 | 1.657 | 4.002190 | 3.280915 |

## Slice/copy ablation

Same production changes except both call shapes use the original 16x16x16.
Five T=1 runs per case, after the final series. Scientific output matches the
baseline exactly; this isolates arithmetic-preserving overhead removal.

| Case | Before (s) | Slice/copy only (s) | Full candidate (s) | SCF / steps / H, slice-only |
|---|---:|---:|---:|---|
| si8-strict-T1-on | 1.817735 | 1.412642 | 1.270901 | 15 / 320 / 2126 |
| si8-standard-T1-on | 0.607639 | 0.448739 | 0.502983 | 15 / 81 / 808 |
| si16-strict-T1-on | 9.551941 | 6.811003 | 5.336214 | 15 / 329 / 3545 |
| si16-standard-T1-on | 3.202838 | 2.008362 | 1.980327 | 18 / 80 / 1714 |

The full candidate adds a strict-mode gain but loses about 12% versus slice-only
on Si8 standard (81 to 89 Davidson steps, 15 to 17 SCF). Si16 standard is
near-neutral versus slice-only. This mixed result is retained explicitly for
dev08 evaluation, not presented as a uniformly accepted optimization.

A failed initial ablation setup assertion left the full candidate unchanged;
those measurements are correctly retained as `after-repeat-runs`, not counted
as slice-only. The corrected ablation has its own binary hash and build log.

## Accuracy and validation

No tolerance or existing test was edited. Expanded kernel tests cover partial
and full panels, dimensions beyond the parallel guard, reduction width 257,
strided input/output views and zero reduction. Maximum relative Frobenius error
against the retained old implementation is 1.0059200815495567e-15, against naive
7.932507683289761e-16. Every output bit matches at T=1,2,4,8. Original real
adapter tests pass too. Fixed-operator tests: 8 passed including ignored release
cases; maximum dense eigenvalue error 2.853273173286652e-14 Ha, below 1e-12 Ha.
The test's timing prints overlapped other correctness work and are not used for
performance claims.

Complex/option-off: all 25 gamma/k222s/k222/k444s/si8 comparisons at
T=1,2,4,8,12 pass `bench_scf --require-identical`: delta energy 0 Ha and identical
density hashes. These timings are correctness-only (may overlap builds).

Strict real versus complex R5 results:

| Fixture | Absolute energy difference (Ha) | SCF complex/real |
|---|---:|---:|
| synthetic 0 | 6.661338147750939e-16 | 10/10 |
| synthetic 1 | 4.440892098500626e-16 | 11/11 |
| synthetic 2 | 3.685940441755520e-14 | 11/11 |
| synthetic 3 (slab) | 8.881784197001252e-16 | 13/10 |
| Si primitive | 1.332267629550188e-13 | 10/12 |
| Si8 | 3.552713678800501e-13 | 15/15 |

All six strict energy differences pass 1e-10 Ha; complete thread comparisons
also pass. The material test nevertheless exits 101 for its pre-existing extra
SCF-count condition, primitive 10 versus 12 (documented in R5/R6). It is not
reported as a passing test. No count threshold was altered.

Standard error against the shared tighter certified complex reference:

| Fixture | Complex error (Ha) | Real error (Ha) | Absolute 5e-9 cap | Relative cap |
|---|---:|---:|---|---|
| synthetic 0 | 1.727831211439934e-10 | 1.022355533564223e-10 | pass | pass |
| synthetic 1 | 2.382233299513814e-10 | 2.382233299513814e-10 | pass | pass |
| synthetic 2 | 4.006044385107543e-10 | 2.874256388452068e-10 | pass | pass |
| synthetic 3 (slab) | 4.598657676879725e-9 | 6.009215569591220e-9 | **FAIL** | pass |
| Si primitive | 1.080682210385930e-10 | 2.336077997711072e-10 | pass | pass |
| Si8 | 1.000891813873750e-9 | 1.265721749632576e-9 | pass | pass |

Relative cap is `max(3*complex_error, 1e-10 Ha)`. Reference plateaus all pass
1e-12 Ha. Standard synthetic suite exits 101 (10 passed, one failed); material
standard test passes. The slab failure matches the documented pre-existing
R6 limitation, now 6.0092155696e-9 versus historical 6.0092151255e-9 Ha.
The architect's R6 disposition treats the absolute cap as a reported observation;
this report still marks it failed and does not claim all fixed gates passed.
No synthetic fixture, tolerance or convergence algorithm is tuned to it.

Both actual before/after CLI sampling attempts also exit 255; the SCF runs
exit 0. Commands and diagnostics are in `profile-status.json` and
`profile-{before,after}.log`. There is no measured local algebra share, and
neither the <30% share target nor named-site sample counts are established.

## Reproduction and dev08 request

All commands run from the ticket worktree with the pinned toolchain and its
single `target/`. No native flags, dependency changes or unsafe were introduced.
The external UPF is `../pseudo/Si.pz-vbc.UPF`, SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Binary/source hashes and machine information: `P32A-evidence/metadata.json`.

- `cargo build --release -p dft-cli`; copy the initial executable to
  `P32A-evidence/before`, the candidate to `P32A-evidence/after`.
- `python3 -I P32A-evidence/measure.py before-runs P32A-evidence/before` and
  the analogous `after-runs .../after` command: 20 rows, five runs each.
- `python3 -I P32A-evidence/identity.py before P32A-evidence/before` and
  `... identity.py after P32A-evidence/after`: five commands per build.
- `cargo test --release -p dft-numerics --test real_algebra panel_alternatives -- --ignored --nocapture --test-threads=1`.
- `python3 -I P32A-evidence/validate.py`: exact constituent commands, exit codes
  and log paths are recorded in `validation.json`, including failures.
- `python3 -I P32A-evidence/profile.py`: both denied attachment attempts.

The initial benchmark-test compile failed with a function-pointer lifetime
inference error (exit 101, `kernels.log`); explicit higher-ranked lifetimes
fixed it. No failing numerical check was weakened.

For dev08, compare the baseline commit above, slice/copy-only ablation, and the
full uncommitted patch after review. Use identical pure/native flags and the
same external UPF, pinned idle cores at T=1,2,4,8,12. Repeat the two-cell strict
and standard timing matrix, record medians/spread, SCF steps, Davidson steps,
H vectors and delivered error. Obtain caller-resolved before/after working-thread
samples at T=1 and T=8 (at least 30 per named site), separating algebra from
FFT/nonlocal/waiting. Compare against QE 7.6's best real-Gamma configuration at
equal accuracy. This task performed no remote runs or server changes.

The performance target's T=1 real/complex margin improves in this Mac series,
but the algebra-share target is unmeasured. General large GEMM remains blocked
by the no-allocation backend API contract. The retained bounded-call change has
a Si8-standard regression against the slice-only ablation; dev08 must decide
whether to retain it or take only the arithmetic-preserving part. There is no
new accepted default or README speed claim. Files remain uncommitted.


## Final verification bookkeeping

The first full `scripts/check.sh` exited 1: build, test and fmt passed, but
Clippy required `as_chunks::<4>()` in a rejected test-only prototype. That
prototype was corrected without changing production arithmetic or weakening
checks. Its release benchmark was rerun (`panel-final-test-code.log`, exit 0)
and the candidate decisions were unchanged. Both the original and corrected
measurements are preserved. README generation and freshness tests exited 0;
the generated README/figures were unchanged, because this is not new accepted
reference evidence. `artifacts/*` was intentionally not modified, per ticket.
No plumbing-only edits outside the named ownership were necessary.

Final `scripts/check.sh` retry: **exit 0**, PASS build, test, fmt, clippy
(`P32A-evidence/check-final.log`). Final disposition: **BLOCKED for full task
acceptance**, with a reviewable uncommitted implementation and completed local
measurements. Missing sampling evidence and the explicitly reported R5/R6
limitations must not be represented as passed gates.
