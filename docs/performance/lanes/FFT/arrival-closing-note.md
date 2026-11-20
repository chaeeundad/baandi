# FFT arrival-skew decision (2026-11-20)

Decision: recommend closing the FFT optimization lane on current evidence.
Do not implement a static ownership change. No demonstrated mechanism supports
saving 0.15 s at si64 MPI8. This is an evidence-limited closure, **not** a finding
that OS jitter or NUMA has been proved. FFT_KEPT_FUSE stays default off, NEUTRAL.

## Evidence and interpretation

Orchestrator relay: source 45340f6, si64 MPI8, seven verified SCF iterations,
776 matching events/rank, complete trace, clock error 684 ns, zero measured
drift. Rank-summed pre-last-entry residence is 2.7990 s: tree_sum 1.3312 s,
alltoall 1.2136 s, allgatherv 0.2542 s. Mean exposure is 0.349875 s/rank,
3.39% of 10.33 s external wall. Fine contributes 2.0283 s (72.5%).
Skew median/mean/p90/max: 504.3/960.4/2409.6/9427.6 us.

Exposure is time inside a collective before the last rank enters; it can include
useful transfer/progress. It is not a direct idle-time measurement. Dividing its
rank sum by eight is a useful opportunity scale, **not a rigorous wall-time
ceiling**. Collective dependencies, release skew and transfer overlap determine
critical-path savings. Saving 0.15 s would be 42.9% of this mean-exposure scale;
we have no evidence that ownership redistribution can recover that fraction.

Rank 3 is last 220 times (28.4%) and experiences only 0.150 s; rank 2 is last
42 times and experiences 0.500 s. Last-arrival frequency says who arrives late,
not whether its immediate kernel is slower. No per-rank stage seconds or matched
per-call release/compute decomposition were supplied with this aggregate relay.

## Ownership audit

`ownership_audit.py` reconstructs the production partition rule without timing.
Assumptions: cubic Gamma sphere, 23871 plane waves, 72^3 grid, eight ranks,
default mirror stick ownership. The integer sphere |G|^2 <= 320 reproduces the
reported basis size exactly. This is an arithmetic reconstruction, not a runtime
manifest from dev08; a different cell or ownership override invalidates it.

| Rank | Sticks | Full coefficients / ORC rows | Transpose planes |
|---|---:|---:|---:|
| 0 | 125 | 2983 | 9 |
| 1 | 124 | 2984 | 9 |
| 2 | 126 | 2986 | 9 |
| 3 | 126 | 2986 | 9 |
| 4 | 126 | 2986 | 9 |
| 5 | 126 | 2982 | 9 |
| 6 | 126 | 2982 | 9 |
| 7 | 126 | 2982 | 9 |

Total 1005 sticks, 35 active z planes. Production `partition` greedily balances
coefficient loads, pairing mirror sticks. Nine transpose planes per rank have
no remainder. Mirror pairs keep both real ORC rows together, so ORC row counts
equal full coefficient counts here. Projector matrices have local ORC rows and
**all** projector channels (`DistributedRealGamma::new` loops over all atoms);
there is no last-rank atom/channel remainder explaining rank 3's GEMM load.

Rank 3 owns exactly the same counts as ranks 2 and 4. Row spread is four rows
(0.134% of mean); stick spread is two (1.59% of mean). Rank 0 is also often last
(167 calls) despite fewer sticks and rows. This rejects a unique excess-count
explanation for rank 3. Address order, cache behavior and placement remain
possible: equal counts do not establish equal kernel cost. Under a simplistic
linear work model, eliminating the full stick spread from the entire 0.7184 s
permutation group would represent only about 0.0114 s, not 0.15 s. This is an
order-of-magnitude diagnostic using earlier bad02ee timing, not a prediction or
an upper bound for the different 45340f6 SCF trajectory.

## Why solver.append appears late

The trace stores the last **completed** instrumented stage. STAGE_APPEND covers
several operations: correction/product column copies, workspace fitting and the
finite-value scan after H. In the lean path, `append_in_place` is timed and then
`border` computes the local overlap/Gram contribution before its tree_sum. The
border stage timer is recorded only after that collective returns. Thus an
arrival labelled solver.append can actually follow substantial unlabelled local
border work. The 166 last arrivals cannot be assigned to append copy cost.

Border reduction dimensions/column widths are shared across ranks; local matrix
work uses the ORC row counts above. No rank-3-only append border or remainder is
visible in this partition. The preceding collective's unequal release times can
also propagate into the next arrival. The existing parser separates
`previous_release` and `local_interval`; these matched values are needed before
blaming the last completed stage. The 120 unattributed events and stage_slot.26
must not be silently attributed to append or FFT.

## Disposition and reopening evidence

No default, numerical operation, collective order or ownership changes. No new
performance run requested. The earlier pass-removal opportunity has no identified
removable pass with supported savings; kept-copy fusion already failed its wall
gate. Payload reduction without wall improvement in NONLOCAL likewise does not
support another generic collective-volume candidate.

If the already queued si32 trace and original si64 results become available,
reopen only if matched events show persistent larger **local compute intervals**
on the late rank after subtracting previous-release differences. Compare
previous_stage_ns, work/shape, site and per-rank FFT stage seconds at identical
map/sequence keys; aggregate last-arriver histograms alone are insufficient.
The existing results fields are documented in wait-profile-request.md. Validate
actual owned rows against nonlocal.projector_gemm previous_stage_shape[0] and
confirm M2_STICKS in the manifest. Placement/NUMA needs affinity and memory
placement evidence; neither can be inferred from rank number. Si32 results are
pending. No numerical equivalence run is needed for this documentation-only
analysis: runtime sources are unchanged, and no Mac timings were taken.

## Local validation

- Arithmetic ownership audit: exit 0; reproducible output ownership-audit.json.
- Existing wait parser regression suite: exit 0, five tests passed.
- `CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast`:
  exit 1. Build/fmt/clippy passed; 620 tests passed, 156 ignored, three known
  CLI pseudopotential-cache failures because TMPDIR is inside the repository.
  The optional process-load query was also denied by the sandbox. Log:
  target/scratch/arrival-closure-check.log. No performance inference from checks.
- `git diff --check`: exit 0. Changes remain uncommitted.
