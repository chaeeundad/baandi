# Phase 2 follow-up: fixed GEMM shapes, packed MPI border payloads

Current decision: UPPER_BORDER and ROTATE_COLUMNS rejected; both remain OFF.
The latest reassessment at the end selects PACKED_REDUCE. Earlier candidate
forecasts below are historical and superseded by the two dev08 decisions.

## Evidence, units and attribution

Measurements are the user's dev08 relay: source **2733278**, instrumented
`on-1`, CPUs 0–7 on socket 0 for MPI8. Phase 0 integrity was 15/15 pairs:
abs(dE)=0, identical science/history/verification/density, no dropped records.
No Mac timing supports this diagnosis. UPPER_BORDER was rejected by the orchestrator for a resolved dev08 loss;
its switch is retained default OFF. Follow-up implementation base: 2bde16a.

The earlier 104.1% exclusive split is **withdrawn**, not evidence of an engine
bug. `phase_profile.buckets[].max` mixes the parent `solve`, children, residual
keys and `.calls` counts. Removing the specified duplicate/non-leaf gives:
MPI8 Si64 100.9% (drop solve.solver_total), one-core Si64 99.1% (drop
solve.other_engine, identical to solve.other), MPI8 Si32 98.9%, one-core Si32
99.1%, MPI8 Si16 98.5%. Remaining closure error is not assigned to this lane.
Serial Si64 lane solve-call envelope is 27.242 s; the phase-profile parent is
32.7484 s. These are different sources/scopes; **do not interchange their
denominators**, add them or claim their difference is a measured algebra cost.

| Model | Cell | Earlier external s | Earlier lane envelope s | Envelope/wall |
|---|---|---:|---:|---:|
| MPI8 | Si64 | 9.7715 | 6.4323 | 65.8% |
| MPI8 | Si32 | 1.9254 | 1.1224 | 58.3% |
| MPI8 | Si16 | 0.6625 | 0.2996 | 45.2% |
| One core | Si64 | 45.5618 | 27.2416 | 59.8% |
| One core | Si32 | 8.1225 | 4.5788 | 56.4% |

These earlier wall/envelope summaries are budget context, not the on-1
per-site measurements below. Models are never compared by a speedup ratio.
Rank max/min was 1.004/1.005/1.002 at Si64/32/16: imbalance is not the lever.
H application and its operator reductions belong to FFT/NONLOCAL; no candidate
here changes H, its products, stopping tolerances or the SCF schedule.

## Actual on-1 counts and times

MPI8's 10,176 records aggregate **all eight ranks**; seconds below are CPU/rank
span sums, not external wall. Dividing by 8 gives a mean rank's accumulated
span time, NOT a critical-rank maximum. Calls are rank-local invocations;
808 EVD records mean 101 replicated solves per rank, not 808 sequential solves.
Every timer is inclusive: eigen is inside eigen_step, tree_sum is inside its
site, and serial residual is inside rotate_residual. Never sum such rows.

| MPI8 Si64 operation | All-rank calls | All-rank s | Mean-rank s | Mean us/call | Dominant (n,m,k) |
|---|---:|---:|---:|---:|---|
| reduced_gram | 1464 | 6.996 | 0.87450 | 4778.8 | 2986,264,132 |
| rotate_residual | 808 | 6.715 | 0.83938 | 8310.4 | 2986,132,132 |
| reduced_eigen_step | 808 | 5.637 | 0.70463 | 6975.9 | 132,0,132 |
| reduced_eigen | 808 | 5.636 | 0.70450 | 6975.3 | 264,264,264 |
| tree_sum | 2576 | 2.325 | 0.29063 | 902.6 | payloads below |
| reduced_transform | 808 | 0.962 | 0.12025 | — | — |
| factor_extension | 896 | 0.886 | 0.11075 | — | — |
| residual | 808 | 0.837 | 0.10463 | 1035.9 | 2986,0,132 |
| return_gram | 152 | 0.388 | 0.04850 | 2552.6 | 376,132,16 |
| ritz_coefficients | 808 | 0.355 | 0.04438 | — | — |
| restart | 88 | 0.183 | 0.02288 | — | — |

| One-core Si64 operation (1025 records total) | Calls | s | Mean us/call | Dominant (n,m,k) |
|---|---:|---:|---:|---|
| rotate_residual | 109 | 5.877 | 53918.5 | 23871,132,132 |
| reduced_gram | 109 | 2.182 | 20015.6 | 23871,132,132 |
| overlap_border | 109 | 2.123 | 19473.7 | 23871,132,132 |
| reduced_eigen_step | 109 | 1.275 | 11700.3 | — |
| reduced_eigen | 102 | 1.162 | 11388.4 | 132,132,132 |
| return_gram | 21 | 0.765 | 36419.6 | 403,132,132 |
| residual | 109 | 0.504 | 4623.2 | — |
| orthogonalize | 2 | 0.174 | — | — |
| bordered_update | 7 | 0.086 | — | 480,479,480 |
| reduced_transform | 109 | 0.133 | — | — |
| factor_extension | 109 | 0.127 | — | — |
| ritz_coefficients | 109 | 0.045 | — | — |

Serial overlap is a separate S border; MPI's reduced_gram contains S/H border
calls. B5 bordered updates occur only on serial in this relay. Do not treat the
missing MPI overlap row as missing S work.

## Arithmetic and bandwidth models

The new relay supplies dominant shapes and means over ALL calls, not the
joint shape/time distribution or sum(n*m*k), sum(m^3). Thus the following is
a **dominant-shape model against the all-call mean**, not an exact per-shape
measurement, an aggregate FLOP total, or a guaranteed removable gap. This
limitation remains even though a numeric useful screening model is now possible.

Measured dev08 MKL single-thread DGEMM rates, GFLOP/s: (23871,264,132)=77.14;
(23871,132,132)=63.96; (2982,264,132)=98.03; 264 cubed=95.40;
396 cubed=94.24. STREAM triad one core=12.97 GB/s. These sustained rates are
**empirical reference rates**, not hardware peak ceilings. GEMM transpose,
packing/layout and eight simultaneously busy cores can alter them. In
particular 98.03 is a nearby 2982-row proxy, not a measured 2986x132x132 rate.

For two rotations (n x m)(m x k), F=4nmk. For one TN Gram border
(m x n)(n x k), F=2nmk. MPI combined S/H trace k=ks+kh, so that formula
already counts both widths; never multiply it by two again. Minimal distinct
array bytes: rotation 8(2nm+mk+2nk); H border 8(nm+nk+mk); S border
8(nm+mk). Residual passes, scratch copies, cache misses, packing and collective
traffic are extra. `max(F/rate, bytes/bandwidth)` is the optimistic comparison
model. A STREAM model of small cache-resident matrices is not a DRAM lower
bound; eight-rank aggregate memory bandwidth was NOT measured here.

### One core, Si64

| Operation/model shape | F GFLOP | Minimum MB | Compute model ms | STREAM model ms | Observed mean ms / compute model |
|---|---:|---:|---:|---:|---:|
| rotation pair, 23871x132x132 | 1.66371 | 100.9705 | 26.0118 | 7.7849 | 53.9185 / 26.0118 = 2.073 |
| H Gram, same shape | 0.83186 | 50.5549 | 13.0059 | 3.8978 | 20.0156 / 13.0059 = 1.539 |
| S Gram, same shape | 0.83186 | 25.3472 | 13.0059 | 1.9543 | 19.4737 / 13.0059 = 1.497 |

Correction to the suggested interpretation: rotation's 0.83186 GFLOP and
13 ms are **one** GEMM, whereas its timer contains **two** plus residuals.
So 53.9/13 is not a fourfold inefficiency. Residual accounts for 0.504/109=
4.624 ms on average, already within rotation; subtracting only that child
leaves approximately 49.295 ms versus the 26.012 ms pair model (1.895x).
Even this is not a measured shape-matched gap. If all calls had the dominant
shape, the pair model would total 2.8353 s versus 5.877 s; that counterfactual
3.0417 s is NOT a claimed saving. At m=264 the measured 77.14 rate gives a
43.135 ms pair model, illustrating how larger-order calls can explain much
of the apparent gap without any kernel defect.

### MPI8, Si64 (each rank separately)

| Operation/model shape | F GFLOP/rank | Minimum MB/rank | Compute model ms | One-core STREAM proxy ms | Observed mean ms / model |
|---|---:|---:|---:|---:|---:|
| rotation pair, 2986x132x132 | 0.20811 | 12.7523 | 2.1229 | 0.9832 | 8.3104 / 2.1229 = 3.915 |
| Gram border, 2986x264x132 | 0.20811 | 6.5852–9.7384 | 2.1229 | 0.5077–0.7508 | 4.7788 / 2.1229 = 2.251 |

The Gram byte range is pure S versus pure H; combined S/H lies between them.
MPI rotation's timer excludes the separately timed residual. MPI Gram includes
its reduction and result fill/symmetrization. At m=264 a rotation pair instead
models 4.2459 ms; no claim that all 808 calls cost the dominant-shape FLOP.
The all-rank measured spans 6.715/6.996 s cannot be compared to a one-rank
model, nor treated as external-wall savings. Per-rank model budgets if all
calls had these shapes: rotation 101*2.1229=0.2144 s; Gram
183*2.1229=0.3885 s. These are clearly labelled counterfactuals, not exact totals.

### Reduced eigensolve: LAPACK, not a GEMM

The numerical adapter calls MKL `dsyevd` for all eigenpairs. LAPACK's driver
uses tridiagonalization, divide-and-conquer tridiagonal eigenvectors and a
back-transform; it is not one cubic GEMM. Its vector-workspace requirement is
1+6m+2m^2 doubles. See [LAPACK DSYEVD documentation](https://netlib.org/lapack/explore-html/d8/d30/group__heevd_ga25b71a69f9921df0a6050aa5883f54f4.html).

Use the requested **9m^3 leading-work class** as a screening approximation,
not an exact dsyevd flop count (deflation/spectrum/backend matter). Conventional
Householder reduction alone has leading work (4/3)m^3 and BLAS2 portions; it
cannot all run at DGEMM rate. Input+output compulsory bytes alone are 16m^2,
which excludes repeated scans, work arrays and reductions. Those bytes give
only a very weak floor. A realistic LAPACK timing floor requires measured
same-order dsyevd/tridiagonalization stages, not 9m^3 divided by peak BLAS alone.

| Order m | Rate proxy GFLOP/s | One GEMM ms (2m^3) | EVD compute model ms (9m^3) | Compulsory MB / STREAM ms |
|---|---:|---:|---:|---:|
| 132 | 95.40 (264 proxy) | 0.04822 | 0.21698 | 0.2788 / 0.02150 |
| 264 | 95.40 | 0.38574 | 1.73582 | 1.1151 / 0.08598 |
| 396 | 94.24 | 1.31789 | 5.93052 | 2.5091 / 0.19346 |

Correction: a **264^3 GEMM is 0.386 ms, not 0.04 ms** at the supplied rate;
0.048 ms corresponds to order 132. MPI average 6.9753 ms / order-264 EVD
model = **4.018x**, not ~175x. At order 396 that same average/model ratio is
1.176x. Serial average 11.3884 ms / order-132 proxy = 52.49x, but it is an
average over 102 different-order EVDs, not the time for the dominant order.
It cannot establish a 52x kernel gap. The replicated MPI EVD is not divided
by 8: reducing PW rows does not reduce this order-m work. At the earlier
Si32/Si16 MPI shares (14.2%/20.4% of solve), even a genuine 25% EVD saving
would represent only ~2.07%/~2.31% of their respective external wall; this
is an Amdahl scenario, not grounds to reopen failed sliced-EVD algorithms.

### MPI reductions, sizes and latency

2576 tree sums across all ranks = 322 calls/rank. Site totals / per-rank:
setup 152/19, apply 672/84, reduced_gram 808/101, lean_border 656/82,
residual_norms 152/19, restart_defect 136/17. Excluding apply leaves 238
non-apply calls/rank; the two border sites total 183/rank, exactly matching
the 1464 Gram calls /8. The two reported payload bins account for:

* 280,896 B *304 rank calls =85.392384 MB (38 calls/rank);
* 270,336 B *200 rank calls =54.067200 MB (25 calls/rank);
* combined 139.459584 MB across ranks, **17.432448 MB/rank**, only 504/2576
  calls. Other payloads remain unknown; this is a partial total, not total
  communication. At 12.97 GB/s merely reading this partial payload once is
  1.344 ms/rank, not a network/collective floor.

Current lean-border payload is 8*m*(ks+kh+1), including the residual-norm
tail; 280896=8*264*133 is consistent with m=264,k=132. For a power-of-two
recursive-doubling sum the model is log2(P)*(L+payload/Bnet), plus local sums
and memory traffic; with P=8 this has three stages. Actual network bandwidth
and latency are not supplied, so no numerical network-floor claim is made.
Mean tree_sum span 0.290625 s/rank includes apply and wait/arrival skew;
it is already inside callers. Packing only retained upper entries is a
possible later bandwidth experiment, but would not remove these 322 latency
payments. No packed-reduction experiment is proposed in this response.

## Candidate reassessment

1. **UPPER_BORDER: REJECTED, default OFF.** The dev08 counter pair removed
   16.38% of traced FLOP (170.792e9 ->142.821e9 across ranks), with exact
   per-rank prediction and bit-identical science, yet the crossed solve
   screen slowed by +0.005679 s (~0.57%), abs(mean)/SE=4.18, all three pairs
   positive. Calls increased 1746 ->1968 (+12.7%). A/A mean -0.000216 s,
   floor0.0029167 s, ratio0.14. The orchestrator's decision is final; no
   repeat with a different narrower tile size. The -17% measured DGEMM
   efficiency for the smaller comparison shape and increased calls invalidate
   the old FLOP-proportional gain scenario. See ledger row
   DAVIDSON-UPPER-BORDER-REJECTED for every supplied number.
2. **SELECTED ALTERNATIVE: serial ROTATE_COLUMNS, default OFF.** The current MKL pair makes
   2*ceil(n/256) GEMMs into scratch and copies them to the output. Reuse the
   numerics layer's existing fixed 32-column `multiply_column_groups` adapter
   twice: each MKL call writes complete contiguous output columns directly.
   At n=23871,k=132: **188 -> 10 MKL calls**, same 4nmk FLOP, and remove
   4nk*8=100.831 MB of scratch write/read traffic per pair (direct output writes
   remain). This is not the rejected fused-residual kernel:
   H X is still produced and kept; no residual policy changes. Candidate
   needs no new FFI, dependency or allocation buffer. Existing fused paths
   take priority and are deliberately outside this candidate.
   Conditional serial Si64 target: 10–25% of the non-residual rotation
   span (5.877-0.504=5.373 s), **0.537–1.343 s** (~1.18–2.95% of earlier
   wall). This is a hypothesis, not a floor-derived promise. Repeated V/HV
   reads across five column groups and tall-skinny GEMM efficiency can erase
   it; only dev08 A/B decides. The supplied rates are for width132, not32.
   Serial Si32 has the same structural hypothesis but no supplied rotation
   seconds for a numeric wall estimate. The distributed rotation code is
   unchanged: expected MPI8 improvement **zero**, used as a control.
3. **Triangular reduced-coordinate products: defer.** Their current budget
   is smaller: transform+coefficient spans 0.178 s serial and mean-rank
   0.164625 s MPI8. A conditional 10–30% reduction yields 0.018–0.053 s
   serial (0.04–0.12% wall), 0.016–0.049 s MPI (0.17–0.51%). These are
   separate operation rows, not eigen_step minus eigen (the earlier proxy
   conflated scopes). No need to implement before the two stronger probes.

The new serial route has a larger controllable call/data-movement reduction
than the rejected arithmetic-pruning attempt, without changing the mathematics.
MPI's EVD row is large, but a corrected LAPACK model does not reveal a new
safe algorithm. REL-SLICES-CLOSED still applies: sliced MRRR failed the
orthogonality gate, inverse-iteration slices cost 1.42–1.64x. Advance,
restart-keep, PR-CG and async alltoall remain rejected; no new premise reopens
them. No residual checks are omitted and no basis capacity is retuned.

## Next measurement and stopping point

Extended DAVIDSON_PROFILE labels each serial rotation pair with shape and
route, enabling an exact structural MKL-call/scatter-byte model per map.
The dev08 request isolates ROTATE_COLUMNS with UPPER_BORDER forced OFF, marks
the old request closed, and pre-registers A/A floor, crossed pairs and
abs(mean)/SE >=3. For a stronger aggregate roofline claim, return shape-grouped
count/FLOP/time rather than applying a mode to every call; existing records
and parser already retain it. Missing joint histograms do not block this
explicitly labelled model or the candidate screen. No adoption yet.

## Rejection response: account for rotation beyond one GEMM

Use unambiguous dimensions: N=PW rows, M=basis columns, K=rotated roots.
Each rotation is (N x M)(M x K). Default MKL calls are (256 x M)(M x K)
into scratch (last row panel shorter); the alternative calls are
(N x M)(M x 32) directly into independent output columns (last group shorter).
It keeps **both N and the entire inner M**, unlike pruning row prefixes in
an upper Gram rectangle. In BLAS notation it changes independent output-N
partitioning, not the reduction-K or the mathematical operation. No product
entry is discarded; 4NMK FLOP are unchanged. This matches the suggested
independent-output split, not another attempt to save arithmetic by narrowing
upper-border tiles. Width32 efficiency is still a risk to be tested, not an
assumption of a guaranteed gain.

Source-level pass inventory for the normal serial non-fused call:

| Work | Arrays / operations | Representative logical bytes (N=23871, K=132) | Candidate |
|---|---|---:|---|
| two GEMMs | X=VC, HX=(HV)C | 4NMK FLOP; two output blocks | unchanged math |
| row scratch materialization | write scratch, read scratch, write X/HX | extra 32NK =100.831104 MB beyond direct output | remove extra scratch traffic |
| finite-output check | read X and HX after GEMMs | 16NK =50.415552 MB | retained |
| residual formation | read X/HX, write R=HX-lambda*X | 24NK =75.623328 MB | retained |
| residual norm | one R scan, four independent sum-of-squares lanes, fixed final sum/sqrt | 8NK =25.207776 MB | retained |
| exceptional norm | repeated hypot pass only for nonfinite/tiny accumulated squares | data-dependent, no rate assigned | retained |
| finite-norm check | read K norms | 8K =1056 B | retained |

Logical bytes are access counts, not measured DRAM traffic. Cache residency
makes the STREAM conversion only a model. Default pair finite checks are in
the pair adapter; the column adapter checks each complete output separately,
covering the same entries. There is no implicit missing third GEMM, H apply
or fresh Ritz check inside this span. C gathering/transformation happens
before it. Row-strided input panels may trigger internal MKL packing; that
cost is not separately measured and is not equated to the scratch-copy model.

Measured residual child is 0.504 s /109 =4.624 ms/call, leaving about49.295 ms
for pair GEMMs, layout/copies, finite scans and glue from the53.9185 ms mean.
At representative M=132, the pair compute model is26.012 ms, so the apparent
non-residual excess is23.283 ms. At M=264, it is only about6.160 ms versus
43.135 ms. Because the observed mean mixes shapes, neither difference is a
measured recoverable budget. The model's extra scratch bytes alone correspond
to7.774 ms at one-core STREAM12.97 GB/s; cache/packing can alter this strongly.
A fused norm-only alternative would touch at most the already measured0.504 s
residual row (~1.1% of prior wall even if free), whereas the selected route
acts on the larger5.373 s pair/glue span. Do not repeat rejected RITZ_FUSED.

## Full EVD frequency per map: no evidence of duplicate full solves

The MPI relay has808 reduced_eigen_step and808 reduced_eigen records across
8 ranks: **101 full EVDs/rank, one full EVD per reduced step**, not808 EVDs on
one critical rank. Serial has109 steps,102 full EVDs and7 bordered attempts,
consistent with B5 replacing seven full solves. Eigen and eigen_step are
nested, not two separate EVDs. Code invokes solve_reduced once after each
new reduced problem; it does not solve again for each rotated root or each
residual norm. Thus the aggregate counts do not expose redundant repeated
full solves that can simply be deleted.

Do not divide101 by11 fine SCF iterations: profiles include coarse,
verification and multiple solve contexts. Setup reductions/return-Gram counts
are also not a reliable count of maps. The original per-map joint records
were not attached in the summary. The parser now emits `eigen_map_work`,
keyed by rank/map/stage/kind/SCF iteration, with solve_ids, reduced_steps,
full_evd_calls, bordered_attempts, full_evd_order_counts, full_evd_sum_m3 and
full_evd_per_solve. This gives exact per-map frequencies from returned files
without introducing more numerical instrumentation. Keep these invariant in
the alternative A/B. Same order on two steps does not establish the same
matrix or justify reusing eigenvectors. Delaying RR changes the algorithm and
is not warranted by these counts or the failed sliced-EVD evidence.

## Reassessment after both shape candidates lost (2026-11-20)

The ROTATE_COLUMNS oracle is decisive: Si64 serial +6.10% (19.20 SE),
Si32 +4.19% (30.46 SE). All six pairs slower. A 256-row pair costs15.50ms,
whereas every tested row-height or column-width alternative costs more.
Default row panels deliver108.52 GFLOP/s versus53.86 for tall32-column groups.
The earlier26.012ms tall-DGEMM floor is therefore not a floor for the actual
panel implementation. It cannot establish23ms of removable non-GEMM work.
The removed51.3MB copy-volume model is only~4ms at STREAM, outweighed by
~22ms GEMM shape loss. Keep both failed switches off and preserve GEMM shapes.

Two remaining candidates were assessed:

1. **Pack distributed border reduction payloads (selected).** Tree sums total
   2.325s over2576 all-rank calls,0.9026ms/call; rank-average0.290625s, not
   critical-rank wall. Border sites account for808+656=1464 calls across ranks
   (183/rank); all sites include672 H-apply reductions belonging to other lanes.
   Full rectangles contain lower entries never consumed. Packing after full
   GEMMs preserves their shape, calls, FLOPs, outputs, and all rank addition
   dependencies. At m264,k132,no pending norms: full(m*(k+1))=35112 doubles;
   upper rows133..264=26202 doubles;280896->209616B (-25.375%). At fresh m132,
   k132:140448->70224B (-50%). Multiple S/H borders and the exact norm tail are
   handled generically; narrow borders need less triangular pruning.
   In-place compaction adds at most one read/write of the retained payload,
   ~0.419MB/representative call (~0.0323ms at the supplied one-core12.97GB/s).
   This is a bandwidth model, not MPI memory bandwidth or measured overhead.
   Latency is unchanged. A conditional10–25% reduction of affected collective
   time is plausible only if bandwidth matters. Even25% of *all* tree time
   is merely0.07266s/rank (~0.74% of prior9.7715s external wall); this is an
   optimistic scale, not an expected measured whole-run gain. Border-only
   timing and exact payloads determine the actual opportunity. Serial gain=0.
2. **Fuse residual formation and four-lane norms (not selected).** Serial
   residual child is0.504s total (~1.1% of45.56s wall), including formation,
   norms and validation. Even deleting all of it cannot recover the old
   apparent rotation gap. Fusion removes one read but risks vectorization
   and arithmetic-order changes; the residual kernel already uses fixed
   groups. MPI residual is0.837/8=0.1046s rank-average, with overlapping parent
   time. A small fraction of those budgets is a hypothesis, not evidence.
   No measured norm-only span identifies a stronger opportunity than payload.

Full EVD remains one per changed reduced problem (101/rank on MPI; serial102
full plus7 bordered steps). The9*m^3 proxy at264 and95.4GF/s is1.736ms versus
6.975ms mean, but GEMM peak is not a LAPACK syevd measurement: tridiagonal
reduction has a different mix, the mean mixes orders, and replicated ranks
share cache/frequency resources. REL-SLICES already tested alternative drivers
and failed the combined speed/orthogonality gate. No identical unchanged
matrix is solved twice; reusing a factorization of a different reduced matrix
is not equivalent. There is no new evidence to reopen that lever, change
restart capacity, or relax tolerance/map work. H application is untouched.

The new `border_payload` span records n=order,m=S width,k=H width,detail=norm
count. It encloses pack/reduce/unpack and nests one `tree_sum` with actual sent
doubles. Profiles remain bounded and disabled without clocks/allocations;
existing exclusive and inclusive timing scopes must still not be added.
