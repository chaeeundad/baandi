# NONLOCAL Phase 2 diagnosis — recommend closing

2026-11-20. Analysis only; no optimization implemented. The latest user task
supersedes the older charter's Mac profiling instruction. All performance inputs
below are **relayed dev08 data**, supplied by the orchestrator for `bad02ee`,
Si64 MPI np8. Local source inspected at
`b7889a7ee461da7784a05dfef7e2236aae2ee03e`. No Mac profiling was performed.
The raw dev08 results, binary/input hashes and exact per-call shapes were not
attached to this relay; these numbers are not independently reproduced here.

## Evidence and denominator

Lane attribution is 0.7825 s / 9.6702 s = 8.09%. Repeats are
0.7825, 0.8125, 0.7486 s (reported spread 8.5%). Timed-1 and counts-1 have
24 identical work keys. These are three profile repetitions, not an A/B gain
or a statistical significance test.

| Stage | Rank-max s | Rank-min s | Calls | Lane share | Asymmetry |
| --- | ---: | ---: | ---: | ---: | ---: |
| Projection | 0.3054 | 0.2584 | 103 | 39.0% | 15.4% |
| Back product | 0.2681 | 0.2178 | 84 | 34.3% | 18.8% |
| Reduction | 0.1178 | 0.0394 | 103 | 15.1% | 66.6% |
| Coupling | 0.0474 | 0.0362 | 84 | 6.1% | 23.6% |
| Image add | 0.0318 | 0.0231 | 84 | 4.1% | 27.2% |
| Audit | 0.0120 | not supplied | 19 | 1.5% | not supplied |

Fine accounts for 0.5425 s over 11 maps (49.32 ms/map), coarse 0.1265 s
over 6 maps (21.08 ms/map), verify_fail 0.0569 s and verify_pass 0.0566 s.
The tiny sum differences are relay rounding. These are sums of bucket/stage
rank maxima: different ranks can supply different maxima. They are attribution
budgets, not a reconstructed critical-path trace. Likewise rates below divide
relayed work by relayed maximum time, not matched per-rank hardware counters.
Do not multiply either work or time by eight in these rate calculations.

The separate H-apply relay is 2.78 s versus QE h_psi 1.91 s, a 0.87 s gap
in the same MPI model; nonlocal is 0.725 s and FFT/local 2.057 s. This H-only
budget differs from the lane's 0.7825 s including energy audits. It does not
establish how much of the gap belongs to nonlocal. No new QE comparison is made.

## Arithmetic and bandwidth floors

Use decimal GFLOP/MB/GB; the supplied single-thread DGEMM reference is
64–108 GFLOP/s and STREAM is 12.97 GB/s. Define the optimistic reference
time as max(F / DGEMM_rate, logical_bytes / STREAM_rate). These are reference
roofline estimates, not guaranteed bounds on a contended MPI rank: STREAM is
not a measurement of these cache-resident panels or eight concurrent ranks.

| Quantity | Projection | Back product |
| --- | ---: | ---: |
| Work, GFLOP | 10.43 | 7.48 |
| Logical traffic, MB | 679.6 | 536.0 |
| Observed aggregate rate, GFLOP/s | 34.15 | 27.90 |
| Mean stage ms/call | 2.965 | 3.192 |
| Logical MB/call | 6.598 | 6.381 |
| Logical arithmetic intensity, FLOP/byte | 15.35 | 13.96 |
| STREAM reference time, s | 0.05240 | 0.04133 |
| STREAM reference ms/call | 0.509 | 0.492 |
| Arithmetic reference at 64 GFLOP/s, s | 0.16297 | 0.11688 |
| Arithmetic reference at 108 GFLOP/s, s | 0.09657 | 0.06926 |

Thus both products together cost 0.5735 s for 17.91 GFLOP, or 31.23 GFLOP/s.
Reaching 64–108 GFLOP/s for the **entire wrappers** would save 0.29366–0.40767 s.
This is theoretical headroom, not an expected gain from a specific change.
A 0.15 s saving solely here requires 42.29 GFLOP/s end-to-end: 26.2% less
stage time or 35.4% greater throughput. Even this is only a necessary lane
budget condition, not proof of 0.15 s whole-run improvement.

Reduction carries 17.3 MB / 103 = 0.168 MB/call. Bytes/STREAM is 0.00133 s
total; STREAM is not a network model. Measured reduction is 0.1178 s, mostly
rank waiting according to the relay, with 66.6% asymmetry. PACKED_REDUCE's
28% payload reduction gave no wall change. No volume-based speedup is assumed.

## Actual algorithm and call-shape model

Read `distributed_real_gamma.rs::projections_reduced` and its H-apply back
product; `real.rs::{multiply_slices_impl,gram_dgemm,dgemm_panels,assemble}`.
Let N be locally owned real rows, P projectors, B active bands, L=256.

* Projection is P^T X: logical output (P,B), reduction N. The MKL wrapper
  makes ceil(N/L) calls (m,n,k)=(P,B,min(L,remaining N)), transposed A,
  lda=ldb=N. Separate P-by-B partials are assembled in fixed ascending order.
* Back product is P W: logical output (N,B), reduction P. Calls are
  (min(L,remaining N),B,P), lda=N, ldb=P. Each writes a contiguous private
  tile, then copies it into the strided destination. The wrapper scans the
  result for nonfinite values. The image is added to local/kinetic H afterward.
* Projector and vector storage is already column-major f64. Complex
  representatives are packed into orthonormal real coordinates once during
  operator construction, not on every product. There is no per-call complex
  GEMM or complex-to-real conversion to eliminate. MKL may internally pack
  strided panels; its cost is not separately timed by Phase 0.
* Distributed dispatch is controlled by M2_DGEMM, not MKL_NONLOCAL. The F11
  dispatch note already identifies these products as MKL / MKL-wide eligible.
  There is no evidence of a current silent fallback needing repair.

For an **illustrative model**, N=ceil(23871/8)=2984, P=256 and B=132 gives
12 calls per product, eleven full (256,132,256) shapes and a tail of 168
reduction/output rows. Each full call is 17.30 MFLOP, with about 1.065 MB of
logical operands/output; the local projector table alone is 6.111 MB. The
global n_pw and band count come from the F11 Si64 table. Equal row ownership
and P=256 are modeling assumptions, not a recovered bad02ee call trace.
Coarse grids, uneven partitions and changing active widths alter this model.
74/103 observed projections have B>=33; only 29 have narrower widths, and
their FLOP/time fraction is unknown. Call fraction is not work fraction.

The reference large-DGEMM rate therefore need not apply: fixed partial
products, short tails, transposed strided reads, repeated internal packing,
assembly, output copying, finite scans and MPI memory/cache contention all
contribute to lower wrapper throughput. Coupling and add account for another
0.0792 s of lane time but are **outside** projection/back timers, so cannot
explain the low rates of those two timers. Attribution among the internal
causes is not established by this profile. The two failed shape-changing
DAVIDSON candidates in the user's relay are further evidence against assuming
that fewer/larger/smaller calls will improve this CPU.

## Candidate assessment (no implementation)

| Candidate | Quantified opportunity | Decision |
| --- | --- | --- |
| Same-shape direct strided back output, fuse image addition | Eliminate at most the measured 0.0318 s add; tile read/write copy model is 16NB = 8F_back/P bytes. At illustrative P=256 this is 233.75 MB, or 0.01802 s at STREAM. Combined optimistic planning estimate about 0.050 s, before finite checking/adapter costs. | Below 0.15 s. beta=1 accumulation also changes rounding; direct output needs a stride-aware adapter. No code. |
| Precompute projector-times-coupling, or improve the small coupling loop without changing GEMM shapes | Even zero-cost coupling saves only 0.0474 s (6.1%). Preweighted projectors still require two GEMMs, extra storage/setup, and can change association. Coupling plus add elimination ceiling is 0.0792 s. | Below 0.15 s. No material-specific diagonal-D assumption. |
| Persistent MKL packing / reduce projection assembly traffic, preserve call dimensions | Full-panel projection tile write+read is approximately 16F_projection/(2L) = 325.94 MB; including destination read/write for assembly gives about 651.88 MB, or 0.05026 s at STREAM, plus tail effects. Packing is opaque and unmeasured; it could matter, but no evidence attributes >=0.15 s to removable packing. | Expected gain cannot be established. Do not convert total arithmetic headroom into a packing speedup. No new packed ABI or shape experiment. |
| Smaller reduction payload / narrow bands | Complete reduction elimination is only 0.1178 s of attributed time; payload cut already failed. Narrow calls are the minority and unweighted. | Do not repeat PACKED_REDUCE, tile narrowing or column groups. |

Copy/assembly estimates are bandwidth-scale estimates, **not hard upper
bounds**: poor loops or contention could exceed them. They also cannot simply
be added to the roofline gap or treated as independently achievable. Likewise
prepacking benefits cannot be ruled out. What is missing is affirmative
evidence for a permitted lever with expected whole-run gain >=0.15 s.

## Ledger reconciliation and disposition

A43 row traversal regressed dev08; A58b/P15B band sharing was slower in all
60 paired rounds despite 3.2x fewer projector reads. Neither is reopened.
P48C's BLAS-3 complex k-point change addressed a different path, already
implemented with a shape gate. B07/B08's historical scalar-loop estimates
do not describe this real-Gamma DGEMM path. F11's fallback opportunity was
classic solver projection, not these nonlocal products.

NONLOCAL-DGEMM-AB initially kept MKL_NONLOCAL opt-in; later ADOPT-DEFAULTS-1
records adoption and this checkout has default `on`. The current user's
instruction forbids promoting it without new evidence. We make **no switch
or default edit**, do not use a toggle as a new optimization, and do not infer
the current MPI route from that serial switch. The historical discrepancy
does not authorize a rollback either.

**Recommend closing LANE-NONLOCAL for this optimization round.** No candidate
has a supported expected saving >=0.15 s at Si64 np8. This is an evidence and
priority decision, not a proof that the lane is optimal. There is no new
switch, Mac numerical run, or dev08 A/B request. Existing Phase 0 evidence
remains historical. Reopen only with same-shape dev08 attribution showing a
removable >=0.15 s cost (or a new explicit task threshold); then pre-register
paired A/A and crossed A/B, |mean|/SE >=3, numerical |dE|<=1e-11 Ha and
unchanged convergence/certificate gates. No other parallel model is inferred
from Si64 np8.
