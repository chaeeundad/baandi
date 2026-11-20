# LANE-FFT Phase 2 diagnosis

Status: candidate selected for dev08 A/B, not adopted. Source profile relayed by
the orchestrator on 2026-11-20: bad02ee, si64 MPI8, CPUs 0–7 socket 0,
72^3 grid, 23871 plane waves, 132 bands. Binary hash/raw results not yet supplied.
This note uses only that dev08 relay and source inspection; no Mac performance
measurement. The follow-up relay supplies all five points below, with zero
count-only/timed differences at every point. No QE comparison is made.

## Integrity and denominators

Timed/count-only controls compare 1152 stage keys with zero count differences;
energy -253.37816772116847 Ha, SCF 11, H 3504, verified in all relayed timed runs.
Rank-0 lane 2.2440 s; rank maximum 2.2649 s (r3), minimum 2.1767 s (r7).
Run wall 8.4265 s. The corrected coverage definition is **lane rank maximum /
run wall**: 2.2649/8.4265 = 26.88%. Rank-0 lane/wall is separately 26.63%;
the discrepancy is resolved, not a run-wise averaging effect. Group totals
below sum to rank-0 2.2440 s. The supplied top-stage
forward/inverse alltoall 0.2677/0.1504 s sum to 0.4181 s, not group 0.3477 s:
these appear to be stage rank maxima, not both rank 0. Raw per-run/rank rows
are required before subtracting stage times. Rates below are relayed rounded
rates, not exact reconstruction of missing byte counters.

| Rank-0 group | Seconds | Lane share | Diagnosis |
|---|---:|---:|---|
| FFT axes | 0.9400 | 41.9% | Existing batched, contiguous lanes; cached plans |
| Permutation/zeroing | 0.7184 | 32.0% | Several necessary layout changes, some redundant writes |
| All-to-all | 0.3477 | 15.5% | Transfer plus arrival/progress/scheduling effects |
| Multiply | 0.1256 | 5.6% | 51.5 GB/s logical, cache-level traffic |
| Kept copy | 0.1123 | 5.0% | Separate copy immediately followed by multiply: removable |

Fine 1.6096 s / 11 maps = 146.3 ms/map; coarse 0.2999 / 6 = 50.0 ms/map;
verify_pass 0.1665 / 1; verify_fail 0.1679 / 1. A verification map is
13.8–14.8% heavier than an average ordinary fine map. Work/kept fraction may
differ; not evidence that its individual FFT kernels are slower. Init,
fine_redo, materialize and between have no lane maps. Enclosing fine FFT and
solve.apply counters have different scope/max aggregation: do not add
pair_unsplit 1.3168, kept_forward 0.3345 or solve.apply 2.4851 to lane time.
Fine h_columns 6025 and h_delta_columns 2244 are not the same counter as H 3504.

## Five-point scale evidence (relayed bad02ee)

Coverage uses lane maximum, whereas group shares use rank 0. The MPI and
one-core tables below are separate models; no cross-model speed ratio is
inferred. All times are instrumented profile times, not candidate A/B results.

| MPI8 case | Run wall s | Lane max s | Coverage | Lane rank spread |
|---|---:|---:|---:|---:|
| si64 | 8.4265 | 2.2649 | 26.88% | 3.89% |
| si32 | 1.9316 | 0.5031 | 26.05% | 0.58% |
| si16 | 0.6726 | 0.1586 | 23.58% | 0.86% |

| One-core case | Run wall s | Lane s | Coverage |
|---|---:|---:|---:|
| si64 | 46.0054 | 11.2682 | 24.49% |
| si32 | 7.9867 | 2.5113 | 31.44% |

| MPI8 case, rank 0 | Axes % | Permutation/zeroing % | Alltoall % | Multiply % | Kept copy % |
|---|---:|---:|---:|---:|---:|
| si64 | 41.9 | 32.0 | 15.5 | 5.6 | 5.0 |
| si32 | 40.1 | 36.0 | 14.2 | 6.6 | 3.1 |
| si16 | 44.9 | 33.4 | 12.5 | 6.1 | 3.1 |

| One-core case | Axes % | Permutation/zeroing % | Alltoall % | Multiply % | Kept copy % |
|---|---:|---:|---:|---:|---:|
| si64 | 52.0 | 34.1 | 0 | 13.9 | 0 |
| si32 | 57.3 | 31.1 | 0 | 11.6 | 0 |

**Primary lever: remove permutation passes and bytes.** Permutation/zeroing
accounts for 31.1–36.0% of the lane at every point. The relay identifies this
group as bandwidth bound under shared socket pressure: MPI8 middle_pack
6.49–7.63 GB/s per rank, roughly 52 GB/s aggregate at the low end, consistent
with socket saturation. The one-core middle_pack observations are
13.72–14.34 GB/s (STREAM triad reference 12.97 GB/s). These are kernel traffic
observations, not a mixed-model wall comparison or hardware bandwidth counters.
Removing traffic is more promising than merely changing traversal order.
Multiply remains cache resident: MPI8 50.4–53.8 GB/s logical, one-core
35.9–40.7 GB/s; do not model it with the DRAM STREAM ceiling.

**Second lever at larger sizes: padded alltoall volume.** MPI8 forward rates
are si64/si32/si16 = 2.38/4.20/4.20 GB/s; inverse = 2.66/2.70/3.06 GB/s.
These use padded collective bytes already counted by the probe. Larger
messages are less efficient here; this is not evidence of a fixed-latency-only
bottleneck. Padding removal still needs the useful-byte counters below.

Rank wait/asymmetry is **scale dependent**: it appears at si64 in the relay,
not as a general MPI8 bottleneck at si32/si16. Their lane spreads are below
1%. The si64 conditional wait analysis below must not be generalized to the
smaller cells or used to promise a general static partitioning gain.

FFT work tracks volume over 1046–60052 points/call (57x), at 1.29–2.15
ns/point, without evidence of a dominant fixed call cost. For example,
forward_fft0 is 1664 points/call at 1.51 ns/point for si16 MPI8; the separate
si64 one-core observation is 60052 points/call at 1.66 ns/point. This supports
retaining existing batching and plan reuse rather than reopening generic FFT
batching/MKL switches. It does not assert mathematically zero call overhead.

## Floors and exact traffic model

Let P=ranks, n0/n1/n2 be chosen FFT axes, S total occupied sticks, s_r local
sticks, Z occupied z planes, l_r local coefficients, p_r owned real planes,
smax=max s_r, pmax=max p_r, b the pair batch width. Complex values are 16 bytes.
Counters count the following **logical** bytes, excluding indices, cache-line
waste, write allocation and backend scratch traffic. Multiply counts 2 real
multiplies/point; gather scaling counts 2 real divisions/coefficient.

| Stage | Bytes per pair (per rank; per plane where stated) |
|---|---|
| scatter | 16 s_r n0 zero fill + 32 l_r coefficient read/write |
| inverse/forward FFT0 | 32 s_r n0 each, minimum lane read/write |
| transpose_pack / transpose_unpack | 32 s_r n0 each |
| middle_unpack, per owned plane | 16 Z n1 clear + 32 S sparse read/write |
| inverse/forward FFT1, per plane | 32 Z n1 each |
| slab_scatter, per plane | 16 n1 n2 clear + 32 Z n1 transpose |
| inverse/forward FFT2, per plane | 32 n1 n2 each |
| multiply, per plane | 40 n1 n2 (complex read/write + real potential read) |
| slab_gather, per plane | 32 Z n1 |
| middle_pack, per plane | 32 S |
| gather_scale | 32 l_r |
| kept_copy, per kept pair/plane | 32 n1 n2 |

A mathematical roofline lower bound is max(F / sustainable compute throughput,
B / sustainable bandwidth at the relevant cache/memory level). Neither
throughput ceiling for these kernels nor an independent socket bandwidth
ceiling is supplied; the follow-up relay estimates saturation from about
52 GB/s aggregate middle-pack logical traffic. The
one-core STREAM triad 12.97 GB/s is a **reference**, not a universal memory
ceiling for eight simultaneous ranks or cache-resident arrays.

For the original si64 MPI8 breakdown, FFT estimate is 5 log2(L) real flops per complex point: 22.925 at L=24,
27.925 at 48, 30.850 at 72. Thus the observed 1.92–2.05 ns/point corresponds
to roughly 11.2–16.1 estimated GFLOP/s over L=24..72 (15–16 applies to L=72,
not every length). Logical arithmetic intensity = 5 log2(L)/32 = 0.716–0.964
flop/byte. A 32-byte STREAM reference costs 2.467 ns/point, which is already
above the observed FFT time: those lane buffers cannot be modeled as full
DRAM streaming at the one-core triad rate. Butterfly loads/stores exceed 32
bytes/point and are mostly cache traffic. For 0.9400 s at the relayed rate,
roughly 459–490 million point-transforms, 14.7–15.7 GB logical traffic; STREAM
reference time 1.13–1.21 s. These ranges are estimates, not recovered counters.
Plan creation is already outside apply; BackendLanes::run uses reusable scratch
and contiguous lane batches. No per-lane plan creation to eliminate.

For si64 MPI8, permutation observed rate divided by STREAM reference: middle_pack 0.50,
transpose_pack 0.71, transpose_unpack 0.86, middle_unpack 1.09,
gather_scale 1.42, slab_gather 1.58, scatter 1.67, slab_scatter 2.46.
For each returned stage with bytes B, reference time is B/12.97e9; measured /
reference is 12.97 / observed_GBps. E.g. middle_pack is 2.00x reference,
transpose_pack 1.40x, transpose_unpack 1.16x. This diagnoses locality/address
work and shared socket bandwidth pressure, not an achievable speedup
guarantee. The five-point evidence favors eliminating bytes over attempting
to make each MPI rank reach the one-core bandwidth reference. Full group byte counters are
missing from the relay, so no invented aggregate permutation floor is given.

Multiply: ~6.47 GB / 0.1256 s, ~161.7 million points, ~323 million flops.
STREAM reference 0.499 s versus observed 0.1256 s demonstrates cache residency.
Kept copy: ~1.359 GB / 0.1123 s, ~42.46 million points; STREAM reference
0.1048 s, close to observed. All estimates use rounded aggregate rates and
must be replaced by exact counter sums when raw JSON arrives.

## Padded exchange and waiting

The code already reports **padded** outgoing alltoall bytes, not useful bytes:
`logical_bytes = 16 * smax * pmax * P * b`. MPI may perform additional internal
copies; that transport traffic is not reported. The claim that padded exchange
volume is larger than this counter is incorrect for bad02ee's instrumentation.
For 72/8, p_r=pmax=9 (no plane-count padding). Inverse useful outgoing bytes on
rank r are 16 b s_r n0; forward useful outgoing bytes are 16 b S p_r.
Across ranks, useful bytes/direction = 16 b S n0; padded bytes/direction =
16 b smax pmax P^2. Utilization = S n0/(smax pmax P^2), or S/(8 smax) here.
23871 coefficients do not determine S or smax: stick support/partition is
needed. Therefore exact padding percentage is not identifiable from the relay.
Counts at the same bucket/rank recover it: inverse useful bytes =
transpose_pack.logical_bytes/2; forward useful bytes = middle_pack.logical_bytes/2.
Compare to each corresponding alltoall.logical_bytes. Include self blocks in
both; off-rank payload additionally excludes rank-specific self traffic.

Relayed outgoing totals per rank: forward 498070272 B, inverse 367670016 B
(sum 865740288 B; eight-rank sum 6925922304 B). Forward is 1.355x inverse
because kept updates need only forward transforms, not a different block size.
At the STREAM reference the per-rank sum is 0.06675 s versus rank-0 0.3477 s
(~5.21x), but this is not transfer time: shared-memory MPI copies, socket
contention, rank arrivals, protocol and synchronization are absent from that
simple reference. At a hypothetical sustained 10 GB/s it is 0.08657 s plus
latency; at 20 GB/s 0.04329 s. Neither bandwidth has been measured here.

The following wait discussion concerns **si64 MPI8 only**. Model each call on rank r as T_r = W_r + alpha(P,algorithm,message) +
B_r/beta(P,message,placement) + local progress/copy costs. The batch widths
vary, so fitting alpha/beta requires per-call sizes and barrier-controlled
microbenchmarks; aggregate calls/bytes cannot separate them. Equal payload
rules out unequal *message volume* as the source of rank timing spread, but
not unequal memory placement/progress/transport costs. If transfer cost and
completion time were identical, max-min would suggest 104.3 ms forward and
39.0 ms inverse additional arrival wait on the longest-waiting rank (39.0%
and 25.9% of their respective maxima). These are conditional proxies, **not
measured wait**, and cannot be added to rank-0 group time. Absolute wait versus
transfer remains unidentifiable (minimum-time rank can still wait).

Before inverse alltoall: coefficient scatter, s_r full n0 FFT lanes, pack.
Partition balances coefficient lengths l_r, not FFT lane counts s_r; mirror
sticks are kept together to avoid Gamma mirror exchanges. Before forward:
p_r planes of equal full-grid FFT2 work, middle sparse work, plus upstream
state/work imbalance. Here all ranks have exactly 9 planes, same potential
point count and S middle entries; rank-independent forward work does not
justify a plane reassignment. Small rank compute differences can accumulate
before the collective; previous collective release, OS scheduling and memory
contention also matter. Static stick reassignment may improve inverse arrivals
but trades off coefficient/nonlocal balance and padding; it cannot be justified
from aggregate rank times alone. Request per-rank FFT0 points/local coefficient
counts from existing scatter/FFT counters before proposing it. A synchronized
transpose microbenchmark can calibrate transfer, but barriers alter the system
and cannot be subtracted as an exact production-wait measurement.

## Candidates and choice (predictions, not acceptance evidence)

1. **Chosen: M2 kept-copy/multiply fusion**, `BAANDI_FFT_KEPT_FUSE=on` (default
   off). Read kept complex value once, perform the same `*= potential`, write
   slab once. 72 -> 40 bytes/kept point; saved bytes = old kept_copy bytes
   (~1.36 GB/rank in this profile), 44.4% of the affected two-pass traffic.
   Existing multiply work/calls stay unchanged; kept_copy becomes zero.
   Conditional expected dev08 saving 0.03–0.10 s/run (~0.4–1.2% of 8.4265 s),
   based on removing a 0.1123 s pass while allowing the fused multiply to cost
   more when it reads the external kept buffer. Do not claim the entire
   0.1123 s is guaranteed or infer wall gain from a relabeled timer.
   At si32/si16 MPI8 kept_copy is only 3.1% of the lane, versus 5.0%
   at si64; expect a smaller affected share, without assigning gains from
   rank-max totals multiplied by rank-0 fractions. The broader permutation
   group remains the larger follow-on opportunity.
   Serial PairFusedPlan already writes source*potential into slab: no serial
   optimization or one-core speedup predicted. No FFT/reduction order change.
2. **Alternative: locality-preserving middle pack traversal**. Precompute
   owner-major sparse copy descriptors and contiguous plane tiles, retaining
   FFT order and alltoall payload. Target middle_pack (relayed max 0.1602 s,
   ~6.49 GB/s); 20–40% stage improvement would be ~0.032–0.064 s on that rank,
   perhaps 0.02–0.06 s critical-path saving (~0.2–0.7% wall). This depends on
   cache/tile costs and must be profiled separately; do not sum rank maxima.
   The five-point socket-saturation evidence weakens this traversal-only
   estimate: without reducing bytes, that predicted gain may not materialize.
   Prioritize pass/byte removal over this alternative in follow-on work.
   Larger redesign could fuse middle layouts with MPI buffers, but middle
   FFTs require dense lanes and MPI needs destination-major blocks: pack and
   unpack straddle a collective/FFT and cannot simply be deleted.
3. Deferred: clear only holes in middle/slab (save zero writes to subsequently
   overwritten entries), or useful-volume variable exchange. All send padding
   is already initialized only at allocation/growth, not zeroed per call;
   there is no repeated padding memset to remove. Sparse clear index traffic
   may cost more (PFFT-4). Variable alltoall needs a different interface and
   actual padding evidence. No gain assigned without the missing raw counts.

Ledger checked: A11 lane staging, A28 many-band 3D, A34 scoped FFT threads,
PFFT-4 indexed transpose, M2-ENTRY-COST-AB-DEV08/M2_FFT neutral, and
MPI-REMAINING-ASYNC-ALLTOALL rejected. This change is a single removed copy,
not new FFT batching, MKL switching, threading or async communication. Existing
M2d batches already reduce transpose call counts. Async remains deferred:
although the new group's share is larger than historical screens, wait/transfer
and overlap opportunity are not established. No physical constants, FFT size,
material, rank-count or band-count special case is introduced.
