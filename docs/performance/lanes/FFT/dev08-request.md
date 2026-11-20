# LANE-FFT Phase 2: kept copy/multiply fusion A/B

Status: ready for orchestrator review; default OFF, not adopted. See
[diagnosis](diagnosis.md) for relayed bad02ee profile, limitations, floors and
alternatives. Candidate source base 10c067e plus the uncommitted patch; record
actual merged commit and binary hashes before running. No Mac timing evidence.
No generic MKL FFT switch or async collective is proposed.

## Updated five-point motivation (bad02ee relay)

All five count-only controls agree with timed counters (zero differences).
Coverage is lane **rank maximum** divided by run wall, resolving the earlier
si64 denominator question: 2.2649/8.4265 = 26.88%. MPI8 lane coverage for
si64/si32/si16 is 26.88/26.05/23.58%; lane rank spread is 3.89/0.58/0.86%.
The visible wait/asymmetry is specific to si64, not a general MPI8 premise.
The separate one-core profiles cover 24.49% (si64) and 31.44% (si32), with
no alltoall or kept_copy. See diagnosis.md for complete per-model tables.

Permutation/zeroing consumes about one third of the lane at every point
(31.1–36.0%). MPI8 middle_pack delivers 6.49–7.63 GB/s per rank, about
52 GB/s aggregate at the low end, consistent with shared socket bandwidth
saturation. **Removing passes/bytes is the primary lever**, with padded
alltoall volume the second lever at larger sizes: forward collective rates
fall from 4.20 GB/s at si16/si32 to 2.38 GB/s at si64. Existing axes track
volume over a 57x points/call range at 1.29–2.15 ns/point; these data do not
justify reopening generic batching, MKL FFT or async alltoall experiments.

The selected kept-copy fusion is a bounded first byte-removal experiment,
not a solution to the entire permutation group. Its affected lane share is
5.0% at si64 MPI8 and 3.1% at si32/si16; serial is already fused. Smaller-cell
A/B therefore matters even if si64 alone shows a resolved gain. Report any
collective-time movement separately: this patch changes neither collective
payload nor partitioning, and cannot claim an alltoall-volume reduction.

The pre-registered counter movement, diagnostic sum, identity gates, A/A
floor, |mean|/SE >= 3 threshold and staged commands below remain unchanged.

## Pre-registration

Same binary and input, A `BAANDI_FFT_KEPT_FUSE=off`, B `on`, MPI8 with one thread
per rank. The only intended effect is eliminating the separate kept-plane
copy immediately before multiplication. Existing serial path is unchanged:
no one-core speed claim, no one-core timing request for this candidate.

Required deterministic counter movement in ordinary `fine` bucket, each rank:

- `phase_profile.fft_lane.per_rank[r].buckets[b].stages[s]` selected by
  `rank`, `bucket == "fine"`, `stage == "kept_copy"`: A.calls > 0 and
  A.logical_bytes > 0; **B.calls = B.points = B.logical_bytes = 0 (100% drop)**.
- Same path for `multiply`: calls, points and logical_bytes unchanged (40
  bytes/point in both arms). Sum over all stages of logical_bytes drops by
  **exactly A.kept_copy.logical_bytes = 32 * A.kept_copy.points**. FFT work,
  collective calls/payload and maps stay identical. Sum of two old passes
  on kept points is 72 bytes/point versus 40 in B. Zero in both arms is
  UNUSABLE, not a success. Reject unexpected changes in other deterministic
  counters. Count-only runs use FFT_PROFILE+PHASE_PROFILE on, STAGE_TIMERS off.
- A diagnostic timed follow-up, only after a counter survivor, must reduce
  **same-rank fine (kept_copy.seconds + multiply.seconds)/maps**. Use the sum
  to avoid crediting mere timer relabeling; report every rank, rank maximum
  of the sum, calls and work-normalized time. No guarantee that time saved
  equals the removed timer. Timed instrumentation is never wall evidence.

End-to-end wall/solve statistics: paired d_i = B_i/A_i - 1; use the process
median warm solve for stage 0, external launch-to-exit wall for SCF stages.
Mean d < 0 and **|mean(d)| / (sample_sd(d)/sqrt(n)) >= 3** are required for
resolved gain. A/A floor f = max(abs(mean(d_AA)), SE(d_AA)); also require
abs(mean(d_AB)) >= 3 f and separated intervals mean_AB +/- 3 SE_AB versus
mean_AA +/- 3 SE_AA. Return raw pairs, mean, SE, intervals, range and ratio;
SE=0 is not infinite confidence (collect more pairs). Three pairs screen,
not strong adoption evidence. If unresolved, return NEUTRAL and request a
larger **same-window** crossed run sized below 600 s, without pooling absolute
walls from different windows. No minimum percentage replaces this criterion.

Identity gate is **bit-identical at fixed rank count**: total energy,
eigenvalues, occupations, scientific SCF history, density bytes/hash and work
counts; verification outcome unchanged on converging runs. The permitted
rounding-class 1e-11 Ha gate is not needed by this change. Different rank
counts retain ADR 0041's strict 1e-10 / standard 5e-9 Ha contracts, not a new
promise of cross-rank SCF bit identity. No coefficient reduction/FFT arithmetic
order is changed. Withdraw for failed identity, mutated kept state, allocations
in the fused loop, nonzero kept_copy in B, unexplained changed work, or a
resolved wall regression on any requested cell. No adoption from Si alone.

## Build once (outside measurement jobs)

Use the same dev08 MPI/MKL environment as the returned profile. No builds or
other lane jobs overlap the timing window. `DFT_SI_UPF` points to the matched
external file; never copy its contents into the repository.

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export FFT_TARGET="$PWD/target/lane-fft-phase2-mpi"
CARGO_TARGET_DIR="$FFT_TARGET" cargo build --release -p dft-cli --features mpi,mkl \
  --bin dft-cli --example davidson_bench
sha256sum "$FFT_TARGET/release/dft-cli" "$FFT_TARGET/release/examples/davidson_bench"
git rev-parse HEAD
```

[dev08-pairs.sh](dev08-pairs.sh) is the exact runner. Its binding command is
the Phase-0 `mpirun -np 8 --map-by ppr:8:socket --bind-to core
--report-bindings`. Inspect reported bindings/topology before accepting results:
eight distinct physical cores, CPUs 0–7, socket 0, matching the returned profile.
Reject a run with different binding; never silently continue unbound. Runtime
implementation and baseline switches stay identical; save the results' full
switch manifest. Use a fresh work directory per job; timeout 580 leaves room
inside the 600 s queue limit.

## Stage 0 first: si32 np8, one fixed-H warm solve

Two separate jobs in order: deterministic counters, then three A/A and three
order-alternated A/B pairs. The benchmark does its standard untimed re-arm,
then warm solves (five repeats/process, process medians). No convergence-loop
speed inference from this stage. Start-density setup is common to both arms;
report `base_density`, warm tolerance and eigenvalue FNV. Phase-profile totals
include setup/re-arm and repeats, so they establish path/counter movement;
benchmark `kinds[kind="warm"]` supplies solve-local identity and timing.

```sh
export FFT_WORK="$PWD/target/fft-phase2-si32-solve-counts"
timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si32 solve counts
```

Counter survivor only:

```sh
export FFT_WORK="$PWD/target/fft-phase2-si32-solve-time"
timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si32 solve time
```

`DAVBENCH_JSON` log record: `phase_profile.fft_lane` has the above counter
paths; `kinds[kind="warm"].solve_ms_median`, `.eigen_fnv`, `.eigenvalues`,
`.h_columns`, `.davidson_iterations`, `.repeats_identical`, `.converged`,
`.counters`. A/A repeats and A/B must match eigenvalue FNV and work. Never
interpret missing/all-zero counters as a neutral measurement.

## Stage 1 survivor: si32 np8, ten-step SCF

Separate counts/time jobs, same mechanism; the runner changes only
`scf.max_iterations` to 10. Exit 3 is expected non-convergence when the cap is
reached; crashes, missing output, divergent history or unequal work invalidate
the screen. Compare all step energies/residuals, coarse and fine maps/counts.

```sh
export FFT_WORK="$PWD/target/fft-phase2-si32-steps-counts"
timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si32 steps counts
```

```sh
export FFT_WORK="$PWD/target/fft-phase2-si32-steps-time"
timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si32 steps time
```

## Stage 2 expansion: converging si64, si32, si16 np8

Only after stage 0/1 survives. One case per job, no multi-case shell grid. For
each case run the counts job first and its time job second. The converging
input keeps default verification on: tail/verification costs are included.
Explicit commands (each `timeout` command is a separate queue job):

```sh
FFT_WORK="$PWD/target/fft-phase2-si64-full-counts" timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si64 full counts
FFT_WORK="$PWD/target/fft-phase2-si64-full-time" timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si64 full time
FFT_WORK="$PWD/target/fft-phase2-si32-full-counts" timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si32 full counts
FFT_WORK="$PWD/target/fft-phase2-si32-full-time" timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si32 full time
FFT_WORK="$PWD/target/fft-phase2-si16-full-counts" timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si16 full counts
FFT_WORK="$PWD/target/fft-phase2-si16-full-time" timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si16 full time
```

Diagnostic timed profiles: a separate job for one converging case, three
alternating A/B pairs, FFT_PROFILE=on, PHASE_PROFILE=on, STAGE_TIMERS=on in both
arms. Repeat si32/si16 as separate jobs only as needed. Do not substitute
these walls for the uninstrumented jobs. Return fine/coarse/verification
buckets separately and every rank.

```sh
FFT_WORK="$PWD/target/fft-phase2-si64-full-profile" timeout 580 bash docs/performance/lanes/FFT/dev08-pairs.sh si64 full profile
```

Return every input/results.json, DAVBENCH_JSON/log, `.exit`, `.wall`, density
hash, binary/UPF hash, load and binding record. Parse CLI counter files with:

```sh
python3 -I scripts/lane_fft_parse.py "$FFT_WORK/counter-A/results.json" \
  "$FFT_WORK/counter-B/results.json" > "$FFT_WORK/lane.csv"
```

For remaining Phase-0 questions, return raw bad02ee rows rather than rounded
summaries: inverse useful bytes = transpose_pack.logical_bytes/2; forward
useful = middle_pack.logical_bytes/2; compare each with already-padded
alltoall.logical_bytes. Separate stage rank maxima from rank-0 totals. No
extra collective barriers or static reassignment are part of this candidate.
