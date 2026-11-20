# INIT Phase 2: local atomic block screen

Status: implementation ready for review, **not adopted**. Base
`c0f2143769f0876931538e0e9d24c5e0c8202de0`; orchestrator records the reviewed
commit and binary hashes. Default-off `BAANDI_ATOMIC_LOCAL_BLOCK=on` reuses
atomic guarded CholQR2 on a single process; no new threaded or MPI algorithm.
Read [diagnosis.md](diagnosis.md) for the two candidates, floor assumptions,
and the A52 distinction. Mac evidence is correctness only.

## First job: stage 1, si32 one core, under 600 seconds

This changes the initial vectors' rounding; use full convergence with
verification, not a Davidson-only kernel screen. Run **only this point**
first. Build outside the job with the existing dev08 toolchain:

```sh
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/lane-init-one \
  cargo build --release -p dft-cli --features mkl
export INIT_BIN="$PWD/target/lane-init-one/release/dft-cli"
# UPF is the master's existing absolute Si.pz-vbc.UPF path.
export INIT_OUT="$PWD/target/init-local-si32"
mkdir -p "$INIT_OUT"
git rev-parse HEAD
sha256sum "$INIT_BIN" "$UPF"
lscpu -e=CPU,CORE,SOCKET,ONLINE
```

One idle physical CPU 0, socket 0; same binding both arms. In a clean shell,
clear inherited BAANDI experiments and record the environment. No change to
cutoff, bands, start, stopping rule, verification or backend between arms.
Fresh output directories, one A/A floor immediately followed by crossed A/B
inside the same window. These commands run on **dev08 only**:

```sh
export BAANDI_INIT_PROFILE=off BAANDI_PHASE_PROFILE=off BAANDI_STAGE_TIMERS=off
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# Three pairs in AB, BA, AB order; A/A establishes the local noise floor.
taskset -c 0 python3 -I scripts/fast_screen.py --binary "$INIT_BIN" --upf "$UPF" \
  --case benchmarks/cases/si32-gamma.toml --threads 1 --seed 601 --pairs 3 \
  --env-a BAANDI_ATOMIC_LOCAL_BLOCK=off --env-b BAANDI_ATOMIC_LOCAL_BLOCK=off \
  --out "$INIT_OUT/aa" --json "$INIT_OUT/aa.json"
taskset -c 0 python3 -I scripts/fast_screen.py --binary "$INIT_BIN" --upf "$UPF" \
  --case benchmarks/cases/si32-gamma.toml --threads 1 --seed 601 --pairs 3 \
  --env-a BAANDI_ATOMIC_LOCAL_BLOCK=off --env-b BAANDI_ATOMIC_LOCAL_BLOCK=on \
  --out "$INIT_OUT/ab" --json "$INIT_OUT/ab.json"
python3 -I docs/performance/lanes/INIT/paired_stats.py \
  "$INIT_OUT/aa.json" "$INIT_OUT/ab.json" > "$INIT_OUT/paired.json"
```

Timing decision is independent of the generic harness verdict. Let d_i be
B wall minus A wall, SE=sample_sd(d)/sqrt(3). Require negative mean,
**|mean|/SE >= 3**, and improvement >3 times max absolute A/A paired difference.
If SE=0, nonzero mean is the limiting infinite statistic (reported as a flag,
not invalid JSON infinity). Report all pairs, medians, within-window range,
mean/SE and A/A floor. A noisy screen is inconclusive, not a rejection or
adoption. No minimum percentage gain. Enforce convergence, verification,
within-arm deterministic density/counters; report any SCF, Davidson,
H-application or verification-work change. Energy gate is 1e-11 Ha cross-arm
for rounding-only results; if exceeded, require complete strict-reference
runs, strict 1e-10 Ha / standard 5e-9 Ha before proceeding.

## Separate diagnostic job: same si32 point

Do not use instrumented wall as acceptance timing. Repeat the exact A/B
screen above with `INIT_OUT` changed to a fresh `target/init-local-si32-profile`,
`BAANDI_INIT_PROFILE=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on`, and
`--pairs 1`. Retain both full results.json files and parse:

```sh
python3 -I scripts/lane_init_parse.py "$INIT_OUT/ab" \
  --json-out "$INIT_OUT/lane.json" --markdown-out "$INIT_OUT/lane.md" \
  --label 'dev08 Phase 2 diagnostic only; local atomic block'
```

Pre-registered mandatory counter movement (JSON paths):

- `run.counters["init_lane.atomic.mgs.calls"]`: positive A -> zero B.
- `run.counters["init_lane.atomic.mgs_projection_macs_upper_bound"]`: positive A -> zero B.
- `run.counters["init_lane.atomic.local_accepted"]`: zero A -> positive B,
  equal to `local_attempts`; `local_fallbacks` must be zero for this screen.
- B has two `block_gram.calls`, two `block_factor.calls`, two
  `block_update.calls` per accepted span. Prefix for all is
  `run.counters["init_lane.atomic.`; each operation has a matching `.ns`.
- `rows`, `columns`, `bands`, `grid_points`, `grid_x`, `grid_y`, `grid_z`,
  `input_bytes`, `block_gram_flops`, `block_update_flops`
  supply the floor inputs. Counts sum across starts, not matrix dimensions
  of a synthetic combined start. Use `effective_start` for the stage and fallback/completion; top-level
  `basis` and `atomic_start` may describe only the fine stage. Grid counters
  above identify the start grid (divide by calls only for identical shapes).

All-zero required counters mean UNUSABLE, never success. The old
`atomic_start.mgs_projection_complex_macs` is a historical proxy even for block paths;
**do not use it** as the work movement gate. `atomic.span.ns` is inclusive;
block fill/Gram/factor/update/gather are component timers, but not an exhaustive
partition. `block_gram.ns` includes reduction on MPI. The new timers do not
add collectives. Existing `phase_profile.buckets` init atomic envelopes retain
their historical meaning; MPI HQ is partly inside Gram. Report shapes and
arithmetic estimates, not measured DRAM bandwidth inferred from logical bytes.

## Conditional widening, separate jobs only after a candidate

1. Repeat the first and diagnostic jobs with `si64-gamma.toml` and fresh
   si64 directories. Require same-window crossed pairs and the measured 3x
   A/A floor; never subtract an earlier table's wall. Confirm the hypothesized
   coarse m=8409, n=256, bands=132 and FFT=54^3 from returned start metadata.
2. Run the existing Phase 0 profile-request MPI8 commands on si16/si32/si64,
   separately, both switch settings, CPUs 0–7 socket 0. These are negative
   controls: `local_attempts=0`, exact same scientific outputs and work counts
   within each fixed-rank model. New block timers describe remaining MPI
   work, but no MPI gain is predicted or claimed for this candidate.
3. Full accuracy/generality before adoption: Si skew/elongated/slab, shifted k,
   small/large bands, other pseudopotential channel structures, strict
   reference and standard certificate, both pure and MKL backends. Run each
   case as its own <=600s job. Guards/fallbacks must remain observable.

Return source/binary/UPF/input hashes, environment/binding, A/A and A/B JSON,
paired.json, diagnostic results.json and parsed lane tables, convergence and
verification, energy deltas, density hashes and work. No QE run is requested;
any future comparison stays within the same parallel model. Stop for the
orchestrator's decision; keep the new switch off pending evidence.
