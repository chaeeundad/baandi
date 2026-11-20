# LANE-GLUE: serial ordered density fusion A/B

Status: implementation ready for relay; **default OFF**, no adoption claim.
The 2026-11-20 follow-up explicitly authorizes this serial density/FFT change,
superseding the historical closing note and original grid/mixer-only scope.
Mac is correctness only. No remote job was executed by this agent.

## Hypothesis and boundaries

New orchestrator evidence: dev08 one-core STREAM triad 12.97 GB/s and
single-thread DGEMM 64–108 GFLOP/s. Source of prior SCF profile: e947a35.
For 15 maps, 128 occupied bands, 72^3 points, conservative accounted traffic
saves 11,466,178,560 bytes / 12.97e9 = **0.884054 s**, about 2.17% of 40.7 s.
This is a hypothesis, not measured removed-pass bandwidth. Cache residency,
FFT cost, write allocation and retained psi storage may reduce the gain.
DGEMM throughput is not an FFT throughput bound.

`BAANDI_SERIAL_DENSITY_FUSED=on` applies only to a standard map, one worker,
Hermitian paired bands, DENSITY_PAIR and DENSITY_REDUCE enabled, a non-MPI
build and no registered communicator. MPI builds, including rank one,
strict reference, complex/unpaired mode and thread budgets >1 retain the old
path. Odd or singly occupied bands retain a split scratch grid and ascending
addition. No communication, GEMM, solver, mixer, stopping or verification changes.
Both ordered pure-Rust and full-grid/MKL FFT adapters implement the same sink.
Each contribution retains its division and addition from +0; bands are added
separately in ascending order. REUSE_PSI writes remain unchanged. Per-band
storage is retained for fallback, so no resident-memory reduction is claimed.

## Correctness and structural counters (separate from timing)

Build the reviewed candidate once; A and B use that same executable. Preserve
source hash, diff, binary hash, UPF hash, compiler/MKL versions, binding, input
and all outputs. Use the same UPF and input generation as the profile window.
Run in an exclusive socket-0 window; no MPI launcher, no concurrent jobs.

```sh
export UPF=/data/sig-dft/pseudo/Si.pz-vbc.UPF  # set to the profile's actual hashed file
export ROOT=/data/sig-dft/lanes/GLUE/serial-fusion
mkdir -p "$ROOT"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/lane-glue-fusion \
  cargo build --release -p dft-cli --features mkl
export BIN="$PWD/target/lane-glue-fusion/release/dft-cli"
sha256sum "$BIN" "$UPF"
git rev-parse HEAD
git diff --stat
```

One case per job; run this block separately with CASE=si64, then si32 (si16
optional). This is the same single-process one-core model in both arms.

```sh
export CASE=si64
python3 -I scripts/fast_screen_dav.py --upf "$UPF" --case "$CASE" \
  --threads 1 --seed 601 --amp 0.001 --emit-input "$ROOT/$CASE.toml"
taskset -c 0 python3 -I docs/performance/lanes/GLUE/fusion_verify.py \
  --binary "$BIN" --work "$ROOT/$CASE-counts" \
  --json-out "$ROOT/$CASE-counts.json" "$ROOT/$CASE.toml"
```

The verifier sets GLUE_PROFILE=on, PHASE_PROFILE=on, STAGE_TIMERS=off and
compares scientific bits including density.bin, SCF history and verification.
Counter paths in each results.json are `run.counters["glue_lane.density.X"]`:

| X | pre-registered result |
|---|---|
| eligible_calls | A=B; includes eligible coarse/fine/verification calls |
| fused_calls | A=0; B=eligible_calls >0 |
| paired_band_points | A=B=sum(2 * executed occupied pairs * local grid points) |
| logical_bytes | A=32*live-band-points; B=A-16*paired-band-points |
| saved_bytes | A=0; B=16*paired-band-points=A.logical_bytes-B.logical_bytes |

These are conservative **logical model** counters, not hardware memory events.
The old 32 counts split clear 8, update 16, fold read 8; the new 16 debits
rho read/write per fused band. Old reduction-output traffic is deliberately
not credited, and retained FFT/psi traffic is excluded. Counters aggregate
actual calls/grids/occupations; do not force a whole-run counter to equal the
15-fine-map model when coarse/verification counts differ. For exactly that
model, paired_band_points=716,636,160 and saving=11,466,178,560 bytes.
Work/FFT call counts, map counts, convergence, verification, final density and
all energy/history bits must agree. Counter mismatch or changed work withdraws
the candidate pending diagnosis. No skipped verification is permitted.

## A/A floor first, then crossed A/B

For each case, first run three crossed A/A pairs (both OFF), inspect the
paired differences, then six A/B pairs. If the A/A mean is resolved at
|mean|/SE>=3, the window is invalid; repair exclusivity/binding and repeat the
whole measurement, without subtracting the drift. Zero empirical SE does
not establish resolution; collect more samples. Separate A/A and A/B jobs
keep Si64 jobs below approximately 10 minutes each. Keep all raw results.

```sh
export BAANDI_GLUE_PROFILE=off BAANDI_PHASE_PROFILE=off BAANDI_STAGE_TIMERS=off
# A/A: run before A/B; do not overlap jobs.
taskset -c 0 python3 -I scripts/fast_screen.py --binary "$BIN" --upf "$UPF" \
  --case "benchmarks/cases/$CASE-gamma.toml" --threads 1 --seed 601 --amp 0.001 \
  --pairs 3 --env-a BAANDI_SERIAL_DENSITY_FUSED=off \
  --env-b BAANDI_SERIAL_DENSITY_FUSED=off --energy-gate 1e-11 \
  --out "$ROOT/$CASE-aa" --json "$ROOT/$CASE-aa.json"
python3 -I docs/performance/lanes/GLUE/fusion_stats.py --aa "$ROOT/$CASE-aa.json"   > "$ROOT/$CASE-aa-floor.json"
# Inspect aa_window_suspect=false and nonzero SE before launching the next job.
# A/B: fixed six pairs, alternating AB / BA. Same input and delivered result.
taskset -c 0 python3 -I scripts/fast_screen.py --binary "$BIN" --upf "$UPF" \
  --case "benchmarks/cases/$CASE-gamma.toml" --threads 1 --seed 601 --amp 0.001 \
  --pairs 6 --env-a BAANDI_SERIAL_DENSITY_FUSED=off \
  --env-b BAANDI_SERIAL_DENSITY_FUSED=on --energy-gate 1e-11 \
  --out "$ROOT/$CASE-ab" --json "$ROOT/$CASE-ab.json"
python3 -I docs/performance/lanes/GLUE/fusion_stats.py \
  --aa "$ROOT/$CASE-aa.json" --ab "$ROOT/$CASE-ab.json" > "$ROOT/$CASE-stats.json"
```

`fast_screen.py`'s old percentage verdict is **not** the adoption rule; use
its raw external-wall observations and fusion_stats.py. Define paired saving
D_i=t_A-t_B, SE=sample_stddev(D)/sqrt(n). Require mean>0 and |mean|/SE>=3,
unchanged accuracy/work, no other-case regression beyond noise. Any resolved
gain counts; 0.8 s is the design estimate screen, not a rejection threshold.
If unresolved, report unresolved; pre-register a fresh larger fixed sample
before another batch, not repeated peeking until significance. Report all
cases, including negative results. No QE comparison is newly measured here.

Expected class is **bit-identical**. Any changed density/history invalidates
this implementation's claim even if |dE|<=1e-11 Ha; investigate first. A proposed
rounding/trajectory variant would need the full strict 1e-10 Ha / standard
5e-9 Ha accuracy and certificate gates, and independent rank-count checks.
Run normal full dev08 validation before adoption, including MPI build/tests;
this Mac cannot exercise MKL or real MPI. Preserve the existing disabled MPI
behavior, including np1 MPI worlds. No MPI timing job is requested for a disabled
candidate.

## MPI8 alternatives requested by the orchestrator

No demonstrated >0.2 s lever is identified. `m2.xc`=0.1559 s includes potential
and audit/FFT/collective work, not just XC; 37.4% spread corresponds to at most
0.0583 s under spread/max, not 0.2 s. Potential=0.1226 s versus QE v_of_rho=0.03 s:
the observed gap is about 0.09 s. Equal fine-grid plane counts do not explain
the timing skew; raw pure-XC spans/arrival data are needed to identify its cause.

Energy=0.2865 s has **no direct QE counterpart**; QE's ~0.05 s electrons rest
is not a matched energy timer. Removing 0.2 s requires ~70% of the whole energy
row; the explicit nonlocal sub-item is only 0.1161 s. The standard Coulomb stop
estimator does not need every reported audit scalar, but current finite-energy
checks, histories, strict gates and verification do consume fresh values.
No acceptance-preserving deferred-audit implementation with >0.2 s attributable
work is established; therefore it is not proposed here. Existing shared
Hartree/local/stop FFT and EWALD/IONIC paths are already adopted, and the prior
ENERGY_IDENTITY substitution lost 0.3–1.2% wall. The STREAM premise changes
serial density's decision, not that energy result.
