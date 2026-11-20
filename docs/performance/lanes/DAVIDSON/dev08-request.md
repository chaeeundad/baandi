# Phase 2 request: upper-border GEMM tiles

Candidate `BAANDI_DAVIDSON_UPPER_BORDER=on`, default **off**. No default change
or adoption claimed. Read diagnosis.md for the relayed profile, candidate
comparison and missing raw shape/hardware inputs. The original Phase 0
request and its integrity evidence remain historical. No Mac timing evidence.

## Preconditions and counter gate

Return the original 2733278 results.json/tables.json (all five points) before
claiming a numeric roofline gap. Also return CPU/binding, effective kernel
rate or bandwidth calibration if available, and clarify why the quoted
exclusive split totals 104.1%. This does not prevent screening the code's
concrete work-elimination hypothesis.

On diagnostic A/B CLI runs, `davidson_profile.operations` now also labels:

* `upper_border_product`: n=local rows, m=full order, k=new width,
  detail=retained prefix. One parent record in BOTH arms for each affected
  product. Logical rectangle FLOP=2*n*m*k. This splits S/H products that the
  old combined Gram record could not distinguish by width.
* `upper_border_tile`: n=local rows, m=tile's last column (exclusive),
  k=tile width, detail=parent prefix. ON only. Actual FLOP=2*n*m*k.
  Nested in upper_border_product, itself inside Gram/overlap. Never sum
  their times. Existing logical solver FLOP counters retain their old
  interpretation; GEMM dispatch counters report executed work.

Mandatory counter, per rank: executed border FLOP must fall by exactly the
shape prediction `sum 2*n*(m*width - tile_end*width)` on equal call shapes.
At least one rank/call must have a positive saving. Zero logical work is
**unusable**, and zero saved work is **no mechanism exercised**, not a win.
No 5% minimum saving is imposed. Parent call counts/shapes, H columns,
Davidson iterations, full EVD counts/orders, rotations, restart/fallback
counts and reduction calls/payloads must stay equal for a rounding-class
claim. Collective payload does NOT shrink in this experiment. Actual GEMM
calls can increase. Reduction-site apply belongs to other lanes.

```sh
python3 -I docs/performance/lanes/DAVIDSON/border-counter.py \
  "$OUT/si64-np8-profile/ab-1-a/results.json" \
  "$OUT/si64-np8-profile/ab-1-b/results.json" > "$OUT/si64-border-work.json"
python3 -I scripts/lane_davidson_parse.py \
  "$OUT/si64-np8-profile/ab-1-a/results.json" \
  "$OUT/si64-np8-profile/ab-1-b/results.json" \
  --json-out "$OUT/si64-lane-tables.json" > "$OUT/si64-lane-tables.md"
```

Counter paths: `ranks[*].{full_rectangle_flop,predicted_tiled_flop,executed_flop}`
in border-work.json; raw records under `davidson_profile.per_rank[*].records`.
Cross-check `convergence.work.gemm_dispatch.rows[*].{site,route,calls,gflop}`
and `parallel.reductions`; do not use a zero serial solver_total as evidence.

## Build once, record provenance

Use the approved dev08 environment and external UPF, no server changes.
The orchestrator supplies the reviewed commit and absolute output paths.

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export CARGO_TARGET_DIR="$PWD/target/lane-davidson-phase2"
cargo build --release -p dft-cli --features mpi,mkl --bin dft-cli --example davidson_bench
export BIN="$CARGO_TARGET_DIR/release/dft-cli"
export BENCH="$CARGO_TARGET_DIR/release/examples/davidson_bench"
export UPF=/absolute/approved/Si.pz-vbc.UPF
export OUT=/absolute/lane-davidson-phase2-results
mkdir -p "$OUT"
git rev-parse HEAD > "$OUT/source-commit.txt"
git diff --binary > "$OUT/source.diff"
lscpu > "$OUT/cpu.txt"
rustc --version --verbose > "$OUT/rustc.txt"
mpirun --version > "$OUT/mpi.txt"
rg '__INTEL_MKL' "$MKLROOT/include/mkl_version.h" > "$OUT/mkl.txt"
```

Each invocation below is a **separate <=600 s job**. Never batch the five
points into one queue job. `timeout` enforces the cap and returns partial logs
on timeout. The shell runner uses physical CPU 0 (one core), or CPUs 0–7
socket 0 (MPI8), `--bind-to core --map-by core`, one MKL thread/rank. Confirm
physical placement with the same master procedure as Phase 0. Do not inherit
unrelated experimental BAANDI switches. Record resolved run config and hashes.

## Stage 0 first: one fixed-H warm solve

Initial job only (A/A floor, then AB/BA/AB; same binary, three pairs each):

```sh
timeout 590s bash docs/performance/lanes/DAVIDSON/phase2-job.sh si64 8 solve
python3 -I docs/performance/lanes/DAVIDSON/paired-summary.py \
  "$OUT/si64-np8-solve" solve > "$OUT/si64-np8-solve-summary.json"
```

The existing davidson_bench uses its start-density warm setup (documented in
fast-loop-protocol), an untimed warm-up, and one timed repeat per process;
repeats across the crossed pairs give the samples. It is real MPI here, NOT
`--sim-ranks`. Timed metric: `DAVBENCH_JSON.kinds[0].solve_ms_median`, external
process time is not the fixed-H metric. Counters and eigenvalues/residuals
are in the same record. Inspect all eigenvalues and occupied band-energy
sum (1e-11 Ha rounding class), independent residual/convergence checks and
work before judging speed. All detailed timers are OFF in timing jobs.

Run one separate diagnostic job to confirm the shape-derived work counter:

```sh
timeout 590s bash docs/performance/lanes/DAVIDSON/phase2-job.sh si64 8 profile
```

This profile job is a single full-convergence OFF/ON pair, detailed timers ON
in both arms. It is NOT acceptance timing. No statistical claim from one pair.

## Stage 1 and conditional expansion

Only after stage 0 has an active counter and no correctness/regression issue:

```sh
timeout 590s bash docs/performance/lanes/DAVIDSON/phase2-job.sh si64 8 screen
```

A screen sets max_iterations=10; exit 3 for expected non-convergence is
accepted only if results.json is complete and all ten history rows are
finite and match the A arm within 1e-11 Ha step by step. Other failures stop.
The diagnostic pair above also checks the converged tail/verification.

Conditional expansion: each command separately, never a combined job:

```sh
timeout 590s bash docs/performance/lanes/DAVIDSON/phase2-job.sh si32 8 screen
timeout 590s bash docs/performance/lanes/DAVIDSON/phase2-job.sh si16 8 screen
timeout 590s bash docs/performance/lanes/DAVIDSON/phase2-job.sh si32 1 screen
timeout 590s bash docs/performance/lanes/DAVIDSON/phase2-job.sh si64 1 screen
```

For each point, use a separate `profile` job for detailed attribution, then
`full` for the full-convergence accuracy/tail gate on survivors. These are
exactly the same runner arguments except the final mode. The one-core Si64
full job may reach the 590 s limit: preserve partial pairs, split follow-up
jobs in fresh output roots, and never silently pool different windows.
Summarize any timing job with `paired-summary.py JOB_DIR screen` (or `full`).
Do not infer serial gains from MPI or vice versa. No threads performance
comparison and no QE comparison are requested.

## Statistical threshold and withdrawal

For crossed pairs, d_i=(B_i-A_i)/A_i; report mean, sample standard error
s(d)/sqrt(N), all pairs and within-window range. Resolved gain requires mean<0
and abs(mean)/se >=3; report an exactly zero SE explicitly rather than divide
by zero. A/A floor is max absolute paired relative difference. On Si64 also
require effect >3 times that measured floor and separated mean +/-3SE
intervals for AA versus AB. Fewer than three complete pairs: no verdict.
An unresolved gain is inconclusive, not a correctness failure or evidence
for adoption. Adopt any resolved improvement, including 1–3%, only through
the orchestrator after generality/accuracy checks.

Withdraw/diagnose if the counter does not move, packing overhead makes wall
regress beyond noise, repeated within-arm eigenvalues differ, or any science
failure occurs. If trajectory/work changes, stop the rounding-class claim:
run full SCFs against the unchanged strict reference (strict 1e-10 Ha,
standard 5e-9 Ha), verification, and the generality set before proceeding.
Generality includes cubic/skew/elongated/slab, small/large bands, synthetic
multiple-projector operators, clustered spectra, simulated rank and pool
independence locally, then real ranks 1/2/4/8 and MKL on dev08. No switch
becomes default from silicon screens alone.

Return inputs/hashes, raw stdout/stderr/exit codes, all results.json and
science/density comparisons, A/A and A/B summaries, counters, load/placement,
and profile tables. Stop here; master decides the next phase.
