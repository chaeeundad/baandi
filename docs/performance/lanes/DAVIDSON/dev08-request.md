# Active request: PACKED_REDUCE (GEMM shapes unchanged)

UPPER_BORDER and ROTATE_COLUMNS are **REJECTED**; keep both OFF. Their
requests below are archived, including headings formerly labelled active.
Do not queue them again. Active commands are exclusively in this section.
No Mac timing is performance evidence.

## Mechanism and pre-registered gates

`BAANDI_DAVIDSON_PACKED_REDUCE=on` compacts consumed entries after the same
full GEMMs in distributed lean border reductions. No serial implementation
changes. No collective is merged or delayed past the pivot test. All buffers
already exist. The lower triangle and padded residual-norm tail are omitted.

On separate diagnostic A/B runs, parse:

```sh
python3 -I docs/performance/lanes/DAVIDSON/packed-counter.py \
  "$OUT/packed-si64-np8-profile/ab-1-a/results.json" \
  "$OUT/packed-si64-np8-profile/ab-1-b/results.json" > "$OUT/packed-si64-counter.json"
python3 -I scripts/lane_davidson_parse.py \
  "$OUT/packed-si64-np8-profile/ab-1-a/results.json" \
  "$OUT/packed-si64-np8-profile/ab-1-b/results.json" \
  --json-out "$OUT/packed-si64-tables.json" > "$OUT/packed-si64-tables.md"
```

Mandatory mechanism: for **every rank**, counter output `[0].ranks[*].sent_bytes`
(OFF) equals `rectangle_bytes`; `[1].ranks[*].sent_bytes` (ON) equals
`predicted_packed_bytes`, with positive `saved_bytes`. Calls and predicted
shape totals must match between arms. Prediction per border is
`8*(sum(j+1,j=m-ks..m-1)+sum(j+1,j=m-kh..m-1)+norm_count)`; OFF is
`8*m*(ks+kh+1)`. Raw `border_payload` records carry n=m,m=ks,k=kh,detail=norms;
the next nested `tree_sum` has n=actual doubles, detail=reduction site.
Check `parallel.reductions` (lean_border/reduced_gram) against those totals.
All other reduction sites must be unchanged; **call count is unchanged**.
Trace dropped=0, run.completed=true. On one-core CLI the path is inactive;
absence of border_payload is expected there, and packed-counter intentionally
rejects inactive evidence. Do not claim a serial mechanism or gain.

Keep exactly equal: GEMM route/call/shape/FLOP counters (including rotate),
full EVD counts/order histograms per map, H applications, restart/fallback
counts, tolerances, SCF history and verification. Target bit identity of
energy/density/scientific history at fixed ranks. Any observed rounding must
be explained by real MPI backend behavior; gate abs(dE)<=1e-11 Ha and
thread/rank independence, standard5e-9/strict1e-10 Ha. A trajectory change
withdraws this communication-only candidate pending investigation.

## Exact execution, each command a separate bounded job

Build once on dev08, source and binaries hashed. Set approved absolute UPF
and OUT paths in the orchestrator job environment. CPU0 for one core; CPUs0–7
socket0 for MPI8. Sequential MKL, one thread per rank. Set BIN/BENCH to the
matching build for each model. Use the same binary for both arms.

```sh
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/packed-mpi cargo build --release -p dft-cli --features mpi,mkl --bin dft-cli --example davidson_bench
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/packed-serial cargo build --release -p dft-cli --features mkl --bin dft-cli --example davidson_bench
export BIN="$PWD/target/packed-mpi/release/dft-cli"
export BENCH="$PWD/target/packed-mpi/release/examples/davidson_bench"
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si64 8 profile
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si64 8 solve
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si64 8 screen
# Expand only after the si64 mechanism passes and screens are promising:
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si32 8 profile
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si32 8 screen
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si16 8 profile
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si16 8 screen
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si64 8 full
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si32 8 full
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si16 8 full
# Serial negative controls; no claimed serial opportunity:
export BIN="$PWD/target/packed-serial/release/dft-cli"
export BENCH="$PWD/target/packed-serial/release/examples/davidson_bench"
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si32 1 screen
timeout 590s bash docs/performance/lanes/DAVIDSON/packed-job.sh si64 1 screen
```

Runner: three A/A pairs first, then three crossed A/B pairs, alternating order;
profile is one untimed-evidence pair with instrumentation enabled. Timed jobs
turn all lane/stage diagnostics off. Same seed601/displacement0.001Bohr,
standard stopping, verification included in full runs; screen max10 steps,
expected non-convergence exit3 only. Preserve hashes, exit codes, stdout/stderr,
results.json and external seconds. Returned profile times are inclusive.

Primary counter must move as above. Primary performance gate: full external
wall paired mean ON-OFF negative with abs(mean)/SE>=3 and larger than A/A
noise floor; screen and solve metrics are triage, never adoption. Also return
border_payload and tree_sum per-map/per-rank times (do not sum parent+child),
rank min/max, EVD work table, physical identity and absolute energies.
Use `paired-summary.py JOBDIR screen|solve|full` for paired statistics; do not mix parallel models. Reject resolved slowdown or numerical
failure; unresolved evidence needs more independent crossed pairs through the
orchestrator, not selective deletion. No fixed percentage minimum gain.
Actual MPI must check binomial/doubling/halving modes if deployed: payload
threshold crossings can change message routes even though elementwise tree
association is designed to agree. No thread-parallel optimization is proposed.

## Archived UPPER_BORDER request — closed, do not execute

## Preconditions and counter gate

The follow-up relay supplies dominant shapes, counts, payload bins and measured
machine rates; `diagnosis.md` now contains a corrected representative-shape
model. The earlier 104.1% split is withdrawn after the parent's duplicate-key
explanation. Exact trajectory-weighted floors still require joint shape/time
histograms; these are not a prerequisite for either controlled A/B. The
UPPER_BORDER commands below remain only as historical reproductions. The independent
ROTATE_COLUMNS experiment is appended at the end; do not stack the switches
for their first screens.

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

## Archived ROTATE_COLUMNS alternative: serial ROTATE_COLUMNS (independent request)

New default-off `BAANDI_DAVIDSON_ROTATE_COLUMNS` applies only to serial real
Ritz rotation pairs, standard and strict. It reuses the existing fixed
32-column GEMM adapter instead of 256-row scratch panels; it keeps X and H X,
residual checks and all work/stop decisions. The MPI distributed path is not
changed. The UPPER_BORDER screen is complete and rejected; do not repeat it.

The proposed gain is fewer internal MKL calls and scratch copies, **not fewer
FLOPs** or a new EVD algorithm. The serial representative-shape pair model is
26.0 ms, not 13.0 ms. The corrected EVD model is in diagnosis.md. Expectations
are conditional: one-core Si64 ~0.54–1.34 s (1.2–3.0% of prior external wall)
if 10–25% of non-residual rotation time is saved; MPI distributed effect zero.
Do not use those percentages as pass/fail minima or cross-model speed ratios.

### Counter contract

`davidson_profile.operations` adds `rotation_pair`; records have n=PW rows,
m=basis order, k=rotated columns, detail=0 (row panels) or 1 (column groups).
It encloses only the two products, inside serial rotate_residual, before
residual construction. Fused residual routes have no such record and retain
priority. The new runner explicitly sets RITZ_FUSED=off in both arms, which
is the current default. UPPER_BORDER=off in both arms.

In the separate diagnostic pair, require:

1. Equal nonzero pair counts, n/m/k sequence, 4nmk FLOP, H columns, EVD order
   distribution, SCF/Davidson iterations, restarts, verification, density and
   eigen residual/energy within the rounding gate. The detail flag changes
   from all 0 to all 1 on eligible calls. An MPI control with zero eligible
   records is labelled **inactive control**, not a measured improvement.
2. Existing actual dispatch counters in
   `convergence.work.gemm_dispatch.rows[*]`: `site=rotate` FLOP must decrease
   by the traced pair FLOP and `site=other` increase by exactly that amount;
   total FLOP is unchanged. Sum across routes; confirm `route=mkl` is active
   on the requested MKL points. A zero or fallback-only route is unusable for
   the MKL mechanism. The high-level wrapper calls still count two products;
   they are NOT the internal MKL call count.
3. The shape-derived internal-call model must fall:
   sum(2*ceil(n/256)) -> sum(2*ceil(k/32)); at 23871 rows/132 outputs,
   **188 -> 10** per pair. The extra scratch write/read model falls from
   sum(32*n*k) bytes to zero (100.831104 MB at that shape). These are exact
   structural counts on the MKL routes, not measured DRAM bytes. Report them
   as models and cross-check the actual dispatch movement from item 2.

```sh
python3 -I docs/performance/lanes/DAVIDSON/rotation-counter.py \
  "$OUT/rotation-si64-np1-profile/ab-1-a/results.json" \
  "$OUT/rotation-si64-np1-profile/ab-1-b/results.json" \
  > "$OUT/rotation-si64-np1-counters.json"
```

Output paths: `[*].ranks[*].total.{pairs,column_pairs,flop,chosen_calls_model,
extra_copy_bytes_model}`, with per-map values in `.maps`. All eligible pairs
must remain represented; dropped=0. Do not sum rotation_pair time with its
rotate_residual parent. New traces remain readable by lane_davidson_parse.py.

### Exact separate jobs

Use the same approved release build/provenance commands above, rebuilt at the
follow-up source, without `dgemm-timers` in acceptance binaries. `BIN`, `BENCH`,
`UPF` and `OUT` have the same meanings. `rotation-job.sh` is a separate runner;
it isolates the active alternative from the closed UPPER_BORDER experiment.

First fixed-H stage-0 job (one core, not MPI):

```sh
timeout 590s bash docs/performance/lanes/DAVIDSON/rotation-job.sh si64 1 solve
python3 -I docs/performance/lanes/DAVIDSON/paired-summary.py \
  "$OUT/rotation-si64-np1-solve" solve > "$OUT/rotation-si64-np1-summary.json"
```

After an active counter and no correctness issue, separate diagnostic and
fixed-ten-step jobs (never concatenate as one queue job):

```sh
timeout 590s bash docs/performance/lanes/DAVIDSON/rotation-job.sh si64 1 profile
timeout 590s bash docs/performance/lanes/DAVIDSON/rotation-job.sh si64 1 screen
```

Conditional expansion, one command per job:

```sh
timeout 590s bash docs/performance/lanes/DAVIDSON/rotation-job.sh si32 1 screen
timeout 590s bash docs/performance/lanes/DAVIDSON/rotation-job.sh si64 8 screen
timeout 590s bash docs/performance/lanes/DAVIDSON/rotation-job.sh si32 8 screen
timeout 590s bash docs/performance/lanes/DAVIDSON/rotation-job.sh si16 8 screen
```

MPI rows are control tests, not expected speedups. If a rank-local serial
solve appears outside the distributed solve, report its context separately
instead of silently attributing it to distributed algebra. For surviving
serial points run `full` to include verification/tail; request converged
skew/slab/synthetic cases and real-MPI rank independence at the master gate.

Same statistics as the UPPER_BORDER protocol: A/A floor first then AB/BA/AB,
at least 3 complete pairs in one window, mean<0 and abs(mean)/SE>=3. For
Si64 also effect>3*A/A floor and separated +/-3SE mean intervals. No minimum
percent gain. Correctness gate |dE|<=1e-11 Ha step by step if work/trajectory
unchanged; otherwise require full SCF with strict 1e-10 / standard 5e-9 Ha
and verification before further performance claims. Failure or regression
beyond noise withdraws this alternative as well; UPPER_BORDER remains rejected.

Risk being tested: full-row, 32-column calls may reread V/HV more often,
may have poorer skinny-GEMM utilization, and MKL grouping may change rounding.
Mac validates the general adapter path and fixed-group worker independence,
but cannot validate MKL speed or rounding. No new threading policy, relaxed
tolerance, residual omission or cached-H assumption is part of this change.

### Rejection-response additions to the acceptance record

UPPER_BORDER A/B was +0.005679 s, abs(mean)/SE4.18,3/3 positive; the master
rejected it. Do not override that decision with its correct work counters or
its small percentage. The next candidate must improve measured solve/wall,
not merely move a model. Reject ROTATE_COLUMNS if calls/copies move as
predicted but wall still worsens. No further upper-tile-size sweeps.

Use unambiguous dimensions in returned tables: PW rows N, basis inner M,
rotated roots K. Default internal products have output256xK with reductionM;
alternative products have outputNx32 with reductionM. Both compute the full
NxK result with unchanged4NMK pair FLOP. This is an independent output split
with fewer calls, not triangular row-prefix pruning. Rate at width32 remains
a risk; the existing width132 throughput is not its measured rate.

Return `.eigen_map_work` from the updated lane parser for each diagnostic
file. It contains rank/map/stage/kind/iteration, solve_ids, reduced_steps,
full_evd_calls, bordered_attempts, full_evd_order_counts and full_evd_sum_m3.
Require identical per-map EVD workload in A/B; never sum the nested eigen_step
and eigen times. Existing relay totals imply101 full EVDs/rank over all MPI
contexts, not101 per fine map. The serial rotation_pair versus residual child
separation must show whether the benefit came from the intended pair/layout
span; residual formula/norm checks stay unchanged.

Mac evidence is reused from the unchanged candidate: saved-base/OFF/ON
Si8/Si32/Si64 and skew/Si8 strict OFF/ON science+density bits identical,
Si32 workers1/3 identical; 13 solver tests and the scalar/tail/nonfinite/pool
unit test passed. The updated parser has a synthetic multi-map EVD/B5 test
and is checked on those existing Mac files without using their times.
