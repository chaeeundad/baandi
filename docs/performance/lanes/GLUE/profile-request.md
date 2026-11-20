# LANE-GLUE — Phase 0 profile request

Status: **PROFILE_REQUEST**. Base `61822bca3b470df1212ee47b47c7eb4e786e22e3`;
changes are uncommitted for the orchestrator to review. Record the reviewed
commit, source diff and binary hashes before running. This is instrumentation
only: no optimization, default change, candidate selection or adoption claim.
The 2026-11-20 charter revision supersedes the original Phase 2 instructions.

## Question and existing context

The user relayed an instrumented dev08 si64 MPI8 profile: density 0.746 s,
energy 0.291 s, potential 0.134 s, mixing 0.032 s, materialization 0.296 s,
between-map work 0.364 s (session setup 0.235 s, handoff 0.123 s), and roughly
0.84 s outside the process. These are supplied context, not measurements from
this run; the underlying raw run identity was not supplied. No QE ratio is
inferred. Determine whether density transforms, grid XC/Hartree, local energy
integration, or mixer arithmetic/communication dominate the remaining GLUE
work on each case/model. Startup and launch remain separate envelopes.

Read ledger rows M2-GRID-MIX-DEFAULT-ON, M2-OWNED-DEFAULT and the existing
MIX_INCREMENTAL/Coulomb-guard history before selecting any Phase 2 candidate.
No rejected optimization is retried here.

## Instrumentation and exact JSON paths

Enable `BAANDI_GLUE_PROFILE=on` **and** `BAANDI_PHASE_PROFILE=on` before
constructing plans/mixers. Both gates are required; GLUE defaults off.
`BAANDI_STAGE_TIMERS=off` counts calls/sizes without reading new clocks.
Tokens are cached at construction. Disabled sites only check an `Option`;
there are no new counter updates, allocations, collectives, clocks, logging
or diagnostic size traversals when disabled. Enabled hooks sit outside
point/band loops and use fixed run-local counters, no trace allocation.
RAII closes timers on error; call counts include entered failing calls.
No numeric expression, summation order, solver branch or tolerance changes.

All new paths are `run.counters["glue_lane.<name>.calls"]` and
`run.counters["glue_lane.<name>.ns"]`. Values are **local to this rank's
process and cumulative over the complete run**, including coarse sessions,
verification, and mixer initialization. Nanoseconds / 1e9 gives seconds;
seconds / entered calls gives the inclusive average call cost.

| name | span |
| --- | --- |
| grid.potential | distributed input potential, including all children |
| grid.audit | distributed output-density energy audit, including all children |
| grid.forward | density validation/packing, forward planes and transpose |
| grid.hartree | per-stick contraction of kernel times squared coefficients |
| grid.inverse | inverse transpose, plane transforms and real extraction |
| grid.xc | per-plane LDA evaluation and XC energy/integral accumulation |
| grid.assemble | local ionic + Hartree + XC effective potential assembly |
| grid.pack | partial-sum/optional field packing and receive sizing |
| grid.gather | partial-sum/optional field allgatherv only |
| grid.reduce | canonical plane/stick sums and flag reduction |
| grid.audit_local | five local charge/double-counting/linear/local integrals |
| dist.mix | distributed block Pulay including children (not final full-density gather) |
| dist.prepare | validation, residual/differences/history and linear trial |
| dist.partial | chunk partials and buffer sizing (Gram/RHS and charge calls pooled) |
| dist.gather | chunk allgatherv only |
| dist.reduce | canonical summation of gathered chunks |
| dist.solve | replicated small spectral solve and coefficients |
| dist.update | application of Pulay coefficients to local points |
| dist.charge | both charge corrections, including nested chunk reductions |
| serial.mix | replicated Pulay, including initialization and all children |
| serial.fit | existing incremental or full fit/update branch, including no-history entry |
| serial.norms | residual/step norms and diagnostic construction |
| serial.charge | final two uniform charge corrections |

Parents and children overlap. In particular `dist.partial/gather/reduce`
include calls inside `dist.charge`; do not add these to charge or total mix.
Grid gather/reduce counters also include the stop estimator when it uses the
same helpers; they are not limited to the potential/audit parents.
Forward and inverse include FFT communication. `grid.hartree` does not
include multiplying the kernel into the inverse-transform coefficients;
that and small unlabelled work remain in the potential parent remainder.
Disabled/zero-call stages are not zero-cost evidence.

Work paths use `run.counters["glue_lane.<suffix>"]`:

| suffix | definition / floor input |
| --- | --- |
| grid.local_points, grid.local_coefficients | sums of rank-owned plane points / complex coefficients over forward calls |
| grid.xc_points | sum of rank-owned XC points over XC calls |
| grid.send_values, grid.recv_values | partial/field gather f64 send/receive payloads; multiply by 8 bytes |
| dist.local_points, dist.history | sums of block points and active differences over successful preparation |
| dist.dot_elements | `(new Gram entries + RHS entries) * local points`, summed; two real flops/element |
| dist.update_elements | differences * local points, summed; four real flops and 32 logical bytes/element for the ordinary update |
| dist.send_values, dist.recv_values | chunk-reduction f64 payloads, including charge; multiply by 8 bytes |
| serial.points, serial.history, serial.update_elements | corresponding replicated mixer sizes; not an exact Gram dot count |

These are logical arithmetic/payload inputs, not measured DRAM traffic or
MPI wire bytes. Repeated scalar coefficient reads/cache traffic and allocation
metadata are excluded. No bandwidth assumption or speedup is inferred.

Existing complementary paths, retained by the parser:

- `convergence.history[].timings.{potential_s,density_s,total_energy_s,mixing_s}`:
  local **fine-iteration** envelopes; no reconstruction of substage traces.
  `fixed_point_certificate.timings` is separately labelled when present.
  Standard verification/coarse maps use the phase buckets below.
- `phase_profile.buckets[]`: coarse, fine, fine_redo, verify_pass, verify_fail,
  init, materialize and between; `.min`/`.max` are already rank-reduced
  envelopes. Never sum maxima from different ranks into a wall time.
- `phase_profile.per_rank[]`: original rank-local bucket closure data.
- `phase_profile.fft_lane.per_rank[].buckets[].stages[]`: existing FFT_PROFILE
  calls, points, logical bytes, estimated FFT flops and stage seconds. These
  pool FFT uses within each bucket; they are **not density-only**. Compare
  them with the density envelope and existing `m2.density`/work FFT counts,
  without assigning all inverse FFT work to density.
- `parallel.m2`, `work`, `basis`, `bands`: existing density/kept-plane/gather
  counts and final sizes. One-core grid stages are outside the owned files;
  they retain existing per-map envelopes and FFT diagnostics, not the new
  distributed-grid subspans. Serial Pulay is newly detailed.

New grid/mixer subspans are cumulative means, not individual map/call traces.
They cannot isolate a coarse versus verification substage by subtraction.
This bounded transport reuses the CLI counter schema within lane ownership;
no shared profiler or CLI writer changes. If Phase 1 needs individual
substage traces, return that requirement with the data before Phase 2.

## Exact dev08 jobs (master relay only)

Use the existing dev08 toolchain, pseudopotential and exclusive socket-0
allocation. Install nothing. One case/model per job, maximum 600 seconds;
return partial outputs and failure codes if exceeded. Build outside jobs:

```bash
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/lane-glue-one \
  cargo build --release -p dft-cli --features mkl
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/lane-glue-mpi \
  cargo build --release -p dft-cli --features mpi,mkl
git rev-parse HEAD
lscpu -e=CPU,CORE,SOCKET,ONLINE
sha256sum target/lane-glue-{one,mpi}/release/dft-cli "$UPF"
```

`UPF` is the master's existing absolute Si.pz-vbc.UPF path. Verify the CPU
IDs below are distinct physical cores on socket 0; replace them if required
by the allocation and record the actual binding. No multithreaded jobs.

| CASE | NP | MODEL | CPUSET |
| --- | ---: | --- | --- |
| si64 | 8 | mpi | 0-7 |
| si32 | 8 | mpi | 0-7 |
| si16 | 8 | mpi | 0-7 |
| si32 | 1 | one | 0 |
| si64 | 1 | one | 0 |

Run this block separately for each row (first row shown). Use a fresh JOB;
retain input, binaries/hashes, all raw results, stdout/stderr and exit codes
under `/data/sig-dft/lanes/GLUE/`. Record all inherited BAANDI overrides;
use reviewed defaults except the explicitly listed diagnostic/guard modes.

```bash
set -euo pipefail
CASE=si64 NP=8 MODEL=mpi CPUSET=0-7
JOB="/data/sig-dft/lanes/GLUE/$CASE-$MODEL"
BIN="$PWD/target/lane-glue-$MODEL/release/dft-cli"
mkdir -p "$JOB"
python3 -I scripts/fast_screen_dav.py --upf "$UPF" --case "$CASE" \
  --threads 1 --seed 601 --amp 0.001 --emit-input "$JOB/input.toml"
git rev-parse HEAD > "$JOB/source.txt"
sha256sum "$BIN" "$UPF" "$JOB/input.toml" > "$JOB/identities.txt"
env | sort | sed -n '/^BAANDI_/p' > "$JOB/inherited-switches.txt"
uptime > "$JOB/load-before.txt"
launcher=(taskset -c "$CPUSET")
if [ "$NP" -gt 1 ]; then
  launcher+=(mpirun -np "$NP" --map-by core --bind-to core --report-bindings)
fi
for MODE in off counts timed1 timed2 timed3; do
  GLUE=on CLOCKS=on FFT=on
  [ "$MODE" != off ] || { GLUE=off; FFT=off; }
  [ "$MODE" != counts ] || CLOCKS=off
  set +e
  env BAANDI_GLUE_PROFILE="$GLUE" BAANDI_PHASE_PROFILE=on \
    BAANDI_STAGE_TIMERS="$CLOCKS" BAANDI_FFT_PROFILE="$FFT" \
    BAANDI_PHASE_WAIT=off BAANDI_MIX_COULOMB_GUARD=off BAANDI_STOP_PREDICTOR=off \
    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    "${launcher[@]}" "$BIN" run "$JOB/input.toml" --output-dir "$JOB/$MODE" \
    > "$JOB/$MODE.stdout" 2> "$JOB/$MODE.stderr"
  RC=$?
  set -e
  echo "$MODE $RC" >> "$JOB/exit-codes.txt"
  [ "$RC" -eq 0 ] || exit "$RC"
done
uptime > "$JOB/load-after.txt"
python3 -I scripts/lane_glue_parse.py "$JOB" \
  --json-out "$JOB/glue-tables.json" --markdown-out "$JOB/glue-tables.md"
```

These are complete SCFs to expose coarse and verification work, three timed
profile replicates and two correctness controls. They are not an optimization
A/B or adoption sample. Do not compare performance across the two models.
No A/A floor, expected gain or adoption threshold is fabricated in Phase 0.
Phase 2/3 will pre-register the selected counter and require resolved gain
`abs(mean)/SE >= 3` under the same-model crossed-pair protocol.

## Required returned evidence and validity gates

Return root `$JOB/$MODE/results.json` and `rank1/results.json` through
`rank7/results.json`, plus every density hash, command/exit code, metadata,
load/binding, stdout/stderr and parser outputs. Keep originals on dev08.
All files must have `run.completed=true`, convergence and verification true.
Check the expected rank set 0..NP-1 exists for every mode, rather than using
only rank 0. Parser accepts individual rank files; completeness of the job
is a relay gate, not inferred from run IDs (which are local per process).

Pre-registered correctness gates: **0 Ha energy delta, identical binary64
scientific values/history/eigenvalues/occupations and identical density bytes**
between off/counts/timed on the same model; unchanged SCF and solver work.
Counts and timed new work/call counters must match exactly. Off counters
are all zero; count-only new `.ns` counters all zero. MPI must exercise
`grid.potential`, `grid.audit`, `dist.mix` with positive calls; one core must
exercise `serial.mix`. Missing/zero coverage is unusable, not an improvement.

Reject the profile as incomplete if a required rank/file/counter is absent,
any run fails, any science/work changes, or clocks appear in new count-only
fields. Existing engine wall clocks can remain nonzero in count-only mode;
they are not the new GLUE timers. Return repeated per-rank values and spread,
never collapse overlapping stages or sum ranks. The largest remaining term
and at least two candidates are decisions for Phase 2 after this relay.
