# LANE-VERIFY — Phase 0 PROFILE_REQUEST

No optimization is proposed or implemented. The revised 2026-11-20 instruction
supersedes the old diagnose/optimize charter. Return the dev08 profile before
Phase 2. Mac execution is for compilation, scientific identity and parser
plumbing only: no Mac duration is evidence for a performance conclusion.

## Measurement contract

New switch: `BAANDI_VERIFY_PROFILE=on` / `[profile].verify_profile = "on"`.
Default **off**. It requires `PHASE_PROFILE=on`; all new clocks additionally
require `STAGE_TIMERS=on`. Off performs no new clock reads, record pushes,
size scans or profile gathers (only the cached enable branches remain).
Count-only mode records sizes, event outcomes and work with zero new seconds.
There is no threshold, judgment, schedule, solver, reuse or arithmetic change.

Scope: scalar fixed-occupation session verification, including single-/multi-k
and owned/lazy MPI paths. Separate spin/finite-T drivers are not instrumented.
Every actually evaluated scheduler call and verification attempt has a record
on each rank, identified by event, iteration and cutoff. Coarse/fine sessions
may reuse iteration indices; retain cutoff and record order. Deferred checks
are scheduler events, not failed or passed verification attempts.

Exact path: `phase_profile.verification_lane` (schema `verify-profile-1`);
MPI also aliases it under `parallel.phases.verification_lane`. Use only the
first path, never count both. `records[]` contains:

- `event`: `attempt` or `schedule`; `rank`, `iteration`, `cutoff_ha`, `timed`.
- `status`: attempt -1 error / 0 rejected or clipped / 1 passed; schedule
  0 proceed / 1 defer. Error-attempt attribution is incomplete.
- `seconds.{owned_input,warm_start,map,coulomb,metrics,decision,record,restore,schedule}`:
  consecutive exclusive intervals on the same rank. Attempt `schedule` is zero;
  a scheduler record has only `schedule`. `owned_input` covers the temporary
  owned-density clone; `warm_start` includes threshold/ratchet and `split_warm`;
  `metrics` includes residual guard, density metrics and occupied eigen shifts;
  `decision` includes the unchanged predictor/clip/gap checks and phase filing;
  `record` includes strict informational metrics/certificate construction;
  `restore` includes return/drop and owned-input restoration. An error's
  unclosed tail is filed in `restore` and must not support phase conclusions.
- `map_seconds.{potential,solve,density,energy}`: existing map phase clocks,
  **nested inside** `seconds.map`; do not add them to exclusive intervals.
- `sizes.{grid_x,grid_y,grid_z,full_density_elements,owned_density_elements,
  warm_complex_bytes,k_blocks,lazy_blocks,bands_sum,plane_waves_sum,
  history_records,previous_attempts}`. History/previous attempts apply to
  scheduler events; the other fields apply to attempts. `warm_complex_bytes`
  is the actual full-matrix payload cloned, zero for lazy placeholders, not
  the size of distributed lazy orbitals. Lazy handles are Arc clones, not
  copies of the distributed coefficients. Band/PW sums are global k-block
  sizes, while owned density elements are rank-local.
- `h_vectors`, `davidson_iterations`: actual attempt work, zero for scheduling.

Floor inputs (not measured DRAM traffic): warm copy reads+writes at least twice
`warm_complex_bytes`; owned input copy reads+writes 16 times
`owned_density_elements`; one density-pair read is 16 times local owned
(or full) elements. These do not account for caches, allocation, other passes,
FFT/algebra or communication and are not a speed prediction.

Use existing `phase_profile.buckets[]` (`verify_pass` / `verify_fail`) to inspect
`max`/`min` solve stages and collective calls; existing `verification.attempts[]`
provides threshold, tau, residual bounds and full outcomes, and
`verification.schedule` provides the forecast/actions. Rank lane records are
gathered with two adapter-only allgathers **after** numerical execution; this
communication and record storage are profiling overhead. No new collectives
occur inside verification. Do not use instrumented whole-run wall as adoption
evidence. Preserve the uninstrumented control and report the profiling overhead
separately. No aggregate may sum maxima from different ranks.

## Exact dev08 protocol (master session only; NOT executed on Mac)

Use the approved native MPI+MKL and one-core MKL builds of this reviewed patch,
separate target directories. Record commit, diff/source hashes, compiler,
binary hashes, input/UPF hashes, binding and complete effective switches. Do
not alter geometry/grid/bands/start/threshold/ladder between the three arms.
Set the paths below to the master session's existing approved converging
inputs (max_iterations sufficient to finish; no ten-map screen). Inputs must
request one thread and the standard path with fixed-point verification.
No OpenMP/threads performance work and no QE comparison in this phase.

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
CARGO_TARGET_DIR=target-verify-serial cargo build --release -p dft-cli --features mkl
CARGO_TARGET_DIR=target-verify-mpi cargo build --release -p dft-cli --features mpi,mkl
export VERIFY_SERIAL="$PWD/target-verify-serial/release/dft-cli"
export VERIFY_MPI="$PWD/target-verify-mpi/release/dft-cli"
# Master supplies its existing approved absolute input, UPF and output paths:
: "${VERIFY_INPUTS:?directory containing si16.toml si32.toml si64.toml}"
: "${VERIFY_UPF:?absolute path to the identical Si UPF used by these inputs}"
: "${VERIFY_WORK:?new output root on dev08}"
sha256sum "$VERIFY_SERIAL" "$VERIFY_MPI" "$VERIFY_UPF" "$VERIFY_INPUTS"/si{16,32,64}.toml
```

Run each following command as **one separate job**, on an otherwise idle node,
with the same master-approved physical-core binding for its controls. The
listed CPU masks require CPUs 0–7 to be distinct physical cores on one socket;
verify topology before running, and record any master-required equivalent mask.
First wait for the master's idle/load gate; no overlapping cases or builds.
Each job runs one discarded warmup per arm, then six retained cycles, alternating
`off/counts/clocks` and `clocks/counts/off` order. The control gives the A/A
repeatability floor; no optimization A/B is being attempted.

```sh
python3 -I docs/performance/lanes/VERIFY/dev08_profile.py --binary "$VERIFY_MPI" --input "$VERIFY_INPUTS/si64.toml" --upf "$VERIFY_UPF" --launcher 'taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core' --repeats 6 --out "$VERIFY_WORK/si64-mpi8"
python3 -I docs/performance/lanes/VERIFY/dev08_profile.py --binary "$VERIFY_MPI" --input "$VERIFY_INPUTS/si32.toml" --upf "$VERIFY_UPF" --launcher 'taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core' --repeats 6 --out "$VERIFY_WORK/si32-mpi8"
python3 -I docs/performance/lanes/VERIFY/dev08_profile.py --binary "$VERIFY_MPI" --input "$VERIFY_INPUTS/si16.toml" --upf "$VERIFY_UPF" --launcher 'taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core' --repeats 6 --out "$VERIFY_WORK/si16-mpi8"
python3 -I docs/performance/lanes/VERIFY/dev08_profile.py --binary "$VERIFY_SERIAL" --input "$VERIFY_INPUTS/si32.toml" --upf "$VERIFY_UPF" --launcher 'taskset -c 0' --repeats 6 --out "$VERIFY_WORK/si32-onecore"
python3 -I docs/performance/lanes/VERIFY/dev08_profile.py --binary "$VERIFY_SERIAL" --input "$VERIFY_INPUTS/si64.toml" --upf "$VERIFY_UPF" --launcher 'taskset -c 0' --repeats 6 --out "$VERIFY_WORK/si64-onecore"
```

Input copies retain their contents, so UPF paths in the approved inputs must
be absolute (or otherwise valid relative to the copied input). The runner
refuses non-Linux hosts and existing output directories. It records commands,
external wall, load, environment, status and SHA-256 in each `manifest.json`;
stdout/stderr and all results/densities are retained. Any failed run stops that
job; return failures too. Timed/profiled modes must deliver identical science,
density bytes, map/attempt sequence and work to the control at the same rank
count; any mismatch stops interpretation and is returned for diagnosis.

Parse one job at a time, excluding warmups and controls:

```sh
python3 -I scripts/lane_verify_parse.py "$VERIFY_WORK/si64-mpi8"/[0-5]-clocks/results.json --json-out "$VERIFY_WORK/si64-mpi8/calls.json" --markdown-out "$VERIFY_WORK/si64-mpi8/table.md"
python3 -I scripts/lane_verify_parse.py "$VERIFY_WORK/si64-mpi8"/[0-5]-counts/results.json --json-out "$VERIFY_WORK/si64-mpi8/counts.json" --markdown-out "$VERIFY_WORK/si64-mpi8/counts.md"
# Repeat these two parser commands for each of the other four job directories.
```

Return raw JSON, manifests, parser tables, exact science/density identity checks,
all rank-local size/count records, and per-repeat critical-rank totals/parts.
Keep models separate. Profile interpretation, candidate design and any
optimization implementation are deferred until the orchestrator relays these
results. No speedup threshold or adoption claim applies to Phase 0.

The parser JSON contains `calls` (raw per-call/rank lane rows), `attempts`
(per-event critical rank, including scheduler events distinguished by `event`),
and `runs` (all events accumulated on EACH rank before selecting that file's
critical rank). `runs[].critical_rank_seconds` is a coherent decomposition of
`runs[].max_lane_total_s`; it is not a sum of per-event or per-part maxima.
Local correctness evidence and limitations: [validation.md](validation.md).
