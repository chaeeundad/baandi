# LANE-INIT — Phase 0 profile request

Status: **PROFILE_REQUEST**. The user's revised instruction of 2026-11-20
supersedes the original lane charter in this worktree. No optimization or
default change is proposed. Stop after the master returns the dev08 profile;
design belongs to Phase 2. Mac runs establish correctness and parser behavior
only. Earlier exploratory Mac measurements are not performance evidence and
are not used to select or justify a change.

Base: `db25e70176bd5d111dab6fe1dc2a110766fcb005`; changes uncommitted for
orchestrator review. The master records the resulting committed source and
binary hashes. Do not time a different source tree under this base hash.

## Scope and instrumentation contract

`BAANDI_INIT_PROFILE=on` (input `[profile] init_profile="on"`) enables the
new run-local INIT counters. Default **off**. `BAANDI_STAGE_TIMERS=off`
retains work/call counters with no new clocks; `on` enables nanoseconds.
`BAANDI_PHASE_PROFILE=on` additionally supplies existing atomic-start,
Ewald, construction and session-setup envelopes and launcher detail.
These switches are independent; INIT_PROFILE alone does not turn on the
global phase profiler. No numerical operations, solver control or defaults
were changed. No new I/O, allocation or collective is performed by the
instrumentation; disabled hooks only check the flag/optional token, without
clocks, metric updates, or diagnostic size traversals.

New values are **local to each process/rank** in `run.counters`. Call counts
include entered operations that fail; work counts are recorded after the
corresponding successful stage. Timers use RAII even on errors. Values are
cumulative across all constructors (full cutoff and coarse ladder included),
not a trace of individual calls. The parser reports inclusive seconds/call.
No additional per-G/atom/radial-sample timer or collective perturbs the loops.
Counters within worker tasks use the run scope; requested jobs use one worker
per rank. No cross-parallel-model comparison is requested.

## Exact jobs (one case and model per job)

On dev08 only, in the reviewed source checkout. Use the master's existing
vendor/MKL/MPI toolchain; install nothing. A profile job is capped at 600 s.
First build (not inside a measured job), retain both build exit codes:

```sh
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/lane-init-one \
  cargo build --release -p dft-cli --features mkl
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/lane-init-mpi \
  cargo build --release -p dft-cli --features mpi,mkl
git rev-parse HEAD
sha256sum target/lane-init-{one,mpi}/release/dft-cli "$UPF"
lscpu -e=CPU,CORE,SOCKET,ONLINE
```

Set `UPF` to the master's existing absolute Si.pz-vbc.UPF path. Run each of
these rows separately; `CPUSET` must contain distinct physical cores on one
socket as verified above. The examples use the existing 0..7 placement;
replace it if that is not the verified allocation, and record the binding.

| CASE | NP | MODEL | CPUSET | binary |
| --- | ---: | --- | --- | --- |
| si64 | 8 | mpi | 0-7 | target/lane-init-mpi/release/dft-cli |
| si32 | 8 | mpi | 0-7 | target/lane-init-mpi/release/dft-cli |
| si16 | 8 | mpi | 0-7 | target/lane-init-mpi/release/dft-cli |
| si32 | 1 | one | 0 | target/lane-init-one/release/dft-cli |
| si64 | 1 | one | 0 | target/lane-init-one/release/dft-cli |

For example, the first independent job is the following bash block; change
only the four row variables for each other job:

```bash
set -euo pipefail
CASE=si64 NP=8 MODEL=mpi CPUSET=0-7
JOB="$PWD/target/lane-init-profile/$CASE-$MODEL"
BIN="$PWD/target/lane-init-$MODEL/release/dft-cli"
mkdir -p "$JOB"
python3 -I scripts/fast_screen_dav.py --upf "$UPF" --case "$CASE" \
  --threads 1 --seed 601 --amp 0.001 --emit-input "$JOB/input.toml"
sha256sum "$BIN" "$UPF" "$JOB/input.toml" > "$JOB/identities.txt"
git rev-parse HEAD >> "$JOB/identities.txt"
uptime > "$JOB/load-before.txt"
launcher=(taskset -c "$CPUSET")
if [ "$NP" -gt 1 ]; then
  launcher+=(mpirun -np "$NP" --map-by core --bind-to core --report-bindings)
fi
for MODE in off counts timed1 timed2 timed3; do
  INIT=on CLOCKS=on
  [ "$MODE" != off ] || INIT=off
  [ "$MODE" != counts ] || CLOCKS=off
  env BAANDI_INIT_PROFILE="$INIT" BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS="$CLOCKS" \
    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    "${launcher[@]}" "$BIN" run "$JOB/input.toml" --output-dir "$JOB/$MODE" \
    > "$JOB/$MODE.stdout" 2> "$JOB/$MODE.stderr"
done
uptime > "$JOB/load-after.txt"
python3 -I scripts/lane_init_parse.py "$JOB" \
  --json-out "$JOB/lane-table.json" --markdown-out "$JOB/lane-table.md"
find "$JOB" -name results.json -o -name density.bin | sort > "$JOB/return-files.txt"
```

Queue wrapper enforces the 600 s deadline; if exceeded, stop the job and
return partial evidence and exit code. Do not run a multi-case batch or
silently discard failed runs. These are complete converging runs so the
cutoff ladder and initialization of later sessions are represented. Three
profile repeats measure reproducibility, not an optimization effect. `off`
and `counts` are correctness controls, not speed comparison arms. Preserve
the default switches of the reviewed commit and report any environment/input
overrides. No A/B adoption threshold or QE speed claim belongs to Phase 0.

## Return files and exact paths

Return all `results.json`, `density.bin` hashes, stdout/stderr, command/exit
codes, source/binary/UPF/input hashes, load, binding and parser JSON/Markdown.
MPI root writes `$JOB/$MODE/results.json`; other ranks write
`$JOB/$MODE/rank1/results.json` through `rank7/results.json`. Do not return
only rank 0. Check numerical histories (excluding their `timings` fields),
all energy components, eigenvalues, verification (excluding `wall_s`) and
density hashes are identical across off/counts/timed within the same model.

For each name below, the exact new JSON paths are
`run.counters["init_lane.<name>.calls"]` and
`run.counters["init_lane.<name>.ns"]`. Nanoseconds are cumulative **inclusive**
times; divide by 1e9 for seconds and by calls for an average call. Do not add
overlapping parent/child spans.

| name | boundary |
| --- | --- |
| preflight | full CLI `preflight::check` |
| structure | periodic structure review, including closest-vector search |
| upf | each preflight UPF read/parse (the later run read remains in launch profile) |
| prepare | `Config::prepare`, including basis/density/grid generation |
| check_basis, check_grid | each per-k validation reconstruction |
| projectors | each concrete eager/owned nonlocal constructor; wrappers do not double count |
| factors | per-species radial/angular table allocation and fill; cache included |
| projector_assembly | per-species coupling, validation, phases and eager packing after factors |
| owned_rows_build | initial selected-row materialization only; excludes later lazy full builds |
| ionic | each concrete local-ionic row constructor |
| phase_tables | per-species table attempt including fallback decision |
| problem | full `ScfProblem` construction, including nested engine/physics stages |
| state_grid | per-k state FFT-grid construction |
| hamiltonian | complex Hamiltonian construction plus worker scratch attachment |
| band_density | per-state band-density workspace construction |
| core_density, start_density | density setup calls, including no-data/disabled checks |
| ionic_fft | initial ionic inverse transform and realness checks, including MPI broadcast |
| worker_ffts | construction of the vector of worker FFT plans (one call per problem) |

Work/size paths use the same `run.counters["init_lane.<suffix>"]` prefix:

| suffix | meaning / floor input |
| --- | --- |
| preflight.wave_rows, preflight.density_rows, preflight.grid_points | summed validated basis sizes; density/grid counted once per preflight |
| factors.rows, factors.shells | G rows and exact squared-norm cache entries, summed over species |
| factors.radial_transforms | shell entries times radial projectors |
| factors.radial_samples | shell entries times sum of projector cutoff mesh lengths; input quadrature samples, not a hardware instruction count |
| factors.harmonics | G rows times magnetic channels; spherical-harmonic calls |
| factors.bytes | logical complex-factor payload written, 16 bytes per coefficient |
| projectors.phase_pairs | atom-row phase evaluations for materialized rows (does not include finite-dot validation) |
| projectors.materialized_bytes | initial packed/owned complex table payload, excluding factor/metadata and lazy builds |
| projectors.owned_calls | number of selected-row constructor calls |
| ionic.rows, ionic.shells | species-row contributions and nonzero-G radial cache entries |
| ionic.radial_samples | nonzero shell entries times retained local radial mesh length; G=0 alpha work excluded |
| ionic.phase_pairs | atom-row structure-factor contributions (table and direct paths both counted) |
| ionic.table_species | successful phase-table species; compare with phase_tables.calls for fallback |
| problem.atoms, problem.bands, problem.kpoints | sum over entered problems (bands after validation) |
| problem.grid_points, problem.wave_rows, problem.density_rows | sum of construction grid / Gamma basis / density sphere sizes |

These are logical work/payload floors, not measured memory traffic or FLOPs.
Use 8/16 bytes times grid points for one real/complex-grid payload. Actual
traffic includes caches, allocator effects, FFT plans and multiple passes;
do not infer effective bandwidth from these floors alone.

Existing detail to return (parser preserves it alongside the new counters):

- `phase_profile.buckets[?bucket=="init"].{min,max}["init.*"]` and
  `["atomic.*"]`, `atomic_start`, `wall`: already reduced min/max across
  ranks. Atomic radial/fill/MGS/Gram/HQ/EVD/rotation and Ewald remain the
  existing instruments. Never reinterpret their maxima as local rank times.
- `phase_profile.buckets[?bucket=="between"].{min,max}["between.session_setup"]`:
  all session setups, not only the first; the parser labels that scope.
- `effective_start[]`, `atomic_start`, `basis`, `bands`: physical shape and
  atomic-start work context; final `atomic_start` may be zero when the fine
  stage starts from prolonged coarse orbitals, so retain `effective_start`.
- `run.counters["owned.*"]`, `["phase_validate.*"]`, `ewald_hits`,
  `ewald_computed`, `problems`: existing lazy/table and validation context.
- `launch.seconds.*`: preflight, second UPF read, MPI launch and construction
  context. Do not combine different ranks' maxima into a synthetic total.

## Phase boundary

After the five profiles and their rank-local tables are relayed, the master
can diagnose dominant terms and request Phase 2. This checkout contains only
instrumentation, its correctness tests, the parser and this request.
No QE source was read, no remote command was executed by this implementer,
and no adoption or speedup is claimed.
