# LANE-NONLOCAL — Phase 0 profile request

Status: **PROFILE_REQUEST** — Phase 0 complete; instrumentation only.
The orchestrator approved `real_gamma.rs` and engine `profile.rs` for wiring.
`NonlocalMark` isolates the engine-side attribution; the existing phase JSON
exporter serializes the appended fields without a CLI behavior change.
The user's revised instruction of 2026-11-20 supersedes the older charter in
this checkout. All performance profiling is on dev08 through the master.
Mac execution is only compilation, correctness and parser validation.

## Questions for the returned profile

Separate P^T X, its rank reduction, D times the projections, P times the
coefficients, addition of the image, and the nonlocal energy audit. Determine
which block widths and phases carry the work before proposing any optimization.
The ordinary fine maps and verification must be separate. Cached projections
may legitimately make projection calls fewer than back-product calls. Projection
counters include overlaps requested by either H application or the nonlocal
energy audit. The audit timer itself covers only its scalar energy contraction;
do not treat all projection time as H-application time.
Do not combine this lane's children with their parent H-apply time.

Read prior attempts A43/A58b (traversal changes rejected), NONLOCAL-DGEMM-AB,
P48C and F11-GEMM-DISPATCH-COUNTERS before interpreting a kernel idea.
`BAANDI_MKL_NONLOCAL` is not changed by this patch. The checked-out registry
already says `on`, despite the older ticket's opt-in wording; record its actual
resolved value in every returned manifest rather than silently changing it.

## One case and model per job, under 600 seconds

Use the master-approved dev08 checkout and external pseudopotential. Build once
before profiling. Record the source commit plus patch hash, binary SHA-256,
compiler version, MKL version, input and pseudopotential SHA-256, CPU binding,
load before and after, exit code, verification outcome and all results.json.
No QE execution or cross-model comparison is part of Phase 0.

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
cargo build --release -p dft-cli --features mpi,mkl
sha256sum target/release/dft-cli
```

Set `UPF` to the approved external Si.pz-vbc.UPF and `OUT` to a fresh dev08
scratch directory. The input generator is already in the repository; it makes
the same fractional displacement amplitude 0.001, seed 601, and one worker.
For each job select exactly one `(CASE, NP)`:
`(si64,8)`, `(si32,8)`, `(si16,8)`, `(si32,1)`, `(si64,1)`.
Do not compare MPI8 with one core. No multithreading sweep.

```sh
mkdir -p "$OUT"
python3 -I scripts/fast_screen_dav.py --upf "$UPF" --case "$CASE" \
  --seed 601 --amp 0.001 --threads 1 --emit-input "$OUT/input.toml"
sha256sum "$UPF" "$OUT/input.toml"
```

Run one discarded warm-up, then three repetitions in each of three modes:
`off` (profile off), `counts` (profile on, stage clocks off), `timed` (both on).
Cross the retained order: `off counts timed`, `timed counts off`,
`counts off timed`. All inputs, numerical switches and binaries stay identical.
These are instrumentation overhead controls, not optimization A/B arms.

```sh
# NP=8 (one case per job)
env BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
  mpirun -np 8 --bind-to core --map-by core \
  target/release/dft-cli run "$OUT/input.toml" --output-dir "$OUT/timed-1"
# NP=1 (same single-rank binary, no launcher; fixed core)
env BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
  taskset -c 0 target/release/dft-cli run "$OUT/input.toml" \
  --output-dir "$OUT/timed-1"
```

Replace the two mode variables with `off/off` or `on/off` for the controls;
replace the output suffix by the mode and retained repetition 1..3. Preserve
stdout/stderr and launcher wall separately. Never take an instrumented wall
as an uninstrumented performance result. Stop and relay any scientific mismatch.

The exact crossed-run driver is `profile.sh` in this directory. Each of the
following is a **separate job**, invoked only after the master reserves the
exclusive socket-0 window (no concurrent jobs):

```sh
timeout 600 bash docs/performance/lanes/NONLOCAL/profile.sh si64 8 /data/sig-dft/lanes/NONLOCAL/si64-np8 "$UPF"
timeout 600 bash docs/performance/lanes/NONLOCAL/profile.sh si32 8 /data/sig-dft/lanes/NONLOCAL/si32-np8 "$UPF"
timeout 600 bash docs/performance/lanes/NONLOCAL/profile.sh si16 8 /data/sig-dft/lanes/NONLOCAL/si16-np8 "$UPF"
timeout 600 bash docs/performance/lanes/NONLOCAL/profile.sh si32 1 /data/sig-dft/lanes/NONLOCAL/si32-onecore "$UPF"
timeout 600 bash docs/performance/lanes/NONLOCAL/profile.sh si64 1 /data/sig-dft/lanes/NONLOCAL/si64-onecore "$UPF"
```

## Data contract and parser

The phase table uses `phase_profile.buckets[]`, with bucket name at `.bucket`.
Per-rank extrema of each named value are `.min[NAME]` and `.max[NAME]`.
Read `phase_profile.ranks` and `.stage_timers` as metadata. Under MPI the same
profile is also exposed at `parallel.phases`. The primary schema remains
`phase-profile-2` with additive fields. For example, select the element with
`bucket == "fine"` and read `.max["nonlocal.projection.calls"]`,
`.max["nonlocal.projection.flops"]`, `.max["nonlocal.projection.bytes"]`
and `.max["nonlocal.projection"]`. Verification uses the separate elements
`bucket == "verify_pass"` and `bucket == "verify_fail"`; retries use
`fine_redo`, and the coarse ladder uses `coarse`.
Names:

- `nonlocal.{projection,reduce,coupling,back,add,audit}`: seconds.
- The same prefixes with `.calls`, `.flops`, `.bytes`: integer work counts.
- `nonlocal.rows_sum`, `.projectors_sum`, `.bands_sum`: projection input-size
  sums, with `.width1`, `.width2`, `.width3_8`, `.width9_32`, `.width33_plus`
  counting projection calls in each band-width bin.

Empty-band or zero-projector samples are omitted (no physical nonlocal work);
zero-row ranks still record the reduction payload. GEMM FLOPs are 2 N P B. Logical bytes are 8(NP+NB+PB) for a real product,
not measured DRAM traffic or a claim about backend packing. Addition is NB
adds and 24 NB logical bytes. Reduction bytes are 8 PB payload bytes, not
wire traffic or the network algorithm's total bytes. Coupling counts include
all dense coupling entries, including stored zeros; no sparse-work assumption.
Serial audit arithmetic is 3 times the number of nonzero D entries times B
(the existing audit skips exact zeros); distributed audit arithmetic is
(2 times the number of stored D entries + 2 P) times B. Coupling itself remains
dense in both paths. Byte figures are one-pass logical lower-bound inputs:
cache reloads, packing, output finiteness scans and cached-projection copies
are not counted as additional traffic. The GEMM wrapper's packing and output
scan remain inside its product timer; serial addition includes its output scan.
The cached-projection copy before the back path is outside these child timers.
Do not interpret their sum as a complete H-apply budget.
The existing GEMM dispatch report supplies actual routes; a wide API request
alone is not evidence that MKL ran. Use optional backend counters separately
if attribution of vendor-internal packing is needed.

```sh
python3 -I scripts/lane_nonlocal_parse.py --machine dev08 \
  "$OUT"/timed-*/results.json > "$OUT/nonlocal.csv"
python3 -I scripts/lane_nonlocal_parse.py --machine dev08 --format json \
  "$OUT"/counts-*/results.json > "$OUT/nonlocal-counts.json"
```

Return raw JSON, CSV and all provenance. Zero calls mean an unusable stage,
not a measured zero cost. Counts-only mode must have zero lane seconds.
Per-call and per-map values are **bounds on bucket averages**: min seconds /
max count to max seconds / min count. Extrema can come from different ranks;
do not label their ratio a particular rank's throughput. These are not an
individual-call trace or an individual-map trace. Phase 2 is held until the
master returns the dev08 data.

## Correctness evidence and scope

Mac pure-Rust release, Si8 and Si32, seed 601, fractional displacement 0.001,
one worker: off / counts-only / timed (six actual CLI runs). Four off/on
comparisons have bit-identical energies, eigenvalues, occupations, scientific
history, work, verification attempts and density bytes. Both cells converge
and pass verification. `equivalence.json` records the binary/input evidence
and density/scientific SHA-256 hashes. No timing supports a performance claim.

`mac-counts-results.json` is the **unmodified actual Si8 counts-only CLI
results.json**, retained as a parser fixture. Reproduce validation:

```sh
TMPDIR="$PWD/target/scratch" python3 -I docs/performance/lanes/NONLOCAL/verify_cli.py \
  --binary target/release/dft-cli --upf /path/to/Si.pz-vbc.UPF \
  --out target/scratch/nonlocal-new-correctness
python3 -I scripts/lane_nonlocal_parse.py --machine mac --format json \
  docs/performance/lanes/NONLOCAL/mac-counts-results.json
TMPDIR="$PWD/target/scratch" python3 -I docs/performance/lanes/NONLOCAL/test_parser.py
```

Mac parser mode excludes all timings/throughput. Counts-only mode has zero lane
seconds and the same lane counts as timed mode. One-core `nonlocal.reduce` is
correctly all zero and must remain labeled unusable, not a free measured stage.
The collection is thread-local and qualified for the requested one worker per
rank / one-core models. Do not use these records to analyze multiworker runs
or `HAPPLY_OVERLAP` on pool workers: that needs separate attribution wiring.
Complex k-point kernels and projector construction are not covered by this
real-Gamma application profile. Real MPI/MKL execution remains for dev08;
Mac simulated-rank tests are correctness checks, not MPI performance evidence.
A two-rank Si8 off/on SCF also matches energy bits, density FNV and all
SCF/Davidson/H work; its coarse, fine and verification buckets have nonzero
projection/reduction counts and zero lane clocks in counts-only mode
(`simulated-equivalence.json`).

Final required Mac gate: `scripts/check.sh fast dft-physics dft-engine dft-cli`
exits 1: build/fmt/clippy pass, 480 tests pass, 3 known CLI cache tests fail
because the required TMPDIR is inside this repository, 131 tests are ignored.
This environmental gate limitation is not hidden by the successful scientific
identity checks. The master should run the normal full gate on dev08 before
merging. No permission or scientific threshold was changed to make a test pass.
