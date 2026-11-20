# LANE-DAVIDSON — Phase 0 profile request

Requested state: **PROFILE_REQUEST**. No optimization or default change.
The user's 2026-11-20 phased instruction supersedes the older lane charter.
Base source: `db25e70176bd5d111dab6fe1dc2a110766fcb005`; changes uncommitted.
Mac output is used **only for correctness and schema verification**. No Mac
clock reading is evidence for a performance conclusion. Earlier exploratory
optimization changes were removed after the scope revision.

## Questions for dev08

Separate per-call cost from call count, per-map work and rank imbalance:

- Measured overlap/reduced Gram rectangles versus factor extension,
  transformed reduced problem, dense EVD/B5 and Ritz coefficients.
- Rotation widths, explicit residual widths and return Gram screening.
- Orthogonalization and solver tree-sum counts, payload sizes and site times.
- Coarse/fine/verification map differences, selected-root widths and solver
  orders. Use the existing restart, fallback and cluster counters alongside
  the new trace, without changing any numerical policy.

No candidates, expected gains or adoption recommendation are proposed in this
phase. REL-SLICES is closed; Advance/bounded-subspace and restart-keep remain
rejected. Historical algebra-gap estimates concern earlier solver paths and
are not substituted for this profile.

## Build and provenance (master session on dev08 only)

Use the approved dev08 toolchain/environment and immutable reviewed source.
No remote installation or configuration changes are requested. From the repo:

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export CARGO_TARGET_DIR="$PWD/target/lane-davidson-profile"
cargo build --release -p dft-cli --features mpi,mkl
export BIN="$CARGO_TARGET_DIR/release/dft-cli"
export UPF=/absolute/path/to/Si.pz-vbc.UPF
export OUT=/absolute/path/to/lane-davidson-profile-results
mkdir -p "$OUT"
git rev-parse HEAD > "$OUT/source-commit.txt"
git diff --binary > "$OUT/source.diff"
sha256sum "$BIN" "$UPF" > "$OUT/hashes.txt"
rustc --version --verbose > "$OUT/rustc.txt"
mpirun --version > "$OUT/mpi.txt"
lscpu > "$OUT/cpu.txt"
rg '__INTEL_MKL' "$MKLROOT/include/mkl_version.h" > "$OUT/mkl-version.txt"
```

`UPF` is the same approved external Si file used by the current master suite;
record its actual absolute path and hash. No UPF is committed. Record compiler,
MKL, MPI versions, CPU/binding and load before/after each job. Build exit 0 is a
prerequisite; record every failed run rather than replacing it silently.

## One case/model per job; at most 600 seconds

Run the following script separately for each point, **not as one queue job**:

```sh
bash docs/performance/lanes/DAVIDSON/profile-job.sh si64 8
bash docs/performance/lanes/DAVIDSON/profile-job.sh si32 8
bash docs/performance/lanes/DAVIDSON/profile-job.sh si16 8
bash docs/performance/lanes/DAVIDSON/profile-job.sh si32 1
bash docs/performance/lanes/DAVIDSON/profile-job.sh si64 1
```

Each point uses the same binary/input for three OFF and three ON runs, crossed
OFF/ON, ON/OFF, OFF/ON. Seed 601, fractional displacement 0.001, one thread per
rank, normal convergence with verification so that coarse and verification
maps are represented. MPI8 is compared only with MPI8; one core only with one
core. No OpenMP/threads experiment and no QE comparison is requested. This is
an instrumentation-overhead control and profile, **not an optimization A/B**.
If a point exceeds the job limit, return the completed repetitions and the
partial/failure logs; split a follow-up by repetition without changing inputs.

## Output contract and interpretation

Switch `BAANDI_DAVIDSON_PROFILE=on` (default off) enables the trace independently
of other timers. Set `PHASE_PROFILE`, `STAGE_TIMERS` and `DAVIDSON_STAGES` on in
both arms to accompany it with the existing exclusive stage tree. OFF creates
no trace buffer, timer reads or additional collectives; hot sites retain only
a thread-local flag branch. This is not a claim of literally zero CPU cycles.

Return all `results.json`, generated input, binary/UPF/input hashes, stdout,
stderr and exit codes. New paths:

- `davidson_profile.schema = "davidson-profile-1"`
- `davidson_profile.fields`, `operations`, `stages`, `kinds`, `reduction_sites`
  describe the compact numeric records (no numeric enum guesses in consumers).
- `davidson_profile.per_rank[R].rank`, `.dropped`, `.records[C]`: each record
  has `call,map,solve,stage,kind,iteration,operation,n,m,k,seconds,detail,tolerance_ha`.
  `map` is a monotonic solve-phase ordinal on that rank; `iteration` is the
  engine SCF iteration. Multiple k/strict-certificate solves may share a map.
  `solve` and `call` are monotonically numbered locally. `tolerance_ha` is the
  actual requested eigen tolerance on solve records, zero on other records. Direct solver APIs
  without engine context use map 0; the requested CLI path has real map labels.
- `dropped` must be zero; bounded capacity is 65536 completed records per rank.
  Overflow invalidates the profile, not the scientific calculation. Recording
  never reallocates inside the solver. Rank gathering occurs only after SCF,
  outside the measured solver, in two length/payload allgathers.
- Timers are **inclusive**. Solve contains all its children; reduced_eigen_step contains full EVD/B5; serial
  `rotate_residual` also contains residual calls; `reduced_gram`/orthogonalize
  contain reductions. Do not sum these rows as exclusive solver wall. Use
  `phase_profile.buckets[*].max/min["solve.*"]` for the existing exclusive tree.

Shape meanings:

| Operation | n, m, k | detail |
| --- | --- | --- |
| solve | local PW rows (global for serial), bands, capacity-per-band | 0 |
| overlap_border / reduced_gram | local PW rows, basis order, new product columns | 1 = combined distributed S/H border; 0 = one rectangle |
| factor_extension / reduced_transform | final order, retained prefix, appended width | 0 |
| ritz_coefficients | order, order, selected columns | 0 |
| return_gram | order, total bands, this rank's checked columns | first checked column |
| reduced_eigen | order, order, output order | 0 = full EVD |
| bordered_update | order, prefix, output order | 0 = attempted B5 update (may decline) |
| reduced_eigen_step | order, prefix, output order | encloses EVD/B5 decision and cache update; 1 = MPI B5 enabled |
| rotate_residual | local PW rows, order, rotated columns | 0 = serial includes residual; 1 = distributed rotation pair only |
| residual | local PW rows, 0, explicit columns | 1 = norm-only; 0 = construction plus norms |
| orthogonalize | PW rows, established basis columns, incoming columns | 0 = serial; 1 = distributed |
| restart | local PW rows, old order, retained bands | 0 = serial contraction including rotation; 1 = MPI contraction copy (rotation is a separate preceding call) |
| tree_sum | payload doubles, ranks, 0 | index in reduction_sites |

The parser supplies GEMM `2mnk` models where exact for the recorded products,
minimum array-access byte models, eigen order-cubed proxy, and reduction payload
bytes. These are **not measured DRAM traffic or network wire bytes**. Factor
extension and orthogonalization preserve raw shape inputs without inventing
exact FLOP counts for their conditional internals. Fused optional kernels can
have different traffic; requested runs retain all defaults. Existing
`convergence.work.gemm_dispatch`, `convergence.work.davidson_lean`,
`parallel.reductions`, `parallel.phases`, and SCF history are returned intact.

## Parser and validity checks

The job script runs:

```sh
python3 -I scripts/lane_davidson_parse.py "$OUT/si64-np8/"on-*/results.json \
  --json-out "$OUT/si64-np8/tables.json" > "$OUT/si64-np8/tables.md"
```

`tables.json[*].calls` retains per-call shape/floor records;
`.map_rank_table` aggregates by map/rank/site; `.map_table` gives per-map rank
minimum/maximum times and counts, never sums of per-rank maxima. The tables
cannot establish a wall saving. Parser rejects missing traces, incomplete run
markers, dropped records, duplicate ranks/calls, malformed/negative/nonfinite
fields and missing solve envelopes. Verify OFF/ON science with the supplied
`verify_identity.py`; all energies, eigenvalues, occupations, scientific SCF
history, verification and density bytes must match exactly. A mismatch stops
profiling and is a correctness issue, not a permissible rounding change.

## Coverage limits

The new detailed trace covers real Gamma serial lean/classic and distributed
workspace paths. The requested standard Gamma inputs take those paths. Complex
k-point algebra and the legacy allocating distributed path do not have complete
new operation traces; their existing counters remain available. MPI orthogonalization has a block-level timer and site-labelled reduction
records, not a per-column local dot-product timer.
No timer is placed inside scalar dot-product loops. Root clusters are observed
through unchanged existing work/results rather than introducing a new policy.

Required next action belongs to the orchestrator: return the dev08 profiles.
This implementer stops at Phase 0. Phase 2 begins only after those profiles.
