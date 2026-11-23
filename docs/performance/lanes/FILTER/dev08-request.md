# FILTER stage 1: optional falsification of the STOP verdict

Status: **do not integrate the filter into SCF**. Mac count-based modeling is
already negative optimistically (best fixed policy -9.8693 s at Si64 MPI8).
This request preserves exact capture/replay commands and a bounded way to
price the missing terms if the master wants independent confirmation. It is
not authorization for a full candidate SCF, new default, server modification,
installation, or a multi-case timing grid.

Two prerequisites: read GitHub issue #24 (unavailable in the implementer sandbox)
and resolve the ticket's mistaken owned bench path. The actual example is
`crates/dft-cli/examples/davidson_bench.rs`; `bench.patch` is prepared/tested but
not installed. The master can apply it after ownership correction. The scratch
builder is for local review only; use the actual Cargo example for MPI/MKL so
native link dependencies are resolved by Cargo.

## Build and capture job (one case, one rank count, <=600 s)

Run from the reviewed checkout on dev08, through the master session. Keep the
UPF outside the repository. `FILTER_UPF` must point to the same Si UPF as the
frozen benchmark; record its SHA-256 and compare it with capture-index.json and
the reference job. The following commands use their own output directory.

```sh
: "${FILTER_UPF:?Set FILTER_UPF to the externally stored verified Si.pz-vbc.UPF}"
export FILTER_WORK="$PWD/target/filter-replay-dev08"
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$FILTER_WORK/si64-np8"
cargo build --release -p dft-cli --features mpi,mkl --bin dft-cli --example davidson_bench
sha256sum target/release/dft-cli target/release/examples/davidson_bench "$FILTER_UPF"
python3 -I scripts/fast_screen_dav.py --upf "$FILTER_UPF" --case si64 \
  --seed 601 --amp 1e-3 --threads 1 --emit-input "$FILTER_WORK/si64.toml"
uptime
timeout 590s env BAANDI_SNAPSHOT_DIR="$FILTER_WORK/si64-np8" \
  mpirun -np 8 --bind-to core --map-by core --report-bindings \
  target/release/dft-cli run "$FILTER_WORK/si64.toml" \
  --output-dir "$FILTER_WORK/si64-np8-scf" \
  > "$FILTER_WORK/si64-np8-capture.log" 2>&1
```

Require physical-core bindings 0–7 on socket 0 in the recorded binding report;
otherwise reject the run and fix the launcher in the master, without using its
timing. Node idle gate: load below 0.30 before each timed process, unchanged
frequency policy, no builds or other jobs overlapping timing. Snapshot capture
wall includes I/O and is never used as a performance result.

Select actual early/mid/late fine maps and their frequencies, rather than
assuming the serial Mac sequence numbers:

```sh
python3 - "$FILTER_WORK/si64-np8" <<'PY'
import hashlib,json,pathlib,sys
root=pathlib.Path(sys.argv[1])
rows=[(p,json.loads(p.read_text())) for p in sorted(root.glob('s*.json'))]
fine=[(p,v) for p,v in rows if v['stage']=='fine' and v['kind']=='ordinary']
assert fine and all(v['ranks']==8 for _,v in fine)
n=len(fine); ends=[0,(n+2)//3,(2*n+2)//3,n]; selected=[]
for i,label in enumerate(['early','mid','late']):
    group=fine[ends[i]:ends[i+1]]
    p,v=group[0] if i==0 else group[-1] if i==2 else group[len(group)//2]
    (root/f'{label}.manifest').write_text(str(p.resolve())+'\n')
    selected.append(dict(phase=label,path=str(p.resolve()),map=v['iteration'],
        weight_actual=len(group),weight_11=11*len(group)/n,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(root/'selected.json').write_text(json.dumps(dict(fine_maps=n,samples=selected),indent=2)+'\n')
print(json.dumps(selected,indent=2))
PY
```

Retain every manifest, rank binary, input-density binary, input TOML, SCF result,
stdout/stderr and SHA-256 manifest on dev08. Return hashes, not the UPF. All
captures should be at the default solver configuration; the capture switch
itself disables the owned-local-potential optimization while serializing.

## Stage 0 job: Si64 np8, middle snapshot, degree 4 only

This is one case/phase/degree and one rank count. Do A/A first, then three
crossed A/B pairs. It must finish within 600 s; the enclosing queue job also
has its hard 600 s timeout. Stop if any single process exceeds its allocated
90 s below; return partial results and split at pair boundaries. Never silently
extend the timeout. All commands use the same binary, input and snapshot.

`BAANDI_M2_LOCAL_FIELDS=off` is required for replay of the captured full potential;
without it the current exact replay rejects the owned/sentinel-potential hash.
Do not switch to legacy replay or ignore an identity failure.

```sh
export BAANDI_M2_LOCAL_FIELDS=off
export FILTER_MANIFEST="$(cat "$FILTER_WORK/si64-np8/mid.manifest")"
mkdir -p "$FILTER_WORK/si64-np8-mid-d4"
# A is exact captured Davidson; B additionally executes the filter and fallback.
filter_a() {
  timeout 90s mpirun -np 8 --bind-to core --map-by core \
    target/release/examples/davidson_bench --input "$FILTER_WORK/si64.toml" \
    --snapshot "$FILTER_MANIFEST" --replay-mode exact --repeats 3
}
filter_b() {
  timeout 90s mpirun -np 8 --bind-to core --map-by core \
    target/release/examples/davidson_bench --input "$FILTER_WORK/si64.toml" \
    --snapshot "$FILTER_MANIFEST" --replay-mode exact --repeats 1 \
    --filter-degrees 4
}
# Job AA: all six runs are A, labelled in order to estimate paired noise.
for pair in 1 2 3; do
  filter_a > "$FILTER_WORK/si64-np8-mid-d4/AA${pair}a.log" 2>&1 || exit $?
  filter_a > "$FILTER_WORK/si64-np8-mid-d4/AA${pair}b.log" 2>&1 || exit $?
done
```

In a **separate <=600 s job**, after validating the A/A floor and stability:

```sh
filter_a > "$FILTER_WORK/si64-np8-mid-d4/AB1-A.log" 2>&1 || exit $?
filter_b > "$FILTER_WORK/si64-np8-mid-d4/AB1-B.log" 2>&1 || exit $?
filter_b > "$FILTER_WORK/si64-np8-mid-d4/AB2-B.log" 2>&1 || exit $?
filter_a > "$FILTER_WORK/si64-np8-mid-d4/AB2-A.log" 2>&1 || exit $?
filter_a > "$FILTER_WORK/si64-np8-mid-d4/AB3-A.log" 2>&1 || exit $?
filter_b > "$FILTER_WORK/si64-np8-mid-d4/AB3-B.log" 2>&1 || exit $?
```

Functions/exports above must be sourced again in a new shell. A metric is
`kinds[0].solve_ms_median / 1000`. B metric is maximum over `per_rank` of
`candidates[0].observation_seconds.total`, **including fallback**. This prototype
allocates its own workspace per solve; do not subtract that cost. The independent
strict oracle, operator reconstruction, density/angle diagnostics and exact
baseline repetitions are outside B's timer and are outside A's solver timer as
well. Original engine replay identity must still pass before B is considered.
Because the reconstructed MPI operator can differ in potential rounding, even
a resolved timing gain is only a lead, not an exact-replay or adoption result.

If the instrumented branch is used, never mix it into these timed pairs.
For A/A require stable work/eigen FNV and no resolved order drift; define the
floor from the paired differences (report samples, mean, SE and spread). For
A/B require |mean(B-A)|/SE >=3 and an effect above that floor. A count failure,
nonfinite metric, timeout, identity exit 3/4, inaccurate result, or a resolved
loss closes the job. A negative optimistic full-cost model closes integration
regardless of a local RR-only gain.

## Separate price/counter job (same single snapshot)

```sh
timeout 590s env BAANDI_M2_LOCAL_FIELDS=off \
  mpirun -np 8 --bind-to core --map-by core \
  target/release/examples/davidson_bench --input "$FILTER_WORK/si64.toml" \
  --snapshot "$FILTER_MANIFEST" --replay-mode exact --repeats 1 \
  --filter-degrees 4 --filter-profile \
  > "$FILTER_WORK/si64-np8-mid-d4/profile.log" 2>&1
```

`profile_seconds_per_rank` gives disjoint setup, bounds, polynomial recurrence,
QR, terminal application, and RR phase totals. H time by exact block width
(including fallback) is **inclusive** within these phases/fallback: never add it
a second time. `total - sum(phases) - fallback` is the unclassified solver
allocation/scalar/control remainder; report it, do not erase it. Report max-rank
phase times as diagnostics only; maxima of individual phases do not sum to one
real rank's wall. Use each rank's own closure before computing maxima. Width-2
bound work, candidate-width H, QR, initial RR/setup, fallback, terminal H and
copies must all receive nonzero count or an explicit reason for absence.

## What to return

One JSON object per job plus all raw `DAVBENCH_JSON` lines:

- commit, source/patch and binary SHA-256, compiler/MKL/MPI versions, build/run
  exit codes, exact command/env, input and UPF hashes, before/after load, binding
  report and elapsed budget;
- capture index: actual map frequencies, stage/kind/map, stopping profile,
  operator fingerprint, potential/density/row hashes, rank-file hashes;
- `snapshot.replay_mode` and `snapshot.captured`; compare captured Davidson/H
  counts and eigen FNV with `kinds[0]`. The reporter prints the work/eigen
  identity booleans in text and exits 4 on mismatch; those booleans are not
  separate JSON fields;
- `snapshot.filter_replay.per_rank[*]`: rebuilt-baseline exact flag,
  potential exact flag and HX difference; do not replace these with the engine
  exact-replay verdict;
- baseline H-width histogram, counters (Gram/rotate flops, EVD calls and sum
  m^3), eigenvalues/residuals; candidate degree/bounds, every attempted sweep,
  fallback reason and executed work, filter H-width histogram, charged histogram,
  bound/terminal columns, QR/Gram/rotate/recurrence flops, RR orders, reduction
  counts and doubles by operation, copied bytes;
- principal angles (rad), occupied eigenvalue error (Ha), projector spectral
  error, density max/L2 errors with units, frozen-band-energy error (Ha), both
  absolute gate flags, strict reference residual. No total-SCF-energy assertion
  follows from these numbers; early maps can fail the band-sum gate even on A;
- all A/A and crossed samples, mean/SE/floor, max-rank full candidate time,
  per-rank optional phase and H-width prices, fallback/terminal cost and closure;
- updated 11-map cost model with actual MPI8 counts and all incremental terms,
  including nonzero accuracy-repair cost if repair is proposed. A no-repair
  inaccurate candidate is a failure, not a faster accepted calculation.

The fallback EVD full order histogram is not instrumented in this patch; only
its count and sum m^3 are available. Do not fabricate it. The current prototype
uses an estimated upper bound and is limited to LDA-PZ without core density.
An exact MPI frozen-operator hook and complete fallback histogram would need a
new, explicitly owned ticket if the cost evidence ever justified reopening.

## Expansion, only after an evidence-supported reason to overturn STOP

One job per phase/degree/rank count. Early and late: substitute their `.manifest`
files in the same commands; degrees 8/12/16/auto are separate jobs. Rank counts
1/2/4/8 require **fresh default captures at that rank count**, then exact replay
at the same count. Substitute `P` in the capture directory and `mpirun -np P`;
select manifests from that directory with the selector's rank assertion updated
to P. For one core omit mpirun and never compare its wall to any MPI wall.
Si32/Si16/skew use their existing benchmark TOMLs and fresh same-rank captures.
Never transplant an np1 rank file into np8 or compare unlike parallel models.
No expansion is presently recommended.
