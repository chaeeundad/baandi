# R-SCF-DAV: fixed-potential dev08 measurement request

Status: standalone solver and installed bench mode; no driver integration,
no changed defaults, no speed/adoption claim. Run only through the master
session. This request does not authorize server changes or installations.
The new premise is a new block-Davidson organization and stopping contract,
not a retry of REL-SLICES, Advance, restart-keep, PR-CG, UPPER_BORDER,
ROTATE_COLUMNS or packed reductions.

## Local evidence and its limits

Source baseline `fb2dc835111d1eb167511a23c772e58f3f2748ce`; changes uncommitted.
[Evidence JSON](../../../../crates/dft-solvers/tests/ref_davidson_evidence.json)
contains source/binary/input hashes, commands, active bands, per-iteration H
columns/orders, collectives, GEMM shape histograms and comparison errors.
Raw logs and rank records are in this worktree's `target/scratch/` (not portable
UPF bundles). The 12 captured physical inputs were copied read-only from
G-FILTER-REPLAY's existing capture index; no new SCF reference was invented.

Mac correctness rerun with the confirmed full-diagonal preconditioner:
48 fixed-potential solves (Si8, Si16, Si32 and skew; early/mid/late;
simulated ranks 1/2/3/8), threshold 1e-12 Ha, capacity 2n, 200-iteration
validation budget. All converge with equal iteration/H counts across ranks.
Maximum occupied rank difference: 4.41313652e-15 Ha;
maximum all-band rank difference: 1.41245904e-11 Ha
(empty threshold 5e-6 Ha). Against the independent legacy residual solve at
1e-12 Ha, maximum occupied eigenvalue difference is
8.92924623e-13 Ha and maximum occupied principal-angle sine
is 7.24731268e-06. Ritz thresholds are not absolute-error
bounds. Full dense skew and synthetic degenerate-edge comparisons also pass.
The degenerate-edge fixture takes 48 iterations / 237 H columns at every
rank count (both root and replicated reduced placement); the nondegenerate
fixture takes 41 / 233.

| Snapshot | Early iterations / H columns | Middle | Late |
| --- | ---: | ---: | ---: |
| Si8 | 28 / 322 | 8 / 94 | 8 / 89 |
| Si16 | 21 / 532 | 10 / 245 | 4 / 121 |
| Si32 | 30 / 1051 | 9 / 492 | 3 / 208 |
| skew | 10 / 61 | 5 / 31 | 2 / 17 |

The unusually tight validation threshold is not the captured early-SCF
threshold. Three early cases need more than the API's default 20 iterations;
the driver must inspect exhaustion and retry. A historical kinetic-only exploratory Si8 early solve
with **all empty bands also forced to 1e-13 Ha** exhausted 200 iterations
(17/20 flags, 1065 H columns), while occupied error was 5.003e-15 Ha. This
negative result is retained; it is not hidden by changing empty tolerances.
Normal first-SCF thresholds are much looser.

The 679-coordinate skew middle snapshot also has a full coordinate-basis
Hamiltonian reference at ranks 1/2/3/8: maximum occupied difference
6.189493362285248e-14 Ha; principal-angle sine <=3.796142102139293e-7.
Synthetic DCT-rotated dense spectra include two degenerate pairs, one ending
at the occupied edge; compare all requested low roots and occupied projectors
with both dense EVD and the existing solver. No low state is missed in these
fixtures. Root/replicated EVD and profile off/on are bit-identical on the
Si8 late three-rank replay (5 iterations, 62 H columns); repeated calls match.
Exact-start, iteration-limit, dependent-start, rank-local failure, and empty
local-row-rank tests pass. No universal completeness proof is claimed.

The legacy bench default, compiled from HEAD and linked to the same unchanged
legacy libraries, matches the edited example's default bit for bit on four
late snapshots (values, residuals, counts and counters). This is default-path
identity, not an end-to-end SCF timing experiment.

Full gate exit 1: build/fmt/clippy pass; 807 tests pass, 3 fail, 181 ignored. Existing repository-local TMPDIR cache-policy
tests fail (three tests demand a cache outside repositories). Scratch remains
inside the worktree as instructed. See the final report for final exit/counts.
README generation/dashboard publication are outside this ticket's OWNED_FILES;
the master owns those integration updates. The frozen production source and
all legacy solver files are unchanged.

## Bench interface

```
cargo build --release -p dft-cli --example davidson_bench
target/release/examples/davidson_bench --input INPUT.toml \
  --snapshot MANIFEST.json --replay-mode ab --ref-arm new --repeats 3
```

New mode only, requires a captured nonempty start:

- `--ref-arm new|legacy|validate`: new solver; existing solver on the same
  reconstructed H with fresh start H products; or new solver plus an untimed
  strict legacy oracle. This legacy arm is a deliberately normalized fresh-H
  comparator using the captured legacy strict-band/buffer policy (recorded),
  not an exact reproduction of captured kept-H work. New and legacy stopping
  definitions can differ; independent delivered-accuracy checks remain required. The original
  `--snapshot` mode without `--ref-arm` remains the exact engine comparator.
- `--ref-tol HA`, `--ref-iterations N`, `--ref-capacity C`, `--ref-first-step`:
  explicit new-solver controls. Without overrides: captured threshold, 20,
  2n, occupied-count policy (first-step flag false).
- `--ref-root`: root EVD instead of replicated. No slices.
- `--ref-profile`: stage clocks; use only diagnostic arms, not timed acceptance.
- `--ref-hx`: include optional returned H X in the timed solve; default false.
  Diagnostic fresh residual H applications occur after timing/counter capture.
- `--ref-dense`: untimed full-H oracle, limited to 1024 coordinates.
- `--sim-ranks P`: correctness only. A serial ORC snapshot can be repartitioned
  to the requested operator rows; marked `repartitioned`, never exact replay.
  A distributed capture must retain its captured rank count and row hashes.

This replay builder currently supports one worker per rank, LDA-PZ without a core density; it checks
physical operator fingerprint, density bytes, layout, dimensions and row hash.
A rebuilt MPI potential may differ in rounding from the capture; `potential_exact`
is recorded and captured H X is **always discarded**. Both solver arms use the
same rebuilt potential. LDA-PW/PBE/core-density reconstruction is not implied
by this bench; the solver itself is operator-agnostic.

`DAVBENCH_JSON.kinds[0]` contains max-rank solve times, science, and
`counters.reference_record` / `reference_per_rank`. `seconds.ref.*` has stage
seconds; setup and finish are explicit. All serialization/rank-record gathering
is outside the timed solve. Exit 5 after JSON means unconverged bands; never
price it as a successful solve. Basis algebra GEMMs exclude opaque operator
kernels. The existing adapter computes full EVD internally even though only
n roots are returned. Startup/output time is not included in this stage-0
solve metric; full SCF timing belongs to R-SCF-1 integration.

## Build and capture: Si64 MPI8, then Si128 MPI8

Use a reviewed checkout, one physical core per MPI rank, one socket, no
oversubscription. No comparison with threads or different rank counts.
Set the external UPF path in the master session; preserve all hashes.

```sh
: "${RSCF_DAV_UPF:?Set to the verified external Si.pz-vbc.UPF}"
export RSCF_DAV_WORK="$PWD/target/ref-davidson-dev08"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$RSCF_DAV_WORK"
cargo build --release -p dft-cli --features mpi --bin dft-cli --example davidson_bench
sha256sum target/release/dft-cli target/release/examples/davidson_bench "$RSCF_DAV_UPF"
git rev-parse HEAD
git diff --stat
for cell in si64 si128; do
  mkdir -p "$RSCF_DAV_WORK/$cell-np8"
  python3 -I scripts/fast_screen_dav.py --upf "$RSCF_DAV_UPF" --case "$cell" \
    --seed 601 --amp 1e-3 --threads 1 --emit-input "$RSCF_DAV_WORK/$cell.toml"
  uptime
  timeout 7200s env BAANDI_M2_LOCAL_FIELDS=off \
    BAANDI_SNAPSHOT_DIR="$RSCF_DAV_WORK/$cell-np8" \
    mpirun -np 8 --map-by core --bind-to core --report-bindings \
    target/release/dft-cli run "$RSCF_DAV_WORK/$cell.toml" \
    --output-dir "$RSCF_DAV_WORK/$cell-scf" \
    > "$RSCF_DAV_WORK/$cell-capture.log" 2>&1
done
```

Check the binding report is physical cores 0–7 on socket 0; the master must
correct the launcher if it is not. Capture timings include I/O and are not
performance evidence. Do not run captures/builds concurrently with timing.
Select actual early/middle/late fine ordinary warm maps and save their paths:

```sh
python3 - <<'PY'
import os,json,pathlib
work=pathlib.Path(os.environ['RSCF_DAV_WORK'])
for cell in ['si64','si128']:
    root=work/f'{cell}-np8'
    fine=[(p,json.loads(p.read_text())) for p in sorted(root.glob('s*.json'))]
    fine=[(p,v) for p,v in fine if v['stage']=='fine' and v['kind']=='ordinary' and v['start']!='none']
    assert fine and all(v['ranks']==8 for _,v in fine)
    for phase,index in [('early',0),('mid',len(fine)//2),('late',len(fine)-1)]:
        (root/f'{phase}.manifest').write_text(str(fine[index][0].resolve())+'\n')
PY
export BAANDI_M2_LOCAL_FIELDS=off
mpirun -np 8 --map-by core --bind-to core \
  target/release/examples/davidson_bench --input "$RSCF_DAV_WORK/si64.toml" \
  --snapshot "$(cat "$RSCF_DAV_WORK/si64-np8/mid.manifest")" \
  --replay-mode exact --repeats 1 > "$RSCF_DAV_WORK/exact-control.log" 2>&1
```

Require exact-control exit 0; do not bypass an identity rejection. Repeat for
Si128. For each selected snapshot run an untimed `--ref-arm validate` at the
captured threshold and then at `--ref-tol 1e-12 --ref-iterations 200`, recording
eigenvalue errors, angles, convergence, and H work. Gate numerical conclusions
against the independent reference; early loose solves are not final-accuracy
certificates. The 2026-11-23 study-note addendum confirms the smooth preconditioner formula.
The replay now supplies kinetic + local average + nonlocal diagonal in Ha;
matched QE iteration identity remains to be measured on dev08.

## A/A floor, crossed A/B, root placement

Start with Si64 middle; only then Si128 middle and early/late controls. Node
load below 0.30, no concurrent jobs, unchanged bindings/frequency policy.
The snippet below executes one cell/phase; change the two environment values
for each later job. Keep `--ref-profile` absent. Each arm discards its first
in-process warm-up. At most 180 s/process, no silent timeout extension.

```sh
export RSCF_DAV_CELL=si64 RSCF_DAV_PHASE=mid
python3 - <<'PY'
import os,json,pathlib,subprocess,time,statistics,math
work=pathlib.Path(os.environ['RSCF_DAV_WORK']); cell=os.environ['RSCF_DAV_CELL']; phase=os.environ['RSCF_DAV_PHASE']
manifest=(work/f'{cell}-np8/{phase}.manifest').read_text().strip()
base=['mpirun','-np','8','--map-by','core','--bind-to','core',
      'target/release/examples/davidson_bench','--input',str(work/f'{cell}.toml'),
      '--snapshot',manifest,'--replay-mode','ab','--repeats','3']
series=[]
for experiment,pairs in [('AA',6),('AB',12),('placement',12)]:
    diffs=[]
    for pair in range(pairs):
        order=['A','B'] if pair%2==0 else ['B','A']; samples={}
        for arm in order:
            if os.getloadavg()[0]>=0.30: raise SystemExit('busy node: stop and reschedule')
            new=experiment=='placement' or (experiment=='AB' and arm=='B')
            cmd=base+['--ref-arm','new' if new else 'legacy']
            if experiment=='placement' and arm=='B': cmd+=['--ref-root']
            path=work/f'{cell}-{phase}-{experiment}-{pair:02}-{arm}.log'
            with path.open('w') as f:
                result=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=180)
            if result.returncode: raise SystemExit(f'failed {path}: {result.returncode}')
            data=json.loads(next(line[14:] for line in path.read_text().splitlines() if line.startswith('DAVBENCH_JSON ')))
            k=data['kinds'][0]
            assert k['converged']==data['bands'] and k['repeats_identical']
            samples[arm]=k['solve_ms_median']/1000
        diffs.append(samples['B']-samples['A'])
    mean=statistics.mean(diffs); se=statistics.stdev(diffs)/math.sqrt(len(diffs))
    row=dict(experiment=experiment,differences_s=diffs,mean_s=mean,se_s=se,
             abs_mean_over_se=abs(mean)/se if se else None)
    series.append(row);print(json.dumps(row),flush=True)
    if experiment=='AA' and se and abs(mean)/se>=3: raise SystemExit('biased AA window: stop')
(work/f'{cell}-{phase}-paired.json').write_text(json.dumps(series,indent=2)+'\n')
PY
```

A resolved gain needs negative paired mean, `|mean|/SE >= 3`, effect at least
three times the measured A/A floor, and no correctness failure. Report the
floor definition explicitly (use `max(abs(AA mean), AA SE)`) and intervals;
otherwise mark inconclusive. Do not infer speed from removed FLOPs. For root
placement choose the faster validated arm only on dev08 evidence; this ticket
does not flip a production default. Repeat a separate diagnostic pair with
`--ref-profile` (and, if needed, `--ref-hx`) to attribute stage cost and compare
H blocks/active bands/logical reduced calls with QE's same-scope counters.
Do not compare the all-stage legacy 183 solves against fine-only 58 iterations.

If MKL is required, repeat as a separate backend experiment with
`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `--features mpi,mkl`, a separate
`CARGO_TARGET_DIR` and new binary/capture hashes. Never pool pure/MKL samples.

Return exact commands/exits, revision plus patch/source hashes, compiler/MPI
versions, binary/UPF/input/capture hashes, bindings/load, all raw JSON/logs,
all eight rank records, A/A and crossed differences/SE, numerical errors,
active bands and H widths per iteration, EVD orders/calls, repairs/restarts,
logical collective bytes, GEMM shapes/flops, and separate diagnostic stage
seconds. Retain snapshots on dev08; do not return or commit UPF contents.
