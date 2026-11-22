# G-FIXTURE-SI128: larger silicon Gamma fixtures

## Scope and geometry

Added `si128-gamma` (2x2x4 conventional cells) and `si216-gamma` (3x3x3),
with exactly the Si64 field set: real potential `si`, 30 Ry wavefunction
cutoff, Gamma, deterministic, `phase = "scf"`, `tier = "long"`, and a
7200 s per-run timeout. Existing fixtures and numerical code are unchanged.
Si512 (4x4x4) is **documentation/preflight only**: no benchmark catalog file,
no Mac SCF, and no authorization to run it before dev08 memory qualification.

The conventional eight fractional sites come from `si8-gamma.toml`:
`(0,0,0), (0,1/2,1/2), (1/2,0,1/2), (1/2,1/2,0),
(1/4,1/4,1/4), (1/4,3/4,3/4), (3/4,1/4,3/4), (3/4,3/4,1/4)`.
Each site is mapped to `(translation + site) / replication`, with translations
in lexicographic x/y/z order, z fastest, followed by the eight sites in their
existing order. Atom order need not match Si64; the regenerated Si64 **set and
lattice match exactly**, so its ordering is left untouched.

The requested conventional edge is 10.261212856717932 Bohr. To preserve the
existing binary64 lattice, generation uses `si64.lattice[0][0] / 2`, whose
short decimal rendering is 10.261212856717933 Bohr. This is the existing
roundoff convention, not a lattice adjustment.

**Ticket distance correction:** this edge gives `sqrt(3)*a/4 =
4.443235503778610 Bohr = 2.351258971274751 Angstrom`, not 4.6 Bohr (2.43
Angstrom). Enforcing 4.6 would change the required Si64 geometry. Periodic
minimum-image pair checks verified all target counts, uniqueness, fractional
coordinates in [0,1), no shorter distances, and exactly four nearest neighbors
per atom. For Si216 the minimum is 4.443235503778607 Bohr (floating roundoff).
No distance tolerance was relaxed to admit a wrong structure: the 1e-12 Bohr
comparison tolerance applies to the analytic diamond bond.

## Engine-derived sizes (Mac preflight, 2026-11-22)

Source HEAD: `58597f9c7b0c8b350aac12df551ee3c6b6bd9624` plus these fixture/docs
changes. `scripts/bench_suite.py::load_case` and `case_input` parse and render
the flat benchmark schema; `target/debug/dft-cli validate INPUT --json` and
`grid INPUT` execute the existing Rust basis/FFT code, **without SCF**.
Each of Si64, Si128, Si216 and the scratch-only Si512 input passed both commands
(exit 0). The fine wavefunction cutoff is 30 Ry = 15 Ha; the density sphere
is at 120 Ry = 60 Ha. The FFT dimensions are automatic, without an override.
"Fine" here means the final 30 Ry stage, not any coarse cutoff-ladder stage.

External read-only `../pseudo/Si.pz-vbc.UPF` has SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Its valence is 4; its radial s and p projectors expand into 1+3 = 4 angular
columns per atom. No potential bytes are added. Defaults give `bands =
electrons/2 + 4`; these band counts must be explicit in the QE arms.

| Case | Diagonal lattice edges (Bohr) | Atoms | Electrons | Bands | Full PW at 30 Ry | Density PW at 120 Ry | Fine FFT grid |
|---|---|---:|---:|---:|---:|---:|---|
| Si64 | 20.522425713435865 / 20.522425713435865 / 20.522425713435865 | 64 | 256 | 132 | 23,871 | 191,735 | 72 x 72 x 72 |
| Si128 | 20.522425713435865 / 20.522425713435865 / 41.04485142687173 | 128 | 512 | 260 | 47,851 | 383,423 | 72 x 72 x 144 |
| Si216 | 30.783638570153798 / 30.783638570153798 / 30.783638570153798 | 216 | 864 | 436 | 80,797 | 646,963 | 108 x 108 x 108 |
| Si512 (documented only) | 41.04485142687173 / 41.04485142687173 / 41.04485142687173 | 512 | 2048 | 1028 | 191,735 | 1,533,851 | 144 x 144 x 144 |

The real-Gamma half sphere has `(G+1)/2` representatives; the ORC vector
still has `G` real coordinates. Do not halve the real payload a second time.

## Memory estimates (payloads, not measured peak RSS)

Let `G` be the full PW count, `B` bands, `Q = 4*atoms` angular projectors,
and `p` MPI ranks. One real-Gamma wavefunction block is `W = 8*G*B` bytes;
one real projector table is `P = 8*G*Q` bytes. Below are the **mean per-rank
payloads** `W/p` and `P/p` (MiB = 2^20 bytes), assuming the M2 owned-row layout.
The actual stick partition is uneven; these are not maximum-rank bounds.
`distributed_real_gamma.rs` partitions projector rows with the ORC rows.

The final two columns give a **full complex replica on each rank**, `2W` and
`2P`. These do not divide by p. Owned/lazy paths can avoid these replicas;
materialization for output, hand-off or fallback can bring them back. A sparse
complex representative projector table adds approximately `P/p` to the real
`P/p` while both are retained (16 bytes per representative versus 8 per ORC
row). Record the actual active path; do not add every alternative as if it
were always simultaneously live.

| Case | np | One real wavefunction block / rank (MiB, mean) | Real projectors / rank (MiB, mean) | Full complex wavefunction replica / rank (MiB) | Full complex projector replica / rank (MiB) |
|---|---:|---:|---:|---:|---:|
| Si64 | 1 | 24.040 | 46.623 | 48.080 | 93.246 |
| Si64 | 8 | 3.005 | 5.828 | 48.080 | 93.246 |
| Si64 | 16 | 1.503 | 2.914 | 48.080 | 93.246 |
| Si128 | 1 | 94.919 | 186.918 | 189.839 | 373.836 |
| Si128 | 8 | 11.865 | 23.365 | 189.839 | 373.836 |
| Si128 | 16 | 5.932 | 11.682 | 189.839 | 373.836 |
| Si216 | 1 | 268.764 | 532.597 | 537.529 | 1065.195 |
| Si216 | 8 | 33.596 | 66.575 | 537.529 | 1065.195 |
| Si216 | 16 | 16.798 | 33.287 | 537.529 | 1065.195 |
| Si512 | 1 | 1503.781 | 2995.859 | 3007.562 | 5991.719 |
| Si512 | 8 | 187.973 | 374.482 | 3007.562 | 5991.719 |
| Si512 | 16 | 93.986 | 187.241 | 3007.562 | 5991.719 |

These two blocks are only a small part of a solve. Davidson V/HV bases,
rotations, residuals, eigensolver and GEMM workspaces, atomic-start arrays,
FFT batches/transposes, mirror exchanges, densities/mixing history, coarse/fine
handoff overlap, MPI buffers and allocator overhead are additional. Dense
subspace matrices are replicated: at capacity `m=6B`, **each** matrix takes
`8*m*m` bytes per rank. CLI `estimated_array_bytes` counts only a complex
orbital block plus one real and one complex FFT grid; it is not a memory gate.

A concrete risk scenario is the serial `RealDavidsonWorkspace::new` in
`crates/dft-solvers/src/real_davidson.rs`: six `G x m` real arrays plus four
`G x B` arrays already total `40W` at m=6B, before its inner workspaces.
Adding one real projector table, a retained full complex table and one full
complex orbital block gives the scenario `42W + 3P`:

| Case | Selected-array scenario (GiB, excludes other workspaces) |
|---|---:|
| Si64 | 1.123 |
| Si128 | 4.441 |
| Si216 | 12.584 |
| Si512 | 70.455 |

This is a documented retained-array scenario, **not a prediction of the
current distributed or lean solver's RSS**, and not a guaranteed upper bound.
For Si512, even `40W` alone is 58.741 GiB; the scenario is 70.455 GiB
(75.651 GB decimal). Its single real orbital/table blocks are 1.469/2.926 GiB.
A full complex projector replica alone is 5.851 GiB per rank, or 46.810/93.621
GiB across np8/np16. Thus "64 GB node" does not establish feasibility.
Keep Si512 out of the runnable catalog. After Si128/Si216, the master must
verify actual MemTotal/available RAM, maximum per-rank and aggregate job RSS,
all simultaneous retained arrays and output/verification peaks on the selected
MPI path. Require aggregate peak below 75% of actual MemTotal and no swap/OOM,
with demonstrated room for the next size; preflight alone cannot qualify it.

## Operation-count scaling relative to Si64

Fixed SCF/Davidson iteration counts, fixed capacity/band multiplier, same
rank count and same kernel choices are assumed. With F FFT points: FFT work
scales as `B*F*log2(F)`, projector products as `G*Q*B`, Gram/rotation as
`G*B^2`, and the replicated subspace eigensolve as `B^3`. These are operation
ratios, not wall predictions or BAANDI/QE ratios. Communication, cache/NUMA,
coarse stages, verification and changing convergence counts are excluded.

| Case | FFT B F log2(F) | Nonlocal G Q B | Gram/rotation G B^2 | Reduced eigensolve B^3 |
|---|---:|---:|---:|---:|
| Si64 | 1.000 | 1.000 | 1.000 | 1.000 |
| Si128 | 4.152 | 7.897 | 7.777 | 7.642 |
| Si216 | 12.205 | 37.732 | 36.927 | 36.036 |
| Si512 | 72.401 | 500.426 | 487.157 | 472.343 |

## dev08 measurement plan (master session only; not run here)

Use the reviewed commit and the existing MPI/MKL toolchain, no installation or
server changes. Build outside timed jobs and record build exit code, source
commit/dirty diff, compiler and MPI/MKL versions, binary SHA-256 and UPF SHA-256.
For example, with the master's existing toolchain environment:

```sh
MKLROOT=/opt/intel/oneapi/mkl/2026.1 CARGO_TARGET_DIR=target/fixtures-si128-mpi \
  cargo build --release -p dft-cli --features mpi,mkl
```

The benchmark catalog is not CLI input. Render each case with the existing
`case_input(case, absolute_upf_path, 7, 1)` as in the reproducer below, then
validate it on dev08. The old fast-screen scripts' case enumerations need not
support these new IDs; do not substitute a smaller built-in case. Explicitly
record the effective standard threshold (default 5e-9 Ha, equivalent to QE
`conv_thr=1e-8` Ry), `verify="fixed_point"`, `nbands`, atomic start, Gamma-real
path, cutoff-ladder settings and all performance overrides. No default or
threshold is changed by this ticket. The timed total includes verification;
also return its time separately. Strict reference checks are separate runs,
not substituted for the standard headline arm.

### Separate MPI-vs-MPI conditions

| Case | np | BAANDI bands = QE nbnd | BAANDI fine grid = QE nr1,nr2,nr3 and nr1s,nr2s,nr3s | Status |
|---|---:|---:|---|---|
| si64-gamma | 8 | 132 | 72,72,72 | Same-window anchor |
| si64-gamma | 16 | 132 | 72,72,72 | Same-window anchor |
| si128-gamma | 8 | 260 | 72,72,144 | New scaling point |
| si128-gamma | 16 | 260 | 72,72,144 | New scaling point |
| si216-gamma | 8 | 436 | 108,108,108 | After Si128 memory/time screen |
| si216-gamma | 16 | 436 | 108,108,108 | After Si128 memory/time screen |

Each row has **BAANDI MPI against QE MPI at the same np**, all library/OpenMP
threads 1 (`OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`).
Use the same distinct physical cores on one socket, verify actual per-rank CPU
binding with the master's launcher and `--report-bindings`/affinity records;
np8 uses eight cores and np16 sixteen. No SMT siblings, competing job or build.
Do not infer affinity from `taskset` alone. Record load, delivered clocks,
MemTotal, MemAvailable, swap and topology before/after each job. np1 appears
in the sizing table only; it is not required as a large-cell timing arm.

For QE 7.6 use the master's fastest qualified build/configuration (P0.3),
`calculation='scf'`, `ibrav=0`, `nat=128/216`, `ntyp=1`, explicit `nbnd`,
`ecutwfc=30`, `ecutrho=120`, `occupations='fixed'`, `input_dft='PZ'`,
`startingpot='atomic'`, `startingwfc='atomic'`, `conv_thr=1e-8`, and
`K_POINTS gamma`. Emit `CELL_PARAMETERS bohr` and `ATOMIC_POSITIONS crystal`
from the exact same rendered geometry, and use the same UPF bytes. Set both
QE dense and smooth FFT grids explicitly to the engine values above; verify
the *reported* grids, electron/band counts and final cutoff in both outputs.
Do not inherit Si64's 72^3 grid for the larger cells. If alternate QE
parallelization settings are tried at the same np, retain every arm and select
its fastest converged/accuracy-qualified configuration before the final
paired window. Publish chosen flags and excluded-arm reasons.

The pristine catalog geometry is the first size/feasibility series. Label its
results "unperturbed size scaling"; Si128 is elongated while Si216 is cubic,
so aspect ratio and FFT smooth dimensions also change. For a later headline
series use the master's approved realistic displacement amplitude and common
geometry seeds in *both* engines; record that amplitude and exact coordinate
hashes. Keep the ten-geometry seeds 801–810 at +/-1e-12 Bohr as the separately
labelled unperturbed-equivalent protocol, and +/-1e-7 Bohr only as a stress
series. Never mix these series into one median. The catalog additions alone
do not establish performance generality beyond silicon.

### A/A floor, then cross-code pairs; one condition per job

1. Perform a single-process preflight first; it is not an SCF memory test.
   Then start with Si128 np8, followed by np16, and use those observed wall/RSS
   values to schedule Si216. A final-comparison qualification job may contain
   one cold run with memory monitoring. Report cold and warmed/cache-hit states
   separately; neither engine may reuse an SCF restart from a previous arm.
2. Before cross-code timing, collect three crossed A/A pairs on the same
   BAANDI binary/input and three on the selected QE binary/input for each
   case/np. Use independent fresh output directories, equal cache conditions,
   and the same binding. Return every external wall sample. For paired
   differences D, report mean, sample SD, SE=SD/sqrt(n), and floor
   `max(abs(mean), 3*SE)`. If `abs(mean)/SE >= 3`, reject that window as drift,
   investigate exclusivity/binding, and repeat without subtracting drift.
   Zero SE is not proof of resolution; take more samples.
3. Only after A/A qualifies, collect at least three crossed BAANDI/QE pairs
   (alternate BAANDI-first/QE-first). Retain raw walls, paired ratios and each
   code's median/min/max, convergence and verification records. If spread is
   above 5% or the gap is near the A/A floor, obtain a fresh predeclared larger
   sample (five or more pairs); do not cherry-pick a favorable window. Report
   BAANDI/QE median-wall ratio for each case/np, with verification included,
   against the target <=0.90. This is a scaling measurement, not an optimization
   adoption claim or a comparison with historical-window walls.
4. One job contains one case, one np and one phase (qualification, BAANDI A/A,
   QE A/A, or paired final comparison); no multi-case loop. Final comparison
   jobs and their required qualification/floor jobs use the `final-` prefix
   and a **7200 s total job cap**, e.g.
   `final-si128-gamma-np8-aa-baandi-01` or
   `final-si216-gamma-np16-compare-01`. The fixture's 7200 s per-run timeout does
   not grant 7200 s to every run inside a job. Stop launching samples when the
   remaining budget cannot fit the observed worst run plus teardown. Split
   the predeclared pairs into numbered jobs/windows if needed, recheck A/A
   when the window changes, and do not pool drifting windows. Non-final stage
   jobs remain subject to the existing 600 s cap. Preserve timeout/OOM/error
   logs and partial evidence; never call an incomplete sample converged.
5. Record maximum rank RSS and aggregate job peak (cgroup/job accounting if
   available, otherwise summed rank RSS sampling with sampling interval and
   limitations). `/usr/bin/time -v mpirun` alone may describe only the launcher;
   it is not a whole-job memory measurement. Account for initialization,
   coarse/fine handoff, verification and final-output materialization. Stop on
   swap/OOM or exhausted headroom. Si512 remains gated as described above.

### Return package

Return through the master: all commands and exit codes, job IDs/time limits,
source/binary/input/UPF hashes and versions, generated BAANDI and QE inputs
(without UPF bytes), each run's stdout/stderr, unedited `results.json` and
QE output, exact grid/PW/electron/band counts, affinity/load/clock records,
per-rank and aggregate RSS/swap evidence, cold/warm state, all A/A and paired
samples/statistics. Return per-run total energy in Ha (convert QE Ry at the
boundary), SCF iterations, Davidson/H applications, verification attempts/status
and wall, final residuals/charge/gap diagnostics, density hashes, and the
available initialization/FFT/nonlocal/algebra/collective/verification phase
breakdown. Use separately instrumented runs for profiling; no profile timings
in the headline median. Preserve raw files remotely with a SHA-256 manifest.

Check deterministic repeats at fixed np (energy/history/density), compare np8
and np16 energies against a separately converged tight reference using standard
5e-9 Ha and strict 1e-10 Ha accuracy classes where applicable, and report the
1e-12 Ha unique-solution cross-path regression check separately. Cross-code
absolute energy differences and internal true-error tests are separate facts;
do not impose a new 1e-12 Ha QE/BAANDI adoption tolerance. No QE reference exists
for these new fixtures until actual commands, version, inputs, UPF hash and
outputs are returned. Failure to converge is a measured outcome, not permission
to loosen thresholds or change bands silently.

## Reproducer and local evidence

All scratch is under `target/fixture-si128/` (ignored). The following Python
3.9+ standard-library program regenerates only the two owned fixture files,
verifies the independent Si64 coordinate set and periodic neighbor counts,
and calls CLI `validate`/`grid` for all sizes. Si512 is generated only as a
scratch CLI input for preflight. Save this block to
`target/fixture-si128/generate_validate.py` and run from the repository root:

```sh
mkdir -p target/fixture-si128/tmp
TMPDIR_990="$PWD/target/fixture-si128/tmp" CHECK_SCCACHE=0 \
  CARGO_NET_OFFLINE=true scripts/check.sh full
python3 -I target/fixture-si128/generate_validate.py
python3 -I -m unittest discover -s scripts/tests -p test_bench_suite.py
```

Adjust only the external read-only UPF path if necessary and record its hash.
The table uses the hash above. The `--generate-only` option skips CLI calls.

```python
import importlib.util, itertools, json, math, pathlib, subprocess, hashlib, sys
root = pathlib.Path.cwd()
out = root / 'target/fixture-si128'
spec = importlib.util.spec_from_file_location('suite', root / 'scripts/bench_suite.py')
suite = importlib.util.module_from_spec(spec)
spec.loader.exec_module(suite)
base = suite.load_case(root / 'benchmarks/cases/si8-gamma.toml')
ref = suite.load_case(root / 'benchmarks/cases/si64-gamma.toml')
a = ref['lattice'][0][0] / 2

def make(n, dims):
    c = dict(ref)
    c.update(id=f'si{n}-gamma', purpose=f'{n}-atom Si Gamma scaling beyond Si64 (G-FIXTURE-SI128); dev08 timing unmeasured',
             lattice=[[a*dims[i] if i == j else 0 for j in range(3)] for i in range(3)],
             atoms=[[(t[i]+b[i])/dims[i] for i in range(3)]
                    for t in itertools.product(*(range(d) for d in dims)) for b in base['atoms']], timeout_s=7200)
    assert len(c['atoms']) == n == len(set(map(tuple, c['atoms'])))
    assert all(0 <= x < 1 for b in c['atoms'] for x in b)
    bond = a*math.sqrt(3)/4
    nearest = float('inf')
    neighbors = [0]*n
    for i, p in enumerate(c['atoms']):
        for j in range(i):
            q = c['atoms'][j]
            distance = math.sqrt(sum((min(abs(p[k]-q[k]), 1-abs(p[k]-q[k]))*a*dims[k])**2 for k in range(3)))
            nearest = min(nearest, distance)
            assert distance >= bond-1e-12
            if abs(distance-bond) < 1e-12:
                neighbors[i] += 1
                neighbors[j] += 1
    assert set(neighbors) == {4}
    return c, {'atoms': n, 'replication': dims, 'minimum_distance_bohr': nearest,
               'minimum_distance_angstrom': nearest*0.529177210903, 'neighbors_per_atom': 4}

check, _ = make(64, (2,2,2))
assert check['lattice'] == ref['lattice']
assert set(map(tuple, check['atoms'])) == set(map(tuple, ref['atoms']))
records = []
upf = (root.parent / 'pseudo/Si.pz-vbc.UPF').resolve()
for n, dims in [(64,(2,2,2)), (128,(2,2,4)), (216,(3,3,3)), (512,(4,4,4))]:
    c, geom = make(n, dims)
    if n in (128,216):
        path = root / f'benchmarks/cases/si{n}-gamma.toml'
        path.write_text('# Lengths in Bohr; fractional atoms replicate the existing Si8 conventional cell.\n'
                        f'# {dims[0]}x{dims[1]}x{dims[2]} conventional cells; Mac preflight only, dev08 SCF.\n'
                        + '\n'.join(f'{k} = {json.dumps(v)}' for k,v in c.items()) + '\n')
        assert suite.load_case(path) == c
    if n == 64:
        c = ref
    path = out / f'si{n}-input.toml'
    path.write_text(suite.case_input(c, upf, 7, 1))
    if '--generate-only' in sys.argv:
        print(json.dumps(geom))
        continue
    result = subprocess.run([str(root/'target/debug/dft-cli'), 'validate', str(path), '--json'], capture_output=True, text=True)
    (out/f'si{n}-preflight.json').write_text(result.stdout)
    (out/f'si{n}-preflight.stderr').write_text(result.stderr)
    assert result.returncode == 0, (n, result.returncode, result.stdout, result.stderr)
    report = json.loads(result.stdout)
    assert report['valid'] and len(report['structure']['atoms']) == n
    assert report['electrons'] == 4*n and report['bands'] == 2*n+4
    grid = subprocess.run([str(root/'target/debug/dft-cli'), 'grid', str(path)], capture_output=True, text=True)
    (out/f'si{n}-grid.txt').write_text(grid.stdout)
    assert grid.returncode == 0, grid.stderr
    record = {**geom, 'lattice': c['lattice'], **{k:report[k] for k in ('electrons','bands','plane_waves_per_k','fft_grid','estimated_array_bytes')},
              'validate_exit_code': result.returncode, 'grid_exit_code': grid.returncode,
              'input_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    records.append(record)
    print(json.dumps(record), flush=True)
if records:
    (out/'sizes.json').write_text(json.dumps({'upf_sha256': hashlib.sha256(upf.read_bytes()).hexdigest(), 'cases':records}, indent=2)+'\n')
```


### Observed check results

- `python3 -I target/fixture-si128/generate_validate.py`: exit 0;
  four structures verified and eight CLI commands (`validate`, `grid`) exit 0.
  JSON/grid output retained under `target/fixture-si128/`; the size tables above
  preserve the results independently of ignored scratch.
- `python3 -I -m unittest discover -s scripts/tests -p test_bench_suite.py`:
  exit 0, **17 tests passed**. The catalog test glob-loads every case and checks
  IDs/rendering; its count assertion is `>=24`, so no test/count edit is needed.
- `TMPDIR_990="$PWD/target/fixture-si128/tmp" CHECK_SCCACHE=0
  CARGO_NET_OFFLINE=true scripts/check.sh full`: exit **1**, build/fmt/clippy
  PASS; **795 passed, 3 failed, 181 ignored**, 87 test binaries, 188 s total
  including the cold build. All three failures are the ticket's known
  pseudopotential-cache policy rejecting an in-repository TMPDIR (the user
  requires scratch inside this worktree). The separate `ps` sandbox diagnostic
  did not prevent the gates from running. This is not a full green gate.
  Failing tests: `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`,
  `pseudo::cache_resolution_and_no_network_side_effects`, and
  `pseudo::cached_fetch_checks_hash_and_provenance`. A direct filtered rerun of
  the first test also exited 101 with the same outside-repositories error
  (`target/fixture-si128/cache-policy-recheck.log`).
- No SCF, MPI execution or QE execution for these new cases on the Mac.
  No energy, convergence or speed result is claimed. No tests or engine code
  changed. No UPF files were copied or added.
- This implementer is restricted to the five ticket-owned paths; registry,
  dashboard uploads and generated root README are left to the master. No new
  SCF validation evidence exists for the README results generator. The local
  sessions registry has no G-FIXTURE-SI128 entry and no local delegate report;
  no new delegated run or recurring orchestrator timer was started here.
