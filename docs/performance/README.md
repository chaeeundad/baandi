# Performance: rules, measurements and the optimization log

## Reconciled status — 2026-11-05

Source checkout: `b5d515a` on `task/P12A-STATUS-RECONCILIATION`. Presence is
checked against merged source and git history, not stale delegate registry notes.
Historical measurements below remain tied to their original builds.

| Work | Current state | Evidence/index and limits |
| --- | --- | --- |
| P10A | Merged opt-in hybrid RMM-DIIS; Davidson default retained | `rmm-diis-benchmark.md`; buffer policies use Davidson fallback |
| P10B | Optional MKL FFT/GEMM merged; remote acceptance recorded | `fft-backend-comparison.md`; final dev08 relay, not a native Mac MKL run |
| P10C | FFT call-structure changes merged | Optimization log P10C final integration; rejected experiments retained in history |
| P10D | Kerker-Pulay and mixer diagnostics merged, opt-in | `scf-mixing-comparison.md`; no general recovery default adopted |
| P10E | GPU architecture/roadmap delivered | `../architecture/gpu-backend.md`; no GPU implementation |
| P10F | No ticket in this checkout | No feature or completion inferred from the alphabetic gap |
| P10G | General benchmark harness and suite delivered | `benchmark-suite.md`, `suite-baseline.json`; synthetic and real-Si coverage distinguished |
| P10H | Superseded by P10J | Git `0cdcf41`; not a separate unmerged backend requirement |
| P10I | Intra-k linear algebra merged; scheduling superseded by J | Optimization log P10I; historical scoped-worker/MKL panel path |
| P10J | Persistent executor and backend-specific products merged | Optimization log P10J; no multiplied outer/inner thread budget |
| P10K | Synchronization/validation glue merged | Optimization log P10K; local timings are not remote scaling acceptance |
| P10L | Fixed pieces/batching and utilization instrumentation merged | Optimization log P10L; utilization/scaling targets not all achieved |
| P10M | Copy and duplicate-GEMM removal merged (`1cae61c`) | Optimization log P10M; CholeskyQR2 rejected after synthetic iteration regression; not an enabled default |

Current cross-code speed claims use only
`tests/reference/perf/dev08_final_relayed.md`, **Part 4**, fdce016, 2026-11-05.
On dev08, shifted 4x4x4 Si at 8–16 cores is 2.0–2.2 times QE at conv_thr
1e-8 and 1.4–1.6 at 1e-12; the ratio changes with QE iteration count.
The 8-atom Gamma pure-Rust time is 5.801 s (T=1) versus 4.151 s (T=8),
only 1.40x scaling. MKL is not uniformly faster: Part 4 reports si8 T=1
6.016 s versus pure 5.801 s. Both are measured, relayed results, not rerun here.
QE internal times below 0.3 s are excluded from README ratios because its
resolution is 0.01 s; the si8 16-core 0.3000 s point is at the inclusion boundary.
Part 1's fe0492b table and R7's pre-P10 tables are historical, not current values.

P11J Part 3 reports all 18 combinations meeting overhead targets with 300 pairs;
that is provisional relayed evidence. Local `tests/reference/runtime/benchmark.json`
retains its own failed/inconclusive targets. The remote raw JSON and independent
reproduction are still TODOs. See [evidence map](../validation/evidence-map.md).

Read this before touching a hot path (anything on the Hamiltonian apply, the eigensolver, FFTs, density accumulation). Physics correctness always comes first; a change that is faster but not equivalent is a bug.

## Rules for every agent
1. Measure first. Profile with `sample` (macOS) or `perf` (Linux) on a release build with debug symbols (`CARGO_PROFILE_RELEASE_DEBUG=true`), read the time at the top of the stack, and write the finding into [optimization-log.md](optimization-log.md) before changing code. No optimization on a guess.
2. Prove equivalence. Use `scripts/bench_scf.py` against [baseline.json](baseline.json) (Si Gamma and k-point cases): the total energy must stay within 1e-12 Ha, and a change that only removes overhead (allocation, copies, skipped zero work) must keep the density file bit-identical (`--require-identical`). A change that alters the algorithm (tolerances, subspace, number of iterations) must say so, keep every convergence gate and the fixed-point certificate, and give the measured energy difference.
3. No allocation inside the per-band or per-iteration loops: allocate work arrays once (the plan, the Hamiltonian, the solver workspace) and reuse them. No logging, I/O or timers in numerical hot loops (macOS `free` and `malloc` even call the timer internally, which shows up as `mach_absolute_time` in profiles).
4. Skip work that is exactly zero (zero lanes in a pruned FFT, empty projector blocks): this is bit-identical and the safest speedup. Prefer general-purpose libraries (faer, rustfft, ndarray) used with their allocation-free entry points (`process_with_scratch`, preallocated buffers).
5. Keep the physics-facing API stable; optimize behind it. A speedup that needs an API change gets an ADR.
6. Record every accepted optimization in [optimization-log.md](optimization-log.md): what the profile showed, what changed, before and after numbers from `scripts/bench_scf.py` on the same machine, equivalence evidence, and what was tried and rejected. Update `baseline.json` only through the decision recorded in the log (never to hide a regression).
7. Report speedups only from measured runs on a quiet machine, with the machine and commit stated; the cross-machine comparison with QE is in `docs/validation.md` R7 and README.

## Tools
- `scripts/bench_scf.py` (standard library only): runs the release `dft-cli` on fixed Si cases (`gamma`, `k222s`, `k222`, `k444s`), reports median wall time, engine and diagonalization time, iterations, energy and the density hash, and compares with the baseline. The pseudopotential path is passed with `--upf` (the file is not in the repository).
- Baselines: [baseline.json](baseline.json) is the state before any optimization (macOS Apple Silicon, release, 1 thread); [baseline-after-o4.json](baseline-after-o4.json) is the state after the Davidson product-reuse change (O4b changed the arithmetic of the residual, so density files differ bit-wise from the first baseline while the energy agrees to 8e-15 Ha): use it for `--require-identical` checks of later overhead-only steps. Timings on the shared Mac are indicative; accepted timings come from the quiet dev08 node (docs/validation.md R7 describes it).

## P10I deterministic intra-k linear algebra

Historical implementation; P10J below supersedes its worker scheduling and MKL path.

The eigensolver obtains the existing intra-k budget from
`LinearOperator::inner_threads()` (default 1); Hamiltonian reports its band-worker
budget and the hybrid counting adapter forwards it. The outer/inner SCF split
is unchanged. GEMM and Hamiltonian work are sequential regions, so their worker
budgets do not multiply. EVD remains sequential.

`GemmWorkspace` uses fixed 256-element panels, including at one thread:
rotations split result rows, while Gram products split the plane-wave reduction
axis and add partial matrices in ascending panel order. The latter is necessary
because subspace Gram results have only a few rows. No reduction grouping depends
on worker count. Panels with fewer than 65,536 complex multiply-adds per
panel execute sequentially to avoid scoped-thread overhead; this scheduling
guard does not change their boundaries or arithmetic. Pure Rust uses borrowed faer submatrix views; MKL uses reusable
packed panels and the existing sequential `gemm_into` adapter. Panel storage is
allocated when the solver workspace is created and reused. Scoped OS-thread
bookkeeping is transient, matching the existing band-worker convention.

`ProjectionWorkspace` replaces projection against the established Davidson basis
with block GEMMs, interleaving projection and internal two-pass MGS twice.
Normalizing the first projected block before reprojecting is essential for
nearly dependent incoming columns. The small incoming block and RMM's internal
orthogonalization retain MGS; thresholds, band policy, certificates and defaults
are unchanged. Projection/GEMM rounding differs from historical builds, so the
acceptance contract is energy agreement within 1e-12 Ha and bit-identical density
across thread budgets within the same backend, not historical density identity.

The allocation-free boundary is `GemmWorkspace::multiply_into`; projection
scratch likewise does not allocate per panel/column. Legacy Davidson still owns
new result matrices and vectors each iteration, as before this change; this is
not a claim of a wholly allocation-free Davidson solver. RMM reuses its existing
outputs. See the P10I entry in [optimization-log.md](optimization-log.md) for
measurements, rejected prototypes, test commands and remote acceptance limits.

## P10J executor and backend-specific decomposition

`ScfProblem` owns an explicit Rayon executor, shared through Hamiltonian and
`LinearOperator::worker_pool()` with both eigensolvers. T=1 runs inline without
creating an OS worker. Nested k, Hamiltonian and GEMM tasks share the total
budget; the global Rayon executor is never configured. Pool construction
counts worker startups once; submission does not update the counter.

Pure Rust retains the P10I fixed 256-plane-wave Gram partials and ordered sum.
Rotation panels write directly through output views. Calls below 262,144 complex
multiply-adds (2,097,152 real flops), or with fewer than four panels, run inline;
the criteria do not depend on requested threads. MKL uses one sequential BLAS
call per logical product at every thread budget, avoiding panel copies and
assembly. Thus bit identity is required across budgets within a backend, while
changed MKL historical rounding still needs the remote energy gate.

Density uses the same SCF pool across k points and bands. Each band's
contribution is independent of scheduling; contributions are reduced in ascending
band order, then k densities in ascending k order. Its output array is reserved
at construction and reused. Hamiltonian tasks split the column-major output
slice directly, removing the old per-band output mutex and copy.

Each complete k-solve or energy-audit phase is installed once in the executor,
including single-k calculations and the entire bounded k-chunk loop. Nested GEMM and band jobs then remain resident in that
pool rather than synchronizing with the external calling thread for each call.
Parallel GEMMs schedule individual fixed tiles instead of tying a multi-tile
chunk to each worker; work stealing can balance uneven completion times without
changing any tile or reduction boundary. No new unsafe code or boxed task is
introduced. Legacy allocating matrix/result ownership APIs remain outside the
allocation-free `multiply_into`, Hamiltonian batch and density transform paths.

Release-only scheduling and traversal experiments:

```
cargo test --release -p dft-numerics dispatch_latency -- --ignored --nocapture
cargo test --release -p dft-numerics granularity_table -- --ignored --nocapture
cargo test --release -p dft-physics nonlocal_traversal -- --ignored --nocapture
cargo test -p dft-engine --lib scf_reuses_workers -- --nocapture
```

The last test runs complete synthetic SCFs and checks executor startup counts
before and after all iterations. See [optimization-log.md](optimization-log.md)
for measured timings, equivalence, outstanding dev08 gates and limitations.

## P10K synchronization and validation glue

The single-k SCF branch now selects serial execution from the **total** thread
budget (`threads <= 1`), not the outer k-point count. With an intra-k budget,
Gamma uses the same resident ordered-solve path as multiple k points. P10J's
claim above that every complete solve was resident did not hold for its separate
single-k branch; that omission is corrected here.

Davidson recurrence rotates the basis and Hamiltonian image through one shared
fixed-panel task tree. Each panel retains the original separate sequential
products; the Gram reduction and all panel boundaries remain unchanged. MKL and
T=1 retain separate rotation calls. Projection reuses the entry validation of
its immutable basis and the output validation of coefficients. Public GEMM still
checks both operands, dimensions and output; newly normalized projection inputs
are checked before use. No unsafe entry point, output-only public validation,
new tolerance, spin policy or thread-dependent arithmetic was introduced.

Panels remain individually stealable. A prototype with only T groups regressed
si8 on the Mac and was removed; a task count cap is not a promise that unused
Rayon workers never wake. Hamiltonian already applies all supplied bands in one
region, splitting them evenly (at most one band difference); setup already has
one region covering k-state construction. Their kernels remain unchanged.

See the P10K entry in [optimization-log.md](optimization-log.md) for the measured
Mac matrix, instrumented ingress/region counts, rejected candidates and remote
acceptance limits. Instrumentation was confined to temporary benchmark builds;
production hot loops contain no new counters, clocks or logging.

## P10L worker analysis and dependency overlap

The production Gram and rotation decomposition remains the P10K 256-row/pw
partition. Finer 64-row rotations and four/eight-column Gram subtiles were
measured and rejected: they changed backend rounding and, for si8, the amount
of Davidson work. Their apparent speedup is not accepted as a scheduling win.

For projection blocks with at least 2048 rows, the first MGS pass updates
independent future columns as each accepted basis vector becomes available.
Each column still receives the same ascending sequence of two-pass scalar
projections. When future work contains at least 16384 complex elements, the
next column's dependent second pass overlaps first-pass updates on later
columns. Small blocks retain the original left-looking loop. Rank decisions,
normalization, tolerances, correction widths and certificates are unchanged.
These are scheduling/copy changes, not a new orthogonalization method.

Hamiltonian application uses one stealable piece per band when the inner
budget owns the whole pool and there are more bands than scratch slots.
Per-worker FFT/nonlocal scratch is preallocated; a slot mutex protects reuse,
with no nested scheduling while held. Outer-k and small-band cases retain
static disjoint groups: mapping a larger pool onto a smaller scratch budget
with modulo alone caused contention in the experiment. T=1 accesses scratch
exclusively without locking. Concatenation reserves its final capacity before
copying, avoiding a reallocating second copy.

### Reproducible analysis build

The ticket owns source directories but not Cargo manifests. This source-only
analysis feature is enabled via rustc cfg, **not** `cargo --features`. Default
builds compile only inline zero-sized no-ops: no clocks, TLS, atomics, trace
allocation or output. Analysis builds allocate/buffer records and perturb
scheduling; use a separate uninstrumented binary for acceptance timings.

```sh
RUSTFLAGS='--cfg feature="worker-analysis"' cargo build --release -p dft-cli
cp target/release/dft-cli /tmp/p10l-analysis
python3 -I crates/dft-numerics/src/pool/measure.py \
  /tmp/p10l-analysis /tmp/p10l-analysis-results /absolute/path/Si.pz-vbc.UPF
RUSTFLAGS='--cfg feature="worker-analysis"' cargo test -p dft-numerics \
  context_follows_stolen_work_and_restores_after_panic
cargo test --release -p dft-engine worker_utilization_synthetic_twenty \
  -- --ignored --nocapture --test-threads=1
```

Run the Python driver from the repository root and select a fresh output
folder. It measures si8 and a skew four-atom cell at T=1,2,4,8,12 and retains
inputs, raw traces, results and density files. `summary.json` reports per-phase
wall/parallel/between-region times and measured worker busy fractions.
`regions.csv` preserves every region's piece count, scheduled slots, individual
worker busy fractions and max/mean ratio (mean includes idle pool workers).
No polling-based busy-worker distribution is inferred from these counters.

Raw `P10L` CSV fields are label, start ns, end ns, scheduled slots, completed
pieces, 256 semicolon-separated worker busy times in ns (empty for plain
inclusive spans), and a nested-busy flag (0/1). The flag marks a serial phase
already charged to its enclosing region; the summary counts it only once. Slot 255 is the external caller, slots 0--254 are Rayon
indices. T=1 therefore reports external-caller activity. `workers=0` denotes a
phase span, not a dispatch region. Regions with one scheduled slot are inline;
`parallel_regions_s` excludes these. Logical regions include both parallel and
inline kernels. Piece timing excludes scratch-lock acquisition and scheduler
waits. Region busy time also includes explicitly timed serial work inside the
region (such as Gram assembly and the overlapping dependent MGS pass). Explicit serial phases cover EVD, residual norms, preconditioning, copies,
projection preparation/update, Gram assembly/validation and dependent MGS. Unclassified glue remains
unmeasured busy time, so busy fractions are lower bounds, not CPU utilization.
Phase spans are inclusive and must not be summed. Region intersections use
interval unions to avoid double-counting nested scopes.

Only one active SCF per process and at most 255 workers are supported by the
analysis feature. Classic `run`/`run_outcome` flushes at the end, including error
returns; direct incremental sessions need their own outer analysis scope.
The SCF scope starts after `ScfProblem` construction and excludes initial setup
and CLI I/O. Whole-pool residency includes serial Davidson code and must not be
substituted for the useful parallel fraction in Amdahl's law.

The synthetic twenty-atom test exercises a skew/slab-like arrangement and both
local-only and multi-projector potentials, comparing full density and energy
bits across thread counts. It deliberately retains six iterations, and is
**not** a converged material calculation or a QE reference. See the optimization
log for retained/rejected experiments, actual gates, and unresolved dev08 goals.

## P10M: retained storage changes and rejected block-QR default

Davidson reserves one live-prefix subspace matrix, two alternating H-subspace
buffers and one correction-product buffer. Appending copies only new columns;
reading the cached H basis borrows its live prefix with the original column
stride. Truncating the owned Ritz-coefficient matrix changes its live width
without copying. A potentially overflowing basis preserves the already computed H rotation
before residual formation; a valid contracted restart reuses that buffer.
The fallback restart GEMM writes directly to the alternate reserved H buffer.
Unused capacity is reserved without initializing columns that will not be read. The accepted default keeps P10L orthogonalization and all
arithmetic/reduction boundaries. An external operator that replaces its output
allocation can force capacity restoration; the engine Hamiltonian preserves it.
This is not a claim that all legacy Davidson allocations have been removed.

The CholeskyQR2 prototype is retained **disabled by default**. To reproduce the
experiment, use `RUSTFLAGS='--cfg feature="block-qr"' cargo build --release -p
dft-cli` (a source-only experiment, not a Cargo feature or SCF default). It uses
faer's sequential Cholesky/triangular solve for the small factor and two Gram /
rotation passes, with persistent scratch. Projection-only 64-row/pw panels give
47 pieces at 2969 plane waves; partials are assembled in ascending panel order.
Small matrices below the existing 2048-row boundary retain MGS. If normal-equation
conditioning is unsuitable, the input is restored and existing MGS decides rank
with the unchanged dependence tolerance. This experiment failed the synthetic
20-atom skew-cell Davidson-iteration gate, despite passing the five Si energy
checks. It must not be enabled as a general default on the Si result alone.

`worker-analysis` now also emits `P10M_COST,start_ns,end_ns,logical_2mnk,
executed_2mnk,explicit_copy_bytes` once per Davidson iteration, at the outer SCF
flush. The iteration context follows stolen/nested work and is distinct from
phase timing, so concurrent k solves do not share a global difference counter.
The FLOP convention here is the requested `2*m*n*k` (complex multiply/add
counting convention; conventional real FLOPs are four times this number).
Executed counts use actual backend-call dimensions, including partial tails;
logical counts use complete requested products. These counters exclude Gram
assembly additions (`panels*m*n` complex adds per Gram), finite checks, scalar
MGS, factorization, and opaque backend microkernel padding/packing. Explicit
copy bytes count source elements copied once, not combined read+write traffic;
they include EVD adapter copies, projection copies, cache/append/restart and
coefficient copies. Scope is the owned numerical/solver copy sites; unmodified
Hamiltonian/FFT copy sites are not instrumented. GEMM counts cover
`GemmWorkspace`, not every GEMM elsewhere in the engine. They exclude memset/zero initialization, implicit allocator
traffic, compiler-generated copies and library-internal memcpy. Thus a ratio of
one rules out extra GEMM dimensions at this adapter boundary, not hidden MKL
implementation work or shape-dependent kernel efficiency.

Reproducible bounded synthetic and paired timing drivers (run from repo root):

```sh
python3 -I crates/dft-numerics/src/pool/orthonormalize_synthetic.py \
  /tmp/before /tmp/new-synthetic-before
python3 -I crates/dft-numerics/src/pool/orthonormalize_synthetic.py \
  /tmp/after /tmp/new-synthetic-after /tmp/new-synthetic-before/summary.json
python3 -I crates/dft-numerics/src/pool/orthonormalize_bench.py \
  /tmp/before /tmp/after /absolute/path/Si.pz-vbc.UPF /tmp/times.json \
  --cases si8 --repeats 30
cargo test --release -p dft-numerics block_qr_si8_shapes -- --ignored --nocapture
RUSTFLAGS='--cfg feature="worker-analysis"' cargo test -p dft-numerics \
  iteration_costs_follow_workers -- --nocapture
```

The synthetic driver exercises 4/8/20 atoms, cubic/skew cells, one/three
projectors, and T=1,2,4,8. Its six-step outcomes deliberately remain unconverged;
they test arithmetic/thread identity and iteration cost, not physical accuracy
or convergence. Timings use alternating pairs after discarded warmups; thirty
samples per binary/budget are required for reported performance ratios. See the
P10M optimization-log entry for the failed candidate, exact counts, gate commands
and remaining remote acceptance work.
