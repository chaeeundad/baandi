# Performance: rules, measurements and the optimization log

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
