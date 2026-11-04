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
