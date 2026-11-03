# ADR 0004: Numerical backend libraries

Status: Accepted (2026-11-02, Claude Code, after review of spike P0C). Applies ADR 0003. Evidence: `docs/backend-spike-results.md`; code preserved at git tag `spike/p0c` (not merged: its debug test takes 707 s).

## Decision
| Need | Library (pinned in Cargo.lock) |
| --- | --- |
| Complex Hermitian eigensolve (dense) and complex GEMM | `faer` 0.24.4 |
| 3D complex FFT | `ndrustfft` 0.6.2 over `rustfft` 6.4.1 (axis APIs and reusable handlers) |
| Arrays, complex numbers | `ndarray` 0.17.2, `num-complex` 0.4.6 |
| Independent test oracle (dev only, optional) | `lapack` 0.20.0 + `blas-src` 0.14.0 / `lapack-src` 0.13.0 with the `accelerate` feature |

The hand-written Accelerate FFI in `dft-numerics` (real ddot/dsyev from P0) is removed in P1; no production `unsafe` for numerics remains. Backends stay behind the narrow traits of ADR 0001.

## Evidence (release build, one thread, median of 3 after warm-up)
- Eigensolve N=1000: faer 0.369 s, LAPACK zheev 0.629 s; max eigenvalue difference 2.1e-13; both meet residual and orthonormality tolerances. Analytic 2x2 [[2,1-i],[1+i,3]] gives 1 and 4 within 1.3e-15 (LAPACK) and 2.2e-16 (faer).
- GEMM N=500: faer 0.0189 s vs naive loop 0.0921 s, relative error 2.5e-15.
- 3D FFT: grids 4x5x6 (vs naive DFT), 15x20x24, 16^3, 27^3, 32^3, 48^3 all pass round trip (worst 1.6e-15), Parseval (worst 1.9e-14) and signed single-mode tests. Measured normalization: forward unscaled, inverse divided by grid size.

## Not established (do not claim)
Accelerate zgemm vs faer; zheevd/zheevr or threaded defaults; multithread scaling; NEON effect; custom rustfft axis loop speed; behavior on nearly degenerate or badly scaled spectra; peak memory. Revisit when profiling real SCF workloads.

## Rules that follow
- Pin exact versions; faer and ndrustfft are 0.x, so wrap them in our traits to contain API churn.
- Tests that take more than about 10 s in debug must be `#[ignore]` with a documented release command; `scripts/check.sh` must stay fast.
- Crates are fetched and vendored by Claude (Codex sandbox has no network); only crates listed here may be used without a new ADR.

## P7P6 notes

Optional standard-library workers parallelize only independent k-point Davidson
solves (including certificate solves). Each faer invocation remains `Par::Seq`
and each FFT remains serial. All k reductions retain their original arithmetic
order for bit identity. Default execution remains one thread; Gamma uses one
thread regardless of the request. Diagonalization timings are region wall time.
