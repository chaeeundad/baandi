# ADR 0003: Prefer general-purpose libraries over hand-written numerics

Status: Accepted (2026-11-02), by user directive. Supersedes the "own 3D FFT loop" and "hand-written Accelerate FFI" parts of ADR 0001 and ADR 0002 where a mature library covers the need.

## Policy
1. Do not hand-write numerical kernels or FFI bindings when a maintained, permissively licensed crate covers them. Hand-write only physics-specific glue (basis mapping, potentials, Hamiltonian, SCF logic).
2. Every adopted crate is checked for license (MIT/Apache-2.0/BSD-class preferred), last release, download base and Apple Silicon support, then recorded in `docs/dependencies.md` with the pinned version.
3. Libraries sit behind the narrow traits from ADR 0001 so backends remain swappable (CPU now, CUDA later). Adapters stay thin.
4. Cross-checks against an independent library are encouraged in tests (e.g. dense Hermitian eigensolve vs the iterative solver).
5. Claude fetches crates (the Codex sandbox has no network) and reviews licenses before they enter a manifest.

## Candidates (crates.io, checked 2026-11-02)
| Need | Candidate | Version / last release | License | Status |
| --- | --- | --- | --- | --- |
| N-d arrays | ndarray | 0.17.2 / 2026-01 | MIT OR Apache-2.0 | adopt (to be confirmed by spike) |
| Complex numbers | num-complex | current | MIT OR Apache-2.0 | adopt |
| Dense linear algebra, pure Rust | faer | 0.24.4 / 2026-06 | MIT | spike candidate (complex Hermitian eig, gemm) |
| Dense LA via system LAPACK | lapack + lapack-src / blas-src (Accelerate feature), ndarray-linalg 0.18.1 | 2025-2026 | MIT OR Apache-2.0 | spike candidate; `accelerate-src` is stale (2018), prefer `blas-src`/`lapack-src` features |
| FFT 1D | rustfft | 6.4.1 / 2025-09 | MIT OR Apache-2.0 | adopt (ADR 0002) |
| FFT N-d | ndrustfft | 0.6.2 / 2026-01 | MIT | replaces our own 3D loop if the spike passes |
| XC functionals | libxc (wrapper, C library) or xcx (pure Rust) | libxc 0.1.4 / xcx 0.3.1 | Apache-2.0 / not yet checked | spike candidate; low adoption for both, so use as cross-check even if LDA-PZ is kept in-tree |
| UPF reading | upf | 0.1.1 / 2026-04 | not yet checked | check license and field coverage before use; our own validation layer on top stays |
| Input/serde | serde, toml | current | MIT OR Apache-2.0 | adopt |

## Spike (ticket P0C) must decide
- faer vs LAPACK-via-Accelerate for complex Hermitian eigenproblems and complex gemm: correctness against analytic cases (2x2 [[2, 1-i],[1+i, 3]] has eigenvalues 1 and 4), build simplicity on Apple Silicon, runtime on a 500-1000 dimension matrix.
- ndrustfft 3D transforms: round trip, Parseval and normalization on non-power-of-two grids (15, 20, 24, 27); NEON effect.
- Whether the existing hand-written `dft-numerics` FFI (real ddot/dsyev) is deleted.

## Consequences
Less code to own and verify, but more dependency review. API churn in 0.x crates (faer, ndrustfft) is contained by the adapter traits and by pinning exact versions in `Cargo.lock`.
