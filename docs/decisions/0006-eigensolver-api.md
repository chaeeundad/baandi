# ADR 0006: Backend-independent Hermitian eigensolver stack

- Date: 2026-11-02
- Task: P4A-EIGENSOLVER
- Status: proposed for reviewer acceptance; implemented

## Context

P3 needs a dense Hermiticity oracle and later plane-wave work needs a matrix-free
solver. Both must be testable without physics. ADRs 0001, 0003 and 0004 isolate
faer in dft-numerics and prefer its established eigendecomposition and GEMM.
No dependencies or manifests change in this task.

## Matrix and numerical adapter

`dft_numerics::linalg::ComplexMatrix` owns private rows, columns and
`Vec<Complex64>` storage, indexed at `i + j*rows`. `new` validates length and
representable storage; `zeros` validates multiplication and byte size before
allocation. Accessors expose dimensions, a read-only storage slice, and mutable
or read-only individual columns. Indexing out of range panics, like Rust slices;
allocation exhaustion is not recoverable. Arithmetic entry points reject
non-finite data. Zero-column blocks are supported.

`hermitian_eigh(&ComplexMatrix)` returns ascending `Vec<f64>` values and matching
eigenvector columns. It rejects nonsquare inputs and checks **both** triangles
including the imaginary diagonal: maximum entrywise `|A_ij-conj(A_ji)|` must be
at most `1e-12 * max(1, max_ij |A_ij|)`. Within this allowance faer uses the lower
triangle and real diagonal; the routine does not silently accept a larger defect.
An empty square matrix returns empty outputs. `LinalgError` distinguishes shape,
non-finite arithmetic, non-Hermiticity (measured defect and allowance), and backend
decomposition failure. No faer types cross the public boundary.

`gemm(a, op_a, b, op_b)` returns `op(a)*op(b)` with each `MatrixOp` either
`Identity` or `Adjoint`. It uses faer GEMM, including rectangular products.
Both GEMM and the low-level faer self-adjoint EVD explicitly use `Par::Seq`;
neither changes process-global parallelism. Scratch and conversions are private.

`orthonormalize(block, against)` returns an `Orthonormalized` result containing
accepted vectors and the original indices of dropped columns. `against` is
optional and must already have orthonormal columns; this is a caller contract.
Two-pass modified Gram-Schmidt projects against both that basis and previously
accepted block columns on each pass. It drops zero columns and columns whose
remaining norm is at most `1e-12` of their original norm. Columns are normalized
before projection to make the relative rank test work even for tiny inputs.
Norms use repeated hypot to avoid unnecessary overflow/underflow. The two-pass
method is chosen over QR because incremental basis extension and reporting
original dropped indices are explicit requirements; it needs no pivot mapping
or backend-specific rank policy. This is the ticket's permitted small rank-filter
algorithm, not a replacement for backend eigensolve or GEMM.

## Operator and solver API

`LinearOperator` exposes `dimension`, `apply(x, y) -> Result<(), SolverError>`
and an overridable `apply_block(x, y)` whose default loops over columns. The
operator is Hermitian by contract and overwrites every output component. Slice
lengths supplied by solvers equal `dimension`. Block implementations must preserve
output dimensions. All solver applications check the returned block for finite
data, including custom overrides. Operator errors propagate.

`Preconditioner::apply(eigenvalues, residuals, output)` transforms a residual
block using its matching Ritz values. It is independent of the operator.
`DiagonalPreconditioner::new(diagonal, denominator_floor)` requires finite
entries and a positive finite floor in eigenvalue units. It applies
`-r_i / safe(d_i-lambda)`, with the denominator magnitude clamped to the floor,
preserving the difference's sign; zero takes positive sign. The floor is caller
selected (0.1 in synthetic tests), not an assumed physical energy scale.

`DenseEigensolver.solve(operator, bands)` applies the operator to unit columns,
then calls `hermitian_eigh` and retains the lowest bands. Its fixed limit is
1024: larger dimensions return `DenseLimit` before allocation/application.
It is an explicit correctness oracle, not an automatically selected fallback.
It returns measured residuals, one iteration and one history entry. Unlike
Davidson, it has no residual tolerance parameter; the caller assesses its
reported dense residuals.

`Davidson` parameters are `bands`, positive finite absolute `tolerance`,
`max_iterations`, `max_subspace`, positive `block_size`, and `seed`.
`solve(operator, initial, preconditioner)` takes optional borrowed initial
columns and an optional borrowed preconditioner. Require `1 <= bands <= n`,
a positive iteration limit, and `min(max_subspace,n) > bands` unless `bands=n`.
Initial guesses must have n rows and no more than the effective subspace capacity.
Dependent guess columns are dropped. Missing rank up to
`min(max(bands,block_size),capacity)` is filled with local SplitMix64 complex
columns; this uses no random crate or global state.

Each iteration applies H to the orthonormal basis, diagonalizes `V^H H V`, forms
the lowest Ritz vectors, then **reapplies H** to compute every actual residual
`H x_j - lambda_j x_j`. Convergence requires every residual 2-norm at most the
absolute tolerance, in the same units as eigenvalues. No eigenvalue-change
criterion, spectral-pattern detection, or known-answer shortcut is used.
Success returns `EigenSolution` with ascending values, matching column-major
vectors, per-band residuals, iteration count and maximum-residual history.

Up to `block_size` unconverged residuals in ascending band order are transformed
and orthogonalized to extend the basis. The number is also bounded by
`capacity-bands`, so a restarted space can accommodate them. Before an extension
would exceed capacity, the solver restarts from all current Ritz vectors,
reorthonormalized. It does not lock individual vectors or remove duplicate
eigenvalues: an exactly degenerate level retains its vector multiplicity.
If all corrections are dependent, it tries fresh seeded independent directions.
Exhausted rank or iteration budget returns `SolverError::NotConverged` with
iterations, maximum residual and all per-band residuals. Bad parameters,
application failures and numerical adapter errors are separate typed variants.
There is no logging, JSON or file I/O in the solver loop.

## Alternatives and limits

LOBPCG, physics operators and performance tuning are deferred. Full-space
Rayleigh-Ritz is not used as a hidden recovery path. Restarting only the requested
Ritz block keeps the API simple but can require many iterations for clustered
spectra. These are residual certificates, not a proof that no lower eigenvalue
exists outside an invariant initial space: a user-supplied invariant guess can
omit lower states. Seeded initial blocks mitigate this generic iterative-solver
limitation. Dense comparison is essential for validation.

Sequential execution and fixed seeds give repeatability on the tested build;
bitwise equivalence across architectures or backend versions is not promised.
The API is allocation-owning, not allocation-free. Very large finite magnitudes
can still overflow intermediate operations and return `NonFinite`.

## Validation and reproduction

All synthetic values and residuals are dimensionless, not Hartree measurements
or QE references. Tests print measured errors with `--nocapture`.

- Numerics: analytic complex 2x2 and Fourier 4x4 eigenvalues (1e-12), seeded N=60
  residual and orthonormality (1e-10), all four rectangular GEMM transpose
  combinations against naive loops (1e-12), dependent column removal and
  orthogonality against an existing basis (1e-12), invalid and non-finite input.
- Solvers: analytic and diagonal cases, exact threefold and 1e-6 near-degeneracy,
  bands=2/3 multiplicity, restart with capacity 9 for n=18, deficient initial
  guesses, same-seed identical solution/history and different-seed agreement,
  explicit nonconvergence diagnostics, invalid parameters, dense limit,
  non-Hermitian input, application failures, and denominator floor signs.
- Ignored release tests: seeded Hermitian N=200, bands=8 (dense eigenvalue
  difference <1e-9, residual <1e-10, orthonormality <1e-10, operator Hermiticity
  relative <1e-12); N=300 diagonal `0.5*k^2` plus 0.02 times seeded Hermitian
  coupling (residual <1e-9, eigenvalue difference <1e-9), asserting fewer
  iterations with the diagonal preconditioner. These are excluded from the
  fast debug gate because repeated large projected decompositions can exceed
  its approximately ten-second per-test budget.

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
# Use the existing offline toolchain alias if needed on this host:
export RUSTUP_HOME=/private/tmp/p0c-rustup
export CARGO_NET_OFFLINE=true
scripts/check.sh
cargo test --offline -p dft-numerics -p dft-solvers -- --nocapture
cargo test --offline --release -p dft-solvers -- --ignored --nocapture
```

Gate JSON and session registry files are outside this ticket's owned paths;
the final structured report supplies measured evidence for the reviewer to
reconcile. Changes remain uncommitted for Claude review.

## P7P3 implementation notes

Davidson now retains unchanged basis products and applies only appended columns.
Reorthonormalized restarts invalidate the cache and receive a fresh block apply.
Bitwise equivalence of this caching step requires column-independent rounding
in `apply_block`, as in the default implementation and the Hamiltonian; a custom
GEMM implementation may select different arithmetic for different block widths.
The additive `solve_with_recurrence(..., bool)` separates cached true-residual
and cached recurrence-residual modes without changing Davidson's fields or the
existing solve signature. Recurrence convergence is always certified by one true
block application; failed certification records the actual maximum in history
and uses actual residual directions to continue. `EigenSolution` additionally
reports vector-column `operator_applications`, including certification (dense
solutions count unit-column construction and residual applications).
Stage 2 is the default after passing the P7P3 SCF energy (1e-12 Ha), fixed-point
certificate and workspace gates; stage 1 remains selectable with `false`.
