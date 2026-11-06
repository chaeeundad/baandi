# dft-solvers

**Purpose:** Eigensolvers (dense baseline, then Davidson/LOBPCG), mixers (linear, then Pulay), convergence checks.

**Invariants to preserve and test:** Residual and orthonormality reported; convergence needs density, energy and eigen residual, never energy alone.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-solvers`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.

## P26D: isolated real Gamma Davidson

`real_davidson.rs` is the real ORC solver family. `RealLinearOperator` is defined
in this crate; physics adapters belong to R5. No dft-physics dependency and no
engine caller are introduced. Keep complex Davidson and its helpers untouched.
Block width is bands; capacity is `min(n, max(max_subspace, 6*bands))`. The kinetic
diagonal is in ORC order and Ha; cold starts sort by kinetic energy with coordinate
order breaking ties. Strict returns require a true residual certificate; the
separate Ritz-change method always returns `residuals_certified = false`.

Use the existing `RealGemmWorkspace` and shared `WorkerPool`, with 256-row tiles
and ordered sums. The real dispatch work estimate is 1,048,576 multiply-adds
(2 flops each), not 262,144 complex multiply-adds (8 flops each). Do not retune
parameters, add scoped threads, use CholeskyQR2, or infer spectral completeness
from a residual certificate. Rank tolerance is the shared 1e-12 dimensionless
constant; restart contraction requires orthogonality defect <=1e-12; the kinetic
denominator floor stays 0.1 Ha.

Fast fixed-operator tests: `cargo test -p dft-solvers --test real_davidson`.
Parallel dispatch and timing observation (ignored in debug):
`cargo test --release -p dft-solvers --test real_davidson -- --include-ignored --nocapture --test-threads=1`.
Evidence and limitations: `docs/performance/real-gamma-r4-notes.md`.
