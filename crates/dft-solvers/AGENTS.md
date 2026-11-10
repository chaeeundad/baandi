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

## Research switches of the real Davidson (2026-11-10)

`RealDavidson::solve_research` takes `ResearchOptions` explicitly; every other
entry point reads them from the environment (`davidson_stages.rs`): stage timers
`BAANDI_DAVIDSON_STAGES=on`, B1 selected-root products `BAANDI_B1=standard|strict|on`
with `BAANDI_B1_SKETCH=count`, B5 bordered eigen-update `BAANDI_B5=on` (`bordered.rs`),
B6 cached projector products for the return certificate `BAANDI_B6=on` (operators opt in
through the `RealLinearOperator` default methods; test `tests/b6_certificate.rs`).
Micro switches read at construction: `BAANDI_WIDE_DOT=on` (eight-lane MGS dot, rounding
class); see docs/performance/requests/micro-stack.md.
All default off and bit-identical when off. Tests: `tests/b1_selected_roots.rs`,
`tests/b5_reduced_eigen.rs`; B4 harness `tests/b4_inner_iterations.rs` (ignored).
Plan and pre-registered rules: docs/performance/requests/b1-b4-b5-davidson.md.
