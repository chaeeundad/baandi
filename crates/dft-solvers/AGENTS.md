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
Algebra-stream switches (docs/performance/requests/algebra-stream.md, test
`tests/algebra_stream.rs`): `BAANDI_RITZ_FUSED=on` (Ritz vectors and residuals in one
pass, H X formed only for a contraction; bit-identical when on), `BAANDI_NARROW_PROJ=on`
and `BAANDI_NARROW_REDUCED=on` (blocks of at most `NARROW_WIDTH` columns by streaming
dot kernels; rounding class), `BAANDI_B5=par` (B5 with root-parallel secular solve and
column-group parallel product; bit-identical to `BAANDI_B5=on` on faer).
All default off and bit-identical when off. Tests: `tests/b1_selected_roots.rs`,
`tests/b5_reduced_eigen.rs`; B4 harness `tests/b4_inner_iterations.rs` (ignored).
Plan and pre-registered rules: docs/performance/requests/b1-b4-b5-davidson.md.

## MPI entry-cost switches (2026-11-14)

`distributed_davidson::DistributedTuning` (read by `from_env`/`current`, forced in
tests with `solve_distributed_rows_tuned` or the per-thread `set_tuning_override`):
`BAANDI_M2_ORTHO=serial` (serial-order orthonormalization: one projection per pass,
DGKS rule, panel Gram-Schmidt with the serial panel rule; rounding class),
`BAANDI_M2_DGEMM=on` (`rotate_local_with`/`project_out_with` wide route; bit-identical
on pure builds), `BAANDI_M2_LEAN=on` (bit-identical slice loops). `BAANDI_M2_FFT=mkl`
selects MKL lanes in `DistributedPairFft` (`with_lanes`). All default off; they apply
to the M1 and M2 distributed solvers. Analysis and dev08 request:
docs/performance/mpi-entry-cost-analysis.md.

## MPI reduced eigenproblem switches (2026-11-15)

Same `DistributedTuning` mechanism: `BAANDI_M2_EIGEN=root` (`eigen_root`: rank 0
solves the reduced eigenproblem, one broadcast of status, the lowest `bands` values
and coefficient columns; bit-identical to the replicated solve; a failure on rank 0
is broadcast so no rank waits forever) and `BAANDI_M2_B5=on` (`bordered`: the serial
B5 secular update from `bordered.rs` when exactly one column was appended, own
capacity x capacity GEMM workspace; rounding class). Both default off.
`reduced_solve_stats()` gives the order distribution (mean, p90, max, sum of m^3,
single-column appends, bordered and broadcast solves), written to results.json
`parallel.reduced_solves`. Analysis: docs/performance/mpi-reduced-eigen-analysis.md.

## Lean Davidson (work order beat-QE phase B, 2026-11-15)

`BAANDI_DAVIDSON=lean|classic` (default `classic`; `ResearchOptions::lean`), module
`src/lean_davidson.rs` (child of `real_davidson`). Acts only on Ritz-value
(standard-profile) solves of the serial real Gamma path; residual (strict) solves, the
certificate, the complex and the MPI distributed solvers are unchanged. Per iteration:
residual vectors of the unconverged roots only; all-band Ritz vectors and H images at
return (the standard residual guard reads their norms) and contraction; unit-norm
corrections appended without projection; overlap S = V^T V kept beside H_r with the
explicit inverse Cholesky factor T (bordered, panels of 64 columns, GEMM workspace);
reduced problem A = T H_r T^T by dsyevd/faer or B5; coefficients T^T Y. Pivot rule:
every new column's squared pivot over its overlap diagonal must be at least
`PIVOT_RATIO_MIN = eps / 1e-12`. Failing appends: stage 1 re-projects the block alone
through T (rank floor `DEPENDENCE_TOLERANCE` on the unit input, then the classic in-block
Gram-Schmidt) and tests again; stage 2 rotates the basis to its Ritz basis once and
projects classically. Counters: `lean_counters()` (fallback and reprojection rates,
pivot statistics) and `flop_counters()` (GEMM flops by site for both paths), written to
results.json `convergence.work.davidson_lean`. Tests `tests/lean_davidson.rs`.
Do not make it the default without the dev08 A/B (docs/performance/attempts-ledger.md
row B-LEAN-DAVIDSON).

## Kept start (REUSE-PSI, 2026-11-16)

`RealDavidsonWorkspace::{keep_returned_products, kept_start, kept_start_mut,
arm_kept_start, invalidate_kept_start}` and `RealDavidson::prepare_workspace`:
a lean return copies X and its H images (the return rotation) into the
workspace; an armed kept start makes the next lean solve start from that X and
the caller-updated H X without the initial application, only when `initial`
equals the kept X to `KEPT_START_TOLERANCE` (1e-14); every other solve consumes
or drops it (classic and residual solves invalidate it). Counters in
`lean_counters()` (`kept_starts`, `kept_declined`, `kept_stored`, `kept_columns`).
Test `tests/lean_davidson.rs::kept_start_skips_initial_application`.

## Lean projection cache (C3a, 2026-11-16)

`BAANDI_LEAN_PROJ_CACHE=on|off` (default off; `ResearchOptions::lean_proj_cache`):
the lean solver keeps P^T V of its basis columns (`lean_davidson::ProjCache`;
loaded through `copy_last_projections` after every apply of basis columns,
rotated by the contraction coefficients, invalidated by the stage-2 fallback
rotation) and, at a return with a complete cache, leaves P^T X = (P^T V) C in
the workspace (`RealDavidsonWorkspace::returned_projections`, cleared at the
start of every solve). The REUSE-PSI kept start carries the kept P^T X. The solve
is bit-identical with the switch on and off; only the caller's nonlocal energy
from these projections is rounding class. Counters in `lean_counters()`
(`proj_cache_*`). Test `tests/lean_proj_cache.rs`; ledger C3a-LEAN-PROJ-CACHE.
