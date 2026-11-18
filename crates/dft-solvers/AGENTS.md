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

Fast fixed-operator tests: `cargo test -p dft-solvers --test all real_davidson::`.
Parallel dispatch and timing observation (ignored in debug):
`cargo test --release -p dft-solvers --test all real_davidson:: -- --include-ignored --nocapture --test-threads=1`.
Evidence and limitations: `docs/performance/real-gamma-r4-notes.md`.

## Research switches of the real Davidson (2026-11-10)

`RealDavidson::solve_research` takes `ResearchOptions` explicitly; every other
entry point reads them from the environment (`davidson_stages.rs`): stage timers
`BAANDI_DAVIDSON_STAGES=on`, B1 selected-root products `BAANDI_B1=standard|strict|on`
with `BAANDI_B1_SKETCH=count`, B5 bordered eigen-update `BAANDI_B5=on` (`bordered.rs`),
B6 cached projector products for the return certificate `BAANDI_B6=on` (operators opt in
through the `RealLinearOperator` default methods; test `tests/all/b6_certificate.rs`).
Micro switches read at construction: `BAANDI_WIDE_DOT=on` (eight-lane MGS dot, rounding
class); see docs/performance/requests/micro-stack.md.
Algebra-stream switches (docs/performance/requests/algebra-stream.md, test
`tests/algebra_stream.rs`): `BAANDI_RITZ_FUSED=on` (Ritz vectors and residuals in one
pass, H X formed only for a contraction; bit-identical when on), `BAANDI_NARROW_PROJ=on`
and `BAANDI_NARROW_REDUCED=on` (blocks of at most `NARROW_WIDTH` columns by streaming
dot kernels; rounding class), `BAANDI_B5=par` (B5 with root-parallel secular solve and
column-group parallel product; bit-identical to `BAANDI_B5=on` on faer).
All default off and bit-identical when off. Tests: `tests/b1_selected_roots.rs`,
`tests/b5_reduced_eigen.rs`; B4 harness `tests/all/b4_inner_iterations.rs` (ignored).
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

## MPI solver workspace (review F06, 2026-11-16)

`BAANDI_M2_WORKSPACE=on|off` (default on; `DistributedTuning.workspace`): the
distributed solver runs on `distributed_davidson/workspace.rs`, a per-thread (per
rank) `DistributedDavidsonWorkspace` reused across iterations and solves: basis and
products with `capacity` columns of storage (appends in place), reduced matrix with
the fixed leading dimension `capacity` and border update (packed m x m copy for the
dense solver), fresh products, lowest coefficients and Gram operands as slices,
persistent Ritz/residual/correction/orthonormalization scratch, GEMM, EVD and B5
workspaces (rebuilt only when (max rows, capacity, pool) changes). Same arithmetic
in the same order: bit-identical to off (tests `tests/distributed_workspace.rs`, ranks
1 to 4, every switch stack; `dft-engine` `m2_workspace_is_bit_identical_in_scf`).
Growth counters `workspace_stats()`; `release_workspace()` frees the calling
thread's workspace; `set_iteration_probe` is a test hook. The classic (non-serial)
orthonormalization and the M1 apply still allocate inside; on pure builds faer's EVD
allocates its own packing (MKL `dsyevd` does not). `dft-mpi`: `tree_sum_f64` and
`allgatherv_f64` keep their receive scratch and count/displacement arrays per thread;
the rank simulator reuses its slots (results unchanged).

## MPI selected-root rotation (review F05, 2026-11-16)

`BAANDI_M2_SELECTED_ROOTS=on|off` (default off until dev08; `DistributedTuning.selected_roots`,
implies the workspace path): on Ritz-change (standard) distributed solves the
unconverged roots are decided from the Ritz-value change before any rotation; V C,
(H V) C and residuals are formed for those columns only (coefficient columns gathered
into `sel_coeff`); all bands are rotated at return (every band's residual reported)
and at a restart (when a restart follows, one full rotation serves both the selected
residuals and the contraction basis). Strict (residual) solves are unchanged.
`rotation_stats()` counts rotation pairs, columns and flops on both paths. Rounding
class in principle (column subsets of one GEMM); bit-identical on pure builds in
`tests/distributed_selected_roots.rs` (ranks 1 to 4, scattered, clustered and
degenerate spectra, capacities 6 and 3). History records the selected maximum per
iteration (as the serial lean solver).

## MPI reduction counters (review F08, 2026-11-16; analysis only)

Both distributed entry points wrap the communicator in `Counted` (forwards every call
unchanged) and count each tree sum by the site last declared on the thread with
`site(Site::..)`: reduced_gram, residual_norms, ortho_normalize, ortho_projection,
panel_cross, panel_dot, panel_norm, restart_defect, apply (the operator's own sums),
setup; payload histogram buckets 1, 2-8, 9-64, 65-512, 513-4096, above
(`REDUCTION_BUCKETS`). Count only (no timer, no allocation); reduction order and
results unchanged. `reduction_stats()` (solves, sites), results.json
`parallel.reductions`, printed by the m2_entry_cost harness.

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

## MPI in-panel Gram-Schmidt reductions (2026-11-16)

`BAANDI_M2_PANEL_GS=classic|pythag|fused` (default classic; `DistributedTuning.panel_gs`,
`PanelGs`): `panel_column` orthogonalises one column against the kept in-panel
columns for both the workspace and the allocating path. classic = the earlier
arithmetic bit for bit (two dot reductions, one norm); pythag = the norm rides with
the second dots (Pythagorean update, explicit norm when it would cancel more than a
factor 4); fused = the input norm rides with the first dots and the second pass runs
only when the per-column DGKS rule (`BAANDI_ORTHO` keep) fails. Rounding class.
Tests: lib `panel_gs_tests`, dft-engine `m2_panel_gs_variants_keep_classes`.
Profile and request: docs/performance/profile-2026-11-16-mpi.md,
docs/performance/requests/mpi-fast-screen-2026-11-16.md. `dft-mpi` also has
`BAANDI_MPI_TREE=rd` (recursive-doubling tree sum, bit-identical for power-of-two
rank counts; set by the CLI with `dft_mpi::set_tree_sum`).

## MPI kept start (MPI-KEPT-START, 2026-11-16)

`RealDavidson::solve_distributed_rows_kept` with `DistributedKeptStart` (workspace
path only): returns the local rows of X in the caller's `x_local`; a Ritz-change
return with `set_keep(true)` stores H X; an armed solve starts from X and the
caller-updated H X (no orthonormalization, no initial apply) when one tree sum
finds `initial` equal to X to `KEPT_START_TOLERANCE` on every rank and the Gram
defect of X is at most 1e-12 (the restart contraction rule). Flags must be set
under replicated conditions. `DistributedPairFft::apply_local_batch_kept` takes
kept real-space slabs for the leading vectors (bit-identical to
`apply_local_batch` for the same slabs). Tests `tests/all/distributed_kept_start.rs`,
lib `kept_planes_batch_equals_plain_batch_bitwise`.

## MPI scaling gap (MPI-SCALING-GAP, 2026-11-17)

`BAANDI_M2_STICKS=greedy|mirror` (default greedy; per-thread override
`distributed_fft::set_mirror_sticks_override`): `partition` assigns a stick and
its Gamma mirror stick to one rank as a unit, so no pair representative of the
real Gamma operator is split and its two mirror alltoalls per apply vanish
(rounding class). `BAANDI_M2_START=ortho|gram` (default ortho; override
`distributed_davidson::set_start_gram_override`, counters `start_gram_stats`):
the workspace solver takes an orthonormal start (Gram defect <= 1e-12 from one
reduction, the restart contraction rule) as its basis without the panel
Gram-Schmidt (rounding class). Solver stages `solve_start` and `solve_return`
time the entry and the N x B gather. Tests: lib
`mirror_partition_keeps_gamma_pairs_on_one_rank`, `tests/distributed_start_gram.rs`.
Analysis docs/performance/mpi-scaling-gap.md.

## MPI lazy orbitals (MPI-LAZY-ORBITALS, 2026-11-17)

`RealDavidson::solve_distributed_rows_local` with `DistributedStart::{Global, Local}`
and `gather: bool`: a `Local` start is this rank's rows in the operator's
`local_rows` order (no replicated block; bit-identical to the same start passed
`Global`); with `gather = false` the workspace path skips the final N x B
allgatherv and returns a 0 x bands eigenvector placeholder, the rows living only in
`x_local`. The non-workspace path gathers a local start and always returns the full
block (callers check `eigenvectors.rows()`). `solve_distributed_rows_kept` is this
entry with `Global` and `gather = true` (unchanged). Test
`tests/all/distributed_lazy_rows.rs` (1 to 4 ranks, workspace full/selected and the
allocating path, bits). Engine design: docs/design/mpi-lazy-orbitals.md.

## Incremental Pulay update (MPI-REMAINING, 2026-11-17)

`BAANDI_MIX_INCREMENTAL=on|off` (default off; per-thread override
`mixing::set_incremental_override`, test hook `PulayMixer::set_incremental`):
`PulayMixer` keeps its difference vectors (and metric images) and Gram entries
across calls, keyed by the serial numbers of the pairs and seeds they are formed
from, and forms only the newest difference and Gram row per step. Reset,
`seed_differences` and `restore` renumber and drop the cache. Same operands,
formulas and orders: bit-identical (test `mixing::incremental_tests`). Counters
`mixing::mixer_counters()`. Analysis docs/performance/mpi-remaining.md.

## Per-core kernels (PERCORE-KERNELS, 2026-11-17)

`BAANDI_LEAN_BORDER=separate|fused` (default separate; `ResearchOptions::lean_border_fused`):
the lean solver applies H to the unit corrections E before the pivot test and forms the
overlap border V^T E and the reduced border V^T (H E) as one product V^T [E | H E]
(`lean_davidson::fused_border`, the reduced step then copies its border in
`reduced_from_fused`). A failing pivot test discards the images (`fused_wasted_columns`),
takes the separate fallback and ends fusing for that solve; only blocks with 4 k <= m fuse
(traffic rule: the copies of E and H E against the saved stream of V). Inactive with `BAANDI_NARROW_REDUCED=on`. Bit-identical
on the aarch64 pure kernel, rounding class on MKL / x86 faer. `BAANDI_RITZ_FUSED=loop`
(`ritz_fused_loop`): the fused Ritz/residual kernel only at the selected-root residuals of
the lean loop (bit-identical). Counters `lean_counters().fused_borders`,
`fused_wasted_columns`. Tests `tests/lean_davidson.rs::percore_*`. Analysis
docs/performance/percore-kernels-2026-11-17.md.

## MPI overlap-matrix Davidson (MPI-OVERLAP-DAVIDSON, 2026-11-17)

`BAANDI_M2_LEAN_ALGEBRA=on|off` (default off; `DistributedTuning.lean_algebra`, implies
the workspace path): Ritz-change (standard) distributed solves run the serial lean
algebra (`distributed_davidson/workspace/lean.rs`, reusing `lean_davidson.rs`'s
`LeanScratch`, `extend_factor`, `reduced_border`, `ritz_coefficients`): raw
corrections appended without Gram-Schmidt, ONE tree sum for the S border (with the
selected residual norms), unit scaling from the S diagonal, pivot test
(`PIVOT_RATIO_MIN`) before H is applied, the serial stage-1/stage-2 fallbacks and
rebuild; H border at the next iteration (`reduced_gram`). Two tree sums per
iteration instead of about 4 + k. Replicated S, T, A, C in the per-rank workspace
(own capacity x capacity GEMM workspace). `M2_EIGEN=root` is ignored on this path;
kept and Gram starts are accepted by the pivot rule (factor kept). Strict solves
unchanged. Counters `lean_algebra_stats()`, reduction site `lean_border`. Tests
`tests/distributed_lean_algebra.rs` (ranks 1 to 4: class against the orthonormal
path, rank decisions equal to the serial lean solver, fallbacks, bits),
dft-engine `tests/m2_lean_algebra.rs`. Design docs/design/mpi-overlap-davidson.md.

## Measurement hooks (BREAKTHROUGH-0-INSTRUMENT, 2026-11-17)

`distributed_davidson::STAGES` has an eleventh stage `restart` (the Gram defect
and contraction copy of the workspace path, the basis replacement of the lean
contraction; a re-orthonormalization at a restart counts as `orthonormalize`,
its application as `apply_total`, so the stages stay exclusive). `stage_add` and
`davidson_stages::Timers::stop` also add to the per-thread mirror
`dft_numerics::stage_clock::local` (one relaxed load when the engine's profile is
off). `DistributedKeptStart::{armed, stored, with_stored}` and
`RealDavidsonWorkspace::{kept_start_armed, load_kept_start}` serve the engine's
snapshot capture and replay only; no solver branch reads them.

## History-preserving restart (BREAKTHROUGH-1-RESTART, 2026-11-17; REJECT, default off)

`BAANDI_RESTART_KEEP=off|1|2|4|8|b/8|b/4|b/2` (`RestartKeep`, `ResearchOptions::restart_keep`,
`DistributedTuning::restart_keep`): at a lean contraction (serial lean solver and the
distributed lean algebra) keep up to r previous Ritz vectors of the unconverged roots
minus their X component, P = V (D - C A) and H P = (H V)(D - C A) from the contraction
rotation (no H application, no extra reduction). History and the coefficient-space
selection live in `LeanScratch` (`KeepScratch`: `keep_record`, `keep_identity`,
`keep_directions`); rank floor `DEPENDENCE_TOLERANCE`, pivot rule `PIVOT_RATIO_MIN`, the
new basis [X | P] is factored with the same rule and falls back to X alone. Off is
bit-identical, and so is any solve without a contraction. Counters
`lean_counters().keep_*`, `lean_algebra_stats().{contractions, keep_*}`. Tests
`tests/restart_keep.rs` (full grids ignored in debug; release command in the file).
Results: ledger RESTART-KEEP-*, design docs/design/restart-keep.md.
## Advance intent of the lean solver (BREAKTHROUGH-3, 2026-11-17)

`SolveIntent::{Converge, Advance(SubspaceBudget)}`, `RealDavidson::advance_reusing`
and `lean_davidson::solve_with_intent`: an `Advance` call runs at most
`correction_blocks` correction blocks / `h_columns` operator columns of the lean
Ritz-change loop from the kept start (or `initial`) and returns an
`AdvancedSubspace` (all-band Ritz vectors, values, residual norms, the true
Ritz-change flags, `SubspaceDiagnostics`); X and H X are kept as at a lean
return. `Converge` is the ordinary `solve`, bit for bit (the budget check never
fires). Classic Davidson or the second-vector experiment decline with
`InvalidParameter`. Test `tests/lean_davidson.rs::advance_respects_its_block_budget_and_converges_when_allowed`.
Engine use and the negative stage-1 result: docs/design/bounded-subspace-scf.md.

## Other labels and the layout cache (BREAKTHROUGH-4-OTHER-AND-INIT, 2026-11-17)

`distributed_davidson::STAGES` has 16 entries: `setup` (row-index exchange and
workspace preparation), `precondition`, `append` (exclusive), `wait`
(`BAANDI_PHASE_WAIT=on`: the `Counted` wrapper puts a timed barrier before every
collective; not exclusive) and `solver_total` (container of the solver call).
The per-thread mirror has exactly 16 distributed-stage slots
(`stage_clock::local::DIST_STAGE..M2_STAGE`): a further stage needs a new slot
range. `BAANDI_M2_LAYOUT_CACHE=on|off` (default off; per-thread override
`set_layout_cache_override`, counters `layout_cache_stats`): `solve_rows_impl`
reuses the previous solve's `Layout` (an `Arc`) when every rank's local rows
are unchanged, decided by one 1-double tree sum (replicated); bit-identical.
