# ADR 0020: Opt-in RMM-DIIS refinement and Davidson recovery

- Task: P10A-RMM-DIIS
- Date: 2026-11-03
- Status: implemented, proposed for reviewer acceptance

`SolverChoice` defaults to `Davidson`; its `Hybrid(HybridRmm)` variant is an
explicit caller choice. `ScfOptions::eigensolver` selects both serial and parallel
k-point paths, including the final fixed-point certificate. Existing Davidson
functions, tolerance schedules, CLI defaults and dependencies do not change.

`RmmDiis::solve` takes an initial full-width block, optional existing
`Preconditioner`, and per-band absolute residual tolerances. It returns the
existing `EigenSolution`. `RmmError` distinguishes upstream numerical/application
errors from rank loss, stagnation and iteration-budget failure; refinement
failures retain a last certified full-rank state when available. Repeated
eigenvalues are not an error; dependent vectors are. A successful return requires
fresh operator residuals for every band and entrywise Gram defect <= 1e-10.

`HybridRmm::solve` screens supplied starts by their true Rayleigh residuals. It
uses tolerance-based Davidson warm-up when needed, then RMM; rank loss or
stagnation/iteration exhaustion permits one Davidson fallback. Terminal Davidson nonconvergence retains its `HybridPhase`, total attempted
operator columns and source residuals. Underlying application errors are never
hidden by fallback. All applications, including
screening, unsuccessful RMM work and certification, are counted. Histories from
completed phases are concatenated; iteration counts count their recorded Ritz
sweeps, not screening or inner trial steps. Standalone hybrid telemetry identifies
warm-up and fallback; the engine adapter returns the existing summary.

The current Davidson cannot return approximate vectors on budget exhaustion,
so warm-up uses a tolerance, not a fixed number of Davidson iterations. It solves
all bands at the strictest caller tolerance on fallback; the RMM API supports
strict/loose bands but this base has no merged P7PC engine policy. No speculative
buffer policy is added. Allocation-free RMM workspace adapters were authorized
additively in dft-numerics; existing adapter functions and tests stay unchanged.
Legacy Davidson warm-up/fallback allocations are outside that refinement claim.

RMM can miss a lower eigenstate absent from an invariant supplied guess, like
Davidson (ADR 0006). Tests against dense references and complete SCF comparisons
remain essential. Near-degenerate wanted boundaries can stagnate explicitly;
this is why fallback and accurate error reporting are part of the API.

See `docs/research/rmm-diis-design.md` for provenance and implementation choices,
and `docs/performance/rmm-diis-benchmark.md` for measured limits. No universal
crossover or production switching recommendation follows from complexity alone.
