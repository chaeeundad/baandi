# P10A: RMM-DIIS design and implementation

Status: implemented, experimental opt-in; default remains Davidson. Base
`0142ebf`, branch `task/P10A-RMM-DIIS`, 2026-11-03. Numerical workspace ownership
was explicitly extended by the orchestrator after the initial design review.

## Approved additive numerical boundary

Existing functions and tests are unchanged. Added in `linalg.rs`:

- `gemm_into`: borrowed column-major faer views and a caller-owned result,
  sequential backend multiplication, shape/finiteness validation.
- `HermitianWorkspace::{new,solve}`: fixed-dimension diagonal storage and faer
  scratch allocated at construction; borrowed matrix and caller-owned values and
  eigenvectors on every solve. Hermiticity checks match the existing adapter.
- `orthonormalize_into`: two-pass MGS with the existing projection arithmetic,
  caller-owned output and scratch, retained-rank return and zeroed trailing
  columns. It does not need an external basis for RMM's complete-band sweeps.
- `reusable_evd_gemm_and_rank_filter` and
  `reusable_adapters_rectangular_and_empty_shapes` in `linalg/tests.rs`: repeated
  agreement with existing adapters, pointer reuse, rank filtering, invalid
  inputs, rectangular adjoints and empty matrices. Empty scratch construction
  bypasses a backend query that assumes positive dimension.

These are sufficient without mutable-width matrix views: the small DIIS matrix
is padded to the fixed history capacity, and inactive indices are excluded from
its affine constraint. All RMM vectors, histories, rotation and DIIS storage,
scalar arrays and history capacity are allocated before iteration. No logging,
timers, I/O or explicit allocations occur inside refinement loops. Existing
Davidson warm-up/fallback retains its allocation-owning API. Backend-private
packing behavior is not measured by a global allocation interceptor (which would
require an unsafe global allocator); pointer reuse tests cover adapter storage.

`lib.rs` only adds the new module and exports. Engine changes are one typed
`ScfOptions::eigensolver` field, its Davidson default, and selection at the two
existing serial/parallel solve sites. No CLI/manifest/default-tolerance changes.
See ADR 0020 and the benchmark document for the exact evidence and limits.

## Sources and provenance

Papers and public algorithm/user documentation were accessed; no QE, VASP or
GPAW implementation source was read.

- [Pulay (1980), Chemical Physics Letters 73, 393–398](https://doi.org/10.1016/0009-2614(80)80396-4):
  DIIS accelerates iterative sequences by combining their error information.
  The publisher full text was unavailable (HTTP 403). This is not a claim to
  have studied the full paper.
- [Kresse and Furthmüller (1996), Physical Review B 54, 11169–11186](https://doi.org/10.1103/PhysRevB.54.11169):
  applies Pulay residual minimization to iterative diagonalization, seeking to
  reduce cubic work. Publisher PDF request returned HTTP 401; an
  [author-uploaded full text](https://www.researchgate.net/publication/225188263_Efficient_iterative_schemes_forab_initiototal-energy_calculations_using_a_plane-wave_basis_set)
  was subsequently found and section III B reviewed. It uses a normalized
  residual minimization with a small generalized residual/overlap problem,
  preconditioned trial updates and inter-sweep orthogonalization/rotation. It
  explicitly warns that incomplete initial eigenspaces can omit states.
  Our implementation choices below are not a transcription of that algorithm.
- [GPAW public residual-minimization description](https://gpaw.readthedocs.io/documentation/rmm-diis.html):
  describes orthogonalization, subspace diagonalization, preconditioned residual
  steps and residual-norm step-length minimization. Our norm-conserving case has
  identity overlap; its PAW details are not imported.
- [QE public input description](https://www.quantum-espresso.org/Doc/INPUT_PW.html#diagonalization):
  the accessed page identifies itself as version 7.5 and documents RMM alternated
  with Davidson or ParO for stability. It is not implementation evidence for
  QE 7.6 and its parameter defaults do not justify BAANDI defaults.
- Internal [Davidson/SCF study note](qe-study/davidson-and-scf-thresholds.md),
  especially measured outcomes: smaller Davidson subspaces, a fuller diagonal
  preconditioner and eigenvalue-change convergence did not improve the measured
  engine cases. This design does not revive those rejected hypotheses.

## Reusable engine contracts

`LinearOperator::apply_block` already gives Hamiltonian batching without exposing
physics to the solver. Check output shape and finiteness even for overrides.
Count actual vector-column applications, including trial directions, residual
certification, warm-up and fallback; do not report only the successful phase.

`Preconditioner::apply` already accepts one Ritz value per residual column and
overwrites a caller-owned output. Its diagonal implementation returns a signed
correction, so the solver must not negate it a second time. Without a
preconditioner use negative residual directions. Residual formation and norms
are physics glue; keep their work buffers persistent.

`EigenSolution` supplies eigenvalues, vectors, per-band residuals, iterations,
maximum-residual history and operator applications. Existing helper routines
return owned allocations and cannot simply be called from the new hot loop.
The existing two-pass MGS dependence threshold is `1e-12`; existing eigensolver
tests require entrywise Gram defect below `1e-10`. Keep these tests and bounds.

Current `scf.rs` builds one Davidson configuration and calls it on both serial
and parallel k-point paths. Its tolerance schedule and final fixed-point
certificate must remain unchanged. P7PC is present as a separate worktree but
its per-band tolerance API is absent from this base. Do not assume that ticket's
session status means it has been merged.

## Published method versus this variant

The 1996 band-history minimization uses a normalization metric and chooses an
initial trial length from a Rayleigh-quotient search. Our variant instead uses
real affine DIIS coefficients (sum one), a regularized residual Gram inverse,
fixed Ritz shifts within a sweep, and a fresh residual least-squares trial length.
This is a residual-minimization/Pulay-DIIS variant, **not an exact reproduction of
the 1996 method**. The choices use the repository's Hermitian workspace adapter
without a new generalized-eigenproblem API. The residual rejection, rank checks,
true certification and Davidson fallback bound failure behavior. Benchmark
conclusions apply to this implementation; they do not establish the performance
of other RMM-DIIS formulations. Published parameter values are not adopted.

## Implemented algorithm (experimental implementation choices)

Let X contain m bands, H be fixed for a solver call, and r_j = Hx_j - e_j x_j.
Use blocked operator applications but separate short histories for each band.
This retains batching while avoiding a single growing Davidson space. A fully
coupled DIIS history would have different storage and conditioning costs and is
not the first implementation.

1. Require a finite, correctly shaped initial block or obtain one from Davidson.
   Two-pass MGS must retain m independent columns. Rotate by the lowest m
   eigenvectors of X* H X before forming band histories.
2. Hold the Ritz shifts fixed during a short inner refinement. For each band,
   form p_j from its preconditioned residual, apply H to the direction block,
   and set q_j = (H - e_j) p_j. The real least-squares trial length is
   alpha_j = -Re(q_j* r_j) / (q_j* q_j). Reject zero/nonfinite denominators or
   nonfinite steps. Any step cap is a measured implementation parameter, not a
   value borrowed from another code.
3. For each band retain pairs of trial vectors and their fixed-shift residuals
   in a fixed-capacity ring. Minimize the norm of the residual combination with
   coefficients summing to one. The residual Gram matrix is Hermitian positive
   semidefinite. A backend Hermitian EVD regularizes inverse metric eigenvalues at
   1e-12 after scaling by its maximum diagonal; real coefficients use the real
   residual Gram matrix. Reject nonfinite coefficients or magnitudes above 10,
   and fall back to a single safe trial when the history is singular.
4. Apply the same coefficients to trial vectors and residuals; evaluate proposed
   states before acceptance. Do not mix histories computed under different
   Hamiltonians or Ritz shifts. Reject residual increases under a documented
   acceptance rule and restart the affected band history.
5. After the inner block, orthonormalize all m bands, rotate by Rayleigh–Ritz
   again and discard histories: band identity and phase may have changed.
   Evaluate fresh H X for certification. Recompute all residuals using the
   returned values and returned vectors, including previously converged bands.
6. Success requires each residual <= its own positive finite absolute tolerance,
   ordered eigenvalues and Gram defect <= `1e-10`. Rank loss, stagnation or the
   iteration budget produces an explicit typed failure with per-band residuals
   from the last certified full-rank state and a flag if none exists.

Cholesky orthogonalization would be efficient on a well-conditioned Gram matrix
but needs a condition/rank guard; Loewdin adds an m-dimensional eigensolve.
Two-pass MGS follows the repository's validated rank policy and is the initial
choice. None removes the O(n m^2) inter-band orthogonalization or O(m^3) rotation.
Reducing the frequency trades these costs against duplicate-state and convergence
risks; benchmark it rather than claiming the cubic cost disappears.

## Degeneracy, initial guesses and hybrid behavior

Equal eigenvalues are legitimate. Detect duplicate *vectors* through rank/Gram
tests, never by deleting repeated eigenvalues. Test a wanted-set boundary inside
an exactly degenerate cluster. Compare eigenvalues and subspace projectors for
complete clusters, rather than arbitrary eigenvector phases. At a partial
cluster boundary, compare distance to the full degenerate eigenspace and rank.

A small residual only certifies an eigenpair, not that it belongs to the lowest
m eigenvalues. An invariant initial block can omit lower states, even with
perfect orthogonality. Expose this limitation as in ADR 0006 and include such
guesses in tests. Warm-up should inject independent exploratory components for
untrusted guesses; it is not a mathematical proof of spectral completeness.

Implemented hybrid policy: Davidson to an explicitly configured warm-up residual
threshold, then RMM refinement; reuse the previous SCF block when its true Rayleigh residuals are below that
threshold. Orthogonality alone is insufficient: the first engine experiment
showed large unnecessary RMM work with coordinate starting vectors. Screening
applications are included in the hybrid operator count.
Current Davidson errors do not return approximate vectors, so a fixed-count
warm-up cannot recover a failed solve's last iterate. A tolerance-based warm-up
is possible without editing Davidson. Record this difference in the API.

On RMM stagnation or rank loss, call Davidson once using the last safe full-rank
block (or the original valid guess if none). Do not mask operator, dimension or
nonfinite errors as numerical stagnation. Propagate warm-up/fallback nonconvergence with
`HybridPhase`, total attempted operator columns and the source per-band residuals. Return telemetry identifying which phases ran and
include all work in totals. Fall back to the strictest requested band tolerance
until Davidson supports a compatible per-band API; this is conservative but can
be expensive. Strict/loose buffer policies are supplied by the caller, and the
SCF certificate must tighten every band.

No implicit environment-variable selection. The additive `SolverChoice::{Davidson,Hybrid}` option affects both k-point
branches and defaults to Davidson (ADR 0020). `HybridRmm::solve` additionally
returns warm-up, fallback and RMM sweep telemetry. The engine uses the existing
`EigenSolution` summary and does not expose phase telemetry in CLI output.

## Validation and limits

Tests cover true residuals, per-band tolerances, exact degeneracy with a partial
wanted cluster, wide spectra, typed stagnation/rank/iteration failures, malformed
operator overrides, fallback work counts and repeatability. A near-degenerate
partial boundary with spacing 1e-7 stagnates around 5.1e-8 residual and is recovered
by Davidson without relaxing the 1e-9 requirement. This is a measured RMM
limitation, not a successful RMM-only convergence claim.

Experimental defaults are history size 4 (validated configurable range 2..=16),
80 outer sweeps, 10 stale sweeps and warm-up residual 1e-3 in operator units.
Each sweep has up to history-size inner steps. Staleness means less than 1%
improvement of the maximum absolute residual; it never authorizes success.
These parameters are conservative trial controls, not measured universal optima.
No production default is changed or recommended. Real and synthetic comparisons
are in [the benchmark document](../performance/rmm-diis-benchmark.md).

At most one fallback is attempted. Operator/backend errors propagate without
retry. A refinement failure contains the last full-rank certified state with
per-band residuals, or explicitly no such state on initial rank loss. The
selection adapter maps terminal errors to the existing engine solver error
boundary. No residual certificate proves that an invariant supplied guess
contains every lowest state; this existing ADR 0006 limitation applies to both
solvers and must not be interpreted as a ground-state completeness guarantee.
