# ADR 0025: Composable residual preconditioning for density mixing

- Task: P10D-SCF-MIXING
- Date: 2026-11-03
- Status: implemented; proposed for reviewer acceptance

## Decision

Keep `DensityMixer`, `LinearMixer`, `PulayMixer` and their default arithmetic
unchanged. Introduce `ResidualPreconditioner`, `DensityTransform`, `Kerker<T>`
and `PreconditionedPulay<P>` in `dft-solvers/src/mixing/kerker.rs`.
The first interface separates a residual operation from history fitting; the
second lets orchestration provide the already adopted FFT without adding an
array-library dependency to the solver manifest. The concrete wrapper satisfies
the existing `DensityMixer` contract. Physical convergence remains exclusively
in SCF, using raw density residuals and the existing fixed-point certificate.

This is Anderson on the filtered map, including a filtered least-squares metric.
It is not Johnson Broyden and not raw-metric Pulay with a filtered final step.
See [the research note](../research/scf-mixing-strategies.md) for details,
conditioning, units, Nyquist convention, diagnostics and literature.

## Authorized integration

The orchestrator approved minimal additional edits to the solver module root
and CLI run wiring. `lib.rs` adds only a path-qualified public module declaration;
no Davidson code is touched. `run.rs` calls config's strategy parser and extends
strategy and history metadata. All other CLI validation and the new parameter
tests live in the originally owned config file.

`Mixing::KerkerPulay { beta, history, q0_inverse_bohr }` selects the new strategy.
Its match branch creates a reusable `FftPlan3d` and two arrays, captured in a
transform closure. This adapter is necessary because the solver manifest has
no ndarray dependency. G comes from the actual reciprocal cell; no geometry
crate dependency is introduced in solvers.

`MixerHistory` extends the existing trait with read-only stored-pair count.
SCF records actual history differences after the proposal, the raw residual
contraction ratio and cumulative feasibility resets. The new variant joins the
existing feasibility branch. There is no change to the fixed-point certificate,
eigensolver policy, raw convergence norm, growth thresholds or legacy mixing
arithmetic. Additional SCF edits are limited to imports, the enum and selection
branch, four record fields and their initialization, history inspection and
the feasibility-branch match/reset counter (including a test record literal).
The larger-than-single-line match branch is FFT adapter construction, not a new
SCF controller or numerical backend.

CLI keys: `mixing="kerker-pulay"`, required `mixing_q0_inverse_bohr` in Bohr^-1,
optional `mixing_history` (default 6), existing `mixing_beta` (default 0.3).
Invalid/inapplicable parameters fail explicitly. No screening length is inferred
from material identity. Existing input forms continue to select the same defaults.
No restart-on-divergence or parameter fitting is hidden in this selection.
The new strategy's deterministic boundary recovery is explicit: with no secant
history, a negative filtered proposal is replaced by a raw convex linear step
with min(beta,1). It is recorded as `mixing_raw_recovery`; the existing controller
still handles other infeasible extrapolations. This addresses measured stagnation
on multiple systems, not a Si-specific threshold. See the research note.

The solver tests import the exported module. The ignored `mixing_scf` test
executes the production CLI from Python, creates synthetic Gaussian-ion UPFs
using the existing fixture formulas, and records inputs, full output histories,
certificates and failures. This allows cross-layer regression without adding
manifest dependencies or changing engine test ownership.

## Consequences

No dependency, lockfile, baseline, ticket, artifact or README changes. Default
mixers remain byte-for-byte untouched. The FFT filter and wrapper scratch reuse
allocations; the existing Pulay history and eigensolve allocations persist.
An allocation-free replacement of that mixer is outside this decision.
No default-parameter recommendation follows from the synthetic response alone.
Real-space local screening, adaptive beta, rollback and potential mixing remain
separate candidates. Material convergence conclusions are restricted to the cases in the comparison
report. No metallic convergence claim is made. Defaults remain unchanged.
