# Bounded subspace updates for intermediate SCF maps (BREAKTHROUGH-3-BOUNDED-SUBSPACE)

Status: design + serial prototype behind `BAANDI_SCF_ADVANCE=off|on` (default
`off`; bit-identical to main when off). Source of the idea: review of
2026-11-17 (`BAANDI_breakthrough_performance_review_2026-11-17.md`, section 7,
R3). Related failures that this design must not repeat: ledger A21 (fixed
tolerance factors 3/10/30/100 on the eigen schedule: SCF count rose, work paid
back), A22 (converged RMM variant: 2 to 10x more H columns), backlog B11
(single-sweep RMM untested because it cannot certify intermediate maps).

## 1. What changes and what does not

Today every SCF map m runs the eigensolver to the Ritz-change tolerance
tau_m of the standard ratchet (`StandardProfile::ratchet`, `session.rs`
`advance`). The proposal is a second, explicit operation contract for
intermediate maps:

- `SolveIntent::Converge` (today's contract, unchanged): return eigenpairs
  whose Ritz values stopped changing below tau_m (standard) or whose true
  residuals are below tau_m (strict).
- `SolveIntent::Advance(SubspaceBudget)` (new): run at most
  `correction_blocks` Davidson correction blocks (and at most `h_columns`
  operator columns) from the kept start and return the **approximate
  occupied subspace**: all-band Ritz vectors X, their Ritz values, their H
  images (kept in the workspace for the next start), the all-band residual
  norms and `SubspaceDiagnostics`. The budget is a deterministic count of
  work, never wall time.

The engine contract (ADR 0043 V1 to V6) is untouched:

- the stop rule, the residual guard, the verification attempt and the
  delivered state are evaluated only on `Converge` maps;
- the final convergence candidate and the verification always run the
  present exact solve;
- `verified = true` keeps its exact meaning; no tolerance, threshold or
  accuracy class changes; `BAANDI_SCF_ADVANCE=off` is bit-identical to main.

This is **not** "lower `max_iterations` and treat `NotConverged` as
success": the solver returns an explicitly typed `AdvancedSubspace` whose
`converged` flags are the true Ritz-change flags (most false), the session
records every such map as approximate (`IterationRecord.advance`), never
stops on it, and decides a completion (exact solve at the same density from
the advanced orbitals) by the rules in section 5. It is not A21 (no
tolerance is scaled: the exact maps keep the ratchet) and not A22 (no
converged RMM; the Davidson correction of the existing lean solver is used
unchanged, only its stopping contract differs).

## 2. Types

```rust
// dft-solvers (real_davidson.rs)
pub enum SolveIntent { Converge, Advance(SubspaceBudget) }
pub struct SubspaceBudget { pub correction_blocks: usize, pub h_columns: usize }
pub enum BudgetStop { Converged, CorrectionBlocks, HColumns }
pub struct SubspaceDiagnostics {
    pub correction_blocks: usize,  // blocks appended in this call
    pub h_columns: usize,          // operator columns applied in this call
    pub stop: BudgetStop,          // why the call returned
    pub kept_start: bool,          // started from the kept X / H X (no initial apply)
    pub unconverged_roots: usize,  // roots whose Ritz change is still above threshold
}
pub struct AdvancedSubspace { pub solution: RealEigenSolution, pub diagnostics: SubspaceDiagnostics }
impl RealDavidson { pub fn advance_reusing(.., budget, ..) -> Result<AdvancedSubspace, SolverError> }
```

```rust
// dft-engine
pub struct AdvanceRecord {          // per IterationRecord (None on exact maps)
    pub budget_blocks: usize,
    pub blocks_used: usize, pub h_columns: usize,
    pub inner_estimate_ha: f64,      // B = sum f_i |r_i|^2 of the advanced state (Ha)
    pub scf_error_ha: f64,           // E_H[rho_out - rho_in] of the advanced map (Ha)
    pub completion: Option<CompletionTrigger>,  // Some = an exact solve followed at the same n
    pub completion_h_columns: usize,
}
pub enum CompletionTrigger { Candidate, InnerNoise, Growth, GapInstability }
```

The approximate state is never a public `ElectronicState` variant: inside the
session it is an ordinary `ScfStep` whose `IterationRecord.advance` is
`Some`, and it never leaves the session (the delivered `ScfResult` is the
exact converged step, V4).

## 3. Consumers: exact eigensolution or approximate state?

| consumer | needs | with Advance |
|---|---|---|
| density build `rho = sum_i f_i |psi_i|^2` | an orthonormal occupied block X with fixed integer f | approximate X is fine: P = X X^T is invariant under X -> X U (U orthogonal); empty buffer bands have f = 0 and never enter |
| energy audit (kinetic, nonlocal, Hartree, XC, double counting) | orbitals and eigenvalues consistent with each other | Ritz pairs satisfy theta_i = x_i^T H x_i exactly, so the band-sum identity `sum f_i theta_i = sum f_i <x_i|H|x_i>` and the explicit audit agree as for an exact solve; the energy of an approximate map is reported in the history only |
| Pulay history (rho_in, F~(rho_in)) | F~ close to F | enters the history only when the inner estimate is small against the SCF residual (section 5, rule R-INNER); otherwise the map is completed and the exact pair enters |
| standard stop estimator E_H[rho_out - rho_in] and residual guard | an exact map | evaluated on exact maps only; an approximate map never stops (rule R-CAND completes it when it would) |
| verification (ADR 0043 V1-V6) | exact warm start from the stop step | unchanged: the stop step is always exact, the verification is the exact solve; `VERIFY_REUSE=half` acts after a failed verification, and the next map is exact by rule R-CAND (previous E below the threshold) |
| forces, stress, checkpoint, observe/agent, spectral audit, NSCF | the delivered exact state | they see only the delivered converged step (exact) or a checkpoint of the session (approximate maps are not special: the mixer history and density are ordinary; the record is additive) |
| REUSE-PSI kept start, H X update | X and H X of the previous return | the Advance return stores X and H X exactly like a lean return; the next solve (Advance or Converge) takes `H_new X = H_old X + (V_new - V_old) X` |

Fractional occupations (Fermi-Dirac, collinear FD), several k points, the
complex solver, the strict reference and the MPI distributed solver fall back
to the exact solve with a recorded reason (`advance_declined` counters with
the reason; the record is `None`). Scope of the prototype: closed-shell fixed
integer occupations, one exact-Gamma k point, real serial lean Davidson
(`BAANDI_DAVIDSON=lean`, the CLI default).

## 4. Budget control (observables only)

Let E_m = E_H[rho_out - rho_in] of map m (Ha, the stop estimator) and
B_m = sum_i f_i |r_i|^2 (Ha, the residual guard's estimator) of the state
that produced rho_out. The exact path accepts, through the ratchet
tau = 0.1 E_prev / N_e, a summed band-energy error of order N_e tau = 0.1
E_prev. The Advance rule uses the same constant:

- R-INNER: an advanced map is accepted when B_m <= ADVANCE_INNER_RATIO * E_m
  with ADVANCE_INNER_RATIO = 0.1 (the ratchet's own factor; not a new
  tuning constant). Otherwise the inner error is comparable to the SCF
  residual (the Pulay secant would learn solver noise): completion.
- Budget start: `BAANDI_SCF_ADVANCE_BLOCKS` (integer 1..=8, default 2)
  correction blocks per map; `h_columns` = blocks x bands (a block corrects at
  most every band).
- Adaptation (deterministic, from the previous advanced map): after a
  completion by R-INNER the budget grows by one block (up to 8); when
  B_m <= 0.1 * ADVANCE_INNER_RATIO * E_m the budget shrinks by one (down to
  1). Both bounds are fixed integers of the prototype.

## 5. Completion and fallback triggers (exact solve at the same density)

A completion is `step_session(n, warm = advanced X, tau_m, ..)`: the exact
solve starts from the advanced orbitals (REUSE-PSI takes the kept X and H X,
the potential is unchanged, so no H update is needed) and its work is
charged to the step like the first-iteration redo.

- R-CAND (convergence candidate): E_m < threshold_ha, i.e. the stop rule
  would fire on an approximate map. Also: the map is run exact (no Advance)
  when the previous step's E was below the threshold (the stop fired and the
  verification failed, or `confirm` territory) or when no history exists
  (iteration 1 and its redo, the initialization of the subspace).
- R-INNER: B_m > 0.1 E_m (section 4).
- R-GROWTH: E_m > E_{m-1} (the SCF residual grew under an approximate map;
  complete to separate solver error from mixer behaviour).
- R-GAP: the highest occupied Ritz value is not below the lowest empty Ritz
  value (occupation/gap instability of the approximate state).
- Any solver error of the Advance call (fallback stage exhaustion, pivot
  failure) falls back to the exact solve with a recorded reason.

## 6. Mixer interaction

e_inner = F~(rho) - F(rho) is unobservable without a second solve. R-INNER
bounds it through B_m (an estimator, not a bound: ADR 0043, review F03) so
that only pairs whose inner error is small against the SCF residual enter the
Pulay history; the record marks them (`advance` is `Some`), and research
runs can compare the per-map ratio B_m / E_m against exact maps (the exact
maps' ratio is recorded too). Periodic diagnostic full solves are a research
option, not part of the prototype. The mixer itself is unchanged (no
re-weighting in the prototype; if the screen shows the mixer learning noise,
the next step is to drop the pair whose inner ratio exceeded the rule).

## 7. Ladder, determinism, MPI

- Cutoff ladder coarse stage: the coarse session is an ordinary session, so
  the switch applies there too (its stop is on the density residual, which
  approximate maps report as every map does).
- Determinism: the budget is a count; the trajectory is reproducible for a
  fixed thread count. Thread-count independence of results is the existing
  rounding class of the lean solver (GEMM reductions), unchanged.
- MPI: the distributed lean algebra (`M2_LEAN_ALGEBRA`) has the same loop
  structure; the Advance contract would be added there after the serial
  screen passes (every budget decision is replicated; the completion trigger
  uses replicated scalars). Not in this prototype: MPI maps decline with the
  reason `distributed solver`.

## 8. Risks and quantitative gates

Risks: (1) SCF step count grows (the mixer sees an inexact map); (2) a
completion after every advanced map doubles the potential/density/audit
work; (3) the first correction block corrects every band (the Ritz-change
flags are NaN at iteration 1), so a 1-block advance costs B columns and may
not reduce the residual enough; (4) gap crossing in approximate states.

Gates (screen, one case, T=1, full SCF to convergence incl. verification,
each run under 10 minutes): continue only if total H columns (recorded +
certificate) fall by at least 5 percent against the default with converged
and verified true, |dE| <= 5e-9 Ha against the default, and the SCF step
count does not grow by more than the H-column saving pays for (wall inside
the spread or better). Then three seeds, then one skewed/slab case. Failure
of the first gate is reported with the per-map ratio and completion counts.

## 9. Result of the serial prototype (2026-11-17, Mac observations, counters exact)

Implemented behind `BAANDI_SCF_ADVANCE=off|on` (default off; bit-identical to
main when off: si16 s601 density SHA and energy equal the main binary's).
Research sub-switches `SCF_ADVANCE_BLOCKS` (1..=8, default 2),
`SCF_ADVANCE_GROWTH`, `SCF_ADVANCE_ADAPT`, `SCF_ADVANCE_SCOPE=all|fine`.
Diagnostics: `IterationRecord.advance` (per map: budget, blocks and H columns
used, kept start, unconverged roots, B, E, ratio, completion trigger and its H
columns) in results.json `convergence.history[].advance`, and the process
counters `convergence.work.scf_advance`.

Stage-1 screen (full SCF to convergence incl. verification, T=1, pure build,
PseudoDojo Si.upf with `xc.functional = "pbe"`, f1e-3 harness inputs; total H
columns = recorded + certificate; all runs converged and verified, |dE| against
the default 2.5e-10 to 7.5e-10 Ha):

| case | default | on (defaults) | growth off | growth off, adapt off, b2 | b3 | b2 fine only |
|---|---|---|---|---|---|---|
| si16 s601 | 981 (SCF 9, ver 2) | 1359 (11, 3) | 1101 (9, 2) | 1237 (11, 3) | 1016 (9, 1) | 1046 (10, 2) |
| si16 s602 | 824 (9, 1) | 1125 (10, 2) | - | 818 (8, 1) | 1208 (10, 2) | 1223 (11, 3) |
| si16 s603 | 1183 (10, 2) | 1172 (10, 2) | - | 1170 (11, 2) | 1107 (10, 1) | 952 (9, 1) |
| si-skew s601 | 127 (6, 2) | 155 (7, 3) | - | 140 (6, 2) | 145 (6, 2) | - |

Verdict: the gate of section 8 (at least 5 percent fewer H columns on one
case, then three seeds) is not met by any variant; the idea is rejected as a
speed change (ledger BT3-ADVANCE-CONTRACT, BT3-GROWTH-TRIGGER,
BT3-BUDGET-VARIANTS, BT3-COARSE, BT3-INNER-RATIO). What was learned:

1. Per advanced map the saving is real but small: 42 to 71 H columns against
   53 to 107 for the exact warm maps at the same SCF stage. The standard
   ratchet already stops exact warm solves after 3 to 7 Davidson iterations,
   and the first all-band block (B columns, Ritz-change flags undefined
   before the second Rayleigh-Ritz) dominates both. Bounding the block count
   cannot remove that first block; review R2 (residual-rank selection of the
   first correction) is the premise that would change this.
2. The completions pay the saving back: a candidate completion costs 36
   columns, an inner-noise completion 60 to 84; after a large Pulay step (E
   dropped three orders) a 2-block advance leaves B / E of 1.5 to 5.9.
3. The growth trigger misreads ordinary Pulay non-monotonicity (the exact path
   itself grows at iteration 3 of si16 s601) as inexactness: 4 of 7
   completions were growth completions of maps with ratio 0.008 to 0.03.
4. The mixer did not visibly learn the inner noise: E per iteration tracks
   the exact path while B <= 0.1 E. What moves is the endgame (one more SCF
   step, one more verification attempt in two of four geometries), i.e. the
   trajectory, not the per-map cost.
5. The coarse ladder stage is 3 to 11 percent cheaper under Advance, but the
   fine stage pays the worse hand-off back (BT3-COARSE).
6. si-slab from the harness inputs does not converge on the default path with
   either UPF (iteration budget exhausted); it is not a usable screen case
   for this idea and was not judged.

The MPI contract (section 7) was not implemented: distributed solves decline
with the reason `distributed solver`. The code stays as an explicit, tested,
switchable contract (`advance_reusing`, `SolveIntent`) for a future first
block that is cheaper than an all-band correction.
