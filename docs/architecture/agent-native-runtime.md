# Agent-native runtime

Status: design only, P11R, 2026-11-03. Code inspected at `425b0f2984be0e05598e2cd6d67fa19394a33553`; no runtime interface described here is implemented. [ADR 0028](../decisions/0028-agent-native-runtime.md) records commitments; [roadmap](../roadmap/runtime.md) assigns gates.

## Placement and resolved overlaps

| Layer | Proposed home and existing seam | Ownership decision |
| --- | --- | --- |
| A numerics | Existing physics/solvers; `DensityMixer::mix/reset`, Davidson through `step_kpoints` | P10D owns mixer algorithms, preconditioning and feasibility; runtime selects validated strategies only |
| B state | Pure `dft-engine::state` types wrapping `IterationRecord`, `ScfResult`, `ScfOutcome` | ESLE `ElectronicState` is a payload in the single P11C/A4 checkpoint envelope, not a second persistence format |
| C step control | Pure `dft-engine::runtime` session wrapping `ScfProblem::run_impl` | No serialization, scheduler, clock, network or agent dependencies in numerical execution |
| C driver | `dft-cli::runtime` drives the session, events, persistence, cancellation and outer coupling | This is ADR 0017's coupling runtime layer, not a competing runtime |
| D integration | CLI-side adapters initially, separately ticketed optional external clients later | ASE, MLIP and external solvers never become numerical dependencies |
| E agent tools | P11F owns one server over the driver | P9A L4 and conversational UI consume that MCP server; no second server or natural-language parsing in engine |

ADR 0016's `save/load` sketch appears beside pure engine types: resolve it explicitly as CLI persistence, never engine file I/O. Its optional exact-trajectory concept is made explicit by the three checkpoint kinds in [state API](scientific-state-api.md). `run_from`/`InitialState` remain additive ESLE warm-start proposals, not an alternative session. ADR 0017's external-field and energy-ledger capabilities remain unimplemented; session extraction does not authorize them. P9A's three design documents are absent in this checkout; reconciliation uses [its ticket](../../tickets/P9A-USER-INTERFACE-DESIGN.md), especially L4. Recheck the delivered P9A documents before implementation.

No new crate before the stepwise API passes equivalence. Consider splitting a pure state crate only with two independent consumers and stable versioned types; consider a runtime adapter crate when CLI and Python/MCP demonstrably duplicate driver logic. Each split needs dependency/license/Apple Silicon review, an ownership update and a new ADR, without upward dependencies from numerics.

## Proposed additive session interface

These are semantic signatures, not Rust declarations or current public APIs:

| Method | Input -> output | Wraps |
| --- | --- | --- |
| `new` | validated `ScfProblem`, driver limits -> Created session | `ScfProblem::new/with_kpoints`, existing option checks |
| `initialize` | optional validated start payload -> Ready snapshot | `run_impl` uniform density, mixer validation/reset, warm orbitals and counters |
| `step` | mutable session -> immutable StepReport or typed failure | One entire `run_impl` iteration using `step_kpoints`, including applicable mixing |
| `snapshot` | shared session -> small read-only StateSnapshot | Already computed `IterationRecord` and session counters |
| `result` | shared session -> `NotReady`, `Unavailable(reason)` or borrowed `ScfOutcome` | Existing `run_outcome` result semantics, no new calculation |
| `apply_policy` | safe-point revision + typed patch -> audit decision | Validated `ScfOptions` subset and `DensityMixer::reset` |

`step_kpoints` currently computes potentials, eigenpairs, output density and energy; it is NOT a complete SCF iteration. Extract the loop locals (`n`, `previous`, mixer, history, best residual, growth count, iteration and schedule state) without reordering arithmetic. One public step performs prechecks, scheduled solve, record creation, divergence checks, optional certificate, terminal-result construction or density mixing/feasibility damping and orbital retention, in that order. The safe point is after this completes, not between solve and mix.

Certificate work is nested in the same public step: call `step_kpoints` on the output density with warm orbitals and `0.1 * eigen_tolerance`. It adds no ordinary iteration, has separate operator counts, and cannot be paused or patched midway. Existing candidate gates are density/energy/eigen/requested-eigen tolerances plus input and output relative charge errors <= 1e-8. Certificate requires density residual <= density tolerance and occupied eigenvalue change <= 1e-7 Ha, with existing gap checks. Only this path sets `converged`.

Preserve the subtle terminal behavior: certified steps and the final retained `run_outcome` step skip mixing; `run` currently mixes its last unsuccessful iteration before returning `NotConverged`. Session construction must retain an explicit legacy completion mode (`strict` or `retain_last`), not silently unify the trajectories. A terminal step reports `mixing_applied=false` when skipped. Certificate failure continues to mixing as today. An eigensolver error remains a typed failure; no invented completed record or partial convergence.

Keep `run` and `run_outcome` unchanged until [equivalence](../validation/runtime-equivalence.md) passes; only then make them thin session drivers. All control and checkpoint operations occur at complete iteration boundaries. Initialization and each step may be expensive; cancellation acknowledges the request promptly in the driver but waits for the boundary. No partial Davidson restart is promised.

The eigensolver seam is `dft-solvers/src/lib.rs::Davidson::solve` (delegating to `solve_with_recurrence/solve_with_restart`); the session retains per-k orbital guesses already held by `run_impl`, not a second Davidson implementation. `dft-solvers/src/mixing.rs::PulayMixer` currently stores private `(input,residual)` pairs and a mean; exact checkpoint support must version these semantics. The runtime beta upper bound is a new control-layer restriction: current solver validation requires only finite positive beta, which is why the legacy beta=100 negative test remains valid without becoming an allowed runtime policy.

## Three-mode extension and implementation inventory

[Execution modes](execution-modes.md) defines Classic/Observed/AgentControlled as policies over this one session. [ControlRequest](runtime-control.md) is the only asynchronous intervention envelope; external agents never call `step_kpoints`, `DensityMixer::mix` or solver methods directly. Compile/runtime capabilities come from the actual merged build, not a ticket label or branch name.

Read-only inventory taken 2026-11-03: this worktree HEAD is `425b0f2`; main was `98f0b50`; P10A/P10D worktree HEADs were `0142ebf`, with feature changes uncommitted. P7PC HEAD `18a5547` is an ancestor of main, but inspection of that commit's source does NOT contain the new `BandTolerances` API; the implementation below exists in P7PC's dirty working tree. A merged commit with a descriptive title is not proof its later working-tree changes were merged. These findings are a source inspection, not rerun validation of those tasks.

| Intervention | Existing implementation inspected | Merge state and runtime decision |
| --- | --- | --- |
| Set beta / reset history | This checkout `dft-solvers/src/mixing.rs`: `LinearMixer`, `PulayMixer::new`, `DensityMixer::reset`; `run_impl` constructs/owns mixer | Present here/main; session mutation wiring is still proposed. Reset the exact active mixer, not a replacement algorithm written by runtime |
| Choose Linear/Pulay | `scf.rs::Mixing::{Linear,Pulay}` and matching `run_impl` branches | Present; same density representation and charge/positivity safeguards required |
| Choose Kerker-Pulay | P10D `scf.rs::Mixing::KerkerPulay{beta,history,q0_inverse_bohr}`; `dft-solvers/src/mixing/kerker.rs::{Kerker,PreconditionedPulay,MixerHistory}`; engine provides FFT/cell/grid to `Kerker::new` | Uncommitted P10D implementation, not in this checkout/main. Advertise only after merge and capability tests; q0 positive Bohr^-1, no filter recreation from agent code |
| Mixer health | P10D `MixingDiagnostics::{residual_norm,contraction,history_used,resets,raw_recoveries,used_raw_recovery}`, `PreconditionedPulay::diagnostics` | Uncommitted. Raw Euclidean residual norm (electron number/Bohr^3, no quadrature) is NOT `density_distance` (Bohr^-3/2); label and convert only with a stated definition |
| Choose compatible eigensolver | P10A `dft-solvers/src/rmm.rs::SolverChoice::{Davidson,Hybrid(HybridRmm)}`, `HybridRmm::solve`, `RmmDiis::solve`; proposed `ScfOptions::eigensolver` dispatches in `step_kpoints` | Uncommitted. Runtime exposes only davidson/hybrid after merge, not a raw-RMM choice merely because `RmmDiis` exists. Preserve warm-orbital shape and strict-band/certificate contracts |
| Eigensolver fallback | P10A `HybridSolution::{phase,fallback}` and `RmmFailureReason`; hybrid implementation owns fallback to Davidson | Uncommitted. Runtime reports existing fallback, does not implement another DIIS or replace solver residual tests |
| Buffer-band accuracy | P7PC `BandTolerances{strict_bands,buffer_tolerance}`, `Davidson::solve_with_band_tolerances_and_certification`, `EigenSolution::buffer`, certificate `buffer_residuals_per_k` | Dirty working-tree implementation only at inspection. Occupied + first empty band strict; code's buffer threshold `max(tolerance,1e-5)` Ha is observed implementation, NOT a new runtime setting or recommendation |
| Tolerance tightening | Existing `ScfOptions::{density_tolerance,energy_tolerance,eigen_tolerance}`, candidate/certificate checks in `run_impl` | Present scalar fields; runtime validator/floors not implemented. Preserve existing certificate and no weakening |
| Pause/resume/cancel/checkpoint | Loop boundary after record/certificate/mixing; retained `ScfOutcome` | No public session/control or checkpoint implementation here; P11A/B/C supply it, no direct kernel hook |

P10A/P7PC compatibility is not automatic: the inspected `SolverChoice::solve` constructs an all-band tolerance vector, whereas P7PC distinguishes strict/buffer bands. Until the merged solver preserves the agreed strict-band and certification semantics, deny live solver selection with `incompatible_solver_policy`; continue with the original solver. Switching algorithms is an intentional trajectory change and needs replay plus final scientific validation, not bit equality to a different algorithm. Reset algorithm-specific cache only after the new choice validates; leave physical state unchanged on rejection. Changes to band count, subspace, physical model or buffer threshold remain outside the autonomous first slice.

Inspection fingerprints (SHA-256 of dirty files, not commits or test evidence):

| Worktree and relative path | SHA-256 |
| --- | --- |
| P10D `crates/dft-solvers/src/mixing/kerker.rs` | `199e6214ff70e2957c840d0d4eceeef77597dfc4b4806ce37bec0325e5dac37f` |
| P10D `crates/dft-engine/src/scf.rs` | `6488ca963c59bff77b42557b2bac8791347795fb0586c7a04d2d0d67c72e3b0f` |
| P10A `crates/dft-solvers/src/rmm.rs` | `69bf6c9b8fbe762324ad78e0c62a5043bca09efdd48d3b919d30bcf16a8df72b` |
| P10A `crates/dft-engine/src/scf.rs` | `3ecbfdd2dce8d0a69b174efe50724e3625fc33faa3c25820e1cc619d19301a94` |
| P7PC `crates/dft-solvers/src/lib.rs` | `f4ade1684c57414bbcd693181fb4524057e858ab486f9dc26d95409716bb7cf7` |
| P7PC `crates/dft-engine/src/scf.rs` | `dfdcb11d98fa87913e5f76b4fd1d9ce0c9fd14930e63690599c1b177ad893de4` |

These parallel worktrees are active and may change immediately after inspection. P11A must re-read merged code and refresh this mapping before implementation. No files in those worktrees were modified by P11R. P11H CLI consumes these capabilities, P11I owns scientific coupling specifications, and P11J measures mode overhead; none should infer readiness from an unmerged implementation.
