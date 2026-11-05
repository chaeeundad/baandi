# Agent-native runtime

Status: reconciled 2026-11-05 at `b5d515a`. P11A–H runtime, checkpoint, diagnostics, CLI, MCP and Python/ASE implementations are merged. The original P11R semantic sketches below are design history, not exact public signatures; current exports are in `crates/dft-engine/src/lib.rs`. [ADR 0028](../decisions/0028-agent-native-runtime.md) records commitments; [roadmap](../roadmap/runtime.md) assigns gates.

## Placement and resolved overlaps

| Layer | Proposed home and existing seam | Ownership decision |
| --- | --- | --- |
| A numerics | Existing physics/solvers; `DensityMixer::mix/reset`, Davidson through `step_kpoints` | P10D owns mixer algorithms, preconditioning and feasibility; runtime selects validated strategies only |
| B state | Pure `dft-engine::state` types wrapping `IterationRecord`, `ScfResult`, `ScfOutcome` | ESLE `ElectronicState` is a payload in the single P11C/A4 checkpoint envelope, not a second persistence format |
| C step control | Pure `dft-engine::runtime` session wrapping `ScfProblem::run_impl` | No serialization, scheduler, clock, network or agent dependencies in numerical execution |
| C driver | `dft-cli::runtime` drives the session, events, persistence, cancellation and outer coupling | This is ADR 0017's coupling runtime layer, not a competing runtime |
| D integration | CLI-side adapters initially, separately ticketed optional external clients later | ASE, MLIP and external solvers never become numerical dependencies |
| E agent tools | P11F owns one server over the driver | P9A L4 and conversational UI consume that MCP server; no second server or natural-language parsing in engine |

ADR 0016's `save/load` sketch appears beside pure engine types: resolve it explicitly as CLI persistence, never engine file I/O. Its optional exact-trajectory concept is made explicit by the three checkpoint kinds in [state API](scientific-state-api.md). `run_from`/`InitialState` remain additive ESLE warm-start proposals, not an alternative session. ADR 0017's external-field and energy-ledger capabilities remain unimplemented; session extraction does not authorize them. P9A documents are now present: `user-interface.md` and `../roadmap/ui.md`; their proposed UI work is not all implemented.

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

Current source inventory (2026-11-05; implementation presence, not fresh numerical acceptance):

| Capability | Merged implementation | Limits |
| --- | --- | --- |
| Complete SCF iteration and three modes | `dft-engine/src/session.rs`, `events.rs`, `state_view.rs` | One session; certificate alone establishes convergence |
| Linear/Pulay/Kerker-Pulay and diagnostics | `scf.rs`, `dft-solvers/src/mixing/kerker.rs`, `diagnostics.rs` | Raw mixer residual and physical density distance have different units |
| Davidson/hybrid RMM and fallback | `dft-solvers/src/rmm.rs`, engine solver dispatch | Opt-in; buffer-band policies retain Davidson fallback; no raw-RMM runtime choice |
| Safe-point beta/reset/solver/tolerance/budget control | `runtime.rs`, `policy.rs` | Host approval or bounded launch authority; no weakening certificate gates |
| Pause/resume/cancel, checkpoint restore | `runtime.rs`, `checkpoint.rs`; CLI checkpoint adapter | Engine payload is in memory; binary persistence belongs to CLI; exact continuation is build/identity constrained |
| MCP | `crates/baandi-mcp` | Observe/Supervised, host-owned inputs and approval; no autonomous MCP authority |
| Python/ASE | `crates/baandi-py`, `python/` | Optional native extension; energy workflows, no forces/stress |

P10A/P10D/P7PC are merged dependencies. Switching algorithms intentionally changes
the trajectory and requires replay and scientific validation. Changes to physical
model, band count or buffer threshold are not authorized runtime policy changes.
See [evidence map](../validation/evidence-map.md) for executed slices and remaining
proof gaps. P11I coupling/MLIP contracts remain specifications, not delivered adapters.
