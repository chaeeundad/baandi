# Three execution modes over one SCF implementation

Status: design, not implemented. Contract source: `tickets/P11-COMMON.md`, Addendum 2026-11-03, read from main at `98f0b50` because this worktree's ticket predates it. The addendum overrides the earlier draft where they differ. [ADR 0022](../decisions/0022-agent-native-runtime.md) commits the boundaries; [control](runtime-control.md) and [events](scientific-state-api.md) specify protocols.

## One numerical implementation

`ExecutionConfig.mode` is a typed enum `Classic | Observed | AgentControlled`, selected at calculation creation. It is not a choice of eigensolver, physics, convergence test or SCF algorithm. All modes ultimately invoke the same extracted `ScfProblem::run_impl` iteration, `step_kpoints`, mixer and fixed-point certificate. Do not maintain a second observable loop or a third agent loop. Conditional observer/control work happens around that iteration boundary, not inside FFT, Hamiltonian or solver kernels.

| Mode | Execution policy | Available externally | Never implied |
| --- | --- | --- | --- |
| Classic | Drive the shared session to completion with observation/control disabled | Existing final `run/run_outcome` outputs and errors | No event stream, telemetry serialization, agent, live status or historical event reconstruction |
| Observed | Same progression plus bounded scalar event publication and on-demand snapshots | Read-only local/async subscribers; final result | No policy mutation, no subscriber-induced wait, no agent reply requirement |
| AgentControlled | Observed plus authenticated safe-point control requests | Supervised proposals by default; explicitly configured bounded autonomy; audit | No LLM-owned convergence, direct arrays/state mutation or inner-loop callback |

Classic may still keep the numerical history that `run_impl` already requires and returns. “No events” does not mean remove existing `IterationRecord` or alter error behavior. Use a disabled observer selected once during setup; do not allocate serialization buffers, queues or timestamp events for Classic. Runtime bookkeeping needed by the shared numerical implementation must remain within the measured overhead target.

## Fast, control and agent paths

| Path | Executes where | Inputs and outputs | Blocking rules |
| --- | --- | --- | --- |
| Fast | Existing numerical crates plus extracted session | Immutable model, current density/orbitals, validated solver/mixer; ordinary record and certificate | No network, LLM, JSON, file I/O, callbacks into external clients or new telemetry allocation |
| Control | Driver at complete iteration boundaries; bounded local diagnostics | Small records, state versions, typed requests; validated atomic changes and audit/events | Bounded CPU work; no waiting for agent/subscriber. Explicit user pause or storage/resource fault is a lifecycle event |
| Agent | Optional external process/service through P11F | Coalesced snapshots/evidence; explanations and proposed ControlRequest | May be slow/unavailable; has no direct numerical reference or mutation lock |

The LLM is outside the loop because reply latency is variable, reasoning may be stale, and model output is not a numerical invariant. If SCF takes 2 ms/iteration while an agent reply takes 3 seconds (illustrative values, not measurements), many iterations or the entire run may finish before the proposal arrives. Keep computing with the last approved policy. Retain only the newest unsent routine observation for that agent, allow one outstanding reasoning job per calculation, and rate-limit new jobs to at most one per second in the example adapter profile. These transport controls do not change the number or order of SCF iterations.

A reply must carry the version it observed and explicit expiry/preconditions; do not reinterpret it as a command for the current state. The driver rejects a terminal/expired request, or performs the bounded revalidation specified in [control](runtime-control.md). Do not pause every iteration to obtain an LLM reply. A user may explicitly pause for review, but that pause and its approval/application are separately audited. Agent timeout does not silently create a pause.

## Mode switches

Mode and numerical policy are different state. Any successful mode switch increments `state_version`, preserves numerical policy revision, density, orbitals, mixer history, budgets, certificate thresholds and RNG state, and emits a never-drop `mode_changed` event. The engine cannot retroactively create events for an unobserved interval.

| Source -> target | Before initialize | After initialize/start |
| --- | --- | --- |
| Same -> same | No-op config | Idempotent response; no fake transition |
| Classic -> Observed or AgentControlled | Allowed creation configuration | Rejected `mode_unavailable_after_start`; start a new calculation if observation is needed |
| Observed -> AgentControlled | Allowed creation configuration | Ready/Paused boundary only; authenticated operator approval, configured audit sink/control quotas/approval policy and available capacity required |
| AgentControlled -> Observed | Allowed creation configuration | Ready/Paused boundary only; operator-approved revocation; reject queued unapplied proposals, revoke autonomy, retain audit and stream |
| Observed or AgentControlled -> Classic | Allowed creation configuration | Rejected; cannot discard an established event/audit contract |
| Any different mode at Converged/Failed/Cancelled/Exhausted | Not applicable | Rejected; historical mode remains part of provenance |

During a Running step a mode request may be queued but not applied. On AgentControlled -> Observed, an already applied change remains part of numerical policy; switching off control does not roll it back. If paused, preserve pause state: explicitly resume under AgentControlled before downgrading if continued execution is wanted. Observed has no agent resume operation. Exact checkpoints retain the mode and stream sequence/version; restore cannot evade the switch table by rewriting metadata. Warm start creates a new calculation with an explicitly chosen mode.

Mode switching itself is an administrative driver operation, not an autonomous agent allow-list action. Observed -> AgentControlled approval uses the same authenticated review mechanism as other control approvals, but does not grant the agent permission to approve itself. Neither switch changes requested physical observables or the physical model.

## Ownership and adoption

P11A owns the single step implementation and mode config; P11B the control validator; P11F the one MCP adapter; P11H-CLI-UX owns CLI commands/sinks and human approval presentation. P11H's proposed `baandi run`, `--observe`, `--agent`, `status` and `explain` must preserve the actual current CLI commands/aliases; this document does not claim they exist. P9A consumes the same server, not another controller.

P11I-COUPLING-API-SPEC owns high/low-level scientific integration contracts; modes do not replace its API or put an agent between external numerical solver calls. P11J-RUNTIME-BENCH owns matched results in [runtime equivalence](../validation/runtime-equivalence.md). P11R supplies requirements and targets only. These three tickets are present on main `98f0b50`, not this checkout's `tickets/`; cross-references use ticket IDs until integration brings those files here.
