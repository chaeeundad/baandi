# Runtime control contracts

Status: P11B control and P11C continuation are merged (2026-11-05); semantic design below wraps `ScfProblem::run_impl`, `ScfOptions`, `Mixing`, `DensityMixer::reset` and `ScfOutcome`. [State/event schemas](scientific-state-api.md) are normative companions.

## Lifecycle

Commands outside this table return `InvalidTransition` without mutation. Only one driver holds the session mutation lock. `Running` includes solve, certificate and mixing; no control operation enters numerics halfway through a step.

| From | Trigger | To | Effect |
| --- | --- | --- | --- |
| Created | initialize succeeds/fails | Ready/Failed | Allocate and validate once |
| Ready | step | Running | Consume next iteration |
| Running | ordinary step completes | Ready | Publish record after mixing |
| Running | certificate passes | Converged | Publish certified outcome |
| Running | iteration budget ends | Exhausted | Retained outcome only in retain_last mode |
| Running | numerical error | Failed | Preserve available history, never certify |
| Ready | pause | Paused | No numerical change |
| Running | pause request | Running | Queue one request; at boundary enter Paused unless terminal |
| Paused | resume | Ready | No new initialization |
| Ready/Paused | cancel | Cancelled | Preserve last complete boundary |
| Running | cancel request | Running | Acknowledge pending; at boundary enter Cancelled unless terminal |
| Ready/Paused | valid policy patch | same | Increment revision, reset mixer if required |
| Exhausted | extend budget with retained continuation | Paused | Explicit extension beyond done count; otherwise new warm-start calc |
| Ready/Paused/terminal | snapshot or export | same | Read only; exact checkpoint needs continuation state |

Terminal = Converged, Exhausted, Failed, Cancelled. Terminal numerical success/failure takes precedence over queued cancel/pause; record both request and disposition. Cancel wins over pause at a nonterminal boundary. Pending requests are not automatically applied: use the ControlRequest version/expiry/precondition checks below. The low-level PolicyPatch still binds the actual commit boundary/revision. Failed or cancelled sessions do not resume; restart is a separate explicit checkpoint request. A strict failure may have history but no `ScfResult`; result access returns `Unavailable`, not fabricated density/energy.

## Policy v1.0

Full `Policy`: `mixing:{kind:linear/pulay/kerker_pulay,beta:f64,history:u64|null,q0_inverse_bohr:f64 (Kerker only)}`, `density_tolerance:f64` (Bohr^-3/2), `energy_tolerance_ha:f64`, `eigen_tolerance_ha:f64`, `max_iterations:u64`, `eigensolver:davidson/hybrid` (hybrid capability-gated). These wrap the matching `ScfOptions` fields; history is required for Pulay/Kerker-Pulay and null for linear. All floats finite, `0 < beta <= 1`, `history >= 1`, tolerances > 0, `max_iterations >= iteration_done`. For a resumable next step require strictly greater than done count. Host quotas additionally cap history at 64 and total iterations at 10000 for the initial adapter profile; these are security limits, not scientific defaults, and may be lowered by deployment.

`PolicyPatch`: `schema_version:string`, `calc_id:string`, `request_id:string`, `expected_revision:u64`, `apply_after_iteration:u64`, `actor_id:string`, `reason_code:stagnation/oscillation/divergence_hint/manual/budget`, `changes:object` containing only the above fields, and `reset_history:bool`. Patches with no changed fields and reset_history=false are rejected; reset_history=true alone is valid. Driver authentication supplies actor identity; callers cannot impersonate it. Revision and completed iteration must match exactly. Unknown fields are rejected, never ignored.

Allowed: beta/strategy/history changes, explicit history reset, budget extension within quota, tolerance tightening (new <= old). Loosening any tolerance is rejected in v1, including through checkpoint import; no `unsafe_relaxed` escape hatch. Immutable for a live session: XC/functional override, pseudopotentials/species/electron count/structure, k mesh, cutoff/grid, occupations/gap threshold, bands/extra bands, seed, threads, eigensolver subspace/schedule and divergence guards. These require a new calculation and recorded lineage. No policy field can set convergence, certificate thresholds, charge safeguards or positivity rules.

Every mixer change resets all mixer history, even beta-only changes, via `DensityMixer::reset` or construction of a validated replacement. No-op changes are rejected. Tolerance/budget-only changes preserve history; explicit reset clears it. P10D owns preconditioner parameters, conservation, feasibility and diagnostic definitions. The proposed v1 tagged mixing schema additionally reserves `kind:kerker_pulay` with `q0_inverse_bohr:f64>0` and history required; only a build advertising the merged compatible P10D capability may accept it. The merged engine implements Kerker-Pulay; adapters must still advertise and validate the actual capability. Reject q0 on linear/pulay. No wire identifier alone grants a capability. Buffer-band recovery is a new calculation, not a hidden live patch.

Validate the complete prospective policy and allocate replacements before commit. Persist audit intent before mutation, then atomically commit policy/revision/reset and outcome record; storage failure prevents applying it. A rejected request leaves numerical state, mixer history and revision unchanged, while audit/event sequence advances. Retrying the same request id and canonical payload returns its recorded decision; different payload with that id is a conflict. Replay uses recorded boundary order, never wall time.

## Audit v1.0

Fields: `schema_version:string`, `calc_id:string`, `request_id:string`, `seq:u64`, `timestamp_utc:string`, `actor_id:string`, `actor_kind:human/rule/agent/replay`, `iteration_done:u64`, `expected_revision:u64|null`, `result_revision:u64`, `reason_code` (patch enum or `invalid_request`), `diagnostic_ids:string[]` (at most 16), `requested_patch:PolicyPatch|null`, `control_request:ControlRequest|null`, `approval_id:string|null`, `observed_state_version:u64|null`, `validated_state_version:u64`, `application_state_version:u64|null`, `request_bytes_sha256:string`, `before:Policy`, `after:Policy`, `decision:accepted/rejected/applied/duplicate`, `rejection_code:null|invalid_schema/invariant/restricted/stale/quota/conflict/transition/io/precondition_failed/expired/stale_policy_revision/future_state_version/terminal/approval_denied/approval_expired/prepare_failed/mode_disallows_control/incompatible_solver_policy/superseded`, `history_action:preserved/reset`, `previous_record_sha256:string|null`, `record_sha256:string`. Malformed requests use null patch/revision, a driver-assigned request id if necessary, and a bounded raw-input digest; never store nonfinite or executable raw payloads in typed audit fields. Rejection has after=before; hash excludes its own field and uses canonical schema field order and round-trip floats. Host authentication and protected append-only storage establish trust; a hash chain alone does not authenticate an attacker-controlled log. Free-text justification is optional external metadata, never policy input.

## Determinism

Promise is scoped to same model, exact build/backend/CPU execution identity, seed, completion mode, requested/effective workers and backend thread counts, same initial arrays and same applied-control log at the same boundaries. Thread setting is fixed at initialization (`ScfOptions::threads`: zero maps to one, capped at k count). Preserve canonical k reduction order. Require bit equality of numerical history, density, energy and certificate; timestamps, timings, process IDs and storage paths are excluded. Existing `tests/threading.rs` cross-thread identity is an additional regression, not a blanket cross-platform guarantee. Record requested and effective threads; changing them on exact restore is incompatible. Live LLM timing is not deterministic; reproducibility means replaying its recorded decisions. See [equivalence tests](../validation/runtime-equivalence.md).

## ControlRequest v1.0: asynchronous proposals, not commands from text

Only AgentControlled admits agent numerical/lifecycle control; read-only queries belong to Observed. Mode transitions use the administrative mechanism in [execution modes](execution-modes.md). The prior PolicyPatch is the internal numerical transaction payload; external callers send this envelope, not an unversioned patch. Types inherit [state API](scientific-state-api.md) finite-number/ID rules.

| Field | Type | Validation |
| --- | --- | --- |
| `schema_version` | literal `1.0` | Unknown version rejected |
| `calc_id`, `request_id` | string | Authenticated scope; idempotency is id plus canonical payload hash |
| `observed_state_version` | u64 | Version of the snapshot used for this decision |
| `observed_policy_revision` | u64 | Must still match for numerical changes; changed policy always rejects stale proposals |
| `action` | `set_beta/reset_history/select_mixer/select_eigensolver/tighten_tolerances/extend_budget/pause/resume/cancel/checkpoint` | Closed allow-list; no arbitrary method name |
| `parameters` | tagged object, table below | Unknown fields rejected; explicit units |
| `preconditions` | at most 8 tagged predicates | Evaluated against current committed state, all must pass |
| `expiry` | `{max_iteration:u64|null, deadline_utc:string|null}` | At least one non-null; reject if done>max_iteration or now>=deadline |
| `stale_policy` | `reject/revalidate` | Default reject; revalidate only for allow-listed actions/predicates |
| `reason_code` | `stagnation/oscillation/divergence_hint/eigensolver/resource/manual/budget` | Explanation only; never sufficient authority |
| `evidence_id` | string or null | Read-only evidence reference; cannot replace runtime checks |

Authenticated actor/session metadata is assigned by the driver, not trusted from the JSON body. Limit control wire body to 64 KiB and nesting depth 8. Wall expiry is converted at receipt to a monotonic deadline using the host clock; store that conversion. On exact resume, reapprove/revalidate pending proposals rather than extending an expired approval. Replays consume recorded application boundaries/dispositions and do not consult current wall time; they still validate physical invariants.

| Action | Exact parameters | Apply semantics |
| --- | --- | --- |
| `set_beta` | `{beta:f64}` | Validate 0<beta<=1 and deployment range; replace mixer configuration and reset history |
| `reset_history` | `{}` | Reset currently selected mixer's history; does not change density |
| `select_mixer` | `{mixing:MixingSpec}` | Tagged Policy mixing object; capability and shape/charge/positivity compatibility checks; reset history |
| `select_eigensolver` | `{solver:davidson/hybrid}` | Merged P10A capability only; preserve model/band counts/tolerances; reset incompatible solver cache, retain compatible orbital guess |
| `tighten_tolerances` | `{density_tolerance:f64|null,energy_tolerance_ha:f64|null,eigen_tolerance_ha:f64|null}` | At least one non-null; floor<=new<=old; no certificate threshold modification |
| `extend_budget` | `{max_iterations:u64}` | done<new<=host maximum; new>current budget |
| `pause/resume/cancel` | `{}` | Lifecycle table; terminal completion wins over pending request |
| `checkpoint` | `{kind:exact_continuation/warm_start/result_export}` | Capture immutable boundary payload under quota; no caller path; report write completion separately |

Hard floors are positive immutable launch fields `tolerance_floors:{density_residual:f64,energy_ha:f64,eigen_ha:f64}` and must satisfy floor<=initial tolerance. They limit excessive tightening/work, not authorize weakening. Values require backend/system evidence and are recorded, not universal scientific constants invented here. The original positive-tolerance and no-relaxation invariants always apply; increasing a floor cannot increase an active tolerance. Without validated floor metadata the adapter rejects tolerance interventions (ordinary configured SCF still runs).

Typed preconditions: `state_is{value:Ready/Paused}`, `policy_revision_is{value:u64}`, `mode_is{value:AgentControlled}`, `diagnostic_active{rule_id:string,rule_version:"1.0",minimum_consecutive:u64}`, `scalar_bound{field:density_residual/charge_error_relative/energy_change_abs_ha,op:lt/le/gt/ge,value:f64}`, `capability_is{value:linear/pulay/kerker_pulay/davidson/hybrid/exact_checkpoint}`. No expression evaluator, string predicate, shell, arbitrary JSONPath or user-supplied code. Missing/insufficient diagnostic evidence makes a predicate false, never presumed true. Predicates cannot read partial Running state or bypass candidate/certificate checks.

## Staleness and worked decisions

A version greater than the current committed version rejects as `future_state_version`. By default observed_state_version must equal the current committed version. `revalidate` permits a different state version only if policy revision is unchanged, the request is unexpired, each current predicate is true, and approval/allow-list explicitly permits revalidation of that action. Recompute diagnostics from current bounded history; a referenced old diagnostic is not enough. The validator records both observed and validated versions and evidence IDs. Prepare a low-level PolicyPatch against the CURRENT exact revision/iteration, then atomically compare-and-commit while holding the mutation lock. Any intervening mutation causes rejection, not silent retry.

1. Oscillation gone: agent observes iteration 12, version 13, revision 0, active `oscillation`, and proposes beta 0.15 from 0.3, expiry iteration 20, `stale_policy=revalidate`, predicate `diagnostic_active(oscillation,1.0,1)`. Reply arrives after iteration 18, version 19; residuals now contract and the last six energy differences no longer alternate. Runtime recomputes the rule, rejects `precondition_failed`, leaves beta/history/revision unchanged and emits `control_rejected`. Do not reset a successful mixer merely because an old alert remains in the agent context.
2. Oscillation persists: same request arrives at iteration 16, version 17 with active predicate and unchanged revision. A preauthorized autonomy rule permits this beta decrease, or an operator approves this exact revalidation scope. Runtime emits `control_accepted(validated_version=17)`, commits at that boundary, resets history and increments revision to 1/version to 18, then emits `control_applied`. Subsequent step 17 sees beta 0.15. These numbers illustrate version accounting, not experimental results.
3. Late or conflicting reply: at iteration 21 the above expiry rejects `expired` even if oscillation persists. If another approved action already advanced revision to 1 at iteration 14, reject `stale_policy_revision` at iteration 16. A proposal received after Converged is rejected `terminal` even when its timestamp deadline has not elapsed. Never revive a finished SCF to apply a late answer.
4. Manual pause: a running calculation may have advanced since the user viewed it. With explicit revalidation, unexpired lifecycle authority and current nonterminal state, pause may be approved for the next boundary; it cannot interrupt Davidson/certification halfway. Cancel received while that step converges is rejected `terminal` and the converged result survives. Paused calculations have no iteration expiry progress; wall deadline still applies.

## Approval and bounded autonomy

Default `approval_mode=Supervised`. Agent proposals remain `pending_approval`, do not reserve the numerical mutation lock and do not pause SCF. Approval is a host-authenticated record `{approval_id,request_digest,approver_id,allowed_revalidation:bool,expires_utc}` bound to one exact payload; another actor/agent cannot approve itself. Pending proposals expire normally; UI approval of a stale request does not skip revalidation. Denial emits rejection with `approval_denied`. Mode changes and autonomy-policy installation/removal always require operator approval. Explicit operator-originated lifecycle requests may carry that authenticated approval in the original request.

Autonomous operation is opt-in via a validated policy file, never “trust the agent”:

```json
{
  "schema_version": "1.0",
  "approval_mode": "PolicyBounded",
  "max_applied_total": 3,
  "fallback": "continue_current",
  "rules": [{
    "id": "reduce-beta-on-oscillation",
    "action": "set_beta",
    "ranges": {"beta_min": 0.05, "beta_max": 0.3},
    "must_not_increase": true,
    "required_diagnostic": "oscillation",
    "max_count": 3,
    "cooldown_iterations": 6,
    "max_state_version_lag": 8,
    "max_request_age_iterations": 8,
    "allow_revalidation": true
  }]
}
```

The full allow-list schema has `schema_version`, `approval_mode:Supervised/PolicyBounded`, `max_applied_total:u64`, `fallback:continue_current/local_rules_v1/pause_at_boundary`, and `rules:Rule[]` (at most 16). Each Rule has the example's count/cooldown/lag/age fields plus `action`, `id`, `required_diagnostic:string|null`, `allow_revalidation:bool`, and an action-specific `ranges` object. `set_beta` uses finite beta_min/beta_max within (0,1] and must_not_increase; `select_mixer` uses `allowed_mixers:enum[]`, history_min/max and positive q0_min/max (Bohr^-1) if Kerker is allowed; `select_eigensolver` uses `allowed_solvers:enum[]`; tolerance rules use the immutable floor and current-value ceilings; budget rules specify max_iterations; checkpoint rules specify allowed kinds/max_bytes. Empty parameter actions require an empty ranges object. No rule means deny autonomous application; the request may still seek supervised approval. Action-specific fields belonging to another action are invalid.

Policy file loading is driver-only and content-hashed in provenance. Restrict counts to host quotas and beta/history ranges to hard invariants; “allowed” cannot enable an unavailable algorithm, changed physics, tolerance weakening or convergence override. Counters/cooldown persist across mode switches and exact restart; a new connection cannot reset them. Track cooldown by completed SCF iterations since last applied action, not LLM message count. Count numerical and lifecycle applications in max_applied_total; read-only snapshots are separately rate-limited. A deny outcome never spends an application count, but request-rate quotas still apply.

## Acceptance, application and failure atomicity

Network receipt acknowledgement is not acceptance. Ordered lifecycle is received -> pending approval (if needed) -> validated -> `control_accepted` -> `control_applied`, or received/pending/validated -> `control_rejected`. Acceptance means semantic validation, approval and queue/journal capacity succeeded at a boundary, not merely that the bytes parsed. Between accepted and applied the mutation lock is retained; no intervening step/control can change the validated state.

Preconstruct/validate replacement mixer/solver state and reserve journal/event capacity first. Write durable intent and stage the new state; publish a durable commit record before making the replacement visible as the next runnable boundary. If persistence/preparation fails before commit, leave original state intact and emit `control_rejected` (possibly following accepted) with `io` or `prepare_failed`. A committed transaction recovered after a crash is applied exactly once from its commit record/checkpoint; an intent without commit is discarded. Never execute the next SCF step until committed state and audit agree. Catastrophic storage loss terminates the calculation rather than asserting durable audit. This is a required transaction protocol, not an existing engine guarantee.

`control_applied` includes before/after state_version and policy_revision and an audit reference. A checkpoint action's applied event means capture at that immutable boundary; later `checkpoint_written` means durable publication, while failure is an error event tied to the request. Mode/lifecycle-only actions change state_version but not numerical policy_revision. Rejections leave numerical state, versions, mixer history and RNG unchanged; event seq and rejected-count telemetry may advance. Idempotent retries return the original disposition and do not emit a second applied event. Requests rejected before admission by authentication/rate limits receive a transport error and bounded aggregate audit; admitted requests always receive individual accepted/rejected/applied dispositions.

## Agent absence and deterministic fallback

Initial default is `continue_current`: no agent connection, slow response, timeout, model service error or invalid proposal does not change numerical parameters. Continue until certificate, numerical failure or configured resource budget. Local diagnostics continue and no live LLM fallback is required. Agent outage is advisory metadata, never evidence of SCF failure.

Optional `local_rules_v1` is an operator-approved deterministic rule program selected at launch, evaluated every boundary regardless of agent connection status. It is not enabled by a wall-clock timeout: that would make trajectory depend on network timing. Its rules/ordering/cooldown are frozen and audited like any other controller. Optional `pause_at_boundary` is an explicit supervised-review policy triggered by a recorded diagnostic/control condition, not by an arbitrary missing heartbeat; record the trigger boundary for replay. Where both local and agent proposals exist, fixed order is operator cancellation, approved operator control, local rule, agent proposal; at most one numerical policy transaction per boundary, lower-priority numerical proposals are rejected as superseded/stale. Certificate/terminal completion precedes all controls. Record this arbitration order in schema/driver version.
