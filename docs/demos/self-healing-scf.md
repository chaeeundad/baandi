# Self-healing SCF demonstration protocol

Status: executable protocol design, not a completed demo or a measured speedup. Future P11E owns the harness and stored results. Wraps the fixture constructors in `crates/dft-engine/tests/scf.rs`, `ScfProblem::run_outcome` and the proposed session. Use [diagnostics](../research/self-healing-scf.md), [safety](../validation/intervention-safety.md) and [equivalence](../validation/runtime-equivalence.md).

## Reproducible baseline and honest scope

Public synthetic fixture: `common::potential(4.0,1.0,true,801)`, `common::cell()` and `common::atoms()` from `crates/dft-engine/tests/common/mod.rs`, species Xx, Gamma `ScfProblem::new`, `ScfOptions::default()` except `max_iterations=1`, threads=1 and seed=2026. Pin the generator and default values by source hash and expanded configuration. The cell helper uses a=10.261212856717933 Bohr and primitive rows (0,a/2,a/2), (a/2,0,a/2), (a/2,a/2,0); atoms have fractional positions (0,0,0) and (1/4,1/4,1/4). Cutoff=5 Ha, Pulay beta=0.3/history=6, density tolerance=1e-7 Bohr^-3/2, energy/eigen tolerances=1e-8 Ha. The existing `divergence_and_iteration_budget` test asserts one-iteration `NotConverged`. Since first-iteration energy change is absent, no first-iteration certificate can pass. This is a guaranteed budget-failure case, not evidence of numerical stagnation.

Existing beta=100 test is a separate negative-control divergence fixture. Reproduce it only through the legacy test command; runtime policy must reject beta>1 before execution. It is not an allowed runtime starting policy or proof of successful live recovery. The historical P7PC Si8 eigensolver stagnation is described in [its ticket](../../tickets/P7PC-BUFFER-BANDS.md). It requires a private UPF/input and may disappear after P7PC merges; do not advertise it as a portable current failure or try to fix Davidson stagnation by changing density mixing.

## Three arms

All arms start from the same immutable synthetic model, policy, seed, one-iteration budget and empty history. Observation and intervention happen at initial Ready boundary `(iteration_done=0,revision=0)`; only the budget differs after the intervention.

| Arm | Deterministic protocol | Expected claim |
| --- | --- | --- |
| Baseline | No intervention; retain_last mode, max_iterations=1 | Exhausted, converged=false |
| Rules | Rule `budget_preflight_v1`: if max_iterations<2, submit one budget extension to 200 at boundary 0 | Tests safe budget repair; successful convergence must still be measured |
| Agent-supervised | Replay the scripted decision below at boundary 0 without a live LLM | Same accepted trajectory as rules arm; recorded actor differs |

Decision file content is one PolicyPatch: `schema_version="1.0"`, matching `calc_id`, `request_id="demo-budget-1"`, `expected_revision=0`, `apply_after_iteration=0`, authenticated `actor_id="demo-agent"`, `reason_code="budget"`, `changes={max_iterations:200}`, `reset_history=false`. Record the complete request and audit result. A live agent is an optional fourth arm and must export replayable decisions; no LLM dependency is required for the three-arm comparison.

Fully converged reference: run the same model with unchanged tolerances, no interventions, Pulay beta=0.3/history=6, max_iterations=1000, same build/thread count; require existing certificate and gap gates. If it fails, record failed preflight and stop the demonstration; do not loosen tolerances or claim recovery. Baseline and reference are fresh runs, not checkpoint imports. Rules and recorded-agent outcomes must agree with reference energy within 1e-10 Ha and with each other bitwise after timing/identity normalization. The reference certificate is necessary but this agreement alone does not prove a global ground state.

## Numerical-instability extension (not yet validated)

To earn a broader “self-healing instability” claim, P11E must additionally find and freeze a valid-policy case: start with the same fixture and Linear beta=1, max_iterations=200, and test the six-sample stagnation/oscillation predicates. Compare no intervention with the documented beta-halving/reset rules and an identical recorded-agent action sequence, using a fully converged reference. If no failure/stagnation occurs, report that and use P10D's independently validated elongated/sloshing fixtures; do not select thresholds after seeing the recovery. Publish the frozen model, failure history and rule version before comparing arms. Until then the guaranteed demo above proves control/replay/budget handling only, not numerical healing or generality.

## Evidence bundle and metrics

For each arm retain expanded input/hash, generator revision, build/backend/thread identity, start-state hash, event JSONL, policy/audit JSONL, outcome/certificate and final array manifest. Report status, ordinary SCF iterations, certificate solves, ordinary plus certificate operator applications, final E (Ha), abs(E-E_reference) (Ha), density residual (Bohr^-3/2), occupied eigenvalue change (Ha), charge errors (relative), interventions accepted/rejected and reasons, wall seconds and observer overhead. Missing final certified energy stays null, not an apparent zero error. No fabricated numbers or QE-reference label. Run each deterministic arm twice; replay must preserve numerical outputs. A release harness command and exact output paths must be supplied by P11E before claiming this protocol executed.

Existing baseline-only reproduction command: `cargo test --release -p dft-engine --test scf divergence_and_iteration_budget -- --ignored --exact --nocapture`. This existing test is ignored in the ordinary debug suite; a debug invocation without `--ignored` does not reproduce its assertions. The release command tests failure behavior, not the future three-arm runtime, and was not executed in this design ticket. Tests over about 10 s debug must be ignored with a documented release command. Extend evaluation across P10G system variety before recommending recovery defaults; this two-atom synthetic demonstration is deliberately narrow.

## Three-mode/control-request extension

Run the three controller arms in AgentControlled with the same frozen execution config and numerical initial policy so mode choice does not confound recovery comparison. Baseline simply sends no control. Fully converged reference may run Classic; P11J separately proves mode equivalence. Rule/recorded-agent arms need an operator-installed grant permitting exactly one `extend_budget` to 200 at boundary 0 (the generic beta allow-list example does not authorize this).

The earlier decision-file PolicyPatch is the low-level intended change, not the public control wire object. Wrap it as ControlRequest v1.0: matching calc/request ids, `observed_state_version=1` after initialization, `observed_policy_revision=0`, `action=extend_budget`, `parameters={max_iterations:200}`, preconditions `state_is(Ready)`, `policy_revision_is(0)`, `mode_is(AgentControlled)`, `expiry={max_iteration:0,deadline_utc:null}`, `stale_policy=reject`, `reason_code=budget`, `evidence_id=null`. The deterministic harness submits it synchronously at initial Ready so the intentionally short expiry is met; this is scripted control, not a claim that a live LLM can answer before the first step. Emit accepted and applied records separately, with version/revision change, and keep the baseline grant unused.

Add a separate stale-agent scenario to the instability extension: freeze an evidence trace where oscillation becomes clear without intervention; replay a delayed proposal tied to the earlier version with `diagnostic_active` precondition. Required outcome is rejection after revalidation, not unnecessary beta reduction. A second scripted trace in which oscillation persists must exercise approved bounded application. If no physical fixture supplies both traces, use labelled control-state unit fixtures and do not present them as measured SCF healing. Capture event coalescing/loss summaries and request latency alongside the original metrics.

The live-agent optional arm never blocks a step. If it misses the expiry or the solve ends first, report no applied intervention and preserve that outcome; do not rerun until a visually convincing recovery is obtained. Explanations remain labelled inference, and the fully converged-reference/certificate gates remain unchanged.
