# Self-healing SCF demonstration protocol

Status: the initial budget-recovery protocol is executed by P11E (results below). Instability recovery and performance claims remain unvalidated. P11E owns the harness and stored results. Wraps the fixture constructors in `crates/dft-engine/tests/scf.rs`, `ScfProblem::run_outcome` and the proposed session. Use [diagnostics](../research/self-healing-scf.md), [safety](../validation/intervention-safety.md) and [equivalence](../validation/runtime-equivalence.md).

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


## Executed P11E demo

Reproduce offline from the repository root:

```sh
scripts/demo_self_healing.sh
python3 -I scripts/make_readme_figures.py
```

Requires the pinned Rust toolchain and the repository's vendored dependencies;
no real UPF, external server, network, or LLM is used. The script runs the ignored
release integration test in `crates/dft-cli/tests/self_healing.rs`. The adapter
`src/heal.rs` is compiled by that test using a path module, deliberately avoiding
changes to the independently owned CLI command dispatcher. There is no new
`dft heal` subcommand. Override `DFT_HEAL_OUTPUT` to write an independent bundle;
the default is `tests/reference/runtime/generated/`. A failure removes the prior
aggregate results file before running, so stale success cannot masquerade as a
successful rerun.

The frozen model is the synthetic fixture specified above, with expanded
`ScfOptions` in each summary, a uniform-density start, seed 2026, one requested
and effective worker, serial faer/FFT, and unchanged engine defaults. The
reference and all three arms use AgentControlled/RetainLast and fresh starts;
only the initial iteration budget differs. A comparison of all initial density
hashes ensures identical starts. The reference has 1000 steps available; the
three arms initially have one. The model is synthetic, not real silicon.

The initialized event supplies the initial observation. `budget_preflight_v1`
recognizes that fewer than two iterations cannot pass a certificate because the
first energy change is absent. The rules arm constructs the budget request;
the agent arm reads the checked-in
`tests/reference/runtime/agent-decision.json`. The host installs a bounded
grant: at most one budget extension, maximum 200, no stale revalidation.
The recorded agent cannot provide its own authority. Both requests execute at
boundary zero, then the same session continues without rebuilding its model.
This is an agent-supervised replay under a host-installed grant, not a live LLM
or a claim of asynchronous human review.

The implementation uses the merged engine schema (`schema_version: 1`, runtime
state `Initialized`) rather than the earlier draft's `"1.0"` and `Ready` spelling.
The compact decision file maps to a complete ControlRequest with fixed typed
state/revision/mode preconditions, exact observed version, expiry at iteration
zero and stale rejection. Unknown fields (including authority) and missing
fields fail closed. Actual iteration events feed P11D diagnostics between steps;
these findings are advisory and never set convergence. Budget repair is a
separate preflight rule, not a fabricated stagnation diagnostic.

| Arm | Terminal status | Ordinary iterations | Certificate solves | Operator applications, including certificates | Certified energy (Ha) | Absolute reference difference (Ha) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Baseline | Exhausted | 1 | 0 | 43 | unavailable | unavailable |
| Deterministic rules | Converged | 30 | 3 | 951 | -4.684561184277879 | 0 |
| Recorded agent | Converged | 30 | 3 | 951 | -4.684561184277879 | 0 |
| Fresh reference | Converged | 30 | 3 | 951 | -4.684561184277879 | 0 |

Both repeats of every arm pass. All final density and complex orbital elements
are compared through SHA-256 of their IEEE bits, along with a hash of the complete
normalized result (including history and certificates). Rules and agent also
match the reference exactly. Certificate density residual is
4.551166513318054e-9 Bohr^-3/2; occupied eigenvalue change is
2.7803310209684184e-8 Ha; final relative charge error is
4.440892098500626e-16. The baseline's retained energy is
-4.2676391604117985 Ha and is explicitly **not certified**. Its certified energy
and reference-error fields are null.

Evidence files in `tests/reference/runtime/generated/`:

- `results.json`: aggregate assertions, units, source SHA-256 values, input/options,
  per-arm status and metrics. `build.txt`: compiler, Cargo, host architecture,
  parent revision and working-tree status. Local uncommitted source hashes bind
  the implementation; the parent commit alone is not its source identity.
- `<arm>-<repeat>/events.jsonl`: ordered complete engine events and typed scalar
  measurements. Every event is drained at each boundary; any dropped-event
  marker fails this demo. `diagnostics.json`: preflight and actual P11D findings.
- `audit.jsonl`: accepted/applied/rejected records, complete debug request and
  authority, typed observation versions, before/after policies, actor and time;
  an ordered SHA-256 chain detects edits. The beta=1.1 negative control has host
  approval but still rejects as `Policy(InvalidBeta)`. Rules have one applied
  request; agent has one rejection and one applied request per repeat.
- `arrays.json`: hash, dimensions and explicit little-endian element ordering;
  `normalized-result.txt`: complete same-build debug result with wall timers
  zeroed. These are evidence, not a portable restart format or volumetric export.

The audit adapter publishes a complete batch by file sync and atomic rename
before runtime mutation, propagating persistence errors. It does not implement
crash recovery, hostile-directory confinement or a cross-build audit parser.
Wall seconds measure runtime initialization, control, stepping and observation;
they exclude model construction and final evidence export. Observer overhead is
null: no paired performance study was performed. There is no speedup claim.

Limits: one synthetic Gamma system, small band count, initial budget repair
only; no new numerical defaults, physical instability intervention, delayed
oscillation replay or global-ground-state proof. The numerical-instability
extension above remains a future experiment, not completed acceptance evidence.
P11J retains ownership of mode overhead and the broader equivalence matrix.
