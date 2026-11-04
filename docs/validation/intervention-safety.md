# Intervention safety

Status: threat model and test plan, with the executed P11E rejection/replay slice below. Enforcement wraps safe boundaries of `ScfProblem::run_impl`, validation of `ScfOptions`, `DensityMixer` reset/replacement, and read-only `ScfOutcome`; the agent never receives direct mutable numerical state. See [policy/audit schema](../architecture/runtime-control.md).

| Threat | Required mitigation | Acceptance test |
| --- | --- | --- |
| Unsafe policy or fake convergence | Strict typed allowlist, finite bounded numbers, no tolerance relaxation, immutable physical model, existing certificate sole authority | NaN/infinity, beta<=0 or >1, history=0, tolerance<=0, max_iterations<done, forged convergence and model changes all reject; compare state/mixer hashes before/after |
| Intervention spam/races | One mutation owner, exact expected revision/boundary, one pending pause/cancel, bounded policy queue | Concurrent patches at same revision: at most one accepted; stale one audited with no state change |
| Replay/conflicting requests | Authenticated actor, idempotency by id plus canonical digest, monotonic sequence and chained audit | Equal replay returns original decision; changed payload same id conflicts; restore does not reuse sequence |
| Resource exhaustion | Preflight memory/work quotas, bounded queues/history/artifacts/events, cooperative cancellation | Huge history/array dimensions, overflow, deep JSON and endless retries reject before allocation; budget exhaustion never certifies |
| Prompt injection through event text | Enum/typed evidence only for automatic actions; untrusted prose displayed as data, never authority | Malicious filenames/error messages cannot become shell commands, tool requests or policy fields |
| Checkpoint tampering/path traversal | Store-confined relative references, length/hash/schema/model validation and trusted audit head | Truncation, symlink escape, ../ paths, altered hash and incompatible model rejected |
| Secrets/UPF exfiltration | Local artifact IDs, authenticated read scopes, no UPF bytes/events, redaction and no arbitrary path read | Agent read cannot escape calc store or obtain credentials/pseudopotential files |
| Audit failure | Durable intent before transactional commit; stop accepting changes if persistence fails | Inject disk-full/write failure: no unaudited mutation; restart resolves incomplete intent as unapplied |

Initial adapter security profile: maximum 64 KiB per control command and 4 KiB per event, nesting depth 8, 16 queued control requests, 10 policy requests/minute/actor, 32 accepted interventions/calculation, 64 mixer vectors, 10000 iterations and 1 GiB artifact export/calculation. Authenticate quotas across connections; reconnect does not reset them. Host sets wall/memory/operator budgets at launch and may only tighten them during execution. Limits are proposed deployment safeguards, not tuned numerical parameters. Quota rejection audit is aggregated after the first rejection per rate-limit interval to avoid log amplification; persist rejection count and interval, without retaining unbounded attacker payloads. Requests accepted into the control queue always receive an individual decision.

Check a hard resource budget before starting the next iteration; in-flight cancellation waits for the safe point, so a wall deadline can be exceeded by one full step. OS-level worker termination is an adapter fallback after a configured grace interval; it cannot emit convergence or an exact checkpoint of partial work. Preserve the last complete checkpoint. Warn in capabilities that per-step memory/time bounds require profiling; no false immediate-cancellation promise.

Trusted boundary: authenticated host driver and immutable model store. Untrusted: agents, imported manifests, external solver messages, diagnostics prose and transports. A compromised trusted host is outside the arithmetic safety guarantee; hashes detect corruption, not host compromise. Agent-supervised execution is reproducible only through recorded decisions. Rejecting malicious control must leave numerical state bit-identical; tests also ensure diagnostics never alter the certificate, positivity or charge gates. Future MCP implementation P11F consumes these checks rather than implementing a parallel validator.

## Three-mode and asynchronous-control threat extensions

| Threat | Enforcement | Required negative test |
| --- | --- | --- |
| Old agent advice applied after recovery | State-version comparison, fresh typed predicates, unchanged policy revision, explicit expiry and bounded revalidation | Oscillation alert at iteration 12 but cleared by 18 rejects without mixer reset |
| Approval laundering | Host-authenticated approval bound to request digest and permitted revalidation; agent cannot self-approve | Change beta/payload after approval, replay another actor's approval or use expired approval: reject |
| Unlimited “autonomous” access | Hashed operator-installed allow-list with per-action ranges/counts/cooldown plus host caps; default Supervised | Omitted action, unsupported solver, fourth application of a three-action grant or early cooldown: reject |
| Mode used to evade audit | Switch table, operator authority, no post-start Classic conversion, revoke unapplied proposals on downgrade | AgentControlled->Classic or rewriting checkpoint mode rejects; downgrade preserves counters/history |
| Slow client stalls SCF | Independent bounded subscriber queues, critical journal replay, no subscriber acknowledgements on fast path | Permanent subscriber stall cannot suspend an otherwise healthy solve; reconnect retrieves critical events |
| Critical-event loss hidden as routine coalescing | Reserved critical credits; named never-drop classes; stop before starting unjournalable work | Fill critical storage, observe explicit terminal telemetry error rather than missing certificate or silent continuation |
| Agent outage changes numerical outcome | Default continue_current; local fallback runs by fixed boundary rule, never heartbeat timing | Same applied-control log with different reply delays gives identical scientific results |

Control requests admitted before a downgrade receive an explicit disposition; proposals received afterward reject `mode_disallows_control`. Observed read-only access is not authority to change settings or request large arbitrary allocations. Checkpoint generation and approved mode switches remain host-authorized operations. Revoking an agent connection does not cancel already committed transactions or erase their audit.

Retain only typed enum/scalar evidence in machine-directed events, with full windows in bounded artifacts. An agent's explanation may contain instructions or misleading hypotheses; treat it as untrusted display text. Never turn it into a predicate, approval, fallback action or new allowed field. Escape displayed HTML/terminal controls in P11H, and preserve raw text only in a separately permissioned artifact if needed.

## P11E executed slice (2026-11-04)

`scripts/demo_self_healing.sh` executes real synthetic SCF with the merged runtime
validator. Each recorded-agent repeat submits an operator-approved beta=1.1
request at boundary zero. Both reject as `Policy(InvalidBeta)`; the boundary
snapshot and initial density bits are unchanged. A subsequent bounded budget
extension is separately accepted and applied (revision 0 to 1), preserving mixer
history. The final result matches the no-intervention reference bit for bit,
including all array elements and normalized history/certificates. This tests
rejection's lack of numerical side effects through the subsequent calculation.

The integration harness additionally rejects unknown decision fields, missing
observation versions, stale observations and future observations. The replay
file cannot set authority. The host installs one budget extension up to 200;
physical model and convergence tolerances are never changed. Exact audit
outcome order is asserted, and every audit hash-chain link is independently
recomputed from the stored JSON. Evidence is under
`tests/reference/runtime/generated/agent-{0,1}/audit.jsonl`.

Scope is the trusted offline demo driver, not a deployed security boundary.
The journal syncs and atomically publishes each batch before application but
has no crash recovery loader; it is not the general checkpoint format. The
broader quota, remote authentication, prompt-injection, resource-failure and
stale-diagnostic threat matrix above is not certified by this demo.
