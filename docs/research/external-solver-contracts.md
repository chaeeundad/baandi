# External solver runtime boundary

Status: design only. P11I-COUPLING-API-SPEC owns the scientific low-level API, not P11R. The first draft's duplicate CoupledRequest/FieldTerm/CoupledResponse schemas are withdrawn. The canonical specification is P11I's `docs/architecture/coupling-api-lowlevel.md`; its high-level companion is `workflow-api.md`, and efficiency/test rationale belongs to `docs/research/external-coupling-efficiency.md`. These files are pending in this checkout; P11I's ticket was read on main `98f0b50`.

“DFT-CES-like” identifies a class of external environment/electronic solver exchange, not a claim of ABI/package compatibility. [ADR 0017](../decisions/0017-coupling-commitments.md) remains the physical-contract authority. Existing seams are `ScfProblem::step_kpoints` potential assembly, `EnergyAudit`, `ScfOutcome` and `dft-solvers::LinearOperator::apply`; a seam is not an implemented external-field or apply-only export API.

## P11R-owned lifecycle handoff

| Runtime operation | Existing type/function wrapped | Scientific contract owner |
| --- | --- | --- |
| Admit validated scientific request | `ScfProblem::new/with_kpoints` plus future external-field validation | P11I request roles, units, basis, boundaries and capabilities |
| Run an inner calculation | Proposed shared session around `run_impl`, one declared execution mode | P11R steps, event/control safety, budgets and certificate |
| Return qualified inner result | `ScfOutcome`, `EnergyAudit`, checkpoint/state artifact references | P11I observables and energy/force accounting; P11R provenance/status |
| Start next environment iteration | New calculation/request with parent provenance | P11I environment identity and outer convergence; never a live physical-model patch |
| Expose future operator evaluation | Existing `LinearOperator::apply` as an internal seam only | P11I ownership/lifetime/transport and mathematical callback validity; not agent-facing arbitrary callback |

The outer numerical solver does not wait for agent approval at every matrix-vector call or inner SCF iteration. A supervisor can propose bounded runtime interventions asynchronously; only the control path may approve/apply them at SCF boundaries. Immutable environment and numerical ownership rules are not bypassed by AgentControlled mode. A changed environment is a new P11I request with lineage, even when the driver reuses compatible allocations or a validated warm start.

Events carry small status/evidence summaries and artifact IDs only. Dense field/operator exchange follows P11I's array transport/ownership contract, not JSON-lines telemetry or LLM messages. P11R neither defines a second array ABI nor chooses shared-memory/IPC implementations. Unsupported fields/observables fail capability negotiation before numerical execution; no force/stress or coupled-energy claim follows merely from runtime availability.

## Integration gates

P11I defines sign, gauge, boundary, energy/force and operator consistency tests before coupling claims. P11R adds runtime tests: separate inner/outer status; exact request identity and conflict handling; cancellation at a complete boundary; no agent dependency in a numerical callback; immutable environment during each inner solve; correct provenance across warm starts; and unchanged zero-field SCF trajectory. P11J measures control/event overhead separately from the transport/coupling costs owned by P11I. P11G/P11F/P11H consume capabilities from that same contract rather than duplicating validation.

No external-solver implementation, coupling experiment, transport benchmark or new numerical result is claimed by this document.
