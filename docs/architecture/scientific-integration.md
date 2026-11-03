# Scientific integration contracts

Status: design only, not implemented. No Python, MCP, MLIP or coupling code is authorized by P11R. Adapter calls ultimately wrap `ScfProblem::new/with_kpoints`, the proposed session, `ScfOutcome` and `EnergyAudit`. [ADR 0017](../decisions/0017-coupling-commitments.md) remains authoritative for physical roles and energy ownership.

## Ownership boundary with P11I (D1/D2/D3)

P11I-COUPLING-API-SPEC owns the canonical high-level workflow/ASE/Python contract (`docs/architecture/workflow-api.md`), low-level coupling contract (`docs/architecture/coupling-api-lowlevel.md`) and MLIP/active-learning contract (`docs/architecture/active-learning-contract.md`). Those files are not present in this checkout when this addendum is written; see the P11I ticket on main `98f0b50`. Do not treat the absent files as implemented features. The first P11R draft's duplicate structure, prediction and field request schemas are withdrawn in favor of that single owner.

P11R supplies the runtime adapter boundary: validate a P11I scientific request into `ScfProblem::new/with_kpoints`, drive the one session under a declared mode, and return qualified `ScfOutcome`/`EnergyAudit` references with provenance. Modes never change P11I's physical interpretation. P11G implements optional Python/ASE binding; P11F implements MCP over the same driver; P11H presents CLI capabilities. None creates a second SCF/control implementation.

Integration invariants retained here, not competing schemas:

- ASE units convert only at I/O boundaries; unavailable forces/stress raise explicit unsupported errors, never zeros. Certified energy is distinct from partial unconverged output. P11I owns shape/layout/cache/observable details and P11G validates binding conventions.
- MLIP uncertainty must distinguish calibrated from uncalibrated/unavailable and retain calibration reference; no silent model mixing or zero-filled DFT force labels. Producer/model identity survives runtime caching, provenance and workflow edges. P11I owns the prediction/dataset contract.
- External numerical solvers use P11I's low-level contract, not agent prompt/event text. The runtime enforces resource bounds, immutable environment per inner calculation, and separate inner certificate versus outer status. No LLM or HTTP dependency enters the tight numerical exchange path.
- E1 records P11I request/response artifact digests; E2 references P11I typed artifact identities. They do not redefine field layouts, force terms, zero-copy lifetimes or active-learning data schemas.

## Provenance E1 proposal

`ProvenanceRecord`: `schema_version:string`, `provenance_id:string` (content digest), `calc_id:string`, `request_id:string`, `parent_ids:string[]`, `input_sha256:string`, `model_digest:string`, `environment_digest:string`, `source_commit:string`, `dirty_patch_sha256:string|null`, `lockfile_sha256:string`, `compiler:string`, `backend:string`, `execution:ExecutionIdentity`, `pseudopotential_sha256_by_species:map<string,string>`, `initialization:{kind:uniform/exact/warm,source_id:string|null}`, `policy_log_sha256:string`, `output_manifests:string[]`, `status:converged/exhausted/failed/cancelled`, `certificate_present:bool`, `created_utc:string`. This fills metadata missing from `ScfResult`; values are adapter-supplied facts, never inferred from energy. Exclude timestamps/host data from content identity; canonical schema ordering and round-trip float representation are versioned. Record exact command as data with secrets redacted, never evaluate it on import. UPF bytes are not embedded or sent to agents.

## Workflow DAG E2 proposal

`Workflow`: `schema_version:string`, `workflow_id:string`, `nodes:Node[]`, `edges:Edge[]`, `budget:{max_nodes:u64,max_total_iterations:u64,max_wall_seconds:u64,max_artifact_bytes:u64}`. `Node`: `id:string`, `kind:scf/inspect/export/mlip/external_solver`, `request_digest:string`, `input_artifact_ids:string[]`, `requested_observables:string[]`, `retry_limit:u64`, `requires_certified_inputs:bool`. `Edge`: `from:string`, `to:string`, `artifact_role:string`. Node status is `pending/ready/running/succeeded/failed/cancelled/blocked`; result references a ProvenanceRecord. Validate unique IDs, acyclicity, typed artifact compatibility and aggregate quotas before scheduling. Unknown node kinds are rejected; no shell or arbitrary-code node.

Execution wraps session creation and results; independent ready nodes may run only under host resource limits. Repeated request id plus equal canonical payload deduplicates; differing payload conflicts. Retries have new attempt ids and lineage, preserve failures, and cannot silently change model/policy. Inner SCF certification is separate from workflow success and outer coupled convergence. DAG caching needs matching full model/execution requirements; an old unconverged result never satisfies a certified dependency. No scheduler currently exists.

External numerical solver D3 uses the separate [contract](../research/external-solver-contracts.md). These adapter proposals refine ADR 0017's deferred runtime without authorizing field insertion or new physical models.
