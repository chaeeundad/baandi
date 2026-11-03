# Scientific workflow API

Status: **Designed**, P11I, 2026-11-03. Semantic schema v1.0, not executable Rust, Python or an ABI. No interface in this document is implemented. Existing building blocks are identified separately below. Extends [scientific coupling](scientific-coupling.md), [ADR 0016](../decisions/0016-esle-layering.md) and [ADR 0017](../decisions/0017-coupling-commitments.md); decision: [ADR 0021](../decisions/0027-coupling-api-split.md).

## Relationship to P11R

P11R documents are not merged into this checkout. Reviewed read-only in the sibling `P11R-RUNTIME-ARCH-DOCS` worktree: `agent-native-runtime.md`, `scientific-state-api.md`, `runtime-control.md`, `scientific-integration.md`, `external-solver-contracts.md` and ADR 0020. These are pending dependencies, not implemented APIs. Their paths below refer to the eventual `docs/architecture/` or `docs/research/` destinations. Reconcile the delivered versions before implementation.

Use P11R's single runtime, lifecycle, `ProvenanceRecord`, `ExecutionIdentity`, `ArrayRef`, workflow DAG and checkpoint envelope. This specification refines D1/D2/E1/E2, not a second scheduler or persistence format. The ESLE `ElectronicState` remains a checkpoint payload; save/load belongs to adapters. P11R's `step` means a full SCF iteration with mixing when applicable and its nested certificate; current `step_kpoints` alone is not that step.

The later [P11 common addendum](../../tickets/P11-COMMON.md) governs execution modes and stale control even where the reviewed P11R draft predates it: `Classic`, `Observed`, `AgentControlled` are policies over one implementation. Classic has no event production; it cannot become Observed after start. Observed/AgentControlled transitions occur at safe points and are audited. Observed never waits for subscribers. AgentControlled adds validated requests, deterministic fallback and approval policy; P11F initially offers Observe and Supervised behavior, not Autonomous. Numerical operations never wait for an LLM.

## Common types and validation

**Designed.** `Id` is an opaque ASCII string of 1–128 characters; `Digest` is lowercase SHA-256 hex (64 characters); counts are u64, physical scalars finite f64. `Version` is `major.minor`. Unknown fields, duplicates, unsupported versions and unknown enum variants fail before execution; migration must be explicit. New optional fields require a negotiated minor version; changed meanings require a major version. No ABI stability is promised.

`Availability<T> = Available{value:T, quality_id:Id} | Unsupported{reason:Reason} | NotComputed{reason:Reason}`. `Reason` includes `missing_capability`, `invalid_combination`, `not_started`, `unconverged`, `failed`, `cancelled`, `not_requested`, `missing_derivative`, `incomplete_ledger`. No absent scalar/vector/tensor is encoded as zero. Unsupported required observables fail preflight; optional ones return an availability tag. Nonfinite results are errors, not available values.

| Record | Required fields and semantics |
| --- | --- |
| `Structure` | `cell_columns_bohr:f64[3][3]`, `positions_fractional:f64[N][3]`, `species:Id[N]`, `atomic_numbers:u16[N]`, `atom_ids:Id[N]`, `periodic:bool[3]`, `origin_bohr:f64[3]`; unique stable atom IDs, lengths equal, finite nonsingular cell, explicit species mapping |
| `ModelContext` | `model_digest:Digest`, `environment_digest:Digest`, `structure:Structure`, `pseudo_sha256_by_species:map<Id,Digest>`, `xc:Id`, `nspin:u64`, `electron_count:f64`, `occupation_policy:Id`, `cutoff_ha:f64`, `grid_shape:u64[3]`, `basis_ordering_id:Id`, `k_fractional:f64[K][3]`, `k_weights:f64[K]`, `k_reduction_rule:Id`; canonical k order and positive normalized weights |
| `ProducerIdentity` | Tagged `DFT{code_version,build_digest,model_digest,environment_digest}` or `MLIP{implementation_version,architecture_id,weights_sha256,training_dataset_digest,training_recipe_digest}`; identity never inferred from a label or energy |
| `WorkflowRequest` | `schema_version:Version`, `request_id:Id`, `payload_sha256:Digest`, `parent_request_id:Id|null`, `model:ModelContext`, `required_observables:Observable[]`, `optional_observables:Observable[]`, `mode:ExecutionMode`, `policy_ref:Id`, `budget:{max_iterations:u64,max_wall_seconds:u64,max_artifact_bytes:u64}`, `initial_checkpoint_ref:Id|null` |
| `WorkflowResult` | `schema_version`, `request_id`, `calc_id:Id`, `model_digest`, `environment_digest`, `producer:ProducerIdentity`, `lifecycle:State`, `inner_status:converged/exhausted/failed/cancelled/not_evaluated`, `outer_status:not_evaluated/not_converged/converged`, `observables:map<Observable,Availability>`, `quality:QualityRecord`, `provenance_ref:Digest`, `state_id:Id|null`, `branch_id:Id|null`, `dataset_refs:DatasetRef[]` |
| `DatasetRef` | `dataset_id:Id`, `version:Id`, `manifest_sha256:Digest`, `record_ids:Id[]`, `store_ref:Id`; immutable content, no arbitrary executable URL |
| `QualityRecord` | `quality_id:Id`, `certificate:CertificateProjection|null`, `validation_refs:Digest[]`, `scope:Id`, `limitations:Reason[]`; certificate absence has a reason; validation references identify evidence, not an assertion of accuracy |
| `CertificateProjection` | `source_calc_id:Id`, `iteration:u64`, `density_residual:f64` (Bohr^-3/2), `occupied_eigenvalue_change_ha:f64`, `density_threshold:f64` (Bohr^-3/2), `eigenvalue_threshold_ha:f64`, `gap_gate:passed/not_applicable/failed`, `passed:bool`; projection of the engine certificate plus its applied criteria, never caller-set convergence |

`State` is P11R's Created/Ready/Running/Paused/Converged/Exhausted/Failed/Cancelled. Creation/setup errors need not fabricate a result. Outer success additionally needs its own certificate (low-level specification); an inner certificate cannot establish it. Failed or exhausted attempts retain partial diagnostics explicitly uncertified.

`Observable` includes `energy`, `free_energy`, `forces`, `stress`, `density`, `eigenvalues`. Energy is `{value_ha:f64, scope:electronic_partial/coupled_total, ledger_ref:Id, reference_id:Id}`. Free energy requires a named thermodynamic functional and temperature, not an alias for energy. Forces are Cartesian `f64[N][3]` in Ha/Bohr, `F=-dE/dR`, keyed by atom IDs and a ledger scope. Stress is Cartesian symmetric `f64[3][3]` in Ha/Bohr^3, tension-positive, defined by `dE=Omega sum_ab sigma_ab d epsilon_ab` for symmetric strain at the stated constraints; pressure is `-trace(sigma)/3`. Derivative outputs must carry the same energy scope and functional identity. Density/eigenvalues use explicit array references, never routine snapshots. Forces and stress remain **Deferred** until complete derivatives pass validation.

## Operations and capability negotiation

**Designed** semantic methods: `capabilities() -> CapabilitySet`, `validate(WorkflowRequest) -> ValidatedRequest|Error`, `submit(ValidatedRequest) -> CalcId`, `status(CalcId) -> StateSnapshot`, `result(CalcId) -> WorkflowResult|NotReady`, `export(CalcId, Observable) -> Availability<ArrayRef>`. Runtime controls reuse P11R; no alternate policy patch API.

`CapabilitySet` carries `schema_version`, `implementation_build:Digest`, `flags:map<Capability,bool>`, and `supported_combinations:Combination[]`. Capability names include energy, density, eigenvalues, forces, stress, external_fields, electrostatic_export, exact_restart, warm_start and step_control. Each `Combination` names supported XC, spin, occupation, boundary, external-term class, observable set and validation references. Flags are necessary, not sufficient: validate the entire requested combination. Design documents never set runtime flags true. A local Hamiltonian setter is not evidence for `external_fields=true`.

Same request ID and canonical payload returns the recorded attempt; different payload is a conflict. Retrying a failure uses a new attempt ID and parent link, never deletes the old attempt. Cache keys cover complete model/environment, settings, initial state, requested properties and execution requirements; positions alone are insufficient. Digest algorithms and canonical float/order rules must pass identity fixtures before content-addressed caching is enabled. Timestamps and hosts are excluded from scientific identity. Environment updates create a new request/model fingerprint, never a live policy change.

Control requests carry calc ID, observed state version, action/parameters, preconditions and expiry (iteration or wall time); accepted, rejected and applied are distinct. Map revision checks to P11R's policy revision and completed iteration; a revision alone does not identify an advancing calculation. Stale requests are rejected or explicitly revalidated. Only allow-listed compatible implemented policies change at complete boundaries; structure/XC/pseudo/electron count/k mesh/cutoff require new calculations. Queue/event guarantees remain runtime responsibilities, including bounded nonblocking delivery and no dropped critical records; no large arrays in events.

## Adapter plans

**Designed; implementation Deferred to P11G.** Python binding wraps validated types, results and eventual session; PyO3/maturin and optional ASE packaging require the ticket's dependency/ABI review. ASE-facing cell rows and Cartesian Angstrom positions convert once to internal columns/fractional Bohr coordinates. Use dft-units, not copied constants. Certified energy converts Ha to eV; future forces convert Ha/Bohr to eV/Angstrom and stress to eV/Angstrom^3 with an explicitly tested tensor/Voigt mapping. Unsupported forces/stress raise a clear NotImplemented error. Unconverged results raise a convergence error with retained diagnostic references. Missing free-energy semantics also raise. Cache invalidates on all model changes. Python buffers need the lifetime rules of the low-level API; releasing the interpreter lock cannot invalidate an owner.

**Designed; implementation Deferred to P11F.** Research agents use the single separate MCP adapter over the same runtime. P9A consumes it. MCP passes small validated requests, snapshots and artifact references; no Python/LLM/HTTP dependency enters the engine, no arbitrary shell/file execution. Agent proposals cannot set certificates or fill missing labels. MLIP and active-learning clients consume the [active-learning contract](active-learning-contract.md), with no mandatory provider.

## Gap table against current code

Inspected checkout `98f0b50` on 2026-11-03; these are source observations, not new numerical verification.

| Interface or building block | Status | Current evidence and missing work |
| --- | --- | --- |
| SCF run/outcome and certificates | Implemented (internal library building block) | `crates/dft-engine/src/scf.rs`: `ScfProblem::run/run_outcome`, `ScfOutcome`, `FixedPointCertificate`; no workflow envelope |
| Energy, density, per-k eigenpairs | Implemented (internal result) | `ScfResult`, `EnergyAudit` in `scf.rs`/`energy.rs`; identity and dataset metadata supplied by future adapter |
| LDA/PBE, insulating k-point model | Implemented (library selections) | `ScfOptions`, `OccupationPolicy`, `dft-physics` XC; support is not blanket validation for all inputs |
| Density binary export | Implemented (CLI artifact) | `crates/dft-cli/src/run.rs`: density payload/manifest; not field import or restart |
| Versioned workflow/provenance/capability envelopes | Designed | This document and P11R; no runtime implementation |
| Full step session, events/control, checkpoints | Designed | P11R/P11A–D; existing `step_kpoints` is only a constituent solve |
| Python/ASE and MCP | Designed | P11G/P11F implementation pending; distinct adapters |
| Active learning and dataset validation | Designed | Contract below; no trainer or deployment service |
| Forces/stress, coupled total/outer certificate | Deferred | No derivative implementation or validated coupling ledger |
| External fields/electrostatic export/operator lease | Designed | Internal numerics exist; public external contracts and transport absent |
| Charged/open/slab boundary models, dynamic embedding | Deferred | No negotiated or validated implementation |
