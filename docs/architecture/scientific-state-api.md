# Scientific state and persistence API

Status: P11A state/events and P11C checkpoint are merged (2026-11-05); broader semantic schemas below remain design contracts, not a claim that every proposed field is exported. Pure projections wrap `IterationRecord`, `FixedPointCertificate`, `ScfResult` and `ScfOutcome` in [scf.rs](../../crates/dft-engine/src/scf.rs). Serialization belongs to dft-cli. [Control](runtime-control.md) defines lifecycle and revisions.

## Common wire rules

Schema v1.0: objects reject unknown fields, duplicate keys and wrong types; integers below are unsigned 64-bit unless stated otherwise; f64 values must be finite. IDs are opaque ASCII strings of 1–128 characters, digests are lowercase 64-character SHA-256 hex. UTC timestamps are RFC3339 strings, assigned by the adapter and excluded from scientific identity. Unsupported major/minor versions are rejected unless a named explicit migration exists. Null means unavailable with a companion reason, never zero. Numerical units are Hartree/Bohr, cell vectors are columns, fractional k coordinates are reciprocal-basis coordinates.

## Read-only StateSnapshot v1.0

| Field | Type and unit | Source |
| --- | --- | --- |
| `schema_version`, `calc_id`, `model_digest`, `environment_digest` | string; environment digest identifies even the no-environment case | Validated immutable model envelope |
| `seq`, `iteration_done`, `state_version`, `policy_revision` | u64 | Driver sequence, completed `IterationRecord` count, accepted revisions |
| `state` | lifecycle enum from control document | Session |
| `mode` | Classic/Observed/AgentControlled | Immutable creation config or approved mode transition |
| `energy_ha`, `energy_change_ha` | f64 or null, Ha | `energy.total_hartree()`, `energy_change` |
| `density_residual` | f64 or null, Bohr^-3/2 | `density_distance`, sqrt(integral(delta n)^2)/N |
| `max_eigen_residual_ha`, `requested_eigen_tolerance_ha` | f64 or null, Ha | Last record |
| `charge_error_relative`, `input_charge_error_relative`, `mixing_fraction` | f64 or null, dimensionless | Last record |
| `min_density_bohr3_inverse` | f64 or null, electron number/Bohr^3 | Last record |
| `operator_applications`, `certificate_operator_applications` | u64 cumulative vector applications, separate sums | Records and nested certificates |
| `davidson_iterations` | u64 cumulative inner iterations, includes certificates | Records and certificates |
| `converged`, `mixing_applied`, `mixer_reset` | bool | Certificate outcome and boundary report |
| `certificate` | null or `{density_residual:f64, occupied_eigenvalue_change_ha:f64, passed:bool}` | Last fixed-point certificate; `passed` also requires gap gate |
| `unavailable` | map of nullable field name to enum `not_started/not_computed/nonfinite/failed` | Explicit absence |

No density, orbital arrays or unbounded history in snapshots. Large k sets are summarized, not copied. Optional array access is an explicit artifact request at a safe point, with byte quota and manifest reference; it is not an event field. Snapshot does not mutate or trigger an eigensolve.

## JSON-lines event v1.0 (three-mode revision)

Events exist only in Observed and AgentControlled. One UTF-8 JSON object per line, maximum 4 KiB after serialization; no array-valued field anywhere in the event, including diagnostics. The draft's nested snapshot/audit payloads and vector evidence are replaced by scalar summaries plus opaque artifact IDs. Full snapshots, audit records and diagnostic windows are separate on-demand store objects. Building/sending them is outside the fast path and subject to byte quotas.

| Required field | Type | Meaning |
| --- | --- | --- |
| `schema_version` | literal string `1.0` | Wire schema, not solver version |
| `calc_id` | bounded opaque string | Calculation identity |
| `seq` | u64 | Strictly increasing produced-event sequence, starting at 1; gaps indicate declared dropping/coalescing |
| `scf_iteration` | u64 | Last completed iteration, 0 before first completion; failure mid-step retains last completed index |
| `state_version` | u64 | Monotonic published boundary-state version, distinct from event seq and policy revision |
| `policy_revision` | u64 | Successful numerical-policy commits only |
| `type` | closed enum below | Event interpretation |
| `timestamp_utc` | RFC3339 string | Adapter wall clock; may move backward, never used for event ordering |
| `scalars` | tagged object of named finite f64/u64/bool/null fields | Values with units defined below; no arbitrary map or serialized prose |
| `details` | tagged object of bounded enum/string/null fields | IDs and reason codes; no arrays or executable text |

`state_version` starts 0 (Created), advances once after initialization, each completed public SCF step, failed step, applied control transaction or completed lifecycle/mode transition. A terminal step and its terminal state are one commit/version; its multiple emitted events share that version. Internal Running entry and pending requests do not mutate the published boundary state. Snapshot reads, rejected/accepted-but-unapplied proposals and subscriber activity do not advance it. `state_changed` events for internal Running notification may carry the preceding committed version and are not usable as a mutable safe point. Overflow of any counter terminates before reuse. Checkpoints persist versions and last assigned sequence.

| `type` | `scalars` | `details` | Delivery class |
| --- | --- | --- | --- |
| `iteration_completed` | `energy_ha`, `energy_change_ha` nullable; `density_residual` (Bohr^-3/2), `max_eigen_residual_ha`, `charge_error_relative`, `operator_applications:u64`, `certificate_operator_applications:u64`, `mixing_applied:bool` | `snapshot_id:string` | Routine, coalescible |
| `diagnostic` | Rule-specific scalar summary (ratios, counts, thresholds) | `rule_id`, `status`, `evidence_id`, `inference_code` | Routine unless error or diagnostic status changes |
| `certificate` | `passed:bool`, `density_residual`, `occupied_eigenvalue_change_ha`, `operator_applications:u64` | `certificate_id:string`, `gap_status:passed/failed/not_applicable` | Never drop, every attempted certificate |
| `state_changed` | empty object | `from`, `to`, `reason` (lifecycle enums) | Never drop |
| `mode_changed` | empty object | `from`, `to` (mode enums), `request_id` | Never drop |
| `control_accepted` | `observed_version:u64`, `validated_version:u64` | `request_id`, `action`, `audit_id`, `approval_id` nullable | Never drop |
| `control_rejected` | `observed_version:u64|null`, `validated_version:u64` | `request_id`, `action` nullable, `reason`, `audit_id` | Never drop for admitted requests |
| `control_applied` | `before_version:u64`, `after_version:u64`, `before_revision:u64`, `after_revision:u64` | `request_id`, `action`, `audit_id` | Never drop |
| `checkpoint_written` | `byte_length:u64` | `artifact_id`, `kind`, `manifest_sha256`, `request_id` nullable | Never drop |
| `failed` | `attempted_iteration:u64`, `record_available:bool` | `code`, `stage`, `evidence_id` nullable | Never drop |
| `events_dropped` | `iteration_events:u64`, `diagnostic_events:u64`, `first_missing_seq:u64`, `last_missing_seq:u64`, `coalesced_count:u64` | `scope:producer/subscriber`, `subscriber_id:string|null` | Never drop; counters describe gaps, not missing critical events |

Reason/error enums extend the earlier draft: `initialized/step_started/step_finished/certificate_passed/budget/pause/resume/cancel/failure` and `invalid_input/nonfinite/diverged/eigensolver/noninsulating/io/quota/incompatible/internal/telemetry_exhausted`. Control reason enums are specified in [control](runtime-control.md). Null scalar values require a typed `unavailable_reason` in details (`not_started/not_computed/nonfinite/failed`); invalid numbers are never written as JSON NaN/Infinity. The full physical record remains in numerical history even if its routine event is coalesced.

## Delivery, capacity and slow consumers

Initial proposed per-calculation limits: 256 routine slots and 128 critical slots, maximum 4 KiB wire payload per slot, plus 2 reserved emergency slots. Thus event payload capacity is at most 1,581,056 bytes; allocator/object overhead and journals are measured separately, not claimed to fit that number. Driver queues use preallocated storage; no event allocation is added inside numerical loops. The numerical producer hands off a bounded typed boundary record non-blockingly. Formatting, wall timestamps and JSON-lines I/O run in the driver/sink worker.

1. Routine `iteration_completed`: coalesce older unsent iterations to the newest complete scalar summary; retain latest energy/residual but do not sum them. Cumulative counters remain cumulative. Increment dropped/coalesced counters for each removed event. If no slot is available, drop this routine event and increment counters; never delay SCF.
2. Routine diagnostics: coalesce repeated unchanged `(rule_id,status)` to latest evidence reference. First activation, clearing, and any error-severity diagnostic are critical, so consumers can track condition transitions. Full six-point evidence is in the bounded diagnostic store, not an event array.
3. Never-drop list: all certificate, lifecycle/mode, accepted/rejected/applied control, checkpoint completion, failure/error, diagnostic-transition and `events_dropped` events. Critical events are retained in an append-only journal and delivered/replayed by `(calc_id,seq)`; clients acknowledge a cursor. No end-to-end delivery promise is made to a disconnected client, but critical history remains queryable for the retention period.
4. Flush one `events_dropped` aggregate before the next critical publication or terminal event, and at least every 64 completed iterations if any loss occurred. Report first/last missing seq as a bounding interval, not a claim that every seq in it was dropped. Receiver sequences are ordered; queue merging never places a lower seq after a higher one. Snapshot response gives the latest state version independently of stream loss.
5. Subscriber queues are independent and bounded (64 slots each in the example profile). Slow consumers lose/coalesce routine items with subscriber-scoped loss reporting. If their critical queue fills, disconnect them with a resume cursor; recover from the journal. Never block other subscribers or SCF for their acknowledgement. Fan-out is capped at 8 subscribers/calculation; others are rejected before allocation.
6. Local subscribers receive borrowed scalar records after a step, never density/orbital references or arbitrary callbacks inside it. Only trusted built-in bounded diagnostics may execute synchronously in the control path. Plugins, UI, sockets and agents always use async queues. Slow/failing optional local observers are detached; their failure cannot change the numerical result. An observer that requires arbitrary blocking work is not a local subscriber.

A bounded, nonblocking queue cannot guarantee unlimited lossless delivery to a permanently stalled sink. Resolve this explicitly: before starting another public step, reserve critical credits for its maximum publication batch (64 slots covers one step's certificate/lifecycle/diagnostic transitions, loss summary and at most 16 admitted control requests). Stop admitting new control requests when credits are unavailable. If the journal cannot drain and the next batch cannot be reserved, do not start that step: emit `failed{code:telemetry_exhausted}` and a terminal state using emergency capacity, retain pending critical records and stop with a telemetry failure. This is a resource fault, not waiting for a subscriber or silently dropping a certificate. Journal byte quota is fixed at launch (example 64 MiB); reserve space for critical termination and reject runs whose configured minimal critical history cannot fit. Optional routine data never consumes the reserved critical journal quota.

The producer never blocks on a full queue. Journal persistence may fail; a bounded emergency record remains in memory and final CLI error reports the failure. A process crash/power loss before persistence cannot guarantee event durability, and must be reported as an incomplete journal on recovery. Exact-continuation checkpoints require all prior critical records durably published. Do not claim that “never drop” defeats disk failure. In Classic none of these streams or queues exists; in Observed/AgentControlled resource faults are explicit and tested. P11J measures saturation and journal failure, separately from the normal overhead target.

## One checkpoint envelope: P11C A4 and ESLE phase B

`CheckpointManifest` fields: `schema_version:string`, `calc_id:string`, `seq:u64`, `timestamp_utc:string`, `kind:exact_continuation/warm_start/result_export`, `iteration_done:u64`, `state_version:u64`, `policy_revision:u64`, `mode:Classic/Observed/AgentControlled`, `completion_mode:strict/retain_last`, `model_digest:string`, `environment_digest:string`, `build_digest:string`, `execution:ExecutionIdentity`, `policy:Policy`, `audit_head_sha256:string|null`, `electronic_state:ElectronicState|null`, `continuation:Continuation|null`, `arrays:ArrayRef[]`, `capabilities:string[]`. ElectronicState is the single ESLE payload described in [electronic-state.md](electronic-state.md); it does not define a second container. Initial checkpoints can lack an evaluated ElectronicState. Classic final result export has seq=0 and null audit head because it has no event stream; this does not fabricate historical events. Live checkpoint controls are AgentControlled-only, and an explicit post-run export is an I/O operation, not telemetry serialization inside Classic execution.

`ExecutionIdentity` contains compiler/version, lockfile hash, source hash including dirty patch hash, target/CPU, numerical backend/version, requested/effective k workers, backend thread counts and seed. All are strings except thread counts and seed (u64). `ArrayRef` contains `role:string`, `path:string` (relative, confined to store), `dtype:f64le/complex128le`, `shape:u64[]`, `layout:C`, `ordering_id:string`, `byte_length:u64`, `sha256:string`, `units:string`. Complex payloads interleave real/imaginary; density follows existing CLI C order. Paths are storage references, never commands. Publish payloads first, immutable manifest last using atomic rename; reject partial files, overflow, symlinks escaping store, hash/length/shape mismatch.

| Kind | Required content and compatibility | Promise |
| --- | --- | --- |
| Exact continuation | Mixed next-input density, previous per-k orbitals, full completed numerical history, best/growth counters, next iteration/schedule state, mixer variant/private history and feasibility state, policy/log position, seed/RNG state if used. Identical model/environment, grid/basis/k ordering, occupations, build, backend, thread topology, completion mode and supported mixer state version | Same trajectory at next boundary; reject if any opaque solver/mixer state cannot be captured |
| Warm start | Valid ESLE density and optional orbitals; validate model and basis for direct reuse; explicit separately validated ESLE transfer plan for changed grid/cell | New calc id, reset mixer, counters and certificate; no trajectory promise |
| Result export | Evaluated `ScfOutcome` plus identity, quality and optional arrays | Inspection only; no continuation claim even if density exists |

Exact continuation occurs only between complete solves, so Davidson's in-flight subspace is unnecessary; persist any future cross-iteration solver cache that affects arithmetic. `DensityMixer` presently has no persistence API: P11C must add separately reviewed state extraction after P10D, never infer history from final density. Reject incompatible exact requests; offer an explicit warm-start request, never automatic downgrade. Recompute finite/nonnegative density, charge integral, occupations, hashes and fingerprint on import; validate charge using existing 1e-8 relative gate. An exported unconverged state remains unconverged. Persistence cannot manufacture a certificate.
