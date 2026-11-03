# dft-engine

**Purpose:** Orchestrates SCF, occupations, total energy and its audit (ion-ion/Ewald, double counting).

**Invariants to preserve and test:** Every energy term has a unit test and appears in results; no hidden constants.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-engine`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.

## Incremental SCF (P11A)

`ScfSession` owns the single SCF iteration implementation (`session.rs`, a child
of `scf` to keep problem internals private, re-exported at the crate root).
Creation consumes a validated `ScfProblem`, `CompletionMode`, and
`ExecutionConfig`; `initialize` performs mixer validation once. `step` executes
one ordinary iteration including its nested certificate and applicable mixing.
It returns `StepReport` with Ready/Converged/NotConverged, or the existing
numerical error. Calls before initialization or after termination reject without
mutation. Only a passed certificate sets convergence. `result` borrows an outcome
when available; `into_result` consumes it. Strict exhaustion has history but no
result; the legacy driver converts its status into `ScfError::NotConverged`.
Retained exhaustion skips final mixing, exactly as the previous implementation.

Execution mode is immutable here. Classic creates no queue, event payload or
telemetry timestamp. Observed and AgentControlled share bounded observation;
AgentControlled does **not** advertise control, pause, replay or checkpoint
capabilities (P11B/C). Assign a distinct `ExecutionConfig.calc_id` for each
observed calculation; the default ID is for standalone use. Snapshots borrow
options and scalar energy terms without copying history/arrays; after a failed
partial step `record_iteration` can exceed the completed `iteration`.
`density`, `history` and `result` are explicit larger-data accessors.

`events.rs` uses preallocated routine/critical partitions. Routine overflow drops
new iteration events, counts each loss, and publishes a loss summary before the
next critical event or every 64 iterations. Pop merges by sequence. Six critical
credits are reserved before each step; two additional emergency slots permit a
loss summary and `TelemetryExhausted` failure without waiting or losing pending
critical events. A consumer drains after steps and owns any persistence/fan-out.
There is no durable journal or background worker in the engine. `Observer` is
invoked explicitly by `EventQueue::drain_to`, never by an eigensolve or step;
adapters must keep slow external subscribers off their numerical producer.
Timestamp values use `SystemTime` at publication; RFC3339/JSON formatting belongs
to adapters. Events contain no array payloads. Sequence gaps have loss summaries;
state versions advance once for initialization and each completed/failed step.

Validation: `scripts/check.sh dft-engine`; the complete equivalence matrix is
ignored in debug to keep gates fast. Run
`cargo test --release -p dft-engine --test stepwise_equivalence -- --include-ignored --nocapture --test-threads=1`.
It covers Gamma, four k points, skew multi-species/nonlocal-plus-local fixtures,
Linear/Pulay/KerkerPulay, threads 1/4, all modes, full array bits and scientific
history/certificates excluding wall-clock durations. Release fingerprints also
allow comparison with the pre-refactor implementation. Timing samples are
interleaved after warmup; shared-host timings do not establish quiet-node targets.
