# dft-engine

## Current ownership (2026-11-05)

Owns SCF orchestration/session, `events`/`state_view`, `runtime`/`policy`,
`diagnostics`, and in-memory `checkpoint` payload/validation; also NSCF and
band-density orchestration. File persistence stays in CLI, transports in adapters.
The P11A/B notes below are historical implementation detail: P11C now supplies
checkpoint/restore APIs, including runtime/audit restoration. No cross-build exact
restart guarantee or engine disk I/O is implied. See `docs/checkpoint-format.md`.


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

Direct sessions expose no public policy/mode mutation. Classic creates no queue,
event payload or telemetry timestamp. Observed and AgentControlled share bounded
observation. Use the P11B `Runtime` wrapper for control, pause and policy replay;
exact checkpoints remain P11C and are explicitly unavailable. Assign a distinct `ExecutionConfig.calc_id` for each
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

## Runtime policy control (P11B)

`runtime::Runtime` exclusively owns a `ScfSession`. Public access to the underlying
session is read-only, apart from draining telemetry. No runtime object, policy,
audit buffer or new per-iteration branch is installed in ordinary `run` or direct
`ScfSession::step`. The original advance/step implementation is unchanged.
The only shared setup refactor extracts construction of the existing mixers.
`tests/stepwise_equivalence.rs` is unchanged.

Lifecycle: Created -> Initialized -> Running -> Initialized, with certificate-only
Converged, budget Exhausted, numerical Failed and controlled Cancelled terminal
states. PauseRequested means an uncommitted request; a boundary first completes
numerics and then validates that request before entering Paused. Terminal numerics
win over pending pause/cancel; approved cancellation precedes other controls.
Resume activates the existing continuation through the logical Restarting state,
without initialization or resetting numerics. Restarting is not checkpoint import.
The public pure transition table does not mutate a runtime or confer convergence.

`begin_step`/`finish_step` allow scheduling and queuing before the synchronous SCF
step. There is one mutation owner; transport threads enqueue to their driver, which
calls `submit` at the next accessible boundary. There is no transport/thread worker
in the engine and no interruption inside an eigensolve. `step` processes controls
before/after an iteration. With no approved request it continues the current
policy; pending human approvals and missing agents never block a running solve.
`safe_point` rejects during Running/PauseRequested/Restarting. All proposals are
validated against the completed boundary, not the transient Running projection.

ControlRequest schema 1 carries calculation/request identity, observed state and
policy revision, a tagged action and parameters, at most 8 typed predicates,
iteration/wall expiry, stale policy and a reason capped at 512 UTF-8 bytes. Labels
are capped at 128 graphic ASCII bytes; patches have at most 10 changes. Wall expiry
and approval expiry are also converted to monotonic deadlines at admission.
Stale state can be revalidated only with unchanged policy revision, explicit
permission and passing current predicates; future versions always reject.
Unavailable diagnostics fail closed. Unknown physical-model interventions are
unrepresentable; checkpoint requests return UnsupportedCapability.

Default authority is supervised: submit with PendingApproval, then let an
AUTHENTICATED HOST call approve for that owned immutable request ID. Authority and
Approval are trusted Rust adapter arguments, not fields to deserialize from an
agent body. Approved requests can also be submitted by the host directly. A changed
payload or requester with a reused ID conflicts; identical retries reuse the prior
disposition without another application. Invalid/oversized input and capacity
failures reject before admission; every admitted, disposed proposal has an audit
record. Pending proposals can be denied or expire and never reserve a numerical
mutation lock. Host transport rate limits and authentication remain adapter duties.

Optional launch-time autonomy is a closed allow-list (beta ranges/nonincrease,
reset, implemented solver selection, tightening, bounded budget and lifecycle).
There are at most 16 rules and 32 total applications, with per-rule counts,
iteration cooldown and version-lag bounds. No rule means deny. Mixer strategy
changes use supervised approval in this initial profile. Administrative mode
changes always need host approval. Downgrading to Observed revokes autonomy and
rejects pending controls; counters persist and upgrading does not silently restore
autonomy. Classic can change mode only before initialization. Observed and
AgentControlled may switch at Initialized/Paused boundaries; a switch preserves
numerical settings, arrays, mixer history and policy revision. Same-mode requests
are idempotent and do not advance the state version.

`policy::prepare_policy` validates immutable options. `ScfSession` has crate-private
prepare/commit APIs used by Runtime only. A mixer change (including beta alone) or
explicit reset constructs and validates an empty replacement before commit; it
preserves density and orbital guesses and increments the reset count. Other changes
preserve mixer history. Hybrid selection retains the merged implementation's
Davidson fallback for buffer-band policies; it does not change band counts.
No schedule or divergence intervention is allowed by P11-COMMON; those enum
variants return named restriction errors. Overlapping mixer assignments reject;
use one complete Mixing value for a combined kind/beta/history change. History is
capped at 64 and the total budget at 10000, and a budget patch must increase it
before exhaustion. Exhausted sessions are not resumed without P11C continuation.
At most one numerical transaction is applied per completed iteration boundary.

Tolerance changes need launch-time, evidence-backed positive lower bounds and
must tighten the current values. They cannot exceed the existing certificate
threshold defaults: density 1e-7 Bohr^-3/2, energy/eigen 1e-8 Ha. Lower bounds limit
excessive tightening, not authorize weakening. There is no unsafe-relaxed escape
hatch; result provenance always records unsafe_relaxed=false. Test floors are
validator fixtures, not scientific defaults. No SCF numerical defaults changed.

Audit is a bounded independent control journal (default 4096, maximum 65536
records), separate from the existing SCF event queue. Each record carries schema,
calculation ID, sequence, SystemTime timestamp, original request (including ID,
reason and observed versions), requester, authority/approval, replay marker and autonomy rule index, complete before/after option
and lifecycle projections, disposition and reset semantics. Accepted and Applied
are separate records; rejected transactions have after=before and preserve state
version, policy revision, density, orbitals and hidden mixer state. Submission is
not semantic acceptance. Journal capacity is reserved for every admitted request.

All validation and replacement preparation precede publication. Optional AuditSink
is a trusted adapter contract: append an entire batch atomically or publish none.
The Accepted/Applied batch is published before an infallible in-memory commit under
exclusive ownership. A sink error commits no numerical change, stores a local
JournalFailed rejection and prevents further driving; RuntimeResult exposes that
failure. The default in-memory journal is NOT durable. No disk I/O, serialization,
hashing/authentication, crash recovery or durable sink implementation is claimed.
Those belong to adapters/P11C; a sink must obey its documented atomic contract.

`RuntimeResult` carries optional ScfOutcome, lifecycle state, policy-change count,
unsafe_relaxed, launch control configuration, journal-failure status and the complete bounded audit. Cancellation
or strict failure never fabricates a scientific result. Only the unchanged SCF
fixed-point certificate creates convergence.

Replay consumes Applied records in recorded sequence at exact iteration boundaries,
checks full before/after control projections, and revalidates the numerical policy.
Wall/approval times are ignored, not extended. Autonomous counters and cooldown
history are reconstructed from the recorded rule indices; replay cannot replenish
a consumed grant. Use the same immutable model,
launch config, build/backend/CPU and thread count; this is policy replay from the
initial state, not exact checkpoint restoration or a cross-platform promise.
Input journal authenticity is the adapter's responsibility. No current predicates
or new agent decisions are substituted for the recorded committed decisions.

Validation commands:
- `scripts/check.sh dft-engine` and `scripts/check.sh`.
- `cargo test --release -p dft-engine --test runtime_policy -- --include-ignored --nocapture --test-threads=1`.
- `cargo test --release -p dft-engine --test stepwise_equivalence modes_mixers_threads_and_species_are_bit_identical -- --ignored --nocapture --test-threads=1`.

The runtime tests cover the complete 10-state/10-trigger matrix, rejection followed
by identical future trajectory (including hidden mixer history), stale/expired
requests, approval, bounded autonomy, journal failure/capacity, mode revocation,
silent-agent fallback, solver/tolerance interventions and full-array replay.
The ignored release replay matrix uses cubic Gamma, four k points and a skewed
multi-species nonlocal-plus-local fixture; it compares energy/density/orbital bits,
full scientific history and certificates, excluding only timings. No silicon-only
default or optimization was adopted, and no new QE or overhead claim is made.

### Verified P11B evidence (2026-11-04)

Release runtime suite: 18 passed, 0 failed (17 ordinary debug tests and one
release-only replay matrix). The transition test checks all 100 state/trigger
pairs. The unchanged stepwise matrix checks 18 fixture/mixer/thread combinations
against all three execution modes (54 comparisons). All 18 complete-state
fingerprints also match an independently built archive of P11A commit
`04167620eca63be4e2a3e60455414db5accbd64f` using the identical matrix command above.
The old and current session code from `advance` through the end of the file is
byte-identical; no control branch/allocation/logging was added to ordinary steps.
These are equivalence checks, not measured overhead targets.

Each replay fixture applies three numerical transactions at boundaries 1, 2 and 3,
including Linear -> KerkerPulay and explicit history reset/budget extension, plus
pause/resume. Energy, complete density/orbital arrays, scientific history and
certificate are bit-identical on replay. Values below are synthetic engine
results, not QE references.

| Fixture | Energy (Ha) | Iterations | Certificate density (Bohr^-3/2) | Occupied eigenvalue change (Ha) |
| --- | ---: | ---: | ---: | ---: |
| Cubic Gamma, nonlocal | -0.62379772639925046 | 14 | 1.614491029e-9 | 1.383107906e-9 |
| Four k points, nonlocal | -0.51412521505579167 | 14 | 4.017338130e-10 | 2.841891056e-9 |
| Skew multi-species, nonlocal + local | -1.52491411485654660 | 15 | 3.049026055e-9 | 9.491289278e-9 |

Numerical defaults, certificate thresholds and existing README reference evidence
are unchanged. README/figure regeneration and its three specification tests pass;
the generated files have no diff. Adding a new README result type or dashboard
JSON entries requires the corresponding orchestrator-owned reporting files.

## Collinear SCF (P13B / MF03)

`scf::collinear` owns the additive `CollinearProblem`, `SpinSector`, result and
fresh-map certificate. Legacy `ScfProblem` constructors pass no spin sector;
no change to their iteration arithmetic or result shape is intended. Shared
construction accepts an internal channel occupied count only for the new API.
There is one solver implementation, called for each channel; density/audit
buffers are reserved once. One common feasibility fraction applies to both
existing mixer proposals. Only raw charge and spin residuals certify convergence;
per-channel gaps apply only to nonempty channels. Output energy/eigenpairs/density
are one state, with fresh certificate metrics separate. No spin runtime or exact
checkpoint support is implied. See `docs/user/collinear-spin.md` for bounds,
measured evidence and limitations.

The geometry/mesh test exceeds the debug budget and is ignored by default:
`cargo test --release -p dft-engine --test collinear -- --include-ignored --nocapture`.
Run `scripts/check.sh` and the unchanged SCF benchmark identity check too.

## MF06 NLCC (P13C)

ScfProblem owns optional immutable core density and reusable XC input storage.
NLCC-off dispatch retains old functional arithmetic. The output-density core XC
is evaluated once before k-weighted audits; Hartree and band double counting use
only valence. `EnergyAudit.nlcc` records core charge, core-potential integral and
valence-potential integral; only `xc` contributes XC energy to total E.
Collinear channels share the spin-neutral convention. NLCC v1 checkpoints reject
on capture and restore; CLI serialization also rejects an active NLCC ledger,
while legacy bytes remain unchanged. Local-potential export uses the same core.
Run the ignored release commands in docs/cli.md; real smoke failures are explicit
and must not be advertised as converged or qualified materials.

## MF04 finite temperature

`finite_temperature.rs` owns prepared occupation arrays and a separate scalar FD
driver; collinear uses the same spin map with a shared or fixed-M number solve.
Do not feed FD into legacy sessions/checkpoints, fixed certificates or a band-only
occupation view. Every retained band is strict and the fresh certificate repeats
the entire map. Direct and band E must use identical fractional weights and
subtract the same TS once. Fixed branches retain their arithmetic. See
`docs/validation/occupations.md` for evidence and limitations.

## MF11 crystal symmetry (P19A)

`enable_symmetry` re-detects the group against stored atoms, checks exact reduced sampling and compatible FFT quadrature, and installs the prepared density projector. Both scalar and spin density builds and fresh certificates use it. Noninvariant input densities and checkpoint identity without symmetry reject.

See `docs/validation/symmetry.md` and `tests/reference/symmetry/` for commands,
observed failures and validation limits. No QE symmetry-on oracle is available.

## Atomic forces (P20A)

`forces(&ForceState)` is post-convergence only. It rejects standard/uncertified
states and reports the complete NC derivative sum. The borrowed cell, atoms,
UPFs, options, density floor and symmetry operations must be the exact original
model. A private setup-time descriptor independently compares these fields and
exact parsed pseudopotential content before force work. File SHA-256 provenance
remains at the adapter boundary; the engine does not receive original file bytes. Do not alter
SCF kernels to calculate forces. The adapter must not reuse a certificate after
mutating result arrays. See `docs/validation/forces.md` for force term ownership,
FD stationarity, independent oracles and versioned force schemas.

Fast tests: `cargo test -p dft-engine --test forces`. Ignored displacement tests
must run in release mode with the commands in that note. Do not weaken the
1e-6 synthetic, 1e-4 real-file or 1e-8 resolved-grid net-force caps.

## Fixed-cell density transfer (P21A / MF07 minimum)

`ScfSession::initialize_density_transfer` is opt-in initialization only. Compare
independent setup model descriptions with ordered atom positions removed; retain
all other identity, sampling, occupation, symmetry and grid checks. Validate
shape, finite/XC-feasible values and relative charge <=1e-8 before mutation.
Do not renormalize or restore orbitals/mixer/history/certification. The caller
records `warm_transfer`, never exact replay. Ordinary initialization and iteration
kernels remain unchanged. Run `cargo test -p dft-engine --test density_transfer`
and the release ignored three-structure warm/fresh energy oracle.
