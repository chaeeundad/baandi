# BAANDI MCP runtime adapter

This is an embeddable synchronous stdio MCP server, not a TOML/UPF importer.
The trusted host constructs `ReferenceInput` from parsed cell, atoms, potentials,
SCF options and a reference `api::Request`, registers it
with `Server::register_input`, then calls `serve_stdio` (or `serve` for another
owned stream). `Server::new(workspace, quotas)` canonicalizes an existing workspace.
No subsequent tool reads or writes any file. Input IDs are opaque single-use
catalog keys, not paths. The host must preflight memory/cost before constructing
inputs and enforce process time/memory limits. The reference engine validates
the start during creation; failed creation consumes the input ID. No
`ScfProblem` or `ScfSession` is constructed. One dedicated thread per calculation
keeps the reference continuation and solver scopes on their owning thread. No network or shell is used.

Initialize with protocol `2025-06-18`, clientInfo and capabilities, then send
`notifications/initialized`. `tools/list` returns the eleven generated schemas in
`schemas/tools.json`. `create_calculation` takes input_id, calc_id, and mode
`observe` or `supervised`; it initializes the runtime. `step` completes one SCF
map (ordinary, verification or retry); `run` takes 1..100 maps and stops on
terminal or paused state. Verification and retry maps can repeat an iteration
number; `completed_maps` counts all successful map returns.
`status` and `get_snapshot` return the same bounded scalar boundary projection.
`get_events` drains up to 256 events and returns the complete bounded runtime
control audit. Schema version 2 snapshots include the complete `last_step`
report with explicit Ha/Ry units, runtime state, safe-point flag and policy
revision. Bounded map events contain a sequence and snapshot; overflow is
reported by `dropped_events`. Control records are retained independently.

`propose_policy_change` accepts one or more of `mixing_beta`, `mixing_history`,
`threshold_ry` (total residual energy, Ry) and `verification` (boolean). Unknown
legacy knobs and empty patches reject with the reference allow-list. Threshold
relaxation is forbidden; numerical ranges and atomic patch application are
validated by the engine. A patch requires a live safe point; there is at most
one successful transaction per map boundary. It and pause/resume/
cancel require request_id, observed_state_version, observed_policy_revision,
max_iteration and reason. They return pending_host_approval until the authenticated
host invokes `Server::approve`, or a runtime rejection disposition. There is no
MCP approval tool or authority field. To integrate an operator UI, drive `handle`
or `serve` in bounded host dispatch intervals and call approve/deny between them;
`serve_stdio` itself blocks on stdin. Approval never bypasses policy validation.
Unsafe proposals remain pending until approval, then receive the engine error in the transport
audit. `approve` returns an error for invalid policies, stale requests and expired
approvals. Pause/cancel approved at an unsafe boundary wait in the engine until
the next safe return; terminal numerical delivery wins. Resume requires paused
state. Control IDs are single-use, and host approval never revalidates stale
state, even if its legacy `allow_revalidation` flag is set.
Observe rejects interventions. Autonomous mode is deferred and unrepresentable.

Default quotas: 8 combined registered inputs/calculations, 4096 frames, 1000 total
map attempts, 65536 bytes per frame. The adapter also bounds
pending proposals (16), request identities and audit records (4096 each), and
map events (256). Control-capacity exhaustion rejects new proposals. Event
overflow drops map events with an explicit counter; audits are not evicted.
Synchronous work cannot be cancelled in the middle of an iteration. Malformed
and unknown input fields reject before dispatch. Oversized stream frames close
the transport; notifications cannot invoke tools. No transport-level retries of
step/create are guaranteed idempotent; control request IDs are single-use within each calculation.

Audit is in-memory, not durable; the host can inspect JSON records through
`runtime(id).audit()`. Audit includes request, action,
host approval identity, outcome and before/after state versions; it is not an
exact policy replay/checkpoint format. Export requires a terminal calculation and returns scalar JSON plus
CSV (`iteration,energy_hartree`). It explicitly omits arrays/restart data; a failed
or cancelled calculation may have no energy. A finite energy never implies
convergence; only the runtime's Converged state does. Verification can be explicitly disabled
by the registered options or an approved policy; Converged alone is not a claim
that verification ran.

No standalone input-loading executable is provided: P9A/CLI hosts reuse this
adapter with their validated input and trusted approval UI. `serve_stdio` is the
production transport entry point. MCP resources, prompts, HTTP, batch requests,
subscriptions, progress, roots and transport cancellation are not advertised.

Regenerate schema documents:

```sh
cargo run --offline -p baandi-mcp --example tool_schemas > crates/baandi-mcp/schemas/tools.json
scripts/check.sh baandi-mcp
```

Required control fields use the shared schema macro; optional policy fields
and their at-least-one constraint are checked against the typed declarations. Semantic
policy ranges remain in the reference runtime; its rejection is recorded by
the adapter after host approval.
Protocol definitions: https://modelcontextprotocol.io/specification/2025-06-18
