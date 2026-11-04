# BAANDI MCP runtime adapter

This is an embeddable synchronous stdio MCP server, not a TOML/UPF importer.
The trusted host constructs `ScfProblem` using existing validation, registers it
with `Server::register_input`, then calls `serve_stdio` (or `serve` for another
owned stream). `Server::new(workspace, quotas)` canonicalizes an existing workspace.
No subsequent tool reads or writes any file. Input IDs are opaque single-use
catalog keys, not paths. The host must preflight memory/cost before constructing
problems and enforce process time/memory limits. No network or shell is used.

Initialize with protocol `2025-06-18`, clientInfo and capabilities, then send
`notifications/initialized`. `tools/list` returns the eleven generated schemas in
`schemas/tools.json`. `create_calculation` takes input_id, calc_id, and mode
`observe` or `supervised`; it initializes the runtime. `step` completes one SCF
iteration; `run` takes 1..100 steps and stops on terminal or paused state.
`status` and `get_snapshot` return the same bounded scalar boundary projection.
`get_events` drains up to 256 events and returns the complete bounded runtime
control audit. Event sequence/timestamps and explicit units are retained.

`propose_policy_change` initially supports mixing_beta only. It and pause/resume/
cancel require request_id, observed_state_version, observed_policy_revision,
max_iteration and reason. They return pending_host_approval until the authenticated
host invokes `Server::approve`, or a runtime rejection disposition. There is no
MCP approval tool or authority field. To integrate an operator UI, drive `handle`
or `serve` in bounded host dispatch intervals and call approve/deny between them;
`serve_stdio` itself blocks on stdin. Approval never bypasses policy validation.
Unsafe proposals remain pending until approval, then receive the same library
rejection audit. Inspect the audit after approval; `approve` returning Ok means
safe-point processing completed, not that the proposed change was applied.
Observe rejects interventions. Autonomous mode is deferred and unrepresentable.

Default quotas: 8 combined registered inputs/calculations, 4096 frames, 1000 total
iteration attempts, 65536 bytes per frame. Runtime defaults separately bound
pending proposals (16), audit (4096 records), routine events (256) and critical
events (128 plus runtime emergency slots). Exhaustion fails closed; no session
is silently evicted. Drain events between runs to avoid critical queue exhaustion.
Synchronous work cannot be cancelled in the middle of an iteration. Malformed
and unknown input fields reject before dispatch. Oversized stream frames close
the transport; notifications cannot invoke tools. No transport-level retries of
step/create are guaranteed idempotent; control request IDs use runtime semantics.

Audit is in-memory, not durable; the host can inspect the original typed records
through `runtime(id).audit()`. Wire audit includes searchable metadata and the
full Rust debug representation of each record (diagnostic, not a stable replay
format). Export requires a terminal calculation and returns scalar JSON plus
CSV (`iteration,energy_hartree`). It explicitly omits arrays/restart data; a failed
or cancelled calculation may have no energy. A finite energy never implies
convergence; only the runtime's certificate-backed Converged state does.

No standalone input-loading executable is provided: P9A/CLI hosts reuse this
adapter with their validated input and trusted approval UI. `serve_stdio` is the
production transport entry point. MCP resources, prompts, HTTP, batch requests,
subscriptions, progress, roots and transport cancellation are not advertised.

Regenerate schema documents:

```sh
cargo run --offline -p baandi-mcp --example tool_schemas > crates/baandi-mcp/schemas/tools.json
scripts/check.sh baandi-mcp
```

The schema macro generates structs and schema properties together. Semantic
policy ranges deliberately remain in Runtime so unsafe beta reaches its audit.
Protocol definitions: https://modelcontextprotocol.io/specification/2025-06-18
