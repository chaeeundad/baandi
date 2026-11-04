# ADR 0029: Minimal supervised MCP stdio adapter

Date: 2026-11-04. Status: accepted architectural direction supplied in P11F;
implementation pending verification/review. Supersedes the SDK review placeholder
in that ticket; extends ADR 0028 without changing numerical contracts.

One `baandi-mcp` adapter depends downward on `dft-engine`. P9A consumes this
server; it does not implement a second controller. The workspace and ownership
entry were pre-registered by the architect. No core crate depends on this crate.

Use stdio JSON-RPC with serde and serde_json, already adopted under MIT OR
Apache-2.0. No MCP SDK, async runtime, network stack or new third-party package.
This keeps the adapter's audit surface small and its ownership synchronous.
A later ADR may adopt the official SDK if broader protocol coverage requires it.

Host registration accepts a constructed, validated ScfProblem under an opaque,
single-use input ID. The host preflights resource requirements before preparing
it. Tool arguments cannot supply UPF paths, shell commands, executable callbacks,
physical model patches or approval credentials. Runtime initialization still
validates the actual session. The adapter performs no filesystem I/O after
canonicalizing the configured workspace at startup; exports are inline JSON/CSV.
File import/export integration belongs to a later CLI host, with its own path
confinement contract. This avoids symlink races and ambient file authority.

Observe maps to Observed, Supervised to AgentControlled with no autonomy policy.
All control proposals, including lifecycle changes, go through Runtime.submit
with PendingApproval. Only the trusted embedding host calls approve/deny. A
pending proposal is not semantically accepted; unsafe beta is rejected with the
same runtime audit after host approval, never applied. Existing runtime expiry,
stale-version checks, quotas and transaction semantics remain authoritative.
Autonomous is not implemented or deserializable in this adapter.

Bound frames, total requests, resident inputs/calculations and total iteration
attempts. Host-selected problem sizes and OS time/memory limits remain necessary:
a synchronous SCF iteration cannot be interrupted by another stdio request.
Runtime events and in-memory audits retain their existing bounded limits and
failure behavior. No durability or crash-recovery claim is made.

Protocol subset: MCP 2025-06-18 initialize/initialized, ping, tools/list and
 tools/call over newline-delimited UTF-8. No HTTP, resources, prompts, subscriptions,
server requests, batching, transport cancellation or progress notifications.
Only tools are advertised. Generated input schemas are checked against the typed
Rust declarations and committed schema document; numerical policy validation is
never duplicated in the transport. Result export is a scalar summary, not an
exact checkpoint, density export or substitute for the CLI scientific writers.

References: [MCP transport specification](https://modelcontextprotocol.io/specification/2025-06-18/basic/transports),
[tools](https://modelcontextprotocol.io/specification/2025-06-18/server/tools),
[P9A interface](../architecture/user-interface.md).
