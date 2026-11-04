TASK_ID: P11F-MCP-SERVER
GOAL: Thin MCP server over the runtime (directive C1, C2): typed tools (create calculation from validated input, step/run N, status, get events, get snapshot, propose policy change, pause/resume/cancel, export result), Observe and Supervised intervention modes only; no shell, no filesystem beyond a configured workspace, quotas, audit log.
OWNED_FILES: crates/baandi-mcp/ (new crate outside the numerical core; needs ADR/ownership entry and dependency review for the MCP/Rust SDK crate: license, maintenance, Apple Silicon)
DEPENDENCIES: P11B, P11D merged; P9A design consulted (one MCP server, P9A consumes it)
ACCEPTANCE_CRITERIA: tool schemas generated and tested; unsafe policy rejected through the tool path with the same audit as the library; Autonomous mode not implemented (documented as later); core crates gain no dependency on the server crate.
GIT: branch task/P11F-MCP-SERVER, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
ARCHITECT DECISION 2026-11-04 (dependency review): no MCP SDK crate and no async runtime. Implement MCP over stdio as minimal JSON-RPC 2.0 with serde and serde_json only (already adopted); tool schemas are written as typed structs plus JSON Schema documents checked by tests. Reason: small audit surface, no tokio, deterministic, license-clean. A later ADR may adopt the official SDK if protocol coverage requires it. Add the crate to ownership.toml and docs/dependencies.md (no new third-party crate). Core crates must not depend on baandi-mcp.
