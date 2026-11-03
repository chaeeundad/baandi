TASK_ID: P11F-MCP-SERVER
GOAL: Thin MCP server over the runtime (directive C1, C2): typed tools (create calculation from validated input, step/run N, status, get events, get snapshot, propose policy change, pause/resume/cancel, export result), Observe and Supervised intervention modes only; no shell, no filesystem beyond a configured workspace, quotas, audit log.
OWNED_FILES: crates/baandi-mcp/ (new crate outside the numerical core; needs ADR/ownership entry and dependency review for the MCP/Rust SDK crate: license, maintenance, Apple Silicon)
DEPENDENCIES: P11B, P11D merged; P9A design consulted (one MCP server, P9A consumes it)
ACCEPTANCE_CRITERIA: tool schemas generated and tested; unsafe policy rejected through the tool path with the same audit as the library; Autonomous mode not implemented (documented as later); core crates gain no dependency on the server crate.
GIT: branch task/P11F-MCP-SERVER, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
