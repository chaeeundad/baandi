# baandi-mcp
Adapter only: stdio JSON-RPC, typed tool inputs, bounded host-owned runtimes.
No numerical implementation, shell, network, arbitrary file paths or autonomous
approval. Host registration accepts an already validated ScfProblem; approvals
are trusted Rust calls, never MCP parameters. Keep all runtime validation/audit
in dft-engine. Verify with `scripts/check.sh baandi-mcp`.
