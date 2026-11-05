# baandi-mcp

Current ownership (2026-11-05): merged P11F adapter, Observe and Supervised only.
Trusted host registers validated problems and approves proposals; a standalone
stdio client cannot grant itself authority. No general file-loading host or MLIP
adapter is implied. ADR 0029 records the boundary; runtime policy remains owned
by dft-engine.

Adapter only: stdio JSON-RPC, typed tool inputs, bounded host-owned runtimes.
No numerical implementation, shell, network, arbitrary file paths or autonomous
approval. Host registration accepts an already validated ScfProblem; approvals
are trusted Rust calls, never MCP parameters. Keep all runtime validation/audit
in dft-engine. Verify with `scripts/check.sh baandi-mcp`.
