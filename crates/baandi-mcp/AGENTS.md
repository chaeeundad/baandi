# baandi-mcp

Current ownership: R8d reference runtime adapter, Observe and Supervised only.
Trusted host registers validated problems and approves proposals; a standalone
stdio client cannot grant itself authority. No general file-loading host or MLIP
adapter is implied. ADR 0029 records the boundary; runtime policy remains owned
by dft-engine.

Adapter only: stdio JSON-RPC, typed tool inputs, bounded host-owned runtimes.
No numerical implementation, shell, network, arbitrary file paths or autonomous
approval. Host registration accepts parsed ReferenceInput (cell, atoms, potentials, options
and reference Request); the engine validates the start during creation. Approvals
are trusted Rust calls, never MCP parameters. Keep numerical policy validation and lifecycle transitions in the reference
engine runtime. Transport approval, stale/expiry checks and bounded JSON audit
are owned by this adapter, as required by the reference runtime contract. Verify with `scripts/check.sh baandi-mcp`.
