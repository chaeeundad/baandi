TASK_ID: P11K-AGENT-PREFLIGHT
GOAL: An optional agent layer on top of the deterministic input preflight (P8G, ADR 0022) that reviews input and structure before a run and proposes fixes. Agent-path only: never in the numerical core, never required, never applies changes by itself.
OWNED_FILES: docs/architecture/agent-preflight.md (new), crates/dft-cli/src/preflight.rs (only to expose findings as typed, stable-code data if not already), crates/dft-cli/tests/preflight_findings.rs (new), docs/cli.md (append), docs/INDEX.md
DEPENDENCIES: P8G merged, P11H CLI merged (for `baandi check --explain`)
DESIGN:
1. Authority split: deterministic preflight rules decide pass/warn/fail. The agent only reads typed findings (stable code, TOML key path, evidence values) plus a read-only structure summary, and returns advisory text and candidate edits as a patch proposal. It cannot change the verdict.
2. Structure review additions that stay deterministic and testable: minimum interatomic distance, overlapping atoms, cell volume and skewness, vacuum size for slabs, electron count parity, k mesh versus cell lengths, cutoff versus pseudopotential suggested cutoff, species/element identity. Each finding has a code, severity, evidence, and a documented fix hint.
3. Agent proposals are diffs on the input file, shown to the user, applied only on explicit approval; restricted fields (pseudopotential, XC, structure, electron count, k mesh, cutoff) are proposed, never auto-applied, and the proposal records the reason.
4. No network call from the CLI; the agent is an external client (MCP/Python) consuming the findings JSON. Without an agent the findings are still printed as plain text.
5. Provenance: the proposal and the user decision are recorded next to results.json (additive).
TESTS: every finding code has a positive and a negative fixture (synthetic UPF only); findings JSON schema stable; proposal format validated; no change in numerical results for valid inputs.
ACCEPTANCE_CRITERIA: findings are machine-readable and deterministic, an agent-free run behaves as before, design doc states what the agent may and may not do.
NON_GOALS: natural-language input (P9A), automatic fixing, pseudopotential download.
GIT: branch task/P11K-AGENT-PREFLIGHT, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
EXTENDED 2026-11-04 by architect: OWNED_FILES also include crates/dft-cli/src/{config.rs,run.rs,main.rs,lib.rs} (minimal wiring, list the lines changed), crates/dft-cli/tests/*, docs/user/cli-modes.md, results.json additive preflight field. Implement the deterministic structure rules, findings JSON on failure, `check --explain`/`validate --json` wiring, proposal validation and provenance. Valid inputs must produce unchanged numerical results.

