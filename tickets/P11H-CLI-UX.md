TASK_ID: P11H-CLI-UX
GOAL: User experience for the three modes: `baandi run` (Classic), `baandi run --observe` (Observed: live human-readable and JSON-lines telemetry to a file/stdout/socket), `baandi run --agent` (Agent-Controlled, attaches the control channel), `baandi status <calc>` (inspect a running or finished calculation from its event log/checkpoint), `baandi explain` (human-readable explanation of a configuration, the preflight result and the diagnostics of a run). Reconcile with the actual CLI first (current binary name and subcommands in dft-cli); keep existing subcommands working (alias, not break).
OWNED_FILES: crates/dft-cli/src/ (run/status/explain modules, event sinks), crates/dft-cli/tests/, docs/user/cli-modes.md
DEPENDENCIES: P11A (mode config, events), P11B (control), P8G (preflight) merged
ACCEPTANCE_CRITERIA: Classic path unchanged in output and bit-identical; observed output readable (iteration, energy, residual, solver status) and machine-readable; errors name the field and the valid range; `explain` never invents physics (reports numbers and labelled hypotheses); natural-language-to-spec translation is a design note only: it must produce an inspectable validated spec, never silently guess physical parameters.
GIT: branch task/P11H-CLI-UX, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
