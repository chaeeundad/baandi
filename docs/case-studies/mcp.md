# MCP: supervised control contract checks

Status: **partial integration case**, not a measured end-to-end scientific case.

## 1. Problem and motivation

An external client must not acquire numerical authority merely by issuing a tool
call. This case checks host approval, stale observations, bounded transport and
read-only observation in the merged MCP adapter.

## 2. Setup

The embedded `baandi-mcp::Server` test host constructs an in-memory synthetic
local NC potential: valence 2, 201 radial points, spacing 0.06 Bohr, no nonlocal
projectors. One Xx ion occupies a 4 Bohr cubic cell at Gamma, cutoff 2.5 Ha,
two bands, budget two iterations. The complete fixture is
`crates/baandi-mcp/tests/server.rs`; its SHA-256 and build date are in
[data/recorded/provenance.json](data/recorded/provenance.json).
The test does not export a generated UPF-byte hash; no real UPF is used.

## 3. Baseline and external path

The host registers a validated problem under an opaque input ID. A JSON-RPC
initialize/initialized handshake uses MCP 2025-06-18. Creation selects `observe`
or `supervised`; `step`/`run` execute through the shared runtime. There is no
standalone TOML-loading MCP executable. These checks do not compare a Classic
run, so energy difference, matched iterations and overhead are unavailable.

## 4. Request sequence and audit

The executable tests create a supervised calculation, submit beta=0.2, assert
policy_revision remains zero while pending, invoke trusted-host `approve`, then
assert revision one and `AuditOutcome::Applied`. A cancellation is denied by the
host and the initialized state remains. Separate lifecycle checks approve pause,
resume and cancel; export after cancellation has null energy. An observation-only
calculation rejects control. A proposal observed before a step is rejected as
StaleState when approval arrives afterward. The optional client cannot send its
own approval credentials. Read [adapter instructions](../../crates/baandi-mcp/README.md).

The new [test log](data/recorded/mcp.log) records executed test names and outcomes.
Audit assertions execute in memory; this log is **not** a serialized wire or audit
transcript. A persisted transcript remains required for full P11L acceptance.

## 5. Results

Nine integration tests passed with none ignored. The bounded-run check asserts
two iterations, one drained event at limit one, and finite exported energy.
It does not assert certified convergence or export the energy value. No numerical
energy, Classic delta or adapter-overhead number is inferred from passing tests.
Process duration in provenance includes Cargo/test overhead, not just SCF.

## 6. Failure cases

Executed checks cover invalid beta, stale approval, observation-only control,
unknown/authority fields, malformed/oversized frames, notification misuse,
registration and iteration quotas. They fail closed. See named checks in the log.

## 7. Limits

One synthetic system, embedded synchronous host, no live agent, network transport,
second system, production approval UI or persisted wire transcript. Host approval
is explicit and optional-agent operation is exercised. A full comparative case
requires a case-owned embedding harness; adding that remains outstanding.

## 8. One-command reproduction

```sh
python3 -I scripts/run_case_studies.py mcp --verify
```

Here verification checks all nine contract tests, not historical energy numbers.
