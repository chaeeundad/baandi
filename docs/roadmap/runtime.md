# Agent-native runtime roadmap

Status (2026-11-05): P11R/A/B/C/D/E/F/G/H merged; P11I specifications delivered; P11K preflight and P11M stagnation studies delivered; P11L case studies partial (MLIP/coupling excluded). P11J local benchmark delivered, dev08 overhead verdict provisional until raw evidence arrives. [ADR 0028](../decisions/0028-agent-native-runtime.md) and [architecture](../architecture/agent-native-runtime.md) define the contract. All phases preserve existing numerical defaults unless a separate measured ticket authorizes changes.

| Ticket/phase | Scope and existing seam | Dependency and exit gate |
| --- | --- | --- |
| P11R | These documents, concrete schemas and overlap decisions | Docs checks and links; no Rust/Python code |
| P11A A1–A3 | Extract session around `ScfProblem::run_impl`, projections of `IterationRecord`, scalar events in driver | After P10A/P10D/P10C/P7PC merge; old versus stepwise bit identity, performance and legacy error semantics |
| P11B B1–B3 | Safe-point controls and diagnostics using `ScfOptions`, `DensityMixer`, certificate observations | P11A; transactional patch/rejection, replay, exact rule tests, no convergence override |
| P11C A4 | Unified checkpoint envelope, ESLE payload, CLI store | P11A and P10D mixer state contract; exact/warm/export round trips and incompatible-load rejection |
| P11D B4 | Recovery rules around session and diagnostics | P11B/P11C; safety adversarial cases, replay and bounded resource use |
| P11E | Three-arm self-healing protocol | P11D, frozen failing fixture/reference; energy difference <=1e-10 Ha, certificates and complete audits; broader instability claim needs additional evidence |
| P11F C1–C2 | One MCP server consuming driver contracts | State/control gates; auth/quotas/schema checks; P9A L4 consumes it |
| P11G D1 | Optional ASE adapter over `ScfOutcome`/session | P11A, dependency review; unit/cache tests, explicit unavailable forces/stress |
| Future D2/D3 | MLIP and external solver adapters | Separately scoped tickets, ADR 0017 physical gates; no implementation in P11R |
| Future E1/E2 | Provenance store and workflow DAG | Stable identities, state/control and adapter contracts; deduplication, lineage and resource tests |

P11H CLI work and independent design work do not authorize concurrent edits to scf.rs. P11A and its P10A/P10D/P10C/P7PC dependencies are merged; P9A documents are present. Keep one validator and one server; do not let each frontend implement its own safety rules.

Implementation acceptance artifacts must include actual commands/exit codes, input/build hashes, thread topology, numerical comparisons and scope limitations. Update README validation through its generator only when new measured implementation evidence exists. No numerical implementation was made in P11R, so no README result or new baseline is justified. Dashboard/session updates belong to the orchestrator's artifact ownership; this ticket is confined to the listed documentation paths.

A crate split requires a later ownership/ADR change, two independent consumers or demonstrated adapter-driver duplication, a stable contract and no upward numerical dependencies. General recovery defaults require P10G's diverse cell/band/k-point/projector/gap evidence. The small synthetic demo alone is insufficient.

## Three-mode addendum ownership and sequencing

Status reconciliation: the current common contract is `tickets/P11-COMMON.md`; the runtime decision is ADR 0028. The earlier draft-number and missing-ticket notes are superseded.

| Ticket | Owns | Consumes from P11R; does not duplicate |
| --- | --- | --- |
| P11H-CLI-UX | Actual CLI aliases/mode flags, status/explain, human and JSON-lines sinks, approval UI; `docs/user/cli-modes.md` | Three mode semantics, typed events, capability/error ranges; no CLI-specific controller or second MCP server |
| P11I-COUPLING-API-SPEC | `workflow-api.md`, `coupling-api-lowlevel.md`, `active-learning-contract.md`, external-coupling efficiency and physical tests | Immutable runtime request lifecycle, state/provenance IDs and mode/control boundary; P11R's earlier duplicate D1/D2/D3 schemas are withdrawn |
| P11J-RUNTIME-BENCH | `scripts/bench_runtime.py`, measured runtime-equivalence results, `tests/reference/runtime/`, optimization-log append | Target definitions and measurement/test matrix; append measured results after implementation, never mark design targets as passed |

Status: P11H and the P11J tooling are merged, and P11I design documents are present. Merged algorithm capabilities are inventoried in [runtime architecture](../architecture/agent-native-runtime.md).

Refined gates: Classic has no event allocations or agent; Observed never waits for subscriber; AgentControlled without intervention stays bit-identical; stale control and approval checks precede mutation; unavailable-agent default preserves policy; no mode transition changes scientific settings; routine loss is counted; critical capacity exhaustion is explicit. P11J measures the <0.5%/<1%/<1% targets on matched normal and short cases, and reports failed/inconclusive targets honestly.
