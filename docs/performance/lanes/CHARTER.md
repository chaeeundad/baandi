# Lane charter (user decision 2026-11-20, revised: profiling runs on dev08, not on the Mac)

## Current state — 2026-11-20..24 (issue #26)

The new reference SCF path (issue #25) is the active work; the threads lane is
dropped. Defaults are unchanged by this reconciliation. dev08 is operated by
`baandi-dev08`; the `dft5090` relay is retired. See the
[current operations and protocol map](../fast-loop-protocol.md#current-state--2026-11-2024-issues-25-and-26).
Use fixed-H replay for kernels, common-input replay for density-update
discrepancies, and full convergence for tail, verification or trajectory changes.
Timing uses [the same median for both engines](../comparison-rules.md#statistics-and-delivered-wall);
fast jobs are limited to 600 s, authorized `final-` jobs to 7200 s.

## Where things run (revision)
The Mac (Apple Silicon, no MKL/MPI parity, shared with other work) is for **compiling, correctness and equivalence only. No timing from the Mac may support a performance conclusion.** All profiling and all A/B timing run on dev08 (exclusive socket-0 window, real MPI, MKL) through the orchestrator and the `baandi-dev08` session. Each lane therefore runs in phases:
- **Phase 0 (agent, Mac):** read the code and ledger, add lane-specific instrumentation (timers/counters, active only under PHASE_PROFILE/STAGE_TIMERS or a default-off switch, no cost when off), and write `docs/performance/lanes/<LANE>/profile-request.md`: exact dev08 commands (cells: si64 np8, si32 np8, si16 np8, si32 and si64 one core as relevant; repeat counts; which counters/timers to return and the JSON paths), plus a script `scripts/lane_<lane>_parse.py` that turns the returned results.json into the lane's tables. Verify on the Mac that instrumentation changes no result (bit-identical energy/density) and the parser works on a Mac-generated file. Stop with status PROFILE_REQUEST.
- **Phase 1 (orchestrator, dev08):** profile run in an exclusive window, raw results preserved under /data/sig-dft/lanes/<LANE>/ and relayed back.
- **Phase 2 (agent, Mac code + dev08 data):** diagnose from the dev08 numbers, write `diagnosis.md`, propose >= 2 candidates, implement one behind a default-off switch, validate equivalence on the Mac, write `dev08-request.md` (A/B). Stop.
- **Phase 3 (orchestrator, dev08):** A/B, approval or rejection; a rejected agent returns to Phase 2 with the data.
Duties 2 and 4 below read accordingly: the diagnosis and all timing evidence come from dev08.
Select experiments by the [live protocol map](../fast-loop-protocol.md#protocol-map-live):
fixed-H replay for kernels, common-input replay for density-update discrepancies,
full convergence for tail/verification/trajectory changes. A ten-step screen
cannot accept or reject an effect outside those steps. Full convergence retains
independent fixed-point verification and charges every failed attempt and
recovery to delivered wall. The [final campaign](../final-measurement-protocol.md)
and [median/statistics rule](../comparison-rules.md#statistics-and-delivered-wall)
govern promotion and timing, with the same median for both engines.
Use 600 s fast jobs or authorized `final-` jobs up to 7200 s as specified in the
map; state the class and diagnostic question in every request.


The SCF computation is split into independent units ("lanes"). Each lane has its own Codex agent that diagnoses, designs and implements improvements on its own judgment. dev08 tests each result; the orchestrator (Claude) approves or rejects. A rejected lane agent brings an alternative.

## Agent duties (every lane)
1. Read AGENTS.md, docs/performance/AGENT-BRIEF.md, attempts-ledger.md (grep your lane's terms) and backlog.md. Never repeat a rejected idea without stating the new evidence.
2. **Diagnose first on dev08**, with instrumentation separate from timing runs.
   Record per-call/map/rank costs and arithmetic/bandwidth/latency limits in
   `docs/performance/lanes/<LANE>/diagnosis.md`, with the raw evidence source.

   > Historical wording, superseded as to the machine: **Diagnose first** (Mac release, small and medium cells; instrument with PHASE_PROFILE/STAGE_TIMERS or your own timers): where does your lane's time go (per call, per map, per rank), what is the floor (arithmetic count, bandwidth, latency)? Write the finding as `docs/performance/lanes/<LANE>/diagnosis.md` with numbers.
3. **Propose >= 2 candidate improvements** with expected gain from the diagnosis, choose one (or more) yourself, and implement behind a switch (default OFF; register in crates/dft-numerics/src/switches.rs, update the switch-count tests in crates/dft-numerics/tests/all/switches.rs and crates/dft-cli/tests/all/run.rs). Prefer new modules over editing hot shared files; keep diffs small so merges with other lanes stay easy.
4. **Validate on the Mac**: energy/density equivalence to the default (rounding-class changes: |dE| <= 1e-11 Ha unless the lane's contract says otherwise; trajectory changes need a full SCF and the strict 1e-10 / standard 5e-9 Ha gates), thread/rank-count independence where relevant, generality over cubic/skewed/slab cells, small/large band counts, real and synthetic PPs. No hard-coded constants tuned to one material.
5. **Write the dev08 request** `docs/performance/lanes/<LANE>/dev08-request.md`: pre-registered (a) the lane counter that must move and by how much, (b) exact A/B command lines per docs/performance/fast-loop-protocol.md (A/A floor first, crossed order-alternated pairs, one case per job, MPI cells si64 np8, si32 np8, si16 np8 and the one-core cells si32, si64 as relevant), (c) the success threshold (resolved gain |mean|/se >= 3), (d) what would make you withdraw it. Same parallel model on both sides; no OpenMP/threads work (user decision 2026-11-19).
   The paired mean/SE criterion is a resolution diagnostic; per-arm and final
   engine comparisons still report medians under the linked statistics rule.
6. Leave changes uncommitted, return the report (docs/schemas/report.schema.json) and stop. The orchestrator commits, merges, requests the dev08 A/B and relays the numbers.
7. If the orchestrator rejects (data attached), diagnose the miss, return with an alternative (new evidence or a different lever), or recommend closing the lane with a ledger row.

## Orchestrator duties
Relay requests to `baandi-dev08`, queue real work under the current job classes,
and stop keepalive when no real work remains. Measurements use socket 0 (np8:
cores 0-7; np16: cores 0-15); the socket-1 gate runs from `/scr/mino/gate` as
`final-a0-gate-<hash>` and is serialized with measurements. Adopt a resolved
same-model median gain only after multi-system full-convergence evidence,
strict 1e-10 Ha / standard 5e-9 Ha, independent verification and applicable
determinism gates; do not promote from a Si-only screen. Record every result
and failure in the ledger and update the dashboard within ownership.

> Historical, superseded runner instruction: Relay requests to dev08 (master session), keep one runner busy, decide adoption by the Adoption rule in AGENTS.md (adopt any resolved same-model gain, reject only for physical hard-coding or correctness failure), record every result in the ledger, update the dashboard.

## Common rules
Scratch inside the worktree (TMPDIR_990="$PWD/target/scratch"); 990 volume only; no QE source; no UPF commits; ledger row per attempt; defaults unchanged until approved. Known: 3 CLI pseudopotential-cache tests fail under the in-repo TMPDIR on the Mac (pass on dev08).
