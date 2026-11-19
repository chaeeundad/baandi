# Lane charter (user decision 2026-11-20)

The SCF computation is split into independent units ("lanes"). Each lane has its own Codex agent that diagnoses, designs and implements improvements on its own judgment. dev08 tests each result; the orchestrator (Claude) approves or rejects. A rejected lane agent brings an alternative.

## Agent duties (every lane)
1. Read AGENTS.md, docs/performance/AGENT-BRIEF.md, attempts-ledger.md (grep your lane's terms) and backlog.md. Never repeat a rejected idea without stating the new evidence.
2. **Diagnose first** (Mac release, small and medium cells; instrument with PHASE_PROFILE/STAGE_TIMERS or your own timers): where does your lane's time go (per call, per map, per rank), what is the floor (arithmetic count, bandwidth, latency)? Write the finding as `docs/performance/lanes/<LANE>/diagnosis.md` with numbers.
3. **Propose >= 2 candidate improvements** with expected gain from the diagnosis, choose one (or more) yourself, and implement behind a switch (default OFF; register in crates/dft-numerics/src/switches.rs, update the switch-count tests in crates/dft-numerics/tests/all/switches.rs and crates/dft-cli/tests/all/run.rs). Prefer new modules over editing hot shared files; keep diffs small so merges with other lanes stay easy.
4. **Validate on the Mac**: energy/density equivalence to the default (rounding-class changes: |dE| <= 1e-11 Ha unless the lane's contract says otherwise; trajectory changes need a full SCF and the strict 1e-10 / standard 5e-9 Ha gates), thread/rank-count independence where relevant, generality over cubic/skewed/slab cells, small/large band counts, real and synthetic PPs. No hard-coded constants tuned to one material.
5. **Write the dev08 request** `docs/performance/lanes/<LANE>/dev08-request.md`: pre-registered (a) the lane counter that must move and by how much, (b) exact A/B command lines per docs/performance/fast-loop-protocol.md (A/A floor first, crossed order-alternated pairs, one case per job, MPI cells si64 np8, si32 np8, si16 np8 and the one-core cells si32, si64 as relevant), (c) the success threshold (resolved gain |mean|/se >= 3), (d) what would make you withdraw it. Same parallel model on both sides; no OpenMP/threads work (user decision 2026-11-19).
6. Leave changes uncommitted, return the report (docs/schemas/report.schema.json) and stop. The orchestrator commits, merges, requests the dev08 A/B and relays the numbers.
7. If the orchestrator rejects (data attached), diagnose the miss, return with an alternative (new evidence or a different lever), or recommend closing the lane with a ledger row.

## Orchestrator duties
Relay requests to dev08 (master session), keep one runner busy, decide adoption by the Adoption rule in AGENTS.md (adopt any resolved same-model gain, reject only for physical hard-coding or correctness failure), record every result in the ledger, update the dashboard.

## Common rules
Scratch inside the worktree (TMPDIR_990="$PWD/target/scratch"); 990 volume only; no QE source; no UPF commits; ledger row per attempt; defaults unchanged until approved. Known: 3 CLI pseudopotential-cache tests fail under the in-repo TMPDIR on the Mac (pass on dev08).
