TASK_ID: G-HANDOFF-OWNED
GOAL: |
  Implement section 5 of docs/performance/plan-2026-11-19.md (read it first). The cutoff ladder hands the coarse result to the fine stage through full matrices: crates/dft-engine/src/cutoff_ladder.rs `warm_state()`/terminal materialization of the lazy block (~L539), `CoarseState` carries a FULL complex orbital matrix, `prolong_orbitals` (~L339) allocates a fine-basis matrix of zeros and copies by Miller index, and the new fine stage consumes it (~L414 constructor). That is mathematically an injection of a smaller index set into a larger one and does not require every rank to hold every coefficient. STEP 1 (price before building): add timers/counters (profile off: zero cost) that give the exclusive cost of just the handoff on the production path: materialization gather, prolongation copy, import check, session setup that consumes the state, bytes per rank, collectives. Run Mac simulated-rank cases (m2_entry_cost with M2_PROFILE_MODE=production M2_PROFILE_LADDER=on, fast_screen_mpi.py) and report the handoff share; the real number comes from dev08 (the orchestrator measures with BAANDI_PHASE_PROFILE=on at si64 np=8). The ceiling of any gain is the measured exclusive handoff cost, NEVER the whole coarse stage (1.2 s at si64 np=8). STEP 2 (only if the handoff share is >= 2 % of the run on the simulated-rank harness or the counters show an avoidable gather): prototype a direct coarse-owner -> fine-owner transfer by integer Miller index behind BAANDI_M2_OWNED_HANDOFF=on|off (default OFF; add to crates/dft-numerics/src/switches.rs, update switch-count tests crates/dft-numerics/tests/all/switches.rs and crates/dft-cli/tests/all/run.rs): precompute the permutation once for this pair of layouts, send only populated rows, fill remaining owned fine rows with zero; check conjugate partners, Gamma packing factors, the zero mode, differing ownership, and identical complex-to-real conventions; keep the full-matrix route as the oracle and for unsupported consumers; preserve import norm/orthogonality checks (LADDER-IMPORT-SPLIT), checkpoint behaviour, potential identity and the error agreement across ranks (a rank-local mismatch must become an agreed error before any collective diverges). A pure copy/permutation in identical storage can be bit-identical; repacking or a different reduction order needs rounding-class validation. Record peak bytes/rank, bytes transmitted, collective calls, gather/materialization count and net handoff seconds. Tests on simulated ranks 1..7 incl. skewed cells, odd/even grids, Gamma and the zero mode, fallbacks, bit identity of the fine-stage result vs the full-matrix route. Do NOT work on the asynchronous fork/join overlap extension (plan section 5 last paragraph): note in the ledger that it needs its own screen.
OWNED_FILES:
  - crates/dft-engine/src/cutoff_ladder.rs
  - crates/dft-engine/src/session.rs   # handoff consumption/timers only
  - crates/dft-engine/src/m2.rs   # counters only
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-engine/tests/owned_handoff.rs   # new
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/performance/plan-2026-11-19.md, crates/dft-engine/AGENTS.md, ledger rows LADDER-IMPORT-SPLIT, LAZY-ASSEMBLE-FUSED, BETWEEN-MATERIALIZE, M2-OWNED-*
IMPLEMENTATION_STEPS:
  1. Add the handoff timers and price it; write the cost in a ledger row with the Mac numbers and the dev08 command.
  2. If the go condition holds, implement the direct transfer behind the switch with the tests above.
  3. scripts/check.sh fast dft-engine dft-numerics dft-cli with TMPDIR_990="$PWD/target/scratch" (never /tmp).
TEST_COMMANDS:
  - scripts/check.sh fast dft-engine dft-numerics dft-cli
ACCEPTANCE_CRITERIA: default off bit-identical; timers show the exclusive handoff cost; if implemented: fine-stage result identical (bits or rounding class stated), error agreement tested, counters show fewer gathered bytes; clippy/fmt clean.
NON_GOALS: async fork/join of coarse and fine preparation, removing the ladder.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
