TASK_ID: G-GUARD-GRAM
GOAL: |
  Implement section 3 of docs/performance/plan-2026-11-19.md (read it first). On the owned M2 route `GridPlan::stop_error` (crates/dft-engine/src/grid_dist.rs ~L1314) already reads both distributed density spectra and forms their difference, but when BAANDI_MIX_COULOMB_GUARD is enabled `DistributedPulay::coulomb_guard` (crates/dft-engine/src/mix_dist.rs ~L195) builds a real trial residual and its own plan transforms 2 fields per check, plus m+1 more per refit (m retained differences), each moving data between ranks. Replace those transforms by Gram algebra: export an epoch-tagged read-only residual spectrum r_hat = rho_out_hat - rho_in_hat on eligible maps; keep the last m+1 ACTUAL residual spectra on the owning ranks; maintain their small Coulomb Gram K_ij = <r_i, r_j>_C with the exact current weights, normalization and canonical fixed-order reductions (one new row + diagonal per new residual); the trial p = r - sum_j gamma_j Delta r_j has coefficient vector a over the raw residual history so <p,p>_C = a^T K a, and the difference Gram / right-hand side for a refit follow from K by row/column differences. This changes only HOW the existing guard is evaluated, no new damping rule. Switch: new BAANDI_M2_GUARD_GRAM=on|off, default OFF (add to crates/dft-numerics/src/switches.rs and update the switch-count tests crates/dft-numerics/tests/all/switches.rs and crates/dft-cli/tests/all/run.rs). Keep the direct path (current transforms) on clamp, stale epoch, unsupported layout, history reseed, an ambiguous decision at the rounding boundary, and recompute near-zero quadratic forms directly when cancellation dominates (state the threshold derivation). The plan's independent algebra check (7 synthetic residuals on an odd 7x9x5 periodic grid, integer-mode Coulomb weights, G=0 excluded, 100 coefficient vectors) matched explicit transforms to 6.9e-16: reproduce it as a unit test, then test production normalization, Nyquist modes, rank-count independent bits and identical decisions on odd/even, skewed and ill-conditioned histories (simulated ranks 1..7), and equality of the guard decisions (check/refit flags, coefficients) with the direct path. Counters (zero cost when timers are off): guard transforms avoided, plan builds, reductions, extra history bytes, check/refit counts. Mac: counts and identity only; real timing (arms A: guard off, B: guard on / GUARD_GRAM off, C: guard on / GUARD_GRAM on, same executable) is run by the orchestrator on dev08: write the dev08 command lines in your report. Credit ZERO seconds against default A until A-to-C wins on dev08; the guard itself stays default off.
OWNED_FILES:
  - crates/dft-engine/src/mix_dist.rs
  - crates/dft-engine/src/grid_dist.rs
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-engine/tests/guard_gram.rs   # new
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/performance/plan-2026-11-19.md, crates/dft-engine/AGENTS.md, ledger rows MIX-COULOMB-GUARD*, SCF-TRAJECTORY, docs/performance/requests/mix-coulomb-guard-dist.md
IMPLEMENTATION_STEPS:
  1. Read the guard/stop code; write the algebra unit test first; then the Gram history and the switch; then the equality tests.
  2. scripts/check.sh fast dft-engine dft-numerics dft-cli with TMPDIR_990="$PWD/target/scratch" (scratch inside your worktree; never /tmp).
TEST_COMMANDS:
  - scripts/check.sh fast dft-engine dft-numerics dft-cli
ACCEPTANCE_CRITERIA: default off bit-identical; GUARD_GRAM on gives identical guard decisions and coefficients to rounding class 1e-12 relative, rank-count independent bits; counters show the avoided transforms; clippy/fmt clean.
NON_GOALS: enabling the guard by default; new damping rules.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
