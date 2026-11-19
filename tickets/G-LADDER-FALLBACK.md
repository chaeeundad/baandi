TASK_ID: G-LADDER-FALLBACK
GOAL: |
  dev08 suite survey (default run, seed 601): on nine synthetic cells (syn-cubic-degenerate, syn-d3, syn-fft-primes, syn-odd-shift, syn-p3, syn-s1, syn-s2, syn-spd, syn-z1) the cutoff ladder is not used; results.json carries fallback_reason = coarse stage: Numerical("k index 0: InvalidInput(\"bands/subspace incompatible with k basis\")"). The synthetic cell syn-gaussian-fcc runs the ladder normally. Find out whether this is a defect (the coarse basis at the ladder cutoff is smaller than the band/subspace count, so the stage errors and the run silently falls back) or a deliberate limit, and why it surfaces as an error string rather than a planned decision.
TASK: (1) Reproduce on the Mac in release (synthetic fixtures live in the test/benchmark fixtures; find them via docs/performance/benchmark-suite.md and tests/). (2) Explain the cause with numbers (coarse plane-wave count vs bands, Davidson block size/subspace size). (3) Propose the general fix behind a default-off switch if it is a defect: e.g. decide before running the coarse stage (plan-time check that the coarse basis can hold bands and the Davidson subspace) and skip the ladder with a clear structured reason, or enlarge the coarse cutoff to the minimum that fits; do not change results when the ladder is not used. Report which of the nine cells would then run a ladder and the Mac total-H effect vs no ladder (3 seeds). Keep |dE| classes strict 1e-10 / standard 5e-9. (4) Write docs/design/ladder-fallback.md and a ledger row. Defaults unchanged; scratch inside worktree (TMPDIR_990="$PWD/target/scratch").
OWNED_FILES:
  - crates/dft-engine/src/cutoff_ladder.rs
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/ladder-fallback.md
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/design/ladder-fraction.md, docs/performance/benchmark-suite.md
TEST_COMMANDS:
  - scripts/check.sh fast dft-engine dft-numerics dft-cli   # 3 CLI cache tests reject in-repo TMPDIR, pass on dev08
ACCEPTANCE_CRITERIA: cause explained with numbers; defect vs design stated; if defect, default-off fix with Mac evidence on the nine cells x 3 seeds; defaults bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
