TASK_ID: G-SCHED-OVERHEAD
GOAL: |
  BAANDI_VERIFY_SCHEDULE is default on (4c3d09d: dev08 one-core si64 -6.81 %), but the one-core si32 control (no failed attempts, bit-identical energy and counts, prediction_evaluations 0) is +0.29 % slower with the switch on (8 of 8 pairs, |mean|/se 6.4): about 23 ms over 10 maps = about 2.3 ms per map, which is NOT prediction arithmetic (counter is 0) but bookkeeping on the schedule path: allocation of a record structure per map, a branch that does extra work every map, cloning of history/spectra, string/JSON building, timers, or a different code path in the SCF loop when the switch is on (crates/dft-engine/src/session.rs, verification schedule code added in F-VERIFY-SCHEDULE; see docs/design/verification-schedule.md). TASK: find and remove it. (1) Micro-profile the SCF loop at si32 and si16 one core on the Mac in release mode with `BAANDI_VERIFY_SCHEDULE=off` vs `on` (use alternating pairs, report the per-map time difference with the pairwise-mean/standard-error method of docs/performance/fast-loop-protocol.md; Mac noise is large, so also use instruction-level evidence: read the code path and count allocations/clones per map with a counting allocator in a test, or dhat-style counters if available offline); (2) list every operation the schedule adds per ordinary map when no attempt failed, and make it zero: no allocation, no clone, no formatting, no extra history copies; keep the data needed after a failed attempt in structures that already exist or are created lazily on the first failure; (3) results must stay bit-identical to the current `on` and to `off` when no attempt fails (tests: the existing schedule tests plus a new test asserting that with no failed attempt the allocation count of the schedule path per map is 0); (4) keep behavior on failed attempts exactly as is (one-core si64 attempts [11,13], trace record). Scratch inside your worktree (TMPDIR_990=$PWD/target/scratch). Report the per-map operations found and removed and Mac numbers; the orchestrator re-measures on dev08 (one-core si32 8 pairs, one-core si64 crossed A/B).
OWNED_FILES:
  - crates/dft-engine/src/session.rs   # schedule bookkeeping only
  - crates/dft-engine/src/verification_tests.rs
  - crates/dft-engine/tests/verify_schedule_alloc.rs   # new
  - docs/design/verification-schedule.md
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/design/verification-schedule.md, ledger rows F-VERIFY-SCHEDULE*
IMPLEMENTATION_STEPS:
  1. Read the schedule code path executed per ordinary map and list its operations.
  2. Remove/lazy them; add the allocation-count test; run scripts/check.sh fast dft-engine dft-numerics dft-cli.
TEST_COMMANDS:
  - scripts/check.sh fast dft-engine dft-numerics dft-cli
ACCEPTANCE_CRITERIA: zero per-map allocations on the schedule path when no attempt fails; bit-identical results; all tests green (known: 3 CLI cache tests reject in-repo TMPDIR and pass on dev08).
NON_GOALS: changing the schedule policy or the certificate.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
