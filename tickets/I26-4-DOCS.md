TASK_ID: I26-4-DOCS
GOAL: |
  Read `gh issue view 26` is not possible in the Codex sandbox (no network), so the issue text relevant to you is quoted below. Issue #26 (user-filed 2026-11-24, reviewed on main 2f91471a): "Align SCF acceptance tests, profiling and benchmark instructions with the current design". Gates for every item: retain independent fixed-point verification and charge failed attempts and recovery to delivered wall; keep the strict 1e-10 Ha and standard 5e-9 Ha classes; multi-system evidence before any default promotion; the independent-source policy (no QE source). Defaults unchanged. Scratch inside the worktree. Do not edit scripts/delegate.sh. Leave changes uncommitted; return docs/schemas/report.schema.json.
  Item: "Reconcile the live instructions for exploration, promotion and timing. Point the agent brief, lane charter and final-table instructions to the applicable protocol. Use fixed-H replay for kernels, common-input replay for a density-update discrepancy, and full convergence for a tail or verification change. An existing A/B showed a 1.2% gain at ten steps and 22% at full length (docs/performance/fast-loop-protocol.md about lines 154-158). Apply the same median statistic to both engines. Keep the 600 s fast-job limit and identify the authorized class for longer diagnostics. Preserve historical requests and raw evidence." Old final-table instruction: docs/performance/final-measurement-protocol.md (about lines 10-18); current comparison rule: docs/performance/comparison-rules.md (about lines 3-8).
  Also reflect operational facts of 2026-11-20..24: dev08 is operated by the baandi-dev08 session (the dft5090 relay is retired), the gate runs on socket 1 from /scr/mino/gate as queue job final-a0-gate-<hash>, keepalive is stopped when there is no real work, measurements on socket 0 (np8 cores 0-7, np16 cores 0-15), the threads lane is dropped, the new reference SCF path (issue 25) is the active work.
TASK: edit the live instruction documents (docs/performance/AGENT-BRIEF.md, docs/performance/lanes/CHARTER.md, docs/performance/final-measurement-protocol.md, docs/performance/comparison-rules.md, docs/performance/fast-loop-protocol.md, docs/performance/README.md if it instructs timing) so they agree with each other and point to one protocol per question type (kernel: fixed-H replay; density update: common-input replay; tail/verification/trajectory change: full convergence with verification charged), one median statistic for both engines, the 600 s fast-job class and the named longer class (final- jobs up to 7200 s), and a short "current state" header. Mark superseded paragraphs as historical instead of deleting them; keep requests and raw evidence untouched. No code changes.
OWNED_FILES:
  - docs/performance/AGENT-BRIEF.md
  - docs/performance/lanes/CHARTER.md
  - docs/performance/final-measurement-protocol.md
  - docs/performance/comparison-rules.md
  - docs/performance/fast-loop-protocol.md
  - docs/performance/README.md
TEST_COMMANDS:
  - git diff --check
ACCEPTANCE_CRITERIA: one consistent protocol map across the live documents, historical text marked not deleted.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
