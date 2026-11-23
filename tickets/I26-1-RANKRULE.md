TASK_ID: I26-1-RANKRULE
GOAL: |
  Read `gh issue view 26` is not possible in the Codex sandbox (no network), so the issue text relevant to you is quoted below. Issue #26 (user-filed 2026-11-24, reviewed on main 2f91471a): "Align SCF acceptance tests, profiling and benchmark instructions with the current design". Gates for every item: retain independent fixed-point verification and charge failed attempts and recovery to delivered wall; keep the strict 1e-10 Ha and standard 5e-9 Ha classes; multi-system evidence before any default promotion; the independent-source policy (no QE source). Defaults unchanged. Scratch inside the worktree. Do not edit scripts/delegate.sh. Leave changes uncommitted; return docs/schemas/report.schema.json.
  Item: "Enforce the accepted cross-rank rounding rule in the release test. On Si16, one versus three simulated ranks straddle a Davidson stopping threshold by a Ritz-change difference of 9.16e-16 Ha, yielding 93 versus 92 H columns. All 24 local case/rank combinations passed the standard energy class and independent verification. The test still requires identical iteration and H-column counts and a blanket 1e-9 e/Bohr^3 density cap. Require agreement up to the first divergence, record both Ritz values and threshold margins, and apply the documented exception only when the rounding evidence and delivered-state checks pass. Include a negative test for an unexplained divergence away from a threshold. Preserve fixed-layout repeatability and publish the release-matrix exit code." The accepted rule is in docs/design/reference-scf.md section "Rank-count acceptance for the reference driver"; the assertions are in crates/dft-engine/tests/reference_driver.rs (about lines 84-155).
TASK: implement exactly this in the test (and any small record fields the driver/solver must expose to identify the divergence point, the two Ritz-change values, the threshold and the margin; no behavior change), add the negative test, keep a fixed-layout repeatability test (same rank count twice: bit-identical), run the release matrix and report its exit code and the table (case, ranks, first divergence map/iteration, margins, delivered-state checks).
OWNED_FILES:
  - crates/dft-engine/tests/reference_driver.rs
  - crates/dft-engine/tests/common/
  - crates/dft-engine/src/refscf/   # record fields only
  - crates/dft-solvers/src/ref_davidson.rs   # record fields only
  - docs/design/reference-scf.md   # rank-count section only
TEST_COMMANDS:
  - scripts/check.sh full
ACCEPTANCE_CRITERIA: rule enforced, negative test present, matrix exit code published, defaults bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
