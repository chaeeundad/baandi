TASK_ID: I26-3-TESTFIX
GOAL: |
  Read `gh issue view 26` is not possible in the Codex sandbox (no network), so the issue text relevant to you is quoted below. Issue #26 (user-filed 2026-11-24, reviewed on main 2f91471a): "Align SCF acceptance tests, profiling and benchmark instructions with the current design". Gates for every item: retain independent fixed-point verification and charge failed attempts and recovery to delivered wall; keep the strict 1e-10 Ha and standard 5e-9 Ha classes; multi-system evidence before any default promotion; the independent-source policy (no QE source). Defaults unchanged. Scratch inside the worktree. Do not edit scripts/delegate.sh. Leave changes uncommitted; return docs/schemas/report.schema.json.
  Item: "Resolve the cache-test and README failures at their source. The R-SCF-2 full gate recorded 825 passes and three cache-test failures with scratch inside the checkout. Give those tests a private writable temp directory on the authorized external volume outside the checkout. Retain the test that rejects unsafe cache paths. Separately, replace the README test's obsolete 'standard profile is unsupported' wording assertion (scripts/tests/test_readme_figures.py about lines 49-58) with checks against the current supported operations and their evidence. Report the three cache-test, README-test and full integration-gate results separately."
  Authorized external volume: /Volumes/990-2TB (the user rule: builds and scratch only there). The Codex sandbox may not be able to write outside the worktree; design the tests so the temp root is taken from an environment variable (for example BAANDI_TEST_TMP) with a safe default, document how the gate sets it (scripts/check.sh already exports TMPDIR_990 for scratch), and verify inside the sandbox with a directory under the worktree that is OUTSIDE the cargo package directories if the volume is not writable from the sandbox; say exactly what you could and could not verify.
OWNED_FILES:
  - crates/dft-cli/tests/   # the three cache tests and helpers only
  - scripts/tests/test_readme_figures.py
  - scripts/check.sh   # only the temp-root export for these tests
TEST_COMMANDS:
  - scripts/check.sh full
  - python3 -I scripts/tests/test_readme_figures.py
ACCEPTANCE_CRITERIA: the three cache tests pass with the external temp root, unsafe-path rejection test kept, README test checks current supported operations, results reported separately.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
