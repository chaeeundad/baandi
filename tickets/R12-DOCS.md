TASK_ID: R12-DOCS
GOAL: |
  Stage R12 of docs/roadmap/refactor-r-plan.md: documentation for the single production SCF path after R0-R11. Users and agents must find the current truth without reading history.
READ FIRST: docs/roadmap/refactor-r-plan.md (whole, incl. sections 7 decisions and 8 evidence), docs/roadmap/reference-generalization.md, docs/design/reference-scf.md, docs/INDEX.md, docs/cli.md, docs/user/*.md, docs/checkpoint-format.md, docs/performance/README.md, docs/performance/switch-retirement.md, crates/*/AGENTS.md, AGENTS.md (do not edit it), README.md (read only), docs/decisions/ (0043, 0044, 0045, 0046).
TASK: |
  1. New ADR docs/decisions/0047-reference-driver-single-path.md (status accepted, user decisions D1-D6 of 2026-11-25): routing (Reference | StrictOracle, G4D to the oracle), the oracle feature and its freeze rule, removed keys (D4), switch rules (D3: rollback switches removed, recreate from git), checkpoint v2, runtime/stepper, evidence table (gen14-gen28 from section 8 of the plan and RESULT paths), known open items (G4D free-moment PBE; restart writes an export container; spin/FD results lack map-boundary clocks; preflight vs results start-fallback reason mismatch noted by dev08 gen27).
  2. Update: docs/INDEX.md; docs/cli.md (drivers, removed keys, checkpoint/restart v2 commands and outputs, observe/agent on reference); docs/user/*.md where they describe legacy-only keys or the legacy driver; docs/checkpoint-format.md (v2 primary, v1 legacy/oracle only); docs/performance/README.md (baselines); crates/dft-engine/AGENTS.md and crates/dft-cli/AGENTS.md (remove or mark historical sections that describe deleted code: say "removed in R10/R11, see git history"); docs/design/reference-scf.md (a top "Status" section: the reference driver is the only production path); docs/roadmap/reference-generalization.md (row R done with evidence).
  3. README.md: do NOT edit it. Write docs/roadmap/readme-corrections.md listing every statement outside the generated block that is no longer true (line, current text, proposed text, reason) for the user to approve (decision D6). The generated block must stay fresh (make_readme_figures + its test).
  4. Mark docs/roadmap/refactor-r-plan.md status "completed except G4D" with links.
OWNED_FILES:
  - docs/
  - crates/dft-engine/AGENTS.md
  - crates/dft-cli/AGENTS.md
  - crates/baandi-py/AGENTS.md
  - crates/baandi-mcp/   # docs/README only, if present
TEST_COMMANDS:
  - python3 -I scripts/make_readme_figures.py && python3 -I scripts/tests/test_readme_figures.py
  - git diff --check
ACCEPTANCE_CRITERIA: ADR 0047 accepted; docs consistent with the code; README correction list ready for the user.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
