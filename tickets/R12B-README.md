TASK_ID: R12B-README
GOAL: |
  The user approved (decision D6, option A, 2026-11-26) all README corrections listed in docs/roadmap/readme-corrections.md. Apply them exactly as proposed.
READ FIRST: docs/roadmap/readme-corrections.md (all entries, including "Generated validation block: provenance follow-up" and "Additional review boundaries"); README.md; scripts/make_readme_figures.py; docs/feature-status.json; scripts/tests/test_readme_figures.py.
TASK: |
  1. README prose outside the generated blocks: apply each proposed text (lines 35, 49, 56, 60, 64, 66, 70) verbatim or with the minimal wording the entry specifies.
  2. Status table: edit docs/feature-status.json for every "Generated status source corrections" entry, then regenerate with python3 -I scripts/make_readme_figures.py; never hand-edit generated text.
  3. Validation block: follow the "provenance follow-up" section: label the drivers of the stored QE/checkpoint/runtime evidence in the generator if the entry asks for it, and add the reference v2 checkpoint and runtime evidence files (tests/reference/checkpoint/reference-v2-results.json, tests/reference/runtime/generated/reference-results.json) as separate rows/paragraphs without replacing stored values. Respect the "Additional review boundaries" (no new speed ratios, keep stated limits, mention restart writes an export container).
  4. Mark each entry in docs/roadmap/readme-corrections.md as "applied (D6 approved 2026-11-26)".
OWNED_FILES:
  - README.md
  - docs/feature-status.json
  - docs/img/   # regenerated figures only
  - scripts/make_readme_figures.py   # driver labels / new evidence rows only
  - scripts/tests/test_readme_figures.py   # only if the generator contract changes
  - docs/roadmap/readme-corrections.md
TEST_COMMANDS:
  - python3 -I scripts/make_readme_figures.py && python3 -I scripts/tests/test_readme_figures.py
  - git diff --check
ACCEPTANCE_CRITERIA: all approved corrections applied; generated blocks regenerated and fresh; no stored numbers replaced.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
