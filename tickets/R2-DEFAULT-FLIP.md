TASK_ID: R2-DEFAULT-FLIP
GOAL: |
  Stage R2 of docs/roadmap/refactor-r-plan.md: the reference driver becomes the default. R1 (driver_select.rs) is merged; legacy-only inputs fall back with recorded reasons; STRICT selects the oracle.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 0, 4 "R2", 5 risks, 7 decisions); crates/dft-cli/src/driver_select.rs; crates/dft-numerics/src/switches.rs REFERENCE_SCF; scripts/bench_scf.py and docs/performance/baseline*.json; the CLI tests.
TASK: |
  1. Flip the REFERENCE_SCF default to "on" in the switch table (meaning text: "production driver; off = legacy standard path, deprecated, removed in R10"); update legacy_defaults_are_unchanged and the switch-count/default tests.
  2. Run the whole workspace test suite. For every CLI or engine test whose expectation is legacy-specific: if it tests a legacy feature, pin it with [scf] driver = "legacy" (or the env) in its input; if it tests generic behaviour, re-baseline the expectation to the reference result and state why the new value is correct (standard class against the oracle). List every touched test in the report under "pinned" and "re-baselined" with one line each.
  3. bench_scf.py: default baseline becomes docs/performance/baseline-reference.json; rename docs/performance/baseline.json to baseline-legacy.json (git mv is done by the orchestrator: just update references and note it); update docs/performance/README.md baseline paragraph.
  4. Proof that the flip changes only the selection: the default run (no env) is bit-identical to the R0 reference baseline (bench_scf --require-identical against baseline-reference.json with no switches set), and reference_identity passes unchanged; the STRICT oracle baseline also unchanged.
OWNED_FILES:
  - crates/dft-numerics/src/switches.rs   # REFERENCE_SCF entry only
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/   # pins and re-baselines only
  - crates/dft-engine/tests/   # pins only if an engine test reads the CLI default (unlikely)
  - crates/dft-cli/src/driver_select.rs   # default only
  - scripts/bench_scf.py
  - scripts/tests/test_bench_scf.py
  - docs/performance/README.md   # baseline paragraph only
  - docs/cli.md   # driver default paragraph only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - scripts/check.sh   (BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical --baseline docs/performance/baseline-reference.json (release, no switches set)
  - python3 -I -m unittest discover -s scripts/tests -p test_bench_scf.py
ACCEPTANCE_CRITERIA: default is reference; full suite green with an explicit pinned/re-baselined list; default run bit-identical to the R0 reference baseline; oracle unchanged.
GIT: leave changes uncommitted (the orchestrator renames baseline.json).
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
