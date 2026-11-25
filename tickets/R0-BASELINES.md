TASK_ID: R0-BASELINES
GOAL: |
  Stage R0 of docs/roadmap/refactor-r-plan.md: freeze identity baselines and a routing corpus before any refactoring. No production code change.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 0, 2, 4 "Gate definitions" and "R0"); scripts/bench_scf.py; docs/performance/baseline.json and docs/performance/README.md; scripts/gen_synthetic_upf.py; crates/dft-cli/src/run.rs routing (options, reference_requested, reference_spin_fd_fallback*); crates/dft-cli/src/config.rs; the legacy-only input keys and UPF requirements listed in plan section 2 (L1-L11) and finding 7.
TASK: |
  1. scripts/bench_scf.py: record scf.driver (and the reference/strict switches in effect) in every baseline and refuse to compare against a baseline whose driver differs. Record docs/performance/baseline-reference.json (BAANDI_REFERENCE_SCF=on) and docs/performance/baseline-oracle.json (BAANDI_REFERENCE_STRICT=on, REFERENCE_SCF off) for the existing bench cases, release build on this Mac (these are identity baselines: energies, density hashes, counts; timings are observations only).
  2. crates/dft-cli/tests/reference_identity.rs: tiny Si2 Gamma LDA with a synthetic UPF (PP_RHOATOM and PP_PSWFC present), in-process, under 10 s in debug; runs the reference path and the strict oracle; asserts density sha256 and energy bits keyed by target_os/target_arch (skip with a note on other targets). Record the current values.
  3. crates/dft-cli/tests/fixtures/routing/*.toml: one minimal input per legacy-only row L1-L11 of plan section 2 (plus plain Gamma LDA, k mesh, FD, fixed spin, symmetry, PBE+NLCC, relax ions that must route to reference), each with a leading comment "# expected_driver = ..." and "# expected_reason = ..." describing TODAY's behaviour with REFERENCE_SCF=on (reference or which rejection/fallback). Add a test that parses every fixture (config validation only, no SCF) so the corpus cannot rot.
  4. Document in docs/performance/README.md (baseline paragraph only) which baseline is which.
OWNED_FILES:
  - scripts/bench_scf.py
  - docs/performance/baseline-reference.json   # new
  - docs/performance/baseline-oracle.json   # new
  - docs/performance/README.md   # baseline paragraph only
  - crates/dft-cli/tests/reference_identity.rs   # new
  - crates/dft-cli/tests/routing_corpus.rs   # new
  - crates/dft-cli/tests/fixtures/routing/   # new
TEST_COMMANDS:
  - cargo test -p dft-cli --test reference_identity
  - cargo test -p dft-cli --test routing_corpus
  - python3 -I scripts/bench_scf.py --help
  - cargo fmt --check && cargo clippy -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: two baselines recorded with driver metadata; identity test and routing corpus pass; no production code changed.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
