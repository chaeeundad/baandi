TASK_ID: R1-DRIVER-SELECT
GOAL: |
  Stage R1 of docs/roadmap/refactor-r-plan.md: one driver selector with recorded reasons for every command; default still legacy. Fixes the precedence bug (finding 3: with REFERENCE_SCF=on and REFERENCE_STRICT=on the reference driver runs and STRICT is silently ignored) and the raw switch reads in config validation (finding 4).
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 0, 1a, 2, 4 "R1", 7 user decisions); crates/dft-cli/src/run.rs (options :123, reference_requested, reference_spin_fd_fallback*, execute_mode_inner header); config.rs (:854, :1227, :1312 raw REFERENCE_SCF reads); relax.rs (:398, :760); bands.rs; checkpoint_io.rs; inspect.rs; refscf_map.rs; the R0 routing corpus crates/dft-cli/tests/fixtures/routing/ and routing_corpus.rs; reference_identity.rs.
TASK: |
  1. crates/dft-cli/src/driver_select.rs: select(config, mode, command, upfs) -> Selection { driver: Reference | LegacyStandard | StrictOracle, reason: Option<&'static str>, source: Default | Env | Input }. Precedence: REFERENCE_STRICT=on -> StrictOracle; then explicit scf.driver input key (add "reference" | "legacy" if absent); then env REFERENCE_SCF; then default (legacy for now). Legacy-only predicates (plan section 2 rows L1-L11, each a declared predicate with a reason constant) are evaluated after UPF parsing, before engine construction; never catch-and-retry a reference error.
  2. Every command (run, relax, bands, dos, checkpoint, restart, explain, refscf-map) and config validation use the selection; no other code reads REFERENCE_SCF/REFERENCE_STRICT directly (add a test that greps the crate sources for direct reads outside driver_select.rs).
  3. results.json for every command gains scf.driver ("reference" | "legacy" | "strict_oracle"), scf.driver_reason, scf.driver_source; keep the existing reference_driver_fallback_reason field for now. One stderr line states the selection.
  4. crates/dft-cli/tests/driver_routing.rs drives the R0 corpus under: default; env on; input reference; input legacy; and STRICT combined with each; asserts driver and reason (update the corpus comments into machine-readable expectations if needed).
  5. Fix the stale mocks in scripts/tests/test_bench_scf.py (R0 found they fail on HEAD).
  Results must be bit-identical for all existing runs (reference_identity, baselines via bench_scf --require-identical for both baselines).
OWNED_FILES:
  - crates/dft-cli/src/driver_select.rs   # new
  - crates/dft-cli/src/lib.rs
  - crates/dft-cli/src/run.rs   # routing and results driver fields only
  - crates/dft-cli/src/config.rs   # routing reads and scf.driver key only
  - crates/dft-cli/src/relax.rs   # routing only
  - crates/dft-cli/src/bands.rs   # routing only
  - crates/dft-cli/src/dos.rs   # routing only
  - crates/dft-cli/src/checkpoint_io.rs   # option resolution only
  - crates/dft-cli/src/inspect.rs
  - crates/dft-cli/src/refscf_map.rs   # routing only
  - crates/dft-cli/tests/driver_routing.rs   # new
  - crates/dft-cli/tests/fixtures/routing/
  - crates/dft-cli/tests/routing_corpus.rs
  - scripts/tests/test_bench_scf.py
  - docs/cli.md   # scf.driver and routing paragraph only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-cli
  - python3 -I -m unittest discover -s scripts/tests -p test_bench_scf.py
  - cargo fmt --check && cargo clippy -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: single selector, STRICT wins, recorded driver/reason/source everywhere, routing matrix passes, results bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
