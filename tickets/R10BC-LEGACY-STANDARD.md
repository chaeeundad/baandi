TASK_ID: R10BC-LEGACY-STANDARD
GOAL: |
  Stages R10b and R10c of docs/roadmap/refactor-r-plan.md. After R6a/R6b/R6c/R7/R8 the only remaining fallback is L1 (free-moment collinear FD with PBE, G4D: negative qualification result) and L9 (legacy-only SCF input keys). Architect decision: route L1 to the strict oracle (driver "strict_oracle", reason G4D, recorded) instead of the legacy standard path; per user decision D4 the legacy-only input keys become ordinary unknown-key errors (no migration message). Then the legacy *standard* path is dead: delete it, keeping only the strict oracle inside src/oracle.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 2 rows L1/L9, 4 "R10b"/"R10c", 7 decisions); crates/dft-cli/src/{driver_select.rs,config.rs,run.rs,finite_run.rs,spin_run.rs,relax.rs,bands.rs,checkpoint_io.rs}; crates/dft-engine/src/oracle/ (session.rs/scf.rs standard-profile branches: stall guard, density forcing, stop predictor, verify schedule, cert screen, stop metric, standard verification attempts); crates/dft-numerics/src/switches.rs (REFERENCE_SCF and the oracle-internal switches); docs/user/*.md and docs/cli.md (input keys).
TASK: |
  1. Routing: driver_select becomes Reference | StrictOracle. L1 inputs select StrictOracle with reason G4D_FALLBACK (record it); remove LegacyStandard, the REFERENCE_SCF switch and the scf.driver = "legacy" value (input error: unknown value). scf.driver = "reference" stays accepted as a no-op.
  2. D4: remove the legacy-only SCF keys listed in L9 (xc.density_floor, scf.mixing linear/kerker-pulay, mixing_q0_inverse_bohr, mixing_metric, spectral_guard, locking, preconditioner, davidson_capacity_per_band, start_density uniform, start_wavefunctions non-atomic, energy_tolerance_ha, density_residual_tolerance, eigen_tolerance, standard_stall_guard, standard_density_forcing, cutoff_ladder "on", cutoff_ladder_fraction) from the config schema, so they fail as unknown keys. If the strict oracle needs any of them internally, hard-code its former default inside oracle/ (the oracle's arithmetic must not change: RG-ORC bit identity).
  3. Delete the legacy standard-path code in oracle/ (all branches only reachable with a standard profile) and the CLI legacy-standard branches (finite_run.rs/spin_run.rs standard routes, relax/bands/checkpoint legacy branches), plus tests that exercise only the legacy standard path (list each). Tests that used driver="legacy" to obtain an oracle comparison switch to BAANDI_REFERENCE_STRICT=on.
  4. Identity: reference baselines and oracle baselines bit-identical (bench_scf --require-identical both), reference_identity, oracle_identity; the full check; report before/after lines, tests, binaries, switches.
OWNED_FILES:
  - crates/dft-cli/src/
  - crates/dft-cli/tests/
  - crates/dft-engine/src/oracle/
  - crates/dft-engine/src/lib.rs   # re-exports only
  - crates/dft-engine/tests/   # deletions/retargeting only
  - crates/dft-numerics/src/switches.rs   # REFERENCE_SCF and orphaned switches only
  - crates/dft-numerics/tests/all/switches.rs
  - crates/baandi-py/   # remove the legacy driver option only
  - crates/baandi-mcp/   # only if it references the legacy driver
  - docs/cli.md   # driver and removed keys
  - docs/user/   # removed keys only
  - docs/performance/attempts-ledger.md   # rows for removed switches
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo check -p dft-engine --no-default-features
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (reference baseline; and BAANDI_REFERENCE_STRICT=on with baseline-oracle.json)
  - python3 -I scripts/make_readme_figures.py && python3 -I scripts/tests/test_readme_figures.py
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: only Reference and StrictOracle remain; legacy standard path and legacy-only keys deleted; both baselines bit-identical; counts reported.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
