TASK_ID: R6A-G4B-DEFAULT
GOAL: |
  Stage R6a of docs/roadmap/refactor-r-plan.md: collinear Fermi-Dirac with fixed or free moment runs on the reference driver by default. Turn REFSCF_SPIN_FD_MOMENT, REFSCF_SPIN_STIFFNESS and REFSCF_HISTORY_CONDITIONING on by default and remove the L2 routing predicate. Evidence: G4B synthetic acceptance 48/48 with the three switches (docs/design/reference-scf.md "Spin and Fermi-Dirac"); dev08 gen15 bcc Fe (PBE+NLCC, 4x4x4, fixed moment 4.4) VERIFIED in 75.8 s on the reference driver while legacy and the strict oracle did not finish in 45 min. Free-moment FD with PBE stays routed to legacy (L1, G4D).
READ FIRST: crates/dft-cli/src/driver_select.rs (l1, l2, spin_fd, reference_mixing); crates/dft-numerics/src/switches.rs (the three entries); docs/decisions/0045-broyden-history-conditioning.md (proposed); tests that pin these switches (grep the three switch names across crates).
TASK: |
  1. Flip the three defaults to on (meaning text cites the evidence); remove l2 and G4B_FALLBACK handling where it becomes dead (keep G4D); update routing fixtures/expectations.
  2. Mark ADR 0045 accepted with the evidence above.
  3. Grep every test that names these switches or asserts the legacy fallback for moment FD; pin or update with one line of justification each; scalar FD, fixed spin and closed-shell results must be bit-identical (they do not read these switches: prove with reference_identity, bench_scf --require-identical against both baselines).
OWNED_FILES:
  - crates/dft-numerics/src/switches.rs   # the three entries only
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/src/driver_select.rs   # l2/G4B only
  - crates/dft-cli/tests/   # pins/expectations only
  - crates/dft-engine/tests/   # pins only
  - docs/decisions/0045-broyden-history-conditioning.md
  - docs/design/reference-scf.md   # "Spin and Fermi-Dirac" status lines only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (reference baseline; and BAANDI_REFERENCE_STRICT=on with baseline-oracle.json)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: moment FD routes to reference by default; ADR 0045 accepted; other results bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
