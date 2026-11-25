TASK_ID: R6B-COMPAT-MODE
GOAL: |
  Stage R6b of docs/roadmap/refactor-r-plan.md (row L10): `compat.mode = "reference"` (ADR 0038: original PW92 prefactor, QE local mesh and radial rule) runs on the reference driver; refstart currently rejects pw92 != Modified, so these inputs fall back to legacy (L10).
READ FIRST: docs/decisions/0038-reference-compat-mode.md; crates/dft-engine/src/refstart.rs (pw92 rejection), refscf/map.rs and electronic_map.rs (XC evaluation), crates/dft-physics (PW92 prefactor variants, local mesh/radial rule options); crates/dft-cli/src/driver_select.rs (l10) and config.rs (compat mode); tests crates/dft-cli/tests/reference_mode_restart.rs and all/reference_mode_*.
TASK: |
  1. Thread the PW92 prefactor choice, the local mesh and the radial integration rule from ScfOptions through refstart and the map/electronic XC evaluators (they are physics constants/choices, not new algorithms). Remove the rejection and the L10 predicate.
  2. Tests: compat-mode Si2 (Gamma and k222) reference vs strict oracle within 5e-9 Ha standard and 1e-10 Ha tight; default (Modified) path bit-identical (reference_identity, bench_scf --require-identical both baselines); the reference_mode tests pass on the reference driver (pin to legacy only those that test legacy-only restart v1).
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs   # compat options only
  - crates/dft-engine/src/refscf/   # XC option plumbing only
  - crates/dft-physics/src/   # only if an option is not reachable from the distributed evaluators
  - crates/dft-cli/src/driver_select.rs   # L10 only
  - crates/dft-cli/tests/   # compat tests and routing expectations only
  - crates/dft-engine/tests/refscf.rs   # new tests appended
TEST_COMMANDS:
  - cargo check --workspace --tests
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (both baselines)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: compat mode on the reference driver agrees with the oracle; default bit-identical; L10 removed.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
