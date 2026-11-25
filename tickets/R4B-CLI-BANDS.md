TASK_ID: R4B-CLI-BANDS
GOAL: |
  Stage R4b of docs/roadmap/refactor-r-plan.md: the CLI bands and dos commands use the reference path: SCF through refscf::api::solve (R3) and the band solve through refscf::nscf (R4a, merged). Remove the L3 routing predicate.
READ FIRST: docs/roadmap/refactor-r-plan.md (section 2 row L3, 4 "R4"); crates/dft-cli/src/bands.rs, dos.rs (bands::Context); crates/dft-engine/src/refscf/{api.rs,nscf.rs}; crates/dft-cli/src/driver_select.rs; output format writers used by bands/dos (all formats must stay unchanged: CSV, dat, agr, bands.dat.gnu, EIGENVAL/DOSCAR, JSON, SVG).
TASK: |
  1. bands/dos: when the selected driver is reference, run the SCF through api::solve and the path/mesh through refscf::nscf with the delivered density (spin channels too); legacy and strict-oracle selections keep the current code. Record scf.driver for both the SCF and the NSCF stage.
  2. Remove the L3 predicate from driver_select.rs and update the routing corpus expectations.
  3. Tests: existing bands_dos format tests pass; reference vs legacy-pinned output of the same input agrees within 1e-7 Ha per eigenvalue (standard-class SCF density difference propagated) on Si2 band path and Si2 DOS; spin case.
OWNED_FILES:
  - crates/dft-cli/src/bands.rs
  - crates/dft-cli/src/dos.rs
  - crates/dft-cli/src/driver_select.rs   # L3 only
  - crates/dft-cli/tests/fixtures/routing/   # L3 expectations only
  - crates/dft-cli/tests/driver_routing.rs   # L3 expectations only
  - crates/dft-cli/tests/all/bands_dos.rs
  - crates/dft-cli/tests/   # new test file only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-cli
  - cargo fmt --check && cargo clippy -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: bands/dos default to the reference path with unchanged formats and agreement with legacy; L3 routing removed.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
