TASK_ID: R7A-FROZEN-STRAIN
GOAL: |
  Stage R7a of docs/roadmap/refactor-r-plan.md: the reference start and delivered state honour ScfOptions::strain_reference (frozen plane-wave basis from the reference cell evaluated at the strained cell), needed for variable-cell relaxation on the reference path. refstart.rs currently rejects strain_reference.
READ FIRST: docs/roadmap/refactor-r-plan.md (section 2 row L4, section 4 "R7"); how legacy implements the frozen Miller basis under strain (crates/dft-engine/src/scf.rs strain_reference, crates/dft-basis); crates/dft-cli/src/relax/cell.rs; crates/dft-engine/src/refstart.rs (rejection ~:398-422); refscf/delivered.rs (stress state); tests/reference/relaxation/cell/validation.json protocol. Never read QE source.
TASK: |
  1. Gamma real and complex k starts: select plane waves on the reference cell's Miller set and evaluate kinetic, structure factors, nonlocal projectors and the density grid at the strained cell, exactly as legacy does; remove the rejection.
  2. delivered::stress_state at the strained cell with the frozen basis (Pulay handling as legacy).
  3. Tests: strained Si2 (Gamma and k222) reference vs strict oracle: energy within the standard class and 1e-10 Ha tight, stress within 1e-8 Ha/Bohr^3 tight; finite-difference check of the frozen-basis energy derivative against the stress; unstrained runs bit-identical.
  If the ticket overruns, finish Gamma completely and report the k part as remaining.
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs   # strain support only
  - crates/dft-engine/src/refscf/delivered.rs
  - crates/dft-engine/src/refscf/k_driver.rs   # strain plumbing only
  - crates/dft-engine/src/refscf/driver.rs   # strain plumbing only
  - crates/dft-engine/tests/refscf_forces.rs   # new tests appended
  - docs/design/reference-scf.md   # "Forces, stress and relaxation" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf_forces
  - cargo test -p dft-engine --test refscf
  - cargo fmt --check && cargo clippy -p dft-engine --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: strained reference runs agree with the oracle; unstrained bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
