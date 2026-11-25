TASK_ID: G4F-STIFFNESS-ROBUST
GOAL: |
  dev08 gen20: O atom (PseudoDojo O.upf PBE+NLCC, 12 bohr box, ecut 84 Ry, 72^3) with Fermi-Dirac theta 0.005 Ha and fixed moment 2.0 on the default reference path (G4B switches now on by default) fails on all 8 ranks after ~6.5 s with `Numerical("reference SCF: spin stiffness evaluator failed on a rank")`. The fixed-sector O case converges (gen19) with ADR 0046 screening (vacuum points below 1e-10 screened). The G4B spin-stiffness estimator term (refscf/electronic_map.rs spin_stiffness, finite differences of the polarized LDA-PW potentials) evidently does not apply the ADR 0046 screening rule and fails on vacuum/negative densities; and its error message is generic.
READ FIRST: docs/decisions/0045 and 0046; docs/design/reference-scf.md "Spin and Fermi-Dirac" (stiffness term) and the vacuum densities paragraph; crates/dft-engine/src/refscf/electronic_map.rs (spin_stiffness, SpinFields) and how the main spin XC evaluation applies screening; crates/dft-physics/src/spin.rs and xc.rs (screening hook).
TASK: |
  1. Apply the ADR 0046 rule inside the stiffness evaluator consistently with the main XC evaluation: points/channels screened in the XC evaluation contribute zero stiffness; the finite-difference step must not push a valid point below zero (use a one-sided or scaled step near the threshold, documented); never fail on vacuum points the main XC accepted.
  2. Error messages: name the stage (start/map/verification), the channel, the minimum density and grid point, and the rank-local reason, as the main XC errors do since G4E.
  3. Tests: a synthetic isolated open-shell atom with core charge in a box with Fermi-Dirac and a fixed moment (and a free-moment LDA-PW variant) that reproduces the failure before the fix and converges after it, reference vs legacy within the standard class where legacy converges; G4B 48/48 synthetic matrix unchanged (bit-identical fingerprints for cases without screened points).
OWNED_FILES:
  - crates/dft-engine/src/refscf/electronic_map.rs   # stiffness evaluator only
  - crates/dft-engine/src/refscf/electronic.rs   # error reporting only
  - crates/dft-physics/src/spin.rs   # screening hook for the response only, if needed
  - crates/dft-engine/tests/refscf.rs   # new tests appended
  - docs/decisions/0045-broyden-history-conditioning.md   # a note on the stiffness term and screening
  - docs/design/reference-scf.md   # "Spin and Fermi-Dirac" stiffness paragraph only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test --release -p dft-engine --test refscf g4 -- --include-ignored --test-threads=1
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-physics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: isolated open-shell FD-moment fixture converges; precise error messages; G4B matrix unchanged.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
