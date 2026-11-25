TASK_ID: G5-FORCES-STRESS
GOAL: |
  Generalization stage G5 (docs/roadmap/reference-generalization.md): forces, stress and relaxation from a delivered reference-SCF state (closed shell, Gamma and k meshes; spin/FD come with G4, design the interface so they plug in). Today forces/stress take the legacy state and accept only certified states; relax.rs drives the legacy SCF.
READ FIRST: docs/design/reference-scf.md (verification, delivered state, K points); crates/dft-engine/src/forces.rs, stress.rs, force_identity.rs; crates/dft-cli/src/relax.rs and run.rs (how forces/stress are requested and written); docs/roadmap/post-performance-tasks.md (band-first nonlocal force contraction note); ADRs on force/stress certificates (grep docs/decisions for force). Never read QE source.
TASK: |
  1. Build ForceState/StressState from the reference driver's delivered state (orbitals per k, occupations, density, potentials, cell, basis) on the distributed layout, without gathering full orbitals on one rank if avoidable; if the legacy kernels need replicated data, gather once at delivery and document memory.
  2. Accept the reference verification (fixed-point verified, ADR 0043/0044) as the certificate for forces/stress; keep rejection for unverified states unless the user sets verify=off explicitly (then record "unverified").
  3. Nonlocal force contraction band-first: one GEMM over bands for <beta|psi> derivatives instead of per-band loops, behind a switch (default off until measured; no timing now), identity test against the current contraction.
  4. relax.rs: run each ionic step with the reference driver; warm start from the previous step's orbitals and density (reuse the existing start seed/prolongation mechanism); per-step records unchanged in format.
  5. Tests (Mac): forces and stress on Si2 displaced, skew cell, a synthetic nonlocal fixture: reference-state vs legacy-state forces within 1e-6 Ha/Bohr and stress within 1e-6 Ha/Bohr^3 (standard class), finite-difference checks of the reference energy (existing helpers), ranks 1/3/8; a 3-step relaxation on Si2 matching the legacy trajectory within the class.
OWNED_FILES:
  - crates/dft-engine/src/forces.rs
  - crates/dft-engine/src/stress.rs
  - crates/dft-engine/src/force_identity.rs
  - crates/dft-engine/src/refscf/delivered.rs   # new: delivered-state adapter (create; do not edit other refscf files except mod.rs declaration)
  - crates/dft-engine/src/refscf/mod.rs   # module declaration only
  - crates/dft-cli/src/relax.rs
  - crates/dft-cli/src/run.rs   # forces/stress request wiring for the reference driver only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-engine/tests/refscf_forces.rs   # new test file
  - docs/design/reference-scf.md   # a "Forces, stress and relaxation" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf_forces
  - cargo test -p dft-engine forces
  - cargo test -p dft-cli relax
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: forces, stress and relaxation from the reference driver agree with legacy; certificate rule documented; band-first contraction behind a switch with identity test.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
