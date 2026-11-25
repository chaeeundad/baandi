TASK_ID: G2B-PBE-NLCC
GOAL: |
  Generalization stage G2B (docs/roadmap/reference-generalization.md): PBE together with NLCC in the reference driver. Today refstart.rs rejects it explicitly ("PBE with core charge (NLCC) is not supported yet") because the GGA evaluator needs the core density's gradient. PseudoDojo files (PBE, mostly with NLCC) need this for any non-Si material.
READ FIRST: docs/design/reference-scf.md "NLCC" and "PBE"; crates/dft-physics/src/gga/ (distributed evaluator); how the legacy path builds PBE with core charge (gradient of rho + rho_core, XC energy, double counting excludes the core). Never read QE source.
TASK: |
  1. Add rho_core(G) to the gradient input (rho_total(G) = rho(G) + rho_core(G)) and rho_core(r) to the XC density on the owned slabs, for the start and every map (cache rho_core(G) once per run); XC potential acts on the valence; double counting and the energy audit exactly as the legacy PBE+NLCC path.
  2. Remove the explicit rejection; keep spin-polarized PBE rejected (G4).
  3. Tests (Mac): synthetic PBE NC fixture with core charge (common::potential_with_nlcc with "PBE"): reference vs legacy delivered energy within the standard class and within 1e-10 Ha tight, ranks 1/3/8; gradient of the total density vs a serial oracle; existing LDA/NLCC and PBE-without-core results unchanged.
  4. Update docs/design/reference-scf.md ("PBE" section, NLCC paragraph).
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs   # NLCC+PBE parts only
  - crates/dft-engine/src/refscf/   # fields (GGA with core) only
  - crates/dft-physics/src/gga/   # only if the evaluator needs a core-density input
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-engine/tests/refstart.rs   # update the rejection test only
  - docs/design/reference-scf.md   # "PBE"/"NLCC" sections only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-engine --test refstart
  - cargo test -p dft-physics
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-physics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: PBE+NLCC runs in the reference driver with legacy agreement on fixtures; other paths unchanged.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
