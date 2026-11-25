TASK_ID: R4A-REF-NSCF
GOAL: |
  Stage R4a of docs/roadmap/refactor-r-plan.md: non-self-consistent band solve on the reference path (engine only), to replace legacy nscf.rs for bands and DOS.
READ FIRST: docs/roadmap/refactor-r-plan.md (sections 2 row L3, 4 "R4"); crates/dft-engine/src/nscf.rs (legacy semantics: fixed density, charge check, at least four empty bands, residual-certified eigenvalues); crates/dft-engine/src/refstart.rs (k start, atomic-orbital start vectors); crates/dft-engine/src/refscf/k_driver.rs, k_map.rs; crates/dft-solvers/src/ref_davidson_complex.rs. Never read QE source.
TASK: |
  1. crates/dft-engine/src/refscf/nscf.rs: nscf(input, density, points, bands, eigen_tolerance) -> per-k eigenvalues (and optionally eigenvectors), with the potential built once from the supplied density (same field code as the SCF map, including NLCC, PBE, spin channels with per-channel density), no mixing, every requested band converged to the tolerance with ref_davidson_complex (Gamma uses complex too unless the real path is trivially reusable), k points processed sequentially on all ranks as in k_driver; same charge check and empty-band rule as legacy.
  2. A new refstart entry that builds the k start from a supplied density (sharing the atomic-orbital start vectors) without running the SCF start.
  3. Tests: on the same input density, reference NSCF eigenvalues equal legacy nscf within 1e-8 Ha (both residual-certified at 1e-8 Ha) on Si2 with a band path (k list incl. Gamma, X, L), a skew two-species synthetic cell, and a PBE NLCC synthetic case; spin two-channel case; ranks 1/3/8.
OWNED_FILES:
  - crates/dft-engine/src/refscf/nscf.rs   # new
  - crates/dft-engine/src/refscf/mod.rs   # module declaration only
  - crates/dft-engine/src/refstart.rs   # new entry only
  - crates/dft-engine/tests/refscf_nscf.rs   # new
  - docs/design/reference-scf.md   # a "Non-self-consistent bands" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf_nscf
  - cargo test -p dft-engine --test refscf
  - cargo fmt --check && cargo clippy -p dft-engine --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: reference NSCF matches legacy within 1e-8 Ha on the fixtures; no change to SCF paths.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
