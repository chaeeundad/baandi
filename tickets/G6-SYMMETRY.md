TASK_ID: G6-SYMMETRY
GOAL: |
  Generalization stage G6 (docs/roadmap/reference-generalization.md): crystal symmetry in the reference SCF driver: irreducible k points with symmetry weights (from the existing KOrbitMap), density symmetrization over the crystal group each map, and symmetrized forces/stress (G5). Today the reference driver rejects symmetry; the legacy path supports it (PreparedSymmetry in crates/dft-engine/src/scf.rs, crates/dft-physics/src/symmetry.rs, crates/dft-geometry/src/symmetry.rs, ADRs 0011/0013).
READ FIRST: docs/design/reference-scf.md (whole, incl. K points, Spin and Fermi-Dirac, Forces sections); crates/dft-engine/src/scf.rs PreparedSymmetry and how the legacy driver symmetrizes the density (in G space or real space, with fractional translations) and where; crates/dft-physics/src/symmetry.rs; ADR 0011, 0013. Never read QE source.
TASK: |
  1. k reduction: the reference k driver takes the irreducible list and weights from the prepared symmetry (same source as legacy); time reversal unchanged.
  2. Density symmetrization on the distributed G layout: apply the group operations to rho(G) of the owned sticks; operations map G to G' owned by other ranks, so build a deterministic exchange plan once per run (one all-to-all per map at most), and symmetrize rho_out before mixing exactly where the legacy driver does (and the start density if legacy does). For spin: both channels (collinear magnetization is a scalar under the group).
  3. Forces and stress from G5: symmetrize as legacy does.
  4. Records: group order, number of irreducible k, symmetrization residual norm per map (cheap diagnostic).
  5. Tests (Mac): Si2 fcc 4x4x4 with symmetry vs without symmetry (full mesh): delivered energy within the standard class and 1e-10 Ha tight, density difference; skew low-symmetry cell (group order small) and a cell with fractional translations (diamond already has them); reference vs legacy with symmetry within the class, ranks 1/3/8; spin and FD case with symmetry; Gamma-only bit-identical to main.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/src/refstart.rs   # symmetric start density only
  - crates/dft-physics/src/symmetry.rs   # distributed exchange plan only, serial behaviour unchanged
  - crates/dft-physics/src/distributed_complex.rs   # symmetrization hook only
  - crates/dft-engine/tests/refscf.rs   # new tests appended at the end only
  - crates/dft-cli/src/run.rs   # reference-driver symmetry wiring and records only
  - docs/design/reference-scf.md   # a "Symmetry" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-physics symmetry
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli reference
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-physics -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: symmetry-reduced reference runs agree with full-mesh and legacy runs; Gamma bit-identical; design documented.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
