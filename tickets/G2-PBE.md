TASK_ID: G2-PBE
GOAL: |
  Stage G2 of docs/roadmap/reference-generalization.md: PBE in the reference SCF driver. refstart.rs rejects PBE ("reference start: PBE requires distributed gradients"). The legacy path has PBE through the xcx library (validated against QE on Si).
READ FIRST: docs/roadmap/reference-generalization.md; how the legacy path builds the PBE potential and energy (gradient of the density by FFT, the xcx call, divergence term, energy audit); docs/design/reference-scf.md (Start, map fields, Per-map energy cost, distributed fields/slabs); docs/performance/AGENT-BRIEF.md. Never read QE source. G1-NLCC runs in parallel and edits the XC/core-charge handling in refstart.rs and refscf/: keep your edits to the gradient/GGA parts and expect a merge.
TASK: |
  1. Distributed density gradient on the owned slabs: forward FFT of the density (already available as rho(G) in the map), multiply by iG per Cartesian component, inverse FFT to the slabs (three transforms), xcx PBE on the owned points, then the divergence term for the potential (three more transforms), all on the existing distributed field plans (no full-grid gathers). Same for the start density.
  2. Energy audit and double counting exactly as the legacy PBE path defines them; LEAN_ENERGY rules unchanged (full audit at delivery/verification).
  3. Remove the PBE rejection; keep spin-polarized PBE rejected (G4).
  4. Tests (Mac): Si (Si.pz-vbc is LDA; use a synthetic pseudopotential or the PseudoDojo Si PBE file only if NLCC (G1) is not needed for it; otherwise test on synthetic fixtures) reference vs legacy PBE delivered energy within the standard class and within 1e-10 Ha tight; gradients vs a serial reference on several grids and simulated ranks 1/3/8.
  5. docs/design/reference-scf.md "PBE" section; docs/performance/lanes/RSCF/g2-pbe-dev08-request.md: legacy vs reference on a PBE Si cell at MPI8 (and after G1, a PseudoDojo PBE cell with NLCC such as Si or MgO), one run each.
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs   # PBE gradient/potential parts only
  - crates/dft-engine/src/refscf/   # fields (GGA) only
  - crates/dft-physics/src/   # only if a distributed gradient helper is added there
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-engine/tests/common/   # new fixture helpers only
  - docs/design/reference-scf.md   # "PBE" section only
  - docs/performance/lanes/RSCF/g2-pbe-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-physics
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-physics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: PBE in the reference driver with legacy agreement on fixtures; distributed gradients without full gathers; LDA results unchanged; dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
