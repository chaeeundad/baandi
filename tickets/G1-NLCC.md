TASK_ID: G1-NLCC
GOAL: |
  First stage of generalizing the accelerated reference SCF path (docs/roadmap/reference-generalization.md, user decision 2026-11-26): support norm-conserving pseudopotentials WITH nonlinear core correction (NLCC) in the reference driver. Today crates/dft-engine/src/refstart.rs rejects them ("requires compatible NC atomic data without core charge").
READ FIRST: docs/roadmap/reference-generalization.md; docs/design/reference-scf.md (Start, map, energy terms, Per-map energy cost); how the legacy path handles NLCC (search the engine for nlcc / core charge: XC potential and energy with rho + rho_core, start density, energy audit, forces excluded here); README "Nonlinear core correction and PseudoDojo NC file support" and its evidence; docs/performance/AGENT-BRIEF.md (hot-path rules). Never read QE source.
TASK: |
  1. Reference start and map: build rho_core on the owned real-space slabs (same structure-factor machinery, reuse the start's tables), add it to the density in the XC potential and XC energy only (never to the Hartree term, the charge, the mixer or dr2), and in the energy audit / double counting exactly as the legacy path defines it. LDA only in this stage (PBE is G2).
  2. Remove the rejection for core charge; keep rejecting what is still unsupported, with clear messages.
  3. Tests (Mac): synthetic NLCC fixture(s) with analytic or legacy-path oracles: reference vs legacy delivered energy within the standard class (and within 1e-10 Ha with tight settings), density and eigenvalues agreement, simulated ranks 1/3/8, rank-count rule; an LDA NC file with core charge if available in the repo test data (never commit UPF files; synthetic only otherwise). Existing Si (no core charge) results bit-identical.
  4. docs/design/reference-scf.md section "NLCC"; a dev08 request (docs/performance/lanes/RSCF/g1-nlcc-dev08-request.md): one run of a PseudoDojo-based cell with NLCC under LDA is NOT physical (PBE files); instead request the legacy-vs-reference comparison on the synthetic fixtures at MPI8 and, if an LDA NLCC UPF exists on dev08, one real cell.
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-engine/tests/common/   # new fixture helpers only
  - docs/design/reference-scf.md   # "NLCC" section only
  - docs/performance/lanes/RSCF/g1-nlcc-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-engine
  - cargo fmt --check && cargo clippy -p dft-engine --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: NLCC supported in the reference driver with legacy agreement on fixtures; no change for non-NLCC inputs; docs and dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
