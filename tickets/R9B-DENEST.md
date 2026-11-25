TASK_ID: R9B-DENEST
GOAL: |
  Stage R9b of docs/roadmap/refactor-r-plan.md: finish the R9a extraction and make refscf a top-level module independent of the legacy scf module, so R9c can put the legacy/oracle files behind a Cargo feature. Pure moves: bit-identical by construction.
READ FIRST: docs/roadmap/refactor-r-plan.md (section 0 finding 1, section 3 layout, section 4 "R9"); tickets/R9A-remaining-scope.patch (prepared by R9a: StartDensity/StartWavefunctions/AtomicStartReport moves out of atomic_start.rs/start_density.rs and path-only substitutions in refscf/driver.rs and electronic.rs; it may no longer apply cleanly after R8c/R10a/R6a, so re-derive it); crates/dft-engine/src/{lib.rs,scf.rs (the #[path = "refscf/mod.rs"] nesting),model/,atomic_start.rs,start_density.rs}; every `crate::scf::`, `super::super::`, `crate::session::`, `crate::collinear::`, `crate::finite_temperature::`, `crate::cutoff_ladder::`, `crate::gamma_real::` use inside refscf/ and refstart.rs.
TASK: |
  1. Apply the R9a remaining scope (re-derived): StartDensity, StartWavefunctions, AtomicStartReport into crate::model with re-exports at the old paths.
  2. De-nest: refscf (and atomic_start if it is shared) become top-level modules declared in lib.rs; remove the #[path] nesting in scf.rs; keep `pub use` at dft_engine::scf::refscf (or wherever external crates import it) so dft-cli, baandi-py and baandi-mcp do not change.
  3. refscf and refstart depend only on crate::model, crate::refscf::*, shared physics modules (energy, forces, stress, functional, grid_dist, parallel, profile, phase_clock), never on legacy driver modules (scf, session, convergence, collinear, finite_temperature, cutoff_ladder, gamma_real, m2, mix_dist, band_density). Add a test (grep-style, under crates/dft-engine/tests) that fails if refscf/*.rs or refstart.rs name one of those legacy modules; list any unavoidable exception with a reason.
  4. Identity: identical test counts before/after (report them), reference_identity, bench_scf --require-identical for both baselines.
OWNED_FILES:
  - crates/dft-engine/src/   # moves, re-exports and use lines only; no logic change
  - crates/dft-engine/tests/refscf_independence.rs   # new
TEST_COMMANDS:
  - cargo check --workspace --tests
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (reference baseline; and BAANDI_REFERENCE_STRICT=on with baseline-oracle.json)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: refscf top-level and independent of legacy driver modules (test enforced); external crates unchanged; identical counts and baselines.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
