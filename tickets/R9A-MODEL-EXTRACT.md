TASK_ID: R9A-MODEL-EXTRACT
GOAL: |
  Stage R9a of docs/roadmap/refactor-r-plan.md (user decision D2: move the strict oracle whole): extract the data types and contracts shared by the reference path and the legacy/oracle code into crate::model, so that later the legacy files can move behind a feature (R9c). Pure moves: bit-identical by construction.
READ FIRST: docs/roadmap/refactor-r-plan.md (section 0 finding 1, section 3 layout, section 4 "R9"); crates/dft-engine/src/{lib.rs,scf.rs,session.rs,convergence.rs,collinear.rs,finite_temperature.rs,cutoff_ladder.rs,gamma_real.rs}; how refscf uses them (grep "crate::scf::", "cutoff_ladder::", "gamma_real::", "collinear::", "finite_temperature::", "session::" under crates/dft-engine/src/refscf and refstart.rs).
TASK: |
  1. Create crates/dft-engine/src/model/{mod.rs,options.rs,result.rs,verify.rs,spin.rs,error.rs} and move into them: ScfOptions and its option enums (Mixing, CompatMode, GammaReal, StartDensity/Wavefunctions, OccupationPolicy...), ScfResult/ScfOutcome/KPointSolution/IterationRecord/PhaseTimings/FixedPointCertificate/AtomicStartReport, Verification* and StandardProfile and assess_fixed_point, SpinSector/ThermodynamicLedger/finite-T validate_options, ScfError. Keep `pub use` re-exports at every old path so no external caller changes.
  2. Move the prolongation helpers used by the reference ladder (cutoff_ladder::prolong_*) into refscf/prolong.rs, and LazyOrbitals/set_delivered_orbitals/lazy_result_allowed (gamma_real) into refscf/lazy.rs, with re-exports at the old paths.
  3. refscf and refstart import from crate::model / refscf::{prolong,lazy} instead of crate::scf::* (keep the #[path] nesting of refscf in scf.rs for now; R9b de-nests).
  4. No logic change. Identity: reference_identity, bench_scf --require-identical against both baselines, full check with identical passed/ignored counts (report before/after counts).
OWNED_FILES:
  - crates/dft-engine/src/model/   # new
  - crates/dft-engine/src/lib.rs
  - crates/dft-engine/src/scf.rs
  - crates/dft-engine/src/session.rs
  - crates/dft-engine/src/convergence.rs
  - crates/dft-engine/src/collinear.rs
  - crates/dft-engine/src/finite_temperature.rs
  - crates/dft-engine/src/cutoff_ladder.rs
  - crates/dft-engine/src/gamma_real.rs
  - crates/dft-engine/src/refstart.rs   # imports only
  - crates/dft-engine/src/refscf/   # imports, new prolong.rs and lazy.rs only (R8c is editing stepper/driver/k_driver/electronic concurrently: touch only their use lines)
TEST_COMMANDS:
  - cargo check --workspace --tests
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (reference baseline; and BAANDI_REFERENCE_STRICT=on with baseline-oracle.json)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: shared surface in crate::model with old paths re-exported; refscf no longer imports crate::scf::*; identical test counts and baselines.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
