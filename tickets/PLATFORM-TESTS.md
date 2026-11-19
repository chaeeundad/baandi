TASK_ID: PLATFORM-TESTS
GOAL: |
  The full gate now also runs on dev08 (x86_64 Linux, MKL 2026.1, OpenMPI) through the baandi-dev08 lane (commands in docs/performance/requests/gate-on-dev08.md). On commit e4186f9 the build, clippy and fmt pass there, but some tests fail because their expected values were recorded on the Mac (aarch64, pure-Rust backend) and assert exact bits or exact iteration counts. Make every test valid on both platforms WITHOUT weakening what it protects:
  - Exact-bit baselines recorded on aarch64 stay exact on aarch64 (use #[cfg(target_arch = "aarch64")] for the exact-bit assertion; on other targets assert the same quantity within the documented rounding class: energy 1e-12 Ha relative to the baseline for a toy cell, or the accuracy class of the test (strict 1e-10 Ha, standard 5e-9 Ha) with the reason in a comment). Never loosen on aarch64.
  - Iteration-count equalities that are really rounding-sensitive (trajectory) are asserted as `abs_diff <= 3` or by an energy-class check ONLY where the test comment already says the comparison is rounding-sensitive; where the test protects a work-count identity (same work with the switch on/off), keep it exact and instead compute the baseline arm and the tested arm on the SAME platform in the test (both arms run in the test) rather than comparing to a stored constant.
  Failures reported by dev08 (default features): (1) crates/dft-engine/tests/spin_share.rs off_preserves_main_bits: energy bits 13829793526668695542 vs BASELINE_ENERGY_BITS 13829793526668695510 (6 ulp) and environment_switch_is_used_without_override (child process of the same check); (2) crates/dft-cli/tests/all/run.rs reuse_psi_switch_keeps_energy_and_reuses_starts: iterations 29 vs 26 (passes in the mkl build). Failures with --features mpi,mkl: (3) crates/dft-numerics fused_pair_plan::tests::mkl_lanes_match_rust_lanes: `MklPairPlan::new([6,6,6]) is not Grid` at fused_pair_plan.rs:669 — investigate whether this is a TEST assumption (MKL lanes unsupported for that shape on x86) or a real bug in MklPairPlan (shape handling, pruned/lane variants); fix the cause if it is a bug; (4) crates/dft-solvers tests algebra_stream rounding_class_switches_keep_roots_and_thread_identity: clustered-96 iterations 39 vs 35; (5) crates/dft-solvers tests all distributed_algebra::gram_matches_serial_within_rounding_and_is_identical_on_every_rank: 1-rank vs serial GEMM differ by 1-2 ulp (the name says 'within rounding': find which assertion demands exact equality and whether it should be exact only on the Rust backend: with MKL the serial GEMM and the distributed 1-rank partial GEMM can order sums differently; keep 'identical on every rank' exact, but serial-vs-distributed within rounding); (6) crates/dft-solvers distributed_selected_roots selected_roots_match_the_full_rotation_path: residual 3.825675439e-7 vs ...441e-7 (an exact-ish comparison of a residual that must be a tolerance 1e-12 relative). MKL_CBWR=COMPATIBLE changes none of these.
OWNED_FILES:
  - crates/dft-engine/tests/spin_share.rs
  - crates/dft-cli/tests/all/run.rs
  - crates/dft-numerics/src/fused_pair_plan.rs   # tests module and, if the cause is a bug, the plan constructor
  - crates/dft-solvers/tests/algebra_stream.rs
  - crates/dft-solvers/tests/all/distributed_algebra.rs
  - crates/dft-solvers/tests/distributed_selected_roots.rs
  - docs/performance/slow-tests.md   # or docs/performance/test-targets.md: add a short 'platform-dependent baselines' note
  - docs/performance/attempts-ledger.md  # one row PLATFORM-TESTS
READ_ONLY_CONTEXT: AGENTS.md, crates/*/AGENTS.md of the touched crates, docs/performance/requests/gate-on-dev08.md
IMPLEMENTATION_STEPS:
  1. Read each failing test and the quoted values; decide per test: platform-gated exact bits, same-platform two-arm comparison, or a documented tolerance. Do not delete any test or lower an aarch64 expectation.
  2. For (3) read MklPairPlan::new and the dispatch that returns the Grid variant; explain in the report why [6,6,6] is not Grid with the mkl feature on x86 and fix the test or the code.
  3. Run on the Mac: scripts/check.sh fast (or `scripts/check.sh fast dft-engine dft-cli dft-solvers dft-numerics`): must stay green on aarch64 with bit-identical expectations there. The x86/MKL behaviour cannot be run in the sandbox: state your reasoning for each x86 expectation so the dev08 lane can confirm.
TEST_COMMANDS:
  - scripts/check.sh fast dft-engine dft-cli dft-solvers dft-numerics
ACCEPTANCE_CRITERIA: all listed tests pass on aarch64 unchanged in strength; on x86 the assertions are rounding-class; no test removed or ignored; clippy/fmt clean; a short list in the report of every assertion changed with its platform condition and class.
NON_GOALS: engine changes, new features.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
