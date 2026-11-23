TASK_ID: I25-FIRST-ATTEMPT
GOAL: |
  In the reference SCF driver, the FIRST fixed-point verification attempt after the stopping estimator fires fails by 10-30x on every measured case (dev08 MPI8: Si64 5.80e-8 Ha, Si128 1.51e-7, elongated similar; pass threshold 5e-9 Ha), costing one full map each time. REFSCF_VERIFY_GATE (default on since 2026-11-24, docs/design/reference-scf.md "Verification cost variants") skips only LATER attempts, because its predictor needs a previous attempt. Goal: a predictor that also covers the first attempt, so a first attempt that will fail is replaced by ordinary SCF work, with the passing criterion unchanged and never delivering without a passing attempt.
TASK: |
  1. Derive a first-attempt prediction from quantities the driver already has, as a general rule with no system-specific constant. Candidates to evaluate (pick by evidence, document the reasoning in the design doc): (a) the modified Broyden history's approximate inverse Jacobian applied to the current residual, giving an estimate of the fixed-point residual of rho_out in the Hartree metric; (b) the observed contraction of the dr2 sequence over the last maps; (c) any combination. Record the prediction for every attempt that runs (predicted versus actual) so the rule can be audited.
  2. Implement it as a new switch REFSCF_VERIFY_GATE_FIRST, default off, active only together with REFSCF_VERIFY_GATE; keep the "no two consecutive skips" safety and the fallback to attempting when undefined.
  3. Fix the label in the si8 release test (crates/dft-engine/tests/refscf.rs stop_to_verified_si8_displaced): the arm with both switches explicitly off must be labelled "off", and a true "default" arm (no explicit switches) must exist.
  4. Mac checks only (user rule 2026-11-24): cargo check, cargo test -p dft-engine --test refscf, cargo test -p dft-cli stop_to_verified, cargo test -p dft-numerics switches, fmt, clippy -D warnings on the touched crates. Show on the existing synthetic fixtures the predicted versus actual error of the first attempt and the map/attempt/skip counts per arm.
  5. Write docs/performance/lanes/RSCF/first-attempt-dev08-request.md: one untimed run per arm (default vs default+GATE_FIRST) on Si64, Si128, elongated, rskew MPI8 with predicted/actual per attempt; then A/A + crossed A/B as in verify-cost-dev08-request.md only if the untimed runs show fewer evaluated maps on Si64 and Si128.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-cli/src/run.rs   # reference_record_json / reference_scf summary fields only
  - crates/dft-cli/src/refscf_check.rs   # new prediction fields only
  - crates/dft-engine/tests/refscf.rs
  - crates/dft-cli/tests/stop_to_verified.rs
  - docs/design/reference-scf.md   # "Verification cost variants" section only
  - docs/performance/lanes/RSCF/first-attempt-dev08-request.md
TEST_COMMANDS:   # Mac rule (user, 2026-11-24): only these; never scripts/check.sh full or release matrices on the Mac
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli stop_to_verified
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
  - (orchestrator, on dev08) scripts/check.sh full, release tests, A/B
ACCEPTANCE_CRITERIA: switch default off, defaults bit-identical, criterion unchanged, predictions recorded, si8 labels fixed, dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
