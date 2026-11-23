TASK_ID: I25-VERIFY-COST
GOAL: |
  Reduce the wall cost of the interval from the first stopping-estimator map to the verified result in the reference SCF driver, WITHOUT changing the passing criterion (fixed-point check of the delivered density: Coulomb error of F(rho_out) - rho_out below the standard 5e-9 Ha at the next-ratchet tolerance, ADR 0043 / docs/design/reference-scf.md). Defaults unchanged; each idea is its own switch, default off, so the orchestrator can A/B it on dev08.
  Measured on dev08 MPI8 (main 5920aee, issue 25 comment): the interval is 25-28% of delivered wall. Si64: stop map 11, pass at 13, 3 attempts, errors 5.80e-8, 6.75e-9, 2.09e-9 Ha; verification 19% of wall, SCF after a failed attempt 9.2%. Si128: stop 15, pass 17, errors 1.51e-7, 3.09e-8, 1.47e-9; verification 17%, recovery 7.7%. Displaced Si8 (Mac fixture): stop 10, pass 12, 3 attempts.
  Observation: every failed attempt already computes one full map F(rho_out) at rho_out, and then an ordinary SCF map follows; the first attempt fails by 10-30x.
TASK: |
  Implement two switchable variants (new entries in crates/dft-numerics/src/switches.rs, default off; update the switch-count tests):
  V1 REFSCF_VERIFY_REUSE: after a failed attempt, the verification map's output is used as SCF work instead of being discarded: take rho_out as the next input density with its already computed output (append the pair to the Broyden history as an ordinary iterate, so the next ordinary map is not needed for that step). Specify in the design doc exactly how the history, ethr ratchet and energy correction treat this iterate. The next attempt still runs on the newly delivered density with the unchanged criterion.
  V2 REFSCF_VERIFY_GATE: skip an attempt (continue ordinary SCF) when a prediction says it will fail. The predictor uses only quantities the driver already has (e.g. the dr2 sequence / observed contraction, the previous attempt's error); it must be a general rule with no system-specific constant, record predicted versus actual for every attempt that runs, and fall back to attempting when the predictor is undefined. A skipped attempt is recorded as skipped. Never deliver without a passing attempt.
  Both variants: delivered state is always the state that passed; record attempts, skipped attempts and the stop-to-verified timestamps as now (refscf-check --stop-to-verified must keep working; add a "skipped" count if V2 is on).
  Mac checks (user rule 2026-11-24: no full gate or release matrices on the Mac): cargo check, the refscf and stop_to_verified tests, and new unit tests on the existing small fixtures (si8 displaced in debug if it fits ~10 s, otherwise #[ignore] with its release command) showing for each variant: same pass criterion, energies within the standard class of the default run, number of maps and attempts. Write docs/performance/lanes/RSCF/verify-cost-dev08-request.md: A/A then crossed A/B per variant against default on Si64 and Si128 MPI8 plus a skewed and a slab cell (generality rule), one case per <=600 s job, verification on, lightweight diagnostics, all output in external wall; report interval parts, maps, attempts, delivered energy error versus strict reference.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-cli/src/refscf_check.rs   # skipped-attempt field only
  - crates/dft-engine/tests/refscf.rs
  - crates/dft-cli/tests/stop_to_verified.rs
  - docs/design/reference-scf.md   # a "Verification cost variants" section only
  - docs/performance/lanes/RSCF/verify-cost-dev08-request.md
TEST_COMMANDS:   # Mac rule (user, 2026-11-24): only these; never scripts/check.sh full or release matrices on the Mac
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli stop_to_verified
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
  - (orchestrator, on dev08) scripts/check.sh full, release tests, A/B
ACCEPTANCE_CRITERIA: both switches default off and bit-identical defaults; criterion unchanged; skipped/reused attempts recorded; tests show map/attempt counts per variant; dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
