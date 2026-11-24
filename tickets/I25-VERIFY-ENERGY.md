TASK_ID: I25-VERIFY-ENERGY
GOAL: |
  User decision 2026-11-25 (option B): keep ONE independent verification map in the reference SCF driver, but judge it by the energy change instead of the density fixed-point (Coulomb) error. Evidence (dev08 MPI8, issue 25): BAANDI's dr2 trajectory matches QE's estimated accuracy map by map and the estimator stops at the same map as QE (si128 15, si64 11); the fixed-point check then fails there by 12-30x (si128 1.51e-7 Ha, si64 5.80e-8 Ha vs 5e-9 Ha) and costs two more maps plus attempts (~25% of si128 wall), while the stop-map total energy is already inside the standard class (1.0e-9 to 3.3e-9 Ha from the strict reference on si64/si128/si216 np8 and si32/si64 np1).
TASK: |
  1. New verification judgment for the reference driver (switch REFSCF_VERIFY_ENERGY, default ON for the reference driver; the legacy path is unchanged): after the estimator stop, run the existing verification map at rho_out (same tolerance rule as now). Pass when |E(verification map) - E(delivered state)| < 5e-9 Ha, where both energies are the same functional form (use the QE-style variational energy the driver already reports; document exactly which terms and which densities). Keep computing and recording the Coulomb fixed-point error, the residual bound and the density residual as diagnostics. If the energy test fails, continue SCF exactly as today (failed attempts and recovery charged) and attempt again at the next stop. Delivered state stays the stop map's state.
  2. With REFSCF_VERIFY_ENERGY on, REFSCF_VERIFY_GATE's predictor (built for the Coulomb error) must not skip the first attempt under the energy test; define the interaction (simplest: the gate is inactive when the energy judgment is on) and record it.
  3. ADR docs/decisions/0044-reference-verification-energy.md: decision, evidence numbers above, what is still guaranteed (an independent map is always run; the energy of the delivered state is checked against a fresh map; the fixed-point error is recorded, not enforced), what is no longer guaranteed (density fixed point to 5e-9 Ha), and how to restore the old rule (switch off).
  4. Tests (Mac, user rule: only cargo check and the touched crates' tests; no check.sh full, no release matrices): on the existing small fixtures the energy judgment passes at the first attempt where the old rule needed retries; delivered energies unchanged within 5e-9 Ha versus the old rule; switch off reproduces the old behaviour bit-identically; refscf-check --stop-to-verified still works.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/src/session.rs   # assess_fixed_point or a sibling energy assessment only
  - crates/dft-engine/src/verification_profile.rs   # if the attempt record needs an energy field
  - crates/dft-cli/src/run.rs   # verification JSON fields only (another ticket edits the timing-table section of this file in parallel)
  - crates/dft-cli/src/refscf_check.rs
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-engine/tests/refscf.rs
  - crates/dft-cli/tests/stop_to_verified.rs
  - docs/decisions/0044-reference-verification-energy.md
  - docs/design/reference-scf.md   # a "Verification by energy" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli stop_to_verified
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: energy judgment default on for the reference driver, one verification map always run, diagnostics kept, switch off bit-identical to the old rule, ADR written.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
