TASK_ID: P3-VERWARM
GOAL: |
  Verification is now the largest gap in the verified headline (dev08 one run at ac6aee9c, si128 MPI8: verified 30.15 s vs verify-off 23.59 s, QE rmm 21.57 s; si64 6.40 vs 4.57). After the estimator stop the driver runs a fixed-point verification map at rho_out (criterion unchanged: Coulomb error < 5e-9 Ha, user decision C), and on failure ordinary SCF maps plus further attempts (si128: 17 maps, 2 attempts, 1 gated skip). Reduce the cost of each verification map and of the recovery maps without changing the passing criterion, the delivered state, or the independence of the check (the verification map recomputes the density from a fresh eigensolve at rho_out).
READ FIRST: docs/design/reference-scf.md (Verification, VERIFY_GATE, Kept HX warm start, stop-to-verified), docs/decisions/0043 and 0044, attempts-ledger rows I25-VERIFY-COST, P3r, KEPTHX2. Never read QE source.
TASK: |
  1. Account a verification map on si128-sized work (Mac fixtures plus the dev08 numbers above): which parts are fresh by design (density at rho_out) and which are fresh only by implementation (the starting H application of the eigensolve, the starting block, the ethr used). P3-KEPTHX deliberately disabled kept-HX for verification maps; argue whether a warm start (previous Ritz block and its H products updated by dV = V[rho_out] - V_old) compromises independence: the check is that F(rho_out) ~ rho_out, i.e. the eigensolve at the potential of rho_out converges to the same density; a warm start changes only the eigensolver's starting guess, not the operator or the convergence tolerance.
  2. Behind a switch REFSCF_VERIFY_WARM (default off): warm-start verification and post-failure recovery maps with kept HX; keep the eigensolve tolerance rule of the verification map unchanged; record per attempt whether it was warm and the H columns spent.
  3. Tests (Mac only): verification decisions and delivered energies identical within rounding to the cold path on the fixtures (cubic, skew, elongated, synthetic nonlocal; simulated ranks 1/3/8); fewer verification H columns; stop-to-verified reader still closes.
  4. docs/performance/lanes/RSCF/verwarm-dev08-request.md: ONE run per arm on Si64 (work-order input) and si128 MPI8, D vs W, verify on; report verified wall, verification seconds, attempts, verification H columns, E - E_strict.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-cli/src/run.rs   # verification record fields only
  - crates/dft-cli/src/refscf_check.rs   # warm flag only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # a "Warm verification maps" section only
  - docs/performance/lanes/RSCF/verwarm-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli stop_to_verified
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: independence argument documented; switch default off; identical verification decisions on fixtures with fewer verification H columns; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
