TASK_ID: G4B-SPIN-MOMENT
GOAL: |
  Qualify collinear Fermi-Dirac with free and fixed moment in the reference driver (G4 left them routed to legacy behind REFSCF_SPIN_FD_MOMENT, default off). Measured failures (docs/design/reference-scf.md "Spin and Fermi-Dirac", G4 continuation review): free moment standard |dF| 1.12e-8 Ha > 5e-9 with the estimator satisfied, moment error 3.4e-4 mu_B; tight runs reject on singular Broyden histories (condition ~1e15-1e16 at iterations 10-13); fixed moment mesh weighted-occupation error 1.41e-4 > 1e-4 electrons. Acceptance bounds stay unchanged (5e-9 / 1e-10 Ha, 1e-4 mu_B, 1e-4 Ha, 1e-4 electrons).
READ FIRST: docs/design/reference-scf.md (Spin and Fermi-Dirac, Broyden), docs/design/onepath3-spin-metals.md (estimator section 2), crates/dft-engine/src/refscf/{electronic.rs,electronic_map.rs,broyden.rs,metric.rs}, the g4b_moment_acceptance test. Never read QE source.
TASK: |
  1. Diagnose before changing anything: for the failing cases record per map the charge, spin and uniform estimator terms, the true error against the strict oracle (energy, moment, mu), the Broyden history condition number, and the occupation residual. Answer: (a) does the run stop while the true error is still above the bound (estimator blind to the spin soft mode) or does it fail to converge; (b) where do singular histories come from (near-duplicate residual differences after convergence stalls, or the fixed-moment/zero-mode projection).
  2. Candidate fixes, each behind its own switch (default off), each documented with its reasoning; no material-tuned constants:
     - Estimator: a spin-stiffness term for the magnetization residual from the local XC spin kernel, E_m ~ 1/2 integral of (d2 e_xc / dm2) d_m^2 (cheap, from the same LDA-PW evaluator), added to or replacing the Coulomb-metric spin term; justify against the onepath3 reasoning (the Coulomb metric under-weights short-wavelength spin errors).
     - Broyden: a general conditioning rule: when the history Gram matrix's reciprocal condition number falls below a machine-precision-derived limit, drop the oldest pairs until it is well conditioned (record every drop), instead of rejecting the run. This changes the existing singular-history contract: write the rule as an ADR draft (docs/decisions/0045-broyden-history-conditioning.md, status proposed) and keep the old rejection when the switch is off.
     - Fixed-moment occupations: check the per-channel number solve tolerance and the weighted-occupation residual definition against legacy; fix if the reference path solves less tightly than legacy.
  3. Re-run the full G4 synthetic acceptance matrix (all 48 comparisons) and the Si2 FD cases with the candidate switches; report which combination passes all bounds. If one does, make REFSCF_SPIN_FD_MOMENT route to the reference driver only together with the passing switches (keep defaults off; the orchestrator flips after dev08 confirms).
OWNED_FILES:
  - crates/dft-engine/src/refscf/electronic.rs
  - crates/dft-engine/src/refscf/electronic_map.rs
  - crates/dft-engine/src/refscf/broyden.rs
  - crates/dft-engine/src/refscf/metric.rs
  - crates/dft-engine/tests/refscf.rs   # g4/g4b tests only, appended
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Spin and Fermi-Dirac" section only
  - docs/decisions/0045-broyden-history-conditioning.md   # new, status proposed
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test --release -p dft-engine --test refscf g4 -- --include-ignored --test-threads=1
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: diagnosis documented; a switch combination that passes all 48 comparisons under unchanged bounds, or a documented negative result with the measured reasons; closed-shell and qualified G4 paths bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
