TASK_ID: G4D-FREE-MOMENT-PBE
GOAL: |
  Qualify free-moment collinear Fermi-Dirac with PBE on the reference driver (today routed to legacy, routing predicate L1 / G4D_FALLBACK). G4C measured: no-core near-full-polarization synthetic case at standard threshold |dF| 9.2e-7 Ha vs the oracle while the fixed-point verification passes (M = 0.99994 vs 1.0; minority-occupation error 2.8e-5 e); with NLCC it passes; tight threshold passes. The G4B spin-stiffness estimator term uses the local LDA-PW spin kernel as a proxy even for PBE (docs/design/reference-scf.md "Spin and Fermi-Dirac", G4C/G4D diagnosis record with per-map estimator terms and true errors). Acceptance bounds unchanged (5e-9 Ha standard, 1e-10 tight, 1e-4 mu_B, 1e-4 Ha mu, 1e-4 e).
READ FIRST: docs/design/reference-scf.md (Spin and Fermi-Dirac incl. the G4D diagnosis record), docs/decisions/0045 and 0046; crates/dft-engine/src/refscf/{electronic.rs,electronic_map.rs} (stiffness term, G4F screening of it), crates/dft-physics/src/gga/distributed.rs (SpinPbe); the g4d_free_moment_pbe_acceptance test. Never read QE source.
TASK: |
  1. Use the diagnosis record (and extend it if needed) to decide why the estimator stops early: (a) the LDA-PW kernel proxy under-weights the PBE spin response near full polarization; (b) the Fermi-level/occupation response of the minority channel near zeta -> 1 is not seen by the estimator; (c) something else (state it).
  2. Implement the fix that follows from the diagnosis, as a general rule (no material- or polarization-specific constant): e.g. a stiffness term from finite differences of the actual functional's spin potentials (PBE when PBE), and/or an occupation-response term for FD. Keep it default-on only if the full G4/G4B/G4C/G4D matrices pass with unchanged bounds; scalar FD, fixed spin, fixed moment and LDA free moment fingerprints must stay bit-identical where the new term is not active.
  3. If it qualifies, remove L1/G4D_FALLBACK (crates/dft-cli/src/driver_select.rs) so free-moment PBE runs on the reference driver; else write the negative result with numbers and keep the fallback.
OWNED_FILES:
  - crates/dft-engine/src/refscf/electronic.rs
  - crates/dft-engine/src/refscf/electronic_map.rs
  - crates/dft-engine/src/refscf/metric.rs
  - crates/dft-physics/src/gga/distributed.rs   # response evaluation only
  - crates/dft-engine/tests/refscf.rs   # g4d tests only
  - crates/dft-cli/src/driver_select.rs   # L1 only
  - crates/dft-cli/tests/   # routing expectations only
  - docs/design/reference-scf.md   # "Spin and Fermi-Dirac" section only
  - docs/decisions/0045-broyden-history-conditioning.md   # note only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test --release -p dft-engine --test refscf g4 -- --include-ignored --test-threads=1
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: free-moment PBE qualified with unchanged bounds and L1 removed, or a documented negative result with the fallback kept.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
