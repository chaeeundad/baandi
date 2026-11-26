TASK_ID: G4D-G4G-RESEARCH
GOAL: |
  Continue the two open research items with the strict oracle as the reference, now that the per-map diagnostics exist (docs/design/reference-scf.md G4D/G4G sections; docs/decisions/0045, 0046, 0047).
  G4G: isolated O atom, fixed collinear sector 4/2 (HGH-like analytic projectors, PBE, vacuum box): the reference driver stalls (dr2 ~1.9e-8 Ha floor, history condition up to 1e14, stall records at maps 68/135 on dev08) while the strict oracle (legacy Pulay, same input) converges in 20 iterations and keeps the axial degeneracy exact (up 2p -0.383429 x2, down LUMO pair exactly degenerate); the reference run splits those pairs by 3e-5/7e-5 Ha and ends 5.4e-5 Ha below the oracle energy without converging. At the start the minority 2p channel is exactly degenerate (gap 3e-15 Ha). 39-70% of the stalled residual lies in the occupied-orbital rotation span (synthetic).
  G4D: free-moment collinear FD with PBE: the reference estimator stops early (PBE spin response 25x the LDA proxy; using the true response fixes the energy but |dmu| 1.6e-3 > 1e-4); on the real O FD case the strict oracle also diverges (growth limit at iterations 13-16).
READ FIRST: docs/design/reference-scf.md (Spin and Fermi-Dirac, G4D/G4G records), crates/dft-engine/src/refscf/{electronic.rs,electronic_map.rs,broyden.rs,stall.rs,metric.rs}, refstart.rs (spin start), crates/dft-engine/src/oracle/ (how the legacy collinear start and Pulay keep the symmetric solution; read only, the oracle is frozen); crates/dft-engine/tests/refscf.rs g4d/g4g tests. Never read QE source.
TASK: |
  G4G first:
  1. Compare reference vs oracle on the synthetic G4G atom map by map from the start: where does the reference first break the axial degeneracy (start orbitals, first Davidson solve, first Broyden step)? The oracle uses a symmetric atomic start and Pulay; the reference uses its own start and Broyden. Test (a) a start that preserves exact degeneracy of degenerate atomic orbitals (orbitals as complete degenerate shells, deterministic ordering), (b) Broyden vs a Pulay-equivalent update on the same trajectory, (c) Davidson degenerate-subspace handling (locking/ordering of exactly degenerate Ritz pairs, rotation within a degenerate occupied/empty boundary).
  2. Fix with a general rule (e.g. degenerate-shell-preserving start; or Ritz vectors of an exactly degenerate cluster rotated to a canonical basis each solve; whichever the data supports), behind a switch default off until dev08; acceptance: the synthetic atom converges and agrees with the oracle within the standard class; all baselines and G4 matrices bit-identical when off.
  G4D next:
  3. With the strict oracle diverging too, first establish a converged reference for the G4D synthetic/real case: tighter mixing on the oracle (smaller beta) or a smaller smearing; record what converges. Then test whether the G4G fix (if any) or the true PBE spin response together with a Fermi-level response term qualifies free-moment PBE under unchanged bounds.
  Report per-map evidence in docs/design/reference-scf.md (G4D/G4G sections), keep every bound unchanged.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-engine/src/refstart.rs
  - crates/dft-solvers/src/ref_davidson_complex.rs   # degenerate-cluster handling only
  - crates/dft-solvers/src/ref_davidson.rs   # degenerate-cluster handling only
  - crates/dft-solvers/src/ref_davidson/
  - crates/dft-numerics/src/switches.rs   # append new switches only
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/   # switch counts only
  - crates/dft-engine/tests/refscf.rs   # g4d/g4g tests appended
  - docs/design/reference-scf.md   # G4D/G4G sections only
  - docs/decisions/0045-broyden-history-conditioning.md   # notes only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test --release -p dft-engine --test refscf g4 -- --include-ignored --test-threads=1
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP=/Volumes/990-2TB/Develop/rust-dft-wt/tmp-<task>, never /private/tmp)
  - python3 -I scripts/bench_scf.py --require-identical (both baselines)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: G4G cause located map by map with a switchable general fix that converges the synthetic atom to the oracle; G4D converged reference established and qualification attempted; everything else bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
