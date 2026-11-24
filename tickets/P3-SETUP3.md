TASK_ID: P3-SETUP3
GOAL: |
  Reference-driver setup is still 2.75 s on dev08 si128 MPI8 vs QE init_run 0.76 s (one run at fa055919); si64 0.74 vs 0.19 s. After REFSTART_MKL and the two REUSE switches, the si128 items are about: atomic_h_apply 0.64, atomic_h_apply_gram 0.73, atomic_orthonormalize 0.53, QR 0.51, local tables + structure factors 0.39, density FFT plan 0.29, field plan 0.17, basis 0.13 (see the P3-SETUP2 setup_s split in results.json).
READ FIRST: docs/performance/AGENT-BRIEF.md; docs/performance/attempts-ledger.md rows REFSTART_MKL, P3-SETUP, P3-SETUP2, P3i; docs/design/reference-scf.md "Start" and "Setup and output cost"; docs/research/qe-study/qe-scf-cycle.md and qe-stage-by-stage.md (how QE builds its starting wavefunctions: atomic wavefunctions then a subspace diagonalization). Never read QE source.
TASK: |
  1. The atomic start applies H to all atomic orbitals, orthonormalizes them twice (CholQR2) and diagonalizes in that subspace. Account rows x columns for each step against the number of bands actually needed, and check whether the orthonormalization is needed at all before a generalized subspace eigensolve that takes the overlap matrix directly (one S and one H Gram, one generalized EVD, one rotation), which removes the QR passes.
  2. Implement the cheaper start behind a switch (default off; rounding class), keeping the start state equal within rounding and the reference trajectory counts unchanged on the fixtures (report any change in the first-map dr2).
  3. Remove remaining setup excess: FFT plan creation (reuse/cache within a run), structure-factor tables (only the G actually used), basis construction.
  4. Mac checks only: cargo check --workspace --tests; cargo test -p dft-engine --test refscf; cargo test -p dft-numerics switches; fmt; clippy -D warnings.
  5. docs/performance/lanes/RSCF/setup3-dev08-request.md: ONE profiled run per arm on si64/si128/elongated/rskew MPI8 with the setup_s split next to QE init_run/wfcinit.
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end (P3-NL also appends)
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Setup and output cost" section only
  - docs/performance/lanes/RSCF/setup3-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: accounting written; cheaper start behind a switch with start state within rounding and unchanged counts; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
