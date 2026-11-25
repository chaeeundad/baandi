TASK_ID: G6B-GRID-SYMMETRY
GOAL: |
  Two findings from dev08 gen8 (G6 symmetry runs on Si2 diamond, Si.pz-vbc, ecut 30):
  1. An explicit basis.fft_grid = [28,28,28] with symmetry on passes preflight but the run uses [32,32,32] (results basis.fft_grid=[32,32,32], fft_grid_overridden=true); nothing reports the replacement. An explicit input value must be honoured exactly or rejected with a clear message; it must never be silently replaced (25^3 is rejected correctly with "symmetry.fft_grid: explicit grid does not preserve translations/rotations"). Find where the grid is changed (config/prepare/symmetry/reference driver), and make the rule one place: explicit grid -> use it if every density operation permutes it, else reject; automatic grid -> choose the smallest symmetric, FFT-friendly grid and record why. Check both drivers (legacy and reference) and the recorded fft_grid_overridden meaning (it should mean "user-specified").
  2. The reference per-map symmetrization residual is flat (5.6e-2 e/bohr^1.5 from first to last map) because it measures ||rho - P rho|| of the irreducible-k sum before projection, which is not expected to vanish. Redefine the record to something informative and document it: report both the pre-projection norm (expected nonzero, normalized by ||rho||) and the post-projection invariance residual ||P(P rho) - P rho|| (expected at rounding), or replace with what legacy records if it has an equivalent.
READ FIRST: crates/dft-cli/src/config.rs (prepare, symmetry), crates/dft-cli/src/preflight.rs (symmetry grid warning), crates/dft-basis/src/lib.rs (FftGrid override), crates/dft-engine/src/refscf/symmetry.rs, docs/design/reference-scf.md "Symmetry". Never read QE source.
TASK: Fix both; tests: explicit 28^3 diamond is used exactly (or rejected with the reason if it does not permute; diamond fractional translations are quarter-cell, 28 is divisible by 4), explicit 25^3 still rejected, automatic grid unchanged for existing inputs (results bit-identical); residual records on the diamond fixture: pre-projection normalized value nonzero, post-projection at rounding. Update docs.
OWNED_FILES:
  - crates/dft-cli/src/config.rs   # grid selection with symmetry only
  - crates/dft-cli/src/preflight.rs   # grid messages only
  - crates/dft-basis/src/lib.rs   # only if the override flag semantics need it
  - crates/dft-engine/src/refscf/symmetry.rs
  - crates/dft-cli/src/run.rs   # symmetry record JSON only
  - crates/dft-cli/tests/   # new tests only
  - crates/dft-engine/tests/refscf.rs   # new tests appended
  - docs/design/reference-scf.md   # "Symmetry" section only
  - docs/cli.md   # fft_grid with symmetry paragraph only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-cli
  - cargo test -p dft-engine --test refscf symmetry
  - cargo fmt --check && cargo clippy -p dft-cli -p dft-engine -p dft-basis --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: explicit grids honoured or rejected, never replaced silently; residual records meaningful; existing automatic-grid results unchanged.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
