TASK_ID: P3-KEPTHX
SOURCE: user work order "BAANDI work order - reduce full-cost SCF work" (2026-11-26), part H; also add the two fast_order_ tests named there.
GOAL: |
  The reference Davidson returns X and HX, but the next ordinary map keeps only X and applies the full Hamiltonian to every starting column (crates/dft-engine/src/refscf/map.rs, crates/dft-solvers/src/ref_davidson.rs setup). On the older MPI path, reusing the starting H products removed 24.8-28.3% of full-H columns and reduced delivered wall by 1.2-5.1% (ledger M2-STICKS-START-REUSE-DEV08-ADOPT, M2_REUSE). On si128 the start block is ~29% of all H columns (260 per solve x 16 solves of 10932 verify-off).
TASK: |
  1. Opt-in switch REFSCF_KEPT_HX (default off). Retain the prior Ritz block X, its H products HX and the previous owned local potential. Use the existing band_density_sum_keeping and apply_local_change (crates/dft-physics/src/distributed_real_gamma.rs, the delta-potential kernel) with bounded slab storage. On a valid warm start, give Davidson H_old X + (V_new - V_old) X instead of its initial full-H application; keep the overlap, rank and convergence checks.
  2. Tie reuse to the exact retained orbital block, fixed operator, basis and row layout. The map currently passes a constant potential epoch of zero: track the actual old/new potentials for the update. Use a fresh H on the first map, cutoff/layout changes (including the REFSCF_LADDER handoff), retries and independent verification maps. All ranks must agree (one collective) before entering the reuse branch.
  3. One local test, named with the prefix fast_order_ in crates/dft-engine/tests/refscf.rs: a tiny two-rank synthetic nonlocal fixture; change the local potential, compare the reused starting product with a direct application (relative Frobenius error <= 1e-13 on a nonzero block) and assert exactly B fewer full-H starting columns for one B-band hit; then invalidate the retained block and assert the next entry takes the full-H path.
  4. Also add the ladder screen test with the same prefix (REFSCF_LADDER already exists on main): one tiny skewed Gamma cell, two simulated ranks: copied coefficients and zero-filled new modes bitwise before any fine-stage reorthogonalization, Gamma conjugate pairs, charge error <= 1e-12 e, max change of an orbital Gram entry <= 1e-12; continue to verified fine delivery and compare with the direct fine solve (<= 1e-10 Ha with tight settings). Reuse the existing ladder test helpers if present.
  5. Run `cargo test -p dft-engine --test refscf fast_order_` and report the test count (must be 2). Mac checks otherwise as usual (cargo check, touched crates, fmt, clippy).
  6. Record per run in results.json (reference_scf): kept-HX cache hits, full-H starting columns replaced, delta-V columns, retained slab bytes per rank.
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-solvers/src/ref_davidson.rs   # accept a supplied initial HX only
  - crates/dft-solvers/src/ref_davidson/
  - crates/dft-physics/src/distributed_real_gamma.rs   # reuse of the delta-potential kernel and kept slabs only
  - crates/dft-engine/tests/refscf.rs
  - crates/dft-cli/src/run.rs   # kept-HX counters JSON only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # a "Kept HX warm start" section only
TEST_COMMANDS:
  - cargo test -p dft-engine --test refscf fast_order_
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-solvers
  - cargo test -p dft-physics
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-solvers -p dft-physics -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: REFSCF_KEPT_HX default off; fast_order_ tests 2/2 pass with the stated tolerances; counters recorded.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
