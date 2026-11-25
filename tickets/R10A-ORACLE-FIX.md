TASK_ID: R10A-ORACLE-FIX
GOAL: |
  dev08 gen20 RG-ORC found that R10A-DELETE-EXPERIMENTS (or a concurrent merge: R6A, R8C, R9A) changed the strict oracle on the Gamma path: si8 cubic Gamma, Si.pz-vbc, ecut 30, grid 36^3, BAANDI_REFERENCE_STRICT=on: energy bits c03f5261b8862aca -> ...2ad4 (dE -3.6e-14 Ha), density differs, SCF iterations 11 -> 9, final-map Davidson iterations 16 -> 7. The displaced si8 forces case also changed (dF 1.8e-7). k-point and FD oracle cases are unchanged. The oracle must be bit-identical to before (user decision D2), so find and restore the arithmetic.
READ FIRST: the R10A commit (git show --stat for "R10A-DELETE-EXPERIMENTS" and "Merge R10A-DELETE-EXPERIMENTS"), docs/performance/attempts-ledger.md rows "retired in R10a", crates/dft-engine/src/{scf.rs,session.rs,gamma_real.rs,cutoff_ladder.rs} and the legacy Gamma strict path; the R6A, R8C and R9A merges only to rule them out.
TASK: |
  1. Reproduce locally: build main before R10A's merge (the parent commit of "Merge R10A-DELETE-EXPERIMENTS" on main) and main HEAD; run a strict si8 Gamma case (or a smaller Gamma case that shows the same change: the bench gamma case did not change, so build a fixture close to si8: cubic 8-atom Si, 36^3, at the Si.pz-vbc in ../pseudo) with BAANDI_REFERENCE_STRICT=on; compare energy bits, iterations and Davidson counts.
  2. Bisect within the R10A diff: which deleted default-on behaviour (a switch whose default was on, a first-block/start step, a ladder stage on Gamma, a Davidson variant) changed the strict Gamma path. Restore exactly that default arithmetic in the oracle (hard-coded at its former default, no switch), keeping the deletions otherwise.
  3. Add a debug-fast regression test (synthetic Gamma fixture that reproduces the divergence if possible) that pins strict-oracle energy bits and iteration counts, so a future deletion cannot change it unnoticed; record the si8 bits from dev08 in the test comments.
OWNED_FILES:
  - crates/dft-engine/src/   # legacy/oracle files only (not refscf/ or refstart.rs)
  - crates/dft-engine/tests/oracle_identity.rs   # new
  - crates/dft-numerics/src/switches.rs   # only if a removed default must return (prefer hard-coding)
  - docs/performance/attempts-ledger.md   # one row describing the fix
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test oracle_identity
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (both baselines)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: strict si8 Gamma bits equal the pre-R10A main locally; cause documented; regression test added.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
