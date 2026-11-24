TASK_ID: P3-RESET
GOAL: |
  At each restart the reference Davidson rotates V and HV and then recomputes a fresh S and H Gram of the contracted basis, factors S and transforms H ("measured refresh"): about 4 r n^2 FLOPs plus a factor and two collectives per restart; si128 MPI8 has 26 restarts per run (restart stage 1.7 s, part of border/factor). After a Rayleigh-Ritz rotation by orthonormal Ritz vectors the contracted basis is, in exact arithmetic, S = I and H = diag(Ritz values); QE's Davidson relies on that. The ideal reset was rejected earlier (docs/design/reference-scf.md "Inner solver cost accounting": a degenerate fixture drifted 1.0154e-11 Ha) under a bit-identity/strict-only premise. Adoption criteria now (user rules 2026-11-10 and 2026-11-25): strict class 1e-10 Ha and standard 5e-9 Ha, verification, convergence and counts; 1e-11 Ha is inside both classes. Re-evaluate with this premise.
READ FIRST: docs/performance/AGENT-BRIEF.md, attempts-ledger.md (P3 rows, the reset rejection), docs/design/reference-scf.md "Inner solver cost accounting", docs/research/qe-study/davidson-iteration-structure-gamma.md. Never read QE source.
TASK: |
  1. Implement REF_DAV_IDEAL_RESET (default off): after a successful contraction set S = I and H = diag(theta) instead of measuring them; keep a cheap guard that detects loss of orthogonality (e.g. a diagonal-only S estimate from the norms already available, or a periodic full refresh every k restarts) and falls back to the measured refresh when the guard trips; record the fallbacks.
  2. Tests (Mac only): fixtures incl. the degenerate one: decisions and H columns, final energies against the measured-refresh run (report the max difference; must stay below 1e-10 Ha), orthogonality error after each restart, rank counts 1/3/8 simulated.
  3. docs/performance/lanes/RSCF/reset-dev08-request.md: ONE profiled run per arm on si64/si128/elongated/rskew MPI8 with restart/border/factor stage seconds, counts and E - E_strict.
OWNED_FILES:
  - crates/dft-solvers/src/ref_davidson.rs
  - crates/dft-solvers/src/ref_davidson/
  - crates/dft-solvers/src/tests.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end (P3-ADDV also appends)
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Inner solver cost accounting" section only
  - docs/performance/lanes/RSCF/reset-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: switch default off with guard and fallback; fixture energies within 1e-10 Ha of the measured refresh; counts reported; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
