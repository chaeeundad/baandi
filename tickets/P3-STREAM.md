TASK_ID: P3-STREAM
GOAL: |
  The reference Davidson's non-GEMM streaming passes cost ~3.4 s on dev08 si128 MPI8 (one profiled run at 0f11a240, rank max, 20 solves) while doing almost no arithmetic: precondition 1.21 s for 0.70 GFLOP (54 M hypot calls, scalar bound), Davidson setup 0.96 s streaming 3.5 GB, normalize 0.31 s (0.2 GFLOP/s), append 0.27 s, apply_check 0.21 s, plus residual formation inside corrections. QE's g_psi + update on the same cell is ~1.6 s total. Make these passes memory-speed or fuse them.
READ FIRST: docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md (P3i row and earlier P3 rows), docs/design/reference-scf.md "Inner solver cost accounting" (incl. the per-stage work counters added by P3-DAVREM), docs/research/qe-study/davidson-iteration-structure-gamma.md (preconditioner addendum). Never read QE source.
TASK: |
  1. Precondition: the denominator d = (1 + x + sqrt(1 + (x-1)^2)) / 2, x = h(G) - e s(G) (Ry), is evaluated per (row, band) with hypot. Replace hypot by an exact-enough sqrt form (document the rounding difference bound), vectorize over rows (contiguous column, chunked loops the compiler can auto-vectorize, no per-element branches), and fuse the residual formation HX - eX, the division by d and the norm accumulation into one pass per column.
  2. Fuse normalize/append/apply_check passes with their producers where the data is already in cache; remove separate full-matrix passes that only scale or copy; avoid zero fills of buffers about to be overwritten.
  3. Each change behind its own switch (default off; bit-identical ones may be proposed default on; the hypot replacement is rounding class). Band decisions and H-column counts unchanged on fixtures; record any rounding branch.
  4. Mac checks only: cargo check --workspace --tests; cargo test -p dft-solvers; cargo test -p dft-engine --test refscf; cargo test -p dft-numerics switches; fmt; clippy -D warnings.
  5. docs/performance/lanes/RSCF/stream-dev08-request.md: ONE profiled run per arm on si64/si128/elongated/rskew MPI8 with the per-stage work table.
OWNED_FILES:
  - crates/dft-solvers/src/ref_davidson.rs
  - crates/dft-solvers/src/ref_davidson/
  - crates/dft-solvers/src/tests.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Inner solver cost accounting" section only
  - docs/performance/lanes/RSCF/stream-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-solvers
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-solvers -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: precondition and streaming passes vectorized/fused behind switches; decisions and H columns unchanged on fixtures; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
