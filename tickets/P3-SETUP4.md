TASK_ID: P3-SETUP4
GOAL: |
  Reference-driver start is 2.30 s on dev08 si128 MPI8 (back-to-back P3s runs at e5ec8ae, verify off) against QE init_run 0.78 s (wfcinit 0.68 s). After REFSTART_MKL, the REUSE switches, REFSTART_GENERALIZED and COMPACT_G, the remaining si128 items (P3l profile) are about: atomic_h_apply 0.64, atomic_h_gram 0.73 (H applied to and projected on all atomic orbitals: 512 columns for 260 bands), density/field FFT plan creation 0.29 + 0.17, local tables + structure factors 0.39, generalized start EVD and rotation ~0.3, basis 0.10.
READ FIRST: docs/performance/AGENT-BRIEF.md; attempts-ledger.md rows REFSTART_MKL, P3-SETUP, P3-SETUP2, P3-SETUP3, P3l, P3i; docs/design/reference-scf.md "Start" and "Setup and output cost"; docs/research/qe-study/qe-scf-cycle.md (QE's starting wavefunctions: H on all atomic wavefunctions then one subspace diagonalization, so the column count matches ours). Never read QE source.
TASK: |
  1. FFT plans: create each distinct FFT plan once per run and share it between the start, the density/field transforms and the H application (the earlier ticket could not because the plan API is in crates/dft-solvers/src/distributed_fft.rs and crates/dft-engine/src/refstart.rs; both are owned here). Prefer plan reuse over re-creation; check rustfft/MKL planner caching. Bit-identical results expected.
  2. Atomic block H application: it applies the full H (local FFT + nonlocal) to all atomic orbitals once. Check that it uses the same tuned paths as the SCF H application (batched pairs, NONLOCAL_SINGLE_GEMM/BACK_ACCUM, wide GEMM routes); route it through them if not. The Gram (S and H borders of the atomic block) must use the single-call MKL route.
  3. Local tables and structure factors: compute only for the owned G rows, vectorize, avoid recomputation across species.
  4. Behind switches where rounding changes (default off); bit-identical changes may be proposed default on with a test. Mac checks only: cargo check --workspace --tests; cargo test -p dft-engine --test refscf; cargo test -p dft-solvers; cargo test -p dft-numerics switches; fmt; clippy -D warnings. Before proposing any default flip, grep the tests for the switch name.
  5. docs/performance/lanes/RSCF/setup4-dev08-request.md: ONE profiled run per arm, four cases, setup_s split next to QE init_run.
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs
  - crates/dft-engine/src/refscf/   # plan sharing only
  - crates/dft-solvers/src/distributed_fft.rs   # plan creation/sharing API only
  - crates/dft-physics/src/distributed_real_gamma.rs   # only to reuse its H-apply path for the atomic block
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Setup and output cost" section only
  - docs/performance/lanes/RSCF/setup4-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-solvers
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-solvers -p dft-physics -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: plans created once and shared; atomic block uses the tuned H/GEMM paths; start state equal (or within rounding behind a switch); dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
