TASK_ID: P0B-COMPLEX-ABI
GOAL: Verify Apple Accelerate complex-valued BLAS/LAPACK calls from Rust FFI, as required by the P0 ticket but not covered yet (only real ddot/dsyev were tested).
OWNED_FILES: crates/dft-numerics/** , docs/environment.md (append a section only). Do not touch other crates, Cargo.toml workspace manifest, or AGENTS.md.
READ_ONLY_CONTEXT: AGENTS.md, docs/environment.md, docs/architecture.md, docs/decisions/0001-numerical-backend-isolation.md, crates/dft-numerics/src
CONSTRAINTS: No network. Zero external crates (std only; define a tiny #[repr(C)] complex struct locally, or reuse the existing one if any). unsafe only inside the existing private FFI adapter, with safety preconditions in comments and tests.
IMPLEMENTATION:
1. Add FFI for complex Hermitian eigensolve (zheev, jobz=V, uplo=U) and complex matrix multiply (cblas_zgemm or zgemm_). Use the same integer width (LP64 / ILP64) discipline as the existing real calls and state which one Accelerate's symbols you link against in docs/environment.md.
2. Tests with analytic answers, all tolerances stated:
   a. 2x2 Hermitian [[2, 1-i],[1+i, 3]] (column-major, interleaved re/im): eigenvalues exactly 1 and 4. Check residual ||Hv - lambda v||, orthonormality V^H V = I, INFO = 0.
   b. 4x4 complex Hermitian built from a fixed non-trivial unitary U and known diagonal D (A = U D U^H); recover D to 1e-12.
   c. zgemm: small A*B with complex entries compared against a naive Rust triple loop, and one case with conjugate-transpose (ConjTrans) on A.
3. Document in docs/environment.md: complex memory layout (interleaved, column-major), leading dimension rules, workspace query for zheev (lwork=-1), and any ABI pitfalls found (e.g. deprecation warnings, LP64 int).
TEST_COMMANDS (PATH needs /opt/homebrew/opt/rustup/bin): cargo build --workspace; cargo test --workspace; cargo fmt --check; cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE: all four commands exit 0; new tests print measured residuals/orthogonality errors (dimensionless) with --nocapture; no new unsafe outside the adapter.
GIT: The sandbox cannot write git metadata for the worktree. Leave changes uncommitted in the worktree and say so in the report; Claude will commit after review.
RETURN_FORMAT: report schema. Do not report an untested command as passing.
