TASK_ID: P3-NL
GOAL: |
  Two remaining per-call gaps on dev08 si128 MPI8 (one run at fa055919, BAANDI verify off vs QE rmm-davidson, rank 0): (a) nonlocal projections calbec 0.0208 vs 0.0172 s per H call (+21%) and back-application add_vuspsi 0.0245 vs 0.0185 s per call (+32%), 83 calls; (b) BAANDI output 0.83 s after convergence (orbital_gather_unpack ~0.6 s, nested collectives ~0.3 s) that QE does not spend in its timed electrons loop.
READ FIRST: docs/performance/AGENT-BRIEF.md; docs/performance/attempts-ledger.md rows M2_DGEMM, M2-OWNED-FIELDS, MPI-LAZY-ORBITALS-FINDING, P3-HAPPLY, P3-QECLOCK, P3k; docs/research/qe-study/qe-hpsi-and-fft.md (calbec/add_vuspsi structure). Never read QE source.
TASK: |
  1. Nonlocal: from the P3-HAPPLY split and gemm_dispatch, give the GEMM shapes, routes and achieved GFLOP/s of the projection (P^T psi) and back-application (P W) per H call, plus the coupling (D matrix) step and the reduction. Remove excess behind a switch (default off): e.g. one GEMM over all atoms/projectors instead of per-species/per-atom blocks, fused coupling, the reduction merged with another collective, projector storage layout that suits the GEMM (column-major, contiguous rows), avoiding re-packing of projectors per call.
  2. Output: deliver only what the run needs: gather orbitals to rank 0 only when the user asked for orbital output (or lazily), write density from the owned layout without an extra gather where the format allows; keep results.json content identical. Behind a switch if any output changes; otherwise default on with a test that results.json and density bytes are identical.
  3. Mac checks only: cargo check --workspace --tests; cargo test -p dft-physics; cargo test -p dft-engine --test refscf; cargo test -p dft-cli reference; fmt; clippy -D warnings.
  4. docs/performance/lanes/RSCF/nl-dev08-request.md: ONE profiled run per arm on si64/si128/elongated/rskew MPI8 with calbec/add_vuspsi per call and output_s.
OWNED_FILES:
  - crates/dft-physics/src/distributed_real_gamma.rs   # nonlocal projection/back-application only
  - crates/dft-physics/src/   # nonlocal operator files reached from it only
  - crates/dft-engine/src/refscf/driver.rs   # final delivery/gather only
  - crates/dft-cli/src/run.rs   # reference-driver output path only
  - crates/dft-numerics/src/switches.rs   # append at the end (P3-SETUP3 also appends)
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/lanes/RSCF/nl-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-physics
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli reference
  - cargo fmt --check && cargo clippy -p dft-physics -p dft-engine -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: nonlocal per-call accounting and a switchable reduction; output path delivers identical results with less gather work; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
