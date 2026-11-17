TASK_ID: LEAN-ALGEBRA-HARDEN
GOAL: |
  An independent Codex review of the MPI generalized Davidson (BAANDI_M2_LEAN_ALGEBRA, default on since 2026-11-17) found three P1 hazards and weak tests. Reproduce each with a FAILING test first, then fix it without changing results on ordinary inputs (bit-identical decisions and bits for the existing fixtures where no hazard triggers) and without adding a collective on the hot path.
  Findings (read .delegate/REVIEW-LEAN-ALGEBRA/report.md for the full text):
  (F1) The pivot ratio test (eps/1e-12) in extend_factor (crates/dft-solvers/src/lean_davidson.rs ~666, distributed_davidson/workspace/lean.rs ~627) does not bound the conditioning of the whole overlap; a block v0=e0, vj=sqrt(0.99) e_{j-1}+0.1 e_j passes every pivot but Q = V T^T has Gram error 8.5e-4 after 8 columns. The design claim 'returned orthogonality ~1e-12' is not enforced; no Gram verification at return.
  (F2) Rank-local error returns (apply_block error, non-finite local output check, local initial row-count error: distributed_davidson.rs ~1413, workspace.rs ~781) return without global agreement while other ranks proceed into the next border reduction: deadlock hazard. (Inherited from the common path.)
  (F3) A factorable but partial Gram start (w < bands) is accepted as if orthonormal; random completion uses I - V V^T where the basis is not orthonormal (needs I - V (V^T V)^-1 V^T); after rebuild only 'zero columns' is checked, so 0 < m < bands can index out of range (workspace.rs ~878, ~909; lean.rs ~388).
  (F4) Tests give false confidence: the residual check compares reported vs actual residual but never bounds the residual; convergence flag is Ritz-change only; fallback comparison reuses the same algebra. Add an independent dense-eigenvalue oracle, an upper bound of the true residual of ALL required bands, non-orthogonal/partial starts, empty local partition, a forced stage-2 after contraction.
OWNED_FILES:
  - crates/dft-solvers/src/distributed_davidson/workspace/lean.rs
  - crates/dft-solvers/src/distributed_davidson/workspace.rs
  - crates/dft-solvers/src/distributed_davidson.rs
  - crates/dft-solvers/src/lean_davidson.rs        # shared extend_factor: serial results must stay bit-identical
  - crates/dft-solvers/tests/distributed_lean_algebra.rs
  - crates/dft-solvers/tests/lean_harden.rs        # new
  - docs/design/mpi-overlap-davidson.md            # correct the orthogonality claim and document the new checks
  - docs/performance/attempts-ledger.md            # one row
READ_ONLY_CONTEXT: AGENTS.md, crates/dft-solvers/AGENTS.md, docs/design/mpi-overlap-davidson.md, .delegate/REVIEW-LEAN-ALGEBRA/report.md
IMPLEMENTATION_STEPS:
  1. Write tests that fail on main for F1 (conditioning example above, serial and 1-4 simulated ranks), F2 (2 ranks: operator completes its collective, rank 1 returns an error or NaN: every rank must return the same error instead of hanging; use a test timeout), F3 (n=4, bands=2, one Gram-start column 1e40 e0 / a partial start: result must equal the normalised start; no panic; rebuild with 0<m<bands handled), F4 (oracle + residual bounds).
  2. F1: add a conditioning guard that does not add a collective: bound the growth of the inverse factor from the quantities already computed (e.g. track the max |T| entry or the pivot minimum and compare with a threshold derived from eps), and verify the Gram defect of the returned/contracted block using data already reduced (or one defect value folded into an existing reduction payload) and fall back to the existing stage 1/2 rules when exceeded. Document the threshold derivation.
  3. F2: carry a local error flag as one extra double inside an EXISTING reduction payload of the next collective on every rank (no new collective on the normal path); all ranks decide together and return the same error. Same for the local initial row-count check at entry (it may use the first reduction).
  4. F3: orthonormalize/normalize the Gram start before use (or project with the Gram inverse), reject/repair partial starts consistently on all ranks, make rebuild require m == bands or refill deterministically.
  5. Run the ordinary fixtures: results and counters bit-identical to main (compare energies/eigenvalue bits/counters recorded from main BEFORE your change); report any decision that changed.
TEST_COMMANDS:
  - cargo test -p dft-solvers --offline
  - scripts/check.sh dft-solvers
ACCEPTANCE_CRITERIA: the new tests fail on main and pass with the fix; ordinary fixtures bit-identical (serial lean and distributed) and reduction/collective counts per solve unchanged on the normal path; no tolerance weakened; clippy/fmt clean.
NON_GOALS: new algorithms, changing capacity rules, performance work.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
