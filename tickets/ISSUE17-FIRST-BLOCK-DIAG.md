TASK_ID: ISSUE17-FIRST-BLOCK-DIAG
GOAL: |
  Implement GitHub issue #17 (full text: tickets/issues/ISSUE-17.md): the first-block / residual diagnostic and snapshot metadata use strict_bands (= n_occ + 1 on the standard path) as the occupied count and write electrons = 2*policy_occupied; it also compares residual norms with a Ritz-change tolerance. Fix: carry actual occupation count / electron count (and occupation vector where integer) from ScfProblem into the Gamma diagnostic and capture context; keep mandatory roots, density-carrying roots and empty buffers as separate fields; carry the convergence-criterion kind; report residuals and Ritz changes separately and rename/remove counts implying a standard-path verdict from a residual comparison; compute the occupied residual spectrum over exactly the occupied block with the true HOMO/LUMO boundary. Production SCF numbers (energy/density bits, solver work) must not change; diagnostics on/off identical apart from explicitly separated diagnostic work.
OWNED_FILES:
  - crates/dft-engine/src/first_block.rs
  - crates/dft-engine/src/gamma_real.rs   # only the diagnostic/capture context plumbing near policy_occupied; another task (ISSUE18) edits the snapshot replay load part
  - crates/dft-engine/src/scf.rs          # only the call sites that build the diagnostic/capture context
  - crates/dft-engine/src/snapshot.rs     # metadata fields (add, keep old readable)
  - crates/dft-engine/tests/first_block_diag.rs   # new
  - docs/performance/first-block-residual-2026-11-17.md   # mark superseded, add a short correction note; do not rewrite history
  - docs/performance/attempts-ledger.md   # one row
READ_ONLY_CONTEXT: AGENTS.md, crates/dft-engine/AGENTS.md, tickets/issues/ISSUE-17.md, docs/performance/fast-loop-protocol.md
IMPLEMENTATION_STEPS:
  1. Read the code at the lines quoted in the issue (they refer to commit 1c52998; use current main). Write the acceptance fixtures as failing tests first: n_occ=2, strict_bands=3, bands=5 -> 4 electrons and gap between indices 1 and 2; a dominant first empty residual must not affect the occupied-only spectrum; separate strict-residual vs standard-Ritz expectations including the last mandatory root.
  2. Implement the plumbing and field split; keep old snapshot manifests readable.
  3. Verify diagnostics on/off leave energy/density bits and solver work counters unchanged (existing fixtures, short fixed-step run).
  4. Do NOT recompute the dev08 tables (the orchestrator does that); only mark the old report as superseded with a pointer to this change.
TEST_COMMANDS:
  - cargo test -p dft-engine --offline
  - scripts/check.sh dft-engine
ACCEPTANCE_CRITERIA: all issue acceptance checks as tests; production bits unchanged; no tolerance relaxed; clippy/fmt clean.
NON_GOALS: deciding whether first-block compression is worth an experiment; performance work.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
