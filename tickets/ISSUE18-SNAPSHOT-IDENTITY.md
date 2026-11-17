TASK_ID: ISSUE18-SNAPSHOT-IDENTITY
GOAL: |
  Implement GitHub issue #18 (full text: tickets/issues/ISSUE-18.md): the measurement-only replay harness (davidson_bench --snapshot, ScfProblem::replay_solve, kept H X load) must validate operator and row-layout identity BEFORE injecting captured products. Separate exact replay (identity gate exits non-zero on mismatch) from an explicitly labeled algorithm A/B on the same physical problem. Record a versioned operator fingerprint in the snapshot manifest (cell, positions, PP content hashes, XC/spin/occupation, cutoff and ordered reciprocal basis, density and potential fingerprints, convention/layout), compare each rank's ordered row mapping with the new operator's, restore or verify tolerance policy / criterion / start state, and turn a rank-local mismatch into an agreed collective error before any rank enters a divergent collective sequence.
OWNED_FILES:
  - crates/dft-cli/examples/davidson_bench.rs
  - crates/dft-engine/src/snapshot.rs      # manifest/fingerprint fields (ISSUE17 adds other metadata fields: keep additions separate and mergeable)
  - crates/dft-engine/src/scf.rs           # replay_solve only
  - crates/dft-engine/src/gamma_real.rs    # kept-product load part only; ISSUE17 edits the diagnostic plumbing elsewhere in the file
  - crates/dft-engine/tests/snapshot_identity.rs   # new
  - docs/performance/fast-loop-protocol.md  # short section on exact replay vs A/B
  - docs/performance/attempts-ledger.md     # one row
READ_ONLY_CONTEXT: AGENTS.md, crates/dft-engine/AGENTS.md, tickets/issues/ISSUE-18.md, docs/performance/fast-loop-protocol.md
IMPLEMENTATION_STEPS:
  1. Write failing tests first: reject a snapshot after (a) moved atom, (b) different same-valence PP, (c) changed XC, (d) changed local row order/ownership, (e) modified density bytes, (f) incompatible stopping profile, while ranks/grid/bands are unchanged; the original snapshot must still reproduce. Add a recompute/invalidation test for a permitted rounding-changing kernel A/B.
  2. Implement fingerprint + checks; old manifests without a fingerprint are accepted only in an explicitly labeled legacy mode with a printed warning.
  3. Keep the new operator-failure agreement guarantees (no rank-local early return before a collective).
TEST_COMMANDS:
  - cargo test -p dft-engine --offline
  - cargo build -p dft-cli --examples --offline
  - scripts/check.sh dft-engine
ACCEPTANCE_CRITERIA: all six rejection cases and the original-reproduces test pass; no change to normal SCF behaviour; clippy/fmt clean.
NON_GOALS: restart-correctness of ordinary SCF; performance work.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
