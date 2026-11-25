TASK_ID: G4G-CLI-RECORDS
GOAL: |
  dev08 gen30: the G4G diagnostics (per-channel gaps, edge residuals, occupied-subspace overlaps, stall records: 3 consecutive maps with history condition > 1e12 and non-decreasing estimate) are computed by the engine but never written by the CLI. On the O HGH fixed-sector run the predicate holds at maps 67 and 134, yet results.json and stderr show nothing. Also the documented workaround needs a correction: Fermi-Dirac 0.005 Ha with fixed moment 2.0 converges (14 maps, verified) but to a different, fractional spherical state (down 2p 1/3 each), 0.0107 Ha above the strict integer-sector energy; it is not a cure for the integer-sector stall.
READ FIRST: crates/dft-engine/src/refscf/ (stall.rs, map report fields stall/channel_gap_ha, the electronic and k drivers' records); crates/dft-cli/src/run.rs (reference_scf record serialization); docs/design/reference-scf.md and docs/decisions/0045 (G4G text).
TASK: |
  1. Serialize per map in results.json reference_scf.iterations[]: channel_gap_ha (per channel), and when present the stall record; add a top-level reference_scf.stalls list; print one stderr line per stall event naming map, condition, estimate and the recommendation.
  2. Correct the recommendation text (engine message and docs): for an open-shell atom whose integer sector stalls, Fermi-Dirac with fixed moment converges to the fractional (spherical) state, which is a different state; for the integer broken-symmetry state use the strict oracle (BAANDI_REFERENCE_STRICT=on) until G4G is resolved.
  3. Tests: a CLI test on the synthetic stalling fixture (short max_iterations) that checks the stall record and gap fields are present; other results unchanged (reference_identity, bench_scf both baselines).
OWNED_FILES:
  - crates/dft-cli/src/run.rs   # reference_scf record fields only
  - crates/dft-cli/tests/   # new test only
  - crates/dft-engine/src/refscf/   # message text and record accessors only
  - docs/design/reference-scf.md   # G4G paragraph only
  - docs/decisions/0045-broyden-history-conditioning.md   # G4G note only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-cli
  - python3 -I scripts/bench_scf.py --require-identical (both baselines)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: stall and gap diagnostics visible in results.json and stderr; corrected recommendation; results unchanged.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
