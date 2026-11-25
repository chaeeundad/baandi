TASK_ID: R6C-NO-ATOMIC-START
GOAL: |
  Row L11 of docs/roadmap/refactor-r-plan.md: pseudopotentials without PP_RHOATOM or PP_PSWFC (or with l > 3 atomic orbitals) currently route to legacy because the reference start requires atomic charge and orbitals. Give the reference start a documented fallback so these files run on the reference driver, then remove L11.
READ FIRST: crates/dft-engine/src/refstart.rs (atomic start, the deterministic kinetic-ordered completion of missing bands), legacy start_density.rs / atomic_start.rs under src/oracle (what legacy does without atomic data: uniform density? superposition? random waves), docs/decisions/0039-atomic-start-wavefunctions.md; crates/dft-cli/src/driver_select.rs l11; scripts/gen_synthetic_upf.py.
TASK: |
  1. Start density without PP_RHOATOM: a neutral uniform valence density (or the same rule legacy uses; state which and why); start orbitals without PP_PSWFC or with l > 3 channels: the existing deterministic kinetic-ordered plane-wave completion for all bands (or only for the missing channels). Record the start kind in the start report and results.
  2. Remove L11 and the explicit rejection; keep a clear error only for truly unusable files (no local potential, etc.).
  3. Tests: synthetic UPFs without PP_RHOATOM, without PP_PSWFC, and with an l=4 channel: reference converges and agrees with the strict oracle within the standard class (and 1e-10 Ha tight); files with atomic data bit-identical (reference_identity, bench_scf both baselines).
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs
  - crates/dft-cli/src/driver_select.rs   # L11 only
  - crates/dft-cli/tests/   # routing expectations and new tests only
  - crates/dft-engine/tests/refscf.rs   # new tests appended
  - scripts/gen_synthetic_upf.py   # new fixture variants only
  - docs/design/reference-scf.md   # start section paragraph only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - scripts/check.sh (CHECK_SCCACHE=0, BAANDI_TEST_TMP outside the checkout)
  - python3 -I scripts/bench_scf.py --require-identical (both baselines)
  - cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: UPFs without atomic data run on the reference driver with oracle agreement; L11 removed; others bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
