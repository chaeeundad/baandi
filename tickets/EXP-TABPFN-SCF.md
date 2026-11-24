TASK_ID: EXP-TABPFN-SCF
KIND: experiment branch (user request 2026-11-25). Never merged into main without a separate decision. Mac-only, small cells.
GOAL: |
  Experiment: an "active" SCF where an external learned controller (TabPFN, run by the orchestrator in Python) observes the reference SCF trajectory map by map and chooses the next mixing action, and we measure how the physics (delivered fixed point, energy, density, verification) and the computation (maps, H columns, Davidson iterations, wall incl. controller overhead) change against the fixed rule (modified Broyden beta 0.7, history 8).
  Your part is the Rust side only: a controller hook, an exploration mode to generate training data, and complete per-map trajectory logs. No ML code in Rust.
TASK: |
  1. Per-map trajectory record (reference driver, every diagnostic level, experiment switch EXP_SCF_TRACE=on, default off): one JSON line per map written by the CLI layer (not inside numerical loops) to <output-dir>/scf_trace.jsonl with: iteration, dr2 (Ry), ethr used, beta used, history occupancy, Broyden diagnostics already in Diagnostics (gamma norms, condition if available), Davidson iterations and H columns of the map, max/mean Ritz change of the last correction, band-energy, Hartree, xc and total energy terms and their change from the previous map, descf, the stopping-estimator decision, verification attempt/skip/error if any, map wall seconds. Keep units explicit in key names.
  2. Controller hook (switch EXP_SCF_CONTROLLER = path to an executable, default empty = off): the driver spawns it once at start; after each map it writes the same JSON line to the child's stdin and reads one JSON line back: {"beta": f64 in (0, 1.5], "reset_history": bool}. The driver applies beta to this step's Broyden update (and clears the history when asked); every other rule (ethr ratchet, stopping estimator, verification criterion and REFSCF_VERIFY_GATE) stays exactly as is. Invalid or missing replies fall back to the fixed rule and are recorded. Controller wall time is measured separately and recorded per map. Never deliver without the unchanged passing verification.
  3. Exploration mode (EXP_SCF_BETA_EXPLORE = seed, default off): beta per map drawn uniformly from [0.2, 1.0] with a deterministic seeded generator (record seed and draws), to create training trajectories.
  4. Tests on the Mac (debug, small): hook round trip with a tiny test controller script (python3 -I, standard library only) that echoes beta 0.7 must reproduce the default trajectory bit-identically; a controller returning garbage falls back; exploration is deterministic per seed. Defaults bit-identical with all switches off.
  5. Write docs/research/exp-tabpfn-scf.md: the hook protocol, the trace schema, and the planned experiment (cells: Si8 displaced, Si16, rskew, elongated 3x1x1, the synthetic nonlocal fixture; how to generate training data with exploration seeds; evaluation metrics: maps, H columns, wall with and without controller time, E - E_strict, density difference to the default run's delivered density, verification attempts).
OWNED_FILES:
  - crates/dft-engine/src/refscf/
  - crates/dft-cli/src/run.rs   # reference-driver trace writing and controller process handling only
  - crates/dft-numerics/src/switches.rs   # append EXP_* switches at the end
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-cli/tests/exp_scf_controller.rs   # new
  - crates/dft-cli/tests/fixtures/exp_echo_controller.py   # new
  - docs/research/exp-tabpfn-scf.md
TEST_COMMANDS:   # Mac rule: no check.sh full, no release matrices
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli --test exp_scf_controller
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: hook, trace and exploration behind default-off switches; echo controller bit-identical to default; fallback tested; design doc written.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
