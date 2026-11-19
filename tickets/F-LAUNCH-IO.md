TASK_ID: F-LAUNCH-IO
GOAL: |
  dev08 time budget (master, si64 MPI np=8, 6c31696 defaults, profile on): external wall B = 8.8835 s, covered_wall 8.3449 s (phase tree closes to 0.14 %), so 0.5386 s (6 % of the run, 14 % of F = all time outside ordinary fine maps) is spent OUTSIDE the engine phases: process launch, MPI_Init, input and pseudopotential parsing, preflight, result/output writing, teardown, binary hashing (the new `run` block hashes the binary at the end, about 5 ms), and rank synchronization at start and end. The other code spends its own startup in its wall too, so every 0.1 s here counts against the 0.90 target (at si64 np=8 the budget 0.90 Q is 4.16 s). TASK: find where those 0.54 s go and remove what can be removed without changing results. (1) Instrument the launcher path with the existing phase/clock machinery or simple wall timestamps printed under BAANDI_PHASE_PROFILE=on in a new `launch` section of results.json (never inside hot loops): time from process start (use the monotonic clock at main entry, and report the gap to the first engine phase), MPI init/finalize, argument and TOML parsing, UPF read + parse, `Config::prepare` (the SCF agent found `run` reuses preflight's prepare result: verify no second preparation happens, and that preflight does not run a heavy check twice on every rank), pseudopotential data transforms done before the engine starts (radial integrals, interpolation tables) that every rank repeats (distribute or share them if replicated and measurable), writing results.json, density.bin and any other outputs (output format conversions, large JSON serialization of histories, synchronous flush, fsync), run-block binary SHA-256 (move off the critical path or compute lazily/only when requested if it exceeds 5 ms), barriers/gathers at the end, and process teardown (freeing large structures explicitly before exit can cost; consider std::process::exit after flush when everything is written, only if safe for MPI finalize ordering). (2) Reduce: for each item > 10 ms propose and implement a fix (parallelize across ranks only where bit-identical; skip redundant work; avoid duplicate preparation; buffered I/O; write outputs on rank 0 only with one gather already available). Output files and results.json content must stay identical (except the new `launch` section when profiling is on); energies and densities bit-identical. (3) Mac: measure with a small case in release mode (`cargo build --release -p dft-cli`) at several rank counts with simulated or real mpirun if available (`mpirun` may not exist on the Mac: then threads/1 rank only) and report the launch section. The real number comes from dev08 (the orchestrator will run `BAANDI_PHASE_PROFILE=on mpirun -np 8 ... si64` and read `launch` + B - covered_wall).
OWNED_FILES:
  - crates/dft-cli/src/main.rs
  - crates/dft-cli/src/run.rs
  - crates/dft-cli/src/run_config.rs
  - crates/dft-cli/src/ (output writers only if a measured item needs it)
  - crates/dft-engine/src/profile.rs   # only to add the launch section/keys
  - crates/dft-physics/src/ (only if a pre-engine table build is replicated and measurable)
  - docs/performance/attempts-ledger.md  # rows F-LAUNCH-IO-*
READ_ONLY_CONTEXT: AGENTS.md, crates/dft-cli/AGENTS.md (if present), docs/performance/run-config.md, docs/performance/fast-loop-protocol.md, docs/performance/plan-2026-11-18-digest.md
IMPLEMENTATION_STEPS:
  1. Add the `launch` timing section (guarded by the existing profile switch, zero cost when off) and run a small case to see the split on the Mac.
  2. Fix the largest items first; keep each fix its own change with its own measurement note.
  3. Run scripts/check.sh fast dft-cli dft-engine dft-physics; keep every existing test green (do not weaken tests); new tests for any changed output path (results.json identical with the new section only when profiling is on).
TEST_COMMANDS:
  - scripts/check.sh fast dft-cli dft-engine
ACCEPTANCE_CRITERIA: outputs and energies bit-identical; launch section present when profiling; every item over 10 ms has either a fix or a documented reason; ledger rows with Mac numbers; dev08 confirmation commands in the report.
NON_GOALS: SCF algorithm changes, engine phase timers.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
