TASK_ID: P3-SETUP
GOAL: |
  Cut the time from process start to the first SCF map of the reference driver. dev08 MPI8 (main 2aee083, median of 3): si128 4.865 s vs QE init_run 0.750 s (wfcinit 0.680); si64 1.091 s vs 0.190 s. Today the setup is not split: results.json init_lane counters are all zero, so first measure, then cut. Also: output before serialization 0.772 s and "unaccounted" 0.624 s (si128), energy evaluation 1.295 s per run (QE has no separate line).
READ FIRST: docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md, docs/performance/backlog.md, docs/design/reference-scf.md ("Start"), docs/research/qe-study/qe-scf-cycle.md and qe-stage-by-stage.md. Never read or request QE source.
TASK: |
  1. Instrument the setup path from main() to the first map with lightweight phase clocks recorded at every diagnostic level (input/UPF parse and hashing, cell/basis/G-vectors, FFT plans, distribution/owned rows, structure factors and local/nonlocal tables, atomic superposition density, atomic-orbital start block and its first orthonormalization/diagonalization, MPI collectives) and serialize them into results.json (reference_scf.setup_s object). Same for the post-convergence output path (gather, density/orbital writes, hashing, JSON).
  2. With the Mac's small cases (debug build, observation only) find the dominant items and remove obvious excess that does not change numbers: replicated work that every rank does in full, repeated hashing of the same bytes, tables built for all G on every rank, O(N^2) loops, avoidable gathers, writes the run did not ask for. Changes that alter rounding go behind a switch (default off) for dev08 A/B; overhead-only changes proven bit-identical by a test may be default on.
  3. Mac checks only (user rule 2026-11-24): cargo check --workspace --tests; cargo test -p dft-engine --test refscf; cargo test -p dft-cli reference; cargo test -p dft-numerics switches; fmt; clippy -D warnings on touched crates.
  4. Write docs/performance/lanes/RSCF/setup-dev08-request.md: si64/si128 MPI8 (same inputs as the verify-cost A/B), one profiled run per arm to get the new setup table, then A/A + crossed A/B for each switch; report the setup table next to QE's init_run lines from the existing matched QE runs.
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs
  - crates/dft-engine/src/refscf/driver.rs   # setup/output clocks and setup-path code only; P3-DAV-DENSE does not touch this file
  - crates/dft-engine/src/refscf/record.rs   # setup clock fields only
  - crates/dft-cli/src/run.rs   # reference-driver setup/output path and reference_scf JSON only
  - crates/dft-cli/src/main.rs   # only if the setup path starts there
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - docs/design/reference-scf.md   # "Setup and output cost" section only
  - docs/performance/lanes/RSCF/setup-dev08-request.md
  - crates/dft-numerics/src/switches.rs   # ONLY if a switch is needed; coordinate: P3-DAV-DENSE also adds switches, so append yours at the end of the list
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-cli reference
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-cli --all-targets -- -D warnings
  - (orchestrator, on dev08) scripts/check.sh full, A/B
ACCEPTANCE_CRITERIA: setup and output split recorded at every diagnostic level; excess removed with identical science; switchable where rounding changes; dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
