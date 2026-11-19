TASK_ID: F-INIT-STATES
GOAL: |
  dev08 (master, si64 MPI np=8, 6c31696 defaults) `init` = 0.879 s = 22 % of the time spent outside ordinary fine maps (3.93 s vs a total budget of 4.16 s = 0.90 x the other code). Earlier work (INIT-SCALE, INIT-SCALE-2, ledger rows INIT-*) distributed ionic_local, Ewald and the atomic-start orthonormalization and halved init, but the per-rank complex nonlocal projector table (the `states` bucket, ~100 MB per rank at si64, memory-bound) is still built in full on every rank, and `cli`, `basis`, `densities`, `rest` also still grow with rank count (replication + memory-bandwidth contention). Ledger INIT-NONLOCAL-TABLE explains why owned-row construction was NOT done: the split atomic-start Gram applies the complex H to column blocks over ALL G rows, so building only owned rows changes the arithmetic (rounding class) and the materialization audits also read the full table. TASK: remove the replicated complex projector table from init with a rounding-class switch BAANDI_INIT_OWNED_STATES=on|off (default OFF; add to crates/dft-numerics/src/switches.rs and update the switch-count tests crates/dft-numerics/tests/all/switches.rs and crates/dft-cli/tests/all/run.rs): (1) read the current construction (crates/dft-physics/src/hamiltonian.rs, nonlocal operator, distributed_real_gamma.rs owned rows, crates/dft-engine/src/scf.rs atomic start, gamma_real.rs) and list every consumer of the full complex table (atomic-start Gram/Rayleigh-Ritz, band energy/nonlocal audits, forces/stress, materialization, certificate screen, fallback paths); (2) with the switch on build only the rank-owned rows of the table (real-Gamma half-sphere rows owned by the rank) plus whatever the atomic start needs expressed as a row-distributed product P^H D P via a tree reduction (use fixed-order sums: rank-count-independent chunks), keep the full table lazily constructible for the consumers that need it (fallback, forces/stress, audits) so nothing breaks, record the lazy constructions in counters (`owned.states_full_builds`); (3) also look for avoidable replication in `cli`, `basis`, `densities`, `rest` (rows with +20..+33 % growth from np2 to np8 at si64): find what each does per rank and remove redundant per-rank work or share one construction per node/process, bit-identical where possible; (4) evidence: unit tests on simulated ranks 1..7 incl. skewed cells (energies in the rounding class, |dE| <= 1e-12 Ha vs the replicated table; atomic-start energies and SCF counts identical), fallbacks exercised, memory per rank counter (bytes of the table) before/after; Mac profile bucket directions with fast_screen_mpi.py (simulated ranks, observation only); dev08 confirmation commands in the report (si64 and si32 np=8, 10 steps, BAANDI_PHASE_PROFILE=on, A: default, B: BAANDI_INIT_OWNED_STATES=on, 3 alternating pairs, rows init, states, atomic, rest, cli, basis). Do not touch crates/dft-solvers. Use TMPDIR_990 under /Volumes/990-2TB/Develop/rust-dft-scratch (never /tmp or /private/tmp).
OWNED_FILES:
  - crates/dft-physics/src/hamiltonian.rs
  - crates/dft-physics/src/distributed_real_gamma.rs
  - crates/dft-engine/src/scf.rs
  - crates/dft-engine/src/gamma_real.rs
  - crates/dft-engine/src/m2.rs
  - crates/dft-numerics/src/switches.rs
  - crates/dft-numerics/tests/all/switches.rs
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - crates/dft-engine/tests/init_owned_states.rs   # new
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, crates/dft-engine/AGENTS.md, crates/dft-physics/AGENTS.md, ledger rows INIT-*, I20-*, docs/wiki/mpi.md
IMPLEMENTATION_STEPS:
  1. Inventory consumers (write it into the ledger row).
  2. Implement the owned table behind the switch; lazy full build for consumers.
  3. Tests (simulated ranks, fallbacks), then `scripts/check.sh fast dft-physics dft-engine dft-numerics dft-cli`.
TEST_COMMANDS:
  - scripts/check.sh fast dft-physics dft-engine dft-numerics dft-cli
ACCEPTANCE_CRITERIA: default off bit-identical; on: rounding class <= 1e-12 Ha and counters showing no full table build on the default production path; fallbacks tested; clippy/fmt clean.
NON_GOALS: SCF algorithm changes, solver changes.
GIT: leave changes uncommitted; Claude commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
