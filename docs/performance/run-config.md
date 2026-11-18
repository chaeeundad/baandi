# Run-owned configuration and validated run report (RUN-CONFIG / RUN-REPORT, 2026-11-18)

Work order: one run-owned configuration instead of process-global switch reads, a test harness that
exploits it, and one validated report. Code: `crates/dft-numerics/src/switches.rs`,
`crates/dft-cli/src/run_config.rs`, `crates/dft-engine/src/{scf,session,cutoff_ladder,profile}.rs`,
`scripts/parse_results.py`. Test targets: `docs/performance/test-targets.md`.

## Configuration

- Sources and precedence: environment (`BAANDI_<NAME>`) > input file > registry default. Input tables:
  `[performance]` (numerical and path switches) and `[profile]` (`phase_profile`, `stage_timers`,
  `davidson_stages`, `snapshot_dir`, `phase_wait`, `first_block`, `dgemm_timers`). Key = registry name
  without `BAANDI_`, lower case; booleans mean `on`/`off`; numbers are their text.
  `cutoff_ladder`, `relax_extrapolation` and `pseudo_dir` are environment only (the first two override
  input-file keys that already exist).
  ```toml
  [performance]
  m2_grid = "on"
  stop_metric = "fft"
  [profile]
  phase_profile = true
  ```
- Environment values keep the legacy rules and error texts exactly (an invalid value falls back as
  before and is recorded as `invalid -> default`/`off`/`error`). Input values are strict: unknown key,
  wrong table, environment-only key or a value outside the allowed set is exit 1 before any work.
- `RunConfig` is resolved once per run and immutable (values); it carries the run's bookkeeping:
  consulted flags (frozen after the SCF so the result writer's own reads do not count), path notes,
  run ID, run-local counters, agreement hash.
- Process and scope: `dft-cli` installs the input's configuration as the process configuration before
  the first switch read (main), so every thread of a CLI run sees it. In-process callers attach a
  configuration per problem (`ScfProblem::set_run_config`) or per thread (`RunConfig::enter`); sessions
  enter their problem's configuration while they initialize and step, and `WorkerPool` carries the
  scope into its workers (every parallel region of the workspace goes through the pool).
- MPI: before any optional collective every rank compares the FNV-1a hash of the effective values of
  every switch except `PSEUDO_DIR`, `SNAPSHOT_DIR`, `DGEMM_TIMERS` (one allgatherv of three doubles per
  rank plus one broadcast of rank 0's run ID). A mismatch is the same error and exit 1 on every rank.

## Migrated (resolved once from the run's configuration, passed explicitly)

| where | switches | resolved at |
|---|---|---|
| `session::SessionSwitches` (was one OnceLock per switch per process) | `STANDARD_STOP`, `STANDARD_RESIDUAL_GUARD` (closed-shell session and verification), `IMPORT_CHECK`, `VERIFY_REUSE`, `CERT_SCREEN`, `POTENTIAL_REUSE`, `MIX_GROWTH_RESET`, `EIGEN_FACTOR` | problem construction, `set_run_config` |
| `ScfProblem` construction | `STOP_METRIC` (also re-resolved by `set_run_config`), `EWALD_CACHE` | construction |
| `CutoffLadderOptions::from_run` | `CUTOFF_LADDER_STANDARD`, `LADDER_STANDARD_SEED`, `LADDER_SECANT`, `LADDER_COARSE_STOP`, `LADDER_COARSE_RULE`, `LADDER_FIRST_TOL` | ladder options; the coarse problem is built and run under the fine problem's configuration |
| run-scoped instead of process-cached (one read per construction, step or audit, through the scope) | `DENSITY_PAIR`, `DENSITY_REDUCE`, `NL_ENERGY`, `ENERGY_IDENTITY`, `EMPTY_LOOSE`, `BUFFER_TOL` | `BandDensity::new`, each audit, each step |
| `dft-cli` | every switch (input tables), `MPI_TREE`, `M2_CAPACITY`, `REFERENCE_STRICT`, `CUTOFF_LADDER` | process start |

Run-local counters (`RunCounter`): stop-rule estimates kept / transformed, Ewald hits / computed,
problems, session steps, verification maps. results.json `run.counters`.

Demonstration: `crates/dft-engine/tests/all/run_config.rs` runs two SCF problems with opposite
`STOP_METRIC` concurrently in one process, without a lock: each is bit-identical to its solo run
(energy, density, eigenvalues, steps), the default run equals the process-configuration run, the
counters are independent (kept estimates only in the kept run, transformed only in the fft run).
The test runs inside the shared engine harness (`--test all run_config::`).

## Remaining on process-global state (list)

Switches still cached in a per-process `OnceLock` (the first read in a process wins; a second run with
a different value in the same process does not see it; correct for the CLI, one run per process):
- `crates/dft-engine/src/gamma_real.rs` (owned by the concurrent M2 work, not refactored here):
  `REUSE_PSI`, `ATOMIC_REAL`, `NL_ENERGY_REAL`, `M2_REUSE`, `M2_LAZY`, `M2_LIGHT_STATE`, `M2_LAZY_SERIAL`;
- `advance.rs`: `SCF_ADVANCE`, `SCF_ADVANCE_BLOCKS`, `SCF_ADVANCE_GROWTH`, `SCF_ADVANCE_ADAPT`,
  `SCF_ADVANCE_SCOPE`;
- measurement: `PHASE_PROFILE` (`profile::enabled`, read on hot paths), `STAGE_TIMERS`
  (`stage_clock`), `FIRST_BLOCK`, `SNAPSHOT_DIR`, `PHASE_WAIT`, `M2_LAYOUT_CACHE`
  (`distributed_davidson.rs`).

Switches read per call through the legacy accessors (correct per run under the scope, but looked up
in deep code rather than handed down): the solver kernels (`dft-solvers`: `B1*`, `B5`, `B6`,
`DAVIDSON`, `ORTHO*`, `MGS`, `WIDE_DOT`, `RITZ_FUSED`, `LEAN_*`, `RESTART_KEEP`, `M2_*` of the
distributed solver, `DAVIDSON_STAGES`), physics (`MKL_NONLOCAL`, `APPLY_LEAN`, `HAPPLY_*`,
`NONLOCAL_BLOCK`, `EWALD_FAST`, `IONIC_FAST`), numerics (`PURE_GEMM`, `PACK_FIXED`, `ROTATE_KERNEL`,
`MKL_*`, `PAIR_FFT`), the M2 engine files (`grid_dist.rs` `M2_GRID`, `mix_dist.rs` `M2_MIX`, `m2.rs`
`EngineTuning`, `parallel.rs` `MPI_M2`), the collinear and Fermi-Dirac drivers
(`STANDARD_RESIDUAL_GUARD`, `SPIN_SHARE`), `REFERENCE_STRICT` in `ScfOptions::from_reference_switch`.

Thread-local test overrides (`set_*_override`, `profile::force`, `stage_clock::force`, ...) keep their
precedence over the configuration. Counters other than the seven run-local ones stay process-global
atomics (`lean_counters`, `reuse_psi_counters`, `grid_dist::counters`, `mix_dist::counters`,
`m2::stats`, `work_timers`, `stage_clock`, `gemm_dispatch`, `phase_clock`, the profile registry).
The results writer derives switch paths from them (one run per CLI process).

Next migration steps (small, mechanical): `advance.rs` into `SessionSwitches`; `gamma_real.rs` and the
M2 files once their owners are done (resolve into `RealSolveState` / `GridPlan` / `EngineTuning` at
creation); solver switches into the Davidson workspace options at workspace creation; then the
corresponding test binaries can drop their locks and join the harness.

## Report

- `scf.switches.<name>`: `value`, `source` (`default`, `input`, `env`, `invalid -> ...`), `default`,
  `allowed`, `meaning`, plus `requested` (environment text, else input text), `consulted`, `path`,
  `fallback_reason` (recorded for `MPI_M2`, `M2_GRID`, `M2_MIX`, `STOP_METRIC`, `REUSE_PSI`,
  `EWALD_CACHE` from the counters after the run).
- `run`: `schema` `run-report-1`, `run_id`, `exit_reason` (`converged`, `not_converged`, `diverged`),
  `source_commit`, `binary.path` and `binary.sha256` (hashed after the timed run), `rank`, `ranks`,
  `switch_agreement_hash`, `counters`, and the completion marker `completed: true`. The outer
  `schema_version` stays 1.
- `phase_profile` (`phase-profile-2`): `per_rank` (each rank's own wall, the sum of its own bucket
  walls, residual) and `closure` (covered max/min over ranks, critical rank, imbalance, wall and
  residual max/min). `covered_wall_s` (sum over buckets of the per-bucket maximum, possibly from
  different ranks) is kept and labelled legacy. The external wall stays the verdict.
- `scripts/parse_results.py RESULTS [--manifest M | --binary B | --expect-sha H] [--json]`: exit 0
  ACCEPT, exit 2 REJECT with one line per problem: truncated/interrupted (invalid JSON or no
  completion marker), stale binary (SHA-256 differs from the manifest or missing), clocks missing or
  zero on an active measured path (profile switch on without a profile; with stage timers on: a
  bucket with maps and zero wall, zero covered time on every rank, Davidson stage timers on with every
  stage zero after iterations, MPI collectives with calls and zero seconds), impossible rows (min >
  max, negative `solve.other`/`other_*`, a rank's covered time above its own wall beyond 1 ms or
  0.1 %). Tests: `scripts/tests/test_parse_results.py` (synthetic, incl. simulated rank imbalance and
  an interrupted write), `crates/dft-cli/tests/all/run_config.rs` (real output: accept, truncated
  copy, stale SHA), `crates/dft-engine/tests/instrument.rs::per_rank_closure_through_the_gather_on_simulated_ranks`
  (the real-MPI gather path on simulated ranks), lib `profile::closure_tests`.
