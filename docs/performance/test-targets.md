# Test targets (TEST-HARNESS, 2026-11-18)

Every Rust integration-test file used to be its own test binary: each one links the whole dependency
tree (engine, solvers, physics, faer), so link time and target-directory size grew with the file count.
With the run-owned configuration (`docs/performance/run-config.md`) the files that hold no
process-global state are compiled into ONE harness binary per crate, `crates/<crate>/tests/all/main.rs`
with one module per former file (`tests/all/<name>.rs`). The files that still hold process-global
state stay separate; each is listed below with the boundary that justifies it.

Running one former target: `cargo test -p <crate> --test all <name>::` (one test:
`--test all <name>::<test>`; `--exact` takes the full path `<name>::<test>`). Live commands in the
crate notes, scripts and current docs were rewritten to this form; historical records (ledger,
optimization log, dev08 requests, tickets, reference evidence) keep the commands as they were run.
`scripts/check.sh` is unchanged (PASS/FAIL per gate; it discovers the binaries from cargo's
compiler-artifact records).

Merge rule (what made a file eligible): no environment mutation, no forced switch or per-thread
override hook or communicator registration, no read of process-global counters or registries, no
static lock, no global allocator, no `crate::` path of its own crate root (also through `#[path]`
includes). Relative `#[path]`, `include*!` and `mod x;` paths were rewritten for the new location.
The engine harness declares the shared fixtures once as `crate::common` (the checkpoint fixture
reaches them there). Classification: a scan of every file for these patterns, then hand checks
(global allocators, crate-root includes).

## Counts (compiler-artifact records of `cargo test --workspace --no-run --message-format=json`)

| | test executables | integration | lib | bin | cdylib |
|---|---:|---:|---:|---:|---:|
| before (1ddcc07) | 184 | 170 | 12 | 1 | 1 |
| after | 81 | 67 (7 harnesses + 60) | 12 | 1 | 1 |

Integration files: 172 (170 at 1ddcc07 plus the two RUN-CONFIG tests); 112 now live in the seven
harnesses (dft-cli 25, dft-engine 43, dft-numerics 12, dft-physics 14, dft-solvers 12, dft-units 2,
dft-upf 4), 60 stay separate.

## Gate walls (Mac, pure build, dev profile opt-level 1, `CHECK_TEST_JOBS=4`; observation O)

Measured 2026-11-18 on the shared Mac (load average 4 to 16 from other sessions), base = an export of
1ddcc07 with its own target directory, new = this branch; both warm before each measurement.

| | base (1ddcc07) | new | how |
|---|---:|---:|---|
| build + link of every test executable after touching `crates/dft-numerics/src/lib.rs` | 10.5, 11.0, 10.1 s | 6.0, 6.4, 6.4 s | `cargo test --workspace --no-run` |
| test-suite wall (all binaries, 4 jobs, as `scripts/check.sh`) | 19.4, 21.6, 23.1 s | 13.6, 14.2, 14.2 s | `scripts/run_test_binaries.py --jobs 4` |
| full `scripts/check.sh` (build, test, fmt, clippy) | 384 s (cold target directory) | 15 s (warm) | not comparable (cold vs warm); all four gates PASS on both |
| SCF wall outside the suite: `dft-cli run` si32 seed 601, T=1 / T=4, release | 5.03 / 2.19 s | 5.19 / 2.26 s | single runs (identity smoke); paired screen: ledger RUN-CONFIG-OVERHEAD, ratio 1.001 +- 2.9 % |

The SCF results are bit-identical between the two binaries (energy bits, density SHA-256, iteration
counts; si8 and si32, T=1 and T=4).

## Remaining separate integration binaries and their boundary

Library unit-test targets (12), the dft-cli binary's own test target and the baandi-py cdylib are
crate boundaries (one per crate, not mergeable across crates). Integration binaries:

| crate | target | boundary |
|---|---|---|
| baandi-mcp | `server` | only integration file of its crate |
| dft-cli | `atomic_start_checkpoint`, `checkpoint_inspect`, `checkpoint_io`, `gamma_real_config`, `pseudo_convention_config`, `reference_mode_restart` | crate root: include the engine checkpoint fixture, which needs the engine `common` as `crate::common` (the CLI harness root would hold the CLI `common`); no process-global state, mergeable once the fixture takes its pseudopotentials as a parameter |
| dft-cli | `forces`, `preflight_findings`, `pw92_prefactor_config`, `spin` | forced switch / per-thread override / communicator registration |
| dft-engine | `audit_shared_identity` | crate root: includes `src/energy.rs`, which needs the engine re-exports at `crate::` |
| dft-engine | `generalized_gram_stage0` | crate root: uses `crate::` paths of its own root |
| dft-engine | `capacity_trajectory`, `fp32_snapshots`, `hybrid_budget`, `preconditioner_study` | reads process-global counters or registries (deltas would mix with concurrent tests) |
| dft-engine | `init_scale`, `m2_coarse_certificate`, `m2_entry_cost`, `m2_kept_start`, `m2_lazy_orbitals`, `m2_lazy_result`, `m2_stage_timers`, `scf_advance`, `stage_timers` | forced switch / per-thread override / communicator registration; reads process-global counters |
| dft-engine | `distributed_scf`, `first_block_diag`, `instrument`, `m2_grid`, `m2_lazy_serial`, `m2_owned` (added on main after the count), `m2_remaining`, `other_init` | forced switch / override / communicator; process-global counters; serialized by a static lock |
| dft-engine | `m2_lean_algebra`, `snapshot_identity` | forced switch / override / communicator; serialized by a static lock |
| dft-engine | `spin_share` | mutates the process environment; forced switch / override |
| dft-geometry | `symmetry` | only integration file of its crate |
| dft-mkl | `in_place`, `linked` | feature boundary (crate-level `cfg(feature = "mkl")`) |
| dft-numerics | `gemm_dispatch`, `work_timers` | reads process-global counters |
| dft-numerics | `pair_fused_plan`, `pure_gemm_alloc` | counting global allocator (process-global) |
| dft-physics | `distributed_real_gamma`, `setup_memo_identity` | forced switch / override (Ewald / ionic scheme overrides, communicator) |
| dft-physics | `work_timers` | reads process-global counters |
| dft-solvers | `agreement_scratch`, `distributed_start_gram`, `lean_harden` | forced switch / override |
| dft-solvers | `algebra_stream`, `b5_reduced_eigen`, `lean_davidson`, `lean_proj_cache` | reads process-global counters |
| dft-solvers | `b1_narrow_margin`, `b1_selected_roots`, `distributed_lean_algebra`, `restart_keep` | reads process-global counters; serialized by a static lock |
| dft-solvers | `distributed_selected_roots` | serialized by a static lock |
| dft-solvers | `distributed_workspace` | counting global allocator (process-global) |
| dft-solvers | `generalized_gram_stage0` | crate root: uses `crate::` paths of its own root |

Joining the harness later: migrate the file's switch overrides to a per-call `RunConfig`
(`ScfProblem::set_run_config`) and its counters to run-local counters, remove the lock, then move the
file into `tests/all/` (rewrite relative paths, add `mod <name>;`).
