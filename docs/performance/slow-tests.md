# Slow tests and the speed of `scripts/check.sh`

Rule (ADR 0004, AGENTS.md): a test that takes more than about 10 s in a debug
build must be `#[ignore]` with a documented release-mode command, and the gate must
stay fast. This file is the single list of such tests and records how the gate
time was measured. Update it whenever a test is ignored for its runtime or the
gate changes.

## Measurement (Mac, 12 cores, 2026-11-18, shared machine)

Method: every test binary built by `cargo test --workspace --no-run` was run once,
niced, with libtest `--report-time` (184 binaries, 695 tests run, 174 already
ignored, 0 failed); the tests above 8 s were then rerun alone
(`--exact <name> --test-threads=1`). Cold gate = empty target directory, the cost
of a new worktree.

| Gate (whole `scripts/check.sh`) | Profile | Cold target | Warm target | Load average |
|---|---|---|---|---|
| Before (serial `cargo test`) | dev opt-level 0 | 458 s (build 29, test compile 29, test run 379, clippy 20) | about 400 s (test run dominates) | about 5 |
| Before, deps optimized only | deps opt 2 | 429 s (test run 309) | | 7 to 10 |
| Before, dev opt 1 + deps opt 2 | dev opt 1 | 256 s (build 80, test compile 52, test run 104, clippy 20) | | 8 to 14 |
| After (this change) | dev opt 1 + deps opt 2, parallel runner | 261 s | 12 s | 37 to 108 (other sessions) |

The after cold number was taken under a load average of 37 to 108; it is an upper
bound. The warm number is the common case (an agent re-running the gate after an
edit); incremental rebuild after touching `dft-numerics` was 2 s build and 9 s
test compile at opt 1 against 4 s and 14 s at opt 0.

The "3+ hours" observed before could not be reproduced on an idle-ish machine
(the old gate took 458 s cold here). Likely contributors, from the measurement:
the serial test run (379 s, multiplied by machine contention when several
worktrees run the gate at once) and a rebuild on every gate run: `dft-upf`'s
`parser` test spawns `cargo run -p dft-cli`, which resolves a different
dependency graph than the workspace build and, outside `cargo test`, the default
rather than the pinned toolchain; both invalidated `faer` and five workspace
crates in the shared target directory (about 40 s per run under load). Both are
fixed in the new gate.

## Tests above 10 s in the old debug profile (opt-level 0)

Measured alone (`--exact --test-threads=1`, niced). With the dev profile at
opt-level 1 (Cargo.toml, this change) all of them run in at most 3.2 s, so none
is ignored; they stay in the default gate. Count ignored for runtime by this
change: 0.

| Test | Old debug s | Debug s at opt 1 | Protects |
|---|---|---|---|
| `dft-solvers --test all distributed_davidson::serial_order_orthonormalization_matches_serial_eigenpairs ` | 26.1 | 2.2 | MPI Davidson orthonormalization order equals serial eigenpairs |
| `dft-physics --lib ewald::tests::serial_and_partitioned_ewald_match_the_previous_evaluator_bit_for_bit` | 21.2 | 3.2 | INIT-SCALE-3 rank-split Ewald sum, bit-identical (changed 2026-11-18) |
| `dft-solvers --test b1_selected_roots b1_results_do_not_depend_on_threads` | 16.6 | 0.8 | selected-roots path thread independence |
| `dft-cli --test all onepath_last_fallbacks::checkpoint_restart_carries_the_verification ` | 15.6 | 1.3 | checkpoint restart keeps the verification step |
| `dft-solvers --test lean_davidson percore_fused_border_and_ritz_loop` | 14.9 | 1.2 | per-core fused border and Ritz loop |
| `dft-solvers --test b1_selected_roots strict_profile_keeps_decisions_and_certifies_every_return` | 13.0 | 0.2 | strict profile decisions and certificate |
| `dft-solvers --test all distributed_davidson::bordered_reduced_update_stays_in_the_rounding_class ` | 12.2 | 0.8 | bordered reduced update rounding class |
| `dft-solvers --test all distributed_davidson::root_eigen_broadcast_is_bit_identical_to_replicated ` | 12.1 | 0.8 | root eigen broadcast bit identity |
| `dft-solvers --test algebra_stream rounding_class_switches_keep_roots_and_thread_identity` | 11.0 | 0.8 | rounding-class switches, thread identity |
| `dft-solvers --test algebra_stream ritz_fused_is_bit_identical` | 10.9 | 0.3 | fused Ritz bit identity |

Run one alone: `cargo test -p <crate> <--test name | --lib> -- --exact <test> --test-threads=1`.
The `dft-engine --test distributed_scf` tests looked slow in a parallel run (9 to
13 s) but take 0.2 to 1.8 s alone: they serialize on a lock, the parallel timer
counted the wait.

## Tests ignored for runtime

None added by this change. The 174 tests that were already `#[ignore]` (release
benchmarks, external files, dev08-only) keep their own documented commands in
their `#[ignore = "..."]` reason or doc comment.

## How the gate runs now

- Cargo.toml: `[profile.dev] opt-level = 1`, `[profile.dev.package."*"] opt-level = 2`.
  Release is unchanged. Debug assertions and overflow checks stay on. rustc does
  not contract or reorder floating-point operations at any opt-level; the same
  695 tests passed, including every bit-comparison and stored-value test.
- `scripts/check.sh`: build; `cargo test --no-run` once; the test binaries run
  concurrently through `scripts/run_test_binaries.py` (`CHECK_TEST_JOBS`, default 4;
  longest first from the timings cached in `target/check-test-times.json`;
  `CHECK_TEST_JOBS=1` runs them one at a time) while clippy runs in
  `target/check-clippy`; doctests through `cargo test --doc`; fmt. Tests that spawn
  cargo get `RUSTUP_TOOLCHAIN` (as under `cargo test`) and their own
  `target/check-nested` directory. Output: one PASS/FAIL line per gate (build, test,
  fmt, clippy), real exit code. Measured warm: 2 jobs 100 s, 4 jobs 93 s, 8 jobs
  91 s before the nested-build fix (all dominated by that rebuild); 4 jobs 12 s
  after it.
- Extra disk per target directory: about 2.6 GB (`check-nested`) and 0.9 GB
  (`check-clippy`); `scripts/wt_cleanup.py` removes them with `target/`.
