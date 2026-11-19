# Slow tests and the speed of `scripts/check.sh`

Rule (ADR 0004, AGENTS.md): a test that takes more than about 10 s in a debug
build must be `#[ignore]` with a documented release-mode command, and the gate must
stay fast. This file is the single list of such tests and records how the gate
time was measured. Update it whenever a test is ignored for its runtime or the
gate changes.

## Gate speed under load (GATE-SPEED, Mac, 2026-11-19)

Problem: after big merges `scripts/check.sh` took 40 to 50 minutes on the shared Mac (log: 83
binaries, 2800 s test wall with 4 jobs; "slowest" `dft-cli pw92_prefactor_config` 164 s,
`dft-engine m2_entry_cost` 164 s, `dft-cli gamma_real_config` 160 s).

Findings (commit 74d62df, 12 cores, 16 GB RAM):
- The three named tests are not slow. Every binary run alone with libtest `--report-time
  --test-threads=1` (load 12 to 30): `m2_entry_cost` 0.04 s (its only test is already
  `#[ignore]`: the binary runs zero tests), `pw92_prefactor_config` 0.42 s, `gamma_real_config`
  0.52 s. All 83 binaries one after another: 29 s, every test passed; no single test above 5 s,
  largest binary `dft-cli test all` 6.9 s. Tests ignored for runtime by this change: 0.
- In the 2800 s log the median binary took about 140 s, including binaries that run no test, so
  the time was not spent in tests. Two machine-wide causes were found: (1) the dft-upf `parser`
  test spawned `cargo run -p dft-cli`, a second cold build of the CLI in `target/check-nested`
  after every merge (145 s of 299 s in a cold gate measured here, plus 2.7 GB per worktree),
  running at the same time as clippy's own cold build and the test binaries; (2) memory: 16 GB
  RAM with 9.1 of 10 GB swap in use and 7.3 M page-outs while many sessions compile, so every
  process start and allocation waits on swap. A binary measured under that pressure is slow
  without being a slow test, so per-test timings are taken alone (method above), never from a
  loaded gate log.
- Fix for (1): the test takes the CLI from `DFT_CLI_BIN`, which `scripts/check.sh` sets to the
  `dft-cli` binary of its build gate when the gate builds dft-cli (the whole workspace, or a
  fast set that contains dft-cli, which every dft-upf change selects). `dft-upf test all`
  went from 145 s to 0.50 s; no `check-nested` directory is created. Plain
  `cargo test -p dft-upf` still works (it falls back to `cargo run`).
- For (2) the gate runs its tests under `nice -n 10` (`CHECK_NICE`, default on) and prints the
  load and the number of running cargo/rustc processes on its first line, so a slow gate can be
  read as load. The remedy for the 40-minute case is fewer concurrent compiles, not fewer tests:
  the `fast` tier on the Mac and the `full` gate on dev08
  (`docs/performance/requests/gate-on-dev08.md`).

Gate walls with the new script (this worktree, warm target, load 13 to 25 from other sessions):

| tier, change | crates | binaries / tests passed | test-run wall | gate wall | load (1 min) |
|---|---|---|---|---|---|
| fast, edit in dft-cli | 1 | 13 / 130 | 3 s | 8 s | 15 |
| fast, interface edit in dft-engine | 4 | 45 / 370 | 7 s | 14 s | 14 |
| fast, interface edit in dft-numerics | 7 | 73 / 679 | 13 s | 25 s | 16 |
| full, nothing changed | 13 | 83 / 732 | 9 s | 13 s | 16 |
| full, interface edit in dft-numerics (merge-like) | 13 | 83 / 732 | 15 s | 29 s | 19 to 25 |
| full, cold (empty target: a new worktree or a rebuild of everything) | 13 | 83 / 732 | 21 s | 191 s | 18 to 36 |
| old script, cold (same commit, before this change) | 13 | 83 / all passed | 145 s in dft-upf alone | 299 s | 2 rising to 98 |

An interface edit is a new public item in the crate, which changes its metadata so every
dependent recompiles. A large merge that rewrites many crates sits between the
merge-like row and the cold row. The 40-minute gates of the user's log came from load
and memory pressure on top of that cold cost, not from the test count: every row
above runs the full test set of its crates.

### Shared build cache: evaluated, rejected for CARGO_TARGET_DIR (2026-11-19)

Experiment: two exports of 74d62df at different paths (as two worktrees), built with
`cargo test --workspace --no-run` into one `CARGO_TARGET_DIR` on the 990-2TB volume
(load 33 to 48 from other sessions).

| step | wall | compiled |
|---|---|---|
| A cold into the empty shared target | 114 s | everything |
| B (fresh mtimes, as a new checkout) into the warm shared target | 35 s | the 13 workspace crates; registry deps reused |
| A, B, A again, nothing edited | 0 to 1 s | nothing: no eviction when the content is equal |
| A and B both edit `dft-numerics` differently, built at the same time | A 48 s, B 47 s | A 9 crates; B **blocked on A's lock, then reused A's `dft-numerics`** |

The last row is a correctness failure: cargo names workspace crates by their
workspace-relative path and judges freshness by mtime, so in a shared target B's
test binaries were linked against A's code and B's own edit was never compiled. A
gate in one worktree could pass on another agent's changes. The lock also
serializes concurrent agent builds. A shared `CARGO_TARGET_DIR` across worktrees is
therefore not used. The safe form of the 114 s to 35 s gain (registry dependencies
compiled once for all worktrees) is a content-addressed compiler cache (sccache as
`RUSTC_WRAPPER`, cache directory on the 990-2TB volume). It is not installed on the
Mac and installing it is a machine change for the user to approve; expected effect
about 80 s off every cold first build, no effect on incremental rebuilds (sccache
does not cache incremental workspace crates).

Disk: the cold target of the new gate is 7.4 GB (no `check-nested`; was about 10 GB
with it). Agent worktree targets under `.claude/worktrees` held 465 GB on 2026-11-19
(40 targets, 8.5 GB of them `check-nested`, 7.1 GB `check-clippy`);
`scripts/wt_cleanup.py --dry-run` lists 577 GB (allocated blocks) removable from
merged, clean, idle worktrees.

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

None (re-measured 2026-11-19: no test above 5 s at opt-level 1 when run alone). The
tests that were already `#[ignore]` (175 at 74d62df: release benchmarks, external
files, dev08-only) keep their own documented commands in their `#[ignore = "..."]`
reason or doc comment.

Rule for a future slow test: measure it alone (`--exact <name> --test-threads=1`,
or the whole binary with `RUSTC_BOOTSTRAP=1 <binary> -Z unstable-options
--report-time --test-threads=1`), never from a loaded gate log. Above about 10 s
alone: `#[ignore = "slow: <release command>"]`, a row here with who runs it (the
dev08 gate, `requests/gate-on-dev08.md` step 4, unless it is the only guard of a
recently changed hot path; then it stays in `full` and the row says so).

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
  (2026-11-19: `check-nested` is no longer created when the gate builds dft-cli, see
  above; `scripts/wt_cleanup.py --check-dirs [--dry-run]` removes both side
  directories in every idle worktree.)
