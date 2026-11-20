# LANE-GLUE Phase 0 report

Status: **PROFILE_REQUEST** — stop here pending the master's dev08 relay.
Base: `61822bca3b470df1212ee47b47c7eb4e786e22e3`. Changes uncommitted.
No optimization designed or implemented. No Mac performance conclusions.

## Changes

- Default-off `GLUE_PROFILE`, additionally gated by `PHASE_PROFILE`, cached
  in grid plans and serial/distributed Pulay mixers. 23 entered-call spans,
  14 work/size counters (60 added fixed run-counter slots including clocks).
  `STAGE_TIMERS=off` makes the new spans count-only. Existing run-counter
  serialization supplies rank-local data without CLI/profile writer changes.
- Existing per-map envelopes and FFT/phase diagnostics complement the new
  cumulative grid and mixer subspans. New spans do not create individual
  substage traces. No clocks/counters/allocations/communication when disabled
  beyond checking the cached optional token at each instrumented boundary.
- `scripts/lane_glue_parse.py` is explicitly authorized by the revised Phase 0
  request beyond the original owned-file list. It preserves rank/run identity,
  validates completed reports/counters, and emits JSON/Markdown tables. Its
  `--correctness-only` mode suppresses Mac timing values.
- [profile-request.md](profile-request.md) contains the five independent dev08
  jobs, three profile replicates per job, controls, paths, counter semantics,
  coverage/identity gates and interpretation limits.

## Correctness evidence

[mac-summary.json](mac-summary.json), [mac-equivalence.json](mac-equivalence.json),
[mac-equivalence-wide.json](mac-equivalence-wide.json), and
[mac-equivalence-synthetic.json](mac-equivalence-synthetic.json) retain real
commands, exit codes, binary/input/density hashes and counters, without timing
values. `verify.py` compares float bit patterns and density bytes; it compares
baseline, disabled, count-only, clock-enabled and phase-gate-disabled modes.

| input | bands | baseline energy (Ha) | outcome in every mode |
| --- | ---: | ---: | --- |
| Si8, displaced cubic | 20 | -31.32177656320248 | converged, verified |
| Si skew, displaced | 8 | -7.2506046203814964 | converged, verified |
| Si32, displaced cubic | 68 | -126.70343204998849 | converged, verified |
| synthetic coupled s projectors | 6 | 0.6354522013776593 | converged, verified |
| Si elongated, displaced | 12 | -14.48485038813729 | existing 100-iteration exhaustion, exit 3 |
| Si slab, displaced | 12 | -15.222453803896748 | existing 100-iteration exhaustion, exit 3 |

30 executions, 24 comparisons against the pristine baseline: maximum energy
delta **0 Ha**, density bytes identical, scientific history/eigenvalues and
solver work identical. Occupations were also compared from all stored files.
No tolerances changed. The two exhausted fixtures are noninterference checks,
not converged physical validation. The first strict-success-only wide harness
stopped at the baseline elongated exit 3; it was rerun with explicit
`--allow-nonconverged`, requiring the same exit/convergence/science in all arms.

Distributed tests cover cubic/skew/slab grids with PZ and PW LDA, owned and
replicated fields, ranks 1/3/8, exact energy/potential bits and work counts.
The plane-mixer trajectory test covers ranks 1..7, seeded histories, sliding
windows and reset, and checks exact trajectory bits across off/counts/timed.
Simulated ranks are correctness evidence only, not real MPI or timing evidence.

## Commands and actual outcomes

Scratch and all builds are under this worktree's `target/` on the 990 volume.
Use `PATH=/opt/homebrew/opt/rustup/bin:$PATH` where cargo is not on PATH.

| command | exit | result |
| --- | ---: | --- |
| `TMPDIR="$PWD/target/scratch" cargo build --release -p dft-cli` before edits | 0 | pristine binary copied to `target/scratch/dft-cli-baseline` |
| same command after instrumentation | 0 | final numerical instrumentation binary; subsequent Rust edits only tests |
| `python3 -I docs/performance/lanes/GLUE/verify.py --binary target/release/dft-cli --baseline target/scratch/dft-cli-baseline --work target/scratch/equivalence --json-out docs/performance/lanes/GLUE/mac-equivalence.json target/scratch/si8.toml target/scratch/si-skew.toml` | 0 | 10 executions, all verified and bit-identical |
| same verifier with `--allow-nonconverged`, `target/scratch/equivalence-wide`, `mac-equivalence-wide.json`, inputs `si-elongated.toml si-slab.toml si32.toml` | 0 | 15 executions; two unchanged baseline exhaustions |
| same verifier with `--allow-nonconverged`, `target/scratch/equivalence-synthetic`, `mac-equivalence-synthetic.json`, input `syn-s2.toml` | 0 | 5 verified synthetic executions |
| `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0 scripts/check.sh fast dft-numerics dft-solvers dft-engine dft-cli` | 1 | build/fmt/clippy PASS; 627 tests passed, 3 known cache tests failed, 161 ignored, 68 binaries |
| `python3 -I docs/performance/lanes/GLUE/test_parser.py` | 0 | 3 parser contract tests (invalid subcases included) |
| `python3 -I scripts/lane_glue_parse.py target/scratch/equivalence/si8/off/results.json target/scratch/equivalence/si8/counts/results.json target/scratch/equivalence/si8/timed/results.json target/scratch/equivalence-wide/si32/counts/results.json target/scratch/equivalence-synthetic/syn-s2/counts/results.json --correctness-only --json-out target/scratch/parser-smoke.json --markdown-out target/scratch/parser-smoke.md` | 0 | 5 Mac-generated files parsed, timing values suppressed |
| `TMPDIR="$PWD/target/scratch" python3 -I scripts/tests/test_readme_figures.py` | 1 | 8 passed, existing `test_structural_scope_and_sources` phrase assertion failed |
| `git diff --check` | 0 | no whitespace errors |

Final gate log: `target/scratch/check-final.log`. Known cache failures:
`pseudo::cached_fetch_checks_hash_and_provenance`,
`pseudo::cache_resolution_and_no_network_side_effects`, and
`pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`.
They reject the in-repo scratch directory on this Mac; no gate exemption or
production fix was introduced. The load-only `ps` probe is sandbox-denied;
the build/test/fmt/clippy gates themselves still execute.

Earlier attempts: first check used `TMPDIR` instead of the script's
`TMPDIR_990` and could not write its default scratch path (exit 1, no gates
verified); corrected as above. An initial instrumentation compile caught a
gather-hook insertion at the wrong call site (exit 101; fixed before any
candidate execution). New grid-counter expectations initially omitted the
second owned-field call and the stop-estimator gather; the first gate and a
focused test failed these assertions. The expectations now reflect the
fixture's actual five gathers per rank; final gate passes the new tests.
A `si-skew` shorthand input-generation attempt failed because the generator
appends `-gamma`; the full `benchmarks/cases/si-skew.toml` path succeeded.

No README or scientific reference data changed: this is diagnostic
infrastructure with no new physics/performance claim; those files are outside
this lane's ownership. The existing README failure asks for the absent phrase
`standard profile is unsupported`, unrelated to this patch. Dashboard/registry
publication remains with the master; the local registry has no LANE-GLUE run
and no `.delegate/LANE-GLUE/report.json` exists in this worktree.

## Remaining limits

Real MPI, MKL and performance measurements run on dev08 through the master,
not here. No QE source read, no QE run, no remote action, no new dependency,
no UPF file committed. Newly detailed grid spans apply to M2; one-core grid
work retains existing envelopes. Rank totals are inclusive and not additive;
new per-call means are not per-map substage traces. All defaults are unchanged.

The required schema's status enum has only DONE/BLOCKED/FAILED; the transport
report uses DONE for completed Phase 0 and carries **PROFILE_REQUEST** as its
phase outcome. No schema edit is made. Phase 2 is not started.
