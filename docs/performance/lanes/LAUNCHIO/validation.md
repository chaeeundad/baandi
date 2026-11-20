# LANE-LAUNCHIO Phase 0 validation

Phase status: **PROFILE_REQUEST** (2026-11-20). Implementation and Mac correctness
checks are complete; remote profiling is not executed. No optimization is
proposed and no default changes. All changes are uncommitted.

Base commit: `61822bca3b470df1212ee47b47c7eb4e786e22e3`. A release executable was
built from unchanged HEAD `run.rs`, confirmed with a second build, and saved as
`target/scratch/lane-launchio/baseline-dft-cli` before restoring the instrumented
file. Both use the pinned pure-Rust toolchain, not MKL or real MPI.
Final binary and UPF SHA-256 values are in [equivalence.json](equivalence.json).
The same external UPF was used by both builds; none was added to the repository.

## Correctness results

`verify_cli.py` performs no external timing. Four inputs at 30 Ry, one thread:
Si2 Gamma, a skewed primitive cell, a 3% anisotropic cell, and two explicit
weighted k points. Each runs the base binary, then the final binary with default
unset, explicit profile off, environment profile on, and input-table profile on.
**20/20 executions** converge and pass fixed-point verification. For every case:

- Total energy, every energy term, eigenvalues and occupations match as f64 bits.
- Density binary and density metadata bytes match (SHA-256 recorded).
- Scientific iteration history and verification reports match, excluding only
  existing `timings` and `_s` wall-time fields.
- `launchio` is absent with unset/off and present with either enable mechanism.
- Eight instrumented outputs pass the parser, including local bucket/launch
  closure, mandatory fields, actual density hash and payload sizes.

| Case | Energy (Ha) | Main history iterations | Delta (Ha) |
| --- | ---: | ---: | ---: |
| Gamma | -7.295530468166501 | 4 | 0 |
| Skew | -7.290036922765172 | 4 | 0 |
| Anisotropic (3%) | -7.305566147681299 | 4 | 0 |
| Two k points | -7.476578066562224 | 8 | 0 |

These are same-code identity oracles, not QE references or convergence-to-exact
physics claims. No Mac timing is summarized or used as performance evidence.
Large bands, multiple species, slabs and real MPI/MKL remain remote scope;
there is no default-adoption or general speed claim.

## Commands and real exits

Run from the worktree, with `/opt/homebrew/opt/rustup/bin` on PATH. All build and
run scratch stays on the worktree's 990 volume.

| Command | Exit | Result |
| --- | ---: | --- |
| `TMPDIR=$PWD/target/scratch cargo build --offline --release -p dft-cli` | 0 | Base, confirmation, candidate and final builds succeeded |
| `TMPDIR_990=$PWD/target/scratch CHECK_NICE=0 scripts/check.sh fast dft-cli` | 1 | Initial environment failure: sandbox rejects sccache; fmt passed |
| `CHECK_SCCACHE=0 TMPDIR_990=$PWD/target/scratch CHECK_NICE=0 scripts/check.sh fast dft-cli` | 1 | build/fmt/clippy PASS; 132 tests pass, 3 fail, 40 ignored; 13 binaries |
| `python3 -I docs/performance/lanes/LAUNCHIO/verify_cli.py --baseline target/scratch/lane-launchio/baseline-dft-cli --candidate target/release/dft-cli --upf ../pseudo/Si.pz-vbc.UPF --output target/scratch/lane-launchio/equivalence-complete` | 0 | 20 executions, exact scientific identity; eight parser smoke checks |
| `python3 -I docs/performance/lanes/LAUNCHIO/test_parser.py target/scratch/lane-launchio/equivalence-final/gamma-on` | 0 | Five tests: real output, optional resource/external records, malformed/duplicate/truncated/missing/NaN data, corrupted density, synthetic two-rank transport and mixed-run rejection |
| `python3 -I scripts/tests/test_readme_figures.py ReadmeIsCurrent.test_readme_validation_block_matches_the_evidence ReadmeIsCurrent.test_figures_match_the_evidence_and_all_runs_pass` | 0 | Two generated-block/figure freshness checks pass |
| `git diff --check` | 0 | No whitespace errors |

Expected negative checks: `profile_dev08.py` invoked on this Mac exits **2**
before reading input or launching a process; parser `--require-external` on a
Mac correctness run without sidecars exits **2**, rather than inventing external
time or RSS. Python files also pass AST syntax checks. The dev08 launcher is not
executed on Linux here; host enforcement and source syntax are the local checks.

## Failures and limits preserved

- Gate `check-no-cache.log` records the ticket's known repository-local TMPDIR
  failures: two `pseudo` integration cache tests and one CLI binary cache test.
  They reject a pseudopotential cache under a repository. Tests/paths were not
  weakened or moved off the authorized volume. The aggregate gate is **not PASS**.
- `ps` load telemetry is denied by the sandbox; the gate still runs all four
  actual checks. `CHECK_SCCACHE=0` avoids the denied sccache process; clippy passes
  with `-D warnings`. Logs remain under `target/scratch/lane-launchio/`.
- The first identity harness's 20%-elongated cell exhausted convergence in the
  **base binary**, exit 3. Its output/log remains under `equivalence/`; no identity
  or performance claim is made for it. Final coverage explicitly uses a smaller
  anisotropic deformation. The initial unused two-k input syntax was corrected
  to weighted rows before that case was executed. No solver tolerance changed.
- An initial README freshness invocation used a wrong unittest class name and
  exited 1; corrected named tests both pass. No README generation is needed:
  Phase 0 supplies no new physics reference result, and freshness is verified.
- MPI init/finalize in Mac output use the single-rank stub. Only dev08 can supply
  real rank timing, page-fault/RSS data, affinity and launch residual evidence.
- Existing preflight internals, per-map glue call counts and allocator traffic
  are not instrumented; see the explicit schema limits in profile-request.md.
- Session/dashboard artifacts are outside ticket ownership and this checkout's
  session registry has no LANE-LAUNCHIO entry. No nested delegate was started.
  Registry/dashboard updates remain the master session's responsibility.

The top-level report schema permits only DONE/BLOCKED/FAILED; report DONE means
Phase-0 deliverables are done. The requested lifecycle state is separately
recorded verbatim in [phase-status.json](phase-status.json): PROFILE_REQUEST.
