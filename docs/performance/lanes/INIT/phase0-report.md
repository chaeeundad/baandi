# LANE-INIT Phase 0 report

Phase status: **PROFILE_REQUEST**. Implementation scope is instrumentation
only. No optimization remains in this diff; the earlier prepared-artifact
reuse experiment was removed when the authoritative revised instruction
arrived. No Mac timing is used as performance evidence. The old Mac timing
scratch is untracked and is not an acceptance record.

## Delivered

- Default-off `INIT_PROFILE`, registered as a `[profile]` switch: 20 call/time
  pairs and 23 work/shape counters. Existing PHASE_PROFILE supplies atomic,
  Ewald and session-setup context. Instrumentation is initialization-only.
- [Exact dev08 profile jobs and metric definitions](profile-request.md): three
  MPI8 cells and two one-core cells, one case per job, three timed profiles
  plus disabled/count-only correctness controls; no QE comparison.
- [Parser](../../../../scripts/lane_init_parse.py): inclusive local seconds per
  call, work and logical payload floors; preserves rank/file boundaries and
  labels existing aggregate phase envelopes separately.
- [Correctness-only runner](verify.py), [parser tests](test_parser.py) and
  [stored Mac equivalence](mac-equivalence.json).

## Correctness evidence

Base source `db25e70176bd5d111dab6fe1dc2a110766fcb005`. The baseline binary
was built before instrumentation. Its SHA-256 and the instrumented binary's
SHA-256, exact inputs/commands and per-run density/history hashes are in
`mac-equivalence.json`. Si.pz-vbc.UPF remains external (no UPF is committed).

Si8 / Si32 (30 Ry, seed 601, fractional displacement +/-0.001, one worker):
all four modes (baseline, instrumented-off, count-only, timed) converge and
verify; all scientific outputs compared are bit-identical. Eight runs,
maximum energy difference **0 Ha**. Tests strip only explicit `timings`
and verification `wall_s` fields. Final total energies are
-31.321776563202498 Ha and -126.70343204998852 Ha respectively. These are
self-equivalence values, not QE references or material accuracy claims.

The first correctness-runner invocation falsely rejected verification wall
clock fields; the diagnostic comparison isolated only `wall_s` differences
with identical density SHA-256. The runner now excludes that explicit timing
field and the complete eight-run matrix was repeated successfully. No
numerical tolerance or comparison of scientific fields was relaxed.

Physics unit test `init_profile_preserves_owned_and_eager_projectors` passes
448 coefficient-bit comparisons across four cell shapes, Gamma/shifted k,
two radial rules, s/p/d/f and repeated p channels, eager and selected rows,
partitions 1/2/4. It also checks byte/phase counts, no lazy materialization
at construction and zero metrics when disabled. This is not an MPI timing
or communicator test.

The parser has three contract tests and was run on actual off/counts/timed
Mac `results.json` files. It distinguishes missing timers from zero measured
work and rejects absent/malformed lane counters. Mac parser tables stay in
`target/scratch`; they are correctness smoke output, not performance evidence.

## Verification commands

All scratch is under this worktree's `target/scratch` on the 990 volume.

| command | exit | outcome |
| --- | ---: | --- |
| `TMPDIR="$PWD/target/scratch" cargo build --release -p dft-cli` | 0 | instrumented release build |
| `TMPDIR="$PWD/target/scratch" python3 -I docs/performance/lanes/INIT/verify.py --binary target/release/dft-cli --baseline target/scratch/dft-cli-baseline --work target/scratch/identity2 --json-out docs/performance/lanes/INIT/mac-equivalence.json target/scratch/si8.toml target/scratch/si32.toml` | 0 | eight scientific identity controls |
| `TMPDIR="$PWD/target/scratch" cargo test --release -p dft-physics --lib init_profile_preserves -- --nocapture --test-threads=1` | 0 | 448 constructor comparisons |
| `TMPDIR="$PWD/target/scratch" cargo test --release -p dft-numerics --test all switches::init_profile_off_count_only_and_scoped -- --nocapture` | 0 | off/count-only/scoped counter contract |
| `TMPDIR="$PWD/target/scratch" python3 -I docs/performance/lanes/INIT/test_parser.py` | 0 | 3 parser contract tests |
| `python3 -I scripts/lane_init_parse.py target/scratch/identity2/si8/off/results.json target/scratch/identity2/si8/counts/results.json target/scratch/identity2/si8/timed/results.json target/scratch/identity2/si32/timed/results.json --json-out target/scratch/parser-mac.json --markdown-out target/scratch/parser-mac.md --label 'Mac parser correctness only; no performance conclusions'` | 0 | 4 real files parsed |
| `python3 -I scripts/tests/test_readme_figures.py` | 1 | 8 pass, existing stale-phrase scope test fails |
| `git diff --check` | 0 | whitespace check |

Required common gate:

```sh
TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0 \
  scripts/check.sh fast dft-numerics dft-physics dft-engine dft-cli
```

Actual exit **1**: **build PASS, test FAIL, fmt PASS, clippy PASS**.
56 test binaries executed: **561 passed, 3 failed, 153 ignored**.
All three failures are the ticket's known in-repository TMPDIR cache-policy
rejections (`pseudopotential cache must be outside repositories`):

- `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`
- `pseudo::cache_resolution_and_no_network_side_effects`
- `pseudo::cached_fetch_checks_hash_and_provenance`

No new test or lint failure was observed. The gate is **not** reported as
passing; the master must run it on dev08 with its approved external scratch.
Raw gate output: `target/scratch/check.log`.
The load preamble cannot invoke `ps` in this sandbox, and `nice` reports a
setpriority warning, but the script continues and actually executes the test
binaries. Do not mistake these diagnostics for skipped build/test gates.

## Limits and handoff

- Actual MPI/MKL execution and all performance diagnosis remain on dev08.
- Counters are per-process cumulative over construction calls, including
  cutoff-ladder stages. Per-call output is a mean, not a call trace. Existing
  phase envelopes contain max/min over ranks and are never summed as local
  times; session-setup context covers all session setups.
- Radial/angular factor time is grouped to avoid per-G clocks. Sizes and
  radial sample counts support later attribution; logical payload floors are
  not measured bandwidth or hardware instruction counts.
- This phase adds no scientific validation result type; README generation
  and dashboard/session changes are outside the ticket's owned files.
  The existing README check has a stale phrase expectation (`standard profile
  is unsupported`) unrelated to this instrumentation; it is not corrected
  here. The orchestrator owns dashboard/session updates per the lane charter.
- Changes are uncommitted. The master reviews/commits and requests the dev08
  profile. This session stops at PROFILE_REQUEST; Phase 2 waits for the relay.
