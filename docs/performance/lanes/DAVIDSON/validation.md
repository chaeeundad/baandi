# LANE-DAVIDSON Phase 0 validation

Outcome: PROFILE_REQUEST. No optimization, numerical policy change, or speed
claim. The user's revised phased instruction is authoritative over the older
charter. The only performance measurements requested are on dev08.

## Code/ledger reading retained

- Standard real Gamma uses the lean generalized small problem; strict uses
  the classic residual path. MPI's current workspace lean path shares factor
  extension, reduced transformation and Ritz coefficients with serial lean.
- Overlap/H borders, normalization and pending residual norms are batched by
  the distributed lean solver. Reduction site labels are available in the
  existing Counted communicator; the new trace measures that same call without
  replacing it or changing the payload.
- Reduced EVD and B5 are distinguished. No range-selected eigensolver,
  tolerance schedule, capacity, restart-keep or cluster-policy change is made.
- Read: AGENT-BRIEF, relevant attempts-ledger/backlog entries, algebra-gap
  analysis, diagonalization plan, fast-loop protocol, ADRs 0003/0004/0041/0043.
  Historical notes do not substitute for a current dev08 profile.

## Scope of the instrumentation wiring

The revised task explicitly requests per-map/per-rank results.json output and
an executable parser. Minimal supporting edits therefore include the engine's
existing profile solve boundary (map context), the CLI report adapter (rank
collection and JSON), distributed workspace helpers, and an engine test, in
addition to the originally named solver/switch files. No dependencies, schema
source files, frozen physics contracts or numerical interfaces are changed.

The trace has 65536 records of 13 f64 fields per rank (6,815,744 bytes of reserved
record payload when enabled). It allocates at the first map/solve boundary;
completed calls append only while capacity remains. Overflow increments a
counter without allocating, and the parser rejects the result. Default OFF
performs no timer reads, trace allocation, counter updates or extra collective;
a cheap flag check remains at instrumented sites. Collection/JSON is after SCF.

## Scientific checks

`correctness.json` contains final measured scientific values and hashes. Raw
CLI files are under `target/scratch/profile-final/`; the pristine binary is
`target/scratch/pristine/dft-cli` built before any source changes from
`db25e70176bd5d111dab6fe1dc2a110766fcb005`. Compare pristine, final OFF, final ON
using `verify_identity.py`: float bits of energy terms, eigenvalues,
occupations, every scientific history/certificate/verification field and exact
density bytes. Only clock fields (`timings`, `wall_s`) are excluded. The new
manifest entry and diagnostic output are expected metadata additions.

`crates/dft-engine/tests/davidson_profile.rs` compares ON/OFF on a skewed
synthetic Gaussian fixture with nonlocal projectors at three simulated ranks,
standard and strict. It checks energy and density bits, history length, each
rank's nonempty trace, map/solve labels and reduction records. Recorder unit
tests cover no allocation while disabled, bounded overflow, context, and exact
collection from unequal/empty rank traces. These are transport and correctness
tests, not MPI performance tests.

Parser tests check aggregation, FLOP/payload models and rejection of nine
malformed/incomplete conditions. A real Mac-generated results.json is parsed
as well. Its clock values are not analyzed or retained as performance evidence.

## Commands and failures

Final command exits and gate summary are appended below when complete.

Exploratory failures retained honestly:

- First fast gate used the script's default sccache; sandbox denied its server
  operation. Retried with `CHECK_SCCACHE=0` (no change to the gate script).
- The gate's optional `ps` load-report subprocess is denied in this sandbox;
  the gate continues. This is not a test success or a performance measurement.
- First CLI build of the instrumentation referenced dft-solvers directly, but
  CLI does not have that normal dependency. Corrected by re-exporting the
  diagnostic module through the existing engine profile adapter; no manifest
  or dependency changes.
- An intermediate gate caught the new switch classified as Performance rather
  than Profile. The classification is corrected and explicitly tested. That
  run was during edits; only the final unchanged-source gate is authoritative.
- A first scientific-comparison assertion excluded only top-level history
  timings and therefore found nested certificate clocks. The delivered
  verifier excludes clock fields recursively and compares all remaining f64
  values by their bits; no numerical tolerance was introduced.
- Three existing CLI pseudopotential-cache tests reject repository-local
  TMPDIR. Scratch remains on the required 990 volume inside this worktree;
  these tests are not changed or skipped to hide the failure.

## Limits

No dev08, real MPI, MKL execution, performance conclusion or adoption decision.
No new scientific README claim: this phase changes observability only, not a
physical result or reference dataset. Dashboard/session registry updates and
the remote job belong to the master session. No commit is created.

## Final unchanged-source results

- `cargo build --release -p dft-cli`: exit **0**. Binary SHA-256
  `8fd32f7c401595030a86c2062dc409e6ae946997ba07056bf8e7d9477db6d43c`.
- `CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast
  dft-solvers dft-numerics dft-engine dft-cli`: exit **1**. Build, fmt and
  clippy **PASS**. 68 test binaries: **613 passed, 3 failed, 160 ignored**.
  Failures are exactly the known repository-TMPDIR restrictions:
  `pseudo::cache_resolution_and_no_network_side_effects`,
  `pseudo::cached_fetch_checks_hash_and_provenance`, and
  `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`.
  All new recorder, simulated-rank and switch tests pass. Log:
  `target/scratch/check-profile-final.log`. The intermediate clippy
  `chunks_exact_to_as_chunks` diagnostic was fixed using `as_chunks::<13>()`;
  no lint allowance or relaxed test was added.
- `cargo test --release -p dft-engine --test davidson_profile`: exit **0**
  during development; the final gate also executes this test successfully.
- `TMPDIR="$PWD/target/scratch" python3 -I
  scripts/tests/test_lane_davidson_parse.py`: exit **0**, two tests, including
  nine rejection cases.
- Six final CLI runs (Si8/Si32 x pristine/OFF/ON): all exit **0**. Two calls
  each to `verify_identity.py` and `lane_davidson_parse.py`: all exit **0**.
  Si8 energy -31.321776563202498 Ha; Si32 -126.70343204998852 Ha; every paired
  energy delta exactly **0 Ha**, density bytes and full scientific histories
  including verification identical. ON traces: **773 / 821 calls**, **18 maps
  each**, **0 dropped records**, all map stage/kind labels known and solve
  tolerance values positive. These are correctness facts, not performance
  measurements. Full hashes and result paths: `correctness.json`.
- `bash -n docs/performance/lanes/DAVIDSON/profile-job.sh`: exit **0**.
- `git diff --check`: exit **0**.

The request's status is recorded verbatim in `phase-status.json` as
`PROFILE_REQUEST`. The existing report schema permits only DONE/BLOCKED/FAILED,
so `report.json` uses DONE for the completed Phase 0 deliverable and states the
phase outcome separately in its design choices. Phase 2 is not started.
