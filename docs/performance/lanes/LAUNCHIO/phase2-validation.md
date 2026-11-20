# LANE-LAUNCHIO Phase 2 local validation

Status: **DEV08_AB_REQUEST**, default off, uncommitted. No dev08 run or Mac
performance conclusion. Relayed profile diagnosis is separate from local
correctness evidence. Base `0413485` (the committed Phase-0 implementation).

## Implementation

Added `BAANDI_BINARY_HASH_CACHE`, default off, registered as a performance switch;
updated both switch-count tests (126 -> 127). On uses a bounded, metadata-keyed
Unix sidecar with corruption checks and atomic best-effort publication; on miss
it still computes and records the ordinary content SHA-256. Off preserves the
old fresh process identity. No density writer, engine calculation, output
suppression, MPI collective or numerical threshold changed. Details and explicit
cache trust limits: [dev08-request.md](dev08-request.md).

## Verified results

[phase2-equivalence.json](phase2-equivalence.json) records hashes and all **32
executions**: Gamma, skew, 3%-anisotropic primitive cells and two weighted k points;
base executable, candidate default unset/off/cold/warm/malformed/input-on and
publication-unavailable modes. Pure Rust, one thread, external Si UPF, 30 Ry.

- Every run converges and passes fixed-point verification.
- **0 Ha** energy difference; exact f64 bits retained throughout scientific JSON.
- Direct byte comparison of `density.bin` and `density.json` passes in all arms.
- Entire results.json is compared after normalizing only volatile/legitimate
  fields: timestamp, input SHA (input-table switch test), run ID, agreement hash,
  independently checked binary path/digest, this switch's manifest/override entry,
  explicit cache diagnostics, all existing `timings` and `_s` timing fields and
  the derived `coarse_wall_fine_step_equivalents`. Remaining JSON is canonically
  serialized and compared byte-for-byte, with f64 values encoded as bit patterns.
  All scientific histories, convergence/certificate work and structural metadata
  remain in this comparison. The code lists exclusions explicitly.
- Every reported binary digest equals an independent Python hashlib scan, even
  when publication fails. **Eight warm/input hits** read zero executable payload
  bytes for hashing. Cold/malformed/unavailable arms each hash the full binary.
- Default/off create no sidecar and add no `run.binary.hash_cache` object.
- Five Rust cache tests pass: cold/hit and same-length overwrite with restored
  mtime; inode/size/mode invalidation; malformed/oversized/corrupted-digest and
  writable records; symlink/directory preservation and missing executable error;
  eight concurrent writers producing one valid cache without temporary leaks.
- Supplemental instrumented cold/hit CLI outputs both pass the Phase-0 parser;
  their provenance independently matches the executable and density matches
  the off result. Those two runs used the preliminary candidate before adding
  test-only source lines; the final 32-run matrix uses the final release digest.
- Python paired-statistics checks reproduce a known mean/SD/SE and leave zero SE
  unresolved. Python sources pass AST syntax checks. The remote runner rejects
  Mac invocation with exit 2 **before any timing run**.

Final release executable length is **15,427,152 bytes on this Mac** (not a dev08
bandwidth input). Binary/UPF SHA-256 values are recorded in the JSON evidence.
The final source/binary digest was independently checked against the evidence.

## Commands and exits

All scratch stays under `target/scratch/lane-launchio-p2/` on the 990 volume.
PATH includes `/opt/homebrew/opt/rustup/bin`.

| Command | Exit | Result |
| --- | ---: | --- |
| `TMPDIR=$PWD/target/scratch cargo build --offline --release -p dft-cli` | 0 | Base, candidate and final builds |
| `CHECK_SCCACHE=0 TMPDIR_990=$PWD/target/scratch CHECK_NICE=0 scripts/check.sh fast dft-cli dft-numerics` | 1 | build/fmt/clippy PASS; 222 passed, 3 known cache failures, 62 ignored; 19 binaries |
| `TMPDIR=$PWD/target/scratch cargo test --offline -p dft-cli --lib binary_hash_cache::tests -- --nocapture` | 0 | 5/5 cache tests pass |
| `python3 -I docs/performance/lanes/LAUNCHIO/verify_phase2.py --baseline target/scratch/lane-launchio-p2/baseline-dft-cli --candidate target/release/dft-cli --upf ../pseudo/Si.pz-vbc.UPF --output target/scratch/lane-launchio-p2/identity-byte-final` | 0 | 32 executions, exact file/scientific identity, independent provenance |
| `python3 -I scripts/lane_launchio_parse.py target/scratch/lane-launchio-p2/profile-cold target/scratch/lane-launchio-p2/profile-warm --output target/scratch/lane-launchio-p2/parser-smoke.json` | 0 | Both instrumented cache paths accepted |
| `python3 -I scripts/tests/test_readme_figures.py ReadmeIsCurrent.test_readme_validation_block_matches_the_evidence ReadmeIsCurrent.test_figures_match_the_evidence_and_all_runs_pass` | 0 | Both generated-block/figure freshness checks pass |
| `git diff --check` | 0 | No whitespace errors |

The aggregate fast gate is **not PASS**: its three existing pseudo-cache tests
reject repository-local TMPDIR, as pre-announced in the ticket. Test/policy
semantics were not weakened and scratch was not relocated. Sandbox denies `ps`
load telemetry, but build/test/fmt/clippy still execute; sccache is disabled.
The initial identity harness correctly found two additional legitimate metadata
fields needing explicit exclusions (derived wall equivalence and the switch's
execution override). Those harness attempts exited 1; no scientific fields or
production behavior changed to resolve them. Their output remains in `identity/`
and `identity-2/`. Final raw files and logs are in `identity-byte-final/`.

## Limits and next owner

Real MPI/MKL and its shared-filesystem behavior are not locally validated. The
Rust eight-writer test is a concurrency/atomic-publication test, not an MPI timing
substitute. Cache coherence/trust assumptions and cold-run cost remain explicit.
No large-band/slab/multi-species performance generality or speedup is claimed;
this implementation never examines physical-system data.

The master runs the four registered A/A + A/B jobs, returns raw profiles/counters
and decides adoption. Until then the switch stays off. The source registry lacks
this lane's live run entry; session/dashboard changes belong to the master and
are outside ticket ownership. README scientific evidence is unchanged and its
generated block/figures remain current. No commit is created.
