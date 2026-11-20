# LANE-INIT Phase 2 implementation report

Default-off candidate implemented; **performance/adoption pending dev08**.
No Mac profiling or performance conclusion. No commits made. Source base:
`c0f2143769f0876931538e0e9d24c5e0c8202de0`.

## Change and rationale

The returned dev08 profile identifies serial atomic-span MGS as 1.21804 s
of the 2.19395 s si64 one-core atomic start. The new
`BAANDI_ATOMIC_LOCAL_BLOCK=on` reuses the already adopted guarded atomic
CholQR2 with a local Gram instead of a collective; it returns directly
without an unnecessary gather. Its numerical guards, positive factor
diagonals, column order and MGS fallback are unchanged. Only atomic-start
construction is affected, never Davidson correction panels. An active
multi-rank communicator excludes the local option even if atomic splitting
is disabled. Default is off. No material/size-specific policy is introduced.

The user's Phase 2 instruction explicitly targets the implementation
“wherever atomic.span/orthonormalize/hq are produced”; this authorizes the
necessary `crates/dft-engine/src/atomic_start.rs` extension to the original
constructor-only file list. No solver/hot-map implementation was edited.
The two candidates, work floors, inclusive-timer caveats and A52 premise
change are in [diagnosis.md](diagnosis.md). The exact staged remote request
and statistical decision are in [dev08-request.md](dev08-request.md).

INIT_PROFILE adds 28 counters (seven call/ns pairs, dimensions/payload,
MGS work upper bound, local attempts/accepts/fallbacks and block FLOP).
Clocks require STAGE_TIMERS; disabled hooks do not allocate, update metrics,
perform collectives or read a clock. Counts are run-local and cumulative;
MPI values retain local rank semantics. The existing parser handles these
labels without a schema change. The switch registry has 119 entries.

## Mac correctness evidence

[phase2-equivalence.json](phase2-equivalence.json) contains hashes, commands,
exit codes, energies, density hashes and counters for 15 complete SCFs:
Si8, Si32 and skew Si, each baseline / new-off / new-on / strict-off /
strict-on. All are pure-backend, one worker, 30 Ry, seed 601, displacement
1e-3 fractional, real Si.pz-vbc. Strict uses BAANDI_REFERENCE_STRICT and
correctly reports no separate standard verification object. Standard runs
all pass the fixed-point verification. The baseline executable was saved
before these edits. Timers are disabled in the correctness harness.

| check | result |
| --- | ---: |
| maximum standard on/off energy difference | 2.842170943040401e-14 Ha |
| maximum strict on/off energy difference | 1.4210854715202004e-14 Ha |
| maximum standard-on vs strict-off energy difference | 1.0450662557559554e-10 Ha |
| new-off vs saved baseline | energy, density and scientific history bit-identical |
| Si32 atomic MGS calls, off -> on | 1 -> 0 |
| Si32 local accepted / fallback, on | 1 / 0 |
| Si32 block Gram / factor / update calls, on | 2 / 2 / 2 |
| Si32 coarse atomic rows / columns | 4189 / 128 |
| Si32 Gram / update real FLOP, on | 549060608 / 549060608 |

The 1e-11 Ha cross-arm rounding check and strict 1e-10 Ha / standard 5e-9 Ha
accuracy gates pass. This is rounding equivalence, not on/off density bit
identity. No tolerance was relaxed. Synthetic atomic-start SCFs with one
cubic atom (six bands, four atomic columns, hence completion) and two skew
atoms pass with the local switch on and simulated MPI ranks 2/3/5/7. The
algebra test includes local and distributed complex/Gamma matrices, unequal
column scales, rank deficiency, zero columns and injected fill failure;
Q vs MGS <1e-12 and orthogonality <1e-14.

## Commands and outcomes

All scratch paths below are on this worktree's external-volume target/.
The stored evidence includes exact per-SCF commands. Compiler/test log
elapsed fields are operational output, never performance evidence.

- `TMPDIR="$PWD/target/scratch" cargo build --release -p dft-cli`: exit 0.
- `TMPDIR="$PWD/target/scratch" cargo test --release -p dft-engine --lib block_cholesky_qr_matches -- --nocapture`:
  exit 0, one test passed; also covered by the final fast gate.
- `BAANDI_ATOMIC_LOCAL_BLOCK=on BAANDI_INIT_PROFILE=on BAANDI_STAGE_TIMERS=off TMPDIR="$PWD/target/scratch" cargo test --release -p dft-engine --test m2_coarse_certificate block_atomic_start_stays -- --nocapture`:
  exit 0, one test passed (two synthetic systems and simulated ranks 2/3/5/7).
- `python3 -I docs/performance/lanes/INIT/phase2_validate.py --binary target/release/dft-cli --baseline target/scratch/phase2/baseline-cli --upf ../pseudo/Si.pz-vbc.UPF --work target/scratch/phase2/scf-final --json-out docs/performance/lanes/INIT/phase2-equivalence.json`:
  exit 0, 15 full SCFs and scientific checks passed.
- `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0 CHECK_NICE=0 scripts/check.sh fast dft-numerics dft-engine dft-cli`:
  **exit 1**, build/fmt/clippy PASS; 51 test binaries, 471 passed,
  3 failed, 135 ignored. The three existing CLI pseudopotential-cache tests
  reject the required in-repository TMPDIR (“pseudopotential cache must be
  outside repositories”). No source change or temp-disk policy bypass was
  made. The script's `ps` load probe also receives a sandbox permission error;
  the gates still execute. Logs: target/scratch/phase2/check-final.log.
- `python3 -I docs/performance/lanes/INIT/test_parser.py`: exit 0, 3 tests.
  The parser also parsed an on-arm Mac results.json with the new labels.
  Five synthetic checks of paired_stats covered significant improvement,
  floor failure, regression, zero difference and zero-SE improvement: passed.
- `git diff --check`: exit 0.

An initial correctness-harness run returned exit 1 because it compared
verification wall fields and required a standard certificate on strict
reference runs. Corrected the harness to strip diagnostic times and require
standard verification only; reran all SCFs. No numerical implementation or
gate tolerance was changed to address that harness mistake.

## Limits and handoff

Actual MPI/MKL and si64 candidate performance remain unmeasured here. The
raw dev08 profile was relayed, not independently reconstructed; shape
assumptions must be checked against the newly added start-grid counters.
Shifted-k full SCF, broader pseudopotentials, thread independence and the
full generality/adoption gate remain with the staged dev08 request. The
local complex algebra and synthetic completion/fallback cases are covered,
but they do not replace those production gates.

No README/reference/dashboard artifacts were changed outside lane ownership;
this is an opt-in candidate, not an adopted benchmark result. The master
owns the eventual full pre-push gate and dashboard publication. Leave all
changes uncommitted and stop at this handoff.
