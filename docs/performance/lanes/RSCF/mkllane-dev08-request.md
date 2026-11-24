# P3-MKLLANE: diagnose the cost around MKL lane kernels

Status: code inspection, descriptor-budget hardening and counters implemented;
**the cause of the measured slowdown and its removal are not yet established**.
No speed claim, default promotion or QE source access. Changes are uncommitted.
Base HEAD: `b69571dbf512a323c17b03732f168c33e2683ae7`.
The master owns dev08 through `baandi-dev08`, the ledger/optimization log,
session/dashboard records and generated README integration. Those shared files
are outside this ticket's OWNED_FILES.

## Evidence and hypothesis order

The ticket reports the following si128 MPI8 stage seconds, rustfft -> MKL:

| Stage | Before (s) | After (s) | Difference (s) |
| --- | ---: | ---: | ---: |
| slab_gather | 0.520 | 0.625 | +0.105 |
| slab_scatter | 0.690 | 0.836 | +0.146 |
| middle_unpack | 0.613 | 0.685 | +0.072 |
| middle_pack | 0.760 | 0.844 | +0.084 |
| multiply | 0.455 | 0.547 | +0.092 |

The differences total 0.499 s; this is a transcription of supplied evidence,
not a new measurement or a critical-path sum of independent rank maxima.
The ledger's Work-order screen records all ten calibrated lane shapes choosing
MKL, 16–40% faster kernels, about 0.9 s kernel saving and wall 33.63 -> 33.46 s.
Raw evidence remains on dev08 at `/data/sig-dft/wo/RESULT.md`; it was not available
to this Mac session. The user reports the same symptom in P3e. The historical
M2_FFT rows already establish that faster FFT calls alone did not justify a
default change. The new premise is explicitly diagnosing the surrounding work.

| Hypothesis | Code evidence | What remains to measure |
| --- | --- | --- |
| (a) Too many MKL threads | `dft-mkl/build.rs` links `mkl_sequential`, not an OpenMP/TBB layer or `mkl_rt`; `src/keep.rs` retains that library under as-needed linking. Previously no descriptor-specific limit was set. Now `Plan1dBatch::new` sets DFTI_THREAD_LIMIT=1 before commit, queries it after commit and rejects any other result. | Resolve the libraries actually loaded on dev08; each profiled call queries the committed limit. An environment variable alone is not evidence of threads used. |
| (b) Extra lane copy/layout | `MklLanes::run` passes the original mutable slice to `Plan1dBatch::transform_in_place`; that plan owns no data buffer. INPLACE, contiguous strides, input/output distance=len, and count transforms. Rust and MKL use the same sticks/middle/slab arrays and pack/scatter/multiply loops. | Explicit adapter copy bytes must be zero. Native MKL scratch/cache traffic is opaque and is not covered by this counter. |
| (c) Commit/allocation in neighbours | Descriptor create/set/commit occur only in construction. Shared descriptors are reused behind a mutex; its lock and dispatch are inside the FFT stage. Batch reserve happens before the per-plane neighbours and is identical across backends. No Rust allocation/copy/commit occurs in the lane transform body. | Native allocator/DFTI stacks and cache misses; first-call lazy initialization must not be mistaken for construction work. |
| (d) FP mode or AVX effects | No existing lane code sets FTZ/DAZ or selects a different ISA per call. Added MXCSR observations do not change FP state. A uniform 12–21% slowdown of these neighbours is compatible with several causes and proves none. | Per-call FP control changes/sticky denormal evidence, delivered cycles/ref-cycles, cache events and sampled ISA/transition activity. |

One real layout exception: `FFT_FUSED_ROWS=on` uses rows only for rustfft;
vendor-capable dispatch falls back to slabs. Force it **off in both arms**.
It is already default off, so this exception alone cannot explain the supplied
results without checking their manifests. No row-fusion change is included.

The descriptor cap closes a missing local budget contract; because the existing
link is sequential, it is **not proof that oversubscription caused the loss**.
No process-global MKL thread/ISA setting, FP mode, arithmetic, layout or default
backend is changed. The cap applies whenever existing M2_FFT=mkl,
FFT_KERNEL_BY_LENGTH or pruned MKL pair dispatch actually constructs MKL lanes.

Intel documents the committed descriptor value as a thread upper bound and
reports 1 in sequential mode: [DFTI_THREAD_LIMIT](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-2/dfti-thread-limit.html).
The budget is expressed as MKL_LONG, consistent with the DFTI C API. Runtime
ABI checks below remain mandatory on Linux x86-64.

## Counters and their limits

Enable the existing `DAVIDSON_STAGES`, or both `FFT_PROFILE` and `PHASE_PROFILE`.
Each MKL lane descriptor keeps bounded local scalar counters, with no clocks,
allocation, locks, formatting or I/O in its diagnostic update. The probe adds
one DftiGetValue query and two MXCSR reads per call; these costs are inside the
FFT stage. Unprofiled calls skip them. Calibration temporarily disables these
probes and resets their counts afterwards, so queries do not bias the selector.
Probe/cache disturbance is asymmetric (only MKL queries DFTI); this instrumented
pair is diagnostic, not a precise estimate of uninstrumented wall savings.

At **workspace teardown**, numerics emits one `P3_MKLLANE` stderr line per
descriptor with PID, run ID, length, batch and `LaneDiagnostics` fields:

- `calls: [forward, backward]`, `points`, `failed_calls`.
- `explicit_copy_bytes`: bytes copied at the owned adapter boundary; zero for
  this direct in-place path. Excludes FFT loads/stores and hidden native copies.
- `thread_limit_min/max`: min/max of DftiGetValue over every profiled call;
  both must be 1 for a nonempty record. Successful calls with that cap use the
  calling thread only. These are descriptor observations, not OS-thread counts
  or the unrelated global `MKL_Get_Max_Threads` setting.
- `mxcsr_samples`, `mxcsr_control_changes`, `ftz_off_calls`, `daz_off_calls`,
  `new_denormal_flag_calls`, `mxcsr_before_or`, `mxcsr_after_or`.
  Control changes exclude sticky exception bits 0..5. The denormal counter
  sees only newly set flags; it cannot count operands or exclude denormals
  when a flag was already set. No flag is cleared. Non-x86 reports no samples.

Records include start/H/density work on the descriptor, excluding calibration;
do not equate their total with H-only FFT_PROFILE counts. Multiple descriptors
can have the same length/batch; sum their calls/points/bytes, take min/max of
limits only for nonempty records. Empty descriptors have no calls and zero
limits. Abnormal termination can lose teardown records; missing records are
unknown, not zero. `Plan1dBatch::diagnostics()` also exposes the typed snapshot.
No native allocation count or measured CPU frequency is invented by this API.

## Requested dev08 job: exactly one profiled off/on pair

Only **si128, MPI8, one execution per arm**, in one exclusive socket-0 window
on physical cores 0–7. No repeated timing campaign, other cell or new QE run.
Use the archived work-order input/UPF/start and identical verification policy;
retain the resolved command, all input hashes and run manifest. Do not rebuild
or run another benchmark concurrently. One arm per fast job, <=600 s; preserve
any failure instead of automatically repeating it.

Build the reviewed patch with the same compiler/native flags as the prior run:

```sh
CARGO_TARGET_DIR=target/p3-mkllane cargo build --release -p dft-cli --features mpi,mkl
cargo test --release -p dft-mkl --features mkl --lib lane_diagnostics -- --nocapture
cargo test --release -p dft-numerics --features mkl --lib mkl_lane_profile -- --nocapture
cargo test --release -p dft-numerics --features mkl --lib mkl_lanes_match_rust_lanes
cargo test --release -p dft-solvers --features dft-numerics/mkl --lib vendor_lanes_match_rust_lanes
```

The unit tests are ABI/correctness checks, not extra si128 runs. Capture
`ldd`/resolved library paths and hashes plus the loaded MKL version. Optionally
set `MKL_VERBOSE=1` for the adapter unit test only to capture its FFT `tLim:1`
and sequential runtime header. Keep verbose **off** for both profiled SCFs:
per-call printing would manufacture a neighbour slowdown. Intel distinguishes
FFT `tLim` from the global `NThr` field in its
[verbose documentation](https://www.intel.com/content/www/us/en/developer/articles/technical/verbose-mode-supported-in-onemkl-112.html).

| Arm | M2_FFT | FFT_KERNEL_BY_LENGTH | FFT_FUSED_ROWS |
| --- | --- | --- | --- |
| off | rust | off | off |
| on | mkl | off | off |

Use explicit backend selection to remove run-local calibration as a confound;
this exercises the same Plan1dBatch used by the ten MKL-selected work-order
shapes. Keep all remaining switches equal to the original manifest. Skeleton
(master resolves actual paths and a separate output destination per arm):

```sh
BAANDI_M2_FFT=rust BAANDI_FFT_KERNEL_BY_LENGTH=off BAANDI_FFT_FUSED_ROWS=off \
BAANDI_DAVIDSON_STAGES=on BAANDI_FFT_PROFILE=on BAANDI_PHASE_PROFILE=on \
BAANDI_STAGE_TIMERS=on MKL_VERBOSE=0 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 8 --map-by core --bind-to core --report-bindings \
  <reviewed-binary> run <frozen-si128-input>
```

Change only M2_FFT to `mkl` and the output destination for on. Capture stderr
separately per rank to avoid interleaving counter records. No environment-only
assertion of thread count: return each nonempty descriptor's queried min/max
and calls, plus the resolved sequential-library evidence.

Collect available counters **during these same two runs**, through the existing
host profiling wrapper, not an additional run: per-rank task-clock, cycles,
ref-cycles, instructions, cache misses, context switches, CPU migrations and
sampled native stacks (DFTI, allocator/memcpy, mutex, neighbours). If permissions
or PMU events are unavailable, return the specific failure; do not install tools
or substitute guessed values. Preserve the exact perf command/event list and
multiplexing coverage. Select model-supported FP-assist/AVX events from the
host's `perf list`; do not hard-code raw event numbers for a different CPU.

Return each rank's stage seconds/calls/points/logical bytes, H columns and maps,
FFT/nonlocal/collective controls, full wall including setup/output, results JSON,
density hash, energy in Ha, verification and exit status. Inspect same-rank
neighbours and normalize by their actual work; never add independent rank
maxima or compare stages with different call counts. Preserve differences in
trajectory/counts caused by backend rounding. Accuracy classes remain strict
1e-10 Ha / standard 5e-9 Ha; do not relax them for this experiment.

## Decision after the pair

1. Unexpected thread limit or non-sequential loaded layer: report the runtime
   mismatch; fail-closed cap prevents oversubscription in the new lane path.
2. Limits 1 and zero explicit bytes: rule out adapter oversubscription/copying,
   not internal MKL scratch traffic. Inspect native allocation/cache samples.
3. FP control changes: isolate the offending native call before adding a
   restoration policy. Zero sticky-flag changes alone do not rule out denormals.
4. Lower delivered clock plus wide-vector samples makes AVX frequency a
   candidate; a rank-wide cycles/ref-cycles ratio alone is not causal proof.
   Any later AVX2 cap must be a separately requested diagnostic: Intel's
   [ISA control](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-2/mkl-enable-instructions.html)
   is process-wide and can change GEMM/3D FFT too. Do not silently apply it to
   this pair or claim it is a lane-only fix.

This pair screens the residual penalty and the hardened path. It cannot by
itself establish a pre/post-patch gain or a general-purpose default. If the
same loss remains, the task is still incomplete; retain existing switches.

## Local validation

Mac outcomes (2026-11-25), logs in ignored `target/mkllane-*.log`:

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace and test targets type-check |
| `cargo check -p dft-solvers -p dft-numerics -p dft-mkl --tests --features dft-numerics/mkl` | 0 | Feature targets type-check; MKLROOT unset |
| `cargo test -p dft-solvers` | 0 | 186 passed, 25 ignored |
| `cargo test -p dft-numerics` | 0 | 91 passed, 22 ignored |
| `cargo test -p dft-mkl` | 0 | 3 passed, no ignored tests; no native MKL execution |
| `cargo fmt --check` | 0 | Clean |
| `cargo clippy -p dft-solvers -p dft-numerics -p dft-mkl --all-targets -- -D warnings` | 0 | Clean |
| Same clippy command with `--features dft-numerics/mkl` | 0 | Clean after correction described below |
| `CHECK_SCCACHE=0 BAANDI_TEST_TMP="$PWD/target/tmp" scripts/check.sh fast dft-solvers dft-numerics dft-mkl` | 0 | All four gates PASS; 27 binaries, 280 passed, 47 ignored |
| `git diff --check` | 0 | Clean |

Direct tests also used worktree-local `BAANDI_TEST_TMP`/`TMPDIR`. The gate's
optional `ps` load inventory was sandbox-denied; build/test/fmt/clippy ran and
passed (gate wall 65 s; not a performance result). Initial feature clippy exited
101 because inline diagnostic storage enlarged `MklPairPlan`; moving optional
counter storage to a construction-only Box fixed it without lint suppression
or a transform-time allocation. Final feature check and clippy both passed.

The new pure adapter tests check counter accumulation and FP control-mask
classification using synthetic register words (120 complex points across two
recorded calls, zero explicit copy bytes); they do not execute native FFTs.
Native tests added for dev08 cover length/batch (1,1), (7,3), (12,5), (25,9),
(36,17), both directions, offset storage, unchanged pointers and output bits,
limit=1, invalid slices, disabled probes and calibration exclusion.
Those tests were **type-checked, not executed on Mac**. The actual x86 MXCSR
instruction, MKL runtime behavior, native allocations, clock/cache effects,
real-system energy/density equivalence and penalty removal all await dev08.
