# P3-RMM: bounded RMM-Davidson implementation and dev08 request

Base: `cab8245bbcbcdd7b1bbd4546f8fd377cd6a17de8`, uncommitted task diff.
Derived from the own-words note
[qe-rmm-davidson.md](../../../research/qe-study/qe-rmm-davidson.md) ONLY.
No QE source was accessed. `REFSCF_RMM` defaults OFF. No dev08 execution,
new QE reference, speed claim or default promotion is made by this task.

## Premise and implementation

The ticket relays approximately 6.5 s of non-H Davidson work on si128 MPI8
and a 2–3% QE RMM gain over QE Davidson. These are motivation, not locally
reproduced numbers or a new reference bundle. A22/A23 converged the inner
RMM solve with repeated sweeps, certification and fallback. This module
instead accepts one sweep of at most four steps, carries H X, and relies
on the SCF iterations plus scheduled Davidson and independent verification.
Thus the earlier negative result does not test this regime.

The driver resolves `[performance] refscf_rmm="on"` / `BAANDI_REFSCF_RMM=on`
once. Davidson runs at iterations 1, 2 and every multiple of five; a first
map redo retains iteration 1. RMM never retries or falls back to Davidson
within a map, regardless of its freeze flags. Verification maps always use
Davidson. Existing failed-verification reuse remains unchanged: when enabled,
a reused map is the already computed Davidson verification, not a new RMM
call. No stopping, energy, mixing or verification formula changes.

The additional run-scoped controls `refscf_rmm_steps` (default 4, integer
1–16) and `refscf_rmm_block` (default 16, integer 1–64) expose the note's
parameters for future diagnostics, not system-specific tuning. All three
controls participate in the usual cross-rank configuration agreement hash.
The default dev08 arms below leave both numerical parameters at defaults.

Algorithm and accounting: [RMM-Davidson](../../../design/reference-scf.md#rmm-davidson).
The CLI emits `solver` on reference iteration records when the RMM switch
is enabled, and full-diagnostic `reference` stage data include `rmm.solves`.
Each RMM record identifies its SCF iteration, exclusive stage seconds, exact
H columns, per-step active counts, per-stage GEMM FLOPs, newest-only DIIS
fallbacks, step clamp counts, QR fallbacks and minimum pre-normalization
GS norm. Davidson's existing work table contains Davidson solves only;
add the separate RMM GEMM totals for whole-run counts. Legacy fields named
`davidson_iterations` count sweep steps on RMM maps; the new `solver` field
is authoritative. RMM does not populate Davidson-only `inner_clocks`.

## Local correctness evidence

The strict comparison fixtures (same 1e-12 Ry estimator in both arms) pass:
all four shapes converge and verify on simulated ranks 1/3/8, with maximum
|E_M-E_D|=1.425526364e-13 Ha. Diagnostics on/off reproduce all scientific
trajectory bits at each fixed rank count. Synthetic fixtures are not QE references.

**Default-estimator acceptance limitation:** an additional development run
applied the same 1e-10 Ha gate at the default 1e-8 Ry estimator. All four
fixtures converged and verified but exceeded that gate (see table below).
That command exited 101; it is not reported as passing. The unchanged
standard stopping/verification contract does not imply a 1e-10 Ha energy
agreement between distinct approximate trajectories. No algorithm constant,
threshold, stopping rule, verification rule or acceptance bound was tuned
to remove this result. The 1e-10 gate remains on the strict fixtures and
in the dev08 request. The default-setting 1e-10 gate fails on these fixtures;
changing the stopping contract is outside this ticket.

Strict-estimator evidence, 1e-12 Ry; 36 runs (D, M, M repeat with diagnostics
off and two workers, for four shapes and ranks 1/3/8). Every run uses the
current energy-based verification switch ON. All 24 D/M runs passed on
one verification attempt; the repeat matches all scientific trajectory bits.
The earlier density-based verification variant also passed the same strict
comparison. Counts include setup H calls within each solve, redo/retry and
verification maps; they exclude atomic-start H work. No timing conclusion
is drawn from these small synthetic cases.

| Shape | Simulated ranks | SCF iterations D / M | Maps D / M | H columns D / M | GEMM FLOPs rank 0 D / M | E_M-E_D (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cubic | 1 | 9 / 10 | 11 / 12 | 127 / 139 | 337868 / 263202 | -1.425526364e-13 |
| cubic | 3 | 9 / 10 | 11 / 12 | 127 / 139 | 155972 / 115538 | -1.418865025e-13 |
| cubic | 8 | 9 / 10 | 11 / 12 | 127 / 139 | 106364 / 75266 | -1.423305918e-13 |
| elongated | 1 | 8 / 11 | 9 / 12 | 106 / 141 | 324296 / 306748 | +1.998401444e-14 |
| elongated | 3 | 8 / 11 | 9 / 12 | 106 / 141 | 145624 / 130520 | +1.820765760e-14 |
| elongated | 8 | 8 / 11 | 9 / 12 | 106 / 141 | 76904 / 62740 | +1.820765760e-14 |
| nonlocal | 1 | 26 / 35 | 27 / 36 | 411 / 631 | 1470904 / 1337866 | -8.326672685e-14 |
| nonlocal | 3 | 26 / 35 | 27 / 36 | 411 / 630 | 674064 / 605234 | -7.882583475e-14 |
| nonlocal | 8 | 26 / 35 | 27 / 36 | 411 / 631 | 456744 / 405578 | -8.459899448e-14 |
| skew | 1 | 11 / 11 | 13 / 13 | 159 / 155 | 573190 / 412610 | +1.110223025e-14 |
| skew | 3 | 11 / 11 | 13 / 13 | 159 / 155 | 267118 / 185214 | +1.243449788e-14 |
| skew | 8 | 11 / 11 | 13 / 13 | 159 / 155 | 149398 / 97754 | +1.221245327e-14 |

Minimum recorded GS norm over the strict fixtures: 0.9999732015. The
nonlocal case changes 26 to 35 outer iterations and 411 to 630/631 H
columns despite fewer dense FLOPs. At three ranks M freezes one extra
band-step; cross-rank-count work identity is not claimed. Energies remain
within the unchanged 1e-10 Ha gate across rank counts.

Default-estimator 1e-8 Ry **failed acceptance** (one rank, both arms converge
and pass the current energy verification on their first attempt):

| Shape | Simulated ranks | SCF iterations D / M | Maps D / M | H columns D / M | GEMM FLOPs rank 0 D / M | E_M-E_D (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cubic-default | 1 | 6 / 6 | 8 / 8 | 96 / 99 | 258842 / 216592 | -2.871041183e-10 |
| elongated-default | 1 | 6 / 7 | 7 / 8 | 85 / 102 | 261756 / 241638 | +3.494387002e-10 |
| nonlocal-default | 1 | 21 / 28 | 22 / 29 | 356 / 519 | 1325616 / 1129472 | +2.626654450e-10 |
| skew-default | 1 | 7 / 8 | 9 / 10 | 119 / 124 | 455790 / 353028 | -6.115339346e-10 |

Maximum |E_M-E_D| is 6.115339346e-10 Ha, exceeding the ticket's 1e-10 Ha
bound. The four `rmm_default_*_strict_energy_gate` tests retain the same
failing assertions and are explicitly ignored as known acceptance failures,
not runtime exclusions. Reproduce with
`cargo test -p dft-engine --test refscf rmm_default_ -- --ignored --nocapture`
(exit 101, four failed). **Default-estimator acceptance remains blocked.**
No standard tolerance was tightened in the production code to manufacture
agreement, and no failing energy gate was relaxed.

Five new solver tests cover sorted low roots of separated/degenerate dense
operators after 30 sweeps, fresh-H consistency of carried H X, orthonormality,
fixed-rank repeatability at 1/3/8 ranks and diagnostics on/off; canonical
DIIS normalization and newest-entry fallback; polynomial scaling and Ry/Ha
thresholds; the fixed schedule; exact constant-state freezing at T=0 with
empty owners; collective invalid-input/dependent-start errors; guarded
panel fallback and H X transformation. Dense-root tolerance is 1e-10 Ha,
orthonormality 2e-12 and carried-H error 2e-13 Ha. No cross-rank bit identity
or dev08/MKL performance is inferred.

Raw local logs (ignored build artifacts) are under `target/tmp/p3-rmm-*`.
A first unit-test compile failed on an ambiguous floating literal (fixed
with an explicit f64 type); the first clippy run found a diagnostic re-export
after a test module (moved before it). Neither was a numerical failure.
The default-estimator failures above persist after those fixes.

The tests print per-arm SCF iterations, total maps including redo/retry and
verification, H columns and rank-local GEMM FLOPs. Different rank counts
change the local FLOP partition, so compare D/M within one rank count only.

Reproduce the new tests and requested gates:

```sh
cargo test -p dft-solvers --lib ref_rmm -- --nocapture
cargo test -p dft-engine --test refscf rmm_verified -- --nocapture
cargo check --workspace --tests
cargo test -p dft-solvers
cargo test -p dft-engine --test refscf
cargo test -p dft-numerics switches
cargo fmt --check
cargo clippy -p dft-solvers -p dft-engine -p dft-numerics -p dft-cli --all-targets -- -D warnings
```

Final local gate results (actual exits, pure Rust on Mac):

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace/test compilation |
| `cargo test -p dft-solvers` | 0 | 185 passed, 25 ignored; five new RMM tests pass |
| `cargo test -p dft-engine --test refscf -- --nocapture` | 0 | 32 passed, 5 ignored (four are the explicit known failures above) |
| `cargo test -p dft-numerics switches` | 0 | 14 passed |
| `cargo test -p dft-cli --test all run::expressions_preserve_energy_and_record_resolutions` | 0 | Switch manifest count and unchanged legacy scientific output |
| `cargo fmt --check && cargo clippy -p dft-solvers -p dft-engine -p dft-numerics -p dft-cli --all-targets -- -D warnings` | 0 | Formatting and touched-crate lints |
| `scripts/check.sh dft-solvers` | 0 | Build/test/fmt/clippy all PASS |
| `cargo test -p dft-engine --test refscf rmm_default_ -- --ignored --nocapture` | 101 | Four known default-estimator energy-gate failures |

For the wrapper, `CHECK_SCCACHE=0 CHECK_NICE=0` and
`BAANDI_TEST_TMP=$PWD/target/tmp TMPDIR_990=$PWD/target/tmp` keep all writes
inside the task sandbox. Its optional load-report `ps` probe was denied by
the sandbox; the four verification gates completed and the wrapper exited 0.
No workspace-wide test run or remote measurement was performed.

The task owns neither the central attempts ledger/optimization log nor the
sessions/dashboard artifacts, README generator/reference JSON and schema.
This document is the integration record for the master to reconcile those
files. The local registry has no P3-RMM entry and the local delegate report
is absent; no recursive delegate process was started inside this task.

## dev08: exactly ONE profiled full run per BAANDI arm/case

The master sends this request to `baandi-dev08`. Use the frozen P3 manifests
for **si64, si128, elongated, rskew, MPI8**, with one worker, OpenMP thread
and MKL thread per rank on physical socket-0 cores 0–7. Preserve the exact
cell/coordinates, pseudopotential bytes, cutoff, FFT grid, band count,
occupations, XC, deterministic atomic start, threshold and all unrelated
adopted switches. Resolve absolute input paths and record SHA-256; case
names alone do not identify input bytes. No installs or server changes.

| Arm | BAANDI_REFSCF_RMM | steps | GS block |
| --- | --- | ---: | ---: |
| D | off | 4 (inactive) | 16 (inactive) |
| M | on | 4 | 16 |

Run both arms with production verification ON for all four cases, plus
both arms with `scf.verify="off"` for si128 only: **10 BAANDI invocations**.
Each is a separate <=600 s fast job, e.g. `fast-p3-rmm-si128-M-verify-on`.
Keep all failures/timeouts and elapsed cost. Do not silently retry, add
warmups, A/A runs, repeats or a median campaign. A job exceeding the limit
is returned to the master, not converted into an unrequested long job.
This one-sample diagnostic screen cannot establish a median wall gain or
promote a default. Verification-off results are explicitly uncertified.

Build reviewed source once in an isolated target with the established flags:

```sh
cargo build --release -p dft-cli --features mpi,mkl
```

Record source/diff hash, build command/exit code, compiler/MPI/MKL versions,
binary SHA-256 and the complete effective switch manifest. Build and launch
paths are filled from the existing approved dev08 manifests by the master.
Launch template (one invocation per job, `rmm_arm=off|on`):

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BAANDI_REFERENCE_SCF=on BAANDI_REFSCF_RMM="$rmm_arm" \
BAANDI_REFSCF_RMM_STEPS=4 BAANDI_REFSCF_RMM_BLOCK=16 \
BAANDI_DAVIDSON_STAGES=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
  /usr/bin/time -p taskset -c 0-7 \
  mpirun -np 8 --bind-to core --map-by core --report-bindings \
  "$baandi_binary" run "$frozen_input" --output-dir "$output_dir"
```

Archive the expanded command, binding report, raw stdout/stderr, external
wall and per-rank results. The verify-off input differs only by that switch;
record both hashes and a semantic diff. Profiling is ON identically for D/M.
No profiled BAANDI/QE wall ratio is permitted.

### Required tables and checks

For every arm, include external delivered wall, setup, ordinary maps,
verification/recovery, output, map count, verification attempts/failures,
H columns, rank-0 and summed-over-ranks GEMM FLOPs, occupied eigenvalues,
minimum GS norm and `E_delivered - E_strict` (Ha). Use the same existing
independently converged strict reference for D/M, with its full provenance.
If no matching strict reference exists, state unavailable; do not substitute
the D result or rounded QE output. Also report `E_M - E_D` (Ha).

| Case / MPI ranks | Verify | D wall (s) | M wall (s) | D/M maps | D/M H columns | D/M GEMM FLOPs | D/M E-E_strict (Ha) | M-D E (Ha) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| si64 / 8 | on | pending | pending | pending | pending | pending | pending | pending |
| si128 / 8 | on | pending | pending | pending | pending | pending | pending | pending |
| si128 / 8 | off | pending | pending | pending | pending | pending | pending | pending |
| elongated / 8 | on | pending | pending | pending | pending | pending | pending | pending |
| rskew / 8 | on | pending | pending | pending | pending | pending | pending | pending |

For each case/verification setting add an exclusive stage table with seconds,
H columns, GEMM FLOPs, rank-0 and rank maxima, and ordinary versus verification
calls separated. D: H apply, borders, factor, reduced EVD, corrections,
restart, repair and return. M: the remaining scheduled/verification Davidson
work plus RMM rotation, dense solve, DIIS, precondition, H apply, line search
and Gram-Schmidt. List per RMM map `n + sum(step_columns)` and compare it
with its recorded total; it must not exceed 5n at the default four steps.
Keep nested operator FFT/nonlocal clocks separate; do not sum nested clocks
or different ranks' maxima as a critical-path wall.

Required correctness: converged; verification passes when enabled; no lost
occupied states; `|E_M-E_D| <= 1e-10 Ha` for this ticket; delivered true error
within the declared strict/standard class. Keep any failed gate visible.
One run is insufficient for repeatability or statistical adoption claims.

## QE-only diagnostic reference (same MPI8 model)

Request **eight additional QE invocations**: one each of Davidson and
`diagonalization='rmm-davidson'` on each of si64/si128/elongated/rskew, using
QE 7.6 and the same frozen physical conditions and MPI8 binding. Only the
solver selection differs within each QE pair; RMM defaults are depth 4,
GS block 16 and unconverged one-sweep acceptance. Keep K_POINTS gamma and
the existing matched Davidson depth/band count/cutoff/FFT/start/conv_thr.
Archive the QE version/build, binary/input/UPF SHA-256, expanded command,
raw outputs, iteration/accuracy lines, stage timers and external wall.
Use existing installation only; no QE source access is needed.

Report QE Davidson versus QE RMM in a **separate QE-only table**, including
SCF iterations, per-iteration average solver iterations, E-E_strict and
wall. Label sample count one, no median or speed conclusion. QE's printed
average iteration count includes nominal GS work and is not a direct H
column count. These QE runs have no independent BAANDI verification step.
Do not form BAANDI/QE ratios from the profiled BAANDI runs. A subsequent
production-output, equal-conditions median campaign is a separate request.
