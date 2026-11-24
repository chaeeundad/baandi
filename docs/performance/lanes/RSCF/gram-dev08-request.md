# P3-GRAM: single-call reference Davidson Gram products

Status: implemented, **default off**, native MPI/MKL and performance acceptance
pending. Run through the master session and `baandi-dev08`; no remote actions,
installs, timing campaign or default promotion were performed on Mac.
Base worktree HEAD: `d7ec4e9565a24d2f88c5b241e4fdd597f50d2e73`.
Changes are uncommitted for review. No QE source was accessed.

## Contract and changed premise

Derived from the own-words [Davidson GEMM accounting note](../../../research/qe-study/davidson-gemm-accounting.md)
and engine inspection, especially [inner solver cost accounting](../../../design/reference-scf.md#inner-solver-cost-accounting).
The ticket supplies a dev08 Si128 MPI8 profile at `66d45e69`, W/R/S on:
border 4.43 s, corrections 2.73 s, `other/mkl` 2353 logical products and
328 GFLOP per rank, local reduction length about 5981 (24 panels). These are
relayed motivation, not measurements reproduced here. Archive the original
commands, raw files, source/binary/input/UPF hashes before citing that profile.

Fixed 256-row AB and 256-term A^T B panels provide task granularity independent
of worker count, ascending-order partial sums, and bounded scratch. The
workspace reserves ceil(long/256) tiles of short^2 doubles. The old pure nano
kernel also bounds library calls to avoid opaque packing allocations. The
wide entry relaxes only the AB whole-output bound; it does not remove TN
reduction panels. `M2_DGEMM` belongs to the other distributed solver and does
not select the reference workspace's TN calls.

Ledger review: A49 rejected finer tiles/cutoff tuning after trajectory changes;
P32A-3 rejected allocating generic pure GEMM; DGEMM-SITES identified call-size
inefficiency; M2_DGEMM/P49A addressed wide dispatch, retaining panels. The new
premise is the measured reference-MPI8 TN bottleneck with one worker per rank,
sequential MKL, and an explicitly authorized rounding-class experiment. This
change does not retune panels, use generic allocating faer GEMM, fuse H/S,
change upper-border work, or revisit the rejected exact-H/S restart reset.

`REF_DAV_SINGLE_GEMM` resolves once from RunConfig at solve entry. Only the
reference workspace's border products (including initial/refresh Grams) and
the two repair overlaps V^T E select the new entry. Rotations and correction
back-application retain their current W route; factor/reduced/coefficient
products and all other callers are unchanged. MKL performs one sequential
`dgemm_at` over all local rows, into the first existing scratch tile, then
scatters into the output view. No new scratch/packing of operands is added.
Full output must fit short^2, as these reference TN rectangles already do.
Empty owners write zero, count as `pure_layout`, and make no vendor call. Native MKL internal packing
is library-owned; this is no new engine allocation, not a promise about opaque
vendor allocations. Pure builds and `MKL_DGEMM=off` keep the old pure arithmetic.

The existing link uses `mkl_sequential`; there is no worker-dependent reduction
or native threading choice in the new call. Same backend/build/layout must
repeat bitwise, including across worker budgets. Across rank counts the
[reference rank acceptance rule](../../../design/reference-scf.md#inner-solver)
allows rounding, with documented threshold straddles only: record the first
branch, both Ritz operands/changes/thresholds and signed margins in Ha;
`32*EPSILON*max(1 Ha, abs(four Ritz operands))` is the existing envelope.
Unexplained earlier count/trajectory differences fail. Every delivered result
must still converge, certify, conserve charge within 1e-10 e and meet strict
1e-10 Ha / standard 5e-9 Ha against independent tight oracles.

## Build, arms and identity

| Arm | BAANDI_REF_DAV_SINGLE_GEMM | Purpose |
| --- | --- | --- |
| D | unset, no input override | Default-off identity control |
| A | off | Existing panels |
| B | on | One sequential MKL call per TN product |

Input key: `[performance] ref_dav_single_gemm = 'on'`. Registry appended at
index 140 (141 entries), with normal environment precedence, consulted state
and MPI agreement hashing. Pin all other effective settings, including the
same W/R/S settings as the supplied profile, and all verification switches.
No change in default main since the profile may be attributed to this switch.

Build profile baseline `66d45e69`, current main, and reviewed candidate in
separate target directories using the established dev08 toolchain:

```sh
CARGO_TARGET_DIR=target/gram cargo build --release -p dft-cli --features mpi,mkl
CARGO_TARGET_DIR=target/gram cargo build --release -p dft-cli --example davidson_bench --features mpi,mkl
```

Record full source commits/diff hash, compiler/MPI/MKL versions, exact build
commands, exit codes and binary SHA-256. Master runs `scripts/check.sh full`
and the reference rank matrix on dev08, serialized with measurements. Require
pre-change main versus candidate D/A bit-identical science/work; old callers
must remain identical even when B is enabled. New manifest entry is excluded
from that comparison, not scientific state.

Use the exact **si64, si128, elongated, rskew** inputs from the existing W/R/S
MPI8 campaign, with frozen start/density snapshots. Master resolves and records
absolute paths and SHA-256 before launch; they are not available locally and
must not be guessed. Match lattice/coordinates, UPF bytes, cutoff, FFT grid,
bands, start, threshold, placement, mixing, verification, outputs and binding.

## Untimed screen, then A/A and crossed A/B

1. Untimed correctness/counter screen D/A versus B on all four cases, MPI8,
   one thread per rank, on socket 0 physical cores 0-7. Include cold, warm,
   late/verification reference fixed-H solves with the same reconstructed
   operator and captured start. Validate A replay first; capture potential and
   row-layout identity. Use the [reference replay interface](dav-dev08-request.md),
   `--ref-arm new --replay-mode ab --ref-profile --ref-hx --repeats 1` with
   explicit captured tolerance/capacity/iteration budget. Never reinterpret
   the harness's legacy arm as the panel baseline. For a snapshot unsupported
   by that interface, use an established reference capture/replay harness or
   report the limitation; do not silently change the functional/input.
2. Run separate profiled full-convergence screens on si64/si128/elongated/rskew
   at MPI8 (D/A versus B) with independent verification retained. These stage
   seconds are diagnostics, not acceptance timings. Retain per-rank records,
   all decisions, and the table below. Repeat B twice for fixed-layout bit
   fingerprints, stripping clocks only. Check the native rank matrix 1/2/3/8
   with B and the thread-budget 1/2/4/8 kernel/solver tests on MKL. Empty owners,
   clustered spectra, repairs and non-multiple-of-256 local lengths matter.
3. If correctness passes and the screen is not slower, start one representative
   case (Si128) in an exclusive window: disposable warmup, three A0/A1 crossed
   pairs for the noise floor, then three A/B crossed pairs (AB, BA, AB, followed
   by the reverse order in the next block). No builds overlap. Widen survivors
   to the other three cases. Close effects require 30 crossed pairs per case
   and at least 20 A/A pairs interleaved across the campaign. Keep raw samples,
   spreads, load, bindings and all failures. No best-sample selection.
4. Profile-free acceptance runs use `BAANDI_DAVIDSON_STAGES=off`; preserve all
   other effective settings. Compare medians of external delivered wall,
   including verification, failures and recovery. Show verification separately.
   Full convergence and multi-system validation are needed for promotion;
   fixed-H/local tests alone cannot promote this default. A reproducible small
   improvement is eligible; there is no percentage rejection threshold.

Launch template (substitute frozen paths and established binding wrapper):

```sh
BAANDI_REFERENCE_SCF=on BAANDI_REF_DAV_SINGLE_GEMM=off \
BAANDI_DAVIDSON_STAGES=off OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 --bind-to core --map-by core --report-bindings \
  /absolute/reviewed/dft-cli run /absolute/frozen/input.toml \
  --output-dir /absolute/unique/A-output
```

Pin the launcher to the established socket-0 physical core set; verify actual
bindings rather than assuming these flags enforce socket selection. Replace
only the new switch and output path for B. One case/rank count per job, fast
class <=600 s; if full convergence cannot fit, request an authorized named
`final-gram-<case>-<arm>` job <=7200 s through the master. Do not truncate tails.

## Required stage and dispatch tables

For each case/arm, report per-rank values and critical-rank closure. Do not sum
maxima belonging to different ranks. Split ordinary/verification and retain
their combined totals.

| Stage / counter | A | B | Unit / expectation |
| --- | --- | --- | --- |
| apply_h, border, factor, reduced_eigen | pending | pending | s, each separate |
| corrections, repair, restart, return | pending | pending | s, each separate |
| external delivered / verification wall | pending | pending | s, separate |
| SCF/maps, verification attempts, solver calls, iterations | pending | pending | counts |
| H calls/columns, restarts, repairs, dropped corrections | pending | pending | same on fixed-H fixtures |
| GEMM FLOPs; collectives and payload | pending | pending | FLOPs, calls, bytes; unchanged shapes/work |
| `other/mkl`, `other/mkl_single` | pending | pending | calls, FLOPs, max_m/max_n/max_k |
| every other `gemm_dispatch` route (including `rotate/mkl_wide`) | pending | pending | retain all nonempty rows |

`gemm_dispatch` counts **logical wrapper products**, not internal panel calls.
Thus B moves TN products from `other/mkl` into `other/mkl_single`, at unchanged
2*m*n*k FLOPs; one nonempty `mkl_single` count is exactly one vendor call.
A's physical panel-call count is sum(ceil(k/256)) over its TN shape trace, not
2353 alone. Retain raw shape records if calculating eliminated backend calls.
Pure Mac `rows_as_if_mkl` is a dispatch prediction, not executed MKL evidence.
Use optional vendor-call timers only in separate diagnostics if needed.

Archive commands, stdout/stderr, results JSON, density/orbitals, scientific
fingerprints, Ritz traces, energy/charge errors, input/UPF/binary hashes and
raw locations. If comparing QE, only MPI8 versus MPI8 with identical conditions
and the same median/external-wall boundary is valid. No new QE run or speed
ratio is claimed by this change.

## Local validation and integration handoff

Executed Mac checks and synthetic numbers are recorded below.
The master owns dashboard/session/README generation, the attempts ledger and
optimization-log integration (outside this ticket's OWNED_FILES), full dev08
gates, review and commit. No recurring session-report timer is active here.

Ledger draft: **P3-GRAM**, reference-only sequential full-local-row MKL TN
entry, default off; derived from the GEMM accounting note and current profile.
New premise: panel overhead persists after wide rotation dispatch was fixed.
No native speed result yet; document accepted/rejected dev08 outcome and lesson
without changing thresholds or suppressing any failed case.


### Executed Mac checks (2026-11-24)

All commands below completed with exit 0 on the final source. Log files are
under `target/tmp/p3-gram-*.log` (local scratch, not committed). No Mac timing
is used as performance evidence.

| Command | Result |
| --- | --- |
| `cargo check --workspace --tests` | Workspace and tests compile, including the final test additions. |
| `cargo test -p dft-numerics` | 89 passed, 22 ignored, 0 failed. |
| `cargo test -p dft-solvers` | 163 passed, 25 ignored, 0 failed (final rerun). |
| `cargo test -p dft-engine --test refscf` | 19 passed, 1 ignored, 0 failed. |
| `cargo fmt --check && cargo clippy -p dft-numerics -p dft-solvers --all-targets -- -D warnings` | Both pass. |
| `cargo check -p dft-numerics -p dft-solvers --tests --features mkl` | MKL branch type-check passes; MKLROOT unset, no link/runtime claim. |
| `cargo test -p dft-numerics --test gemm_dispatch -- --nocapture` | Dispatch counts, scalar oracle, padded/strided bounds and repeat/thread checks pass. |
| `BAANDI_REF_DAV_SINGLE_GEMM=on target/debug/deps/ref_davidson-4f96658c9cc927aa --nocapture` | All 5 existing independent Davidson oracle tests pass. |
| `BAANDI_REF_DAV_SINGLE_GEMM=on target/debug/deps/refscf-601b960de0af45e8 --nocapture` | 19 pass, 1 ignored with the switch enabled. |
| `git diff --check` | Pass. |

The executable suffixes are build-specific. Cargo equivalents for the three
binary reruns are `cargo test -p dft-numerics --test gemm_dispatch -- --nocapture`,
`BAANDI_REF_DAV_SINGLE_GEMM=on cargo test -p dft-solvers --test ref_davidson -- --nocapture`
and `BAANDI_REF_DAV_SINGLE_GEMM=on cargo test -p dft-engine --test refscf -- --nocapture`.
One initial `cargo test -p dft-solvers` exited 101 because the new test attempted
to clone a non-Clone Record. Fixed by cloning its existing Clone iteration
records and comparing all work/decision fields with clocks removed. This was
a test compilation error, not a numerical failure; final full rerun passes.

Synthetic evidence:

- Dispatch: **40 nonempty single-entry products, 18,142,176 nominal FLOPs**;
  max scalar-oracle error **1.3642420526593924e-12**, dimensionless (gate 1e-11).
  Local reduction lengths 0/1/255/256/257/1025, padded leading dimensions,
  strided output sentinels, invalid storage/stride, NaN rejection and the
  MKL_DGEMM opt-out are covered. Both repetitions at T=1/2/4/8 are bit-identical;
  the 1025x40x27 case exceeds the original parallel dispatch threshold.
  The Mac `mkl_single` row is **as-if-MKL**, not executed vendor calls.
- Solver: **24 layout/placement configurations, 96 rank-local off/on pairs**:
  dimensions/bands 12/3, 17/5, 7/4, 1025/3; cyclic ranks 1/3/8; root and
  replicated EVD. The small cases run up to 40 corrections, the large shape
  uses an explicit two-iteration budget. Every pair preserves all decisions,
  H calls/columns, GEMM shapes/FLOPs, collectives, repairs/drops/restarts and
  reduced orders. Pure off/on eigenvalues, vectors, HX and Ritz history are
  bit-identical (**0 Ha eigenvalue difference**). Repeated on solves and
  T=1 versus T=4 have identical scientific bits and decision records.
- Explicit repair fixture (257 rows): two full-block projection products,
  one retained and one dropped correction, **0 H calls / 0 H columns**; old/on
  and repeated-on correction bits and work agree.
- Existing dense/legacy oracle with on: unclustered **41 iterations / 233 H
  columns**, clustered **48 / 237**, at ranks 1/2/3/8 and both EVD placements.
  Maximum eigenvalue error **7.327e-14 Ha**, occupied projector error
  **2.392e-9**, dimensionless, under the unchanged oracle tolerances.

These establish pure fallback isolation and synthetic correctness, not MKL
rounding or native MPI repeatability. The same MKL-enabled tests and full
rank/accuracy matrix are required on dev08 before adoption. The default is
still off. No full check.sh, release rank campaign, remote measurement,
dashboard/README/ledger mutation or commit was performed in this worktree.
