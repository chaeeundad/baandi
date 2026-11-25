# P4-DAVREM: Davidson remainder, one-run dev08 screen

Status: implementation candidates, default off; native MKL and timings pending.
Base checkout: `50333c72a8f2061359cad7a80b0f014dc9eed7eb`; uncommitted task diff.
The master reviews/commits and relays this request to `baandi-dev08`.
No implementer server access, installs, new QE reference or QE source access.

## Question and prior evidence

Ticket-relayed fa4cfad MPI8 snapshot: si216 Davidson remainder 25.69 s versus
matched QE regterg minus h_psi minus rdiaghg 19.2 s; si128 6.46 versus 5.21 s.
These are historical diagnostic boundaries, not new measured speed claims.
Stage seconds (rank maxima) for si216: restart 5.12, S border 4.94, H border
4.60, return 3.89, corrections 3.41, setup 2.11, transform 1.23,
coefficients 0.83, factor 0.78, normalize 0.76, precondition 0.55,
apply_check 0.52. Obtain the original command/input/UPF/binary hashes before
reusing those numbers in a comparison. New work is explained in
[Inner solver cost accounting](../../../design/reference-scf.md#inner-solver-cost-accounting).

The old upper-border split lost 0.57% despite 16.38% fewer FLOPs; it added
12.7% calls and reduced DGEMM throughput. New premise explicitly authorized
by this ticket: a square upper-only GEMMT kernel, not narrower DGEMM tiles.
Its rectangle plus triangle still adds a call and can regress. The stacked
rotation adds packing traffic; measure it, do not predict its benefit from
call count alone. S rank validation remains before H, including failed appends.

## Arms: exactly ten profiled SCF executions

Use frozen **si128 and si216 MPI8**, **verify off**, one thread/rank. Freeze
identical cell/coordinates, cutoff, FFT grid, UPF, bands, start, threshold,
placement, all adopted switches and binding across arms. Do not revert the
adopted P3 SINGLE_GEMM, WIDE_ROTATE, IDEAL_RESET or STREAM switches. The only
arm differences are these three new switches:

| Arm | REF_DAV_TRIANGLE | REF_DAV_STACK_ROTATE | REF_DAV_SETUP_FUSED |
| --- | --- | --- | --- |
| D | off | off | off |
| T | on | off | off |
| R | off | on | off |
| S | off | off | on |
| TRS | on | on | on |

**ONE run per arm/case**, ten total, no warmups, repeated pairs, A/A series,
automatic retries or new QE runs. Each job <=600 s; retain failed jobs and
elapsed cost. Follow the full verify-off stopping condition, not an iteration
truncation. This is the user-requested development profile, not a verified
production-speed or multi-system adoption campaign.

```sh
cargo build --release -p dft-cli --features mpi,mkl
cargo test -p dft-mkl --features mkl --test dgemmt
cargo test -p dft-solvers --features dft-numerics/mkl p4_ -- --nocapture
```

The linked tests are correctness gates, not additional SCF screen repetitions.
Use the same supported MKLROOT and toolchain as the frozen snapshot; check the
GEMMT symbol and adapter contracts before starting the ten runs. Freeze absolute
binary, verify-off input and unique output paths, then execute each arm:

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BAANDI_REFERENCE_SCF=on BAANDI_MKL_DGEMM=on \
BAANDI_REF_DAV_TRIANGLE="$triangle" \
BAANDI_REF_DAV_STACK_ROTATE="$stack_rotate" \
BAANDI_REF_DAV_SETUP_FUSED="$setup_fused" \
BAANDI_DAVIDSON_STAGES=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
  /usr/bin/time -p taskset -c 0-7 \
  mpirun -np 8 --map-by core --bind-to core --report-bindings \
  "$baandi_binary" run "$baandi_input" --output-dir "$unique_output"
```

Verify the actual socket-0 physical cores and one-thread-per-rank binding;
record it, not merely the intended affinity. Preserve source commit/diff SHA256,
binary/input/UPF SHA256, versions/build flags, effective switch manifests,
load before/after, completion/exit codes, all-rank RSS, console, results.json,
density.bin, SCF and Davidson histories. Do not silently reuse a mismatched
historical QE setup. Only MPI8 against MPI8 may appear in a QE comparison.

## Required tables and gates

For each cell and arm, publish every exclusive reference stage, including
setup, apply_h, apply_check, append, border_s, border_h, factor, transform,
reduced_eigen, coefficients, corrections, precondition, normalize, repair,
repair_panel, restart, return, decisions and glue:

| Case/arm/stage | Max seconds (rank ID) | That rank's GFLOP | GFLOP/s | GEMM calls | Streamed GB | GEMM operand GB | Logical GB/s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Pending dev08 | pending | pending | pending | pending | pending | pending | pending |

Read `convergence.work.davidson_stages.work.per_rank`; compute rates on each
rank, then select the maximum-seconds rank's entire row. Do not divide max
FLOPs by an unrelated rank's max seconds. Keep the all-rank table. GFLOP is
(gemm_flops + scalar_flops)/1e9; streamed GB excludes GEMM operands; logical
GB/s includes both. EVD/H arithmetic is opaque, so zero modeled FLOPs are
not a claim of zero work. Stage rank maxima need not add to critical wall.
Report wall, total remainder, H columns/calls, iterations, reduced solves,
restarts, repair/drop counts, active-root decisions, energy Ha and density
identity/difference next to D. Show both old QE-matched clocks and new exclusive
stage sum with their distinct boundaries. No speedup claim from Mac clocks.

Expected movement (executed shapes, no assumptions about Si band counts):

- T: per border pair save `2*r*k*(k-1)` FLOPs; two calls become four for
  nonempty old basis, remain two on a full square. `_gemmt` shape entries
  use `m*(m+1)*k` FLOPs, not `2*m*n*k`. Stage/iteration FLOPs already use this
  convention; legacy dense dispatch omits GEMMT, adapter timer includes it.
  Adapter scatter: 32*m*k versus baseline 32*(m+k)*k bytes/pair (both excluded
  from stage bytes; report this shape-derived supplement separately). Counted
  lower zero fills: 8*k*(k-1) bytes/pair.
- R: each restart and HX-return pair keeps `4*r*m*n` FLOPs; two logical calls
  become one sequential MKL call (baseline wide kernels may make many panel
  calls). Packing/unpacking adds 32*r*(m+n) streamed bytes. Setup adds
  16*r*(capacity+n) initialized bytes/rank. Report peak RSS and whether call
  savings outweigh this extra traffic. X-only returns stay unchanged.
- S: remove 8*r*n + 8*r setup reads per solve with SQRT on, plus another
  8*r*n when supplied HX is present. Off accounting now includes its prior
  uncounted HX finite scan. Normalize/append/apply-check already contain P3
  fusion; list their remaining GB rather than claiming a duplicate saving.

Require equal active sets, H columns/calls, restart/repair/drop decisions and
convergence outcomes on fixtures and screen cases. R/T are rounding class;
use strict 1e-10 Ha / standard 5e-9 Ha classes against the existing independent
reference, report 1e-12 Ha internal agreement separately. S requires bits.
For any first divergent rank/stop branch, retain previous/current Ritz values,
threshold and margin, accepted/dropped correction indices, factor panel/rank
and pivot ratio versus `PIVOT_RATIO_MIN`; use the existing rank-rule diagnostic
replay through the master. Do not relax tolerances or adopt based on final
energy alone. Verify-off results do not provide a new fixed-point certificate.

## Local validation and limits

Local command results and numerical fixture evidence are recorded below after
execution. Native MKL execution is unavailable on this Mac; a successful
feature type-check is not a linked-kernel correctness result. Default-off
paths, pure fallback, cyclic partitions (including empty ranks), clustered
spectra, capacity exhaustion/repair, root/replicated reduced solves, cold/warm
starts, repeatability and 1/4-worker identity are covered by task tests.

The shared ledger, dashboard/session JSON and README generator are outside
OWNED_FILES. This document is the task-local attempt/evidence record for the
master to reconcile after review; no measured performance or external reference
result is being added to README by this ticket.

### Executed Mac evidence (2026-11-25)

Pure Rust task tests (`cargo test -p dft-solvers p4_ -- --nocapture`): **4 pass**.
The trajectory matrix covers 12 shape/rank combinations, two reduced placements,
D against T/R/S/TRS, repeated solves and worker counts 1/4. Maximum observed
candidate/off eigenvalue difference: **0 Ha**. Active sets, H calls/columns,
repair/drop/restart decisions and communication are equal. This does not test
native GEMMT: T follows the documented pure fallback on Mac. The 513-row
fixture exercises multiple row panels and an unconverged iteration-limit
return; the 7-row/4-band fixture forces capacity repair and includes an empty
rank at MPI-size 8. These are simulated collectives, not a real MPI run.

A separate integer scalar oracle checks stacked V/HV rotations, both returned
fields and X-only returns. For r=257, m=5, n=3 the tested stacked path has
**15,420 FLOPs**, **1 logical product** and **65,792 explicit packing/unpacking
bytes** (6.5792e-5 GB); FLOPs equal the two-product baseline. Setup traffic
assertions verify 8*r*(n+1) bytes saved without supplied HX and 8*r*(2*n+1)
with supplied HX, at SQRT on. Nonfinite start, supplied HX, diagonal and bad
HX shape on one rank are rejected collectively, with no completed-work record.
Padded distinct-operand GEMMT, untouched lower/padding sentinel, zero sizes,
LP64/stride/short-storage rejection tests type-check but require dev08 execution.

Verified command exits:

| Command | Exit | Result |
| --- | --- | --- |
| `cargo check --workspace --tests` | 0 | workspace/test type check |
| `cargo check -p dft-solvers -p dft-numerics -p dft-mkl --tests --features dft-numerics/mkl` | 0 | MKL feature type check, no native link |
| `cargo test -p dft-solvers` | 0 | 191 passed, 25 ignored on final source |
| `cargo test -p dft-engine --test refscf` | 0 | 60 passed, 7 ignored |
| `cargo test -p dft-mkl` | 0 | 3 non-native tests passed; linked adapter tests disabled |
| `cargo fmt --check && cargo clippy -p dft-solvers -p dft-numerics -p dft-mkl --all-targets -- -D warnings` | 0 | formatting and requested lints |

Tests use `BAANDI_TEST_TMP=$PWD/target/tmp TMPDIR=$PWD/target/tmp` where needed.
The full gate `BAANDI_TEST_TMP="$PWD/target/tmp" CHECK_SCCACHE=0 scripts/check.sh`
returned **1**: build/fmt/clippy PASS, test FAIL (**947 passed, 3 failed,
198 ignored**, 96 binaries). Its unchanged pseudo-cache tests reject a cache
under a repository (`pseudopotential cache must be outside repositories`):
`cache_resolution_and_no_network_side_effects`,
`cached_fetch_checks_hash_and_provenance`, and
`fetch_boundary_verifies_before_publication_and_never_overwrites`.
The writable external-volume scratch is inside this worktree; the normal
external non-repository scratch is outside sandbox write permission. No cache
policy, tests or tolerances were changed to conceal these failures. The gate's
optional load-report `ps` probe is also denied by the sandbox; the actual four
gates still execute. Whole-gate wall 352 s is not a numerical speed measurement.
After the final accounting-only adjustment (exclude adapter scatter equally
on both arms), workspace/MKL type checks, the 191-test solver suite and requested
clippy were repeated successfully.

Final changed-crate gate:
`BAANDI_TEST_TMP="$PWD/target/tmp" CHECK_SCCACHE=0 scripts/check.sh fast dft-solvers dft-numerics dft-mkl`
returned **0**, build/test/fmt/clippy all PASS: **287 passed, 0 failed,
47 ignored** in 28 binaries. Extra
`cargo clippy -p dft-solvers -p dft-numerics -p dft-mkl --all-targets --features dft-numerics/mkl -- -D warnings`
returned **0** on final source (MKLROOT unset; type/lint only).

To isolate the unrelated full-gate failure, the existing linked CLI test binary
`target/debug/deps/all-b0a793b818630f7a pseudo:: --nocapture` was run with the
same `BAANDI_TEST_TMP`/`TMPDIR`: exit **101**, 4 pass and the two named cache
success tests fail with exactly `pseudopotential cache must be outside
repositories`. The source-local third failure and same message are in the
full-gate output. Master should rerun the full gate with its permitted
external non-repository scratch. No worktree files outside OWNED_FILES were
modified, and all changes remain uncommitted.
