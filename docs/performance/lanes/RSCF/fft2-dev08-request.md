# P3-FFT2: fused last-axis rows, development screen

Status: default-off implementation ready for review; requested Mac checks pass.
The additional numerics gate fails a pre-existing GEMM-default assertion
(documented below). No speed claim. Base `8d1f3584eaa634f24b2bd2d11a1c1c341489919f`, changes uncommitted.
Master execution through `baandi-dev08`; no implementer remote changes.
Shared dashboard/session records, README generation, ledger/log integration
and commit remain the master's work, outside this ticket's OWNED_FILES.

## Mechanism and changed premise

[Pass-by-pass accounting](../../../design/reference-scf.md#p3-fft2-pass-ledger-and-row-fused-experiment)
separates every buffer read/write and the ideal one-transpose interface from
what this patch actually removes. The supplied `9e99bd3` Si128 MPI8 profile
attributes approximately 4.9 s of 10.2 s FFT time to kernels, with the rest
largely movement. Already-batched alltoall is not batched again. P3-FFT's
neutral direct-pack/sparse-clear result and PFFT-4's per-lane overhead remain
counterevidence: removed bytes are not a measured speedup.

`BAANDI_FFT_FUSED_ROWS=on` performs last-axis inverse -> V multiply -> forward
on one reusable row before gathering back into middle. Immutable inverse
input retains zero holes indefinitely. No slab traffic on ordinary apply;
no repeated slab clear. The first two axes, partition, normalization,
coefficient arithmetic, collective counts and padded payload are unchanged.
Kept-plane and standalone plane/density operations retain the slab route.
The new row route applies to all ordinary distributed local applications,
including reference SCF; it is not a reference-driver-only algorithm.

Exact removed interface bytes per pair/rank: `16*p_r*n1*n2` (sparse clear off),
or `16*p_r*(n1*n2-Z*n1)` relative to sparse clear on. Remaining middle/row
copies still move the same logical elements; reuse of short rows can reduce
cache misses, which needs measurement. The slab allocation is retained.
Candidate is conservatively rounding class: immutable and in-place RustFFT
interfaces are not guaranteed to choose identical arithmetic on every CPU.
No tolerance, stopping rule, occupation, grid or default is changed.

Derived from own-words notes fft-call-structure, qe-hpsi-and-fft and
qe-gamma-real-path, and BAANDI's implementation. No QE source access.
A26/PFFT-3/P39A-1..4 motivate local turnaround fusion; this reduces its
working set from a plane to a row. A11/PFFT-4 warn that smaller lane calls
can lose. M2_FFT is not revived. MKL lanes explicitly selected by M2_FFT=mkl
fall back to slab and record the reason in the switch manifest; an MKL build
with **M2_FFT=rust**, the current default, runs the candidate.

## Exactly eight profiled executions

The ticket's **development protocol** governs this request: **ONE profiled
run per arm per case**, no warmup executions, repeats, A/A runs, crossed
pairs, medians or precise timing campaign. Use one binary for both arms.
This overrides the older request's repeated-run recipe. Any later precise
stacked-default campaign is separate work, not part of this request.

| Case, model | A: FFT_FUSED_ROWS=off | B: FFT_FUSED_ROWS=on |
| --- | --- | --- |
| si64, MPI8, one thread/rank | 1 profiled run | 1 profiled run |
| si128, MPI8, one thread/rank | 1 profiled run | 1 profiled run |
| elongated, MPI8, one thread/rank | 1 profiled run | 1 profiled run |
| rskew, MPI8, one thread/rank | 1 profiled run | 1 profiled run |

Use the original converging reference-driver input manifests from
`/data/sig-dft/vc/` / the current P3 campaign, including elongated and rskew.
Record their actual absolute paths and hashes before execution. Labels are
not input specifications. Preserve cell/coordinates, cutoff, FFT grid, UPF,
band count, start, threshold, verification, output policy and all other
switches. In both arms force M2_FFT=rust, FFT_DIRECT_PACK=off and
FFT_SPARSE_CLEAR=off to isolate this mechanism. Preserve the adopted
M2_STICKS partition and reference solver defaults identically.

Master build, after review, with the same established native flags as the
profiled baseline:

```sh
CARGO_TARGET_DIR=target/p3-fft2 cargo build --release -p dft-cli --features mpi,mkl
```

Record build exit, compiler/MPI/MKL versions, complete source commit and
patch hash, binary SHA-256, exact command/environment, and input/UPF/start
hashes. Use separate output directories, no overwriting prior evidence.
Execute the same reference-driver command used by the current P3 screen;
only change the candidate switch and output directory. Invocation skeleton
(the master resolves the existing manifest paths, not invented inputs):

```sh
BAANDI_FFT_FUSED_ROWS=off BAANDI_M2_FFT=rust \
BAANDI_FFT_DIRECT_PACK=off BAANDI_FFT_SPARSE_CLEAR=off \
BAANDI_DAVIDSON_STAGES=on OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 8 --map-by core --bind-to core --report-bindings \
  <reviewed-binary> run <frozen-reference-driver-input-for-arm-A>
```

Arm B uses `BAANDI_FFT_FUSED_ROWS=on` and the separate output destination.
Use the established socket-0 physical-core 0-7 placement and confirm actual
binding; `--map-by core` alone is not proof of that placement. No overlapping
builds or measurements. One case/arm per fast job, at most 600 s. Preserve
verification and recovery inside delivered wall; never truncate a run to
manufacture a convergence result. A timeout/failure is retained with its
exit code and reason, without an automatic replacement run.

## Required outputs and interpretation

Save raw per-rank results, merged maxima/minima, stdout/stderr and exit codes.
For each case/arm report:

- External delivered wall, convergence and independent verification, energy
  in Ha, density/eigenvalue fingerprints, strict-reference error with its
  existing provenance, and off/on energy difference. Strict 1e-10 Ha and
  standard 5e-9 Ha classes remain; unique-solution 1e-12 Ha agreement remains
  an internal check. If bits differ, report rounding explicitly.
- SCF maps, Davidson iterations, complete H columns-per-call sequence,
  two-band pairs, odd tails, retry/verification calls and attempts. Expect
  unchanged counts on bit-identical cases; any difference must be reported,
  not attributed to a faster FFT.
- Outer H total and all FFT children: scatter; inverse first axis; transpose
  pack; inverse alltoall; middle unpack; inverse second axis; slab scatter;
  inverse last axis; multiply; forward last axis; slab gather; forward second
  axis; middle pack; forward alltoall; transpose unpack; forward first axis;
  gather/scale. Retain nonlocal/kinetic/ORC/mirror controls and unassigned time.
- Per-stage points, logical bytes and calls. For stages 6..10, B reports row
  regions (n1 times the calls of A); point totals and transform counts remain
  equal. Only stage 6 bytes decrease, by the accounting formula. Per-call
  rounding of estimated FFT FLOPs can change the integer estimate slightly;
  use transform points/axis lengths for exact work comparisons.
- S, Z, s_r, p_r, n0/n1/n2, batch widths, useful/padded payload, and rank-local
  closure. Do not sum separately selected rank maxima. Keep the FFT and
  whole-H totals beside the individual children.

The finer row split adds diagnostic timer/counter calls. Its overhead can
increase unassigned H time and distort a profiled speed difference; disclose
this when interpreting the single-run result. With profiling off there are
no new clocks/TLS updates. The test verifies point/byte accounting separately
with clocks disabled. The observed speed direction is a development screen,
not a confidence interval or final performance claim. Keep the switch off
until the master evaluates the four cases; record negative/neutral results
as well as wins. No QE execution is requested. Any later cross-code table
must use MPI8 against MPI8 and complete reference provenance.

## Mac verification

Mac, pinned Rust 1.99.1, pure backend. All outputs below are correctness/work
counts, not performance acceptance. Raw logs are in `target/tmp/fft2/` of the
task worktree, not committed external-reference data.

| Command | Exit | Result |
| --- | --- | --- |
| `cargo check --workspace --tests` | 0 | Workspace/tests type-check |
| `cargo test -p dft-solvers` | 0 | 165 passed, 25 ignored (before the additional backend-selection test below) |
| `cargo test -p dft-physics` | 0 | 95 passed, 18 ignored |
| `cargo test -p dft-engine --test refscf -- --nocapture` | 0 | 19 passed, 1 ignored |
| `cargo fmt --check` | 0 | Format gate |
| `cargo clippy -p dft-solvers -p dft-physics --all-targets -- -D warnings` | 0 | No warnings |
| `cargo test -p dft-solvers --lib fused_rows -- --nocapture` | 0 | Both new tests pass, including backend selection |
| `cargo test -p dft-physics --lib fft_pack_switches -- --nocapture` | 0 | 4,200 rank-local H-bit comparisons (1,680 exercise the new switch), maximum difference 0 Ha |
| `cargo test -p dft-cli --test all run::expressions_preserve_energy_and_record_resolutions -- --exact` | 0 | Manifest count 144 and output record pass |
| `CHECK_SCCACHE=0 BAANDI_TEST_TMP="$PWD/target/tmp/fft2" CHECK_TEST_JOBS=2 scripts/check.sh dft-numerics` | **1** | Build/fmt/clippy PASS; tests 88 passed, 1 failed, 22 ignored: pre-existing stale GEMM-default assertion |

The numerics gate failure is `tests/gemm_dispatch.rs:167`, which asserts
`!default.on(Id::REF_DAV_SINGLE_GEMM)` although that switch's default is on.
Both the assertion and the default-on registry entry are present unchanged
in base HEAD `8d1f3584`; `git show HEAD:<path>` confirmed them. This ticket
only appends the FFT switch and updates the two permitted switch counts.
The failing GEMM test is outside ownership and is left unchanged for the
master; **the broad numerics gate is not reported as passing**. The gate's
optional `ps` load diagnostic and `nice` priority adjustment were denied by
the sandbox; build, tests, fmt and clippy still executed. The new registry
checks passed within the numerics test suite.


The solver test makes **1,260 rank-local paired-transform comparisons**:
9 supports/grids (cubic 6^3, elongated 5x7x15, shifted/asymmetric 5x7x11,
forced y-first 3x5x7, empty 1x2x3, sparse 3x3x17/36/72 and full 3^3), ranks
1/2/3/8, sparse clear on/off, and batch widths 3/1/5/0/2. It separately
compares inverse rows with the original public inverse-plane API, checks
V=1 round trip, signed-zero/random potentials, direct-input fusion, dirty
reuse, stable scratch pointers, untouched slab sentinels and immutable zero
holes. Observed maximum on/off relative difference **0**, zero non-bitwise
arrays/rows. Stage point counts match; removed stage-6 bytes equal the exact
formula for both sparse-clear settings. The rounding cap (2e-14 relative)
was fixed before executing the test; it was not relaxed.

The physics test adds FFT_FUSED_ROWS alone and stacked with the P3-FFT
switches to the existing H-bit oracle: cubic/skewed/elongated cells, ranks
1/2/3/8, greedy/mirror ownership, batched/unbatched, band widths 5/1/2/3/0.
Production physics code and the concurrently owned profile counters are
untouched; only this existing test's FFT configuration arms change.

A standalone scratch test built against the current engine rlibs reuses the
unchanged `refscf.rs` Gaussian fixtures and `stop_to_verified_science`
comparison (energy, density/orbital bits, eigenvalues, histories, verification,
map/solve records and call counts with clocks removed). Both FFT arms are
explicit RunConfig entries, not environment overrides hidden by fixture
scopes. It asserts switch consultation, convergence and passed verification
on every rank. Reproducer and exact rustc link command are retained in
`target/tmp/fft2/build_scf_probe.py` and `scf_probe_build.log`; execution:
`target/tmp/fft2/refscf_fft2 fft2_reference_fixture_identity --exact --nocapture`
(exit **0**, output `scf_on_off.log`). No engine source/test file was edited.

| Synthetic fixture | Ranks tested | Maps / solves / verification attempts in both arms | Maximum on/off energy difference |
| --- | --- | --- | --- |
| Local Gaussian (z=2, rc=1.8 Bohr, skew length 7 Bohr, 2 bands) | 1/2/3/8 | 6 / 7 / 1 | 0 Ha, all compared science/count bits equal |
| Nonlocal Gaussian (z=4, rc=1.2 Bohr, skew length 11 Bohr, 3 bands) | 1/2/3/8 | 28 / 30 / 2 | 0 Ha, all compared science/count bits equal |

At rank 1 the totals are respectively -1.10503701485504457 Ha and
-2.61337340593422773 Ha. These synthetic energies are **not QE references**.
Different rank counts retain the existing reduction-rounding differences;
identity here is on/off at each fixed rank count. The real Si external-file
SCF test stays ignored on the Mac. Actual MPI, MKL lanes, x86 RustFFT rounding
and the four material performance screens remain remote validation work.
