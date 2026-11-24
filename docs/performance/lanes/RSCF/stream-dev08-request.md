# P3-STREAM: reference Davidson streaming passes

Status: four default-off experiments; Mac correctness and counter checks only.
Master reviews and commits; `baandi-dev08` executes the request. No remote
changes or runs were made here. Derived from the preconditioner addendum of
[the own-words iteration note](../../../research/qe-study/davidson-iteration-structure-gamma.md);
no QE source was read. The implementation uses this engine's storage and
collective contracts. Defaults are not proposed for promotion before dev08.

## Motivation and independent arms

Ticket-relayed P3i profile, 0f11a240, si128 MPI8, componentwise rank maxima,
20 solves: precondition 1.21 s (54 million hypot calls, 0.70 GFLOP), setup
0.96 s (3.5 GB logical traffic), normalize 0.31 s, append 0.27 s,
apply_check 0.21 s. These are diagnostic observations, not new measurements
or a sum of simultaneous rank-local critical paths. No profiled BAANDI/QE
speed ratio is inferred. The new premise is the measured scalar/streaming
floor after the earlier GEMM and restart-copy improvements, not a retry of
rejected capacity, Gram approximation, or reduction-order experiments.

All names in the table have prefix `BAANDI_REF_DAV_` in the environment;
lower-case `ref_dav_*` names work under `[performance]`. Explicitly set all
four in every arm. All default **off**.

| Arm | SQRT | FUSED_APPEND | CHECK_APPEND | REUSE_ROWS |
| --- | --- | --- | --- | --- |
| A: baseline | off | off | off | off |
| B: vectorizable sqrt/residual/norm | on | off | off | off |
| C: normalize and append | off | on | off | off |
| D: H finite check and append | off | off | on | off |
| E: initialized row storage reuse | off | off | off | on |
| F: combined | on | on | on | on |

B/F are rounding class; C/D/E preserve arithmetic and row order bit for bit.
B keeps the rationalized denominator in Ha with x=2(h-e), rather than using
the cancellation-prone literal expression. The exact formula, rounding
bound, extreme-column fallback, traffic counts and remaining passes are in
[the design contract](../../../design/reference-scf.md#p3-stream-contiguous-kernels-and-fused-stores).
There is no threshold, band-policy, reduction-tree or H-work change.

## Exactly one profiled execution per arm and case

| Case | Model | Executions |
| --- | --- | --- |
| si64 | MPI8, one thread/rank | A through F, once each |
| si128 | MPI8, one thread/rank | A through F, once each |
| elongated | MPI8, one thread/rank | A through F, once each |
| rskew | MPI8, one thread/rank | A through F, once each |

**24 executions total.** No warmups, repeated samples, A/A, crossed pairs,
medians or automatic replacement of failures in this ticket's development
request. Preserve failures/timeouts and their exit codes. One case/arm per
600 s fast job, `p3-stream-<case>-<arm>-mpi8`. Do not truncate convergence or
extend a job silently; an overrun is recorded for the master. A subsequent
precision/adoption campaign is separate work, using the current protocol.

Use full reference-SCF convergence with independent verification and all
retry/recovery costs inside delivered wall. This both reaches repeated solves
for E and detects rounding-induced trajectories in B/F. Frozen-H replay is
the follow-up if a discrepancy needs localization, not an extra execution
authorized by this one-run request. Keep full diagnostics on in every arm.

Use the exact frozen P3i cell/input/UPF/start and strict-reference manifests
(`/data/sig-dft/p3i/REPORT.md`, and its referenced `/data/sig-dft/vc/` and
`/data/sig-dft/p3g/` records). Resolve and record actual paths and SHA-256 on
dev08; case labels are not input specifications. Preserve coordinates,
cutoff, FFT grid, pseudopotential bytes, bands, occupied count, atomic start,
threshold, independent verification and output policy. Keep adopted
WIDE_ROTATE, SINGLE_GEMM, RESTART_CORRECTIONS, SWAP_BUFFERS and LEAN_RESET on;
keep FFT_FUSED_ROWS and REFSCF_VERIFY_ENERGY off. Retain all other effective
settings exactly and save the complete switch manifest.

Build one reviewed candidate with the established MPI/MKL environment/native
flags and no installs:

```sh
CARGO_TARGET_DIR=target/p3-stream cargo build --release -p dft-cli --features mpi,mkl
```

Record full source commit/patch hash, compiler/MPI/MKL versions, build flags
and exit code, binary SHA-256, input/UPF/start hashes, exact command lines,
node load, and actual binding to socket-0 physical cores 0–7. Serialize
builds, master gates and all measurement jobs. Invocation skeleton (master
resolves the frozen reference-driver input and unique output destination):

```sh
BAANDI_REF_DAV_SQRT=off BAANDI_REF_DAV_FUSED_APPEND=off \
BAANDI_REF_DAV_CHECK_APPEND=off BAANDI_REF_DAV_REUSE_ROWS=off \
BAANDI_DAVIDSON_STAGES=on OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 8 --map-by core --bind-to core --report-bindings \
  <reviewed-binary> run <frozen-input-with-unique-output-path>
```

Change only the four candidate flags and output destination per arm. Capture
the existing caller-resolved CPU profile during that same execution: hypot,
sqrt/vector sqrt, copies/zeros, norm loops, finite scans and collective waits.
Record tool/rate/sample counts. Fewer than 30 relevant samples cannot support
a quantitative attribution. If profiling fails, retain the single execution
and report the limitation without rerunning. Inspect optimized x86 code for
packed sqrt/divide and absence of per-row calls in the ordinary B kernel;
Mac code generation is not evidence of dev08 compiler/backend behavior.

## Required evidence and per-stage work table

Retain stdout/stderr, exit status, external delivered wall, results.json,
density/eigenvector fingerprints and tight-reference provenance. Compare all
band flags/active sequences, Ritz values and threshold operands, dropped
corrections, repairs, restarts, H calls/columns, GEMM shapes, maps, verification
attempts and outcomes. C/D/E must be bit-identical to A. B/F must retain band
decisions and H columns in this screen: list the first changed threshold
operand if any, even when delivered accuracy still passes. Report absolute
energy error in Ha against the independent tight reference (strict 1e-10 Ha,
standard 5e-9 Ha) and fixed-point certificate. Never mask a changed branch
with a final rounded energy or adjust a tolerance to pass.

For **every** stage in `convergence.work.davidson_stages.work`, including
setup, precondition, normalize, append, apply_check and repair_panel, emit:

| case/arm/rank/stage | solves/calls | exclusive s | GEMM/scalar FLOPs | sqrt/hypot calls | copy/zero bytes | logical streamed/GEMM bytes | collective calls/bytes/s | nominal GFLOP/s | logical GB/s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| filled from raw records | | | | | | | | | |

Keep rank0, all per-rank records and componentwise maxima; compute rates on
each rank before maxima. Do not divide max FLOPs by another rank's max time
or add maxima as if they came from one rank. Collective seconds are included
in stage seconds. Logical traffic is not measured DRAM bandwidth. State
exact cold/reused setup counts: E reuses four row matrices after successful
same-shape solves, while two returned row matrices and adapter scratch remain
fresh. Initial H validation and exceptional repair normalization remain.

Expected identities: B ordinary precondition hypot calls become zero and
sqrt calls equal its baseline hypot count; a nonzero B hypot count exposes
an extreme-column fallback. C saves 8*r*k logical bytes at ordinary accepted
normalization/appends; D saves 8*r*k at each correction H check/append; E
saves 16*r*(capacity+bands) zero bytes per same-shape reused solve. Count
moves between normalize/apply_check/append stages must not be mistaken for
removed arithmetic. Report the sum across affected stages as well.

## Local attempt record

Implementation and tests are recorded here because the ledger, optimization
log, dashboard and generated README evidence lie outside ticket ownership.
No new reference dataset or performance acceptance is claimed. The master
can transfer this record to the central ledger after review.

Baseline source: `2717b5bbdee8559c3d5ab1b139a83841c691310a` plus this
uncommitted patch. Mac compiler: `rustc 1.99.1 (58cf8f9c3 2026-10-13)`.

| Command | Exit | Result |
| --- | --- | --- |
| `cargo check --workspace --tests` | 0 | workspace test compilation |
| `cargo test -p dft-solvers` (initial) | 101 | new traffic oracle incorrectly equated normalized columns with H columns; corrected as described below |
| `cargo test -p dft-solvers` (final) | 0 | 176 passed, 25 ignored, no failures |
| `cargo test -p dft-engine --test refscf` | 0 | 24 passed, 1 ignored |
| `cargo test -p dft-numerics switches` | 0 | 14 passed |
| `cargo fmt --check` | 0 | clean |
| `cargo clippy -p dft-solvers -p dft-numerics --all-targets -- -D warnings` | 0 | clean |
| `bash scripts/check.sh dft-solvers` | 0 | build/test/fmt/clippy PASS; 17 test binaries, 176 passed, 25 ignored |
| `cargo test -p dft-solvers --lib stream -- --nocapture` | 0 | 5 focused tests, including fallback and failure agreement |

Gates used worktree-local `BAANDI_TEST_TMP`/`TMPDIR=target/tmp` and
`CHECK_SCCACHE=0`. The check script's optional `ps` load probe was denied by
the sandbox; the actual four gates all completed and passed. The full
workspace check was rerun after the test-only traffic-oracle correction and
also exited 0. The final solver test, fmt and clippy gates include that
correction. No production code changed afterward.

The first traffic assertion expected 192 bytes saved; the counters correctly
reported 384 bytes. Three directions had been normalized/appended but only
one survived factor repair to receive H. The corrected independent formula
is `8 * (normalized_row_entries + H_correction_row_entries)`, using the
normalization division count and H column count separately. No numerical
tolerance or solver decision changed to fix the test.

The main equivalence grid comprises dimensions/bands/capacity
(12,3,6), (17,5,15), (7,4,7), cyclic rank partitions at 1/3/8 ranks (including
an empty rank), clustered spectra in the latter two fixtures, restart reuse
off/on, and arms B/C/D/E/C+D+E/F. Each arm is executed twice, exercising
retained storage: **108 configurations, 216 comparisons** against 18
baselines. Every active-band sequence, decision, repair/drop count, GEMM
shape, communication count and H count is unchanged. C/D/E/C+D+E preserve
eigenvalues, X, HX and complete Ritz histories bit for bit.

| Synthetic dimension/bands/capacity | H columns | Correction iterations | Maximum eigenvalue difference, Ha |
| --- | --- | --- | --- |
| 12 / 3 / 6 | 13 | 4 | 2.220446049250313e-16 |
| 17 / 5 / 15 | 20 | 3 | 0 |
| 7 / 4 / 7 | 7 | 2 | 0 |

These are Hamiltonian eigenvalue differences, not SCF total-energy accuracy
or real-cell validation. The rounding arms changed no branch in this grid.
The denominator sweep covers 87,537 values (dense ordinary range, signed
binary exponents, deterministic mantissa samples): 31 output denominators
differ, maximum relative difference **2.195943771839758e-16**, below the
8*EPSILON comparison bound. Chunk tails at 0/1/31/32/33/65/257 rows retain
ascending norm bits. A finite extreme diagonal uses exactly 36 hypot calls
and zero precondition sqrt calls; the ordinary fixture uses zero hypot and
36 sqrt calls. Exact-zero compaction, failed H application and a NaN from
only one simulated rank pass, with a successful solve after each failure.
The cache test verifies exactly `16*r*(capacity+bands)` fewer zero bytes on
second same-shape solves.

Code-generation observation only: an opt-level=3 compile of the exact
`column::<true>`/denominator bodies, with a public wrapper to instantiate
them, emits `fsqrt.2d`, `fdiv.2d`, and `fmaxnm.2d` in the ordinary row loop on
Apple Silicon, without a hypot call in that loop. This is neither a timing
nor an x86/MKL result. Reproduce from the checked-in bodies:

```sh
python3 - <<'PYCODE'
from pathlib import Path
s = Path('crates/dft-solvers/src/ref_davidson/stream.rs').read_text()
s = s[s.index('const SQRT_LIMIT'):s.index('pub(super) fn count')]
s = s.replace('pub(super)', 'pub')
b = Path('crates/dft-solvers/src/ref_davidson.rs').read_text()
b = b[b.index('pub fn preconditioner_denominator'):b.index('\n#[cfg(test)]')]
s += b + '\npub fn kernel(h: &[f64], x: &[f64], hx: &[f64], out: &mut [f64], e: f64) -> f64 { column::<true>(h,x,hx,out,e) }\n'
Path('target/tmp/stream_codegen.rs').write_text(s)
PYCODE
rustc --crate-type lib --edition 2024 -C opt-level=3 --emit asm \
  target/tmp/stream_codegen.rs -o target/tmp/stream_codegen.s
```

The source-generation and rustc commands both exited 0. Local raw logs and assembly remain in
`target/tmp/p3-stream-*` and `target/tmp/stream_codegen.s`; they are disposable
build artifacts, not committed reference datasets. dev08 timings, actual MPI,
x86 vectorization, MKL equivalence, full real-cell convergence and promotion
remain unmeasured. No QE comparison or performance improvement is claimed.
