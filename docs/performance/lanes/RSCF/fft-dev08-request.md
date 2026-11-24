# P3-FFT: reference H FFT movement candidates

Status: implementation ready for review, both switches default off. Mac
correctness only; no speed or adoption claim. Base
`d7ec4e9565a24d2f88c5b241e4fdd597f50d2e73`; changes left uncommitted.
The master owns dev08 execution through `baandi-dev08`, full gates,
ledger/log/backlog, dashboard/README integration and commit. Those shared
artifacts are outside this ticket's OWNED_FILES. No QE source was accessed;
`crates/dft-numerics/src/real.rs` remains untouched (P3-GRAM owner).

## Mechanisms and prior evidence

[Exact byte accounting](../../../design/reference-scf.md#h-apply-fft-cost-accounting)
records the ticket's `66d45e69` profile and all baseline stages. Existing
pruning, reuse and batched alltoall are retained. New mechanisms:

1. `FFT_DIRECT_PACK`: ORC-to-stick scatter avoids intermediate hbuf
   write/read, saving 32*L_r bytes per transformed pair on rank r. Arithmetic,
   returned layout, forward scaling and mirror exchange are identical.
2. `FFT_SPARSE_CLEAR`: cache contiguous complement ranges at construction;
   clear only destinations not overwritten by scatter. Saves
   16*(L_r + p_r*(S+Z*n1)) bytes per pair for unique G indices. Additional
   range traversal can lose to a full contiguous memset, so keep an
   independent arm. All three clear stages report their changed logical bytes.

No revivals of A11 lane staging, A28 many-band 3D, A34 FFT threads,
M2_FFT backend switching, async alltoall, kept-plane copy fusion, or stick
ownership balancing. PFFT-4's indexed-transpose failure and the neutral
FFT-KEPT-FUSE/closed arrival study are explicit counterevidence against
assuming traffic reduction guarantees gain. The new premise is the user's
reference-driver Si64/Si128 `66d45e69` profile and P3-FFT authorization; the
sparse candidate uses contiguous hole spans rather than an indexed
replacement transpose. These are unmeasured candidates, not adoption.

## Arms and identity

| Arm | BAANDI_FFT_DIRECT_PACK | BAANDI_FFT_SPARSE_CLEAR |
| --- | --- | --- |
| D | unset | unset |
| A | off | off |
| B1 | on | off |
| B2 | off | on |
| B12 | on | on |

Input keys under `[performance]` are `fft_direct_pack` and
`fft_sparse_clear`; normal environment precedence, consulted manifest and
cross-rank agreement hash apply. Set the same values on all ranks. D/A must
match bitwise; do not time D as a third optimization. Other effective
switches, including `M2_FFT`, `M2_STICKS`, verification and P3-GRAM controls,
remain fixed. Direct packing affects ordinary distributed H application,
batched or unbatched; density and kept-local-change still use their existing
packing. Sparse clear affects every use of this distributed FFT plan.

Use frozen Si64/Si128 inputs and UPFs from the original profile
(`/data/sig-dft/vc/` manifests), plus the established **elongated** and
**rskew** fixtures, all MPI8. Resolve their actual absolute paths from the
master's manifests; cell labels do not define grids, coordinates or bands.
Record input/UPF/start/snapshot SHA-256, full base/main/candidate commits,
source diff hash, compiler, MPI/MKL versions, build flags and binary SHA-256.
No reconstructed input, changed band count/cutoff/grid, tolerance, mixing,
start or output policy between arms. Charge independent verification,
failed attempts, recovery, initialization and output in delivered wall.
QE comparisons, if refreshed, are MPI8 against MPI8 under identical
conditions and require actual command/version/input/UPF provenance.

Build base, current main and candidate in separate target directories;
main drift is separate from candidate effect. Follow established native
build flags, identically in all arms. Example build (master only):

```sh
CARGO_TARGET_DIR=target/p3-fft cargo build --release -p dft-cli --features mpi,mkl
scripts/check.sh full
```

Use the existing MPI8 binding on socket 0 physical cores 0-7, one thread per
rank (`OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`), preserve `--report-bindings`
and load data. No builds overlap measurement. Fast jobs <=600 s including
output; split the campaign into bounded jobs. Full convergence that cannot
fit requires an authorized named `final-` job <=7200 s; never truncate it.

## Stage 0: untimed screen and work accounting

First run Si64, then Si128, elongated and rskew for A/B1/B2/B12 with
`BAANDI_DAVIDSON_STAGES=on`. These are diagnostic runs, excluded from the
acceptance sample. Retain per-rank JSON, stdout/stderr, exit codes and merged
rank min/max output. Start with identical fixed-H cold/warm snapshots where
the reference harness supports them, then the reference driver's complete
screen. A legacy solver replay is supplementary, never a substitute for the
reference-driver split. Use the existing invocation/input manifests from
[happly-dev08-request.md](happly-dev08-request.md); do not invent a new
CLI or change a production SCF budget to force completion.

Required rows per case/arm/rank, split ordinary/retry/verification:

- H calls, full columns-per-call sequence, complex pairs, odd tails;
  check columns = 2*pairs - tails and both directions' transform counts.
- apply_h total, outer pair_pack/unpack, mirror exchange, kinetic, local
  1D/2D transforms, both alltoalls, transpose pack/unpack, grid
  scatter/gather/scale and multiply. Preserve nonlocal rows as controls.
- Each FFT stage's calls, point-transforms, logical bytes; L_r, s_r, p_r,
  S, Z, grid, maximum stick/plane counts and batch widths from the actual
  layout or equivalent counters. Useful and padded payload, no double
  counting self blocks or incoming bytes. Do not add separate rank maxima
  to form a wall total.
- Output H/eigen/density bits, energy Ha, convergence, verification outcome,
  work counts, H-apply critical-rank closure and unassigned remainder.

Expected: B2 reduces only zero-store bytes in FFT stages 0/4/6, by the
formula; FFT point counts and MPI payload/calls stay equal. B1 removes the
outer input hbuf round trip, but **moves ORC arithmetic into FFT scatter**.
The probe's scatter bytes remain abstract coefficient traffic. Inspect
combined outer pair_pack + scatter and full apply, not pair_pack alone.
Claim neither an FFT flop reduction nor a smaller alltoall payload.
The combined arm should equal the sum of logical byte savings, while wall
benefits need not add. Density-only work responds only to B2.

Reject any correctness mismatch; do not loosen tolerances. If a switch has
a clear whole-apply regression, retain off and record it; a timing label
migration is not a regression or gain by itself. Only screen survivors
without a resolved slowdown advance. An unresolved noisy screen requires
A/A + crossed pairs rather than a speed verdict from one instrumented run.

## Stage 1: A/A floor, then crossed A/B per surviving switch

Disable diagnostics for acceptance timing: DAVIDSON_STAGES, FFT_PROFILE,
PHASE_PROFILE, PHASE_WAIT and COLLECTIVE_TRACE off; no new profiling or
barriers. Use external delivered wall through process/MPI exit. Hold all
other effective switches and resources constant.

1. In one exclusive window establish A/A noise on Si64 MPI8 (at least
   six pairs, equal warm-up policy). Preserve every raw sample and failure.
2. For each surviving B1 and B2 independently take at least six crossed
   pairs alternating A-B / B-A in the same window, interleaving controls.
   Report A and B medians using the same estimator, paired differences,
   SE, |mean|/SE, gain/floor, order effects and work/identity checks.
   Extend a small unresolved effect only under the live fast-loop policy;
   no arbitrary percent-gain cutoff. A resolved gain requires >=3x floor
   and |mean|/SE >=3; an unresolved sample is not an adoption claim.
3. Repeat survivors on Si128, elongated and rskew MPI8 with their own A/A
   floor. Measure B12 only after individual attribution; compare it against
   A and the best surviving individual switch to check interaction.
4. Full reference SCF plus independent verification, strict 1e-10 Ha and
   standard 5e-9 Ha accuracy, deterministic repeats and applicable rank
   matrix on the generality set precede promotion. No new default from a
   Si-only or fixed-H result. A resolved gain of 1-3% is eligible; a
   correctness failure or resolved regression is not.

Return exact commands/versions/hashes, raw tables, stderr/exit codes,
case/arm/rank work and split timing, confidence/noise diagnostics and a
separate verdict per switch. No mixed threads/MPI comparator columns.

## Mac validation

The local commands and actual results are recorded below after execution.
No Mac timings are acceptance evidence. Simulated ranks exercise placement
and collective semantics, not native MPI transport; pure builds do not
validate MKL kernels. Full gates and native MPI/MKL remain with the master.

Executed 2026-11-24, pure Rust Mac build, pinned toolchain:

| Command | Exit | Result |
| --- | --- | --- |
| `cargo check --workspace --tests` | 0 | Workspace and test targets checked |
| `cargo test -p dft-solvers` | 0 | 162 passed, 25 ignored |
| `cargo test -p dft-physics` | 0 | 95 passed, 18 ignored |
| `cargo test -p dft-engine --test refscf` | 0 | 19 passed, 1 ignored |
| `cargo fmt --check` | 0 | Formatting clean |
| `cargo clippy -p dft-solvers -p dft-physics -p dft-numerics -p dft-cli --all-targets -- -D warnings` | 0 | All touched crates clean |
| `cargo test -p dft-numerics --test all switches::` | 0 | 6 passed |

New transform matrix: 6 support/grid fixtures, ranks 1/2/3/8, two lane
backend requests (both rustfft on this pure build), five batch widths
including zero, three nonbaseline arms: **2520 rank-local comparisons**,
bit differences zero, sparse-clear logical-byte deltas exactly equal to
16*np*(L_r+p_r*(S+Z*n1)). Includes full and empty support, first-axis y,
odd grids, growing/shrinking batches and dirty buffer reuse. Existing
roundtrip/full-grid and serial pair-plan oracles also pass.

New H matrix: three actual cell geometries, ranks 1/2/3/8, greedy/mirror
ownership, batched/per-pair modes, five odd/even/empty widths and three
nonbaseline arms: **2520 rank-local comparisons**, H differences **0 Ha**
bitwise. This fixture isolates local-potential/kinetic packing (no nonlocal
projectors); existing distributed H/nonlocal/reference tests supply the
unchanged integration controls. The physics focused test additionally
printed this count in `target/p3-fft-focused-physics.log`.

Both new tests are fast (well under the debug 10 s limit); no benchmark or
performance conclusion is drawn from their runtime. Local gate transcripts
are `target/p3-fft-gate-{0..6}.log`, machine-readable exit records
`target/p3-fft-gates.json` (scratch, not committed). All 282 executed gate
tests pass; 44 existing tests remain ignored. `git diff --check` and the
owned-file audit pass; real.rs and pre-existing design sections are unchanged.
Only the two appended switches and count 140 -> 142 change the registry;
reconcile the additive count when merging P3-GRAM's separate appended entry.

No `scripts/check.sh full`, native MPI/MKL run, new QE comparison or dev08
speed measurement was performed locally. Default-off retention is
intentional. Shared evidence/dashboard/README generation and attempt-ledger
integration remain master-owned; this file records the candidate attempts,
results and limits within ticket ownership.
