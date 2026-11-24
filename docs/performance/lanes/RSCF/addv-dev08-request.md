# P3-ADDV: nonlocal coupling and direct back-application

Status: B and C are implemented, default off. The orchestrator approved the
beta=1-only adapter/test ownership extension. Final Mac validation is recorded
below; native MKL execution and dev08 performance remain pending.

Source base: `941a67fd68f9a5619a8476ea9206b747dfbaf887`. Changes are
uncommitted. Derived from the matrix-product structure in the own-words note
`docs/research/qe-study/qe-hpsi-and-fft.md`; no QE source read. Defaults stay
unchanged. Mac checks establish correctness, not speed. This request is for
ONE profiled run per arm and cell, not a promotion or precision campaign.

## Per-H-call accounting

Use r = owned ORC rows on the particular rank, p = total projector columns,
b = this H block's active width, C = allocated band capacity, s = max(p,C).
Do not infer these values from a material name. P is column-major r by p;
X and HX are r by b; Q and W are p by b. All are already contiguous with
leading dimensions r or p. Atom a has c_a coupling channels and a row-major
c_a by c_a D block. Mixed species and off-diagonal D are supported.

| Stage | Logical shape (m,n,k) / operations | Baseline route and storage |
| --- | --- | --- |
| Local Q = P^T X | (p,b,r); 2rpb FLOPs | `other/mkl_single` with NONLOCAL_SINGLE_GEMM and M2_DGEMM on; one sequential MKL call into existing p*b tile, then copy tile to local projection |
| Global Q | deterministic tree sum of p*b doubles | reads local projection, writes `reduced`; 8pb bytes payload, unchanged |
| W = DQ | 2b sum_a(c_a^2) nominal FLOPs | scalar atom/band/j/k loop including zero D entries; writes separate W directly, no W copy |
| Image = P W | (r,b,p); 2rpb FLOPs | wide entry, ceil(r/256) output-row panels; one sequential DGEMM per panel, beta=0, copy each tile column into image |
| HX += image | rb additions | separate contiguous pass, 24rb nominal bytes (read HX/image, write HX) |

Back route: `other/mkl` if r*b <= s*s, `other/mkl_wide` otherwise,
provided min(256,r)*b <= s*s. The latter panel bound can reject small-short
shapes (`pure_panel_bound`); the long/short/layout bounds also apply. These
are **logical wrapper** counts: a single `mkl_wide` row call may contain many
native DGEMM calls. With M2_DGEMM off the caller uses `unlabelled/pure_direct`;
without MKL the slice entry reports `pure_build`; MKL_DGEMM off reports
`pure_switch`. Empty owned ranks skip local products but join the tree sum.

There is **no per-H Rust zero fill of W or image** in the ordinary nonempty
path. Construction zero-initializes storage once. DGEMM beta=0 replaces its
tile; the pure product replaces its first reduction chunk and adds subsequent
chunks. HX is first overwritten by local+kinetic unpacking. Empty ranks zero
only their local p*b projection before the collective. Zero-width H returns
early. W is borrowed by every back panel without a Rust copy; BLAS-internal
packing is not an observed Rust copy. Removing an imagined W-copy/zero-fill
cannot be counted as a saving.

### Why 103 / 206 is not half of the projections falling back

`crates/dft-cli/src/run.rs` builds the P3-HAPPLY display using
`single = s.wide_calls / 2` when NONLOCAL_SINGLE_GEMM is on. The solver's old
counter records TWO entries per H: one projection and one back product. This display produces
103 single + 103 wide out of 206 entries for 103 projection entries and
103 **back** entries. It is not a measured projection shape fallback.
`DistributedRealHamiltonian::projections_reduced` selects the single entry
for every nonempty rank with M2_DGEMM on, including the nonlocal-energy audit.
The p*b output always fits max(p,C)^2 because b <= C. No remaining nonempty
projection shape in this path needs another switch. Serial `RealNonlocal`
uses its own MKL_NONLOCAL path and is outside this MPI8 mechanism.

This algebraic explanation does not reconcile different profiling windows:
the ticket also quotes 93 H calls. Raw per-stage/per-rank records are needed
to reconcile 93 with 103 (startup, audit, verification and solve boundaries).
The 31b23e11 raw run is not in this checkout; do not invent numerical shapes,
widths, exclusive coupling clocks or achieved GFLOP/s. Existing CLI derived
entry labels are estimates, not the authoritative backend trace; with the
new direct-back arm they will still describe the old entry arrangement.
Use actual `gemm_dispatch` routes/shapes and stage-delimited records.

The ticket's motivating rank-0 BAANDI rates are 0.0193 s/calbec call and
0.0246 s/add_vuspsi call (93 calls, one relayed run). The rounded products
are 1.7949 s and 2.2878 s, **not recovered raw totals**. QE rates supplied
in the ticket are motivation only: no new QE reference or speed comparison
is established without the original command/version/input/UPF hashes.

## Switches and numerical contract

C: `BAANDI_NONLOCAL_FUSE_COUPLING=on` (default off). After the deterministic
rank reduction the local projection buffer is dead. Write DQ there directly
and pass it as W to the back product; omit the separate p*C-double weights
allocation. The global reduced Q remains untouched. The full D loop and
its order are unchanged, so this is bit-identical even for off-diagonal D.
There is no extra W copy, no in-place overwrite hazard and no allocation in
H. This is storage/lifetime fusion, not moving D through the rank reduction
(which would change rounding). It saves 8pC bytes per operator; any cache
or wall benefit is unmeasured and may be zero.

B: `BAANDI_NONLOCAL_BACK_ACCUM=on` (default off), conditional on
M2_DGEMM and MKL_DGEMM and an MKL build. One sequential DGEMM with
(m,n,k)=(r,b,p), lda=r, ldb=p, ldc=r, alpha=1, beta=1 writes directly into HX.
No whole-output scratch bound, image tile, image copy or separate HX-add
pass. Keep the ordinary panel route as the pure/MKL-disabled fallback.
This is rounding class: BLAS grouping and beta accumulation can differ even
when the reduction length p is unchanged. No material-specific dimension
threshold or changed rank reduction. Dispatch reports `other/mkl_single`;
shape distinguishes it from projection. The existing per-call span must
include the direct GEMM and output finite check; no standalone Add span.

The new `dft-mkl::dgemm_add_at` shares checked dispatch and native timers
with the unchanged public beta=0 `dgemm_at` API. Zero reduction preserves C
for accumulation; replacement still clears it. Both validate LP64 dimensions,
C length and (for nonempty products) input strides/addressed lengths before
FFI. The unsafe block documents initialization, exclusive output borrowing,
layout, LP64 bounds and synchronous pointer lifetime. No duplicate CBLAS
binding is introduced outside the established adapter (ADR 0019).

`RealGemmWorkspace::try_multiply_slices_single_add` validates shapes in every
build. On pure builds or MKL_DGEMM=off it returns false without touching HX or
recording a product; the original fallback then records its effective route.
Thus Mac as-if-MKL counters describe the executed panel fallback, not the
unexecuted beta=1 candidate. The image allocation remains available for the
fallback but is not read or written on the direct route. Existing phase byte
estimates remain logical-product estimates, not measured direct-HX traffic.
The direct Back span includes accumulation; there is no separate Add span.

Contract tests cover nonzero initial C, repeated accumulation, padded inputs,
transpose, 257-row products, zero dimensions/reduction (including signed zero
and NaN payload preservation), invalid lengths/strides and LP64 overflow.
The numerics test covers output larger than the old scratch bound, pool
budgets 1/4, unchanged fallback output/storage and native nonfinite rejection.
The H fixture additionally poisons the image buffer: on native MKL with B on,
it must remain untouched. No MKL numerical execution is claimed on the Mac.

## dev08 request (master / baandi-dev08 only)

Four frozen campaign cells: si64, si128, elongated, rskew. MPI8 only, one
thread per rank, physical cores 0-7 on socket 0 with recorded actual binding.
Reuse exact input/UPF/cutoff/grid/band/start/threshold/verification settings.
Record source commit + patch SHA-256, binary/input/UPF SHA-256, compiler and
MPI/MKL versions, full command, exit status and load. No installs or remote
execution by the implementer.

| Arm | NONLOCAL_BACK_ACCUM | NONLOCAL_FUSE_COUPLING |
| --- | --- | --- |
| A | off | off |
| C | off | on |
| B | on | off |
| BC | on | on |

Exactly **16 profiled runs**, ONE per arm per cell, no warmups or repeats.
Run the native correctness prerequisites below before profiling. Hold
NONLOCAL_SINGLE_GEMM=on, M2_DGEMM=on, MKL_DGEMM=on and all unrelated switches
fixed. Use fixed-H replay for attribution where the frozen campaign provides
it; use its full reference-SCF profile to recover calbec/add_vuspsi clocks.
Do not mix replay clocks and full-run clocks. One case/arm per 600 s fast job;
a longer job requires the master's named authorized final job (up to 7200 s).
Verification remains at the identical campaign setting and its work is
separated; no setting is weakened to obtain a performance result.

```sh
CARGO_TARGET_DIR=target/p3-addv cargo build --release -p dft-cli \
  --features mpi,mkl,dft-numerics/dgemm-timers

BAANDI_NONLOCAL_BACK_ACCUM=off BAANDI_NONLOCAL_FUSE_COUPLING=off \
BAANDI_NONLOCAL_SINGLE_GEMM=on BAANDI_M2_DGEMM=on BAANDI_MKL_DGEMM=on \
BAANDI_DAVIDSON_STAGES=on BAANDI_STAGE_TIMERS=on BAANDI_PHASE_PROFILE=on \
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/usr/bin/time -p mpirun -np 8 --map-by core --bind-to core --report-bindings \
/absolute/reviewed/dft-cli run /absolute/frozen/input.toml \
--output-dir /absolute/unique/p3-addv/case/arm
```

Resolve exact paths and ensure socket-0 core binding with the master before
launch. Keep complete rank outputs and stdout/stderr. Detailed profiling is
for BAANDI A/B hypothesis tests only, never a QE-comparison table.

Return per case/arm/rank: r,p,C, exact columns_per_call, H calls/columns;
projection/reduction/coupling/back-and-add seconds; calbec and add_vuspsi
seconds/call; nonlocal `gemm_dispatch` rows (site,route,calls,FLOPs,max_m/n/k)
and native DGEMM calls/time. Attribute startup and energy audit separately;
max shapes alone are not a shape histogram. For a rank and one common window,
F = sum_H(2rpb), projection GFLOP/s = F/(1e9*projection_s), back effective
GFLOP/s = F/(1e9*back_and_add_s). Coupling work is sum_H(2b sum(c_a^2)).
calbec/call = (projection_s+reduction_s)/H_calls;
add_vuspsi/call = (coupling_s+back_and_add_s)/H_calls. Never divide rank-0
work by a max-over-ranks clock. Byte counters are nominal logical traffic,
not measured memory bandwidth; direct beta=1 also reads existing HX.

C must retain H/energy/density bits and all work counts. B/BC must pass the
H component relative/absolute 1e-13 gate and report energy differences in Ha,
density differences, SCF work and verification outcomes, with strict 1e-10 Ha
and standard 5e-9 Ha classes against the independent oracle. Native MKL
thread/rank determinism and mixed-projector tests remain prerequisites for
promotion. A single profiled run has no uncertainty estimate; no default
promotion or claimed speedup follows from this request.

## Native correctness prerequisites (dev08, not extra profiled arms)

```sh
cargo test -p dft-mkl --features mkl --test dgemm_accumulate
cargo test -p dft-numerics --features mkl single_add_preserves_fallback_and_accumulates_beyond_tile_bounds
cargo test -p dft-physics --features dft-numerics/mkl --lib nonlocal_single_projection_matches_panel_route -- --nocapture
```

These must pass before the 16 requested profiles. The synthetic H fixture
uses pool budgets 1/3, single projection on/off and all four B/C combinations;
it asserts exact coupling/storage equivalence and rounding-class back-product
equality, including zero-owned-row ranks. Also retain the campaign's native
MPI rank reproducibility and independent SCF verification requirements.

## Local validation

All requested Mac commands passed. Logs: `target/tmp/addv-final-*.log`;
command/exit-code record: `target/tmp/addv-final-commands.json` (local,
untracked evidence). Builds used CARGO_BUILD_JOBS=3; the check script also
used CHECK_TEST_JOBS=2, CHECK_SCCACHE=0 and BAANDI_TEST_TMP/TMPDIR under this
worktree's target/tmp. No Mac timing is speed evidence.

| Command | Exit |
| --- | ---: |
| `cargo check --workspace --tests` | 0 |
| `cargo test -p dft-physics` | 0 |
| `cargo test -p dft-engine --test refscf` | 0 |
| `cargo test -p dft-numerics` | 0 |
| `cargo fmt --check` | 0 |
| `cargo clippy -p dft-physics -p dft-numerics --all-targets -- -D warnings` | 0 |
| `cargo check -p dft-mkl --features mkl` | 0 |
| `cargo check -p dft-mkl --features mkl --tests` | 0 |
| `cargo check -p dft-physics --features dft-numerics/mkl --tests` | 0 |
| `cargo test -p dft-physics --lib nonlocal_single_projection_matches_panel_route -- --nocapture` | 0 |
| `cargo clippy -p dft-mkl --features mkl --all-targets -- -D warnings` | 0 |
| `BAANDI_NONLOCAL_BACK_ACCUM=on BAANDI_NONLOCAL_FUSE_COUPLING=on cargo test -p dft-engine --test refscf` | 0 |
| `scripts/check.sh fast dft-physics dft-numerics dft-mkl` | 0 |

- Physics: 97 passed, 18 ignored. Numerics: 90 passed, 22 ignored.
- refscf: 28 passed, 1 ignored both at defaults and with B/C enabled in the
  environment. These Mac runs exercise the pure fallback, not native beta=1.
- The final check script passed build/test/fmt/clippy: **188 passed, 0 failed,
  40 ignored**, 15 test binaries. Its advisory load probe could not execute
  `ps` in the sandbox; the four actual gates completed successfully.
- **350 rank/block fixtures**, maximum H component difference **0 Ha** on
  pure Rust. Five cells (including an origin-only cell with empty owned
  ranks), P=1/2/3/8, widths 0/1/2/5/17, mixed 0/1/3/5-channel atom blocks,
  off-diagonal D and pool budgets 1/3. All four B/C arms plus the panel
  projection baseline are covered. Coupling and repeated H after an energy
  audit agree bitwise; nonlocal band energies agree exactly.
- Registry count: 160; both new switches default off. Existing defaults and
  the beta=0 signature/semantics are retained.
- MKL adapter, adapter tests and physics MKL tests **type-check** on Mac;
  MKL-feature adapter clippy also passes. MKLROOT is unset, and Cargo reports
  that MKL executables cannot link. No native MKL numerical test ran here;
  the three dev08 prerequisite commands above remain unexecuted.
- `git diff --check`: exit 0. Changes remain uncommitted.

No speedup, QE agreement, native MPI run or default promotion is claimed.
Registry/dashboard/README/ledger changes remain outside this ticket's ownership;
this document supplies the orchestrator's artifact reconciliation handoff.
