# P3-SETUP4: shared FFT plans, atomic Gram and owned local tables

Base: `07671b7afa9155aa6d3f7c99d9d061b86c65f874`, uncommitted task diff.
All three new switches default **off**. No dev08 run, speed claim, default
promotion, new QE reference, dependency or QE source access in this task.
Derived from the initialization principles in the own-words note
`docs/research/qe-study/qe-scf-cycle.md`; implementations use BAANDI's APIs.

## Diagnosis and implementation record

The ticket relays si128 MPI8 reference setup 2.30 s at e5ec8ae, verify off,
versus QE init_run 0.78 s (nested wfcinit 0.68 s). These are motivation, not
new measurements. P3l's nested H/Gram clocks must not be added to their parent.

- `REFSTART_FFT_SHARE`: a constructor scope owns a cache for one run/rank.
  rustfft plans are keyed by length/direction through one FftPlanner. MKL
  descriptors are keyed by length AND batch count; zero batches allocate no
  descriptor. Different sparse/dense supports and first-axis choices may share
  identical 1D kernels, never scatter maps or distributions. The start and H
  constructors participate in the same scope. DensityPlan gets a fresh FFT
  workspace from StartFields, copying its already computed layout and cloning
  kernel handles without re-planning or repartitioning. All coefficient,
  transpose and scratch buffers are private. The scope is restored on errors;
  kernel handles outlive it through the returned state. No cross-run cache.
- Backend audit: `RustLanes::plan` previously made a fresh rustfft planner per
  consumer (caching only inside that constructor). `MklLanes::plan` constructs
  three independent descriptors and the existing DFTI adapter creates/commits
  them explicitly; there is no application descriptor cache. The new solver
  API uses that public adapter with one nonempty axis per cached descriptor.
  MKL descriptor access is mutex protected because its API requires mutable
  access; the reference consumers execute sequentially. Measure this lock
  overhead in SCF as well as setup. Rust plans use immutable Arc handles and
  private scratch, with no execution lock. FFT_FUSED_ROWS shares its rustfft
  plans too; it remains an unrelated default-off switch.
- `REFSTART_SINGLE_GEMM`: S=X^t X and every H border X^t HX use
  `multiply_slices_single_into` when REFSTART_MKL is on. Previously these
  used the 256-row reduction-panel entry. Each border retains its dimensions,
  B-wide H chunking, rank reduction and symmetrization. On MKL this changes
  reduction rounding; pure Rust/MKL_DGEMM=off use the existing pure arithmetic.
  Rotations already use `multiply_slices_wide_into` and are unchanged.
- Atomic H audit: `projected_h` already calls the SAME
  `DistributedRealHamiltonian::apply_block` as the SCF adapter, with the same
  operator object. It already uses paired batched FFTs, M2_DGEMM,
  NONLOCAL_SINGLE_GEMM, NONLOCAL_FUSE_COUPLING, NONLOCAL_BACK_ACCUM and the
  wide GEMM fallback. No physics kernel edit, extra H route or reduced atomic
  column count is justified. All J columns still receive H, in B-wide chunks.
- `REFSTART_LOCAL_SHELLS`: both arms already compute radial/phase values only
  for used owned G. The candidate computes exact-bit shell IDs and sqrt(G^2)
  once for all species, then indexes contiguous per-species radial tables.
  Direct angles are formed in a contiguous independent-row pass; phase and
  charge/factor accumulation preserve each row's species/atom order. No
  interpolation, factored phases, quadrature or tolerance changes. This exposes
  independent arithmetic to vectorization; no measured SIMD/speed claim.
  Species-specific radial integrals remain distinct. Extra transient memory is
  O(owned density rows + owned shells), not grid-by-atom storage.

Density basis construction, basis sorting, metadata copies, buffer allocation,
EVD/factors and the full-indexed operator potential remain. A density-plan
clock is not zero after sharing: it still includes basis and metadata work.

## ONE profiled run per arm, four cases

The master relays to `baandi-dev08`; no implementer remote access. Reuse the
frozen P3s case inputs for **si64, si128, elongated, rskew, MPI8** with the
same matched QE MPI8 provenance. One worker/OpenMP/MKL thread per rank,
physical socket-0 cores 0-7, identical binding and all physical inputs.
Record exact input/UPF hashes, geometry, cutoff, FFT grid, bands, atomic start,
threshold, source/diff and binary hashes. Freeze all adopted switches at the
reviewed base, especially REFSTART_MKL, REUSE, REFSTART_GENERALIZED,
REFSTART_COMPACT_G and the nonlocal switches. Check effective manifests;
input-table overrides must not silently mask arm settings.

| Arm | REFSTART_FFT_SHARE | REFSTART_SINGLE_GEMM | REFSTART_LOCAL_SHELLS |
| --- | --- | --- | --- |
| D | off | off | off |
| F | on | off | off |
| G | off | on | off |
| L | off | off | on |
| FGL | on | on | on |

Exactly **20 BAANDI invocations**, one profiled full-convergence run per arm
per case, each a <=600 s fast job. No warmup, repetitions, A/A or precision
campaign. Preserve failures/timeouts and elapsed cost; no automatic retries.
Use the P3s **verify=off** inputs in all arms to isolate the reported start
cost and QE-comparable SCF interval. Label every row verify-off; these runs do
not establish independent verification or production delivered-wall acceptance.
Keep switches off pending dev08 review and the separate verification/accuracy
gates; this request does not authorize a default flip.

Build once with the established toolchain and MKL/MPI installation:

```sh
cargo build --release -p dft-cli --features mpi,mkl
```

One job invocation (master fills frozen absolute paths and arm values):

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BAANDI_REFERENCE_SCF=on BAANDI_REFSTART_MKL=on \
BAANDI_REFSTART_FFT_SHARE="$fft_share" \
BAANDI_REFSTART_SINGLE_GEMM="$single_gram" \
BAANDI_REFSTART_LOCAL_SHELLS="$local_shells" \
BAANDI_DAVIDSON_STAGES=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
  /usr/bin/time -p taskset -c 0-7 \
  mpirun -np 8 --map-by core --bind-to core --report-bindings \
  "$baandi_binary" run "$baandi_input" --output-dir "$unique_output"
```

Record MKL_DGEMM/M2_DGEMM/M2_FFT effective settings, build flags and compiler,
MPI/MKL versions. Retain raw logs, results/density, start/end, exit, bindings,
node load, RSS and external wall through required output and MPI exit.

### Required adjacent setup table

This is clock attribution, not a speed ratio: profiled BAANDI setup and QE
init_run have different boundaries. Retrieve existing matched QE **MPI8** raw
init_run/wfcinit with command/version/input/UPF provenance and fastest eligible
diagonalization layout. Missing provenance means unavailable, never substitute
the ticket's numbers. No QE run is requested. No cross-parallel-model column.

| Case, MPI8 | D setup_s (s) | F (s) | G (s) | L (s) | FGL (s) | QE init_run (s) | QE wfcinit nested (s) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| si64 | pending | pending | pending | pending | pending | retrieve | retrieve |
| si128 | pending | pending | pending | pending | pending | retrieve | retrieve |
| elongated | pending | pending | pending | pending | pending | retrieve | retrieve |
| rskew | pending | pending | pending | pending | pending | retrieve | retrieve |

Next to EACH case row include the complete setup_s split for every arm:
engine and CLI-inclusive setup, main_to_first_map, field_fft_plan,
operator_fft_plans_and_tables, density_fft_plan, basis_gvectors,
atomic_local_tables_structure_factors, atomic_overlap_gram, factor/metric,
atomic_h_apply, atomic_h_gram, atomic_h_apply_gram (inclusive), raw fill,
QR/fallback, generalized EVD/rotation and atomic_block_total (inclusive).
Report rank 0 and min/max/argmax where captured. Do not add nested phases,
maxima from different ranks, or wfcinit to init_run. Also retain per-map H/FFT
and total SCF costs to detect descriptor-lock regressions.

Report actual G, each L_r, J, B, grid, projector count and chunk widths. S
costs 2 L_r J^2 FLOPs; all H borders together cost 2 L_r J^2. The change is
one MKL call per rectangle, not fewer mathematical FLOPs or H columns.
For si128's nominal J=512/B=260, confirm two H chunks (260,252), one S Gram
and two H Grams, with the same rank-reduction payloads as D.

F/L must preserve every start/trajectory/delivered scientific bit. G/FGL may
change rounding: preserve start eigenvalues within 1e-12 Ha, orthonormality
and linear selected-subspace residual within 1e-12, and final unique-solution
energy within 1e-12 Ha as an internal regression. Retain strict 1e-10 Ha /
standard 5e-9 Ha independent-reference classes for later adoption. Report
first-map dr2 signed deltas (Ry), every map/redo, Davidson iteration/H column,
restart/repair/drop, reduced dimension and convergence status. CLI does not
emit raw start orbitals; the new start-oracle test on the reviewed MKL build
checks those properties. Do not infer start identity from delivered density.

## Mac checks and integration limits

Local checks use pure Rust and simulated collectives, never real MPI/MKL or
QE. Final gate results are recorded below. All changes remain
uncommitted. Central attempts/log/dashboard/session files, README generator
and reference evidence are outside OWNED_FILES; this document is the master's
integration record. The local sessions registry has no P3-SETUP4 row and no
`.delegate/P3-SETUP4/report.json`; no nested delegate was launched. No dashboard
upload or recurring orchestrator timer was started by this implementer.


### Executed local gates

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | All workspace test targets checked |
| `cargo test -p dft-engine --test refscf` | 0 | 36 passed, 5 existing ignored |
| `cargo test -p dft-solvers` | 0 | 185 passed, 25 existing ignored across targets |
| `cargo test -p dft-numerics switches` | 0 | 8 unit + 6 integration passed |
| `cargo fmt --check` | 0 | Format gate |
| `cargo clippy -p dft-engine -p dft-solvers -p dft-physics -p dft-numerics --all-targets -- -D warnings` | 0 | No warnings |
| `cargo test -p dft-engine --test refscf setup4 -- --nocapture` | 0 | 3 new tests, 0.22 s combined (debug-budget observation only) |
| `BAANDI_TEST_TMP="$PWD/target/tmp" CHECK_SCCACHE=0 scripts/check.sh dft-numerics` | 0 | build/test/fmt/clippy PASS; 90 passed, 22 ignored |
| `git diff --check` | 0 | No whitespace errors |

The check.sh load diagnostic could not run sandbox-restricted `ps`; this was
not a gate failure, and all four actual gates completed. Early development
checks returned 101 for an incorrect FftPlanner Default derive and then an
incorrect unwrap of simulate's Vec return; both were fixed before these gates.
No failing test tolerance was changed. Existing tests were untouched except
the two authorized switch counts (165 to 168); new tests were appended.

- FFT cache/workspace test: full and sparse supports on 4x4x4 and 3x5x7
  grids, ranks 1/8, both vendor selections (pure fallback on this Mac),
  sparse-clear and fused rows enabled. Repeated constructors and independent
  workspace copies request no new kernels. Unique rustfft length/direction
  counts are **2** and **6**, respectively, and a new cache starts at zero.
  Alternating applications on separate buffers match private-plan bits.
- Start test: 12 shape/rank configurations, 108 logical starts including
  **96 candidate comparisons** across F/G/L/FGL at 1/2 workers. Cubic,
  skew two-potential, elongated and slab; ranks 1/3/8; nonlocal/local-only,
  J>B and J<B, odd column tails; the cubic serial fixture crosses 256 rows.
  Every field, Ritz value and orbital bit matches the baseline on pure Rust;
  measured Ritz difference is **0 Ha**. Orthogonality and linear subspace
  residual satisfy the unchanged **1e-12 dimensionless** assertions (bounds,
  not measured maxima). Results match between worker counts at fixed ranks.
- Four verified trajectory pairs (eight complete SCFs): cubic and elongated
  at ranks 1/3. F+L versus D gives **delta E = 0 Ha**, identical delivered
  density bits, every dr2 bit and per-map H work/call count, with passed
  independent verification. Map counts are **21/21** for cubic at ranks 1/3
  and **9/9** for elongated. No across-rank bit-identity claim.
- Grepped the test tree for all three switch names; the added test explicitly
  asserts each default is off. No default flip is proposed. MKL single-call
  rounding, descriptor behavior/performance and real MPI remain dev08 gates.

Raw local logs for this session: `/tmp/setup4-{check,refscf,solvers,switches,fmt,clippy,focused,gate}.log`.
These are local scratch, not committed benchmark evidence; the summary and
reproduction commands above are the durable handoff. No numbers are labelled
as new QE reference or Mac performance acceptance.
