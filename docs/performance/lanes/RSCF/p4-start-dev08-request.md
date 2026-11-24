# P4-START: atomic Ritz and FFT construction cost

Base: `1a3fc4347ce65e0933a6f03c189b8621ae0cebfa`, plus the reviewed uncommitted
P4-START diff. No default changes, remote measurements, new dependencies or QE
source access. Derived from the initialization principles in the own-words
note `docs/research/qe-study/qe-scf-cycle.md`, implemented with BAANDI APIs.

## Audit and attempt record

The ticket supplies a **relayed motivation snapshot**, fa4cfad, dev08 MPI8,
verify off: si128 engine setup 2.125 s and QE init_run 0.76 s; si216 5.834 s
and 2.90 s. These numbers have not been reproduced here and are not a new QE
reference. Retrieve the original commands/version/inputs/UPF hashes before
putting them in the diagnostic table below.

1. **Atomic H:** `projected_h` calls the same retained
   `DistributedRealHamiltonian::apply_block` object as the SCF adapter. It uses
   batched pairs (default true), the same FFT plans, M2_DGEMM,
   NONLOCAL_SINGLE_GEMM, NONLOCAL_FUSE_COUPLING and NONLOCAL_BACK_ACCUM.
   Every J atomic column is needed to preserve the Ritz space; retain B-wide
   chunks rather than a second J-wide operator. The supplied si216 1.764 s /
   1152 columns is 1.53125 ms/column, versus the supplied SCF approximation
   0.38 ms/column. This difference does **not** identify a slower code route:
   cold allocation/cache effects, chunk shapes, collective wait and the SCF
   kept-HX denominator need separation. Kept-HX cannot replace the initial
   fresh H: no earlier orbital density/slabs/HX exist. No new H kernel is added.
2. **Dense factor and EVD:** `REFSTART_CHOLESKY=on` replaces the hand-written
   complex Cholesky plus explicit triangular inverse with sequential faer
   real Cholesky and triangular solves. For D_ii=sqrt(S_ii), factor
   D^-1 S D^-1 = L L^t once; whiten A with two triangular solves, use one
   existing symmetric EVD (MKL dsyevd on the enabled MKL path), solve
   L^t Y=U_B, return X D^-1 Y. No inverse, no second factor and no second EVD
   on success. The unit-column 1e-6 pivot guard and 1e-12 metric-defect bound
   are unchanged; rejected factors/metrics take the existing QR/MGS fallback.
   J<B keeps the original completion path. The metric check also uses solves.
   Small work remains replicated once per rank, with **no MPI sum** of small
   matrices. Root-only work/broadcast is not introduced. Faer factor/solves
   are used on both builds because the existing MKL adapter has no dpotrf or
   dtrsm entry and that adapter is outside ownership. This is a rounding-class
   experiment, not a claim that faer beats MKL dsygvd.
3. **Tables/fill:** field coefficients already evaluate only represented owned
   density G, once per species, with exact-bit radial-shell caching and shared
   direct atom phases for charge/ionic sums. P3t LOCAL_SHELLS had no gain
   (0.415 -> 0.457 s on si128); leave it off rather than repeat that experiment.
   `REFSTART_SPECIES_FILL=on` addresses different duplicated work: build one
   amplitude table for a species, then evaluate direct phases once per atom
   and owned representative across *all* radial channels and m. Adjacent ORC
   partners reuse identical amplitudes/phases. Species/radial/atom/m output
   ordering and scalar expression order stay unchanged. Contiguous angle and
   fill passes expose independent rows to vectorization; no SIMD speed claim.
   Scratch: 8 L M_s + 32 L bytes (M_s=sum_c(2l_c+1) for one species), excluding
   the existing radial map and output. No G-by-atom phase table. Nonlocal
   diagonal already uses the same sparse complex table supplied to the real
   operator, reuses each adjacent ORC representative, and drops the source
   afterwards; no second projector quadrature or global-valued table is added.
4. **FFT plans:** REFSTART_FFT_SHARE already surrounds *field and operator*
   construction and forks the field workspace for DensityPlan. Equal kernels
   share handles; full-field and sparse wave supports may have different axes,
   batches and owners, so their layouts cannot be interchanged. First field
   planning is not removable by borrowing a later, different-support plan.
   A remaining duplicate is wave partitioning: owned_projector_indices,
   operator position_owners, then operator FFT each compute it.
   `REFSTART_PARTITION_CACHE=on` keeps the last exact partition in the existing
   constructor scope, avoiding the latter two constructions. Key: dimensions,
   ordered indices, rank count, mirror policy (including test override), balance
   policy. Partition is rank-independent; each plan still selects its rank and
   retains private work buffers. Backend is not part of this geometry-only key;
   kernel caching remains separately keyed by backend/length/batch. Only one
   partition is retained, released when start returns; metadata copying and
   zeroed work allocation remain. Cache scratch is O(G + sticks + ranks),
   excluded from the existing principal numerical-array envelope. No SCF
   apply-loop cache lookup, global cache, or shared mutable buffer is added.
5. **Preflight:** CLI validation and basis enumeration are outside this ticket's
   edits; retain their cost in the table. No setup-at-QE target is claimed met
   without the requested measurements.

## ONE profiled run per arm and case

The master relays to `baandi-dev08`; implementer does not access the server.
Use the frozen fa4cfad-snapshot **si64, si128, si216, MPI8** inputs, verify off,
all physical conditions and adopted switches unchanged. Confirm actual J/B/G,
not the nominal counts in the ticket. One OpenMP/MKL/engine thread per rank,
physical socket-0 cores 0-7, identical binding for BAANDI and the matched QE
MPI8 evidence. Keep default build flags (no new target-cpu experiment).

| Arm | REFSTART_CHOLESKY | REFSTART_SPECIES_FILL | REFSTART_PARTITION_CACHE |
| --- | --- | --- | --- |
| D | off | off | off |
| C | on | off | off |
| O | off | on | off |
| P | off | off | on |
| COP | on | on | on |

Exactly **15 profiled BAANDI invocations**, ONE per arm/case, <=600 s per job.
No warmups, A/A runs, repeats, automatic retries or precise campaign. Preserve
failed/timed-out jobs and elapsed cost. All arms run to the verify-off stopping
condition; they do not establish independent verification or default adoption.
New switches remain off until a separate review. No additional QE run requested.

```sh
cargo build --release -p dft-cli --features mpi,mkl
```

Freeze absolute binary/input/output paths and arm variables before launching:

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BAANDI_REFERENCE_SCF=on BAANDI_REFSTART_MKL=on \
BAANDI_REFSTART_GENERALIZED=on BAANDI_REFSTART_SINGLE_GEMM=on \
BAANDI_REFSTART_FFT_SHARE=on BAANDI_REFSTART_LOCAL_SHELLS=off \
BAANDI_REFSTART_CHOLESKY="$cholesky" \
BAANDI_REFSTART_SPECIES_FILL="$species_fill" \
BAANDI_REFSTART_PARTITION_CACHE="$partition_cache" \
BAANDI_DAVIDSON_STAGES=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
  /usr/bin/time -p taskset -c 0-7 \
  mpirun -np 8 --map-by core --bind-to core --report-bindings \
  "$baandi_binary" run "$baandi_input" --output-dir "$unique_output"
```

Record source commit/diff SHA-256, binary and input/UPF SHA-256, versions,
build flags, effective switch manifests, exit/completion status, start/end,
load, physical bindings, all-rank RSS, results, density.bin and raw console.
Ensure input tables do not silently override the frozen adopted configuration.
The listed verify-off input selection is essential; no verification toggle is
inferred from profiling switches.

## Adjacent setup clock attribution (seconds; MPI8 only)

Per explicit ticket request, put profiled setup next to existing matched QE
init_run as **diagnostic clock attribution**, not a production speed ratio or
proof of target attainment. Their boundaries/profiling differ. A missing QE
command/version/input/UPF hash means unavailable. Never use a mixed-model column.

| Case, MPI8 | D engine setup | C | O | P | COP | QE init_run | QE wfcinit (nested) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| si64 | pending | pending | pending | pending | pending | retrieve | retrieve |
| si128 | pending | pending | pending | pending | pending | retrieve | retrieve |
| si216 | pending | pending | pending | pending | pending | retrieve | retrieve |

Next to each case, include `setup_s` split for every arm, rank 0 plus min/max
and critical rank from all-rank records: CLI-inclusive and engine totals,
preflight, basis_gvectors, field_fft_plan, distribution_owned_rows,
nonlocal_tables, operator_fft_plans_and_tables, nonlocal_diagonal,
atomic_local_tables_structure_factors, atomic_orbital_fill,
atomic_overlap_gram, atomic_overlap_factor_inverse (D/O/P) or
atomic_overlap_factor_solve (C/COP), atomic_overlap_metric_check,
atomic_generalized_diagonalize_rotate, all QR/fallback rows, atomic_h_apply,
atomic_h_gram, atomic_h_apply_gram (inclusive), atomic_block_total (inclusive),
density_fft_plan and main_to_first_map. Never add parents/children,
collectives_nested, maxima of different ranks, or QE wfcinit to init_run.

For atomic H report J, B, L_r, chunk widths, actual pair counts, full H calls
and columns; split FFT/mirror/transpose, nonlocal projection/back accumulation,
and collective wait where existing profiling permits. Compare SCF **full-H**
columns at matching width/diagnostics; report kept-HX columns separately.
Do not divide a full-H time by a count including cheap kept-HX updates. The
setup_h_apply clock excludes outer X/HX allocation but includes lazy allocation
inside the FFT/operator on its first call; identify this before calling it a
per-column kernel regression. No new hot-loop clocks are added here.

O/P require exact start/trajectory bits for each build/rank count. C/COP require
start Ritz energies within 1e-12 Ha, orthogonality and selected-subspace defect
within 1e-12 on nondegenerate fixtures; equivalent degenerate subspaces may
rotate. Run the P4 start oracle on the reviewed MKL build as a correctness
check, separate from the 15 timed invocations. Preserve map/redo counts, H
columns, eigenvalue/density differences and final energy against the existing
independent strict reference; standard class remains 5e-9 Ha, strict 1e-10 Ha.
A verify-off sample is not independently certified. Later default promotion
requires the separate verification/generality gates; no new campaign here.

## Local validation and integration boundary

Mac pure-Rust correctness and simulated ranks only. No MKL/MPI runtime or QE
measurement. Final local results are appended after execution. The tests cover
cubic/skew/elongated/slab, species/projectors, odd bands, completion, and ranks
1/3/8. Analytic matrix tests independently check generalized eigenpair residuals
and the S metric. Plan tests check every cache-key component, nested scopes and
private buffers. All original tests remain unchanged except registry counts.

Central ledger/log, session/dashboard JSON, README generator/evidence are outside
OWNED_FILES. This request is the attempt/integration record for the master.
The checkout registry has no P4-START row and `.delegate/P4-START/report.json`
is absent; no nested delegate or dashboard upload was launched. No recurring
orchestrator timer is active in this implementer session. Leave changes
uncommitted; master records the gates and remote outcomes in central artifacts.

### Executed Mac gates

All exits below are actual process exits on the final code. Scratch/logs are
under `target/tmp/p4-start/` (untracked, not durable reference artifacts).

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | All workspace test targets checked |
| `cargo test -p dft-engine --test refscf -- --nocapture` | 0 | 55 passed, 5 existing ignored |
| `cargo test -p dft-solvers` | 0 | 187 passed, 25 existing ignored across targets |
| `cargo test -p dft-numerics` | 0 | 93 passed, 22 existing ignored across targets |
| `cargo test -p dft-numerics generalized_tests -- --nocapture` | 0 | 2 analytic/invalid-input tests passed |
| `cargo fmt --check` | 0 | Final format gate |
| `cargo clippy -p dft-engine -p dft-solvers -p dft-numerics --all-targets -- -D warnings` | 0 | Final lint gate |
| `BAANDI_TEST_TMP="$PWD/target/tmp" TMPDIR="$PWD/target/tmp" CHECK_SCCACHE=0 scripts/check.sh dft-numerics` | 0 | Build/test/fmt/clippy PASS, 93 passed and 22 ignored |
| `git diff --check` | 0 | Whitespace gate |

The gate's load-report `ps` call is sandbox-blocked (PermissionError); all
four actual check.sh gates ran and passed. No test threshold was changed.
No production/default flip, benchmark claim or new stored reference result was
made; README figures were not regenerated because their evidence is unchanged
and their generator/output are outside this ticket's ownership.

- P4 start comparison: four shapes x ranks 1/3/8 x four nonbaseline arms x
  thread budgets 1/2 = **96 candidate starts** plus 12 baselines. All field,
  density-G and diagonal bits match; O/P also preserve orbital/Ritz bits.
  C/COP maximum Ritz difference **8.881784197001252e-15 Ha**. Each candidate
  passes the unchanged 1e-12 orthogonality/selected-subspace bounds; each arm
  is bit-identical across thread budgets. The test asserts the new factor is
  actually used without QR fallback for J>=B, the old factor clock is absent,
  and the operator remains batched and usable after cache-scope destruction.
- Independent analytic generalized pencil tests, n=1/4/13/65: maximum entrywise
  eigenpair residual **9.925387274759139e-14**, S-metric defect
  **7.105427357601002e-15**, factor metric defect **5.551115123125783e-16**
  (all dimensionless). Known lowest eigenvalues agree within 1e-12.
  Near-dependent and indefinite metrics reject; malformed/nonfinite inputs
  and incompatible requested dimensions return errors.
- Two synthetic verified trajectories (cubic/skew), ranks 1/3/8, D versus COP:
  all 12 runs converge with independent verification passed. Largest final
  energy difference **3.19300141882195e-13 Ha** (cubic, 8 simulated ranks).
  COP takes 21 maps for cubic and 13 for skew. Duplicate same-position atomic
  functions reject collectively in both arms at every rank count. These are
  synthetic correctness tests, not Si/QE delivered-accuracy evidence.
- Partition-cache test covers dimensions, ordered support, rank count,
  mirror/greedy, balanced/unbalanced, private mutable workspaces, independent
  kernel-share-off operation, and nested scope restoration. Start construction
  reuses a wave partition twice; first field planning, copies and allocation
  remain deliberately visible in the existing setup clocks.
