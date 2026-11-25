# P4-MEM: allocation census and lifetime trimming

Status: scoped implementation and Mac numerical checks passed; **full allocation
coverage is blocked by ticket ownership**. No remote run or performance
acceptance is claimed. All changes remain uncommitted. No QE source was read.

## Finding and implementation

The existing `reference_scf.peak_bytes_per_rank` is an analytic envelope, not
an allocation measurement. It remains labelled as such. It cannot explain the
campaign RSS difference by subtraction.

The inspected allocation owners are:

| Owner | Large storage and lifetime | This patch |
| --- | --- | --- |
| refstart atomic block | raw atomic span, Gram/factor/rotation arrays and its GEMM workspace, only during start | Already dropped on return; no claim of removing a retained atomic block |
| StartState | density_g, first potential, initial Ritz energies survive all maps although unused by this driver | Release before the first map; ionic and static diagonal stay live |
| reference Davidson | V/HV/X/HX at capacity, E/HE at band width; reduced matrices and private LeanScratch | Preserve repair capacity; return band-sized storage without a copy or a different GEMM |
| Davidson row cache | V/HV and E/HE survive between solves and formerly after delivery | Cache four capacity-wide arrays in the candidate; return two band-wide arrays. Release cache before fine delivery/coarse handoff |
| RealGemmWorkspace | one short-by-short tile for **every 256 long rows**, plus packing storage | Largest unmeasured candidate; private to dft-numerics, outside ownership |
| Hamiltonian | real projectors, image, projected/reduced blocks, mirror buffers, hbuf, GEMM workspace, full-index potential | Private to dft-physics, outside ownership |
| distributed FFT | sticks, middle, slab, growing batched send/recv, lane scratch, kernels and metadata | Field and density plans have distinct workspaces; kernels are already shared. Private to dft-solvers/distributed_fft.rs, outside ownership |
| kept HX | cloned active X/HX, owned potential, retained complex pair slabs; temporary supplied HX/delta | Release supplied HX after solve, preserve slab work/reuse, release all cached storage before delivery |
| Broyden | previous x/F plus dx/dF history | **Already owned half-sphere**, not full-grid history. No representation or summation-order change |
| map/verification | input fields, field cache, output density, old orbital guess, candidate/check maps | Release old guess and dead potential/diagonal immediately after solve; operator density/plane map after redistribution |
| output | two replicated real gather buffers, unpacked complex block, two public complex copies | Release gather input before unpacking; release mixer/kept/cache storage before delivery. Public copies remain API-owned |

For diagnosis only, the GEMM allocation rule in `real.rs::with_pool` is
`8 * ceil(long/256) * short^2` bytes of initialized tiles. The Davidson uses
`long=max(local_rows,capacity), short=capacity`. The Hamiltonian uses
`long=local_rows, short=max(projector_count,bands)`. Atomic setup has its own
workspace with `short=max(atomic_columns,bands)`. These are **source formulas,
not measured Si128/Si216 numbers** and are never serialized as measurements.
Single-GEMM dispatch does not reduce these constructor allocations. Shrinking
this shared workspace safely needs an adapter change respecting its fallback
and repair paths; it is not possible by merely changing the Davidson bound.

`BAANDI_REFSCF_MEMORY_TRIM=on` / `[performance] refscf_memory_trim='on'`
enables the lifetime and return-storage changes above. Default **off** until
dev08 validation. No arithmetic, tolerance, capacity, FFT pair count, history
length or slab-retention default changes. `REFSCF_KEPT_RECOMPUTE` remains off:
KEPTHX2 adopted retained slabs; the P4 recompute screen saved memory but
increased Si128 solve/wall time. This patch does not repeat that experiment.

The constructor still has four capacity-wide row arrays plus two band-wide
arrays. The compact arm reuses all four capacity-wide arrays across solves;
its two band-wide arrays become returned X/HX after the last repair. This
reduces new row initialization on a warm solve from two capacity-wide arrays
to two band-wide arrays. The cache itself grows by the complementary amount;
do **not** present the return-storage reduction as equal whole-process RSS
savings. It reduces overlap with later solves, and the cache is explicitly
released before output. Empty ranks and changes of cache configuration are
handled without a numerical branch.

## New results.json block and coverage boundary

`reference_scf.memory_breakdown`, schema `reference-memory-1`, contains:

- `component_names` and rank-ordered `components_per_rank`, in bytes:
  actual Vec capacities observed at constructor/map/delivery boundaries;
  known newly allocated matrix storage uses its initialized backing length.
- `davidson_names` and rank-ordered `davidson_per_rank`: V, HV, X, HX, E, HE,
  accessible reduced arrays, vectors, returned backing and cached rows.
  Capacity/return ownership is tracked at construction and recycling, not
  inferred from the final contracted active matrix shape.
- `complete=false`, with an explicit `unmeasured` list. Components are
  independent maxima; **do not sum them into a simultaneous peak**, subtract
  that sum from RSS, or call it whole-run allocator accounting. Coarse-stage
  storage, inaccessible scratch and allocator overhead are not covered.
- `start_vectors_released` is released payload, not additional live storage.
  Returned/cached Davidson bytes overlap the workspace categories in time.

No I/O, JSON, global allocator hooks, per-element logging or numerical-storage
allocations were added to the hot loops. Fixed-size thread-local counters are
scoped to the run stage and restored on errors/nesting. Rank gathering occurs
once at delivery.

To meet the ticket's **every large allocation** criterion, a follow-up needs
read-only storage introspection in these currently unowned locations:
`dft-numerics/src/real.rs` (matrix backing, GEMM/EVD),
`dft-solvers/src/lean_davidson.rs`,
`dft-solvers/src/distributed_fft.rs`,
`dft-physics/src/distributed_real_gamma.rs`, and backend FFT adapters as needed.
Atomic-start instrumentation also exceeds refstart.rs's "freeing arrays only"
permission. Library-internal plan allocation must remain unavailable unless
measured by an approved allocator/profiler; no size formula can substitute.
The scope request is pending. The present patch is useful partial attribution,
not completion of that criterion.

## Exactly ONE run per arm and case on dev08

Master session relays this request through `baandi-dev08`. Use the frozen
Si128/Si216 MPI8 inputs from `/data/sig-dft/precise-2026-11-26/`, with the
current adopted P4 start switches in **both** arms. Select the verified input
variant (independent fixed-point verification on, its cost included). Freeze
exact input/UPF/binary hashes and the complete switch manifest. No new QE run,
no one-core or thread/MPI comparison, no warmup, repeats or automatic retries.

| Case/model | Baseline arm D | Candidate arm M | Count |
| --- | --- | --- | --- |
| si128 MPI8 | refscf_memory_trim=off | refscf_memory_trim=on | 1 each |
| si216 MPI8 | refscf_memory_trim=off | refscf_memory_trim=on | 1 each |

Exactly **four BAANDI invocations**, one case/arm per job, 600 s timeout each.
Keep all other settings identical: cell, coordinates, cutoff/grid, bands,
projectors, occupations, start, convergence threshold, verification, output,
kept-HX/recompute/ladder settings, toolchain/build flags and physical binding.
Use production/light diagnostics in both arms, including this lightweight
census. Keep REFSCF_KEPT_RECOMPUTE=off and REFSCF_LADDER=off. No silent input
TOML override of an environment arm. Reject an arm whose results manifest
reports the wrong effective value.

```sh
cargo build --release -p dft-cli --features mpi,mkl
```

Use the campaign's existing per-rank GNU time wrapper and process-wall timer;
do not substitute launcher RSS for rank RSS. Bind eight MPI ranks to the same
eight physical cores as the frozen campaign, with one engine/OpenMP/MKL thread
per rank. For a new output directory and frozen paths in each job:

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BAANDI_REFSCF_MEMORY_TRIM="$arm" \
  timeout 600s /usr/bin/time -p taskset -c 0-7 \
  mpirun -np 8 --map-by core --bind-to core --report-bindings \
  "$rank_time_wrapper" "$baandi_binary" run "$baandi_input" \
  --output-dir "$unique_output"
```

Here `rank_time_wrapper` is the existing campaign wrapper, writing GNU
`time -v` maximum RSS into **distinct rank files** through process exit.
Record units explicitly: GNU time KiB multiplied by 1024 for bytes; decimal
MB and MiB must not be mixed. Preserve output-preparation VmHWM from results
as a separate boundary from the complete-process external peak.

Return raw results.json, per-rank RSS, total wall, input/binary/UPF/source/diff
SHA-256, launcher/binding, switch manifests, all scientific output hashes,
SCF and verification counts and failed-attempt cost. Require identical energy,
density and orbital bits where materialized, per-map dr2, solver iterations,
H columns, correction decisions when available, verification attempts and
convergence. A successful one-run screen is provisional, not a measured median
or a multi-system adoption result. If wall regresses beyond known run noise,
retain off; no additional runs are authorized here.

| Case/model | Arm | Peak RSS bytes by rank | Component census | Wall s | E/density/work identity |
| --- | --- | --- | --- | --- | --- |
| si128 MPI8 | D | pending dev08 | pending dev08 | pending | pending |
| si128 MPI8 | M | pending dev08 | pending dev08 | pending | pending |
| si216 MPI8 | D | pending dev08 | pending dev08 | pending | pending |
| si216 MPI8 | M | pending dev08 | pending dev08 | pending | pending |

These are intentionally unfilled: Mac synthetic fixtures cannot supply actual
allocation/RSS numbers for the remote Si workloads. If full attribution is
required before the screen, resolve the ownership expansion above first.

## Mac checks and attempt evidence

Base commit: `f9cfd1a93798cf125d69b752921be7e34c5bedc7`. No Mac timing
supports a speed claim. New tests compare both arms at identical rank/thread
counts; the candidate also preserves thread-count identity. All four full SCFs
pass independent fixed-point verification at a 1e-12 Ry SCF threshold. All
energy differences are **0 Ha**; density, orbitals, eigenvalues, per-map dr2,
Davidson corrections/work and verification decisions match bit for bit.
The rollback arm also restores the baseline orbitals after cache-key changes.

| Synthetic fixture | Simulated ranks | Delivered E (Ha), both arms | Returned X/HX backing bytes per rank, D -> M | Released start vector bytes per rank |
| --- | ---: | ---: | --- | --- |
| cubic, local | 1 | -1.10301826312977003 | 1824 -> 912 | 12320 |
| skew, nonlocal | 3 | -0.839575255872317294 | 672,576,576 -> 336,288,288 | 4384,4256,3744 |
| elongated, nonlocal | 8 | -0.765046225087195220 | 288,576,192,192,192,192,0,0 -> 144,288,96,96,96,96,0,0 | 792,872,872,872,872,512,512,352 |
| slab, local | 3 | -1.07153420657288700 | 864,960,768 -> 432,480,384 | 3272,3416,3128 |

Each row was checked with candidate engine thread budgets 1 and 2. These are
small synthetic systems, not Si measurements or independent QE references.
The zero-owned-row ranks in the elongated fixture are deliberate. A separate
solver test covers capacity smaller than 2n, repeated cached solves, cache
disabling/rollback and empty owners. It probes actual initialized backing with
`reshape_within_storage`: the measured width succeeds and one more column
fails, independently of the census formula. Solver repairs retain full scratch.
Lazy delivery tests also cover trimming with deferred materialization and an
explicit orbital request, at two and three simulated ranks.

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace/test type check |
| `cargo test -p dft-engine --test refscf` (final rerun) | 0 | 62 passed, 7 ignored |
| `cargo test -p dft-engine --test refscf memory_trim_ -- --include-ignored --nocapture` | 0 | 2 passed, four shapes, 1/2 threads, numerical evidence above; 0.18 s debug test execution |
| `cargo test -p dft-solvers` | 0 | 188 passed, 25 ignored across targets |
| `cargo fmt --check && cargo clippy -p dft-engine -p dft-solvers -p dft-cli --all-targets -- -D warnings` | 0 | Clean |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | 9 passed; shared generated evidence remains current |
| `git diff --check` | 0 | Clean |

The broader memory test initially carried an ignore attribute for runtime;
a measured 0.18 s for both new tests justified enabling it in the normal gate.
No tolerance or production algorithm was adjusted to obtain a pass.

The first standard gate exited **1** (786 passed, 1 failed, 179 ignored): the
new CLI assertion caught a JSON placement mismatch. The census was serialized
inside the pre-existing `reference_scf.run` summary while the documented API
and test required `reference_scf.memory_breakdown`. The adapter now moves that
object to the documented location, without duplicating arrays. Build/fmt/clippy
passed; the failing test was `run::tests::reference_ladder_records_and_fine_stop_clock`.
The targeted command
`BAANDI_TEST_TMP=/private/tmp/p4-mem-json TMPDIR=/private/tmp cargo test -p dft-cli --lib run::tests::reference_ladder_records_and_fine_stop_clock`
then exited **0** (1 passed), checking real written results.json for both
verified and explicitly unverified ladder runs.

Final standard gate:

`BAANDI_TEST_TMP=/private/tmp/p4-mem-gate CHECK_SCCACHE=0 CHECK_NICE=0 scripts/check.sh fast dft-engine dft-solvers dft-numerics dft-cli`

exited **0**, with **PASS build/test/fmt/clippy**, **75 binaries,
787 passed, 0 failed, 179 ignored**, 48 s gate wall. This includes the final
lazy-delivery and actual JSON serialization tests. The auxiliary `ps` load
report was denied by the sandbox on both gate invocations; it did not prevent
any gate from running. Raw local logs: `/private/tmp/p4-gate-final.log`,
`p4-engine-final.log`, `p4-solvers.log`, `p4-memory-tests.log`, `p4-json-test.log`.

Final owned-file audit: 16 changed/new files, all in the ticket's allow-list;
`git diff --check` exit 0. Changes are intentionally uncommitted. No existing
ledger, dashboard, sessions registry, README or generated figure was edited.
This checkout has no P4-MEM entry in `artifacts/sessions.json` and no local
`.delegate/P4-MEM/report.json`; registry/dashboard reconciliation belongs to
the master session. No recurring orchestrator timer was installed here.

Disposition: **BLOCKED for full acceptance**, with tested partial attribution
and an opt-in reduction ready for review. Complete large-allocation coverage
and measured Si128/Si216 MPI8 breakdowns remain outstanding; no peak-RSS
reduction or absence of runtime regression has been established on dev08.

Shared README/ledger/dashboard files are outside this ticket's ownership;
this owned request is the attempt record for the orchestrator to reconcile.
