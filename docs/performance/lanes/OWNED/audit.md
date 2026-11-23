# G-OWNED-ROWS: construction and map-boundary ownership audit

Base: `f526588` (2026-11-23). Static source audit plus **previously relayed
 dev08** measurements; no Mac profiling or speed claim. GitHub CLI reads of
issues 20 and 24 failed (network); the connected GitHub read API supplied both
issues. The work order's line references describe an older implementation.

## Reconcile the premise before changing code

- `DistributedRealHamiltonian::new` already constructs local kinetic rows and
  directly packs local projector rows. `real_gamma::pack_projector_rows` retains
  global validation. The I20-OWNED-ROWS ledger entry records the always-on change.
- `RealSolveState` already separates its `Arc<HalfSphere>` from `Option<Operator>`.
  `ensure_serial_op` creates the serial real Hamiltonian on demand. The
  I20-LAZY-SERIAL-DEFAULT-ON entry records adoption on 2026-11-18: Si64 MPI8
  -1.22%, fine `solve.real_state` 0.1286 -> 0.0072 s in a 10-step screen.
  These are historical dev08 results, not this patch's measurements.
- `INIT_OWNED_STATES`, `M2_LAZY_SERIAL`, `M2_LOCAL_FIELDS`, `M2_LOCAL_DENSITY`,
  `M2_LAZY`, and `M2_LAZY_RESULT` are already on. Do not add duplicate switches,
  revert their defaults, or credit their savings to this task.
- On stable maps `ensure_state_with`, `DistributedRealHamiltonian::same_operator`
  and `grid_plan` already reuse geometry and FFT plans. Potential epochs change;
  geometry does not. No evidence of repeated geometry/plan construction on every
  map justifies a new cross-problem cache. The fine and coarse ladder problems
  have different bases/plans and coexist by design. Item (c) is deferred.

## Payload notation and frequency

All sizes below are **logical live array payload bytes per rank**, excluding
Vec/Arc headers, alignment, allocator capacity, FFT library opaque plans and MPI
buffers. A replica of X bytes costs P*X bytes across P ranks of one cell. For
owned blocks sum actual local lengths (never round G/P and call it a peak).
G = full wavefunction PW/ORC rows; H=(G+1)/2 representatives; D = density PW;
F = Nx*Ny*Nz grid points; L_r = local ORC rows; F_r = local plane points;
Q = angular projector columns; B = requested bands; C = solver capacity;
A = atoms; J = atomic-span columns; h = mixer history. usize=8, f64=8,
Complex64=16 bytes on both target machines. GVector/Momentum payloads use
`size_of::<T>()`, not an assumed packed C layout. Length sums below are exact
payload formulas even where a backend's actual capacity is unavailable.

Construction frequency is per problem/k-state (normally fine + coarse once per
run), unless stated otherwise. A solve can occur on an ordinary map, a standard
redo, or a verification attempt. Thus 'per solve' must not be multiplied only
by the final history length. Allocation totals are not simultaneous peak RSS.

| Object / source site | Per-rank payload and lifetime/frequency | Ownership finding |
|---|---|---|
| CLI input/UPF parsing, preparation and preflight | input bytes + radial tables O(radial points/projectors); basis/grid preflight arrays O(G+D); per rank/run, sometimes repeated in preparation | Required metadata is replicated; full radial input is not a row block. SHA cache already on; no new I/O work here. |
| `ScfProblem::construct`: Gamma, density, k-state bases and force-model descriptors | sum of G/D vector lengths times Rust element sizes; grid index map 24G per copy; per problem; `gamma.clone`, `k_grid.clone`, force-model copies | Full metadata still replicated; useful for validation, ionic construction, fallback and output. Gamma construction reuse already implemented. |
| ionic G coefficients and ionic real grid | 16F + 8F (ionic reciprocal data is a full FFT Array3), per problem; coarse/fine separate | Full arrays on all ranks. Distributed ionic construction is a separate experiment; do not alter phase arithmetic here. |
| Ewald keys/ions and atom/species/PP descriptors | O(A) positions/charges; cached scalar Ewald, per problem | Cache already checks cell/charge/position. Not a per-map Ewald evaluation. |
| initial density, optional core/XC density, stop-metric arrays | 8F per real field, 16F per reciprocal stop array; per problem with map updates; `starting_density` clones 8F | Full input/read APIs and fallbacks remain. `StopKept`/Coulomb buffers are FFT-grid sized; replacing them needs consumer contracts. |
| complex Hamiltonian kinetic, potential and basis metadata | 8G + 8F + basis metadata; one per state, potential updated on gathered maps/consumer request | Still replicated even with sparse nonlocal data; complex atomic HQ/fallback can use them. |
| complex Hamiltonian serial FFT scratch | sum of plan Vec lengths times 8/16, plus opaque plan storage; one per worker/state | Full-grid serial plans remain needed for atomic HQ/fallback, not the lazily removed **real** operator. FusedPlan backend controls exact scratch; do not claim 16F is its total. |
| complex nonlocal source | owned representative rows 16*H_r*Q, local real operator 8*L_r*Q; full fallback materialization 16GQ; coupling 8*sum(q_atom^2) | Sparse source already default. Full table appears only at consumers needing it; validation/radial tables remain replicated. |
| `BandDensity::new` in SCF setup | **8F(B+1)** for B contribution grids and total density, plus serial FFT plans and optional 16G per worker phi; per state/problem | Largest identified idle serial payload on successful M2 density path. Requires lazy wrapper/fallback work in `band_density.rs` and consumers (outside OWNED_FILES); follow-up priority. Zero-filled virtual allocation is not proof of resident RSS. |
| `HalfSphere` / ownership setup in distributed H | 4*(2H+G) descriptor; temporary owners 8G, slots 8H, position lookup 8G, representative tuples `size_of::<(usize,bool)>()*G`, links total <=8H plus P² Vec headers | Created once per operator, not per map. Direct owned kinetic/projectors already done; temporary metadata is O(G), not O(GQ). |
| distributed H kinetic/projector rows | 8L_r + 8L_r Q + coupling above | Already owned; do not repeat original issue #20 fix. |
| distributed H potential | **8F**, with unowned entries NaN in owned mode; one allocation/operator; local plane writes each map | Layout/API requires full indexing (`potential`, snapshot, kept delta, FFT apply). Candidate physical storage 8F_r, but requires FFT/snapshot API work outside scope. |
| distributed FFT and mirror metadata/scratch | per-field sum(size_of(element)*Vec length); H image 8L_r C; hbuf 16*ceil(C/2)*local_FFT_positions; mirror x send+recv 32P*max_link*C, h send+recv 32P*max_link*ceil(C/2); P^T X/weights 24QC | Mostly owned or required communication buffers, sized once/operator/capacity. Dense Q*C values legitimately replicated for reduction. Metadata/cache rewrite is not established as a time lever. |
| `solve_real_gamma` diagonal bridge | full 8G plus selected 8L_r; each solve | Still builds full ORC diagonal although distributed H has local kinetic. Small candidate: Si128 full temporary 382,808 B/rank; not selected alone. |
| **full complex initial -> full ORC packed -> selected rows** in `solve_real_gamma` | source 16GB remains; **8GB + 8L_r B** bridge matrices per supplied full start, including cold guess / materialized handoff / fallback from mismatched lazy rows | **Selected:** default-off M2_OWNED_INPUT builds only 8L_r B for workspace standard M2. Removes 8GB per eligible entry/rank. No new geometry cache. |
| ordinary lazy warm start / returned local blocks | O(8L_r B) copies and kept X/HX; per map/attempt; local density and psi planes owned | Already avoids global ORC packing. New switch must report zero calls on a purely lazy input; no invented saving on those maps. |
| kept potential/delta in `prepare_m2_kept_start` | old potential 8F_r in owned mode (else 8F), **delta 8F** with owned plane values, temporary current/old planes <=16F_r, kept image 8L_r B | Delta full indexing remains. Independent local-storage follow-up; potential epoch checked before apply. |
| session initialize trial mixer | full starting n 8F, trial residual/input/result and distributed mixer gather scratch O(F); validated then reset; per session including coarse/fine | 0.235 s envelope is not proof all of it is avoidable. `initialize_owned` has a separate handoff validation path already. Must preserve validation errors, normalization and collective agreement. |
| ladder materialization/prolongation | source/target complex 16G_coarse B and 16G_fine B, gathered real rows up to 8G_coarse B; density 8F_coarse/8F_fine; optional h pairs of difference fields | Once per handoff. Default owned handoff remains off. This patch only removes target ORC packing, not source/target complex blocks or interpolation. |
| atomic span Q/HQ and Ritz | 16G_start J per complex block; local row blocks ~16G_start J/P; 16J² reduced matrix; rotated guess 16G_start B | Full Q gather/complex HQ remain. Separate INIT experiment; no orthogonalization change. |
| grid field plan and density/mixer | owned work/planes/sticks, full-grid layout/weights; mixer state ~16hF_r plus own scratch; per plan and map | Numerical grid/field split already exists. Do not treat replicated small metadata or canonical reduction buffers as removable density replicas. |
| terminal density, orbitals, serialization | density allgather 8F/rank, byte conversion 8F/rank, P*8F file output per cell; optional full orbital 16GB plus gather scratch; per delivery/consumer | Lazy orbital result already default for CLI consumers needing no orbitals. Density still delivered on all ranks, output contract belongs to CLI. No root-only output change in this patch. |

## Rank by measured seconds, not by byte size alone

Numbers are relayed, different commits/windows, inclusive and **not additive**.
Source: `requests/si64-gap-attribution-2026-11-20.md` (279278d), INIT diagnosis
(c0f2143), GLUE/LAUNCHIO diagnoses (e947a35). No new QE run.

| Priority/envelope, Si64 MPI8 | Seconds | Interpretation for this lane |
|---|---:|---|
| atomic start excess / INIT atomic total | +0.41 / 0.52558 | Large, separate arithmetic/representation follow-up; not warm packing's budget. |
| other init + between + unassigned excess | +0.39 | Residual bucket, not uniquely avoidable ownership work. |
| session setup / ladder handoff | 0.235 / 0.123 | Highest boundary envelopes; full fields and complex imports persist. Do not add these to residual bucket. |
| io/materialize excess; LAUNCHIO materialize / density write | +0.30; 0.23333 / 0.20749 | Mostly all-rank output; outside this core patch. |
| solve outside solver envelope excess | +0.77 | Broad budget, not all construction. Components setup 0.122, start 0.170, return 0.051, finish 0.040, restart 0.025, real_state 0.014, entry 0.009 are nested/partial, not a closure. |
| complex-to-ORC input packing | not isolated; inside prepare/setup and entry | Cannot justify a 0.15 s forecast. Selected by memory criterion, not by claiming the whole 0.122 s. |
| INIT projector assembly rank max | 0.00442 | Already owned construction. Not a plausible new 0.15 s lever. |

**No remaining scoped item has a substantiated expected >=0.15 s gain.**
Do not close on that basis alone: removing the full ORC input block saves
99,530,080 B = 94.919 MiB per rank at Si128 (G=47,851, B=260), and
25,207,776 B = 24.040 MiB at Si64 (G=23,871, B=132). These are allocation
payload reductions per eligible entry; live/RSS improvement depends on allocator
and peak overlap. Across MPI8 the Si128 removed payload is 796,240,640 B.
The idle BandDensity grids are larger still: Si128 1,558,683,648 B/rank
(8*72*72*144*(260+1)), a separately owned follow-up, not a claimed saving here.

## Implementation and invariants

`M2_OWNED_INPUT` defaults OFF independently of all adopted switches. It only
routes a full complex input through local packing when the existing standard
M2 workspace path can take `DistributedStart::Local`. It keeps the same input,
rows, bands, solver, tolerances, reduction order and output. Unsupported strict,
M1/serial, or non-workspace paths retain global packing. Lazy and replay-local
starts already local remain untouched. Snapshot capture still receives the
same local rows. No new allocation inside H/FFT/solver numerical loops.

`pack_orbital_rows` deliberately follows **HalfSphere::pack**, not projector
packing (which has a different zero-mode tolerance). It scans every coefficient,
zero-mode imaginary component, every conjugate pair, and scaling overflow on
every rank, including empty owners. The output is written only after validation.
Only selected rows receive sqrt(2)*Re/Im. Tests compare independent full packing,
including signed zero, subnormals, threshold equality, NaN/Inf and overflow.

No geometry/charge/position cache is introduced: each call borrows the current
operator's rows and half sphere. Existing operator identity invalidates new
Hamiltonians (including cell, atoms, charge/PP and basis), grid and solver key
changes. The lifetime is a solve, so there is no new invalidation policy. Active
session communicator replacement is **not** a supported remapping API: existing
GridPlan checks dimensions, not communicator identity. Changing rank count/layout
requires a new problem/session. Do not claim safe mid-session rank migration.
The tests reconstruct across cell shapes, cutoffs, ranks and synthetic PP types;
existing lazy-serial tests exercise same-operator potential reuse and fallback.

## Counter contract and validation

Under `INIT_PROFILE=on`, `run.counters["owned_lane.input.*"]` reports per-run,
per-rank successful standard M2 full-start bridge payloads (entered calls also
count if validation fails): calls, direct_calls, allocated_bytes, baseline_bytes,
saved_bytes. Allocated bytes count actual requested f64 matrix elements at the
pack and extracted-diagnostic-start allocation sites, **not allocator/RSS samples**.
For successful entries with shapes (G_i,L_ri,B_i):

- baseline_bytes = sum 8*(G_i+L_ri)*B_i in both arms;
- A allocated_bytes = baseline_bytes, direct_calls=saved_bytes=0;
- B allocated_bytes = sum 8*L_ri*B_i; saved_bytes=sum 8*G_i*B_i;
- B allocated_bytes + saved_bytes = baseline_bytes; direct_calls=calls>0.

Aborted invalid packing counts the entered allocation but not a later extraction;
therefore the successful-entry identity is intentionally not asserted on failure.
No clocks or counter atomics run when INIT_PROFILE is off. Totals exclude the
solver's independent persistent start copy, complex source, and all other objects.
See `validation.md` for executed gates and `dev08-request.md` for pending timing.

### Serial FFT payload detail

On dev08's MKL build `FusedPlan` is `GridFusedPlan`: two complex FFT grids
are **32F bytes**, plus 24G index bytes and the backend plan/scratch, per serial
worker. The pure Rust `OrderedFusedPlan` has 24G index bytes plus 8 times the
sum of scatter/planes/plane_sticks/active_x/forward_z/gather lengths and 16
times the sum of sticks/partial/forward/real_plane/kernel scratch/lane lengths.
Those support-dependent lengths need backend counters before an exact RSS
prediction. The unused density contribution grids alone already dominate the
Si128 follow-up; no opaque plan-memory saving is included in its 8F(B+1) estimate.

The new fast invalidation test explicitly replaces cell, position, PP charge,
basis cutoff and FFT dimensions, requires a new descriptor identity, and checks
potential-only reuse. For each variant reconstructed layouts at ranks 1/2/3/8
cover every ORC row once. This is new-session/layout coverage, not active-session
migration support.
