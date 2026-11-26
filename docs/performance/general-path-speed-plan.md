# General-path speed plan (GSPEED-PLAN, 2026-11-26)

Scope: the complex k-point reference path (`refscf/k_driver.rs`, `k_map.rs`,
`k_kept.rs`) and the spin/Fermi-Dirac path (`electronic.rs`,
`electronic_map.rs`), their start (`refstart.rs` k section), the complex
Davidson (`dft-solvers/src/ref_davidson_complex.rs`) and the k operator
(`dft-physics/src/distributed_complex.rs`). Analysis only; no code changed.
Base: main `d5938242`. User decision 2026-11-26: speed on the general path is
the top priority.

Labels: **C** = read from code (exact), **R** = relayed dev08 number,
**O** = Mac observation, **E** = estimate (assumptions stated). No timing in
this document is acceptance evidence. QE behaviour is taken only from the
own-words notes in `docs/research/qe-study/`.

Relayed starting point (R, dev08 MPI8, same inputs as QE 7.6):

| case | atoms | irreducible k | QE MPI8 total | BAANDI MPI8 |
|---|---|---|---|---|
| IrI3 | ? | ? | 132 s | 1194 s (SCF + forces + stress, gen41; Pulay assessment probably included) |
| GeTe | 2 | 260 (QE nsym 1) | 43 s | SCF 378 s, forces 23 s, stress 102 s |
| B12 | 12 | 172, grid 54^3, 22 bands | 84 s | SCF 560 s (12 maps) |

Pending from dev08 (kprof): BAANDI `BAANDI_PHASE_PROFILE` and QE clock
reports for the three cases. Section 7 lists what must be read from them
before any ticket is launched.

## 1. Path map (C)

| stage | file:line | what it does |
|---|---|---|
| dispatch | `refscf/api.rs:381-457` | policy present -> `electronic::run_continuation`; Gamma + `gamma_real` -> real driver; everything else -> `k_driver::run_continuation` |
| k start | `refstart.rs:1726-1778`, `1782-2040` | k bases, union grid, common fields; per k: owned rows, `NonlocalOperator::new_owned_with_rule`, `DistributedComplexHamiltonian::new_owned`, atomic columns, two-pass MGS, one H apply, Rayleigh-Ritz |
| k SCF loop | `k_driver.rs:378-1079` | ethr schedule (`next_ethr_ry`, QE-like), first-iteration redo, estimator stop, verification gate, Broyden mixing, permap fused sums |
| k map | `k_map.rs:53-288` | sequential loop over k (line 92); every k uses the whole communicator |
| per-k kept HX | `k_kept.rs:43-112` | exact X/HX/V cache; hit costs a both-way FFT of all bands (`add_local_change`) |
| complex Davidson | `ref_davidson_complex.rs:289-685` | owned-row block Davidson, generalized reduced problem, ideal reset |
| k operator | `distributed_complex.rs:23-473` | stick-distributed pair FFT, kinetic, nonlocal by two GEMMs + one tree sum |
| density symmetry | `refscf/symmetry.rs:166-200` | G-space projection, applied twice per map (second = diagnostic) |
| spin/FD map | `electronic_map.rs:500-922` | same per-(spin,k) loop, occupations after all solves, output fields every map |
| symmetry / k list | `dft-geometry/src/symmetry.rs:67`, `:247`; `dft-cli/src/config.rs:139-165`, `:639-642` | operation detection (tolerance 1e-7 Bohr, `config.rs:105`), orbit reduction with optional time reversal |

## 2. Coverage of the Gamma-path optimizations (C)

Has = present on that path, Partial = weaker form, No = absent, n/a = Gamma only by nature.

| Gamma item (ledger row) | k path | spin/FD path | note |
|---|---|---|---|
| Real wavefunctions, pair FFT, half sphere (B02/R-series) | n/a | n/a | only TRIM k could use it (backlog B18 G-H8, speculative) |
| Border-only S/H Gram + triangle (P3m, GEN-2432530 `REF_DAV_TRIANGLE`) | **No**: full m x m S and H Gram after every append (`ref_davidson_complex.rs:477`, `:552`) | No | rounding class to port |
| Ideal reset (P3m) | Has (`:594-644`) | Has | |
| sqrt preconditioner (P3j `REF_DAV_SQRT`) | Partial: per-element branch inside the row loop (`:441-452`) | Partial | |
| Fused append / check append / row cache / swap buffers (P3i, P3j) | **No**: eight row matrices and a GEMM workspace allocated per call (`:331-360`); `validate_finite` on both GEMM operands each call (`dft-numerics/src/linalg/tiled.rs:383-384`) | No | |
| Restart corrections, wide-rotate dispatch (P3) | No restart corrections; rotations go to MKL zgemm on mkl builds (`tiled.rs:432-437`) | same | |
| Error agreement per GEMM | one 1-double tree sum after every rotate/gram/apply (`:48-56`, `:138`, `distributed_complex.rs:248-264`) | same | about 10 latency-only collectives per Davidson iteration |
| Nonlocal single GEMM / fused coupling / back accumulate (P3l, P3m) | Partial: single projection GEMM, scalar coupling loop, separate back GEMM then an add pass (`distributed_complex.rs:266-354`) | same | |
| Atomic start: species fill, generalized Ritz, Cholesky, compact G, single GEMM, FFT share (P3l, P3t, P4) | **No**: radial quadrature per (G row, species channel, atom, m) (`refstart.rs:1861-1897`), column MGS with one tree sum per column pair (`:1947-1975`), no plan-cache scope, FFT layout built twice per k (`:1827`, `:1838`) | No, and the whole k start runs once per spin channel (`electronic.rs:747-759`) | largest code-visible fixed cost, section 4 |
| Lean energy (P3q) | Has (`k_map.rs:244-250`) but kinetic and nonlocal band energies still computed per k per map (`k_map.rs:168-182`) | **No**: output fields every map (`electronic_map.rs:809`) | |
| Permap fused sums (P3) | Has (`k_driver.rs:692-709`) | No | |
| Verification gate V2 (I25-VERIFY-COST) | Has (`k_driver.rs:472-479`) | No | |
| Verification warm start (P3-VERWARM) | Partial (kept HX, `k_driver.rs:528-551`) | Partial | |
| Kept HX (KEPTHX2) | Partial: no retained psi(r) slabs, hit = full both-way FFT of all bands (`k_kept.rs:75-82`, `distributed_complex.rs:361-400`) | Partial | saves only kinetic + nonlocal of the start block |
| Lazy output (P3l) | No: delivery allreduces zero-padded full orbitals per k (`k_driver.rs:1141-1166`) | No | once per run |
| Shared FFT workspaces / potential (P3t `FFT_SHARE`) | **No**: each k operator owns a full-grid `potential` and `delta` (`distributed_complex.rs:124-125`) and private FFT buffers | No, doubled for spin | B12 E: about 0.43 GB per rank of per-k full-grid arrays |
| Atomic start reuse at ionic steps | k seed exists (`k_driver.rs:78-86`) | No (`api.rs:304-312`) | |

## 3. MPI work distribution on the k path

BAANDI (C): no pools. Every k point is solved by all ranks with the plane-wave
(stick) distribution; k points run one after another (`k_map.rs:92`,
`ref_davidson_complex.rs:3`). Per Davidson iteration on one k the code issues
about 17 collectives (C, count from the code): rotate x/hx agreements (2),
correction norms (1), S Gram agreement + reduce (2), apply_block (entry
agreement, projection agreement, projection tree sum, back-GEMM agreement,
exit agreement, two FFT all-to-alls: 7), post-apply agreement (1), H Gram
(2), reduced solve agreement (1), plus the restart Grams. Per k per map
further: potential plane map, residual tree sum, nonlocal energy projection
and tree sum, density plane map (`k_map.rs:93-193`).

QE (note `qe-stage-by-stage.md:101`, `:134`): pools are used by default only
when the rank count exceeds half the number of z planes and nk > 1; otherwise
plane-wave decomposition, as BAANDI. With `-nk` each pool holds whole k points
and needs no communication inside a k solve; QE's best configuration for k
meshes is therefore expected to be `-nk` = ranks (the note marks it
unmeasured). The comparison rules (`comparison-rules.md`, item 2) require the
fastest eligible QE configuration of the same MPI model, so QE `-nk 8` is a
legitimate and probably faster reference than the 43/84/132 s numbers.

Expected cost of BAANDI's choice for 172-260 k on 8 ranks (E):
- Pure latency: about 17 collectives x 5-15 us (ADR 0041 tree-sum table:
  14 us for a 16 x 68 block at P=8) = 0.1-0.25 ms per Davidson iteration;
  GeTe 260 k x about 15 maps x about 4 iterations = 16k iterations ->
  2-4 s. Not the main term.
- Efficiency of small distributed kernels: GeTe and B12 have a few hundred to
  about a thousand rows per rank per k; zgemm, preconditioner and FFT slab
  work on such blocks run far from peak and every collective is a barrier.
  Plane imbalance: owned planes `ceil(n0/8)/(n0/8)` = 1.037 at n0 = 54,
  1.11 at 36. QE's own plane-wave-only k444s scaled 1.7x on 8 ranks
  (`qe-stage-by-stage.md:101`, relayed). A pool layout replaces all of this by
  serial per-k work plus one density reduction per map; load imbalance is
  `ceil(nk/8)/(nk/8)` = 1.015 (260 k), 1.023 (172 k) plus per-k cost spread.
  Estimate for the k-loop part at MPI8: 1.5-2.5x faster (E, to be replaced by
  the kprof phase split and a QE `-nk 1/2/4/8` scan).
- Memory: pools hold only nk/8 operators per rank (full rows each), about the
  same bytes as today's nk operators with 1/8 of the rows, minus the per-k
  full-grid arrays.

## 4. Per-k fixed costs (C, with estimates E)

1. **Atomic start radial quadrature per (row, atom, m)** (`refstart.rs:1861-1897`).
   `projector_transform_with_rule` is recomputed for every atom and every m of
   a channel; it depends only on (species channel, |k+G|). O: 12.7 / 17.0 /
   18.2 us per call for l = 0 / 1 / 2 on a 1534-point mesh (Mac M-series,
   release, scratch bench under `target/scratch-bench`, 20000 calls each; exit 0).
   B12 (PseudoDojo B: mesh 1534, 2S + 2P -> 48 columns for 12 atoms; N_pw per k
   assumed 6000-10000 from the 54^3 grid, 750-1250 rows per rank):
   172 k x 750-1250 x 48 x 13-25 us = **80-260 s** of the 560 s (E). With a
   per-(channel, |q| bits) memo per k the count falls by 24x (only 2 radial
   channels): 3-11 s. GeTe (Ge/Te: 3D, 4S/5S, 4P/5P -> 18 columns): 20-80 s (E).
   Bit-identical: the same call with the same q returns the same bits; the
   per-atom phase and the angular factor keep their product order.
2. **Projector table per k** (`hamiltonian.rs:2492-2495`): exact |q|^2 memo
   per species; at generic k every row has its own |k+G|, so one quadrature per
   (row, radial projector). B12 (4 radial projectors): 7-21 s; GeTe (6+6):
   10-30 s (E). An interpolation table (QE's approach, note
   `qe-hpsi-and-fft.md:14`, `:54`) was rejected in backlog B18 (H-6) because
   its 1e-8..1e-9 error exceeded the then 1e-12 Ha gate. New premise: the
   adoption gate is now strict 1e-10 / standard 5e-9 Ha and the per-k cost now
   matters (nk up to 260). A fine-table study is a class (b) item, not a
   default change.
3. **Start orthonormalization** (`refstart.rs:1947-1975`): column MGS, two
   passes, one tree sum per column pair: B12 about 2 x 48 x 47 / 2 + 48 = 2300
   collectives per k, 400k per run, 2-6 s latency (E) plus scalar
   O(c^2 rows) loops. The Gamma path's generalized Ritz start (P3l) avoids it.
   FD band growth (`refstart.rs:1642-1725`) has the same per-candidate pattern.
4. **Spin: the whole k start twice** (`electronic.rs:747-759`): two identical
   `start_electronic` calls build identical operators, projectors, FFT plans
   and start vectors; only the potential differs per channel.
5. **Davidson workspace per call** (`ref_davidson_complex.rs:331-360`): eight
   rows x cap matrices plus a GEMM workspace allocated and zeroed for every
   call (B12: about 6 MB per call, about 2.6k calls per run, E 1-3 s).
6. **Full S/H Gram per iteration** (`:477`, `:552`): O(m^2 rows) instead of
   O(m k rows), plus 2 m^2 doubles reduced per Gram instead of the border.
7. **Per-k per-map energy terms** (`k_map.rs:158-182`): residual pass + tree
   sum, kinetic pass, `nonlocal_band_energies` (one projection GEMM + tree sum)
   for every k in every map although lean maps use only band sums.
8. **Per-k plane maps** (`k_map.rs:93-101`, `:183-193`;
   `electronic_map.rs:541-549`, `:734-742`): a new `PlaneMap` and a new
   potential vector per k; an all-to-all of the full potential and of the
   density per k whenever the operator's first axis differs from the field's
   (`grid_dist.rs:404-407` is a copy otherwise). The potential is the same for
   every k of a channel.
9. **Kept HX on k** (`k_kept.rs:58-82`): bitwise compare of X, clone of HX,
   then a both-way FFT of every band; the saving is only kinetic + nonlocal of
   the start block. The FFT count per map equals QE's start h_psi.
10. **Symmetrization twice per map** (`symmetry.rs:178-181`): the second
    projection is a diagnostic only.
11. **Plan reuse**: no `PlanCache::enter` scope in `finish_kpoints`; each k
    builds private FFT lanes, transposes and slabs (`distributed_fft.rs:159-166`,
    `:660-685`).
12. **Eigen-tolerance and iteration counts**: the k path uses the Gamma
    schedule (`k_driver.rs:387-393`) and QE's per-band rule (empty bands
    `max(5 thr, 5e-6)`, `ref_davidson.rs:70-76`; QE note
    `qe-solvers-and-subspace.md:39`), retry rule max(5, bands/4) as QE. The FD
    path requires every retained band at the full threshold
    (`electronic_map.rs:563-564`), while QE loosens bands whose occupation is
    below 1 percent of the k weight (`qe-solvers-and-subspace.md:90`).
    Davidson iterations per k per map against QE's "average # of iterations"
    (typically 2-3) are not yet known: read them from kprof (BAANDI
    `records[].calls[].iterations`, QE output per SCF iteration).

## 5. k-point symmetry (C)

BAANDI: `detect` finds lattice-preserving integer rotations with fractional
translations at tolerance 1e-7 Bohr (`config.rs:105`), `reduce_mesh` keeps
the operations that preserve the (shifted) mesh and adds time-reversal
partners when `ktime_reversal` is true (default for k grids,
`config.rs:294-296`). Fractional translations incommensurate with the FFT
grid are not dropped: the grid is enlarged (`symmetry.rs:79-100`), or an
explicit grid is rejected.

Differences that can change the irreducible count against QE's nks:
- Tolerance: 1e-7 Bohr is tighter than QE's acceptance (the note does not give
  QE's value; to be checked from QE output "N Sym. Ops." on the same file).
  Materials Project coordinates rounded at 1e-6 can lose operations here and
  produce MORE k points than QE.
- **Collinear spin runs disable time reversal and symmetry reduction**
  (`config.rs:639-642`, `electronic.rs:634-638`): a spin-polarized run solves
  the full mesh, up to 2x (time reversal) times the group order more k than
  QE's LSDA run, which keeps both.
- QE drops incommensurate fractional translations; BAANDI keeps them (fewer k
  for BAANDI when it happens; the comparison then needs QE's grid).
- `compare_materials_qe.py:122` already requires equal k count and weights, so
  any timing pair that passed it is on equal k sets. GeTe has one operation on
  both sides (relayed); B12 and IrI3 must be checked from `symmetry.irreducible_k_count`
  against QE "number of k points" and from the spin mode of IrI3.

## 6. Ranked candidates

Gain columns are E for the three cases at MPI8, verify-off target metric,
before the kprof numbers. "Class": BI = bit-identical, RC = rounding class
(accuracy classes, fixed-rank reproducibility), TR = trajectory change (full
convergence on the generality set, strict 1e-10 / standard 5e-9 Ha, verification).

### (a) Overhead removals

| rank | ID | change | cases | gain (E) | effort | risk | class |
|---|---|---|---|---|---|---|---|
| 1 | GK1 | atomic start radial memo per (species channel, |q| bits) per k; per-atom phases only | B12, IrI3 (many atoms per species), GeTe | B12 70-250 s, GeTe 15-70 s | S | low | BI |
| 2 | GK2 | spin: build the k start once, share operators between channels (potential per channel set before each solve) | spin runs | half of spin setup and operator memory | S-M | low | BI |
| 3 | GK3 | complex Davidson: workspace reused across calls and k (row cache as `ref_davidson::work`), drop repeated `validate_finite` on known-finite operands, merge per-GEMM agreements into one per iteration | all | 3-8 % of solve | S | low | BI |
| 4 | GK4 | per-k energy terms only at audit (kinetic, nonlocal band energies deferred under lean), residual norms from the solver's own last residual where available | all k | 2-5 % | S | low | BI |
| 5 | GK5 | potential per channel redistributed once per map per distinct operator axis; density accumulated per axis in operator layout and moved once (k order kept) | all k | 1-3 % | S | low | BI |
| 6 | GK6 | shared full-grid potential/delta and FFT lanes/buffers across k operators (plan-cache scope in `finish_kpoints`) | all k | memory -0.4 GB/rank (B12), 0-3 % | M | low | BI |
| 7 | GK7 | second (diagnostic) symmetry projection opt-in; delivery orbitals by allgatherv | symmetric cases | <2 % | S | low | BI (diagnostic field changes) |
| 8 | GK8 | electronic path: lean energy, permap, verification gate as on the k path | FD/spin | 5-15 % of FD/spin maps (P3q showed per-map energy 1.0 -> 0.17 s at si128) | M | low | BI |

### (b) Algorithmic or layout changes

| rank | ID | change | cases | gain (E) | effort | risk | class |
|---|---|---|---|---|---|---|---|
| 1 | GK9 | **k-point pools**: stage 1 pool size 1 (each rank owns whole k points, `Serial` collective inside a k, operators with full rows), one density reduction per map, allgatherv of eigenvalues/residual sums for FD occupations, stop and verification; stage 2 general `npool` with a sub-communicator (MPI_Comm_split in the FFI adapter) for nk < ranks. ADR 0041 option D | all k cases | k loop 1.5-2.5x at MPI8 | L | medium (start, map, verification, checkpoints, derivatives on one layout) | RC across rank counts (ADR 0041 contract), fixed-rank reproducible |
| 2 | GK10 | border-only S/H Gram with triangle in the complex Davidson (port of the real path's border) | all k | 5-15 % of solve | M | low | RC |
| 3 | GK11 | generalized Ritz start for k (Gram + Cholesky/overlap factor + one Rayleigh-Ritz, no column MGS); FD band growth through the same block rule | all k | B12 2-6 s + loops | M | low | RC, start changes |
| 4 | GK12 | collinear spin: time reversal within each channel, crystal operations that map atoms onto atoms of the same initial moment; ADR with an independent full-mesh oracle | spin runs (IrI3 if spin-polarized) | k count / 2 or more | M | medium (physics policy) | TR + policy |
| 5 | GK13 | FD empty-band threshold as QE (bands with occupation below 1 % of the k weight use the empty-band threshold) | FD runs | unknown, H columns | S | medium (density accuracy of partially occupied tails) | TR |
| 6 | GK14 | fine radial interpolation tables for projectors and atomic orbitals (premise of B18 H-6 changed by the 1e-10 / 5e-9 classes) | many-k | 10-30 s per case | M | medium (accuracy) | RC/TR |
| 7 | GK15 | kept psi(r) slabs for k (bounded by a memory budget) | small grids | <5 % | M | memory | BI |

Not proposed (ledger): capacity sweeps (A08/A13/P31A), RMM default (P3r),
ladder (2026-11-26 retry), build flags (2026-11-26), locking (P24A).
GK14 repeats B18 H-6 only with the stated new premise.

## 7. First step: read kprof before code (no code)

Read from the dev08 phase profiles and QE clocks, per case:
1. BAANDI setup split: `complex_k_atomic_starts`, `nonlocal_tables`,
   `operator_fft_plans_and_tables` against `init_run` (QE). This confirms or
   kills GK1/GK11/GK14 (section 4 items 1-3).
2. Per map: `solve_s`, `density_s`, `potential_s`, `energy_s`, verification
   share; Davidson iterations and H columns per k per map against QE
   `c_bands` calls and average iterations; h_psi time per call.
3. IrI3: spin mode, smearing, Pulay assessment present or not (PULAY-OPTIN),
   BAANDI irreducible k against QE nks, forces/stress split.
4. QE layout line ("R & G space division", "K-points division: npool") of the
   relayed runs, and a QE `-nk 1/2/4/8` scan on GeTe and B12 (QE's best
   configuration is the target).

## 8. Staged plan (tickets)

Order: GK0 -> GK1 + GK3 + GK4 (one bit-identical batch) -> GK9 stage 1 ->
GK10 + GK11 -> GK2/GK8/GK12 for spin/FD -> GK9 stage 2 -> GK13/GK14 only if
kprof shows the terms. Each ticket: one mechanism, own switch for A/B and
rollback (adoption rule 2026-11-10), dev08 development protocol: one run per
arm, `BAANDI_PHASE_PROFILE=on` in both arms, same binary hash pair recorded,
cases GeTe, B12, IrI3 plus one Gamma check (si64) that must stay unchanged.

| ticket | owned files | acceptance |
|---|---|---|
| GK0 measurement | none (docs/performance only) | kprof table per case (section 7); QE `-nk` scan; decision which of GK1/GK9 leads |
| GK1 atomic start memo | `crates/dft-engine/src/refstart.rs` (k section) | k start orbitals bit-identical (test comparing column bits on a 2-species, 3-k synthetic fixture and on Si k222s); `complex_k_atomic_starts` drops on dev08; delivered density.bin identical |
| GK3 complex Davidson workspace and agreements | `crates/dft-solvers/src/ref_davidson_complex.rs`, `crates/dft-solvers/src/ref_davidson/work.rs` | bit-identical eigenpairs on existing solver tests and k222s/GeTe-size synthetic; error agreement tests (`complex_reference_agrees_nonfinite_operator_output`) still pass |
| GK4 deferred per-k energy terms | `crates/dft-engine/src/refscf/k_map.rs`, `k_driver.rs` | delivered energies and density bit-identical; intermediate audit fields null as in P3q |
| GK5/GK6 layout and sharing | `k_map.rs`, `electronic_map.rs`, `crates/dft-physics/src/distributed_complex.rs`, `refstart.rs` | bit-identical density; peak RSS per rank reported |
| GK9 k pools stage 1 | new `crates/dft-engine/src/refscf/k_pool.rs`, `k_map.rs`, `k_driver.rs`, `electronic_map.rs`, `refstart.rs` (k section), `forces.rs`/`stress.rs` k ownership (reuse MPI-DERIV ranges), checkpoint k state; switch `REFSCF_K_POOLS` | np=1 bit-identical to today; np = 2, 3, 8 against np=1 within strict 1e-10 / standard 5e-9 Ha, same k count, same map count or justified; fixed np bit-reproducible (two runs); nk < np falls back to today's layout; verification unchanged; generality set (cubic, skewed, elongated, slab, Gamma-in-list, shifted mesh, FD, spin) |
| GK10 border Gram | `ref_davidson_complex.rs` | RC: eigenvalue agreement 1e-12 Ha on solver fixtures, SCF classes on the k set, H columns within run-to-run equality or reported |
| GK11 generalized k start | `refstart.rs` (k section) | RC + full convergence on the k generality set; first-map H columns reported |
| GK2/GK8 spin and FD sharing | `electronic.rs`, `electronic_map.rs`, `refstart.rs` | bit-identical delivered results (spin and FD fixtures, FD band growth) |
| GK12 spin symmetry | `dft-cli/src/config.rs`, `electronic.rs`, `refscf/symmetry.rs`, ADR | full-mesh oracle agreement within strict class on a ferromagnet and an antiferromagnet fixture with symmetry-breaking moments rejected |

dev08 request template (per ticket): commits and binary hashes of both arms,
`mpirun -np 8` on socket 0 cores 0-7, one thread per rank, verify off for the
target metric and on for one check run, `BAANDI_PHASE_PROFILE=on` in both
arms, one run per arm per case, results.json + phase table + QE clock table
stored with sha256, accuracy against the case's strict reference.

## 9. Limits

- No dev08 profile of these cases was available; every gain above is E and
  must be replaced by kprof numbers before ranking is final.
- N_pw per k for B12 and GeTe is inferred from the grid, not read from inputs.
- QE's pool default and tolerance come from the notes; QE's actual layout in
  the relayed runs is unknown.
- The Mac quadrature timing is an observation on a different CPU.
