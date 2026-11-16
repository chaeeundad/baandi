# MPI lazy orbitals: local stick rows as the SCF state (ticket MPI-LAZY-ORBITALS, 2026-11-17)

Switch: `BAANDI_M2_LAZY=on|off` (default **off** until dev08). Review item F07 (second half; the first half, the local rows of X from the solver, was MPI-KEPT-START section 7). Ledger rows MPI-LAZY-ORBITALS-*. dev08 request: docs/performance/requests/mpi-lazy-orbitals.md.

## 1. Problem
Every standard M2 solve ends with an allgatherv of the full real eigenvector block (N x B doubles; every rank receives all of it) and a complex Hermitian completion of it (N x B complex per rank). The node receive volume is P x N x B x 8 bytes per solve, so this step grows with the rank count while everything around it shrinks. dev08 microbenchmark (relayed, 2026-11-17, real MPI): si32 (11877 x 68) 557 / 1434 / 2050 us per call at np 2 / 4 / 8, node volume 12.9 / 25.8 / 51.7 MB; si64 (23871 x 132) 4551 / 7823 / 9994 us per call, 50 / 101 / 202 MB. With about 30 solves per si32 SCF run (fine steps, redo, verification) that is 60 ms at np=8 for si32 and 300 ms for si64, plus the unpack and the per-band copies of the complex block that follow every solve, and the ORC pack of the warm start before the next solve.

## 2. Consumer inventory (every reader of the gathered N x B block after an M2 standard solve)
Read from the code at 4080fa0. "Local" = can use this rank's stick rows directly; "reduction" = local rows plus a small collective; "full" = needs the whole block (gathered lazily at that point, once, cached in the handle).

| # | consumer | where | when (production standard path) | needs | lazy treatment |
|---|---|---|---|---|---|
| 1 | Complex unpack into `EigenSolution.eigenvectors` | `gamma_real::solve_real_gamma` | every solve | full | skipped: the solution carries a 0 x B placeholder and the k point carries the handle (`LazyOrbitals`) |
| 2 | Warm start (`initial`) of the next solve: ORC pack, row extraction, kept-start equality check | `solve_real_gamma` -> `solve_distributed_rows_*`, `workspace::solve_on` (`local_rows`) | every solve after the first | local | the handle's own rows enter the solver as the local start (`DistributedStart::Local`), after the elementwise ORC round trip `sqrt2 * (x / sqrt2)` (row 0 unchanged) that pack(unpack(x)) applies in the gathered path, so the start is bit-identical; one 1-double tree sum checks that every rank's handle matches its operator rows (else: materialize and take the gathered path) |
| 3 | Kept start X (MPI-KEPT-START) | `RealSolveState::m2_local`, `prepare_m2_kept_start` | every solve with `BAANDI_M2_REUSE=on` | local | already local (MPI-KEPT-START); unchanged |
| 4 | M2d density (band densities on the planes) | `ScfProblem::m2_density` | every map | local | already local (`m2_local`); unchanged |
| 5 | Replicated band density (`BandDensity::compute_keep`) | `step_with_occupations` | only when the M2 density declines (finite T, several k, grid mismatch) | full | materialize first (replicated condition) |
| 6 | Audit: per-band complex copies (`orbitals` vector) | `step_with_occupations` | every map | none for lazy | not built for lazy k points (empty columns of the right count; length checks only) |
| 7 | Audit: kinetic energy `sum_b f_b sum_G |G|^2/2 |c|^2` | `energy::evaluate_shared_with_partials` | every map | reduction | `sum_b f_b sum_r d_r x_rb^2` over the local rows with the ORC kinetic diagonal d (the preconditioner diagonal), one tree sum of the occupied band values: rounding class (different summation order; the complex form sums |c(G)|^2 and |c(-G)|^2 = (a^2 + b^2)/2 each, the ORC form a^2 + b^2) |
| 8 | Audit: nonlocal energy | `ScfProblem::m2_nonlocal_energy` | every map without `BAANDI_ENERGY_IDENTITY` | reduction | already distributed (one tree sum); if it declines, materialize |
| 9 | Audit: `prepare_shared` | `energy::prepare_shared` | every map | none (length check) | unchanged |
| 10 | Band-sum identity nonlocal (B2, `BAANDI_ENERGY_IDENTITY=on`) | `apply_band_sum_identity` | non-final maps | none | uses the distributed kinetic (rounding class) |
| 11 | Explicit nonlocal audit of the final state (B2) | `ScfProblem::explicit_nonlocal_energy` | final state with B2 | full | materialized with the final state (row 14) |
| 12 | Standard redo of the first step (warm start) | `session::advance` | first step when the estimate is already below the tolerance | local | handle passed with the placeholder |
| 13 | Verification map (ADR 0043) warm start, VERIFY-REUSE orbitals | `session::verify_fixed_point` | at every stop | local | handle passed; the verification's own handle is returned for VERIFY-REUSE |
| 14 | Final state: `ScfResult.orbitals`, `k_points[].solution.eigenvectors` (forces, stress, results/bands output, volumetric writers, Python/MCP results) | `session::advance` outcome block | once per run | full | materialized (one allgatherv per k block) before the result is built; result blocks are full complex matrices as before |
| 15 | Spectral audit (`spectral_guard`, off by default) | `audit_reported_spectrum` | at a stop when requested | full | materialize first |
| 16 | Certificate screen (strict only) | `certificate_shift_estimate` | strict path | full / local (`M2_SCREEN=dist`) | not lazy: strict solves always gather |
| 17 | Strict certificate re-solve | `session::advance` strict branch | strict path | full | not lazy |
| 18 | Ladder hand-off (coarse -> fine prolongation by Miller index) | `cutoff_ladder` via `ScfSession::warm_state` | once per ladder run | full | `warm_state` materializes (coarse block, once) |
| 19 | Checkpoint capture (`ScfSession::checkpoint`, exact continuation, result export) | `session.rs`, `runtime.rs`, CLI `checkpoint_io` | on request | full | materialized through the handle (cached; every rank runs the same CLI code, so the collective is entered by all ranks) |
| 20 | Checkpoint import / restore / Track A import | `import_electronic` | on request | (writes full) | produces full blocks; no handle |
| 21 | Snapshots / events / StepReport | `snapshot`, events | every step | none | carry no orbitals |
| 22 | Public `ScfProblem::step_kpoints*` (tests, stage-0 bench, external callers) | `scf.rs` | callers outside the session | full | never lazy (unchanged); the bench uses the hidden `step_kpoints_lazy` / `materialize_orbitals` |
| 23 | Collinear, finite-T, NSCF, direct minimization, M1 path, threads-only path | other drivers | | | never reach the M2 standard real path or never lazy |

Result: in the production standard loop (rows 1 to 13) nothing needs the full block; it is needed once at the end (row 14), at a ladder hand-off (row 18, coarse size) and on explicit requests (rows 11, 15, 19).

## 3. Design
- `gamma_real::LazyOrbitals` (engine-internal state of one k block, cheap to clone: `Arc`): this rank's stick rows of the returned ORC eigenvectors (a copy of the solver's `x_local`, n_loc x B), the global ORC row indices of these rows, the global length N, the half sphere (shared `Arc` per real state), the communicator, and a `OnceLock` cache of the materialized complex block. `materialize()` = three allgathervs (row counts, row indices, the rows; every rank enters them at the same session point) + the same Hermitian completion as the gathered path; cached until the handle is dropped (the next solve makes a new handle).
- Solver (`dft-solvers`): `RealDavidson::solve_distributed_rows_local` with `DistributedStart::{Global, Local}` and `gather: bool`; on the workspace path `gather = false` returns a 0 x B eigenvector placeholder and the local rows in `x_local`. `solve_distributed_rows_kept` is that entry with `Global` and `gather = true` (unchanged behaviour). On the non-workspace path a local start is gathered and the block is always returned (the engine checks the row count).
- Engine: `KPointSolution.lazy` (crate-private) carries the handle with the 0 x B placeholder; `ScfProblem::step_session` (session entry) passes warm handles and requests lazy returns when the switch is on; the public `step_kpoints*` never return lazy blocks. The session's continuation keeps `previous_lazy` beside `previous`; every boundary in rows 5, 8, 11, 14, 15, 18, 19 calls `materialize`.
- Eligibility (replicated): switch on, M2 active (two or more ranks, `BAANDI_MPI_M2` on), a distributed operator exists, standard profile (Ritz-change solve), workspace path (`BAANDI_M2_WORKSPACE` or `BAANDI_M2_SELECTED_ROOTS` on; default on), single k point, the M2 density and nonlocal energy available.

## 4. Identity contract
- Switch off: no code path changes (the session passes `lazy = false`, no handle exists).
- Switch on, fixed rank count: bit-reproducible (every decision replicated; the kinetic sum is a fixed-order local sum plus the deterministic tree sum).
- Switch on against off at the same rank count: solves, densities, SCF / Davidson / H counts, verification decisions and the materialized orbitals are bit-identical (the warm start is the same bits, gathers are exact data movement, the standard stop rule and the verification use the density metric and residuals only). The reported energies differ in the kinetic term by rounding only (distributed ORC sum; class 1e-14 relative), so `total_hartree` and the energy columns are rounding class, never the density.
- Accuracy classes unchanged (strict solves never lazy; standard 5e-9 Ha unaffected by a rounding-class kinetic sum).

## 5. Counters
`parallel.m2` (results.json, `m2_stats`): `lazy_returns` (solves that returned local rows only), `lazy_starts` (solves started from a handle's rows), `lazy_materializations` (gathers on demand) and stage `lazy_materialize`; `collectives.allgatherv.{calls,bytes}` shows the drop. Expected in the production standard loop: allgatherv N x B gathers 0 per step, one per k block at the end (plus one per ladder hand-off).

## 6. Screens (Mac, pure build, simulated ranks; counters exact C, wall observation O)
Labels as in docs/performance/mpi-scaling-gap.md. Simulated collectives are memcpy and thread barriers, so the time the real allgatherv costs on dev08 (section 1) is not in these walls; the counters and the identity are exact.

### Stage 0 (one solve at a fixed Hamiltonian, `fast_screen_mpi.py --level solve`, si32/si64 f1e-3 seed 601, 3 alternating pairs, 5 repeats, davidson_bench sha256 de8596ae940667e7...)
| cell | allgatherv bytes per rank per map (C) | allgatherv calls | cold solve B/A (O) | warm solve B/A (O) | identity |
|---|---|---|---|---|---|
| si32 np=2 | 11.04 -> 4.57 MB (-58.5 %) | 6 -> 5 | 0.991 (spread 1 %) | 0.994 (1 %) | eigenvalue bits equal, 68/68 converged |
| si32 np=4 | 11.04 -> 4.57 MB | 6 -> 5 | 0.971 (1 %) | 0.991 (1 %) | bits equal |
| si32 np=8 | 11.04 -> 4.57 MB | 6 -> 5 | 0.927 (5.4 %) | 0.971 (1 %) | bits equal |
| si64 np=4 | 34.36 -> 9.15 MB (-73.4 %) | 6 -> 5 | 0.969 (1.6 %) | 0.988 (1.5 %) | bits equal, 132/132 |

The removed volume per rank is the N x B block (si32 11877 x 68 x 8 = 6.46 MB, si64 23871 x 132 x 8 = 25.2 MB); what remains per map is the density plane gather and the XC potential gather (full grid each) plus the solver's per-solve row-index exchange (N doubles). H columns, Davidson iterations and GEMM FLOP identical; `m2.lazy_returns` 1 per map, `m2.lazy_starts` 1 per warm map.

### Stage 1 (production SCF: standard + verification + C1 ladder + capacity rule, `fast_screen_mpi.py`, si32/si64 f1e-3 seed 601, 3 alternating pairs, m2_entry_cost sha256 48bbdca3...)
| cell | allgatherv bytes (C, per rank, whole run) | calls | SCF / Davidson / H | density hash | dE (Ha) | wall A -> B (O, medians) |
|---|---|---|---|---|---|---|
| si32 np=4 | 224.4 -> 102.9 MB (-54.1 %) | 157 -> 137 | 13 / 74 / 2884 both | equal | 1.4e-14 | 4.052 -> 3.940 s (-2.8 %, spread 0.07 s) |
| si32 np=8 | 224.4 -> 102.9 MB (-54.1 %) | 157 -> 137 | 13 / 74 / 2884 both | equal | 2.8e-14 | 4.208 -> 3.857 s (-8.3 %, spread 0.18 s) |
| si64 np=4 | 455.3 -> 159.1 MB (-65.1 %) | 103 -> 92 | 10 / 51 / 4353 both | equal | 2.8e-14 | 12.968 -> 12.313 s (-5.1 %, spread 0.75 s) |
| si32 np=4, `M2_REUSE=on` both arms | 224.4 -> 102.9 MB | 157 -> 137 | 13 / 74 / 2272 both | equal | 1.4e-14 | 3.771 -> 3.601 s (-4.5 %, spread 0.07 s) |

Counters per run (si32): `lazy_returns` 26 (every standard solve of the coarse and fine stages and the verification), `lazy_starts` 24 (every warm solve), `lazy_materializations` 2 (ladder hand-off and result), `lazy_start_mismatches` 0, `solve_unpack` calls 26 -> 0; si64 17 / 15 / 2. Stage seconds per rank (O, si32 np=8): `solve_unpack` 0.052 -> 0 s, `solve_prepare` 0.053 -> 0.017 s (no warm-start pack), simulated allgatherv 0.074 -> 0.047 s, `lazy_materialize` 0.012 s, `lazy_kinetic` 0.003 s. Tree sums +50 per si32 run (one 1-double agreement check per lazy start, one B-value kinetic reduction per map).

Expected on dev08 (E, from the relayed per-call times, about 30 solves per run at si32 and 20 at si64, gather only): si32 17 / 43 / 62 ms per run at np 2 / 4 / 8, si64 91 / 156 / 200 ms, plus the unpack and pack (Mac O: 0.1 to 0.25 s per run per rank on si64) and the removed per-rank memory traffic of the N x B complex block.

## 7. Engine tests
- dft-solvers `tests/distributed_lazy_rows.rs`: local start bit-identical to the same start given replicated, lazy return placeholder 0 x B, `x_local` equal to the gathered rows bit for bit, eigenvalues and work equal, ranks agree (1 to 4 ranks, workspace full and selected roots, allocating path).
- dft-engine `tests/m2_lazy_orbitals.rs`: production SCF (verification on) cubic and skewed cells, P = 2 and 3, kept start off and on: density, eigenvalues, `ScfResult.orbitals`, `k_points[0]` eigenvectors and SCF / Davidson / H / verification counts bit-identical to the switch off, |dE| <= 1e-12 Ha (observed 0 to 1.8e-15), nonlocal energy bits equal, ranks agree, fixed rank count reproducible, one gather per rank; a checkpoint taken after three lazy steps holds orbitals bit-identical to the switch-off capture and the continued run delivers the same state.
