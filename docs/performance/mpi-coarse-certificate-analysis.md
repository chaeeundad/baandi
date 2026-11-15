# MPI coarse stage and certificate: the non-scaling part of the M2 wall

Work order: docs/performance/requests/work-order-beat-qe-2026-11-14.md, phase E2. MPI against MPI only; no threads number appears here.
Evidence tags: (R) relayed from dev08, (O) Mac observation (pure build, simulated ranks = threads of one process, single runs, not dev08 evidence), (C) read from the code or counted exactly, (E) estimate.
Status 2026-11-15: analysis complete; fixes 1 and 2 implemented behind switches that are OFF by default (`BAANDI_M2_SCREEN=dist`, `BAANDI_M2_ATOMIC=split`); dev08 measurement pending (request in section 7).

## 1. Starting numbers (R)
si64 strict f1e-3 seed 601, MKL, one thread per rank, arm 1 (M2 default switches), closure of results.json:

| seconds | np=2 | np=4 | np=8 | np=2 -> 8 |
|---|---:|---:|---:|---:|
| total | 73.75 | 51.52 | 46.07 | 1.60x |
| phases (fine session) | 52.84 | 33.04 | 25.14 | 2.10x |
| ladder_coarse | 15.92 | 12.49 | 12.18 | 1.31x |
| certificate outside phases | 3.58 | 4.49 | 6.94 | 0.52x (slower) |
| of which certificate screen | 3.553 | 4.460 | 6.904 | 0.51x |
| coarse potential (atomic start inside) | 2.348 | 2.491 | 2.773 | 0.85x |
| coarse atomic start | 2.263 | 2.420 | 2.700 | 0.84x |
| coarse diagonalization | 12.4 (E, from 1.48x) | | 8.4 | 1.48x |
| setup | 0.87 to 1.00 | | | |

si32 strict screen 0.429 / 0.506 / 0.747 s against outside-phases 0.436 / 0.514 / 0.759 s. The certificate block count is fixed (si64 3, si32 2), so the growth is per-block screen cost; the certificate solve itself (phases) scales 2.47x. At np=8 the non-scaling items are coarse diagonalization 8.4 s + screen 6.9 s + coarse atomic start 2.7 s (45 percent of the wall).

## 2. What runs where on the M2 path (C)
One thread per rank (`threads = 1`), so every replicated item runs single-threaded on every rank at the same time.

| item | code | M2 (np >= 2) | collectives |
|---|---|---|---|
| coarse setup (problem at E1) | `cutoff_ladder::coarse_problem` -> `ScfProblem::new` | replicated | none |
| coarse atomic start: span (radial transforms, columns, two-pass Gram-Schmidt) | `atomic_start::span` | replicated, serial | none |
| coarse atomic start: H Q | `atomic_start::ritz` (complex `Hamiltonian::apply_block`) | replicated, serial: ATOMIC-REAL is excluded when a communicator is active (`gamma_real::atomic_apply_real`, `SharedComm::active().is_some()`) | none |
| coarse atomic start: Gram Q^H HQ, EVD, rotation | `atomic_start::ritz` | replicated | none |
| potential: Hartree | `input_fields` -> `hartree::evaluate` | replicated (full-grid FFT pair) | none |
| potential: LDA XC | `m2::xc` | distributed (M2d) | 1 tree_sum + 1 allgatherv |
| Davidson (coarse, fine, certificate solve) | `solve_real_gamma` -> `solve_distributed_rows` | distributed (M2c), reduced eigenproblem replicated | see mpi-reduced-eigen-analysis.md |
| density, nonlocal energy | `m2_density`, `m2_nonlocal_energy` | distributed (M2d) | alltoall, allgatherv, tree_sum |
| mixing, prolongation | session, `cutoff_ladder::hand_off` | replicated | none |
| certificate screen: fields | `certificate_shift_estimate` -> `input_fields` | Hartree replicated, XC distributed; reused by the certificate solve (POTENTIAL-REUSE) | 1 tree_sum + 1 allgatherv |
| certificate screen: single-band shifts | `M2Density::weighted_integrals` | distributed | 1 tree_sum |
| certificate screen: degenerate clusters | `certificate_shift_estimate` (cluster loop) | replicated, serial: per cluster band a scatter and a complex 3D FFT on the full fine grid, then s(s+1)/2 full-grid complex products, then a small EVD | none |
| certificate checks | `density_distance`, `check_gap` | replicated, cheap | none |

Nothing is root-only. No rank waits on a root.

## 3. Why these items do not scale or get slower
### 3.1 Screen (C, O)
The screen issues 3 to 4 collectives per block (XC tree_sum and allgatherv, the integrals tree_sum); their per-rank time falls with ranks (Mac `screen_integrals` stage 0.045 s at np=2, 0.017 s at np=8). The time sits in the replicated cluster loop. Mac si64 (O, per rank, 3 blocks, 127 cluster bands and 692 entries per block):

| Mac si64 per rank | np=2 | np=8 |
|---|---:|---:|
| cluster band FFTs (381 full-grid complex FFTs) | 0.94 s | 1.95 s |
| cluster entries (2076 full-grid products) | 0.46 s | 1.23 s |
| fields + integrals | 0.07 s | 0.07 s |
| screen total | 1.48 s | 3.27 s |

The same single-thread work takes 2.1x longer per rank when 8 ranks run it concurrently: the working set (127 complex full-grid fields) is streamed from memory by every rank, so the copies share the socket's L3 and memory bandwidth (the mechanism of MPI-SCALING-BEYOND-8 and MPI-REDUCED-EIGEN). dev08 shows the same shape (3.55 -> 6.90 s, 1.94x). It is not collective latency.

### 3.2 Coarse atomic start (C, O)
Entirely replicated and single-threaded (table 2); the H Q falls back to the complex apply because ATOMIC-REAL is threads-only. Mac si64 per rank: span 0.55 / 0.78 s (two-pass Gram-Schmidt of 256 atomic columns is 90 percent of it), H Q 0.50 / 0.72 s, Gram 0.09 / 0.12 s, rotation and EVD 0.06 s, at np=2 / np=8: the same contention growth (1.4x). The coarse potential phase is the atomic start plus 0.07 to 0.08 s (R: 2.348 - 2.263 s).

### 3.3 Coarse diagonalization (C, O, E)
Mac si64: 80 of the 102 Davidson iterations of the run are coarse (6 coarse SCF steps). Per Davidson iteration the distributed work scales with the coarse grid (54^3 coarse grid: about a third of the fine grid's points at half the cutoff), but the per-iteration costs that do not depend on the cutoff stay: the replicated dense reduced eigenproblem (order up to 6 x bands = 792 for si64, the same at both cutoffs), the per-column tree_sum reductions of the orthonormalization (count fixed by the trajectory, latency growing with log2 P) and the per-iteration allgathers. With 78 percent of the solves in the coarse stage, about 0.78 x 7.20 = 5.6 s of the dev08 np=8 reduced_eigen (R, MPI-REDUCED-EIGEN) is coarse (E): most of the 8.4 s coarse diagonalization at np=8. This is the strong-scaling limit of the replicated reduced problem, made more visible because the coarse stage has a third of the distributed work per iteration.

## 4. Fixes, ranked by expected seconds at si64 strict (E unless marked)

| rank | fix (switch) | rounding | np=2 | np=4 | np=8 | status |
|---|---|---|---:|---:|---:|---|
| 1 | screen clusters from the distributed band values: `band_values` on each rank's FFT planes (the stick FFT of the operator, one batched transpose), plane-local products, one tree_sum per cluster (`BAANDI_M2_SCREEN=dist`) | screen estimate only, relative 1e-16 level; the certificate meaning is unchanged (the screen can only reject, a pass still needs the full solve) | -2.7 | -4.0 | -6.4 | implemented, off |
| 2 | atomic start split over ranks: span columns built cyclically, two-pass Gram-Schmidt distributed with the serial arithmetic (one broadcast per column, first passes right-looking on the owners), H Q and the Gram columns in contiguous blocks with one allgatherv of the small matrix (`BAANDI_M2_ATOMIC=split`) | span bit-identical; Gram block GEMM rounding class | -0.9 | -1.4 | -1.9 | implemented, off |
| 3 | rank-0 reduced solve and B5 (`BAANDI_M2_EIGEN=root`, `BAANDI_M2_B5=on`, already implemented): their effect is concentrated in the coarse stage (78 percent of the solves) | exact / rounding class | 0 | | -0.5 to -1.5 | implemented, off; measure the coarse diagonalization column |
| 4 | distributed dense EVD (ScaLAPACK `pdsyevd` / ELPA) | rounding class, fixed per P | | | -3 to -4 (coarse share of -3 to -5) | needs the dev08 microbenchmark and library adoption (MPI-REDUCED-EIGEN) |
| 5 | recursive-doubling tree sum (power-of-two P) | exact | | | -0.3 | not done (MPI-REDUCED-EIGEN item 3) |
| 6 | H Q of the atomic start through the distributed real Hamiltonian (ATOMIC-REAL on M2) | rounding class | | | -0.1 to -0.2 after fix 2 | not done (small after fix 2) |
| 7 | distributed Hartree, prolongation, coarse setup | | | | under -0.3 together | not proposed (Mac: coarse setup 0.15 s, prolongation 0.02 s per rank) |

Fixes 1 and 2 at np=8: about -8.3 s of 46.07 s (E); the certificate block and the atomic start would then scale with ranks. The fixed residue (Hartree, fields, the serial second-pass chain of the Gram-Schmidt, which is at most half of its work) is under 1 s at np=8.

Expected numbers derive from the Mac per-rank ratios applied to the dev08 seconds: screen (O) np=2 1.48 -> 0.40 s, np=8 3.27 -> 0.24 s; atomic start (O) np=2 1.20 -> 0.72 s, np=8 1.69 -> 0.51 s.

## 5. Mac observations (O, pure build, single runs, energies and densities bit-identical in every row)
`crates/dft-engine/tests/m2_entry_cost.rs` with `M2_PROFILE_LADDER=on` (C1 ladder, Si.pz-vbc, 15 Ha, f1e-3 seed-601 harness displacements); "on" means both switches.

| case | wall off -> on | screen per rank off -> on | atomic start per rank off -> on | certificate wall per rank off -> on |
|---|---:|---:|---:|---:|
| si32 np=2 | 9.49 -> 9.11 | 0.365 -> 0.091 | 0.187 -> 0.119 | 1.64 -> 1.26 |
| si32 np=4 | 7.29 -> 6.97 | 0.409 -> 0.070 | 0.218 -> 0.089 | 1.21 -> 0.86 |
| si32 np=8 | 7.43 -> 6.72 | 0.557 -> 0.060 | 0.271 -> 0.079 | 1.30 -> 0.76 |
| si64 np=2 | 57.93 -> 55.35 | 1.477 -> 0.403 | 1.201 -> 0.723 | 10.92 -> 9.89 |
| si64 np=8 | 40.57 -> 36.34 | 3.399 -> 0.218 | 1.688 -> 0.512 | 8.28 -> 5.20 |

Bit identity (O): energy bits and a density FNV-1a hash of si32 np=2 and np=3 are equal for the HEAD build (git archive of e554622) and for this branch with every arm (off, screen, atomic, both): 0xc05fad02a3081918 / 0x2d32453b6b90a5c8 (np=2), 0xc05fad02a30818eb / 0xce90111edf820c14 (np=3); si64 np=8 off and both: 0xc06fac198a263f80 / 0xeccc7ea398bec6af. Work counts identical in every row (si64: SCF 11, Davidson 102, H 9023; coarse Davidson 80).

## 6. Determinism and accuracy contract (ADR 0041)
- Switches off: the code paths are the previous ones (the replicated screen code is unchanged inside an `else`; `band_densities` now shares its body with `band_values` through a closure that evaluates the same expression `w * v * v / volume`, tested bit for bit). Pure-build SCF bits equal to HEAD (section 5).
- `M2_SCREEN=dist`: the conditions (M2 active, distributed operator and M2 density present, cluster boundaries from the replicated eigenvalues) are replicated, so every rank enters the same transforms and the same tree_sum. Plane-local sums in fixed order plus the binomial tree_sum: bit-reproducible at a fixed rank count. The estimate differs from the replicated one by rounding only; the screen can only reject (margin 1.02), so the certificate meaning is unchanged.
- `M2_ATOMIC=split`: the Gram-Schmidt gives each column exactly the serial sequence of operations (bit-identical to `orthonormalize` at every P, tested); a nonfinite column travels as a status word so all ranks return the same error and no rank waits. The Gram matrix is a GEMM of column blocks (rounding class; observed one ulp of the total energy at P=5 on the pure build). Fixed P: bit-reproducible.
- Neither switch changes anything on one rank (threads-only path) or with `BAANDI_MPI_M2=off`.

## 7. Tests (Mac, pure build)
- `crates/dft-physics/tests/distributed_real_gamma.rs` `distributed_band_values_match_serial_and_band_densities`: plane-local psi equal to the complex full-grid transform within 2.8e-16 to 3.5e-16 relative (limit 1e-13), weight psi^2/volume equal to `band_densities` bit for bit, values independent of the rank count; cubic, skewed and elongated cells, ranks 1 to 5.
- `crates/dft-engine/src/atomic_start.rs` `split_two_pass_gram_schmidt_is_bit_identical_to_the_serial_routine`: ranks 2 to 5, with a zero and an exactly dependent column.
- `crates/dft-engine/tests/m2_coarse_certificate.rs` (simulated communicator, switches forced per rank by `m2::set_tuning_override`): `distributed_screen_matches_the_replicated_screen` (one-atom cubic cell with a degenerate occupied cluster, two-atom skewed cell; ranks 2 and 3: screen estimates equal to the replicated ones within 1e-10 relative, SCF |dE| against the threads-only run 8.9e-16 to 3.6e-15 Ha (strict limit 1e-10), energy bits equal to the replicated-screen arm, ranks bit-equal, reproducible, distributed clusters taken) and `split_atomic_start_matches_the_replicated_start` (s and p atomic orbitals; ranks 2, 3, 5: |dE| against the threads-only run at most 5.3e-15 Ha, against the replicated M2 start 0 to 2.7e-15 Ha, same atomic orbital and H counts, ranks bit-equal, reproducible).
- `crates/dft-numerics/src/switches.rs` `mpi_entry_cost_switches_defaults` (both new switches); `crates/dft-cli/tests/run.rs` switch count 58 -> 60.
- New results.json fields (`parallel.m2`): stages `screen_clusters`, `atomic_start`, `op.band_values`; counters `screen_cluster_bands`, `screen_dist_clusters`, `atomic_split_starts`, `band_value_builds`.
- Profile: `DFT_SI_UPF=<Si.pz-vbc.UPF> M2_PROFILE_LADDER=on M2_PROFILE_CASE=si64 M2_PROFILE_RANKS=2,8 M2_PROFILE_SERIAL=off cargo test --release -p dft-engine --test m2_entry_cost -- --ignored --nocapture` (prints energy bits and density hash).

## 8. dev08 request (MPI against MPI only)
Build: this branch with `--features mkl,mpi` (MKLROOT=/opt/intel/oneapi/mkl/2026.1), one binary, sha256 recorded; MKL_NUM_THREADS=1, one thread per rank, `--map-by core --bind-to core` on socket 0, load gated as in BAANDI-MPI-RANK-SWEEP.
Cells: si32 and si64, strict and standard, f1e-3 seed 601, ranks 2, 4, 8. Baseline (arm 1) measured in the same session three times, interleaved with the other arms.
Arms (same binary, environment A/B):
1. defaults (both switches off),
2. `BAANDI_M2_SCREEN=dist`,
3. `BAANDI_M2_ATOMIC=split`,
4. both,
5. both plus `BAANDI_M2_EIGEN=root BAANDI_M2_B5=on` (coarse diagonalization column).
Report per row: wall; the closure table (setup, phases, certificate_outside_phases, step_overhead, ladder_coarse, remainder) with `certificate.screen_s`, `certificate.wall_s`, `certificate.blocks`, `ladder_coarse.{setup_s, session_s, prolong_s, atomic_start_s, phases}`; `parallel.stages`; `parallel.m2.stages` (`screen_clusters`, `atomic_start`, `op.band_values`) and counters (`screen_dist_clusters`, `atomic_split_starts`); `parallel.collectives` (bcast calls appear with arm 3); SCF / Davidson / H counts; total energy |dE| against arm 1 at the same rank count (strict 1e-10 Ha, standard 5e-9 Ha); every rank's energy equal to rank 0; density SHA of two repeats of arm 4 at np=4 (bit-reproducibility at fixed rank count). Standard runs have no screen (the screen runs before strict certificates only) and their verification step is unaffected by arm 2: they measure arm 3 and must equal arm 1 otherwise.
Also run `cargo test --features mkl -p dft-engine --test m2_coarse_certificate` (MKL rounding of the Gram block GEMM).
Adoption (AGENTS.md adoption rule): each switch whose arm improves the median wall at np=2, 4 and 8 without a correctness failure becomes the default with its variable kept for rollback.
