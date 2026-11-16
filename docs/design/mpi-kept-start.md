# MPI kept start on the stick rows (ticket MPI-KEPT-START, 2026-11-16)

Switch: `BAANDI_M2_REUSE=on|off` (default **off** until dev08 stage 2; rounding class). Ledger rows MPI-FS-KEPT-START (finding) and MPI-KEPT-START (this change). Serial precedent: C2 REUSE-PSI (ledger C2-REUSE-PSI, adopted on dev08 as C2-REUSE-PSI-DEV08-ADOPT). Review item absorbed: F07 (local rows of the returned X).

## 1. Problem
At equal SCF and Davidson counts the M2 path applied H to about 29 percent more columns than the serial path (si32 production: 2699 against 2087), because every standard M2 solve re-orthonormalized the warm start and applied the full distributed H to it (pair FFTs both ways, mirror exchanges, projector GEMMs with a tree sum, kinetic pass), and then the engine re-extracted the local rows of the gathered eigenvectors (`extract_local`) for the distributed density and energy. The serial path has, since C2, a kept start that replaces the initial apply by a potential change on kept real-space orbitals.

## 2. What is kept (per rank, per k state)
| item | owner | shape | lifetime |
|---|---|---|---|
| X: local rows of the returned eigenvectors | engine `RealSolveState::m2_local` (unchanged field; now filled by the solver) | local rows x bands | until the next solve of the state |
| H X: return rotation of the H images, (H V) C | `DistributedKeptStart` in `M2Reuse` | local rows x bands | stored by a Ritz-change return, consumed by the next solve |
| V_old: potential of the stored H X | `M2Reuse::v_old` | full grid (replicated) | set at every keeping return |
| plane slabs psi(r) of the leading fully occupied pairs | `M2Reuse::planes` | kept pairs x plane points (complex, re = band 2t, im = band 2t+1, unscaled inverse transform) | written by the M2 density build after the keeping return |
| (V_new - V_old) X scratch | `M2Reuse::image` | local rows x bands | allocated once per state |

Memory per rank: X and H X add 2 x 8 x n_loc x bands bytes (n_loc about n / P), the slabs 16 x kept pairs x plane points (plane points about grid / P). si64 at np=4 (grid 60^3 = 216000 points as a sizing example, 64 kept pairs): slabs about 16 x 64 x 54000 = 55 MB per rank, X and H X a few MB. Memory gate: `PSI_KEEP_MAX_BYTES` (2 GiB, the serial bound) applied to the replicated estimate kept pairs x ceil(grid / P) x 16 bytes, so every rank takes the same decision; above it the kept H X is still used and every pair is transformed both ways (counter `reuse_gated`).

## 3. The start
Before a standard M2 solve of the same operator (`prepare_m2_kept_start`): if the previous solve stored X and H X, the density build ran after it (`collect`), the solve has a warm start and the grid matches, the engine computes delta = V_new - V_old on the full grid, and if any entry is nonzero forms (V_new - V_old) X with `DistributedRealHamiltonian::apply_local_change`: the kept pairs take their slab from `planes` (multiply, r -> G only), the remaining pairs (empty bands, an odd tail, a gated or missing slab) are packed and transformed both ways; one mirror-row alltoall only if a pair is not kept, one transpose to r only if a pair is not kept, one transpose back, the mirror-h alltoall. No kinetic or nonlocal work. H X += that image; the start is armed.

The solver (`solve_distributed_rows_kept`, workspace path) decides with replicated data: (a) one tree sum of the local mismatch flags of `initial` against X (shape, or any entry differing by more than `KEPT_START_TOLERANCE` = 1e-14, the serial value), (b) if all ranks match, the Gram defect of X (one bands x bands Gram, the `restart_defect` site) must be at most 1e-12, the contraction rule the solver already uses at a restart. Taken: basis = X, products = H X, no orthonormalization and no initial apply (bands columns saved). Declined: the plain start (orthonormalize the local rows of `initial`, apply H). Unlike the serial lean solver, the distributed solver needs an orthonormal basis, hence (b) instead of a factored overlap.

## 4. Invalidation (everything is dropped or not armed)
- New operator: `ensure_state` rebuilds the state (different key, pool or `update_from_complex` false: new atoms or projectors, basis, grid), or the distributed operator is rebuilt (`state.dist` re-created): `m2_reuse = None`.
- Ladder coarse <-> fine: the stages are different problems with different states; nothing crosses.
- Failed solve: the stored return is consumed at entry (`stored = false`) and set again only by a successful keeping return.
- Strict (residual) or certified solve, M1 path, threads-only path, non-workspace distributed solver (`BAANDI_M2_WORKSPACE=off` and `BAANDI_M2_SELECTED_ROOTS=off`): `m2_reuse = None` (the switch has no effect there).
- Exact checkpoint: the existing `forget_reuse` flag (set by an exact-continuation capture) also drops `m2_reuse`, so a continued and a restored run agree bit for bit, as for REUSE-PSI.
- `initial` absent or not equal to X (restored session, density transfer, other warm start): declined (counter `reuse_declined`) or not armed (`reuse_not_kept`).
- Two solves without a density build in between: X and H X are still valid, the slabs are not (kept pairs reset to 0 at every return): every pair is transformed both ways.

## 5. Reproducibility contract
- Fixed rank count: bit-reproducible. Every decision (armed, taken, kept pairs, potential changed, memory gate) is a function of replicated data or of one tree sum, so all ranks enter the same collectives in the same order; every arithmetic step has a fixed order (the kept slab is the slab the inverse transform produces, `apply_local_batch_kept` has the bits of `apply_local_batch` for the same slab, tested at 1 to 5 ranks).
- Against the switch off: rounding class (basis X instead of the re-orthonormalized pack/unpack round trip of X; H_old X + dV X instead of H_new X, equal to 1.2e-16 to 2.4e-16 relative in the physics unit test). Accuracy classes unchanged (strict 1e-10 Ha is not affected: strict solves never take the kept start; standard 5e-9 Ha).
- Different rank counts: rounding class, as every M2 change.
- Density bits do not depend on keeping slabs (tested).

## 6. Interaction with other switches
- `BAANDI_M2_WORKSPACE` / `BAANDI_M2_SELECTED_ROOTS`: the kept start lives in the workspace solver only (both return branches store X and H X; with selected roots the return already rotates all bands). The workspace keeps its own buffers; X and H X are copied into it at the start (no aliasing between the engine's matrices and the per-thread workspace).
- `BAANDI_M2_DENSITY=fused|bands`: both density forms keep slabs (`band_density_sum_keeping`, `band_densities_keeping`).
- `BAANDI_M2_PANEL_GS`, `BAANDI_MPI_TREE`, `BAANDI_M2_BATCH`, `BAANDI_M2_DGEMM`: independent (the potential change is always batched; the extra collectives are one 1-double tree sum and one bands x bands Gram tree sum per armed solve).
- `BAANDI_M2_SCREEN`: strict only, never meets the kept start.

## 7. F07 scope (what changed, what remains)
Done: the local rows of X come from the solver (`x_local`), no `extract_local` pass over the gathered n x bands block. Remaining (deferred, not in this ticket): the N x B allgatherv at the return and the complex unpack stay, because the engine's SCF state carries complex full orbitals that several consumers read every step (session warm start, checkpoint capture, ladder prolongation, explicit audits, output). Making them lazy is an engine-wide type change (a `DistributedRealState` with explicit materialization at those boundaries); it needs its own ticket. The allgatherv volume is measured in the screen (`collectives.allgatherv.bytes`) so that the remaining share is known (profile: allgatherv 1 to 4 percent of the Mac wall).

## 8. Counters
Engine (`parallel::m2_stats`, results.json `parallel.m2`): `reuse_armed`, `reuse_taken`, `reuse_declined`, `reuse_unchanged_potential`, `reuse_not_kept`, `reuse_columns_saved`, `reuse_plane_keeps`, `reuse_gated`, `reuse_plane_bytes_max` (a maximum, not a sum; the screen harness divides it by the rank count, read it raw). Operator: `inverse_pair_transforms`, `forward_pair_transforms` (apply, density/values, potential change), `kept_plane_pairs`, `potential_change_columns`, `applied_columns` (full applies only).

## 9. Stage-1 screen (exact counters; Mac pure, simulated ranks, production si32/si64, f1e-3 seed 601)
| cell | np | H columns | inverse pair FFTs | tree sums | alltoall bytes | allgatherv bytes | taken / armed |
|---|---:|---|---|---|---|---|---|
| si32 | 2 | 2699 -> 2087 (-22.7 %) | 3402 -> 2766 (-18.7 %) | 17263 -> 13225 (-23.4 %) | -5.5 % | 0 % | 20 / 22 |
| si32 | 4 | 2884 -> 2272 (-21.2 %) | 3687 -> 3021 (-18.1 %) | 17955 -> 13716 (-23.6 %) | -4.3 % | 0 % | 21 / 24 |
| si32 | 8 | 2884 -> 2272 (-21.2 %) | 3687 -> 3021 (-18.1 %) | 17955 -> 13716 (-23.6 %) | -3.6 % | 0 % | 21 / 24 |
| si64 | 4 | 4353 -> 3165 (-27.3 %) | 4730 -> 3770 (-20.3 %) | 22170 -> 16260 (-26.7 %) | -7.7 % | 0 % | 15 / 15 |

SCF and Davidson counts identical in every cell; |dE| 0 to 5.7e-14 Ha. Declined starts are solves whose `initial` was not the kept X. Wall: observation only (section 5 of the request); the walls of the one-pair cells overlapped a concurrent build and are void. Clean reruns (O): si32 np=4, 3 pairs at load 9: medians 6.71 -> 4.75 s (spread 2.7 s); si64 np=4, 3 pairs at load 7.6 to 8.7: 15.32 -> 14.37 s (-6.2 %, 3 of 3 pairs faster, spread 1.9 s), apply stage 5.08 -> 3.69 s per rank.

## 10. Tests
- dft-solvers `distributed_fft::tests::kept_planes_batch_equals_plain_batch_bitwise` (1 to 5 ranks, every kept count, bits).
- dft-physics `potential_change_from_kept_planes_matches_plain_apply` (two cells, 5 and 6 columns, 1 to 4 ranks: kept and unkept changes bit-identical, H_old X + dV X against H_new X at most 2.4e-16 relative, gate 1e-13; density bits unchanged by keeping).
- dft-solvers `tests/distributed_kept_start.rs` (1 to 4 ranks, selected roots on and off: taken, exactly `bands` H columns saved, eigenvalues within 3.3e-15 of the plain start, ranks agree, `x_local` equals the returned rows bit for bit, a 1e-10 perturbed `initial` declined with the plain H count).
- dft-engine `distributed_scf::m2_kept_start_keeps_the_standard_class_and_saves_columns` (cubic and skewed, 2 and 3 ranks: dE against threads-only and against off at most 2.1e-14 Ha, SCF and Davidson counts identical, H columns -34 / -32 percent, ranks bit-identical, fixed rank count bit-reproducible in energy and density).
