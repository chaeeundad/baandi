# MPI entry cost: why the M2 path at np=2 is slower than the serial path

Work order: docs/performance/requests/work-order-beat-qe-2026-11-14.md, phase E2 (MPI against MPI only; no threads number appears here).
Evidence tags: (R) relayed from dev08, (O) Mac observation (pure build, simulated ranks = threads of one process, not dev08 evidence), (C) counted from shapes and code, (E) estimate.
Status 2026-11-14: analysis complete; fixes 1 to 4 below implemented behind switches that are OFF by default; dev08 measurement pending (request at the end).

## 1. The starting numbers (R, ledger MPI-ENTRY-COST-DIAGNOSIS)
si64 strict f1e-3 seed 601, MKL, one thread per rank. np=1 runs the serial threads-only path (no `parallel` block), np>=2 the distributed M2 path (`parallel.model = mpi_stick_rows`), so np=1 against np=2 compares two different code paths, not one path at two rank counts.

| | np=1 | np=2 | np=4 | np=8 |
|---|---:|---:|---:|---:|
| total_s | 109.3 | 160.7 | 108.4 | 73.7 |
| diagonalization | 80.2 | 128.5 | 81.3 | 43.7 |
| op.apply_fft_transposes | (serial h_apply 28 s in all) | 13.60 | 7.41 | 4.62 |
| op.apply_nonlocal | | 9.85 | 5.74 | 3.46 |
| op.apply_pack_unpack_kinetic | | 1.83 | 1.10 | 0.89 |
| collectives (all) | 0 | 3.37 | 4.62 | 4.57 |

Work counts are identical at every rank count (Davidson 97, H 9008, certificate H 1096).

Reading that changes the diagnosis: the apply stages explain only part of the np=2 penalty. At np=2 the apply is 25.3 s against an ideal 14 s (half of the serial 28 s), an excess of about 11 s. The rest of the diagonalization is 128.5 - 25.3 = 103 s against a serial 80.2 - 28 = 52 s for the whole algebra (ideal at np=2: about 26 s plus the replicated reduced eigenproblem). The distributed Davidson ALGEBRA (rotation, orthonormalization, residuals) carries an excess of roughly 70 s, six times the apply excess. The dev08 split of that 103 s exists in results.json `parallel.stages` (reduced_gram, reduced_eigen, rotate, residual, orthonormalize, apply_total) but was not relayed; it is the first item of the dev08 request.

## 2. Mac decomposition at equal kernels (O)
`crates/dft-engine/tests/m2_entry_cost.rs` (ignored test, release) runs the same Si SCF (Si.pz-vbc, 15 Ha, f1e-3 deterministic displacements, strict defaults) on the threads-only path at T=1 and on the M2 path with simulated ranks, and prints the serial Davidson stages (`BAANDI_DAVIDSON_STAGES=on`), the distributed stages, the M2 operator stages and the collectives, per rank. The pure build uses the same FFT kernels (rustfft lanes) and the same faer GEMM on both paths, so differences are structural, not kernel differences.

si32 (2x2x1 conventional; SCF 16, Davidson 279, H 6732 on every path), seconds per rank:

| stage | serial T=1 | M2 np=2 | M2 np=2 per unit work (x serial) |
|---|---:|---:|---:|
| wall | 19.32 | 18.33 | |
| h_apply (local 4.39 + nonlocal 1.73) | 6.09 | 3.47 | 1.14 |
| of which FFT incl. transposes / nonlocal | 4.39 / 1.73 | 2.35 / 0.95 | 1.07 / 1.10 |
| rotate | 4.91 | 3.56 | 1.45 |
| orthonormalize | 3.44 | 6.26 | 3.64 |
| residual | 0.06 | 0.54 | about 17 (includes wait in its reduction) |
| reduced Gram + reduced eigen | 1.58 + 1.10 | 0.99 + 1.19 | eigen replicated |
| tree_sum calls / s | 0 | 39363 / 0.64 | |

On equal kernels the distributed apply is at its ideal (1.07 to 1.14 per unit work); the orthonormalization costs 3.6 times its serial counterpart per unit work and is the largest single stage at np=2.

## 3. Root causes (C unless marked)

### 3.1 Orthonormalization: the classic distributed order does about four times the serial work, half of it as scalar loops
`distributed_davidson::orthonormalize` (M1 design, reused by M2): normalise, project twice against the basis, renormalise, `mgs_panel` over the whole block, project a third time, `mgs_panel` again. `mgs_panel` is a two-pass classical Gram-Schmidt of every column against ALL earlier kept columns of the block with three reductions per column (two passes of dot products, one norm) and per-column allocations. The serial solver (`Projection::run_into`) projects once, repeats only when the DGKS rule fails (`BAANDI_ORTHO=dgks:0.1`, default since MICRO-ORTHO), and runs the in-block Gram-Schmidt in panels (MGS-B, 16 columns, on MKL builds; ORTHO-BLOCK shape rule on pure builds) with the cross-panel part as GEMM.
- Scalar in-block work per rank: classic 8 nl sum(w^2) flops; serial order about 60 nl sum(w) in-panel plus 4 nl sum(w^2) as GEMM. si64 np=2 (nl about 12000, 97 iterations, mean correction width about 77 from (9008 - 11 x 140) / 97, sum(w^2) about 7.5e5): classic 72 Gflop of streaming dot/axpy on a 13 MB working set (memory bound, 1 to 3 Gflop/s on one Xeon core: 25 to 70 s, E); serial order 5 Gflop scalar in L2-sized panels plus 36 Gflop GEMM.
- Basis projections: 4 nl m w flops each (m the basis width, about 250 on average, E); classic 3 per orthonormalization, serial order 1 when DGKS holds: 320 against 107 Gflop per rank at si64 np=2.
- Reductions: classic about 6 w + 5 tree_sums per Davidson iteration; this is the source of the 76258 tree_sum calls (97 iterations x about 786, i.e. w near 130 in the wide early iterations). The calls cost only 1.23 s at np=2 (16 us each); the problem is the O(nl w^2) scalar work they accompany, not their latency.

### 3.2 Distributed rotations and projections fall off dgemm on MKL builds
`distributed::rotate_local` and the rotation inside `project_out` call `multiply_slices_into`, whose dgemm route requires the whole output (rows x cols) to fit the square workspace bound capacity^2. Rotation outputs are nl x bands (si64 np=2: 12000 x 140 = 1.7e6 against 280^2 = 7.8e4), so on MKL builds every rotation and projection-apply runs the pure-Rust nano kernel. The serial solver uses `rotate_pair_slices` and (P49A) the wide route, both per-panel dgemm. Measured rates (R, NONLOCAL-DGEMM-AB): dgemm 27 GFLOP/s against the pure tiled kernel 9.7 GFLOP/s on one dev08 core. Rotation work per rank at si64 np=2: 2 x 2 nl m bands per iteration, about 163 Gflop: 16.8 s pure against 6.0 s dgemm (E).

### 3.3 Distributed nonlocal GEMMs never use dgemm
`DistributedRealHamiltonian::projections_reduced` and the back product call `RealGemmWorkspace::multiply_into` (pure tiled kernel) and the module documents that `BAANDI_MKL_NONLOCAL` is not honoured; the serial `RealNonlocal` takes both products through `multiply_slices_wide_into` (per-panel dgemm, default ON on MKL since ADOPT-DEFAULTS-1). Same 2.8x rate gap: op.apply_nonlocal 9.85 s at np=2 would be about 3.5 s (E). The tree_sum of P^T X (count x m) is one call per apply and stays.

### 3.4 Distributed FFT uses rustfft lanes; the serial MKL path uses the full-grid MKL 3D FFT
`DistributedPairFft` is the mirror-slab pair plan of `PairFusedPlan<RustLanes>` (rustfft 1D batches along the stick, middle and slab buffers) with one alltoall per direction for the whole block. On MKL builds the serial path uses `MklPairPlan::Grid` (full-grid MKL 3D FFT, `BAANDI_PAIR_FFT=grid`), not the pruned pair plan. The pure-build comparison (section 2) shows the distributed layout itself at its ideal; on dev08 the kernel is the difference: PFFT-2 (O) puts 64 to 66 percent of the pair apply in the 1D kernels, and PFFT-5 (E) gives about 213 / 487 us per transform for the rustfft mirror plan against about 97 / 235 us for an MKL 1D stick pipeline at 36^3 / 72x36x36. Transposes: 2 alltoalls per apply (plus 2 mirror exchanges); their cost is inside op.apply_fft_transposes (alltoall total 1.97 s at np=2 including mirrors). Packing into the f64 send buffers is per element (scalar re/im copies with index arithmetic); no per-call allocation (buffers grow to the widest batch once).

### 3.5 Smaller items
- Residual and correction blocks use `RealMatrix` element indexing (bounds-checked, column-major index per element) and allocate per iteration; serial computes them fused. Mac residual stage 0.54 s per rank at si32 np=2 (part of it is waiting in its norm reduction).
- Per-iteration allocations and copies: `append_columns` copies the whole basis and products every iteration (2 nl m doubles; at si64 np=2 about 48 MB per iteration, about 1 s per run, E); `apply` allocates its output block per call.
- The distributed solver lacks the serial solver's algorithmic stack (B1 selected roots, B5 bordered reduced eigen update, B6 projector-product reuse in the certificate, RITZ_FUSED): the reduced eigenproblem is replicated and solved from scratch every iteration (1.19 s per rank at si32 np=2 against 1.10 serial; it does not shrink with ranks).
- Collectives do not shrink with ranks (3.4 to 4.6 s): tree reductions grow with log2 P and the alltoall block size is set by the largest (rank, plane range) pair.

## 4. Fixes, ranked by expected seconds at si64 strict (E, per rank, from the counts above)

| rank | fix (switch) | touches | np=2 expected | np=8 expected | status |
|---|---|---|---:|---:|---|
| 1 | serial-order orthonormalization: one projection per pass, DGKS second-pass rule, panel Gram-Schmidt with the serial panel rule (`BAANDI_M2_ORTHO=serial`) | distributed_davidson.rs | -20 to -40 s | -5 to -10 s | implemented, off |
| 2 | rotations and basis projections through the wide per-panel dgemm route (`BAANDI_M2_DGEMM=on`) | distributed.rs, distributed_davidson.rs | about -11 s (rotation) plus the projection applies | about -3 s | implemented, off |
| 3 | distributed nonlocal on the dgemm route (`BAANDI_M2_DGEMM=on`, same switch) | distributed_real_gamma.rs | about -6 s | about -2 s | implemented, off |
| 4 | MKL DFTI 1D batch lanes in the distributed FFT (`BAANDI_M2_FFT=mkl`) | distributed_fft.rs, fused_pair_plan.rs (`BackendLanes`) | about -6 s | about -2 s | implemented, off |
| 5 | slice loops for residual and correction blocks (`BAANDI_M2_LEAN=on`, bit-identical) | distributed_davidson.rs | under -1 s | small | implemented, off |
| 6 | fuse the in-panel norm with the second-pass dot products (Pythagorean norm, 3 to 2 reductions per column) | distributed_davidson.rs | -0.3 s | -0.5 s | not done (rounding at the dependence floor needs its own test) |
| 7 | serial algorithm stack in the distributed solver (B5 bordered eigen update, B6 certificate reuse, B1) and reuse of basis/product storage instead of `append_columns` | distributed_davidson.rs | several s | grows in relative terms | not done (E2 step 2, after Phase B) |
| 8 | replicated stages (Hartree, mixing, reduced eigenproblem) | m2.rs, scf | | the np>=8 ceiling (MPI-SCALING-BEYOND-8) | not done; mixing was deliberately replicated in M2d |

Sum of 1 to 5 at np=2: about -45 to -65 s of 160.7 s, which brings np=2 near the np=1 serial path (109.3 s); at np=8 about -12 to -17 s of 73.7 s (E). Fixes 1 and 2 address the 70 s algebra excess of section 1, fixes 3 and 4 the 11 s apply excess.

Mac result of fixes 1, 2 and 5 (O, pure build: fix 2 is bit-identical there, so the gain is fix 1 plus 5), wall seconds, work counts identical on every row:

| case | default | M2_ORTHO=serial, M2_LEAN, M2_DGEMM | orthonormalize stage per rank | tree_sum calls per rank |
|---|---:|---:|---:|---:|
| si32 np=2 | 18.33 | 13.74 (-25 %) | 6.26 -> 2.10 | 39363 -> 23507 |
| si32 np=4 | 12.89 | 10.22 (-21 %) | 4.11 -> 1.57 | 39363 -> 23507 |
| si16 np=2 | 2.48 | 1.90 (-23 %) | 0.70 -> 0.26 | 17055 -> 10398 |
| si16 np=4 | 1.97 | 1.70 (-14 %) | 0.58 -> 0.30 | 17055 -> 10398 |

si32 np=2 with `BAANDI_M2_ORTHO=serial` alone: 13.96 s (the orthonormalization order is the whole Mac gain). |dE| against the serial T=1 run 2.8e-14 to 8.5e-14 Ha in every row. Single runs on a shared Mac: observations, not acceptance.

## 5. Ledger cross-check (nothing rejected is repeated)
- CholeskyQR2 (A52) and the M1 Cholesky panel that dropped near-dependent corrections (MPI-M1-GATE1) are not used: fix 1 keeps the per-column relative floor (`DEPENDENCE_TOLERANCE`) and two in-panel passes, exactly the serial MGS-B/ORTHO-BLOCK structure (adopted: MGS-B R2 on MKL, P48D-ORTHO-BLOCK on pure) and the adopted DGKS rule (MICRO-ORTHO).
- Fix 2 is the distributed counterpart of P49A-PROJECT-WIDE (serial, opt-in, dev08 pending) and of the adopted MKL rotation dgemm; fix 3 of NONLOCAL-DGEMM-AB/ADOPT-DEFAULTS-1; fix 4 reuses the HAPPLY-PAIR-MKL lanes (serial opt-in, no dev08 number yet). None was rejected.
- M2d decisions are untouched: batched transposes (`BAANDI_M2_BATCH`), replicated mixing and Hartree, distributed density/XC/energy.
- FFT lane batching (A11) and band-grouped 3D FFTs (A28) are not revived: fix 4 changes only the 1D kernel library inside the existing lane buffers.

## 6. Determinism and accuracy contract (ADR 0041)
- All switches off: bit-identical to the current M2 code (no arithmetic change; the FFT keeps the same rustfft plans through `BackendLanes::Rust`).
- `M2_LEAN`: bit-identical (same element expressions). `M2_DGEMM` on pure builds: bit-identical (both routes are the same product); on MKL builds a rounding change with fixed 256-row panels, so results do not depend on the thread count. `M2_FFT=mkl`: rounding class on MKL builds, one DFTI descriptor per lane buffer; pure builds have no vendor lanes and stay bit-identical. `M2_ORTHO=serial`: rounding class (different summation order, same rank decisions rule).
- Every branch depends on reduced (replicated) values only, so all ranks take the same branches and enter the same collectives; a fixed rank count is bit-reproducible; one rank never enters this code (threads-only path).

## 7. Tests (Mac, pure build)
- `crates/dft-solvers/src/distributed_fft.rs` `vendor_lanes_match_rust_lanes`: vendor lanes within 1e-14 of rustfft lanes for ranks 1, 2, 3, 5 on cubic, elongated and skewed grids; bit-identical on pure builds.
- `crates/dft-physics/tests/distributed_real_gamma.rs` `nonlocal_dgemm_route_matches_the_default_route`: H x with the dgemm route within 1e-14 relative of the default (bit-identical on pure), ranks 1, 2, 3, 5, m = 1, 3, 16, three cell shapes.
- `crates/dft-solvers/tests/distributed_davidson.rs` `lean_and_wide_tuning_keep_the_classic_bits` (bit identity of eigenvectors and work counts, ranks 2 to 5) and `serial_order_orthonormalization_matches_serial_eigenpairs` (panel rules whole block, 4 and shape rule, DGKS on and off: eigenvalues within 1e-9 of the serial solver, residuals below 1e-8, ranks bit-equal, reproducible per rank count, iterations within 10 percent plus 2 of the classic order; scattered row partition and clustered spectra).
- `crates/dft-engine/tests/distributed_scf.rs` `m2_tuned_solver_matches_threads_only`: SCF with all Davidson switches forced per simulated rank (`set_tuning_override`), cubic and skewed cells, strict and standard, 2 and 3 ranks: |dE| 0 to 2.2e-16 Ha strict, at most 6.6e-11 Ha standard; ranks bit-equal; reproducible.
- `crates/dft-numerics/src/switches.rs` `mpi_entry_cost_switches_default_off`; `crates/dft-cli/tests/run.rs` switch count 52 -> 56.
- Profile harness: `DFT_SI_UPF=<Si.pz-vbc.UPF> BAANDI_DAVIDSON_STAGES=on M2_PROFILE_CASE=si32 M2_PROFILE_RANKS=2,4 cargo test --release -p dft-engine --test m2_entry_cost -- --ignored --nocapture`.

## 8. dev08 request (MPI against MPI only)
Build: the branch HEAD with `--features mkl,mpi` (MKLROOT=/opt/intel/oneapi/mkl/2026.1), one binary, sha256 recorded; MKL_NUM_THREADS=1, one thread per rank, `--map-by core --bind-to core` on socket 0, load gated as in BAANDI-MPI-RANK-SWEEP.
Cells: si32 and si64, strict and standard, f1e-3 seed 601 (add 602, 603 if time allows), ranks 2, 4, 8 (np=1 is the serial path and is not affected; report it once as the reference of the same binary).
Arms (same binary, environment A/B, interleaved per seed):
1. all off (current default),
2. `BAANDI_M2_ORTHO=serial`,
3. `BAANDI_M2_DGEMM=on`,
4. `BAANDI_M2_FFT=mkl`,
5. `BAANDI_M2_LEAN=on`,
6. all four together.
Report per row: wall, results.json `parallel.stages` (reduced_gram, reduced_eigen, rotate, residual, orthonormalize, apply_total, apply_operator), `parallel.m2.stages` (op.* and engine stages), `parallel.collectives` (calls and seconds), SCF / Davidson / H counts, total energy (|dE| against arm 1 at the same rank count: strict 1e-10 Ha, standard 5e-9 Ha), density SHA on two repeats of arm 6 at np=4 (bit-reproducibility at fixed rank count) and every rank's energy equal to rank 0. Also run `cargo test --features mkl -p dft-solvers --lib vendor_lanes_match_rust_lanes`, `cargo test --features mkl -p dft-physics --test distributed_real_gamma nonlocal_dgemm_route` and `cargo test --features mkl -p dft-numerics --lib mkl_lanes_match_rust_lanes` (MKL rounding checks that the Mac cannot run).
Adoption (AGENTS.md adoption rule): each switch whose arm improves the median wall at np=2, 4 and 8 without a correctness failure becomes the default with its variable kept for rollback.
