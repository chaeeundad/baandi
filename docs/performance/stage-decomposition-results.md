# Stage decomposition results (2026-11-09): where the QE gap is

Data: [dev08_stage_decomposition_relayed.md](../../tests/reference/perf/dev08_stage_decomposition_relayed.md) (A QE timers, B/C BAANDI stage profile and counters at ee86f07, D ceilings). Plan and stage list: [stage-decomposition.md](stage-decomposition.md). One core, real Gamma, pure native build (the MKL build moves algebra into an unattributable BLAS bucket). Seconds of BAANDI per stage = call-site share x non-sampled wall; QE stages from its timers (FFT = vloc_psi, nonlocal = calbec + add_vuspsi, Davidson algebra = regterg total - h_psi, other = wall - the three). Observation: the grouping (FFT includes the pair transform of the density, S8 driver and H dispatch are put into algebra) is the orchestrator's, not the master's.

## 1. Gap by stage group (BAANDI / QE seconds, ratio, share of the total gap)
| case | profile | wall B / Q | FFT + local | nonlocal | Davidson algebra | other (setup, density, potential, mixing, audit) |
|---|---|---|---|---|---|---|
| si8 | standard | 1.10 / 0.47 = 2.4x | 0.26 / 0.16 = 1.6x (15% of gap) | 0.12 / 0.01 = 12x (18%) | 0.35 / 0.05 = 7x (48%) | 0.37 / 0.25 (19%) |
| si8 | strict | 2.33 / 0.55 = 4.2x | 0.41 / 0.21 = 2.0x (11%) | 0.11 / 0.01 = 11x (6%) | 1.52 / 0.06 = 25x (82%) | 0.28 / 0.27 (0%) |
| si16 | standard | 3.06 / 1.41 = 2.2x | 0.58 / 0.66 = 0.9x (-5%) | 0.47 / 0.08 = 5.9x (24%) | 1.41 / 0.24 = 5.9x (71%) | 0.60 / 0.43 (10%) |
| si16 | strict | 9.73 / 1.83 = 5.3x | 1.73 / 0.87 = 2.0x (11%) | 0.82 / 0.12 = 6.8x (9%) | 6.46 / 0.34 = 19x (77%) | 0.72 / 0.50 (3%) |
| si32 | standard | 16.18 / 5.04 = 3.2x | 2.93 / 2.48 = 1.2x (4%) | 3.46 / 0.45 = 7.7x (27%) | 8.06 / 0.98 = 8.2x (64%) | 1.73 / 1.13 (5%) |
| si32 | strict | 49.53 / 6.70 = 7.4x | 5.15 / 3.25 = 1.6x (4%) | 3.96 / 0.63 = 6.3x (8%) | 38.24 / 1.46 = 26x (86%) | 2.18 / 1.36 (2%) |
QE comparison class: conv_thr 1e-8 for standard, 1e-10 for strict (strict BAANDI certifies, QE does not).

## 2. What the table says
1. **Davidson algebra is 64 to 86 percent of the whole gap** at every size and in both classes (standard 5.9 to 8.2x its QE counterpart, strict 19 to 26x). It is the only group whose ratio grows with the size in strict.
2. **FFT is near parity** (0.9 to 2.0x; the pair transform itself is faster than QE per call, section D4) and is 4 to 15 percent of the gap, mostly the S2 share in the strict runs.
3. **Nonlocal is 6 to 12 times QE's calbec + add_vuspsi** and 6 to 27 percent of the gap, the largest group after algebra, and the largest single share of the standard class at si32 (21 percent of the wall). QE does it with BLAS; the MKL dgemm ceiling at those shapes is 80 to 93 GFLOP/s (D1).
4. Setup, density, potential, mixing, energy audit together are 0 to 19 percent of the gap, within 1.1 to 1.6x of QE's; the phase timers show total_energy at 3.0 to 3.7 s on si32 (6 percent of the wall), which the sampling does not see (master's limit 2): worth a look, not a priority.
5. Inside the algebra of the pure build (si32 strict, share of wall): rotate 35.0, gram/MGS products 15.2, **driver 21.4**, eigen 3.5, reduced 1.0, residual 0.7. After MKL-R1/R2 (MKL build, ee86f07): rotate 15.5, MKL BLAS 15.9, **driver 24.2**, gram 3.3, eigen 3.6, reduced 1.8. The driver bucket is the self time of `real_davidson::run_into`: the orthonormalization of the correction block against the basis and within itself is done with scalar column dot products and axpys (two passes, O(n k^2) per step, no GEMM). It does not shrink with MKL and is now the largest single item of the MKL strict profile.

## 3. Optimization points, ranked by (share of gap) x (plausible reduction), with the ledger check
| # | Point | Why | Expected gain (to be measured on dev08, distributions) | Ledger/rejection check |
|---|---|---|---|---|
| 1 | **Blocked orthonormalization of the correction block** (replace scalar two-pass MGS in `Projection::run_into` by a GEMM-based scheme: block classical Gram-Schmidt against the basis already uses GEMM; within the block use Cholesky-QR with a second pass, or blocked MGS with panel GEMM updates, and keep the orthogonality defect gate) | S8 driver is 15 to 24 percent of the wall in every configuration, does not respond to MKL, and is serial scalar code at about 1 to 2 GFLOP/s | up to 10 to 20 percent of the wall in strict and 8 to 20 percent in standard (the share itself) | Not tried before: P37A generalized (non-orthonormal) Davidson and P24A locking changed the algorithm; this keeps an orthonormal basis and only changes how the block is orthonormalized (rounding-level trajectory change: judge by perturbed-start distributions as for MKL-R1) |
| 2 | **Nonlocal operator through GEMM** (projections P^T X and back-application as dgemm over band blocks, real path `real_gamma::add_block`) | 6 to 12x QE, 21 percent of the standard si32 wall, ceiling 80 to 93 GFLOP/s | up to 10 to 18 percent of the wall in standard, 6 to 10 percent in strict | Ledger A43 (P10J cache-blocked traversal) and P15B (band-blocked nonlocal) were slower on dev08; the new premise is the measured MKL dgemm ceiling at the exact shapes (D1) and that the current add/add_block code is a scalar loop, so this must be an MKL-feature path (like MKL-R1/R2) or a real SIMD kernel, not another traversal order |
| 3 | **Pure-Rust build GEMM** (rotate and gram are 50 percent of pure si32 strict; legacy rotate ran at 17 GFLOP/s against 85 for MKL; the auto-vectorised blocked kernel ROT1 was slower) | the default non-MKL build keeps this cost | up to 25 to 35 percent on non-MKL builds | ROT1 rejected; needs explicit SIMD micro-kernels (pulp is already in the dependency tree through faer; policy check required) or shipping MKL as the supported x86 backend |
| 4 | Reduce Davidson steps and H vectors (QE uses 0.1x the vectors per SCF step and a Ritz-change stop; BAANDI standard already does) | step-count part of the algebra ratio | unknown | locking, generalized, band-count rules, preconditioners, mixing all rejected (ledger P24A, P37A, P38A, P40A, P42A, P43A, P47A) |
| 5 | Energy audit / total_energy phase (3 to 3.7 s on si32) and density (1.2 to 1.4 s on si32 vs QE sum_band 0.47) | small groups with a visible ratio | 3 to 7 percent on si32 | audit shared terms already done (P28B); the density path was not profiled at call-site level |

## 4. Next steps
1. Implement point 1 (blocked orthonormalization) in `dft-solvers`: pure-Rust blocked version first (works for every build), validate with `bench_scf --require-identical` for the untouched default if the old path is kept as the default until dev08 distributions pass the frozen rule.
2. Point 2 as an MKL-feature path after point 1, using the same same-binary on/off switch method.
3. Ask the master for the finer split of the driver bucket (confirm that the leaf `run_into` is the MGS loops, not copies) before coding point 1, and for call-site attribution of the MKL BLAS bucket by a wrapper build (the master's offer).
