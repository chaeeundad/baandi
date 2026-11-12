# Work order P48: BLAS-3 and size scaling without any accuracy trade

Written 2026-11-13 for the performance owner and the implementers. Tickets: `tickets/P48A-SCALING-PROFILE.md` to `tickets/P48F-THREAD-PHASES.md`. Read `docs/performance/AGENT-BRIEF.md`, the ledger and the backlog first (AGENTS.md, Performance owner). Labels: M = measured on dev08 and relayed, O = Mac observation, D = arithmetic on recorded numbers, E = estimate.

## 1. Why: the gap grows with size because of the prefactor of the O(N^2) and O(N^3) terms

Source: `qe-vs-baandi-2026-11-11.md` (M, 208b226, one core, unperturbed).

| | si8 -> si32 wall | fitted exponent (D) | fit T = c1 x + c2 x^2 + c3 x^3, x = N/8 (D, three points only) |
|---|---|---|---|
| QE standard | 0.46 -> 2.99 s (6.5x) | N^1.35 | c2 0.056, c3 0.008 |
| BAANDI standard | 0.52 -> 7.25 s (13.9x) | N^1.90 | c2 0.229, c3 0.040 |
| BAANDI strict | 0.99 -> 13.3 s (13.5x) | N^1.87 | c2 0.257, c3 0.104 |

BAANDI's quadratic and cubic coefficients are 4 to 5 times QE's while the linear terms are at or below QE's, so small cells hide the gap. No operation of BAANDI has a higher asymptotic order than QE's (checked operation by operation, 2026-11-13); the difference is the constant of the same orders. Three code facts explain most of it:

1. **Pure-Rust GEMM runs in tiny tiles.** `dft-numerics/src/real.rs::bounded_matmul` calls faer `matmul` on 16x16x16 (Gram) or 4x4x256 (rotation) sub-blocks to stay in faer's allocation-free nano-GEMM branch. No packing, no cache blocking: the legacy rotation ran at about 17 GFLOP/s against 85 for MKL at the same shapes (stage-decomposition-results.md section 3, M). Every BLAS-3 site of the pure build (Gram, rotation, reduced matrix, real-path nonlocal products) inherits this.
2. **The complex (k-point) nonlocal operator is BLAS-1 per band.** `dft-physics/src/hamiltonian.rs::NonlocalOperator::add` forms projections by dot products and back-applies by axpys per projector channel, one band at a time (`apply_block` loops over columns). The real Gamma path already uses GEMM (`real_gamma.rs::RealNonlocal`, MKL dgemm default since NONLOCAL-DGEMM-AB, -6.4/-4.2 percent at T=1, M).
3. **Narrow Davidson steps are BLAS-2.** On si32 strict 66.4 percent of Davidson apply calls carry one correction column (ledger MPI-M1-CEILING, M); the projection, reduced-matrix border and in-block MGS of such a step are matrix-vector work bound by memory bandwidth, and the paired real FFT computes and discards a second column.

At eight cores a fourth fact multiplies the gap: the Amdahl serial fraction derived from the 8-core table is 50 percent for BAANDI si32 against 13 percent for QE (D: S8 = 1.77 vs 4.15). Memory-bound BLAS-2 work does not scale with cores, so items 1 and 3 also cap the thread speed-up.

BLAS-3 does not change the asymptotic order; it raises the flop rate of the same flops (arithmetic intensity) and lets that work scale with cores. Its effect on the measured exponent is that the crossover to N^3 dominance moves to larger N.

## 2. Accuracy contract (common to every P48 ticket; no exception)

Nothing in P48 may change what is computed, only how. Specifically:

- Unchanged: every stopping predicate and tolerance (strict residual, standard estimator, eigen schedule, buffer floors), the fixed-point certificate and its gap check, band counts, cutoff and grids, the pseudopotential treatment, mixing parameters, Davidson capacity.
- A change is either **bit-identical** (overhead only) or **rounding class** (different summation order or blocking). Rounding-class gates, from AGENTS.md and the 2026-11-11 decisions: strict class |E - E_ref| <= 1e-10 Ha against the strict certified reference of the same coordinates; standard class <= 5e-9 Ha; every strict run certifies; the 1e-12 Ha cross-path check on unique-solution fixtures is run and reported; density SHA-256 identical across T = 1, 2, 4, 8 within a backend (pure and MKL separately); SCF and Davidson counts reported, no new convergence failure on any fixture.
- Orthonormality of every Davidson basis stays inside the existing defect gate; rank decisions (`DEPENDENCE_TOLERANCE`) keep their per-column meaning.
- No physical hard-coding (adoption rule (a)): block widths, panel sizes and thresholds come from data size and shape rules, never from a material.
- Generality set for every adoption: `bench_suite.py` default set plus si16, si32, si64 Gamma (strict and standard, seeds 601-610 at +-1e-7 Bohr and fractional +-0.001), k444s, mgo-40, syn-p3, si-skew, syn-slab, si-elongated.
- Hot-path rules (`docs/performance/README.md`): no allocation, logging or I/O in numerical loops; `unsafe` only in FFI adapters.

**Excluded from P48 because each one changes the physics or the delivered accuracy:** real-space nonlocal projectors (aliasing of non-band-limited projectors), fp32 or mixed precision (P47B), radial projector interpolation tables (backlog H-6), looser eigen tolerances or estimator changes, smaller Davidson capacity (P31A), locking (P24A).

## 3. Tickets, order and dependencies

| Ticket | Content | Depends on | Expected effect (E unless marked) |
|---|---|---|---|
| P48A | Post-adoption scaling profile si8 to si64, T=1 and 8, both backends; per-site GFLOP/s and BLAS level; width histogram; serial fraction per phase | none (start now) | decides the order of B to F; adds the si64 case (file added with this work order) |
| P48B | Packed BLAS-3 GEMM for the pure build (replace the nano tiles) | A (baseline numbers) | pure build: rotation, Gram, reduced and nonlocal sites were about 50 percent of pure si32 strict wall before the adoption round (M, stage decomposition); 25 to 35 percent of pure wall if per-site GFLOP/s reaches 60 percent of MKL |
| P48C | Complex k-path nonlocal as GEMM over band blocks | A | k-mesh cases (k444s, mgo-40): the real-path analogue gave 4 to 6 percent at T=1 (M) starting from a faster kernel; 5 to 15 percent |
| P48D | Blocked orthonormalization inside the correction block (panel MGS + GEMM trailing update) | A shows the in-block MGS share; B for the pure kernel | up to the measured share of the S8 driver bucket (15 to 24 percent of wall before DGKS, M; current share unknown) |
| P48E | Width floor for narrow Davidson steps (generalizes SECOND-VECTOR v4) | v4 dev08 result; B and D change the cost of an extra column | si32 strict: v3 gave -9.3 percent wall (M) and the gain is the narrow-step share |
| P48F | Thread the phases that P48A shows serial at T=8; shape-based parallel thresholds | A | T=8 only; target serial fraction at si32 below 25 percent |

Start P48A immediately. B and C are independent of each other and can run in parallel worktrees after A's baseline exists. D follows B for the pure build (MKL build can start after A). E waits for the v4 verdict. F follows A.

## 4. Decisions for the user (one recommendation each)

1. **si64 QE reference on dev08** (one QE run per class at 1 and 8 ranks, about 15 minutes of node time): recommend yes; without it the scaling exponent rests on three points.
2. **`BAANDI_ENERGY_IDENTITY` default ON** (exact identity for Rayleigh-quotient eigenvalues; the certificate solve and the final state keep the explicit audit; ledger ENERGY-IDENTITY: final energy and density bit-identical on si16, final difference at most 1.5e-12 Ha on si32, O): recommend yes after a dev08 A/B; it is not an accuracy trade, but the ledger marks the change of history energies as a semantic change that needs the user.

## 5. dev08 requests (through the master session; protocol of AGENT-BRIEF.md)

Each ticket writes its request into its report. Common items: three-point comparison (baseline, main, branch) with binary SHA-256; si16, si32, si64 Gamma strict and standard, k444s and mgo-40; T=1 core 0 and T=8 socket 0; seeds 601-610; paired same-binary env A/B where the change is a switch; perf profiles with at least 30 samples on the changed function; per-site GFLOP/s from the GEMM accounting build; energy, certificate and density-SHA tables.
