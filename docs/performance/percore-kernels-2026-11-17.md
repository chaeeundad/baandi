# Per-core kernels of one Davidson solve (ticket PERCORE-KERNELS, 2026-11-17)

Question: one core BAANDI is 1.78 to 2.60 times slower than QE (dev08, 2026-11-17, new default stack; T=8 1.46 to 2.43; MPI 2.6 to 3.6, part of it this per-core cost). Which kernels of the diagonalization run far from the rate the build can reach, what does each one cost per H column against QE, and what is the cheapest fix for the worst three. Work at the diagonalization level (stage 0, [fast-loop-protocol.md](fast-loop-protocol.md)): one Davidson solve at a fixed Hamiltonian, counters first, then time.

Labels: **O** Mac observation (Apple M6 core, pure build: no MKL, every GEMM on the pure-Rust kernel; shared host, load 3 to 21 during this session; never acceptance evidence); **R** relayed dev08 number (source named); **E** estimate; **C** exact counter (machine-independent). Main e09067d, branch `worktree-agent-aa09292b818f370fe`.

## 1. Method
- Harness: `crates/dft-cli/examples/davidson_bench.rs`, `scripts/fast_screen_dav.py`. Inputs `fast_screen_dav.py --emit-input` (si32-gamma / si64-gamma, seed 601, fractional amplitude 1e-3, T=1), base density = converged `dft-cli run` density.bin of the same input. Cold = atomic start at 5e-3 Ha; warm = kept start (REUSE-PSI) at the ratchet tolerance (si32 3.59e-11 Ha), i.e. the typical mid-SCF solve.
- Stage seconds: `BAANDI_DAVIDSON_STAGES=on` (3 repeats; timer runs are not used for A/B timing). FLOP and calls: `gemm_dispatch` rows (C). Function shares: macOS `sample` 1 ms, 12 s of repeated warm maps, 6995 samples below `lean::solve` (O; includes the untimed re-arm solves of the harness).
- Per-iteration shapes (C): a temporary trace (not committed). si32 warm: (basis m after the append, block k) = (68, 68), (136, 64), (200, 10), (210, 7), (217, 7), (224, 7), (231, 6), (237, 5), (242, 2), return at m = 244; cold: (68, 68), (136, 54), return at m = 190. So a warm solve has two wide appends and seven narrow ones (k <= 10).
- Kernel microbenchmark at exactly these shapes: `cargo test --release -p dft-numerics --test percore_kernels -- --ignored --nocapture percore_kernel_microbench` (new, this ticket) and the existing `pure_gemm_microbench` (achievable rate of the build's library kernel, faer `matmul` with full packing = "library" column).

## 2. Per-kernel table, si32, one core (Mac pure build, O; counters C)

Cold solve 264 ms (3 Davidson iterations, 190 H columns); warm solve 387 ms (10 iterations, 176 H columns). Rates are 2mnk / stage seconds.

| kernel | stage | cold ms | warm ms | warm share | GFLOP (warm) | calls (warm) | shapes (rows x out x reduction) | achieved GFLOP/s (warm) | achievable on this build (O) | AI (flop/byte) |
|---|---|---|---|---|---|---|---|---|---|---|
| local pair FFT (H apply) | h_apply.local | 110.6 | 106.6 | 27.5 % | about 3.0 nominal (92 pair transforms x 32.7 MF, 5 N log2 N) | 92 pair transforms | grid 72 x 72 x 36, pruned mirror-slab plan | 28 nominal (16 counting only the pruned lanes) | rustfft kernel-only floor 0.84 of the pair time (PFFT-5, E) | 1 to 2 |
| nonlocal P^T X and P W (H apply) | h_apply.nonlocal | 41.7 | 48.2 | 12.4 % | 1.07 | 18 (9 applies x 2) | 11877 x 128 x k (TN) and 11877 x k x 128 (NN) | 22 | 46 to 58 (packed / library at these shapes) | 2 k / 8 |
| reduced border V^T (H E) | reduced | 21.9 | 40.5 | 10.5 % | 0.869 | 10 | m x k x 11877, m 68..244, k 68..2 (TN) | 21.5 | 38 to 46 (library TN) | k / 4 |
| overlap border V^T E (+ normalize, append copy) | orthonormalize | 21.9 | 38.3 | 9.9 % | 0.869 | 10 | as the reduced border (TN) | 22.7 | 38 to 46 | k / 4 |
| Ritz X_s = V C_s, H X_s (+ return rotation) | rotate | 39.9 | 78.3 | 20.2 % | 1.878 | 20 | 11877 x k x m (NN), narrow k 2..10 in 7 of 10 iterations | 24.0 (narrow: 4.5 to 10, microbenchmark) | 45 to 48 (library NN) | k / 4 |
| factor, T H_r T^T, coefficients, EVD | eigen | 4.1 | 26.7 | 6.9 % | 0.27 GEMM + sum m^3 = 95.9e6 (10 full EVDs) | | m <= 244 | | | |
| precondition + residual | | 1.8 | 2.5 | 0.6 % | | | | | | |

Sample shares below `lean::solve` (warm and re-arm maps, 6995 samples, O): H apply 54.4 % (FFT 2254, GEMM kernel 907, copy 160), reduced border 10.2 %, Ritz loop rotation 10.1 % (of which `memmove` 22 %), overlap border 8.5 %, return rotation 5.8 %, EVD 4.9 %, start projection 2.4 %, factor_fresh 1.9 %. Parents of every `memmove` sample: the pure GEMM kernel 297 (4.2 % of the solve: one variable-length copy per eight doubles when a strip is packed and when a tile column is stored), rustfft out-of-place butterflies 186, the pair plan's `inverse_slab` 162.

Achievable rates on the Mac (pure_gemm_microbench, T=1, GFLOP/s, O): per-panel packed kernel / faer library: rotate shapes 34 to 42 / 45 to 48, TN gram shapes 23 to 25 / 38 to 46, narrow TN (k = 4) 13 / 16.5, nonlocal 28 to 47 / 57 to 58. An Apple core's f64 NEON peak is about 64 to 72 GFLOP/s (E: four 128-bit FMA pipes).

### dev08, MKL build (R; no new dev08 numbers in this ticket)
- dgemm ceiling one core: 98.5 GFLOP/s square; production shapes standalone 76 to 101 (D1, tests/reference/perf/dev08_stage_decomposition_relayed.md).
- In process (DGEMM-SITES, strict, 2026-11-09, classic Davidson): NN rotate 58.3 GFLOP/s (70 % of its shape ceiling), TN overlap 27.7 (27 %), TN reduced 30.7 (34 %), projection apply 9.7 (12 %). The TN products with a 256-row reduction panel are the farthest from the ceiling (DGEMM-CENSUS: TN at k = 256 runs 42 to 59 standalone; efficiency follows the work per call; 25 % of the classic projection calls had one column).
- FFT: the mirror pair plan is at or above the library per transform (72 x 72 x 36: 647 us per band against MKL in place 727, FFTW3 667, QE fftw timer 762 per transform; D addendum). The local FFT is not where the per-core gap is.
- The lean default changed the shapes (borders instead of projections; narrow appends); the MKL per-kernel rates of the lean solve are not measured yet: the dev08 request asks for them (stage timers plus `BAANDI_DGEMM_TIMERS`).

### Per H column against QE (si32, one core)
| | QE 7.6 dev08 (R, E derivation) | BAANDI Mac warm (O) | BAANDI Mac cold (O) |
|---|---|---|---|
| diagonalization per H column | 1.87 ms (c_bands 6.26 s / about 3340 columns) | 2.20 ms | 1.39 ms |
| H apply per column | 1.29 ms (h_psi 4.31 s): local 1.05, nonlocal 0.23 | 0.89 ms: local 0.61, nonlocal 0.27 | 0.81 ms |
| Davidson algebra per column | 0.54 ms (overlap + update + last + rdiaghg 1.79 s) | 1.06 ms | 0.47 ms |
| algebra GEMM work per column | about 10 to 12 MFLOP (E: capacity 2 x bands, mean basis about 1.5 x bands) | 22 MFLOP (C: 3.89 GFLOP / 176) | 14.6 MFLOP (C) |

QE's column count is derived from its fftw calls (3914, minus sum_band's 15 x 34 inverse transforms and about 60 potential transforms, two per pair of bands; E). The machines differ, so only the structure is compared: the H apply per column is at QE's level; the algebra per column is about twice QE's because each column carries about twice the GEMM work (the basis is up to 6 x bands, mean about 3 x bands in the warm solve, against QE's 2 x bands; capacity is settled, ledger A13 / SERIAL-CAPACITY-SCREEN) and because the narrow products run at a small fraction of the GEMM rate.

## 3. Ranking: wall share x gap to the achievable rate (warm si32, Mac O; dev08 direction R)
| rank | kernel | share | achieved / achievable | excess share of the solve | dev08 direction |
|---|---|---|---|---|---|
| 1 | TN borders (overlap + reduced) | 20.4 % | 22 / 40 | about 9 % | worse: TN in process 28 to 31 of 76 to 101 |
| 2 | NN Ritz rotation, narrow blocks | 20.2 % | 24 / 45 (narrow 4.5 to 10) | about 9 % | NN 58 of 84 |
| 3 | nonlocal products | 12.4 % | 22 / 46 to 58 | about 6 % | MKL dgemm wide; per-call rate of the lean shapes unmeasured |
| 4 | local pair FFT | 27.5 % | at library rate (dev08) | about 4 % (kernel floor) | at or above MKL/FFTW/QE per transform |
| 5 | eigen (factor, EVD) | 6.9 % | faer EVD | small | MKL dsyevd |

## 4. Fixes implemented (switches default off) and stage-0 screens

### 4.1 Rank 1: fused border, `BAANDI_LEAN_BORDER=fused`
Mechanism. The lean solver formed the overlap border V^T E before the H apply (for the pivot test) and the reduced border V^T (H E) one iteration later: two transposed products over the same basis V (n x m), each with k output columns. Fused: H is applied to the unit corrections first, then one product V^T [E | H E] (k doubled, V streamed once) gives both borders; the pivot test follows. A failing pivot test (measured rate 0 on the Si cases) discards the images (counter `lean.fused_wasted_columns`) and takes the existing fallback, which applies H again to its own columns. Fresh bases (start, contraction) are unchanged. Cost: copies of E and H E side by side (about 4 n k of traffic against the saved n m stream of V), hence the shape rule: only blocks with 4 k <= m fuse (the first, wide append of every solve keeps the separate products). A failed pivot test under fused wastes the k H columns of the block; after the first one the solve continues separately (at most one block per solve). The synthetic fixtures of `tests/lean_davidson.rs` show why both rules are needed: without them the near-degenerate fixture discarded 439 of 1296 H columns (capacity 3) and 98 of 403 (capacity 6), the cubic one 20 of 64; with the per-solve rule 10 of 867 / 315 / 54; with the shape rule added only the near-degenerate fixture discards one narrow block (10 columns). On the Si stage-0 cases the pivot test never failed (`fused_wasted_columns` 0).
Identity. On the aarch64 pure kernel every output entry is one ascending FMA chain inside its fixed 256-row panel, independent of the output width, so the fused product is bit-identical (unit test `percore_kernels::fused_border_matches_two_products`, solver tests `lean_davidson::percore_*` assert the whole solve bit for bit, T = 1, 2, 4, 8). MKL dgemm and x86 faer may round differently with the width: rounding class there.
Why it should matter more on dev08: the TN products are the farthest from the ceiling there, and DGEMM-SITES found the dgemm efficiency growing with the work per call; the fused call has twice the columns and half the call count.

Kernel microbenchmark (O, two TN products of width k against one of width 2k, copies included; two runs, load 3 to 7): narrow k = 1 to 12: 0.49 to 0.88 of the separate time (si32 k = 10: 0.73 to 0.79; k = 2: 0.51; si64 k = 6: 0.62 to 0.65); wide k = 64 to 132: 0.93 to 1.17 (copy cost, within noise). Rates: si32 k = 2 TN 11 -> 22 GFLOP/s.

Stage 0 (O, si32, 3 pairs x 5 repeats, pure, T=1): cold ratio 1.004 (spread 1.0 %) and 0.978 (4.0 %) in two runs, NEUTRAL (two wide appends only); warm 0.964 (1.75 %) **CANDIDATE** and 0.949 (5.2 %, NEUTRAL by the spread); eigenvalues bit-identical; counters: reduced-site GEMM calls 10 -> 1 (warm), GEMM FLOP unchanged, H columns and Davidson iterations unchanged, `lean.fused_borders` 9 per warm solve, `fused_wasted_columns` 0. si64 (load 5 to 21): cold 0.921 (5.0 %) CANDIDATE, warm 1.001 (1.7 %) NEUTRAL. A five-pair si32 rerun at load 7.5 had a 59 % spread (void). These runs fused every append. Final version with the shape rule (binary 190ec3e2, si32, load 11, time void: spread 20 %): 7 of 9 warm appends fuse (reduced-site calls 10 -> 3), cold appends are all wide and stay separate (cold path identical to the default), bit-identical.

### 4.2 Rank 2: narrow NN rotation; `BAANDI_PACK_FIXED=on` (aarch64 pure kernel) and `BAANDI_RITZ_FUSED=loop`
(a) PACK_FIXED. The narrow Ritz rotations (k = 2 to 10 output columns, reduction m = 200 to 460) run at 4.5 to 10 GFLOP/s on the pure kernel: each eight-row strip of V is packed with one variable-length `copy_from_slice` (a `memmove` call) per reduction term, and with few output columns the packing is not amortized. The switch packs full strips and stores full tile columns with constant-length copies, and packs transposed strips (TN: borders, P^T X) in gather order. Pure data movement: bit-identical (same eigenvalue FNV, solver tests). Acts only on the aarch64 pure kernel (`PureGemm::Simd`); x86 pure runs faer's packed GEMM and MKL builds run dgemm, so dev08 is not affected.
Stage 0 (O, si32, bit-identical): first version (contiguous branch only, load 7): cold 0.992 (spread 1.9 %) NEUTRAL, warm 0.976 (1.9 %) CANDIDATE; final version with the transposed gather branch (load 2 to 3): cold 0.963 (2.6 %), warm 0.947 (1.0 %), **CANDIDATE** in both kinds. si64 (first version, load 6 to 21): cold 1.001 (7.0 %), warm 0.999 (2.2 %), NEUTRAL.
(b) RITZ_FUSED=loop. New premise against RITZ-FUSED-LEAN-SCREEN (rejected because the fused return does not form H X, which REUSE-PSI needs): the fused one-pass Ritz and residual kernel only at the loop site, where H X_s is never used again; the return keeps the pair rotation. Microbenchmark (O): fused/pair 0.80 to 1.69, mostly slower on the pure kernel (the fused path's per-panel tile combine costs more than the saved n x k streams). Stage 0 (O, si32): cold 0.999 NEUTRAL, warm 1.019 (spread 1.2 %) **REJECT on the pure build**; bit-identical. The saved streams can only pay where the GEMM itself is fast (MKL): one MKL stage-0 screen in the dev08 request, no further work otherwise.

### 4.3 Rank 3: nonlocal products
No new kernel: on the pure aarch64 kernel the projector products are the same NN (P W) and TN (P^T X) packed products as above, so PACK_FIXED covers them; on MKL builds they are dgemm `mkl_wide` already. The dev08 request asks for their in-process MKL rate at the lean shapes before a nonlocal-specific change (candidates kept in the backlog: fp32 projector storage was rejected by B34; real-space projectors; one GEMV-shaped call for single columns, HAPPLY-NOT-DONE).

### 4.4 Also screened
`BAANDI_NARROW_REDUCED=on` (existing switch, premise changed: lean borders): fires in 1 of 10 warm iterations (only k <= 4), warm 0.980 (spread 4.5 %) NEUTRAL, cold 0.999 NEUTRAL; rounding class (eigenvalues not bit-identical). No further work: the fused border covers the narrow borders without a separate kernel.

## 5. si64 and the combined stack
si64 warm solve shapes (C): (132, 132), (264, 127), (391, 127), (518, 23), then six narrow appends k = 4 to 1 up to m = 556; 11 Davidson iterations, 425 H columns (cold: 3 iterations, 374 H columns). The si64 solve takes 1.4 s cold and 2.5 to 2.9 s warm on the Mac (O).

Combined stack `LEAN_BORDER=fused` + `PACK_FIXED=on` (O, bit-identical in every run): si32 (load 6 to 8) cold 0.989 (3.0 %), warm 0.947 (4.3 %), CANDIDATE; si32 (load 3) cold 1.063 (5.3 %, harness REJECT on the cold kind; pairs 1.085 / 0.979 / 1.063), warm 0.937 (4.7 %); si64 (load 3 to 7) cold 0.987 (1.1 %) CANDIDATE, warm 0.983 (15.9 %, one outlier pair). The Mac was shared with sibling sessions (load 2 to 21). The cold kind of the stack is unresolved: its two appends are wide, and wide fused borders cost 0.93 to 1.17 of the separate products in the kernel microbenchmark (the copies of E and H E are the extra work). These stack runs predate the shape rule; with it the cold solve takes no fused append, which removes the cold question.

## 6. What each fix could give (E)
- Fused border: Mac pure -3.6 % of a warm si32 solve (O); on dev08 MKL the TN borders are a larger share of the algebra and run at a third of the ceiling, so the halved basis stream and the doubled call width are expected to give more there (the narrow TN calls are at 4 to 6 GFLOP/s on MKL for one column, DGEMM-CENSUS); bounded by the border share of the solve (about 20 %).
- PACK_FIXED: Mac pure only, -3.7 % (cold) and -5.3 % (warm) of a si32 solve at low load (O); no dev08 effect.
- RITZ_FUSED=loop: none on the pure build; MKL unknown (one stage-0 screen).
- Not addressed and larger: the algebra work per H column is about twice QE's (basis size), and on the threads/SCF level the extra maps (coarse stage, verification) of profile-2026-11-16-threads.md.

## 7. Reproduction
```
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si32 --emit-input <dir>/si32.toml
target/release/dft-cli run <dir>/si32.toml --output-dir <dir>/si32-cli
BAANDI_DAVIDSON_STAGES=on target/release/examples/davidson_bench --input <dir>/si32.toml --density <dir>/si32-cli/density.bin --repeats 3
python3 -I scripts/fast_screen_dav.py --upf <UPF> --case si32 --density <dir>/si32-cli/density.bin --b BAANDI_LEAN_BORDER=fused --expect lean.fused_borders:up
python3 -I scripts/fast_screen_dav.py --upf <UPF> --case si32 --density <dir>/si32-cli/density.bin --b BAANDI_PACK_FIXED=on
cargo test --release -p dft-numerics --test percore_kernels -- --ignored --nocapture percore_kernel_microbench
cargo test --release -p dft-numerics --test pure_gemm -- --ignored --nocapture pure_gemm_microbench
```
dev08 request: [requests/percore-kernels.md](requests/percore-kernels.md).
