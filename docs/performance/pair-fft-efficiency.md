# Pair FFT efficiency: why the pruned pair path is no faster than the unpruned MKL transform, and what would be

Date: 2026-11-07. Owner: perf-optimizer (measurement and design only; no repository code changed). Labels: **O** = Mac observation (not acceptance), **R** = relayed dev08 number (D08 parts in `tests/reference/perf/dev08_materials_qe_relayed.md`), **E** = estimate from arithmetic stated here, **C** = exact count. Nothing in this document is dev08 acceptance. The scratch crate that produced every O number is `/private/tmp/claude-501/-Users-mino-Develop-rust-dft/6774a210-70b0-490c-9ba0-b3db1cceafa1/scratchpad/fftlab` (not in the repository; session-scoped, so the algorithm is described in section 5 in prose; it links `dft-numerics` by path and uses rustfft 6.4.1). Mac state: Apple M-series, release build, no extra flags, one thread, load average 3.5 to 5.9 during the runs (not idle; all timing claims below use alternating pairs or medians of 40 samples and are comparisons inside one process, but absolute microseconds carry a few percent of load noise).

## 0. Result in ten lines

1. The pair apply does 57 percent of the unpruned flops (C: 4.156 of 7.236 Mflop at 36^3, 8.833 of 15.406 at 72x36x36). The last pass before the potential and the first after it are full-grid and cannot be pruned, and they are 58 percent of the pruned flops; occupation 6.5 percent of the grid buys only a 43 percent flop saving.
2. Where the time goes (O, 40 samples, kernel-skip ablation of a bit-identical clone of `OrderedFusedPlan`): 1D FFT kernel calls 64 to 66 percent, everything else 34 to 36 percent (zero-fill, gathers and scatters of lanes, four scale passes, potential multiply). Pack and unpack are 1.5 percent (2.8 + 2.9 us at 36^3, 5.9 + 6.0 us at 72x36x36).
3. The kernel itself is the floor: rustfft runs the 36-point transform at 15.0 Gflop/s (62 ns per lane) and the 72-point one at 21.6 Gflop/s on the Mac (O). Even with zero data movement the pair apply could not go below 274 us at 36^3 and 510 us at 72x36x36 (O), against measured 373 and 704 for the library.
4. Best restructuring found, "mirror x-slab" (inverse x on sticks, then per x slab y inverse, z inverse, V multiply, z forward, y forward in a 31 KB working set, then forward x on sticks): bit-different but equal to the library within 5e-16 relative; paired, 40 rounds: 0.867 of the library time at 36^3 (323.5 vs 373.1 us) and 0.911 at 72x36x36 (641.8 vs 704.2 us) (O). Across nine shapes (cubic 20 to 48, elongated, two slabs, odd lengths, skewed) the ratio is 0.66 to 0.90 (O, section 4).
5. Everything else tried is neutral or worse (section 4): dropping scale passes alone 0.96, batching rustfft calls 0.99, lane-vectorised prototype kernels 0.92 to 0.96 (worse than the plain mirror on NEON), the existing opt-in `PlaneFusedPlan` 0.98 to 1.02, index-table transposes 0.90 to 0.94 (worse than the plain mirror).
6. On dev08 per band (= per transform equivalent) the mirror would read about 185 us at 36^3 and 444 us at 72x36x36 (E: Mac ratios applied to the relayed 213.6 and 486.9), which beats the unpruned MKL adapter (195.0 and 488.5, R) but is still 1.37 and 1.30 times slower than QE (135 and 341, R).
7. Closing the rest of the QE gap needs a faster kernel, not a better layout: at the measured Mac rate the kernel alone is 137 us per transform equivalent at 36^3, equal to QE's whole call. Only an MKL 1D batched stick pipeline (E: 97 and 235 us, 0.72 and 0.69 of QE) or SIMD lane kernels on dev08 (unmeasured; on the Mac the first prototype only reaches rustfft parity) can do that.
8. New evidence that reopens the MKL route: the "MKL full 195.0 us" figure is the repository adapter, which copies the 373 KB array in and out around an in-place descriptor; QE's raw call at the same grid is 135 us. The two copies plausibly explain a large part of the 60 us difference (E, unverified; cheap dev08 check in section 7).
9. Complex path: the proposed work adds a new plan type and touches no line of `fused.rs`; the complex path stays bit-identical by construction (section 6).
10. Recommended ticket (section 5): `PairFusedPlan` (mirror x-slab, pure Rust), real path only, gate energy 1e-12 Ha, expected 9 to 13 percent of the pair FFT on the Mac, perhaps more on dev08 because the library's two 0.4 to 0.8 MB intermediates spill the 1 MB L2 there and the new working set does not.

## 1. Work per stage (C: counted from the actual sphere)

Fixtures: the `pair_microbenchmark` ellipsoid (q <= 9, scale 12/d): 36^3 has 3071 plane waves (si8 production 2969), 72x36x36 has 6097 (matches si16), grid occupation 6.6 and 6.5 percent. Occupied x lines (one per (y,z) with a plane wave): 253 at both grids. Active z planes 19, active x columns 19 (36^3) and 37 (72x36x36), occupied (x,y) z-sticks 253 and 503. Flop model 5 N log2 N per complex 1D transform.

### 1.1 Current implementation (`OrderedFusedPlan` called through `PairPlan::apply_local_pair`)

Order: pack, inverse x on the 253 compact sticks, inverse y on active z planes (all x), inverse z on all (x,y) lines fused with the potential multiply per y plane, forward x on all (y,z) lines (strided, out of an xz plane buffer), forward y on active x columns, forward z on occupied sticks, gather, unpack.

| stage | 36^3 lanes x length | Mflop | 72x36x36 lanes x length | Mflop |
|---|---|---:|---|---:|
| inverse x (sticks) | 253 x 36 | 0.235 | 253 x 72 | 0.562 |
| inverse y | 684 x 36 | 0.637 | 1368 x 36 | 1.273 |
| inverse z (+ V multiply) | 1296 x 36 | 1.206 | 2592 x 36 | 2.412 |
| forward x (all lines) | 1296 x 36 | 1.206 | 1296 x 72 | 2.879 |
| forward y | 684 x 36 | 0.637 | 1332 x 36 | 1.240 |
| forward z (sticks) | 253 x 36 | 0.235 | 503 x 36 | 0.468 |
| total (lanes) | 4466 | 4.156 | 7344 | 8.833 |
| unpruned pair apply (two 3D transforms) | | 7.236 | | 15.406 |
| ratio | | 0.574 | | 0.573 |

Per lane the code also executes: a gather of the lane into a contiguous buffer (strided except inverse x and forward z), a zero-fill of the missing entries, the transform, a scale by 1/n (inverse) and, for the z lanes, a second scale by the grid size, and a scatter back. Complex element moves outside the kernels (loads plus stores, C): 575,476 at 36^3 (9.2 MB) and 1,146,704 at 72x36x36 (18.4 MB); 53 percent of them are in the inverse-z stage (zero fill, gather of 19 values, two scale passes, multiply, store). Array footprints touched per apply: sticks 142 / 285 KB, `partial` 385 / 770 KB, `forward` 385 / 749 KB, potential 364 / 729 KB (36^3 / 72x36x36): the total 1.3 and 2.5 MB exceeds the 1 MB per-core L2 of the dev08 Xeon Gold 6230 at both grids (and fits the Mac cache, which is one reason the Mac under-reads data-movement cost). Passes over the full-grid-sized arrays: `partial` written once and read once, `forward` written once, read and rewritten by the y pass, rewritten by the z pass and read by the gather (about five touches), the potential read once, sticks four touches; the real grid is never formed.

### 1.2 Ideal stick-pruned algorithm (QE note: z sticks, then y on active x columns, then x on all lines)

Lane count and flops depend on the order of axes: with first axis a on the stick count A_a (projection of the sphere along a), the second pass is over N_a x extent_c lines and the third over N_a x N_b full lines. C table (inverse, one direction; the forward is the mirror image):

| grid | order | lanes (first, second, third) | Mflop per transform |
|---|---|---|---:|
| 36^3 | any | 253, 684, 1296 | 2.078 |
| 72x36x36 | x, y, z (ours, current inverse) | 253, 1368, 2592 | 4.247 |
| 72x36x36 | y, x, z or z, x, y | 503, 684, 2592 | 4.399 |
| 72x36x36 | z, y, x (QE order, FFTW pruned) | 503, 1332, 1296 | 4.586 |

So the QE order is 7 percent more expensive than starting with the 253-stick axis on this elongated cell; the axis with the fewest sticks should come first, and the last inverse pass should be the one that runs along contiguous potential lines (ours: z, because the potential is stored z-fastest). The ideal count is 2 x 2.078 = 4.156 Mflop at 36^3 (equal to the current count: at cubic grids our current order already is the ideal flop count) and 2 x 4.247 = 8.494 at 72x36x36 (3.8 percent below the current 8.833, because the current forward starts with the long x axis on all 1296 lines). In flops the current plan is therefore within 4 percent of the ideal pruned algorithm; the whole difference is data movement and per-lane overhead.

The half sphere does not reduce the sticks: the packed function a + i b has plane waves on both G and -G. The second (-G) stick holds conj(a)+i conj(b) while the first holds a+i b; four real numbers per point are carried by the two sticks, so no stick is redundant, and pairing two real bands in one complex transform is already exactly as efficient per band as a real-to-complex transform of one band would be (half the flops per band). No further halving exists in the pair path (E, argued, not measured).

The unavoidable floor: all pruned algorithms still do the two full-grid passes (58 percent of the flops) because the potential is applied on every real-space point; only an approximation (truncating the potential or coarse grids, QE note 5 section 2) would change that, and it violates the 1e-12 Ha gate.

## 2. Measured stage breakdown (O; kernel-skip ablation of a bit-identical clone of `OrderedFusedPlan`)

The clone (same arithmetic, same operation order, checked bit-identical to the library output at both grids) carries stage timers. A second run replaces every rustfft call by nothing, which leaves exactly the data movement, scaling and multiply. Medians of 40 samples, each the mean of 100 applies; microseconds per pair apply (two bands). The clone is about 11 percent slower than the library in absolute terms (415 vs 373 us at 36^3 and 797 vs 710 at 72x36x36) because of the ablation hooks; read the columns as proportions.

| stage | 36^3 full | movement only | kernel (difference) | 72x36x36 full | movement only | kernel |
|---|---:|---:|---:|---:|---:|---:|
| zero and scatter sticks | 2.3 | 2.2 | 0.1 | 4.9 | 4.5 | 0.4 |
| inverse x on sticks | 17.0 | 1.7 | 15.3 | 28.7 | 3.8 | 24.9 |
| inverse y | 71.7 | 32.6 | 39.1 | 143.5 | 66.6 | 76.9 |
| inverse z, scales, V multiply | 134.6 | 52.8 | 81.8 | 269.7 | 107.6 | 162.1 |
| forward x (gather, scatter active) | 112.3 | 32.1 | 80.2 | 195.9 | 62.4 | 133.5 |
| forward y | 63.1 | 21.3 | 41.8 | 122.0 | 41.5 | 80.5 |
| forward z on sticks | 16.0 | 0.0 | 16.0 | 31.5 | 0.0 | 31.5 |
| final gather | 1.1 | 1.3 | 0.0 | 2.5 | 2.4 | 0.1 |
| sum | 418.2 | 143.9 (34 %) | 274.3 (66 %) | 798.6 | 288.9 (36 %) | 509.7 (64 %) |
| pack plus unpack (outside the plan) | 5.7 | | | 11.8 | | |

1D kernel micro-measurements (O, 256 lanes contiguous, in place, 40 samples of 50 calls):

| length | batched call ns per lane | Gflop/s (5 N log2 N) | per-lane calls | out of place | scale pass alone |
|---:|---:|---:|---:|---:|---:|
| 25 | 48.2 | 12.1 | 50.2 | 50.2 | 6.3 |
| 36 | 62.0 | 15.0 | 65.3 | 64.5 | 8.2 |
| 48 | 83.6 | 16.0 | 88.1 | 85.3 | 10.5 |
| 72 | 102.7 | 21.6 | 109.1 | 106.9 | 17.1 |

Cross-check: the lane counts of section 1 times these per-lane times give 277 us (36^3) and 519 us (72x36x36), equal to the measured kernel columns (274, 510); the kernel stage cost is entirely the 1D transforms, not call overhead (batched versus per-lane calls differ by 5 percent per lane). Per stage the movement is largest where the lane is gathered from a strided source and scattered back: inverse z 39 percent of its stage, inverse y 45 percent, forward y 34 percent (the zero fills and the two scale passes add to the z stage). So: kernel 64 to 66 percent; gather, scatter, zero-fill and scale 31 to 33 percent; potential multiply plus final gather about 3 percent; pack plus unpack 1.5 percent of the apply.

## 3. Why pruned equals unpruned MKL (R numbers, per transform equivalent = half a pair apply)

| grid | BAANDI pruned pair (R) | MKL full adapter (R) | QE per call (R) | flops of the BAANDI pruned transform | rate | flops of the unpruned transform | MKL rate | QE rate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 36^3 | 213.6 us | 195.0 | 135 | 2.078 Mflop | 9.7 Gflop/s | 3.618 | 18.6 | 26.8 |
| 72x36x36 | 486.9 | 488.5 | 341 | 4.417 | 9.1 | 7.703 | 15.8 | 22.6 |

The pruned path does 57 percent of the flops at about half the flop rate of the vendor transform, so the two cancel. Half of that rate loss is the data movement of section 2 (it would take the rate from 9.7 to about 14.6 if the movement were free); the other half is that rustfft's 36-point kernel itself (15 Gflop/s on the Mac) is slower than the vendor's whole-3D rate (26.8 at QE). The MKL adapter's own 195 us contains two full-array copies around the in-place descriptor (`dft_mkl::Plan3d::transform`), which QE's raw call does not pay (E: two 373 KB copies are of the order of 20 to 40 us at cache bandwidth); this accounts for part of the 195 versus 135 difference and is cheap to check (section 7, item 3). The MKL unpruned baseline is therefore itself beatable by the adapter's overhead alone.

## 4. Candidates (O; paired alternating rounds, 40 rounds of 50 applies for the table, medians; microseconds per pair apply; all outputs checked against the library, max relative difference 3e-16 to 8e-16)

| ID | candidate | 36^3 | 72x36x36 | ratio to library | verdict |
|---|---|---:|---:|---|---|
| L | library `OrderedFusedPlan` through `apply_local` (reference) | 373.1 | 704.2 | 1.000 | |
| M | **mirror x-slab**, batched rustfft calls | 323.5 | 641.8 | **0.867 / 0.911** | best, adopt |
| M-pl | same with per-lane rustfft calls instead of batched | 327.7 | 653.9 | 0.876 / 0.921 (separate run) | batching worth 1 percent |
| M-ix | M with an index-table transpose instead of zero-fill plus scatter | 334.4 | 664.0 | 0.896 / 0.943 | worse than M |
| ML | M with the z and y slab transforms replaced by a lane-vectorised SoA prototype kernel (four adjacent x slabs per group) | 341.9 | 676.5 | 0.916 / 0.961 | worse on NEON; dev08 AVX unmeasured |
| PF | existing `PlaneFusedPlan` (z first, plane-fused, opt-in A27) | 381.9 | 692.0 | 1.024 / 0.983 | neutral; confirms A27 |
| NS | clone with all per-lane scale passes dropped (same order; rounding-only change) | 400.3 vs clone 415.3 | 761.1 vs clone 797.0 | 0.964 / 0.955 of the clone | about 4 percent of the old structure, not alternated |

Stage view of M (O, timers on, 40 samples of 100 applies): at 72x36x36 kernel calls 540 us (y 169, z 319, x 52), slab gathers, zero fills and transposes 87, potential multiply 18, stick scatter and final gather 9: sum 653 us, kernel 83 percent, everything else 17 percent (36^3: kernel 274 of 330 us, 83 percent). The movement share fell from 34 to 36 percent to 17 percent, the kernel time is unchanged (lane counts equal at 36^3; 8.494 vs 8.833 Mflop at 72x36x36).

Mechanism of M (why it wins): the middle of the apply never forms a full real-space grid or a transposed full-grid intermediate; after the inverse x pass every x slab (y,z plane) is independent, so y inverse, z inverse, potential multiply, z forward and y forward run on one 31 KB slab that stays in L1; the forward z pass runs on the same lanes right after the inverse z pass, so the real-space potential multiply reads the potential z-contiguously with no transposition; the forward x pass becomes a contiguous batched call over only the 253 occupied sticks (the library's forward x over all 1296 lines of length 36, strided, is gone: 112 us of 418 at 36^3). The stick array (292 KB at 72x36x36) is updated in place by the slabs, so the whole apply keeps one 292 KB array plus the slab and the streamed potential (5.3 and 10.7 MB of element moves against 9.2 and 18.4 MB, C).

Generality (O; paired 30 rounds of 30 applies, library vs M; npw = sphere size, occupation = npw / grid):

| shape | npw | occupation | library us | M us | ratio |
|---|---:|---:|---:|---:|---:|
| cubic 20x20x20 | 515 | 0.064 | 58.1 | 47.2 | 0.812 |
| cubic 25^3 | 1021 | 0.065 | 137.1 | 117.7 | 0.859 |
| cubic 36^3 | 3071 | 0.066 | 371.9 | 323.4 | 0.869 |
| cubic 48^3 | 7153 | 0.065 | 874.7 | 773.2 | 0.884 |
| elongated 72x36x36 | 6097 | 0.065 | 700.1 | 631.4 | 0.902 |
| slab 36x36x18 | 1565 | 0.067 | 194.0 | 166.9 | 0.860 |
| slab with vacuum 36x36x72 (sphere short in z) | 1985 | 0.021 | 595.8 | 395.4 | 0.664 |
| odd lengths 15x20x24 | 465 | 0.065 | 35.8 | 27.6 | 0.769 |
| skewed 36^3 (shear 0.35 in the metric) | 3053 | 0.065 | 367.0 | 310.7 | 0.847 |

M is never slower; the gain is largest where the sphere occupies few planes of a long axis (the vacuum slab, which is exactly the case where the current inverse-x-first ordering leaves a full-length strided forward x pass). Maximum relative difference to the library 3e-16 to 6e-16 in every row.

Things that did not help and the reason (each a candidate from the task): (a) a layout where every transform is contiguous by full out-of-place transposes: M is better, because the slab working set makes the y transposition L1-resident and a full transposition adds a full-array pass; (b) fusing the potential multiply into the last pass: already fused in the library (the multiply is inside the z-lane loop), in M it is a separate 31 KB pass costing 2.7 percent, so there is nothing to win; (c) skipping empty columns in the y and x passes: the library already skips them (684 of 1296, 1368 of 2592 y lanes), M keeps it; (d) half sphere to halve sticks: no redundancy exists (section 1.2); (e) batched rustfft calls: 1 percent, kept because it is free; (f) plan-per-axis reuse: both implementations already plan each axis once; (g) larger groups of bands or planes: A28 stays rejected, M works per pair; (h) A29 plane pruning rule: M prunes without a threshold (empty planes and columns are removed in the first pass, and no strided 2D alternative exists in pure Rust, so the 60 percent rule does not arise); (i) A11 four-lane staging: ML is the lane-vectorised version of that idea and loses on the Mac.

Lane-vectorised kernel prototype (O, E3): a plain Stockham mixed-radix 4x3x3 kernel written with fixed-width arrays that the compiler vectorises across lanes reaches 58.3 ns per lane (16.0 Gflop/s) at length 36 with four lanes, against rustfft 62.0; at length 72 it is slower (147 against 103 ns). It is untuned (no special first-twiddle case, generic index arithmetic) and the Mac has 128-bit vectors; on dev08's 256 and 512-bit vectors a lane-vectorised kernel has more to gain than rustfft's already wide-vector kernels, but that is unmeasured and must come from a dev08 run (section 7).

## 5. Which ordering reaches or beats unpruned MKL, and what the faithful implementation reaches

Per transform equivalent (half a pair apply), dev08 reference numbers (R) and projections (E = relayed number times the Mac ratio; the Mac ratio understates the benefit on a CPU with a 1 MB L2, E):

| grid | current (R) | M projected | unpruned MKL adapter (R) | QE (R) | pure floor: kernel only at Mac rate scaled by the dev08/Mac ratio | MKL 1D stick pipeline at QE's flop rate (+25 percent) |
|---|---:|---:|---:|---:|---:|---:|
| 36^3 | 213.6 | 185 | 195.0 | 135 | 137 x 1.14 = 156 | 2.078 Mflop / 26.8 = 78 -> 97 |
| 72x36x36 | 486.9 | 444 | 488.5 | 341 | 270 x 1.38 = 373 | 4.247 / 22.6 = 188 -> 235 |

- M reaches and beats the unpruned MKL adapter at both grids (185 < 195.0 and 444 < 488.5, E), by 5 and 9 percent. It does not reach QE (1.37 and 1.30 times slower, E).
- Even with perfect data movement the pure rustfft pair path is at best at parity with QE: the kernel-only floor is 156 and 373 us (E), 1.16 and 1.09 times QE. Exceeding QE in pure Rust therefore needs SIMD kernels faster than rustfft's on the target CPU (an unmeasured chance, dev08 microbenchmark first) and cannot come from layout.
- A faithful stick pipeline on MKL 1D batched descriptors (z or x on sticks with distance and stride, y on active columns, x/z fused with the potential; descriptors owned by the plan, sequential layer, as the fastqe and QE notes describe) is expected at 78 to 97 us (36^3) and 188 to 235 us (72x36x36) per transform equivalent if the strided batched 1D rate reaches the vendor 3D rate that QE measured, i.e. 0.55 to 0.72 times QE's per-call time and 2.1 to 2.7 times faster than the current path. The assumption (strided 1D batches at 3D efficiency, per-descriptor commit cost already paid at plan time) is untested here; the fastqe report measured 5x against its unpruned baseline on other hardware. This stays a dev08-only measurement (Mac has no MKL). Because the current MKL build routes the real path through the copying adapter and full-grid zero-fills (`FusedPlan` alias under the `mkl` feature), the first MKL step is just a raw in-place pair path (no copies), which should already reach QE's 135 us per transform at 36^3 (E).

### Ticket-ready design (T1, pure Rust, real path only)

Name: `PairFusedPlan` in a new file `crates/dft-numerics/src/fused_pair_plan.rs` (new file; `fused.rs` and `fused_pair.rs` untouched except one `pub mod` and one re-export line in `lib.rs`).

Interface: `PairFusedPlan::new(shape, ordered_indices) -> Result<Self, FftError>`; implements `PairDensityOperator` and, because the trait requires it, `FusedOperator` (`apply_local` as the pair kernel; the single-band `accumulate_density` computed from the same inverse half, unused by the real path). Plugged in by the closure at `crates/dft-physics/src/real_gamma.rs` line 194 (`PairLocalWorkspace::new(map, capacity, pool, || FusedPlan::new(dims, g_to_fft))`), under `#[cfg(not(feature = "mkl"))]`; the density workspace (`PairDensityWorkspace`) takes the same plan type because it is generic over `PairDensityOperator`.

Algorithm (all arithmetic in `f64`, one plan per fixed pair piece, no allocation after construction, no threads inside):
1. Plan construction (C, once): choose the axis order automatically: the first axis is the one with the fewest occupied lines (ties: lowest index); the last axis is the one along which the potential is contiguous (z); the middle axis is the third. All counts below use these axis roles (first = x here); for the production grids with z last this is x, y, z. Build: stick slots numbered plane-major (plane by active last-axis index, then y within the plane); `scatter[i] = slot * n_first + coordinate`; per plane the list of (y, slot); active-plane list for the last axis; rustfft plans (inverse and forward, per axis length) and one scratch sized to the maximum in-place requirement; buffers: sticks (slots x n_first), a y-batch buffer (active planes x n_mid), a z slab (n_mid x n_last).
2. Apply: (a) zero the stick array and scatter the packed coefficients into it (rustfft out-of-place calls destroy their input, so there is no persistent-zero buffer; the 292 KB zero fill costs 1 to 2 percent); (b) one batched in-place inverse transform over all sticks; (c) for each slab of the first axis: zero the y-batch buffer and gather the occupied entries of the slab from the sticks (stride n_first), one batched inverse transform over the active planes, transpose into the z slab (zero fill, then scatter of the active last-axis indices), one batched inverse transform over the n_mid lines, multiply by the potential slab (contiguous because the potential is stored with the last axis fastest and the first axis slowest), one batched forward transform, extract the active indices into the y-batch buffer, one batched forward transform, scatter only the occupied entries back into the sticks in place; (d) one batched forward transform over the sticks; (e) gather with the final factor 1/(n_x n_y n_z) (all unnormalised transforms: the per-lane scale passes of `OrderedFusedPlan` are dropped; the product of its four scale factors is one, so the result is the same operator).
3. Density: the same inverse half (a to c up to and including the z inverse) followed by accumulating the weighted squared real and imaginary parts per grid point instead of the multiply and the forward half, so the pair density path no longer needs the library inverse (this is the second user of the same code; gate it by the existing pair density tests).
4. Odd tail band: `PairMap::pack` with the zero partner is unchanged.
5. Determinism: one slab order, one stick order, no reductions across pieces; results are bit-identical across thread counts and across repeated calls, and independent of the persistent worker pool.

Validation gates (all must pass; none relaxed): `scripts/check.sh` exit 0; a new test comparing `PairFusedPlan::apply_local` against `OrderedFusedPlan::apply_local` on at least nine fixtures (the shapes of the table in section 4, including a skewed metric and an odd length such as 15x20x24, and a sphere with no plane waves at some planes) with max relative difference below 1e-13; the same for pair density against the existing `OrderedFusedPlan` pair density; stable buffer pointers across calls; bit-identical outputs across 1, 2, 4, 8, 12 threads through `PairLocalWorkspace`; `scripts/bench_scf.py --cases gamma,k222s,k222,k444s,si8 --require-identical` for the COMPLEX path: energy and density bit-identical to the parent commit (the complex path does not use the new type); for the real path: energy within 1e-12 Ha of the stored baseline and the same SCF iteration and certificate behaviour on gamma and si8 and on the R5 synthetic cases; generality set in `scripts/bench_suite.py` (cubic, skewed, elongated, slab, 4/8/20 atoms, synthetic projectors); microbenchmark extension of `pair_microbenchmark` with the new type next to the old one on 36^3 and 72x36x36 (30 or more samples, alternating).

Expected gain: Mac 9 to 13 percent of the pair FFT (13 percent at 36^3, 9 percent at 72x36x36), 13 to 34 percent on other shapes; FFT is 35 to 39 percent of the real+standard working time (R, D0), so the SCF gain is about 3 to 5 percent on the Mac scale; on dev08 plausibly more (cache) but that is a hypothesis. Honest reading: this closes about a quarter of the pair-FFT gap to QE, not the gap.

Ticket T2 (needs a dev08 microbenchmark first; MKL builds): raw in-place MKL pair path (no adapter copies, no full-grid zero fill beyond what the transform needs), then the 1D batched stick pipeline of section 5. Ticket T3 (needs a dev08 microbenchmark first): lane-vectorised kernels (SIMD across slabs, four to eight lanes, radix 2/3/4/5 butterflies) behind the same `PairFusedPlan` interface; expected only if the dev08 prototype beats rustfft at lengths 25 to 72.

## 6. Bit identity and what the complex path shares

- `OrderedFusedPlan` (fused.rs lines 91 to 295) is used by the complex path for both `apply_local` and `accumulate_density`, and also by the real pair path today: `PairPlan::apply_local_pair` calls the plan through `FusedOperator`, and `PairDensityOperator for OrderedFusedPlan` (fused.rs lines 582 to 607) reuses its `inverse_xy` and `inverse_z`. The `Kernels` struct (per-lane transform plus 1/n scale) is private to it and to `PlaneFusedPlan`. The complex band density (`dft-engine/src/band_density.rs`) uses `FusedPlan::new`, the same alias.
- A new type in a new file shares none of these private functions; the complex path therefore stays bit-identical without testing change, as long as `fused.rs` and the `FusedPlan` alias are not edited. The ticket must not alias `PairFusedPlan` into `FusedPlan`.
- The real path has no bit-identity contract with the complex path (the real-path energies already differ in the last digits by design); the change of rounding (sum order of the transforms, removed scale passes) is a stated rounding change under the P14B rule: energy within 1e-12 Ha, iteration counts and certificate unchanged on realistic fixtures, reported if not. Determinism across thread counts within a backend is kept (no thread-dependent reduction).
- Under the `mkl` feature `FusedPlan` is `GridFusedPlan<FftPlan3d>`; T1 does not apply there (cfg gate), T2 addresses it.

## 7. dev08 request (ask the orchestrator to relay; acceptance only from dev08)

1. Microbenchmark on dev08 native build (`scripts/build_native.sh` flags, node idle, load below 0.30): rustfft 1D kernel rates at lengths 25, 36, 48, 72 (the E1 table of section 2), the library pair apply, M and ML at 36^3 and 72x36x36 and 25^3 (the scratch crate runs all of it: `fftlab e1`, `fftlab paired`, `fftlab e3`, `fftlab gen`; it needs its manifest path to `dft-numerics` edited to the dev08 checkout). Report medians of at least 30 samples and the cache-spill hypothesis (library versus M ratio at 72x36x36 on the Mac 0.91; if dev08 shows 0.8 or lower the L2 effect is confirmed).
2. In the same session the lane-vectorised prototype (`fftlab e3`) with AVX2 or AVX-512 enabled by `target-cpu=native`: if it exceeds rustfft by 1.3 times at lengths 36 and 72, T3 is justified.
3. MKL checks: raw in-place DFTI 3D at 36^3 and 72x36x36 without the adapter's two copies, with the same descriptor settings as `Plan3d` (expected near QE's 135 and 341 us if the copies explain the 195 and 488 us), and one 1D batched strided inverse of length 36 over 253 contiguous sticks plus 684 strided lines, to measure the rate T2 depends on.
4. After T1 lands: the standard three-point alternation (parent, main, branch), si8 and si16 real strict and standard T=1 and T=8, bit-identical iteration counts expected within 1 iteration, `perf` profile with at least 30 samples on `apply_local`.

## 8. Ledger check (premises)

- A29 (plane pruning only below 60 percent active columns): not revived; M has no strided 2D alternative to choose from, and the 60 percent rule belongs to vendor 2D transforms. A27/PF measured again: neutral (1.02 and 0.98), consistent with its row.
- A28 (grouping bands or planes), A34 (packed-lane FFT threads), A35 (scoped threads), A11 (four-lane staging): untouched. M processes one pair per plan, uses no threads, and its lane group is a slab, not an adjacent-lane staging of strided gathers; ML is the A11-like variant and lost on the Mac, which agrees with A11. New evidence only in two places: (1) the A11 premise "strided gather is not the bottleneck on the Mac" is correct for full 3D transforms of 24^3 to 36^3 but here the movement is 34 to 36 percent of the pruned pair apply, because pruned lanes are fewer and cheaper per lane while their gather and scatter cost is unchanged (the kernel is 62 ns per lane, the lane movement plus scales about 33 ns); (2) the dev08 L2 (1 MB) is smaller than the library's working set at both grids (1.3 and 2.5 MB), which the Mac L2 hides.
- A26 (P10C fused ordered arithmetic): its bit-identity design is the reason the library has four scale passes and an asymmetric order; M deliberately gives that up for the real path only.
- A24/A25/A41 (MKL backend): the copy cost inside `Plan3d::transform` is new evidence for B08/B28; it does not reopen the MKL default (A57: native pure beats MKL on si8).
- A52 and the QE note's S-A to S-E: T1 is the pure-Rust counterpart of S-A/S-B restricted to the pair path; S-C (density) is included as item 3 of the design.

## 9. Limits

Mac observations under a non-idle machine (load 3.5 to 5.9; ratios from alternating rounds are more reliable than absolute times); one machine type; no dev08 data for M, ML, or the 1D kernel rates; the stage breakdown was taken on the clone (11 percent slower than the library in absolute terms, proportions used); the kernel-minus-movement split is by ablation (movement-only run), so cache-state differences between the two runs are included in "kernel"; the lane-vectorised prototype is untuned; the MKL stick-pipeline numbers are arithmetic, not measurements; the half-sphere argument is analytic. No QE source was needed beyond the existing notes for this report; the stick structure (z sticks, y on active x columns, x on all lines, FFTW-pruned, MKL-unpruned) is as in `docs/research/qe-study/qe-hpsi-and-fft.md` section 1.
