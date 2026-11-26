# Performance attempts ledger

Compiled 2026-11-05 from the repository documents named in each row. Read this before proposing or implementing any performance change. Companion files: [AGENT-BRIEF.md](AGENT-BRIEF.md) (process), [backlog.md](backlog.md) (ideas not yet tried), [optimization-log.md](optimization-log.md) (raw evidence).

## A. How to read this, and the rule

- One row per attempt: adopted, rejected, superseded, or opt-in. "Mac" = shared Apple Silicon Mac (12 logical CPUs, indicative only). "dev08" = 2 x Xeon Gold 6230, Rocky 9.7, pinned, quiet; dev08 numbers are RELAYED by the remote session (not reproduced in this repository) unless a row says otherwise.
- Status words: **adopted** (in main as default), **opt-in** (in main, off by default), **rejected** (not merged or removed), **superseded** (kept at the time, replaced by a later row), **gated** (not rejected, but its precondition measurement failed or was not run).
- "not recorded" means the number was searched for in the sources and not found. Numbers marked "derived" are arithmetic on recorded numbers (shown in the row).
- Source column: `OL` = optimization-log.md, `R` = README.md in this folder, `N1..N4` = docs/research/qe-study notes (N1 davidson-and-scf-thresholds, N2 fft-call-structure, N3 davidson-work-accounting, N4 davidson-gemm-accounting), `D08` = tests/reference/perf/dev08_final_relayed.md (parts 1 to 5), `D08B` = dev08_baseline.json, `D08T` = dev08_threads.json.

**Rule: never repeat a rejected, superseded or gated idea without writing down the new evidence that changes its premise** (a different system class, a different backend, a different measured bottleneck, or a variant that differs in the mechanism that failed). "It was tried on the Mac" or "it was only tried on Si" is a premise-change only if the row says so. Rejection reasons below are specific: read the mechanism before arguing the idea is new. When you finish an investigation, add rows here and update [backlog.md](backlog.md).

## B. Attempts

Notation: `->` before to after. Wall seconds unless stated. SCF = self-consistent iteration. 2mnk = complex GEMM work counted as 2*m*n*k.

### B1. Early FFT, Davidson and SCF work (Nov 3, 2026)

| ID | Ticket / idea | Status | Measured effect (where) | Mechanism: why it worked or failed | Lesson | Source |
|---|---|---|---|---|---|---|
| A01 | P7P1 (O1): plan FFT once, persistent scratch, `process_with_scratch`, no per-call allocation | adopted | Mac: Gamma 1.304 -> 0.960 (1.36x), k222s 6.549 -> 4.623 (1.42x), energy and density bit-identical. Profile before (6445 samples): FFT 48.0%, allocation about 30% (malloc/free/memset 15.3% plus mach_absolute_time 15.0% inside macOS free); after (5020 samples): FFT 75.4%, allocation 2.5%. dev08: measured together with A02 | ndrustfft allocated lane and process scratch on every call; macOS free calls the timer | Allocation in hot loops is the first thing to remove; bit-identical | OL "O1 / P7P1-FFT-ALLOC" |
| A02 | P7P2 (O3): pruned FFT, skip lanes and planes that are exactly zero for the wavefunction sphere | adopted | Mac: Gamma 0.907 -> 0.643 (1.41x), k222s 4.499 -> 3.007 (1.50x), k222 8.871 -> 5.799 (1.53x); 1001 of 1728 lanes executed on 24^3 radius 6 (-42.07%). dev08 (commit f6e3679, A01+A02 together, relayed): Gamma 2.278 -> 1.219 (1.87x), k222s 10.972 -> 5.602 (1.96x), k222 21.236 -> 10.735 (1.98x), k444s 98.975 -> 46.975 (2.11x), identical to baseline | Zero lanes do no work; dependencies cannot cross absent coordinates | Skipping exact-zero work is bit-identical and safest | OL "O3 / P7P2"; D08B bench_steps |
| A03 | P7P5 (O2): hypot-free norm (4 accumulators, sqrt only when safe), 4-accumulator projection dots | adopted | Mac, shared host: SCF Gamma 0.157 -> 0.162, k222s 0.572 -> 0.511, k222 1.066 -> 0.955 (norm only: 0.180 / 0.531 / 1.736); ignored kernel test (2048 entries, 20000 calls) norm 0.689 -> 0.0147 s (overlapped Clippy, indicative). Profile before: hypot 12.2%, ortho 8.3%, GEMM 8.7% (2330 samples). dev08: not recorded | Complex `hypot` per element was 12% of the profile; rounding changes at 1e-16 so not bit-identical (energy deltas <= 4.4e-15 Ha) | A 47x kernel gain did not show as a measurable SCF gain on a noisy host; no SCF claim was made. Rejected inside: block GEMM for ortho (copies through adapter), skipping 2nd MGS pass (no error bound) | OL "P7P5 O2" |
| A04 | P7P3 stage 1 (O4a): Davidson keeps H*v products of unchanged basis columns instead of re-applying H to the whole basis every iteration | adopted | Mac: vector applications 3281 -> 1497 (2.19x; k222s 2.32x, k222 2.34x); wall Gamma 0.938 -> 0.490, k222s 4.828 -> 2.224, k222 10.036 -> 4.273; bit-identical to baseline. dev08 (3b94d53, with A05): Gamma 1.219 -> 0.710, k222s 5.602 -> 2.808, k222 10.735 -> 5.341, k444s 46.975 -> 22.409 | Old code applied H about 40 times per iteration where 16 suffice (basis 48, block 8) | Count operator applications first; it was a 2x redundancy | OL "P7P3 / O4a and O4b"; D08B |
| A05 | P7P3 stage 2 (O4b): residual from stored products, `(H V) C`, one true block application only to certify | adopted (default) | Mac: applications 1497 -> 833 (Gamma), wall 0.490 -> 0.338, k222s 2.224 -> 1.427, k222 4.273 -> 2.733; energy within 8.0e-15 Ha of stage 1; density hashes differ (rounding) | Saves one H application per residual; failed certification falls back to true residual | Rounding change breaks historical density hashes: re-baseline via `baseline-after-o4.json` | OL "P7P3" |
| A06 | P7P4 (O4c) a: SCF eigen-tolerance schedule `min(0.1 * L2 density residual, 1e-2 Ha)`, floor `0.1 * eigen_tolerance` | adopted | Counts alone (before -> schedule): Gamma apps 833 -> 592, k222s 3863 -> 2604, k222 7437 -> 4924, k444s 31963 -> 21078; SCF steps 11 -> 10 on k222s and k222 | Early SCF steps do not need tight eigenvectors | Schedule constant c=0.1 not swept; it is empirical | OL "P7P4 / O4c" |
| A07 | P7P4 b: cold start from the `bands` lowest-kinetic-energy plane waves instead of random vectors | adopted | Alone: Gamma 833 -> 690, k222 7437 -> 6145. Combined A06+A07 (apps before -> after): Gamma 833 -> 422, k222s 3863 -> 1883, k222 7437 -> 3661, k444s 31963 -> 15015 (about 2.0x). Mac wall Gamma 0.266 -> 0.156, k222s 1.104 -> 0.568, k222 2.090 -> 1.062. dev08 (ec77e7e, relayed): Gamma 0.710 -> 0.427, k222s 2.808 -> 1.457, k222 5.341 -> 2.666, k444s 22.409 -> 10.606 | First-step Davidson iterations per k point 29 -> 5 (Gamma) | Initial subspace quality matters most in the first SCF step only | OL "P7P4"; D08B |
| A08 | P7P4: subspace capacity 32 / 64, block size 4 | rejected | cap 32: apps Gamma 462 vs 422, k222 4086 vs 3661; cap 64: identical to 48; block 4: iterations Gamma 98 vs 50, wall not consistently better | Smaller capacity forces restarts; block 4 needs more iterations | Keep capacity 48, block = bands **Premise changed (V 6.1, 2026-11-05):** capacity is `max(48, 6*bands)`, so 8-band cases have 48 but si8 (20 bands) has 120; these sweeps used 8 bands only. | OL "P7P4" |
| A09 | P7P6: one scoped thread group per balanced block of k points (outer parallelism), ordered reductions | adopted | Mac: k444s T=1/2/4 5.425 / 2.717 / 1.882 (2.0x, 2.9x); k222s 1.997x / 2.687x at T=2/4. dev08 (1decf9c, relayed): k444s T=1 9.346, T=4 4.270, T=8 3.450, T=16 3.057 (3.06x); Gamma flat 0.403 -> 0.397; bit-identical densities across T | Independent k solves; reductions kept sequential in k order | Single-k (Gamma) is not helped by k parallelism. Rejected: parallel density/energy reductions (order), nested FFT or faer threads, default all cores, global mutex | OL "P7P6"; D08T |
| A10 | P7PA: parallel per-k setup, density partials (folded in ascending k), energy audit, certificate | adopted | Mac T=8, time outside diagonalization: k222s 0.129 -> 0.067, k222 0.219 -> 0.080, k444s 0.783 -> 0.209 (-73%). dev08 (20cb416, relayed): k444s T=16 engine 3.043 -> 1.211 (diag 0.839), wall 1.226 (7.56x vs T=1); 4 and 8 k points best at T=8, slightly worse at T=16 (not investigated) | At T=16 the non-diagonalization part was 2.19 s of 3.043 s engine (derived from D08T phases), all sequential | After parallelising the main loop, Amdahl moves to setup/density/audit | OL "P7PA"; D08T |
| A11 | P7P9 (H8): four adjacent strided FFT lanes gathered into one staging block | rejected | Transform alone: 1.01-1.04x on 24^3 to 36^3, 0.85-0.96x tiny grids; end to end A/B five runs min 1.023 vs 1.019, median 1.069 vs 1.051 (1.00-1.02x); W=8, 16 not better. dev08: not run | Strided gather is not the bottleneck on the Mac | Reconsider only if FFT is dominant again or on a CPU where strided access is much costlier **Premise changed (H 4b):** about 3800 gather/scatter lane operations per band at 36^3; dev08 was never run; reopens only as backlog B18. | OL "P7P9"; N1 section 6 |

### B2. Davidson hypotheses from the first QE study note (P7P7, P7PC, P14A)

| ID | Ticket / idea | Status | Measured effect (where) | Mechanism | Lesson | Source |
|---|---|---|---|---|---|---|
| A12 | P7P7 S1 (H3a): restart by contracting stored basis and H products instead of re-orthonormalising and re-applying H | adopted | Default cases: no change (they never restart at capacity 48). Test problems: applications 165 -> 114 and 92 -> 64. Later relevant: si8 has 14 restarts (P14A) | Saves one H application per retained Ritz column per restart | Matters for larger systems that do restart | OL "P7P7"; N1 section 6 |
| A13 | P7P7 S2 (H3b): subspace capacity 2x, 3x, 4x band count (instead of 6x) | rejected | Wall medians G/k222s/k222 vs 0.1404/0.4997/0.9269: 2x 0.2430/0.5944/1.1100; 3x 0.1540/0.4926/0.9399; 4x 0.1400/0.4710/1.1011; applications 2x 452/2018/4076 vs 373/1730/3398 | Our basis is orthonormalised with corrections appended; a larger basis converges in fewer applications | QE parameter choices do not transfer where the algorithm differs **Premise changed (V 7, 6.1):** run with 8 bands and with all roots always rotated and orthonormalised; at 20 bands capacity is 120 (mean basis about 70 vs QE about 30) and the sweep was never repeated; to be retested only after locking or a non-orthonormal basis (backlog B17). | OL "P7P7"; N1 |
| A14 | P7P7 S3 (H5): fuller diagonal preconditioner (kinetic + mean local potential + nonlocal diagonal, smooth positive denominator) | rejected | Gamma needed 11 SCF steps instead of 10; applications 427/1799/3461 vs 373/1730/3398 | No iteration reduction on these cases | Do not retry without a case where the kinetic diagonal is known to fail | OL "P7P7"; N1 |
| A15 | P7P7 S4 (H1): Ritz-value-change convergence test inside the loop, with strict residual refinement before return | rejected | Iterations 45/221/441 vs 43/212/422; k222 saved 38 applications but added 19 iterations | The strict refinement pass re-did the work the cheaper test saved | A different variant (no refinement pass, certificate only at the end) was not tried: see backlog B04 **Premise changed (V 7):** the cost-oriented use (small problem decides which roots get work, certificate kept) is untested; QE's test also tolerates up to 5 unconverged roots (S 5.4). **Update 2026-11-06:** the Ritz-change rule inside an uncertified profile is P17A-3; the P15A rows below show it explains most of the step gap. | OL "P7P7"; N1 |
| A16 | P7P7 S5 (H4): generalized non-orthonormal reduced problem | gated | Not attempted: prerequisite stages were negative, highest numerical risk | Would remove orthonormalisation and hypot cost (about 20% of the post-O4c profile) | Not rejected **Premise changed (V 7, backlog B13):** now a proposed route with a mandatory rank/conditioning diagnostic; A52's failure mechanism (rank decision vs rounding) is undetermined. | OL "P7P7"; N1 |
| A17 | P7PC: buffer bands. Strict tolerance for occupied bands plus the first empty band; buffer bands floor-only `max(tol, 1e-5)`, capacity at least 6 columns per band (48), recurrence-only acceptance of intermediate tolerances | adopted | si8 before: fails (exit 4, 400 iterations, residual 4.24e-6 Ha); after: 15 SCF steps, 283 iterations, 2060 main applications, E = -31.321803601031377 Ha. Gamma/k222s/k222 main applications 373 -> 281, 1730 -> 1266, 3398 -> 2543 (certificate 49 -> 38, 153 -> 113, 263 -> 208); wall 0.140 -> 0.123, 0.500 -> 0.410, 0.935 -> 0.788 (Mac) | Band 20 cut a six-member near-degenerate cluster (gaps about 1.5e-6 Ha), so the top band never reached the strict tolerance | A hard failure on the 8-atom cell was found only after widening beyond 2-atom cases: test shapes broadly | OL "P7PC" |
| A18 | P7PC: looser buffer factor 100 x tolerance, capacity 60 or 120 | rejected | Converges only with capacity 60/120 but more SCF steps (Gamma 10 -> 11, si8 15 -> 17: 701 / 311 iterations); floor-only at capacity 48 with strict selection fails (max 1.55e-4 Ha) | Loosening only buffers cannot certify the cluster | Keep floor-only target | OL "P7PC" |
| A19 | P7PC / H7 / P14A: fewer or more extra bands as default (1 or 2 instead of 4) | rejected | P7PC: si8 2198 / 2567 apps vs 2060, wall 5.597 / 6.461 vs 4.630. P14A: 1 extra: k444s H vectors 11181 -> 8310 (10 SCF), si8 2416 -> 2444, skewed 20-atom fails to converge in 100 SCF; 2 extra: k444s 9068, si8 17 SCF and 3249, skewed converges to an energy +29.3839772182878 Ha away | Saves on many-k, loses on 8 atoms; skewed cell changes fixed point | Default 4 stays; QE comparison inputs use nbnd=8 explicitly | OL "P7PC", "P14A"; N3 |
| A20 | H2 (looser threshold for empty bands) | adopted (via A17) | Buffer-origin share of H vectors (P14A): si8 24.92%, k444s 34.39%, skewed 10.98%, elongated 28.65% | Realised as buffer-band policy; lowest empty band kept strict for the gap | Closed | N3 |
| A21 | P14A (H6 premise): SCF eigen-tolerance schedule multiplied by factor 3/10/30/100 on the density-residual term (cap 1e-2 and floor 1e-9 unchanged) | rejected (this family only) | si8 SCF count 15 -> 16/16/21/57 (H vectors 2416 -> 2072/1794/1722/2280); k444s 10 -> 12/15/35/58 (11181 -> 10582/10896/16781/22046); 20-atom skewed diverges at 3 and 10; elongated 9 -> 17/19/27/32 | Looser eigenvectors raise the number of SCF iterations; fewer Davidson steps are paid back | H6 itself (energy-weighted estimate, first-step redo) was not implemented; only this fixed-factor family is excluded. Equal SCF count is also not equal trajectory **Premise changed (S 8.2-8.3):** QE's rule differs (Hartree-metric E_res, ratchet, gap mapping, first-iteration redo) and the mixer metric is L2 vs QE's Hartree metric; the failure may be mixer noise sensitivity, so a different family is not excluded (backlog B04, B16). | OL "P14A"; N3; tests/reference/perf/p14a/work-tables.md |

### B3. Other solvers, backends and FFT structure (P10A to P10E, P10H)

| ID | Ticket / idea | Status | Measured effect (where) | Mechanism | Lesson | Source |
|---|---|---|---|---|---|---|
| A22 | P10A: hybrid RMM-DIIS eigensolver (Davidson warm-up/fallback) | opt-in (Davidson stays default) | Mac complete SCF (Davidson vs hybrid s): Gamma 0.208 vs 0.374; k222 1.620 vs 4.928; k222s 0.822 vs 2.294; k444s 4.837 vs 28.853; 16-atom 0.843 vs 2.681; 32 bands 11.895 vs 32.573; 64 bands 61.689 vs 229.351; max pair dE 5.9e-13 Ha. Fixed-H solver: hybrid wins only cold wide-quadratic spectra at 16/64 bands (0.0057 -> 0.0018, 0.220 -> 0.094) and reverses at 128 (7.89 vs 1.03). RSS (128 bands): davidson 70.6 MB, hybrid 77.1, pure RMM 28.7 | Hybrid needs 2-10x more H columns (k444s 17628 vs 136022); no monotone crossover with band count | Fixed-operator wins do not transfer to SCF. Revisit only with a real SCF case above 64 bands **Premise changed (V 3.3):** QE's RMM is one non-converged sweep per SCF iteration with Davidson at iterations 1, 2 and every fifth; only the converged variant was tested here. | rmm-diis-benchmark.md |
| A23 | P10A: first hybrid variant that skipped warm-up for an orthonormal start | rejected | Si Gamma 2750 applications vs Davidson 470; four-cell Si 46480 vs 2679 | Coordinate starts leave RMM far from convergence | Final policy screens true residual and warms up | rmm-diis-benchmark.md |
| A24 | P10B: optional MKL FFT and GEMM backend (feature `mkl`, sequential layer) | opt-in | dev08 T=1 (fe0492b, relayed), pure / MKL: Gamma 0.345 / 0.301, k222s 1.003 / 0.807, k222 1.850 / 1.453, k444s 7.116 / 5.493 (1.30x), si8 9.784 / 6.334 (1.54x). Later (fdce016): k444s 6.127 / 5.034 (1.22x), si8 5.801 / 6.016 (MKL slower). Kernel microbench (ADR 0019, 1 core): MKL/rustfft 3.9x (24^3), 3.0x (25^3), 2.8x (36^3), GEMM 2.5-4.2x. Energy pure vs MKL <= 4.6e-13 Ha (si8); density hashes identical across T within a backend, never across backends | Faster vendor kernels, but SCF gain is far smaller than kernel gain (call path, copies, small matrices, EVD not routed) | MKL is not uniformly faster at SCF level; SCF counts on si8 differ between backends and builds (fe0492b: pure 16, MKL 15; fdce016 and main: pure 15, MKL 16; P14B: both 15), unresolved cause (OL P10K/P10L) **Premise changed (H 1.9, 4d):** QE with MKL-DFTI does not prune FFTs, and the repository has no `[profile.release]` or target-cpu flags (dev08 builds may lack AVX2), which also bears on MKL-vs-pure comparisons. | fft-backend-comparison.md; D08 parts 1, 4, 5 |
| A25 | P10H: close the MKL call-path gap (fft3d 1.3-1.7x slower than the C reference, zgemm 750xN 1.6x vs 4.1x, zheev not routed through MKL) | superseded by P10J | Ticket text only; no numbers beyond those in the ticket | Gap partly packing/panel glue, fixed in A41 | EVD still not routed: see backlog B07 | tickets/P10H-MKL-CALL-PATH.md; R |
| A26 | P10C: fused wavefunction FFT, ordered arithmetic (compact x sticks, active-z, per-y xz planes with inverse z, potential, forward x fused); fused density; zero-weight bands skipped | adopted | Mac, five alternating pairs: Gamma 0.160 -> 0.148, k222s 0.527 -> 0.486, k222 0.992 -> 0.904; on main 9d6fef1: Gamma 0.1254 -> 0.1138, k222s 0.405 -> 0.359, k222 0.757 -> 0.671, si8 3.725 -> 3.634; density phase Gamma 0.0106 -> 0.0036, si8 0.104 -> 0.0526; per call 36^3 local 399.0 -> 373.8 us, density 343.8 -> 196.4 us. Bit-identical. Slab fixture wall rose 0.00441 -> 0.00479. dev08: not recorded separately | Fewer full-grid passes; density skips FFT for zero occupations | Gain is small for the apply path (6%), large for density (1.75x); QE-note gains were measured on other hardware | OL "P10C"; batched-fft-study.md |
| A27 | P10C: symmetric plane-fused variant (inverse z/y/x order, `PlaneFusedPlan`) | opt-in (explicit) | Energy within 1e-12 Ha, but density hashes change; Gamma/k222s/k222 dE <= 2.7e-15 Ha | Reordered axes change rounding | Kept off by default to preserve bit-identity | batched-fft-study.md |
| A28 | P10C: group 1 to 32 bands per plane pass (potential read once), many-band 3D mode | rejected (test-only) | 25^3: 168.8 -> 164.6 us/band at group 8, 176.7 at 32; 36^3: 441.6 -> 433.6 at 8, 532.0 at 32; 15x20x24 up to +33% | Multiply is compute-bound; large groups spill cache | Matches the QE-side report (grouping +-1%) | batched-fft-study.md; N2 section 2 |
| A29 | P10C: plane pruning only when active columns <= 60% | opt-in (explicit), no default threshold | 36^3 at 1/36 active: 426.8 -> 242.0 us; at 21/36: 543.4 -> 429.5; helps sparse supports only | Strided y passes less efficient than native 2D | No tuned threshold adopted | batched-fft-study.md |
| A30 | P10C: MKL fused path | not done | `FusedPlan` is the generic grid adapter under the MKL feature; no MKL fusion speedup claimed | Batched vendor stick transforms not implemented | Open: backlog B08 | OL "P10C final integration" |
| A31 | P10D: Kerker-type preconditioned Pulay mixing (`KerkerPulay`) | opt-in (default Pulay stays) | Analytic dielectric model 13/17/21 -> 3 map evaluations; complete SCF q0=0.5: fcc 24 -> 21, elongated nonlocal 22 -> 21; most others need more iterations; Si four cases slower; stretched supercell and slab fail for all strategies | Helps only long-wavelength sloshing | No universal q0; failures are reported, not hidden | scf-mixing-comparison.md |
| A32 | P10E: GPU architecture / roadmap | n/a (design only) | No code, no numbers | - | See docs/architecture/gpu-backend.md | R |

### B4. Threading and scheduling at single k point (P7PB, P10I to P10M)

| ID | Ticket / idea | Status | Measured effect (where) | Mechanism | Lesson | Source |
|---|---|---|---|---|---|---|
| A33 | P7PB: scoped band workers inside one k point (H apply, density), deterministic budget split outer=min(N,nk), inner=N/outer | superseded (mechanism) by A40 | Mac, B(N)/A(N): Gamma 1.31 / 1.46 / 1.49 (N=2/4/8), diagonalization 0.080 -> 0.039; synthetic block timing 1.90x / 3.01x / 4.55x; k cases 0.97-1.00 at fixed N (outer k uses the budget). dev08: not recorded | Independent bands; std scoped OS threads plus a mutex for completed-column copy | Thread creation is cheap on the Mac and was costly on dev08 (see A34, A44): never accept scoped-thread scaling from the Mac | OL "P7PB" |
| A34 | P7PB: packed-lane FFT threads (scoped) | rejected | Speedup 25^3: 0.77 / 0.67 / 0.44 (2/4/8 threads); 36^3: 1.11 / 0.95 / 0.82; 48^3: 0.97 / 0.94 / 0.95 | Spawn plus gather cost exceeds work per FFT (0.1-0.8 ms) | Rejects this packing and spawning strategy, not every parallel FFT; band-level parallelism is the better axis | OL "P7PB" |
| A35 | P10I: deterministic tiled GEMM (fixed 256-plane-wave tiles, ordered partial sums, borrowed faer views) with scoped threads | adopted; scheduling superseded by A40 | Mac si8 wall T=1/2/4/8 before -> after: 3.537 -> 3.754 (+6.1% slower), 3.100 -> 2.622, 2.920 -> 2.026, 2.838 -> 1.811 (1.57x at T=8); T1/T8 scaling 1.247 -> 2.072. dev08 (494bece, ticket-supplied) pure si8 T=1/4/8/16: 6.89 / 6.20 / 7.61 / 10.15 (T1/T8 = 0.91x, derived); 10 of 22 T=8 samples inside `pthread_clockjoin`; T=16 had one live thread in 12 of 29 samples | Recursive OS thread spawn per GEMM and per H block; fine on Mac, fatal on dev08 | A 2x Mac gain became 0.91x on dev08. Mac scaling is not evidence | OL "P10I", "P10J" |
| A36 | P10I: blocked projection against the established basis (two GEMM passes, normalise between passes) | adopted | Architect-supplied Mac profile before (about 1550 samples, si8 T=1): GEMM microkernel 553 (36%), subtract_projection 264 (17%), FFT about 270 (18%), nonlocal add 93, fused local 60. dev08 MKL T=1: projection 1.354 -> 0.443 s | Replaces scalar projection (single-threaded, four lanes) by GEMM | Interleave normalisation and reprojection; see A37 | OL "P10I", "P10J" |
| A37 | P10I: project twice before normalising | rejected | RMM 64-band quadratic-spectrum Davidson test failed at 800 iterations (residual 5.93e-9 vs threshold 1e-9) | Near-rank-deficient incoming block amplifies round-off | Reproject only after normalising | OL "P10I" |
| A38 | P10I: pack Rust operands into contiguous panels | rejected | 5-8% many-k overhead; borrowed views removed it; MKL still packs (adapter needs contiguous operands) | Extra copies | Prefer borrowed views | OL "P10I" |
| A39 | P10I: scheduling guard (panels below 65,536 complex multiply-adds run inline) | adopted, superseded by A41 guard | k222s T=8 median 0.1291 before the guard | Worker launch cost on tiny panels | Shape-based guard, independent of T | OL "P10I" |
| A40 | P10J: SCF-owned persistent Rayon pool shared by k solves, GEMM, projection, Hamiltonian bands, density; whole k-solve and audit phases installed once; fixed tiles as individual stealable jobs; no mutex or copy in H output | adopted | Mac si8 before -> final: T=1 3.783 -> 3.743, T=2 2.556 -> 2.471, T=4 2.032 -> 1.943, T=8 1.826 -> 1.771 (alternating 7-pair T=8: 1.830 -> 1.761, faster in 7/7); T1/T8 = 2.11x. First review version regressed Mac T=8 by 4.5%; fixed by per-tile jobs and resident phases. dev08: acceptance run was not isolated; the next relayed state (fdce016) gives pure si8 T=1 5.801, T=8 4.151 (1.40x) | Removes thread creation; but regions are small (mean 294-327 us) and utilisation low (A47) | Local T=1 changes within +-2.5%: needs paired runs to resolve | OL "P10J" |
| A41 | P10J: MKL uses one sequential logical product per GEMM (no panel packing, copies, assembly, per-panel scans) | adopted | dev08 (ticket-supplied) MKL T=1 profile before: tile glue 0.325 -> 3.874 s, memcpy/memset 0.433 -> 1.107 s, kernel 1.137 -> 0.664 s. After: not isolated; MKL si8 T=1 6.334 (fe0492b) vs 6.016 (fdce016) | Glue overhead was larger than the kernel | MKL intra-GEMM parallelism is given up; k solves and H bands still parallel | OL "P10J"; D08 parts 1, 4 |
| A42 | P10J: dispatch threshold 262,144 complex multiply-adds and at least 4 panels, independent of T | adopted | Mac: inline 14.9 us vs forced pool T=8 21.5 us at 1025x8; crossover near 262,400 (53.4 vs 31.6 us); dispatch ingress 4-9 us external vs 0.2-1.6 us resident (T=2..16) | Dispatch cost versus work | Mac-measured; not a portable optimum | OL "P10J" |
| A43 | P10J: cache-blocked Nonlocal traversal (fixed 256-row output blocks, precomputed weights) | rejected on dev08 (reverted in 9fcb356, bit-identical) | Mac microbench 8 atoms x 4097 rows x 3 channels 197.4 -> 117.5 us; 20 atoms x 16385 x 9: 2883 -> 2925 (neutral). dev08: `NonlocalOperator::add` 0.474 s (17 of 227 samples) at fe0492b vs 0.929 s after P10I (about 2.0x; the first claim of 2.9x used 0.323 s from only 6 of 117 samples); the revert gave 0.800 s (239 samples) but wall unchanged beyond noise (MKL T=1 6.612 -> 6.589, pure T=8 4.567 -> 4.540) | Mac shows shape-dependent mixed results; remote function time doubled | Microbench, Mac and low-count profiles all misled. See C **Premise changed (H 3, 4):** the traversal tried was over output rows, not over bands; the 1.52 MB projector matrix does not fit dev08's 1 MB L2 (hypothesis, not proof); band-blocked traversal is a different change (backlog B07). | OL P10J, correction note |
| A44 | P10J: serial intra-k density as an interim version | superseded | Mac si8 T=8 density 0.016815 -> 0.051618 s; restored parallel on the pool after architect authorisation (density 0.017271 vs 0.017021) | Removing scoped threads exposed serial density | Interim versions are not acceptance | OL "P10J resumed" |
| A45 | P10K: validation reuse (single finite scan per new block), paired basis/H rotation in one task tree, Gamma uses the resident ordered-solve path (`threads <= 1` is the only serial branch) | adopted | Mac si8 T=1/8/12: 3.775 -> 3.666 (1.030x), 1.727 -> 1.647 (1.048x), 1.707 -> 1.590 (1.073x); ablation: validation alone -0.047 s (T=1) / -0.027 s (T=8); resident alone no Mac win; installed regions 2463 -> 2055 (-1 per Davidson iteration), external ingress inside solves 2463 -> 0. dev08 (supplied before): main-thread latch wait 32%, multiply_into 18.8%, subtract_projection 13.8%; 8044 FUTEX_WAKE, 1686 blocking waits | Gamma had bypassed the resident path (P10J claim was wrong) | Check the call path of the single-k branch | OL "P10K" |
| A46 | P10K: limit recursive leaves to min(T, panels) | rejected | si8 T=8 1.815 vs 1.727 (about 5% slower) | Fewer stealable pieces worsen tails; Rayon workers wake anyway | Do not cap pieces | OL "P10K" |
| A47 | P10K: finite-scan as one Boolean reduction; spin/yield policy | rejected | Scan microbench 1024-524288 elements: at most 2.6% gain (0.328 vs 0.321 us to 176.5 vs 171.8 us); spin policy not adopted | Too small to explain the dev08 gap | Do not burn cores on spin | OL "P10K" |
| A48 | P10L: right-looking first-pass MGS for tall blocks (>= 2048 rows), overlap of the dependent second pass, per-band stealable H pieces (Gamma), reserved concatenation | adopted | Mac si8 T=8 1.664 -> 1.580 (1.053x, 5.0%), T=1 3.687 -> 3.668, k444s T=8 0.577 -> 0.571; 20-atom skew synthetic T=8 5.0% (local-only) and 5.6% (projectors) faster. dev08 si8 T=8 worker utilisation (relayed): 19.5% before (293 samples), 22.5% (P10L), 23.0% at fdce016 (250 samples, mean 1.84 of 8 busy) vs 33% needed for 2x and 55% for 3x | Dependency chain in MGS; small stealable pieces | Utilisation goal NOT met. Do not read the 5% as a scaling fix | OL "P10L"; D08 part 4 |
| A49 | P10L: finer GEMM subtiles (4-column Gram, 64-row rotation, 8-column) and lower dispatch cutoff 65,536 | rejected | Fine tiles: si8 T=8 1.659 -> 1.518 but Davidson count 340+68 -> 293+69 (workload changed), T=1 3.679 -> 3.820; 8-column: about 2% only; cutoff change: gains and losses vary with T and shape | Smaller faer calls change rounding, hence Davidson trajectory | Finer GEMM needs an accepted rounding policy (see C) | OL "P10L" |
| A50 | P10L: dynamic H pieces under outer-k parallelism | rejected | First modulo-scratch prototype regressed k222s | Scratch slot collisions when inner budget < pool | Static groups for outer-k | OL "P10L" |
| A51 | P10M: storage and duplicate-GEMM removal in Davidson (reserved live-prefix basis, in-place truncation, new-column append, saved restart product reused) | adopted | si8 explicit copied bytes 3083884288 -> 736959264 over 408 iterations; GEMM 2mnk 21834221954 -> 21636961594 (14 duplicate restart GEMMs removed); Mac 30 pairs si8: T=1 3.667 -> 3.604 (1.018x), T=8 1.627 -> 1.585 (1.026x), T=12 1.580 -> 1.539; k444s T=1/8/12 1.0018/1.0047/1.0099. First version regressed k444s T=1 (2.549 -> 2.569, 28 of 30 pairs slower); fixed by removing zero-fill and duplicate restart GEMM | Less copying and recomputation | Run the 30-pair gate on a many-k case too | OL "P10M" |
| A52 | P10M: CholeskyQR2 block orthonormalisation with 64-row projection panels (47 pieces) | rejected (`--cfg feature="block-qr"` only) | Mac exploratory: si8 Davidson 408 -> 374, energy within 2.4e-13 Ha, T=8 1.574 -> 1.447 (3 runs); but 20-atom skewed synthetic iterations 233 -> 329 (64-row) or 311 (256-row) | Conditioning of normal equations changes Davidson trajectory | Passed 5 Si cases, failed a different cell: confirms the generality rule **Premise changed (V 7):** CholeskyQR2 yields the same Q as MGS up to rounding, so the skewed-cell regression was rounding amplification or a different rank decision at the pivot test; no record shows which. | OL "P10M" |

### B5. Recent work accounting (P14A, P14B)

| ID | Ticket / idea | Status | Measured effect (where) | Mechanism | Lesson | Source |
|---|---|---|---|---|---|---|
| A53 | P14A: work accounting (steps, H blocks, vectors, restarts, buffer origin) | adopted as instrumentation only | Local si8: 15 SCF, 340 ordinary + 68 certificate steps (408), 410 H block calls, 2416 H vectors, 14 restarts; k444s 1686 steps, 1891 blocks, 11181 vectors. Relayed dev08: QE 83 steps / 98 H blocks vs BAANDI 290 steps / 2070 vectors | Units differ (QE h_psi per block, ours per vector) | Compare steps, not vectors; QE ethr is a Ritz-value change in Ry, ours a residual in Ha **Premise changed (S 5.4-5.5, V 2.5):** QE's definition at conv_thr 1e-8 is weaker (Ritz-change test, up to 5 unconverged roots, no certificate), so part of the 290 vs 83 step gap is definition; QE steps equal H blocks minus (SCF iterations + 1). **Update 2026-11-06:** P15A replay: 283 -> 39 rows at identical tolerances; QE-ladder 89 rows vs relayed 83: the gap is mostly the stopping rule. | OL "P14A"; N3 |
| A54 | P10G general benchmark suite, P11A/P11J observer-overhead measurement | n/a (measurement infrastructure, no optimization) | Suite: 19 cases, synthetic and Si, three persistent non-convergers (si-elongated, syn-s3, syn-z14). P11A: Gamma 0.112248 vs 0.111934 s, four k 0.359221 vs 0.359967 s (no resolvable regression). P11J: see C.6 | - | Use the suite for the generality gate | benchmark-suite.md; OL P11A/P11J; D08 part 3 |
| A55 | P14B: incremental reduced Hamiltonian (keep prefix, two rectangular products) | adopted (default, pure and MKL) | si8 2mnk 21636961594 -> 9489226838 (-56.14%), k444s 3300404692 -> 2258007134 (-31.58%); skewed 20-atom -14.75%, elongated 4-atom -18.84%. dev08 (7fe2e9a vs 65036b7, relayed): si8 MKL T=1/4/8/16 6.044 -> 4.814, 4.794 -> 3.738, 4.798 -> 3.718, 4.984 -> 3.899 (1.26-1.29x); si8 pure T=1/4/8/16 5.800 -> 5.463 (1.06x), 4.010 -> 3.877, 4.135 -> 4.026, 4.578 -> 4.510 (1.015-1.034x); four k-point cases 0.958-1.032 (inside spread). Work: si8 MKL 338 -> 302 steps and 16 -> 15 SCF; GEMM kernel samples only 1.286 -> 1.111 s (0.175 s) | Reduced matrix was recomputed fully (46% to 57% of summed GEMM time on k444s and si8); the relay attributes most of the MKL gain to less work (one SCF iteration and 36 Davidson steps fewer), not faster GEMM; the SCF-count change coincides with changed GEMM rounding | Mac predicted 1.42x for pure si8 T=1; dev08 gave 1.06x. Rounding changes fixed points: 4 toy cases changed iteration count (14 -> 16, 8 -> 10, 13 -> 14, 17 -> 16), 3 failures persisted, 12 Bohr cube H vectors rose 1996 -> 2107 **Premise changed (dev08 relay, V 6.2):** the 56% explicit GEMM flop cut gave only 1.06x on pure dev08 (1.26-1.29x MKL), so flops do not map to time for pure Rust; one-rectangle reduced H (V-5) reverses the 'no conjugation fill' choice (backlog B01). | acceptance.md; D08 part 5; N4 |
| A56 | P14B: lock converged roots / skip Ritz and residual for converged roots | gated | Safe removable work under the current rule: zero (a root below tolerance can fail after the next rotation); ceiling = Ritz/residual 5386241040 of 9489226838 2mnk on si8 (56.8%), 946782528 of 2258007134 on k444s | Unlocked stopping rule needs all residuals | Needs a stability argument and measurement: backlog B01 **Premise changed (V 2.1-2.3, 7):** the zero bound holds only for unlocked roots; QE decides convergence from the small problem and builds no N-sized object for converged roots; hard locking with a cluster rule and certificate unlock is an open route (backlog B01). **Update 2026-11-06:** locking was tried (P24A-1, P24A-2): MACs -6.0% / -12.5% only and failed gates; the ceiling is not realisable by prefix locking. | N4; p14b/adopted/acceptance.md |

Row count: 56 (A01 to A56).

### B6. Rows added 2026-11-06 (symmetry, forces/stress)

| ID | Ticket / idea | Status | Measured effect (where) | Mechanism | Lesson | Source |
|---|---|---|---|---|---|---|
| A63 | P19A: crystal symmetry (k reduction, density symmetrisation), `symmetry.mode = "auto"` | opt-in (default off) | Mac observation, shifted 4x4x4 Si: 64 -> 32 -> 10 k points, k444s T=1 3.328 -> 1.106 s (3.01x), applications 12973 -> 4040, full-vs-reduced dE 1.24e-14 Ha; QE symmetry-on count unconfirmed (no stored QE output) | Fewer k points, density projector | Needs commensurate FFT grid; dev08 and QE symmetry-on bundle pending | docs/validation/symmetry.md |
| A64 | P20A-P23A: forces, stress, relaxation, cell relaxation | adopted (features) | Mac observation, one sample: forces 28.3% of SCF+force on si8, 10.9% on k444s; stress 58.4% of SCF+pass on si8, 44.6% on k444s (Pulay scans extra) | New phases, not optimised | Force and stress cost does not shrink with the standard profile (uncertified); add them to the benchmark set | tests/reference/forces/cost.json; stress/cost.json |

Rows A63-A64 were added 2026-11-06 (the P15A, P15C, P17A, P24A and native-flag results are recorded in the append-only sections below by their own tickets; A57-A59 numbers are used twice there).


Remarks beginning "Premise changed" (2026-11-05) were added after the four QE study notes (G, H, S, V in docs/research/qe-study/); the rows themselves are unchanged. Plan: [exceed-qe-plan.md](exceed-qe-plan.md).

## C. Measurement methodology and pitfalls (learned the hard way)

1. **Mac and dev08 behave differently; the Mac is not an acceptance machine.** P10I measured 1.57x at T=8 and 2.07x T1/T8 scaling on the Mac and 0.91x T1/T8 (6.89 s to 7.61 s) on dev08: thread creation is cheap on macOS and costly on Linux with 80 logical CPUs (A35). P10A to P10M local timings are "observations". Mac wins predicted for dev08 failed twice more: P14B pure si8 T=1 1.42x on the Mac became 1.06x on dev08 (A55), and the blocked Nonlocal traversal helped Mac microbenchmarks and doubled the function on dev08 (A43). Acceptance numbers come from dev08 through the orchestrator.
2. **Three-point comparison.** When a change lands after other merges, measure the baseline commit, the previous main and the branch (the P14B relay used fdce016, 65036b7 and 7fe2e9a with separate target directories and sha256 per binary) so an effect is attributed to the change and not to the 26 commits in between.
3. **Profile sample counts.** A ratio from fewer than 30 samples is directional only. The 2.9x "Nonlocal regression" came from a 6-of-117-sample reference; the corrected value from 17 of 227 samples was 2.0x. Rule in the log: never quote a profile-derived ratio as firm below about 30 samples. mkl T=1 169-sample profile: figures above 30 samples firm, others directional (D08 part 5).
4. **Build hash embeds the git commit.** `crates/dft-cli/build.rs` writes `DFT_GIT_COMMIT` into the binary, so a new commit gives a new binary sha256 with identical code (P11J p100 vs p300). Record the source hashes next to the binary hash; do not read a hash change as a code change.
5. **Shared target directory contamination.** `scripts/bench_runtime.py` built candidate and reference revisions into one `target`, so both `dft-cli` packages overwrote each other's artifacts (fixed in 8be0c05 with `target/bench-reference`). Always build each commit under test into its own target directory and record exit codes.
6. **Pairs needed to resolve small overheads.** 20 pairs gave 95% intervals up to 3.9 percentage points and the null control (classic) was inconclusive in 11 of 18 combinations in each of three runs (met/inconclusive 8/10, 10/8, 9/9). 100 pairs: 15 of 18 met. 300 pairs: 18 of 18 met, classic -0.207 to +0.214%, observed -0.048 to +0.388%, agent +0.012 to +0.552% (relayed, D08 part 3). Local Mac (load 11.6 to 14.5) gave not_met 7, inconclusive 10, met 1 over 18 cells. Use 30 alternating pairs for a 1-3% wall claim (P10M) and 7 alternating pairs only to settle a single outlier (P10K k444s T=2: 1.0068x; P10L k222s T=2: 0.18%).
7. **QE timer resolution is 0.01 s.** Entries below about 0.3 s are not valid ratios (one tick up to 20% at 0.05-0.27 s). Only k444s and si8 qualify; the si8 16-core 0.3000 s point is at the boundary. QE internal WALL excludes mpirun start-up; BAANDI wall includes process start; `/usr/bin/time` and `bench_scf.py` differ by up to about 5% (k444s 94.3 vs 98.98 s).
8. **QE ratios depend on conv_thr.** k444s N=8: 2.21 (1e-8) vs 1.61 (1e-12); si8 N=1: 8.48 vs 6.36. QE energies are equal to 8 digits for both; 1e-12 only adds iterations (QE 14 vs 20 SCF on si8). State the conv_thr with every ratio. QE spread above 5% remained in 6 of 80 cases after 9 repeats; a batch with a loaded node (load 1.9) gave 30.8% and 18.6% spreads and was not used.
9. **Relayed versus reproduced numbers.** dev08 numbers reach this repository as text; `dev08_final.json` (sha256 8fb74fb0...) never transferred (base64 corruption). Every row here marks relayed numbers; do not write relayed numbers into README as measured here (README uses D08 part 4 only, with its limits). Do not treat an unsynchronised thread-count, commit or backend as comparable.
10. **Density bit-identity versus energy gates.** Overhead-only changes must keep density bytes identical. Algorithm or rounding changes (A05, A35, A55) are accepted on energy <= 1e-12 Ha plus certificate, and density identity across thread budgets within one backend; density hashes across backends always differ and historical hashes differ from main. Do not rewrite baselines; record same-HEAD pairs. Harness exit 5 is a hash/fingerprint mismatch (thread count is part of the input fingerprint), not necessarily a physics failure.
11. **Fixed-point multiplicity makes toy iteration counts unreliable.** In syn-s2 (7 plane waves) a displacement of 1e-14 reaches different converged energies (0.6304095283, 0.6348260419, 0.6354522014 Ha) and the iteration count moves 13 to 24; so pure vs MKL can converge to different solutions. The 1e-12 Ha cross-backend gate holds only for unique solutions (docs/validation/scf-fixed-point-multiplicity.md; narrative-only evidence). Consequence: the architect superseded identical-count rules for toy fixtures in P14B (realistic fixtures keep SCF count and energy within 1e-12 Ha; converged toys keep energy tolerance; iteration changes reported).
12. **Generality.** Test cubic, skewed, elongated and slab cells, 2 to 20 atoms, Gamma and k meshes, synthetic projector structures. Things that passed Si and failed elsewhere: CholeskyQR2 (A52), extra bands (A19), tolerance factors (A21), fine GEMM tiles (A49), the slab wall-time rise under fused FFT (A26).
13. **Sandbox cannot profile.** `sample` returned exit 255 in nearly every implementer run, so implementers rely on architect or dev08 profiles. Counter builds (`worker-analysis`, P14A/P14B accounting) perturb timing and are never used for acceptance; use a separate uninstrumented binary. Region busy fractions are lower bounds; pool residency is not parallel fraction.
14. **Mixed units.** QE steps are `h_psi` per block (nbnd), BAANDI counts vector applications; compare Davidson steps (QE 83 vs 290) not vectors. QE `ethr` is a Ritz-value change in Ry, ours a residual norm in Ha. "Equal SCF count" is weaker than equal trajectory.
15. **Instrumentation and microbench traps.** Adding a case to a standalone benchmark changed compiler-generated timings while the CLI binary was byte-identical (P10J); a 47x kernel gain (A03) did not appear in SCF; stale per-binary baselines (`baseline.json` vs `baseline-after-o2.json`, `-o4`, `-o4c`) must match the step under test.

## D. Bottleneck map today (as of fdce016 plus P14B, 2026-11-05)

Sources: D08 parts 4 and 5 (relayed dev08), OL P10I/P10L/P10M, N3/N4, tests/reference/perf/p14a/work-tables.md, p14b/gemm-tables.md. "Derived" = arithmetic on recorded numbers.

**Where the wall time goes**

| Item | si8 Gamma (8 atoms, 20 bands, 2969 PW, 36^3) | k444s (32 k, 8 bands, 25^3) |
|---|---|---|
| Davidson share of SCF (Mac, ordinary phases, P14A) | 93.0% to 97.2% per SCF iteration; density 0.9-2.8%, energy audit 1.0-2.9%, potential 0.4-1.0%, mixing 0.0-1.3% | 80.2% to 89.5% per SCF iteration; energy audit 6.6-12.9%, density 3.7-6.2% |
| Mac si8 T=1 phases (P10I final): wall 3.754, diag 3.472, density 0.051, potential 0.019, mixing 0.013, energy audit 0.059 | density scales with T (0.017 at T=8); potential, mixing and audit stay 0.019/0.014/0.059 at every T (serial) | dev08 k444s T=16: engine 1.211, diagonalization 0.839 (after P7PA) |
| Mac profile before P10I (about 1550 samples, si8 T=1, architect-supplied, stale) | GEMM microkernel 36%, subtract_projection 17%, FFT about 18%, nonlocal add 6%, fused local 4%, EVD 1.7% | not recorded |
| GEMM site shares (P14B accounting, Mac, si8, before A55; share of summed GEMM logical-call time, 1.78 s in total) | reduced Hamiltonian 56.8% (ordinary 48.3 + certificate 8.4), Ritz/residual 27.2%, four projection sites 16.0% (derived by summing the table) | reduced Hamiltonian 46.4%, Ritz/residual 24.9%, projection 28.8% |
| GEMM work after A55 (2mnk, si8) | total 9489226838; Ritz/residual 5386241040 (56.8%); second projection pass 1296609804 (13.7%) | total 2258007134; Ritz/residual 946782528 (41.9%) |

**Davidson steps and per-step cost**

| Quantity | QE 7.6 (dev08, relayed) | BAANDI |
|---|---|---|
| si8 N=1 Davidson steps, conv_thr 1e-8 | 83 steps (98 H block calls), 14 SCF, 5.93 steps per SCF | dev08 (fdce016, T=1): 290 steps, 2070 H vectors, 15 SCF, 19.33 per SCF; P14B build 279 steps / 2044 vectors; Mac local: 340 + 68 certificate = 408 steps, 410 H blocks, 2416 vectors, 14 restarts |
| si8 N=1 wall | 0.684 s (PWSCF; cegterg 0.45 s over 14 calls; h_psi 0.36 s over 98 calls) | pure 5.801 s, MKL 6.016 s (P14B: pure 5.463, MKL 4.814) |
| Per-step wall (derived: total wall / steps) | 0.684 / 83 = 8.2 ms (cegterg alone 5.4 ms) | 5.801 / 290 = 20.0 ms; ratio 2.4x per step times 3.49x more steps = 8.5x observed ratio (pure / QE = 8.48) |
| k444s steps | about 25 per k point over 6 SCF (4.7 + 3.8 + 1.2 + 3.5 + 4.5 + 7.0), 1010 H blocks | 1686 ordinary steps = 52.7 per k point (derived), 1891 H blocks incl. certificate |
| Tail effect (P14A si8, Mac) | QE steps spike to 12, 11, 20 right after ethr tightens, else 1-2 | BAANDI steps per SCF 9, 10, 10, 24, 11, 13, 15, 16, 33, 34, 18, 52, 45, 23, 27: SCF 4, 9, 10, 12, 13 hold 188 of 340 steps (55%, derived) with 3.7-6.4 H vectors per step (e.g. SCF 12: 52 steps for 192 vectors); one late root per SCF (early/late 19/1) |
| First SCF iteration share of H vectors | not recorded | si8 104 of 2416 (4.3%, derived); k444s 1410 of 11181 (12.6%, derived) |

**Parallel scaling (dev08, relayed)**

- si8 pure T=1 5.801 s to T=8 4.151 s: 1.40x (P14B build 5.463 to 4.026: 1.36x, derived). Mac T1/T8 is 2.32x (P10L), so Mac scaling overstates dev08.
- Worker utilisation si8 pure T=8 (12 runs, 250 samples): busy workers 0: 9.6%, 1: 29.6%, 2: 34.8%, 3: 20.4%, 4: 4.4%, 5: 1.2%, 6 or more: 0.0%; mean 1.84 of 8 busy workers (23.0%; P10L 22.5%); max/mean imbalance 3.01 (P10L 3.29). Serial part 9.2% (0.382 s), parallel region 3.769 s at 26.2% utilisation; at 100% utilisation 1.369 s (4.24x); actual 4.151 s. Targets: 33% utilisation for 2x, 55% for 3x (P10L ticket text).
- Mac si8 T=8 (P10L, instrumented, SCF 1.411 s): scheduled parallel regions 0.752 s; between logical regions 0.556 s; Davidson wall 1.302 s of which 0.550 s (42.2%) outside scheduled regions, 0.459 s worker means, 0.200 s max-minus-mean imbalance, 0.093 s unaccounted. Inclusive phases: Hamiltonian 0.310, Gram 0.288, projection 0.260, rotation pair 0.140, EVD 0.139 (sequential), MGS second pass 0.074 and first pass 0.063. Worker 0 busy 0.685, others 0.28-0.35. Imbalance max/mean: Gram 1.29 median (worst 2.70), MGS first pass 2.24 (worst 5.55), Hamiltonian 1.20.
- Dispatch is not the bottleneck on the Mac: resident empty region 0.19 us at T=8, external 7 us; mean installed region 327 us; 2055 regions per solve (137 per SCF iteration) after A45.
- Pure vs MKL on dev08: MKL is not uniformly faster (si8 T=1 6.016 vs 5.801 at fdce016, MKL faster for k444s 5.034 vs 6.127). After A55: si8 MKL 4.814 vs pure 5.463 at T=1.
- k-point parallelism (fe0492b): k444s pure T=1 7.116, T=16 1.003 (7.1x), MKL 5.493 to 0.853 (6.4x); Gamma is flat (0.345 at T=1 to 0.361 at T=16), k222s stops scaling beyond T=4 (0.428 to 0.458 at T=16) and k222 beyond T=8 (0.498 to 0.528).
- QE ratios (BAANDI / QE internal time, fdce016, relayed): k444s 1e-8 N=8 2.21 pure, 1.99 MKL; N=16 2.17 / 2.05; 1e-12 N=8 1.61 / 1.44; si8 N=1 8.48 / 8.80 (1e-8), 6.36 / 6.60 (1e-12); si8 N=16 at the timer limit 15.2 / 16.8. After P14B: si8 N=1 1e-8 7.99 / 7.04, 1e-12 5.99 / 5.28.

**Update 2026-11-05 (notes H, V).** QE h_psi costs about 3.7 ms per Davidson step against 3.5-4.5 ms of BAANDI H time per step: the H application is not the per-step gap. The gap is Davidson algebra (V cost model: about 5600N vs 1100N multiply-adds per step, matching the P14B GEMM total within 20%) and the vector count (about 2.2x); Davidson capacity is `max(48, 6*bands)` = 120 for 20 bands. Full decomposition: exceed-qe-plan.md section a.

**Reading the map.** (1) The step count is 3.5x QE's on si8 and the per-step cost 2.4x; neither of the five tolerance/band/preconditioner ideas tried so far reduced steps without raising SCF count (A13-A15, A19, A21). (2) After A55, half of the remaining GEMM work is Ritz/residual formation for all roots, even in tail steps that carry about four vectors. (3) Parallel utilisation (23%) is the dominant loss on dev08 at T>=4, caused by small regions with serial dependency chains (MGS, EVD, residual norms), not by dispatch latency. (4) Non-Davidson phases (potential, mixing, energy audit) do not scale at Gamma: about 0.09 s of 1.8 s at Mac T=8 (derived, 5%). (5) k-point parallelism is already good for 32 k points and saturates for 4 and 8.


## P15A measurement additions (2026-11-05; append-only)

The historical row count above describes A01–A56. These new rows extend it to
A59 without rewriting prior attempts. Derived from solver note V-1 and SCF
note S-E0–E2. Evidence: `tests/reference/perf/p15a/`,
`docs/research/qe-study/replay-stopping-results.md`.

| ID | Experiment | Status | Measured effect | Lesson / limits |
| --- | --- | --- | --- | --- |
| A57 | P15A frozen-H stopping and tolerance replay | measurement only | si8 native 283 rows / 2063 H vectors; Ritz predicate with identical per-root numeric tolerances 39 / 749; QE-empty policy 38 / 729; Hartree ratchet on native history 54 / 1000. k444s 1686 / 10290 to 740 / 5631 under ratcheted Ritz. | Stopping strength can explain the size of the step gap on these fixed histories. No actual SCF policy changed; no equal-accuracy speedup. QE per-root rate traces remain missing. A15/A21 are not repeated or overturned. |
| A58 | P15A E0 and offline schedule/first-redo diagnostic | measurement only | E_res < 5e-9 Ha first at si8 11, k444s 7, skew 4, elongated 5, 16-atom Gamma 7; no crossing in slab/16-atom mesh at 100 iterations. First-redo condition fires on skew and 16-atom mesh. | Conditional stored-history ratchet is not a new SCF trajectory. Actual redo work requires its changed start/residual. Two new cases exhaust identically before/after. H6 remains unproved. |
| A59 | P15A independent dense spectral-completeness oracle | measurement only; follow-up needed | si8 initial all-root 1e-10 Ha residual replay misses a low root: occupied error 0.110723486 Ha. Initial synthetic all-root discrepancies 2.958860388 / 0.380087693 / 8.225605511 / 0.523220374 Ha (skew/elongated/16-atom/slab). Native final skew occupied error 1.514951753 Ha versus its frozen dense spectrum. | True residuals do not certify lowest-root completeness. Baseline outputs are bit-identical, so this is an existing finding. Resolve completeness before claiming equal accuracy or adopting locking/schedule changes. No repair in this ticket. |
### B6. Added 2026-11-05 after the table above

| ID | Ticket / idea | Status | Measured effect (where) | Mechanism | Lesson | Source |
| --- | --- | --- | --- | --- | --- | --- |
| A57 | Build flags: `-C target-cpu=native`, LTO fat with codegen-units 1, both (no source change; commit 84bfbf3) | native: adopted for benchmarking via scripts/build_native.sh, not the portable default; LTO: rejected | dev08 relayed: native 1.12 to 1.40x on pure Rust (si8 T=8 1.399, si8 T=1 1.316, k444s T=1 1.229), 1.05 to 1.07x on MKL builds; LTO alone 1.008 to 1.057; native plus LTO within 1 to 4 percent of native; energy and density hashes bit-identical in all 24 comparisons | hand-written loops (subtract_projection above all: 0.822 to 0.224 s, 46 percent of the saving) compile for 128-bit x86-64 without the flag; library kernels (rustfft, faer, MKL) already dispatch at run time; no cross-crate inlining bottleneck | check build flags before any hand-written-loop optimization; native makes pure Rust beat the default MKL build on si8 (4.122 versus 4.814 s T=1); portability policy needed for release builds (benchmarks record both) | D08 part 7 |

## P15C spectral completeness (2026-11-05; append-only)

| ID | Experiment | Status | Measured effect | Lesson / limits |
| --- | --- | --- | --- | --- |
| P15C-1 | Final-state full-space spectrum | measurement | si8 occupied max 1.39e-14 Ha; sampled k444s 9.55e-15 Ha; skew final 1.51495 Ha and occupied band-sum excess 27.9580 Ha | Realistic samples pass; skew legacy convergence is incomplete. Sampling is not all-k proof. |
| P15C-2 | Cold 2*bands independent certificate solve | rejected as default | extra H si8 1035 (45.70%), k444s 11463 (102.52%) | Recovers measured frozen-H roots but expensive. |
| P15C-3 | One-root orthogonal-complement screen, full-space small fallback | experimental, not adopted | extra H si8 89 (3.93%), k444s 2618 (23.41%); complete realistic energy/density bits and 15/10 iterations unchanged | Large-operator detection heuristic, not inertia certificate. Skew recovery exhausts 100 SCF steps; target unmet. |
| P15C-4 | Coordinate plus random initial space of width 2*bands with guard | experimental, not adopted | si8/k444s final energy changes 2.70e-13 / 3.55e-15 Ha; 15/10 steps unchanged; synthetic Gamma 10→9; skew still exhausts | No default accepted; not a repeat of capacity-only optimization. Mac walls are observations, no dev08 measurement. |

| P15C-5 | Runtime final-certificate-only spectral audit, all reported k/spin states | explicit option accepted; default **off** | si8 +90/2265 H (3.9735%), k444s +2618/11181 (23.4147%); 15/10 steps unchanged, energy change 0 Ha and identical density | Target <3% unmet even at final pass only. No early trigger added. Cumulative retries/work/errors retained; failed repairs reject certification. Mac wall observations only; dev08 request in validation/spectral-dev08-request.md. |
| P15C-6 | Dense eigensolver at every skew SCF step, same initial density/mixing | diagnostic, no defaults tuned | 100 steps exhausted; max eigen residual 2.889e-14 Ha; final density residual 2.079e-3 vs 1e-7 target; initial occupied gap 1.776e-15 Ha | Corrected SCF difficulty is not solely rescue history. Soft synthetic near-degenerate integer-occupation stress case, not evidence against realistic insulating final samples. Final audit catches legacy 27.958 Ha band-sum excess, retries once and refuses false convergence. |

## P17A preparation (2026-11-05; append-only)

| ID | Experiment | Status | Measured effect | Lesson / limits |
| --- | --- | --- | --- | --- |
| P17A-1 | Two-level strict references across nine fixtures | partial measurement | Six plateaus <=1e-12 Ha; strict default energy differences si8 4.97380e-13 Ha, k444s -1.77636e-15 Ha, Si 3.37508e-14 Ha, MgO 3.83693e-13 Ha | Skew tight/tighter diverge (exit 4); slab and mesh16 exhaust (exit 3). Accuracy stays null without a usable reference. `perf/p17a/strict-preparation-report.json`; no spectral-completeness claim. |
| P17A-2 | E1 calibration evidence inventory; standard/E4 readiness | blocked, no algorithm tried | No raw requested QE si8/k444s 1e-8/1e-14 pair; four requested standard/Hartree configurations rejected by unchanged CLI (exit 1) | No default or speedup adopted; profile and metric implementation still required. Do not substitute rounded relay energies for delivered error. A15/A21 are neither repeated nor overturned. |

## P17A implementation after architect follow-up (2026-11-05)

| ID | Experiment | Status | Measured effect | Lesson / limits |
| --- | --- | --- | --- | --- |
| P17A-3 | Ratcheted Ritz standard profile, 3 thresholds | opt-in implemented; default strict | Frozen 5e-9 Ha: si8 Davidson 345->69, H 2265->742, energy error 1.00089e-9 Ha; k444s 1827->866, H 11181->5864, error 1.33280e-10 Ha. Six reference-backed systems. | Different explicitly uncertified contract, not an A15/A21 strict-path loosening. Stored QE estimates calibrate scale; equal printed energies do not prove zero error. No dev08 speed acceptance. |
| P17A-4 | Standard generality stress fixtures | negative results retained | Elongated H 171->172; skew standard 100 SCF / 4168 H exhausted; slab 13 SCF / 297 H diverged; mesh16 100 SCF / 8182 H exhausted. | Three cases lack tight references: work-only, no delivered accuracy. Skew strict convergence is already spectrally suspect. |
| P17A-5 | Independent Hartree Pulay metric, strict E4 | experimental switch, L2 default retained | H si8 2265->2190, k444s 11181->11206, MgO 12202->12382; MgO energy delta -1.40687e-12 Ha exceeds 1e-12 target. | Not a uniformly successful optimization. SCF counts unchanged on six reference cases; mesh16 strict/Hartree diverges. |
| P17A-6 | Strict preservation, standard determinism and restart | verified | Five strict cases at T=1/T=8: 0 Ha delta and identical density versus pre-feature. Nine standard cases x two metrics: 18 T=1/T=8 pairs identical including failures. Exact restart density/orbital bytes match for both metrics. | Timers excluded; standard is never certified. Evidence: `perf/p17a/final-verification.json`, `thread-identity.json`. |


## P24A locking and Hermitian fill (2026-11-06)

Derived from `qe-solvers-and-subspace.md` hypotheses 2/5, notes only. Baseline
46dfeddd370f5aaf37c898aca41735269c36cc05; uncommitted experiment. These rows
change A56's unlocked-root premise, not the A08/A13/A15/A21 algorithms.
[Full results and commands](../research/qe-study/p24a-locking-results.md).

| ID | Experiment | Status | Measured effect | Lesson / limits |
|---|---|---|---|---|
| P24A-1 | Whole-cluster lowest-prefix hard locking, 1/sqrt(m) margin, final rotation/certificate and auto-unlock; one-rectangle reduced H | implemented, experimental default off | Si8 MAC 4,744,613,419 -> 4,459,794,280 (-6.00%), H 2265 -> 2809, steps 345 -> 386; k444s MAC 1,129,003,567 -> 987,797,005 (-12.51%), H 11181 -> 13500, steps 1827 -> 1861. Mean tail active widths 19.717/20 and 7.729/8. | Misses 25% MAC and <=2% H-increase targets. New lock contracts active space; all-root audits add work. No time acceptance. |
| P24A-2 | Component and margin ablations | not adopted | Rectangle-only Si8 MAC 4,855,582,763, SCF 15 -> 16; locking/two-rectangle MAC 5,212,958,324; doubled margin divisor MAC 4,651,704,502 vs combined 4,459,794,280. k444s respectively 1,029,000,526 / 1,079,833,280 / 1,045,720,012 vs 987,797,005. | Reduced single-call work does not guarantee reduced trajectory work. Twice the margin does not fix syn-s2; keep exploratory 1/sqrt(m), default off. |
| P24A-3 | Generality, final-state dense and thread oracles | adoption gates failed, negative evidence retained | Five Si energies within 1.14e-13 Ha and counts unchanged. Six Si/MgO energies within 8.25e-13 Ha, sampled dense roots pass; MgO60 SCF 19 -> 18 fails unchanged-count gate. Syn-s2 energy +0.004416513597413463 Ha; skew iteration 8 divergence/missing root 12; Si8 iteration 15 buffers 17/19 max error 3.72e-9 Ha versus 1e-9 cap. P15C detects/repairs Si8 at iteration 16. | 16/19 suite cases converge, but convergence alone is insufficient. Pure T=1,2,4,8 identical states; off is baseline-bit-identical. Existing slab/mesh failures are not equivalent final states. MKL execution/dev08 acceptance outstanding. |

### Setup/non-diagonalization profile (2026-11-07, measurement only, Mac observation)
| ID | Experiment | Status | Measured effect | Lesson / limits |
|---|---|---|---|---|
| PROF-1 | Where the 24% outside the five phases goes, si8 real Gamma + standard (30 runs JSON, 30 `sample` runs) | measurement only | wall 0.6148 s, engine setup 0.1269 s (20.6%), five phases 0.4698 s, remainder about 0.011 s engine glue and about 0.008 s process/output; inside setup `ionic_local` 86% (erf alone 42%), `NonlocalOperator::new` 10% | The gap is engine setup, not first-iteration or final-step cost. Setup is serial and unchanged at T=8 (0.127 s of 0.570 s). Cost scales with cell volume (si2 Gamma setup 46% of wall). No code changed. See profile-real-standard-si8.md. |
| PROF-2 | Energy phase breakdown si8 (30 `sample` runs) and k444s symmetry-on | measurement only | si8 0.0535 s: nonlocal 51%, Hartree 21%, XC 12%, local 9%; k444s (10 k) energy 0.098 s = 18% of wall, about 86% k-independent and recomputed per k | Per-k audit repeats common terms; candidates C2 to C4 in the profile note. Samples attach caveat: first runs lost 0.1 s of setup until a suspended-start attach was used. |


## P28A setup memoization (2026-11-07, append-only)

| ID | Experiment | Status | Measured effect | Lesson / limits |
|---|---|---|---|---|
| P28A-C1 | Precompute local short-range integrand; exact squared-G bit-key cache, unchanged public quadrature and G=0 branch | implemented, local identity verified | Mac engine_setup_s medians, 30 runs per variant: gamma 0.032649792 -> 0.005270334 s; si8 0.151553292 -> 0.025463604 s; k444s symmetry auto 0.084302459 -> 0.051255021 s (C1 only). | Extends B18 exact-shell idea with PROF-1 evidence; no interpolation or reassociation. Keep rab after Bessel multiplication. Shared Mac load, no remote acceptance. |
| P28A-C5 | Exact squared-q memo of unchanged projector transforms, preallocated flat cache storage | implemented with C1, bit-identity retained | Incremental C1 -> C1+C5 setup medians: gamma 0.005270334 -> 0.001974292 s; si8 0.025463604 -> 0.008021958 s; k444s 0.051255021 -> 0.028515375 s. All 25 strict case/thread pairs: delta E = 0 Ha, identical density, SCF and work counts; 96 multi-shape/atom/k component cases plus ten mesh cases pass. | Removes original per-G Vec allocation; angular/phase/D arithmetic unchanged. Extra scratch O(N_G * N_projectors). No C5 rejection needed. Exact inputs, binary hashes, raw samples and identity fields in optimization-log.md P28A. No QE source read or remote/QE timing claim. |
## P28B shared audit and energy-only Hartree (2026-11-07; append-only)

Baseline/current main `534595c4d19346dc1fe6cdd5db87415635668f86`; candidate
uncommitted. Mac portable release observations, 30 samples per variant/case,
T=1 standard; si8 real Gamma without symmetry, k444s complex with symmetry.
Full commands, binary hashes, raw timing samples and limitations:
[optimization log, P28B](optimization-log.md#p28b-audit-shared-terms--shared-audit-and-energy-only-hartree-2026-11-07).

| ID | Experiment | Status | Measured effect | Lesson / limits |
| --- | --- | --- | --- | --- |
| P28B-C2 | Evaluate common audit terms once per output density, preserve every per-k expression and weighted summation | implemented, local identity verified; dev08 pending | Energy phase k444s 0.1022862485 -> 0.0153458960 s (-85.00%); si8 0.0553770835 -> 0.0557856665 s (+0.74%, no Gamma saving claimed) | Removes K-fold common work; does not skip audits or repeat A43. Gamma has K=1, so C2 is a multi-k benefit. |
| P28B-C3 | Add energy-only Hartree and share its forward density FFT with the local integral | implemented, local identity verified; dev08 pending | C2-only -> final energy phase si8 0.0557856665 -> 0.0459679575 s (-17.60%); k444s 0.0153458960 -> 0.0117628755 s (-23.35%). Combined baseline -> final -16.99% / -88.50% | Same frequency arithmetic and sums; retain the raw local FFT when Hartree's tiny-negative clamp changes input bits. Existing Hartree evaluator unchanged. |
| P28B-validation | Legacy audit oracle, complete supported CLI matrix, spin/FD, generality | verified locally | 8 audit fixtures x 3 k points plus 4 Hartree grid fixtures: every field identical; 140 CLI combinations: 0 Ha delta, identical density/audits/history/counts; 60 unchanged gamma_real=on multi-k rejections; 3 spin/FD CLI pairs identical; 19 generality outcomes identical (16 converged, 3 pre-existing failures); workspace gates exit 0 | 40 final bench commands exit 0. Explicit optional historical-Si8-anchor bypass supports symmetry/real/standard baselines without changing baseline checks. Local raw outputs under target/p28b-evidence require preservation before target cleanup. No QE or remote acceptance claim. |

## P29A standard-profile stagnation (2026-11-07; Mac observations)

| ID | Experiment | Status | Measured effect | Lesson / limits |
| --- | --- | --- | --- | --- |
| P29A-1 | Ten mean-centered relative-density perturbations <=1e-10 on ten benchmark/material cases | study complete | 100 runs, 80 converged; 32 density-stagnation events / 30 MKL-like trajectory events. Si8 SCF 20/25.5/32 min/median/max; Si16 18/21.5/27. Eight strict-reference pairs agree within 1e-12 Ha; max standard energy error 2.245319e-9 Ha. | Mac reproduces a flat requested-tolerance/noisy-density trajectory without MKL. Elongated/slab also fail strict references; no accuracy claim for them. |
| P29A-2 | Five-interval, factor-two error-progress guard; trial caps 0.5 and 0.1 of window-entry tolerance | experimental, default off; performance selection failed | 0.5 leaves Si8 trajectories unchanged. 0.1 changes Si8 SCF 20/25.5/32→20/24.5/27, median H 1167.5→1179 (+0.99%). All 60 nonstagnating histories identical. | Neither factor passes the frozen work-saving selection rule. Retain 0.1 only as a documented opt-in tail-robustness experiment; no default recommendation. A21's strict schedule is unchanged. |
| P29A-3 | Held-out k222s/slab/Si16/Si-PBE/MgO-PBE, same ten perturbations | negative generality result retained | Si16 SCF/median H unchanged, max H 2318→2389. Slab median H 5420→13080 (+141.33%), still 0/10 converged; training elongated median H 5165.5→11803.5 (+128.51%), still 0/10. All-set H 275581→412309 (+49.61%); successful-only H 172409→172309 (-0.058%). Max guarded energy error 3.164040e-9 Ha. | Guard costs much more than it saves on the whole set. No parameter tuning on held-out results, no timing or QE claim. |
| P29A-4 | Default-path and persistence gates | verified locally | Workspace build/test/fmt/clippy exit 0; five cases × T=1,2,4,8,12 × strict/standard guard-off: 50 final comparisons, every command exit 0, zero energy delta and identical density. Guard checkpoint continuation/option mismatch and legacy wire omission tests pass. | Original strict/certificate thresholds unchanged; standard stays uncertified. Additional checkpoint/session edits are provenance plumbing. |

Protocol, failure treatment, commands and two full example histories:
[standard-stall-study.md](standard-stall-study.md). Raw evidence is in the
architect-approved worktree `P29A-evidence/`, never under `target/`; archive it
before cleanup. This is not an accepted performance optimization or dev08 result.

## P29B standard-profile robust gain (2026-11-07; Mac observations)

| ID | Experiment | Status | Measured effect | Lesson / limits |
| --- | --- | --- | --- | --- |
| P29B-1 | Default strict and standard, same ten P29A seeds on ten benchmark/material cases, plus twenty repeated reference runs | measured | 80/100 converge in each arm; all 120 reference/standard histories match P29A. Si8 paired median total-H strict/standard 1.7849x [95% bootstrap 1.6503, 2.0044], worst pair 1.4036x; strict H median 2059.5, standard 1167.5. | Historical 2.65x H / 3.0x dev08 wall is one trajectory, not robust gain. Strict certificate work included; main-only counts separately retained. No gain/accuracy assigned to elongated/slab failures. |
| P29B-2 | Residual-squared cap on ordinary Ritz ratchet, coefficient 1, four training cases x five seeds | research switch, default off | Training H 8784→7359 (-16.22%), all 20 converge, max energy error 1.213515e-9 Ha. Gamma/Si8 improve, skew/low regress. | One inexact-solve-inspired mechanism; not the rejected P29A delayed stall guard or A21 strict-schedule loosening. Conditional coefficient 0.1 not tested after training criterion passed; freeze before holdout. |
| P29B-3 | Frozen density cap, held-out seeds 6–10; six whole-case holdouts and four training-case new-seed holdouts | general performance selection failed | All-set H 139332→146296 (+5.00%); successful-only 85631→86270 (+0.75%); whole-case holdouts +6.31%, new seeds of training cases -15.16%. Si8/Si16 H -23.08%/-13.55%; material Si/MgO +3.69%/+5.80%. All 40 prior successes remain; ten stress attempts still fail. Max error 8.633592e-10 Ha. | No rule beats standard on the held-out aggregate criterion. Improvements on large Si cells do not generalize; do not recommend as default. No parameter retuning on holdout. |
| P29B-4 | Default-off identity and switch provenance/continuation | locally verified | Workspace build/test/fmt/clippy exit 0; 50 five-case/thread/profile bench comparisons exit 0, energy delta 0 Ha, identical density/SCF counts; exact continuation/mismatch and legacy wire tests pass. | Strict/certificate arithmetic untouched; standard uncertified. session.rs dispatcher and checkpoint/result/inspection edits are additive plumbing; no extra target directory or artifacts edits. |

Full distributions, paired bootstrap, main-versus-certificate counts, mechanism,
commands and limitations: [standard-robust-gain.md](standard-robust-gain.md).
Raw evidence: `P29B-evidence/` (orchestrator archives before cleanup). These are
local pure-Rust work observations, not a robust dev08 wall or QE speed claim.

## P30A deterministic real Gamma pair scheduling (2026-11-07; Mac observations)

Derived from `real-gamma-design.md` sections 6/12 (R7); no QE source read.
Full protocol, work counts, raw-evidence locations and remaining-site profile:
[real-gamma-r7-results.md](real-gamma-r7-results.md).

| ID | Experiment | Status | Measured effect | Lesson / limits |
| --- | --- | --- | --- | --- |
| P30A-1 | One fixed pair per persistent-pool piece, private preallocated FFT/pack/unpack scratch; same pair partition for engine ORC boundary conversions | implemented; bit identity verified | Old serial pair output equals pooled output bit for bit on four shapes, both FFT backends, odd/even/empty/reused blocks and T=1,2,4,8,12. R5 four synthetic plus two Si fixtures retain pre-change energy/density/work bits at all five budgets; option-off 25 case/budget comparisons have 0 Ha difference and identical density/work. Workspace build/test/fmt/clippy exit 0. | No arithmetic, threshold, default, density grouping or complex-path change. Paired density is already pooled but is not used by the engine. No pool accessor added. Scratch scales with capacity, not worker count. |
| P30A-2 | Same-HEAD real-on T=8 before/after; five retained runs per row | Mac observation; dev08 pending | Si8 strict 1.382801 -> 1.108805 s, standard 0.446245 -> 0.371749 s; Si16 strict 8.580309 -> 5.581573 s, standard 2.128872 -> 1.793700 s. Energy/density/main and certificate counts identical before/after. | Sequential before/after batches on shared host, not quiet-node acceptance. Specified Si16 is 1x1x2, not the relayed dev08 2x1x1 orientation. |
| P30A-3 | On/off strict/standard, Si8/Si16, T=1,2,4,8, five samples after warmup | all four overall targets missed | Off vs on T1/T8: Si8 strict 2.192x/1.422x, standard 2.178x/1.424x; Si16 strict 2.417x/1.459x, standard 2.408x/1.441x. On/off T8 wall ratios 1.0001/1.1983/1.4976/1.3665, respectively. | Pair parallelism improves the real path but does not establish parity with complex or equally good scaling. All on/off work trajectories differ; no tolerance/algorithm change to conceal that. Remaining-site timing is reported separately; no guessed serial fractions. Default stays off. |
| P30A-4 | Temporary named-site timing after missed target; three SCFs per case/profile, all sites >=45 calls, source restored afterward | diagnostic only | Si8 strict per SCF: H 0.293172 s of solver 1.036358 s; local pair 0.207341 s, H construction 0.039418 s (pair scratch 0.033034 s), kinetic 0.004962 s, coupling 0.000414 s, output addition/validation 0.002773 s. Full four-condition site table in report. | Largest measured remainder is solver outside H; it includes parallel algebra and serial work, not a measured all-serial fraction. MGS/EVD attribution remains unresolved; no guess substituted. `sample` unavailable (exit 255), named timing fallback used; diagnostic times are not acceptance timings. |
## P32A real Gamma algebra (2026-11-07; append-only, Mac observations)

| ID | Experiment | Status | Measured effect | Lesson / limits |
|---|---|---|---|---|
| P32A-1 | Contiguous finite/Gram/projection/residual/correction loops, Ritz prefix view, projection-buffer swap | candidate, uncommitted | T=1 slice-only medians: Si16 strict 9.552 -> 6.811 s, standard 3.203 -> 2.008 s; old energy/density bits retained | Removes bounds-checked indexing and copies without reassociation; shared Mac load, remote acceptance pending. |
| P32A-2 | Rotation 4x4x256 calls instead of 16x16x16; unchanged Gram shape and outer 256-row ordered tiles | candidate, not accepted | Full T=1 Si8 strict 1.818 -> 1.271 s, steps 320 -> 253; Si16 strict 9.552 -> 5.336 s, steps 329 -> 264. Si8 standard 0.608 -> 0.503 s but slice-only 0.449 s, SCF 15 -> 17. | 4096-MAC calls stay allocation-free; kernel gain is mixed with changed trajectories. All 25 complex-off comparisons identical; kernel max old/naive relative errors 1.01e-15/7.94e-16; dense eigen max 2.86e-14 Ha. R6 slab absolute cap still fails, primitive extra SCF-count condition still fails. |
| P32A-3 | Full-panel column SIMD, straight slice products, alternative bounded shapes, general large GEMM | rejected alternatives | 256x120 times 120x20: old 44.07 us; column SIMD 65.60 us, slices 66.94 us, chosen bounded shape 40.55 us, general GEMM 19.91 us | Full GEMM is faster but allocates packing memory with no public caller-owned scratch API. Small-call alternatives do not consistently improve Gram. New evidence is D0 real-path algebra dominance, not a repeat of complex A49 tuning. |
| P32A-4 | Before/after sampling in task sandbox | blocked | Both actual CLI sampling attachments exit 255, SCF exit 0 | No <30% algebra-share claim. dev08 sampling and performance acceptance requested through orchestrator. |

Full matrix, all candidates, commands, failures and raw-evidence locations:
[real-gamma-r8-algebra-results.md](real-gamma-r8-algebra-results.md).
## P31A — Davidson capacity (2026-11-07, Mac observation; default unchanged)

Reopens A08/A13 only on the new 20–68-band premise in D1a/B17, derived from the
prose note `docs/research/qe-study/qe-solvers-and-subspace.md`; no QE source read.
Adds research multiplier 2..=12, default 6, to all solve paths and restart metadata.
The absolute subspace floor remains 48. No solver/tolerance/default retuning.

The uninstrumented T=1 sweep has 152 runs including paired tight references;
137 converge and 15 fail, all failures at multiplier 2 in strict Si8/Si16.
Si8 complex and real each converge 0/5 versus 5/5 at 6; Si16 complex converges
1/5, real 4/5 versus 5/5 at 6. All failures report the unchanged 400-step Davidson
limit. Median-of-condition wall ratios for 2/3/4 versus 6 are
1.17475/1.05999/1.00160 (failed attempts are censored costs, not speedups).
The last ratio is unresolved noise, not evidence of a 0.16% regression.
Standard median delivered errors breach the prespecified 1.5x rule in three
conditions at 2 and one at 3. Converged strict error is at most 1.16245e-11 Ha.

Si32 strict 2/3/6: wall 143.9721/68.7006/64.9406 s, SCF counts 15/13/13;
shrinking to 3 increases H work despite the same SCF count. **No smaller multiplier
qualifies as a robust default; keep 6.** Small cases with the same effective
capacity expose timing noise and must not be counted as algorithmic speedups.
See [the study](davidson-capacity-study.md) for every row, perturbation distributions,
restart/work accounting, named-site algebra attribution, actual verification logs,
and scope limits. Raw data are in `P31A-evidence/`, not target or tests/reference;
no new QE or dev08 performance claim. Future work needs a different algorithmic
premise or independent dev08 evidence, not another unqualified capacity reduction.

P31A mechanism check (same first perturbation, instrumented T=1 Si16): named
algebra share falls 60.81% -> 54.07% (strict) and 44.42% -> 39.54% (standard) from
6 to 2. Strict EVD steps rise 224 -> 2905, and absolute named algebra time rises
5.6340 -> 20.1646 s. Coverage includes the incremental reduced-matrix extension;
224/2905 EVD samples and 410/5026 extension calls exceed the 30-sample requirement.
This supports repeated-work cost overwhelming cheaper per-step algebra, not a
recommendation to reduce capacity. All 70 default comparisons are bit-identical;
all 70 accounting histories agree; final workspace gates pass.

## P33A — real Davidson workspace and column scheduling (2026-11-07)

Starting main `73f59e347e5b6eb74d8cfdbbe8ce87d55f7132c2`; candidate uncommitted.
Motivation is the relayed post-P32A profile (48.3% algebra, 12.6% memset/memcpy),
not a repeat of A49/P32A arithmetic-changing call shapes or A52 factorization.

| ID | Experiment | Status | Measured effect | Lesson / limits |
|---|---|---|---|---|
| P33A-1 | Reuse initialized V/HV, reduced matrix, Ritz/correction/residual/projection and EVD buffers; preserve compact EVD layout | implemented, whole-task acceptance blocked | 30/30 before/after R5 snapshots byte-identical across T=1,2,4,8,12; all 25 off comparisons retain 0 Ha delta, density and work counts; 8 original fixed-operator tests pass; workspace build/test/fmt/clippy exit 0 | Eliminates solver-owned step allocations; does not establish an allocation-free EVD backend. |
| P33A-2 | Fixed four-column residual/correction/normalization/subtraction and ordered Gram assembly | local observations, not accepted performance | Si8 strict T8 0.762807 -> 0.746050 s workspace-only -> 0.751955 s final. Si16 strict T8 2.710891 -> 2.650844 s workspace-only -> 2.614851 s final, same 14 SCF / 264 Davidson / 3254 H vectors | Five measured runs after warmup per row; effects are small/mixed and do not prove a general gain on a shared Mac. No dev08 utilization claim. |
| P33A-3 | Counting allocator enabled after initial H through solver return; independent EVD/GEMM isolation | BLOCKED: gate 5 fails, test exit 101 | T1 cold solves 2969/20: 91 allocations; 11877/68: 1499. Preallocated EVD sizes 20/40/88/160: 0/32/88/121; direct real AB and A^T B: 0 each | Caller-owned EVD scratch does not cover internal general-GEMM packing/bookkeeping. A library scratch API is needed; do not claim no allocation from a source-only buffer audit or silently change arithmetic. Zero-allocation assertion remains unchanged. |
| P33A-4 | Resolve `matmul_1_4_dyn` interpretation and remaining serial work | investigation complete within stated scope | AVX f64 first dimension counts four-lane row vectors; normal 4x4xK rotation calls select this symbol. Synthetic final T8 measured EVD+MGS shares: 11.79% (2969/20), 8.26% (11877/68) | Not evidence of a missed per-column GEMM loop. Serial shares are lower bounds from wrapper timings, not SCF Amdahl fractions; restart/validation remainder is not separately attributed. |

All 16 Si8/Si16 strict/standard timing rows, workspace-only ablation, routine
counts, exact failed gate, commands and backend follow-up scope:
[real-gamma-r9-workspace-results.md](real-gamma-r9-workspace-results.md).
Raw evidence stays in `P33A-evidence/`; no artifacts or QE source changes.

## P34A — full crystal-group density projection (2026-11-07)

Research option, default unchanged; no performance optimization claimed. New
premise: relayed QE shifted 4x4x4 symmetry-on/off difference survives matched
FFT grids, time-reversal controls and tightened convergence, whereas BAANDI's
mesh-group projection agrees with off. No QE source was accessed.

`[symmetry] density_group="crystal"` keeps the 12-operation k reduction (10
representatives) but projects density with all 48 crystal operations. All must
permute the FFT grid; incompatible grids and spin/derivatives fail closed.
On Si shifted 4x4x4, 32^3, strict T=1, tightened crystal-minus-off is
**-1.5538989465824216e-5 Ha** versus relayed -1.554e-5 Ha: **0.00650279%** error.
Frozen support/partial bounds were 5%/20%; **hypothesis supported**. This is
consistent with different effective k sampling, not a correctness ranking or
proof of the reference implementation's internals. No convention offset fitted.

Independent real-space permutations: unsymmetrized density differs by
2.1268811698e-3 electrons/Bohr^3 under each of the 36 non-mesh operations;
crystal density's all-operation error <=6.2450045136e-17. Unshifted mesh/crystal
energy and density identical. Default path: 100 case/thread/profile/symmetry
combinations retain energy/density bits and work counts. Raw evidence lives in
`P34A-evidence/`; [study and reproduction](../validation/symmetry-full-group.md).
No dev08 timing or default recommendation follows from this experiment.

## P15B dev08 verdict (2026-11-07)

| ID | Attempt | Result | Numbers and lesson | Source |
| --- | --- | --- | --- | --- |
| A58b | dev08 isolation of A58 (band-shared nonlocal, branch `task/P15B-NONLOCAL-BLOCK`, commit 963f94c, base ce64b57), pure/MKL x default/native | rejected; NOT in main | bit-identical (8/8) but slower in every configuration, candidate faster in 0 of 60 paired rounds: si8 T=1 +6.05 %, T=4 +8.46 %, T=8 +7.20 %, T=16 +7.55 %, k444s T=1 +4.74 %, T=8 +2.77 %; sampling: nonlocal bucket 12.4 % to 17.0 %, explains 93 % of the increase. Main never contained this code (branch only), so no revert was needed; ticket P35A closed without a code change. Lesson: fewer projector reads (3.2x) do not pay on this Xeon; judge nonlocal changes on dev08 with paired alternation. | tests/reference/perf/dev08_materials_qe_relayed.md section (6) |

## THR thread-width analysis (2026-11-07; Mac observations, measurement and scratch experiments only)

No repository source changed; every build lives in the session scratchpad. Full tables, commands and limits: [thread-width-analysis.md](thread-width-analysis.md). Mac = Apple M6 with 2 Super, 4 Performance and 6 Efficiency cores, so thread counts above 2 or 6 change core type; none of these is an acceptance number.

| ID | Experiment | Status | Measured effect | Lesson / limits |
|---|---|---|---|---|
| THR-1 | Wall T=1,2,4,8 (5 to 7 runs each, alternating order, spread 0.2 to 2.8 %) and instrumented per-phase wall time, si8 and si16-gamma strict, gamma_real off and on; plus k222s, k222, k444s with symmetry and k444s full at T up to 16 | measurement only | si8 real 1.195 / 0.854 / 0.751 / 0.751 s (1.59x), si8 complex 2.396 to 1.096 (2.19x), si16 real 5.239 to 2.624 (2.00x), si16 complex 8.969 to 3.240 (2.77x); Amdahl fit f = 0.545 / 0.362 / 0.414 / 0.260. si8 real T=8 phase budget: EVD 9.0 %, in-block MGS 6.8 %, nonlocal GEMM 6.0 %, real-Hamiltonian rebuild per SCF iteration 4.9 % (si16 4.4 %), audit + potential + mixing + setup + process 13.2 %, 12-panel GEMM phases 34 % at 1.5x, local FFT pairs 18.6 % at 2.8x. H calls with at most 4 vectors: 177 of 293 | The serial fraction is 0.29 plainly serial plus about 0.26 inefficiency of parallel regions; the width limit is structural (pairs per narrow call, 12 fixed panels, 1M multiply-add guard, serial EVD/MGS/nonlocal). Mac scaling overstates dev08 (relayed: real 1.20x, complex 1.54x on si8) |
| THR-2 | `sample` at 1 ms with suspended start, per thread and per phase, T=8 (si8 and si16, real and complex: 6/6/4/4 runs; si8 T=1: 3 runs) | measurement only | mean busy threads 2.21 (si8 real, 8,042 active samples), 2.99 (si8 complex), 2.69 (si16 real, 22,532), 3.65 (si16 complex); owner thread 82 to 91 % active, other workers 16 to 41 %; EVD 310 to 1,390 samples (100 % owner thread), in-block MGS 563 to 3,263 samples (60 to 98 % owner), real-H build 218 samples | Reproduces dev08 (about two busy). Stolen GEMM tiles carry no caller frames, so GEMM per phase comes from the clock guards |
| THR-3 | Patched `rayon-core` 1.13.0 (scratch copy): spin 4000 rounds instead of 32 before sleeping; plus wake microbenchmark | rejected (negative) | si8 real strict T=4 0.7548 to 0.7766 s, T=8 0.7526 to 0.7913 s; complex T=8 1.1154 to 1.2005 s; si16 real T=8 2.639 to 2.808 s; k444s sym T=8 0.2487 to 0.2654 s. Microbenchmark: 12 pieces after at least 20 µs of serial work cost 35 to 87 µs extra at T=8, 7 to 10 µs at T=4 | Refutes sleep/wake latency as the cause of the narrow width on the Mac; consistent with A47. Not a new idea to retry; dev08 wake data still unknown |
| THR-4 | `scf.rs:804` inner budget = pool size instead of T/outer (scratch edit), k444s and k222s with symmetry, 5 runs | positive, small, not proposed as done | k444s sym wall T=4 0.3083 to 0.2933 s (-4.9 %), T=8 0.2487 to 0.2215 (-10.9 %), T=12 0.1968 to 0.1818 (-7.6 %), T=16 0.1983 to 0.1864 (-6.0 %); k222s sym T=4 -1.8 %, T=8 -1.9 %; density sha256 identical to the T=1 default at T=8 and 16 (one case) | Revives A50 with new evidence: per-worker scratch (inner = pool size) removes the slot collisions that regressed the modulo prototype. Gate set not run |
| THR-5 | k-point scaling: per-k solve time at T=1 (max/mean 1.04, 10 k), sum over SCF of the slowest k solve 0.066 / 0.107 / 0.119 s at T=1 / 8 / 10 (efficiency 0.93 / 0.75 / 0.77), serial remainder 38 to 41 % of the T>=10 wall (k444s sym) and 60 % at T=8 (k222 sym); CLI `symmetry_map` evaluated at seven call sites (71 of 246 main-thread samples) | measurement only | k444s sym T=1 0.7519 s, T=10 0.1965 (3.83x), T=12 0.1968, T=16 0.1983; k222 sym (3 k) T=4 equals T=3 (0.153 s) | Outer width is nk and inner is floor(T/outer) = 1 above T=10; for few k points the serial setup/potential/CLI remainder, not the solver, is the floor |

## QPAR QE stage-by-stage study (2026-11-07; analysis only, no code, derived rows are arithmetic on relayed numbers)

Notes and plan: [qe-stage-by-stage.md](../research/qe-study/qe-stage-by-stage.md), [qe-parity-plan.md](qe-parity-plan.md). No source or measurement was produced; D = derived.

| ID | Experiment | Status | Measured effect | Lesson / limits |
|---|---|---|---|---|
| QPAR-1 | Per-stage and per-vector accounting, QE real Gamma si8/si16 timers (M) against dev08 D0 shares of BAANDI real strict T=1 (M) | analysis | QE: h_psi 49 / 58 %, algebra 11 / 18 %, init 11 / 7 % of PWSCF; about 800 H vectors on si8 from 964 FFT calls (D), 0.22 ms per vector. BAANDI: local FFT 31.9 %, nonlocal 3.0 %, algebra 48.3 %; 0.54 ms per H vector (D, 2.4x); algebra about 1.5 s against QE 0.05 s (D, 30x); complex path 0.61 ms per vector so pairing gained only about 1.2x per band instead of 2x (D, unconfirmed) | The H-vector count is at parity on the unperturbed standard trajectory (772 to 864 against about 800, M); the gap is per-vector FFT cost and algebra. The pair-path FFT figure needs a dev08 microbenchmark before any action (plan Q1) |
| QPAR-2 | Core-count scaling QE np=8 against BAANDI T=8 (relayed batch-7 C/D tables) | analysis | QE 1 to 8 ranks: si16 2.54x, si32 2.97x. BAANDI best-of standard: 1.69x, 1.49x. Gap at one core 3.75 / 4.40, at eight cores 5.6 / 8.8 (D) | Gap grows with cores and size: decomposition over the plane-wave index (QE) against bands/pairs (BAANDI). QE's default layout for k meshes below half the z planes is plane-wave, not pools; its best k configuration (`-nk`) is unmeasured |
| QPAR-3 | Reopen test for B13 (non-orthonormal generalized basis) | premise changed, not yet implemented | Backlog condition "algebra above 40 % of wall on a dev08 profile after native and standard": real strict D0 gives 48.3 % (T=1), 57.1 % (T=8); standard-profile share not measured | Reopens A16/B13 only as an isolated real-path module with the stage-0 rank diagnostic (A52 premise) and the Q1 standard D0; P24A, A13, P31A are not repeated (different structure: unconverged-only work, one rectangle, cheap steps) |
## P36A persistent real state and paired rotation (2026-11-07; Mac observations)

Baseline `4bce6acbe563fb7aec04dd3db876765b05b0110f`; candidate uncommitted on
`task/P36A-REAL-PERSISTENT-STATE`. Full raw samples, ablation, commands and
limits: [p36a-results.md](p36a-results.md), `P36A-evidence/`. No dev08 acceptance
or new QE run. The unchanged zero-allocation requirement remains a blocker.

| ID | Attempt | Status | Measured evidence | Lesson / limits |
|---|---|---|---|---|
| P36A-1 | B21: retain real Hamiltonian/pair plans/nonlocal image and Davidson GEMM/projection/EVD/matrix scratch per k state; refresh potential, invalidate on operator identity or solve-storage settings | numerical identity verified; allocation gate blocked | 70/70 final snapshots match baseline (seven systems, strict/standard, T=1/2/4/8/12); 25/25 complex/off comparisons preserve energy/density/work; full SCF forced-rebuild oracle passes; check.sh exit 0. Alternating five-triple Mac Si8 strict T=8: baseline 0.752295 s, workspace-only 0.711260 s; standard T=8: 0.337523 -> 0.300506 s. Reused-workspace allocation probe still counts 70/72 allocations at 2969x20 T=1/8 and 1499/1502 at 11877x68 (exit 101); isolated EVD sizes 20/40/88/160 count 0/32/88/121, real GEMMs 0 | Per-SCF ownership removes repeat numerical-storage construction but cannot provide the backend's opaque EVD packing storage. Do not claim an allocation-free loop or acceptance. Shared-host outliers are retained; dev08 must decide speed |
| P36A-2 | B22 V/HV only: one real Ritz rotation region containing the two existing fixed-panel trees; no arithmetic/tile/guard change | bit-identical; performance provisional and mixed | Si8 strict workspace-only -> paired: T=4 0.719752 -> 0.727408 s, T=8 0.711260 -> 0.699368 s; five alternating measured triples. Pair-kernel bit oracle passes all five budgets, zero/boundary/Si8 shapes and poisoned output. 2969 rows retain 12 panels each, now 24 independent pieces in one region | Ports A45 without repeating A49's finer tiles. Small individual effects are not established by five samples; no standalone rotation speedup or default acceptance is claimed. Audit/nonlocal overlap and P36B budget changes remain outside this ticket |
## P36B k-point hybrid budget and shared symmetry preparation (2026-11-07)

Reopens A50 only on THR-4's new full-pool scratch evidence, not the rejected
modulo-slot prototype. No QE source was read. Candidate is uncommitted;
[full report](p36b-results.md), raw evidence in `P36B-evidence/`.

| ID | Experiment | Status | Measured effect | Lesson / limits |
|---|---|---|---|---|
| P36B-1 | inner=T on the existing pool, outer=min(T,nk), per-worker Hamiltonian scratch; unchanged fixed arithmetic and ordered density reductions | identity verified; performance not accepted | 162 mesh/profile/symmetry/budget conditions per arm, 648 records across baseline/budget/cache/final: all energy/eigenvalue/history bits and density SHA-256 agree. 25 final bench comparisons pass; all nine-budget spin/FD snapshots match baseline, including aggregate Davidson/H counts 52/189 and 26/434 respectively. | No duplicated budget logic in session.rs/collinear.rs. Existing Hamiltonian scratch supports the full-pool rule; no numerical kernel or real-path changes. Memory grows with T scratch slots per k. |
| P36B-2 | Immutable verified symmetry preparation cached by exact geometry, mesh, time reversal and tolerance; passed to ordinary SCF | ordinary SCF implemented; whole-task BLOCKED on ownership | Cache reuse/invalidation and wrong-engine-geometry rejection tests pass; workspace build/test/fmt/clippy exit 0. Crystal density operations reuse the same detection. | Spin/FD wrappers still call legacy preparation: forwarding needs spin_run.rs/finite_run.rs/finite_temperature.rs plus collinear.rs beyond its budget-only allowance; scope approval pending. No artifact edits. |
| P36B-3 | Alternating before/budget/final, three symmetric meshes x six budgets, warmup plus five retained samples each | noisy Mac observations; no general gain established | k444s T8 0.313165 -> 0.290331 -> 0.270459 s; T16 0.242049 -> 0.247617 -> 0.234434 s. Final T1 regresses k222 +7.72%, k444s +2.26%; final T16 changes k222/k222s/k444s -1.63/-1.81/-3.15%. All 324 timing outputs retain scientific bits. | One-minute load 12.27–18.26, worst cell spread 115.05%; heterogeneous Mac and uncontrolled host activity prevent attribution or acceptance. No task-owned build overlap. Full phase/remainder tables retained; five samples cannot establish small effects. dev08 decides. |

### P36B completion after architect scope extension (2026-11-07)

This supersedes P36B-2's ownership blocker; the earlier timing evidence is retained.

| ID | Experiment | Status | Measured effect | Lesson / limits |
|---|---|---|---|---|
| P36B-4 | Forward the shared preparation through fixed-spin and both FD dispatches; legacy collinear preparation shared by both channels | implemented; all local correctness gates pass | Additional production files: engine collinear.rs/finite_temperature.rs and CLI spin_run.rs/finite_run.rs. 144 original-main/candidate CLI runs (four policies, off/auto, nine budgets) retain scientific JSON and density bytes. Final 162-condition matrix, 25 bench comparisons and 18 accounting snapshots all match baseline; workspace build/test/fmt/clippy exit 0. | Additive methods, existing public signatures/checkpoint formats unchanged. Antiunitary collinear preparation remains rejected. No real-path/artifact changes. |
| P36B-5 | Repeat all 18 three-arm timing conditions after checking pre-batch one-minute load <8; retain five samples plus warmup per arm | protocol complete; no general performance acceptance | Admission load 6.89–7.98 vs original run load 12.27–18.26. All 324 repeated outputs preserve scientific bits. k444s T8 median before/budget/complete 0.313631/0.278012/0.265787 s, ranges [0.294928,0.320271]/[0.265403,0.278791]/[0.262681,0.272319] s, spreads 8.08/4.82/3.63%. | Original batch has 9/18 rows >10%; lower-load repeat still has 12/18. Every such row is marked UNUSABLE, all samples retained. Load admission is insufficient to guarantee isolation; do not claim gains from noisy medians. Full wall/phase min/max tables in p36b-results.md; dev08 decides. |

All task-owned executable copies outside target are removed after verification;
raw evidence, hashes and cleanup manifest remain. No extra target directory.

## P37A — generalized Gram stage 0 (2026-11-07; measurement only)

Derived from the own-words `qe-stage-by-stage.md` and
`qe-solvers-and-subspace.md`; no QE source read. Reopens B13/A16 only on Q2's
relayed real-algebra premise; does not repeat A13 capacity-only tuning, A52
CholeskyQR2, or P24A hard locking. Baseline
`6c34eed6125c3e05338ddef486edcf1a8bb30c9b`, Mac pure Rust, T=1.

| ID | Experiment | Status | Measured result | Lesson / limits |
|---|---|---|---|---|
| P37A-1 | Nonorthogonal unit-norm corrections; capacities 2/3/6; Cholesky and canonical fallback diagnostic on exact original SCF operators/starts | recommendation rejected by prespecified safety rule | 28 histories, 2,410 original solves, 67,257 diagnostic steps. Strict Si8 fallback at relative 1e-8: 64.60/64.57/34.83%; minimum relative Gram eigenvalues -2.429e-16/-1.175e-15/-1.388e-16. Strict Si16 2/3: fallback 41.39/36.98%, minima 3.537e-15/4.074e-15 | Fails <5% fallback and >=1e-12 strict minimum. Standard-only safety cannot override an all-profile rule. Diagnostic uses an independently rebuilt Gram; a cached-only production solver has no safety guarantee here. |
| P37A-2 | A52 twenty-atom skew geometry, uniform and perturbed starts; current real strict path | no general module recommended | 2/3/6 fallback 0/1.08/2.03%; minimum relative eigenvalue 3.727e-6/5.768e-11/4.760e-10. Counted large-product/MGS MAC per step savings 20.97/44.99/18.99%; total MAC ratios 5.70/3.59/0.72. Diagnostic exhaustion 12/8/0 of 12 solves | A52 safety passes this assay, but cheap steps do not ensure convergence or lower total work. Does not explain the historical complex QR2 iteration regression. Six SCF steps remain a bounded, unconverged fixture. |
| P37A-3 | Actual-call real solver GEMM/MGS accounting, generalized rectangle/packed-rotation/contraction dimensions; separate small-solve model | measured subset; full EVD count not certified | Strict capacity-6 MAC/step savings Si8/Si16 18.80/18.62%; standard Si8 45.54/50.28/50.32% at 2/3/6. All si-skew and si-slab rows miss 40%. 1,176 original/control scientific records identical; 25 production benchmark comparisons at T=1/2/4/8/12 give 0 Ha and identical densities | Small EVD is explicitly modeled, not an exact faer instruction count. No wall-time or dev08 claim. Production sources/defaults/artifacts/UPFs unchanged. |

Per-case/profile/capacity spectra, percentiles, Cholesky failures, rank deficits,
work tables, censored outcomes, fixed decision rule, exact commands and gate
results: [generalized-gram-stage0.md](generalized-gram-stage0.md). Raw evidence:
`P37A-evidence/`, including the superseded full-rebuild pilot. No module or
default is adopted; no tolerance is relaxed.

## PFFT pair FFT efficiency (2026-11-07; measurement and scratch experiments only, Mac observations; details [pair-fft-efficiency.md](pair-fft-efficiency.md))

| ID | Experiment | Status | Measured result | Lesson / limits |
|---|---|---|---|---|
| PFFT-1 | Count of 1D transforms, flops and element moves of the pair apply (pack, inverse, V multiply, forward, unpack) at 36^3 and 72x36x36 against the ideal stick-pruned algorithm and the unpruned count | analysis (C) | Pruned count 4.156 / 8.833 Mflop against unpruned 7.236 / 15.406 (57.4 %); ideal order x,y,z with mirrored forward 4.156 / 8.494 (QE order z,y,x 4.586 per transform at 72x36x36, 7 % worse); 4466 / 7344 lanes; 9.2 / 18.4 MB of lane element moves outside the kernels; library working arrays 1.3 / 2.5 MB against a 1 MB dev08 L2; pack plus unpack 5.7 / 11.8 us (1.5 %) | The flop count is already within 4 % of the ideal pruned count; the two unprunable full-grid passes are 58 % of the flops. The half sphere gives no stick reduction in the pair path (pairing is already the r2c-equivalent). The loss is data movement plus kernel rate |
| PFFT-2 | Stage breakdown by kernel-skip ablation on a bit-identical clone of `OrderedFusedPlan`, 40 samples of 100 applies; 1D kernel microbench lengths 25, 36, 48, 72 | measured (O) | 1D kernel calls 66 % (36^3, 274 of 418 us) and 64 % (72x36x36, 510 of 799 us); zero-fill, gather, scatter and scale passes 34 / 36 %; inverse z stage alone 32 % of the total. rustfft 62.0 ns per lane (15.0 Gflop/s) at length 36, 102.7 ns (21.6) at 72, 48.2 (12.1) at 25, 83.6 (16.0) at 48; per-lane calls only 5 % above batched; scale pass 8.2 ns per lane at length 36 | The kernel cost equals lanes x ns per lane (277 / 519 us predicted, 274 / 510 measured); a perfect layout cannot go below the kernel floor 274 / 510 us (library 373 / 704) |
| PFFT-3 | Mirror x-slab pair plan (scratch): inverse x on sticks, per x slab y inverse, z inverse, V multiply, z forward, y forward in one 31 KB buffer, forward x on sticks; unnormalised transforms | candidate: recommended for ticket T1 | Paired, 40 rounds: 373.1 -> 323.5 us (0.867) at 36^3, 704.2 -> 641.8 us (0.911) at 72x36x36; nine shapes 0.66 to 0.90 (largest for a vacuum slab 36x36x72: 0.664); max relative difference to the library 3e-16 to 8e-16; movement share 34 to 36 % -> 17 % | Gain is from not forming the full-grid intermediates and from the contiguous potential multiply and forward passes; it is a rounding change for the real path only (new plan type; `fused.rs` and the complex path untouched). Dev08 gain may be larger because the old working set exceeds the 1 MB L2 |
| PFFT-4 | Variants of PFFT-3 and of the old structure | rejected or neutral (O) | Index-table transpose instead of zero-fill plus scatter 0.896 / 0.943 (worse than PFFT-3); per-lane instead of batched rustfft calls +1 percent; dropping all per-lane scale passes in the old order 0.964 / 0.955 of the clone (not alternated); existing opt-in `PlaneFusedPlan` (A27) 1.024 / 0.983; lane-vectorised SoA slab kernels (four adjacent x slabs per group, plain Stockham 4x3x3 and 4x2x3x3) 0.916 / 0.961 and kernel at length 36 58.3 ns per lane against rustfft 62.0, at 72 147 against 103 | Layout is exhausted at about 0.87 to 0.91 on the Mac; further gain needs a faster kernel (vendor MKL 1D batches or SIMD kernels on dev08). A27 stays neutral, A29 not revived, A11 reconfirmed on the Mac (ML is its lane-vectorised form) |
| PFFT-5 | Arithmetic projection of the QE gap from the relayed per-call times and the measured flop rates | analysis (E) | Per transform equivalent at 36^3 / 72x36x36: current 213.6 / 486.9 us (R), mirror 185 / 444 (E), unpruned MKL adapter 195.0 / 488.5 (R), QE 135 / 341 (R); kernel-only floor of pure rustfft 156 / 373 (E); MKL 1D stick pipeline at QE's flop rate 97 / 235 (E, +25 percent allowance) | Pure Rust layout reaches and beats the unpruned MKL adapter by 5 to 9 percent but stays 1.3 to 1.37x behind QE; exceeding QE needs MKL sticks or faster SIMD kernels, both needing dev08 microbenchmarks first |
| PFFT-6 | Reading of the MKL adapter (`dft_mkl::Plan3d::transform`) | new evidence (code reading, unmeasured) | The adapter copies the input into its buffer and the result out around an in-place descriptor (two full-array copies per transform); the relayed MKL full figure 195.0 us at 36^3 includes them, QE's raw call is 135 us | Part of the 60 us difference between "MKL full" and QE is probably adapter overhead; check with a raw in-place call on dev08 before any MKL pair-path work (backlog B33) |

## P38A band-count study (2026-11-07; Mac observations)

Source `9d32045155007d893bfd9ccf86585167e7c2ab4d`; no production changes.
Protocol, full condition tables and limitations: [band-count-policy.md](band-count-policy.md).
Raw evidence: `P38A-evidence/`, external UPFs only. No dev08 speed acceptance.

| ID | Attempt | Result | Numbers and lesson |
| --- | --- | --- | --- |
| P38A-1 | A19 reopened as a conditional multiplet study, extras 1–12, strict/standard, T=1, three repetitions | simple first-empty-multiplet oracle rejected; default remains four | Si8 strict main H: 2063 at four extras, 1253 at the six-member multiplet completion, 1024 at seven. Held-out Si16 first-empty multiplet has two members, but two extras cost 3205 H versus default 2657 (+20.62%, exceeding the 10% per-case gate), with 14 versus 12 SCF steps. Threshold sensitivity 1e-7 through 1e-4 Ha does not change either multiplicity. Converged-spectrum correlation is not a validated pre-SCF/first-iteration predictor. No auto option implemented. |
| P38A-2 | Pre-SCF free-electron first-empty shell completion on the Si8 training case | rejected diagnostic | Shells end at bands 1/7/19, selecting three extras above 16 occupied; measured H is 2321 versus default 2063 (+12.51%). The free-electron shell differs from the physical first-empty group. No fitted policy or independent held-out success is claimed. |

Completed 2026-11-08 (local date): 24 cases × two profiles × twelve counts ×
three runs = 1728 samples; 1095 converge, 294 reject incompatible basis/band
requests, 306 exhaust SCF budgets, and 33 diverge. Held-out syn-s2 strict extras
1/2 reach energies +4.4165136e-3 Ha above default; completing its first-empty
group does not restore the default fixed point. Five ambiguous spectrum rows
were replayed with identical energy bits/density SHA-256 and checked at fixed
density (maximum all-root residual 9.9977605473e-6 Ha, grouping threshold 1e-4 Ha).
All 543 successful strict scans retain certificates (maximum certified residual
9.9881473239e-10 Ha). All 70 supported default comparisons and cross-thread
checks are bit-identical; 30 forced-real multi-k combinations reject unsupported
inputs. Final workspace gates, release probe and evidence audit exit 0.
Wall medians are loaded-Mac observations (load 2.6733–29.6738, including two
concurrent T=1 collectors), not speed acceptance. Integer four stays unchanged;
Phase B is not entered and no recommended automatic option is introduced.
## P39A mirror-slab pair plan (2026-11-07; Mac observations, acceptance blocked)

Source: [p39a-results.md](p39a-results.md), starting HEAD
`f66e3966eceb5cedaca6411910f6738e3c8f4866`; raw evidence in `P39A-evidence/`.
Derived from T1/PFFT-3 and its own-words FFT study-note context, not QE source.

| ID | Attempt | Status | Evidence | Lesson |
|---|---|---|---|---|
| P39A-1 | Separate pure-Rust `PairFusedPlan`: plane-major sticks, mirror slab, batched rustfft, final 1/N; local and pair-density traits | implemented, uncommitted; not accepted | Ten shapes: maximum local relative error 2.703408120714201e-16, pair-density 4.821305842995627e-16; repeated/thread bits agree at T=1,2,4,8,12; six cold allocation probes count 0; stable numeric-buffer addresses. Complex/off 25 comparisons retain 0 Ha delta and identical density. Workspace build/test/fmt/clippy exit 0. | Backend cfg belongs in numerics (physics has no `mkl` feature). Complex/MKL aliases and shared fused modules stay unchanged. Engine density still uses completed complex orbitals. |
| P39A-2 | Pair microbenchmark, 30 alternating samples at 25^3/36^3/48^3/72x36x36 | Mac observation only | Library to mirror us/band: 70.401→60.326, 223.668→180.353, 504.943→419.628, 389.631→334.583; ratios 0.857/0.806/0.831/0.859. No concurrent task builds; shared-host load, all raw samples retained. | Confirms layout benefit on this host; does not establish dev08 cache effects or QE parity. No SIMD/MKL/threading variant revived. |
| P39A-3 | Strict/standard Si8/Si16 alternating end-to-end runs, T=1/8, five retained pairs; R5 snapshots across five budgets | **BLOCKED: strict Si16 benchmark energy gate fails** | Seven-system snapshot strict maximum delta 7.815970093361102e-13 Ha and all 14 profile groups deterministic. Exact benchmark-suite Si16 strict delta -1.438849039914203e-11 Ha at T=1/8 exceeds 1e-12 Ha. Snapshot Si16 and benchmark Si16 differ in seed/atom ordering. Main work and trajectories change; complete timing/work table in report. | Kernel norm agreement at roundoff does not guarantee finite-tolerance SCF energy agreement at 1e-12 Ha. Do not accept faster trajectories while the accuracy gate fails; no tolerance/seed/default tuning. Existing primitive count and standard slab absolute-cap failures are separately retained. |

### P39A architect disposition and completed validation (2026-11-08)

This supersedes P39A-3's acceptance block, without deleting its measured
inter-trajectory delta. The architect replaced that gate with delivered accuracy
against the same fixture's complex-tight reference: `abs(new-tight) <=
max(2*abs(previous_real-tight), 2e-11 Ha)`, tight energy tolerance 1e-12 Ha,
density residual 1e-10, gamma_real off. Production FFT/solver code and all
measured real inputs are unchanged; component gates are unchanged.

| ID | Follow-up | Status | Evidence | Lesson |
|---|---|---|---|---|
| P39A-4 | Recompute seven exact snapshot references and two exact CLI benchmark references, using the existing target tree; apply revised strict gate and unchanged standard relative criterion | **DONE, local ticket gates pass; dev08 performance pending** | All 11 strict rows (seven snapshot fixtures plus Si8/Si16 T1/T8) and six standard rows pass. Exact benchmark Si16 tight -6.31292444216179334e1 Ha; previous/new errors 1.49071865962469019e-11 / 5.18696197104873134e-13 Ha, cap 2.98143731924938037e-11 Ha. Largest new snapshot error 3.09086090055643581e-12 Ha, cap 2e-11 Ha. Full-precision energies/deltas and hex energies in p39a-results.md / delivered-accuracy.json. Revised summarize.py exit 0; workspace check exit 0; five-pair Mac tables already complete. | A reference-based delivered-accuracy gate distinguishes SCF trajectory spread from accuracy loss. Historical slab absolute-cap exceedance and primitive 10/12 count are pre-existing and explicitly excluded by the architect. Mac time ratios remain changed-work observations, not isolated FFT speedups or remote acceptance. |

## P39B MKL caller-owned in-place FFT (2026-11-08)

| ID | Attempt | Status | Evidence | Lesson |
|---|---|---|---|---|
| P39B-1 | Add checked Plan3d in-place entry, allocation-free Fft3d fallback with caller scratch, and GridFusedPlan dispatch | Local implementation and gates complete; dev08 runtime/performance pending | Baseline 010b4e73; pure Rust 25/25 SCFs at T=1,2,4,8,12: energy delta 0 Ha, energy bits/density hashes/work counts equal. Workspace gates and MKL feature all-target check/clippy exit 0. Per MKL FFT: two copies removed, 32N bytes copied payload (1492992 bytes at 36^3). [p39b-results.md](p39b-results.md) contains commands, hashes and durable comparison table. | New premise from relayed raw FFT rates: wrapper/raw 191.13/147.25 us at 36^3, 486.30/323.29 us at 72x36x36; 25^3 neutral. These motivate the patch but are not patch timings. Does not repeat A34/A35; no algorithm/thread/default change. Vendor alignment-dependent bits still require remote tests. |

## P40A mixing distribution survey (2026-11-08; append-only)

Source HEAD `49877d03a0ef8476cffc2fafd45021c6fe1d3d43`; Mac pure-Rust
observations, measurement-only test and documents. Derived from the own-words
stage 2/5 note `docs/research/qe-study/qe-stage-by-stage.md`, never QE source.
Premise change against P17A-5: P29A/P29B's 20–34-SCF spread makes a single
trajectory an invalid evidence unit. A21's tolerance factors and A22's RMM
algorithm are not retried. Protocol and decision rule were archived before
measuring; no default, tolerance, production source or accessor changes.

| ID | Attempt | Outcome | Measured evidence | Lesson / source |
|---|---|---|---|---|
| P40A-1 | A default; B beta 0.5; C beta 0.7/history 8; D beta 0.5/history 8; E existing Hartree metric; F C+Hartree selected by frozen B–D distribution ranking | All alternatives case-dependent; no candidate universal default | 2,244 SCFs, 10 cases, complex plus real on seven Gamma cases, ten paired starts per profile plus two tight references per arm/path. 1,716 converge; 528 failures remain confined to elongated/slab, no new failures. C standard improves 8/10 complex and 5/7 real case medians, but strict only 6/10 and 2/7, with >10% regressions. | Uniform-default rule requires >=70% improved cases, no >10% regression, no new failures and <=1.5x delivered error in both profiles/paths. Every alternative fails available-case criteria as well as lacking stress-case references. `mixing-survey.md`; `P40A-evidence/analysis.json`. |
| P40A-2 | F distribution and paired hierarchical bootstrap (20,000 case/seed resamples) | Standard work gain is conditional; strict regressions and accuracy limits prevent adoption | Standard paired SCF complex -28.205% [-35.282,-22.222], real -22.222% [-40.000,0.000]; total H including failures complex -20.166% [-26.786,-10.673], real -12.429% [-21.455,-0.803]. Si8 standard medians complex/real 25.5/27.5→16/16, strict 14/14→16/16 (+14.29%). E alone Si8 complex standard max error 2.245319e-9→1.666601e-8 Ha. | Do not infer strict improvement or accuracy from standard SCF savings. F is adaptively selected, intervals exploratory, no wall-time or dev08 acceptance. Tight cross-arm max delta 7.247536e-13 Ha; strict default error can exceed 1e-12 Ha despite certification. All workspace gates exit 0; 25/25 same-source benchmark repeat pairs preserve energy, density and work. |

## P41A trajectory divergence (2026-11-08; measurement only)

Premise: P29A/P29B and P40A show a trajectory distribution, not a universal
mixer/default improvement. This isolates the loose-state density channel;
it does not repeat A21's rejected fixed-factor loosening or propose a tuning.
Derived hypothesis inventory uses the own-words QE stage note, sections 2/4/5;
no QE source was read. Local Apple Silicon/pure-Rust evidence only.

| ID | Measurement | Result | Numbers and lesson |
|---|---|---|---|
| P41A-1 | First-map loose/tight pairs, ten starts × Si8/Si16 × complex/real × strict/standard (80 pairs) | Initial occupied-subspace sensitivity directly observed before mixing | Input density separation <=1.246e-12 Bohr^-3/2 produces loose output separation up to 5.124e-3; tight output <=5.919e-13. Loose occupied-subspace defect reaches 1.393 (dimensionless). Strict first-map occupied/empty gaps are roundoff-scale, versus 0.11072348615 Ha in tight solves; standard loose gaps about 0.00163–0.00170 Ha. Nearly equal loose eigenvalues/residual norms do not imply the same occupied density. This is more than an occupied-space gauge rotation. The particular EVD/expansion operation is not isolated. |
| P41A-2 | Paired default / every ordinary solve tight at 1e-10 Ha / same-potential tight density only; unchanged mixing and outer profile predicates | Loose density is the dominant observed large-Gamma standard SCF-count variance channel | Si8 complex SCF 20/25.5/32 -> 9/9/9 -> 14/14/15; real 21/27.5/30 -> 9/9/9 -> 13/14/15. Si16 complex 18/21.5/27 -> 7/7/7 -> 13/13/13; real 15/19/25 -> 7/7/7 -> 7/7/7. Density-only removes 98.4%/97.4% of Si8 count variance, all Si16 variance. Work rises: Si8 complex median total H 1167.5 -> 2834.5 / 5854.5; Si16 complex 1854.5 -> 4654 / 9444. These are expensive oracles, not optimizations. |
| P41A-3 | Strict profile and generality controls | No universal iteration-variance cure or default recommendation | Si8 strict is already 14 SCFs for every seed on both paths; density-only takes 17–19 complex and 18–19 real. k222 standard is 10/10/10 SCFs and 1600 H for every seed; PseudoDojo Si standard is 11/11/11 and 6974 H, strict 15/15/15 and 16526 H. Stable multi-k/material controls prevent generalizing large-Gamma sensitivity to every input. Tight strict final certificates use 1e-10 Ha too; three earlier logs with ordinary final certificates were superseded and retained. |

Full protocol, seed-pair amplification, work accounting, threshold histograms,
reference comparison limits and source-only trace method:
[trajectory-divergence.md](trajectory-divergence.md). Raw logs, analysis and
instrumentation diff: `P41A-evidence/` (archive before worktree cleanup).
Only a doc-hidden visibility accessor changes production source; numerical
bodies/defaults and artifacts/* are unchanged. No performance acceptance or
new QE measurement is claimed.

P41A stopping attribution (default 1e-10 ensembles): all 360 Gamma seed-pair
first step-count divergences have actual passing/failing band comparisons.
Standard crossings: 254 occupied, 4 first-empty, 26 higher-buffer bands;
strict: 102 occupied, 90 first-empty, 97 higher-buffer (multiple bands can
cross per pair). The first empty is mandatory but uses the empty threshold in
standard; nominal required-prefix metadata must be converted before comparing.
Si8 complex standard seeds 1/2 branch at SCF 6, 2 versus 5 Davidson steps,
band 15 passing 1.65444e-4 < 1.67869e-4 Ha versus failing
1.73930e-4 > 1.51819e-4 Ha. Their density-residual norms already differ by
2.51407e-5 Bohr^-3/2 at SCF 1, before internal stopping flags differ.
Discrete stopping is a later amplifier, not established as the initiating event.
Noise amplitudes 1e-11/1e-10/1e-9 retain broad Gamma distributions, but medians
change: this is not proof of exact amplitude independence. No default is tuned.

P41A verification complete: 60 groups / 600 converged reported trajectories,
2700 seed pairs, 120 instrumented default histories identical field-for-field,
80 first-map pairs, and zero unresolved default branching pairs. Workspace
build/test/fmt/clippy exit 0; 25 five-case/thread benchmark comparisons under
`--require-identical` have energy delta 0 Ha and identical density. README
freshness: nine tests, exit 0. `P41A-evidence/verify.py` exits 0 and confirms
that the production function body and artifacts/* are unchanged. Three
superseded certificate-control groups and one duplicate parallel control are
kept/accounted for without inflating the ten-seed sample size. Work remains
uncommitted; archive raw evidence before cleanup.

## P42A truncation-boundary study (2026-11-08; Mac count observations)

A19/P38A reopened on the P41A perturbation-distribution premise, specifically
last-requested-band closure rather than first-empty multiplet completion.
No production/default/solver/tolerance change; no QE source read.

| ID | Attempt | Decision | Evidence and lesson |
|---|---|---|---|
| P42A-1 | Ten P29A starts per extra count, two strict references, complex/real; own-density N+1 and wider spectral probes | Converged boundary closure is neither necessary nor sufficient for stable SCF counts; no auto option justified | Si8 complex extra=4 reproduces 20/25.5/32 SCFs, median H 1167.5. Extra=8 gives 10/10/10 and H 629 despite a wide-probe gap 8.48e-12 Ha inside a triplet. Extra=6 has gap 0.0719762645 Ha and complex 10/10/11, but real 24/25/30. Same final spectral gap does not determine both arithmetic-path trajectories. |
| P42A-2 | Frozen descriptive oracle: first resolved gap >=4 extras, minimum over k relative gap >=0.05; fallback four, no probe cost charged | Fails per-case work bound even as an optimistic converged-spectrum diagnostic | Training cubic synthetic: H 98 ->112 (+14.29%). Held-out k222s: H 731 ->1083 (+48.15%), SCF 9 ->11; strict reference H +17.40%/+23.37%. Threshold not retuned. No first-loose predictor or recommended option follows. |
| P42A-3 | Audit spectral completeness and equal-trajectory claims | Residual certification alone does not certify lowest-root completeness; equal counts are weaker than equal scalar histories | Narrow cold probes can omit roots despite <=1e-5 Ha residuals. Wide own-density probes retain disagreement/unknown flags. Si8 extra=8 has one integer-work sequence but ten distinct exact scalar histories. Do not rank unresolved roundoff splittings as physical gap ordering without the resolution sensitivity check. |

Protocol, full coverage status, correlation/accuracy/failure tables and commands:
[band-boundary-policy.md](band-boundary-policy.md), raw `P42A-evidence/`.
All numbers above are local deterministic count observations, not dev08 wall
acceptance. Phase B and ten-start ordinary strict candidate comparisons remain
conditional on a justified rule; none was selected. The integer default stays four.

P42A final coverage: 25 cases, 44 applicable case/path groups, 676 conditions,
8,112 attempted starts. Standard: 4044 converged / 2716 unsuccessful; strict
references: 810 certified / 542 unsuccessful. Of 404 converged reference pairs,
396 meet the fixed 1e-12 Ha energy plateau. Wide gaps support 335 fully converged,
reference-qualified condition correlations; missing/incomplete spectra remain
unknown. Held-out diagnostic median H ratio = 1.000 and median spread change = 0
on both paths; MgO H +37.87%, error ratios 5.9995/2.5067/2.8943, large synthetic
real loses previously convergent seed 2. All recommendation conditions are not
satisfied. Full evidence audit exit 0; workspace build/test/fmt/clippy exit 0;
70 exact default repeat pairs plus 14 cross-thread groups pass, 30 unsupported
forced-real multi-k inputs reject. The unchanged P29A harness exactly reproduces
all ROW/END records for Si16 real extras 8 and 10. No production option added.

## P43A fixed occupied-band fraction (2026-11-08; stopped negative study)

A19/P38A reopened on P41A loose-density sensitivity and the P42A ensemble
premise, with a frozen f=0.375 before measurement. No solver/default change.

| ID | Attempt | Decision | Evidence and lesson |
|---|---|---|---|
| P43A-1 | New skewed 16-atom Si Gamma cell, 32 occupied; extra=max(4,ceil(3*occupied/8))=12 versus 4; new P29A density seeds 101–110, real and complex, two strict references per count/path | Reject; stop before further systems and Phase B | Complex standard SCF 18/18/18 ->17/17/17, median total H 1612 ->1815 (+12.593%, beyond frozen +10% gate). Real SCF 19/19/19 ->17/17/17, H 1753 ->1814.5 (+3.508%). All 40 standard starts converge, all eight strict references certify; all accuracy ratios <=1.5. Increasing buffers reduces steps but increases work on a new already-stable cell. No retuning or path-specific exception. |

Full protocol, units, accuracy, strict work, commands and scope:
[band-fraction-rule.md](band-fraction-rule.md). Raw evidence in
`P43A-evidence/`. Si32, elongated Si12, MgO and Si8/Si16 consistency reruns
were not run after the mandatory stop; they are not passes. No existing
8-atom synthetic benchmark fixture was found. These are Mac count observations,
not dev08 wall acceptance. Workspace build/test/fmt/clippy exit 0; 70 exact
default repeat pairs, 14 cross-thread groups and 30 expected forced-real
multi-k rejections; evidence audit exit 0. Integer four remains unchanged,
and no auto option is implemented or recommended. No QE source was read.

## P45A atomic starting subspace (2026-11-08; stopped negative prototype)

A07/B27 reopened on P41A's first-map sensitivity premise. Derived from stage 1
of the own-words `qe-stage-by-stage.md`; no QE source read. Test-only UPF atomic
orbitals, projector-rule radial transforms, real harmonics, rank-checked MGS
and first-H Rayleigh–Ritz; no production/parser/default/tolerance changes.

| ID | Attempt | Decision | Evidence and lesson |
|---|---|---|---|
| P45A-1 | Paired seeds 201–210, standard plus two tight references, Si8/Si16 complex and real | Strong cubic-cell improvement, insufficient for adoption | Si8 complex SCF 20/25/29 ->13/13/13, real 19/25.5/31 ->13/13/13; Si16 complex 17/21/26 ->9/9/9, real 15/20/21 ->9/9/9. H medians including initial H block: 1138.5/1197.5 ->658; 1863.5/1721.5 ->948. First-map pair density separation falls from 1.59e-3–3.86e-3 to <=4.13e-14 Bohr^-3/2. |
| P45A-2 | Same frozen rule on externally archived skewed Si16, complex | Reject; mandatory stop before further controls or Phase B | Free-electron SCF 18/18/18 ->atomic 11/12/12: spread 0 ->1 violates the required >=50% reduction. Maximum energy error 9.283383e-10 ->2.094893e-9 Ha (2.2566x); density error 1.281516e-7 ->2.010910e-7 Bohr^-3/2 (1.5692x), both violate <=1.5x. First-map separation 1.147978e-13 ->2.598881e-5 Bohr^-3/2. Atomic starts remove the cubic pathology but introduce sensitivity in an already stable skewed cell. |

Primary evidence: 120 completed SCFs (100 standard, 20 tight strict), all
converged; strict reference energy plateaus <=8.243e-13 Ha. Initialization
charges 32/64 H applications and separately records radial evaluations and
Gram/rotation MACs. Quoted H totals are lower bounds on full setup work:
MGS/EVD and radial overhead are not converted into H equivalents. No positive
full-cost gate or performance acceptance is claimed; independent spread and
accuracy failures already require stopping. One already-started real-skew
control was interrupted; its one completed strict reference is excluded from
the primary matrix. PseudoDojo Si/MgO and k222 controls were not run after the
stop. Full protocol, arithmetic, failure logs and commands:
[atomic-start-prototype.md](atomic-start-prototype.md), `P45A-evidence/`.

## P46A atomic-start option (2026-11-08; R2 evaluation pending)

A07/P45A reopened by explicit architect authorization for a default-off option
and a new preregistered rule on new systems. This does not overturn the P45A
negative verdict or select a default from silicon. Derived from P45A and its
own-words stage-1 note; no QE source read.

| ID | Attempt | Current result | Lesson |
|---|---|---|---|
| P46A-1 | Production PP_PSWFC data, first-density atomic Ritz start, both existing Davidson initial-vector paths | P45A consistency: ten standard seeds 201–210 each give Si8 13/13/13 and Si16 9/9/9 on complex and real Gamma | The test-only span ports without changing density, mixing or tolerances. Initial H count 32/64; Gram/rotation MACs 4,940,416/37,856,000; MGS MACs 5,890,496/47,698,560. These are separate work units, not a total-work acceptance. |
| P46A-2 | Synthetic deficient span, plane-wave completion, missing PSWFC, spin, FD and restart contracts | Completion exposed unreserved append capacity; fixed by reserving requested bands before appending. Gamma completion uses conjugate cosine/sine pairs. | P45A full-rank silicon did not qualify completion. Rank deficiency rejects; supplied warm vectors bypass construction. Defaults remain numerically unchanged. |

[Protocol and implementation](atomic-start-option.md); raw `P46A-evidence/`.
New-system R2 evaluation has not begun: the ticket's si32 sample-count conflict
and dimensional conversion of MACs into the total H-work bound require a
preregistered resolution. No recommended-option or dev08 speed claim is made.

P46A local gates: workspace build/test/fmt/clippy exit 0; 70 pre-hook/option-build
default identity pairs (0 Ha difference, identical density SHA-256), 30 expected
forced-real multi-k rejections, and 48 atomic whole-SCF cross-thread density
fingerprint matches. All 48 P45A consistency SCFs converge; reference plateaus
<=3.980e-13 Ha. Evidence verifier exits 0. These are implementation/consistency
checks, not new-system R2 acceptance. Changes remain uncommitted.

### P46A R2 clarification and completed stop (2026-11-08)

Architect clarification before data: si32 uses five standard + one tight strict
per arm; other cases ten + two. Cost is median single-thread whole-SCF wall,
paired alternating order, same Mac/build; H/MAC diagnostics remain separate.

| ID | System and observed result | R2 outcome |
|---|---|---|
| P46A-3 | si32 real: SCF 19/29/47 -> 11/11/11, median wall 14.786158 -> 6.205336042 s; maximum atomic standard energy error 1.681144112808397e-11 Ha | No failed evaluable condition; both strict references certify, single-reference strict-error ratio unavailable |
| P46A-4 | k444s: SCF 10 -> 9, median wall 1.233611104 -> 1.255154646 s (+1.7464%); maximum atomic energy error 2.2710011648996442e-10 Ha, 1.01965 times default | Observed per-system gates pass; timing increase is a local observation |
| P46A-5 | k222: SCF 10 -> 9, median wall 0.3089391875 -> 0.2814508955 s; maximum energy error 5.6632920575339085e-11 -> 8.8816420884541e-10 Ha (15.68283 times) | Reject recommendation under frozen factor-3 accuracy bound; mandatory stop |

60 native SCFs completed (50 standard, 10 certified strict references). No
convergence loss. No later MgO/Si-material/low-symmetry data collected, and no
retuning. Whole-SCF wall includes startup: fewer H applications do not remove
radial/MGS/EVD costs. Non-isolated Mac observations only; no dev08 acceptance.
Full records and automatic stop marker: `P46A-evidence/r2/`; protocol and
limitations: [atomic-start-option.md](atomic-start-option.md). Default unchanged;
the architect-authorized opt-in code remains uncommitted for review.

## P47A — switchable kinetic preconditioners (2026-11-08)

| ID | Mechanism | Status | Measured result | Lesson |
|---|---|---|---|---|
| P47A-1 | Per-band mean-kinetic TPA damping; positive-clipped kinetic Ritz denominator, floor 0.1 Ha; default unchanged | Screening stopped, no recommendation | Si8 complex strict, including certificate: default/TPA/ritz_shift steps 318/268/266, H vectors 2022/1944/1843, SCF 14/14/14. All certify. Tight-reference energy differences 0/4.973799150320701e-14/4.618527782440651e-14 Ha. | The frozen 1.5-times-default strict error test cannot pass for either variant when the computed baseline difference is zero. This is a reference-resolution limitation, not evidence of physical accuracy degradation. Stop before the standard ensembles or broader performance suite; do not add a post-hoc error floor. |

A14's potential/nonlocal diagonal augmentation is not repeated: only the damping
or denominator sign clipping changes, retaining kinetic h_diag. The default
already subtracts the Ritz value. Constants and rule R were documented before
measurement in [preconditioner-variants.md](preconditioner-variants.md).
Counts are local reproducible observations, not dev08 acceptance or QE results.
Raw commands and observations: task worktree `P47A-evidence/`.

### P47A architect amendment and resumed result (R1)

The architect explicitly amended only strict accuracy to
max(1.5 * default error, 1e-12 Ha), after only Si8 complex strict selection
evidence had been seen. P47A-1's zero-error ratio stop is superseded, not a
physics failure. Identity/determinism outputs existed separately; no standard
five-start performance ensemble had been evaluated for selection.

| ID | Mechanism | Status | Measured result | Lesson |
|---|---|---|---|---|
| P47A-2 | Same fixed TPA and positive-clipped Ritz options, unchanged standard accuracy and SCF gates | Both rejected for recommendation under R1; default unchanged | Independent tight run: 20 SCF, -31.3218036010317320 Ha. Both strict variants pass R1. Si8 standard TPA seed 402: SCF 21 -> 28, steps 95 -> 87, H 1011 -> 1276, errors 1.454594666938647e-9 -> 3.293294525974488e-10 Ha. Ritz_shift seed 401: SCF 28 -> 29, steps 100 -> 90, H 1253 -> 1283, errors 8.097416070995678e-10 -> 1.508126956650813e-10 Ha. | All starts converge and accuracy passes; SCF counts increase. Reduced Davidson steps alone do not establish reduced total work. Stop each variant at its first failure; later seeds and systems intentionally unrun. No retuning. |

Raw: `P47A-evidence/si8-amended-screen.log`; protocol, timing limitations and
unreached generality gates: [preconditioner-variants.md](preconditioner-variants.md).
## P47B — fp32 H feasibility (2026-11-08; measurement only, Mac observations)

| ID | Attempt | Verdict | Evidence and lesson |
|---|---|---|---|
| P47B-1 | f32/f64 rustfft lanes, real-pair H and sequential faer nonlocal/Davidson shapes; ordinary and native release; physical early/converged Si8/Si16/Si32 snapshots; standard histories with P14A counts | **Do not recommend implementation under the frozen 10% median wall-saving rule** | The actual fused plan is f64-only, so the explicitly allowed lane fallback measures precision ratios, not a production fused speedup. The fp64 experiment matches production H to about 3e-15 Ha; fp32 orthogonal residual noise reaches about 1.2e-6 Ha and projected Ritz errors about 3.1e-7 Ha. Use 1e-5 Ha only as an empirical standard/Ritz eligibility screen; 1e-7 Ha is already below measured root errors. All nine current standard histories converge and independent counts match exactly. Even crediting all eligible diagonalization time to H leaves the median wall ceiling below 10% in both builds. The implemented f32 lane H does not beat production fused f64 H. Full shapes, paired times, exact upper-ceiling percentages, commands and limitations: [fp32-feasibility.md](fp32-feasibility.md), raw `P47B-evidence/`. No production option or source change, no fp32 certificate, no QE or dev08 run. |

This precision experiment does not reopen A11/A28/A34/A43 traversal/threading
proposals. Revisiting it needs a new measured fused-f32 or hardware premise,
not a tighter eligibility threshold chosen to make the bound pass. Davidson
GEMM ratios are reported but not credited to H-only savings. The ceiling is a
fixed-trajectory precision-transfer model, not a theorem about an unimplemented
kernel or a mixed-precision SCF trajectory.
| STR1 | Orchestrator (Codex out of credits): shell and species memoization of the radial pseudopotential transforms in the stress and force passes (`dft-physics` stress local_or_core, nonlocal metric, forces local_or_core): values and derivatives depend only on (potential, bit pattern of |G|), so each is computed once per shell and reused across the atoms of a species; the per-species nonlocal strain metric is built once per potential instead of once per atom | adopted (bit-identical) | Mac observation, si8 strict real Gamma, `[output] forces+stress`: 8.07 s -> 5.04 s total wall (stress alone 8.07 -> 5.63 s after local-term cache, frozen-stress elapsed 3.20 -> 1.52 s); forces and stress results (all terms, ledger) identical to the baseline binary including every digit (JSON comparison without elapsed fields); `scripts/check.sh` all PASS. Profile before: radial transforms (Bessel, erf, exp, mesh integrals) about 52 percent of samples of the whole stress run, 87 percent of frozen_stress. The remaining stress cost is the cutoff-derivative Pulay estimate (three extra strict SCF runs, about 73 percent of the run after this change) | Lesson: the dev08 stress pass (9.7x SCF, relayed 2026-11-08) was dominated by repeated per-atom, per-G radial integrals and by the three scan SCFs, not by band algebra; dev08 confirmation pending. Next: warm-start the two upper-cutoff scan SCFs from the converged density on the shared scan grid (not bit-identical; needs a tolerance gate) and skip the centre scan when the scan grid equals the production grid |
| REG1 | Suspected dev08 default-SCF regression 88cb469 -> 2dda886 (+8.35 percent, si8 strict real, 10 of 10 paired rounds, relayed 2026-11-08): the P47A preconditioner plumbing put a three-way `match` and an extra `finite` check inside the per-element correction loop of the real Davidson; the default damping now takes the original tight loop again (research variants have their own loops) | adopted (default loop restored, bit-identical) | Mac observation: pre-fix and fixed binaries alternate 8 rounds on si8 strict real, medians 1.010 s and 1.008 s (paired difference 0.002 s): the regression does not reproduce on the Mac, so the cause is NOT established; `bench_scf --require-identical` gamma and si8 identical; check.sh PASS. dev08 measurement requested (88cb469 vs this commit, and bisect of the commits in between if the +8 percent remains) | Lesson: a research option added to a hot loop must not change the default loop's instruction stream; verify on the acceptance machine, the Mac did not show it. **dev08 bisect (relayed 2026-11-08, 5 commits, 10 alternating rounds): the first bad commit is b08dc90 (Merge P47A), +7.74 percent paired median, 10 of 10 rounds; P47B and STR1 are noise; REG1 (bffbba2) recovers 6.04 points (0 of 10 slower); residual +1.38 percent vs 88cb469 (8 of 10 rounds, min -0.86), weak evidence, not chased (possible x86 code-layout sensitivity, as in P30A where codegen-units=1 removed it).** |
| STR2 | Orchestrator: cutoff-derivative Pulay scan of the stress pass runs its three SCF energies in the standard profile (residual threshold 1e-9 Ha) instead of strict, with an automatic strict repeat if a scan fails to converge. Derivation (before measuring): the estimate is 2c/(3V) dE/dc with dc = 0.75 and 1.5 Ha, so a scan energy error of 1e-9 Ha moves it by about 1e-11 Ha/Bohr^3 against a reported uncertainty of 1e-7 and a warning threshold of 1e-5 | adopted | Mac observation, strict production run with `[output] stress`: si8 Pulay assessment 3.898 s -> 1.497 s, whole run 5.04 s -> 2.93 s; Pulay estimate -3.933884e-6 -> -3.933873e-6 (difference 1.1e-11), uncertainty 2.8832e-7 -> 2.8831e-7; gamma (2 atoms) estimate difference 2.5e-13, wall 0.33 -> 0.25 s; k222 -3.179039704e-6 -> -3.179039950e-6, wall 2.39 -> 1.61 s; k222s -4.78346977e-6 -> -4.78346958e-6, wall 1.28 -> 0.86 s; frozen stress tensor and pressure bit-identical; `scripts/check.sh` PASS. Rejected in the same step: density warm start of the second and third scan from the previous converged density (Pulay time 3.90 -> 3.76 s only, -4 percent: the strict cold Davidson and certificate dominate, not the density iterations) and skipping the centre scan when the scan grid equals the production grid (never applies: si8 production 36^3, scan grid larger) | The estimate is not bit-identical any more (about 1e-11 Ha/Bohr^3, four orders below the stated uncertainty); the frozen stress and forces are untouched. Skewed, slab and spin fixtures were not run on the Mac: the strict fallback covers non-convergence, dev08 confirmation pending |
| PRT1 | Orchestrator: partial rotation in the real Davidson (rotate the Ritz pairs and form residuals only for roots that can still change a decision, complete the rest on demand before finish, certificate and restart; env-gated prototype, not merged) | rejected | Mac observation si8 real: strict profile slower (1.061 s vs 0.998 s) and a different trajectory (first SCF density residual 0.026204 vs 0.026189, final energy differs at 2.5e-14 Ha): converged roots do not stay converged (their residuals rise again as directions are added), so stale flags change the selection; standard profile (Ritz-change flags are exact) same Davidson iteration counts per SCF step but energy differs by 2e-14 Ha (compact GEMM groups columns differently) and wall only 0.386 -> 0.370 s (-4 percent) | Lesson: the all-roots rotation is not the dominant per-step cost on these shapes (si16 strict real profile on the Mac: GEMM about 45 percent in total, nano_gemm 4x4 micro-kernels 25 + 9 percent, FFT 17 percent, projection 7 percent, nonlocal 5 percent), and exact-trajectory preservation needs fresh residuals for every root in strict mode |
| ROT1 | Orchestrator: blocked rotation kernel for the real Davidson (Q C and HQ C): 8 x 4 register tile, right operand packed in blocks of four, one ascending-p accumulation chain per output element (independent of row partition and thread budget), fused multiply-add where the build has it; opt-in switch `RotateKernel::Blocked` (`BAANDI_ROTATE_KERNEL=blocked`, default Legacy) | rejected on dev08 (code kept as an off-by-default switch) | Trigger: dev08 call-site profile (2026-11-08, df4a485): `real::rotate` is 37.2 percent of si32 strict and 27.5 percent of si16 strict samples, the leaf is the 1x4 faer nano-GEMM micro-kernel, derived rate about 2.3 GFLOP/s against 85 GFLOP/s of single-core MKL dgemm at the same shapes (m = 5939, k = 408, n = 68). Mac observation (aarch64, where the legacy NEON micro-kernel is good): legacy 29.6 to 31.0 GFLOP/s, blocked 16.7 to 17.5 GFLOP/s at the si16 and si32 shapes; so the Mac cannot show the gain and the blocked kernel is slower there. Gates passed: blocked equals legacy within 4 k eps relative, blocked results bit-identical across thread budgets 1, 2, 4, 8, 12 and row partitions; default path `bench_scf --require-identical` gamma and si8 identical; check.sh PASS | Blocked changes the summation order, so SCF trajectories differ at rounding level (Mac si8: 14 vs 15 SCF steps in an earlier prototype); adoption on dev08 needs energy within the strict tolerance, equal or lower SCF and Davidson step counts and a wall-time gain on si16 and si32 | **dev08 result (relayed 2026-11-08, ba4310a, same binary, env switch, alternating rounds): microbenchmark legacy 11.1 to 16.9 GFLOP/s, blocked 6.1 to 7.6 GFLOP/s (blocked 0.45 to 0.64x); full SCF blocked slower in 26 of 26 paired comparisons (si8 strict +10.8, standard +3.8; si16 strict +16.9, standard +8.5; si32 standard +10.0 percent); si32 strict +107.6 percent because the trajectory changed (SCF 11 to 14, Davidson 217 to 411, H 4669 to 6440) from rounding-level differences; energies within the strict tolerance everywhere. Correction of the profile derivation: legacy rotate runs at about 17 GFLOP/s at the si32 shape (measured), not 2.3 (derived from a one-rotation-per-step flop model; the real rotate work is about 7 times that model), so the gap to single-core MKL dgemm (85 GFLOP/s) is about 5x, not 37x; the sampled shares (rotate 37.2 percent, subspace algebra 57.4 percent of si32 strict) are unaffected.** | Lesson: auto-vectorised plain Rust does not beat the faer/gemm nano micro-kernels on either NEON or AVX-512; a faster rotation needs real SIMD micro-kernels or MKL dgemm (backlog design R8, whose trigger of 25 percent algebra share is now met: 50 to 57 percent on strict) |
| MKL-R1 | Orchestrator: real-path rotation (Q C and HQ C of the real Davidson) through single-thread MKL dgemm under the `mkl` feature (`dft_mkl::dgemm`, new FFI adapter function; `RealGemmWorkspace::rotate_pair_slices`, 256-row panels each one dgemm into a private tile copied into place; panel grouping fixed, so the result does not depend on the thread budget); without the feature the slice entry point is the ordinary `rotate_pair` (bit-preserving) | opt-in by build feature, dev08 evidence pending | Trigger (backlog design R8, trigger 25 percent met): dev08 profile 2026-11-08, subspace algebra 50 to 57 percent of strict, `real::rotate` 37.2 percent of si32 strict; measured legacy rotate about 17 GFLOP/s vs 85 GFLOP/s single-core MKL dgemm at the same shape. Mac: no MKL, type check only (`cargo check --features dft-numerics/mkl`, clippy with the feature `-D warnings`), non-MKL default path `bench_scf --require-identical` gamma and si8 identical, check.sh PASS | MKL dgemm results are not bit-identical to the faer path and are only reproducible for fixed shapes and alignment (no CNR setting yet): adoption needs energy within the strict tolerance, SCF and Davidson step counts not above the default, a wall gain on si16 and si32, and run-to-run identical density hashes |
| MKL-R1 dev08 | dev08 pure-vs-MKL A/B at 6f4bc31 (relayed 2026-11-08): determinism 4/4 (density SHA identical run-to-run and T=1 vs T=8 for both builds), energies within the tolerance and MKL error never larger except si16 standard; but Davidson/SCF counts differ in all six rows (MKL FFT and BLAS take a different trajectory already before MKL-R1) and wall is +26.6 percent on si32 strict (Davidson 217 -> 418, per step MKL 1.5x faster) while si16 strict is -4.85 percent (Davidson 242 -> 175, per step 1.3x slower): the pure-vs-MKL comparison cannot isolate the kernel | not decided; same-build switch added: `BAANDI_MKL_DGEMM=off` disables the dgemm rotation inside an MKL build (default on), so the next A/B compares identical trajectories except for rotation rounding | Lesson: compare a kernel inside one build; a cross-backend A/B mixes trajectory length with per-step cost |
| MKL-R1 A/B | dev08 same-binary A/B at 560c73b (relayed 2026-11-08, one MKL native binary, `BAANDI_MKL_DGEMM` off vs on, alternating rounds, T=1, n=5/5/3): dgemm rotation faster in 26 of 26 paired rounds: si8 strict -8.25, standard -4.28; si16 strict -9.55, standard -5.04; si32 strict -11.87, standard -6.37 percent; five rows bit-identical (energy, SCF/Davidson/H counts, density SHA); density SHA stable run-to-run and T=1 vs T=8 for both settings (4/4); the exception is si32 strict where the trajectory diverges (Davidson 342 -> 418, H 5857 -> 5945, SCF 13 = 13, dE 3.456e-11 Ha, |dE| vs tight 1.663e-11 -> 5.119e-11; per Davidson step -27.7 percent) | frozen rule not met on one row (counts above the off run, and dE beyond 2e-12 Ha), therefore NOT recommended as default: the dgemm rotation is now opt-in (`BAANDI_MKL_DGEMM=on`, MKL builds only) until a distribution check shows whether the si32 divergence is systematic or one rounding-level trajectory (the known 1e-10 sensitivity of the standard and strict Davidson paths) | The kernel gain is clear (4 to 12 percent, bit-identical results where trajectories agree); only trajectory chaos stands between it and adoption |
| MKL-R1 dist | dev08 perturbed-start distribution at 91aa2b3 (relayed 2026-11-08, one MKL native binary, dgemm rotation off vs on, same perturbed inputs: atoms shifted by +-1e-7 Bohr, seeds 501 to 510 (si16) and 501 to 505 (si32), T=1): si16 strict 10 starts Davidson 281/303/367 on both sides, all energies identical, median wall 12.931 -> 11.330 s (-12.38 percent); si16 standard 10 starts Davidson 69/79.5/94 identical, wall 4.365 -> 4.128 s (-5.43 percent); si32 strict 5 starts: Davidson more with on in 2 (+25, +31), equal in 2, fewer in 1 (-1), median 582 -> 581, SCF equal in all, median wall 97.631 -> 79.351 s (-18.72 percent); all 25 runs converged, max |dE| vs tight 1.468e-11 (strict), 4.399e-09 (standard) | adopted as default for MKL builds (frozen follow-up rule 4 of 4 met); `BAANDI_MKL_DGEMM=off` is the switch back | The single si32 strict divergence (342 -> 418) was one tail of a mixed-direction distribution, not a systematic effect; gain grows with trajectory length (-11.9 percent unperturbed, -18.7 percent perturbed) |
| MKL-R2 | Orchestrator: the same single-thread MKL dgemm panels for the rest of the real Davidson algebra under the `mkl` feature: Gram Q^T X (reduction panels of 256, tiles shaped m x n summed in the fixed `assemble` order), projection products Q (Q^T X) (256-row panels) and the reduced-matrix update V^T HV (new `RealGemmWorkspace::multiply_slices_into`, slice inputs with leading dimensions; non-MKL builds forward to `multiply_into` on strided views, bit-preserving); same switch `BAANDI_MKL_DGEMM=off` | adopted for evidence collection (default on in MKL builds together with MKL-R1; dev08 distribution evidence pending) | Trigger: the dev08 profile buckets `real::gram` 15.1 percent, `real::rotate` 37.2 percent (it also contains the projection Identity products, which MKL-R1 had not moved) and `reduced_matrix_into` 1.6 percent of si32 strict. Mac: type check only for the MKL path (cargo check and clippy with `--features dft-numerics/mkl`, `-D warnings`), default path `bench_scf --require-identical` gamma and si8 identical, check.sh PASS | Needs the same dev08 same-binary off/on distribution (si16 strict and standard, si32 strict, perturbed starts) and the frozen rule: median Davidson steps not above off, all converged, |dE| within class, median wall gain; fall back to off if it fails |
| MKL-R2 dist | dev08 MKL-R2 measurement at 863b47a (relayed 2026-11-08, one MKL native binary, `BAANDI_MKL_DGEMM=off` disables R1 and R2 together, perturbed inputs shared, T=1): combined R1+R2 vs neither, perturbed distribution: si16 strict 10 starts Davidson 281/303/367 identical on both sides, median wall 12.900 -> 10.550 s (-18.22 percent); si16 standard identical counts, 4.357 -> 4.021 s (-7.70); si32 strict 5 starts Davidson median 582 = 582 (more with on in 1 start, fewer in 1), wall 98.213 -> 71.845 s (-26.85); all converged, max |dE| 3.011e-11 (strict) and 4.399e-09 (standard); determinism: density SHA identical run-to-run and T=1 vs T=8 for both settings. R2-only increment (this run's on vs the R1 distribution run, same seeds): si16 strict -6.9, si16 standard -2.6, si32 strict -9.5 percent. Unperturbed single starts (n=5/5/3): wall -7.9 (si8 strict), -6.5 (si8 std), -1.5 (si16 std), -25.5 (si32 strict), -6.8 (si32 std) percent but si16 strict +9.75 percent because that single trajectory lengthens (Davidson 175 -> 270, SCF 13 -> 15), and the unperturbed si32 strict |dE| 1.238e-10 is just outside the 1e-10 class | adopted: R1+R2 default for MKL builds (frozen rule 5 of 5 met on the perturbed distributions; the rule was set on distributions); the unperturbed si16 strict reversal and the 1.24e-10 single value are recorded as the known tail behaviour | Lesson: single unperturbed starts are one draw of a chaotic trajectory distribution (known 1e-10 sensitivity); judge algebra changes by distributions; the combined kernel gain is 8 to 27 percent on MKL builds |
| FIN1 | Orchestrator: remove the redundant input finiteness scans of the real Davidson products (`rotate_pair`, `rotate_pair_slices`, `multiply_into`, `multiply_slices_into`: every call read all of A and B once before the product, about 20 MB per si32 rotation); a non-finite input propagates to every dependent entry of C (NaN and infinity survive products and sums), so the existing output scans still reject it with the same `NonFinite` error | adopted (bit-identical results) | Trigger: dev08 S8 split (relayed 2026-11-09, MKL strict, line-level sampling with inlined frames): dgemm input scans `real.rs:330/331` 7.6 percent of the si32 wall (5.0 percent si16), output scans 2.7 percent, finiteness checks in total 10.3 percent of the si32 wall; MGS dot 5.6 and subtraction 2.6 percent of the wall. Mac observation: si16 strict real, alternating 3 rounds, 4.108 s -> 3.859 s (-6.1 percent), total energy bit-identical; `bench_scf --require-identical` gamma and si8 identical; dft-numerics and dft-solvers test suites pass | dev08 confirmation requested; the MGS dot/axpy loops (about 8 percent of the wall) are the next item (blocked orthonormalization) |
| FIN1 dev08 | dev08 FIN1 verification (relayed 2026-11-09, c000d11 against ee86f07, same flags, pure and MKL, T=1, alternating, n=5/5/3): energy, SCF/Davidson/H counts and density SHA identical in 12 of 12 rows; wall pure: si8 strict -10.7, std -4.0; si16 strict -8.1, std -2.2; si32 strict -3.6, std +1.26 percent (0 of 3 rounds faster); MKL: si8 strict -15.7, std -6.4; si16 strict -14.4, std -6.9; si32 strict -13.9, std -7.3 percent. Frozen rule (no case slower) missed by the single row pure si32 standard (+1.26 percent, identical numerics, so a code-placement effect of the P30A kind); the other 11 rows gain 2.2 to 15.7 percent | kept (already on main; numerics unchanged), codegen-units check of the +1.26 percent row requested | The dgemm call-site timers (MKL, strict, FIN1): si32 rotate 40984 calls 5.108 s (11.4 percent of the wall), projection_overlap 39574 / 1.956 s (4.4), reduced 20492 / 1.227 s (2.7), projection_apply 31396 / 0.908 s (2.0): total 20.5 percent of the wall, equal within 0.1 point to the sampled 20.6 percent; si16 15.9 percent; the timers carry no gram calls because the Gram products are tagged by their caller (`projection_overlap`, `reduced`) |
| FIN1 cgu | dev08 check of the single regression row of FIN1 (relayed 2026-11-09): pure si32 standard, ee86f07 vs c000d11, both with `-C codegen-units=1`: 17.071 s vs 16.710 s (-2.11 percent, 3 of 3 rounds faster) against +1.26 percent with the default codegen units: the sign flips, so the +1.26 percent was a code-placement effect and not a FIN1 effect (same family as P30A; codegen-units=1 itself is slower than the default build for both commits, 16.71 vs 16.19 s) | FIN1 adopted without reservation | The frozen rule's 'no case slower' was triggered only in the default build by placement noise; future rules for bit-identical changes should be judged on the median over all rows and re-checked with codegen-units=1 when a single row flips |
| DGEMM-SITES | dev08 per-site dgemm rates with `gflops_total` (relayed 2026-11-09, MKL, strict, T=1): si32 rotate 40984 calls 5.126 s 298.8 GFLOP 58.3 GFLOP/s (70 percent of the shape ceiling), projection_overlap 39574 / 1.941 s / 53.7 GFLOP 27.7 GFLOP/s (27 percent), reduced 20492 / 1.224 s / 37.5 GFLOP 30.7 (34 percent), projection_apply 31396 / 0.899 s / 8.72 GFLOP 9.7 (12 percent); total 132446 calls 9.19 s, 20.5 percent of the wall, ceiling pace 4.62 s (10.2 percent of the wall of slack, 6.4 percent in the three small sites); si16 rotate 57.6 (75 percent), overlap 21.3, reduced 23.0, apply 13.1 GFLOP/s, total 15.8 percent of the wall; efficiency follows the work per call monotonically (7.29 MFLOP per call 70 percent, 1.83 34, 1.36 27, 0.28 12) | open: larger work per dgemm call (bigger fixed panels) is the candidate; needs the panel-size dgemm sweep and a T=8 check because fixed panels set the thread-independent task granularity | The 256-row panel is a thread-independence constant, not a performance optimum |
| DGEMM-CENSUS | dev08 census of the real dgemm calls by an LD_PRELOAD interposer on cblas_dgemm (relayed 2026-11-09, gemmspy, outside the repository): exactly 62232 (si16) and 132446 (si32) calls, equal call for call to the dgemm-timer counts; two families, beta = 0: NN m = 256 (rotate 40984 + projection_apply 31396 calls, 307.5 GFLOP, 51.4 GFLOP/s) and TN with reduction panel 256 (projection_overlap 39574 + reduced 20492, 91.2 GFLOP, 29.0 GFLOP/s); lda = 5915 (si16) and 11877 (si32) because the real basis is 2 ngw - 1 rows (the earlier microbenchmark halved it by mistake: the panel-size sweep at d = 2958/5939 is void). MKL cliff in the Trans-A case at exactly k = m: for m = 204, n = 68 TN 84.3 GFLOP/s at k = m, 58.2 at k = m + 1 (m = 136: 87.0 / 59.8; m = 144, n = 36: 77.2 / 43.5); every TN call with output rows m < 256 and reduction panel 256 sits on the slow side; k panel = largest multiple of 8 not above m gives x1.22 to x1.48 for m = 72 to 204, ~1.0 for m <= 68 or m >= 263; weighted gain 0.451 s on si32 (1.0 percent of the wall), 0.118 s on si16 (1.26 percent); suggested rule panel_k = m >= 128 ? round_down_8(min(m, 256)) : 256. NN panel: flat within 5 percent from 64 to 4096 rows, 256 near the best: keep. Three regimes at the exact production shapes (si32 strict, dgemm 9.136 s): packed 4.566 s, hot (real lda) 5.043, walk (cycling panels as in production) 5.677, production 9.136: 3.459 s (7.7 percent of the wall) is neither shape, stride nor panel cache but in-process (operand eviction and per-call cost outside MKL); 25 percent of the si32 projection_apply calls have n = 1 (single-band, latency-bound; only batching helps) | corrections: the 98.5 GFLOP/s square ceiling is not reachable at these shapes (NN 76 to 88, TN 42 to 59 at k = 256, n = 1 calls 4 to 6 GFLOP/s), realistic si32 headroom 3.46 s (7.7 percent of the wall), not 4.57 s; the removed finiteness input scan acted as a prefetch for the dgemm that followed (walk+scan 82.9 vs plain walk 77.4 GFLOP/s at m = 256, n = 68, k = 407), so FIN1 nets 4.15 s of 4.62 s on si32 (9.3 percent of the wall, not 10.3) | candidates: (1) TN reduction panel rule (about 1 percent), (2) batching single-band projection_apply calls, (3) the 7.7 percent in-process residual; panel-size end-to-end builds (128/256/512/1024, T=1 and T=8) are running on dev08 | The 'ceiling' used so far was a square-matrix number; use per-shape standalone rates at the production lda |
| PANEL-INSITU | dev08 in-situ panel builds from a /scr clone of 12fa159 (relayed 2026-11-09, strict, interleaved, 3 rounds; PANEL 128/256/512/1024 and an adaptive Transpose reduction panel `m >= 128 ? round_down_8(min(m,256)) : 256`): counts differ and density SHA differs in all 16 comparisons (the panel changes the summation order and therefore the trajectory; si16 PANEL=512 -25 percent is 13/213 iterations instead of 15/270); per H application vs 256: si16 128 +1.4, 512 -11.9, 1024 -6.6, adaptive -14.9 percent; si32 128 +9.1, 512 +0.2, 1024 +3.8, adaptive +9.0 percent; dgemm GFLOP/s flat across all variants (si16 37.4 to 39.0, si32 41.6 to 44.0), dgemm share of the wall flat (13.9 to 15.7, 19.1 to 21.4 percent); the adaptive rule moves the bulk of the TN time to the fast kernel (si16 +23 percent TN GFLOP/s, upper bound) but on si32 the call count rises 60066 -> 101507 and the tile scratch doubles (TN -3 percent) because 35 percent of si32 TN calls already have m >= 256; T=8/T=1 speedup is independent of the panel (si16 1.59 to 1.64, si32 1.90 to 1.97), so panel count is not what limits thread scaling | rejected (no panel change; no variant meets the frozen rule, none supports a default change) | Lesson: a kernel-level cliff does not translate into a gain when the call count and tile traffic grow; the real residual is in-process data flow (3.46 s, 7.7 percent of the si32 wall) and the 1.6x/1.9x T=8 scaling |
| MGS-B | Orchestrator: blocked Gram-Schmidt for the correction block of the real Davidson (`BAANDI_MGS=blocked`, default off): panels of 16 columns, two block projections against the kept columns by GEMM (`multiply_slices_into`, site Other) and the unchanged two-pass MGS and rank floor inside the panel; replaces the O(n k^2) scalar dot/axpy chain that is 5.6 + 2.6 percent of the si32 wall and is serial (so it also limits T=8 scaling) | opt-in, dev08 evidence pending | Mac observation (pure faer GEMM path, not representative of MKL): slower, si8 0.916 -> 1.035 s, si16 3.763 -> 4.266 s, and the trajectory changes (SCF 14 -> 15, rounding level); dft-solvers tests pass with the switch on; the case rests on MKL dgemm at n = 16 columns | Needs the dev08 same-binary switch comparison on MKL and pure builds (T=1 and T=8, perturbed distributions); default stays old |
| MGS-B dev08 | dev08 same-binary comparison of `BAANDI_MGS=blocked` at e6d9b9f (relayed 2026-11-09; default path bit-identical to 12fa159, si16 strict MKL 15/270/3504 SHA c7a98de81903f466): unperturbed wall off -> blocked, MKL T=1: si16 strict -18.4, std -11.7; si32 strict +15.7, std -9.3 percent; T=8: si16 strict -21.0, std -12.9; si32 strict +2.2, std -16.1; pure builds lose in 5 of 8 cells (si32 strict T=1 +44.7); perturbed distributions MKL T=1: si16 strict (seeds 501-510) wall -1.09 percent, 8 of 10 faster, Davidson median 303.0 -> 303.5; si16 standard -0.81 percent, 79.5 -> 79.5; si32 strict (501-505) -3.49 percent, 5 of 5, 582 -> 585; all converged; cost per H application falls with blocked in all 8 T=8 cases (MKL si16 strict -6.0, std -3.4, si32 strict -9.4, std -12.7; pure si32 strict -11.0, std -19.1 percent) and the T=8/T=1 scaling improves in all 8 (MKL si32 strict 1.91 -> 2.16, pure si32 strict 2.34 -> 2.82); |dE| vs tight smaller with blocked in 6 of 8 cells (MKL si32 strict 1.238e-10 -> 1.160e-11); density SHA stable run-to-run and identical T=1 vs T=8 in all blocked cells | rejected as default under the frozen rule for both builds (Davidson step medians above off, T=1 wall criterion contradictory between unperturbed and distribution measures, si32 strict +2.2 percent at T=8, pure losses); kept as the opt-in switch | The mechanism is confirmed (the serial MGS is a T=8 bottleneck: cost per step -3 to -19 percent at T=8, scaling +13 to +21 percent on si32) but rounding-level trajectory movement is larger noise than the gain at T=1; the T=8 distributions are the open measurement |
| MGS-B T=8 | dev08 MGS-B T=8 distributions (relayed 2026-11-09, MKL, same binary, env unset vs `BAANDI_MGS=blocked`, taskset 0-7, 30 starts x 2 interleaved rounds): frozen rule FAILS on one of six criteria, the Davidson median (si16 strict 303.0 -> 303.5, si32 strict 582 -> 585; si16 standard 79.5 = 79.5, si32 standard 95 = 95), passes the other five: 60 of 60 runs converged, |dE| max 1.123e-12 (strict) and 5.684e-13 (standard), median wall gain si16 strict -3.04, si16 standard -1.83, si32 strict -12.40, si32 standard -11.15 percent, no case worse than +1 percent, density SHA stable for all 30 blocked seeds in both rounds; Davidson counts identical on 25 of 30 seeds (differences: si16 501 367 -> 389, 502 286 -> 302; si32 503 476 -> 457, 504 582 -> 585, 505 586 -> 651); restricted to the identical-work seeds the wall gain is si16 strict -3.22 (8 of 10 seeds), si16 standard -1.82, si32 strict -12.42 (2 of 5), si32 standard -11.15 percent, spreads of a few tenths of a percent | frozen rule failed: no adoption on this data; a new rule on fresh seeds is pre-registered (see the next row) | The cost-per-step effect is large and clean (about -11 to -12.5 percent of the si32 wall at T=8 and identical work); only the trajectory tail (5 of 30 starts) breaks a median rule that was written as 'not above off' |
| MGS-B R2 | Pre-registered re-test of MGS-B on FRESH seeds (set by the orchestrator on 2026-11-09 before any new data; the first rule failed and is not changed retroactively, this is a new test on independent data): MKL build, same binary, env unset vs blocked, perturbed starts seeds 601-610 (si16 strict and standard) and 601-605 (si32 strict and standard), T=8 (taskset 0-7) and T=1 for the two strict cases; adopt blocked as the default for MKL builds only if ALL hold: median Davidson steps of blocked at most 1.03 x off on each case; blocked has more Davidson steps on at most 40 percent of the seeds of each case; all runs converge; every |dE| within the class (strict 1e-10, standard 5e-9); median wall gain at T=8 at least 3 percent on si16 strict and si32 strict and no case median worse than +1 percent at T=8 or T=1; density SHA stable run-to-run (two rounds) and identical between T=1 and T=8 for blocked | pending dev08 data | The 1.03 tolerance and the 40 percent seed bound encode that medians over 5 to 10 seeds cannot resolve half a Davidson step; they are fixed now and not to be adjusted after seeing the data |
| MGS-B R2 result | dev08 MGS-B R2 on fresh seeds (relayed 2026-11-09, MKL binary cad43dfb at e6d9b9f, env unset vs `BAANDI_MGS=blocked`, two interleaved rounds, 120 runs): pre-registered rule PASSES 7 of 7: median Davidson blocked/off 1.0000 in all cases except si32 strict 1.0103 (T=8 and T=1); blocked has more Davidson steps on 0 percent of seeds except si32 strict 40 percent (2 of 5, exactly the bar); 120 of 120 converged; max off-vs-blocked |dE| 2.612e-11 (strict), 3.268e-13 (standard); median wall at T=8: si16 strict 5.767 -> 5.566 s (-3.10 percent, 10 of 10 seeds), si16 standard -2.27, si32 strict 31.432 -> 27.428 (-12.68, 5 of 5), si32 standard -11.47; at T=1: si16 strict -0.68, si32 strict -4.07 percent; density SHA stable across rounds for all 30 blocked runs and identical T=1 vs T=8 in all 15 pairs; Davidson and H counts bit-identical on every seed of the other four cases (so those wall numbers are pure cost per step); cross-round agreement with R1 within 0.5 points. Master's flags: the 40 percent criterion has no resolution at n = 5 (2 of 5 passes, 3 of 5 fails), si16 strict at T=8 clears its 3 percent bar by 0.10 points, and against the tight reference si32 strict seed 604 is outside 1e-10 on both paths (off 1.359e-10, blocked 1.097e-10) | adopted: blocked Gram-Schmidt is the default for MKL builds only (`dft_numerics::MKL_BACKEND`), `BAANDI_MGS=classic` switches back, pure builds keep the classic loop (blocked loses there) | The gain on the equal-core relevant case is large and replicated (si32 -11 to -12.7 percent at T=8, -4.1 at T=1); the n = 5 resolution limit of the 40 percent criterion and the 0.10 point margin of si16 strict at T=8 are recorded as the weak points |
| MGS-B post-adoption | dev08 confirmation and production check at 43bed40 (relayed 2026-11-09, MKL): switch mapping verified (si16 strict T=1 unset = blocked 13/210/2942 SHA 3483c9fb7613348a, classic 15/270/3504 SHA c7a98de81903f466); confirmation si32 strict T=8 fresh seeds 606-615, classic vs default: median wall -21.64 percent (10 of 10 faster, per-seed -20.6 to -22.7), Davidson median 615.0 -> 612.5 (ratio 0.996), default has more steps on 1 of 10 seeds (the 40 percent criterion resolved at n = 10 in favour of the default), all converged, SHA stable, max |dE| 2.558e-12; unperturbed production check within the commit, classic -> default: si16 strict T=1 -22.87, T=8 -26.04; si16 standard -14.85, -17.93; si32 strict +8.38 (0 of 3; that start takes 14/464/6512 instead of 13/330/5777, its pre-adoption trajectory), T=8 -8.48; si32 standard -14.15, -23.43 percent. Layout effect of the adoption commit on the now-off classic path: e6d9b9f classic vs 43bed40 classic +4.77 (si16 strict T=1), +7.53 (T=8), +6.32 (si32 strict T=1), +11.83 (T=8), +6.42 and +9.56 percent (si32 standard), counts, energies and SHA identical (dE 0), and the default path itself is unchanged (e6d9b9f blocked vs 43bed40 default -0.68 to +0.98 percent) | adopted, not reverted. The honest user-facing gain is pre-adoption classic (e6d9b9f) to the new default (43bed40): si16 strict -18.6 (T=1) and -21.2 (T=8), si32 strict +15.1 (T=1, single start) and +1.9 (T=8, single start), si32 standard -8.7 and -16.3 percent; on the perturbed-start distributions the gain is -21.6 percent (si32 strict T=8) | The within-commit percentages overstate the gain because the baseline path slowed by 5 to 12 percent through code placement; unperturbed single starts of si32 strict are one draw of a chaotic trajectory distribution and flat or worse, the distributions are the reliable measure; a codegen-units=1 rebuild of both commits would show whether the classic gap is placement only (low priority, the classic path is off by default) |
| AUD1 | Orchestrator: band-parallel nonlocal term of the energy audit (`energy::evaluate_shared` and `evaluate`: per-band Re<c|V_nl|c> split over the worker pool with `WorkerPool::join`, each band's own arithmetic unchanged, accumulation in ascending band order with the occupation weights, so results are bit-identical to the serial loop) | adopted (bit-identical), dev08 evidence requested | Trigger: dev08 stage profile at 43bed40 (relayed 2026-11-09, MKL, si32): `total_energy` does not scale with threads (strict 4.197 s at T=1, 4.075 s at T=8, 17.4 percent of the T=8 wall; standard 3.289 -> 3.106 s, 37.9 percent), average working threads 2.77 (strict) and 1.90 (standard) of 8. Cause found by Mac sampling (si16 strict): the time inside `evaluate_shared` is almost entirely the serial per-band `NonlocalOperator::apply` of the nonlocal energy (the dev08 call-site classifier had filed it under FFT because the inlined closure sits in `MklFftPlan3d::EnergyAudit::step_with_occupations`). Mac observation si16: T=8 `total_energy` phase 0.238 -> 0.069 s, wall 1.744 -> 1.575 s (-9.7 percent), total energy and input hash bit-identical, T=1 unchanged (3.757 vs 3.743 s); `bench_scf --require-identical` gamma and si8 at T=8 identical; check.sh PASS | dev08: same-binary comparison not available (no switch), compare against 43bed40 MKL and pure at T=8 and T=1, si32/si16, strict and standard, with density SHA identical expected | The remaining audit work (shared terms, Hartree/XC FFTs) is small; the next serial or idle part at T=8 is inside the diagonalization (2.77 of 8 threads working) |
| AUD1 dev08 | dev08 AUD1 verification at 886a0fe against 43bed40 (relayed 2026-11-09, MKL and pure, paired alternating, n = 5/5/3): 16 of 16 rows bit-identical (counts, dE 0, density SHA, converged); T=8 total_energy: mkl si16 strict 0.557 -> 0.257 s (-53.9 percent), standard 0.444 -> 0.199 (-55.3), si32 strict 4.102 -> 1.064 (-74.1), standard 3.101 -> 0.809 (-73.9), pure si32 strict 2.924 -> 0.777 (-73.4), standard 3.195 -> 0.837 (-73.8); wall T=8: mkl si16 strict -6.20, std -11.44, si32 strict -12.31, std -27.64; pure -7.09, -13.24, -11.98, -25.33 percent; diagonalization unchanged; T=1: six of eight rows flat to better (-0.23 to -0.40 percent), two MKL rows nominally slower (si16 strict +0.24, standard +0.36 percent); re-measured at n = 15 paired: -0.19, +0.08 and, for the si32 standard control, +0.34 percent, signs do not reproduce, wall sd 0.006 to 0.031 s (noise floor 0.2 to 0.3 percent), T=1 total_energy phase unchanged (-0.95, +0.07, +0.16 percent) | frozen rule: criteria 1 and 3 pass, criterion 2 ('median wall not slower on any row') reads FAIL as written on rows that change from sample to sample at or below the noise floor; AUD1 is kept (already on main) under the lesson recorded after FIN1: bit-identical changes are judged on all rows together and a flipping row is re-measured; no rule edit | Classifier corrections from the master: the audit is not FFT (leaf analysis: scalar arithmetic 35 to 38 percent, memcpy/memset 17 to 18, FFT kernels 4.6 to 5.1, `hamiltonian::add` named) because the generic `ordered_chunks::<.., MklFftPlan3d, ..>` frame matched the FFT rule; `callsite.py` now has `c5_dgemm_panels` and skips that wrapper, re-bucketed 43bed40 stacks: si32 standard T=1 a_fft_local 49.2 -> 30.5, e_audit 7.3 -> 26.0 percent (phase timer 23.4); at T=8 the audit share must be read from the phase timer, not the sampler |
| GR-AUTO | Orchestrator: `gamma_real` default changed from off to auto (ADR 0040, user approval 2026-11-09); engine and CLI defaults, seven test files updated (default expectation, legacy checkpoint header pinned to off, result-key counts 24 -> 26, spin and spectral-guard comparisons pinned to the complex path) | adopted (user-approved default change) | check.sh PASS on the Mac after the test updates; dev08 default-vs-off measurement requested | Users running defaults now get the real path on qualifying Gamma cells; baselines for default-path identity use `gamma_real = "off"` explicitly |
| GR-AUTO dev08 | dev08 verification of the new default (relayed 2026-11-09, d836731): default (auto) vs explicit off, 24 runs: real path used everywhere, speedups 1.1 to 3.2x (MKL si32 strict T=8 66.1 -> 20.6 s, 3.22x), all in the accuracy classes; multi-k falls back with a stated reason and bit-identical results to off; spin fallback pending (input error on the measurement side) | confirmed; results.json now always records gamma_real, gamma_real_used and gamma_real_fallback_reason (the `!= Off` guards hid the explicit-off case once the default became auto) | The default is now the fast path for single-Gamma users |
| MKL-EVD | Orchestrator: MKL `dsyevd` (all eigenpairs, lower triangle, workspace queried once) for the reduced eigenproblem of the real Davidson (`dft_mkl::SymmetricEigen`, `RealEighWorkspace`; the sign normalisation of `symmetric_eigh` is applied afterwards, so the eigenvector convention is unchanged); dev08 microbenchmark (relayed 2026-11-09): dsyevd 0.54 to 0.62 x faer at orders 96 to 408 (si32 order 408: 18861 -> 11679 us), dsyevr with only the lowest b pairs slower than dsyevd for every order this solver uses (crossover order/b about 9 to 10; subspace is at most 6b) so dsyevr was not built | adopted as the MKL-build default after dev08 pre-registered distribution test (opt-in at d836731, fresh seeds 701-710 si16, 701-705 si32, two interleaved rounds, same MKL binary, `BAANDI_MKL_EVD` unset vs on): 7 of 7 criteria pass: median wall si16 strict T=8 -5.06 (10/10), T=1 -2.90 (10/10), si16 standard T=8 -1.26 (10/10), si32 strict T=8 -5.66 (5/5), T=1 -2.38 (5/5), si32 standard T=8 -1.65 (5/5); Davidson counts bit-identical on 28 of 30 seed-cases, si32 702 652 -> 650 and 703 578 -> 547 (median 578 -> 560, ratio 0.969); all converged, max |dE| 9.8e-13 (strict) and 1.0e-12 (standard), density SHA stable across rounds and identical T=1 vs T=8 on all 15 strict seeds; measured gain matches the microbenchmark-based projection (5.5 to 6.1 percent at T=8, 2.2 to 2.5 at T=1) | `BAANDI_MKL_EVD=off` switches back to faer; the non-MKL build is unchanged; the sizing method (microbenchmark plus profile share) predicted the result, usable for the next candidate |
| OFFSET-G0 | Orchestrator, from the dev08 diagnosis (relayed 2026-11-09): the QE vs BAANDI energy offset of the Si.pz-vbc cases (+1.273e-6 Ha per atom, si8 +1.0184e-5 Ha) is one-electron (97.8 percent of the total; Ewald exact, Hartree +8.2e-7 Ry, xc -3.7e-7 Ry for si8), proportional to 1/volume (2-atom primitive cell at a = 4.90, 5.43, 6.00 Angstrom: dE/atom x volume constant to 0.78 percent) and a fixed relative error 1.7473e-5 of the G=0 local term (diff / local_zero_mode constant to 0.034 percent); cause: the log mesh starts at r_1 = 0.001308 and the reference G=0 constant omits [0, r_1] (about 4 r_1^2 = 6.8e-6 against the integral -0.3913: 1.75e-5), BAANDI integrates the Coulomb part analytically; new option `local_mesh = "qe_like_g0"` (default off) | adopted as a research option; Mac: si8 strict total energy -31.32181378329 Ha against QE -31.321813785 (difference 1.7e-9 Ha, QE print resolution) instead of 1.0184e-5; dev08 confirmation of si16/si32 and the lattice scan requested | The earlier 'per atom constant' reading was an artefact of three cells with one density; the radial keys could not move it because it sits in the G=0 constant, not in the quadrature; the speed-table ratios can now be given at equal target (QE-like convention) with the option on |
| QE-START | dev08 QE start-option benchmark (relayed 2026-11-09, QE 7.6, 1 core, same inputs): QE logs show `atomic+random` start wavefunctions (more atomic orbitals than bands: 32/64/128 for si8/16/32) and `superposition of free atoms` starting potential (charge 31.9960 -> renormalised 32.0000); startingwfc: default atomic+random vs atomic vs random: si32 1e-10 default 13 SCF / 186 h_psi / 6.623 s, atomic 8 / 35 / 3.550 s (-46 percent), random 15 / 219 / 8.028 s; si32 1e-8 default 4.981 s, atomic 3.022 s; si16 1e-10 1.826 -> 1.322 s; si8 about -4 percent; mixing beta 0.5 ndim 6 changes si32 1e-10 by -6 percent only (so the start, not the mixing); QE sensitivity: SCF counts invariant at 1e-7 and 1e-3 Bohr perturbations (inner work +15 to 20 percent at 1e-3, wall +8 to 10 percent); QE per-iteration ethr spans eleven orders of magnitude (si32 1e-10: iterations 1 to 3 use 2.0, 1.0, 2.0 average Davidson iterations at ethr 7.6e-4); QE on the hard cells of benchmarks/cases (fixed occupations, no smearing, nosym, 1 core): si-skew 51 / 180 SCF iterations (1e-8 / 1e-10), si-low 11 / 15, si-elongated 137 / 200 (1e-10 not converged), si-slab 200 / 200 (neither converged); BAANDI SCF counts equal QE's on si32 strict (13 vs 13); BAANDI start density test (Mac, opt-in `start_density = "atomic"`): si8 14 -> 15 SCF, si16 14 -> 14, wall slightly worse, i.e. the starting DENSITY does not move the SCF count | the QE baseline for 'best configuration' is startingwfc = atomic; BAANDI's atomic start wavefunctions (P46A) are now the relevant comparison; start density superposition stays opt-in, no gain | QE's default start is not its best (randomisation harmful on large cells); QE's low SCF count does not survive hard geometries (180 iterations on a sheared cell, no convergence on the slab and, at 1e-10, the elongated cell) |
| ATOMIC-DEFAULT | P46A atomic start wavefunctions (`start_wavefunctions = "atomic"`) made the default (ADR 0042) after the dev08 pre-registered rule R3 passed 8 of 8 (relayed 2026-11-09, f75a817, MKL and pure, fresh seeds 801-810 / 801-805, T=1 and T=8): perturbed median wall gains -15 to -42 percent, unperturbed si32 strict 40.640 -> 19.889 s (13 -> 10 SCF), si16 strict 7.485 -> 4.014 s, si8 standard 0.998 -> 0.624 s; atomic more SCF on 0 percent of seeds; all converged and certified; |dE| at most 1.65e-10 Ha (standard) and 3.4e-13 (strict); SCF spread over perturbed seeds 7 -> 0, 2 -> 0, 2 -> 1, 1 -> 0; SHA stable and T=1 equal T=8 | adopted as default (user-requested equal-footing evaluation against QE's best start) | Mac: si8 strict real, default now 10 SCF iterations (14 before), same energy | The start was the largest single lever left on the SCF side; the QE gap at one core falls from median 4.0x to about 2.9x and si32 strict from 11.6x to about 5.7x |
| MPI-M1B | Orchestrator: row-block distributed dense algebra primitives (`dft_solvers::distributed`: `RowBlock` partition with scatter/gather, `gram` = local A^T B GEMM plus deterministic tree sum, `rotate_local`, `project_out` = X - Q (Q^T X), `column_norms_squared`), generic over `dft_mpi::Collective`, no MPI needed to test | adopted (M1b of ADR 0041, plan docs/roadmap/mpi-m1-plan.md) | tests on 1, 2, 3, 4 and 8 simulated ranks (threads): gram within 8 n eps relative of the serial GEMM, bit-identical on every rank and run to run at fixed rank count, bit-identical to the serial GEMM at one rank; projection, rotation and norms equal serial within 1e-12; Mac observation only (correctness, no scaling claim) | The orthonormalisation of the correction block cannot reuse the scalar MGS under distribution (every dot would be a reduction: about 1000 small reductions per Davidson step); M1c uses a panel Gram plus Cholesky (CholeskyQR2) with a rank floor and a fallback |
| STR3 | Orchestrator: the three cutoff-scan SCF runs of the stress pass (centre, +5 and +10 percent cutoff) run concurrently in scoped threads with ceil(T/3) threads each when the thread budget is at least 3 (sequential otherwise); each run is an independent deterministic SCF whose result does not depend on its thread count, so every energy and the whole stress block are unchanged bit for bit | adopted | Trigger: dev08 derivative-pass measurement at the atomic default (relayed 2026-11-09): the Pulay assessment is 84 to 87 percent of the added cost of forces plus stress (si16 T=1 7.926 of the added 9.44 s; T=8 4.650 s), threads it only 1.70x, while forces (0.25 s) and stress proper (1.2 s) are small and serial (1.02x, 1.03x). Mac observation, si8 strict with stress, clean A/B of two binaries differing only by this change: T=1 wall 1.412 -> 1.436 s (unchanged), T=8 Pulay assessment 0.443 -> 0.229 s (-48 percent), wall 0.760 -> 0.552 s (-27 percent), stress block identical | dev08 confirmation requested; forces and stress proper stay serial (about 15 percent of the added cost) | Independent runs inside one call are the easy parallelism the threaded SCF cannot reach |
| HARD-CELLS | dev08 hard cells at the atomic default (relayed 2026-11-09, MKL, T=1, 1 run each, fixed occupations, nosym; QE at conv_thr 1e-8 / 1e-10): si-skew standard 16 SCF (free electron 16, QE 51), strict 13 (17, QE 180, wall 0.190 s against QE 0.696 s); si-low standard 8 (9, QE 11), strict 11 (11, QE 15); si-elongated standard now aborts with Diverged ('density residual exceeded growth limit') at iteration 123 where free electron exhausted its 200 iterations, strict exhausted (QE standard converged in 137, strict not); si-slab exhausted in both classes (QE not converged either); multi-k with the atomic default: k222s and k444s converge in 8 SCF (9 before), gamma_real_used false with reason 'multiple k points or shifted mesh'. Mac probe on the two cells (standard, 200 iterations, atomic start): mixing variants (beta 0.5 / 0.7, Hartree metric, Kerker-Pulay q0 = 1, history 8) all end not converged, final density residuals 1.6e-3 to 1.1e-2 (elongated) and 1.2e-3 to 8.3e-3 (slab) | recorded; no mixing setting rescues the two cells; si-skew and si-low are where BAANDI beats QE (3.7x faster, 13.8x fewer SCF on si-skew strict); si-elongated standard is the one cell where QE converges and BAANDI does not (open) | Root cause of the elongated/slab failure is not identified (charge sloshing or near-degenerate fixed occupation are the candidates; QE's `nbnd = N/2` input differs from BAANDI's extra bands) |
| MPI-M1CDE | Orchestrator: distributed real Davidson (`RealDavidson::solve_distributed`, CholeskyQR2 panel orthonormalisation with pivot floor, replicated reduced eigenproblem, H on owned column pairs after allgatherv), engine branch for two or more ranks (`dft_engine::parallel` thread-local communicator, one rank keeps the existing path), CLI `--features mpi` (MPI_Init in main, non-root ranks write `<out>/rank<r>/`, `parallel` key in results only when ranks>1) | adopted pending the dev08 scaling gate (ADR 0041) | simulated ranks P=2,3,4,5: eigenvalues within 1e-9 of serial, residuals < 1e-8 by independent operator application, all ranks bit-equal, bit-reproducible at fixed P; engine SCF on 2 and 3 simulated ranks within 5e-9 Ha of threads-only on two cell shapes, distributed path confirmed by counter; Mac correctness only, no scaling claim | Not yet measured with real MPI: needs dft5090 to build `--features mpi` on dev08 under /scr and run the pre-registered gate (si32 strict >= 2.2x at 4 ranks, >= 3.5x at 8 ranks vs 1 rank). Non-Davidson SCF stages are replicated per rank in M1 (no scaling there), so the gate measures the Davidson share only |
| MIX-GROWTH-RESET | Orchestrator: opt-in `BAANDI_MIX_GROWTH_RESET=<factor>` resets the Pulay history when the density residual grows by more than the factor in one iteration (default off) | rejected as a fix, kept as opt-in (inert) | dev08 (relayed 2026-11-09, master): inert at factor 3 on si-skew, si-low, si32, si16 (iterations and energies bit-identical, standard and strict); si-elongated standard, 50 iterations: max growth 15.6 to 4.2 and resets fire, but the best residual gets worse (free_electron 9.25e-04 to 2.94e-03; atomic 2.95e-03 to 2.87e-03 at 3, 4.0e-03 at 2.5), nothing converges at factor 2.5, 3 or 5 | The growth spikes (up to 16.9x, never above 20 so the factor-20 guard never fires) are a symptom of the 8e-04 plateau, not its cause: resetting removes the extrapolation that reached the plateau. The floor is identical across mixing variants, starts, band counts and grid (15x15x45 as QE); next hypothesis is a fixed density-potential inconsistency of about 3.6e-05 Ha against QE -15.138214 (term-by-term energy diff on this cell at QE's density). si-slab shows the same pattern (growth up to 9, residual floor 1.0e-03) |
| MPI-M1-GATE1 | dev08 first scaling measurement of the M1 distributed Davidson (relayed 2026-11-09, master; cli-mpi 54f4c28, MKL, Open MPI 5.0.8, si32 strict, np=1 with the same MPI binary as baseline) | gate not met; distribution correct | wall 20.409 s (1 rank), 17.133 s (4 ranks, 1.19x), 15.638 s (8 ranks, 1.31x = 16% efficiency); energies agree to 2.1e-13 Ha; SCF iterations 10 at all rank counts; thread path on the same 8 cores 2.43x (8.179 s); Davidson iterations 91 -> 121 (4 ranks) and 117 (8 ranks) while H columns stay flat (3397 -> 3432/3437); np=2 smoke: rank1 output bit-identical to rank 0 (density.bin SHA e6fef2f314e08103), np=2 vs np=1 dE 1.059e-12 Ha | Cause of the extra Davidson iterations found on the Mac (same session): the panel orthonormalisation by Gram-matrix Cholesky (CholeskyQR2) dropped nearly dependent correction columns that the serial in-block MGS keeps (n=700, 32 bands: 21 of 32 kept at the first iteration, 8 iterations against 6); replaced by a distributed MGS panel (two passes per column, drop only at relative norm <= 1e-12, extra re-projection against the basis); simulated 2 ranks now match the serial iteration counts (6 and 6, 207 against 206 H columns). Amdahl is the remaining limit: density, potential, mixing and nonlocal setup are replicated per rank in M1; to be re-measured on dev08 |
| MPI-M1-THREADBUG | Orchestrator: found from the ac8cc91 gate (relayed 2026-11-09, master): with `scf.threads >= 2` the SCF solves k points in pool threads (`ordered_solves`), where the thread-local communicator is not visible, so every rank silently ran the serial solver (np=4 x T=2 bit-identical to np=1, Davidson 91, dE exactly 0); the hybrid rows of that table (np=2 x T=4 1.83x, np=4 x T=2 1.24x) therefore measured replicated thread-path runs, not distributed ones. Fix: with an active communicator the k points are solved sequentially on the caller thread (also required by MPI thread level single); test `distributed_with_threads_takes_the_distributed_path_and_keeps_iterations` (3 simulated ranks, threads 2, 24 bands: distributed path counter increments, Davidson iterations 38 against 38 serial) | adopted; hybrid and pure-MPI rows must be re-measured | Open: pure MPI at threads=1 still showed 121 Davidson iterations on si32 at ac8cc91 although simulated ranks match serial on synthetic Hamiltonians, so a second cause specific to the real Hamiltonian or the atomic start is not excluded | Lesson: every scaling row needs a proof the distributed branch ran (counter or per-rank collective count in results.json), not only a rank count |
| MPI-M1-GATE2 | dev08 re-run of the M1 scaling table at a12ad8f (relayed 2026-11-09, master; cli-mpi c9a315819c2c7ecf, MKL, si32 strict, 8 physical cores of socket 0, median of 2) with proof the distributed branch ran (ranks>1, rank dirs, dE vs np=1 nonzero at 1e-13) | gate not met; MPI slower than threads on one node | wall: np=1 20.234 s; np=4 19.107 (1.06x); np=8 16.723 (1.21x); np=2xT=4 15.422 (1.31x); np=4xT=2 14.957 (1.35x, best); thread path np=1xT=8 8.131 s (2.49x); Davidson iterations 91 at np=1 and 117 to 121 distributed, H columns 3397 to 3439, dE vs np=1 up to 7.7e-13 Ha | The extra Davidson iterations are all in SCF step 10 (33 at np=1, 59 to 63 distributed; steps 1-9 and 11 identical), a single strict band (index 64, first empty band of 68, one member of a degenerate cluster cut by the strict prefix) plateaus at about 1.03e-6. Mac reproduction (ledger-only, `clustered_*` ignored tests, n=800, 68 bands, strict prefix 65): serial and 4 simulated ranks give equal iteration counts on 8 of 9 clustered cases (18/18, 31/31, 26/26, 26/27, 35/35, 15/15, 17/17, 27/27) and one outlier 16 vs 24, so the distributed algorithm matches the serial distribution and the step-10 excess looks like rounding-level sensitivity of a cluster-cut tail, not a defect. To be confirmed on dev08 with the perturbed-start protocol (distribution of total Davidson iterations at np=1 against np=4). Per-collective call counts and seconds are now in results.json `parallel.collectives` for the next measurement |
| MPI-M1-GATE3 | dev08 perturbed-start protocol and collective accounting at 337f134 (relayed 2026-11-09, master; cli-mpi 39c456ac8ba9b614, si32 strict, seeds 601-610, +-1e-7 Bohr, T=1) | adopted as the standard scaling case; Davidson excess is sensitivity, not a defect | total Davidson np=1 and np=4: median 283, min 221, max 308 (identical distributions, equal on 9 of 10 seeds, seed 606 +1); np=4 wall median 35.315 s against 41.598 s at np=1 (1.18x) with identical work; unperturbed si32 is anomalously easy (91 Davidson against perturbed median 283); np=4 x T=2 unperturbed wall 15.269 s: 261 allgatherv and 19092 tree_sum calls per rank, collectives 7.6 to 13.0 percent of wall (1.16 to 1.98 s per rank, 0.82 s spread = load imbalance, bottleneck rank differs per collective) | Communication is not the limit: at least 87 percent of the wall is compute that does not scale, so the gap to the thread path (8.1 s) is replicated or non-distributed work. Next measurement: per-phase timers per rank at np=1, np=4 T=1 and np=4 T=2 (diagonalization should be about a quarter if the distributed Davidson scales) |
| MPI-M1-GATE4 | dev08 per-phase timers per rank at 337f134 (relayed 2026-11-09, master; si32 strict, perturbed seed 601, Davidson 308 in every configuration): np=1 T=1 total 41.3 s (diagonalization 33.9, potential 1.69, density 1.47, total_energy 4.09, mixing 0.12); np=4 T=1 total 35.2 s (diagonalization 26.5, potential 2.0, density 1.6, total_energy 4.9, mixing 0.13, ranks within 0.2 s); np=4 T=2 total 27.3 s (diagonalization 21.6); dgemm sites (dgemm-timers build, wall identical to 0.4 percent) | finding: the distributed Davidson does not scale; the replicated stages are the smaller term | diagonalization falls only 1.28x at 4 ranks (ideal 4x); dgemm is 12.8 percent of it (3.4 s of 26.5) and collectives 1.28 s, leaving 21.8 s non-BLAS non-communication against 24.1 s serial; per-rank GEMM flops 173 of 414 GFLOP serial (41.8 percent, four ranks execute 167 percent of the serial work); the replicated stages are 8.6 s = 24 percent of the 35.5 s wall (cap 4.8x); rank 0 is consistently last to arrive (allgatherv wait 0.38 s against 3.0 to 3.4 s on the others) | Mac findings, same session: (1) the operator columns are divided among ranks (test `operator_columns_are_divided_among_ranks`: 116 apply calls serial, 116 in total over 4 simulated ranks), so H-psi itself is not replicated by construction; (2) the distributed solver recomputed the full reduced matrix every iteration (m x m over local rows) where the serial solver forms only the new columns, which explains the 67 percent extra GEMM work; now incremental. (3) per-stage wall timers (reduced_gram, reduced_eigen, rotate, residual, orthonormalize, apply_total, apply_operator) are in results.json `parallel.stages` to attribute the 21.8 s; the replicated dsyevd of the reduced matrix (m up to 408, every iteration, every rank) is a prime suspect for the non-scaling part |
| MPI-M1-GATE5 | dev08 stage attribution at d149999 (relayed 2026-11-10, master; cli-mpi e2b895b3829cafe2, si32 strict, perturbed seed 601, T=1, 308 Davidson and 4909 H applications in every configuration) | M1 distributed Davidson limited by column-pair granularity; MPI stays slower than threads on one node | wall 41.968 s (np=1), 33.670 (np=4, 1.25x), 27.250 (np=8, 1.54x); incremental reduced matrix bought 6 percent at np=4 and cut per-rank dgemm GFLOP at np=4 from 173 to 30.4; rank-0 stages np=4 / np=8: reduced_gram 0.81/0.60, reduced_eigen 1.91/2.29 (the only stage that gets worse with P), rotate 5.75/1.38, residual 1.20/0.47, orthonormalize 4.98/2.60, apply_total 6.78/6.23 (operator 6.0/5.2); apply_operator on rank 0 is 1.95x (np=4) and 2.26x (np=8) the mean of the other ranks; replicated stages (potential, density, total_energy, mixing) 7.43 s serial, 8.72 at np=4, 10.28 at np=8 (38 percent of the np=8 wall); per-rank dgemm GFLOP doubles from np=4 (30.4) to np=8 (60.8), unexplained | Cause of the rank-0 asymmetry (code reading, `block_counts`): the H application is divided over column PAIRS and the first `pairs % P` ranks take the remainder, so a block of w columns gives ceil(w/2/P) pairs to rank 0 first; with the few unconverged columns of late Davidson iterations (4909 applications over 308 iterations, about 16 per call, many calls with fewer pairs than ranks) rank 0 or a few ranks do all of it and the others wait. Rotating the remainder does not help (the wall of a call is the maximum over ranks); a pair cannot be split without splitting the FFT itself, so the next real step is a G-vector distributed 3D FFT with transposes (plan M2), which also distributes the replicated stages | Decision for the user: stop MPI at M1 (opt-in feature, correct, 1.54x at np=8 against 2.49x for threads on the same 8 cores; dev08 has no multi-node evidence) and return to single-node work, or invest in M2 |
| MPI-M1-CEILING | dev08 trace of apply-call widths (relayed 2026-11-10, master; instrumented /scr clone, si32 strict perturbed seed 601, 458 Davidson iterations, 5232 correction columns) and a reporting finding | decides M1 against M2; reporting fixed | columns per apply call: mean 11.4, median 1; 66.4 percent of calls carry exactly one column, 11.8 percent carry 65 to 68 columns and 69 percent of the column work; with cost per call = ceil(pairs/P) the ceiling of a perfectly balanced pair split is 1.76x at P=2, 2.79x at P=4, 3.91x at P=8, 5.61x at P=32 (75 to 88 percent of calls have fewer pairs than ranks), so only splitting within a column (distributed FFT, M2) removes it. Reporting: 150 of 458 traced Davidson iterations (33 percent) were not in `history[].davidson_iterations`: they are the fixed-point certificate re-solves (54 and 78 iterations at two failed certificates and 18 at the final passing one), stored in each history entry's `fixed_point_certificate.davidson_iterations`, which every Davidson total taken from the history missed | adopted: `convergence.work` in results.json now reports recorded Davidson iterations and H vectors, certificate iterations, H vectors and solves, so the total work is the sum | The MPI gate comparisons stay valid (both sides measured alike) but absolute Davidson counts quoted from the history understate work by about a third in the strict class; the hard late SCF steps are exactly where certificates fail and re-solve |
| NONLOCAL-DGEMM | Orchestrator: real-path nonlocal projector products (P^T X and P W) through dgemm panels on MKL builds, opt-in `BAANDI_MKL_NONLOCAL=on` (default off, bit-identical when off; non-MKL builds forward to the existing kernel). The earlier note that the real nonlocal is a scalar loop was outdated: `RealNonlocal::add_block` already used the pure-Rust tiled GEMM, not MKL; the new `multiply_slices_wide_into` also admits the long-output non-transposed product (per-panel scratch bound only) | implemented, awaiting dev08 A/B | Mac: clippy and tests clean, `cargo check --features dft-numerics/mkl` compiles; no speed claim | Needs dev08: same-binary env A/B (off against on) on si32 strict and standard, perturbed seeds 601-610, T=1 and T=8; accept only if energies stay in the accuracy classes, density SHA stable per setting, and the nonlocal share (stage decomposition) drops; the coupling loop between the two products stays scalar |
| NONLOCAL-DGEMM-AB | dev08 same-binary env A/B of BAANDI_MKL_NONLOCAL (relayed 2026-11-10, master; cli-mkl-dc9dc03 d07207338af24876, si32 perturbed seeds 601-610, interleaved off/on per seed, T=1 taskset core 0, T=8 cores 0-7) | NOT ADOPTED as default by the pre-registered rule (>= 3 percent median wall in both T=1 and T=8 on both profiles); stays opt-in | median wall off to on: strict T=1 41.329 to 38.691 s (-6.38 percent, worst seed -5.67), strict T=8 16.609 to 16.270 (-2.04, worst -1.57; fails), standard T=1 13.615 to 13.039 (-4.23, worst -3.24), standard T=8 6.043 to 5.955 (-1.46, worst -0.38; fails); all 40 seed pairs faster; energy differences 1.4e-14 to 7.7e-13 Ha (inside both classes); total Davidson identical on 38 of 40 pairs (seed 605: 567 to 569, seed 606: 491 to 493); density SHA reproducible per setting and equal at T=1 and T=8; dgemm sites: only `other` moves (72474 to 116090 calls, 13.2 to 53.7 GFLOP; the 40.5 GFLOP moved cost 1.48 s = 27 GFLOP/s, the pure-Rust tiled GEMM did the same work at about 9.7 GFLOP/s) | Lesson: a single-threshold pre-registration hid a clean 4 to 6 percent T=1 win with zero regressions; the T=8 gain shrinks plausibly because the Rust tiled path is already spread over 8 pool threads (inference, not measured; a T=2/T=4 sweep would test it). Decision rule kept as registered: no default change. Future pre-registrations should state the thread budgets the default is meant for (a single-thread default is a legitimate scope for a numerical-kernel switch). Proposal for the user: adopt as default only for MKL builds at T=1 if the scope is declared in a new pre-registration, or keep opt-in |
| ELONGATED-GAP | dev08 trace of the non-converging si-elongated cell (relayed 2026-11-10, master; BAANDI_EIGTRACE in a /scr clone only; 16 electrons, 8 occupied of 12 bands, cell 8x8x24 Bohr, standard profile, 400 iterations) | diagnosed: ill-conditioned fixed-occupation problem, not a mixing or solver defect | HOMO-LUMO gap median 1.98e-4 Ha (min 6.1e-7, max 2.2e-2) over SCF iterations 20-400 with free-electron start, below 1e-3 Ha on 79 percent and below 1e-2 Ha on 96.6 percent of iterations, final 4.1e-5 Ha; it jumps by three orders between consecutive iterations (9.7e-5 at 31, 2.5e-3 at 33, 7.5e-3 at 46, 1.1e-5 at 47); bands 7, 8, 9 are a near-degenerate triplet (17 meV) and the 8-occupied cut falls inside it; atomic start median gap 7.0e-3 Ha and diverges on the growth limit instead; no level swap or crossing over 380 iterations; QE ran nbnd = 8 (no empty state), mixing beta 0.7, converged in 137 iterations but its estimated scf accuracy oscillates between 0.02 and 1.06 Ry from iteration 7 to 129 (limit cycle, escapes at 130), highest occupied level differs by 29 meV from BAANDI | The eigen solve exits 110 to 230 times (max 650 to 1030) above the requested eigen tolerance with zero unconverged bands: the standard profile converges bands on Ritz-value change (`ritz`), by design uncertified, which is too loose to fix the orientation inside a near-degenerate cut and feeds noise into the density (density floor about 380x the requested eigen tolerance); the strict class also fails on this cell (floor 8.75e-4). Remedy is physics-side: occupations.mode = fermi_dirac with a smearing width above the triplet splitting (the engine supports it), not more mixing variants; next: run the cell with smearing 0.02, 0.05, 0.1 eV and with extra_bands 0 (QE-like eigenproblem). Proposal: a not-converged hint when the gap falls below 1e-3 Ha at fixed occupations (needs per-iteration gap in the history records) |
| ELONGATED-SMEARING | dev08 remedy test on the si-elongated cell (relayed 2026-11-10, master; cli-mkl-dc9dc03, strict profile because the standard profile is rejected under Fermi-Dirac, free-electron start) and extra_bands 0 control | smearing fixes the cell; extra_bands 0 refuted as a remedy; si-slab is the same mechanism | fermi_dirac strict: 0.02 eV converged in 27 iterations (E -15.138036632 Ha, bands 7/8/9 occupations 1.99964/1.60271/0.25583), 0.05 eV in 22 (E -15.137042296), 0.1 eV in 20 (E -15.134665449), band caps 16 and 32 agree to 1e-10 Ha; the smeared E approaches QE's fixed-occupation -15.138214450 Ha as the width shrinks (3.55e-3, 1.17e-3, 1.78e-4 Ha away); extra_bands 0 fixed occupations: standard diverges at iteration 40 (min density residual 7.2e-5, below the 8.3e-4 floor of the 12-band run, then lost), strict diverges too; si-slab (standard, 200 iterations): gap median 2.77e-4 Ha, 77.3 percent of iterations below 1e-3 Ha, ordering changes 0, eigen residual / requested tolerance median 111, same mechanism within a factor of 1.4 of si-elongated | Config findings: the standard profile is rejected under finite temperature with a message that named unpolarized SCF although the run is unpolarized (the condition is the occupation channel path); message fixed to say closed-shell fixed occupations are required and to suggest strict with smearing; band_cap must be at least 13 for 8 occupied bands (empty-state coverage rule). Remedy for near-degenerate cuts at the Fermi level is smearing, not band count or mixing | Open: repeat the sweep on si-slab; a not-converged hint when the gap is below 1e-3 Ha at fixed occupations (needs the gap in the history record); whether standard profile support under smearing is wanted |
| SLAB-SMEARING | dev08 Fermi-Dirac strict sweep on si-slab (relayed 2026-11-10, master; cli-mkl-77a41d7 422757fa663b8eec, band_cap 16, free-electron start, 400-iteration bound) and the QE reference run | remedy transfers; si-slab is a cell QE does not converge either | smearing 0.02 eV converged in 29 iterations (E -15.232295396 Ha, mu -0.07126582), 0.05 eV in 35 (E -15.231672422), 0.10 eV in 32 (E -15.229428001), certificate passed in all three, occupations sum to 16.000000000; fixed occupations floor at 1.045e-3 density residual for 200 iterations; bands 4/5 and 8/9 are EXACTLY degenerate pairs (identical to 8 digits at every width), band 7 sits 9.9e-4 Ha below mu and the empty pair 1.6e-3 Ha above it; occupation of band 7 at 0.10/0.05/0.02 eV: 1.297/1.436/1.587, of bands 8 and 9: 0.397/0.284/0.207 each (trend toward integer filling, not reached at 0.02 eV); QE (orig-mpi-dfti, nbnd 8) prints `convergence NOT achieved after 200 iterations` with a best estimated scf accuracy 2.26e-4 Ry at iteration 61 and 0.178 Ry at the end | Both hard cells now read as one mechanism: a fixed integer occupation cut inside a (near-)degenerate cluster at the Fermi level; QE needs 137 iterations (123 in a limit cycle) on si-elongated and fails on si-slab, BAANDI with smearing converges both in 20 to 35 iterations. The benchmark set must label these cells accordingly: fixed-occupation comparisons on them are not evidence of a BAANDI convergence defect; compare smeared runs, or QE with smearing. A zero-width extrapolation (widths down to 0.005 eV) would test whether the fixed integer occupation is a valid ground state of si-slab at all; not requested |
| GAP-HINT | Orchestrator: a not-converged fixed-occupation run now prints a hint when the final indirect HOMO/LUMO gap is below 1e-3 Ha (use Fermi-Dirac smearing, strict profile, band_cap above the occupied count); `docs/cli.md` gains the rule of thumb with the measured evidence | adopted (message only, no numerical change) | check.sh PASS; no history or checkpoint format change (the gap comes from the final eigenvalues already printed) | The per-iteration gap was deliberately not added to the iteration record: the record is part of the checkpoint wire format |
| SECOND-VECTOR | Orchestrator: opt-in `BAANDI_SECOND_VECTOR=v1|v2|v3` fills the free slot of an odd-sized correction block of the serial real Davidson (v1 previous correction of the selected band with the largest residual, v2 preconditioned residual of the converged band with the largest residual, v3 unpreconditioned residual of the worst strict band); default off, bit-identical; plan and acceptance rules in docs/roadmap/second-vector-plan.md | implemented, Mac signal unfavourable, dev08 A/B decides | Mac synthetic clustered fixture (n=800, 68 bands, strict prefix 65, 12 cases, serial iterations/H applications): total iterations over the 12 cases: off 280, v1 296, v2 285, v3 286; example g4s1 26 -> 34 (v1, v2, v3), g3s0 18 -> 17 (v2); H applications rise with every variant; the existing real_davidson tests pass with the variant off and v1, `buffer_and_ritz_semantics` (real against complex to 1e-12 on buffer bands) fails with v2 and v3 as expected for a changed trajectory | The fixture applies mostly wide blocks, so it does not represent the single-column tail of si32; the pre-registered dev08 A/B (si16/si32 strict and standard, seeds 601-610, three non-silicon fixtures, T=1 and T=8) is the decision; a variant that fails any of the five acceptance rules is rejected |
| SECOND-VECTOR-AB1 | dev08 same-binary env A/B of the first registration (relayed 2026-11-10, master; cli-mkl-2290b01 8fd6b1b6534fecb9, T=1 core 0, perturbed seeds 601-610) | v1, v2, v3 rejected as registered; conditional gain found on the premise cell; second registration v4 | si16 standard (first cell of rule 2): median Davidson off 43, v1 43 (0.0 percent), v2 42 (-2.3, worst seed +18.2), v3 52 (+20.9, worst +63.2); wall +0.75/+0.46/+6.63 percent; energies inside class (v3 up to 3.3e-9 Ha, 66 percent of the standard class); si32 strict T=1 (measured separately, not in the verdict): off 432 / 40.513 s, v1 416.5 (-3.6) / -0.5 percent wall, v2 402 (-6.9) but H +6.9 percent and +3.0 percent wall, v3 384 (-11.1), H -6.6 percent, wall 36.740 s (-9.3 percent), worst seed +8.3; all runs converged, energies <= 1.2e-12 Ha | The idea pays only where a single-column tail exists (si32 strict) and costs where it does not (si16 standard, Mac fixture), so v4 = v3 gated to blocks with exactly one column, with a measurable gate counter (`convergence.work.second_vector_columns`); rule 2 reworded to remove the per-seed against per-cell ambiguity (plan section 'Second registration'); T=8, density SHA reproducibility, certificates and the non-silicon fixtures were not measured in the first run |
| CERT-SCREEN | Perf owner (single-node sub-session, scope decision U10), opt-in `BAANDI_CERT_SCREEN=on` (strict class; default off, bit-identical when off): before the fixed-point certificate's tight re-solve, a first-order estimate of the occupied-eigenvalue movement it will measure (limit 1e-7 Ha) is computed from the step's stored per-band densities (single bands: diagonal element of the potential change) and, for bands within 1e-3 Ha of each other, from the small degenerate-perturbation matrix of the cluster (one FFT per band of the cluster); above 1.02 x 1e-7 Ha the certificate is declared failed without the solve. A failed certificate's re-solve is discarded by the SCF loop anyway, so the trajectory cannot change; a pass still needs the full solve. `BAANDI_CERT_SCREEN=trace` prints estimate and actual. Code: `ScfProblem::certificate_shift_estimate` (scf.rs), session.rs certificate branch, `BandDensity::weighted_integrals` | implemented; Mac observation strongly positive, dev08 A/B pre-registered (docs/performance/requests/single-node-speed.md) | Measurements (Mac, pure Rust, observations): (1) the estimate equals the real re-solve change to 4 digits on 52 certificates (22 on the small set: si8/si16 seeds 601-605 incl. two failed ones, and 30 on si32 seeds 601-610 incl. 20 failed ones; skew, k222, k222s, low, high, synthetic s2/p3/d3/spd/cubic-degenerate; worst relative error 6e-4, absolute 5e-11 Ha); a plain cluster-mean version without the cluster matrix underestimated by 2 to 3x on cubic Gamma cells and was discarded. (2) Every one of the ten si32 seeds 601-610 fails its first TWO certificates on the eigenvalue criterion (changes 4.1e-7 to 9.6e-7 Ha against 1e-7; the density criterion passes at 2.5e-8 to 6.4e-8 against 1e-7) and passes the third: failed re-solves cost 52 to 80 Davidson iterations each (107 to 155 per run), the passing one 14 to 20; certificate work is 26 percent of all H vectors (1643 to 1795 of about 6600) and 30 percent of Davidson iterations. si16: one seed of five fails (1.15e-7), si8: one of five (1.30e-7); the remaining si8/si16 passes sit at 0.5 to 0.8 of the limit. (3) A/B, same binary, one process at a time, T=1: si32 seed 601 wall 30.40 -> 23.97 s (-21.2 percent), seed 602 27.91 -> 23.18 (-16.9); Davidson 458 -> 326 and 412 -> 290 iterations, H vectors 6660 -> 5352 (-19.6 percent) and 6563 -> 5273; T=8 seed 601 11.10 -> 8.90 s (-19.8); si16 seed 601 (one failed certificate) Davidson 227 -> 165 (-27 percent), H vectors -11.9 percent; seeds without a failure unchanged in work (si16-602 wall +1.2 percent, noise, screen cost about 0.02 s). Density SHA-256 and total energy bit-identical off against on in every compared run (17 of the 20-case set, base binary against the screen build; si32 601/602 at T=1 in three modes; si32 601 at T=8; si16 601/602 at T=8) | The strict class pays about a quarter of its H work for certificate re-solves that cannot pass because the SCF first reaches the density criterion several iterations before the eigenvalue criterion (Hartree long-wavelength amplification: eigenvalue movement is 10 to 20 times the density residual on si32). The first-order estimate is exact enough to predict the verdict with margin 2 percent and costs one potential build plus a few FFTs. A failed certificate carries zero solver work in results.json (`davidson_iterations: 0`, step residuals, estimate as `occupied_eigenvalue_change_ha`); no checkpoint format change. Not yet covered: the metallic/Fermi-Dirac and spin drivers (their own certificate code), standard class (no certificate). Needs dev08 per the pre-registered rules (identity bit-exact, no extra work, median wall -3 percent at T=1 and -2 percent at T=8 on si32 strict) |
| NL-ENERGY | Perf owner, opt-in `BAANDI_NL_ENERGY=fast` (default off, bit-identical when off): the nonlocal energy of the energy audit (`total_energy` phase) from the projector products alone (no back-application, no output vector), occupied bands only (unoccupied bands contribute exactly zero to the occupation-weighted sum), four bands per pass over a projector column with separate accumulators per band, so a band's value does not depend on grouping or thread count. Code: `NonlocalOperator::expectation_values` (hamiltonian.rs), energy.rs `band_nonlocal_partials` | implemented; Mac energy-phase -55 to -60 percent, wall effect inside noise, dev08 A/B pre-registered | Profile first (Mac sample, si16, about 3000 samples, directional): the complex `NonlocalOperator::add` from the energy audit is 107 self samples (3.5 percent of the run) and `band_nonlocal` sits under 97 of the 112 audit samples, so about 85 percent of the `total_energy` phase on si16 was the per-band legacy apply plus dot product (a latency-bound complex reduction followed by a second pass). Result (Mac, same binary): `total_energy` phase si16 0.210 -> 0.093 s, si32 T=1 1.64 -> 0.66 s and 1.62 -> 0.69 s, T=8 0.35 -> 0.21 s; total-energy difference 1.3e-13 (si16), 2.6e-13 and 2.8e-13 Ha (si32), eigenvalues bit-identical, density SHA unchanged, Davidson totals unchanged; unit test `expectation_blocks_match_apply_and_ignore_grouping` (equals apply plus dot to 1e-12 relative, grouping never changes bits) | On the Mac the saving is 1 s of 30 s at si32 (3 percent; wall A/B deltas -5.0 and +2.0 percent are inside the 1.5 s run-to-run spread of the Davidson phase); on dev08 the phase is 4.09 s of 42 s (MPI-M1-GATE4) with the legacy nonlocal as scalar Rust, so the expected saving is 2 to 3.5 s (5 to 8 percent), an estimate until the pre-registered A/B (phase timer -40 percent, wall -2 percent at T=1 and T=8, energy within 1e-12 Ha, work unchanged) |
| DENSITY-PAIR | Perf owner idea, not implemented: density accumulation for the real Gamma path with two bands per complex FFT (`PairDensityOperator::accumulate_density_pair` exists in dft-numerics and is tested; eigenvectors from the real path are exactly Hermitian-symmetric so each band is a real function) | untested | Mac profile: density phase 0.156 s of 3.5 s on si16 (4.4 percent), dev08 si32 1.47 s of 42 s (3.5 percent); halving the transforms would save at most 1.7 percent of wall. Not started because (a) it changes the per-band contribution arrays that CERT-SCREEN reads (a pair array merges two bands; the screen would need its own FFTs, about the cost it saves), (b) worker ranges must be rounded to pairs to keep the reduction order independent of the thread count, (c) rounding change | Record only; reopen if the density phase grows relative to Davidson (large k meshes) or the screen is dropped |
| POTENTIAL-REUSE | Perf owner note, not implemented: the screen builds the Hartree/XC potential of the output density and, when it predicts a pass, the certificate solve builds it again | untested | cost about one `potential` phase step (0.1 s on si32 Mac, 0.4 percent of wall) once per run | below the noise of any A/B; reopen only if the screen is made to run more often |
| SINGLE-NODE-AB1 | dev08 A/B of the single-node switches (relayed 2026-11-10, master; f1f50f9, binary 5065e753055774f8ea4e30ea4146882687e39f2cca514ec3c896a691eb07b19f, interleaved off/on, perturbed seeds 601-610, T=1 core 0 and T=8 cores 0-7, per-run records at /scr/mino/rustdft/sn/runs.jsonl) | BAANDI_CERT_SCREEN and BAANDI_NL_ENERGY=fast each REJECTED as registered by one clause (margins under one point); stay opt-in; a new registration with fresh seeds is proposed to the user | CERT_SCREEN: identity PASS (energy and density SHA equal in 40 of 40 pairs), no extra work PASS; wall si32 strict T=1 40.526 to 32.760 s (-19.16 percent), T=8 16.121 to 13.587 (-15.72), worst seed -15.24; si16 strict T=1 +0.78 percent (limit +1, pass), T=8 +1.36 (limit +1, FAIL; five repeats +1.49/+0.81/+1.10/+1.77/+3.09, mean +1.65: real small overhead); screens: si32 10 of 10 seeds had one (1305 Davidson iterations avoided per thread budget), si16 1 of 10 (46 avoided); the estimate is essentially exact (si32 seed 601: bound 6.5096e-7 vs actual 6.5095e-7 at iteration 10, 5.6136e-7 both at 11, 3.3184e-9 both at 12). NL_ENERGY=fast: energy within 1e-12 Ha in all eight cells (max 6.3e-13), Davidson and SCF counts unchanged in 80 of 80 pairs, wall PASS (si32 strict T=1 -5.67, T=8 -2.46, standard T=1 -12.84, T=8 -4.21, worst seed +0.88), total_energy phase si32 strict T=1 -57.23 (pass), T=8 -39.45 against -40 required (FAIL; three repeats -39.46/-39.28/-39.54, firm). Both together on si32 strict: wall T=1 40.607 to 30.802 s (-24.15 percent), T=8 16.143 to 13.286 (-17.70), total_energy phase -62.9/-48.3 percent, energies within 3.7e-13 Ha, density SHA equal in 20 of 20 pairs | Verdicts stand as registered (rules are not moved after data). Observation for a new registration: the failed clauses are a cost bound on a cell where the screen has nothing to screen (about +1.6 percent of a 2.9 s run) and a timer proxy clause (phase -40 percent at T=8) whose real target, wall, passes with margin; a new registration on fresh seeds (701-710) with the same large-cell gains and a small-cell cost bound stated up front is a user decision because it changes defaults. Rule 5 of the second-vector plan is dropped (restart counts not reported) |
| POTENTIAL-REUSE-IMPL | Implementer (proposal B3), with `BAANDI_CERT_SCREEN=on|trace`: the Hartree/XC fields the screen builds for the output density are kept and consumed by the next map only on a bit-exact density match (the certificate solve), dropped otherwise and after every certificate branch; `BAANDI_POTENTIAL_REUSE=off` disables it for A/B; nothing kept with the screen off. Same change: the screen's cluster band transforms run one per pool worker with their own FFT plan (`ordered_chunks`) and the cluster matrix entries are split over the pool (`cluster_entries`), each entry a serial sum, so the estimate is bit-identical for any thread count; trace lines gain `screen_s=`. Code: scf.rs `ReusedFields`, `certificate_shift_estimate(.., keep_fields)`, session.rs `potential_reuse()` | implemented, opt-in through the screen; part of the second registration (docs/performance/requests/default-candidates.md) | Mac observations (pure Rust, shared machine, wall noise +-10 percent per run): instrumented screen on si16 T=8 seeds 601-602: total 30 ms = potential build 2.5 ms + band transforms 20.5 ms (31 bands in clusters of 2 to 6) + cluster entries 4.8 ms; after the change screen_s 15.2 to 17.0 ms on seeds 601-605 (was 33.7 to 35.5 ms), certificate `potential_s` 2.5 ms -> 0.08 ms (si16) and 5.1 ms -> 0.14 ms (si32 T=8); bounds identical to the printed digits; density SHA and energy bit-identical off / screen-no-reuse / screen in all 45 si16 T=8 runs (seeds 601-605, three repeats) and in the si32 T=8 runs (seeds 601-603); si16 T=8 median wall vs off: screen without reuse +2.65 percent, with reuse +0.86 percent (observation, inside the noise); si32 T=8 screen -21.0 percent median (Davidson 1284 -> 873 over three seeds). Unit test `certificate_potential_reuse_is_bit_identical` (cubic cell with a threefold cluster: estimate, eigenvalues, density and energy bits equal for reuse on/off and 1/4 threads) | The proposal's premise (overhead = one potential build) was wrong on the Mac: the build is under a tenth of the screen; the per-band transforms were the cost and they parallelise. Whether the +1.4 to +1.7 percent si16 T=8 overhead of SINGLE-NODE-AB1 disappears is for dev08 under the new registration; at T=1 the transforms stay serial (the AB1 T=1 cost passed its bound) |
| ENERGY-IDENTITY | Implementer (proposal B2, parity-plan row 2.3 / Q6), opt-in `BAANDI_ENERGY_IDENTITY=on` (default off, bit-identical when off): in ordinary SCF maps of the integer-occupation spin-unpolarized driver the nonlocal energy is E_nl = sum f eps - E_kin - integral V_loc n_out - integral (V_H + V_xc)[n_in] n_out (exact for Rayleigh-quotient eigenvalues of H[V_in]), so no projector work; the history `energy` of those iterations equals the band-structure form (`band_form`). The certificate solve keeps the explicit audit (flag cleared before it); the final state's nonlocal energy is recomputed explicitly with the audit's arithmetic (`ScfProblem::explicit_nonlocal_energy`, honours `BAANDI_NL_ENERGY`) and the identity-minus-audit difference is reported as `convergence.work.energy_identity` (results.json, nested). Fermi-Dirac maps fall back automatically; the collinear and finite-temperature drivers have their own energy code and are not touched | implemented, opt-in; SEMANTIC CHANGE of history energies needs the user's decision before any default | Unit test `energy_identity_matches_explicit_audit` (skewed one-atom cell with s+p projectors): one map identity vs explicit nonlocal differs by 3.1e-15 Ha (all other audit fields bit-identical); full SCF on/off: same iteration count, final energy and density bit-identical, final-state difference 2.7e-15 Ha. Mac (observations): si16 T=8 seeds 601-605 x2: final energy and density SHA bit-identical to off in all runs, `total_energy` phase median 55 -> 32 ms; final-state |difference| at most 8.8e-13 Ha; si32 T=8 seeds 601-603: phase 0.396 -> 0.151 s (identity alone), NL_ENERGY=fast 0.239 s, both 0.125 s, all three switches 0.086 s; final difference at most 1.5e-12 Ha; si32 T=1 seed 601: phase 1.739 s off, 0.776 NL_ENERGY=fast, 0.506 identity, 0.232 both (wall of that series not usable: shared-machine load moved the Davidson phase by +-20 percent) | With NL_ENERGY=fast already taking most of the phase, the identity adds about 0.54 s at T=1 and 0.11 s at T=8 on si32 (Mac), i.e. 1 to 2 percent of wall, not the 8 to 10 percent the proposal expected against the old audit; dev08 phase share decides. The identity difference at the final state (1e-12 Ha class) is a consistency check of eigenpairs, Hamiltonian and density, worth reporting regardless of speed |
| B6-PHYSICS | Implementer (proposal B6): `RealNonlocal::{projector_count, last_projections, add_block_with_projections}` and `RealHamiltonian::{projector_count, last_projections, apply_block_with_projections}` in crates/dft-physics/src/real_gamma.rs; the default apply path keeps its arithmetic (coupling and back product moved to a private helper). Solver half (cache P^T V per basis column, rotate with restarts, certificate from (P^T V) C) NOT implemented: needs more than a hook in real_davidson.rs (owned by another agent); described in docs/performance/requests/b6.md; no `BAANDI_B6` switch exists yet | physics half implemented and tested; solver half deferred (recommended together with Track A) | Unit test `cached_projections_reproduce_apply` (skewed two-species cell, 104 projector columns): H X from cached projections vs plain apply relative 1.1e-16 (rounding class); direct projections reproduce apply bit for bit. Expected saving (E): one projector product per certificate apply, about 0.3 to 0.6 percent of a si32 strict run | Small alone; worth it only where every apply is full width (Track A) |
| C1-LADDER | Implementer (proposal C1, design docs/roadmap/c1-cutoff-ladder.md): `scf::cutoff_ladder::run_ladder` (engine, default off) and `[scf] cutoff_ladder = "on"` / `BAANDI_CUTOFF_LADDER=on|off|<fraction>` (CLI run path only): coarse SCF at E1 = 0.5 E2 on the E1 grid, stopped at density residual 1e-4 or 30 iterations, exact prolongation (orbitals by Miller index, density by Fourier zero padding), unchanged strict session at E2 from the prolonged density and warm orbitals with an empty mixer; eligible only for one exact Gamma point on the real path with fixed uniform occupations, strict profile, no symmetry/strain reference; everything else falls back with a recorded reason | implemented, default off; dev08 A/B requested (docs/performance/requests/c1-cutoff-ladder.md) | Mac T=1 strict, one run each (observations): full-cutoff SCF iterations plain -> ladder si8 11 -> 10, seeds 601/602 19 -> 12 and 17 -> 11; si16 seeds 601-603 20/23/22 -> 12/16/14; si32 seed 601 24 -> 15; si-skew 21 -> 11; synthetic fcc/skewed/slab 20/12/18 -> 12/12/12; full-cutoff H vectors si32 11122 -> 6747 plus 3426 coarse (1388 full-grid equivalents at grid ratio 0.42); wall si8 0.97 -> 0.73 (seed 601), si16 5.06/5.44/5.35 -> 3.77/4.47/4.38, si32 39.7 -> 28.7 s, si-skew 0.194 -> 0.109; worse only on a 24 ms synthetic cell (0.024 -> 0.034, coarse setup); every run certified at the full cutoff, max |dE| 8e-13 Ha; fractions 0.6/0.7 not better than 0.5 on si16 | The coarse stage replaces the expensive early full-cutoff iterations (far from self-consistency) at about 35 to 42 percent of the per-column cost; the full-cutoff part still needs 10 to 16 iterations from a 1e-4 hand-off (the coarse/fine fixed points differ at that level), so the gain is in the removed early iterations, not in a shorter tail. Standard profile and k meshes are not eligible yet (real path only) |
| C3-EXTRAPOLATION | Implementer (proposal C3, design docs/roadmap/c3-ionic-extrapolation.md): `[relax] start_extrapolation = "none" | "first_order" | "second_order"` / `BAANDI_RELAX_EXTRAPOLATION` (ions mode, warm start only; default none, bit-identical): start n_at(R_new) + d with d = n_scf - n_at carried (first order) or extrapolated by the least-squares fit of the new displacement on the last two steps (second order); n_at is the engine's atomic-charge superposition (new accessor `ScfProblem::atomic_superposition_density`, built with `start_density = Atomic`, which only sets the cold start a warm transfer replaces); negative samples clipped with the charge restored; history of three d fields in restart.json (optional field, resume exact); per-frame `scf_h_vectors`, start kind, clipped fraction | implemented, default none; dev08 A/B requested (docs/performance/requests/c3-ionic-extrapolation.md) | Mac T=1, six BFGS steps, SCF iterations steps 1-6: si8 Gamma 30 Ry none 84, first_order 76, second_order 70 (H vectors 9543 / 8668 / 8129; wall 5.98 / 5.54 / 5.27 s; cold every step 98, 6.21 s); synthetic dimer slab 117 / 100 / 107; synthetic fcc 121 / 112 / 123; synthetic H dimer with a Gaussian charge not matched to its potential 76 / 87 / 82 (worse); step-1 energies (same positions) equal to 1.5e-13 Ha; along the path the differences (si8 max dE 5.3e-8 Ha second order) are below the cold-versus-warm difference of the existing driver (3.0e-7 Ha), final positions within 5e-6 Bohr | The atomic part carries most of an ionic step's density change only when PP_RHOATOM matches the pseudopotential; second order helps on smooth BFGS paths (si8) and can lose on toy fixtures with irregular steps (synthetic fcc). Orbital extrapolation (design stage 2) not implemented: the transfer contract carries no orbitals and the density start alone already gives the measured gain |
| C3-DESIGN | Proposal C3 / parity-plan 7.2: ionic-step density (and optional orbital) extrapolation for relaxation | design note only: docs/roadmap/c3-ionic-extrapolation.md | none | Not a modest change: needs a public atomic-superposition accessor, an orbital-carrying transfer entry and a restart format bump; acceptance test pre-drafted there |
| C1-DESIGN | Proposal C1: cutoff ladder (nested bases, prolongation, full-cutoff hand-off and certificate) | design note only: docs/roadmap/c1-cutoff-ladder.md | none | Experiment pre-drafted for after Track A stage 1 |
| DAVIDSON-STAGES | Perf sub-session (proposal step 0 item 1): env-gated stage timers in the serial real Davidson, `BAANDI_DAVIDSON_STAGES=on` (default off; clock reads and relaxed atomic adds only, no I/O or allocation; default path bit-identical): iteration wall time and counts by correction width (0 = returning, 1, 2-8, 9-64, 65+) and stage seconds (reduced, eigen, rotate, residual, orthonormalize, h_apply, sketch, precondition), exposed by `real_davidson::stage_seconds()` and nested in `results.json` under `convergence.work.davidson_stages` | adopted as an opt-in measurement tool; dev08 part A requested (docs/performance/requests/b1-b4-b5-davidson.md) | Mac observation (pure Rust, T=1, load 5-8, seed 601): si32 strict iterations by width 15/265/83/29/54, ms per iteration 75/32/53/160/189; shares reduced 6.3, eigen 5.6, rotate 39.6, residual 0.5, orthonormalize 26.7, H apply 24.0 percent; at T=8 eigen 14.3 and orthonormalize 38.1 percent; si32 standard H apply 48.9, rotate 23.2, eigen 2.3; si16 strict eigen 8.6, rotate 31.8. Same energy (17 digits) and density SHA as the parent-commit binary with timers on and off on 8 cells | A width-1 iteration costs about 20 H columns on si32 strict (32 ms vs 1.5 ms per column): the per-iteration fixed cost, not H, dominates the tail (proposal F1 confirmed on the Mac); rotation and orthonormalization are the two big stages, the serial EVD grows to 14 percent at T=8. The MKL split needs dev08 part A |
| B1-SELECTED-ROOTS | Perf sub-session (proposal B1): Ritz vectors, H images and residual vectors only for the roots that can receive a correction in narrow iterations (fewer than half the bands; restart, return, limit and wide iterations keep the all-band path), `BAANDI_B1=standard|strict|on` (default off, bit-identical when off; inactive with `BAANDI_SECOND_VECTOR`). Standard: exact Ritz-change flags. Strict: fixed-seed sign sketch (s = min(128, capacity), A = S V and B = S HV per appended column, rotated on contraction restarts), est_j = ||B c_j - theta_j A c_j||, est below half the threshold flags the root converged for that iteration, every other root gets an explicit residual; returns decided on the all-band path with the unchanged certification. Second sketch form `BAANDI_B1_SKETCH=count` (CountSketch, N adds instead of s N per column); S = identity when the space is not larger than the sketch. Counters in `convergence.work.b1`. Code: real_davidson.rs `solve_research`, `Sketch`; davidson_stages.rs | implemented, opt-in; Mac mixed (dense sketch roughly neutral, count sketch positive); dev08 A/B pre-registered (requests/b1-b4-b5-davidson.md parts B and addendum) | Ledger PRT1 premise change: flags are no longer stale (sketch plus explicit residuals; 0 sketch disagreements in every Mac run), and the target is si32, not si8. Mac si32 strict T=1 seed 601: Davidson 446 / H 6633 unchanged (one dense run 6634), energy within 3.4e-13 Ha, 317 (dense) and 292 (count) of 446 iterations narrow, columns formed 16675 against 30328; rotate + residual + sketch: off 10.6-14.0 s, dense 9.6-10.6 (sketch 2.2-2.6), count 6.6-9.4 (sketch 0.16-0.24); summed iteration time off median 28.7 s, count 21.5/25.3/30.7. si16 strict: count 0.64 s vs off 0.77-0.82, dense 0.95-0.98 (worse). Standard: si32 rotate 0.81-1.04 -> 0.57-0.63 s (15 of 45 narrow), energies within 8.5e-14 Ha. Tests: `tests/b1_selected_roots.rs` (decisions, H work, certification, thread identity on kinetic and clustered fixtures; generality set `--ignored` in release); the existing real Davidson tests pass with `BAANDI_B1=on` | In strict mode the explicit set is much larger than the corrected set (on average 25 of 68 roots are formed in a narrow iteration, while 265 of 446 iterations correct a single root): roots that sit between half and the full threshold need explicit residuals. The dense sketch's maintenance (2 s N flops per appended column, about 40 GFLOP per si32 run) cancels most of the rotation saving on the Mac; the count form removes that cost. Generality fixtures have tiny bases and see few narrow iterations |
| B4-INNER | Perf sub-session (proposal B4): sequential preconditioned directions with a small per-root Rayleigh-Ritz for narrow blocks (at most 2 selected roots, up to 4 inner directions), evaluated without engine code in the `#[ignore]` harness `tests/b4_inner_iterations.rs` (one Davidson written for the test, baseline and variant on the same kinetic, clustered and edge-gap fixtures; the frozen-Hamiltonian replay captures only the complex path, so no real si32 tail capture exists) | rejected on the Mac fixture as pre-registered | No fixture reached the 30 percent single-column filter of the rule (baseline single-column share 0 to 10 percent on 10 fixtures, including 4 tail fixtures with a 1e-4 to 2e-3 Ha gap at the strict edge). Where narrow steps exist the variant cut full steps 3 to 9 percent (inner = 4: 90 -> 82, 82 -> 77, 78 -> 73, 94 -> 84) at +0 to +2 percent H columns; modelled H + 33 x steps at most -6 percent (rule: steps -25 percent, cost -15 percent) | Synthetic dense fixtures do not reproduce the si32 single-column tail (265 of 446 width-1 iterations there); revisit only with a captured real Gamma tail solve (a real-path capture in the replay module would be needed) |
| B5-BORDERED | Perf sub-session (proposal B5): bordered (secular-equation) update of the reduced eigenproblem when exactly one column was appended, `BAANDI_B5=on` (default off, bit-identical when off): arrowhead in the previous Ritz basis, deflation of negligible weights and close pole pairs (Givens, tolerance 8 eps scale), roots relative to the nearer pole by safeguarded Newton, eigenvectors from Loewner-recomputed weights, C' = [[C, 0], [0, 1]] Y by the fixed-order GEMM, sign convention of `symmetric_eigh`; full solve on any doubt (no convergence, non-finite, broken interlacing, recomputed weights off by more than 1e-9 ||z||) and after 8 consecutive updates. Counters `convergence.work.b5`. Code: bordered.rs | implemented (gate passed: eigen share above 5 percent on strict cells), Mac effect inside noise; dev08 A/B pre-registered (part C) | Gate: in-situ eigen share 5.6 (si32 strict T=1), 8.6 (si16 strict), 10.1 (si8 strict), 14.3 percent (si32 strict T=8); 2.3 to 2.4 standard. Microbenchmark (faer): bordered 0.48 to 0.68 of the full solve at orders 96 to 408 (0.60 vs 1.24 ms, 7.64 vs 11.25 ms); the m x m GEMM dominates the update. In the solver: si32 strict 242 of 446 eigenproblems bordered, 0 declined; energies within 1.7e-13 Ha; work unchanged on si8-si32, one Davidson iteration moved on three synthetic fixtures (energies within 3e-15 Ha); eigen stage 1.53/1.58 s vs off 1.42/1.64 s (no visible change under load). Tests: unit test against the full solve on clustered, degenerate, tiny- and zero-coupling cases (eigen-residual and orthonormality 1e-13 x m relative); integration test in `tests/b5_reduced_eigen.rs` (roots 1e-10 Ha, certification, thread identity) | Expected saving about 1 percent of Davidson time at T=1 and 3 percent at T=8 on the Mac; on dev08 the MKL dgemm makes the update relatively cheaper against dsyevd, so the A/B there decides. Exact degeneracies (cubic Si clusters) are handled by deflation, never declined in these runs |
| TRACK-A0 | Track A stage 0/1 skeleton (speed-restructure proposal section 3): all-band direct minimization of the KS energy over the occupied subspace (variant A0: real Gamma pair density, existing Hartree/XC, `RealHamiltonian::apply_block`, projected gradient, TPA-preconditioned Polak-Ribiere CG with reset, Loewdin retraction, energy-decrease line search) as a pre-converger with hand-off to the unchanged strict session (`scf::direct_min`, `ScfMethod::DirectThenScf`, default off, no CLI key) | STOPPED by the pre-registered kill rule (Mac counts); module kept default-off for research | Mac release T=1, gradient tol 1e-6 Ha. Unperturbed total H columns A (minimizer + hand-off rec+cert) vs strict (rec+cert): si-skew 322 vs 614 (0.52), si8 1144 vs 1319 (0.87), si16 3157 vs 1953 (1.62), si32 6011 vs 3934 (1.53), mgo8 3647 vs 2346 (1.56), syn-skew2 259 vs 187 (1.39), syn-fcc 1312 vs 685 (1.92), syn-lowsym6 1460 vs 835 (1.75). Seeds 601-603 medians: si32 6556 vs 5204 (1.26), si16 3395 vs 2261 (1.50), syn-skew2 252 vs 187, syn-fcc 1335 vs 707, syn-lowsym6 1440 vs 835, mgo8 3257 vs 2346 (1.35 to 1.89). Minimizer alone 46-80 H-eq per occupied band on Si/MgO (atomic start), 100-280 on synthetic (free-electron start). Certified energies of both paths within 2.6e-12 Ha on si32, up to 3.6e-11 on si16 seed 601 and 8.8e-12 on mgo8 (both certified). Determinism: minimizer energy bits and densities identical at T=1/4/8. Wall (Mac observation, not evidence): si32 seeds 17.8/33.8/20.6 s vs strict 29.4/29.1/27.5 s | A0 needs about 2 full-width applies per iteration and 26-40 iterations on real fixtures (linear rate about 0.9 per iteration from a free-electron start); the line search hits the energy rounding floor near 1e-6 Ha gradient on 32-64 bands; a minimizer density error above the 1e-7 certificate window makes the hand-off certificate fail repeatedly (up to 2582 certificate columns on si32 seed 602). The dense algebra is cheap, hence the wall observation | Do not re-open A0. A new registration needs a variant with about one apply per iteration (A1/A2) and a cheaper hand-off (buffer bands), and must halve the minimizer cost on si32 to pass the count rule; the 1e-12 cross-path energy gate would also have to be met. Fixed-occupation minimizers reach non-aufbau states on small-gap cells (si-elongated, si-slab) | docs/performance/track-a-stage0.md, tests/reference/perf/track-a/, docs/performance/requests/track-a-stage1.md |
| ADOPT-DEFAULTS-1 | Orchestrator, first application of the user's adoption rule (AGENTS.md 'Adoption rule', 2026-11-10): BAANDI_CERT_SCREEN, BAANDI_NL_ENERGY=fast and BAANDI_MKL_NONLOCAL are now ON by default (rollback: BAANDI_CERT_SCREEN=off, BAANDI_NL_ENERGY=exact, BAANDI_MKL_NONLOCAL=off) | adopted | dev08 evidence in SINGLE-NODE-AB1 and NONLOCAL-DGEMM-AB: screen bit-identical in 40 of 40 pairs and -19.2 percent (T=1) / -15.7 (T=8) on si32 strict; NL_ENERGY energies within 6.3e-13 Ha, wall better in all eight cells, both together -24.2 / -17.7 percent; MKL nonlocal 40 of 40 pairs faster; check.sh PASS (four) with the new defaults | The earlier percentage-threshold rejections are superseded by the adoption rule; the small-cell screen cost was halved by the parallel screen (POTENTIAL-REUSE-IMPL); a fresh-seed dev08 confirmation is still requested to catch regressions, not to veto |
| DEFAULTS-CONFIRM | dev08 fresh-seed confirmation of the new defaults (relayed 2026-11-10, master; 0bce368, binary 173cb73941f722b210e4d9e81e6b7d0769b9521d897a1ee1d5848425fce17956, seeds 701-710, interleaved off/on, T=1 core 0, T=8 cores 0-7) | confirmed: CERT_SCREEN + NL_ENERGY=fast adopted (default ON since a50247d) | median wall off to on: si32 strict T=1 39.191 to 29.720 s (-24.16 percent, worst seed -20.61), T=8 15.668 to 12.640 (-19.33, worst -15.90); si16 strict T=1 5.085 to 4.901 (-3.61), T=8 2.903 to 2.908 (+0.16, worst +0.72; the old +1.36 percent clause is fixed); si32 standard T=1 -12.72, T=8 -4.31; si16 standard T=1 -7.16, T=8 -1.33 (worst +0.45); 160 runs, zero identity violations: density SHA equal off/on in every run, energies within 1e-12 Ha, SCF counts equal, Davidson work never above off; density SHA equal between T=1 and T=8 on 10 of 10 seeds in all eight cells; attribution (si32 strict, seeds 701-705): CERT_SCREEN alone -18.05 percent (T=1) / -17.07 (T=8), NL_ENERGY alone -5.79 / -1.85, both -23.11 / -18.99, total_energy phase -62.7 / -48.0 percent; POTENTIAL_REUSE isolation: +0.24 percent (si16 strict T=8) and +0.10 (si32 strict T=1), bit-equal in 10 of 10 pairs | The parallel screen FFTs, not the potential reuse, removed the small-cell cost (the screen costs 0.19 to 0.20 s per certificate on si32, three certificates per run; the certificate's own potential build is 0.0007 s); the screen's bound matches the actual change to five digits | Both defaults stay ON. BAANDI_POTENTIAL_REUSE stays on (bit-identical, small positive) |
| MPI-M2D | M2d sub-session: (1) batched FFT transposes of the stick-distributed H apply (one alltoall per direction for all column pairs of a block, `DistributedPairFft::*_batch`, default on, `BAANDI_M2_BATCH=off` for A/B); (2) distributed SCF stages on the M2 path (`BAANDI_MPI_M2=on`, 2+ ranks): density from the distributed eigenvectors on FFT planes (one batched transpose plus one allgatherv), LDA XC by grid blocks (one tree_sum plus one allgatherv), audit nonlocal energy from distributed projections (one tree_sum), certificate-screen band integrals (one tree_sum); mixing and Hartree stay replicated; (3) `parallel.m2` counters and stage timers in results.json | implemented; Mac correctness only; dev08 request docs/performance/requests/mpi-m2d.md (adoption rule: any improvement over M1 is adopted) | alltoall per apply of m columns 2*ceil(m/2)+2 -> 4 (m=16: 18 -> 4; Davidson gate n=81: 7.8 -> 4.2 per iteration), batched bit-identical to per-pair on 3 cells x m in {1,3,16} x ranks {1,2,3,5,8}; band densities vs serial full-grid FFT <= 7.4e-16 relative and bit-identical across 1..6 ranks; nonlocal band energies <= 7.7e-16 relative; engine SCF cubic and skewed, 2-4 simulated ranks: strict dE <= 6.7e-16 Ha, standard dE <= 6.2e-11 Ha vs threads-only, bit-reproducible at fixed rank count, all ranks bit-equal; H-apply and Davidson gates of M2b/M2c unchanged | Mixing was not distributed on purpose (0.12 s of 41 s at np=1 in MPI-M1-GATE4; a grid allgatherv per step would cost about what it saves); GGA, NLCC, density floors, finite temperature and k meshes keep the replicated stages. Expected gain on dev08: the replicated stages were 7.4 s (np=1) growing to 10.3 s (np=8) in MPI-M1-GATE5; density and XC now scale with ranks, transposes per apply fall from up to 70 to 4 at 68 columns |

## Micro-stack (2026-11-10; implementer, Mac observations; dev08 request docs/performance/requests/micro-stack.md)

Mac = shared Apple Silicon, release pure Rust, perturbed seeds 601-602 (+-1e-7 Bohr), other agents' builds running (wall noise +-20 percent); phase timers and counts are the evidence, wall is not. Profile (sample, si32 strict T=1, 23 s): nano-gemm/faer GEMM about 55 percent, FFT about 13, `Projection::run_into` self time 8.7 (scalar MGS dot/axpy), `NonlocalOperator::expectation_values` 2.2, atomic start about 2.7 (first-step `potential` phase; the per-step potential build is only about 5 ms), density 2.8.

| ID | Idea | Status | Measured (Mac) | Lesson / next |
|---|---|---|---|---|
| MICRO-DENSITY-REDUCE | `BAANDI_DENSITY_REDUCE=on`: density band sum over fixed 4096-point blocks, pool-parallel, ascending band order per point; zero-occupation grids cleared once and skipped | opt-in, **bit-identical** (unit test `pair_and_reduce_switches`, T=1/2/3/4/8) | si16 T=1 and T=8 seeds 601-602: energy and density SHA equal; density phase T=8 0.041 -> 0.035-0.039 s (noise) | Candidate for immediate default if dev08 wall is not worse |
| MICRO-APPLY-LEAN | `BAANDI_APPLY_LEAN=on`: real H apply kinetic term without a per-element integer modulo; nonlocal add without repeated input scans (X scanned once, Y once at the end) | opt-in, **bit-identical** (test `apply_lean_is_bit_identical`, incl. NaN rejection) | no Mac-visible change (expected 0.5 percent of H time or less) | Same as DENSITY-REDUCE |
| MICRO-B6 | `BAANDI_B6=on` (`ResearchOptions::b6`): solver half of B6-PHYSICS: cache P^T V per basis column (appends, contractions, rebuilds), return certificate from (P^T V) C via `apply_block_with_projections` | opt-in; eigenpairs, density, energy and counts bit-identical, certificate residual digits move | si16 seed 601: residuals relative 4e-7 of 1e-10 values, H counts 2162+253 both; test `b6_certificate.rs` (every certificate cached, contracted restarts covered, thread-independent) | Saving is one projector product per strict solve (B6 estimate 0.3-0.6 percent); dev08 decides |
| MICRO-NL-ENERGY-REAL | `BAANDI_NL_ENERGY_REAL=on`: audit nonlocal energy on the real Gamma path from one real projector GEMM of the returned ORC vectors (`RealNonlocal::expectation_values`) instead of complex projector sums | opt-in, rounding class | `total_energy` phase si32 seed 601 0.586 -> 0.078 s (-87 percent, about 2 percent of wall), si16 0.092 -> 0.027 s; energy 3e-14 Ha, density SHA and counts identical; test `real_nonlocal_expectation_matches_complex` | Largest item of the stack on the Mac; stacks with NL_ENERGY=fast (replaces it on the real path) |
| MICRO-ATOMIC-REAL | `BAANDI_ATOMIC_REAL=on`: atomic-start Ritz H Q through the real Gamma Hamiltonian; the real state built there is reused by the first solve (shared `ensure_state`) | opt-in, rounding class (trajectory moves at rounding) | `potential` phase si32 0.677 -> 0.493 s, si16 0.160 -> 0.102 s; si32 energy 1.8e-13 Ha, SCF/Davidson/H counts identical (308, 4909+443); test `atomic_apply_real_matches_complex` | The remaining atomic-start cost is the complex MGS and complex Gram (real versions are a follow-up) |
| MICRO-DENSITY-PAIR | `BAANDI_DENSITY_PAIR=on`: DENSITY-PAIR implemented with a split kernel (`SplitPairDensity`) so each band keeps its contribution grid (CERT_SCREEN integrals unchanged); fixed pairs (2k, 2k+1), even-aligned worker ranges | opt-in, rounding class, thread-count independent (test, T=1..8 bit-equal) | si16 seed 601 energy 7e-15 Ha; density timer inside noise | Premise (a) of the DENSITY-PAIR row removed by the split kernel, (b) by fixed pairs; gain is half the density FFTs (about 1 percent) |
| MICRO-WIDE-DOT | `BAANDI_WIDE_DOT=on`: eight-lane dot in the in-block MGS | opt-in, rounding class | microbench 5 -> 10 elements/ns, but si32 `run_into` self share unchanged (8.7 -> 9.6 percent incl. dot, noise): MGS is bound by L2 traffic on the M-series core | Kept only for the x86 measurement; drop if dev08 shows nothing |
| MICRO-VEC-DOT / MGS-FUSED | bit-identical lane-exact rewrites of `dot`/`norm` (`chunks_exact`) and a fused update-plus-next-dot MGS sweep | rejected (not committed) | slower on the Mac: the loop vectorizer turns the four ordered lane sums into scalar in-order reductions (run_into+dot 1119 -> 1725 samples) | A four-lane ordered sum is latency bound; only more lanes (rounding change) can be faster |
| MICRO-ORTHO-DGKS | `BAANDI_ORTHO=dgks` / `dgks:<ratio>` (coordinator priority after the dev08 stage profile: orthonormalize 33-45 percent of non-H Davidson time): skip the second projection pass of an append when every column kept at least `ratio` of its unit norm after the first pass (DGKS); contraction orthogonality gate unchanged | opt-in, rounding class | si32 seed 601 T=1: ratio 1/sqrt2 skips 10 of 313 appends (diag 18.88 -> 18.50 s); ratio 0.1 skips 256 of 312 (diag 18.76 -> 16.05 s, -14 percent; blocked MGS 17.36 -> 15.55 s), energy within 4.4e-13 Ha, counts unchanged; test `dgks_projection_stays_orthonormal` (defect <= 1e-13, classic and blocked, repeat path) | Davidson corrections usually lose more than 30 percent of their norm to the basis, so the textbook 1/sqrt2 rarely fires; the loss-of-orthogonality bound eps sqrt(m)/ratio allows 0.1 with four orders of margin to the 1e-12 gate |
| MICRO-ARM1 | dev08 arm 1 of the micro-stack request (relayed 2026-11-11, master; fa629c1, binary 05e998150bd7c70c968c2ce0e76947d705d73f9a2b4427c922c2f5d99c15e575, seeds 601-610, interleaved off/on, C1 ladder and the other three switches at their defaults in both arms): BAANDI_DENSITY_REDUCE + BAANDI_APPLY_LEAN + BAANDI_B6 | adopted: all three now default ON (rollback `=off`) | bit-identical in all 160 runs (energy bits, density SHA, Davidson and H counts, certificate verdicts, SCF counts; density SHA equal between T=1 and T=8); median wall off to on: si32 strict T=1 -0.25 percent, T=8 -2.35; si32 standard T=1 -0.68, T=8 -3.05; si16 strict T=1 +0.12 (worst +2.21), T=8 -1.45; si16 standard T=1 -1.31, T=8 -2.20; seven of eight cells improve, worst seed over all cells +2.21 percent at T=1 (off-arm spread 13.8 to 15.5 percent on the strict cells); density phase alone -6.6 to -7.9 percent per seed at T=1 on si32 strict | The gain is thread-dependent (DENSITY_REDUCE parallelises a serial band sum); the strict-cell seed spread with the C1 ladder on is larger than the effects, so only bit-identical or paired-per-seed comparisons are resolvable for small changes; the rounding-class switches and ORTHO need the paired per-seed deltas reported alongside the spread |
| MICRO-ARM2 | dev08 rounding-class switches alone (relayed 2026-11-11, master; fa629c1, MKL, si32 strict and standard, T=1 and T=8, seeds 601-605, interleaved, paired per-seed sign tests; every run converged, every |dE| inside its class, zero violations in 160 runs) | adopted: BAANDI_NL_ENERGY_REAL, BAANDI_ATOMIC_REAL, BAANDI_DENSITY_PAIR, BAANDI_WIDE_DOT default ON (rollback `=off`) | cell medians off to on, strict T=1 / strict T=8 / standard T=1 / standard T=8: NL_ENERGY_REAL -3.82 / -1.27 / -8.32 / -1.93 percent; ATOMIC_REAL -1.04 / -1.18 / -4.29 / -2.01; DENSITY_PAIR -0.90 / -0.76 / -3.60 / -1.60; WIDE_DOT -0.82 / -2.09 / -0.80 / -1.67 (16 of 16 cell medians improve); paired: NL_ENERGY_REAL 5 of 5 seeds in all four cells (standard T=1 -8.64 percent, the largest single effect), ATOMIC_REAL and DENSITY_PAIR 5 of 5 in six cells, 4 of 5 in strict T=8 (unresolved but negative); WIDE_DOT 5 of 5 on strict T=8 (-2.06, the best T=8 switch) and on three more cells, 4 of 5 on standard T=1; also arm-1 sign test: bit-identical stack 9 or 10 of 10 seeds in seven cells, si16 strict T=1 unresolved (3 of 10, +0.43 percent) | WIDE_DOT does pay on x86 although the Mac saw nothing (four-lane accumulators); with the ladder on, strict-cell seed spread is 11 to 13 percent, so only paired comparisons resolve small effects | No regression in any cell; unresolved cells are still negative |
| MICRO-ORTHO | dev08 BAANDI_ORTHO A/B (relayed 2026-11-11, master; fa629c1, MKL seeds 601-610, pure seeds 601-605, interleaved, load gated; stage timers on for the mechanism) | adopted: `dgks:0.1` is now the default (`BAANDI_ORTHO=off` restores two passes; plain `dgks` at 1/sqrt(2) is marginal and not the default) | dgks:0.1 median wall off to on (worst seed): si32 strict T=1 27.654 to 24.990 s -9.64 percent (-9.13), T=8 -9.69 (-8.69); si32 standard T=1 -6.55, T=8 -6.57; si16 strict T=1 -8.07, T=8 -8.61; si16 standard T=1 -4.87, T=8 -5.53 (worst seed negative in all eight cells); pure si32 strict T=1 42.703 to 37.712 s -11.69 percent, T=8 16.245 to 13.669 -15.86; plain dgks: -0.2 to -3.1 percent, worst seed positive in four cells; mechanism: orthonormalize stage -43 percent (T=1) and -44 (T=8), second passes skipped 208-235 and repeated 60-118 per run at 0.1 against 8-9 skipped at 1/sqrt(2); Davidson and H counts identical across all arms on all three seeds (310/4749, 240/4662, 317/4675), so the gain is removed work, not a changed path; the part-A profile (orthonormalize 33.1 percent) predicts -14 percent of Davidson time, measured -11.0 | The next target by the same profile is h_apply (38 percent of Davidson time on si32 strict T=1); the orthonormalize stage still carries about 19 percent after the cut |

## H application (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/h-apply.md)

Trigger: dev08 stage profile (MKL, si32 strict, new defaults): `h_apply` 38 percent of Davidson time at T=1, 33 at T=8; MPI-M1-CEILING: 66 percent of apply calls carry one column, 65-68-column calls carry 69 percent of the column work. On MKL builds the pair local part is `GridFusedPlan` (full-grid 3D DFTI, full zero fill), so a single column pays a whole complex transform for one real field.

| ID | Experiment | Status | Measured result | Lesson / limits |
|---|---|---|---|---|
| HAPPLY-OVERLAP | `BAANDI_HAPPLY_OVERLAP=on`: local pair FFTs and the nonlocal image (P^T X, coupling, P W into a private buffer) as the two branches of one pool join, then Y += kinetic and Y += image in the serial order (plain and B6 cached-projection apply) | opt-in, **bit-identical** (test `apply_overlap_is_bit_identical`, widths 1/2/5/8, T=1/2/4) | Mac microbench (60^3, 12053 rows, 256 projectors): T=4 width 1 0.724, width 2 0.695, width 16 1.037, width 68 1.006 of base; T=1 0.995 to 1.011. SCF si8/si16 strict/standard T=1/T=4: energy bits, density SHA and counts equal in all 32 comparisons; si8 strict T=4 wall 0.954/0.987 | For 1-2 columns the local part is one serial pair piece while the nonlocal GEMM panels are parallel: running them side by side hides most of the nonlocal time at T>1. Wide blocks already fill the pool; neutral there. Not a repeat of A33/A35 (no new threads; the persistent pool) |
| HAPPLY-FOLD | `BAANDI_HAPPLY_FOLD=on`: odd tail column (every single-column call) through a half-length complex 3D transform of the z-folded real field (`Z = psi(2m) + i psi(2m+1)`, potential split into even/odd z samples, unfold by the two-real-halves identity); pair path kept when the last grid axis is odd | opt-in, rounding class (tail relative 1.0e-16 against the zero-partner pair; thread-count independent) | Same-FFT full-grid proxy (rustfft in both arms, stand-in for the MKL GridFusedPlan): 48^3 0.381, 60^3 0.459, 72^3 0.520 of the pair time. Pure build (pruned pair plan as base): neutral (width 1 1.034 T=1, 1.005 T=4). SCF |dE| 2.8e-14 to 5.5e-13 Ha, counts identical, density SHA T=1 = T=4 | The real-to-complex saving is only available to a single real field; pairing (PFFT-1) already gives it to two. Gains only where the pair path is a full-grid transform (MKL build); a pruned folded plan for the pure build is not implemented |
| HAPPLY-PAIR-MKL | `BAANDI_PAIR_FFT=pruned` (MKL builds): `PairFusedPlan` mirror-slab layout on MKL DFTI 1D batch descriptors (`dft_mkl::Plan1dBatch`, one per lane buffer, contiguous, unscaled) instead of `GridFusedPlan`; backlog B33 step 2 | opt-in, rounding class; MKL code type-checked only (no MKL on the Mac); unit test `mkl_lanes_match_rust_lanes` for dev08 | no Mac number possible; estimate PFFT-5 (MKL stick pipeline 97 / 235 us per transform equivalent against MKL raw full 147 / 323 us, relayed P39B) | Pluggable `PairLanes` keeps the pure path bit-identical (same rustfft calls). Applies to all widths, so it is the only item here aimed at the wide 65-68-column calls |
| HAPPLY-NOT-DONE | (2) fusing pack with the scatter / unpack with the gather of the pair map, (3) a GEMV-shaped nonlocal product for one column, (5) per-size choice between MKL 3D and rustfft | not implemented | PFFT-1: pack plus unpack 1.5 percent of the pair apply; microbench: one column's nonlocal (256 projectors x 12053 rows, two products) is about 0.6 ms, the cost of streaming the 24.7 MB projector matrix twice (memory bound, a GEMV shape reads the same bytes); P39B/PFFT-5: MKL raw 3D beats the rustfft mirror plan on x86, so (5) reduces to HAPPLY-PAIR-MKL | Remaining nonlocal lever for single columns is fewer bytes (real-space projectors or fp32 storage), not a different BLAS shape |
## Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md)

Mac = shared Apple Silicon (load 4 to 11), release pure Rust (faer, no MKL), si32-gamma perturbed seed 601; stage seconds normalised by the same run's h_apply stage (untouched by every switch); counts and SHA exact. Kernel microbenchmark: `cargo test --release -p dft-numerics --test algebra_stream stream_timing -- --ignored --nocapture`. Tests: `dft-numerics/tests/algebra_stream.rs`, `dft-solvers/tests/algebra_stream.rs` (roots within 1e-11 Ha of the default path and 1e-10 of known spectra, certified strict returns, orthogonality defect 1e-12, thread identity at 1/2/4/8 workers).

| ID | Idea | Status | Measured (Mac) | Lesson / next |
|---|---|---|---|---|
| ALGEBRA-STREAM-FUSED | `BAANDI_RITZ_FUSED=on`: Ritz vectors and residuals in one pass over V and HV (each panel of HV C lives only in its tile, combined at once with the panel of X); H X formed late (single-side rotation, same kernel) only for a contraction restart; also used for the B1 subset | opt-in, **bit-identical** (faer tested; MKL same dgemm calls and tiles, dev08 confirms) | microbench 68 bands, m 408: 43.2 -> 43.1 ms T=1, 10.7 -> 10.3 ms T=8; in solver rotate+residual / H 1.421 -> 1.467 (T=1), 1.136 -> 1.204 (T=8), inside the noise; SHA equal to OFF on si16/si32 strict and standard | faer rotation is compute bound (about 28 GFLOP/s); the saved 3 n x bands streams (about 20 MB per si32 iteration) matter only where the GEMM is fast (MKL); dev08 decides |
| ALGEBRA-STREAM-NARROW-PROJ | `BAANDI_NARROW_PROJ=on`: correction blocks of at most 4 columns projected by streaming eight-lane dots (S = Q^T X, no tiles) and one fused update pass X -= Q S (fixed 1024-row blocks, no FMA) | opt-in, rounding class, thread-independent | microbench k=1: 3.64 -> 1.60 ms (T=1), 0.99 -> 0.52 (T=8); k=2 5.33 -> 2.19; k=4 5.15 -> 4.01; k=8 7.85 -> 8.03 (hence width 4); in solver orthonormalize / H -7 percent (T=1), -5 (T=8); energies within 6e-13 Ha, counts identical | CGS needs Q^T x before the update, so two reads of Q per pass are the floor; fusing the next reduced border into the projection is impossible (it needs H v_new); the symmetric rewrite moves the stream to HV and adds a 1/l cancellation |
| ALGEBRA-STREAM-NARROW-REDUCED | `BAANDI_NARROW_REDUCED=on`: reduced-matrix border for at most 4 new columns by the same streaming dot kernel | opt-in, rounding class | microbench k=1 1.11 -> 0.56 ms (T=1), 0.40 -> 0.24 (T=8); in solver reduced / H -1 to -5 percent | small stage (4 percent dev08) |
| ALGEBRA-STREAM-B5-PAR | `BAANDI_B5=par`: B5 with secular roots, Loewner weights and eigenvector columns per root on the pool and the m x m coefficient product in fixed 32-column groups | opt-in; identical to `B5=on` bit for bit on faer (test), MKL a fixed different dgemm grouping | `bordered_timing` order 408: full 14.3, bordered serial 10.3, par T=4 4.0, par T=8 3.2 ms; in solver si32 strict T=8 eigen / H 0.654 -> 0.589 (B5=on) -> 0.492 (par), 224 of 402 eigenproblems bordered, 0 declined | B5 part C still pending on dev08; `par` is the T=8 form; range-selected dsyevr (MKL-EVD) and threaded dsyevd (thread-dependent) stay excluded |
| ALGEBRA-STREAM-B1-MAC | Mac re-check of `BAANDI_B1=on BAANDI_B1_SKETCH=count` against the new defaults (DGKS 0.1, B6, ladder) | stays opt-in (dev08 part B pending; adoption needs dev08 wall) | si32 strict non-H Davidson time / H 2.305 -> 1.904 (T=1, -17 percent), 2.907 -> 2.680 (T=8, -8); standard 0.668 -> 0.468 (T=1), 1.387 -> 1.187 (T=8); 245 of 362 strict and 15 of 45 standard iterations narrow; counts identical, energies within 3.4e-13 Ha | Largest Mac effect of the set; with the stack (fused, narrow, B5=par) T=8 strict 2.535 (-13 percent) |
| PRE-EXISTING TEST FAILURE | `dft-engine scf::tests::energy_identity_matches_explicit_audit` fails on main bc33efe (final total 1 ulp apart between identity on and off); passes at dff2fe8 and with `BAANDI_ORTHO=off` | reported, not fixed (outside this area) | deterministic; caused by the DGKS default commit, not by the algebra-stream branch (the branch's default path is bit-identical to bc33efe) | owner of scf.rs / DGKS default to decide: bitwise assertion vs. rounding class |
| CAMPAIGN-1 | dev08 cumulative campaign of the 2026-11-10 to 11-11 adoption round (relayed 2026-11-11, master; fa629c1, one binary, three arms: OLD = all new defaults switched off, NOW = defaults at fa629c1 (C1 ladder, CERT_SCREEN, NL_ENERGY fast, MKL nonlocal), BEST = defaults plus the eight winners; seeds 601-605; winners-together A/B on seeds 601-610) | adopted: the eight winners (DENSITY_REDUCE, APPLY_LEAN, B6, NL_ENERGY_REAL, ATOMIC_REAL, DENSITY_PAIR, WIDE_DOT, ORTHO=dgks:0.1) are default ON since bc33efe | BEST against OLD, median wall: si32 strict T=1 40.652 to 22.768 s (-43.99 percent, 1.79x), T=8 16.086 to 10.662 (-33.72); si32 standard T=1 13.402 to 8.438 (-37.04), T=8 5.781 to 4.589 (-20.60); si16 strict T=1 4.946 to 3.795 (-23.26), T=8 2.783 to 2.412 (-13.33); si16 standard T=1 2.612 to 1.934 (-25.94), T=8 1.552 to 1.356 (-12.61); BEST against NOW -10.97 to -24.32 percent in all eight cells; winners together on seeds 601-610: -11.25 to -24.25 percent MKL, pure si32 strict -15.67 (T=1) / -21.20 (T=8), every seed in every cell; accuracy over the chain: max |E(BEST)-E(OLD)| 9.1e-13 (si32 strict), 5.7e-13 (si32 standard), 3.5e-12 (si16 strict), 1.3e-13 (si16 standard) Ha, inside the strict class; zero violations; density SHA equal between T=1 and T=8 on 10 of 10 seeds in all eight cells; the stack is close to additive (individual si32 strict T=1 effects sum to about -16.5 percent against a measured -17.38) | Only cell where NOW was not faster than OLD: si16 strict T=8 (2.808 against 2.783, the residual small-cell screen cost), repaid by the winners; phase attribution on pure si32 strict T=1: density -36 to -37 percent, total_energy -79 to -80, diagonalization -12 to -15, potential unchanged; next target by the stage profile: h_apply (38 percent of Davidson time on si32 strict T=1); nothing in this stack changed a trajectory (Davidson and H counts equal), so the 14 to 16 percent off-arm seed spread only matters for trajectory-changing changes (C1, v4, B1) | si32 strict T=1 is 1.79x faster than at the start of the round; changes still pending dev08: H apply (OVERLAP, FOLD, PAIR_FFT), algebra streaming (NARROW_PROJ, RITZ_FUSED, B5 par, NARROW_REDUCED), B1, MPI M2/M2d, Track A wall, v4, C3 |
| B1-PART-B | dev08 B1 A/B (relayed 2026-11-11, master; fa629c1, current defaults, winners not set, three settings off/dense/count interleaved per seed, seeds 601-610, 240 runs; generality cells built with scripts/bench_suite.py helpers) | adopted: B1 default ON with the count sketch (`BAANDI_B1=off` disables, `BAANDI_B1_SKETCH=dense` selects the dense sketch, which is a regression on small cells) | total Davidson and H vectors identical to off in all 240 runs, zero accuracy violations, sketch disagreements 0, all converged; narrow iterations 68.0 percent (count) / 70.9 (dense) on si32 strict T=1, 32 to 38 percent elsewhere; paired wall off to count: si32 strict T=1 -4.24 percent (10/10), T=8 -2.00 (10/10); si32 standard T=1 -1.74 (10/10), T=8 -0.75 (10/10); si16 standard T=1 -1.40 (9/10), T=8 -1.27 (8/10 unresolved); si16 strict T=1 +1.13 (4/10), T=8 +0.89 (3/10) both unresolved; dense: si16 strict T=1 +6.61 and T=8 +2.70 with 0/10 seeds improving (a regression), si32 standard -1.89/-1.18 (10/10) | Generality: syn-spd, syn-p3, syn-d3 and si-skew are identical to off (energies bit-equal except si-skew standard 2.66e-15 Ha); si-slab does not converge in any arm and takes a different non-converged path with B1 on, not a B1 failure | si16 strict is the only cell where count is not clearly good (unresolved, small); keep watching with 30 seeds if it matters |
| B5-PART-C | dev08 B5 A/B (relayed 2026-11-11, master; fa629c1, current defaults, interleaved per seed, seeds 601-610, 160 runs) | adopted by the adoption rule (nothing regresses, one cell resolved better): B5 default ON (`BAANDI_B5=off` disables; `par` is the thread-parallel variant, pending its own A/B); the registered rule-4 stage clause (-25 percent) FAILED and is reported only for the mechanism | zero accuracy violations, all converged; bordered update fires on 53.4 percent of Davidson eigenproblems on si32 strict (36.0 si16 strict, 6.1 si32 standard, 2.4 si16 standard), declined 0.00 percent; eigen stage si32 strict T=1 1.1369 to 1.1088 s (-2.5 percent), T=8 1.2999 to 1.2509 (-3.8); paired wall: si32 strict T=1 -0.77 percent 10/10 (resolved), T=8 -0.16, si32 standard +0.33/-0.32, si16 strict -0.28/-0.58, si16 standard +0.05/-0.14, all unresolved | Key observation: the eigen stage is larger at T=8 (1.30 s) than at T=1 (1.14 s) on the same work because the reduced eigensolve does not thread, so it is a serial bottleneck (6.1 percent of Davidson at T=1, 16.3 at T=8); the lever is making that solve use the threads (BAANDI_B5=par; Mac: si32 strict T=8 eigen stage -25 percent), not doing less work; B5=par goes to dev08 next |
| MPI-M2D-AB | dev08 MPI M2d scaling (relayed 2026-11-11, master; fa629c1, MPI binary sha256 114a3613a562ffc28434a5a71a962d5b3da70ba91bd5275f3b9bfe68b07b8503, --features mkl,mpi, Open MPI 5.0.8, MKL/OMP threads 1, 8 physical cores of socket 0, si32 strict perturbed seeds 601-610, all configurations interleaved within each seed, load gated) | adopted: BAANDI_MPI_M2 default ON with two or more ranks (`BAANDI_MPI_M2=off` selects M1); the superseded M2e percentage gates are not met but the adoption rule decides | median wall: np=1 T=1 28.080 s; M1 np=4 21.977 (1.28x), np=8 18.366 (1.53x); M2d np=4 T=1 18.784 (1.49x, -14.53 percent vs M1), np=8 T=1 12.782 (2.20x, -30.40 percent vs M1); M2d np=2 x T=4 22.835 (1.23x); np=4 x T=2 17.204 (1.63x); threads np=1 x T=8 12.585 (2.23x); M2 branch proof 10/10 on every M2d row (parallel.model mpi_stick_rows, m2.enabled, counters positive; np=4 example rows_solves 16, batched_applies 358, unbatched 0, density_builds 16, xc_evaluations 34, nonlocal_energy_audits 16, screen_integrals 3, applied_columns 6366); total Davidson 319 in every configuration (trajectory unchanged); max |dE| vs np=1 6.3e-13 to 8.5e-13 Ha; every run converged | Distributed MPI catches the thread path on one node for the first time (np=8 1.6 percent behind threads at the same 8 cores, M1 was 33 percent behind); pure ranks beat hybrids; the earlier gates (np=4 2.2x, np=8 3.5x) are not met (1.49x, 2.20x); harness note: a first sweep without --map-by slot:PE=8 --bind-to core pinned the threads-only rank to one core (29.454 s) and was discarded; still owed: M2_BATCH=off A/B at np=8, si32 standard cells, determinism repeat at np=4, stage/collective breakdown; next MPI questions: np=16 across both sockets against T=16 threads (NUMA), and the remaining replicated stages (Hartree, mixing) |
| EMPTY-LOOSE | dev08 first-empty-band hypothesis and A/B (relayed 2026-11-11, master; 1f449ad, DAVTRACE serial patch in a /scr clone for the trace, binary cli-mkl-1f449ad 1da9563ec97c387f for the A/B, paired off/on, perturbed seeds 601-610, unperturbed median of 3) | adopted: BAANDI_EMPTY_LOOSE default ON (`=off` restores the strict first empty band in every step) | trace, si32 strict T=1: in 100 percent of single-unconverged Davidson iterations the lone band is an EMPTY band, band 64 (the first empty) in 90 to 100 percent, never an occupied band; unperturbed 16 of 99 iterations (16.2 percent) vs perturbed 120-197 of 284-369 (42-53 percent): the symmetry-breaking penalty is this tail (single-unconverged iterations x10, Davidson total 64 to 240-317); standard is insensitive because it converges on Ritz-value change; A/B: si32 strict perturbed T=1 21.767 to 18.923 s -13.06 percent (10/10 seeds, Davidson 319 to 125.5, -60.66), T=8 10.432 to 9.007 -13.66; si32 strict unperturbed T=1 13.354 to 12.845 -3.81 (Davidson 64 to 56), T=8 -3.64; si16 strict perturbed T=1 -3.65 (Davidson 122.5 to 96.5), T=8 -5.42; si16 strict unperturbed -4.71 / -5.28; zero violations: every |dE| inside the strict class, certificate verdicts identical seed by seed, density SHA equal between T=1 and T=8; identity cells: syn-spd, syn-p3, syn-d3 bit-identical, si-skew Davidson 61 to 72 with H 390 to 380, dE 8.9e-15 Ha, verdicts equal | Davidson work falls 60.7 percent but wall only 13.1 because a width-1 iteration costs 28 ms against 208 ms for a full one (stage profile); the QE ratio on the identical displaced geometry improves from 6.08x to 5.28x; the remaining perturbed-over-unperturbed gap is still +47 percent (18.923 against 12.845 s): next question is what the WIDE iterations look like in the perturbed run |
| EIGEN-FACTOR-AB | dev08 A/B of the eigen tolerance schedule coefficient (relayed 2026-11-11, master; 92805b8, binary cli-mkl-ef bcbb1e430c8a3e29, 0.03 and 0.01 against the 0.1 default, paired per run, T=1, load gated; zero correctness violations) | NOT adopted: default stays 0.1 (`BAANDI_EIGEN_FACTOR` stays an experiment switch) | si32 strict perturbed wall 18.996 s to 21.242 (+11.82 percent, 0.03) and 24.396 (+28.43, 0.01), SCF steps 10/10/10, Davidson 125.5/175.5/270; si32 unperturbed 12.824 to 15.034 (+17.24) and 16.125 (+25.74), SCF 8/8/8; si16 strict perturbed 3.543 to 3.753 (+5.93) and 4.012 (+13.24), si16 unperturbed +16.45 / +24.20; si32 0 of 10 seeds improved in both arms, si16 perturbed 1 of 10 and 0 of 10 (p 0.021, 0.002); identity cells go the other way: si-skew at 0.03 SCF 14 to 10, Davidson 72 to 66, wall 0.237 to 0.196 s (-17 percent), syn-spd SCF 11 to 10, wall 0.015 to 0.013 | The stall behind the two extra SCF steps of the perturbed run is NOT caused by a loose eigen tolerance (schedule 10 times tighter at every step, still 10 SCF steps), so cause 1 lies in the density mixing, not in the eigensolve; a constant coefficient helps where the SCF count is the bottleneck (si-skew) and hurts where Davidson is (si32, si16), so no single value fits both; a regime-reading rule (tighten only while the SCF count is still falling) is the shape the data suggests, not a design; next lead: the mixing history/step at the stalled mid-schedule steps of the perturbed run (density residual 1.855e-6 to 1.732e-6) |
| BUFFER-TOL | dev08 A/B of the empty-buffer residual floor in intermediate steps (relayed 2026-11-11, master; f52b79b, binary cli-mkl-bt, paired per run, T=1, load gated, EMPTY_LOOSE and the other defaults on in all arms) | adopted: default 1e-3 (`BAANDI_BUFFER_TOL=1e-5` restores the old floor); T=8 and the SHA cross-check requested | wall 1e-5 to 1e-4 / 1e-3 / 1e-2: si32 strict perturbed 19.034 to 17.883 (-6.05) / 17.598 (-7.54) / 18.041 (-5.21), 10 of 10 seeds in all arms, Davidson 125.5 / 103.5 / 100 / 142; si32 unperturbed 12.802 to +0.24 / +1.22 / +0.57; si16 strict perturbed 3.576 to -2.46 / -2.81 / -1.30 (8, 7, 6 of 10 seeds: unresolved); si16 unperturbed +0.65 / +2.33 / +1.85 | The optimum is interior (1e-3 beats 1e-2, where Davidson rises again to 142); certificates stay strict by construction and by measurement: every run converged, every final performed certificate inside 1e-7 (6.9e-9 to 9.2e-8), final strict residuals 7e-10 to 1e-9, buffer residuals 4.6e-6 to 1e-5, zero failures even at 1e-2; the first-pass 'violations' were a harness artefact (a nonexistent key); NOT verified: the gap check (band_gap is null in results.json for these cells); identity cells: syn-spd/p3/d3 bit-identical, si-skew Davidson 72 to 54 at 1e-3; the perturbed si32 strict T=1 chain: 21.759 (before EMPTY_LOOSE) to 18.923 to 17.598 s, QE ratio 6.08x to 5.28x to 4.91x | Unperturbed cells show +1 to +2 percent at 1e-3 (3 seeds, unresolved), perturbed (realistic) cells gain; chosen for the realistic case |
| C1-LADDER-DEV08 | dev08 A/B of the C1 cutoff ladder (relayed 2026-11-11, master; c70de1f, binary sha256 9fa16ed0544c604af126f58c5f6879979fc78428a034d0c6d0347219c1bda51a, MKL, one binary env A/B, interleaved off/on per seed, taskset core 0 / 0-7, load gated, seeds 601-610; per-run detail /scr/mino/rustdft/c1/runs.jsonl) | adopted: BAANDI_CUTOFF_LADDER default ON (strict real Gamma fixed occupation; standard and k meshes fall back) | at the project protocol (+-1e-7 Bohr): si32 strict T=1 28.901 to 27.717 s (-4.10 percent, worst seed +7.25), T=8 12.694 to 12.742 (+0.38), si16 strict T=1 -4.00 / T=8 +1.59, standard cells bit-identical (fallback); at realistic displacements, si32 strict T=1 seeds 601-603, new defaults: fractional +-0.01 (about 0.21 Bohr) 66.357 to 47.773 s (-28.01 percent, SCF 26-29 to 16-17), fractional +-0.001 (about 0.021 Bohr) 62.244 to 39.787 (-36.08, SCF 23-24 to 11-13); against the old baseline (CERT_SCREEN=off NL_ENERGY=exact MKL_NONLOCAL=off) the project-protocol gain is -3.58 percent, so the defaults do not explain the small gain; ladder block si32: coarse 5 SCF iterations, 40-52 coarse Davidson, 679-685 full-grid-equivalent H, coarse stage 3.59-3.83 s, coarse residual 2.9e-5, grid 54x54x25 against 72x72x36; correctness: zero violations in 160 runs, every ladder run passed the full-cutoff certificate, |dE| 2.8e-14 to 8.4e-12 Ha, density SHA equal between T=1 and T=8 | The +-1e-7 Bohr protocol starts almost converged (12 SCF steps) and hides the ladder; realistic displacements cut 28 to 36 percent; later chain at 4b83503: ladder default on in all later arms (convergence.work excludes the coarse stage: cutoff_ladder.coarse_davidson_iterations, 32 on seed 601) | This row backs the comment in dft-cli config.rs (cutoff_ladder); the ladder strict-only restriction is the open item S1 of the 2026-11-11 work order |
| LADDER-STANDARD | Work order 2026-11-11 S1: C1 ladder for the standard profile, `BAANDI_CUTOFF_LADDER_STANDARD=on` (engine `CutoffLadderOptions.standard`), default OFF. Coarse stage in the standard profile at E1 = 0.5 E2 (uncertified), unchanged stop rule and prolongation, unchanged full-cutoff standard session; first fine eigen tolerance = what the coarse ratchet would use next (`BAANDI_LADDER_STANDARD_SEED=coarse`, default of the switch) or the cap (`=cap`) | implemented, default off; kill rule NOT triggered (S1 continues); dev08 A/B requested (docs/performance/requests/ladder-standard-secant.md); STANDARD-CLASS RISK found on the Mac (see numbers) | Mac release, pure Rust, T=4/8 (thread-independent counts), seeds 601-603, +-1e-7 Bohr (b) and +-0.001 fractional (f); median total H work in full-grid equivalents (fine recorded + fine atomic start + coarse scaled by the grid ratio), plain standard -> standard ladder (coarse seed): si32 b 2059 -> 1761, f 2691 -> 2462 (fine SCF 11 -> 6, 14 -> 9); si16 b 948 -> 990, f 1292 -> 1120; si8 b 670 -> 446, f 799 -> 605; si-skew 227 -> 228 / 227 -> 228; si-low 178 -> 202 / 178 -> 205; mgo8 (PBE, NLCC, floor 1e-10 as in the CLI) 1314 -> 1257 / 1223 -> 1237; syn-skew2 152 -> 130 / 156 -> 136; syn-fcc 593 -> 417 / 630 -> 421; syn-lowsym6 523 -> 645 / 493 -> 444; syn-slab 112 -> 131 / 114 -> 120. Every run converged. |E - E_strict plain| max: plain standard 2.1e-9 Ha; ladder coarse seed 1.18e-8 Ha (syn-slab f seeds 602, 603: 9.9e-9 and 1.18e-8; every other ladder-arm run <= 1.9e-9), ladder cap seed 6.2e-9 Ha (si8 f, one seed; syn-fcc b 4.7e-9) | Kill rule (pre-registered, unchanged): median total H lower on the largest silicon (si32, both amplitudes) AND on >= 2 non-silicon cells (syn-skew2, syn-fcc both amplitudes, mgo8 b, syn-lowsym6 f) -> S1 survives. But the standard estimator (Coulomb residual energy of one step) is not a bound: shorter ladder trajectories end with true errors above the 5e-9 class on a toy slab (coarse seed) and on si8 (cap seed); neither seed rule is clean on the Mac. dev08 must check clause 2 per seed; any violation rejects the switch. Evidence: tests/reference/perf/ladder-standard-secant/ (raw rows and summaries) |
| LADDER-SECANT | Work order 2026-11-11 S2: secant pairs across the ladder hand-off, `BAANDI_LADDER_SECANT=on|<h>` (on = 4; engine `CutoffLadderOptions.secant_pairs`), default OFF. `PulayMixer` stores seeded difference pairs as a closed block (never differenced against absolute pairs), `seed_differences` API, `PulayState.seed_dx/seed_dr` (checkpoint/restore exact; wire fields omitted when empty), last h coarse differences prolonged with `prolong_density`, aged out first; spectral cut unchanged; no mixer parameter changed (P40A-1 not repeated) | implemented, default off; Mac counts mixed (positive only on the larger silicon cells at the realistic amplitude); dev08 A/B requested with h=4 | Mac, same set, strict fine SCF steps / total H eq, ladder -> secant h=2/4/6: si32 f 11/5580 -> 11/5528, 10/5237, 11/5699; si32 b 9/4286 -> 9/4472, 9/4431, 9/4431; si16 f 11/2686 -> 10/2486, 11/2649, 11/2689; si16 b 9/2023 -> 9/2076, 9/2138, 9/2138; si8 f 10/1176 -> 9/1138, 10/1175, 10/1175; si8 b 10/1188 -> 9/1116, 9/1106, 9/1132; si-skew b 453 -> 390/390/416, f 389 -> 468/390/480; mgo8 b 11/1905 -> 12/2056, 14/2324, 13/2152, f 11/2037 -> 14/2456, 14/2470, 16/2734; syn-fcc b 16/531 -> 19/738, 18/647, 19/661; syn-lowsym6, syn-skew2, syn-slab, si-low: equal or up to +2 steps. Standard (with S1): si32 f 2462 -> 2193/2339/2328, b 1761 -> 2569/1752/2246; mixed elsewhere. Strict |dE| vs no-ladder max 3.1e-11 Ha (si16), all certified; standard max |E - E_strict| 5.0e-9 (syn-lowsym6 f h=6) | The prolonged secants help where the dielectric (low-G) response dominates the mixing stall (larger Si cells at the realistic amplitude: si32 strict -6 percent H at h=4, si16 -7 percent at h=2) and hurt on hard/short-period cells (MgO PBE +17 to +34 percent, syn-fcc), where the coarse map's differences misrepresent the fine response (E1 = 0.5 E2 is far below the O dataset's needs). Expected gain 1 to 3 steps (work order E) not seen on the Mac: at most -1 step. No h is uniformly best. A default needs dev08 wall on si16/si32 and a generality check; on this evidence it would be a silicon-only result |
| STANDARD-STOP | Offline fit of the standard-class stopping rule to the dev08 table of 60 trajectories (master relay 2026-11-11, /scr/mino/rustdft/m6tab/steps.csv, 582 rows, decoded here, gzip sha256 511fff2bb235d65168f3c7c87beb4e91dbc569157752dc4d68a62b67f54c6449): si8/si16/si32 standard, ten geometries at +-1e-12 Bohr (seeds 801-810) and ten at +-1e-7 (601-610), each with its own strict certified reference | implemented as `BAANDI_STANDARD_STOP=confirm` (default `plain`); adoption waits for the non-silicon generality table | baseline: every trajectory finishes below 5e-9 Ha (7.6e-11 to 1.7e-10), the first step at or below 5e-9 is 8 (si8, si16 1e-12, si32 1e-12: 10 steps), 9 (si32 1e-7) or 10-11 (si16 1e-7); a plain threshold increase stops early: estimator threshold 1e-8 alone makes 21 of 60 trajectories stop above 5e-9 Ha (step 7, true error 6.4e-8 Ha); the estimator under-estimates only at step 7 (si8 10 of 10, worst ratio 0.117; si32 10 of 10, worst 0.245; si16 never); best rule with zero early stops: stop when estimator <= 2 x threshold (1e-8) AND the previous step's estimator <= 20 x threshold AND the estimator contracted by 0.3 or better; saves 22 steps over 60 trajectories (si32 1e-12: 1 each, si32 1e-7: 0 to 2, si16 1e-7: 0 to 1, si8: 0), worst delivered true error 3.7e-9 Ha; the Aitken-type predicted-next-error rule stops early on all 60 | The estimator is 4 to 9 times conservative on most steps and not stable enough to scale; the savings are limited to si32 by the estimator itself (at si8 step 8 the estimator is 3.4e-8 while the true error is 4.3e-9): a better error estimate is the lever for si8/si16; the rule is fitted on silicon standard only, so mgo-40 and a synthetic multi-projector cell must be checked before the default changes |
| STANDARD-STOP-GENERALITY | dev08 step tables for the non-silicon standard cells (relayed 2026-11-11, master; aeb8665, binary cli-mkl-new 6fdb03136c8512c0): mgo-40 (Mg and O, 3x3x3 k, PBE, kerker-pulay 0.3 q0 1.0 as the fixture sets) and syn-p3 (synthetic multi-projector), ten geometries each at +-1e-12 and +-1e-7 Bohr (fractional shift = amplitude / lattice-vector length), each with its own strict certified reference, 360 rows, gzip sha256 cf329eb6cf23a71668e2a5da921d17926291b4c995ad6c6b57e7e8b2cae8c493 | `BAANDI_STANDARD_STOP=confirm` is SAFE but INERT on both cells and stays OFF by default | both cells are amplitude-insensitive (every trajectory has the same SCF count: 13 on mgo-40, 5 on syn-p3; first step at or below 5e-9 is 12 and 3), so the 40 trajectories are 20 distinct conditions; final true errors max 1.5e-10 (mgo-40) and 5.6e-12 (syn-p3), 0 of 40 finish above 5e-9; estimator/true ratio min 2.30, median 3.89, max 15.96 on mgo-40 and min 20.5, median 23.5 on syn-p3, no under-estimating step; replay of the confirm rule: zero early stops and ZERO steps saved in 40 of 40 (it fires at the final step: contraction <= 0.3 is not met earlier on mgo-40) | The saving of the rule is limited to silicon (si32, 1 to 2 steps), the rule is fitted on silicon, and its safety on mgo-40 rests on the 2x clause with a 15 percent margin (estimator ratio as low as 2.30 admits a true error up to 4.3e-9 against 5e-9); a material with a slightly less conservative estimator would breach it. Decision: no default change; a lever that generalises is a better error estimate (energy-difference based) instead of a threshold on the Coulomb-metric estimator; mgo-40 and syn-p3 are insensitive to amplitude so a second condition would need fractional +-0.001 (offered, not requested) |
| MPI-SCALING-BEYOND-8 | dev08 rank-count and socket scan of the M2d path (relayed 2026-11-11, master; aeb8665, mpi binary 87ae00eeb7ab075e, si32 strict perturbed seeds 601-605, median, M2d; CPUs 0-19 socket 0, 20-39 socket 1, 40-79 SMT) | MPI development is frozen at M2d on this node; the ceiling is replicated work plus memory-bandwidth contention, not NUMA or communication | wall (np=8 socket 0 = 10.400 s, 131 Davidson, 4362 H in every configuration): np=8 split across sockets 10.392 (-0.08 percent, boundary free), np=10 11.122 (+6.94), np=12 11.480 (+10.38), np=16 socket 0 16.628 (+59.89), np=16 two sockets 16.705 (+60.62); rank spread 0.0 to 0.2 percent (no imbalance); np=16 vs np=8 rank-0 phases: potential 0.165 to 0.563 s (+241 percent), mixing 0.153 to 0.286 (+87), density 0.236 to 0.477 (+102), total_energy 0.219 to 0.437 (+100), diagonalization 5.514 to 7.101 (+29), collectives max over ranks 0.829 to 2.105 s (allgatherv 0.099 to 0.442, alltoall 0.256 to 0.569, tree_sum 0.474 to 1.094); earlier same session: np=16 two sockets vs np=8 socket 0 collectives explain 0.54 s of a 6.4 s regression (8 percent); owed M2d items: M2_BATCH=off at np=8 batched 10.325/10.414/10.508 vs unbatched 10.325/10.396/10.377 s (batching buys nothing measurable now, energies bit-identical), si32 STANDARD np=1 8.48 s, np=4 7.25, np=8 6.47 (1.31x against 1.76x strict), determinism at np=4 bit-identical (E -126.70399809612987, density SHA equal on rank 0 and 1), M2 stages at np=8: op.apply_fft_transposes 38.7 percent of the 1.744 s M2 layer (17 percent of the wall), nonlocal 21.3, density 11.9, band_densities 6.6, mirror exchange 5.8, pack/unpack/kinetic 5.5, xc 5.1; threads references: T=8 socket 0 8.531 s, T=16 two sockets 8.781 s (the two T=16 rows of the rank-scan script were invalid: scf.threads hard-coded 1, discarded) | The replicated stages (potential, mixing, density, energy: every rank does the full work) roughly double per rank when the rank count doubles because all ranks share the memory bandwidth of the node; the distributed stages (Davidson +29 percent) suffer the same; np=8 is the best MPI configuration and the threads-only path (8.5 s) still beats it; the first-reading hypotheses NUMA and message cost were refuted by the two controlled comparisons; only a multi-node run (separate memory bandwidth per rank group) can show the benefit MPI exists for, and dev08 is one node; remaining MPI levers if multi-node appears: distribute the Hartree/mixing/energy stages (plane slabs), a smaller replicated footprint per rank |
| LADDER-S2-SECANT | dev08 A/B of BAANDI_LADDER_SECANT off->on (h=4), relayed 2026-11-12 (master; 49abc57, binary cli-s1 efb0d7ad5335c743), si16/si32, seeds 601-610, T=1 and T=8, paired. Strict: regression in 8 of 8 cells (+6.6 to +23.3 %, six with p <= 0.021; si32 b1e-7 T=1 17.618 -> 21.728 s); secant was used (seeded median 4, fallback 0) and still lost; clause-1 violation in 3 si32 b1e-7 seeds (601, 606, 610; dE 1.08e-10, 1.55e-10, 1.05e-10 Ha vs 1e-10), reproduced at T=1 and T=8. Standard (S1 on), T=1: si32 b1e-7 +17.4 %, si16 b1e-7 +11.8 %, f1e-3 cells unresolved. | REJECTED: slower everywhere and breaks the strict class; stays off. Lesson: secant extrapolation of the ladder seed does not return its cost in SCF steps. |
| LADDER-S1-STANDARD | dev08 A/B of BAANDI_CUTOFF_LADDER_STANDARD off->on, relayed 2026-11-12 (master; 49abc57), si16/si32 standard 160 pairs: all within 5e-9 Ha (8.2e-11 to 3.9e-9), density SHA equal T=1 vs T=8. Wall: f1e-3 -6 to -16 % (si32 T=1 -15.8 %, T=8 -10.4 %); b1e-7 T=8 +2.9 % (si32) and +6.9 % (si16) because only SCF steps drop (11 -> 7) while total H work is equal. si8 standard T=1: -12.5 % at both amplitudes, but one pair violates (seed 604, f1e-3, 1.80e-8 Ha vs strict of the same coordinates; Mac saw 6.2e-9 on the same cell). | PENDING: default stays off until the off-arm error of seed 604 and the syn-slab/si8/mgo8 harness are in. |
| STANDARD-RESIDUAL-GUARD | Correctness guard for the standard-class stop, motivated by LADDER-S1-STANDARD si8 seed 604 f1e-3 (stop at Coulomb estimate 4.88e-9 Ha, true error 1.80e-8 Ha; final step requested eigen tolerance 4.93e-10, max occupied residual 3.65e-5 against 1.9e-6 to 1.4e-5 on the 18 other runs). Hypothesis: the density estimator misses the second-order band-energy error of a loose eigensolve. Switch BAANDI_STANDARD_RESIDUAL_GUARD=on (default off) also requires B = sum_k w_k sum_occupied f_i |r_i|^2 (Ha, c = 1/Ha, Kato-Temple second-order bound with the gap scale set to 1 Ha, no fitted constant) below the class threshold; it can only delay a stop | opt-in, PENDING dev08 A/B (si8 standard f1e-3 seeds 601-610, ladder on, guard off vs on) | Mac: unit tests of the arithmetic and of delay-only; off is bit-identical (the switch is read only after the ordinary rule already stops). Limit: residual norms in standard solves may be recurrence norms, and c = 1/Ha is a scale choice, so a small-gap material can still exceed the bound; a large-residual stop without energy error would cost extra steps |
| LADDER-S1-CLAUSE2-HARNESS | dev08 clause-2 harness (relayed 2026-11-12, master; 49abc57; cutoff_ladder.rs bench, 720 rows, no errors; Rust-RNG geometries, differ from the CLI python-RNG ones), standard ladder arm vs strict plain of the same coordinates, 5e-9 Ha, seeds 601-610, T=1 and 8. si8: 0/20 over limit (max 1.92e-9), wall -14 to -33 %, work -29 to -33 %. mgo8: 0/20 (max 1.81e-9), wall -6 to -20 %, work -4 % (b1e-7) / +7.8 % (f1e-3). syn-slab b1e-7: 0/10 but wall +22/+38 % and work +17 %; syn-slab f1e-3: 3/10 over (602 9.89e-9, 603 1.18e-8, 606 1.32e-8 Ha), wall +17/+24 %, work +9.6 %. All three violators stop at fine_scf = 7, the other seven at 8 (plain 11-12 steps, safe; secant 9 steps, safe). fine_scf = 7 is not a general predictor (mgo8 b1e-7 all 7, no violation). Strict clause 1 passes (max dE 4e-15 / 6.7e-13 / 6.1e-12). | S1 stays off by default; combined with the CLI si8 seed 604 case (1.8e-8) the failure is one-step-early stopping of the standard estimator under the ladder. STANDARD_RESIDUAL_GUARD (c82de94) is the candidate fix; syn-slab f1e-3 is the discriminating cell. Lesson: the coarse stage costs more than it saves on tiny cells (syn-slab). |
| STANDARD-RESIDUAL-GUARD-DEV08 | dev08 A/B of BAANDI_STANDARD_RESIDUAL_GUARD off/on with S1 on (relayed 2026-11-12/13, master; c82de94, cli-rg 27fd92c8e7035d60). si8 f1e-3 seeds 601-610 T=1: seed 604 1.801e-8 -> 6.87e-10 Ha (guard adds 5 steps), other nine seeds identical to the last digit, wall median +0.06 %, worst 2.86e-9. si16/si32: guard never fires, wall -0.02 / -0.21 % (noise), worst 4.0e-9 / 1.74e-9. syn-slab f1e-3: violations 3/10 -> 0/10 (max 1.32e-8 -> 2.87e-9), fine_scf {7:3, 8:7} -> {8:10}, extra work only in the 3 firing seeds (+13-14 H each). mgo8 f1e-3: off already 0/10 violations; guard fires in 5/10 seeds (+2 SCF, +216 H), wall +9.9 % (T=1) / +8.3 % (T=8), worst dE 1.64e-9 -> 1.83e-10. | ADOPTED: S1 (CUTOFF_LADDER_STANDARD) and the guard default ON together (standard class only). Cost: mgo8 over-fires (about +9 % wall there, no accuracy needed); candidate refinement: a smaller c or a gap-aware c (see mgo8 note). Opt out with BAANDI_STANDARD_RESIDUAL_GUARD=off / BAANDI_CUTOFF_LADDER_STANDARD=off. |
| MIXING-BETA-STAGE1 | dev08 beta 0.3 vs 0.5 (relayed 2026-11-13, master), five geometries per cell: si8 strict -11.7 % for 0.5 (5/5), si16 strict +5.7 % (0/5), si8 standard -18.8 % (5/5), si16 standard -0.9 % (unresolved). No nonconvergence, |dE| max 2.1e-12 strict / 2.1e-10 standard. | The class rule (strict single-k 0.3, else 0.5) is not supported: direction depends on the cell, not the class. Stage 2 (si32, si-k222s, mgo-40, ten geometries) pending; default unchanged until then. |
| P48-QUEUE | P48 work order adopted 2026-11-13 (user approved order and both decisions). Order: A (scaling profile on current main, not 208b226) -> F, C, D (MKL) in parallel -> E (after v4 verdict; stop predicate and residual guard unchanged) -> B (pure-Rust GEMM, last: dev08 and the QE comparison use the MKL build). Also tracked: si64 QE reference run, ENERGY_IDENTITY dev08 A/B then default, mixing beta stage 2 (si32, si-k222s, mgo-40), M3/M2/M4/M5, displacement-amplitude scan then QE table update, guard over-firing on mgo8 (+9 % wall) refinement, standard estimator recalibration to QE accuracy. | OPEN | Nothing is dropped except rejected rows (S2 secant, STANDARD_STOP=confirm, direct minimisation). |
| MIXING-BETA-STAGE2 | dev08 beta 0.3 -> 0.5, ten geometries (seeds 801-810, +-1e-12 Bohr), T=1, cli-s1 efb0d7ad5335c743 (relayed 2026-11-13, master). Wall change for 0.5: si8 strict -6.6 % (8/10, p 0.109), si8 std -19.0 %, si16 strict +4.7 % (1/10), si16 std -1.0 %, si32 strict +12.4 % (0/10), si32 std -3.5 %, si-k222s strict -5.5 %, std -16.6 %, mgo-40 strict -11.8 %, std -10.9 % (kerker-pulay, beta only). No nonconvergence, no mixer reset, max |dE| 2.1e-10. Mechanism: where 0.3 wins (si16/si32 strict) SCF counts are equal (8 -> 8) and 0.5 only adds Davidson work (+8 to +10 % H); where 0.5 wins the SCF count drops. | Default kept (strict single-k 0.3, else 0.5): it is right in 9 of 10 cells; si8 strict loses 7-12 % of 0.9 s. No atom-count threshold from two points. Open: si16 with 2x2x2 k and si64 strict (queued) to tell atom count from Gamma-only. |
| P48C-NONLOCAL-BLOCK | Complex (k-point) nonlocal of a band block by two GEMMs per fixed 32-band panel (Pk = B^H X, image B (D Pk); MKL zgemm / faer sequential), packed B per k built at construction, switch `BAANDI_NONLOCAL_BLOCK` (default off); branch on 3fea9e8, uncommitted | implemented, opt-in, dev08 pending | Mac O: switch off bit-identical (bench_scf 5 cases x T=1,8); on: dE <= 2.7e-15 Ha, counts identical, density T=1 = T=8, mgo-40 energy identical to the last digit, robustness suite statuses/counts identical, QE k-point comparisons PASS; kernel 35-50 vs 16 GFLOP/s (2.3-3.1x); SCF mgo-40 T=1 -3.8 %, T=8 -2 % (2 pairs) | New premise vs A43/A58b: BLAS-1 replaced by BLAS-3 (B read once per panel at GEMM speed), not a traversal change. Memory: one extra copy of the projector columns per k (n p 16 B). dev08 decides (k444s, mgo-40, si-k222s) | OL P48C |
| P48A-SCALING | dev08 post-adoption scaling profile (relayed 2026-11-13, master; main 3fea9e8, mkl and pure builds, seeds 601-610, +-1e-7 Bohr, bench_suite beta removed so class defaults apply). Wall medians T=1 (mkl/pure): standard si8 0.459/0.414, si16 1.964/1.991, si32 8.168/9.457, si64 38.935/52.551 s; strict si8 0.939/1.021, si16 3.302/3.803, si32 15.140/20.780, si64 96.969/233.380 s. log2 ratio si32->si64: mkl std 2.25, mkl strict 2.68, pure std 2.47, pure strict 3.49. Accuracy contract: T=1 vs T=8 density SHA equal 16/16 combos 10/10; pure vs mkl max 1.2e-12 Ha; standard vs strict max 3.5e-9 Ha, none above 5e-9. si64 strict T=1 stages (pure/mkl): rotate 119.2/20.6 s, orthonormalize 34.6/18.3, h_apply 31.3/28.1, reduced 11.4/4.0, eigen 12.1/5.6. Width histogram si64 strict: 194 steps, 92 with width 1 (47 %), 62 with width 65+. T=8: rotate speeds up 7.76x in pure, orthonormalize only 1.75-2.20x, dense eigen 0.97-1.06x (serial). | ADOPTED as baseline. Targets: B (pure rotate, T=1), D (orthonormalize, T=8), E (width-1 steps), F (see P48A-OTHER). Cubic-fit coefficients are ill-conditioned (negative c2): quote the log ratios only. |
| P48A-OTHER | The "other" time (total minus summed phases, 21-37 % of the T=8 wall) is fully explained (relayed 2026-11-13, master, saved results.json seeds 601-603): other = setup_s + ladder coarse (setup+stage) + certificate solves + result writing (below measurement limit). Standard residual -2.3..+1.8 %; strict residual follows the certificate-solve count (about 0.18 s/0.052 s per solve at si32 T=1/T=8, 1.24/0.46 s at si64). Coarse stage thread speed-up only 1.44-2.33 against fine diagonalization 2.03-3.81; coarse is about 30 % of the T=8 wall in standard (mkl si32 30.6 %, si64 29.5 %), 17 % in strict si64. | F's first target: thread efficiency of the ladder coarse stage and of certificate solves; also report both inside phases. May also remove the S1 T=8 regression at tiny amplitudes. |
| QE-SI64-REFERENCE | QE 7.6 K_POINTS gamma, orig-mpi-dfti, unperturbed, same cell/cutoff/grid/bands (relayed 2026-11-13, master). si32: conv_thr 1e-8 np=1 5.00 s (10 it), np=8 1.70 s; 1e-10 np=1 6.59 s, np=8 1.97 s. si64: 1e-8 np=1 21.39 s (11 it), np=8 4.32 s; 1e-10 np=1 27.68 s; np=8 pending. BAANDI mkl / QE at T=1 vs np=1: si32 std 1.63x, strict 2.30x; si64 std 1.82x, strict 3.50x; at 8 cores: si32 std 2.80x, strict 3.82x, si64 std 4.69x. Growth si32->si64 (log2): QE std 2.10, strict 2.07; BAANDI std 2.25, strict 2.68. Energy offset per atom 1.27e-6 Ha at both sizes. | Standard has nearly caught QE's growth rate (S1); strict is the open gap. |
| P48D-ORTHO-BLOCK | Panel Gram-Schmidt of the real Davidson correction block (`Projection::run_into`): per panel two GEMM pairs against all earlier kept block columns, then the unchanged two-sweep scalar MGS and `DEPENDENCE_TOLERANCE` inside the panel; panel width from a shape rule (largest multiple of 4 with n b doubles in half of a fixed 1 MiB L2, at least 4: si8 20, si16 8, si32/si64 4); switch `BAANDI_ORTHO_BLOCK` (default off, overrides `BAANDI_MGS`); DGKS untouched; branch on 92b26ae, uncommitted | implemented, opt-in, dev08 pending | Mac O (pure, shared host): unit tests: kept-column counts identical classic/MGS-B/block on random, near-dependent 1e-9/3e-12/1e-12/1e-13 and exactly dependent blocks with and without basis, DGKS off and 0.1; defect <= 1.55e-15; bits identical T=1,2,4,8. SCF strict: small cells bit-identical (b >= block width); si16 dE 3.6e-14 Ha, si32 8.3e-13, si64 1.2e-12 Ha; SCF counts equal; H/Davidson equal at si16/si32, si64 6747/176 -> 6751/178; density SHA identical T=1,2,4,8 (si16, si32). Stage orthonormalize si32 T=8 0.67-0.76 -> 0.40-0.49 s (3 rounds), T=1 1.30-1.96 -> 1.23-1.34 s; si64 T=8 7.4-7.6 -> 4.0-4.6 s; `BAANDI_MGS=blocked` (b 16) on the Mac matches ON at si32 | On MKL this is MGS-B (adopted, b 16) with a different width (8/4) and the wide dgemm gate; on pure it parallelizes the serial in-block MGS through the existing (nano-tile) GEMM, inherits P48B. No Cholesky (A52), no change of passes (A37); optimization log "P48D"; dev08 request in the P48D report |
| ENERGY-IDENTITY-DEV08 | dev08 BAANDI_ENERGY_IDENTITY off vs on (relayed 2026-11-13, master; main 3fea9e8, cli-main-mkl 3e6a9ea476d77548, T=1 core 0, seeds 601-610, +-1e-7 Bohr, paired). Wall for on: si16 std +0.50 %, si16 strict +0.41 %, si32 std +1.23 % (0/10, p 0.002), si32 strict +0.33 %. Density SHA identical 40/40 pairs, max |dE| 4.3e-14 Ha, SCF steps unchanged in every pair, no nonconvergence. | NOT ADOPTED (stays off): no speed gain, a 0.3-1.2 % cost; accuracy is unaffected. No user decision needed any more. |
| QE-SI64-CORRECTION | Correction of QE-SI64-REFERENCE (master, 2026-11-13): all eight QE runs completed. si64 conv_thr 1e-10 np=8: 5.81 s, 14 SCF, converged, E -253.379764520 Ha. BAANDI mkl / QE at 8 cores: si64 strict 43.309 / 5.81 = 7.45x (single-core gap 3.50x times thread-scaling gap 2.13x: BAANDI 2.24x vs QE 4.76x from T=1 to T=8). QE np=1 to np=8 scaling: si32 std 2.94, si32 strict 3.34, si64 std 4.95, si64 strict 4.76; BAANDI si64 std 1.92, strict 2.24. | Two axes: single-core kernels (P48B/D/E) and thread scaling (P48F). |
| P48C-NONLOCAL-BLOCK-DEV08 | dev08 BAANDI_NONLOCAL_BLOCK off -> on (relayed 2026-11-13, master; main 92b26ae, cli-p48c 37070d8b16041af92d4a0cb06d19570a494a42c9a7ffaaec1ff7182aa9f919ee, seeds 601-610, paired). Wall T=1: k444s std -1.98 %, strict -1.79 %; mgo-40 std -4.16 %, strict -5.43 %; si-k222s std -1.89 %, strict -1.11 %. T=8: mgo-40 -0.95 % / -1.30 % (10/10), k444s std +0.84 %, strict -0.30 %, si-k222s std +1.53 %, strict +1.30 % (0/10, p 0.002; wall 0.27 s). SCF and Davidson counts identical in all 12 combinations, max |dE| 5.7e-14 Ha, density SHA equal T=1 vs T=8 in 24/24, Gamma si8 bit-identical. | ADOPT with a shape gate: the block path pays where the projector matrix is large (mgo-40) and loses 1-1.5 % on tiny cells at T=8 (fixed GEMM cost). Gate on packed-projector bytes (cache rule), default ON above the gate; below it the old per-band path (bit-identical to off). **Final rule (2026-11-13, uncommitted branch):** `BAANDI_NONLOCAL_BLOCK=auto` (new default; `on` always, `off` never): block path per k when n_pw x n_projectors x 16 B >= 262,144 B (a quarter of a fixed 1 MiB L2, P48D-style hardware-shape constant, shapes only, thread-independent). si2 30 Ry (si-k222s, k444s): 749 x 8 x 16 = 95,872 B -> per-band path; mgo-40: 1495 x 21 x 16 = 502,320 B -> block. Unit tests: below-gate auto bit-identical to off (T=1,2,4,8), gate decision from shapes only. dev08 re-check auto vs off pending. |
| MIXING-BETA-BOUNDARY | dev08 beta 0.3 vs 0.5 boundary cells (relayed 2026-11-13, master): si16 with 2x2x2 k strict -6.55 % for 0.5 (5 geometries, SCF 13 -> 12); si64 Gamma strict +15.91 % (3 geometries, SCF 8 -> 8, H +12 %). Gamma strict by atom count, effect of 0.5: si8 -6.6 %, si16 +4.7 %, si32 +12.4 %, si64 +15.9 %. | Atom count is the main axis inside Gamma-only strict, k points remove the trend. Default unchanged (strict single-k 0.3, else 0.5): right for si16/si32/si64 Gamma, wrong only for si8 strict (0.06 s). |
| P48B-PACKED-GEMM | Real pure-Rust GEMM sites (Gram, reduced, rotation pair, Ritz/residual, column groups, real Gamma nonlocal) through a packed kernel per shape-only panel: x86_64 faer `matmul` Par::Seq (faer's packed AVX-512 GEMM); aarch64 own 8x4 NEON register tile via pulp safe wrappers, one ascending FMA chain per entry; switch `BAANDI_PURE_GEMM=nano|packed` (default nano); branch f70deda + fixes, uncommitted | implemented, opt-in, dev08 pending | Mac O: switch off bit-identical to 92b26ae (gamma, si8, T=1,2,4,8); on: gamma 8.9e-16, si8 5.5e-13 Ha (SCF 9 -> 10), si16 4e-15, si32 1.17e-12 Ha (counts equal), density identical T=1..8; robustness suite statuses identical, unique-solution max dE 5.4e-13, syn-s2 other fixed point (5.0e-3, known multiplicity); 0 steady-state allocations (counting allocator, T=1,4); microbench 41-49 vs 21-30 GFLOP/s nano (faer library 47-59); SCF si32-gamma T=1 -7 %, T=8 -5 % (3 alternating runs, loaded host) | New premise vs ROT1/A49: real SIMD intrinsics / packed library GEMM instead of auto-vectorised loops or finer nano tiles. faer's general GEMM allocates on aarch64 (6520 per measured phase), hence the own kernel there. First version panicked on one-row transposed views (shape-based contiguity test); robustness suite caught it (syn-s1, syn-fft-primes, syn-cubic-degenerate), fixed by stride tests and regression shapes | Contiguity must be tested on strides, never inferred from shape; run the robustness suite before any SCF timing. dev08 decides: x86 path never run (allocation, T identity, GFLOP/s vs MKL), si32 strict trajectory (ROT1 lost there by trajectory change), seeds 601-610 | OL P48B |
| P48F-PHASE-CLOCK | Diagnostics only (no numerical change): results.json `timings` gains `certificate` (blocks, wall, screen, phases inside, outside phases), `atomic_start_s`, `steps_s`, `ladder_coarse` (setup, session, prolong, total, coarse phases, coarse atomic start, steps, certificate) and `closure` (total = setup + phases + certificate outside phases + step overhead + ladder coarse + remainder). Mac O si32 standard T=8: setup 0.079, phases 2.119, step overhead 0.039, ladder coarse 0.649, remainder 0.025 s (sum = total 2.911 s exactly); si32 strict T=8: certificate 2 blocks, wall 0.871 s of which 0.774 s inside phases, 0.095 s screen. | ADOPTED (diagnostic, always on) | Mac si64 standard (O): coarse stage 6.69 s T=1 / 3.34 s T=8; its atomic start (span + two-pass MGS + H Q + Ritz, charged to the first coarse potential phase) about 1.1 s T=1 / 0.9 s T=8, nearly serial: MGS 0.50/0.52 s, Q^H HQ 0.11/0.11 s, raw columns 0.05 s, H Q 0.34/0.14 s. The fine session never runs the atomic start (it starts from the prolonged state). | OL "P48F" |
| P48F-PHASE-THREADS | `BAANDI_PHASE_THREADS=on` (default off): atomic-start span columns built in parallel (fixed one-column chunks) and the two-pass MGS scheduled right-looking (`linalg::orthonormalize_pooled`: after each group of 8 accepted columns every later column takes its first-pass projections in parallel; second pass sequential). Bit-identical to off by construction (each column receives the serial operation sequence; no reduction split); unit test at pool sizes 1/2/3/8 with zero, dependent and near-dependent columns. Mac O: density SHA identical off/on and across T for si32 standard T=1/2/4/8, si32 strict T=1/8, si16 strict T=1/2/4/8, si16 standard, si8-k222, si-skew, si-slab (si-slab variant 601 does not converge in 100 iterations in both arms, pre-existing). si64 standard T=8 (2 pairs, Mac load about 8): atomic start 0.906 -> 0.687 s, coarse stage 3.197 -> 3.043 s, total within noise. | PENDING dev08 A/B | Bound: at most half of the MGS projections can run in parallel bit-identically (the second pass of column j needs all earlier columns final). The remaining serial atomic-start pieces (Q^H HQ and Q rot GEMMs, about 0.16 s at si64 on the Mac) need column-group GEMMs whose bits differ from the off path (rounding class); not done. The coarse Davidson itself also scales worse than the fine one: on the smaller coarse basis width-1 steps fall below `PARALLEL_WORK` (real.rs, owned by P48B; si32 coarse npw 4189 x 1 x ~70 bands < 1,048,576 MAC) and the serial reduced EVD has the same order as at the fine level. | OL "P48F" |
| PROTOCOL-FRACTIONAL-AMPLITUDE-BUG | Master found (2026-11-13) a bug in the new shared module p48lib.py (written 2026-11-12): the fractional +-0.001 displacement was divided by the lattice length, i.e. f<amp> acted as a Cartesian Bohr displacement. Corrected: b<amp> is Bohr per component, f<amp> a fractional component added directly; QE and BAANDI generators are now shared (p48lib displace()). Source check (master, corrected report): every earlier script keeps its own generator and distinguishes b and f correctly (gen_pert.py b only: C1 ladder, EMPTY_LOOSE, BUFFER_TOL, cl/knife; s1.py and s2.py and rg.py correct b and f: S1, S2, residual guard; mb.py/mb2.py b only +-1e-12 Bohr). The p48lib callers used b amplitudes only (1e-7, 1e-12); p48f.py (the only f user) had not run. | NO delivered measurement was affected; amplitude labels in earlier rows are correct. Lesson: the first summary widened the scope without reading the sources; check the generator code before declaring old rows suspect. |
| QE-STARTINGWFC-CORRECTION | Master (2026-11-13): the QE si32 5.00 s of the first equal-band run came from an input without `startingwfc = 'atomic'` (QE default start, 10 SCF); with atomic start si32 1e-8 np=1 is 2.967 s (7 SCF, 28 h_psi), equal to the 11-11 value 2.993 s. QE PHASE A (si32, atomic, medians of 3): nbnd 68 vs QE default: 1e-8 np=1 2.967 vs 2.798 s, np=8 0.724 vs 0.690; 1e-10 np=1 3.487 vs 3.233, np=8 0.831 vs 0.799 (band-count effect 3.9-7.3 %, so PHASE B ran both). | The 2026-11-13 table QE columns are invalid (QE too slow); BAANDI/QE gaps in that table are underestimates (standard si32 1.26x becomes about 2.1x). Replaced when PHASE B is relayed. |
| P48D-ORTHO-BLOCK-DEV08 | dev08 BAANDI_ORTHO_BLOCK off vs on (and MKL classic arm), b7cb5d5, seeds 601-610, both arms DAVIDSON_STAGES=on (walls not comparable with P48A). pure: on is faster: si64 T=8 -18.8 % (65.5 -> 53.2 s; orthonormalize 16.4 -> 6.5 s), si32 T=8 -14.5 %, si64 T=1 -7.4 %, si16 -0.5/-3.1 %. MKL: no gain: si64 T=1 +4.2 %, si32 T=1 +1.7 %, si64 T=8 -0.6 %, si16 -1.4 %/+0.6 %; MKL classic arm (MGS-B off) is +11 to +32 % worse, so MKL already gets the blocking from MGS-B. |dE| max 3.1e-13 Ha, k444s and mgo-40 bit-identical, density SHA equal T=1 vs T=8 in 15/15 combinations. | ADOPTED for pure builds (default on when not MKL_BACKEND), off on MKL builds (backend-dependent default like MGS). |
| M2-DISPLACEMENT-SCAN | dev08 displacement-amplitude scan, si32 Gamma, BAANDI T=1 and QE np=1, same geometries (p48lib displace), seeds 601-605 (relayed 2026-11-14, master): sensitivity ratio BAANDI/QE strict 29.1 (b1e-7), 0.89, 0.68, 0.58 (f1e-4), 0.54 (f1e-3), 0.82 (f1e-2); standard 50.7, 1.51, 0.57, 0.55, 0.61, 0.77. BAANDI increase at b1e-7 is 23.3 % strict / 30.4 % standard (was 63 % before the adoption round). Ladder on/off: absolute wall lower with the ladder at all amplitudes except b1e-5 (+2.6 %); its percentage growth is larger because it lowers the unperturbed baseline. BAANDI growth is work growth (narrower correction blocks), not per-column cost. No nonconvergence, standard error max 3.879e-9 Ha. | Headline amplitude for the QE table: f1e-3 (realistic), b1e-7 kept as stress column. si32 f1e-3: standard 1.71x, strict 2.71x of QE at one core. Lead for E (width floor): the amplitude growth is mainly narrow blocks. |
| M2C-AND-LADDER-F1E-2 | dev08 (relayed 2026-11-14, master): si32 T=8 / np=8: strict unperturbed 6.334 / 0.845 s, b1e-7 +19.1 % / +1.7 % (sensitivity 11.2), f1e-3 +51.4 % / +152.6 % (0.34); standard unperturbed 3.968 / 0.740, b1e-7 +20.0 % / +1.1 % (18.1), f1e-3 +44.3 % / +141.9 % (0.31). Ladder f1e-2 standard T=1: on 11.894 s (SCF 11) vs off 15.233 s (SCF 22), -21.9 %. Incident: dev08 sat idle 28 minutes (01:53-02:21) because a running chain script was edited in place with sed (bash reads by byte offset), so the m2s step was skipped and later chains waited for M2S_DONE; recovered by launching m2s directly. | Rule for the master session: never edit a running chain script, attach new steps as separate new chains. At realistic amplitudes BAANDI is the more stable code at T=8 as well. |
| M2-STAGE-WIDTH-PASS | dev08 stage timers and width histogram, si32, seeds 601-603, unperturbed vs b1e-5 vs f1e-4 (relayed 2026-11-14, master; DAVIDSON_STAGES=on): width-1 steps are NOT the cause of the displacement cost (standard 0 -> 7-9 calls costing 0.08-0.09 s of a 2.0 s increase; strict 13 -> 12 calls); the increase sits in wide blocks (9-64: 0.671 -> 1.388 s standard, 4.245 -> 6.472 s strict; 65+: 1.625 -> 2.404 s standard) and in eigen (strict +140 %), rotate (+104 %), orthonormalize (+98 %) against h_apply +57 %: more Davidson iterations, not narrower blocks. | The width-floor idea (P48E) is not supported by the displaced-geometry data; it keeps its target only for large unperturbed strict (si64 strict, 47 % width-1 steps). |
| HEADLINE-F1E-3-ALL-CELLS | dev08 BAANDI f1e-3 on all cells (38048b9, cli-f3-mkl), against QE re-measurement at the same coordinates; BAANDI/QE one core: si8 std 1.38x strict 1.73x; si16 1.98x / 2.58x; si32 1.74x / 2.69x; si64 1.74x / 3.62x; T=8 against np=8: 2.09/2.58, 3.65/4.23, 3.31/4.48, 3.73/6.55. Density SHA T=1 equals T=8 everywhere, standard error max 3.747e-9 Ha, no nonconvergence. | Headline of docs/performance/qe-vs-baandi-2026-11-13.md. |
| P48B-PACKED-GEMM-DEV08 | dev08 BAANDI_PURE_GEMM nano vs packed (relayed 2026-11-14, master; 38048b9 pure native, cli-p48b-pure). x86 path first run clean: pure_gemm tests 4/4, steady-state allocations 0, off path bit-identical (gamma, si8), robustness suite same statuses (si-elongated, syn-s3, syn-z14 unconverged in both). Microbench GFLOP/s nano -> packed (equals faer library): rotate si64 9.5 -> 61.3, gram 15.7-16.4 -> 53.7-57.0, reduced 16.4-17.2 -> 62.3-63.0, nonlocal 14.7-16.2 -> 54.2-64.8 (target 60 % of MKL met; narrow gram n=4 7.7 -> 9.9). Wall T=1 packed vs nano: si16 std -15.6 %, strict -20.2 %; si32 std -25.3 %, strict -31.8 %; si64 std -38.0 % (51.07 -> 31.66 s). T=8: -0.2 to -5.1 %. SCF/Davidson/H counts identical, max |dE| 8.5e-13 Ha, density SHA identical at T=1,2,4,8, certificate identical. (Microbench T=8 header printed threads 1: the T=8 microbench numbers are T=1 repeats.) | ADOPTED for pure builds (default packed when not MKL_BACKEND, nano on MKL builds where dgemm replaces these sites). si64 strict and T=8 rest pending but no regression anywhere so far. |
| P48B-FINAL-TABLE | Completed P48B A/B (relayed 2026-11-14, master; 38048b9 pure, seeds 601-610): wall nano -> packed: si16 std T1 1.995 -> 1.684 (-15.6 %), T8 -0.24 %; strict T1 3.762 -> 3.002 (-20.2 %), T8 +0.08 %; si32 std T1 9.320 -> 6.958 (-25.3 %), T8 3.904 -> 3.784 (-3.1 %); strict T1 20.445 -> 13.953 (-31.8 %), T8 7.146 -> 6.779 (-5.1 %); si64 std T1 51.065 -> 31.659 (-38.0 %), T8 16.226 -> 14.171 (-12.7 %); strict T1 231.354 -> 92.447 (-60.0 %, 2.50x), T8 53.665 -> 34.874 (-35.0 %). SCF/Davidson/H counts identical in all 12 combinations, max |dE| 5.97e-13 Ha, certificate values identical, density SHA identical at T=1,2,4,8, no nonconvergence. Microbench T=8 (PURE_GEMM_BENCH_THREADS=8; the earlier "T=8" run was T=1): packed 1.1-4.2x nano (si64 rotate 73.1 -> 308.1 GFLOP/s), equal to or slightly above faer library (si32 reduced 253.5 vs 240.0). | Confirms adoption (pure builds default packed, commit 2dbac32). The gain grows with size, also at T=8 for si64. Pure stack baselines of P48A, P48D, P48B and P48F tables are not comparable (defaults changed twice). |
| P48C-AUTO-RECHECK-DEV08 | dev08 BAANDI_NONLOCAL_BLOCK auto (default) vs off (relayed 2026-11-14, master; 2dbac32 MKL, seeds 601-610): k444s and si-k222s (below the gate) bit-identical (dE 0, density SHA 10/10, wall -0.42..+0.58 %); mgo-40 (above the gate) T=1 std -4.40 %, strict -6.42 %, T=8 -1.05 % / -0.87 %, dE 5.7e-14 Ha, density SHA differs (rounding class, intended); SCF/Davidson/H identical in 24 rows; scf.switches reports NONLOCAL_BLOCK = auto, source default (key BAANDI_NONLOCAL_BLOCK). | CONFIRMED: auto default stays; the gate removes the 1-1.5 % loss on small k-point cells. |
| P48F-PHASE-CLOCK-AND-THREADS-DEV08 | dev08 closure and BAANDI_PHASE_THREADS off/on (relayed 2026-11-14, master; 2dbac32 MKL, seeds 601-610). Closure sums exactly to total_s (si64 strict T=8: setup 0.897 + phases 32.667 + certificate outside phases 0.847 + step overhead 0.097 + ladder coarse 6.821 + remainder 0.354 = 41.683 s). The coarse stage of the ladder spends 43 % (2.944 s of 6.821 s) in the atomic start (si64 strict T=8): 7.1 % of the whole wall. PHASE_THREADS on: bit-identical (dE 0, density SHA 10/10, counts identical); T=8 atomic start 2.861 -> 2.252 s (-21 %) at si64, coarse -11.7 %, wall -2.50 % (std), -1.73 % (std f1e-3); si32 T=8 -1.93 % / -0.63 %; T=1 +0.16..+0.31 % (noise, no thread effect). | ADOPTED (default on): bit-identical, T=8 gains 0.6-2.5 %, no regression beyond noise. Remaining unthreaded atomic start 2.25 s (11.6 % of si64 T=8 wall): second Gram-Schmidt pass and the products Q^H HQ and Q rot. |
| P49A-PROJECT-WIDE | Gap-analysis rank 1 (davidson-algebra-gap-analysis.md 3.1): `Projection::run_into` apply X -= Q S (site ProjectionApply) and the MGS-B cross-panel apply (site Other) through `multiply_slices_wide_into` (fixed 256-row panel dgemm, only the per-panel tile bound) instead of `multiply_slices_into` (whole-output bound rows x cols <= capacity^2: c <= 26 at si64, 14 at si32, 7 at si16, wider blocks fell to the nano pure kernel); switch `BAANDI_PROJECT_WIDE` (default off until dev08; ORTHO_BLOCK already took the wide route for its panel apply); branch on 06589a4, uncommitted | implemented, opt-in, dev08 pending | Mac O (pure build: both routes forward to the same product, so on = off bitwise by construction): unit tests GEMM level n 4096 x k 40 x c 30, capacity 96 (wide/square within 2 k eps sum abs(a)abs(b), bitwise T=1,2,4,8) and projection level (classic and MGS-B, random / near-dependent / exact blocks: kept counts equal, defect <= 1e-12, T=1,2,4,8 bitwise); MKL path type-checked only (clippy --features dft-numerics/mkl); SCF si16/si32 gamma base and seed 601 f1e-3: main / off / on density SHA, energy and iterations identical at T=1 and T=8 | Premise: MKL-R2 meant both projection products on dgemm and P48B kept nano on MKL builds on that premise; DGEMM-SITES (projection_apply 8.72 vs overlap 53.7 GFLOP) shows the gap. dev08 decides: rounding class only on MKL builds |
| WORK-ORDER-BEAT-QE | User work order 2026-11-14 (docs/performance/requests/work-order-beat-qe-2026-11-14.md): BAANDI wall <= 0.90 x the fastest QE configuration per cell and per parallel model (one core, threads vs OpenMP, MPI vs MPI); single standard path with verification (certificate) on by default; GPU out of scope; MPI reopened. Phases: P0 (single path, timers that close, honest QE reference incl. no diago_full_acc, rmm-davidson, ppcg, OMP and MPI 1/2/4/8), B (lean Davidson algebra: B1 never rotate all bands, B2 capacity, B3 generalized problem vs one-pass orthonormal, B4 reduced solve), C (C1 pruned FFT on MKL, C2 no double transform of converged orbitals, C4 PGO/huge pages), D (SCF steps: D0 baseline and QE-equivalent mixer, D1 ladder hand-off), E (threads, MPI). Checkpoints: after B step 4 ratio above 1.25 on si32 and si64 -> re-plan; after D step 6 above 1.0 -> scope decision. Overlap noted: P49A PROJECT_WIDE and P49B B1 narrow/margin (in dev08 queue) are prerequisites/variants of B1; A13/P31A (capacity) and P37A (generalized problem) are reopened by the work order with stated new premises. | OPEN | Decisions recorded in AGENTS.md; mixed-model columns removed from qe-vs-baandi-2026-11-11 and 2026-11-13. |
| P0.2-WORK-TIMERS | Work order beat-QE P0.2, timers that close (instrumentation, no numerical change): under the existing switch `BAANDI_DAVIDSON_STAGES=on` (no new switch), `convergence.work.davidson_stages` gains `fft` (paired orbital transforms of the real Gamma h_apply split into inverse / potential multiply / forward, density transforms, other wave transforms; transform and column counts; `wave_fft_s`, `wave_fft_calls` = 2 x pairs + density + 2 x other, the counterpart of QE's `fftw` wave-transform count), `h_apply_parts` (apply, local, nonlocal projection P^T X, nonlocal back product, remainder) and `by_run_stage` (fine, ladder_coarse, verification = certificate block) with every block; pooled keys are the sums and keep their meaning. Mechanism: one generic body per FFT plan (`apply_local_with<K: Laps>`, `NoLaps` compiles away), clocks per pair (per slab in the mirror-slab plan), per density transform and per block apply, relaxed atomics per run stage (index set at serial boundaries in `cutoff_ladder::coarse_stage` and the session certificate block). Mac si16 (2x1x1 conventional, 30 Ry, f1e-3 seed 601, strict, loaded host 5-10): energy -63.129243717904316 Ha and density sha 005531472e8c26a6 identical for base 30442a9 vs branch (off) and off vs on at T=1 and T=8; wall base vs branch-off 2.667 / 2.639 s (T=1, 12 pairs), 1.146 / 1.129 s (T=8, 8 pairs); off vs on 2.617 / 2.718 s (T=1, +3.9 percent, noise of the loaded host), 1.120 / 1.121 s (T=8). si16 T=1 on: wave FFT 1.05 s of 2.57 s total (pairs 1553 = 3055 columns, density 256 transforms), h_apply 1.17 s = local 0.99 + nonlocal 0.09 + 0.08 + 0.01. | ADOPTED (diagnostic, default off) | Times are thread-seconds (summed over workers), counts are thread-independent; at T>1 divide by busy threads before setting beside wall. QE's fftw timer covers wave transforms only, so compare with `wave_fft_s` minus or including `density_s` as QE's sum_band share requires (QE fftw includes sum_band transforms). Mac numbers are observations; dev08 request in the P0.2 report. |
| P48F-PHASE-THREADS-FINAL | Completed P48F A/B BAANDI threads, BAANDI_PHASE_THREADS off vs on (relayed 2026-11-14, master; 2dbac32, seeds 601-610, T=8): mkl si64 strict f1e-3 49.246 -> 45.375 s (-7.86 %), pure si64 standard b1e-7 -4.05 %, mkl si64 standard b1e-7 -2.50 %, pure si64 strict b1e-7 -2.25 %, mkl si64 strict b1e-7 -2.08 %, mkl si32 standard b1e-7 -1.93 %, mkl si64 standard f1e-3 -1.73 %, pure si64 strict f1e-3 -1.09 %, pure si32 strict -1.05 / -1.01 %, mkl si32 strict f1e-3 -0.63 %: all T=8 combinations improve, no regression; T=1 -0.84..+1.08 % (noise, no gain claimed). Bit-identical (dE 0, density SHA equal, SCF/Davidson/H identical, closure closes). Mechanism: ladder_coarse atomic start -21 % at si64 T=8. | CONFIRMS the default-on adoption (40b6b84 era). Remaining: atomic start 2.25 s (about 11 % of si64 T=8), unthreaded second Gram-Schmidt pass. |
| P0.3-QE-MPI-INVALID | Incident (relayed 2026-11-14, master): all 153 QE MPI rows of the P0.3 sweep were invalid (wall 0.003 s, empty output, recorded as nonconverged) because mpirun was called with `--cpu-set`, which Open MPI 5.0.8 rejects as deprecated, so pw.x never ran; the harness did not distinguish a failed launch from nonconvergence. Fix: `--map-by core --bind-to core` (ranks on cores 0..n-1, socket 0 for n <= 8; verified si16 1e-8 np=2: 11 iterations, 0.54 s); OpenMP rows (207) unaffected and kept; MPI rows are rerun; BAANDI MPI leg started with the corrected command. Interim QE OpenMP results (si8, si16, f1e-3, seeds 601-603 medians): diago_full_acc=false is fastest in nearly all combinations (one si8 point reversed); si16 1e-8 1/2/4/8 threads 1.063/0.958/0.871/0.816 s (1.30x), 1e-10 1.293/1.164/1.047/0.979 s; si8 1e-8 0.324/0.290/0.281/0.297 s (slower again at 8 threads). | Rule for every harness: exit code check, empty output or wall below 0.1 s is INVALID, the chain stops and reports on the first INVALID cell. The sweep order is now BAANDI MPI ranks, QE MPI rerun, si64 QE leg (reduced). |
| BAANDI-MPI-RANK-SWEEP | dev08 BAANDI MPI (M2d, `--features mpi`, ranks 1/2/4/8, one thread per rank, `--map-by core --bind-to core`, mpi build 6f437a92b0f235b7), relayed 2026-11-14, master, si8/si16/si32 (si64 running); wall medians (unperturbed 1 run, f1e-3 seeds 601-603). si32 strict f1e-3: np=1 20.169, np=2 23.839, np=4 15.394, np=8 13.677 s (best 1.47x); si32 strict unperturbed 12.532 / 14.485 / 9.819 / 8.305 (1.51x); si16 strict f1e-3 4.373 / 4.672 / 3.351 / 2.868 (1.52x); si8 standard unperturbed np=1 fastest. np=2 is slower than np=1 in 11 of 12 combinations (up to +18 %). Work counts are rank independent on displaced geometries (si32 strict f1e-3 SCF 11, Davidson 132, H 4656 at all ranks); on the symmetric unperturbed si8/si16 strict the Davidson path varies with ranks. Guard check: 54 rows, 0 FAILED, 0 INVALID (the earlier run died on the rejected `--cpu-set`). | Single-model measurement for the MPI against MPI table (QE MPI rerun pending). BAANDI MPI scales to at most about 1.5x at 8 ranks on these cells: the E2 work (distribution of the replicated stages) is needed before the MPI target is reachable. |
| BAANDI-MPI-RANK-SWEEP-SI64 | dev08 BAANDI MPI ranks 1/2/4/8, si64 (relayed 2026-11-14, master; 70 rows, 0 FAILED, 0 INVALID): wall np=1/2/4/8: standard unperturbed 27.290 / 41.553 / 29.664 / 24.805 s (best 1.10x), standard f1e-3 46.757 / 67.589 / 48.946 / 35.608 (1.31x), strict unperturbed 85.567 / 144.700 / 103.795 / 62.133 (1.38x), strict f1e-3 112.163 / 163.838 / 111.305 / 77.314 (1.45x). np=2 slower than np=1 by +45 % (std f1e-3) and +69 % (strict unperturbed); np=4 still 21 % slower at strict unperturbed. Best 8-rank scaling over cells: si8 1.03-1.19x, si16 1.34-1.52x, si32 1.37-1.47x, si64 1.31-1.45x. Work counts identical across ranks on displaced geometries (si64 strict f1e-3 SCF 11, Davidson 176, H 10068). | The np=2 penalty grows with problem size (si8 +9 %, si16 +7 %, si32 +18 %, si64 +45-69 %): fixed partition/communication cost proportional to the problem, not to the rank count. Lead for the MPI target (work order E2): the np=1 run may not use the same code path as np>=2 (verify), and the replicated stages and collectives need the per-stage table (parallel.stages). |
| MPI-ENTRY-COST-DIAGNOSIS | dev08 stored results.json, si64 strict f1e-3 seed 601 (relayed 2026-11-14, master; mpidiag.py): np=1 has no `parallel` block (serial path, no M2); np>=2 runs the distributed M2 path (`mpi_stick_rows`), so np=1 against np=2 compares the serial path with the distributed path. total_s np=1/2/4/8 = 109.3 / 160.7 / 108.4 / 73.7; diagonalization 80.2 / 128.5 / 81.3 / 43.7; collectives 0 / 3.37 / 4.62 / 4.57 s (6.6 % at np=2: alltoall 1.97 s in 980 calls, tree_sum 1.23 s in 76258 calls, allgatherv 0.17 s); m2 stages np=2/4/8: op.apply_fft_transposes 13.60 / 7.41 / 4.62 s, op.apply_nonlocal 9.85 / 5.74 / 3.46, pack/unpack/kinetic 1.83 / 1.10 / 0.89, band_densities 1.61 / 1.04 / 0.89; density is faster distributed (3.195 -> 2.047 s). Work counts identical (Davidson 97, H 9008, certificate H 1096, 3 solves). Per-rank energies equal rank 0 exactly. | Cause of the np=2 penalty: the distributed apply kernels, not communication. Per-rank apply work at np=2 (25.3 s for transposes + nonlocal + pack) is about the cost of the whole serial h_apply (28 s) although each rank handles half the rows: the distributed FFT and nonlocal kernels are about 1.6-1.8x slower per unit of work than the serial MKL pair-FFT and dgemm path. Collectives do not shrink with ranks (3.4 -> 4.6 s) and become the next ceiling. |
| MPI-ENTRY-COST-FIXES | Analysis and opt-in fixes of the M2 entry cost (docs/performance/mpi-entry-cost-analysis.md; branch on e7cb6aa, uncommitted). Reading of MPI-ENTRY-COST-DIAGNOSIS corrected: the apply excess at si64 np=2 is about 11 s (25.3 s against an ideal 14), the distributed Davidson algebra excess about 70 s (103 s against an ideal about 26 plus the replicated eigenproblem); the 76258 tree_sum calls come from the classic distributed orthonormalization (3 reductions per column, two whole-block sweeps, three basis projections), whose O(nl w^2) scalar work is the cost, not the reduction latency. Root causes (C): classic order does about four times the serial orthonormalization work; distributed rotations/projections and the distributed nonlocal miss dgemm on MKL builds (whole-output bound, `multiply_into`); distributed FFT uses rustfft lanes while the serial MKL path uses the full-grid MKL 3D FFT. Switches (all default off): `BAANDI_M2_ORTHO=serial` (serial order: one projection per pass, DGKS, panel Gram-Schmidt with the serial panel rule), `BAANDI_M2_DGEMM=on` (wide per-panel dgemm route for rotations, projections and the M2 nonlocal), `BAANDI_M2_FFT=mkl` (MKL DFTI batch lanes, `BackendLanes`), `BAANDI_M2_LEAN=on` (slice loops, bit-identical) | implemented, opt-in, dev08 pending | Mac O (pure, simulated ranks, Si.pz-vbc 15 Ha f1e-3, work counts identical in every row, |dE| vs serial <= 8.5e-14 Ha): si32 np=2 18.33 -> 13.74 s (-25 %; ORTHO alone 13.96), np=4 12.89 -> 10.22 (-21 %); si16 np=2 2.48 -> 1.90 (-23 %), np=4 1.97 -> 1.70 (-14 %); orthonormalize stage per rank si32 np=2 6.26 -> 2.10 s; tree_sum calls per rank 39363 -> 23507 (si32), 17055 -> 10398 (si16); on equal kernels the distributed apply is at its ideal (FFT 1.07, nonlocal 1.10 per unit work against serial). Tests: vendor lanes 1e-14 (bit-identical pure), nonlocal dgemm route 1e-14 (bit-identical pure), lean and wide bit-identical (pure), serial order within 1e-9 of serial eigenvalues for three panel rules x DGKS on/off, engine SCF strict |dE| <= 2.2e-16 / standard <= 6.6e-11 Ha, ranks bit-equal, reproducible. Expected at si64 strict (E): np=2 -45 to -65 s of 160.7, np=8 -12 to -17 s of 73.7 | The np=2 penalty is mostly the distributed algebra, which was written for M1 before MGS-B, DGKS and the wide dgemm route existed in the serial solver; bringing the serial solver's structure into the distributed one is the first E2 step. No rejected idea repeated (no Cholesky, A52/MPI-M1-GATE1; M2d replicated stages and mixing untouched). Next: dev08 A/B per switch (request in the analysis doc), then fuse the in-panel norm, then the serial B1/B5/B6 stack in the distributed solver |
| QE-PPCG-REMOVED | Master (2026-11-14): QE 7.6 removed the PPCG diagonalization in December 2024 (`Error in routine iosys (1): PPCG diagonalization not supported anymore (Dec. 2024)`, exit at input parsing). All 42 ppcg rows of the earlier OpenMP sweep and the MPI ppcg rows failed (wall 0.04 s, no convergence); the new run guard (qe_check) stopped the QE MPI rerun on the first one. Only rmm-davidson remains as an alternative QE diagonalization (works: si32 1e-8 unperturbed full_acc=true MPI np=2/4/8 1.710 / 1.002 / 0.796 s, 8 SCF). si64 reduced matrix is now (true, david), (false, david), (false, rmm-davidson): 84 combinations. | The work order's P0.3 ppcg arm is impossible with this QE; the question "is a QE alternative faster than david" is open only for rmm-davidson; if rmm-davidson is faster for a cell and model, that number is the QE target (decision 6). |
| P0.3-QE-FASTEST-PARTIAL | dev08 QE fastest configurations so far (relayed 2026-11-14, master; atomic start, nbnd = BAANDI value, f1e-3 seed 601 unless stated; INVALID guard clean). si8 and si16 (david only), full_acc=false is fastest everywhere: wall omp1/2/4/8 and mpi2/4/8: si8 1e-8 0.322 / 0.289 / 0.277 / 0.294 and 0.266 / 0.222 / 0.214; si8 1e-10 0.406 / 0.364 / 0.346 / 0.371 and 0.314 / 0.297 / 0.242; si16 1e-8 1.071 / 0.963 / 0.871 / 0.816 and 0.646 / 0.454 / 0.397; si16 1e-10 1.293 / 1.164 / 1.047 / 0.979 and 0.765 / 0.566 / 0.462. QE MPI is 2.06x faster than QE OpenMP at 8 cores for si16 1e-8 (0.397 vs 0.816 s); QE OpenMP saturates at 4-8 threads (si8 slower at 8). si32 1e-8 f1e-3 (3 seeds, same pattern): rmm-davidson beats david in 6 of 7 settings: mpi2 2.583 vs 2.954 s (-12.6 %), mpi4 1.555 vs 1.848 (-15.9 %), mpi8 1.215 vs 1.513 (-19.7 %), omp1 4.899 vs 5.447 (-10.1 %), omp2 3.454 vs 3.831, omp4 3.010 vs 3.245, omp8 3.157 vs 3.119 (david +1.2 % faster). | QE target for si32 = rmm-davidson with full_acc=false (omp8 excepted): the 0.90 x QE target for si32 moves down by 7-20 %. Same-model pairs (BAANDI f1e-3 medians): si16 threads 8 vs QE OpenMP 8: 1.841 / 0.816 = 2.26x; MPI 8 vs MPI 8: 1.973 / 0.397 = 4.97x; si8 threads 1.94x, MPI 2.93x. The MPI axis is the furthest from the 0.90 target. |
| P0.3-SI32-CORRECTED | Correction of P0.3-QE-FASTEST-PARTIAL (master, 2026-11-14): the david vs rmm-davidson comparison mixed diago_full_acc true and false in its medians. With full_acc fixed to false (si32, 39 comparisons): rmm wins 30, david 9, and all david wins are on the unperturbed geometry. Displaced f1e-3: 1e-8 s601 rmm vs david mpi2 -22.9 %, mpi4 -26.8 %, mpi8 -29.8 %, omp1 -19.8 %, omp8 -10.3 %; 1e-10 s602 mpi2 -31.0 %, mpi4 -36.6 %, mpi8 -41.8 %. Unperturbed: david wins at 1e-8 (mpi2 +8.8 %, mpi8 +5.5 %, omp8 +24.5 % for rmm), mixed at 1e-10. si32 overall fastest QE (all full_acc=false; david and rmm both): 1e-8 s601 mpi2 2.497 david, mpi4 1.517 rmm, mpi8 1.169 rmm, omp1 4.715 david, omp8 2.702 david; 1e-8 unperturbed mpi2 1.521 rmm, mpi4 0.911 rmm, mpi8 0.696 rmm, omp1 2.834 rmm, omp8 1.675 david; 1e-10 s601 mpi2 2.845, mpi4 1.723, mpi8 1.343 all david; 1e-10 unperturbed mpi2 1.811 david, mpi4 1.081 rmm, mpi8 0.825 rmm; si32 OpenMP 1e-10 still missing (64 combinations queued after the si64 leg). | The QE target is the per-(cell, threshold, geometry, model) minimum over david and rmm-davidson with full_acc=false (david is much faster with full_acc off, rmm gains less). The earlier "rmm is 7-20 % faster" statement is retracted: always fix full_acc before comparing diagonalizations. |
| M2-ENTRY-COST-AB-DEV08 | dev08 A/B of the four MPI entry-cost switches (relayed 2026-11-15, master; a81b69d, cli-m2ab 1887b267..., BAANDI MPI, one thread per rank, `--map-by core --bind-to core`, f1e-3 seed 601 (si32: 601-603 medians), 0 INVALID in 65 rows so far; table headline "BAANDI MPI, switch arms"). si64 standard np=2: arm1 (off) 68.540 s vs serial path 46.923 s; M2_ORTHO=serial 47.874 (-30.2 %), M2_DGEMM=on 53.640 (-21.7 %), M2_FFT=mkl 67.518 (-1.5 %), M2_LEAN=on 66.475 (-3.0 %), all four 34.689 s (-49.4 %, 26 % faster than the serial path); np=4 arm2 -19.9 %. si64 strict np=2: off 158.981 vs serial 110.216 s; ORTHO -30.6 %, DGEMM -25.1 %, FFT +0.85 %, LEAN -0.16 %, all four 73.614 s (-53.7 %, 33 % faster than serial); stages ortho 60.9 -> 10.6 s (-83 %), rotate 30.5 -> 9.2 (-70 %), apply 26.4 -> 19.3, collectives 3.65 -> 2.58. si32 standard np=2/4/8 all-four -31.5 / -23.3 %/..., si32 strict np=4 -31.3 %. Work counts identical (si64 strict SCF 11, H 10104; one Davidson count 161 vs 162 from rounding), max |dE| 1.99e-13 Ha, every rank energy equals rank 0; M2_LEAN always bit-identical; reduced eigen unchanged by the switches (1.0-1.4 s at si32, 5.5 s at si64 strict, replicated, grows with ranks). MKL rounding tests (--features dft-numerics/mkl) vendor_lanes_match_rust_lanes, nonlocal_dgemm_route_matches_the_default_route, mkl_lanes_match_rust_lanes all pass. | ADOPTED: BAANDI_M2_ORTHO=serial, BAANDI_M2_DGEMM=on, BAANDI_M2_LEAN=on default (2026-11-15); M2_FFT stays rust (neutral). The np=2 entry penalty is gone. Next ceiling: the replicated dense reduced eigenproblem (5.5 s at si64 strict, rising with ranks) and collectives. |
| MPI-REDUCED-EIGEN | Analysis of the replicated dense reduced eigenproblem of the distributed Davidson (docs/performance/mpi-reduced-eigen-analysis.md; branch on 87d3b7b, uncommitted). Starting point (R): reduced_eigen si32 strict np=8 1.37 of 10.0 s, si64 strict np=2 5.45 of 73.6 s, np=8 7.20 of 47.7 s, unchanged by the entry-cost switches. Counts (C, Mac analogues, exact): one all-pairs dense solve per Davidson iteration on every rank, no B5/B1; si32 293 solves, order mean 224.6, p90 393, max 408, sum m^3 6.5e9, single-column appends 25 %; si64 228 solves, mean 460.8, p90 743, max 787, sum m^3 3.87e10, single-column 18 %; tree_sum 23507 (si32) and 46722 (si64) calls per rank, about 85 % per-column in-panel reductions, 14 % warm-start orthonormalization; no material exact fusion (MGS dependence chain). Switches (default off): `BAANDI_M2_EIGEN=root` (rank 0 solves, one broadcast of status + lowest pairs; bit-identical), `BAANDI_M2_B5=on` (serial B5 bordered update on single-column iterations; rounding class); results.json `parallel.reduced_solves` | implemented, opt-in, dev08 pending | Mac O (pure, simulated ranks, single runs on a loaded Mac, work counts identical in every row): reduced_eigen si32 np=4 1.400 -> 1.194 s (root), si32 np=2 1.202 -> 1.130 (B5), si64 np=4 6.737 -> 5.481 (root + B5); B5 fired on 70/74 (si32) and 42/42 (si64) single-column solves. Capacity per band (Mac O si32 np=2): 4 -> Davidson +27 %, H +14.8 %, eigen -39 %, other Davidson stages +16 %, wall +10 %; 3 -> H +34 %, SCF +1: net loss, not proposed. Tests: root bit-identical to replicated (solver P=2-5, engine SCF energy and density bits); B5 within 1e-9 eigenvalues, SCF |dE| <= 4.4e-16 Ha strict, <= 6.6e-11 Ha standard; ranks bit-equal; reproducible. Expected si64 strict (E): np=8 root 0 to -1.7 s, B5 -0.2 to -0.4 s | The stage is replicated single-thread work that grows with ranks (same work 32 % slower per rank at np=8: concurrent copies share L3 and memory bandwidth, as in MPI-SCALING-BEYOND-8); a root solve removes the contention but not the critical path. Only a distributed dense EVD (ScaLAPACK pdsyevd / ELPA) shortens the path: next a dev08 microbenchmark. Collectives: exact gain only from recursive-doubling tree sum (bit-identical for power-of-two P, about half the latency); count reductions (Pythagorean norm, skip orthonormal warm start) are rounding class. Not repeated: dsyevr (MKL-EVD), threaded dsyevd, CholeskyQR (A52, MPI-M1-GATE1), smaller capacity (A13, confirmed on the MPI path) |
| M2-AB-SI64-STRICT-NP8-CORRECTION | Correction of M2-ENTRY-COST-AB-DEV08 (master, 2026-11-15): the si64 strict np=8 row of the A/B is withdrawn: its arm-1 baseline was 76.769 s while three repeats now give 61.246-61.389 s (median 61.276 s), a node-state artefact (load 2.38 right after the np=4 arms), so every arm of that row was compared with a wrong baseline. Rerun, 3 repeats per arm: arm 1 61.276 s, arm 4 (FFT) 62.132 (+1.40 %), arm 5 (LEAN) 61.732 (+0.74 %, dE 0, so the former -20.4 % was noise), arm 6 (all) 47.443 s (-22.58 %, spread 0.25 s, dE 8.5e-14). Arms 2 and 3 not remeasured (extrapolation from the old numbers would be about -16 % and -1 %, an estimate, not a measurement). Reduced eigen unchanged by every arm (7.18-7.35 s, 15.2 % of arm 6). Gain of all switches shrinks with ranks: si64 strict np=2 -53.70 %, np=4 -41.55 %, np=8 -22.58 %. Other nine combinations unaffected. | The adoption decision stands (ORTHO serial and DGEMM effective, FFT none, LEAN neutral and bit-identical). Rule: when a baseline row looks off, rerun the baseline three times before using any arm of that row; single runs with load above 1 are not baselines. |
| HARNESS-SHADOW-INCIDENT | Incident (relayed 2026-11-15, master): a harness script named `re.py` shadowed Python's `re` module (statistics -> fractions -> re.compile circular import), so three jobs died at their first line (arm 2/3 rerun, reduced-eigen A/B, P0.1); the chains printed their completion tokens without checking the return code and released the next steps over empty outputs. About 13 minutes lost (05:28-05:41), no measurement contaminated; the jobs were rerun (redochain 05:41). si32 OpenMP remainder completed normally (65 rows, 0 INVALID). | Rules for every dev08 chain: (1) never name a script like a standard library module (checked all names in the directory), run harness Python with `python3 -I`; (2) every chain step checks the return code and stops on non-zero; (3) a chain must not print a DONE token unless the step produced its expected output file with rows. |
| M2-AB-SI64-STRICT-NP8-ARMS-REMEASURED | Same-session remeasure (relayed 2026-11-15, master; a81b69d, cli-m2ab, si64 strict f1e-3 seed 601, np=8, 3 repeats per arm, load 0.02, 0 INVALID): arm1 61.449 / 61.149 / 61.858 s (spread 0.709); arm2 M2_ORTHO=serial 49.970 s (-18.68 %), arm3 M2_DGEMM=on 59.492 s (-3.18 %), arm4 M2_FFT=mkl +1.40 %, arm5 M2_LEAN +0.74 %, arm6 all -22.58 %; dE 0 for arms 2 and 3. Withdrawn single-run figures (-33.15 / -20.95 / -17.46 / -20.38 / -37.85 %) came from an outlier baseline row (76.769 s); the arm wall times agree between sessions within 1 s, only the baseline was off. np=2 (-53.70 %) and np=4 (-41.55 %) rows stay valid. Gains are nearly independent (arm2 + arm3 -21.9 % vs arm6 -22.6 %); reduced_eigen unchanged (7.20-7.26 s, 15.2 % of arm 6). Interim reduced-eigen A/B (si32 standard np=2/4): M2_EIGEN=root +0.06 %, B5 +1.35 %, both +0.52 % (noise, dE 0 or 2.8e-14); si32 np=8 baseline row flagged BASELINE UNSTABLE (load 1.27) and to be remeasured. | Adoption of ORTHO serial / DGEMM / LEAN unchanged; the np=8 gain at si64 strict is -22.6 % not -37.9 %. The reduced-eigen switches show no effect at si32 standard (order too small); si64 decides. |
| MPI-COARSE-CERT | Non-scaling coarse stage and certificate of the M2 path (docs/performance/mpi-coarse-certificate-analysis.md; branch on e554622, uncommitted). Starting point (R, si64 strict np=2/4/8): certificate screen 3.553 / 4.460 / 6.904 s (= certificate outside phases), coarse atomic start 2.263 / 2.420 / 2.700 s, coarse diagonalization 1.48x from np=2 to 8 (8.4 s at np=8). Code (C): the screen's degenerate-cluster loop (full-grid complex FFT per cluster band, s(s+1)/2 full-grid products) and the whole atomic start (span, two-pass Gram-Schmidt, complex H Q because ATOMIC-REAL is threads-only, Gram) are replicated single-thread work on every rank; the screen's 3-4 collectives per block are not the cause. Switches (default off): `BAANDI_M2_SCREEN=dist` (cluster entries from `band_values` on each rank's FFT planes, one tree_sum per cluster; rounding class of the screen estimate only), `BAANDI_M2_ATOMIC=split` (cyclic span columns, distributed two-pass Gram-Schmidt with the serial arithmetic and one broadcast per column, H Q and Gram in column blocks plus one allgatherv; Gram GEMM rounding class) | implemented, opt-in, dev08 pending | Mac O (pure, simulated ranks, single runs; energy and density bits equal to the HEAD build in every arm at si32 np=2/3 and si64 np=8): screen per rank si64 np=2 1.477 -> 0.403 s, np=8 3.399 -> 0.218 s; atomic start 1.201 -> 0.723 s, 1.688 -> 0.512 s; wall si64 np=2 57.93 -> 55.35 s, np=8 40.57 -> 36.34 s; si32 np=8 7.43 -> 6.72 s. Tests: band values vs full-grid FFT <= 3.5e-16 rel, densities bit-equal; split Gram-Schmidt bit-identical P=2-5; SCF |dE| vs threads-only <= 5.3e-15 Ha strict, ranks bit-equal, reproducible. Expected dev08 si64 strict (E): -2.7 / -4.0 / -6.4 s (screen) and -0.9 / -1.4 / -1.9 s (atomic) at np=2/4/8 | Replicated single-thread full-grid work gets slower per rank as more ranks stream it concurrently (shared L3 and memory bandwidth: the MPI-SCALING-BEYOND-8 / MPI-REDUCED-EIGEN mechanism); distributing it removes both the replication and the contention. The coarse diagonalization is a different limit: 80 of 102 Davidson iterations (Mac si64) are coarse, and their replicated reduced eigenproblem and per-column reductions do not shrink with the cutoff; the levers are M2_EIGEN=root/B5 (measure the coarse column) and a distributed dense EVD | Not repeated: CholeskyQR (A52, MPI-M1-GATE1) is not used (the split Gram-Schmidt keeps the serial per-column rule); no capacity change (A13) | this branch |
| MPI-REDUCED-EIGEN-DEV08 | dev08 A/B of BAANDI_M2_EIGEN=root and BAANDI_M2_B5=on (relayed 2026-11-15, master; main 70500a9, base = M2_ORTHO serial + M2_DGEMM + M2_LEAN, f1e-3 seed 601, 3 repeats per arm, si32 five seed-repeats; 0 INVALID). si64 strict wall np=2/4/8 vs 1-base: EIGEN=root +0.15 / -0.51 / -1.18 %, B5 +0.29 / +0.70 / +0.75 %, both +0.09 / +1.20 / -1.74 % (np=8 both: 46.970 -> 46.154 s, spread 0.039); si64 standard all within -1.2..+0.7 %; si32 -0.5..+2.5 % (noise). reduced_eigen stage: np=2/4/8 5.42 / 5.94 / 7.16 s in the baseline, np=8 6.69 (root) / 6.45 (both); it grows with ranks while every other solver stage halves per doubling. EIGEN=root is bit-identical to the baseline (SHA 3/3, dE 0); B5 and both are rounding class (dE <= 1.4e-13 Ha), reproducible per rank count (SHA repeat EQUAL at np=4 in all four cells). Capacity probe (si64 strict np=8, single run): scf.davidson_capacity_per_band=4 reduced_eigen 7.18 -> 4.09 s (-43 %) but Davidson 162 -> 223, H 10104 -> 11490, ortho 3.55 -> 4.57 s, wall 46.875 -> 45.334 s (-3.3 %, needs repeats). Closure of the np=8 run (arm 1): phases scale 2.10x np=2 to np=8 but ladder_coarse 1.31x, certificate screen 0.51x (gets slower), reduced_eigen 15 % of the wall. | M2_EIGEN and M2_B5 NOT adopted (no consistent gain; np=8 both -1.7 % but np=4 +1.2 %); kept as opt-in. The next MPI ceilings are the replicated coarse stage and the certificate screen (MPI-COARSE-CERT) and the capacity lever. A capacity sweep 3/4/5/6 with repeats is queued (new premise against A13: replicated reduced eigen on the MPI path). |
| THREAD-INDEPENDENCE-FALSE-ALARM | Alarm and resolution (2026-11-15): the P0.1 harness reported density SHA T=1 != T=8 on all 24+ T=8 rows including the strict reference. Cause: p01.py reset its `store` dictionary inside the T loop, so the T=8 pass compared with an empty dictionary (`None == sha` is False). Recomputing from the stored density.bin files: 40 of 40 pairs (si32 and si64 x unperturbed and seeds x REF, S8, S10, S8NV) are bit-identical between T=1 and T=8. sha1t8 stage A (si16 strict f1e-3 seed 601, no switches, 9 binaries across 3fea9e8, 38048b9, 2dbac32, cfaa1bc, 4a41feb, mkl and pure) and stage A2 (MKL forced single sequential thread) all give T=1 == T=8 in density SHA and energy bits. The pure density SHA changes once between commits (3fea9e8 d74df9ad, from 2dbac32 7a128a9f: the P48B/P48D rounding-class defaults), not between thread counts. Mac pure build (7e70c0d) also T=1/2/4/8 identical (gamma, si8, k222s, si16 seed 601). Also found: p01 died on si-elongated seed 901 (KeyError total_energy after 100 SCF stall, now NaN-safe); si-elongated (unperturbed and 901) and syn-p3 seed 903 stall at 100 iterations in all arms including REF (cell property, not the single path); the chain guard `grep INVALID RUN` matched the summary line `INVALID RUNS: 0` (fixed to `INVALID RUN |`). | No thread-independence violation. Lessons: a comparison must be validated against a known-equal pair before it is reported as a failure; harness statistics variables must live outside loops that reset them; guards must not match their own summary lines. Mac vs dev08 energy for the "same" si16 f1e-3 seed 601 differ by 1.5e-4 Ha: the two inputs are different geometries (the P48 agent's generator against the harness inp()), not comparable. |
| MPI-COARSE-CERT-AB-PARTIAL | dev08 A/B of BAANDI_M2_SCREEN=dist and BAANDI_M2_ATOMIC=split (relayed 2026-11-15, master; main 8852dc3, cli-mc 7df44c89, BAANDI MPI, f1e-3 seed 601, 47 of 60 table rows, 0 INVALID; si64 strict np=4/8 pending). Both switches are bit-identical to the baseline in every row (SHA 2/2, dE 0). si32 strict screen_s np=2/4/8: replicated 0.431 / 0.504 / 0.738 s (grows with ranks), dist 0.142 / 0.084 / 0.064 (-67 / -83 / -91 %); wall si32 strict np=2/4/8: -1.88 / -5.08 / -9.97 % with both (np=8 spread 0.049); atomic_s (coarse atomic start) replicated 0.381 / 0.418 / 0.483 s vs split 0.257 / 0.209 / 0.206; si64 standard atomic_s 2.267 / 2.415 / 2.709 vs 1.550 / 1.225 / 1.125 (-32 / -49 / -58 %), wall -1.67 / -5.36 / -7.41 %; si64 strict np=2 screen 3.597 -> 1.025 s (-71.5 %), wall -3.18 %. Standard runs have no screen: arm 2 equals arm 1 in wall and SHA, arm 4 equals arm 3. The B5 and EIGEN arm (5) is not consistently better than arm 4. | ADOPTED (2026-11-15): BAANDI_M2_SCREEN=dist and BAANDI_M2_ATOMIC=split default (bit-identical, the rank penalty of both stages removed); si64 strict np=8, the largest expected effect (screen 14.7 % of the wall), pending. |
| MPI-COARSE-CERT-AB-FINAL | Final dev08 A/B of M2_SCREEN=dist and M2_ATOMIC=split (relayed 2026-11-15, master; 8852dc3, 60/60 rows, 0 INVALID; header "BAANDI MPI, coarse/certificate arms"). si64 strict np=8 f1e-3 seed 601: arm1 47.026 s; screen=dist 40.484 (-13.91 %, spread 0.006); atomic=split 45.256 (-3.76 %); both 39.223 s (-16.59 %, -7.80 s); both + EIGEN=root + B5 37.830 (-19.56 %). Measured against estimate: screen -6.55 vs -6.4 s (6.919 -> 0.372 s, -94.6 %), atomic start -1.59 vs -1.9 s, total -7.80 vs -8.3 s; new broadcast cost of the distributed screen 256 bcast calls 0.662 s (benefit/cost about 10:1, may grow with ranks). Coarse potential_s 2.788 -> 1.197 (follows the atomic start). si64 strict np=2/4 both-arm: -4.02 / -9.12 %; si32 strict np=2/4/8 -1.88 / -5.08 / -9.97 %. Arms 2, 3, 4 bit-identical in all 60 rows (SHA, dE 0); arm 4 np=4 SHA repeat EQUAL in four cells; standard runs have no screen (arm 2 = arm 1, arm 4 = arm 3). Remaining ceilings at np=8 (arm 4, 39.05 s): ladder coarse diagonalization 8.38 s (21.4 %), reduced_eigen 7.43 s (19.0 %), apply 6.98 s (17.9 %), orthonormalize 3.74 (9.6 %), rotate 3.10 (7.9 %); screen plus atomic start now 4.1 %. QE MPI np=8 si64 strict f1e-3 is 5.307 s: BAANDI 39.2 s is 7.4x (MPI against MPI), and slower than QE on ONE core (24.9 s). | ADOPTED, defaults M2_SCREEN=dist and M2_ATOMIC=split (commit 226736f). Conclusion for planning: the MPI gap is dominated by the single-core gap (4.4x) times a weaker scaling (BAANDI 8 ranks 2.9x of its serial path against QE 4.7x); single-core algebra (work order phase B) comes first, MPI engineering has diminishing returns until then. |
| B-LEAN-DAVIDSON | Work order beat-QE phase B (B1 at all widths, B3a, B2/B4 documented) on the serial real Gamma path: `BAANDI_DAVIDSON=lean` (default classic; module `crates/dft-solvers/src/lean_davidson.rs`, Ritz-value solves only; strict residual solves, certificate, complex and MPI solvers unchanged). Per iteration residuals of the unconverged roots only; all-band X and H X only at return (the standard residual guard reads the norms) and contraction (QE's update/last split); unit-norm corrections appended without projection; S = V^T V kept beside H_r with the explicit inverse Cholesky factor T (bordered, 64-column panels, GEMM workspace), A = T H_r T^T by dsyevd/faer or B5, coefficients T^T Y. Pivot rule (fixed numerical rule, no material constant): every new column's squared pivot over its overlap diagonal >= eps/1e-12 = 2.2e-4 (S-orthonormality loss eps/ratio held at the classic 1e-12 contraction gate). Failing append: stage 1 re-projects the block alone through T (one GEMM n m k, rank floor DEPENDENCE_TOLERANCE on the unit input, classic in-block Gram-Schmidt) and retests; stage 2 rotates the basis to its full Ritz basis once (no H application) and projects classically. Counters `convergence.work.davidson_lean` (reprojection and fallback rates, pivot statistics, GEMM GFLOP by site for both paths, full eigensolves and sum m^3). Capacity per band stays the input `scf.davidson_capacity_per_band` (default 6). Ledger premises: P37A-1's rejection rule covered all profiles at strict-level Gram eigenvalues; with one standard path and Ritz-value solves only it no longer applies (P37A-3: 45-50 % fewer multiply-adds per standard step). A13/P31A assumed all roots rotated and an orthonormal basis; B1 at all widths and B3a remove both (capacity sweep B2 is the next step, counts first). A52 is not repeated: the Cholesky acts on the small overlap with a pivot bound and an MGS fallback; the classic MGS is unchanged. PRT1 concerned stale residual flags; Ritz-change flags are exact each iteration | implemented, opt-in; dev08 pending | Mac O, pure build, f1e-3 seed 601 (bench_suite variant, mixing beta line removed), counts exact: si16 std T=1 SCF 7 / Dav 38 / H 886 identical, dE 8e-15 Ha, GEMM 5.84 -> 4.93 GFLOP; si32 std SCF 10 / Dav 52 / H 2381 identical, dE 2.8e-14, GEMM 57.08 -> 47.86 (projection 13.50 + in-block GS 8.25 -> overlap 10.78 + GS 2.70 + small 0.97; Ritz 24.54 -> loop 13.09 + return 8.23 + contraction 1.31); si64 std T=8 SCF 10 / Dav 57 / H 4470 identical, dE 2.8e-14, GEMM 386.8 -> 322.9 (projection 85.7 + GS 58.6 -> overlap 72.3 + GS 20.0 + small 5.8; Ritz 170.1 -> 85.5 + return 67.0); reprojections and fallbacks 0, pivot minimum 0.010 / 0.0055 / 0.0049; full eigensolves and sum m^3 identical (B5 unchanged). Strict profiles: no Ritz-value solves, lean never engages, density SHA identical to classic. Classic default bit-identical to HEAD 2acbce5 (energy, density SHA, SCF/Davidson/H) on si16 std/strict and si32 std/strict at T=1 and si16 std/strict, si32 std at T=8. Lean density SHA identical at T=1/2/4/8 (si16, si32 std). Generality (std, f1e-3 or unperturbed): syn-spd, syn-p3, syn-d3, si-skew, syn-fft-primes, syn-cubic-degenerate converge with identical SCF/Davidson/H counts and |dE| <= 4.4e-15 Ha; stage-2 fallbacks 0 on every cell; syn-cubic-degenerate needs stage 1 in 5 of 7 appends (in-block near-dependence of symmetry-related corrections of degenerate triplets, pivots 2e-6 to 1e-4 against basis ratios 0.6 to 1), dropping 2 columns, H count equal to classic; syn-100-bands and si-elongated do not converge in either arm (100 steps, different non-converged trajectories; syn-100-bands reprojection rate 0.6 %). Fixture tests (cubic, skew, near-degenerate, synthetic multi-projector; capacity 6 and 3): eigenvalues within 1e-9 Ha of classic, identical iteration and H counts, returned residuals equal true residuals to 1e-11, orthonormality defect <= 1e-11 (measured 2-4e-15); both fallback stages exercised by an ill-conditioned fixture; bit-identical at 1/2/4/8 workers | An earlier draft reused the pre-contraction eigenvectors in a fallback after a contraction in the same iteration (non-orthogonal Ritz basis, rebuilds); fixed (after a contraction the basis already is the Ritz basis). A pure Cholesky pivot bound on its own sends symmetric degenerate cells to the expensive basis rotation; the block re-projection keeps the rotation for genuine basis ill-conditioning. Expected (E, gap-analysis rates): lean has no effect on si64 strict (residual solves); si64 standard T=1 MKL about -6 to -9 s of 46 s with PROJECT_WIDE off (the removed projection apply ran on the nano kernel), about -2.5 to -4 s if PROJECT_WIDE is adopted first. Kill rule for dev08: stage-2 fallback rate < 5 %, total H columns +5 % at most (si32 and two non-silicon cells) |
| MPI-CAPACITY-SWEEP-DEV08 | dev08 sweep of scf.davidson_capacity_per_band 3/4/5/6 on the BAANDI MPI path (relayed 2026-11-15, master; 8852dc3 + M2_SCREEN=dist + M2_ATOMIC=split, f1e-3 seed 601, baseline measured in the same session, 180 runs, 60 rows, 0 INVALID). cap6 equals the default everywhere (si32 default 6, si64 default 5 or more: capacity default scales with the band count). cap3 vs default: si32 standard np=2/4/8 -15.8 / -19.3 / -21.9 % (dE 2.04e-9 Ha, 41 % of the 5e-9 bound, SCF 10 -> 9); si64 standard -14.6 / -16.4 / -19.8 % (dE 9.8e-10, 20 %, SCF 11 -> 10, H -5 %, all of red_eigen, rotate and ortho fall); si32 strict np=2/4/8 +7.3 / -0.9 / -9.8 % (dE 2e-13); si64 strict np=2/4/8 -0.09 / -6.83 / -13.94 % (dE 1.0e-11, 10 % of 1e-10; np=8 red_eigen 7.19 -> 2.57 s, ortho 3.56 -> 4.63 s, H +14 %). cap4 non-monotonic: standard -5 to -9 %, strict +7.6 % (si64 np=2) / +4.9 % / -3.6 %, Davidson iterations 162 -> 223. cap5 -1.7 to -9 %, dE <= 1.2e-11, safe but smaller. Mechanism: a smaller subspace shrinks the replicated dense reduced eigenproblem (red_eigen -43 to -64 %), costing more Davidson iterations in strict (ortho +31 % at np=2) while ortho scales with ranks and red_eigen does not; so the gain grows with rank count and cell size and is larger for the standard threshold. | Candidate: capacity 3 (generality check on non-Si cells queued, gencap). Not yet adopted. Reopens A13/P31A with the new premise (replicated reduced eigenproblem, rank-dependent). |
| B-LEAN-DAVIDSON-DEV08 | dev08 A/B of BAANDI_DAVIDSON=lean against classic (relayed 2026-11-15/16, master; main 622a75d, cli-lean-mkl 889c5588, cli-lean-pure 1bb731f8, 564 runs, 84 table rows, 0 INVALID; configurations: explicit standard 5e-9, single path S8 5e-9, S10 5e-11; f1e-3 seeds 601-610 (si16, si64 601-603)). Wall B-lean vs A-classic: MKL T=1 si16 std -6.4 %, S8 -7.0 %, S10 -13.0 %; si32 std -8.4 %, S8 -2.7 % (spread 44 %), S10 -9.3 %; si64 std -12.1 %, S8 -20.6 % (spread 19 %), S10 -14.2 %; T=8 similar (-6.7 .. -17.9 %); pure T=8 larger than T=1 (si16 S10 -16.1 %, si32 std -10.4 %, si64 S8 -20.9 %). Work identical in all 84 rows (SCF, Davidson, H equal integers), |dE| equal except one 3 % difference at 5e-12. Kill rules: stage-2 fallback rate 0.0 and stage-1 reprojections 0 in all 36 blocks (pivot_min 2e-3 .. 1.7e-2), total H increase 0 %. Density SHA T=1 == T=8 in 42 of 42 rows. Kill-rule cells (syn-spd, syn-p3, syn-d3, syn-cubic-degenerate) stall in both arms; si-skew identical (6 / 24 / 144); syn-cubic-degenerate A 204 / 1608 vs B 203 / 1607 (one Davidson iteration, both stalled). Arm C (classic + PROJECT_WIDE) -12.3 % .. +4.5 %, sign depends on configuration: not the source of the lean gain. | ADOPTED: BAANDI_DAVIDSON=lean default (Ritz-criterion solves only; residual-criterion strict solves unchanged). Non-Si converging cells (mgo-40, k444s, k222s) in the generality queue. Side result: pure builds beat MKL builds at T=8 (si64 std-explicit pure 17.78 s vs MKL 21.26 s). |
| D1-LADDER-HANDOFF-KNOBS | Work order beat-QE 2026-11-14 phase D1: ladder hand-off scan switches, defaults bit-identical to the adopted ladder. `BAANDI_LADDER_COARSE_STOP=<x in (0,1)>` (default 1e-4), `BAANDI_LADDER_COARSE_RULE=fixed|atomic` (atomic: max(level, distance between full-cutoff and prolonged coarse atomic start densities), needs scf.start_density=atomic), `BAANDI_LADDER_FIRST_TOL=auto|cap|coarse|next` (strict seeds via the new `seed_strict_tolerance` on the residual schedule; standard auto already = next); fraction stays `BAANDI_CUTOFF_LADDER=<f>`. results.json cutoff_ladder block adds coarse_fine_step_equivalents, coarse_wall_fine_step_equivalents, total_step_equivalents, total_h_full_grid_equivalents, fine_scf_steps/h, atomic_truncation_distance, handoff_to_final_density_distance. LADDER-S2-SECANT not repeated. | implemented, defaults unchanged; dev08 scan requested | Mac counts (O), f1e-3 seed 601, uniform start, T=8 engine bench `ladder_d1_scan`; fine SCF / total step eq (fine + coarse H eq per fine step) / total H eq. si16 standard: off 13/13/1290; E1 0.5 stop 1e-4 auto 8/10.55/1134 (default), coarse 8/10.67/1096, stop 1e-3 auto 8/9.71/1141, stop 3e-3 or 1e-2 11/12.2/1227; E1 0.7 stop 1e-3 coarse 8/10.72/1282. si32 standard: off 14/14/2668; 0.5/1e-4 auto 9/11.84/2359 (default), 0.5/1e-3 10/11.69/2376, 0.5/1e-2 12/12.91/2585, 0.7/1e-3 auto 9/11.65/2254 (best). si16 strict: off 15/15/3242; 0.5/1e-4 auto 11/13.32/2671, first=coarse 10/12.36/2512 (best). si32 strict: off 16/16/7135; 0.5/1e-4 auto 11/13.33/5580 (best H), coarse 11/13.23/5820; looser stops cost 3 to 4 fine steps (14-15). All converged; strict |dE| <= 5.3e-12, standard |E - E_strict| <= 3.9e-9 Ha. Atomic rule: d = 1.3e-5 to 1.8e-5 while the measured coarse-to-fine density distance is 3.1e-4 to 4.6e-4 (20 to 30x larger), so the rule never binds at 1e-4: no effect, the atomic start density is a poor proxy for the coarse-to-fine difference. | Mac: no arm reaches fine steps <= QE minus two at si32 (QE 11 at f1e-3: needs 9); stopping the coarse stage earlier than 1e-4 trades about 1 coarse step equivalent for 1 to 3 fine steps (net loss on strict, neutral on standard); fraction 0.7 costs 1.3 to 1.9 more coarse step equivalents than 0.5 without fewer fine steps except si32 standard stop 1e-3. first=coarse helps strict si16 by one step, mixed elsewhere. The measured handoff distance (3e-4 to 5e-4 at E1=0.5) is the quantity a previous-run rule would use; it is above 1e-4, which supports a coarse stop near 3e-4 to 1e-3 only if dev08 confirms (si32 standard 0.5/1e-3 equals the default within 0.15 step eq). Mac numbers are counts; dev08 decides. |
| C2-REUSE-PSI | Work order beat-QE C2 (no double transform of the same orbitals): `BAANDI_REUSE_PSI=on|off` (default off until dev08; rounding class). A lean return keeps X and H X (its return rotation of the stored H images, workspace `KeptStart`); the density build keeps psi(r) (unscaled inverse pair transform, `SplitPairDensity::accumulate_density_pair_split_keep`) of the leading occupied pairs; the next threads-only standard solve of the same k state starts from H_new X = H_old X + (V_new - V_old) X with the potential change applied to the kept psi(r) by one forward transform per pair (`ForwardProduct`, `PairLocalWorkspace::apply_block_kept`), no inverse transform, nonlocal product or kinetic pass; unkept pairs (empty bands, odd tail) take the full pair transform with the potential change. The solver takes the start only when `initial` equals the kept X to 1e-14 and uses X itself as the start basis (S = X^T X factored as any fresh lean basis). Fallback to the plain start: first solve of a state, rebuilt operator (new atoms/projectors/basis/grid, other ladder stage), MPI path, strict (residual) or classic solve, failed return, initial differing from X; memory gate 2 GiB of kept grids per k state (above it H X is still reused, every pair transformed both ways). No allocation per step (grids and X/H X slots allocated once per state). Counters `convergence.work.reuse_psi`, FFT `kept_*` work timers | implemented, opt-in; dev08 pending | Stage 0 (Mac pure, lean, standard, f1e-3 seed 601, counts exact): si32 bands 68 (34 H pairs, 32 occupied pairs), 10 fine + 7 coarse solves, initial applies 680 of 2381 fine H columns (28.6 %) and 476 of 1553 coarse; si64 bands 132 (66 / 64 pairs), 10 fine + 6 coarse, initial applies 1320 of 4470 fine (29.5 %) and 792 of 2646 coarse. Si64 T=1 per-op cost: pair inverse 1.08 ms, forward 0.98 ms, nonlocal 0.82 ms per column. With the switch: reused starts si16 11 / si32 15 / si64 14 of 13 / 17 / 16 solves (first coarse and first fine solve not reusable), 0 declined; inverse transforms saved 176 / 480 / 896, columns without nonlocal and kinetic work 396 / 1020 / 1848; fine H columns 886->670, 2381->1769, 4470->3282; kept grid bytes 24 MB / 96 MB / 382 MB. SCF and Davidson counts identical (7/38, 10/52, 10/57), |dE| 4e-15 / 0 / 8.5e-14 Ha; density SHA identical at T=1/2/4/8 (si16, si32). Mac T=1 wall, 7 alternating pairs, medians: si32 5.382 -> 4.884 s (-9.3 %, 7/7), si64 26.906 -> 24.265 s (-9.8 %, 7/7). Generality (standard): si-skew, syn-spd, syn-p3, syn-d3, syn-fft-primes, syn-cubic-degenerate converge with identical SCF/Davidson counts and |dE| <= 4e-15 Ha; k222s not Gamma (no reuse); si-elongated(-u,-b05), si-slab-b05, syn-100-bands-u fail in 100 steps in both arms; syn-100-bands (= -b05) converges in 89 steps off but not on, and fails with three other rounding perturbations of the off arm (DENSITY_PAIR=off, WIDE_DOT=off, classic Davidson): a trajectory-sensitive stall, not a reuse defect. Synthetic toy (CLI test, 26 steps): reuse moves E by 1.6e-9 Ha, DENSITY_PAIR=off by 2.2e-9 (both in the standard class) | Unit gates: H_new X against a plain apply relative 1.2e-16 worst (gate 1e-13; three cells, s-f projectors, even/odd blocks, kept and unkept pairs, bit-identical at 1/2/4 threads); kept grid + forward_product = apply_local to 1e-14 for the pair and grid plans, density bits unchanged by keeping. Measured Mac gain exceeds the per-op estimate (about 7 %): fewer passes over the 382 MB working set. dev08 request in the C2 report |
| C3-NL-ENERGY-FROM-PROJECTIONS | Analysis only (work order C3, with C2): can the audit's nonlocal energy come from cached projections at no cost? Today (NL_ENERGY_REAL default on) every map computes P^T X once (2 n n_proj bands flop: si64 1.6 GFLOP per solve, about 60 ms at T=1 Mac, about 0.75 s or 2.7 % of the si64 Mac T=1 wall over 16 solves, E). The lean solver keeps no projections (B6 caches them on the classic path only), and C2 removes the initial apply, so no free P^T X exists yet | not implemented; candidate C3a | Premise change by C2: P doesn't depend on V, so P^T X of the kept start is the previous return's P^T X, and P^T X at return = (P^T V) C (n_proj m bands, about 2 % of the n-length product at si64) if the lean solver keeps P^T V per basis column through the existing `copy_last_projections` hook and rotates it with H V at contraction and return. Expected (E) 1.5 to 2.7 % on si64; not the ENERGY-IDENTITY mechanism (that removed the nonlocal term via the band-sum identity and was neutral) |
| C1-PAIR-FFT-PRUNED-DEV08 | dev08 BAANDI_PAIR_FFT default (grid) vs pruned, MKL build (relayed 2026-11-16, master; 4a41feb, cli-p01-mkl; f1e-3; si32 seeds 601-610, si64 601-603; arms measured with the timers off). Wall pruned vs default: si32 standard T=1 -13.94 % (10/10 seeds, p 0.002), T=8 -20.26 %; si32 strict T=1 -10.64 %, T=8 -14.89 %; si64 standard T=1 -16.19 %, T=8 -21.56 %; si64 strict T=1 -11.69 %, T=8 -15.57 % (3/3 seeds). SCF, Davidson and H identical, max |dE| 1.7e-13 Ha, density SHA differs by design (different transform path). Timers (STAGES on, seed 601, thread seconds): wave FFT is 33.5 % (si32) and 33.3 % (si64) of the T=1 wall, 13 % and 18 % at T=8 per thread; pair_inverse 1.6-1.7x pair_forward; timer overhead +0.13 % at si32 T=1, energy and SHA identical. | ADOPTED: BAANDI_PAIR_FFT=pruned default on MKL builds (pure builds keep grid; no pure measurement). Lead for C2: the inverse transform is the expensive half. |
| AGGREGATION-BIAS-RULE | Harness bias found by master (2026-11-16) in the MgO supercell capacity table: the default arm had 3 repeats on the unperturbed geometry and 1 on three displaced ones (6 values), the other arms 1 per geometry (4 values); the raw median weighted the default arm toward the fast unperturbed geometry and printed capacity 3 at +9.7 % while the per-geometry ratios give -2.6 % (T=1 standard), -1.8 % (S8), -5.2 % (np=8, all four geometries negative). Demo with the same numbers: raw median +15.8 %, per-seed aggregation -0.5 %. Per-geometry T=1 signs are mixed (+3.2, -4.0, -7.0, -1.2 %); displaced geometries are single runs. closure.total_s (process start excluded) differs from the table wall. | RULE: compare arms by the median over seeds of per-seed medians (med_by_seed), never by a raw median over runs, unless every arm has the same seeds and repeat counts. Harnesses mixer.py, c2.py, d1.py and the rerun of mgsc use it; earlier tables with one seed per block (cap, M2) are unaffected; lean and C1 tables are being recomputed from stored files. |
| B-LEAN-DAVIDSON-RECOUNT | Correction of B-LEAN-DAVIDSON-DEV08 by per-seed aggregation (relayed 2026-11-16, master; stored closure.total_s, no new runs). The raw-median aggregation weighted the classic arm (first seed repeated 3 times) toward its first seed. Corrected B-lean vs A-classic (median over seeds of per-seed medians), std-explicit / S8 / S10, T=1 and T=8: mkl si16 -6.5/-6.7, -6.3/-6.9, -6.2/-6.3; si32 -8.6/-8.0, -9.0/-8.5, -9.4/-9.2; si64 -12.2/-9.3, -12.9/-10.3, -12.2/-9.9; pure si16 -6.3/-9.5, -4.9/-9.1, -5.5/-9.4; si32 -7.5/-10.8, -7.1/-11.3, -7.3/-12.0; si64 -10.4/-11.4, -11.2/-12.5, -11.1/-12.1 (range -4.9 to -12.9 %). Withdrawn: "maximum -20.6 %" (mkl si64 S8 T=1 is -12.91 %), "si32 S8 not significant" (now -9.01 %), "arm C loses in S8" (arm C PROJECT_WIDE is -1.0 to -4.3 % everywhere; two sign flips mkl si32 S8 +3.59 -> -3.09 and +4.63 -> -1.29). Unchanged: work identical, fallback 0, density SHA T=1 == T=8 in 42 of 42, gain grows with cell size, pure gains more at T=8 and mkl at T=1. Other tables audited without bias: C1 pair-fft (raw equals by-seed), p01, m2ab, m2rep, cap, QE sweep. | Lean stays adopted; the quoted gain is -5 to -13 %, not -3 to -21 %. All later tables use med_by_seed. |
| MPI-CAPACITY-RULE | Adopted rule (2026-11-16): with four or more MPI ranks and no explicit scf.davidson_capacity_per_band the CLI uses 3 Davidson basis columns per band (switch BAANDI_M2_CAPACITY=auto|off, default auto). Evidence: capacity sweep (si64 strict np=4/8 -6.8 / -13.9 %, si32 strict np=8 -9.8 %, np=4 noise, si64 standard np=4/8 -16.4 / -19.8 %, si32 standard -19 / -22 %) and the 132-band MgO supercell (np=8 standard -5.2 %, S8 -6.2 %, per-geometry all negative but one; serial T=1 inconclusive, mixed signs). np=2 and serial runs keep the default 6 (si32 strict np=2 +7.3 %). Mechanism: the replicated dense reduced eigenproblem is a ceiling that grows with ranks; a smaller basis shrinks it at the price of more Davidson iterations (strict) or none (standard, si64). Accuracy: standard dE up to 2-3e-9 of the 5e-9 bound at si32 (61 % worst at cap4), MgO cap3 22 % of the bound, strict dE <= 1e-11. | ADOPTED as a rank-count rule (not a material constant); explicit input always wins; off by BAANDI_M2_CAPACITY=off. Serial-path capacity stays an open question (perturbed MgO geometries as 3 repeats would settle it). |
| MIXER-HARTREE-GENERALITY | dev08 mixer arm B (pulay, beta 0.7, history 8, mixing_metric hartree; approximation of QE's mixer) against the default mixer, seven cells, standard and S8, f1e-3 seeds 601-603, T=1, med_by_seed, strict references 3/3 (relayed 2026-11-16, master; cli-lean-mkl 622a75d; 0 INVALID). Wall: si16 std -9.7 %, S8 -7.2 %; si32 std -8.2 %, S8 +14.50 % (SCF 10 -> 11, H 2598 -> 3142; the only clear loss); si64 std -3.3 %, S8 -0.9 %; mgo-40 std -21.3 %, S8 -19.4 % (SCF 13 -> 9; default there is kerker-pulay so arm B also removes the Kerker preconditioner: not separable); si-k444s -4.6 / -4.9 %; si-k222s -2.1 / -2.7 %; si-skew +0.9 / -6.7 %. Accuracy margin shrinks (fewer SCF steps): si64 std |dE| 37 -> 50 % of the 5e-9 bound, si16 std 14 -> 38 %, mgo-40 3 -> 24 %, si-k444s 14 -> 0.8 %. | NOT adopted as a default yet: 12 of 14 blocks gain, but si32 S8 loses 14.5 % and the accuracy margin approaches half the bound at si64. Open: add a kerker-pulay + hartree arm on mgo-40 to separate the metric from the Kerker removal; a rule for mixer choice (class and kind of system) stays an open design item. |
| C3a-LEAN-PROJ-CACHE | Review F10 / backlog C3a (2026-11-16, on a4c7a8f): `BAANDI_LEAN_PROJ_CACHE=on|off` (default off until dev08). The serial lean Davidson keeps P^T V of its basis columns (`lean_davidson::ProjCache`, count x capacity): copied from the operator (`copy_last_projections`) after the start apply, every appended block and a rebuild; rotated by the contraction coefficients through the fixed-panel GEMM workspace; the stage-2 fallback rotation (m columns) invalidates it until a full apply. At the return it forms P^T X = (P^T V) C (count x m x bands GEMM) and the workspace reports it (`returned_projections`, cleared at the start of every solve, only for lean returns); the REUSE-PSI kept start carries the kept P^T X (same operator: the engine rebuilds the workspace with the operator). The engine's NL-ENERGY-REAL audit reads it through `RealHamiltonian::nonlocal_expectation_values_from_projections` (same coupling-sum code) instead of a fresh n x count x bands product; any miss falls back to the direct product. Counters: results.json `davidson_lean.proj_cache` (returns, misses, kept_loads, columns_copied, audits_from_cache, audits_direct). | implemented, default off; dev08 A/B requested | Mac O: unit tests (4 atoms with 3/2/4/3 projectors, off-diagonal D, bands 7 and 9, capacity per band 2 (contractions 2-3) and 6, T=1 and 3): max abs P^T X difference against a direct product 1.8e-15 to 2.9e-15, band nonlocal energies 2.7e-15 to 7.1e-15 Ha (gate 1e-12), kept start 2.2e-15 / 6.2e-15 Ha; projections NOT bit-identical to the direct product (different summation order; rounding class), bit-identical across thread counts; the solve (eigenvalues, vectors, iterations, H count) bit-identical with the cache on and off; physics test: energies from the latest apply's projections equal `nonlocal_expectation_values` bit for bit. SCF (si8, one atom displaced by 1e-4 fractional, standard profile, Mac release, T=1 and T=4): total energy on-off 3.6e-15 Ha (REUSE_PSI=on: -7.1e-15 Ha, 9 kept loads), SCF 5, Davidson 28, H 351 (271 with REUSE_PSI) identical, density SHA identical (the audit energy does not feed the density), 11/11 audits from the cache, 0 misses. Expected gain (E, from C3-NL-ENERGY-FROM-PROJECTIONS): one n x n_proj x bands GEMM per map, si64 1.5-2.7 % T=1; Mac wall on si8 is noise. si32 f1e-3 standard T=1 Mac pure, 5 alternating pairs, load 5-7 (busy host, O): median 4.988 -> 4.968 s (-0.4 %, within noise), dE exactly 0 Ha, 13/13 audits from the cache. Lesson: reusing B6's scalar `rotate` loop would cost n_proj x m x bands serial multiply-adds per return (si64 about 5e7), comparable to the GEMM it saves at MKL speed, so the cache rotates through the GEMM workspace. |
| F11-GEMM-DISPATCH-COUNTERS | Review F11 (2026-11-16): count-only route counters for every real GEMM wrapper of `RealGemmWorkspace` (`dft_numerics::gemm_dispatch`: per site x route calls, 2mnk flops, max m/n/k; routes mkl, mkl_wide, pure_build, pure_switch, pure_output_bound, pure_panel_bound, pure_long, pure_short, pure_layout, pure_direct, narrow; pure builds also count the route an MKL build would take, `rows_as_if_mkl`), results.json `convergence.work.gemm_dispatch`, accounting build records `real.<site>.<route>` and results.json `gemm_accounting`. The dgemm decision is a clause-by-clause decomposition of the old condition (MKL unit test over 7,200 combinations; type-checked on the Mac). Analysis docs/performance/gemm-dispatch-analysis-2026-11-16.md | implemented (no behaviour change); fixes proposed, defaults unchanged | Mac O (counts exact) + D (MKL route of the same shapes, MKL orthogonalization defaults MGS=blocked, ORTHO_BLOCK=off), si16/si32/si64 f1e-3: standard (lean) pure fallback 0.7 / 1.3 / 0 % of the real GEMM flops (MGS-B cross-panel apply n x 16 only, si16/si32); strict 10.2 / 11.9 / 10.4 % (projection apply X -= Q S: 79/156, 73/116, 105/249 calls carrying 1.90 / 15.26 / 178.79 GFLOP, plus MGS-B 0.31 / 3.97 GFLOP at si16/si32). Rule: an A B output n_pw x c falls back when c > 36 bands^2 / n_pw, so full-width projections fall back whenever n_pw > 36 bands (always for plane waves; Si 30 Ry has 165-181 PW per band). Certificate solves (classic residual Davidson) carry the same site. si8 energy and density SHA bit-identical with the counters. Proposal: dev08 A/B of `BAANDI_PROJECT_WIDE=on` on strict and S8 runs, then a shape rule (whole-output bound only for A^T B) instead of a per-site switch; E si64 strict T=1 about 13 s of single-core work at assumed 10 vs 40 GFLOP/s. |
| F21-STAGE-TIMERS-SWITCH | Review F21 (2026-11-16): `BAANDI_STAGE_TIMERS=on|off` (default on = unchanged behaviour; `dft_numerics::stage_clock`). Off = count-only: the MPI collective wrapper (`SharedComm::timed`), the distributed Davidson stages, the M2 operator stages and the engine M2d stages count every entry but read no clock (`stage_clock::start` returns None; no `Instant::now`/`elapsed`), stage and collective seconds stay 0. New in both modes: per-stage entry counts (results.json `parallel.stages[].calls`, `parallel.m2.counters["calls.<stage>"]`) and collective payload bytes (`parallel.collectives[].payload_bytes`: f64 values delivered on this rank x 8: allgatherv output, tree_sum vector, bcast buffer, alltoall send buffer); `parallel.stage_timers` records the mode. The atomic-start clock shared with the phase clock is still read (phase timings unchanged) and only kept for the M2 stage in time mode. | implemented; default unchanged; dev08 on/off A/B requested | Mac O (simulated 2 ranks, one-atom Gaussian NC fixture, debug): off and on arms count identical stage entries, collective calls (allgatherv 206, tree_sum 2532, alltoall 764, bcast 0, barrier 0) and payload bytes (588,832 / 167,888 / 813,184 B), M2 counters equal; every second counter stays exactly 0 in the off arm and is positive in the on arm; energy, history length and density bits identical between modes. Overhead not measured on the Mac (no MPI runtime here); the off mode exists so dev08 can time it (timers-on minus timers-off on the same binary) and run acceptance with the instrumentation removed. Lesson: tree_sum dominates the call count on small cells (2532 calls for 168 kB), i.e. latency, not volume; the bytes column makes that visible next to the seconds. |
| SINGLE-PATH-DEFAULT-FLIP | User decision 2026-11-16 (review 2026-11-16 F01-F03): `scf.convergence_profile` removed (input error naming `standard_threshold_ha`, `verify`, `BAANDI_REFERENCE_STRICT`); every classic run uses the standard path (5e-9 Ha) with the fixed-point verification; `BAANDI_SINGLE_PATH` replaced by `BAANDI_REFERENCE_STRICT=on|off` (default off; switch count stays 66); strict fallbacks (relax, stress, spin, Fermi-Dirac, locking auto, runtime modes, non-run commands) record `single_path.status = not_applicable` with the reason; results carry `execution` (requested/effective profile, verification, solver lean/classic/real/complex, backend, parallel model, non-default switches); MPI multi-k/spin/FD runs report `parallel.model = mpi_replicated` (F02); the residual guard is documented and named as an empirical estimator (F03, 2x2 counter-example in ADR 0043). Mac only: `scripts/check.sh` gates; no timing measured. | Adopted by user decision, not by a measured speedup. Default Si results change only because the default path changed (standard class instead of strict: SCF count, energy within 5e-9 Ha, Pulay beta 0.5, density hash, verification cost inside every timed run); `BAANDI_REFERENCE_STRICT=on` is the former strict default bit for bit and `verify = "off"` the former explicit standard profile. Strict baselines (`docs/performance/baseline.json`, QE validations) need the reference switch; `bench_scf.py` refuses a profile mismatch. Owed on dev08: the p01-single-path list as measurements of the new default (delivered error, verification pass rate and cost, forces). |
| ONEPATH-1 | User decision 2026-11-16 (strict leaves production in stages; stage 1 forces and relaxation, ticket `tickets/ONEPATH-1-FORCES-RELAX.md`, ADR 0043 amendment ONEPATH-1): fixed-cell relaxation runs every ionic step on the standard path (`scf.standard_threshold_ha`, default 5e-9 Ha) with the fixed-point verification (failed verification continues that step's SCF; warm starts and extrapolation unchanged); the `relaxation` strict fallback is removed; variable-cell modes and `output.stress` stay strict (reasons `stress (variable-cell relaxation)`, `stress`); `verify = off` with `[relax]` is an input error; `forces.preconditions` records threshold, attempts and the last Coulomb residual; relaxation frames carry a deterministic `verification` record (no wall time, resume stays bit-identical) | adopted (default; rollback `BAANDI_REFERENCE_STRICT=on`) | Mac, synthetic fixtures, forces standard vs strict reference: cubic two-H (local ions) max abs dF 1.54e-7 Ha/Bohr at 5e-9 (5.8e-6 of max abs F 2.64e-2; dE 2.1e-11 Ha; 1 attempt) and 2.69e-8 at 5e-11 (1.0e-6); skewed four-atom two-species nonlocal fixture 4.81e-6 (debug) / 5.28e-6 (release) Ha/Bohr at 5e-9 (3.3e-5 / 3.7e-5 of 0.144; dE 4.9e-10 / 5.9e-10 Ha; 3 attempts, V5 continuation exercised) and 6.5e-7 / 5.9e-7 at 5e-11 (4.5e-6 / 4.1e-6; 2 attempts); bounds fixed before measuring 1e-3 (5e-9) and 1e-4 (5e-11) relative. Relaxation, compact skewed two-H fixture, force_tol 1e-4: final positions 1.6e-6 Bohr apart (bound 1e-3), 17 ionic steps both, SCF iterations 41 standard (18 verification attempts) vs 122 strict, H vectors 539 vs 846, wall 1.07 vs 0.92 s (noise level). First relaxation fixture (cubic two-H, flat minimum at the cell-symmetric point, curvature about 2e-4 Ha/Bohr^2): both converged, 5.4e-3 Bohr apart, final max force 9.1e-7 vs 2.1e-9; there the stopping force sets the position, not the SCF path; fixture replaced by a stiff one rather than widening the bound. Existing oracles at the default 5e-9: C3 extrapolation identity max abs dF 1.07e-6 (bound 1e-6, strict class) and the 1e-12 Ha energy-scan relaxation exhausts 200 steps; they now ask for 5e-11 (dF 9.8e-8 / 6.5e-8) and 5e-13 (distance error 3.3e-5 Bohr, bound 1e-4) | the default threshold resolves forces to about 1e-5 relative on these fixtures; ionic energy tests far below the threshold need a tighter `standard_threshold_ha` (documented in docs/cli.md, docs/user/relaxation.md). dev08: `docs/performance/requests/onepath1-forces-relax.md` (Si8/Si16 forces seeds 601-610 at 5e-9 and 5e-11; Si8, slab, skew relaxations). Next stages: stress, spin/metals, engine default |
| M2-WORKSPACE-F06 | Review 2026-11-16 F06: persistent `DistributedDavidsonWorkspace` (`BAANDI_M2_WORKSPACE=on|off`, crates/dft-solvers/src/distributed_davidson/workspace.rs). Basis and H products in capacity-column storage (in-place append instead of copying the whole basis per iteration), reduced matrix with fixed leading dimension and border update (no per-iteration zero matrix and no `reduced_prev` clone; packed m x m copy for the dense solver), fresh products / lowest coefficients / Gram operands as slices, persistent Ritz, residual, correction and serial-order orthonormalization scratch, GEMM/EVD/B5 workspaces rebuilt only when (max rows, capacity, pool) changes; one workspace per thread (= per rank), reused across solves. dft-mpi: `tree_sum_f64` receive scratch and `allgatherv_f64` count/displacement arrays kept per thread; simulator slots reused. Same arithmetic order. (C) bit-identical to off: tests/distributed_workspace.rs at ranks 1, 2, 3, 4 on the entry-cost stack, root, B5, root+B5, classic and whole-block serial order, strict and standard, capacities 6 and 3, cold and warm starts, repeated solves of different shapes on one workspace (all eigenvector, eigenvalue, residual, history bits and work counts equal); dft-engine `m2_workspace_is_bit_identical_in_scf` (cubic and skewed, strict and standard, P=2/3: energy and density bits equal; 195 of 232 solves warm). (O, counting allocator, toy n=240, 8 bands, steady-state iterations): heap bytes per iteration on ranks that do not run the dense solve 0 (off 177-285 KB, largest single allocation 51-77 KB); the solving rank keeps 33 KB per iteration, all inside faer's EVD packing (pure build only; MKL dsyevd uses preallocated storage). (O, Mac si32 strict f1e-3 seed 601, pure build, simulated ranks, same binary sha256 f382796b..., load 9.5-12, alternating off/on, single runs): energy bits and density FNV identical in all 12 runs (np=2 0xc05fad02a30818d2, np=4 0xc05fad02a30818d0); wall np=2 off 21.99/15.95/23.96 on 16.60/19.74/21.99 s (no signal under this load); np=4 off 17.03/11.48/14.11 on 13.76/9.83/13.07 s (3/3 pairs faster, -7 to -19 %, load-confounded); stage residual 0.926/0.328/0.493 -> 0.503/0.236/0.414 s and reduced_gram 1.208/0.713/1.009 -> 0.821/0.589/0.812 s per rank at np=4. Max RSS inconclusive (1.2-2.5 GB, varies with load). | implemented, default ON (bit-identical, storage only; rollback `BAANDI_M2_WORKSPACE=off`); dev08 confirmation requested | Storage changes in the distributed solver can be proved exact by construction (same GEMM calls on views, same reduction calls); the remaining per-iteration heap traffic on pure builds is the faer EVD's internal packing. Classic (non-serial) orthonormalization and the M1 apply still allocate inside (not default). | this branch |
| M2-SELECTED-ROOTS-F05 | Review 2026-11-16 F05 step (1): selected-root rotation in the MPI distributed Davidson, standard (Ritz-change) solves only (`BAANDI_M2_SELECTED_ROOTS=on|off`, implies the workspace path). Converged roots decided from the Ritz-value change before any rotation; V C, (H V) C and residuals only for the unconverged roots (coefficient columns gathered); all bands rotated at return (every band's residual reported and returned) and at restarts (when a restart follows, one full rotation serves the selected residuals and the contraction basis; without that, the first version did 11 % MORE rotation flops than the full path on capacity-3 cases). Strict residual path unchanged (bit-identical, test). Step (2) of F05 (distributed overlap S / bordered Cholesky / generalized reduced solve) not started. (C) tests/distributed_selected_roots.rs, ranks 1-4, scattered / clustered / exact-degenerate (groups of 3 and 4) spectra, capacities 6 and 3: iterations, H columns, restarts, contractions, max basis, convergence flags and root ordering identical; eigenvalues and every returned residual bit-identical on the pure build (rounding class in principle: column subsets of one GEMM; MKL dgemm may block differently by width); returned residuals equal an independent residual from the full matrix to 1e-9 relative; tree-sum count unchanged, payload -0.1 to -0.5 %; rotation flops -8 to -27 % on these toys. (O, Mac si32 standard f1e-3 seed 601, pure build, simulated ranks, binary sha256 6478e21e..., load 4.7-5.2, alternating, 2 pairs per np): SCF/Davidson/H 15/69/2723 both arms, energy bits and density FNV identical (np=2 0xc05fad02a30717d8, np=4 0xc05fad02a30717e8); rotation columns per rank 4692 -> 2723 (-42 %), rotation flops 1.728e10 -> 8.457e9 (np=2) and 8.641e9 -> 4.229e9 (np=4) (-51 %); stage rotate 0.607 -> 0.344 s (np=2), 0.316 -> 0.178 s (np=4); residual 0.05 -> 0.03 s; wall np=2 4.613/4.546 -> 4.340/4.265 s (-6 %), np=4 3.365/3.357 -> 3.276/3.208 s (-3 to -4 %); tree_sum calls 8406 per rank both arms. | opt-in (default OFF until the dev08 A/B; adopt if np=2/4/8 median wall improves with standard-class accuracy) | The F05 rotation saving is real but bounded by the rotate stage share (Mac si32 standard: about 13 % of the wall at np=2); H columns and reductions are unchanged by construction, so step (2) (orthogonalization removal) is the larger lever for reductions. When a restart follows, rotate all bands once instead of selected-then-all. | this branch |
| M2-REDUCTION-COUNTS-F08 | Review 2026-11-16 F08 (analysis only, no algorithm change): per-site counters of every tree sum in the distributed Davidson (communicator wrapper `Counted`, count only, sites reduced_gram, residual_norms, ortho_normalize, ortho_projection, panel_cross, panel_dot, panel_norm, restart_defect, apply, setup; payload buckets 1 / 2-8 / 9-64 / 65-512 / 513-4096 / >4096 doubles); `reduction_stats()`, results.json `parallel.reductions`. (C) test: per-site sums equal the communicator's own count and payload exactly (P=2, 3). Results unchanged: Mac si32 energy bits and density FNV identical to the F05/F06 runs. (O/C, Mac si32 f1e-3 seed 601, pure build, simulated ranks, counts per rank; counts are trajectory properties, not machine properties): strict np=2 23451 solver sums (17 solves): panel_dot 12758 (54 %, payload 1 or 2-8 doubles in 100 %), panel_norm 7555 (32 %, 1 double), panel_cross 1582, ortho_projection 368, reduced_gram 293, residual_norms 297, ortho_normalize 293, apply 297, restart_defect 8; strict np=4 same total, panel_dot 13726 / panel_cross 614 (the pure-build panel shape rule depends on the local row count, so the split differs with P); standard np=2/4 8391 sums (15 solves): panel_dot 4692/5084, panel_norm 2723. In-panel dot+norm sums are 87 % (strict) and 88 % (standard) of all solver reductions and carry at most 64 doubles each; their payload is under 2 % of the reduced doubles (the Gram/projection/apply sums carry the bytes). Mac tree_sum time 0.40 s (np=2) and 0.73 s (np=4) per rank at si32 strict (4 and 7 % of the wall; latency-bound). | analysis counters in place (always on, count only) | Confirms the analysis table of mpi-reduced-eigen-analysis.md section 6 with exact counts: the call count is the per-column modified Gram-Schmidt (two dots and one norm per column), so latency per call (recursive doubling, item 3) or fewer in-panel reductions (Pythagorean second-pass norm, item 4; F05 step 2 removing the projection) are the levers; payload aggregation does not help. dev08 should read `parallel.reductions` beside `parallel.collectives` to get seconds per call by P. | this branch |

| C2-REUSE-PSI-DEV08-ADOPT | 2026-11-16 | dev08 (master), f42d7e3, 90 rows / 687 runs, INVALID 0, MKL and pure, T=1 | `BAANDI_REUSE_PSI` off vs on, paired same binary: 36 of 36 blocks negative (mkl -1.4..-4.8 %, pure -0.9..-6.3 %; si64 standard mkl -4.8 %, pure -6.2 %); H columns -25..-29 % but wall only -3..-6 % (the reuse removes the initial apply only); SCF and Davidson counts unchanged, dE unchanged (rounding), density SHA identical across T=1,2,4,8 in 54 of 54 comparisons; 21 of 22 armed starts kept. bench (33 cases): syn-100-bands unchanged (both arms unconverged at 100 SCF steps, identical counters); single-run rows (si-gamma -24.7 %, syn-odd-shift -15.2 %, si-elongated -11.7 %, syn-fft-primes +63.7 %) are single runs, repeat requested (benchxchain). | **Adopted as default (on) 2026-11-16** by the adoption rule; `BAANDI_REUSE_PSI=off` stays as rollback. The syn-fft-primes repeat decides whether a shape rule is needed. |
| D1-LADDER-SCAN-DEV08-VOID | 2026-11-16 | dev08, f42d7e3 | Baseline arm spread 83 % of its median (standard) so every percentage of the scan is void (arms with identical SCF/dav/H showed -50 %). | Not used. Re-measure at low priority with >= 3 baseline runs and a 2 % spread gate. |
| MIXER-BETA-DECOMPOSITION-DEV08 | 2026-11-16 | dev08, mgo-40 standard S8 seeds 601-603 | A default kerker 0.3: 7.913 s; B pulay+hartree beta 0.7: -21.87 %; C kerker+hartree 0.3 (metric only): +6.78 %; D kerker+hartree 0.7: -16.76 %. Metric alone is a loss; the gain is beta 0.3 -> 0.7 (-22 points), Kerker removal -5 points. si32 S8 +14.5 % was seed variation in the number of fixed-point verification attempts (seeds 601/602: 1 -> 3 attempts; seed 603: 4 -> 2). | Hartree metric not adopted. Next single-variable test: kerker-pulay beta 0.7 without metric over the general set (beta.py queued); a fixed beta is adopted only if the hold-out set passes. |

| MPI-CAPACITY-RULE-DEV08-CORRECTED | 2026-11-16 | dev08 (master), a4c7a8f, 24 rows, INVALID 0 (an earlier relayed table was withdrawn: master had not read the output file) | `BAANDI_M2_CAPACITY` auto vs off, same binary. np=2: no difference (-0.60 / -0.05 / -0.37 / +0.80 %). np>=4 off vs auto: standard +15.24 / +18.01 (si32 np4/np8), +16.13 / +17.48 % (si64 np4/np8); strict si64 np8 +12.46 %, si32 np8 +2.91 %, si64 np4 +3.39 % (baseline spread >2 %, held back). Exception: si32 strict np4 off is -4.52 % faster (auto raises Davidson iterations 132 -> 219 and the reduced-eigen saving does not cover it). auto also changes the work (si32 standard np4 SCF 9/dav 44/H 2225 vs 10/48/2329; dE up to 2.04e-9 Ha, inside the standard class; strict dE <= 2.3e-13). | Rule stays adopted (net gain in 7 of 8 np>=4 cells). The si32 strict np4 exception is recorded as a regression to explain; repeat x5 requested. |
| MIXER-BETA07-DEV08-CANDIDATE | 2026-11-16 | dev08, f42d7e3, 21 rows, INVALID 0, standard unless noted | `kerker-pulay beta 0.7` without metric (arm E) vs the case default (A): mgo-40 -16.11 %, si32 -24.11 %, si64 -25.51 %, si-skew -5.35 %, si-k444s +4.80 % (dav 839, H 5526), si64 S8 -19.22 %, si32 S8 (10 seeds) -16.76 %; B (pulay + Hartree metric) +10.41 % on si32 S8 (10 seeds): rejected. Verification attempts min-max A 1-5 -> E 1-2 (part of the gain). |dE|/bound grows (mgo-40 3.0 -> 20.9 %, si32 10.7 -> 14.2 %) but stays inside the class. si-slab, si-elongated, syn-p3 fail in both arms (no evidence). | Candidate default, not adopted yet: still missing the A-arm beta per cell, hold-out cases that converge (slab, elongated, small gap), and the si-k444s +4.8 % explanation. |
| ONEPATH-2 | User decision 2026-11-16 (strict leaves production in stages; stage 2 stress and variable-cell relaxation, ticket `tickets/ONEPATH-2-STRESS.md`, ADR 0043 amendment ONEPATH-2): `output.stress` and relax modes `cell`/`ions_cell` run the standard path with the fixed-point verification (strict fallbacks `stress` and `stress (variable-cell relaxation)` removed; `verify = off` with stress is an input error); the engine accepts a standard state for stress exactly when forces do; the three Pulay cutoff-scan SCFs of a standard caller run standard WITH verification at min(1e-9 Ha, caller threshold), retried once at the caller threshold when looser, no strict fallback (strict callers unchanged bit for bit); `stress.preconditions` standard branch and `pulay_assessment.scan_*` keys | adopted (default; rollback `BAANDI_REFERENCE_STRICT=on`); pre-registered stress criterion FAILED at the default threshold on one fixture (below) | Mac release, synthetic fixtures, stress standard vs strict reference (bound fixed before measuring: 1e-7 Ha/Bohr^3, the synthetic stress oracle cap): cubic two-H local, max abs dsigma 2.06e-9 Ha/Bohr^3 at 5e-9 (9.4e-6 of max abs sigma 2.19e-4; dP 5.0e-10, dPulay 2.3e-13, dE 2.1e-11 Ha, 1 attempt) and 6.47e-10 at 5e-11 (3.0e-6); skewed four-atom two-species s,s,p fixture **1.305e-7 at 5e-9 (FAIL, 2.7e-5 of 4.90e-3; dP 2.4e-8, dE 5.9e-10 Ha, 3 attempts)** and 6.13e-8 at 5e-11 (pass, 1.25e-5; dE 4.5e-11 Ha). Decomposition against a tightened strict reference (energy 1e-12, density 1e-10, eigen 1e-11): default strict 2.1e-11 / 8.6e-11 (cubic / skew), so the exceedance is the standard path's own error; skew 5e-10 -> 5.9e-8, 5e-11 -> 6.1e-8 (same stopping step), 5e-13 -> 5.0e-9 Ha/Bohr^3; SCF iterations skew 21 (5e-9), 44 (5e-10, 5e-11), 75 (5e-13) vs strict 17. Variable-cell relaxations (bound 1e-3 Bohr): BCC two-H cell isotropic stress_tol 1e-7, 7 frames both, cell 1.6e-8 Bohr apart, SCF iterations 12 standard vs 22 strict, H vectors 171 vs 299, wall 0.78 vs 0.74 s; compact skewed two-H ions_cell full force_tol 1e-4 stress_tol 1e-7, 53 frames both, cell 5.4e-6 and positions 3.2e-6 Bohr apart, SCF 98 vs 220, H vectors 1080 vs 3271, wall 5.10 vs 5.01 s (noise level) | stress is first order in the density error (energy second order), so its error falls like the square root of the delivered energy error; the default 5e-9 Ha does not guarantee the 1e-7 Ha/Bohr^3 stress class (one of two fixtures misses by 1.3x); recommendation 5e-11 Ha for stress-sensitive runs, default unchanged pending the user decision and dev08 `docs/performance/requests/onepath2-stress.md`. Next stages: spin/metals, locking auto, runtime modes, engine default |
| ONEPATH-2-STRESS-DEFAULT | Orchestrator decision 2026-11-16 (adoption rule: the default must meet the accuracy class), after ONEPATH-2 criterion a failed at 5e-9 Ha: an input with `output.stress` or `relax.mode = cell/ions_cell` and no `scf.standard_threshold_ha` runs at `STRESS_DEFAULT_THRESHOLD_HA` = 5e-11 Ha (rule by observable, not by material); forces, ions relaxation and plain runs keep 5e-9; an explicit value always wins; `single_path.standard_threshold_ha` / `single_path.threshold_source` (`user`, `stress_default`, `default`) and the relaxation `scf` block record it | adopted, PROVISIONAL until dev08 `docs/performance/requests/onepath2-stress.md` (si8, si16, si-skew, si-slab, strained si8 at 5e-9/5e-10/5e-11 vs a tightened strict reference) | Mac release, synthetic, at the new default (stress bound 1e-7 Ha/Bohr^3 unchanged): cubic two-H max abs dsigma 6.47e-10 (dE 7.8e-13 Ha), skewed nonlocal 6.13e-8 (margin 1.6x; dE 4.5e-11 Ha), scans at [5e-11]x3, 1 verification attempt each; the former KNOWN FAILURE test now passes because the default changed, not the bound. Variable-cell relaxations at the stress default: BCC cell isotropic, cell 9.1e-9 Bohr from strict, SCF iterations 17 vs 22 strict (were 12 at 5e-9), H vectors 189 vs 299; compact skewed ions_cell full, cell 8.9e-7 and positions 3.3e-7 Bohr from strict (were 5.4e-6 / 3.2e-6 at 5e-9), SCF 148 vs 220 (were 98 at 5e-9), H vectors 1280 vs 3271; tetragonal shape oracle relative 6.4e-8 (bound 1e-5). Cost on the skewed stress fixture: 44 SCF iterations at 5e-11 vs 21 at 5e-9 (strict 17) | the tighter default costs SCF iterations where stress is asked for and nowhere else; the margin of 1.6x is thin, so dev08 decides whether 5e-11 holds on the generality set (revisit the default, never the bound) |
| ONEPATH-3 | User decision 2026-11-16 (strict leaves production in stages; stage 3 spin and metals, ticket `tickets/ONEPATH-3-SPIN-METALS.md`, design `docs/design/onepath3-spin-metals.md`, ADR 0043 amendment ONEPATH-3): collinear spin (fixed sector, FD) and scalar Fermi-Dirac run the standard path with the fixed-point verification (strict fallback reasons `collinear spin` / `Fermi-Dirac occupations` removed). Collinear stop estimate E_H[d_n] + E_H[d_m] + E_0[d_m] (magnetization in the same Coulomb metric, uniform moment change weighted like the longest wavelength); residual guard sum_s,k w f |r|^2 with actual (fractional) occupations; ratchet, first-iteration redo and Ritz-change Davidson as closed shell; FD every retained band mandatory, ratchet instead of 0.1 x eigen tolerance per step; verification repeats the whole map incl. Fermi-level solve, passes on the same estimate + number/moment/edge (+ fixed-sector gaps); forces/stress accept a verified spin/FD state; strict-caller stress scans of spin/FD go straight to strict (reference bit-identical) | adopted (default; rollback `BAANDI_REFERENCE_STRICT=on`); all pre-registered bounds met | Mac release, synthetic, bounds fixed before measuring. (a) zero-moment collinear vs closed shell, both standard+verification, skewed cell, LDA-PW and PBE: dE 1.0e-15 / 1.1e-16 Ha, dn 1.4e-16 / 5.5e-16, 7 vs 7 SCF steps (bound 1e-10 Ha, 1e-9, +-1). (b) standard 5e-9 Ha vs strict reference: ferro fixed sector (He nl + H, skewed, N_up 2 N_down 1) dE 2.32e-10 Ha, dM 2.2e-16 mu_B, SCF 15 vs 45; scalar FD metal (H, skewed, 2x2x2, theta 0.01 Ha, 12 bands) dF 1.98e-12 Ha, dmu 8.85e-7 Ha, sum w|df| 3.83e-6 e, SCF 5 vs 9; collinear FD free moment (H, seed 0.5) dF 5.50e-11 Ha, dM 6.70e-6 mu_B, dmu 6.64e-6 Ha, sum w|df| 6.70e-6 e, SCF 12 vs 39 (bounds 5e-9 Ha, 1e-4 mu_B, 1e-4 Ha, 1e-4 e); one verification attempt each, last E_H 8.9e-11 / 2.1e-13 / 1.2e-11 Ha. (c) forced failed verification (lib test): spin fixed and spin FD 2 attempts, +1 SCF step, dE 4.2e-12 / 8.9e-11 Ha; scalar FD 2 attempts; exhausted budget not converged. No wall-time claim (synthetic sizes) | the magnetization term is an estimator without an energy-bound proof; first-order quantities (mu, M, f) differ by up to 7e-6 at the 5e-9 class, consistent with sqrt scaling; dev08 request `docs/performance/requests/onepath3-spin-metals.md` (real Al/Fe/Ni-like NC fixtures, delivered error and cost vs strict). Next stages: locking auto, runtime modes, engine default |

| M2-WORKSPACE-SELECTED-ROOTS-DEV08 | 2026-11-16 | dev08 (master), 69e13f1, MPI vs MPI (no QE column), si32/si64 x standard/strict x np 2/4/8, f1e-3 seeds 601-603, arms interleaved per seed, med_by_seed, INVALID 0 | A = `M2_WORKSPACE=off`, B = default (workspace on), C = B + `M2_SELECTED_ROOTS=on`. B vs A: energy bits and density SHA identical at 12 of 12 points, wall -4.34..-12.50 %. C vs B: strict dE = 0 and SHA identical (3/3 at every point), wall within +-0.7 points of B except si32 strict np4 (-2.2 points); standard dE 3e-14..3e-13 Ha (class 5e-9), SHA differs (rounding class), wall same or faster at all six cells (-0.04..-3 points); C repeated at np=4 gives the same density SHA 4 of 4. Standard-arm spreads are large (up to 33 % of the median at np=2): standard differences are inside the spread. Rotation FLOPs -51 % (Mac count). | Workspace stays default on. **Selected-roots adopted as default on (2026-11-16)** under the adoption rule (no correctness failure, never slower beyond noise, strict bit-identical); `BAANDI_M2_SELECTED_ROOTS=off` is the rollback. |
| ONEPATH-4 | User decision 2026-11-16 (strict leaves production; stage 4, the last one, ticket `tickets/ONEPATH-4-LAST-FALLBACKS.md`, ADR 0043 amendment ONEPATH-4, single path complete): no strict fallback remains. `scf.locking = "auto"` is an input error on the default path naming `BAANDI_REFERENCE_STRICT=on` (locking is hard deflation of the strict residual Davidson; the standard Ritz-change Davidson retests every root, so they cannot be combined without a new eigensolver); observe/agent runtime modes and `bands`, `dos`, `checkpoint`, `restart`, `explain` run the standard path with the verification (caller exclusion removed); checkpoint v1 carries `verification`, `options.verification`, `continuation/result.verification_attempts` (additive keys omitted at defaults, old files unchanged; exact restore validates attempts against the history's verification records); runtime policy rejects a hybrid eigensolver on a standard session (`StandardRequiresDavidson`); `nscf` clears the SCF profile (one residual solve); engine default flip: `ScfOptions::default()` = standard 5e-9 Ha + verification + Pulay 0.5, `ScfOptions::strict_reference()` = former default bit for bit, `from_reference_switch()` mirrors the switch; Python `Options()` = standard + verification, `reference_strict` mirrors the switch, `stress()` refuses a threshold above 5e-11 Ha (`StandardProfile::STRESS_THRESHOLD_HA`); every engine/CLI/MCP test and fixture that relied on the strict default asks for `strict_reference()` explicitly (141 sites in 63 files, no tolerance changed) | adopted (default; rollback `BAANDI_REFERENCE_STRICT=on`, engine `strict_reference()`); all pre-registered bounds met | Mac debug, synthetic. Engine: exact checkpoint after a forced failed verification resumes bit for bit (energy, density, history, 2 attempts without wall), dropped/passed-flag/mode-mismatch tampering rejects (History/History/Options). CLI: observe mode vs classic (ladder off) identical energy bits, density hash and attempt count; agent mode (ignored socket test, passes on this Mac) verified; standard vs strict reference in observe mode dE 8.0e-12 Ha (bound 5e-9); checkpoint/restart with verification bit-identical to the uninterrupted run with REUSE_PSI off, dE 1.1e-16 Ha with REUSE_PSI on (bound 1e-10 Ha, pre-registered); verification option/attempts round-trip bit for bit, strict files gain no key; bands/dos SCF verified and recorded; locking auto rejects by default and runs under the switch. Python (installed debug wheel, 9/9 tests): native, CLI (ladder off) and ASE bit-identical on the default path; default vs reference within 5e-9 Ha; evidence energies moved 1.1e-12 to 7.6e-11 Ha from the former strict values; ASE stress oracle at 5e-11 Ha max FD error 8.8e-8 eV/Angstrom^3 (bound 1e-5); ASE BFGS vs CLI relax at 5e-13 Ha: positions 3.9e-6 Bohr, dE 3.4e-14 Ha (bounds 1e-4, 1e-10). No wall-time claim | restart exactness is defined against the checkpointed session (it forgets kept REUSE_PSI state), so a run that never checkpointed matches bit for bit only with REUSE_PSI off; runtime tolerance interventions act on strict-only options and do not move the standard stop rule (documented limitation); Python `stress()` at the 5e-9 default would leave the stress class, so it refuses rather than weakening. dev08 request `docs/performance/requests/onepath4-last-fallbacks.md` (locking, runtime modes, checkpoint exactness, bands, Python default, strict-reference identity) |
| ONEPATH-3B | dev08 finding 2026-11-16 (006baaf): Si8 Gamma PBE zero-moment collinear 22 SCF / 2 attempts / 6.33 s vs unpolarized 9 / 1 / 2.67 s (ticket `tickets/ONEPATH-3B-SPIN-COST.md`, ADR 0042 amendment ONEPATH-3B). Root cause: the spin driver never used the atomic start (construct recorded "spin-polarized occupations", the collinear map started cold solves from free-electron plane waves); the ONEPATH-3 Mac fixtures had no PP_PSWFC so both sides fell back there. Mixing (per-channel Pulay is exactly half the closed-shell step at equal channels), estimator, ratchet, redo, band count, preconditioner and spin XC ruled out by the free-electron control (equal SCF counts). Fix: cold fixed-sector solves build the atomic span once per k and run Rayleigh-Ritz per channel (closed-shell construction); FD (scalar and collinear) keeps free electrons with reason "Fermi-Dirac occupations" | adopted (default; `scf.start_wavefunctions = "free_electron"` restores the old start); unpolarized code path untouched | Mac release, pure, T=1, PseudoDojo Si.upf (sha256 3982...5282), 30 Ry: spin 18 SCF / 3 attempts / 3.20 s -> 9 / 1 / 1.58 s, E -33.47297201392545 vs unpolarized complex -33.47297201392599 (|dE| 5.4e-13 Ha, 9 / 1 / 0.74 s); default unpolarized (real Gamma + ladder) 6 / 2 / 0.47 s. Synthetic dimer (s+p orbitals): free-electron start dE <= 8.3e-16 Ha; atomic start SCF 7 vs 7, attempts equal, dE 5.8e-10 (LDA-PW) / 2.6e-10 Ha (PBE): pre-registered 1e-10 bound NOT met there, inside the closed-shell atomic-start rounding envelope 1.5e-9 / 3.0e-10 Ha (1e-13 fractional displacements) | an atomic start amplifies rounding-level potential differences at the 5e-3 Ha first tolerance (closed shell too); spin wall stays about 2x unpolarized (two channel solves; an equal-channel shortcut is a separate proposal); dev08 request `docs/performance/requests/onepath3b-spin-cost.md` |
| ONEPATH-2B | dev08 bug (binary 7a58d83, real Si8 PBE, relax mode cell isotropic and ions_cell): standard arms failed at the first step that expanded the cell with Numerical("Miller index [0, 0, 18] is outside the canonical frequencies for [36, 36, 36]") in both attempts; strict converged. Root cause (backtrace on the Mac, real PseudoDojo Si PBE, Si8, 30 Ry): `CoulombMetric::new`, the Hartree-metric estimator of the standard stop rule (built for every standard problem, never on the strict path) and of the opt-in Hartree mixing metric, built its density sphere from the current cell, ignoring the frozen Miller reference (`strain_reference`) on the frozen grid; after a cell expansion the sphere left the grid. Not the ladder, REUSE_PSI, POTENTIAL_REUSE, LEAN_PROJ_CACHE, APPLY_LEAN or real Gamma (each off: failure unchanged). Fix: the metric uses the problem's own density list (`options.basis(cell, 0, 1).density_basis(cell)`), identical sphere when no reference (ordinary runs bit-identical by construction); three call sites | fixed (no tolerance, threshold or strict-path change) | Regression `crates/dft-engine/tests/frozen_cell_metric.rs` (synthetic, cell +15 % on the frozen [9,9,9] grid, own grid [12,12,12]): fails before the fix with "Miller index [-5, 0, 0] ... [9, 9, 9]", passes after; standard vs strict dE 9.0e-11 Ha. Real Si8 PBE (Mac release, lattice x 0.98, mode cell isotropic, bfgs, force_tol 1e-3, stress default 5e-11): before the fix the standard arm failed at step 3 like dev08; after: converged in 8 steps like strict, final cell 3.2e-8 Bohr and positions 2.4e-8 Bohr from strict, SCF iterations 73 vs 50 strict, wall 11.2-13.2 vs 13.5-21.7 s (noise level). FINDING (separate, path-independent): Si8 ions_cell with one atom displaced (lattice x 1.01, displacement 0.005 isotropic or 0.01 full cell, force_tol 1e-3) fails on the STRICT reference too (step 20 "SCF unconverged", step 23 "density residual exceeded growth limit"), and the standard arm follows the same trajectory to the same failure (24 steps, no Miller-index reason); trajectory energy -33.530 Ha vs the relaxed diamond -33.494 Ha, i.e. the walk leaves the diamond basin before the SCF fails; not investigated here | Also: `cell_relaxation::combined_scf_enthalpy_derivative` (analytic vs FD 1e-7) failed after main ONEPATH-4 made `run::options` standard (2.5e-5 at 5e-9 Ha, 2.4e-7 at 5e-13: stress is first order in the residual error); it now uses the strict reference options explicitly (5.0e-10, bound unchanged). Lesson: the standard stop-rule metric is a density-space object and must follow the same frozen basis protocol as the density; any new density-space object must use `ScfOptions::basis` rather than the current cell. dev08: `docs/performance/requests/onepath2b-cell-relax.md`. Open: the ions_cell divergence on real Si8 (strict and standard) needs its own ticket |
| MPI-FS-PANEL-GS | 2026-11-16 MPI fast screen (stage 1, docs/performance/profile-2026-11-16-mpi.md): `BAANDI_M2_PANEL_GS=classic|pythag|fused` (default classic) in the distributed panel Gram-Schmidt (`distributed_davidson::panel_column`, both the workspace and the allocating path). pythag: the norm of the once-projected column rides with the second dot reduction, final norm by the Pythagorean update |v1|^2 - |s2|^2 (explicit norm reduction when the update cancels more than a factor 4); fused: the input norm rides with the first dot reduction and the second pass runs only when the per-column DGKS rule of `BAANDI_ORTHO` (keep 0.1, the rule the serial solver applies against the basis) fails, then as pythag. Per-column dependence decisions (`DEPENDENCE_TOLERANCE`) kept; no Cholesky (A52, MPI-M1-GATE1 respected) | CANDIDATE (stage 2 requested: docs/performance/requests/mpi-fast-screen-2026-11-16.md); pythag dominated by fused (not proposed) | (C/O) si32 production np=4 simulated, same binary, 3 alternating pairs, load 6.0-6.4: tree_sum calls per rank 17955 -> 7172 (-60.1 %), panel_dot 10968 -> 5484 (every column with k > 0 needed one pass: the DGKS rule never failed), panel_norm 5824 -> 525, tree_sum payload -3 %; SCF/Davidson/H 13/74/2884 both arms; energy bits equal, density hash differs (rounding class), reproducible; wall 4.954 -> 4.505 s (-9.1 %, spread 0.25 s), orthonormalize stage 0.661 -> 0.450 s, tree_sum seconds 0.459 -> 0.247 s. np=8 (1 pair): 17955 -> 7056 (-60.7 %), orthonormalize 0.922 -> 0.522 s, wall 5.717 -> 4.836 s. pythag (2 pairs): 17955 -> 12471 (-30.5 %), wall neutral, dE 1.4e-14 Ha. Engine SCF test (cubic and skewed, strict and standard, P=2/3): classic/pythag/fused dE vs threads-only <= 7.8e-16 Ha strict, <= 8.8e-11 Ha standard, ranks bit-equal, reproducible; tree sums e.g. 1104/938/796 | The in-panel dot+norm reductions (88 % of all solver reductions, M2-REDUCTION-COUNTS-F08) are latency bound; folding the norm into the dot reduction and skipping the second pass when DGKS holds removes 60 % of all tree sums without changing any work count. The serial solver's DGKS rule transfers to the in-panel pass because the in-panel coefficients of Davidson corrections are small. Mac wall is an observation (simulated ranks do not model MPI latency); dev08 decides |
| MPI-FS-TREE-RD | 2026-11-16: `BAANDI_MPI_TREE=binomial|rd` (default binomial; read by the CLI at MPI start, `dft_mpi::set_tree_sum`). rd = recursive doubling with `MPI_Sendrecv` (partner rank XOR 2^s, log2 P rounds, every rank ends with the sum) instead of the binomial reduce to rank 0 plus broadcast (2 log2 P rounds); power-of-two P only, other P keep binomial. Backlog MPI item (a), mpi-reduced-eigen-analysis.md item 3 | implemented, opt-in, FFI not runnable on the Mac: dev08 microbenchmark + bit check first (`examples/mpi_tree_bench.rs`), then stage 2 | (C) bit identity proven arithmetically and by the unit test `recursive_doubling_equals_the_binomial_tree` (P = 1, 2, 4, 8, 16, values where association matters): every rank's result equals `tree_associate` bit for bit (IEEE addition is commutative; each round adds the lower block's partial and the higher block's). Counter: message rounds per tree sum 2 log2 P -> log2 P (np=8: 6 -> 3); call counts unchanged. FFI diff: one extern declaration (`MPI_Sendrecv`) and `doubling_steps` in crates/dft-mpi/src/ffi.rs | Not screenable on simulated ranks (a simulated tree sum is two thread barriers whatever the pattern). The bit identity means the switch can be chosen by rank count without violating the ADR 0041 rule that one rank count has one fixed reduction function |
| MPI-FS-DENSITY-FUSED | 2026-11-16: `BAANDI_M2_DENSITY=bands|fused` (default bands). On the standard path the M2 density accumulates sum_b f_b psi_b^2 / volume on each rank's planes in ascending band order (`DistributedRealHamiltonian::band_density_sum`, sink form of `band_planes`) instead of writing the bands x plane-points per-band buffer and summing it in a second pass; the per-band buffer only serves the certificate screen, which runs on the strict path only (fused is not taken there) | NEUTRAL at stage 1 (bit-identical; rider arm D of the stage-2 request) | (C) bit-identical: same per-point sum order (pairs ascending, column 2t before 2t+1), test `distributed_band_values_match_serial_and_band_densities` (fused sum equals the ascending-band sum of `band_densities` bit for bit, ranks 1-5, three cell shapes). (O) si32 production np=4, 3 pairs: energy bits and density hash equal, work counts identical, m2 density stage 0.283 -> 0.241 s (-15 %), op.band_densities 0.222 -> 0.215 s, wall 5.045 -> 4.721 s (inside the 0.54 s spread). Memory: removes bands x plane points doubles per rank (si64: about 132 x grid/P values) | Only a memory-traffic item; worth measuring on dev08 where 8 ranks share one socket's bandwidth (MPI-SCALING-BEYOND-8 mechanism), not on the Mac |
| MPI-FS-VERIFY-COST | 2026-11-16 (coordinator request: verification cost under MPI): harness knob `M2_PROFILE_VERIFY=off` (sizing only) | finding | (O/C) si32 production np=4 (cap 3): verification on -> off SCF 13 -> 10, Davidson 74 -> 62, H 2884 -> 2322, reduced solves 142 -> 110, wall 4.66 -> 3.16 s (-32 %), dE 1.97e-9 Ha; same at cap 6: MPI np=4 12/73/2699 -> 10/65/2282, serial T=1 12/73/2087 -> 10/65/1670 (wall 7.01 -> 5.04 s). The verification is distributed exactly like an SCF map (Davidson, density and XC distributed; Hartree, Coulomb metric and density distance replicated); its cost is the failed attempts (two extra maps plus the passing one), identical on the serial and the MPI path | Not an MPI distribution problem: the verification costs the same share on the serial path. The lever is the standard stop rule (estimator calibration, decision 3 of 2026-11-11) or a cheaper verification map, not MPI engineering. Belongs to the single-path work order |
| MPI-FS-KEPT-START | 2026-11-16 finding from the stage-1 counters | implemented as MPI-KEPT-START (next row; opt-in BAANDI_M2_REUSE) | (C) at equal SCF/Davidson counts the MPI path applies H to 29 % more columns than the serial path (si32 production cap 6: 2699 vs 2087; cap 3 at np=4: 2884): the serial REUSE-PSI kept start (C2: H_new X = H_old X + (V_new - V_old) X from the kept psi(r)) and the lean solver do not exist on the M2 path, which re-applies H to the warm start every solve and gathers N x B to the replicated interface at every return | Next structural MPI item (review F07 plus C2 on the stick rows): keep X and H X in local rows across SCF steps (no N x B allgatherv, no extract_local), take psi(r) on the planes from the density build, and form (V_new - V_old) psi with one batched forward transpose. Expected (E): -20 % H columns, about -5 to -8 % wall at np=4 |
| MPI-KEPT-START | 2026-11-16, ticket MPI-KEPT-START, design docs/design/mpi-kept-start.md: `BAANDI_M2_REUSE=on|off` (default off until dev08; rounding class). Standard M2 solves keep local X (filled by the solver: no `extract_local`) and H X (return rotation of the H images, `DistributedKeptStart`); the M2 density build keeps the plane slabs of the leading fully occupied pairs (`band_densities_keeping`/`band_density_sum_keeping`, replicated memory gate 2 GiB on pairs x ceil(grid/P)); the next standard solve of the same operator adds (V_new - V_old) X (`apply_local_change`: kept pairs r -> G only, others both ways, no kinetic/nonlocal) and starts from basis X, products H X (no orthonormalization, no initial apply) when one tree sum finds initial = X to 1e-14 on every rank and the Gram defect of X <= 1e-12 (restart contraction rule). Invalidation: operator rebuild, strict/certified solve, failed solve, `forget_reuse` (exact checkpoint), non-workspace solver | implemented, opt-in; dev08 stage 2 requested (docs/performance/requests/mpi-kept-start.md) | Stage 1 (C, exact counters; Mac pure, simulated ranks, production = standard + verification + ladder + capacity rule, f1e-3 seed 601, one binary sha256 efe65d9b...): si32 np=4: SCF/Davidson 13/74 identical, H columns 2884 -> 2272 (-21.2 %), operator applied columns 5670 -> 4242 (-25.2 %), inverse pair transforms 3687 -> 3021 (-18.1 %), forward 2855 -> 2957 (+3.6 %, unkept pairs of the potential change), tree sums 17955 -> 13716 (-23.6 %: the warm-start Gram-Schmidt and the initial apply's projector sums), alltoall bytes -4.3 %, allgatherv bytes unchanged (N x B gather at return kept, F07 remainder), 24 armed / 21 taken / 3 declined; si32 np=2 2699 -> 2087 (-22.7 %), np=8 2884 -> 2272; si64 np=4 4353 -> 3165 (-27.3 %), inverse -20.3 %, tree sums -26.7 %, 15/15 taken. |dE| 0 to 5.7e-14 Ha, fixed-np density hash reproducible. Wall (O): si32 np=4 3 pairs at load 9, medians 6.71 -> 4.75 s with 2.7 s spread (CANDIDATE by counters); the single-pair si64 np=4, si32 np=2/8 walls overlapped a concurrent `check.sh` build (load 9-21) and are void; clean rerun si64 np=4, 3 pairs at load 7.6-8.7: medians 15.32 -> 14.37 s (-6.2 %, B faster in 3 of 3 pairs, spread 1.9 s), apply stage 5.08 -> 3.69 s per rank. Unit gates: kept-slab batch bit-identical to the plain batch (1-5 ranks); H_old X + dV X vs H_new X <= 2.4e-16 relative (gate 1e-13); solver test saves exactly `bands` columns at 1-4 ranks, perturbed initial (1e-10) declined; SCF test (cubic/skewed, P=2/3) dE <= 2.1e-14 Ha, SCF/Davidson identical, H -34 %, bit-reproducible | The serial lesson (C2 on dev08: H -25..-29 % gave wall -3..-6 %) bounds the expectation; on MPI the removed tree sums (-24 %) and the warm-start orthonormalization add latency savings that grow with np. F07 is only partly absorbed: the return allgatherv and complex unpack remain because the engine keeps complex full orbitals (warm start, checkpoint, prolongation); lazy materialization is a separate engine-wide ticket |

## Fast-loop screens (2026-11-16; threads/one-core owner, Mac observations O, `scripts/fast_screen.py`, si32-gamma f1e-3 seed 601 T=1 unless stated, main 62c821f pure build; profile docs/performance/profile-2026-11-16-threads.md; dev08 request docs/performance/requests/verify-reuse.md)

| ID | Idea | Status | Measured (Mac O, counters exact) | Lesson / next |
|---|---|---|---|---|
| FAST-SCREEN-HARNESS | `scripts/fast_screen.py`: one case, env/binary/input A/B, alternating pairs, counters with the all-zero self-check, verdict REJECT/NEUTRAL/CANDIDATE | adopted (tool) | one si32 T=1 run 5 to 6 s on the Mac, a three-pair screen under a minute | Pick a seed with failing verifications (603, 605, 610) when the idea touches the verification |
| SERIAL-CAPACITY-SCREEN | Serial lean path `scf.davidson_capacity_per_band` 3 and 4 against the default 6 (reopens A13/P31A with the lean premise, MPI-CAPACITY-RULE "serial open question") | REJECT (stage 1) | cap3: H 1918 -> 2244 (+17 %), verification attempts 1 -> 2, Davidson 53 -> 57, GEMM +9.8 %, wall +6.9 % (spread 1.8 %), dE 3.6e-11 Ha; cap4: H +13.9 %, attempts 1 -> 2, Davidson +17 %, wall +9.0 % (spread 1.5 %) | On one core the smaller basis costs more Davidson iterations and pushes the verification over the threshold; the MPI gain came from the replicated reduced eigenproblem, which the serial path does not have. Do not retry on the serial path without a new premise |
| C3a-SCREEN | `BAANDI_LEAN_PROJ_CACHE=on` (C3a, existing request) | CANDIDATE (stage 1) | GEMM 64.91 -> 62.20 GFLOP (-4.2 %), all other counters and the density SHA identical, energy identical to the printed digits; wall median ratio 0.957 (spread 2.5 %, 3 pairs) | Matches the profile share of `expectation_values` (2.1 %) plus GEMM; already pending dev08 under C3a-LEAN-PROJ-CACHE, add it to the stage-2 list of verify-reuse.md |
| PROJECT-WIDE-SCREEN | `BAANDI_PROJECT_WIDE=on` | not screenable on the Mac | pure build: bit-identical by construction (counters, energy, SHA identical; wall ratio 0.964 inside the 6 % spread) | The fallback it removes exists only on MKL builds (gemm-dispatch-analysis-2026-11-16.md); dev08 only |
| RITZ-FUSED-LEAN-SCREEN | `BAANDI_RITZ_FUSED=on` on the lean default path | REJECT | H vectors 1918 -> 2598 (+35 %), GEMM calls +61 %, GFLOP +11.5 %, wall +5.4 % (spread 4.3 %), dE 2.8e-14 Ha | Interaction: the fused return does not form H X, so REUSE-PSI keeps nothing and every solve pays its initial apply. Never stack RITZ_FUSED with REUSE_PSI |
| VERIFY-REUSE-WARM | `BAANDI_VERIFY_REUSE=warm` (new, default off): after a FAILED verification the next step starts from the verification solve's orbitals (so its REUSE-PSI kept start is taken) | NEUTRAL (stage 1), kept opt-in | seed 603 (4 attempts): H -1.5 %, Davidson +5.6 %, GEMM +3.5 %, wall ratio 0.951 inside an 11 % spread, dE 1.0e-11 Ha | The saved initial apply is offset by the warm start of a different Hamiltonian; the gain needs the mixer to use the map (feed) |
| VERIFY-REUSE-FEED | `BAANDI_VERIFY_REUSE=feed` (new, default off): warm plus the pair (rho_out, F(rho_out)) of a failed verification enters the mixer after the step's own pair, so the verification map counts as an SCF map instead of being discarded. Only acts after a failed attempt (bit-identical otherwise: si32 seeds 601/602/604/606/607, si-k444s) | CANDIDATE (stage 1), mixed in stage-2-like widening; dev08 stage 2 requested | si32 seed 603: H 2984 -> 2465 (-17.4 %), attempts 4 -> 3, fine SCF 13 -> 12, GEMM -11.9 %, wall ratio 0.849 (spread 5.1 %); seed 605: H -15.7 %, attempts 3 -> 2, wall 0.913 (spread 1.1 %); seed 610: H 3328 -> 2546 (-23.5 %), attempts 5 -> 3, wall 0.853 (spread 9.7 %); seeds 608/609 (3 and 2 attempts both arms): H -2.8 / -1.7 %, Davidson +1 / +3 %; si16 seeds 601/602/605: H +13 / +16 / +12 %, attempts 2 -> 3 (seed 601: second verification 5.09e-9 against 5e-9, a near-threshold miss), wall seed 601 +4.5 % (3 pairs, spread 4.1 %); si16 603/604 neutral; si-skew H -20 %, attempts 3 -> 2; si8 seed 603 H -15.5 %, attempts 3 -> 2; si8 601 -1.8 %. Delivered error against the strict reference (`BAANDI_REFERENCE_STRICT=on`, same binary): feed 1.4e-12 / 3.8e-11 / 3.3e-12 / 1.4e-10 Ha (si32 603, si32 610, si16 601, si8 603) against off 4.5e-11 / 4.8e-11 / 7.0e-11 / 1.3e-10 Ha; all far inside 5e-9. Engine test `verification_tests::failed_verification_reuse_keeps_the_fixed_point` (forced failures; off = default bit for bit; warm and feed within 1.4e-11 Ha) | Failing verifications are the per-seed variance driver (si32 seeds: 1 to 5 attempts; the stop estimator was 10 to 20 times below the verification residual at every failing stop). Feed removes most of the waste on long tails, but changes which near-threshold stop passes, so small cells can lose one attempt: stage 2 on dev08 with many seeds decides (med_by_seed) |
| VERIFY-TAU-QUESTION | Cost of a PASSED verification (si32 8 %, si64 16 % R): 5 Davidson iterations at tau_v = ratchet(tau, E, N_e) | not implemented (needs a user decision) | si32: tau_v 8.4e-14 Ha against the last step tolerance 5.9e-12 Ha; the Ritz-change rule needs at least two iterations and a full correction block in the first | tau_v is fixed by ADR 0043 V1; any looser tau_v is an ADR change, not a performance switch. Listed as a question in verify-reuse.md |
| VERIFY-FORESIGHT-ANALYSIS | Can the standard stop rule foresee a failing verification cheaply (coordinator question 2026-11-16)? Offline over the Mac arm-A runs of this session (si32 seeds 601-610, si16 601-605, si8 601/603, si-skew, si-k444s; 55 verification attempts, 26 passed, 29 failed) | analysis only, nothing implemented | Ratio verification residual / stop estimator at the stop: passed median 19.4 (range 0.1 to 30.4), failed median 21.1 (6.1 to 140.7); density contraction of the stop step: passed median 0.34 (0.15 to 0.71), failed 0.44 (0.19 to 0.82). The ratio is about 20 in both groups, so failures are stops whose estimator lies between about threshold/20 and threshold; contraction separates the groups only weakly | A foresight gate (verify only when the estimator is below threshold/k) is the same as a tighter stop threshold: it trades a failed attempt for an ordinary step and needs a tuned k. With VERIFY-REUSE=feed a failed attempt is already an ordinary map (its output enters the mixer), so feed removes the waste without a tuned constant; no gate proposed. tau_v stays as ADR 0043 fixes it |

| MPI-TREE-RD-DEV08-STAGE0-ADOPT | 2026-11-16 | dev08 (master), eab020c, `mpi_smoke` P=2,3,4,8 OK; `mpi_tree_bench` microbenchmark (worst rank, us per call), binomial -> recursive doubling | Bit check enforced at P=2,4,8 (OK), skipped at P=3 (not a power of two, falls back to binomial). P=2: payload 1 double 0.48 -> 0.29 (-39.6 %), 64: 1.20 -> 0.70 (-41.7 %), 27744: 78.15 -> 42.69 (-45.4 %); P=4: -37.0 / -31.6 / -44.9 %; P=8: 1.43 -> 1.03 (-28.0 %), 3.01 -> 2.47 (-17.9 %), 215.42 -> 133.80 (-37.9 %); P=3 identical (fallback). First dev08 run of the new FFI (`MPI_Sendrecv`, `doubling_steps`). | **`BAANDI_MPI_TREE=rd` adopted as default (2026-11-16)** (bit-identical for power-of-two ranks, never slower, fixed-order reproducibility kept); rollback `BAANDI_MPI_TREE=binomial`. Whole-solve effect to be read in the MPI stage-1 screen (arm C vs A). |

| VERIFY-REUSE-FEED-DEV08-STAGE1 | 2026-11-16 | dev08 (master), 82209bb MKL, T=1 core 0, 3 alternating pairs, stage 1 (one case per row) | `BAANDI_VERIFY_REUSE=feed` vs off. si32 s603 (A has 4 verification attempts): wall -11.12 %, H columns 2984 -> 2466 (-17.4 %), attempts 4 -> 3, SCF 13 -> 12, GEMM GFLOP -11.8 %, dE 4.3e-11 Ha. si32 s601 (A: 1 attempt): all nine counters identical, dE = 0, SHA identical (wall +0.48 % inside the spread: NEUTRAL). si16 s601/s602/s605 (A: 2 attempts): REGRESSION +10.31 / +11.93 / +9.88 % wall, attempts 2 -> 3, SCF 10 -> 11, H +13.0 / +15.9 / +12.3 %. All runs converge, |dE| <= 6.6e-11 Ha. The seeds were chosen where the agent had seen an effect (selection bias): stage 2 over seeds 601-610 decides. | `feed` only pays when failures repeat; with exactly two attempts it creates an extra attempt. Not adopted; the agent explores variants (feed from the 2nd failed attempt on, rule from the measured ratio) under the same switch. |
| LEAN-PROJ-CACHE-DEV08-STAGE1 | 2026-11-16 | dev08, 82209bb MKL, si32 T=1 | `BAANDI_LEAN_PROJ_CACHE` off vs on: wall -0.91 % (A 8.246 s, spread 0.014; B 8.171 s, spread 0.004), all counters identical, dE 0, SHA identical. Same order as the A-vs-A floor (0.2-3.9 %). | Not decidable at SCF level; re-measure at the diagonalization level (stage 0, FASTBENCH-DAVIDSON). |
| VERIFY-REUSE-FEED-DEV08-S1 | dev08 stage 1 of `feed` (relayed 2026-11-16, coordinator; 82209bb MKL, T=1, 3 alternating pairs) | superseded by VERIFY-REUSE-HALF | si32 s603 -11.12 % wall, H -17.4 %, attempts 4 -> 3, SCF 13 -> 12; si32 s601 identical; si16 s601/s602/s605 +10.3 / +11.9 / +9.9 %, attempts 2 -> 3, SCF 10 -> 11, H +12..16 % | Full feed only pays when failures repeat; after exactly one failure it causes an extra attempt |
| VERIFY-REUSE-VARIANTS | Variants behind the same switch (Mac O, counters, one pair per geometry; si16 seeds 601/602/605 (2 attempts each on off) and si32 seeds 603/605/610 (4/3/5 attempts)); totals over the six geometries, off -> variant | feed2, replace REJECT; half CANDIDATE | `feed`: attempts 18 -> 17, SCF 69 -> 68, H -10.1 % (si16 +12..16 % each). `feed2` (first failure keeps the default continuation, repeated failures feed): si16 identical, si32 603 -14.4 %, 605 0 %, 610 -14.6 %; total attempts 18 -> 16, SCF 69 -> 67, H -7.5 %. `replace` (only the verification pair enters the mixer, the step's pair is not recorded): attempts 18 -> 31, SCF 69 -> 84, H +39.3 % (si32 610: 5 -> 9 attempts), REJECT. `half` (both pairs enter the mixer, the next input is the midpoint of the step's own Pulay proposal and the fed proposal): si16 601/602/605 attempts 2 -> 2, SCF 10 -> 10, H -1.2 / 0.0 / +0.6 %; si32 603 4 -> 2 attempts, SCF 13 -> 11, H -26.7 %; 605 3 -> 2, 12 -> 11, -15.8 %; 610 5 -> 2, 14 -> 11, -35.0 %; total attempts 18 -> 12, SCF 69 -> 63, H -19.4 % | A ratio rule (feed only when the attempt's own verification/threshold ratio is in a range) was not built: the first failing attempts of the regressing si16 seeds (ratio 3.1) and of the winning si32 seeds (2.1) do not separate, so any cut would be a tuned constant. The midpoint is the unweighted mean of two proposals (no fitted factor) |
| VERIFY-REUSE-HALF | `BAANDI_VERIFY_REUSE=half` widened (Mac O): si32 seeds 601/602/604/606/607 (1 attempt) bit-identical; 608 attempts 3 -> 2, H -16.5 %; 609 H -3.1 %; si16 603/604 H -0.2 / 0.0 %; si8 601 -1.8 %, si8 603 3 -> 2 attempts, H -15.5 %; si-skew H -8.7 %; si-k444s identical. Total over these 13 geometries H -2.7 %, attempts 23 -> 21, never an extra attempt or SCF step in 19 geometries. Wall (3 pairs): si32 610 ratio 0.755 (spread 9 %, CANDIDATE), si16 601 0.967 (spread 8.5 %, neutral). Delivered error against the strict reference: si32 610 1.6e-10 Ha, 605 2.1e-10 Ha (gate 5e-9); against off <= 1.9e-10 Ha everywhere | CANDIDATE for dev08 stage 2 (the one variant picked); engine test covers all modes (forced failures, within 1.4e-11 Ha of the default fixed point) | Damping the fed step keeps the first-failure trajectory close to the default (no extra attempt on si16) while still using the verification map on long tails |
## Stage-0 Davidson screens (FASTBENCH-DAVIDSON, 2026-11-16; Mac observations O, counters exact; `scripts/fast_screen_dav.py`, one solve at a fixed Hamiltonian, si32-gamma f1e-3 seed 601 T=1, base density = converged `dft-cli run` density.bin, cold = atomic start at 5e-3 Ha, warm = kept start at the ratchet tolerance 3.6e-11 Ha, 3 alternating pairs x 5 repeats, pure build, load 3 to 5; protocol docs/performance/fast-loop-protocol.md stage 0)

| ID | Idea | Status | Measured (solve-time ratio B/A median, spread; counters per solve) | Lesson / next |
|---|---|---|---|---|
| FASTBENCH-DAVIDSON | `crates/dft-cli/examples/davidson_bench.rs` + `scripts/fast_screen_dav.py` (+ `fast_screen_mpi.py --level solve`): cold and warm solve, counters first, eigenvalue identity | adopted (tool) | A/A: cold 0.992 (+-1.0 %), warm 0.989 (+-4.6 %), bit-identical; si32 screen 44 to 59 s, si8 2.5 s, si32 np=4 simulated 27 s; one si32 process: cold solve 0.26 s (3 Dav, 190 H columns), warm 0.36 s (10 Dav, 176 H) | Run every kernel idea here before stage 1. Engine refactor: atomic start moved into `atomic_guesses` (CLI density SHA identical si8, si32) |
| DAV0-LEAN-PROJ-CACHE | `BAANDI_LEAN_PROJ_CACHE=on` | CANDIDATE (stage 0; matches C3a-SCREEN stage 1) | cold 0.973 (+-1.4 %), warm 1.031 (+-4.6 %); map GEMM FLOP -5.2 % / -4.1 % (the audit's P^T X is read from the cache: `nl_audit.direct` 1 -> 0), Davidson GEMM +0.1 %, H and Davidson counts identical, eigenvalues bit-identical | The saving is in the map's audit, not the solve; pending dev08 C3a |
| DAV0-DAVIDSON-CLASSIC | `BAANDI_DAVIDSON=classic` against lean | REJECT (confirms the lean default) | cold 1.086, warm 1.392 (+-1 %); warm H columns 176 -> 244 (+39 %, no kept start), GEMM FLOP +54 %, GEMM calls x3.4; rounding class (max d eps 1.9e-15 Ha) | Do not repeat |
| DAV0-ORTHO-BLOCK | `BAANDI_ORTHO_BLOCK=off` (pure default on) | NEUTRAL on time; on-default confirmed by counters | cold GEMM calls 108 -> 44, FLOP -5.3 %, time 1.010 (+-6.7 %); warm identical counters, 0.995 | Off saves GEMM calls only in the cold solve; no time gain; default unchanged |
| DAV0-PURE-GEMM-NANO | `BAANDI_PURE_GEMM=nano` (pure default packed) | REJECT (confirms packed) | cold 1.20, warm 1.14 (+-2 to 4 %), same FLOP | Do not repeat on pure builds |
| DAV0-PROJECT-WIDE | `BAANDI_PROJECT_WIDE=on` | not screenable on the Mac (pure build bit-identical, 1.006 / 0.998) | counters identical | MKL builds only: run stage 0 on dev08 (`--features mkl`) |
| DAV0-NONLOCAL-BLOCK | `BAANDI_NONLOCAL_BLOCK=on` | not reachable (si8 Gamma bit-identical, counters identical) | the switch acts on the complex k-point nonlocal only | Stage 0 is Gamma-only; screen at stage 1 on a k mesh |
| DAV0-M2-PANEL-GS | `BAANDI_M2_PANEL_GS=fused` / `pythag`, si32 simulated np=4 | CANDIDATE (counters) | fused: tree_sum calls per solve 581 -> 223 (cold), 1476 -> 681 (warm), bytes -2 %, time 0.963 (+-4.8 %) / 0.958 (+-6.5 %); pythag: 581 -> 402, 1476 -> 1060, time 0.982 / 0.987 (+-1 %); rounding class (d eps <= 1.1e-14 Ha) | Latency not modelled on simulated ranks: time with real MPI on dev08 (stage-0 `mpirun` commands in the protocol) |
| DAV0-M2-SELECTED-ROOTS | `BAANDI_M2_SELECTED_ROOTS=off` (A) against the default on (B), np=4 | default confirmed | warm rotation FLOP 2.30e9 -> 4.93e8 (-79 %), map GEMM FLOP -49 %, solve time 0.797 (+-1 %); cold 0.993; eigenvalues bit-identical on this case | Do not repeat |
| DAV0-MPI-TREE-RD | `BAANDI_MPI_TREE=rd`, np=4 simulated, 1 pair | counters only (not screenable) | all counters and eigenvalues identical (the simulated tree sum ignores the pattern) | dev08 real MPI only |
| DAV0-MPI-WARM-GAP | finding from the counters | open | warm solve on np=4 (capacity 3, no kept start on M2): 23 Davidson iterations, 272 H columns against 10 / 176 on the threads path at the same density | Supports MPI-FS-KEPT-START (kept start on the stick rows) and the capacity question; stage-0 measurable once implemented |

## MPI scaling gap (MPI-SCALING-GAP, 2026-11-17; Mac observations O on simulated ranks, counters exact C, estimates E on relayed dev08 numbers D; analysis docs/performance/mpi-scaling-gap.md, dev08 request docs/performance/requests/mpi-scaling-gap.md)

| ID | Idea | Status | Measured | Lesson / next |
|---|---|---|---|---|
| MPI-SCALING-GAP-MODEL | Why the M2 solve scales 2.2x from np 2 to 8 (D: si32 cold 235/141/109 ms, warm 428/264/200 ms) | finding | (E) Amdahl T = s + w/np: cold (same work at every np) s = 61 ms, w = 342 ms (15 % of the work, 56 % of the np=8 solve); warm np 4/8 (capacity 3) s = 136 ms (68 % of np=8). Split T = a H + b D + c: a = 0.63 / 0.47 ms per H column at np 4 / 8, b = 3.6 / 2.6 ms per Davidson iteration; a = 0.30 + 1.31/P ms, b = 1.7 + 7.4/P ms: at np=8 warm the per-column constant is 83 ms (41 %), the per-iteration constant 39 ms (19 %), per-solve interface about 12 ms. (O, stage 0 si32 capacity 3, np 2 -> 8): rotate 3.3-3.5x, apply 2.4x (transposes 2.6x), reduced_gram 1.8-2.0x, orthonormalize 0.95-1.03x, reduced_eigen 0.81-0.86x, solve_start 0.67x, mirror exchange 0.51-0.61x, solve_return 0.18x, pack+unpack 0.48-0.52x; stick/row imbalance at most 1 % (C). SCF level (O si32 production): replicated potential/audit/mixing 21-29 ms per fine step per rank, flat or growing (about 10 % of the np=4 wall); ladder coarse session flat from np 4 (0.96 -> 1.03 s) | The non-scaling part is mainly a constant per H column in the distributed apply (four alltoalls per apply with a node volume independent of np, plus about three in-panel reductions per appended column), then the replicated reduced eigenproblem; Part A of the dev08 request decides between transposes and reductions. New instrumentation: solver stages `solve_start`, `solve_return`, engine stages `m2.solve_prepare`, `m2.solve_unpack`, davidson_bench `m2.*` stages, m2_entry_cost SCF phases; `examples/mpi_transpose_bench.rs` (alltoall and N x B allgatherv per call at the solve payloads) |
| MPI-SCALING-GAP-STICKS | `BAANDI_M2_STICKS=greedy|mirror` (default greedy): stick (z, o) and its mirror ((-z) mod n2, (-o) mod n1) assigned to one rank as a unit in `distributed_fft::partition`; no ORC pair representative is split, `max_items` = 0, the two mirror alltoalls per apply and their pack/unpack vanish | CANDIDATE (stage 0 and stage 1; dev08 Part B/C) | (C/O) stage 0 si32 np=4, 3 pairs: alltoall calls 14 -> 7 (cold), 94 -> 47 (warm), bytes -28.5 / -27.6 %, solve time 0.977 / 0.981 (spread 1 %), max d eps 9.6e-15 Ha. SCF si32 np=4 (3 pairs): alltoall calls -50 %, mirror stage 0.068 -> 0 s, wall 4.365 -> 4.379 s (spread 0.25 s), dE 1.4e-14 Ha, reproducible. Stick imbalance (C): si32 np=8 1.010 -> 1.058, otherwise <= 1.010 | Rounding class (rows summed in another order). The cost is a slightly coarser balance unit (two sticks); a balance rule on stick count could recover it if dev08 shows the imbalance |
| MPI-SCALING-GAP-START | `BAANDI_M2_START=ortho|gram` (default ortho): workspace solver takes the start block as its basis when one Gram reduction finds a defect <= 1e-12 (restart contraction rule), else the panel Gram-Schmidt; counters `start_gram_stats` | CANDIDATE (stage 0 and stage 1; dev08 Part B/C) | (C/O) stage 0 si32 np=4: tree sums 581 -> 379 (cold, -35 %), 1476 -> 1274 (warm), solve time 0.944 / 0.971 (spread 1 %), max d eps 3.8e-14 Ha, every start taken. SCF si32 np=4: tree sums 17955 -> 12703 (-29.3 %), 26 of 26 starts taken (atomic, warm, coarse, prolonged, verification), solve_start 0.288 -> 0.040 s per rank, wall inside the 1.0 s spread (load 4.5-6.7), dE 5.7e-14 Ha | Every start the SCF hands to the M2 solver is already orthonormal (Ritz vectors, atomic Ritz vectors, zero-padded prolongation); the Gram-Schmidt chain only re-did that at about 200 reductions per solve |
| MPI-SCALING-GAP-BOTH | M2_STICKS=mirror + M2_START=gram, SCF level | CANDIDATE | (O) si32 np=8 production, 2 pairs: wall 4.623 -> 4.234 s (-8.4 %, spread 0.13 s), alltoall calls 620 -> 310, bytes -33.5 %, tree sums -29.3 %, energy bits equal, reproducible; si64 np=4, 2 pairs: 14.004 -> 13.471 s (-3.8 %, spread 0.47 s), alltoall bytes -33.0 %, tree sums -30.2 %, dE 5.7e-14 Ha; SCF / Davidson / H identical in every pair. Engine test (cubic, skewed, P 2/3): strict dE <= 5.2e-15 Ha, standard <= 3.7e-11 Ha against threads-only | Simulated ranks do not model MPI latency or bandwidth: dev08 Part B (stage 0, np 2/4/8) and Part C (SCF) decide |
| MPI-SCALING-GAP-CAP2 | `scf.davidson_capacity_per_band = 2` on the M2 path (new premise against A13 / SERIAL-CAPACITY-SCREEN: replicated reduced eigenproblem) | REJECT (stage 0) | (O/C) si32 np=4 / np=8 against capacity 3: warm Davidson 23 -> 19, H 272 -> 334 (+23 %), sum m^3 -47 %, reduced_eigen 23 -> 13 ms, solve 285 -> 298 ms (+4.5 %) / 264 -> 281 ms (+6.5 %); cold 135 -> 125 / 124 -> 117 ms (3 iterations only) | Extra H columns cost more than the halved eigenproblem saves; with the per-column constant on dev08 even more so. Do not repeat without a cheaper apply |
| MPI-SCALING-GAP-NOT-DONE | Third term (replicated reduced eigenproblem; replicated N x B interface and grid stages) | no cheap fix | (C) root and B5 measured (MPI-REDUCED-EIGEN-DEV08), dsyevr premise unchanged (order/bands <= 3), threading needs idle cores and the equal-core comparison binds np ranks to np cores (none idle at any np); the N x B gather feeds the audit kinetic energy, warm start, prolongation and checkpoint every step (engine-wide lazy orbitals, mpi-kept-start.md section 7); fusing the residual-norm and correction-norm reductions saves one reduction per iteration (about 0.2 ms per warm solve): not worth a switch | Next levers, sized in the analysis: distributed full-grid FFT for Hartree/audit/Coulomb metric (the `DistributedPairFft` plane API over the full box), lazy orbitals, distributed back-transform of the reduced eigenproblem (tridiagonalization stays replicated) |

## MPI lazy orbitals (MPI-LAZY-ORBITALS, 2026-11-17; Mac observations O on simulated ranks, counters exact C, relayed dev08 D, estimates E; design docs/design/mpi-lazy-orbitals.md, dev08 request docs/performance/requests/mpi-lazy-orbitals.md; review F07 second half)

| ID | Idea | Status | Measured | Lesson / next |
|---|---|---|---|---|
| MPI-LAZY-ORBITALS-FINDING | The N x B allgatherv at every M2 solve return grows with np (every rank receives the whole block) | finding (D, coordinator relay 2026-11-17, real MPI microbenchmark) | si32 (11877 x 68) 557 / 1434 / 2050 us per call at np 2 / 4 / 8 (+268 %), node volume 12.9 / 25.8 / 51.7 MB; si64 (23871 x 132) 4551 / 7823 / 9994 us, 50 / 101 / 202 MB; the alltoall transposes of the same solve scale (si32 34-pair block 1593 / 889 / 521 us) | The gather and the complex unpack sit in the np-independent part of the solve (Amdahl s of MPI-SCALING-GAP-MODEL); the consumer inventory (mpi-lazy-orbitals.md section 2) shows that nothing in the standard loop needs the full block |
| MPI-LAZY-ORBITALS | `BAANDI_M2_LAZY=on|off` (default off until dev08): standard M2 session solves return local rows only (`solve_distributed_rows_local`, `gather = false`, 0 x B placeholder plus a `LazyOrbitals` handle on the k point); the next solve starts from the rows (ORC round trip applied, so the start bits equal the gathered path's); audit kinetic energy from the rows (one tree sum of B values), nonlocal from the existing M2 projections; the full block is gathered once on demand (result, ladder hand-off, checkpoint, spectral and explicit audits, fallbacks) and cached | CANDIDATE (stage 0 and stage 1); dev08 stage 0 and SCF level requested | (C/O) stage 0 (one solve, 3 pairs): eigenvalue bits equal in every cell; allgatherv bytes per rank per map si32 11.04 -> 4.57 MB (-58.5 %, np 2/4/8), si64 np=4 34.36 -> 9.15 MB (-73.4 %); solve B/A cold/warm si32 np2 0.991/0.994 (neutral), np4 0.971/0.991, np8 0.927/0.971, si64 np4 0.969/0.988. Stage 1 (production = standard + verification + ladder + capacity rule, 3 pairs): density hash and SCF / Davidson / H equal in every pair, dE 1.4e-14 to 2.8e-14 Ha; allgatherv bytes si32 -54.1 % (224.4 -> 102.9 MB per rank and run), si64 -65.1 %; `solve_unpack` 26 -> 0 calls, `solve_prepare` -60 to -68 %; walls (O) si32 np4 -2.8 %, np8 -8.3 %, si64 np4 -5.1 % (inside the 0.75 s spread), si32 np4 with M2_REUSE on both arms -4.5 %. Engine test (cubic, skewed, P 2/3, reuse off/on): state bit-identical to off, dE <= 1.8e-15 Ha | The SCF trajectory is bit-identical (only the reported kinetic sum is rounding class), so the dev08 decision is pure wall. Remaining allgathervs per map: density planes and XC potential (full grid each; the distributed full-grid backlog item) and the solver's per-solve row-index exchange (N doubles, cacheable per operator). E on dev08 from the gather alone: si32 -17 / -43 / -62 ms per run (np 2/4/8), si64 -91 / -156 / -200 ms, plus unpack and pack |

| VERIFY-REUSE-HALF-DEV08-ADOPT | 2026-11-17 | dev08 (master), 6b78f06 MKL, T=1, stage 1 (si32/si16 seeds, 3 pairs) and stage 2 (arm F, 1 pair per seed alternating) | `BAANDI_VERIFY_REUSE=half` (keep a failed verification's output, mix it after the step's own result, start the next step from the midpoint of the two Pulay proposals). si16 10 seeds -2.98 % median, si32 10 seeds -3.01 % (T=8 -2.98 %), seeds with >= 3 failed attempts -10..-25 % (si32 603 -19.2 %, 605 -10.0 %, 608 -10.6 %, 610 -25.1 %), seeds with one attempt bit-identical. si64 seeds 601-610: median -16.29 %, summed wall -10.36 % (445.1 -> 399.0 s), attempts fell on 8/10 seeds, rose on 1/10 (602: 2 -> 3 attempts, +12.1 %), all converged and verified, max |dE| 1.28e-10 Ha. Pre-registered 4 criteria PASS (the first 3-seed si64 look, +0.47 %, was a small-sample artifact). | **Adopted as default (half), 2026-11-17**; rollback `BAANDI_VERIFY_REUSE=off`. |
| M2-LAZY-DEV08-STAGE0 | 2026-11-17 | dev08, 05e55d5 mkl,mpi, davidson_bench real MPI, si32/si64 x np 2/4/8, 3 pairs | `BAANDI_M2_LAZY` off vs on: eigenvalue fingerprint, H columns and Davidson iterations identical at 12 of 12 points; allgatherv bytes -58.6 % (si32, predicted -58.5 %) and -73.4 % (si64, predicted -73 %); solve time B/A grows better with ranks: si32 cold 0.992 / 0.924 / 0.816 (np 2/4/8), si64 cold 0.968 / 0.925 / 0.816, warm 0.980 / 0.949 / 0.903. The only MPI candidate with a double-digit gain at np=8. | Strong candidate; SCF-level real-MPI confirmation queued (seeds 601-603, strict cell), then default on. |

| M2-LAZY-PANELGS-DEV08-ADOPT | 2026-11-17 | dev08 (master), 05e55d5 mkl,mpi, SCF level real MPI, si32/si64 production np 2/4/8 seeds 601-603 (19 rows, INVALID 0) plus stage 0 | `BAANDI_M2_LAZY` off/on: density SHA identical 19/19, work counters identical 19/19, ranks equal 19/19, energy differs by at most 2.27e-13 Ha (about 1 ulp at 250 Ha; reduction order of the audit sum), strict cell (si64 np=8 `BAANDI_REFERENCE_STRICT=on`) allgatherv bytes ratio 1.0000, dE 0. Wall B vs A: si32 np2/4/8 -1.64 / -7.09 / -13.18 %, si64 -2.50 / -4.13 / -15.35 %; allgatherv bytes -52 % (si32) and -66 % (si64). (The harness verdict FAIL came from the requested 'energy bit-identical' clause and from requiring `verified` on strict rows, which carry no such key: harness defects, recorded.) `BAANDI_M2_PANEL_GS=fused`: stage 0 real MPI -3.8..-5.0 % solve time in 6 of 6 cells, tree_sum calls -55..-60 %, rounding class (eigenvalue bits differ). | **Both adopted as defaults (2026-11-17)** (1 ulp energy difference accepted as rounding class); rollbacks `BAANDI_M2_LAZY=off`, `BAANDI_M2_PANEL_GS=classic`. |

## MPI remaining terms (MPI-REMAINING, 2026-11-17; Mac observations O on simulated ranks, counters exact C, relayed dev08 D, estimates E; analysis docs/performance/mpi-remaining.md, dev08 request docs/performance/requests/mpi-remaining.md)

| ID | Idea | Status | Measured | Lesson / next |
|---|---|---|---|---|
| MPI-REMAINING-MAP | Per-stage scaling of the new default stack (lazy, fused, selected roots, verify half, rd tree, workspace) and the SCF phases outside the solve | finding | (O, stage 0, ms per rank, warm si32 np 4 -> 8) solve 209 -> 180 (1.16x), apply 93.5 -> 67.6, rotate 22.4 -> 11.9, orthonormalize 49.8 -> 48.0 (flat), reduced_eigen 18.5 -> 22.4 (grows), residual sync 4.2 -> 8.8 (grows), solve_return / solve_unpack 0 (lazy); non-shrinking group 104 of 180 ms (58 %). si64 warm np 4 -> 8: 1630 -> 1284, orthonormalize 536 -> 392, reduced_eigen 201 -> 219. (O, SCF si32 production np 2/4/8) wall 5.29 / 3.26 / 3.05 s; ladder coarse session 1.06 / 0.72 / 0.76 s (25 % at np=8, flat); verification 0.70 / 0.37 / 0.29 s; replicated fine potential + audit + mixing 0.21 / 0.20 / 0.27 s (grows; 9 % at np=8) plus the untimed stop-rule FFT pair per map | Biggest np-independent terms: ladder coarse stage, panel Gram-Schmidt chain, verification, replicated grid stages, replicated reduced eigenproblem. Against QE at one core the Davidson algebra (orthonormalize 7.5x, rotate 12x, eigen 2.9x) is the stage family furthest away, and it is the family that does not scale under MPI: next ticket = overlap-matrix (generalized reduced problem) Davidson on M2 (backlog). QE MPI timers requested (Part A) for an MPI-against-MPI stage table |
| MPI-REMAINING-MIX-INCREMENTAL | `BAANDI_MIX_INCREMENTAL=on|off` (default off): Pulay mixer keeps differences, metric images and Gram entries across calls (keyed by pair/seed serial numbers), forms only the newest difference and Gram row; reset/seed/restore drop the cache | CANDIDATE (stage 1; dev08 Part C arm I) | (C) bit-identical: unit test (window, checkpoint/restore, reset, seeds, l2 and metric, capacities 2/6), engine test threads-only and M2 P=2 (energy, density, counts equal). (O, si32 production, 3 pairs) differences formed 78 -> 18 (-77 %), dot products 476 -> 216 (-55 %), mixing phase np=4 0.054 -> 0.028 s (-48 %), np=8 0.083 -> 0.043 s; wall np=4 3.024 -> 2.991 s (-1.1 %, spread 0.12 s), np=8 2.921 -> 2.834 s (-3.0 %, spread 0.17 s); si64 np=4 (with STOP_METRIC kept) mixing 0.105 -> 0.054 s, differences 48 -> 13, dots 266 -> 131 | Overhead-only change on every parallel model (the mixer is replicated per rank and on the threads path); the saving per rank is np-independent, its share grows with np |
| MPI-REMAINING-STOP-KEPT | `BAANDI_STOP_METRIC=fft|kept` (default fft): the standard stop-rule estimate E_H[rho_out - rho_in] by Parseval from the map's own coefficients (rho_in kept by `hartree::evaluate_keeping_density`, rho_out moved out of the audit), used only when both densities are bit-equal to the kept ones and unclamped | CANDIDATE (stage 1; dev08 Part C arm K) | (C) estimate within 1e-12 relative of the transform form (unit test, skewed cell); SCF state bit-identical in every Mac cell and engine fixture (the estimate only gates the stop and the first ratchet). (O) si32 np=4/8: stop-rule FFT pairs 22 -> 0 per run (every estimate kept, no fallback), si64 np=4 17 -> 0; first build copied the audit coefficients (np=8 wall 3.223 -> 3.223 s, spread 0.27 s; potential/energy phases +10 to +15 ms under load 5), after moving them (5 pairs) np=4 3.327 -> 3.292 s (-1.1 %, spread 0.05 s), np=8 3.239 -> 3.186 s (-1.6 %, spread 0.13 s), bit-identical | Rounding class of the estimate only; strict runs have no standard stop rule; collinear and finite-T estimators unchanged |
| MPI-REMAINING-STACK | Combined stack S = `M2_STICKS=mirror M2_START=gram M2_REUSE=on` on the adopted defaults (fused panel Gram-Schmidt, lazy orbitals), with and without I and K | CANDIDATE (stage 0 and 1; dev08 Parts B and C) | (C) test `m2_remaining::combined_mpi_stack_keeps_the_standard_class` (cubic, skewed, P 2/3): |dE| vs threads-only 2.7e-15 to 1.6e-14 Ha, counts equal to threads-only, grid switches change the stacked energy by 0 Ha, ranks bit-equal, reproducible. (O, stage 1, 3 pairs) si32 np=4 3.318 -> 2.992 s (-9.8 %, spread 0.03 s), np=8 3.192 -> 2.960 s (-7.3 %, spread 0.20 s), si64 np=4 12.80 -> 11.88 s (-7.2 %, spread 0.35 s); H columns -26.5 % (si32) / -27.3 % (si64), alltoall calls -50 %, bytes -38 to -41 %, tree sums -25 to -29 %, solve_start 0.085 -> 0.028 s, dE <= 5.7e-14 Ha; S + I + K si32 np=8 3.183 -> 2.802 s (-12.0 %, spread 0.05 s). Stage 0 (si32, 3 pairs, cold / warm solve B/A): np=2 0.951 / 0.897, np=4 0.956 / 0.913, np=8 0.940 / 0.932; warm H 272 -> 204, alltoall calls halved, tree sums cold -32 %, warm -11 % | First test of all pending MPI switches together; correctness holds, the counters add up (the kept start and the mirror/gram switches act on different terms) |
| MPI-REMAINING-REDUCED-EIGEN | (b) root solve / B5 / partial EVD for the replicated reduced eigenproblem, premise re-checked | not done (premise unchanged) | (O) share grew after the lazy gather removal (si32 np=8 warm 12 %, si64 17 %), but MPI-REDUCED-EIGEN-DEV08 (D) measured root and B5 within +-1.7 % and the root solve does not shorten the critical path; B5 acts on 11 of 127 solves | Fewer reduced solves (overlap-matrix Davidson) or a distributed dense EVD are the only changed premises |
| MPI-REMAINING-ASYNC-ALLTOALL | (c) nonblocking alltoall overlapping the FFT transposes with the next pair batch | REJECT (analysis, no code) | (D) 34-pair si32 transpose 1593 / 889 / 521 us at np 2/4/8 (scales 3.1x); warm si32 solve 47 transposes of about 6 pairs, about 5 ms at np=8 against about 200 ms (2 to 3 %); the fft_transposes stage is mostly 1D FFT compute that already scales | Best case below 3 % of the solve for a new FFI surface (MPI_Ialltoall, request handles, unsafe); mirror sticks already halve the calls. Do not repeat without a profile showing the alltoall itself above about 10 % of the solve |

| M2-STICKS-START-REUSE-DEV08-ADOPT | 2026-11-17 | dev08 (master), e09067d mkl,mpi, SCF level real MPI, seeds 601-603, 24 + 36 rows, all converged, INVALID 0; stage 0 earlier (40 combinations) | `BAANDI_M2_STICKS=mirror` (M), `BAANDI_M2_START=gram` (G), `BAANDI_M2_REUSE=on`. Stage 0 solve ratios: M 0.94-0.99, G 0.93-0.98, MG 0.88-0.94, MGF (+fused panel) 0.82-0.92 in every cell (counters: alltoall calls exactly halved, bytes -27.6 %, tree_sum -25..-35 %). SCF level vs the default stack: si32 np4 M -3.80 %, G -1.57 %, MG -5.15 %; si64 np8 M -3.25 %, G -1.59 %, MG -3.18 %; SCF counts and H columns identical, max |dE| 5.1e-13 Ha. M2_REUSE (A off, B on): H columns -24.8..-28.3 % in all six cells, SCF counts identical, wall -1.2 % (si32) .. -5.1 % (si64 np4), max |dE| 1.14e-13 Ha. | **All three adopted as defaults (2026-11-17)** (never slower, rounding class, no work added); rollbacks `BAANDI_M2_STICKS=greedy`, `BAANDI_M2_START=ortho`, `BAANDI_M2_REUSE=off`. |

## Per-core kernels (PERCORE-KERNELS, 2026-11-17; Mac observations O, pure build, counters exact C, relayed dev08 R; analysis docs/performance/percore-kernels-2026-11-17.md, dev08 request docs/performance/requests/percore-kernels.md; stage 0 = `scripts/fast_screen_dav.py`, si32/si64-gamma f1e-3 seed 601 T=1, converged base density, 3 alternating pairs x 5 repeats; load 2 to 21 during the session, spreads as printed)

| ID | Idea | Status | Measured | Lesson / next |
|---|---|---|---|---|
| PERCORE-TABLE | Per-kernel table of one lean solve at one core (stage timers, gemm_dispatch, `sample` 6995 samples, microbenchmarks at the traced shapes) | finding (O, C) | si32 warm 387 ms / 176 H columns: local pair FFT 27.5 % (1.16 ms per pair transform), Ritz rotations 20.2 % at 24 GFLOP/s (narrow k <= 10: 4.5 to 10), reduced border 10.5 % at 21.5, overlap border 9.9 % at 22.7, nonlocal 12.4 % at 22, eigen 6.9 %; build's library rates 38 to 58; `memmove` from the pure GEMM packing 4.2 % of the solve. Per H column (structure only, different machines): H apply 0.89 ms (QE dev08 1.29), algebra 1.06 ms (QE 0.54), algebra GEMM work 22 MFLOP per column (QE about 10 to 12, E: basis size) | Ranking by share x gap: TN borders, narrow NN rotation, nonlocal products; the local FFT is at the library rate on dev08 (D addendum) and is not the per-core gap |
| PERCORE-LEAN-BORDER | `BAANDI_LEAN_BORDER=fused`: H applied to the unit corrections before the pivot test, one TN product V^T [E | H E] gives the overlap and the reduced border (basis streamed once, twice the width); failing pivot test discards the images and takes the separate fallback; final form: only blocks with 4 k <= m fuse, and fusing ends for the solve after the first failed pivot test | CANDIDATE (stage 0, Mac); dev08 Part A/B requested | (C) reduced-site GEMM calls 10 -> 1 per warm solve, GEMM FLOP, H columns and Davidson iterations unchanged, `fused_borders` 9 (si32 warm) / 10 (si64 warm), `fused_wasted_columns` 0; eigenvalues bit-identical (aarch64 pure: one FMA chain per entry, any width). Kernel microbenchmark (O): narrow k = 1..12 0.49 to 0.88 of two products, wide 0.93 to 1.17. Stage 0 (every append fused) si32: warm 0.964 (spread 1.75 %) and 0.949 (5.2 %), cold 1.004 (1.0 %) and 0.978 (4.0 %); si64: cold 0.921 (5.0 %), warm 1.001 (1.7 %); final form si32: 7 of 9 warm appends fused, cold unchanged, time void (load 11). Synthetic fixtures without the rules: near-degenerate 439 of 1296 H columns discarded (capacity 3), cubic 20 of 64; with both rules at most one narrow block (10 columns) per solve | Gain is in the narrow appends; a fused append bets that the pivot test passes (H is applied before it), so the rules bound the loss where it fails often; on MKL the TN borders are further from the ceiling (DGEMM-SITES), so dev08 decides; rounding class on MKL/x86 |
| PERCORE-PACK-FIXED | `BAANDI_PACK_FIXED=on`: aarch64 pure kernel packs full strips and stores full tile columns with constant-length copies (no `memmove` per eight doubles), transposed strips in gather order | CANDIDATE (stage 0, Mac; aarch64 pure only) | bit-identical (eigenvalue FNV, counters); first version (contiguous branch only): si32 cold 0.992 (1.9 %), warm 0.976 (1.9 %), si64 cold 1.001 (7.0 %), warm 0.999 (2.2 %); with the gather branch, load 2 to 3: si32 cold 0.963 (2.6 %), warm 0.947 (1.0 %) | No effect on dev08 (x86 pure runs faer, MKL runs dgemm); only for aarch64 hosts; a pure-build default needs the Mac generality set, not dev08 |
| PERCORE-RITZ-LOOP | `BAANDI_RITZ_FUSED=loop`: fused Ritz and residual kernel at the selected-root residuals only (return keeps the pair rotation for REUSE-PSI; new premise against RITZ-FUSED-LEAN-SCREEN) | REJECT on the pure build (stage 0); one MKL stage-0 row requested | bit-identical; si32 cold 0.999 (2.4 %), warm 1.019 (1.2 %); microbenchmark fused/pair 0.80 to 1.69 (mostly slower: the per-panel tile combine costs more than the saved n x k streams) | Do not repeat on pure builds; MKL only if Part A shows a gain |
| PERCORE-NARROW-REDUCED | `BAANDI_NARROW_REDUCED=on` on the lean path (premise changed: lean borders) | NEUTRAL (stage 0) | fires in 1 of 10 warm iterations (k <= 4 only); si32 cold 0.999 (1.0 %), warm 0.980 (4.5 %); rounding class | Superseded by the fused border, which covers every narrow border without a second kernel |
| PERCORE-STACK | `LEAN_BORDER=fused` + `PACK_FIXED=on` | mixed (stage 0, Mac): harness CANDIDATE in two runs (warm), REJECT in one (cold, si32 at load 3) | bit-identical; si32 (load 6 to 8) cold 0.989 (3.0 %), warm 0.947 (4.3 %); si32 (load 3) cold 1.063 (5.3 %, pairs 1.085 / 0.979 / 1.063; border alone cold 1.004 / 0.978 and pack alone 0.963, so the wide cold appends are the open question: kernel microbenchmark wide fused 0.93 to 1.17), warm 0.937 (4.7 %); si64 (load 3 to 7) cold 0.987 (1.1 %), warm 0.983 (15.9 %, one outlier pair). A five-pair border-only rerun at load 7.5 had a 59 % spread (void) | bit-identical on the Mac; dev08 stack row only for LEAN_BORDER (PACK_FIXED is a no-op there) |

## MPI overlap-matrix Davidson (MPI-OVERLAP-DAVIDSON, 2026-11-17; Mac observations O on simulated ranks, counters exact C, estimates E; design docs/design/mpi-overlap-davidson.md, dev08 request docs/performance/requests/mpi-overlap-davidson.md; review F05 step 2)

| ID | Idea | Status | Measured | Lesson / next |
|---|---|---|---|---|
| MPI-OVERLAP-DAVIDSON-LEAN | `BAANDI_M2_LEAN_ALGEBRA=on|off` (default off): standard (Ritz-change) M2 solves run the serial lean algebra on the stick rows: raw corrections appended without Gram-Schmidt, one tree sum for the S border (with the selected residual norms), unit scaling from the S diagonal, pivot test (eps/1e-12) before H, serial stage-1/stage-2 fallbacks, H border at the next iteration; generalized reduced problem A = T H_r T^T replicated; kept and Gram starts accepted by the pivot rule | CANDIDATE (stage 0 and 1; dev08 Parts A/B) | (C) tests: eigenvalues within 1e-10 Ha of the orthonormal path (cubic, degenerate, skewed; cap 3 and 6; ranks 1-4; cold and warm), rank decisions, iterations and H columns equal to the serial lean solver (incl. the ill-conditioned fixture: stage 1, stage 2 to the random block), ranks bit-equal, fixed-np reproducible; SCF test standard class 0 to 1.2e-14 Ha vs threads-only, verified. (O, stage 0, warm solve B/A) si32 np 2/4/8 0.747 / 0.876 / 0.815, si64 np 4/8 0.862 / 0.853; cold 0.86 to 0.94; tree sums warm si32 -79 / -62 / -61 %, si64 -73 / -72 %; H columns and iterations equal; stage clocks si32 np=8 warm: orthonormalize 67 -> 7 ms, algebra stages 120 -> 84 ms. (O, stage 1 production) si32 np=4 -7.5 % (spread 0.9 %), np=8 -10.5 % (within the 13 % spread), si64 np=4 -6.1 % (spread 4 %), tree sums -66 %, counts unchanged; with the pending stack S -> S+L si32 np=8 -10.3 % (spread 7 %), tree sums 4494 -> 511 per rank | Two reductions per iteration instead of about 4 + k; the remaining chains are the start orthonormalization (use with M2_START=gram) and stage 1 at capacity 3 (backlog) |
| MPI-OVERLAP-DAVIDSON-ONESUM | Same algebra with ONE tree sum per iteration: S border, H border and residual norms reduced together after H is applied to the raw corrections (pivot test after the apply) | REJECT (stage 0) | (O, si32 np=4 warm solve) 15 of 22 appends fail the pivot test, 42 H columns per solve applied to rejected blocks: H columns 272 -> 314 (+15 %), alltoall calls 94 -> 154 (+64 %); tree sums 683 -> 254 (vs 261 for the two-sum form), solve 0.952 (vs 0.876) | Do not trade H applications (the stage that scales) for one latency-bound reduction; the pivot test must precede the apply. The serial counterpart PERCORE-LEAN-BORDER (main 650eed2) fuses only while no pivot test has failed in the solve; on the M2 path the second reduction is the S border before H, which costs no H column, so that guard is not needed here |
| MPI-OVERLAP-DAVIDSON-PIVOT-RATE | Observation: late warm solves at capacity 3 (np >= 4) fail the pivot test on most appends (si32 15/22, si64 14/23; capacity 6 at np=2: 0/9; production SCF screens: 0 in every run) | finding | (O, C) stage 1 costs one norm reduction, the in-block Gram-Schmidt of the block and a second S border; orthonormalize stage 7 ms of the 194 ms warm si32 np=8 solve | Check the serial lean rate at capacity 3 on the threads path before designing a cheaper stage 1 (backlog) |

| MIX-INCREMENTAL-STOP-METRIC-DEV08-ADOPT | 2026-11-17 | dev08 (master), 2963c86 mkl,mpi, SCF level real MPI, seeds 601-603, 30 rows, all converged, SCF counts identical | `BAANDI_MIX_INCREMENTAL=on` (I) and `BAANDI_STOP_METRIC=kept` (K) vs default (A): density SHA identical and |dE| = 0 in every row (bit-identical in these cells). Alone: si32 np4 I +0.22 %, K +0.45 % (inside the spread), si64 np8 I -2.11 %, K -1.85 %. On top of the sticks+start+reuse stack S: si32 np4 S -5.58 % -> SIK -8.26 %, si64 np8 S -6.86 % -> SIK -9.70 % (2.7-2.8 points added). | **Both adopted as defaults (2026-11-17)** (bit-identical in the screened cells, never slower beyond noise, additive on top of the MPI stack); rollbacks `BAANDI_MIX_INCREMENTAL=off`, `BAANDI_STOP_METRIC=fft`. |

| LEAN-BORDER-DEV08-NEUTRAL | 2026-11-17 | dev08 (master), 6c4354c mkl, si32 seed 601 T=1; stage 0 (one solve) and SCF level (to convergence, then 10 fixed steps) | `BAANDI_LEAN_BORDER=fused`: stage 0 warm solve ratio 0.9758 (fused borders 0 -> 7, wasted H columns 0, H columns equal). SCF level to convergence: B vs A +2.31 % inside the 6.9 % spread, SCF 10/10, H columns 1717/1717. 10 fixed steps: only the GEMM call count moved (1250 -> 1218, -2.56 %), GEMM GFLOP identical (64.873), H columns equal, |dE| 1.4e-14 Ha, wall inside the spread. | **Closed as neutral; stays default off** (the fusion removes calls, not flops, so it does not show in the SCF wall). |
| M2-LEAN-ALGEBRA-DEV08-STAGE0 | 2026-11-17 | dev08 (master), 6c4354c mkl,mpi, davidson_bench real MPI np=8, 108 JSON, 3 pairs | `BAANDI_M2_LEAN_ALGEBRA=on` (L) vs default (A): si32 cold -6.96 %, warm -9.41 %; si64 cold -7.39 %, warm -12.77 %; tree_sum calls 155 -> 15, 616 -> 187, 289 -> 15, 1255 -> 265; H columns and Davidson iterations identical (190/204/374/478, 3/23/3/25); band energy difference 5e-15..2e-14 Ha (rounding class, eigen FNV differs as declared). Arm O (`M2_LAZY=off`) is +11..24 % slower than A and OL (lazy off + L) is not better than A: L adds on top of LAZY, not instead of it. | Candidate; 10-step fixed SCF confirmation (si32 np=8, <1 min) queued, then default on. |

## Breakthrough restart (BREAKTHROUGH-1-RESTART, 2026-11-17; review 2026-11-17 R1 / E1 / E2, backlog B05; Mac observations O on simulated ranks and T=1, counters exact C; design docs/design/restart-keep.md, dev08 request docs/performance/requests/breakthrough-restart.md)

| ID | Idea | Status | Measured | Lesson / next |
|---|---|---|---|---|
| R1A-SAME-H-CONTINUE | Same-H continuation of the first-step standard redo (continue the first solve's basis, H V, reduced matrices and last Ritz change at the tighter tolerance instead of a warm restart) | CLOSED before code: ceiling zero | (C) redo flag `standard_first_iteration_redo` set in 0 of 9 default runs (si8/si16/si32 seeds 601-603, ladder on) and 0 of 6 with `BAANDI_CUTOFF_LADDER=off` (si8/si32); kept-start counter `unchanged_potential` 0 in all 15 | The first map is never within tolerance x N_e of its fixed point (ladder seed or 5e-3 cap from the atomic start). Reopen only for a workload with redos (restart from a converged density, relaxation with density transfer); the contract is written in the design, section 1 |
| RESTART-KEEP-R1B | `BAANDI_RESTART_KEEP=off|1|2|4|8|b/8|b/4|b/2` (default off): at a lean contraction keep up to r previous Ritz vectors of the unconverged roots minus their X component, P = V (D - C A), H P = (H V)(D - C A) from the same rotation as X (no H application), rank floor DEPENDENCE_TOLERANCE, pivot rule eps/1e-12, fallback to X alone; serial lean and the distributed lean algebra (`M2_LEAN_ALGEBRA=on`), zero extra reductions per restart | REJECT (stage 0 and stage 1) | Stage 0 si32 np=8 warm solve from the start density (4 contractions in A; the converged-density warm and both cold solves have none): H 617 -> 619 / 618 / 657 / 623 (r = 1 / B/8 / B/4 / B/2), rotation -2 / +6.5 / +26 / +45 %, GEMM +3.5 to +41 %, contractions 4 -> 6-7, occupied eigenvalue moved 1.3-2.0e-6 Ha (lower Ritz values: the plain restart stalls at this rough potential); 3-pair screen r = 1 REJECT (H 824 -> 826, warm 1.087, spread 6.5 %). Stage 1 si32 np=8 production: r = 1 H 1888 -> 1877 (-0.6 %), Davidson 67 -> 70, rotation -1.2 %, sum m^3 +2.7 %; r = B/4 H +0.6 %, rotation +8.4 %, sum m^3 +5.9 %, contractions 11 -> 13; r = B/2 H -1.6 %, rotation +16 %, sum m^3 +15 %; |dE| <= 3.2e-11 Ha. Threads path T=1 (capacity 6): 0 contractions in si8/si16/si32 SCF, bit-identical | At capacity 3 restarts are rare (11 in 22 solves of a si32 np=8 SCF) and the kept columns take the room of the next correction block, so the basis fills sooner and restarts more often; no H column is saved. Synthetic skew fixture (many restarts): H +16 / +8 / -1 (r = 1 / B/4 / B/2) at 1e-11; the Ritz-change stop fires elsewhere (one setting 9e-8 against 1.1e-8 Ha plain error at tolerance 1e-9). Do not retry at capacity 3 or 6 |
| RESTART-KEEP-CAP2-LEAD | Same switch at capacity 2 per band (`M2_PROFILE_CAPACITY=2`), new premise against MPI-SCALING-GAP-CAP2 | lead (stage 1 Mac counters only; dev08 job 1 of the request) | (C, si32 np=8 production, seed 601) capacity 2: SCF 13, Davidson 116, H 2720, 51 contractions, sum m^3 2.70e8; + r = B/4: SCF 11, Davidson 91, H 2206 (-19 %), sum m^3 2.28e8, rotation -1.5 %, tree sums -26 %, |dE| 1.3e-10 Ha; + r = B/2: H 2189. Against capacity 3 (H 1888, sum m^3 4.21e8, rotation 3.12e9): +17 % H, +25 % rotation, -46 % sum m^3 | Kept directions remove 62 % of the capacity-2 H excess (two SCF steps fewer); still +17 % H against the current rule, so a dev08 loss is expected (CAP2 lost at +23 %). One seed, one cell; a capacity change needs stage 2/3 under the generality rule |
## Bounded subspace SCF (BREAKTHROUGH-3-BOUNDED-SUBSPACE, 2026-11-17; Mac observations O, counters exact C; design docs/design/bounded-subspace-scf.md, dev08 request docs/performance/requests/breakthrough-advance.md; review 2026-11-17 section 7, R3). Pure build, T=1, full SCF to convergence incl. verification, PseudoDojo Si.upf (PBE, never committed) with `xc.functional = "pbe"`, benchmark inputs from `scripts/fast_screen.py --xc pbe` (f1e-3), one run per arm (counters are deterministic; walls at load 4 to 6 are not judged). Switch off is bit-identical to main d5d04b1 (si16 s601 density SHA 085b9c28..., energy -67.42404477459324 Ha, baseline binary 73ac5dba...).

| ID | Idea | Status | Measured | Lesson / next |
|---|---|---|---|---|
| BT3-ADVANCE-CONTRACT | `SolveIntent::Advance(SubspaceBudget)` on the serial lean Davidson (`RealDavidson::advance_reusing`): at most N correction blocks / H columns from the kept start, returns the all-band Ritz state typed `AdvancedSubspace` with true Ritz-change flags; session (`BAANDI_SCF_ADVANCE=on`, default off) builds density, energy and the Pulay pair from it, never stops on it, completes it by the exact solve at the same density (warm = advanced X, kept start taken, potential unchanged) on candidate (E < threshold), inner noise (B > 0.1 E, the ratchet's factor), growth (E > E_prev) or gap crossing; first map, redo and every completion exact; verification unchanged (ADR 0043 V1-V6 intact) | implemented, default off; REJECT as a speed idea (stage 1) | (C) si16 s601 total H columns 981 -> 1359 (+38 %; SCF 9 -> 11, verification attempts 2 -> 3, 16 advanced maps, 7 completions: 4 growth, 2 inner noise, 1 candidate); s602 824 -> 1125 (+37 %); s603 1183 -> 1172 (-1 %); si-skew s601 127 -> 155 (+22 %). |dE| vs default 2.5e-10 to 7.5e-10 Ha (standard class kept), all converged and verified | Per advanced map the saving is real but small (42 to 71 H columns against 53 to 107 for the exact warm maps at the same SCF stage: the ratchet already stops exact warm solves after 3 to 7 Davidson iterations and the first all-band block dominates both). The saving is eaten by completions (36 to 84 columns each), one more SCF step and one more verification attempt in two of four geometries. The mixer does not visibly learn the inner noise (E per iteration tracks the exact path while B <= 0.1 E); the endgame is what moves. Not A21 (no tolerance scaled) and not A22 (no converged RMM). |
| BT3-GROWTH-TRIGGER | Growth completion (`BAANDI_SCF_ADVANCE_GROWTH=off` removes it) | REJECT the trigger | (C) si16 s601: the exact path's own E grows at iteration 3 (1.8e-4 -> 3.5e-3 Ha, Pulay); with the trigger 4 of 7 completions were growth completions on maps whose inner ratio was 0.008 to 0.03. Without it: 1359 -> 1101 H (still +12 % vs default 981); s602 1125 -> 818 (-1 %); s603 1172 -> 1170 (-1 %); si-skew 155 -> 140 (+10 %) | Residual growth is a mixer property, not an inexactness signal; the inner ratio already separates the two. Keep the switch for research, default on stays (it was the pre-registered rule) only because the whole idea is off. |
| BT3-BUDGET-VARIANTS | Blocks 1 / 2 / 3 fixed (`SCF_ADVANCE_ADAPT=off`), adaptation on/off, fine stage only (`SCF_ADVANCE_SCOPE=fine`) | REJECT (no consistent winner) | (C) total H columns, default -> variant: si16 s601 981 -> b1 1352 / b2 1237 / b3 1016 / b2-fine 1046; s602 824 -> b2 818 / b3 1208 / b2-fine 1223; s603 1183 -> b2 1170 / b3 1107 / b2-fine 952; si-skew 127 -> b2 140 / b3 145. 1 block: 10 of 11 advanced maps complete by inner noise (a single all-band block after a large mixer step leaves B >> 0.1 E) | The only win (s603 fine-only -19.5 %) is one fewer verification attempt (certificate H 197 -> 97) and one fewer SCF step; the same variant loses 48 % on s602. Trajectory luck, not a mechanism. The gate (>= 5 % fewer H columns on one case, then three seeds) is not met on any variant. |
| BT3-COARSE | Advance inside the cutoff-ladder coarse stage (scope all vs fine) | finding | (C) coarse H columns 659 -> 586 to 638 (-3 to -11 %) on si16 s601 with 7 to 9 coarse steps; the fine hand-off then differs and the fine stage pays it back (+1 to +3 fine H columns per saved coarse column in the runs above) | A cheaper coarse stage that delivers a worse hand-off is not a saving; R5 of the review (ladder re-evaluation) should measure hand-off quality, not coarse cost alone. |
| BT3-INNER-RATIO | Observation: inner/SCF ratio B / E of accepted 2-block advanced maps 0.008 to 0.09; after a large Pulay step (E dropped 3 orders) a 2-block advance leaves ratio 1.5 to 5.9 and completes; near the stop (E ~ 1e-9) ratio 0.3 to 1.3 | finding | (C) si16 s601/603 histories in the request document | The budget a map needs is set by how far the mixer moved the density, which is known only after the solve: a budget rule on observables before the solve (E_prev, previous ratio) cannot see the step size of the current mixing. A posteriori completion is then the only safe rule, and it costs a second potential/density/audit. |

## Measurement infrastructure (BREAKTHROUGH-0-INSTRUMENT, 2026-11-17; review 2026-11-17 sections 13.1 to 13.3 and E3; Mac observations O on simulated ranks np=4 and T=1, counts exact C; notes docs/performance/profile-2026-11-17-phases.md, docs/performance/first-block-residual-2026-11-17.md; dev08 request docs/performance/requests/breakthrough-instrument.md). No speed claim; three default-off switches whose on/off state leaves energy bits and density bits unchanged (main d5d04b1 harness against this branch, si8/si32 np=4 production + ladder; dft-cli si32 T=1 density SHA 9491950d...; test `dft-engine/tests/instrument.rs`).

| ID | Idea | Status | Measured | Lesson / next |
|---|---|---|---|---|
| BT0-PHASE-PROFILE | `BAANDI_PHASE_PROFILE=on`: exclusive nested phase profile per rank thread (init / coarse / fine / fine_redo / verify_pass / verify_fail / materialize; map phases; exclusive solve stages incl. the new `restart` stage; collectives and operator counts beside them; max and min over ranks, never a sum); results.json `phase_profile` / `parallel.phases`, `m2_entry_cost`, `davidson_bench`; count-only with `STAGE_TIMERS=off` | adopted as a tool (default off) | (O, si32 np=4 production, wall 3.10 s, covered 97 %) init 5.8 % (atomic start 3.3 %), coarse 18.3 %, fine 61.3 % (solve 50.5 % of the wall: apply 43 % of it, orthonormalize 15 %, rotate 12 %, reduced border 9 %, EVD 6 %, residual 3 %, uncovered loop 10 %), verification 11.4 % (one failed + one passed attempt, each 5.7 %), remainder 3.0 %; fine H 2568 = 1888 full + 680 delta-V columns, FFT pairs 1295 fwd / 1327 inv; rank imbalance of exclusive stages 2 to 5 % except residual / other (waits) | the earlier "outside the recorded fine phases 1.28 s" is now coarse 0.57 + verification 0.35 + init 0.18 + mixing 0.05 + 0.09 remainder; use these buckets as the denominators of every R1 to R7 claim (MPI against MPI) |
| BT0-SNAPSHOT-REPLAY | `BAANDI_SNAPSHOT_DIR=<dir>`: every real Gamma solve writes density (rank 0), this rank's start rows with the armed kept H X, manifest with context and outcome, outside the stage clocks (IO subtracted from the map's diagonalization); `davidson_bench --snapshot <manifest.json>` replays it (`ScfProblem::replay_solve`, injection keyed by rank) | adopted as a tool (default off) | (C) replay reproduces Davidson iterations, H columns and the eigenvalue FNV of the captured solve: threads si32 T=1 cold / kept it06 / verification (4-203, 4-112, 5-201 H), 4 simulated ranks si32 warm (30 iterations, 617 H), 2 simulated ranks in the engine test (cold, kept, verification); the captured X of a kept start must be the solver's kept block, not the round-tripped start (one ulp changes the eigenvalue bits) | stage-0 A/B of real mid-SCF states (coarse middle, hand-off, fine early/middle/late, verification before and after a failure) is now a replay, not the synthetic warm map; the delta-V update itself is not replayed (its H X is stored); the same-H redo never occurred in these runs (nothing to capture; see R1A-SAME-H-CONTINUE) |
| ISSUE18-SNAPSHOT-IDENTITY | Replay identity gate (issue #18): versioned operator fingerprint in the snapshot manifest, per-rank row-hash, potential, stopping-profile and kept-H X checks before injection, agreed over ranks; `davidson_bench --replay-mode exact|ab|legacy` (exact mismatch exit 3, unreproduced exact replay exit 4; A/B recomputes H X on a changed potential rounding) | adopted as a tool (measurement only, default SCF unchanged) | (C) `tests/snapshot_identity.rs`: original kept snapshot reproduces iterations, H columns and eigenvalue FNV (1 rank and 2 simulated ranks); rejects moved atom, same-valence PP, XC, swapped and rank-local row order (agreed on both ranks), one flipped density bit, strict vs standard profile and strict-band change; A/B with tampered potential FNV drops and recomputes H X | Replays recorded before this change (BT0-SNAPSHOT-REPLAY evidence) were unchecked; re-run them in exact mode before using them as decision evidence. |
| BT0-FIRST-BLOCK | `BAANDI_FIRST_BLOCK=on` (review E3, read-only): per warm solve a small Rayleigh-Ritz on the kept B vectors, true residuals of the B Ritz pairs, singular values of R_o = H X_o - X_o A_o, effective ranks at 1e-1..1e-4, counts of unconverged roots by the Ritz-change rule and by the true residual | finding (no solver change) | (O/C, si32 np=4 and T=1, si16, si8) late fine maps (from fine it 6 to 7 on): rank 1 at 10 % of sigma_max, 32 to 60 of 65 at 1 %, 65 at 0.1 %; sigma_min 3e-7 to 1e-5 Ha against tolerances 1e-12 to 1e-9: every occupied direction stays above the threshold (`r>tol` = n_o always); `unc_all` = 68 = the rule's first block at every warm solve; ladder hand-off is the worst start (sigma_max 1.3e-1, full rank) against 2e-3 for the next kept start; verification maps start as far from convergence as ordinary maps (sigma_max 2 to 3e-4 against tolerance 2e-12); boundary gap stable 3.4 to 4.5e-4 Ha | the Ritz-change rule does not over-correct the first block on this trajectory (the review's R2 premise is not confirmed for the FIRST block); a relative-rank compression would drop directions still 4 to 6 orders above tolerance, so R2 can only be judged on a replayed late-step solve by H columns to the same converged result; the hand-off quality (R5) is the larger first-block lever |

| LADDER-REEVAL-DEV08 | 2026-11-17 | dev08 (master), current default stack (081ff8a+) mkl,mpi, si32 np=8 production, seed 601, to convergence, 3 alternating pairs | cutoff ladder on (default) 3.480 s (spread 0.035) vs off 4.083 s (spread 0.032): off is +17.31 % slower; the ladder earns 14.8 % on the current stack; cost condition T_coarse + T_transfer < T_fine(no ladder) - T_fine(ladder) holds (difference 0.603 s). | **Ladder stays on by default** (review R5 / E6 closed). |
| M2-DENSITY-FUSED-DEV08-HOLD | 2026-11-17 | dev08, si32 np=8, 10 fixed steps, 3 pairs | `BAANDI_M2_DENSITY=fused` -3.04 % (spread 0.71 %) but no usable counter moved. | Held: judged again with the exclusive phase profile (density phase seconds) before any adoption. |
| M2-LEAN-ALGEBRA-DEV08-10STEPS | 2026-11-17 | dev08, si32 np=8, 10 fixed steps, 3 pairs | `BAANDI_M2_LEAN_ALGEBRA=on`: per-step 0.3490 -> 0.3423 s (-1.91 %, larger spread 1.53 %); collective calls 4809 -> 853 (-82.3 %), reductions 4273 -> 317 (-92.6 %); H columns 2254 = 2254, Davidson 49 = 49, SCF 10 = 10. The stage-0 solve gain (-9.4 %) mostly disappears at SCF level (latency terms are a small share of the step). | Weak candidate: si64 np=8 10-step and a repeat queued; adopted if si64 >= -3 % or the repeat agrees. |
| M2-LEAN-ALGEBRA-DEV08-ADOPT | 2026-11-17 | dev08 (master), a220a85 mkl,mpi, np=8, 10 fixed steps, 3 alternating pairs each | `BAANDI_M2_LEAN_ALGEBRA=on` vs default: si64 per-step 1.2225 -> 1.1756 s (-3.83 %, spread 0.56 %), collective calls -86.6 %, reductions -95.3 %, GEMM calls -13.9 %, H columns/Davidson/SCF identical, |dE| 2.8e-14 Ha; si32 repeat -2.59 % (spread 0.53 %; first look -1.91 % with a 1.53 % spread), collective calls -82.3 %, reductions -92.6 %, work identical, |dE| 2.8e-14 Ha. Stage 0: solve -7..-13 %. | **Adopted as default on (2026-11-17)** (pre-registered criteria met: si64 >= -3 %, si32 repeat same direction; work and accuracy unchanged, rounding class); rollback `BAANDI_M2_LEAN_ALGEBRA=off`. |

## Other and init (BREAKTHROUGH-4-OTHER-AND-INIT, 2026-11-17; Mac observations O on simulated ranks np=8/4/2 and T=1, counts exact C; note docs/performance/profile-2026-11-17-other-init.md, dev08 request docs/performance/requests/breakthrough-other-init.md). Starting point (relayed dev08, a220a85, si32 np=8): fine `solve.other` 0.2601 s = 24.7 % of the fine solve, remainder 0.2552 s, init 0.5606 s.

| ID | Idea | Status | Measured | Lesson / next |
|---|---|---|---|---|
| BT4-LABELS | Finer exclusive labels of the phase profile: solver `setup`, `precondition`, `append`, `solver_total` (other split into `other_solver` / `other_engine`), engine `real_state`, `prepare` and `unpack` now subtracted, `entry`, `kept_update`, `finish`; init construction parts and atomic-start parts; new bucket `between` (ladder hand-off and truncation, session set-up, stop rule, verification glue); the coarse construction filed under `init`; per-rank derivation of `other` before max/min | adopted as a tool (profile off by default, results bit-identical: energy bits 0xc05fad02a3080b2a, density FNV 0xa6d327ec57f20279 with and without) | (O, si32 np=8) fine `other` 0.188 -> 0.0005 s, coarse 0.095 -> 0.0003 s, remainder 3.0-5.6 % -> 0.3 %; CLI threads si32 T=1 remainder 0.2 % | `solve.other` was NOT the solver's loop body (other_solver ~0): it was the engine's solve interface, mainly the serial real Hamiltonian rebuild (`real_state`) and the kept-start delta-V update (`kept_update`) |
| BT4-PHASE-WAIT | `BAANDI_PHASE_WAIT=on` (diagnostic, default off): timed barrier before every collective of the distributed solver (`solve.wait`, not exclusive) and at every M2 solve entry (`solve.entry_wait`) | tool (timing only, never an acceptance arm) | (O, si32 np=8, loaded Mac) fine wait 0.249 s in 2547 barriers, coarse 0.144 s, entry_wait fine 0.022 s; bits identical | Real-MPI split of wait from transfer is Job 1 of the request; simulated-rank waits are CPU contention, not imbalance evidence |
| BT4-LIGHT-STATE | `BAANDI_M2_LIGHT_STATE=on` (default off): the replicated serial real Gamma Hamiltonian of an M2 solve state (rebuilt at the first map of every ladder stage) is built with block capacity 2 (one pair plan) instead of the Davidson capacity (capacity/2 = 102 full-grid pair plans + real projector scratch per rank); it only supplies the ORC half sphere there; full build if the distributed operator is unavailable | CANDIDATE (stage 1 O; dev08 Job 2/3) | (O) `real_state` fine 0.1095 -> 0.0107 s, coarse 0.0548 -> 0.0026 s; si32 np=8 wall B/A 0.918 (5 pairs: 0.866 1.073 0.915 0.961 0.918); bit-identical (energy bits, density FNV, eigenvalues; test `tests/other_init.rs` cubic and skewed, 2 ranks); peak RSS not reduced (2.2-2.6 GB both arms: the untouched scratch was not resident) | The hidden term of `solve.other` on dev08 (0.26 s) is most likely this rebuild (plan creation and allocation per rank, 8 processes on one node); confirm with Job 1 |
| BT4-EWALD-CACHE | `BAANDI_EWALD_CACHE=on` (default off): a construction on the same thread with bit-equal lattice, charges and fractional positions reuses the Ewald energy (ladder coarse problem after the full-cutoff one) | CANDIDATE (counter; wall inside the noise on the Mac) | (C) Ewald evaluations 2 -> 1 per ladder run; (O) `init.ewald` 0.067 -> 0.043 s (si32 np=8, max over ranks), 0.041 s for two on the threads path; wall B/A 0.985 (5 pairs, spread larger) ; bit-identical | One Ewald is about 1 % of the si32 np=8 wall; acts on every path with the ladder (threads too) |
| BT4-LAYOUT-CACHE | `BAANDI_M2_LAYOUT_CACHE=on` (default off): reuse the solver's row-index layout of the previous solve when all ranks' local rows are unchanged (one 1-double tree sum decides, replicated) instead of two allgathervs (one of N doubles) and the index rebuild per solve | CANDIDATE (counter; latency item, dev08 decides) | (C, si32 np=8) allgatherv calls 117 -> 77 per rank, payload -1.4 MB per rank, tree sums +22; (O) wall B/A 1.000 (5 pairs; simulated ranks do not model latency); bit-identical | Backlog item "cacheable per operator" (MPI scaling gap follow-ups) done |
| BT4-STACK | LIGHT_STATE + EWALD_CACHE + LAYOUT_CACHE together | CANDIDATE (stage 1 O; dev08 Job 2) | (O) si32 np=8 B/A 0.936 (7 pairs; A 3.097 s, B 2.934 s), si16 np=4 0.959 (3 pairs), si8 np=2 0.930 (3 pairs); identical energy bits and density FNV in every run of all three cases | Stacked candidate for dev08 stage 1 (one case si32 np=8, 3 pairs); widen only after |
| BT4-KEPT-UPDATE | Observation: `kept_update` (the MPI-KEPT-START delta-V columns) is the second hidden term of `solve.other` (fine 0.072 s, coarse 0.026 s, verification 0.014 s, O) | finding (no change) | (C) 680 fine / 544 coarse delta-V columns; 2 of 34 pairs not kept (partly occupied + empty) force the mirror exchange and an inverse transform | Inherent transform work; only the partly occupied pair's slab could still be kept (backlog) |
| BT4-ATOMIC-SPLIT-CHAIN | Observation: with `M2_ATOMIC=split` the distributed two-pass Gram-Schmidt of the 128 span columns is a chain of 128 broadcasts | finding (no change) | (O) np=8 orthonormalize 0.058 s against 0.069 s serial (threads T=1); Ritz 0.047 s = gram 0.030 + eigh 0.004 + rotate 0.013 | The split saves little; a block orthonormalization changes the span's rounding (trajectory, rounding class, stage 1 to convergence); the rotation can be split by band columns (backlog) |

### LEAN-ALGEBRA-HARDEN (2026-11-17; correctness, no speed claim)

| ID | Status | Evidence | Mechanism and lesson |
|---|---|---|---|
| LEAN-ALGEBRA-HARDEN | Implemented; review pending | Original `0fb09c0`: 8/8 adversarial integration tests fail (Gram defect, 3 s error-agreement timeouts, partial-start panic, empty-rank reduction mismatch); shared inverse-factor unit test also fails. Patched: 9/9 integration cases pass; forced stage 2 after contraction and nonempty rank-loss rebuild unit tests pass. 24 distributed cold/warm fixture records and 6 serial records retain all eigenvalue/vector/residual bit hashes and decision counts; normal tree-sum calls unchanged. `cargo test -p dft-solvers --offline`: exit 0 (133 passed, 25 ignored); `scripts/check.sh dft-solvers`: exit 0, all four gates PASS. Apple Silicon pure Rust only; no dev08 timing or SCF-energy claim. | A good sequence of Cholesky pivot ratios does not bound the inverse factor. Check inverse growth and actual return/contraction Gram defects; carry operator failure flags in existing sums; orthonormalize partial starts before completion; retain replicated panel width on empty ranks. Details and baseline hashes: `docs/design/mpi-overlap-davidson.md`, section 9. |
## ONEPATH-5-SPIN-SHARE (2026-11-17)

| ID | Mechanism | Status | Evidence | Limits / lesson |
| --- | --- | --- | --- | --- |
| ONEPATH-5-SPIN-SHARE | `BAANDI_SPIN_SHARE=on` (default off): one complex Davidson solve per k point for exactly equal fixed-sector standard spin channels; copy spectral data and keep both density/audit paths. Compare full effective potentials bitwise every map (zero-bit threshold), including verification. Decline disabled, strict, MPI, FD, nonzero seed, unequal sector, unequal/fractional occupations or band counts, different potentials, or different warm starts. Per-run count-only report; copied channel retains logical iterations but records zero H columns. | Implemented, default-off candidate. Ownership extension approved by the orchestrator: both collinear CLI writers expose `spin.spin_share` (shared/fallback maps, reason counts, actual Davidson H columns, logical Davidson iterations). Count-only; physical state schema and digests unchanged. No default adoption. | Baseline measured before engine edits on `b1c95a61282e66e741e98d1ef0199545faf1cae2`: synthetic skewed 5-Bohr cell, Gaussian nonlocal fixture, LDA-PW, 2 electrons, 4 bands, cutoff 2.5 Ha, standard+verification: E = -0.91551991195296512 Ha, bits 13829793526668695510. Off retains those bits. On T=1 and T=4: dE = 0 Ha, density/orbitals/eigenvalues/residuals bit-identical, M = 0 mu_B, 6 SCF steps, 46 logical Davidson iterations, 7 shared maps including verification, actual Davidson H columns 144 -> 72. Atomic-start dimer: LDA-PW H columns 262 -> 131, PBE 266 -> 133, atomic-start H columns 16 -> 8 in both; off/on spectral states, densities and energies bit-identical. Finite-T 0.01 Ha, unequal sector and nonzero-seed runs converge with only the expected fallback reason. Unit tests also inject zero-integral local magnetization, unequal/fractional occupations, band counts, unequal starts and a one-bit potential difference. | Mac pure-Rust correctness/work-count observation, no wall-time or QE speed claim. New tests: `cargo test -p dft-engine --test spin_share --offline -- --nocapture`; `cargo test -p dft-engine spin_share --offline --lib`. Five integration tests, two guard unit tests, three switch tests and the CLI manifest test pass. CLI JSON test covers on/off, nonzero seed and finite-T fallback (all verified), exact halving of H columns, and exclusion from physical state schemas; independent spin schema/cube validation and engine+CLI clippy also pass. Final engine check: build/fmt/clippy PASS, test FAIL for the baseline MPI issue below. Stage 1/2 dev08 through master still required: paired off/on fixed-sector zero-seed Si and non-Si cells, skewed/elongated and k meshes, pure/MKL T=1/4/8, verification included, H/atomic-start counts and wall separately. MPI and finite-T sharing excluded. Full engine check has an unrelated existing `distributed_scf::m2_panel_gs_variants_keep_classes` failure (tree sums [160,160,160]); reproduced with HEAD's original collinear source, exit 101. Orchestrator reports this is fixed on main at `0fb09c0` and authorizes ignoring it on this older base; no edit to distributed_scf.rs. The remaining complete engine suite passes with `cargo test -p dft-engine --offline -- --skip m2_panel_gs_variants_keep_classes` (exit 0; 189 passed, 70 ignored, 1 filtered). README generator produces no diff; its existing structural wording test expects obsolete `standard profile is unsupported`, exit 1 (known stale-text check, waived by the orchestrator). No new stored README evidence outside ticket ownership. |

| OTHER-INIT-DENSITY-DEV08-ADOPT | 2026-11-17 | dev08 (master), 0fb09c0 mkl,mpi, np=8, 10 fixed steps, session noise floor measured first (first pair discarded) | si32: `BAANDI_M2_LIGHT_STATE` alone -4.47 % (floor 0.21 %, 21x), `BAANDI_EWALD_CACHE` alone -1.98 % (floor 0.19 %), `BAANDI_M2_LAYOUT_CACHE` alone +0.23 % (floor 0.18 %: no time gain, 15 collectives became 19 reductions), all three -6.40 % (floor 0.26 %, 24.6x, energy and density SHA identical); si64: LIGHT_STATE alone -5.19 %, all three -7.83 / -7.96 %. LIGHT_STATE and EWALD_CACHE leave every usable counter identical (pure overhead removal). `BAANDI_M2_DENSITY=fused`: -3.04 / -3.42 / -3.21 % in three measurements (floor 0.39 %, 8x), exclusive phase profile shows the density phase 0.2727 -> 0.1690 s (-38 %), wall -0.1015 s matching the density saving -0.1037 s; init and solve unchanged. | **Adopted as defaults: M2_LIGHT_STATE=on, EWALD_CACHE=on, M2_DENSITY=fused (2026-11-17)**; `M2_LAYOUT_CACHE` stays off (no time gain). Rollbacks `BAANDI_M2_LIGHT_STATE=off`, `BAANDI_EWALD_CACHE=off`, `BAANDI_M2_DENSITY=bands`. |

| LEAN-ALGEBRA-HARDEN-OVERHEAD | Implemented; dev08 remeasurement pending | Relayed dev08 hardening cost: GEMM FLOP +2.15% si32 / +2.36% si64 np8, calls +1.1–1.4%, payload +5%, si64 wall +0.65%. Local stage-0 si32 fixed converged density np4/8 cold/warm: follow-up matches **every** `0fb09c0` counter and eigenvalue/residual/band-energy bits. Per-rank warm GEMM calls 397 -> 396; GFLOP np4 1.468470286 -> 1.441010662, np8 0.759092707 -> 0.745362895; reduction doubles 124524 -> 119740. Synthetic 24 work records and 30 bit/decision records also match original. | Replicated compensated C^T S C replaces the unconditional N-sized Gram; cancellation/threshold-sensitive cases retain the direct check at 1e-12. Failure-only NaN propagation reuses existing reduction payload. Naive m-scaled screening was overconservative and retained the ordinary warm Gram; compensated summation removed that source of uncertainty. Uniform-density warm stress cases still pay for the necessary direct check; np4 retains the pre-existing hardening rebuild. No local timing adoption claim. Evidence and all before/after counts: `docs/design/mpi-overlap-davidson.md`, section 10. |

| ISSUE19-AGREEMENT-SCRATCH | Workspace-owned agreement, previous-Ritz and cold diagonal scratch; dev08 timing pending | Against hardened HEAD `9c6d52d9eb40d1bdddec8977fed5b82eb2f5fc65`, 114 sorted fixture records match: 24 work records, 30 eigenvalue/vector/residual bit and decision records, 60 per-rank collective-order/payload and unit-occupation band-energy (Ha) records. SHA-256 of the combined records: `7219e2729bd8cf80bc7faf9b59dc249cf432614cba3ed962b511f4c1d6bd8271`. Commands: `cargo test -p dft-solvers --offline --test distributed_lean_algebra -- --nocapture`, before with original HEAD sources, after with patched sources; logs `/tmp/issue19-{before,after-final}-trace.log`. `agreement_scratch` checks four simulated ranks, empty partitions, growth/shrink, error/NaN then healthy recovery on lean and root-eigen workspace paths: agreement cold/growth requests 200/1448 bytes, fixed-shape reuse 0 allocations / 0 bytes for agreement and previous-Ritz scratch. Existing nine `lean_harden` tests retain factor/returned-vector guards. Final `cargo test -p dft-solvers --offline`: exit 0, 135 passed / 25 ignored; `CARGO_NET_OFFLINE=true scripts/check.sh dft-solvers`: exit 0, build/test/fmt/clippy PASS. | Current HEAD already forwards healthy reduction buffers directly using failure-only NaN propagation: no status slot, copy, payload or collective change is needed. The agreement buffer is lent from the existing rank workspace and returned by Drop even on errors; failure flags remain per-solve. Both workspace paths reuse previous Ritz values; cold diagonal gather/reorder buffers now persist too. The cited root-broadcast allocation is in the allocating reference path; its production workspace counterpart was already persistent. Output-owned history/flags, layout creation, initial-guess and numerical-library allocations remain: this is not a whole-solver zero-allocation claim. Enable per-rank named-site profiling with `set_scratch_allocation_profile(true)` and read `scratch_allocation_stats()` after each solve; disabled costs one entry check, no clocks/atomics/per-collective accounting. No timing or total SCF energy claim; master must measure dev08. |
| ISSUE17-FIRST-BLOCK-DIAG | Diagnostic correction verified; consumer scope extension approved | Synthetic occupied-block fixture: 2 occupied / 4 electrons / 3 mandatory / 5 bands, HOMO–LUMO gap 3 Ha; occupied singular values 0.25 and 0.125 Ha unchanged when first empty residual increases from 1 to 100 Ha. Strict and standard 3-step synthetic SCF: energy/density bits and complete solver results (including work counters) identical with diagnostics off/on; corrected snapshot metadata checked. Focused tests 2 passed; library/focused-test clippy passed. Orchestrator approved the three consumer files; serializer and test references updated without weakening numerical checks. Full offline engine tests passed (exit 0; 195 passed, 70 ignored), including threads/simulated-M2 instrumentation identity and snapshot replay; CLI tests passed (exit 0; 126 passed, 40 ignored). `scripts/check.sh dft-engine` passed build/test/fmt/clippy (exit 0). README regeneration produced no tracked diff. | Physical occupations now come from ScfProblem, independent of mandatory roots. Residual verdicts and absolute spectral-threshold counts are absent for Ritz-change solves; entry Ritz changes are unavailable. Historical first-block report marked superseded; no dev08 table recomputation, optimization decision or performance claim. |
| SPIN-SHARE-DEFAULT-ON | 2026-11-17 | BAANDI_SPIN_SHARE default on. dev08 (master, 10-step screens): si8 -18.1 %, si16 -23.8 % (operator applications exactly halved, bit-identical energies and density); fallbacks verified clean: mpi (all maps fall back, no gain and no harm for rank runs), strict_reference, unequal_sector (total_moment 2.0). | ADOPTED | Gain only for single-process runs; MPI users get none until sharing is implemented for ranks. |
| RETURN-GRAM-LOCKSTEP | Implemented (Mac stage 0); dev08 confirmation pending | Cause of the +6..11% residual of `7320cd8` (grows with B: si16 +2.8%, si32 +5..8%, si64 +11% on dev08): `return_gram_defect` formed C^T S C as one compensated (Kahan) scalar sum at a time, m^2 B latency-bound work on one thread, replicated on every rank. Mac si32 T=8 (m 136..244, B 68): 2.5..8.8 ms per lean return vs about 5 ms for the N x B GEMM it replaced; 8.8 ms = 7% of the 118 ms warm solve. Fix: same per-entry operation sequence, 16 rows of P = S C in lockstep (reads row k of the bitwise-mirrored S, vectorizes), 8 columns of C^T P in lockstep, columns split over the solver's WorkerPool (rank pool in the distributed lean state). Screen 8.8 -> 1.55 ms (T=1) -> 0.29 ms (T=8); test `return_gram_lockstep_matches_scalar_reference_bitwise` (defect bits and decision equal to the scalar form, T 1/3/4). Mac si32 T=8 warm solve, median of 6 per-process floors (load 9..11): a430f8b 118.95 ms, 405b729 122.42 (1.029), fix 113.90 (0.958); eigenvalue FNV and Davidson/H counters identical in all three. T=1 inconclusive (spread 15%). Lesson: an "O(m^2 B) replicated" check is not free when it is scalar and serial: keep reduced-space checks vectorized and on the pool. |
| RETURN-GRAM-DISTRIBUTED | Implemented (Mac, simulated ranks); dev08 confirmation pending | Residual after RETURN-GRAM-LOCKSTEP (dev08, master): T=8 single process none (eecd57f -2.48 % vs a430f8b), MPI si64 np=8 +1.92 % (si32 np=8 +0.75 %) with less counted work: the lean-return screen C^T S C is O(m^2 B) replicated on every rank on one core and does not shrink with np, while the direct N x B Gram it replaced was O(N B^2 / np). Mac probe (simulated np=8, contended): 1.5..5.7 ms per return per rank at si64 (m 146..364, B 132). Fix (distributed lean path only): each rank forms the unchanged lockstep column sums for its own contiguous column range (`return_gram_range`, equal work m^2 + m (j+1) per column) and sends (defect, amplification) in its own two slots of the existing residual-norm reduction (payload B + 2 np; zeros elsewhere are exact): maxima, decision and defect bit-identical to the replicated screen (test `return_gram_column_ranges_reproduce_replicated_screen_bitwise`, ranks 1..7 incl. non-finite S). Flagged returns take one extra reduction (combined, then the unchanged direct Gram); serial path unchanged. Mac simulated np=8, stage `residual` (contains the check), median per kind, eecd57f -> new: si64 cold 4.65 -> 3.64 ms, warm 14.20 -> 11.49 ms; si32 cold 2.24 -> 1.68, warm 3.31 -> 3.07 (load 6..30: walls are noise, not judged). Eigenvalue FNV, iterations, H columns, GEMM calls/FLOP identical to eecd57f in all cells; reduction doubles +2 np per return (+B+2np on a flagged return). Rejected first (same row, not repeated without new evidence): C^T S C through the GEMM workspace. Plain GEMM needs an m eps a allowance (a = max |C|^T|S||C| measured 10..1e5), which sent ordinary si32 warm returns (m 244, a 11) to the direct Gram (margin 1.31); a blocked form (GEMM inside 32-term blocks + Kahan over blocks, rigorous allowance 1.01 (32+6) eps a) kept every probed decision but was no faster than the lockstep kernel: per call 0.74..1.32 x (si64 np=8 354 vs 372 ms per run, si32 np=8 106 vs 104, si32 T=8 2.55 vs 2.51). Lesson: a replicated reduced-space check must be split over ranks like the N-sized work it replaces; the vectorized compensated kernel is already at GEMM speed for these shapes. Design: docs/design/mpi-overlap-davidson.md "Return check distributed over the ranks". |

## ISSUE-20 replicated M2 construction (2026-11-17; Mac observations O on simulated ranks, counters exact C; branch commits 1aefb43 part 1, ef81fd8 part 2, base main eecd57f)

| ID | Change | Status | Numbers | Lesson / next |
|---|---|---|---|---|
| I20-LAZY-SERIAL | `BAANDI_M2_LAZY_SERIAL=on` (default off until dev08): an M2 real Gamma solve state owns only the ORC half-sphere descriptor (`Arc<HalfSphere>`) and the distributed operator; `RealSolveState::op` is `Option` and the serial `RealHamiltonian` (pair plan(s), real projector table, kinetic vector) is built in place only when a branch applies it (distributed operator unavailable -> M1 fallback; threads path). Supersedes the light state on the M2 path when on | IMPLEMENTED (stage 1 O); dev08 confirmation pending | (C) serial operator builds in the lazy arm: 0 (test `dft-engine/tests/m2_lazy_serial.rs`, 2 ranks cubic and skewed; the M1 fallback via a test hook builds it lazily and matches the eager arm bit for bit). (O) si16 np=4 production + ladder, 3 pairs: bit-identical (energy bits, density FNV, every counter), `m2.real_state` 0.002 -> 0.000 s per rank, wall B/A 0.981 inside spread; si32 np=8, 3 pairs: bit-identical (energy 0xc05fad02a3080b2a, 11 SCF, 1888 H columns), `real_state` 0.010 -> 0.002 s, wall 2.063 -> 2.034 s (-1.4 %, spread 0.21 s, load 3). Binary 4666427e... | Removes what the light state left (per-rank real projector conversion, one pair plan, kinetic vector, potential copy per map). Expected dev08 move: `solve.real_state` of the fine and coarse buckets (light state left about 0.011 s fine O) and a little `init`; small, adopt by the adoption rule if dev08 median improves |
| I20-OWNED-ROWS | `DistributedRealHamiltonian::new` packs only its owned kinetic and projector rows (`pack_projector_rows`) instead of a full ORC kinetic vector and a full-length temporary per projector; validation over all representatives kept (rank-independent error) | DONE (always on, bit-identical, no switch) | (C) unit test: owned rows equal the full packing bit for bit on 5 row subsets, same verdict on invalid columns; end-to-end si16 np=4 old/new code: energy bits 0xc04f9085028f55bb, density FNV 0xd3f261f6d6905c70, 9 SCF / 52 Davidson / 907 H columns identical. Time not measured separately (inside `solve.prepare`, 0.014 s si32 np=8 O) | Per-rank memory traffic of the construction falls by (projectors + 1) x n doubles; the validation scan (one norm per representative and projector) stays replicated by design. Parts 3 (geometry/potential/plan descriptor split with invalidation) and 4 (distributed ionic-local) of ISSUE-20 not started |

## Replicated reduced eigenproblem on the lean MPI path (REDUCED-EIGEN-LEAN, 2026-11-18; Mac observations O on simulated ranks, pure build (faer EVD), load 4 to 6, counters exact C, relayed dev08 D; base main 0c4d652)

| ID | Change | Status | Numbers | Lesson / next |
|---|---|---|---|---|
| REL-FINDING | Why `solve.reduced_eigen` grows np4 -> np8 (D, master eecd57f, 10-step fine screens: reduced_eigen +21 % si32, +22 % si64 at constant work) | finding | (D) per-unit-work slowdown np4 -> np8 from the same relayed rows: reduced_eigen 1.21 (work constant), rotate 1.32 (0.66 for an ideal 0.5), reduced_gram 1.36, kept_update 1.24, apply 1.18: the EVD is NOT anomalous, every core runs about 1.2 to 1.36 x slower per unit of work with 8 busy cores (Xeon Gold 6230: AVX-512/all-core frequency and shared L3/memory bandwidth); it only stands out because its work does not shrink. dft-mkl links `mkl_sequential`, so the EVD is single-threaded per rank whatever `MKL_NUM_THREADS` says (no oversubscription). The lean path (default) ignores `BAANDI_M2_EIGEN=root` (every rank solves). (O, si32 np=4 warm, temporary timers) stage 1 = 83 % dense EVD (26.1 ms of 31.3), 7 % `extend_factor`, 9 % `reduced_border`; inside the EVD call 7 % was validation/sign overhead (1.51 of 22.27 ms). (O, faer) eigenvalues-only (tridiagonalization + values) is 52 to 62 % of the full EVD at n = 136..600, so a distributed eigenvector/back-transform stage can remove at most about 40 % x (1 - 1/P); MKL-EVD measured dsyevr(b roots) slower than dsyevd, so on MKL the replicated reduction is probably a larger share (microbenchmark REL-SLICES). A distributed Householder tridiagonalization needs about 2 m small collectives per EVD (400 to 800 at m = 200..400): more latency than the EVD itself at these orders, not pursued | Shrink or distribute the EVD; do not chase the growth as a contention bug. Optional dev08 probe of the frequency part: the same np4/np8 rows with `MKL_ENABLE_INSTRUCTIONS=AVX2` (rounding class, diagnostic only) |
| REL-CHECKS | `RealEighWorkspace::solve_into` validation in one tiled pass (finiteness, max |a|, max asymmetry over 32 x 32 tile pairs instead of three sweeps with a full-stride mirror walk) and one read sweep for output finiteness plus the sign pivots (flips only after the whole block is finite) | DONE (always on, bit-identical, no switch) | (C) unit test `real::tests::reduced_eigen_checks_match_former_sweeps`: 90 cases (n 1..100 around the tile size, asymmetric, NaN/inf anywhere, ties, signed zeros, 1e300 scale) give the same error kinds, defect/allowance bits and output bits as the former code kept verbatim in the test; si32 np=4 bench eigenvalue FNV unchanged (cold 0x4a8a8c3b366b6701, warm 0x35da8f77685b506e). (O) overhead around the raw EVD per warm solve 1.51 -> 0.60 ms, cold 0.25 -> 0.09 ms (about -4 % of stage 1 on faer) | MKL EVD is about 0.55 x faer, the overhead is the same absolute work, so the expected dev08 share is about 5 to 8 % of the stage before, 2 to 3 % after (E). Applies to every caller of `RealEighWorkspace` (serial Davidson too) |
| REL-SLICES | dev08 microbenchmark `crates/dft-mkl/examples/evd_slices.rs` (MKL, one core): `dsyevd` against the sliced pipeline replicated `dsytrd` + `dstebz` (k lowest values) + per-rank slice of k/P eigenvectors (`dstein` or `dstemr` by index) + `dormtr` on the slice; m = 136..816, k = ceil(0.34 m), P = 1/4/8, matrix models dense and Davidson-like (diagonal Ritz block + dense m/3 border); accuracy of the independently computed, gap-aligned slices assembled at P = 8 (max |Z^T Z - I|, max residual / max |A|) | READY (type-checked with `--features mkl`, not runnable on the Mac); dev08 pending | Decision rule (pre-registered): implement `BAANDI_M2_EIGEN=slices` (default off, rounding class) only if sliced(P=8) / dsyevd <= 0.80 at m = 204 and 408 in both models with the m x k allgather (D collective timers) added, and orthogonality <= 1e-13 for the chosen variant; otherwise close the distributed-EVD lever for orders below 1000 | (E) faer analogue: at most -40 % x 7/8 of the EVD; the MKL ratio decides |
| REL-TREE-RHD | `BAANDI_MPI_TREE=rhd` (new value, default stays `rd`): payloads of at least `dft_mpi::HALVING_MIN_DOUBLES` = 2048 doubles by a recursive-halving reduce-scatter (partners rank XOR 1, 2, 4, ... as rd; the rank with the round's bit clear keeps the lower half of the shared segment) plus a recursive-doubling allgather of the kept segments (copies), smaller payloads rd; power-of-two P only, else binomial | IMPLEMENTED (type-checked with `--features mpi`; FFI not runnable on the Mac); dev08 microbenchmark pending | (C) why tree_sum grows np4 -> np8 at si32: the lean path issues three reductions per Davidson iteration (apply projections, S border, H border; Mac si32 np=8 warm: 15 + 15 + 16 of 50 solver reductions, 2570 doubles per call on average) whose payloads are replicated m x c blocks that do not shrink with P, and rd moves the whole payload in each of log2 P rounds (2 -> 3 rounds: up to +50 % per call at equal calls; the D +28 % at si32 is that, the -18 % at si64 is the wait share shrinking). rhd moves 2 (P-1)/P of the payload (1.75 n at P = 8 instead of 3 n) in 2 log2 P rounds. (C) bit identity: every element is associated as in rd, hence as in the binomial tree (pure model `dft_mpi::halving_doubling_associate`, test `halving_doubling_equals_the_binomial_tree`: P = 1..16, lengths 0..100 incl. empty segments, cancellation-sensitive values); `mpi_tree_bench` bit-checks binomial / rd / rhd / rhd with threshold 0 on every rank. (E) from the relayed rd timings (P = 8, 27744 doubles 133.8 us, about 5 GB/s and 1 us per round): 64 KB border -30 %, 222 KB -37 %, below about 1000 doubles a loss (hence the threshold) | Calibrate `HALVING_MIN_DOUBLES` from the `rhd-forced` column; adopt as default only if the fine-solve tree_sum seconds drop on the si32/si64 np=8 screens with identical bits (no rounding change possible) |

## INIT-SCALE: making `init` scale with MPI ranks (2026-11-18; Mac observations O on 8 simulated ranks, load 3 to 5, pure build; counters exact C; branch base main 0c4d652, commits 03fd272, 5262fcf, a7205c9, 529e485)

Trigger (dev08, master eecd57f, 10-step fixed screens): si32 `init` 0.629 / 0.603 / 0.637 s at np 2 / 4 / 8 while the rest scales about -50 %; si64 atomic start 1.127 s of which orthonormalize 0.849 s (-14 % np2 to np8), ionic_local 0.449 s (+9 %), states 0.314 s (+20 %), ewald 0.274 s (+12 %). Harness: `m2_entry_cost` production + ladder, `BAANDI_PHASE_PROFILE=on`, numbers are max over ranks [min]. Baseline (O, same binary family): si32 np8 init 0.207 to 0.227 s (ionic_local 0.043 to 0.048, states 0.033 to 0.041, ewald 0.034 to 0.040, atomic start 0.085 to 0.092 = orthonormalize 0.043 to 0.047 + Ritz 0.040); si64 np8 init 1.249 s (ionic_local 0.123, states 0.226, ewald 0.205, atomic start 0.654 = orthonormalize 0.328 + Ritz 0.315 of which rotate 0.142).

| ID | Change | Status | Numbers | Lesson / next |
|---|---|---|---|---|
| INIT-ATOMIC-BLOCK | `BAANDI_M2_ATOMIC_BLOCK=on` (default off, rounding class): with `M2_ATOMIC=split`, every rank fills its contiguous row block of the atomic span (same elementwise arithmetic), orthonormalizes by row-distributed CholQR2 (local Gram, one tree_sum, replicated Cholesky + triangular inverse, local update; twice; real GEMMs on the interleaved data at Gamma, where overlaps and factors are real) and gathers the rows once; the Ritz rotation Q C runs by row blocks + one allgatherv. Falls back to the split Gram-Schmidt (counter `atomic_block_fallbacks`) when a pivot is below 1e-6 (CholQR validity bound, sqrt(eps) is 1.5e-8), a rank cannot fill, or data are nonfinite | IMPLEMENTED (stage 1 O); dev08 A/B pending | Cause (code reading, C): in `orthonormalize_split` the second Gram-Schmidt pass of column j runs on its owner before the broadcast, a serial chain of about m n^2 complex MACs (m = span rows, n = atomic orbitals) independent of P; CholQR2 costs about 4 m n^2 / P (2 m n^2 / P real at Gamma). (O) si64 np8: orthonormalize 0.328 -> 0.18 to 0.23 s, rotate 0.142 -> 0.035 s, atomic start 0.654 -> 0.39 to 0.41 s; si32 np8 atomic start 0.092 -> 0.053 s. SCF iterations and Davidson counts unchanged (si32 11 / 67, si64 10 / 51); energy bits differ by 1 to 2 ulp (si32 0x...0b2a -> 0x...0b29, si64 equal), density FNV differs (rounding class). The row-split rotation alone left the density FNV unchanged (faer); MKL not checked, so it stays under the switch. Tests: lib `block_cholesky_qr_matches_gram_schmidt_and_falls_back_together` (P 2..7, complex and Gamma, \|Q_block - Q_gs\| < 1e-12, joint fallback); `m2_coarse_certificate::block_atomic_start_stays_in_the_class_of_the_split_start` (P 2/3/5/7, \|dE\| <= 3.6e-14 Ha vs threads-only, ranks equal, reproducible, no fallback) | Rounding class: adoption needs dev08 time plus the stage-1 class check (and the generality set for a default). Remaining serial pieces of the block path: replicated n x n Cholesky / inverse (n^3/3 each, n = 256 at si64) and the full-row allgatherv (m n complex, same volume as the split broadcasts). Ritz eigh (n x n, replicated) untouched |
| INIT-IONIC-SPLIT | ISSUE-20 part 4: ionic local coefficients by contiguous density-G blocks per rank (`ionic_local_rows` / `ionic_local_with_conventions_rows`, per-G arithmetic of the full call) + one allgatherv with a per-rank failure flag (any failure repeats the replicated call and its error); always on with a communicator | DONE (always on, bit-identical, no switch) | (C) 168 partition cases (1..7 parts; cubic, skewed, elongated, slab; Gamma and general k; Baandi/Simpson 3/8, QE-like/quadratic end, QE-like G0) bit-identical; engine ionic potential bit-identical on 2..7 simulated ranks (`tests/init_scale.rs`); si32 / si64 np8 production + ladder energy bits and density FNV unchanged (0xc05fad02a3080b2a / 0x92de843e24a85b81, 0xc06fac198a25f9b6 / 0x269ac897397fc0f5). (O) `init.ionic_local` np8: si64 0.123 -> 0.031 to 0.067 s, si32 0.043 -> 0.013 to 0.016 s (max over ranks; min 0.008 to 0.024) | Basis order is by \|G\|, so equal G counts per rank also split the radial shells, but low-\|G\| blocks hold more distinct shells (more radial integrals): the max includes waiting for those ranks. A shell-weighted split would balance it; the FFT of the assembled coefficients (`densities`, replicated full grid) is small |
| INIT-EWALD-SPLIT | Ewald sum split over the ranks, bit-identical for every rank count: `ewald::EwaldPlan` numbers the real-space candidates (pair i, j, then the translations in the serial enumeration order) and the reciprocal candidates; every rank evaluates the accepted terms of its contiguous block; the two running sums pass down the ranks (rank r folds its terms into the sums of rank r - 1, P broadcasts of three values incl. a failure flag), which is the serial association exactly. Serial `energy()` refactored onto the same prelude and term functions | DONE (always on, bit-identical, no switch) | (C) verbatim copy of the previous evaluator vs `energy()`, concatenated and chained plan, 1..7 parts, 16 cell/ion cases (cubic, skewed, elongated, slab; 1..9 ions, mixed charges, positions outside [0,1)): all nine output fields bit-identical; coincident ions fail in both; engine Ewald fields bit-identical on 2..7 simulated ranks; si32 / si64 np8 end-to-end bits unchanged. A first version gathered every accepted term (allgatherv): si64 0.205 -> 0.124 s only (about 2M real-space terms, 16 MB per rank, r_cut about 100 Bohr at tolerance 1e-13) - replaced by the chained sums. (O) `init.ewald` np8: si64 0.147 to 0.205 -> 0.040 s, si32 0.034 to 0.040 -> 0.013 s | The EWALD_CACHE (default on) already removes the coarse ladder's second evaluation; this removes the replication of the remaining one. The serial chain is P broadcasts of 3 doubles (latency only) |
| INIT-NONLOCAL-TABLE | (not implemented; analysis) `states` = per-k complex `NonlocalOperator` table (all atoms x channels x full sphere) + `Hamiltonian::new` + band-density workspace | OPEN (backlog) | (O, temporary probe, si64 np8): NonlocalOperator 0.07 to 0.18 s fine (npw 23871) and 0.075 to 0.095 s coarse (npw 8409) per rank, Hamiltonian::new 0.01 to 0.05 s, band density 0.001 to 0.006 s. The spread for identical work across ranks points to memory contention (table about 100 MB per rank at si64 fine) | The table is per-G independent like the ionic local (bit-identical split possible), but a split + allgatherv writes the full table on every rank plus a gather buffer (about 3x the memory traffic), so the gain is doubtful. The scalable fix is not materializing the full complex table on M2 ranks (owned rows only; consumers: the distributed operator packing, split_gram's complex H apply of the atomic start, audits and fallbacks): ISSUE-20 part 3 territory |
| I20-LAZY-SERIAL-DEFAULT-ON | 2026-11-18 | BAANDI_M2_LAZY_SERIAL default on. dev08 (master, 95e1242 family, 10-step fixed): si64 np=8 -1.22 % (floor 0.09 %), fine solve -4.5 %, solve.real_state 0.1286 -> 0.0072 s (-94.4 %); si32 np=8 -0.8..-1.4 % (inside the noisy floor, real_state -89.9 %, only 1 % of the run); counters, energy and density SHA identical. | ADOPTED | Gain is larger for larger cells; no harm on si32. |
| REL-SLICES-CLOSED | 2026-11-18 | M2_EIGEN=slices not built. Pre-registered rule (<= 0.80 x dsyevd at m=204 and 408, both matrix models, orthogonality <= 1e-13): dev08 evd_slices reps 20, 8 cores busy. sliced-mrrr time 0.43-0.53 x but orthogonality 1.7e-13 / 6.9e-13 / 2.7e-13 in 3 of 4 cells (worse with m); sliced-stein orthogonality 1e-14 but time 1.42-1.64 x. | REJECTED | Distributed dense EVD lever closed below m ~ 1000. |
| REL-TREE-RHD-THRESHOLD | 2026-11-18 | BAANDI_MPI_TREE=rhd threshold made communicator-size dependent (55000 doubles at <= 2 ranks, 4096 at 3-4, 8192 at >= 5) from the dev08 mpi_tree_bench crossovers; the fixed 2048 was 28-59 % slower than rd. rd stays the default until an A/B of rhd with the new thresholds shows a gain (tree_sum is only ~1 % of the run, so the lever is <= 0.4 %). | PENDING A/B | Bit-identical to rd for any threshold. |

## INIT-SCALE-2: replicated init pieces that grew with ranks (2026-11-18; base main e902853, commits 29bfc48, 65b78ad, 6220dbc, c639337; Mac observations O on simulated ranks, load 5 to 7; bits C)

Trigger (dev08 master e7c94c8, 10-step screens, np2 -> np8): init scales (si32 -20 %, si64 -29.5 %) but its replicated parts grow: states +20 % (si64 np8 0.3215 s), cli +33 %, rest +29 %, basis +24 %, densities +21 %; distributed Ewald at np=2 83 to 102 % slower than replicated (si32 0.0595 -> 0.1089 s, si64 0.2451 -> 0.4949 s), 48 % faster at np=8. Probe (O, si64 np8, per construction, max over ranks): `states` = nonlocal table 0.04 to 0.09 s (of which the P48C packed copy about 1/3), Hamiltonian::new 0.02 to 0.04 s (a second density sphere), band density 0.002 s; `rest` = Coulomb metric 0.02 to 0.04 s (gamma basis + a third density sphere); `basis` = density sphere 0.009 to 0.019 s. Identity reference (base and head, production + ladder harness): si16 np3 0xc04f9085028f4f47 / 0xa1df51b1440229b2, si32 np2 0xc05fad02a3080e0e / 0xa34b5d1bae989d97, si32 np4 0xc05fad02a3080b24 / 0x561c9f9b5d8dbdf7, si64 np2 0xc06fac198a263f2c / 0x01c4b7936315215c, si64 np8 0xc06fac198a25f9b6 / 0x269ac897397fc0f5 (energy bits / density FNV).

| ID | Change | Status | Numbers | Lesson / next |
|---|---|---|---|---|
| INIT-EWALD-GATE | `ewald_evaluate` splits only for communicators of at least `EWALD_SPLIT_MIN_RANKS` = 4 ranks (named constant with the dev08 numbers); below, every rank evaluates the serial sum, no collective | DONE (always on, bit-identical either side) | (C) `tests/init_scale.rs`: P = 2..7, split counted exactly for P >= 4, Ewald and ionic bits equal to the serial construction. (O, Mac simulated ranks, 3 pairs) si64 np2 `init.ewald` 0.089 -> 0.182 s, i.e. the OPPOSITE of dev08 | Suspicion to check on dev08: the split's chained broadcasts are the first collective after the replicated `states`, so `init.ewald` (max over ranks) absorbs the rank skew of `states` (dev08 probe: identical work 0.07 to 0.18 s across ranks). Judge the gate by `init` total and the next collective bucket (`atomic.*`), not `init.ewald` alone; if init total is not better at np=2, set the constant to 2 |
| INIT-PACKED | `NonlocalOperator` stores all atoms' projector columns once, in one column-major table shared by `Arc`; `AtomProjectors::columns()` returns a borrowed view; the P48C block path multiplies the table itself (`set_block` sets a flag, no second packed copy); constructors fill the table directly (same products) | DONE (always on, bit-identical, no switch) | (C) `setup_memo_identity` (verbatim old constructor) and the P48C block / traversal tests unchanged; five end-to-end bit pairs above identical. Memory: one projector table less per rank and k state (si64 fine about 98 MB) when the block gate is on | The remaining `states` cost is the per-atom phase factors exp(-iG.tau) (64 x 23871 at si64), compute-bound. A split over atoms + allgatherv would remove that compute but write the table twice more (gather buffer + conversion, no safe f64 view of Complex64 without a new dependency): not done. Removing the full complex table from M2 ranks needs the atomic-start Gram from row-distributed projections (P^H D P, rounding class) and lazy full tables for audits (ISSUE-20 part 3) |
| INIT-REUSE | `Hamiltonian::new_with_grid` takes the engine's k grid (checked in O(N_G) against `map_miller`) instead of building the density sphere again; `CoulombMetric::with_sphere` takes the construction's density list (same vectors as `options.basis(..).density_basis(..)`, strain reference included); the exact Gamma point reuses the Gamma basis | DONE (always on, bit-identical) | (C) lib test `new_with_grid_matches_new_and_checks_the_grid` (Gamma and general k: apply bits equal; other dims or another basis rejected); five bit pairs identical. (O, si64 np8, 4 alternating pairs vs e902853, INIT-PACKED + INIT-REUSE) `init.states` 0.121 -> 0.070 s (-42 %), `init.rest` 0.029 -> 0.006 s (-80 %), init -10.6 % | Density spheres per construction 3 -> 1. Hamiltonian::new now 2 to 4 ms |
| INIT-CLI | `dft-cli run` reuses the preflight's `Config::prepare` result (`Report::prepared`, not serialized) instead of preparing again | DONE (always on, bit-identical) | (C) CLI tests preflight / preflight_findings / run pass. Saves one density sphere, two wavefunction bases and one grid per rank (si64 about 0.01 to 0.02 s O, not timed end to end on the Mac) | The engine still builds its own bases (API change needed to pass the CLI's); `structure::review` and UPF parsing remain replicated and were not profiled |
| INIT-DENSITY-SPLIT | Opt-in atomic-superposition start density (`start_density = atomic`, also the relax predictor): coefficients by contiguous G blocks (per-G arithmetic of the full evaluation) + one allgatherv, failure flag -> replicated call; counter `density_splits` | DONE (always on with a communicator, bit-identical) | (C) `tests/init_scale.rs`: start-density grid bits equal on P = 2..7 (cubic and skewed, two local-mesh conventions), split counted | Default runs use the uniform start: there `densities` is the FFT plan + one transform of the ionic potential (O 0.011 s si64 np8), left replicated |
| INIT-SCALE-2-NOT-DONE | Analysis only: density sphere enumeration (`PlaneWaveBasis::at_k`, 191735 vectors at si64 4 x 15 Ha) 7.1 ms warm, 18 ms cold on one core; `sort_unstable_by` 7.7 ms (no gain, reverted). `between.session_setup` (O 0.23 s si64 np8, outside init) was not examined | - | (O) | Pruning the box enumeration (42 % of the box is inside the sphere) could save about 3 ms per sphere; not worth the risk |
| INIT-SCALE-2-MAC | Mac O, si64 production + ladder, 3 alternating pairs, base e902853 vs head c639337 (bits identical) | O only | np8: init 0.789 -> 0.687 s (-13.0 %), states 0.105 -> 0.058 (-44.9 %), rest 0.025 -> 0.005 (-79 %), basis -13.9 %, ewald +7.7 % (noise; still split), total wall -3.9 %. np2: states -26 %, rest -80 %, ewald +104 % (gate, see INIT-EWALD-GATE), init +17 %, total +7.5 % (load 7.4) | dev08 decides (docs/performance/requests/init-scale-2.md) |
| M2-STAGE-TIMERS | 2026-11-18 | Instrumentation only (no speed claim). `BAANDI_DAVIDSON_STAGES=on` now fills `convergence.work.davidson_stages` on the MPI M2 path: `stage_s` from the distributed solver stages (reduced_gram -> reduced, reduced_eigen -> eigen, rotate, residual, orthonormalize, apply_total -> h_apply, precondition; solver_total counts timed_solves), `h_apply_parts` from the distributed operator (`local_s` = mirror exchange + pack + distributed pair FFT with transposes + mirror-h exchange + unpack + kinetic; `nonlocal_projection_s` = local P^T X + tree sum; `nonlocal_back_s` = coupling + P W + image), `fft.pair_*` / `fft.density_*` / `fft.kept_*` counts and seconds, and new `collectives.{all,in_h_apply}.<op>.{calls,s}` per run stage. Clocks need `BAANDI_STAGE_TIMERS=on` too (off: counts only). MPI: every number is the max over ranks, `min_over_ranks` the min (one allgatherv at the end, only with the switch on). Mac test `m2_stage_timers` (2 simulated ranks): 174 applies, apply 0.0591 s = local 0.0494 + projection 0.0072 + back 0.0025 (remainder 8e-6 s), 352 FFT pairs, 348 alltoall and 174 tree sums inside H; energy and density bits equal on and off. Off cost: one relaxed load per block apply / stage entry / collective, no clock. | DONE (diagnostic) | dev08 si32 np=8 run with STAGES+STAGE_TIMERS+PHASE_PROFILE gives the local / projection / back split per run stage. Simulated ranks share the process counters (and the process-wide run-stage index), so per-rank values and the stage split are exact only under real MPI. |
## Between, verification glue and materialize (BETWEEN-MATERIALIZE, 2026-11-18; Mac observations O on 8 simulated ranks, pure build, loaded Mac: buckets only, walls not comparable; counters exact C; base main e902853, commits 5fe2796, b3f8509, 3f7d82b)
Trigger (dev08, master e902853, 10 fixed steps, np=8): terms with no QE counterpart are 13-15 % of the run: si32 verification 0.264 + materialize 0.038 + between 0.119 s; si64 verification 0.532 + materialize 0.139 + between 0.652 s, and between grows +46 % from np=2 to np=8 for si64. What the buckets hold (O, si64 np=8, base): `between` 0.271 s = session_setup 0.214 (fine-session import of the prolonged ladder state; 0.19 s of it is the orthonormality check of the imported N x B block, a scalar O(N B^2/2) complex Gram loop replicated on every rank: N 23871, B 132), stop_rule 0.036 (density residual, charges, min density, kept-coefficient Coulomb estimate: full-grid passes, about 3 ms per step), ladder_handoff 0.025 (orbital prolongation, density zero padding with two FFTs, ElectronicState copies), verify_glue 0.003; `materialize` 0.252 s = two lazy gathers (coarse block at the hand-off, N 8409; the delivered fine block, N 23871): allgatherv 0.034 s, scatter into an N x B ORC array 0.05-0.12 s, per-column unpack 0.10-0.14 s per rank. `verify_pass` is one ordinary map at the delivered density (ADR 0043), unchanged here.

| ID | Change | Status | Numbers | Lesson / next |
|---|---|---|---|---|
| LADDER-IMPORT-SPLIT | `BAANDI_IMPORT_CHECK=split` (default; `serial` = former loop): the import check of `initialize_electronic` computes the same per-entry overlaps (sequential row sums from zero, bit-identical) eight columns i in lockstep against one column j, and with a registered communicator only the columns j mod P of the rank, one 1-double tree sum of the failure flags (same verdict on every rank) | DONE (default on, bit-identical; dev08 confirmation pending) | (O, si64 np=8, same binary A/B) check 0.19 s -> 0.02-0.03 s per rank; `between.session_setup` 0.209 -> 0.048 s; `between` 0.266 -> 0.106 s. si32 np=8 session_setup 0.032 -> 0.011-0.017 s. Energy bits and density FNV identical (si32 0xc05fad02a3080b2a / 0x92de843e24a85b81, si64 0xc06fac198a25f9b6 / 0x269ac897397fc0f5). Lib test `split_import_check_matches_serial_verdict` (orthonormal, defects above/below 1e-8, NaN, inf; 1-5 simulated ranks) | The replicated, latency-bound validation of a hand-off was the biggest `between` item and explains its growth with np (8 ranks running the same scalar loop share the memory system). Any replicated O(N B^2) check is a split candidate |
| LAZY-ASSEMBLE-FUSED | `LazyOrbitals::materialize` writes each gathered ORC coordinate straight into its two complex slots (x / sqrt2, mirror imaginary negated: the arithmetic of `HalfSphere::unpack`), column by column across the rank segments, without the N x B ORC intermediate; same errors (length, nonfinite) | DONE (always on, bit-identical, no switch) | (O, si64 np=8) `materialize` 0.31 -> 0.18 s (with the delivered gather still on); lib test `fused_assembly_matches_scatter_then_unpack` (cubic, skewed, elongated; 1-5 uneven segments, permuted rows, signed zero, subnormal; same error strings) | One pass and one 25 MB intermediate less per gather; the allgatherv and the output block's first touch remain |
| LAZY-RESULT | `gamma_real::set_delivered_orbitals(false)` + `BAANDI_M2_LAZY_RESULT=on` (default on; off = always gather): an M2 result keeps its lazy blocks when the caller reads no orbitals. dft-cli `run` declares it without `output.forces`, `output.stress`, a volumetric quantity other than density/hartree/local_potential, observe/agent mode (replicated: same input on every rank); `ScfResult::materialize_orbitals` (collective) gathers on demand, `orbitals_delivered()`; forces and stress refuse a lazy result; the coarse-terminal ladder branch gathers itself; the energy-identity audit still gathers. Counter `lazy_results` | DONE (default on, bit-identical outputs; dev08 confirmation pending) | (O, si64 np=8) `materialize` 0.18 -> 0.018 s (only the coarse hand-off gather left), si32 0.007 -> 0.003 s; energy bits and density FNV identical. Test `tests/m2_lazy_result.rs` (cubic and skewed, 2 and 3 ranks: energy, density, eigenvalues, history identical; materialized orbitals bit-equal to the gathered run; forces refuse). Consumer inventory: default `run` output reads eigenvalues and density only | Library callers keep the old behaviour (default true). Relax, bands, checkpoint never set the flag |
| BETWEEN-REST (not done) | stop_rule (about 3 ms per step: separate full-grid passes for residual, two charges, min density) and ladder_handoff (prolongation + ElectronicState round trip: 2 x 50 MB copies at si64) | OPEN (backlog) | (O, si64 np=8) stop_rule 0.031-0.039 s, ladder_handoff 0.023-0.032 s per run | Fusing the stop-rule passes is bit-identical (separate accumulators, same order) but worth about 1-2 ms per step; a matrix-form hand-off (no coefficient Vec) saves two copies. Both under 0.5 % of the run |
## Distributed full-grid stages of the M2 path (M2-GRID, 2026-11-18; Mac observations O on simulated ranks, pure build, load 5 to 7, counters exact C; base main e902853, commits 414c5e2, f0f45ce, ecc1ad2; dev08 request docs/performance/requests/m2-grid.md)

Trigger (master, dev08, e902853, per SCF step, np=8, MPI against MPI): outside the solve `potential` 0.0228 / 0.0522 s, `mixing` 0.0124 / 0.0262 s, `energy` 0.0208 / 0.0482 s (si32 / si64), replicated full-grid work that does not shrink with ranks (mpi-remaining.md section 2).

| ID | Change | Status | Numbers | Lesson / next |
|---|---|---|---|---|
| M2-GRID-MAP | Which replicated loops dominate (temporary timers, not committed; si32 production np=8, median per map) | finding (O) | potential 7.2 ms: Hartree FFT pair of the input density (replicated complex 3D FFT plus a full-grid G^2 rebuild) about 6 ms, input XC (distributed + allgatherv) about 1 ms, effective sum 0.4 ms; energy 7.3 ms: output-density Hartree FFT 3.0 ms, full-grid audit sums 1.1 ms, output XC 0.7 ms (its potential allgatherv is never used by the audit), lazy kinetic + distributed nonlocal energy 2.2 ms; mixing 2.0 ms (Pulay window, incremental) | The FFTs, not the XC, are the replicated cost; the output XC potential gather was pure waste |
| M2-GRID-POTENTIAL-AUDIT | `BAANDI_M2_GRID=on` (default off, rounding class): slab plan = M2a `DistributedPairFft` with every grid point as an index; potential (forward transform of the owned planes, 4 pi/G^2 on owned sticks, inverse, XC and effective on the planes, ONE allgatherv of the effective potential plus partial sums), audit (output transform, charge, double counting, Bregman linear terms, local energy in real-space form, output XC energy; ONE allgatherv of partials, no XC potential gather), stop rule from kept owned coefficients (one small allgatherv). Canonical per-plane / per-stick sums | CANDIDATE (Mac; dev08 Stage 1) | (C) lib tests ranks 1..7, cubic and skewed, odd/even grids, PZ/PW: identical bits for every rank count, <= 1e-12 relative to the replicated physics routines; engine production SCF cubic/skewed P 2/3: abs dE <= 4.4e-15 Ha, same SCF counts. (O) si32 np=8 alone: wall 2.388 -> 2.221 s (-7.0 %, spread 0.13 s), energy 0.105 -> 0.043 s, potential 0.108 -> 0.072 s; residual per map 2.6 ms (forward 0.50, inverse 0.45, XC+effective 0.28, gather 0.30 ms) | Rank-count independent bits are possible here because the transforms never change per-line arithmetic; only the sums needed a canonical order. Not done: keeping V_eff on the planes for the distributed operator (it still copies the full potential from the complex Hamiltonian; one allgatherv of N doubles per map remains) |
| M2-GRID-MIX | `BAANDI_M2_MIX=on` (default off, rounding class): engine-side `DistributedPulay` (dft-solvers untouched) for plain Pulay on the M2 path: window on the rank's block of fixed 4096-point chunks, newest Gram row + rhs as chunk partials (one allgatherv), replicated eigenproblem, two charge corrections (one allgatherv each), one allgatherv of the new density | CANDIDATE (Mac; dev08 Stage 1) | (C) lib test ranks 1..7 against PulayMixer with seeds, window sliding, reset: <= 1.3e-12 on densities of order 1 (no growth over 10 fed-back steps), identical bits for every rank count; engine SCF: abs dE <= 4.4e-15 Ha. (O) with GRID: si32 np=8 wall 2.463 -> 2.241 s (-9.0 %, spread 0.08 s), mixing 0.053 -> 0.015 s; si32 np=4 2.735 -> 2.593 s (-5.2 %), potential 0.077 -> 0.044, energy 0.086 -> 0.049, mixing 0.036 -> 0.015 s; si64 np=8 9.320 -> 8.606 s (-7.7 %, spread 0.23 s), potential 0.242 -> 0.136, energy 0.355 -> 0.174, mixing 0.094 -> 0.026 s; abs dE 2.8e-12 / 9.3e-12 Ha, SCF/Davidson/H identical | Exact checkpoints reject the mixer (opaque state); Hartree-metric and Kerker mixers stay replicated |
| M2-GRID-PLAN | Plan construction: stick owners from two small allgathervs of owned stick ids instead of `position_owners` (a second full partition) | done (always on inside the switch, identical bits) | (O) si32 np=8 per rank: coarse plan 11 -> 6 ms, fine 23 -> 11..22 ms (two plans per run) | The remaining cost is `DistributedPairFft::new`'s own partition of the full index set (dft-solvers, not touched) |
| M2-GRID-MIX-DEFAULT-ON | 2026-11-18 | BAANDI_M2_GRID and BAANDI_M2_MIX default on. dev08 (master, 6be8f09, 10-step fixed, 3 pairs, A defaults vs G=GRID vs GM=GRID+MIX): wall GM si32 np8 -15.12 %, si64 np8 -9.17 %, si32 np4 -7.83 %, si64 np4 -4.24 % (G alone -10.51 / -6.51 / -5.83 / -3.34 %); potential -28..-45 %, energy -40..-66 %, mixing -62..-78 %; |dE| <= 8e-12 Ha (standard class 5e-9), SCF counts identical, verification passed, all 8 cells improve. | ADOPTED | Rounding class, fixed-order sums (bits identical for every rank count); EWALD_SPLIT_MIN_RANKS lowered to 2 in the same push (np=2 gated serial path measured slower: si64 init.ewald 0.4977 split vs 0.8178 gated). |

## Ewald and ionic local construction (INIT-EWALD-FAST, 2026-11-18; Mac observations O, release, one thread, examples `ewald_bench` / `ionic_bench`; counts exact C; dev08 request docs/performance/requests/init-ewald-fast.md; base main 9d942a0)
| ID | Change | Status | Numbers | Lesson / next |
|---|---|---|---|---|
| INIT-EWALD-SLOWDOWN | Cause of the dev08 `init.ewald` regression (eecd57f replicated about 0.245 s -> 0.813 s at one core, si64) | analysed | (C) si64 reference scheme: 28.1M box candidates (19^3 per pair x 4096 pairs), 3,865,600 accepted erfc terms, the same before and after INIT-SCALE-3 (529e485). (O) eecd57f evaluator 0.069 to 0.080 s, 529e485 serial `energy()` 0.139 s (2.0x; dev08 3.3x), bits equal | The refactor routed every one of the 28M candidates through `real_term` (bounds-checked position fetch and delta per candidate, `Result<Option<f64>>`, closure through `enumerate`); the erfc work was unchanged. A per-candidate refactor of a 10^7-candidate loop is not free: measure the serial path after any "pure refactor" of a hot loop |
| INIT-EWALD-PRUNE | Reference scheme (default): every ion pair's translations are enumerated only near the cutoff sphere (`ImageWalker`: necessary conditions per axis from b_k . A x = 2 pi x_k, per-row quadratic interval for m2, each widened by one index and a 1e-9 relative radius slack), the exact old predicate on every visited candidate, same order; plan candidates are (pair, m0 row); forces and stress use the same walker | DONE (always on, bit-identical, no switch) | (C) lib `ewald::tests::reference_scheme_matches_the_previous_evaluator_bit_for_bit`: energy (9 fields), every partition 1..7 of the plan, forces and stress bit-identical to a verbatim copy of eecd57f on 16 cell/ion cases; engine `tests/init_scale.rs` P = 2..7. (O) si64 energy 0.069 s (eecd57f 0.069 to 0.080, 529e485 0.139); forces 0.151 s | The reference scheme is now erfc-bound (3.9M calls, about 17 ns each); no further bit-identical gain is available because eta = sqrt(pi)/Omega^(1/3) with the 2 rho padding puts r_cut at 125 Bohr for 64 ions |
| INIT-EWALD-FAST | `BAANDI_EWALD_FAST=on` (default off, rounding class): `Scheme::Balanced`: eta = sqrt(pi) (N c)^(1/6) / Omega^(1/3) from the real/reciprocal cost balance (c = `BALANCE_COST_RATIO` = 16, a machine cost ratio, scan flat within 10 % for 12..64 on cubic, skewed, elongated and slab cells with 1..216 ions); cutoffs from rigorous lattice-counting tail bounds (circumradius, shell term plus integral) with budget min(tol (sum |Z|)^2/Omega^(1/3), tol N), t stepped by 0.01; ion pairs i <= j (one term per pair); reciprocal half space with per-ion, per-axis phase tables; forces and stress take the same eta and cutoffs; plan split bit-identical for any rank count | opt-in (dev08 pending) | (C) lib tests: against a Neumaier-compensated reference at tolerance 1e-16 (its own eta dependence < 1e-13 Ha/atom): energy <= 8.0e-15 Ha/atom (caps pre-registered: 1e-12 Ha/atom), forces <= 4.9e-14 Ha/Bohr (cap 1e-10), stress <= 8.8e-13 Ha/Bohr^3 (cap 1e-12; limited by the uncompensated stress reference) on cubic, skewed, elongated, slab cells with 1, 2, 8 (debug) and 64 ions (ignored, release); plan partitions 1..7 bit-identical; engine P = 2..7 bit-identical with the switch on. (O) si64: eta 0.274, r_cut 21.7 Bohr (reference 124.9), 2080 pair terms + 2520 G (reference 3.87M + 948); energy 0.069 -> 0.00064 s (108x), forces 0.151 -> 0.0042 s; E_F - E_R = 1.9e-9 Ha (the reference scheme's own summation rounding, see compensated test) | Adopt on dev08 if `init.ewald` drops (expected to about 1 ms) and nothing regresses; the total energy moves by about 2e-9 Ha at si64 (inside the standard class 5e-9 Ha, outside the 1e-12 cross-path check, so strict reference fixtures that store Ewald bits stay on the reference scheme) |
| INIT-IONIC-SHELL | `ionic_local*`: the density list is sorted by |G|^2, so the previous shell's form factor is checked before the HashMap (same cached value) | DONE (always on, bit-identical) | (C) `setup_memo_identity` (verbatim old constructor, 168 partition cases) unchanged. (O) si64 direct path 0.0611 -> 0.0549 s | Small; the direct path's cost is the per-ion per-G sin/cos (12.3M pairs at si64) |
| INIT-IONIC-FAST | `BAANDI_IONIC_FAST=on` (default off, rounding class): structure factor from per-species phase tables exp(-i m_k b_k.tau) built from the full density list (identical on every rank), lazily cached row products T0 T1, one complex product per ion and G; falls back to the direct phase when a unit Miller vector is outside the list or G is not the Miller combination (1e-14 (1 + |G|)) | opt-in (dev08 pending) | (C) `setup_memo_identity::ionic_phase_tables_match_direct_phases_and_partitions`: partitions 1..7 bit-identical, deviation from direct phases 2.5e-17 relative (cap 1e-13) on cubic, skewed, elongated, slab, three mesh conventions, two species; fallback bit-identical on a 1.5 Bohr axis; engine P = 2..7 bit-identical with the switch on. (O) si64 (191735 G, real Si.pz-vbc): 0.0549 -> 0.0134 s (4.1x), max |dV| 2.6e-16 Ha | dev08 reported 0.390 s for this row at one core (Mac 0.055 s), so dev08 sin/cos is the larger share and the gain should be larger there. Same tables would serve `core_density` and the atomic start density (not changed here) |
| INIT-EWALD-INDEPENDENT-REFERENCE | Which Ewald scheme is right at si64 (dev08: E_F - E_R = 1.865e-9 Ha Ewald, 1.794e-9 Ha total)? New test-only evaluator `src/ewald/tests/independent.rs`, sharing no code with either engine scheme or `old_reference.rs`: own lattice algebra (volume in double-double), own image enumeration and reciprocal half space, cutoffs eta r_c = g_c/(2 eta) = 7.5 (analytic tail bounds <= 9e-21 Ha), phases reduced mod 1 in double-double, every term (images, S(G), |S|^2, self, background) accumulated in double-double; same convention as the engine (G = 0 dropped, background -pi Q^2/(2 Omega eta^2)) | analysed; no default change (proposal below) | (C) validation: NaCl and CsCl Madelung constants rel <= 3.8e-16 at eta 0.6/1/2 eta0; one-component lattice with background: sc 2x2x2, fcc conventional, bcc conventional = 8/4/2 x primitive to <= 1.9e-16 rel, Wigner-Seitz constants 1.760118884 / 1.791747230 / 1.791858511 (literature 1.76011888 / 1.79174723 / 1.79185851); forces and stress vs 4th-order finite differences of its own energy 3.1e-13 Ha/Bohr, 1.9e-11 Ha (Omega sigma). Benchmark geometry (fast_screen_dav seed 601, amplitude 1e-3, Z = 4), eta = 1, 2, 4, 8 x eta0: si64 E = -268.7771109045238828 Ha (2 eta0, hi+lo), spread 1.31e-13 Ha; si32 E = -134.3895366548956972, spread 7.9e-14. Default scheme: si64 -268.7771109063182848 (ref - default = +1.794e-9 Ha), si32 -134.3895366549529058 (+5.721e-11). FAST: si64 -268.7771109045238518 (ref - fast = -3.1e-14), si32 -134.3895366548957213 (+2.4e-14). Cause, measured: the default's own term list summed in double-double equals the reference to 2.6e-13 (si64) / 6.8e-14 (si32), and the reference at the default's own (eta, r_cut, g_cut) gives the same +1.794e-9, so it is NOT truncation: the naive running sum of the 3,861,408 real-space terms (accumulator about 1378 Ha, half ulp 1.5e-13) drops 1.795e-9 Ha; 3,469,094 of the terms lie below that half ulp (their sum 3.2e-9 Ha) and are partly lost, one-sided because every term is positive for equal-sign charges. Forces and stress are not affected: default 2.5e-14 Ha/Bohr, 1.3e-13 Ha/Bohr^3; FAST 1.5e-15, 1.4e-16 against the reference | FAST is the correct one (about 5.8e4 times closer at si64: 3.1e-14 vs 1.794e-9 Ha). The default scheme's Ewald energy, hence every default total energy, is low by a geometry-dependent, systematic amount that grows faster than N (si32 5.7e-11, si64 1.8e-9 Ha); comparisons between runs of the same scheme cancel it, but a "tightened reference" or any comparison with another code or with FAST at the 1e-12 level carries it. Proposal (not applied): either (a) make `BAANDI_EWALD_FAST=on` the default (dev08 request init-ewald-fast.md, already 108x faster and accurate), or (b) a switch `BAANDI_EWALD_COMPENSATED=on` that accumulates the reference scheme's real-space terms with Neumaier compensation (cost measured on the Mac: +0.31 ms over a 58.5 ms si64 evaluation, +0.5 %; the rank-split chain then has to pass the compensation word as a fourth double). Tightening the default's cutoffs does not help: the error is summation, not truncation. Test runs in the debug gate (1.9 s): `cargo test --release -p dft-physics --lib ewald::tests::independent -- --nocapture --test-threads=1` |
| SCF-START-DENSITY-WIRING | Is `scf.start_density = "atomic"` a dead option (dev08 master 6be8f09, si64 MPI8: `atomic_start.atomic_orbitals_used = 0`, "supplied initial vectors retained" in both arms)? | FINDING: by design, reporting fixed | (C, Mac) The coarse problem clones the run's options, so the coarse stage (the first stage) builds and starts from the atomic superposition; `atomic_truncation_distance` (si16 1.84e-5) exists only then; si16 coarse maps 10 -> 4. The top-level `atomic_start` block is the FINE session's (start_wavefunctions, not density), whose initial vectors are the prolonged coarse orbitals: "retained" is correct and the coarse atomic start is in `cutoff_ladder.coarse_start_h_vectors` (si16 64). New: `effective_start` (per stage: density, wavefunctions, atomic orbitals used) and `cutoff_ladder.coarse_history` in results.json; engine `CutoffLadderReport.coarse_start_density`, `coarse_atomic_start`, `coarse_history`; test `cutoff_ladder::ladder_reports_the_start_of_each_stage` | Two separate starts (density, wavefunctions) and two stages: report each; a 10-fixed-map screen cannot see a start-density effect (it shortens the coarse stage and the run, not a map) | this session (SCF-numerics owner) |
| START-DENSITY-ATOMIC-LADDER | `start_density = "atomic"` vs uniform on the current stack (ladder on, standard 5e-9, f1e-3, seeds 601-603, T=4, counts) | CANDIDATE for dev08 (supersedes the "no gain" of QE-START, measured before the ladder) | (C, Mac) coarse maps -> 4 on si8/si16/si32 (from 5-10), 7 on si-skew (from 11). Total H (coarse+fine+certificate) 17069 -> 15905 (-6.8 %): si8 -16.0/-27.6/-6.9 %, si16 -14.6/-19.3/-4.1 %, si32 +1.4/+2.0/-2.5 % (verification attempts 1 -> 2 on 601/602), si-skew -20.0/-20.6/-19.8 %; mgo8 (rocksalt 8 atoms, PseudoDojo PBE, Gamma, 60 Ry) -15.9/-16.7 %. |dE| <= 5e-10 Ha. H to the 5e-9 band -11.2 % | The start density pays inside the coarse stage; earlier "no gain" had no ladder. Default unchanged until dev08 to-convergence runs (request scf-trajectory.md) | this session |
| SCF-BOUNCE-ATTRIBUTION | dev08 si64 MPI8 step 3 standard error up 14.8x (0.00180 -> 0.02660 Ry) then down two orders. Reproduce and attribute on the Mac with the new read-only mixer record (`history[].pulay`: differences, kept, Gram condition, coefficients, predicted residual ratio, step length) | FINDING | (C) Reproduces at fine map 3 on si16 (x23, 5.5e-4 -> 1.4e-2 Ha energy error, H 79 -> 119 -> 116), si32 seeds 601-603 (x12/x14/x18), si-skew (x86-x109); coarse map 4 on si16/si32 (x8-x57, m = 2); absent on si8, si-low, si-k222s, mgo8. Every bounce map is evaluated at the density of the FIRST Pulay step (m = 1) after the hand-off: gamma -0.79 to -1.02, step 0.9-1.0 x |R| (plain beta step 0.5); the next step (m = 2, e.g. [-0.73, 0.67], 0.17 x |R|) undoes it. Attribution: (i) inner solve: NO (eigen tolerance 4.85e-7 -> 5e-14 floor in every map: same bounce x19, same gamma -0.988, 2.3x H); (iii) hand-off: NO (ladder off: same at the first Pulay steps of the cold run, map 4 3.1e-2 -> 2.5 Ha); (iv) verification: NO (runs only at the end); (ii) mixer: YES, the L2 least squares sees the screened short-wavelength modes (r2 ~ r1/2, so gamma ~ -1 doubles the step), and the long-wavelength modes that linear mixing at 0.5 already over-corrects are over-extrapolated: their L2 weight is small, their Hartree energy large (stdErr/drho^2 grows 2.5e3 -> 1.5e5 across the bounce). Elongated/slab cells do not converge in either case (MIX-GROWTH-RESET plateau, separate) | BT3's "residual growth is a Pulay trajectory property" confirmed and localised: the first one-difference Pulay step after a (re)start, in the L2 metric | this session |
| MIX-COULOMB-GUARD | `BAANDI_MIX_COULOMB_GUARD=off|trace|on` (new, default off; replicated `PulayMixer`, unweighted fit only; not the M2 distributed mixer, kerker-pulay or `mixing_metric = hartree`): after the L2 fit, measure the predicted residual P = R - dR gamma in the Coulomb metric (`CoulombMetric`, the standard stop's metric); on: if P^T C P > R^T C R refit THIS step's coefficients in the Coulomb metric (history unchanged). No tuned constant: the rule is "the fit must not predict a larger Hartree-weighted residual" | CANDIDATE for dev08 (default off) | (C, Mac, counts, f1e-3, T=4) trace is bit-identical and predicts exactly the bounce steps (ratio 1.46/1.97 si16, 1.23 si-skew; all other steps 0.01-0.93). On, standard 5e-9, seeds 601-605: H to reach the 5e-9 band 25446 -> 22469 (-11.7 %): si16 -19.0/-22.1/-9.8/0/-24.0 %, si32 -9.4/-14.0/-21.6/-12.7/-12.1 %, si-skew -5.5/-5.9/-1.0/-1.0/+2.6 %, si8/si-low/si-k222s/mgo8 0 refits, bit-identical. Full runs 29868 -> 29265 (-2.0 %): si16 -8.6 to -19.8 %, si32 -8.9/+9.1/-1.4/+11.7/-0.3 %, si-skew +5.0 to +6.9 %: the losses are extra verification attempts (si32 s604: stop estimate 3.6e-9 at map 9, true 8.8e-8; si-skew s602: estimate 2.8e-11, true 5.0e-9, and arm A also fails once there). At threshold 1e-10 (si16/si32/si-skew x3): full H 16347 -> 15153 (-7.3 %), to target -11.9 %, |dE| <= 1.5e-12 Ha. With S (atomic start density): to target -16.8 %, full -5.6 % (si32 s601 +19.8 % from 4 verification attempts). Slab-like si8slab: unconverged in both arms. Unit test `coulomb_guard_trace_is_inert_and_on_refits_only_a_predicted_growth` (two-mode sloshing model, plain and incremental update) | The guard fixes the trajectory; the measured full-run gain is diluted by the standard stop estimate under-reading the true error near the stop on some trajectories (a stop-rule/verification issue, VERIFY-FORESIGHT territory, not a mixer one). Not a repeat of the Hartree-metric rows (P40A-1, MIXER-HARTREE-*): those change every step's fit; this refits only steps whose L2 fit predicts Hartree-weighted growth (2-4 per run). Related to ISSUE-21 (response problem) as a cheap detector, not a response model | this session |

## Owned potential and density of the M2 path (M2-OWNED, 2026-11-18; base 6be8f09; Mac counters exact C, wall O; dev08 request docs/performance/requests/m2-owned.md)

Trigger: with GRID+MIX every map still made three N-double allgathervs per rank (effective potential, output density, mixed density). Finding: the distributed H apply multiplies the potential only on the rank's own FFT planes (`multiply_slab` at `grid_start`), and the M2_GRID slab plan produces it on planes of the same rule; the axes coincide unless the wavefunction sphere and the full grid choose different first FFT axes (only seen with a forced grid on an elongated cell).

| ID | Change | Status | Numbers | Lesson / next |
|---|---|---|---|---|
| M2-OWNED-FIELDS | `BAANDI_M2_LOCAL_FIELDS=on` (default off): warm standard maps keep V_eff on the slab planes, moved to the operator layout (`PlaneMap`, alltoall only for crossed axes), `set_owned_potential` with an epoch (checked at the solve: stale = error), kept-start delta V on owned planes (+1 one-double tree sum), complex Hamiltonian full potential gathered only for a consumer (atomic start, screen, complex fallback; M1 fallback errors). Declines counted per reason (cold map, no operator yet, ...) | CANDIDATE (bit-identical; dev08 stage 1) | (C) `tests/m2_owned.rs`: cubic, skewed, elongated with crossed axes, P 2/3: energy, density, eigenvalues, H columns bit-identical to GRID alone; potential gathers per rank 18 -> 1 (cubic), one redistribution per owned map for crossed axes, 0 consumer gathers | The potential gather was the smallest of the three; the value is in the next row |
| M2-OWNED-DENSITY | `BAANDI_M2_LOCAL_DENSITY=on` (with FIELDS, MIX): owned input/output density (epoch-tagged slab planes) through potential, audit, stop rule (epoch match), session metrics (per-plane canonical sums, one allgatherv), plane-mode `DistributedPulay` (`mix_block`, `seed_block`), convex safeguard (exact min); full density only at delivery, ladder hand-off, `materialize_density`; `density()` None at an owned boundary (exact checkpoint rejects) | CANDIDATE (rounding class; dev08 stage 1) | (C) tests: |dE| <= 3.6e-15 Ha vs GRID+MIX and threads-only, same SCF counts, verification passed, ranks equal, reproducible; plane-mode mixer within 1e-11 of PulayMixer at P 1..7, rank-count independent bits. (C) si32 np=8 production+ladder, per rank, arm GRID+MIX -> +FIELDS+DENSITY: allgatherv bytes 100.3 -> 33.5 MB (-66.6 %), calls 231 -> 252 (small reductions), alltoall unchanged, tree sums 470 -> 490; potential/density/mixer gathers 22/22/22 -> 2/0/2 (the remaining ones: cold maps, one-time mixer validation, 2 materializations: hand-off + result), full-buffer copies 42 -> 4; E bits equal, SCF/Dav/H 11/67/1888 both. (O, load 26, 2 pairs) wall 2.270 -> 2.254 s (inside spread 0.14 s); potential phase 0.048 -> 0.022 s, mixing 0.017 -> 0.009 s, density 0.095 -> 0.082 s | Mac cannot price it: the removed gathers cost about 37 ms of 2.27 s here; dev08 at si64 (N doubles per gather larger, allgatherv grows with np) decides. Not done: owned density for spin/finite T/k points/symmetry/NLCC/GGA (declined by the static eligibility, gathered path kept); runtime mixer replacement refused in the owned loop |
| M2-OWNED-DECLINE | Safe by construction before the default flip: consumer inventory of every reachable reader of active-boundary state (only the `checkpoint` subcommand (exact/warm, `Boundary`/`Mixer` reject), an exact restore, `run --agent` mixer replacement (`PrepareFailed`), `heal::run` (panicked on `density()` None) and the hidden `atomic_start_vectors` (stale complex potential overwrote its potential) broke; forces, stress, volumetric/density output, relax predictors, ladder, NSCF read the DELIVERED result, gathered at delivery). New `ScfProblem::decline_owned(OwnedDecline::{Checkpoint, RuntimeControl})`: both owned modes take the gathered path, `Checkpoint` also the replicated Pulay mixer (checkpoint v1 carries it); declared by `checkpoint_io::create` (exact/warm), automatically by `restore_exact` and by every `AgentControlled` session (`--agent`, heal). Reasons counted per map / per session: `owned.fields_declined.<r>` and new `owned.density_declined.<r>` (r adds checkpoint, runtime_control, spin (collinear driver), functional (non-LDA, e.g. PBE)); results.json `scf.switches.BAANDI_M2_LOCAL_{FIELDS,DENSITY}.path` owned/gathered and `fallback_reason` (comma list of reasons). `atomic_start_vectors` clears the stale flag | DONE (C) | (C) `tests/m2_owned_decline.rs`, P=2, all owned switches requested: checkpoint declared -> warm and exact checkpoints of an active boundary succeed, exact restore continues bit for bit, result bit-identical to the gathered replicated-mixer flow, density_declined.checkpoint 4 (2 sessions x 2 ranks), fields 44, owned maps 0; agent-controlled -> boundary density every step, bits = gathered; forces/stress/density output keep the owned mode: FIELDS bit-identical forces and stress, FIELDS+DENSITY |dE| 3.6e-15 Ha, max |dF| 6.0e-14 Ha/Bohr, max |dsigma| 4.2e-17 Ha/Bohr^3, charge 8 - 3e-15; PBE and spin (LDA-PW, 5/3) bit-identical to the gathered path with the reason counted | Forces/stress/density output deliberately NOT declined (they never see owned state; declining would cost relax runs the gain). New consumer of active-boundary state: declare it or extend the automatic declines |
| M2-OWNED-DEFAULT | `BAANDI_M2_LOCAL_FIELDS` and `BAANDI_M2_LOCAL_DENSITY` default ON (switches.rs; `=off` restores the gathered path for A/B and rollback). Tests that compare instrumented or checkpointed runs bit for bit: `instrument.rs` pins the owned switches off (snapshot capture declines them, reason snapshot), `m2_lazy_orbitals.rs` checkpoint test declares `Checkpoint`; no tolerance changed | ADOPTED (dev08 master, np=8, on GRID+MIX: F wall si64 -1.83 %, si32 -2.54 %; FD si64 -4.16 %, si32 -6.03 %; energies bit-identical to gathered, F density SHA identical, FD rounding class; declines only cold atomic start and operator consumers) | (C) `scripts/check.sh` dft-numerics/dft-engine/dft-cli PASS (build, test, fmt, clippy) after the two test pins | Confirm on dev08 with the default binary (no switch env): arm A must equal the earlier FD arm |
| MIX-COULOMB-GUARD-DIST | The SCF-TRAJECTORY guard ported to `DistributedPulay` (chunk mode and the plane mode of the owned density): same rule (after the L2 fit, if P^T C P > R^T C R with P = R - dR gamma, refit THIS step in the Coulomb metric; no constant). C from the owned planes: own `GridPlan` (built only when the switch is not off) with the stop-rule weights V/(2N) 4 pi/G^2, `GridPlan::coulomb_forms` = raw forward transform of signed plane fields, N sum_G w Re(conj(a) b) as per-stick partials, one allgatherv, canonical stick order (rank-count independent). Chunk mode: R planes from the replicated densities, P (and dR on a refit) by one allgatherv each. Per-map record `history[].pulay` now filled on both distributed variants (residual/step norms ride the existing first charge reduction: no extra collective with the guard off) | CANDIDATE for dev08 (default off, bit-identical off) | (C) lib `mix_dist::tests::coulomb_guard_matches_the_replicated_guard_for_every_rank_count`: P 1..7, chunk and plane mode, cubic and skewed cells, odd/even grids, sloshing model with 2 refits: trajectories within 4.2e-11 (off, trace) and 2.9e-11 (on) of `PulayMixer`+`CoulombMetric`, refit flags equal, ratios to 1e-6, trace bit-identical to off, rank-count independent bits. (C) m2_entry_cost production+ladder, Si pz-vbc, guard off vs base 635eb9a: energy bits, density FNV and SCF/Dav/H identical on si16 P2/3 and si32 P2/3/4; trace = off bit for bit. On (M2_MIX + owned density) vs off, SCF/Dav/H: si16 P2/P3 10/48/916 -> 9/50/770 (H -15.9 %), P4 9/52/907 -> 10/61/868 (-4.3 %); si32 P2/P3 11/70/1850 -> 10/53/1649 (-10.9 %), P4 11/67/1888 -> 11/47/1757 (-6.9 %); every count equals the replicated-mixer guard (M2_MIX=off) and, at P2/3, the serial T=1 guard; energies equal to 1-3 ulp; on vs off |dE| <= 9.2e-11 Ha | The MPI arms of docs/performance/requests/scf-trajectory.md no longer need M2_MIX=off; request `docs/performance/requests/mix-coulomb-guard-dist.md` |

## Run-owned configuration, test harness, validated report (RUN-CONFIG / TEST-HARNESS / RUN-REPORT, 2026-11-18; Mac observations O, pure release build, si32-gamma seed 601 amp 1e-3, T=1, load 3 to 16; bits C; base main 1ddcc07; docs/performance/run-config.md, docs/performance/test-targets.md)
| ID | Change | Status | Numbers | Lesson / next |
|---|---|---|---|---|
| RUN-CONFIG | Switch registry resolved once per run into an immutable `RunConfig` (environment > input `[performance]`/`[profile]` > default); session switches resolved per problem; legacy accessors read the thread scope (propagated by `WorkerPool`) or the process configuration installed by dft-cli; MPI agreement hash | DONE (default, bit-identical) | (C) si8 and si32, T=1 and T=4: energy bits, density SHA-256 and iteration counts equal to 1ddcc07. Concurrency test: two runs with opposite `STOP_METRIC` in one process, no lock, each bit-identical to its solo run, independent run counters | 24 switches were cached per process in OnceLocks (first read wins); 6 converted to run-scoped reads here, the rest listed in run-config.md (gamma_real.rs, advance.rs, measurement switches) |
| RUN-CONFIG-OVERHEAD | New binary against 1ddcc07, default path (`fast_screen.py --binary/--binary-b`, 3 and 4 pairs) | measured | (O) median ratio 1.037 +- 8.7 % (load 13 to 16), 1.001 +- 2.9 % (load 4); results bit-identical. Binary SHA-256 at the end of a run about 5 ms (14.8 MB) | Neutral: accessors add one thread-local read; nothing is looked up in a hot loop |
| PROFILE-OVERHEAD | Uninstrumented arm A against `BAANDI_PHASE_PROFILE=on` (same binary, 4 pairs); `BAANDI_DAVIDSON_STAGES=on` (4 pairs); stage 0 `fast_screen_dav.py --b BAANDI_PHASE_PROFILE=on` (3 pairs) | measured | (O) phase profile 0.994 +- 2.8 % (bit-identical); Davidson stage timers 1.012 +- 1.0 %; stage 0 cold solve 1.020 +- 7.3 %, warm 1.003 +- 5.2 % (same Davidson iterations and H columns) | Phase profile cost is below the noise; Davidson stage timers cost about 1 % (keep both arms on when they are on, as before) |
| TEST-HARNESS | 112 of 172 integration-test files with no process-global state moved into one `tests/all` binary per crate (7 harnesses); 60 stay separate with their boundary listed | DONE | (O) test executables 184 -> 81; build + link after touching dft-numerics 10.1 to 11.0 s -> 6.0 to 6.4 s; test-suite wall (4 jobs) 19.4 to 23.1 s -> 13.6 to 14.2 s; full check.sh PASS on both trees | The remaining binaries are dominated by forced switches and process-global counters; they join once their overrides become per-call RunConfig and their counters run-local |
| RUN-REPORT | results.json `run` block (run ID, exit reason, binary SHA-256, agreement hash, run counters, completion marker), switch manifest with requested/consulted/path/fallback, phase profile closed per rank (`phase-profile-2`), `scripts/parse_results.py` | DONE | (C) parser rejects: truncated and interrupted files, stale binary SHA, profile switch on without profile, buckets with maps and zero wall, zero covered time, Davidson stage timers all zero, collectives with calls and zero seconds, min > max, negative other, covered above a rank's wall; simulated rank imbalance: legacy sum of maxima 5.0 s above both walls (4.5, 4.0 s), per-rank closure 4.0 / 3.5 s, residual 0.5 s | dev08 jobs should run the parser on every results.json with the job's binary SHA before using a number |
| SCF-TAIL-ACCOUNTING | Where the last maps and verification attempts go (master question after dev08 si64 np=8: maps 6-9 2.3e-7 -> 7.4e-8, a failed verification at map 10 with true error 7.7e-10). Mac release T=4, f1e-3, production defaults (ladder, guard, VERIFY_REUSE=half), 10 systems (si8, si16, si32, si-skew, si-low, si8-k222, mgo8 PBE two species, syn-gaussian-fcc, syn-p3, syn-spd) x seeds 601-605 x thresholds 5e-9 / 1e-10 = 100 runs, true error per map against a strict reference run of the same geometry (`BAANDI_REFERENCE_STRICT=on`) | measured, finding | (C) all 100 converged and verified, max delivered true error below class in 100/100. Failed verification attempts: 64, every one of them verified a state whose true error was ALREADY below the class threshold (0 of 64 above). Maps after the true error entered (and stayed in) the class, medians over seeds: si16 1-2, si-skew 2, si8 1-2, si8-k222 1-2, syn-gaussian-fcc 1-2, si32 0 (5e-9) / 2 (1e-10), mgo8 0-1, si-low 0-1; their H plus failed attempts after that point = 14.5 % of all H columns (15771/108866). Verification maps = 14.4 % of all H: failed attempts 6.7 % (median 1.33 x the H of their stop map, 0.53-2.80), passed 1.02 x. Estimator / true error at the attempts (T >= 1e-11 Ha): stop estimate S min 0.33 (syn-gaussian-fcc at T~1e-12) / 0.73 (si32), median 1.96; verification V = E_H[F(rho_out) - rho_out] min 0.19 (syn-spd), median 14, max 298 (si-skew); V/T medians si-skew 83, si16 37, si32 21, si8 13, si8-k222 5.5, mgo8 2.5, si-low 2.5, syn-p3 0.74, syn-spd 0.19 | The stop estimate does track the true error (S/T about 1-3); the verification criterion reads high on screened cells and LOW on weakly screened ones, exactly as the linear Hartree model says (V/T = A (1 + A) per mode, S/T = 1 + 1/A). The pair rule (S < thr AND V < thr) never under-reads: max(S, V)/T min 1.48. Failed attempts are therefore false failures of a conservative pair, not under-reading stops | tests/reference: none (Mac scratch runs; scripts in the report) |
| SCF-TAIL-PREDICTOR | `BAANDI_STOP_PREDICTOR=off|trace|on` (new, default off; `session::verification_cross`, `StopPredictor`, `CROSS_SAFETY` = 4): cross estimate -X = -E_H-form <rho_out - rho_in, F(rho_out) - rho_out> by polarization from one extra Coulomb-metric evaluation of E_H[F(rho_out) - rho_in] in the verification map (no solve). Linear Hartree model: -X / T = 1 + A >= 1 (V reads A (1 + A)). trace prints S, V, -X and the verification energy change; on also passes an attempt whose old criterion fails when 0 < -X < thr / 4, S < thr and the residual guard passes; replicated densities only (owned-plane loop: old criterion) | trace kept as an instrument; `on` NOT supported by the evidence (default stays off); dev08 trace on si32/si64 requested to close it | (C) trace = off bit for bit in 100/100 runs (energy bits, SCF count). -X/T at the attempts (T >= 1e-11): calibration set (si8, si16, si-skew, mgo8, syn-gaussian-fcc) n=58 min 0.654 (syn-gaussian-fcc s603 5e-9: T 1.06e-9, -X 6.9e-10) p05 0.66 median 5.0; held-out (si32, si8-k222, si-low, syn-p3, syn-spd) n=30 min 1.57 median 3.1; X itself is negative in 164/164 attempts as the model predicts. One-mode fit T ~ S(-X)/(S - X) under-reads (min 0.22): rejected. Offline replay (trajectory after the stop removed): pass rule at thr/f: f=1 34 runs earlier, H -8.2 %, worst delivered T 0.37 thr; f=2 19 runs, -1.6 %; f=4 10 runs, -0.36 %; f=8 2 runs; f=16 none. Measured `on` arm (f=4; si8, si16, si32, si-skew, si8-k222, syn-gaussian-fcc x 601-605 x 5e-9/1e-10 = 60 runs): exactly the replayed 10 runs change (si-skew 5, si16 2, syn-gaussian-fcc 3; each one map and one attempt less, -9 to -14 % H), total H 97982 -> 97590 (-0.40 %), all verified, worst delivered true error 0.275 x threshold; si32, si8, si8-k222 unchanged. Model: -X/T = 1 + A grows with the cell (si-skew median 9.4, si16 6.5), so at dev08 si64 map 10 (V 1.0e-8, T 7.7e-10) -X is expected near 3e-9, above thr/4 | The verification failures are structural (V/T = A(1+A)); -X removes one factor of A but not the second, and one measured point (0.654) forbids using -X without a safety factor that eats the gain. No estimate built from these three maps is both tighter and conservative on the Mac set. A real gain would need a different certified quantity (for example delivering F(rho_out) after a failed attempt, an ADR 0043 V4 change and a user decision) | docs/performance/requests/scf-tail.md; engine tests `verification_cross_matches_the_linear_model`, `stop_predictor_trace_is_inert_and_on_keeps_the_class` |
| SCF-TAIL-STALL | Attribution of slow tail maps (Mac set of SCF-TAIL-ACCOUNTING, maps >= 5, 487 maps) | finding, no change tested | (C) stall maps (density contraction > 0.8, energy error flat or up) 47 of 487 = 6.4 % of their H; only on silicon: si8 32 %, si32 13 %, si16 12 %, none on mgo8, si-skew, si-low, si8-k222, synthetic cells. Pulay window state at the stall: sliding (window full at 6, oldest pair dropped) on si32 10/38 vs 0/40 while growing, but si8 14/40 growing vs 11/39 sliding and si16 6/40 vs 2/29; median contraction sliding 0.44 vs growing 0.31; Gram condition up to 5e5 at the stalls. Inner eigen solve at the stall maps is TIGHTER than elsewhere: max occupied residual / density residual median 1.43 (stall) vs 6.69 (other maps), max residual grew against the previous map in 6 % vs 23 %, tolerance unchanged from the previous map 30 % vs 18 % | Not attributable to one general cause: the inner schedule is not it, and the window effect is a silicon-only signal (si32). Under the generality rule no mixing-window change was tested |
| PLATFORM-TESTS | Platform-aware test baselines: preserve aarch64 assertions; distinguish rounding trajectories from exact same-platform work/rank/thread identities; correct stale MKL pair-plan default assumption | Mac fast gate PASS; dev08 confirmation pending | Mac fast gate: 64 binaries, 583 passed, 0 failed, 157 existing ignored; build/test/fmt/clippy exit 0, 13 s wall (observation only); MKL-feature test clippy/type-check exit 0 on Mac (no linking or MKL execution). First run exit 1: three unrelated pseudo-cache tests rejected the repository-local TMPDIR; rerun with /private/tmp scratch passed. Relayed e4186f9 failures: spin energy bit patterns differ by 32 ULP = 3.552713678800501e-15 Ha (derived from the supplied integers, not the ticket's 6-ULP label); reuse SCF 29 vs 26; clustered Davidson 39 vs 35; MKL one-rank Gram differs by 1-2 ULP; selected residual difference about 2e-16 Ha. No speed measurement | MKL PAIR_FFT defaults to pruned, independent of shape. Toy energy 1e-12 Ha, standard energy 5e-9 Ha, roots 1e-11 Ha; residual uses 1e-12 relative to a unit-Ha floor on non-aarch64 (pure relative 1e-12 would reject the supplied values). No production changes or new ignored tests | docs/performance/slow-tests.md, platform-dependent baselines; tickets/PLATFORM-TESTS.md |

## F-LAUNCH-IO (2026-11-19; base ff0c83caaccdfc26b83598cc7b073d980f055d17; Mac observations, dev08 pending)

| ID | Mechanism / decision | Evidence and limits |
|---|---|---|
| F-LAUNCH-IO-CLOCKS | Add `launch-1` to completed closed-shell `run` results only when the resolved phase-profile switch is on (environment or input). Main-entry bootstrap, MPI init/finalize, argument/switch/input parsing, preflight, prepared-value reuse, UPF read/parse/hash, construction, materialization/output writers, phase gather/report, run-block hash, JSON serialization/write, Rust teardown and console write are measured. No engine phase clocks or numerical code changed. | Each rank retains its existing files and appends its own small launch object after dropping the engine communicator and finalizing MPI. No extra collective, barrier, full-result rewrite, fsync or process::exit. The original JSON bytes before the final brace are preserved. `closure.before_init_s + after_profile_s + within_profile_residual_s = outside_covered_s = main_wall_s - covered_s`, using THIS rank's covered wall. Nested detail clocks are not additive. `main_wall_s` excludes the final launch append/diagnostic and post-main runtime/OS teardown; external B also includes the launcher/loader before main. No claim that a local main timestamp measures those intervals. Four bootstrap Instant reads remain with profiling off because the input profile table is not known at main entry; subsequent launch hooks take no clocks or allocations when off. Library callers and spin/FD drivers do not produce the process-lifecycle section. |
| F-LAUNCH-IO-SPLIT | Small real Si diamond, two atoms, 30 Ry, default standard path with verification/ladder, pure release; 1 rank and T=1/2/4, three runs each. mpirun is unavailable. Initial nine-run instrumented medians: binary SHA **0.034410167 s**, preflight including prepare/first UPF **0.000860542 s**, second UPF read/parse/hash **0.000500917 s**, construction **0.002567958 s**, materialization including outputs **0.001192750 s**, profile report **0.000416166 s**, JSON serialize **0.000254958 s**, JSON write **0.000164375 s**, density.bin write **0.000171083 s**, teardown **0.000099875 s**. | All except binary hashing are below 10 ms in this small case. MPI init/finalize use the feature-off single-rank stand-in, NOT a real MPI measurement. Original/instrumented external wall medians (s, T=1/2/4): 0.081796458 / 0.078557083 / 0.081632833 and 0.077764334 / 0.074299875 / 0.072510334; different loaded intervals, NOT evidence of a profiling speedup. One-minute load 27–35 during these batches. Raw outputs and logs: worktree `target/f-launch-io/{baseline,instrumented}.json`, corresponding output directories. |
| F-LAUNCH-IO-HASH | Candidate: stream the executable through the existing SHA-256 implementation with 64 KiB scratch, retaining the exact digest/path, error behavior and once-per-process identity. Removes the executable-sized allocation (baseline binary 15,103,888 bytes). Keeps synchronous hashing and every provenance field. | Initial nine-run streaming median **0.029459042 s** versus **0.034410167 s**; longer alternating AB/BA screen (31 pairs, discard first = 30 measured pairs) did NOT establish a speed gain: hash medians whole-file **0.042655750 s**, stream **0.042205125 s**; external walls **0.1095370205 / 0.1082419790 s**. Median paired hash saving **-0.0003696045 s**, range **[-0.021707832, +0.009223291] s**, only 14/30 positive; load 15.25/20.31/17.63. Keep as a bounded-memory candidate for dev08 review, not an adopted speedup. Residual >10 ms is mostly the required exact SHA scan; the vendored sha2 library uses its software backend on this Mac unless an additional feature/dependency is enabled (manifests are outside ownership). Lazy/omitted/cached hashes would change provenance or risk stale identity; background hashing would consume an extra runnable thread and needs an equal-core A/B before adoption. No such change made. Raw pairs: `target/f-launch-io/hash-pairs.json`. |
| F-LAUNCH-IO-PREFLIGHT | Verified `preflight::check` calls Config::prepare once and returns it; `run` takes that value (new `prepare_reused` proof field). It does NOT prepare again on this successful path. | Source review: preflight does structure review once, reads/parses each UPF, builds prepared bases, and builds a per-k basis again for subspace validation; run reads/parses the UPF again. Thus some work IS redundant, although it is already inside `init.cli`, not wholly in B minus covered. The measured inclusive preflight is <1 ms here; separate internal prepare/UPF clocks and elimination of this duplication require `preflight.rs`/`config.rs`, outside this ticket's owned files. Do not re-run preparation just to time it. No radial integration/interpolation table build occurs in the CLI before ScfProblem construction (preflight only inspects PP_PSWFC counts); engine construction is already detailed by existing init/atomic profile keys. No replicated numerical tables changed without a measured dev08 cost. |
| F-LAUNCH-IO-OUTPUT | Keep existing per-rank evidence files, one whole-buffer write per JSON/density file, and normal destruction/finalization. No writer rewrite beyond the small late launch append. | Small-case writes and teardown are sub-ms. Suppressing rank directories would change the existing output contract; a new gather solely for launch timing would itself add latency. Finalize cannot be included in a pre-finalize gather. Rank-local results allow the master to compare min/max after the run, without collectives after MPI_Finalize. Large-system end-of-run work and real MPI startup/finalize remain dev08 measurements. |
| F-LAUNCH-IO-VALIDATION | Initial gate used `TMPDIR_990=$PWD/target/tmp`: build/fmt passed, 448 tests passed, 3 cache tests failed because the cache is intentionally forbidden inside a repository. Tests were not weakened; rerun uses the sandbox-writable external temporary directory. | `ps` load telemetry and `nice` are sandbox-denied; neither is a correctness gate. Initial clippy also caught a byte-array style warning in the new append helper; corrected to a byte string. The first supplemental identity harness had missing explicit k-point weights and was corrected after the baseline rejected that input (not an engine regression). README test suite has an unrelated existing stale-phrase assertion (`standard profile is unsupported`); the focused README block and figure freshness checks both pass. README/evidence schemas, session/dashboard artifacts, and generator are outside ticket ownership and unchanged. |
| F-LAUNCH-IO-GATES | `TMPDIR_990=/private/tmp/f-launch-io-check CHECK_NICE=0 scripts/check.sh fast dft-cli dft-engine dft-physics` **exit 0**, build/test/fmt/clippy PASS: **451 passed, 0 failed, 127 ignored**, 47 binaries, 19 s gate wall (loaded Mac). Four new unit tests cover late append byte preservation/truncation, inactive clock behavior, local closure/state cleanup and streamed SHA chunk boundaries/errors. | Supplemental baseline-versus-candidate process checks: Gamma, skewed, elongated and two explicit k points; each with profile off, environment on and input-table on (**12 combinations**). Exact energy bits, density.bin and density.json bytes, work/history/scientific JSON agree; only existing nondeterministic timing/timestamp/run-ID/binary-identity fields and the new launch object are excluded. Actual executable hashes independently checked by Python hashlib. `prepare_reused=true`; launch absent off; local launch closure sum agrees within 1e-12 s. `target/f-launch-io/identity.json`, `identity-check.log`. |

F-LAUNCH-IO reproduction and handoff (no QE run or cross-code speed claim):

- Baseline executable SHA-256: `f083c539433d362aec9293e73f12d5c63926375f141335bdc3cfbec321482c0d`; instrument-only `bcb52af858b0e2b4f81c8e746aa44d48428ab6fb40f16245ad6b96f7bf92323c`; streaming screen `4075fa0d8a781a8d1f5659dc39946d9b2bbd30a97e516637db6fe9888f1bcdc4`. These are saved binaries of successive changes, not commits. Final build contains later closure/prepare-proof metadata; its identity is recorded in the final measurement JSON.
- External Si.pz-vbc.UPF SHA-256: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`. No UPF added to the repository. Si2 energy **-7.2955304681665085 Ha**, delta from baseline **0 Ha**, density SHA-256 **32c29dcbd8a8433311caf83c30c6ed3250df883786d2e5eedb5d362da65d4967**, identical for every T=1/2/4 sample and all 62 hash-screen executions.
- Local: `cargo build --release -p dft-cli`; inputs derive from `examples/si_lda.toml`, replacing the external UPF path and appending `[scf] threads=1`, 2 or 4. Scripts and raw evidence are retained under ignored `target/f-launch-io/`: `measure.py`, `hash_pairs.py`, `identity_check.py`. Timings are observations only; no CPU affinity or quiet-host acceptance on this Mac.
- dev08, master only: build the base and reviewed candidate in separate target directories with the SAME features/toolchain as the original si64 np=8 measurement. Record both binary hashes and preserve the input/UPF, switches, binding and start. Example candidate command (substitute the master's existing paths):
  ```sh
  BAANDI_PHASE_PROFILE=on /usr/bin/time -p mpirun -np 8 --map-by core --bind-to core "$BIN" run "$SI64_INPUT" --output-dir "$OUT"
  ```
  Run alternating baseline/candidate pairs at np=1/2/4/8 (same parallel model in each pair), check byte-identical density, exact energies/work and independent binary hashes before comparing external B. Read `$OUT/results.json` and `$OUT/rank*/results.json`: `launch.seconds.{mpi_init,mpi_finalize,binary_sha256,preflight_including_prepare_and_upf,results_serialize,results_write,run_teardown,profile_gather}`, milestones and closure. Compare each rank's local outside-covered time; report the range over ranks. Do not sum maxima from different ranks. Also retain B minus `phase_profile.closure.covered_max_s` and the legacy B minus `covered_wall_s` for comparison with the original 0.5386 s budget; the latter can combine maxima from different ranks. External launcher/loader/post-main work remains an external-clock residual, not an invented in-process measurement.
- Pending: attribution of the actual dev08 0.5386 s, real MPI feature build/lifecycle check, si64 costs of replicated preflight and output, and quiet-node adoption evidence for streaming SHA. No optimization default is claimed accepted from Mac timings. README block/figures remain current; the separate stale-phrase Python test needs its owner to reconcile it with the already-merged standard-force documentation. Dashboard/session updates belong to the master (registry confirms this delegated run is `running`; worktree artifacts were not rewritten).

F-LAUNCH-IO final release observation: `cargo build --release -p dft-cli` exit 0;
SHA-256 `0c257198ef83e4777a91731f91ecce918e15b02a5368eaf644a6cf34c6d8f048`.
`python3 target/f-launch-io/measure.py final target/release/dft-cli 3` produced
nine completed runs (`final.json`), then the twelve identity combinations were
repeated against that final binary (`identity.json`). `scripts/parse_results.py`
accepts `final-t1-0/results.json` with the independent binary check (exit 0).
Final T=1/2/4 medians, seconds: external wall **0.086323250 / 0.076963375 /
0.075183041**, main wall **0.080064834 / 0.072605500 / 0.070469708**, local
outside-covered **0.038107583 / 0.036962535 / 0.036992086**. The nine-run median
hash is **0.035835917 s**, preflight **0.000983083 s**, UPF **0.000489708 s**,
JSON serialize/write **0.000248458 / 0.000239750 s**, teardown **0.000107625 s**.
Load at completion: **14.76 / 37.76 / 29.82**; these observations do not establish
a speedup over the earlier batches. The remaining ~36 ms is dominated by exact
hashing on the Mac; no inference about dev08 MPI startup/finalize is made.
| PLATFORM-TESTS | Platform-aware test baselines: preserve aarch64 assertions; distinguish rounding trajectories from exact same-platform work/rank/thread identities; correct stale MKL pair-plan default assumption | Mac fast gate PASS; dev08 confirmation pending | Mac fast gate: 64 binaries, 583 passed, 0 failed, 157 existing ignored; build/test/fmt/clippy exit 0, 13 s wall (observation only); MKL-feature test clippy/type-check exit 0 on Mac (no linking or MKL execution). First run exit 1: three unrelated pseudo-cache tests rejected the repository-local TMPDIR; rerun with /private/tmp scratch passed. Relayed e4186f9 failures: spin energy bit patterns differ by 32 ULP = 3.552713678800501e-15 Ha (derived from the supplied integers, not the ticket's 6-ULP label); reuse SCF 29 vs 26; clustered Davidson 39 vs 35; MKL one-rank Gram differs by 1-2 ULP; selected residual difference about 2e-16 Ha. No speed measurement | MKL PAIR_FFT defaults to pruned, independent of shape. Toy energy 1e-12 Ha, standard energy 5e-9 Ha, roots 1e-11 Ha; residual uses 1e-12 relative to a unit-Ha floor on non-aarch64 (pure relative 1e-12 would reject the supplied values). No production changes or new ignored tests. Follow-up: relayed 7a8ffe2 dev08 gate passed 6/7 failures; degenerate4 residual delta 1.386e-11 Ha exposed per-band basis dependence. Compare cluster residual squared norms at 1e-12 with unit-scale floor; cluster gaps < 1e-8 times common spectrum scale, partitions agree between paths. Original aarch64 bounds retained; independent true residual checked FIRST for both arms. Mac detected singleton scattered/clustered, quartets [0..4,4..8], triplets [0..3,3..6] plus singleton [6..7] at capacities 3/6 and P=2/3; targeted tests 4 passed, 1 existing ignored; all 140 solver unit/integration tests passed, 25 existing ignored; clippy/fmt exit 0. Full cargo test exit 1 at rustdoc scratch-directory creation. Requested fast gate blocked before execution by sandbox denial of mktemp in the required external-volume scratch directory | docs/performance/slow-tests.md, platform-dependent baselines; tickets/PLATFORM-TESTS.md |

## F-LADDER-HANDOFF — closed negative experiment (2026-11-19)

- Source `7243075`, Mac pure Rust, one core, no new QE run. Current solved
  coarse orbitals are prolonged to the unchanged final basis/grid; one occupied
  fine H block supplies the kinetic-inverse high-G indicator (shifts 0.5/1/2
  Ha, policy shift 1 Ha). This is an empirical indicator, not an LDA/PBE bound.
  The earlier sandbox blocker was resolved by the user's worktree-local
  `TMPDIR_990=$PWD/target/scratch` correction; no wider permissions were used.
- **Scientific status CLOSED / policy NO-GO.** 126 fixed-stop candidates on
  si8, si16, si32, si-skew, MgO8, synthetic skew and synthetic slab, seeds
  601–603; all converged and verified. Six stops from 1e-2 to 3e-5, fraction
  0.5. Counts include atomic start, handoff, all fine/redo solves, every
  verification attempt, and delta-potential updates. H columns and block calls
  are separate; every requested probe is charged (one block, 2–64 fine columns).
- Registered sparse-trigger energy-balance rule, q in {0.001, 0.01, 0.1, 1},
  trained on seed 601: best q=0.1 already has mean work ratio 1.004606.
  Held-out seeds 602/603, median paired total-work proxy / adopted rule:
  MgO **0.906958**, Si8 **0.980086**, syn-slab **0.853172**; si16 **1.074884**,
  si32 **1.074247**, syn-skew2 **1.000143**, si-skew **1.000000**. Only **3/7**
  wins and no skewed win, below the four-system gate. The global fixed-stop
  training optimum is 1e-2 (ratio 0.904042), also not uniformly better on the
  held-out systems. No indicator switch or production default was added.
- Calibration: high-G variation over unique states is only **0.58–5.08 %**
  within a system, while remaining fine/verification work varies substantially;
  its association with recovery work is not transferable. All 21 coarse
  trajectory prefixes agree across stops. Analytic weighted-residual test,
  real-versus-complex H probes (1e-11 Ha entry tolerance), 21 exact ordinary
  ladder controls and 21 capture-on/off energy/work comparisons pass. An
  independent reader validates **2,202** snapshot manifests and binary shapes.
- Additional **105** capture-disabled configuration runs, all converged and
  verified. Fraction 0.35 reduces the work proxy on four systems, but every
  one needs more fine work; 0.4 improves three. Minimum valid coarse grids
  improve the proxy on five, but observed Mac walls regress on MgO (1.361),
  syn-skew2 (1.144) and syn-slab (1.063): grid-volume weighting misses FFT
  factorization costs. A coarse-only 4x looser ratchet helps MgO/si16 only;
  its factor=1 driver matches density/orbitals/H exactly on all 21 geometries.
  No configuration default changes; no Mac wall observation is adoption evidence.
- Accuracy against accurately completed strict fine references of the same
  geometries: max |dE| **1.30466659698e-9 Ha** (handoff scan),
  **4.31099822507e-9 Ha** (configuration scan), both below the unchanged
  5e-9 Ha standard class. All fine cutoffs, grids, PPs and verification unchanged.
- **Required gate not green:** `scripts/check.sh fast dft-engine dft-numerics
  dft-cli` exits **1**. Build/fmt/clippy pass; **446 passed, 3 failed,
  133 ignored**, 48 binaries. Three unchanged CLI pseudo-cache tests require
  caches outside repositories and reject the mandated in-worktree TMPDIR.
  The experiment's initial duplicate-fixture-module lint was fixed; the final
  clippy gate passes. Scope isolation `scripts/check.sh fast dft-engine
  dft-numerics` exits **0** (315 passed, 0 failed, 93 ignored; all four gates
  pass); it does not waive the required gate. Cache-policy/test changes are
  outside this ticket's scope.
- [Complete design, per-candidate counts, coarse phase/tolerance tables,
  commands, hashes and limitations](../design/ladder-handoff.md).
  Reproducer: `scripts/ladder_handoff.py`; raw JSON/CSV and snapshots:
  `target/scratch/handoff-{main,config,analysis}`. Preserve/export those raw
  artifacts before worktree cleanup. The implementation is test-only and the
  long experiment is ignored by normal gates. **Lesson:** an almost constant
  high-G discretization floor does not select a cheaper handoff reliably;
  reduced coarse work must earn back its fine recovery and probe costs.

## F-INIT-STATES: owned projector initialization (2026-11-19)

| Attempt | Mechanism | Status | Evidence and limits |
| --- | --- | --- | --- |
| F-INIT-STATES | `BAANDI_INIT_OWNED_STATES=on` (default off): owned Gamma complex projector representatives, lazy full table, fixed-chunk atomic overlaps; root-only initial ionic FFT | Implemented; aggregate gate blocked by cache-policy tests; dev08 adoption pending | Base `7243075d009e51eed0e8b59786f915ef0048dcb1`. Tests and local observations below; no solver changes. |

The ticket explicitly authorizes revisiting INIT-NONLOCAL-TABLE / INIT-PACKED
as a rounding-class change (cap 1e-12 Ha), retaining the replicated default.
This changes the premise of the earlier deferral. The local class and count
checks are recorded below; default adoption still requires dev08 evidence.

### Full complex table consumers at this revision

* `hamiltonian.rs`: `NonlocalOperator::new_with_rule` creates every full-sphere
  column, shared through one `Arc<Vec<Complex64>>` by the operator and atom
  views. Scalar `add` / `apply` / `apply_into`, `expectation_values`, and the
  block-panel application all read it. `Hamiltonian` complex applications
  therefore require it, including complex Davidson and fallback paths.
* `distributed_real_gamma.rs`: `DistributedRealHamiltonian::new` calls
  `source.nonlocal().atoms()` and `pack_projector_rows` for each column. Its
  retained real table is already owned-row only, but the packing validation
  reads the full complex column. Replacing just construction without replacing
  that validation would still materialize the table on every rank.
* `real_gamma.rs`: `RealNonlocal::new` reads every column for the serial real
  operator. `gamma_real.rs` can construct that operator lazily for the M1 or
  serial fallback; the lazy full complex table must remain available there.
* `scf.rs::atomic_guesses` calls `atomic_start::ritz`; `split_gram` distributes
  **columns**, applies complex H over **all G rows**, then gathers the small
  Gram matrix. The unsplit branch also applies complex H. Supplying `hq` to
  the existing Ritz API avoids both complex applications, but still computes
  the Gram from that image. `gamma_real::atomic_apply_real` explicitly declines
  an active communicator, so enabling ATOMIC_REAL alone cannot fix MPI init.
  `atomic_start.rs` is outside this ticket's owned files; this implementation
  uses the existing interface through scoped projector weights.
* `scf.rs::explicit_nonlocal_energy` and `energy.rs`'s scalar/projector
  expectation audits require full columns. The M2 ordinary audit can instead
  use `m2_nonlocal_energy` / `nonlocal_band_energies`, or the band-sum identity;
  passing a `&NonlocalOperator` into an audit must not itself trigger a build.
* `scf.rs::operator_fingerprint` traverses all columns for snapshot capture
  and replay. It is an explicit full-table consumer outside plain SCF.
  Orbital materialization itself is a gather/unpack; subsequent explicit
  audits or fingerprints are what require the complex projector data.
* The certificate pre-screen reads the complex **local potential**, densities
  and delta-V integrals, not projector columns directly. A resulting real M2
  solve can remain owned; a complex fallback or explicit final energy audit
  cannot. Potential materialization and projector materialization need
  separate counters and eligibility decisions.
* Post-SCF forces construct a separate `NonlocalOperator::new` in engine
  `forces.rs`; physics `forces::nonlocal` consumes full columns. Physics
  `stress::nonlocal` likewise constructs its own operator and contracts full
  columns with strain derivatives. Keep these constructors eager by default;
  their correctness must not depend on SCF's owned representation.
* Collinear workspaces call `scratch_len` and `apply_into` on their operators.
  Non-Gamma, spin, incompatible real starts and other excluded configurations
  must retain the ordinary constructor or an explicit lazy full fallback.

### Implementation and arithmetic

`NonlocalOperator` retains only each rank's positive Gamma representatives,
selected with the same `position_owners` rule as the distributed operator.
Species radial/angular factors (full G rows, shared across that species' atoms)
and positions form an owned reconstruction recipe. This is O(N_G times species
channels), not the replicated O(N_G times atom channels) table. Full phase
arguments and factors are checked at construction; the distributed operator
packs directly from the owned representatives, without calling `columns()`.
A shared `OnceLock` materializes full complex columns only for a consumer that
reads them. Both radial quadrature rules reproduce the eager columns bitwise.

At atomic start, P = B^H Q is formed in globally numbered 1024-G tiles. Batches
of eight tile slots use one tree reduction; each slot has one contributing
rank, then every rank adds the slots in fixed global order. Thus the projection
association is independent of rank count (unlike a sum of rank-local partial
sums). D P is small and replicated. A scoped Hamiltonian cache supplies each
rank's column subset of these weights to the existing `split_gram`: local and
kinetic FFTs and Q^H H Q stay column-split, and B D P reconstructs only a
bounded row tile at a time. The usual small Gram gather, Hermitian cleanup,
eigendecomposition and rotation remain unchanged. The cache clears on success,
error and unwind; subsequent complex applications lazily build the full table.
This computes the nonlocal reduced term P^H D P without gathering a full H Q
or adding a replicated full Gram product. An initial H-Q-gather implementation
was replaced during review before retained performance observations.

Counters are exposed as `owned.states_full_builds`, `owned.states_table_bytes`,
`owned.states_replicated_bytes`, and `owned.states_recipe_bytes` in run counters
and the existing M2 counters (`m2.owned.*` in fast_screen JSON). Full-build counts
exclude eager construction. Byte counters accumulate allocations across operator
generations (including the coarse ladder); they are **not peak RSS**. Full
fallback storage increments table bytes. They exclude the distributed real
table, temporary atomic tiles/weights and basis metadata. A zero full-build
counter in both A and B alone is not evidence: tests also assert sparse retained
bytes, execute the owned path and force a counted full fallback.

### Other init replication (source inspection, no new timing)

| Bucket | Per-rank work | Remaining constraint |
| --- | --- | --- |
| `cli` | Input/UPF reading, preflight and prepared bases/grid | INIT-CLI already removed the second prepare call. Production CLI edits are outside this ticket (CLI test ownership is switch-count only). |
| `basis` | Gamma basis, density basis, two FFT grid specifications | The engine already reuses Gamma and the density sphere under INIT-REUSE. Immutable sharing inside a process could help simulated ranks, but is not shared memory across MPI processes. Basis enumeration changes are outside ownership. |
| `densities` | FFT plan, optional core and atomic-superposition density, ionic inverse FFT and real extraction | With the switch on, rank 0 alone performs the initial ionic inverse FFT and broadcasts its real values. Other ranks avoid the complex output grid allocation and transform; failure is broadcast before the field. Tests compare the field bits. Mutable plans remain per rank. |
| `rest` | Validation, worker FFT plans, occupations, Coulomb metric and force-model identity | Coulomb metric already reuses the density sphere. Removing validation or sharing mutable FFT workspaces is not an equivalent optimization. No additional removal was implemented or measured. |

### Environment and verification

The first preflight could not write the external scratch directory (exit 1).
The user's correction resolved that restriction: all subsequent builds, tests
and experiment outputs use `TMPDIR_990="$PWD/target/scratch"` and
`TMPDIR="$TMPDIR_990"`, with `/opt/homebrew/opt/rustup/bin` prepended to PATH.
No `/tmp` fallback or wider permission was used. The local registry has no
F-INIT-STATES entry; registry/dashboard/README changes are outside the ticket's
explicit owned list. Changes are uncommitted for master review.

Required gate (executed; outcome below):

```sh
scripts/check.sh fast dft-physics dft-engine dft-numerics dft-cli
```

### Verification results (C)

`cargo test --release -p dft-engine --test init_owned_states -- --nocapture`
returned **0**: **6 passed**, 0 failed/ignored, 0.33 s test execution. Evidence
log: `target/scratch/test-acceptance.log` (scratch, not committed).

* Cubic and skewed cells, simulated ranks 1..7: fixed-chunk D P weights are
  bit-identical across rank counts for the same Q; maximum H Q entry difference
  **2.7755575615628914e-17 Ha** against the full-table complex operator. The
  fixture crosses the 1024-G tile boundary at appreciable projector values.
* Initial atomic Ritz vectors are independently evaluated with a full-table
  Hamiltonian at the uniform input density (zero Hartree G=0, constant PZ XC,
  no NLCC). Maximum energy difference **9.7699626167013776e-15 Ha** across
  cubic/skewed cells and ranks 1..7; no tolerance relaxed.
* Complete production SCFs over the same shapes/ranks: maximum energy
  difference **6.2172489379008766e-15 Ha**, **12 SCF steps in both arms**,
  identical atomic H-application counts, passed fixed-point verification.
  Default owned fields/density are enabled. M1 fallback is separately forced
  with gathered fields: full table lazily built exactly once per rank, versus
  zero on the ordinary owned M2 path.
* One versus two intra-rank threads: matching energy bits, density and step
  count, sparse construction in both cases. MPI identity is captured before
  pool dispatch, not read from a worker's thread-local communicator registry.
* Both radial rules (even mesh), coupled s channels, p channels, multiple
  species and a species with no projectors: sparse values and lazily rebuilt
  columns match the eager columns bitwise. Repeated complex apply constructs
  once; an injected callback error clears the atomic-weights scope.
* Initial ionic-potential arrays match bitwise with root-only transformation.
  At seven ranks the cubic fixture retains **24,000..24,720 bytes/rank** of
  complex representatives versus **340,560 bytes/rank** for the full table;
  the skewed fixture **26,400..27,360** versus **374,640 bytes/rank**. These
  exclude species reconstruction factors and the distributed real table.

Final required `scripts/check.sh fast dft-physics dft-engine dft-numerics dft-cli`
returned **1**: **build PASS, fmt PASS, clippy PASS**; **537 tests passed,
3 failed, 149 ignored**, 54 test binaries, 12 s gate wall. The failures are all
unchanged CLI cache-policy tests:

* `pseudo::cache_resolution_and_no_network_side_effects`
* `pseudo::cached_fetch_checks_hash_and_provenance`
* `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`

They create their cache under `std::env::temp_dir()`. With the required
worktree-local TMPDIR, `dft-cli/src/pseudo.rs::outside_repository` correctly
rejects that location (`pseudopotential cache must be outside repositories`).
Those files are not owned by this ticket. No guard was weakened, test skipped,
or outside-worktree scratch used to conceal the failures. Gate log:
`target/scratch/gates-acceptance.log`. Its optional `ps` load diagnostic also
hits sandbox `Operation not permitted`; the script continues through all gates.
The first compile attempts exposed ordinary type/import errors (corrected);
the first gate additionally caught constant-size chunk Clippy warnings (fixed).

### Mac observation and retained design (O, not adoption)

Both final screens returned **0** using the existing `fast_screen_mpi.py`,
production mode, ladder on, 8 simulated ranks, 3 alternating pairs, phase
profile on, delivered orbitals off. These are converged runs, not fixed ten-step
runs. A separate JSON check (exit 0) enforces the ticket's **1e-12 Ha** cap
(the harness's standard-class cap is looser), matching SCF/Davidson/H counts,
repeatability within each arm, sparse byte counts and zero B full builds.

| Cell | A/B SCF, Davidson, H columns | Maximum abs dE (Ha) | Complex table bytes/rank A -> B | Additional B recipe bytes/rank |
| --- | --- | --- | --- | --- |
| si32 | 11, 67, 1888 (identical) | 2.842170943040401e-14 | 32,903,168 -> 2,056,704 | 1,028,224 |
| si64 | 10, 51, 3165 (identical) | 2.842170943040401e-14 | 132,218,880 -> 8,264,192 | 2,065,920 |

Bytes are construction totals across the fine and coarse operators, averaged
across ranks by the simulated-rank harness. B has zero lazy full builds on all
six retained runs. Density hashes differ A/B (rounding class), but repeat
identically within each arm. Neither an RSS claim nor a bit-identical A/B
trajectory claim is made.

Bucket medians, seconds (max-over-ranks within each run):

| Bucket | si32 A | si32 B | si64 A | si64 B |
| --- | ---: | ---: | ---: | ---: |
| init total | 0.110691 | 0.121854 | 0.972810 | 0.463038 |
| states | 0.023665 | 0.017899 | 0.092870 | 0.035792 |
| atomic start | 0.060914 | 0.077682 | 0.795783 | 0.378280 |
| densities | 0.005392 | 0.004369 | 0.013599 | 0.007474 |
| rest | 0.002642 | 0.002815 | 0.007860 | 0.005093 |
| basis | 0.013194 | 0.014511 | 0.036464 | 0.030530 |
| cli | unavailable | unavailable | unavailable | unavailable |

**No speedup established.** The harness bypasses the CLI (`init.cli` is zero,
not a measured zero cost). si32 total-wall samples A = 2.175988 / 2.128672 /
2.297577 s, B = 2.231019 / 2.143808 / 2.332140 s, load 3.70..5.45. si64 samples
A = 7.610771 / 18.702073 / 18.766105 s, B = 7.329880 / 13.016739 / 16.623432 s,
load 4.48..17.81: strongly nonstationary, so the si64 median directions cannot
be attributed causally to the change. No builds from this session overlapped
these retained runs. The switch stays OFF; dev08 decides adoption.

The first 256-G/eight-slot variant reduced storage but increased atomic-start
cost: si32 `atomic.hq` 0.066642 s, si64 0.466868 s (Mac observations). Retained
1024-G fixed tiles and omission of unused final slots reduce the **additional**
projection tree payload: si32 **6,291,480 -> 1,310,728 bytes/rank**, si64
**41,943,080 -> 9,437,200 bytes/rank** (including failure-flag payloads; global
shape, not material or rank count, sets the arithmetic). Final `atomic.hq`
medians are 0.018392 s and 0.100040 s, respectively, with the load caveat above.
Back application still reconstructs tiles on each rank for its band subset;
there is a memory/compute tradeoff and this is not proven to improve init wall.
The initial exploratory si32 screen overlapped a diagnostic rebuild at the end
and is excluded; the retained screens use `--binary` after builds finish.

Binary SHA-256: `c8a1af409695d53d8b88020eaf31c51c23a6c90d2d037fd6957e7c5048da1469`.
External Si UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Raw scratch records (not committed):
`target/scratch/si32-1024.json` SHA-256
`32154f37fddc7fd22e0e6ea43413af57daf363426e2b269bfa6988e439945ee0`;
`target/scratch/si64-1024.json` SHA-256
`c16b3dd6b08d811d39fbe19c6d0c36164efc581780aeabb40d8d0f85df446573`.
`target/scratch/acceptance-summary.json` holds the explicit post-check summary.
No QE calculation was run or claimed.

Mac observation recipe (external UPF, simulated ranks, observation only):

```sh
BAANDI_PHASE_PROFILE=on M2_PROFILE_ORBITALS=off \
python3 -I scripts/fast_screen_mpi.py --binary "$SCREEN_BINARY" --upf "$DFT_SI_UPF" \
  --case si64 --ranks 8 --pairs 3 --mode production \
  --a BAANDI_INIT_OWNED_STATES=off --b BAANDI_INIT_OWNED_STATES=on \
  --json-out "$TMPDIR_990/si64-owned-states.json"
```

Dev08 confirmation is for the master session after implementation and review;
the following is a command template, not an executed benchmark. Set `BIN` to
the reviewed release MPI binary, `INPUTS` to matching si32/si64 inputs with
`scf.max_iterations = 10`, and `OUT` to a fresh result directory. Preserve the
6c31696-default configuration, identical input/UPF/FFT/bands/start/thresholds,
one thread per rank, and record binary/input/UPF hashes and command exits.

```sh
for cell in si32 si64; do
  for pair in 1 2 3; do
    for arm in off on; do
      env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BAANDI_PHASE_PROFILE=on \
        BAANDI_INIT_OWNED_STATES="$arm" \
        mpirun -np 8 --map-by core --bind-to core \
        "$BIN" run "$INPUTS/$cell-10step.toml" \
        --output-dir "$OUT/$cell-$pair-$arm"
    done
  done
done
```

Report max/min-over-ranks and three-pair medians/spread for total `init`,
`init.states`, atomic-start total and parts, `init.rest`, `init.cli`,
`init.basis` (also `init.densities` to detect shifted work); energy differences
in Ha, SCF/Davidson/H counts, certificate status, full-build counts and per-rank
table bytes. Require no full build on the default production path in arm B.
These dev08 commands were not executed here. The master must use matched
10-step inputs; the local screen runs converged production SCF with the ladder
and therefore is not a reproduction of the fixed-step dev08 timing protocol.
| F-INIT-STATES-DEFAULT-ON | 2026-11-19 | BAANDI_INIT_OWNED_STATES default on. dev08 (master, 9cc6a7a, np=8, crossed A/B, floor first): si64 ten-step -11.14 % (floor 0.1118 s, effect 8.7x floor), si64 converging input -11.36 % (floor 0.0274 s, effect 43x floor, intervals separate), tail included; per-rank complex projector table 132.2 MB -> 8.3 MB (-93.7 %); si16/si32 np4/si-skew np4 neutral (+0.09/+0.23/-0.27 %, unresolved), no regressions; syn-spd, syn-p3, MgO(2-atom PBE smoke), si32 np=1 bit-identical; si64 np8 |dE| 2.84e-14 Ha, si32 np4 1.4e-14, si16 np8 7.1e-15, si-skew np4 3.6e-15 (standard 5e-9, strict 1e-10); SCF/Davidson/H counts identical, all converged. | ADOPTED | Gain appears at si64 size (memory/cache); small cells neutral. Rounding class: bits can differ at the 1e-14 level. |

## Distributed cutoff handoff pricing (G-HANDOFF-OWNED, 2026-11-19)

Initial pricing status: **implemented; prototype initially blocked on ownership**.
The orchestrator subsequently approved the narrow Gamma/session interface extension;
the completed default-off prototype and its measurements are recorded below.
Base `999f5fcb695e7b1bef8ad58d0963e604e5f526cc`, uncommitted worktree.
This implements the price-first part of plan-2026-11-19 section 5. It does not
repeat LADDER-IMPORT-SPLIT or LAZY-ASSEMBLE-FUSED and does not change the default
handoff or an SCF numerical rule.

`BAANDI_PHASE_PROFILE=on` exports `handoff.*` through the existing engine M2
counter/stage report (CLI `parallel.m2`, simulated harness `m2_counters` and
`stages.m2.handoff.*`). Off: no new clock reads, counter updates, snapshots or
allocations; only the disabled checks remain. The nine exclusive stage entries
are materialize (orbital gather, unpack and returned clone), orbital_copy,
density_copy (FFT prolongation and charge/positivity check), encode (electronic
state), setup (fine session construction/initialization and seeding), import
(density checks and coefficient decoding), check (the unchanged distributed
norm/overlap check), export_density (owned input gather), export_copy (remaining
coarse-state copies). Nested spans subtract child time. The additional
`materialize_collective` entry is **inside** materialize: never add it again.
Coarse SCF maps, coarse session setup, final delivery and atomic truncation
comparison are excluded. Fine session setup uses the same path even without a
tolerance seed; checkpoint and direct-min imports retain their entry point.

`cutoff_ladder::handoff_profile()` exposes rank-local cumulative times/counts.
M2 exports process sums, which the simulated harness divides by P: the following
are **rank means**, not a maximum/critical-path measurement. Counts are exact.
The gather-byte counter is the orbital allgatherv's delivered buffers
`8 * (P + Nc + Nc * bands)` per materialization and rank; it is not physical
wire traffic, does not include density/setup gathers, and is zero for cached
materialization. Collective counts include the entire handoff. Full-matrix bytes
are `16 * (Nc + Nf) * bands`. The coefficient-only live-payload upper bound is
`max(40*Nc*bands, 16*Nc*bands + 32*Nf*bands)` bytes per rank; it excludes density,
index metadata, allocator overhead, operators and solver workspaces and is
**not a measured RSS peak**. The first term deliberately bounds unpack plus
clone conservatively. A terminal coarse result already materialized by its last
map has zero handoff materializations; no map cost is relabelled as handoff.

Mac observations (O), pure release, 8 simulated ranks, production standard,
verification on, ladder on, seed 601 / fractional displacement 1e-3, **one A/B
pair per cell**. A = phase profile off, B = on. Shared-host load varied; wall
changes between A and B are not speedup claims. Source fixture and all remaining
controls are the existing `m2_entry_cost` production harness. UPF SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Measured executable SHA-256
`36d83663a0af93f18753b1a8abcb85fbd3470e294e20241baa9b97f9e54d7f29`.

| Observation / count, per rank | Si32 P=8 | Si64 P=8 |
| --- | ---: | ---: |
| Profile-on whole-run wall (s) | 1.742104 | 7.882185 |
| Exclusive handoff rank mean (s) | 0.032431 | 0.107617 |
| Rank-mean handoff / whole-run wall | 1.8616% | 1.3653% |
| Materialize incl. returned clone (s) | 0.001966 | 0.009546 |
| Collective time inside materialize (s; not additive) | 0.000357 | 0.001650 |
| Orbital copy (s) | 0.002002 | 0.006551 |
| Density prolong/check (s) | 0.003572 | 0.008731 |
| Encode (s) | 0.001313 | 0.006452 |
| Fine setup excluding import (s) | 0.016498 | 0.036032 |
| Import excluding overlap check (s) | 0.002960 | 0.010684 |
| Norm/overlap check (s) | 0.003881 | 0.028825 |
| Export density / remaining copy (s) | 0.000175 / 0.000064 | 0.000665 / 0.000131 |
| Orbital materializations | 1 | 1 |
| Orbital gather delivered bytes | 2,312,392 | 8,947,240 |
| Handoff allgatherv / tree_sum calls | 9 / 1 | 9 / 1 |
| Handoff alltoall / bcast / barrier calls | 0 / 0 / 0 | 0 / 0 / 0 |
| Full matrix payload bytes | 17,479,808 | 68,175,360 |
| Coefficient live-payload upper bound bytes | 30,401,984 | 118,590,912 |
| Fine SCF / Davidson / H columns, both arms | 11 / 67 / 1888 | 10 / 51 / 3165 |
| Energy (Ha), both arms | -126.703285939970 | -253.378117630305 |
| Profile-on/off energy difference (Ha) | 0, identical bits | 0, identical bits |
| Profile-on/off density FNV | identical | identical |

**Go condition holds through the avoidable gather**, despite both observed
rank-mean fractions being below 2%. The ceiling is the measured exclusive
handoff, never the historical 1.2007 s coarse stage. Even that ceiling includes
necessary fine preparation and density interpolation that an orbital permutation
alone will not remove. Real MPI latency and critical-path pricing remain dev08
work for the master session.

Reproduction (scratch only on the worktree volume; no UPF copied/committed):

```sh
mkdir -p target/scratch
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export TMPDIR_990="$PWD/target/scratch" TMPDIR="$PWD/target/scratch"
cargo test --offline --release -p dft-engine --test m2_entry_cost --no-run
python3 -I scripts/fast_screen_mpi.py --binary target/release/deps/m2_entry_cost-27ed8b949e92fddf --upf ../pseudo/Si.pz-vbc.UPF --case si32 --ranks 8 --pairs 1 --a BAANDI_PHASE_PROFILE=off --b BAANDI_PHASE_PROFILE=on --json-out target/scratch/price-si32.json
python3 -I scripts/fast_screen_mpi.py --binary target/release/deps/m2_entry_cost-27ed8b949e92fddf --upf ../pseudo/Si.pz-vbc.UPF --case si64 --ranks 8 --pairs 1 --a BAANDI_PHASE_PROFILE=off --b BAANDI_PHASE_PROFILE=on --json-out target/scratch/price-si64.json
# Equivalent direct harness controls (production/ladder are also script defaults):
DFT_SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF" M2_PROFILE_SERIAL=off M2_PROFILE_CASE=si64 M2_PROFILE_RANKS=8 M2_PROFILE_MODE=production M2_PROFILE_LADDER=on BAANDI_PHASE_PROFILE=on cargo test --offline --release -p dft-engine --test m2_entry_cost -- --ignored --nocapture
```

Master-session dev08 command (not executed here; use the current approved Si64
production input, matching PP/grid/bands/start/threshold and one thread/rank):

```sh
env MKL_NUM_THREADS=1 BAANDI_PHASE_PROFILE=on \
  mpirun -np 8 --map-by core --bind-to core \
  target/release/dft-cli run ~/fastbench/si64.toml \
  --output-dir ~/fastbench/G-HANDOFF-OWNED-si64-p8-profile
```

Record external monotonic wall, input/binary/UPF hashes, actual starts,
verification work, per-rank closure, and the new exclusive handoff fields.
Current CLI M2 fields are rank-local (root exports its rank), not a cross-rank
maximum; preserve that qualification and use the per-rank accessor if collecting
all ranks. Do not sum unrelated rank maxima. No QE run/comparison was performed.

Initial implementation boundary (resolved by the subsequent scope approval): `LazyOrbitals` and its private `LazyInner` live in
`crates/dft-engine/src/gamma_real.rs`, outside this ticket's OWNED_FILES. There is
no local-row accessor or constructor, only full materialization. The fine solve
also consumes its private matching/round-trip interface there. A direct route
requires a minimal, reviewed extension in that file to export/construct owned
ORC storage and match the fine layout. Accessing private storage from the owned
files or disguising a gather as a transfer is not a valid prototype. Scope
extension was requested; pending it, no `BAANDI_M2_OWNED_HANDOFF` switch or partial
transfer is installed. The full-matrix route stays unchanged. Follow-up must
preserve sqrt(2)*(x/sqrt(2)) on nonzero ORC rows, G=0 and signed-zero conventions,
conjugate pairing, agreed validation failures before redistribution, import
checks, unsupported-consumer fallback and checkpoint behavior.

Asynchronous coarse/fine preparation and fork/join are **not implemented**;
section 5's last paragraph needs its own eight-physical-core economic screen.

Validation at the final source (scratch `target/scratch`, uncommitted):

- `cargo test --offline --release -p dft-engine --test owned_handoff -- --include-ignored --nocapture`: exit 0, 2 tests passed (0.19 s). Profile off/on agree in every energy, density, eigenvalue and delivered complex-orbital bit, history H counts and across ranks; ranks 1 through 7, cubic/skewed cells, even/odd fine grids, active and terminal coarse boundaries. Profile off produces zero new metrics. Active lazy boundaries record one gather; the terminal retained audit already materializes in its map, so its handoff records zero. Initial test expectation incorrectly required one gather at that terminal boundary; it was corrected to match that existing consumer behavior, not by changing the engine.
- `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli`: exit 1; build, fmt and clippy passed; 453 passed, 3 failed, 134 ignored in 50 test binaries. All three failures are existing CLI pseudopotential-cache tests: `pseudo::cache_resolution_and_no_network_side_effects`, `pseudo::cached_fetch_checks_hash_and_provenance`, and unit `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`. They put their cache in `std::env::temp_dir()`, which the required scratch path places inside the repository; the production cache policy then rejects with `pseudopotential cache must be outside repositories`. No scratch relocation, test bypass, production-policy relaxation or unowned test edit was made. Reproduced the unit separately: `TMPDIR="$PWD/target/scratch" cargo test --offline -p dft-cli --bin dft-cli`, exit 101, same rejection.
- `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics`: exit 0, all four gates pass, 322 passed / 0 failed / 94 ignored across 37 binaries, 42 s wall. Both check invocations printed a sandbox `ps` permission error in their informational load probe; this did not prevent the gates from running.
- Both `fast_screen_mpi.py` commands above: exit 0; exact energy bits, density FNV and SCF/Davidson/H work agree between profile off/on. Timing differences between those arms are observations only.
- `python3 -I scripts/tests/test_readme_figures.py`: exit 1, 8 passed / 1 failed. Unchanged `test_structural_scope_and_sources` still requires the old phrase `standard profile is unsupported`, absent from current generated structural text. README/source evidence/scripts were not edited; no new physical validation result or speedup is claimed. README generation and dashboard/session JSON updates are outside this ticket's ownership.
- `git diff --check`: exit 0. No commit, remote command, QE run, async fork/join, default flip or new dependency.

Final direct `m2_entry_cost` production+ladder Si64 P=8 run: exit 0, wall 12.201468 s, rank-mean exclusive handoff 0.124560 s (1.0209%); the same 1 orbital materialization, 8,947,240 gathered orbital bytes/rank and 9 allgatherv + 1 tree_sum calls/rank. This run overlapped other local verification work and is not a timing comparison with the earlier A/B pair. The first direct invocation used a relative UPF path and failed before SCF (exit 101); Cargo runs the test from the crate directory, so the documented direct command now uses an absolute UPF path.

### G-HANDOFF-OWNED prototype after approved scope extension

Status: **CANDIDATE, implemented, default OFF; dev08 measurement pending**.
The orchestrator approved narrow `gamma_real.rs` interfaces and session plumbing.
`BAANDI_M2_OWNED_HANDOFF=on` retains the coarse lazy handle, gathers only its
integer row ownership metadata, computes one coarse/fine Miller permutation,
and moves ORC row payload directly to the owning fine ranks. New fine rows are
zero-filled locally. `off` is the full-matrix oracle and remains the default.
The registry has 111 switches. No `scf.rs` change was needed.

Implementation and validation contracts:

- All transfer/import policy lives in `cutoff_ladder.rs`. Gamma-real supplies a
  read-only local matrix/row/half-sphere/communicator view and a constructor for
  the existing lazy handle. It does not expose or modify the Hamiltonian or its
  potential. The existing first-solve row-layout match and potential identity
  checks remain in place.
- Coarse/fine representative and conjugate Miller indices, the zero mode, shape,
  finite values, source ownership coverage/uniqueness and layout identities are
  checked before the coefficient exchange. Local failures are reduced to the
  same error before ranks can enter differently shaped collectives. Eligibility
  and lazy-handle availability must agree too; partial availability rejects.
- Source storage can be in any local order and ownership can change between
  cutoffs. One plan enumerates each source/destination edge in global coarse ORC
  order. One all-to-all exchanges the packed row payload. The current Collective
  API has **fixed-size all-to-all only**, so messages are padded to the largest
  edge. Only populated coarse rows are logical payload; fine zero rows are never
  sent. The reported sent bytes include **all padding and self messages**, not a
  claim of sparse alltoallv or physical network wire traffic. An alltoallv adapter
  would require a separate transport change and measurement.
- Raw ORC values are copied without arithmetic. The existing warm-start
  `sqrt(2) * (x / sqrt(2))` convention is applied once, exactly as in the full
  complex oracle; G=0 is not scaled. The lazy materializer remembers newly added
  pairs so they return `(+0,+0)` complex zeros, matching `prolong_orbitals`, not
  the `-0` mirror imaginary value that generic Hermitian unpack would introduce.
- The density interpolation, charge/positivity/occupation checks and owned
  density installation are unchanged. Norm/orthogonality admission is still
  **1e-8**. It sums complex-unpack-equivalent local products with a deterministic
  tree. A conservative `64 * eps * N * max(1, sqrt(norm_i*norm_j))` rounding
  envelope detects proximity to the threshold; ambiguous cases materialize and
  use the original sequential-complex-row/split-column oracle. The envelope is
  only a fallback criterion, never a relaxed acceptance tolerance. The Gram
  reduction is rounding-class; the accepted fine start and delivered results
  tested here are bit-identical.
- A single rank, strict SCF, disabled lazy/distributed operator, declared
  checkpoint/runtime consumer, snapshot capture, or an already materialized
  terminal coarse block uses the existing full route. Exact checkpoint code and
  formats are unchanged. Full materialization of the transferred block is also
  bit-identical to the complex injection oracle (including signed zeros).
- Profile-off does not update the new metrics. Profile-on adds exclusive
  `handoff.transfer` (plan, metadata, exchange and eligibility check) to the
  previous nine parts. `materialize_collective` remains a **nonadditive subset**.
  The actual lazy materializer now records gathered-byte dimensions at its call
  site, so an ambiguous import fallback is counted correctly too.
  `owned_transfers`, `sent_bytes_with_padding`, `populated_payload_bytes` and
  `check_fallbacks` are counters. `owned_peak_bytes_sum` is the process sum of
  rank-local high-water **coefficient/message/Gram buffer payloads**; the harness
  divides by P. It excludes layout metadata, allocator overhead, density and
  solver/operator storage and is not process RSS. Use the per-rank
  `handoff_profile()` accessor for individual rank values. No cross-rank maximum
  is fabricated from the process-average harness or root-only CLI report.

Final-code validation:

- Release `owned_handoff` integration tests: 4 passed. Ranks 1..7, cubic/skewed,
  odd/even fine grids, active and terminal boundaries, production fine SCF to
  convergence with its certificate, exact energy/density/eigenvalue/orbital bits
  and H counts against the off route. Explicit fallback tests cover a declared
  checkpoint consumer, lazy off and strict SCF. A consumer declaration on only
  one rank yields the same error on every rank before the coarse solve.
- Release library `owned_handoff_tests`: 3 passed. Pure permutation/materialized
  complex bits and packed warm-start bits on ranks 2..7 and both grid parities,
  reversed rows and unrelated source/fine ownership, G=0, conjugate partners,
  signed zero and subnormal values; rank-local duplicate rows, NaN and mismatched
  fine layout all agree on an error; defects below, at and above 1e-8 match the
  original import verdict. Rank 1 takes the integration-tested full fallback.
- Final `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli`: exit 1, build/fmt/clippy pass, **458 passed / 3 failed / 134 ignored** in 50 binaries, 20 s. The same three known CLI cache-policy failures described above are the only failures; the orchestrator confirmed they are unrelated and pass on dev08. The informational `ps` load probe remains sandbox-denied. No tests were skipped or altered to hide these failures.
- The initial low-level fault test accidentally assigned no rows to rank 0 for
  P=3, so no fault was injected. Its ownership fixture was corrected to guarantee
  a rank-0 row; all error-agreement tests then pass. No production behavior or
  tolerance was changed to satisfy that test.

No asynchronous fork/join was implemented. It still needs its own economic
screen on the same eight physical cores. Nothing here changes the ladder,
convergence rules, solver defaults or the performance acceptance authority.

Final Mac screen (O), pure release, 8 simulated ranks, **three alternating pairs
per cell**, both arms phase-profile on; A=OWNED_HANDOFF off, B=on. Same production
standard/verification/ladder controls, UPF and seed as the pricing screen above.
Binary SHA-256 `2cc600c80069aacd71d86bd391df3f96058de8ea55fdace7aa51dd060047a416`.
Raw records remain in `target/scratch/owned-{si32,si64}-final.json` (uncommitted
scratch); commands below reproduce them. Times below are medians of rank-mean
exclusive handoff sums, **not sums of independently median-selected parts**.

| Metric (per rank unless stated) | Si32 A -> B | Si64 A -> B |
| --- | ---: | ---: |
| Exclusive handoff (s) | 0.060205 -> 0.034637 | 0.155803 -> 0.106210 |
| Whole-run wall median (s; observation only) | 2.028417 -> 2.022573 | 7.304849 -> 7.245540 |
| Orbital materializations | 1 -> 0 | 1 -> 0 |
| Orbital gather delivered bytes | 2312392 -> 0 | 8947240 -> 0 |
| Whole-run allgatherv delivered bytes | 33528112 -> 31249616 | 109060088 -> 100180504 |
| Whole-run allgatherv calls | 252 -> 252 | 199 -> 199 |
| Handoff allgatherv / tree_sum / alltoall calls | 9 / 1 / 0 -> 9 / 13 / 1 | 9 / 1 / 0 -> 9 / 13 / 1 |
| All-to-all sent bytes, including padding | 0 -> 940032 | 0 -> 2433024 |
| Populated row payload bytes | 0 -> 284852 | 0 -> 1109988 |
| Full-route coefficient payload peak bound / owned-route tracked buffer peak | 30401984 -> 2972552 | 118590912 -> 9127008 |
| Ambiguous overlap-check fallbacks | 0 -> 0 | 0 -> 0 |

The observed exclusive savings are **25.568 ms (Si32)** and **49.593 ms (Si64)**.
Both are below the original price-first budgets of 32.431 / 107.617 ms. The
original approximately 1.4–1.9% handoff ceiling remains the planning constraint;
**none of the historical 1.2007 s coarse stage is credited**. The new off-arm
handoff samples vary (Si32 31.971–61.009 ms; Si64 125.192–195.494 ms), so these
Mac medians are not stable production estimates. Whole-wall within-arm spreads
are 0.485851 / 0.646316 s, much larger than the A/B median wall differences
(0.005844 / 0.059309 s). No whole-run speedup or default adoption is claimed.
The distributed check still costs work, and agreement adds reductions; the first
fine solve also performs its existing lazy layout check. Only the explicitly
measured handoff difference is reported as a handoff observation.

All six runs per cell have identical energy bits, density FNV and SCF/Davidson/H
counts: Si32 11/67/1888, Si64 10/51/3165. Energy difference is **0 Ha**. The final
off-arm Si32/Si64 fingerprints also match the earlier price-first off-arm records,
so the unchanged default has an independently saved pre-prototype oracle. The
full-state synthetic tests additionally compare eigenvalues and all delivered
complex coefficient bits. No QE execution or cross-model comparison is involved.

Final reproduction (both screen commands exit 0):

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export TMPDIR_990="$PWD/target/scratch" TMPDIR="$PWD/target/scratch"
cargo test --offline --release -p dft-engine --test m2_entry_cost --test owned_handoff --no-run
python3 -I scripts/fast_screen_mpi.py --binary target/release/deps/m2_entry_cost-27ed8b949e92fddf --upf ../pseudo/Si.pz-vbc.UPF --case si64 --ranks 8 --pairs 3 --a BAANDI_PHASE_PROFILE=on,BAANDI_M2_OWNED_HANDOFF=off --b BAANDI_PHASE_PROFILE=on,BAANDI_M2_OWNED_HANDOFF=on --json-out target/scratch/owned-si64-final.json
python3 -I scripts/fast_screen_mpi.py --binary target/release/deps/m2_entry_cost-27ed8b949e92fddf --upf ../pseudo/Si.pz-vbc.UPF --case si32 --ranks 8 --pairs 3 --a BAANDI_PHASE_PROFILE=on,BAANDI_M2_OWNED_HANDOFF=off --b BAANDI_PHASE_PROFILE=on,BAANDI_M2_OWNED_HANDOFF=on --expect m2.handoff.materializations:down --expect collectives.allgatherv.bytes:down --json-out target/scratch/owned-si32-final.json
cargo test --offline --release -p dft-engine --test owned_handoff -- --include-ignored --nocapture
cargo test --offline --release -p dft-engine --lib owned_handoff_tests -- --nocapture
```

Master-session dev08 follow-up (commands **not executed here**): use the approved
MPI+MKL environment and the same current production inputs; first price Si64 P=8,
then repeat Si32 P=8. Keep all starts, PP/grid/bands/threshold, binding and thread
counts identical. Preserve input/UPF/binary hashes, verification verdict/work,
energy bits and density hashes, external monotonic wall and the raw results.
Use the normal load gate before **each** launch. The following is a single-arm
command template; execute off/on then on/off then off/on as three pairs in unique
output directories. No remote installation or input changes are requested.

```sh
cargo build --offline --release -p dft-cli --features mpi,mkl
sha256sum target/release/dft-cli ~/fastbench/si64.toml ~/fastbench/si32.toml
env MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 BAANDI_PHASE_PROFILE=on BAANDI_M2_OWNED_HANDOFF=off \
  mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si64.toml \
  --output-dir ~/fastbench/G-HANDOFF-OWNED-si64-p8-off-1
env MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 BAANDI_PHASE_PROFILE=on BAANDI_M2_OWNED_HANDOFF=on \
  mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si64.toml \
  --output-dir ~/fastbench/G-HANDOFF-OWNED-si64-p8-on-1
```

Read `parallel.m2.counters["handoff.*"]` and the named entries of
`parallel.m2.stages`. Sum exactly `materialize`, `orbital_copy`, `density_copy`,
`encode`, `setup`, `import`, `check`, `export_density`, `export_copy`, `transfer`;
do not add `materialize_collective`. Root CLI exports root-local handoff values;
retain this qualification or collect the public per-rank snapshots in the master
harness. Report bytes including padding and the additional reduction latency.
Then repeat the paired launches with `BAANDI_PHASE_PROFILE=off` for uninstrumented
launch-to-exit acceptance. A local timing benefit does not flip the default.

| F-VERIFY-SCHEDULE | 2026-11-19, Mac pure Rust T=1, HEAD dba7d824 plus uncommitted patch; 8 systems (si8/16/32/skew, MgO8 PBE, syn-gaussian-fcc/p3/spd) x seeds 601–603 x thresholds 5e-9 / 1e-10 Ha. Default-off off/trace/on retry scheduling; criterion, stop rule and `verify_fixed_point` byte-identical. Cost proxy m per map, 1.3m per attempt, wait above failure probability 1/2.3. | **NEGATIVE / NO-GO: only 3 systems benefit, fewer than required 4; do not promote.** Explicit on retained for reproduction, default off. | Initial pre-failure geometric predictor: 14 proposed waits, 8 would defer a pass, rejected before on. Revised forecast V_last * S_now / S_at_failure, future contraction from current S ratio and (when available) last-attempt ratio: 5/5 positive trace votes are actual failures; predicted extra maps exact/short/long = 3/1/1. Real final A/B, 48 pairs: H including coarse start 48,944 -> 48,676 (-0.548%); attempts 85 -> 81, failures 37 -> 33, maps after first attempt 37 -> 37. Four improved runs, 44 identical work, zero H/extra-map regressions; all 5e-9 rows unchanged. Changes: si8 s602 tight H 928 -> 860, si16 s601/s602 tight 1959 -> 1870 / 1943 -> 1835, syn-gaussian-fcc s602 tight 417 -> 414. Every system median paired H change 0%. Pristine/off and off/trace scientific identities 48/48 each. All results pass a real attempt; max error vs same-geometry strict reference 1.4631e-10 Ha (5e-9 class), 8.1855e-12 Ha (1e-10 class). | Failed-attempt half-reuse changes the continuation, so pre-failure rates overpredict waiting and saved attempt work is not net benefit. Re-anchor after actual maps and re-evaluate each boundary; no first-attempt forecast without a calibrated scale. Positive sample count is only five; no probability-bound claim. Mac counts only, no dev08 timing. Final check: build/fmt/clippy pass, 460 tests pass / 3 known in-repo TMPDIR pseudo-cache failures / 134 ignored. README check has one pre-existing stale-phrase failure. | docs/design/verification-schedule.md (complete 48-row table, calibration, hashes, scratch harness, crossed dev08 commands); raw uncommitted target/scratch/ab and attempts.csv. |
### G-HANDOFF-OWNED: dev08 rejection and collective-budget audit

The orchestrator relayed the following dev08 result for `70f01fd`: 10 crossed
pairs in a clean window, Si64 MPI8, all 10 pairs slower with owned handoff,
mean whole-run increase **0.278 s / 3.12%**. `between.ladder_handoff` fell from
0.0967 to 0.0352 s; orbital gathered bytes fell from approximately 8.9 MB to
zero and materializations from one to zero. The mechanism worked but the
candidate **failed the production speed screen**. These are relayed numbers,
not runs performed in this worktree; raw launch/input/binary hashes were not
provided in this message. Keep `BAANDI_M2_OWNED_HANDOFF` default OFF. The previous
Mac results do not establish adoption. The orchestrator limits a revision to
at most two handoff collectives and will close the ticket on another regression.
The available whole-run saving ceiling is only about 0.06–0.10 s.

Source audit on worktree HEAD `84643cf` (the implemented handoff, preceding the
relayed dev08 integration commit): the 13 tree sums are **not per-column or
per-chunk loops**. All are fixed transaction-level calls:

| Site | tree_sum calls |
| --- | ---: |
| `run_ladder`: eligibility agreement | 1 |
| `warm_owned_handoff`: availability, density-present agreement | 2 |
| `prolong_owned`: local shape/finite check, header identity, partition size, ownership/plan validity, Miller mapping | 5 |
| `initialize_owned_ladder`: imported shape, overlap input validity, Gram reduction, overlap verdict, density-import verdict | 5 |
| Total | 13 |

The same measured scope also contains **9 allgatherv calls** and **1 all-to-all**:
three orbital ownership/header exchanges, one coarse owned-density gather, three
fine-session dummy mixer-validation gathers (the mixed value is discarded and
history reset), and two fine `GridPlan::new` stick-ownership exchanges. Thus the
current normal-path handoff is **23 collectives**, not 13. Batching scalar flags
alone cannot meet a two-collective *whole-handoff* limit. Even eliminating all
orbital metadata/validation calls leaves six pre-existing density/setup gathers.
Do not move them before/after timer boundaries or relabel them to pass the gate.

A fresh Mac audit on the unchanged code reproduced the exact 9/13/1 counts:
`M2_PROFILE_MODE=production M2_PROFILE_LADDER=on`, Si64, 8 simulated ranks,
`BAANDI_PHASE_PROFILE=on BAANDI_M2_OWNED_HANDOFF=on`. One observation, no A/B claim:
wall 11.024827 s, exclusive rank-mean handoff 0.157631 s, energy
-253.37811763030507 Ha, zero orbital gathered bytes and zero materializations.
This is a baseline reproduction, **not a revised implementation or improvement**.
Command (exit 0; raw output `target/scratch/batch-audit-si64.log`):

```sh
TMPDIR="$PWD/target/scratch" DFT_SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF" M2_PROFILE_SERIAL=off M2_PROFILE_CASE=si64 M2_PROFILE_RANKS=8 M2_PROFILE_MODE=production M2_PROFILE_LADDER=on BAANDI_PHASE_PROFILE=on BAANDI_M2_OWNED_HANDOFF=on cargo test --offline --release -p dft-engine --test m2_entry_cost -- --ignored --nocapture
```

Required revision design (not yet implemented): a fixed-size descriptor/status
exchange establishes safe payload dimensions and shared eligibility; a single
owner-directed payload exchange carries coefficients together with repeated
validation/Gram contributions and the coarse density blocks. Every rank sees
the same validation contributions, including flags from sources whose orbital
rows it does not own. In particular, receiver-local row-coverage checks cannot
return early before peers agree. The density gather must disappear into that
payload, the discarded initial mixer validation must be performed locally with
the same parameter/data predicates, and the fine-grid stick ownership must be
computed locally from the deterministic partition rather than communicated.
Existing norm-threshold ambiguity fallback must remain collectively agreed;
its extra oracle communication cannot be silently excluded from the counters.
The requested count assertion must sum *all five* operation counters over the
handoff on ranks 2..7, not just `tree_sum` or the transfer function.

Scope dependency: `GridPlan::new` is in `crates/dft-engine/src/grid_dist.rs`, which
is not in OWNED_FILES or the previously granted gamma_real/scf/session extension.
Its two stick-ownership allgathervs must be removed for the fine-handoff
construction to fit this budget. A narrow constructor scope extension was
requested, specifically to calculate those identical metadata arrays from the
existing deterministic partition. No change to grid arithmetic, FFTs or ordinary
map behavior is proposed. The scope question is pending; no numerical code,
collective counters, default, test expectation or timing boundary was changed
in this audit, and no second dev08 candidate is ready.

Audit verification: `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0
scripts/check.sh fast dft-engine dft-numerics dft-cli` exited 1: 458 passed,
3 failed, 134 ignored; build, fmt and clippy passed. The three failures are
the acknowledged CLI cache tests rejecting an in-repository temporary cache;
raw output is `target/scratch/batch-audit-check.log`. No revised collective-count
test is claimed: the unchanged implementation still performs 23 calls.

### G-HANDOFF-OWNED: batched transaction after the dev08 rejection (2026-11-19)

**Status: implemented, default OFF; awaiting the orchestrator's second dev08
screen. No speed adoption claim.** The user approved the narrow `grid_dist.rs`
constructor extension. The old owned implementation used 13 tree reductions,
9 allgathervs and 1 all-to-all (23 calls). The revised normal handoff uses
**one allgatherv and one all-to-all, zero tree reductions, two total calls**.
All counts include export, transfer, import and fine session initialization;
no initialization work or collective was moved outside the measured spans.

Transaction and error contract:

- The first exchange combines eligibility/availability, integer-layout and
  consumer fingerprints, local row/payload lengths, finite/shape flags and all
  overlap partial sums. Its extent uses the *completed coarse SCF's* common
  band count, not the local lazy handle's potentially invalid column count.
  A malformed local handle, including a wrong band count, cannot change the
  descriptor exchange size. The coarse solve already requires common bands
  and numerical inputs; this does not add a new cross-rank input contract.
- Reconstruct the existing binomial-tree sum locally from those contributions.
  The second exchange sends populated coefficient rows directly to their fine
  owners and repeats only ownership metadata and coarse density planes to all
  peers. Every rank checks the same complete row partition, density lengths and
  Gram data before consuming its destination-only coefficients. Missing blocks,
  duplicate/missing rows, nonfinite coefficients, local consumer/layout/band
  mismatches and malformed density blocks produce an agreed error.
- A norm/overlap within the existing rounding envelope of the 1e-8 admission
  boundary selects the original sequential complex oracle on **every rank**.
  That decision is made after the first exchange: only in this rare case the
  second payload replicates coefficients as well. It still takes **two calls**,
  preserves the original verdict, and records fallback/materialization bytes.
  The ordinary path has zero orbital materializations and gathered bytes.
  This is an oracle fallback, not a relaxed tolerance or changed norm convention.
- The initial mixer trial is discarded by the original constructor. The owned
  import constructor validates its plain Pulay parameters and starting density
  locally; no dummy distributed charge reductions are needed. Fine GridPlan
  metadata uses the existing deterministic partition and the same sorted stick
  order. This constructor branch is scoped to owned handoff setup and the ON
  switch, independently of profiling: coarse initialization and default OFF
  retain their original communication. Grid arithmetic, potential state,
  checkpoint fallbacks and the ordinary session constructor are unchanged.

Validation (all scratch below `target/scratch`): release `owned_handoff` **5
passed**, including the 1..7 rank/cubic/skew/odd/even matrix, completed fine
certificates, profile on/off and full-route bit identity, unsupported consumer
fallbacks and a named whole-handoff <=2-calls test on ranks 2..7. Release
`owned_handoff_tests` **4 passed**: arbitrary/reversed source ownership, Gamma
packing/G=0/signed zeros/subnormals, rank-local faults and missing blocks/bad
density, and norm-boundary oracle verdicts. The boundary test uses an independent
transport wrapper to assert exactly two calls, including boundary-oracle cases,
on every rank count 2..7. Wrong-band descriptors are also tested. Final
`scripts/check.sh fast dft-engine dft-numerics dft-cli` exited **1**: **460 passed,
3 failed, 134 ignored**; build/fmt/clippy pass. The same three acknowledged CLI
cache tests reject the in-repository scratch location. No unrelated cache fix.
Raw logs: `batch-owned-final.log`, `batch-unit-final.log`, `batch-final-check.log`.

Mac observations, three crossed pairs per case, production ladder, eight
simulated ranks, same binary SHA-256
`a9ba16717189f4cab15ecd1f4f20adc44fd933f6dd767d671e444291fadb2e86`.
A is full route (switch OFF); B is revised owned route (ON). Handoff seconds
are the median of each run's **sum of the ten exclusive handoff parts**, with
`materialize_collective` excluded from that sum. These harness numbers are
rank means, not sums of per-part maxima.

| Observation | Si32 OFF -> ON | Si64 OFF -> ON |
| --- | ---: | ---: |
| Exclusive handoff median (s) | 0.034584 -> 0.033829 | 0.179371 -> 0.095521 |
| Median per-run handoff share (%) | 1.944 -> 1.872 | 2.329 -> 1.214 |
| Orbital gathered bytes/rank (B) | 2,312,392 -> 0 | 8,947,240 -> 0 |
| Orbital materializations | 1 -> 0 | 1 -> 0 |
| Handoff allgatherv/tree_sum/alltoall | 9/1/0 -> 1/0/1 | 9/1/0 -> 1/0/1 |
| ON submitted transfer bytes/rank, padding included (B) | 1,615,728 | 3,946,288 |
| ON populated coefficient bytes/rank (B) | 284,852 | 1,109,988 |
| ON mean rank-local peak buffer bound (B) | 4,287,624 | 12,014,816 |
| Whole-run allgatherv delivered bytes/rank (B) | 33,528,112 -> 29,414,088 | 109,060,088 -> 96,938,704 |
| Whole-run median (s) | 1.795946 -> 1.802706 | 7.691011 -> 7.847321 |
| Whole-run median change (%) | +0.376 | +2.032 |
| Harness wall spread (s) | 0.149924 | 0.243944 |

The transfer byte counter includes the submitted descriptor/Gram and padded
all-to-all send buffer. It is an API payload count, not a measurement of MPI's
network algorithm. Peak bounds include coefficient and transaction buffers,
exclude allocator/RSS and other SCF storage, and the simulated harness reports
an average of rank-local high-water bounds, **not the worst rank's RSS**.
No boundary oracle occurred in these screens. Energy bits, density hashes,
SCF/Davidson/H work are identical in both arms (Si32 -126.70328593997024 Ha,
density `0x928361d517f75f05`; Si64 -253.37811763030507 Ha, density
`0x37808e2041c875bd`). The OFF fingerprints match the previous implementation's
saved screen records too.

For comparison with the rejected owned implementation, the fresh pre-revision
Si64 ON audit above recorded **0.157631 s**, 13 tree reductions and 23 total
calls; the revised ON median is **0.095521 s**, zero tree reductions and two
calls, with gathered orbital bytes **0 -> 0**. These are different Mac windows
(single old observation versus three new pairs), not a paired speed estimate.
The earlier three-pair owned observations were 0.034637 s (Si32) and 0.106210 s
(Si64); they likewise are not a current timing control. The same-binary current
A/B whole-run medians are slower on the ON arm, within the observed spreads.
The harness's `CANDIDATE` label reflects its counter screen; it is **not** a
whole-run speed/adoption result. The dev08 ceiling remains the measured
approximately 0.06–0.10 s, never the 1.2 s coarse stage. A second dev08 regression
closes this attempt under the orchestrator's decision.

Reproduction (each release test and screen exited 0; the final check exited 1
for the known cache failures):

```sh
export TMPDIR_990="$PWD/target/scratch" TMPDIR="$PWD/target/scratch"
cargo test --offline --release -p dft-engine --lib owned_handoff_tests -- --nocapture
cargo test --offline --release -p dft-engine --test owned_handoff -- --include-ignored
cargo test --offline --release -p dft-engine --test m2_entry_cost --no-run
M2_PROFILE_MODE=production M2_PROFILE_LADDER=on python3 -I scripts/fast_screen_mpi.py --binary target/release/deps/m2_entry_cost-27ed8b949e92fddf --upf ../pseudo/Si.pz-vbc.UPF --case si64 --ranks 8 --pairs 3 --a BAANDI_PHASE_PROFILE=on,BAANDI_M2_OWNED_HANDOFF=off --b BAANDI_PHASE_PROFILE=on,BAANDI_M2_OWNED_HANDOFF=on --expect m2.handoff.tree_sum:down --expect m2.handoff.materializations:down --json-out target/scratch/batched-si64-final.json
M2_PROFILE_MODE=production M2_PROFILE_LADDER=on python3 -I scripts/fast_screen_mpi.py --binary target/release/deps/m2_entry_cost-27ed8b949e92fddf --upf ../pseudo/Si.pz-vbc.UPF --case si32 --ranks 8 --pairs 3 --a BAANDI_PHASE_PROFILE=on,BAANDI_M2_OWNED_HANDOFF=off --b BAANDI_PHASE_PROFILE=on,BAANDI_M2_OWNED_HANDOFF=on --expect m2.handoff.tree_sum:down --expect m2.handoff.materializations:down --json-out target/scratch/batched-si32-final.json
CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli
```

Master-session dev08 command templates (**not executed here**). Use the existing
approved MPI+MKL build/input and load gate; repeat ten crossed pairs, unique
output directories, with identical binding/start/PP/grid/bands/threshold and
certificate policy. Keep raw external wall, binary/input/PP hashes, scientific
fingerprints, handoff counters and exclusive part timings for every arm. The
acceptance conditions are zero handoff tree reductions, at most two handoff
collectives on the supported owned route, zero ordinary orbital gathers,
unchanged accuracy/work/certificate, and no repeat of the whole-run regression.

```sh
export TMPDIR_990="$PWD/target/scratch" TMPDIR="$PWD/target/scratch"
cargo build --offline --release -p dft-cli --features mpi,mkl
sha256sum target/release/dft-cli ~/fastbench/si64.toml
# One crossed pair; orchestrator repeats/reverses order under its load gate.
env MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 BAANDI_PHASE_PROFILE=on BAANDI_M2_OWNED_HANDOFF=off \
  mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si64.toml \
  --output-dir "$PWD/target/scratch/owned-handoff-si64-p8-off-1"
env MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 BAANDI_PHASE_PROFILE=on BAANDI_M2_OWNED_HANDOFF=on \
  mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run ~/fastbench/si64.toml \
  --output-dir "$PWD/target/scratch/owned-handoff-si64-p8-on-1"
```

Async coarse/fine fork/join remains unimplemented and needs its own economic
screen. Changes remain uncommitted for Claude's review.
| G-HANDOFF-OWNED-CLOSED | 2026-11-19 | BAANDI_M2_OWNED_HANDOFF closed (default off stays). dev08 crossed pairs, clean window: 70f01fd si64 np8 +3.12 % (10/10 pairs slower, handoff.tree_sum 1 -> 13); cd33d3f (collectives <= 2: tree_sum 13 -> 0, allgatherv 9 -> 1) still +1.27 % (9/10 pairs slower, |mean|/se 4.31); si32 np8 -0.27 % unresolved. Handoff itself falls 66 % (0.097 -> 0.032 s) and gathered bytes 8.9 MB -> 0, but sent bytes with padding are 3.55x the populated payload (3.95 MB vs 1.11 MB) and owned peak bytes +32 %, so the whole run is slower; ceiling was only 0.06-0.10 s of 9 s. | REJECTED | Rule: a remaining resolved regression closes the ticket. Possible future premise: padding-free exchange (alltoallv with exact counts) would have to beat the plain gather in whole-run seconds. |

| F-VERIFY-SCHEDULE-LAZY-TRACE | Follow-up to dev08 41d2d79 (user-relayed crossed A/B): one core Si64 -7.06%, all 4 pairs faster, abs(mean)/SE 218, attempts [11,12,13] -> [11,13], H 4008 -> 3871, Davidson 71 -> 69, abs(dE) 4.11e-11 Ha; one core Si32 +0.28%, all 4 pairs slower, abs(mean)/SE 4.87, no failures and identical bits/counts; MPI np=8 Si64 unchanged. | Default OFF; adoption pending master's post-fix dev08 measurements. Prior Mac NO-GO is historical, not a veto on the new dev08 evidence. | Scheduler moved entirely inside eligible verification branch; first attempt bypasses prediction, later predictions require a failed real attempt. Cold helper, scalar allocation-free predictor, audit journal allocated once at first attempt and preallocated for the fixed budget; no ordinary-map schedule branch/call/counter/logging. Optional `verification.schedule` in JSON for trace/on, absent off: map, predicted n, actual action, would_wait, measured error, actual-attempt ratio and pass; waits have null measured/pass fields. Mac 48 off/on/trace triples (144 SCFs): each arm numerically identical to prior implementation, 48/48 off/trace identities. No-failure subset 19 geometries / 120 fine maps: zero predictor calls, identical science/counts, 14,017 total H per arm. Small syn-p3 seeds 601–603: standard maps 5/5, H 13/13, one attempt, zero predictions; tight maps 6/6, H 14/14. Full on journal: 81 attempts + 4 waits, 37 post-failure predictor calls. | Rounding/path class, not path-identical: omitting a failed solve changes reuse/continuation; dev08 abs(dE) 4.11e-11 Ha is accepted inside strict 1e-10. Trace remains numerically inert. Zero per-map work is not zero reporting cost; JSON construction and journal allocation still occur at attempt/output boundaries. No local wall-time adoption claim. Check build/fmt/clippy pass; 463 tests pass, 3 known pseudo-cache TMPDIR failures, 134 ignored. | docs/design/verification-schedule.md (field contract, Mac counters, relayed dev08 numbers, crossed commands plus schedule.csv); target/scratch/followup-ab and followup-summary.json. |
| F-VERIFY-SCHEDULE-DEFAULT-ON | 2026-11-19 | BAANDI_VERIFY_SCHEDULE default on. dev08 crossed pairs, clean windows: one core si64 -6.81 % (6/6 pairs negative, |mean|/se 177.5; earlier 41d2d79 -7.06 %): attempts [11,12,13] -> [11,13], failures 2 -> 1, Davidson 71 -> 69, H 4008 -> 3871, |dE| 4.11e-11 Ha (path class within 1e-10; not path-identical); MPI np8 si64 no effect (schedule unchanged, -0.18 % unresolved); one core si32 control +0.29 % (8/8 pairs positive, |mean|/se 6.4, bits identical, prediction_evaluations 0: the cost is bookkeeping about 2.3 ms per map, not prediction arithmetic); trace records verification.schedule correctly (predicted n=1, waited, pass at map 13). | ADOPTED | Net positive at the target sizes; small cells pay about 0.3 %. Follow-up: locate the 2.3 ms/map bookkeeping. Delivered state can differ at the 1e-11 Ha level when an attempt that would have failed is skipped. |

## G-SCHED-OVERHEAD (2026-11-19): ordinary-map allocation audit

Baseline `b58a5a35523419a2ae829fd39c36a1a326307e23`, default-on schedule;
uncommitted candidate, Mac pure release, one thread. **Implemented, dev08
performance confirmation pending.** The relayed Si32 +0.29% (23 ms / ten maps)
does not locate that time in ordinary-map bookkeeping: both source inspection
and a test-only counting allocator find **zero added allocations per nonterminal
ordinary map in the baseline already**. No schedule-specific history/spectrum
clone, formatting, JSON, timer, environment read or branch runs there.

The actual redundant work was at the first eligible verification and output:
reserve `max_iterations` diagnostic rows, push/fill the first row, then clone
the journal into the result. The candidate bypasses scheduling before any
failure, creates retry storage on the first failure, creates only one output
row when the first attempt passes, and moves an existing journal into the
terminal result. No policy, certificate, arithmetic, tolerance or default change.

Counting allocator, real seven-map skew/nonlocal SCF: off counts
`[440,59,60,45,58,57,96]`; baseline on/trace
`[440,59,60,45,58,57,98]`; candidate on/trace
`[440,59,60,45,58,57,97]`. Ordinary added allocations stay **0**; final reporting
allocations fall **2 -> 1**, with first-success capacity **100 -> 1**. The
required public diagnostic Vec/JSON is not allocation-free. The regression
test checks this distinction, full scientific identity, and a live allocator
probe. A helper test proves the first boundary needs neither history access
nor journal creation. The forced-failure trace rows and predictor fallbacks
remain covered; the synthetic wait test now clears its real fixture's
inherited default-on journal before replacing its history.

Eight alternating off/on pairs collected for each build/case; conservatively
exclude the potentially overlapping last baseline Si32 / first baseline Si16
pairs. Results below are **on minus off mean ± sample SE, ms per fine map**,
where ordinary time is steps minus certificate wall; ten fine maps in each run.

| Build | Case | Pairs | Ordinary mean ± SE (ms/map) | abs(mean)/SE | Shared-host load |
| --- | --- | ---: | ---: | ---: | --- |
| Before | Si32 | 7 | -0.566346 ± 3.405864 | 0.17 | 4.25–6.25 |
| Before | Si16 | 7 | +0.811029 ± 3.256482 | 0.25 | 4.46–4.95 |
| Candidate | Si32 | 8 | -5.159820 ± 3.574622 | 1.44 | 8.12–20.04 |
| Candidate | Si16 | 8 | -1.092127 ± 2.031926 | 0.54 | 6.43–8.12 |

All deltas are unresolved. No before/after-window subtraction or speedup claim.
Si32 is the no-failure control (attempt 10, zero predictions); Si16 has one
failure, attempts [9,10], one prediction and no wait. All 64 timing runs have
identical scientific results across builds and off/on (timing/journal metadata
excluded), including density bytes and work. Si32 E = -126.70343204998852 Ha;
Si16 E = -63.129095464190456 Ha. These are internal identity checks, not QE
references. The predictor and complete verification function regions are
byte-identical to baseline.

Validation: `CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli`
exit 1: build/fmt/clippy PASS, **466 passed / 3 known CLI in-repository TMPDIR
cache failures / 134 ignored**, 51 binaries. Final release schedule tests:
4 passed; allocation gate: 1 passed. README freshness suite: 8 passed / 1
pre-existing scope-text expectation failure (`standard profile is unsupported`);
README and figure freshness assertions themselves pass. No out-of-scope fix.

Lesson: allocating at a terminal boundary is not per-map allocation, and dividing
whole-run wall by map count is not attribution. This patch removes measured
redundant storage but does **not** establish the origin or elimination of the
relayed 23 ms. Master remeasures one-core Si32 eight alternating pairs and
one-core Si64 crossed A/B, preserving [11,13] and trace records. If the regression
persists, inspect verification/allocator-cache/output effects on dev08.
Detailed operations, hashes, phase-bucket timings, commands and raw evidence
locations: `docs/design/verification-schedule.md`, `target/scratch/summary.json`.

Failed-attempt same-mode baseline/candidate identity: six off/on/trace pairs
(12 SCFs), Si16 at 1e-10 Ha and Si64 at 5e-9 Ha, seed 601. Si16 on [10,12]
with a wait at 11; Si64 on **[11,13]** with a wait at 12. Off/trace remain
[10,11,12] and [11,12,13], respectively. Both predictor counters remain 2;
all trace rows, science and density bytes are identical to the corresponding
baseline arm. Scratch `identity-summary.json` SHA-256
`8ad80c4b33b825993f131d1cba60548ca0a789942622df380aaaea11632a102c`.
| G-LADDER-MINGRID-CLOSED | 2026-11-19 | Codex implemented BAANDI_LADDER_COARSE_GRID=default|min (minimum valid smooth coarse FFT grid) and measured 7 systems x 3 seeds plus simulated ranks 1-7: the existing default already selects the smallest valid 2/3/5-smooth grid, so `min` is bit-identical to default (42 rank/standard/strict conditions, 2,430 captured files identical, all grids/maps/H counts equal, W ratio 1.0, max strict-reference difference 8.5e-10 Ha < 5e-9). The earlier F-LADDER-HANDOFF offline "minimum valid grid" gain used non-smooth sizes that are slower in wall time. | CLOSED, not merged | No code merged (the switch would be a no-op). Coarse-stage cost needs a different lever. Branch task/G-LADDER-MINGRID kept in ../rust-dft-wt for reference. |

## F-LADDER-FRACTION — closed negative rule study (2026-11-19)

- **Scientific result: NO-GO for the tested general rules; defaults unchanged.**
  Mac pure-Rust one-core release, seven systems (si8/16/32, si-skew, MgO8,
  synthetic skew and slab), seeds 601–603. All **126/126** matrix arms
  converged and passed unchanged verification. No new QE run or dev08 timing.
- Mechanism: 32 fraction-screen arms, eight systems including si64, seed 601,
  fractions 0.5/0.45/0.4/0.35; **288 fine-map records** include standard error,
  actual eigen tolerance, H, Davidson work, density residual and Pulay records.
  Handoff diagnostics include low/high-G density error, Hartree/XC response,
  variational orbital energy error and occupied-subspace principal angles.
  Si32/0.4 crosses the 1e-4 residual stop after six coarse maps (9.94947e-5),
  versus seven at 0.5/0.45 (~2.6e-5). Auto first tolerance loosens about 27x
  (9.59e-8 to 2.61e-6 Ha), and fine maps rise 10 to 13. Si64 takes six
  coarse maps at every fraction, keeps its seed near 4e-6 Ha, and fine maps
  fall 13 to 12 at 0.4 despite comparable representation-error growth.
- Causal two-by-two si32/601 ablation: coarse state 6/7 maps and tolerance
  seed from state 6/7 yield fine maps **13/11/11/10** (state6/seed6,
  state6/seed7, state7/seed6, state7/seed7); fine+verification H
  **2233/2471/1966/2017**. Tightening the seed alone reduces maps but can
  increase work. The Ritz-floor ablation takes ten maps but 5206 fine H.
  Shared first-Pulay bounce is not a loss classifier; history transfer was
  not changed or repeated. High-G omission is real but cannot alone explain
  the opposite si32/si64 map-count changes.
- Independent default-off run-config switches:
  `BAANDI_LADDER_FRACTION_RULE=grid` fills the requested fraction's FFT
  plateau by the first excluded density Miller vector (no fitted/material
  constants); `BAANDI_LADDER_STOP_RULE=modes` scales the requested residual
  by sqrt(nonzero coarse density modes / half-cutoff modes), capped at one
  (Parseval per-mode budget, explicitly not a discretization-error bound);
  `BAANDI_LADDER_WATCH=on` records early error above the first-two-map
  contraction forecast after completion and never discards or repeats work.
- Against 0.5, paired median **H-eq** wins: plain 0.4 **4/7**, grid **5/7**,
  modes **4/7**, combined **5/7**. Systems with at least one fine-H or
  fine-map increase: **5/7, 5/7, 5/7, 6/7**, respectively. Raw unweighted
  total-H wins: **2/7, 2/7, 1/7, 2/7**. No rule meets the generality gate.
  Grid H-eq ratios by system: si8 0.908663, si16 0.999826, si32 0.980393,
  si-skew 0.873109, MgO8 1.046156, syn-skew2 0.891416, syn-slab 1.123627.
  Si32 seed 603 still takes twelve fine maps versus eleven at 0.5.
- Supplemental si64 seed 601: baseline H-eq **5474.453125**; grid
  **4560.699219** (-16.6912%), modes **4552.455078** (-16.8418%), each
  twelve versus thirteen fine maps. The large-system win survives in this
  geometry, but does not remove the multi-system counterexamples. H-eq is
  grid-volume weighted; extra basis setup/algebra need actual dev08 wall.
- Worst |E-Estrict| **4.310998447110137e-9 Ha**, syn-slab/603 grid and
  combined, below 5e-9 with only **6.89001552889863e-10 Ha (13.78%)** margin.
  Baseline worst 8.463689749049763e-10 Ha. No accuracy or verification
  criterion was changed; the narrow margin supports no default adoption.
- Identity: 21 exact baseline offline/production controls and 84 candidate
  production-driver controls; 21 watch/off numerical/history/handoff-array
  pairs; 14 trace/off pairs; four pre-rule si32 fraction histories/array
  comparisons. New unit tests cover grid preservation/idempotence/Parseval
  budgets and exact watch-off delivered density/orbitals. The diagnostic
  forecast flags 11/21 lower-fraction runs and is not predictive enough to
  select a corrective policy.
- Required check exit **1**: build/fmt/clippy pass, **468 passed, 3 failed,
  134 ignored**, 51 binaries. The three existing CLI cache tests reject
  the required in-repository TMPDIR; isolated bin-test rerun exits 101 with
  the same reason. No new test failure. Python analysis validates all matrix
  and mechanism records; scripts and complete per-geometry tables in
  [ladder-fraction.md](../design/ladder-fraction.md). Raw evidence under
  `target/scratch/fraction-*` must be preserved before cleanup. No timings
  accepted from the loaded Mac, no remote changes, no default changes.
- **Lesson:** coarse-stop crossing and its tolerance seed interact. Neither
  geometric basis size nor a residual budget alone predicts nonlinear Pulay
  recovery; removing fine maps is insufficient when it increases fine H.
  Dev08 crossed-pair commands are documented for characterization only;
  the master decides further work and adoption.

## G-LADDER-FALLBACK — predictable coarse basis rejection (2026-11-19)

| ID | Change / status | Measured result | Lesson / evidence |
|---|---|---|---|
| G-LADDER-FALLBACK | Default-off `BAANDI_LADDER_BASIS_CHECK=on`: decide dimensional incompatibility before coarse problem construction; typed engine rejection plus stable CLI switch path `skipped_coarse_basis`, numeric fallback reason. No cutoff, bands, solver, tolerance or default changes. | Mac pure release T=1, base a8347ddf: eight reported errors reproduce (coarse PW/bands 1/6 on d3/fft-primes/p3/s1/s2/spd, 7/8 on cubic-degenerate/z1); original odd-shift is already skipped for multiple k points, so its reported dev08 error needs original input. Requested nine x seeds 601–603 x strict/standard: 54/54 converge; on/off/unset/no-ladder energy difference 0 Ha, density hashes and total H identical; no newly enabled ladders. Gaussian-fcc positive control retains 27 coarse PW / 8 bands and the same ladder. Additional independently displaced nine: 53/54 converge, existing syn-s2/602 strict budget exhaustion unchanged in every arm. Combined 240 pristine/candidate scientific identities and 120 on/off/unset identities. | Planning/diagnostic defect, not a 48-PW minimum: Davidson capacity is clipped to Npw and B=Npw is legal. Fix avoids futile setup, saves zero H; no wall or adoption claim. Required gate exit 1: build/fmt/clippy pass, 470 tests pass, 3 known in-repo TMPDIR cache failures, 134 ignored. README test has existing stale-phrase failure. Syn-s2 cross-profile multiplicity precludes absolute-reference claim; paired absolute-dE classes strict 1e-10 / standard 5e-9 Ha remain unchanged. Full numbers, exploratory failures, commands, hashes and limits: [ladder-fallback.md](../design/ladder-fallback.md); raw scratch `target/scratch/{baseline,candidate,displaced-nine-baseline,displaced-nine}`. No QE run, no dev08 acceptance, leave uncommitted. |

### LADDER-FRACTION-DEV08-SCREEN (2026-11-19, dev08 via master; ledger by orchestrator)
Fraction 0.4 (default 0.5): crossed pairs wall -8.43 % si8 np4, -7.41 % si16 np8, -2.88 % si-skew np4, -2.33 % si32 np8, **+9.71 % si32 1core**, -16.85 % si64 np8, -12.83 % si64 1core. Mechanism (si32): coarse stage margin m = stop/final residual = 4.0 at default; at 0.4 one coarse map disappears (m 1.03), last coarse std_err 6.8e-5 -> 3.15e-3, first_tol=auto loosens 46x, fine 10 -> 13. m >= 3.7 separates 9/9 default-arm pairs where fine steps rise (risk predictor, not a gain predictor). Raising BAANDI_LADDER_COARSE_STOP alone is refuted (si32 np8 2e-4 -0.94 %, 1core +6.43 %, 4e-4 +4.91/+18.14 %; si64 4e-4 +6.21 %; si16 -3.32/-4.54 %). syn-gaussian-fcc (m 7.7) -1.8 %/-2.4 % unresolved and si-low (m 3.4) +0.4..+1.1 % unresolved: the coarse basis does not shrink there (27 plane waves at both fractions), so nothing changes. Survey: only 12 of 34 suite cases converge with a running ladder (k-point cells and nine synthetic cells fall back; the synthetic ones are a plan-time basis shortage, G-LADDER-FALLBACK). Mac rule study (F-LADDER-FRACTION, 126 arms): no rule meets the generality criterion; worst |dE| 4.31e-9 Ha (13.8 % margin to standard). Verdict: default stays 0.5; candidate two-condition gate (coarse basis shrinks at the candidate fraction, known before SCF) + risk (margin) is a hypothesis, not adopted. Lesson: ten-step and Heq indicators predict magnitude only where the basis actually shrinks.

### SCOREBOARD-2026-11-19 and LADDER-PLATEAU (dev08 via master; ledger by orchestrator)
Scoreboard at d6be1d0 (default equal to f78fb30; unmeasured-gate note: two default-off commits only), uninstrumented, first iteration discarded, certificate verified, one case/model per job: si64 MPI8 B 10.6081 s (in-window spread 1.0 %), Q 4.3578 s (rmm), B/Q 2.4343, B-0.90Q 6.6861 s; si32 MPI8 B 2.1998 s (3.0 %), Q 1.2789 s (rmm), B/Q 1.7201, B-0.90Q 1.0488 s; si32 1core B 8.1631 s (0.1 %), Q 4.7251 s (david), B/Q 1.7276, B-0.90Q 3.9105 s. si64 absolute is window-bound (earlier window 8.8801 s for the same input and binary family); never compare table-to-table.
Plateau pairs (same coarse FFT grid): si32 np8 0.40 -> 0.41 wall -2.01 % (coarse maps 6 -> 7, fine 13 -> 11, H 1900 -> 1993); si32 0.40 -> 0.45 -6.03 %; si64 np8 0.40 -> 0.41 **+28.86 %** (coarse 6 -> 6, fine 11 -> 13, H 3044 -> 3710, |dE| 2.05e-10). Verdict: no plateau; cutoff changes where the coarse stage stops and how many fine maps follow. Heq sign failed again when coarse and fine map counts move in opposite directions (si32 0.41: Heq +4.2 %, wall -2.0 %). F-VERIFY-COST A/B closed (superseded by F-VERIFY-SCHEDULE).

## G-PHASE-VALIDATE (2026-11-19): exact finite-dot range proof

| ID | Change / status | Measured result | Lesson / evidence |
| --- | --- | --- | --- |
| G-PHASE-VALIDATE | Default-off `BAANDI_PHASE_VALIDATE=on`: owned-projector validation scans finite component maxima once (G per constructor, positions per species); maxima <= 2^500 prove every three-product dot finite. Otherwise original loop/errors. Factor/norm checks, actual phases, lazy table and collectives unchanged. Rank/run-local constructor/call/pair/accept/fallback counters; exclusive nanoseconds only under PHASE_PROFILE. | Mac arm64 pure release, Si64 15/7.5 Ha: exactly 2 constructors and 2,065,920 checks per owned rank; on 0 checks / 2 fast accepts / 0 fallbacks. One-thread constructor A0,A1,B0,B1,A2, 12 retained groups: off median 0.000780146 s, on 0.0000717085 s; paired median saving 0.0007075205 s, off/off SD 0.0000703289222 s, 10.06x noise (resolved >=3x). Off scan 3.77028% of constructor wall, not full engine init. Actual one-core engine path: zero calls, zero saving. Two simulated ranks separately confirm counts; scan/init fractions 1.21008%/1.21761% off -> 0.107290%/0.099016% on; simulations are not MPI timing. | Implemented; no default or whole-run speed claim until dev08. 30 numerical edge cases, 160 off/on constructor pairs across 4 cell shapes / 2 k types / 2 radial rules; 13,360 complex coefficients bit-identical, exact errors and lazy-once behavior; existing eager coefficient oracle also passes with switch on. Required gate exit 1: build/fmt/clippy PASS, 559 passed / 3 known in-repo TMPDIR cache failures / 152 ignored. Real MPI unavailable; per-rank critical wall and complete init denominator remain for master dev08 A/B. Proof depends only on binary64 range, never material/cutoff. Raw spans, commands, hashes, reproducible harness and scope: [phase-validate.md](../design/phase-validate.md). No QE source or QE run. |
## F-LADDER-2X2 — captured-state/seed causal test (2026-11-19 UTC)

**Closed for a universal tight-seed rule; defaults unchanged.** Unlike the earlier
0.40 six-/seven-map state ablation, this crosses the actual natural 0.50/0.40
handoff density AND warm orbitals with both observed seeds. Mac pure release,
T=1, si32 and held-out size si64, fractional displacement ±0.001, seeds 601–603.
The source base is `d6be1d0cab87cb0654d6b6de2c45f8c9a1fe6cfb`; no QE/dev08 run.

- **72/72** continuations and **18/18** unsplit natural controls converge;
  unchanged standard verification passes. **48** diagonal/repeated-control
  comparisons are byte-identical in energy/eigenvalues, density, orbitals and
  recorded scientific history. A 0.45 natural cross-check is kept for every
  geometry. The complete record has **824** maps and **30** unique cells.
- Si32/601, cells `(state50,seed50)/(state50,seed40)/(state40,seed50)/(state40,seed40)`:
  maps **10/11/11/13**, fine+verification H **1918/1987/2455/2233**. Tightening
  state40 removes two maps but adds **222 H (+9.94%)**. Loosening state50 adds
  one map and 69 H. The map/H interactions are **+1 / −291**.
- Si32/602: all four cells take ten maps and seed40 is already tighter than
  seed50. Si32/603: *loosening* state40 removes a map (12 → 11), H 2211 → 2064.
  Si64/601: state40 takes 12 maps/4061 H at either seed. Si64/602: tightening
  state50 worsens 12 → 13 maps, H 4556 → 4690. Si64/603 has a 40.3x seed jump
  at the 7 → 6 coarse-map crossing, yet tightening state40 leaves 12 maps and
  raises H **4154 → 4861 (+17.02%)**. These are geometry-dependent interactions,
  not support for a new controller. No default-off production rule was added;
  the conditional seven-system rule-validation stage was not triggered.
- Si64/603 state50+seed40 is the only unique first-redo cell: injected
  6.152031800255166e-6 Ha becomes effective **5.426896967241076e-7 Ha**; first-map
  cost 706 H / 6 blocks includes both solves. All resets/refits are zero.
  Failed verification attempts, per-map H/blocks, Coulomb error, charge and
  complete measured setup/fine/verified times are retained. Common-low-G density
  error and complete occupied principal-angle spectra compare each state to its
  geometry's strict reference with a common 0.40 PW domain.
- Largest absolute strict-reference difference **3.001474624397815e-10 Ha**,
  below the standard 5e-9 Ha gate; **17/30** unique cells also lie within strict
  1e-10 Ha. No assertion that all standard results meet strict accuracy. Maximum
  output charge error **2.2737367544323206e-12 electrons**.
- Test-only harness: reuse exact captures, explicit `first-seed:<Ha>` override,
  native-host seed discovery, alternating full-construction pairs and A/A floor.
  Pristine-parent si32 natural trajectories reproduce at 0.50/0.40/0.45. The
  minimal-grid arm now plans its old grid before one constructor, eliminating
  its discarded coarse problem; si16 default/minimal before-after arrays,
  histories, energies and H counters match. The si64 alias was already present
  at this base; a real-fixture regression test now protects it. Production source
  before the test module is byte-identical; no session/switch/default edit.
- **14/14** full-construction smoke arms pass, each rebuilt handoff matching the
  capture. One retained si32/601 pair has full seconds 4.3911/4.6240 for state50
  and 4.7793/4.4706 for state40 (seed50/seed40). Retained A/A is +1.836%; one pair
  is not acceptance. Matrix A/A drift reaches 39.8% on si32/602. All Mac wall
  values are observations; no small timing effect is claimed.

Full cells, all unique per-map records, commands, hashes, diagnostics and dev08
order-alternated/A/A recipe: [ladder-2x2.md](../design/ladder-2x2.md). Raw evidence:
`target/scratch/{2x2,baseline-natural,minimal-before,minimal-after,full-pairs-smoke}`;
archive before cleanup. No production optimization is adopted. **Lesson:** map
recovery can buy more H work, the seed ordering itself changes with geometry,
and the ordinary redo can change the delivered initial tolerance. Keep the
default and distinguish these effects before fitting any handoff policy.

Validation: final required gate exit **1**, build/fmt/clippy PASS, **470 passed /
3 failed / 135 ignored**, 51 binaries, 17 s. The three failures are the known
in-repository TMPDIR pseudo-cache restrictions. The initial gate also found an
ignored-driver `useless_vec`, fixed to a six-element array; a fresh final-build
si32/601 six-arm run reproduces **49** binary/history files byte for byte.
Python fixture/pair tests **3 passed**; complete matrix analysis and all extra
harness checks exit 0. README check separately has its existing stale-phrase
failure (8 pass / 1 fail). Source, logs, exact failures and final build hashes
are in the design note. The overall required gate remains red. All edits stay
inside OWNED_FILES, and changes remain uncommitted.

## LANE-INIT Phase 0 — profile request (2026-11-20)

User's revised phased workflow supersedes the earlier Mac-diagnose/implement
charter for this task. **PROFILE_REQUEST**, instrumentation only; no candidate
optimization, default change, adoption decision or speed claim.

- Read startup code, INIT-CLI/INIT-PACKED/F-INIT-STATES/G-PHASE-VALIDATE records
  and the backlog. Added default-off `BAANDI_INIT_PROFILE`: 20 nested call/timer
  pairs and 23 work/shape counters, local to the run/rank. `STAGE_TIMERS=off`
  retains counts without clocks. No per-G/atom timer, I/O or collective.
- Mac correctness only: Si8 and Si32, baseline / new-off / counts / timed,
  eight converged and verified runs; all energy components, scientific
  iteration histories, eigenvalues, verification and density bytes identical;
  maximum |delta E| = **0 Ha**. Off records zero in all new counters; counts
  records zero nanoseconds; timed and counts record identical work.
- Synthetic projector validation: **448** eager/owned constructor comparisons,
  cubic/skew/slab/elongated, Gamma/shifted k, both radial rules, s/p/d/f and
  repeated p projectors, selected-row partitions of size 1/2/4; coefficient
  bits, lazy state and counters checked. These are row partitions, not real
  MPI runs. No Mac timings support a performance conclusion.
- Parser checks missing/malformed records, nanosecond units, rank identity,
  payload floors, and missing-time semantics; parses actual Mac-generated
  results for all three switch modes. Existing global phase min/max remain
  separate from new rank-local times.
- Detailed verification, commands, limitations and pending dev08 jobs:
  [Phase 0 report](lanes/INIT/phase0-report.md),
  [profile request](lanes/INIT/profile-request.md).
  Remote measurement and any Phase 2 design are the orchestrator's next step.

Phase 0 gate result: `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0
scripts/check.sh fast dft-numerics dft-physics dft-engine dft-cli` exit **1**;
build/fmt/clippy PASS; 561 passed, 3 known CLI in-repo TMPDIR cache-policy
failures, 153 ignored. The `ps`/`nice` sandbox diagnostics did not prevent
the actual test binaries from running. README script exit 1: 8 tests pass,
existing stale structural-scope phrase test fails. No new scientific result
or README claim is added. Instrumented release build, eight-run identity
runner, owned/eager constructor test, switch contract test and parser tests
all exit 0. Changes remain uncommitted; stop at PROFILE_REQUEST.
## LANE-DAVIDSON Phase 0 — requested 2026-11-20

Phase outcome: **PROFILE_REQUEST**, not an optimization candidate. The revised
user instruction supersedes the earlier diagnose/optimize lane run. Exploratory
triangular-product code was removed; no change to arithmetic, tolerance,
subspace size, root selection, restarts or defaults is retained. No Mac timing
supports a performance conclusion.

Default-off `BAANDI_DAVIDSON_PROFILE` records bounded per-call real-Davidson
shapes, requested solve tolerance, inclusive durations, map/solve IDs and
rank-local tree-sum payload/site records. It distinguishes lean factor extension,
reduced transformation, coefficients, return Gram, full EVD, B5 attempt and
whole reduced-eigen step. Buffers are prepared outside the numerical loop;
overflow is explicit and makes the parser reject the profile. OFF has no clock
reads, trace storage or extra collectives (a flag branch remains). The CLI
collects traces only after SCF; JSON stays outside the solver.

Mac role: compilation, scientific bit-identity and parser transport checks only.
Evidence and final gate outcomes: `lanes/DAVIDSON/correctness.json` and
`lanes/DAVIDSON/validation.md`. Dev08 protocol and exact one-job-per-point
commands: `lanes/DAVIDSON/profile-request.md`. Parser:
`scripts/lane_davidson_parse.py`. No adoption recommendation; await dev08 data
before Phase 2. REL-SLICES, Advance and restart-keep remain closed/rejected.
## LANE-FFT Phase 0 (2026-11-20; revised user instruction)

PROFILE_REQUEST: instrumentation only, no optimization or performance conclusion.
The superseded optimization and Mac timing interpretation were removed. Added
construction-cached, default-off FFT_PROFILE under PHASE_PROFILE: 18 separate
local-potential stages, calls/points/logical bytes/conventional FFT work,
count-only support, no clocks/TLS writes at disabled kernel sites. Serial and
M2 plans keep numerical operations unchanged; shared phase/JSON wiring was
explicitly approved. Per-rank bucket records use a separate struct, and u64
metrics use exact limb encoding through the existing profile gather.
Mac correctness: Si8 OFF/timed/count-only energy difference 0 Ha, eigenvalues,
occupations, scientific histories and density bytes identical; all verified,
11 SCF iterations. 144 parser rows per enabled actual CLI file; integer
counters equal. Rank-1/2/3/8 synthetic M2 identity, serial apply/kept/density
exclusion and 3-rank wire-format tests pass. No Mac timing supports a speed
claim. [Profile request](lanes/FFT/profile-request.md),
[equivalence evidence](lanes/FFT/cli-equivalence.json),
[validation](lanes/FFT/validation.md). dev08 supplies the profile before Phase 2.
## LANE-VERIFY-PHASE0 (2026-11-20): PROFILE_REQUEST

User's revised phase instruction supersedes the optimization charter for this
run. No optimization/candidate selection retained. Default-off
`BAANDI_VERIFY_PROFILE` adds per-attempt/per-rank size, work and consecutive
phase records, plus actual scheduler-call records; PHASE_PROFILE is required,
STAGE_TIMERS gates every new clock. Post-run adapter gathering only; judgment,
ratchet, schedule and reuse are unchanged. The parser preserves per-call and
whole-run critical-rank closure rather than summing unrelated rank maxima.

Mac evidence is **correctness/plumbing only**, never performance: 48 synthetic
SCFs / 30 mode comparisons (cubic/skew/slab, 3/6 bands, 1/2 k blocks, simulated
2/3 ranks, forced failed attempts); five real-PP inputs x five arms (20 identity
comparisons), maximum dE **0 Ha**, density/scientific history identical. Four
real inputs converge; real slab exhausts identically in every arm, while the
synthetic slab converges. Parser 5 tests pass and actual Mac JSON parses; final
release Si8 five-arm smoke passes. Final required gate exit **1**: build/fmt/
clippy pass, **471 passed / 3 known in-repo TMPDIR cache failures / 136 ignored**.
Initial build/API, timer-granularity test, driver baseline-convergence assumption
and clippy failures were fixed and are documented, not omitted.

Status: **PROFILE_REQUEST**, no dev08 execution or speed claim. Exact five
case/model jobs, six retained profile/control cycles, schemas and parser paths:
[profile request](lanes/VERIFY/profile-request.md). Scope, hashes, numerical
identity and real exit codes: [validation](lanes/VERIFY/validation.md) and
[correctness record](lanes/VERIFY/correctness.json). Await the orchestrator's
dev08 profile before Phase 2. Changes remain uncommitted.
## LANE-NONLOCAL Phase 0 — PROFILE_REQUEST (2026-11-20)

Instrumentation only, no optimization/default change. The user's revised
instruction moved all profiling and performance conclusions to dev08; the
previous in-session candidate and probe source were removed. The orchestrator
subsequently approved instrumentation wiring in `real_gamma.rs` and engine
`profile.rs`. Six real-Gamma nonlocal child stages (projection, rank reduction,
coupling, back product, addition, audit) now reach the actual CLI phase JSON,
with per-bucket rank extrema, logical FLOP/byte inputs and width counts.
`NonlocalMark` is isolated from shared solver snapshot slots. OFF has no clock,
TLS/counter arithmetic, allocation, lock or I/O; only the enable branch remains.

Mac correctness only: Si8 and Si32, seed 601, fractional displacement 0.001,
one worker, off/counts-only/timed: **6/6** CLI runs converge and verify, **4/4**
off/on pairs have identical scientific f64 bits, density bytes, eigenvalues,
occupations, history, work and verification attempts; **|dE| = 0 Ha**.
Counts-only lane clocks are zero, timed and count-only work counts match.
Four actual CLI JSON files parse successfully. An unmodified counts-only CLI
result is retained as a regression fixture; **5/5** parser tests pass.
Simulated Si8 at two ranks additionally has identical energy bits, density FNV,
SCF/Davidson/H work off/on, with nonzero projection/reduction counts in coarse,
fine and verification buckets. This is not real MPI timing. Earlier synthetic
coverage: 45 operator/energy pairs bit-identical over cubic/skewed/slab shapes,
0/1/7 projector channels and 0/1/2/3/17 bands; 11 rank integration tests pass.

Required final `scripts/check.sh fast dft-physics dft-engine dft-cli`: **exit 1**,
build/fmt/clippy PASS, **480 passed / 3 failed / 131 ignored**. The failures are
the known CLI cache tests rejecting the in-repository TMPDIR. Optional `ps`
load reporting and `nice` also encounter sandbox restrictions; no numerical
failure is observed. The earlier physics-only gate exited 0 (92/0/18).
No Mac timing supports a performance conclusion, no QE/dev08 run, no adoption.

Exact dev08 commands, field semantics, limits and correctness reproducers:
[profile-request.md](lanes/NONLOCAL/profile-request.md). Evidence:
[equivalence.json](lanes/NONLOCAL/equivalence.json),
[simulated-equivalence.json](lanes/NONLOCAL/simulated-equivalence.json),
[report.json](lanes/NONLOCAL/report.json). Collection is qualified for one
worker per rank / one core; averages are per-bucket bounds, not event traces.
Stop at PROFILE_REQUEST until the master relays the dev08 profile.

## LANE-INIT Phase 2 — local atomic block (2026-11-20)

- **Implemented, default off, dev08 screen pending.**
  `BAANDI_ATOMIC_LOCAL_BLOCK=on` reuses the adopted guarded atomic CholQR2
  locally, with no multi-rank communicator. Two Gram GEMMs, two factor/inverse
  stages and two update GEMMs replace serial two-pass MGS. Existing guards
  and MGS fallback, completion, HQ and Ritz remain. MPI arithmetic is unchanged.
- New premise relative to A52: that rejection concerned Davidson correction
  panels; this is the fixed atomic input span, already validated/adopted on
  MPI under INIT-ATOMIC-BLOCK. The returned dev08 profile identifies serial
  MGS as the remaining target: si64 one-core atomic start 2.19395 s,
  orthonormalize 1.21804 s of external 42.9552 s. MPI8 start 0.52558 s of
  9.3551 s is separately reported, not used as a cross-model speed ratio.
  Inclusive times are not summed. Raw dev08 records were not supplied locally.
- Hypothesis: one-core si64 net saving 0.5–0.9 s if sequential GEMM sustains
  20–40 GFLOP/s; not a measurement. Alternative Gamma ORC-native span/HQ/Ritz
  deferred because it requires a broader representation change. No MPI gain
  predicted for the chosen switch. New INIT_PROFILE counts distinguish MGS,
  accepted/fallback local blocks, block Gram/factor/update/gather and floor sizes.
- Mac correctness (pure backend): full Si8, Si32, skew Si standard and strict
  SCFs; max cross-arm delta 2.842170943040401e-14 Ha (standard),
  1.4210854715202004e-14 Ha (strict). Standard vs strict reference max
  1.0450662557559554e-10 Ha, below 5e-9 Ha. Disabled candidate retains baseline
  scientific bits, including density/history; strict reference has no separate
  verification object, while all standard runs passed verification. Synthetic
  cubic/skew atomic-start SCFs and simulated P=2/3/5/7 pass with the switch on;
  local algebra tests include complex/Gamma and rank-deficient/fill failures.
  No Mac performance conclusion. Evidence and final command status are in
  [phase2-report.md](lanes/INIT/phase2-report.md) and the linked JSON.
- No default adoption or QE comparison. Required next evidence: positive MGS
  count -> zero, accepted local blocks with no fallback; same-window dev08 A/A
  then crossed A/B, negative mean wall difference with |mean|/SE >= 3 and >3x
  measured floor, correctness and later generality gates. See
  [diagnosis](lanes/INIT/diagnosis.md) and
  [dev08 request](lanes/INIT/dev08-request.md). Leave uncommitted and stop.

## LANE-DAVIDSON Phase 2 — upper-border GEMM tiles (2026-11-20)

Default-off `BAANDI_DAVIDSON_UPPER_BORDER`: compute upper-border rectangles in
64-column tiles and omit entries discarded by symmetric cleanup. Serial
real overlap/reduced and MPI lean S/H borders; no H, EVD, tolerance, capacity,
restart, collective count/payload change. Candidate 2 (deferred): omit zeros
in triangular reduced-coordinate products. REL-SLICES, Advance, restart-keep,
async alltoall and PR-CG remain rejected/closed.

Dev08 source 2733278 relay: Si64 MPI8 solve envelope 6.4323/9.7715 s,
reduced Gram 14.4%, rotation 13.8%, eigen 11.0% (inclusive); one-core
27.2416/45.5618 s, rotation 21.6%, reduced Gram 8.0%, overlap 7.8%.
No inclusive rows added. Numeric floor/gap diagnosis is **incomplete**:
raw block/count/payload traces and rate ceilings were not relayed; requested
from the orchestrator. The quoted exclusive shares total 104.1%, needing
aggregation clarification. Conditional Amdahl scenarios, not measured gains,
are in `lanes/DAVIDSON/diagnosis.md`.

Mac correctness only: Si8/Si32 base/OFF/ON and skew/Si8 strict OFF/ON science
and density bytes identical, abs(dE)=0 Ha; Si32 workers 1/3 identical.
Slab both exit 3 after 100 iterations with identical science (not a converged
generality pass). Direct m=137 border test across simulated ranks 1/2/3/4/8
agrees within 1e-12; 13 ON solver oracle tests pass. Final fast gate:
615 passed, 3 known TMPDIR-cache failures, 160 ignored; build/fmt/clippy PASS.
Counter/parser schema tested on Mac output; **no Mac performance conclusion**.
Status: implemented for dev08 screening, default OFF, floor diagnosis awaiting
raw data. Request pre-registers exact tile FLOP reduction, A/A floor then
crossed pairs, abs(mean)/SE >=3; all five requested points in separate jobs.
Evidence: `lanes/DAVIDSON/{phase2-correctness.json,phase2-validation.md,dev08-request.md}`.

### LANE-INIT-ATOMIC-LOCAL-BLOCK-DEFAULT-ON (2026-11-19/20, dev08 via master; ledger by orchestrator)
BAANDI_ATOMIC_LOCAL_BLOCK (guarded CholQR2 for the single-process atomic start; no effect with MPI ranks) default flipped on. dev08 crossed pairs, one core, 3 pairs, A/A floor first: si64 -1.98 % (|mean|/SE 17.99, floor 0.00707 s; gain in ladder_coarse 7.416 -> 6.584 s), si32 -1.09 % (14.91, 0.088 s vs 3x floor 0.0348 s); si16, si8, si-skew neutral (|mean|/SE 1.42/0.49/0.62, all negative sign); no resolved loss; |dE| <= 2.8e-14 Ha (si64 bit identical); SCF/H/verified identical in every pair. Counters moved as pre-registered (mgs.calls 1 -> 0, local_accepted 0 -> 1, fallbacks 0). Mac fast gate: 572 passed, 0 failed. Lesson: atomic start was 63-86 % of the init bucket (init 6-15 % of wall); init profile problem span covers only 12-45 % of init.
## LANE-DAVIDSON Phase 2 follow-up — corrected models / ROTATE_COLUMNS

User-relayed dev08 source2733278 on-1 counts and rates resolve the earlier
missing-data screen prerequisite. Withdraw the earlier 104.1% split: flat
phase-profile keys contained a non-leaf duplicate; corrected closures are
98.5–100.9%. Serial lane envelope27.242s and phase parent32.7484s are separate
sources. MPI1464 Gram calls/6.996s and808 rotations/6.715s aggregate eight
ranks (mean-rank0.8745/0.8394s), not wall. Serial rotation5.877s includes
residual0.504s. Detailed numbers/models: `lanes/DAVIDSON/diagnosis.md` and
`dev08-relay-model.json`.

Important arithmetic corrections: 23871x132x132 rotation is TWO GEMMs,
1.6637GFLOP ->26.01ms at63.96GFLOP/s, not13ms; mean53.92ms/model=2.073.
264-cubed GEMM at95.4GFLOP/s is0.386ms, not0.04ms; 9m^3 EVD model1.736ms,
MPI mean6.975ms/model=4.018. Dominant shapes and all-call means are not a
joint distribution, so these are screening models, not exact removable gaps.
REL-SLICES, bounded-subspace, restart-keep, PR-CG and async alltoall stay closed.

UPPER_BORDER is retained unchanged, A/B already queued. New default-OFF
DAVIDSON_ROTATE_COLUMNS routes serial non-fused Ritz pairs through the existing
32-column GEMM adapter: same FLOP/X/HX/residuals, whole-row direct output on
MKL, no new allocation buffer/FFI. At representative rows23871, outputs132:
internal-call model188 ->10 and scratch write/read model100.831MB ->0 per
pair. Hypothesis: save10–25% of serial non-residual rotation (0.54–1.34s,
about1.2–3.0% of prior Si64 wall); MPI distributed effect zero. These are
conditional estimates, not Mac or dev08 candidate performance results.

Mac correctness: Si8/Si32/Si64 saved-base/OFF/ON, skew OFF/ON and Si8 strict
OFF/ON science/density bit-identical, abs(dE)=0Ha. Si32 workers1/3 identical.
Slab both exit3 after100 steps with identical science, not a converged pass.
13 ON solver oracle tests and new scalar-reference/tail/nonfinite/pool test
pass. Fast gate616 passed/3 known TMPDIR cache failures/160 ignored,
build/fmt/clippy PASS. No Mac timing conclusion. Separate dev08 request keeps
UPPER_BORDER off, requires route/call/copy counter movement with unchanged
FLOP/work, A/A floor then crossed pairs, abs(mean)/SE>=3. Defaults unchanged,
no adoption. Evidence: `lanes/DAVIDSON/rotation-{correctness.json,validation.md}`.

## DAVIDSON-UPPER-BORDER-REJECTED — orchestrator decision (2026-11-20)

| Attempt | dev08 evidence | Decision / lesson |
|---|---|---|
| DAVIDSON-UPPER-BORDER-REJECTED | Si64 MPI8, Stage 0, timed metric `DAVBENCH kinds[0].solve_ms_median`; A/A first, three crossed pairs. A/A mean -0.000216 s, floor 0.0029167 s, abs(mean)/SE 0.14. A/B differences +0.005178 / +0.003617 / +0.008241 s, mean +0.005679 s, abs(mean)/SE 4.18, about +0.57%, 3/3 slower. Separate OFF/ON counter pair: predicted/executed tile FLOP equal on 8/8 ranks; saving 3.494e9–3.499e9 FLOP/rank; 183 products both arms, tiles 0 -> 405/rank; all-rank traced FLOP 170.792e9 -> 142.821e9 (-16.38%). Reported gemm_dispatch GFLOP 69.247 -> 65.751; calls 1746 -> 1968 (+12.7%). E=-253.37816772116847 Ha both arms, abs(dE)=0 Ha, science bit-identical, SCF 11, H 3504, verified. | **REJECTED: resolved slowdown.** Default remains OFF; implementation retained, no resubmission of narrower upper tiles. Relayed DGEMM rates 23871x264x132=77.14 versus 23871x132x132=63.96 GFLOP/s (-17%) explain why smaller products can lose efficiency despite fewer FLOPs; additional calls compound the loss. Eliminate copies/calls while retaining product dimensions/work instead of assuming every removed FLOP saves time. |

Source: user/orchestrator relay; timed candidate commit was not supplied in
the decision (do not invent one). These counter-pair and timed-screen records
have different scopes and are not added. The previous “queued” status is
superseded by this decision. No Mac timing used.

Alternative selected: already implemented, uncommitted default-off
DAVIDSON_ROTATE_COLUMNS. Retain both PW rows and full inner basis dimension;
partition independent output columns using the existing adapter. It removes
no arithmetic and instead reduces row-scratch calls/copies (188 -> 10 modeled
MKL calls at the serial representative shape). It is not another upper-border
arithmetic-pruning variant. Reuse the unchanged implementation's Mac evidence
in `lanes/DAVIDSON/rotation-correctness.json`; only documentation and the
per-map EVD parser are refined in this rejection-response turn. Independent
serial A/B request is now active; UPPER_BORDER stays OFF in every arm.

## LANE-FFT Phase 2: M2 kept-plane copy fusion (2026-11-20)

Implemented default-off `BAANDI_FFT_KEPT_FUSE`, awaiting dev08 A/B, not adopted.
Relayed bad02ee si64 MPI8 profile: rank-0 lane 2.2440 s, including 0.1123 s
kept_copy, 0.1256 multiply, 0.7184 permutation/zeroing. This motivates removing
the kept copy followed immediately by a slab read/multiply/write: 72 -> 40
logical bytes/kept point, unchanged scalar arithmetic. Conditional prediction
0.03–0.10 s/run (0.4–1.2% of relayed 8.4265 s wall); **not a measurement**.
Serial already fuses this operation and is unchanged. Alternative middle-pack
locality and unresolved padding/wait questions are in
[diagnosis](lanes/FFT/diagnosis.md); no revival of generic MKL FFT, async
alltoall, lane threading or many-band 3D FFT. Relayed rank-max versus rank-0
and wall-ratio inconsistencies are explicitly recorded, not silently combined.
Mac correctness only: five support/grid fixtures, ranks 1/2/3/8 and pairs
1/3/8 bit-identical; real Si8 simulated MPI2/3/8 OFF/ON dE 0 Ha, density FNV
and work identical. Real-UPF warm fixed-H harness actual JSON: each rank's
fine kept_copy 288 -> 0 calls, 11943936 -> 0 bytes; every other stage's
counts unchanged. CLI serial control scientific fields/density bits identical.
Fast gate: 358 pass, 3 known TMPDIR cache failures, 87 ignored; build/fmt/clippy
pass. See [validation](lanes/FFT/phase2-validation.md),
[dev08 request](lanes/FFT/dev08-request.md): A/A floor first, three crossed
A/B pairs, fixed-H si32 then SCF, expansion si64/si32/si16 MPI8; resolved gain
requires negative paired mean with |mean|/SE >= 3 and floor/interval checks.
No performance conclusion from the Mac and no default change.
## DAVIDSON-ROTATE-COLUMNS-REJECTED — orchestrator decision (2026-11-20)

| Attempt | dev08 evidence | Decision / lesson |
|---|---|---|
| DAVIDSON-ROTATE-COLUMNS-REJECTED | One core, A/A floor then three crossed pairs, all positive: Si64 +6.10%, abs(mean)/SE=19.20, 4.74x floor; Si32 +4.19%, 30.46, 29x floor. Si64 pairs109 unchanged, FLOP207.00e9 unchanged, column_pairs0->109, chosen_calls_model16954->720, extra_copy_bytes_model5.59e9->0. Dispatch rotate -207.0018 GFLOP, other +207.0018, total580.791 unchanged, MKL. Si32 calls6138->382, copy1.18e9->0. Si32 bit-identical; Si64 abs(dE)=1.705e-13 Ha. | **REJECTED: resolved losses on both cells.** Retain default OFF. Reducing calls/copies changed GEMM efficiency adversely; do not repeat GEMM-shape changes. |

Relayed single-thread MKL core0: default256x132x132=108.52 GFLOP/s,
82.21 us/call; candidate23871x132x32=53.86 GFLOP/s,3744 us/call.
Pair model188*82.21us=15.46ms versus10*3744us=37.44ms. Removed
scratch51.3MB/pair (~4ms at12.97GB/s, relay's copy-volume convention) is
smaller than22ms shape loss. Height128/256/512/1024/2048/4096/8192/23871
pairms15.78/15.50/17.46/20.03/40.16/40.01/39.93/25.21; the2048–8192
L2 cliff is about42.7GFLOP/s. Column widths32/64/132 give40.49/30.69/25.11ms.
None beats existing256-row15.50ms. Source is orchestrator relay, not Mac.
This supersedes the previous active ROTATE_COLUMNS request.

## DAVIDSON-PACKED-REDUCE — alternative pending dev08 (2026-11-20)

Default-off DAVIDSON_PACKED_REDUCE compacts only the consumed upper-triangle
entries and norm tail after unchanged distributed-lean GEMMs. Same collective
count and per-element rank summation order, no new allocations, no tolerance
change. Representative264-order132-column payload280896->209616 B
(-25.37%); exact trajectory saving is returned by packed-counter.py. Both
rejected switches remain OFF. Serial path unchanged. No performance adoption
or Mac performance claim; see packed-validation.md and dev08-request.md.

## LANE-GLUE Phase 0 (2026-11-20): PROFILE_REQUEST

Instrumentation only, no optimization/default change or Mac performance claim.
`GLUE_PROFILE=on` plus `PHASE_PROFILE=on` enables 23 grid/Pulay spans and 14
work/size counters, rank-local in `run.counters`; clocks obey STAGE_TIMERS.
Existing fine-map/phase/FFT data supplies density and coarse/verification
envelopes. Counters are cumulative, subspans inclusive, not individual traces.
Mac correctness: 30 runs over six inputs (8..68 bands plus synthetic coupled
s projectors), 24 comparisons with pristine baseline, max abs(dE) **0 Ha**,
density bytes/scientific history/solver work identical. Four inputs converge
and verify; elongated/slab preserve the baseline's iteration exhaustion.
Simulated grid ranks 1/3/8 and seeded/sliding/reset plane Pulay ranks 1..7
preserve bits with instrumentation off/counts/timed. Final named fast gate:
627 passed / 3 known in-repo TMPDIR cache failures / 161 ignored;
build/fmt/clippy PASS. Parser contracts and five real Mac-file smoke cases
pass. README checker retains its unrelated stale-phrase failure (8/9 pass).
No speedup or largest-term diagnosis inferred; request si64/si32/si16 MPI8
and si32/si64 one-core profiles through master, then stop. Evidence and
limitations: [phase0-report.md](lanes/GLUE/phase0-report.md);
exact jobs/counter gates: [profile-request.md](lanes/GLUE/profile-request.md).
## LANE-LAUNCHIO Phase 0 (2026-11-20): PROFILE_REQUEST

Base `61822bca3b470df1212ee47b47c7eb4e786e22e3`; instrumentation only, no candidate
optimization or default change. Reuses F-LAUNCH-IO's existing PHASE_PROFILE-gated
lifecycle clocks (including real MPI init/finalize on dev08). Adds input read/parse,
SCF boundary, density encoding/SHA and history conversion clocks; per-species
second-pass UPF read/parse/hash seconds and byte/projector sizes; rank-local
existing between/materialize bucket clocks and map counts. No added collective
or numerical-loop I/O. Preflight internals, allocation counts and per-map glue
call counts remain unmeasured. No rejected handoff or hashing idea is repeated.

Mac correctness only: **20 executions**, four cells (Gamma, skew, 3% anisotropic,
two weighted k points), base plus default/off/environment-on/input-on: **0 Ha**
energy difference, exact density/metadata and scientific-history bits, all
converged and verified. Eight parser smoke checks and five parser tests pass.
Fast CLI gate: **build/fmt/clippy PASS; 132 passed, 3 known in-repo TMPDIR cache
failures, 40 ignored**, exit **1**. Initial sccache sandbox failure was avoided
with CHECK_SCCACHE=0; no tolerance/test changes. A 20%-elongated preliminary cell
failed to converge in the base binary (exit 3), not counted as equivalence.
No Mac timing is evidence. Exact hashes, commands and limitations:
`docs/performance/lanes/LAUNCHIO/validation.md`, `equivalence.json`.

Next: master-only dev08 jobs for Si64/32/16 MPI8 and Si32/64 one-core, five crossed
instrumentation off/on pairs after discarded warmups, rank resource sidecars,
complete raw output and strict parser checks. See
`docs/performance/lanes/LAUNCHIO/profile-request.md` and
`scripts/lane_launchio_parse.py`. This measures attribution/perturbation only;
optimization design and adoption wait for Phase 2. Lifecycle status
**PROFILE_REQUEST**; changes uncommitted.
## DAVIDSON-PACKED-REDUCE-NEUTRAL — orchestrator decision

| Attempt | Relayed dev08 evidence | Decision / lesson |
|---|---|---|
| DAVIDSON-PACKED-REDUCE-NEUTRAL | Si64 MPI8, source61822bc. Mechanism exact on8/8 ranks: calls183 unchanged; sent_bytes282797568 ->203221312, equal to ON prediction; saving79576256 B (-28.1%). A/A mean -0.004163 s, floor0.004163 s, abs(mean)/SE1.60. A/B mean -0.000172 s, abs(mean)/SE0.05,0.04x floor. Byte totals retain the relay's scope; do not multiply them by8. | **NEUTRAL, not adopted; not a resolved loss.** Default OFF retained. A substantial payload reduction produced no resolved gain on this machine; no further volume-only candidate without new evidence. |

## DAVIDSON-LANE-CLOSE — recommendation after three closed candidates

UPPER_BORDER: resolved+0.57% Si64MPI8 loss. ROTATE_COLUMNS: resolved+6.10%
Si64/+4.19%Si32 serial losses. PACKED_REDUCE: neutral despite28.1% fewer
communicated bytes. All switches stay default OFF, no numerical code changed
in this closing review. Recommend closing the current lane, not declaring the
algebra gap solved or the implementation universally optimal.

Latest relayed Si64MPI8 attribution: subspace algebra gap to same-model QE
+1.91s; BAANDI reduced_gram1.002s, reduced_eigen0.988s, rotate0.963s,
kept_update0.573s. These are BAANDI component times, not four independent
QE differences. Fine iterations58 versus QE whole-SCF58; relayed cost99ms
versus69ms/iteration, with unequal coarse/verification scopes explicitly
retained. Outside-envelope gap+0.77s; listed BAANDI components total0.431s,
not0.77s. No unsupported assignment of the0.339s numerical difference.

A0.3s candidate needs30.4% of EVD,31.2% of rotation,52.4% of kept_update,
or69.6% of the listed outside components. Existing dsyevd scratch is reused;
selective dsyevr was slower at these orders (MKL-EVD); REL-SLICES failed its
speed/orthogonality gate. No duplicate unchanged reduced problem is identified.
Kept_update includes delta-potential H application, not just copy overhead.
The code and measured copy/shape experiments do not establish an eligible
non-GEMM saving of0.3s. Full reasoning and reopening criteria: diagnosis.md,
"Closing assessment". No new switch, no new dev08 queue, no Mac timings.
## LANE-NONLOCAL Phase 2 — recommend closing (2026-11-20)

Analysis only, no implementation or default change. Relayed by the user:
dev08 `bad02ee`, Si64 np8 lane 0.7825 s / wall 9.6702 s (8.09%); repeats
0.7825/0.8125/0.7486 s, spread 8.5%, timed/counts 24 work keys identical.
Projection 10.43 GFLOP / 0.3054 s = 34.15 GFLOP/s; back 7.48 GFLOP /
0.2681 s = 27.90 GFLOP/s. At supplied 64–108 GFLOP/s reference, combined
theoretical headroom is 0.294–0.408 s, not a measured removable overhead.
Logical bytes/STREAM (12.97 GB/s) gives 0.0524/0.0413 s reference times.
Existing real-Gamma 256-panel DGEMMs have no per-call complex conversion.

Assessed same-shape direct-output/add fusion (about 0.050 s optimistic
bandwidth-scale estimate at modeled P=256), coupling elimination (entire
stage only 0.0474 s), and persistent packing (unattributed, no supported
>=0.15 s expectation). Reduction is 0.1178 s with 66.6% asymmetry; the user's
PACKED_REDUCE relay reports no wall gain despite 28% less payload. Do not
repeat A43/A58b traversal, tile narrowing or column groups. No candidate
meets this round's evidence-based expected >=0.15 s screen. Recommend closing,
not claiming optimality. No Mac timing, numerical code edit or new dev08 job.
Detailed assumptions, ledger/default reconciliation and reopening conditions:
[diagnose.md](lanes/NONLOCAL/diagnose.md),
[closing-note.md](lanes/NONLOCAL/closing-note.md). Raw dev08 files/hashes were
not attached; the relay was not independently reproduced locally.

## FFT-KEPT-FUSE-NEUTRAL (orchestrator relay, 2026-11-20)

CLOSED NEUTRAL; FFT_KEPT_FUSE stays default off. Extended si32 MPI8 crossed
run, socket 0 cores 0–7: A/A 6 pairs, A/B 8 pairs, first pair discarded as
reported; 32 runs identical in eigen_fnv, h_columns and iterations. Relayed
A/A paired mean +0.003227, SE 0.004108; A/B +0.008744, SE 0.005753,
|mean|/SE 1.52, 2.13 times floor, intervals overlap; 5/8 A/B pairs slower.
Counter gate exact: kept_copy 25344 -> 0 calls, 2102132736 -> 0 bytes,
multiply unchanged. No resolved wall gain; estimated 96 pairs to resolve
-0.5% is not worth pursuing. Do not repeat the fusion or infer speedup from
its removed traffic. Pair counts/discard convention are retained as relayed;
raw samples/binary hash were not supplied here.

Independent NONLOCAL bad02ee si64 MPI8 relay: lane 0.7825 s, 8.09% of wall;
reduce (zero flops) max/min 0.1178/0.0394 s, 66.6% spread, versus compute
15–19%. PACKED_REDUCE removes 28.1% payload without wall change. Together
with FFT evidence this prioritizes arrival-skew attribution, not another
volume-only change. It does not prove all collective volume irrelevant at all
sizes, or identify the last rank from aggregate durations. Follow-up:
[wait/arrival diagnostic design](lanes/FFT/wait-diagnostic.md) and
[profile request](lanes/FFT/wait-profile-request.md). No new optimization:
0.3 s si64 MPI8 saving from a specific permutation-pass removal is not supported
by current counters (it would require about 42% of the 0.7184 s group).

## LANE-FFT arrival trace implementation (2026-11-20)

Instrumentation-only scope approved after FFT-KEPT-FUSE-NEUTRAL. Added default-off
COLLECTIVE_TRACE with PHASE_PROFILE: bounded per-rank map/collective entry/exit,
previous completed stage/work/shape, FFT/NONLOCAL/solver reduction sites.
No new application collective or arithmetic change; two diagnostic gathers
only after SCF. Count-only reads no trace clocks; OFF allocates no trace buffer.
Same-host clock anchors report uncertainty/drift; parser identifies unresolved
last-rank ties and separates previous-release skew from local interval skew.
Exposure is not exact idle wait; parent stage durations remain inclusive.

Mac correctness: OFF/count/timed Si8 simulated MPI2/3/8 energy bits/density FNV/
work identical (dE 0 Ha); serial CLI scientific fields and density bytes identical,
verified. Actual CLI profile serializer on warm MPI3: 408 events/rank, zero
count/timed signature differences, 408 parsed fine rows. Five parser tests and
buffer/transport/all-five-operation unit tests pass. Final fast gate: 705 pass,
3 known TMPDIR cache failures, 178 ignored; build/fmt/clippy pass. No Mac speed
claim, no optimization. See [request](lanes/FFT/wait-profile-request.md),
[validation](lanes/FFT/arrival-validation.md), [evidence](lanes/FFT/arrival-equivalence.json).

## 2026-11-20 — LANE-VERIFY Phase 2: closed without a qualifying lever

**CLOSED — no implementation.** Six dev08 Si64 MPI8 instrumented runs at
`bad02ee`, as relayed by the user: median external wall **8.5466 s**, median
verification sum **1.0475 s (12.26%)**, two attempts and 11 ordinary maps in
all runs. Attempt 0 fails (6 Davidson iterations, 380 H columns); attempt 1
passes (5 iterations, 388 H columns). No Mac performance measurements were used.

Default HALF already carries failed verification orbitals into the next map;
M2_REUSE can update retained H X for a changed potential. Actual cache-hit counts
were not relayed. The next Hamiltonian still needs a new solve and convergence
test; the full Davidson search space is not retained. Reimplementing the adopted
reuse is not a new saving. Energy's reported 0.0499 s and the combined non-solve
phase envelope of 0.1676 s do not establish a >=0.2 s lever (independent phase
maxima are not a coherent wall bound). Full search-space recycling is an
unqualified, separately owned solver hypothesis, not a predicted gain.

The judgment rule and attempt timing are unchanged. No source edits, switch,
A/B request or adoption claim. One-core profile remains pending; no MPI-to-one-core
extrapolation. Details and correctness-check scope:
[diagnosis](lanes/VERIFY/diagnosis.md), [closing note](lanes/VERIFY/closing-note.md).

Phase 2 Mac checks: `scripts/check.sh fast` (CHECK_BASE=HEAD, sccache off,
in-repository scratch) exit **1**: build/fmt/clippy pass; **617 passed /
3 known CLI cache TMPDIR failures / 157 ignored**. `git diff --check` exit **0**.
Documentation only; no new numerical equivalence or performance claim.

## LANE-GLUE Phase 2 — recommend closing (2026-11-20)

Analysis only, from user-relayed dev08 `e947a35` (instrumented, socket 0).
Si64 GLUE 1.1690 s MPI8 / 8.2797 s one core; density 0.7331 / 6.3032 s.
No candidate substantiates the task's >=0.2 s MPI8 or >=0.8 s serial expected
saving. Pair FFTs and MPI fused/owned density already adopted; no rank-wide
density sum remains on the eligible owned path. Serial MKL density still uses
a full-grid adapter and split contribution arrays: fixture model (72^3,
128 occupied bands, 15 maps) has 960 paired inverse transforms, 0.420213 s/map
normalization; direct accumulation could remove 11.466 GB logical traffic,
but >=0.8 s requires <=14.333 GB/s effective removed-pass throughput, absent
from the relay. This is a break-even condition, not a gain prediction.
`m2.xc` includes potential/audit FFTs and collectives; 37.4% asymmetry cannot
be attributed to pointwise XC without the pure-XC counters. ENERGY_IDENTITY
previously lost; dropping fresh verification density would invalidate the
residual. Five candidates and their limits are recorded in
[diagnosis.md](lanes/GLUE/diagnosis.md); [closing-note.md](lanes/GLUE/closing-note.md)
records reopening evidence. No switch, optimization, A/B job, default or Rust
source change. No Mac profiling. Named engine gate exit 0: build/test/fmt/clippy
PASS, 260 passed, 0 failed, 74 ignored. Changes uncommitted; stop here.
The task's screening threshold does not supersede the standing adoption rule.

### LANE-GLUE serial density fusion — 2026-11-20, default-off, dev08 request

Supersedes the Phase-2 closing recommendation above. **New premise:** the
orchestrator relayed one-core STREAM triad 12.97 GB/s (DGEMM 64–108 GFLOP/s)
and explicitly authorized implementing the serial density/FFT candidate.
11,466,178,560 conservative logical bytes / 12.97e9 predicts **0.884054 s**
at Si64 one core, about 2.17% of 40.7 s; this is a hypothesis, not measured gain.
No rejected GEMM/ENERGY_IDENTITY experiment is retried.

Implemented `BAANDI_SERIAL_DENSITY_FUSED=on` (default off): ordered pair sink
avoids per-band density-grid traffic, keeps each band addition and FFT/psi
arithmetic unchanged. Standard, one-core, no registered communicator and
non-MPI build only; strict/complex/multi-worker/MPI paths unchanged. Arrays
remain allocated for fallback; no resident-memory claim. Counters under
GLUE_PROFILE+PHASE_PROFILE charge old 32 versus new 16 bytes/paired band-point.
Mac final baseline/OFF/ON scientific bits match: 12 converged+verified standard
runs (Si8, skew, synthetic S2, Si32), strict skew 3 arms, thread fallback 2 arms;
6 elongated/slab arms retain baseline exhaustion and identical histories.
All |dE|=0 Ha. Si32 accounted bytes 5,249,359,872 → 2,624,679,936; exact predicted
saving 2,624,679,936 bytes over actual coarse/fine/verification calls.

Final fast gate: build/fmt/clippy PASS, 728 passed / known 3 Mac TMPDIR
pseudo-cache failures / 179 ignored, exit 1. Switch-on simulated M2 suite
3/3 passes, dispatch tests exclude ranks 1/2/3/8. Generic full-grid and pure
FFT sink tests pass; actual MKL/real MPI remain dev08 validation. README checker
retains its unrelated structural wording failure (8/9 pass). No Mac timing,
QE run, default adoption or commit. Energy/XC still have no demonstrated
acceptance-preserving >0.2 s MPI8 lever; m2.xc is inclusive, energy has no
matched QE timer. **Lesson:** STREAM changes candidate selection, but logical
bytes are not measured DRAM bytes; only crossed same-model wall decides adoption.

Evidence: [fusion-validation.md](lanes/GLUE/fusion-validation.md),
[fusion-mac-equivalence.json](lanes/GLUE/fusion-mac-equivalence.json),
[dev08-request.md](lanes/GLUE/dev08-request.md). A/A first, crossed six A/B pairs
per case, |mean|/SE>=3; any resolved gain counts, 0.8 s is not an adoption floor.
## 2026-11-20 — LANE-VERIFY reopened: first-attempt contraction scheduling

**Experimental, default unchanged; no adoption recommendation.** The
orchestrator now permits timing changes while preserving ADR 0043 judgment.
Source correction: the adopted scheduler predicts only retries, not the first
attempt. Its np8 prediction count of one is not an initial false-pass forecast.
Relayed Si64 one-core median verification **4.8284 s / 41.5704 s (11.62%)**;
MPI8 **1.0475 s / 8.5466 s (12.26%)**. No Mac profiling.

Added opt-in `VERIFY_SCHEDULE=first` and `first-trace` values (registry count
unchanged, default `on`). The first uses a threshold-relative geometric residual
amplitude tail, the slower Coulomb/L2-squared contraction, and a two-map bound
from the first crossing; it never marks a skip as a pass. First-trace preserves
baseline numerical decisions; both retain adopted retry scheduling. No fitted
material constants or certificate/tolerance changes.

Final 8 systems x 3 seeds x 2 accuracy classes: **48/48** certified and within
class; precursor default/parent **48** and strict-path **24** identities;
final trace/parent **48** scientific/density identities plus **4** fresh direct
default/strict identities. Only the opt-in predictor changed between variants. First baseline failures **29**, predicted failure
skips **28**, wrong skips **16**. Actual attempts **81->51**, maps **419->444**,
total H **48,676->47,588**, failed-attempt H **1,814->255**, failed Davidson
iterations **151->16**. These sums are work counts, not wall savings. Maximum
energy error against internal strict references: standard **1.4630785472036223e-10
Ha**, tighter class **3.5704772471945034e-12 Ha**.

**Negative target evidence:** supplemental Si64 seed 601 retained two attempts.
The lag-only precursor moved failure from map 11 to 12 (total H 6,644->6,622).
The final tail model waits to map 13, which still fails, then passes at map 14:
maps **13->14**, total H **6,644->6,661**, final strict-reference error
**2.3021584638627246e-12 Ha**. Omitting the failed HALF product changes the
subsequent state; the baseline's later maps cannot simply be assumed identical.
Neither model establishes removal of the relayed failed-attempt cost. Retain
negative evidence and the switch for orchestrator review, not a default flip.

Release scheduler tests: **6 passed**, including simulated 2/3-rank identity.
Final fast gate exit **1**: build/fmt/clippy pass, **473 passed / 3 known CLI
in-repository TMPDIR failures / 137 ignored**. Initial test-type/clippy and
strict-JSON assumptions were fixed; final checks include the fixes.
[Design and historical features](lanes/VERIFY/first-schedule.md),
[final validation](lanes/VERIFY/first-schedule-validation.md),
[conditional dev08 A/A and A/B request](lanes/VERIFY/dev08-request.md).
Changes uncommitted; no remote run or speed claim.

## LANE-LAUNCHIO Phase 2 (2026-11-20): metadata-keyed executable digest cache

**Default OFF; DEV08_AB_REQUEST; no speed/adoption claim.** User-relayed dev08
profile e947a35, median six instrumented on-runs: fresh binary SHA-256 costs
0.10446/0.10474/0.10456 s on Si64/32/16 MPI8, 0.09387/0.09424 s on Si64/32 one
core. Source inspection proves density_bin_write is already-replicated bytes
written by all ranks (P*D output payload), not a gather inside that timer. The
terminal density allgatherv happens before CLI serialization. Raw byte sizes and
bandwidth/resource data were not in the relay; floors are explicit formulas,
not invented device rates. See `lanes/LAUNCHIO/diagnosis.md` for alternatives.

New premise versus F-LAUNCH-IO's earlier caching deferral: real dev08 fixed-cost
evidence and user authorization to cache while retaining provenance. Implemented
`BAANDI_BINARY_HASH_CACHE=on` behind a default-off registry switch. Unix sidecar
keys device/inode/length/mtime+ns/ctime+ns/owner/mode; validates opened-file/path
identity and bounded records, checks record corruption, atomically publishes a
mode-0600 sidecar. All misses/publication failures still produce a fresh digest.
No omitted provenance, numerical work, output suppression or new collectives.
Trusted coherent build-directory metadata required; malicious same-UID sidecar
forgery is outside the contract, and master independent pre/post SHA is mandatory.
Expected warm saving 0.05–0.10 s is a planning estimate, not measured here.

Mac correctness only, base 0413485: **32 executions across four inputs**, base
and unset/off/cold/warm/malformed/input-on/unavailable arms, direct density file
byte identity and entire normalized results JSON identity; **0 Ha** energy delta,
all converged/verified, every binary digest matches independent Python hashlib.
Eight warm/input hits report zero hashed executable bytes. Five cache unit tests
cover invalidation, corruption, unsuitable paths and eight concurrent writers.
Fast gate **exit 1**: build/fmt/clippy pass, **222 passed / 3 known in-repo TMPDIR
pseudo-cache failures / 62 ignored**, 19 binaries. No tolerance/test weakening.
README generated-block/figure freshness passes. Exact commands, failures, hashes
and exclusions: `lanes/LAUNCHIO/phase2-validation.md` and `phase2-equivalence.json`.

Master request `lanes/LAUNCHIO/dev08-request.md`: Si16/32/64 MPI8, Si64 one-core;
retain cold cost, six A/A pairs first, twelve crossed warm A/B pairs with phase
clocks off, then three separate diagnostic pairs. Require warm hits/zero hash
payload, independent provenance, exact outputs, diagnostic hash-span attribution,
and positive resolved whole-run gain |mean|/SE >=3 without regressions beyond
noise. The >=75% hash-span drop is a diagnostic forecast, not an adoption
percentage gate; any resolved whole-run gain qualifies. Default remains off;
changes uncommitted until review and dev08 approval.
## FFT-ARRIVAL-OWNERSHIP-CLOSED (2026-11-20)

Recommendation: close FFT optimization on available evidence; no new switch or
ownership change. Relayed 45340f6 si64 MPI8 trace: 776 matched events/rank,
684 ns clock error, rank-summed exposure 2.7990 s (0.349875 s/rank), median/p90
arrival skew 504.3/2409.6 us. Rank 3 last 220 times but reconstructed default
mirror ownership is identical to ranks 2/4: 126 sticks, 2986 ORC rows, nine
transpose planes. Across ranks rows differ by four, sticks by two. No supported
0.15 s balancing gain. solver.append is a last-completed-stage label and can
precede untimed local border work, not causal attribution. Jitter/NUMA remains
unproved; raw matched interval/stage data and queued si32 are missing. Mean rank
exposure is not a critical-path wall bound. See [closing note](lanes/FFT/arrival-closing-note.md)
and reproducible [ownership audit](lanes/FFT/ownership_audit.py). No Mac profiling,
no runtime changes; FFT_KEPT_FUSE remains OFF/NEUTRAL.

### LAUNCHIO-BINARY-HASH-CACHE-DEFAULT-ON (2026-11-20, dev08 via master; ledger by orchestrator)
BAANDI_BINARY_HASH_CACHE (cache of the executable SHA-256 in a sidecar keyed by size/mtime/inode; run.binary.sha256 is still produced; invalidated when the file changes; hit reads 0 bytes) default flipped on. dev08 crossed pairs, A/A floor first: si16 np8 -0.104577 s (|mean|/SE 19.05, 12/12 pairs positive, -10.68 % of 0.9788 s median A), si32 np8 -0.095911 s (20.24, 12/12, -4.15 %); si64 np8 and one core pending at adoption time (the gain is a fixed 0.1 s, cannot lose). Mechanism gate: hit x16 with hashed_bytes 0 and one miss that fills the cache (22,106,160 B), sha256 identical in every arm; energy difference 0, density bytes identical. Predicted 0.1045 s (LAUNCHIO phase 0 binary_sha256) matched. Mac fast gate 492 passed, 0 failed. Side effect to remember: a small cache file is written next to the executable when that directory is writable (no write otherwise; a failed write is only a miss).

### GLUE-SERIAL-DENSITY-FUSED-DEFAULT-ON (2026-11-20, dev08 via master; ledger by orchestrator)
BAANDI_SERIAL_DENSITY_FUSED (standard path, single process, one-core lane: the per-band density accumulation avoids the intermediate density array; band order and FFTs unchanged) default flipped on. dev08 crossed pairs, 6 per cell, A/A floor first, 1 core: si64 -2.0801 s (|mean|/SE 62.44, 6/6 positive, -4.77 % of the same-window A median 43.6417 s), si32 -0.2949 s (119.19, 6/6, -3.73 % of 7.9103 s), si16 -0.085344 s (7.83, 6/6, -3.65 % of 2.3397 s), si-skew -0.001713 s (3.43, 5/6, -1.03 %, resolution edge, A/A ratio 2.78); si8 first attempt aborted by the A/A guard (window suspect), retest queued. No resolved loss. Energies bit identical in every cell (|dE| 0), SCF/H/verified identical. Counters exact: logical bytes halved (si32 5.25 GB -> 2.62 GB, si64 26.8 GB -> 13.4 GB), eligible calls all fused. Prediction for si64 0.884 s, measured 2.080 s (2.35x). The gain applies to the one-core lane only (no effect with MPI ranks). Mac full gate 798 passed, 0 failed. Lesson: the winning lane candidate removed intermediate-array memory traffic; those changing GEMM shapes lost, those reducing communication volume were neutral.


## G-FIXTURE-SI128 — larger-cell measurement surface (2026-11-22)

Status: fixtures and dev08 plan added; **no optimization or speed claim**.
The new premise is scale beyond the existing 64-atom catalog, not a retry of
any rejected numerical/kernel idea. Added only Si128 (2x2x4) and Si216 (3x3x3),
30 Ry Gamma, long tier, 7200 s per-run cap, same Si64 field set and potential.
No code, defaults, tolerances or existing fixtures changed.

Mac engine preflight/grid exits 0: Si128 512 electrons / 260 bands / 47,851 PW /
72x72x144 FFT; Si216 864 electrons / 436 bands / 80,797 PW / 108^3 FFT.
Regenerated Si64 coordinates match exactly as a set. Every generated atom has
four periodic nearest neighbors at 4.443235504 Bohr (the preserved lattice does
not imply the ticket's 4.6 Bohr). No overlapping/missing atoms.
Si512 is documented only: 191,735 PW, 144^3 FFT, real orbital/projector blocks
1.469/2.926 GiB globally; a specified retained-array scenario is 70.455 GiB,
not measured RSS. No runnable Si512 catalog entry; 64 GB feasibility unverified.

Validation: `scripts/check.sh full` with worktree-local TMPDIR and sccache off
exit **1**, build/fmt/clippy PASS; 795 tests passed, 3 known cache-policy tests
failed because scratch is inside the repository, 181 ignored (87 binaries).
Benchmark suite Python tests: **17 passed**, exit 0. Geometry + four sizes'
CLI validate/grid: exit 0. Mac performed no 128+ atom SCF; no dev08/QE runs.
Lesson: basis preflight and block payloads establish workload size, not peak
memory or time. Measure all-rank RSS and actual convergence before admitting
Si512 or inferring a scaling law. Full size tables, formula assumptions,
reproducer and same-model np8/np16 A/A-first requests:
[fixtures-si128.md](lanes/fixtures-si128.md). Changes left uncommitted.
## G-FILTER-REPLAY — resident block Chebyshev, stage 1 (2026-11-23)

**STOP integration; replay artifact only.** New premise versus the older
Chebyshev warning: price a resident fixed-width block against the relayed Si64
MPI8 2.953 s Gram/eigen/rotation budget (183 EVDs), rather than reopen PR-CG
TRACK-A0, REL-SLICES, Advance, restart-keep or rejected GEMM-shape changes.
No production caller/default/dependency changed; existing paired Gamma H and
wide rotation kernels are used. Degree grid 4/8/12/16 plus spectrum/tolerance
auto, 24 counted bound H columns, scaled recurrence, existing CGS2 QR, one RR
per sweep plus setup RR, exact captured stopping profile, recorded fallback to
Davidson from the original snapshot. All failed work is charged.

Mac counts (not timing acceptance): 15 snapshots, 75 arms on Si8/16/32/64/skew;
actual fine map counts 11/10/10/13/8. Serial captured solves reproduce H counts
and eigenvalue FNV. Si64 early/mid/late Davidson H = 593/289/235; d4 =
5436/1608/1608, d8 = 4380/3192/2136. Frequencies 5/4/4 scaled to 11 maps give
33883.38 H (d4), 36564 (d8), 53094.46 (d12), 70518.46 (d16), 56208.31 (auto,
including fallback), versus the relayed **MPI8** fine baseline 3504 columns.
This is a count projection across machines, never a mixed-model timing ratio.

Relayed prices: 2.782 s / 6025 H columns, DGEMM 64–108 GFLOP/s,
STREAM 12.97 GB/s/core. Over-crediting all 2.953 s to the fine replacement and
allowing ideal cached QR, faster nonlocal work and maximal baseline EVD orders
still gives net savings d4/d8/d12/d16/auto =
**-9.8693/-10.8026/-17.6246/-24.8206/-20.5793 s**. Conditional pessimistic d4
is -35.5381 s. These are cost estimates, not dev08 measurements; omitted opaque
allocation/accuracy-repair work is nonnegative and cannot rescue the verdict.
Stop under the issue's nonpositive-optimistic exit rule. No degree or threshold
was retuned to Si to evade rejection.

Accuracy limits are explicit: energy is a fixed-H occupied band sum, not total
SCF energy. Early/mid/late standard 5e-9 Ha gates pass 5/25, 20/25, 25/25 arms;
strict 1e-10 Ha gates pass 2/25, 4/25, 22/25. Captured early Davidson also fails
these absolute frozen-energy gates; matching its profile does not make early
maps certified. Principal angles are degeneracy-safe; density errors/units and
all failures are retained. No full candidate SCF or total-energy acceptance.
Simulated np2/4 exercise rank-local payloads/fallback; exact engine replay needs
M2_LOCAL_FIELDS=off (capture disables owned fields). Reconstructed MPI operator
rounding differs (HX <=3.89e-16 Ha); its baseline is marked non-exact and stale
HX is dropped. No claim of exact MPI filter identity or rigorous spectral bound.

Implementation delivery is **blocked at the bench ownership boundary**: the
ticket names nonexistent `dft-solvers/examples/davidson_bench.rs`, while the real
harness is in dft-cli. Correction requested, no answer yet; the actual CLI file
was not edited. The tested patch and scratch-only builder are reviewable in
`lanes/FILTER/`; they do not claim an installed CLI mode. `gh issue view 24`
failed (network), so issue scope is taken only from the supplied ticket.

Tests: solver 4/4; prepared bench principal-angle oracle 1/1 and clippy pass;
profile off/on has identical counts/accuracy/residuals/sweeps. Five default
SCFs retain identical energy values and density SHA-256, all verified. Full
check exit 1: build/fmt/clippy pass, 799 pass, 3 known in-repo TMPDIR pseudo-cache
failures, 181 ignored. No commit. Details, reproducible tables/model and a
bounded dev08 falsification request: [filtered-subspace.md](../design/filtered-subspace.md),
[summary.json](lanes/FILTER/summary.json), [dev08-request.md](lanes/FILTER/dev08-request.md).

## G-MIXER-DIVERGENCE (2026-11-23; Mac counts/correctness only)

| Attempt | Mechanism and changed premise | Evidence | Disposition / limits |
|---|---|---|---|
| MIX-CONTROL-AUDIT | Five control questions: beta already multiplies R and Delta R in both Pulay implementations; mixing_fraction is only extra positivity damping. Fatal divergence uses 3 consecutive normalized L2 residuals >20x best, not the Coulomb error and not energy jumps. The projected Pulay ratio is read-only. The standard eigen schedule retains its historical minimum when error grows. | Relayed si128 map 4 L2/best = 2.85016, despite 272x Coulomb growth; serial beta 0.2/0.5/0.8 one-secant oracle plus distributed chunks/planes P=1/2/3 added. Exact source anchors in `lanes/MIXER/source-anchors.txt`. | Finding, not a missing-beta fix. No double beta or acceptance change. |
| MIX-CELL-KERKER | New default-off `BAANDI_MIX_CELL_KERKER`: q0^2 = max(0, (2*pi/Omega^(1/3))^2 - gmin^2) from the actual density basis; no species/size threshold or fitted q0. Reuses Kerker-Pulay. Cubic cells are exactly inert. Unlike A31's fixed-q sweep, the premise is the long-wavelength excess relative to an equal-volume cube and a normally convergent 1x1x4 stress geometry. | 9 systems x seeds 601/602/603, f1e-3. 1x1x4 fine maps 26/79/37 -> 8/9/9; total recorded H including coarse/certificate/start 8533/21876/10125 -> 2578/2475/2532. Si8 maps 11/9/9 unchanged; Si16 10/10/10 -> 7/6/7; Si32 10/10/11 -> 9/8/8. Skew/Gaussian-fcc/p3 unchanged. All 21 formerly successful runs still verify; max standard error vs tight strict reference 1.6735910435272672e-9 Ha; strict-C vs strict-A <=6.252776074688882e-13 Ha. | Candidate for dev08, NOT a timing result or adopted default. Old si-elongated/slab remain 100-map failures in all 3 seeds; no oracle there. MPI uses replicated Kerker and loses owned-density mixing; real MPI8 wall must decide. |
| MIX-CELL-COULOMB-COMBINATION | Independently run existing Coulomb guard alone and together with C, without changing STOP_PREDICTOR or verification. | Guard alone on 1x1x4: 34/39/29 maps (baseline 26/79/37); Si32 10/11/12 vs 10/10/11, skew 9/9/9 vs 8/8/8. C+guard is identical to C on active Kerker cases (existing guard does not apply to Kerker) and regresses skew to 9/9/9; Gaussian-fcc 9/8/9 vs 8/9/8. | C alone is the locally supported combination. No revival/default promotion of MIX_COULOMB_GUARD from these data. |
| MIX-CELL-IDENTITY | Pre-change executable vs candidate-off; candidate C threads 1 vs 4; final binary off confirmation. | 27/27 pre/post scientific histories, energy/density and H identical; 27/27 thread pairs identical; final-binary off: 9/9 systems at seed 601 identical. Simulated M2 si16 P2/P3: 8 maps, 36 Davidson iterations, 619 fine H each; energy vs serial <=5.471e-13 Ha, P2 vs P3 differs 1 ulp (not a rank-bit identity claim). | Full gate exit 1: build/fmt/clippy PASS; 798 tests pass, 3 existing pseudo-cache tests fail because writable scratch is inside the repository. README checks: 8/9 pass; existing structural-scope text assertion fails. No unrelated tests or tolerances changed. |
| MIX-RECOVERY-EIGEN-DESIGNS | Recovery from measured standard-error/model mismatch with saved best pair; independently a dimensionless growth-reactive eigen cap. | Algorithm contracts, failure cases and preregistered tests in `docs/design/mixer-divergence.md`. Session wiring is in unowned `session.rs`; scope request pending. | Designs only; no fictitious guard/eigen switches registered or measured. Fixed-point acceptance unchanged. Wrapper diagnostic forwarding also requires unowned `mixing/kerker.rs`. |

Raw traces, input/binary/UPF hashes, command exits, failed inputs and reproducer:
`docs/performance/lanes/MIXER/`. No UPF stored. Remote commands, A/A first,
crossed AB/BA jobs and 600/7200 s caps: `lanes/MIXER/dev08-request.md`.

## G-OWNED-ROWS — owned input packing (2026-11-23)

**Default OFF; dev08 A/B requested; no speed/adoption claim.** Base f526588.
Issue #20's direct kinetic/projector rows and lazy serial real operator are
already implemented/adopted; do not re-credit them. Current audit finds no
new substantiated >=0.15 s scoped lever, but the remaining full complex -> full
ORC -> local input bridge allocates 8GB extra bytes per eligible entry/rank.
`M2_OWNED_INPUT` packs the selected rows directly on the standard M2 workspace
path, preserving the full HalfSphere validation and coefficient arithmetic.
Strict/M1/serial/non-workspace and already-local inputs retain their paths.
No new geometry cache; no ionic or atomic orthogonalization arithmetic change.

Mac correctness: Si8/16/32/si-skew x simulated ranks 1/2/3/8, plus six synthetic
local/nonlocal standard/fallback/strict pairs: **22 A/B pairs, 68 rank records,
0 Ha energy difference**, exact density/scientific history/eigen-solutions;
standard verification passes. Rank one is the unchanged serial negative control.
Calls/allocation-request counters agree exactly with independently reconstructed
G/L_r/B shapes. Active Si32 saves 6,461,088 B/rank/entry. Static Si128 prediction
99,530,080 B/rank/entry; not RSS. Existing idle serial BandDensity grids total
1,558,683,648 B/rank at Si128, a higher-memory follow-up requiring another owner.

Final full gate exit **1**: build/fmt/clippy PASS, 797 passed / 3 known CLI
in-repository scratch-cache failures / 183 ignored. Scoped engine gate exit 0:
265 passed / 0 failed / 77 ignored, all four gates PASS. README regeneration leaves
no diff; its separate test suite has an unrelated stale expected-phrase failure
(8 passed / 1 failed). No tolerances or failing tests weakened. Dev08 request:
Si32/64/16/128 MPI8, A/A first, crossed pairs, exact counters, all-rank RSS and
whole-run |mean|/SE>=3; defaults unchanged. Lesson: audit the current HEAD before
reimplementing historical hotspots; remaining input-boundary memory is distinct
from already-resident warm maps and previously adopted construction savings.
[Audit](lanes/OWNED/audit.md), [validation](lanes/OWNED/validation.md),
[evidence](lanes/OWNED/mac-evidence.json), [dev08](lanes/OWNED/dev08-request.md).

## R-SCF-START — final-cutoff owned initialization (2026-11-23)

**Explicit library entry only; defaults unchanged; dev08 timing/RSS pending.**
Base d97d86430fd943ed6abfd1cfbacdb580846f6cef, uncommitted worktree.
Derived from initialization principles in the own-words note
`docs/research/qe-study/qe-scf-cycle.md`; no QE source read. This is a new
reference-driver start, not a revival of the rejected coarse owned handoff.
It has no coarse/fine transfer or padded coarse-to-fine exchange.

Motivation is the ticket's relayed startup/memory decomposition and the
INIT/LAUNCHIO/OWNED audits. The additional saving is eliminating replicated
serial problem/operator/BandDensity state and full atomic Q/HQ/start blocks;
previously adopted sparse projectors/lazy serial real H are not re-credited.
`M2_OWNED_INPUT` remains default off and is not used by the new entry. The new
function directly builds local projector/ORC rows, owned charge and potential,
uses the existing guarded CholQR2 factorization and distributed H with chunked
atomic HQ/Ritz, and retains final-cutoff state. No randomization; deterministic
kinetic completion. The existing Ewald split/order is reused before allocating
large fields. Density metadata and sparse projector source are explicitly
released before the atomic block peak. Detailed byte/frequency audit and
numerical/unsupported-input contracts: `docs/design/reference-scf.md#start`.

### Measured local correctness (no Mac speed claim)

Stored evidence: `crates/dft-engine/tests/refstart-evidence.json`, pure Rust,
Rust 1.99.1, six release tests, **126 per-rank comparison records**. Si8/16/32
and skew Si8 at final 15 Ha, deterministic 1e-4 fractional displacements;
local/nonlocal synthetic PPs, library LDA-PW, missing-band completion and empty
owners; 1/2/3/8 simulated ranks. Dependence/unsupported-input rejection and
same-seed/different-unused-seed determinism also pass. External Si UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
The reference is existing BAANDI at the same final cutoff, **not QE**.

| Quantity | Maximum absolute discrepancy | Gate |
| --- | ---: | ---: |
| Starting density | 8.326673e-16 electrons/Bohr^3 | 2e-12 |
| First potential | 1.632028e-14 Ha | 2e-11 |
| Gram minus identity | 1.709743e-14, dimensionless | 2e-12 |
| Selected atomic-Ritz span column projection residual | 6.282379e-12, dimensionless | 2e-8 |
| Initial Ritz energies | 9.769963e-15 Ha | 2e-11 |

Zero mode is bit-equal to N_el/Omega; normalized charge meets relative 2e-12.
Eigenvector signs/degenerate rotations are not physical discrepancies; tests
compare the subspace and Ritz spectrum rather than requiring arbitrary gauges.
No cap was relaxed after observing results.

Legacy standard-driver continuation from the rank-3 owned start, imported
through its unchanged warm-checkpoint API: **8/8 SCF steps**, both fixed-point
verifications pass. Legacy energy -0.623797726375388772 Ha, new-start energy
-0.623797726375388439 Ha, |delta| **3.330669e-16 Ha** (internal cap 1e-12 Ha).
This end-to-end check is synthetic, not a new Si or external-code SCF oracle.

### Per-rank bytes: named-array bounds, not RSS

Actual rank maxima at P=8 from the final-cutoff correctness fixtures:

| Case | Owned output block, bytes | New principal envelope, bytes | Legacy principal envelope, bytes |
| --- | ---: | ---: | ---: |
| Si8 | 60,160 | 1,514,144 | 15,547,616 |
| Si16 | 216,576 | 4,176,768 | 51,848,992 |
| Si32 | 810,560 | 13,033,792 | 187,091,936 |
| skew Si8 | 60,320 | 1,522,784 | 15,554,400 |

The output block counter is exactly 8*L_r*B. The envelopes are conservative
named-array overlap formulas, **not measured allocator peaks**. FFT/GEMM/XC/EVD
scratch, communication buffers, metadata, radial caches, allocator and MPI
are excluded; the full-indexed 8F operator potential remains. These numbers
must not be reported as process memory or RSS ratios. Full high-water RSS is
requested on dev08, per rank. The large serial `8F(B+1)` idle density payload
is absent, independently of whether zero pages would have been resident.

### Commands and unresolved gates

All scratch stayed in `target/scratch/refstart` inside the worktree.

- `DFT_SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF" cargo test --offline --release -p dft-engine --test refstart -- --include-ignored --nocapture --test-threads=1`: **exit 0**, 6 passed, 0 failed. Record contains no timing claim.
- `CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh dft-engine`: **exit 0**, build/test/fmt/clippy PASS; 271 passed, 0 failed, 79 ignored.
- `CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh full`: **exit 1**, build/fmt/clippy PASS; **808 passed, 3 failed, 185 ignored**. Existing CLI tests `pseudo::cached_fetch_checks_hash_and_provenance`, `pseudo::cache_resolution_and_no_network_side_effects`, and `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites` reject their in-repository temporary cache with "pseudopotential cache must be outside repositories". Same pre-existing conflict recorded in G-OWNED-ROWS. No test/path guard weakened; files outside ticket ownership left unchanged. The optional load report also cannot run `ps` in this sandbox; gates still execute.
- `python3 -I scripts/make_readme_figures.py`: **exit 0**, generated files unchanged. New evidence is in the ticket-owned tests directory; extending the shared README evidence generator is outside this ticket's ownership.
- `python3 -I scripts/tests/test_readme_figures.py`: **exit 1**, 8 passed / 1 failed; existing `test_structural_scope_and_sources` still expects the outdated phrase "standard profile is unsupported". No README/test wording changed to mask it.
- Initial compilation exposed test sharing/type errors, corrected before the numerical runs; first full Clippy found the new helper placed after a test module, corrected by moving it before the module. Final Clippy passes; all failed attempts remain in worktree scratch logs.

[Dev08 request](lanes/RSCF/start-dev08-request.md): Si64/128/216 MPI8, one case
per job, A/A floor then crossed A/B, same final-cutoff atomic start, first
Davidson boundary including workspace setup and all-rank RSS. A host adapter
must consume `StartState` directly before legacy problem construction; this
ticket does not add a CLI flag or edit the R-SCF-1 driver. LDA without NLCC is
implemented; PBE, NLCC and frozen strain reject explicitly. No new default,
speed claim, full allocation-peak claim, or QE result is adopted.
## R-SCF-1 — independent reciprocal-density reference driver (2026-11-23)

Derived from note `docs/research/qe-study/qe-density-mixing-spec.md`; no QE
source was accessed. This is issue #25 patch 1, an opt-in verification surface,
not an adopted optimization and not a revival of a rejected mixing heuristic.
The motivating dev08 map-count/overhead decomposition is recorded in
`docs/design/reference-scf.md`; no new Mac timing claim is made here.

| Attempt | Measured correctness/work | Decision and lesson |
| --- | --- | --- |
| Reciprocal half-sphere Broyden, unnormalized secants, beta=0.7, m=8, spectral symmetric solve | Closed-form 2D linear fixed point and ring eviction pass. Three-step synthetic file replay has maximum next-density error 2.168404344971009e-19 e/bohr^3; largest relative shell error 2.168404344971009e-16; charge difference 0 electrons. Serial and simulated MPI3 replay agree. Fourier/metric tests cover simulated ranks 1/2/3/8; owned driver exercises the cold-map verification boundary and two maps at MPI2. | Default off. Singular/numerically rank-deficient A is reported, never regularized or silently reset. Real QE data and measured tolerance are still pending dev08. |
| Final-cutoff atomic-start driver, same rule in Ry and Ha | Seed 801, fractional displacement amplitude 1e-3, 1 thread: Si8 stops in 10 iterations, E=-31.32178356464003 Ha, dr2=9.025138772492297e-9 Ry; Si16 9 iterations, E=-63.12910143369001 Ha, dr2=6.809305771391534e-9 Ry. `verify=off` exits 0. Default requested verification exits 3 with an explicit patch-2 TODO and no convergence/certification claim. | Runnable skeleton; not delivered-accuracy or speed evidence. Fixed-point verification is too coupled to session control to expose unchanged without importing the legacy switches. No default promotion. |
| Legacy default isolation, base c534359c7ca6a7b2ee18694f905c453839bcb764 | Si8/Si16/Si32/skew: dE=0 Ha, density bytes, eigenvalues and all scientific iteration history fields bit-identical (wall fields excluded); 11/10/12/8 iterations respectively. Separate `bench_scf.py --cases si8 --require-identical` also passes, dE=0 Ha. | No default numerical change. Four fingerprints and reproduction inputs in `lanes/RSCF/local-validation.json`. |
| Synthetic low-cutoff negative-ringing probe | A narrow Gaussian atomic charge (rc=0.7 Bohr, cutoff 2 Ha, deliberately small test cell) generated -0.0014546618092563957 e/bohr^3 after mixing and correctly failed the existing density validator. A resolved broad Gaussian (rc=1.8 Bohr) tests the two-map serial/MPI bridge. | No clipping or tolerance change was added to make the unresolved synthetic model run. Such inputs remain explicit errors; reference mixing is not a positivity-preserving method. |

Native MPI/MKL and QE 7.6 density comparison are requested in
`lanes/RSCF/dev08-request.md`. The new mixer keeps only owned numerical G rows
and owned FFT planes during MPI iterations; inherited problem construction and
final delivery still use the existing engine's full-array boundary. MPI PBE,
NLCC/floor and symmetry bridges are not implemented and fail explicitly.
Kernel internals retain their existing backend settings; the new orchestration
uses no `ScfSession`, cutoff ladder, or legacy mixer/stop scheduling.

Whole-workspace validation with worktree-local TMPDIR: build/fmt/clippy pass;
813 tests pass, 3 existing pseudopotential-cache tests fail because they require
a cache outside every repository while the ticket requires scratch inside this
worktree. Those files are outside this ticket's ownership. README regeneration
has no diff; its existing specification test also fails on the stale expected
phrase `standard profile is unsupported` (current generated text supports the
standard path after ADR 0043). These failures are retained in the local
validation record, not reported as passing. Dashboard/session/README result-type
updates require orchestrator-owned files; no numerical gate is claimed adopted.

### R-SCF-1 delivered reference follow-up (2026-11-23)

Production Broyden replay of genuine QE 7.6 Si8 k1..3 passes: max density error 3.46945e-18 e/bohr^3; all ten lowest shells, charge and rounded dr2 pass. Initial beta 0.7/0.3 recovery and beta 1 output pass. See [measured evidence](lanes/RSCF/qe-si8-validation.json). No timing claim. Full gate retains three accepted local TMPDIR cache failures; patch-2 verification TODO remains.
## R-SCF-DAV — standalone reference-SCF inner solver (2026-11-23)

**Implemented, not wired into the driver, no adopted default or speed claim.**
New premise: the user-authorized reference-SCF block organization and explicit
inner stopping semantics. Derived from note
`docs/research/qe-study/davidson-iteration-structure-gamma.md` only; no QE source
read. REL-SLICES, Advance/bounded-subspace, restart-keep, PR-CG, UPPER_BORDER,
ROTATE_COLUMNS and packed reductions stay closed. The supplied motivation's
183 reduced calls and 58 iterations have different stage scopes; no per-fine-map
redundant-solve inference is made from their quotient.

Fresh all-band H at entry, active-band corrections, bordered H/S, generalized
reduced solve via the existing guarded Cholesky/faer adapter, capacity 2n,
GEMM-only contraction of V/HV, and explicit occupied/empty Ritz-change thresholds.
No retained H input. Root and replicated EVD are exposed for a same-problem
measurement. Per-call/per-iteration H widths, active indices, EVD orders/calls,
repairs/restarts, collective payloads and GEMM shapes/flops are returned; stage
clocks are optional. Ordinary driver/legacy source remains untouched.

**Correctness finding, not a speed optimization:** the ideal identity/diagonal
reset after contraction/whitening accumulated 1.0154099783221682e-11 Ha error
on a synthetic degenerate spectrum. Re-measuring H/S from rotated V/HV fixes
it with no new H or extra EVD; the extra borders/factor work is charged. The
original test tolerance was retained. Synthetic roots now differ from full EVD
by at most 6.54e-14 Ha in the degenerate fixture, including ranks 1/2/3/8 and
root/replicated placement; occupied projectors agree with dense and legacy.
This refinement may affect speed and must be priced, not omitted from the model.

Mac fixed-H evidence (not timing acceptance): 48/48 converged replay solves,
Si8/16/32/skew x early/middle/late x simulated ranks 1/2/3/8; solver ethr 1e-12 Ha,
200-iteration validation budget. Maximum occupied difference from strict legacy
residual 1e-12 Ha reference **1.1123879595231756e-12 Ha**; maximum occupied
principal-angle sine **7.747457492838181e-6**. The stopping ethr is not an absolute
error bound (Si16 middle slightly exceeds ethr). Iteration/H pairs, identical
across these rank counts: Si8 30/326,8/93,8/88; Si16 21/530,9/241,4/121;
Si32 28/1045,9/487,3/209; skew 9/60,5/31,2/17. Full 679-coordinate skew middle
oracle: <=6.189493362285248e-14 Ha, ranks 1/2/3/8. Complete degenerate occupied
edge synthetic test passes; arbitrary spectrally incomplete invariant starts
still require an independent driver-level audit.

Negative stress retained: Si8 early, forcing ALL bands including empties to
1e-13 Ha, exhausted 200 iterations, 17/20 flags, 1065 H columns (occupied error
5.003e-15 Ha). No empty threshold or tolerance was relaxed to hide it. The
normal first-step setting is loose; the default 20-call budget explicitly
returns unconverged flags and the bench exits 5 on exhaustion.

Si8 late simulated rank 3: profile off/on and root/replicated have identical
bits, 5 iterations/62 H columns; default replay matches the HEAD example on
four systems. Full check build/fmt/clippy pass; three pre-existing pseudo-cache
tests fail with repository-local TMPDIR, as in preceding tasks. No workaround
moves scratch outside the owned worktree. The allowed note does not contain
the preconditioner formula; the implemented smooth formula is explicit in the
design and needs architect confirmation before exact QE iteration claims.

Evidence and exact dev08 capture/A-A/crossed A-B commands (Si64/Si128 MPI8):
[request](lanes/RSCF/dav-dev08-request.md),
[inner contract](../design/reference-scf.md#inner-solver),
[stored Mac evidence](../../crates/dft-solvers/tests/ref_davidson_evidence.json).
The dev08 stage-0 result, root-placement choice, and full-SCF integration are
pending through the master; no Mac timing is used to choose a default.


### R-SCF-DAV preconditioner addendum follow-up (2026-11-23)

The main study-note addendum resolves the formula question above. The smooth
Ha denominator was already algebraically correct; the physical replay supplied
only kinetic energy. It now passes kinetic + mean total local potential + the
real nonlocal projector-pair diagonal, with NC overlap 1, shared across Gamma
pair coordinates without extra weights. The API now names this full diagonal
explicitly. Existing solver comparison arms retain their kinetic inputs.
Derived from note `davidson-iteration-structure-gamma.md`, addendum dated
2026-11-23; no QE source consulted and no production default changed.

All 48 physical snapshot solves converge on simulated ranks 1/2/3/8, with
identical iteration and H-column counts across ranks. Early/mid/late iterations:
Si8 28/8/8; Si16 21/10/4; Si32 30/9/3; skew 10/5/2.
Maximum occupied error against strict legacy: 8.92924623e-13 Ha;
maximum occupied principal-angle sine: 7.24731268e-6;
maximum occupied rank difference: 4.41313652e-15 Ha.
Synthetic tests: 5 pass, including degenerate occupied edge (48 iterations,
237 H columns at every rank count) and the independent Ry formula/shift test.
Dense skew and root/replicated/profile comparisons pass. These are correctness
observations, not speed acceptance or measured QE iteration identity.

`CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh full`
exits 1: build/fmt/clippy pass, 807 tests pass, 181 ignored, the same three
pseudo-cache tests fail because their temporary cache is inside the repository.
README generation exits 0 and produces no tracked changes; its Python test
exits 1 on the existing missing phrase `standard profile is unsupported`.
No unrelated test, tolerance or scope change was made. Evidence JSON and the
request table above are refreshed; historical kinetic-only negative experiments
remain explicitly labelled as historical. No dev08 measurements performed.

## R-SCF-MAP — same-input density diagnostic (2026-11-23)

Diagnostic infrastructure, no optimization/default adoption. Derived from the
own-words density-mixing / SCF-cycle notes and the existing reference-Davidson
contract; no QE source consulted. New `refscf-map` teacher-forces QE's input
at every step, calls existing Hamiltonian/density kernels with legacy Davidson
(default) or the opt-in reference solver, and retains all mismatch rows before
exit 4. The reference bridge is isolated from production dispatch. G matching,
charge, ten shells, mixer metric, printed energy/bands/occupations, analytic
Hartree and available potential-response checks are recorded. No direct QE
XC/local potential grid is present; missing comparisons are null.

Mac Si8 one-rank k1..k3, 20 bands, 30/120 Ry, 36^3, printed ethr
1e-2 / 1.37e-3 / 6.11e-5 Ry. Legacy max density errors:
1.412541591e-4 / 4.016489166e-5 / 1.180480939e-5 e/bohr^3;
reference: 1.411562359e-5 / 1.026520345e-5 / 4.814860351e-6 e/bohr^3.
Mixer distances respectively 1.387939016e-3 / 1.068642280e-4 /
2.627257499e-5 Ry and 1.966403362e-5 / 1.282576976e-5 /
5.596930412e-6 Ry. Both first fail at **k=1 rho_out all G** against the
2.72e-16 e/bohr^3 equality tolerance (4 times the user/dev08-relayed QE rank
floor upper end, not independently remeasured on Mac). The disagreement is a
finding, not a relaxed correctness gate or an attribution to mixing.

Input FFT roundtrip <=4.423e-18 e/bohr^3, charge error <=3.553e-14 electrons,
analytic input Hartree difference <=8.442e-17 Ha. Hartree/XC energies computed
on QE's output density agree within QE's 5e-9 Ry print rounding; Ewald differs
by 1.203e-8 Ry. At common ethr=1e-12 Ry the two solvers' occupied eigenvalues
agree to 1.016e-11 eV. Simulated reference ranks 2/3 retain iteration/H counts
and differ from serial density by <=2.914e-16 / 2.498e-16 e/bohr^3 in real
space. QE start orbitals and unrounded thresholds are absent, so loose-map
mismatch cannot by itself distinguish the potential from inner-solve/start
error. Neither the smaller reference error nor these Mac numbers select a
production solver default.

Release map tests exit 0 (4 passed); parser tests exit 0 (2 passed). Full
`CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh full` exits 1:
build/fmt/clippy pass; 827 tests pass, 188 ignored, the same 3 repository-local
pseudo-cache failures remain. README test exits 1 on the pre-existing missing
phrase `standard profile is unsupported` (8 pass, 1 fail). Default A/B against
pre-edit HEAD: Gamma/k222s/Si8 energy delta exactly 0 Ha, density hashes identical,
SCF iterations unchanged 4/8/4, verification enabled. No Mac speed claim.

[Complete stored evidence](../../crates/dft-cli/tests/fixtures/refscf-map-si8.json),
[reproduction/table/dev08 request](lanes/RSCF/map-dev08-request.md). The request
covers Si128 np1/np8 and Si8 np1, both solver arms, per-step stdout ethr,
provenance and attribution follow-ups. Native MPI/MKL and Si128 measurements
remain with the master/dev08 session. Changes left uncommitted.
### R-SCF-2 — independent full reference SCF (2026-11-23)

Integration only; **not adopted as default, acceptance incomplete**. Derived
from the own-words notes `qe-density-mixing-spec.md`,
`davidson-iteration-structure-gamma.md` and `qe-scf-cycle.md`; no QE source read.
The opt-in driver now consumes parsed input directly and uses refstart +
owned-row ref_davidson + patch-1 Broyden, with orbitals retained between maps,
Ry/Ha threshold mapping, first-step redo and a fresh fixed-point verification.
A failed verification continues the original history (ADR 0043 V5); Si8 first
stops at 10, fails verification at 3.06723442e-8 Ha, and passes at step 12.
Records expose every call/active width/order/H column, energies, history
conditioning and verification cost. Memory payload bounds are distinguished
from Linux VmHWM (Mac null) and external whole-process measurements.

Mac correctness matrix: Si8/16/32, established skew Si2, elongated Si24 and
closed-shell synthetic nonlocal PP, each at simulated ranks 1/2/3/8. All 24
converge and pass the fixed-point/5e-9 Ha class checks. Maximum strict-oracle
energy differences by case: 1.43884904e-10, 1.09466214e-10, 3.80993015e-11,
7.59472485e-11, 6.06945605e-11 and 1.69247072e-9 Ha. Exploratory open-shell
single-atom z=4 and sheared Si8 strict oracles did not converge; no result is
claimed for them. The accepted fixture set and all caveats are in
[patch2-validation.md](lanes/RSCF/patch2-validation.md).

**Unresolved:** Si16 rank 3 differs at SCF step 12. A 9.15933995e-16 Ha Ritz
change difference straddles ethr: 2 versus 3 corrections, 92 versus 93 H
columns; dr2 differs 3.15849833e-11 Ry and final density 8.17472951e-8 e/Bohr^3.
Other five cases pass trajectory tests. No threshold/tolerance widening or
forced extra iteration was introduced. Lesson: a passed delivered accuracy
class does not establish rank-identical threshold branches; preserve the
failing test and per-band Ritz evidence before a determinism decision.

Genuine committed QE Si8 stdout: first divergence at step 1 (dr2
0.43855533 versus 0.4347506747835862 Ry); ethr/correction counts agree initially
(0.01 Ry, 3). Already before solving, recovered initial rho(G) differs by
8.64045255e-7 e/Bohr^3 max. Direct radial integration/start clipping versus
QE's documented table/no-clipping start has not been separated further;
no claim that Broyden caused the full-trajectory mismatch. The user's relayed
patch-1 pass is retained with measured-floor density-replay limits
1e-17 + 1e-12*amplitude (same rank), 1e-16 + 1e-12*amplitude (cross-QE-rank).

Pristine same-HEAD comparisons: four default CLI cases have bit-identical
density/energy/eigenvalues and timing-stripped histories/verification records;
`bench_scf.py` gamma same-HEAD identity passes. The historical baseline rejects
standard-versus-strict (exit 3); strict mode differs in density with energy
3.55e-15 Ha (exit 5 with require-identical), so it cannot establish identity of
this patch. No timing ratio is used. `scripts/check.sh full` exits 1:
build/fmt/clippy pass, 825 tests pass, 188 ignored, the same three worktree-cache
tests fail. New release CLI tests 2/2, Broyden tests 10/10, refstart 6/6 and
ref_davidson 5/5 pass. Full matrix exits 101 on the retained Si16 rank failure.
README generation exits 0/no tracked change; its old scope-text test fails.

[Full dev08 request](lanes/RSCF/run-dev08-request.md): Si128/64/216 np8 and
Si32/64 np1, matching full QE references, one case per `final-` job, 7200 s,
full trajectory/verification/rank-memory evidence. No remote run performed;
no default or speed decision is justified yet.

### I26-2-DIAGLEVEL dev08 A/B (2026-11-24): diagnostic level, no measurable wall effect

Hypothesis: forced full diagnostics (per-band Ritz traces, inner clocks, duplicated JSON) inflate
delivered wall of the reference driver. Change: lightweight level is the default since main 6cdda2e;
full via BAANDI_DAVIDSON_STAGES=on. dev08, release mpi,mkl binary sha 6c599b30..., MPI np8 bound to
cores, one thread per rank, verification on in both arms, wall around the whole mpirun incl. output.
A/A 20 pairs then crossed A/B 30 pairs per case; 206/206 deliveries verified, science identical
(density.bin byte-identical on every rank).

| case | A/A floor F | median full / light (s) | mean log(light/full) +- SE | abs/SE | verdict |
|---|---|---|---|---|---|
| si64 np8 | 0.00536 | 9.491 / 9.540 | +0.00293 +- 0.00479 | 0.61 | inconclusive |
| si128 np8 | 0.00560 | 56.469 / 56.367 | -0.00141 +- 0.00447 | 0.31 | inconclusive |

Bytes: light removes ~6.9 MB (si64) and ~13.9 MB (si128) of output per run (results.json 1.16 -> 0.30
MB and 2.13 -> 0.40 MB). Lesson: at MPI8 on dev08 the diagnostic payload costs below the 0.5 % floor;
the QE gap is not explained by tracing. The lightweight default stays (smaller output, no cost).
Evidence: /data/sig-dft/i26/diaglevel/{si64,si128}/ on dev08 (TABLE.md, table.json, samples-*.tsv).

### I25-VERIFY-COST dev08 A/B (2026-11-24): V2 predicted-fail skip ADOPTED, V1 reuse kept off

Reference driver, standard 5e-9 Ha, fixed-point verification criterion unchanged. dev08 MPI8, release
mpi,mkl 0ce77fdc (sha 12011068...), one thread per rank, lightweight diagnostics. 20 A/A + 30 crossed
pairs per variant and case, 672 samples, all VERIFIED, |E - E_strict| <= 2.4e-10 Ha, science identical
to each arm's setup run. Beyond noise = |mean|/SE >= 3 and |mean| >= 3F.

| case | A/A floor F | V1 reuse log(B/A) | V2 gate log(B/A) | median wall A -> V1 / V2 (s) |
|---|---|---|---|---|
| rskew (2 atoms) | 0.0229 | +0.0062 +- 0.0170 inconclusive | +0.0057 +- 0.0157 inconclusive | 0.176 |
| elongated (24) | 0.0095 | +0.0600 +- 0.0090 slower | -0.0106 +- 0.0081 inconclusive | 1.511 -> 1.621 / 1.504 |
| si64 | 0.0058 | -0.0927 +- 0.0044 faster | -0.0557 +- 0.0057 faster | 9.448 -> 8.637 / 8.929 |
| si128 | 0.0066 | +0.0443 +- 0.0039 slower | -0.0670 +- 0.0043 faster | 56.63 -> 59.47 / 53.06 |

V2 skips one predicted-failing attempt (prediction = previous attempt error scaled by the dr2 ratio; no
constant tuned to a system) and never delivers without a pass: faster on si64 and si128, neutral elsewhere
-> default on (REFSCF_VERIFY_GATE). V1 helps only when reuse removes an attempt (si64); on si128 and the
elongated cell it adds 2 to 3 logical iterations and recovery time -> stays off. Limitation: all cells are
Si with one pseudopotential; the slab cell needs smearing and was excluded. Lesson: the first attempt
after the estimator stop fails by 10-30x; a predictor for the FIRST attempt is the next lever.
Evidence: dev08 /data/sig-dft/vc/<case>/{TABLE.md,table.json}, setup.md.

### I25-FIRST-ATTEMPT dev08 screen (2026-11-24): first-attempt predictor REJECTED (stays off)

REFSCF_VERIFY_GATE_FIRST (secant estimate of the fixed-point residual from the Broyden history) at
de2f1d29, untimed screen, MPI8 mpi,mkl, four cells, all VERIFIED, |E - E_strict| <= 1.2e-10 Ha.
Evaluated maps default -> first: si64 16 -> 16, si128 20 -> 20, elongated 18 -> 16, rskew 9 -> 10.
The doomed first attempt is skipped correctly (si64 predicted 5.17e-8 vs counterfactual 5.80e-8), but the
next candidate still misses (6.75e-9 / 3.09e-8 vs 5e-9), so the failure only moves one map later. On
first-pass cells the secant overpredicts 50-700x (rskew 6.48e-8 vs actual 9.69e-11) and costs a map.
Pre-registered stop: no timing run. Lesson: with the criterion fixed, skipping attempts is exhausted;
the remaining tail cost is the 2 extra maps needed to bring the delivered fixed-point error from 10-30x
the threshold to below it. Also found: maps after a skipped attempt are counted as "unassigned" by the
stop-to-verified reader (si8 0.041 s); to be fixed with the next reader change.
Evidence: dev08 /data/sig-dft/fa/SCREEN.md, screen.json.

### P3 timed (2026-11-24): REF_DAV_WIDE_ROTATE, REF_DAV_RESTART_CORRECTIONS, REFSCF_SETUP_PHASE_REUSE ADOPTED

Cause found from the dev08 gap profile: the reference Davidson's rotations (V C, HV C; si128 184 calls,
186 GFLOP per rank) went through the plain GEMM entry, whose whole-output bound sent every real-cell
rotation to the pure-Rust kernel (~15 GFLOP/s per core). W routes them through the wide entry (MKL
panels). R omits the identity rotations right after a restart (4rnu FLOPs per affected correction).
S reuses structure factors in the atomic start. dev08 MPI8, mpi,mkl 0e8b08da (sha 5e8b1a3e...),
lightweight diagnostics, verify on; 20 A/A + 30 crossed pairs per comparison; 842 samples all
VERIFIED, |E - E_strict| < 5e-9, maps/attempts identical across arms.

| case | W vs default log ratio (x) | X=W+R+S vs W log ratio (x) | median wall default -> W -> X (s) |
|---|---|---|---|
| rskew | -0.0005 +- 0.0221 inconclusive | -0.0162 +- 0.0201 inconclusive | 0.177 |
| elongated | -0.0009 +- 0.0064 inconclusive | -0.0162 +- 0.0049 (1.8F) inconclusive | 1.484 -> 1.503 |
| si64 | -0.1048 +- 0.0063 faster (x0.900) | -0.0089 +- 0.0047 inconclusive | 8.915 -> 8.000 |
| si128 | -0.2515 +- 0.0040 faster (x0.778) | -0.0290 +- 0.0042 faster (x0.971) | 52.719 -> 41.126 -> 40.056 |

Adoption rule (2026-11-10): faster where measurable, no regression beyond noise anywhere -> all three
default on (W is rounding class: H columns 14204 vs 14203 on si128, energies unchanged). Lesson: check
gemm_dispatch routes first whenever a dense stage is slow; the pure fallbacks are silent.
Evidence: dev08 /data/sig-dft/p3/<case>/TABLE.md, table.json, SCREEN.md.

### REFSTART_MKL (2026-11-24): ADOPTED

Atomic-start Gram and rotation products moved from the direct pure-Rust kernel (si128: 7 calls, 17 GFLOP
per rank, route pure_direct) to the dispatching GEMM entries. dev08 MPI8 mpi,mkl 66d45e6, lightweight
diagnostics, 20 A/A + 30 crossed pairs per case, all samples VERIFIED and identical to their arm's setup.
si64 log ratio -0.0178 +- 0.0036 (x0.982, 7.951 -> 7.826 s); si128 -0.0347 +- 0.0034 (x0.966,
40.501 -> 39.075 s); elongated and rskew inconclusive. Not slower anywhere -> default on. Setup on si128
4.83 -> 3.36 s (profiled screen). Evidence: dev08 /data/sig-dft/p3b/<case>/TABLE.md.

### P3c timed (2026-11-25): REF_DAV_SINGLE_GEMM and FFT_DIRECT_PACK/FFT_SPARSE_CLEAR not adopted (inconclusive)

dev08 MPI8 mpi,mkl e62f85a6, lightweight, 20 A/A + 30 crossed pairs per arm and case, all samples verified.
GF (all three) vs default: si128 -0.0181 +- 0.0049 (z 3.7 but 2.7F < 3F; medians 40.488 -> 39.726 s),
si64 -0.0103 +- 0.0056, elongated/rskew inconclusive. F alone (bit-identical): inconclusive everywhere
(si128 -0.0047 +- 0.0037). Profiled screen showed border 5.07 -> 3.59 s with G, so a G-only campaign
with more pairs may resolve the si128 trend. Switches stay default off. Evidence: dev08 /data/sig-dft/p3c/.

### REF_DAV_SINGLE_GEMM provisionally ADOPTED (2026-11-25, development protocol)

User rule 2026-11-25: during development single-run screens decide; precise campaigns later for the stacked
defaults. Evidence: P3c screen border 5.07 -> 3.59 s (si128), P3c GF vs D si128 -1.8% (z 3.7, below 3F),
partial P3d G-only si128 (9e99bd3, 7 A/A + 20 crossed pairs): medians 39.126 -> 38.518 s, mean log ratio
-0.0120 +- 0.0045. All samples verified, counts identical, E - E_strict unchanged. Default on (rounding
class); to be re-measured in the consolidated precise campaign. Evidence: dev08 /data/sig-dft/p3d/si128/TABLE.md.

### P3f/P3g screens (2026-11-25, development protocol, one run each, dev08 MPI8)

- REFSCF_VERIFY_ENERGY (ADR 0044): energy change ~0.85x the Coulomb error at every attempt; si64 7.62 ->
  7.28 s, si128 38.39 -> 40.64 s (gate inactive -> 3 attempts), elongated 1.51 -> 1.69 s. Default off
  (user decision C); QE comparisons now carry verify-off and verified BAANDI walls.
- FFT_FUSED_ROWS (P3-FFT2, bit-identical): si128 42.30 -> 41.53 s, si64 7.52 -> 7.92, elongated
  1.56 -> 1.67, rskew 0.170 -> 0.250; slab_gather and multiply slower. Not adopted, stays off.
- QE-matched clock (P3-QECLOCK, si128 one run each): h_psi 0.154 vs 0.137 s/call (93 vs 82 calls),
  Davidson remainder 9.25 vs 5.29 s, reduced eigensolve 2.25 vs QE serial rdiaghg 4.17 s, setup 3.38 vs
  init_run 0.83 s, density 2.12 vs sum_band 1.84 s. QE ran -ndiag 1 (serial), not its earlier best 2x2.
Evidence: dev08 /data/sig-dft/{p3f,p3g}/SCREEN.md, /data/sig-dft/qeclock/out/.

### P3i (2026-11-25, development protocol): four bit-identical switches ADOPTED; first efficiency table

REF_DAV_SWAP_BUFFERS, REF_DAV_LEAN_RESET (P3-DAVREM), REFSTART_DIAGONAL_REUSE, REFSTART_ORBITAL_REUSE
(P3-SETUP2): dev08 MPI8 one profiled run each at 0f11a240, all bit-identical on four cases. Wall
si64 7.75 -> 7.41, si128 38.15 -> 37.05, elongated 1.61 -> 1.43, rskew equal; si128 setup 3.39 -> 2.75 s
(nonlocal_diagonal 0.549 -> 0.080, atomic_orbital_fill 0.377 -> 0.114). Default on.
si128 per-stage efficiency (rank max, 20 solves): apply_h 15.78 s (alltoall 4.1 GB); restart 1.94 s at
44 GFLOP/s per rank; border S/H 1.85/1.81 s at ~50; corrections 1.11 s at 33; return 1.03 s at 47;
precondition 1.21 s for 0.7 GFLOP (scalar hypot bound); Davidson setup 0.96 s streaming 3.5 GB; normalize
0.31 s at 0.2 GFLOP/s. GEMMs are now at 60-80% of the ~60 GFLOP/s per-core peak; the non-GEMM passes
(precondition, setup, normalize, append, apply_check: ~3 s) are the cheap targets. Evidence: dev08
/data/sig-dft/p3i/REPORT.md.

### P3j (2026-11-25, development protocol): P3-STREAM switches ADOPTED

REF_DAV_SQRT (vectorized residual/precondition/norm, sqrt instead of hypot; rounding class), REF_DAV_FUSED_APPEND,
REF_DAV_CHECK_APPEND, REF_DAV_REUSE_ROWS. dev08 MPI8 one profiled run each at ec553c93, all verified, counts
identical (si128 H columns 14204 -> 14203), |E - E_D| <= 1.7e-13 Ha. Wall si64 7.48 -> 7.18, si128 37.23 ->
36.06, elongated 1.55 -> 1.45, rskew equal. si128 precondition 1.251 -> 0.252 s (0.6 -> 3.3 GFLOP/s),
setup 0.934 -> 0.699 s, append copy moved into normalize/apply_check (net -0.06 s). Davidson stage sum
29.21 -> 27.96 s. Default on. Evidence: dev08 /data/sig-dft/p3j/REPORT.md.

### P3l (2026-11-25, development protocol): REFSTART_GENERALIZED, REFSTART_COMPACT_G, NONLOCAL_SINGLE_GEMM, REFSCF_LAZY_OUTPUT ADOPTED

dev08 MPI8 one profiled run each at 31b23e11, all verified, counts identical, |E - E_D| <= 2.8e-13 Ha.
Wall si64 7.20 -> 7.13, si128 35.67 -> 34.60, elongated 1.49 -> 1.42, rskew equal. si128 setup 2.647 ->
2.297 s (generalized Ritz start replaces two QR passes), output 0.666 -> 0.009 s (orbitals stay owned),
nonlocal projection 1.986 -> 1.708 s (124 projection GEMMs from 256-row panels to one call); add_vuspsi per
call unchanged. Evidence: dev08 /data/sig-dft/p3l/REPORT.md.

### P3m (2026-11-25, development protocol): REF_DAV_IDEAL_RESET, NONLOCAL_FUSE_COUPLING, NONLOCAL_BACK_ACCUM ADOPTED

The ideal reset (S = I, H = diag(theta) after a restart, orthogonality guard, full refresh every 4th
restart) was rejected earlier under a bit-identity premise (1.0e-11 Ha drift); re-evaluated under the
strict 1e-10 / standard 5e-9 classes. dev08 MPI8 one profiled run per arm at c9abca1b, all verified,
counts identical, |E - E_D| <= 2.3e-13 Ha, E - E_strict <= 1.2e-10. Wall D/R/N/RN si64 6.89/6.83/6.97/6.79,
si128 35.57/33.47/33.70/34.33 (RN single-run noise; stage savings additive), elongated 1.44/1.31/1.28/1.26.
si128 R: border S+H 3.75 -> 2.80 s (26 fewer GEMM pairs), factor/transform -0.16 s; N: add_vuspsi per call
0.0250 -> 0.0203 s (back_and_add 2.37 -> 1.97 s). Default on. Evidence: dev08 /data/sig-dft/p3m/REPORT.md.

### P3q (2026-11-26, development protocol): REFSCF_LEAN_ENERGY ADOPTED

Per-map energy only where consumed, fields of the mixed density reused as the next map's input fields,
full audit at delivery and verification. dev08 MPI8 one profiled run per arm at b7e6cfe8: delivered
density.bin and energies bit-identical on all cases (intermediate audit fields null), counts identical.
si128 per-map energy 1.168 -> 0.275 s (verified) and 1.025 -> 0.173 s (verify off); si64 0.235 -> 0.059;
elongated wall 1.48 -> 1.21 s. energy_correction (~0.25 s) untouched. Default on.
Evidence: dev08 /data/sig-dft/p3q/REPORT.md.

### P3r (2026-11-26, development protocol): REFSCF_RMM (QE-style rmm-davidson) not adopted, stays off

Implemented from the own-words note docs/research/qe-study/qe-rmm-davidson.md (Davidson at SCF iterations
1, 2 and every 5th, one RMM-DIIS sweep of at most 4 steps otherwise, verification maps Davidson). dev08
MPI8 one run per arm at e5ec8ae: slower everywhere because it adds maps and H columns. si128 verified
35.49 -> 43.56 s (17 -> 20 maps, H columns 14203 -> 17405), verify off 28.56 -> 31.57 s (15 -> 16 maps);
si64 6.70 -> 7.62, elongated 1.29 -> 1.34. RMM sweeps on si128 cost 21.1 s for 14 sweeps (apply_h 11.7,
Gram-Schmidt 3.9, DIIS 2.3) against QE rrmmdiagg 6.6 s for 10 calls; QE rmm itself is only 2-3% faster than
QE Davidson on these cells (si128 21.89 vs 22.40 s). Lesson: on Si up to 128 atoms the Davidson trajectory
already matches QE's; RMM pays only if its sweep overheads (GS, DIIS) are brought to QE's level and the map
count does not grow. Kept as a switchable module (generality rule). Evidence: dev08 /scr/mino/rustdft/p3r/.

### P3t (2026-11-26): REFSTART_FFT_SHARE and REFSTART_SINGLE_GEMM ADOPTED, REFSTART_LOCAL_SHELLS kept off

dev08 MPI8 one profiled run per arm at 79cf8703 (all three on vs off): all verified, counts identical,
|E - E_D| <= 2.3e-13 Ha; not bit-identical on dev08 MKL (Mac pure build was), so rounding class. Wall
si64 6.56 -> 6.66 (noise), si128 32.47 -> 32.15, elongated 1.35 -> 1.17. si128 setup 2.343 -> 2.254 s:
density FFT plan 0.252 -> 0.115 (FFT_SHARE), atomic overlap Gram 0.086 -> 0.062 (SINGLE_GEMM); LOCAL_SHELLS
showed no gain (local tables 0.415 -> 0.457) and stays off. Evidence: dev08 /data/sig-dft/p3t/REPORT.md.

### Work-order screen (2026-11-26, user work order "reduce full-cost SCF work"), dev08 Si64 MPI8, one run per arm at bb416c73

Archived vc Si64 input (sha af97b69b...), strict -253.3781677212637 Ha, socket 0 cores 0-7, verify on, 120 s limit.
| arm | status | wall s | E - E_strict | fine maps | fine full-H cols | attempts | peak RSS/rank |
|---|---|---|---|---|---|---|---|
| D | verified | 6.41 | +1.15e-10 | 13 | 6094 | 2 | 206 MB |
| H (REFSCF_KEPT_HX) | verified | 6.60 | +1.15e-10 | 13 | 4642 (-24%) | 2 | 259 MB |
| L (REFSCF_LADDER) | verified | 6.73 | +4.86e-11 | 10 (+9 coarse) | 5166 | 2 | 219 MB |
Decisions (work-order wording): H "work reduced, speed unresolved": 11 cache hits, 1452 starting full-H columns
replaced, solve -0.15 s but density +0.48 s (retained-block density pass) and +53 MB/rank. L "not a follow-up
candidate as configured": fine saving ~1.22 s < coarse stage 1.72 s (9 coarse maps to 1e-6 Ry). Both stay off.
Also si128 FFT_KERNEL_BY_LENGTH (now wired): all 10 lane shapes choose MKL (16-40% faster per batch in
calibration); wall 33.63 -> 33.46 s; kernels -0.9 s but slab gather/scatter, middle pack/unpack and multiply +0.5 s
(same pattern as M2_FFT=mkl in P3e). Evidence: dev08 /data/sig-dft/wo/RESULT.md.

### PULAY-OPTIN (2026-11-26, user decision): stress Pulay assessment OFF by default, opt-in

Trigger: dev08 stress runs where the isotropic cutoff-derivative Pulay assessment (three extra cold SCF runs at
c, 1.05c, 1.10c) was 51% of the wall on GeTe and 71% on B12 (1861 s of 2626 s); QE does no such scan and the
assessment never changes the stress tensor (diagnostic only). Change: off by default, opt-in with
`[output] stress_pulay_assessment = true` (CLI, also for the nominal-cutoff frames of cell relaxation) or
`Session.stress(pulay_assessment=True)` (Python); engine switch `stress::PulayAssessmentMode` on
`refscf::api::Derivatives` and `stress::stress_with`. Off: `stress.pulay_assessment =
{"status":"not_evaluated",...}`, `cutoff_scf_runs = 0`. Stress tensor bit-identical on/off (CLI test
`all::stress::stress_output_identity_schema_and_paths`, spin/FD in `forces.rs`, facade deferred vs eager in
`api_identity`). Status: adopted (default off). Lesson: a diagnostic that costs more than the observable belongs
behind an explicit switch.

### MPI-DERIV (2026-11-26, Mac correctness; dev08 acceptance pending): k-distributed post-SCF derivatives

Trigger: dev08 MPQE entry 17 (B12, 172 k, MPI8) spent >36 min after a ~9 min SCF+verification in forces/stress.
Change: the per-k nonlocal force loop and the per-k kinetic+nonlocal stress loop run through
`forces::accumulate_blocks`: with an active communicator each rank evaluates a contiguous k range, one allgatherv
of the per-k parts (3*atoms or 18 values per k plus one status value per rank), accumulation in global k order,
symmetrization after. Threads-only path is the former loop (bit-identical). Mac simulated ranks 2/3/7 on a 5-k
two-atom fixture: forces and frozen stress bit-identical to the replicated evaluation (difference 0). Also
found (not changed, ticket non-goal): the full stress always runs the Pulay assessment, three extra cold SCF
runs at c, 1.05c, 1.10c with the scan threshold min(1e-9 Ha, caller) and verification; on the Mac fixture
forces 0.003 s, frozen stress 0.009 s (it re-runs the force pass), scans 0.062 s. On B12 the scans are
expected to be the larger part of the 36 min; dev08 must report `stress.pulay_assessment.elapsed_seconds`
separately from `stress.elapsed_seconds` and `forces.elapsed_seconds`. Early record: results.json with
`status = derivatives_pending` is written before derivatives. Lesson: before distributing a post-SCF stage,
split its wall into kernels and hidden SCF runs.

### KEPTHX2 (2026-11-26, development protocol): REFSCF_KEPT_HX ADOPTED

Kept slabs now written inside the density accumulation loop with buffers reused across maps. dev08 MPI8 one run
per arm at e9250e4c, all verified: si64 6.37 -> 6.30 s (density 0.389 -> 0.529, solve 4.017 -> 3.829; 1452 of
5044 ordinary full-H columns replaced); si128 32.82 -> 31.40 s (solve 23.20 -> 21.06, density 2.12 -> 2.80; 3900
of 12261 columns replaced, -27%). Cost: peak RSS +52 MB/rank (si64), +271 MB/rank (si128; slab 194 MB, the rest
to be explained). Default on. Also found: dev08 mpi,mkl builds use baseline x86-64 (no target-cpu); QE is built
for the host. Evidence: dev08 /scr/mino/rustdft/kh2/runs/.

### Build flags (2026-11-26, dev08 si128 MPI8, one run per arm, two passes): no change

RUSTFLAGS none / -C target-cpu=x86-64-v3 / -C target-cpu=native at fa454494, results bit-identical across builds.
Light pass walls 31.34 / 30.67 / 32.06 s (rustfft), with MKL lanes 30.21 (v3) / 30.27 (native); profiled pass
30.42 / 30.49 / 31.65, MKL 31.88 / 30.20. No consistent winner; native auto-vectorized AVX-512 slows the
pure-Rust path (slab gather/scatter, kernels), consistent with the AVX-512 clock penalty seen in the MKL lane
study (MKL_ENABLE_INSTRUCTIONS=AVX2 removed the neighbour slowdown but made MKL kernels 40-55% slower).
Decision: keep the default build (no target-cpu). Single-run noise on si128 is about +-1 s.
Evidence: dev08 /data/sig-dft/bflag/si128/, /scr/mino/rustdft/mkllane/runs/si128/.

### Ladder retry with looser coarse stop (2026-11-26, dev08 MPI8, one run per arm at ac6aee9c): rejected

REFSCF_LADDER with REFSCF_LADDER_RATIO 1e4 / 1e5 (coarse stop at 1e-4 / 1e-3 Ry), KEPT_HX on, verify on.
si64: D 6.18 s, L1e4 7.71, L1e5 7.22; si128: D 31.79, L1e4 38.11, L1e5 37.48. Fine maps do not fall (si64 13 -> 14,
si128 17 -> 16) and fine full-H columns rise because the kept-HX warm start is lost at the handoff; coarse stage
costs 0.9-1.1 s (si64) and 4.9-6.7 s (si128). The reference driver's trajectory already matches QE's map count,
so a coarse pre-solve has little to remove. Ladder stays off.
Headline snapshot at the same binary (one run each): si64 verify off 4.57 s / verified 6.40 s / QE rmm 4.30 s
(1.06x / 1.49x); si128 23.59 / 30.15 / 21.57 s (1.09x / 1.40x). Evidence: dev08 /scr/mino/rustdft/lad2/runs/.

### P3-VERWARM (2026-11-26, development protocol): REFSCF_VERIFY_WARM ADOPTED

Verification and recovery maps start from the kept HX updated by dV; independent density, tolerances and the
Coulomb criterion unchanged. dev08 MPI8 one run per arm at 185489c: verification H columns si64 1050 -> 786,
si128 1942 -> 1422 (one band block per attempt), Coulomb errors identical, maps/attempts identical.
si128 verification 4.03 -> 3.72 s, wall 32.21 -> 30.36 s (most of it single-run spread); si64 6.14 -> 6.35 (noise).
Default on. Evidence: dev08 /scr/mino/rustdft/vw/runs/.

### PRECISE-2026-11-26: consolidated campaign at fa4cfadc

Full table: docs/performance/precise-2026-11-26.md. Verify-off / QE best 1.014-1.189 (MPI8 si216 1.014,
si128 1.098; 1 core si128 1.189); verified / QE best 1.074-1.648. 180 samples, no failures, accuracy classes held.

### P4 screen (2026-11-26, dev08 MPI8 verify off, one run per arm at db29dd7)

S = REFSTART_CHOLESKY + REFSTART_SPECIES_FILL + REFSTART_PARTITION_CACHE: setup si64 0.681 -> 0.682, si128 2.247 ->
2.126, si216 6.075 -> 5.374 s (atomic block 3.59 -> 2.84); wall si128 23.76 -> 23.50, si216 82.91 -> 81.55; rounding
class (|E - E_D| <= 1.1e-13). ADOPTED. K = REFSCF_KEPT_RECOMPUTE (bit-identical): peak RSS -45/-191/-505 MB per
rank, but solve +0.18/+1.57/+1.36 s; wall si64 +0.07, si128 +1.13, si216 -1.09 (density drop, unconfirmed). Kept
off as a memory-saving option. Evidence: dev08 /data/sig-dft/p4start/RESULT.md.

### GEN-2432530 (2026-11-25, dev08 one run per arm, MPI8 verify off; last timing screen before generalization)
REF_DAV_TRIANGLE alone (T) vs D: si128 wall 23.92 -> 23.66 s, si216 82.93 -> 81.11 s; si216 profiled pair: border_s
5.31 -> 4.50, border_h 4.64 -> 3.85 s (254 -> 180 GFLOP each), other stages within noise; maps/iterations/H vectors
identical; |E_T - E_D| <= 4.6e-13 Ha. ADOPTED (default on). STACK_ROTATE and SETUP_FUSED stay off. Memory breakdown
(D): si128 VmHWM 726 MB/rank, measured 506 MB (kept slabs 185, Davidson cached rows 71 + returned backing 48);
si216 1.91 GB/rank, measured 1.38 GB (kept slabs ~500, cached rows 202, returned 134). Generalization agreement:
G2 si32 PBE E_ref - E_legacy = -1.2e-10 Ha, G1 Mg NLCC +3.6e-10 Ha, both verified, both within the standard class
of the strict oracle. Evidence: dev08 /data/sig-dft/gen/RESULT.md.


## R10a legacy experiment retirement (2026-11-25)

User decision D3 and R10A-DELETE-EXPERIMENTS; removal, not a new speed claim.
Stored evidence remains under `tests/reference`. The reference ladder, live
fixed-map diagnostic, floor physics tests, and snapshot capture/replay remain.

| Removed switch | Prior evidence | Disposition |
| --- | --- | --- |
| `BAANDI_CUTOFF_LADDER_STANDARD` | LADDER-STANDARD / LADDER-S1-STANDARD | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_LADDER_STANDARD_SEED` | LADDER-STANDARD / LADDER-S1-STANDARD | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_LADDER_SECANT` | LADDER-SECANT / LADDER-S2-SECANT | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_LADDER_COARSE_STOP` | D1-LADDER-HANDOFF-KNOBS / F-LADDER-FRACTION | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_LADDER_COARSE_RULE` | D1-LADDER-HANDOFF-KNOBS | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_LADDER_FRACTION_RULE` | F-LADDER-FRACTION | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_LADDER_STOP_RULE` | F-LADDER-FRACTION | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_LADDER_WATCH` | F-LADDER-FRACTION | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_LADDER_BASIS_CHECK` | G-LADDER-FALLBACK | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_LADDER_FIRST_TOL` | D1-LADDER-HANDOFF-KNOBS | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_SCF_ADVANCE` | BT3-ADVANCE-CONTRACT / BT3-BUDGET-VARIANTS | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_SCF_ADVANCE_BLOCKS` | BT3-ADVANCE-CONTRACT / BT3-BUDGET-VARIANTS | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_SCF_ADVANCE_GROWTH` | BT3-GROWTH-TRIGGER | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_SCF_ADVANCE_ADAPT` | BT3-ADVANCE-CONTRACT / BT3-BUDGET-VARIANTS | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_SCF_ADVANCE_SCOPE` | BT3-ADVANCE-CONTRACT / BT3-BUDGET-VARIANTS | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_M2_OWNED_HANDOFF` | G-HANDOFF-OWNED-CLOSED | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_CUTOFF_LADDER` | C1-LADDER / C1-LADDER-DEV08 | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_VERIFY_PROFILE` | LANE-VERIFY Phase 0 | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |
| `BAANDI_FIRST_BLOCK` | BT0-FIRST-BLOCK | retired in R10a, code deleted; recreate from git history `d049b09a756c7da9f0b076698bd90aec6faec8bc` if needed |

### R10a deleted-test inventory

Surviving prolongation, snapshot, fixed-map and floor tests are retained.

| Source | Deleted test | Reason |
| --- | --- | --- |
| `crates/dft-cli/tests/all/run.rs` | `final_spectral_guard_is_additive_and_audits_only_reported_state` | tests deleted code |
| `crates/dft-engine/tests/owned_handoff.rs` | `profile_is_observational_at_active_and_terminal_handoffs` | tests deleted code |
| `crates/dft-engine/tests/owned_handoff.rs` | `profile_matrix_one_through_seven_ranks` | tests deleted code |
| `crates/dft-engine/tests/owned_handoff.rs` | `whole_handoff_uses_at_most_two_collectives_on_two_through_seven_ranks` | tests deleted code |
| `crates/dft-engine/tests/owned_handoff.rs` | `unsupported_consumers_keep_the_full_route` | tests deleted code |
| `crates/dft-engine/tests/owned_handoff.rs` | `rank_local_consumer_mismatch_is_an_agreed_error` | tests deleted code |
| `crates/dft-engine/tests/scf_advance.rs` | `scf_advance_keeps_the_standard_class_and_records_every_advanced_map` | tests deleted code |
| `crates/dft-engine/tests/first_block_diag.rs` | `physical_boundary_spectrum_and_criterion` | tests deleted code |
| `crates/dft-engine/tests/first_block_diag.rs` | `fixed_steps_preserve_production_bits_and_work` | tests deleted code |
| `crates/dft-engine/tests/all/direct_min.rs` | `default_method_is_the_ordinary_scf` | tests deleted code |
| `crates/dft-engine/tests/all/direct_min.rs` | `minimizer_descends_to_the_scf_state_and_hands_off_certified` | tests deleted code |
| `crates/dft-engine/tests/all/direct_min.rs` | `minimizer_is_thread_count_independent` | tests deleted code |
| `crates/dft-engine/tests/all/direct_min.rs` | `ineligible_problems_fall_back_with_a_reason` | tests deleted code |
| `crates/dft-engine/tests/all/direct_min.rs` | `track_a_stage0_bench` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `ladder_reaches_the_full_cutoff_fixed_point_skewed` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `ladder_reaches_the_full_cutoff_fixed_point_slab` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `ineligible_problem_falls_back_bit_identically` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `ladder_result_is_independent_of_thread_count` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `standard_ladder_is_opt_in_and_stays_in_the_standard_class` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `secant_seeds_keep_the_strict_fixed_point` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `standard_and_secant_ladder_is_independent_of_thread_count` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `d1_defaults_are_the_adopted_ladder` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `d1_knobs_change_only_their_quantity_and_the_fine_stage_converges` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `atomic_coarse_rule_only_raises_the_stop_level` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `ladder_reports_the_start_of_each_stage` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `d1_knobs_on_a_k_point_problem_fall_back_bit_identically` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `d1_knobs_are_independent_of_thread_count` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `ladder_d1_scan` | tests deleted code |
| `crates/dft-engine/tests/all/cutoff_ladder.rs` | `ladder_standard_secant_bench` | tests deleted code |
| `crates/dft-engine/src/gamma_real.rs` | `owned_permutation_matches_complex_injection_and_packing` | tests deleted code |
| `crates/dft-engine/src/gamma_real.rs` | `owned_rank_local_faults_are_agreed_before_transfer` | tests deleted code |
| `crates/dft-engine/src/gamma_real.rs` | `owned_missing_block_and_bad_density_are_agreed` | tests deleted code |
| `crates/dft-engine/src/gamma_real.rs` | `owned_overlap_boundary_uses_original_verdict` | tests deleted code |
| `crates/dft-engine/src/cutoff_ladder.rs` | `basis_plan_matches_constructor_at_dimension_boundaries` | tests deleted code |
| `crates/dft-engine/src/cutoff_ladder.rs` | `planned_basis_skip_preserves_full_cutoff_result` | tests deleted code |
| `crates/dft-engine/src/cutoff_ladder.rs` | `high_g_indicator_has_the_known_diagonal_value` | tests deleted code |
| `crates/dft-engine/src/cutoff_ladder.rs` | `fraction_rules_preserve_grid_and_parseval_budget` | tests deleted code |
| `crates/dft-engine/src/cutoff_ladder.rs` | `fraction_watch_is_bit_identical_and_forecast_is_diagnostic` | tests deleted code |
| `crates/dft-engine/src/cutoff_ladder.rs` | `measure` | tests deleted code |
| `crates/dft-engine/src/cutoff_ladder.rs` | `measure_2x2` | tests deleted code |
| `crates/dft-engine/src/replay_capture.rs` | `signed_mode_metric_has_analytic_units_and_sphere_filter` | tests deleted code |
| `crates/dft-engine/src/verification_tests.rs` | `verification_lane_profile_identity` | tests deleted code |
| `crates/dft-engine/src/verification_profile.rs` | `gather_preserves_unequal_rank_records_and_count_only_payloads` | tests deleted code |
| `crates/dft-engine/src/advance.rs` | `completion_rules_in_order` | tests deleted code |
| `crates/dft-engine/src/advance.rs` | `budget_adaptation_is_bounded` | tests deleted code |

### R10a strict-oracle exception (D2)

The stored strict Gamma baseline uses the legacy ladder by default. Removing
that default changed density bits despite an energy difference of only
9.77e-14 Ha. The user confirmed D2: keep the oracle's default arithmetic whole.
`cutoff_ladder.rs` therefore retains only that fixed-default initialization:
fraction 0.5, coarse density residual 1e-4, at most 30 coarse maps, no secant
memory, no tolerance seed, no fraction/basis/owned-handoff experiment branches.
It has no ladder switch reads. The switch retirement rows refer to the deleted
configuration and variant code, not this retained oracle default. The live
reference driver still uses its independent ladder. No baseline was changed.

### R10a final verification

Rust line counts include every `.rs` file recursively in the three affected
crates; inline tests are counted in `src`, and stored reference evidence is
excluded.

| Directory | Before files / lines | After files / lines |
| --- | ---: | ---: |
| `crates/dft-engine/src` | 73 / 49421 | 67 / 43196 |
| `crates/dft-engine/tests` | 87 / 38333 | 83 / 35917 |
| `crates/dft-cli/src` | 25 / 20846 | 25 / 20661 |
| `crates/dft-cli/tests` | 49 / 13957 | 49 / 13883 |
| `crates/dft-numerics/src` | 19 / 14774 | 19 / 14701 |
| `crates/dft-numerics/tests` | 18 / 4814 | 18 / 4792 |

- Switches: **187 → 168** (19 removed; retirement rows above).
- Test binaries: **112 → 109**. Named test inventory: **1274 → 1229**
  (45 removed, zero added); ignored inventory: **237 → 230**.
- Full runner summaries: before **1038 passed / 0 failed / 237 ignored**;
  after **1000 passed / 0 failed / 230 ignored**. Runner summary totals exceed the named
  inventory by one before and after; both measurements are reported separately.
- `cargo check --workspace --tests`: exit 0.
- `CHECK_SCCACHE=0 BAANDI_TEST_TMP=/tmp/r10a-test scripts/check.sh`: exit 0;
  build, tests, fmt, clippy all PASS. Sandbox denies the optional `ps` load
  probe; no verification gate is skipped or fails.
- `cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings`:
  exit 0.
- `cargo test -p dft-cli --test reference_identity -- --nocapture`: exit 0,
  one passed; both reference and strict golden values unchanged.
- `python3 -I scripts/make_readme_figures.py` and
  `python3 -I scripts/tests/test_readme_figures.py`: exit 0, nine tests passed.
  README block, figures, baseline files and `tests/reference` have no diff.
- The mode-comparison fixture now pins `gamma_real='off'` to compare direct
  sessions without the fixed-default oracle ladder; this replaces its retired
  ladder-off switch pin. The retired ladder routing fixture was deleted.

Both identity commands use `target/release/dft-cli`, three repeats, one thread,
`--cases gamma,k222s,k222 --require-identical`, and UPF SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Reference: `--baseline docs/performance/baseline-reference.json`, exit 0.
Oracle: `BAANDI_REFERENCE_STRICT=on` and
`--baseline docs/performance/baseline-oracle.json`, exit 0. Every energy bit and
density hash matches; all energy differences are exactly 0 Ha.

| Path | Case | Energy bits (Ha, binary64 hex) | Fine iterations | Density SHA-256 |
| --- | --- | --- | ---: | --- |
| reference | gamma | `c01d2e9f89e0349a` | 7 | `a0c76c8c7741125ad2cf59dee1534a0e4c29bc306002e6de22e56b92c99573e2` |
| reference | k222s | `c01fad7d54bc623e` | 7 | `67120c8ad65d1b86021aec1aa39c82e5e3dfd39f5dc758473589a9cd77c29d5c` |
| reference | k222 | `c01f532df8c633b9` | 7 | `680dd07126157dadfe5d2cb2cf036a3a05e5778e4afa8c99357651d7c3b3cdef` |
| strict oracle | gamma | `c01d2e9f89ffbffc` | 8 | `2fb29221af964b9ab95d635c1b23a0fade6f9b4676ab6a8d8fe65f290b29fcb5` |
| strict oracle | k222s | `c01fad7d54beed04` | 10 | `913aa5b4c39d9a49d0afc17d26efbc6f07a89cefd67c67fd1521287339baa564` |
| strict oracle | k222 | `c01f532df8db2672` | 10 | `4edeedb4469f5324df3da08c9d1e5810f7711367e7971c4cacb2921615a110fd` |

Raw command logs and inventories: `/tmp/r10a-evidence` (local run artifacts).
No timing acceptance or remote performance claim is made. Changes are uncommitted.

| Attempt | Evidence and result | Decision / lesson |
| --- | --- | --- |
| R10A-ORACLE-FIX (2026-11-25; explicit CLI off pin restored) | Root cause in `5b807dce`: `config.rs` retained/accepted `scf.cutoff_ladder="off"`, but `run.rs` replaced the old `None => problem.run_outcome()` dispatch with unconditional `cutoff_ladder::run_oracle`. The strict run therefore gained a half-cutoff coarse start; its eigen tolerance and Davidson algorithm were unchanged. Initial probes omitted this key and correctly found identical default-ladder arithmetic, missing the CLI routing regression. Fresh build of requested pre-R10A `a5ce9cf4` (isolated archive/target; git metadata is read-only) versus HEAD `10dd123a`, pure Rust 1.99.1 macOS/aarch64, strict CLI with `driver=legacy`, atomic density, Pulay beta 0.7/history 8, 20 bands, 200 maps, one thread, cubic Si8 a=10.261212856717933 Bohr, Gamma, 30 Ry, 36^3, Si.pz-vbc SHA-256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`: pre energy `c03f5261b88629b1` (-31.321803601030926 Ha), 10 maps, Davidson `[5,5,4,6,10,9,8,8,10,10]`; broken HEAD `c03f5261b8862afa`, 9 maps, `[3,5,5,6,5,6,6,6,6]`, dE=-1.1688428003253648e-12 Ha. Removing the pin on the old CLI reproduces the broken density exactly. Fixed off-pin energy, every density byte, eigen-tolerance history and Davidson list equal pre-R10A; density SHA-256 `05371fad92da2d2a0405953484c2e55bc6c9fb1966b30acb0a6a8f40b93e55eb`. Omitted-key result also remains identical. Atom 1 displaced by fractional `[0.001,-0.0005,0.0007]`: pre/fixed energy `c03f5261389bd4e2`, 11 maps, identical density and force bits (max dF=0 Ha/Bohr). Full CLI regression `run::tests::strict_cli_ladder_off_preserves_oracle` fails on the broken dispatch and passes after the fix; it checks both off and omitted inputs against engine entry points on every backend, plus pre-R10A macOS/aarch64 goldens: `bfe3f626c809494f` Ha bits, 11 maps, Davidson `[3,2,3,2,2,3,1,3,3,3,2]`, density SHA-256 `310f2763cbc5f728e40f510460c39eb8ed19326a4f2484ee0b7351285fe84a1e`; debug runtime 0.16 s with golden checks. Final verification: `cargo check --workspace --tests`, CLI regression, engine `oracle_identity`, `CHECK_SCCACHE=0 BAANDI_TEST_TMP=/tmp/r10a-oracle-tests scripts/check.sh`, and `cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings` all exit 0; full gate 1014 passed / 0 failed / 230 ignored (111 binaries). Both `bench_scf.py --require-identical` baselines (Gamma/k222s/k222, three repeats, one thread) exit 0, dE=0 Ha and identical density hashes. README regeneration and nine freshness tests exit 0 with no generated-file diff. Raw inputs, build/run logs and comparison JSON: `/tmp/r10a-oracle-resume`; invocation `BAANDI_REFERENCE_STRICT=on <pre/head/fixed>-cli run evidence/si8.toml --output-dir evidence/<arm>`; input SHA-256 `01665256347f3b347c808e416ad726f9aa13446ca0004d6ed10e2ce99141271c`. Binary SHA-256 pre/head/fixed: `478f9ff8b62be1ddd47d21d0ea528435250e595ab75eaa8ef233ca5e1754a0e9` / `f3f0b1d66e5b4660e3d461b434363964ac05e13f818092a87595c6c834114329` / `44d1ee726f501097b8d8bfaa58559224cea91f43f5860f10ca261440f8229b96`. | Honor the existing off-only input pin; when omitted, keep the fixed-default oracle ladder. No retired experiment, environment switch, tolerance or baseline is restored/changed. SCF_ADVANCE and FIRST_BLOCK were off and are not the cause. R6A/R8C/R9A need no arithmetic change. Keep the earlier engine ladder golden; its coverage is the omitted-key path. Relayed dev08 evidence `/data/sig-dft/bisect20/RESULT.md` identifies R10A with and without MKL, including np1; local bit patterns are platform-specific, and the fixed binary still needs dev08 RG-ORC confirmation through the orchestrator. |

## R11a reference switch retirement (2026-11-25)

D3 retirement only; no new timing/adoption claim. Default arithmetic retained.
Complete reader/retention inventory and validation: [switch-retirement.md](switch-retirement.md).

| Switch | Previous evidence | Disposition |
| --- | --- | --- |
| `REF_DAV_RESTART_CORRECTIONS` | P3 timed (line 3539) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_MKL` | REFSTART_MKL (2026-11-24): ADOPTED (line 3561) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_WIDE_ROTATE` | P3 timed (line 3539) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_VERIFY_GATE` | I25-VERIFY-COST (line 3503) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_SETUP_PHASE_REUSE` | P3 timed (line 3539) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_SINGLE_GEMM` | REF_DAV_SINGLE_GEMM provisionally ADOPTED (2026-11-25, development protocol) (line 3578) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_SWAP_BUFFERS` | P3i (line 3598) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_LEAN_RESET` | P3i (line 3598) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_DIAGONAL_REUSE` | P3i (line 3598) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_ORBITAL_REUSE` | P3i (line 3598) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_SQRT` | P3j (line 3611) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_FUSED_APPEND` | P3j (line 3611) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_CHECK_APPEND` | P3j (line 3611) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_REUSE_ROWS` | P3j (line 3611) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_GENERALIZED` | P3l (line 3620) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_COMPACT_G` | P3l (line 3620) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_LAZY_OUTPUT` | P3l (line 3620) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_IDEAL_RESET` | P3m (line 3628) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_LEAN_ENERGY` | P3q (line 3638) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_FFT_SHARE` | P3t (line 3659) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_SINGLE_GEMM` | P3t (line 3659) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_KEPT_HX` | KEPTHX2 (2026-11-26, development protocol): REFSCF_KEPT_HX ADOPTED (line 3682) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_VERIFY_WARM` | P3-VERWARM (2026-11-26, development protocol): REFSCF_VERIFY_WARM ADOPTED (line 3711) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_CHOLESKY` | P4 screen (line 3724) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_SPECIES_FILL` | P4 screen (line 3724) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_PARTITION_CACHE` | P4 screen (line 3724) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_TRIANGLE` | GEN-2432530 (line 3732) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_SPIN_STIFFNESS` | R6a / G4B 48/48, dev08 gen15; ADR 0045 (manifest adoption record; no separate ledger row) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_HISTORY_CONDITIONING` | R6a / G4B 48/48, dev08 gen15; ADR 0045 (manifest adoption record; no separate ledger row) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_VERIFY_REUSE` | I25-VERIFY-COST (line 3503) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_VERIFY_GATE_FIRST` | I25-FIRST-ATTEMPT (line 3525) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSCF_VERIFY_ENERGY` | P3f/P3g screens (2026-11-25, development protocol, one run each, dev08 MPI8) (line 3588) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REFSTART_LOCAL_SHELLS` | P3t (line 3659) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_STACK_ROTATE` | GEN-2432530 (line 3732) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |
| `REF_DAV_SETUP_FUSED` | GEN-2432530 (line 3732) | retired in R11a, recreate from git `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` if needed |


## R10BC-LEGACY-STANDARD (2026-11-25): legacy standard switch retirement

The two CLI routes are Reference and StrictOracle; the negative G4D qualification
routes to StrictOracle with the reason recorded. L9 input keys are unknown keys.
This is deletion work, not a performance adoption claim. No dev08 timing was run.

| Removed registry switch | Reason |
| --- | --- |
| `REFERENCE_SCF` | Reference is the production route; strict oracle has its own switch. |
| `STANDARD_STOP` | The legacy standard stop rule has been deleted from the oracle session. |
| `VERIFY_REUSE` | Legacy standard verification retry/reuse has been deleted. |
| `STOP_METRIC` | Legacy standard Coulomb stopping estimates have been deleted. |
| `VERIFY_SCHEDULE` | Legacy standard verification scheduling has been deleted. |
| `STOP_PREDICTOR` | Legacy standard cross-estimate verification has been deleted. |
| `REUSE_PSI` | Standard-only retained real-space start and H-X update deleted. |
| `M2_REUSE` | Standard-only distributed kept-start solver branch deleted. |
| `M2_LAZY` | Standard-only lazy solve return/start branch deleted; strict solves retain full delivered orbitals. |
| `M2_OWNED_INPUT` | Standard-only direct owned-row input packing branch deleted. |
| `SPIN_SHARE` | Standard-only shared spin-channel solve deleted; strict channels retain independent solves. |
| `SERIAL_DENSITY_FUSED` | Standard-only direct density sink deleted; strict per-band certificate grids retained. |

The registry shrinks from 133 to 121 entries (127 at the continuation start).
`STANDARD_RESIDUAL_GUARD` stays because `refscf/electronic.rs` consumes it.
`CERT_SCREEN` stays because the strict oracle consumes it; deleting that screen
would change the protected arithmetic. Shared reference lazy-result and owned-grid
switches remain because the reference driver consumes them.

Mac identity evidence from the checked-in benchmark script: release
Gamma/k222s/k222 against `baseline-reference.json` and `baseline-oracle.json`
each passes `--require-identical`: delta energy 0 Ha and identical density hashes,
reference iteration counts 7/7/7 and oracle counts 8/10/10. The script now requires
only the retained strict switch metadata and recognizes historical oracle labels
without changing any stored numbers. Ten benchmark metadata tests pass.
Local records: `/tmp/r10bc-cont/{reference,oracle}-approved.json` and
`/tmp/r10bc-cont/bench-{reference,oracle}-approved.log`. No tolerance was relaxed.

Final continuation test pass: 103 binaries, 953 passed, 0 failed, 217 ignored.
Seven native Python binding tests pass, including same-build Si2 Gamma/k222
CLI identity. The README generator and its nine tests pass without generated
changes. After the approved shared-orphan cleanup, the full gate passes:
build, tests, formatting, workspace clippy, production-only check and
production-only clippy. Command: `CHECK_SCCACHE=0 BAANDI_TEST_TMP=/tmp/r10bc-cont scripts/check.sh`; exit 0, log `/tmp/r10bc-cont/full-check-approved.log`.
The remaining clippy cleanup simplified the mixing-beta validation and the
strict-only Davidson profile test, and removed detached legacy test comments.

Continuation inventory (base `ca846fd6f3371052b5c8b8399518d5c8db87365a`):
all `crates/**/*.rs` excluding target, 205595 -> 196470 physical lines;
1231 -> 1180 literal test attributes. Nine obsolete engine test targets and
the 1914-line legacy Davidson example are removed. Each retired test and its
reason is listed in `/tmp/r10bc-cont/test-migration.md`; no comparison tolerance
was relaxed. New changes remain uncommitted for review.

## R11B-ORACLE-SWITCHES — D3 oracle switch retirement (2026-11-25)

Code retirement at frozen defaults, not an optimization or timing claim. Source
commit `1535810f4d6abffbe36483dd82c174ece6664b61`; readers, retained boundaries and verification are recorded in
[switch-retirement.md](switch-retirement.md#r11b-freeze-oracle-only-switches-d3).
No QE source or new QE comparison was used.

| Switch | Disposition | Evidence / lesson |
| --- | --- | --- |
| `INIT_OWNED_STATES` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `CERT_SCREEN` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `POTENTIAL_REUSE` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `MIX_GROWTH_RESET` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `MIX_CELL_KERKER` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `EMPTY_LOOSE` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `BUFFER_TOL` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `DENSITY_PAIR` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `DENSITY_REDUCE` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `NL_ENERGY_REAL` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `ATOMIC_LOCAL_BLOCK` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `ATOMIC_REAL` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `PHASE_THREADS` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `MIX_COULOMB_GUARD` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_LIGHT_STATE` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_LAZY_SERIAL` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `IMPORT_CHECK` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_DENSITY` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_ORTHO` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_LEAN` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_EIGEN` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_B5` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_WORKSPACE` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_SELECTED_ROOTS` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_PANEL_GS` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_START` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_LEAN_ALGEBRA` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `M2_LAYOUT_CACHE` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `B1_NARROW` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `B1_MARGIN` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `B1_ROWS` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `B6` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `SECOND_VECTOR` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `MGS` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `ORTHO_BLOCK` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `WIDE_DOT` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `PROJECT_WIDE` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `DAVIDSON_ROTATE_COLUMNS` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `DAVIDSON_UPPER_BORDER` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `LEAN_BORDER` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `DAVIDSON_PACKED_REDUCE` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `MIX_INCREMENTAL` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `MKL_NONLOCAL` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `APPLY_LEAN` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `HAPPLY_OVERLAP` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |
| `HAPPLY_FOLD` | Retired in R11b; former default frozen. | Oracle-only caller audit; see R11b verification. Preserve shared consumers and numerical fallbacks. |

### R11b part 2 — frozen oracle controls (2026-11-25)

No speed claim; local oracle/reference baseline identity passes 3 cases each, all seven recorded fields exact. See [reader audit and gates](switch-retirement.md#r11b-part-2-remaining-oracle-arithmetic-controls-2026-11-25).

| Switch | Outcome | Evidence and lesson |
| --- | --- | --- |
| `EIGEN_FACTOR` | Retired; fixed 0.1. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `MPI_M2` | Retired; fixed on. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `M2_BATCH` | Retired; fixed on. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `M2_CAPACITY` | Retired; fixed auto (3 bands at >=4 ranks; otherwise 6). | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `M2_MIX` | Retired; fixed on. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `B1` | Retired; fixed on. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `B1_SKETCH` | Retired; fixed count. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `B5` | Retired; fixed on (serial). | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `RITZ_FUSED` | Retired; fixed off. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `NARROW_PROJ` | Retired; fixed off. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `NARROW_REDUCED` | Retired; fixed off. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `ORTHO` | Retired; fixed DGKS 0.1. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `DAVIDSON` | Retired; fixed lean for Ritz-change; orthonormal for strict. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `LEAN_PROJ_CACHE` | Retired; fixed off. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `RESTART_KEEP` | Retired; fixed off. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |
| `EWALD_CACHE` | Retired; fixed on. | CLI metadata readers removed or made constant; oracle-only arithmetic frozen. Both baselines exact (delta energy 0 Ha). |

## GSPEED-PLAN — general (k, spin/FD) path speed analysis (2026-11-26; analysis only, no code)

Trigger (relayed dev08 MPI8, same inputs as QE 7.6): IrI3 QE 132 s vs BAANDI 1194 s (SCF+forces+stress, gen41),
GeTe 260 k QE 43 s vs BAANDI SCF 378 s, B12 172 k QE 84 s vs BAANDI SCF 560 s. Plan and file:line map:
docs/performance/general-path-speed-plan.md. Code findings (exact): (1) no k pools, every k is G-distributed over all
ranks and solved sequentially (~17 collectives per complex Davidson iteration, ~10 of them 1-double error agreements);
(2) k atomic start evaluates the radial quadrature per (row, atom, m) instead of per (species channel, |k+G|), and
orthonormalizes by column MGS with one tree sum per column pair; spin runs build the whole k start twice;
(3) complex Davidson recomputes full S/H Grams every iteration (real path: border + triangle) and allocates its
workspace per call; (4) per-k per-map kinetic/nonlocal energy terms and per-k potential/density plane maps;
(5) collinear spin disables time reversal and symmetry reduction (full mesh). Mac observation (O, release, 20000
calls, 1534-point mesh): radial quadrature 12.7/17.0/18.2 us per call for l = 0/1/2; estimate (E) B12 atomic-start
fill 80-260 s of 560 s, removable by a bit-identical memo. No speed claim; ranking waits for the kprof phase profiles
and a QE `-nk 1/2/4/8` scan. Lesson: before generalizing a fast path, list which Gamma optimizations the new path
inherits; the k start and the complex solver inherited none of the P3/P4 start and Davidson work.
