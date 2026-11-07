# Thread width analysis: why T=8 uses about two of eight threads (2026-11-07)

Owner: perf-optimizer. Mac observations only (Apple M6, 12 logical CPUs: 2 Super + 4 Performance + 6 Efficiency cores, macOS 27.0.1; a thread count above 2 or 6 starts using slower or different cores, so Mac scaling is not an Amdahl experiment on identical cores and is never acceptance). dev08 numbers are quoted as relayed (`tests/reference/perf/dev08_d0_profile_relayed.md`, `dev08_materials_qe_relayed.md`). No repository source was edited: all instrumented and patched builds live in the session scratchpad (paths in section 7). Ledger rows THR-1 to THR-5 and backlog items B21 to B25 record what was run.

Build under test: release `dft-cli` from the working tree at HEAD `fbd802b` (crates unchanged at `4df5358`), sha256 `2f9caa40...d703956`, Si UPF sha256 `e8d93375...5d68d5217`, inputs from `scripts/bench_scf.py case_input` (si8, k222, k222s, k444s) and `benchmarks/cases/si16-gamma.toml` through `scripts/bench_suite.py case_input`, strict profile, `gamma_real` off or on, `threads = T`, `[symmetry] mode = "auto"` where stated.

## 1. Result in twelve lines

1. Mac wall T=1 to T=8 (5 to 7 runs each, spread at most 2.8 %): si8 real 1.195 to 0.751 s (1.59x, T=4 already 0.751), si8 complex 2.396 to 1.096 (2.19x), si16 real 5.239 to 2.624 (2.00x), si16 complex 8.969 to 3.240 (2.77x). dev08 (relayed) is lower in every cell (si8 real 1.20x, complex 1.54x; si16 1.59 / 2.04; si32 1.93 / 2.36).
2. Amdahl fit `wall = a + b/T` over T=1,2,4,8: effective serial fraction f = 0.545 (si8 real), 0.362 (si8 complex), 0.414 (si16 real), 0.260 (si16 complex); implied ceilings 1.83x, 2.76x, 2.42x, 3.84x.
3. The phase-sum (instrumented, section 3) splits f for si8 real into 0.29 of work that is plainly serial (EVD, in-block MGS, nonlocal GEMM, per-SCF real-Hamiltonian build, audit, potential, mixing, setup) and 0.26 more that is parallel code with poor efficiency (the 12-panel GEMMs and narrow H blocks).
4. At T=8 the owner thread of the solve is 82 to 91 % active and the other workers 16 to 41 %: mean busy threads per run 2.2 (si8 real), 2.7 (si16 real), 3.0 (si8 complex), 3.65 (si16 complex), from 6 / 4 / 6 / 4 `sample` runs at 1 ms (5.6k to 22k active samples per cell); this reproduces dev08 (about two busy of 17).
5. The width limit is structural, not dispatch latency: in the real path 177 of 293 H calls carry at most 4 vectors (1 or 2 FFT pairs), the GEMMs have 12 fixed 256-row panels and run serially below 1,048,576 multiply-adds, and EVD, in-block MGS and the nonlocal GEMM are single-threaded.
6. Refuted hypothesis (THR-3): Rayon workers going to sleep between regions (sleep after about 33 yields) costs the time. A patched `rayon-core` that spins 4000 rounds instead of 32 gave no gain (si8 real T=8 0.7526 to 0.7913 s, T=4 0.7548 to 0.7766 s, complex T=8 1.115 to 1.200 s). A wake microbenchmark does show 35 to 87 µs extra per 12-piece region at T=8 when it follows at least 20 µs of serial work (7 to 10 µs at T=4), so it is not zero, but removing it by spinning does not help here (consistent with A47).
7. Largest serial items at T=8 in si8 real (0.7526 s): EVD 9.0 %, outside-solve serial phases 13.2 % (audit 6.4 %, potential 2.6 %, mixing 1.9 %, setup 1.3 %, process 1.0 %), in-block MGS 6.8 %, nonlocal GEMM 6.0 %, real-Hamiltonian rebuild per SCF iteration 4.9 % (si16: 4.4 %, plus 1.8 % solver setup), projection overlap GEMMs 8.2 %, Ritz rotation 13.9 %.
8. Real Hamiltonian is rebuilt at every SCF iteration (`from_complex`, 16 times): 0.036 s of 1.17 s at T=1 and 0.037 s of 0.753 s at T=8, serial; it allocates ceil(capacity/2) = 60 FFT pair plans per solve (memset share on dev08 is 7 %).
9. k points: k444s with symmetry (10 irreducible): T=1 0.752 s, T=10 0.1965 s (3.83x), T=12 0.1968, T=16 0.1983: no gain above T=10 because outer width is min(T, nk) and inner = floor(T/outer) = 1; 38 to 41 % of the T>=10 wall is serial (CLI symmetry detection 12 %, setup, potential, energy, mixing).
10. Giving every k solve the whole pool as inner budget (one line, scratch build, THR-4) cut k444s symmetric wall by 4.9 % (T=4), 10.9 % (T=8), 7.6 % (T=12), 6.0 % (T=16) with density sha256 identical to the T=1 default at T=8 and T=16.
11. k222 with symmetry (3 irreducible points) is already limited by the serial remainder (60 % of the T=8 wall: setup 0.022, potential 0.020, energy 0.009, mixing 0.005, outside-engine 0.023 s); `symmetry_map` is evaluated at seven call sites (config.rs 85, 147, 155, 1020, 1041, preflight.rs 518, run.rs 866).
12. Top three proposals (section 5): P1 persistent real Hamiltonian and solver state across SCF iterations (bit-identical, E: 4 to 6 % at T=8, 3 to 4 % at T=1); P2 overlap of independent phases (nonlocal image concurrent with local FFT pairs, audit of iteration i concurrent with iteration i+1, one region for the V and HV Ritz rotations) (bit-identical, E: 10 to 20 % at T=8, 0 at T=1); P3 full-pool inner budget for k solves plus one symmetry map (bit-identical, measured 6 to 11 % at 10 k).

## 2. Wall-time scaling and Amdahl (item 2 of the task)

Protocol: alternating T order inside each repetition (T=1,2,4,8 repeated), one discarded warm-up run per T, wall = `perf_counter` around the CLI process. Host load average (1 min) was 2.5 to 9.7, mostly the runs themselves; the machine had other users at about 1.4 before the session. Speed-up = median(T=1) / median(T). Karp-Flatt f = (1/S - 1/T) / (1 - 1/T).

| Case (strict) | runs | T=1 (s) | T=2 | T=4 | T=8 | S(8) | Karp-Flatt f at T=2 / 4 / 8 | fit f | 1/f |
|---|---|---|---|---|---|---|---|---|---|
| si8 real | 7 | 1.1949 | 0.8541 | 0.7510 | 0.7508 | 1.59 | 0.430 / 0.505 / 0.575 | 0.545 | 1.83 |
| si8 complex | 7 | 2.3958 | 1.5635 | 1.2132 | 1.0963 | 2.19 | 0.305 / 0.342 / 0.380 | 0.362 | 2.76 |
| si16 real | 5 | 5.2393 | 3.4885 | 2.8813 | 2.6238 | 2.00 | 0.332 / 0.400 / 0.429 | 0.414 | 2.42 |
| si16 complex | 5 | 8.9690 | 5.4634 | 3.9500 | 3.2400 | 2.77 | 0.218 / 0.254 / 0.270 | 0.260 | 3.84 |

Min-max spread of every cell: 0.2 to 2.8 % (largest: si8 real T=8, 0.7323 to 0.7537 s). A fit through T<=4 only underpredicts T=8 wall by 0.6 to 13 % (si8 real: predicted 0.656, measured 0.751), so the Mac loses further above 4 threads (Efficiency cores; section 1 caveat).

Serial constants outside the diagonalization (identical at every T, medians, seconds): si8 potential 0.0190, mixing 0.0141, energy audit 0.0483, setup 0.0085, process start/output 0.0076 (sum 0.0975 = 8.2 % of T=1 real wall and 13.0 % of T=8 real wall); si16 0.0352 + 0.0258 + 0.2380 + 0.0218 + 0.0125 = 0.333 s (6.4 % of T=1 and 12.7 % of T=8 real wall). Density is parallel (si8 0.0509 to 0.0177 s).

Phase-sum Amdahl (si8 real; times from section 3): plainly serial work at T=1 is 0.34 s of 1.170 s (f_hard = 0.29); if all parallel phases were perfect at T=8 the wall would be 0.34 + 0.837/8 = 0.445 s (2.6x), against 0.753 measured. The remaining 0.31 s is the inefficiency of the parallel regions. The ceiling with the present serial work is 1/f_hard = 3.4x even with infinite perfect threads.

## 3. Per-phase wall time (instrumented scratch build, items 1 and 3)

Method: wall-clock guards (`Instant::now`, 40 ns) at every Davidson phase and in `solve_real_gamma`, accumulated over the whole SCF in static atomics on the owning thread (so each phase value is wall time of that phase seen by the solve), printed at the end. Overhead check: instrumented walls 1.166 / 0.857 / 0.755 / 0.753 s at T=1/2/4/8 against 1.195 / 0.854 / 0.751 / 0.751 uninstrumented. 5 runs per T for si8 and 3 for si16 (counts below are per SCF, identical across runs; medians shown). Counts: 16 real solves, 291 reduced-matrix calls (Davidson steps including the certificate pass), 293 H calls, 275 projection calls (si16: 15, 318, 320, 303).

### 3.1 si8 real strict (seconds; T1/T8 = speed-up of the phase)

| Phase | code | T=1 | T=2 | T=4 | T=8 | T1/T8 | per iteration µs (T=1 to T=8) |
|---|---|---:|---:|---:|---:|---:|---|
| solve total | `gamma_real.rs:86` | 1.0217 | 0.7261 | 0.6339 | 0.6358 | 1.61 | |
| H local FFT pairs | `fused_pair.rs:397-437` | 0.3940 | 0.2283 | 0.1653 | 0.1401 | 2.81 | 1354 to 481 |
| H nonlocal (GEMM, couplings) | `real_gamma.rs:120-164` | 0.0564 | 0.0469 | 0.0456 | 0.0453 | 1.24 | 194 to 156 |
| H kinetic add | `real_gamma.rs:234` | 0.0031 | 0.0032 | 0.0036 | 0.0041 | 0.75 | |
| reduced matrix (Gram) | `real_davidson.rs:523` | 0.0560 | 0.0403 | 0.0366 | 0.0365 | 1.54 | 192 to 125 |
| EVD | `real.rs:395` | 0.0652 | 0.0655 | 0.0663 | 0.0676 | 0.96 | 224 to 232 |
| Ritz rotation (V and HV) | `real_davidson.rs:493` | 0.1736 | 0.1100 | 0.0914 | 0.1049 | 1.65 | 597 to 361 |
| residual | `real_davidson.rs:579,627` | 0.0064 | 0.0047 | 0.0049 | 0.0080 | 0.81 | |
| corrections (preconditioner) | `real_davidson.rs:558` | 0.0037 | 0.0028 | 0.0024 | 0.0036 | 1.05 | |
| projection total | `real_davidson.rs:700` | 0.2102 | 0.1710 | 0.1645 | 0.1692 | 1.24 | 722 to 581 |
| of which overlap Gram | | 0.0927 | 0.0668 | 0.0623 | 0.0617 | 1.50 | 319 to 212 |
| of which projection rotate | | 0.0702 | 0.0562 | 0.0532 | 0.0553 | 1.27 | 241 to 190 |
| of which in-block MGS (serial) | `real_davidson.rs:760-779` | 0.0482 | 0.0485 | 0.0496 | 0.0514 | 0.94 | 166 to 177 |
| real Hamiltonian build per SCF | `gamma_real.rs:126` | 0.0359 | 0.0363 | 0.0363 | 0.0371 | 0.97 | |
| solver setup per solve | `real_davidson.rs:111-` | 0.0097 | 0.0095 | 0.0094 | 0.0106 | 0.91 | |
| append, restart, unpack | | 0.0030 | 0.0028 | 0.0028 | 0.0033 | | |
| **engine wall** | | 1.1698 | 0.8574 | 0.7548 | 0.7526 | 1.55 | |

Outside the solve at T=8: potential 0.0193, mixing 0.0144, energy audit 0.0485, setup 0.0095, process 0.0077, density 0.0177 (parallel).

T=8 composition of 0.7526 s: hard serial inside the solve 0.237 s (31.5 %: EVD, MGS, nonlocal, `from_complex`, setup, small items), serial phases outside the solve 0.099 s (13.2 %), density 0.018 (2.4 %), 12-panel GEMM phases 0.259 s (34.3 %, 1.5x on 8 threads), local FFT pairs 0.140 s (18.6 %, 2.8x).

### 3.2 si16 real strict (3 runs; wall 5.068 / 3.504 / 2.882 / 2.639 s at T=1/2/4/8)

solve 4.556 to 2.253; H local 1.264 to 0.394 (3.2x); H nonlocal 0.305 to 0.175 (1.75x); reduced 0.292 to 0.126; EVD 0.254 to 0.257 (serial, 9.7 % of the T=8 wall); Ritz rotation 1.100 to 0.375; projection 1.102 to 0.701 (Gram 0.445 to 0.206, rotate 0.354 to 0.189, in-block MGS 0.320 to 0.323: 12.2 % of the T=8 wall, serial); `from_complex` 0.116 constant (4.4 %); solver setup 0.048 (1.8 %); energy audit 0.238 constant (9.0 % of the T=8 wall).

### 3.3 Complex path (counts 345 steps for si8, 232 for si16; wall seconds)

| Phase | si8 T=1 | si8 T=8 | T1/T8 | si16 T=1 | si16 T=8 | T1/T8 |
|---|---:|---:|---:|---:|---:|---:|
| H apply | 0.9598 | 0.2783 | 3.45 | 3.0025 | 0.6852 | 4.38 |
| reduced H | 0.2054 | 0.1082 | 1.90 | 0.7949 | 0.2840 | 2.80 |
| EVD (serial) | 0.1239 | 0.1284 | 0.96 | 0.4246 | 0.4271 | 0.99 |
| Ritz rotation | 0.4589 | 0.1686 | 2.72 | 2.1232 | 0.5259 | 4.04 |
| residual loops (serial, `lib.rs:957`) | 0.0336 | 0.0341 | 0.99 | 0.0815 | 0.0816 | 1.00 |
| preconditioner (serial) | 0.0129 | 0.0130 | 0.99 | 0.0333 | 0.0333 | 1.00 |
| orthonormalize (MGS + projection) | 0.4123 | 0.2328 | 1.77 | 1.8418 | 0.7507 | 2.45 |
| restart | 0.0045 | 0.0046 | 0.99 | 0.0131 | 0.0131 | 1.00 |
| wall | 2.4045 | 1.1154 | 2.16 | 8.9649 | 3.2548 | 2.75 |

The complex path has wider H blocks (single bands, no pairing) and no per-SCF rebuild, so its serial fraction is lower; its serial or poorly parallel items at T=8 are EVD (11.5 % of si8 wall, 13.1 % of si16), orthonormalize (20.9 %, 23.1 %), residual plus preconditioner plus restart (4.7 %, 3.9 %).

### 3.4 H application by block width (real path; seconds, T=1 to T=8; calls per solve set)

| vectors per call | pairs (pieces) | calls | si8 T=1 | si8 T=8 | T1/T8 | si16 T1/T8 |
|---|---|---:|---:|---:|---:|---:|
| 1 to 2 | 1 | 94 | 0.0405 | 0.0409 | 0.99 | 0.99 |
| 3 to 4 | 2 | 83 | 0.0676 | 0.0417 | 1.62 | 1.59 |
| 5 to 8 | 3 to 4 | 31 | 0.0440 | 0.0211 | 2.09 | 2.02 |
| 9 to 16 | 5 to 8 | 16 | 0.0429 | 0.0136 | 3.16 | 2.97 |
| 17 to 20 | 9 to 10 | 69 | 0.2614 | 0.0751 | 3.48 | 3.75 |

At T=8 the calls with at most 8 vectors cost 0.104 s of 0.192 s of H time (54 %, 14 % of the engine wall) and speed up at most 1x to 2x: that is the width limit (60 % of the H calls are tail steps with one to four unconverged bands). The complex path has the same shape (calls of at most 4 vectors: 220 of 345; T1/T8 1.2x to 2.6x). The wide calls reach 3.5x on 8 threads because 10 pairs on 8 workers need two rounds (ceil(10/8) = 2: limit 5x).

## 4. Per-thread activity from `sample` (item 1)

Method: process started suspended through `posix_spawn(START_SUSPENDED)`, `sample <pid> 120 1 -mayDie` attached, then continued (the earlier attempt with `sh -c 'kill -STOP $$; exec ...'` attached to the shell and was discarded; the first runs of the earlier PROF-1 lost 0.1 s). 1 ms interval requested; `sample` records about 0.65 samples per requested millisecond (si8 real T=8 run of 0.91 s wall gives about 600 samples per thread, one sample per 1.5 ms), so shares are shares of sampled time and the sample counts below are not milliseconds. Thread totals include 9 threads (8 workers and the queue-main thread that blocks in the install latch). A bucket is quoted only with at least 30 samples; `*` marks fewer. Mangled Rust v0 names were bucketed by substring of the innermost phase frame; work stolen by another worker carries no caller frames, so stolen GEMM tiles form their own bucket ("caller n/a").

| Config (T=8) | runs | thread-samples | idle share | active samples | mean busy threads | owner active | other workers active |
|---|---:|---:|---:|---:|---:|---:|---|
| si8 real | 6 | 32,747 | 75.4 % | 8,042 | 2.21 | 0.82 | 0.16 to 0.23 |
| si8 complex | 6 | 48,010 | 66.8 % | 15,933 | 2.99 | 0.83 | 0.27 to 0.33 |
| si16 real | 4 | 75,537 | 70.2 % | 22,532 | 2.69 | 0.89 | 0.24 to 0.26 |
| si16 complex | 4 | 94,584 | 59.5 % | 38,349 | 3.65 | 0.91 | 0.37 to 0.41 |

Shares of active samples (T=8), same configs in order: H apply (all workers) 33.7 / 40.2 / 28.4 / 37.3 %; stolen GEMM tiles 25.4 / 14.4 / 39.4 / 18.1 %; projection incl. MGS 9.7 / 8.7 / 10.1 / 10.8 % (MGS and glue 563 / 985 / 1719 / 3263 samples, 98 / 74 / 98 / 60 % on the owner thread); EVD 3.9 / 3.7 / 3.7 / 3.6 % (310 / 597 / 836 / 1390 samples, 100 % owner); Ritz rotation 4.4 / 20.5 / 4.8 / 22.0 %; reduced matrix 2.2 / 2.7 / 1.7 / 2.2 %; real-Hamiltonian build 2.7 % (218 samples, 92 % owner) and 1.9 % for si16 real; spinning in `find_work` and `steal` (no phase) 5.8 / 3.7 / 2.1 / 1.8 %; density 4.6 / 1.7 / 3.5 / 1.8 %; potential, mixing, setup, audit, CLI each below 2 % (potential 156 samples, mixing 66, audit 22 *). T=1 si8 real (3 runs, 2,867 samples): H apply 36.4 %, projection 17.6 % (MGS 8.6 %), Ritz rotation 14.5 %, EVD 6.4 %, reduced matrix 4.6 %, real-H build 3.7 %, density 4.2 %. These agree with the dev08 D0 profile (local FFT 31.9 %, algebra 48.3 % at T=1; T=8 local FFT 30.7 %, algebra 57.1 %, waiting share 84.4 %).

## 5. Serial or poorly parallel regions: location, limit, proposal

Per Davidson iteration (T=8, si8 real) the solve passes 12 to 16 regions (`pool.install` or `WorkerPool::join` trees): H local pairs, reduced Gram (2), Ritz rotation (2), residual, corrections, projection (1 + 2 passes x 4), nonlocal GEMMs (usually serial). Total regions per run about 3000 to 4700; a mean region carries 20 to 300 µs of work.

| Region | Code | What limits parallelism | Proposal | E: gain at T=8 | Bit-identity risk |
|---|---|---|---|---|---|
| In-block MGS (6.8 % si8, 12.2 % si16) | `real_davidson.rs:760-779` | columns of one block are serial: each does two passes against all kept columns, scalar dot and axpy over n = 2969 | right-looking first pass with dependent second pass (the complex path got this in A48, the real path never did); region costs 20 to 60 µs per column against about 6 µs of work, so only worthwhile with a column-pair or panel-dot variant | 2 to 4 % (E) | right-looking keeps each column's operation order: bit-identical; panel-split dots change association |
| EVD (9.0 % si8, 9.7 % si16 real; 11.5 % / 13.1 % complex) | `real.rs:395`, `lib.rs:703` | single dense solve of m = 20 to 120 per step, 232 µs | none cheap: range-selected solve is B06 (1 to 4 %), overlapping is impossible (everything after depends on it) | 0 to 4 % | library change, rounding |
| Nonlocal GEMM and coupling (6.0 % si8, 6.6 % si16) | `real_gamma.rs:120-164`, `real.rs:209` | 32 projectors x 2969 x w multiply-adds is below the 1,048,576 guard for w < 11: serial; the kinetic add is serial | **run the nonlocal image concurrently with the local FFT pairs** (independent: image goes to its own buffer, then `y = (local + kinetic) + image` in the same order as now) | 3 to 5 % (E: hide 55 to 80 % of 0.045 s) | none if the sum order is kept |
| Real Hamiltonian build per SCF iteration (4.9 % / 4.4 %) and solver setup (1.4 % / 1.8 %) | `gamma_real.rs:126`, `real_gamma.rs:181`, `real_davidson.rs:111-` | rebuilt for every eigensolve: copies projectors to the real layout, allocates ceil(capacity/2) = 60 pair plans, zero-fills workspaces; the potential is the only per-iteration change | **persist `RealHamiltonian`, `PairLocalWorkspace`, `RealGemmWorkspace`, `Projection` and `RealEighWorkspace` in `KState`; refresh only the potential** | 5 to 6 % (E: 0.040 s of 0.753; si16 0.164 s of 2.64), 3 to 4 % at T=1 | none (same operations, cached state; allocation order only) |
| Energy audit (6.4 % si8, 9.0 % si16), potential 2.6 %, mixing 1.9 % | `scf.rs:1308-1337`, `session.rs`, `mixing.rs` | strictly after density, single thread at Gamma | **run the audit of iteration i concurrently with iteration i+1** (needs only output density, orbitals and its own FFT plan; the last audit and any convergence decision stay synchronous); potential and XC point work could be parallel (bit-identical) | 5 to 9 % (E: 14 of 15 audits hidden) | none if the arithmetic is unchanged; needs a snapshot of orbitals and a private FFT plan; changes only when a result is available |
| Ritz rotation (13.9 % si8) | `real_davidson.rs:493`, `real.rs:177-355` | two separate 12-panel regions (V C, then HV C), each about 300 µs at T=1: T=8 slower than T=4 (0.105 vs 0.091 s) | one region with 24 pieces (the complex path has `rotate_pair`, A45) | 2 to 4 % (E) | none: same tiles |
| Projection overlap and rotation GEMMs (8.2 % + 7.3 %) | `real.rs:177-355`, `real_davidson.rs:700-` | 0.7 to 1.4 million multiply-adds per call: around the 1M guard, 12 panels of 12 µs | fuse overlap, rotation and subtraction per 256-row panel is impossible (overlap needs a reduction across all panels): two barriers per pass remain | 0 to 2 % | n/a |
| Narrow H calls (54 % of H time at T=8, 14 % of the wall) | `fused_pair.rs:416-437` | one piece per FFT pair; 1 to 2 pieces for 60 % of calls | split one pair's inverse and forward transforms into plane groups (z lines, xy planes) as stealable pieces: arithmetic per line unchanged | 4 to 7 % (E, 2x on narrow calls); needs dev08 wake data | none if the plane partition is fixed |
| H local pairs, wide (10 pairs on 8 workers) | same | two rounds | split the pair of 2 FFTs further or accept; limit 5x | 1 to 2 % | none |
| Density, potential FFT | `band_density.rs`, `hartree.rs` | density parallel; potential serial Hartree plus XC | parallel FFT lines for Hartree (lines are independent) | 1 to 2 % | none |
| k outer width | `scf.rs:802-806`, `scf.rs:440-470` | outer = min(T, nk), inner = floor(T/outer) | see section 6 | | |

What does not limit the width on the Mac: dispatch (0.19 to 1.6 µs resident, A42), thread creation (persistent pool), sleep after 32 rounds (THR-3), `inner` budget at one k point (it is T), per-band loop count (the pair count is fixed by the block, not the thread budget).

Check against the ledger. P1 is not in the ledger: P30A text says "the Hamiltonian snapshot is still rebuilt per eigensolve" and P33A removed only solver-owned step allocations; this is a different object (H and plans). P2a (nonlocal concurrent with local) differs from A43 and A58b, which changed the nonlocal traversal order and lost on dev08: here the traversal is unchanged and only its time is hidden, so the A43/A58b mechanism (projector streaming, L2 size) is not repeated; A43's revert note says to judge on dev08 with paired alternation, which the request below does. P2c is the A45 mechanism ported to the real path. P2b (speculative audit) is new (backlog B20 C7 skips the audit, a history change; this keeps it). The narrow-block FFT split revives A34 (packed-lane FFT threads, rejected on spawn plus gather cost with scoped threads) with the stated new evidence: persistent pool, resident join 0.19 µs, and 54 % of H time at T=8 sitting in one- to two-piece calls; A28 (band-grouped 3D FFT) is not revived (no grouping across bands). A33 and A35 (scoped OS threads) are not revived. A46 (leaf caps) and A49 (finer GEMM tiles, changes rounding) are not touched. A47 is supported (spin policy gives nothing, THR-3). A50 (dynamic H pieces under outer k, rejected on scratch-slot collisions with an inner budget below the pool) is revived in P3 with new evidence: scratch slots per pool worker (inner = pool size) remove the collisions that caused the regression, and the measured k444s gain is 6 to 11 %.

## 6. k-point cases (item 4)

Measured on the Mac with `[symmetry] mode = "auto"` (5 runs, spread at most 3 %): k222s has 2 irreducible points (not 3), k222 has 3, k444s has 10; k444s without symmetry has 32.

| Case | T=1 (s) | T=2 | T=3 | T=4 | T=8 | T=10 | T=12 | T=16 | fit f | note |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| k222s sym (2 k) | 0.1864 | 0.1218 | 0.1207 | 0.1048 | 0.0956 | | | | 0.438 | outer 2, inner T/2 |
| k222 sym (3 k) | 0.2656 | 0.2120 | 0.1526 | 0.1532 | 0.1320 | | | | 0.419 | T=2: two k on one worker |
| k444s sym (10 k) | 0.7519 | 0.4191 | | 0.3083 | 0.2487 | 0.1965 | 0.1968 | 0.1983 | 0.207 | flat above 10 |
| k444s full (32 k) | 2.1693 | | | 0.6812 | 0.4776 | | | 0.3891 | 0.110 | 5.58x at T=16 |

Why T=8 and T=16 scale poorly with few k points:

1. **Outer width is the k count.** `scf.rs:802-806` sets outer = min(T, nk) and inner = max(1, T/outer). For 10 k at T=16 outer is 10 and inner is 1 (`inner = 16/10`), so six threads never get work: T=10, 12 and 16 give the same wall (0.1965, 0.1968, 0.1983 s). For 3 k at T=4 inner is 1 (T=3 and T=4 equal, 0.1526 and 0.1532 s); at T=8 inner is 2 and the inner parallelism (complex path, H bands and GEMM tiles) scales 1.3x to 1.5x (k222s sym: per-k solve 0.062 s at T=1 and 0.049 s at T=4 with inner 2; diagonalization phase 0.1252 s at T=1 to 0.0402 s at T=8).
2. **Imbalance under the halving partition** (`ordered_solves_at`, `scf.rs:440-470`, splits k and workers in half recursively): 10 k on 8 workers gives at most 2 k on a worker (efficiency 10/16 for equal costs); 3 k on 2 workers gives 2 + 1 (k222 sym T=2 1.25x). Per-k solve cost at T=1 is uniform (max/mean 1.04 for 10 k), so the loss is the partition, not unequal k.
3. **Memory and core effects.** The per-k solve is 0.061 s at T=1 but 0.092 s at T=10 and T=16 (k444s sym, sum of 11 SCF iterations); the sum over iterations of the slowest k is 0.066 s at T=1, 0.107 s at T=8, 0.119 s at T=10 and T=16, an efficiency of 0.93 / 0.75 / 0.77 (instrumented run, one sample each; the Mac heterogeneous cores are a confounder).
4. **Serial remainder.** At T>=10 the diagonalization is 0.119 s of 0.1965 s; the rest is serial or per-k setup: process and CLI 0.022 (symmetry detection and preflight: 71 of 246 main-thread active samples sit in `Config::symmetry_map`, which is evaluated at seven call sites; `reduce_mesh` 59 samples, `reciprocal_action` 39), setup 0.022, potential 0.016, density 0.018, energy 0.010, mixing 0.005. Floor for 10 k with one thread per k: 0.061 + 0.076 = 0.137 s (5.5x); measured 3.8x. For k222 sym at T=8 the serial remainder is 0.079 s of 0.132 s (60 %). The dev08 observation that symmetric k222 is slower than the full mesh at T=8 and T=16 (0.390 vs 0.299 s relayed) has the same cause: 8 independent solves became 3.

What a hybrid outer x inner decomposition needs:

- A single budget: every k solve forks pieces regardless of its own inner share and the pool decides who helps whom. That is what `WorkerPool` (one Rayon pool, nested `join`) was designed for; the only obstacle is the `threads <= 1` serial guards and the scratch count in `Hamiltonian::with_band_pool`. With inner = T (one scratch set per pool worker, so no slot is shared between workers) the H apply and GEMM paths already work: the band partition, 256-row tiles and ordered reductions are independent of T.
- Measured THR-4 (scratch edit of `scf.rs:804` to `inner = options.threads.max(1)`, five runs each, instrumented build): k444s sym wall T=4 0.3083 to 0.2933 s (-4.9 %), T=8 0.2487 to 0.2215 s (-10.9 %), T=12 0.1968 to 0.1818 s (-7.6 %), T=16 0.1983 to 0.1864 s (-6.0 %); k222s sym T=4 0.1048 to 0.1029 s (-1.8 %), T=8 0.0956 to 0.0938 s (-1.9 %). density.bin sha256 `d609b9dd44c63506...` for the unchanged build at T=1, 8, 16 and for the patched build at T=8 and 16 (one case only: k444s sym; the full gate set is not run). Cost: T scratch sets per k (10 k x 16 x one FFT scratch: tens of MB).
- Then the serial remainder is the limit (38 to 41 % at T>=10): compute `symmetry_map` once and pass it through (E: 5 % of the T>=10 wall, 71 samples directional), parallel per-k setup is already used (setup 0.040 s at T=1 to 0.022 at T=10), and the audit/energy plus potential costs shrink only through B20 C2/C3 style or overlap work.

## 7. Ranked proposals (gain / risk) and what to ask dev08

Ranking is expected gain divided by risk, for the Mac numbers above; none is accepted without dev08.

| Rank | Proposal | Expected gain (E) | Risk | Ledger check |
|---|---|---|---|---|
| 1 | P1 persistent real Hamiltonian, pair plans and solver workspaces across SCF iterations (section 5, row 4) | T=8: 5 to 6 % (si8 0.040 s, si16 0.164 s); T=1: 3 to 4 %; may be larger on dev08 where memset is 7 % of samples | very low: bit-identical, same operations | new object; not A51/P33A; no rejected row |
| 2 | P2 overlap of independent phases: (a) nonlocal image concurrent with local FFT pairs, (b) audit of iteration i concurrent with iteration i+1, (c) V and HV rotation in one region | (a) 3 to 5 %, (b) 5 to 9 %, (c) 2 to 4 %; together 10 to 20 % at T=8, 0 at T=1; Amdahl: si8 real 0.753 to about 0.62 s (T1/T8 1.59x to about 1.9x) | low to medium: bit-identical by construction; (b) needs a private FFT plan and an orbital snapshot; (a) must be judged on dev08 with paired alternation (A43, A58b) | (a) not A43/A58b (traversal unchanged); (c) = A45 for the real path; (b) new |
| 3 | P3 full-pool inner budget for k solves (per-worker scratch) plus one `symmetry_map` | k444s sym T=8 to 16: 6 to 11 % measured here (density identical, one case), plus about 5 % from the single symmetry map (E) | low: scheduling only, but rerun the whole bit-identity gate across T | revives A50 with new evidence (no slot collisions); A33/A35 not revived |
| 4 | P4 narrow-block FFT split into plane groups (stealable pieces) | 4 to 7 % (E) at T=8 si8 real | medium: wake cost on dev08 unknown; A34 premise stated | revives A34 |
| 5 | P5 right-looking in-block MGS for the real path | 2 to 4 % (E), si16 more | medium: only worthwhile with batching | A48 ported; rounding unchanged |
| 6 | parallel Hartree FFT lines and XC points, residual/precondition (1 to 2 % each) | 1 to 3 % total | low | none |

Not recommended: spin or sleep policy in Rayon (THR-3), another capacity or tolerance change (A13, A21, P31A), band-grouped FFT (A28), scoped threads (A33, A35), finer GEMM tiles (A49).

dev08 requests (through the orchestrator): (1) the wake microbenchmark (`cargo test --release -p dft-numerics wake_latency -- --ignored --nocapture` in the scratch tree; add it as an ignored test in the ticket that needs it) at T=2, 4, 8, 16 with `taskset`, to decide whether P4 is viable and whether wake costs are 7 to 85 µs there; (2) per-phase wall times for si8 real strict T=1 and T=8, which needs a `cfg`-gated phase-timer module in the ticket (the scratch module is `crates/dft-numerics/src/prof.rs` in the scratchpad tree, 40 lines; a ticket should add it behind a feature, never default); (3) three-point alternation per the brief for P1 and P2, T=1, 4, 8, 16, 10 rounds; (4) for P3, k444s and k222 with symmetry at T=8 and T=16, density hash across T per backend.

## 8. Limits of this study

- Mac only, heterogeneous cores, shared host; `sample` under-reports wall time (about 0.65 of requested) and gives shares; GEMM tile samples carry no caller frames (stolen pieces) so GEMM per phase comes from the wall-clock guards, not from `sample`.
- The instrumented build counts 291 reduced-matrix calls for si8 real strict (relayed dev08 b022b61: 258 Davidson steps; trajectories differ by rounding and by whether the certificate pass is counted); the phase values are wall time per call site, not an exclusive profile (EVD is also timed inside nothing else).
- The Amdahl phase-sum assumes perfect scaling of the parallel phases as the optimistic bound; estimates marked E are arithmetic on the tables above, not measurements of an implementation.
- THR-4 checked identity on one case (k444s sym) at two thread counts; not the generality set.
- No QE source was read; no remote numbers were produced.

## 9. Reproduction (scratchpad `/private/tmp/claude-501/-Users-mino-Develop-rust-dft/6774a210-70b0-490c-9ba0-b3db1cceafa1/scratchpad/`)

- Build: `CARGO_TARGET_DIR=$S/target cargo build --release --offline -p dft-cli` (exit 0, 45.6 s). Instrumented tree `$S/inst` (copy of `crates`, plus `dft-numerics/src/prof.rs` and edits in `real_davidson.rs`, `lib.rs` (solvers), `real_gamma.rs`, `gamma_real.rs`, `scf.rs`), target `$S/target_inst`. Patched Rayon: `$S/rayon-core-patched` (`ROUNDS_UNTIL_SLEEPY` 32 to 4000), tree `$S/inst2`, target `$S/target_inst2`. Inner-budget experiment: tree `$S/inst3`, target `$S/target_inst3`.
- Wall and Amdahl: `amdahl.py <case> <real> <std> <reps> <T list> [sym]`; fits `fit.py`; raw results `$S/w/*.json`, summaries `$S/amdahl_out.txt`, `$S/amdahl_k.txt`.
- Phase timing: `iprof.py` (`isum2.py`, `wid.py` tables); sampling: `susp.py`, `samplerun.py`, `bucket2.py` (raw `sample` files `$S/w/sm/*.txt`).
- Wake microbenchmark: `cargo test --release --offline -p dft-numerics wake_latency -- --ignored --nocapture` in `$S/inst` (6.6 s, exit 0). Median region times (12 pieces of 20 µs, T=8): extra 14.5 µs with gap 0, 35 to 59 µs with gaps 20 to 100 µs, 53 µs at 300 µs, 21 µs at 1 ms; (12 x 60 µs): 19 to 87 µs; T=4: 1.8 µs at gap 0, 7 µs at gaps 20 to 300 µs; T=2: 8 to 10 µs.
