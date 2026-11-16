# MPI remaining terms after the 2026-11-17 adoptions (ticket MPI-REMAINING, 2026-11-17)

MPI against MPI only (comparison-rules.md). No threads number, no one-core number and no QE number of another parallel model appears in any table or ratio below; the one-core QE stage decomposition is used only for a one-core statement (section 3.2), and Part A of the dev08 request asks for QE's MPI timers so that the per-stage comparison can be made MPI against MPI. Process: fast-loop-protocol.md (stage 0 = one solve at a fixed Hamiltonian, stage 1 = one SCF). Labels as in mpi-scaling-gap.md: **D** dev08 relay (coordinator), **O** Mac observation (pure build, simulated ranks = threads of one process on a shared 12-CPU Mac at load 3 to 6; a simulated tree sum is two thread barriers and an alltoall a memcpy, so collective seconds are waits, not MPI latency; np=8 shares the Mac with other sessions), **C** exact count or read from the code, **E** estimate.

Starting point (D, coordinator, 2026-11-17): defaults M2_LAZY on, M2_PANEL_GS fused, VERIFY_REUSE half, M2_SELECTED_ROOTS on, MPI_TREE rd, workspace on; BAANDI MPI / QE MPI wall at equal ranks si32 np 2/4/8 = 2.60 / 2.92 / 3.20, si64 np 4/8 = 3.11 / 3.64 (target 0.90). Unadopted dev08 stage-0 candidates: M2_STICKS=mirror (alltoall calls halve, bytes -27.6 percent, solve -3 to -5 percent), M2_START=gram (tree sums -25 to -35 percent, solve -2 to -7 percent), the three together with PANEL_GS fused (MGF) solve -10 to -15 percent in every cell (SCF-level confirmation pending), M2_REUSE (Mac: -21 to -27 percent H columns; dev08 pending).

## 0. Plain summary
- What limits scaling now (np 4 -> 8, new default stack): inside the solve, the stages that do not shrink are the panel Gram-Schmidt orthonormalization (a chain of small reductions per appended column), the replicated reduced eigenproblem and the waits at the first synchronization after the rotation; together about 55 to 60 percent of the np=8 warm si32 solve on the Mac (section 1). Outside the solve, the ladder coarse stage (about a quarter of the np=8 SCF wall, flat from np 4 to 8), the fixed-point verification (about a tenth) and the replicated full-grid stages (potential, energy audit, mixing, stop rule: 6 to 9 percent of the wall and growing with np because every rank streams the same full grid) (section 2).
- Against QE (one core against one core, the only stored decomposition): the Davidson algebra is the stage family furthest from QE (5.8x at si64), while the H application is 1.6x. The same family is what does not scale under MPI. QE builds the reduced overlap with one matrix product per iteration and solves a generalized reduced problem; BAANDI's distributed solver orthonormalizes every appended column by Gram-Schmidt with its own reductions. That is the largest single lever left and the next ticket (section 6), not a switch of this one.
- Implemented here (default off, stage-0/1 screened, dev08 request written): `BAANDI_MIX_INCREMENTAL=on` (the Pulay mixer keeps its difference vectors and Gram entries across steps; bit-identical; mixing phase -48 percent on the Mac, Gram dot products -55 percent, difference vectors formed -77 percent) and `BAANDI_STOP_METRIC=kept` (the stop-rule Coulomb estimate from the Fourier coefficients the Hartree potential and the energy audit already computed: no FFT pair per SCF map, 22 full-grid FFT pairs per si32 run removed; the estimate is rounding class, the SCF state was bit-identical in every screened cell). Both remove replicated work from every rank, so their absolute saving is the same at every np and their share grows with np. Expected (E): together 1 to 3 percent of the np=8 wall on dev08.
- Checked here: the pending solver candidates together (mirror sticks, Gram start, fused panel Gram-Schmidt, kept start, lazy orbitals) and with the two new switches keep the standard class (|dE| <= 1.6e-14 Ha against threads-only on cubic and skewed cells at 2 and 3 ranks; the two grid switches change the stacked energy by 0 Ha), ranks agree, a fixed rank count is bit-reproducible (new test `tests/m2_remaining.rs`). Mac screens of the solver stack (mirror + Gram start + kept start on the adopted defaults): stage 0 si32 solve -5 to -10 percent at np 2/4/8, stage 1 SCF wall -7 to -10 percent (si32 np 4/8, si64 np 4), -12 percent at si32 np=8 with the two grid switches added; H columns -26 to -27 percent, alltoall calls -50 percent, tree sums -25 to -29 percent (section 4.3). The stack is the biggest immediate dev08 item.
- Not done, with the reason: (b) the replicated reduced eigenproblem (root solve, B5) has no changed premise beyond its larger share; (c) overlapping the transposes with nonblocking alltoall: dev08 shows the transposes scale (521 us per 34-pair call at np=8, 3.1x from np 2), the alltoall part is about 2 to 3 percent of the np=8 solve, so the best case is below the noise of a new FFI path: REJECT before code.

## 1. Stage 0: per-stage scaling of the new default stack (O, exact counters)
`davidson_bench --sim-ranks P`, si32 / si64 f1e-3 seed 601 at the converged `dft-cli run` density, cold = atomic start at 5e-3 Ha, warm = kept start at the ratchet tolerance (4.3e-11 Ha), default capacity rule (6 per band at np=2, 3 at np >= 4: np=2 warm is a different algorithm, compare np 4 and 8), 3 repeats, `BAANDI_DAVIDSON_STAGES=on` (stage clocks), medians, ms per rank, main e09067d. With M2_LAZY the per-solve interface is gone: `solve_return` 0.0 to 0.3 ms and `solve_unpack` 0 at every np (MPI-SCALING-GAP measured 5.8 / 33 ms at si32 / si64 np=8 before).

| stage (ms per rank) | si32 warm np2 | np4 | np8 | np4->8 | si64 warm np2 | np4 | np8 | np4->8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| solve | 286.5 | 209.2 | 180.0 | 1.16x | 2519 | 1630 | 1284 | 1.27x |
| apply_total | 120.4 | 93.5 | 67.6 | 1.38x | 742 | 551 | 400 | 1.38x |
| - fft_transposes (1D FFTs + alltoall) | 81.2 | 62.2 | 41.2 | 1.51x | 424 | 286 | 212 | 1.35x |
| - mirror exchange | 2.6 | 4.2 | 6.1 | 0.69x | 15.7 | 18.9 | 23.9 | 0.79x |
| - nonlocal | 33.3 | 24.5 | 18.4 | 1.33x | 280 | 231 | 151 | 1.53x |
| orthonormalize (panel Gram-Schmidt, fused) | 68.7 | 49.8 | 48.0 | 1.04x | 891 | 536 | 392 | 1.37x |
| rotate | 40.8 | 22.4 | 11.9 | 1.88x | 304 | 209 | 111 | 1.88x |
| reduced_gram | 22.6 | 15.6 | 14.0 | 1.11x | 245 | 168 | 115 | 1.46x |
| reduced_eigen (replicated) | 23.3 | 18.5 | 22.4 | 0.83x | 241 | 201 | 219 | 0.92x |
| residual (first sync after rotate) | 1.7 | 4.2 | 8.8 | 0.48x | 30.5 | 24.5 | 75.5 | 0.32x |
| solve_start | 5.8 | 4.6 | 4.8 | 0.96x | 58.5 | 29.4 | 26.6 | 1.11x |
| tree_sum waits (simulated) | 7.5 | 22.4 | 49.6 | | 141 | 106 | 256 | |
| H columns / Davidson iterations | 244 / 10 | 272 / 23 | 272 / 23 | | 557 / 11 | 603 / 24 | 603 / 24 | |
| tree_sum calls | 580 | 683 | 664 | | 1516 | 1489 | 1354 | |
| alltoall calls / MB per rank | 42 / 110.8 | 94 / 63.8 | 94 / 35.8 | | 46 / 501.7 | 98 / 293.8 | 98 / 150.9 | |
| sum of m^3 (reduced solves) | 9.6e7 | 4.9e7 | 4.9e7 | | 1.39e9 | 7.08e8 | 7.08e8 | |

Cold (same work at every np: si32 3 Davidson / 190 H columns, si64 3 / 374): si32 solve 187.9 / 101.7 / 90.1 ms (apply 109.5 / 58.4 / 40.2, orthonormalize 25.0 / 14.5 / 12.9, reduced_eigen 3.2 / 3.3 / 4.2); si64 974 / 731 / 465 ms (np=4 si64 cold had a load spike in solve_start, 102 ms).

Reading: per rank, the H application and the rotation scale (1.4 to 1.9x per doubling); the orthonormalization (si32 flat from np 4 to 8), the replicated reduced eigenproblem (grows, 0.83 to 0.92x) and the residual synchronization (waits for the slowest rank after the rotation: grows) do not. At np=8 warm si32 the non-shrinking group (orthonormalize 48.0 + reduced_eigen 22.4 + reduced_gram 14.0 + residual 8.8 + mirror 6.1 + start 4.8 ms) is 104 ms of 180 (58 percent). On the Mac the orthonormalize flatness is partly the simulated barrier (two thread barriers per reduction); on dev08 the PANEL_GS fused result (-4 to -5 percent of the solve for 55 to 60 percent fewer tree sums, D) gives the effective cost of about 10 us per reduction including waits, and the remaining compute of the stage (two passes of V^T w against the whole basis per appended column) is local and scales.

## 2. Stage 1: SCF phases of the new default stack (O, `m2_entry_cost` production = standard + verification + ladder + capacity rule, si32 f1e-3 seed 601, per rank, max over ranks)
| phase (s per rank) | np2 | np4 | np8 | np4->8 |
|---|---:|---:|---:|---|
| wall | 5.292 | 3.255 | 3.047 | 1.07x |
| SCF / Davidson / H | 11 / 70 / 2530 | 11 / 67 / 2568 | 11 / 67 / 2568 | |
| fine diagonalization | 2.939 | 1.643 | 1.376 | 1.19x |
| fine potential (Hartree FFT pair replicated, XC distributed + allgatherv) | 0.061 | 0.068 | 0.096 | grows |
| fine density (distributed + plane gather) | 0.207 | 0.128 | 0.099 | 1.29x |
| fine energy audit (Hartree FFT replicated, sums over the full grid) | 0.100 | 0.078 | 0.095 | grows |
| fine mixing (Pulay, l2, full grid, replicated) | 0.049 | 0.058 | 0.082 | grows |
| ladder coarse session (54x54x25, 4189 PW, 9 steps) | 1.062 | 0.716 | 0.756 | flat |
| - coarse diagonalization | 0.787 | 0.487 | 0.520 | flat |
| verification (certificate clock) | 0.700 | 0.374 | 0.293 | 1.28x |
| outside the recorded fine phases | 1.936 | 1.279 | 1.300 | flat |

VERIFY_REUSE=half took the run from 13 to 11 SCF steps against the MPI-SCALING-GAP table (same case). The replicated fine grid stages (potential + energy + mixing) are 0.21 / 0.20 / 0.27 s (4.0 / 6.2 / 8.9 percent of the wall) and grow with np on the Mac; MPI-SCALING-BEYOND-8 (D) showed the same stages doubling per rank from np 8 to 16 on dev08 (one socket's memory bandwidth shared by every rank streaming the same full grid). Not in any phase clock: the stop-rule Coulomb estimate (one forward and one inverse full-grid FFT per map, about 22 per si32 run including the redo, verification and coarse maps).

## 3. The biggest remaining terms
### 3.1 np-independent (do not shrink from np 4 to 8; O, si32 np=8)
| term | where | per rank at np=8 | share of the wall | class |
|---|---|---:|---:|---|
| ladder coarse session | outside the fine loop | 0.756 s | 25 % | per-iteration costs of a 4189-PW problem (replicated reduced eigenproblem of the same order as the fine one, reductions); the serial coarse stage costs 1.19 s, so running it replicated is no fix |
| panel Gram-Schmidt orthonormalization (all solves) | solve | 0.372 s | 12 % | latency chain + basis projections |
| verification | outside | 0.293 s | 10 % | an SCF map; its cost is the attempts (MPI-FS-VERIFY-COST), not MPI |
| replicated grid stages (fine potential, audit, mixing; plus the untimed stop rule) | per SCF map | 0.27 s + stop rule | 9 % + | replicated full-grid work, bandwidth-bound |
| reduced eigenproblem (all solves) | solve | 0.190 s | 6 % | replicated dense EVD |

### 3.2 Against QE (one core against one core only; qe-stage-decomposition-2026-11-14.md)
Per stage at si64 strict-equivalent f1e-3: BAANDI / QE h_apply 1.6x, rotate 12x (regterg:last), orthonormalize 7.5x (regterg:overlap), dense eigen 2.9x (rdiaghg), reduced 1.8x (regterg:update); the Davidson algebra together 5.8x, potential 1.7x, density 1.2x, mixing 0.75 to 0.91x. The stages furthest from QE at one core (orthonormalization, dense eigen, rotation) are the stages that scale worst under MPI (section 1): QE's regterg never orthonormalizes the expanded basis by Gram-Schmidt; it forms the reduced overlap with one product per iteration and solves the generalized reduced problem (rdiaghg). The MPI-against-MPI version of this statement needs QE's MPI timers (request Part A).

## 4. Candidates of this ticket (screens; O on the Mac, counters exact)
### 4.1 (a) `BAANDI_MIX_INCREMENTAL=on` (dft-solvers `PulayMixer`)
The plain Pulay update forms all m difference vectors (2 m grid vectors allocated and written) and all m^2 Gram dot products at every step, although only the newest difference is new. The switch keeps the differences (and, with a metric, their metric images) and the Gram entries across calls, keyed by the serial numbers of the pairs and seeds they come from, and forms only the newest difference and Gram row; reset, seeding and restore renumber and drop the cache. Same operands, formulas and summation orders: bit-identical (unit test over window filling and sliding, checkpoint/restore with an unchanged pair count, reset, seeds that age out, l2 and a metric, capacities 2 and 6; engine test on the threads-only and the M2 path: energy and density bits and every count equal). Applies to every parallel model (the mixer is replicated on every rank and on the threads path). Counters `dft_solvers::mixing::mixer_counters()` (mix calls, differences formed, dot products, metric applications), printed by m2_entry_cost as `grid_counters`.

| cell (O) | differences formed | dot products | mixing phase per rank | wall A -> B (3 pairs) | identity |
|---|---|---|---|---|---|
| si32 np=4 | 78 -> 18 (-77 %) | 476 -> 216 (-55 %) | 0.054 -> 0.028 s (-48 %) | 3.024 -> 2.991 s (-1.1 %, spread 0.12 s) | bit-identical |
| si32 np=8 | 78 -> 18 | 476 -> 216 | 0.083 -> 0.043 s (-48 %) | 2.921 -> 2.834 s (-3.0 %, spread 0.17 s) | bit-identical |

### 4.2 (a) `BAANDI_STOP_METRIC=kept` (dft-engine stop rule)
The standard stop rule transforms the residual rho_out - rho_in forward and back on the full grid at every map, although the Hartree potential has just transformed rho_in and the energy audit rho_out (forward, same normalization). The switch keeps both coefficient sets (the input copy inside `hartree::evaluate_keeping_density`, bit-identical potential; the audit's array is moved, not copied) together with the two densities they belong to, and evaluates E_H = N sum_G w(G) |rho_out(G) - rho_in(G)|^2 (Parseval) when the densities handed to the estimate are bit-equal to the kept ones and neither needed clamping; otherwise it transforms as before. Rounding class of the estimate (unit test: 1e-12 relative to the transform form on a skewed cell); the estimate only gates the stop and the first-step ratchet, so the SCF state changes only if a decision sits on a rounding edge (bit-identical in every Mac cell and engine fixture). Counters `dft_engine::stop_metric_counters()` (kept, transformed).

| cell (O) | stop-rule FFT pairs per run | kept estimates | wall A -> B (3 pairs) | identity |
|---|---|---|---|---|
| si32 np=4 | 22 -> 0 | 22 | 3.327 -> 3.292 s (-1.1 %, 5 pairs, spread 0.05 s; audit coefficients moved) | bit-identical |
| si32 np=8 | 22 -> 0 | 22 | 3.239 -> 3.186 s (-1.6 %, 5 pairs, spread 0.13 s; first build with copied audit coefficients: 3.223 -> 3.223 s at load 5) | bit-identical |
| si64 np=4 (with MIX_INCREMENTAL) | 17 -> 0 | 17 | 12.70 -> 11.83 s (load 5 to 9, one A run at 17.9 s: spread 6 s, not judged); mixing phase 0.105 -> 0.054 s | bit-identical |

### 4.3 (d) the combined stack: mirror sticks + Gram start + fused panel Gram-Schmidt + kept start + lazy orbitals (+ the two grid switches)
Correctness (C, test `tests/m2_remaining.rs::combined_mpi_stack_keeps_the_standard_class`, production standard with verification, cubic and skewed cells, P = 2 and 3): |dE| against threads-only 2.7e-15 to 1.6e-14 Ha, SCF / Davidson / H / verification counts equal to the threads-only run (the kept start brings the M2 H count down to the threads count, 173 against 269 without it on the cubic cell at P=2), the grid switches change the stacked energy by 0 Ha, ranks bit-equal, fixed rank count bit-reproducible. The kept start had been tested with lazy orbitals (MPI-LAZY-ORBITALS) and the mirror/Gram switches without them (MPI-SCALING-GAP); this is the first test of all of them together.
Screens (O; arm S = `M2_STICKS=mirror M2_START=gram M2_REUSE=on` against the defaults, which already carry fused and lazy):

| stage 0, si32, 3 pairs (davidson_bench sha256 f8decef703fc40bc) | cold solve B/A | warm solve B/A | alltoall calls cold / warm | tree sums cold / warm | warm H columns |
|---|---|---|---|---|---|
| np=2 | 0.951 (spread 1 %) | 0.897 (8.5 %) | 14 -> 7 / 42 -> 21 | 252 -> 169 / 580 -> 496 | 244 -> 176 |
| np=4 | 0.956 (3.9 %) | 0.913 (1 %) | 14 -> 7 / 94 -> 47 | 224 -> 151 / 683 -> 609 | 272 -> 204 |
| np=8 | 0.940 (1 %) | 0.932 (6.3 %) | 14 -> 7 / 94 -> 47 | 214 -> 145 / 664 -> 594 | 272 -> 204 |

| stage 1 (production SCF, 3 pairs, m2_entry_cost) | wall A -> B | H columns | alltoall calls / bytes | tree sums | dE |
|---|---|---|---|---|---|
| si32 np=4, S | 3.318 -> 2.992 s (-9.8 %, spread 0.03 s) | 2568 -> 1888 (-26.5 %) | -50 % / -37.6 % | -26.0 % | 0 (energy bits equal) |
| si32 np=8, S | 3.192 -> 2.960 s (-7.3 %, spread 0.20 s) | 2568 -> 1888 | -50 % / -41.0 % | -25.3 % | 1.4e-14 Ha |
| si64 np=4, S | 12.80 -> 11.88 s (-7.2 %, spread 0.35 s) | 4353 -> 3165 (-27.3 %) | -50 % / -40.7 % | -29.0 % | 5.7e-14 Ha |
| si32 np=8, S + MIX_INCREMENTAL + STOP_METRIC=kept | 3.183 -> 2.802 s (-12.0 %, spread 0.05 s) | 2568 -> 1888 | -50 % | | 1.4e-14 Ha |

Every arm reproducible at fixed np (density hash equal between repeats), SCF and Davidson counts unchanged; the kept start was taken 20 times per si32 run (1360 H columns saved per rank).

### 4.4 (b) replicated reduced eigenproblem and small GEMMs: not done
MPI-REDUCED-EIGEN-DEV08 (D): root solve +0.15 / -0.51 / -1.18 percent, B5 +0.29 / +0.70 / +0.75 percent of the si64 strict wall at np 2/4/8. Premise check: the eigenproblem's share of the warm solve grew (si32 np=8 warm 12 percent, si64 17 percent, O) because the gather and unpack are gone, but the mechanism is unchanged: a root solve removes contention, not the critical path (every rank waits for the broadcast), and B5 only acts on single-column iterations (11 of 127 solves in the si32 run). Capacity 2 was rejected in MPI-SCALING-GAP-CAP2. A changed premise needs either a distributed dense eigensolver (no library adopted) or fewer reduced solves; the second is the overlap-matrix Davidson of section 6.

### 4.5 (c) nonblocking alltoall overlapping the transposes with the next pair batch: REJECT before code
Premise: the transposes are on the critical path. dev08 (D, mpi_transpose_bench, MPI-LAZY-ORBITALS-FINDING): a 34-pair si32 transpose takes 1593 / 889 / 521 us at np 2 / 4 / 8 (it scales 3.1x); a warm si32 solve has 47 transposes of about 6 pairs each, about 5 ms at np=8 against a solve of about 200 ms (2 to 3 percent). The `fft_transposes` stage of section 1 is mostly the 1D FFT compute, which already scales and which an overlap would compete with for the same memory bandwidth. Best case below 3 percent of the solve, for a new FFI surface (MPI_Ialltoall, request handles, unsafe) and a restructured batch pipeline: not worth a switch. Mirror sticks (pending dev08) already halve the alltoall calls.

## 5. What each fix could give on dev08 (E)
| fix | mechanism | expected at np=8 | scaling statement |
|---|---|---|---|
| MIX_INCREMENTAL | O(m^2 N) -> O(m N) replicated grid work per step | 0.5 to 1.5 % of the wall | saving per rank independent of np; share grows with np and with socket bandwidth sharing |
| STOP_METRIC=kept | one full-grid FFT pair and four grid passes per map removed | 0.5 to 1.5 % | same |
| stack S (mirror + gram + kept start, with fused and lazy adopted) | half the alltoalls, fewer reductions, 21 to 27 % fewer H columns | 10 to 20 % of the solve (D stage-0 MGF -10 to -15 %, kept start not yet on dev08) | acts on the distributed part and the latency chain; the H-column saving scales like the apply |
| overlap-matrix (QE-style) distributed Davidson (next ticket) | no per-column Gram-Schmidt; S and H borders in one reduction per iteration | removes most of the orthonormalize stage (12 % of the np=8 wall, O) and part of the reduced_gram | removes a latency chain whose length does not shrink with np |
| distributed full-grid FFT (backlog) | Hartree potential, audit and stop rule on planes | up to the remaining replicated grid share (5 to 8 % at np=8) | turns replicated work into 1/np work plus two transposes |

## 6. Next ticket (sized, not started here): overlap-matrix Davidson on the M2 path
The serial lean solver (B-LEAN-DAVIDSON) already does it on one rank: unit-norm corrections appended without projection, the overlap S = V^T V kept beside H_r, an inverse Cholesky factor bordered per block, the reduced problem T H_r T^T, pivot rule and two fallback stages. On the distributed path the borders V^T W, W^T W, V^T H W and W^T H W are local GEMMs plus one tree sum per iteration (they can share the reduction that `reduced_gram` already makes), which replaces the panel Gram-Schmidt chain (about one reduction per appended column with fused, plus two basis projections per column). Work: port the lean bordering, pivot rule and fallbacks into `distributed_davidson/workspace.rs`, behind `BAANDI_M2_DAVIDSON=lean`; stage-0 counters: tree sums per solve, orthonormalize seconds, H columns (the lean fallback rate decides whether H columns grow).

## 7. Determinism and accuracy
- MIX_INCREMENTAL: bit-identical by construction and test; no tolerance involved.
- STOP_METRIC=kept: the decision to use kept coefficients is replicated (bit comparison of replicated densities), every rank takes the same branch; fixed rank count bit-reproducible; rounding class of the estimate only. Strict runs have no standard stop rule and never use it.
- Default path unchanged: `dft-cli run` si32 f1e-3 seed 601 T=1 (threads path, pure build) density.bin SHA-256 9491950d470503bf... with the main e09067d binary, with this branch, and with this branch under `BAANDI_MIX_INCREMENTAL=on BAANDI_STOP_METRIC=kept` (C).
- Accuracy classes unchanged (strict 1e-10 Ha, standard 5e-9 Ha, verification unchanged).
- No `unsafe`, no FFI change, no allocation added in a hot loop (the incremental mixer reuses dropped difference buffers; the kept coefficients reuse two grid buffers per problem).

## 8. dev08 request
docs/performance/requests/mpi-remaining.md: QE MPI stage timers (Part A), stage 0 of the solver stack (Part B), stage 1 of I, K, S and SIK (Part C).
