# MPI scaling gap: why the M2 solve scales poorly (ticket MPI-SCALING-GAP, 2026-11-17)

MPI against MPI only (comparison-rules.md). No threads number and no QE number in any table below; the QE ratios of the problem statement are the coordinator's dev08 relay and are quoted only as the motivation. Process: docs/performance/fast-loop-protocol.md (stage 0 = one solve, then the SCF level).
Labels: **D** dev08 measurement relayed by the coordinator, **O** Mac observation (pure build, simulated ranks = threads of one process on a shared 12-CPU Mac, load 3 to 6; a simulated tree sum is two thread barriers, so collective seconds are wait times, not MPI latency), **C** exact count or read from the code, **E** estimate (model on D numbers, stated assumptions).

## 0. Plain summary
- The dev08 solve (si32, one Davidson solve) speeds up only 2.2x from np=2 to np=8. Fitting time = s + w/np to the cold solve (same work at every np) gives a part that does not shrink with np of about 61 ms (15 percent of the one-rank-equivalent work, but 56 percent of the np=8 solve); for the warm solve (np 4 and 8, same work) 136 ms (68 percent of the np=8 solve).
- Splitting the dev08 times into a cost per H column and a cost per Davidson iteration (section 3) points at three non-scaling terms: (1) a constant of about 0.3 ms per H column inside the distributed H apply (FFT transposes and mirror exchanges: four alltoalls per apply whose node-level volume does not shrink with np), about 40 percent of the np=8 warm solve; (2) about 1.7 ms per Davidson iteration that does not shrink (the replicated dense reduced eigenproblem plus the chain of latency-bound reductions), about 20 percent; (3) per-solve and per-SCF-step replicated work (start orthonormalization chain, N x B pack, gather and unpack, Hartree, mixing and audit on the full grid), about 6 percent of a solve plus about 10 percent of the SCF wall.
- Fixes behind switches (default off, rounding class, fixed-np bit-reproducible): `BAANDI_M2_STICKS=mirror` (term 1: a stick and its Gamma mirror stick go to the same rank, so the two mirror alltoalls of every apply disappear: alltoall calls -50 percent, alltoall bytes -28 to -34 percent) and `BAANDI_M2_START=gram` (terms 2/3: an orthonormal start is taken after one Gram reduction instead of a panel Gram-Schmidt: tree sums -29 percent at the SCF level). Together on the Mac: si32 np=8 production wall -8.4 percent, si64 np=4 -3.8 percent (observations). The existing `BAANDI_M2_PANEL_GS=fused` (dev08 -4 to -5 percent) is the other lever on term 2. The replicated reduced eigenproblem has no cheap fix left (capacity 2 screened and rejected here; root, B5, dsyevr, threaded dsyevd rejected earlier); the replicated grid stages need a distributed full-grid FFT (next ticket).
- What each fix could give on dev08 (E, from the model; the dev08 request decides): mirror sticks 5 to 15 percent of the np=8 solve if the transposes are the per-column constant (up to a third of their volume and half of their calls go away); gram start 2 to 6 percent of the solve (it removes about 200 reductions and a B-column Gram-Schmidt per solve; at the about 10 us effective cost per reduction inferred from the PANEL_GS dev08 result that is 2 ms per solve plus the Gram-Schmidt compute).

## 1. Stage-0 per-stage scaling (O, exact counters)
`crates/dft-cli/examples/davidson_bench.rs` with simulated ranks, si32-gamma f1e-3 seed 601 at its converged density (dft-cli run), cold = atomic start at 5e-3 Ha, warm = previous orbitals at the ratchet tolerance 4.3e-11 Ha, 3 repeats per process, medians, ms per rank. New in this ticket: solver stages `solve_start` (entry to the first iteration: row extraction, kept-start check, initial orthonormalization) and `solve_return` (kept return and the N x B gather), engine stages `m2.solve_prepare` (operator update, diagonal, ORC packing of the N x B start) and `m2.solve_unpack` (ORC unpacking of the N x B result); davidson_bench reports the engine M2 stages (`m2.*`).

Capacity fixed at 3 per band (`scf.davidson_capacity_per_band = 3`) so that the work is identical at every np: cold 3 Davidson / 190 H columns, warm 23 / 272, sum of m^3 9.7e6 / 4.87e7, tree sums 581 / 1476, alltoall calls 14 / 94 per rank.

| stage (ms per rank) | cold np2 | np4 | np8 | np2->8 | warm np2 | np4 | np8 | np2->8 | class |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| solve (diagonalization phase) | 202.3 | 135.3 | 123.5 | 1.64x | 406.0 | 285.3 | 263.6 | 1.54x | |
| apply_total | 110.5 | 69.1 | 45.9 | 2.41x | 182.2 | 110.3 | 77.6 | 2.35x | distributed + transposes |
| - of which fft_transposes (1D FFTs + alltoall) | 81.0 | 50.0 | 30.7 | 2.64x | 126.2 | 74.3 | 48.0 | 2.63x | |
| - nonlocal (local GEMMs + 1 tree sum) | 23.7 | 14.1 | 10.1 | 2.34x | 47.4 | 27.6 | 20.3 | 2.34x | |
| - mirror exchange (2 alltoalls per apply) | 2.3 | 2.7 | 3.8 | 0.61x | 3.8 | 5.5 | 7.5 | 0.51x | anti-scaling |
| orthonormalize (panel Gram-Schmidt) | 28.2 | 25.5 | 29.7 | 0.95x | 99.8 | 85.1 | 96.9 | 1.03x | latency chain |
| rotate | 22.1 | 11.0 | 6.4 | 3.45x | 47.3 | 25.0 | 14.2 | 3.33x | scales |
| reduced_gram | 14.0 | 8.4 | 7.1 | 1.97x | 29.3 | 18.3 | 15.9 | 1.84x | |
| reduced_eigen (replicated dsyevd/faer) | 3.8 | 4.0 | 4.6 | 0.81x | 22.1 | 23.1 | 25.6 | 0.86x | replicated |
| residual (first sync after rotate) | 1.4 | 1.6 | 1.9 | | 3.6 | 4.2 | 7.3 | | imbalance |
| solve_start | 9.6 | 11.2 | 14.3 | 0.67x | 9.8 | 10.8 | 14.7 | 0.67x | latency chain |
| solve_return (N x B allgatherv + reorder) | 1.0 | 2.2 | 5.8 | 0.18x | 1.1 | 2.3 | 6.2 | 0.18x | anti-scaling |
| m2.solve_prepare + solve_unpack | 2.6 | 4.1 | 5.0 | 0.52x | 2.3 | 3.1 | 4.8 | 0.48x | anti-scaling |
| collectives: tree_sum seconds (waits) | 7.0 | 20.4 | 41.0 | | 19.9 | 50.1 | 96.7 | | simulated |
| collectives: alltoall seconds | 6.3 | 7.6 | 11.8 | | 10.6 | 14.0 | 19.7 | | |
| alltoall bytes per rank (MB) | 88.3 | 45.1 | 25.4 | | 125.1 | 63.8 | 35.8 | | node total cold 177 / 180 / 203, warm 250 / 255 / 286 MB (C) |
| allgatherv bytes per rank (MB) | 11.0 | 11.0 | 11.0 | | 11.0 | 11.0 | 11.0 | | node total grows with np (C) |

si64 (cap rule: 6 at np=2, 3 at np >= 4; np 4 and 8 comparable; np=8 is oversubscribed on the Mac): cold solve 956 / 679 / 690 ms, warm 2441 / 1718 / 1470 ms; warm np4 -> np8: apply 533 -> 380, orthonormalize 541 -> 457, reduced_eigen 193 -> 225 (replicated, 13 to 15 percent of the warm solve), rotate 192 -> 107, reduced_gram 155 -> 112, solve_return 13 -> 32, mirror exchange 20 -> 26 ms; sum of m^3 7.45e8 (warm, cap 3), 1.22e9 (cap 6 at np=2).

Default capacity rule (6 at np=2): the warm solve at np=2 takes 10 Davidson iterations and 244 H columns against 23 and 272 at np >= 4; the dev08 np=2 -> 4 warm comparison therefore mixes two algorithms (cap 6 has the larger replicated reduced eigenproblem, cap 3 more iterations). The threads path at the same density (T=1, kept start): 10 / 176.

Reading (C/O): the distributed compute stages scale (rotate 3.3 to 3.5x from np 2 to 8, apply 2.4x), the stages that do not shrink are the replicated reduced eigenproblem, the latency chains (initial and per-iteration panel Gram-Schmidt, about 3 reductions per appended column in the classic scheme) and the replicated N x B interface (allgatherv node volume P x N x B, pack and unpack O(N x B) per rank). On the Mac the orthonormalize flatness is mostly the simulated barrier (thread wake-up); on dev08 a recursive-doubling tree sum costs 0.3 to 2.5 us (MPI-TREE-RD-DEV08-STAGE0-ADOPT), but the PANEL_GS dev08 result (-4 to -5 percent of the wall for about 60 percent fewer tree sums) implies an effective cost of about 10 us per reduction including the waits.

### Stick partition (C, `sticks.py` replica of `distributed_fft::partition`)
| case | sticks | np | sticks per rank (greedy) | max/mean greedy | max/mean mirror | padding of the transpose blocks |
|---|---:|---:|---|---:|---:|---:|
| si32 (72x72x36, 11877 G) | 499 | 2/4/8 | 249-250 / 124-125 / 62-63 | 1.002 / 1.002 / 1.010 | 1.002 / 1.010 / 1.058 | equal to max/mean |
| si64 (72x72x72, 23871 G) | 1005 | 2/4/8 | 502-503 / 250-252 / 125-126 | 1.001 / 1.003 / 1.003 | 1.003 / 1.003 / 1.003 | equal to max/mean |

Load imbalance of sticks, planes and rows is not a cause (rows max/mean at most 1.004). The transposed volume per pair of columns is sticks x n0 complex values (si32 35928, three times the 11877 G coefficients), independent of np at the node level. The mirror rule costs 5.8 percent stick imbalance at si32 np=8 (66 against 62 sticks on the worst rank) and nothing elsewhere.

## 2. Amdahl estimate from the dev08 solves (D numbers, E fit)
Real MPI on dev08 (coordinator relay, si32, one solve): cold 235 / 141 / 109 ms, warm 428 / 264 / 200 ms at np 2 / 4 / 8.
- Cold (identical work at every np: 3 Davidson, 190 H columns, the same reduced orders): least squares of T = s + w/np over the three points gives s = 61 ms, w = 342 ms (residuals +3 / -5 / +5 ms). Serial fraction s / (s + w) = 15 percent of the one-rank-equivalent work (Amdahl ceiling 6.6x), but s is 26 / 43 / 56 percent of the np = 2 / 4 / 8 solve. From np 4 and 8 alone: s = 77 ms.
- Warm (np 4 and 8 only, same work at capacity 3): s = 136 ms, w = 512 ms: 21 percent of the work, 68 percent of the np=8 solve. np=2 runs capacity 6 (10 Davidson iterations, larger reduced orders) and is not on the same curve.
- A term that grows with np (latency log P, memory-bandwidth sharing between ranks of one socket) appears in s; a three-parameter fit with a log2 P term is ill-conditioned on three points (negative s) and is not used.
- Which stages must account for s: the replicated reduced eigenproblem is only about 2 ms of the cold s (sum of m^3 9.7e6; dev08 rate about 0.19 ns per m^3 unit from MPI-REDUCED-EIGEN-DEV08 si64 np=8, contention included) and about 9 ms of the warm s. The interface (pack, gather, unpack) is at most about 10 ms. The rest must be a cost per H column (next section).

## 3. Cost per H column and per Davidson iteration (E, dev08 numbers)
Model T = a(P) H + b(P) D + c with the H columns and Davidson iterations of section 1 (cold 190 / 3, warm 272 / 23 at np >= 4) and a per-solve fixed cost c of 10 to 12 ms (start, return, prepare, unpack: the Mac np=4/8 values):
- np=4: a = 0.63 ms per H column, b = 3.6 ms per iteration; np=8: a = 0.47 ms, b = 2.6 ms. Both improve only 1.35x from np 4 to 8.
- With a = alpha + beta/P and b = b0 + b1/P: alpha = 0.30 ms per H column, beta = 1.31 ms; b0 = 1.7 ms per iteration, b1 = 7.4 ms.
- At np=8 the warm solve (200 ms) splits into alpha H = 83 ms (41 percent), b0 D = 39 ms (19 percent), c = 12 ms (6 percent): together 134 ms, the Amdahl s of section 2. The cold s (61 to 77 ms) is alpha x 190 = 58 ms plus c.
- b0 (1.7 ms per iteration) matches the replicated reduced eigenproblem (warm mean order about 128, 1 to 1.5 ms per dsyevd at the dev08 rate) plus the per-iteration reductions.
- alpha (0.3 ms per H column) is the main serial term. Per appended column the distributed path does: its share of four alltoalls per apply (two FFT transposes, two mirror exchanges: si32 warm 94 calls, 63.8 MB per rank at np=4, the node-level volume 250 to 286 MB per warm solve at every np), its share of one nonlocal tree sum, and about three in-panel reductions (classic Gram-Schmidt). Which of these dominates is not decidable on the Mac (simulated alltoall is a memcpy, a simulated reduction a barrier): the dev08 request measures the alltoall and allgatherv time per call at np 2/4/8 (`mpi_transpose_bench`) and the per-stage seconds of the real-MPI stage-0 JSON.

## 4. Outside the solve: SCF level (O, `m2_entry_cost` production = standard + verification + ladder + capacity rule, si32 f1e-3 seed 601, per rank, max over ranks)
| phase (s per rank) | np2 | np4 | np8 | scaling |
|---|---:|---:|---:|---|
| wall | 6.51 | 4.81 | 4.91 | |
| SCF / Davidson / H | 12 / 73 / 2699 | 13 / 74 / 2884 | 13 / 74 / 2884 | capacity rule at np >= 4 |
| fine diagonalization | 3.45 | 2.22 | 2.14 | |
| fine potential (Hartree replicated, XC distributed + allgatherv) | 0.075 | 0.083 | 0.121 | grows |
| fine density (distributed band densities + plane gather) | 0.252 | 0.192 | 0.180 | |
| fine energy audit (Hartree FFT, local, double counting, kinetic from full orbitals: replicated) | 0.124 | 0.116 | 0.148 | grows |
| fine mixing (Pulay with the Coulomb metric: full-grid FFT, replicated) | 0.062 | 0.076 | 0.102 | grows |
| ladder coarse session (54x54x25, 4189 PW, 9 steps) | 1.23 | 0.96 | 1.03 | flat from np 4 |
| - coarse diagonalization | 0.90 | 0.70 | 0.76 | per-iteration costs dominate |
| verification (certificate clock) | 1.11 | 0.93 | 0.90 | |
| outside the recorded fine phases (setup, coarse stage, verification) | 2.54 | 2.13 | 2.22 | |

Reading: the replicated full-grid stages (potential, audit, mixing: about 21 to 29 ms per fine SCF step per rank, flat or growing with np) are about 10 percent of the np=4 wall on the Mac; MPI-SCALING-BEYOND-8 (D) showed the same stages doubling per rank from np 8 to 16 on dev08 (shared memory bandwidth). The Pulay mixer and the Hartree/metric FFTs are replicated on every rank over the full density grid; the XC is distributed (contiguous blocks plus one allgatherv), the density is distributed. The ladder coarse stage does not scale beyond np=4: at a third of the plane waves the per-iteration costs (replicated reduced eigenproblem of the same order as the fine stage, reductions) dominate. The verification is distributed like an SCF map (MPI-FS-VERIFY-COST): its cost is the failed attempts, not MPI.

## 5. Top three non-scaling terms and the cheapest fix for each
| term | share (E, np=8) | fix (switch, default off) | stage-0 (O, si32 np=4, 3 pairs) | SCF level (O) | identity |
|---|---|---|---|---|---|
| 1. per-H-column constant in the apply: 2 FFT transposes + 2 mirror exchanges per apply, node volume independent of np | about 40 percent of the warm solve | `BAANDI_M2_STICKS=mirror`: stick (z, o) and its mirror (-z, -o) assigned to one rank as a unit (`distributed_fft::partition`), so no ORC representative is split, `max_items` is 0 on every rank and both mirror alltoalls (and their pack/unpack) are skipped | alltoall calls 14 -> 7 (cold), 94 -> 47 (warm); alltoall bytes -28.5 / -27.6 percent; solve time 0.977 / 0.981 (spread 1 percent); max d eps 9.6e-15 Ha | si32 np=4: alltoall calls -50 percent, mirror stage 0.068 -> 0 s, wall neutral (spread 0.25 s); dE 1.4e-14 Ha | rounding class (rows summed in another order) |
| 2. per-iteration and per-solve latency chains (panel Gram-Schmidt reductions) and the replicated reduced eigenproblem | about 20 percent | `BAANDI_M2_START=gram`: the start block is taken as the basis when one Gram reduction finds a defect <= 1e-12 (the restart contraction rule), else the panel Gram-Schmidt; plus the existing `BAANDI_M2_PANEL_GS=fused` | tree sums 581 -> 379 (cold, -35 percent), 1476 -> 1274 (warm, -14 percent); solve time 0.944 / 0.971; max d eps 3.8e-14 Ha | si32 np=4: tree sums 17955 -> 12703 (-29 percent), 26 of 26 starts taken, solve_start 0.288 -> 0.040 s per rank; dE 5.7e-14 Ha | rounding class |
| 3. replicated per-solve interface (N x B pack, gather, unpack) and per-step grid stages (Hartree, audit, mixing) | about 6 percent of the solve plus about 10 percent of the SCF wall | no cheap fix: the full orbitals are read every step (audit kinetic energy, warm start, ladder prolongation, checkpoint), so dropping the gather is the engine-wide lazy-orbital change of docs/design/mpi-kept-start.md section 7; the grid stages need a distributed full-grid FFT (the `DistributedPairFft` plane API over the full box is the primitive) | measured and sized here; `mpi_transpose_bench` times the N x B allgatherv at np 2/4/8 on dev08 | | |

Combined (O, SCF level, both new switches): si32 np=8 production wall 4.623 -> 4.234 s (-8.4 percent, 2 pairs, spread 0.13 s), alltoall calls 620 -> 310, alltoall bytes -33.5 percent, tree sums -29.3 percent, energy bits equal; si64 np=4 14.004 -> 13.471 s (-3.8 percent, spread 0.47 s), alltoall bytes -33.0 percent, tree sums -30.2 percent, dE 5.7e-14 Ha. SCF / Davidson / H counts identical in every pair.

Screened and rejected here: capacity 2 per band (`scf.davidson_capacity_per_band = 2`, new premise: replicated reduced eigenproblem on MPI, A13 and SERIAL-CAPACITY-SCREEN were serial): stage 0 np=4 / np=8 warm H columns 272 -> 334 (+23 percent), Davidson 23 -> 19, sum m^3 -47 percent, reduced_eigen 23 -> 13 ms, solve 285 -> 298 ms (+4.5 percent) / 264 -> 281 ms (+6.5 percent); cold -8 / -5 percent (only 3 iterations). The saved eigenproblem time is smaller than the cost of the extra H columns at every np; on dev08 the per-column constant alpha makes extra columns more expensive still.

Not done, with the reason:
- Root solve plus broadcast of the reduced eigenproblem, B5 bordered update: measured on dev08 (MPI-REDUCED-EIGEN-DEV08), no consistent gain; premise unchanged.
- Partial eigensolver (dsyevr/MRRR): MKL-EVD measured it slower than dsyevd below order/bands about 9; with capacity 3 the ratio is at most 3: premise unchanged. A changed premise exists only for a distributed back-transform (each rank computes bands/P eigenvectors of the tridiagonal problem and back-transforms them, the tridiagonalization of about 4/3 m^3 of the about 9 m^3 stays replicated): needs the LAPACK pieces (dsytrd, dstemr, dormtr) on both builds; backlog.
- Threading the replicated parts with idle cores: the comparison runs np ranks on np cores (comparison-rules.md, identical cores and binding), so there is no idle core at any np (at np=8 on an 8-core binding certainly none, and np=2 runs on 2 cores). Not applicable.
- Reduction batching beyond PANEL_GS: fusing the residual-norm reduction with the correction-norm reduction of the next orthonormalization is bit-identical but saves one reduction per iteration (about 23 per warm solve, 0.2 ms at 10 us each): too small to carry a switch.

## 6. Determinism and accuracy
- `M2_STICKS=mirror`: the partition is a deterministic function of shape, indices and rank count; fixed-np bit-reproducible; rounding class against greedy (row order changes the local sums). Tests: dft-solvers lib `mirror_partition_keeps_gamma_pairs_on_one_rank` (exact partition, +G and -G on one rank for symmetric spheres, greedy splits some pairs at P >= 2), dft-physics `mirror_sticks_match_serial_without_mirror_exchanges` (H X against the serial operator below 1e-12 relative at P = 1 to 6, exactly 2 alltoalls and 1 tree sum per apply, bit-reproducible), dft-engine `m2_scaling_gap_switches_keep_classes` (cubic and skewed cells, P = 2, 3, strict dE <= 5.2e-15 Ha, standard <= 3.7e-11 Ha against threads-only; ranks bit-equal; reproducible).
- `M2_START=gram`: the decision uses one reduced Gram matrix (identical on every rank); taken starts skip the panel Gram-Schmidt (rounding class); fixed-np bit-reproducible. Tests: dft-solvers `tests/distributed_start_gram.rs` (P = 1 to 4, selected and full rotation: cold and warm starts taken, a start scaled by 1 + 1e-6 declined, eigenvalues within 1e-10 Ha of the plain start, ranks bit-equal, reproducible) and the engine test above.
- Accuracy gates unchanged (strict 1e-10 Ha, standard 5e-9 Ha); no tolerance changed.

## 7. dev08 request
docs/performance/requests/mpi-scaling-gap.md: real-MPI stage-0 per-stage table at np 2/4/8 (the open question of section 3), the transpose and gather microbenchmark, stage-0 and SCF-level A/B of the two switches.
