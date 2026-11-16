# MPI profile and fast-screen candidates (2026-11-16)

MPI against MPI only; no threads number and no QE ratio in this document (comparison-rules.md). Process: docs/performance/fast-loop-protocol.md (stage 1 on the Mac, harness section "Harness (MPI)").
Labels: **M/D** measured on dev08 and relayed (ledger row named), **O** Mac observation (pure build, simulated ranks = threads of one process on a shared 12-CPU Mac, load 5 to 7; simulated ranks do not model MPI latency or socket bandwidth), **C** exact count or read from the code, **E** estimate.

## 1. dev08 stage shares (M/D)
The latest relayed full stage table is si64 STRICT np=8 at 8852dc3 (ledger MPI-COARSE-CERT-AB-FINAL, arm 4 = the adopted screen and atomic-start switches, before the capacity rule, workspace and selected roots), wall 39.05 s:

| stage | s | share |
|---|---:|---:|
| ladder coarse diagonalization | 8.38 | 21.4 % |
| reduced_eigen (replicated dense EVD, fine + coarse) | 7.43 | 19.0 % |
| apply (distributed H) | 6.98 | 17.9 % |
| orthonormalize | 3.74 | 9.6 % |
| rotate | 3.10 | 7.9 % |
| screen + atomic start | 1.6 | 4.1 % |

Changes since (M/D): capacity rule 3 per band at np >= 4 (MPI-CAPACITY-SWEEP-DEV08: si64 strict np=8 reduced_eigen 7.19 -> 2.57 s, orthonormalize 3.56 -> 4.63 s, wall -13.9 %); workspace -4 to -12.5 % wall, selected roots neutral to -3 % (M2-WORKSPACE-SELECTED-ROOTS-DEV08). reduced_eigen grows with ranks at fixed work (si64 strict 5.42 / 5.94 / 7.16 s at np 2 / 4 / 8, MPI-REDUCED-EIGEN-DEV08). Collectives (older stack, si64 strict, MPI-ENTRY-COST-DIAGNOSIS): 3.4 / 4.6 / 4.6 s at np 2 / 4 / 8, tree_sum 16 us per call at np=2.
**Not available:** a relayed stage table of the current production default (standard 5e-9 Ha with verification, ladder, capacity rule, workspace, selected roots) at any np; the coordinator reports that the latest dev08 MPI comparison table (8f4d921) is invalid (built without `--features mpi`) and is being re-measured. The dev08 request of this round asks for `parallel.stages`, `parallel.reductions` and `parallel.collectives` on every row, which closes this gap.

## 2. Mac stage shares of the production default (O)
`crates/dft-engine/tests/m2_entry_cost.rs`, `M2_PROFILE_MODE=production M2_PROFILE_LADDER=on` (CLI default path: standard + verification + ladder + capacity 3 at np >= 4), Si.pz-vbc 15 Ha, harness f1e-3 seed-601 displacements, single runs, per rank, seconds and share of the wall:

| stage | si32 np=4 | si32 np=8 | si64 np=4 | si64 np=8 |
|---|---:|---:|---:|---:|
| wall | 4.79 | 5.25 | 16.21 | 19.45 |
| SCF / Davidson / H columns | 13 / 74 / 2884 | 13 / 74 / 2884 | 10 / 51 / 4353 | 10 / 51 / 4353 |
| apply_total | 1.77 (37 %) | 1.33 (25 %) | 5.10 (31 %) | 3.85 (20 %) |
| orthonormalize | 0.62 (13 %) | 0.86 (16 %) | 2.09 (13 %) | 2.07 (11 %) |
| rotate | 0.36 (7 %) | 0.20 (4 %) | 1.45 (9 %) | 0.84 (4 %) |
| reduced_gram | 0.24 (5 %) | 0.22 (4 %) | 1.00 (6 %) | 0.83 (4 %) |
| reduced_eigen | 0.20 (4 %) | 0.26 (5 %) | 0.57 (4 %) | 0.69 (4 %) |
| M2 density | 0.30 (6 %) | 0.25 (5 %) | 0.99 (6 %) | 1.28 (7 %) |
| atomic start (split) | 0.12 | 0.11 | 0.71 (4 %) | 0.63 (3 %) |
| ladder coarse session | 0.91 (19 %) | 1.07 (20 %) | 3.62 (22 %) | 3.64 (19 %) |
| verification block (certificate clock) | 0.92 (19 %) | 0.95 (18 %) | 1.26 (8 %) | 1.76 (9 %) |
| tree_sum calls / s | 17955 / 0.37 (8 %) | 17955 / 1.19 (23 %) | 22170 / 1.10 (7 %) | 22170 / 2.76 (14 %) |
| alltoall calls / s | 620 / 0.22 | 620 / 0.34 | 394 / 0.77 | 394 / 1.17 |
| solver reductions: in-panel dot + norm | 16792 of 17877 (94 %) | 17062 of 17877 | 20372 of 22119 (92 %) | 21172 of 22119 |

Reading (C/O): the call count of the tree sums is a trajectory property (identical at every np); 92 to 94 % of the solver reductions are the in-panel dot and norm reductions of the panel Gram-Schmidt, each carrying at most 64 doubles. Their time grows with np (latency per call grows with log2 P; on the Mac also with oversubscription): the clearest MPI-specific lever. The ladder coarse session is about a fifth of the wall at both sizes; it runs the same distributed code at a third of the grid, so its fixed per-iteration costs (reductions, replicated EVD) weigh more (mpi-coarse-certificate-analysis.md 3.3).

## 3. Verification and coarse stage under MPI (coordinator request, C/O)
- Distribution (C): the verification (`session.rs` `verify_fixed_point`) is one more ordinary map: `step_kpoints_with_certification` with the warm start, i.e. the distributed Davidson (M2c), the distributed density and LDA XC (M2d), the replicated Hartree, the replicated Coulomb metric and density distance. The coarse ladder stage is an ordinary session on the coarse problem through the same M2 path; the atomic start is split over ranks (M2_ATOMIC=split, adopted); the reduced EVD of the coarse stage is replicated with the same order as the fine stage.
- Cost (O, `M2_PROFILE_VERIFY=off` sizing knob, ledger MPI-FS-VERIFY-COST): si32 np=4: SCF 13 -> 10, H 2884 -> 2322, wall 4.66 -> 3.16 s (-32 %), dE 1.97e-9 Ha. The same on the serial path at T=1 (12 -> 10 SCF, 7.01 -> 5.04 s). The verification cost is its failed attempts (two extra maps before the passing one), not its distribution: an algorithmic item for the standard stop rule (decision 3 of 2026-11-11), not an MPI item.
- Not MPI-specific either way: the verification share at np=4 (19 % O) is close to the serial T=1 share (28 % O here; 16 % at si64 T=1 on dev08, coordinator note).

## 4. Candidates (top 5 by expected share x plausibility) and stage-1 results
| rank | idea (switch) | premise / ledger | stage 1 (O/C, si32 production np=4) | verdict |
|---:|---|---|---|---|
| 1 | in-panel Gram-Schmidt reductions fused: input norm with the first dots, second pass only when the per-column DGKS rule fails, Pythagorean norm (`BAANDI_M2_PANEL_GS=fused`; `pythag` variant) | review F08, M2-REDUCTION-COUNTS-F08; not A52/MPI-M1-GATE1 (per-column decisions kept, no Cholesky) | tree_sum 17955 -> 7172 (-60 %), work counts identical, energy bits equal, wall -9 % (spread 0.25 s); np=8 -61 %, orthonormalize -43 %; pythag -31 % calls, wall neutral | **CANDIDATE** (fused); pythag dominated |
| 2 | recursive-doubling tree sum, bit-identical for power-of-two P (`BAANDI_MPI_TREE=rd`) | backlog MPI (a), mpi-reduced-eigen-analysis item 3 | not screenable on simulated ranks; rounds 2 log2 P -> log2 P (C); bit identity unit-tested | **CANDIDATE for the dev08 microbenchmark** (FFI first run) |
| 3 | kept start on the stick rows (REUSE-PSI C2 + F07: X and H X kept in local rows across SCF steps, H_new X = H_old X + dV X) | new: serial path has it, M2 does not | H columns MPI 2699 vs serial 2087 at equal Davidson counts (+29 %, C) | candidate, not implemented (exceeds the 2-iteration budget); expected -20 % H columns (E) |
| 4 | fused distributed band density on the standard path (`BAANDI_M2_DENSITY=fused`) | backlog item (f); the per-band buffer only serves the strict screen | bit-identical; density stage -15 %, wall inside the spread | NEUTRAL (rider arm on dev08: bandwidth item) |
| 5 | verification solve cost (coordinator) | section 3 | -32 % wall without verification, same on the serial path | finding: not an MPI lever; to the standard stop-rule work |

Not taken, with the ledger reason: (d) partial-eigenpair `dsyevr` for the reduced problem: MKL-EVD measured it slower than `dsyevd` below order/b about 9; with capacity 3 the order/b is at most 3, the premise is weaker than before (not repeated). (b) distributing or shrinking the coarse-stage EVD: the only shortening is a distributed dense EVD (ScaLAPACK/ELPA, library adoption, dev08 microbenchmark first, MPI-REDUCED-EIGEN item 6); root solve and B5 were measured and not adopted (MPI-REDUCED-EIGEN-DEV08); capacity already moved to 3 (MPI-CAPACITY-RULE). (e) certificate screen distribution beyond the merged one: the screen does not run on the standard production path (strict only), nothing left to distribute there. (c) the N x B gather at solve return alone: allgatherv is 1 to 4 % of the wall (O) and is subsumed by candidate 3.

## 5. Determinism and accuracy (per change)
- `M2_PANEL_GS=fused|pythag`: rounding class; every branch depends on reduced (replicated) values, so all ranks take the same branches; fixed P bit-reproducible; engine test `m2_panel_gs_variants_keep_classes` (strict dE <= 7.8e-16 Ha, standard <= 8.8e-11 Ha against threads-only), unit test `panel_gs_tests` (dependent column dropped, nearly dependent column takes the second pass, reductions per column 3 / 2 / 1-2).
- `MPI_TREE=rd`: bit-identical for P = 2, 4, 8, 16 (unit test `recursive_doubling_equals_the_binomial_tree`); P not a power of two keeps the binomial code; FFI runs on dev08 only (`examples/mpi_tree_bench.rs` asserts the bit identity on every rank).
- `M2_DENSITY=fused`: bit-identical (same per-point sum order; physics test), standard path only (the strict screen keeps the per-band data).

## 6. dev08 request
docs/performance/requests/mpi-fast-screen-2026-11-16.md (stage 2: si32/si64 production x np 2/4/8, seeds 601-603, arms default / fused GS / rd tree / fused density / all three, FFI first run).
