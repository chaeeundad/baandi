# dev08 request: history-preserving restart (ticket BREAKTHROUGH-1-RESTART, 2026-11-17)

MPI against MPI only (BAANDI `--features mpi,mkl` at np ranks; comparison-rules.md: no threads, one-core or QE number of another parallel model in any table of this request). Design and Mac results: docs/design/restart-keep.md. Ledger rows R1A-SAME-H-CONTINUE, RESTART-KEEP-R1B, RESTART-KEEP-CAP2-LEAD. Protocol: docs/performance/fast-loop-protocol.md (hard rule: every job below finishes in 10 minutes; one case, one np).

## Status in one line
R1-a (same-H continuation) was not implemented: the redo never fires in the default runs (0 of 15), nothing to measure. R1-b (`BAANDI_RESTART_KEEP`) is REJECT on exact counters at the current capacity rules (machine independent: H columns do not fall, rotation and GEMM FLOPs rise); **no dev08 timing is requested for it at capacity 3**. One optional stage-1 job tests the only lead (capacity 2 per band with kept directions against the current capacity 3 rule).

## Build (branch HEAD of this ticket)
- `MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-cli --features mpi,mkl` and the stage-1 harness binary as in requests/mpi-overlap-davidson.md part B. Record sha256 and build exit codes.
- **FFI diff: none.** No `unsafe`, no new MPI call (the kept directions use the existing rotation and border; zero extra reductions per restart in the normal case).
- `MKL_NUM_THREADS=1`, one thread per rank, `mpirun -np 8 --map-by core --bind-to core` on socket 0, load below 0.30 before each process, uptime before and after, the usual INVALID rule.

## Job 1 (optional, stage 1, one case, np=8, three alternating pairs, six production runs of a few seconds each; well inside 10 minutes)
si32-gamma f1e-3 seed 601, production standard with verification and ladder (the CLI default path), `BAANDI_M2_LEAN_ALGEBRA=on` on both arms (the candidate baseline):
| arm | environment / input |
|---|---|
| A | defaults (capacity rule: 3 per band at np >= 4) |
| B | `scf.davidson_capacity_per_band = 2` in the input, `BAANDI_RESTART_KEEP=b/4` |
Report per arm: wall (median of 3, spread), SCF steps, Davidson iterations, H columns (`convergence.work.h_vectors_recorded` + certificate), `parallel.reduced_solves` (sum m^3), rotation FLOPs, `parallel.collectives` tree_sum calls (the distributed keep_* counters are not in results.json; the Mac harness counts them), energy (|dE| against A, gate 5e-9 Ha), verified.
Mac expectation (O, simulated ranks, counters exact): B has +17 % H columns, +25 % rotation FLOPs, -46 % sum m^3, +7 % tree sums, equal SCF steps (11). MPI-SCALING-GAP-CAP2 found on one warm solve that +23 % H columns were not paid back by the halved reduced eigenproblem; B is expected to be slower. A B faster than A beyond the spread reopens a capacity study (stage 2: si64 np=8, a non-Si cell); otherwise the lead is closed.

## Not requested (and why)
- R1-b at capacity 3 (np 4/8) or on the threads path: counters already decide (section 4/5 of the design: H -0.6 % to +6.5 %, rotation up to +45 %; the threads path has no contraction in si8/si16/si32 SCF runs, so the switch is a no-op there).
- R1-a: no redo in the measured runs; reopen only with a workload that has redos.
