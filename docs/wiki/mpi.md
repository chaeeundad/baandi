# MPI in one page

Approved by the user on 2026-11-09 (AGENTS.md scope decision U9). Contract: docs/decisions/0041-mpi-parallel-layer.md; plan and gates: docs/roadmap/mpi-m1-plan.md. Crate `dft-mpi`: FFI for the Open MPI ABI behind feature `mpi` (no MPI on the Mac; type check only), single-rank stand-in without it, `Collective` trait with an in-process simulator, deterministic `tree_sum_f64` (fixed binomial association, bit-identical on every rank for a fixed rank count; faster than allgather plus ordered sum from P = 4, one rule for every P).
Measured (dev08, shared memory, one socket, us per call, P = 2 / 4 / 8): Allgatherv 11877 x 19 columns 82 / 254 / 400; Allreduce 408 x 68 44 / 80 / 129; tree sum 408 x 68 94 / 154 / 228; Alltoall of the correction-block transpose 39 / 37 / 26. Staying inside one socket matters (1.15 to 1.9x on reductions). Multi-node scaling cannot be measured on dev08.
Design: row-block basis with local algebra and tree reductions, H on column pairs with two all-to-all transposes per Davidson step, reduced eigenproblem replicated; Amdahl bound of a band-parallel-only design is about 1.8x at 8 ranks, below the 2.43x of threads, which is why the algebra is distributed too. Gate: si32 strict at least 2.2x at 4 ranks and 3.5x at 8 ranks over one rank.

## Status 2026-11-09 (update)
M1a to M1e are implemented and pass on simulated ranks (Mac, correctness only): `dft_solvers::distributed_davidson` (`RealDavidson::solve_distributed`), `dft_engine::parallel` (thread-local communicator; the engine takes the distributed branch only when the registered communicator has two or more ranks and the preconditioner is the default), and `dft-cli --features mpi`. Non-root ranks write to `<output>/rank<r>/`. Open: dev08 build and smoke run, then the pre-registered scaling gate (ledger row MPI-M1CDE). Non-Davidson stages are replicated per rank in M1.

First dev08 gate (ledger MPI-M1-GATE1): 1.19x at 4 ranks and 1.31x at 8 ranks on si32 strict against the thread path's 2.43x. The extra Davidson iterations (91 to 121) came from Cholesky panel orthonormalisation dropping near-dependent corrections; fixed with a distributed MGS panel (matches serial iteration counts on simulated ranks). Re-measure; the replicated non-Davidson stages are the next Amdahl limit.

Gate re-run (ledger MPI-M1-GATE2): best MPI 1.35x (np=4 x T=2) against threads 2.49x on 8 cores; extra Davidson iterations traced to one cluster-cut empty band at SCF step 10; Mac reproduction suggests sensitivity, not a bug. Results.json now carries per-collective time (`parallel.collectives`).

Ledger MPI-M1-GATE3: perturbed seeds 601-610 are the standard MPI scaling case (unperturbed si32 is anomalously easy); the extra Davidson iterations were sensitivity; communication is at most 13 percent of wall, so the remaining gap is non-scaling compute.

Ledger MPI-M1-GATE4: the distributed Davidson scales only 1.28x at 4 ranks; the replicated stages (24 percent of wall) are the smaller Amdahl term. Reduced matrix is now incremental and per-stage timers are in results.json `parallel.stages`; next measurement attributes the non-BLAS 21.8 s.

Ledger MPI-M1-GATE5: M1 granularity limit found (column pairs divided over ranks with the remainder on the first ranks; late iterations apply few columns). Further scaling needs a distributed FFT (M2). Decision pending with the user: stop at M1 or invest in M2.

Ledger MPI-SCALING-BEYOND-8: M2d reaches parity with threads at np=8 on one socket; beyond 8 ranks the replicated stages double per rank (memory-bandwidth sharing) and np=16 is 60 percent slower, on one socket as on two; the socket boundary and communication are not the cause. MPI work is frozen at M2d until a multi-node resource exists.
