# dev08 request: MPI fast-screen candidates, stage 2 (2026-11-16)

MPI against MPI only (BAANDI `--features mkl,mpi` at np ranks; no threads and no QE number in any table of this request). Stage 2 of docs/performance/fast-loop-protocol.md; stage 3 (generality set, accuracy classes) is decided on dev08 after this. Screening evidence and profile: docs/performance/profile-2026-11-16-mpi.md; ledger rows MPI-FS-PANEL-GS, MPI-FS-TREE-RD, MPI-FS-DENSITY-FUSED.

## Candidates (all default off; one binary, environment A/B)
| arm | switch | identity | Mac stage-1 (O, si32 production np=4 simulated) |
|---|---|---|---|
| A | none (defaults) | | |
| B | `BAANDI_M2_PANEL_GS=fused` | rounding class | tree_sum calls per rank 17955 -> 7172 (-60 %), panel_norm 5824 -> 525, work counts identical, energy bits equal; np=8: 17955 -> 7056, orthonormalize stage -43 % |
| C | `BAANDI_MPI_TREE=rd` | bit-identical for np = 2, 4, 8 (other np stay binomial) | not screenable on simulated ranks (message pattern only): 2 log2 P -> log2 P rounds per tree sum |
| D | `BAANDI_M2_DENSITY=fused` | bit-identical | m2 density stage -15 %, wall within spread (NEUTRAL); removes the bands x plane-points buffer on the standard path |
| E | B + C + D | | |

## Step 0: first run of the changed FFI (cannot run on the Mac)
The only FFI change is the declaration of `MPI_Sendrecv` and its use in `World::tree_sum_f64` when `BAANDI_MPI_TREE=rd` and the rank count is a power of two (crates/dft-mpi/src/ffi.rs, `doubling_steps`); the CLI sets the pattern once at MPI start (crates/dft-cli/src/main.rs). Before any SCF arm:
1. `cargo build --release --features mpi -p dft-mpi --examples` and `mpirun -np P target/release/examples/mpi_smoke` for P = 2, 3, 4, 8 (must print OK).
2. `mpirun -np P --map-by core --bind-to core target/release/examples/mpi_tree_bench` for P = 2, 4, 8 on socket 0, and P = 3 (bit check skipped, must still print OK). It asserts that binomial and recursive doubling return the same bits on every rank for payloads 1, 8, 64, 408 and 27744 doubles, then prints microseconds per call of both patterns (worst rank). Report the full output. If the bit check fails, stop and do not run arm C or E.

## Build and protocol
- This branch HEAD, `cargo build --release --features mkl,mpi` with `MKLROOT=/opt/intel/oneapi/mkl/2026.1`, one binary, sha256 recorded, build exit code.
- `MKL_NUM_THREADS=1`, one thread per rank, `--map-by core --bind-to core` on socket 0, load gated as in BAANDI-MPI-RANK-SWEEP (load below 0.30 before each run), uptime before and after; INVALID rule of P0.3-QE-MPI-INVALID (exit code, empty output, wall below 0.1 s).
- Cells: si32 and si64, production default (standard 5e-9 Ha with verification, ladder, capacity rule), f1e-3 seeds 601, 602, 603, np = 2, 4, 8. Arms A to E interleaved per seed (A B C D E, then E D C B A for the next seed, ...). One strict reference cell for the class check of arm B: si64 `BAANDI_REFERENCE_STRICT=on` np=8 seed 601, arms A and B.
- Report per row: wall; SCF / Davidson / H counts; `parallel.reductions` (calls and payload per site); `parallel.collectives` (calls, payload bytes, seconds per operation); `parallel.stages` (reduced_gram, reduced_eigen, rotate, residual, orthonormalize, apply_total); `parallel.m2.stages` (density, op.band_densities, op.apply_*); the closure table (setup, phases, certificate/verification, ladder_coarse); verification attempts; total energy; density SHA.
- Medians by seed, spread per cell; any cell whose arm-A spread exceeds 5 % is repeated.

## Correctness conditions (per arm against arm A at the same np)
- C and D: energy bits and density SHA identical in every row (bit-identical by construction; any difference is a defect, stop).
- B and E: |dE| within the standard class 5e-9 Ha (strict reference cell: 1e-10 Ha), converged, verification passed, Davidson count within the run-to-run variation of A.
- Every arm: every rank's energy equal to rank 0; arm E density SHA of two repeats at np=4 identical (bit-reproducible at a fixed rank count).
- Also run `cargo test --features mkl -p dft-solvers --lib panel_gs_tests` and `cargo test --features mkl -p dft-engine --test distributed_scf m2_panel_gs_variants_keep_classes` (MKL rounding of the fused Gram-Schmidt).

## Decision (AGENTS.md adoption rule)
Each arm that improves the median wall at np = 2, 4 and 8 on both cells without a correctness failure and without a regression beyond noise becomes a stage-3 candidate (full generality set and accuracy classes before the default changes). Expected (E, from the counters): B removes about 60 % of the solver tree sums; at the dev08 latency per call (16 to 40 us, MPI-ENTRY-COST-DIAGNOSIS, wiki/mpi.md) that is 0.2 to 1 s at si64 np=8 plus the synchronisation waits it removes; C halves the rounds of every remaining tree sum; D is a memory-traffic item whose effect grows with ranks sharing one socket's bandwidth.
