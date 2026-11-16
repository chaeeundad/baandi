# dev08 request: MPI kept start on the stick rows, stage 2 (ticket MPI-KEPT-START, 2026-11-16)

MPI against MPI only (BAANDI `--features mkl,mpi` at np ranks; no threads and no QE number in any table of this request). Stage 2 of docs/performance/fast-loop-protocol.md. Design: docs/design/mpi-kept-start.md. Ledger rows MPI-FS-KEPT-START, MPI-KEPT-START.

## Candidate (default off; one binary, environment A/B)
| arm | switch | identity | Mac stage 1 (O, simulated ranks, production si32 np=4, same binary) |
|---|---|---|---|
| A | none (defaults) | | |
| B | `BAANDI_M2_REUSE=on` | rounding class | SCF / Davidson 13 / 74 identical, H columns 2884 -> 2272 (-21.2 %), operator applied columns 5670 -> 4242, inverse pair transforms 3687 -> 3021 (-18.1 %), tree sums 17955 -> 13716 (-23.6 %), alltoall bytes -4.3 %, energy bits equal, density hash differs (rounding), fixed-np reproducible; 21 of 24 armed starts taken; wall median 6.71 -> 4.75 s with a 2.7 s spread at load 9 (not evidence) |

## Build and protocol
- Branch HEAD of this ticket, `cargo build --release --features mkl,mpi` with `MKLROOT=/opt/intel/oneapi/mkl/2026.1`, one binary, sha256 and build exit code recorded. No FFI change in this ticket.
- `MKL_NUM_THREADS=1`, one thread per rank, `--map-by core --bind-to core` on socket 0, load below 0.30 before each run, uptime before and after; INVALID rule of P0.3-QE-MPI-INVALID (exit code, empty output, wall below 0.1 s).
- Cells: si32 and si64, production default (standard 5e-9 Ha with verification, ladder, capacity rule, workspace, selected roots), f1e-3 seeds 601, 602, 603, np = 2, 4, 8. Arms A, B interleaved per seed (A B, then B A, ...).
- Report per row: wall; SCF / Davidson / H counts; `parallel.m2` counters (`applied_columns`, `inverse_pair_transforms`, `forward_pair_transforms`, `kept_plane_pairs`, `potential_change_columns`, `reuse_armed`, `reuse_taken`, `reuse_declined`, `reuse_not_kept`, `reuse_gated`, `reuse_plane_bytes_max`); `parallel.collectives` (calls, payload bytes, seconds per operation); `parallel.reductions`; `parallel.stages` and `parallel.m2.stages`; verification attempts; total energy; density SHA; peak RSS per rank (`/usr/bin/time -v` max resident) for the memory statement.
- Medians by seed, spread per cell; any cell whose arm-A spread exceeds 5 % is repeated.

## Correctness conditions (B against A at the same np)
- |dE| within the standard class 5e-9 Ha in every row, converged, verification passed, SCF and Davidson counts within the run-to-run variation of A.
- Every rank's energy equal to rank 0; arm B density SHA of two repeats at np=4 identical (bit-reproducible at a fixed rank count).
- `reuse_taken` > 0 in every row and `reuse_gated` = 0 (otherwise report the grid and kept pairs).
- Also run `cargo test --release --features mkl -p dft-solvers --test distributed_kept_start`, `cargo test --release --features mkl -p dft-physics --test distributed_real_gamma potential_change`, `cargo test --release --features mkl -p dft-engine --test m2_kept_start` (MKL rounding of the kept start).
- One strict cell (si64 `BAANDI_REFERENCE_STRICT=on` np=8 seed 601, arms A and B): energy and density SHA bit-identical (strict solves never take the kept start; the switch must be inert there).

## Decision (AGENTS.md adoption rule)
B becomes a stage-3 candidate if the median wall improves at np = 2, 4 and 8 on both cells with no correctness failure and no regression beyond noise. Expected (E, from the counters and the serial C2 result on dev08: H columns -25..-29 % gave wall -3..-6 % at T=1): the saved columns remove the initial apply's nonlocal GEMMs with their tree sum, the inverse transform and the mirror exchange of the kept pairs, the orthonormalization of the warm start and the `extract_local` pass; a few percent of the wall, more at higher np where the removed tree sums (-24 %) weigh more. The N x B allgatherv at return is unchanged (F07 remainder, design section 7).
