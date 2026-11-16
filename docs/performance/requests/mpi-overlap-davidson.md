# dev08 request: overlap-matrix Davidson on the M2 path, stage 0 then stage 1 (ticket MPI-OVERLAP-DAVIDSON, 2026-11-17)

MPI against MPI only (BAANDI `--features mpi,mkl` at np ranks; comparison-rules.md: no threads number, no one-core number and no QE number of another parallel model in any table of this request). Design: docs/design/mpi-overlap-davidson.md. Ledger rows MPI-OVERLAP-DAVIDSON-*. Protocol: docs/performance/fast-loop-protocol.md (stage 0 = one solve, stage 1 = one SCF per cell).

## Why
From np 4 to 8 the panel Gram-Schmidt of the correction block does not shrink (a chain of about one tree sum per appended column) and is about a quarter of the np=8 warm si32 solve (MPI-REMAINING-MAP, O). `BAANDI_M2_LEAN_ALGEBRA=on` replaces it on standard (Ritz-change) solves by the serial lean solver's algebra: unit-norm corrections appended without Gram-Schmidt, the overlap S kept beside H_r with its inverse Cholesky factor (pivot rule eps/1e-12, the serial fallbacks), generalized reduced problem. Two tree sums per iteration instead of about 4 + k. Strict solves unchanged. Mac stage 0 (O, si32 warm solve): tree sums -61 to -79 percent, H columns and Davidson iterations unchanged, solve -12 to -25 percent.

## Build (branch HEAD of this ticket)
- `MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-cli --features mpi,mkl --example davidson_bench` and `cargo build --release -p dft-cli --features mpi,mkl` (production runs); `cargo build --release -p dft-cli --features mkl` (base densities). Record sha256 and build exit codes.
- **FFI diff for this run: none.** No `unsafe`, no new MPI call: the new path uses the existing deterministic tree sum only (fewer calls). Replicated small products go through the existing fixed-panel GEMM workspace.
- `MKL_NUM_THREADS=1`, one thread per rank, `mpirun -np P --map-by core --bind-to core` on socket 0 (cores 0..P-1), `BAANDI_MPI_TREE=rd` (default), load below 0.30 before each process, uptime before and after, INVALID rule of P0.3-QE-MPI-INVALID (exit code, empty output, wall below 0.1 s).

## Part A: stage 0 (one solve), si32 and si64, np 2/4/8, three alternating pairs per np
Inputs and base densities as in requests/mpi-scaling-gap.md Part A2 (`python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case $C --emit-input ~/ovl/$C.toml`, then `target/release/dft-cli run ~/ovl/$C.toml --output-dir ~/ovl/$C-cli` once per cell with the MKL build).
Defaults since 0d87d6e include mirror sticks, the Gram start and the kept start (the Mac stack arm S of the design's screens); arm L adds the switch. Arm O (the old start, `BAANDI_M2_START=ortho`) is a side check that the lean path also pays when the start is orthonormalized by Gram-Schmidt.
| arm | environment | identity |
|---|---|---|
| A | defaults (mirror, gram, reuse on) | |
| L | `BAANDI_M2_LEAN_ALGEBRA=on` | rounding class |
| O | `BAANDI_M2_START=ortho` | rounding class |
| OL | `BAANDI_M2_START=ortho BAANDI_M2_LEAN_ALGEBRA=on` | rounding class |
```
for C in si32 si64; do for P in 2 4 8; do for PAIR in 1 2 3; do for ARM in A L O OL; do   # rotate the arm order per pair
  case $ARM in A) E=;; L) E=BAANDI_M2_LEAN_ALGEBRA=on;;
    O) E=BAANDI_M2_START=ortho;; OL) E="BAANDI_M2_START=ortho BAANDI_M2_LEAN_ALGEBRA=on";; esac
  env $E BAANDI_DAVIDSON_STAGES=off mpirun -np $P --map-by core --bind-to core target/release/examples/davidson_bench \
    --input ~/ovl/$C.toml --density ~/ovl/$C-cli/density.bin --repeats 5 \
    | grep DAVBENCH_JSON > ~/ovl/A-$C-np$P-$ARM-p$PAIR.json
done; done; done; done
```
Then once per (cell, np, arm) with `BAANDI_DAVIDSON_STAGES=on` (stage table only, never timed): `seconds` of `dist.orthonormalize`, `dist.reduced_gram`, `dist.reduced_eigen`, `dist.rotate`, `dist.residual`, `dist.apply_total`, `dist.solve_start`, `coll.tree_sum`, `coll.alltoall`. On the lean arms `dist.reduced_gram` contains both border reductions (S border with the residual norms, H border) and therefore the wait that `dist.residual` held before; `dist.reduced_eigen` contains the factor extension and the border of A; `dist.orthonormalize` only the fallback stages.

Counters to report per arm: `coll.tree_sum.calls`, `coll.tree_sum.bytes`, `coll.alltoall.calls`, `h_columns`, `davidson_iterations`, `dist.lean_solves`, `dist.lean_border_reductions`, `dist.lean_appends`, `dist.lean_reprojections`, `dist.lean_fallbacks`, `dist.lean_rebuilds`, `dist.lean_wasted_columns`, `dist.lean_factored_starts`, `gemm_flop`.

Pass/fail per (cell, np), arm L against A and OL against O:
1. Correctness: occupied eigenvalues within max(10 tol, 1e-10) Ha, same converged roots, eigenvalue FNV identical between the processes of one arm (fixed-np reproducibility; the binary prints rank 0's JSON with the maximum over ranks of every repeat; rank agreement is covered by the tests below).
2. Counters: `h_columns` and `davidson_iterations` equal (the Mac shows equality in every cell; a difference is a finding, not a failure, if 1 holds), `coll.tree_sum.calls` lower, `dist.lean_rebuilds` = 0, `dist.lean_wasted_columns` = 0.
3. Time: median solve time of L not worse than A beyond the pair spread (and OL against O).
Mac expectation (O, section 8 of the design, measured on the earlier defaults with the ortho start, i.e. arms O -> OL): si32 warm tree sums 580 -> 120 (np=2), 683 -> 261 (np=4), 664 -> 257 (np=8); solve B/A cold 0.86 / 0.94 / 0.90, warm 0.75 / 0.88 / 0.82 at np 2/4/8. With the Gram start (arm A -> L) the Mac SCF screen gave tree sums 4494 -> 511 per rank and -10.3 percent wall at si32 np=8.

## Part B: stage 1 (SCF), si32 and si64, production default, seeds 601 602 603, np 2/4/8
Only if Part A passes 1 and 2 in every cell. Inputs `python3 -I scripts/fast_screen_dav.py --upf <UPF> --case $C --seed $SEED --emit-input ~/ovl/$C-$SEED.toml`; arms A and L (O and OL optional), interleaved per seed:
```
for C in si32 si64; do for SEED in 601 602 603; do for P in 2 4 8; do for ARM in A L; do   # rotate the order per seed
  (env as in Part A) mpirun -np $P --map-by core --bind-to core target/release/dft-cli run ~/ovl/$C-$SEED.toml \
    --output-dir ~/ovl/B-$C-$SEED-np$P-$ARM
done; done; done; done
```
Report per row: wall; SCF / Davidson / H; verification attempts; total energy; density SHA; `timings` (diagonalization, ladder coarse); `parallel.stages`, `parallel.collectives` (calls, bytes, seconds), `parallel.reductions` (sites `lean_border` and `reduced_gram` on the lean arms).
Correctness: |dE| of L against A within the standard class 5e-9 Ha, converged, verification passed, every rank's energy equal to rank 0, two repeats at np=4 with identical density SHA. One strict cell: si64 `BAANDI_REFERENCE_STRICT=on` np=8 seed 601, arms A and L: energy bits equal (strict solves never take the lean path; this checks that nothing leaks).

Also run (MKL rounding): `cargo test --release --features mkl -p dft-solvers --test distributed_lean_algebra`, `cargo test --release --features mkl -p dft-engine --test m2_lean_algebra`.

## Decision (AGENTS.md adoption rule)
L becomes the default when its median wall improves at np 2, 4 and 8 on both cells with no correctness failure and no regression beyond noise. Expected (E): the orthonormalize stage and most of its waits leave the solve (dev08 effective cost about 10 us per reduction including waits, MPI-SCALING-GAP: about 400 fewer tree sums per warm si32 solve, about 4 ms, plus the local projection compute of the stage, two passes of V^T w against the basis per column); together 5 to 15 percent of the np=8 solve, more at si64 (larger k per iteration), partly offset by the factor and A border (replicated, O(m^2 k) per iteration). Risk: the pivot test fails on most appends of late warm solves at capacity 3 (np >= 4; 15 of 22 on the Mac si32 warm solve), where stage 1 re-runs a short in-block Gram-Schmidt; if that dominates on dev08, the follow-up is a cheaper stage 1 (backlog).
