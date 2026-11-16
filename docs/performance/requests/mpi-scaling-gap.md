# dev08 request: MPI scaling gap, stage 0 then SCF level (ticket MPI-SCALING-GAP, 2026-11-17)

MPI against MPI only (BAANDI `--features mkl,mpi` at np ranks; no threads number and no QE number in any table of this request). Analysis: docs/performance/mpi-scaling-gap.md. Ledger rows MPI-SCALING-GAP-*. Protocol: docs/performance/fast-loop-protocol.md (stage 0 = one solve, stage 1/2 = SCF level).

## Why
The si32 solve speeds up 2.2x from np=2 to np=8 (cold 235 / 141 / 109 ms, warm 428 / 264 / 200 ms). The model of the analysis (section 3) attributes the non-scaling part mainly to a constant of about 0.3 ms per H column in the distributed apply (FFT transposes and mirror exchanges) and about 1.7 ms per Davidson iteration (replicated reduced eigenproblem and reductions). Part A measures which stage carries it (the Mac cannot: simulated collectives are memcpy and barriers). Parts B and C time the two new switches.

## Build
- Branch HEAD of this ticket. `MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-cli --features mpi,mkl --example davidson_bench`, `cargo build --release -p dft-mpi --features mpi --example mpi_transpose_bench`, `cargo build --release -p dft-cli --features mpi,mkl` (production runs) and `cargo build --release -p dft-cli --features mkl` (base density). Record sha256 and build exit codes. No FFI change in this ticket (the new example uses `alltoall_f64`, `allgatherv_f64` and `barrier` only).
- `MKL_NUM_THREADS=1`, one thread per rank, `mpirun -np P --map-by core --bind-to core` on socket 0 (cores 0..P-1), load below 0.30 before each process, uptime before and after; INVALID rule of P0.3-QE-MPI-INVALID (exit code, empty output, wall below 0.1 s).

## Part A: where the time goes (no switch)
A1. Transpose and gather microbenchmark, np = 2, 4, 8, three runs each:
```
for P in 2 4 8; do for R in 1 2 3; do
  mpirun -np $P --map-by core --bind-to core target/release/examples/mpi_transpose_bench
done; done
```
Report every line (us per call worst rank, node GB/s) for si32 and si64, alltoall at 1, 6 and B/2 pairs and the N x B allgatherv.

A2. Real-MPI stage-0 per-stage table, si32 and si64 (f1e-3 seed 601, base density from `dft-cli run`), np = 2, 4, 8, three processes each, default capacity rule and once more with the capacity fixed at 3 (identical work at every np):
```
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si32 --emit-input ~/scal/si32.toml
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si64 --emit-input ~/scal/si64.toml
target/release/dft-cli run ~/scal/si32.toml --output-dir ~/scal/si32-cli      # mkl build, once
target/release/dft-cli run ~/scal/si64.toml --output-dir ~/scal/si64-cli
# capacity-3 copies: add "davidson_capacity_per_band = 3" under [scf] -> si32c3.toml, si64c3.toml
for C in si32 si32c3 si64 si64c3; do B=${C%c3}; for P in 2 4 8; do for R in 1 2 3; do
  mpirun -np $P --map-by core --bind-to core target/release/examples/davidson_bench \
    --input ~/scal/$C.toml --density ~/scal/$B-cli/density.bin --repeats 5 \
    | grep DAVBENCH_JSON > ~/scal/A2-$C-np$P-r$R.json
done; done; done
```
Report per (case, np): `solve_ms_median` and `map_ms_median` of cold and warm, and from `seconds`: `dist.apply_total`, `m2.op.apply_fft_transposes`, `m2.op.apply_mirror_exchange`, `m2.op.apply_nonlocal`, `m2.op.apply_pack_unpack_kinetic`, `dist.orthonormalize`, `dist.rotate`, `dist.reduced_gram`, `dist.reduced_eigen`, `dist.residual`, `dist.solve_start`, `dist.solve_return`, `m2.solve_prepare`, `m2.solve_unpack`, `coll.alltoall`, `coll.tree_sum`, `coll.allgatherv`; from `counters`: `coll.*.calls`, `coll.*.bytes`, `dist.reduced_sum_m3`, `h_columns`, `davidson_iterations`. These are rank-0 values per repeat (stage seconds of rank 0; solve_ms is the maximum over ranks).

Pass/fail of Part A (measurement, no adoption): valid when every process exits 0, the eigenvalue FNV is identical between the three processes of one (case, np) and the counters are deterministic. The decision it feeds: if `coll.alltoall` plus `m2.op.apply_mirror_exchange` per warm solve does not shrink from np=4 to np=8 (A1 us per call flat or growing), the transposes are the per-column constant and Part B arm M is the lever; if `dist.reduced_eigen` plus `dist.orthonormalize` carry most of the non-shrinking time, term 2 is.

## Part B: stage 0 A/B of the switches (si32 and si64, capacity rule default, np 2/4/8, three alternating pairs per np)
| arm | environment | identity | Mac stage 0 (O, si32 np=4, simulated) |
|---|---|---|---|
| A | defaults | | |
| M | `BAANDI_M2_STICKS=mirror` | rounding class | alltoall calls 14 -> 7 / 94 -> 47 (cold / warm), bytes -28.5 / -27.6 %, solve 0.977 / 0.981, max d eps 9.6e-15 Ha |
| G | `BAANDI_M2_START=gram` | rounding class | tree sums 581 -> 379 / 1476 -> 1274, solve 0.944 / 0.971, max d eps 3.8e-14 Ha, every start taken |
| MG | both | rounding class | |
| MGF | both + `BAANDI_M2_PANEL_GS=fused` | rounding class | (PANEL_GS fused alone: dev08 -4 to -5 %) |
```
for P in 2 4 8; do for PAIR in 1 2 3; do for ARM in A M G MG MGF; do   # rotate the arm order per pair
  case $ARM in A) E=;; M) E=BAANDI_M2_STICKS=mirror;; G) E=BAANDI_M2_START=gram;;
    MG) E="BAANDI_M2_STICKS=mirror BAANDI_M2_START=gram";;
    MGF) E="BAANDI_M2_STICKS=mirror BAANDI_M2_START=gram BAANDI_M2_PANEL_GS=fused";; esac
  env $E mpirun -np $P --map-by core --bind-to core target/release/examples/davidson_bench \
    --input ~/scal/si32.toml --density ~/scal/si32-cli/density.bin --repeats 5 \
    | grep DAVBENCH_JSON > ~/scal/B-si32-np$P-$ARM-p$PAIR.json
done; done; done
```
(same for si64). Pass/fail per arm against A at the same np: occupied eigenvalues within max(10 x tolerance, 1e-10) Ha, the same number of converged roots, eigenvalue FNV identical between the processes of one arm (fixed-np reproducibility), and the expected counters moved (M: `coll.alltoall.calls` halved and `m2.op.apply_mirror_exchange` 0; G: `coll.tree_sum.calls` down and `dist.start_gram_taken` equal to 1 per solve). An arm whose median solve time is worse than A beyond the pair spread at any np is a stage-0 REJECT for that np.

## Part C: SCF level (stage 2; only arms that pass Part B)
Production default (standard 5e-9 Ha with verification, ladder, capacity rule, workspace, selected roots, rd tree), si32 and si64, f1e-3 seeds 601, 602, 603, np 2, 4, 8, arms A and the passing arms interleaved per seed, the CLI (`dft-cli run`) with the environment switches. Report per row: wall; SCF / Davidson / H; verification attempts; total energy; density SHA; `parallel.stages`, `parallel.m2.stages`, `parallel.collectives` (calls, bytes, seconds), `parallel.reductions`.
- Correctness: |dE| within the standard class 5e-9 Ha against A in every row, converged, verification passed; every rank's energy equal to rank 0; arm density SHA of two repeats at np=4 identical. One strict cell (si64 `BAANDI_REFERENCE_STRICT=on` np=8 seed 601, arms A and MG): |dE| <= 1e-10 Ha, ranks equal, reproducible.
- Also run `cargo test --release --features mkl -p dft-solvers --test distributed_start_gram`, `cargo test --release --features mkl -p dft-solvers --lib mirror_partition`, `cargo test --release --features mkl -p dft-physics --test distributed_real_gamma mirror_sticks`, `cargo test --release --features mkl -p dft-engine --test distributed_scf m2_scaling_gap` (MKL rounding).

## Decision (AGENTS.md adoption rule)
An arm becomes the default when its median wall improves at np = 2, 4 and 8 on both cells with no correctness failure and no regression beyond noise (M only changes anything at np >= 2; G acts at every np). Expected (E, analysis section 5): M removes half of the alltoall calls and 28 to 34 percent of the alltoall bytes of every apply, worth 5 to 15 percent of the np=8 solve if the transposes carry the per-column constant; G removes about 200 reductions and one B-column Gram-Schmidt per solve (2 to 6 percent of the solve; about 30 solves per SCF run). The Part A table decides the next lever: distributed full-grid FFT for the replicated Hartree / audit / mixing stages, the engine-wide lazy orbitals (no N x B gather), or a distributed back-transform of the reduced eigenproblem.
