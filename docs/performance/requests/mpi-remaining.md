# dev08 request: MPI remaining terms, stage 0 then stage 1 (ticket MPI-REMAINING, 2026-11-17)

MPI against MPI only (BAANDI `--features mkl,mpi` at np ranks, QE MPI at the same np; comparison-rules.md: no threads number, no one-core number and no QE number of another parallel model in any table of this request). Analysis: docs/performance/mpi-remaining.md. Ledger rows MPI-REMAINING-*. Protocol: docs/performance/fast-loop-protocol.md (stage 0 = one solve, stage 1 = one SCF per cell).

## Why
After the 2026-11-17 adoptions (M2_LAZY, M2_PANEL_GS=fused, VERIFY_REUSE=half, selected roots, rd tree, workspace) BAANDI MPI is 2.6 to 3.6x QE MPI at equal ranks. Three switches are still waiting for dev08 (M2_STICKS=mirror, M2_START=gram, M2_REUSE=on); this ticket adds two replicated-grid switches (MIX_INCREMENTAL, STOP_METRIC=kept) and checks the combined stack. Part A also asks for QE's own MPI stage timers, so that the per-stage comparison can be made MPI against MPI.

## Build (branch HEAD of this ticket)
- `MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-cli --features mpi,mkl --example davidson_bench` and `cargo build --release -p dft-cli --features mpi,mkl` (production runs); `cargo build --release -p dft-cli --features mkl` (base density). Record sha256 and build exit codes.
- **FFI diff for this run: none.** No `unsafe`, no new MPI call (the two new switches are replicated arithmetic only; the stack switches use alltoall, allgatherv, tree sums and broadcasts that exist since MPI-SCALING-GAP).
- `MKL_NUM_THREADS=1`, one thread per rank, `mpirun -np P --map-by core --bind-to core` on socket 0 (cores 0..P-1), load below 0.30 before each process, uptime before and after, INVALID rule of P0.3-QE-MPI-INVALID (exit code, empty output, wall below 0.1 s).

## Part A: QE MPI stage timers (no BAANDI change; MPI against MPI)
The stored QE stage decomposition (qe-stage-decomposition-2026-11-14.md) is one core only. For the MPI comparison we need QE's own timers at np 2, 4, 8 for si32 and si64 (f1e-3 seed 601, the P0.3 MPI inputs, conv_thr 1e-8, the fastest QE MPI configuration of P0.3): `c_bands`, `h_psi`, `vloc_psi`, `calbec`, `add_vuspsi`, `regterg`, `rdiaghg`, `regterg:overlap`, `regterg:update`, `regterg:last`, `fftw`, `sum_band`, `v_of_rho`, `mix_rho`, `PWSCF` (cpu and wall), three runs, all lines of the timing block. BAANDI side at the same np: `parallel.stages`, `parallel.m2.stages`, `parallel.collectives`, `timings` of results.json of the default production run (Part C arm A gives it).

## Part B: stage 0 (one solve), si32 and si64, np 2/4/8, three alternating pairs per np
Inputs and base densities as in requests/mpi-scaling-gap.md Part A2 (`fast_screen_dav.py --emit-input`, `dft-cli run` once per cell).
| arm | environment | identity |
|---|---|---|
| A | defaults | |
| S | `BAANDI_M2_STICKS=mirror BAANDI_M2_START=gram BAANDI_M2_REUSE=on` | rounding class |
```
for C in si32 si64; do for P in 2 4 8; do for PAIR in 1 2 3; do for ARM in A S; do   # alternate the order per pair
  E=; [ $ARM = S ] && E="BAANDI_M2_STICKS=mirror BAANDI_M2_START=gram BAANDI_M2_REUSE=on"
  env $E BAANDI_DAVIDSON_STAGES=off mpirun -np $P --map-by core --bind-to core target/release/examples/davidson_bench \
    --input ~/rem/$C.toml --density ~/rem/$C-cli/density.bin --repeats 5 \
    | grep DAVBENCH_JSON > ~/rem/B-$C-np$P-$ARM-p$PAIR.json
done; done; done; done
```
Then once per (cell, np, arm) with `BAANDI_DAVIDSON_STAGES=on` (stage table only, never timed): report `seconds` of `dist.orthonormalize`, `dist.reduced_eigen`, `dist.reduced_gram`, `dist.rotate`, `dist.residual`, `dist.apply_total`, `m2.op.apply_fft_transposes`, `m2.op.apply_mirror_exchange`, `m2.op.apply_nonlocal`, `dist.solve_start`, `coll.tree_sum`, `coll.alltoall`, and the counters `coll.*.calls`, `coll.*.bytes`, `h_columns`, `davidson_iterations`, `reuse.armed`, `reuse.columns_saved`.
Pass/fail per np: occupied eigenvalues within max(10 tol, 1e-10) Ha of A, same converged roots, eigenvalue FNV identical between the processes of one arm, `coll.alltoall.calls` halved, `coll.tree_sum.calls` down; median solve time of S not worse than A beyond the pair spread.

## Part C: stage 1 (SCF), si32 and si64, production default, seeds 601 602 603, np 2/4/8
Inputs: `python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case $C --seed $SEED --emit-input ~/rem/$C-$SEED.toml` (the stage-1 displacement stream, amplitude 1e-3). Arms interleaved per seed (`dft-cli run` with the environment switches; the CLI default path: standard 5e-9 Ha with verification, ladder, capacity rule):
| arm | environment | identity | Mac stage 1 (O, simulated ranks, si32) |
|---|---|---|---|
| A | defaults | | |
| I | `BAANDI_MIX_INCREMENTAL=on` | bit-identical (energy, density SHA, every counter) | mixing phase -48 % (np=4), Gram dot products -55 %, difference vectors formed -77 % |
| K | `BAANDI_STOP_METRIC=kept` | rounding class of the stop estimate only (bit-identical on the Mac cells) | stop-rule full-grid FFT pairs 22 -> 0 per run |
| S | `BAANDI_M2_STICKS=mirror BAANDI_M2_START=gram BAANDI_M2_REUSE=on` | rounding class | see analysis section 4 |
| SIK | S + I + K | rounding class | |
```
for C in si32 si64; do for SEED in 601 602 603; do for P in 2 4 8; do for ARM in A I K S SIK; do   # rotate the order per seed
  case $ARM in A) E=;; I) E=BAANDI_MIX_INCREMENTAL=on;; K) E=BAANDI_STOP_METRIC=kept;;
    S) E="BAANDI_M2_STICKS=mirror BAANDI_M2_START=gram BAANDI_M2_REUSE=on";;
    SIK) E="BAANDI_M2_STICKS=mirror BAANDI_M2_START=gram BAANDI_M2_REUSE=on BAANDI_MIX_INCREMENTAL=on BAANDI_STOP_METRIC=kept";; esac
  env $E mpirun -np $P --map-by core --bind-to core target/release/dft-cli run ~/rem/$C-$SEED.toml \
    --output-dir ~/rem/C-$C-$SEED-np$P-$ARM
done; done; done; done
```
Report per row: wall; SCF / Davidson / H; verification attempts; total energy; density SHA; `timings` (potential, density, energy, mixing, diagonalization, ladder coarse); `parallel.stages`, `parallel.m2.stages`, `parallel.collectives` (calls, bytes, seconds). Correctness: arm I energy bits and density SHA equal to A in every row (bit-identical claim); K, S, SIK |dE| within the standard class 5e-9 Ha against A, converged, verification passed, every rank's energy equal to rank 0, two repeats at np=4 with identical density SHA. One strict cell: si64 `BAANDI_REFERENCE_STRICT=on` np=8 seed 601, arms A and SIK: |dE| <= 1e-10 Ha, ranks equal (strict never uses the standard stop rule or lazy orbitals; this checks that nothing leaks).

Also run (MKL rounding): `cargo test --release --features mkl -p dft-engine --test m2_remaining`, `cargo test --release --features mkl -p dft-solvers --lib incremental_tests`, `cargo test --release --features mkl -p dft-engine --lib coulomb_estimate_from_kept`.

## Decision (AGENTS.md adoption rule)
An arm becomes the default when its median wall improves at np 2, 4 and 8 on both cells with no correctness failure and no regression beyond noise. I is bit-identical and acts on every parallel model (also threads and one core), so a neutral MPI result with a measured mixing-phase drop is enough to adopt it (overhead-only change); K changes only the stop estimate's rounding. Expected (E, analysis section 5): I and K together 1 to 3 percent of the np=8 wall (replicated grid work removed from every rank, more where 8 ranks share one socket's bandwidth); S 10 to 20 percent of the np=8 solve (stage-0 dev08 relays: MGF -10 to -15 percent; kept start -21 to -27 percent H columns on the Mac).
