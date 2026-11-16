# dev08 request: per-core kernels, stage 0 then stage 1 (ticket PERCORE-KERNELS, 2026-11-17)

Threads model only (BAANDI threads at T=1 and T=8; no QE and no MPI number in any table of this request). Analysis: [../percore-kernels-2026-11-17.md](../percore-kernels-2026-11-17.md); ledger rows PERCORE-*. Protocol: [../fast-loop-protocol.md](../fast-loop-protocol.md) (stage 0 = one Davidson solve at a fixed Hamiltonian; stage 1 = one SCF).

## Why
On the Mac pure build the transposed border products of the lean Davidson (overlap V^T E and reduced V^T (H E)) are 20 % of a warm si32 solve at about half the build's GEMM rate, and the narrow Ritz rotations run at 4.5 to 10 GFLOP/s. On dev08 MKL the transposed products were the farthest from the dgemm ceiling (DGEMM-SITES: 28 to 31 of 76 to 101 GFLOP/s, efficiency growing with the work per call), but the lean default's per-kernel MKL rates were never measured. Two switches act on MKL builds:
- `BAANDI_LEAN_BORDER=fused`: H applied to the corrections first, then one TN product V^T [E | H E] (twice the columns, one stream of V) gives both borders. MKL: rounding class (dgemm rounding can depend on the width). Mac pure: bit-identical, warm si32 stage 0 -3.6 % (spread 1.75 %).
- `BAANDI_RITZ_FUSED=loop`: the existing fused Ritz and residual kernel at the selected-root residuals only (H X_s not stored; the return keeps the pair rotation for REUSE-PSI). Bit-identical. Mac pure: REJECT (+1.9 %); asked only because the saved streams can pay where dgemm is fast.
`BAANDI_PACK_FIXED=on` acts only on the aarch64 pure kernel and has no effect on dev08 (both x86 builds: run it once in Part A on the pure build as a no-op check: identical counters and eigenvalue FNV).

## Build
Branch HEAD of this ticket. `MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-cli --features mkl --example davidson_bench` and `cargo build --release -p dft-cli --example davidson_bench` (pure, separate target directories), `cargo build --release -p dft-cli [--features mkl]` for the base densities and stage 1. Record sha256 of every binary and the build exit codes. `MKL_NUM_THREADS=1`, `taskset -c 0` (T=1) or `taskset -c 0-7` (T=8), load below 0.30 before each process, uptime before and after.

## Part A0: the per-kernel table of the lean default on MKL (no A/B)
```
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si32 --emit-input ~/percore/si32.toml
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si64 --emit-input ~/percore/si64.toml
target/release/dft-cli run ~/percore/si32.toml --output-dir ~/percore/si32-cli      # mkl build
target/release/dft-cli run ~/percore/si64.toml --output-dir ~/percore/si64-cli
for C in si32 si64; do for T in 1 8; do
  BAANDI_DAVIDSON_STAGES=on BAANDI_DGEMM_TIMERS=on taskset -c 0-$((T-1)) \
    target/release/examples/davidson_bench --input ~/percore/$C.toml \
    --density ~/percore/$C-cli/density.bin --repeats 3 | grep DAVBENCH_JSON > ~/percore/A0-$C-T$T.json
done; done
```
(the input files carry `threads = 1`; for T=8 edit `threads = 8` in a copy of the input.) Report per case and T: cold and warm `seconds.stage.*` (h_apply, local, nonlocal, reduced, orthonormalize, rotate, eigen), `counters.gemm.*` (calls, flop per site and route) and, if `BAANDI_DGEMM_TIMERS` fills them, the per-site dgemm seconds; derived GFLOP/s per site. This replaces the Mac pure table of the analysis by the build that is compared with QE.

## Part A: stage 0 A/B (same binary, environment A/B, three alternating pairs)
```
for C in si32 si64; do for T in 1 8; do for X in BAANDI_LEAN_BORDER=fused BAANDI_RITZ_FUSED=loop "BAANDI_LEAN_BORDER=fused,BAANDI_RITZ_FUSED=loop"; do
  taskset -c 0-$((T-1)) python3 -I scripts/fast_screen_dav.py --upf <UPF> --case $C --threads $T \
    --binary target/release/examples/davidson_bench --no-build --density ~/percore/$C-cli/density.bin \
    --work ~/percore/work-$C --b "$X" --json-out ~/percore/A-$C-T$T-$(echo $X | tr ',=' '__').json
done; done; done
```
MKL build all rows; the pure build (x86 faer path) for `LEAN_BORDER=fused` at T=1 and T=8 on si32, plus `BAANDI_PACK_FIXED=on` once (si32 T=1, must be identical in counters and eigenvalue FNV).
Report per row: verdict, cold and warm median ratio and spread, identity line (bit-identical or the max occupied eigenvalue difference), the counters `lean.fused_borders`, `lean.fused_wasted_columns`, `gemm.*` calls/flop, H columns, Davidson iterations.
Pass (per row): no REJECT; `lean.fused_wasted_columns` = 0; identity bit-identical or max |d eps_occ| below the harness gate; H columns and Davidson iterations equal to A or within one iteration (rounding class on MKL). Thread independence: the eigenvalue FNV of the B arm equal at T=1 and T=8 on each backend.

## Part B: stage 1 (one SCF; only for the Part A survivors)
`python3 -I scripts/fast_screen.py --binary target/release/dft-cli --upf <UPF> --env-b <switch>` on si32 f1e-3 seed 601 T=1 and T=8, MKL build, three pairs; si64 seed 601 T=1 one pair per arm if the si32 rows pass. Report wall medians, ratios and spread, energy difference (gate: standard 5e-9 Ha; strict reference run `BAANDI_REFERENCE_STRICT=on` once with each arm: |dE| <= 1e-10 Ha), density SHA identical across T=1 and T=8 within an arm, verification attempts, SCF / Davidson / H counts.

## Decision
Stage-0/1 results are screening, not adoption (fast-loop rule). A switch that improves the MKL stage-1 wall beyond the spread with no correctness failure goes to stage 2 (si16, si32, si64, k444s, mgo-40, one skewed or slab cell) and stage 3 before any default change (AGENTS.md adoption rule).
