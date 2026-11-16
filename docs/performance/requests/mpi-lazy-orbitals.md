# dev08 request: MPI lazy orbitals, stage 0 then SCF level (ticket MPI-LAZY-ORBITALS, 2026-11-17)

MPI against MPI only (BAANDI `--features mpi,mkl` at np ranks; no threads number and no QE number in any table of this request). Design and consumer inventory: docs/design/mpi-lazy-orbitals.md. Ledger rows MPI-LAZY-ORBITALS-*. Protocol: docs/performance/fast-loop-protocol.md (stage 0 = one solve, stage 1/2 = SCF level).

## Why
The dev08 transpose microbenchmark (2026-11-17, relayed) showed that the N x B allgatherv at every M2 solve return grows with the rank count (every rank receives the whole block): si32 (11877 x 68) 557 / 1434 / 2050 us per call at np 2 / 4 / 8, si64 (23871 x 132) 4551 / 7823 / 9994 us. `BAANDI_M2_LAZY=on` removes it from the SCF loop: the solve returns local rows only, the next solve starts from them, the audit takes the kinetic energy from them (one reduction of B values), and the full complex block is gathered once at the end (and at a ladder hand-off, checkpoint or explicit audit). Expected per solve: the gather (above), the complex unpack of N x B per rank, the ORC pack of the warm start and its row extraction; per SCF run about 30 solves.

## Build
- Branch HEAD of this ticket. `MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-cli --features mpi,mkl --example davidson_bench`, `cargo build --release -p dft-cli --features mpi,mkl` (SCF runs), `cargo build --release -p dft-cli --features mkl` (base density). Record sha256 and build exit codes. No FFI change in this ticket.
- `MKL_NUM_THREADS=1`, one thread per rank, `mpirun -np P --map-by core --bind-to core` on socket 0, load below 0.30 before each process, uptime before and after; INVALID rule of P0.3-QE-MPI-INVALID (exit code, empty output, wall below 0.1 s).

## Part A: stage 0, one solve (si32 and si64, np 2 / 4 / 8, A/B on the same binary)
```
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si32 --emit-input ~/lazy/si32.toml
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si64 --emit-input ~/lazy/si64.toml
target/release/dft-cli run ~/lazy/si32.toml --output-dir ~/lazy/si32-cli      # mkl build, once
target/release/dft-cli run ~/lazy/si64.toml --output-dir ~/lazy/si64-cli
for C in si32 si64; do for P in 2 4 8; do for PAIR in 1 2 3; do
  for ARM in $( [ $((PAIR % 2)) = 1 ] && echo "A B" || echo "B A" ); do
    env $( [ $ARM = B ] && echo BAANDI_M2_LAZY=on ) mpirun -np $P --map-by core --bind-to core \
      target/release/examples/davidson_bench --input ~/lazy/$C.toml \
      --density ~/lazy/$C-cli/density.bin --repeats 5 \
      | grep DAVBENCH_JSON > ~/lazy/A-$C-np$P-$ARM-p$PAIR.json
  done
done; done; done
```
Report per (case, np, arm): `solve_ms_median` and `map_ms_median` of cold and warm; `coll.allgatherv.calls`, `coll.allgatherv.bytes`, `m2.lazy_returns`, `m2.lazy_starts`; seconds `coll.allgatherv`, `dist.solve_return`, `m2.solve_prepare`, `m2.solve_unpack`, `dist.apply_total`; the eigenvalue FNV and the converged roots.

Pass/fail (per case and np, against A):
- Identity: the eigenvalue FNV of B equals A's (the start is the same bits and the gathered block is exact data movement: the solve is bit-identical); converged roots equal; the FNV identical between the processes of one arm (fixed-np reproducibility). Any difference is a FAIL (report it, do not adopt).
- Counters: `m2.lazy_returns` = 1 per timed map in B (0 in A); `m2.lazy_starts` = 1 per warm map in B; `coll.allgatherv.bytes` per map down by at least N x B x 8 bytes per rank (si32 6.46 MB, si64 25.2 MB; Mac np=4: si32 11.0 -> 4.57 MB per rank per map).
- Time: B's median `solve_ms` not worse than A beyond the pair spread at any np (a worse median is a stage-0 REJECT for that np). Expected (E, from the relayed per-call times): si32 -0.6 / -1.4 / -2.1 ms per solve at np 2 / 4 / 8 plus the unpack, si64 -4.6 / -7.8 / -10 ms.

## Part B: SCF level (stage 2; only if Part A passes)
Production default (standard 5e-9 Ha with verification, ladder, capacity rule, workspace, rd tree; the CLI `dft-cli run`), si32 and si64, f1e-3 seeds 601, 602, 603, np 2, 4, 8, arms A (defaults) and B (`BAANDI_M2_LAZY=on`) interleaved per seed, three pairs per cell. Also one cell of each with `BAANDI_M2_REUSE=on` on both arms (the kept start reads the same local rows). Report per row: wall; SCF / Davidson / H; verification attempts; total energy (bits); density SHA; `parallel.m2` (`lazy_returns`, `lazy_starts`, `lazy_materializations`, `lazy_start_mismatches`, stages `lazy_materialize`, `lazy_kinetic`, `solve_unpack`, `solve_prepare`); `parallel.collectives` (calls, bytes, seconds).

Pass/fail:
- Correctness: density SHA of B equals A's in every row (bit-identical SCF trajectory); SCF / Davidson / H / verification counts equal; |E_B - E_A| <= 1e-12 Ha (the kinetic energy is a distributed sum: rounding class only; the standard class 5e-9 Ha is the hard gate); converged and verified; every rank's energy equal to rank 0's; two repeats of B at np=4 with the same density SHA and energy bits. `lazy_start_mismatches` = 0. `lazy_materializations` = 1 per run without the ladder, 2 with it (coarse hand-off and result).
- One strict cell (si64 `BAANDI_REFERENCE_STRICT=on`, np=8, seed 601, arms A and B): identical bits (strict solves never go lazy).
- Unit tests on MKL: `cargo test --release --features mkl -p dft-solvers --test distributed_lazy_rows`, `cargo test --release --features mkl -p dft-engine --test m2_lazy_orbitals`.

## Decision (AGENTS.md adoption rule)
B becomes the default when its median wall improves at np = 2, 4 and 8 on both cells with no correctness failure and no regression beyond noise. The gain is expected to grow with np (the removed gather grows as P x N x B per node) and with N x B (si64 more than si32).
