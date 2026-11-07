# dev08 D0 profile: si8 strict, real Gamma path on versus off, T=1 and T=8 (relayed, 2026-11-07)

Relayed by the dft5090/master session. eu-stack sampling every 20 ms, 5 runs accumulated per configuration; si8, strict, pure native, frame-pointer build (sha256 `53706fddf6a2b038ff023c9841b0a75efa7cda196f80f8a77e5df78f224e17a3`, `RUSTFLAGS="-C target-cpu=native -C force-frame-pointers=yes"`, main `8894b0b`), `taskset -c 0` or `0-7`. Symbols are Rust v0 mangled (dev08 elfutils cannot demangle them): paths were reconstructed from length-prefixed identifiers. Sampling overhead (wall time while sampling versus unsampled batch-7 B medians): T=1 off +19 %, T=1 on +15 %, T=8 off +14 %, T=8 on +9 %; used for shares only. Raw stacks on dev08: `/scr/mino/rustdft/d0/stacks_si8_T{1,8}_{off,on}_strict.txt`.

## Share of working-thread samples (waiting excluded); `*` = fewer than 30 samples, do not quote

| Bucket | T=1 off | T=1 on | T=8 off | T=8 on |
| --- | ---: | ---: | ---: | ---: |
| a local FFT | 32.7 % | 17.8 % | 27.2 % | 20.9 % |
| b nonlocal | 12.8 % | 2.0 %* | 14.9 % | 2.1 %* |
| c subspace algebra | 20.9 % | 60.2 % | 21.1 % | 66.4 % |
| d residual and correction | 3.5 %* | 2.8 %* | 4.9 %* | 4.1 %* |
| f unattributable complex GEMM | 19.4 % | 0.0 % | 23.7 % | 0.0 % |
| e other | 10.6 % | 17.2 % | 8.2 % | 6.5 %* |
| working samples | 602 | 495 | 389 | 292 |
| all-thread stack samples | 606 | 500 | 1601 | 2267 |
| waiting share (all threads) | 0.7 % | 1.0 % | 75.7 % | 87.1 % |

Bucket f: the complex path's `gemm` crate micro-kernels carry no unwind information (96 of 98 samples at stack depth 1), so those samples cannot be assigned to nonlocal or algebra; b + c + f is 53.1 % (T=1) and 59.7 % (T=8) of working samples for the complex path (upper bound for nonlocal plus algebra). The real path has no such samples (0 of 500 and 0 of 2267): it uses `nano_gemm::copy_millikernel`, `x86::f64::avx::matmul_2_4_16` and `linalg::matmul::matmul_imp`, i.e. a different GEMM backend.

## Top symbols (self share of all-thread samples)

T=1 off: 12.54 % `FusedOperator::apply_local`, 12.54 % `NonlocalOperator::add`, 11.39 % `perform_fft_f64`, 6.60 % `linalg::subtract_projection`, 5.61 % + 4.13 % + 1.82 % + 2.31 % complex `gemm` micro-kernels and packing, 3.63 % `fused`, 2.64 % memset, 2.64 % + 2.48 % `linalg::tiled::Gem(multiply)`, 2.48 % `tridiag_evd::qr_algorithm`.
T=1 on: 13.60 % `dft_numerics::real`, 11.80 % `real_davidson::run`, **10.60 % `ops::index::IndexMut::index_mut` (53 samples)**, 7.20 % `matmul_2_4_16`, 6.20 % `apply_local`, 5.80 % memcpy, 5.60 % `perform_fft_f64`, 4.80 % `nano_gemm::copy_millikernel`, 3.20 % memset, 2.80 % `cbrt`, 2.80 % `qr_algorithm`.
T=8: waiting dominates (off: futex 70.96 %, sched_yield 3.81 %; on: futex 83.86 %, sched_yield 2.16 %).

Threads at T=8: both configurations keep 17 threads; the CLI reports `threads_requested 8, threads_used 1, outer_threads 1, inner_threads 8` (single k point: no outer parallelism, band direction only, resident pool of 16 workers). Waiting share per thread: off 30.6 to 94.7 % (median about 80 %); on 25.5 to 98.0 % (all but three threads 72 to 98 %), i.e. effectively three or four threads work in the real path.

## Orchestrator reading

With the real path on, the subspace algebra (real_davidson::run, dft_numerics::real, element-wise indexing with bounds checks, memcpy, small-tile GEMM calls) takes about 60 % of the working time at one thread and 66 % at eight, while FFT plus nonlocal fall to about 20 %. The planned D0 gate (FFT plus nonlocal below 60 percent moves algebra work forward) is met for both paths; for the real path algebra clearly dominates. Single-thread real versus complex GEMM backends differ.
