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

# dev08 profile of the real path after P30A/P32A (main 72574f9, relayed 2026-11-07)

Relayed by the dft5090/master session. eu-stack 20 ms sampling, 20 runs accumulated per configuration; si8 strict, gamma_real on, pure native frame-pointer build (sha256 `a7381682f97c05b369c3de527556ccb07450fbde3a14e66e5e522153a4b22e3f`, `-C target-cpu=native -C force-frame-pointers=yes`), `taskset -c 0` / `0-7`; sampled wall 3.14 to 3.17 s versus 3.134 s unsampled (no sampling overhead this time). Raw stacks: `/scr/mino/rustdft/d0b/stacks_si8_T{1,8}_on_strict.txt` on dev08. `*` = fewer than 30 samples, do not quote.

| Bucket (share of working samples) | T=1 on, 8894b0b | T=1 on, 72574f9 | T=8 on, 8894b0b | T=8 on, 72574f9 |
| --- | ---: | ---: | ---: | ---: |
| a local FFT | 17.8 % | 31.9 % | 20.9 % | 30.7 % |
| b nonlocal | 2.0 %* | 3.0 % | 2.1 %* | 3.7 % |
| c subspace algebra | 60.2 % | 48.3 % | 66.4 % | 57.1 % |
| d residual and correction | 2.8 %* | 1.0 %* | 4.1 %* | 0.4 %* |
| f unattributable GEMM | 0.0 % | 0.2 %* | 0.0 % | 0.3 %* |
| e other | 17.2 % | 15.5 % | 6.5 %* | 7.8 % |
| working samples | 495 | 1713 | 292 | 1014 |
| all-thread samples | 500 | 1733 | 2267 | 6481 |
| waiting share | 1.0 % | 1.2 % | 87.1 % | 84.4 % |

Top symbols T=1 (self share of all samples): 12.00 % `dft_numerics::real`, 10.10 % `perform_fft_f64`, 9.98 % `x86::f64::avx::matmul_1_4_dyn`, 8.25 % `FusedOperator::apply_local`, 6.98 % memset, 6.00 % `real_davidson::run`, 5.65 % memcpy, 3.69 % `nano_gemm::copy_millikernel`, 3.58 % `fused`, 2.48 % `matmul_imp`, 2.19 % `qr_algorithm`, 2.02 % `hamiltonian::add`, 1.90 % `cbrt`. T=8: futex 81.39 %, then `matmul_1_4_dyn` 1.99 %, `real` 1.76 %, FFT 1.48 %.
Changed by P32A: `IndexMut::index_mut` (10.6 %) left the top 15; `matmul_2_4_16` 7.2 % to 1.7 % while `matmul_1_4_dyn` rose to 10 %; memset 3.2 % to 7.0 % and memcpy 5.8 % to 5.65 % together 12.6 % (larger than the FFT kernel symbol); cbrt (XC) 1.9 %.
Parallel width at T=8 (per-thread waiting share over 20 runs): 8894b0b on strict median 95.9 %, 72574f9 on strict median 89.5 % (180 thread records), 8894b0b off strict median 80.6 %; threads with under 50 % waiting are 11 % in all three, about two of 17. The CLI reports `outer_threads 1, inner_threads 8, threads_used 1`.
Reading: algebra is still about half of the real path's working time (FFT plus nonlocal 34.9 %), zero-fill and copies take 12.6 %, and only about two worker threads are busy at T=8 in every configuration, which limits the thread scaling to about 1.24x.
