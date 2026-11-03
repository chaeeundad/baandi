# ADR 0019: Optional Intel MKL backend for FFT and GEMM, selected at compile time

- Date: 2026-11-03
- Status: accepted for implementation (user suggested compile-time selection on 2026-11-03); default behavior unchanged
- Context: profile after O4c: FFT 54 percent, GEMM about 9 percent, eigensolver small. Micro-benchmark on the dev08 node (2 x Xeon Gold 6230, one core, `tools/lib_bench`, `crates/dft-numerics/examples/lib_bench.rs`, three runs each, spread 1 to 3 percent): MKL 3D complex FFT is faster than rustfft by 3.9 (24^3), 3.0 (25^3), 2.1 (27^3), 2.8 (36^3), 4.1 (48^3) times; MKL zgemm (750 x m adjoint times 750 x 8) is faster than faer by 2.5 (m = 8), 3.5 (16), 3.9 (32), 4.2 (48) times; the Hermitian eigendecomposition is close (faer 10 percent faster at 16, MKL 1.12 and 1.28 times faster at 32 and 48). Estimate for the whole run, single thread: about 1.5 to 1.7 times faster. The reference code uses the same MKL.

## Decision
1. Backend selection is a COMPILE-TIME choice by Cargo feature of `dft-numerics`: default pure Rust (rustfft, faer); `--features mkl` selects MKL for the 3D FFT and for GEMM. The eigendecomposition stays on faer (no gain). No runtime switching, no dynamic loading, no cost in the default build.
2. FFI lives in a new crate `dft-mkl` (ownership.toml: `unsafe_allowed`), an optional dependency of `dft-numerics`. Unsafe blocks document their preconditions and are tested; every other crate stays `unsafe_code = "deny"`. The MKL headers and libraries are not copied or vendored: the adapter declares the few C functions it needs itself (DFTI descriptor functions, `cblas_zgemm`) and the build script links `MKLROOT/lib/intel64` (sequential threading layer; the engine owns its threads).
3. Licensing: MKL is proprietary (Intel Simplified Software License); it is never committed or redistributed by this project, the user supplies it (`MKLROOT`), and the default build does not need it. The MIT license of this repository is unaffected. MKL is x86-64 only: the Apple Silicon development machine keeps the pure Rust default (an Accelerate backend is a possible later addition behind another feature).
4. Numerics: backends differ in rounding, so results are bit-identical only within one backend. Within a backend the existing guarantees stay (identical for every thread count). Across backends the gates are the energy difference (1e-12 Ha for the benchmark cases) and the unchanged QE agreement. Normalization conventions stay as documented (forward divides by N, inverse unscaled): MKL scaling is set through the descriptor.
5. The pruned transforms (zero-lane skipping) keep their API. For MKL they may be implemented by batched 1D transforms over the needed lanes or by the full 3D transform followed by zeroing the outputs outside the pattern (the documented contract); the choice is made by measurement (full MKL 3D at 25^3 is already about 1.8 times faster than the pruned rustfft transform on that node).
6. Testing: default CI (`scripts/check.sh`) runs without MKL. MKL builds are type-checked on any machine with `cargo check --features mkl` (declarations only) and built and tested on a machine that has MKL (dev08) with `MKLROOT` set.

## Alternatives
- Runtime-selected backends or dynamic loading: rejected (complexity, ABI risk, no benefit).
- FFTW: GPL, rejected by the project license policy.
- Linking MKL through a third-party `-src` crate that downloads or builds it: rejected (license and offline-build policy).

## Consequences
New crate `dft-mkl`; a Cargo feature; documentation of the build; a benchmark protocol for the two backends (`tools/lib_bench`, `scripts/bench_scf.py`). Open: Accelerate for macOS, MKL for the eigensolver, call-structure changes that profit from batched transforms (ticket P7PD).
