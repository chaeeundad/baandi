# FFT and GEMM backend boundary (P10B stage 1)

`dft_numerics::FftPlan3d` is a compile-time alias: `RustFftPlan3d`
by default, `mkl_backend::MklFftPlan3d` under `dft-numerics/mkl`.
The explicit Rust type remains available for comparison. The existing Rust
implementation and its implementation-specific lane-count tests are unchanged;
those unit tests explicitly name the Rust type even in an MKL build. MKL API
contract tests are in `mkl_backend.rs`, FFI tests in `dft-mkl/tests/linked.rs`.
All downstream physics and engine callers use the selected alias. Eigensolves
remain on faer; only `linalg::gemm` selects CBLAS.

Arrays are standard C layout `(i,j,k)`, k contiguous. MKL owns a persistent
in-place double-complex descriptor and a full-grid buffer. Calls copy input to
that buffer, transform, then copy to output. Creation/configuration/commit and
buffer allocation happen once per plan. Forward scale is `1/(nx*ny*nz)` in the
descriptor; backward scale is the default 1. No Rust allocation or logging occurs
in transforms. The initial placement follows ADR 0019's in-place measurements;
out-of-place has not been measured in this adapter and is not claimed inferior.

Pruned inverse requires canonical zero input outside reciprocal support and
returns the complete real-space transform. Pruned forward computes the complete
transform and sets output outside support to canonical `(0,0)`. Empty support,
duplicate indices, layout and shape errors retain their contracts. In particular,
zeroing real-space inverse outputs would be incorrect. No batched pruning or
fused stick algorithm is included; P10C owns that trait and stage 2 integration.
The call-structure note is design context, not evidence for an unmeasured kernel.

`dft-mkl` confines every unsafe operation, checks all dimensions and buffer lengths,
and uses `num_complex::Complex64`'s C-compatible pair-of-f64 representation.
DFTI uses C `long`; LP64 CBLAS dimensions use C `int` (not ILP64). DFTI errors
carry their status through `FftError::Mkl`; error text is fetched only when formatted.
GEMM validates finite inputs and outputs at the established numerics boundary.
Its four identity/adjoint combinations pass checked column-major buffers directly
without converting them to faer matrices. Zero inner dimension returns exact zero.

A descriptor is `Send`, not `Sync`. Each worker owns a separate plan, and every
transform requires an exclusive mutable borrow. Moving a committed plan to a
worker is tested. Only `mkl_intel_lp64`, `mkl_sequential` and `mkl_core` are linked;
no OpenMP runtime or process-global thread-setting functions are used. Bitwise
repeatability is required within a fixed backend/build, not between vendors or
between CPU generations. Foreign-library internal allocations are not measured
by Rust allocation assertions and need profiling on dev08.

The handwritten declarations follow Intel's public
[DFTI interface documentation](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-0/dfticreatedescriptor.html).
The dev08 header probe checks the enum values and ABI sizes before execution.
MKL headers, binaries and source are never vendored.

## Fused operator API

P10C adds `dft_numerics::fused::FusedOperator`: `apply_local` maps ordered sphere
coefficients and C-order real potential samples to P V psi; `accumulate_density`
adds `weight * |psi|^2 / volume` with the weight and volume kept separate to
preserve rounding. Shape, coefficient order, last-write-wins duplicate scatter,
normalization and validation-before-write are documented in rustdoc. Plans own
all descriptors and scratch, use no internal threads and allocate only at
construction. Each existing band worker owns one plan, including a single
worker when the thread budget is one.

The compile-time `fused::FusedPlan` alias selects `OrderedFusedPlan` by default.
It uses compact x sticks and active-z inverse intermediates, then one xz plane
per real y to fuse inverse z / local potential / forward x. Compact requested
x columns feed forward y/z. Both historical x/y/z orders and per-axis inverse
scales are preserved, so the default is bit-identical to the Rust oracle.
`PlaneFusedPlan` is the explicit z-stick/xy-plane experiment; its reversed
inverse order changes bits and is not selected by the alias.

With the `mkl` feature, `FusedPlan` selects
`GridFusedPlan<dft_numerics::FftPlan3d>`, a reusable compatibility adapter around
the existing selected MKL 3D plan. This keeps P10B's native backend and numerical
sequence instead of silently replacing it with Rust FFTs. The adapter retains
full-grid scratch and is not advertised as an MKL fused optimization. The generic adapter is tested with the explicit Rust backend. MKL-feature
checking currently stops at unchanged upstream `linalg.rs` imports (`Accum`
and `matmul`); full native checking and linked validation remain pending.
P10B stage 2 can replace this alias with a measured vendor implementation of
the same trait using batched 1D stick descriptors with its own strides/distances.
No FFI, process-global thread setting or MKL descriptor code is added by P10C.

See [ADR 0026](../decisions/0026-fused-wavefunction-operators.md) and the
[measurement study](../performance/batched-fft-study.md) for the pure-Rust
oracle, pointer-stability, grouping and SCF evidence.
