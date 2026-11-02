# P0 architecture and P1 proposal

Nine independent crates reserve ownership boundaries. P0 has no inter-crate
or external crate dependencies; no placeholder physical types are exported.

| Crate | Future responsibility |
| --- | --- |
| dft-units | Atomic units and explicit I/O conversions |
| dft-geometry | Cell, positions, reciprocal lattice |
| dft-numerics | Safe numerical interfaces and private native adapters |
| dft-basis | Plane-wave indexing and cutoff |
| dft-upf | Independently implemented public UPF specification reader |
| dft-physics | Norm-conserving, unpolarized LDA operators |
| dft-solvers | Eigensolvers and SCF iteration |
| dft-engine | Calculation orchestration |
| dft-cli | Input/output and invocation |

Expected dependency direction: CLI → engine → solvers/physics → basis/geometry/
units/numerics; UPF feeds physics. Introduce actual edges only with their tickets.
Unsafe is denied workspace-wide; only the private Accelerate adapter allows it.
The P0 adapter is test-only, fixed-size, and exports no production backend API.
Other skeleton tests are compilation sentinels, not scientific validation.

## P1 sketch (proposal only, not implemented or type-checked)

```rust,ignore
// Complex64, error types and direction need explicit P1 definitions.
trait NumericalBackend {
    // Full n*n column-major Hermitian input; overwrite with column eigenvectors.
    // Eigenvalues ascending, length n. Validate lengths/finite entries first.
    fn eigh(&self, n: usize, a: &mut [Complex64], eigenvalues: &mut [f64])
        -> Result<(), NumericalError>;
    fn fft3(&self, shape: [usize; 3], data: &mut [Complex64], direction: Direction)
        -> Result<(), NumericalError>;
}
```

FFT proposal: contiguous x-fastest data, forward exp(-i G·r) unscaled,
inverse exp(+i G·r) divided by grid size. This discrete pair does not replace
the physical Ω^-1/2 wavefunction normalization. Units stay Bohr/Hartree;
B=2π(A⁻¹)ᵀ, lattice vectors are columns, Σ|c|²=1. P1 must document grid
mapping and verify round-trip, impulse and Parseval identities before adoption.
Boundary checks must reject dimension overflow before FFI integer conversion.
No logging, JSON or file I/O belongs in numerical kernels.

## FFT candidates — 미검증 / unverified offline recollection

No network research was performed or package selected.

| Candidate | License recollection (미검증) | Maintenance / Apple Silicon questions (미검증) |
| --- | --- | --- |
| FFTW | GPL; alternative commercial licensing may exist | Native build, SIMD/planner behavior, wrapper upkeep and distribution obligations require review |
| rustfft | MIT OR Apache-2.0 | Pure Rust portability attractive; current release, aarch64 SIMD and performance require verification |
| Accelerate vDSP | Apple SDK terms | OS-managed native support; complex packing, scaling and dimensional API require tests |

Complex BLAS/LAPACK and FFT remain P1 work. No performance or cross-platform
claim follows from these P0 real-valued smoke tests.
