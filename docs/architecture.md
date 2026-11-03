# Current architecture

The workspace implements a CPU, full-complex Bloch-sampled (Gamma and Monkhorst-Pack
k points, time-reversal reduction only), spin-unpolarized, insulating
norm-conserving LDA-PZ SCF. [Physics conventions](physics-conventions.md)
are authoritative for units, reciprocal vectors, FFT normalization, zero modes,
projector normalization and energy bookkeeping. Internally the engine uses
Bohr/Hartree; lattice vectors are matrix columns.

## Crate boundaries

The allowed workspace edges below come from [ownership.toml](../ownership.toml).
Names in the dependency columns omit the `dft-` prefix. Actual edges are read
from each crate's `Cargo.toml`; external dependencies are inventoried in
[dependencies](dependencies.md).

| Crate | Implemented responsibility | Allowed workspace dependencies | Declared workspace dependencies |
| --- | --- | --- | --- |
| dft-units | Unit newtypes and conversions | none | none |
| dft-geometry | Cell, atoms, reciprocal lattice and coordinates | units | none |
| dft-numerics | Complex dense algebra and reusable FFT adapter | none | none |
| dft-basis | G spheres, cutoff, grid choice and scatter/gather maps | units, geometry, numerics | units, geometry; numerics (dev) |
| dft-upf | Validated UPF reading and radial transforms | units | units |
| dft-physics | Density, kinetic, Hartree, XC, Ewald and Hamiltonian | units, geometry, numerics, basis, upf, solvers | all allowed edges |
| dft-solvers | Dense oracle, Davidson, preconditioning and mixing | numerics | numerics |
| dft-engine | SCF orchestration, convergence certificate and energy audit | physics, solvers, basis, geometry, upf | all allowed edges |
| dft-cli | TOML, command dispatch, provenance and output | engine, units, geometry, basis, upf, physics, solvers | engine, units, geometry, basis, upf; physics (dev) |

The physics-to-solvers production edge is for `LinearOperator`, as specified
by ownership. `dft_engine::input` re-exports `Atom`, `Cell`, `Pseudopotential`
and `parse_upf` for clients. `dft_physics::data` re-exports the array, complex,
unit and numerical adapter types used by the engine, avoiding a direct engine
backend dependency. No numerical kernel performs JSON serialization or file I/O.

## One SCF run

```text
TOML + user UPF files
  -> CLI Config / PreparedInput + parsed Pseudopotential
  -> Cell / PlaneWaveBasis / FftGridSpec
  -> ionic local potential + NonlocalOperator + Ewald energy
  -> ScfProblem
       input density -> Hartree + LDA-PZ -> Hamiltonian
       Hamiltonian (LinearOperator) -> Davidson -> orbitals -> output density
       -> EnergyAudit + convergence checks
       -> density mixing -> next iteration, or fixed-point certificate
  -> ScfOutcome / ScfResult -> CLI RunReport
  -> results.json + density.json + density.bin
```

1. `dft_cli::config::Config::parse` rejects unknown structured fields;
   `Config::prepare` constructs `PreparedInput` with geometry, wavefunction and
   density bases and an FFT grid. Input lattice-vector rows become internal
   columns; the input cutoff is converted from Ry to Ha. K points come from the
   legacy Gamma list, a Monkhorst-Pack `kgrid` with `kshift`, or an explicit
   weighted list ([CLI](cli.md), ADR 0013); numeric fields may be written as
   expressions or unit strings (ADRs 0012, 0014, 0015). Only `lda-pz` is accepted. `validate` and `grid` stop without reading pseudopotentials.
2. `dft_cli::run::execute` reads the supplied UPF bytes, records SHA-256 hashes,
   and parses through `dft_engine::input::parse_upf`. `dft-upf` wraps the `upf`
   reader with supported-feature validation and unit conversion. See
   [UPF support](upf-support.md) for restrictions, including projector limits
   and accepted functional labels. Relative UPF paths use the TOML directory.
3. `ScfProblem::new` validates species, the even valence-electron count, bands
   and solver controls. It builds `PlaneWaveBasis`, its density basis, and
   `FftGridSpec`, respecting an optional grid override. `ionic_local` assembles
   reciprocal local potentials using radial transforms and structure factors;
   `FftPlan3d` transforms them to the real grid. `NonlocalOperator` assembles
   atom-centered Kleinman–Bylander projectors. Ewald ion-ion energy is prepared.
4. Starting from uniform density, `ScfProblem::step` evaluates Hartree and
   LDA-PZ potentials and updates `Hamiltonian`. Its `LinearOperator`
   implementation applies kinetic, FFT local convolution and nonlocal terms
   without constructing a dense Hamiltonian. `Davidson` uses a
   `DiagonalPreconditioner` and reuses previous orbitals as initial guesses.
5. Occupations are fixed at two electrons per occupied band and zero for empty
   bands. `density::from_orbitals` forms the output density. `energy::evaluate`
   returns `EnergyAudit`: kinetic, local, nonlocal, Hartree, XC and Ewald sum to
   the total; local zero mode, band sum, double counting, band form,
   Harris–Foulkes and the functional difference bound are audit diagnostics.
6. `LinearMixer` or `PulayMixer` proposes the next density. Charge and positivity
   safeguards, density/energy/occupied-eigen-residual gates, and an additional
   fixed-point certificate govern convergence. Energy change alone is not enough.
   `ScfProblem::run_outcome` returns `ScfOutcome` on success or budget exhaustion;
   `run` returns a success-only `ScfResult` or `ScfError`.
7. The CLI serializes the outcome and provenance after SCF returns. Its
   `RunReport` carries console text and convergence status; result and density schemas,
   output-path protection, overwrite behavior and exit statuses are specified in
   [CLI documentation](cli.md). Budget exhaustion is not reported as success.

## Public numerical interfaces

These are implemented APIs, not proposed signatures:

| Boundary | Public traits and types |
| --- | --- |
| Units/geometry | `Bohr`, `Angstrom`, `Hartree`, `Rydberg`, `ElectronVolt`; `Cell`, `Atom` |
| Basis | `GVector`, `PlaneWaveBasis`, `FftGridSpec` |
| UPF | `Pseudopotential`, `Projector`, `UpfError`; `radial::RadialMesh`, `radial::LocalTransform` |
| FFT | `Fft3d`, `FftPlan3d`, `FftError` |
| Dense algebra (`dft_numerics::linalg`) | `ComplexMatrix`, `MatrixOp`, `Orthonormalized`, `LinalgError`; functions `hermitian_eigh`, `gemm`, `orthonormalize` |
| Solvers | `LinearOperator`, `Preconditioner`, `DiagonalPreconditioner`, `DenseEigensolver`, `Davidson`, `EigenSolution`, `SolverError` |
| Mixing (`dft_solvers::mixing`) | `DensityMixer`, `LinearMixer`, `PulayMixer` |
| Physics | `hamiltonian::Species`, `AtomProjectors`, `NonlocalOperator`, `Hamiltonian` (all in `hamiltonian`); `hartree::HartreeGrid`, `xc::XcPoint`, `xc::XcGrid`, `ewald::Ion`, `ewald::EwaldOptions`, `ewald::EwaldEnergy` |
| Engine | `ScfOptions`, `Mixing`, `SpeciesPotential`, `ScfProblem`, `ScfStep`, `ScfResult`, `ScfOutcome`, `ScfError`, `IterationRecord`, `FixedPointCertificate`, `PhaseTimings`; `energy::AuditInput`, `energy::EnergyAudit` |

`LinearOperator` exposes `dimension`, `apply` and `apply_block`.
`Preconditioner` transforms residual blocks; `DensityMixer` manages density
updates and history. `DenseEigensolver` is an explicit small-system correctness
oracle, not an automatic fallback for Davidson.

faer calls are confined to `dft-numerics/src/linalg.rs`, behind owned
`ComplexMatrix` storage and functions; no faer type crosses that public boundary.
Both eigendecomposition and GEMM explicitly select `Par::Seq`.
`FftPlan3d` implements `Fft3d` with reusable ndrustfft handlers over rustfft and
serial axis transforms. Dense algebra currently has adapter functions, not a
unified `NumericalBackend` trait. A search of all workspace Rust files under
`crates/` found no `unsafe` keyword, `extern "C"`, or Accelerate reference;
the root manifest denies unsafe code. This says nothing about dependency internals.

## Verification and future work

[scripts/check.sh](../scripts/check.sh) runs workspace build, tests, formatting
and clippy. Unit and analytic tests cover unit conversions, reciprocal geometry,
basis mapping, FFT round trips/Parseval, complex eigensolves, radial transforms,
Hartree modes, XC derivatives and Ewald behavior. Independent cross-checks include
Davidson versus the dense solver, local convolution versus explicit matrices,
exchange versus test-only xcx, and the Gaussian electrostatic identity in
[engine tests](../crates/dft-engine/tests/electrostatic.rs). Synthetic fixtures
are not real-material references.

Long tests are `#[ignore]` in the debug gate. Release invocations include:

```sh
cargo test --offline --release -p dft-solvers -- --ignored --nocapture
cargo test --offline --release -p dft-engine -- --ignored --nocapture
cargo test --offline --release -p dft-cli --test run -- --ignored --nocapture
python3 -I -m unittest discover -s scripts/tests -v
```

See [ADR 0006](decisions/0006-eigensolver-api.md),
[ADR 0009](decisions/0009-scf-and-energy.md) and [CLI tests](cli.md).
The Python suite covers interchange refusals and stored real-bundle regressions.
[Validation](validation.md) records the Si/QE comparison and its scope;
[reference bundles](../tests/reference/README.md) preserve provenance. Real
Fourier-density evidence is sampled/aggregate only, and the comparison cannot
assess orthogonality from the exported files.

Later phases intend to add a library-backed PBE functional (k-point sampling, units and expression inputs are implemented, see ADRs 0011 to 0015)
([ADR 0007](decisions/0007-xc-and-special-functions.md)), and GPU backends behind
numerical interfaces ([ADR 0003](decisions/0003-prefer-general-purpose-libraries.md)).
These are intentions requiring new implementation and validation, not current
features; the existing FFT/operator/preconditioner traits provide separation,
not a ready GPU backend. MVP approval is required before expanding scope.
