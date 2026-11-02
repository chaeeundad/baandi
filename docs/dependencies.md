# Direct dependencies

Exact resolved versions below were read from Cargo.lock; license expressions
were checked in the corresponding vendored Cargo.toml files. The lockfile,
workspace manifests, vendor, and Cargo source configuration are unchanged by
P1A. Use `CARGO_NET_OFFLINE=true` for verification. Rust 1.99.1 is pinned in
rust-toolchain.toml (edition 2024, resolver 3).

| External direct dependency | Locked version | License | Workspace consumers |
| --- | --- | --- | --- |
| nalgebra | 0.35.0 | Apache-2.0 | dft-geometry, dft-basis |
| thiserror | 2.0.21 | MIT OR Apache-2.0 | dft-geometry, dft-numerics, dft-basis, dft-cli |
| ndarray | 0.17.2 | MIT OR Apache-2.0 | dft-numerics |
| num-complex | 0.4.6 | MIT OR Apache-2.0 | dft-numerics, dft-basis |
| rustfft | 6.4.1 | MIT OR Apache-2.0 | dft-numerics (also used through ndrustfft) |
| ndrustfft | 0.6.2 | MIT | dft-numerics |
| clap | 4.6.7 | MIT OR Apache-2.0 | dft-cli |
| serde | 1.0.229 | MIT OR Apache-2.0 | dft-cli |
| toml | 1.1.8+spec-1.1.0 | MIT OR Apache-2.0 | dft-cli |

There are currently no inter-workspace path dependencies or build dependencies.
dft-units uses std only. faer was selected in ADR 0004 for future dense linear
algebra but is not declared in the current workspace; it is not claimed as a
current dependency. LAPACK and the Accelerate source selectors were spike-only.
The hand-written Accelerate adapter and native tests were removed in P1A.

Adoption/maintenance review: ADR 0003 records general library policy and release
checks; ADR 0004 records the measured ndrustfft/rustfft backend choice and Apple
Silicon results; ADR 0005 records nalgebra and utility-library checks. Current
foundation crates build and test on the installed aarch64-apple-darwin toolchain;
no new releases, dependencies, native installs, or vendor changes are introduced.
The historical backend spike documents internal allocation and performance
limits. This direct-dependency inventory is not a fresh transitive license audit.
The project's own LICENSE remains undecided.

## P3A dft-physics usage

No dependencies were added. Production physics glue uses the existing
path dependencies dft-units, dft-geometry, dft-basis and dft-numerics;
ndarray/num-complex hold grids and coefficients, nalgebra handles fixed
three-dimensional Ewald geometry, thiserror supplies typed errors, and
libm::erfc evaluates Ewald screening (ADR 0007). The existing xcx 0.3.1
dev-dependency cross-checks LDA_X and supplies a loose LDA_C_PW sanity
comparison; it remains test-only and is not a PZ81 oracle. License and
backend decisions remain those already recorded in ADRs 0004/0005/0007.

## P2A dft-upf usage

No external dependency was added. The predeclared `upf` 0.1.1 reader
(MIT OR Apache-2.0; ADR 0008) supplies typed UPF arrays; a supplemental
metadata scanner recovers functional and cutoff indices omitted by its model.
The reader enumerates only three numbered beta/chi entries; this limit is
explicitly documented and rejected rather than silently truncating data.
`thiserror` supplies typed validation errors, `libm` 0.2.16 supplies erf,
and `dft-units` performs Ry-to-Ha conversion. Elementary l<=3 angular
functions and index-mesh Simpson quadrature follow ADR 0008. The only
manifest addition is the ownership-permitted dft-cli -> dft-upf path edge.
No SHA-256 dependency or implementation was introduced. See `upf-support.md`
for the nominal/alternative UPF unit distinction and synthetic verification.

## P2B Hamiltonian path dependencies

`dft-physics` adds only the ownership-approved workspace path dependencies
`dft-upf` (validated radial transforms) and `dft-solvers` (LinearOperator;
solver comparisons in tests). No external crates, versions or licenses change.
The test imports the existing synthetic UPF text generator with a selectable
801-point resolution; the default UPF fixture remains 601 points. Production
numerics perform no parsing or I/O during application.

Cargo automatically adds these two path edges to its lockfile while running
checks. The ticket forbids retaining Cargo.lock edits, so the original lockfile
is restored after verification. The reviewer must include those two dependency
names in the dft-physics lockfile entry when committing the manifest change.
Artifact/session registry paths are outside this ticket's ownership; the final
structured report supplies gate evidence for the orchestrator to reconcile.

## P5A SCF and mixing

No external dependency, backend version, license or vendor changes. dft-engine
adds only the ownership-permitted path edges to dft-physics, dft-solvers,
dft-basis, dft-geometry and dft-upf. dft-solvers uses its existing dft-numerics
edge for the small Anderson Gram eigensolve (sequential faer), with a spectral
pseudoinverse cutoff. No hand-written dense linear solver is introduced.

The additive `dft_physics::data` facade re-exports existing physics API types
(ndarray arrays, complex numbers, Hartree units, FFT plans, matrix/error types)
and existing libm erf/erfc, so engine code requires no forbidden direct backend
edge. New physics methods update a Hamiltonian potential transactionally and
expose read-only nonlocal data. All previous signatures and tests remain intact.

Cargo regenerates the dft-engine lockfile dependency list during verification;
Cargo.lock is restored afterward because it is outside P5A ownership. The
reviewer must reconcile those five path dependency names when committing. No
new package versions need resolving or downloading. All verification is offline
on the existing Apple Silicon toolchain. Artifact/session JSON updates remain
with the orchestrator because those files are outside this ticket's ownership.

## P5B CLI output

No external dependency added. The CLI uses the already adopted sha2 and
serde_json (ADR 0010) and adds the permitted dft-engine path dependency. CLI tests additionally use
the permitted dft-physics path dependency for the existing libm erf re-export
in the previously verified synthetic Gaussian fixture.
UTC civil-date formatting and compile-time git/rustc capture use std only.
The engine input facade re-exports existing geometry/UPF types; sequential
faer/FFT execution remains explicit in the existing numerical backend.
