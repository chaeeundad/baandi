# UPF v2 NC support (P2A)

`dft-upf::parse_upf(&str)` and `from_file(path)` read through the vendored
`upf` 0.1.1 crate and return an immutable `Pseudopotential`. No real Si file,
other real pseudopotential, or QE calculation has been used. No QE source
was read. The format reference is the public specification distributed at
`../rust-dft-vendor/upf/docs/reference/upf-spec.html`; the reader's MIT/Apache
Rust model was inspected to establish its field coverage.

## Supported subset and reader limits

- UPF versions 2.0, 2.0.0, 2.0.1; exactly one element label, NC fully nonlocal
  form; scalar or nonrelativistic. Fake `Xx` is allowed in synthetic tests.
- Non-PAW, non-ultrasoft, no spin-orbit, no nonlinear core correction, no bare
  Coulomb dataset, no GIPAW or all-electron reconstruction. True flags are
  explicitly rejected, as are their unsupported sections even with false flags.
- At most three projectors, with l=0..3. The vendored reader enumerates only
  PP_BETA.1 through .3 (and PP_CHI.1 through .3). Larger files are explicitly
  unsupported; this is a reader limit, not a physical restriction of NC.
- `functional` must be nonempty and is stored/reported without classifying it.
  SUPPORTED means format support; it does **not** certify LDA compatibility.
  The engine must select/reject the functional before using the dataset.
- PP_LOCAL, PP_MESH, PP_NONLOCAL and PP_RHOATOM are required by the reader.
  The reader also requires its documented header and mesh attributes; absence
  of such reader-required metadata may produce a Parse error. Projectors need
  `index`, `angular_momentum`, `cutoff_radius`, `cutoff_radius_index`.
- A supplemental balanced-markup scanner reads `functional` and
  `cutoff_radius_index`, which the reader omits. It checks data-tag nesting,
  duplicate attributes, truncation, multiple roots and unsupported sections.
  It skips comments, processing instructions and informational PP_INFO contents,
  handles quoted attributes and XML character references, and rejects DTDs.
  It is not a replacement numerical-array/XML deserializer.
- The original reader model is available via `raw()`, in original UPF units;
  atomic density, optional pseudo-wavefunctions and generation metadata are
  preserved there. It is provenance/auxiliary input, not a numerical-kernel API.
  The reader itself does not retain every informational XML attribute. The own
  type separately retains the functional and cutoff counts that affect physics.

The input is expected to use the nominal UPF convention: stored beta is
`u(r)=r*beta(r)` in Bohr^-1/2 and D is in Ry. PP_LOCAL and D are converted
once with `dft_units::Hartree::from(Rydberg(x))`; u, r and rab are unchanged.
Explicit `units` attributes must match `Ry`, `Bohr`, or `Bohr^-1/2` as
appropriate. The specification warns of converters using energy-scaled beta
and reciprocal-energy D. There is no reliable discriminator for unlabelled
alternative scaling: no heuristic inference is attempted. Such producers
require a separate audited adapter before engine use. Non-nominal explicit
unit annotations are rejected.

No checksum implementation or dependency was added. Compute SHA-256 of the
original input bytes in the calling script (for example `shasum -a 256 file`);
`inspect-upf` does not claim or print a checksum.

## Validation and typed errors

`UpfError` variants are:

| Variant | Meaning / identifying payload |
| --- | --- |
| Io | File access or UTF-8 reading error |
| Parse | Reader XML/number error, malformed or truncated document |
| Version(version) | Unsupported UPF version |
| Unsupported(feature) | PAW, ultrasoft, spin-orbit, core correction, bare Coulomb, GIPAW, all-electron wavefunctions, paw_as_gipaw, non-NC pseudo_type, relativistic mode, unsupported section, l>3, projector count>3, or explicit non-nominal units |
| Missing(field) | Missing supplemental metadata, PP_LOCAL, required numbered PP_BETA, or other explicitly inspected field |
| Invalid { field, reason } | Invalid flag, label, charge, dimensions, samples, indices, section placement, attribute duplication, or cross-l coupling |
| NonSymmetricD | max absolute pair discrepancy exceeds 1e-12 Ha |

Validation checks finite positive z_valence; at least two mesh points;
finite nonnegative strictly increasing r; finite positive rab; equal mesh
array lengths matching both header and PP_MESH declarations; finite local
potential and atomic density; optional pseudo-wavefunction lengths, finiteness
and nonnegative finite occupations; and matching orbital count. Optional array
size attributes for mesh, local, density, D and beta are checked against actual
lengths. The vendor's mesh length validation errors are mapped to Invalid with
PP_R, PP_RAB, PP_LOCAL or PP_RHOATOM as field identifiers.

Projector tag and index must agree and follow consecutive one-based order;
l must not exceed header l_max; cutoff radius must be finite/nonnegative;
cutoff index is the count of valid points (1..mesh_size). Arrays may have this
length or the full mesh length. In the latter case samples after the cutoff
must be exactly zero. All samples must be finite; at r=0, u must be exactly
zero for a regular beta. The validated type stores only the valid prefix and
retains l and cutoff count. D has exactly n_proj^2 finite real entries, in
row-major order, with symmetry checked at **absolute 1e-12 Ha**; it is not
silently symmetrized. Couplings between different l must be exactly zero.

`RadialError::Mesh` identifies invalid mesh inputs and
`RadialError::Argument` identifies invalid arguments, dimensions, nonfinite
integrands, or unrepresentable results. G=0 is rejected by `at_g`; use
`alpha_z`. Negative wave numbers, nonpositive volume/charge, out-of-range l,m,
and zero harmonic directions are errors.

## Numerical API and conventions

See sections 14–18 of [physics-conventions.md](physics-conventions.md) for
formula derivations, angular phases, G=0 bookkeeping and the exact Parseval
normalization. The public radial module supplies `RadialMesh`, Simpson
integration, `spherical_bessel`, `real_spherical_harmonic`, `LocalTransform`
and `projector_transform`. Convenience methods on the validated potential
bind transforms to its own arrays. Direct evaluation is used; no tables or
interpolation are introduced. Kernels perform no I/O and contain no unsafe.

`LocalTransform::new` uses an optional r_max, default min(mesh end, 10 Bohr),
and reports `effective_r_max()`: the last mesh point <= r_max. It does not
extrapolate an off-mesh endpoint. Explicit r_max must lie within the mesh and
include at least two points. Quadrature begins at the first mesh radius;
if it is positive, the unsampled origin interval is omitted, not invented.
The caller must check convergence of this omission and of the short-range
upper tail. For bounded V, its local short-range omission is O(r_min^3).
Tests use r_min=0. No universal error bound for arbitrary supplied meshes or
poorly resolved potentials is claimed. The analytic Coulomb subtraction
extends to infinity regardless of the numerical cutoff.

## Synthetic verification and measured errors

`tests/common/mod.rs` generates a 601-point logarithmic mesh and an in-memory
UPF string labelled SYNTHETIC: Xx, Z=4, two s and one p projector, a symmetric
D with nonzero s–s coupling, and an erf-screened Coulomb local potential in Ry.
Round-trip assertions check every mesh/projector sample and exact Ry/Ha
conversion. Negative tests assert the error variant **and identifying payload**
for the requested unsupported features, missing fields, invalid meshes,
nonfinite local data, D dimensions/asymmetry, version and truncation.
The CLI test writes its own 401-point synthetic string to a temporary
`.xml` file and removes it; no UPF dataset is checked in.

Commands (PATH includes `/opt/homebrew/opt/rustup/bin`):

```sh
scripts/check.sh
cargo test --offline -p dft-upf -p dft-cli -- --nocapture
```

All tests are fast debug tests; the entire radial test binary took about
0.83 seconds in the measured run, so no ignored/release-only test is needed.
Measured errors below are mathematical synthetic comparisons, **not QE references**.

| Check | Required bound | Measured error |
| --- | --- | --- |
| Simpson exponential, n=101 / 201 / 401 / 801 / 802 | <2e-7 Bohr; n>=801 <1e-10 Bohr | 2.209e-9 / 1.383e-10 / 8.644e-12 / 5.403e-13 / 5.331e-13 Bohr |
| Simpson Gaussian, same n | same bounds | 2.518e-9 / 1.574e-10 / 9.840e-12 / 6.150e-13 / 6.118e-13 Bohr |
| Polynomial integration, index grid, n=2..9 | absolute <1e-12 in test units | asserted; cubic for n>=3, linear for n=2 |
| Bessel closed forms, x=0.2..10000.1, away from zeros | relative <1e-13 | 7.524e-15 |
| Bessel recurrence / small-x asymptote | normalized residual / relative <1e-13 | asserted; recurrence scale is sum of absolute operands |
| Harmonic orthonormality, 12 Gauss-Legendre by 32 azimuthal nodes | absolute <1e-12 | 2.998e-15 |
| Addition theorem, 100 fixed-seed direction pairs | absolute <1e-13 | 9.992e-16 |
| Cartesian p-harmonic signs | absolute <1e-15 | asserted |
| Local Fourier, 12001-point log mesh, rc=0.3 and 0.5 Bohr, A=0 and -2 Ha, s=0.4 Bohr, Z=4, Omega=120 Bohr^3, G=0.1..12 Bohr^-1 (120 points) | relative <1e-8 | max 1.967e-12 |
| alpha-Z closed forms and direct Coulomb-subtracted integral | absolute <1e-11 Ha | closed form max 1.388e-17 Ha |
| Gaussian beta, l=0..2, s=0.4/0.8 Bohr, q=0/0.1/0.7/1/3/6 Bohr^-1, 4001 radial points | relative <1e-9 (exact zero at q=0 for l>0) | max 8.089e-13 |
| Projector Parseval, 401 linear q nodes to 12/s | relative <1e-9 | max 2.220e-16 |

Local cutoff sensitivity uses rc=1.2 Bohr, A=-0.3 Ha, s=0.7 Bohr,
G=1 Bohr^-1, Z=4, Omega=120 Bohr^3. Requested r_max 3/5/8/10 Bohr
selects effective radii 2.99926634/4.99974568/7.99956021/9.99538900 Bohr.
Absolute Fourier errors are 2.017318e-6, 1.974317e-10, 5.551115e-17,
and 5.551115e-17 Ha, respectively (final bound <1e-11 Ha).
At high G, true transforms can be exponentially tiny: relative accuracy is
conditioned by cancellation/roundoff. The reported relative bounds apply to
the stated analytic cases, not all potentials and all wave numbers.

## CLI

`dft-cli inspect-upf <file>` prints element, valence charge, functional label,
NC type, mesh count/range in Bohr, each projector's l and cutoff index, and
the full small D matrix in Ha. It reports SUPPORTED on successful validation.
Failure exits nonzero with the path and the first explicit unsupported feature
or invalid field. Error reporting is fail-fast; it does not claim to enumerate
all problems in a malformed input. Existing `validate` and `grid` are intact.

The added CLI path dependency requires Cargo's corresponding dft-cli lockfile
entry. Cargo generates that entry during tests; the ticket forbids editing
Cargo.lock, so it is restored afterward. The reviewer must include the
one-line dependency entry when committing the manifest change. No dependency
version changes are required. Artifact JSON/session registry files are outside
this ticket's owned files and are left to the orchestrator.
