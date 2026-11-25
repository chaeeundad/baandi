# UPF v2 NC support (P2A, P2C)

`dft-upf::parse_upf(&str)` and `from_file(path)` read through the vendored
`upf` 0.1.1 crate and return an immutable `Pseudopotential`. The reader and the
radial mathematics were developed on synthetic fixtures only; the real
`Si.pz-vbc.UPF` (provenance unverified, see candidate-si-upf.md) was later read
and used for validation (docs/validation.md R1 to R5), which exposed the three
real-world deviations handled in P2C and P5C. The file was never given to the
implementation agent. No QE source was read. The format reference is the public specification distributed at
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
  The shared `dft-upf::functional_is_lda_pz` helper is used by the engine and CLI
  to accept exactly the case-sensitive token sequences `SLA PZ` and
  `SLA PZ NOGX NOGC`. It splits on ASCII whitespace; leading, trailing and repeated
  whitespace is ignored for compatibility only. `NOGX NOGC` explicitly means no
  exchange or correlation gradient correction. All other sequences are rejected
  by SCF validation. The original label, including whitespace, is retained in
  the parsed potential and CLI `results.json`.
- PP_LOCAL, PP_MESH, PP_NONLOCAL and PP_RHOATOM are required by the reader.
  The reader also requires its documented header and mesh attributes; absence
  of such reader-required metadata may produce a Parse error. Projectors need
  `index`, `angular_momentum`, and either `cutoff_radius_index` or
  `cutoff_radius` (see the support/tail policy below).
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

## Real-world generator compatibility (P2C)

Two producer conventions were missing from the original synthetic coverage:

- The specification lists `scalar | full | nonrelativistic`, but older
  generators also emit `relativistic="no"`. The reader accepts `no` as a
  synonym of `nonrelativistic`, while retaining the original spelling in
  `raw().header.relativistic`. `scalar` remains supported. `full` and every
  unknown value fail with `unsupported feature: relativistic=<value>`.
- Real numeric sections wrap rows, commonly with `columns="4"`. The vendored
  reader's serde decoder for mesh/local/density arrays separates tokens on
  spaces; literal newlines can join two numbers into an invalid token.
  Before calling the reader, we replace each ASCII whitespace character
  (space, tab, LF, CR, vertical tab, form feed) with one space. This supports
  LF and CRLF, mixed separators and runs of whitespace in PP_R, PP_RAB,
  PP_LOCAL, PP_BETA.*, PP_DIJ, PP_CHI.* and PP_RHOATOM. `columns` is a layout
  hint, not a size constraint; the final row may be short.

This is safe for supported numerical data because whitespace only separates
numbers: no non-whitespace character, numeric spelling, sign, exponent or unit
is changed, and no values are inferred. Non-ASCII characters are unchanged.
Supplemental metadata validation uses the original text. The reader receives
normalized whitespace throughout the document, including informational text
and attribute values, so `raw()` does not preserve their original whitespace
formatting. Each replacement occupies one byte: byte offsets remain aligned
with the original input, but reader line/column locations refer to the
normalized single-line text and must not be treated as original locations.
Typed validation errors retain their field names and identifying payloads.

Synthetic tests compare all arrays bit for bit (zero tolerance, including
Ry-to-Ha converted fields), with three optional PP_CHI arrays, 601 mesh
samples and nine D entries; wrapped final rows contain one entry. Separate
variants cover the accepted relativistic spellings, LF, CRLF, tabs, double
spaces, CR, form feed, vertical tab and mixed separators. Existing negative
and radial-math tests remain in place. A process-boundary test under
`dft-upf/tests/parser.rs` invokes `cargo run --offline -p dft-cli -- inspect-upf`
on a temporary wrapped CRLF, `relativistic="no"` fixture and checks SUPPORTED.
No real UPF dataset or independent real-file result was available for this
implementation; these are synthetic compatibility checks only.

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
l must not exceed header l_max; a supplied cutoff radius must be finite and
nonnegative. The declared cutoff index is metadata (1..mesh_size).
Arrays may have this length or the full mesh length. All samples, including
tails beyond the declared index, must be finite. A finite residual at r=0 is retained: it has
zero weight in the projector integral. The validated type stores only the
effective support and retains l, the declared index, effective support index,
and largest absolute tail sample. D has exactly n_proj^2 finite real entries, in
row-major order, with symmetry checked at **absolute 1e-12 Ha**; it is not
silently symmetrized. Couplings between different l must be exactly zero.

### Full stored beta support (G7B-UPF-FULL-SUPPORT)

History: G7-UPF-BETA-TAIL introduced metadata-based truncation with a printed
precision tolerance. G7B supersedes that rule: dev08 gen26 found real SG15
Si_ONCV_PBE-1.0 and O.pbe-hgh tails that it rejected. No tolerance or
printed-precision estimate now determines acceptance or support.

The [public UPF description](https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format)
defines the stored function `u(r)=r*beta(r)`; the
[UPF v2 schema documentation](https://orange-box.cc/qeschemas/qe_pp-v2.html)
lists the optional cutoff metadata. BAANDI's explicit NC reader policy is:

- Record `cutoff_radius_index` when present (one-based). If only a radius is
  given, infer the declared index from the last mesh point with
  `r <= cutoff_radius`, without interpolation. An index takes precedence over
  a conflicting radius. Missing both, an index outside `1..mesh_size`, or an
  array of neither declared-index nor full-mesh length remains an error.
- Reject non-finite samples anywhere, including beyond the declared index.
  All finite tail values are accepted, irrespective of amplitude, sign,
  printed precision, generator, or whether the declared prefix is all zero.
- The effective support index is exactly
  `max(declared_index, last_nonzero_sample_index)`. An all-zero array uses
  the declared index. Positive and negative exact zeros are both zero;
  even the smallest nonzero subnormal extends support.
- Store and integrate the unmodified samples through this effective index.
  No nonzero sample is truncated, thresholded, interpolated, or renormalized.
  `raw()` retains the entire original array, including zero padding.
  Files without nonzero tails keep their previous support and numerical bits.
- `Projector::cutoff_index()` and `effective_support_index()` both return
  the effective length; `declared_cutoff_index()` records the metadata-derived
  index. `largest_tail_magnitude()` is `max(abs(u[declared_index..]))` in
  Bohr^-1/2 (zero for absent or exactly zero tails).
  `inspect-upf --json` exposes an ordered `projector_support` array with
  `index`, `l`, `declared_cutoff_index`, `effective_support_index`, and
  `largest_tail_magnitude_bohr_minus_half` for each validated projector.

Both production and strict drivers, and the nonlocal and derivative setup,
consume the same effective sample slice through the existing radial APIs.
No physics-kernel or quadrature formula changes are needed. Integration ends
at the effective index; extending across extra zero padding is not required
and could change endpoint quadrature weights on earlier nonzero samples.

Synthetic tests in `crates/dft-upf/tests/all/beta_tail.rs` cover ONCV-like
smooth sign-changing continuations and HGH-like Gaussian tails, odd and even
full-support lengths, and both quadrature rules. An independent full-mesh
weight-vector quadrature with `j0(x)=sin(x)/x` must agree within absolute
1e-14 Bohr^3/2 at q=0, 0.1, 1, and 8 Bohr^-1. Separate tests check unchanged
bits for zero-tailed versus shortened files, radius/index precedence, exact
zero padding, subnormal tails, zero prefixes, and non-finite rejection at
every sample position. These fixtures are synthetic, not copied UPF data.

G7B local validation (2026-11-25, macOS/aarch64): the maximum independent
full-mesh quadrature difference was `5.32907051820075139e-15 Bohr^3/2`.
Both stored Si.pz-vbc baselines (reference and strict oracle,
Gamma/k222s/k222, three repeats each) passed `bench_scf.py
--require-identical`: all six energy differences were exactly zero Ha and
all density SHA-256 hashes matched. A synthetic three-projector CLI JSON
check independently verified declared indices `[90, 90, 90]`, effective
indices `[101, 95, 90]`, and tail maxima
`[2.8818320503339954e-17, 1.3494495090618492e-8, 0] Bohr^-1/2`.
`reference_identity` passed for both drivers. `cargo test -p dft-upf`
passed 25 tests with one external-file test ignored; the workspace test gate
passed 1,023 tests with zero failures and 231 ignored.
`CHECK_SCCACHE=0 BAANDI_TEST_TMP=/private/tmp/g7b-tests scripts/check.sh`
passed build, test, fmt, Clippy, production-check and production-clippy
(exit 0). `cargo check --workspace --tests` and workspace all-target Clippy
with warnings denied also passed. The sandbox denied the optional `ps` load
probe and `nice` priority adjustment; neither affected the gate results.
Generated README
figures were refreshed without a diff and all nine freshness tests passed.
Real SG15/HGH files were unavailable locally; their dev08 qualification is
not claimed by these synthetic checks. No new QE reference run was performed.

An optional external SG15 Si2 smoke compares the production reference driver
with the strict oracle at the same 10 Ry cutoff and PBE functional, requiring
verification and `abs(delta E) < 5e-9 Ha` (standard class). It is ignored in
debug gates, requires an external file, and is not a QE comparison or a
cutoff-converged physical reference:

```sh
DFT_ONCV_UPF=/absolute/Si_ONCV_PBE-1.0.upf cargo test --release \
  -p dft-engine --test all oncv_beta_tail -- --ignored --nocapture
```

Without the variable it looks under the worktree's `../pseudo`. Missing files
fail explicitly when this ignored test is requested; no UPF is committed.

Historical G7 validation before the G7B support change
(2026-11-25, macOS/aarch64, Rust 1.99.1): the same ignored
test passed with `DFT_ONCV_UPF=$PWD/../pseudo/pseudodojo/ex/Si.upf`, SHA-256
`39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282`.
The reference energy was `-7.80444120941249153 Ha`, strict oracle energy
`-7.80444120970207322 Ha`, difference `2.89581691959028831e-10 Ha`.
This available PseudoDojo file has six projectors with exactly zero tails;
it establishes existing-file compatibility and exercises the SCF smoke,
not SG15/HGH nonzero-tail qualification. Those files were not available.
The nonzero-tail evidence is synthetic. Both stored Si.pz-vbc baseline
comparisons (reference and strict oracle, Gamma/k222s/k222, three repeats)
passed `bench_scf.py --require-identical`: all six energy differences were
zero Ha and all density SHA-256 hashes matched. Generated README figures
were refreshed without a diff; the nine figure/evidence tests passed.
`cargo test -p dft-upf` passed 20 tests (one external-file test ignored).
`CHECK_SCCACHE=0 BAANDI_TEST_TMP=/private/tmp/g7-upf-beta-tail scripts/check.sh`
passed all gates: 1,022 tests passed, none failed, 231 ignored, including
the production-only checks. `cargo check --workspace --tests`, formatting
and workspace all-target Clippy with warnings denied also passed. The sandbox
denied the gate's optional `ps` load probe; this did not affect any gate.

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

### P7B functional compatibility

The engine selects LDA-PZ or PBE explicitly. PBE recognizes `PBE` and the
four-token `SLA PW PBX PBC` label and its older spelling `SLA PW PBE PBE` (for example Si.pbe-rrkj.UPF), normalized for whitespace and case.
Other PBE-family labels (PBE0, PBEsol) are not PBE. LDA matching is unchanged.
`xc.allow_functional_override = true` permits a mismatch and records both the
requested functional and unchanged UPF header, plus the permission/applied
flags in results.json and the console. It does not bypass any UPF structural
restriction (NC, no NLCC, spin-unpolarized). Parser acceptance alone is not
functional compatibility.
