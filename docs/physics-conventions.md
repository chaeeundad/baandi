# Physics conventions

## 1. Units

Internally use Hartree atomic units: lengths in Bohr, energies in Hartree
(Ha), with electron mass, elementary charge, reduced Planck constant and
4*pi*vacuum permittivity equal to one. Convert only at I/O boundaries.

`dft-units` supplies distinct `Bohr`, `Angstrom`, `Hartree`, `Rydberg`, and
`ElectronVolt` newtypes with explicit `From` conversions:

- `RYDBERG_PER_HARTREE = 2.0`: 1 Ha = 2 Ry exactly.
- `BOHR_IN_ANGSTROM = 0.529177210903`: multiply Bohr by this to obtain Angstrom.
- `HARTREE_IN_EV = 27.211386245988`: multiply Ha by this to obtain eV.

The latter two values are the CODATA 2018 central values from the
[NIST 2018 constants table](https://physics.nist.gov/cuu/pdf/all_2018.pdf).
Reverse conversions divide by the same constants. These measured physical
constants are not exact definitions. P6 must confirm which constants QE
uses before comparisons; equality with QE is not established.

Geometry uses explicitly unit-named methods with plain arrays to avoid
exposing matrix-library types. Its test's fcc lattice parameter is a fixed
converted Bohr fixture; production conversion constants live only in dft-units.

## 2. Lattice and positions

A has lattice vectors a1, a2, a3 as **columns**, in Bohr. Array matrix entries
are `[Cartesian component][vector index]`. User-facing `lattice_vectors`
has one vector per **row**. `Cell::from_user_rows_bohr` transposes that
representation into internal columns; `from_matrix_columns_bohr` instead
accepts conventional matrix entries directly. Neither constructor parses TOML.

Volume Omega = |det A| in Bohr^3. `CellError::NonFinite` rejects non-finite
entries. With s = max_ij |A_ij|, `NearSingular` rejects s = 0 or
|det(A/s)| <= `MIN_SCALED_DETERMINANT = 1e-12` (dimensionless).
`Unrepresentable` rejects non-finite or zero computed volume, or a non-finite
inverse or reciprocal matrix. This scale-independent determinant threshold
is not a condition-number bound. Left-handed cells (det A < 0) are accepted.

B = 2*pi*inverse(A)^T in Bohr^-1, so a_i dot b_j = 2*pi*delta_ij even for
left-handed cells. `reciprocal_columns_inverse_bohr` returns B with the same
matrix indexing as A. Cartesian position r = A*s_frac; fractional position
s_frac = inverse(A)*r. `fractional_to_cartesian_bohr` and
`cartesian_bohr_to_fractional` do not wrap positions. `Atom` stores a species
label and a fractional position without interpreting or validating the label;
position values are not validated or wrapped.

## 3. Plane waves

For the Gamma-only complex representation:

```text
psi(r) = Omega^(-1/2) * sum_G c(G) exp(+i G.r)
sum_G |c(G)|^2 = 1
```

Here G is an integer combination of reciprocal lattice columns. The FFT
adapter does not apply the wavefunction's Omega^(-1/2) prefactor.

## 4. FFT grid and convention

Dimensions are (n1,n2,n3), N = n1*n2*n3. Arrays are standard row-major
(C order) `ndarray::Array3<Complex64>` with index (i,j,k), flat index
`(i*n2 + j)*n3 + k`, and grid point
r = (i/n1)*a1 + (j/n2)*a2 + (k/n3)*a3.

```text
to_real:       f(r) = sum_G f(G) exp(+i G.r)          (no scale factor)
to_reciprocal: f(G) = (1/N) sum_r f(r) exp(-i G.r)
Parseval:     sum_r |f(r)|^2 = N * sum_G |f(G)|^2
```

For an axis of length n, index m represents frequency m if
`m <= (n-1)/2` (integer division), and m-n otherwise. For even n, the
Nyquist index n/2 represents **-n/2**. Positive and negative Nyquist waves
are indistinguishable on the sampled grid; the negative sign is the chosen
label. P1B will implement basis enumeration/index mapping.

The measured ndrustfft default is forward unscaled with exp(-i...), inverse
scaled by 1/N with exp(+i...) (see `backend-spike-results.md`). The adapter
calls `ndfft` along all three axes and divides by N for `to_reciprocal`,
and calls `ndifft` along all axes and multiplies by N for `to_real`.
Constant-function results are guaranteed within the stated floating-point
tolerances, not bit-exact. All inputs use the three-axis FFT and documented
scaling.

## 5. Plan lifecycle and verification tolerances

`Fft3d` is the small backend trait. `FftPlan3d` owns three ndrustfft handlers
and two full-grid scratch arrays allocated once. Calls reuse them and write
into caller-owned output. No per-call plan construction occurs. ndrustfft
may allocate internal lane scratch; allocation-free execution is not claimed.
Zero dimensions, size overflow, and sizes whose Complex64 storage exceeds
isize::MAX return `InvalidGrid`; allocation failure is not recoverable.
Shape mismatches or nonstandard input/output layouts return
`IncompatibleArray` before modifying output.

All reported errors are dimensionless unless specified. No tolerances are
adjusted to accommodate a failing result:

| Check | Bound |
| --- | --- |
| Unit round trips (both directions) | relative error < 1e-15 |
| 1 Ha = 2 Ry | exact |
| Lattice reciprocal duality | absolute error < 1e-12 |
| Fractional position round trip | max absolute error < 1e-12 |
| Hand-computed cell volumes | relative error < 1e-12 |
| Both FFT round-trip directions | max entry error / max input magnitude < 1e-12 |
| Parseval | abs(sum_r/N - sum_G) / sum_G < 1e-12 |
| Unit-amplitude signed modes (coefficients and sampled values) | max absolute error < 1e-12 |
| Conjugate symmetry for real input | max absolute error < 1e-12 |
| Linearity, both transforms | max entry error / max expected magnitude < 1e-12 |
| Band-limited coefficients on two grids | max entry error / max expected magnitude < 1e-12; each shared peak relative difference < 1e-12 |
| Constants 1 and 0.3-0.7i on 4x5x6 and 15x20x24: G=0 value and inverse values | relative error < 1e-14 |
| Constant-function non-origin coefficients | maximum absolute magnitude < 1e-14 |
| Repeating either transform with one plan | identical arrays |

Geometry cases are cubic (10 Bohr), fcc primitive (5.43 Angstrom converted
to Bohr, volume a^3/4), hexagonal, strongly skewed triclinic, and left-handed.
Tests also reject singular, near-singular, non-finite, and unrepresentable
cells and distinguish user rows from matrix rows.

FFT round trips and Parseval use 4x5x6, 15x20x24, 16^3, and 27^3 grids.
Signed-mode tests on 4x5x6 use frequencies (1,-1,2), index (1,4,2), and
(-2,1,-3), index (2,1,3), including even-axis Nyquist indices. A two-mode
band-limited function is checked on 4x5x6 and 15x20x24. These are mathematical
fixtures, not QE reference results. All tests are intended for the fast debug
gate; run `scripts/check.sh` and
`cargo test --offline -p dft-geometry -p dft-units -p dft-numerics -- --nocapture`.

## Not yet defined

- G-vector cutoff, enumeration, and density grid rule (P1B).
- Hartree G=0 handling, ionic terms, and energy double counting.
- UPF projector conventions, XC parameterization, and SCF convergence.
- QE constant agreement and physical reference validation (P6).

## 6. P1B Gamma-only basis and density sphere

This section defines the P1B items previously listed as not yet defined.
Miller indices (m1,m2,m3) in Z^3 label G = m1*b1 + m2*b2 + m3*b3,
where the reciprocal columns B come from dft-geometry in Bohr^-1.
`PlaneWaveBasis` stores the cutoff in Hartree and each vector's Miller
indices, Cartesian components (Bohr^-1), and squared norm (Bohr^-2).
It includes the full complex Gamma sphere, including G=0 and both G and
-G; it is not the real-Gamma half sphere.

The wavefunction inclusion condition is
`0.5*|G|^2 <= Ecut_wfc*(1 + 1e-12)`. The density basis uses cutoff
`4*Ecut_wfc` with the identical relative boundary allowance, hence radius
2*Gmax_wfc. Sums and differences of two wavefunction vectors lie within
this density sphere, including the same relative radial allowance.
Cutoffs must be finite and positive; nonrepresentable bounds are errors.

Enumeration uses the duality identity m_i = a_i dot G / (2*pi).
Cauchy-Schwarz gives |m_i| <= |a_i|*Gmax/(2*pi) even for skewed cells.
The implementation uses the tolerance-expanded radius, rounds the bound
up, and adds one integer for floating-point rounding protection before
looping over the full box. Tests independently enumerate a box two larger
on every side for cubic, fcc primitive, hexagonal and strongly skewed
triclinic cells. Vectors are sorted by ascending computed |G|^2, without a
comparison tolerance, then by (m1,m2,m3) lexicographically for exact ties.
G=0 is first. A cubic reciprocal unit lattice tests the exact 0.5 Ha shell
and both sides of the boundary allowance.

## 7. P1B FFT grid and index tables

For the density sphere let M_i = max |m_i|. The required axis lengths are
2*M_i + 1, ensuring distinct sampled frequencies for all density vectors.
Automatic grid selection rounds each required length up to the next
integer having only prime factors 2, 3 and 5 (including 1). This is an
independent implementation, not a claim of agreement with QE; P6 compares.
An explicit override must meet every required length; otherwise
`GridTooSmall` is returned. Larger overrides need not be smooth.
Zero dimensions and unrepresentable storage sizes are rejected.

An axis frequency m maps to m for m >= 0 and n+m for m < 0.
The canonical allowed interval is `-floor(n/2) <= m <= floor((n-1)/2)`;
in particular even grids allow -n/2 and reject +n/2, matching section 4.
Out-of-range indices return `IndexOutOfRange`. Grid construction stores
coefficient-to-grid indices and a C-order flat grid-to-coefficient table;
entries outside the selected coefficient basis are None. Grid dimensions
and tables are immutable after validation. Supply the density basis from
the same cell and four-times wavefunction cutoff when constructing a grid.

Coefficient vectors use the G-list order. Scatter fills an
`Array3<Complex64>` with zeros and writes each c(G) to its mapped index.
Gather reads these indices in G-list order. Wrong coefficient counts and
array shapes are typed errors. No data-dependent transform shortcuts are
used. As with the FFT backend, allocation failure is not recoverable;
large but representable inputs can demand substantial time and memory.

## 8. P1B normalization verification

With sum_G |c(G)|^2 = 1, apply P1A `to_real` to the scattered coefficients
and divide the result by sqrt(Omega) to obtain psi(r). Then
`(Omega/N) * sum_r |psi(r)|^2 = 1`. Fixed-seed normalized complex test
coefficients verify scatter -> to_real -> to_reciprocal -> gather with
maximum absolute coefficient error < 1e-12 and the normalization identity
with absolute dimensionless error < 1e-12. Repeating on an enlarged grid
must preserve gathered coefficients with maximum absolute error < 1e-12.

The large-cutoff count test compares the discrete full-sphere count to
`Omega*Gmax^3/(6*pi^2)` (sphere volume divided by reciprocal-cell volume).
A cubic fixture has Miller-space radius 20; a 5% relative bound allows
finite-shell lattice-point fluctuations and is only an asymptotic sanity
check. The Si primitive fixture at a=5.43 Angstrom and Ecut_wfc=30 Ry prints
its counts and FFT dimensions and checks internal consistency only. These
are not external reference numbers. All P1B tests belong to the fast debug
gate and can be inspected with
`cargo test --offline -p dft-basis -p dft-cli -- --nocapture`.
