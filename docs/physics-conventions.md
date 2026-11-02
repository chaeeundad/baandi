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

## 9. P3A density and kinetic energy

`dft-physics` uses `Array3<f64>` for real-space number density n in
Bohr^-3: this enforces reality and avoids storing meaningless imaginary
components. Transform inputs are explicitly converted to complex arrays.
The reciprocal density is `to_reciprocal(n)`, so n(r) = sum_G n(G) exp(iG.r)
and n(0) is the mean density, not the integrated electron count.
For each orbital, psi_b = `to_real(scatter(c_b))` / sqrt(Omega), with
sum_G |c_b|^2 = 1. Occupations are finite, nonnegative electron counts
(2 for a doubly occupied orbital). The implementation evaluates
n(r) = sum_b f_b |psi_b(r)|^2 without normalizing supplied coefficients.
The caller supplies normalized coefficients and a compatible density grid
from the same cell. Thus (Omega/N) sum_r n(r) = sum_b f_b.

T = sum_b f_b sum_G (|G|^2/2) |c_b(G)|^2, in Ha. Grid and coefficient
count mismatches and invalid numerical inputs are typed errors. Pointwise
routines perform no I/O, and pointwise XC allocates no memory.

## 10. P3A Hartree and full-grid frequencies

V_H(G) = 4*pi*n(G)/|G|^2 for nonzero G, and V_H(0) = 0 explicitly.
E_H = (Omega/2) sum_{G != 0} 4*pi*|n(G)|^2/|G|^2, in Ha.
V_H(r) = `to_real(V_H(G))`. Every represented FFT frequency participates;
there is no density-sphere filter. Even-grid Nyquist indices use the
negative frequency rule in section 4. Components outside the density
sphere are the caller's responsibility. On skew cells with occupied mixed
Nyquist modes this canonical kernel need not preserve conjugate symmetry;
therefore the potential is returned as a complex array without discarding
imaginary components or silently symmetrizing the kernel. Alias-free
physical densities supported inside the density sphere avoid this issue.

The divergent G=0 contributions of Hartree, local pseudopotential and
Ewald must later be handled together through the alpha-Z term. No such
assembly is implemented in P3A. The real-space energy test compares E_H
with (Omega/(2N)) sum_r [n(r)-mean(n)] V_H(r), including a check on the
imaginary residual. Subtracting the mean explicitly removes the G=0 term;
since V_H has zero mean this equals (1/2) integral n V_H to roundoff.
The Poisson test applies the exact spectral multiplier |G|^2 and compares
with 4*pi*[n(r)-mean(n)], with no finite-difference Laplacian.

## 11. P3A unpolarized LDA-PZ

This is Slater exchange (alpha = 2/3) plus PZ81 correlation, the combination
called `sla+pz` by QE. All energies and potentials below are in Hartree;
eps is energy per electron, whereas n*eps is energy per volume.

```text
rs = (3/(4*pi*n))^(1/3)
eps_x = -(3/4)*(3/pi)^(1/3)*n^(1/3)
v_x = (4/3)*eps_x
rs >= 1: eps_c = gamma / d, d = 1 + beta1*sqrt(rs) + beta2*rs
         gamma=-0.1423, beta1=1.0529, beta2=0.3334
         deps_c/drs = -gamma*(beta1/(2*sqrt(rs)) + beta2)/d^2
rs < 1:  eps_c = A*ln(rs) + B + C*rs*ln(rs) + D*rs
         A=0.0311, B=-0.048, C=0.0020, D=-0.0116
         deps_c/drs = A/rs + C*(ln(rs)+1) + D
```

Since drs/dn = -rs/(3n), differentiating n*eps_c gives
v_c = eps_c - (rs/3)*deps_c/drs. Hence v_xc = d[n*(eps_x+eps_c)]/dn.
For n <= 1e-10 Bohr^-3, all eps and v components are set to zero.
This is our explicit vanishing-density choice intended to match QE;
P6 must verify it. Non-finite n and n < -1e-14 Bohr^-3 are typed errors.
The closed interval [-1e-14,0] Bohr^-3 is treated as zero roundoff density
in XC and Hartree. No positive-density threshold is applied in Hartree.

Grid evaluation returns v_xc, E_xc = (Omega/N) sum_r n*eps_xc, and the
separate double-counting integral (Omega/N) sum_r n*v_xc, in Ha.
Rounded PZ81 constants leave small finite jumps at rs=1; they are not
smoothed. Tests calculate the one-sided limits directly from the listed
parameters, report measured jumps, and check their size below 1e-4 Ha.
Finite differences of n*eps_xc use a relative density step of 1e-5,
including rs=0.999 and 1.001; stencils do not cross the discontinuous join
or the density threshold, where differentiability is not claimed.

## 12. P3A Ewald point charges and error control

For charges Z_I at positions r_I = A*s_I, with a uniform compensating
background and eta in Bohr^-1, the energy in Ha is

```text
E_real = (1/2) sum_IJ Z_I Z_J sum_R' erfc(eta*|r_I-r_J+R|)/|r_I-r_J+R|
E_recip = (2*pi/Omega) sum_{G != 0} exp(-G^2/(4*eta^2))/G^2
          * |sum_I Z_I exp(iG.r_I)|^2
E_self = -eta/sqrt(pi) * sum_I Z_I^2
E_background = -pi*(sum_I Z_I)^2/(2*Omega*eta^2)
E = E_real + E_recip + E_self + E_background
```

The prime excludes I=J,R=0. Distinct ions at identical wrapped positions
are errors. Fractional positions are wrapped into [0,1) before summation;
charges may have either sign and must be finite. `libm::erfc` supplies the
special function. Default eta = sqrt(pi)/Omega^(1/3); callers may set any
finite positive eta. Real and reciprocal contributions are returned
separately for audit. **P6 note:** compare this combined energy and its
background/zero-mode convention with QE's reported "ewald contribution"
using real run provenance; no agreement with QE is claimed here.

The default absolute truncation budget is 1e-13*Qabs^2/L Ha, with
Qabs=sum_I |Z_I| and L=Omega^(1/3). It is relative to an electrostatic
scale, not to E, which may vanish by cancellation. The returned tail bound
excludes floating-point summation error. Let rho be half the sum of the
three lattice-column lengths, so a centered primitive parallelepiped
covers every point within rho of its lattice site. Let rho_g be the same
quantity for reciprocal columns and Omega_g=(2*pi)^3/Omega.
For spherical cutoffs R=u+2*rho and Gmax=v+2*rho_g, positive u,v give
conservative omitted-energy bounds per Qabs^2:

```text
B_real = sqrt(pi)/(Omega*eta^3*u) * (1+rho/u)^2 * exp(-eta^2*u^2)
B_recip = (2*pi/Omega) * (8*pi*eta^2/(Omega_g*v))
          * (1+rho_g/v)^2 * exp(-v^2/(4*eta^2))
```

To derive these bounds, tile space by centered primitive cells (shifted
by r_I-r_J for real space). A cell associated with an omitted site lies
outside cutoff-rho, and its site's distance is at least |x|-rho.
Replace the monotone radial kernel by its value there and integrate over
that exterior region. Substitution w=|x|-rho gives lower limit
cutoff-2*rho and a radial factor (w+rho)^2. Use
erfc(eta*w) <= exp(-eta^2*w^2)/(sqrt(pi)*eta*w), bound
(1+rho/w)^2 by its lower-limit value, and use
integral_u^infinity exp(-a*w^2) dw <= exp(-a*u^2)/(2*a*u).
For reciprocal space also use |structure factor| <= Qabs.

Set t=sqrt(-ln(tolerance)), u=t/eta, v=2*eta*t, then increase t by 1
until B_real+B_recip <= tolerance/L. Enumeration uses duality bounds:
|m_i+delta_s_i| <= |b_i|*R/(2*pi) in real space and
|m_i| <= |a_i|*Gmax/(2*pi) in reciprocal space. Wrapped differences have
|delta_s_i| < 1, so real bounds include 1 for the displacement. Ceil each
bound and add one integer for rounding protection, enumerate the box,
and retain only vectors inside the spherical cutoff. This holds for
skew and left-handed cells. Unrepresentable bounds are typed errors;
extremely skew cells or extreme eta can be expensive. No special
configuration-dependent energy shortcuts are used.

## 13. P3A validation tolerances

All P3A tests run in the debug gate and print values with `--nocapture`.
The deterministic pseudo-random fixture uses four normalized complex
orbitals on fcc Si (a=5.43 Angstrom, Ecut=10 Ry).

| Check | Bound |
| --- | --- |
| Density electron count | relative < 1e-12 |
| Single plane wave and weighted-superposition kinetic energy | relative < 1e-14 |
| Uniform Hartree energy and maximum potential | absolute < 1e-14 Ha |
| Cosine Hartree energy, potential, integral, linearity, spectral Poisson | relative/max-amplitude < 1e-12 |
| Exchange closed form and xcx LDA_X (n=1e-6 through 10) | relative < 1e-12 |
| XC energy-density central derivative (rs=0.1 through 20) | relative < 1e-7 |
| PZ81 branch jumps vs parameter-derived one-sided limits | absolute < 1e-10 Ha |
| PZ81 vs xcx PW92 correlation | absolute < 0.005 Ha/electron (different fits; sanity only) |
| Uniform XC energy and integral n*v_xc | relative < 1e-12 |
| Ewald eta variation, translation, scaling (0.8 and 1.3) | relative < 1e-10 |
| Ewald doubled-cell replication | relative < 1e-9 |
| Ewald cutoff refinement (1e-13 to 1e-16) | sum of tail bounds + 1e-13 Ha roundoff |

The Si Ewald energy print is informational, not a reference benchmark.
Run `scripts/check.sh dft-physics`, `scripts/check.sh`, and
`cargo test --offline -p dft-physics -- --nocapture`.

Additional tests exercise odd-grid signed modes and even-grid mixed
Nyquist modes on the fcc cell (spectral coefficients and Poisson residual
< 1e-12 relative), and a neutral pair in equivalent cubic and skew,
left-handed lattice representations (< 1e-10 relative). Doubling all
charges checks quadratic Ewald scaling (< 1e-12 relative).
