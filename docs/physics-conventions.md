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

## 14. P2A UPF boundary and radial storage

The public UPF specification shipped with `upf` 0.1.1 defines Bohr lengths,
Ry local potentials and nominal Ry D_ij. The stored projector is
u_l(r)=r*beta_l(r), with u in Bohr^-1/2, hence beta in Bohr^-3/2.
Our validated representation retains u, converting only PP_LOCAL and D
with the shared Ry-to-Ha conversion at parse time. The nonlocal energy
bilinear beta*D*beta then has Hartree units. The specification mentions
alternative converter scalings (energy-valued beta and reciprocal-energy D);
these cannot be inferred from unlabelled numbers. P2A assumes the nominal
convention, rejects explicit alternative units, and requires producer auditing
before real-data use. See `upf-support.md` for reader and validation limits.

The projector cutoff index is a count: samples 0..cutoff_index are valid
(exclusive Rust upper bound), equivalent to UPF's one-based final index.
At r>0 beta=u/r. At r=0 its finite limit is not reconstructed from a single
zero u sample; the transform uses r*u directly, whose origin value is zero.
This avoids interpolation or an arbitrary beta(0). No volume factor, Y_lm,
(-i)^l or atomic structure factor is included in the radial projector result.

## 15. P2A mesh quadrature and elementary angular functions

Treat the mesh index t as a continuous coordinate with unit step and
rab_i=dr/dt at node i. Then integral f(r)dr = integral f(r(t))rab(t)dt.
Apply Simpson's 1/3 rule to successive two-interval panels of this product.
For odd n, all n-1 intervals are covered by 1/3 panels. For even n>=4,
the last three intervals use the 3/8 rule, with weights (3/8)*(1,3,3,1);
all preceding intervals use 1/3 panels. For n=2 use the index-space
trapezoid. A one-point projector prefix has no interval and integrates to
zero. Panel sums use compensated summation. This is fourth-order for smooth
index-space integrands, except for the explicitly two-point trapezoid case.
The tests use r_i=r0*(exp(i*dx)-1), rab_i=dx*(r_i+r0), r0=0.02 Bohr;
convergence is measured rather than assuming every arbitrary mesh is resolved.

Spherical Bessel j_l follows its regular power series:

```text
j_l(x) = x^l/(2l+1)!! * sum_{k>=0} a_k,
a_0=1,  a_k/a_(k-1) = -x^2/[2k*(2l+2k+1)].
```

For x<4 the series stops once a term is below machine epsilon times the
partial sum (80 terms is a conservative cap for this range). For x>=4,
use j_0=sin(x)/x, j_1=(j_0-cos(x))/x and
j_(l+1)=(2l+1)*j_l/x-j_(l-1), through l=3. At zero the regular series gives
j_0=1 and the other j_l=0. Tests compare independently evaluated closed forms
away from zeros at relative 1e-13, and the x^l/(2l+1)!! limit for
x=1e-8,1e-20,1e-90 at relative 1e-13. Recurrence residuals are normalized
by the sum of absolute operands and bounded by 1e-13; normalizing by the
small cancelled sum at large x would measure conditioning rather than
function accuracy. No irregular y_l function is implemented.

Let complex harmonics use Condon-Shortley associated Legendre polynomials:
Y_l^m=N_lm*P_l^m(cos(theta))*exp(i*m*phi),
N_lm=sqrt((2l+1)/(4*pi)*(l-m)!/(l+m)!). Our real basis is

```text
m=0:  N_l0 P_l(cos(theta))
m>0:  sqrt(2)*(-1)^m*N_lm*P_l^m(cos(theta))*cos(m*phi)
m<0:  sqrt(2)*(-1)^|m|*N_l|m|*P_l^|m|(cos(theta))*sin(|m|*phi).
```

Within each l the order is m=-l,...,0,...,+l; flattening successive l blocks
would use index l*l+l+m. Thus l=1 is proportional to (y,z,x).
P_m^m=(-1)^m*(2m-1)!!*(1-z^2)^(m/2), P_(m+1)^m=(2m+1)zP_m^m,
and (l-m)P_l^m=(2l-1)zP_(l-1)^m-(l+m-1)P_(l-2)^m give all l<=3.
The real basis is an orthogonal change from complex harmonics, preserving
unit angular norm and the addition theorem. Quadrature uses 12 Gauss-Legendre
nodes in cos(theta) and 32 uniform azimuths, checking all 16x16 overlaps to
absolute 1e-12. The addition theorem uses 100 fixed-seed unit-vector pairs
and absolute 1e-13. Directions supplied to the API are normalized internally.

## 16. P2A local Fourier transform and finite zero mode

Angular integration of exp(-i G.r) gives 4*pi*j_0(Gr). Set
W(r)=V(r)+Z*erf(r)/r. The smooth unit-width reference Coulomb potential
has Fourier transform -4*pi*Z*exp(-G^2/4)/G^2: applying Poisson's equation
reduces the result to the Fourier transform of the normalized Gaussian
charge pi^(-3/2)*exp(-r^2). Thus for G>0,

```text
Vloc(G) = (4*pi/Omega) integral r^2 W(r) j_0(Gr) dr
          -4*pi*Z*exp(-G^2/4)/(Omega*G^2).
```

The numerical integrand is formed as r^2*V+Z*r*erf(r), finite at r=0.
libm supplies erf. The analytic tail extends to infinity. The numerical
upper bound defaults to min(mesh end,10 Bohr), rounded down to a mesh point;
`effective_r_max` exposes that choice. No table or interpolation is used.
Mesh integration starts at r_min; if positive, the omitted origin interval
requires convergence checking. Short-range tail truncation also requires
convergence; it is not the truncation of the long-range Coulomb tail.

Define the returned alpha-Z quantity C (in Ha) precisely by

```text
C = (4*pi/Omega) integral_0^infinity r^2 [V(r)+Z/r] dr
  = (4*pi/Omega) integral_0^infinity r^2 W(r) dr + pi*Z/Omega.
```

The last equality follows from integral_0^infinity r*erfc(r)dr=1/4:
insert erfc(r)=(2/sqrt(pi))*integral_r^infinity exp(-t^2)dt, interchange
integrals, and use the Gaussian second moment. Expanding the analytic tail
also gives -4*pi*Z/(Omega*G^2)+pi*Z/Omega+O(G^2), establishing the same
finite remainder. `alpha_z()` therefore adds pi*Z/Omega to the numerically
truncated W integral, including the reference's finite remainder exactly.
It equals the infinite-domain definition up to the stated numerical
quadrature/origin/short-range-tail omissions. It is not merely integral W.

This defines the local zero-mode constant consistent with dropping the
bare Coulomb 1/G^2 divergence. Section 10 sets V_H(0)=0; section 12 specifies
the Ewald compensating-background convention. In later assembly, retain
C for each ion (the local G=0 potential is their sum); its electron–ion
energy contribution is N_e*sum_I C_I. Combine with the existing Hartree/Ewald
zero-mode/background terms without adding the erf reference constant a
second time. P2A implements no electrostatic energy assembly or structure
factors and makes no QE zero-mode agreement claim.

Independent closed forms follow from the Gaussian Fourier integral and
Poisson's equation. For V=-Z*erf(r/rc)/r+A*exp(-r^2/s^2),

```text
Vloc(G) = [-4*pi*Z*exp(-G^2*rc^2/4)/G^2
           + A*pi^(3/2)*s^3*exp(-G^2*s^2/4)]/Omega
C = [pi*Z*rc^2 + A*pi^(3/2)*s^3]/Omega.
```

The tests evaluate these closed forms separately from the radial integrator,
including A=0, and compare C to direct quadrature of r^2*V+Z*r.
Local Fourier relative tolerance is 1e-8 over G=0.1..12 Bohr^-1;
zero-mode and final cutoff convergence absolute tolerance is 1e-11 Ha.
Parameters, measured errors and cutoff sensitivity are in `upf-support.md`.

## 17. P2A projector Fourier-Bessel normalization

Define B_l(q)=4*pi*integral r^2*beta(r)*j_l(qr)dr
=4*pi*integral r*u(r)*j_l(qr)dr. It has units Bohr^(3/2) under the
nominal UPF convention. Angular integration of a single beta(r)*Y_lm(rhat)
against a plane wave gives (-i)^l*B_l(q)*Y_lm(qhat). Later P2B must multiply
by Omega^(-1/2) and the atomic phase; those are not part of this routine.

For beta(r)=r^l*exp(-r^2/(2*s^2)), the Gaussian Cartesian transform is
(2*pi)^(3/2)*s^3*exp(-s^2*q^2/2). Applying the corresponding harmonic
polynomial in derivatives with respect to q yields, without radial
quadrature,

```text
B_l(q)=4*pi*sqrt(pi/2)*s^(2l+3)*q^l*exp(-s^2*q^2/2).
```

These analytic test functions have arbitrary amplitudes; dimensional units
can be restored by multiplying the r^l amplitude by Bohr^(-l-3/2).
The comparison uses l=0,1,2 and relative tolerance 1e-9.

Spherical-Bessel completeness reads
integral_0^infinity q^2*j_l(qr)*j_l(qr')dq=pi*delta(r-r')/(2*r^2).
Substituting the definition of B gives the exact Parseval statement

```text
(2/pi) integral q^2 |B_l(q)|^2 dq
    = (4*pi)^2 integral r^2 |beta(r)|^2 dr.
```

Equivalently, with integral |Y_lm|^2 dOmega=1, three-dimensional Fourier
Parseval is integral q^2|B_l|^2dq/(2*pi)^3 = integral r^2|beta|^2dr.
There is no extra angular 4*pi for a normalized single Y_lm. For the
Gaussian the right radial integral is
s^(2l+3)*Gamma(l+3/2)/2, evaluated from half-integer Gaussian moments.
Tests integrate the **computed** B_l on an independent 401-point linear
q mesh up to 12/s and compare at relative 1e-9.

## 18. P2A numerical evidence boundaries

Only synthetic strings and analytic functions have been verified. No UPF
file is committed or taken from a real pseudopotential library, and no QE
comparison is made. Parser fixtures have 601 points; separate analytic
quadrature meshes are refined to 12001 points for local transforms and
4001 for projectors. Absolute D symmetry tolerance is 1e-12 Ha. Exact zeros
beyond a declared projector cutoff and between different-l D channels are
structural requirements, not adjustable tolerances. All P2A test tolerances,
mesh sizes and measured errors are listed in `upf-support.md`. Unsampled
origin behavior, producer-dependent beta scaling, and exponentially small
large-G signals remain explicit numerical/input limitations.

## 19. P2B ionic and Kleinman–Bylander operators

`dft_physics::hamiltonian` uses the fixed full-complex basis order of section 6.
`Species` supplies a validated potential and Cartesian atomic positions in Bohr.
For each species s, S_s(G)=sum_I exp(-i G.tau_I). `ionic_local` returns
V_ion(G)=sum_s v_s(|G|) S_s(G), in Ha, **only on the supplied density sphere**;
all other Fourier components are zero. Supply the density basis of the same
cell and four-times wavefunction cutoff. The local transform already contains
1/Omega and the Coulomb tail of section 16. For G=0 use alpha_z(), including
its finite reference remainder, exactly once per atom. Nonzero modes use
at_g(), without an additional alpha-Z shift. Scatter on the density mapping
and apply to_real to obtain the grid potential; its imaginary roundoff is
retained for inspection by the caller before conversion to a real array.

Rayleigh expansion and translation of beta_i(r)Y_lm(rhat) give

```text
P_{I,i,m}(G) = <G|beta_{I,i,l,m}>
            = Omega^(-1/2) exp(-i G.tau_I) (-i)^l Y_lm(Ghat) B_i(|G|)
V_NL = sum_{I,i,j,m} |beta_{I,i,m}> D_ij <beta_{I,j,m}|
p_{I,j,m} = sum_G conj(P_{I,j,m}(G)) c(G)
(V_NL c)(G) = sum_{I,i,j,m} P_{I,i,m}(G) D_ij p_{I,j,m}.
```

B is the radial transform of section 17, with no additional 4*pi. The real
harmonics use section 15, ordered m=-l..l for each radial projector. At G=0
we evaluate the harmonic along the positive z axis: all l>0 radial transforms
vanish exactly, and Y_00=1/sqrt(4*pi), so this choice has no physical effect.
Nonzero D entries between unequal l are explicitly rejected (also rejected
by the parser); coupling is repeated only for matching l and m. The parser's
symmetry allowance of 1e-12 Ha is retained without symmetrizing D.

Construction evaluates each radial transform once per species, radial index
and G, and each harmonic per channel and G; atoms reuse these tables with
individual translation phases. Stored atom columns and small channel coupling
matrices are exposed read-only. Application first projects, then couples and
expands; no plane-wave-square matrix is stored. P is dimensionless under the
nominal UPF units and D and V_NL are in Ha. The public standalone nonlocal
apply allocates its small projection buffer; Hamiltonian application reuses it.

## 20. P2B local convolution and Hamiltonian lifecycle

`Hamiltonian::new` takes a real Array3<f64> V_eff(r) in Ha, the wavefunction
basis and cell, and its nonlocal operator. The caller forms V_eff from ionic,
Hartree and XC contributions; no SCF or energy assembly occurs here. The cell,
basis and nonlocal volume must describe the same physical cell. The constructor
checks the nonlocal G list against the wavefunction basis and derives the
minimum density grid from the supplied cell and cutoff, accepting any larger
valid dimensions supplied by the potential. All values must be finite.

```text
(V_eff c)(G) = gather(to_reciprocal(V_eff(r)*to_real(scatter(c))))
(H c)(G) = (|G|^2/2)c(G) + (V_eff c)(G) + (V_NL c)(G).
```

The wavefunction normalization factors sqrt(Omega) cancel between the two
transforms. Section 4's FFT scaling therefore gives v0*c for constant V_eff.
The reciprocal and real scratch arrays, FFT plan and projection scratch are
allocated at construction. `LinearOperator::apply` and `apply_block` reuse
these across calls; the block method borrows scratch once and streams columns.
Interior mutability is needed for the trait's shared receiver. The operator is
not Sync; a conflicting scratch borrow returns an application error. Backend
lane allocation remains possible as described in section 5. No per-call plan
creation, numerical I/O, unsafe code or input-dependent FFT shortcuts occur.

For Fourier support inside the density sphere the local matrix element is
exactly V_eff(G-G') up to roundoff. Indeed every difference G-G' is inside
that sphere, with axis bounds M_i, and every supported potential mode also
has bounds M_i. A wrapped equality q=G-G' modulo n_i would require a nonzero
multiple of n_i with magnitude at most 2*M_i, whereas n_i>=2*M_i+1. Thus none
is possible, even though the full product may have modes outside the grid.

Outside-sphere components require care: represented canonical modes outside
the sphere have no projected matrix elements, since G-G' is inside it. They
do not automatically cause aliasing. True unresolved high frequencies can
alias into the sphere on sampling; for example a cosine with Miller index
(n1+1,0,0) samples identically to (1,0,0), causing a spurious 0.5 Ha matrix
element for a unit-Ha cosine. `truncate_to_density_sphere` FFTs the input,
zeros represented coefficients outside the supplied density sphere, and
inverts. It returns real values, discarding imaginary roundoff (the retained
sphere contains both signs and excludes Nyquist modes). It cannot remove
already-aliased in-sphere coefficients. Resolve the original function on a
sufficiently fine grid before truncating if such frequencies matter.

## 21. P2B independent synthetic validation

The in-memory UPF fixture is SYNTHETIC Xx, Z=4, with two s and one p Gaussian
projector and a symmetric D with off-diagonal s coupling. It uses 801 logarithmic
radial points to 12 Bohr and an erf-screened Coulomb local potential. The
existing UPF test generator now exposes `fixture_with_points`; its original
`fixture()` still uses 601 points. No parser behavior or real UPF data changed.
All ordinary tests use the fcc diamond primitive cell, a=5.43 Angstrom,
two atoms, Ecut_wfc=6 Ry, 65 plane waves and a 12x12x12 FFT grid.

Reference matrices are built entrywise: local entries look up independently
indexed Fourier coefficients at the signed Miller difference (no FFT), kinetic
entries use |G|^2/2, and nonlocal entries explicitly sum atom and radial pairs.
Reference projector columns use separate Cartesian s,p harmonics and sine/cosine
phases and the trusted P2A radial transform; they do not read production tables.
The full Hamiltonian test includes ionic Fourier coefficients plus small real
band-limited terms. Grid samples are synthesized by a direct Fourier sum.
Fixed-seed complex vectors exercise Hermiticity. No tolerance is inferred from
the observed result or loosened after testing.

| Check | Bound |
| --- | --- |
| Ionic imaginary part / maximum real magnitude | <1e-12 |
| Ionic translation and direct coefficients | max absolute <1e-12 Ha |
| One-atom structure-factor modulus error | <1e-12 |
| Projector columns and translation | max absolute <1e-12, dimensionless |
| Nonlocal dense action | max absolute <1e-12 Ha |
| Nonlocal and H Hermiticity | relative <1e-12 |
| p-channel angular sum vs (3/(4*pi*Omega))*B_1(q)^2 | relative <1e-12 |
| Local convolution and constant potential action | max absolute <1e-12 Ha |
| H entries vs independent dense matrix | max absolute <1e-11 Ha |
| Free-electron dense eigenvalues | max absolute <1e-12 Ha |
| Full dense eigenvalues vs independent matrix | max absolute <1e-11 Ha |
| Lowest eight Davidson eigenvalues vs independent dense matrix | max absolute <1e-9 Ha |
| Davidson absolute residual 2-norm | <1e-10 Ha |
| Block vs separate applications | max absolute <1e-14 Ha |
| Consecutive applications | identical arrays |
| Filter retention/removal; analytic alias example | max absolute <1e-12 Ha |

The large-cutoff projector norm check is ignored in debug. The exact identity
for a periodized localized projector b(r)=beta_i(r)Y_lm(rhat) is

```text
sum_G |P(G)|^2 = integral_cell |sum_R b(r-R)|^2 dr
              = sum_R integral_R3 conj(b(r)) b(r-R) dr.
```

The R=0 term is integral_0^infinity |u_i(r)|^2 dr, dimensionless. Equality to
that isolated radial norm requires nonoverlapping periodic images (or the
large-cell limit), not merely increasing cutoff. The ignored test reports
finite sums at 16 and 64 Ry for all five channels, alongside the radial norm
and relative differences; there is no unjustified convergence tolerance.
Gaussian image overlaps and radial quadrature error remain at fixed cell.

Reproduction (offline toolchain alias on this host):

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export RUSTUP_HOME=/private/tmp/p0c-rustup
scripts/check.sh
cargo test --offline -p dft-physics -p dft-upf -- --nocapture
cargo test --offline --release -p dft-physics --test hamiltonian projector_norm_convergence -- --ignored --nocapture
```

These comparisons are mathematical synthetic tests, not QE references.

Measured P2B results on this host (the commands above exited 0): ionic
imaginary/real ratio 0; translation error 5.721958498152797e-17 Ha;
alpha-Z 0.046523648459421495 Ha per atom. Nonlocal action error
1.1487223069389168e-15 Ha, relative Hermiticity 3.3647837479443684e-16,
and p rotation relative error 6.249593393505841e-16. Local convolution
error 9.485749680535094e-16 Ha; full H entry error 8.883093345818638e-16 Ha;
relative H Hermiticity 6.864949379572692e-16; block error 0. Davidson took
88 iterations with maximum residual 8.749042190309645e-11 Ha and eigenvalue
error 9.769962616701378e-15 Ha. The complete fast P2B test binary took
5.88 seconds; the release norm test took 0.04 seconds, excluding compilation.

| Cutoff | Plane waves | Relative norm difference, s0 | s1 | each p channel |
| --- | --- | --- | --- | --- |
| 16 Ry | 283 | 1.3107218632658238e-3 | 2.2330470096246616e-5 | 9.001145971986984e-3 |
| 64 Ry | 2277 | 4.423816356522835e-11 | 2.3077474408052375e-5 | 9.001139925153145e-3 |

Radial norms are 0.15666426716443535, 0.44311346272637686 and
1.831615675600615, respectively (dimensionless). The broader projectors'
nonvanishing fixed-cell differences illustrate the periodic-image caveat;
no isolated-norm equality at fixed cell size is asserted.

## 22. P2B review: isolated normalization and independent image evidence

`projector_norm_isolation_and_images` closes the fixed-cell ambiguity in
section 21. It is ignored in debug and must be run with:

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export RUSTUP_HOME=/private/tmp/p0c-rustup
cargo test --offline --release -p dft-physics --test hamiltonian projector_norm_isolation_and_images -- --ignored --nocapture
scripts/check.sh
```

There is one atom at the origin in a cubic cell of side 40 Bohr. The actual
last radial mesh point is 12.000000000000012 Bohr. The nearest-image
**support-to-support gap** is 40-2*12 = 16 Bohr (15.999999999999975 in floating
point), exceeding the required 10 Bohr. The stored cutoff makes each
projector zero beyond this radius, so distinct periodic supports do not overlap.
The exact complete-basis identity is

```text
b(r) = beta_i(r) Y_lm(rhat)
P(G) = Omega^(-1/2) integral_R3 b(r) exp(-i G.r) dr
sum_G |P(G)|^2 = integral_cell |sum_R b(r-R)|^2 dr
              = integral_R3 |b(r)|^2 dr
                + sum_{R != 0} integral_R3 conj(b(r)) b(r-R) dr.
```

The isolated term is integral |u_i(r)|^2 dr, because u=r*beta and the
harmonic has unit angular norm. All these norms and overlaps are dimensionless.
The test evaluates the production projector columns, including every m, at
four cutoffs. The largest-cutoff relative error must be below **1e-8 for
each channel**. Errors must decrease while above that bound; strict
monotonicity is not imposed once quadrature and summation roundoff dominate.
No pre-existing tolerance changed.

The initial 32 Ry endpoint failed for s0 with error 5.334481956165149e-7.
This was investigated using independent Cartesian Gaussian Fourier integrals,
not accepted by relaxing the bound. For beta=exp(-a*r^2),

```text
|P_s(G)|^2 = pi^2/(4*Omega*a^3) * exp(-G^2/(2*a))
|P_p_axis(G)|^2 = |P_s(G)|^2 * 3*G_axis^2/(4*a^2).
```

Here a=1, 1/2, 1/3 Bohr^-2 for s0, s1, p respectively; p uses the second
formula with a=1/3. The independently summed s0 deficit at 32 Ry is
5.334481676388947e-7, and numerical-versus-analytic norm error is
2.799222200852258e-14 relative to the radial norm. Thus the failed endpoint
is a resolved finite reciprocal-cutoff tail, not an unexplained normalization
error. The test now includes 64 Ry and asserts numerical-versus-analytic
agreement below 1e-10 at every cutoff and for every channel.

Measured absolute relative differences from isolated radial norms:

| Cutoff (Ry) | Plane waves | s0, m=0 | s1, m=0 | p, m=-1 | p, m=0 | p, m=+1 |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 24405 | 4.624596137854697e-2 | 1.146562961354958e-3 | 2.205771232657305e-4 | 2.205771232665077e-4 | 2.205771232669518e-4 |
| 16 | 68999 | 1.148225308833961e-3 | 5.370953475303253e-7 | 3.689525795991244e-9 | 3.689526462125059e-9 | 3.689527017236571e-9 |
| 32 | 195269 | 5.334481956165149e-7 | 5.007105841059456e-14 | 3.169686735304822e-13 | 3.177458296477198e-13 | 3.181899188575699e-13 |
| 64 | 552629 | 2.709721336202620e-12 | 5.007105841059456e-14 | 3.169686735304822e-13 | 3.177458296477198e-13 | 3.181899188575699e-13 |

To predict the Si-cell image contribution independently, the test integrates
the Gaussian translated overlap analytically in real space. Completing the
square with t=r-R/2 gives

```text
I_s(R)/I_s(0) = exp(-a*R^2/2)
I_p_axis(R)/I_p_axis(0) = (1-a*R_axis^2) exp(-a*R^2/2).
```

The p factor follows from the Gaussian moment
< t_axis^2 > = 1/(4*a), since
r_axis*(r_axis-R_axis) = t_axis^2-R_axis^2/4. Its negative contributions
explain why the periodic p norm is smaller than the isolated norm. The test
sums these real-space overlaps over nonzero fcc translations with Miller
coordinates in [-4,4]^3 and verifies agreement with [-5,5]^3 to 1e-15
absolute. It does not use production radial transforms, Bessel functions,
harmonics or projector tables to compute the predicted image corrections.

The Gaussian formulas extend to infinity, whereas the fixture stops at
12 Bohr. This difference is negligible with an explicit bound: the largest
omitted radial norm fraction, for p, is 1.4597160698215574e-39, obtained
from integral_12^infinity r^4 exp(-2*r^2/3) dr using integration by parts
and erfc. Cauchy-Schwarz bounds the overlap perturbation per image relative
to the isolated norm by 2*sqrt(f)+f. Across the 728 nonzero translations
in [-4,4]^3 this is below 5.563e-17. Outside that box,
R^2 >= a_lattice^2*|m|^2/4, so the closest possible omitted image is farther
than 25.65 Bohr and the compact supports cannot overlap. The remaining
infinite-Gaussian image tail is bounded by

```text
sum_{k=5}^infinity (2*k+1)^3 * (1+3*a_lattice^2*k^2)
                 * exp(-a_lattice^2*k^2/24) < 2.446e-41.
```

Here k=max|m_i|; the overestimate counts an entire cube per shell and uses
the slowest exponent a=1/3. Successive terms are bounded by a geometric
ratio (6/5)^5*exp(-11*a_lattice^2/24). These errors are far smaller than
the required image-correction comparison tolerance.

At 64 Ry in the original Si cell, define the signed correction
(sum_G |P|^2 / radial_norm)-1. The independent prediction and measured result:

| Channel | Measured signed correction | Predicted signed image correction | Relative disagreement |
| --- | --- | --- | --- |
| s1 | +2.307747440810637e-5 | +2.307747440887537e-5 | 3.332257552956900e-11 |
| p, m=-1 | -9.001139925153168e-3 | -9.001139925151924e-3 | 1.381821349982378e-13 |
| p, m=0 | -9.001139925153168e-3 | -9.001139925151926e-3 | 1.379894123552834e-13 |
| p, m=+1 | -9.001139925153168e-3 | -9.001139925151922e-3 | 1.383748576411921e-13 |

Each of these four comparisons asserts relative disagreement below **1e-3**.
The corresponding absolute image terms are +1.022593959629611e-5 for s1
and approximately -1.648662898518281e-2 for each p channel. For completeness,
s0's measured correction is 4.423816868381891e-11 versus prediction
4.438073340755477e-11; its tiny signal is limited by finite-cutoff and
floating-point errors, so it is checked with the isolated 1e-8 absolute
relative-norm bound, not the s1/p image-signal-relative certificate.

The release test passed in 13.80 seconds excluding compilation. This evidence
separates normalization correctness in the isolated cell from the independently
predicted periodic-image contribution in Si; no production code changed.
