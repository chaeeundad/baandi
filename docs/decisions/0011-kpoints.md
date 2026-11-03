# ADR 0011: Weighted Bloch bases and insulating k-point SCF

- Date: 2026-11-03
- Task: P7A-KPOINTS
- Status: implemented; proposed for reviewer acceptance

## Convention and API

`psi_nk(r) = Omega^(-1/2) sum_G c_nk(G) exp(+i (k+G).r)`, with
`sum_G |c_nk(G)|^2 = 1`. Library callers supply dimensionless fractional
reciprocal coordinates; internally `k = B*k_fractional` is Cartesian Bohr^-1.
`PlaneWaveBasis::at_k(cell, cutoff, fractional)` canonicalizes fractional k
modulo integers to `[-1/2,1/2)`, then enumerates `|k+G|^2/2 <= Ecut` using the
existing relative cutoff allowance 1e-12. The dual-lattice bounding box has an
extra integer for the fractional shift. Sorting retains the original G-norm,
then Miller-index order to preserve Gamma rounding and solver initialization.
`vectors()` still means G; `momenta()` means k+G with corresponding G labels.
`k_fractional()` and `k_cartesian()` expose the canonical coordinates.
Reciprocal translations use the same canonical basis and coefficients, rather
than returning coefficients in an unwrapped gauge. Coordinates must be finite;
floating-point fractional information lost in very large input coordinates
cannot be reconstructed.

`PlaneWaveBasis::new` calls `at_k` at zero. Kinetic energy and Hamiltonian
kinetic diagonals consume `momenta()`. NonlocalOperator uses the same code at
Gamma and nonzero k; no Gamma-specific physics implementation exists.
The density basis remains the unshifted sphere at 4*Ecut. Its grid rule and
optional oversampling are unchanged and independent of the quadrature.
In a skewed cell a shifted sphere's absolute G can lie outside canonical FFT
frequencies even when all G-G' differences fit the density sphere. Such G
indices are mapped modulo the fixed grid dimensions, with collision rejection;
Gamma's existing canonical-index checks remain unchanged. The dedicated skew
fixture has q=(0.01,0,0), G=(0,10,0) Bohr^-1 on a one-point density grid.
The FFT transforms the periodic part u_nk. The Bloch phase cancels in
`|psi_nk|^2`; local Hamiltonian application is the usual periodic FFT product.

`dft_geometry::kpoints` is an additive module:

- `KPoint { fractional, weight }` supplies a quadrature element.
- `monkhorst_pack([n1,n2,n3], [s1,s2,s3])` uses
  `k_i=(m_i+s_i)/n_i`, `m_i=0..n_i-1`, with each shift exactly 0 or 1/2 in
  grid-step units, wrapped to `[-1/2,1/2)`. This explicitly names the
  Gamma-centered unshifted convention, not the original centered MP formula.
- `reduce_time_reversal` merges k with -k modulo reciprocal integers, retaining
  the first representative and summing weights. It also coalesces duplicates;
  it is idempotent on its output. Equivalence allowance is 1e-12 fractional.
  No atomic or space-group symmetry is consulted. For uniform input, weights
  equal multiplicity divided by the full count.
- `validate` rejects empty sets, nonfinite coordinates, nonpositive/nonfinite
  weights, and absolute pairwise weight-sum error above 1e-14. It never renormalizes.

Independent integer enumeration establishes unshifted 2^3 -> **8**, half-shifted
2^3 -> **4**, unshifted 4^3 -> **36**, half-shifted 4^3 -> **32**. All eight
unshifted 2^3 points are self-inverse. The ticket's unshifted 2^3 -> 4 example
therefore applies to the half-shifted convention instead. Odd 3^3 gives 14
representatives for either shift. All measured weight-sum errors for grids 2, 3 and 4 are zero. A 12^3
grid also checks summation stability with independent chunked summation.

## Nonlocal convention

For q=k+G, projector columns are

```text
P_I,i,l,m(G) = Omega^(-1/2) (-i)^l Y_real_lm(q_hat)
               * beta_transform_i(|q|) * exp(-i q.tau_I)
V_NL = sum_I P_I D_I P_I^H.
```

The radial transform includes the established 4*pi radial integral and the
UPF conversion already in dft-upf. D couples only equal l,m channels.
The special direction at q=0 is chosen only when its squared norm is zero,
not when the Miller index is zero (G=0 generally has q=k). This matters for
nonzero-k p projectors. Atom positions are Cartesian Bohr. No UPF parser,
harmonic, radial-transform or unit-conversion implementation changes.

## SCF and occupations

`ScfProblem::with_kpoints(cell, atoms, potentials, options, points, policy)` is
the new entry point. `OccupationPolicy::Insulator { minimum_gap_hartree }`
requires a positive finite threshold and at least one empty band at every k.
`PartialOccupations` is explicitly rejected, as are odd/nonintegral electron
counts. Every k has f_n=2 on the lowest N_el/2 bands and zero otherwise.
Intermediate fixed-point maps use these integer occupations; accepting a
converged result also requires the sampled indirect gap
`min_k eps_LUMO(k) - max_k eps_HOMO(k)` to exceed both the supplied threshold
and twice the eigen tolerance. Both the candidate and its fresh fixed-point
certificate are checked. A failing gap returns `ScfError::NonInsulating`.
A finite quadrature certifies only its sampled gap, not unsampled regions.

The legacy `new` constructor preserves the historical Gamma occupation contract
(including occupied-only band counts) and numerical results. It has no added
gap certificate. This compatibility choice is explicit; the new k API cannot
silently opt out. `run_outcome` may retain an unconverged state as before; its
reported gap is then a diagnostic, not a certificate.

Each k retains its own basis, FFT mapping, Hamiltonian and previous Davidson
vectors. Davidson uses that k's diagonal kinetic preconditioner (0.1 Ha floor).
No solver trait or eigensolver implementation changes. `step_kpoints` accepts
one initial matrix per point. Legacy `step` delegates for one point and rejects
multi-point calls. All points share the same effective potential and density.

```text
n(r) = sum_k w_k sum_n f_n |psi_nk(r)|^2
T = sum_k w_k sum_n f_n sum_G |c_nk(G)|^2 |k+G|^2/2
E_NL = sum_k w_k sum_n f_n <c_nk|V_NL(k)|c_nk>
B = sum_k w_k sum_n f_n eps_nk
```

Energy audits for each orbital block use the **same total weighted density**,
then are weight-averaged. Thus kinetic/nonlocal/band terms are k-weighted and
Hartree, XC, ionic-local and Ewald terms occur once (up to weight-sum roundoff).
In particular this does not average Hartree/XC of individual k densities.
Direct, band-form and Harris-Foulkes audits retain the ADR 0009 definitions;
the residual bound uses weighted eigen residuals. Repeated common-density
functional evaluation is accepted here; performance work is out of scope.

Existing density, charge, energy-change and eigen-residual convergence thresholds
are unchanged. The maximum occupied residual and the fixed-point eigenvalue
change cover **all** k points. Davidson iteration counts sum the per-k work.
`ScfStep::k_points` and `ScfResult::k_points` report every eigenpair block with
its coordinates and weight; the old solution/eigenvalues/orbitals fields are
first-k compatibility views. `ScfResult::gap_hartree` reports the sampled gap,
or None for the legacy constructor. Orbitals use the canonical basis ordering.

## Validation and reproduction

No dependencies, CLI/TOML, real UPF files, or QE results are added. Synthetic
pseudopotentials are generated in memory. The 2-species SCF motif has generic
positions and different Gaussian widths, so inversion cannot exchange its
atoms. The full/reduced comparison tests actual time reversal without space-
group reduction. Tests use tight SCF controls (density 1e-12, energy 1e-13 Ha,
eigen residual 1e-12 Ha) for the 1e-10 comparison, not relaxed acceptance bounds.

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export CARGO_NET_OFFLINE=true
scripts/check.sh
cargo test --offline -p dft-geometry kpoints -- --nocapture
cargo test --offline -p dft-engine --test kpoints -- --nocapture
cargo test --offline --release -p dft-engine --test kpoints -- --ignored --nocapture
cargo test --offline --release -p dft-engine -p dft-physics -- --ignored --nocapture
```

Longer SCF regression and invariance tests are ignored in the debug gate and
carry their release command. The small Gamma reference was measured before
editing with `cargo test --offline -p dft-engine --test scf
functional_labels_preserve_energy_and_reject_other_functionals -- --nocapture`:
-1.8482289692971974 Ha. The second Gamma reference uses the existing synthetic
nonlocal acceptance fixture and the measurements recorded in ADR 0009.
Neither number is a QE reference.

Measured evidence is recorded below after final verification. Session/gate JSON
files are outside this ticket's owned paths; the reviewer must reconcile the
final report into those artifacts. Changes remain uncommitted as required.

Final numerical evidence (Hartree unless otherwise specified):

| Check | Measured result | Required bound |
| --- | --- | --- |
| Free electron, skew cell, nonzero k, 7 plane waves | max eigenvalue error 2.2204460492503131e-16 Ha | 1e-12 Ha |
| Nonlocal H(k) Hermiticity, release | 8.8817841970012523e-16 Ha | 1e-12 Ha |
| Nonlocal direct projector dense construction | max error 9.9340186018569745e-16 Ha | 1e-10 Ha |
| Reciprocal translation k -> k+(1,-2,3) | energy/eigenvalue/density differences all 0 | 1e-10 Ha / electron Bohr^-3 |
| Full/reduced half-shifted 2^3 SCF | 2.2050857162549464 / 2.2050857162549451 Ha | energy difference 1e-10 Ha |
| Full/reduced energy difference | 1.3322676295501878e-15 Ha | 1e-10 Ha |
| Full/reduced maximum density difference | 1.5681900222830336e-15 electron Bohr^-3 | 1e-10 electron Bohr^-3 |
| Unequal quadrature weights, full/reduced energy difference | 1.7763568394002505e-15 Ha | 1e-10 Ha |
| Unequal weights, maximum density difference | 5.6898930012039273e-16 electron Bohr^-3 | 1e-10 electron Bohr^-3 |
| Small Gamma, legacy and k API | both -1.8482289692971974 Ha, pre-P7 difference 0 | 1e-10 Ha |
| Nonlocal Gamma, k API | -4.6845611842778592 Ha, ADR 0009 difference 0 | 1e-10 Ha |
| Metallic 3x1x1 sample | rejected, indirect gap -0.06694262433472065 Ha | no insulating success |
| 12^3 quadrature, independent chunk sum | weight error 1.4432899320127035e-15 | 1e-14 |
| Skewed shifted basis on 1x1x1 FFT | kinetic error 0 Ha, scatter/gather exact | 1e-12 Ha |

The uniform full/reduced run's sampled indirect gap is 0.07127723463067315 Ha;
the unequal-weight run's is 0.0716851885221701 Ha. Both charge and independent
band-form energy assertions pass. Five k-point release tests pass in 0.70 s
excluding compilation. All eight ignored engine/physics acceptance tests pass,
including the existing electrostatic, SCF and projector-normalization suites.
The whole-workspace gate passes build, tests, format and clippy. An intermediate
clippy failure on redundant decimal digits in the fixed reference literal was
corrected without changing its represented f64 value or any tolerance.
