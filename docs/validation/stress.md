# NC stress and pressure (P22A / MF14)

The engine computes the complete post-convergence frozen-Miller strain
derivative for its supported NC models. `stress(&StressState)` adds a separate
isotropic cutoff-derivative Pulay estimate. `frozen_stress` exposes the analytic
pass alone for independent validation; its `pulay=None` explicitly means no
cutoff assessment. Neither API asserts cutoff convergence.

## Sign, units and thermodynamic quantity

With lattice columns A, A(eta)=(I+eta)A, fixed fractional ions and k points,
sigma_ab=(1/Omega) dQ/d eta_ab. Q=E for fixed insulating occupations and
Q=F=E-TS for Fermi-Dirac occupations at fixed theta,N (and fixed M when
requested). The Cartesian 3x3 tensor is in Ha/Bohr^3, tensile positive.
P=-trace(sigma)/3 is positive for compression. A cell stretched from stable
equilibrium has positive tensile stress; length alone does not determine sign.
For symmetric shear eta_ab=eta_ba=h, dQ/dh=2 Omega sigma_ab. There is no
implicit engineering-shear factor in the tensor or ASE stress vector.

Explicit conversions: 1 Ha/Bohr^3 = 2 Ry/Bohr^3, approximately 29421.0157 GPa
= 294210.157 kbar. One kbar is 0.1 GPa. Opposite-sign conventions require an
explicit multiplication by -1 before comparison. No external-code comparison
or silent sign/unit conversion is included in this evidence.

## Analytic contributions

All derivatives are evaluated after convergence; no SCF iteration requests
stress. The engine returns an eight-entry derivative ledger, per-term tensors,
pressure, symmetric tensor, raw antisymmetric part/norm, timing, and retained
band counts. Inactive core/entropy entries are distinguished from unavailable
active derivatives. A missing active derivative must reject the total, not
supply a zero. The supported closed NC ledger is checked through the force
preconditions and exhaustive energy audit. New physical energy fields require
a derivative-coverage review.

- Kinetic: -sum_kbG w_k f_kb |c_kbG|^2 q_a q_b/Omega, q=G+k.
  Fractional k, Miller indices, and normalized coefficients are fixed.
- Hartree: differentiate the full FFT Coulomb metric and volume. Each nonzero
  Fourier component contributes 2*pi*|n_G|^2/G^2 times
  (2 G_a G_b/G^2-delta_ab). G=0 remains excluded exactly as in energy.
- Local: differentiate the Coulomb-subtracted radial quadrature analytically
  using j0'=-j1, its analytic Coulomb tail, and 1/Omega. The alpha-Z zero
  mode contributes its nonzero isotropic volume derivative. Integer ionic
  structure-factor phases are invariant at fixed fractional positions.
- Nonlocal: differentiate radial Bessel transforms and real harmonic angular
  factors, plus Omega^-1/2. Polynomial solid-harmonic recurrence derivatives
  avoid polar-axis singularities for l=0..3, the supported UPF range. Contract
  both sides of the KB expectation with the full Hermitian D matrix, including
  multiple radial projectors with the same l. Fractional ionic/k phases do not
  differentiate. No finite-difference projector derivative is used.
- XC: LDA-PZ, LDA-PW and PBE include (E_xc-integral n_val V_xc)/Omega on the
  diagonal. PBE additionally differentiates the spectral gradient metric using
  library sigma derivatives. Collinear channels retain their separate gradients
  and the single cross-sigma derivative; no hidden factor of two. The same
  screening and Nyquist conventions as the energy are used. A configured
  unpolarized density floor uses the existing masked potential chain rule and
  the gradient of the floored combined density.
- NLCC: contract the mean spin XC potential against the analytic core Fourier
  derivative, including its radial, G=0 and 1/Omega terms. This supplies the
  core density response excluded from the valence integral above, not a second
  core energy. Core is split equally between collinear spin channels.
- Ewald: real-distance, reciprocal-metric and volume derivatives, with
  -E_background*delta_ab/Omega. The fixed-eta self derivative is zero. Eta
  response cancels between the four infinite sums; tests vary eta to measure
  truncation consistency. The energy tail bound is not a stress error bound.
- Entropy: zero explicit strain derivative at fixed fractions and theta. At a
  stationary FD solution the occupation response in dE cancels theta*dS up to
  chemical-potential constraints, sum w df=0 in each constrained sector.
  The oracle differentiates converged F, not E, and uses a fixture with resolved
  fractional occupations and TS>1e-6 Ha.

## Frozen and changing basis protocols

`PlaneWaveBasis::in_cell` is additive: it preserves the exact ordered Miller
list and fractional k, recomputes Cartesian G/G+k, and retains the reference
cutoff as metadata. It does not re-sort or reapply a cutoff sphere. Legacy
constructors and their arithmetic are unchanged.

`ScfOptions.strain_reference: Option<Box<Cell>>` selects all orbital and
ionic/core density reciprocal lists in the reference cell; `grid_override`
is mandatory for this protocol and must hold the reference FFT dimensions.
This opt-in path is for strain oracles, including spin and FD. Default None
uses production cutoff construction. Identity checks include the reference
cell. Existing checkpoint v1 continuation export rejects a frozen reference
instead of silently discarding it. Strain scans use the full input quadrature;
this work does not establish a broken-symmetry strain scan from a reduced mesh.

The reported stress is always `basis_protocol=frozen_miller_set`. The separate
estimate uses the standard isotropic approximation

    sigma_fixed_cutoff - sigma_frozen ~= (2 Ecut / (3 Omega)) dQ/dEcut * I.

Three independent strict SCFs at Ecut, 1.05 Ecut, 1.10 Ecut rebuild their bases.
They share an FFT grid at least as large as both the original grid and the
largest-cutoff requirement; the center is recomputed on that grid too. The
reported correction uses the 5% slope. Its uncertainty diagnostic is the
absolute difference between the 5%/10% estimates plus the configured SCF energy
threshold contribution. It is a two-step spread, **not a rigorous error bound**;
cutoff shell discreteness and unresolved anisotropic errors require additional
cutoff/k/grid scans. If no orbital shell is crossed, `shells_resolved=false`
and the estimate warns even if both measured slopes are zero.

Warn when abs(correction)+uncertainty exceeds 1e-5 Ha/Bohr^3 or shell sampling
is unresolved. `cutoff_converged=false` is always explicit. The estimate is
reported separately and never added silently to the frozen tensor or pressure.
Its three SCFs have their own elapsed time and are excluded from analytic-pass
timing. A scan failure propagates an error rather than fabricating an estimate.

## Preconditions, symmetry and interfaces

The independent setup-model identity (cell, atom order/positions, potential
contents, functional/floor, reciprocal model, occupations and symmetry), strict
profile/tolerance caps, and fresh passed stationary certificate are the same
as forces. Standard profile, unconverged or mismatched states reject. Validation
currently reuses the force pass; its work is included in analytic elapsed time
and reported separately as `precondition_force_band_passes`. This conservative
implementation has overhead; no performance optimization is claimed.

The raw antisymmetric Frobenius norm must be below 1e-10 Ha/Bohr^3. It is
reported before averaging, so group averaging cannot conceal a failed check.
For symmetry auto, the exact spatial operations retained by SCF average each
term as R_cart sigma R_cart^T with R_cart=A R_frac A^-1. Returned terms and
symmetric total are symmetric. Raw physical stresses are not optimizer forces.

CLI `[output] stress=true` writes a `stress` block, including the ledger,
per-term tensor, pressure, preconditions, Pulay estimate, warnings and timings.
Fixed unpolarized runs print a table; spin/FD also write `stress.txt`.
The extension-version-2 proposal schema has available/unavailable branches;
unavailable values must be absent. Version 1 cannot advertise stress.

The existing Gamma/unpolarized ASE adapter supports `atoms.get_stress()` in
ASE order (xx,yy,zz,yz,xz,xy), tensile-positive eV/Angstrom^3. Conversion uses
BAANDI's unit constants. ASE 3.26.0's own `calculate_numerical_stress` is the
independent sign/shear oracle. Native `Session.stress()` returns the tensor in
Ha/Bohr^3 and emits a Python warning for a flagged Pulay assessment. Spin and FD
are engine/CLI capabilities; this ticket does not extend the existing Python
session to those models or alias energy as ASE free_energy.

## Reproduction and evidence

All synthetic tests contain only generated fixtures, not real-element claims.
Strain steps are 0.002, 0.001, 0.0005 in three diagonal and three symmetric shear
directions. Fixed-Miller SCFs retain fractional atoms/k and FFT dimensions.
Tests require the synthetic 1e-7 Ha/Bohr^3 cap and second-order error reduction
or a resolved small numerical plateau, without relaxing tolerances. Cubic
isotropy is checked at 1e-9 Ha/Bohr^3; all supported totals enforce the
antisymmetric norm gate.

```sh
cargo test --release -p dft-physics --test stress_terms -- --include-ignored --nocapture
cargo test --release -p dft-engine --test stress synthetic_local_nonlocal -- --ignored --nocapture
cargo test --release -p dft-engine --test stress synthetic_models -- --ignored --nocapture
cargo test --release -p dft-engine --test stress symmetry_preconditions -- --ignored --nocapture
DFT_NC_DIR=/path/to/NC cargo test --release -p dft-engine --test stress real_files -- --ignored --nocapture
cargo test -p dft-cli --test stress -- --nocapture
python3 -I tests/reference/stress/measure_cost.py --binary target/release/dft-cli --upf /path/Si.pz-vbc.UPF --output tests/reference/stress/cost.json
```

The extended release sweeps are ignored in debug. `DFT_STRESS_CASE` selects
skew, slab, nlcc, pbe, spin, spin-pbe or fd for diagnosis, and
`DFT_STRESS_MATERIAL` selects Si or MgO. Slab-like tests remain 3D periodic;
no 2D Coulomb cutoff or isolated-boundary stress is claimed. The initial NLCC
1 Ha fixture failed the existing positivity rule; it was rerun at 2 Ha with a
48^3 grid. The initial 0.01 Ha FD fixture had negligible entropy; the final
0.05 Ha/48-band fixture resolves fractions and passes the fresh certificate.
Neither failure was hidden by clipping density or weakening a tolerance.

External readable PseudoDojo files were used for Si at a=10.2612 and 10.6 Bohr,
and MgO at a=7.96 Bohr, with 20/40 Ha cutoffs respectively. The test repeats at
1.5 times cutoff, twice the grid dimensions, and a full 2x2x2 mesh. All three
are **inconclusive as cutoff/k/grid-converged material results**: measured
maximum uncertainty exceeds 1e-5/3 Ha/Bohr^3. Their frozen-set strain and
energy-versus-volume comparisons are retained as finite-basis observations,
not certified material acceptance. No UPF is committed.

Stored logs, summaries, command exit codes, UPF hashes and Mac cost observations
live in `tests/reference/stress/`. `measure_cost.py` compares stress on/off
energy and density bits, SCF iterations and H/operator work before reporting
analytic-pass share on si8 and k444s. Timing samples are Mac observations only,
not dev08 acceptance or a speedup claim. README regeneration remains assigned
to the orchestrator's later ticket.


### Measured results

The ten synthetic total-stress cases have maximum final-step error
1.3006765445505025e-9 Ha/Bohr^3, below the frozen 1e-7 cap. The largest recorded
raw antisymmetric norm is 8.571370772892337e-20 Ha/Bohr^3. The polarized
PBE+NLCC channel-swap pair differs by 6.358195220323992e-15 Ha/Bohr^3.
The FD metal has theta=0.05 Ha, TS=0.05634887876816284 Ha,
E=-0.35539139181363266 Ha and F=-0.4117402705817955 Ha; the free-energy
strain error is at most 2.9410121018223767e-10 Ha/Bohr^3. The fixed-moment
spin-FD case gives at most 2.96504169741002e-10 Ha/Bohr^3.

| Real finite-basis observation | Maximum final strain error, Ha/Bohr^3 | Cutoff/k/grid uncertainty, Ha/Bohr^3 | Qualification |
| --- | ---: | ---: | --- |
| Si, a=10.2612 Bohr | 1.0647756456920454e-9 | 1.1386220514642813e-3 | inconclusive |
| Si, a=10.6 Bohr | 7.454107108038205e-10 | 9.236788447450177e-4 | inconclusive |
| MgO, a=7.96 Bohr | 5.933744163123117e-9 | 2.219785600683727e-3 | inconclusive |

ASE's own strain-energy oracle differs by at most
3.6657227253655833e-10 eV/Angstrom^3. Eight installed-wheel Python tests pass
with zero skips, including native/CLI equality and existing relaxation tests.
The explicit CLI stress schema test additionally covers fixed, fixed-collinear,
FD and collinear-FD output. Its spin/FD group takes over 10 s in debug, so run
`cargo test --release -p dft-cli --test stress -- --include-ignored --nocapture`.

| Mac T=1 observation (one pair) | SCF, s | Analytic pass, s | Pass / (SCF+pass) | Separate Pulay scans, s | Retained bands / force-validation bands |
| --- | ---: | ---: | ---: | ---: | ---: |
| si8 | 2.552709833 | 3.584649583 | 58.407% | 8.745434875 | 20 / 20 |
| k444s | 2.460176 | 1.981509667 | 44.612% | 8.533295792 | 256 / 256 |

The pass counters count retained band records per derivative sweep; they are
not a count of all floating-point operations. There are three additional SCFs
per Pulay assessment. Stress on/off has identical energy, density SHA-256,
SCF iteration count, Davidson iteration count and operator applications on both
cost cases. Five legacy cases also pass `bench_scf.py --require-identical`
against an unmodified same-HEAD archive (zero energy difference, identical
density hashes). This final identity run uses the default thread budget; no
all-thread performance acceptance is claimed. The historical repository baseline
still differs in density and is retained as an initial failed comparison, not
used to assert bit identity for the new feature.

The initial isolated-term evidence and same-HEAD baseline records predate the
full feature; `validation.json`, `scf-final.json`, `cost.json` and `gates.json`
are the current evidence entry points. A performance-owner ledger entry and
README/dashboard publication are orchestrator-owned follow-up artifacts;
this ticket makes no optimization acceptance or publication claim.
