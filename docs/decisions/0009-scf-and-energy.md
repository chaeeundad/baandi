# ADR 0009: Charge-conserving SCF and independently audited energy

- Date: 2026-11-02
- Task: P5A-SCF-ENERGY
- Status: implemented; proposed for reviewer acceptance

## API and ownership

`dft_solvers::mixing::DensityMixer` mixes flat, real, uniform-grid number-density
samples in electrons/Bohr^3. `mix(input, output, mean)` takes the density used
in the Hamiltonian, the orbital density returned by that Hamiltonian, and the
prescribed mean N_e/Omega. The residual is **output - input**. It returns the
next input density. It does not mix potentials, orbitals, or charges per voxel.
`reset()` discards extrapolation history. Invalid dimensions, nonfinite values,
or nonpositive beta/mean/history return errors. Linear beta may exceed one;
unstable parameters are not silently capped.

`dft_engine::ScfProblem::new(cell, atoms, species_potentials, options)` prepares
ionic local/nonlocal operators and Ewald once. Atom labels resolve against
`SpeciesPotential { label, potential }`. The parser's arbitrary functional label
is not a physics guarantee: this driver accepts only the whitespace-separated
label `SLA PZ`. Valence charges determine the electron count; nonintegral counts
(discrepancy >1e-10 electron), odd counts, and insufficient bands are rejected.
Occupations are 2 for the lowest N_e/2 bands and zero for the rest. Bands default
to occupied + 4 and may not exceed the plane-wave dimension. Ecut is in Ha.
There is no real-pseudopotential or chemical accuracy claim.

`run()` starts from uniform density. `step(density, initial_vectors, tolerance)`
is an explicit single fixed-point map, also available to independent validators.
The result includes eigenvalues, occupations, orbitals, density, charge, every
energy term, Ewald subterms and tail bound, iteration records and phase timings.
There is no separate public Boolean that can be set to claim convergence:
`Ok(ScfResult)` is the success certificate. `NotConverged { history }` and
`Diverged { reason, history }` are distinct errors. Odd electrons and insufficient
bands have their own variants. Backend failures are reported separately.

Only the five permitted path edges from dft-engine to physics, solvers, basis,
geometry, and UPF are added. No external dependency is added. Because ownership
does not permit a direct engine-to-numerics edge, `dft_physics::data` re-exports
the array, complex, unit, FFT and matrix types already used in physics APIs,
and the adopted libm erf/erfc functions used by analytic fixtures. It implements
no new numerical backend. The only other physics additions are transactional
`Hamiltonian::set_local_potential` and read-only `Hamiltonian::nonlocal`.
The eigensolver implementation and existing tests are unchanged.

## Mixing and safeguards

Linear mixing is x_next = x + beta R. Pulay uses type-II Anderson with up to
m=6 differences (m+1 input/residual pairs). With columns DX and DR,

```text
gamma = argmin ||R - DR gamma||_2
x_next = x + beta R - (DX + beta DR) gamma.
```

The uniform quadrature factor would multiply the entire least-squares objective,
so Euclidean grid dot products suffice. Solve the small Gram system DR^T DR
through `dft_numerics::hermitian_eigh`, using its sequential faer backend.
Eigenvalues at or below 1e-12 times the largest are discarded: this is a spectral
pseudoinverse regularization, not a hand-written dense solver. Zero/singular
history remains finite; no ridge bias is added to well-resolved directions.
History is FIFO, deterministic, and automatically reset if shape or prescribed
mean changes, as well as on explicit reset. A positive representable history
capacity is required.

After either mixer, two uniform mean corrections enforce the requested charge
to floating-point summation accuracy. They alter only the zero Fourier mode,
without choosing a privileged spatial sample. No pointwise density clipping
or renormalization of orbitals is used.

Anderson is unconstrained and may propose a negative density. The driver then
uses x + theta*(proposal-x), with theta the minimum of 1 and 0.9*x_i/(x_i-proposal_i)
over negative proposed samples, and resets Pulay history. This is a convex
feasibility safeguard, preserving charge and Fourier support. Its fraction and
reset are recorded. Linear mixing is not damped this way: an invalid linear
proposal produces a divergence error. No synthetic-value recognition exists.

The original beta=0.3 linear trial was unstable for the small-gap synthetic
fixture. The acceptance run uses beta=0.05 and a 600-iteration budget; Pulay
uses beta=0.3, m=6. These are solver controls, not changed acceptance tolerances.
The second seed required the Pulay positivity safeguard; it is not skipped.

## Iteration and convergence

`V_eff = V_ion + V_H[n_in] + V_xc[n_in]`. The Hamiltonian keeps its FFT and
projector scratch across potential updates. Density and Hartree use a separate
reusable FFT plan. Davidson is warm-started with all previous eigenvectors and
uses kinetic diagonal preconditioning with a 0.1 Ha denominator floor. Default
subspace capacity is 48, Davidson iteration limit 400, SCF limit 200, seed 2026.
No dense eigensolver fallback is hidden in the SCF driver.

The default eigensolver schedule is
`max(0.1*eigen_tolerance, 1e-5 * 0.5^(iteration-1))` Ha. Both initial tolerance
and decay are configurable. Final convergence requires a solve requested at
least as tightly as eigen_tolerance, as well as its actual occupied residuals.
On **one and the same iteration**, require:

1. sqrt(Omega/Ngrid * sum (n_out-n_in)^2)/N_e <= 1e-7 by default.
2. Absolute change of the direct total energy <= 1e-8 Ha by default; the first
   iteration cannot meet this condition because it has no previous energy.
3. Maximum occupied eigen residual <= 1e-8 Ha by default.
4. Relative output charge error <= 1e-8; input charge is checked too.

The density residual has units Bohr^(-3/2) when electron counts are dimensionless.
The numeric tolerance is expressed in atomic units, not claimed dimensionless.

A small-gap fixed-point map can amplify its input residual. Therefore a candidate
success is additionally certified by rebuilding from n_out and re-diagonalizing
at 0.1*eigen_tolerance. This must change density by no more than density_tolerance
and occupied eigenvalues by no more than 1e-7 Ha. The original candidate orbitals,
density and audited energy are returned, with this independent certificate.
A failed certificate continues ordinary SCF iteration. Certificate Davidson work
and phase timings are included in the record. This extra solve was necessary:
the initial four-criterion candidate's external density residual was 1.78e-7.
No tolerance was loosened to accept it.

Nonfinite mixed densities/arithmetic/observables are divergence. Density below
-1e-14 Bohr^-3 is outside the existing P3A domain and is also divergence. Residual
growth beyond 20 times max(best residual, density_tolerance) for three consecutive
iterations is a distinct divergence route, configurable by factor and patience.
Exhausting the SCF budget returns NotConverged with all completed records.
Convergence is never inferred from energy change alone. Tests independently
invalidate each gate and verify that it prevents success.

Records hold direct/Harris energies, dE, density/eigen residuals, requested eigen
tolerance, Davidson iterations, input/output charge errors, minimum output
density, mixing safeguards, and timings for potential, diagonalization, density,
energy and mixing. All are in memory. The loop has no logging, JSON or file I/O.
Bitwise determinism concerns numerical outputs on the same build/backend; elapsed
timings are necessarily excluded. Cross-architecture bitwise identity is not
promised.

## Energy and independent electrostatics

Sections 23–25 of `docs/physics-conventions.md` give the exact formulas, zero-mode
cancellation, Gaussian correction and functional-difference bound. Production
energy evaluates every direct term using final orbitals and n_out. It also
reports the band identity with input potentials, the Harris-Foulkes functional
at n_in, and the separate local G=0 contribution. No Gaussian correction is
added to production SCF energy: that correction belongs only to the independent
comparison of different ionic conventions.

The key test assembles production Hartree, production P2A radial-transform local
energy including alpha-Z, and production point-ion Ewald. It compares to a
separate analytic reciprocal sum for neutral Gaussian ions and prescribed
periodic Gaussian electron densities. Its correction is derived from real-space
Gaussian pair interactions and self-energy, never from differences of computed
energies. It covers fcc/two Z=4 atoms and a skew triclinic/three Z=2,3,1 atoms;
rc=0.8,1.0,1.2 Bohr; two electron shapes; rigid translations; reciprocal cutoffs
72 and 90 Ha (48 cases). No nonlocal part is present in that test. Every case
must agree at relative 1e-9, independently of SCF convergence.

The parser backend cannot deserialize an empty explicit PP_DIJ text field.
The zero-projector synthetic fixture omits that optional tag, which the adopted
reader correctly defaults to an empty matrix. No parser change or dummy nonlocal
operator is needed.

## Verification and limits

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export RUSTUP_HOME=/private/tmp/p0c-rustup
export CARGO_NET_OFFLINE=true
scripts/check.sh
cargo test --offline -p dft-solvers -p dft-engine -- --nocapture
cargo test --offline --release -p dft-engine -p dft-solvers -- --ignored --nocapture
```

The debug gate includes a seven-plane-wave synthetic SCF, input validation,
convergence-gate tests, the additive Hamiltonian update contract, mixer linear
rate/finite termination/charge/reset/singular tests, and all unchanged existing
tests. The larger SCF sweep and electrostatic cross-check are ignored in debug
to honor the approximately ten-second per-test budget. The release command also
reruns the existing ignored Davidson tests.

SCF acceptance uses fcc a=5.43 Angstrom, two synthetic Z=4 ions, Gaussian local
potential and two s/one p projector with off-diagonal s coupling, 1201 radial
points, 10 Ry cutoff. It tests both mixers, same-seed bitwise outputs, another
seed, every-iteration charge conservation <1e-10, output minimum >=-1e-14
Bohr^-3, full energy audit, and external self-consistency. Controlled input
perturbations of amplitudes 0.04,0.02,0.01,0.005 test the functional-difference
power-law slope in [1.7,2.3], away from PZ branch discontinuities. The direct and
band forms agree within sum f*r plus 1e-12 Ha summation allowance. The direct and
Harris forms agree within the derived remainder bound plus the same allowance.

XC is nonlinear and is not bandlimited by the wavefunction density sphere.
The default alias-free density grid is not a promise of converged XC quadrature.
The translation/grid test resolves XC on 32^3 and 36^3 grids and asserts 1e-8 Ha
for both; the 8/10/12 Ry sweep uses 32^3. The default grid's measured XC quadrature
shift relative to 32^3 must be reported, rather than hidden or treated as a
G=0 discrepancy. No reference-energy values are asserted.

PZ81 retains its existing rounded branch discontinuity and low-density threshold.
Quadratic remainder asymptotics are asserted only away from those crossings.
General metallic occupations, smearing, real UPFs, CLI output and QE comparisons
are outside scope. Path dependencies require Cargo's regenerated engine lockfile
entry; the ticket prohibits retaining Cargo.lock edits, so the original is
restored after checks and the reviewer must reconcile it. Artifact/session JSON
paths are also outside ownership, and are left for the orchestrator with the
measured final report. Changes remain uncommitted.

## Measured acceptance evidence

On the existing aarch64 offline toolchain, the three commands above exited **0**.
`scripts/check.sh` reported PASS for build, test, fmt and clippy. The targeted
debug command passed 1 engine unit test, 2 engine integration tests and 8 solver
tests; 4 engine and 2 existing solver tests were ignored for release. The release
command passed all 6 ignored tests. The electrostatic binary took 18.55 s and
the SCF binary 22.50 s in release (including its concurrent tests); the small
SCF integration binary took 0.10 s in debug. These are observations, not timing
guarantees.

The Gaussian check's maximum relative discrepancy across all 48 cases was
**3.293860069749e-14**. The analytic constants and pair corrections are:

| Cell | rc (Bohr) | C=S-P (Ha) | P (Ha) |
| --- | --- | --- | --- |
| fcc, two Z=4 | 0.8 | 15.95769081399 | 4.020631e-7 |
| fcc, two Z=4 | 1.0 | 12.76602533096 | 1.276419e-4 |
| fcc, two Z=4 | 1.2 | 10.63538812443 | 3.072686e-3 |
| triclinic, Z=2,3,1 | 0.8 | 6.981381403713 | 1.085033e-4 |
| triclinic, Z=2,3,1 | 1.0 | 5.583187881788 | 2.004044e-3 |
| triclinic, Z=2,3,1 | 1.2 | 4.643742400825 | 1.058420e-2 |

At the largest width the pair correction is not negligible at the acceptance
accuracy; explicitly retaining it is essential. For example fcc rc=1, shape 0
has H=5.110931839806, local=-8.174945859835 (zero mode=0.7443783753507),
Ewald=-8.399471866510, and independent direct energy=1.302539444426 Ha.
The engine energy plus the analytically derived C agrees with this direct value.

The 10 Ry SCF with the default 15^3 grid gave:

| Quantity | Linear beta=0.05 | Pulay beta=0.3, m=6 |
| --- | --- | --- |
| Iterations | 252 | 15 |
| Total energy (Ha) | -4.6845611842778414 | -4.6845611842778592 |
| Density residual (Bohr^-3/2) | 9.672845558e-8 | 2.474588155e-8 |
| Max occupied eigen residual (Ha) | 8.929907996e-10 | 6.028917784e-10 |
| Minimum density (Bohr^-3) | 6.408853085e-3 | 6.408849945e-3 |

The energy difference is 1.776356839e-14 Ha. Same-seed energy, density, eigenvalues,
orbitals, iteration count and tested history observables are bitwise identical.
Seed 997 converged in 14 iterations, with one recorded positivity/history reset,
and differs by 8.881784197e-15 Ha. All printed SCF runs' maximum input/output
charge error over their histories was at most 2.620126338e-14 relative.

The Pulay term audit (Ha) is:

```text
T              4.6708948362521312
E_loc         -1.4485109028860925   (G=0: 0.74437837535074669)
E_NL           2.5014035784107116
E_H            0.31766640481016623
E_xc          -2.3265432343551438
E_ewald       -8.3994718665096322
sum           -4.6845611842778592  = -9.3691223685557183 Ry
band sum       3.3260985338528712
double count  -2.3976889779238766
band form     -4.6845611842778627
Harris        -4.6845611842778307
bound          2.885005432e-9
```

Ewald subterms are real=1.375261756529, reciprocal=0.1260128817421,
self=-4.950373252391, background=-4.950373252391 Ha, with a 4.829e-13 Ha
conservative truncation bound. External re-diagonalization changes occupied
eigenvalues by 1.315682785e-9 Ha and density by 1.597075436e-8 Bohr^-3/2.
The controlled residuals 1.955931422e-3, 9.797890577e-4, 4.903490424e-4 and
2.452876899e-4 Bohr^-3/2 give |direct-Harris| of 5.758434398e-4,
1.436030211e-4, 3.585025831e-5 and 8.955848996e-6 Ha. Successive exponents
are 2.008989, 2.004713 and 2.002415.

At 32^3, rigid translation changes energy by 1.060485033e-12 Ha;
32^3 -> 36^3 changes it by 1.208810829e-12 Ha. The **default 15^3 -> 32^3
change is 9.663798294923e-8 Ha**, larger than 1e-8. This is measured nonlinear
XC quadrature error on the minimal density grid, not an accepted relaxation of
the 32^3/36^3 invariance test. Density is bandlimited; its fractional powers in
LDA are not. The expected size for this fixture's default-grid error is about
1e-7 Ha; no universal bound is asserted for other systems. Oversampling resolves
it to the measured 1e-12 Ha scale.

| Cutoff (Ry), grid 32^3 | Total energy (Ha) | LUMO-HOMO (Ha) |
| --- | --- | --- |
| 8 | -4.6830233107047423 | 1.275672457215e-2 |
| 10 | -4.6845610876398762 | 7.947850095614e-3 |
| 12 | -4.6953191199626243 | 8.568851506765e-3 |

No external benchmark assertion is attached to this trend. The oversized linear
beta=100 case returns Diverged after one completed iteration; the one-iteration
budget returns NotConverged. Mixer tests measured the exact 0.8 linear contraction
rate (error 0.027487790694399994 after 12 steps), two-mode Anderson error
4.1522341120980855e-14 after three steps, and zero measured charge defect after
correcting deliberately wrong-charge output arrays. Singular/collinear histories
and history resets pass. Existing ignored Davidson results remain N=200 residual
7.788161106254981e-11 and N=300 preconditioned/unpreconditioned iterations 52/2286.

Additional unit-test tolerances are explicit: linear-rate error <2e-15 in
number-density fixture units, prescribed-charge defect <1e-14 electrons for
unit voxel volume, two-mode Anderson termination error <1e-12 in density units,
and collinear-history fixed-point error <1e-13. Reset/fresh-mixer agreement and
transactional Hamiltonian update/fresh-operator agreement use exact equality.
Ionic inverse-FFT imaginary components must be <=1e-12*max(1,max|Re V_ion|) Ha.
The existing density/Hartree/XC validation tolerances remain unchanged.

## P7P4 notes

O4c adds `ScfOptions::residual_eigen_schedule` (default true): first request
1e-2 Ha, subsequently `max(0.1*eigen_tolerance,
min(1e-2, 0.1*previous_density_residual))` Ha. The residual is the existing
normalized density L2 norm evaluated numerically in atomic units; the factor
0.1 is an empirical scheduling coefficient, not an equality of physical units.
False selects the original geometric schedule with its unchanged controls.
`free_electron_initial_guess` (default true) supplies lowest-kinetic coordinate
vectors for cold solves, with deterministic basis-index tie order. Existing
Davidson rank completion and Rayleigh–Ritz are reused. False restores the random
start. No solver parameter, convergence gate or fixed-point certificate changes.
The required final requested tolerance remains <= eigen_tolerance. Additive
in-memory application/per-k counters separate main and certificate work.
Measurements, rejected parameter settings and limitations are in the P7P4 entry
of `docs/performance/optimization-log.md`.

## P7PC notes

SCF uses a strict prefix of occupied bands plus the first empty band, capped
at the number of bands. Only occupied states enter the zero-temperature density
and energy; the lowest empty state is additionally needed by the unchanged
sampled-gap test. Higher empty states are subspace buffers. The engine's buffer
residual target is `max(requested_tolerance, 1e-5)` Ha. This measured local
choice is derived from study-note hypothesis H2, not copied from QE; the
suggested 100-times factor increased Gamma and si8 SCF iterations and was
rejected. There are three buffers in an eight-band/four-occupied calculation:
the fourth empty band count includes the strict lowest empty band.

A fixed-point certificate still solves the strict prefix at
`0.1 * eigen_tolerance` and retains the density and occupied-eigenvalue-change
tests. It records strict and buffer residual arrays separately per k point.
Buffer residuals are informational and are not certificate claims. The gap
criterion is unchanged. The default number of extra bands remains four.

Tolerance splitting alone did not fix si8: later SCF iterates stalled in the
strict lowest-empty member of the near-degenerate cluster. Engine working
capacity is therefore `max(davidson_subspace, 6 * bands)`, capped at dimension
by Davidson. The engine option now supplies a minimum capacity; the standalone
Davidson `max_subspace` contract is unchanged. This preserves the 48-dimensional
eight-band benchmark space and gives si8 120 dimensions. A 60-dimensional
trial converged but required 701 iterations versus 311 for 120 with the same
100-times policy; the final floor-only policy further reduced si8 to 283.
See the optimization log for timings, rejected trials and verification.

After the buffer policy passed energy and QE gates, intermediate SCF maps were
measured without the extra final Davidson true-residual block application.
This applies only while the scheduled tolerance is greater than
`eigen_tolerance`. Public `step`/`step_kpoints`, all final-accuracy maps and every
fixed-point certificate continue to recompute true residuals. Iteration-limit
acceptance also always recomputes them. The measured experiment preserves the
floor-policy density hashes and energies while saving 72/288/576/280 vector
applications on Gamma/k222s/k222/si8. No certificate work is skipped.
