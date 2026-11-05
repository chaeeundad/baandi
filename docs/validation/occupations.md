# P17B finite-temperature occupations

Implementation and evidence for MF04, against ADRs 0031–0034. This is a
synthetic numerical validation and an external-file smoke, not a QE-validated
metal calculation. The acceptance caps were taken from the MF04 row and ADR
0032 before measurements; no cap was relaxed.

## Numerical boundary

`dft_physics::occupations::FermiDirac::evaluate` accepts sorted spin/k/band
energies, explicit channel capacities (2 or 1,1), positive normalized k weights,
N, theta in Ha, and a total-number or fixed-moment constraint. Output storage
belongs to the caller. It allocates nothing and performs no I/O. Failed calls
must discard their output buffer. The engine clears any prior evaluation before
attempting a new root.

For each constrained sector, start at the minimum and maximum eigenvalues with
one theta of padding. Double that padding until the monotone number function
brackets the target, bounded by 256 attempts. Bisect with the overflow-safe
midpoint `(lo/2 + hi/2)`. Accept only a number residual no greater than
`1e-13 * sector_target`; a small bracket by itself never passes. Stalled or
unrepresentable brackets and exhausted iterations fail explicitly. Endpoint
sector targets have no finite chemical potential and reject. Equal-energy
states in a sector receive identical fractions. Root brackets and iteration
counts are exported; there is no universal mu-change convergence gate.

The probability and entropy use t=exp(-abs((epsilon-mu)/theta)) and
q=t/(1+t). The entropy of one capacity-one state is
log1p(t)+abs((epsilon-mu)/theta)*q, with the infinite-argument limit zero.
This preserves small hole entropy when the reported probability rounds to one.
Occupations include spin capacity but exclude k weight. Entropy includes both.
Independent tests use a tanh logistic and explicit binary entropy, not these
production functions.

## SCF and retained-band policy

The scalar driver reuses the existing Hamiltonian/eigensolver, density and
energy-audit path with an optional prepared occupation workspace. The spin map
solves every spin/k block before the shared root, then accumulates density and
energy in the original spin/k/band order. Fixed-policy summation order is
unchanged. Fixed-policy certificates, buffer-band tolerances and gap checks
remain on their existing branches. A finite-T spin result has no legacy
band-only occupation view or fixed-sector certificate; use `finite_occupations`
and `finite_certificate`. Its `sector` carries the setup axis/seed, not the
unconstrained final channel populations.

At preparation, the first implementation automatically fills the entire
explicit `band_cap` (CLI default 32), treating `scf.nbands` as a requested lower
bound. This conservative policy allocates occupation scratch once, does not
resize inside SCF, and records `prepare_explicit_cap` and the retained count.
A cap larger than the plane-wave basis rejects; decrease it or increase cutoff.
A cap without sufficient empty-state coverage rejects; increase it and repeat
the band study. Every retained finite-T eigenpair is strict. No tuning of the
legacy eigensolver buffer policy or default mixing is made.

Both maximum top-band fraction <=1e-8 and weighted top-band electrons
<=1e-10*max(1,N) are required. These are retained-basis diagnostics, **not an
infinite omitted-state bound**. The b,b+4,b+8 study remains necessary. Root and
occupation arrays are prepared once; the existing solver/map/mixer owned-result
allocation lifecycle remains. This is not a claim that the entire legacy SCF
implementation is allocation-free.

E is evaluated from the same output density, orbitals and fractional occupations
in both direct and band-sum audits. Positive TS is subtracted exactly once in
each free-energy audit. The CLI writes E, -TS and F ledger lines and labels F as
the variational quantity. All force/stress capabilities remain unavailable.

A candidate is followed by a fresh Hamiltonian and complete finite-T map at
0.1 times the eigen tolerance. Its certificate checks raw charge density and,
for spin, raw spin density; all retained eigen residuals; weighted absolute
occupation change; number and optional moment residuals; delta F; and both
edge tests. Density integrals are checked independently of root sums. Thresholds
cannot exceed 1e-7 Bohr^-3/2 for density or 1e-8 Ha for F/eigen residuals.
Occupation change is <=1e-8*max(1,N), number/moment <=1e-10*max(1,N).
The reported state remains the ordinary candidate, as in the existing driver.
The CLI adds a canonical policy digest to the certificate.

Unconstrained collinear FD uses one mu and allows channel populations to change.
Fixed M uses two roots at (N+M)/2 and (N-M)/2, including fractional targets.
An optional initial moment seed selects the initial density in the unconstrained
FD path; it is provenance, not model identity or a moment constraint. A symmetric
seed may retain a symmetric stationary branch. No global magnetic-minimum claim
is made. Existing fixed-sector seeds and behavior are unchanged.

## Reproduction and measured evidence

Run from the repository root (add the pinned rustup bin directory to PATH):

```sh
cargo test -p dft-physics --lib occupations -- --nocapture
cargo test -p dft-engine --lib finite_temperature -- --nocapture
cargo test -p dft-cli --lib finite_cli -- --nocapture
cargo test -p dft-engine --test collinear -- --nocapture
scripts/check.sh
python3 -I scripts/tests/test_readme_figures.py
DFT_NC_DIR=/path/to/nc/files cargo test --release -p dft-engine --lib real_al_finite -- --ignored --nocapture
```

Stored evidence: `tests/reference/occupations/synthetic.txt`, `al-smoke.txt`,
`summary.json`, `cli-unpolarized.json`, and the two legacy benchmark JSON files.
The synthetic tests use one soft-Coulomb NC atom, cubic 5 Bohr cell, Gamma,
cutoff 2.5 Ha (4.5 Ha for the N=2 low-T study), and no nonlocal projector.
They are not element predictions. The pure occupation oracle also covers
unequal k weights and two fixed-moment roots.

- Nine symmetric +/-a cases: a=0,0.01,1 Ha and theta=0.001,0.01,0.1 Ha;
  mu=0, complementary fractions, analytic occupations <=1e-12 and entropy
  <=1e-10 dimensionless. Degenerate capacity-weighted entropy is g ln(2).
- Synthetic N=1, theta=0.01 Ha: E=-0.46189330421395736 Ha,
  F=-0.47575624782515624 Ha, mu=-0.23127105844711310 Ha, 11 iterations.
- 12->16->20 bands: max delta F=1.5543122344752192e-14 Ha,
  max density distance=6.1760083645929322e-15 Bohr^-3/2.
- Stationary derivative at theta=0.05 Ha, h=0.004/0.002/0.001 Ha:
  abs(dF/dtheta+S/k_B)=5.6727159210012701e-6,
  1.3911603249372462e-6, 3.4608951726688986e-7 dimensionless, an O(h^2) trend.
- Gapped N=2 theta=0.1/0.01/0.001 Ha: delta E from fixed is
  6.6208358088497565e-2, 2.2204460492503131e-15,
  4.4408920985006262e-16 Ha. Low-T density distances are below 3e-16 Bohr^-3/2.
- Real Al external-file smoke: primitive fcc 3.8 Bohr off-diagonal columns,
  half-shifted 2x2x2 full mesh, PBE, 20 Ha cutoff, 20 bands, theta=0.01 Ha.
  E=-2.3089573748074912 Ha, mu=0.27591212552271216 Ha, 8 iterations.
  UPF SHA-256 b02eaa07c5d98f5eeae1bec4155854a7a94ea3635ca6b42e914f3f7ebd6912fb.
  This is not the R2 k/cutoff/band convergence study or a QE comparison.

The legacy gamma/k222s/k222/k444s benchmark used the initial gate's original
**debug** executable and the changed debug executable, one repetition, T=1,
against the same external Si UPF. All energies differ by 0 Ha, density SHA-256
values match, and all iteration counts remain 10. Timings are incidental,
not performance acceptance. Baseline source HEAD is
84f06a478549ce39151ba6095f5980bb652d8d26. P13B tests were not edited.

Earlier development failures were retained in the run transcript: formatting
and clippy duplicate-test-module findings; too many bands for a small basis;
and an insufficient high-temperature empty-state cap. The latter two were
fixed by increasing the test basis/retained coverage, not changing thresholds.
Inspection found an intermediate spin audit incorrectly using setup weights;
it was corrected and the direct/band identity regression assertion was added.
No failed development output is presented as accepted evidence.

## Boundaries and integration handoff

FD v1 supports the existing scalar functionals and collinear LDA-PW/PBE.
Unvalidated combinations (legacy band/DOS/volumetric requests, XC floor,
runtime control, final spectral guard, and v1 checkpoints) reject explicitly.
A fixed v1 checkpoint cannot be resumed under FD; the test creates a real fixed
checkpoint and checks this rejection. No v2 persistence is implemented.
Results use the proposed `state_contracts` extension with canonical physical
model/basis/state records, explicit spin/k/band layout and unavailable derivatives.
A full external JSON Schema validator was unavailable in the offline environment;
structural fields and numerical identities are tested in Rust.

Required edits outside this ticket's OWNED_FILES: extend
`scripts/make_readme_figures.py` for the new occupations evidence and regenerate
README.md (and a figure if selected). The existing generated README gate passes,
but does not yet include this new result family. The orchestrator owns dashboard
publication and the main checkout's run registry. No frozen schema or git
metadata was edited. Changes are uncommitted for review.
