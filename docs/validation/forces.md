# P20A force implementation and validation

This is fixed-cell, norm-conserving Hellmann–Feynman differentiation. Forces
are Cartesian, input atom order, Ha/Bohr. `dft_engine::forces(&ForceState)`
returns Ewald, local, nonlocal and core arrays, their complete sum, net force,
maximum absolute component, a ledger with active flags, energies and audit roles, elapsed seconds and the number of band visits.
The engine independently retains a setup-time model description in each result.
Before force work, it compares the cell, ordered atom positions/species, cutoff,
functional, actual grid, k points/weights, density floor, symmetry operations and
occupation policy against that record. The occupation record includes fixed
occupations or FD theta/constraint/band count and the collinear sector/axis.
Mismatches fail `observables.forces.model: ...; use its original model or rerun
SCF with the requested model`. The record is private and is not derived from
mutable result arrays at query time.

Original UPF file SHA-256 values are not passed to the engine. Instead, species
identity compares all parsed physics content exactly: element, functional,
valence, radial mesh/weights, local potential, core density, every projector's
l/cutoff/radial samples and the full D matrix. Independently parsed identical
potentials match; changed contents fail without hashing or pointer identity.
File-byte/comment provenance remains the adapter's SHA-256 responsibility.
The retained PP arrays are shared across result clones; comparisons run once
per force query, outside SCF numerical kernels. This is a model identity check,
not a cryptographic attestation of manually modified orbitals or certificates.

## Formula and ownership

The engine's closed NC ledger names kinetic, local, nonlocal, Hartree, XC,
core, Ewald and entropy owners. `require_complete` collects every unavailable
active derivative before returning an error; there is no partial ForceResult.
New energy owners must extend this ledger before exposing their model. Stress
is unavailable throughout. Hubbard and external coordinate-dependent fields
are not represented by this engine model and are not supported force models.

With forward-normalized Fourier coefficients and volume Omega, a local or
core contribution has force `-Omega sum_G G Im(conj(field_G) v_I,G)`.
The field is valence density for the local pseudopotential and V_xc for the
moving core. G=0 contributes no displacement derivative. Core Fourier
coefficients use the same truncated density sphere and radial transform as
SCF. Collinear core derivatives use `(V_up + V_down)/2`, since half of the
spin-neutral core enters each channel. PBE uses the full discrete variational
potential, including its divergence, with the existing Nyquist convention.
The XC density-floor chain-rule mask, when present, is retained.

For each projector overlap p=<beta|psi>, dp/dR=i<q beta|psi>, q=k+G.
The nonlocal force is `-2 sum_kn w_k f_kn Re(dp^dagger D p)`. The full D matrix,
all magnetic channels, all radial projectors and both spin channels participate.
There is one derivative traversal of each band's projector coefficients; setup
reconstructs projectors, and storage stays linear in plane waves and atoms.
No dense plane-wave nonlocal matrix is formed.

Ewald uses the energy evaluator's eta and real/reciprocal truncation spheres.
Self and background terms have no positional derivative at fixed cell. The
force sign is negative energy gradient. Kinetic, Hartree and valence XC have
no additional explicit ionic derivative at a stationary fixed-basis state.
For FD the objective is F=E-TS at fixed theta,N and any fixed moment; the
stationary occupation response cancels, so no extra df or entropy force is
added. Reported occupations include spin capacity and normalized k weights.

The exact spatial subgroup used by SCF averages each term as a polar vector
using A R A^-1 and its atom permutation. It does not subtract net force or
project optimizer constraints. Symmetry-off applies no averaging.

## Preconditions and interface limits

Standard-profile and uncertified states reject force queries. Tolerances may
not exceed density 1e-7 and energy/eigen 1e-8 Ha. Fixed, collinear and FD results
have separate views; density-only state does not suffice. The caller must not
modify a public result's arrays/certificate and then reuse certification.

CLI `[output] forces=true` writes a force block and text table. Fixed-spin and
FD output also write `forces.txt`. ASE `get_forces` converts Ha/Bohr to eV/Angstrom
using dft-units conversion functions, and rejects unconverged states. Its
existing Gamma/unpolarized model limits remain; stress is unsupported.

`docs/schemas/results-state.proposal.json` now accepts extension versions 1
and 2. Version 1 retains its original unavailable-only contract. Version 2 adds
`forces-available.proposal.json`, with strict/converged/certificate-passed
preconditions, Cartesian arrays, full term breakdown, net force and maximum
component. Available forces cannot appear under version 1. Quantity is E for
fixed occupations and F for FD. Unavailable branches prohibit all numeric
payloads and preserve the old `not_implemented` reason while allowing explicit
`not_requested`, failure and incomplete-coverage reasons. CLI spin/FD output
emits extension_version 2; the outer schema_version stays unchanged. The
legacy scalar profile has no state_contracts envelope, so its additive force
block validates directly against the same available-force schema.

The dependency-free validator used by existing spin output checks now also
checks force schemas. Seventeen accepted/rejected full-envelope fixtures cover
v1/v2 compatibility, fabricated nonzero unavailable forces, missing terms,
wrong units/shape/quantity, standard/unconverged/uncertified requests and boolean
preconditions. Cross-array lengths and force sums are independently checked on
actual CLI output. No jsonschema package is required.

## Reproduction and acceptance

Fast tests:

```
cargo test -p dft-engine --test forces
cargo test -p dft-cli --test forces
```

Release-only scans (each displacement holds cell, basis sphere and FFT grid
fixed; every atom/Cartesian component uses h=0.01,0.005,0.0025 Bohr):

```
cargo test --release -p dft-engine --test forces converged -- --ignored --nocapture
cargo test --release -p dft-engine --test forces spin_core_and_finite -- --ignored --nocapture
DFT_NC_DIR=/external/pseudodojo cargo test --release -p dft-engine --test forces real_files -- --ignored --nocapture
```

Synthetic cap is 1e-6 Ha/Bohr, net force 1e-8 Ha/Bohr on a resolved grid.
The trend check requires second-order improvement or a roundoff-scale plateau;
nonzero components explicitly check the sign. Ewald tests all axes, two eta
values and cubic/skew cells. Frozen local/nonlocal tests include shifted k,
two s projectors with off-diagonal D and p magnetic channels. A separate
molecule-like cubic test reimplements scalar Ewald via fixed integer boxes
and explicit pair cosines and local energy via an analytic Gaussian transform.
It never calls production Ewald/local energies to generate its oracle.

Converged tests include local-only and coupled nonlocal, NLCC, unpolarized
LDA-PW, fixed-sector spin and FD free energy on a skew cell. PBE+NLCC has a
separate grid refinement because the nonlinear quadrature's translation
residual does not vanish merely by converging SCF. Real Si and MgO use
PseudoDojo inputs with O/Si displaced by (0.03,0.02,0.01) Bohr. Real cap is
1e-4 Ha/Bohr only after the sum of independent SCF, cutoff and FFT-grid
force changes is below one third of the cap. Failing that condition, or an
unresolved translation residual, is INCONCLUSIVE and fails the test.
No QE force reference was used or generated.

Initial failures are retained in evidence: coarse synthetic grids gave net
forces up to 2.78e-6 Ha/Bohr; 48-cubed resolved the 2 Ha local/nonlocal tests; the 2/3 Ha scan uses 48/96-cubed respectively.
PBE+NLCC required further refinement. MgO cutoff 20->30 Ha differed by
3.8624e-3 Ha/Bohr, so its scan was refined to 40->60 Ha. No cap was enlarged.
Raw command outputs and on/off benchmarks are in `tests/reference/forces/`.

## Cost and legacy path

No SCF kernel, iteration, default or tolerance was changed. Forces execute
only on explicit post-convergence requests. The four-case force-off/on
benchmark records identical energies and density hashes. These are same-build
comparisons against the stored force-off record, not a claim that today's
engine matches the historical pre-optimization baseline byte for byte.

`measure_cost.py` records the input, UPF SHA-256, per-term forces, band visits,
force-pass time and SCF time for si8 and k444s. One shared-host Mac observation
was 1.007871 s / 20 band visits for si8 and 0.297878 s / 256 visits for k444s;
SCF times were 2.553888 s and 2.439059 s. Force share of SCF+force was therefore
28.30% and 10.88%. These include reconstruction/setup and are observations,
not dev08 acceptance or a performance optimization claim.

## Remaining owner work outside ticket ownership

The architect authorized and this follow-up implements the production schema
change and independent model checks. The previous schema integration blocker
is resolved. The README owner still extends `scripts/make_readme_figures.py`
for force records and regenerates README/images. Session/dashboard artifacts
remain outside OWNED_FILES; the orchestrator owns their registration/upload.

## Recorded result summary

The checked command/exit-code list, fixture hashes, external UPF hashes and
numerical observations are in `tests/reference/forces/manifest.json`.

| Check | Maximum observed discrepancy (Ha/Bohr) |
| --- | ---: |
| Converged synthetic local/nonlocal, 2/3 Ha | 3.05415e-8 |
| Frozen-orbital coupled nonlocal | 5.36098e-7 |
| Independent molecule-like scalar functional | 5.44023e-8 |
| NLCC LDA-PW / PBE | 5.76727e-8 / 5.72958e-8 |
| Spin/unpolarized reduction | 9.92262e-16 |
| FD free-energy derivative, theta=0.05 Ha | 2.69674e-8 |
| Distorted Si / MgO fixed-basis derivative | 5.42589e-9 / 2.80724e-8 |

The Si and MgO cutoff+grid+SCF uncertainty sums are 5.55474e-7 and
3.11319e-5 Ha/Bohr, respectively, below 1e-4/3. MgO's minimum grid has a
7.42e-6 Ha/Bohr net component; the explicitly refined grid reduces it to
3.40e-9. Thus the net-force qualification is for the resolved grid, not a
claim that every default grid has that accuracy.

The earlier composite NLCC/PBE/spin/FD run executed every PBE component
successfully but exited 101 because the 24-band FD retained space had too few
empty bands: the occupation certificate rejected its non-negligible top-band
tail. The documented FD force configuration is **theta=0.05 Ha, 48 bands,
Davidson subspace 96**. Subspace 48 cannot support a 48-band Davidson solve and
rejects at setup, so it was enlarged as well. This increases representation and
solver capacity, not a tolerance adjustment or force-fitting choice. All force,
SCF and finite-temperature certificate caps are unchanged. The corrected scan
has first occupation 1.8727364 electrons and TS=0.0563488778 Ha, distinguishing
F from a zero-entropy limit. `DFT_FORCE_CASE=nlcc|pbe|spin|fd` selects a subcase.
The architect accepted the composite-failure explanation; the unchanged long
PBE derivative scan is not repeated for this FD-only correction.

Follow-up checks re-run workspace gates, the corrected FD force scan, real-file
force/uncertainty scans, independent model rejection tests (including changed
FD temperature, constraint and band cap), actual CLI schema validation, native
Python/ASE tests, legacy force-on/off bit identity and README freshness. Results
and exit codes are recorded in `tests/reference/forces/followup.json`.
