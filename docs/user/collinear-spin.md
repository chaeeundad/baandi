# Fixed-sector collinear spin (MF03 / P13B)

The engine supports collinear valence-only NC SCF with `lda-pw` (Slater exchange
plus standard Perdew–Wang 1992 correlation) and `pbe` through xcx 0.3.1.
Both functionals also support unpolarized calculations. Existing `lda-pz` remains
unpolarized only; requesting it with collinear spin returns:

```
xc.spin: required polarized capability unavailable for lda-pz; select pbe or lda-pw
```

This is a fixed integer sector, not an unconstrained search for a magnetic ground
state. Supply nonnegative integer `n_up` and `n_down` whose sum is the neutral UPF
valence count, or supply `total_moment` such that `(N+M)/2` and `(N-M)/2` are
nonnegative integers. Occupation capacity is one per spin/k/band; occupations
exclude k weights. Each nonempty channel needs an empty band and a positive
sampled indirect gap above `max(1e-7 Ha, 2*eigen_tolerance)`. An empty channel
still has a computed retained spectrum, but its occupied-edge gap is null.

```toml
[xc]
functional = "lda-pw"  # or "pbe"; UPF XC must match, unless explicitly overridden

[spin]
mode = "collinear"
n_up = 3
n_down = 2
# Alternatively replace both counts with: total_moment = 1.0
axis = "z"             # x, y, z; default z
# Alternatively: axis_cartesian = [0.0, 0.0, 1.0] (unit vector)
initial_moment_seed = 1.0  # provenance only; initialization uses uniform fixed counts

[basis]
# Existing cutoff/mesh entries go here.
kgrid = [2, 2, 2]
kshift = [0.5, 0.5, 0.5]
ktime_reversal = false   # required for collinear meshes; no symmetry reduction
```

The example sector presupposes five valence electrons; it is not a material
input or a pseudopotential recommendation. The axis is a labeling convention in
this SOC-free model. `M=N_up-N_down` is reported in conventional spin-polarization
mu_B, not signed electron angular magnetic moment. Seeds are recorded but do not
change the sector or model identity. Omitted `[spin]` or `mode="unpolarized"`
without additional spin keys follows the legacy scalar path.

## Calculation and certificate

`dft_engine::collinear::CollinearProblem` is an additive API. It accepts a
`SpinSector`, `ScfOptions`, and a full k quadrature. The caller must not submit a
symmetry-reduced quadrature; the CLI explicitly rejects time reversal reduction.
All common solver controls use the established eigensolver implementation. New
SCF density, gradient, FFT, potential, orbital-copy, nonlocal-projection and
history storage is reserved at setup and reused. The two spin channels use the
same solver and the existing Linear, Pulay or KerkerPulay mixers. A single
feasibility fraction applies to both channel proposals; no output density or
magnetization is clipped or rescaled after diagonalization.

The Hartree field uses total density. XC uses the two valence channels and
three GGA invariants `(sigma_uu,sigma_ud,sigma_dd)`. xcx returns energy per
electron; the adapter multiplies by total density. The up-channel divergence
uses `2*e_sigma_uu*grad(n_up) + e_sigma_ud*grad(n_down)` (one cross factor).
Direct and band-sum energy audits subtract the potential integrated against
valence output density, and include Hartree, XC and Ewald once.

Convergence requires raw total-density and magnetization residuals, energy,
number/moment conservation and strict eigen residuals in both channels. A fresh
map at `0.1*eigen_tolerance` must also pass both raw residuals, energy change,
occupied eigenvalue change (`1e-7 Ha`, the legacy certificate cap), strict
occupied/first-empty residuals, channel/total/moment sums and channel gaps.
Residual norms are `sqrt(integral delta^2)/N` in Bohr^-3/2, normalized by total N.
Higher empty states retain the existing buffer tolerance. The returned candidate
orbitals, density and energy all refer to the same state; fresh-map metrics are
reported separately. Budget exhaustion returns a retained unconverged state and
CLI exit 3. It never fabricates a passed certificate.

The architect approved xcx result-vector allocations at density-field evaluation
boundaries, outside orbital, k and Hamiltonian inner loops. Output-state energy
assembly and fresh certification each need their own XC evaluation as in the
legacy energy audit. Existing solver and mixer result/history allocations are
unchanged. No fully allocation-free solver claim or performance gain is made.

## Output and deliberately unavailable combinations

Collinear runs write `results.json` with `schema_version=2` and profile
`collinear_fixed_v1`; unpolarized output remains schema 1 without a new spin
block. The spin block follows `spin.proposal.json`; `state_contracts` follows
`results-state.proposal.json` including fixed occupations, zero entropy, E=F and
unavailable forces/stress. `proposal_status="proposed"` is retained because it is
a constant in that frozen extension schema, not a claim that the engine emitted
no result. Model/basis/state digests use the versioned binary record encoding
`baandi-mf03-record-v1` (sorted named fields, type tags, length-prefixed UTF-8,
little-endian integers/f64, normalized signed zero). Records accompany digests;
initial seeds and actual converged state are separate from model identity. These
are result identifiers, not support for MF07/08 warm/exact spin restart.

Band arrays are spin-major, then input k order, then band; capacity is `[1,1]`.
The total moment, both raw density residuals, measured channel electron counts,
per-channel gaps and strict residuals appear in the certificate. Channel densities
are included in JSON in C order (z fastest). The legacy scalar orbital view is
ambiguous for collinear results and explicitly errors.

Four Gaussian cubes are written first: `density-up.cube`, `density-down.cube`,
`density-total.cube`, `magnetization.cube`. Coordinates are Bohr and scalars are
respectively electrons/Bohr^3 and signed electron difference/Bohr^3. The origin
is zero and z is contiguous. Known atomic numbers come from UPF element metadata;
if synthetic atomic numbers are unavailable, grid-only cubes explicitly omit
atoms and JSON retains the structure and omission flag. `density.bin` and
`density.json` additionally retain the documented legacy total-density binary
layout. No spin information is inferred from that total-density file.

Spin-aware symmetry, smearing, noncollinear/SOC, NLCC, U, forces/stress, runtime
control/telemetry, spin restart and spin band/DOS/legacy volumetric postprocessing
are not supplied here. Explicit requests for unsupported CLI combinations reject
with a key, reason and remedy. Classic collinear SCF still writes its standard
channel cubes. Real-material comparisons require separately matched, licensed
pseudopotentials and a requested reference run; none is claimed by this ticket.

## Reproduction and evidence

```sh
scripts/check.sh
cargo test -p dft-physics --test spin_oracles -- --nocapture
cargo test -p dft-physics --test spin_grid -- --nocapture
cargo test --release -p dft-engine --test collinear -- --include-ignored --nocapture
cargo test -p dft-cli --test spin -- --nocapture
python3 -I scripts/tests/test_readme_figures.py
```

The full geometry/mesh matrix is ignored in debug (about 8 seconds in release).
Its fixture uses `tests/common/mod.rs`'s analytic synthetic Gaussian NC data,
401 radial points, cutoff 2.5 Ha, widths 1 Bohr, charges 1 or 2, zero or three
projectors, 2/4 atoms, 8/16 bands, the four approved cubic/skew/elongated/slab
cells, Gamma/full half-shifted 2x2x2/unequal 0.2,0.3,0.5 quadratures. Precise
positions, generation parameters and solver settings are in the test. Fixture
source hashes and measurements are retained in
[`collinear-evidence.json`](../../crates/dft-engine/tests/collinear-evidence.json).
The reduction iteration cap was frozen at a difference of at most one; all
23 converged pairs had exactly equal iteration counts. One skew unequal-weight
PW92 case is noninsulating: both paths reject with gap approximately
`-4.90921617e-4 Ha`; this is rejection evidence, not a converged energy result.
The 23 accepted pairs have maximum energy difference `6.66133815e-15 Ha` and
normalized density distance `6.03926407e-14 Bohr^-3/2`.

The odd one-electron Gamma fixture converges with M=+1 and M=-1, with exact
energy/channel-swap equality for both functionals. Its energies are
`-0.48409468993977545 Ha` (PW92) and `-0.48439633882319222 Ha` (PBE).
A skew three-electron fixture exercises both occupied channels (2/1) and all
three mixers: energies agree within `1e-10 Ha`, with both gaps certified.

The PW92 correlation oracle is independently transcribed from the
architect-supplied Perdew and Wang, Phys. Rev. B 45, 13244 (1992) equations,
including standard rounded coefficients and f''(0)=1.709921. Its maximum combined
exchange/correlation error is `7.07712968e-15 Ha/electron`; interior potential
error is `1.38815626e-11 Ha`. PBE correlation absolute values are validated only
through spin reduction/scaling limits and derivative consistency, **not** an
independent transcription or matched external-code run. Exact exchange scaling,
spin swap, three-step functional perturbations and the cross-gradient factor
are separately tested. Full-grid functional derivatives on skew even/odd grids
have maximum error `6.58009064e-12 Ha`. No acceptance tolerance was relaxed.

The unchanged `scripts/bench_scf.py` compares against a same-machine pre-change
binary built from `fa5de50197a7a10a5abef37aa3914d8e3a681a6b`. Gamma, k222s,
k222 and k444s all have `delta E=0 Ha`, bit-identical density hashes and ten
iterations. Raw before/after benchmark records and exact command are retained
in the evidence JSON. Single shared-host timings do not establish a speedup or
an allocation performance effect. The historical performance baseline predates
later accepted arithmetic changes; this ticket tests identity to its actual
parent implementation instead of silently replacing that historical baseline.

README figure-generator/reference integration and dashboard artifacts are outside
this ticket's owned paths. The orchestrator can import the retained evidence into
`tests/reference`, extend `scripts/make_readme_figures.py`, and regenerate the
README validation block; no hand-entered README claim is added here.
