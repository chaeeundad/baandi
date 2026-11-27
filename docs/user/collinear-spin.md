# Fixed-sector collinear spin (MF03 / P13B)

The engine supports collinear NC SCF, including optional NLCC, with `lda-pw` (Slater exchange
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
initial_moment_seed = 1.0  # does not change the fixed electron sector
# or per atom in [structure]: `Fe 0 0 0 magmom=1.0` (the sum becomes the seed)

[basis]
# Existing cutoff/mesh entries go here.
kgrid = [2, 2, 2]
kshift = [0.5, 0.5, 0.5]
# ktime_reversal defaults to true for a grid: k ~ -k within each spin channel
# (ADR 0050); set false for the full mesh.
```

The example sector presupposes five valence electrons; it is not a material
input or a pseudopotential recommendation. The axis is a labeling convention in
this SOC-free model. `M=N_up-N_down` is reported in conventional spin-polarization
mu_B, not signed electron angular magnetic moment. Seeds are recorded but do not
change the sector or model identity. Omitted `[spin]` or `mode="unpolarized"`
without additional spin keys follows the unpolarized reference path.

## Calculation and certificate

The production entry point is `dft_engine::refscf::api`, which dispatches to
`refscf::electronic` for spin/FD. The two density channels use the reference
Broyden/plain update (beta 0.7/history 8 by default) and the distributed complex
solver. `CollinearProblem` is a feature-gated strict-oracle compatibility API.
Opt-in identity-spin spatial symmetry is supported; magnetic spin swaps and
magnetic time reversal are not inferred.

The Hartree field uses total density. XC uses the two valence channels and
three GGA invariants `(sigma_uu,sigma_ud,sigma_dd)`. xcx returns energy per
electron; the adapter multiplies by total density. The up-channel divergence
uses `2*e_sigma_uu*grad(n_up) + e_sigma_ud*grad(n_down)` (one cross factor).
Direct and band-sum energy audits subtract the potential integrated against
valence output density, and include Hartree, XC and Ewald once.

**Default path (Stage R).** Fixed-sector LDA-PW/PBE and qualified FD inputs use
Reference with independent fixed-point verification, threshold 5e-9 Ha by default.
Free-moment collinear FD PBE remains unqualified (G4D) and routes to StrictOracle.
`BAANDI_REFERENCE_STRICT=on` explicitly selects that oracle with precedence;
there is no automatic retry after a numerical failure. Current spin metric,
stiffness, number/moment and verification contracts are in
[reference SCF](../design/reference-scf.md#spin-and-fermi-dirac) and
[ADR 0047](../decisions/0047-reference-driver-single-path.md).
The strict certificate described next belongs to the oracle, not the production
reference driver's estimator. Its tolerances are not accepted legacy TOML keys.

Strict reference: convergence requires raw total-density and magnetization residuals, energy,
number/moment conservation and strict eigen residuals in both channels. A fresh
map at `0.1*eigen_tolerance` must also pass both raw residuals, energy change,
occupied eigenvalue change (`1e-7 Ha`, the legacy certificate cap), strict
occupied/first-empty residuals, channel/total/moment sums and channel gaps.
Residual norms are `sqrt(integral delta^2)/N` in Bohr^-3/2, normalized by total N.
Higher empty states retain the existing buffer tolerance. The returned candidate
orbitals, density and energy all refer to the same state; fresh-map metrics are
reported separately. Budget exhaustion returns a retained unconverged state and
CLI exit 3. It never fabricates a passed certificate.

Reference initialization uses atomic charge/orbitals when available. Missing
atomic data use the documented neutral uniform/deterministic kinetic completion,
with executed start provenance. Explicit start values accept only `"atomic"`;
`"uniform"`/random starts are not current CLI options. Gamma-real and the
reference cutoff ladder remain restricted to their supported scalar models.

The architect approved xcx result-vector allocations at density-field evaluation
boundaries, outside orbital, k and Hamiltonian inner loops. Output-state energy
assembly and fresh certification each need their own XC evaluation as in the
legacy energy audit. Existing solver and mixer result/history allocations are
unchanged. No fully allocation-free solver claim or performance gain is made.

## Output and deliberately unavailable combinations

Reference collinear/FD runs write schema version 2,
`output_profile="reference_spin_fd_v1"` and `reference_scf.schema="reference-spin-fd-1"`.
They retain spin/k/band axes, occupations and E/TS/F data as applicable, plus
routing and verification. The older `collinear_fixed_v1` profile and
`baandi-mf03-record-v1` identifiers describe oracle/historical evidence.
Spin/FD output currently lacks the scalar driver's map-boundary clocks.

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

Reference spin/FD supports NLCC, eligible spatial symmetry, verified forces/stress,
observe/agent runtime and exact checkpoint v2. See [CLI](../cli.md) and
[checkpoint limits](../checkpoint-format.md). Free-moment PBE still uses the
G4D oracle, whose spin runtime/restart restrictions remain. CLI spin/FD relaxation,
spin band/DOS postprocessing, noncollinear/SOC and U remain unsupported.
Unsupported outputs reject rather than fabricating data.

## Historical P13B reproduction and evidence

The measurements below retain their original oracle/legacy implementation and
are not a new Stage R material qualification. Current evidence is indexed in ADR 0047.

```sh
scripts/check.sh
cargo test -p dft-physics --test all spin_oracles:: -- --nocapture
cargo test -p dft-physics --test all spin_grid:: -- --nocapture
cargo test --release -p dft-engine --test all collinear:: -- --include-ignored --nocapture
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
