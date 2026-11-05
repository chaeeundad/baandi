# dft-physics

**Purpose:** Kinetic, local ionic, nonlocal KB, Hartree, LDA XC, density; matrix-free H*psi.

**Invariants to preserve and test:** Hartree G=0 convention; XC v = d(n*eps)/dn finite-difference test; H Hermitian; density integrates to N_valence.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-physics`, then `scripts/check.sh` for the whole workspace.

## MF03 polarized XC (P13B)

`tests/spin_oracles.rs` validates the vendored xcx 0.3.1 polarized LDA exchange,
standard PW92 correlation, and PBE exchange/correlation boundary. Production `spin.rs` supplies tagged density, LDA-PW/PBE spin XC, and reusable
spectral/Poisson buffers. The engine owns spin SCF. Run:

```sh
cargo test -p dft-physics --test spin_oracles -- --nocapture
```

The independent PW92 oracle is transcribed from the architect-supplied equations
of J. P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992), including the standard
rounded coefficients and f''(0)=1.709921. It does not call xcx. Caps were frozen
before measurement: 1e-12 Ha/electron and potential relative 1e-8 or absolute
1e-10 Ha near zero. The test covers total densities 0.01, 0.1, 1 Bohr^-3 and
polarizations 0, 0.3, 0.9, 1. Central derivatives exclude an empty channel and
use three steps with a second-order trend or a resolved sub-cap plateau.
The initial sigma difference steps suffered roundoff; increasing the perturbation
resolved the signal without changing the oracle, library or acceptance caps.

Observed maxima on 2026-11-05: combined PW92/exchange energy error
7.07712968089913e-15 Ha/electron and interior potential error
1.3881562566098182e-11 Ha; PBE exchange scaling error
2.0122792321330962e-14 Ha/electron, density derivative error
5.2889137514000595e-11 Ha and sigma derivative error
8.586198418925051e-10 Ha Bohr^5 (within the relative cap).
Spin-swap potential differences were zero. Zero-polarization reduction here is
only the point functional boundary, not a complete SCF reduction certificate.
PBE correlation has no independently transcribed absolute-value oracle in this
test; only reduction, swap and derivative consistency are established. The
gradient perturbation explicitly tests the single cross-sigma flux factor.
These are synthetic functional checks, not matched external-code calculations.

The architect approved xcx's owned result-vector allocation once per density
field evaluation (input map, output energy audit, and fresh certificate), outside
orbital/k/Hamiltonian loops. All new adapter grid, spectral, local-potential,
orbital-copy and nonlocal-projection buffers are prepared once and reused.
Existing solver and mixer ownership APIs retain their existing allocations;
this is not a claim of a wholly allocation-free SCF. No production XC formula
is transcribed. PZ remains unpolarized only. No measured optimization is claimed.

`tests/spin_grid.rs` checks functional perturbations on 12x14x16 and 15x12x10
skew grids, zero-spin grid reduction, tagged density validation and retained
output storage. Maximum extrapolated grid-derivative error was
6.5800906390300895e-12 Ha. It also establishes bit-identical Poisson energy
and real potential against the existing full-grid evaluator on those grids.
See `docs/user/collinear-spin.md` for full SCF scope and evidence.

**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.

## MF06 NLCC (P13C)

`density::core_density` prepares immutable core Fourier coefficients at setup,
with shell caching, density-sphere truncation, 1/volume and ionic phases. It never
renormalizes or clips core density. Spin XC uses half the core in each channel
and differentiates the combined density; its valence potential integral excludes
core and its separate core integral is audit-only. Combined-density screens and
negative-density policies remain those of the selected functional. See the
engine `tests/nlcc.rs` for localized three-step differences, translation, analytic
multiple-projector assembly and converged spin reduction. No force/stress
capability is established and no new per-orbital allocation is introduced.

## MF04 occupation roots

`occupations.rs` is a pure, allocation-free number/entropy boundary with explicit
spin/k/band layout and caller-owned output. Keep the independent tanh/binary
entropy tests independent of the production stable-tail formula. Failed output
buffers are not states. See `docs/validation/occupations.md` for the bracket,
relative electron tolerance, sector endpoints and frozen oracle caps.

## MF11 crystal symmetry (P19A)

`symmetry::DensitySymmetrizer` prepares source indices, translation phases and persistent FFT buffers. Apply rho_sym(G)=mean exp(-2*pi*i*G.tau) rho(R^T G). Do not allocate inside application. The independent real-space field oracle must remain independent.

See `docs/validation/symmetry.md` and `tests/reference/symmetry/` for commands,
observed failures and validation limits. No QE symmetry-on oracle is available.

## Explicit ionic derivatives (P20A)

`ewald::forces` uses the energy's eta/truncation protocol. `forces::local_or_core`
contracts phase derivatives against a forward-normalized Fourier field, and
`forces::nonlocal` differentiates both sides of the full Hermitian D coupling.
These APIs are fixed-electronic-state partial derivatives; only the engine may
claim a complete stationary force. Species-grouped atom order must be restored
by callers (the engine supplies one species entry per input atom). Core forces
are the explicit XC derivative, not an additional core energy. Historical P13C
statements about absent force capability describe that earlier ticket's scope.

## Analytic stress (P22A)

`kinetic::stress`, `ewald::stress`, and `stress` implement tensile-positive
frozen-Miller strain derivatives in Ha/Bohr^3. XC gradient metrics live beside
their functional adapters. The engine alone certifies a complete total. Core
response is contracted against the XC potential, not added as a core energy.
All supported l=0..3 channels have a two-projector coupled-D fixture; solid
harmonic derivatives avoid polar-axis singularities. No numerical projector
finite differences or SCF hot-loop calls are used. See
`docs/validation/stress.md` for equations, conventions and scope.
Run `cargo test --release -p dft-physics --test stress_terms -- --include-ignored --nocapture`.
