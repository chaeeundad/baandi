# dft-physics

**Purpose:** Kinetic, local ionic, nonlocal KB, Hartree, LDA XC, density; matrix-free H*psi.

**Invariants to preserve and test:** Hartree G=0 convention; XC v = d(n*eps)/dn finite-difference test; H Hermitian; density integrates to N_valence.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-physics`, then `scripts/check.sh` for the whole workspace.

## MF03 polarized XC (P13B)

`tests/all/spin_oracles.rs` validates the vendored xcx 0.3.1 polarized LDA exchange,
standard PW92 correlation, and PBE exchange/correlation boundary. Production `spin.rs` supplies tagged density, LDA-PW/PBE spin XC, and reusable
spectral/Poisson buffers. The engine owns spin SCF. Run:

```sh
cargo test -p dft-physics --test all spin_oracles:: -- --nocapture
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

`tests/all/spin_grid.rs` checks functional perturbations on 12x14x16 and 15x12x10
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
Run `cargo test --release -p dft-physics --test all stress_terms:: -- --include-ignored --nocapture`.

## P26C real Gamma Hamiltonian (R3, library only)

`real_gamma::{RealNonlocal, RealHamiltonian}` uses R1 ORC and deterministic real
GEMM, R2 paired local transforms, and a duplicated kinetic diagonal. Construct
`RealHamiltonian::from_complex` with the exact source basis, block capacity and
persistent worker pool. It snapshots the real potential and projectors; later
updates of the source Hamiltonian do not change the snapshot. Applications
borrow caller-owned `RealMatrix` blocks and reuse construction-sized scratch.
No engine dispatch or existing complex arithmetic changes. Direct
`RealNonlocal::new` requires a matching Gamma `HalfSphere` from the caller.

Projectors are converted from public complex representative entries only. The
G=0 imaginary allowance is 1e-14 times the column maximum (zero columns accepted).
Per-atom row-major coupling blocks remain separate between P^T X and P W.
The 1e-13 relative component and Hermiticity gates allow FFT pairing and GEMM
summation rounding; they are not SCF or external-reference tolerances. Pair
norms must satisfy R2's factor-1e3 debug assertion. Local pairs are serial;
real GEMM uses the persistent pool with its existing fixed panel decomposition.

```sh
cargo test -p dft-physics --test all real_gamma_hamiltonian:: -- --nocapture
cargo test --release -p dft-physics --test all real_gamma_hamiltonian::large_parallel -- --ignored --nocapture
cargo test --release -p dft-physics --test all real_gamma_hamiltonian::timing_blocks -- --ignored --nocapture --test-threads=1
```

The latter two tests are ignored in the fast debug gate. Ten small synthetic
shapes test 1–20 atoms, mixed species, repeated s/p/d/f radial channels, off-diagonal
couplings and odd/even bands. The large test crosses the actual parallel GEMM
threshold and compares every output bit at T=1,2,4,8. See
`docs/performance/real-gamma-r3-timing.md` for observations and scope limits.

Opt-in apply switches (2026-11-11, docs/performance/requests/h-apply.md):
`BAANDI_HAPPLY_OVERLAP=on` (local and nonlocal image in one pool join,
bit-identical) and `BAANDI_HAPPLY_FOLD=on` (odd tail column through the
half-length folded transform of `dft_numerics::fused_fold`, rounding class);
tests `apply_overlap_is_bit_identical`, `apply_fold_matches_pairs`, ignored
microbenchmark `happly_widths`.

## M2 fused band density (2026-11-16)

`DistributedRealHamiltonian::band_density_sum` accumulates the plane density in
ascending band order through the sink form of `band_planes` (bit-identical to
summing `band_densities`); the engine takes it on the standard path with
`BAANDI_M2_DENSITY=fused` (default bands). Test
`distributed_band_values_match_serial_and_band_densities`.

## M2 kept start (MPI-KEPT-START, 2026-11-16)

`band_densities_keeping` / `band_density_sum_keeping` copy the pair slabs of the
leading pairs (density bits unchanged); `apply_local_change` applies a potential
change only (no kinetic, no nonlocal), kept pairs forward-transformed from their
slabs, the rest both ways; `potential()` accessor. Counters
`inverse_pair_transforms`, `forward_pair_transforms`, `kept_plane_pairs`,
`potential_change_columns`. Test `potential_change_from_kept_planes_matches_plain_apply`.

## Packed projector table and grid reuse (INIT-PACKED / INIT-REUSE, 2026-11-18)

`NonlocalOperator` stores every atom's projector columns once, in one
column-major table (`rows` x channels, atom order) shared through an `Arc`;
`AtomProjectors::columns()` returns a borrowed `ProjectorColumns` view
(`len`, `iter`, `column(j)`, `[j]`). The P48C block path multiplies that table
directly, so `set_block` only sets a flag and no second packed copy exists
(`block_bytes` reports the shared table when the block path is on). The
physical constructors fill the table directly (same products, same bits:
`tests/setup_memo_identity.rs` compares against a verbatim old constructor);
fixtures use `NonlocalOperator::from_columns`. `Hamiltonian::new_with_grid`
takes the caller's k grid (checked in O(N_G) against `map_miller`) instead of
building the density sphere again; test
`new_with_grid_matches_new_and_checks_the_grid`.

## Ewald schemes and ionic phase tables (INIT-EWALD-FAST, 2026-11-18)

`ewald::Scheme` (`BAANDI_EWALD_FAST`, per-thread `set_scheme_override`): the
default `Reference` keeps eta sqrt(pi)/Omega^(1/3) and the historical cutoffs;
its real-space images go through `ImageWalker` (conservative pruning of the old
box, exact old predicate, same order), bit-identical to eecd57f for energy,
forces, stress and every plan partition (verbatim copy in
`src/ewald/tests/old_reference.rs`). `Balanced` (rounding class) picks eta from
the cost balance (`BALANCE_COST_RATIO`), cutoffs from rigorous tail bounds,
pairs i <= j and the reciprocal half space with phase tables. Do not change the
reference scheme's arithmetic or order. `BAANDI_IONIC_FAST`
(`hamiltonian::set_ionic_fast_override`): `PhaseTables` structure factor in
`ionic_local*`, built from the full density list so row blocks stay
bit-identical. Release accuracy check:
`cargo test --release -p dft-physics --lib ewald::tests -- --include-ignored --nocapture`;
benches `examples/ewald_bench`, `examples/ionic_bench <UPF>`.

Independent accuracy reference (INIT-EWALD-INDEPENDENT-REFERENCE, 2026-11-18):
`src/ewald/tests/independent.rs` is a test-only double-double evaluator that
shares no code with either scheme (Madelung, supercell, finite-difference and
eta-independence checks; si32/si64 benchmark geometry in
`benchmark_positions.rs`). It shows that the default scheme's naive real-space
sum loses 1.8e-9 Ha at si64 (one-sided rounding of millions of sub-ulp terms)
while `Balanced` is within 3e-14 Ha. Keep it independent: do not route it
through engine helpers. Runs in the debug gate (1.9 s).
