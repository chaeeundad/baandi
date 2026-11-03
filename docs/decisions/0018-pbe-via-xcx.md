# ADR 0018: PBE through xcx

- Date: 2026-11-03
- Status: accepted 2026-11-03. The user approved moving `xcx` from a dev-dependency to a normal dependency of `dft-physics` (it brings a file-level MPL-2.0 component into linked binaries). Manifest changed by Claude; the PBE implementation is ticket P7B.
- Applies ADR 0003 (prefer general-purpose libraries) and supersedes the "re-review before it moves" note of ADR 0007.

## Decision
Use `xcx` 0.3.1 (pure Rust, no C, already vendored) as the source of the GGA exchange-correlation values and derivatives for PBE (`GgaXPbe` and `GgaCPbe`: energy per particle, `vrho`, `vsigma`). Project code supplies only the physics glue that a plane-wave GGA needs and that no XC library provides: the density gradient by FFT, `sigma = |grad n|^2`, the potential `V_xc = vrho - div(2 vsigma grad n)` evaluated in reciprocal space, density thresholds, and the energy integral. LDA-PZ stays in-tree because no surveyed library has PZ81.

## Evidence (spike, scratch program, not in the repository; to be turned into repository tests by ticket P7B)
Unpolarized, Hartree units, seven (n, sigma) points from n = 0.01 to 1.0 and sigma from 0 to 3, compared with closed forms typed independently from the published PBE96 and PW92 formulas (kappa 0.804, mu 0.2195149727645171, beta 0.06672455060314922, gamma (1 - ln 2)/pi^2, PW92 parameters):
- maximum |exc_x(xcx) - exc_x(closed form)| = 2.2e-16 Ha; maximum |exc_c difference| = 1.4e-17 Ha.
- `vrho` and `vsigma` against central finite differences of `n * exc`: maximum relative difference 3.9e-11 (exchange) and 8.2e-11 (correlation) for `vrho`; `vsigma` agrees to the nine printed digits.
The original seven points, confirmed by Claude for P7B, are (n in Bohr^-3,
sigma in Bohr^-8): `(0.01,1e-6)`, `(0.03,4e-4)`, `(0.03,0)`, `(0.1,0.05)`,
`(0.2,2)`, `(0.3,0.2)`, `(1.0,3)`.
This checks the library against the published formulas; it does not check a plane-wave GGA SCF, which needs its own validation (below).

## License and dependency review
`xcx` is `(MIT OR Apache-2.0) AND MPL-2.0`: original code MIT or Apache-2.0, libxc-derived functionals MPL-2.0 (file-level copyleft; source of those files must stay available; our own code is unaffected and stays MIT). Linking it into released binaries therefore adds an MPL-2.0 component to the notice set. Rust edition 2021, MSRV 1.87 (our toolchain 1.99.1), no native build, pure Rust so Apple Silicon is covered by the same tests. Dependencies pulled by `xcx` beyond what is already vendored: none (all vendored). Maintenance: single upstream author, young crate; mitigation: keep it behind a small `dft-physics` function interface and keep the closed-form oracle tests, so it can be replaced by another library or by in-tree code.

## Alternatives
- libxc through a C binding: native build and FFI (`unsafe` only in an adapter), MPL-2.0 as well, more platform risk: rejected.
- Hand-written PBE: rejected by ADR 0003 (physics glue only); kept as a fallback because the closed-form oracle already exists.

## Validation required before the feature counts (ticket P7B)
Closed-form oracle tests (as the spike), finite-difference check of `E_xc[n]` against `V_xc` on the FFT grid including the divergence term, rigid-translation invariance, `sigma = 0` limit, density-threshold behavior documented, then a QE 7.6 comparison with the same pseudopotential for both codes using QE's functional override (a code-to-code consistency test; the validation Si file was generated for LDA, so the physical meaning of the PBE result is not claimed) through the dft5090 session, with the tolerances of R3 and R6.
## Consequences
`docs/dependencies.md` must list the MPL-2.0 component for the linked set; the UPF functional label check needs an explicit, recorded override for the code-to-code test (never silent).

## P7B implementation results (2026-11-03)

The new GGA adapter batches each xcx component over the grid and uses a real,
skew-adjoint spectral derivative with each even-axis Nyquist component zeroed
before Cartesian conversion. The negative-density and 1e-10 Bohr^-3 screening
rules are stated in physics-conventions.md. LDA numerical source is untouched.
The exact discrete XC remainder extends the existing Harris-Foulkes difference
bound to GGA; a local-density quadratic bound is not asserted for PBE.

Independent PBE96/PW92 oracle tests use the exact seven original spike points
listed above. Maximum exc errors: exchange 2.220446049250313e-16 Ha,
correlation 2.0816681711721685e-17 Ha (limit 1e-14 Ha).
Maximum relative central-difference errors for (vrho,vsigma): exchange
(2.6960762705971615e-11,6.000294732823276e-11), correlation
(9.937766006440036e-11,4.520355414753177e-10), each below 1e-8.
The sigma step is 1e-5 times max(sigma,4*k_F^2*n^2), balancing cancellation
against central-difference truncation at zero and small sigma; the analytic
extension is used on the negative side at sigma=0. A fixed 1e-10 sigma step
at (0.03,0) lost relative precision (1.655e-8) through cancellation; no
acceptance tolerance was changed.

For periodic Gaussian sums, grids 12x14x16 (orthogonal cell) and 15x12x10
(skew cell) have variational derivative errors 8.138585509e-12 and
8.235417123e-12 Ha (limit 1e-8 Ha), translation errors 3.552713679e-15 and
8.881784197e-16 Ha (limit 1e-12 Ha), and uniform-limit energy errors
1.278976924e-13 and 1.421085472e-14 Ha (limit 1e-12 Ha).
Commands: `cargo test --offline -p dft-physics --test pbe -- --nocapture`.
These are mathematical fixtures, not QE comparisons.

Synthetic two-species SCF with PBE labels, 2 electrons and 2.5 Ha cutoff:
Gamma converges in 11 iterations to 1.2198550926939631 Ha; shifted 2x2x2
converges in 12 to 2.1981939378887070 Ha. Charge is
2.0000000000000031 electrons in both. Band-form errors are 8.881784197e-16
and 4.440892099e-16 Ha; Harris-Foulkes errors are 1.776356839e-15 and
3.996802889e-15 Ha, below exact-remainder/eigen residual bounds
5.412823312e-9 and 1.149250706e-14 Ha plus 1e-12 Ha scalar-roundoff allowance.
Both pass fresh fixed-point certification with unchanged default tolerances.
Command: `cargo test --offline --release -p dft-engine --test pbe -- --ignored --nocapture`.
The real-Si PBE/QE comparison remains deferred to Claude's dft5090 session.

The workspace build/test/fmt/clippy gate passes. The Python suite's existing
106 tests pass, including the stored real-Si density regression: rho(0)
difference 0 electrons/Bohr^3, sum-|rho(G)|^2 relative difference
1.1078135768010342e-8, largest sampled coefficient difference
4.008255174783675e-9 electrons/Bohr^3. These are existing LDA/QE data, not PBE.
A temporary snapshot harness run on baseline commit
`59a4b1cd4cb411fd431ff86548f5d143e458aa90` and the P7B implementation compared
all energy-audit scalar terms selected by the harness, charge, every density
grid value, and every k eigenvalue for synthetic Gamma and shifted 2x2x2 LDA
runs: all f64 bit patterns match exactly. SHA-256 of the printed snapshot lines:
`cf11c51c2567da3456c69200d3329e7f2a586558bd95ec59db79d41bf2aa0413` (Gamma),
`c804e7d67e6b476e208a3d3e53a8e617f5c981856a3d5a501a8736f5dfae14d7` (2x2x2).
This is a same-toolchain/machine comparison, not a cross-platform bit guarantee.

The authorized shared QE reader now accepts PBE results and checks explicit
boolean override permission, the applied flag, and their consistency with the
UPF label. The generator emits `input_dft = 'pbe'` and comments containing the
original UPF label and both flags. Missing permission for a mismatched label
is rejected. A byte-for-byte test preserves the existing LDA golden input;
all pre-existing LDA tests remain unchanged. No pseudopotential is opened or
modified. The actual PBE QE run and README evidence integration belong to
Claude; README.md and tests/reference are unchanged.

Final verification: `scripts/check.sh` exits 0 (build, test, fmt, clippy).
`cargo test --offline --release -p dft-engine -p dft-physics -p dft-cli -- --ignored --nocapture`
exits 0, including the original LDA ignored tests and documentation tests.
`python3 -I -m unittest discover -s scripts/tests` exits 0: 113 tests,
including seven new PBE interchange and LDA byte-identity tests. The initial uniform-energy test
failure was resolved by ndarray multi-accumulator summation, and the initial
legacy-output failures by preserving the LDA JSON shape; no tolerance changed.
