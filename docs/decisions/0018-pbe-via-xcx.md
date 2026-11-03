# ADR 0018: PBE through xcx (proposal, awaiting approval to promote the dependency)

- Date: 2026-11-03
- Status: proposed. The only part that needs the user's approval is moving `xcx` from a dev-dependency to a normal dependency of `dft-physics` (it brings a file-level MPL-2.0 component into linked binaries).
- Applies ADR 0003 (prefer general-purpose libraries) and supersedes the "re-review before it moves" note of ADR 0007.

## Decision
Use `xcx` 0.3.1 (pure Rust, no C, already vendored) as the source of the GGA exchange-correlation values and derivatives for PBE (`GgaXPbe` and `GgaCPbe`: energy per particle, `vrho`, `vsigma`). Project code supplies only the physics glue that a plane-wave GGA needs and that no XC library provides: the density gradient by FFT, `sigma = |grad n|^2`, the potential `V_xc = vrho - div(2 vsigma grad n)` evaluated in reciprocal space, density thresholds, and the energy integral. LDA-PZ stays in-tree because no surveyed library has PZ81.

## Evidence (spike, scratch program, not in the repository; to be turned into repository tests by ticket P7B)
Unpolarized, Hartree units, seven (n, sigma) points from n = 0.01 to 1.0 and sigma from 0 to 3, compared with closed forms typed independently from the published PBE96 and PW92 formulas (kappa 0.804, mu 0.2195149727645171, beta 0.06672455060314922, gamma (1 - ln 2)/pi^2, PW92 parameters):
- maximum |exc_x(xcx) - exc_x(closed form)| = 2.2e-16 Ha; maximum |exc_c difference| = 1.4e-17 Ha.
- `vrho` and `vsigma` against central finite differences of `n * exc`: maximum relative difference 3.9e-11 (exchange) and 8.2e-11 (correlation) for `vrho`; `vsigma` agrees to the nine printed digits.
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
