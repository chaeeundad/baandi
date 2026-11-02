# ADR 0007: XC functional implementation and special functions

Status: Accepted (2026-11-02, Claude Code). Applies ADR 0003.

## Decision
- LDA-PZ (Slater exchange with alpha = 2/3 plus Perdew-Zunger 1981 correlation, the combination QE calls `sla+pz`) is implemented in `dft-physics` as closed-form functions. It is physics glue that must match one specific parameterization exactly; no surveyed library offers PZ81: `xcx` 0.3.1 has LDA_X, LDA_C_PW and VWN variants but no LDA_C_PZ (checked in the vendored source).
- `xcx` is a dev-dependency only, used to cross-check the exchange part to tight tolerance and to sanity-check correlation against PW92 (different parameterization, so only a loose bound). Its license is `(MIT OR Apache-2.0) AND MPL-2.0` (MPL part derived from libxc, file-level copyleft); as a test-only dependency it is not linked into released binaries. Re-review before it moves to `[dependencies]`.
- `libm` (MIT OR Apache-2.0) provides `erfc` for the Ewald sums; no hand-written special functions.

## Consequences
XC correctness rests on finite-difference derivative tests, continuity of the PZ81 branches at rs = 1, the exchange cross-check against xcx, and later the P6 comparison with QE. A future PBE phase should adopt a library (libxc binding or xcx) as the primary source.
