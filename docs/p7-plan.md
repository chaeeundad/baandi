# P7 plan (started 2026-11-03, after MVP approval)

## Current status — 2026-11-05

P7A/P7A2 k-point engine/CLI and P7B PBE through xcx are merged. ADR 0018
accepted xcx 0.3.1; library selection is no longer open. R6 covers unshifted and
shifted 2x2x2 plus shifted 4x4x4, not unshifted 4x4x4. R8 covers Gamma and
shifted 2x2x2 PBE with an explicit LDA-UPF override. P7C/P10E produced GPU
design documents only; there is no GPU backend. The original plan below records
intended gates, not evidence that every proposed comparison was performed.
See [evidence map](validation/evidence-map.md) for commands and evidence limits.


Goal: extend the verified Gamma-only LDA engine without weakening any MVP gate. Each step keeps the Gamma-only result bit-for-bit comparable (E = -7.295530468210927 Ha at 30 Ry, see docs/validation.md) and is validated against QE 7.6 through the dft5090 session on dev08 (never a local pw.x).

| Ticket | Content | Gate |
|---|---|---|
| P7A | Multiple k-points: Monkhorst-Pack grid, time-reversal reduction only (k and -k), per-k plane-wave basis |k+G|^2/2 <= Ecut, per-k kinetic and nonlocal terms, weighted density and energy sums | Gamma-only unchanged to 1e-10 Ha; independent checks below; QE comparison on Si with a 2x2x2 and a 4x4x4 shifted and unshifted grid |
| P7B | PBE exchange-correlation through a maintained library (candidate: a pure-Rust crate or libxc FFI adapter; license, maintenance and Apple Silicon checks first, ADR required) | Library values against the in-tree LDA path on the LDA limit; QE PBE comparison |
| P7C | GPU backend design note only (numerical traits that a Metal or CUDA backend would implement; no GPU code) | ADR reviewed by the user |

Order: P7A, then P7B, then P7C. P7B needs a library decision that Claude prepares and the user approves if it adds FFI or a new license family.

## P7A design constraints
- Symmetry: no space-group reduction in P7A. Only time reversal (rho_k = rho_-k for a real density) is used, so the reduced set is exact without density symmetrization.
- The Gamma point must reproduce the existing path exactly (same basis, same energies).
- Independent checks (no QE needed): (1) a one-k-point run at k = Gamma equals the MVP result; (2) shifting k by a reciprocal lattice vector G0 gives the same energy and eigenvalues (basis relabeling invariance); (3) a free-electron-like synthetic potential has eigenvalues |k+G|^2/2 exactly; (4) the full unreduced grid and the time-reversal-reduced grid give the same energy and density to 1e-10.
- Convergence in the number of k-points is reported, never tuned to match QE.
- k-point weights sum to 1; occupations remain spin-unpolarized insulator-like (no smearing).

## Remaining provenance question
- PBE library choice is resolved by ADR 0018; xcx is merged.
- Whether P7 also covers a second pseudopotential; provenance-verified Si UPF would help (see docs/mvp-review.md).
