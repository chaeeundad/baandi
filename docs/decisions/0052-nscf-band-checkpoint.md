# ADR 0052: Many-band NSCF solvers and the `nscf_bands_v1` checkpoint unit (R1)

- Status: **Accepted** 2026-11-28 (ticket RPA-R1-MANY-BANDS, proposed 2026-11-27). Evidence: QE 7.6 nscf on dev08, five cases within 8.1e-6 Ha (gate 2e-5), `tests/reference/response/r1/dev08-2026-11-27/summary.json`.
- Scope: `refscf::nscf::nscf_bands_with` (engine), `dft_solvers::chebyshev_complex`,
  `DistributedComplexHamiltonian::{dense_matrix, spectral_upper_bound}`,
  `baandi nscf`, `docs/checkpoint-format.md` section "NSCF band units".
  Plan: `docs/roadmap/rpa-crpa-plan.md` 3.2, 3.5, 3.6 and ticket R1 (user
  decisions D1-D10 of 2026-11-26).

## Decision

1. Three band solvers behind one fixed-density entry point. `Davidson` is the
   unchanged reference path (`nscf_with` delegates to it bit for bit).
   `Dense` builds the explicit H_k replicated on every rank from gathered
   owned data (kinetic, projector rows) and one serial DFT of the gathered
   potential, indexing V_N(G - G') modulo the grid exactly as the FFT apply
   convolves, and runs the complete faer Hermitian EVD. The matrix has the
   same bits for every rank count, so each rank keeps its rows of identical
   eigenvectors without communication. `Chebyshev` is filtered subspace
   iteration (degree 10, SVQB orthonormalization twice, Rayleigh-Ritz,
   contiguous locking, buffer columns max(8, bands/10) never certified) with
   a rigorous upper bound max T + max V + max(0, max_atom lambda_max(D))
   lambda_max(P^H P).
2. Certificate: every returned band of every solver carries a fresh
   `||H psi - e psi||` from an independent `apply_block`; a band above the
   requested tolerance fails the k point.
3. Checkpoint: one BAANDICP container per completed (channel, k) with
   rank-indexed sections and a fingerprint map, manifest written last;
   resume re-certifies loaded bands and rejects any mismatch.
4. The minimal CLI surface is `baandi nscf` with an `[nscf]` table
   (`bands`, `band_solver`, `eigen_tolerance_ha`, `k_mesh`, Chebyshev
   controls); the plan's `[response]` keys follow with R2.

## Evidence (Mac, 2026-11-27; observations, not dev08 acceptance)

`cargo test -p dft-engine --test refscf_nscf_many` (debug, synthetic skewed
two-species, 3 k) and the ignored release matrix with Si.pz-vbc (Si2 4x4x4
full mesh, skewed Si 2x2x2) and PseudoDojo PBE (LiF 2x2x2, wurtzite AlN):
dense complete-basis residuals <= 5.2e-14 Ha, completeness
|X X^H - I|, |X^H X - I| <= 1.1e-14, Davidson vs dense <= 3.0e-14 Ha,
Chebyshev vs dense <= 3.8e-12 Ha with residuals <= 1e-6 Ha, rank counts
1/2/4: dense bit-identical, Chebyshev <= 1.4e-14 Ha, np=1 bit-identical to
the serial adapter. CLI Si2 30 Ry: Chebyshev 300 bands vs dense 1.4e-12 Ha.
QE 7.6 oracle runs: `tests/reference/response/r1/README.md` (dev08, pending).

## Consequences and limits

- Dense memory is about 3 N_PW^2 x 16 B per rank and every rank solves every
  k (no k or matrix distribution yet); the EVD is sequential faer.
- Chebyshev defaults are not calibrated on the generality set; its cost on
  the R1 cases is recorded but no speed target applies (D10).
- The unit gather is an allgatherv (each rank briefly holds a full unit).
