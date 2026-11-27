# ADR 0054: Response basis, chi0 and optical-limit conventions (R2)

- Status: **Proposed** (ticket RPA-R2-CHI0, 2026-11-27; orchestrator review pending).
- Scope: new crate `dft-response`, `dft_engine::response`,
  `dft_physics::stress::projector_k_gradient`, `NscfResult::{bases, grid_dims}`,
  `baandi dielectric`. Plan: `docs/roadmap/rpa-crpa-plan.md` 2.1-2.4, 3.1,
  3.4, 3.6 and ticket R2 (user decisions D1-D10 of 2026-11-26, option A each).

## Decision

1. Response basis: `{G : |q+G|^2/2 <= E_cut^chi}` from
   `PlaneWaveBasis::at_k(cell, E_cut^chi, q)` (q wrapped to [-1/2, 1/2)), so
   G = 0 is row 0 at q = 0. `E_cut^chi <= 4 E_cut` (density sphere) is
   enforced; with the density grid every pair density is alias-free.
2. Pair density `rho_nm(q+G) = <psi_nk| e^{-i(q+G).r} |psi_m,k+q>
   = sum_G1 conj(c_nk(G1)) c_m,k+q(G1+G)` (cell-normalized orbitals). An
   umklapp `k + q = k' + G0` is applied by placing row M of k' at grid point
   M - G0. Formed on full real-space grids with one-rank batched stick FFTs
   (`DistributedPairFft::new_direct`), forward with 1/N.
3. chi0 on the imaginary axis (sum over states, D1):
   `chi0_GG'(q, i w) = sum_k s_k sum_{n occ, m emp} rho rho'^* 2 D / (D^2 + w^2)`,
   `D = e_nk - e_m,k+q < 0`, `s_k = 2 w_k / Omega`; valid for a
   time-reversal-symmetric k set (checked by the driver). Each frequency
   accumulates `-(Q Q^H)`, `Q = P diag(sqrt(-f))`: Hermitian negative
   semidefinite by construction.
4. Optical limit: `rho_nm(q) = q . a_nm`, `a_nm = <u_n| dH_k/dk |u_m> / (e_m - e_n)`
   with `dH_k/dk = (k+G) + sum_atom (dB D B^H + B D dB^H)`; head
   `chi0_00 = q^2 q^.L.q^`, wings `chi0_G0 = q q^.W_G`, body G = 0 row and
   column zero. `projector_k_gradient` returns the projector table of
   `NonlocalOperator::new` (same channel order and phase) and the k
   gradient of its radial-angular factor (quadrature slopes, no
   interpolation); the phase gradient cancels in every atom-diagonal product.
5. Dielectric: symmetrized `eps~ = 1 - v^1/2 chi0 v^1/2`,
   `v^1/2(p) = sqrt(4 pi)/|p|`; optical tensor
   `E = 1 - 4 pi L - (4 pi)^2 W^H [eps~_BB]^-1 W / (|G||G'|)` (local fields,
   Schur complement over G != 0) and `1 - 4 pi L` (no local fields);
   `eps_M(q^) = q^.E.q^`. Reported: real symmetric parts and the dropped
   asymmetry. Finite q: `eps_M = 1/[eps~^-1]_00`.
6. Frequencies (D5): w = 0 plus `n` Gauss-Legendre nodes on t in (0,1),
   `w = w0 t/(1-t)`; `w0` = mean of the lowest direct transition and the
   largest included transition unless given.
7. MPI level 1 (D7): NSCF in the owned-row layout; per k the eigenvectors
   are gathered (permutation, identical bits on all ranks); units
   (k, block of empty bands) are dealt round robin in k order; each rank
   accumulates a replicated partial chi0 for all frequencies; one
   `tree_sum_f64` combines them. One rank is bit-identical to the
   threads-only run; frequencies are split over threads (each frequency
   sequential), so results do not depend on the thread count.
8. Guards: metals (indirect gap <= 1e-6 Ha, a numerical zero) reject (D6);
   spin-polarized densities reject (R3); non-time-reversal-symmetric k sets
   reject.

## Validation (internal oracles; `crates/dft-response/tests/oracles.rs`)

Pair FFT vs convolution with umklapp; (a) free-electron gas vs the exact
discrete sum (1e-12) and the analytic Lindhard function; (b) static columns
vs supercell finite differences (1e-6 gate); (c) Hermitian NSD; (d)
high-frequency closure over every empty state; (e) k.p head and wings with
the nonlocal velocity vs finite-q extrapolation; projector gradient vs
finite differences of `NonlocalOperator` rows; band-count monotonicity;
rank-reduction contract. Engine: `crates/dft-engine/tests/response_r2.rs`.
External oracles (Yambo, BerkeleyGW, QE ph.x) are dev08 items:
`tests/reference/response/r2/README.md`.

## Consequences

No checkpoint of chi0 units yet (`chi0_q_v1` later); q -> 0 only in the CLI;
the finite-q path exists in the library (tests) for R3. Chebyshev bands may
be selected; defaults unchanged.
