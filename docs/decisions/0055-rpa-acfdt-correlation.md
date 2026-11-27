# ADR 0055: RPA correlation energy (ACFDT) conventions and per-q checkpoint (R3)

- Status: **Accepted** 2026-11-28 (ticket RPA-R3-ACFDT, proposed 2026-11-27). Evidence: finite-q E_c vs ABINIT 10.8.3 on dev08 within 0.27 meV/atom (gate 1), q = 0 convention difference falling as 1/N_q, `tests/reference/response/r3/dev08-2026-11-28/summary.json`.
- Scope: `dft_response::acfdt`, `dft_response::transitions::Transitions::block`,
  `dft_numerics::linalg::hermitian_cholesky`, `dft_engine::response::rpa`,
  `baandi rpa` / `[task] type = "rpa"` / `[rpa]`. Plan:
  `docs/roadmap/rpa-crpa-plan.md` 2.1, 2.3-2.7, 3.5, 3.6 and ticket R3 (user
  decisions D1-D10 of 2026-11-26, option A each). Builds on ADR 0054.

## Decision

1. Energy. `E_c = sum_q w_q (1/2pi) int_0^inf dw Tr[ln(1 - X) + X]`,
   `X = v^1/2 chi0 v^1/2` (ADR 0054 conventions, `v^1/2(p) = sqrt(4 pi)/|p|`),
   Hartree per cell, q weights summing to one. X is Hermitized (the dropped
   defect is reported); `ln det(1 - X)` comes from the Cholesky factor of
   `1 - X` (`hermitian_cholesky`, faer, sequential). A failed factorization
   is an error (chi0 lost its sign), never clipped.
2. q sampling. q runs over the Gamma-centred k mesh; NSCF on the FULL mesh
   (no wavefunction rotation, plan 2.6 stage 1). `k + q = k' + G0` with
   G0 from the wrapped Bloch vectors of the bases; the k' rows are placed at
   M - G0 (ADR 0054 item 2). q reduction: `reduce_mesh` with time reversal
   and either the identity (`time-reversal`) or the detected crystal group
   (`crystal`, default; only mesh-preserving operations, exact orbits),
   or none. E_c(q) is invariant under these operations (unitary similarity
   of chi0 with a permutation of the response sphere).
   G-sphere mode (`[rpa] response_sphere = "G"`, one Miller set
   `|G|^2/2 <= E_cut^chi` for every q; amendment R3-QIMAGE, 2026-11-27).
   There the set {q + G} is not invariant under q -> q - G0, so E_c(q) is a
   function of the q VECTOR. Before the amendment the engine evaluated the
   [-1/2, 1/2) reduced-coordinate wrap of the mesh value, which is not
   rotation covariant (on fcc it is often a longer image than the first-zone
   one) and is not what the record printed (the [0, 1) mesh value), so
   symmetry stars split (Si dojo 4x4x4, 6 Ha, dev08: star (1/2,1/4,0)
   -0.422000 / -0.421084 / -0.420111 Ha) and crystal differed from none by
   1.1 meV/atom. Convention now: each mesh q is evaluated at its
   first-Brillouin-zone image, the shortest q - G0 in the cell metric
   (`first_zone_images`, ties at relative 1e-10 in |q|^2); a mesh q on a
   zone face, edge or corner with t tied images contributes the mean over
   them (t units of weight 1/(N_q t)). This definition is invariant under
   every lattice rotation and time reversal, so the reductions group units
   whose image vectors are related by S or -S without G0 (mesh-preserving
   subgroup of `reduce_mesh`), and none, time reversal and crystal agree
   to rounding. The stored q of a G-mode record is the evaluated vector.
   The q+G mode (default) is unchanged bit for bit.
3. Alias guard at finite q: product components `P = (M - G0) - G1` and the
   response sphere must not differ by a nonzero grid-superlattice vector; a
   candidate n is excluded by the Miller box ranges or by
   `|n N . B| > r_left + r_right + r_response`. Otherwise the run rejects
   (larger grid needed).
4. q = 0. `X(d)` is finite and direction dependent: `X_00 = 4 pi d.L.d`,
   `X_G0 = sqrt(4 pi) v^1/2(G) d.W_G`. Schur complement over the body
   `B = 1 - X_BB`: `Tr[ln(1-X(d)) + X(d)] = ln det B + ln eps_M^LF(d) + Tr X_BB + 4 pi d.L.d`
   with `eps_M^LF(d) = 1 - 4 pi d.L.d - d.Q.d`, `Q = U^H B^-1 U`,
   `U_Gb = sqrt(4 pi) v^1/2(G) W_Gb` (one Cholesky per frequency; every
   direction costs a 3x3 form). The q = 0 mesh term is the weighted average
   of this function over a Lebedev set (6, 14, 26 default, 50; exactness
   checked against sphere moments); every run also records the sets 6-50
   and a 24 x 48 product rule (the direction-convergence check) and the
   body-only value (G = 0 dropped, diagnostic). The |q| dependence inside
   the mini zone is not integrated: an O(1/N_q) error measured by the k-mesh
   series (plan 2.3).
5. Cutoff series. One pass over the pair densities of the largest sphere
   feeds one accumulator per cutoff (rows gathered by Miller index; the
   spheres are nested); band blocks are split at every cutoff's band count,
   so cutoff c sees exactly bands `< N_b(c)`. `bands = auto` ties
   `N_b(c) = N_G(c, q = 0)` (clamped to the smallest k basis, plan 2.5); an
   explicit count applies to every cutoff. Extrapolation: least squares
   `E(c) = E_inf + A c^{-3/2}` with residuals and the plan-2.5 linearity
   check (last three residuals <= 0.5 meV/atom); the same energies fitted in
   `1/N_G(q = 0)` are reported as a diagnostic (discrete shells).
6. Frequencies: the R2 mapped Gauss-Legendre grid (no w = 0 node needed);
   `omega0 = auto` as in ADR 0054 item 6, over every occupied channel.
7. Spin: one density (factor 2) or two collinear densities with integer
   electron counts per channel (from the density integral, 1e-6
   tolerance); `chi0 = chi0_up + chi0_down` with weight `w_k / Omega`; each
   occupied channel must be gapped; an empty channel (fully polarized
   moment) contributes nothing. The NSCF accepts collinear densities only
   with lda-pw or pbe.
8. MPI level 1 (D7): units (q, channel, k, band block) are dealt round robin
   in that order with one global counter; per q every cutoff accumulator is
   reduced by `tree_sum_f64`; the log det runs replicated on identical bits.
   One rank is bit-identical to the threads-only run; frequencies are split
   over threads (each frequency sequential).
9. Checkpoint per q (`rpa_q_v1`, CLI): one JSON record per irreducible q
   (`rpa_q_<index>.json`, written by rank 0 through a temporary file and
   rename) with every float as IEEE-754 bits plus readable copies, and a
   fingerprint `sha256("rpa_q_v1|" + input SHA-256 + "|" + SCF density SHA-256)`.
   The engine re-checks q, weight, optical flag, omega0 bits, cutoff count
   and frequency count of a loaded record; any mismatch rejects. chi0 is
   never stored (plan 3.5: recompute, never store pair densities; a
   `chi0_q_v1` payload is not needed for small systems). The NSCF runs on
   resume (its own `nscf_bands_v1` store is not wired here).

## Validation (internal; no external number yet)

- Analytic: the plasmon formula `E_c = (1/2)(sum Omega - sum Delta - Tr K)`,
  `Omega^2 = eig(D^2 + 2 D^1/2 K D^1/2)`, `K = s R^H v R`, for random
  multi-transition models through the accumulator, Cholesky and the
  frequency grid (`crates/dft-response/tests/acfdt.rs`): 32 points
  within 7e-13 Ha (relative <= 5e-13).
- Cholesky log det vs eigenvalue path (`finite_q_term_eigen`) <= 1e-12
  relative; Schur q = 0 form vs the direct determinant of the full X(d)
  <= 1e-12; Lebedev exactness <= 1e-14.
- Engine (`crates/dft-engine/tests/response_r3.rs`, synthetic skewed cell):
  np = 1 bit-identical, np = 2, 3 within 2e-15 Ha; q reduction (none vs
  time reversal) <= 1.4e-14 Ha; frequency 12 -> 24 points 4e-7 meV/atom;
  Lebedev-26 vs product rule 1e-10 meV; band-count monotonic; spin channels
  of an unpolarized density equal the unpolarized path to 7e-18 Ha;
  checkpoint resume bit-identical, stale record rejects.
- CLI (`crates/dft-cli/tests/rpa_r3.rs`): CSV views equal the JSON, the
  q-weighted sum and the frequency integral reproduce E_c, resume, single
  entry, guards.
- Si.pz-vbc observations on the Mac: `tests/reference/response/r3/mac-2026-11-27/`.
- External oracle (ABINIT `gwrpacorr`, plan R3 gates 1 meV/atom matched and
  2 meV/atom extrapolated): `tests/reference/response/r3/README.md`, dev08.

## Consequences

No EXX (R3b), no IBZ wavefunction rotation, no Sternheimer path; the q = 0
term is the direction-averaged limit only. Memory: one accumulator per
cutoff (N_w N_G(c)^2 complex each) plus thread workspaces. Metals still
reject (D6).
