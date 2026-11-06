# Design note: real-wavefunction Gamma path (backlog B02, plan stage S4a)

Written 2026-11-06 by the performance owner as ticket input for Codex implementers. Documentation only: no code was written or run for this note. It is built on the analysis note `docs/research/qe-study/qe-gamma-real-path.md` (G) and `qe-hpsi-and-fft.md` (H), both in our own words, and on a read-only reading of the current BAANDI code (`dft-solvers/src/lib.rs`, `dft-physics/src/hamiltonian.rs` and `density.rs`, `dft-numerics/src/fused.rs`, `dft-basis/src/lib.rs`, `dft-engine/src/scf.rs`, `collinear.rs`). No QE source is used here and implementers need none. Labels: **M** = dev08 measured and relayed, **O** = Mac observation, **E** = estimate, **der.** = arithmetic on recorded numbers. Nothing here is a new measurement.

## 1. Goal, scope, and what is NOT claimed

Goal: at the Gamma point (one k point, k = 0, unshifted) run the eigensolver, the Hamiltonian application and the density accumulation on real wavefunctions, so that two bands share one complex FFT and all dense subspace algebra is real arithmetic on half-length data. QE's real Gamma path does exactly this and is 1.87x faster than QE's own complex path on si8 at one core (M, dev08 part 6: 0.366 s vs 0.684 s at conv_thr 1e-8; 1.3 to 1.9x across cases), with the h_psi call count nearly unchanged (94 vs 98). The relevant target is QE's best configuration, so this path is the largest remaining factor on Gamma cases (exceed-qe-plan.md, stage S4a: time factor 0.48, range 0.43 to 0.58, E).

Scope of the first release of this feature: unpolarized and collinear-spin SCF at Gamma, norm-conserving, real local potential (LDA, GGA, NLCC), Fermi-Dirac or insulating occupations. Out of scope: any k != 0, shifted meshes, spin-orbit, noncollinear magnetism, complex fields, other time-reversal-invariant k points (G-H8, speculative), locking, RMM and the spectral guard (these fall back to the complex path with a logged reason). The complex path stays the default and stays bit-identical; the real path is an isolated, switchable second solver family (a new ADR is required: section 14).

It is not a way to close the whole gap to QE: QE gains the same factor from the same trick, so this prevents the gap from widening (G section 9, item 5).

## 2. Why it is exact (one paragraph an implementer needs)

At k = 0 the operator H = T + V_loc(r) + V_NL commutes with complex conjugation of the periodic function when V_loc(r) is real and every projector obeys beta(-G) = conj(beta(G)). BAANDI's projectors do (real spherical harmonics, factor (-i)^l, structure factor exp(-iG.tau); `NonlocalOperator::new`). Then the complexified spectrum equals the spectrum of H restricted to the real-structured subspace {c : c(-G) = conj c(G)}, which has REAL dimension N_full (the full-sphere plane-wave count, odd at Gamma). Every eigenspace has a real-structured basis, so solving the real symmetric problem on that subspace gives the same eigenvalues, the same density and the same energies as the complex solver, with half the numbers per vector. The Davidson algorithm (blocks, restarts, tolerances, preconditioner) is unchanged; only the arithmetic differs.

## 3. Where it plugs in (code map)

| Layer | Existing | New (separate files; no edit of complex kernels) |
|---|---|---|
| `dft-basis` | `PlaneWaveBasis` (full sphere sorted by |G|^2 then Miller index), `FftGridSpec::g_to_fft` | `HalfSphere` map built from a Gamma basis: representatives, plus/minus positions in the full list, validity check |
| `dft-numerics` | `ComplexMatrix`, faer GEMM and `hermitian_eigh`, tiled deterministic GEMM, `FusedOperator` plans | `RealMatrix` (f64, column-major), real GEMM with the same fixed-tile ordered reduction, real symmetric eigen solve (faer `self_adjoint_evd` on f64), pair methods on the fused plan |
| `dft-physics` | `Hamiltonian`, `NonlocalOperator` (public `atoms()`, `columns()`, `coupling_hartree()`), `density::from_orbitals_fused` | `RealHamiltonian` (kinetic, paired local, real nonlocal), pair density accumulation |
| `dft-solvers` | `LinearOperator`, `Davidson` (complex) | `RealLinearOperator` trait, `RealDavidson` module with the same modes |
| `dft-engine` | `scf.rs` per-k solve (line about 963), `collinear.rs` per-channel solve | dispatcher, option `gamma_real = off | auto | on` (default `off`), real states in `ElectronicState`, expansion to full complex vectors at the output boundary |

Check `ownership.toml` and each crate's AGENTS.md for owned files in every ticket (section 12). The existing `NonlocalOperator` is read through its public accessors, so it is not edited.

## 4. State layout: orthonormal real coordinates (ORC), my recommendation

QE stores the half sphere as complex numbers and repairs every product with a rank-one G=0 correction (G section 1 item 3, section 6). I recommend a cleaner representation, derived here (QE does not do exactly this; G-H9 is the nearest relative and is smaller).

**ORC layout.** A real-structured state with full-sphere coefficients c(G) is stored as a real vector x of length N_full:

- x[0] = c(0) (real by symmetry);
- for each half-sphere representative number i = 1 .. n_half-1 (the pair {G_i, -G_i}): x[2i-1] = sqrt(2) Re c(G_i), x[2i] = sqrt(2) Im c(G_i).

Length 1 + 2 (n_half - 1) = N_full. Properties (all testable):

1. Isometry. Sum over the full sphere |c(G)|^2 = c(0)^2 + 2 sum_half |c|^2 = x.x. So the norm, every inner product <a|b> (it is real), every overlap and reduced-matrix entry are PLAIN real dot products of ORC vectors, with no factor 2 and no G=0 correction. The existing algebra (orthonormalise, project, Gram, Ritz rotation) carries over with `f64` in place of `Complex64`.
2. The invariant c(-G) = conj c(G) holds by construction: there is no stored imaginary G=0 coefficient that could drift (G risk 1 disappears), and no operation can leave the real-structured subspace.
3. Every operator that is real-structured is a real N_full x N_full symmetric matrix in these coordinates: kinetic is diagonal (value T(G_i) copied to both slots of i), the nonlocal term is P D P^T with a real P (section 7), the local potential acts through the pair FFT (section 6).
4. Memory: a vector is N_full doubles, i.e. half of the complex vector (N_full complex numbers). Davidson basis and H*basis storage halves.

**Compatibility rule.** Anything outside the real path sees full-sphere complex coefficients obtained by Hermitian completion: c(G_i) = (x[2i-1] + i x[2i]) / sqrt(2), c(-G_i) = conj c(G_i), c(0) = x[0], in the existing full-sphere ordering. Checkpoints, volumetric and band-density writers, forces, stress, the energy audit and any other consumer keep using the complex representation through this completion (cheap, O(N)). The checkpoint stores the completed complex orbitals; the reader of a complex checkpoint into a real run is the optional conversion in section 10.

**Alternative kept for reference** (not recommended as the primary layout): QE-style half-sphere complex arrays with the factor 2 and a rank-one correction after every product. It reuses the complex storage but puts a special case on every kernel and every unit test. ORC trades an O(N) scale (sqrt 2) inside pack and unpack (which are written anyway) for no special cases anywhere else.

## 5. Half sphere and mirror table (`HalfSphere`)

Input: the Gamma `PlaneWaveBasis` (momenta sorted by |q|^2 then Miller). Construction:

1. Build a map Miller -> full index. For every full index f find mirror(f) with Miller -m. If any mirror is missing the basis is not inversion symmetric and the dispatcher falls back (log reason). At k = 0 the sphere test uses |G|^2 computed from exactly negated Cartesian components, so +G and -G always pass together; a unit test asserts this for cubic, skewed, elongated and slab cells.
2. Representatives: G = 0 (index 0 in the sorted list; assert norm exactly 0) and the member of each pair whose Miller triple is lexicographically positive (first nonzero component > 0). Keep them in the existing full-list order, so the half list inherits the prefix-by-|G|^2 property used by the density sphere truncation and the free-electron cold start.
3. Store `plus: Vec<u32>` (index of representative i in the full list) and `minus: Vec<u32>` (index of its mirror; equal to plus for i = 0). Also an FFT-grid mirror check: `g_to_fft()[minus]` must be the grid position of -G (it is, because the plan is built from Miller triples; test it).

The fused FFT plan is built from the FULL-sphere index list exactly as today and is reused unchanged by the first implementation (pack into a full-sphere vector, call the existing routine, unpack). Pruning, plane fusion and the potential multiply fusion are therefore untouched (G section 8 item 11). Do not add real-to-complex transforms (same factor 2, new code path).

## 6. FFT pairing of bands

Let a, b be two ORC vectors (bands) and let ca, cb be their complex half-sphere coefficients (section 4). Pack into the full-sphere complex vector phi:

- phi[0] = x_a[0] + i x_b[0];
- for i >= 1: phi[plus_i] = ca + i cb, and phi[minus_i] = conj(ca) + i conj(cb).

phi is the full-sphere coefficient vector of psi_a(r) + i psi_b(r), both real functions. The existing local apply (inverse, multiply by real V(r), forward, divide by grid size) gives h = FT[V psi_a + i V psi_b]. Unpack:

- (V psi_a)(G_i) = (h[plus_i] + conj(h[minus_i])) / 2, (V psi_b)(G_i) = (h[plus_i] - conj(h[minus_i])) / (2i);
- at G = 0: V psi_a = Re h[0], V psi_b = Im h[0] (the imaginary leftovers, which are rounding noise, are discarded here, which is what keeps G=0 real);
- store into ORC with the sqrt(2) factor.

Cost: ceil(m/2) inverse plus ceil(m/2) forward transforms per block of m bands, i.e. half of today. Pack and unpack are about 2 N_full complex words of copy plus about 6 flops per G per pair (about 20 kflop on si8), negligible against about 8.6 Mflop of transform work per pair (G section 7.1).

**Odd number of bands.** The leftover band is packed with a zero partner (b = 0): the same pair routine is used, the extra transform work is the same as one complex transform today (no loss compared with the complex path, no special code path).

**Pairing rule (determinism).** Pairs are (0,1), (2,3), ... in the column order of the block passed to `apply_block`. The rule depends only on column order, never on thread count, pool state or which bands converged earlier. Because H is only applied to normalised vectors (Davidson applies H to orthonormalised corrections and to Ritz vectors; verify in the solver), the partner contamination is of the order of one rounding of the combined function, about 1e-16 relative to unit norm. Do not apply the paired routine to vectors of very different norm (a tiny residual next to a unit vector would lose relative accuracy); the real Davidson must keep this property and a debug assertion must check norms are within a factor 1e3 of each other.

**Density.** For each pair the inverse transform alone gives psi_a(r) = Re and psi_b(r) = Im; add w_a Re^2 + w_b Im^2 divided by the cell volume to the density grid. A pair is skipped when both weights are exactly zero (occupied bands are first, so only the trailing buffer pairs are skipped; the complex path skips each zero-weight band). Accumulation order: fixed pieces of pairs, each piece summing into its own grid buffer, buffers added in ascending piece order, so the result does not depend on the thread count. The downstream density (real grid array) is identical in type to today's, so mixing, Hartree, XC, symmetrisation (P19A, density only) and NLCC are untouched.

New methods needed on the plan (pure and generic `GridFusedPlan`, so the MKL build keeps working through the full-grid fallback): `apply_local_pair` (or the pack/call/unpack composition in `dft-physics`) and `accumulate_density_pair(phi, w_a, w_b, volume, density)`. The existing single-vector methods stay.

## 7. Nonlocal projectors in real arithmetic

The complex operator is sum over atoms and channels of beta_j D_jk beta_k^dagger with columns beta_j(G) (including 1/sqrt(volume) and the phase exp(-iG.tau)), coupling block D per atom (real, symmetric).

With ORC vectors, define for each column the real vector B_j: B_j[0] = Re beta_j(0), B_j[2i-1] = sqrt(2) Re beta_j(G_i), B_j[2i] = sqrt(2) Im beta_j(G_i), taken from the representative G_i = plus_i. Then:

- projection p_j = <beta_j|x> = sum over the full sphere conj(beta_j) c = B_j . x, a real number, because beta_j(-G) = conj beta_j(G);
- the back-application adds sum_jk B_j D_jk p_k to y (still real-structured).

So the nonlocal part for a block X (N_full x m) is two real GEMMs with the same real matrix P = [B_1 ... B_nkb] (N_full x nkb) and a tiny per-atom real product in between: proj = P^T X, w = D_blockdiag proj, Y += P w.

Facts to build and test:

1. P is built from the public complex columns of `NonlocalOperator` (no edit of that file). Only the representative entries are used. Validation at construction: |Im beta_j(0)| must be below 1e-14 times max|beta| (it is zero: the factor (-i)^l with odd l multiplies a radial transform that vanishes at |G| = 0, and the l = 0 and even-l values at G = 0 are real); if not, the dispatcher falls back.
2. The inner dimension is N_full (about 2970 for si8), real: 2 flops per multiply-add, a quarter of the 8-flop complex count (G section 7.3).
3. Storage is halved (8 bytes instead of 16 per element) and the block traversal reads P once per block per GEMM. The ledger declined a "per-block BLAS rewrite" of the COMPLEX nonlocal operator because it changes scalar association and A43 (row-blocked traversal) lost on dev08; this is a different change: a new real operator with no bit-identity contract against the complex one, and the reuse axis is the band block (see backlog B07 for the bit-identical complex alternative, which remains a separate ticket). State both in the ticket report.
4. Multi-projector fixtures (several radial projectors per l, d and f channels, species with different channel counts) must be in the tests: the per-atom coupling block is applied per atom, exactly as the complex path does.

## 8. Real Davidson (symmetric eigenproblem, real GEMM)

`RealDavidson` implements the same algorithm as the complex `Davidson` in `dft-solvers/src/lib.rs`, with the same parameters and meanings: block size = bands, capacity `max(davidson_subspace, 6 * bands)`, residuals in Ha, preconditioner correction -r/(d - lambda) with the sign-preserving floor 0.1 on the diagonal d, strict prefix plus buffer tolerance, certification by a true residual, and the standard-profile Ritz-value-change mode (ADR 0035). Everything below is the real image of the complex step:

1. Basis V (N_full x m) and stored products HV, ORC columns, orthonormal in the plain real dot product.
2. Reduced matrix T = V^T H V is real symmetric. Incremental update (the adopted P14B structure): keep the prefix, compute the rectangle between old and new columns and the new square block with real GEMMs; fill the other triangle by symmetry. In the real path the one-rectangle form is free (no conjugation, an exactly symmetric result), which was the open choice B01/V-5 in the complex path; the real path has no bit-identity contract with the complex one (G section 8 item 6).
3. Small problem: real symmetric eigen solve (faer `self_adjoint_evd` for f64, same library policy as today). The dense oracle for tests is the existing `DenseEigensolver` (dimension limit 1024) on complexified operators or a real dense equivalent.
4. Ritz vectors, Ritz products and residuals R = HV C - V C Lambda: real GEMMs with real coefficient matrices, N_full x m times m x q.
5. Orthonormalisation of corrections against the basis: blocked projection (two passes with normalisation between, A36/A37 rule kept: do not project twice before normalising) and MGS inside the block, all with real dots. CholeskyQR2 stays rejected (A52).
6. Restart by contracting stored basis and H products (A12), certification with one true application at the end.
7. Real GEMM kernels must follow the existing determinism design: fixed tiles (256 rows) with ordered partial sums, independent of the worker count, scheduled on the shared persistent pool (A35, A40, A42 dispatch threshold in real multiply-adds scaled so that the work per tile is comparable: re-derive, do not copy the complex 262,144 figure). Scoped OS threads per GEMM are forbidden (A33, A35).
8. Cold start (`free_electron_initial_guess`): the lowest-kinetic ORC unit vectors in the order of the half list: coordinate 0 (G = 0) and then the cos-type and sin-type vectors of the representatives with the smallest |G|^2. This is a real-structured start. It can differ from the complex start inside degenerate shells that are cut at the band count (the complex start takes single plane waves, which are not real-structured), so the first SCF iteration can have a different Davidson step count; compare totals, not per-iteration counts.
9. Modes in the first release: strict with buffer tolerances and certification, and Ritz-change (standard profile). Locking, spectral guard, RMM and hybrid are not implemented; the dispatcher falls back to the complex path when those options are on.
10. Convergence semantics are identical: the residual 2-norm of an ORC vector equals the complex residual norm (isometry), Ritz values are the same Ha numbers, so `eigen_tolerance`, buffer tolerance and the ratcheted standard thresholds are used unchanged and not retuned.

Iteration counts: for real-structured starts the complex Davidson IS the real Davidson up to rounding, so equal counts are expected; the only systematic difference is the start in item 8. This is a hypothesis to check (risk R2), not an assumption.

## 9. Interaction with SCF features

- **Standard profile (ADR 0035).** The threshold ladder, the Ritz-change stopping rule and the uncertified flag are properties of the solver mode and the SCF driver; the real solver implements the same mode (ticket R6) and reports `certified: false` identically. Strict path: certificate and spectral guard semantics as in section 8 item 9.
- **Fermi-Dirac / smearing (MF04).** Occupations enter only as density weights and in the energy; they are computed from the Ritz values, which are identical numbers in both paths. Partially occupied bands raise the number of bands to be paired in the density, no other change. Entropy term and Fermi level search untouched. Test with a metallic fixture.
- **Collinear spin.** Each spin channel has its own potential and its own Davidson solve (`collinear.rs`); the real path runs per channel and pairs bands only inside a channel (different potentials cannot share a transform). Staged after the unpolarized path (ticket R5b).
- **Symmetry (P19A).** Acts on the density grid after accumulation; at Gamma there is one k point and nothing to reduce. Unchanged.
- **Forces, stress, relaxation.** They consume orbitals through `state_view`; they receive the Hermitian completion (section 4) and need no change. Forces and stress keep the strict path (ADR 0035 forbids mixing with standard).
- **PBE/GGA, NLCC.** Only real-space potential and density; real-structured. No change.
- **Other operators later (DFT+U, hybrids).** Any operator must provide a real-structured ORC form or the dispatcher falls back.
- **Dispatcher.** `gamma_real = auto` selects the real path only if all hold: one k point equal to 0 after wrapping and an unshifted mesh, collinear or unpolarized, norm-conserving, no complex field, no locking/spectral guard/RMM, HalfSphere construction succeeds, projector G=0 check passes. Otherwise the complex path runs and the reason is logged once. `on` fails with an error instead of falling back (for tests). Default `off` until the gates of section 11 pass on dev08; then the user decides whether `auto` becomes the default.

## 10. Checkpoint, restart, outputs

The checkpoint keeps the full-sphere complex format (Hermitian completion on write), so readers, cube/XSF writers and band-density output are unchanged; the checkpoint records `gamma_real: true` in its metadata. Reading a complex-path checkpoint into a real run: first release cold-starts the orbitals (density restart works unchanged) and logs it. Optional later conversion: take Re psi(r) and Im psi(r) (each is real-structured), orthonormalise rank-revealingly and keep the lowest `bands`. Wavefunction phases and signs are gauge: tests compare gauge-invariant quantities only (section 11).

## 11. Equivalence tests against the complex Gamma path

**Tolerances and why.**

| Quantity | Gate | Rationale |
|---|---|---|
| Component: pair application vs two single complex applications | 1e-14 relative (normwise) | one transform rounding; random real-structured unit vectors |
| Component: real nonlocal projections and application vs complex | 1e-13 relative | 2970-term sums, a few ulp |
| Component: H_real x vs H_complex x on random real-structured vectors | 1e-13 relative, 10 cell shapes | sum of the above |
| Fixed-H eigenvalues vs the dense oracle | 1e-12 Ha | solver tolerance chosen tight in the test |
| Total energy, SCF at default strict settings, vs complex path | **1e-10 Ha** | the engine's own strict default stays within 5e-13 Ha of its tight reference on every fixture with a plateau (P17A-1: si8 4.97e-13, Si 3.4e-14, MgO 3.8e-13), so two valid trajectories to the same fixed point differ by about 1e-12; 1e-10 leaves two orders of margin for different trajectories and still rejects every symmetry or G=0 mistake, which show up at 1e-6 Ha or larger. The gate is not relaxed if a fixture fails: the cause is found. |
| Same, tightened SCF (density_tolerance 1e-10) | reported, expected below 1e-12 Ha | shows the floor; not a gate |
| Occupied eigenvalues, density L2 difference | gate fixed before the real run from the observed complex-vs-complex spread (different seed, T = 1 vs 8) and written into the ticket report; not changed afterwards | density_tolerance is 1e-7; the floor is set by the complex path's own noise, not guessed |
| SCF iterations | equal within one on realistic fixtures; toy fixtures judged by energy only (ledger C.11) | |
| Davidson steps and H vectors | reported per fixture; mean over the suite not more than 10 percent above the complex path (E), otherwise investigate before any default | iteration-count risk R2 |
| Determinism | energy and density hash bit-identical for T = 1, 2, 4, 8, 12, pure and MKL separately | existing contract |
| Option off | bench_scf `--require-identical` on gamma, k222s, k222, k444s, si8 unchanged | isolation |

Never compare wavefunction coefficients; compare eigenvalues, energies, density, and sums over degenerate subspaces (cubic Si has triple degeneracies). Bit-identity with the complex path is not a contract; density hashes change.

**Fixture set (generality rule; each shape listed, not Si alone).** Si2 (R1, stored energy -7.295530468210927 Ha in G; the suite has -7.295530468210772, a build difference of 1.5e-13), si8 (20 bands, buffer cluster), 4-atom elongated cell, skewed 20-atom synthetic, slab-like cell, 12-Bohr cubic cell, 16-atom Gamma case, MgO (two species), si-20 and mgo-40 for large band counts; synthetic projector fixtures with 0, 3 and 5 channels and several radial functions per l; an odd band count; a small-gap and a large-gap insulator; a metallic Fermi-Dirac fixture with fractional occupations; a collinear polarized fixture; a PBE fixture; NLCC. The three persistent suite non-convergers (si-elongated, syn-s3, syn-z14) are not required to converge but the behaviour (fail or pass) is reported for both paths.

**Component tests in the first tickets** (each against a straightforward complex reference on random real-structured vectors): ORC round trip and isometry (norm and dot equal the full-sphere complex values to 1e-15); pair pack/unpack round trip; pair vs two single applications; odd band; mirror table on four cell shapes; real nonlocal vs complex; real GEMM vs complex GEMM on real data; real Davidson on dense real matrices with clustered spectra; G=0 imaginary leftovers do not appear anywhere (they cannot be stored).

**QE evidence.** The validation reference for performance is a dev08 QE run with `K_POINTS gamma` (and the complex k = 0 run as second reference), recorded as a bundle with command, version, input and UPF SHA-256 (AGENTS.md rule; backlog M0 part 1 is relayed numbers only, the bundle is still missing). The old R1 to R8 QE energy comparisons stay valid as independent evidence; energies of the real path are compared with the complex BAANDI path first, then with QE with unchanged tolerances.

## 12. Ticket breakdown (each independent, testable, switchable, default off)

Order follows G (H1, H2 foundation, H3 independent, H5 safety) with ORC, and does not start before the dev08 profile R0 refreshes the fractions. Every ticket: `scripts/check.sh` exit 0, `scripts/check.sh <crate>`, new tests with stated tolerance, no change to complex kernels, bench_scf `--require-identical` for the complex path, ticket report with the sections of AGENT-BRIEF.md, and one ledger row from the orchestrator or me.

| Ticket | Content | Files (new, check ownership.toml) | Gate | Switch |
|---|---|---|---|---|
| R0 (no code) | dev08: QE `K_POINTS gamma` bundle; post-native, post-standard profile of si8 Gamma (at least 30 samples per named site, T = 1) giving the current FFT, GEMM-family, nonlocal and other shares; pair-FFT microbenchmark request. ADR "second solver family: real Gamma path" drafted by the architect | docs only | profile shows FFT + algebra + nonlocal at least 60 percent of wall (else expected gain is re-estimated before R2 onward) | n/a |
| R1 | `HalfSphere`, ORC pack/unpack (full complex <-> ORC), `RealMatrix`, real deterministic tiled GEMM, real symmetric eigen wrapper | `dft-basis` half.rs, `dft-numerics` real.rs | isometry 1e-15; GEMM vs complex-on-real 1e-14; mirror tests on cubic, skewed, elongated, slab | library only |
| R2 | Paired local application and paired density accumulation on the unchanged full-sphere plan (pure and generic fallback); microbenchmark binary (test-only) | `dft-numerics` fused_pair.rs (or additions behind new methods) | pair vs two singles 1e-14; microbenchmark at 25^3, 36^3, 48^3 at least 1.8x per band over two single applications (O first, dev08 decides) | library only |
| R3 | `RealNonlocal` (P from public columns, per-atom coupling), `RealHamiltonian::apply_block` (kinetic, paired local, nonlocal) | `dft-physics` real_gamma.rs | H_real vs H_complex 1e-13 on 10 shapes, multi-projector, odd count; timing per block vs current traversal on 4, 8, 20 atoms (O; 257, 2969, 16385 rows) | library only |
| R4 | `RealDavidson`: strict, buffer tolerances, certification, restart contraction, incremental reduced matrix; real GEMM scheduling on the shared pool | `dft-solvers` real_davidson.rs | fixed-H eigenvalues vs dense oracle 1e-12 Ha; iteration counts vs complex Davidson on fixed-H test problems (same operator, real start); determinism over thread budgets | library only |
| R5 | Engine integration (unpolarized): `RealStates`, dispatcher with logged fallback, density via pair accumulation, energy audit, Hermitian completion at the boundary, checkpoint metadata, option `gamma_real = off | auto | on` | `dft-engine` gamma_real.rs plus minimal hooks in scf.rs/session.rs (architect approves hook points) | section 11 system gates on the fixture set; option off bit-identical; `on` + fallback cases tested (k != 0, shifted mesh, locking on) | default `off` |
| R5b | Collinear spin and Fermi-Dirac through the real path | `dft-engine` collinear hooks | same gates on polarized and metallic fixtures | same switch |
| R6 | Standard-profile (Ritz-change) mode in `RealDavidson`; forces and stress checked through completion | `dft-solvers`, tests in `dft-engine` | frozen-threshold P17A cases: energy within 1e-10 Ha of the complex standard run, same uncertified flag; force and stress match the complex path to the existing identity-test tolerances | same switch |
| R7 | Pair work unit for the worker pool (pairs per piece, fixed piece size), deterministic density reduction | `dft-physics`, `dft-numerics` | bit-identical energy and density across T = 1, 2, 4, 8, 12; dev08 request: si8 T = 1, 4, 8, 16 pure and MKL, 3-run medians | same switch |
| R8 (optional, after R5 and dev08) | Half-sphere density and potential bookkeeping (G-H6), Hartree and mixing on N_half entries; MKL dgemm/dsyevd adapters for the real kernels; decision on `auto` as default | various | SCF counts unchanged, energy within 1e-12 Ha of R5 | same switch |

Each ticket may be launched alone with the stated tests. R1 to R4 add no engine behaviour. R5 is the first ticket that changes any reachable code path and it is guarded by the option.

## 13. Expected gain (all E unless stated; replace by the R0 profile)

Anchors: QE's own real/complex ratio 1.87x on si8 at one core, 1.3 to 1.9x across cases, with the same number of h_psi calls (M, relayed). G's derivation for BAANDI on si8 (2969 plane waves, 36^3 grid, 20 bands, 32 projector columns; G section 7, assumed time fractions from the ledger):

| Component | Mechanism | Factor (E) | Share of Gamma SCF assumed in G (E, to be replaced by R0) |
|---|---|---|---|
| FFT (local apply, density) | half the transforms; the pruned plan is unchanged | 1.95 to 2.0 | 0.35 |
| Dense products (reduced matrix, projections, Ritz, rotation) | 8 flops -> 2 per multiply-add, half the bytes | 3 to 4 | 0.40 |
| Nonlocal | real block GEMMs instead of per-vector loops | 4 to 8 (about 4x flops, plus reuse) | 0.065 |
| Rest (density FFT counted above, XC, Hartree, mixing, small solves, setup) | unchanged | 1.0 to 1.2 | 0.185 |

Whole SCF: about 2.1x (1.7 to 2.4), i.e. time factor 0.43 to 0.58 (exceed-qe-plan S4a). Staged: R2 plus R3 about 1.2 to 1.3x alone; R4 brings it to about 1.9x. Caveats: (1) after the standard profile (O: 3x fewer steps) and native flags (M: 1.3x) the share of Davidson algebra falls and the shares above change, so the whole-run factor must be re-read from the R0 profile; gains on the GEMM family shrink with its share. (2) Real kernels on dev08 may not reach the complex flop efficiency (the dev08 pure run gained only 1.06x from a 56 percent GEMM-flop cut, A55): the 3 to 4x on dense products is the least safe number. (3) Multi-core scaling is unproven (utilisation 23 percent, A48); R7 targets it, no number promised. (4) Where BAANDI can exceed QE and not only match: our pruned fused FFT per transform is ahead of QE's MKL-DFTI build (no pruning; H section 3), so the FFT term of the real path can be lower than QE's per pair; the other terms only reach parity.

Memory: Davidson basis and H products halve (about 5.7 MB complex to 2.9 MB real at the full capacity on si8, G 7.2), which helps large band counts where the basis leaves the cache (up to 2x bandwidth advantage, E).

## 14. Decisions that belong to the user

1. Approve the real Gamma path as a second, isolated solver family (new ADR, permanent maintenance cost) with default `off`; recommended yes (largest remaining Gamma factor, QE measured 1.3 to 1.9x).
2. After the dev08 gates: make `gamma_real = auto` the default for Gamma runs? Recommended yes only after R5/R6/R7 gates pass on dev08 and the results (energy, determinism, iteration counts) are shown; the complex path stays available.
3. Gate wording: this note uses energy 1e-10 Ha against the complex path (the research note G proposed 1e-12). Confirm 1e-10 as the fixed gate, with 1e-12 reported as an observation.
Not user decisions: ORC versus QE-style layout (technical, section 4), the dev08 measurement requests (orchestrator through the remote session).

## 15. Risks (ordered by how likely they are to cost time)

1. **Gain below the estimate on dev08.** Real kernels, memory traffic and Amdahl; A55 shows flop cuts do not map to pure-Rust time. Mitigation: R0 profile first, R2 microbenchmark gate, dev08 acceptance after each of R3, R4, R5.
2. **Different Davidson iteration counts** (start vectors, rounding, buffer cluster at the band cut). Mitigation: report counts per fixture, suite-mean gate, capacity retune only through B17 and not on silicon alone.
3. **Pair partner contamination** at rounding level when norms differ; mitigated by applying H only to normalised vectors, a debug assertion, and a flicker test on certificates with thresholds near 1e-14.
4. **Degenerate manifolds and gauge**: band-resolved outputs are not comparable coefficient by coefficient; tests use invariants only.
5. **Fixed-point multiplicity on toy fixtures** (ledger C.11): iteration counts and even energies can change; judge realistic fixtures by the gates and toys by energy tolerance and the dense oracle.
6. **Hidden complex-only assumptions** in consumers of orbitals (forces, stress, energy audit, checkpoint, band density). Mitigation: Hermitian completion at one boundary function with its own round-trip test; no consumer sees ORC.
7. **Projector G=0 check** fails for exotic UPF data (non-vanishing odd-l radial transform at zero): the dispatcher falls back; tested with a synthetic fixture.
8. **Scope creep**: k != 0 and other TRIM points are not part of this work (G-H8); locking and spectral guard are not supported in the real path in the first release and the dispatcher must say so when they are on.
9. **Second solver family** doubles the number of Davidson paths to keep correct; mitigated by isolation, shared test fixtures, and the complex path as permanent oracle. A later unification as one generic solver over a scalar type is possible but touches frozen-hash complex code and is not proposed here.

## 16. Ledger conformance

No rejected item is repeated: band-grouped FFT (A28) groups bands to reuse the potential; here two bands share one transform, a different mechanism with the factor 2 fixed by arithmetic. Lane staging (A11), threaded FFT (A34), scoped threads (A33, A35), row-blocked nonlocal traversal (A43), CholeskyQR2 (A52), capacity changes (A08, A13), hybrid RMM (A22) are not used. The incremental reduced matrix (A55) and restart contraction (A12) carry over. After any ticket completes, add ledger rows (R-prefixed IDs) with numbers, and update backlog B02.
