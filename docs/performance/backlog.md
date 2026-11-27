# Performance backlog: ideas not yet tried, or only partly tried

Revised 2026-11-05 after the four QE study notes of the same date (G = `qe-gamma-real-path.md`, H = `qe-hpsi-and-fft.md`, S = `qe-scf-cycle.md`, V = `qe-solvers-and-subspace.md`; numbering inside a note as "H1" means that note's hypothesis 1, written G-H1, H-3, S-9, V-2). Ranked by (justified gain x chance it works) / risk. The goal behind the ranking is [exceed-qe-plan.md](exceed-qe-plan.md): match, then exceed QE 7.6 at equal cores and equal accuracy against QE's best configuration. Ledger IDs (`A..`) refer to [attempts-ledger.md](attempts-ledger.md). Gains are ESTIMATES or CEILINGS unless stated measured; the basis is given. Acceptance is always dev08 via the orchestrator ([AGENT-BRIEF.md](AGENT-BRIEF.md)). Labels: **ready**, **needs-measurement**, **blocked** (needs a user decision, a case, or a prerequisite).

## What the notes changed (read first)
1. The gap to QE is NOT the Hamiltonian application: QE h_psi costs about 3.7 ms per Davidson step, BAANDI about 3.5-4.5 ms (H). The gap is Davidson algebra around H (all-root Ritz and residual work, orthonormalisation, large basis) and the H vector count (about 2.2x). By the V cost model one BAANDI step costs about 5600N multiply-adds against about 1100N for QE and reproduces the P14B explicit GEMM total within 20% (1.1e10 vs 9.49e9 in 2mnk). Old ranks that treated FFT/nonlocal as the main lever are lowered; Davidson algebra is raised to rank 2.
2. Davidson capacity is `max(48, 6*bands)` (`crates/dft-engine/src/scf.rs`): 120 for the 20-band si8 case, not 48. The mean basis is about 70 (QE about 30). Capacity experiments A08/A13 were run with 8 bands (capacity 48) and say nothing about 20 bands with all-root algebra.
3. QE's convergence definition at conv_thr 1e-8 is weaker (Ritz-value change test, up to 5 unconverged roots accepted in SCF, no certificate, no gap check). Part of the step-count gap is definition, not speed (S section 5.4-5.5, V 2.5). A fair comparison needs an explicit equal-accuracy definition (plan section b).
4. QE's MKL (DFTI) FFT path does not prune (H 1.9): our pruned pipeline is already ahead of it per transform; an FFTW3-built QE prunes, so QE's best configuration may be an FFTW3 build (to be asked of the remote session).
5. The repository has no `[profile.release]` and no target-cpu flags (H section 4d): hand-written loops are compiled for baseline x86-64 on dev08. A build-flag experiment on dev08 is requested.
6. QE's real-Gamma path (`K_POINTS gamma`) is the correct reference for Gamma-only runs and is now MEASURED (dev08 part 6): si8 N=1 0.366 s at 1e-8 (0.500 s at 1e-12) vs 0.684 s complex, 1.87x (1.3-1.9x across cases), h_psi calls 94 vs 98 (the saving is per call, not fewer calls), QE algebra only about 0.05 s. The real target on si8 is BAANDI/QE = 13-15, not 7-8. The real Gamma path is therefore the largest single factor for Gamma cases. QE and BAANDI on dev08 are both default builds (QE -O3 plus MKL dispatch; no RUSTFLAGS for BAANDI); the build-flag experiment is being run by the remote session. The P14B "safe removable bound is zero" holds only for unlocked roots (V 7).
7. dev08 pure si8 T=1 gained only 1.06x from a 56% cut in explicit GEMM flops (A55), so the flop model under-predicts time effects on pure Rust; a post-P14B dev08 profile is the first step before any algebra work.

## Update 2026-11-06 (after P17A, P19A, P24A, P15A/C, native flags; labels M = dev08 measured, O = Mac observation, E = estimate)
- **Done or decided:** native flags measured (M: pure 1.12-1.40x, bit-identical; pure native beats default MKL on si8); standard profile P17A implemented (O: si8 2.51 -> 0.80 s, 345 -> 69 steps, 2265 -> 742 H vectors, dE 1.0e-9 Ha; uncertified; dev08 equal-accuracy run pending); symmetry P19A implemented, default off (O: k444s 32 -> 10 k, 3.33 -> 1.11 s, dE 1.2e-14 Ha; QE symmetry-on reference pending); P15A replay (stopping rule explains most of the step-count gap: 283 -> 39 rows at identical tolerances, QE-ladder 89 rows vs QE 83); P15C spectral guard optional (+4% si8, +23% k444s H columns); forces, stress, relaxation, cell relaxation exist.
- **B01 locking and one rectangle: tried (P24A), not adopted.** O: complex multiply-adds -6.0% si8, -12.5% k444s (model predicted 30-40%), H applications +24% / +21%, failed gates (syn-s2 4.4e-3 Ha; si8 buffer roots 3.7e-9 Ha; skew diverges); option `locking = "auto"`, default off. Lowered to option; do not repeat without a new mechanism (ledger P24A-1).
- **B13 non-orthonormal basis: priority lowered.** The serial scalar MGS it targeted shrank by itself under native flags (`subtract_projection` 0.822 -> 0.224 s) and the standard profile cut steps 5x; revised estimate 5-12% of wall (E), high risk on skew fixtures; park unless a dev08 profile after native + standard shows algebra above 40% of wall.
- **B04 and B14 resolved in part:** the profile exists (P17A); remaining: dev08 equal-accuracy comparison, repair of skew/slab/mesh16 under standard, Hartree-metric mixing not adopted (MgO energy gate).
- **B02 real Gamma path: rank 2, start now.** QE gets 1.3-1.9x (1.87x on si8) from it with the h_psi call count unchanged; after standard + native the H application and phases are a larger share, so it is the largest remaining Gamma factor (E: time factor 0.43-0.58).
- **B07 band-blocked nonlocal:** ticket P15B exists (not merged when written); `NonlocalOperator::add` is unchanged by native flags (0.6 s of 4.1 s, directional, below 30 samples).
- **B08 MKL pruned FFT:** lowered; pure native beats default MKL on si8.
- **B15 symmetry:** implemented; remaining: dev08 and a QE symmetry-on bundle, default decision (plan section d).
- **B19 new: forces, stress and relaxation cost.** Measured shares (O): forces 28.3% of SCF+force on si8, 10.9% on k444s; stress 58.4% / 44.6% of SCF+pass (Pulay cutoff scans extra). The standard profile does not apply to forces (uncertified), so after SCF gets faster these phases dominate force/stress/relaxation runs. First experiment: profile the force and stress phases on si8 and k444s on dev08 and compare with QE `forces`/`tstress` timings; add the workloads to the benchmark set (plan section b).
- Re-ranking after this update: 1 M0 (still open: QE symmetry-on, equal-accuracy run, post-native profile), 2 B02 real Gamma, 3 B04/B14 follow-ups and equal-accuracy evidence, 4 B07 band-blocked nonlocal, 5 B15 default decision plus k-point work, 6 B19 forces/stress, 7 B03 parallel decomposition (needed for multi-core), then B10, B05, B16, B09, B17; B01 (option), B13, B08, B11 lowered.

## Ranking (as of 2026-11-05; see the update above for 2026-11-06)

| Rank | ID | Idea (source notes) | Label | Gain basis |
|---|---|---|---|---|
| 1 | M0 | dev08 measurement package: QE real Gamma, QE variants (FFTW3 build, nbnd, ndim, symmetry), BAANDI post-P14B profile, build-flag experiment, E0/E1/E2 and frozen replay (G-H0, V-1, S-E0..E2, H-4) | ready | decides every rank below |
| 2 | B01 | Davidson algebra: hard locking with margin and cluster rule; one-rectangle reduced H (V-2, V-5; ledger B01) | needs-measurement | V model: 30-40% of subspace multiply-adds; ceiling 56.8% of remaining si8 GEMM 2mnk |
| 3 | B13 | Non-orthonormal basis with conditioned small generalised problem (V-3, QE-style) | needs-measurement | V model: 25-30% of subspace work and the serial MGS (4-26% of dev08 T=1 wall) |
| 4 | B04 | Step-count reduction: Hartree-metric error E_res, energy-weighted ratcheted schedule with first-iteration redo, replay first (S-2, S-1, V-1; old H6) | needs-measurement | ceiling 3.49x steps; 15-40% if a quarter to half is closed (S-9.2) |
| 5 | B14 | Additive "standard" convergence profile and equal-accuracy harness (S-6, S-E1) | blocked (user decision) | definition share unknown; QE moving conv_thr 1e-8 to 1e-12 costs 33% |
| 6 | B07 | Cheap bit-identical H-side items: band-blocked nonlocal, Hartree/kernel caches, QE-style cheap energy in non-final iterations, release flags (H-1, H-4, S-7, S-8) | ready | nonlocal 5-6% (H-1), flags 0-10% of non-library part (H-4), audit 3-5% (S-7) |
| 7 | B02 | Real-valued Gamma path (G-H1..H6; two bands per FFT, half sphere, real GEMM, real nonlocal) | needs-measurement (after rank 2/3 structure); raised in value | QE gains 1.87x on si8 N=1 from it (measured), so BAANDI must follow to stop falling behind; G 7.4 estimate 2.1x (1.7-2.4) for BAANDI |
| 8 | B03 | Band-level decomposition, fewer larger parallel regions; pair work unit in real path (ledger B03, G-H4) | needs-measurement | utilisation 23% vs 33% for 2x |
| 9 | B15 | k-point symmetry reduction and density symmetrisation (S-9, S section 3.1) | blocked (user decision, ADR) | QE default k444s about 10 irreducible k vs 32 (3.2x, estimate, standard fcc result not recorded here) |
| 10 | B08 | MKL pruned stick FFT pipeline; matrix-product nonlocal (H-3, H-2) | needs-measurement | per vector FFT 0.3-0.4 -> 0.2-0.25 ms (estimate); nonlocal 0.2 -> 0.05 ms |
| 11 | B10 | Starting density and wavefunctions (atomic or free-electron subspace) (S-3, S-4, V-8, S-E6..E8) | needs-measurement | first SCF iteration 4.3% (si8) to 12.6% (k444s) of H vectors; 3-5% (S-9.4) |
| 12 | B05 | Thick restart with previous directions (V-4); tail-root acceleration | needs-measurement | tail iterations hold 55% of si8 steps |
| 13 | B16 | Hartree-metric mixing metric (S-1, S-E4), mixer resets survey (S-5, S-E9) | needs-measurement | 0-3 SCF iterations (one is 6.7% of si8) |
| 14 | B06 | Serial remainder at Gamma; parallel audit/potential; EVD routing (ledger B06) | ready (audit), needs-measurement (EVD) | 5-10% of si8 T=8 wall |
| 15 | B09 | k-point parallel efficiency (ledger B09) | needs-measurement | 7.1x of 16 on k444s |
| 16 | B17 | Capacity and restart retune after ranks 2-3 (V-6) | needs-measurement | unknown |
| 17 | B11 | RMM single-sweep hybrid as switchable module (V-7); converged RMM at >=128 bands | blocked | A22 negative; regime untested |
| 18 | B18 | Low priority or exploratory: lane staging re-test on dev08 (H-5), radial table/shell cache (H-6), truncated-potential box (H-7), stick-major order (H-8), scaled G=0 (G-H9), other TRIM points (G-H8) | needs-measurement | each below 5% or unproven |
| 19 | B12 | MPI, GPU | blocked | out of scope (AGENTS.md) |

Common gates (details in the brief): energy within 1e-12 Ha of the stored baseline on gamma, k222s, k222, k444s, si8 and the generality set (cubic, skewed, elongated, slab; 4/8/20 atoms; synthetic projectors); SCF count and certificate unchanged for realistic fixtures (P14B rule for toys); density bit-identical across thread counts within a backend; no tolerance relaxed; no new convergence failure.

## Items

### M0 dev08 measurement package (no code; ready)
- Content: (1) DONE (part 6): QE `K_POINTS gamma`, si8 0.366 s / 0.500 s, gamma 0.070 s / 0.084 s at N=1; remaining: store as a reference bundle with command and UPF SHA-256 (G-H0); (2) QE variants: FFTW3 build if one exists on dev08 (H 1.8), `nbnd` default (occupied only), `diago_david_ndim` 2 and 4, symmetry on for k444s and k222, pools versus plane-wave MPI at N=8 and 16; (3) BAANDI post-P14B `perf` profiles, at least 30 samples on each named site, si8 T=1 and T=8, pure and MKL; (4) build-flag experiment: `lto="fat"`, `codegen-units=1`, `-C target-cpu=x86-64-v3` (H-4), same protocol, `--require-identical`; (5) accuracy: E0 (E_res, Hartree metric, stored every SCF iteration), E1 (BAANDI energy error and iteration when E_res < 5e-9 Ha), E2 (offline replay of stopping rules), frozen-Hamiltonian replay (V-1).
- Gain: none directly; decides ranks 2-5 and replaces estimates by numbers. Risk: remote time; QE timer 0.01 s (use cases above 0.2 s).
- Builds on / avoids: ledger C.7, C.8; avoids mixing denominators (C.14).

### B01 Hard locking and one-rectangle reduced H (V-2, V-5; my B01)
- Gain: V model: B term (Ritz vectors and products for all roots, 2 q m) is half of a step and dominates tail steps (k=2, q=70: B=2800 vs A+C about 600); locking removes most of B, shrinks q in A, C and the EVD (cubic). Estimate 30-40% of subspace multiply-adds, 40-60% of EVD time, 10-20% of si8 T=1 wall (more at T=8). Ceiling from P14B: Ritz/residual 5386241040 of 9489226838 2mnk (56.8%). One rectangle (V-5): 6-8% of subspace multiply-adds, nearly free, reverses the explicit P14B "no conjugation fill" choice.
- Risk: a locked vector degraded by cluster rotation (lock lowest roots first, whole clusters, margin c about sqrt(m), unlock if the certificate fails); trajectory change (rounding policy needed); dev08 pure showed weak flop-to-time coupling (A55).
- Prerequisites: M0 profile; re-crossing count from the P14A collector (does a root fall below tolerance and later exceed it?).
- Validate: gates; certificate on all roots; H vectors not above baseline by more than 2%; P14B counters show at least 25% fewer subspace multiply-adds on si8 and k444s; skewed, elongated, 12-Bohr cubic, 19-case suite.
- Builds on A55, A56, A53; avoids A15 (strict refinement) and A52 (CholeskyQR2). Conflicts: B13 (do locking first, then non-orthonormal, never together; V-3).
- First experiment: the re-crossing count (already planned in B01 of the previous backlog) plus the V model per-step accounting with the P14B site counters on si8 tail steps.

### B13 Non-orthonormal Davidson (V-3, the old H4)
- Gain: removes term C (about one third of multiply-adds plus the sequential scalar MGS, 4-26% of dev08 T=1 wall); adds about 5% for the Gram rectangle; net 25-30% of subspace work (V estimate). QE structure: normalise corrections only, solve `H c = e S c` by canonical orthogonalisation with a Cholesky fast path.
- Risk: accuracy eps*cond(S_red); nearly dependent corrections are exactly the 20-atom skewed case where CholeskyQR2 regressed (A52): mandatory stage 0: log retained rank and smallest Gram eigenvalue under MGS and CholeskyQR2 to learn whether that regression was a rank decision or rounding. QE keeps accuracy only through frequent restarts and a weak stopping rule (V 1.3).
- Gates: as B01 plus min eigenvalue of S_red logged, no cond above 1e8 without fallback.
- Builds on: A16 (gated), A52. Conflicts: B01 (order), B17.
- First experiment: the stage-0 rank/conditioning diagnostic on the skewed fixture.

### B04 Step-count reduction (S-2, S-1, V-1; old H6, B04)
- Gain: ceiling 3.49x steps (290 vs 83); S estimates 15-40% of Davidson steps if a quarter to half of the gap closes; part of the gap is definition (see B14) and part QE start quality (iterations 1-2: 3-4 steps vs 9-10). BAANDI residual tolerance may be 100-1000x tighter than QE's late Ritz-change threshold (V 6.3, estimate), worth 4-8 steps per late SCF iteration; but P14A (A21) shows uniform loosening raises SCF count, so QE's looseness is tolerated by its mixer and not by ours (unexplained).
- Risk: high for the schedule; none for the replay.
- Prerequisites: frozen-Hamiltonian replay (V-1) separating criterion strength from convergence rate; E0, E2, E3 (calibrate kappa), E5 (ratchet).
- Validate: Davidson steps lower AND SCF count equal on at least 8 systems of different shape, gap, band count and k mesh; certificate unchanged.
- Builds on A21 (different family: metric, ratchet, gap mapping tau = sqrt(g epsilon/N_e)); avoids A15 (no strict refinement) and A14.
- First experiment (offline): E0/E2 replay from stored histories before any solver change.

### B14 Standard convergence profile and equal-accuracy harness (S-6; user decision)
- Additive, labelled profile: convergence by Hartree-metric E_res threshold, energy from output density with first-order correction, certificate only on request; strict default unchanged. Without it the comparison with QE stays definition-limited. See plan sections b and d.

### B07 Cheap H-side items (H-1, H-4, S-7, S-8; was B06/B07 partly)
- Band-blocked nonlocal (H-1): process g bands per projector column, bit-identical per band, removes L3 streaming (3.04 MB per band) and the latency-bound dependent chain; nonlocal 0.2 -> 0.05-0.08 ms per vector, about 5-6% of T=1 wall on dev08 (estimate); judge on dev08 only; must cite A43 and show the matrix-reuse effect separately from row blocking.
- Release/target flags (H-4): 0-10% of the non-library part, unknown; native artifact or run-time dispatch, not a silent floor raise.
- QE-style cheap energy in non-final iterations (S-7): 3-5% on Gamma (audit 0.007 of 0.137 s), the full audit kept for the converged iteration and certificate; Hartree kernel cache and loop over the density sphere, real fields r2c (S-8): small today, needed at QE speed (QE non-eigensolver phases 0.21 s of 0.68 s).
- Validate: bit-identical for blocked traversal and caches; energy and certificate unchanged for the audit change.

### B02 Real-valued Gamma path (G-H1..H6)
- Design note and ticket breakdown (R0 to R8, orthonormal real coordinates layout, 1e-10 Ha equivalence gate): [real-gamma-design.md](real-gamma-design.md) (2026-11-06).
- Gain: G 7.4 estimate 2.1x (1.7-2.4) on the Gamma SCF: FFT 2.0x, dense products 3-4x, nonlocal 4-8x; re-profile first (the note's fractions are assumed). QE's best configuration gets the same factor, so this prevents the gap from widening rather than closing it (G 9.5).
- Staging: G-H1 paired local application and density (FFT 2x, about 1.2-1.3x total), G-H2 real Davidson dense products (about 1.9x with H1), G-H3 real block nonlocal, G-H4 pair work unit, G-H5 dispatcher with logged fallback, G-H6 half-sphere density and potential.
- Risk: new wavefunction representation and a second Davidson family (permanent maintenance); G=0 bookkeeping; needs an ADR; bit-identity with the complex path not a contract. Dependency: build it on the rank 2-3 solver structure, not before, to avoid rewriting the algebra twice.
- Validate: G section 8 (energy and eigenvalues to 1e-12 Ha against the complex path on Si2, si8, skewed, elongated, multi-projector, odd band count, gapped and buffered); QE `K_POINTS gamma` reference bundle (M0).
- 2026-11-13: the pure-build dense products now have a packed kernel behind `BAANDI_PURE_GEMM=packed` (P48B, ledger P48B-PACKED-GEMM; Mac 41-49 vs 21-30 GFLOP/s nano at the si16-si64 shapes, O). Open follow-ups: dev08 verdict on the x86 path, packed-B reuse across the column groups of one rotation (B columns are read in place today), P48D blocked orthonormalization on this kernel.

### B03 Band-level decomposition (ledger B03, G-H4)
- Gain: utilisation 23.0% (1.84 of 8 workers busy) vs 33% needed for 2x and 55% for 3x; scheduling-only attempts gave about 5% each (A40, A45, A48). Structurally different decomposition (per-band kernels for residual, preconditioning, projection; fewer larger regions) not tried; estimate unknown.
- Risk: rounding and trajectory change (A49), T=1 efficiency (A35). Avoids A33, A34, A46, A49, A50.
- First experiment: serial chain length per Davidson iteration from the `worker-analysis` build on dev08 (Mac 42.2% outside scheduled regions).

### B15 k-point symmetry reduction (S-9; user decision, ADR)
- QE reduces k points by crystal symmetry by default; repository QE references use `nosym`. For k444s the default is about 10 irreducible points against 32 (3.2x, estimate). Without it BAANDI is behind QE's best by that factor on every k-mesh case whatever else is done. Needs density symmetrisation in G space, symmetry detection and an ADR; physics-scope decision (see plan section d).

### B08 MKL pruned stick pipeline and matrix-product nonlocal (H-3, H-2; old B08, B07)
- Pruned stick FFT with plan-owned descriptors, sequential layer: per vector FFT about 0.3-0.4 -> 0.2-0.25 ms (estimate; QE with MKL pays the full 7.2 Mflop per band at 136 us); brings H per vector to or below QE's 0.38 ms with H-1. Nonlocal as one ZGEMM per block (QE 42 us per band vs ours about 0.2 ms) after H-1 only if a gap remains; changes rounding (A55 rule).
- Avoids A29 (keep the 60% rule), A34 (threads in FFT), A11/A28.
- 2026-11-13: the k-path part (nonlocal as one GEMM pair per band block) is implemented as P48C behind `BAANDI_NONLOCAL_BLOCK` (ledger P48C-NONLOCAL-BLOCK); open follow-ups: reuse Pk in `expectation_values`, avoid the packed copy by storing the columns contiguously (touches every `columns()` user), dev08 verdict.

### B10 Starting density and wavefunctions (S-3, S-4, V-8)
- Atomic charge superposition start (needs UPF atomic-density section; BAANDI starts from a uniform density), free-electron-subspace or atomic-plus-random Rayleigh-Ritz start: first iteration 9/10 steps vs QE 3/4; 3-5% (S-9.4), more on many-k (first iteration 12.6% of k444s H vectors). QE also skips the first small solve when the start is Ritz (V 2.4).
- Risk low; synthetic fixtures need a synthetic atomic density.

### B05 Thick restart with previous directions (V-4)
- Keep previous Ritz differences of unconverged roots across the 14 si8 restarts; unknown gain; do not combine with smaller capacity before measuring (H3b premise, B17).
- 2026-11-17 (BREAKTHROUGH-1-RESTART): implemented as `BAANDI_RESTART_KEEP` on the lean paths (serial and distributed lean algebra), REJECT at the current capacities (ledger RESTART-KEEP-R1B: H -0.6 to +0.6 % at SCF level, rotation up; the threads path has no contraction at capacity 6). The only open part is the capacity-2 pair (ledger RESTART-KEEP-CAP2-LEAD: -19 % H against plain capacity 2, still +17 % against capacity 3); dev08 job in requests/breakthrough-restart.md. Same-H continuation of the redo (R1-a) closed: no redo in the default runs.

### B16 Hartree-metric Pulay least squares and mixer reset survey (S-1, S-5)
- Same update, inner product in the Hartree kernel over the density sphere (G=0 excluded); 0-3 SCF iterations and, more importantly, noise tolerance that P14A showed missing; QE history 8, beta 0.7 untested here. Survey `mixing_resets` before any change. Mixing itself is 0.1-1% of the run.

### B06, B09 (unchanged from the previous backlog)
- B06 serial remainder: potential 0.019, mixing 0.014, audit 0.059 s constant over T on the Mac (about 5% of si8 T=8), EVD 0.139 s (10%); fixed-chunk reductions only; EVD via library only (V 6.2: range-selected path worth 1-4% of wall, hand-written LAPACK-style path rejected by ADR 0003).
- B09 k-point efficiency: k444s 7.1x of 16 (pure), 6.4x (MKL); largest-first ordering by per-k cost; first experiment the per-k time distribution on dev08.

### B17 Capacity and restart retune (V-6)
- Only after B01/B13: sweep multipliers 2, 3, 4, 6 on the whole set with capacity defined on active vectors; also on 20+ band cases (capacity 120 today). No default change from silicon alone.

### B11 RMM (V-7; blocked)
- QE's `rmm-davidson` is one non-converged sweep per SCF iteration with Davidson at iterations 1, 2 and every fifth; BAANDI's converged variant (A22) lost on every complete SCF; the single-sweep regime is untested but cannot certify intermediate maps; default change unlikely under the generality rule. Converged RMM at 128 or more bands remains blocked on a real case.

### B18 Exploratory or low priority
- Lane staging re-test on dev08 only (H-5; reopens A11 on its own clause); radial projector table (H-6, interpolation error 1e-8..1e-9 exceeds the 1e-12 Ha gate unless finer) [taken up by GK14-RADIAL-TABLES 2026-11-26 for the k start and k projectors: 8-point, 0.005 Bohr^-1 table, error at quadrature roundoff, see ledger] or per-shell caching (exact, only symmetric k); truncated-potential wavefunction box for early SCF iterations (H-7, same family as A21, high risk); stick-major coefficient order (H-8, microseconds); scaled G=0 (G-H9, below 1%); real-structured states at other time-reversal-invariant k points (G-H8, speculative, QE does not do it); ParO, band-by-band CG, PPCG are not recommended (V section 3; no reduction of H applications).

### B20 Setup and audit redundancy (profile-real-standard-si8.md, 2026-11-07; label O = Mac observation, E = estimate)
- C1 `ionic_local`: 86% of engine setup (setup is 21% of si8 real/standard wall, 46% of si2 Gamma wall): G-independent `erf(r)` evaluated per (G, mesh point) and the radial integral repeated for every G of the density sphere although only about 525 of 23,871 |G|^2 are distinct (cubic si8). Precompute the per-species integrand part and memoize by exact bit pattern of |G|^2. Bit-identical, E: 0.10 s of 0.615 s. Extends the exact per-shell caching idea of B18 (H-6) with measured evidence; check the NLCC core-density transform for the same pattern.
- C2 audit k-independent terms (Hartree, XC, local energy, charge, double counting, remainders) are recomputed for every k: k444s symmetry-on, 10 k, energy phase is 18% of wall, about 86% of it shared work (O: 0.076 s redundant of 0.539 s). Compute once, weight as now. Bit-identical.
- C3 energy-only Hartree in the audit and reuse of its forward FFT for `local_energy` (3 transforms to 1): E 1.3% of si8. Bit-identical.
- C4 nonlocal energy from projections only (about 3% of si8, rounding change 1e-16). C5 exact |q| memo in `NonlocalOperator::new` (2% at Gamma). C6 pair density for gamma_real (2-3%, rounding change, directional). C7 skip audit in non-converged standard iterations (8.7% si8, 18% k444s; changes history output, needs a decision). Details and ranking: [profile-real-standard-si8.md](profile-real-standard-si8.md).

### B12 MPI and GPU (blocked)
- Out of scope (AGENTS.md); GPU design note exists (A32).

## Not on the backlog (closed, see ledger)
Capacity changes at 8 bands (A08, A13; premise re-opened only for 20+ bands, B17), fuller preconditioner (A14), Ritz-change test with strict refinement (A15), multiplied tolerance factors (A21), fewer or more extra bands as default (A19), FFT lane batching/threads (A11 reopened only as B18, A34), scoped-thread schemes (A33, A35), CholeskyQR2 (A52), finer GEMM tiles (A49), blocked row-wise Nonlocal traversal (A43), band-grouped FFT (A28), hybrid RMM as default (A22). Re-open only with new evidence stated as in ledger section A.

## B17 follow-up — P31A (2026-11-07)

Capacity 2/3/4 versus default 6 has now been measured at 20, 36 and 68 bands,
with Gamma complex/real, strict/standard, five Si8/Si16 density perturbations,
a primitive cell, a k mesh, and PseudoDojo Si/MgO. Smaller capacity is
**case-dependent, not a robust default**: 2 loses strict convergence on 15
perturbations, 2/3 breach standard accuracy-ratio conditions, and none reduces
the equally weighted median wall ratio. The 4 result is indistinguishable from
noise, with only single-perturbation coverage. Keep default 6; the additive
research option remains available for controlled studies.

See [P31A study](davidson-capacity-study.md) and its ledger entry. Do not repeat
the same capacity-only sweep as an optimization proposal without new evidence.
Dense-algebra improvements that preserve the useful larger subspace remain a
separate premise; any acceptance decision still requires dev08 evidence.

## B21 to B25 thread-width follow-ups (thread-width-analysis.md, 2026-11-07; label O = Mac observation, E = estimate)
Ranking by expected gain / risk; acceptance only on dev08. Premises and ledger checks are in the analysis (section 5 and 7).
- **B21 persistent real Hamiltonian and solver workspaces (rank 1, ready).** `RealHamiltonian::from_complex` (0.036 to 0.037 s per run of 16 solves, 4.9 % of si8 and 4.4 % of si16 real T=8 wall, serial) and per-solve workspaces (1.4 % / 1.8 %) are rebuilt at every SCF iteration, including ceil(capacity/2) = 60 FFT pair plans. Keep them in the k-point state, refresh the potential only. E: 5 to 6 % at T=8, 3 to 4 % at T=1, bit-identical; may be larger on dev08 (memset 7 % of samples). Gate: bit-identical density across T; dev08 three-point alternation.
- **B22 overlap of independent phases (rank 2, ready).** (a) nonlocal image concurrent with local FFT pairs (separate buffer, same sum order; not A43/A58b, traversal unchanged), (b) energy audit of SCF iteration i concurrent with iteration i+1 (private FFT plan, orbital snapshot, final audit and convergence decision synchronous), (c) V and HV Ritz rotation in one region (A45 for the real path). E: 10 to 20 % at T=8, 0 at T=1; bit-identical by construction. First experiment: (b) alone, since the audit is 6.4 % (si8) and 9.0 % (si16) of the T=8 wall.
- **B23 k-point hybrid budget (rank 3, ready).** Inner budget = pool size with per-worker scratch, one `symmetry_map` per run. O: k444s sym -4.9 / -10.9 / -7.6 / -6.0 % at T=4/8/12/16 with identical density (one case). Remaining floor is the serial setup, potential, energy and CLI remainder (38 to 41 % at T>=10). Revives A50 with new evidence.
- **B24 narrow-block FFT split (rank 4, needs-measurement).** 54 % of H time at T=8 is in calls with at most 8 vectors (1 to 4 FFT pairs; 60 % of all H calls carry at most 4 vectors). Plane-group pieces inside one pair (revives A34 on a persistent pool); needs the dev08 wake microbenchmark first (Mac: 7 to 87 µs per region). E: 4 to 7 %.
- **B25 real-path in-block MGS (rank 5, needs-measurement).** In-block MGS is serial: 6.8 % (si8) and 12.2 % (si16) of the T=8 wall; the complex path got a right-looking first pass in A48, the real path did not. Panel-split dots change rounding; right-looking keeps it. E: 2 to 4 %. **2026-11-13:** the left-looking panel form with an L2 shape-rule width is implemented as P48D (`BAANDI_ORTHO_BLOCK`, opt-in, ledger P48D-ORTHO-BLOCK); the right-looking trailing-update form is still open if P48D's dev08 A/B shows the in-panel part dominating.
- Rejected here: longer Rayon spin before sleep (THR-3); wake latency is not the main cause on the Mac. Not recommended: spin or yield policies, scoped threads, finer GEMM tiles, band-grouped FFT.

## B26 to B30 from the QE stage-by-stage study (2026-11-07; plan [qe-parity-plan.md](qe-parity-plan.md); label E = estimate)
- **B26 QE-structured real Davidson (reopens B13 with new evidence; plan Q2, needs-measurement then isolated module).** Generalized small problem, unit-norm corrections, packed unconverged-only columns, one rectangle per matrix, capacity measured at 2, 3, 6 after steps are cheap, Ritz-change stop (standard) and contraction-based all-root residual plus the existing certificate (strict). Evidence: D0 algebra 48.3 % / 57.1 % (T=1 / T=8), QE 11 to 18 %. E wall factor 0.6 to 0.8 at T=1. Risk high on skew; stage 0 mandatory (rank and smallest Gram eigenvalue on skew, elongated, slab, Si8 replays). Supersedes the "park" decision of B13 for the real path only.
- **B27 Start and mixer package (plan Q3, extends B10 and B16).** Atomic starting density (needs UPF atomic charge parser; synthetic analytic charge for fixtures), starting rotation, beta 0.7 and history 8 survey, Hartree metric, all judged on the P29A perturbation distribution (not single trajectories). E 1 to 3 SCF iterations of 15 strict; standard median 24.5 against QE 12 is the target.
- **B28 Pair FFT parity (plan Q5, revives B08 on the pair path).** Derived 0.54 ms per H vector against QE 0.22 ms and a real-path FFT gain of only about 1.2x per band; first the dev08 microbenchmark, then pure-Rust tuning or the MKL pruned stick pipeline (only route below QE-MKL per band). E 10 to 25 % of strict wall.
- **B29 QE reference completion (plan Q1).** QE `nbnd` occupied-only, `diago_david_ndim`, mixing beta, `-nk` pools for k meshes, QE force/stress clocks, h_psi vector counts from verbose output, FFT backend of the QE build; BAANDI real+standard D0 and wake-latency microbenchmark.
- **B19/B30 update 2026-11-27 (GD1-DERIV-SPEED):** nonlocal force/stress kernel, radial tables, shared shells and the dropped duplicate force pass are done (ledger GD1). Next candidates if dev08 still shows derivatives above QE: share one XC/GGA evaluation between the force core term and the stress XC term when both are requested; reuse the SCF k bases instead of rebuilding them per k; compute forces and stress in one k pass when both are requested (one O/W/Y per atom serves both); distribute the replicated PBE grid work over ranks.
- **B19/B30 update 2026-11-27 (GD2-DERIV-XC):** shared XC between the force core term and the stress XC term, one k pass for forces and stress (joint pass, `DERIV_JOINT`), split pointwise PBE (`DERIV_XC_SPLIT`) and tabulated local/core shells (`DERIV_SHELL_TABLES`) are done (ledger GD2). Remaining derivative candidates: distribute the 10 full-grid FFTs of the derivative PBE and the local/core transforms (rounding class: plane-distributed FFT as in `gga::distributed`); split the spin PBE (`SpinXc`) the same way; joint pass for the relax adapter (`evaluate_forces_with_operations`) and the oracle/CLI legacy callers; reuse SCF k bases. Do not reuse the SCF's last XC potential for derivatives (it belongs to the mixed input density, not the delivered one). Small-wall metals: read the new spin/FD launch record (PHASE_PROFILE) on dev08 before changing fixed costs; Mac shows input parse 51 ms and binary hash 53 ms per run as the only non-SCF fixed items above 10 ms.
- **B30 Relaxation and stress workload (plan Q8, extends B19).** Make the cutoff-derivative estimate opt-in, carry orbitals and an atomic-density delta across ionic steps, tight first threshold; compare with QE clocks.

## B31 to B33 from the pair FFT efficiency study (2026-11-07; [pair-fft-efficiency.md](pair-fft-efficiency.md); label O = Mac observation, E = estimate)
- **B31 `PairFusedPlan` mirror x-slab pair plan, pure Rust, real path only (rank 1 of the FFT items, ready; ticket T1 in the document).** New plan type, no edit of `fused.rs`; unnormalised transforms, per-slab fused inverse/V/forward in a 31 KB buffer, batched rustfft calls on sticks, pair density from the same inverse half. O: 0.867 (36^3) and 0.911 (72x36x36) of the library pair apply, nine shapes 0.66 to 0.90, difference to the library 5e-16 relative; E: SCF 3 to 5 percent. Gate: energy 1e-12 Ha, complex path bit-identical, determinism across T, generality set. Refines B28.
- **B32 Lane-vectorised FFT kernels for the pair path (needs-measurement on dev08).** SIMD across slabs, four to eight lanes, radix 2/3/4/5; the Mac prototype only reaches rustfft parity at length 36 (58.3 vs 62.0 ns) and is slower at 72. Run first on dev08 (`target-cpu=native`); adopt only if at least 1.3x faster than rustfft at lengths 25 to 72. Pure-Rust kernel-only floor is 156 / 373 us per transform equivalent against QE 135 / 341 (E), so this is the only pure-Rust route to exceed QE.
- **B33 MKL raw in-place pair path, then 1D batched stick pipeline (revives B08 with new evidence; needs a dev08 microbenchmark).** The relayed "MKL full 195 us" includes two adapter copies; QE's raw call is 135 us. Step 1: raw in-place DFTI pair path without copies or full-grid zero fill (E: near QE per call); step 2: 1D batched strided stick pipeline with plan-owned descriptors (E: 97 / 235 us per transform equivalent, 0.72 / 0.69 of QE). The `mkl` feature currently routes the real path through `GridFusedPlan`. **2026-11-11: step 2 implemented as opt-in `BAANDI_PAIR_FFT=pruned` (contiguous batched lanes of the mirror-slab plan, ledger HAPPLY-PAIR-MKL; dev08 request requests/h-apply.md).** Remaining pure-build item from the same study: a pruned folded single-column plan (HAPPLY-FOLD helps only the full-grid MKL path).

## P47A preconditioner follow-up (2026-11-08, amended rule R1)

TPA mean-kinetic damping and positive-clipped kinetic Ritz denominators have
switchable implementations. The architect explicitly superseded the original
zero-baseline strict-ratio stop with a 1e-12 Ha strict accuracy floor. Both
variants pass that amended accuracy condition against a separate tighter run.

Neither is recommended: resumed Si8 standard screening rejects TPA at seed 402
(SCF 21 -> 28; H vectors 1011 -> 1276) and ritz_shift at seed 401 (SCF 28 -> 29;
H vectors 1253 -> 1283). All starts converge and satisfy unchanged standard
accuracy bounds; increased SCF count is the decisive failure. Stop before later
seeds/systems as required. No floor, tolerance, capacity or default was retuned.
Do not repeat this mechanism without a new premise addressing its increased
SCF work. See [derivation, amendment and evidence](preconditioner-variants.md).
## B34 — early-SCF fp32 H (P47B, 2026-11-08; measurement completed)

**No implementation ticket recommended by the frozen 10% median wall-saving
rule.** [The feasibility study](fp32-feasibility.md) measures ordinary/native
single-thread rustfft and faer, physical-state H errors and standard-history
eligibility. Even an intentionally generous ceiling that credits all eligible
Davidson time to H remains below 10% on the equally weighted Si8/Si16/Si32
median. The f32 lane fallback itself is slower than production fused f64 H.
The standard Ritz screen is 1e-5 Ha; moving it to 1e-7 Ha conflicts with measured
projected-root errors, and strict certificates must remain entirely fp64.

Reopen only with new evidence from a generic fused-f32 precision experiment or
a different measured hardware/workload premise. Include conversion and cached-H
rebuild costs, physical-state/correction-vector accuracy, complete fp64
certificate isolation and trajectory changes. No default, production option or
optimization was adopted. Standalone Davidson GEMM speedups are not H savings.

## Free second vector for single-column applies (idea, 2026-11-10, not tried)

Evidence: on si32 strict (perturbed seed 601) 66.4 percent of Davidson apply calls carry exactly one correction column (median 1, mean 11.4; ledger MPI-M1-CEILING). The paired real FFT transforms two columns for the price of one, so a single-column call pays for a second column and discards it. Idea: when the correction block has an odd number of columns, fill the free slot with a useful second direction (for example the previous step's correction for the same band, a conjugate-direction term as in LOBPCG, or the next unconverged-but-buffer band) so the subspace grows faster per H application. Cost per call is unchanged; benefit is fewer Davidson iterations in the long single-band tail. Needs: a switchable module, the perturbed-start distribution protocol on several systems (generality rule), the certificate and accuracy classes unchanged, and a check that it does not reduce the orthogonality margin. Not an MPI feature; applies to the threads-only path.

## Single-node sub-session findings (2026-11-10, ledger CERT-SCREEN, NL-ENERGY, DENSITY-PAIR, POTENTIAL-REUSE)
- Implemented, opt-in, awaiting the pre-registered dev08 A/B (docs/performance/requests/single-node-speed.md): `BAANDI_CERT_SCREEN=on` (skip certificate re-solves that cannot pass; Mac -17 to -21 percent wall at si32, bit-identical density) and `BAANDI_NL_ENERGY=fast` (nonlocal energy of the audit from projections only, energy phase -55 to -60 percent).
- Next ideas this opens, not tried: (a) make the SCF itself aware of the eigenvalue criterion: on si32 every seed reaches the density criterion two iterations before the eigenvalue criterion (changes 4e-7 to 1e-6 against 1e-7), so the iterations 10 and 11 are spent for a certificate that cannot pass; the estimate could steer the last iterations (tighter eigen tolerance only where it helps) instead of just skipping the re-solve; any such change alters the trajectory and needs the generality set. (b) The same first-order estimate extended to the density criterion of the certificate (needs first-order wavefunction response, which Davidson would compute anyway). (c) Warm start of the certificate solve from first-order corrected vectors (the passing re-solve still costs 14 to 20 iterations on si32 and 49 to 69 on si8). (d) The screen for the Fermi-Dirac and spin drivers (they carry their own certificate code). (e) DENSITY-PAIR (below 2 percent, rounding change).

## P48F follow-ups (2026-11-13; labels O = Mac observation, E = estimate)

- Atomic-start GEMMs Q^H HQ and Q rot over fixed column groups on the SCF pool (rounding class against the off path, thread-count independent); about 0.16 s at si64 on the Mac (O).
- Coarse-level width-1 Davidson steps fall below `PARALLEL_WORK` in real.rs (P48B owns it): a bytes-and-MAC rule per worker would let the coarse stage thread them.
- Serial reduced EVD: a fixed-decomposition parallel tridiagonalization (reduction order independent of T) or keeping B5 bordered updates; its share is larger at the coarse level because the order is the same as at the fine level.
- setup_s (0.8 s at si64 on dev08, M): ask dev08 for `engine_setup_s` versus `setup_s` first (CLI preflight versus `ScfProblem::new`); on the Mac the engine part is dominated by Ewald 0.07 s and the ionic local potential 0.05 s (O).
- MPI status 2026-11-16 (profile-2026-11-16-mpi.md): (a) is implemented as `BAANDI_MPI_TREE=rd` (MPI-FS-TREE-RD) and a stronger form of (b) as `BAANDI_M2_PANEL_GS=fused` (MPI-FS-PANEL-GS), both waiting for dev08 stage 2. New: kept start on the stick rows (MPI-FS-KEPT-START: X and H X kept in local rows across SCF steps, H_new X = H_old X + dV X with psi(r) from the density build; MPI path applies H to 29 % more columns than the serial path at equal Davidson counts).
- MPI (docs/performance/mpi-reduced-eigen-analysis.md, ledger MPI-REDUCED-EIGEN): (a) recursive-doubling `tree_sum_f64` in the FFI adapter, bit-identical to the binomial tree for power-of-two P, about half the latency per call (E si64 np=8 -0.3 to -0.5 s); (b) Pythagorean second-pass norm in the distributed panel Gram-Schmidt, 3 -> 2 reductions per column with an explicit-norm fallback (rounding class); (c) skip re-orthonormalizing a warm start whose Gram defect is at most 1e-12 (rounding class, -14 percent of the reductions); (d) dev08 microbenchmark of MKL ScaLAPACK `pdsyevd` against `dsyevd` at orders 400/600/792, P = 2/4/8, before any library adoption (the only change that shortens the reduced-eigen critical path).

- C3a (2026-11-16, from C2-REUSE-PSI / C3-NL-ENERGY-FROM-PROJECTIONS): lean solver keeps P^T V beside H V (copy after each apply via `copy_last_projections`, rotate at contraction and return) so the audit's nonlocal energy reads (P^T V) C and the C2 kept start carries P^T X; removes one n x n_proj x bands GEMM per map (si64 E 1.5 to 2.7 %). Own switch, counts first. IMPLEMENTED 2026-11-16 as `BAANDI_LEAN_PROJ_CACHE` (default off; ledger C3a-LEAN-PROJ-CACHE); dev08 A/B pending.

## MPI scaling gap follow-ups (2026-11-17, MPI-SCALING-GAP; docs/performance/mpi-scaling-gap.md)
- Distributed full-grid FFT for the replicated Hartree potential, audit Hartree energy and Coulomb metric of the M2 path (`DistributedPairFft` over the full box: fill planes from the replicated density, `to_reciprocal_plane`, `finish_to_reciprocal`, local 4 pi / G^2, one tree sum for the energy, inverse and `gather_planes`); replicated potential/audit/mixing are about 10 percent of the Mac np=4 wall and grow with np on dev08 (MPI-SCALING-BEYOND-8). Rounding class.
- Lazy orbitals on the M2 path (no N x B allgatherv and complex unpack per solve; audit kinetic energy from local rows): IMPLEMENTED 2026-11-17 as `BAANDI_M2_LAZY` (default off; ledger MPI-LAZY-ORBITALS-*, design docs/design/mpi-lazy-orbitals.md, dev08 request docs/performance/requests/mpi-lazy-orbitals.md). Remaining allgathervs per map: the solver's per-solve row-index exchange (N doubles, 1/B of the removed block; cacheable per operator), the M2 density plane gather and the XC potential gather (full grid each; they go away only with the distributed full-grid stages above).
- Distributed back-transform of the reduced eigenproblem: tridiagonalize replicated (about 4/3 m^3), each rank computes bands/P eigenvectors of the tridiagonal problem and back-transforms them, one allgather of the m x bands coefficients; needs dsytrd/dstemr/dormtr on MKL and a pure-build route.
- Stick balance rule on stick count when `M2_STICKS=mirror` (si32 np=8 imbalance 1.058) if dev08 shows it.

## Fast-loop follow-ups (2026-11-16, threads/one-core owner; profile-2026-11-16-threads.md)
- VERIFY-REUSE: feed regressed si16 on dev08 stage 1; variants screened (feed2 and replace rejected); `half` is the stage-2 candidate (requests/verify-reuse.md).
- Passed-verification cost (8 to 16 %): tau_v is fixed by ADR 0043; a user decision, not a switch.
- Stop-estimator underestimation (10 to 20x below the verification residual at failing stops, si32 seeds 603/605/610): a better predictor of F(rho_out) - rho_out (for example from the Pulay subspace) would avoid futile verifications; analysis only so far.

## MPI remaining terms (2026-11-17, MPI-REMAINING; docs/performance/mpi-remaining.md)
- Overlap-matrix (QE-style generalized reduced problem) Davidson on the M2 path: IMPLEMENTED 2026-11-17 as `BAANDI_M2_LEAN_ALGEBRA=on` (default off; ticket MPI-OVERLAP-DAVIDSON, docs/design/mpi-overlap-davidson.md, ledger MPI-OVERLAP-DAVIDSON-*); two tree sums per iteration (S border before H, H border after; the one-sum variant was rejected because it applies H to rejected blocks); dev08 request docs/performance/requests/mpi-overlap-davidson.md.

## Overlap-matrix Davidson follow-ups (2026-11-17, MPI-OVERLAP-DAVIDSON)
- Cheaper stage 1 at capacity 3: on the late warm solves at np >= 4 (capacity 3 per band) 14 to 15 of 22 to 23 appends fail the pivot test (si32, si64, O), and stage 1 runs a norm reduction, the in-block Gram-Schmidt chain of the block and a second S border. Check first whether the serial lean solver at capacity 3 shows the same rate (threads path, same density); then candidates that keep the rank-decision meaning: in-block dependence from the k x k Gram of E' (one reduction, Cholesky with the same per-column floor, MGS only for the columns it rejects), or a second projection pass of E' through the factor before the in-block step. Counters: `dist.lean_reprojections`, `reductions.ortho_*`, `panel_*`.
- Start with the factor test by default on the lean path (`M2_START=gram` semantics without the switch): the start orthonormalization is most of the remaining tree sums with the switch on (si32 np=8 production run: 2033 tree sums per rank, of which panel_dot 1426; with the Gram start 511, O); with the pivot rule any start block with a factorable Gram matrix can be the basis.
- Replicated reduced work after the switch: the factor extension and the border of A add O(m^2 k) per iteration on every rank (reduced_eigen stage +6 ms of 31 at the warm si32 np=8 solve, O); a root-only factor plus a broadcast of the coefficients would need all m eigenvectors for stage 2 (or a broadcast of T on fallback only).
- Ladder coarse stage with the lean algebra: the coarse M2 solves are standard solves, so the switch acts there too; measure the coarse-session share at stage 1 (si32 np=8, 25 percent of the wall before).
- M2 output XC energy without the potential allgatherv: the audit's output-density XC needs the energy only (no core density on the M2 XC path), but `m2::xc` gathers the full potential (N doubles per map); bit-identical, small (one allgatherv of 1.5 MB at si32 per map).
- Distributed Pulay mixer (history on grid blocks, one tree sum of the m x m Gram and m right-hand sides, one allgatherv of the mixed density): after MIX_INCREMENTAL the remaining replicated mixer work is O(m N); checkpoints need the full history (gather on capture).
- Ladder coarse stage at np >= 4 (25 percent of the si32 np=8 wall, flat from np 4 to 8; serial coarse stage 1.19 s is slower than distributed 0.72 s on the Mac): per-iteration costs of a 4189-PW problem dominate; the overlap-matrix Davidson acts here too.

## Per-core kernel follow-ups (2026-11-17, PERCORE-KERNELS; docs/performance/percore-kernels-2026-11-17.md)
- `BAANDI_LEAN_BORDER=fused` and `BAANDI_PACK_FIXED=on` are stage-0 candidates on the Mac (bit-identical there); dev08 stage 0/1 (MKL and x86 pure) in requests/percore-kernels.md. The MKL per-kernel table of the lean default (Part A0) decides the next kernel item.
- Fresh-basis borders (start, contraction): V^T V and V^T (H V) could also be one product V^T [V | H V], and their diagonal blocks are symmetric (half the flops of the k x k block, dsyrk on MKL); about 2 % of a warm solve (factor_fresh), more in solves with contractions. Not implemented.
- Narrow NN rotations stay far below the GEMM rate even after PACK_FIXED (k <= 4: GEMV shaped); a streaming kernel for X_s = V C_s with k <= 4 (one pass over V, like NARROW_PROJ for the TN side) is the remaining pure-build lever; on MKL check the per-call rate first (Part A0).
- Algebra work per H column is about twice QE's (22 against about 10 to 12 MFLOP at si32, basis size): any capacity change must pass the SCF-level screens (A13 / SERIAL-CAPACITY-SCREEN premises), not stage 0 alone.

## Bounded subspace SCF follow-ups (2026-11-17, BREAKTHROUGH-3; docs/design/bounded-subspace-scf.md, ledger BT3-*)
- The `Advance` contract exists (`RealDavidson::advance_reusing`, `BAANDI_SCF_ADVANCE`, default off) and lost at stage 1 on si16 (three seeds) and si-skew: the standard ratchet already makes exact warm maps cheap (3 to 7 Davidson iterations), so bounding the block count saves 20 to 30 percent of a map at best and the completions, extra SCF steps and extra verification attempts pay it back. Do not re-run the same screen; the premise that changes it is a cheaper FIRST block (review R2: residual-rank selection of the first correction, so a 1-block advance is not an all-band block) or a solver whose per-block cost is far below the exact solve's (R4 periodic RR). With such a first block the Advance contract can be re-screened in one run (`BAANDI_SCF_ADVANCE=on`, `SCF_ADVANCE_GROWTH=off`).
- MPI path: not implemented (declines with `distributed solver`); only worth adding after a serial win.
- Mixer re-weighting of approximate pairs (design section 6) was not needed: no accepted advanced map showed the mixer learning noise (E tracked the exact path). Keep as a note, not a task.
- Ladder (review R5): judge the coarse stage by the fine hand-off quality (fine H columns after hand-off), not by coarse cost; BT3-COARSE shows a cheaper coarse stage paying back in the fine stage.

## Other and init follow-ups (2026-11-17, BREAKTHROUGH-4-OTHER-AND-INIT; docs/performance/profile-2026-11-17-other-init.md, ledger BT4-*)
- dev08 decides the three bit-identical removal switches (`M2_LIGHT_STATE`, `EWALD_CACHE`, `M2_LAYOUT_CACHE`; request breakthrough-other-init.md). If light state wins, the next step is to stop building the serial real Hamiltonian at all on the M2 path (a half-sphere-only state; the light build still converts every projector to real, `real_state` 0.011 s fine O).
- Keep the slab of the partly occupied pair in the M2 density build so that only the empty pair is transformed both ways at `kept_update` (fine 0.072 s O); the mirror exchange stays while any pair is not kept.
- Replicated construction work per rank (O, two constructions at si32 np=8): `ionic_local` 0.044 s (per-G independent: G-range split plus one allgatherv, bit-identical), the remaining Ewald evaluation (reciprocal G split: rounding class of E only, density untouched), `DistributedRealHamiltonian::new` packing every projector over all N before keeping the local rows.
- Atomic start at np=8: block orthonormalization of the 128 span columns (Gram + Cholesky with the dependence rule) instead of the 128-broadcast chain (rounding class of the trajectory: stage 1 to convergence); split the Ritz rotation Q C by band columns (one allgatherv).
- `between` (O, si32 np=8): session set-up 0.036 s, stop rule 0.019 s, ladder hand-off 0.028 s: replicated full-grid passes; distribute only if dev08 shows them above 1 %.

## General k / spin / FD path (GSPEED-PLAN, 2026-11-26; plan docs/performance/general-path-speed-plan.md)
Labels E = estimate before the dev08 kprof profiles. Bit-identical batch first: GK1 k atomic-start radial memo
(B12 E 70-250 s), GK3 complex Davidson workspace reuse and merged error agreements, GK4 per-k energy terms at the
audit only, GK5/GK6 one potential move per axis and shared full-grid/FFT buffers, GK7 diagnostic symmetry
projection opt-in, GK8 lean energy/permap/verification gate on the spin/FD path, GK2 one k start shared by both
spin channels. Layout/algorithm: GK9 k-point pools (ADR 0041 option D; E 1.5-2.5x of the k loop at MPI8), GK10
border Gram in the complex Davidson, GK11 generalized Ritz k start, GK12 time reversal and moment-preserving
symmetry for collinear spin (policy + ADR), GK13 QE-like FD empty-band threshold, GK14 fine radial tables (reopens
B18 H-6 with the new 1e-10/5e-9 classes as premise). First step GK0: read kprof (setup split, per-map phases,
Davidson iterations per k vs QE, IrI3 spin/Pulay/k count, QE layout) and request the QE -nk scan.

## After GK16 (GK17-SCF-GAP, 2026-11-27; plan section 10)
Gap is cost per unit of work, not H columns (IrI3: 0.50x QE's H columns, 0.81x its FFTs, same core-s per k-iteration).
- dev08 A/B first: `K_KEPT_SLABS=on` at `K_KEPT_SLAB_MB` 1024 and 4096 (peak RSS per rank), `FFT_ONE_RANK_DIRECT=off`
  vs on, `M2_FFT=mkl` (MKL lanes on the k operators; disables the direct sticks), `K_FIRST_EMPTY_RULE=on` (GeTe map 1).
- GK9 stage 2: sub-communicators for the remainder k (IrI3 51 k on 8 ranks: critical rank 7 k vs mean 6.375, ~9 %).
- Kept slabs on the spin/FD path (`electronic_map`, same mechanism; FD/spin cases only).
- Middle-stage FFT movement (gathers/scatters/clears 12-15 % on the Mac k profiles): fuse the y-gather with the
  z-row transform only where it keeps bits (FFT_FUSED_ROWS was rounding class and slower on Gamma, P3f).
- If dev08 shows the rustfft kernels slower than QE's per transform: calibrated per-length kernels for private
  (k) lanes (today `FFT_KERNEL_BY_LENGTH` acts only inside a plan-cache scope, which the k start does not open).

## After GK18 (GK18-SETUP-SMALL, 2026-11-27)
- Lazy delivery on the spin/FD k driver (`electronic.rs` BandBlock; IrI3: 51 k x 42 bands x ~15k rows broadcast to
  every rank in energy-only runs) with materialize hooks in forces/stress/checkpoint/NSCF for electronic results.
- Shared per-rank scratch for the k operators (FFT stick blocks, `back`, GEMM workspace) instead of one set per k
  (plan GK6): memory per rank and first-touch page faults in the start H apply on Linux.
- Per-k `FftGridSpec::new` (1.77 MB `fft_to_g` per k at 48^3) in the k start: a k-only grid spec without the inverse
  map (the k operators read `g_to_fft` only).
- If the GLIBC_TUNABLES arm of the GK18 dev08 request wins: set the allocator thresholds at process start (mallopt
  in the MPI/CLI adapter) or evaluate a maintained allocator crate (ADR 0003 checks).
- Start H apply and SCF H apply: FFT kernel speed per transform on dev08 (`M2_FFT=mkl` arm of GK17).

## After GF2 (GF2-FD-SCF-COST, 2026-11-27; ledger GF2-FD-SCF-COST)
- dev08 A/B of `FD_EMPTY_THRESHOLD` (weighted default, percent, off) and `FD_EMPTY_VERIFY=on` on the MP metals (Au,
  Ru, Sr, Hg, Ta6Zn7, Hf, TiGa), verify off and on; decide FD_EMPTY_VERIFY (verification-map eigen thresholds only,
  criterion unchanged) with the orchestrator/user.
- Ta6Zn7 k count: BAANDI 34 irreducible k against QE 26 (nsym 12); read BAANDI's `reference_scf.symmetry`
  (group_order, irreducible_k_count) and preflight warnings; 34 / 26 = 1.31 more k, 5 against 4 k on the critical rank
  at MPI8 (1.25x), the largest single factor if the group or mesh reduction differs.
- GK9 stage 2 remainder layout (few k, large npw): Ta6Zn7 34 k on 8 ranks = 4 whole k per rank plus 2 remainder k
  each split over 4 ranks: about 76k rows per rank against 87k (`auto` npool 8; E -12 % of the k loop); Hg 76 k:
  10 against 9.5 k (E -4 %).
- Reduced problem of the complex Davidson on MKL builds: Cholesky, two triangular solves and the full Hermitian EVD
  (all n eigenvectors of an order up to 2n) run on faer even with `mkl`; an MKL route (zpotrf/ztrsm + zheevd, or
  zheevr for the lowest n of up to 2n) behind a switch; Mac stage share 2-3 % of the Davidson time (TaZn13, Au),
  more on dev08 where the GEMMs are MKL. Needs an MKL machine for the FFI tests (dev08).
- A cheaper guard than the FD buffer against missed top states: the buffer costs 4 of 18 columns (Au) and 14 of 145
  (Ta6Zn7) per Davidson call, but removing it lets the top retained bands converge to wrong eigenpairs (GF2); a
  smaller buffer for large n or a buffer only in the first and verification maps are untested options.
- Forced first Davidson iteration (GF2 counts, Mac, maps >= 2, weighted arm): every call corrects every solved
  column once before any convergence test, 62-93 % of the SCF H columns (Au 0.73, Ru 0.74, Sr 0.93, Al slab 0.62,
  TaZn13 0.86); the buffer is 4 of those 18 columns on Au. A residual-screened first iteration (kept HX makes
  r = HX - eps X free before any H application: skip a column whose predicted Ritz change from its residual is far
  below its threshold, keep the Ritz-change exit test) would cut this for every path; it changes what the first
  iteration does, so it needs a safety argument (no stale buffer columns, no missed states) and the class matrix.

### GS1D-B (2026-11-27): reciprocal-space density symmetrization for grid-incompatible rotations
QE keeps rotations that do not map a non-cubic explicit FFT grid onto itself (Cl2 on [40,36,80]: QE 8 ops / 42 k, BAANDI 4 ops / 68 k). BAANDI requires grid-mapping rotations because XC is evaluated on the real-space grid. Option: a switch that keeps such rotations and symmetrizes the density in reciprocal space; adopt only if reduced-vs-full-mesh energies stay within strict 1e-10 / standard 5e-9 Ha on several cells. Deferred: affects explicit non-cubic grids only.

## After GM1 (GM1-FD-MIXING, 2026-11-27; ledger GM1-FD-MIXING)
- dev08 A/B of `SCF_ETHR_GUARD` off / forecast / redo (verify on) on Cs2TcI6, Ta6Zn7, TiGa, Au, Ru, Fe; wall and maps.
  Insulators on the fixed-occupation k driver are not affected by the switch (control rows only).
- Sensitivity on the QE side: QE Cs2TcI6 with mixing_beta 0.69 / 0.71 or diago_thr_init changed by 2x; if QE's map
  count also moves by 2-4, the 15 vs 11 gap is trajectory noise, not a BAANDI deficit.
- Cs2TcI6 it1 dr2 0.405 (BAANDI) vs 0.253 Ry (QE): check the start density (QE "starting charge" line, negative
  atomic charge, radial cutoff of PP_RHOATOM) and the band count/cutoff identity; the Mac model says map noise is
  not the reason.
- Overshoot after large drops (dr2 x3-30 on new strongly responding directions): candidates are an adaptive scalar
  initial inverse Jacobian from the newest secant (beta_k = min(beta, -<dx,dF>/<dF,dF>)) or an energy-based adaptive
  damping (Herbst and Levitt 2022); both are research items that need the full FD and insulator class matrix.
- Start-symmetry check with clipped atomic densities (side finding GM1): compare within the density sphere or after
  the same sphere projection, or clip after projection; the model Cs2TcI6 at 30 Ry is refused today.
