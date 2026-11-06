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

### B03 Band-level decomposition (ledger B03, G-H4)
- Gain: utilisation 23.0% (1.84 of 8 workers busy) vs 33% needed for 2x and 55% for 3x; scheduling-only attempts gave about 5% each (A40, A45, A48). Structurally different decomposition (per-band kernels for residual, preconditioning, projection; fewer larger regions) not tried; estimate unknown.
- Risk: rounding and trajectory change (A49), T=1 efficiency (A35). Avoids A33, A34, A46, A49, A50.
- First experiment: serial chain length per Davidson iteration from the `worker-analysis` build on dev08 (Mac 42.2% outside scheduled regions).

### B15 k-point symmetry reduction (S-9; user decision, ADR)
- QE reduces k points by crystal symmetry by default; repository QE references use `nosym`. For k444s the default is about 10 irreducible points against 32 (3.2x, estimate). Without it BAANDI is behind QE's best by that factor on every k-mesh case whatever else is done. Needs density symmetrisation in G space, symmetry detection and an ADR; physics-scope decision (see plan section d).

### B08 MKL pruned stick pipeline and matrix-product nonlocal (H-3, H-2; old B08, B07)
- Pruned stick FFT with plan-owned descriptors, sequential layer: per vector FFT about 0.3-0.4 -> 0.2-0.25 ms (estimate; QE with MKL pays the full 7.2 Mflop per band at 136 us); brings H per vector to or below QE's 0.38 ms with H-1. Nonlocal as one ZGEMM per block (QE 42 us per band vs ours about 0.2 ms) after H-1 only if a gap remains; changes rounding (A55 rule).
- Avoids A29 (keep the 60% rule), A34 (threads in FFT), A11/A28.

### B10 Starting density and wavefunctions (S-3, S-4, V-8)
- Atomic charge superposition start (needs UPF atomic-density section; BAANDI starts from a uniform density), free-electron-subspace or atomic-plus-random Rayleigh-Ritz start: first iteration 9/10 steps vs QE 3/4; 3-5% (S-9.4), more on many-k (first iteration 12.6% of k444s H vectors). QE also skips the first small solve when the start is Ritz (V 2.4).
- Risk low; synthetic fixtures need a synthetic atomic density.

### B05 Thick restart with previous directions (V-4)
- Keep previous Ritz differences of unconverged roots across the 14 si8 restarts; unknown gain; do not combine with smaller capacity before measuring (H3b premise, B17).

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
- Lane staging re-test on dev08 only (H-5; reopens A11 on its own clause); radial projector table (H-6, interpolation error 1e-8..1e-9 exceeds the 1e-12 Ha gate unless finer) or per-shell caching (exact, only symmetric k); truncated-potential wavefunction box for early SCF iterations (H-7, same family as A21, high risk); stick-major coefficient order (H-8, microseconds); scaled G=0 (G-H9, below 1%); real-structured states at other time-reversal-invariant k points (G-H8, speculative, QE does not do it); ParO, band-by-band CG, PPCG are not recommended (V section 3; no reduction of H applications).

### B12 MPI and GPU (blocked)
- Out of scope (AGENTS.md); GPU design note exists (A32).

## Not on the backlog (closed, see ledger)
Capacity changes at 8 bands (A08, A13; premise re-opened only for 20+ bands, B17), fuller preconditioner (A14), Ritz-change test with strict refinement (A15), multiplied tolerance factors (A21), fewer or more extra bands as default (A19), FFT lane batching/threads (A11 reopened only as B18, A34), scoped-thread schemes (A33, A35), CholeskyQR2 (A52), finer GEMM tiles (A49), blocked row-wise Nonlocal traversal (A43), band-grouped FFT (A28), hybrid RMM as default (A22). Re-open only with new evidence stated as in ledger section A.
