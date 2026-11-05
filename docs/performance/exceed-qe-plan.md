# Plan: match and then exceed QE 7.6 at equal cores and equal accuracy

Written 2026-11-05 by the performance owner. Target (AGENTS.md, Performance target): BAANDI must match and then exceed QE 7.6 on the benchmark set at equal core count and equal accuracy, against QE's BEST configuration per case. Sources: relayed dev08 data (`tests/reference/perf/dev08_final_relayed.md` parts 4, 5 and 6, `qe_iteration_data_relayed.md`), the P14A/P14B tables, the four QE study notes G, H, S, V of 2026-11-05 (see [backlog.md](backlog.md) for the key), the ledger. Every number marked (est.) is an estimate, (der.) is arithmetic on recorded numbers, (rel.) is relayed from dev08. Nothing here is a new measurement. Items needing the remote session are collected as M0 in the backlog.

## Update 2026-11-06: what is measured now (read first)
Labels: **M** = measured on dev08 and relayed (part 7 of dev08_final_relayed.md); **O** = Mac observation (stored in the repository, not acceptance); **E** = estimate.
1. **Native build flags (M, part 7, commit 84bfbf3).** `RUSTFLAGS="-C target-cpu=native"` on the pure build: si8 T=1 5.423 -> 4.122 s (1.316x), T=8 4.019 -> 2.756 s with LTO added (1.458x); k444s T=1 1.229x, T=8 1.121x; across the eight case/backend/thread combinations pure 1.12-1.40x, MKL 1.05-1.10x; energy, density hashes, SCF/Davidson/H counts bit-identical (24 comparisons). 46% of the saving is `subtract_projection` (0.822 -> 0.224 s, the only entry above 30 samples); `NonlocalOperator::add` (0.603 -> 0.645 s, directional) is unchanged. Pure native beats the default MKL build on si8 (4.122 s vs 4.814 s at T=1; 2.756 vs 3.370 s at T=8). LTO alone 1.01-1.06x. Against QE best (0.366 s): 11.3 (was 14.9).
2. **Standard convergence profile P17A (ADR 0035, accepted; Mac O).** si8: strict 2.512 s, 345 Davidson steps, 2265 H vectors, 15 SCF; standard (threshold 5e-9 Ha) 0.797 s, 69 steps, 742 H vectors, 14 SCF: wall factor 0.317 (3.15x), work factor 0.33 (H vectors); delivered energy error 1.0e-9 Ha against BAANDI's own tight reference (QE estimate for si8 2.7e-9 Ha, stored, not an exact QE guarantee), eigenvalue error 1.7e-6 Ha, density L2 1.2e-5. k444s: 2.456 -> 1.401 s (0.57), 1827 -> 866 steps, 11181 -> 5864 H vectors, energy error 1.3e-10 Ha. Result is uncertified (`certified: false`), Hartree-metric mixing not adopted (E4 failed a 1e-12 Ha gate on MgO), skew and 16-atom mesh fixtures exhaust the iteration limit, slab diverges, three stress fixtures have no tight reference. The dev08 comparison against QE at QE's measured delivered accuracy is requested and not yet received.
3. **P15A replay (O/work counters).** On fixed si8 histories, replacing the residual stopping predicate by the Ritz-value change at identical numeric tolerances: 283 -> 39 rows, 2063 -> 749 H vectors; QE-like empty thresholds 38/729; for the first 14 SCF iterations the relayed QE threshold ladder gives 89 rows (75 expansions) against 257 (243) for BAANDI native, QE relays 83. The stopping rule therefore explains most of the 290 vs 83 step-count gap; the per-step convergence rate is not the main difference. The replay also found that small residuals do not imply spectral completeness (dense oracle: low roots missed, si8 tight replay occupied discrepancy 0.11 Ha; skew fixture 1.51 Ha).
4. **P15C spectral guard.** Final-only option `spectral_guard = "final"` (default off): +90 H columns on si8 (3.97%) and +2618 on k444s (23.41%), SCF counts unchanged; realistic final states (six Si/MgO fixtures) are complete (max occupied discrepancy 3.6e-14 Ha); the 3% target is not met, so it stays optional; the skewed fixture still cannot converge.
5. **k-point symmetry P19A (O).** Shifted 4x4x4 Si: 64 -> 32 (time reversal) -> 10 irreducible points; Mac k444s T=1: 3.328 s (TR only) -> 1.106 s (crystal), ratio 3.01, operator applications 12973 -> 4040; reduced vs full-mesh energy difference 1.24e-14 Ha. The 10-point count is confirmed only by our own detection; QE with symmetry on is the reference for k-point cases and the dev08 request is pending (no stored QE symmetry-on k444s output).
6. **P24A locking (O, work counters; option `locking = "auto"`, default off).** Complex multiply-adds -6.0% on si8 and -12.5% on k444s against the 30-40% the solver-note model predicted; H applications rose (2265 -> 2809, 11181 -> 13500); failed the accuracy gates (syn-s2 energy changes 4.4e-3 Ha, si8 buffer roots 3.7e-9 Ha above the 1e-9 Ha dense cap, skew diverges); not an optimisation. The V cost model over-predicts removable work, and the scalar MGS it was meant to attack shrank by itself under native flags (point 1).
7. **Forces, stress, relaxation, cell relaxation exist (O).** Mac cost observations (one sample): forces are 28.3% of SCF-plus-force time on si8 (force time about 0.39 x SCF) and 10.9% on k444s; the analytic stress pass is 58.4% of SCF-plus-pass on si8 and 44.6% on k444s (Pulay cutoff scans are extra). Once SCF is 3x faster (standard profile), forces and stress become a larger share than SCF itself.

## (a) Gap model

### a.1 si8 Gamma, one core, conv_thr 1e-8 (fdce016, rel.)
BAANDI pure 5.801 s, MKL 6.016 s; QE internal 0.684 s (complex Gamma, MKL-DFTI, nbnd 20, nosym). Ratio 8.48 pure, 8.80 MKL. With P14B (rel.): pure 5.463 s (7.99), MKL 4.814 s (7.04).

**Measured QE best configuration (part 6, rel.): QE real Gamma path (`K_POINTS gamma`), si8, one core, 0.366 s at conv 1e-8 (12 SCF iterations, 81 Davidson steps, h_psi 94 calls 0.18 s, fftw 964 calls 0.13 s, real Davidson 0.23 s) and 0.500 s at 1e-12 (18 iterations); same total energy; 1.87x faster than the complex path at N=1 (1.3-1.9x across cases). The real target is therefore BAANDI/QE-best = 14.9 (pure) and 13.2 (MKL) at 1e-8 with fdce016 numbers, and 5.463/0.366 = 14.9, 4.814/0.366 = 13.2 with P14B (the relay quotes 14.9 and 13.2); at 1e-12 10.9 and 9.6. The h_psi call count is nearly equal in the two QE paths (94 vs 98): the real path halves the cost per call (two real bands per complex FFT, real symmetric Davidson), it does not reduce calls. QE build: plain -O3, MKL run-time dispatch, no -march=native; BAANDI on dev08 was always a default release build (no RUSTFLAGS, no `.cargo/config`), so both are default builds; a target-cpu/LTO experiment is being measured by the remote session.**

| Factor | QE | BAANDI | Ratio | Source / label |
|---|---|---|---|---|
| SCF iterations | 14 | 15 | 1.07 | qe_iteration_data_relayed.md |
| Davidson steps (total) | 83 | 290 (Mac local: 340 + 68 certificate) | 3.49 | same; per SCF iteration 5.93 vs 19.33 (3.26) |
| H vectors per step | about 11.4 (about 950 / 83, est.; upper bound 12.9) | 7.14 (2070 / 290) | 0.63 | H section 3 (est.) |
| H vectors total | about 950 (est.) | 2070 | about 2.2 | H section 3 (est.) |
| Cost per H vector | 0.38 ms (h_psi 0.36 s / 950, der.) | 0.5-0.65 ms (est.) | 1.3-1.7 | H section 3; FFT at parity (MKL), nonlocal about 5x QE (0.22 vs 0.042 ms) |
| H application time | 0.36 s | 1.0-1.6 s (der.) | 2.9-4.4 | = vectors 2.2 x cost 1.3-1.7 |
| Davidson algebra around H | 0.11 s (c_bands 0.47 minus h_psi 0.36) | 4.0-4.7 s (der. by subtraction) | 36-43 | the dominant factor |
| Non-eigensolver phases | 0.21 s (init 0.05, sum_band 0.07, v_of_rho 0.03, mix_rho 0.04, rest) | 0.1-0.2 s (user quote, Mac phases) | 0.5-1.0 | S section 1, 7 |
| Sum check | 0.36 + 0.11 + 0.21 = 0.68 | 1.0-1.6 + 4.0-4.7 + 0.1-0.2 = 5.8 | | by construction |

Decomposition of the algebra factor 36-43: step count 3.49 x per-step multiply-adds (V model: about 5600N for BAANDI with k about 7, q about 70, m 20; about 1100N for QE with k about 9, q about 30; ratio 5.1, est.; model total 1.1e10 vs recorded 9.49e9 2mnk, within 20%) = 17.8x; the remaining 2.0-2.4x is lower efficiency per flop (scalar sequential MGS, copies, large EVD, parallel dependency chains; est.). Capacity is `max(48, 6*bands)` = 120 for 20 bands, so the mean basis is about 70 against QE's 2x bands (about 30).

Remaining factors not visible in the time split:

| Factor | Facts | Size |
|---|---|---|
| Accuracy definition | QE: Ritz-value-change test, up to 5 unconverged roots accepted, no certificate, no gap check, last threshold 5.9e-11 Ry (residual about 1e-5, est.). BAANDI: residual in Ha, certificate solve (Mac local: 68 of 408 steps = 16.7%, 204 of 2416 H vectors = 8.4%, der.), 4 extra bands (24.9% of H vectors), tolerance 100-1000x tighter late (V 6.3, est.). QE alone moving conv_thr 1e-8 -> 1e-12 costs 0.684 -> 0.912 s (+33%) | not measured; plausibly 1.3-2x of the ratio (est.); E1/E2/replay measure it |
| QE configuration | the 8.48 compares against QE's complex Gamma, nbnd 20, MKL-DFTI (no FFT pruning, 0.29 s of 0.68 s is FFT). Measured real Gamma path: 0.366 s (part 6), FFT 0.13 s (2.22x fewer FFT calls), real Davidson 0.23 s of which h_psi 0.18 s, i.e. QE algebra only about 0.05 s. Not yet measured: occupied bands only (nbnd 16, est. 1.1-1.25x), FFTW3 build (prunes, est. a further 0.05 s of the 0.13 s FFT part), symmetry | QE best on si8 = 0.366 s measured (est. 0.28-0.36 s if nbnd 16 and FFTW3 help); ratio vs best with P14B pure 14.9, MKL 13.2 (measured/rel.) |

### a.2 k444s, one core, conv_thr 1e-8 (fdce016, rel.)
BAANDI pure 6.127 s, MKL 5.034 s; QE 1.204 s; ratios 5.09 and 4.18. At N=16: 0.834 / 0.788 s vs 0.384 s: 2.17 / 2.05.

| Factor | QE | BAANDI | Ratio | Source / label |
|---|---|---|---|---|
| SCF iterations | 5 (8 at 1e-12) | 10 | 2.0 (1.25 at 1e-12) | D08 part 4 |
| Davidson steps | about 790 (24.7 per k x 32, der. from the six relayed averages; 5 vs 6 iterations ambiguous) | 1686 (+141 certificate) | 2.1-2.3 | qe_iteration_data_relayed.md, p14a |
| Time per step | 1.204 / 790 = 1.5 ms | 6.127 / 1827 = 3.4 ms | about 2.2-2.4 (der.) | mixes H, algebra, phases |
| Non-eigensolver phases | about 12% (electrons 1.09 vs c_bands 0.95) plus init 0.10 s | density 3.7-6.2% + audit 6.6-12.9% + potential/mixing | about parity | p14a, S section 1 |
| QE configuration | default symmetry: about 10 irreducible k for the shifted 4x4x4 fcc mesh (standard result, not recorded in this repo, to be verified) vs 32 with time reversal only | 32 (nosym) | QE best about 1.204 x 10/32 = 0.38 s plus symmetrisation overhead (est.) | S section 3.1; repository QE references use nosym |

Ratio vs QE best on k444s, N=1: 6.127 / (0.38-0.45) = 14-16 (est.).

### a.3 What the model says
(1) Hamiltonian application is not the main gap (V and H agree: 3.7 ms vs 3.5-4.5 ms per step). (2) Davidson algebra is: about 70-78% of BAANDI's time vs 16% of QE's. (3) Step and vector counts are about 2.2-3.5x and partly definition. (4) QE's best configuration lowers the target by 1.3-1.9x (real Gamma, measured, 1.87x on si8 N=1) on Gamma-only and by about 3x (symmetry, est.) on k meshes; our earlier ratios understated the gap: on si8 the real target is 13-15, not 7-8. For QE's real path, algebra around h_psi is only about 0.05 s (0.23 - 0.18): the BAANDI algebra factor against QE best is about 80-90 (der.). (5) Caveat on the model: dev08 pure si8 T=1 gained only 1.06x from a 56% cut in explicit GEMM flops (A55) while MKL gained 1.26x; the flop model over-predicts time gains for pure Rust until the post-P14B profile is known.

## (b) What the target means operationally

### b.1 Benchmark set
Protocol cases gamma, k222s, k222, k444s, si8 (QE references exist, complex Gamma, nosym); the 19-case suite of ticket P10G for generality (BAANDI only today); additional cases still missing and required for any "exceed" claim: a 16-atom Gamma case and the 32/64-band supercells (`super4`, `super8`), a 128-or-more-band case, a slab and an elongated cell with QE references, and Gamma-only cases above 0.2 s QE time (timer resolution 0.01 s; entries below about 0.3 s are excluded). Additional workloads recommended for the set (second phase, after SCF parity): single-point forces, the stress tensor and a short BFGS relaxation of a displaced cell, each with QE's own `forces`/`tstress`/`relax` runs as reference (cost shares measured on the Mac: forces 28.3% of SCF+force on si8, 10.9% on k444s; stress 58.4% and 44.6% of SCF+pass). Layouts: 1, 4, 8, 16 cores; QE with the same number of MPI processes (or the best of pools, plane-wave and OpenMP at that core count).

### b.2 Equal accuracy: options
1. **Status quo (equal conv_thr label).** QE conv_thr 1e-8 vs BAANDI strict defaults with certificate. Simple, but compares different guarantees; BAANDI is penalised by its certificate and buffer bands. Already approved; keep as the "strict" column.
2. **Delivered-accuracy comparison (recommended).** For each case run each code to its own tight reference (QE 1e-12, BAANDI strict) and require, for the timed run, the same bounds on: total-energy error against that code's own reference (QE at 1e-8 is already below the printed 8 digits in Ry, i.e. under about 5e-9 Ry, rel.), occupied-eigenvalue error, and density error in the Hartree metric E_res. Numbers fixed after E1 shows what QE actually delivers. Needs a harness and QE runs; no change of BAANDI semantics.
3. **Additive "standard" convergence profile (needs user approval).** A documented, labelled profile in which convergence is declared by the Hartree-metric E_res threshold, the energy is the output-density energy with the first-order correction (descf), Davidson uses an energy-matched test, the fixed-point certificate and gap check run only on request, and optionally buffer bands are off. The strict default is unchanged; no tolerance is relaxed silently; every headline number states the profile. A "standard" result counts as equal-accuracy only if option 2 shows its delivered error is no worse than QE's at conv_thr 1e-8.
4. **Cross-code energy match.** Not usable as a gate: BAANDI and QE totals differ by 1.28e-6 Ha/atom (existing comparisons) for convention reasons; each code is compared with itself.

### b.3 QE best configuration (each to be measured on dev08 and stored as a reference bundle with command, version, input and UPF SHA-256)
Real Gamma path (`K_POINTS gamma`) for Gamma-only; default symmetry for k meshes; default `nbnd` (occupied bands) for insulators, with the buffer-band variant as second reference; `diago_david_ndim` 2 and 4; `rmm-davidson` for many bands; FFTW3 build if one exists (QE's MKL-DFTI path does not prune); MPI layout (pools versus plane-wave) and OpenMP at equal core count; mixing beta default 0.7; conv_thr chosen by option 2. Build flags: asking the remote session which `-O`/`-march` the QE build used; likewise for BAANDI (no `[profile.release]` and no target-cpu flags in the repository today).

## (c) Staged roadmap (dependency order), recomputed 2026-11-06

Reference: QE best (real Gamma path, rel., part 6) 0.366 s at conv 1e-8, QE complex 0.684 s. Start: BAANDI default pure Rust build on dev08, si8 Gamma N=1: 5.423 s (M, part 7 build (a)). Every stage factor carries its label; E factors are not independent and later stages are optimistic about Amdahl. Standard-profile factor transferred from Mac (O) to dev08 is an assumption (E), supported by the work factor 0.33.

| Stage | Content | Factor (label) | Time (s) | Ratio to QE best 0.366 s | Ratio to QE complex |
|---|---|---|---:|---:|---:|
| start | default pure build | 1 (M) | 5.423 | 14.8 (M) | 7.9 |
| S1a | native build flags (done as measurement; needs a native artifact decision) | 0.760 (M, bit-identical) | 4.122 | 11.3 (M) | 6.0 |
| S1b | band-blocked nonlocal (P15B ticket), cheap energy and kernel caches | 0.90 (0.85-0.95, E) | 3.71 | 10.1 | 5.4 |
| S2 | locking + one rectangle (P24A) | 1.0 (0.97-1.02; O: MACs -6%, H applications +24%, failed gates; not enabled) | 3.71 | 10.1 | 5.4 |
| S3 | standard profile P17A (uncertified) | 0.317 (O wall; work 0.33; range 0.28-0.40, E on dev08) | 1.18 | 3.2 (2.8-4.0) | 1.7 |
| S4a | real Gamma path | 0.48 (0.43-0.58, E; QE's own gain 1.87x measured) | 0.57 | 1.55 (1.4-1.9) | 0.83 |
| S4c | remaining H-side items (pure FFT/nonlocal parity with QE per call) | 0.95 (E) | 0.54 | 1.5 | 0.79 |

Range at the end for si8 N=1: optimistic 4.122 x 0.85 x 0.28 x 0.43 x 0.92 = 0.39 s (ratio 1.06); mid 0.54 s (1.5); pessimistic 4.122 x 0.95 x 0.40 x 0.58 x 1.0 = 0.91 s (2.5). So parity on Gamma at one core is now inside the plausible range (it was not before part 6/7 and P17A), exceeding is unproven.

k444s (N=1, pure): default 6.260 s (M) -> native 5.093 s (M, 1.229x) -> standard 0.57 (O) 2.90 s -> k-point symmetry 0.33 (O: 1.106/3.328) 0.96 s -> S1b 0.90: 0.87 s. QE best with symmetry is not measured (est. 1.204 x 10/32 = 0.38 s plus overhead, 0.38-0.45 s): ratio 2.0-2.3 (E); optimistic 5.093 x 0.5 x 0.30 x 0.85 = 0.65 s (1.5); pessimistic 5.093 x 0.65 x 0.36 = 1.19 s (2.8). vs QE complex nosym 1.204 s: 0.72. P15C spectral guard if enabled adds 23% H on k444s.

Multi-core (si8): native + LTO pure T=8 2.756 s (M); QE real N=8 0.144 s and N=16 0.160 s are below QE's 0.3 s timer limit (not valid), so equal-core claims need larger cases (see (b)). With S1b x S3 x S4a (0.9 x 0.32 x 0.48 = 0.14) the T=8 time becomes about 0.38 s (E) before any scaling gain, i.e. 2.6 x QE's 0.144 s: unmeasurable until a larger case exists.

Order, dependencies, alternatives (changed from the previous version):
1. **S3 is now the largest single stage and is already implemented** (P17A); what remains is the dev08 comparison against QE at its measured delivered accuracy (ADR 0035 requires it before any claim). Keep `certified:false` visible.
2. **S2 shrinks.** Locking is optional and weak (O: -6% / -12.5% MACs, more H work, failed gates); the one-rectangle component alone raised si8 SCF 15 -> 16 in the ablation. **B13 (non-orthonormal basis) loses most of its premise:** its main benefit was the scalar sequential MGS (4-26% of dev08 T=1 wall) and native flags already cut `subtract_projection` from 0.822 s to 0.224 s (-0.60 s, 46% of the saving); P17A cut tail steps by 5x (345 -> 69 steps), the other target of locking; the V model over-predicted locking by 5x (30-40% predicted vs 6-12.5% measured). Revised B13 estimate: 5-12% of wall (E, was 25-30% of subspace work), high risk (skew fixtures), priority lowered.
3. **S4a real Gamma is now the largest remaining factor for Gamma cases** (QE gets 1.3-1.9x, measured, with the h_psi call count unchanged; mechanism per-call cost). After S3 the H application and phases are a larger share of the run, so the real path matters more, not less; it should start now (the dependency on S2 structure is removed because S2 is not coming). Real Gamma versus the k-point path remain separate modules; symmetry (P19A, done, default off) applies only to the k path.
4. S1a needs a user decision for a native artifact (bit-identical on dev08); S1b is the P15B ticket.
5. MKL pruned FFT (old S5) drops in priority: the pure native build beats the default MKL build on si8; keep for the MKL build only if M0 shows MKL wins on k-point cases (k444s MKL 1.05-1.08x vs pure native 1.23x; native pure T=1 5.093 vs MKL default 5.050 s, i.e. parity).
6. Standard-profile follow-ups: skew/slab/16-atom-mesh fixtures do not converge under standard; they limit the "default" claim and are the next accuracy work.
7. Forces/stress/relaxation (benchmark extension, see (b)): force time is not reduced by the standard profile; a new backlog item B19.

## (d) Decisions the user must make (recomputed; earlier approvals are not repeated)
Already settled: standard profile accepted (ADR 0035); symmetry implemented (P19A, default off); locking kept as an option (default off).
1. **Switch symmetry on by default for k-point runs?** Recommend yes for the benchmark comparison and the speed claims only after dev08 confirms the 10-point count and energy equality; keep `off` as the library default until then.
2. **Native CPU build artifact** (bit-identical, 1.12-1.40x pure): recommend yes, a separate "native" artifact for benchmarking and for users on known hardware; the default portable build stays.
3. **Which build counts in the comparison with QE:** recommend the native pure build (QE also uses plain -O3 with MKL dispatch, so the comparison stays honest only if both have default or both optimised flags; state it).
4. **Standard profile in headline claims:** recommend a comparison table with both columns, strict and standard, each labelled, standard shown with delivered accuracy; no claim of "faster than QE" from standard alone until the dev08 equal-accuracy run arrives.
5. **Real Gamma path (ADR, second solver family):** recommend yes, start now as an isolated module with the complex path as oracle (largest remaining Gamma factor, QE measured 1.3-1.9x).
6. **Drop B13 (non-orthonormal basis) and keep locking off:** recommend yes: park both; reopen only if the dev08 profile after native + standard shows algebra above 40% of wall.
7. **Add relaxation and stress workloads to the benchmark set:** recommend yes as a second phase (after SCF parity), with QE's own relax/stress runs as references; needs QE inputs and dev08 time.
8. **Larger benchmark cases for equal-core claims** (QE above 0.3 s at N=8 and 16): recommend yes (16-atom and 32-atom Si, MgO 40-atom already exists for BAANDI); needs QE runs.
9. **Spectral guard:** recommend leave `final` optional (default off); enable automatically only in strict certified production runs if the user wants completeness evidence (cost 4% si8, 23% k444s).
10. **Rounding-level changes standing rule and buffer bands in standard profile:** unchanged recommendations from the previous version (yes; keep 4 buffers in strict, evaluate fewer in standard).

## (e) Honest assessment (2026-11-06)

| Case class | Assessment | Arithmetic (labels as above) |
|---|---|---|
| Gamma-only single core (si8) vs QE real Gamma | **Parity now plausible**, exceeding unproven (was: out of reach) | 14.8 now (M) -> 11.3 native (M) -> about 3.2 with standard (O/E) -> about 1.5 with real Gamma (E); range 1.06-2.5. Needs the dev08 equal-accuracy run: if standard is not accepted as equal accuracy, the strict path ends at about 4.9 (4.122 x 0.9 x 0.48 = 1.78 s vs 0.366 s, E) |
| k-point cases at one core (k444s) vs QE with symmetry | Parity not reached even optimistically; within about 2x plausible | 2.0-2.3 mid (E), 1.5 optimistic, 2.8 pessimistic; QE symmetry reference not yet measured; k-point symmetry factor 0.33 (O) |
| Multi-core (8-16) single-k and k-point cases | Not assessable: QE times are below its 0.3 s timer limit on si8 and k444s | needs larger cases; BAANDI scaling si8 T=1 -> T=8 is 4.122 -> about 2.8 s (1.5x, native + LTO, M); QE's own N=1 -> 8 scaling at this size is 0.366 -> 0.144 s (2.5x, below the timer limit) |
| Large systems (20+ atoms, 100+ bands) | Unknown | no QE reference; P17A standard factors on si-20 (strict 5.97 s -> 2.95 s, 0.49) and mgo-40 (7.52 -> 4.18, 0.56) show smaller gains than si8 (0.32) |
| Forces, stress, relaxation | Behind and not yet measured against QE; will dominate after SCF speedups | forces 28.3% (si8), 10.9% (k444s) of SCF+force; stress 58.4%, 44.6% of SCF+pass (O); standard-profile forces are not validated (ADR 0035 forbids mixing), so force workloads stay on the strict path: SCF factor 1 |
| Tiny cases (gamma, k222s) | Not measurable (QE timer) | QE 0.06-0.2 s |


What would change this assessment:
- The dev08 equal-accuracy run shows QE's delivered error is much smaller than 1e-9 Ha (then standard needs a tighter threshold and the factor rises towards 0.5; ratio about 4 after real Gamma) or larger (the factor falls, parity easier).
- The Mac standard factor 0.317 does not transfer to dev08 (it is an observation on one machine and the work factor 0.33 is only a proxy); si-20 and mgo-40 factors are 0.49 and 0.56.
- Real Gamma delivers less than 1.7x on BAANDI (QE got 1.87x on si8, 1.29x on the tiny Gamma case).
- A dev08 profile of native + standard shows the remaining time dominated by a few sites (algebra, nonlocal 0.6 s, memcpy) rather than spread out.
- Failure of the standard profile on skew/slab/mesh fixtures to be repaired would restrict the claim to well-behaved insulators.
