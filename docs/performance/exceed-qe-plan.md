# Plan: match and then exceed QE 7.6 at equal cores and equal accuracy

Written 2026-11-05 by the performance owner. Target (AGENTS.md, Performance target): BAANDI must match and then exceed QE 7.6 on the benchmark set at equal core count and equal accuracy, against QE's BEST configuration per case. Sources: relayed dev08 data (`tests/reference/perf/dev08_final_relayed.md` parts 4 and 5, `qe_iteration_data_relayed.md`), the P14A/P14B tables, the four QE study notes G, H, S, V of 2026-11-05 (see [backlog.md](backlog.md) for the key), the ledger. Every number marked (est.) is an estimate, (der.) is arithmetic on recorded numbers, (rel.) is relayed from dev08. Nothing here is a new measurement. Items needing the remote session are collected as M0 in the backlog.

## (a) Gap model

### a.1 si8 Gamma, one core, conv_thr 1e-8 (fdce016, rel.)
BAANDI pure 5.801 s, MKL 6.016 s; QE internal 0.684 s (complex Gamma, MKL-DFTI, nbnd 20, nosym). Ratio 8.48 pure, 8.80 MKL. With P14B (rel.): pure 5.463 s (7.99), MKL 4.814 s (7.04).

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
| QE configuration | the 8.48 compares against QE's complex Gamma, nbnd 20, MKL-DFTI (no FFT pruning, 0.29 s of 0.68 s is FFT). Best configuration: real Gamma (est. 1.7-2.4x faster, G 7.4), occupied bands only (nbnd 16, est. 1.1-1.25x), FFTW3 build prunes (saves about 0.13 s of the FFT part, est.) | QE best on si8 about 0.22-0.40 s (est.); ratio vs best with P14B pure 14-25 (der.) |

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
(1) Hamiltonian application is not the main gap (V and H agree: 3.7 ms vs 3.5-4.5 ms per step). (2) Davidson algebra is: about 70-78% of BAANDI's time vs 16% of QE's. (3) Step and vector counts are about 2.2-3.5x and partly definition. (4) QE's best configuration lowers the target by about 2x (real Gamma) on Gamma-only and about 3x (symmetry) on k meshes; our current ratios understate the gap. (5) Caveat on the model: dev08 pure si8 T=1 gained only 1.06x from a 56% cut in explicit GEMM flops (A55) while MKL gained 1.26x; the flop model over-predicts time gains for pure Rust until the post-P14B profile is known.

## (b) What the target means operationally

### b.1 Benchmark set
Protocol cases gamma, k222s, k222, k444s, si8 (QE references exist, complex Gamma, nosym); the 19-case suite of ticket P10G for generality (BAANDI only today); additional cases still missing and required for any "exceed" claim: a 16-atom Gamma case and the 32/64-band supercells (`super4`, `super8`), a 128-or-more-band case, a slab and an elongated cell with QE references, and Gamma-only cases above 0.2 s QE time (timer resolution 0.01 s; entries below about 0.3 s are excluded). Layouts: 1, 4, 8, 16 cores; QE with the same number of MPI processes (or the best of pools, plane-wave and OpenMP at that core count).

### b.2 Equal accuracy: options
1. **Status quo (equal conv_thr label).** QE conv_thr 1e-8 vs BAANDI strict defaults with certificate. Simple, but compares different guarantees; BAANDI is penalised by its certificate and buffer bands. Already approved; keep as the "strict" column.
2. **Delivered-accuracy comparison (recommended).** For each case run each code to its own tight reference (QE 1e-12, BAANDI strict) and require, for the timed run, the same bounds on: total-energy error against that code's own reference (QE at 1e-8 is already below the printed 8 digits in Ry, i.e. under about 5e-9 Ry, rel.), occupied-eigenvalue error, and density error in the Hartree metric E_res. Numbers fixed after E1 shows what QE actually delivers. Needs a harness and QE runs; no change of BAANDI semantics.
3. **Additive "standard" convergence profile (needs user approval).** A documented, labelled profile in which convergence is declared by the Hartree-metric E_res threshold, the energy is the output-density energy with the first-order correction (descf), Davidson uses an energy-matched test, the fixed-point certificate and gap check run only on request, and optionally buffer bands are off. The strict default is unchanged; no tolerance is relaxed silently; every headline number states the profile. A "standard" result counts as equal-accuracy only if option 2 shows its delivered error is no worse than QE's at conv_thr 1e-8.
4. **Cross-code energy match.** Not usable as a gate: BAANDI and QE totals differ by 1.28e-6 Ha/atom (existing comparisons) for convention reasons; each code is compared with itself.

### b.3 QE best configuration (each to be measured on dev08 and stored as a reference bundle with command, version, input and UPF SHA-256)
Real Gamma path (`K_POINTS gamma`) for Gamma-only; default symmetry for k meshes; default `nbnd` (occupied bands) for insulators, with the buffer-band variant as second reference; `diago_david_ndim` 2 and 4; `rmm-davidson` for many bands; FFTW3 build if one exists (QE's MKL-DFTI path does not prune); MPI layout (pools versus plane-wave) and OpenMP at equal core count; mixing beta default 0.7; conv_thr chosen by option 2. Build flags: asking the remote session which `-O`/`-march` the QE build used; likewise for BAANDI (no `[profile.release]` and no target-cpu flags in the repository today).

## (c) Staged roadmap (dependency order)

Predicted time and ratio are estimates; stage gains are multiplied for a mid case and are not independent (Amdahl lowers later gains; ranges in the text). Start state: P14B, pure Rust, si8 Gamma N=1 = 5.46 s (rel.), QE complex 0.684 s, QE best central 0.30 s (range 0.22-0.40, est.). k444s start 6.127 s, QE complex 1.204 s, QE best with symmetry central 0.41 s (est.).

| Stage | Content (backlog id) | Expected time factor (est.) | si8 N=1: ratio to QE complex / to QE best | k444s N=1: ratio to QE complex (nosym) / to QE best (sym) |
|---|---|---|---|---|
| S0 | dev08 measurements M0; no code | 1.0 | 8.0 / 18 | 5.1 / 15 |
| S1 | bit-identical and cheap: band-blocked nonlocal, release flags, cheap energy in non-final iterations, kernel caches (B07) | 0.85-0.93, mid 0.89 | 7.1 / 16 | 4.5 / 13 |
| S2 | Davidson algebra: locking + one rectangle (B01), then non-orthonormal basis (B13); capacity retune afterwards (B17) | 0.55-0.80, mid 0.67 (basis: non-H machinery 65-70% of time, cut by 30-65%) | 4.8 / 11 | 3.0 / 9 |
| S3 | step count and accuracy definition: standard profile (B14, needs approval), Hartree-metric mixing (B16), energy-weighted schedule with redo (B04), atomic starts (B10) | 0.55-0.85, mid 0.70 | 3.3 / 7.6 | 2.1 / 6.3 |
| S4a | real Gamma path (B02), Gamma-only | 0.43-0.58, mid 0.50 (QE best already includes it) | 1.7 / 3.8 | not applicable |
| S4b | k-point symmetry (B15) on the k path | 0.31-0.40 (10 of 32 k points, plus symmetrisation) | not applicable | 0.7 / 6.3 (vs complex nosym); vs best (symmetric QE) 2.1 / 2.1 |
| S5 | H per vector to QE parity or below: MKL pruned FFT, matrix-product nonlocal (B08); parallel decomposition (B03) for N>1 | 0.88-0.96 at N=1; larger gain only at N>1 | 1.55 / 3.5 | 2.0 / 2.0 |

Range of the cumulative product at the end (si8 N=1): optimistic 0.85 x 0.55 x 0.55 x 0.43 x 0.88 = 0.097 (time 0.53 s, ratio 0.78 vs complex, 1.8 vs best); pessimistic 0.93 x 0.80 x 0.85 x 0.58 x 0.96 = 0.35 (1.9 s, ratio 2.8 vs complex, 6.4 vs best). k444s N=1 after S1-S4b and S5: optimistic 0.85 x 0.55 x 0.5 x 0.31 x 0.9 = 0.065 (0.40 s, ratio 1.0 vs best); mid 2.0 vs best; pessimistic 0.93 x 0.8 x 0.8 x 0.4 x 0.96 = 0.23 (1.4 s, 3.5 vs best).

Order, dependencies, alternatives:
1. S0 first (decides everything; the post-P14B dev08 profile may reorder S1/S2).
2. S1 can run in parallel with S0 reading; it conflicts with nothing and is bit-identical except the energy change.
3. S2: locking and non-orthonormal basis are ALTERNATIVES in structure and must not be done together (V-3); order locking then non-orthonormal on the same fixtures; stage 0 of the non-orthonormal route is the rank/conditioning diagnostic on the 20-atom skewed case (A52). Capacity retune only after S2.
4. S3 is gated by the user decision on the standard profile; E0-E2 and the frozen-Hamiltonian replay come first (offline, no solver change); the Hartree-metric mixer (B16) precedes the schedule (S: mixer must accept the noise before looser tolerances). Alternative to the energy-weighted schedule: single-sweep RMM in the SCF loop (B11): conflicts with locking (uncertified intermediate maps), low priority.
5. S4a real Gamma versus the k-point path: separate modules; real Gamma must be built on the S2 solver structure (otherwise the algebra is rewritten twice and a second Davidson family doubles maintenance); S4b symmetry applies only to the k path; neither helps the other's case class.
6. S5: matrix-product nonlocal only if band-blocked nonlocal (S1) leaves a measured gap; MKL pruned FFT only for the MKL build; B03 is independent but changes rounding like S2, so it follows S2.

M0 (must be measured first on dev08): QE real Gamma and variants; post-P14B profiles; build-flag experiment; E1 (how much of the 8.5x is definition); frozen-Hamiltonian replay (criterion strength versus convergence rate); re-crossing count for locking.

## (d) Decisions the user must make (each with one recommendation)
1. **Approve an additive "standard" convergence profile** (b.2 option 3)? Recommend YES, strictly additive and labelled, strict default untouched, and no "exceed QE" claim without the option-2 accuracy evidence.
2. **Equal-accuracy definition:** recommend option 2 (delivered accuracy against each code's own reference), with numerical bounds set after E1.
3. **QE best-configuration scope:** recommend all items in b.3 (real Gamma, symmetry, nbnd default, ndim 2 and 4, FFTW3 build if available, best MPI layout) so the target means what AGENTS.md says.
4. **Remote actions:** permission for the remote session to build one extra QE (FFTW3) only if no such build exists in `~/fq/builds` (AGENTS.md forbids unrequested remote changes and bulk installs); recommend YES, one build, no installs elsewhere.
5. **k-point symmetry reduction in BAANDI** (ADR, density symmetrisation): recommend YES, after S2, because without it every k-mesh case stays about 3x behind QE's default.
6. **Real-valued Gamma path** (ADR, second solver family): recommend YES, after S2 and after the M0 real-Gamma QE number, built as an isolated switchable module with the complex path as oracle.
7. **Rounding-level changes as a standing rule** (P14B precedent: realistic fixtures keep SCF count and energy within 1e-12 Ha, converged toys keep energy tolerance and report iteration changes): recommend YES as standing rule for S1-S5, because every structural item changes rounding (fixed-point multiplicity makes counts on toys unreliable).
8. **Buffer bands in the standard profile:** recommend keep 4 extra bands in strict, evaluate zero in standard only after E2, because the 8-atom cluster failure (A17) is the reason they exist.
9. **Native-CPU build artifact** (`x86-64-v3` or dispatch) for benchmarking: recommend YES as a separate artifact, default build unchanged.
10. **Scope of "equal cores":** recommend QE allowed its best MPI/OpenMP layout at each core count, BAANDI threads only.

## (e) Honest assessment

| Case class | Assessment | Arithmetic |
|---|---|---|
| Many-k (k222, k444s) at 8-16 cores, both codes symmetric | Parity plausible, exceeding possible but unproven | k444s N=16 is 2.17 now vs QE nosym; stage factors 0.89 x 0.67 x 0.70 = 0.42 give about 0.9 vs QE nosym; with symmetry on both sides about 1-2 vs best, limited by k-parallel saturation at about 10 k points (k222-like saturation at T=8 seen); exceeding needs intra-k scaling better than today |
| Many-k at 1 core | Parity only in the optimistic case | optimistic 1.0 vs best, mid 2.0, pessimistic 3.5 (c) |
| Gamma-only small cells (si8) at 1 core vs QE real Gamma | Out of reach with current ideas; matching QE complex Gamma is plausible | best-case 0.53 s vs QE best 0.22-0.40 s; floor: QE spends 0.11 s on algebra, BAANDI 4.0-4.7 s, and S2 gives at most 1.5-3x of that; H alone at vector-count parity is about 0.18-0.25 s (real path) so total below QE's 0.30 s requires algebra plus phases below 0.1 s, a 40x cut from today, and S2 delivers 1.5-3x |
| Gamma/single-k at 8-16 cores (si8) | Out of reach | QE 0.21-0.30 s (at the timer limit) vs BAANDI T=8 4.15 s; single-core stages give about 2.4x (0.42 factor), threading 1.40x now; even 100% utilisation of 8 workers gives 4.24x on the parallel part (ideal 1.369 s): remaining gap about 4-6x vs QE complex and 8-12x vs QE real |
| Large systems (20+ atoms, 100+ bands) | Unknown; parity plausible after S2 because both codes then scale as N m^2 in the algebra and BAANDI has no MPI communication | no QE reference or large BAANDI case exists; P10G has `syn-100-bands` only for BAANDI |
| Cases below QE's timer limit (gamma, k222s) | Not measurable | QE 0.06-0.2 s, ticks up to 20% |

What would change this assessment:
- E1/E2/replay show the definition share is 2x or more of the 8.5x (then S3 may close the step gap cheaply and Gamma parity moves into reach).
- M0 shows QE real Gamma is less than 1.7x faster than complex (the best-config target drops) or that no FFTW3 QE build is competitive.
- The post-P14B dev08 profile shows the algebra concentrated in a few removable sites (scalar MGS 4-26% of wall, memcpy 0.17-0.66 s) rather than spread over the step: then S2 would exceed the flop-model 1.5-3x.
- A measured vector count at or below QE's about 950 with the cost per vector at 0.38 ms or less (MKL pruning plus nonlocal blocking): then H time falls to QE parity and exceeding becomes possible through pruned FFT (QE's MKL build does not prune) and equal-or-better parallel scaling.
- Evidence to the contrary: if locking or the non-orthonormal basis fail the skewed-cell gates (as CholeskyQR2 did), S2 shrinks to the cheap one-rectangle change and the gap stays above 5x on Gamma.
