# Where the non-diagonalization time goes: si8 real Gamma + standard profile (Mac observation)

Written 2026-11-07 by the performance owner. Measurement only; no source was changed. Everything here is a Mac observation (shared machine, load average about 1.4 during the runs), not an acceptance number. Binary: `target/release/dft-cli` (sha256 prefix 558bd73dd6c3345d), HEAD b59d10b, UPF `Si.pz-vbc.UPF`.

## Method

- Wall and phase data: the CLI report JSON (`timings.setup_s`, `engine_setup_s`, `engine_total_s`, `total_s`, `phases`), 30 runs per case through the case inputs of `scripts/bench_scf.py` (`--threads 1`; `gamma_real = "on"`, `convergence_profile = "standard"`; k444s with `[symmetry] mode = "auto"`). Spread over 30 runs is under 2% for every number below.
- Site attribution: macOS `sample` (1 ms) of the release binary, one sample file per run, 30 runs per case (si8: about 14.9k stack samples, k444s: about 11.4k), aggregated by inclusive frame counts. The first attempt attached `sample -wait` after launch and lost the first 0.1 to 0.15 s (it under-counted setup: 13.7% against the clock's 21%). The final si8 run starts the process suspended, attaches `sample`, then continues it; setup then reads 22.2% against 20.9% by clock. The tool still misses about 15 to 20% of the wall time overall, so fractions are shares of sampled time; seconds below are share times the JSON clock total and are labelled "(s)" as estimates. Aggregation script and raw files are in the session scratchpad, not in the repository.
- Code read: `crates/dft-engine/src/scf.rs` (`construct`, `step_with_occupations`), `energy.rs`, `session.rs`, `dft-physics/src/hamiltonian.rs`, `hartree.rs`, `dft-upf/src/radial.rs`.

## 1. si8 real Gamma standard: where the 0.146 s goes (30 runs, medians)

| Item | Seconds | Share of wall |
|---|---:|---:|
| Wall (process) | 0.6148 | 100% |
| CLI total (`total_s`) | 0.6111 | 99.4% |
| engine total | 0.6079 | 98.9% |
| engine setup (`ScfProblem::construct`, `engine_setup_s`) | 0.1269 | 20.6% |
| five phases (diag 0.3385, energy 0.0535, density 0.0469, potential 0.0173, mixing 0.0136) | 0.4698 | 76.4% |
| engine remainder (iteration glue, history records, mixer-side clones; engine_total - setup - phases) | about 0.011 | 1.8% |
| CLI before the engine (`setup_s - engine_setup_s`) | 0.0024 | 0.4% |
| CLI after the engine plus process start and output (wall - engine_total - 0.0024) | about 0.005 | 0.8% |

So the 0.146 s are: **0.127 s in engine setup (87%)**, about 0.011 s engine glue, about 0.008 s process and output. There is no first-iteration or final-step cost outside the phases (the standard profile has no certificate; JSON writing is about 6 of 14.9k samples; `WorkerPool::new` at T=1 and FFT plan creation do not appear in the samples).

### Engine setup breakdown (30 runs, 3319 samples inside `construct`)

| Site | Share of construct | Seconds (est.) | Share of wall |
|---|---:|---:|---:|
| `hamiltonian::ionic_local` (local pseudopotential on the density sphere) | 86% (2853 samples) | about 0.109 | about 18% |
| of which `erf` evaluated inside the per-G radial integral | 42% of construct (1388 samples) | about 0.053 | about 8.6% |
| `NonlocalOperator::new` (projector radial transforms per G, structure phases) | 10% (338) | about 0.013 | 2.1% |
| `RealNonlocal::new` | 2.7% (90) | about 0.003 | 0.5% |
| Ewald energy | 1.2% (40) | about 0.002 | 0.3% |
| `Hamiltonian::new`, band density, FFT plans, remainder | about 1% | about 0.001 | 0.2% |

Mechanism (read from `dft-upf/src/radial.rs`, `dft-physics/src/hamiltonian.rs`): `ionic_local` calls `LocalTransform::at_g` once for every G of the density sphere (23,871 vectors for si8 at 30 Ry, density cutoff 4x), and each call integrates over the whole radial mesh (hundreds of points) with one `erf` and one spherical Bessel (sin) per point. The `erf(r)` term does not depend on G at all, and V(|G|) depends only on |G|: the sphere has only about 525 distinct |G|^2 values (computed offline for the cubic cell with the same floating-point order; 661 shells by Miller class), a factor of 45. Cost therefore scales with cell volume (and is independent of atom count): the primitive si2 Gamma case has engine setup 0.031 s of 0.067 s wall (46%, 30 runs), 4.1x smaller than si8 at 4x the volume, consistent with the volume law.

At T=8 on the Mac (10 runs, directional) setup is unchanged at 0.127 s (single k, serial) while wall is 0.570 s; serial non-diagonalization work (setup 0.127 + energy 0.053 + potential 0.017 + mixing 0.014) is 0.21 s, 37% of the T=8 wall. On dev08 it will be a larger share if real-Gamma diagonalization scales (not measured; Mac scaling is not evidence for dev08).

## 2. The energy phase (0.0535 s, 8.7% of si8 wall)

Inclusive shares inside `energy::evaluate_with_core_xc` (30 runs, 1275 samples; one k point):

| Part | Share of the phase | Seconds (est.) | Note |
|---|---:|---:|---|
| `<c\|V_NL\|c>`: `NonlocalOperator::apply` per orbital (complex path, 3 passes: projections, weights, y accumulate, then dot) | 51% | about 0.027 | allocates a scratch vector per call |
| Hartree energy of the output density (`hartree::evaluate`: forward FFT, kernel loop, inverse FFT, alloc of squared frequencies) | 21% | about 0.011 | the inverse FFT and the potential are never read by the audit (only `energy`) |
| XC energy (`xc::evaluate`, LDA-PZ) | 12% | about 0.0066 | needed |
| `local_energy` (another forward FFT of the same density) | 9% | about 0.0048 | recomputes the transform Hartree already did |
| kinetic, charge, double counting, Bregman remainders, sums | about 7% | about 0.004 | |
| Ewald | 0 | 0 | computed once in `construct` |

The audit is part of the convergence logic only in strict mode: `converged()` (scf.rs) needs `energy_change`, from `total_hartree()` = kinetic + local + nonlocal + Hartree + XC + Ewald. The other audit fields (band form, Harris-Foulkes, double counting, functional-difference bound) are reported but do not enter the stopping test. In the standard profile (`session.rs`, converged = `standard_error_ha < threshold`) the energy does not decide convergence either; the audit then only fills the history record (per-iteration energies and `energy_change`), the final reported energy, and the non-finite check (`step.energy.total_hartree().is_finite()`).

Potential phase (0.0173 s): Hartree potential (2 FFTs), XC potential, effective potential build and a clone per state. Density phase (0.0469 s): 95% in `band_density::transform` / fused `accumulate_density`, 63% of it rustfft butterflies (one inverse 3-D FFT per band, 16 bands); the real-Gamma pair trick (R2) is used in the Hamiltonian apply (`real_gamma.rs`) but not in the density, which still takes the complex orbitals one band at a time.

## 3. k444s with symmetry (10 irreducible k of 64, 9 SCF iterations, 30 runs)

Wall 0.539 s; engine setup 0.067 s (12%, `construct` is only 2% of samples here, 224 of 11.4k, plus `enable_symmetry` 1.1%); CLI before the engine 0.013 s (symmetry detection and mesh reduction, directional); diag 0.296, energy 0.098, density 0.039, potential 0.0125, mixing 0.004.

**The energy phase is 18% of the wall (20.7% of samples) and most of it is repeated work.** The audit runs once per k point (`ordered_chunks` over states in `step_with_occupations`), and each call recomputes terms that do not depend on k: Hartree energy of the output density (43% of the audit samples), XC energy (25%), `local_energy` (18%), charge and the double-counting loops. These are the same numbers for every k (same density, same input potentials); `weighted_audits` then averages them with weights summing to one. Only kinetic, nonlocal, band sum and the eigen-residual bound are per k (about 10% of the audit here). So about 86% of 0.098 s = 0.084 s is shared work done 10 times, 0.076 s redundant, about 14% of the wall at T=1 (sample counts 2361 in the audit, 30 runs).

si16-gamma was not profiled: the stored suite records eigensolver failures for it (`benchmark-suite.md`), so it does not run to convergence.

## 4. Candidates for the non-diagonalization share, ranked by saving / risk

Ledger check: none of these is a ledger row. Closest relations: A43 (blocked nonlocal traversal, rejected, different mechanism: this is projection-only energy, not a traversal change), B18 per-shell caching of radial transforms ("exact, only symmetric k": extended here to the local potential at Gamma with new measured evidence: 86% of construct), backlog B07 (cheap energy in non-final iterations, S-7), A10 (parallel audit, adopted; this removes redundancy rather than parallelizing it). Not revived: lane batching, band-grouped FFT, hybrid RMM, nothing from the do-not-repeat list.

| Rank | Candidate | Saving (si8 real standard, T=1, Mac est.) | Other cases | Result change | Risk |
|---|---|---|---|---|---|
| 1 | **C1: precompute the G-independent part of the local-potential radial integrand (r^2 V + Z r erf(r) and rab) once per species, and memoize `at_g` by the exact bit pattern of |G|^2 within one `ionic_local` call** | about 0.10 s of 0.615 s (16%) | si2 Gamma: about 0.027 s of 0.067 s (40%); any cell size: the cost is N_G x mesh independent of atoms, +/-G pairs alone give 2x, cubic/high symmetry up to 45x; k meshes pay once, so about 2-10% | none: same operations in the same order on the same values; bit-identical energy and density | very low |
| 2 | **C2: compute the k-independent audit terms (Hartree, XC, local energy, charge, double counting, remainders) once per iteration and reuse them in every k audit** | none at Gamma; k444s (10 k): about 0.076 s of 0.539 s (14%); k222s/k222 proportionally | grows with k count; at T>1 it cuts CPU, wall gain smaller because audits ran in parallel | none if the shared value is computed by the same function and then weighted exactly as now (each k still multiplies the same f64 by its weight) | low to medium (touches `scf.rs` audit loop and `energy.rs`; ownership of both files) |
| 3 | **C3: energy-only Hartree in the audit and reuse of its forward FFT for `local_energy`** (drop the unused inverse FFT, the unused potential array and the per-call squared-frequency allocation; 3 transforms become 1) | about 0.008 s (1.3%) | per audit, so also helps C2's single computation | none: the forward transform of the same real array is identical | low |
| 4 | **C4: nonlocal energy from projections only** (the energy needs sum_f <c\|beta_i> D_ij <beta_j\|c>: skip the weight-times-column accumulation and the dot with c, optionally as one GEMM) | about 0.018 s (3%) with the same loops, up to about 0.025 s with a GEMM | k444s: nonlocal is 4.5% of its audit, small | rounding change at about 1e-16 relative in E_NL (energy within 1e-12 Ha, density bit-identical, SCF counts to be checked) | low to medium (must be declared a rounding change) |
| 5 | C5: bit-exact |q| memo and shared Bessel evaluation in `NonlocalOperator::new` | about 0.012 s (2%) at Gamma | near zero at generic k (k+G norms rarely repeat) | none | low |
| 6 | C6: pair density (two real bands per complex FFT) for gamma_real, extending R2 to `BandDensity` | about 0.012 to 0.016 s (2-3%; estimate from the FFT share of the density phase, directional) | Gamma only | rounding change (density no longer bit-identical to the complex path); real-Gamma acceptance pending | medium |
| 7 | C7: skip the audit in non-converged iterations of the standard profile | 0.054 s (8.7%) si8, 0.098 s (18%) k444s | standard profile only | **changes the report**: per-iteration energies and `energy_change` in the history disappear, and the finite-energy divergence check would move; convergence itself is unaffected (it uses `standard_error_ha`); strict mode needs the energy every iteration | needs a decision; also needs `step_with_occupations` split so the audit runs after the convergence test |

Not worth doing: thread-pool startup (not visible in the samples), structure-factor reuse (about 190k sin/cos in `ionic_local`, about 3-5 ms, 0.5%; becomes the dominant remainder of `ionic_local` after C1 only), CLI JSON output (6 of 14.9k samples).

Combined C1 to C3 and C5 (all bit-identical) would take about 0.13 s off 0.615 s on si8 (about 21%) and about 0.04 of 0.067 s (60%) on si2 Gamma; with C2 the k444s wall drops by about 14%. A rough bound only; effects were not tested together and none was implemented.

### Cautions

- C1 must keep the zero-G branch (`alpha_z`) separate and key on exact bit patterns, not tolerances; the same pattern probably exists in the NLCC core density (`spherical_density_transform`, not profiled here since Si.pz-vbc has no core correction) and should be checked for fixtures that have one.
- C2 changes only how often a function is called; keep the per-k weighting expression so `weighted_audits` yields the identical bits, and keep the symmetric (`symmetry.apply`) density as the input.
- Acceptance: three-point comparison and `--require-identical` on gamma, k222s, k222, k444s, si8 at T=1,2,4,8,12 for C1 to C3, C5; energy within 1e-12 Ha and density identical across T for C4, C6; dev08 numbers through the orchestrator (setup fractions there may differ because the machine and libm differ).

## Limits

Mac only, shared machine, load about 1.4; `sample` under-reports total time by about 15-20% and was attached from process start only in the final si8 run (k444s, 30 runs, used `sample -wait`, so its setup share is understated; its energy-phase shares are unaffected). Seconds from sample fractions are estimates. The 13 ms CLI pre-engine time in k444s is only partly explained (`enable_symmetry` about 5 ms). Distinct-|G| count (525) was computed offline for the cubic cell and the Rust order of floating-point operations was assumed equal; the real count must be measured by the implementer (a counter, not in this study).
