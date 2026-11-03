# Electronic-state landscape: motivation, formulation, limits, hypotheses

Status: design document, no implementation. Directive record: [esle-directive.md](esle-directive.md). Facts about the current code are marked "verified" with a file reference; everything else is a proposal or a hypothesis.

## 1. Motivation
A Kohn-Sham SCF run returns one self-consistent density. The Kohn-Sham map `F` (input density -> output density through V_eff[n], diagonalization and occupation) can have several fixed points `n = F[n]`. Which one a run reaches depends on the starting density, the mixer and (when the physics allows it) spin and occupation freedom. Reusing a converged state as the next start is cheap but can keep every later calculation on the same metastable branch. ESLE treats the set of fixed points as the object of study instead of a single number.

Two objectives, kept apart:
1. Candidate search at fixed structure: collect distinct converged solutions, report the lowest energy found. Never a proof of the global minimum.
2. Tracking: follow a branch `n_alpha(r; lambda)` along a path `R(lambda)`. Continuity of the SCF is not evidence that the branch is the lowest.

## 2. Formulation
Fixed model M = (cell, atoms, pseudopotentials with hashes, XC functional, cutoff, FFT grid, k set with weights, electron count, occupation policy, spin treatment, tolerances). Within one M:
- Solution set `S(M) = { n : F_M[n] = n }`. A solution is stable if the SCF Jacobian `dF/dn - 1` has no eigenvalue of the wrong sign (energy Hessian positive); saddle points are also fixed points and can be reached only by special starts. Whether a converged SCF is a minimum is a separate check (stability analysis) that the current engine does not do.
- Each solution has an energy `E_alpha = E[n_alpha]` that is only comparable to another solution's energy inside the same M. Energies from different structures, cutoffs, k sets or functionals are never ranked against each other as "better states".
- Gauge freedom: orbitals are defined up to a unitary within the occupied subspace at each k (and a phase). Physical content is the density and the occupied-subspace projector `P_k = sum_occ |psi_nk><psi_nk|`.

## 3. What the current engine can and cannot show (verified in code)
- Physics: spin-unpolarized, norm-conserving, LDA-PZ, Bloch k points, insulating occupations only; partial occupations and smearing are rejected (`crates/dft-engine/src/scf.rs`, `OccupationPolicy`, `PartialOccupations`). No spin, no DFT+U.
- SCF always starts from a uniform density `N_el / Omega` (`ScfProblem::uniform_density`, used in `run_impl`). The single-step functions `step` and `step_kpoints` accept an input density and optional initial orbitals, but the loop that drives them is internal.
- The UPF parser checks `PP_RHOATOM` values but exposes no atomic-density field; an atomic-superposition start therefore needs a parser extension first.
- Results: `ScfResult` holds energy audit, eigenpairs per k, density, history; the CLI writes `results.json`, `density.bin`, `density.json`. There is no checkpoint, restart or state serialization and no wavefunction file.
Consequence: for the present physics (an insulating Si-like cell) we expect essentially one solution. Multi-start there is mainly a test of infrastructure (the same state reached from different starts, reproducibility, tolerance sensitivity), not a source of competing minima. Competing solutions need spin polarization, occupation freedom (metals with smearing), or correlation corrections.

## 4. State distance (proposal; tolerances measured later)
All metrics need an identical model M except where stated; incompatible pairs are reported "not comparable", never forced.
- Density: `d_n = sqrt(integral (n_A - n_B)^2 d^3r) / N_el`, same definition as the existing SCF residual (`density_distance`), unit Bohr^(-3/2) per electron. Same cell and FFT grid required; otherwise compare through G-space on the common G sphere.
- Occupied subspace (gauge invariant): per k, `S_k = C_A^dagger C_B` over occupied bands in the common plane-wave basis; `d_sub(k) = sqrt(1 - ||S_k||_F^2 / N_occ)` (RMS sine of principal angles, dimensionless, in [0,1]); aggregate `max_k` and the weighted mean. Invariant to unitary rotations and phases; cost O(N_pw N_occ^2) per k.
- Spectrum: per k, occupied (and a stated number of empty) eigenvalues, maximum absolute difference in Hartree; reported, not used as a sole identity criterion (band crossings reorder).
- Energy: `|E_A - E_B|` in Hartree, informational. Similar energies never merge states.
- Local moments, spin density, occupation matrices: future fields (need spin / DFT+U).
Classification is three-way: SAME if all active metrics are below `c_same * noise`, DISTINCT if any is above `c_diff * noise`, UNRESOLVED in between. `noise` is measured, per metric, from repeated runs of the same model with different seeds, thread counts and tolerance schedules (Phase B reproducibility tests, starting from the floors recorded in `docs/validation.md`: density symmetry spread 2.1e-11 ours, SCF tolerance sensitivity 1e-12 Ha). The factors `c_same < c_diff` are chosen and documented from that data, not assumed here.
False merge: a metric that is blind to a real difference (energy-only, density-only for states differing in a localized occupation). Mitigation: require all applicable metrics, keep UNRESOLVED. False split: noise from loose tolerances or different SCF paths; mitigation: tight tolerance recomputation before declaring DISTINCT.

## 5. Branch tracking (proposal)
For a path lambda_0..lambda_n, a branch is a chain of states, each solved from a transferred predecessor. Records: parent link, state distance to the parent after transfer and after convergence, energy and density-localization measures versus lambda, convergence failures, and the comparison against an independent fresh solve at the same lambda (the only evidence of inherited-state bias). Branch crossing is flagged when two chains' states become SAME or swap energy order; disappearance is flagged as a non-convergence event, never silently re-seeded.

## 6. Hypotheses (falsifiable; outcomes reported either way)
H1. In the present insulating model, independent starts (uniform, perturbed, previous) converge to states with all metrics within noise; i.e. one solution. (Falsified if a distinct converged state appears.)
H2. Density transfer from a neighbouring structure reduces SCF iterations versus a uniform start for small distortions, with identical final state within noise.
H3. Continuation through a structural change reaches the same state as a fresh solve in the insulating model; differences first appear only once spin or occupation freedom exists.
H4. The cost of the extra independent solves per accepted distinct solution is measurable and can be budgeted; an adaptive policy is only worth building if a simple fixed multi-start misses solutions in a benchmark where solutions exist.

## 7. Limits stated up front
No finite search proves global optimality. Lower energy found by another start proves the earlier candidate was not the lowest known. Stability (minimum versus saddle) is not assessed by SCF convergence. Results depend on the pseudopotential: the validation UPF has unverified provenance (`docs/mvp-review.md`). Learned strategies, if ever built, must beat deterministic baselines (Phase F).
