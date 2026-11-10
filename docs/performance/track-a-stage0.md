# Track A stage 0: all-band direct minimization (variant A0) as a pre-converger

Date 2026-11-10. Proposal: speed-restructure proposal section 3 (Track A). Ledger row: TRACK-A0. Code: `crates/dft-engine/src/direct_min.rs`, bench and gates `crates/dft-engine/tests/direct_min.rs`. Raw outputs: `tests/reference/perf/track-a/stage0-mac-t1.txt` (unperturbed, all fixtures) and `tests/reference/perf/track-a/stage0-mac-seeds.txt` (kill-rule cells, seeds 601-603). All numbers are Mac observations (Apple Silicon, release, pure Rust, T=1); counts are exact on any machine, wall is not evidence.

## Verdict
The pre-registered kill rule fires on the Mac fixtures: the median total H columns of A0 (minimizer + hand-off recorded + hand-off certificate) exceed the strict-only total (recorded + certificate) on si32 (ratio 1.26) and on all four non-silicon fixtures (1.35 to 1.89). **Track A stops for variant A0.** The module stays in the tree, default off, as an isolated switchable research module (generality rule); no default changes.

## What was implemented
- Variant A0 (proposal 3.4): energy E[C] = sum f c.(H c) - E_H + E_xc - int n v_xc + E_ewald, which is the audit's direct functional; density from the real pair path (`PairDensityWorkspace`, first engine use of `accumulate_density_pair`), potential from the existing Hartree/XC code, G = H C from `RealHamiltonian::apply_block`; projected gradient R = G - C (C^T G); per-band TPA preconditioner (`dft_solvers::tpa_weight`, `kinetic_scale`, now public, arithmetic unchanged) followed by projection; Polak-Ribiere (non-negative) with reset every 50 iterations and on a non-descent direction; Loewdin retraction orth(C + tD); line search through E(0), E'(0) = 2 f tr(D^T G) and E(t), accepting only a strictly lower evaluated energy, otherwise shrinking; orthonormality defect checked every iteration (abort above 1e-10). Stops: gradient tolerance (default 1e-6 Ha), stagnation (gradient not below 0.9x its value 30 iterations earlier), failed line search, max iterations. None of them declares convergence.
- Start: the atomic Rayleigh-Ritz block at the starting density when the problem uses the atomic start (counted as start H columns); otherwise lowest-kinetic ORC unit vectors with a seeded, kinetically damped component on every row.
- Hand-off: one Ritz rotation of the occupied block, buffer bands from the start's unoccupied atomic Ritz vectors then lowest-kinetic unit vectors (orthogonalized twice), Hermitian completion, then `ScfSession::initialize_electronic` (new crate-private hook in session.rs: ordinary `initialize` followed by the same validation/import as `warm_start`; `warm_start` itself needs host digests and rejects NLCC, `initialize_density_transfer` drops the orbitals). The unchanged strict steps and fixed-point certificate produce the result.
- Switch: `direct_min::run_with_method(problem, ScfMethod, &DirectMinOptions, retain_last)`; `ScfMethod::Scf` (default) is exactly `run_outcome`. `MethodOutcome.start_path` records `Scf`, `DirectMinimization` or `FallbackScf(reason)`. Fallback reasons: several k points, spin, non-Gamma k, `gamma_real = off`, density symmetrization, XC density floor, non-uniform occupations (Fermi-Dirac uses a separate driver and never reaches this path). No CLI key was added.
- Gates (debug, in check.sh): default method bit-identical to `run`; minimizer energy equals the strict SCF energy within 1e-10 Ha at gradient 1e-7 (synthetic skew two-species fixture, measured 2e-15 to 1e-12 on the real fixtures); energies strictly decreasing; hand-off certified within 1e-9 Ha; minimizer energy, H counts, orbitals and density bit-identical for 1, 2 and 4 threads; ineligible problem falls back with its reason.

## Unperturbed fixtures (T=1, gradient tolerance 1e-6 Ha)
H columns are real columns applied (one H-equivalent per band = occupied band count columns). Strict = ordinary SCF from its default start (atomic where available).

| case | occ | min iter | min H cols (H-eq/band) | stop | hand-off iter, cols rec+cert | total A | strict iter, cols rec+cert | ratio A/S | dE cert vs strict (Ha) |
|---|---:|---:|---:|---|---|---:|---|---:|---:|
| syn-skew2 | 2 | 53 | 202 (101) | gradient | 2, 39+18 | 259 | 13, 171+16 | 1.39 | 5e-15 |
| syn-fcc | 4 | 141 | 1112 (278) | gradient | 5, 164+36 | 1312 | 22, 647+38 | 1.92 | -5e-14 |
| syn-lowsym6 | 6 | 120 | 1224 (204) | gradient | 3, 181+55 | 1460 | 16, 777+58 | 1.75 | -2e-14 |
| si-skew | 4 | 27 | 204 (51) | gradient | 2, 87+31 | 322 | 21, 550+64 | 0.52 | -2e-14 |
| si8 | 16 | 26 | 736 (46) | gradient | 2, 262+146 | 1144 | 11, 1134+185 | 0.87 | -4e-14 |
| si16 | 32 | 37 | 2560 (80) | line search | 2, 418+179 | 3157 | 10, 1746+207 | 1.62 | -2e-13 |
| mgo8 | 32 | 36 | 2308 (72) | line search | 5, 660+679 | 3647 | 15, 2053+293 | 1.56 | 5e-12 |
| si32 | 64 | 38 | 4992 (78) | line search | 2, 581+438 | 6011 | 10, 3540+394 | 1.53 | 5e-13 |
| si-elongated | 8 | 35 | 536 (67) | stagnated | not converged (100) | - | not converged (100) | - | - |
| si-slab | 8 | 32 | 496 (62) | gradient | not converged (100) | - | not converged (100) | - | - |

si-elongated and si-slab do not converge on either path (known suite limitation, benchmark-suite.md); the minimizer's fixed-occupation state lies 0.12 and 0.03 Ha above the SCF's last energy there, i.e. it is not the aufbau state (see caveats).

## Kill-rule cells, perturbed seeds 601-603 (T=1)
| case | total A per seed | median A | strict total per seed | median S | median ratio |
|---|---|---:|---|---:|---:|
| si32 | 6556, 10150, 6109 | 6556 | 5223, 5184, 5204 | 5204 | 1.26 |
| si16 | 3474, 3395, 2767 | 3395 | 2163, 2261, 2268 | 2261 | 1.50 |
| syn-skew2 | 250, 253, 252 | 252 | 187, 187, 187 | 187 | 1.35 |
| syn-fcc | 1318, 1605, 1335 | 1335 | 707, 800, 703 | 707 | 1.89 |
| syn-lowsym6 | 1432, 1442, 1440 | 1440 | 835, 835, 835 | 835 | 1.72 |
| mgo8 | 3165, 3257, 3382 | 3257 | 2346, 2346, 2346 | 2346 | 1.39 |

Limits: Mac only, three seeds instead of ten, pure Rust only, si32 is the only si32-like cell; the rule is applied as registered to these.

## Where the work goes (mechanism)
- The minimizer alone costs 46 to 80 H-equivalents per occupied band on real silicon and MgO (proposal estimate 25 to 60): A0 spends about 2 energy evaluations per iteration (each a full-width apply) and needs 26 to 40 iterations. On si32 the minimizer alone (4992 columns) already exceeds the whole strict SCF (3934 unperturbed). The strict SCF needs about 60 H columns per occupied band including buffer bands and certificate.
- Synthetic fixtures have no atomic orbitals; from the free-electron start the minimizer needs 120 to 140 iterations (linear convergence about 0.9 per iteration).
- Line searches reach the energy rounding floor near 1e-6 Ha gradient on 32 and 64 bands (si16, si32, mgo8 stop as "line search failed", gradient 1.3e-6 to 4.6e-6); every failed search costs up to 7 evaluations.
- Hand-off: when the minimizer state is accurate, the strict session needs 2 iterations (the first at tolerance 1e-2 returns the warm vectors almost unchanged, the second goes straight to 1e-9 plus one certificate). When the minimizer density error is above the certificate's 1e-7 window, the certificate fails repeatedly (si16 at gradient 1e-5: 5 iterations and 1364 to 1438 certificate columns; mgo8 at 1e-6: 4 to 5 iterations, 480 to 680 certificate columns; si32 seed 602: 6 iterations, 2582 certificate columns).
- Gradient tolerance sweep (1e-5, 1e-6, 1e-7 on six fixtures, Si and synthetic): 1e-5 makes hand-off certificates fail repeatedly (si8 9 iterations, si16 5), 1e-7 hits the rounding floor; 1e-6 was chosen as default from all six fixtures, not silicon alone.

## Observations that do not change the verdict
- Wall (Mac, T=1, observation only): si32 seeds 601/602/603 minimizer + hand-off 17.8 / 33.8 / 20.6 s against strict 29.4 / 29.1 / 27.5 s; unperturbed si32 15.0 s against 16.7 s. The per-iteration dense algebra of A0 is much cheaper than Davidson's, so wall can be lower with more H work. The registered rule is on counts and stands.
- Cross-path energy: certified energies of the two paths differ by up to 3.6e-11 Ha (si16 seed 601), 8.8e-12 (mgo8), 2.6e-12 (si32 seed 603): both satisfy the strict certificate, but stage 1's 1e-12 Ha cross-path gate would fail on these cells.
- Determinism: minimizer energy bits, H counts and densities identical at T=1, 4 and 8 on all six unperturbed fixtures tested.

## Caveats (documented, not defects)
- Several fixed points: the minimizer descends to a local minimum of the fixed-occupation energy, the mixed SCF can stop at another self-consistent state (docs/validation/scf-fixed-point-multiplicity.md). On small-gap or metallic-like cells (si-elongated, si-slab) the minimizer converges to a non-aufbau state; starting the strict session from such a state either fails to converge or (observed on syn-fcc with an unperturbed unit-vector start, before the seeded start was added) trips the divergence guard when the SCF leaves the state (density residual 7e-7 to 3e-2). `start_path` records which path produced the start.
- A pure unit-vector start can stay inside one symmetry sector (syn-fcc reached a state 0.049 Ha above the certified SCF one); the seeded start component removes this.
- The Hartree/XC evaluation and `set_local_potential` are the engine's existing allocating APIs (as in every SCF step); all orbital-sized buffers are allocated before the loop.

## What would change the premise
A variant that needs about one apply per iteration (A1/A2 in the proposal) and a hand-off that does not re-solve the buffer bands from scratch would have to halve the minimizer cost to pass the count rule on si32; that is a new registration, not a continuation of this one.
