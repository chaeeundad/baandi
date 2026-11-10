# C3: start-state extrapolation between ionic steps (design note, 2026-11-10)

Status: design only, no code. Source: speed-restructure proposal C3, `docs/performance/qe-parity-plan.md` row 7.2 (Q8). Written from the physics; no QE source was read for this note.

## Why no code in this round
The relaxation driver (`crates/dft-cli/src/relax.rs`) has a single start-state hook, `Context::evaluate(.., warm: Option<(&ScfProblem, &[f64])>, ..)`, which hands the previous converged density to `ScfSession::initialize_density_transfer` (ions) or `initialize_cell_density_transfer` (cell). An extrapolating predictor needs three things this path does not have:
1. the atomic-superposition density at the OLD and the NEW positions: the engine builds it only when `ScfOptions.start_density = Atomic`, privately inside `ScfProblem::new` (`start_density::atomic_superposition`), with no public accessor;
2. orbitals: the transfer contract explicitly carries "no interpolation, orbitals, mixer history or certificate" (dft-engine AGENTS.md, P21A/P23A); warm orbitals need a new session entry point and a gauge/subspace-alignment step;
3. history of two or three previous ionic steps in `restart.json` (`Restart.density` holds one density), i.e. a restart format version bump with its own compatibility tests.
Each is a contract change in a frozen-by-test area (density transfer, checkpoint/restart), so it is more than a modest change and gets its own ticket.

## Predictor (principle)
Write the converged density as n_scf(R) = n_at(R) + d(R), where n_at is the superposition of neutral-atom densities (moves rigidly with the atoms, carries most of the change of an ionic step) and d is the smooth bonding part.
- First order ("atomic"): n_start(R_new) = n_at(R_new) + d(R_old). The atomic part follows the atoms exactly; the bonding part is carried over.
- Second order: d_start = d_k + a (d_k - d_{k-1}) + b (d_{k-1} - d_{k-2}), with a, b from a least-squares fit of the new positions to the previous ones (R_new - R_k ~ a (R_k - R_{k-1}) + b (R_{k-1} - R_{k-2})); fall back to first order when fewer than three previous steps exist or the fit is ill-conditioned (the fitted step deviates from the actual step by more than a fixed fraction, chosen by measurement).
- Charge: n_at and every d integrate to N and 0 respectively, so the start charge is exact to rounding; the existing transfer check (relative charge 1e-8, finite, XC-feasible) stays. If the extrapolated density has negative samples beyond the feasibility rule, use the first-order (or the plain previous-density) start and record the fallback reason.
- Orbitals (optional second stage): linear extrapolation psi_start = 2 psi_k - psi_{k-1} after aligning psi_{k-1} to psi_k by the unitary that maximises their overlap (Procrustes/Loewdin on the occupied subspace, degenerate-safe), then re-orthonormalise; fed as the warm start of the first solve. Orbitals are a convergence aid only; the certificate is unchanged.

## Interface (proposed)
- Engine: `ScfProblem::atomic_superposition_density(&self) -> Result<Array3<f64>, ScfError>` (built on demand from the stored species data, same quadrature as the Atomic start); `ScfSession::initialize_extrapolated(source, density, orbitals: Option<&[ComplexMatrix]>)` with the same identity checks as density transfer plus a per-k orbital shape check; provenance string `extrapolated_first_order` / `extrapolated_second_order`.
- CLI: `[relax] start_extrapolation = "none" | "atomic" | "second_order"` (default `none`, i.e. today's `warm_start` behaviour, bit-identical), `extrapolate_orbitals = false`. `restart.json` version +1 stores up to three previous d fields and positions; older versions load with history length 1.
- results: per ionic step, `scf_iterations`, `davidson_iterations`, `start_kind`, `fallback_reason`.

## Acceptance test (pre-register before running)
Cases: `si16-k222-relax1` and `si32-gamma-relax1` extended to 6 BFGS steps from a perturbed start, plus one skewed and one slab-like synthetic cell (generality rule), strict class, T=1 and T=8 on dev08.
1. Identity of physics: the certified energy and forces of every ionic step within the strict class (energy 1e-8 Ha, forces 1e-6 Ha/Bohr against `start_extrapolation = none`); identical BFGS step count and final structure within 1e-5 Bohr.
2. Gain: median SCF iterations per ionic step (steps 2 to 6) reduced by at least 20 percent for `atomic` and at least 30 percent for `second_order`; total wall of the relaxation reduced by at least 10 percent on si32 at T=1 and T=8.
3. No regression: no ionic step with more SCF iterations than the baseline by more than 2; every fallback reported with its reason.
Unit tests: charge exactness of the extrapolated start on a moved two-species cell; second-order predictor reproduces a density that is exactly quadratic in a scalar path parameter; rejection paths (shape, negative charge, missing history).

## Risks
Large steps (BFGS line search, trust-radius resets) make extrapolation worse than the plain previous density: the fit check and fallback cover it. Symmetry-reduced runs need the extrapolated density to stay invariant (n_at is invariant for symmetric positions; d differences inherit invariance). Cell relaxation changes the grid metric: start with ions-only mode.
