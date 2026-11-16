# ADR 0043: One SCF path (standard profile) with a fixed-point verification after convergence

- Status (amended 2026-11-16): **Default flipped by user decision 2026-11-16; the input key `scf.convergence_profile` is removed.** The single standard path with the fixed-point verification is the only SCF path of a classic run; `BAANDI_SINGLE_PATH` is removed and the strict path is a reference tool behind `BAANDI_REFERENCE_STRICT=on`. See the amendment at the end; a second amendment (ONEPATH-1, same day) moves forces and fixed-cell relaxation onto the single path, a third (ONEPATH-2) stress and variable-cell relaxation, a fourth (ONEPATH-3) collinear spin and Fermi-Dirac occupations. The original status line follows for the record.
- Status: **Accepted in direction by the user (decisions 1 and 5 of 2026-11-14, `docs/performance/requests/work-order-beat-qe-2026-11-14.md` section 1)**; this ADR fixes the definition the work order left open (P0.1 item c). The implementation is behind `BAANDI_SINGLE_PATH` (default `off`) until the dev08 check list `docs/performance/requests/p01-single-path.md` shows convergence and delivered error in class on the generality set. Switching the default is a separate, one-line change recorded as an amendment here.
- Supersedes, when the default flips: the "strict is the default" clause of ADR 0035. ADR 0035's labelling rules stay (a standard result says so); the "never certified" rule is replaced by the verification below, which is reported under its own name and is not the strict certificate.

## Context

BAANDI has two SCF paths. The strict path (default) stops on density L2 <= 1e-7, energy change <= 1e-8 Ha and residuals <= 1e-8 Ha, then runs the fixed-point certificate: a re-solve at the output density with eigen tolerance 0.1 x eigen_tolerance (1e-9 Ha), accepting when the density changes by at most density_tolerance and no occupied eigenvalue moves by more than 1e-7 Ha. The certificate also drives the per-step tolerances (the strict ratchet works toward it). The standard path (ADR 0035, P17A) stops when the Coulomb (Hartree-metric) estimator E_H[rho_out - rho_in] falls below `threshold_ha` (5e-9 Ha default, the standard accuracy class) with the residual guard (sum of occupied f_i |r_i|^2 below the same threshold, default on); it is uncertified, so forces, stress and relaxation reject it.

The user decided on 2026-11-14 that the standard path becomes the only SCF path, that accuracy is a threshold as in QE, and that the certificate becomes a verification step after convergence, on by default, with its cost inside every timed run.

## Decision

1. **One path.** A run without an explicit `scf.convergence_profile` uses the standard profile. The threshold is `scf.standard_threshold_ha` when given, otherwise `StandardProfile::default().threshold_ha` (5e-9 Ha, the standard class today; `FLOOR_HA` 5e-14 Ha and `CAP_HA` 5e-3 Ha unchanged). Mixing defaults follow the standard class (Pulay beta 0.5), as they do today for an explicit standard run. No numerical kernel changes.
2. **Strict stays callable.** `convergence_profile = "strict"` keeps today's strict path bit for bit, including its certificate (which keeps driving the strict tolerances). All accuracy gates and references keep using strict. `convergence_profile = "standard"` written explicitly keeps today's uncertified standard behaviour unless `verify = "fixed_point"` is also given (bit-identical old inputs).
3. **Verification is on by default** on the single path; `scf.verify = "off"` skips it, `scf.verify = "fixed_point"` requests it (also on an explicit standard profile). `verify` together with an explicit strict profile is an input error (strict always runs its certificate). The verification does not drive any per-step tolerance.
4. **What the verification establishes** (the definition P0.1 item c asked for). Let rho_out be the output density of the step at which the standard stop rule fired, with orbitals psi_out and step eigen tolerance tau.
   - V1, the map. One more application of the ordinary Kohn-Sham map F at rho_out: potential from rho_out, Davidson warm-started from psi_out, eigen tolerance tau_v = the tolerance the standard ratchet assigns to the next step (`StandardProfile::ratchet(tau, E, N_e)`, never looser than tau), same band-tolerance policy (occupied and first empty bands mandatory) as the ordinary steps.
   - V2, the pass criterion. The fixed-point residual of the delivered density, measured in the same metric and at the same threshold as the stop rule: E_H[F(rho_out) - rho_out] < threshold_ha. When the residual guard is on (default), also sum_k w_k sum_occupied f_i |r_i|^2 of the verification solve < threshold_ha. The sampled gap check and, with a density floor, the clipped-charge check of the strict certificate are applied to the verification map as well.
   - V3, what it does not require. Eigen-residuals at 1e-9 Ha and the 1e-7 Ha occupied-eigenvalue change (strict certificate criteria) are not required. Both are computed from the same re-solve and reported as information (`strict_criteria_met`), so the strict criteria can be compared on every verified run at no extra cost.
   - V4, what is delivered. The delivered state (density, orbitals, energy, eigenvalues) is the converged step's own state, not the verification map's output. Verification therefore never changes a delivered number; a verified run and an unverified run of the same input (`verify = "off"`) deliver bit-identical results when both stop at the same step.
   - V5, failure. A failed verification does not end the run: the session continues with the ordinary mixing step from the converged step (the verification output is discarded, the mixer history is kept), and the next time the stop rule fires it verifies again. Every attempt is recorded. The iteration budget bounds the loop; budget exhaustion is reported as not converged, exactly as today.
   - V6, cost and reporting. The verification is one warm-started Davidson solve, one density build and one Coulomb metric evaluation. Its wall time is reported as its own line (`timings.verification_s`, sum over attempts) and inside the existing certificate closure (`timings.certificate`), and its Davidson iterations and H vectors are reported under `convergence.work` (certificate counters) and in `verification.attempts`. `results.json` carries `verification` (mode, passed, attempts with the estimator, threshold, tau_v, residual bound, density L2 change, occupied eigenvalue change, work, wall) and `verified`.
5. **Downstream features.**
   - Forces (`output.forces`): require a passed verification on the single path (or the strict certificate, unchanged). A standard state without a passed verification is rejected by the engine with the existing error class (`observables.forces: strict convergence profile required ...`). On the single path with `verify = "off"`, `output.forces = true` is an input error. If the verification never passes within the budget the run ends not converged and forces error as today.
   - Stress (`output.stress`): its internal cutoff scans re-run SCF with the caller's options and fall back to "the caller's strict settings"; that fallback has no meaning when the caller is standard. Until the stress driver is adapted, a run that requests stress stays on the strict path (recorded as `single_path.status = "not_applicable"`, reason `stress`). Target behaviour: same as forces.
   - Relaxation: every ionic step needs forces from a verified state. Until the relaxation driver is adapted and the force difference standard against strict is measured (dev08 list item F), a relaxation input stays on the strict path (`not_applicable`, reason `relaxation`). Target behaviour: each ionic step verifies; a failed verification continues that step's SCF; the user chooses the threshold for relaxation and MD from the measured force differences (work order P0.1, last bullet).
   - Paths the standard profile does not support today (collinear spin, Fermi-Dirac occupations, `scf.locking = "auto"`) and the observe/agent runtime modes (checkpoint v1 does not carry the verification state) stay on the strict path, recorded as `not_applicable` with the reason.
   - Checkpoint v1 cannot carry the verification option or its attempts; encoding a continuation with verification on is rejected with a message, and exact restore of a standard history that holds a verification record is rejected (`Options`), as for other non-portable options.
6. **Switch.** `BAANDI_SINGLE_PATH=on|off` (default `off`) selects the single path for plain inputs. Every `run` with the switch on records `single_path` in `results.json` (`status`: `applied`, `not_applicable` with reason; `verification`: `fixed_point` or `off`). The switch is the rollback after the default flips.

## Why this definition

- The question a user asks of a converged standard run is whether the density is a fixed point of the Kohn-Sham map to the accuracy they asked for. The stop rule measures E_H[rho_out - rho_in], the mismatch of the *input*; V2 measures the mismatch of the *delivered* density under one more map, in the same metric and at the same threshold. For a contracting SCF it is about the square of the contraction smaller than the stop estimate, so it passes in the normal case and fails where the stop fired on a lucky step, a stall or an under-converged eigensolve (the residual guard catches the last one).
- The strict certificate's 1e-9 Ha residuals and 1e-7 Ha eigenvalue bound belong to the strict accuracy class (1e-10 Ha). Requiring them on the standard path would make the verification stricter than the threshold the user chose and would cost a tight solve on every run; that is what the strict profile is for.
- Delivering the converged step's state (V4) keeps verification an observation: it can only delay convergence, never change a delivered number, which keeps verified and unverified timings comparable and makes "switch off is bit-identical" testable.

## Consequences

- With the switch on, plain runs change path (standard instead of strict): SCF counts, energies (within the standard class, 5e-9 Ha), mixing beta (0.5) and density hashes change. With the switch off nothing changes (tested bit for bit).
- The verification cost is inside every timed BAANDI run and reported separately (work order decision 5).
- Open, measured on dev08 before the default flips (request `p01-single-path.md`): convergence and delivered error at thresholds 1e-8 and 1e-10 on si16/si32/si64 f1e-3 and the generality set; verification pass rate and cost; force differences standard against strict.

## Amendment: default flipped, profile key removed (user decision 2026-11-16)

Decision (user, 2026-11-16, in response to the code and performance review of
2026-11-16, items F01 to F04): remove `scf.convergence_profile`; the single
standard path with verification is the only path. This amendment replaces items
1, 2 and 6 above; items 3 to 5 stand with the changes noted.

1. **Input.** `scf.convergence_profile` (any value, including `"strict"` and
   `"standard"`) is an input error whose message names the replacement:
   accuracy by `scf.standard_threshold_ha` (default 5e-9 Ha, unchanged),
   verification by `scf.verify` (`fixed_point` default, `off`), strict only as a
   reference tool. Every classic `run` uses the standard path with the
   verification (status `applied`). `scf.verify = "off"` reproduces the former
   explicit `convergence_profile = "standard"` bit for bit (same options).
2. **Reference strict.** `BAANDI_REFERENCE_STRICT=on|off` (default `off`)
   replaces `BAANDI_SINGLE_PATH` in the switch table (count unchanged, 66).
   With `on` the same input runs the strict path with its certificate (status
   `reference_strict`); standard-only keys are ignored there and the result says
   so. Accuracy gates, QE validations (`scripts/verify_qe.py` rejects non-strict
   results) and strict baselines use it; `scripts/equal_accuracy.py` sets it per
   row, `scripts/bench_scf.py` refuses to compare a run against a baseline of
   the other profile and names the switch. The engine API is unchanged:
   `ScfOptions::default()` is still strict (`standard = None`), so engine tests
   and the strict reference keep their arithmetic.
3. **Strict fallbacks (item 5, unchanged in substance).** Inputs with
   relaxation, stress, collinear spin, Fermi-Dirac occupations or
   `scf.locking = "auto"` keep the strict path and record
   `single_path.status = "not_applicable"` with the reason; in these inputs the
   standard-only keys (`standard_threshold_ha`, `verify`, `standard_stall_guard`,
   `standard_density_forcing`) are input errors instead of being silently
   dropped. The observe/agent runtime modes and the non-`run` commands (bands,
   checkpoint, restart, inspection) keep a plain input on the strict path
   (checkpoint v1 cannot carry the verification); an input that writes a
   standard-only key asks for the standard path explicitly and gets it there
   (with `verify = "off"`, the former explicit standard profile that checkpoints
   carry). Forces use the passed verification (item 5); `output.forces` with
   `verify = "off"` is an input error.
4. **Recording (review F01).** Every result carries `single_path` (status,
   reason, effective profile, verification mode) and `execution`: requested
   controls, effective profile, threshold, verification, effective solver
   (`real_gamma_lean_davidson`, `real_gamma_classic_davidson`,
   `complex_davidson`, or the MPI distributed solver), real-Gamma use and its
   fallback reason, Davidson mode and lean solve count, backend (MKL or pure,
   MPI feature, ranks, threads), parallel model and every switch whose value is
   not the default. A stored baseline therefore identifies its combination.
5. **MPI replication (review F02).** With more than one k point, spin or
   Fermi-Dirac there are no k-point pools: every rank solves every k point.
   Results report `parallel.model = "mpi_replicated"`, the reason and
   `distributed_solves = 0` (the model follows the solves that ran, not the
   rank count or the M2 switch), and rank 0 warns. Replicated execution is
   deterministic: every rank follows the same branches and enters the same
   collectives; the engine test `multiple_k_points_run_replicated_on_every_rank`
   (two k points, two simulated ranks) checks bit-equal ranks and agreement with
   the threads-only run within 1e-10 Ha. Such runs are excluded from MPI
   against MPI comparisons. k-point pools are not built in this change.

Consequences for default runs: a plain Si input now runs the standard path with
the verification instead of the strict path. Its SCF count, energy (within the
standard class, 5e-9 Ha), mixing beta (0.5 instead of 0.3 for single-k strict
runs), density hash and work counts change **only because the default path
changed**; no kernel, tolerance or threshold changed. The strict path, reached
by `BAANDI_REFERENCE_STRICT=on`, is bit-identical to the former default. The
dev08 check list `docs/performance/requests/p01-single-path.md` was not run
before this user decision; its measurements (delivered error, verification
pass rate and cost, force differences) remain owed and are now measurements of
the default rather than gates of the flip.

### Residual guard is an estimator (review F03)

The residual guard (`BAANDI_STANDARD_RESIDUAL_GUARD`, default on, used by the
standard stop rule and by V2) computes B = 1 Ha^-1 * sum_k w_k sum_i f_i
||r_ik||^2. It is an **empirical estimator**, not an error bound: the
Kato-Temple form ||r||^2 / gap needs a spectral separation, and choosing the
unit gap 1 Ha does not guarantee that the actual separation is at least 1 Ha;
an SCF fixed-point residual is not the total-energy error either.
Counter-example (2 x 2, arithmetic, Python): H = diag(0, 0.01 Ha),
x = (sqrt(1 - 1e-6), 1e-3): the Rayleigh-quotient error is 1e-8 Ha while
||r||^2 * 1 Ha^-1 = 9.99999e-11 Ha, below the 5e-9 Ha threshold. This refutes
the guard as a general guarantee; it is not a run of the BAANDI SCF failing on
that input. The code now names it so (`StandardProfile::residual_energy_estimate`,
comments, switch meaning); the result key `residual_bound_ha` keeps its name for
compatibility. No number, threshold or tolerance changed. The adequacy of the
default rests on measured delivered error against strict references (dev08
list above); a true bound would need an occupied/unoccupied spectral
separation estimate that handles degeneracy and finite occupations, which is
not designed here.

## Amendment: forces and fixed-cell relaxation on the single path (ONEPATH-1, user decision 2026-11-16)

Decision (user, 2026-11-16): the strict path disappears from production use in
stages: (1) forces and relaxation (this amendment, ticket
`tickets/ONEPATH-1-FORCES-RELAX.md`), then stress, then spin and metals, then
the engine default; strict stays only behind `BAANDI_REFERENCE_STRICT=on` as an
accuracy reference tool. This amendment replaces the relaxation bullet of item 5
and the word "Relaxation" in item 3 of the first amendment.

1. **Forces.** Unchanged in substance (item 5): a passed verification makes the
   standard state eligible, and the engine computes the forces with the same
   formulas from the converged step's own state (V4). `forces.preconditions`
   now also records `standard_threshold_ha`, `verification_attempts` and the
   last attempt's `verification_coulomb_error_ha` and
   `verification_residual_bound_ha` (additive optional keys in
   `docs/schemas/forces-available.proposal.json`).
2. **Fixed-cell relaxation** (`relax.mode = "ions"`, the default) runs every
   ionic step on the standard path with the verification. A failed
   verification continues that step's SCF (V5); density transfer, warm starts
   and `start_extrapolation` are unchanged; the half-step retry and every
   optimizer rule are unchanged. The `relaxation => strict not_applicable`
   fallback is removed: `single_path.status = "applied"`. `scf.verify = "off"`
   with `[relax]` is an input error (every step needs forces).
3. **Threshold.** The relaxation threshold is `scf.standard_threshold_ha`
   (default 5e-9 Ha, unchanged). Users tighten it for force-sensitive runs
   (tight `force_tol_ha_bohr`, small `energy_tol_ha`): the SCF energy is
   resolved to about the threshold, so an ionic energy-change test far below
   it cannot pass reliably. Measured force differences against the strict
   reference are in ledger row ONEPATH-1 (Mac, synthetic fixtures) and are
   requested on dev08 (`docs/performance/requests/onepath1-forces-relax.md`).
4. **Still strict.** `output.stress` (reason `stress`) and the variable-cell
   relaxation modes `cell` and `ions_cell` (reason `stress (variable-cell
   relaxation)`) keep the strict path, and standard-only keys in those inputs
   stay input errors. Spin, Fermi-Dirac, `locking = "auto"`, the observe/agent
   modes and the non-`run` commands are unchanged.
5. **Recording.** The relaxation `results.json` carries `single_path`, an `scf`
   block (effective profile, threshold, verification mode) and per frame a
   `verification` record (passed, threshold, attempts, last Coulomb residual,
   residual estimate, `strict_criteria_met`, H vectors). It holds no wall time,
   because frames are stored in `restart.json` and a resumed relaxation must
   reproduce the uninterrupted files bit for bit.
6. **Reference identity.** With `BAANDI_REFERENCE_STRICT=on` a relaxation runs
   the strict path with the same options as before this amendment (frames carry
   no verification record); only the added top-level `single_path` and `scf`
   keys of `results.json` differ from the former output.

## Amendment: stress and variable-cell relaxation on the single path (ONEPATH-2, user decision 2026-11-16)

Decision (user, 2026-11-16): stage 2 of removing the strict path from
production use (ticket `tickets/ONEPATH-2-STRESS.md`). This amendment replaces
item 4 of the ONEPATH-1 amendment and the stress bullet of item 5.

1. **Stress** (`output.stress`) runs the standard path with the verification
   (`single_path.status = "applied"`). The engine accepts a standard state for
   stress exactly when it accepts it for forces (a passed verification; the
   stress-specific rejection of any standard state is removed). The tensor is
   computed with the same formulas from the converged step's own state (V4).
   `stress.preconditions` carries the forces' standard keys (threshold,
   verification passed, attempts, last Coulomb residual and residual estimate);
   the results proposal schema gives `stress.preconditions` a strict branch
   (unchanged) and a verified-standard branch. `verify = "off"` with
   `output.stress` is an input error.
2. **Cutoff scans.** The isotropic Pulay estimate re-runs three SCFs (Ecut,
   1.05 Ecut, 1.10 Ecut). For a standard caller every scan runs the standard
   path with the verification at `min(1e-9 Ha, scf.standard_threshold_ha)`:
   never looser than the former scan threshold and following a caller who
   tightens. A scan that fails there is repeated once at the caller's own
   threshold when that is looser; there is no strict fallback. The SCF term of
   the uncertainty uses max(`energy_tolerance`, loosest scan threshold used).
   For a strict caller (`BAANDI_REFERENCE_STRICT=on`) the scans are unchanged
   bit for bit (standard at 1e-9 Ha unverified, strict fallback).
   `pulay_assessment` records `scan_profile`, `scan_thresholds_ha` and
   `scan_energy_resolution_ha`.
3. **Variable-cell relaxation** (`relax.mode = "cell"`, `"ions_cell"`) runs every
   step, including basis rebuilds and the final nominal-cutoff point, on the
   standard path with the verification; frames carry the ONEPATH-1
   `verification` record; the frozen Miller chart, rebuilds and Pulay handling
   are unchanged. The reasons `stress` and `stress (variable-cell relaxation)`
   are removed. Still strict: collinear spin, Fermi-Dirac, `locking = "auto"`
   (also when they request stress), the observe/agent modes and the non-`run`
   commands; the Python session's `stress()` uses the engine API (strict
   default) and is unchanged.
4. **Stress sensitivity to SCF convergence (measured, Mac, synthetic).** Stress,
   like forces, is first order in the density error while the energy is second
   order, so the stress error falls roughly as the square root of the delivered
   energy error. Pre-registered bound (ticket criterion a, the synthetic stress
   oracle cap 1e-7 Ha/Bohr^3 of `docs/validation/stress.md`): at 5e-11 Ha both
   fixtures pass (cubic 6.5e-10, skewed nonlocal 6.1e-8 Ha/Bohr^3); at the
   default 5e-9 Ha the cubic fixture passes (2.1e-9) and **the skewed nonlocal
   fixture fails (1.305e-7 Ha/Bohr^3, relative 2.7e-5)**. Against a tightened
   strict reference the default strict path is within 8.6e-11, so the
   exceedance is the standard path's own error, not the reference's
   (decomposition: 5e-10 -> 5.9e-8, 5e-13 -> 5.0e-9 Ha/Bohr^3 on the skewed
   fixture). The bound is not widened; the failing test is kept and marked.
5. **Stress default threshold (orchestrator decision 2026-11-16, PROVISIONAL).**
   Under the adoption rule "the default must meet the accuracy class", an input
   that asks for stress (`output.stress`, or `relax.mode = "cell"` or
   `"ions_cell"`) and does not set `scf.standard_threshold_ha` runs at
   `STRESS_DEFAULT_THRESHOLD_HA` = 5e-11 Ha. This is a rule by observable, not by
   material. Forces, `relax.mode = "ions"` and plain runs keep 5e-9 Ha; an
   explicit user value always wins, tighter or looser. `results.json` records
   the effective threshold and its origin in `single_path.standard_threshold_ha`
   and `single_path.threshold_source` (`user`, `stress_default` or `default`;
   null off the standard path); the relaxation `scf` block carries the same two
   keys. Evidence: two synthetic fixtures on the Mac (item 4), with a margin of
   only 1.6x on the skewed nonlocal fixture (6.1e-8 against 1e-7 Ha/Bohr^3). The
   rule is provisional until the dev08 list
   `docs/performance/requests/onepath2-stress.md` (si8, si16, si-skew, si-slab,
   strained si8 at 5e-9, 5e-10 and 5e-11 Ha against a tightened strict
   reference) confirms it; if dev08 shows 5e-11 missing the class anywhere, the
   default is revisited, not the bound. Cost: on the skewed fixture 44 SCF
   iterations at 5e-11 against 21 at 5e-9 (strict 17). Variable-cell
   relaxations at 5e-9 reached the strict reference's cell within 1.6e-8 Bohr
   (BCC, isotropic) and 5.4e-6 Bohr (compact skewed, full cell and ions), well
   inside the pre-registered 1e-3 Bohr; the numbers at the new default are in
   ledger row ONEPATH-2-STRESS-DEFAULT.
6. **Reference identity.** With `BAANDI_REFERENCE_STRICT=on` stress runs and
   variable-cell relaxations run the strict path with the same options and scan
   routing as before; the stress block gains only the additive
   `pulay_assessment.scan_*` keys.

## Amendment: collinear spin and Fermi-Dirac occupations on the single path (ONEPATH-3, user decision 2026-11-16)

Decision (user, 2026-11-16): stage 3 of removing the strict path from
production use (ticket `tickets/ONEPATH-3-SPIN-METALS.md`, design
`docs/design/onepath3-spin-metals.md`). This amendment replaces the spin and
Fermi-Dirac parts of item 5 ("Paths the standard profile does not support"),
of item 3 of the first amendment and of item 3 of the ONEPATH-2 amendment.

1. **Routing.** Collinear spin (fixed sector and Fermi-Dirac) and scalar
   Fermi-Dirac inputs run the standard path with the verification
   (`single_path.status = "applied"`); the reasons `collinear spin` and
   `Fermi-Dirac occupations` are removed; standard-only keys are accepted.
   `scf.locking = "auto"` keeps the strict fallback. The engine no longer
   rejects a standard profile for spin channels; it rejects the Hartree mixing
   metric (unchanged) and the experimental stall guard / density forcing
   (closed-shell experiments) for them.
2. **Estimator.** Collinear: E_spin = E_H[d_n] + E_H[d_m] + E_0[d_m], the
   charge and magnetization residuals in the same Coulomb metric plus the
   uniform magnetization (total-moment change, possible only under an
   unconstrained FD moment) weighted like the longest wavelength of the
   density sphere. It reduces to the closed-shell estimator when both channels
   are equal. Like the residual guard it is an estimator, not an energy bound
   (the magnetization's energy cost is set by the spin stiffness, not by
   4 pi / G^2). Residual guard: sum over spin channels and k of
   w f |r|^2 with the actual occupations (capacity 1 per channel, 2 closed
   shell; fractional under FD). Ratchet, first-iteration redo, Ritz-change
   Davidson and band policy (all retained FD bands mandatory) as on the
   closed-shell standard path; complex solver (real Gamma excludes spin and FD).
3. **Verification.** V1 repeats the whole map at the delivered densities,
   including the Fermi-level solve; V2 passes on the same estimator and the
   guard of the verification solve below the threshold plus the number, moment
   and band-edge checks of the strict certificate, and (fixed sector) the
   channel gaps; V3 records charge and spin density changes, occupation change,
   free-energy change, chemical-potential change and moment change as
   information (`strict_criteria_met`); V4 delivers the stop step's own state
   (energies, F, mu, moment, occupations); V5 and V6 unchanged.
4. **Downstream.** Forces and stress accept a spin or FD standard state with a
   passed verification (same rule as ONEPATH-1/2). Stress cutoff scans of a
   standard caller run standard with verification; for a strict caller with a
   spin or FD state the scans go straight to strict (before this amendment the
   standard first attempt always failed for these problems), so the reference
   stays bit-identical. Relaxation with spin or FD and the observe/agent modes
   with spin or FD stay input errors.
5. **Measured (Mac, synthetic, ledger ONEPATH-3; bounds fixed before
   measuring).** (a) zero-moment collinear against closed shell on the standard
   path: |dE| 1e-15 Ha, same SCF count (bound 1e-10 Ha, +-1 iteration). (b)
   standard 5e-9 Ha against the strict reference: ferromagnetic fixed sector
   |dE| 2.3e-10 Ha, |dM| 2e-16 mu_B; scalar FD metal (2x2x2 mesh) |dF| 2.0e-12
   Ha, |dmu| 8.9e-7 Ha, sum w|df| 3.8e-6 e; collinear FD with free moment |dF|
   5.5e-11 Ha, |dM| 6.7e-6 mu_B, |dmu| 6.6e-6 Ha, sum w|df| 6.7e-6 e (bounds
   5e-9 Ha, 1e-4 mu_B, 1e-4 Ha, 1e-4 e). SCF counts 15/5/12 against strict
   45/9/39. dev08 request: `docs/performance/requests/onepath3-spin-metals.md`.
6. **Reference identity.** With `BAANDI_REFERENCE_STRICT=on` spin and FD runs
   take the strict drivers with the same options as before; results gain only
   the `single_path`/`execution` keys they already had.
