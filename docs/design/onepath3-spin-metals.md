# ONEPATH-3: collinear spin and Fermi-Dirac occupations on the single standard path

Status: design for ticket `tickets/ONEPATH-3-SPIN-METALS.md` (user decision
2026-11-16, stage 3 of removing the strict path from production use). Written
before the implementation; the acceptance bounds in section 7 were fixed before
any measurement and are not changed afterwards.

## 1. Starting point

- The standard path (ADR 0035, ADR 0043) exists only in `ScfSession` (closed
  shell, integer occupations). Its stop rule is the Coulomb estimate
  E_H[rho_out - rho_in] < threshold (default 5e-9 Ha) together with the residual
  guard B = 1 Ha^-1 sum_k w_k sum_i f_i |r_ik|^2 < threshold; the per-step eigen
  tolerance follows the ratchet tau_{n+1} = clamp(min(tau_n, 0.1 E_n / N));
  after the stop the fixed-point verification (V1 to V6) runs one more map at the
  delivered density.
- Collinear spin (`scf::collinear`, fixed sector and FD) and scalar Fermi-Dirac
  (`scf::finite_temperature`) have their own drivers with the strict stop rule
  and the strict fresh-map certificate. `ScfProblem::construct` rejects a
  standard profile for any spin channel (also the scalar FD problem, which is
  built as channel 0), the CLI routes both to strict with reasons
  `collinear spin` / `Fermi-Dirac occupations`, and forces/stress of such a state
  accept only the strict certificate.

## 2. Stop estimator with two spin densities

Notation: n = n_up + n_down (charge), m = n_up - n_down (magnetization),
residuals d_s = n_s,out - n_s,in, d_n = d_up + d_down, d_m = d_up - d_down. K is
the Coulomb metric of `CoulombMetric` (r.K.r = E_H[r], G = 0 omitted, density
sphere only).

Estimator:

    E_spin = E_H[d_n] + E_H[d_m] + E_0[d_m]
           = 2 (E_H[d_up] + E_H[d_down]) + E_0[d_m]

with the uniform term E_0[d_m] = (Omega/2) (4 pi / G_min^2) <d_m>^2, where
<d_m> = Delta M / Omega is the cell average of d_m and G_min the shortest
nonzero reciprocal vector of the density sphere (the largest Coulomb weight of
the metric).

Reasons:

- Charge part: identical to the closed-shell estimator, so the standard
  threshold keeps its QE-comparable meaning (Hartree-metric error estimate).
- Magnetization part: the spin residual is invisible to the Hartree term (it
  produces no electrostatic potential), so it needs its own term. It is measured
  in the same positive metric so that one threshold governs both and the two
  parts are commensurate in Ha. This is a measuring rule, not an energy bound:
  the true second-order energy cost of a magnetization error is set by the spin
  stiffness (XC and kinetic response), which is local and not 4 pi / G^2. Like
  the residual guard (ADR 0043, review F03) it is an empirical estimator whose
  adequacy rests on measured delivered errors against the strict reference
  (section 7) and on the verification.
- Uniform mode: the charge residual has no G = 0 component (both densities hold
  N electrons). The magnetization residual has none in a fixed sector or under a
  fixed-moment FD constraint (per-channel number solves), but under an
  unconstrained FD moment the total moment of the output can differ from the
  input. A positive metric must give that mode a weight; it is weighted like the
  longest wavelength the cell resolves. The term vanishes identically in the
  other cases and in the unpolarized limit, so it is always evaluated.
- Unpolarized limit: d_up = d_down gives d_m = 0 and E_spin = E_H[d_n], the
  closed-shell estimator exactly (test a).

## 3. Residual guard with spin and occupations

    B = 1 Ha^-1 sum_s sum_k w_k sum_i f_{s,k,i} |r_{s,k,i}|^2

with f the occupation of the band in electrons including the channel capacity
(1 per spin channel, 2 for the scalar FD problem; fractional under FD), k
weights excluded from f, bands with f = 0 excluded (exact zeros only; every FD
band with a positive fraction counts with its fraction). Unpolarized limit: two
channels with f = 1 give the closed-shell 2 |r|^2. The guard stays an
estimator, not a bound (ADR 0043 F03); under FD its Kato-Temple motivation is
weaker still, because partially occupied bands near mu have no gap to the
neighbouring states. It can only delay a stop.

## 4. Eigen tolerances and band policy

- Ratchet unchanged: first step 5e-3 Ha (`StandardProfile::CAP_HA`), then
  `StandardProfile::ratchet(tau, E_spin, N)` with N the total electron count.
  The first-iteration redo of the closed-shell path (estimate below tau N at
  iteration 1: repeat the map at the ratcheted tolerance) is reproduced, so the
  unpolarized limit follows the same schedule.
- Eigensolver: the standard Ritz-change Davidson (`solve_ritz_changes`), as the
  closed-shell standard path on the complex solver. Real Gamma stays excluded
  for spin and FD (ADR 0040 exclusions unchanged); `results.json`
  `execution.effective_solver` reports `complex_davidson`.
- Band policy: fixed sector as today (occupied bands and the first empty band
  mandatory, standard empty threshold for higher buffers). FD: every retained
  band stays mandatory (the drivers mark all bands occupied for the policy),
  because every band can carry a fractional occupation.
- The scalar FD driver today solves every step at 0.1 x eigen_tolerance; on the
  standard path it follows the ratchet instead.

## 5. Stop rule

Standard stop at iteration n: E_spin,n < threshold, B_n < threshold (guard on,
the default), and the number conditions that the strict driver already
requires (channel electron counts within 1e-10 N for a fixed sector; FD number
and moment residuals within 1e-10 N; FD edge coverage, which is an input error
as today). The opt-in `BAANDI_STANDARD_STOP=confirm` variant and the
experimental stall guard / density forcing are not carried over to these
drivers (documented; the switches keep their meaning on the closed-shell path).

## 6. Verification map (V1 to V6 with spin and occupations)

- V1, the map. At the delivered input pair (the stop step's output densities
  of both channels): potentials, warm-started Davidson for every spin channel
  and k point at tau_v = ratchet(tau, E_spin, N), then for FD a fresh number
  solve (new mu, per-spin mu under a fixed moment, new occupations, new moment
  under an unconstrained moment), then the output densities. This is the whole
  Kohn-Sham map including the occupation step; nothing of the stop step's
  occupations is reused.
- V2, pass. E_spin[F(rho) - rho] < threshold, and with the guard on the B of
  the verification solve (with the verification's own occupations) <
  threshold. Validity checks of the strict certificate that are not
  convergence tolerances also apply: channel electron counts (fixed sector),
  FD number/moment residuals (1e-10 N) and edge coverage; per-channel sampled
  gaps of a fixed sector on the stop step and on the verification solve
  (required max(1e-7 Ha, 2 x eigen_tolerance), same rule as strict: a closed
  gap is a NonInsulating error, not a failed verification).
- V3, information only. Charge and spin density L2 changes, occupied
  eigenvalue change, and for FD the weighted occupation change
  sum w |Delta f|, the free-energy change, the change of each chemical
  potential and of the total moment. `strict_criteria_met` reports whether the
  strict certificate criteria (density <= density_tolerance for charge and spin,
  energy or free-energy change <= energy_tolerance, residuals <= eigen_tolerance,
  occupied eigenvalue change <= 1e-7 Ha for a fixed sector, occupation change
  <= 1e-8 N for FD) would have held.
- V4, delivered. The stop step's own state: densities, orbitals, eigenvalues,
  energies, occupations, mu, entropy term and moment. The verification's
  occupations and mu are discarded (the existing restore of the FD workspace),
  so verification never changes a delivered number.
- V5, failure. A failed attempt continues with the ordinary mixing step from
  the stop step (mixer histories of both channels kept); the next stop verifies
  again; the iteration budget bounds the loop and exhaustion is "not converged".
- V6, recording. `CollinearResult` and `FiniteTResult` carry `standard`,
  `verification` (the shared `VerificationReport`) and per-iteration
  `standard_error_ha`; `VerificationAttempt` gains optional
  `spin_density_residual`, `occupation_change_electrons`,
  `chemical_potential_change_ha` and `moment_change` (None on the closed-shell
  path). The CLI writes `verification`, `verified`, `timings.verification_s`,
  `single_path` (status `applied`) and `execution` into the spin and FD
  `results.json`.

## 7. Acceptance (bounds fixed before measuring)

- (a) Unpolarized limit, engine level: collinear fixed sector n_up = n_down on
  the standard path against the closed-shell standard path with the same
  options (complex solver on both): |Delta E| <= 1e-10 Ha, density distance
  <= 1e-9 and SCF counts within 1 (the bounds of the existing strict test
  `zero_moment_reduces_to_same_unpolarized_functional`). Fallback reading if the
  two trajectories cannot be made identical: not applied; a failure is reported.
- (b) Fixtures, standard (default 5e-9 Ha, verification on) against
  `BAANDI_REFERENCE_STRICT=on` through the CLI:
  - ferromagnetic synthetic fixture (fixed sector N_up > N_down, nonlocal
    projector, skewed cell);
  - smeared metallic synthetic fixture (scalar FD, theta = 0.01 Ha, k mesh);
  - spin-polarized smeared fixture (collinear FD with unconstrained moment).
  Bounds: |Delta E| (free energy F under FD) <= 5e-9 Ha (standard class);
  |Delta M| <= 1e-4 mu_B; |Delta mu| <= 1e-4 Ha per chemical potential; weighted
  occupation difference sum w |Delta f| <= 1e-4 electrons. The last three are
  first-order quantities (error about the square root of the energy error), so
  they get a first-order bound; they are stated with the measured values.
- (c) Verification passes on every converged standard run of (b); a forced
  failure (test hook) continues the SCF and the run still converges with two or
  more attempts.
- (d) `single_path.status = "applied"` with no reason for spin and FD runs;
  `verify = "off"` gives an unverified standard run; forces with
  `verify = "off"` stay an input error; relaxation of spin/FD stays an input
  error (unchanged); standard-only keys are accepted.
- (e) All existing spin and smearing tests pass; tests that test the strict
  reference set `BAANDI_REFERENCE_STRICT=on` (CLI) or keep the engine API
  default (`ScfOptions::default()` is strict).

## 8. Downstream features

- Forces: a spin or FD standard state with a passed verification is accepted
  exactly like a closed-shell one (ONEPATH-1 rule); unverified standard states
  are rejected.
- Stress: accepted on the same rule. Cutoff scans of a standard caller run the
  standard path with verification (ONEPATH-2 rule, now possible for spin/FD
  problems). For a strict caller with a spin or FD state the scans go straight to
  the strict settings: before this change the standard first attempt was always
  rejected by the engine for these problems, so going straight to strict keeps
  the reference bit for bit.
- Relaxation with spin or FD stays an input error ("relax currently requires
  unpolarized fixed occupations"), unchanged.
- Runtime modes and non-run commands keep their behaviour (spin/FD already
  reject runtime control; bands with FD go through `options()` and stay strict
  for plain inputs).

## 9. Still unsupported (input errors, documented)

`scf.locking = "auto"` (stays strict `not_applicable`, as before), the Hartree
mixing metric with spin or FD (engine error as today), spin with `lda-pz`,
spin/FD with observe/agent modes and with relaxation, FD with a nondefault
preconditioner. No new error is introduced.

## 10. Risks

- The magnetization term is an estimator without an energy-bound proof
  (section 2). Mitigation: verification at the same metric and the measured
  comparison (b); the dev08 request extends the measurement.
- FD near a gap: mu is ill-conditioned (ADR 0032); its difference is a
  diagnostic there. The metallic fixture is not near a gap.
- Charge sloshing in metals with an unconstrained moment can make the
  verification fail repeatedly; V5 then continues and the budget decides.
- Results of spin/FD runs change from the strict path (SCF counts, energies
  within the class, mixing beta 0.5 for single-k runs as on the closed-shell
  standard path); the strict reference stays bit-identical behind the switch.
