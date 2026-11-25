# ADR 0045: Oldest-first conditioning of Broyden histories

Status: accepted (R6A-G4B-DEFAULT, 2026-11-25). Enabled by default for collinear FD.

## Motivation and scope

Free-moment collinear FD reaches nearly dependent, widely scaled residual
secants while approaching the strict oracle. The current unnormalized Gram
solve rejects at conditions near 1e15–1e16, even where the true energy already
meets the tight accuracy class. This is distinct from premature standard
stopping; history recovery alone cannot certify delivered accuracy.

`REFSCF_HISTORY_CONDITIONING=on` selects conditioning for the collinear FD reference entry
(default on). Scalar FD, fixed collinear sectors, fixed unpolarized drivers
and legacy drivers retain the old arithmetic even with the switch on. The
Broyden adapter exposes explicit setup-time opt-in for manufactured tests.
Default adoption is supported by the synthetic and dev08 evidence below.

## Accepted rule

For p retained unnormalized pairs, diagonalize the same metric Gram matrix
A with the existing symmetric EVD adapter. Define rcond=lambda_min/lambda_max.
Require rcond > sqrt(p*epsilon_f64), positive finite spectrum maximum and a
finite condition. The factor p*epsilon is the scale of accumulated floating
point perturbations; its square root conservatively limits inverse-Gram
amplification to approximately half the available digits. This is a numerical
conditioning policy, not a material-specific convergence tolerance. It is
stricter than the legacy numerical-singularity criterion p*epsilon.

When the condition fails, remove the oldest retained (dx,dF) pair, rebuild
the Gram matrix and RHS, and repeat. Record each removal's pre-drop condition
in chronological order. A zero Gram spectrum also removes its oldest pair.
If no pair remains, take the ordinary beta-preconditioned linear step, saving
the current x and F as the next secant's origin. Nonfinite data and EVD/backend
errors still fail; this rule does not mask them. Ordinary full-history FIFO
eviction is not a conditioning drop. Redos and convergence/verification maps
never insert pairs. Failed verification retains the ordinary input pair.

Do not normalize secants, add a ridge, use a pseudoinverse, adjust beta or
project an unconstrained uniform spin mode. The unchanged linear solve acts
on the retained history. Scratch storage is allocated at construction. Gram
reuse is disabled in the conditioned FIFO branch to avoid stale ring indices;
the existing ring/cache remains untouched with the switch off. The optional
ring-based secant predictor returns unavailable for a conditioned history;
the electronic FD entry does not use it.

## Compatibility and acceptance

Off keeps the existing min <= p*epsilon*max rejection, insertion order and
arithmetic. Every drop is explicit in `Diagnostics` and optional per-map traces.
Conditioning changes the singular-history contract only when selected.
Acceptance: all 48 G4 synthetic comparisons pass with the three switches on
at unchanged energy, moment, chemical-potential and weighted-occupation bounds;
Si2 FD and already qualified paths retain their recorded identity. Detailed
measurements live in [Spin and Fermi-Dirac](../design/reference-scf.md#spin-and-fermi-dirac).

The R6a work order records dev08 gen15 bcc Fe (PBE+NLCC, 4x4x4 k mesh,
fixed moment 4.4 mu_B) VERIFIED in 75.8 s on the reference driver; legacy
and the strict oracle did not finish within 45 min. This is convergence
and verification evidence, not a completed strict-oracle accuracy comparison
or a QE timing comparison. It was supplied by the orchestrator, not rerun locally.

R6a enables REFSCF_SPIN_FD_MOMENT, REFSCF_SPIN_STIFFNESS and
REFSCF_HISTORY_CONDITIONING by default and removes the L2 routing predicate.
Explicit off values remain available for A/B diagnostics without changing
driver selection. Free-moment FD with PBE remains on legacy under L1/G4D.
Local rank simulation is not real-MPI or MKL validation.


## G4F stiffness screening correction

The local LDA-PW stiffness proxy follows ADR 0046 on the combined valence plus
half-core channels. Freeze the accepted XC mask for both finite-difference
probes: screened channels have zero response and fully screened points have
zero stiffness. Shorten a probe near the threshold to keep the mask fixed,
using its actual span (possibly one-sided); do not introduce a density floor.
The globally audited input is not subjected to a second rank-local charge gate
on derivative scratch. Above the threshold with unchanged probe spans, and
with screening off, the original arithmetic is retained. This corrects the
stiffness evaluator only; no convergence threshold or screening budget changes.
See the stiffness paragraph in `docs/design/reference-scf.md` for the step rule.

Local regression fixture: synthetic one-electron ion in a 12 Bohr cubic box,
2.6 Ha wavefunction cutoff, Gaussian NLCC and a short-range Gaussian well,
12 bands, FD theta=0.005 Ha. The fixed moment is 0.9999999 mu_B (strictly
inside the finite-chemical-potential domain). The original evaluator passes
one rank but fails at three simulated ranks with the reported generic
stiffness error; the main XC charge gate has passed. With the correction,
fixed-moment PBE and free-moment LDA-PW converge and verify at simulated ranks
1/3/8. Maximum free-energy differences from converged strict legacy runs are
1.6056489471338864e-11 and 1.1717293801893902e-11 Ha, respectively, below the
unchanged 5e-9 Ha standard class. The free variant starts fully polarized to
compare the same stationary spin branch; a partial seed reached different
stationary solutions in the two drivers and is not an accuracy oracle.

Reproduce with `cargo test --release -p dft-engine --test refscf g4f_vacuum_fd
-- --include-ignored --nocapture`. Both tests are ignored in the debug gate.
`cargo test --release -p dft-engine --lib g4f -- --nocapture` covers frozen
masks, threshold-adjacent positive probes, zero vacuum response, an independent
energy-curvature stencil, unscreened on/off bits, and collective errors from
one rank to an empty owner for start/map/verification. The 48 existing G4
synthetic fingerprints are byte-identical before and after this correction.
These are local synthetic checks, not a dev08 O rerun, real MPI validation,
a QE comparison, or a performance claim.

The unfiltered release `g4` command also runs two unrelated tests that already
fail on the baseline checkout: `g4_si_gap` requires the external `DFT_SI_UPF`,
and `g4d_free_moment_pbe_acceptance` rejects the unqualified legacy free-moment
PBE no-core case. Neither is hidden or relaxed by this correction; they are
separate from the passing G4F regressions and unchanged G4 synthetic matrix.

Final local gates: `cargo check --workspace --tests`, the debug `refscf`
suite (84 passed, 26 ignored), and targeted fmt/clippy all exit 0.
`CHECK_NICE=0 CHECK_SCCACHE=0 BAANDI_TEST_TMP=/private/tmp/g4f-gates
scripts/check.sh fast dft-engine dft-physics` exits 0: build/test/fmt/clippy
pass, 45 binaries, 462 passed, 129 ignored. Its informational `ps` probe is
denied by the sandbox. Excluding only the two baseline failures named above,
the release G4 suite passes all 12 tests. README generation and freshness
checks pass with no generated changes; reporting JSON and new README evidence
are outside this ticket's owned files. Changes are left uncommitted for review.
