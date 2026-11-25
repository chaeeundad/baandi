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
