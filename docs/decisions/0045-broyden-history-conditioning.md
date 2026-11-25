# ADR 0045: Oldest-first conditioning of Broyden histories

Status: proposed (G4B-SPIN-MOMENT, 2026-11-25). Defaults remain off.

## Motivation and scope

Free-moment collinear FD reaches nearly dependent, widely scaled residual
secants while approaching the strict oracle. The current unnormalized Gram
solve rejects at conditions near 1e15–1e16, even where the true energy already
meets the tight accuracy class. This is distinct from premature standard
stopping; history recovery alone cannot certify delivered accuracy.

`REFSCF_HISTORY_CONDITIONING=on` opts the collinear FD reference entry into
the candidate. Scalar FD, fixed collinear sectors, fixed unpolarized drivers
and legacy drivers retain the old arithmetic even with the switch on. The
Broyden adapter exposes explicit setup-time opt-in for manufactured tests.
No performance/default promotion is established by local tests.

## Proposed rule

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
reuse is disabled in the candidate FIFO branch to avoid stale ring indices;
the existing ring/cache remains untouched with the switch off. The optional
ring-based secant predictor returns unavailable for a conditioned history;
the electronic FD entry does not use it.

## Compatibility and acceptance

Off keeps the existing min <= p*epsilon*max rejection, insertion order and
arithmetic. Every drop is explicit in `Diagnostics` and optional per-map traces.
The candidate changes the singular-history contract only when selected.
Qualification requires all 48 G4 comparisons at unchanged energy, moment,
chemical-potential and weighted-occupation bounds, Si2 FD, and bit identity
of already qualified paths. A negative result leaves moment policies
unqualified and defaults off. Local rank simulation is not dev08 real MPI
or MKL validation. Detailed measurements live in the Spin and Fermi-Dirac
section of `docs/design/reference-scf.md`.
