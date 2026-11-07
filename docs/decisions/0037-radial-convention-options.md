# ADR 0037: Research radial convention options

Date: 2026-11-07. Status: **Proposed**. Defaults remain unchanged.
Task: P27D-RADIAL-CONVENTION-SCORE.

## Intent and implemented contract

Expose two independent research choices under `[pseudo]`:
`local_mesh = "baandi" | "qe_like"` and
`radial_even_rule = "simpson_3_8" | "simpson_quadratic_end"`.
Omission selects the first variant of each. Unknown values are rejected.
No convention is claimed to be more correct; a default change would require
its own decision and evidence beyond this diagnostic.

The alternate local count is the 1-based index of the first sample strictly
beyond 10 Bohr, or the full count if absent, rounded down to odd. Both the
local transform and ionic-local reconstruction use that count. Core density
keeps its current full mesh; projectors keep their stored cutoff prefixes.

The alternate even-count quadrature uses Simpson 1/3 on the first n-1 points
and weights (-1, 8, 5)/12 on the last three samples for the final interval.
It reaches local, core and projector transforms through additive APIs.
Default entry points and numerical bodies are unchanged. Odd counts use
the existing rule. The two-point fallback remains trapezoidal.

These definitions come only from the supplied own-words note
`../research/qe-study/qe-radial-endpoint-convention.md`. No QE source was read.

## Architect corrections

The original gate requiring cubic exactness from both even-count rules
conflicted with the prescribed weights. The architect approved the corrected
gate: Simpson 3/8 is cubic-exact, the quadratic end interval is quadratic-exact,
and its preceding Simpson panels remain cubic-exact. Tests retain the specified
weights and explicitly verify the +h^4/4 monic-cubic error. The four-point
unit-step result is 20.5 versus the exact 20.25. Odd-count identity is tested.

The architect also extended ownership to the CLI checkpoint wire implementation
and additive option plumbing. Missing fields in old checkpoints default to the
legacy choices; nondefault strings are serialized and both knobs participate
in exact-restart comparison. Unknown wire values reject. No checkpoint version
or existing public function signature changed.

## Metadata and capability limits

SCF options and scalar results retain both enum choices. CLI output records
nondefault choices under `scf.pseudo` and, for the scalar result, `pseudo`.
Default omission preserves old output shape. The shared construction path
carries the knobs through other SCF modes and NSCF; the evidence in this ticket
is unpolarized PBE only.

Derivative model identity includes both options. The separate analytic force
and stress radial derivative routines are not extended by this diagnostic;
nondefault conventions explicitly reject derivative requests. This prevents
mixing SCF energies with derivatives of the legacy Hamiltonian. Default
force/stress behavior is unchanged.

## Measured consequence

The full 2 x 2 x 6 SCF table is in
`../validation/radial-convention-score.md`. One common combination, qe_like
plus simpson_quadratic_end, reduces all six residuals to the 1e-8 Ry scale
or below. The other combinations do not explain Si and MgO together at that
scale. This is a candidate explanation under this input protocol, not universal
cross-code equivalence. QE's relayed total energies are printed to eight
Ry decimal places (approximately +/-5e-9 Ry rounding uncertainty).
No per-species choices or tolerance tuning were used. The default remains
baandi plus simpson_3_8 regardless of the diagnostic outcome.
