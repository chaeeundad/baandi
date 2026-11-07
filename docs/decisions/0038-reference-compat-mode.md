# ADR 0038: named reference compatibility mode

Status: Accepted by the user on 2026-11-08.
Ticket: P44A-REFERENCE-MODE.

## Decision

`[compat] mode = "reference"` selects four existing conventions for unpolarized
PBE energies and densities. This is configuration resolution, not new numerics:

| Setting | Resolved value | Existing contract |
| --- | --- | --- |
| `xc.pw92_prefactor` | `original` | [ADR 0036](0036-pw92-prefactor-convention.md) |
| `pseudo.local_mesh` | `qe_like` (reference-style local mesh) | [ADR 0037](0037-radial-convention-options.md) |
| `pseudo.radial_even_rule` | `simpson_quadratic_end` | [ADR 0037](0037-radial-convention-options.md) |
| `symmetry.density_group` | `crystal` when `symmetry.mode="auto"`; otherwise `mesh` | [Full-group contract and evidence](../validation/symmetry-full-group.md) |

The symmetry convention predates this ADR but has no separate numbered ADR;
this decision records its inclusion, without changing k-point reduction or
symmetry activation. The switch does not turn symmetry on.

An explicit selected key with the same value is allowed. A differing value
rejects with both `compat.mode = "reference"` and the conflicting key named.
With symmetry off or absent, an explicit `density_group="crystal"` conflicts;
the resolved density group remains the default `mesh`.

Spin-polarized PBE, LDA functionals, analytic forces, stress and relaxation are
unsupported. The mode identifies its complete scope in one diagnostic:
`compat.mode = "reference" is for unpolarized PBE energies and densities`.
Existing restrictions on UPFs, occupations, symmetry, checkpoint capture and
other capabilities remain applicable. This decision does not extend them.

Keys, messages and result metadata use the neutral name `reference`, never a
reference implementation's product name. The established `qe_like` enum value
is retained unchanged because existing inputs and evidence depend on it.
External code names occur only in evidence documents citing measurements.

## Evidence and limits

The [evidence and user guide](../validation/reference-mode.md) records the six
Si/MgO totals and the exact reproduction commands. With the original PW92
prefactor and the two radial options, the previously measured residuals range
from -9.3061e-9 to +3.3657e-8 Ry. The shifted 4x4x4 crystal-minus-off shift is
-1.5538989465824216e-5 Ha, versus the relayed rounded -1.554e-5 Ha:
absolute discrepancy 1.010534175784e-9 Ha. The symmetry evidence used LDA;
it validates the existing projector separately, not an LDA capability of
this PBE-only mode. The material totals used symmetry off. These are distinct
studies, not six new calculations combining all four active conventions.
Neither convention set is claimed to be more correct. No general cross-code
identity or derivative equivalence follows from these measurements.

## Defaults, metadata and restart

Absent `compat` or `mode="none"` keeps every default and arithmetic path
unchanged. Default metadata and checkpoint wire omission remain unchanged.
Mode results, CLI preflight and explain output identify `compat_mode` and all
four resolved values, including inactive/default density-group selection.

The mode label travels through SCF options and results. Exact continuation
compares it independently of the numerical options, so a manual equivalent
input without the mode cannot resume a mode-labelled exact checkpoint, or
vice versa. Explicit WarmStart is permitted subject to existing compatibility
checks. An old checkpoint missing the field loads as `none`; unknown values
reject. Symmetry checkpoint capture remains unsupported by the v1 envelope.
