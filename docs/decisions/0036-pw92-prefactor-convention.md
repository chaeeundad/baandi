# ADR 0036: explicit PW92 prefactor convention for unpolarized PBE

- Date: 2026-11-07
- Status: **Accepted (option 1), 2026-11-07**. The user accepted the recommendation:
  `modified` stays the default and `original` stays an explicitly named option for
  QE-equivalence comparisons.
- Ticket: P27B-PW92-ORIGINAL-OPTION.

## Facts and sources

PW92, [Phys. Rev. B 45, 13244 (1992)](https://doi.org/10.1103/PhysRevB.45.13244),
Eq. (10) and Table I, prints the unpolarized correlation prefactor A=0.031091.
The existing xcx-backed implementation uses A=0.0310907 (modified convention).
The other parameters are 0.21370, 7.5957, 3.5876, 1.6382 and 0.49294.
Neither convention is claimed to be more correct. The modified value is closer
to the high-density coefficient, but the original printed fit used rounded
parameters; a comparison must identify its convention.

[PBE, Phys. Rev. Lett. 77, 3865 (1996)](https://dft.uci.edu/pubs/PBE96.pdf),
Eqs. (7) and (8), defines the correlation gradient correction using an
independent gamma=(1-ln(2))/pi^2. That gamma is unchanged. Full-precision beta
0.06672455060314922 and exchange parameters remain those of ADR 0018.
The production scalar implementation and its analytic chain rule were written
from these formulas and the existing independent P27A diagnostic, without QE
source access or copying. The architect's own-word study note supplied the
convention hypothesis, not source code. This ticket authorizes the narrowly
scoped in-tree correlation exception to ADR 0018; exchange still uses xcx.

P27A predicted fixed-density changes -2.1002508994e-6 Ry (si-20) and
-6.8705824638e-6 Ry (mgo-40). New self-consistent changes are
-2.1002509278e-6 Ry and -6.8705848548e-6 Ry. All six material calculations
converged. See the [dated validation record](../validation/xc-prefactor-diagnostic.md#p27b-2026-11-07-explicit-convention-and-self-consistent-results)
for all energies, gates, commands and limitations. No new QE run was made.

## Options and recommendation

1. Keep modified as the default and expose explicitly named original for
   QE-equivalence comparisons. **Recommended:** preserves existing results and
   makes a convention difference visible and reproducible.
2. Change the default to original. This would change existing energies and
   densities and requires a separate user decision and regression migration.
3. Keep modified only. This preserves results but leaves an avoidable convention
   mismatch in comparisons.

## Implemented contract

```toml
[xc]
functional = "pbe"
pw92_prefactor = "original" # optional; omission means "modified"
```

Only the literal strings `modified` and `original` are accepted. Original is
supported only for unpolarized PBE; LDA and collinear PBE reject it at the CLI
and engine construction boundaries. This does not override UPF functional
compatibility or alter the density threshold, floor, NLCC, exchange, gamma or
any convergence tolerance. No spin-PBE extension is implied.

`ScfOptions.pw92_prefactor` defaults to Modified. Default evaluation delegates
to the unchanged functional/library path. Original has separate correlation
assembly with the same spectral derivative, screening and quadrature. The
potential differentiates n*(ec+H) analytically, including both density and
sigma dependence. Forces use this XC potential; stress also selects its
explicit gradient-metric derivative. The private derivative model identity
includes the convention to reject stale certificates with changed options.

SCF results carry the enum; CLI `scf.pw92_prefactor` and scalar top-level
`pw92_prefactor` record original. Omission denotes modified and preserves
legacy JSON shape. Finite-temperature outputs retain the convention in their
SCF options metadata. Existing checkpoint JSON omits modified and decodes an
absent field as modified. Original is serialized explicitly. Exact restart
compares the option in both directions; changing it requires explicit WarmStart
with a new calculation identity. Existing NLCC/floor checkpoint restrictions
remain applicable.

## Decision still open

User confirmation of recommendation 1, or a separately scoped change of
convention policy. No decision is required to retain the currently specified
modified default. Validation here establishes consistency and the measured
Si/MgO effect, not universal equality with QE or a more accurate functional.
