# Crystal symmetry (P19A / MF11)

The opt-in `symmetry.mode="auto"` path detects operations in the original input
basis, selects the mesh-preserving subgroup, builds irreducible weights and
projects each newly built density into its invariant subspace. The default is
`off`; it retains the existing time-reversal setting and numerical path.
No external symmetry library or QE source was used.

## Contracts and algorithms

`dft_geometry::symmetry::detect` searches integer columns with lengths bounded
by the inverse lattice row norms. This is not a fixed [-1,1] matrix search:
non-reduced skew bases can require larger integers. A candidate must have
determinant +/-1 and preserve the lattice metric. Cartesian atomic tolerance
is `tolerance_bohr` (default 1e-7 Bohr); each fractional component is bounded by
`tolerance_bohr * norm(row_i(A^-1))`. Metric product allowance is
`tol*(norm(a_i)+norm(a_j))+tol^2` Bohr^2. Atom images are enumerated with these
bounds, rather than component-rounding a Cartesian nearest image in skew cells.
Species labels must match. Ambiguous mappings reject; accepted permutations
must be bijections. Rotations and wrapped translations have deterministic
lexicographic order. Closure includes translations and atom permutations.
The bounded search rejects more than 2,000,000 integer vector candidates per
column and tolerances with fractional component bounds >=0.1; it never returns
a silently incomplete operation list when that search bound is exceeded.

`reduce_mesh` takes a duplicate-free uniform full mesh. It selects only
operations that map every point onto that mesh (fractional tolerance 1e-12).
Its strict variant rejects any excluded operation. Time-reversal partners are
explicit antiunitary entries for unpolarized runs only. Representatives use the
first full-mesh index; weights are integer multiplicity/full-count and must sum
to one within 1e-14. Each full index has a representative, operation, integer
shift and multiplicity. The convention is
`k_full = (+/- R^-T) k_rep + reciprocal_shift`.

For `x' = R*x + tau`, the density projector is
`rho_sym(G) = mean_op[exp(-2*pi*i*G.tau) rho(R^T G)]`.
The corresponding Bloch phase is `exp(-i*(k+G).tau)` with physical vectors, or
`exp(-2*pi*i*(k_frac+G_miller).tau_frac)`. Density coefficients outside the
alias-free density sphere are zero. Missing rotated sphere entries reject.
The source indices, phases and three complex FFT buffers are prepared once;
application allocates no buffers, logs nothing and performs no file I/O.
Input densities that are not invariant within 1e-12 electrons/Bohr^3 reject,
including externally supplied `step_kpoints` densities. Both collinear channels
use the same unitary identity-spin projector. No magnetic antiunitary or
spin-swapping operations are enabled. Uniform fixed total moment is invariant;
local constraints, external fields and Hubbard terms are not supported and
unknown configuration tables reject instead of being ignored.

### Nonlinear quadrature is part of the invariant model

A real-space FFT quadrature must be permuted by the operations:
`N_i*tau_i` and `N_i*R_ij/N_j` must be integers. The first has an absolute
integer-distance allowance of 1e-8; the second is checked in integer arithmetic.
Otherwise nonlinear XC can break symmetry even when Fourier density phases are
correct. Automatic grids choose 2/3/5-smooth dimensions satisfying these rules,
first preserving independent axis lengths and then trying a common dimension.
Search stops at 512 per axis and fails explicitly for incommensurate origins.
An explicit noninvariant `basis.fft_grid` rejects. This can enlarge the automatic
grid only with symmetry enabled. Full/reduced comparisons MUST use the same grid.
The projector itself supports arbitrary fractional translations; the nonlinear
SCF quadrature restriction is a separate engine requirement.

An initial actual-Si test using the legacy 25^3 grid failed the density oracle:
1.3012983775327992e-8 Bohr^-3/2, despite delta E=3.7898573168604344e-12 Ha.
That grid cannot represent diamond's quarter translation as a permutation.
Selecting 32^3 on **both** sides resolves the failure without changing the SCF
or oracle tolerances. This is not a claim that 25^3 off and 32^3 auto have
identical discretization. An initial MgO test without the existing NLCC floor
also failed on negative core-density ringing; the final test uses the existing
1e-10 XC floor policy on both sides, with its original clipping-charge checks.

## Validation and provenance

Commands and raw output are stored under `tests/reference/symmetry/`. The
release engine test uses tightened density=1e-11 Bohr^-3/2,
energy=1e-12 Ha and eigen=1e-11 Ha tolerances. It tests full versus reduced SCF
and reduced T=1 versus T=4 exact energy/density bits for the shape/material set.
The geometry tests independently invert k maps with R^T, reconstruct Bloch
phases in real coordinates, check closure, multiplicities, shifted mesh rejection,
species separation, origin shifts and tolerance-edge ambiguous atoms.
The physics oracle directly evaluates a translated analytic field at
`R^-1(x-tau)` without using the production Fourier mapping; its maximum observed
error was 2.220446049250313e-15 electrons/Bohr^3. It also tests idempotence and
constant-mode preservation.

Actual material evidence uses external UPFs, never repository copies:
Si diamond 2 atoms, shifted and unshifted 4^3; Si cubic 8 atoms, shifted 2^3;
MgO rocksalt, PseudoDojo PBE/NLCC, shifted 2^3. These are finite-cutoff symmetry
oracles, not cutoff convergence or transferability claims. Synthetic cases cover
cubic, skew and slab shapes, local/nonlocal potentials, collinear moment and FD.
For Si 4^3: 48 crystal operations, 48 mesh-preserving operations unshifted and
12 shifted; irreducible counts are respectively 8 and 10. The 64 full points
have zero weight-sum error. Stored QE outputs explicitly say `No symmetry found`
and give k222=8 and k222s=4 (time reversal only); there is no stored symmetry-on
k444s QE output to independently confirm our 10-point count. No QE reference
number or QE speedup is claimed.

## Mac work and timing observation

`tests/reference/symmetry/mac-k444s.json` records complete input text, binary and
UPF SHA-256, three alternating full/TR-only/crystal runs, load and raw timings.
All use 30 Ry, 32^3 FFT, T=1 and the tightened tolerances above. There was no
concurrent build during timing. This is a Mac observation, not dev08 acceptance.

| Sampling | k points | Median wall s | Iteration operator applications | Certificate applications field | SCF steps |
|---|---:|---:|---:|---:|---:|
| Full | 64 | 6.647655625 | 25946 | 1962 | 14 |
| Time reversal only | 32 | 3.327696416 | 12973 | 981 | 14 |
| Crystal + time reversal | 10 | 1.106044584 | 4040 | 307 | 14 |

The observed wall ratio TR-only/crystal is 3.00865. Relative to the full mesh,
the reduced energy difference is 1.2434497875801753e-14 Ha and density distance
4.665499908802415e-16 Bohr^-3/2. The work-count reduction is the portable evidence;
these timings do not establish QE-relative performance. No default was changed.

## Limits and owner follow-up

No space-group naming, primitive-cell standardization, magnetic time reversal,
spin swap, SOC, fields, Hubbard terms, or local constraints are implemented.
The engine rejects symmetry checkpoint capture/import because the existing
checkpoint identity has no symmetry record; CLI auto is Classic-only. No
checkpoint schema is silently extended. Arbitrary incommensurate origins may
require `off` or a symmetry-compatible input origin; large integer searches fail
closed. Geometry detection alone does not prove a future added Hamiltonian term
invariant: adding one must extend the model-invariance checks.

Files outside this ticket's OWNED_FILES need owner edits: extend
`scripts/make_readme_figures.py` for this evidence and regenerate the README
validation block; add this attempt and its grid lesson to
`docs/performance/attempts-ledger.md`, `optimization-log.md` and `backlog.md`;
update the orchestrator's artifacts/dashboard using its authoritative session
registry. Existing README freshness tests can pass without including this new
result type. No permission to modify those files is inferred here.

For dev08, request the same stored inputs, T=1/4/8/16 on equal cores, pure and
MKL, binary/UPF hashes, density/energy equality, operator counts, and quiet-node
alternating runs. Request a new QE 7.6 symmetry-on k444s bundle with actual
operation/k counts, command, input and UPF hash; compare QE's best configuration
at equal delivered accuracy. Mac numbers must not be substituted for that run.
