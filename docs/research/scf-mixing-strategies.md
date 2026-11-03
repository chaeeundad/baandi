# Density mixing strategies (P10D)

Status: optional solver strategy integrated into SCF and CLI; default unchanged.
The orchestrator approved minimal module-export and CLI wiring edits.
No QE source was read. No metallic convergence is claimed; integer occupations
and the engine's insulating-gap test remain mandatory.

## Audit of the current implementation

`DensityMixer::mix(x, F(x), mean)` consumes real, uniform-grid electron density
in electrons/Bohr^3; R=F(x)-x. `LinearMixer` proposes x+beta R. Beta must be
positive and finite, but may exceed one. `PulayMixer` uses type-II Anderson:

    gamma = argmin ||R - DR gamma||_2
    next = x + beta R - (DX + beta DR) gamma.

It stores a FIFO of m+1 input/residual pairs, permitting m differences; defaults
are beta=0.3 and m=6. Differences and Gram matrices are rebuilt each call.
The faer-backed Hermitian Gram eigensolve discards eigenvalues <=1e-12 times
the largest. This is a spectral pseudoinverse, not ridge regularization; forming
normal equations squares the condition number. Singular histories are supported.
History resets explicitly or when shape/mean changes. Existing storage, cloning,
small eigensolve allocations and result Vec allocations are unchanged.

Both mixers apply two uniform mean shifts to enforce prescribed Ne/Omega to
summation roundoff. They neither clip samples nor enforce positivity themselves.
SCF damps a negative Pulay proposal along its displacement from the old input:
theta=min(1, min_negative 0.9*x/(x-proposal)); it then resets history. Linear
mixing has no feasibility safeguard. This convex operation conserves charge
and Fourier support. Density below -1e-14 Bohr^-3 is rejected by the controller.

SCF reports sqrt(Omega/Ngrid * sum R^2)/Ne, with units Bohr^-3/2 when electron
count is dimensionless. It uses the **unfiltered** residual. Three successive
residuals above divergence_factor*max(best_residual,density_tolerance) give
`Diverged`; defaults are factor=20, patience=3. Nonfinite densities, arithmetic,
and mixer errors also fail explicitly. A finite iteration limit gives
`NotConverged`, not success. No restart currently tries to rescue that divergence.

On the same iteration, the density, energy change, occupied eigen residual,
requested eigensolve tolerance and charge gates must all pass. A second map
from the candidate output, solved at 0.1*eigen_tolerance, must also pass the
fixed-point density and occupied eigenvalue certificate. Mixing must not weaken
these tests or use a small filtered norm as evidence of convergence. See ADR 0009.

## Published candidates and interpretation

* **Kerker:** damp long wavelengths with P(G)=G^2/(G^2+q0^2). q0 is an inverse
  screening length, not a dimensionless beta. This removes the model Coulomb
  singularity and leaves high-frequency components close to unchanged.
  [Kerker, Phys. Rev. B 23, 3082 (1981)](https://doi.org/10.1103/PhysRevB.23.3082).
* **Modified Broyden:** Johnson's limited-memory inverse-response updates use
  information from earlier iterations to accelerate a nonlinear fixed point.
  This is a distinct candidate, not an alternative name for the implementation
  below. It needs separate history weighting and conditioning validation.
  [Johnson, Phys. Rev. B 38, 12807 (1988)](https://doi.org/10.1103/PhysRevB.38.12807).
* **Preconditioned Pulay:** combine a response model with multisecant history.
  The norm used to fit history and the preconditioner applied to a step are
  separate design choices. We explicitly choose Anderson on the filtered map;
  fitting raw R then filtering only the correction would be a different method.
  **Density versus potential:** density mixing gives a direct charge constraint;
  potential mixing requires fixing the constant-potential gauge and reconstructing
  density before applying charge and fixed-point tests. It would require more
  SCF integration than this ticket permits. A broad comparison and discussion of
  failure cases is [Woods, Payne and Hasnip, J. Phys.: Condens. Matter 31,
  453001 (2019)](https://doi.org/10.1088/1361-648X/ab31c0).
* **Local response:** spatially varying elliptic preconditioners can distinguish
  vacuum from screened regions. A single global q0 cannot represent both well;
  uniform Kerker can suppress useful long-wavelength insulating modes. A local
  model adds coefficient selection, an elliptic solve and boundary/null-space
  handling. This is a promising separate strategy for heterogeneous slabs,
  not implemented here. [Lin and Yang, SIAM J. Sci. Comput. 35, S277–S298
  (2013)](https://doi.org/10.1137/120880604).
* **Adaptive steps and restart:** measured residual growth can inform a trust
  region or history rejection, but density feasibility alone does not prove
  residual reduction. Multisecant safeguards are studied by [Marks and Luke,
  Phys. Rev. B 78, 075114 (2008)](https://doi.org/10.1103/PhysRevB.78.075114).
  Our strategy retains the existing feasibility reset and explicit divergence
  error and adds a recorded raw convex recovery for an infeasible history-free
  filtered step; it does not retry the SCF with fitted parameters. Adaptive beta, rollback and
  Johnson updates are deferred until a broader measured test set exists.

## Implementation contract

`mixing/kerker.rs` defines a `ResidualPreconditioner` and `PreconditionedPulay<P>`
implementing the existing `DensityMixer` interface. The wrapper evaluates raw R,
filters it, then passes x and x+PR to the unchanged Pulay implementation.
Consequently its least-squares metric is ||PR-D(PR)gamma||_2 and its history
contains the filtered residual. Adding x then subtracting it inside legacy
Pulay incurs ordinary cancellation roundoff; no exact arithmetic claim is made.
Raw residual diagnostics are computed before this operation.

`Kerker<T>` accepts a C-order grid, reciprocal columns B (including 2*pi), q0
in Bohr^-1 and a reusable transform adapter. The adapter contract is forward
1/N, inverse unnormalized. Keeping the adapter at the caller boundary avoids
adding an ndarray dependency to solvers. Weights and complex scratch allocate
only at construction; filtering itself allocates no buffers. Wrapper residual
and mapped buffers are reused. The inherited Pulay allocations remain and are
a known limitation, not an accepted allocation-free SCF optimization.

At zero frequency the filtered coefficient is assigned exactly zero. Legacy
mean correction removes summation error after extrapolation. Skew cells use
|B*m|, not a diagonal-cell approximation. On even grids a Nyquist coefficient
aliases positive/negative indices: weights average over all signs of Nyquist
coordinates. This makes conjugate-paired weights equal in a skew cell; simply
using a single signed Nyquist vector can create a complex inverse transform.
Non-Nyquist modes use the exact Kerker expression. No spectral support is added.
Invalid shapes, nonfinite/singular reciprocal matrices and invalid q0 are errors.

Diagnostics expose raw Euclidean residual norm, raw contraction ratio, number of
history differences presented to the fit, and reset count. History size is not
numerical rank. The contraction ratio is absent on the first call, after reset,
or after a zero residual. The solver does not certify convergence: its caller
must retain physical gates. Mean changes reset history; a shape change requires
a new grid-specific preconditioner. Errors propagate to the controller.

## Recommendation and scope

Keep the current default. Trial Kerker only as an explicit choice when residuals
show long-wavelength sloshing; compare raw residual histories and energies at
unchanged tolerances across q0 and beta. Larger q0 damps more wavelengths and
can slow insulators and vacuum regions. No universal parameter is inferred from
a silicon result or from the analytic dielectric model. Use the existing linear
strategy for transparent diagnostic experiments; it is not guaranteed stable.
The synthetic model isolates a known response and alone cannot establish
material performance. The separate SCF matrix tests synthetic pseudopotentials,
several k meshes, Si and elongated/vacuum cells. See the comparison report for
case-specific successes and failures; no metallic conclusion follows.

## Controller and CLI integration

`MixerHistory` is an additive extension trait implemented for the unchanged
legacy mixers and the new wrapper. It exposes stored-pair count, avoiding a
change to `DensityMixer` itself. SCF reads the actual post-proposal history size,
records raw residual/current-to-previous contraction and cumulative feasibility
resets, and includes Kerker–Pulay in the existing positivity safeguard. Final
iterations which need no proposal have a null history-used field. No eigen
schedule, divergence threshold, energy gate or certificate arithmetic changes.
The closure in the new selection branch captures an FFT plan and two arrays;
it uses the adopted numerics adapter, with no per-filter allocation. FFT errors
propagate through the existing mixer-to-divergence path.

The CLI accepts `mixing="kerker-pulay"`, required `mixing_q0_inverse_bohr`, and
optional `mixing_history` (six by default). No global screening length is chosen
silently. Effective parameters and diagnostics are included in results metadata.
See [CLI documentation](../cli.md) and ADR 0025. The detailed material matrix is
an ignored cross-layer process test, so the solver crate gains no engine, array
or serialization dependency. Generated synthetic UPFs are temporary evidence,
not committed data; real Si UPFs are read in place.

## Measured boundary-stagnation failure and recovery

Initial material runs exposed a failure of merely wrapping filtered Pulay with
the existing feasibility reset. In Si Gamma and synthetic elongated/vacuum cells,
a history-free filtered step pointed negative at a nearly zero-density sample.
Convex feasibility fractions shrank by about a factor of ten per iteration;
restarting filtered history repeated the same outward direction. At q0=1 the
Si raw residual stalled at 0.02932 Bohr^-3/2 (98 resets in 100 iterations), and
the synthetic elongated/slab cases also exhausted 150 iterations. This was
not an FFT-normalization or charge error and was not fixed by tolerances.

The new strategy now detects a negative proposal with exactly one stored pair
(no differences) and substitutes x+min(beta,1)*R using the existing charge-corrected
LinearMixer. For nonnegative input/output this is a convex density combination,
so no density clipping, privileged sample or material-specific threshold is
needed. It preserves the valid sampled (x,PR) history pair. Subsequent history
extrapolation still uses the unchanged controller safeguard. Recovery is recorded
in solver diagnostics and each SCF record as `mixing_raw_recovery`; it is not
counted as a reset unless the controller actually resets history. This behavior
is restricted to the new strategy; legacy Pulay arithmetic is unchanged.

The decision is motivated by the same observed failure across independent Si,
synthetic nonlocal elongated and synthetic local vacuum cases, and by the convex
feasibility argument. A dedicated solver test starts on the density boundary and
compares the recovery to an independent invocation of legacy linear mixing.
No claim of universal convergence follows. Very large q0 still produces honest
stagnation/NotConverged because convergence is judged on the raw residual.
