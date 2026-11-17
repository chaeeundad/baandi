# GitHub issue #21: Research: response-based SCF preconditioning with a measured Sternheimer cost budget

Reviewed at `1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4` (2026-11-17). This is a static review of the implementation and committed experiment records; no new BAANDI/QE timing run was performed in this review. Related tracking issue: #13.

**Suggested priority:** P3 — default-off algorithm research.

### Question

Can a small, refreshed approximation to the physical density response remove enough expensive SCF/verification maps to repay its construction? This changes the **outer nonlinear solve** rather than only speeding Davidson kernels.

For `R(rho)=F(rho)-rho`, the local dielectric operator is `D=I-chi0*K_Hxc`. A Newton-like density correction solves `D delta_rho = R`. The needed action of `chi0` can be computed without a large explicit dielectric matrix.

### Research evidence and limits

- Das & Gavini, **Phys. Rev. B 107, 125133 (2023)**, [DOI](https://doi.org/10.1103/PhysRevB.107.125133), [author manuscript](https://arxiv.org/html/2211.07894v2): LRDM builds/accumulates a low-rank dielectric approximation. Their efficient response construction is coupled to a finite-element Chebyshev-filtered subspace. It is not a drop-in plane-wave Davidson implementation, and its difficult heterogeneous systems do not establish a benefit for small Si cells.
- Herbst & Sun, **SIAM J. Sci. Comput. 48, B571–B597 (2026; online July 9)**, [DOI/full text](https://epubs.siam.org/doi/full/10.1137/25M1759033): inexact GMRES controls nested Sternheimer accuracy through density-response error estimates. Its reported speedup concerns **linear-response calculations**, not ground-state SCF. This is a basis for controlling the extra response cost, not a transferable SCF speedup claim.

### Why this is a new premise

The ledger's A31 Kerker/Pulay experiment and LADDER-SECANT history transfer are not this method. Uniform metallic screening and coarse-grid secants can misrepresent semiconductor/heterogeneous response. BT3 also showed that residual growth can be a Pulay trajectory property rather than inner-solver noise. None of those results licenses another blind mixing-parameter sweep.

### Bounded experiment

1. **Offline response audit first.** At a few captured early/mid/late Si32 and non-Si states, implement a matrix-free `chi0 delta_V` action with charge conservation. For fixed integer occupations, solve projected Sternheimer equations in the complement of the occupied subspace. Keep occupied-occupied gauge rotations out of the density response. Include Hartree and the supported XC derivative; start with validated LDA, record unsupported cases.
2. Check response against central finite differences of a tightly converged SCF map over several perturbation amplitudes. Finite differences are a validation oracle, not the production preconditioner. Ensure eigen-solve noise is below the finite-difference signal.
3. Test one small response subspace aligned with actual residual history, with an explicit refresh/invalidation criterion. Use regularization and a fallback to the present mixer if the model predicts the wrong response or becomes unstable.
4. Choose inner response accuracy by its contribution to the density-step error. Do not compare a residual norm directly with a Ritz-change tolerance; do not treat a heuristic gap estimate as a proved bound.
5. Track response H applications, FFTs, reductions, added memory, rejected steps, SCF maps, failed verification maps, and **total** time. A lower outer iteration count alone is not success.

Break-even test: if the response costs C and removes k maps of measured cost M, require k*M > C plus all added control overhead. Si already converges in relatively few maps, so the method may primarily benefit heterogeneous/small-gap cases. Report that limitation if observed; do not force it on the default Si path.

### Deliverable and acceptance

A short design note, validated response primitive, and one default-off CPU experiment; no immediate default change. Reuse the same response API later for verification prediction or derivative calculations only if separately justified. Do not re-enable the rejected fixed-block Advance strategy on the strength of a paper alone.

### Validation and adoption

Use the repository's current two-lane fast-loop protocol. Start with one short case and an identity/correctness smoke; use full-SCF Stage 1 for any trajectory change. Final acceptance uses the designated MPI benchmark host, the current default stack, fixed-point verification included, the existing accuracy classes and same-model/same-rank QE comparisons. Measurable small improvements remain eligible under the adoption rule; do not impose a new minimum percentage. Report failed or negative experiments as well as wins.
