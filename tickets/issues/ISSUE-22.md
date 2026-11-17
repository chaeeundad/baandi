# GitHub issue #22: Research: prototype occupied-subspace energy minimization as an alternative to nested SCF/Davidson

Reviewed at `1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4` (2026-11-17). This is a static review of the implementation and committed experiment records; no new BAANDI/QE timing run was performed in this review. Related tracking issue: #13.

**Suggested priority:** P3 — default-off architectural research.

### Hypothesis

The larger architectural opportunity is to optimize the occupied subspace and density together, avoiding repeated approximate diagonalization of successive potentials. A direct minimizer is a separate experimental driver, not a retuning of the existing fixed-block Advance path.

Start with the currently supported CPU real-Gamma, fixed doubly occupied, norm-conserving setting. Optimize the true orbital energy over orthonormal occupied X. The tangent residual direction is proportional to

```text
G = H[rho(X)] X - X (X^T H[rho(X)] X)
```

and density is invariant under occupied-space orthogonal rotations. Keeping this invariance avoids solving for irrelevant rotations within equally occupied states. Spin/occupation factors in the actual energy gradient must follow BAANDI's conventions and be verified, not inferred from this schematic residual.

### Research basis, with the negative evidence included

Altmann, Peterseim & Stykel, **J. Sci. Comput. 101, 6 (2024)**, [DOI/full text](https://doi.org/10.1007/s10915-024-02612-3), derive inexact Riemannian Newton methods and test plane-wave Kohn–Sham problems. Their pentacene example uses fewer Newton outer iterations but similar total times across methods; the difficult graphene example benefits more. Thus fewer iterations do not imply lower wall time.

A more recent [CPC 312, 109596 (2025)](https://doi.org/10.1016/j.cpc.2025.109596) study formulates direct minimization on a complex Stiefel manifold, but in localized non-orthogonal bases. Its abstract supports investigating the algorithm family, not transferring a benchmark or implementation to BAANDI. The 2024 full text is the primary algorithm reference here.

### Proposed prototype

1. Expose an **energy + tangent-gradient** primitive for arbitrary orthonormal occupied orbitals. Reuse the explicit orbital functional in [energy.rs](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-engine/src/energy.rs#L1), density construction and H application. Do not use an SCF band-sum shortcut with stale Ritz values as the energy of a nonstationary trial state.
2. Validate directional derivatives of the complete energy, including kinetic, ionic-local, nonlocal, Hartree and XC terms, along retracted tangent directions. Validate invariance under occupied rotations and normalization in the packed real-Gamma representation.
3. Start with a preconditioned limited-memory quasi-Newton or Riemannian conjugate-gradient pilot plus a safeguarded line search. This tests the new dataflow with only energy/gradient evaluations. Then evaluate an inexact Newton variant only if a matrix-free Hessian-vector action is affordable and independently checked.
4. Reuse geometry, projectors and FFT plans, and retain valid H products at each accepted state. Count every trial energy, line-search rejection, H application, reduction, orthogonalization and final canonicalization.
5. Obtain empty states and canonical eigenvalues in an explicit terminal stage required by the existing output/verification contract. Include its full cost. For unequal/fractional occupations, the simple occupied-rotation invariance no longer suffices; decline until an occupation-aware formulation is separately validated.

### Why BT3 does not settle this question

BT3's Advance still uses an inexact SCF map and pays for completing solves, mixer behavior and failed verifications. Here the accepted step is selected using the orbital energy and its consistent gradient. That is a changed mathematical problem formulation, not “try fewer Davidson steps again.”

### Go/no-go

Compare from the same atomic/prolonged start and at the same full-cutoff final accuracy. Use a small Si case plus a nonsymmetric/non-Si case before an MPI implementation. Keep #15's spectral-completeness concern: a small stationary residual does not prove the lowest occupied subspace was found. Final ground-state/verification checks and the existing strict reference remain necessary.

The pilot is a research branch. Abandon or narrow its scope if line searches, Hessian solves or terminal empty-state recovery consume the saved work. No QE parity claim precedes an equal-condition measurement.

### Validation and adoption

Use the repository's current two-lane fast-loop protocol. Start with one short case and an identity/correctness smoke; use full-SCF Stage 1 for any trajectory change. Final acceptance uses the designated MPI benchmark host, the current default stack, fixed-point verification included, the existing accuracy classes and same-model/same-rank QE comparisons. Measurable small improvements remain eligible under the adoption rule; do not impose a new minimum percentage. Report failed or negative experiments as well as wins.
