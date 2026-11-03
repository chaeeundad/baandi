# Additional Directive: General-Purpose Scientific Coupling Architecture

**Project:** Rust DFT electronic-structure engine  
**Audience:** Claude Code (technical lead/reviewer) and Codex via MCP (implementation and independent review)  
**Status:** Architecture/design directive; do not interrupt active k-point validation  
**Priority:** Preserve future coupling flexibility without expanding the DFT solver's immediate implementation scope.

## 1. Mission

Design a **general-purpose, physically rigorous, extensible coupling interface** for the Rust DFT engine. The engine should serve as a reliable electronic-structure component in larger scientific workflows while remaining a focused Kohn-Sham DFT codebase.

**DFT-CES and DFT-CES2 are reference use cases, not architectural dependencies or the sole target.** Additional use cases include classical and polarizable QM/MM, electrostatic embedding, explicit/implicit/hybrid solvation, electrochemical interfaces, applied fields, charge-constrained calculations, external self-consistency loops, ML/MLIP data generation, and continuum/FEM parameterization. Treat these as design stress tests, not features to implement now.

The objective is **physics-aware interoperability**, not a generic plugin system or a monolithic multiscale simulator.

## 2. Non-negotiable Scope and Scheduling

- **Continue existing k-point implementation and validation unchanged.** Do not introduce refactors that destabilize it.
- Finish the baseline, independently validated plane-wave Kohn-Sham SCF solver before integrating coupling functionality.
- For now: inspect repository, design interfaces, document physical contracts, create a phased roadmap and focused test plan. **No major feature implementation without explicit approval.**
- Do not implement MD, MLIP inference/training, force fields, FEM solvers, solvent models, or external orchestration inside `dft-core`.
- Do not make DFT-CES(2), LAMMPS, ASE, Python, REST, gRPC, or any specific framework a required dependency of the core.
- Keep mathematical/numerical kernels independent of serialization, networking, file formats, and external workflow scheduling.

## 3. Architectural Separation

Propose logical boundaries, adapted to the repository's existing structure (do not create crates prematurely):

1. **DFT Core:** structures, reciprocal/real-space grids, Hamiltonian, eigensolver, SCF, energy, forces, stress, electronic-state representations.
2. **Coupling Contracts:** typed external physical inputs and observable outputs, immutable context, units, frames, validity metadata, compatibility checks.
3. **Coupling Runtime (future):** explicit outer iteration/handshake, state exchange, checkpoint/restart, synchronization, convergence policy and failure handling.
4. **Adapters (future, out of core):** file interchange, Python/ASE, QM/MM driver, DFT-CES-family driver, MLIP datasets, FEM material-property exporters, remote execution.
5. **ESLE (separate):** multiple electronic solutions, branch identity, state comparison, exploration policy and provenance.

A minimal `EnergyForceCalculator` interface is useful for consumers but **must not dictate internal SCF or Hamiltonian architecture**. Consider separate interfaces for calculation requests/results, field injection, observables and state continuation rather than one oversized trait.

## 4. Required Physical Input Contracts

Design **capability-gated**, typed inputs, with explicit availability and unsupported-feature errors:

- External scalar potential on a 3D grid, clearly distinguishing **electrostatic potential (volts / energy per charge)** from **electron potential energy (Hartree or eV)**, including electron-charge sign.
- External charge distributions (point charges, continuous density, or gridded charge) as separate representations; specify whether and where Poisson solving occurs.
- Optional external electric fields, time-independent perturbations and future generalized embedding operators. A uniform electric field under periodic boundary conditions is **not** equivalent to an arbitrary periodic scalar field; document Berry-phase/slab-specific restrictions rather than claiming generic support.
- Cell and grid geometry: lattice vectors, origin, ordering/strides, reciprocal-space convention, periodic axes, interpolation, FFT normalization and aliasing constraints.
- Electron count, total charge, spin treatment, occupations, XC/pseudopotential/basis/k-point metadata and boundary-condition compatibility.
- Optional model-specific terms through explicit extension points **only when a consistent energy functional and derivative contracts can be defined**. Do not force nonlocal or density-dependent interactions into a scalar-potential abstraction.

Clarify whether an external contribution is fixed, density-dependent, geometry-dependent, time-dependent or updated by an external driver. These distinctions affect variational consistency and forces.

## 5. Required Output Contracts

Specify versioned outputs with quality and provenance, when actually supported by the underlying solver:

- Total energy and **decomposed energy terms**, with a precise accounting policy for external coupling terms and double counting.
- Electron density, spin density where supported, electrostatic potential and fields, with clear distinction between electronic-only, ionic-only and total quantities.
- Forces and stress only when corresponding derivatives are implemented and validated; include whether external-field and coupling contributions are present.
- Eigenvalues, occupations, Fermi level, SCF residuals and convergence diagnostics where available.
- `ElectronicStateId`, `BranchId` (ESLE), initialization lineage and compatibility metadata; never imply one SCF solution is the guaranteed ground state.
- Units, coordinate frame, reference potential/gauge, boundary conditions, numerical tolerances, provenance, and warnings/validity domain.
- Large arrays via a versioned binary or array-capable format, with lightweight JSON metadata as needed. Do not mandate JSON for dense 3D fields or wavefunctions.

**FEM-facing outputs:** report microscopic properties and conditions of validity. A 0 K elastic tensor is not a general finite-temperature constitutive law; homogenization, microstructure, finite-temperature effects and constitutive modeling belong to separate workflows.

**MLIP-facing outputs:** preserve state labels and physical conditions so distinct electronic solutions at identical atomic coordinates are not silently merged into a single-valued training target.

## 6. Coupling Semantics and Lifecycle

Specify a driver-neutral sequence:

1. Negotiate capabilities and validate physical context.
2. Receive structure, coupling terms and optional compatible electronic-state checkpoint.
3. Solve the requested electronic problem to declared convergence criteria.
4. Return observables, accounting terms, state identity and quality metadata.
5. Let an external driver update the environment and determine whether to repeat.

The DFT core must **not** assume the outer loop is DFT-CES, MD, a continuum solver or ML optimization.

Support both **one-way** (external environment fixed) and **two-way** (environment responds to DFT outputs) workflows in the future. Keep inner electronic SCF convergence distinct from outer coupled-loop convergence. Specify deterministic restarts, idempotent request identity, partial/failure statuses, cancellation and resource limits as future runtime requirements.

For density-dependent embedding, clearly document functional derivatives, variational stationarity, free-energy versus total-energy interpretation, and how consistent forces/stress would be obtained. Do not report an arbitrary external potential as though it came from a known total-energy functional.

## 7. ESLE Integration

- A coupling request can optionally identify an initial electronic state or target branch.
- Track branches across external-environment updates while preserving independent multi-start exploration.
- Record external-environment identity in state compatibility and provenance.
- Distinguish **electronic metastability** from multiple self-consistent **coupled QM/environment solutions**.
- Never rank coupled states using DFT energy alone when environmental free-energy terms differ.
- Do not let ESLE orchestration leak into the Hamiltonian or SCF numerical kernels.

## 8. Reference Use Cases to Test the Design (Not Implement Now)

A. **Fixed electrostatic embedding:** inject a well-defined periodic scalar electron-energy potential and export electron density; verify energy accounting and finite-difference response.

B. **Mean-field explicit-solvent coupling (DFT-CES/DFT-CES2 family):** external driver exchanges averaged solvent environment and QM electrostatic observables; demonstrate how a richer embedding model would be accommodated without core modifications.

C. **Polarizable or density-dependent QM/MM:** environment responds to QM density; require a documented coupled energy/free-energy functional and derivative consistency.

D. **Implicit solvation:** density-dependent reaction field with its own solver and energy functional, connected through an adapter.

E. **MLIP/active-learning workflow:** export structures, energies, forces, stress, state labels and numerical quality, without embedding ML inference.

F. **FEM workflow:** export properly framed, unit-checked, provenance-rich microscopic material properties, not unqualified macroscopic constitutive laws.

G. **ESLE under an external environment:** track competing electronic branches while external fields or solvent conditions change; distinguish fixed-environment energies from coupled free energies.

For each use case, document data exchanged, required core capability, external responsibilities, physical consistency constraints and current gaps.

## 9. Validation and Acceptance Criteria

Before any implementation, propose tests covering:

- Units and electron-charge sign (including an analytically solvable constant potential shift).
- Real/reciprocal grid transforms, origin, interpolation and periodicity.
- Electron-number conservation and density normalization.
- External-potential energy accounting and double-counting prevention.
- Finite-difference energy derivatives versus forces/stress **only for well-defined, supported functionals**.
- Consistency across restart and serialization round-trips.
- Clear rejection of incompatible k-point, cell, spin, charge or pseudopotential metadata.
- Separation of inner SCF convergence and outer coupled convergence.
- Multiple-state provenance and branch tracking without assuming global minimum.
- CPU/GPU backend neutrality for physical contracts.

Use established scientific conventions and interoperable formats where practical; evaluate existing ecosystem approaches before proposing custom standards. Do not copy incompatible licenses or code.

## 10. Agent Workflow and Deliverables

**Claude Code:** remain technical lead; inspect current code, coordinate, review scientific assumptions and preserve existing milestones.

**Codex via MCP:** conduct an independent architecture/design review and, if useful, delegate focused reviews of physical consistency, serialization and interoperability. **Do not begin large implementation changes.**

Produce repository-appropriate documents (adjust paths to existing conventions):

1. `docs/architecture/scientific-coupling.md` — module boundaries, contracts, lifecycle and example typed interfaces.
2. `docs/research/coupling-physics.md` — energy/free-energy, charge sign, double counting, boundary conditions, forces/stress and limitations.
3. `docs/architecture/coupling-use-cases.md` — capability matrix for cases A–G, including DFT-CES(2) as examples only.
4. `docs/roadmap/coupling.md` — gated roadmap: design now → baseline SCF → fixed-field validation → observables/derivatives → external drivers → coupled state exploration.
5. An ADR describing which interface commitments must be made now and which should be deferred.

**First response required:** concise repository-grounded gap analysis, proposed minimum stable interfaces, compatibility with current k-point work, high-risk physics assumptions, and explicit decisions requiring approval.

## 11. Definition of Success

The architecture succeeds if new scientific coupling methods can be implemented **primarily in external drivers/adapters**, while the DFT solver remains physically correct, independently testable and focused. Future integration must not require introducing MLIP, MD, FEM or any single solvent methodology into the numerical core.

**Guiding principle: Make the DFT engine composable, not monolithic; physics-aware, not merely format-compatible.**
