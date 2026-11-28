# Extensible BAANDI Research Program — Decision Log

- **Program**: Extensible Electronic Structure and Multiscale Simulation Engine (Tracks A, B, C)
- **Location**: `docs/research/extensible-dft/`
- **Governance (Directive Section 10.1)**:
  - Literature review, mathematical derivation, isolated numerical prototypes, research branches, benchmarks, profiling, and failure-driven alternative exploration run **autonomously** within an approved research scope.
  - **Research-lead approval / review is required ONLY for**: (1) initial track selection (Gate 0 → Gate 1 entry), (2) major scope changes, (3) major core-architecture changes, and (4) production integration (Gate 5).
  - Within an approved scope, Gate-to-Gate progression does **not** require separate human sign-off at every Gate boundary; every Gate decision (`Go / Revise / Hold / Stop`), tolerance, and failure analysis is logged here transparently.

---

## Standing Program Rules

1. **No Unverified Citations (Section 2.1)**: Only papers and DOIs directly verified via primary lookup may be cited; any unverified placeholder must be marked `[Citation needed]`.
2. **Zero Pure-PW Performance Regression (Sections 2.3 & 6.4)**: Extensions must not add allocations, I/O, branches, or extra MPI collectives inside `DistributedRealHamiltonian::apply_block`, `DistributedComplexHamiltonian::apply_block`, `ref_davidson`, or `ref_davidson_complex`. When disabled (`Option::None`), pure-PW execution must remain bit-identical.
3. **Production Contract Preservation (ADR 0047)**: All production extensions target `refscf/` and `dft-physics` (or dedicated opt-in crates). Frozen `oracle/` default arithmetic is never altered, and `cargo check/clippy -p dft-engine --no-default-features` must stay green.
4. **No Premature Universal Abstractions (Section 6.6)**: Do not introduce a monolithic multi-track "universal engine trait" or large core refactoring before concrete vertical implementations prove the minimal shared boundary.
5. **Independent Verification & Licensing Hygiene (Sections 7 & 9)**: Independent verification (never relying solely on the implementing agent's self-check); never relax tolerances post-hoc to manufacture a pass; never copy GPL/LGPL source code (e.g., QE, Environ, Wannier90, ABINIT, Yambo [GPL]; ASE [LGPL-2.1+]). Implement physics kernels cleanly from published mathematical formulations under `[workspace.lints.rust] unsafe_code = "deny"`.

---

## Decision Entries

### D-001 — Gate 0 Scope Boundary (2026-11-28)
- **Status**: Executed (independently audited by `Gate 0 Quality Auditor` subagent)
- **Decision**: Execute Gate 0 strictly as repository assessment, codebase/literature triage answering all 10 Directive Section 12 questions, and track prioritization without production Rust code changes, large refactoring, or premature 12-item theory specifications across all three tracks.
- **Evidence**: [repository-assessment.md](repository-assessment.md), [triage-report.md](triage-report.md).

### D-002 — Anchor New Seams in `refscf/` Rather Than Frozen `oracle/` (2026-11-28)
- **Status**: Adopted (consequence of ADR 0047)
- **Context**: Nov 2026 coupling documents (`docs/architecture/scientific-coupling.md`, `coupling-codex-review.md`, ADR 0017) referenced `ScfProblem` in `scf.rs`, which is now frozen under `oracle/` (`strict-oracle` feature).
- **Decision**: All future local-potential, solvation, external-field, and grand-canonical hooks target the production `refscf/{api, driver, k_driver, electronic, map, k_map, electronic_map}` and `refstart.rs` seams, preserving ADR 0017's sign/role contracts while obeying ADR 0047's single-path rule.

### D-003 — Retain Typed Electrostatic Semantics from ADR 0017 as Shared B/C Foundation (2026-11-28)
- **Status**: Proposed for Gate 1
- **Decision**: Implement explicit semantic types distinguishing electrostatic potential $\phi$ (Ha/e), electron potential energy $U_{\mathrm{el}} = -\phi$ (Ha), positive electron density $n$, and signed charge density $\rho$ before wiring either Track B solvation or Track C electrostatic embedding.

### D-004 — Defer Track A (Mixed Basis / PW+AO / ALB) Behind Track B & Track C (2026-11-28)
- **Status**: Recommended (pending user approval)
- **Rationale**: Track A requires replacing BAANDI's orthonormal PW eigenproblem ($S = I$) and `KeptHx` warm-start architecture with a generalized $H\psi = \epsilon S\psi$ solver, cross-density quadrature, and Pulay forces/stress. Doing Track A first carries the highest core-intrusion and performance-regression risk with the lowest cross-track reuse.

### D-005 — Separate Mechanical vs Electrostatic QM/MLIP Layers (2026-11-28)
- **Status**: Recommended (pending user approval)
- **Rationale**: Subtractive/ONIOM/$\Delta$-ML QM/MLIP (e.g., BAANDI + SevenNet/MACE) requires zero Rust SCF core changes and belongs in the adapter/workflow layer (`baandi-py` / `python/baandi/ase.py`), whereas electrostatic QM/MM and polarizable QM/MLIP reuse the exact same typed `ExternalEnvironment` seam as Track B.

### D-006 — Gate 1 Entry Track Selection (2026-11-28)
- **Status**: **Pending User Approval**
- **Recommended Option (Option 1 — Track B First: Implicit Solvation + Typed Electrostatic Seam → GC-DFT & Track C Reuse)**:
  - **Gate 1 (Autonomous after approval)**: Full 12-item theory document `implicit-solvation-theory.md` (with GC-DFT readiness hooks), exact variational energy/potential/force derivation (including ionic Gaussian pseudo-charge reconciliation and NLCC core-cavity force analysis), SCCS vs. VASPsol/LPB comparison, generalized Poisson PCG vs. direct-minimization algorithm comparison (`algorithm-comparison.md`), and `dev08` external-oracle audit.
  - **Gate 2 (Autonomous after approval)**: Isolated numerical prototype of the distributed FFT-preconditioned generalized Poisson / LPB solver + independent verification suite (Level-1 MMS $< 10^{-12}$, Born/Debye–Hückel analytic limits, exact $\epsilon_b = 1$ vacuum reduction, MPI rank independence) → `feasibility-report.md` with Go/Revise/Hold/Stop decision.
- **Alternative Option (Option 2 — Track C First)**:
  - Broaden `baandi-py`/`ase.py` ($k$-points, spin, FD) + typed fixed-$U_{\mathrm{ext}}$ / GEEP point-charge QM/MM seam + subtractive QM/MLIP prototype (SevenNet/MACE).
- **Alternative Option (Option 3 — Track A First)**:
  - Isolated mixed-basis (PW+AO) or adaptive orthogonalized local enrichment prototype crate outside `refscf` hot paths.
