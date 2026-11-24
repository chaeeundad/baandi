# Tasks after the performance work (user, 2026-11-20)

The performance campaign (BAANDI <= 0.90 x fastest QE, one-core and MPI lanes) comes first. Items below start only after it is declared complete by the user.

1. **Apply TabPFN to the structure-optimization process** (added 2026-11-20, user request). Scope idea only, no design yet: use a tabular foundation model (TabPFN) inside ionic/cell relaxation, for example to predict the next step or the Hessian/curvature correction from the history of positions, forces and stresses, and to decide when to stop. Needs, before any code: a design contract approved by the user (what is predicted, what stays exact), an independent validation oracle (the existing finite-difference and energy-scan checks of ionic and cell relaxation, `docs/validation/relaxation.md` and `cell-relaxation.md`), a license/maintenance/Apple Silicon/Linux check for TabPFN and its weights (new dependency rules, ADR 0003), and a rule that every relaxation result remains certified by the engine's own forces/stress, so the model can only propose steps, never replace the physics. Ticket and ADR to be written at kickoff.

2. **Generalize the accelerated reference SCF path, then refactor** (user decision 2026-11-25). Trigger: the
   speed target is near or met (si128 MPI8 verify-off was 1.19x the fastest QE on 2026-11-25). The reference
   driver (issue 25) is today Gamma-only, LDA-only, closed-shell fixed occupations, NC pseudopotentials without
   core charge, no density symmetrization, and not wired to forces/stress/relaxation. Extend it, each through its
   own ticket and the existing validation oracles: k-point meshes (complex path), PBE (distributed gradients),
   NLCC/PseudoDojo files, collinear spin and Fermi-Dirac occupations, symmetry, forces/stress/relaxation. Then a
   refactor that makes it the single SCF path and retires duplicated legacy code, with the README feature table
   kept true at every step.
