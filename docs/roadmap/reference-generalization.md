# Generalizing the accelerated reference SCF path (user decision A, 2026-11-26)

The reference driver (issue 25) is today Gamma-only, LDA, closed-shell fixed occupations, norm-conserving
pseudopotentials without core charge, no density symmetrization, and not wired to forces, stress or relaxation.
User decisions: 2026-11-25 (generalize when the speed target is near or met, then refactor) and 2026-11-26 (A:
start now, in parallel with the speed work). Each stage is its own ticket; the oracle for every stage is the
legacy path on the same input (already validated against QE where noted) plus the existing analytic/synthetic
tests; the reference path keeps its accuracy classes, verification and rank-count rule.

| stage | capability | oracle | notes |
|---|---|---|---|
| G1 | NLCC (core charge in XC, start density) | legacy path, synthetic NLCC oracles, Si/MgO PseudoDojo smoke | refstart currently rejects core charge |
| G2 | PBE (distributed density gradients) | legacy PBE (validated vs QE on Si) | refstart rejects PBE |
| G3 | k-point meshes (complex path, k distribution) | legacy k-mesh path (validated vs QE) | largest stage; keep the Gamma-real path as a special case |
| G4 | Fermi-Dirac occupations, collinear spin | legacy finite-T and collinear oracles | needs occupation update in the map and Broyden on two spins |
| G5 | forces, stress, relaxation from the reference state | legacy forces/stress tests, finite differences | band-first nonlocal contraction (post-performance-tasks.md) |
| G6 | symmetry (k reduction, density symmetrization) | legacy symmetry evidence | |
| R | refactor: reference path becomes the single SCF path, legacy duplicates retired | full suite, README table | after G1-G6 |

Generality evidence for speed defaults (AGENTS.md "Generality rule"): from G1/G2 on, every speed screen adds
at least one non-Si NLCC/PBE system (e.g. MgO, GaAs, SiC with PseudoDojo files; never committed).
