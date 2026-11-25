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

## Evidence (dev08, MPI8, one run per arm; legacy pulay vs reference broyden vs strict legacy oracle)

| Stage | Case | E_ref - E_legacy (Ha) | E_ref - E_oracle (Ha) | Notes |
|---|---|---|---|---|
| G1 | Mg atom, LDA-PZ NLCC (Mg.pz-n-vbc) | +3.6e-10 | +4.4e-10 | both verified; core electrons identical (gen 2432530) |
| G2 | si32 PBE (Si.pbe-rrkj, old label accepted at e3a75a5) | -1.2e-10 | +3.3e-11 | rerun without override differs by -2.8e-14 |
| G2B | MgO 8 atoms, PseudoDojo PBE, O with NLCC | +6.1e-10 | +7.3e-10 | occupied eigenvalues 3.8e-6 Ha |
| G3 | Si2 4x4x4 unshifted (36 k) | +1.3e-10 | +1.9e-10 | per-k occupied eigenvalues <= 3.7e-6 Ha |
| G3 | Si2 2x2x2 shifted (4 k) | +1.3e-10 | +2.1e-10 | E_ref - E_QE(k222s_c) = +2.56e-6 Ha, the archived cross-code offset of this input |
| G5 | si8 Gamma, one atom +0.05 bohr, forces+stress, default standard | +2.97e-10 (vs oracle) | | max dF vs oracle 8.2e-6 Ha/bohr (legacy 4.9e-6), dsigma 1.1e-9; labels "verified" |
| G5 | same, standard threshold 5e-13 Ha | | +2.0e-13 | dF 2.2e-7 Ha/bohr, dsigma 4.2e-11; 3-step relax vs legacy: dE <= 4.1e-9 Ha, dF <= 3.1e-7, dr <= 9.3e-7 bohr (gates 5e-9 / 1e-6 / 1e-6 pass) |

All within the standard class (5e-9 Ha). Evidence: dev08 /data/sig-dft/gen/RESULT.md /data/sig-dft/gen4/RESULT.md and /data/sig-dft/gen5/RESULT.md (binary e3a75a5, sha256 183df05a...5a49).
