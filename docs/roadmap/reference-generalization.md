# Generalizing the accelerated reference SCF path (user decision A, 2026-11-26)

Status (Stage R, 2026-11-25): the reference driver is the only production path.
NLCC, PBE, complex k points, fixed spin/eligible FD, symmetry, forces/stress,
relaxation, NSCF and runtime/checkpoint v2 are implemented. Free-moment PBE (G4D)
remains unqualified and routes to the strict oracle. Current contracts and scoped
evidence: [ADR 0047](../decisions/0047-reference-driver-single-path.md).
The G-stage table and early measurements below are the historical rollout record;
their original rejection notes do not describe current refstart support.

| stage | capability | oracle | notes |
|---|---|---|---|
| G1 | NLCC (core charge in XC, start density) | legacy path, synthetic NLCC oracles, Si/MgO PseudoDojo smoke | refstart currently rejects core charge |
| G2 | PBE (distributed density gradients) | legacy PBE (validated vs QE on Si) | refstart rejects PBE |
| G3 | k-point meshes (complex path, k distribution) | legacy k-mesh path (validated vs QE) | largest stage; keep the Gamma-real path as a special case |
| G4 | Fermi-Dirac occupations, collinear spin | legacy finite-T and collinear oracles | needs occupation update in the map and Broyden on two spins |
| G5 | forces, stress, relaxation from the reference state | legacy forces/stress tests, finite differences | band-first nonlocal contraction (post-performance-tasks.md) |
| G6 | symmetry (k reduction, density symmetrization) | legacy symmetry evidence | |
| R | Done except G4D: single production reference path, oracle isolated | [ADR 0047 evidence](../decisions/0047-reference-driver-single-path.md#evidence-gen14gen28), [R plan](refactor-r-plan.md#8-evidence-log) | R11b implementation and gen22–gen28 evidence require integration reconciliation; README corrections await D6 approval |

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
| G4 | Si2 Fermi-Dirac theta 0.01 Ha, 4x4x4 (36 k, 32 bands) | -2.7e-10 (F) | +2.2e-11 (F) | mu_ref - mu_oracle -1.5e-7 Ha (legacy -3.4e-6); verified |

All within the standard class (5e-9 Ha). Evidence: dev08 /data/sig-dft/gen/RESULT.md /data/sig-dft/gen4/RESULT.md and /data/sig-dft/gen5/RESULT.md (binary e3a75a5, sha256 183df05a...5a49).

Historical gaps found by the early runs (G4C is now implemented; current open items are in ADR 0047): spin with PBE needs a two-channel distributed GGA evaluator (stage G4C; real PseudoDojo files are PBE, so this blocks spin on real materials in the reference driver); `refscf-check --stop-to-verified` does not read the k/FD results schema (fold into G5B). Evidence for G4: dev08 /data/sig-dft/gen6/RESULT.md.

### bcc Fe (dev08 gen15, 2026-11-25)

PseudoDojo Fe (PBE, NLCC), 2-atom cubic cell a = 5.42 bohr, ecut 80 Ry, 4x4x4 (64 k per spin, time reversal off), Fermi-Dirac 0.01 Ha, fixed total moment 4.4 mu_B, MPI8. Reference driver with the G4B switches on: VERIFIED in 75.8 s, 18 maps, F = -250.682665063 Ha, M = 4.4 mu_B, RSS 578 MB. The legacy driver and the strict legacy oracle each failed to finish within 45 min (2700 s timeouts; legacy prints no progress), so no oracle agreement exists for this case. Evidence: dev08 /data/sig-dft/gen15/RESULT.md.
