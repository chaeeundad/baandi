# C1: cutoff ladder (design note, 2026-11-10)

Status: design only, no code. Source: speed-restructure proposal section 5, C1. Run the experiment only after Track A stage 1 (the all-band minimiser), as the proposal orders.

## Idea
Plane-wave bases are nested: the sphere |k+G|^2/2 <= E1 is a subset of the sphere for E2 > E1, and an FFT grid built for E1 is a sub-grid (in reciprocal space) of the grid for E2. So an orbital or density converged at a reduced cutoff embeds EXACTLY into the full basis by zero padding in G. The ladder runs the first part of the solve (density far from self-consistency, orbitals far from converged) at a reduced cutoff and grid, where every H apply, FFT and GEMM is cheaper by roughly (E1/E2)^{3/2}, prolongs the state, and finishes at the full cutoff. Hand-off and certificate are at the full cutoff only, so the certified result is the result of today's full-cutoff fixed point; the coarse stage only supplies a start.

## Engine entry points
- Basis and grid: `ScfOptions.cutoff_hartree`, `ScfOptions::basis(cell, k, factor)`, `PlaneWaveBasis::new`, `FftGridSpec::new` (density sphere and grid follow the cutoff; a fixed `fft_grid` override must be dropped for the coarse stage).
- Problem: a second `ScfProblem` built with the reduced cutoff (same atoms, pseudopotentials, options otherwise); the coarse SCF runs with the ordinary session (strict tolerances loosened ONLY in the coarse stage, never at the full cutoff).
- Prolongation, new helper in dft-basis: coefficients by Miller index from the coarse to the fine basis (exact injection; coarse coefficients are re-orthonormalised in the fine metric, which is the identity for an injection); density by zero padding in G on the fine grid (Fourier interpolation; the density sphere of the coarse grid is a subset of the fine one).
- Hand-off: `ScfSession::initialize_density_transfer`-like entry with a different basis (new; the existing transfer requires identical model and grid), carrying the prolonged density and the prolonged orbitals as the warm start of the first full-cutoff solve (`step_kpoints(.., initial: Some(..), ..)` already accepts warm matrices). The mixer history starts empty at the full cutoff (histories are not translatable between bases).
- With Track A: the coarse stage is a minimisation at E1 to a loose gradient norm; prolong; continue the minimisation at E2; certificate unchanged.

## Pre-registered experiment (to run after Track A stage 1)
Cells si16, si32 (strict and standard), plus a skewed and a slab-like synthetic cell and one case with NLCC (generality rule); T=1 and T=8 on dev08; fresh seeds.
- Arms: no ladder; ladder with E1 = 0.5 E2; ladder with E1 = 0.7 E2. Coarse-stage stop: density residual 1e-4 (fixed now).
- Identity: certified energy within the class of the no-ladder run (strict: 1e-8 Ha energy, same certificate verdict), eigenvalues within the class, density SHA NOT required (different start).
- Gain: median wall at least 15 percent better on si32 strict at T=1 and T=8 for the better of the two E1 choices, that choice stated before the run as E1 = 0.5 E2 (the 0.7 arm is reported only); no cell worse than +2 percent median.
- Reported: coarse and fine H-vector counts (coarse counted with the grid-size ratio), full-cutoff SCF iterations against the no-ladder run.

## Risks
- Coarse states can order differently (a band that is unoccupied at E1 occupied at E2): the full-cutoff certificate catches it; the gain is lost, the result is not wrong.
- Hard pseudopotentials: at E1 far below the dataset's needs the coarse density is qualitatively wrong and the fine stage needs as many iterations as a cold start; E1 must respect a per-dataset minimum (e.g. the UPF's suggested cutoff), else the ladder is skipped.
- Aliasing: the coarse grid must still hold the coarse density sphere (twice the coarse wave-function radius), otherwise the coarse XC is aliased; GGA gradients and NLCC core densities are evaluated on the coarse grid and change the coarse fixed point (start quality only).
- Fixed FFT grids, symmetry with grid compatibility, Fermi-Dirac, spin and checkpoints: first version only for the integer-occupation spin-unpolarized driver without a fixed grid; everything else falls back to no ladder with a recorded reason.
- Measurement trap: on small cells the coarse stage's setup (second problem, plans) can exceed the saving; the +2 percent bound covers it.
