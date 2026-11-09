# ADR 0042: scf.start_wavefunctions defaults to atomic

Status: Accepted by the pre-registered rule R3 (dev08, f75a817, 2026-11-09); the user asked on 2026-11-09 for the QE comparison to be made at QE's best configuration and for the BAANDI start to be evaluated on equal terms.

## Context
The dev08 comparison with QE 7.6 showed that QE starts from randomised atomic wavefunctions and that its default start is not its best: `startingwfc = 'atomic'` makes QE 1.0 to 1.9x faster on the larger cells and freezes its SCF count. BAANDI started from the lowest free-electron plane waves, needed 13 to 20 SCF iterations and moved with 1e-7 Bohr perturbations (SCF count spread 2 to 7). The opt-in atomic start (P46A: orbitals from the UPF `PP_PSWFC`, Rayleigh-Ritz in the atomic span as Davidson initial vectors) failed an earlier rule (R2) only through a relative accuracy bound that penalised fewer iterations.

## Evidence (R3 stage 2, relayed 2026-11-09, f75a817, both builds, fresh seeds)
8 of 8 pre-registered criteria pass. Perturbed medians (wall, atomic against free electron, T=1): MKL si16 strict -40.7, si16 standard -28.2, si32 strict -31.2, si32 standard -20.1 percent; pure si16 strict -42.4, standard -31.3, si32 strict -32.1, standard -18.6 percent; T=8 si32 strict -29.6 (MKL) and -30.2 (pure), standard -15.0 and -18.6 percent. Unperturbed T=1 MKL: si8 standard -37.5, strict -21.2; si16 standard -25.5, strict -46.4; si32 standard -17.3, strict -51.1 percent (si32 strict 40.64 s to 19.89 s, 13 to 10 SCF iterations). Atomic needs more SCF iterations on 0 percent of the seeds; every run converged; strict runs certified; |dE| vs tight at most 1.65e-10 Ha (standard class 5e-9), strict 5.7e-14 to 3.4e-13; the SCF-iteration spread over perturbed seeds falls to 0 or 1 in every case (7 to 0, 2 to 0, 2 to 1, 1 to 0); density SHA stable run to run and identical T=1 versus T=8.

## Decision
`ScfOptions::default().start_wavefunctions` and the CLI key `scf.start_wavefunctions` default to `atomic`. When the start does not apply (a species without `PP_PSWFC`, collinear spin, Fermi-Dirac occupations, finite-difference paths) the free-electron start is used and the reason is recorded in `results.json` (`atomic_start.fallback_reason`); `free_electron` remains selectable. `results.json` always records `scf.start_wavefunctions`. Old checkpoints without the key are read as `free_electron`, so exact continuation of old runs keeps its start.

## Consequences
- Default SCF counts, trajectories and density hashes change (the start differs); energies agree with the free-electron start within the accuracy classes (not bit-identical).
- Baselines for default-path identity checks must pin `start_wavefunctions = "free_electron"` (and `gamma_real = "off"` for the legacy path).
- The QE comparison at the like-for-like 1-core point (QE best start): si32 strict 11.6x to about 5.7x, median over the six cases 4.0x to 2.9x (master estimate, MKL, unperturbed).
- This closes the part of the gap that was the start; the remaining gap is per-step cost (Davidson algebra) and parallel efficiency.
