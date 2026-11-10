# dev08 request: C3 ionic-step start extrapolation (`BAANDI_RELAX_EXTRAPOLATION`), 2026-11-10

Status: implemented on branch `worktree-agent-a719eab3f9d2f2305` (commit in the hand-off report), default `none` (bit-identical to the existing warm transfer). Design: `docs/roadmap/c3-ionic-extrapolation.md`. Per the adoption rule (AGENTS.md, 2026-11-10) there is no percentage threshold: a mode is adopted when the median relaxation wall improves and nothing regresses beyond noise; it stays switchable by the environment variable.

## What the switch does
`[relax] start_extrapolation = "none" | "first_order" | "second_order"`, overridden by `BAANDI_RELAX_EXTRAPOLATION` (same values). Ions mode with `warm_start = true` only (cell mode and cold starts are unchanged). The converged density is split as n = n_at(R) + d, where n_at is the superposition of the UPF atomic charges (`PP_RHOATOM`, the engine's own atomic start density, clipped at zero and normalised to the electron count) at the same positions:
- first_order: start at n_at(R_new) + d_k;
- second_order: d_start = d_k + a (d_k - d_{k-1}) + b (d_{k-1} - d_{k-2}), a and b the least-squares fit of the new displacement R_new - R_k by R_k - R_{k-1} and R_{k-1} - R_{k-2}; with two previous steps only the a term; with nearly parallel previous steps (Gram determinant below 1e-10 |u|^2 |v|^2) the b term is dropped; first order with one previous step.
Negative samples of the predicted start are clipped and the charge restored (the projection the atomic start itself uses); the clipped charge fraction is recorded per frame (`start_clipped_fraction`). Without atomic charge data the plain warm density is used and `start_fallback` says why. The certificate, forces and optimizer are unchanged; the predictor is a start only. The history (up to three d fields and positions) is stored in `restart.json` (optional field, absent when off), so `--resume` reproduces an uninterrupted run exactly (test `relax_extrapolation`).
Per frame, results now carry `scf_h_vectors` (solve plus certificate H vectors) in every mode, and `provenance` names the start kind (`warm_transfer`, `extrapolated_first_order`, `extrapolated_second_order`).

## Mac observations (release, T=1, one run each, not acceptance)
Six BFGS steps from a perturbed start (atom 0 moved by (0.02, 0.01, 0) fractional), force tolerance 1e-9 so all six steps run; default strict relax tolerances. SCF iterations per ionic step (step 0 is the cold start, identical in all modes):

| case | none (warm transfer) | first_order | second_order | cold every step |
| --- | --- | --- | --- | --- |
| si8 Gamma 30 Ry (Si.pz-vbc) | 16, 12, 13, 12, 15, 16, 16 (sum 1-6: 84; 5.98 s) | 16, 11, 11, 12, 13, 14, 15 (76; 5.54 s) | 16, 11, 11, 11, 12, 12, 13 (70; 5.27 s) | 16, 15, 15, 15, 16, 18, 19 (98; 6.21 s) |
| synthetic dimer slab 10 Ry | 19, 12, 14, 15, 28, 22, 26 (117) | 19, 13, 14, 16, 19, 19, 19 (100) | 19, 13, 13, 14, 28, 20, 19 (107) | 19, 23, 28, 26, 32, 36, 34 |
| synthetic fcc 8 Ry | 20, 15, 20, 18, 25, 22, 21 (121) | 20, 14, 15, 19, 19, 23, 22 (112) | 20, 14, 16, 26, 26, 21, 20 (123) | 20, 20, 21, 20, 21, 18, 21 |
| synthetic H dimer 3 Ry (CI test, Gaussian charge not matched to the potential) | 21, 16, 15, 15, 15, 15 | 21, 18, 17, 18, 17, 17 | 21, 18, 16, 16, 15, 17 | n/a |

H vectors follow the iterations (si8 steps 1-6: none 9543, first_order 8668, second_order 8129). Identity: the energy at step 1 (same positions in every mode) agrees to 1.5e-13 Ha (si8) and 3e-14 (synthetic); later steps differ only because the trajectory itself moves by the SCF force tolerance: the existing cold-start relax differs from the warm transfer by the same amount (si8: max |dE| along the path 3.0e-7 Ha cold, 1.4e-7 first_order, 5.3e-8 second_order; max |dF| 6.7e-7, 6.6e-7, 1.27e-6 Ha/Bohr; final max |dx| 5.3e-6, 4.6e-6, 3.6e-6 Bohr). With the tight tolerances of the CI test the modes agree to 6e-13 Ha, 7e-12 Ha/Bohr and 6e-11 Bohr. Lesson: the predictor depends on the quality of `PP_RHOATOM`; with a charge that does not match the potential (the synthetic H test) n_at adds error and the extrapolated start is worse than the plain previous density.

## Protocol
One MKL build of the branch HEAD (`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`); the A/B is the environment only: `BAANDI_RELAX_EXTRAPOLATION=none`, `first_order`, `second_order`. Cases: `benchmarks/cases/si16-k222-relax1` and `si32-gamma-relax1` with `[relax] max_steps = 6` and `force_tol_ha_bohr = 1e-9` (so all six steps run), starting from the case geometry plus perturbed seeds 601-610 (the usual per-atom perturbation of the seed protocol, so every seed has a distinct BFGS path); strict SCF class only (the relax driver rejects the standard profile and any SCF tolerance looser than energy 1e-8 Ha / density 1e-7; run the cases once with the default relax tolerances and once with `energy_tolerance_ha = 1e-10`). T=1 (taskset core 0) and T=8 (cores 0-7). Interleave the three modes per seed; node idle (load below 0.30); uptime before and after.

## Correctness clauses (any violation rejects that mode)
1. Every ionic step's SCF converges with a passed certificate; the relaxation exits with the same status as `none` and the same number of frames.
2. Step 1 (identical positions in every mode): |dE| below 1e-10 Ha (strict class) and forces within 1e-8 Ha/Bohr against `none`.
3. Along the path: |dE| per step within 1e-8 Ha, |dF| within 1e-6 Ha/Bohr and final structure within 1e-5 Bohr of `none`, OR, where the trajectory has moved, no larger than the difference between `none` and a cold-start run (`warm_start = false`) of the same seed (the path's own SCF-tolerance noise; report both).
4. Thread-count independence: per mode, energies and frames bit-identical between T=1 and T=8.

## Adoption
A mode is adopted as the default (rollback `BAANDI_RELAX_EXTRAPOLATION=none`) when the median relaxation wall over the seeds improves on si16-k222 and si32 at T=1 and T=8 and no case regresses beyond the noise of the `none` arm. If both modes pass, the one with the lower median wall is proposed.

## What to send back
Per run: wall, per frame `scf_iterations`, `scf_h_vectors`, `provenance`, `start_clipped_fraction`, `start_fallback`, `energy_ha`, `max_force_ha_bohr`, final positions, exit code; binary sha256, commit, build command, uptime before and after, interleave order.
