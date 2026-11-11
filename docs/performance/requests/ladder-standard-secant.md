# dev08 request: standard-profile ladder (S1) and secant seeds across the hand-off (S2), 2026-11-11

Status: implemented, both switches default OFF (bit-identical defaults). Work order 2026-11-11 sections S1 and S2. Per the adoption rule (AGENTS.md, 2026-11-10) there is no percentage threshold: each switch is adopted as the default when the median wall improves on dev08 and no correctness clause below is violated; it stays switchable by its environment variable.

## What the switches do

### S1 `BAANDI_CUTOFF_LADDER_STANDARD=on` (engine `CutoffLadderOptions.standard`)
The C1 ladder (`BAANDI_CUTOFF_LADDER`, default on) was strict-only (`ineligibility` returned "standard profile"). With the switch on, a standard-profile run is eligible:
1. Coarse stage at E1 = 0.5 E2 with the **standard** profile (the same ordinary session, uncertified), stopped by the unchanged rule (density residual 1e-4 Bohr^-3/2 or 30 coarse iterations).
2. Unchanged exact prolongation (orbitals by Miller index, density by Fourier zero padding).
3. Unchanged full-cutoff standard session (same Coulomb-metric convergence test, same 5e-9 Ha threshold, same ratchet), started from the prolonged density and warm orbitals. The only difference to a cold standard session is the eigen tolerance of its FIRST step:
   - `BAANDI_LADDER_STANDARD_SEED=coarse` (default of the switch): the tolerance the coarse standard session would have used next, `ratchet(t_last, E_err_last, N_e) = clamp(min(t_last, 0.1 E_err_last / N_e), 5e-14, 5e-3)`. Reason: a cold session starts at the cap 5e-3 Ha and relies on the first-iteration redo rule (discard the first solve when its error is already below `tolerance * N_e`); after a hand-off the start is already good, so starting at the cap would almost always trigger that redo and pay a discarded solve. The coarse ratchet value is the session's own rule applied to the last measured error; from the first fine record on, the fine ratchet takes over. No constant was added.
   - `BAANDI_LADDER_STANDARD_SEED=cap`: the cold-start rule (cap plus redo), for comparison.
The delivered accuracy class is unchanged by construction: only the full-cutoff session decides convergence, with the unchanged estimator and threshold.

### S2 `BAANDI_LADDER_SECANT=on|<h>` (engine `CutoffLadderOptions.secant_pairs`, on = 4)
`PulayMixer` now keeps an optional block of seeded difference pairs (dx, dr), stored explicitly and closed: no difference is ever formed between a seed and an absolute (input, residual) pair, so the coarse/fine map mismatch (the two fixed points differ at the 1e-4 level) never enters the Gram matrix. API `seed_differences(dx, dr, mean)` (empty history only; `MixerHistory::seed_differences` rejects for other mixers, recorded as `secant_fallback`). `PulayState` carries the seeds (`seed_dx`, `seed_dr`), so checkpoint/restore stays exact; checkpoint v1 serializes them only when non-empty (legacy bytes unchanged). After the coarse stage, the last h coarse difference pairs are prolonged with `prolong_density` (linear, exact for the band-limited coarse fields) and seeded into the fine mixer. The first fine step evaluates the map at the prolonged density; its mix is the seeded least squares with the existing spectral cut (eigenvalues below 1e-12 of the largest dropped). Seeds age out oldest first: seeds + window differences <= history capacity. No mixer parameter changed (P40A-1 not repeated).

## Mac counts (release, pure Rust, T=4, results thread-count independent; observations, not acceptance)
See the ledger rows LADDER-STANDARD and LADDER-SECANT for the numbers; raw rows and medians in `tests/reference/perf/ladder-standard-secant/`. Short form:
- S1 kill rule (median total H in full-grid equivalents lower on si32 and on >= 2 non-silicon cells): passed (si32 2059 -> 1761 at 1e-7 Bohr, 2691 -> 2462 at 0.001 fractional; syn-skew2, syn-fcc, mgo8 at 1e-7 Bohr lower). Regressions in H: si16 at 1e-7 Bohr, si-low, syn-slab, syn-lowsym6 at 1e-7 Bohr.
- S1 accuracy risk: the standard estimator is not a bound. With the coarse seed rule two of three syn-slab seeds at 0.001 fractional ended 9.9e-9 and 1.18e-8 Ha from the strict energy; with the cap rule one si8 seed ended 6.2e-9 Ha. Plain standard max on the same set: 2.1e-9 Ha. Clause 2 below is therefore the decisive check on dev08.
- S2: at most one SCF step fewer (si8, si16, si32 at 0.001 fractional), more steps on mgo8 and syn-fcc; strict energies within 3.1e-11 Ha of the no-ladder run, all certified.
Command:
`LADDER_CASES=... LADDER_SEEDS=601,602,603 LADDER_AMPS=b1e-7,f1e-3 cargo test --release -p dft-engine --test cutoff_ladder ladder_standard_secant_bench -- --ignored --nocapture --test-threads=1`
(`b1e-7` = +-1e-7 Bohr Cartesian uniform per component, `f1e-3` = +-0.001 fractional uniform per component, deterministic xorshift per seed, written in the test.)

## Protocol
One MKL build of the branch HEAD (`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`), the same binary for all arms; arms differ by environment only. Defaults of fa629c1-and-later stay on in every arm (C1 ladder on, EMPTY_LOOSE, BUFFER_TOL 1e-3, B5, ...).
- Cells: `benchmarks/cases/si16-gamma`, `si32-gamma`; profiles standard and strict.
- Coordinates: seeds 601-610 at BOTH amplitudes: (a) +-1e-7 Bohr (the MPI-M1-GATE3 / C1 generator), (b) +-0.001 fractional (`scripts/bench_suite.py --variants` generator with `--displacement 0.001`, or the same per-seed generator as (a) with the fractional amplitude; state which).
- Threads: T=1 (taskset core 0) and T=8 (cores 0-7). Node idle (load below 0.30), uptime before and after; interleave the arms per seed; paired per seed.
- Arms:
  - S1 (standard only): `BAANDI_CUTOFF_LADDER_STANDARD=off` (baseline) vs `=on` (seed rule coarse). Optional third arm `=on BAANDI_LADDER_STANDARD_SEED=cap`.
  - S2 (strict; and standard together with S1=on): `BAANDI_LADDER_SECANT=off` vs `=on` (h=4); optional `=2` and `=6`.
- Report per run: wall, `convergence.iterations`, `convergence.work` (recorded and certificate H vectors), the `cutoff_ladder` block (used, fallback_reason, coarse iterations and H, `coarse_h_full_grid_equivalents`, `standard_seed_ha`, `secant_seeded`, `secant_fallback`), energy (17 digits), density SHA, exit code; binary sha256, commit, build command.

## Correctness clauses (any violation rejects the switch)
1. Strict: every run converges with a passed full-cutoff certificate; |E(arm) - E(no-ladder strict)| below 1e-10 Ha for every seed and T.
2. Standard: every run converges; |E(arm) - E(strict reference of the same coordinates)| below 5e-9 Ha (the standard class), for every seed and T.
3. Thread-count independence: density SHA identical between T=1 and T=8 within each arm.
4. No new non-convergence on any seed. si-elongated / si-slab are known non-convergent and are reported, not counted.

## Adoption
Per switch: adopted as the default when the median wall over seeds 601-610 improves in its cells (si16 and si32, T=1 and T=8, both amplitudes) and no cell regresses beyond the run-to-run noise of the off arm; correctness clauses above hold. No percentage threshold.
