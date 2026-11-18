# dev08 request: SCF trajectory (start density, first-Pulay bounce, Coulomb guard), 2026-11-18

Ledger rows SCF-START-DENSITY-WIRING, SCF-BOUNCE-ATTRIBUTION, MIX-COULOMB-GUARD, START-DENSITY-ATOMIC-LADDER. MPI against MPI and threads against threads only (comparison-rules.md); no QE column in this request. Load below 0.30, `--map-by core --bind-to core`, one socket (CPUs 0-7), the master's start/done/failed/idle signal files on the shared directory. Report binary sha256, raw stdout, load before and after.

## Why
dev08 (master, 6be8f09, si64 MPI8, 10 fixed maps) saw (1) `start_density = "atomic"` change nothing and (2) a non-monotone standard error (step 3 up 14.8x). Mac findings (counts, O):
- (1) is NOT a dead option. The ladder's coarse stage is the first stage and takes the requested start density (`atomic_truncation_distance` exists only then); the fine stage starts from the prolonged coarse state, and the top-level `atomic_start` block is the FINE session's ("supplied initial vectors retained" by design). results.json now has `effective_start` (per stage) and `convergence.work.cutoff_ladder.coarse_history`. A 10-fixed-map screen cannot show the gain: atomic start shortens the coarse stage (si8/si16/si32 to 4 coarse maps) and the run to convergence.
- (2) reproduces on si16, si32 and si-skew (not si8, si-low, si-k222s, mgo8): fine map 3 is evaluated at the density made by the FIRST Pulay step after the hand-off (one difference, gamma close to -1, step 2x the plain beta step). Same pattern without the ladder (at the first Pulay steps of the cold run) and with a near-exact inner solve (tolerance 5e-14): not inner-solver noise, not the hand-off, not verification. It is the L2 fit over-extrapolating long-wavelength (Hartree-heavy) components.

## Arms (to convergence, production defaults, f1e-3, seeds 601-603, 3 alternating pairs)
Mixer diagnostics and the guard live in the replicated `PulayMixer`; the M2 distributed mixer (`BAANDI_M2_MIX`, default on at np >= 2) has neither. MPI arms therefore pin `BAANDI_M2_MIX=off` on BOTH sides of the guard A/B.
- A: defaults.
- S: `[scf] start_density = "atomic"` (input line; no switch).
- MA: `BAANDI_M2_MIX=off` (MPI only; the guard's baseline).
- G: `BAANDI_M2_MIX=off BAANDI_MIX_COULOMB_GUARD=on` (MPI); `BAANDI_MIX_COULOMB_GUARD=on` (T=1, T=8).
- SG: S plus G.
```
# si64 np=8 per-step table (the original finding), to convergence:
for ARM in A S MA G SG A S MA G SG A S MA G SG; do
  case $ARM in A) E=""; X="";; S) E=""; X=atomic;; MA) E="BAANDI_M2_MIX=off"; X="";;
    G) E="BAANDI_M2_MIX=off BAANDI_MIX_COULOMB_GUARD=on"; X="";; SG) E="BAANDI_M2_MIX=off BAANDI_MIX_COULOMB_GUARD=on"; X=atomic;; esac
  IN=<si64-production.toml>; [ -n "$X" ] && { cp $IN /tmp/in-$ARM.toml; echo 'start_density = "atomic"' >> /tmp/in-$ARM.toml; IN=/tmp/in-$ARM.toml; }   # the line goes in [scf]
  T=traj-si64-np8-$ARM-$RANDOM
  env $E mpirun -np 8 --map-by core --bind-to core target/release/dft-cli run $IN --output-dir ~/fastbench/$T > ~/fastbench/$T.stdout 2>&1
done            # then si32 np=8 and np=4; threads: si32 T=1 and T=8 with arms A S G SG (no BAANDI_M2_MIX)
```
Per run report: wall, coarse maps + fine maps, H columns (coarse, fine, certificate), Davidson iterations, verification attempts, |dE| against arm A, density SHA. Per fine and coarse map (from `history` and `cutoff_ladder.coarse_history`): energy error against the final energy (the like-for-like column with QE), `standard_error_ha`, `density_residual`, `operator_applications`, `pulay.differences`, `pulay.coefficients`, `pulay.coulomb_predicted_ratio`, `pulay.coulomb_refit`.

## Expected (Mac counts, E)
- S: total H -6.8 % (12 runs si8/si16/si32/si-skew), mgo8 -16 %; si32 -2.5 to +2.0 % (one extra verification attempt).
- G: H to reach the 5e-9 Ha band -11.7 % (si16 -10 to -24 %, si32 -9 to -22 %, si-skew -1 to -6 %); full runs -2.0 % at 5e-9 and -7.3 % at threshold 1e-10, because the standard stop estimate under-reads the true error near the stop on some trajectories (extra verification attempts; si-skew +5 to +11 %). Inert (bit-identical) on si8, si-low, si-k222s, mgo8.
- SG: H to the band -16.8 %, full runs -5.6 % (si32 seeds 601/602 +7 to +20 % from 3 extra verification attempts).

## Decision (adoption rule of 2026-11-10)
Each of S (a default change of `start_density`) and G (`BAANDI_MIX_COULOMB_GUARD` default) is adopted when the median wall to convergence improves on si32 and si64 at np 4/8 and T=1/8 and nothing regresses beyond noise; correctness: standard class |dE| <= 5e-9 Ha against A, verification passed. The 10-fixed-map protocol is not valid for either (both change the trajectory length, not the cost per map).
