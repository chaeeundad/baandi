# P3-ENERGY: dev08 development screen

Master-session request only; no remote action from this worktree. Derived from
`docs/research/qe-study/qe-scf-cycle.md` (own words), hypothesis 7.

Run exactly ONE fully profiled run per arm/case: si64, si128, elongated, rskew,
MPI8, verification on; additionally si128 verification off. Ten runs total.
Use the current lane's identical input/UPF/grid/bands/start/threshold, binary,
8 ranks, one thread per rank and binding in both arms. Record binary commit and
SHA, full command, resolved switches, input and UPF SHA, rank placement.

- A: `BAANDI_REFSCF_LEAN_ENERGY=off`.
- B: `BAANDI_REFSCF_LEAN_ENERGY=on`.
- All other switches identical, including verification gate/reuse/energy mode.
- Full diagnostics on in both; cache reuse is present in both arms.

Collect results JSON, density bytes/hash, stdout and phase/QE-clock profiles.
Report wall, per-map potential/energy/energy_correction, density_transform,
mix, output, verification time, first-stop-to-verified interval; SCF/maps,
Davidson iterations and H columns (including verification), every attempted
verification outcome. Compare delivered complete energy terms (Ha), density
max difference (electron/Bohr^3) and hash, dr2/ethr histories and discrete counts.
Intermediate missing audit terms must be null; mixed variational energies remain
finite. Require unchanged counts/verdicts and delivered differences within
rounding; investigate any discrepancy before adoption. Default remains off.

This is a development screen, not statistical speed evidence. Do not add repeat
runs or compare parallel models. No new QE reference is requested. The supplied
aa6a062 si128 profile (energy 1.08 s, correction 0.25 s, potential 0.29 s across
15 maps, verify off) motivates the work; it is not this candidate's measurement.
