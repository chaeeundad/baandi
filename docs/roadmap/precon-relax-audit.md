# Ionic relaxation audit (PRELAX-PRECON Phase 0, 2026-11-27)

Scope: fixed-cell ionic relaxation (`[relax] mode = "ions"`) as of main
`2efb09be`. Variable-cell (`relax/cell.rs`) is out of scope except where shared.
Line numbers refer to that commit. Mac timings are observations, not acceptance
numbers. Specification: [precon-relax-spec.md](precon-relax-spec.md); design:
[precon-relax-design.md](precon-relax-design.md).

Status 2026-11-27: bottlenecks 1 (stop rule) and 3 (native FIRE units) are
addressed by RELAX-STOP-RULE / ADR 0056 (max atomic force norm, energy test
optional, default 0.02 eV/A, FIRE in ASE units); native frames now record
`max_force_norm_ha_bohr`. Line numbers below still refer to `2efb09be`.

## 1. Code map

| item | where | behaviour |
|---|---|---|
| settings | `crates/dft-cli/src/relax.rs:143-185` | method bfgs/fire, force_tol 1e-3 Ha/Bohr, energy_tol 1e-8 Ha, max_steps 200, max_displacement 0.2 Bohr, trust_radius 0.2 Bohr, warm_start true, start_extrapolation none, fixed list |
| optimizer state | `crates/dft-cli/src/relax/optimizer.rs:5-35` | dense inverse Hessian `h` (3N x 3N), FIRE velocity/dt/alpha, previous x and force, trust radius |
| BFGS | `optimizer.rs:63-92` | see 2.1 |
| FIRE | `optimizer.rs:40-62` | see 2.2 |
| step cap | `optimizer.rs:94-100` | uniform scaling so that the largest per-atom Euclidean displacement <= min(trust radius, max_displacement) |
| constraint projection | `relax.rs:240-248`, applied to forces and step at `relax.rs:1092-1100` | Cartesian components zeroed; raw physical forces kept in frames |
| ionic loop | `relax.rs:1083-1188` | propose, try once, retry once fresh at half displacement, accept, update radius/stop counters, save restart + outputs |
| evaluation | `relax.rs:398-783` (`evaluate_with`) | reference SCF through `refscf::api`, verified forces (+ D3 if configured), projected max component, frame record, warm seed for the next step |
| stop rule | `relax.rs:1158-1179` | see 2.4 |
| restart | `relax.rs:284-298` (payload), `835-878` (save/load) | see 2.7 |
| outputs | `relax.rs:1201-1400` | trajectory.extxyz, final.xyz/.cif/POSCAR/final.toml, results.json, relax.log, one stdout line per step |

## 2. Findings

### 2.1 BFGS (`optimizer.rs:63-92`)
- Dense inverse-Hessian BFGS. Initial H^-1 = I Bohr^2/Ha, i.e. isotropic springs of
  1 Ha/Bohr^2 = 97.2 eV/A^2 on every Cartesian coordinate (ASE BFGS default
  alpha = 70 eV/A^2): no element, bond or connectivity information.
- Update: standard inverse BFGS with s = x_k - x_{k-1}, y = -(F_k - F_{k-1});
  skipped when s.y <= 1e-12 |s||y| (or <= 1e-30). The pair is formed between
  consecutive ACCEPTED frames (a half-step retry just changes s).
- Step p = H^-1 F; a non-finite or uphill (p.F < 0) direction RESETS the whole
  Hessian to the spring model and uses p = F.
- No line search and no step rejection: a step that raises the energy is
  accepted; the only reaction is halving the trust radius (`relax.rs:1166-1170`),
  which otherwise grows by 1.1 per step up to `trust_radius_bohr`. This is a
  capped BFGS step, not a trust-region subproblem.
- Memory and work are O((3N)^2) per step (72 MB at 1000 atoms, ~7 GB at 10^4):
  a scalability ceiling independent of the DFT cost.

### 2.2 FIRE (`optimizer.rs:40-62`)
- Unit fictitious masses, dt0 = 0.1, dt_max = 1, alpha0 = 0.1, N_min = 5, f_inc 1.1,
  f_dec 0.5, f_alpha 0.99, all in Hartree atomic units. The displacement per step
  is dt^2 F (Bohr per Ha/Bohr). ASE FIRE uses the same dimensionless numbers in
  eV/A and A, i.e. a displacement dt^2 F(eV/A) in A. For the same physical force
  and dt, the native step is smaller by the factor
  (51.42 eV/A per Ha/Bohr) x (1.8897 Bohr per A) = 97.2
  (equivalently the native FIRE is ASE FIRE with masses 97 times heavier). With
  forces of 1e-2 Ha/Bohr and dt = 0.1 the first steps are 1e-4 Bohr
  (observed: `exp-tabpfn/relax/runs/si8-a0.05-s1-fire/out/relax.log`, step
  lengths 4e-5 .. 1.5e-3 Bohr over the first 14 steps). dt only grows after five
  positive-power steps and saturates at 1, so the native FIRE is very
  conservative by construction (units, not physics).
- Nonpositive power: velocity zeroed, dt halved, alpha reset. Retries also zero
  the velocity (`relax.rs:1171-1173`).

### 2.3 Warm starts, extrapolation, SCF per step
- Every accepted structure starts a new reference SCF (`refscf::api`). With
  `warm_start = true` the `IonicSeed` of the previous successful evaluation
  (density and compatible orbitals; `relax.rs:422-445`, stored at `762-780`) is
  used; potentials and mixer history are rebuilt. `warm_start = false` uses fresh
  atomic starts. A half-step retry always starts fresh (`relax.rs:1128`).
- Optional C3 extrapolation (`start_extrapolation`, env override
  `BAANDI_RELAX_EXTRAPOLATION`, `relax.rs:38-141`): first/second order on the
  bonding density n_scf - n_at; default none.
- SCF per step: the single standard path with fixed-point verification on every
  step (ONEPATH-1); ions mode default threshold 5e-9 Ha; looser energy/eigen/
  density tolerances reject at preflight (`relax.rs:961-967`); `verify = "off"`
  rejects. Forces only from a verified state (`relax.rs:603-610`).
- Symmetry auto: density and raw forces symmetrised; constraints must commute
  with the group (`relax.rs:905-930`).

### 2.4 Stop rule (`relax.rs:1158-1179`)
- Force measure: **maximum absolute free Cartesian component** of the projected
  force (`relax.rs:742`), not the atomic norm. Max component <= norm <= sqrt(3)
  max component, so the native 1e-3 Ha/Bohr component test corresponds to a norm
  between 1e-3 and 1.73e-3 Ha/Bohr (0.051..0.089 eV/A), 2.6..4.5 times looser
  than the 0.02 eV/A norm target of the benchmark.
- Energy: |E_k - E_{k-1}| < energy_tol (default 1e-8 Ha) must hold together
  with the force test on **two consecutive accepted steps**. Consequences: at
  least one evaluation beyond the first one that meets the force test; and the
  energy test, with a default only 2x the SCF threshold (5e-9 Ha), keeps the
  relaxation running long after the forces are converged. Example
  (`exp-tabpfn/relax/runs/si8-a0.05-s1-bfgs`, Si8 LDA k222): component < 1e-3
  at step 4, energy test met at steps 15-16, i.e. 12 extra evaluations (75 % of
  the run).
- No stop rule uses the atomic force norm; the benchmark therefore uses an
  external evaluator (design Section 5).

### 2.5 Fixed atoms, PBC
- Fixed atoms / components (`[relax] fixed`) zero the projected force and the
  step; unchanged fractional coordinates of frozen atoms are copied bit-exactly
  (`relax.rs:1110-1121`).
- Positions are unwrapped Cartesian; fractional coordinates are recomputed
  without wrapping; the cap is per atom in Cartesian Bohr. The optimizer itself
  does not need minimum-image logic (it works on displacements), but a neighbour
  graph for a preconditioner will (design 6.2).

### 2.6 Failure policy (`relax.rs:1102-1182`)
- A failed or unconverged candidate SCF is retried once, fresh, at half the
  proposed displacement; if that fails the run ends `failed` (exit 4) with the
  last good structure. Initial SCF failure: one fresh retry at the same
  structure. Failures are recorded with coordinates and reason; each failed
  attempt is a full SCF cost.

### 2.7 Restart
- `restart.json` after every accepted step (atomic publication, f64 as IEEE
  bits, checksum): optimizer state (dense H^-1, FIRE state), frames, last
  density, stop counter, failures, extrapolation history. The orbital seed
  (`IonicSeed`) is NOT persisted: the first SCF after `--resume` is a fresh start.

### 2.8 Records
- relax.log: step, energy (Ha), projected max component (Ha/Bohr), step length
  (Bohr), SCF iterations, provenance (fresh/warm_transfer/extrapolated).
- results.json frames: positions, raw forces, energy, projected max component,
  step length, SCF iterations, `scf_h_vectors` (solves + certificates),
  verification summary, start fallback/clipping; failure records. **No wall time**
  per step (deliberately: frames are part of the bit-identical restart), no count
  of attempted evaluations as such (it is frames + failures), no atomic force norm.

## 3. Mac baseline (observations)

Executable: `target/release/baandi` built from `2efb09be` on the PRELAX-P0
worktree (pure build, no MKL/MPI). Raw results:
`tests/reference/relax-precon/results/mac-baseline-2026-11-27/`.

Native optimizers as configured (own stop rule: max component < 1e-3 Ha/Bohr and
|dE| < energy_tol on two consecutive steps; energy_tol 1e-8 Ha for Si8, 1e-7 Ha
for Si8-AlP), warm starts, one process per run. "target at" is the first
evaluation whose max atomic force norm is below 0.02 eV/A (3.889e-4 Ha/Bohr),
recorded by the harness without stopping the run.

| case (threads) | seed | method | evals at native stop | 0.1 / 0.05 / 0.02 eV/A reached at | first max component < 1e-3 | SCF iterations | wall (s) |
|---|---|---|---|---|---|---|---|
| Si8 LDA k222, 0.2 Bohr (1) | 1 | BFGS | 22 | 8 / 9 / 18 | 8 | 133 | 41 |
| Si8 LDA k222, 0.2 Bohr (1) | 2 | BFGS | 28 | 13 / 14 / 16 | 14 | 174 | 65 |
| Si8 LDA k222, 0.2 Bohr (1) | 1 | FIRE | 77 | 46 / 53 / 67 | 51 | 432 | 224 |
| Si8 LDA k222, 0.2 Bohr (1) | 2 | FIRE | 72 | 48 / 54 / 59 | 54 | 375 | 177 |
| Si8-AlP PBE k222 (4) | 1 | BFGS | 22 | 6 / 15 / 16 | 14 | 160 | 77 |
| Si8-AlP PBE k222 (4) | 1 | FIRE | 62 | 40 / 52 / 58 | 52 | 394 | 346 |

No SCF failure, no retry and no energy rise in any of these runs.

Earlier run kept for context (exp-tabpfn, older executable sha256 d49da41d...,
Mac, 4 threads; observation only): 63-atom Si59Al3P1 + vacancy
(`relax/runs/vac2-s1-native`, PBE 48 Ry Gamma), native BFGS warm: 0.1 / 0.05 /
0.02 eV/A reached at 21 / 36 / 38 evaluations, native stop at 61 (energy test
1e-7 Ha). The energy above the final value decays slowly over steps 20-60
(1.6e-3, 8e-4, 4e-4, 2e-4 Ha every five steps): a soft collective relaxation
around the defect. ASE BFGS (cold, max-component stop at 1e-3) stopped after 21
calls 5.1e-4 Ha (14 meV) ABOVE the native final energy: the max-component
1e-3 stop is premature for this defect.

Early/late stages and SCF noise (Si8 seed 1, native BFGS, per-step trace in
`results/mac-smoke-2026-11-27/si8-lda-k222-a0.2-s1-baandi-bfgs-warm.json`):
- early (evaluations 1-6): large capped steps (0.07-0.15 Bohr), energy falls
  37 mHa, force norm 6.8e-2 -> 2.5e-3 Ha/Bohr;
- late (7-18): steps of 3e-3 to 1e-2 Bohr, far below the 0.2 Bohr cap, force
  norm stalls near 1e-3 Ha/Bohr for eight evaluations while the energy still
  falls by about 1e-5 Ha per step. The per-step energy changes are three orders
  of magnitude above the SCF threshold (5e-9 Ha) and the forces are verified, so
  this tail is **ill-conditioning (soft collective mode not represented by the
  isotropic 97 eV/A^2 start Hessian), not SCF noise**. A3 (7 evaluations) and
  PreconLBFGS (9 including its probe) do not show the stall;
- SCF noise would matter only for the native energy stop test (1e-8 Ha, 2x
  the SCF threshold), not for the 0.02 eV/A force target: force differences
  attributable to the SCF threshold were measured at about 1e-6 Ha/Bohr
  (ONEPATH-1 note in docs/validation/relaxation.md), 400 times below the target.
- Warm starts: late steps take 4-5 SCF iterations warm against 7-8 cold;
  over the whole BFGS run 120.5 versus 136 iterations (median of two seeds,
  -11 %), FIRE 366 versus 526 (-30 %), because FIRE's many tiny steps profit most.

## 4. Bottlenecks (ranked by expected effect on evaluations)

1. **Native stop rule** (2.4): with the default energy tolerances the native
   run continues 4 to 12 evaluations (+22 % to +75 %) past the point where the
   0.02 eV/A norm target is met (Si8 BFGS 22/18 and 28/16, AlP 22/16; 63-atom
   vacancy 61/38), and its max-component force test is looser than the norm
   target. This is the largest and cheapest lever, but it is a stop-rule
   question (user-visible semantics), separate from preconditioning; it must be
   decided explicitly, not changed silently.
2. **Ill-conditioned late stage** (BFGS start Hessian = isotropic 97 eV/A^2):
   long tails of small steps with nearly constant force (Si8 evals 7-18; the
   63-atom vacancy steps 20-60). This is what a physical preconditioner targets.
   Smoke evidence (Si8, 2 seeds, evaluations to 0.02 eV/A, median): native BFGS
   17, ASE BFGS/LBFGS 15, ASE PreconLBFGS(Exp) 8.5 (including one probe), A3 7.
3. **Native FIRE is ~97x too conservative by units** (2.2): 59-67 evaluations
   on Si8 against 38-45 for ASE FIRE and 16-18 for BFGS; FIRE is not a
   competitive baseline as implemented. A unit fix (masses or dt in eV/A-like
   scaling) is a separate small change, out of this task's scope.
4. **Conservative steps**: not limiting in these cases; the cap (0.2-0.3 Bohr)
   binds only in the first 2-3 steps; the trust radius never halved (no energy
   rises).
5. **Unnecessary evaluations**: none in the optimizer itself (every native
   evaluation is an accepted step; no line search, no probes); the waste is in
   the stop rule (1) and in the conditioning (2).
6. **SCF noise**: not limiting for a 0.02 eV/A force target at the default
   5e-9 Ha threshold; it is close to limiting for the native energy test
   (energy_tol 1e-8 Ha = 2x the threshold).
7. **Per-evaluation cost**: warm starts already cut SCF iterations (-11 % BFGS,
   -30 % FIRE). The restart does not persist the orbital seed (first SCF after
   resume is cold). Dense O(N^2) BFGS memory is a scalability limit beyond
   ~10^3 atoms.

## 5. Phase 1 harness smoke results (Mac observations, 2026-11-27)

Harness: `tests/reference/relax-precon/` (design Sections 5, 8). Common stop:
first candidate with max atomic force norm < 0.02 eV/A; every attempted
evaluation counted. Medians over seeds; raw JSON and `summary.json` in
`results/mac-smoke-2026-11-27/` (Si8, 1 thread) and `results/mac-alp-2026-11-27/`
(Si8-AlP, 4 threads). No SCF failure in any run.

| arm | Si8 LDA k222 0.2 Bohr, 2 seeds: evals to 0.1 / 0.05 / 0.02 eV/A | Si8-AlP PBE k222, 3 starts: evals to 0.1 / 0.05 / 0.02 eV/A |
|---|---|---|
| baandi-bfgs:warm | 10.5 / 11.5 / 17 | 6 / 15 / 16 (per start 16, 16, 35) |
| baandi-bfgs:cold | 10.5 / 11.5 / 17 | 6 / 15 / 16 |
| baandi-fire:warm | 47 / 53.5 / 63 | (not run in the matrix; own-stop baseline: 58) |
| ase-bfgs:cold | 9.5 / 11.5 / 15 | 5 / 13 / 15 |
| ase-lbfgs:cold | 9.5 / 11.5 / 15 | 5 / 13 / 15 (per start 14, 15, 30) |
| ase-fire:cold | 20.5 / 27 / 41.5 | 11 / 23 / 35 |
| ase-precon-lbfgs:cold | 5.5 / 7 / 8.5 (incl. 1 probe) | 5 / 8 / 9 (per start 9, 8, 14; incl. 1 probe, 2 rejected trials in total) |
| a3-lbfgs:cold | 4.5 / 5.5 / 7 | 3 / 5 / 6 (per start 6, 6, 9) |

Paired ratio of evaluations to 0.02 eV/A against native BFGS (same seed; all
pairs won): A3 0.41 (Si8) and 0.38 (AlP); PreconLBFGS 0.50 and 0.50; ASE
BFGS/LBFGS 0.88 and 0.88; ASE FIRE 2.44 and 1.62; native FIRE 3.70 (Si8).
Cold-versus-cold wall ratios against native BFGS cold: A3 0.43 (Si8) / 0.38
(AlP); PreconLBFGS 0.52 / 1.10. The Mac walls are noisy (per-evaluation wall
of identical 11-iteration SCFs varied 7-24 s during the AlP runs, and the warm
native arm was slower than the cold one there), so these walls are not
evidence; the evaluation counts are.

Same-minimum check: all Si8 arms agree (|dE| < 1e-4 Ha, final structures within
0.05 Bohr after removing the rigid translation). Si8-AlP start 3: A3 and ASE
FIRE end 0.15 Bohr away from native BFGS with |dE| = 5e-6 and 1.7e-6 Ha:
a very soft mode (all three are below the force target; the energy differences
are 3 orders below a meaningful structural difference), flagged and not counted
as a different minimum. Start 3 is also the case with the long native/LBFGS
tail (35 and 30 evaluations) where the preconditioned methods need 9-14.

Reading: on these small bulk-like cells the physics-based metric halves
(PreconLBFGS, including its probe) or more than halves (A3) the evaluations
to 0.02 eV/A; the gain is concentrated in the late stage (0.05 -> 0.02 eV/A).
This is an 8-atom observation with 2-3 seeds on the Mac; it motivates, but does
not replace, the dev08 benchmark with the 63-atom defect cases.

## 6. Reproduction

```sh
PY=/Volumes/990-2TB/Develop/rust-dft-wt/venv-convert/bin/python   # Python 3.9.6, ASE 3.26.0, NumPy 2.0.2, SciPy 1.13.1
cd tests/reference/relax-precon
$PY -I run_bench.py --case cases/si8_AlP-pbe-k222.json --seed 1 --method baandi-bfgs --start warm \
   --native-own-stop --exe ../../../target/release/baandi \
   --pseudo-root /Volumes/990-2TB/Develop/rust-dft-wt/pseudo --threads 4 --max-evals 120 \
   --workdir /scratch/rb/alp1-bfgs --out results/mac-baseline-2026-11-27/si8_AlP-pbe-k222-s1-baandi-bfgs-warm-own.json
python3 -I summarize.py results/mac-baseline-2026-11-27
```
