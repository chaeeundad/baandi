# ADR 0058: Opt-in preconditioned L-BFGS for fixed-cell ionic relaxation

- Status: **Proposed** (ticket PRELAX-P2, 2026-11-28; orchestrator review
  pending; dev08 rerun of stages 1-3 and the generality cases not yet made).
- Scope: `crates/dft-cli/src/relax.rs` (method value, records, restart fields,
  dispatch), `relax/precon.rs` (model Hessian, sparse solve, L-BFGS memory,
  line search), `relax/precon_driver.rs` (ionic loop). No change to SCF,
  forces, BFGS or FIRE; the default optimizer stays BFGS.
- Relates to: ADR 0056 (stop rule on the atomic force norm), ADR 0057
  (force-tied SCF threshold, used as the Armijo noise allowance),
  [design note](../roadmap/precon-relax-design.md) (Phase 0/1, Section 10 for
  the dev08 rerun), [audit](../roadmap/precon-relax-audit.md).

## Context

Phase 1 (dev08, 2026-11-27/28; one core; external stop at max atomic force norm
0.02 eV/A, every attempted evaluation counted) measured the experimental A3
optimizer (L-BFGS whose H0 is the inverse of a fixed bond-spring model Hessian,
Armijo line search) through the BAANDI `run` calculator: median evaluations Si8
LDA 0.05 Bohr 5, Si8 0.2 Bohr 7, Si8-AlP 6, the three 63-atom Si-Al-P vacancy
cases 7/9/9, i.e. 0.21-0.41 times native BFGS (ASE PreconLBFGS(Exp) 0.39-0.50,
ASE LBFGS about 0.87), with no extra failures and no different minima. That
meets the pre-registered rule (at least 20 percent fewer than native BFGS and
ASE LBFGS). A3 had two properties that a production optimizer must not keep: a
fixed absolute stiffness (5 eV/A^2) and a dense inverse.

## Decision

`[relax] method = "precon_lbfgs"` (ions mode; `cell`/`ions_cell` reject at
preflight; BFGS remains the default; no new key, no environment variable).

1. **Model Hessian (shape).** P1 = sum over neighbour pairs (i, j, lattice image)
   of w B (e_i - e_j)(e_i - e_j)^T + c I on the free Cartesian coordinates, with
   B = (1 - s) I + s n n^T, s = 0.7, w = exp(-A (d/r0 - 1)), A = 2.5,
   r0 = r_cov,i + r_cov,j (Cordero et al., Dalton Trans. 2008, 2832; Z <= 96,
   otherwise an input error), pairs closer than 1.3 r0 (1.8 r0 between a
   molecular element H/C/N/O and another element), c = 0.01. These are the A3
   values that Phase 1 measured. Reasons: a bond stretch dominates the Hessian
   of a covalent or ionic bond (for a central pair force the transverse part is
   phi'/r, near zero at equilibrium); the transverse 0.3 represents angular
   stiffness that a pair graph cannot express; the exponential makes
   compressed bonds stiffer and stretched contacts softer (Packwood et al.,
   J. Chem. Phys. 144, 164109 (2016), use the same form with A = 3 and
   isotropic blocks, which Phase 1 measured worse: 0.39-0.50 against
   0.21-0.41); the stabiliser pins the translation null space, isolated atoms
   and molecules in vacuum (c relative to a unit bond weight; Packwood use 0.1,
   A3 0.01). Fixed components are eliminated (rows and columns removed), so P1
   is SPD on the free subspace: v.P1 v >= c |v|^2 (tested).
2. **Scale mu without a fitted constant or extra evaluation.** P = mu P1.
   Before the first secant pair mu = 1 Ha/Bohr^2 / mean free diagonal of P1:
   the BFGS starting spring (H0 = I Bohr^2/Ha, the existing prior), distributed
   by the model. After every accepted step with positive curvature
   (s.y > 1e-12 |s||y|) mu = y.P1^-1 y / s.y from that pair, the P-metric form of
   the standard L-BFGS scaling gamma = s.y / y.y (Nocedal and Wright eq. 7.20).
   Alternatives measured on the Mac smoke (table below): s.y / s.P1 s (the
   Rayleigh quotient of design 6.3(a)), mu from the first pair only, and a first
   step scaled to the displacement cap instead of the prior.
3. **L-BFGS**: memory 10 (Phase 1 value), two-loop recursion with
   H0 = P^-1; a non-descent or non-finite direction resets the history and uses
   -P^-1 g. P1 solves: block-Jacobi preconditioned conjugate gradients to a
   relative residual 1e-13 (no new dependency: `faer` is not a dependency of
   `dft-cli`, and the solve costs O(bonds) per iteration, negligible against an
   SCF). Neighbour search: cell lists on wrapped fractional coordinates with
   bins at least one cutoff wide along each plane normal (thin cells scan the
   needed images), O(N).
4. **Rebuild** P1 at the current positions when an atom moved more than 0.1 of
   the shortest bond since the last build (weights then change by at most
   about exp(0.25) = 1.28); the L-BFGS pairs are kept (they are secant data of
   E; the two-loop recursion accepts a different H0 at every step). Measured
   on a Morse fcc vacancy model (31 atoms, +-0.4 Bohr starts, 3 seeds,
   evaluations to 1e-5 Ha/Bohr; `relax::precon::tests`): rebuild fraction 0.5
   (Packwood's rule) 34/18/11, 0.1 22/17/10, every step with pairs kept
   20/18/10, every step with pairs reset 110/26/11; native BFGS 74/64/44.
   No rebuild happens in the Phase 1 DFT cases (moves of 0.2-0.3 Bohr against
   0.44 Bohr thresholds).
5. **Line search** per ionic step along p (capped so that no atom moves more
   than min(`max_displacement_bohr`, `trust_radius_bohr`); no trust-radius
   adaptation): accept at t when E(t) <= E0 + 1e-4 t g.p + eps,
   eps = 2 x the SCF energy accuracy (the standard threshold, ADR 0057 for
   relaxations; 9.5e-11 Ha at the default tolerance), or when the trial already
   passes the force-only stop rule (it is a converged structure on verified
   forces). Otherwise t <- the zero of the linearly interpolated directional
   derivative g(t).p (exact on a quadratic; we have the full gradient at every
   trial), or the energy parabola if the curvature is not positive, within
   [0.1 t, 0.5 t]. Three trials per direction; then one round of three on the
   preconditioned steepest-descent direction with a history reset; then status
   `line_search_failed`, exit 3, last good structure kept. Failed SCF: one retry
   fresh at half the step (counted), a second consecutive failure ends the run
   with exit 4 (BFGS semantics).
6. **Warm starts** as BFGS: the next trial starts from the seed of the last
   successful evaluation (reference driver; accepted or rejected), the strict
   oracle path from the last accepted density.
7. **Records and restart.** Every attempted evaluation gets a global
   `evaluation` index (frames `optimizer.evaluation`, failures, rejected
   trials). Frames gain an `optimizer` record (mu, mu source, t, trials,
   rejections, failures, fallback, pair stored, history length, rebuild and
   graph size); `relaxation.rejected_trials` holds every evaluated but rejected
   trial (positions, energy, forces, SCF iterations and H vectors);
   `relaxation.optimizer` the summary and the model constants. `relax.log`
   appends five columns; a rejected trial republishes the outputs at once and
   prints `trial step ...`. restart.json gains `precon` (mu, pairs, build
   positions, counters; IEEE bits like the rest) and `rejected`; the dense
   BFGS/FIRE state is empty for this method. P1 is rebuilt from the build
   positions on resume (deterministic; round trip tested bit for bit).

## Evidence in this change (Mac, observations)

Mac (Apple Silicon), release build of this change, 4 threads, harness
`tests/reference/relax-precon` with `--scf-threshold-ha force-tied`
(4.73e-11 Ha for every arm), external stop at max atomic force norm
0.02 eV/A, every attempted evaluation counted. Raw JSON and summaries:
`tests/reference/relax-precon/results/mac-p2-smoke-2026-11-28/`.

| case | start | native precon (warm) | native BFGS (warm) | Python A3 (cold) | final E precon / BFGS / A3 (Ha) |
|---|---|---|---|---|---|
| Si8 LDA k222, 0.2 Bohr | seed 1 | 8 | 17 | 7 | -31.6697197 / -31.6697201 / -31.6697191 |
| Si8 LDA k222, 0.2 Bohr | seed 2 | 9 (1 rejected trial) | 16 | 7 | -31.6697209 / -31.6697201 / -31.6697209 |
| Si8 LDA k222, 0.2 Bohr | seed 3 | 8 | 16 | 7 | -31.6697205 / -31.6697197 / -31.6697207 |
| Si8-AlP PBE k222 | start 1 | 8 | 16 | 6 | -34.6021424 / -34.6021398 / -34.6021419 |

Evaluations to 0.02 eV/A. Paired ratio native precon / native BFGS: 0.50
(Si8, 3/3 wins) and 0.50 (AlP); warm-vs-warm wall ratio 0.49 and 0.55 (Mac
observation). Native precon / A3: 1.14 and 1.33; the extra evaluation is the
first step under the stiff BFGS prior (A3's fixed 5 eV/A^2 happens to be close
to the Si scale), after which mu settles at 0.07-0.10 Ha/Bohr^2. No failures;
all end energies agree within 3e-6 Ha (all stopped at the same force
tolerance). The ignored CLI test
`relaxation::precon_lbfgs_needs_fewer_evaluations_than_bfgs` (Si8, own stop
rule, fixed displacement pattern): 9 against 17 evaluations, same minimum
(dE 4.1e-7 Ha, max dx 2.9e-3 Bohr).

Variants of item 2 on Si8 seeds 1-3 (evaluations to 0.02 eV/A): prior +
y.P1^-1 y / s.y 8/9/8 (chosen); prior + s.y / s.P1 s 10/8/9; prior + mu from
the first pair only 9/8/9; first step scaled to the cap + y.P1^-1 y / s.y
8/7/8; cap + s.y / s.P1 s 9/8/8. The spread is one evaluation and not
decisive; on the Morse model the prior was better (22/17/10 against 38/15/9
with the cap-scaled first step), so the documented prior stays. The A3 gap
(one evaluation) suggests a better first-step scale; candidates for the dev08
rerun are recorded here, not tuned on silicon (generality rule).

## Generality and open items

- Shape constants are the Phase 1 A3 values, not fitted in this change; mu is
  estimated per run. Phase 1 covered Si (LDA, PBE) and Si-Al-P; the ticket
  requires an ionic oxide and a molecule/surface case on dev08 before any
  default discussion (design Section 10).
- Not covered: variable cell (`ions_cell`), internal coordinates for floppy
  molecules (torsions), the FIRE metric variant.
- The A3 harness arm evaluates through `baandi run` (cold); the native arm is
  warm or cold; wall comparisons only within a start class.
