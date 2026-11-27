# Preconditioned ionic relaxation: design note (PRELAX-PRECON, Phase 0/1)

Status: design input for Phase 1 (benchmark) and, only if Phase 1 shows a
measurable DFT benefit, Phase 2 (native Rust optimizer). Binding requirements:
[user specification](precon-relax-spec.md). Current-code facts: [audit](precon-relax-audit.md).
Nothing here changes production code.

## 1. Problem statement in one paragraph

Ionic relaxation minimises E(R) with forces F = -dE/dR from a verified SCF per
geometry. Each evaluation costs one full SCF (seconds to hours), so the figure of
merit is the number of attempted evaluations to reach max_i |F_i| < 0.02 eV/A
(3.889e-4 Ha/Bohr), and then end-to-end wall. Quasi-Newton methods (BFGS,
L-BFGS) converge in a number of steps that grows with the condition number of
the Hessian *as seen through the initial inverse Hessian H0*. With H0 = alpha^-1 I
(the current native BFGS and ASE BFGS/LBFGS), that condition number is the
ratio of stiffest bond-stretch to softest collective mode (long-wavelength
acoustic-like modes, soft molecule/surface modes, rotations of adsorbates),
which grows with system size and with heterogeneity. A preconditioner P ~ H that
is cheap, sparse and SPD replaces H0 by P^-1 and removes most of that spread.

## 2. Literature (what is used, in our own words)

- **Packwood et al., J. Chem. Phys. 144, 164109 (2016), doi:10.1063/1.4947024**
  ("universal preconditioner for condensed phase materials"). P is a weighted
  graph Laplacian over atom pairs within r_cut, with off-diagonal weight
  -mu exp(-A (r_ij / r_nn - 1)) (A = 3, r_nn the nearest-neighbour distance,
  r_cut = 2 r_nn), Kronecker-expanded to 3 Cartesian components (isotropic
  blocks), plus a diagonal stabilisation mu c_stab (c_stab = 0.1) that makes it
  SPD despite the Laplacian's translation null space. The energy scale mu is fitted
  once by a long-wavelength test displacement v (a sine of the Cartesian
  coordinates over the cell): mu = <grad E(x+v) - grad E(x), v> / <P_1 v, v>,
  with P_1 the matrix for mu = 1; it costs one extra force evaluation. It is used
  as the initial inverse metric of L-BFGS (two-loop recursion with P^-1) or as the
  metric of FIRE, with an Armijo line search for L-BFGS. The reported benefit grows
  with system size and with soft modes (defects, surfaces, interfaces, vacuum);
  for small stiff bulk cells the gain is small. Variants in the same family: a
  constant-weight Laplacian ("C1"), force-field based P (bond/angle/dihedral
  Hessians of a model potential) for molecules (Mones, Ortner, Csanyi, Sci. Rep. 8,
  13991 (2018)), and P for path methods (Makri, Ortner, Kermode, J. Chem. Phys.
  150, 094109 (2019)). Bibliographic details other than Packwood's DOI are to be
  verified before they enter user documentation.
- **Model Hessians / internal coordinates for molecules**: Lindh et al. model
  Hessian (Chem. Phys. Lett. 241, 423 (1995)): bond, angle and torsion force
  constants with exponential distance dependence from element-row parameters;
  redundant/delocalised internal coordinates (Pulay and Fogarasi; Baker; Bakken
  and Helgaker) and translation-rotation internal coordinates (geomeTRIC, Wang and
  Song, J. Chem. Phys. 144, 214108 (2016)). These are strongest for isolated
  molecules and weakly bound adsorbates (rotations, torsions); in periodic solids a
  Cartesian sparse P is simpler and sufficient. They are a later option, not part of
  Phase 2's first cut.
- **Globalisation**: line search (Armijo sufficient decrease; Wolfe adds a
  curvature condition that guarantees s.y > 0 for BFGS updates) versus trust
  region (step limited to a radius updated from the ratio of actual to predicted
  decrease). Both need energies consistent with forces; both must tolerate SCF
  noise in E of order the SCF threshold (Nocedal and Wright, Numerical
  Optimization, ch. 3, 4, 7).
- **Pfrommer et al., J. Comput. Phys. 131, 233 (1997)**: physically scaled
  initial Hessian for solids (bulk modulus / phonon frequency); ASE ships it as
  a preconditioner. Mentioned for completeness; it is a scalar H0.

## 3. ASE PreconLBFGS (LGPL-2.1): behaviour, studied, not copied

Read in the installed ASE 3.26 (`ase/optimize/precon/{lbfgs,precon}.py`,
`ase/utils/linesearcharmijo.py`) to describe behaviour; no code enters this MIT
tree. ASE is used only as an external benchmark through the Phase 1 harness.

- Default preconditioner `Exp(A=3)`: r_nn estimated from the structure,
  r_cut = 2 r_nn, c_stab = 0.1, sparse assembly; mu is estimated at the first
  preconditioner build by the sinusoidal test displacement described above, which
  calls the calculator on a perturbed COPY of the atoms: one extra DFT evaluation
  (the harness labels it `auxiliary` and counts it). mu is floored at 1 eV/A^2.
- P is rebuilt only when some atom moved more than 0.5 r_nn since the last build;
  mu is not re-estimated unless `reinitialize` is set.
- L-BFGS: memory 100, two-loop recursion with z = P^-1 q in the middle.
- Armijo line search with c1 = 0.23 (a strong sufficient-decrease requirement):
  initial step length from a quadratic interpolation of the previous energy
  decrease (falls back to 1), rescaled so that no atom moves more than `maxstep`
  (ASE default 0.04 A if not given; the harness passes the case's cap), then
  backtracking by quadratic interpolation, never below a1/10. Each trial is one
  DFT evaluation. A failed line search resets the history once, then raises.
- Convergence: max atomic force norm (and stress for variable cell).
- With a single-state calculator, the main structure is forgotten after the mu
  probe, so ASE would re-evaluate an identical geometry; the harness calculator
  memoises exact geometries so that no identical SCF is ever run twice (hits are
  reported as `memo_hits`, not as evaluations).

## 4. The experimental A3 optimizer (exp-tabpfn)

Source (outside the repo): `exp-tabpfn/tools/icl/run_precon.py` (driver, mode
`A3`), `atomistic.py` (base metric), `preconditioner.py` (metric algebra; the
other modes A4..A9 add learned or multisecant corrections and are not part of
this task). A3 is:

- a FIXED model Hessian B_p assembled once from the starting geometry over
  element-aware neighbour pairs: cutoff 1.3 (r_cov,i + r_cov,j), 1.8 for pairs
  of a molecular element (H, C, N, O) with a non-molecular one (adsorbate to
  surface), periodic minimum image; stiffness k = 5 exp(-2.5 (d/r0 - 1)) eV/A^2
  with r0 the covalent-radius sum; anisotropic 3x3 bond block
  k (0.3 I + 0.7 n n^T) (stretch-dominant, unlike Exp's isotropic blocks);
  graph-Laplacian assembly; fixed atoms removed; plus 0.05 I eV/A^2;
- its exact dense inverse (eigendecomposition) as the L-BFGS H0, memory 10;
- Armijo with c1 = 1e-4, halving, at most 25 trials; non-descent fallback to
  -B_p^-1 g; curvature-pair filter s.y > 1e-12 |s||y|;
- no mu estimation (absolute scale fixed by the constant 5 eV/A^2), never rebuilt.

EMT evidence (`exp-tabpfn/cl/stage2_emt.json`, `stage2_preconlbfgs.json`,
12 seeds each, fmax norm 0.02 eV/A, all runs converged; recomputed here):

| system | ASE BFGS | ASE LBFGS | ASE FIRE | ASE PreconLBFGS(Exp) | A3 |
|---|---|---|---|---|---|
| Cu(111)+CH3OH | 172 | 158.5 | 270 | 88 | 69 |
| Ni(111)+C2H4 | 185.5 | 184 | 300.5 | 99 | 99.5 |
| Au/Cu alloy (32 atoms) | 21 | 21 | 49 | 12 | 13 |

(median evaluations). Two caveats found while recomputing: (1) on both
adsorbate systems most seeds end in **different minima** depending on the method
(|E_A3 - E_BFGS| median 0.025 eV, up to 1.29 eV; PreconLBFGS similar), so
evaluation counts there compare different end states and must be paired with a
same-minimum check; on the alloy all methods agree within 2e-4 eV. (2) EMT is a
cheap pair-functional surrogate; it says nothing about SCF noise, failures or
warm-start effects. The spec already states that this does not establish a DFT
benefit.

Generality note (AGENTS.md generality rule): A3's fixed 5 eV/A^2 and the
1.3/1.8 cutoffs are model constants, not fitted to one material, but the absolute
scale is not estimated; a native version must estimate the scale (Section 6.3).

## 5. Phase 1 benchmark design (implemented: `tests/reference/relax-precon/`)

Methods (spec): BAANDI BFGS, BAANDI FIRE (`baandi relax`), ASE BFGS, ASE LBFGS,
ASE FIRE, ASE PreconLBFGS(Exp), A3 + L-BFGS (harness re-expression of the
experimental A3 in `relaxbench/a3.py`; own code, no ASE code copied).

Fairness (spec section 6), as implemented:
- identical structure (seeded start), UPFs (SHA-256 recorded), XC, cutoff, k
  points, bands, SCF settings (the same `[basis]/[xc]/[scf]` text for every
  arm), threads (`--threads`, the case's own thread count is stripped), fixed
  atoms and the same per-atom displacement cap (the case's `max_displacement_bohr`);
- **one common external evaluator**: after every attempted evaluation, a
  candidate geometry (not an auxiliary probe) whose max atomic force NORM over
  free atoms is below 3.889e-4 Ha/Bohr stops the run; intermediate thresholds
  0.1 and 0.05 eV/A are recorded as first-hit evaluation counts. The native stop
  rule (max Cartesian component + energy change, twice) is disabled by a
  1e-9 Ha/Bohr force tolerance and the process is terminated by the harness; the
  native energy tolerance is kept because it drives the BFGS trust radius;
- every attempted evaluation is counted: failed SCFs (native retry included),
  rejected line-search trials (geometry never accepted by the optimizer),
  auxiliary probes (ASE mu estimation);
- **warm and cold are never mixed in wall comparisons**. `baandi run` has no
  start-density input, so all ASE/A3 arms are cold (fresh atomic start per
  evaluation). Native arms run both cold (`warm_start=false`) and warm. The
  summary prints wall ratios only within one start class; evaluation-count
  ratios are reported across classes (they do not depend on the start except
  through SCF noise and failures);
- per evaluation: role, accepted/rejected, failed + reason, wall, SCF
  iterations, H vectors, energy, max atomic force norm and max component;
  per run: status, attempted/failed/rejected/auxiliary counts, totals, end-to-end
  wall, final energy/force/positions, executable and UPF SHA-256, host, ASE version.
- ASE optimizers stop only through the evaluator (`fmax` set to 1e-12).
- Paired summary (`summarize.py`): medians, per-seed paired ratio of
  evaluations to 0.02 eV/A against a baseline arm, wins/ties/losses, and a
  different-minimum flag (|dE| > 1e-4 Ha or max final displacement > 0.05 Bohr).

Known asymmetries (recorded, not hidden): the native relax retries a failed SCF
once at half displacement, ASE arms stop on the first failure; the native warm
arm reuses density and orbitals, which a Phase 2 native preconditioned optimizer
would also do, so the decisive wall comparison is native-warm versus a native
preconditioned optimizer (Phase 2), while Phase 1 decides on evaluation counts
(cold vs cold for wall). Peak memory and preconditioner build/apply times are not
recorded in Phase 1 (ASE internals); they are Phase 2 metrics.

Cases (existing inputs only, `cases/*.json`):

| case | atoms | settings (from) | starts |
|---|---|---|---|
| si8-lda-k222-a0.2 | 8 | 24 Ry LDA-PZ, k 2x2x2, 20 bands (`relax_cases.py`) | seeded, 0.2 Bohr |
| si8-lda-k222-a0.05 | 8 | same | seeded, 0.05 Bohr (trivial, ~5 evals) |
| si8_AlP-pbe-k222 | 8 | 36 Ry PBE, k 2x2x2, PseudoDojo (`relax/sys/si8_AlP-s{1,2,3}`) | 3 pre-perturbed (0.1 Bohr) |
| si_alp_vac-pbe-g48 | 63 | 48 Ry PBE, Gamma, 130 bands (`relax/sys/si_alp_vac-s{1,2}`) | 2 configurations (0.2 Bohr) |
| si_alp_vac2-pbe-g48 | 63 | 48 Ry PBE, Gamma, 129 bands (`relax/sys/si_alp_vac2-s{1,2}`) | 2 configurations |

Surfaces/adsorbates (TiO2(110)+CH3OH, `exp-tabpfn/tools/tio2.py`) are added only
after a measurable benefit on the bulk/defect cases (spec Phase 1 order).

## 6. Phase 2 design (conditional on Phase 1 evidence)

### 6.1 Placement and interface
`crates/dft-cli/src/relax/` (driver layer, per crates/dft-cli/AGENTS.md); a new
`method = "precon_lbfgs"` value in the existing `[relax]` table, no new
environment variable, defaults unchanged until established. Reuses the existing
evaluation, warm-start (`IonicSeed`), retry, constraint projection, restart and
output paths, so warm starts are preserved by construction.

### 6.2 Preconditioner
P = mu (L_w (x) blocks) + mu c_stab I, restricted to the active (free)
coordinates:
- graph: element-aware neighbour pairs (covalent-radius cutoffs as A3, or
  r_cut = 2 r_nn per species pair as Exp; Phase 1 decides which), periodic
  minimum image with cell lists (O(N)), sparse CSR;
- weights exp(-A (r/r0 - 1)), A = 3; blocks isotropic (Exp) or stretch-dominant
  (A3) — again decided by Phase 1 (A3 vs PreconLBFGS);
- SPD on the active subspace: fixed atoms are eliminated (Dirichlet rows and
  columns removed), which pins every component containing a fixed atom; c_stab > 0
  pins every component without one (bulk translation, disconnected molecule in
  vacuum, isolated atom). Component-wise constraints (fixed Cartesian components)
  are eliminated per coordinate the same way;
- solve: sparse Cholesky of the active block (faer sparse LLT with fill-reducing
  ordering; adopted crate) factorised at each rebuild; apply = two triangular
  solves, O(nnz). For very large N, PCG with incomplete Cholesky as fallback;
- rebuild when any atom moved more than 0.5 r_nn since the last build or the
  neighbour graph changed; L-BFGS pairs are kept across rebuilds because the
  two-loop recursion accepts a different H0 at every iteration (the pairs are
  secant data of E, not of P); the curvature filter keeps the implicit inverse
  SPD. A reset of the history is triggered only by a failed line search.

### 6.3 Scale mu without extra evaluations (to be validated)
Packwood's probe costs one evaluation. Alternatives to test: (a) first step
scaled to the trust radius (as BFGS today), then mu from the Rayleigh quotient
of the first accepted secant pair, mu = (s.y)/(s.P_1 s), bounded to a factor
[0.1, 10] change per update; (b) Packwood's probe on the first SCF with a warm
start (cheap SCF). Adoption needs measured evidence across the case set
(generality rule), never a constant tuned on silicon.

### 6.4 L-BFGS and globalisation
Two-loop recursion with z = P^-1 q, memory m = 20 (Phase 1 shows whether 10/100
matters); descent check g.p < 0 else p = -P^-1 g and history reset; per-atom
displacement cap = min(trust radius, max_displacement_bohr) (unchanged
semantics); Armijo sufficient decrease with a noise tolerance:
E(x + t p) <= E(x) + c1 t g.p + eps_E, eps_E = a small multiple of the SCF energy
threshold (the standard threshold is 5e-9 Ha); backtracking by safeguarded
quadratic interpolation; at most 3 trials before a history reset and fallback to
the preconditioned steepest-descent step; curvature check s.y > 0 before storing
a pair. Rejected trials are evaluations: the warm seed used for the next trial is
the nearest evaluated geometry (the existing seed is from the last successful
evaluation).

### 6.5 Robustness (spec section 8)
SCF failure: the existing single retry at half displacement from a fresh start,
then fail with last good structure; a failed trial counts as an attempted
evaluation and shrinks the trust radius; convergence only on real verified forces;
final report flags a different minimum (energy/structure) against a reference run
when benchmarking.

### 6.6 Restart/checkpoint
Restart payload adds the L-BFGS pairs (m x 2 x 3N floats), mu, c_stab, A and
the positions at the last P build; P itself is rebuilt on resume (deterministic
from positions), so the restart stays small. Same IEEE-bits encoding and
checksum as today.

### 6.7 Tests (spec section 12), planned
SPD and inverse accuracy of P on random structures (|P P^-1 v - v| < 1e-10),
translation/rotation/permutation invariance, PBC images (small cells where
neighbours are periodic images), fixed atoms and fixed components (SPD, zero
steps on fixed coordinates), disconnected fragments (molecule in vacuum), weakly
bound adsorbate (large condition number), sparse build/solve scaling O(N) to
1e4 atoms, L-BFGS history across rebuilds, SCF-failure recovery, restart
consistency (byte-identical restart from a committed boundary), analytic
quadratic/Rosenbrock oracles, BFGS/FIRE regression unchanged.

## 7. Phase 1 decision rule

Proceed to Phase 2 only if, on the DFT cases with at least 3 seeds each, A3 or
PreconLBFGS(Exp) reduces median attempted evaluations to 0.02 eV/A by at least
20 percent against both native BFGS and ASE LBFGS, with no additional failures
and no unflagged different minima; otherwise document the negative result and
stop (spec acceptance clause).

## 8. Harness environment and calculator choice

- Python: any external venv with ASE >= 3.23 and SciPy (Mac: the existing
  `/Volumes/990-2TB/Develop/rust-dft-wt/venv-convert`, Python 3.9.6, ASE 3.26.0,
  NumPy 2.0.2, SciPy 1.13.1; no installation was made). Nothing in the engine or
  its build depends on it.
- The BAANDI Python binding (`crates/baandi-py`, `python/baandi/ase.py`) was
  not used: its ASE calculator builds a new `Session` per calculation (cold
  start, same as the CLI path, so no warm-start advantage), does not forward
  `Options.kgrid` (Gamma only through ASE), and has no XC selection (LDA
  workflows), so it cannot run the PBE k-point cases. Building the wheel
  (maturin is only in `../pyenv`) would add work without changing the cold/warm
  situation. The harness therefore uses `baandi run` per evaluation
  (`relaxbench/clicalc.py`): every BAANDI input is accepted, the cost includes
  process start, UPF parsing and setup per evaluation (not separated by the
  harness; a whole Si8 Gamma `baandi run` took 0.9 s on the Mac, so the
  overhead is small against 2-20 s SCFs but it is a cold-path cost that the
  native warm path does not pay).
- `--launcher` puts the same command prefix (for example `taskset -c N`) in
  front of every BAANDI process of every arm.

## 9. dev08 Phase 1 benchmark plan (proposal; acceptance numbers come from dev08)

- Executable: release build of the merge commit of this branch (no engine
  change), one core per run (`--threads 1`, pinned with
  `--launcher "taskset -c <core>"`); independent runs fan out over distinct cores,
  at most one run per physical core, no SMT siblings, so evaluation counts are
  unaffected and walls are comparable within a lane. MPI is not used in Phase 1
  (the relaxation arms would have to agree on one rank layout through
  `mpirun`; possible with `--launcher`, not needed for evaluation counts).
- Arms (all cases): baandi-bfgs:warm, baandi-bfgs:cold, baandi-fire:warm,
  ase-bfgs:cold, ase-lbfgs:cold, ase-fire:cold, ase-precon-lbfgs:cold,
  a3-lbfgs:cold (8 arms); plus one `--native-own-stop` baandi-bfgs:warm run per
  seed for the stop-rule audit.
- Cases and seeds:
  1. si8-lda-k222-a0.2: seeds 1-10 (cheap; distribution of ratios);
  2. si8_AlP-pbe-k222: the 3 existing starts;
  3. si_alp_vac-pbe-g48 and si_alp_vac2-pbe-g48: the 2 existing starts each
     (63 atoms; the case where conditioning should matter most);
  4. after a measured benefit only: TiO2(110)+CH3OH (fixed bottom layer;
     `exp-tabpfn/tools/tio2.py`), 2-3 azimuth seeds.
- Budget: max 200 evaluations per run. Expected cost (Mac per-evaluation
  observations, one-core dev08 assumed 2-4x slower than 4 Mac threads):
  Si8 ~2-3 s per SCF, 80 runs x ~20 evaluations, about 1-2 core-hours;
  Si8-AlP ~10-20 s, 24 runs x ~20, about 2-3 core-hours; 63-atom cases
  several minutes per cold SCF (280 s on 4 Mac threads with an older build),
  32 runs x ~30-60 evaluations, roughly 150-400 core-hours, i.e. one to two
  days on 16-32 cores in parallel. Run the 63-atom set last and only for the
  arms that remain competitive after cases 1-2 if the budget is tight.
- Outputs: `tests/reference/relax-precon/results/dev08-1core-<date>/` raw JSON,
  `summarize.py --baseline baandi-bfgs:warm` and `--baseline baandi-bfgs:cold`
  summaries; the decision rule of Section 7.
