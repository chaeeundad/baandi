# Reference SCF path (issue #25): design

Status: design draft by the orchestrator, 2026-11-23. User decisions of the same day: all lane/A-B work stopped; a NEW SCF driver replaces the legacy one; numerical kernels are reused only where they are good, and bold rewrites are allowed wherever they promise speed (the old code stays in git history). Implementation goes through Codex tickets; QE conventions reach implementers only through the Claude note `docs/research/qe-study/qe-density-mixing-spec.md` and the published literature.

## Why a new driver
The dev08 decomposition (docs/performance/requests/si64-gap-attribution-2026-11-20.md, architecture brief addendum) shows the ratio to QE as per-map cost x map count x non-SCF overhead: si32 1.51 x 0.83 x 1.66, si64 1.56 x 1.00 x 1.49, si128 1.35 x 1.93 x 1.31 (A/D 3.586 measured), si216 1.21 x 1.50 x 1.47. The per-map cost is closing with size; the map count (si128: 29 maps vs QE 15, maps 3-13 diverge with the standard error growing 272x while no control reacts) and the overhead outside the SCF loop (coarse ladder 18-24 % of wall at every size, verification 9-11 %) are the levers. The legacy driver (`ScfSession::advance`, about 1100 lines, plus the cutoff ladder and verification scheduling, more than 130 switches) is where both live. A new driver that follows a published, QE-equivalent fixed-point scheme lets us (1) compare every iteration with QE term by term and (2) add BAANDI's own accelerations back one at a time, each with a measured gain.

## Scope of the new path
Real-Gamma, integer occupations, spin-unpolarized, norm-conserving, LDA/PBE: the cases of the benchmark suite (Si8...Si216, skew, slab, elongated, MgO Gamma, synthetic PPs). Everything else (k points, spin, smearing, DFT+U, forces/stress outputs) keeps using the legacy driver until the new one covers it; the CLI chooses by profile.

## Architecture
A new module `crates/dft-engine/src/refscf/` (no dependency on `session.rs`, `cutoff_ladder.rs` or the switch registry beyond one profile switch):

1. `profile`: default-off execution profile `BAANDI_REFERENCE_SCF=on` (input `[scf] driver = "reference"`), resolved once, recorded in results.json (`scf.driver`). The CLI branch sits where `run.rs` chooses between `run_ladder` and `run_outcome`; output schema stays the same (`ScfOutcome`), plus a `reference_scf` block per iteration.
2. `density`: the mixed quantity in reciprocal space exactly as defined in the note (half-sphere storage for real Gamma, G=0 handling, which G shells are mixed and what happens to the rest, normalization, units), distributed by G vectors under MPI (owned G slabs from the existing grid plan; no full-grid replication).
3. `metric`: the Hartree-like inner product used for the residual norm (the SCF accuracy estimate) and for the history overlap matrix, with its exact weight in |G|^2, the G=0 rule and the half-sphere weights; one reduction per inner product, batched.
4. `broyden`: the modified Broyden/Anderson update (Johnson, PRB 38, 12807 (1988)) with the startup step, beta, history length, insertion/eviction order and the linear solve, as specified in the note; deterministic reductions (fixed order) so results are rank-count independent within rounding.
5. `driver`: the loop. Input n_1 = atomic superposition at the FINAL cutoff (no coarse ladder in the reference arm). For each iteration k: inner threshold from the QE-equivalent rule (initial value, update from the previous accuracy estimate and the electron count, floor), one fixed-point map `ScfProblem::step_session(n_k, warm, threshold)` (the existing distributed Hamiltonian/Davidson/density kernels), accuracy estimate dr2 = metric(n_out - n_in), the first-iteration redo rule, stop when dr2 is below the threshold that the standard accuracy class maps to (decision 3 of 2026-11-11: match QE's conv_thr 1e-8 Ry equivalent), otherwise mix. After the stop: BAANDI's independent fixed-point verification (certificate), counted in the delivered wall.
6. `record`: per iteration: dr2, threshold, Davidson iterations and H columns, beta and history occupancy, step norm, energy terms, wall per phase; plus run totals. This is the comparison surface with QE (estimated scf accuracy, ethr, avg # of iterations).

## Patches (issue #25)
- Patch 1 (runnable density update): modules 1-4 and the driver skeleton; a reproducible test that calls the production update routine on (a) synthetic arrays with closed-form answers and (b) QE reference arrays (input density, output density, history) generated on dev08 from a recorded QE 7.6 Si128 input; compare the next density, total charge and low-G components in common units against a documented tolerance; report the first failing term if any. Deliver: commit, CLI command, comparison table.
- Patch 2 (common start convergence): same atomic start at the final cutoff for both codes; compare output densities first (a mismatch points at potential, occupations or the inner solve), then next inputs (a mismatch points at mixing); show whether Si128's early excursion disappears and the map count at equal accuracy (QE's 15 is a comparator, never a stopping rule). Record the reference-style stop and BAANDI's certificate separately, both charged to the delivered wall.
- Patch 3 (speed): dev08 pure MPI8 matched runs, Si128 with Si64 and Si216 controls, A/A floor first, crossed A/B, effect >= 3x floor and |mean|/SE >= 3, uninstrumented external wall including startup, certificate and output; skewed/slab and non-Si PP coverage; strict 1e-10 / standard 5e-9 Ha before any default change. Then add BAANDI accelerations back one at a time (kept orbitals first, then the ladder only if it pays at equal accuracy), and rewrite kernels where the per-iteration profile shows waste (per-iteration non-column cost was 99 ms vs QE 69 ms at si64; Gram/rotate/update 2.53 s vs QE 0.78 s of non-eigen algebra).

## Rules
No QE source to implementers; equations and conventions from the note and the literature only, the note's INFERRED items verified against QE reference data. Same parallel model on both sides. No constant tuned to silicon. The legacy driver stays the default until the new path passes the gates on the generality set; then the old driver and its switches are removed in one step.

## Implemented conventions (R-SCF-1, patch 1)

The independent implementation is `dft_engine::refscf`. It is a child module
of `scf` solely to reuse private map/setup data without widening the engine API;
it neither constructs a session nor calls the cutoff ladder. The CLI selects it
with `BAANDI_REFERENCE_SCF=on` or `[scf] driver = "reference"`. Both are opt-in;
`driver = "legacy"`/an absent key with the switch off preserves the existing
numerical path. Environment `on` also selects reference with `driver="legacy"`;
explicit `driver="reference"` remains a request even with the switch off.
`[scf] mixing_beta` and `mixing_history` configure beta and m; their reference
defaults are **0.7** and **8** (m in 1..25). `reference_preconditioning="tf"`
enables Thomas-Fermi screening; `plain` is the default. The optional mixing
names `plain`, `broyden`, and `pulay` denote this one difference-history method.
There is no legacy mixer selection or damping/positivity recovery in this arm.

Only the own-words [density specification](../research/qe-study/qe-density-mixing-spec.md)
was used for QE conventions; no QE source was consulted. The method references
are D. D. Johnson, *Phys. Rev. B* **38**, 12807 (1988),
doi:10.1103/PhysRevB.38.12807; P. Pulay, *Chem. Phys. Lett.* **73**, 393–398
(1980), doi:10.1016/0009-2614(80)80396-4; G. Kresse and J. Furthmueller,
*Phys. Rev. B* **54**, 11169 (1996), doi:10.1103/PhysRevB.54.11169. The
implementation uses the specified difference-space least-squares equations;
it does not claim to reproduce every variant in these papers.

### Density and metric

For valence density only, in electrons/bohr^3,

\[
 \rho_G=N_r^{-1}\sum_r\rho(r)e^{-iG\cdot r},\qquad
 \rho(r)=\sum_G\rho_G e^{iG\cdot r},\qquad \rho_0=N_{el}/\Omega.
\]

Reciprocal vectors are columns of B, including 2 pi, with G=B m and integer
Miller triple m. Keep m1>0, or m1=0,m2>0, or m1=m2=0,m3>=0. This is a
*density* half sphere at four times the wavefunction energy cutoff, not the
smaller orbital half sphere. BAANDI basis order is retained locally. The
reference reader uses exact Miller keys and validates complete membership;
no positional correspondence or inferred QE shell-tie order is assumed.
A sphere reaching an FFT Nyquist frequency is rejected rather than aliased.

The smooth mask selects ngms history entries. All G=0 linear combinations are
carried, with zero metric weight. For ngms<ngm, the remaining dense entries
receive x+beta F and never history subtraction or TF screening. The current
CLI has only the 4x density cutoff, hence ngms=ngm; the general mask is tested.
NLCC core charge never enters the mixed vector. No mixing-step normalization,
clipping, or feasibility damping is performed. The existing atomic-start
kernel normalizes its superposition once (and clips initial Fourier ringing).

For the half sphere, excluding zero and high-G entries,

\[
 \langle a,b\rangle_{Ry}=8\pi\Omega
  \sum_{G\ne0}^{half,smooth}\frac{\Re(a_G^*b_G)}{|G|^2},
 \qquad \langle a,b\rangle_{Ha}=\tfrac12\langle a,b\rangle_{Ry}.
\]

Here Omega is bohr^3 and |G|^2 is bohr^-2. The factor two from the implicit
negative half is included. dr2 is this form on F=y-x, before any mixing or
history mutation. Each batch uses compensated local sums and **one** allgather
of scalar partials, summed in ascending rank order. Fixed rank counts receive
identical scalars; rank-count changes are permitted to change rounding. No
full density is gathered by the metric or Broyden module.

### Update and solve

At accepted mixing call k>1, slot (k-2) mod m receives the **unnormalized**
differences dF=F_previous-F and dx=x_previous-x; occupancy is min(k-1,m).
The previous input and residual are saved before the current update. k=1
has no differences and performs x_next=x+beta F. A redo or stop never calls
the mixer, hence never inserts/evicts history.

\[
 A_{ij}=\langle dF_i,dF_j\rangle,\quad c_i=\langle dF_i,F\rangle,\quad
 A\gamma=c,\qquad x_{next}=x+\beta F-\sum_i\gamma_i(dx_i+\beta dF_i).
\]

A's triangle and c share one reduction batch. The existing symmetric EVD
adapter (faer; MKL dsyevd on the optional MKL backend) supplies
A=Q Lambda Q^T and gamma=Q Lambda^-1 Q^T c. This is a symmetric spectral
factorization, different from the note's indefinite-factorization route but
solving the same system. cond2=lambda_max/lambda_min is recorded. A numerical
rank deficiency lambda_min<=p*epsilon*lambda_max, nonfinite result, or backend
failure stops with a named error; **no regularization, normalized secants,
pseudoinverse, automatic restart, or silently reduced history** is used.
Near-rank-deficiency rejection is deliberately stricter than a zero-pivot-only
check. Storage and numerical work arrays are reserved at construction.

TF acts only on the history-corrected residual, using G^2/(G^2+a), with
r_s=(3 Omega/(4 pi N_el))^(1/3), a=(12/pi)^(2/3)/r_s in bohr^-2. It does not
act on dx. The plain profile carries G=0 without special projection; TF's
zero factor applies to its residual, as in the stated formula.

### Driver, energy, and explicit patch-1 limits

The start is the atomic valence superposition at the final cutoff. Missing
PP_RHOATOM is an input error, not a uniform fallback. Each iteration uses one
existing `ScfProblem::step_session` map, with the single first-iteration redo
allowed by the specification. All ethr quantities below are in Ry:

* k=1 starts at 1e-2. If dr2<ethr*max(1,N_el), redo once at
  max(0.1*dr2/max(1,N_el),1e-13), with unchanged density and empty history.
* k=2 resets the previous cap to 1e-2, then applies the update. Thereafter,
  ethr=max(min(ethr,0.1*dr2_previous/max(1,N_el)),1e-13).
* The existing map receives ethr/2 in Ha. Stop when dr2<2*standard_threshold_ha,
  thus default 5e-9 Ha maps to conv_thr=1e-8 Ry. No extra legacy residual guard,
  ladder, stall schedule, mixer reset, or positivity damping is used by the
  new driver. Kernel internals still use their existing backend configurations;
  this patch does not rewrite the Hamiltonian/Davidson internals.

On mixing, the recorded first-order correction is

\[
 dE_{scf}=-\frac{\Omega}{N_r}\sum_r(x_{next}-y)
 [V_H[x_{next}]+V_{xc}[x_{next}]],
\]

in Ha (Ry is twice this). The recorded variational energy is
E_band - integral y(V_H[x]+V_xc[x]) + E_H[x_next]+E_xc[x_next]+E_ewald+dE_scf.
The ionic potential is density independent and cancels from this correction.
At the stop, dE_scf=0. Existing `energies` remain the direct orbital functional
of the last output y; these intermediate variational values are additive
`reference_scf` diagnostics. Even on exhaustion, the ordinary delivered
`density.bin` remains the last **output** density to keep it consistent with
that energy/orbital schema; it is not the note's unconverged QE restart file.
The computed mixed-next-input is used for every subsequent map.

**Verification TODO (issue #25 patch 2):** the existing fixed-point verifier
is coupled to session switch state, scheduling, predictor and reuse behavior.
It is not called or silently substituted here. With the default requested
verification, reaching the reference stop returns `converged=false`, an empty
(unpassed) verification report, a clear console/results reason, and CLI exit 3.
With explicit `verify="off"`, the estimator stop returns success but is never
certified. No production accuracy or speed acceptance follows from these runs.

Real Gamma, integer occupation and unpolarized NC are checked. Runtime control,
forces/stress/relaxation and a complex-Gamma request reject. Serial LDA/PBE reuse
the existing map. MPI patch 1 additionally requires the existing owned LDA
map's eligibility (no NLCC/core floor, symmetry, snapshot; M2 grid/local-field/
local-density kernels available). Unsupported MPI inputs reject explicitly;
there is no replicated-density fallback. Mixing owns G rows from the same
full-grid `DistributedPairFft` plan as `grid_dist`; transforms own real planes.
The inverse fills only the positive half and forms 2 Re(sum_half)-rho0, avoiding
conjugate full-grid replication. Existing construction/atomic start may still
hold full arrays; a full density is materialized at final output only.

Per-iteration JSON records include dr2 and ethr in Ry/Ha, beta, history slot
and occupancy, gamma norms, A condition number, Davidson/H counts **including
redo work**, and map/FFT/metric/mixing/energy-correction seconds. Stop records
have no solve condition/gamma because no history solve was performed.

### Reference files and tolerance

`refscf-check INPUT --qe DIRECTORY --steps 6` reads `k{k}-out` and `k{k}-mix`,
recursing at most four directory levels to find exactly one charge-density.dat
per snapshot. It checks four-byte Fortran markers, both endiannesses, lengths,
Gamma/nspin, finite values, reciprocal basis, unique Miller entries and exact
membership. The dev08 probe supplied to this task reports Si128 ngm=191712,
5,368,052 bytes (=116+28*ngm), and rho0=512/Omega; those are **relayed probe
facts**, not newly measured here.

The initial input is (mix1-beta*out1)/(1-beta); beta=1 recovery rejects.
Optional `k1-in` cross-checks the recovered atomic input. Optional `k1-beta1`,
`k1-beta-1`, or `beta1` cross-checks out1 (a beta=1 output alone cannot uniquely
recover the initial density). Every later input is QE's preceding mixed file,
not BAANDI's replay result. The reader recognizes stdout, stdout.txt, pw.out,
qe.out, output.out, run.out in each k-mix directory, or k{k}-mix.out/.stdout
at the root. Missing stdout estimates are null, not invented. Printed '<'
is an upper bound; '=' uses half the last printed decimal's resolution.

Checks are ordered: file/cell/membership, optional initial checks, per-k dr2
when available, A/gamma solve, all G, zero mode, lowest ten distinct nonzero
shells, total charge. It exits at the first failed term and names it. Output
contains max/RMS errors in e/bohr^3, relative errors, charge in electrons,
shell G^2 in bohr^-2 and dr2 in Ry/Ha. Absolute/relative defaults (1e-12
e/bohr^3, 1e-9) are **provisional diagnostic tolerances**, adjustable by `--atol`
and `--rtol`, not a measured QE agreement claim. Real QE tolerances await the
rank-count/repeated-run evidence requested in [the dev08 protocol](../performance/lanes/RSCF/dev08-request.md).

Synthetic three-step replay (including reversed file order and simulated MPI3)
measured max density difference 2.168404344971009e-19 e/bohr^3, maximum relative
shell error 2.168404344971009e-16, and zero charge error. This is a closed-form
linear fixed-point oracle, **not QE data**. The real Si8 fixture test emits
`SKIP QE si8 reference` while the delivery branch/files are absent. Detailed
Mac evidence and default-path fingerprints are in
[local-validation.json](../performance/lanes/RSCF/local-validation.json).

## Delivered QE regression (2026-11-23 follow-up)

Supersedes the earlier absent-fixture status. Copied all 46 files (2,806,532
bytes) byte-for-byte from the main worktree; no UPF is included. Provenance,
QE version, commands and hashes remain in the unmodified fixture README,
MANIFEST and SHA256SUMS. This is density-update correctness evidence, not
Mac performance evidence or independent full-SCF trajectory agreement.

`cargo test -p dft-engine --test refscf -- --nocapture` passes all 10 tests,
including the real production Broyden replay. Initial densities recovered with
beta 0.7 and 0.3 agree to 3.46945e-18 e/bohr^3 (relative max 1.17140e-16);
beta 1 output agrees to 3.83323e-20 e/bohr^3. Each iteration uses QE's own x/F;
history occupancy is 0, 1, 2. All 11,936 Miller indices are matched in BAANDI
norm/Miller ordering. The ten lowest nonzero shells are checked separately.

| k | max density error (e/bohr^3) | RMS density error (e/bohr^3) | G0 error (e/bohr^3) | max error over lowest 10 shells (e/bohr^3) | dr2 (Ry) | dr2 minus printed (Ry) |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 3.46944695e-18 | 3.18984922e-20 | 3.46944695e-18 | 1.53329342e-19 | 0.438555331969 | 1.96919914e-09 |
| 2 | 3.46944695e-18 | 3.35435183e-20 | 3.46944695e-18 | 8.67361738e-19 | 0.019564039712 | -2.88015854e-10 |
| 3 | 8.67361738e-19 | 1.10392361e-20 | 0.00000000e+00 | 8.67361738e-19 | 0.000853988046877 | -1.95312339e-09 |

Charge differences versus QE are -3.74849e-15, -3.74849e-15, 0 electrons;
absolute total charge agrees with 32 electrons within 1.43e-14 electrons.
Maximum shell relative max/RMS errors are 6.4000e-16 / 9.6936e-16.
Regression limits: max and RMS density errors <= 1e-19 + 1e-13 times the
reference group's maximum amplitude, and charge error < 1e-12 electrons.
These give rounding headroom over the measured errors without accepting
physics-scale discrepancies. dr2 must match the printed rounded value within
half its last decimal place plus 1e-14 Ry, even though QE prints '<'.
CLI general diagnostic defaults remain unchanged; this fixture is tested at
the tighter measured limits. Native MPI/MKL Si128 remains a dev08 check.

CLI command (replace the replay TOML's UPF placeholder with the local file):

```sh
target/release/dft-cli refscf-check target/scratch/si8-qe-replay.toml \
  --qe tests/reference/qe/i25-si8/si8-np1 --steps 3 --atol 1e-19 --rtol 1e-13
```

Exit 0; no failing term. Machine-readable per-shell errors and diagnostics:
[qe-si8-validation.json](../performance/lanes/RSCF/qe-si8-validation.json), [qe-si8-cli.jsonl](../performance/lanes/RSCF/qe-si8-cli.jsonl).
Final full gate: exit 1, 813 passed / 3 known cache failures / 181 ignored;
build, fmt and clippy pass. The orchestrator explicitly accepts those TMPDIR
failures as an environment limitation. Fixed-point verification remains the
accepted patch-2 TODO. No commit was created.
