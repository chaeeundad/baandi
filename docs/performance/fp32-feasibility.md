# P47B: early-SCF single-precision feasibility

Measurement only, 2026-11-08. Production sources, engine options, tolerances,
certificates, manifests, README and artifacts are unchanged. All timings are
single-thread Mac observations, not dev08 acceptance or QE comparisons.

**Decision: do not recommend an implementation ticket.** At the empirical
standard/Ritz eligibility threshold of 1e-5 Ha, the generous median wall-saving
ceiling is **6.78% ordinary / 6.18% native**, below the frozen 10% rule.

## Question and frozen decision

Recommend an implementation ticket only if the upper bound on the median
end-to-end wall saving over equally weighted Si8, Si16 and Si32 is at least
10%. A kernel speedup is not a wall-time saving: if `s = t64/t32`, the fraction
of eligible work removed is `1 - 1/s`, not `s` or `s - 1`. For eligible H time
`H_e` and total wall time `W`, the fixed-trajectory estimate is
`(H_e/W) * (1 - 1/s)`. Changed iteration counts, conversion overhead, rebuilding
cached H images and additional fp64 verification can only reduce this idealized
saving if the fp64 trajectory is held fixed. A changed SCF trajectory can go
either way and is not predicted by this model.

The premise is literature table row 1 and section 6 in
[diagonalization-literature.md](../research/diagonalization-literature.md).
This is not a repeat of rejected FFT lane batching (A11), band grouping (A28),
thread spawning (A34/A35), or nonlocal traversal (A43): the independent variable
is scalar precision at identical shapes and algorithms. No external-code source
was read and no new literature performance claim is made.

## Reproduction and scope

Two ignored integration-test harnesses contain all Rust experiment code:

- `crates/dft-engine/tests/fp32_snapshots.rs`: unmodified standard SCFs at T=1,
  uniform start plus two small deterministic density perturbations; the existing
  P14A collector independently records every solve's tolerance, Davidson steps,
  H vectors and block calls. Saves early and converged physical Si orbitals and
  H components through public APIs. The converged snapshots are **standard,
  uncertified**, not strict certificate outputs.
- `crates/dft-numerics/tests/fp32_feasibility.rs`: rustfft f32/f64 full-grid lanes,
  real pairs, faer sequential nonlocal products and Davidson GEMM shapes. Plans,
  inputs and output buffers are constructed outside the timed loop. Library
  internal GEMM packing is included. Storage conversions are excluded. Each
  condition has 30 alternating-order pairs after warmup; raw timings are kept.

`PairFusedPlan` is hard-coded to `Complex64` and `Fft<f64>`; it is not generic.
The ticket explicitly permits rustfft lanes in this case. Consequently these
FFT ratios are precision observations on the fallback algorithm, **not measured
speedups of the production pruned mirror-slab plan**. The experiment's fp64
complete H output is checked against the current production `RealHamiltonian`
on every saved column (absolute norm difference < 1e-11 Ha). The fp32 experiment
rounds inputs, potential, kinetic diagonal, projectors and coupling to f32;
FFT, products and H summation are all f32. Residual/error analysis is f64.

The standard schedule is a Ritz-energy-change threshold, not a residual norm:
`tau_1 = 5e-3 Ha`, then
`max(5e-14 Ha, min(tau_previous, 0.1 E_Hartree_error / N_e))`.
Buffer thresholds have the existing `max(5 tau, 5e-6 Ha)` rule. This study
keeps the occupied/first-empty threshold as the eligibility criterion; the
buffer rule never authorizes reduced precision in a certificate.

Raw binary snapshots, logs, build failures and corrections, source provenance,
recorded historical histories, CSV summaries and analysis script are in
`P47B-evidence/` in this worktree. UPF files are external and are not copied.
The snapshot binary is an internal measurement exchange: little-endian u64
`nx,ny,nz,N,b,p`, two `ceil(N/2)` flattened FFT-index arrays (plus, minus), then
little-endian f64 arrays `V[Ngrid], kinetic[N], X[N*b], P[N*p], D[p*p], H64X[N*b]`.
Matrices are column-major; X and projectors use orthonormal real Gamma
coordinates. Potential, kinetic, coupling and H outputs are in Hartree units.
No claim of a user-facing scientific output format is made.

Commands (run from this worktree, without simultaneous builds and timings):

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
P14A_WORK_ACCOUNTING=1 cargo test --offline --release -p dft-engine \
  --test fp32_snapshots -- --ignored --nocapture --test-threads=1 \
  > P47B-evidence/engine-final.log 2>&1
# Optional counting/timing rerun without regenerating snapshots:
P14A_WORK_ACCOUNTING=1 P47B_SKIP_SNAPSHOTS=1 cargo test --offline --release \
  -p dft-engine --test fp32_snapshots -- --ignored --nocapture --test-threads=1 \
  > P47B-evidence/engine-wall.log 2>&1
P14A_WORK_ACCOUNTING=1 P47B_SI32_N=11877 cargo test --offline --release \
  -p dft-numerics --test fp32_feasibility -- --ignored --nocapture --test-threads=1 \
  > P47B-evidence/ordinary.log 2>&1
CARGO_TARGET_DIR=target/p47b-native RUSTFLAGS='-C target-cpu=native' \
  P14A_WORK_ACCOUNTING=1 P47B_SI32_N=11877 cargo test --offline --release \
  -p dft-numerics --test fp32_feasibility -- --ignored --nocapture --test-threads=1 \
  > P47B-evidence/native.log 2>&1
python3 -I P47B-evidence/analyze.py
python3 -I P47B-evidence/decision.py
python3 -I P47B-evidence/tables.py
scripts/check.sh
```

The UPF defaults to `../pseudo/Si.pz-vbc.UPF`; `P47B_UPF` overrides that path.
The repository check script accepts a crate name, not the literal argument
`all`; running it without arguments is the all-workspace gate.

<!-- P47B measured tables follow; regenerated from P47B-evidence. -->
## Measured speed

Times below are medians in microseconds; speedup is the median paired `t64/t32`.
FFT rows contain one inverse plus one forward full-grid transform (one real pair),
including lane gathering/scattering and normalization, without potential or pair packing.
Complete local/H rows additionally include pair packing, potential multiplication,
kinetic and nonlocal work as appropriate. Raw CSVs include every measured shape and pair.

| Grid | ordinary f64 us | ordinary f32 us | speedup | native f64 us | native f32 us | speedup |
|---|---|---|---|---|---|---|
| 25x25x25 | 303.30 | 220.71 | 1.374 | 230.35 | 196.83 | 1.202 |
| 32x32x32 | 436.11 | 322.11 | 1.345 | 342.51 | 247.85 | 1.367 |
| 36x36x36 | 803.08 | 616.45 | 1.315 | 642.92 | 526.34 | 1.277 |
| 36x36x72 | 1615.64 | 1320.60 | 1.281 | 1175.25 | 924.38 | 1.276 |
| 48x48x48 | 2701.40 | 2054.65 | 1.376 | 1545.75 | 1126.48 | 1.390 |
| 72x36x36 | 1674.58 | 1237.98 | 1.345 | 1200.63 | 954.03 | 1.266 |
| 72x72x36 | 5839.56 | 2803.03 | 1.865 | 2318.40 | 1740.27 | 1.334 |

Actual 30 Ry benchmark shapes are Si8: N=2969, b=20, p=32, 36³;
Si16: N=5915, b=36, p=64, 36×36×72; Si32: N=11877, b=68, p=128,
72×72×36. N is the full-sphere/ORC dimension, b the band count and p the
expanded projector count. Requested 25³, 32³ and 48³ are additional FFT
shape probes, not mislabeled production grids; Si8 at this cutoff cannot use
25³ or 32³ as the alias-free density grid. 72×36×36 tests the transposed elongated shape.

Complete H timings use actual saved orbitals; nonlocal is `P (D (P^T X))` with
faer `Par::Seq`. D is stored as a dense block-diagonal matrix in this experiment,
so the coupling call also visits off-atom zeros. The production adapter uses
small atom blocks. These are precision comparisons of the test implementation.

| Operation/block | ordinary f64 us | ordinary f32 us | speedup | native f64 us | native f32 us | speedup |
|---|---|---|---|---|---|---|
| h-si8-early | 9768.48 | 6853.72 | 1.353 | 6276.61 | 5057.97 | 1.240 |
| local-si8-early | 6196.01 | 5083.36 | 1.201 | 6076.33 | 4907.70 | 1.241 |
| nonlocal-si8-early | 146.59 | 80.20 | 1.828 | 148.93 | 81.13 | 1.831 |
| h-si8-converged | 6275.91 | 5129.51 | 1.212 | 6234.83 | 5042.49 | 1.237 |
| local-si8-converged | 6253.64 | 5072.82 | 1.235 | 6073.49 | 4880.52 | 1.244 |
| nonlocal-si8-converged | 170.79 | 93.65 | 1.824 | 149.77 | 81.65 | 1.834 |
| h-si16-early | 23179.55 | 18740.74 | 1.233 | 22803.80 | 17628.43 | 1.296 |
| local-si16-early | 21245.46 | 17358.47 | 1.221 | 21401.16 | 16704.61 | 1.278 |
| nonlocal-si16-early | 1061.36 | 622.68 | 1.697 | 1024.47 | 567.01 | 1.791 |
| h-si16-converged | 22411.28 | 17747.04 | 1.245 | 22675.94 | 17524.90 | 1.294 |
| local-si16-converged | 21295.94 | 17230.84 | 1.230 | 21447.99 | 16706.09 | 1.282 |
| nonlocal-si16-converged | 1001.36 | 574.14 | 1.748 | 1056.55 | 587.21 | 1.806 |
| h-si32-early | 88246.56 | 66763.72 | 1.332 | 87853.93 | 64875.18 | 1.359 |
| local-si32-early | 79561.21 | 62738.70 | 1.295 | 79256.95 | 59749.07 | 1.323 |
| nonlocal-si32-early | 7220.72 | 3942.69 | 1.822 | 7605.97 | 3973.65 | 1.900 |
| h-si32-converged | 88574.34 | 68149.51 | 1.303 | 84002.49 | 62630.81 | 1.342 |
| local-si32-converged | 80360.20 | 62952.58 | 1.286 | 76399.38 | 58286.55 | 1.313 |
| nonlocal-si32-converged | 7261.64 | 3981.91 | 1.835 | 6998.80 | 3820.80 | 1.826 |

Davidson representative faer shapes below use q=2b. The raw CSVs additionally
cover q=b and q=6b, correction widths 2 and b, rotation, Gram/projection products,
and nonlocal projection/reconstruction widths 2 and b. Complete H also covers
block widths 1, 2, 4, 8 and b. Rotation shapes repeated in the two correction-width
sweeps have 60 paired samples; other shapes have 30. Davidson products are reported
for feasibility but **not counted as H-only wall savings**. No fp32 orthogonalizer
or eigensolver is proposed by this ticket.

| Case | Product | m×n×k | ordinary f64 us | ordinary f32 us | speedup | native f64 us | native f32 us | speedup |
|---|---|---|---|---|---|---|---|---|
| si8 | gram | 40x20x2969 | 110.70 | 79.86 | 1.395 | 107.69 | 76.28 | 1.397 |
| si8 | rotation | 2969x20x40 | 76.60 | 37.96 | 1.976 | 73.44 | 37.10 | 1.976 |
| si8 | projection | 2969x2x40 | 11.88 | 5.99 | 1.994 | 11.52 | 5.97 | 1.941 |
| si16 | gram | 72x36x5915 | 590.98 | 381.65 | 1.547 | 591.40 | 380.51 | 1.554 |
| si16 | rotation | 5915x36x72 | 516.87 | 262.74 | 1.966 | 515.04 | 260.59 | 1.972 |
| si16 | projection | 5915x2x72 | 46.19 | 21.83 | 2.134 | 46.06 | 21.90 | 2.104 |
| si32 | gram | 136x68x11877 | 3874.30 | 2360.03 | 1.638 | 3798.28 | 2256.29 | 1.683 |
| si32 | rotation | 11877x68x136 | 3575.28 | 1857.12 | 1.919 | 3504.53 | 1829.13 | 1.915 |
| si32 | projection | 11877x2x136 | 171.44 | 97.36 | 1.725 | 168.52 | 96.40 | 1.751 |

## Accuracy and usable tolerance

For a normalized saved column x, let δ = H32(round32(x)) − H64(x). Report
`||δ||₂` (Ha), `||δ||₂ / ||H64 x||₂` (dimensionless), and the orthogonal
noise `||(I−xxᵀ)δ||₂` (Ha); its relative version divides by `||H64 x||₂`. The latter is the residual-noise floor after a
Rayleigh-value update, not a claim that a mixed-precision SCF has converged.
Every entry below is a maximum over the saved bands; maxima in different columns
need not belong to the same band. A near-zero eigenvalue can make the relative
metric large, so engine tolerances must be compared in Ha.

| Build | Snapshot | Bands | ||δ|| Ha | Relative residual floor | Orthogonal floor Ha | |xᵀδ| Ha | f64 vs production Ha |
|---|---|---|---|---|---|---|---|
| ordinary | si8-early | 20 | 1.254e-06 | 9.950e-06 | 1.220e-06 | 2.903e-07 | 1.546e-15 |
| ordinary | si8-converged | 20 | 6.717e-07 | 1.239e-05 | 6.670e-07 | 8.450e-08 | 1.356e-15 |
| ordinary | si16-early | 36 | 7.702e-07 | 7.629e-06 | 7.519e-07 | 2.061e-07 | 2.898e-15 |
| ordinary | si16-converged | 36 | 6.403e-07 | 1.090e-05 | 6.374e-07 | 1.554e-07 | 1.632e-15 |
| ordinary | si32-early | 68 | 6.887e-07 | 8.184e-06 | 6.756e-07 | 1.675e-07 | 1.850e-15 |
| ordinary | si32-converged | 68 | 1.134e-06 | 2.003e-05 | 1.133e-06 | 1.597e-07 | 1.911e-15 |
| native | si8-early | 20 | 1.254e-06 | 9.950e-06 | 1.220e-06 | 2.903e-07 | 1.546e-15 |
| native | si8-converged | 20 | 6.717e-07 | 1.239e-05 | 6.670e-07 | 8.450e-08 | 1.356e-15 |
| native | si16-early | 36 | 7.702e-07 | 7.629e-06 | 7.519e-07 | 2.061e-07 | 2.898e-15 |
| native | si16-converged | 36 | 6.403e-07 | 1.090e-05 | 6.374e-07 | 1.554e-07 | 1.632e-15 |
| native | si32-early | 68 | 6.887e-07 | 8.184e-06 | 6.756e-07 | 1.675e-07 | 1.850e-15 |
| native | si32-converged | 68 | 1.134e-06 | 2.003e-05 | 1.133e-06 | 1.597e-07 | 1.911e-15 |

The frozen projected Ritz experiment forms symmetrized `XᵀH32X` and
`XᵀH64X`, solves both small problems in fp64, and compares sorted roots.
Symmetrization is diagnostic, not a proposed repair for the production solver.
The unsymmetrized fp32 projected-H Hermitian defect is also retained.

| Build | Snapshot | Max Ritz error Ha | Max Hermitian defect Ha |
|---|---|---|---|
| ordinary | si8-early | 3.064e-07 | 1.597e-07 |
| ordinary | si8-converged | 8.856e-08 | 1.605e-07 |
| ordinary | si16-early | 2.292e-07 | 1.656e-07 |
| ordinary | si16-converged | 1.568e-07 | 1.064e-07 |
| ordinary | si32-early | 2.012e-07 | 1.764e-07 |
| ordinary | si32-converged | 1.705e-07 | 1.715e-07 |
| native | si8-early | 3.064e-07 | 1.597e-07 |
| native | si8-converged | 8.856e-08 | 1.605e-07 |
| native | si16-early | 2.292e-07 | 1.656e-07 |
| native | si16-converged | 1.568e-07 | 1.064e-07 |
| native | si32-early | 2.012e-07 | 1.764e-07 |
| native | si32-converged | 1.705e-07 | 1.715e-07 |

For the standard Ritz-change schedule, use **tau >= 1e-5 Ha** as the empirical
loose-step eligibility screen: it exceeds twenty times the largest measured
projected-root error (twenty allows two noisy successive Ritz values and a
tenfold margin). This is a screening margin, not a perturbation theorem or
a certified safe switch. At 1e-7 Ha the ordinary experiment already has
7/20 Si8 early, 6/36 Si16 early and 2/36 Si16 converged roots with larger
projected errors; at 1e-8 Ha all snapshot groups have failing roots.
A residual-based solver has a different restriction: the observed absolute
orthogonal floor is about 1e-6 Ha, so 1e-6 and tighter residual solves cannot
reliably use this H; a tenfold margin rounded upward gives 1e-4 Ha as a
conservative residual-screen threshold. Do not compare the dimensionless
relative error directly with a threshold in Ha.

These probes establish a noise floor on physical states, not convergence of
an fp32 Davidson recurrence. They do not bound error on every intermediate
correction direction or under near-degeneracy. Threshold sensitivity below
makes the recommendation auditable without selecting a threshold to meet 10%.

## Standard-history eligibility and wall ceiling

Current histories are three starts per system (uniform plus two deterministic
relative density perturbations <=5e-11); all nine standard runs converge.
The existing P14A counting collector is enabled at build time and independently
agrees exactly with every run’s sum of history Davidson steps and H vectors.
No standard first-step redo occurred in these histories. Timing includes one
production model construction and SCF, but excludes the extra source model
used solely by the test to initialize density transfer, snapshot extraction
and log output. Times from counting builds are observations, not acceptance.

| Case | Seed | SCF steps | Wall s | Eligible/all Davidson steps | Eligible/all H vectors | Eligible diag / wall % |
|---|---|---|---|---|---|---|
| si8 | 0 | 19 | 0.470027 | 38/82 | 444/923 | 30.34 |
| si8 | 1 | 24 | 0.552961 | 49/101 | 634/1136 | 32.38 |
| si8 | 2 | 23 | 0.532923 | 43/102 | 557/1076 | 29.39 |
| si16 | 0 | 17 | 1.680381 | 40/75 | 889/1558 | 37.16 |
| si16 | 1 | 15 | 1.610284 | 35/78 | 662/1448 | 31.53 |
| si16 | 2 | 18 | 1.822406 | 43/90 | 868/1641 | 33.87 |
| si32 | 0 | 33 | 16.259555 | 50/128 | 1800/5242 | 24.27 |
| si32 | 1 | 23 | 12.206765 | 32/104 | 1101/3910 | 21.07 |
| si32 | 2 | 29 | 15.316155 | 38/133 | 1371/4815 | 20.03 |

Previously recorded P31A real/standard, capacity=6 histories supply five
additional starts each for Si8 and Si16 (historical checkout, not pooled with
current timing). No equivalent Si32 history file was present at that path;
the current Si32 counting runs provide its evidence. Source paths and SHA-256
are in `P47B-evidence/recorded-manifest.json`.

| Historical case | Starts | Eligible step median % | Eligible vector median % |
|---|---|---|---|
| si8-gamma | 5 | 47.78 | 56.84 |
| si16-gamma | 5 | 45.33 | 47.42 |

We do not have H-region clocks for each production SCF call. Instead use
a deliberately **more generous ceiling**: `H_e <= D_e`, where D_e is the
entire diagonalization time in eligible iterations. It includes EVD,
orthogonalization, reductions, and other work that an fp32 H cannot accelerate.
For each case and build, use the **largest** median paired complete-H speedup
over early/converged snapshots and widths 1,2,4,8,b; then compute
`U = (D_e/W) * (1 - 1/s_max)`. Take the median over the three starts per
case, then the equally weighted median over Si8/Si16/Si32. This upper estimate
does not require the incorrect assumption that equal H vector counts imply
equal call costs. It is conditional on transferring the measured lane precision
ratio to production H and holding the fp64 trajectory fixed; it is not a
mathematical bound on an unimplemented fused f32 kernel or changed SCF trajectory.

| Build | Case | Eligible H vectors % | Eligible steps % | s_max | H-only saving % | D_e / wall % | Wall ceiling % |
|---|---|---|---|---|---|---|---|
| ordinary | si8 | 51.77 | 46.34 | 1.353 | 13.50 | 30.34 | 7.91 |
| ordinary | si16 | 52.89 | 47.78 | 1.250 | 10.58 | 33.87 | 6.78 |
| ordinary | si32 | 28.47 | 30.77 | 1.372 | 7.71 | 21.07 | 5.71 |
| native | si8 | 51.77 | 46.34 | 1.249 | 10.33 | 30.34 | 6.05 |
| native | si16 | 52.89 | 47.78 | 1.296 | 12.07 | 33.87 | 7.73 |
| native | si32 | 28.47 | 30.77 | 1.415 | 8.35 | 21.07 | 6.18 |

The H-only column is the fraction of total H work nominally removed under
equal cost/vector; it is **not** end-to-end saving. Treating all wall time as H
would ignore measured setup, density, mixing, audit and Davidson algebra costs.
The D_e ceiling already overcredits all eligible Davidson work to H.

| Eligibility threshold Ha | Ordinary median wall ceiling % | Native median wall ceiling % |
|---|---|---|
| 1e-04 | 5.79 | 5.60 |
| 1e-05 | 6.78 | 6.18 |
| 1e-06 | 8.97 | 9.10 |
| 1e-07 | 10.54 | 10.30 |


## Certification boundary and recommendation

Do **not** open an implementation ticket on this evidence. The empirical
standard-eligible threshold is 1e-5 Ha; neither release build reaches the frozen
10% median wall-saving rule even after charging all eligible diagonalization
work to H. Lowering eligibility to 1e-7 Ha to cross 10% would select a regime
where the measured projected-root errors already exceed the requested tolerance.
The unmodified production fused f64 H is also faster than the implemented f32
lane fallback on all six snapshot blocks in the ordinary build. A future generic
fused precision study could change the premise, but that speedup is not measured
here and must not be treated as established upside.

Strict certification requires a hard fp64 boundary, not just a tighter tolerance:

- Keep the full fixed-point certificate, H construction, energy/density audit,
  occupations, scalar reductions and independent residual evaluation fp64.
  A standard run remains uncertified regardless of the apparent residual.
- Discard and recompute **all** cached fp32 H images and reduced Hamiltonians
  before a fp64 solve/certificate. Casting an fp32 image to f64 does not restore
  its lost accuracy. Restarted/contracted Davidson bases can carry old images
  long after a nominal precision switch.
- Preserve fp64 Gram matrices, orthogonalization, eigenproblem solves, mixing
  state and convergence logic. The standalone Davidson GEMM speedups in this
  note do not authorize reduced precision in those operations.
- Independently recompute fp64 residuals on the current Hamiltonian; do not
  certify recurrence residuals assembled from rounded cached products. Ensure
  precision-specific buffers cannot alias certificate workspaces.
- Account for projector cancellation, Gamma pair separation, odd tails,
  Hermitian defects, degeneracies and correction-vector amplification. These
  Si snapshots cannot establish generality for skewed/slab cells, metals,
  different projector structures or near-zero gaps.

This is a rejection of an implementation recommendation under the frozen
measurement rule, not proof that mixed precision can never help. A revival needs
new measured precision ratios for the production fused structure (or a different
machine/workload), physical-state accuracy beyond this silicon set, explicit
fp64 cache invalidation and certificate isolation, and an end-to-end trajectory
comparison including conversion/refinement costs. dev08 evidence would still be
required for performance acceptance. No remote machine or QE calculation was run.

## Verification and limitations

The raw evidence retains failed draft compiles (wrong public field name and a
missing generic `Default` bound) and their corrected successful builds. The
f64 experimental H/production comparison is a numerical oracle, while the
projected Ritz test is a diagnostic, not a mixed-precision SCF acceptance test.
No total-energy equivalence claim is made for an unimplemented fp32 engine.

The host is macOS 27.0.1 arm64, rustc 1.99.1, LLVM 23.1.1, rustfft 6.4.1,
faer 0.24.4, source HEAD `0676a294da8405689a5497415179b4d6e1b38185`.
UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
CPU model and process listing were denied by the sandbox. Load averages and
all paired samples are retained; shared-host noise means these medians are
observations, not confidence bounds or remote acceptance. No local timing run
was overlapped with this task's builds. Cross-build timing comparisons do not
replace paired f32/f64 comparisons within each build.

Production files and `artifacts/*` remain unchanged. README generation is not
applicable: this ticket is measurement-only, owns no README/reference fixture
files and verifies no production implementation. All study tests are ignored
by the debug gates. Changes remain uncommitted on
`task/P47B-FP32-FEASIBILITY` for review.
