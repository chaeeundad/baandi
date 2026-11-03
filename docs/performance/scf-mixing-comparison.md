# P10D: optional density mixing implementation and measured comparison

Status: implemented and verified; defaults unchanged. The new strategy is
selectable in the library and CLI, with recorded parameters, raw-residual
contraction, history usage, resets and boundary-recovery diagnostics.
See [design and literature](../research/scf-mixing-strategies.md),
[ADR 0025](../decisions/0025-preconditioned-density-mixing.md) and
[CLI keys](../cli.md#optional-kerkerpulay-density-mixing-p10d).

## Recommendation and limits of the result

**Keep default Pulay.** The exact analytic dielectric model improves from
13/17/21 map evaluations to 3 with its known inverse preconditioner. In complete
SCFs, q0=0.5 Bohr^-1 reduces the synthetic fcc case from 24 to 21 iterations
and the synthetic elongated nonlocal case from 22 to 21; these modest reductions
do not establish general acceleration. Most other cases need more iterations.
The synthetic vacuum case is substantially worse, and higher q0 can fail the
fixed-point certificate within the budget. Si's four existing regression cases
are slower in iteration count with the tested Kerker parameters.

All strategies fail on the tested stretched Si supercell and Si vacuum slab:
these are reported failures, not metallic or material-convergence claims.
A wider-band supercell variant also fails within the budget. No thresholds,
certificate, gap policy or eigensolver controls were relaxed to accept a case.
Use the optional strategy for diagnosed long-wavelength sloshing and compare
raw residuals at unchanged tolerances; no universal q0 is recommended.

## Implementation scope and safety contracts

The new solver module composes an FFT residual filter with existing Pulay and
has a deterministic raw convex recovery for an infeasible first correction.
Existing `mixing.rs` is untouched. `lib.rs` changes only two module-export lines.
CLI `run.rs` changes only the strategy-resolution call and strategy/iteration
metadata. Parameter parsing and its tests are in config.rs.
SCF additions are the enum branch, FFT adapter construction, history introspection,
four diagnostic fields and their initialization, and extending the existing
feasibility match/reset counter. The controller and convergence arithmetic are
not rewritten. The adapter's two arrays and FFT plan are allocated once.
The existing Pulay history/eigensolve/result allocations remain a limitation;
this is not an allocation-free mixer rewrite or a performance optimization of
the default path. No unsafe code or dependencies were added.

## Environment and reproducibility

Starting commit `0142ebf`, uncommitted task changes, arm64 macOS,
rustc 1.99.1 (58cf8f9c3 2026-10-13), release profile, one engine thread.
The CPU-model query was denied (sysctl exit 1). The shared host runs other jobs;
three repetitions below characterize observed variability, not quiet-node speed.
Default-path measurements were taken before integration and after the final
recovery implementation. `sample PID 1` on the pre-change k444s run exited 255
because process examination was denied; the profiled CLI run itself exited 0.
Consequently no stack-profile conclusion or accepted engine speedup is asserted.
Claude should repeat full-SCF timing/profiling on quiet dev08 before accepting
any performance claim. The task excludes baseline and optimization-log edits;
this document retains the algorithm experiment, rejected behavior and timings.

Validation UPF SHA-256: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`. It was read in place, never
copied into the repository. Its provenance remains unverified as documented in
the existing reference bundles. Final CLI binary SHA-256:
`4bebfdabf7e5ac495fcede09d81ff6122263a952965ebb710f1f6b187547a3c4`.

## Complete SCF systems and measured gaps

The ignored cross-layer test `crates/dft-solvers/tests/mixing_scf.rs` contains the
complete input generator. It invokes the actual CLI without introducing an
engine dependency into solvers. Temporary Gaussian-ion UPFs use the formulas
from `dft-engine/tests/common/mod.rs`, evaluated by Python math.erf, with 1201
logarithmic radial samples, rc=1 Bohr. Projector cases have two s and one p
projector with off-diagonal s coupling; local cases have none. They are numerical
fixtures, not real-element pseudopotentials. Shapes include cubic, fcc, skew,
elongated and slab-like vacuum cells; 1,2,4,8 atoms and 5,8,12,20 bands occur.
Gamma and shifted 2x1x1, 2x2x1, 2x2x2 and 4x4x4 meshes are tested. The even skew
case explicitly uses 24x24x32 FFT samples; other cases use the existing grid rule.

Synthetic cases use 10 Ry, beta=0.3 and 150 iterations. Existing Si regression
cases use the benchmark inputs unchanged (30 Ry, beta=0.5, 100 iterations).
The Si supercell is 1x1x2 primitive replicated along the third lattice vector,
then stretched by 1.15, at 20 Ry, beta=0.3 and 150 iterations, with 12 and 20
bands separately. The Si slab is an 8-atom conventional cell plus two cell
lengths of vacuum in z, at 12 Ry, beta=0.3 and 150 iterations. These two cutoffs
are characterization settings, not converged material energies. All use history
6 and unchanged density 1e-7 Bohr^-3/2, energy 1e-8 Ha and eigen residual 1e-8 Ha
controls. q0 values are 0.5,1,2 Bohr^-1; no fitting by species is performed.

| System | Bands | Sampled insulating gap of default converged result (Ha) |
| --- | ---: | ---: |
| synthetic-cubic-z2-p0 | 5 | 3.0138379331e-01 |
| synthetic-skew-z2-p3-k211 | 8 | 2.1433863096e-02 |
| synthetic-slab-z2-p0 | 8 | 2.3076481889e-01 |
| synthetic-fcc-z4-p3 | 8 | 7.9480698311e-03 |
| synthetic-even-skew-z2-p0-k221 | 12 | 2.4967981241e-01 |
| si-gamma | 8 | 7.8512231322e-02 |
| si-k222s | 8 | 7.9948786686e-02 |
| si-k222 | 8 | 1.7213826483e-02 |
| si-k444s | 8 | 4.1614159566e-02 |

These sampled gaps span small to large gaps in the tested discretizations; they
are not a proof of a gap throughout continuous reciprocal space. Failed Si
supercell/slab states are not labeled insulating or metallic.

## SCF convergence matrix

Numbers are evaluated SCF iterations including the final certificate iteration.
`N/150` means NotConverged at the budget; `D` means explicit residual-growth
Diverged. Every row/strategy was repeated three times with identical exits;
where results were returned, energy, density hash, counts and complete numerical
iteration histories (excluding timings) were identical across repeats.

| Case | Default Pulay | Kerker q0=0.5 | Kerker q0=1 | Kerker q0=2 |
| --- | ---: | ---: | ---: | ---: |
| synthetic-cubic-z2-p0 | 14 | 15 | 17 | 37 |
| synthetic-skew-z2-p3-k211 | 22 | 21 | 30 | 57 |
| synthetic-slab-z2-p0 | 15 | 51 | N/150 | N/150 |
| synthetic-fcc-z4-p3 | 24 | 21 | 24 | 28 |
| synthetic-even-skew-z2-p0-k221 | 13 | 14 | 22 | 57 |
| si-gamma | 10 | 12 | 14 | 17 |
| si-k222s | 10 | 11 | 13 | 19 |
| si-k222 | 10 | 11 | 13 | 16 |
| si-k444s | 10 | 11 | 13 | 18 |
| si-supercell-stretch | D | D | D | D |
| si-slab-vacuum | N/150 | N/150 | N/150 | N/150 |
| si-supercell-stretch-wide | N/150 | N/150 | N/150 | N/150 |

For successful paired SCFs the maximum absolute total-energy difference from
same-case default Pulay is **1.332267629550e-13 Ha**. The maximum input/output
charge defect across all returned histories, including budget failures, is
**1.771027768882e-12 electrons** (summation roundoff). Each returned density minimum
was >=-1e-14 electrons/Bohr^3. The harness enforces charge error <1e-10 electrons,
the unchanged physical convergence gates and the independent certificate;
its latest version also asserts successful paired energies differ <=1e-12 Ha.

The synthetic slab at q0=1 is a useful counterexample to energy-only convergence:
after 150 iterations its raw residual is 9.6191408141e-9 Bohr^-3/2 and dE is
-1.4210854715e-13 Ha, but the certificate's occupied eigenvalue change is
2.2653225645e-7 Ha, exceeding the unchanged 1e-7 Ha threshold. It correctly
returns NotConverged. No failed trajectory contributes a converged-energy claim.

The 12-band stretched Si runs reach small raw residuals then grow to about
0.014 Bohr^-3/2 and trigger the existing three-step growth detector. They fail
in 15/15/16/19 iterations for default/q0=0.5/1/2. Increasing to 20 bands avoids
that particular early growth error but all four runs exhaust 150 iterations.
This is not evidence of a successful material SCF or a diagnosed eigensolver
root cause. The Si slab similarly exhausts 150 iterations for all strategies.

### Three-run complete-SCF timings and final numerical values

CSV below is a retained, machine-readable summary of real measured results.
Time is CLI wall time in seconds (includes setup and output), with min/median/max
from three serial repetitions of each case. Different cases were measured on a
shared host; small differences cannot establish a speedup. Empty energy/residual
fields indicate an error without a final state, not zero values. Failed final
energies are informational only. Reset and recovery counts refer to the first
repeat (identical across repeats).

```csv
case,strategy,exit,iterations,wall_min_s,wall_median_s,wall_max_s,energy_ha,raw_residual_bohr_m3_over_2,max_charge_error_electrons,resets,raw_recoveries
synthetic-cubic-z2-p0,pulay,0,14,0.077774,0.080547,0.373518,-1.3969438183603748e+00,1.3723820442608110e-08,2.0872192862952943e-14,1,0
synthetic-cubic-z2-p0,kerker-0.5,0,15,0.089930,0.092093,0.112345,-1.3969438183603771e+00,3.4649415740205291e-09,1.2656542480726785e-14,1,0
synthetic-cubic-z2-p0,kerker-1,0,17,0.092469,0.093257,0.099849,-1.3969438183603513e+00,5.0279087462904499e-08,2.3536728122053319e-14,1,0
synthetic-cubic-z2-p0,kerker-2,0,37,0.130767,0.147640,0.152325,-1.3969438183603575e+00,4.2514512840794752e-08,3.2862601528904634e-14,3,0
synthetic-skew-z2-p3-k211,pulay,0,22,0.201147,0.248998,0.254351,-7.8702655543461109e-01,1.3127926025677674e-08,8.6153306710912148e-14,0,0
synthetic-skew-z2-p3-k211,kerker-0.5,0,21,0.320768,0.386800,0.591063,-7.8702655543458311e-01,2.7043392393350367e-08,8.6153306710912148e-14,3,3
synthetic-skew-z2-p3-k211,kerker-1,0,30,0.261949,0.309323,0.401803,-7.8702655543459432e-01,1.7745764779119751e-08,8.6153306710912148e-14,5,4
synthetic-skew-z2-p3-k211,kerker-2,0,57,0.356279,0.388182,0.391581,-7.8702655543460809e-01,2.5819154722097987e-08,8.6153306710912148e-14,9,9
synthetic-slab-z2-p0,pulay,0,15,0.157590,0.160087,0.166860,-2.8432473120033350e+00,3.1015693763477976e-09,6.7501559897209518e-14,2,0
synthetic-slab-z2-p0,kerker-0.5,0,51,0.365690,0.372071,0.500826,-2.8432473120033670e+00,2.0209054643006413e-09,6.7501559897209518e-14,16,12
synthetic-slab-z2-p0,kerker-1,3,150,0.957255,1.057821,1.116271,-2.8432473120032249e+00,9.6191408140615445e-09,6.7501559897209518e-14,62,56
synthetic-slab-z2-p0,kerker-2,3,150,1.033672,1.143176,1.674721,-2.8419732886612028e+00,3.1013432914125702e-03,6.7501559897209518e-14,73,51
synthetic-fcc-z4-p3,pulay,0,24,0.082337,0.085891,0.093833,-4.6845611842778565e+00,7.3317694355227850e-09,7.1054273576010019e-14,0,0
synthetic-fcc-z4-p3,kerker-0.5,0,21,0.068642,0.074290,0.077952,-4.6845611842778023e+00,3.9342173570756040e-08,7.1054273576010019e-14,0,0
synthetic-fcc-z4-p3,kerker-1,0,24,0.078199,0.082370,0.085167,-4.6845611842778503e+00,1.5646578625203987e-08,7.1054273576010019e-14,0,0
synthetic-fcc-z4-p3,kerker-2,0,28,0.080147,0.081273,0.087789,-4.6845611842778530e+00,1.4591049789179721e-08,7.1054273576010019e-14,0,0
synthetic-even-skew-z2-p0-k221,pulay,0,13,0.216309,0.232826,0.241150,-1.3669704386157968e+00,2.6495161804603902e-09,1.1324274851176597e-13,1,0
synthetic-even-skew-z2-p0-k221,kerker-0.5,0,14,0.248048,0.304952,0.315800,-1.3669704386158004e+00,3.4313997227903851e-08,1.1324274851176597e-13,1,0
synthetic-even-skew-z2-p0-k221,kerker-1,0,22,0.378717,0.426772,0.449168,-1.3669704386157948e+00,1.4667459683687663e-08,1.1324274851176597e-13,3,1
synthetic-even-skew-z2-p0-k221,kerker-2,0,57,0.776395,0.879679,1.004501,-1.3669704386157955e+00,2.2955177793139916e-08,1.1324274851176597e-13,15,7
si-gamma,pulay,0,10,0.151977,0.160164,0.170967,-7.2955304682107718e+00,4.1634599768829370e-08,1.5631940186722204e-13,1,0
si-gamma,kerker-0.5,0,12,0.171203,0.171850,0.172197,-7.2955304682108935e+00,5.0929702002991409e-09,1.5631940186722204e-13,2,1
si-gamma,kerker-1,0,14,0.187511,0.187833,0.189609,-7.2955304682108828e+00,3.3936959630700281e-09,1.5631940186722204e-13,2,2
si-gamma,kerker-2,0,17,0.211109,0.212601,0.260246,-7.2955304682109050e+00,2.3710423238840313e-08,1.5631940186722204e-13,2,2
si-k222s,pulay,0,10,0.554005,0.620753,0.626563,-7.9194234124690670e+00,6.8382774262092927e-09,1.5631940186722204e-13,0,0
si-k222s,kerker-0.5,0,11,0.557917,0.580029,0.594083,-7.9194234124690714e+00,1.5264316877273961e-09,1.5631940186722204e-13,0,0
si-k222s,kerker-1,0,13,0.639072,0.647494,0.697892,-7.9194234124690874e+00,4.9074373688657822e-09,1.5631940186722204e-13,1,1
si-k222s,kerker-2,0,19,0.968693,1.025188,1.035727,-7.9194234124690661e+00,2.0854651462593440e-08,1.5631940186722204e-13,2,2
si-k222,pulay,0,10,0.967604,0.978164,0.982982,-7.8312300571211448e+00,5.2367282940062214e-09,1.5631940186722204e-13,0,0
si-k222,kerker-0.5,0,11,1.019456,1.024529,1.031128,-7.8312300571211431e+00,1.6041061361055807e-09,1.5631940186722204e-13,1,0
si-k222,kerker-1,0,13,1.138382,1.165852,1.193972,-7.8312300571211262e+00,5.0786492330533508e-10,1.5631940186722204e-13,1,1
si-k222,kerker-2,0,16,1.410974,1.414303,1.460555,-7.8312300571210596e+00,3.9564135492572505e-08,1.5631940186722204e-13,2,2
si-k444s,pulay,0,10,4.999052,6.101537,6.975186,-7.9259761266480808e+00,5.5647399177958846e-09,1.5631940186722204e-13,0,0
si-k444s,kerker-0.5,0,11,5.165012,5.326394,6.419094,-7.9259761266480755e+00,1.4574765148661685e-09,1.5631940186722204e-13,0,0
si-k444s,kerker-1,0,13,5.690457,6.159338,6.773117,-7.9259761266480648e+00,2.6876617270903912e-09,1.5631940186722204e-13,1,1
si-k444s,kerker-2,0,18,6.917449,7.771365,9.837779,-7.9259761266480613e+00,3.3026632244558751e-08,1.5631940186722204e-13,2,2
si-supercell-stretch,pulay,4,,0.381456,0.382244,0.393783,,,,,
si-supercell-stretch,kerker-0.5,4,,0.343695,0.346445,0.351222,,,,,
si-supercell-stretch,kerker-1,4,,0.361653,0.366652,0.372306,,,,,
si-supercell-stretch,kerker-2,4,,0.389692,0.389738,0.395058,,,,,
si-slab-vacuum,pulay,3,150,7.247047,7.281945,7.678691,-3.1228401698483410e+01,1.8536300861782054e-03,9.3081098384573124e-13,3,0
si-slab-vacuum,kerker-0.5,3,150,9.828073,10.356783,10.663548,-3.1231150327838051e+01,1.7884116263108935e-03,9.3081098384573124e-13,19,15
si-slab-vacuum,kerker-1,3,150,9.469166,9.730529,11.031577,-3.1227571983520889e+01,1.9355505364752031e-03,9.3081098384573124e-13,28,25
si-slab-vacuum,kerker-2,3,150,11.604533,11.761604,13.769068,-3.1225453522278272e+01,2.2073777399339955e-03,9.3081098384573124e-13,44,45
si-supercell-stretch-wide,pulay,3,150,4.051136,4.109828,4.258191,-1.5095183439815381e+01,5.3480485881486879e-03,1.0427214647279470e-12,0,0
si-supercell-stretch-wide,kerker-0.5,3,150,3.191844,3.465869,5.813593,-1.5091817184941597e+01,5.8729263976185853e-03,1.0427214647279470e-12,1,0
si-supercell-stretch-wide,kerker-1,3,150,5.032803,6.246808,9.359231,-1.5101214672029728e+01,3.8375677007263000e-03,1.0427214647279470e-12,1,1
si-supercell-stretch-wide,kerker-2,3,150,3.828402,3.985103,4.078238,-1.5096579292992669e+01,5.2469901404085788e-03,1.0427214647279470e-12,8,6
```

### Residual and energy history comparisons

Selected samples of complete recorded trajectories are below. Total energies
are in Ha and raw SCF density residuals in Bohr^-3/2. Full per-iteration raw
histories, certificates, diagnostics, parameters, input TOMLs, stdout/stderr
and density files remain in `/private/tmp/p10d-final/` and `/private/tmp/p10d-wide/`;
they are reproducible using the retained test driver. The JSON record is
`scf-matrix.json` in each directory. Temporary evidence is not committed.

| Case / strategy | Iteration | Raw residual | Total energy (Ha) |
| --- | ---: | ---: | ---: |
| synthetic-fcc-z4-p3 / pulay | 1 | 7.479408767e-02 | -4.267639160411770 |
| synthetic-fcc-z4-p3 / pulay | 2 | 4.690448129e-02 | -4.388322187532433 |
| synthetic-fcc-z4-p3 / pulay | 5 | 3.074975402e-02 | -4.678211218547236 |
| synthetic-fcc-z4-p3 / pulay | 10 | 5.558595109e-03 | -4.684229560119530 |
| synthetic-fcc-z4-p3 / pulay | 24 | 7.331769436e-09 | -4.684561184277857 |
| synthetic-fcc-z4-p3 / kerker-0.5 | 1 | 7.479408767e-02 | -4.267639160411770 |
| synthetic-fcc-z4-p3 / kerker-0.5 | 2 | 5.012374404e-02 | -4.370448422105659 |
| synthetic-fcc-z4-p3 / kerker-0.5 | 5 | 3.142729418e-02 | -4.679066757083177 |
| synthetic-fcc-z4-p3 / kerker-0.5 | 10 | 7.005151531e-03 | -4.683960982669417 |
| synthetic-fcc-z4-p3 / kerker-0.5 | 21 | 3.934217357e-08 | -4.684561184277802 |
| synthetic-skew-z2-p3-k211 / kerker-0.5 | 1 | 6.280160106e-02 | -0.410047947595276 |
| synthetic-skew-z2-p3-k211 / kerker-0.5 | 2 | 4.790498557e-02 | -0.478799928049573 |
| synthetic-skew-z2-p3-k211 / kerker-0.5 | 5 | 1.809033174e-02 | -0.720894820383247 |
| synthetic-skew-z2-p3-k211 / kerker-0.5 | 10 | 1.902021989e-03 | -0.787004101738887 |
| synthetic-skew-z2-p3-k211 / kerker-0.5 | 21 | 2.704339239e-08 | -0.787026555434583 |
| si-gamma / pulay | 1 | 6.826656818e-02 | -7.172342783970571 |
| si-gamma / pulay | 2 | 2.760641616e-02 | -7.271614631586712 |
| si-gamma / pulay | 5 | 1.382313393e-04 | -7.295529729291241 |
| si-gamma / pulay | 10 | 4.163459977e-08 | -7.295530468210772 |
| si-gamma / kerker-1 | 1 | 6.826656818e-02 | -7.172342783970571 |
| si-gamma / kerker-1 | 2 | 4.397966638e-02 | -7.229410852002419 |
| si-gamma / kerker-1 | 5 | 2.897387215e-03 | -7.295475814556433 |
| si-gamma / kerker-1 | 10 | 3.409604820e-06 | -7.295530467936500 |
| si-gamma / kerker-1 | 14 | 3.393695963e-09 | -7.295530468210883 |
| synthetic-slab-z2-p0 / kerker-1 | 1 | 1.123913245e-01 | -2.533064634405786 |
| synthetic-slab-z2-p0 / kerker-1 | 2 | 9.822219296e-02 | -2.524759825371158 |
| synthetic-slab-z2-p0 / kerker-1 | 5 | 2.813023850e-02 | -2.777907534641514 |
| synthetic-slab-z2-p0 / kerker-1 | 10 | 1.622913902e-02 | -2.822366644429525 |
| synthetic-slab-z2-p0 / kerker-1 | 150 | 9.619140814e-09 | -2.843247312003225 |
| si-slab-vacuum / kerker-1 | 1 | 4.998991055e-02 | -7.227269630895378 |
| si-slab-vacuum / kerker-1 | 2 | 4.361158797e-02 | -9.133497072349755 |
| si-slab-vacuum / kerker-1 | 5 | 2.650221504e-02 | -19.600143724131733 |
| si-slab-vacuum / kerker-1 | 10 | 2.127002859e-02 | -28.940146739159488 |
| si-slab-vacuum / kerker-1 | 150 | 1.935550536e-03 | -31.227571983520889 |

## Parameter variation and rejected behavior

Additional three-repeat runs at beta=0.1 and 0.8 on the synthetic fcc and Si
Gamma cases preserve all gates. Counts below use the same q0 grid and history=6.
All listed parameter-variation runs converge and repeat identically; this is a
limited stability result for these systems, not a guaranteed safe beta range.

| Case | Beta | Pulay | q0=0.5 | q0=1 | q0=2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| synthetic-fcc-z4-p3 | 0.1 | 32 | 22 | 42 | 61 |
| si-gamma | 0.1 | 15 | 17 | 20 | 28 |
| synthetic-fcc-z4-p3 | 0.8 | 20 | 19 | 22 | 26 |
| si-gamma | 0.8 | 9 | 11 | 11 | 13 |

Before recovery was implemented, filtered Pulay plus the unchanged controller's
feasibility reset became trapped on a near-zero density boundary in Si Gamma,
the synthetic elongated nonlocal cell and the synthetic slab. At q0=1, Si
stalled at 0.0293218823 Bohr^-3/2 with 98 resets in 100 iterations; synthetic
elongated/slab cases exhausted 150 iterations. Failure was across several
systems, not Si alone. The initial strict Si characterization therefore failed
(exit 1); the initial synthetic test also failed (cargo exit 101) because its
negative-error assertion looked only at the truncated tail of long stderr.
Both are retained as development failures, not passing evidence.

The accepted deterministic recovery uses raw convex linear mixing with
min(beta,1) only when a negative filtered proposal has no history differences.
It uses no tuned density threshold or per-material recognition. This avoids
repeating an outward filtered direction, preserves charge without clipping,
and is recorded as `mixing_raw_recovery`. With this change, the synthetic
elongated q0=0.5 case converges in 21 iterations and Si Gamma in 12. Synthetic
slab q0=1/2 still fails in the allowed budget and is not hidden.
Full failed pre-recovery records remain in `/private/tmp/p10d-synthetic/` and
`/private/tmp/p10d-real/`. The research note explains the feasibility argument.

## Failure detection and raw-residual certification

- q0=1e100 Bohr^-1 suppresses the filtered update in Si Gamma, but the raw
  residual stays 0.0682755072 Bohr^-3/2; a three-iteration budget correctly returns
  exit 3 / NotConverged. No filtered norm can bypass the convergence gates.
- Linear beta=100 returns exit 4 / Diverged explicitly.
- Kerker–Pulay beta=1e308 on the synthetic fcc case returns exit 4 / Diverged,
  reporting nonfinite arithmetic, not silently reducing parameters or declaring
  convergence. The unchanged controller captures mixer failure in its history.
- Invalid q0, beta, zero history, unknown strategy and inapplicable keys fail
  config tests; expression-valued parameters are tested. Solver tests also
  cover transform failure, restart, mean change and nearly stationary maps.

## Analytic response model

Grid 5x3x5, 75 samples, uniform target density 1 electron/Bohr^3. Initial
zero-mean perturbation is 1e-4*sin(17*i+3) with its sampled mean removed.
B rows are [[0.7,0.12,0],[0,0.9,0.13],[0,0,0.6/stretch]] Bohr^-1;
stretch=1,4,12 changes the skew reciprocal metric and corresponding cell.
These are synthetic dielectric response problems, not pseudopotential cells.
The physical-cell volume is not used in the harness norm or charge metric.

Each nonzero mode obeys R_G=-(1+1/G^2)*(x_G-target_G), evaluated from independently
constructed analytic response coefficients and a direct DFT. The unaccelerated
beta=0.3 map can amplify long wavelengths. At q0=1 Bohr^-1, Kerker exactly
inverts this model dielectric on nonzero modes. All strategies use beta=0.3,
history=6, the same start, and stop on raw Euclidean residual <1e-10
in electron/Bohr^3 sample units. The count includes the initial and final map
evaluations; 3 evaluations mean 2 mixed updates. No energy is defined for this
fixed-point harness. Its success is not an SCF fixed-point certificate.

| Stretch | Strategy | Map evaluations | Time (s), one release run | Max sum-density defect (unit voxel: electrons) | Final raw residual |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | default | 13 | 0.000834 | 1.421e-14 | 3.262716e-11 |
| 1 | q0=0.5 | 11 | 0.001015 | 2.842e-14 | 1.136951e-11 |
| 1 | q0=1 | 3 | 0.000234 | 0.000e0 | 2.025967e-15 |
| 1 | q0=2 | 14 | 0.001359 | 2.842e-14 | 7.553645e-11 |
| 4 | default | 17 | 0.000925 | 1.421e-14 | 6.953316e-11 |
| 4 | q0=0.5 | 12 | 0.001339 | 2.842e-14 | 3.575911e-11 |
| 4 | q0=1 | 3 | 0.000278 | 0.000e0 | 3.557914e-15 |
| 4 | q0=2 | 15 | 0.002114 | 2.842e-14 | 4.963576e-11 |
| 12 | default | 21 | 0.001597 | 4.263e-14 | 7.425447e-11 |
| 12 | q0=0.5 | 13 | 0.001277 | 2.842e-14 | 3.905271e-11 |
| 12 | q0=1 | 3 | 0.000232 | 0.000e0 | 3.126416e-15 |
| 12 | q0=2 | 21 | 0.002673 | 7.105e-14 | 2.592747e-11 |

For the exact model preconditioner all three cases take 3 evaluations versus
13,17,21 for default Pulay. q0=0.5 and 2 remain convergent here, but q0=2 has
no iteration advantage at stretch=12. This supports an optional response-based
strategy, not a universal q0/default recommendation. The worst measured charge
sum defect across these runs is 7.105e-14 electrons for unit voxel volume.
The actual Fourier zero-mode multiplier is exactly zero; real-space charge
is conserved to summation roundoff, not claimed bitwise exact.

Additional tests cover independent cosine-mode attenuation on odd/even skew
grids, conjugate Nyquist weights, invalid q0/cell/grid/densities, changed-mean
and explicit restart, contraction diagnostics, injected transform failure,
and nearly stationary maps with derivative 0.999 and beta=0.1,0.3,0.8.
The latter require final pointwise error <2e-8 in density units and charge-sum
error <2e-14 electrons. They do not establish robustness to a general nonlinear
map or replace SCF divergence tests. The original four debug tests finished in 0.30 s total; the final suite adds
a fifth explicit boundary-recovery test, and all five pass in debug/release.

### Complete release residual traces

Each line records raw residual norms at every map evaluation, in the same
sample-density units as above. This retained evidence is generated from actual
test output; no timing or numerical result is inferred from a source paper.

```text
stretch=1 default: [0.0007888384169352465, 0.00048206932666466137, 3.0131630561394637e-5, 1.490959400325882e-5, 2.5998053202205406e-6, 5.958744637947082e-7, 1.9802617712735933e-7, 3.039656132157059e-8, 1.4078779975365462e-8, 1.8021701479146424e-9, 3.201644540156861e-10, 1.0840456020167026e-10, 3.262715577641586e-11]
stretch=1 q0=0.5: [0.0007888384169352465, 0.00050334159463228, 2.2718999001509705e-5, 4.736704619917562e-6, 6.991505043391476e-7, 1.0164120313232879e-7, 1.7991479902193218e-8, 2.335258516872764e-9, 1.2860719369510377e-9, 1.132537128569847e-10, 1.1369505934682133e-11]
stretch=1 q0=1: [0.0007888384169352465, 0.0005521868918546471, 2.025966904339278e-15]
stretch=1 q0=2: [0.0007888384169352465, 0.0006470276411431983, 5.495711808969515e-5, 8.102765671248046e-6, 2.430985894092113e-6, 5.070530632844439e-7, 4.4372320882181825e-7, 9.315597330212645e-8, 7.962939280029003e-8, 1.1101353150202714e-8, 1.449434020527091e-9, 3.1631748831137616e-10, 2.715085410596159e-10, 7.553645144778366e-11]
stretch=4 default: [0.0008537425459847741, 0.0012110302357858929, 0.0014766001009823605, 0.00019678913562786045, 2.878584656541245e-5, 9.674225142971467e-6, 1.372316353608032e-6, 2.3827573489866865e-7, 9.719909673298429e-8, 2.581994625528102e-8, 8.251648353911814e-8, 7.100342969160453e-9, 5.2013384401354444e-9, 7.0046917675189135e-9, 3.929044498738942e-9, 3.1499560002904794e-10, 6.953316033763185e-11]
stretch=4 q0=0.5: [0.0008537425459847741, 0.0005209649102737615, 3.889004937530742e-5, 1.571846513628533e-5, 3.7894591194465197e-6, 1.0165464556618044e-6, 1.2514983382858265e-7, 1.722034000663281e-8, 9.593203365590968e-9, 1.270988061389946e-9, 1.493907946641659e-10, 3.5759111346320856e-11]
stretch=4 q0=1: [0.0008537425459847741, 0.0005976197821894365, 3.5579140431534355e-15]
stretch=4 q0=2: [0.0008537425459847741, 0.0007147066957853837, 0.00010660677132467111, 4.500202309390465e-5, 1.808128039543983e-5, 5.922441191474513e-6, 5.349627029293809e-6, 1.1573911975425776e-6, 2.0237229430829992e-7, 1.762137821476116e-7, 3.217671496100159e-8, 7.455783955117573e-9, 8.67494923381224e-10, 7.349971679963302e-10, 4.9635758172911673e-11]
stretch=12 default: [0.0009057079836752709, 0.00473779173102185, 0.0010726243563422987, 0.0011744036162507773, 3.3818636087555526e-5, 1.390560505889264e-5, 2.540652097457114e-6, 4.571445194690054e-7, 1.1343845125073053e-7, 2.3137140134121186e-7, 2.5806762654759367e-7, 5.333121566750958e-7, 3.3347861650640054e-8, 2.155405061126663e-8, 2.880269440792048e-9, 9.615693699100575e-10, 2.8780031118136567e-10, 1.7788969880405077e-10, 1.088917001901046e-9, 5.70646747590037e-10, 7.425447429203609e-11]
stretch=12 q0=0.5: [0.0009057079836752709, 0.000524010599310338, 9.032421464754476e-5, 2.376629943564602e-5, 6.220295322937743e-6, 1.6331217559780735e-6, 1.88602997751268e-7, 2.7259909049338625e-8, 1.342724103561707e-8, 2.3298075219813932e-9, 3.015474462653858e-10, 1.5921510563086556e-10, 3.9052708005030185e-11]
stretch=12 q0=1: [0.0009057079836752709, 0.0006339955885728598, 3.1264164235997256e-15]
stretch=12 q0=2: [0.0009057079836752709, 0.0007674071516340834, 0.0001837381542654623, 8.934616767419142e-5, 3.5939047031634707e-5, 1.1713608622130954e-5, 1.0574149181146383e-5, 2.2658817228150636e-6, 3.662711586931621e-7, 3.1750738778311773e-7, 5.350516231778486e-8, 1.5512118736739647e-8, 3.43350530199374e-9, 3.018967233113267e-9, 3.3951715246961995e-10, 2.7190231594794116e-10, 2.814340177901386e-10, 2.2960284915657627e-10, 1.1881621788804926e-10, 1.0571734382870002e-10, 2.592746736535847e-11]
```


### Three serial release repetitions

All three invocations of the solver test command above, adding
`--test-threads=1`, exited 0 (4 tests each). Counts, charge defects and complete
residual traces were identical across repetitions. Wall durations below include
the direct-DFT harness; no production performance inference is warranted.

| Stretch | Strategy | Median time (s) | Min–max (s) |
| ---: | --- | ---: | --- |
| 1 | q0=0.5 | 0.000996 | 0.000995–0.001098 |
| 1 | q0=1 | 0.000233 | 0.000230–0.000242 |
| 1 | q0=2 | 0.001275 | 0.001274–0.001279 |
| 4 | default | 0.000846 | 0.000827–0.000867 |
| 4 | q0=0.5 | 0.001127 | 0.001096–0.001382 |
| 4 | q0=1 | 0.000240 | 0.000234–0.000286 |
| 4 | q0=2 | 0.001551 | 0.001374–0.001602 |
| 12 | default | 0.001120 | 0.001016–0.001161 |
| 12 | q0=0.5 | 0.001302 | 0.001190–0.001362 |
| 12 | q0=1 | 0.000256 | 0.000231–0.000267 |
| 12 | q0=2 | 0.002063 | 0.001960–0.002125 |

## Default Si equivalence and stored QE comparisons

`bench_scf.py --require-identical` against `baseline-after-o2.json` passed for
Gamma, k222s, k222 and k444s both before and after final integration. All remain
10 iterations, bit-identical density hashes and zero energy difference (Ha).
A three-repeat integrated run and a final single-repeat run passed; timing
ratios printed by the benchmark are not accepted speedup claims on this host.

Existing real QE 7.6 output bundles were compared to fresh default-mixer outputs
without changing tolerances. These are **comparisons to stored runs**, not new
QE executions. Commands/version/input/UPF hashes remain in the original bundle
provenance. No pw.x was installed or run locally and no QE source was read.

| Stored reference | Comparison | Measured difference | Exit |
| --- | --- | ---: | ---: |
| si30_ours | Maximum eigenvalue difference | 1.283668146834e-6 Ha | 0 |
| si30_ours | Maximum compared energy-term difference (one-electron) | 2.594966043379e-6 Ha | 0 |
| k222_c | Total energy per atom / maximum eigenvalue | 1.281439427459e-6 Ha/atom / 1.990419700128e-6 Ha | 0 |
| k222s_c | Total energy per atom / maximum eigenvalue | 1.281265466613e-6 Ha/atom / 1.772334107711e-6 Ha | 0 |
| k444s_c | Total energy per atom / maximum eigenvalue | 1.281675959586e-6 Ha/atom / 2.100133025743e-6 Ha | 0 |

## Commands, exits and reproduction

Every command below was executed with exit 0. The full physical matrix ran
11 cases x 4 strategies x 3 repetitions plus negative tests (265.60 s); the
20-band supercell and beta variation sweeps were separate runs. A final source
revision added explicit cross-strategy energy assertions and the new-strategy
nonfinite test; the final targeted release command executes those assertions
on cubic, fcc and Si Gamma. The retained complete matrix also satisfies the
same energy bound (measured above), and numerical histories were independently
compared across all three repetitions with timing fields removed.

```sh
cargo build --release --offline -p dft-cli
cargo test --offline -p dft-solvers -p dft-cli
CARGO_NET_OFFLINE=true scripts/check.sh
cargo test --release --offline -p dft-engine -- --ignored
DFT_MIXING_REPEATS=3 DFT_MIXING_EVIDENCE=/private/tmp/p10d-final \
  cargo test --release --offline -p dft-solvers --test mixing_scf -- --ignored --nocapture
DFT_MIXING_REPEATS=1 DFT_MIXING_CASES=synthetic-cubic-z2-p0,synthetic-fcc-z4-p3,si-gamma \
  DFT_MIXING_EVIDENCE=/private/tmp/p10d-acceptance \
  cargo test --release --offline -p dft-solvers --test mixing_kerker --test mixing_scf \
  -- --include-ignored --nocapture
python3 -I scripts/bench_scf.py --binary target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222,k222s,k444s --repeats 1 --threads 1 \
  --baseline docs/performance/baseline-after-o2.json --require-identical \
  --record /private/tmp/p10d-final-default.json
python3 -I scripts/verify_qe.py \
  --ours /private/tmp/p10d-final/si-gamma-pulay-0/out/results.json \
  --bundle tests/reference/si30_ours \
  --json /private/tmp/p10d-qe-gamma.json --markdown /private/tmp/p10d-qe-gamma.md
python3 -I scripts/compare_kpoints_qe.py \
  --ours /private/tmp/p10d-final/si-k222-pulay-0/out/results.json \
  --qe-out tests/reference/kpoints/k222_c/qe.out
python3 -I scripts/compare_kpoints_qe.py \
  --ours /private/tmp/p10d-final/si-k222s-pulay-0/out/results.json \
  --qe-out tests/reference/kpoints/k222s_c/qe.out
python3 -I scripts/compare_kpoints_qe.py \
  --ours /private/tmp/p10d-final/si-k444s-pulay-0/out/results.json \
  --qe-summary tests/reference/kpoints/k444s_c/qe_summary.json
python3 -I -m unittest discover -s scripts/tests -p 'test_compare_kpoints.py'
git diff --check
```

The current full matrix also includes `si-supercell-stretch-wide` (12 cases
rather than 11); to reproduce the originally measured main sweep, set
`DFT_MIXING_CASES` to the first 11 names in the convergence table. Set it to
`si-supercell-stretch-wide` for the separate wider-band sweep. Set
`DFT_MIXING_CASES=synthetic-fcc-z4-p3,si-gamma` and `DFT_MIXING_BETA=0.1` or `0.8`
for the parameter runs, with distinct evidence directories. `DFT_CLI_BINARY`
and `DFT_SI_UPF` override paths. Never commit generated UPFs. The model suite's
three-repeat timing records are retained in `/private/tmp/p10d-repeat-{0,1,2}.log`.

Final gates: build/test/fmt/clippy all PASS; five solver tests passed; material
characterization passed its invariants and explicitly retained failures; eight
existing ignored engine tests passed. Reference-comparison Python tests passed.
Initial missing-test-target launch (before creating the tests directory) exited
101; the two pre-recovery characterization failures and profile denial are
recorded above. No numerical tolerance was relaxed after a failure.

## Acceptance and remaining limitations

Implementation, published-candidate survey, strategy selection, reproducible
parameters, charge tests, explicit failures, fixed-point certificate retention,
model improvement, material/shape/k/band diversity, default bit identity and
existing QE comparisons are verified. Physical successes and failures are
limited to the exact systems and settings above; no general metallic, slab,
chemical-accuracy or speedup claim is made. Further work includes heterogeneous
local screening, potential mixing, Johnson Broyden, adaptive trust-region logic,
and independent investigation of the failed stretched Si states. None is hidden
in the new defaults. Quiet-dev08 performance acceptance remains with Claude.

Artifacts/dashboard, README, baselines, manifests and Cargo.lock were explicitly
excluded by the ticket and are untouched. This report replaces the earlier
ownership-blocked prototype status after the orchestrator authorized the two
minimal connection files. Changes remain uncommitted for review, as required.
