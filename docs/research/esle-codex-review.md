# Independent ESLE architecture review

Task: P8R-ESLE-REVIEW. Reviewed 2026-11-03. This assessment was developed from the directive and repository implementation without consulting another ESLE design document. Only this report is changed; source files, gate records and session records are outside this ticket's ownership. No build, test or new scientific calculation was run. Measurements below are explicitly attributed to existing repository records, not independently reproduced here.

Evidence labels: **Verified in code** means static inspection of the cited implementation, not execution. **Recorded evidence** means an existing document reports a measurement or requirement. **Inferred/proposed** identifies my scientific interpretation or a future contract. All paths are repository-relative; line references identify the reviewed snapshot.

## Assessment

**Inferred/proposed:** ESLE can begin after baseline validation with additive initial-density support, an immutable state envelope and conservative comparisons. It needs no new crate initially. The current system already returns structured orbital results; replacing the SCF architecture would be unnecessary. The immediate scientific value for insulating nonmagnetic Si is reproducibility, initialization robustness and continuation validation, not an expectation of many metastable solutions. Neither an SCF convergence certificate nor a collection of starts establishes a global minimum.

**Verified in code:** weighted Bloch SCF exists (`crates/dft-engine/src/scf.rs:291`, `:444`, `:577`), despite the Gamma-only introduction and future-k-point wording in `docs/architecture.md:3`. **Recorded evidence:** ADR 0011 reports synthetic nonzero-k, reciprocal-translation and full/reduced quadrature checks (`docs/decisions/0011-kpoints.md:156–179`). The validation record ends with a Gamma-only Fourier-density comparison (`docs/validation.md:201–206`); `docs/p7-plan.md:7` still requires real Si shifted/unshifted k-grid QE comparisons. **Inferred:** implementation and recorded synthetic validation must not be represented as completed real-material k-point validation. This review does not authorize crossing that gate or changing interfaces during it (`docs/research/esle-directive.md:9`).

## 1. Implementation inventory

### SCF inputs, outputs and numerical seams

**Verified in code:**

- `ScfProblem::new(cell, atoms, potentials, options)` selects legacy Gamma; `with_kpoints` adds a point slice and `OccupationPolicy` (`crates/dft-engine/src/scf.rs:269–319`). Species resolve atom labels to borrowed parsed pseudopotentials. Options include cutoff, requested/extra bands, FFT override, linear/Pulay mixing, density/energy/eigen tolerances, eigensolver tolerance schedule, iteration/subspace limits, seed and divergence controls (`:103–141`). Functional identity is checked against the UPF label (`:363–370`). There is no independent charge input: electron count is the sum of ionic valences, required positive, integral and even (`:391–417`).
- `run` returns success-only `ScfResult` or `ScfError`; `run_outcome` also returns a final evaluated result on iteration-budget exhaustion with `converged=false` (`:625–640`, `:795–859`). Numerical errors still return errors. Some errors contain iteration history, others only a cause (`:24–46`); not every failure yields a recoverable state.
- `ScfResult` holds energy audit, first-k compatibility eigenvalues/orbitals, common occupations, all k-point eigenpairs, optional sampled gap, final output density, integrated charge, iteration history, Ewald breakdown and timings (`:197–215`). `EigenSolution` holds eigenvalues, coefficient matrix, per-band residuals, iteration count and residual history (`crates/dft-solvers/src/lib.rs:117–124`). Geometry, pseudopotential hashes and basis labels are not in `ScfResult`.
- Energy is the sum of kinetic, ionic-local, nonlocal, Hartree, XC and Ewald terms. Local zero mode, band sum, double counting, band-form and Harris–Foulkes values are diagnostics, not extra summands (`crates/dft-engine/src/energy.rs:8–31`, `:97–159`). Per-k audits use the same total density and are weighted (`:177–198`; `scf.rs:591–614`).
- `LinearOperator` exposes dimension and vector/block application; `Hamiltonian` applies kinetic, FFT local potential and nonlocal projectors (`crates/dft-solvers/src/lib.rs:30–44`; `crates/dft-physics/src/hamiltonian.rs:328–389`). Davidson accepts optional coefficient guesses and a preconditioner (`crates/dft-solvers/src/lib.rs:175–215`). These interfaces already separate physics from eigensolving.

**Verified in code:** the initial density is explicitly uniform, `N_e/Omega` (`scf.rs:498–500`, `:647`), and the initial orbital guesses are `None` (`:655`). Subsequent iterations reuse each k block (`:852–856`). `step_kpoints(n, initial, tolerance)` already accepts density plus optional matrices in quadrature order, but performs one map, not the certified mixed SCF driver (`:517–622`). `step` only supports a single k point (`:501–513`).

**Inferred/proposed:** add a full-run entry point at `run_impl` initialization, preserving the existing loop and certification. A density-only start initializes `n`; orbitals-only must first construct the weighted density with the fixed occupations; both should have an explicit meaning: density determines the first potential and orbitals are eigensolver guesses. If both are declared a reconstructed physical state, check their density agreement; if deliberately inconsistent, record them as an initialization recipe, not a certified state. Validate charge, shape, finite values, layout, positivity, basis metadata and orbital rank before entering the loop. Shape checks alone in `step_kpoints` are not restart compatibility validation.

### K points, density, orbitals, occupations, spin

**Verified in code:** `KPoint` is fractional reciprocal coordinates plus positive normalized weight; validation checks weights without normalization (`crates/dft-geometry/src/kpoints.rs:4–57`). Grid generation and time-reversal reduction are explicit; reduction keeps the first representative and merges equivalent k and -k, with no space-group symmetry (`:60–116`). Input order survives in the engine's `Vec<KState>`, each containing its own basis, FFT map and Hamiltonian (`scf.rs:246–251`, `:444–471`).

**Verified in code:** `PlaneWaveBasis::at_k` canonicalizes k, selects by kinetic momentum k+G, and sorts by G norm and Miller labels (`crates/dft-basis/src/lib.rs:78–150`). Coefficients are complex, dimensionless, stored by band columns with entry `i+j*rows` (`crates/dft-numerics/src/linalg.rs:26–60`). Density is `Array3<f64>` in electrons/Bohr^3, obtained from normalized orbital coefficients with occupation and volume factors (`crates/dft-physics/src/density.rs:1–31`). The engine sums these densities with k weights (`scf.rs:568–587`). In formulas,

```text
psi_nk(r) = Omega^(-1/2) sum_G c_nk(G) exp(i(k+G).r)
n(r) = sum_k w_k sum_n f_n |psi_nk(r)|^2
sum_G |c_nk(G)|^2 = 1;  sum_k w_k = 1;  integral n dr = N_e.
```

The Bloch convention is also recorded in `docs/decisions/0011-kpoints.md:9–24`.

**Verified in code:** the same `Vec<f64>` occupations is used at every k: two electrons in each of the lowest `N_e/2` bands, zero in empty bands (`scf.rs:474–476`). There is no spin channel, spinor, magnetization or chemical-potential solve; XC is explicitly unpolarized PZ81 (`crates/dft-physics/src/xc.rs:1–42`). Partial occupations are explicitly rejected (`scf.rs:177–185`, `:299–306`). The k-point API requires an empty band and checks the sampled indirect gap against its policy and eigen tolerance; legacy Gamma has no gap certificate (`:414–417`, `:861–885`). A sampled positive gap says nothing about unsampled k regions.

### Persistence and reproducibility

**Verified in code:** `density.bin` is the final output density as little-endian f64; `density.json` contains schema version, dimensions, cell, volume, units, ordering, expected/integrated electrons and density checksum (`crates/dft-cli/src/run.rs:279–283`). `results.json` contains input and UPF hashes, tool/compiler identity, atoms/cell, cutoff/grid summary, occupations/eigenvalues, energies, controls, history, convergence, timings and threading; new k inputs add point/weight/eigenvalue records (`:300–322`). Legacy Gamma omits that block. Writes occur after SCF and are separate operations (`:439–441`); `docs/cli.md:49–57` documents overwrite and partial-output risks.

**Verified in code, bounded absence finding:** inspection/search of `crates/*/src` finds no SCF checkpoint reader or restart command using those density files. The serialization writer does not persist orbital coefficients, per-k Miller lists, mixer history, current mixed input density, Davidson subspaces or random-generator progress. The solver's “restart” is an internal subspace restart, not disk recovery (`crates/dft-solvers/src/lib.rs:175–215`; `docs/decisions/0006-eigensolver-api.md`, Operator and solver API). Existing density data can support a future density warm start, but cannot reconstruct the exported occupied projector or resume an interrupted trajectory exactly. Reconstruction also needs accessible matching UPF bytes, not merely their digests. The three files lack a transaction manifest binding all payloads into one immutable attempt.

**Verified in code:** randomness is local SplitMix64 with configured seed (`crates/dft-solvers/src/lib.rs:196`, `:347–365`); SCF gives every k solve that seed (`scf.rs:544–559`). Density and audit accumulation follow vector order (`:577–614`; `energy.rs:179–198`). Basis sorting is explicit (`dft-basis/src/lib.rs:110–115`); quadrature validation uses pairwise summation but production density reduction is an ordered accumulation (`dft-geometry/src/kpoints.rs:34–54`). faer uses `Par::Seq` (`dft-numerics/src/linalg.rs:135–143`, `:186–200`); FFT axes are called serially (`dft-numerics/src/lib.rs:74–84`). The engine reports one execution thread (`dft-engine/src/lib.rs:12–13`).

**Inferred:** these choices support same-build repeatability, not portable bitwise identity. Changed point/species/atom order, basis cutoff boundary rounding, compiler/CPU/backend versions, floating-point instruction selection, reduction trees, degenerate eigenvector gauges and future parallelism can alter bits. Timestamps and durations already prevent byte-identical result JSON (`run.rs:301`, `:310`). Preserve execution order separately from canonical comparison order. A state content hash identifies an artifact, not a physical equivalence class; exclude volatile metadata from reproducibility comparisons. Record source-tree/build digest and dependency/toolchain versions as well as HEAD, since a commit alone does not identify uncommitted code (`docs/validation.md:115–119`).

## 2. ElectronicState identity and compatibility

**Inferred/proposed**, grounded in the gaps in `ScfResult` (`scf.rs:197–215`) and the CLI envelope (`run.rs:279–322`): separate a versioned persisted state from an attempt record and from an exact runtime checkpoint.

Mandatory for a current-model state:

- Schema/convention version, artifact ID, payload digests, explicit units and layout, and capability flags (density present; occupied orbitals present or absent).
- Cell columns in Bohr, ordered species/atom identities and fractional positions; explicit neutral electron count; pseudopotential content hashes with algorithm, valences and NC model. Paths/names are locators, not identity.
- Exact XC model/parameterization and implementation/convention version; explicit unpolarized, fixed-integer occupation policy; no implicit default for spin, charge or smearing. Record whether the legacy no-gap contract or the insulating certificate applies.
- Actual canonical k list and weights, execution order, full-grid/shift information where applicable, time-reversal reduction and representative mapping; explicit absence of space-group reduction. Store the actual quadrature even if generation parameters exist.
- Wavefunction cutoff in Ha, density-support convention, FFT dimensions/origin/layout, canonical Bloch convention, and per-k ordered Miller labels (or a versioned deterministic recipe plus basis digest, verified on load). A plane-wave count is insufficient.
- Real density and its integral; per-k occupations and occupied rank; computed band counts and eigenvalue arrays in Ha for a solver-produced state. Validate sum_k w_k sum_n f_n=N_e.
- State status, final residuals and certificate, requested thresholds, energy audit and energy-definition version for evaluated states; solver/build provenance, initialization recipe, seed, and optional parent represented explicitly as absent for independent starts.

An arbitrary starting density is an `InitialGuess`, not a fabricated evaluated state: it has no invented energy or eigenvalues. A failed attempt can lack a state entirely. Evaluated but unconverged states remain storable but ineligible for accepted-state ranking.

Optional: full occupied coefficient blocks (mandatory for projector comparison), empty-band coefficients, full SCF/eigensolver histories, raw input, timings, host metadata, parent links and path parameter when applicable. An occupied-factor representation is already a compact projector representation; a dense plane-wave projector is unnecessary. A lossy compressed representation must carry an error bound used in comparisons. The minimum density-only record supports density restart but cannot claim complete orbital-state reconstruction. For scientific ESLE archives I recommend retaining occupied blocks by default. Fermi level is optional and convention-labelled: a gapped fixed-occupation system has no uniquely selected value inside its gap, and the current solver does not compute one.

Future, gated fields: spin densities/spinors and spin convention, local moments plus integration/projector definition, Hubbard subspaces and occupation matrices with local gauges, U/J and double-counting definitions, fractional occupations, chemical potential, smearing scheme/width, entropy and the choice of internal/free/extrapolated energy. These fields must not silently activate unsupported physics.

### Compatibility is an operation, not one Boolean

**Inferred/proposed:** return `Exact`, `Convertible(plan)`, or `Incompatible(reasons)` for a named operation; never silently convert.

| Field | Same-problem restart / state comparison | Explicit transfer or conversion |
| --- | --- | --- |
| UPF hash, XC/model, spin/occupation policy | Mismatch is an error | Cross-model guesses are a separate future experiment, never the same branch identity |
| Electron count, atom species/positions, cell | Must match for same-geometry state identity and energy ranking | Deformation/displacement can seed a new problem with parent provenance; different charge is outside current engine |
| Cutoff and ordered per-k basis | Exact reuse requires matching basis or verified permutation | Same-cell zero embedding/truncation by Miller label is possible; record discarded norm, reorthonormalize and recertify |
| FFT dimensions/conventions | Exact density reuse requires match | Bandlimited Fourier resampling can be valid; nonlinear XC quadrature changes require energy reevaluation |
| K set and weights | Require equivalent weighted set with verified coordinate/band mapping | Reordering is exact semantically, not necessarily bitwise; reduced/full conversion requires valid time-reversal transformation; different sampling allows density seeding only initially |
| Number of empty bands | Occupied rank/policy must match | Different empty-band counts can still permit occupied-projector comparison; compare only common spectral windows and regenerate missing solver guesses |
| Solver controls/build | Must match for an exact replay claim | Warm restart may change controls; accepted comparison needs both states certified under a common accuracy protocol |
| Schema, payload integrity, units | Unknown convention, corrupt payload, nonfinite data or inconsistent dimensions is an error | Only explicit tested schema/unit migration is permitted |

Equal cell volume is not equal cell geometry. Periodic translations, atom permutations or origin shifts require an explicit consistent mapping of atoms, density and orbital phases; do not fit them away implicitly. Energy ranking uses identical physical and numerical settings, or reevaluates both candidates at common settings. A mathematical cross-cutoff density/projector distance can be reported as a discretization study, not automatic same-state merging. Evidence for these distinctions: basis/cell dependence (`dft-basis/src/lib.rs:78–154`), weighted nonlinear audit (`energy.rs:177–198`), and observed grid dependence (`docs/validation.md:139–149`).

## 3. State comparison and tolerance calibration

All formulas below are **inferred/proposed**, derived for compatible current-model states. They operate on certified results and return the metric vector, compatibility decision, calibration identity and an `EquivalentWithinResolution`, `Distinct`, or `Undetermined` classification. No direct coefficient equality is used. Fix the spatial origin and k mapping first. Do not merge records destructively; keep an equivalence association with its evidence and policy version.

### Density distances

On a common cell/grid, let M be the number of points and delta n_j=n_A(r_j)-n_B(r_j). Use

```text
d_n2 = sqrt((Omega/M) sum_j (delta n_j)^2) / N_e   [Bohr^(-3/2)]
d_ninf = max_j |delta n_j|                        [electron/Bohr^3]
d_n2_scaled = sqrt(Omega) * d_n2                  [dimensionless]
d_ninf_scaled = (Omega/N_e) * d_ninf              [dimensionless].
```

The first formula exactly matches `density_distance` (`crates/dft-engine/src/scf.rs:888–899`); the infinity metric is a proposal. Report both raw and scaled definitions, never an unlabeled “relative L2.” Parseval gives `d_n2 = sqrt(Omega sum_G |delta n_G|^2)/N_e` for the documented FFT normalization (`crates/dft-numerics/src/lib.rs:18–28`, `:74–84`). Different grids in the same cell require common Fourier support and a checked resampling operation, including Nyquist treatment; discard loss must be reported. Grid infinity distance is a sampled maximum, not a guaranteed continuum supremum.

False merge: L2 averages away localized changes and scales with cell/electron count; coarse grids miss features; distinct orbital projectors can share a scalar density. Infinity distance helps local sensitivity but still misses unsampled structure. False split: inconsistent origins, grid interpolation, truncation ringing or different convergence accuracy. Both distances require context-specific calibration and cannot establish orbital identity alone.

### Occupied-subspace distance

For each matched k, embed occupied coefficient matrices A_k and B_k into the union of the same-cell Miller basis, inserting zeros for absent rows. Verify column orthonormality first. Let ranks be r_A,k and r_B,k, S_k=A_k^H B_k and P_A,k=A_k A_k^H. Then

```text
||P_A,k - P_B,k||_F^2 = r_A,k + r_B,k - 2 ||S_k||_F^2
D_P^2 = sum_k w_k ||P_A,k - P_B,k||_F^2
D_P_normalized = sqrt(D_P^2 / sum_k w_k (r_A,k+r_B,k)).
```

All are dimensionless. Under the current occupation contract ranks agree, and electron density matrices are twice these projectors; keep the convention explicit. The formula follows by expanding the trace of the squared projector difference. It is invariant to phases, band permutations and occupied-unitary rotations. A rank mismatch is an occupation-compatibility failure, not a distance to quietly threshold. Singular values of S additionally expose principal-angle outliers that a weighted aggregate could hide.

Algebraic cost, not a measured performance claim: forming overlaps takes O(sum_k N_G,k r_A,k r_B,k) arithmetic and O(r_A,k r_B,k) temporary overlap storage per k, besides coefficient storage. It avoids forming N_G,k-by-N_G,k projectors. The repository's adjoint GEMM adapter is suitable (`crates/dft-numerics/src/linalg.rs:161–200`). Roundoff may make the trace expression slightly negative; clamp only within a documented measured arithmetic error envelope, otherwise reject. Do not use square-root cancellation to justify an invented projector tolerance.

Different cells have different Hilbert-space domains. For tracking only, choose and record a unitary pullback to fractional coordinates, `u_tilde(s)=sqrt(Omega) u(A s)`. Equal Miller labels then identify Fourier functions on the unit torus; compare union-basis coefficients under that explicit affine identification. This is a chosen continuation similarity, not an invariant physical distance between unrelated cells. Changed lattice orientation, cell representation or origin needs a mapping; unrelated/supercell k folds should initially be rejected. Projection onto only the intersection can erase distinguishing components and must not be used as an unqualified equivalence test.

False merge: averaging over k hides a local change, missing/filtered basis components erase differences, or an invalid time-reversal mapping identifies different spaces. Report maximum per-k discrepancy as well. False split: occupation-boundary changes near gap closure, nonorthogonal input, truncation, or inappropriate cell mapping. Degeneracy wholly inside the occupied manifold is harmless; a degeneracy crossing its boundary invalidates this fixed-rank insulating comparison.

### Spectrum and energy

Compare sorted eigenvalues at each matched k over the same specified band window:

```text
d_eps_inf = max_k,n |eps_A,nk - eps_B,nk|                    [Ha]
d_eps_rms = sqrt(sum_k w_k sum_n (eps_A,nk-eps_B,nk)^2 / b) [Ha]
d_E = |E_A - E_B|                                           [Ha/cell].
```

Here b is the common number of compared bands. Also report occupied and empty windows separately, sampled direct/indirect gaps, and per-atom energy only with a declared atom count. The existing energy and spectrum storage supports this (`scf.rs:197–218`, `:861–871`; `energy.rs:26–31`). Require the same potential-zero convention; any aligned-spectrum diagnostic must record the shift and cannot replace raw energies. Cross-geometry spectra are tracking diagnostics only.

False merge: isospectral or near-degenerate-energy states can have different densities/projectors; energy is stationary near a solution and can hide appreciable state error. False split: changed potential zero, empty-band windows, quadrature, cutoff or solver accuracy. An energy window is a ranking/consistency filter, never sufficient identity. Preserve a lower-energy state even when the classifier thinks it equivalent; investigate disagreement under tighter common settings.

### Measured starting evidence and missing calibration

**Recorded evidence, not new measurements:**

| Repository observation | Value and scope | Implication |
| --- | --- | --- |
| Real Si tolerance sweep | Difference to tightest energy is 1.4210854715202004e-14 Ha at density tolerance 1e-7 Bohr^(-3/2); the loosest listed point differs by 1.3873346915715956e-12 Ha (`docs/validation.md:143–147`) | Energy sensitivity for that fixture, not a universal noise bound |
| Real Si grid sweep | Grid energies differ by at most 7.8e-8 Ha; gap changes by 1e-6 eV (`docs/validation.md:139–149`) | Grid error exceeds that fixture's SCF energy sensitivity and is not a branch difference |
| Sampled Fourier-density QE comparison | Largest sampled difference 4.0e-9 electron/Bohr^3; cubic-equivalent coefficient spreads 2.1e-11 (ours) and 6.8e-9 (QE) electron/Bohr^3 (`docs/validation.md:205`) | Fourier sample noise, not full-density L2 or real-grid infinity noise |
| Synthetic k full/reduced comparison | Energy difference 1.3322676295501878e-15 Ha; maximum density difference 1.5681900222830336e-15 electron/Bohr^3 (`docs/decisions/0011-kpoints.md:165–166`) | Quadrature-equivalence evidence on that synthetic fixture only |

The real density evidence is sampled/aggregate only (`docs/validation.md:206`). No measured projector-distance noise distribution or ESLE classifier false-merge rate is established by these records. The QE acceptance bounds in `docs/validation.md:50–57` are validation gates, not same-engine clustering tolerances; the cross-code energy offset discussed at `:197–198` is not random numerical noise.

**Inferred/proposed calibration method:** for each fixed physical/numerical problem, measure all metrics on repeated same-seed runs, distinct orbital seeds, tighter SCF/eigen tolerances, fresh post-SCF rediagonalization, serialized reloads, and equivalent full/reduced k representations. Keep grid/cutoff studies separate as discretization error. Use gauge rotations/permutations as known-equivalent controls and deliberately perturbed, independently reconverged candidates as challenge cases. Establish a measured same-state envelope for every metric, select and document thresholds above that envelope only when separated from resolved candidate differences, then test on held-out repetitions. If envelopes overlap, label unresolved and retain both. Record threshold units, normalization, case coverage and version; do not borrow Si thresholds for a larger cell or smaller gap. No numerical merge tolerances are proposed before those measurements. Pairwise thresholding is not transitive: avoid chain-merging a cluster whose extremes are distinguishable.

## 4. Legitimate multi-start under current physics

**Verified in code:** uniform density and seeded orbital initialization are implemented (`scf.rs:498–500`, `:647–655`; `dft-solvers/src/lib.rs:196–215`). Changing the solver seed changes orbital guesses, not the initial density.

**Verified in code:** `Pseudopotential::raw()` exposes the underlying reader model (`crates/dft-upf/src/parser.rs:49–64`). `PP_RHOATOM` is parsed, size-checked where declared, and checked for finite samples (`:385–446`). It is therefore incorrect to say it is unavailable. There is no normalized atomic-density accessor or production superposition initializer in the inspected source: SCF still uses `uniform_density`. Raw access is explicitly documented as unsuitable for direct numerical-kernel use (`:61`).

**Inferred/proposed:** atomic superposition is legitimate after a validated adapter establishes the UPF radial-density convention, radial weights, integrated valence, endpoint behavior and periodic Fourier transform. Do not mistake the radial distribution for pointwise three-dimensional density or for the local potential. Test its integral against valence before use; finite samples alone do not establish normalization or positivity. Construct periodic atom-centered densities with structure factors and check the total on the actual grid. Do not silently substitute guessed atomic densities when required data are absent.

Other legitimate proposed starts are a compatible previous density, orbitals with the same fixed occupations, and smooth real density perturbations with zero integral. For a perturbation, impose conjugate Fourier symmetry, omit its zero mode, restrict support to the target density representation, and bound its amplitude to keep density nonnegative. Save seed, generator/version, modes and amplitude. Use convex damping rather than pointwise clipping, which changes Fourier support and charge; the existing Pulay safeguard follows this principle (`scf.rs:830–848`). An arbitrary positive normalized density is a potential-building guess, not proof of finite-basis orbital representability. The next orbital solve constructs an admissible output density.

**Inferred:** multiple fixed points are possible for nonlinear KS equations, but the present nonmagnetic, integer-filled Si-like problem likely has one accessible solution for the chosen structure/model. Symmetry-breaking charge patterns or small-gap nonlinear response could yield alternatives; numerical failures, underconvergence and basis errors can imitate them. Do not prescribe magnetic moments, local orbital occupations or promoted-band “states” through this solver's fixed Aufbau occupations. Multi-start validates basin robustness, seed independence, restart integrity and false-split resistance even if every accepted run coincides. It cannot validate false-merge resistance on genuinely distinct states without an appropriate additional fixture.

A particular risk is an invariant warm-start subspace that omits lower eigenstates: small eigen residuals certify eigenpairs, not that they are the lowest. Davidson accepts user guesses after orthonormalization and stops on residuals (`crates/dft-solvers/src/lib.rs:197–215`, `docs/decisions/0006-eigensolver-api.md`, Alternatives and limits). Independent fresh guesses and small-system dense checks are needed; a gap between the returned bands does not cure omitted lower states.

## 5. Continuation and density transfer

**Inferred/proposed**, using the existing FFT normalization, basis labels and density construction (`dft-numerics/src/lib.rs:18–28`; `dft-basis/src/lib.rs:149–154`; `dft-physics/src/density.rs:24–26`):

For unchanged cell and displaced atoms, the old grid density can seed a rebuilt problem. Rebuild ionic local/nonlocal and Ewald terms from the new atoms; never retain the old cached `ScfProblem` Hamiltonian (`scf.rs:373–473`). This is a new state with a parent link, not a restart of identical conditions.

For A_old to A_new, use a declared affine fractional-coordinate transfer:

```text
n_guess(A_new s) = (Omega_old/Omega_new) n_old(A_old s).
```

This preserves electron count by change of variables. Transform to reciprocal coefficients, match integer Miller labels under the chosen cell map, multiply density coefficients by the volume ratio, truncate modes unrepresentable on the target grid and zero-pad missing ones. Set/check G=0 as N_e/Omega_new. Respect negative Nyquist conventions, conjugate pairs and collision checks; do not interpolate flattened arrays or match Cartesian G vectors from different cells as if they were identical. Check and report integral, imaginary leakage, discarded spectral norm, positivity and round-trip loss. A final positive multiplicative charge correction may address numerical integration loss, with the factor recorded; reject substantial defects instead of disguising a bad map as renormalization. Any positivity repair should preserve charge/support and be recorded as a changed guess.

Once atomic densities are independently validated, a displacement-aware extrapolation can use

```text
n_guess = T[n_old - n_atomic_old] + n_atomic_new,
```

where T is the same charge-preserving transfer and both atomic sums use consistent pseudopotentials and normalization. The residual should integrate to zero; correct/check it and damp the update to avoid negative density. This avoids merely dragging bonding density with atoms without accounting for atomic motion, but is a heuristic, not an adiabatic theorem or a branch certificate. No such extrapolator exists in the inspected SCF path (`scf.rs:647–655`).

Risks: cutoff membership changes, grid ringing, XC quadrature changes, arbitrary cell/origin mapping, closing sampled gap, and loss of localization can all mimic branch changes. Even density transfer to a different k mesh is only a guess; its energy and occupations require a fresh solve. Periodic independent starts, reverse-path checks and smaller path increments can diagnose hysteresis and branch-assignment ambiguity. Failure to converge is a recorded failure, not proof that a physical branch disappeared. Tracking should allow ambiguous links/crossings and retain failed attempts rather than force one label through every geometry.

## 6. Interfaces, crate ownership and least invasive sequence

**Inferred/proposed sketches only:** these are contract sketches, not compiled code. Named envelope/metadata/error types below denote new owned types implementing the rules above. Keep numerical values in the existing `Array3`, `ComplexMatrix` and Hartree-named f64 style; do not introduce backend or serde types into solver traits. Evidence for that style is `scf.rs:103–215`, `dft-solvers/src/lib.rs:30–53` and `dft-physics/src/data.rs:1`.

```rust
// dft-engine: validated metadata + owned numerical initial guesses.
pub enum InitialState {
    Density { density: Array3<f64>, basis: StateBasis },
    Orbitals { blocks: Vec<OrbitalBlock>, basis: StateBasis },
    DensityAndOrbitalGuess {
        density: Array3<f64>,
        blocks: Vec<OrbitalBlock>,
        basis: StateBasis,
    },
}
pub struct ScfRequest {
    pub initial: Option<InitialState>, // None preserves existing uniform start.
}
impl ScfProblem {
    pub fn run_from(self, request: ScfRequest) -> Result<ScfOutcome, ScfError>;
}

// Pure envelope; metadata describes the problem actually evaluated.
pub struct EvaluatedState {
    pub metadata: StateMetadata,
    pub outcome: ScfOutcome,
    pub initialization: InitializationRecord,
}
pub struct AttemptRecord {
    pub id: AttemptId,
    pub parent: Option<StateId>,
    pub initialization: InitializationRecord,
    pub result: Result<StateId, FailureRecord>,
}
// ElectronicState is the versioned, validated archive representation.
// It may omit orbital payloads, but exposes that loss of capability explicitly.
```

`ScfOutcome` already supplies structured convergence and result data; do not introduce a replacement for its numerical fields. Failures must retain the recipe and controls even when the engine cannot return a final state. The archive should distinguish initial, evaluated-unconverged and certified states in its type/status validation.

```rust
// Pure comparison in dft-engine; no filesystem access.
pub trait StateComparison {
    fn compare(
        &self,
        a: &ElectronicState,
        b: &ElectronicState,
        context: &ComparisonContext,
    ) -> Result<ComparisonReport, ComparisonError>;
}
pub trait ExplorationPolicy {
    fn next(
        &mut self,
        context: &ExplorationContext,
        attempts: &[AttemptRecord],
    ) -> Option<InitializationPlan>;
}

// dft-cli boundary: persistence implementation + driving loop live here.
pub trait StateStore {
    fn load(&self, id: &StateId) -> Result<ElectronicState, StoreError>;
    fn save(&mut self, state: &ElectronicState) -> Result<StateId, StoreError>;
    fn record_attempt(&mut self, attempt: &AttemptRecord) -> Result<(), StoreError>;
}
```

The CLI resolves load/save operations between solver invocations. Policies emit plans referencing IDs; they do not read files in a Hamiltonian call. `save` validates metadata/payload consistency and commits an immutable manifest after payload writes; partial writes must never appear as complete states. Unknown schema/model capabilities, malformed dimensions and digest mismatches fail closed. Classification and branch registries preserve attempts, failed outcomes, ties and ambiguous matches. Start with a deterministic fixed-budget policy; no learned policy is needed to exercise these seams.

**API impact:** adding `run_from` and new types is additive. Existing `new`, `with_kpoints`, `run`, `run_outcome`, `step`, `step_kpoints`, `LinearOperator` and `Davidson::solve` need no public signature changes. `run_impl` would change privately to accept validated guesses. Adding a required field to public `ScfOptions`/`ScfResult` breaks downstream struct literals; replacing first-k fields or changing `run` return type is breaking. Avoid both via the request and envelope. A read-only basis/occupation descriptor accessor is an additive extension worth considering so clients do not reverse-engineer private `KState` order. Archive schema changes are separate from preserving CLI schema-1 output. Any eventual removal of compatibility views needs its own ADR and migration.

**Crate count:** zero new crates initially. Use pure state/compatibility/comparison/exploration modules in `dft-engine`, persistence and command driving in `dft-cli`, density transfer/atomic-density physics in existing `dft-physics`/`dft-upf` modules when separately ticketed. Keep the existing eigensolver and numerical adapters. This respects `ownership.toml:18–25`: engine can depend on physics/solvers/basis/geometry/upf, while CLI owns I/O. Engine backend types can use the existing physics data facade; do not add a forbidden direct engine-to-numerics dependency or put numerical comparison in the I/O-only CLI. Atomic-density support would touch existing crates, not require more crates. Extract a landscape crate only if independent clients later justify an ownership/dependency change; a generic `dft-core` replacement would duplicate the existing decomposition.

Sequence after the baseline gate: (a) specify model/basis compatibility and additive initialization ADR; (b) density-only `run_from` and envelope, preserving default regressions; (c) immutable density-state persistence and strict reload, recording all failures; (d) optional orbital persistence/restart and projector comparison with calibration; (e) deterministic multi-start and conservative grouping; (f) explicit density transfer and branch tracking against independent starts. Atomic superposition can join after its radial-density validation. Each step is independently reviewable and does not require a Hamiltonian redesign.

## 7. Scientific risks and independent-review checklist

The following is **inferred/proposed** for every future ESLE ticket, grounded in current convergence gates (`scf.rs:739–817`, `:901–908`), fixed occupations (`:474–476`), energy audits (`energy.rs:97–159`) and the directive's review requirement (`docs/research/esle-directive.md:33`).

- **Electron conservation:** check density integral and zero Fourier mode before/after load, perturbation, resampling and mixing; test changed volume, grid and atom displacement. Do not mask wrong units with normalization.
- **Occupation consistency:** verify weighted occupation sum, occupied rank, orthonormality and band ordering at every k; test deficient/invariant guesses and independent lower-spectrum checks. Reject unsupported spin/fractional occupations; retain the distinction between legacy Gamma and sampled-gap certificates.
- **Basis compatibility:** verify cell, k coordinates, weights, Miller ordering, Bloch phase, coefficient normalization, Nyquist and FFT layout. Test reordered equivalent bases, truncation loss, corrupted shapes and incompatible cells. Equal dimensions do not imply equal bases.
- **Energy validity:** compare only common model/structure/charge/sampling/discretization; preserve raw components and G=0 convention. The functional-difference audit is not a bound on unknown global-minimum error. Reevaluate transferred guesses before ranking.
- **Convergence interpretation:** preserve density, energy, occupied-residual, charge and fixed-point gates together; the returned state is the candidate output, while the extra map supplies a certificate (`scf.rs:795–815`). A lower unconverged energy is not an accepted branch. Record failures without treating them as branch annihilation.
- **False merging/splitting:** gauge and band-rotation controls, localized/k-local differences, near-gap cases, sampling changes, nontransitive clusters and missing orbitals must be exercised. Keep unresolved states separate; persist metric/calibration versions and reasons.
- **Reproducibility:** repeat same build/seed, vary seeds independently, load-save-load exactly, reject corrupt/partial archives, and check numerical equivalence across supported builds separately from bitwise replay. Record all attempted starts, including failures and budget exhaustion. Changing warm guesses must not silently change the physical model.
- **Continuation bias:** compare inherited and fresh starts at the same geometry; reverse paths and independent branches must be distinguishable in provenance. A parent link is computational ancestry, not proof of physical branch continuity.

**Unsupported targets and dependency order:** Al requires fractional occupations, a chemical-potential solve, smearing/entropy and consistent energy definitions, plus metallic k-convergence validation before ESLE uses it as a metallic benchmark (`scf.rs:299–306`). Spin requires spin-resolved density/occupations/XC and Hamiltonian/energy validation before magnetic initializations; Fe is not a current-model target (`xc.rs:1–42`). DFT+U requires defined correlated local projectors/overlap metric, occupation matrices, interaction parameters, double counting and variational potentials/energy before orbital-pattern exploration; ordinary NC pseudopotential projectors are not automatically Hubbard subspaces. The inspected Hamiltonian contains kinetic/local/NC nonlocal terms only (`hamiltonian.rs:328–356`). For the directive's magnetic correlated targets, validate spin before spin-dependent DFT+U; metallic occupations are a separate prerequisite for metallic systems, not intrinsically for every insulating +U case. PBE is planned, not implemented by the current XC path (`docs/p7-plan.md:8–11`; `xc.rs:39`). New materials also need supported, validated pseudopotentials and independent reference evidence; MgO support is not established merely because multiple species are representable.

## 8. Disagreements and unresolved semantics in the directive

**Inferred/proposed:** the directive is scientifically sound in separating lowest-found search from tracking and refusing global-optimum claims (`docs/research/esle-directive.md:6`, `:29`). I would refine these points before implementation:

1. “ElectronicState” currently mixes a physical solution, a trial guess, a reproducible artifact and a runtime checkpoint (`:16`). Give these distinct contracts. Density-only restart is feasible without exact trajectory resumption; neither should imply recovery of orbitals that were never saved.
2. Fermi level, local moments and occupation matrices must not be mandatory populated fields for the present insulating scalar model (`:16`). Use explicit unsupported/absent capabilities. A single scalar similarity is premature; retain several metrics and an unresolved classification.
3. “Occupation control needs DFT+U or similar” (`:15`) is a useful scope guard, not a general physical necessity: a separately defined constrained/excited-state DFT method could also control occupations. Such a method would change the problem and is not enabled by forcing the present solver to keep arbitrary occupied bands.
4. No universal physical distance exists between orbitals in different cells without an identification map. Branch IDs are hypotheses supported by continuation/comparison, not intrinsic hashes (`:18`). Do not equate solver failure with disappearance, or smooth ancestry with adiabatic continuity through a crossing.
5. “A lower-energy solution found elsewhere” (`:29`) only invalidates a lowest-known claim when evaluated for the same physical problem under compatible accuracy/settings. A lower energy at a different geometry, charge, XC or pseudopotential does not do so.
6. Al is a later metallic-physics validation target, not an immediate fixed-occupation ESLE smoke test (`:29`). The larger magnetic/correlated materials need independent method validation first. A Si multi-start result with one cluster is useful but cannot establish discovery capability for those targets.
7. Keep the proposed logical layers without creating the named collection of crates (`:22`). The current structure already provides the relevant boundaries. Learned exploration and claims based on distinct-solution counts are premature until classifier resolution, failures and baseline policies are measured (`:26`, `:29`).

Open decisions for architecture review: archive orbital retention policy; the exact-replay versus warm-restart promise; the initial supported cell-identification maps; calibration/uncertainty governance; and whether symmetry-related broken-symmetry states should be distinct records grouped into a separate symmetry orbit. Until those decisions are made, default to immutable records, explicit transfer and conservative non-merging.
