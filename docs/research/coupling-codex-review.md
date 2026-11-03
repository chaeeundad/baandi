# Independent scientific coupling review

Task: P9R-COUPLING-REVIEW. Scope: design review only, against this checkout. No source, schema, dependency, session-registry or validation-artifact changes. This assessment starts from the directive and implementation, independently of any future coupling design. ESLE documents are used only for state identity and continuation. No new physical calculation or QE run is claimed.

Evidence labels below distinguish **verified in code**, **recorded evidence** (documentation or upstream documentation), and **inferred** (physical deductions and proposed contracts). Repository paths and line numbers refer to the reviewed checkout. Interface sketches are proposals, not compilable additions or promises of implemented capabilities.

## Assessment

**Inferred:** endorse the directive's separation of numerical core, physical contracts and external drivers. The smallest useful future increment is a fixed, real, cell-periodic electron-energy field with an independently audited energy contribution. A universal energy/force calculator is premature. The current engine can supply baseline energies, densities, eigenpairs and diagnostics, but cannot yet execute even case A through its SCF API. Neither serialization nor a trait can supply missing variational physics.

**Verified in code:** the engine is weighted Bloch, spin-unpolarized, fixed-occupation LDA-PZ (`crates/dft-engine/src/lib.rs:1`, `scf.rs:180–185`, `scf.rs:367–370`, `scf.rs:391–417`). All k blocks share one effective potential (`scf.rs:529–537`); results expose no forces or stress (`scf.rs:198–215`). External-field inputs are absent from `ScfProblem` and its constructors (`scf.rs:232–245`, `scf.rs:269–330`). This is a favorable seam for an additive future feature, not an existing coupling capability.

**Inferred:** commit now to semantics: typed field meaning, explicit unsupported errors, separate inner/outer convergence, complete energy ownership, immutable environment identity and unavailable observables. Defer concrete ABI, runtime implementation, restart guarantees, charged-cell models and derivatives. Preserve every k-point gate and legacy entry point. User decisions are listed in section 9; this review does not request implementation approval.

## 1. Local potential, energy and zero modes

**Verified in code:** the exact data path is:

| Quantity | Storage and assembly | Evidence |
| --- | --- | --- |
| Ionic electron potential energy | `ScfProblem.ionic_g: Array3<Complex64>`; `ionic_real: Array3<f64>` after inverse FFT | `crates/dft-engine/src/scf.rs:236–237`, `:431–443` |
| Hartree electron potential energy | local `h: HartreeGrid`, with complex `potential` and `reciprocal_potential` | `scf.rs:530`; `crates/dft-physics/src/hartree.rs:10–14`, `:63–84` |
| XC electron potential energy | local `x: XcGrid`, real `potential`, functional energy and potential integral | `scf.rs:531`; `crates/dft-physics/src/xc.rs:80–99` |
| Effective local field | real `effective = ionic_real.clone()`, adding `h.re + x` at each sample, copied into each k Hamiltonian | `scf.rs:532–537` |
| Hamiltonian multiplication | private real `potential`, multiplied with the FFT-transformed periodic orbital, then transformed back | `crates/dft-physics/src/hamiltonian.rs:263–268`, `:340–356` |

**Inferred:** add validated `U_ext(r)` in Ha to `effective` before `set_local_potential`, identically at every k. Retain its separate identity and array; do not disguise it as a pseudopotential or change `ionic_g`. Require matching cell/grid/origin and real, finite values. Public `set_local_potential` already permits operator-level tests (`hamiltonian.rs:310–319`), but the SCF rebuilds the field each step: an external caller cannot inject a persistent field by that method alone.

**Verified in code:** `EnergyAudit` contains `kinetic`, `local`, `local_zero_mode`, `nonlocal`, `hartree`, `xc`, `ewald`, `band_sum`, `double_counting`, `band_form`, `harris_foulkes`, and `functional_difference_bound` (`crates/dft-engine/src/energy.rs:8–24`). The total sums only kinetic + local + nonlocal + Hartree + XC + Ewald (`:26–31`); `local_zero_mode` is already inside `local`, not an extra energy. `local_energy` computes the full-grid Parseval integral of output density and ionic coefficients, including zero mode (`:35–55`, `:102`, `:150`). Double counting is integral n_out times (V_H[n_in] + V_xc[n_in]); the band and Harris forms are distinct away from self-consistency (`:116–159`). All k audits use the same total density and common nonlinear terms occur once under normalized weights (`:177–198`, `scf.rs:576–614`).

**Inferred:** for a fixed field, the direct functional needs a separate `external_electron = (Omega/M) sum_j n_out,j U_j` in Ha/cell. In the band form that contribution is already in eigenvalues: do not add it a second time or subtract it as Hartree/XC double counting. The current band/Harris expressions acquire the fixed-field contribution through the band sum; the direct expression must be extended to agree. Audit bounds must continue using the same input field and output orbitals, including in the fresh certificate. For general embedding functional F_emb, the band reconstruction instead needs `- integral n_out U_emb[n_in] + F_emb[n_out]`, plus any separately owned environment/ion terms. A scalar field alone supplies no F_emb.

**Verified in code / recorded evidence:** zero-mode policies are different physical conventions, not one universal zeroing rule:

* Ionic local: `ionic_local` selects `alpha_z()` at Miller zero (`hamiltonian.rs:70–77`). Its finite Coulomb-subtracted constant is `4*pi/Omega * short_integral(0) + pi*Z/Omega`, with an erf-subtracted integrand (`crates/dft-upf/src/radial.rs:213–233`). `local_zero_mode = N * ionic_g[0].re` (`energy.rs:103`, `:150`). This finite remainder must remain.
* Hartree: `V_H(G)=4*pi*n(G)/G²` for G nonzero, `V_H(0)=0`, energy with the one-half factor (`hartree.rs:63–76`). It solves for the mean-subtracted electron number density, not a neutral total signed charge distribution (`docs/physics-conventions.md:232–252`).
* XC: pointwise derivative and its mean are retained; no zero-mode subtraction (`xc.rs:80–99`). XC is not an electrostatic potential.
* Ewald: reciprocal zero excluded (`crates/dft-physics/src/ewald.rs:159–174`), point self term and uniform background `-pi Q_ion²/(2 Omega eta²)` included (`:180–193`). This is the ionic component's compensated convention, not an implemented user-selectable charged-cell correction.

**Inferred:** a finite supplied periodic U has a legitimate mean U0. Preserve it and record its reference. A zero-mean convention may be requested explicitly, with the removed constant recorded; silently applying the Hartree rule to U would change energies. For fixed N electrons and electron-only constant U0, H' = H + U0 I, every eigenvalue shifts by U0 Ha, and total energy shifts by N U0 Ha/cell. Density, orbitals up to gauge/degenerate rotations, occupations and gaps are unchanged provided the same insulating state is followed. There is no one-half factor. Positive U0 raises electron energies.

For an electrostatic reference shift phi0, U0 = -phi0 in atomic units, so the electronic shift is -N phi0; ionic electrostatic energy shifts by +sum Z_I phi0. The combined QM shift is (sum Z_I - N) phi0, zero for a neutral QM cell when both contributions are included. Any additional external sources/background must be accounted for consistently. The electron-only N U0 shift is therefore not a gauge-invariant total charged-system observable. A periodic net charge has no unmodified periodic Poisson solution; choosing a background and reference defines a different physical model, not merely a numerical stabilization. Current SCF fixes N=sum Z_I (`scf.rs:386–403`), so charged-cell results cannot be advertised today.

## 2. Units, signs and geometry contract

**Verified in code:** density is positive electron *number* per Bohr cubed (`crates/dft-physics/src/density.rs:1–26`), not signed charge. Hartree's positive kernel is repulsive electron potential energy (`hartree.rs:71–76`); ionic local has the attractive Coulomb tail (`radial.rs:220–227`). **Recorded evidence:** atomic units set the magnitude of elementary charge to one, not the electron charge to +1 (`docs/physics-conventions.md:5–7`). **Inferred:** declare q_e=-1, rho_el=-n, and U_el=q_e phi=-phi in atomic units. For volts, U_el in eV is numerically minus phi in volts, then convert eV to Ha once. An electric field is E=-grad phi, not -grad U_el.

**Inferred:** validated types must distinguish dimensions *and physical role*. Two Ha arrays can represent entirely different observables. Use private constructors and immutable grid metadata; raw f64 values remain acceptable inside established numerical loops.

```rust
struct ElectronEnergyHa(f64);
struct ElectrostaticPotentialHaPerE(f64);
struct SignedChargeE(f64);
struct ElectronNumberDensityPerBohr3(f64);
struct SignedChargeDensityEPerBohr3(f64);
struct ScalarGrid<Q> { context: GridContext, values: Vec<Q> }
struct GridContext {
    cell_columns_bohr: [[f64; 3]; 3], origin_bohr: [f64; 3],
    dims: [usize; 3], layout: Layout, periodic_axes: [bool; 3],
    fft_convention: ConventionId, reference: PotentialReference,
}
fn electron_energy(phi: ElectrostaticPotentialHaPerE) -> ElectronEnergyHa;
```

**Verified in code / recorded evidence:** `dft-units` currently imports uom Energy, Length, Pressure and Ratio, not charge or potential (`crates/dft-units/src/lib.rs:55–60`). Its energy/length conversions are available (`:154–175`); ADR 0012 keeps internal f64 and places conversion at boundaries (`docs/decisions/0012-units-library.md:7–27`). **Inferred:** extend dft-units later with uom electric-potential/charge dimensional conversions and explicit semantic newtypes. Do not claim existing `Hartree` prevents phi/U confusion. Adapters convert volts/eV/Angstrom and charges per voxel versus density before validation; physics accepts atomic units only. Generated constants and existing conversions avoid new hardcoded physical constants. UPF Ry conversions stay solely in dft-upf (ADR 0008).

**Verified in code:** lattice columns, reciprocal transpose and coordinate transformations are in `crates/dft-geometry/src/lib.rs:44–52`, `:86–101`; FFT normalizations in `crates/dft-numerics/src/lib.rs:18–31`, `:74–84`; canonical negative Nyquist and C-order maps in `crates/dft-basis/src/lib.rs:241–267`. **Recorded evidence:** origin is implicitly zero and samples are r=A(i/n1,j/n2,k/n3) (`docs/physics-conventions.md:63–76`). **Inferred:** initial field support should require this exact grid, with all three axes periodic. Later adapters may resample with explicit Fourier support, origin phase and discarded-weight reports. For samples at r0+A s, the sampled Fourier coefficient gains exp(+iG.r0); conversion back to the global-origin coefficient removes that phase. A real field with arbitrary high-frequency samples can alias into matrix elements; the density-sphere rule does not license unrestricted imported fields. Project/validate representable support, reject ambiguous Nyquist cases, and converge quadrature. Do not silently interpolate or strip imaginary parts. The existing Hartree complex-Nyquist caveat is documented at `docs/physics-conventions.md:234–243`.

## 3. External charges and boundary conditions

**Verified in code:** Hartree calls `checked_density` on each input sample (`hartree.rs:59–62`), so it is not a general signed-charge Poisson API. It removes the zero coefficient without testing total neutrality (`:68–76`). `scf.rs:391–403` establishes neutral, positive even valence count separately. `ewald::Ion` accepts signed point charges, but its current SCF caller supplies ionic valences (`ewald.rs:9–13`, `:89–93`; `scf.rs:386–391`).

**Inferred:** distinguish these future inputs:

| Input | Conversion and owner | Required physical policy |
| --- | --- | --- |
| Prescribed U_el grid | No Poisson solve in core; adapter supplies electron energy | Declared reference, origin, periodicity and fixed-field energy |
| Signed rho_ext grid | Future pure `dft-physics` Poisson operator, or external solver with declared identical convention | phi_G=4*pi*rho_ext,G/G², U_G=-phi_G; neutrality/background and zero mode explicit |
| Point charges q_a,R_a | Separate regularized/Ewald representation; not delta samples on one voxel | Singularity, self subtraction, periodic images, source-source and ion-source terms; reject coincident singular sources |

For normalized Fourier coefficients, rho_ext,G=(1/Omega) sum_a q_a exp(-iG.R_a). Point charges are not bandlimited. Grid deposition requires a specified shape/width, charge conservation and convergence checks; a bare point attractive potential also risks unphysical electron spill into the classical region. Environment pseudization, repulsion and QM/MM boundary treatment belong to a documented model, not an automatic conversion.

Prefer separate electron Hartree and external-source solves initially. Simply putting rho_ext into the current positive-density Hartree routine gets sign and meaning wrong. A total signed-charge solver would also require replacing/reconciling existing electron-ion/Ewald decomposition to avoid counting Coulomb terms twice. The electron-source cross energy is integral n U_ext, without one-half; source self/environment energy and ion-source energy require explicit ownership. If rho_ext changes in response to n, section 4 applies.

A nonzero cell integral requires an explicit compensating background or rejection; a fixed-N neutral QM subsystem plus charged environment is not automatically neutral. Periodic images remain even for net-neutral sources. Three-dimensional periodic Coulomb electrostatics does not implement isolated molecules, 2D slab Coulomb truncation, dipole corrections or open-boundary FEM electrostatics. Vacuum supercells are periodic approximations requiring convergence, not isolated boundary conditions. A uniform field has nonperiodic scalar potential -E.r and cannot be encoded as a smooth periodic scalar field; Berry-phase electric enthalpy and slab sawtooth/dipole methods require distinct capabilities. Neither appears in the current Hamiltonian/SCF path (`hamiltonian.rs:329–356`, `scf.rs:529–537`). Do not infer support from a `periodic_axes` metadata field.

## 4. Variational consistency and observables

**Inferred from the energy path** (`energy.rs:99–159`): at fixed geometry a converged eigenproblem is not by itself proof of a coupled variational energy. For stationary E[n,R], F_I=-partial E/partial R_I at stationary n; orbital-response terms cancel only with correct functional derivatives and sufficient convergence. Include nonlocal pseudopotential derivatives and Ewald derivatives in the ordinary force as well. Fixed-cell plane waves avoid atom-centered basis Pulay terms, but stress/cell changes need basis/cutoff, metric, quadrature and field-frame treatment; this is not an argument that stress is trivial.

| Class (inferred contract) | Energy and derivative conditions | Missing terms / restriction |
| --- | --- | --- |
| Fixed U independent of n,R | E_DFT + integral n U is variational; delta E/delta n adds U | Electron-only fixed lab field has no explicit dU/dR. Physical electrostatic embedding also needs E_ion-ext=sum Z_I phi_ext(R_I), giving +Z_I E_ext(R_I). |
| Density-dependent U[n] | Require named F_emb with U=delta F_emb/delta n, including any auxiliary variables minimized to stationarity | Integral n U is generally not F_emb. Linear reciprocal response can yield a one-half factor; nonlinear reaction fields cannot borrow that factor. Dielectric/cavity derivatives and polarization self energy matter. |
| Geometry-dependent U[n,R] | Same functional condition, plus explicit geometry derivative | Force includes -integral n partial U/partial R for fixed-in-n U, ion-field terms, environment parameter/position derivatives and moving-grid/interpolation derivatives. External point charges receive reaction forces in the driver. |
| Time-dependent U(t) | Instantaneous ground-state surfaces may be evaluated quasistatically if individually stationary | Explicit work term integral n partial U/partial t; no conserved static total energy claim. No real-time propagation, electronic friction or nonadiabatic forces exists here. |
| Outer-driver updates | Freeze environment for each inner solve; separately converge joint variational equations | Inner certificate does not certify outer stationarity or solvent equilibrium. Coupled free energy may include entropy and reservoir work; never rank with E_DFT alone. |

The classes overlap: a field can depend on density, geometry and time simultaneously. Declare dependencies separately from consistency/derivative availability. An empirical potential with no energy primitive can at most return electronic diagnostics and an explicitly partial/frozen-field energy; it must not claim a conservative coupled total, forces, stress or thermodynamic ranking.

**Verified in code:** inspection of engine/physics source APIs finds no force/stress evaluator; `ScfResult` has none (`scf.rs:198–215`), and CLI reports energy/density/eigenpairs instead (`crates/dft-cli/src/run.rs:279–313`). **Recorded evidence:** pressure conversion explicitly does not enable stress (`docs/cli.md:264`). **Inferred:** unavailable must be a typed result with reason, never a zero vector/tensor. MLIP energy-only data generation is possible for current neutral unpolarized insulating LDA models; force/stress training targets and MD forces are unavailable. A finite-difference script would be a separately validated derived observable, not an implemented analytic force. Preserve state/environment labels, cutoff/k sampling, convergence, UPF hashes and the energy's accounting scope in every training record. Do not merge different environments or electronic states at identical atom coordinates.

Electrostatic observables need another distinction: V_eff includes XC and pseudopotential short-range physics and is not phi_total. Electronic phi is -V_H in this convention; ionic point-charge phi is not simply minus the entire local pseudopotential, and nonlocal projectors have no scalar electrostatic equivalent. Define valence-only versus reconstructed/all-electron source and core conventions before exporting a QM electrostatic field to a solvent solver. Spin density and a unique metallic Fermi level cannot be fabricated; the current result supplies a sampled gap, not a Fermi-level field (`scf.rs:198–215`).

## 5. Minimum separate interfaces and ownership

**Inferred proposals**, grounded in `ownership.toml:7–25`, current `ScfOutcome` (`scf.rs:247–257`) and ESLE's proposal status (`docs/architecture/electronic-state.md:1–21`). These sketches use abstract value types to express contracts, not a dependency shopping list.

```rust
trait CapabilityProvider {
    fn capabilities(&self) -> CapabilitySet;
    fn validate(&self, request: &CalculationRequest) -> Result<ValidatedRequest, Unsupported>;
}
struct CalculationRequest {
    request_id: RequestId, canonical_payload_digest: Digest,
    model: ModelContext, environment: EnvironmentContext,
    terms: Vec<TermSpec>, initial: Option<ContinuationRequest>,
    requested: ObservableSet, convergence: InnerCriteria, budget: ResourceBudget,
}
struct CalculationResult {
    request_id: RequestId, model_digest: Digest, environment_digest: Digest,
    status: SolveStatus, energy: Availability<EnergyLedger>,
    observables: ObservableBundle, quality: QualityRecord,
    state_id: Option<ElectronicStateId>, lineage: InitializationLineage,
}
enum Availability<T> { Available(T), Unsupported(Reason), NotComputed(Reason) }
enum SolveStatus { Converged, BudgetExhausted, Cancelled, Failed(Failure) }
```

Capabilities must be negotiated for combinations, not just Boolean flags: current neutral/even-N, nspin=1, LDA-PZ, NC subset, fixed occupations and sampled insulating gap, 3D periodic Bloch model. Today advertise no external term, no charged-cell control, no spin, PBE, force, stress, restart or time dynamics. Requested unsupported observables can be rejected before solving or explicitly returned unavailable according to a declared strictness policy; never silently ignored.

```rust
struct Dependencies { density: bool, geometry: bool, time: bool }
enum Consistency {
    FixedLinearEnergy,
    Variational { functional: FunctionalId, quantity: EnergyOrFreeEnergy },
    PotentialOnly,
}
trait ExternalFieldTerm {
    fn identity(&self) -> TermIdentity;
    fn dependencies(&self) -> Dependencies;
    fn consistency(&self) -> Consistency;
    fn evaluate(&self, context: &ImmutableContext, n: &ElectronDensity)
        -> Result<TermEvaluation, TermError>;
}
struct TermEvaluation {
    electron_potential: ScalarGrid<ElectronEnergyHa>,
    functional_energy: Availability<EnergyContribution>,
    explicit_derivatives: Availability<GeometryDerivatives>,
    accounting: ContributionOwnership,
}
trait Observables {
    fn get(&self, observable: ObservableKind) -> Availability<QualifiedObservable>;
}
```

The first field implementation supports only fixed linear energy, with the numerical integral audited independently. The trait is a mathematical provider contract, not permission for remote I/O inside SCF. Each term declares which electron, ion, source-self and environment energies it owns, including duplicate contribution IDs. Report electronic partial and coupled total/free energy separately. Reject a claimed variational implementation whose derivative/energy validation is absent. Nonlocal embedding needs a separate Hermitian-operator interface and trace accounting; do not stretch this scalar trait to cover it. Density-dependent implementations may require new inner/outer algorithms; the directive's “primarily external” goal cannot promise zero core changes for all future functionals.

```rust
trait StateContinuation {
    fn assess(&self, state: &StateDescriptor, target: &ModelEnvironment,
              purpose: ContinuationPurpose) -> Compatibility;
    fn prepare(&self, state: &StateDescriptor, plan: ApprovedTransfer)
        -> Result<InitialGuess, StateError>;
}
enum Compatibility { Exact, TransferRequired(TransferPlan), Incompatible(Vec<Reason>) }
struct ContinuationRequest {
    parent: ElectronicStateId, target_branch: Option<BranchId>,
    purpose: ContinuationPurpose,
}
trait CalculationRuntime {
    fn submit(&mut self, request: CalculationRequest) -> Result<JobId, SubmissionError>;
    fn poll(&self, job: JobId) -> JobStatus;
    fn cancel(&mut self, job: JobId) -> CancelAcknowledgement;
    fn result(&self, job: JobId) -> Availability<CalculationResult>;
}
```

**Recorded evidence:** ESLE distinguishes initial guesses, evaluated states and attempts (`docs/architecture/electronic-state.md:20–21`), proposes model fingerprints and payload hashes (`:28–49`), and restricts comparisons/transfers by operation (`:51–58`). Its independent review distinguishes exact replay and warm restart (`docs/research/esle-codex-review.md:87`) and treats branch continuity as a hypothesis (`:280–292`). **Inferred:** add environment content digest, gauge, source/functional parameters, boundary/background convention and accounting version to model compatibility. A field change produces a new model/state; it may permit a warm transfer, not equal-model energy ranking. BranchId is optional external tracking metadata, not an intrinsic wavefunction hash or proof of ground state. Compare branches at the *same* environment; cross-environment curves need explicit thermodynamic accounting. Preserve independent fresh starts alongside continuation. Follow ESLE's prose/ADR 0016 placing persistence in CLI, not its ambiguous `save/load` functions adjacent to an engine sketch (`electronic-state.md:71–75`).

**Inferred lifecycle:** negotiate and validate, resolve immutable inputs, submit, solve, seal result, then let the external driver update the environment. Same request ID plus same canonical payload returns the existing job/result; same ID with different content is a conflict. Include input arrays, units, physical model, initial-state/seed choice and tolerances in identity; exclude timestamps/host display metadata. Each environment update is a new request with parent/outer-iteration lineage. A request ID is not an electronic-state ID. Specify canonical float/ordering rules before promising content-addressed identity.

Cancellation is cooperative at future safe boundaries and acknowledged separately from completion. Retain only completed, internally consistent steps; partial output carries unconverged status, actual residuals and missing-observable reasons. Failed attempts may have no state. Atomic publication belongs to adapters/runtime. Resource limits distinguish inner iterations, outer iterations and wall/resource budgets. A repeated ID deduplicates work; it does not guarantee bitwise replay across builds/backends. Current `run_outcome` already distinguishes budget exhaustion (`scf.rs:637–643`), but cancellation/checkpoint recovery are future work.

**Inferred ownership mapping:** physical units in dft-units; geometry in dft-geometry; grid support in dft-basis; mathematical field/Poisson/derivative operations in dft-physics using dft-numerics; request validation, energy accounting, state compatibility and outcome types in dft-engine; files/checksums/network/process lifecycle in adapters, initially dft-cli. Engine already obtains unit/array/backend types through `dft_physics::data`, preserving existing permitted edges (`docs/architecture.md:29–33`; `ownership.toml:22–25`). No `dft-core` crate exists in this topology. No new crate or direct engine dependency is justified now. Frozen ownership/schema changes require later tickets.

Commit to these separations now as design constraints. Implementing even fixed U needs SCF assembly, both energy paths and certificate coverage later. Defer general term ABI, derivatives, nonlinear embedding, asynchronous runtime, array import and state transfer until independently scoped gates. Keep backend-specific arrays and serialization out of the physical contract; this is backend neutrality, not implemented GPU capability.

## 6. Serialization and interoperability

**Verified in code:** CLI writes headerless little-endian f64 density in C order, with grid, lattice, volume, electron-number units, normalization and SHA-256 in JSON (`crates/dft-cli/src/run.rs:279–283`, `:442–443`; **recorded evidence:** `docs/cli.md:168–175`). This is density export, not a generic potential or restart reader. **Inferred:** retain it now. A future envelope must add array role, explicit origin, periodic axes, reference/gauge, model/environment digests, payload length, schema/version and convergence status without changing legacy outputs. A checksum detects changed bytes; it does not validate physical meaning or provide authenticity. Readers must reject wrong lengths/shapes, overflow, nonfinite data, unsupported versions and path escapes before allocation/use. Write payloads first and publish an immutable manifest last; do not reinterpret an incomplete directory as a checkpoint.

**Recorded evidence (upstream sources, consulted for this review); inferred adoption costs:** file-format interoperability does not require copying an implementation or inheriting its code license. Binding/library/transitive licenses still need review at the pinned version.

| Option | Evidence and role | Assessment for this repository |
| --- | --- | --- |
| HDF5 | Multidimensional self-describing storage; upstream BSD-style/3-clause license: [HDF Group](https://www.hdfgroup.org/solutions/hdf5/), [license](https://github.com/HDFGroup/hdf5/blob/develop/LICENSE) | Good later multi-array archive candidate. Typical C-library/Rust binding path adds native build/link and Apple Silicon validation. Keep optional in an adapter, with explicit physics metadata. |
| NetCDF | Named scientific array data and portable formats: [Unidata FAQ](https://docs.unidata.ucar.edu/netcdf-c/current/faq.html); [BSD-3 licensing](https://www.unidata.ucar.edu/software/netcdf/licensing) | Useful if a consumer already uses its dimension conventions. NetCDF-4/HDF5 adds layered native dependencies; classic variants have different capabilities. No reason to adopt both HDF5 and NetCDF now. |
| Zarr | Chunked compressed arrays; Python implementation MIT: [upstream](https://github.com/zarr-developers/zarr-python) | Candidate when chunk/object-store exchange is required. Pin format version/codecs, storage transaction semantics and checksums; review the actual Rust implementation independently. Python is not required in core. More moving parts than current single payload. |
| Cube | ASE documents cube reading/writing, including origin: [format documentation](https://docs.ase-lib.org/ase/io/formatoptions.html) | Useful visualization/interchange export. Text precision, unit conventions and limited metadata make it a poor authoritative checkpoint. Format use does not require Gaussian code; this review does not approve any proprietary implementation dependency. |
| openPMD | Metadata schema over storage such as HDF5/ADIOS2/JSON, not a competing single binary encoding: [reference API](https://github.com/openPMD/openPMD-api) | Worth a later mesh/unit metadata mapping. Reference API carries LGPL/GPL license texts and a C++/backend build stack; review exact component obligations rather than assuming permissive core linkage. Defer adoption until a consumer needs it. |
| ASE/extxyz | [ASE](https://docs.ase-lib.org/about.html), [extxyz specification/parser](https://github.com/libAtoms/extxyz) (MIT parser; native C/libcleri/PCRE2 build path) | Prefer later external structure/energy dataset adapter. ASE is a Python ecosystem with LGPL-2.1-or-later licensing to review for the selected release; keep it optional. extxyz is not a dense-field archive; link fields and complete provenance by digest. Explicitly map lattice ordering, eV/Angstrom, force and stress sign/order conventions. |

**Inferred recommendation:** preserve LE-f64 + JSON as the authoritative near-term representation; first optional structure adapter can use extxyz, and first richer field/archive adapter can use HDF5 after dependency approval. Reconsider Zarr for an actual chunked-storage requirement and openPMD for an actual metadata consumer. No format alone resolves field sign, total-potential meaning, gauge, branch identity or double counting. No new packages installed, no upstream code copied, and no performance claims or binding compatibility certification are made here.

## 7. Use cases A–G

Each row is **inferred** from the verified core limits in sections 1–4 and the directive's use cases (`docs/research/coupling-directive.md:91–107`). Present feasibility always means neutral, even-electron, spin-unpolarized, insulating NC LDA-PZ Bloch calculations; PBE, metals, smearing, forces and stress are not available. Legacy Gamma does not provide the new sampled-gap certificate (ADR 0011:90–94).

| Case | Data exchanged | Required core capability | External responsibilities and physical constraints | Present feasibility / gap |
| --- | --- | --- | --- | --- |
| A Fixed electrostatic embedding | U_el grid + cell/reference/environment ID in; density, energy ledger, eigenpairs/residuals out | Fixed scalar injection and audited integral nU; optional electrostatic outputs | Build compatible periodic field; distinguish electron-only perturbation from ion-inclusive electrostatics; preserve gauge and neutrality | Baseline density export and operator primitives exist; SCF injection/ledger absent. First bounded extension, not usable today. |
| B Mean-field explicit solvent / DFT-CES family | Averaged solvent charge/potential or model-specific terms; QM density and qualified electrostatic observables | A plus validated electrostatic source/output definitions; richer terms only through supported operators/functionals | Solvent sampling/MD, averaging, thermodynamic functional, outer mixing/convergence and environment free energy | Baseline QM snapshots possible; no coupled loop or source reconstruction. Mean field of averaged density is not automatically an average free energy. No assertion that all CES-family terms reduce to scalar U. |
| C Polarizable QM/MM | Permanent charges, induced response/parameters, U; density and possibly reaction forces back | Functional/derivative-consistent embedding; forces only after independent implementation | Polarization solver, self-energy, exclusions, link/boundary treatment and force partition; joint stationarity | None of the coupled workflow is implemented. Fixed-field snapshots could follow A; no conservative dynamics or force training. |
| D Implicit solvation | Reaction field, dielectric/cavity state, functional components; density/cavity response | Functional derivative plus matching energy, later cavity/geometry derivatives | Continuum PDE solver, solvent free energy, non-electrostatic terms, boundary matching | Baseline density may feed a one-way external analysis; no self-consistent solvent. Constant-dielectric linear-response one-half rule cannot be assumed for density-dependent cavities. |
| E MLIP / active learning | Structures, energy ledger, state/environment IDs and quality; requested forces/stress | Qualified observables and eventually validated derivatives | Sampling, labeling, training/inference, uncertainty selection and provenance | Energy-only baseline datasets possible now. Force/stress columns unavailable, not zero; branches/environments must remain distinct. Electronic-state IDs still a proposal. |
| F FEM | Microscopic cell/energy/property data with orientation, units and validity domain | Energy scans available; later validated stress/response/elastic derivatives | Strain sampling, homogenization, temperature/microstructure, constitutive fitting and continuum BCs | External energy-versus-strain study can be designed now; no built-in stress/elastic tensor. A fitted 0 K tensor requires basis/k convergence and is not a finite-temperature constitutive law. |
| G ESLE in an environment | Initial state/branch hint and environment; evaluated states, lineage and qualified energies | A plus additive continuation and compatibility/state archive | Multi-start policy, branch hypotheses, common-environment comparison and total free-energy ranking | Seeded baseline reruns possible, but run starts uniform (`scf.rs:647–655`). ESLE and field APIs unimplemented; multiple electronic or coupled solutions are not established by current Si evidence. |

## 8. Validation specifications without feature implementation

No new tests or source are written in this task. **Inferred:** the following specifications can be prepared now; some can run against existing public operator/FFT routines, others must remain deferred integration gates. A mock passing is not a coupled SCF validation.

**Recorded evidence:** relevant measured scales in `docs/validation.md` are: direct/Harris difference about 9e-14 Ha for R1 (`:14`); tolerance-sweep energy differences to the tightest point of 1.4210854715202004e-14 Ha at density tolerance 1e-7 and up to 1.3873346915715956e-12 Ha at the loosest listed setting (`:143–147`); grid spread about 7.8e-8 Ha (`:139–149`); R3 occupied eigen residual 4.3152503285188984e-9 Ha (`:177`); R5 sampled density difference 4.0e-9 electron/Bohr^3 and our cubic-equivalent spread 2.1e-11 electron/Bohr^3 (`:205`). These are fixture-specific observations, not universal noise floors. The systematic alpha-Z offset (`:197–198`) and cross-code acceptance bounds (`:50–57`) are not same-engine perturbation noise. No measured full-density response, force or stress noise floor is given.

**Inferred tolerance rule:** calibrate same-grid, same-build, same-model repeats/tighter solves for each future fixture. Let epsilon_E be its observed energy uncertainty envelope and epsilon_n its measured density-integral uncertainty; record their construction and coverage. The R2 numbers motivate the starting calibration scale only. For central response D_h=[E(+h)-E(-h)]/(2h), with dimensionless amplitude h multiplying fixed W in Ha, energy uncertainty contributes at most (epsilon_E,+ + epsilon_E,-)/(2|h|); density uncertainty contributes |integral delta n W|, bounded by ||delta n||_2 ||W||_2. Add the measured stencil truncation envelope from an h sweep. Choose h on a resolved plateau; do not invent a universal decimal tolerance or loosen existing gates. Existing FFT tolerances below are recorded contracts, not measurements claimed by this review.

| Proposed test | Can be done before SCF field implementation? | Oracle and acceptance rationale |
| --- | --- | --- |
| Constant shift | Yes, existing Hamiltonian setter + dense eigensolver and `local_energy`; full SCF later | H(U0)c-H(0)c=U0 c, all eigenvalues shift U0, unchanged occupied projector/density, energy delta N U0. Both signs, multiple k, and explicit phi-to-U negative sign. Operator/FFT roundoff plus measured eigensolver residual sets bound; future SCF bound includes calibrated epsilon_E and U0 times charge uncertainty. |
| FFT/grid/origin | Yes for existing transforms; future importer refusal deferred | Single Fourier modes, skew and left-handed cells, even negative Nyquist, constant and translated fields; verify phase, Parseval, integral and cell rows/columns. Existing bounds are 1e-12 scaled FFT error and 1e-14 constant-mode error (`docs/physics-conventions.md:102–117`). No measured interpolation envelope yet: bandlimited exact fixtures first, then convergence for general interpolation. |
| Charge sign/Poisson | Yes as analytic fixtures and existing positive-density Hartree tests; signed API deferred | Positive source produces positive phi and negative U; negative source reverses it. Constant n has zero Hartree field under mean subtraction. Test Fourier multiplier and q/Omega normalization; signed external input must not pass as electron number density. Test net-charge rejection/background metadata separately. |
| Weak-field response | Frozen-operator identity now; reoptimized DFT response deferred | At stationary fixed-N branch dE/dh=integral n_h W. Compare central derivative with density integral on identical grid/k set. Use the uncertainty formula above and h/tolerance sweep. Avoid gap closure and PZ threshold/branch crossings; density and energy stationarity must be certified. |
| Double counting | Algebraic/direct-local fixtures now; complete external audit later | Fixed field counted once in direct and band forms; local_zero_mode not added twice. Artificial linear response functional with known half factor contrasts with fixed field. Add known density-independent ion/environment constant and verify it affects total once but not electron Hamiltonian. Existing bound is `functional_difference_bound` (`energy.rs:139–159`), to be rederived for new functional terms. |
| Serialization | Independent reading/re-encoding existing output can be tested now; restart later | LE-f64 byte-identical round trip, hash/length/order checks, recomputed N. Exact bytes need no floating tolerance; physical integral uses independently calibrated summation/charge allowance, not a QE energy bound. Reject corrupted/truncated payloads, mismatched origin/units and unsupported versions. |
| Metadata rejection | Existing constructor restrictions now; coupling validator later | Wrong cell/grid, k weights/order, spin, charge, PP digest, occupation/XC model, potential role/reference and environment identity. Exact digest/schema decisions; no numerical tolerance for incompatible physics. Distinguish approved transfer from exact restart. |
| Lifecycle/ESLE | Specification/mock state machine now, real cancellation/restart later | Repeated ID deduplication, changed payload conflict, cancellation before/after a completed step, partial results never marked converged, independent inner/outer statuses; parent does not certify branch identity. |
| Forces/stress | Deferred | Finite differences only after a complete stationary functional and derivative evaluator exist. Calibrate displaced/strained energy noise and basis/grid errors; no current force/stress tolerances can honestly be inferred from one Si energy sweep. |

**Inferred:** constant-shift tests should include no-field regression and direct/band/Harris agreement to catch the tempting Hamiltonian-only implementation. Keep long SCF checks ignored in debug with release commands as required by `AGENTS.md:26`; no tests are authorized for addition by this ticket. Preserve zero-field Gamma compatibility, full/reduced-k equivalence and original convergence thresholds. R6 records multi-k energy/eigenvalue comparisons but no k-run density comparison (`docs/validation.md:208–222`); do not upgrade this evidence to full state or response validation.

## 9. Risks, disagreements and decisions

**Inferred, high priority:**

1. “General purpose” must not imply every model fits a scalar field, has an energy primitive, or can be implemented with no core changes. Nonlocal, density-gradient and time-dependent terms need separate physics. The directive appropriately mentions this, but its flexibility goal needs this explicit limit (`docs/research/coupling-directive.md:43–48`, `:144–146`).
2. A field API without direct-energy integration, ion/environment partition and gauge metadata can produce plausible eigenpairs and incorrect total energies. The current independent direct/band audits are an asset to preserve (`energy.rs:99–159`).
3. Current “neutrality” is fixed electron count plus component zero-mode conventions, not a general charge/environment/background policy (`scf.rs:391–403`, `hartree.rs:68–76`, `ewald.rs:180–184`). Reject net-charged requests initially; do not introduce jellium silently.
4. V_eff is not a reusable electrostatic observable. Solvent/FEM consumers need explicitly defined sources, pseudopotential core treatment and compatible boundary conditions (sections 3–4).
5. Declaring a functional ID is not evidence of variational consistency. Density-dependent response, moving cavities and environment relaxation need derivative tests, and thermodynamic ranking needs the complete energy/free energy.
6. Coupling metadata belongs in ESLE compatibility, but branch identity is not request identity. Independent multistart and unresolved classifications must survive continuation; a failed solve is not a disappeared state (`docs/research/esle-codex-review.md:280–292`).
7. Documentation has historical stale statements: `docs/architecture.md:52–54` says only Gamma, while its introduction and code support Bloch k; its future-work paragraph also still lists k support. Use current code plus ADR 0011 and R6, without editing those files in this ticket. ESLE serialization interfaces remain proposals, not current features.

**Recorded evidence / inferred protection:** the directive requires unchanged active k validation (`docs/research/coupling-directive.md:17–20`). R6 already records three Si k comparisons with limited scope (`docs/validation.md:208–222`); it does not authorize altering that work. Do not change k canonicalization, weights, shifted basis ordering, grid selection, real scalar time-reversal assumptions, nonlocal projectors, occupation/gap gates, alpha-Z, XC or legacy outputs. Future zero-field requests must follow the existing path; no global refactor justified. A real cell-periodic scalar field preserves time reversal in the present spinless model, but a future magnetic/vector/nonlocal extension may not. Re-negotiate that capability instead of reusing reduced k sets blindly.

**Inferred decisions requiring user/architect approval before implementation:**

* Approve first scope as neutral fixed-N, exact-grid, real periodic electron-energy fields; decide whether phase one is explicitly electron-only or also includes a supplied electrostatic ion contribution.
* Choose the potential reference policy (preserved mean by default is recommended), allowed source regularization and charged-cell/background policy (initial rejection recommended).
* Set energy ownership and output names: electronic partial, interaction, source-self/environment, coupled total or free energy; decide whether potential-only models are rejected or diagnostic-only.
* Confirm state fingerprint extension and warm-transfer versus replay promises jointly with ESLE; decide branch-label governance without assuming metastable Si solutions.
* Approve observable order: fixed-field energy/density first; full ordinary and coupling forces/stress require separate physics tickets before MLIP force data or dynamics.
* Approve optional adapter dependencies only for a concrete consumer, with pinned licenses/maintenance/Apple Silicon checks. No new crate now.
* Approve calibration fixtures and measured response/derivative acceptance envelopes after the current k milestone. Existing tolerances remain fixed.

Suggested gates: document semantics now; complete/preserve baseline k validation; implement and validate fixed-field accounting; then observable/derivative contracts; then external responsive drivers; finally coupled-state exploration. A gate advances only on evidence for its physics, not because a transport adapter works.

## Review verification

`scripts/check.sh` exited 0: PASS build, PASS test, PASS fmt, PASS clippy. `git diff --check` exited 0. The new file was also checked explicitly with `git diff --no-index --check /dev/null docs/research/coupling-codex-review.md` (exit 1, expected for a new file compared with /dev/null; no whitespace diagnostics), since ordinary git diff does not inspect untracked files. No newly proposed coupling test is claimed as run. The only intended tracked-path addition is this document. Session/artifact JSON updates and all other directive deliverables are outside this ticket's owned files and remain for the orchestrator. Leave this file uncommitted for review.
