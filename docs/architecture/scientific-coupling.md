# Scientific coupling architecture (design; nothing here is implemented)

Status: design reconciled with the independent Codex review ([coupling-codex-review.md](../research/coupling-codex-review.md)). Directive: [coupling-directive.md](../research/coupling-directive.md). Physics: [coupling-physics.md](../research/coupling-physics.md). Use cases: [coupling-use-cases.md](coupling-use-cases.md). Roadmap: [../roadmap/coupling.md](../roadmap/coupling.md). ADR: [0017](../decisions/0017-coupling-commitments.md). ESLE state identity: [electronic-state.md](electronic-state.md).
The engine stays a focused Kohn-Sham component. Coupling methods (DFT-CES family, QM/MM, solvation, MLIP workflows, FEM exports) are external drivers and adapters; none of them, and no framework (LAMMPS, ASE, Python, REST, gRPC), becomes a dependency of the numerical core.

## 1. Gap analysis against the repository (verified in code)
| Need | Present state | Where |
|---|---|---|
| Insert an external field | none; the local potential of each SCF step is `ionic_real + V_H + V_xc`, rebuilt every step and copied to every k Hamiltonian; the only seam is that assembly | `crates/dft-engine/src/scf.rs` (~532), `crates/dft-physics/src/hamiltonian.rs` |
| Energy ledger | `EnergyAudit` (six physical terms plus diagnostics); public fields, so adding a field breaks struct literals | `crates/dft-engine/src/energy.rs` |
| Forces, stress | not implemented anywhere | searched `crates/*/src` |
| Spin, smearing, charged cells | none; electron count must be neutral and even | `scf.rs` |
| Signed-charge Poisson | none; Hartree accepts the electron number density only and drops G=0 | `crates/dft-physics/src/hartree.rs` |
| State identity, restart | none (ESLE proposal only) | [electronic-state.md](electronic-state.md) |
| Array export | little-endian f64 `density.bin` + JSON with SHA-256; no potential or checkpoint reader | `crates/dft-cli/src/run.rs`, `docs/cli.md` |
| Units | `dft-units` over uom, closed unit lists, strict parser | ADR 0012 |
Compatibility with the running k-point work: every item below is additive; k canonicalization, weights, shifted basis ordering, grid selection, nonlocal projectors, occupation and gap gates, alpha-Z, XC and the legacy outputs stay as they are; a zero-field request must follow the existing path exactly (regression: the real-Si Gamma density stays bit-identical).

## 2. Layers and ownership (no new crates now)
| Layer | Content | Home |
|---|---|---|
| Core | structure, grids, Hamiltonian, eigensolver, SCF, energy audit, (future) forces and stress, state types | existing crates; units in `dft-units`, geometry in `dft-geometry`, grids in `dft-basis`, field and Poisson math in `dft-physics`, request validation, energy accounting, state and outcome types in `dft-engine` |
| Coupling contracts | typed external inputs and outputs, immutable context, capability negotiation, validity metadata | pure modules in `dft-engine` (no I/O) and typed quantities in `dft-units` |
| Runtime (future) | outer iteration, request identity, checkpoint, cancellation | adapters, first in `dft-cli`; not in the numerical crates |
| Adapters (future) | files, Python/ASE, QM/MM and CES drivers, MLIP datasets, FEM exporters | outside the core, optional dependencies only for a concrete consumer |
| ESLE | states, branches, comparison, exploration | separate (ADR 0016) |
`ownership.toml` already allows these edges; a crate is split off only when an independent client needs it.

## 3. Interfaces (sketches, not code)
Several small interfaces instead of one calculator trait; none of them dictates the SCF or Hamiltonian internals.
```rust
// capability negotiation: combinations, not booleans
trait CapabilityProvider {
    fn capabilities(&self) -> CapabilitySet;
    fn validate(&self, r: &CalculationRequest) -> Result<ValidatedRequest, Unsupported>;
}
struct CalculationRequest { request_id: RequestId, payload_digest: Digest, model: ModelContext,
    environment: EnvironmentContext, terms: Vec<TermSpec>, initial: Option<ContinuationRequest>,
    requested: ObservableSet, convergence: InnerCriteria, budget: ResourceBudget }
struct CalculationResult { request_id: RequestId, model_digest: Digest, environment_digest: Digest,
    status: SolveStatus, energy: Availability<EnergyLedger>, observables: ObservableBundle,
    quality: QualityRecord, state_id: Option<ElectronicStateId>, lineage: InitializationLineage }
enum Availability<T> { Available(T), Unsupported(Reason), NotComputed(Reason) } // never a zero vector or tensor
enum SolveStatus { Converged, BudgetExhausted, Cancelled, Failed(Failure) }

// external terms declare how they depend and what they promise
struct Dependencies { density: bool, geometry: bool, time: bool }
enum Consistency { FixedLinearEnergy, Variational { functional: FunctionalId, quantity: EnergyOrFreeEnergy }, PotentialOnly }
trait ExternalFieldTerm {
    fn identity(&self) -> TermIdentity; fn dependencies(&self) -> Dependencies; fn consistency(&self) -> Consistency;
    fn evaluate(&self, ctx: &ImmutableContext, n: &ElectronDensity) -> Result<TermEvaluation, TermError>;
}
trait Observables { fn get(&self, k: ObservableKind) -> Availability<QualifiedObservable>; }
trait StateContinuation { /* Exact | TransferRequired(plan) | Incompatible(reasons), see ESLE */ }
```
Typed quantities distinguish physical role, not only dimension: electron potential energy (Ha), electrostatic potential (Ha per e), signed charge, electron number density, signed charge density; grids carry cell columns, origin, dimensions, layout, periodic axes, FFT convention, potential reference and boundary conditions. Raw f64 stays inside numerical loops.

## 4. Lifecycle (driver neutral)
1. Negotiate capabilities and validate the physical context (reject unsupported combinations, never ignore them).
2. Receive structure, coupling terms and an optional compatible state checkpoint.
3. Solve the electronic problem to the declared inner criteria.
4. Return observables, the energy ledger, state identity and quality, with unavailable items typed as unavailable.
5. The external driver updates the environment and decides whether to repeat.
Inner SCF convergence and outer coupled convergence are separate statuses; the inner certificate never certifies outer stationarity. One-way (fixed environment) and two-way (environment responds) workflows share the same request/result types. Request identity is idempotent: same id plus same canonical payload returns the existing result, same id with different content is a conflict; the request id is not an electronic-state id; every environment update is a new request with parent lineage; cancellation is cooperative and acknowledged separately; partial results carry unconverged status and the reason for each missing observable. Identity excludes timestamps and host metadata; canonical float and ordering rules must be fixed before content-addressed identity is promised. A repeated id deduplicates work but does not promise bitwise replay across builds.

## 5. Energy ledger and compatibility
The ledger names each contribution and its owner (electron-external, ion-external, source self and environment energy), so nothing is counted twice and an electronic partial energy is never reported as a coupled total or free energy. Coupled states are never ranked by DFT energy alone when the environment free-energy terms differ. Environment content digest, potential reference, background convention and accounting version are added to the ESLE model fingerprint: a field change is a new model (a warm transfer is allowed, equal-model energy ranking is not). Electronic metastability and multiple coupled QM/environment solutions are different phenomena and are labeled separately.

## 6. Serialization
Keep little-endian f64 plus JSON manifest with SHA-256 as the authoritative near-term representation; a future envelope adds array role, origin, periodic axes, reference, model and environment digests, payload length and schema version, writes payloads first and publishes an immutable manifest last. A checksum detects changed bytes, not physical validity. Candidate later adapters (from the review, licenses to be rechecked at the pinned version): extxyz for structure and label datasets, HDF5 for a multi-array archive, openPMD only if a consumer needs its mesh metadata; Zarr and NetCDF only for a concrete need; cube files for visualization only. No code is copied from any of them.

## 7. What may be committed now versus deferred
Commit now (semantics only): typed field meaning and sign, explicit unsupported errors, separate inner and outer convergence, complete energy ownership, immutable environment identity, typed unavailable observables. Defer: any concrete ABI, runtime, restart guarantee, charged-cell model, derivative implementation, nonlocal or density-dependent terms, array import, transport adapters. See ADR 0017.

## 8. Test plan prepared before any implementation
1. Constant potential shift: all eigenvalues shift by `U0`, density and occupied projector unchanged, energy `+N U0` for an electron-only shift and `(sum Z - N) phi0` for an electrostatic reference shift (zero for a neutral cell); both signs, several k, explicit phi to U sign.
2. FFT, grid and origin: single Fourier modes, skew and left-handed cells, negative Nyquist, translated fields, Parseval, cell rows versus columns.
3. Charge sign and Poisson: a positive source gives positive phi and negative U; constant density gives no Hartree field; q/Omega normalization; net charge rejected or background stated; signed input never accepted as electron number density.
4. Weak-field response: central derivative of the stationary energy equals `integral n W` on the same grid and k set; uncertainty from measured energy noise and density-integral noise; avoid gap closure.
5. Double counting: a fixed field counted once in the direct and band forms; `local_zero_mode` not added twice; a known density-independent ion/environment constant counted in the total once and not in the Hamiltonian.
6. Serialization round trip: byte-identical payloads, hash, length and order checks, recomputed electron count; corrupted or truncated payloads rejected.
7. Metadata rejection: cell, grid, k weights and order, spin, charge, pseudopotential digest, occupation and XC model, potential role and reference, environment identity.
8. Lifecycle state machine (mock): id deduplication, conflicting payload, cancellation before and after a completed step, partial results never marked converged, separate inner and outer status.
9. Forces and stress: finite differences only after a complete stationary functional and a derivative evaluator exist; tolerances calibrated from displaced or strained energy noise and basis errors.
10. CPU and GPU neutrality: physical contracts expressed independently of backend array types.
Tolerances are derived from measured noise floors in `docs/validation.md` (fixture specific; recalibrated per fixture), never copied between systems; long runs are `#[ignore]` with release commands (AGENTS.md).
