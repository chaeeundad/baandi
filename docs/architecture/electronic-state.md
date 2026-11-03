# ElectronicState: representation, serialization, compatibility (proposal)

Status: proposal; nothing here is implemented and no existing interface changes until Phase B (see [../roadmap/esle.md](../roadmap/esle.md), ADR [0016](../decisions/0016-esle-layering.md)). Statements about current code are marked "verified" with a file.

## 1. Current state of the code (verified)
| Topic | Fact | Where |
|---|---|---|
| SCF input | `ScfProblem::with_kpoints(cell, atoms, potentials, points, options)`; `new` is the Gamma wrapper | `crates/dft-engine/src/scf.rs` |
| Start density | always uniform `N_el/Omega`; no way to pass another start to `run` | `uniform_density`, `run_impl` |
| Seam for an initial state | `step(n, initial_orbitals, tol)` and `step_kpoints` already take an input density and optional per-k orbitals; the driving loop is internal | `scf.rs` |
| Result | `ScfResult`: `EnergyAudit`, per-k `KPointSolution` (point, `EigenSolution`), `occupations` (first k, legacy), `density: Array3<f64>`, `history: Vec<IterationRecord>`, timings | `scf.rs` |
| Occupations | fixed insulator rule (f = 2 for lowest N_el/2 bands at every k); partial occupations rejected by `OccupationPolicy` | `scf.rs` |
| Spin | none (`nspin = 1` implied) | no spin field anywhere |
| Wavefunctions | complex coefficients per k in that k's plane-wave order (`ComplexMatrix`); kept in memory only | `scf.rs`, `dft-basis` |
| Persistence | CLI writes `results.json`, `density.bin` (f64 LE, C order) with `density.json` (grid, units, SHA-256); no wavefunction file, no reader, no restart | `crates/dft-cli/src/run.rs`, `docs/cli.md` |
| Atomic density | UPF `PP_RHOATOM` is validated but not exposed | `crates/dft-upf/src/parser.rs` |
| Determinism | fixed seed in options; single execution thread recorded by the CLI; no other source of nondeterminism found yet (to be certified by Phase B tests, not assumed) | `ScfOptions::seed`, `docs/cli.md` |

## 2. Layering
Numerical kernels (basis, Hamiltonian, eigensolver, mixer, SCF loop) stay free of persistence, registries and policies. New code lives in two new layers added only when their phase starts:
- `dft-state` (Phase B): data types, versioned serialization, compatibility checks, provenance. Dependencies: dft-units, dft-geometry, serde stack only (no solver crates), so any producer and consumer can use it.
- `dft-landscape` (Phase C): comparison, registry, initialization strategies, multi-start, tracking, exploration policy. Depends on dft-state and dft-engine.
Two crates at most; `dft-validation` is not a crate (tests and scripts already exist: `tests/reference`, `scripts/`). Both are proposals subject to `ownership.toml` changes made when the phase starts.

## 3. ElectronicState v1 fields
M = mandatory in v1 (a state without it is invalid); O = optional; F = future (reserved name, rejected if present before the physics exists).
| Group | Field | Level | Notes |
|---|---|---|---|
| Header | `schema_version` (major.minor), `state_id`, `created_utc` | M | `state_id` = SHA-256 of the canonical manifest without `created_utc` and without provenance host fields |
| Model (M) | cell lattice (Bohr), atoms (species, fractional positions) | M | |
| | pseudopotential ids and SHA-256 per species | M | hash of the exact file bytes; the file itself is never stored |
| | functional label, `nspin` (1), occupation policy, electron count | M | |
| | cutoff (Ha), FFT grid, basis ordering rule id | M | |
| | k set: fractional coordinates, weights, reduction rule (e.g. time reversal only) | M | |
| | `model_fingerprint` = SHA-256 of the model fields above | M | energies are comparable only for equal fingerprints |
| Solution | total energy and the full `EnergyAudit` components (Ha) | M | |
| | charge density (array reference + shape + SHA-256) | M | |
| | per-k eigenvalues (Ha) and occupations | M | |
| | convergence record: final residuals, iteration count, tolerances used, converged flag | M | history may be summarized; full history optional |
| | wavefunction coefficients per k | O | large; needed for subspace metrics and orbital warm starts |
| | Fermi level / gap | O | gap exists for insulators |
| Provenance | solver version (git commit, dirty flag), toolchain, lockfile hash, thread count, seed, command line | M | |
| | parent `state_id`, initialization method with parameters | M (null parent allowed) | `Uniform`, `AtomicSuperposition`, `FromState`, `RandomPerturbation{seed, amplitude}`, ... only methods valid for the model are accepted |
| Future | spin density, local moments, occupation matrices, smearing/temperature, Hubbard parameters, symmetry operations | F | adding them bumps the minor version; models using them get a different fingerprint |

## 4. Serialization
- Manifest: JSON (`serde_json`, ADR 0010 stack), `deny_unknown_fields`, strict types, finite floats only, keys in a fixed order, floats written with round-trip digits.
- Arrays: little-endian f64 (complex as re,im pairs) in sibling files with a SHA-256 in the manifest; layout and ordering stated per array (density: C order, k fastest, same as `density.json`).
- Validation on load: schema version supported, hashes match, shapes consistent with the model (grid, k count, band count), electron count from the density integral equals the stated count within a stated tolerance, occupations sum to the electron count, `state_id` and `model_fingerprint` recomputed and equal. Any failure is a typed error; nothing is repaired silently.
- Round trip test: save then load gives bit-identical arrays and manifest; this is part of the Phase B gate.

## 5. Compatibility rules
| Use | Must be equal (else error) | May differ (explicit conversion step, recorded) |
|---|---|---|
| Compare two states (distance metrics) | `model_fingerprint` for energy and spectrum metrics | none; density metrics may also compare across grids through the common G sphere if the cell is equal, stated in the report |
| Restart from density, same model | `model_fingerprint` | none |
| Restart from wavefunctions | `model_fingerprint`, basis ordering | none (a different cutoff or cell needs the conversion below) |
| Density transfer to another structure | pseudopotentials, functional, `nspin`, electron count, occupation policy | cell, positions, grid, cutoff, k set: the density is resampled in G space, renormalized to the electron count, and the new state records the parent and `FromState{transfer}` |
Energies of states with different fingerprints are never ranked; a report may show them side by side with a visible "different model" marker.

## 6. Interface sketches (not code; additive only)
```rust
// dft-engine (additive, Phase B): the existing run() is unchanged.
pub enum InitialDensity { Uniform, Given(Array3<f64>) }
pub struct StartState { pub density: InitialDensity, pub orbitals: Option<Vec<ComplexMatrix>> }
impl ScfProblem { pub fn run_from(self, start: StartState) -> Result<ScfOutcome, ScfError>; }

// dft-state
pub struct ElectronicState { /* fields of section 3 */ }
pub fn from_outcome(problem_model: &ModelSpec, outcome: &ScfOutcome, prov: Provenance) -> ElectronicState;
pub fn save(state: &ElectronicState, dir: &Path) -> Result<(), StateError>;
pub fn load(dir: &Path) -> Result<ElectronicState, StateError>;
pub fn check_compatible(a: &ElectronicState, b: &ElectronicState, use_: Use) -> Result<(), Incompatible>;

// dft-landscape
pub trait StateDistance { fn distance(&self, a: &ElectronicState, b: &ElectronicState) -> Result<DistanceReport, Incompatible>; }
pub trait Initializer { fn start(&self, model: &ModelSpec, parent: Option<&ElectronicState>) -> Result<StartState, InitError>; }
pub trait ExplorationPolicy { fn next(&mut self, history: &Registry, budget: &mut Budget) -> Option<Attempt>; }
```
Registry records every attempt (initialization, outcome, failure, cost) and never discards non-converged runs. The Hamiltonian and eigensolver never see these types.

## 7. Breaking versus additive
All Phase B changes are additive (`run_from`, new crates). The one place that needs care is `ScfResult.occupations`/`orbitals` (legacy first-k fields): they stay as they are. Nothing in the current k-point validation is touched.
