# ElectronicState: representation, serialization, compatibility (proposal)

Status: proposal, reconciled with the independent Codex review (`../research/esle-codex-review.md`); nothing here is implemented and no existing interface changes until Phase B (see [../roadmap/esle.md](../roadmap/esle.md), ADR [0016](../decisions/0016-esle-layering.md)). Statements about current code are marked "verified" with a file.

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
| Atomic density | UPF `PP_RHOATOM` is parsed and checked and reachable through `Pseudopotential::raw()`; no normalized accessor or superposition initializer exists | `crates/dft-upf/src/parser.rs` |
| Determinism | fixed seed in options; single execution thread recorded by the CLI; no other source of nondeterminism found yet (to be certified by Phase B tests, not assumed) | `ScfOptions::seed`, `docs/cli.md` |

## 2. Layering (revised after the independent review)
Numerical kernels (basis, Hamiltonian, eigensolver, mixer, SCF loop) stay free of persistence, registries and policies. No new crates initially (ADR 0016): pure state types, compatibility checks, comparison and exploration policy are modules in `dft-engine` (no file access, allowed by `ownership.toml`); state persistence (immutable manifest committed after payload writes) and the driving loop live in `dft-cli`; density transfer and atomic-density physics are separate tickets in `dft-physics` and `dft-upf`. A crate is extracted only when an independent client justifies an ownership change. The first draft proposed `dft-state` and `dft-landscape`; that was dropped because the repository already has the needed boundaries.
Four distinct contracts, not one: (1) `InitialGuess` (a start: density and/or orbital guesses plus its recipe; it has no energy or eigenvalues), (2) `EvaluatedState` (the outcome of a run: converged or not; only converged and certified states are eligible for ranking), (3) `AttemptRecord` (every attempt, including failures with no state), (4) the archive `ElectronicState` (versioned, validated; may omit orbital payloads and then states that capability loss explicitly). An exact trajectory checkpoint (mixer history, subspaces) is a fifth, optional thing and is not promised.
Gaps in the current outcome that the envelope must fill: `ScfResult` carries no geometry, pseudopotential hashes, basis labels or build identity (verified `scf.rs`, `ScfResult`).

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
| Compare two states (distance metrics; result is Exact, Convertible(plan) or Incompatible(reasons) per operation) | `model_fingerprint` for energy and spectrum metrics | none; density metrics may also compare across grids through the common G sphere if the cell is equal, stated in the report |
| Restart from density, same model | `model_fingerprint` | none |
| Restart from wavefunctions | `model_fingerprint`, basis ordering | none (a different cutoff or cell needs the conversion below) |
| Density transfer to another structure | pseudopotentials, functional, `nspin`, electron count, occupation policy | cell, positions, grid, cutoff, k set: the density is resampled in G space, renormalized to the electron count, and the new state records the parent and `FromState{transfer}` |
Energies of states with different fingerprints are never ranked; a report may show them side by side with a visible "different model" marker.

## 6. Interface sketches (not code; additive only)
```rust
// dft-engine (additive, Phase B): the existing run() is unchanged; None keeps the uniform start.
pub enum InitialState {
    Density { density: Array3<f64>, basis: StateBasis },
    Orbitals { blocks: Vec<OrbitalBlock>, basis: StateBasis },
    DensityAndOrbitalGuess { density: Array3<f64>, blocks: Vec<OrbitalBlock>, basis: StateBasis },
}
pub struct ScfRequest { pub initial: Option<InitialState> }
impl ScfProblem { pub fn run_from(self, request: ScfRequest) -> Result<ScfOutcome, ScfError>; }

// dft-engine state module (pure, no file access)
pub struct ElectronicState { /* fields of section 3 */ }
pub fn from_outcome(problem_model: &ModelSpec, outcome: &ScfOutcome, prov: Provenance) -> ElectronicState;
pub fn save(state: &ElectronicState, dir: &Path) -> Result<(), StateError>;
pub fn load(dir: &Path) -> Result<ElectronicState, StateError>;
pub fn check_compatible(a: &ElectronicState, b: &ElectronicState, use_: Use) -> Result<(), Incompatible>;

// dft-engine landscape module (pure) and a dft-cli StateStore for persistence
pub trait StateDistance { fn distance(&self, a: &ElectronicState, b: &ElectronicState) -> Result<DistanceReport, Incompatible>; }
pub trait Initializer { fn start(&self, model: &ModelSpec, parent: Option<&ElectronicState>) -> Result<StartState, InitError>; }
pub trait ExplorationPolicy { fn next(&mut self, history: &Registry, budget: &mut Budget) -> Option<Attempt>; }
```
Registry records every attempt (initialization, outcome, failure, cost) and never discards non-converged runs. The Hamiltonian and eigensolver never see these types.

## 7. Breaking versus additive
All Phase B changes are additive (`run_from` and new modules); adding a required field to public `ScfOptions` or `ScfResult`, replacing the first-k fields or changing the `run` return type would be breaking and is avoided. The one place that needs care is `ScfResult.occupations`/`orbitals` (legacy first-k fields): they stay as they are. Nothing in the current k-point validation is touched.
