# ADR 0016: ESLE layering and additive-only state interfaces

- Date: 2026-11-03
- Status: proposed (design only; no code changes)
- Context: the user's ESLE directive (docs/research/esle-directive.md) asks for multiple-solution search and branch tracking on top of the validated SCF. The engine currently has uniform-density starts only, no state serialization, no spin and no smearing (see docs/architecture/electronic-state.md section 1).

## Decisions
1. The numerical layers (basis, Hamiltonian, eigensolver, mixer, SCF loop) stay independent of persistence, registries and exploration policy. ESLE code sits above them.
2. Phase B interface changes are additive: a new `ScfProblem::run_from(StartState)` beside the unchanged `run`; legacy first-k fields of `ScfResult` stay. No breaking change during the k-point validation.
3. State types go into a new `dft-state` crate (data, serialization, compatibility, provenance) depending only on units, geometry and the serde stack; exploration code goes into `dft-landscape` (Phase C). At most two new crates; both are added to `ownership.toml` when their phase starts.
4. Energies and spectra are comparable only between states with an identical `model_fingerprint`; incompatible comparisons return a typed error, never a silent result.
5. State classification is three-way (SAME, DISTINCT, UNRESOLVED) with tolerances derived from measured noise floors; energy similarity alone never merges states; no search result is described as a global minimum.
6. Initialization strategies are enabled only if the implemented physics supports them (no spin starts without spin-polarized DFT, no occupation control without DFT+U). Spin, smearing and DFT+U are separate, separately validated features (Phase E) and need user approval and ADRs.

## Alternatives rejected
- Embedding checkpointing inside `ScfProblem`: couples kernels to I/O and blocks GPU backends sharing state semantics.
- One crate for everything: harder to enforce the dependency direction in `ownership.toml`.
- Identity by wavefunction equality or by energy: not gauge invariant / merges distinct states.

## Consequences
Phase B starts only after Phase A gates. A scientific review checklist applies to each ESLE ticket (docs/roadmap/esle.md). Open items: tolerance factors (measured in Phase B), exposure of `PP_RHOATOM` in the UPF reader, the second pseudopotential for MgO.
