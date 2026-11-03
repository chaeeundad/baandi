# ADR 0016: ESLE layering and additive-only state interfaces

- Date: 2026-11-03
- Status: proposed (design only; no code changes)
- Context: the user's ESLE directive (docs/research/esle-directive.md) asks for multiple-solution search and branch tracking on top of the validated SCF. The engine currently has uniform-density starts only, no state serialization, no spin and no smearing (see docs/architecture/electronic-state.md section 1).

## Decisions
1. The numerical layers (basis, Hamiltonian, eigensolver, mixer, SCF loop) stay independent of persistence, registries and exploration policy. ESLE code sits above them.
2. Phase B interface changes are additive: a new `ScfProblem::run_from(StartState)` beside the unchanged `run`; legacy first-k fields of `ScfResult` stay. No breaking change during the k-point validation.
3. No new crates initially (revised after the independent review, which found that `ownership.toml` already provides the needed boundaries): pure state types, compatibility, comparison and exploration policy are new modules in `dft-engine` (no file access); persistence and the driving loop live in `dft-cli`; density transfer and atomic-density physics go into the existing `dft-physics` and `dft-upf` when separately ticketed. A `dft-state` or `dft-landscape` crate is extracted only when an independent client justifies an ownership change (a later ADR).
4. Energies and spectra are comparable only between states with an identical `model_fingerprint`; incompatible comparisons return a typed error, never a silent result.
5. State classification is three-way (SAME, DISTINCT, UNRESOLVED) with tolerances derived from measured noise floors; energy similarity alone never merges states; no search result is described as a global minimum.
6. Initialization strategies are enabled only if the implemented physics supports them (no spin starts without spin-polarized DFT, no occupation control without DFT+U). Spin, smearing and DFT+U are separate, separately validated features (Phase E) and need user approval and ADRs.

## Alternatives rejected
- Embedding checkpointing inside `ScfProblem`: couples kernels to I/O and blocks GPU backends sharing state semantics.
- New `dft-state` and `dft-landscape` crates now (the first draft proposed them): rejected because the directive asks to avoid crate proliferation and `ownership.toml` already separates pure engine logic from CLI I/O.
- Identity by wavefunction equality or by energy: not gauge invariant / merges distinct states.

## Consequences
Phase B starts only after Phase A gates. A scientific review checklist applies to each ESLE ticket (docs/roadmap/esle.md). Open items: tolerance factors (measured in Phase B), a validated normalized atomic-density adapter over `Pseudopotential::raw()`, the second pseudopotential for MgO, archive orbital-retention policy, the exact-replay versus warm-restart promise, the initial supported cell-identification maps. Independent review: docs/research/esle-codex-review.md (Codex run P8R-ESLE-REVIEW); its main corrections are in docs/research/electronic-state-landscape.md section 7.
