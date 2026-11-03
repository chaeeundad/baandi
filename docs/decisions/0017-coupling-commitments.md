# ADR 0017: Coupling commitments now and deferrals

- Date: 2026-11-03
- Status: proposed (design only)
- Context: the user's general-purpose coupling directive (docs/research/coupling-directive.md); independent review docs/research/coupling-codex-review.md; design docs/architecture/scientific-coupling.md.

## Commit now (semantics, no code)
1. External fields are typed by physical role (electron potential energy in Ha versus electrostatic potential in Ha per e versus signed charge versus electron number density), with sign convention q_e = -1 and `U_el = -phi`.
2. Unsupported capabilities, observables and combinations are explicit typed errors or typed "unavailable" results, never ignored and never zero-filled.
3. Inner SCF convergence and outer coupled convergence are separate statuses.
4. Every contribution to an energy has a named owner; an electronic partial energy is never reported as a coupled total or free energy; coupled states are not ranked by DFT energy alone.
5. Environment identity is immutable and part of the ESLE model fingerprint; request identity is idempotent and distinct from state and branch identity.
6. The numerical core stays free of serialization, networking, file formats, workflow scheduling, MD, MLIP, FEM, solvent models and any specific framework; coupling lives in contracts inside existing crates and in external adapters.
7. No new crates now (ownership.toml already provides the boundaries); additive changes only; zero-field requests follow the existing path unchanged.

## Deferred
Concrete ABI and runtime, restart guarantees, charged-cell and background models, nonlocal and density-dependent embedding, derivatives (forces, stress) and their coupling terms, array import, transport adapters, optional format dependencies (extxyz, HDF5, openPMD), GPU contracts beyond backend-neutral semantics.

## Rejected
One oversized `EnergyForceCalculator` trait (would dictate SCF internals); forcing every environment into a scalar potential; using `V_eff` as an exported electrostatic potential; treating a periodic scalar field as a uniform electric field; JSON for dense fields.

## Consequences
Gate 2 of docs/roadmap/coupling.md needs the user's decisions listed there. The first implementation will touch the SCF assembly, both energy paths and the fixed-point certificate and will have to keep every existing k-point result unchanged.
