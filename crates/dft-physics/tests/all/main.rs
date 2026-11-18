//! One test harness binary for the integration tests of `dft-physics` that hold
//! no process-global state (no environment mutation, forced switches,
//! per-thread override hooks, process-global counters or registries, static
//! locks, global allocators). Each former test target is a module here: run
//! one with `cargo test -p dft-physics --test all <module>::`. Targets that do hold
//! process-global state stay separate binaries; the list with the boundary
//! that justifies each is docs/performance/test-targets.md (TEST-HARNESS,
//! 2026-11-18).
// Each former test binary declared its own fixture modules (`common`, `#[path]`
// includes); kept as they were, so a fixture file may be loaded by several modules.
#![allow(clippy::duplicate_mod)]
mod hamiltonian;
mod pbe;
mod pbe_pw92_materials;
mod pbe_pw92_oracles;
mod radial_convention_diag;
mod radial_convention_score;
mod real_gamma_hamiltonian;
mod scf_helpers;
mod spin_grid;
mod spin_oracles;
mod stress_terms;
mod symmetry;
mod volumetric;
mod xc_prefactor_diag;
