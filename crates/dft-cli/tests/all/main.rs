//! One test harness binary for the integration tests of `dft-cli` that hold
//! no process-global state (no environment mutation, forced switches,
//! per-thread override hooks, process-global counters or registries, static
//! locks, global allocators). Each former test target is a module here: run
//! one with `cargo test -p dft-cli --test all <module>::`. Targets that do hold
//! process-global state stay separate binaries; the list with the boundary
//! that justifies each is docs/performance/test-targets.md (TEST-HARNESS,
//! 2026-11-18).
// Each former test binary declared its own fixture modules (`common`, `#[path]`
// includes); kept as they were, so a fixture file may be loaded by several modules.
#![allow(clippy::duplicate_mod)]
mod atomic_start;
mod bands_dos;
mod buffer_bands;
mod cell_relaxation;
mod commands;
mod inspect_upf;
mod materials;
mod preflight;
mod pseudo;
mod reference_mode_config;
mod reference_mode_equivalence;
mod reference_mode_relaxation;
mod relax_extrapolation;
mod relaxation;
mod run;
mod run_config;
mod self_healing;
mod stress;
mod symmetry;
mod volumetric;
