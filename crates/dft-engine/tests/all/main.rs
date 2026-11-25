//! One test harness binary for the integration tests of `dft-engine` that hold
//! no process-global state (no environment mutation, forced switches,
//! per-thread override hooks, process-global counters or registries, static
//! locks, global allocators). Each former test target is a module here: run
//! one with `cargo test -p dft-engine --test all <module>::`. Targets that do hold
//! process-global state stay separate binaries; the list with the boundary
//! that justifies each is docs/performance/test-targets.md (TEST-HARNESS,
//! 2026-11-18).
// Each former test binary declared its own fixture modules (`common`, `#[path]`
// includes); kept as they were, so a fixture file may be loaded by several modules.
#![allow(clippy::duplicate_mod)]
mod atomic_start_contract;
mod atomic_start_option;
mod atomic_start_prototype;
mod band_boundary_scan;
mod band_boundary_spectrum;
mod band_count_spectrum;
mod band_fraction_scan;
mod capacity_contract;
mod checkpoint;
mod collinear;
// `common/checkpoint.rs` (included by several modules) reaches the shared
// fixtures as `crate::common`.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
mod cutoff_ladder;
mod density_transfer;
mod diagnostics;
mod electrostatic;
mod force_identity;
mod forces;
mod frozen_cell_metric;
mod gamma_real_diagnosis;
mod gamma_real_materials;
mod gamma_real_scf;
mod gamma_workspace_snapshot;
mod kpoints;
mod mixing_survey;
mod nlcc;
mod nscf;
mod pbe;
mod preconditioner_contract;
mod real_nc;
mod run_config;
mod runtime_policy;
mod scf;
mod spin_atomic_start;
mod stepwise_equivalence;
mod stress;
mod symmetry;
mod symmetry_full_group;
mod threading;
mod trajectory_divergence;
