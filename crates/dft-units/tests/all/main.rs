//! One test harness binary for the integration tests of `dft-units` that hold
//! no process-global state (no environment mutation, forced switches,
//! per-thread override hooks, process-global counters or registries, static
//! locks, global allocators). Each former test target is a module here: run
//! one with `cargo test -p dft-units --test all <module>::`. Targets that do hold
//! process-global state stay separate binaries; the list with the boundary
//! that justifies each is docs/performance/test-targets.md (TEST-HARNESS,
//! 2026-11-18).
// Each former test binary declared its own fixture modules (`common`, `#[path]`
// includes); kept as they were, so a fixture file may be loaded by several modules.
#![allow(clippy::duplicate_mod)]
mod expressions;
mod units;
