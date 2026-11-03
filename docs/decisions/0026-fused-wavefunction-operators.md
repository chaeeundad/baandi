# ADR 0026: Sphere-local FFT operators with worker-owned scratch

Date: 2026-11-03. Status: implemented for review; performance acceptance on the
quiet node remains with the orchestrator. Scope: P10C-FUSED-FFT, pure Rust.

## Context

The local Hamiltonian application scatters to a full reciprocal grid, copies
that grid for each transform, and multiplies the potential in a separate pass.
Density also scatters and transforms unoccupied bands. The prose study note
`docs/research/qe-study/fft-call-structure.md` motivates compact sticks and
plane-local work; no QE source was read for this change.

The current inverse and forward both execute x/y/z, with reciprocal-axis
scaling after each inverse lane and total-size multiplication at the end.
A symmetric z/y/x inverse followed by x/y/z forward changes floating-point
results. The measured alternative meets the normwise and SCF energy bounds
but changes density hashes. Bit identity remains the default contract.

## Decision

Add `dft_numerics::fused::FusedOperator`, with plan-owned scratch and two
operations: projected local-potential application and weighted density
accumulation. The plan is constructed from dimensions and ordered sphere
indices. C-order real arrays and the input/output coefficient order are
explicit. Repeated indices use last-write-wins scatter semantics. Weight and
volume remain separate arguments to preserve multiplication/division order.
Caller-owned outputs and validation-before-write avoid per-call allocation.
The physics layer validates finite orbitals and potentials.

`OrderedFusedPlan` is the pure-Rust default. Its pipeline is:

1. Scatter to compact x sticks and inverse x, retaining per-lane scaling.
2. Inverse y on active z planes, storing a compact intermediate.
3. For each real y, inverse z into one xz plane, multiply the potential, and
   immediately transform x forward. Gather only requested x columns.
4. Forward y and z on the compact requested support; gather and divide by N.

This asymmetric fusion preserves each value's x/y/z arithmetic in both
transforms. No full real-space wavefunction grid, full-grid copy or full-grid
zero-fill is needed. At full support, compact intermediates necessarily become
dense. Density stops after inverse z and accumulates in the original arithmetic
order. Exactly zero occupations skip the FFT after caller validation.

`PlaneFusedPlan` provides the requested symmetric z-stick/xy-plane experiment
behind the same trait. It batches contiguous z sticks and processes one xy
plane at a time. The x pass writes a separate plane; unrequested x columns
remain zero. Y pruning is requested only at <=60% active columns. The dense
control uses separable 1D transforms because rustfft has no native 2D plan.
This backend is explicit and is not wired as the default.

One plan belongs to each existing band worker. No FFT threads, global caches,
new dependencies, unsafe code or production timers are introduced. The MKL
stage-2 backend can implement the trait with cached batched 1D descriptors and
its own stride/distance choices. The compile-time `FusedPlan` alias preserves P10B selection: default Rust uses
`OrderedFusedPlan`; `mkl` uses `GridFusedPlan<FftPlan3d>` around the selected
native 3D backend until stage 2 is measured. This compatibility adapter owns
reusable full-grid scratch; it is not a vendor fused optimization. No MKL
implementation file or FFI is changed.

## Evidence and consequences

See `docs/performance/batched-fft-study.md` for exact commands, oracle coverage,
performance samples and rejected alternatives. The existing generic `Fft3d`
density API remains available; `from_orbitals_fused` is additive, and the engine
uses persistent fused plans. The only crate-root change is `pub mod fused;`,
authorized by the orchestrator. An initial ordered prototype retained one full
real grid; it was superseded by the asymmetric plane pipeline above.

Group sizes 1/2/4/8/16/32 are measured with actual stage grouping and shared
potential loads in test-only code. No grouping default is adopted from noisy
single-machine measurements. This ADR adds no CPU-specific size threshold or
Si-specific numerical parameter.

## Integration on current main

Main `9d6fef14ce86dce821bd359c6f58fef0aecad8ab` includes the optional MKL
backend, buffer bands, RMM-DIIS and preconditioned mixing. The authorized
constructor-only change in `dft-engine/src/scf.rs` creates `BandDensity` for
every inner thread budget, including one, without changing any solver, mixing
or thread policy. This ADR was renumbered from the uncommitted 0021 draft to
0026 to avoid the volumetric-output ADR imported from main.
